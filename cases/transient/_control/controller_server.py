"""Boundary-control worker for SIREN's `*-capi` transient references.

One process, one SIREN run. Speaks the newline-delimited protocol of
`siren/src/transient_albedo.c` over a Unix-domain socket:

    recv "<time> <albedo_coeff> <pow_lhs> <pow_rhs>\n"   ->   send "<albedo>\n"

(a three-field line from a pre-two-sided binary is still accepted, with
pow_rhs as its last field.)

All controller modes run through this *same* bridge, first-step handling, log
schema and saturation accounting, so a difference between two modes is
attributable to the control law and not to the harness:

    echo        return the albedo SIREN passed in (open-loop transparency gate)
    pid-legacy  faithful reimplementation of SIREN's own C PID
                (kp=0.004, ki=0, kd=0.001) -- the as-shipped baseline
    pid-tuned   same velocity-form structure, gains supplied per case
    pi-only     AdaptiveBoundaryController(backstepping=False): the PI trim,
                observer and identifier, with the field feedback held at its
                design anchor -- the decisive ablation
    full        the method: PI trim + adaptive backstepping field feedback

Two-sided modes (need the 4-field bridge; refuse a 3-field binary):

    observe-2s  pid-tuned drives the plant; the two-output observer and
                swapping identifier run passively on (pow_lhs, pow_rhs), so
                identification is scored independently of any new control law
    pid-mean    pid-tuned's structure on 0.5*(pow_lhs/pl0 + pow_rhs/p0)*p0 --
                a two-sensor total-power proxy with NO observer: the arm a
                total-power win has to beat to be credited to the observer
    pi-total    two-output observer; PI trim on its total-power estimate,
                field feedback held at the design anchor
    full-total  the same plus adaptive backstepping field feedback

Why the ablation matters: `documents/phase5_findings.md` records that on the
spatially-flat case the regulation metrics were insensitive to the backstepping
gain, i.e. the PI trim was doing the regulating. `pi-only` vs `full` is the
control that separates method from tuning.
"""

import argparse
import json
import os
import socket
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "..", "..", "src"))

MODES = ("echo", "pid-legacy", "pid-tuned", "pi-only", "full",
         "observe-2s", "pid-mean", "pi-total", "full-total")
TWO_SIDED = ("observe-2s", "pid-mean", "pi-total", "full-total")
OUTPUTS_2S = ("lhs", "rhs")


def _envf(name, default):
    v = os.environ.get(name)
    return float(v) if v is not None and v != "" else default


def parse_args(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", required=True, choices=MODES)
    ap.add_argument("--dt", type=float, required=True,
                    help="must match the deck's deltat")
    ap.add_argument("--sock", default="controller.sock")
    ap.add_argument("--out-dir", default=".")
    ap.add_argument("--snap-times", default="0.0,0.5,1.0,1.5,2.0",
                    help="times at which to dump phihat/lam_hat, matching tedit")
    ap.add_argument("--actuator-dx", type=float, default=None,
                    help="SIREN's dx at the actuated face; default = design grid h")
    ap.add_argument("--n-design", type=int, default=40)
    return ap.parse_args(argv)


# ----------------------------------------------------------------- PID modes
class VelocityPID:
    """SIREN's own controller, ported verbatim from `transient_albedo.c`.

    Velocity form: the output is added to the *previous* albedo each step, so
    the accumulation of kp*error already supplies integral action -- which is
    why the shipped gains set ki=0. Reproduced exactly, including the
    uncontrolled first step, so `pid-legacy` is the real shipped baseline.
    """

    def __init__(self, kp=0.004, ki=0.0, kd=0.001):
        self.kp, self.ki, self.kd = float(kp), float(ki), float(kd)
        self.raw = None      # pre-clamp command, for honest saturation counting
        self.p0 = None
        self.prev_error = 0.0
        self.prev_time = 0.0
        self.interror = 0.0

    def step(self, time_s, albedo_in, pboundary):
        alb = albedo_in
        if self.p0 is None:                       # first call: latch, no action
            self.p0 = pboundary
            self.raw = alb
            error = 0.0
        else:
            error = self.p0 - pboundary
            dt = time_s - self.prev_time
            derror_dt = (error - self.prev_error) / dt if dt != 0.0 else 0.0
            self.interror += error * dt
            alb += (self.kp * error + self.kd * derror_dt
                    + self.ki * self.interror)
            self.raw = alb
            alb = min(max(alb, -1.0), 1.0)
        self.prev_time = time_s
        self.prev_error = error
        return alb


def build_adaptive(args, backstepping, outputs=("rhs",)):
    from backstep1d.controllers.adaptive import AdaptiveBoundaryController
    from backstep1d.plants.siren_onegroup import (
        SirenOneGroupPlant, SirenOneGroupDesignModel)

    basis_env = os.environ.get("BS_BASIS", "hat")
    n_basis = int(_envf("BS_N_BASIS", 8))
    if basis_env == "physics":
        # Single psi(x)=(H-x)^2 basis function (see `observers/passive.py`);
        # BS_N_BASIS is meaningless here, so it is overridden rather than
        # left to silently mis-size theta0 below.
        n_basis = 1
    theta0_env = os.environ.get("BS_THETA0")
    theta0 = (np.full(n_basis, float(theta0_env))
              if theta0_env not in (None, "") else None)

    plant = SirenOneGroupPlant(N=args.n_design, thermal_bc="robin")
    dm = SirenOneGroupDesignModel(plant)
    ctrl = AdaptiveBoundaryController(
        plant, dm, dt=args.dt,
        c=_envf("BS_C", 200.0),
        n_basis=n_basis,
        c_obs=_envf("BS_C_OBS", 100.0),
        gamma=_envf("BS_GAMMA", 200.0),
        resolve_every=int(_envf("BS_RESOLVE_EVERY", 100)),
        dtheta_tol=_envf("BS_DTHETA_TOL", 1.0),
        gain_ramp=int(_envf("BS_GAIN_RAMP", 0)),
        lam_bounds=(_envf("BS_LAM_LO", -150.0), _envf("BS_LAM_HI", 20.0)),
        robin_alpha=0.0, theta0=theta0,
        ki=_envf("BS_KI", 5.0), kp=_envf("BS_KP", 0.02),
        backstepping=backstepping,
        actuator_dx=args.actuator_dx,
        # "direct" is the as-shipped adaptive law; "swapping" swaps its
        # regressor for the exact output sensitivity (see
        # `observers/passive.py`). Selected by environment rather than by
        # `--mode` so that every existing mode string, and therefore every
        # cached run directory, keeps its meaning.
        # Two outputs are only offered with the swapping law (see
        # `observers/passive.py`), so the two-sided modes force it.
        regressor=("swapping" if len(outputs) > 1
                   else os.environ.get("BS_REGRESSOR", "direct")),
        outputs=outputs,
        # "hat" (default) is the as-shipped n_basis-dimensional partition of
        # unity; "physics" is the single matched-shape diagnostic basis (see
        # `observers/passive.py`, `scripts/bc8_physics_basis.py`).
        basis=basis_env,
        # "term" (default) is the as-shipped per-term clamp; "sum" clamps only
        # the summed albedo command (see `controllers/adaptive.py`). Selected
        # by environment for the same cache-meaning reason as BS_REGRESSOR.
        bs_clip=os.environ.get("BS_CLIP", "term"),
        **drum_kw(args))
    return ctrl


def drum_kw(args):
    """DRUM_TAU [s] lag and DRUM_T180 [s] for a full 180-deg sweep (A from 1
    to -1, so |dA/dt| <= 2/T180). Unset = ideal actuator."""
    t180 = _envf("DRUM_T180", 0.0)
    return {"drum_tau": _envf("DRUM_TAU", 0.0),
            "drum_rate": 2.0 / t180 if t180 > 0 else float("inf")}


def guard(obj, name, counters, key):
    """Make an online kernel re-solve non-fatal.

    `ObserverGainKernel.gain` raises on a complex spectrum or a failed pole
    placement -- a real risk once lam_hat(x) is discontinuous. Losing a gain
    refresh degrades the observer; letting the exception escape aborts SIREN
    and throws away the whole run's log. Keep the previous gain and count it,
    so the count appears in the diagnostics rather than being silent.
    """
    orig = getattr(obj, name)

    def wrapped(*a, **k):
        try:
            return orig(*a, **k)
        except (RuntimeError, np.linalg.LinAlgError) as exc:
            counters[key] += 1
            counters.setdefault("last_error", str(exc))
            return None
    setattr(obj, name, wrapped)


def main(argv=None):
    args = parse_args(argv)
    out = os.path.abspath(args.out_dir)
    os.makedirs(out, exist_ok=True)
    snap_times = [float(s) for s in args.snap_times.split(",") if s.strip()]

    if os.path.exists(args.sock):
        os.remove(args.sock)
    srv = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
    srv.bind(args.sock)
    srv.listen(1)
    conn, _ = srv.accept()
    stream = conn.makefile("rw")

    counters = {"steps": 0, "nonfinite_meas": 0, "nonfinite_cmd": 0,
                "gain_fail": 0, "kernel_fail": 0, "sat_hi": 0, "sat_lo": 0}
    pid = None
    pid_cmd = None
    pid_drum = None
    ctrl = None
    phi_scale = None
    p0 = None
    pl0 = None
    alb_prev = None
    snaps = []
    next_snap = 0
    theta_trace = []
    # The identifier's update DIRECTION, traced on the same cadence as theta.
    # theta's own trajectory only shows the subspace it ended up in; the generic
    # regressor trace shows the subspace it was ever offered. For "direct" this
    # is <psi_i, phihat>/<psi_i>; for "swapping" it is xi_i = d yhat/d theta_i.
    # Keep xi_trace too for compatibility with existing swapping diagnostics.
    regressor_trace = []
    xi_trace = []
    # z(t) is the observer's own output error (measured boundary minus its
    # reconstruction). It is the single number that says whether the interior
    # estimate driving the backstepping term is any good, so it is traced
    # rather than only summarized -- bc3/bc5 gate on it.
    z_trace = []
    # What the bridge actually delivered, every step and in every mode: bc9
    # checks it bit-for-bit against SIREN's own controller_log.csv.
    recv_trace = []
    # Two-sided modes: the observer's total-power estimate, every step.
    phat_trace = []
    split_trace = []

    for line in stream:
        parts = line.split()
        # 4 fields since the two-sided bridge (time, A, pow_lhs, pow_rhs); 3
        # from a pre-lhs_observer binary. pboundary is always the last field,
        # so every existing mode sees exactly the value it saw before.
        if len(parts) == 4:
            time_s, albedo_in, pow_lhs, pboundary = (float(x) for x in parts)
        elif len(parts) == 3:
            time_s, albedo_in, pboundary = (float(x) for x in parts)
            pow_lhs = None
        else:
            continue
        recv_trace.append((time_s, albedo_in,
                           np.nan if pow_lhs is None else pow_lhs, pboundary))

        # ---- first sample: latch the operating point, build the controller
        if counters["steps"] == 0:
            p0 = pboundary
            pl0 = pow_lhs
            phi_scale = p0
            with open(os.path.join(out, "handshake.json"), "w") as f:
                json.dump({"mode": args.mode, "p0": p0, "pl0": pl0,
                           "phi_scale": phi_scale,
                           "A0": albedo_in, "t0": time_s, "dt": args.dt,
                           "actuator_dx": args.actuator_dx,
                           "n_design": args.n_design,
                           "regressor": os.environ.get("BS_REGRESSOR",
                                                       "direct")}, f, indent=2)
            if args.mode == "pid-legacy":
                pid = VelocityPID()
            elif args.mode == "pid-tuned":
                pid = VelocityPID(kp=_envf("PID_KP", 0.004),
                                  ki=_envf("PID_KI", 0.0),
                                  kd=_envf("PID_KD", 0.001))
            elif args.mode in ("pid-mean", "observe-2s"):
                pid = VelocityPID(kp=_envf("PID_KP", 0.004),
                                  ki=_envf("PID_KI", 0.0),
                                  kd=_envf("PID_KD", 0.001))
            if pid is not None:
                from backstep1d.actuator import DrumActuator
                pid_drum = DrumActuator(args.dt, A0=albedo_in, **{
                    {"drum_tau": "tau", "drum_rate": "rate"}[k]: v
                    for k, v in drum_kw(args).items()})
            if args.mode in TWO_SIDED and pl0 is None:
                raise SystemExit("controller_server: mode %s needs pow_lhs; "
                                 "SIREN sent the 3-field protocol" % args.mode)
            if args.mode in ("pi-only", "full"):
                ctrl = build_adaptive(args, backstepping=(args.mode == "full"))
            elif args.mode in ("observe-2s", "pi-total", "full-total"):
                ctrl = build_adaptive(args,
                                      backstepping=(args.mode == "full-total"),
                                      outputs=OUTPUTS_2S)
            if ctrl is not None:
                guard(ctrl.identifier, "_resolve_gain", counters, "gain_fail")
                guard(ctrl, "_resolve", counters, "kernel_fail")

        # ---- a diverging plant must not kill the run before the log is written
        if not np.isfinite(pboundary):
            counters["nonfinite_meas"] += 1
            alb = -1.0                       # maximum leakage; nothing else left
            counters["steps"] += 1
            stream.write("%.17g\n" % alb)
            stream.flush()
            alb_prev = alb
            continue

        # ---- controller
        # Velocity-form PID accumulates on its own previous COMMAND; with an
        # ideal drum that is exactly the albedo SIREN echoes back.
        pid_prev = albedo_in if pid_cmd is None else pid_cmd
        if args.mode == "echo":
            alb = albedo_in
        elif args.mode == "pid-mean":
            alb = pid_cmd = pid.step(time_s, pid_prev,
                                     0.5 * (pow_lhs / pl0 + pboundary / p0) * p0)
            alb = pid_drum.step(alb)
        elif pid is not None:
            alb = pid_cmd = pid.step(time_s, pid_prev, pboundary)
            alb = pid_drum.step(alb)
            if ctrl is not None and counters["steps"] > 0:
                # observe-2s: the identifier sees the albedo actually applied
                # (post-clamp), exactly as SIREN will apply it.
                ctrl.observe_two_sided(pow_lhs / phi_scale,
                                       pboundary / phi_scale,
                                       min(max(alb, -1.0), 1.0))
        else:
            if counters["steps"] == 0:
                alb = albedo_in              # same uncontrolled first step as PID
            elif args.mode in TWO_SIDED:
                ctrl.step_two_sided(pow_lhs / phi_scale, pboundary / phi_scale,
                                    objective="total")
                alb = ctrl.albedo_applied
            else:
                ctrl.step_collocated(pboundary / phi_scale)
                alb = ctrl.albedo_applied
        if ctrl is not None and ctrl.identifier.n_out > 1:
            phat_trace.append(ctrl.total_power_estimate())

        if not np.isfinite(alb):
            counters["nonfinite_cmd"] += 1
            alb = alb_prev if alb_prev is not None else -1.0

        # Requested (pre-clamp) command: the clamp is only binding if this left
        # [-1, 1]. Distinguishes a genuinely saturated actuator from a
        # controller sitting idle at the mirror value A = 1.
        raw = pid.raw if pid is not None else (
            getattr(ctrl, "A_raw_last", None) if ctrl is not None else None)
        if raw is not None and np.isfinite(raw):
            if raw > 1.0:
                counters["sat_hi"] += 1
            elif raw < -1.0:
                counters["sat_lo"] += 1
        alb = min(max(alb, -1.0), 1.0)

        # ---- interior-field snapshots for bc5 validation
        if ctrl is not None and next_snap < len(snap_times) \
                and time_s >= snap_times[next_snap] - 1e-12:
            snaps.append({"t": time_s,
                          "phihat": np.asarray(ctrl.identifier.phihat).copy(),
                          "lam_hat": np.asarray(ctrl.identifier.lam_hat()).copy()})
            next_snap += 1
        if ctrl is not None and counters["steps"] % 100 == 0:
            ident = ctrl.identifier
            theta_trace.append([time_s] + list(ident.theta))
            if ident.swapping:
                # two outputs: the 2 x M regressor, flattened row-major
                regressor = np.ravel(ident.xi)
                xi_trace.append([time_s] + list(regressor))
            else:
                overlap = ident.Psi[:, :ident.grid.N] @ (
                    ident.phihat * ident.grid.wq[:ident.grid.N])
                regressor = overlap / ident.psi_mass
            regressor_trace.append([time_s] + list(regressor))
        if ctrl is not None:
            z_trace.append(ctrl.identifier.z_last)
            split_trace.append((getattr(ctrl, "A_bs_last", np.nan),
                                getattr(ctrl, "A_trim_last", np.nan)))

        counters["steps"] += 1
        alb_prev = alb
        stream.write("%.17g\n" % alb)
        stream.flush()

    # ---- SIREN closed the socket: the run is over, persist diagnostics
    diag = dict(counters)
    diag.update({"mode": args.mode, "p0": p0, "dt": args.dt})
    if ctrl is not None:
        zt = np.asarray(z_trace, dtype=float)
        zf = zt[np.isfinite(zt)]
        ntail = max(1, zf.size // 4)
        diag.update({"regressor": ctrl.identifier.regressor,
                     "basis": ctrl.identifier.basis,
                     "n_resolves": ctrl.n_resolves,
                     "n_saturated_ctrl": ctrl.n_saturated,
                     # Fraction of steps on which the backstepping command was
                     # unrealizable (asked the actuator for less leakage than a
                     # perfect mirror) and was therefore clipped to the design
                     # anchor -- contributing exactly nothing. At ~1 the method
                     # is inert and any regulation is the PI trim's.
                     "n_bs_clipped": int(ctrl.n_bs_clipped),
                     "bs_clip_frac": float(ctrl.n_bs_clipped
                                           / max(1, counters["steps"])),
                     "bs_clip_mode": ctrl.bs_clip,
                     # bs_clip="sum": how often the term asked for less
                     # leakage than the anchor (and was KEPT, not clipped)
                     "bs_over_anchor_frac": float(ctrl.n_bs_over_anchor
                                                  / max(1, counters["steps"])),
                     "z_rms": float(np.sqrt(np.mean(zf ** 2))) if zf.size else None,
                     "z_rms_tail": float(np.sqrt(np.mean(zf[-ntail:] ** 2)))
                                   if zf.size else None,
                     "n_theta_reject": int(getattr(ctrl.identifier,
                                                   "n_theta_reject", 0)),
                     "grid_h": float(ctrl.grid.h),
                     "actuator_dx_used": float(ctrl.actuator_dx),
                     "theta_final": [float(v) for v in ctrl.identifier.theta],
                     "A_int_final": float(ctrl.A_int),
                     "albedo_final": float(ctrl.albedo_cmd)})
    with open(os.path.join(out, "diag.json"), "w") as f:
        json.dump(diag, f, indent=2)
    np.save(os.path.join(out, "recv_trace.npy"),
            np.asarray(recv_trace, dtype=float).reshape(-1, 4))
    if snaps:
        np.savez(os.path.join(out, "snapshots.npz"),
                 t=np.array([s["t"] for s in snaps]),
                 phihat=np.array([s["phihat"] for s in snaps]),
                 lam_hat=np.array([s["lam_hat"] for s in snaps]),
                 x=np.asarray(ctrl.grid.x),
                 theta_trace=np.array(theta_trace) if theta_trace else np.zeros((0, 0)),
                 regressor_trace=(np.array(regressor_trace)
                                  if regressor_trace else np.zeros((0, 0))),
                 xi_trace=np.array(xi_trace) if xi_trace else np.zeros((0, 0)),
                 z_trace=np.asarray(z_trace, dtype=float),
                 split_trace=np.asarray(split_trace, dtype=float),
                 phat_trace=np.asarray(phat_trace, dtype=float))
    print("controller_server: mode=%s steps=%d %s"
          % (args.mode, counters["steps"],
             {k: v for k, v in counters.items() if k != "steps" and v}),
          file=sys.stderr)


if __name__ == "__main__":
    main()
