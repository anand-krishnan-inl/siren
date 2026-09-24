"""Boundary-control worker for SIREN's `*-capi` transient references.

One process, one SIREN run. Speaks the newline-delimited protocol of
`siren/src/transient_albedo.c` over a Unix-domain socket:

    recv "<time> <albedo_coeff> <pboundary>\n"   ->   send "<albedo>\n"

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

MODES = ("echo", "pid-legacy", "pid-tuned", "pi-only", "full")


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


def build_adaptive(args, backstepping):
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
        regressor=os.environ.get("BS_REGRESSOR", "direct"),
        # "hat" (default) is the as-shipped n_basis-dimensional partition of
        # unity; "physics" is the single matched-shape diagnostic basis (see
        # `observers/passive.py`, `scripts/bc8_physics_basis.py`).
        basis=basis_env)
    return ctrl


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
    ctrl = None
    phi_scale = None
    p0 = None
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

    for line in stream:
        parts = line.split()
        if len(parts) != 3:
            continue
        time_s, albedo_in, pboundary = (float(x) for x in parts)

        # ---- first sample: latch the operating point, build the controller
        if counters["steps"] == 0:
            p0 = pboundary
            phi_scale = p0
            with open(os.path.join(out, "handshake.json"), "w") as f:
                json.dump({"mode": args.mode, "p0": p0, "phi_scale": phi_scale,
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
            elif args.mode in ("pi-only", "full"):
                ctrl = build_adaptive(args, backstepping=(args.mode == "full"))
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
        if args.mode == "echo":
            alb = albedo_in
        elif pid is not None:
            alb = pid.step(time_s, albedo_in, pboundary)
        else:
            if counters["steps"] == 0:
                alb = albedo_in              # same uncontrolled first step as PID
            else:
                ctrl.step_collocated(pboundary / phi_scale)
                alb = ctrl.albedo_cmd

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
                regressor = ident.xi
                xi_trace.append([time_s] + list(ident.xi))
            else:
                overlap = ident.Psi[:, :ident.grid.N] @ (
                    ident.phihat * ident.grid.wq[:ident.grid.N])
                regressor = overlap / ident.psi_mass
            regressor_trace.append([time_s] + list(regressor))
        if ctrl is not None:
            z_trace.append(ctrl.identifier.z_last)

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
                 z_trace=np.asarray(z_trace, dtype=float))
    print("controller_server: mode=%s steps=%d %s"
          % (args.mode, counters["steps"],
             {k: v for k, v in counters.items() if k != "steps" and v}),
          file=sys.stderr)


if __name__ == "__main__":
    main()
