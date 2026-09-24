/* Boundary-control bridge: forwards the albedo actuation decision to an
 * external controller process over a Unix-domain socket.
 *
 * Protocol (newline-delimited text, one exchange per timestep):
 *   send    "<time> <albedo_coeff> <pboundary>\n"   with %.17g fields
 *   receive "<albedo>\n"
 *
 * %.17g round-trips an IEEE double exactly; anything coarser makes last-bit
 * sensor noise indistinguishable from a control bug.
 *
 * Paths are taken from the environment so several runs can share a working
 * tree without colliding:
 *   SIREN_CTRL_SOCK  socket to connect to    (default "controller.sock")
 *   SIREN_CTRL_LOG   per-step CSV trace      (default "controller_log.csv")
 *
 * pboundary above is pow_rhs, the power at the actuated face x=L.
 */

/* snprintf and the POSIX socket declarations are not visible under -ansi
 * unless a feature macro is defined before the first include. */
#define _POSIX_C_SOURCE 200809L

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/socket.h>
#include <sys/un.h>
#include "transient_albedo.h"

#define SOCK_PATH_DEFAULT "controller.sock"
#define LOG_PATH_DEFAULT  "controller_log.csv"

static FILE *log_file = NULL;
static int sock_fd = -1;
static int first = 1;
static double log_p0;

double max(const double a, const double b){ return (a > b) ? a : b; }
double min(const double a, const double b){ return (a < b) ? a : b; }

double transient_albedo(const double time, const double albedo_coeff, 
    const double pow_lhs, const double pow_rhs)
{
  char sbuf[128];
  char rbuf[64];
  int len;
  double alb;
  double error;
  /* the controlled output is the actuated face, x=L; pow_lhs is logged for
   * the left-hand-side observer but not sent to the controller */
  const double pboundary = pow_rhs;

  if (first)
  {
    struct sockaddr_un addr;
    const char *sock_path;
    const char *log_path;

    sock_path = getenv("SIREN_CTRL_SOCK");
    if (sock_path == NULL || sock_path[0] == '\0')
      sock_path = SOCK_PATH_DEFAULT;

    log_path = getenv("SIREN_CTRL_LOG");
    if (log_path == NULL || log_path[0] == '\0')
      log_path = LOG_PATH_DEFAULT;

    sock_fd = socket(AF_UNIX, SOCK_STREAM, 0);
    if (sock_fd < 0)
    {
      perror("transient_albedo: socket");
      exit(1);
    }

    memset(&addr, 0, sizeof(addr));
    addr.sun_family = AF_UNIX;
    strncpy(addr.sun_path, sock_path, sizeof(addr.sun_path) - 1);

    if (connect(sock_fd, (struct sockaddr *)&addr, sizeof(addr)) < 0)
    {
      fprintf(stderr, "transient_albedo: connect to '%s': ", sock_path);
      perror(NULL);
      exit(1);
    }

    log_p0 = pboundary;
    error = 0.0;
    first = 0;

    log_file = fopen(log_path, "w");
    if (log_file == NULL)
    {
      fprintf(stderr, "transient_albedo: cannot open log '%s': ", log_path);
      perror(NULL);
      exit(1);
    }
    fprintf(log_file, "time,pboundary,albedo,error,pow_lhs\n");
  }
  else
  {
    error = log_p0 - pboundary;
  }

  len = snprintf(sbuf, sizeof(sbuf), "%.17g %.17g %.17g\n",
                 time, albedo_coeff, pboundary);
  if (write(sock_fd, sbuf, (size_t)len) != (ssize_t)len)
  {
    perror("transient_albedo: write");
    exit(1);
  }

  len = 0;
  while (len < (int)sizeof(rbuf) - 1)
  {
    char c;
    if (read(sock_fd, &c, 1) != 1)
    {
      fprintf(stderr, "transient_albedo: read failed or connection closed\n");
      exit(1);
    }
    if (c == '\n')
      break;
    rbuf[len++] = c;
  }
  rbuf[len] = '\0';

  alb = strtod(rbuf, NULL);
  /* clamp: safety net regardless of what the server returns */
  alb = min(max(alb, -1.0), 1.0);

  /* %.17e, not %.6e: at pboundary ~ 2e-2 a 6-decimal field quantizes to
   * 1e-8 absolute, i.e. 5e-7 relative -- coarser than the boundary
   * regulation error the controlled runs achieve, so M1 below that was
   * reading the log format rather than the plant. */
  fprintf(log_file, "%.17e,%.17e,%.17e,%.17e,%.17e\n",
          time, pboundary, alb, error, pow_lhs);
  fflush(log_file);

  return alb;
}
