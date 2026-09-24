#ifndef TRANSIENT_ALBEDO_H
#define TRANSIENT_ALBEDO_H

double max(const double a, const double b);
double min(const double a, const double b);

double transient_albedo(const double time, const double albedo_coeff,
    const double pow_lhs, const double pow_rhs);

#endif
