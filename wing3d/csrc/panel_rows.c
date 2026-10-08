/* Fused panel-method kernels (Morino Dirichlet + Kutta).
 *
 * Mirrors wing3d/solver.py numpy semantics bit-closely (same formulas,
 * same clamps/cutoffs). Self-terms are left for the Python caller
 * (analytic overrides), so these compute raw far+near sums for ALL pairs.
 * Pure C99, no dependencies. Layouts are C-contiguous row-major.
 */
#include <math.h>

#define PI 3.14159265358979323846

/* Signed solid angle of fan triangles at target points.
 * P[np,3], fan[n,F,3,3], fmask[n,F] -> om[np,n] (summed, unscaled). */
void solid_angle_block(const double *P, int np_, const double *fan,
                       const unsigned char *fmask, int n, int F, double *om) {
    for (int t = 0; t < np_; t++) {
        double px = P[3 * t], py = P[3 * t + 1], pz = P[3 * t + 2];
        for (int j = 0; j < n; j++) {
            double sum = 0.0;
            for (int f = 0; f < F; f++) {
                if (!fmask[j * F + f]) continue;
                const double *T = fan + ((j * F + f) * 3 * 3);
                double ax = T[0] - px, ay = T[1] - py, az = T[2] - pz;
                double bx = T[3] - px, by = T[4] - py, bz = T[5] - pz;
                double cx = T[6] - px, cy = T[7] - py, cz = T[8] - pz;
                double la = sqrt(ax * ax + ay * ay + az * az);
                double lb = sqrt(bx * bx + by * by + bz * bz);
                double lc = sqrt(cx * cx + cy * cy + cz * cz);
                /* b x c */
                double qx = by * cz - bz * cy;
                double qy = bz * cx - bx * cz;
                double qz = bx * cy - by * cx;
                double num = ax * qx + ay * qy + az * qz;
                double ab = ax * bx + ay * by + az * bz;
                double bc = bx * cx + by * cy + bz * cz;
                double ca = cx * ax + cy * ay + cz * az;
                double den = la * lb * lc + ab * lc + bc * la + ca * lb;
                sum += 2.0 * atan2(num, den);
            }
            om[t * n + j] = sum;
        }
    }
}

/* (1/4pi) int 1/r dS per panel. Far: monopole; near: quadrature.
 * cent[n,3] (hoisted centroids), area[n], size[n]. */
void source_pot_block(const double *P, int np_, const double *qp,
                      const double *qw, int n, int Q, const double *area,
                      const double *size, const double *cent, double *out) {
    const double fpi = 4.0 * PI;
    for (int t = 0; t < np_; t++) {
        double px = P[3 * t], py = P[3 * t + 1], pz = P[3 * t + 2];
        for (int j = 0; j < n; j++) {
            double dx = cent[3 * j] - px;
            double dy = cent[3 * j + 1] - py;
            double dz = cent[3 * j + 2] - pz;
            double r = sqrt(dx * dx + dy * dy + dz * dz);
            double v;
            if (r > 4.0 * size[j]) {
                v = area[j] / (fpi * r);
            } else {
                double s = 0.0;
                for (int q = 0; q < Q; q++) {
                    double qx = qp[(j * Q + q) * 3] - px;
                    double qy = qp[(j * Q + q) * 3 + 1] - py;
                    double qz = qp[(j * Q + q) * 3 + 2] - pz;
                    double d = sqrt(qx * qx + qy * qy + qz * qz);
                    if (d < 1e-12) d = 1e-12;
                    s += qw[j * Q + q] / d;
                }
                v = s / fpi;
            }
            out[t * n + j] = v;
        }
    }
}

/* -(1/4pi) int r_vec/r^3 dS per panel -> out[t*n+j, 3].
 * Far uses +area*rv/4pi d^3 with rv = p - c (matches numpy branch). */
void source_vel_block(const double *P, int np_, const double *qp,
                      const double *qw, int n, int Q, const double *area,
                      const double *size, const double *cent, double *out) {
    const double fpi = 4.0 * PI;
    for (int t = 0; t < np_; t++) {
        double px = P[3 * t], py = P[3 * t + 1], pz = P[3 * t + 2];
        for (int j = 0; j < n; j++) {
            double dx = cent[3 * j] - px;
            double dy = cent[3 * j + 1] - py;
            double dz = cent[3 * j + 2] - pz;
            double r = sqrt(dx * dx + dy * dy + dz * dz);
            double ox = 0.0, oy = 0.0, oz = 0.0;
            if (r > 4.0 * size[j]) {
                double d = r < 1e-14 ? 1e-14 : r;
                double s = area[j] / (fpi * d * d * d);
                /* rv = p - c */
                ox = s * (-dx);
                oy = s * (-dy);
                oz = s * (-dz);
            } else {
                double sx = 0.0, sy = 0.0, sz = 0.0;
                for (int q = 0; q < Q; q++) {
                    /* d = qp - p ; summand (qw/dist^3)*d, negated after */
                    double qx = qp[(j * Q + q) * 3] - px;
                    double qy = qp[(j * Q + q) * 3 + 1] - py;
                    double qz = qp[(j * Q + q) * 3 + 2] - pz;
                    double dist = sqrt(qx * qx + qy * qy + qz * qz);
                    if (dist < 1e-12) dist = 1e-12;
                    double w = qw[j * Q + q] /
                               (dist * dist * dist);
                    sx += w * qx;
                    sy += w * qy;
                    sz += w * qz;
                }
                ox = -sx / fpi;
                oy = -sy / fpi;
                oz = -sz / fpi;
            }
            out[(t * n + j) * 3] = ox;
            out[(t * n + j) * 3 + 1] = oy;
            out[(t * n + j) * 3 + 2] = oz;
        }
    }
}

/* Vortex-loop (Biot-Savart) velocity per panel, unit strength.
 * ed[n,E,2,3], emask[n,E]; cutoff 1e-10, d clamp 1e-14. */
void loop_vel_block(const double *P, int np_, const double *ed,
                    const unsigned char *emask, int n, int E, double *out) {
    const double fpi = 4.0 * PI;
    for (int t = 0; t < np_; t++) {
        double px = P[3 * t], py = P[3 * t + 1], pz = P[3 * t + 2];
        for (int j = 0; j < n; j++) {
            double sx = 0.0, sy = 0.0, sz = 0.0;
            for (int e = 0; e < E; e++) {
                if (!emask[j * E + e]) continue;
                const double *E0 = ed + ((j * E + e) * 2 * 3);
                const double *E1 = E0 + 3;
                double r0x = E1[0] - E0[0];
                double r0y = E1[1] - E0[1];
                double r0z = E1[2] - E0[2];
                double r1x = px - E0[0], r1y = py - E0[1], r1z = pz - E0[2];
                double r2x = px - E1[0], r2y = py - E1[1], r2z = pz - E1[2];
                double crx = r1y * r2z - r1z * r2y;
                double cry = r1z * r2x - r1x * r2z;
                double crz = r1x * r2y - r1y * r2x;
                double denom = crx * crx + cry * cry + crz * crz + 1e-10;
                double d1 = sqrt(r1x * r1x + r1y * r1y + r1z * r1z);
                double d2 = sqrt(r2x * r2x + r2y * r2y + r2z * r2z);
                if (d1 < 1e-14) d1 = 1e-14;
                if (d2 < 1e-14) d2 = 1e-14;
                double cos_term =
                    r0x * (r1x / d1 - r2x / d2) +
                    r0y * (r1y / d1 - r2y / d2) +
                    r0z * (r1z / d1 - r2z / d2);
                double s = cos_term / denom / fpi;
                sx += crx * s;
                sy += cry * s;
                sz += crz * s;
            }
            out[(t * n + j) * 3] = sx;
            out[(t * n + j) * 3 + 1] = sy;
            out[(t * n + j) * 3 + 2] = sz;
        }
    }
}
