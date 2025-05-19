#ifndef PBC_H
#define PBC_H

#include <math.h>

/**
 * Apply periodic boundary condition wrapping to a coordinate difference.
 * @param x Difference in coordinate.
 * @param box Box length.
 * @return The wrapped coordinate difference.
 */
double pbc(double x, double box);
#define PBC pbc  // Create a macro for backward compatibility

/**
 * Same as pbc but with a clearer name for validation code.
 * @param x Difference in coordinate.
 * @param box Box length.
 * @return The minimum distance considering periodic boundaries.
 */
static inline double pbc_distance(double x, double box) {
    return pbc(x, box);
}

/**
 * Wrap a position (x,y) into the periodic box.
 * @param x Pointer to x coordinate.
 * @param y Pointer to y coordinate.
 * @param boxx Box length in x.
 * @param boxy Box length in y.
 */
void pbc_wrap_position(double* x, double* y, double boxx, double boxy);

/**
 * Calculate minimum image distance between two points.
 * @param v1 First coordinate.
 * @param v2 Second coordinate.
 * @param L Box length.
 * @return Minimum image displacement.
 */
double pbc_min_image_delta(double v1, double v2, double L);

/**
 * Minimum image convention for a coordinate difference (dx,dy).
 * @param dx Pointer to x difference.
 * @param dy Pointer to y difference.
 * @param boxx Box length in x.
 * @param boxy Box length in y.
 */
void pbc_apply_min_image(double* dx, double* dy, double boxx, double boxy);

#endif // PBC_H
