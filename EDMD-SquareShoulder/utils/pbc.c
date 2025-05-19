#include <math.h>
#include "pbc.h"

double pbc(double x, double box) {
    return x - box * round(x / box);
}

// Wrap a 2D position into the box
void pbc_wrap_position(double* x, double* y, double boxx, double boxy) {
    *x = *x - boxx * floor(*x / boxx);
    *y = *y - boxy * floor(*y / boxy);
}

// Implementation for new declaration
double pbc_min_image_delta(double v1, double v2, double L) {
    double dx = v2 - v1;
    return dx - L * round(dx / L);
}

// Rename the existing function to avoid conflicts
void pbc_apply_min_image(double* dx, double* dy, double boxx, double boxy) {
    *dx -= boxx * round(*dx / boxx);
    *dy -= boxy * round(*dy / boxy);
}
