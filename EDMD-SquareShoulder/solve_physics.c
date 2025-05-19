#include "physics/physics.h" // Include the header that declares these functions
#include "core/sim_context.h"
#include "utils/pbc.h"
#include "physics/cell_interactions.h"
#include <math.h>
#include <float.h>
#include <stdio.h>
#include <assert.h>

// This file contains implementations of the core physics functions
// to resolve the linker errors

/**
 * @brief Numerically stable quadratic solver for event time calculation
 */
double solve_collision_time(double dx, double dy, double dvx, double dvy, double sigma_sum) {
    double a = dvx*dvx + dvy*dvy;
    double b = 2*(dx*dvx + dy*dvy);
    double c = dx*dx + dy*dy - sigma_sum*sigma_sum;
    double discriminant = b*b - 4*a*c;
    if (discriminant < 0) return DBL_MAX;
    double sqrt_disc = sqrt(discriminant);
    double t1 = (-b - sqrt_disc) / (2*a);
    double t2 = (-b + sqrt_disc) / (2*a);
    if (t1 > 1e-12) return t1;
    if (t2 > 1e-12) return t2;
    return DBL_MAX;
}

/**
 * @brief Calculate time for particles to reach shoulder entry distance
 */
double calculate_shoulder_entry_time(SimContext *ctx, Particle *p1, Particle *p2) {
    if (!ctx || !p1 || !p2) return DBL_MAX;
    // Similar to calculate_collision_time, but with shoulder_radius
    double shoulder_radius = ctx->params.sigma * ctx->params.lambda_shoulder;
    
    double dx = pbc_min_image_delta(p1->x, p2->x, ctx->xsize);
    double dy = pbc_min_image_delta(p1->y, p2->y, ctx->ysize);
    double dvx = p2->vx - p1->vx; // Relative velocity p2 relative to p1
    double dvy = p2->vy - p1->vy;

    // Time to reach shoulder_radius
    // We are looking for entry, so particles should be approaching this radius.
    // If currently dist > shoulder_radius, they must be approaching (dr.dv < 0).
    // If currently dist < shoulder_radius (i.e. inside core), this is not an entry event from outside.
    // This function should calculate time to reach shoulder_radius *from outside*.
    
    double dist_sq = dx*dx + dy*dy;
    if (dist_sq < shoulder_radius * shoulder_radius - 1e-9) { // Already inside or at shoulder
        return DBL_MAX; 
    }

    double t = solve_collision_time(dx, dy, dvx, dvy, shoulder_radius);
    // solve_collision_time returns smallest positive t for reaching separation 'sigma_sum'
    // We need to ensure they are approaching: b = dx*dvx + dy*dvy should be < 0.
    // solve_collision_time itself should handle this by returning DBL_MAX if no approach.
    
    return (t < DBL_MAX && t > 1e-9) ? (ctx->current_time + t) : DBL_MAX;
}

/**
 * @brief Calculate time for particles to reach shoulder exit distance
 */
double calculate_shoulder_exit_time(SimContext *ctx, Particle *p1, Particle *p2) {
    if (!ctx || !p1 || !p2) return DBL_MAX;
    // This function calculates time to reach shoulder_radius *from inside*, moving apart.
    double shoulder_radius = ctx->params.sigma * ctx->params.lambda_shoulder;
    
    double dx = pbc_min_image_delta(p1->x, p2->x, ctx->xsize);
    double dy = pbc_min_image_delta(p1->y, p2->y, ctx->ysize);
    double dvx = p2->vx - p1->vx; 
    double dvy = p2->vy - p1->vy;

    double dist_sq = dx*dx + dy*dy;
    // Must be inside shoulder (dist < shoulder_radius) and moving apart (dr.dv > 0)
    if (dist_sq >= shoulder_radius * shoulder_radius - 1e-9) { // Already outside or at shoulder boundary
        return DBL_MAX;
    }
    if (dx*dvx + dy*dvy <= 0) { // Not moving apart
        return DBL_MAX;
    }

    // Use solve_collision_time to find time to reach shoulder_radius.
    // It finds the smallest positive t. If they are inside and moving apart,
    // one root of quadratic for reaching distance R will be positive.
    double t = solve_collision_time(dx, dy, dvx, dvy, shoulder_radius);

    return (t < DBL_MAX && t > 1e-9) ? (ctx->current_time + t) : DBL_MAX;
}

/**
 * @brief Predict cell crossing events
 */
void predict_cell_crossing_events(SimContext* ctx) {
    for (int i = 0; i < ctx->num_particles; ++i) {
        predict_cell_crossing(ctx, &ctx->particles[i]);
    }
}

#if 0
/**
 * @brief Predict all events for a given particle
 */
void predict_all_events_for_particle(SimContext *ctx, Particle *p) {
    // Implementation specific to the simulation
}
#endif

// Only the solver functions should be defined here.
