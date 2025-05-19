#include "collision_validation.h"
#include "../utils/pbc.h"
#include <math.h>

#define EPS 1e-12  // Epsilon for numerical precision

/**
 * @brief Validates a collision event to ensure numerical precision.
 * 
 * This function checks if the particles are actually at the correct distance
 * for a collision, within a small epsilon tolerance. This helps detect numerical 
 * drift that can occur in simulations.
 * 
 * @param ctx Pointer to the simulation context
 * @param ev Pointer to the collision event to validate
 * @return 1 if the collision is valid, 0 otherwise
 */
int validate_collision(SimContext* ctx, Event* ev) {
    if (!ctx || !ev) return 0;
    
    // Only validate collision events
    if (ev->type != EVENT_COLLISION) return 1;
    
    // Get the particles involved
    Particle *p1 = &ctx->particles[ev->p1_idx];
    Particle *p2 = &ctx->particles[ev->p2_idx];
    
    // Calculate the current distance between particles
    double dx = pbc_min_image_delta(p1->x, p2->x, ctx->xsize);
    double dy = pbc_min_image_delta(p1->y, p2->y, ctx->ysize);
    double actual_dist = sqrt(dx*dx + dy*dy);
    
    // Expected distance at collision is the sum of radii
    double expected_dist = p1->radius + p2->radius;
    
    // 0.01% tolerance for numerical precision
    if (fabs(actual_dist - expected_dist)/expected_dist > 1e-4) {
        // Log the error if a logging hook is available
        if (ctx->hooks.log_error) {
            ctx->hooks.log_error("Collision distance violation: %.12f != %.12f", actual_dist, expected_dist);
        }
        return 0;  // Invalid collision
    }
    
    return 1;  // Valid collision
}

// This file should contain the only definition of validate_collision.
