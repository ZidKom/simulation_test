#include "hooks.h"
#include "energy_audit.h"
#include "../core/sim_context.h"
#include "../utils/pbc.h" // Added for periodic boundary conditions
#include <math.h> // For M_PI and pow
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>

/**
 * @brief Default error logging function
 * 
 * @param format Format string for the error message
 * @param ... Variable arguments for the format string
 */
void default_log_error(const char* format, ...) {
    va_list args;
    va_start(args, format);
    fprintf(stderr, "[ERROR] ");
    vfprintf(stderr, format, args);
    fprintf(stderr, "\n");
    va_end(args);
}

/**
 * @brief Validates collision distance between particles
 * 
 * @param ctx Simulation context
 * @param ev Collision event
 */
void validate_collision_distance(SimContext* ctx, Event* ev) {
    if (!ctx || !ev || ev->type != EVENT_COLLISION || 
        ev->p1_idx < 0 || ev->p1_idx >= ctx->num_particles ||
        ev->p2_idx < 0 || ev->p2_idx >= ctx->num_particles) {
        return;
    }
    
    Particle* p1 = &ctx->particles[ev->p1_idx];
    Particle* p2 = &ctx->particles[ev->p2_idx];
    
    // Apply periodic boundary conditions when calculating distance
    double dx = p2->x - p1->x;
    double dy = p2->y - p1->y;
    // Apply minimum image convention for periodic boundaries
    dx = pbc_min_image_delta(p1->x, p2->x, ctx->xsize);
    dy = pbc_min_image_delta(p1->y, p2->y, ctx->ysize);
    
    double dist = sqrt(dx*dx + dy*dy);
    double expected = p1->radius + p2->radius;
    
    if(fabs(dist - expected) > 1e-9) {
        if (ctx->hooks.log_error) {
            ctx->hooks.log_error("Collision distance violation: %.12f != %.12f", dist, expected);
        } else {
            fprintf(stderr, "[COLLISION] Distance violation: %.12f != %.12f\n", dist, expected);
        }
        
        #ifdef VALIDATION_MODE
        abort(); // Terminate on validation errors in validation mode
        #endif
    }
}

/**
 * @brief Validates cell indices for all particles
 * 
 * @param ctx Simulation context
 */
void validate_cell_indices(SimContext* ctx) {
    if (!ctx) return;
    
    for(int i=0; i < ctx->num_particles; i++) {
        Particle* p = &ctx->particles[i];
        if (p->cellx < 0 || p->cellx >= ctx->n_cells_x || 
            p->celly < 0 || p->celly >= ctx->n_cells_y) {
            
            if (ctx->hooks.log_error) {
                ctx->hooks.log_error("Cell index violation: particle %d at cell (%d,%d), bounds [0,%d)x[0,%d)",
                    i, p->cellx, p->celly, ctx->n_cells_x, ctx->n_cells_y);
            } else {
                fprintf(stderr, "[CELL] Index violation: particle %d at cell (%d,%d), bounds [0,%d)x[0,%d)\n",
                    i, p->cellx, p->celly, ctx->n_cells_x, ctx->n_cells_y);
            }
            
            #ifdef VALIDATION_MODE
            abort(); // Terminate on validation errors in validation mode
            #endif
        }
    }
}

/**
 * @brief Default cell crossing validation
 * 
 * @param ctx Simulation context
 * @param ev Cell crossing event
 */
void validate_cell_crossing(SimContext* ctx, Event* ev) {
    (void)ev; // Unused parameter
    // Just validate cell indices after every crossing
    validate_cell_indices(ctx);
}

/**
 * @brief Creates a default set of validation hooks
 * 
 * @return SimulationHooks Default validation hooks
 */
SimulationHooks create_default_validation_hooks() {
    SimulationHooks hooks;
    
    hooks.log_warning = NULL;
    hooks.log_error = default_log_error;
    hooks.validate_collision = validate_collision_distance;
    
    return hooks;
}

/**
 * @brief Validates initial conditions of the simulation
 * 
 * @param ctx Simulation context
 */
void validate_initial_conditions(SimContext* ctx) {
    // Check box size matches particle positions
    for(int i=0; i < ctx->num_particles; i++) {
        if(i >= ctx->num_particles) {
            if(ctx->hooks.log_error)
                ctx->hooks.log_error("Particle index %d exceeds allocated count", i);
            break;
        }
        // Check if particles are within the primary box [0, xsize) and [0, ysize)
        // Note: Particles could be slightly outside due to PBC wrap-around logic
        // if they were loaded that way before wrapping.
        // A more robust check might consider particles within [-epsilon, xsize+epsilon]
        // For now, a strict check based on typical expectations after pbc_wrap_position.
        if(ctx->particles[i].x < 0.0 || ctx->particles[i].x >= ctx->xsize ||
           ctx->particles[i].y < 0.0 || ctx->particles[i].y >= ctx->ysize) {
            // Log an error if a particle is detected outside the defined simulation box boundaries.
            // This check is crucial for ensuring that all particles are correctly initialized within the simulation domain.
            ctx->hooks.log_error("Particle %d out of bounds: (%.2f,%.2f) in box (%.2f,%.2f)",
                               i, ctx->particles[i].x, ctx->particles[i].y,
                               ctx->xsize, ctx->ysize);
        }
    }
    
    // Check particle density
    if (ctx->xsize > 0 && ctx->ysize > 0 && ctx->num_particles > 0 && ctx->params.default_particle_radius > 0) {
        double area = ctx->xsize * ctx->ysize;
        double occupied_area_sum = 0.0;
        for(int i=0; i < ctx->num_particles; i++) {
            occupied_area_sum += M_PI * pow(ctx->particles[i].radius, 2);
        }
        double packing_fraction = occupied_area_sum / area;
        if(packing_fraction > 0.7) {
            if(ctx->hooks.log_error) {
                ctx->hooks.log_error("High packing fraction: %.3f", packing_fraction);
            }
        }
    }
}
