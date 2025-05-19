#include "cell_interactions.h"
#include "physics.h"                 // For predict_all_events_for_particle
#include "../core/event_system/hybrid.h" // For schedule_event
#include "../core/particle_pool.h"      // For allocate_event_from_pool
#include "../utils/cell_list.h"       // For update_particle_cell_membership
#include "../utils/pbc.h"            // For pbc_wrap_position
#include <float.h>                   // For DBL_MAX
#include <math.h>                    // For fabs, fmin
#include <stdio.h>                   // For error logging

/**
 * @brief Predicts cell crossing events for a particle
 * 
 * Calculates when a particle will cross each boundary of its current cell
 * and schedules the corresponding events.
 * 
 * @param ctx Simulation context
 * @param p Particle to predict cell crossings for
 */
void predict_cell_crossing(SimContext* ctx, Particle* p) {
    if (!ctx || !p) return;
    
    double event_abs_time;

    // Calculate x-boundary crossing times
    double tx_rel = DBL_MAX; // Relative time to crossing
    double x_bound;
    EventType x_event_type = EVENT_NONE;

    if (fabs(p->vx) > 1e-12) {
        if (p->vx > 0) {
            // Moving right - will cross right boundary of current cell
            x_bound = (p->cellx + 1) * ctx->cell_size;
            // Ensure particle is not already past or exactly at the boundary in the direction of travel
            if (p->x < x_bound - 1e-12) { 
                tx_rel = (x_bound - p->x) / p->vx;
                x_event_type = EVENT_CELL_CROSS_X_POS;
            }
        } else { // p->vx < 0
            // Moving left - will cross left boundary of current cell
            x_bound = p->cellx * ctx->cell_size;
            // Ensure particle is not already past or exactly at the boundary in the direction of travel
            if (p->x > x_bound + 1e-12) {
                tx_rel = (x_bound - p->x) / p->vx; // tx_rel will be positive
                x_event_type = EVENT_CELL_CROSS_X_NEG;
            }
        }

        if (tx_rel > 0 && tx_rel != DBL_MAX && x_event_type != EVENT_NONE) {
            event_abs_time = ctx->current_time + tx_rel;
            if (event_abs_time < ctx->current_time + 1e-9) event_abs_time = ctx->current_time + 1e-9;
            
            Event* ev = allocate_event_from_pool(ctx);
            if (ev) {
                ev->time = event_abs_time;
                ev->type = x_event_type;
                ev->p1_idx = p->id;
                ev->p2_idx = -1; // Not a pair event
                ev->cross_dir = (p->vx > 0) ? 0 : 1; // 0 for +x, 1 for -x (example)
                schedule_event(ctx, ev);
            }
        }
    } else { 
        // Velocity is near zero, effectively no crossing in x planned
        tx_rel = DBL_MAX;
    }

    // Calculate y-boundary crossing times
    double ty_rel = DBL_MAX; // Relative time to crossing
    double y_bound;
    EventType y_event_type = EVENT_NONE;

    if (fabs(p->vy) > 1e-12) {
        if (p->vy > 0) {
            y_bound = (p->celly + 1) * ctx->cell_size;
            if (p->y < y_bound - 1e-12) {
                ty_rel = (y_bound - p->y) / p->vy;
                y_event_type = EVENT_CELL_CROSS_Y_POS;
            }
        } else { // p->vy < 0
            y_bound = p->celly * ctx->cell_size;
            if (p->y > y_bound + 1e-12) {
                ty_rel = (y_bound - p->y) / p->vy;
                y_event_type = EVENT_CELL_CROSS_Y_NEG;
            }
        }

        if (ty_rel > 0 && ty_rel != DBL_MAX && y_event_type != EVENT_NONE) {
            event_abs_time = ctx->current_time + ty_rel;
            if (event_abs_time < ctx->current_time + 1e-9) event_abs_time = ctx->current_time + 1e-9;

            Event* ev = allocate_event_from_pool(ctx);
            if (ev) {
                ev->time = event_abs_time;
                ev->type = y_event_type;
                ev->p1_idx = p->id;
                ev->p2_idx = -1; // Not a pair event
                ev->cross_dir = (p->vy > 0) ? 2 : 3; // 2 for +y, 3 for -y (example)
                schedule_event(ctx, ev);
            }
        }
    } else {
        ty_rel = DBL_MAX;
    }
}

/**
 * @brief Handles a cell crossing event in the positive X direction
 * 
 * @param ctx Simulation context
 * @param p Particle crossing the cell boundary
 */
void handle_cell_cross_x_pos(SimContext* ctx, Particle* p) {
    if (!ctx || !p) {
        if (ctx && ctx->hooks.log_error) ctx->hooks.log_error("Null context or particle in handle_cell_cross_x_pos");
        return;
    }

    // Expected boundary particle is crossing
    double expected_boundary = (double)(p->cellx + 1) * ctx->cell_size;

    // Correct position to be slightly past the boundary
    // This ensures it's in the new cell for subsequent calculations.
    p->x = expected_boundary + 1e-9; // Small epsilon past boundary

    // Apply PBC wrapping to ensure particle position is within the primary simulation box
    // This is crucial if the cell boundary is also a box boundary.
    pbc_wrap_position(&p->x, &p->y, ctx->xsize, ctx->ysize);


    int old_cellx = p->cellx;
    int old_celly = p->celly;

    // Update cell index
    p->cellx = (p->cellx + 1);
    if (p->cellx >= ctx->n_cells_x) {
        p->cellx = 0; 
        // Position already wrapped by pbc_wrap_position if it crossed box boundary
    }
    
    // Update particle's cell membership in the cell list data structure
    update_particle_cell_membership(ctx, p, old_cellx, old_celly);

    // Invalidate old events and predict new ones for this particle
    invalidate_events_for_particle(ctx, p->id);
    predict_all_events_for_particle(ctx, p);
}

/**
 * @brief Handles a cell crossing event in the negative X direction
 * 
 * @param ctx Simulation context
 * @param p Particle crossing the cell boundary
 */
void handle_cell_cross_x_neg(SimContext* ctx, Particle* p) {
    if (!ctx || !p) {
        if (ctx && ctx->hooks.log_error) ctx->hooks.log_error("Null context or particle in handle_cell_cross_x_neg");
        return;
    }
    
    double expected_boundary = (double)p->cellx * ctx->cell_size;
    p->x = expected_boundary - 1e-9; // Small epsilon before boundary

    pbc_wrap_position(&p->x, &p->y, ctx->xsize, ctx->ysize);

    int old_cellx = p->cellx;
    int old_celly = p->celly;

    p->cellx = (p->cellx - 1);
    if (p->cellx < 0) {
        p->cellx = ctx->n_cells_x - 1;
        // Position already wrapped
    }

    update_particle_cell_membership(ctx, p, old_cellx, old_celly);
    invalidate_events_for_particle(ctx, p->id);
    predict_all_events_for_particle(ctx, p);
}

/**
 * @brief Handles a cell crossing event in the positive Y direction
 * 
 * @param ctx Simulation context
 * @param p Particle crossing the cell boundary
 */
void handle_cell_cross_y_pos(SimContext* ctx, Particle* p) {
    if (!ctx || !p) {
        if (ctx && ctx->hooks.log_error) ctx->hooks.log_error("Null context or particle in handle_cell_cross_y_pos");
        return;
    }

    double expected_boundary = (double)(p->celly + 1) * ctx->cell_size;
    p->y = expected_boundary + 1e-9;

    pbc_wrap_position(&p->x, &p->y, ctx->xsize, ctx->ysize);
    
    int old_cellx = p->cellx;
    int old_celly = p->celly;

    p->celly = (p->celly + 1);
    if (p->celly >= ctx->n_cells_y) {
        p->celly = 0;
        // Position already wrapped
    }

    update_particle_cell_membership(ctx, p, old_cellx, old_celly);
    invalidate_events_for_particle(ctx, p->id);
    predict_all_events_for_particle(ctx, p);
}

/**
 * @brief Handles a cell crossing event in the negative Y direction
 * 
 * @param ctx Simulation context
 * @param p Particle crossing the cell boundary
 */
void handle_cell_cross_y_neg(SimContext* ctx, Particle* p) {
    if (!ctx || !p) {
        if (ctx && ctx->hooks.log_error) ctx->hooks.log_error("Null context or particle in handle_cell_cross_y_neg");
        return;
    }

    double expected_boundary = (double)p->celly * ctx->cell_size;
    p->y = expected_boundary - 1e-9;

    pbc_wrap_position(&p->x, &p->y, ctx->xsize, ctx->ysize);

    int old_cellx = p->cellx;
    int old_celly = p->celly;

    p->celly = (p->celly - 1);
    if (p->celly < 0) {
        p->celly = ctx->n_cells_y - 1;
        // Position already wrapped
    }
    
    update_particle_cell_membership(ctx, p, old_cellx, old_celly);
    invalidate_events_for_particle(ctx, p->id);
    predict_all_events_for_particle(ctx, p);
}

/**
 * @brief Generic handler for cell crossing events
 * 
 * Dispatches to the specific handler based on event type
 * 
 * @param ctx Simulation context
 * @param ev Cell crossing event
 */
void handle_cell_crossing(SimContext* ctx, Event* ev) {
    if (!ctx || !ev) return;
    
    Particle* p = &ctx->particles[ev->p1_idx];
    if (!p) return;
    
    switch (ev->type) {
        case EVENT_CELL_CROSS_X_POS:
            handle_cell_cross_x_pos(ctx, p);
            break;
        case EVENT_CELL_CROSS_X_NEG:
            handle_cell_cross_x_neg(ctx, p);
            break;
        case EVENT_CELL_CROSS_Y_POS:
            handle_cell_cross_y_pos(ctx, p);
            break;
        case EVENT_CELL_CROSS_Y_NEG:
            handle_cell_cross_y_neg(ctx, p);
            break;
        default:
            fprintf(stderr, "Error: Invalid cell crossing event type %d\n", ev->type);
            break;
    }
}
