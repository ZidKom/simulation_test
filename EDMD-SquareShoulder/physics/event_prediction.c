#include "event_prediction.h" // Should be specific header if exists, or general physics.h
#include "physics.h"          // For declarations like calculate_collision_time, solve_collision_time, NEVER
#include "../core/sim_context.h"
#include "../core/event_system/hybrid.h" // For schedule_event
#include "../core/particle_pool.h"     // For allocate_event_from_pool
#include "cell_interactions.h" // For predict_cell_crossing
#include <math.h>             // For DBL_MAX, fabs
#include <assert.h>
#include <stdio.h>            // For debugging, if any

// NOTE: Solver functions like solve_collision_time, calculate_shoulder_entry_time,
// calculate_shoulder_exit_time are defined in solve_physics.c and declared in physics.h.
// This file (event_prediction.c) should *call* them, not define them.

// Helper function to predict and schedule events between two particles
void predict_pair_events(SimContext* ctx, Particle* p1, Particle* p2) {
    if (!ctx || !p1 || !p2 || p1->id == p2->id) return;

    // Standard collision (core collision)
    double t_coll_abs = calculate_collision_time(ctx, p1, p2);
    if (t_coll_abs != DBL_MAX && t_coll_abs > ctx->current_time + 1e-12) {
        Event* ev = allocate_event_from_pool(ctx);
        if (!ev) { 
            if(ctx->hooks.log_error) ctx->hooks.log_error("Failed to allocate event for PP collision p%d-p%d", p1->id, p2->id);
            return; 
        }
        ev->type = EVENT_COLLISION;
        ev->time = t_coll_abs;
        ev->p1_idx = p1->id;
        ev->p2_idx = p2->id;
        schedule_event(ctx, ev);
    }

    // Shoulder entry
    double t_shoulder_entry_abs = calculate_shoulder_entry_time(ctx, p1, p2);
    if (t_shoulder_entry_abs != DBL_MAX && t_shoulder_entry_abs > ctx->current_time + 1e-12) {
        Event* ev = allocate_event_from_pool(ctx);
        if (!ev) { 
            if(ctx->hooks.log_error) ctx->hooks.log_error("Failed to allocate event for shoulder entry p%d-p%d", p1->id, p2->id);
            return; 
        }
        ev->type = EVENT_SHOULDER_ENTRY;
        ev->time = t_shoulder_entry_abs;
        ev->p1_idx = p1->id;
        ev->p2_idx = p2->id;
        schedule_event(ctx, ev);
    }
    
    // Shoulder exit
    // Only predict exit if they are currently in a shoulder interaction together
    if (p1->in_shoulder && p2->in_shoulder && 
        ((p1->shoulder_partner_idx == p2->id) && (p2->shoulder_partner_idx == p1->id)) ) {
        double t_shoulder_exit_abs = calculate_shoulder_exit_time(ctx, p1, p2);
        if (t_shoulder_exit_abs != DBL_MAX && t_shoulder_exit_abs > ctx->current_time + 1e-12) {
            Event* ev = allocate_event_from_pool(ctx);
            if (!ev) { 
                if(ctx->hooks.log_error) ctx->hooks.log_error("Failed to allocate event for shoulder exit p%d-p%d", p1->id, p2->id);
                return; 
            }
            ev->type = EVENT_SHOULDER_EXIT;
            ev->time = t_shoulder_exit_abs;
            ev->p1_idx = p1->id;
            ev->p2_idx = p2->id;
            schedule_event(ctx, ev);
        }
    }
}

void predict_all_events_for_particle(SimContext* ctx, Particle* p) {
    assert(ctx && p);

    // Get current cell of particle p
    // Assuming p->cellx and p->celly are up-to-date
    int current_cell_x = p->cellx;
    int current_cell_y = p->celly;

    // Iterate over 3x3 neighborhood of cells (including current cell)
    for (int dx_cell = -1; dx_cell <= 1; dx_cell++) {
        for (int dy_cell = -1; dy_cell <= 1; dy_cell++) {
            int neighbor_cell_x = (current_cell_x + dx_cell + ctx->n_cells_x) % ctx->n_cells_x;
            int neighbor_cell_y = (current_cell_y + dy_cell + ctx->n_cells_y) % ctx->n_cells_y;

            // Iterate through particles in this neighbor cell
            Particle* neighbor = ctx->cells[neighbor_cell_x][neighbor_cell_y].head;
            while (neighbor) {
                if (p->id != neighbor->id) { // Don't predict with self
                    // Predict pair events (collision, shoulder entry, shoulder exit)
                    predict_pair_events(ctx, p, neighbor);
                }
                neighbor = neighbor->next_in_cell;
            }
        }
    }

    // Predict cell crossings for this particle p
    predict_cell_crossing(ctx, p);
}