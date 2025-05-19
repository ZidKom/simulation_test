#include "shoulders.h"
#include "../core/sim_context.h"
#include "physics.h" // For NEVER and event prediction functions
#include "event_prediction.h" // For predict_all_events_for_particle
#include "particle_interactions.h" // For tracking shoulder state and partner invalidation
#include "../core/event_system/hybrid.h" // For schedule_event
#include "../core/particle_pool.h" // For allocate_event_from_pool
#include "../utils/pbc.h" // For pbc_min_image_delta
#include <assert.h>
#include "../validation/energy_audit.h"

#include <math.h>
#include <stdio.h> // For error logging if needed
#include <float.h> // For DBL_MAX

/**
 * @brief Handles the consequences of two particles entering a square shoulder potential.
 *
 * When two particles enter a shoulder, their interaction potential changes.
 * This function is called when an EVENT_SHOULDER_ENTRY occurs.
 *
 * @param ctx Pointer to the simulation context.
 * @param event Pointer to the shoulder entry event.
 */
void handle_shoulder_entry(SimContext *ctx, Event *event) {
    if (!ctx || !event || event->type != EVENT_SHOULDER_ENTRY) {
        if (ctx && ctx->hooks.log_error) {
            ctx->hooks.log_error("Invalid call to handle_shoulder_entry: null ctx/event or wrong event type.");
        }
        return;
    }
    if(event->p1_idx < 0 || event->p1_idx >= ctx->num_particles ||
       event->p2_idx < 0 || event->p2_idx >= ctx->num_particles ||
       event->p1_idx == event->p2_idx) {
        if (ctx && ctx->hooks.log_error) {
            ctx->hooks.log_error("Invalid particle indices in shoulder entry event: p1=%d, p2=%d", event->p1_idx, event->p2_idx);
        }
        return;
    }
    Particle *p1 = &ctx->particles[event->p1_idx];
    Particle *p2 = &ctx->particles[event->p2_idx];
    double shoulder_radius = ctx->params.sigma * ctx->params.lambda_shoulder;

    // --- Position Correction ---
    double dx_current = pbc_min_image_delta(p1->x, p2->x, ctx->xsize);
    double dy_current = pbc_min_image_delta(p1->y, p2->y, ctx->ysize);
    double dist_current = sqrt(dx_current*dx_current + dy_current*dy_current);

    if (fabs(dist_current - shoulder_radius) > 1e-9) {
        if (dist_current > 1e-12) {
            double correction_dist_half = 0.5 * (dist_current - shoulder_radius);
            double u_nx = dx_current / dist_current;
            double u_ny = dy_current / dist_current;

            // Move toward the boundary (correct sign)
            p1->x -= correction_dist_half * u_nx;
            p1->y -= correction_dist_half * u_ny;
            p2->x += correction_dist_half * u_nx;
            p2->y += correction_dist_half * u_ny;

            pbc_wrap_position(&p1->x, &p1->y, ctx->xsize, ctx->ysize);
            pbc_wrap_position(&p2->x, &p2->y, ctx->xsize, ctx->ysize);
        }
    }
    // Debug output for post-correction distance
    dx_current = pbc_min_image_delta(p1->x, p2->x, ctx->xsize);
    dy_current = pbc_min_image_delta(p1->y, p2->y, ctx->ysize);
    dist_current = sqrt(dx_current*dx_current + dy_current*dy_current);
    if (fabs(dist_current - shoulder_radius) > 1e-8) {
        printf("[DEBUG] Shoulder ENTRY: Post-correction distance violation: %.12f (expected %.12f)\n", dist_current, shoulder_radius);
    }

    // --- Physics: KE/PE exchange ---
    double nx = dx_current / dist_current;
    double ny = dy_current / dist_current;
    double dvx = p2->vx - p1->vx;
    double dvy = p2->vy - p1->vy;
    double v_rel_n = dvx * nx + dvy * ny;
    double m1 = p1->m;
    double m2 = p2->m;
    double mu = (m1 * m2) / (m1 + m2);
    double KE_normal_relative = 0.5 * mu * v_rel_n * v_rel_n;

    if (v_rel_n < -1e-9) { // Approaching
        if (KE_normal_relative > ctx->params.U + 1e-12) {
            double v_prime_rel_n_sq = (KE_normal_relative - ctx->params.U) / (0.5 * mu);
            if (v_prime_rel_n_sq < 0) v_prime_rel_n_sq = 0; // Clamp for safety
            double v_prime_rel_n = -sqrt(v_prime_rel_n_sq); // Keep approaching direction
            double delta_v_rel_n = v_prime_rel_n - v_rel_n;

            p1->vx -= (mu / m1) * delta_v_rel_n * nx;
            p1->vy -= (mu / m1) * delta_v_rel_n * ny;
            p2->vx += (mu / m2) * delta_v_rel_n * nx;
            p2->vy += (mu / m2) * delta_v_rel_n * ny;

            p1->in_shoulder = 1; p1->shoulder_partner_idx = p2->id;
            p2->in_shoulder = 1; p2->shoulder_partner_idx = p1->id;
            if (ctx->hooks.validate_collision) track_energy_injection(ctx, -ctx->params.U);
        } else {
            // Reflect elastically
            double delta_v_rel_n = -2 * v_rel_n;

            p1->vx -= (mu / m1) * delta_v_rel_n * nx;
            p1->vy -= (mu / m1) * delta_v_rel_n * ny;
            p2->vx += (mu / m2) * delta_v_rel_n * nx;
            p2->vy += (mu / m2) * delta_v_rel_n * ny;

            p1->in_shoulder = 0; p1->shoulder_partner_idx = -1;
            p2->in_shoulder = 0; p2->shoulder_partner_idx = -1;
        }
    } else {
        // Not approaching: treat as stale event, just re-predict
        if (ctx && ctx->hooks.log_warning) {
            ctx->hooks.log_warning("Shoulder entry event for non-approaching particles (v_rel_n=%.2e). P%d-P%d @ t=%.6f. No velocity change.", v_rel_n, p1->id, p2->id, event->time);
        }
        invalidate_events_for_particle(ctx, p1->id);
        invalidate_events_for_particle(ctx, p2->id);
        predict_all_events_for_particle(ctx, p1);
        predict_all_events_for_particle(ctx, p2);
        return;
    }

    // Invalidate all events for both particles and partners
    invalidate_events_for_particle(ctx, p1->id);
    invalidate_events_for_particle(ctx, p2->id);
    invalidate_partner_events(ctx, p1->id);
    invalidate_partner_events(ctx, p2->id);
    predict_all_events_for_particle(ctx, p1);
    predict_all_events_for_particle(ctx, p2);
}

/**
 * @brief Handles the consequences of two particles exiting a square shoulder potential.
 *
 * When two particles exit a shoulder, their interaction potential changes back to normal.
 * This function is called when an EVENT_SHOULDER_EXIT occurs.
 *
 * @param ctx Pointer to the simulation context.
 * @param event Pointer to the shoulder exit event.
 */
void process_shoulder_exit(SimContext *ctx, Event *event) {
    if (!ctx || !event || event->type != EVENT_SHOULDER_EXIT) {
        if (ctx && ctx->hooks.log_error) {
            ctx->hooks.log_error("Error: process_shoulder_exit called with invalid parameters.");
        }
        return;
    }

    if (event->p1_idx < 0 || event->p1_idx >= ctx->num_particles ||
        event->p2_idx < 0 || event->p2_idx >= ctx->num_particles) {
        if (ctx && ctx->hooks.log_error) {
            ctx->hooks.log_error("Error: process_shoulder_exit called with invalid particle indices.");
        }
        return;
    }

    Particle *p1 = &ctx->particles[event->p1_idx];
    Particle *p2 = &ctx->particles[event->p2_idx];

    // Advance particles to the exact time of the event (should be done by process_event)
    // double dt = event->time - ctx->current_time;
    // if (dt > 1e-12) {
    //     advance_particle_position(p1, dt);
    //     advance_particle_position(p2, dt);
    // }
    // ctx->current_time = event->time;

    // Get current separation with periodic boundary conditions
    double dx = pbc_min_image_delta(p1->x, p2->x, ctx->xsize);
    double dy = pbc_min_image_delta(p1->y, p2->y, ctx->ysize);
    double dist = sqrt(dx*dx + dy*dy);
    double shoulder_radius = ctx->params.sigma * ctx->params.lambda_shoulder;

    // Move particles to be exactly at the shoulder boundary
    if (fabs(dist - shoulder_radius) > 1e-9) {
        // if (ctx->hooks.log_warning) {
        //    ctx->hooks.log_warning("Particles not exactly at shoulder distance (%.12f vs %.12f) at exit. Adjusting.", dist, shoulder_radius);
        // }
        if (dist > 1e-12) { // Avoid division by zero
            double unit_dx = dx / dist;
            double unit_dy = dy / dist;
            double m_tot = p1->m + p2->m;

            double p1_target_x = p1->x - (dist - shoulder_radius) * unit_dx * (p2->m / m_tot);
            double p1_target_y = p1->y - (dist - shoulder_radius) * unit_dy * (p2->m / m_tot);
            double p2_target_x = p2->x + (dist - shoulder_radius) * unit_dx * (p1->m / m_tot);
            double p2_target_y = p2->y + (dist - shoulder_radius) * unit_dy * (p1->m / m_tot);
            
            p1->x = p1_target_x;
            p1->y = p1_target_y;
            p2->x = p2_target_x;
            p2->y = p2_target_y;

            pbc_wrap_position(&p1->x, &p1->y, ctx->xsize, ctx->ysize);
            pbc_wrap_position(&p2->x, &p2->y, ctx->xsize, ctx->ysize);

            dx = pbc_min_image_delta(p1->x, p2->x, ctx->xsize); // Recalculate for nx, ny
            dy = pbc_min_image_delta(p1->y, p2->y, ctx->ysize);
            dist = sqrt(dx*dx + dy*dy); 
            // Debug output for validation
            if (fabs(dist - shoulder_radius) > 1e-8) {
                printf("[DEBUG] Shoulder EXIT: Post-correction distance violation: %.12f (expected %.12f)\n", dist, shoulder_radius);
            }
        }
    }
    
    // Normal vector from p1 to p2
    double nx = dx/dist;
    double ny = dy/dist;

    // Relative velocity
    double dvx = p2->vx - p1->vx;
    double dvy = p2->vy - p1->vy;
    // Normal component of relative velocity (v_rel_n = (v2-v1) . n)
    // Should be positive if they are indeed exiting (moving apart)
    double v_rel_n_initial = dvx * nx + dvy * ny;

    if (v_rel_n_initial <= 1e-9) {
        // If they are not separating at the exit boundary, something is wrong.
        // This could be a stale event or a logic error.
        // if (ctx->hooks.log_warning) {
        //     ctx->hooks.log_warning("Shoulder exit called for non-separating particles (v_rel_n=%.2e). No velocity change applied by U.", v_rel_n_initial);
        // }
        // Mark as out of shoulder and re-predict.
        track_shoulder_state(ctx, event->p1_idx, event->p2_idx, 0);
        invalidate_partner_events(ctx, p1->id);
        invalidate_partner_events(ctx, p2->id);
        predict_all_events_for_particle(ctx, p1);
        predict_all_events_for_particle(ctx, p2);
        return;
    }

    // Increase normal KE by U.
    // KE_normal_final = KE_normal_initial + U
    // 0.5 * mu * v_rel_n_final^2 = 0.5 * mu * v_rel_n_initial^2 + U
    // v_rel_n_final^2 = v_rel_n_initial^2 + 2*U/mu
    double m1 = p1->m;
    double m2 = p2->m;
    double mu = (m1 * m2) / (m1 + m2);
    
    double v_rel_n_initial_sq = v_rel_n_initial * v_rel_n_initial;
    double v_rel_n_final_sq = v_rel_n_initial_sq + (2.0 * ctx->params.U / mu);
    // v_rel_n_final must be positive as they are exiting.
    double v_rel_n_final = sqrt(v_rel_n_final_sq); 

    // Change in relative normal velocity
    double delta_v_rel_n = v_rel_n_final - v_rel_n_initial;
    
    // Impulse J = mu * delta_v_rel_n
    double impulse_j = mu * delta_v_rel_n;

    // Apply impulse
    // p1 velocity change: + (J/m1)*n  (n points from p1 to p2, J is positive, so p1 speeds up along n)
    // p2 velocity change: - (J/m2)*n  (p2 speeds up along -n)
    p1->vx += (impulse_j / m1) * nx;
    p1->vy += (impulse_j / m1) * ny;
    p2->vx -= (impulse_j / m2) * nx;
    p2->vy -= (impulse_j / m2) * ny;

    // Mark particles as no longer in a shoulder interaction
    track_shoulder_state(ctx, event->p1_idx, event->p2_idx, 0);

    // Invalidate other events for these particles as their trajectories changed
    invalidate_partner_events(ctx, p1->id);
    invalidate_partner_events(ctx, p2->id);
    
    // Re-predict all events for both particles
    predict_all_events_for_particle(ctx, p1);
    predict_all_events_for_particle(ctx, p2);
}
