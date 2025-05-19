#include "../core/sim_context.h"
#include "collisions.h"
#include "physics.h"  // For NEVER definition
#include "particle_interactions.h" // For invalidate_partner_events
#include "../validation/energy_audit.h" // For track_energy_injection
#include <math.h>
#include <immintrin.h>
#include "../utils/pbc.h"
#include "../core/event_system/hybrid.h" // For invalidate_events_for_particle
#include "event_prediction.h"

// Vectorized, impulse-based energy-injected collision with validation
void process_collision(SimContext* ctx, Event* ev) {
    if (!ctx || !ev) return;
    
    // Get particle indices from the event
    int p1_idx = ev->p1_idx;
    int p2_idx = ev->p2_idx;
    
    if (p1_idx < 0 || p1_idx >= ctx->num_particles || 
        p2_idx < 0 || p2_idx >= ctx->num_particles) {
        // It's good practice to log this error if a logging mechanism is available
        // For example: if (ctx->hooks.log_error) ctx->hooks.log_error("Invalid particle indices in process_collision");
        return;
    }
    
    Particle *p1 = &ctx->particles[p1_idx];
    Particle *p2 = &ctx->particles[p2_idx];
    
    // Use correct periodic boundary condition functions for distance calculation
    double dx = pbc_min_image_delta(p1->x, p2->x, ctx->xsize);
    double dy = pbc_min_image_delta(p1->y, p2->y, ctx->ysize);
    double dist_sq = dx*dx + dy*dy; // Work with squared distance as long as possible
    
    // It's possible particles have slightly overlapped due to floating point precision
    // or if the event time was not exact.
    // For a robust simulation, we should handle this, but for now, assume they are at contact.
    // A check like: if (dist_sq < (p1->radius + p2->radius)*(p1->radius + p2->radius) - 1e-9) might be needed.

    double dist = sqrt(dist_sq);
    if (dist < 1e-12) { // Avoid division by zero if particles are at the exact same spot
        if (ctx && ctx->hooks.log_error) ctx->hooks.log_error("Collision between overlapping particles p%d and p%d at time %.12f", p1_idx, p2_idx, ctx->current_time);
        // Invalidate and re-predict to be safe, then return
        invalidate_events_for_particle(ctx, p1_idx);
        invalidate_events_for_particle(ctx, p2_idx);
        predict_all_events_for_particle(ctx, p1);
        predict_all_events_for_particle(ctx, p2);
        return;
    }
    
    // Calculate normal vector at point of impact
    double nx = dx/dist;
    double ny = dy/dist;
    
    // Corrected relative velocity calculation for v_rel_n_approach
    // dvx and dvy should represent (v1 - v2) if v_rel_n_approach is (v1-v2).n
    // Or, if v_rel_n_approach is (v_particle1 - v_particle2) . normal_from_1_to_2
    double dvx_p1_minus_p2 = p1->vx - p2->vx;
    double dvy_p1_minus_p2 = p1->vy - p2->vy;
    double v_rel_n_approach = (dvx_p1_minus_p2 * nx + dvy_p1_minus_p2 * ny); // This is (v1-v2).n

    // If particles are not approaching (v_rel_n_approach >= 0), they are moving apart or tangentially.
    // A collision event should only be processed if they are approaching.
    // A small positive value can occur due to floating point inaccuracies.
    if (v_rel_n_approach >= -1e-9) { // If separating or moving tangentially (use -1e-9 to allow for small numerical noise if they are perfectly at contact and should collide)
                                     // If v_rel_n_approach is positive, they are separating.
                                     // If v_rel_n_approach is zero, they are moving tangentially.
                                     // If v_rel_n_approach is negative, they are approaching.
                                     // The impulse calculation uses -(1+e) * v_rel_n_approach. If v_rel_n_approach is positive (separating), impulse would be negative, pushing them together. This is wrong.
                                     // Collision should occur if v_rel_n_approach < 0.
        if (ctx && ctx->hooks.log_warning && v_rel_n_approach > 1e-9) { // Log if clearly separating
             ctx->hooks.log_warning("Collision event for separating particles p%d-p%d at t=%.12f, v_rel_n=%.3e. Skipping impulse.", p1_idx, p2_idx, ev->time, v_rel_n_approach);
        }
        // Invalidate and re-predict to clear this potentially stale event.
        invalidate_events_for_particle(ctx, p1_idx);
        invalidate_events_for_particle(ctx, p2_idx);
        predict_all_events_for_particle(ctx, p1);
        predict_all_events_for_particle(ctx, p2);
        return;
    }
    
    // Calculate masses and reduced mass
    double m1 = p1->m;
    double m2 = p2->m;
    double mu = (m1*m2)/(m1+m2);
    
    // The prompt's fix for v_rel_n_approach was:
    // double dvx = p1->vx - p2->vx;
    // double dvy = p1->vy - p2->vy;
    // double v_rel_n_approach = (dvx * nx + dvy * ny); // approaching
    // This definition means v_rel_n_approach is negative if approaching.
    // The impulse J = (-(1 + restitution_coef) * v_rel_n_approach) / (inv_mass1 + inv_mass2);
    // If v_rel_n_approach is negative, then -v_rel_n_approach is positive. So J is positive.
    // p1->vx += J * inv_mass1 * nx; (nx is from 1 to 2, so p1 moves along n)
    // p2->vx -= J * inv_mass2 * nx; (p2 moves against n)
    // This seems correct. My previous `v_rel_n_approach_magnitude = -(dvx * nx + dvy * ny);` was to make it positive.
    // Let's stick to the prompt's `v_rel_n_approach` which is negative for approach.

    // Re-evaluating the impulse logic from the prompt's context:
    // dvx = p1->vx - p2->vx; dvy = p1->vy - p2->vy;
    // v_rel_n_approach = (dvx * nx + dvy * ny); // Negative if approaching
    
    // Standard elastic collision restitution_coef = 1.0
    double restitution_coef = 1.0; // Assuming elastic for now, can be param
    double inv_mass1 = 1.0 / m1; // Assuming m1, m2 > 0
    double inv_mass2 = 1.0 / m2;


    // Energy injection logic from prompt's context (e_eff)
    double e_eff = restitution_coef;
    if (ctx->params.deltaE > 0 && v_rel_n_approach < -1e-9) { // v_rel_n_approach is negative
        double v_rel_n_approach_sq = v_rel_n_approach * v_rel_n_approach;
        if (mu * v_rel_n_approach_sq < 1e-12) {
            // KE too small
        } else {
            double term_under_sqrt = 1.0 + (2.0 * ctx->params.deltaE) / (mu * v_rel_n_approach_sq);
            if (term_under_sqrt < 0) term_under_sqrt = 0; // Should not happen with deltaE > 0
            e_eff = sqrt(term_under_sqrt); // This e_eff is the "effective coefficient of restitution"
            if (ctx->hooks.validate_collision) track_energy_injection(ctx, ctx->params.deltaE); // Assuming track_energy_injection is part of validation hooks or similar
        }
    }
    
    // Impulse J = mu * (1 + e_eff) * (-v_rel_n_approach)
    // Or, J = (-(1+e_eff) * v_rel_n_approach) / (inv_mass1 + inv_mass2) if using the common formula structure.
    // Let J_scalar = mu * (restitution_coef + e_eff) * (-v_rel_n_approach); where v_rel_n_approach is negative.
    // The formula for J is typically: J = (-(1+e) * m1*m2/(m1+m2) * v_rel_n_initial)
    // where v_rel_n_initial = (v1-v2).n. If v_rel_n_initial is negative (approaching), J is positive.
    double J_scalar = mu * (1.0 + e_eff) * (-v_rel_n_approach); // -v_rel_n_approach is positive approach speed

    p1->vx += J_scalar * nx * inv_mass1;
    p1->vy += J_scalar * ny * inv_mass1;
    p2->vx -= J_scalar * nx * inv_mass2;
    p2->vy -= J_scalar * ny * inv_mass2;

    // Remove the SIMD block for velocity updates.
    
    // Invalidate old events for both particles and their interaction partners
    invalidate_partner_events(ctx, p1_idx);
    invalidate_partner_events(ctx, p2_idx);
    
    // Predict new events for both particles
    predict_all_events_for_particle(ctx, p1);
    predict_all_events_for_particle(ctx, p2);
    
    // Update collision counters
    p1->coll++;
    p2->coll++;
    
    // Validation hook for collision if defined
    if (ctx->hooks.validate_collision) {
        ctx->hooks.validate_collision(ctx, ev);
    }

    // ======== POSITION CORRECTION ==========
    // Calculate exact contact point positions
    dx = pbc_min_image_delta(p1->x, p2->x, ctx->xsize);
    dy = pbc_min_image_delta(p1->y, p2->y, ctx->ysize);
    dist = sqrt(dx*dx + dy*dy);
    double required_dist = p1->radius + p2->radius;
    if (fabs(dist - required_dist) > 1e-12) {
        // Calculate correction vector
        double correction = (dist - required_dist)/dist;
        double mass_ratio_p1 = p2->m/(p1->m + p2->m);
        double mass_ratio_p2 = p1->m/(p1->m + p2->m);
        // Move exactly to contact point
        p1->x -= dx * correction * mass_ratio_p1;
        p1->y -= dy * correction * mass_ratio_p1;
        p2->x += dx * correction * mass_ratio_p2;
        p2->y += dy * correction * mass_ratio_p2;
        // PBC consistency
        pbc_wrap_position(&p1->x, &p1->y, ctx->xsize, ctx->ysize);
        pbc_wrap_position(&p2->x, &p2->y, ctx->xsize, ctx->ysize);
        // Debug output for validation
        dx = pbc_min_image_delta(p1->x, p2->x, ctx->xsize);
        dy = pbc_min_image_delta(p1->y, p2->y, ctx->ysize);
        dist = sqrt(dx*dx + dy*dy);
        if (fabs(dist - required_dist) > 1e-8) {
            printf("[DEBUG] Collision: Post-correction distance violation: %.12f (expected %.12f)\n", dist, required_dist);
        }
    }
}

// Calculate collision time using correct box size from context
double calculate_collision_time(SimContext* ctx, Particle* p1, Particle* p2) {
    if (!ctx || !p1 || !p2) return NEVER;
    
    // Use correct periodic boundary condition functions
    double dx = pbc_min_image_delta(p1->x, p2->x, ctx->xsize);
    double dy = pbc_min_image_delta(p1->y, p2->y, ctx->ysize);
    double dvx = p2->vx - p1->vx;
    double dvy = p2->vy - p1->vy;
    double r = p1->radius + p2->radius;
    
    // Use the numerically stable solver from event_prediction.c
    double t = solve_collision_time(dx, dy, dvx, dvy, r);
    
    return (t != DBL_MAX) ? t : NEVER;
}
