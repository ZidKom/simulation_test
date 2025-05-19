#include "../core/sim_context.h"
#include "../utils/pbc.h"
#include <math.h>
#include <float.h>

#define COLLISION_TOL 1e-9
#define SHOULDER_TOL 1e-8

/**
 * @brief Validates if a shoulder event occurs at the expected distance
 * 
 * @param ctx Simulation context
 * @param ev Shoulder event to validate
 * @return int 1 if valid, 0 if invalid
 */
int validate_shoulder_event(SimContext* ctx, Event* ev) {
    // Get particle pointers from indices
    if (!ctx || !ev || ev->p1_idx < 0 || ev->p2_idx < 0 || 
        ev->p1_idx >= ctx->num_particles || ev->p2_idx >= ctx->num_particles) {
        return 0;
    }
    
    Particle* p1 = &ctx->particles[ev->p1_idx];
    Particle* p2 = &ctx->particles[ev->p2_idx];
    
    double shoulder_radius = ctx->params.sigma * ctx->params.lambda_shoulder;
    
    // Calculate distance using periodic boundary conditions
    double dx = pbc_min_image_delta(p1->x, p2->x, ctx->xsize);
    double dy = pbc_min_image_delta(p1->y, p2->y, ctx->ysize);
    double dist = sqrt(dx*dx + dy*dy);
    
    int valid = 0;
    
    if(ev->type == EVENT_SHOULDER_ENTRY) {
        // Entering shoulder - distance should be exactly shoulder_radius
        valid = (fabs(dist - shoulder_radius) <= SHOULDER_TOL);
        
        // Additional check: particles should be approaching each other
        double dvx = p1->vx - p2->vx;
        double dvy = p1->vy - p2->vy;
        double dr_dv = dx*dvx + dy*dvy;
        
        if (dr_dv > 0) { // Particles moving apart
            valid = 0; // Should not be entering shoulder if moving apart
            if (ctx->hooks.log_error) {
                ctx->hooks.log_error("Shoulder entry validation: particles moving apart (dr·dv=%.12f)", dr_dv);
            }
        }
    } else if(ev->type == EVENT_SHOULDER_EXIT) {
        // Exiting shoulder - distance should be exactly shoulder_radius
        valid = (fabs(dist - shoulder_radius) <= SHOULDER_TOL);
        
        // Additional check: particles should be moving away from each other
        double dvx = p1->vx - p2->vx;
        double dvy = p1->vy - p2->vy;
        double dr_dv = dx*dvx + dy*dvy;
        
        if (dr_dv < 0) { // Particles approaching
            valid = 0; // Should not be exiting shoulder if approaching
            if (ctx->hooks.log_error) {
                ctx->hooks.log_error("Shoulder exit validation: particles approaching (dr·dv=%.12f)", dr_dv);
            }
        }
    }
    
    if(!valid && ctx->hooks.log_error) {
        ctx->hooks.log_error("Shoulder event %s invalid: dist=%.12f (expected:%.12f), dx=%.12f, dy=%.12f",
                       (ev->type == EVENT_SHOULDER_ENTRY) ? "entry" : "exit",
                       dist, shoulder_radius, dx, dy);
        // Additional diagnostics
        double dvx = p1->vx - p2->vx;
        double dvy = p1->vy - p2->vy;
        double dr_dv = dx*dvx + dy*dvy;
        ctx->hooks.log_error("Particles: p1(id:%d, %.12f,%.12f,%.12f,%.12f) p2(id:%d, %.12f,%.12f,%.12f,%.12f), dr·dv=%.12f",
                       p1->id, p1->x, p1->y, p1->vx, p1->vy, 
                       p2->id, p2->x, p2->y, p2->vx, p2->vy, dr_dv);
    }
    
    return valid;
}
