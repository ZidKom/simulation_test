#include "hooks.h"
#include "../core/sim_context.h"
#include <math.h>
#include <stdio.h>

/**
 * @brief Calculates the current total energy of the system
 * 
 * @param ctx Simulation context
 * @return double Total energy (kinetic + potential energy from shoulders)
 */
double calculate_total_energy(SimContext* ctx) {
    if (!ctx) return 0.0;
    
    // Calculate kinetic energy
    double ke = 0.0;
    for (int i=0; i < ctx->num_particles; i++) {
        Particle* p = &ctx->particles[i];
        ke += 0.5 * p->m * (p->vx*p->vx + p->vy*p->vy);
    }
    
    // Calculate shoulder potential energy
    double pe = 0.0;
    // Iterate through all particles that are in a shoulder interaction
    for (int i=0; i < ctx->num_particles; i++) {
        Particle* p = &ctx->particles[i];
        if (p->in_shoulder) {
            // Add half of the shoulder potential energy U for each particle
            // (will be added twice total, once for each particle in the pair)
            pe += 0.5 * ctx->params.U;
        }
    }
    
    return ke + pe;
}

/**
 * @brief Validates system energy against tracked energy
 * 
 * @param ctx Simulation context
 */
void validate_energy(SimContext* ctx) {
    if (!ctx) return;
    
    // Calculate kinetic energy
    double ke = 0.0;
    for(int i=0; i < ctx->num_particles; i++) {
        ke += 0.5 * ctx->particles[i].m * 
             (ctx->particles[i].vx*ctx->particles[i].vx + 
              ctx->particles[i].vy*ctx->particles[i].vy);
    }
    
    // Calculate shoulder potential energy
    double pe = 0.0;
    for (int i=0; i < ctx->num_particles; i++) {
        Particle* p = &ctx->particles[i];
        if (p->in_shoulder) {
            pe += 0.5 * ctx->params.U;
        }
    }
    
    double total_energy = ke + pe;
    
    if(fabs(total_energy - ctx->energy.tracked) > ctx->energy.tolerance) {
        if (ctx->hooks.log_error) {
            ctx->hooks.log_error("Energy violation: Δ=%.6f (current: KE=%.6f + PE=%.6f = %.6f, tracked=%.6f)",
                               total_energy - ctx->energy.tracked, ke, pe, total_energy, ctx->energy.tracked);
        } else {
            fprintf(stderr, "[ENERGY] Violation: Δ=%.6f (current: KE=%.6f + PE=%.6f = %.6f, tracked=%.6f)\n", 
                   total_energy - ctx->energy.tracked, ke, pe, total_energy, ctx->energy.tracked);
        }
    }
}

/**
 * @brief Updates the tracked energy after a system state change
 * 
 * @param ctx Simulation context
 */
void update_tracked_energy(SimContext* ctx) {
    if (!ctx) return;
    
    // Calculate current total energy
    ctx->energy.tracked = calculate_total_energy(ctx);
}

/**
 * @brief Updates tracked energy after energy injection during a collision
 * 
 * @param ctx Simulation context
 * @param energy_injected Amount of energy added to the system
 */
void track_energy_injection(SimContext* ctx, double energy_injected) {
    if (!ctx) return;
    
    ctx->energy.tracked += energy_injected;
}

/**
 * @brief Initialize energy tracking system
 * 
 * @param ctx Simulation context
 */
void init_energy_tracking(SimContext* ctx) {
    if (!ctx) return;
    
    // Set initial tracked energy
    ctx->energy.tracked = calculate_total_energy(ctx);
    
    // Set default tolerance if not already set
    if (ctx->energy.tolerance <= 0.0) {
        ctx->energy.tolerance = 1e-6; // Default tolerance
    }
}

/**
 * @brief Calculates the total kinetic energy of the system
 * 
 * @param ctx Simulation context
 * @return double Total kinetic energy
 */
double calculate_total_kinetic_energy(SimContext* ctx) {
    double ke = 0.0;
    for (int i = 0; i < ctx->num_particles; i++) {
        Particle* p = &ctx->particles[i];
        ke += 0.5 * p->m * (p->vx*p->vx + p->vy*p->vy);
    }
    return ke;
}
