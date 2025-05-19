#include <stdlib.h>
#include <math.h>
#include "velocity_init.h"

void initialize_particle_velocities(SimContext* ctx, double max_abs_velocity_component) {
    if (!ctx || ctx->num_particles == 0) {
        if (ctx && ctx->hooks.log_warning) {
            ctx->hooks.log_warning("initialize_particle_velocities: No particles to initialize or null context.");
        }
        return;
    }
    double total_mass = 0.0;
    double sum_px = 0.0;
    double sum_py = 0.0;
    // Assign random velocities and accumulate total momentum
    for (int i = 0; i < ctx->num_particles; ++i) {
        Particle* p = &ctx->particles[i];
        p->vx = ((double)rand() / RAND_MAX) * 2.0 * max_abs_velocity_component - max_abs_velocity_component;
        p->vy = ((double)rand() / RAND_MAX) * 2.0 * max_abs_velocity_component - max_abs_velocity_component;
        sum_px += p->m * p->vx;
        sum_py += p->m * p->vy;
        total_mass += p->m;
    }
    if (fabs(total_mass) < 1e-12) {
        if (ctx->hooks.log_warning) {
            ctx->hooks.log_warning("initialize_particle_velocities: Total mass of particles is near zero. Cannot ensure zero COM velocity.");
        }
        return;
    }
    double com_vx = sum_px / total_mass;
    double com_vy = sum_py / total_mass;
    // Subtract COM velocity from each particle
    for (int i = 0; i < ctx->num_particles; ++i) {
        Particle* p = &ctx->particles[i];
        p->vx -= com_vx;
        p->vy -= com_vy;
    }
    // Optionally log final COM momentum
    if (ctx->hooks.log_info) {
        double final_sum_px = 0.0, final_sum_py = 0.0;
        for (int i = 0; i < ctx->num_particles; ++i) {
            Particle* p = &ctx->particles[i];
            final_sum_px += p->m * p->vx;
            final_sum_py += p->m * p->vy;
        }
        ctx->hooks.log_info("Particle velocities initialized. Max component: %.2f. COM momentum Px=%.3e, Py=%.3e.",
                           max_abs_velocity_component, final_sum_px, final_sum_py);
    }
}
