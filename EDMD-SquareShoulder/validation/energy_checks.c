#include "../core/sim_context.h"
#include "energy_audit.h"
#include <math.h>
#include <stdio.h>
#include <stdatomic.h>  // For atomic operations

#define ENERGY_TOLERANCE 1e-6

void validate_energy(SimContext* ctx) {
    double current_energy = 0.0;
    for(int i=0; i<ctx->num_particles; i++) {
        current_energy += 0.5 * ctx->particles[i].m *
                         (ctx->particles[i].vx*ctx->particles[i].vx +
                          ctx->particles[i].vy*ctx->particles[i].vy);
    }
    double delta = fabs(current_energy - ctx->energy.tracked);
    if(delta > ENERGY_TOLERANCE) {
        fprintf(stderr, "Energy violation: Δ=%.6f\n", delta);
        // Optionally: recovery or abort
    }
    ctx->energy.tracked = current_energy;
}
