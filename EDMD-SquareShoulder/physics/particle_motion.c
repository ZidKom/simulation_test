#include "event_prediction.h"
#include "../utils/pbc.h"
#include "../utils/cell_list.h" // For update_particle_cell
#include "../core/sim_context.h" // For SimContext definition (though already in .h)

/**
 * @brief Advances a particle's position based on its velocity over a specified time
 * 
 * This function linearly extrapolates the particle's position using its current velocity
 * for the specified time duration.
 * 
 * @param p Pointer to the particle to advance
 * @param dt Time duration to advance
 */
void advance_particle_position(SimContext* ctx, Particle* p, double dt) {
    if (!ctx || !p || dt < 0) return; // Safety check
    
    // Simple linear motion update
    p->x += p->vx * dt;
    p->y += p->vy * dt;
    
    // Apply periodic boundary conditions
    pbc_wrap_position(&p->x, &p->y, ctx->xsize, ctx->ysize);
    
    // Update the particle's cell based on its new position
    update_particle_cell(ctx, p);
}
