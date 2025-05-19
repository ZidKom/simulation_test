#pragma once

#include "../core/sim_context.h"

/**
 * @brief Advances a particle's position based on its velocity over a specified time
 * 
 * This function linearly extrapolates the particle's position using its current velocity
 * for the specified time duration. It also applies periodic boundary conditions
 * and updates the particle's cell membership.
 * 
 * @param ctx Pointer to the simulation context (for box size and cell updates)
 * @param p Pointer to the particle to advance
 * @param dt Time duration to advance
 */
void advance_particle_position(SimContext* ctx, Particle* p, double dt);
