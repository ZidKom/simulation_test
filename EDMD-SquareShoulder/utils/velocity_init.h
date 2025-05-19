#ifndef VELOCITY_INIT_H
#define VELOCITY_INIT_H
#include "../core/sim_context.h"
/**
 * Initializes particle velocities from a uniform distribution with zero COM velocity.
 * Assigns random vx and vy components to each particle, drawn from a uniform
 * distribution between -max_abs_velocity_component and +max_abs_velocity_component.
 * After initial assignment, the velocities are adjusted so that the center of mass
 * velocity of the entire system is zero.
 *
 * @param ctx Pointer to the simulation context, containing particles and their properties.
 * @param max_abs_velocity_component The maximum absolute value for each velocity component (vx, vy).
 */
void initialize_particle_velocities(SimContext* ctx, double max_abs_velocity_component);
#endif
