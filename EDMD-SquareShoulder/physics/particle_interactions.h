#ifndef PARTICLE_INTERACTIONS_H
#define PARTICLE_INTERACTIONS_H

#include "../core/sim_context.h"

/**
 * @brief Tracks shoulder interaction states between particle pairs
 * 
 * @param ctx Simulation context
 * @param p1_idx Index of first particle
 * @param p2_idx Index of second particle
 * @param in_shoulder Flag indicating if particles are now in shoulder interaction (1) or not (0)
 */
void track_shoulder_state(SimContext* ctx, int p1_idx, int p2_idx, int in_shoulder);

/**
 * @brief Invalidates all events for a particle and its interaction partners
 * 
 * @param ctx Simulation context
 * @param p_idx Index of the particle whose events should be invalidated
 */
void invalidate_partner_events(SimContext* ctx, int p_idx);

#endif // PARTICLE_INTERACTIONS_H
