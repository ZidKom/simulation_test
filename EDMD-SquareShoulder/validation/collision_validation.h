#ifndef COLLISION_VALIDATION_H
#define COLLISION_VALIDATION_H

#include "../core/sim_context.h"

/**
 * @brief Validates a collision event to ensure numerical precision.
 * 
 * This function checks if the particles are actually at the correct distance
 * for a collision, within a small epsilon tolerance. This helps detect numerical 
 * drift that can occur in simulations.
 * 
 * @param ctx Pointer to the simulation context
 * @param ev Pointer to the collision event to validate
 * @return 1 if the collision is valid, 0 otherwise
 */
int validate_collision(SimContext* ctx, Event* ev);

#endif // COLLISION_VALIDATION_H
