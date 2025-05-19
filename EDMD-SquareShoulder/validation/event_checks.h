#pragma once

#include "../core/sim_context.h"

// Validation tolerance constants
#define COLLISION_TOL 1e-9
#define SHOULDER_TOL 1e-8

/**
 * @brief Validates if a collision event occurs at the expected distance
 * 
 * @param ctx Simulation context
 * @param ev Collision event to validate
 * @return int 1 if valid, 0 if invalid
 */
int validate_collision(SimContext* ctx, Event* ev);

/**
 * @brief Validates if a shoulder event occurs at the expected distance
 * 
 * @param ctx Simulation context
 * @param ev Shoulder event to validate
 * @return int 1 if valid, 0 if invalid
 */
int validate_shoulder_event(SimContext* ctx, Event* ev);
