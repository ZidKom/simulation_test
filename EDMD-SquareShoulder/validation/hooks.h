#pragma once

#include "../core/sim_context.h"
#include <stdio.h>
#include "../simulation_hooks.h" // Use canonical SimulationHooks definition

// Callback function types for the simulation hooks
typedef void (*LogFunctionPtr)(const char* format, ...);
typedef void (*EventHandlerPtr)(SimContext* ctx, Event* event);
typedef void (*ParticleUpdateHandlerPtr)(SimContext* ctx, int particle_idx);

// Function declarations
void default_log_error(const char* format, ...);
void validate_collision_distance(SimContext* ctx, Event* ev);
void validate_cell_indices(SimContext* ctx);
void validate_cell_crossing(SimContext* ctx, Event* ev);
void validate_initial_conditions(SimContext* ctx);
SimulationHooks create_default_validation_hooks();

#ifndef VALIDATION_HOOKS_H
#define VALIDATION_HOOKS_H

#endif // VALIDATION_HOOKS_H
