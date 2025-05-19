#ifndef SIMULATION_HOOKS_H
#define SIMULATION_HOOKS_H

#include "core/sim_context.h"

// Define the simulation hooks structure
typedef struct SimulationHooks {
    void (*log_error)(const char* fmt, ...);
    void (*log_warning)(const char* fmt, ...);
    void (*log_info)(const char* fmt, ...);
    void (*validate_collision)(SimContext* ctx, Event* ev);
} SimulationHooks;

// Function to create default validation hooks
SimulationHooks create_default_validation_hooks();

#endif // SIMULATION_HOOKS_H
