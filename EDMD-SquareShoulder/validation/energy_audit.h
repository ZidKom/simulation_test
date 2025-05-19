#pragma once

#include "../core/sim_context.h"

// Energy calculation and validation functions
double calculate_total_energy(SimContext* ctx);
double calculate_total_kinetic_energy(SimContext* ctx);
void validate_energy(SimContext* ctx);
void update_tracked_energy(SimContext* ctx);
void track_energy_injection(SimContext* ctx, double energy_injected);
void init_energy_tracking(SimContext* ctx);
