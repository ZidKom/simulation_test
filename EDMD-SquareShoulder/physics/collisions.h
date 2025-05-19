#pragma once

#include "../core/sim_context.h"

// Forward declarations
double calculate_collision_time(SimContext *ctx, Particle *p1, Particle *p2);
void process_collision(SimContext *ctx, Event *event);
