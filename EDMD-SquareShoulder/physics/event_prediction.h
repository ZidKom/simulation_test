#pragma once

#include "../core/sim_context.h"

// Forward declarations
void predict_all_events_for_particle(SimContext* ctx, Particle* p);
void predict_pair_events(SimContext* ctx, Particle* p1, Particle* p2);
void predict_cell_crossing_events(SimContext* ctx);
double calculate_shoulder_entry_time(SimContext *ctx, Particle *p1, Particle *p2);
double calculate_shoulder_exit_time(SimContext *ctx, Particle *p1, Particle *p2);
double solve_collision_time(double dx, double dy, double dvx, double dvy, double sigma_sum);
