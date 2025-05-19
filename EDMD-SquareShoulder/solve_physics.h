#pragma once
#include "core/sim_context.h"
#include "physics/physics.h"

double solve_collision_time(double dx, double dy, double dvx, double dvy, double sigma_sum);
double calculate_shoulder_entry_time(SimContext *ctx, Particle *p1, Particle *p2);
double calculate_shoulder_exit_time(SimContext *ctx, Particle *p1, Particle *p2);
void predict_cell_crossing_events(SimContext* ctx);
