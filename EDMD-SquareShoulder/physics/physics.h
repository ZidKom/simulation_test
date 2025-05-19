#pragma once

#include <float.h> // For DBL_MAX
#include "../core/sim_context.h" // For SimContext, Particle, Event types

#define NEVER DBL_MAX

// Function prototypes from physics module files

// From collisions.c
// Calculates the time until the next collision between two particles.
double calculate_collision_time(SimContext *ctx, Particle *p1, Particle *p2);
// Processes a collision event, updating particle velocities and predicting new events.
void process_collision(SimContext *ctx, Event *event);

// From shoulders.c
// Handles a shoulder entry event, adjusts velocities, and predicts subsequent events.
void handle_shoulder_entry(SimContext *ctx, Event *event);
// Processes a shoulder exit event
void process_shoulder_exit(SimContext *ctx, Event *event);

// From event_prediction.c
// Predicts all potential future events (collisions, cell crossings, shoulder entries) for a single particle.
void predict_all_events_for_particle(SimContext* ctx, Particle* p);
// Predicts cell crossing events for all particles.
void predict_cell_crossing_events(SimContext* ctx);
// Numerically stable solver for quadratic equation to find event times
double solve_collision_time(double dx, double dy, double dvx, double dvy, double sigma_sum);
// Calculate time until particles reach the shoulder entry boundary
double calculate_shoulder_entry_time(SimContext *ctx, Particle *p1, Particle *p2);
// Calculate time until particles reach the shoulder exit boundary
double calculate_shoulder_exit_time(SimContext *ctx, Particle *p1, Particle *p2);
