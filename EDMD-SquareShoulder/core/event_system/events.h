#pragma once
#include "../sim_context.h" // SimContext now defines Event and Particle

// Event system interface
void init_event_system(SimContext *ctx); // Initializes event pool, Paul lists, BST within SimContext
void schedule_event(SimContext *ctx, Event *ev);
Event *get_next_master_event(SimContext *ctx); // Use consistent name with hybrid.h

// Event pool management specific functions (if not solely handled by init_event_system)
Event* allocate_event_from_pool(SimContext *ctx);
void release_event_to_pool(SimContext *ctx, Event *ev); // Renamed from release_event for clarity

// Event creation functions
// These now take SimContext to access event pool and current_time
Event *create_collision_event(SimContext *ctx, Particle *p1, Particle *p2, double event_abs_time);
Event *create_shoulder_event(SimContext *ctx, Particle *p1, Particle *p2, double event_abs_time, EventType type);
Event *create_cell_crossing_event(SimContext *ctx, Particle *p1, double event_abs_time);
Event *create_system_event(SimContext* ctx, double event_abs_time, EventType type);
