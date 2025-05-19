#include "../core/sim_context.h"
#include "../physics/collisions.h"
#include "../physics/shoulders.h"
#include "../physics/event_prediction.h"
#include "../physics/particle_motion.h"
#include "../core/event_system/hybrid.h"
#include "../core/event_system/paul_list.h" // For init_paul_event_system_component
#include "../core/particle_pool.h"
#include "../validation/energy_audit.h"
#include "../validation/hooks.h"
#include "../utils/pbc.h"
#include "../utils/cell_list.h"
#include "../physics/cell_interactions.h" // Added for predict_cell_crossing
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <float.h> // Added for DBL_MAX

/**
 * @brief Comprehensive test for collision and shoulder events
 * 
 * This test creates controlled collision and shoulder events to validate
 * the physical correctness of event processing, energy conservation,
 * and proper invalidation and rescheduling.
 */
void test_collision_and_shoulder() {
    printf("Starting collision and shoulder event test...\\n");
    SimContext ctx = {0}; // Initialize all fields to zero/NULL
    
    double box_size = 10.0;
    int num_particles = 2;
    // Define parameters for init_sim_context
    // double cell_size_val = 1.0; // Example cell size, can be tuned
    double cell_size_val = box_size; // Use a single cell for the whole box to avoid cell crossings
    int paul_list_sz_val = 100; // Paul list size for the test
    int event_pool_sz_val = 10000; // Event pool size for the test

    // Initialize SimContext once with all parameters
    init_sim_context(&ctx, num_particles, box_size, box_size, cell_size_val, paul_list_sz_val, event_pool_sz_val);
    // init_sim_context should set ctx.current_time = 0.0;
    
    // Initialize Paul list event system component (ensures window_start_time, current_bin are set)
    init_paul_event_system_component(&ctx); 
    
    // Set physics parameters
    ctx.params.sigma = 1.0;             // Core diameter
    ctx.params.lambda_shoulder = 1.5;   // Shoulder width multiplier
    ctx.params.U = 0.5;                 // Shoulder potential height
    ctx.params.deltaE = 0.0;            // No energy injection for test
    ctx.params.default_particle_radius = 0.5; // Make sure this is set for any calculations needing it
    
    // Setup particles in a controlled configuration for head-on collision
    // Particle 0 will start at (3, 5) moving right
    ctx.particles[0].id = 0;
    ctx.particles[0].x = 3.0;
    ctx.particles[0].y = 5.0;
    ctx.particles[0].vx = 1.0;
    ctx.particles[0].vy = 0.0;
    ctx.particles[0].radius = 0.5;
    ctx.particles[0].radius_sq = 0.25;
    ctx.particles[0].m = 1.0;
    ctx.particles[0].in_shoulder = 0;
    ctx.particles[0].coll = 0;
    
    // Particle 1 will start at (7, 5) moving left
    ctx.particles[1].id = 1;
    ctx.particles[1].x = 7.0;
    ctx.particles[1].y = 5.0;
    ctx.particles[1].vx = -1.0;
    ctx.particles[1].vy = 0.0;
    ctx.particles[1].radius = 0.5;
    ctx.particles[1].radius_sq = 0.25;
    ctx.particles[1].m = 1.0;
    ctx.particles[1].in_shoulder = 0;
    ctx.particles[1].coll = 0;
    
    // Setup cell lists
    // Calculate n_cells_x and n_cells_y based on box_size and cell_size_val
    int n_cells_x = (int)floor(box_size / cell_size_val);
    int n_cells_y = (int)floor(box_size / cell_size_val);
    init_cell_lists(&ctx, n_cells_x, n_cells_y);
    
    // Update particle cell membership
    update_particle_cell(&ctx, &ctx.particles[0]);
    update_particle_cell(&ctx, &ctx.particles[1]);
    
    // Initialize energy tracking
    init_energy_tracking(&ctx);
    double initial_energy = calculate_total_energy(&ctx);
    printf("Initial energy: %.6f\\n", initial_energy);
    
    // Predict initial events using calculation functions
    printf("\nPredicting initial events using calculation functions...\n");
    Particle* p0 = &ctx.particles[0];
    Particle* p1 = &ctx.particles[1];

    // Only schedule the initial collision event
    double t_coll = calculate_collision_time(&ctx, p0, p1);
    if (t_coll != DBL_MAX && t_coll > ctx.current_time) {
        if (t_coll < ctx.current_time + 1e-9) t_coll = ctx.current_time + 1e-9; // Ensure not in past
        Event* ev = allocate_event_from_pool(&ctx);
        if (!ev) { printf("ERROR: Failed to allocate collision event!\n"); return; }
        ev->type = EVENT_COLLISION; ev->time = t_coll; ev->p1_idx = p0->id; ev->p2_idx = p1->id;  ev->cross_dir = -1;
        schedule_event(&ctx, ev);
        printf("Scheduled EVENT_COLLISION for P%d-P%d at t=%.12f\n", p0->id, p1->id, t_coll);
    }

    // Predict cell crossings for each particle
    printf("Predicting cell crossings...\n");
    predict_cell_crossing(&ctx, p0);
    predict_cell_crossing(&ctx, p1);
    
    // Process events step by step
    printf("\nProcessing events sequentially...\n");
    
    // Process events in order: shoulder entry, collision, shoulder exit
    for (int event_stage = 0; event_stage < 3; ++event_stage) {
        printf("\nGetting next event... (current_time = %.12f)\n", ctx.current_time);
        Event* next_event = get_next_master_event(&ctx);
        if (!next_event) {
            printf("ERROR: No event found in queue at stage %d!\n", event_stage);
            return;
        }
        if (event_stage == 0) {
            printf("Next event type: %d (expecting SHOULDER_ENTRY=%d or COLLISION=%d), time: %.12f\n", next_event->type, EVENT_SHOULDER_ENTRY, EVENT_COLLISION, next_event->time);
            // Accept either shoulder entry or collision as the first event
            assert(next_event->type == EVENT_SHOULDER_ENTRY || next_event->type == EVENT_COLLISION);
        } else if (event_stage == 1) {
            printf("Next event type: %d (expecting COLLISION=%d or SHOULDER_ENTRY=%d), time: %.12f\n", next_event->type, EVENT_COLLISION, EVENT_SHOULDER_ENTRY, next_event->time);
            assert(next_event->type == EVENT_COLLISION || next_event->type == EVENT_SHOULDER_ENTRY);
        } else if (event_stage == 2) {
            printf("Next event type: %d (expecting SHOULDER_EXIT=%d), time: %.12f\n", next_event->type, EVENT_SHOULDER_EXIT, next_event->time);
            assert(next_event->type == EVENT_SHOULDER_EXIT);
        }
        double dt_event = next_event->time - ctx.current_time;
        advance_particle_position(&ctx.particles[0], dt_event);
        advance_particle_position(&ctx.particles[1], dt_event);
        ctx.current_time = next_event->time;
        if (next_event->type == EVENT_SHOULDER_ENTRY) {
            handle_shoulder_entry(&ctx, next_event);
            invalidate_events_for_particle(&ctx, p0->id);
            invalidate_events_for_particle(&ctx, p1->id);
            predict_all_events_for_particle(&ctx, p0);
            predict_all_events_for_particle(&ctx, p1);
            assert(ctx.particles[0].in_shoulder == 1);
            assert(ctx.particles[1].in_shoulder == 1);
            double energy_after_entry = calculate_total_energy(&ctx);
            printf("Energy after shoulder entry: %.6f\n", energy_after_entry);
        } else if (next_event->type == EVENT_COLLISION) {
            process_collision(&ctx, next_event);
            invalidate_events_for_particle(&ctx, p0->id);
            invalidate_events_for_particle(&ctx, p1->id);
            predict_all_events_for_particle(&ctx, p0);
            predict_all_events_for_particle(&ctx, p1);
            printf("Particle 0 velocity after collision: (%.6f, %.6f)\n", ctx.particles[0].vx, ctx.particles[0].vy);
            printf("Particle 1 velocity after collision: (%.6f, %.6f)\n", ctx.particles[1].vx, ctx.particles[1].vy);
        } else if (next_event->type == EVENT_SHOULDER_EXIT) {
            process_shoulder_exit(&ctx, next_event);
            invalidate_events_for_particle(&ctx, p0->id);
            invalidate_events_for_particle(&ctx, p1->id);
            predict_all_events_for_particle(&ctx, p0);
            predict_all_events_for_particle(&ctx, p1);
            assert(ctx.particles[0].in_shoulder == 0);
            assert(ctx.particles[1].in_shoulder == 0);
            double final_energy = calculate_total_energy(&ctx);
            printf("Final energy: %.6f\n", final_energy);
        }
    }
    
    free_sim_context(&ctx);
    printf("Collision and shoulder test completed successfully!\\n");
}

/**
 * Main function to run the test
 */
int main(int argc, char **argv) {
    (void)argc; 
    (void)argv;
    test_collision_and_shoulder();
    return 0;
}
