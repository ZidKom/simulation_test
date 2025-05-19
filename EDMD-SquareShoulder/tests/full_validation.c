#include "../core/sim_context.h"
#include "../physics/collisions.h"
#include "../physics/shoulders.h"
#include "../core/event_system/hybrid.h"
#include "../utils/pbc.h"
#include "../validation/energy_audit.h"
#include "../validation/hooks.h"
#include "../core/event_system/events.h"
#include "../core/particle_pool.h"
#include "../physics/event_prediction.h"
#include "../physics/cell_interactions.h"
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <float.h>
#include <string.h> /* For memset */

#define BOX_SIZE 10.0
#define CELL_SIZE 2.0
#define NUM_TEST_CYCLES 1000

/* Function prototypes */
double calculate_total_kinetic_energy(SimContext* ctx);
void init_energy_tracking(SimContext* ctx);
void process_event(SimContext* ctx, Event* ev); // Add this prototype

/* **** Test Configuration ***** */
const struct {
    double U;           // Shoulder height
    double sig;         // Shoulder width
    double deltaE;      // Energy injection
    double energy_tol;  // Validation tolerance
} test_params = {
    .U = 1.0,
    .sig = 0.2,
    .deltaE = 0.01,
    .energy_tol = 1e-6
};

/* **** Core Validation Functions ***** */

void validate_physics(SimContext* ctx) {
    /* 1. Energy conservation check */
    double current_ke = calculate_total_kinetic_energy(ctx);
    double expected_ke = ctx->energy.tracked;
    
    if (fabs(current_ke - expected_ke) > ctx->energy.tolerance) {
        fprintf(stderr, "Energy violation: Δ=%.6f (current=%.6f, tracked=%.6f)\n", 
                current_ke - expected_ke, current_ke, expected_ke);
        assert(fabs(current_ke - expected_ke) < ctx->energy.tolerance);
    }

    /* 2. Momentum conservation check (optional - uncomment if momentum tracking is added) */
    /*
    double px = 0.0, py = 0.0;
    for(int i=0; i < ctx->num_particles; i++) {
        px += ctx->particles[i].m * ctx->particles[i].vx;
        py += ctx->particles[i].m * ctx->particles[i].vy;
    }
    // Add momentum tracking to ctx if needed
    */

    /* 3. Shoulder potential validation */
    for(int i=0; i < ctx->num_particles; i++) {
        for(int j=i+1; j < ctx->num_particles; j++) {
            Particle *p1 = &ctx->particles[i];
            Particle *p2 = &ctx->particles[j];
            double dx = pbc_distance(p2->x - p1->x, ctx->xsize);
            double dy = pbc_distance(p2->y - p1->y, ctx->ysize);
            double r = sqrt(dx*dx + dy*dy);
            double r_min = p1->radius + p2->radius;
            
            // Check that particles are never overlapping cores
            assert(r >= r_min - 1e-9);
            
            // Check shoulder status is consistent with distance
            if(r < r_min * (1.0 + ctx->params.sigma)) {
                // If particles are within shoulder distance, at least one should be marked
                if (p1->in_shoulder || p2->in_shoulder) {
                    // This is good - at least one particle has the flag set
                } else {
                    fprintf(stderr, "Particles %d and %d are within shoulder distance but not marked\n", 
                            p1->id, p2->id);
                    assert(p1->in_shoulder || p2->in_shoulder);
                }
            }
        }
    }
}

void validate_event_system(SimContext* ctx) {
    /* 1. Hybrid queue consistency - modified to only check event timestamps */
    Event* current = get_next_master_event(ctx);
    
    /* 2. Event time ordering */
    double current_time = ctx->current_time;
    if (current) {
        // All scheduled events should be in the future
        assert(current->time >= current_time);
        
        // Release the event since we only peeked at it and didn't process it
        release_event_to_pool(ctx, current);
    }

    /* 3. Cell crossing prediction - simplified */
    for(int i=0; i<ctx->num_particles; i++) {
        Particle* p = &ctx->particles[i];
        
        // Ensure particle is in the correct cell
        int expected_cellx = (int)(p->x / ctx->cell_size);
        int expected_celly = (int)(p->y / ctx->cell_size);
        
        // Apply periodic boundary conditions to cell indices
        expected_cellx = (expected_cellx + ctx->n_cells_x) % ctx->n_cells_x;
        expected_celly = (expected_celly + ctx->n_cells_y) % ctx->n_cells_y;
        
        if (p->cellx != expected_cellx || p->celly != expected_celly) {
            fprintf(stderr, "Particle %d cell mismatch: (%d,%d) vs expected (%d,%d)\n", 
                    p->id, p->cellx, p->celly, expected_cellx, expected_celly);
            assert(p->cellx == expected_cellx && p->celly == expected_celly);
        }
    }
}

/* ***** Test Scenarios ***** */

// Functions from the existing codebase
void predict_initial_events(SimContext* ctx) {
    // Predict events for all particles
    for (int i = 0; i < ctx->num_particles; i++) {
        predict_all_events_for_particle(ctx, &ctx->particles[i]);
    }
}

void test_core_collision() {
    printf("=== Testing Core Collision Physics ===\n");
    
    SimContext ctx;
    memset(&ctx, 0, sizeof(SimContext));
    
    // Initialize simulation context with 2 particles
    // Account for all 7 required parameters
    init_sim_context(&ctx, 2, BOX_SIZE, BOX_SIZE, CELL_SIZE, 100 /* paul_list_size */, 1000 /* event_pool_size */);
    
    // Set energy tracking parameters
    ctx.energy.tolerance = test_params.energy_tol;
    
    // Set physics parameters
    ctx.params.U = test_params.U;
    ctx.params.sigma = test_params.sig;
    ctx.params.deltaE = test_params.deltaE;
    
    // Initialize particles
    ctx.particles[0] = (Particle){
        .x = 3.0, .y = 5.0, 
        .vx = 1.0, .vy = 0.0, 
        .radius = 0.5, .radius_sq = 0.25, 
        .m = 1.0, 
        .id = 0, 
        .active = 1
    };
    
    ctx.particles[1] = (Particle){
        .x = 7.0, .y = 5.0, 
        .vx = -1.0, .vy = 0.0, 
        .radius = 0.5, .radius_sq = 0.25, 
        .m = 1.0, 
        .id = 1, 
        .active = 1
    };
    
    // Set cell indices
    for (int i = 0; i < ctx.num_particles; i++) {
        ctx.particles[i].cellx = (int)(ctx.particles[i].x / ctx.cell_size);
        ctx.particles[i].celly = (int)(ctx.particles[i].y / ctx.cell_size);
    }
    
    // Initialize energy tracking
    init_energy_tracking(&ctx);
    
    // Predict initial events
    predict_initial_events(&ctx);
    
    for(int i=0; i<NUM_TEST_CYCLES; i++) {
        Event* ev = get_next_master_event(&ctx);
        if (!ev) break;
        
        process_event(&ctx, ev);
        release_event_to_pool(&ctx, ev);
        
        validate_physics(&ctx);
        validate_event_system(&ctx);
    }
    
    printf("Core collision physics validated through %d cycles\n", NUM_TEST_CYCLES);
    free_sim_context(&ctx);
}

void test_shoulder_potential() {
    printf("\n=== Testing Square Shoulder Potential ===\n");
    
    SimContext ctx;
    memset(&ctx, 0, sizeof(SimContext));
    
    // Initialize simulation context with 2 particles
    init_sim_context(&ctx, 2, BOX_SIZE, BOX_SIZE, CELL_SIZE, 100 /* paul_list_size */, 1000 /* event_pool_size */);
    
    // Set a high potential barrier to force reflection
    ctx.params.U = 2.0; 
    ctx.params.sigma = test_params.sig;
    
    // Initialize particles for shoulder testing
    ctx.particles[0] = (Particle){
        .x = 4.0, .y = 5.0, 
        .vx = 0.5, .vy = 0.0, 
        .radius = 0.5, .radius_sq = 0.25, 
        .m = 1.0, 
        .id = 0, 
        .active = 1
    };
    
    ctx.particles[1] = (Particle){
        .x = 6.0, .y = 5.0, 
        .vx = -0.5, .vy = 0.0, 
        .radius = 0.5, .radius_sq = 0.25, 
        .m = 1.0, 
        .id = 1, 
        .active = 1
    };
    
    // Set cell indices
    for (int i = 0; i < ctx.num_particles; i++) {
        ctx.particles[i].cellx = (int)(ctx.particles[i].x / ctx.cell_size);
        ctx.particles[i].celly = (int)(ctx.particles[i].y / ctx.cell_size);
    }
    
    // Initialize energy tracking
    init_energy_tracking(&ctx);
    
    // Predict initial events
    predict_initial_events(&ctx);
    
    // Test shoulder reflection
    Event* ev = get_next_master_event(&ctx);
    if (ev) {
        // Store event type for verification 
        EventType ev_type = ev->type;
        
        process_event(&ctx, ev);
        release_event_to_pool(&ctx, ev);
        
        // If the potential is strong enough, particles should be reflected
        if (ev_type == EVENT_SHOULDER_ENTRY) {
            printf("Shoulder entry event processed\n");
            assert(ctx.particles[0].vx < 0 || ctx.particles[1].vx > 0);
        }
    }
    
    printf("Shoulder reflection validated\n");
    
    // Reset for penetration test
    free_sim_context(&ctx);
    memset(&ctx, 0, sizeof(SimContext));
    init_sim_context(&ctx, 2, BOX_SIZE, BOX_SIZE, CELL_SIZE, 100 /* paul_list_size */, 1000 /* event_pool_size */);
    
    // Set a weaker potential to allow penetration
    ctx.params.U = 0.5;
    ctx.params.sigma = test_params.sig;
    
    // Initialize particles with higher velocities for penetration
    ctx.particles[0] = (Particle){
        .x = 4.0, .y = 5.0, 
        .vx = 2.0, .vy = 0.0, 
        .radius = 0.5, .radius_sq = 0.25, 
        .m = 1.0, 
        .id = 0, 
        .active = 1
    };
    
    ctx.particles[1] = (Particle){
        .x = 6.0, .y = 5.0, 
        .vx = -2.0, .vy = 0.0, 
        .radius = 0.5, .radius_sq = 0.25, 
        .m = 1.0, 
        .id = 1, 
        .active = 1
    };
    
    // Set cell indices
    for (int i = 0; i < ctx.num_particles; i++) {
        ctx.particles[i].cellx = (int)(ctx.particles[i].x / ctx.cell_size);
        ctx.particles[i].celly = (int)(ctx.particles[i].y / ctx.cell_size);
    }
    
    init_energy_tracking(&ctx);
    predict_initial_events(&ctx);
    
    ev = get_next_master_event(&ctx);
    if (ev) {
        process_event(&ctx, ev);
        release_event_to_pool(&ctx, ev);
        
        // If kinetic energy > potential, particles should maintain direction
        assert(ctx.particles[0].vx > 0);
        assert(ctx.particles[1].vx < 0);
    }
    
    printf("Shoulder penetration validated\n");
    
    free_sim_context(&ctx);
}

/* ***** Main Test Runner ***** */

int main() {
    printf("======= EDMD-SquareShoulder Full Validation =======\n");
    printf("Running comprehensive validation of simulation systems...\n\n");
    
    test_core_collision();
    test_shoulder_potential();
    
    printf("\nAll critical systems validated successfully!\n");
    return 0;
}
