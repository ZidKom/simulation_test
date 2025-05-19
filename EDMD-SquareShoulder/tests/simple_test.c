#include "../core/sim_context.h"
#include "../physics/collisions.h"
#include "../physics/shoulders.h"
#include "../physics/event_prediction.h"
#include "../physics/particle_motion.h"
#include "../physics/particle_interactions.h"
#include "../core/event_system/hybrid.h"
#include "../core/particle_pool.h"
#include "../validation/energy_audit.h"
#include "../utils/pbc.h"
#include "../utils/cell_list.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <assert.h>

/**
 * @brief Simple test to validate the basic functionality
 */
int main() {
    printf("Starting simple test...\n");
    
    // Initialize a simple simulation context with just 2 particles
    SimContext ctx;
    double box_size = 10.0;
    int num_particles = 2;
    
    // Initialize context with default parameters
    init_sim_context(&ctx, num_particles, box_size, box_size, 1.0, 0, 10000); // Use a larger event pool
    
    // Setup validation hooks
    ctx.hooks.log_error = printf;  // Simple error logging
    
    // Set physics parameters
    ctx.params.sigma = 1.0;             // Core diameter
    ctx.params.lambda_shoulder = 1.5;   // Shoulder width multiplier
    ctx.params.U = 0.5;                 // Shoulder potential height
    ctx.params.deltaE = 0.0;            // No energy injection for test
    
    // Setup particles in a simple configuration
    ctx.particles[0].id = 0;
    ctx.particles[0].x = 3.0;
    ctx.particles[0].y = 5.0;
    ctx.particles[0].vx = 1.0;
    ctx.particles[0].vy = 0.0;
    ctx.particles[0].radius = 0.5;
    ctx.particles[0].radius_sq = 0.25;
    ctx.particles[0].m = 1.0;
    ctx.particles[0].in_shoulder = 0;
    
    ctx.particles[1].id = 1;
    ctx.particles[1].x = 7.0;
    ctx.particles[1].y = 5.0;
    ctx.particles[1].vx = -1.0;
    ctx.particles[1].vy = 0.0;
    ctx.particles[1].radius = 0.5;
    ctx.particles[1].radius_sq = 0.25;
    ctx.particles[1].m = 1.0;
    ctx.particles[1].in_shoulder = 0;
    
    // Setup cell lists
    init_cell_lists(&ctx, 10, 10);
    
    // Update particle cell membership
    update_particle_cell(&ctx, &ctx.particles[0]);
    update_particle_cell(&ctx, &ctx.particles[1]);
    
    // Initialize energy tracking 
    init_energy_tracking(&ctx);
    
    printf("Initialization complete\n");
    
    // Test particle_interactions.c functionality
    printf("Testing track_shoulder_state...\n");
    track_shoulder_state(&ctx, 0, 1, 1);
    
    printf("Testing invalidate_partner_events...\n");
    invalidate_partner_events(&ctx, 0);
    
    printf("Test completed successfully!\n");
    
    // Clean up
    free_sim_context(&ctx);
    
    return 0;
}
