#include "../core/sim_context.h"
#include "../core/particle_pool.h"
#include "../core/event_system/hybrid.h"
#include <stdio.h>
#include <stdlib.h>

int main() {
    printf("Starting minimal test...\n");
    
    // Initialize a simple simulation context with just 2 particles
    SimContext ctx;
    double box_size = 10.0;
    int num_particles = 2;
    
    // Initialize context with default parameters
    init_sim_context(&ctx, num_particles, box_size, box_size, 1.0, 0, 1000);
    
    // Setup particles
    ctx.particles[0].id = 0;
    ctx.particles[0].x = 3.0;
    ctx.particles[0].y = 5.0;
    ctx.particles[1].id = 1;
    ctx.particles[1].x = 7.0;
    ctx.particles[1].y = 5.0;
    
    // Test event system
    printf("Testing event allocation...\n");
    Event* event1 = allocate_event_from_pool(&ctx);
    if (!event1) {
        printf("Failed to allocate event1\n");
        return 1;
    }
    
    event1->time = 1.0;
    event1->type = 1; // Arbitrary type for test
    event1->p1_idx = 0;
    event1->p2_idx = 1;
    
    printf("Scheduling event...\n");
    schedule_event(&ctx, event1);
    
    printf("Getting event back...\n");
    Event* retrieved = get_next_master_event(&ctx);
    if (!retrieved) {
        printf("Failed to retrieve event\n");
        return 1;
    }
    
    printf("Retrieved event: time=%.2f, type=%d, p1=%d, p2=%d\n", 
           retrieved->time, retrieved->type, retrieved->p1_idx, retrieved->p2_idx);
    
    printf("Test completed successfully!\n");
    
    // Clean up
    free_sim_context(&ctx);
    
    return 0;
}
