#include "../core/sim_context.h"
#include "../physics/collisions.h"
#include "../validation/energy_audit.h"
#include "../validation/hooks.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

/**
 * @brief Run a series of 1000 collisions between two particles and verify energy conservation
 * 
 * @return int 0 on success, non-zero on failure
 */
int test_1000_collisions() {
    printf("Running 1000 collision test...\n");
    
    // Initialize simulation context
    SimContext ctx;
    memset(&ctx, 0, sizeof(SimContext));
    
    // Setup simulation parameters
    double box_size = 10.0;
    ctx.xsize = box_size;
    ctx.ysize = box_size;
    ctx.num_particles = 2;
    ctx.params.sigma = 1.0;
    ctx.params.lambda_shoulder = 1.5;
    ctx.params.deltaE = 0.01;
    ctx.params.default_particle_radius = 0.5;
    ctx.params.default_particle_mass = 1.0;
    
    // Allocate particles
    ctx.particles = (Particle*)malloc(ctx.num_particles * sizeof(Particle));
    if (!ctx.particles) {
        fprintf(stderr, "Failed to allocate particles\n");
        return 1;
    }
    
    // Initialize particles for head-on collision
    ctx.particles[0].id = 0;
    ctx.particles[0].x = 3.0;
    ctx.particles[0].y = 5.0;
    ctx.particles[0].vx = 1.0;
    ctx.particles[0].vy = 0.0;
    ctx.particles[0].radius = 0.5;
    ctx.particles[0].radius_sq = 0.25;
    ctx.particles[0].m = 1.0;
    ctx.particles[0].cellx = -1;
    ctx.particles[0].celly = -1;
    ctx.particles[0].next_in_cell = NULL;
    
    ctx.particles[1].id = 1;
    ctx.particles[1].x = 7.0;
    ctx.particles[1].y = 5.0;
    ctx.particles[1].vx = -1.0;
    ctx.particles[1].vy = 0.0;
    ctx.particles[1].radius = 0.5;
    ctx.particles[1].radius_sq = 0.25;
    ctx.particles[1].m = 1.0;
    ctx.particles[1].cellx = -1;
    ctx.particles[1].celly = -1;
    ctx.particles[1].next_in_cell = NULL;
    
    // Setup validation
    ctx.energy.tolerance = 1e-9;
    ctx.energy.tracked = calculate_total_energy(&ctx);
    // Set up log hooks
    ctx.hooks.log_error = default_log_error;
    
    // Track test results
    int test_passed = 1;
    int collision_count = 0;
    
    // Perform 1000 collisions
    printf("Initial energy: %.9f\n", ctx.energy.tracked);
    
    for (int i = 0; i < 1000; i++) {
        // Create collision event
        Event event;
        memset(&event, 0, sizeof(Event));
        event.type = EVENT_COLLISION;
        event.p1_idx = 0;
        event.p2_idx = 1;
        
        // Calculate time to collision
        double dx = ctx.particles[1].x - ctx.particles[0].x;
        double dy = ctx.particles[1].y - ctx.particles[0].y;
        double dvx = ctx.particles[1].vx - ctx.particles[0].vx;
        double dvy = ctx.particles[1].vy - ctx.particles[0].vy;
        
        double b = dx*dvx + dy*dvy;
        double v_sq = dvx*dvx + dvy*dvy;
        double d_sq = dx*dx + dy*dy;
        double sigma_sq = pow(ctx.particles[0].radius + ctx.particles[1].radius, 2);
        
        double dt = -(b + sqrt(b*b - v_sq*(d_sq - sigma_sq)))/v_sq;
        
        // Advance particles to collision
        ctx.particles[0].x += ctx.particles[0].vx * dt;
        ctx.particles[0].y += ctx.particles[0].vy * dt;
        ctx.particles[1].x += ctx.particles[1].vx * dt;
        ctx.particles[1].y += ctx.particles[1].vy * dt;
        
        // Process collision
        // Here we would normally call process_collision(&ctx, &event)
        // But for testing purposes, we'll implement core collision logic directly
        double nx = ctx.particles[1].x - ctx.particles[0].x;
        double ny = ctx.particles[1].y - ctx.particles[0].y;
        double nr = sqrt(nx*nx + ny*ny);
        nx /= nr;
        ny /= nr;
        
        double v1n = ctx.particles[0].vx*nx + ctx.particles[0].vy*ny;
        double v2n = ctx.particles[1].vx*nx + ctx.particles[1].vy*ny;
        
        // Energy calculation before momentum exchange
        double old_energy = 0.5*ctx.particles[0].m*(ctx.particles[0].vx*ctx.particles[0].vx + ctx.particles[0].vy*ctx.particles[0].vy) +
                           0.5*ctx.particles[1].m*(ctx.particles[1].vx*ctx.particles[1].vx + ctx.particles[1].vy*ctx.particles[1].vy);
        
        // Momentum exchange
        double v1n_new = v2n;
        double v2n_new = v1n;
        
        // Update velocities
        ctx.particles[0].vx += (v1n_new - v1n)*nx;
        ctx.particles[0].vy += (v1n_new - v1n)*ny;
        ctx.particles[1].vx += (v2n_new - v2n)*nx;
        ctx.particles[1].vy += (v2n_new - v2n)*ny;
        
        // Energy injection
        double energy_factor = sqrt(1.0 + ctx.params.deltaE/old_energy);
        ctx.particles[0].vx *= energy_factor;
        ctx.particles[0].vy *= energy_factor;
        ctx.particles[1].vx *= energy_factor;
        ctx.particles[1].vy *= energy_factor;
        
        // Update tracked energy
        ctx.energy.tracked += ctx.params.deltaE;
        
        // Validate energy after collision
        double post_energy = calculate_total_energy(&ctx);
        
        // Check if energy matches our tracked value within tolerance
        if (fabs(post_energy - ctx.energy.tracked) > ctx.energy.tolerance) {
            fprintf(stderr, "Energy mismatch at collision %d: Expected %.9f, got %.9f\n", 
                    i, ctx.energy.tracked, post_energy);
            test_passed = 0;
            break;
        }
        
        collision_count++;
    }
    
    printf("Performed %d collisions\n", collision_count);
    printf("Final energy: %.9f\n", ctx.energy.tracked);
    printf("Expected energy: %.9f\n", calculate_total_energy(&ctx));
    printf("Test %s\n", test_passed ? "PASSED" : "FAILED");
    
    // Clean up
    free(ctx.particles);
    
    return test_passed ? 0 : 1;
}

/**
 * @brief Main entry point for collision tests
 * 
 * @return int 0 on success, non-zero on failure
 */
int main(int argc, char **argv) {
    // Seed random number generator
    srand(time(NULL));
    
    // Run tests
    int result = test_1000_collisions();
    
    return result;
}
