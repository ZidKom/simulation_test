#include "../core/sim_context.h"
#include "../physics/collisions.h"
#include "../physics/shoulders.h"
#include "../physics/event_prediction.h" // Corrected include for predict_all_events_for_particle
#include "../physics/particle_motion.h"
#include "../core/event_system/hybrid.h"
#include "../core/particle_pool.h" // For init_particle_pool, free_particle_pool
#include "../validation/energy_audit.h" // For calculate_total_kinetic_energy
#include "../utils/pbc.h" // For pbc_min_image_delta (though not strictly needed for these unit tests if box is large)
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <float.h> // For DBL_MAX

// Helper function to initialize a basic SimContext for testing
void setup_test_context(SimContext* ctx, int num_particles, double box_size, double deltaE_val, double U_val) {
    init_sim_context(ctx, num_particles, box_size, box_size, box_size, 100, 10000); // Use a larger event pool
    ctx->params.sigma = 1.0;
    ctx->params.lambda_shoulder = 1.5; // Shoulder is sigma * lambda_shoulder = 1.5
    ctx->params.deltaE = deltaE_val;
    ctx->params.U = U_val;
    ctx->params.default_particle_radius = ctx->params.sigma / 2.0;
    ctx->params.default_particle_mass = 1.0;
    ctx->current_time = 0.0;

    // Ensure particles are initialized within the pool
    for (int i = 0; i < num_particles; ++i) {
        ctx->particles[i].id = i;
        ctx->particles[i].radius = ctx->params.default_particle_radius;
        ctx->particles[i].radius_sq = ctx->particles[i].radius * ctx->particles[i].radius;
        ctx->particles[i].active = 1; // Mark as active
        ctx->particles[i].in_shoulder = 0;
        ctx->particles[i].shoulder_partner_idx = -1;
        ctx->particles[i].coll = 0;
        // Cell info might not be strictly necessary if not testing cell interactions directly
        ctx->particles[i].cellx = 0;
        ctx->particles[i].celly = 0;
        ctx->particles[i].next_in_cell = NULL;
    }
}

void test_elastic_collision_unequal_masses() {
    printf("Starting test_elastic_collision_unequal_masses...\n");
    SimContext ctx = {0};
    setup_test_context(&ctx, 2, 100.0, 0.0, 0.0); // deltaE = 0, U = 0

    Particle* p1 = &ctx.particles[0];
    Particle* p2 = &ctx.particles[1];

    // Setup: p1 (m=1) moving right, p2 (m=2) stationary
    // p1: x=0, y=0, vx=1, vy=0, m=1
    // p2: x=1, y=0, vx=0, vy=0, m=2 (collision at sigma/2 + sigma/2 = sigma = 1.0)
    p1->x = 0.0; p1->y = 0.0; p1->vx = 1.0; p1->vy = 0.0; p1->m = 1.0;
    p2->x = ctx.params.sigma - 1e-9; p2->y = 0.0; p2->vx = 0.0; p2->vy = 0.0; p2->m = 2.0;
    // Ensure radii are consistent with sigma for contact
    p1->radius = ctx.params.sigma / 2.0;
    p2->radius = ctx.params.sigma / 2.0;


    double initial_ke = calculate_total_kinetic_energy(&ctx);
    double initial_px = p1->m * p1->vx + p2->m * p2->vx;
    double initial_py = p1->m * p1->vy + p2->m * p2->vy;

    printf("  Initial state: P0(x=%.2f,y=%.2f,vx=%.2f,vy=%.2f,m=%.1f), P1(x=%.2f,y=%.2f,vx=%.2f,vy=%.2f,m=%.1f)\n",
           p1->x, p1->y, p1->vx, p1->vy, p1->m, p2->x, p2->y, p2->vx, p2->vy, p2->m);
    printf("  Initial KE: %.6f, Px: %.6f, Py: %.6f\n", initial_ke, initial_px, initial_py);

    // Create a dummy event for process_collision
    Event collision_event = {0};
    collision_event.p1_idx = p1->id;
    collision_event.p2_idx = p2->id;
    collision_event.type = EVENT_COLLISION; // Ensure type is set
    // process_collision doesn't use event.time, but good practice
    collision_event.time = ctx.current_time; 

    // Manually trigger collision processing
    // Need to ensure particles are at contact distance for process_collision logic
    // The `process_collision` function expects particles to be exactly at `sigma` distance.
    // We adjust p2->x slightly to simulate this, assuming p1 is at 0.
    // p2->x = p1->x + ctx.params.sigma; // This would be after advancing to collision time.
    // For this unit test, we assume they are already at contact.
    // The positions set above (p2->x = sigma - epsilon) are *just before* contact.
    // process_collision itself doesn't advance particles.
    
    // To correctly use process_collision, we should simulate them being at contact.
    // The function calculates relative velocity at approach.
    // Let's assume they are exactly at contact for the purpose of the test call.
    // The function internally calculates nx, ny based on current positions.
    // For a head-on collision along x-axis: nx = -1 if p1 is on left, p2 on right.
    // (p2->x - p1->x) / dist. If p1 at 0, p2 at 1, then (1-0)/1 = 1. So nx=1.
    // This means impulse is applied along positive x for p1, negative x for p2.

    process_collision(&ctx, &collision_event);

    double final_ke = calculate_total_kinetic_energy(&ctx);
    double final_px = p1->m * p1->vx + p2->m * p2->vx;
    double final_py = p1->m * p1->vy + p2->m * p2->vy;

    printf("  Final state: P0(vx=%.6f,vy=%.6f), P1(vx=%.6f,vy=%.6f)\n",
           p1->vx, p1->vy, p2->vx, p2->vy);
    printf("  Final KE: %.6f, Px: %.6f, Py: %.6f\n", final_ke, final_px, final_py);

    // Theoretical velocities after 1D elastic collision:
    // v1_final = (m1-m2)/(m1+m2)*v1_initial + (2m2)/(m1+m2)*v2_initial
    // v2_final = (2m1)/(m1+m2)*v1_initial + (m2-m1)/(m1+m2)*v2_initial
    // Here, v1_initial=1, v2_initial=0, m1=1, m2=2
    // v1_final = (1-2)/(1+2)*1 + 0 = -1/3 * 1 = -0.333333
    // v2_final = (2*1)/(1+2)*1 + 0 = 2/3 * 1 = 0.666667
    double expected_v1x = (p1->m - p2->m) / (p1->m + p2->m) * 1.0 + (2 * p2->m) / (p1->m + p2->m) * 0.0;
    double expected_v2x = (2 * p1->m) / (p1->m + p2->m) * 1.0 + (p2->m - p1->m) / (p1->m + p2->m) * 0.0;
    
    printf("  Expected v1x: %.6f, v2x: %.6f\n", expected_v1x, expected_v2x);

    assert(fabs(p1->vx - expected_v1x) < 1e-9);
    assert(fabs(p2->vx - expected_v2x) < 1e-9);
    assert(fabs(p1->vy - 0.0) < 1e-9); // No change in y-velocity
    assert(fabs(p2->vy - 0.0) < 1e-9); // No change in y-velocity
    assert(fabs(final_ke - initial_ke) < 1e-9); // Conserve KE
    assert(fabs(final_px - initial_px) < 1e-9); // Conserve Px
    assert(fabs(final_py - initial_py) < 1e-9); // Conserve Py

    free_sim_context(&ctx);
    printf("test_elastic_collision_unequal_masses PASSED\n\n");
}

void test_energy_injected_collision() {
    printf("Starting test_energy_injected_collision...\n");
    SimContext ctx = {0};
    double deltaE_val = 0.5; // Inject 0.5 units of energy
    setup_test_context(&ctx, 2, 100.0, deltaE_val, 0.0); // U = 0

    Particle* p1 = &ctx.particles[0];
    Particle* p2 = &ctx.particles[1];

    // Setup: p1 (m=1) moving right, p2 (m=1) stationary. Head-on collision.
    p1->x = 0.0; p1->y = 0.0; p1->vx = 1.0; p1->vy = 0.0; p1->m = 1.0;
    p2->x = ctx.params.sigma - 1e-12; p2->y = 0.0; p2->vx = 0.0; p2->vy = 0.0; p2->m = 1.0;
    p1->radius = ctx.params.sigma / 2.0;
    p2->radius = ctx.params.sigma / 2.0;


    double initial_ke = calculate_total_kinetic_energy(&ctx);
    double initial_px = p1->m * p1->vx + p2->m * p2->vx;
    printf("  Initial state: P0(x=%.2f,y=%.2f,vx=%.2f,vy=%.2f,m=%.1f), P1(x=%.2f,y=%.2f,vx=%.2f,vy=%.2f,m=%.1f)\n",
           p1->x, p1->y, p1->vx, p1->vy, p1->m, p2->x, p2->y, p2->vx, p2->vy, p2->m);
    printf("  Initial KE: %.6f, Px: %.6f, deltaE: %.6f\n", initial_ke, initial_px, deltaE_val);

    Event collision_event = {0};
    collision_event.p1_idx = p1->id;
    collision_event.p2_idx = p2->id;
    collision_event.type = EVENT_COLLISION;
    collision_event.time = ctx.current_time;

    process_collision(&ctx, &collision_event);

    double final_ke = calculate_total_kinetic_energy(&ctx);
    double final_px = p1->m * p1->vx + p2->m * p2->vx;

    printf("  Final state: P0(vx=%.6f,vy=%.6f), P1(vx=%.6f,vy=%.6f)\n",
           p1->vx, p1->vy, p2->vx, p2->vy);
    printf("  Final KE: %.6f, Px: %.6f\n", final_ke, final_px);

    // For equal masses, head-on, p1 hits stationary p2:
    // Elastic: p1 stops, p2 gets p1's velocity.
    // With deltaE:
    // mu = m1*m2/(m1+m2) = 1*1/(1+1) = 0.5
    // v_rel_approach = v1x - v2x = 1.0 - 0.0 = 1.0
    // e_eff = sqrt(1 + 2*deltaE / (mu * v_rel_approach^2))
    // e_eff = sqrt(1 + 2*0.5 / (0.5 * 1^2)) = sqrt(1 + 1 / 0.5) = sqrt(1 + 2) = sqrt(3)
    // J = mu * (1 + e_eff) * v_rel_approach_magnitude
    // J = 0.5 * (1 + sqrt(3)) * 1.0
    // v1x_final = v1x_initial - J/m1 = 1.0 - 0.5 * (1 + sqrt(3)) / 1.0 = 1.0 - 0.5 - 0.5*sqrt(3) = 0.5 - 0.5*sqrt(3)
    // v2x_final = v2x_initial + J/m2 = 0.0 + 0.5 * (1 + sqrt(3)) / 1.0 = 0.5 + 0.5*sqrt(3)
    double mu = (p1->m * p2->m) / (p1->m + p2->m);
    double v_rel_n_approach_sq = (p1->vx - p2->vx)*(p1->vx - p2->vx); // Assuming head on, nx=1, so v_rel_n = v_rel_x
    double e_eff = sqrt(1.0 + (2.0 * deltaE_val) / (mu * v_rel_n_approach_sq));
    double J = mu * (1.0 + e_eff) * fabs(p1->vx - p2->vx); // fabs for magnitude

    // Initial velocities for calculation: v1x_init = 1.0, v2x_init = 0.0
    double expected_v1x = 1.0 - J / p1->m; 
    double expected_v2x = 0.0 + J / p2->m;

    printf("  Expected v1x: %.6f, v2x: %.6f\n", expected_v1x, expected_v2x);

    assert(fabs(p1->vx - expected_v1x) < 1e-9);
    assert(fabs(p2->vx - expected_v2x) < 1e-9);
    assert(fabs(final_ke - (initial_ke + deltaE_val)) < 1e-9); // KE should increase by deltaE
    assert(fabs(final_px - initial_px) < 1e-9); // Momentum conserved

    free_sim_context(&ctx);
    printf("test_energy_injected_collision PASSED\n\n");
}

void test_shoulder_entry_reflection() {
    printf("Starting test_shoulder_entry_reflection...\n");
    SimContext ctx = {0};
    double U_val = 1.0; // Shoulder height
    setup_test_context(&ctx, 2, 100.0, 0.0, U_val); // deltaE = 0

    Particle* p1 = &ctx.particles[0];
    Particle* p2 = &ctx.particles[1];

    // Setup: p1 (m=1) moving towards stationary p2 (m=1).
    // KE_normal_approach < U for reflection.
    // Shoulder width = sigma * lambda_shoulder = 1.0 * 1.5 = 1.5
    // Particles start at distance > 1.5, e.g., 2.0.
    // p1: x=0, y=0, vx=0.5, vy=0, m=1
    // p2: x=1.5-eps, y=0, vx=0, vy=0, m=1 (at shoulder boundary)
    p1->x = 0.0; p1->y = 0.0; p1->vx = 0.5; p1->vy = 0.0; p1->m = 1.0;
    p2->x = ctx.params.sigma * ctx.params.lambda_shoulder - 1e-9; // Just at shoulder boundary for p1 to enter
    p2->y = 0.0; p2->vx = 0.0; p2->vy = 0.0; p2->m = 1.0;
    // Radii are sigma/2, shoulder boundary is at sigma * lambda_shoulder center-to-center
    p1->radius = ctx.params.sigma / 2.0;
    p2->radius = ctx.params.sigma / 2.0;


    // Relative velocity normal component: v_rel_n = (v1-v2) . n
    // Here, n is from p1 to p2, so n = (1,0) if p1 is left of p2.
    // v_rel_x = p1->vx - p2->vx = 0.5 - 0 = 0.5
    // mu = 0.5 for m1=m2=1.
    // KE_normal_approach = 0.5 * mu * v_rel_n^2 = 0.5 * 0.5 * (0.5)^2 = 0.25 * 0.25 = 0.0625
    // U_val = 1.0. So KE_normal_approach (0.0625) < U_val (1.0) -> reflection.
    double mu_shoulder = (p1->m * p2->m) / (p1->m + p2->m);
    double v_rel_x_approach = p1->vx - p2->vx; // Assuming collision along x
    double ke_normal_approach = 0.5 * mu_shoulder * v_rel_x_approach * v_rel_x_approach;
    printf("  U_val=%.6f, KE_normal_approach=%.6f (Should reflect)\n", U_val, ke_normal_approach);
    assert(ke_normal_approach < U_val);

    double initial_ke = calculate_total_kinetic_energy(&ctx);
    double initial_px = p1->m * p1->vx + p2->m * p2->vx;

    Event shoulder_event = {0};
    shoulder_event.p1_idx = p1->id;
    shoulder_event.p2_idx = p2->id;
    shoulder_event.type = EVENT_SHOULDER_ENTRY;
    shoulder_event.time = ctx.current_time; // Process immediately

    // Manually set particles to be at the shoulder boundary for handle_shoulder_entry
    // The function itself adjusts positions to be *exactly* at boundary.
    // For this test, we assume they are brought to the boundary by event time.
    // p1->x = p2->x - ctx.params.sigma * ctx.params.lambda_shoulder; // if p2 is at origin
    // Let p1 be at 0, p2 be at shoulder_radius
    p1->x = 0;
    p2->x = ctx.params.sigma * ctx.params.lambda_shoulder;


    handle_shoulder_entry(&ctx, &shoulder_event);

    double final_ke = calculate_total_kinetic_energy(&ctx);
    double final_px = p1->m * p1->vx + p2->m * p2->vx;

    printf("  Final state: P0(vx=%.6f,vy=%.6f, in_shoulder=%d), P1(vx=%.6f,vy=%.6f, in_shoulder=%d)\n",
           p1->vx, p1->vy, p1->in_shoulder, p2->vx, p2->vy, p2->in_shoulder);
    printf("  Final KE: %.6f, Px: %.6f\n", final_ke, final_px);

    // Elastic reflection: v1_final = -v1_initial, v2_final = v2_initial (if m1=m2, v2_initial=0)
    // For reflection, impulse J = -2 * mu * v_rel_n_approach
    // v_rel_n_approach = v1x - v2x = 0.5 (as calculated before, assuming nx points from p1 to p2, so v_rel . n > 0 if approaching)
    // The function uses v_rel_n = (v2-v1).n, so if p1 approaches p2 from left, n=(1,0), v_rel_n = (v2x-v1x)*1 = -0.5
    // Impulse J = -2 * mu * (-0.5) = mu. (This J is in the direction of n)
    // p1->vx_new = p1->vx_old + J_normal_component_for_p1 / m1
    // p2->vx_new = p2->vx_old - J_normal_component_for_p2 / m2
    // In handle_shoulder_entry, for reflection:
    // impulse_J = -2.0 * mu * v_rel_n; (where v_rel_n is (v2-v1).n, so it's negative if approaching)
    // p1->vx += impulse_J * nx / p1->m;
    // p2->vx -= impulse_J * nx / p2->m;
    // If p1 at left, p2 at right, nx = (p2x-p1x)/dist = 1.0
    // v_rel_n = (p2->vx - p1->vx)*nx + (p2->vy - p1->vy)*ny = (0 - 0.5)*1 = -0.5
    // impulse_J = -2.0 * 0.5 * (-0.5) = 0.5
    // p1->vx_new = 0.5 + 0.5 * 1.0 / 1.0 = 1.0  -- This is wrong.
    // The impulse in the code is defined as J = mu * (1+e) * v_rel_approach_abs.
    // For elastic reflection e=1. J = 2 * mu * v_rel_approach_abs.
    // v_rel_approach_abs = 0.5. J = 2 * 0.5 * 0.5 = 0.5.
    // p1_vx_final = p1_vx_initial - J/m1 = 0.5 - 0.5/1 = 0.0
    // p2_vx_final = p2_vx_initial + J/m2 = 0.0 + 0.5/1 = 0.5
    // This is for p1 hitting a wall.
    // For two particles, elastic head-on, equal mass: they swap velocities.
    // So p1 should have 0 vx, p2 should have 0.5 vx.
    // Let's re-check the `handle_shoulder_entry` reflection logic:
    // impulse_J = -2.0 * mu * v_rel_n; (v_rel_n is (v2-v1).n)
    // v_rel_n = (p2->vx - p1->vx)*nx = (0 - 0.5)*1 = -0.5 for nx=1 (p1 left, p2 right)
    // impulse_J = -2.0 * 0.5 * (-0.5) = 0.5.
    // p1->vx += impulse_J * nx / p1->m = 0.5 + 0.5*1/1 = 1.0. Still not right.
    // The velocities should be swapped for equal mass elastic collision.
    // p1.vx_final = 0, p2.vx_final = 0.5
    // The issue might be in my interpretation of `v_rel_n` or `nx` in the test vs code.
    // Code: `v_rel_n = (p2->vx - p1->vx) * nx + (p2->vy - p1->vy) * ny;`
    // `nx = (p2->x - p1->x) / dist;` `ny = (p2->y - p1->y) / dist;`
    // If p1=(0,0) p2=(1.5,0), then nx=1, ny=0.
    // v_rel_n = (0 - 0.5)*1 + (0-0)*0 = -0.5
    // impulse_J = -2.0 * mu_shoulder * v_rel_n = -2.0 * 0.5 * (-0.5) = 0.5.
    // p1->vx += impulse_J * nx / p1->m; => p1->vx = 0.5 + 0.5 * 1.0 / 1.0 = 1.0.
    // p2->vx -= impulse_J * nx / p2->m; => p2->vx = 0.0 - 0.5 * 1.0 / 1.0 = -0.5.
    // This means p1 speeds up, p2 moves left. This is not a simple reflection or swap.
    // Ah, the `handle_shoulder_entry` reflection is like particle hitting a fixed potential barrier represented by the other particle.
    // It's not a two-body elastic collision in the standard sense for this reflection part.
    // It's an elastic reflection of the *relative normal velocity*.
    // So, v_rel_n_final = -v_rel_n_initial.
    // (v2_f - v1_f).n = - (v2_i - v1_i).n
    // (v2x_f - v1x_f) = - (v2x_i - v1x_i) = - (0 - 0.5) = 0.5
    // And momentum conservation: m1*v1x_f + m2*v2x_f = m1*v1x_i + m2*v2x_i
    // v1x_f + v2x_f = 0.5 (since m1=m2=1)
    // Solving: v2x_f - v1x_f = 0.5 and v2x_f + v1x_f = 0.5
    // 2*v2x_f = 1.0 => v2x_f = 0.5
    // v1x_f = 0.0
    // So, expected_v1x = 0.0, expected_v2x = 0.5. This is velocity swap.

    assert(fabs(p1->vx - 0.0) < 1e-9); // p1 stops
    assert(fabs(p2->vx - 0.5) < 1e-9); // p2 moves with p1's initial velocity
    assert(p1->in_shoulder == 0); // Should not enter shoulder on reflection
    assert(p2->in_shoulder == 0);
    assert(fabs(final_ke - initial_ke) < 1e-9); // KE conserved in elastic reflection
    assert(fabs(final_px - initial_px) < 1e-9); // Momentum conserved

    free_sim_context(&ctx);
    printf("test_shoulder_entry_reflection PASSED\n\n");
}

void test_shoulder_entry_penetration() {
    printf("Starting test_shoulder_entry_penetration...\n");
    SimContext ctx = {0};
    double U_val = 0.1; // Shoulder height
    setup_test_context(&ctx, 2, 100.0, 0.0, U_val); // deltaE = 0

    Particle* p1 = &ctx.particles[0];
    Particle* p2 = &ctx.particles[1];

    // Setup: p1 (m=1) moving towards stationary p2 (m=1).
    // KE_normal_approach > U for penetration.
    // p1: x=0, y=0, vx=1.0, vy=0, m=1
    // p2: x=1.5-eps, y=0, vx=0, vy=0, m=1
    p1->x = 0.0; p1->y = 0.0; p1->vx = 1.0; p1->vy = 0.0; p1->m = 1.0;
    p2->x = ctx.params.sigma * ctx.params.lambda_shoulder - 1e-9; 
    p2->y = 0.0; p2->vx = 0.0; p2->vy = 0.0; p2->m = 1.0;
    p1->radius = ctx.params.sigma / 2.0;
    p2->radius = ctx.params.sigma / 2.0;

    double mu_shoulder = (p1->m * p2->m) / (p1->m + p2->m); // 0.5
    double v_rel_x_approach = p1->vx - p2->vx; // 1.0
    double ke_normal_approach_initial = 0.5 * mu_shoulder * v_rel_x_approach * v_rel_x_approach; // 0.5 * 0.5 * 1^2 = 0.25
    printf("  U_val=%.6f, KE_normal_approach_initial=%.6f (Should penetrate)\n", U_val, ke_normal_approach_initial);
    assert(ke_normal_approach_initial >= U_val);

    double initial_total_ke = calculate_total_kinetic_energy(&ctx); // KE of system = 0.5 * 1 * 1^2 = 0.5
    double initial_px = p1->m * p1->vx + p2->m * p2->vx; // 1.0

    Event shoulder_event = {0};
    shoulder_event.p1_idx = p1->id;
    shoulder_event.p2_idx = p2->id;
    shoulder_event.type = EVENT_SHOULDER_ENTRY;
    shoulder_event.time = ctx.current_time;

    // Set positions for handle_shoulder_entry
    p1->x = 0;
    p2->x = ctx.params.sigma * ctx.params.lambda_shoulder;

    handle_shoulder_entry(&ctx, &shoulder_event);

    double final_total_ke = calculate_total_kinetic_energy(&ctx);
    double final_px = p1->m * p1->vx + p2->m * p2->vx;

    printf("  Final state: P0(vx=%.6f,vy=%.6f, in_shoulder=%d), P1(vx=%.6f,vy=%.6f, in_shoulder=%d)\n",
           p1->vx, p1->vy, p1->in_shoulder, p2->vx, p2->vy, p2->in_shoulder);
    printf("  Final total KE: %.6f, Px: %.6f\n", final_total_ke, final_px);

    // Expected: Normal component of KE of relative motion is reduced by U.
    // KE_normal_final = KE_normal_initial - U
    // 0.5 * mu * (v_rel_n_final)^2 = 0.5 * mu * (v_rel_n_initial)^2 - U
    // (v_rel_n_final)^2 = (v_rel_n_initial)^2 - 2U/mu
    // v_rel_n_final = sqrt( (v_rel_n_initial)^2 - 2U/mu )
    // v_rel_x_initial = 1.0. (v_rel_x_initial)^2 = 1.0
    // v_rel_x_final_sq = 1.0 - 2*0.1/0.5 = 1.0 - 0.4 = 0.6
    // v_rel_x_final = sqrt(0.6)
    // v_rel_x_final = p1->vx_final - p2->vx_final (if still approaching, or could be separating if U is large)
    // The code sets (v2_f - v1_f) = sign((v2_i-v1_i)) * sqrt( |(v2_i-v1_i)|^2 - 2U/mu )
    // (v2x_f - v1x_f) = - sqrt( (v1x_i - v2x_i)^2 - 2U/mu ) = -sqrt(1.0 - 2*0.1/0.5) = -sqrt(0.6)
    // Momentum conservation: p1->vx_final + p2->vx_final = p1->vx_initial + p2->vx_initial = 1.0
    // Let v1f = p1->vx_final, v2f = p2->vx_final
    // v2f - v1f = -sqrt(0.6)
    // v2f + v1f = 1.0
    // Adding: 2*v2f = 1.0 - sqrt(0.6) => v2f = (1.0 - sqrt(0.6))/2.0
    // Subtracting: -2*v1f = -sqrt(0.6) - 1.0 => v1f = (1.0 + sqrt(0.6))/2.0
    double expected_v_rel_x_final_sq = v_rel_x_approach * v_rel_x_approach - 2 * U_val / mu_shoulder;
    double expected_v_rel_x_final = -sqrt(expected_v_rel_x_final_sq); // Negative as still approaching or just passed through

    double expected_v1x = ( (p1->vx + p2->vx) - expected_v_rel_x_final ) / 2.0; // (Sum_v - Vrel_final)/2
    double expected_v2x = ( (p1->vx + p2->vx) + expected_v_rel_x_final ) / 2.0; // (Sum_v + Vrel_final)/2
    // Re-deriving from code logic:
    // v_rel_n_new_magnitude = sqrt(v_rel_n_approach_sq - (2.0 * ctx->params.U) / mu);
    // delta_v_rel_n = v_rel_n_new_magnitude - v_rel_n_approach_magnitude; (v_rel_n_approach_magnitude is positive)
    // impulse_J_normal = mu * delta_v_rel_n; (This impulse is to change v_rel_n)
    // This is not how the code does it. Code:
    // v_rel_n_final_magnitude = sqrt(KE_normal_approach * 2.0 / mu - 2.0 * U / mu) -- this is wrong.
    // It should be: v_rel_n_final_magnitude = sqrt( (v_rel_n_initial_magnitude^2) - 2*U/mu )
    // The code calculates new v1,v2 based on changing KE_normal by -U.
    // Total KE of system should be initial_total_ke - U
    // (The U is "lost" to potential energy, so kinetic energy decreases by U)

    printf("  Expected v1x: %.6f, v2x: %.6f\n", expected_v1x, expected_v2x);
    assert(fabs(p1->vx - expected_v1x) < 1e-9);
    assert(fabs(p2->vx - expected_v2x) < 1e-9);
    assert(p1->in_shoulder == 1);
    assert(p2->in_shoulder == 1);
    assert(fabs(final_total_ke - (initial_total_ke - U_val)) < 1e-9); // Total KE reduced by U
    assert(fabs(final_px - initial_px) < 1e-9); // Momentum conserved

    free_sim_context(&ctx);
    printf("test_shoulder_entry_penetration PASSED\n\n");
}

void test_shoulder_exit() {
    printf("Starting test_shoulder_exit...\n");
    SimContext ctx = {0};
    double U_val = 0.5; // Shoulder height
    setup_test_context(&ctx, 2, 100.0, 0.0, U_val); // deltaE = 0

    Particle* p1 = &ctx.particles[0];
    Particle* p2 = &ctx.particles[1];

    // Setup: Particles are inside shoulder and moving apart.
    // p1: x=0.6, y=0, vx=-0.5, vy=0, m=1, in_shoulder=1
    // p2: x=2.0, y=0, vx=0.5, vy=0, m=1, in_shoulder=1
    // Shoulder boundary is 1.5. They are exiting.
    // Let's place them exactly at the shoulder boundary, moving apart.
    p1->x = 0.0; p1->y = 0.0; p1->vx = -0.5; p1->vy = 0.0; p1->m = 1.0;
    p2->x = ctx.params.sigma * ctx.params.lambda_shoulder; 
    p2->y = 0.0; p2->vx = 0.5; p2->vy = 0.0; p2->m = 1.0;
    
    p1->in_shoulder = 1; p1->shoulder_partner_idx = p2->id;
    p2->in_shoulder = 1; p2->shoulder_partner_idx = p1->id;
    p1->radius = ctx.params.sigma / 2.0;
    p2->radius = ctx.params.sigma / 2.0;


    double initial_total_ke = calculate_total_kinetic_energy(&ctx); // 0.5*1*(-0.5)^2 + 0.5*1*(0.5)^2 = 0.125 + 0.125 = 0.25
    double initial_px = p1->m * p1->vx + p2->m * p2->vx; // 0.0

    printf("  Initial state: P0(vx=%.2f, in_s=%d), P1(vx=%.2f, in_s=%d)\n", p1->vx, p1->in_shoulder, p2->vx, p2->in_shoulder);
    printf("  Initial total KE: %.6f, Px: %.6f, U_val: %.6f\n", initial_total_ke, initial_px, U_val);

    Event shoulder_event = {0};
    shoulder_event.p1_idx = p1->id;
    shoulder_event.p2_idx = p2->id;
    shoulder_event.type = EVENT_SHOULDER_EXIT;
    shoulder_event.time = ctx.current_time;

    process_shoulder_exit(&ctx, &shoulder_event);

    double final_total_ke = calculate_total_kinetic_energy(&ctx);
    double final_px = p1->m * p1->vx + p2->m * p2->vx;

    printf("  Final state: P0(vx=%.6f,vy=%.6f, in_shoulder=%d), P1(vx=%.6f,vy=%.6f, in_shoulder=%d)\n",
           p1->vx, p1->vy, p1->in_shoulder, p2->vx, p2->vy, p2->in_shoulder);
    printf("  Final total KE: %.6f, Px: %.6f\n", final_total_ke, final_px);

    // Expected: Normal component of KE of relative motion is increased by U.
    // KE_normal_final = KE_normal_initial + U
    // 0.5 * mu * (v_rel_n_final)^2 = 0.5 * mu * (v_rel_n_initial)^2 + U
    // (v_rel_n_final)^2 = (v_rel_n_initial)^2 + 2U/mu
    // v_rel_n_final = sqrt( (v_rel_n_initial)^2 + 2U/mu )
    // v_rel_x_initial = -1.0. (v_rel_x_initial)^2 = 1.0
    // v_rel_x_final_sq = 1.0 + 2*0.5/0.5 = 1.0 + 2.0 = 3.0
    // v_rel_x_final = sqrt(3.0)
    // Since they are separating, (v1_f - v2_f) should be -sqrt(3.0)
    // Or (v2_f - v1_f) = sqrt(3.0)
    // Momentum conservation: p1->vx_final + p2->vx_final = initial_px = 0.0
    // So, v1f + v2f = 0 => v1f = -v2f
    // v2f - (-v2f) = sqrt(3.0) => 2*v2f = sqrt(3.0) => v2f = sqrt(3.0)/2.0
    // v1f = -sqrt(3.0)/2.0
    double expected_v1x = -sqrt(3.0) / 2.0;
    double expected_v2x = sqrt(3.0) / 2.0;
    
    printf("  Expected v1x: %.6f, v2x: %.6f\n", expected_v1x, expected_v2x);

    assert(fabs(p1->vx - expected_v1x) < 1e-9);
    assert(fabs(p2->vx - expected_v2x) < 1e-9);
    assert(p1->in_shoulder == 0);
    assert(p2->in_shoulder == 0);
    // Total KE of system should be initial_total_ke + U (as potential energy is converted back to KE)
    assert(fabs(final_total_ke - (initial_total_ke + U_val)) < 1e-9);
    assert(fabs(final_px - initial_px) < 1e-9); // Momentum conserved

    free_sim_context(&ctx);
    printf("test_shoulder_exit PASSED\n\n");
}


int main() {
    test_elastic_collision_unequal_masses();
    test_energy_injected_collision();
    test_shoulder_entry_reflection();
    test_shoulder_entry_penetration();
    test_shoulder_exit();
    // Add calls to other tests here
    return 0;
}

