#include "../core/sim_context.h"
#include "../core/event_system/hybrid.h"
#include "../physics/cell_interactions.h"
#include "../physics/collisions.h"
#include "../utils/pbc.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

// Dummy validation hooks
void dummy_validate_collision(SimContext* ctx, Event* ev) { (void)ctx; (void)ev; }
void dummy_validate_energy(SimContext* ctx) { (void)ctx; }
void dummy_validate_momentum(SimContext* ctx) { (void)ctx; }
void dummy_log_error(const char* fmt, ...) { (void)fmt; }

void test_cell_crossing() {
    printf("Running cell crossing test...\n");
    SimContext ctx = {0};
    ctx.xsize = ctx.ysize = 10.0;
    ctx.cell_size = 2.0;
    ctx.n_cells_x = ctx.n_cells_y = 5;
    ctx.num_particles = 1;
    ctx.particles = malloc(sizeof(Particle));
    ctx.particles[0].x = 1.0;
    ctx.particles[0].y = 1.0;
    ctx.particles[0].vx = 1.0;
    ctx.particles[0].vy = 0.0;
    ctx.particles[0].cellx = 0;
    ctx.particles[0].celly = 0;
    ctx.particles[0].id = 0;
    ctx.cells = calloc(ctx.n_cells_x, sizeof(Cell*));
    for (int i = 0; i < ctx.n_cells_x; i++) {
        ctx.cells[i] = calloc(ctx.n_cells_y, sizeof(Cell));
    }
    
    // Ensure event pool is initialized with a large enough size for tests
    int event_pool_size = 10000;
    init_sim_context(&ctx, ctx.num_particles, ctx.xsize, ctx.ysize, ctx.cell_size, ctx.event_system.paul_list_size, event_pool_size);

    // Test cell crossing event prediction
    predict_cell_crossing(&ctx, &ctx.particles[0]);
    
    // Verify the crossing time is correct (should be 1.0 seconds)
    // Process cell crossing event (X+)
    Event ev = {.type = EVENT_CELL_CROSS_X_POS, .p1_idx = 0};
    handle_cell_crossing(&ctx, &ev);
    assert(ctx.particles[0].cellx == 1);
    printf("Cell crossing test passed.\n");
    
    // Clean up
    for (int i = 0; i < ctx.n_cells_x; i++) {
        free(ctx.cells[i]);
    }
    free(ctx.cells);
    free(ctx.particles);
}

void test_energy_injection() {
    printf("Running energy injection test...\n");
    SimContext ctx = {0};
    ctx.num_particles = 2;
    ctx.particles = malloc(2 * sizeof(Particle));
    ctx.particles[0].x = 0.0; ctx.particles[0].y = 0.0;
    ctx.particles[0].vx = 1.0; ctx.particles[0].vy = 0.0;
    ctx.particles[0].radius = 0.5;
    ctx.particles[0].m = 1.0;
    ctx.particles[0].id = 0;
    ctx.particles[1].x = 2.0; ctx.particles[1].y = 0.0;
    ctx.particles[1].vx = -1.0; ctx.particles[1].vy = 0.0;
    ctx.particles[1].radius = 0.5;
    ctx.particles[1].m = 1.0;
    ctx.particles[1].id = 1;
    ctx.params.U = 0.0;
    ctx.params.sigma = 0.0;
    ctx.params.deltaE = 0.1;
    ctx.current_time = 0.0;

    double ke0 = 0.5 * ctx.particles[0].m * (ctx.particles[0].vx * ctx.particles[0].vx) +
                 0.5 * ctx.particles[1].m * (ctx.particles[1].vx * ctx.particles[1].vx);

    for (int i = 0; i < 1000; ++i) {
        Event ev = {.type = EVENT_COLLISION, .p1_idx = 0, .p2_idx = 1};
        process_collision(&ctx, &ev);
    }
    
    double ke1 = 0.5 * ctx.particles[0].m * (ctx.particles[0].vx * ctx.particles[0].vx) +
                 0.5 * ctx.particles[1].m * (ctx.particles[1].vx * ctx.particles[1].vx);

    double expected = ke0 + 1000 * ctx.params.deltaE;
    assert(fabs(ke1 - expected) < 1e-6);
    printf("Energy injection test passed.\n");
    free(ctx.particles);
}

int main(void) {
    test_cell_crossing();
    test_energy_injection();
    printf("All event tests passed.\n");
    return 0;
}
