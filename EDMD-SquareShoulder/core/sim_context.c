#include "sim_context.h"
#include "particle_pool.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <stdint.h>
#include <math.h> // For fmin
#include "../utils/cell_list.h" // For init_cell_lists declaration
                             

#if defined(_MSC_VER) || defined(__MINGW32__)
#include <malloc.h> // For _aligned_malloc and _aligned_free
#else
#include <stdalign.h>  // (though not directly used by aligned_calloc here)
#endif

// Helper for aligned allocation
static void* aligned_calloc_internal(size_t n, size_t size, size_t alignment) {
    void* ptr = NULL;
    size_t total_size = n * size;
    if (total_size == 0) return NULL; // Avoid issues with zero-size allocations
    if (size > 0 && n > SIZE_MAX / size) {
        fprintf(stderr, "Overflow detected in aligned_calloc_internal: n=%zu, size=%zu\n", n, size);
        return NULL;
    }

#if defined(_MSC_VER) || defined(__MINGW32__)
    ptr = _aligned_malloc(total_size, alignment);
    if (ptr) {
        memset(ptr, 0, total_size);
    }
#else
    // posix_memalign requires alignment to be a power of two multiple of sizeof(void*)
    // and total_size to be a multiple of alignment.
    // For simplicity, we ensure alignment is at least sizeof(void*) and a power of 2.
    // A robust solution would handle size not being a multiple of alignment if required by specific posix_memalign versions.
    if (alignment < sizeof(void*)) alignment = sizeof(void*);
    // Ensure alignment is a power of two (simplified check)
    if ((alignment & (alignment - 1)) != 0) {
        // Find next power of two if not already (crude way)
        size_t p2 = 1;
        while(p2 < alignment) p2 <<= 1;
        alignment = p2;
    }
    
    if (posix_memalign(&ptr, alignment, total_size) == 0) {
        memset(ptr, 0, total_size);
    } else {
        ptr = NULL;
    }
#endif
    return ptr;
}

static void aligned_free_internal(void *ptr) {
    if (!ptr) return;
#if defined(_MSC_VER) || defined(__MINGW32__)
    _aligned_free(ptr);
#else
    free(ptr);
#endif
}


void init_sim_context(SimContext *ctx, int num_particles_val, double x_size_val, double y_size_val, 
                     double cell_size_val_param, int paul_list_size_val, int event_pool_sz_val) {
    if (!ctx) return;
    memset(ctx, 0, sizeof(SimContext));

    ctx->num_particles = num_particles_val;
    ctx->xsize = x_size_val;
    ctx->ysize = y_size_val;

    // Default physics parameters (can be overridden by config file later)
    ctx->params.sigma = 1.0; // Default core diameter
    ctx->params.lambda_shoulder = 1.5; // Default shoulder width factor (e.g., shoulder is 1.5x core)
    ctx->params.U = 1.0; // Default shoulder height
    ctx->params.deltaE = 0.0; // Default energy injection
    ctx->params.default_particle_radius = ctx->params.sigma / 2.0;
    ctx->params.default_particle_mass = 1.0;

    // Cell size calculation based on max interaction range
    // Max interaction range: sigma_core * lambda_shoulder (diameter of shoulder)
    // Or, if using prompt's formula: sigma_core * (1 + lambda_shoulder)
    // Using prompt's formula:
    double max_interaction_range = ctx->params.sigma * (1.0 + ctx->params.lambda_shoulder); 
                                     // If sigma is diameter, this is diameter * (1+factor)
                                     // If lambda_shoulder is the ratio of shoulder diameter to core diameter,
                                     // then max_interaction_range = ctx->params.sigma * ctx->params.lambda_shoulder.
                                     // Sticking to prompt: sigma * (1 + lambda_shoulder)
    
    // Cell size should be at least the maximum interaction range to ensure
    // interacting particles are in adjacent cells.
    ctx->cell_size = max_interaction_range; 
    // If cell_size_val_param was provided and is larger, it could be used, but typically it's derived.
    // For simplicity, we derive it here. The parameter cell_size_val_param is ignored for this fix.

    if (ctx->cell_size <= 1e-9) { // Prevent division by zero if range is tiny
        ctx->cell_size = fmin(ctx->xsize, ctx->ysize) / 10.0; // Fallback
        if (ctx->cell_size <= 1e-9) ctx->cell_size = 1.0;
    }
    
    ctx->n_cells_x = (int)(ctx->xsize / ctx->cell_size);
    ctx->n_cells_y = (int)(ctx->ysize / ctx->cell_size);

    if (ctx->n_cells_x < 1) ctx->n_cells_x = 1;
    if (ctx->n_cells_y < 1) ctx->n_cells_y = 1;
    
    // Adjust cell_size to make cells fit perfectly into the box dimensions
    // This can be important for some PBC and cell logic.
    // However, the prompt implies setting cell_size first, then n_cells.
    // If n_cells_x/y are calculated by floor(box/cell_size), then cell_size might not divide box evenly.
    // For now, following the prompt's implied order.

    init_cell_lists(ctx, ctx->n_cells_x, ctx->n_cells_y);

    // Initialize particles
    ctx->particles = (Particle*)aligned_calloc_internal(ctx->num_particles, sizeof(Particle), 64);
    if (!ctx->particles && ctx->num_particles > 0) {
        fprintf(stderr, "[EDMD] Failed to allocate particle memory.\n");
        exit(EXIT_FAILURE);
    }
    
    // Initialize particles with basic values
    for(int i=0; i < ctx->num_particles; ++i) {
        ctx->particles[i].id = i;
        ctx->particles[i].radius_sq = 0.0; // Will be set properly when radius is assigned
        ctx->particles[i].next_in_cell = NULL;
    }

    // Initialize event system components
    ctx->event_system.paul_list_size = paul_list_size_val;
    ctx->event_system.paul_dt = 0.1; // Default value, should be configurable
    ctx->event_system.current_paul_bin = 0;
    
    // Allocate array of pointers to Event (heads of Paul bins)
    ctx->event_system.paul_lists = (Event**)aligned_calloc_internal(ctx->event_system.paul_list_size, sizeof(Event*), 64); 
    if (!ctx->event_system.paul_lists && ctx->event_system.paul_list_size > 0) {
        fprintf(stderr, "[EDMD] Failed to allocate Paul lists.\n");
        aligned_free_internal(ctx->particles);
        exit(EXIT_FAILURE);
    }
    // After allocating ctx->event_system.paul_lists, loop through and initialize each element to NULL
    if (ctx->event_system.paul_lists) {
        for (int i = 0; i < ctx->event_system.paul_list_size; ++i) {
            ctx->event_system.paul_lists[i] = NULL;
        }
    }
    
    // BST root initialization
    ctx->event_system.event_tree_root = NULL;

    // Initialize cell lists
    ctx->cell_size = cell_size_val_param;
    ctx->n_cells_x = (int)(ctx->xsize / ctx->cell_size);
    ctx->n_cells_y = (int)(ctx->ysize / ctx->cell_size);
    
    if (ctx->n_cells_x == 0) ctx->n_cells_x = 1; // Ensure at least one cell
    if (ctx->n_cells_y == 0) ctx->n_cells_y = 1;

    ctx->cells = (Cell**)aligned_calloc_internal(ctx->n_cells_x, sizeof(Cell*), 64);
    if (!ctx->cells && ctx->n_cells_x > 0) {
        fprintf(stderr, "[EDMD] Failed to allocate cell list (rows).\n");
        aligned_free_internal(ctx->particles);
        aligned_free_internal(ctx->event_system.paul_lists);
        exit(EXIT_FAILURE);
    }
    
    for(int i = 0; i < ctx->n_cells_x; i++) {
        ctx->cells[i] = (Cell*)aligned_calloc_internal(ctx->n_cells_y, sizeof(Cell), 64);
        if (!ctx->cells[i] && ctx->n_cells_y > 0) {
            fprintf(stderr, "[EDMD] Failed to allocate cell list (cols for row %d).\n", i);
            // Free previously allocated rows
            for(int k=0; k<i; ++k) aligned_free_internal(ctx->cells[k]);
            aligned_free_internal(ctx->cells);
            aligned_free_internal(ctx->particles);
            aligned_free_internal(ctx->event_system.paul_lists);
            exit(EXIT_FAILURE);
        }
        // All cell heads are already zeroed by aligned_calloc_internal
    }
    
    // Initialize event pool with 10x safety margin
    // ctx->event_pool_size = event_pool_sz_val * 10; // Previous logic, potentially too large or small
    ctx->event_pool_size = event_pool_sz_val; // Use the carefully calculated event_pool_sz_val directly
    
    // Ensure ctx->event_pool_storage is allocated before calling init_event_pool
    ctx->event_pool_storage = (Event*)aligned_calloc_internal(ctx->event_pool_size, sizeof(Event), 64);
    if (!ctx->event_pool_storage && ctx->event_pool_size > 0) {
        fprintf(stderr, "[EDMD] Failed to allocate event pool. Requested size: %d\n", ctx->event_pool_size);
        for(int i=0; i < ctx->n_cells_x; ++i) aligned_free_internal(ctx->cells[i]);
        aligned_free_internal(ctx->cells);
        aligned_free_internal(ctx->particles);
        aligned_free_internal(ctx->event_system.paul_lists);
        exit(EXIT_FAILURE);
    }

    // Initialize event pool after allocation
    if (ctx->event_pool_storage) {
        // init_event_pool is in core/particle_pool.c and should be called after event_pool_storage is allocated.
        init_event_pool(ctx, ctx->event_pool_size); 
    } else if (ctx->event_pool_size > 0) {
        // This case should ideally be caught by the allocation check above, but as a safeguard:
        fprintf(stderr, "[EDMD] Event pool storage is NULL despite positive requested size: %d\n", ctx->event_pool_size);
        for(int i=0; i < ctx->n_cells_x; ++i) aligned_free_internal(ctx->cells[i]);
        aligned_free_internal(ctx->cells);
        aligned_free_internal(ctx->particles);
        aligned_free_internal(ctx->event_system.paul_lists);
        exit(EXIT_FAILURE);
    }
    
    // Default physics parameters
    ctx->params.U = 0.5;
    ctx->params.sigma = 0.2;
    ctx->params.lambda_shoulder = 1.5; // Default shoulder width factor
    ctx->params.deltaE = 0.01;
    ctx->params.default_particle_radius = 0.5;
    ctx->params.default_particle_mass = 1.0;
    
    // Energy tracking initialization
    ctx->energy.tracked = 0.0;
    ctx->energy.tolerance = 1e-6;

    // Log initialization
    ctx->log_file = NULL;
    ctx->log_verbosity = 0;

    // Validation hooks to NULL by default
    ctx->hooks.validate_collision = NULL;
    ctx->hooks.log_error = NULL;
    ctx->hooks.on_cell_crossing = NULL;
    
    // Default simulation control
    ctx->sim_end_time = 10.0; // Default value

    // Use box sizes from file if already set by particle loader
    if(ctx->xsize == 0 || ctx->ysize == 0) { 
        ctx->xsize = x_size_val;
        ctx->ysize = y_size_val;
    }
    
    // Cell list initialization after box size is set
    // Ensure params are initialized before this, especially default_particle_radius and sigma
    if (ctx->params.default_particle_radius == 0) {
        // Provide a default or handle error if params not set
        // For now, let's assume it might be set later or use a placeholder if critical for init
        // Or, ensure init_simulation_parameters(&ctx->params); is called before init_sim_context.
    }
    
    // Check for zero cell_size to prevent division by zero if params are not ready
    if (ctx->params.default_particle_radius > 0 && (1 + ctx->params.sigma) > 0) {
        ctx->cell_size = 2.0 * ctx->params.default_particle_radius * (1 + ctx->params.sigma);
        if (ctx->cell_size > 0) {
            ctx->n_cells_x = (int)(ctx->xsize / ctx->cell_size);
            ctx->n_cells_y = (int)(ctx->ysize / ctx->cell_size);

            // Ensure at least one cell in each dimension
            if (ctx->n_cells_x == 0) ctx->n_cells_x = 1;
            if (ctx->n_cells_y == 0) ctx->n_cells_y = 1;
        } else {
            // Handle error: cell_size is not positive
            // Fallback or error logging
            ctx->n_cells_x = 1; // Fallback
            ctx->n_cells_y = 1; // Fallback
        }
    } else {
        // Handle error: params not properly initialized for cell calculation
        // Fallback or error logging
        ctx->n_cells_x = 1; // Fallback
        ctx->n_cells_y = 1; // Fallback
        ctx->cell_size = ctx->xsize; // Fallback to a single cell if dimensions are known
        if (ctx->ysize < ctx->cell_size && ctx->ysize > 0) ctx->cell_size = ctx->ysize;
        if (ctx->cell_size == 0 && ctx->xsize > 0) ctx->cell_size = ctx->xsize; // if ysize was 0
        else if (ctx->cell_size == 0 && ctx->ysize > 0) ctx->cell_size = ctx->ysize; // if xsize was 0
        else if (ctx->cell_size == 0) ctx->cell_size = 1.0; // Absolute fallback
    }
}

void free_sim_context(SimContext *ctx) {
    if (!ctx) return;

    // Free particles array
    aligned_free_internal(ctx->particles);
    ctx->particles = NULL;

    // Free paul lists array (not the individual events)
    if (ctx->event_system.paul_lists) {
        aligned_free_internal(ctx->event_system.paul_lists);
        ctx->event_system.paul_lists = NULL;
    }

    // Free cell lists
    if (ctx->cells) {
        for (int i = 0; i < ctx->n_cells_x; i++) {
            free(ctx->cells[i]);
        }
        free(ctx->cells);
        ctx->cells = NULL;
    }
    // Free event pool (which contains all events)
    aligned_free_internal(ctx->event_pool_storage);
    ctx->event_pool_storage = NULL;
    ctx->free_event_list = NULL; // No longer valid after pool is freed
    ctx->event_system.event_tree_root = NULL; // Events in BST are from the pool

    // Close log file if open
    if (ctx->log_file) {
        fclose(ctx->log_file);
        ctx->log_file = NULL;
    }
    
    // Clear struct to prevent use after free
    memset(ctx, 0, sizeof(SimContext));
}
    
// Definitions for aligned_alloc and aligned_free if they were declared in sim_context.h
// and not defined as static. For now, they are static helpers.
// void* aligned_alloc(size_t alignment, size_t size) {
//     return aligned_calloc_internal(1, size, alignment); // Assuming n=1 for a single block
// }
// void aligned_free(void *ptr) {
//    aligned_free_internal(ptr);
// }
