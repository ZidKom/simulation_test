#include "sim_context.h"
#include <stdlib.h>
#include <string.h>
#include <xmmintrin.h>
#include <float.h> // For DBL_MAX
#include <stdio.h> // For fprintf, exit, EXIT_FAILURE

void init_particle_pool(SimContext* ctx, int N) {
    if (!ctx || N <= 0) return;
    ctx->particles = (Particle*)_mm_malloc(N * sizeof(Particle), 64);
    ctx->num_particles = N;
    for (int i = 0; i < N; ++i) {
        Particle* p = &ctx->particles[i];
        p->x = p->y = p->vx = p->vy = 0.0;
        p->radius = 0.0;
        p->cellx = p->celly = -1;
        p->active = 0;
    }
}

void release_particle(SimContext* ctx, Particle* p) {
    (void)ctx; // Unused parameter
    if(!p || !p->active) return;
    p->active = 0;
    // Optionally: zero memory for safety
    memset(p, 0, sizeof(Particle));
}

// Event pool management (see md_event_functions.c:26-53)
int init_event_pool(SimContext *ctx, int pool_size_param) {
    if (!ctx) return -1;

    // New heuristic based on num_particles if pool_size_param is a suggestion or 0
    // If pool_size_param is explicitly given (e.g. from config), use it.
    // The prompt implies event_pool_sz_val = ctx.num_particles * (ctx.num_particles + 5) * 20;
    // This suggests init_event_pool is called *after* ctx->num_particles is known.
    
    int actual_pool_size;
    if (pool_size_param > 0) { // If a specific size is requested
        actual_pool_size = pool_size_param;
    } else { // Use heuristic based on num_particles
        if (ctx->num_particles > 0) {
            actual_pool_size = ctx->num_particles * (ctx->num_particles + 5) * 20;
        } else {
            actual_pool_size = 1000; // Default if num_particles not yet set or is 0
        }
    }
    if (actual_pool_size < 1000) actual_pool_size = 1000; // Ensure a minimum size

    ctx->event_pool_storage = (Event*)aligned_alloc(64, actual_pool_size * sizeof(Event));
    if (!ctx->event_pool_storage) {
        if (ctx->hooks.log_error) ctx->hooks.log_error("Failed to allocate event pool storage.");
        ctx->event_pool_size = 0;
        return -1;
    }
    ctx->event_pool_size = actual_pool_size;
    memset(ctx->event_pool_storage, 0, actual_pool_size * sizeof(Event));

    ctx->free_event_list = NULL;
    for (int i = 0; i < ctx->event_pool_size; ++i) {
        ctx->event_pool_storage[i].next = ctx->free_event_list;
        ctx->free_event_list = &ctx->event_pool_storage[i];
    }
    return 0;
}

Event* allocate_event_from_pool(SimContext *ctx) {
    if (!ctx) return NULL;

    if (!ctx->free_event_list) {
        // Pool is exhausted, try to resize
        size_t old_size = ctx->event_pool_size;
        size_t new_size = (old_size > 0) ? (old_size * 2) : 1000; // Double the size, or start with 1000
        
        if (ctx->hooks.log_warning) ctx->hooks.log_warning("Event pool exhausted. Attempting to resize from %zu to %zu.", old_size, new_size);

        Event* new_pool_storage = (Event*)realloc(ctx->event_pool_storage, new_size * sizeof(Event));
        if (!new_pool_storage) {
            if (ctx->hooks.log_error) ctx->hooks.log_error("Failed to reallocate event pool to size %zu.", new_size);
            return NULL; // Out of memory
        }
        ctx->event_pool_storage = new_pool_storage;

        // Link the new chunk of events into the free list
        // The new events are from index old_size to new_size - 1
        // Important: The existing free_event_list is now invalid if realloc moved the memory.
        // We should rebuild it for the new part.
        // The prompt's fix:
        // for (int i = ctx->event_pool_size; i < new_size; ++i) { // This should be size_t
        //    new_pool[i].next_free = ctx->free_event_list; // next_free is 'next' in current Event struct
        //    ctx->free_event_list = &new_pool[i];
        // }
        // This prepends the new items.

        // Initialize and link the newly allocated part
        ctx->free_event_list = NULL; // Reset and build, or carefully append. Prompt implies reset for new part.
                                     // To be safe, let's assume the old free_event_list is invalid.
                                     // We are adding new_size - old_size new events.
        for (size_t i = old_size; i < new_size; ++i) {
            memset(&ctx->event_pool_storage[i], 0, sizeof(Event)); // Initialize new events
            ctx->event_pool_storage[i].next = ctx->free_event_list;
            ctx->free_event_list = &ctx->event_pool_storage[i];
        }
        
        ctx->event_pool_size = new_size;

        if (!ctx->free_event_list) { // Should not happen if resize succeeded and old_size < new_size
             if (ctx->hooks.log_error) ctx->hooks.log_error("Event pool still exhausted after resize attempt.");
            return NULL; 
        }
    }

    Event* new_event = ctx->free_event_list;
    ctx->free_event_list = new_event->next; // 'next' is the free list pointer in struct Event
    
    // Initialize the event (memset already done for new parts during resize)
    // For events taken from existing free list, ensure they are clean.
    // The memset in schedule_event or a dedicated init_event function is better.
    // For now, let's ensure it's clean here.
    memset(new_event, 0, sizeof(Event)); 
    // new_event->id = ctx->next_event_id++; // Assuming SimContext has next_event_id
    return new_event;
}

void release_event_to_pool(SimContext *ctx, Event *event) {
    if (!ctx || !event) {
        // fprintf(stderr, "Warning: Attempted to release a NULL event or with NULL context.\n");
        return;
    }

    // Add the event to the head of the free list
    event->next = ctx->free_event_list;
    ctx->free_event_list = event;
}

void free_event_pool(SimContext *ctx) {
    if (ctx) {
        // The actual memory (ctx->event_pool_storage) is freed in free_sim_context.
        // Here, we just null out the list head as a good measure, though it's redundant
        // if free_sim_context is called right after.
        ctx->free_event_list = NULL;
        // ctx->event_pool_storage will be freed by the caller (free_sim_context)
    }
}
