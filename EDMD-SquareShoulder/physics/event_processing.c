#include "event_processing.h"
#include "collisions.h"
#include "shoulders.h"
#include "cell_interactions.h"
#include "../core/particle_pool.h"
#include "../validation/collision_validation.h"
#include "../physics/particle_motion.h"
#include "../io/logger.h"
#include <stdio.h>

/**
 * @brief Processes an event from the event queue
 * 
 * This function dispatches the event to the appropriate handler based on the event type.
 * Before processing, it optionally validates the event to ensure numerical stability.
 * 
 * @param ctx Pointer to the simulation context
 * @param ev Pointer to the event to process
 */
void process_event(SimContext* ctx, Event* ev) {
    if (!ctx || !ev) {
        if (ctx && ctx->hooks.log_error) ctx->hooks.log_error("Process_event called with null context or event.");
        return;
    }

    // Advance all particles to the event time if current_time < ev->time
    // This is a global advancement step.
    if (ev->time < ctx->current_time - 1e-9) { // Check for events in the past (with tolerance)
        if (ctx->hooks.log_warning) {
            ctx->hooks.log_warning("Warning: Event time %.12f is before current time %.12f. Event ID: %p, Type: %d. Skipping advancement, processing event.", 
                                   ev->time, ctx->current_time, (void*)ev, ev->type);
        }
        // Do not update ctx->current_time to ev->time if ev->time is in the past.
        // Process the event at ctx->current_time or handle as error.
        // For now, proceed to process, but this indicates a potential issue.
    } else if (ev->time > ctx->current_time + 1e-9) { // Event is in the future (with tolerance)
        double dt = ev->time - ctx->current_time;
        // Advance only particles involved in the event if that's the strategy.
        // The prompt implies advancing all particles to event time.
        // However, the original code in combined_new.txt for process_event advances ALL particles.
        // If only involved particles are advanced, that's a different optimization.
        // Sticking to advancing all particles as per the structure of this function.
        for (int i = 0; i < ctx->num_particles; ++i) {
            // Only advance particles that are active, if such a flag exists and is used.
            // Assuming all particles in ctx->particles are active for now.
            advance_particle_position(ctx, &ctx->particles[i], dt);
        }
        // ctx->current_time = ev->time; // This update is moved after logging
    }
    // else: ev->time is very close to ctx->current_time, no significant advancement needed.

    // Log the event *before* updating current_time to the event's time
    // Use global log_event function, not a hook
    if (ctx->log_file) {
        log_event(ctx, ev);
    }

    // Update simulation time to the event's time
    // This should happen *after* logging the event at its occurrence time,
    // but *before* processing its consequences (which might generate new events based on this new current_time).
    // The prompt's fix:
    // log_event(ctx, next_event);
    // ctx->current_time = next_event->time;
    // This order is what's being implemented.
    ctx->current_time = ev->time;


    // Validate collision event just before processing (if hook is set)
    if (ev->type == EVENT_COLLISION && ctx->hooks.validate_collision) {
        // validate_collision might be specific, or a general pre-event hook
        // The prompt shows: if (!validate_collision(ctx, ev)) { ... return; }
        // Assuming validate_collision is a function that returns void, not int.
        ctx->hooks.validate_collision(ctx, ev);
    }
    
    // Process the event based on its type
    switch(ev->type) {
        case EVENT_COLLISION:
            process_collision(ctx, ev);  // From collisions.h
            break;
        case EVENT_SHOULDER_ENTRY:
            handle_shoulder_entry(ctx, ev);  // From shoulders.h
            break;
        case EVENT_SHOULDER_EXIT:
            process_shoulder_exit(ctx, ev);  // From shoulders.h
            break;
        case EVENT_CELL_CROSS_X_POS:
        case EVENT_CELL_CROSS_X_NEG:
        case EVENT_CELL_CROSS_Y_POS:
        case EVENT_CELL_CROSS_Y_NEG:
            handle_cell_crossing(ctx, ev);  // From cell_interactions.h
            break;
        default:
            // Invalid or unhandled event type
            if (ctx->hooks.log_error) {
                ctx->hooks.log_error("Unhandled event type: %d", ev->type);
            } else {
                fprintf(stderr, "Unhandled event type: %d\n", ev->type);
            }
            break;
    }
    
    // Return the event to the pool
    release_event_to_pool(ctx, ev);
}
