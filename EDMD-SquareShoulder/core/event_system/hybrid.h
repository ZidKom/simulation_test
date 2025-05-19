#ifndef HYBRID_EVENT_SYSTEM_H
#define HYBRID_EVENT_SYSTEM_H

// Forward declare SimContext and Event to avoid circular dependencies
// Actual definitions will be included in .c file via sim_context.h -> events.h
struct SimContext;
struct Event;

/**
 * @brief Schedules an event into the hybrid event system (Paul list or BST).
 * 
 * The decision of where to place the event (Paul list or BST) is based on the event's time
 * relative to the Paul list's time horizon.
 * Event's list pointers (next, prev) are nullified before scheduling.
 * 
 * @param ctx Pointer to the simulation context.
 * @param event Pointer to the event to be scheduled. The event should be allocated from the event pool.
 *              The event's time and other properties must be set before calling this function.
 */
void schedule_event(struct SimContext* ctx, struct Event* event);

/**
 * @brief Selectively rebuilds the event queue for only the specified particles.
 *
 * Instead of rebuilding the entire event queue, this function only recalculates
 * events for the specified particles and their neighbors, significantly improving
 * performance for localized changes.
 *
 * @param ctx Pointer to the simulation context
 * @param modified_particles Array of particle indices that have been modified
 * @param count Number of particles in the modified_particles array
 */
void smart_queue_rebuild(struct SimContext* ctx, int* modified_particles, int count);

/**
 * @brief Retrieves and removes the next master event (earliest event) from the hybrid system.
 * 
 * This function considers events from matured Paul list bins (handled by advance_paul_list), 
 * the current Paul list bin, and the BST to find the event with the globally minimum time.
 * Events from matured Paul list bins that are not yet due (i.e. their time > current_time despite being in an old bin) 
 * are rescheduled. Due events from matured bins that are not the absolute earliest are also rescheduled.
 * The returned event is removed from the event system.
 * 
 * @param ctx Pointer to the simulation context. The ctx->current_time should be up-to-date.
 * @return Pointer to the next master event, or NULL if no events are pending.
 *         The caller is responsible for processing the event and eventually releasing it
 *         back to the event pool if it's no longer needed. The event's list pointers
 *         (next, prev, next_in_bin, prev_in_bin) will be appropriately managed or nullified.
 */
struct Event* get_next_master_event(struct SimContext* ctx);

/**
 * @brief Invalidates all events related to a specific particle.
 * 
 * This function is used when a particle's state changes significantly (e.g., after a collision,
 * after crossing cell boundaries) and all its future events need to be invalidated.
 * 
 * @param ctx Pointer to the simulation context.
 * @param p_idx Index of the particle whose events should be invalidated.
 */
void invalidate_events_for_particle(struct SimContext* ctx, int p_idx);

#endif // HYBRID_EVENT_SYSTEM_H
