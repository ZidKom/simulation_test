#include "paul_list.h"
#include "../sim_context.h"
#include "hybrid.h"
#include "bst_queue.h"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

/**
 * @brief Initializes the Paul list component of the event system.
 *
 * This function is called by init_sim_context after the main SimContext fields
 * (like num_paul_bins, paul_dt, and the paul_lists array itself) have been allocated and set.
 * It ensures the Paul list bins are empty (NULL pointers) and sets initial time markers.
 *
 * @param ctx Pointer to the simulation context.
 */
void init_paul_event_system_component(SimContext *ctx) {
    assert(ctx != NULL);
    assert(ctx->event_system.paul_lists != NULL);
    assert(ctx->event_system.paul_list_size > 0);
    // The paul_lists array (size paul_list_size + 1 for overflow) should be allocated by init_sim_context.
    // Initialize all bin heads to NULL.
    for (int i = 0; i <= ctx->event_system.paul_list_size; ++i) { // Iterate up to and including the overflow bin
        ctx->event_system.paul_lists[i] = NULL;
    }
    ctx->event_system.current_paul_bin = 0;
    ctx->event_system.paul_list_window_start_time = ctx->current_time; // Initialize Paul window to current sim time
}

/**
 * @brief Adds an event to a specific Paul list bin.
 *
 * The event is added to the head of the singly linked list for that bin.
 * The event's assigned_paul_bin field is updated.
 *
 * @param ctx Pointer to the simulation context.
 * @param ev Pointer to the event to be added.
 * @param target_bin_index The pre-calculated index of the Paul list bin (0 to num_paul_bins, inclusive of overflow).
 */
void add_to_paul(SimContext *ctx, Event *ev, int target_bin_index) {
    assert(ctx != NULL && ev != NULL);
    assert(target_bin_index >= 0 && target_bin_index <= ctx->event_system.paul_list_size);

    // Link event to the head of the target bin's list
    ev->next_event_in_paul_bin = ctx->event_system.paul_lists[target_bin_index];
    ctx->event_system.paul_lists[target_bin_index] = ev;
    ev->assigned_paul_bin = target_bin_index;
}

/**
 * @brief Removes an event from its Paul list bin.
 *
 * The event's assigned_paul_bin field is used to locate the bin.
 * The event must be present in that bin's list.
 *
 * @param ctx Pointer to the simulation context.
 * @param ev Pointer to the event to be removed.
 */
void remove_from_paul(SimContext *ctx, Event *ev) {
    if (!ctx || !ev) return;
    int bin = ev->assigned_paul_bin;
    if (bin < 0 || bin >= ctx->event_system.paul_list_size) return;
    Event **head = &ctx->event_system.paul_lists[bin];
    Event *cur = *head, *prev = NULL;
    while (cur && cur != ev) { prev = cur; cur = cur->next_event_in_paul_bin; }
    if (cur == ev) {
        if (prev) prev->next_event_in_paul_bin = cur->next_event_in_paul_bin;
        else *head = cur->next_event_in_paul_bin;
        cur->next_event_in_paul_bin = NULL;
        cur->assigned_paul_bin = -1;
    }
}

/**
 * @brief Advances the Paul list system by one time step (paul_dt).
 *
 * This function is called when the simulation logic determines it's time to
 * move to the next Paul list bin, typically because the current bin is empty or
 * its time window has been processed.
 *
 * It involves:
 * 1. Identifying the current Paul bin that is being "left".
 * 2. Advancing the `paul_list_window_start_time` by `paul_dt`.
 * 3. Updating `current_paul_bin` to the next bin in the circular buffer.
 * 4. Taking all events from the bin that was just "left".
 * 5. Clearing the "left" bin (setting its head to NULL).
 * 6. Re-scheduling each event from the "left" bin back into the hybrid event system
 *    (using `schedule_event`). These events will be placed according to the new
 *    `paul_list_window_start_time` and `current_paul_bin`.
 *
 * @param ctx Pointer to the simulation context.
 */
void advance_paul_list(SimContext *ctx) {
    if (!ctx || ctx->event_system.paul_list_size == 0) return;

    // The bin that is now "in the past" relative to the window moving forward
    int old_bin_idx = ctx->event_system.current_paul_bin;

    // Advance current bin index and window start time
    ctx->event_system.current_paul_bin = (ctx->event_system.current_paul_bin + 1) % ctx->event_system.paul_list_size;
    ctx->event_system.paul_list_window_start_time += ctx->event_system.paul_dt;

    Event* current_ev = ctx->event_system.paul_lists[old_bin_idx];
    ctx->event_system.paul_lists[old_bin_idx] = NULL; // Clear the old bin

    while (current_ev) {
        Event* next_ev = current_ev->next_event_in_paul_bin;
        
        // Reset Paul list specific flags before rescheduling
        current_ev->next_event_in_paul_bin = NULL;
        // current_ev->prev_in_bin = NULL; // If using doubly linked lists for bins
        current_ev->assigned_paul_bin = -1;
        // current_ev->is_in_paul_list = 0; // If such a flag existed

        // As per prompt: adjust time relative to new window
        current_ev->time -= ctx->event_system.paul_dt * ctx->event_system.paul_list_size;


        schedule_event(ctx, current_ev); // Reschedule the event
        current_ev = next_ev;
    }
}

/**
 * @brief Retrieves the earliest event from the current Paul list bin.
 *
 * This function scans only the current active Paul bin (ctx->event_system.current_paul_bin)
 * to find the event with the smallest time. It does not remove the event.
 * Returns NULL if the current bin is empty.
 *
 * @param ctx Pointer to the simulation context.
 * @return Pointer to the earliest event in the current Paul bin, or NULL if none.
 */
Event* paul_get_earliest(SimContext *ctx) {
    if (!ctx) return NULL;
    int bin = ctx->event_system.current_paul_bin;
    if (bin < 0 || bin >= ctx->event_system.paul_list_size) return NULL;
    Event* e = ctx->event_system.paul_lists[bin];
    Event* min_ev = NULL;
    while (e) {
        if (!min_ev || e->time < min_ev->time) min_ev = e;
        e = e->next_event_in_paul_bin;
    }
    return min_ev;
}
