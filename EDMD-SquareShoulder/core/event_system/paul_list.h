#ifndef PAUL_LIST_H
#define PAUL_LIST_H

#include "../sim_context.h"
#include "events.h"

/**
 * @brief Initializes the Paul list component of the event system.
 *
 * Assumes that SimContext fields like num_paul_bins, paul_dt, and paul_lists array
 * have been primarily set up by init_sim_context. This function can perform
 * additional checks or minor setup if needed.
 *
 * @param ctx Pointer to the simulation context.
 */
void init_paul_event_system_component(SimContext *ctx);

/**
 * @brief Adds an event to a specific Paul list bin.
 *
 * This function is typically called by schedule_event in hybrid.c when an event
 * is determined to fall within the Paul list time horizon and a target bin
 * has been calculated.
 *
 * @param ctx Pointer to the simulation context.
 * @param ev Pointer to the event to be added.
 * @param target_bin_index The pre-calculated index of the Paul list bin where the event should be placed.
 */
void add_to_paul(SimContext *ctx, Event *ev, int target_bin_index);

/**
 * @brief Removes an event from its Paul list bin.
 *
 * The event's assigned_paul_bin field is used to locate the bin.
 *
 * @param ctx Pointer to the simulation context.
 * @param ev Pointer to the event to be removed.
 */
void remove_from_paul(SimContext *ctx, Event *ev);

/**
 * @brief Advances the Paul list system by one time step (paul_dt).
 *
 * This function is called when the simulation time progresses beyond the current
 * Paul list bin's coverage. It involves:
 * 1. Advancing the paul_list_window_start_time.
 * 2. Moving the current_paul_bin_idx to the next bin.
 * 3. Taking all events from the bin that has just become "stale" (the previous
 *    current_paul_bin_idx) and rescheduling them using schedule_event().
 *    These events might go into new Paul list bins or the BST.
 *
 * @param ctx Pointer to the simulation context.
 */
void advance_paul_list(SimContext *ctx);

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
Event* paul_get_earliest(SimContext *ctx);

#endif // PAUL_LIST_H
