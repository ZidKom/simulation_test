#ifndef PARTICLE_POOL_H
#define PARTICLE_POOL_H

#include "sim_context.h" // For Event and SimContext (specifically event_pool and free_event_list)

/**
 * @brief Initializes the event pool and free list within the SimContext.
 *
 * This function is typically called by init_sim_context.
 * It pre-allocates a block of Event structures and links them into a free list.
 *
 * @param ctx Pointer to the simulation context.
 * @param pool_size The number of Event objects to pre-allocate in the pool.
 * @return 0 on success, -1 on failure (e.g., memory allocation failed).
 */
int init_event_pool(SimContext *ctx, int pool_size);

/**
 * @brief Allocates an Event object from the SimContext's event pool.
 *
 * Retrieves an event from the free list. If the free list is empty,
 * it may attempt to expand the pool (not implemented here) or return NULL.
 *
 * @param ctx Pointer to the simulation context.
 * @return Pointer to an allocated Event, or NULL if the pool is exhausted.
 */
Event* allocate_event_from_pool(SimContext *ctx);

/**
 * @brief Releases an Event object back to the SimContext's event pool.
 *
 * Adds the event to the head of the free list, making it available for reuse.
 *
 * @param ctx Pointer to the simulation context.
 * @param event Pointer to the Event object to release.
 */
void release_event_to_pool(SimContext *ctx, Event *event);

/**
 * @brief Frees the memory allocated for the event pool in SimContext.
 *
 * This function is typically called by free_sim_context.
 *
 * @param ctx Pointer to the simulation context.
 */
void free_event_pool(SimContext *ctx);


#endif // PARTICLE_POOL_H
