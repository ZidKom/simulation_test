#ifndef EVENT_PROCESSING_H
#define EVENT_PROCESSING_H

#include "../core/sim_context.h"

/**
 * @brief Processes an event from the event queue
 * 
 * This function dispatches the event to the appropriate handler based on the event type.
 * 
 * @param ctx Pointer to the simulation context
 * @param ev Pointer to the event to process
 */
void process_event(SimContext* ctx, Event* ev);

#endif // EVENT_PROCESSING_H
