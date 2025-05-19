#ifndef SHOULDERS_H
#define SHOULDERS_H

#include "../core/sim_context.h" // For SimContext and Event structures

/**
 * @brief Handles the consequences of two particles entering a square shoulder potential.
 *
 * This function is called when an EVENT_SHOULDER_ENTRY occurs.
 * It predicts the corresponding EVENT_SHOULDER_EXIT and re-predicts other events
 * for the involved particles.
 *
 * @param ctx Pointer to the simulation context.
 * @param event Pointer to the shoulder entry event.
 */
void handle_shoulder_entry(SimContext *ctx, Event *event);

/**
 * @brief Handles the consequences of two particles exiting a square shoulder potential.
 *
 * This function is called when an EVENT_SHOULDER_EXIT occurs.
 * It re-predicts events for the involved particles.
 *
 * @param ctx Pointer to the simulation context.
 * @param event Pointer to the shoulder exit event.
 */
void process_shoulder_exit(SimContext *ctx, Event *event);

#endif // SHOULDERS_H
