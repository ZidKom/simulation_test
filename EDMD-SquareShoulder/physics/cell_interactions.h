#ifndef CELL_INTERACTIONS_H
#define CELL_INTERACTIONS_H

#include "../core/sim_context.h"

/**
 * @brief Predicts cell crossing events for a particle.
 *
 * Calculates when a particle will cross each boundary of its current cell
 * and schedules corresponding events: EVENT_CELL_CROSS_X_POS, EVENT_CELL_CROSS_X_NEG,
 * EVENT_CELL_CROSS_Y_POS, or EVENT_CELL_CROSS_Y_NEG.
 *
 * @param ctx The simulation context.
 * @param p The particle for which to predict the cell crossings.
 */
void predict_cell_crossing(SimContext *ctx, Particle *p);

/**
 * @brief Generic handler for cell crossing events
 * 
 * Dispatches to the specific handler based on event type
 * 
 * @param ctx Simulation context
 * @param ev Cell crossing event (must be one of the cell crossing event types)
 */
void handle_cell_crossing(SimContext *ctx, Event *ev);

/**
 * @brief Handles a cell crossing event in the positive X direction
 * 
 * @param ctx Simulation context
 * @param p Particle crossing the cell boundary
 */
void handle_cell_cross_x_pos(SimContext *ctx, Particle *p);

/**
 * @brief Handles a cell crossing event in the negative X direction
 * 
 * @param ctx Simulation context
 * @param p Particle crossing the cell boundary
 */
void handle_cell_cross_x_neg(SimContext *ctx, Particle *p);

/**
 * @brief Handles a cell crossing event in the positive Y direction
 * 
 * @param ctx Simulation context
 * @param p Particle crossing the cell boundary
 */
void handle_cell_cross_y_pos(SimContext *ctx, Particle *p);

/**
 * @brief Handles a cell crossing event in the negative Y direction
 * 
 * @param ctx Simulation context
 * @param p Particle crossing the cell boundary
 */
void handle_cell_cross_y_neg(SimContext *ctx, Particle *p);

#endif // CELL_INTERACTIONS_H
