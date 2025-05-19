#ifndef CELL_LIST_H
#define CELL_LIST_H

#include "../core/sim_context.h"

/**
 * @brief Initialize cell lists for neighbor search
 * 
 * @param ctx Simulation context
 * @param n_cells_x Number of cells in x dimension
 * @param n_cells_y Number of cells in y dimension
 */
void init_cell_lists(SimContext* ctx, int n_cells_x, int n_cells_y);

/**
 * @brief Update a particle's cell assignment
 * 
 * @param ctx Simulation context
 * @param p Particle to update
 */
void update_particle_cell(SimContext* ctx, Particle* p);

/**
 * @brief Update a particle's membership in the cell list data structure
 * 
 * @param ctx Simulation context
 * @param p Particle to update
 * @param old_cellx Previous cell x index
 * @param old_celly Previous cell y index
 */
void update_particle_cell_membership(SimContext* ctx, Particle* p, int old_cellx, int old_celly);

/**
 * @brief Get list of particles in neighboring cells
 * 
 * @param ctx Simulation context
 * @param cellx Cell x index
 * @param celly Cell y index
 * @param neighbors Array to store particle pointers
 * @param max_neighbors Maximum number of neighbors to store
 * @return int Number of neighbors found
 */
int get_cell_neighbors(SimContext* ctx, int cellx, int celly, Particle** neighbors, int max_neighbors);

/**
 * @brief Clear cell list data structures and free memory
 * 
 * @param ctx Simulation context
 */
void free_cell_lists(SimContext* ctx);

#endif // CELL_LIST_H
