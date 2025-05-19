#include "../core/sim_context.h"
#include "cell_list.h"
#include <stdlib.h>
#include <string.h>
#include <assert.h>

/**
 * @brief Initialize cell lists for neighbor search
 * 
 * @param ctx Simulation context
 * @param n_cells_x Number of cells in x dimension
 * @param n_cells_y Number of cells in y dimension
 */
void init_cell_lists(SimContext* ctx, int n_cells_x, int n_cells_y) {
    assert(ctx && n_cells_x > 0 && n_cells_y > 0);
    
    ctx->n_cells_x = n_cells_x;
    ctx->n_cells_y = n_cells_y;
    ctx->cell_size = ctx->xsize / n_cells_x; // Assuming square cells
    
    // Allocate 2D array of cells
    ctx->cells = (Cell**)aligned_alloc(64, n_cells_x * sizeof(Cell*));
    for (int i = 0; i < n_cells_x; ++i) {
        ctx->cells[i] = (Cell*)aligned_alloc(64, n_cells_y * sizeof(Cell));
        for (int j = 0; j < n_cells_y; ++j) {
            ctx->cells[i][j].head = NULL;
        }
    }
}

/**
 * @brief Update a particle's cell assignment
 * 
 * @param ctx Simulation context
 * @param p Particle to update
 */
void update_particle_cell(SimContext* ctx, Particle* p) {
    assert(ctx && p);
    
    // Calculate cell indices based on position
    int cx = (int)(p->x / ctx->cell_size);
    int cy = (int)(p->y / ctx->cell_size);
    
    // Apply periodic boundary conditions
    cx = (cx + ctx->n_cells_x) % ctx->n_cells_x;
    cy = (cy + ctx->n_cells_y) % ctx->n_cells_y;
    
    // Update particle's cell indices
    int old_cx = p->cellx;
    int old_cy = p->celly;
    
    // Only update if cell changed
    if (cx != old_cx || cy != old_cy) {
        update_particle_cell_membership(ctx, p, old_cx, old_cy);
    }
}

/**
 * @brief Update a particle's membership in the cell list data structure
 * 
 * @param ctx Simulation context
 * @param p Particle to update
 * @param old_cellx Previous cell x index
 * @param old_celly Previous cell y index
 */
void update_particle_cell_membership(SimContext* ctx, Particle* p, int old_cellx, int old_celly) {
    // Remove from old cell if valid indices
    if(old_cellx >= 0 && old_celly >= 0) {
        Particle** head = &ctx->cells[old_cellx][old_celly].head;
        while(*head && *head != p) head = &(*head)->next_in_cell;
        if(*head == p) *head = p->next_in_cell;
    }
    
    // Apply PBC to new indices
    p->cellx = (p->cellx + ctx->n_cells_x) % ctx->n_cells_x;
    p->celly = (p->celly + ctx->n_cells_y) % ctx->n_cells_y;
    
    // Add to new cell
    p->next_in_cell = ctx->cells[p->cellx][p->celly].head;
    ctx->cells[p->cellx][p->celly].head = p;
}

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
int get_cell_neighbors(SimContext* ctx, int cellx, int celly, Particle** neighbors, int max_neighbors) {
    int count = 0;
    
    // Loop over cell and its 8 neighbors (assuming 2D)
    for (int i = -1; i <= 1; i++) {
        for (int j = -1; j <= 1; j++) {
            // Apply PBC to neighbor cell indices
            int nx = (cellx + i + ctx->n_cells_x) % ctx->n_cells_x;
            int ny = (celly + j + ctx->n_cells_y) % ctx->n_cells_y;
            
            // Add particles from this cell to neighbors list
            Particle* p = ctx->cells[nx][ny].head;
            while (p && count < max_neighbors) {
                neighbors[count++] = p;
                p = p->next_in_cell;
            }
        }
    }
    
    return count;
}

/**
 * @brief Clear cell list data structures and free memory
 * 
 * @param ctx Simulation context
 */
void free_cell_lists(SimContext* ctx) {
    if (!ctx || !ctx->cells) return;
    
    // Free the memory for each row
    for (int i = 0; i < ctx->n_cells_x; i++) {
        if (ctx->cells[i]) {
            free(ctx->cells[i]);
        }
    }
    
    // Free the row pointers array
    free(ctx->cells);
    ctx->cells = NULL;
    
    // Reset cell count
    ctx->n_cells_x = 0;
    ctx->n_cells_y = 0;
}
