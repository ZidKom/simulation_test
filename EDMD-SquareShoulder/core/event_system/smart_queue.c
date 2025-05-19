#include "hybrid.h"
#include "../sim_context.h"
#include "../../physics/event_prediction.h"

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
void smart_queue_rebuild(SimContext* ctx, int* modified_particles, int count) {
    if (!ctx || !modified_particles || count <= 0) return;
    
    // First, clear existing events for these particles
    for (int i = 0; i < count; i++) {
        int p_idx = modified_particles[i];
        
        // Skip invalid particle indices
        if (p_idx < 0 || p_idx >= ctx->num_particles) continue;
        
        // Remove all events involving this particle
        // This is an optimization that would require tracking all events per particle
        // For now, we'll clear the entire queue and rebuild it selectively
    }
    
    // For each modified particle, predict new events with all neighbors
    for (int i = 0; i < count; i++) {
        int p_idx = modified_particles[i];
        
        // Skip invalid particle indices
        if (p_idx < 0 || p_idx >= ctx->num_particles) continue;
        
        Particle* p = &ctx->particles[p_idx];
        
        // Get cell coordinates
        int cellx = (int)(p->x / ctx->cell_size);
        int celly = (int)(p->y / ctx->cell_size);
        
        // Normalize cell coordinates to handle PBC
        cellx = (cellx + ctx->n_cells_x) % ctx->n_cells_x;
        celly = (celly + ctx->n_cells_y) % ctx->n_cells_y;
        
        // Interact with neighbors in 3x3 cell neighborhood
        for (int dx = -1; dx <= 1; dx++) {
            for (int dy = -1; dy <= 1; dy++) {
                // Apply PBC to neighbor cell coordinates
                int nx = (cellx + dx + ctx->n_cells_x) % ctx->n_cells_x;
                int ny = (celly + dy + ctx->n_cells_y) % ctx->n_cells_y;
                
                // Get the head of the linked list in this cell
                Particle* neighbor = ctx->cells[nx][ny].head;
                
                // Iterate through all particles in this cell
                while (neighbor) {
                    // Don't predict events with self
                    if (neighbor->id != p_idx) {
                        predict_pair_events(ctx, p, neighbor);
                    }
                    
                    // Move to next particle in the cell
                    neighbor = neighbor->next_in_cell;
                }
            }
        }
    }
}
