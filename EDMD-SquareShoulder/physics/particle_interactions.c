#include "particle_interactions.h"
#include "event_prediction.h"
#include "../core/event_system/hybrid.h" // For invalidate_events_for_particle
#include "../utils/pbc.h" // For pbc_min_image_delta

/**
 * @brief Tracks shoulder interaction states between particle pairs
 * 
 * This function updates the shoulder interaction state between two particles
 * and tracks which pairs of particles are currently within shoulder distance.
 * 
 * @param ctx Simulation context
 * @param p1_idx Index of first particle
 * @param p2_idx Index of second particle
 * @param in_shoulder Flag indicating if particles are now in shoulder interaction (1) or not (0)
 */
void track_shoulder_state(SimContext* ctx, int p1_idx, int p2_idx, int in_shoulder) {
    if (!ctx || p1_idx < 0 || p1_idx >= ctx->num_particles || 
        p2_idx < 0 || p2_idx >= ctx->num_particles) {
        return;
    }
    
    Particle* p1 = &ctx->particles[p1_idx];
    Particle* p2 = &ctx->particles[p2_idx];
    
    // Update in_shoulder flag for both particles
    p1->in_shoulder = in_shoulder;
    p2->in_shoulder = in_shoulder;
    
    // Here we could extend to track specific particle pairs in shoulder interactions
    // if we wanted to allow multiple simultaneous shoulder interactions per particle
}

/**
 * @brief Invalidates all events for a particle and its interaction partners
 * 
 * This comprehensive invalidation ensures that when a particle's state changes,
 * all events involving any of its interaction partners are also invalidated.
 * This is essential for physical correctness in complex many-body interactions.
 * 
 * @param ctx Simulation context
 * @param p_idx Index of the particle whose events should be invalidated
 */
void invalidate_partner_events(SimContext* ctx, int p_idx) {
    if (!ctx || p_idx < 0 || p_idx >= ctx->num_particles) {
        return;
    }
    
    Particle* p = &ctx->particles[p_idx];
    
    // First invalidate this particle's events
    invalidate_events_for_particle(ctx, p_idx);
    
    // Find all particles in the same and neighboring cells
    int cell_x = p->cellx;
    int cell_y = p->celly;
    
    // Loop over 3x3 neighborhood of cells
    for (int dx = -1; dx <= 1; dx++) {
        for (int dy = -1; dy <= 1; dy++) {
            // Calculate cell indices with proper PBC
            int cx = (cell_x + dx + ctx->n_cells_x) % ctx->n_cells_x;
            int cy = (cell_y + dy + ctx->n_cells_y) % ctx->n_cells_y;
            
            // Get the cell
            Cell* cell = &ctx->cells[cx][cy];
            
            // Loop through particles in this cell
            Particle* other = cell->head;
            while (other) {
                // Skip the original particle
                if (other->id != p_idx) {
                    // Check if there's an interaction (currently simplified to just distance check)
                    double dx = pbc_min_image_delta(p->x, other->x, ctx->xsize);
                    double dy = pbc_min_image_delta(p->y, other->y, ctx->ysize);
                    double r_squared = dx*dx + dy*dy;
                    double shoulder_radius = ctx->params.sigma * ctx->params.lambda_shoulder;
                    
                    // If particles are close enough, invalidate the other particle's events too
                    if (r_squared < shoulder_radius * shoulder_radius * 4.0) { // Conservative radius
                        invalidate_events_for_particle(ctx, other->id);
                        
                        // Predict new events for this partner - wrap in null check for test environment
                        if (other) {
                            predict_all_events_for_particle(ctx, other);
                        }
                    }
                }
                other = other->next_in_cell;
            }
        }
    }
}
