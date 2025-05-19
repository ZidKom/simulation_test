#include "../core/sim_context.h"
#include "../core/event_system/events.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <time.h>
#include <omp.h>

// Validation output stub

static const char* event_type_str(EventType type) {
    switch(type) {
        case EVENT_NONE: return "NONE";
        case EVENT_COLLISION: return "COLLISION";
        case EVENT_SHOULDER_ENTRY: return "SHOULDER_ENTRY";
        case EVENT_SHOULDER_EXIT: return "SHOULDER_EXIT";
        case EVENT_CELL_CROSS_X_POS: return "CELL_CROSS_X_POS";
        case EVENT_CELL_CROSS_X_NEG: return "CELL_CROSS_X_NEG";
        case EVENT_CELL_CROSS_Y_POS: return "CELL_CROSS_Y_POS";
        case EVENT_CELL_CROSS_Y_NEG: return "CELL_CROSS_Y_NEG";
        case EVENT_INVALID: return "INVALID";
        default: return "UNKNOWN";
    }
}

void init_logger(SimContext* ctx, const char* log_filename_base, int rank) {
    if (!ctx) return;
    
    char log_path[256];
    if (rank > 0) {
        // For potential MPI support
        snprintf(log_path, sizeof(log_path), "%s.%d.log", log_filename_base, rank);
    } else {
        snprintf(log_path, sizeof(log_path), "%s.log", log_filename_base);
    }
    
    ctx->log_file = fopen(log_path, "a");
    if (!ctx->log_file) {
        fprintf(stderr, "[LOGGER] Failed to open log file: %s\n", log_path);
        exit(EXIT_FAILURE);
    }
    ctx->log_verbosity = 1;
}

void log_event(SimContext* ctx, Event* ev) {
    if (!ctx || !ctx->log_file) return;
    flockfile(ctx->log_file);
    
    if (!ev) {
        fprintf(ctx->log_file, "[TIME %.4f] EVENT NULL\n", ctx->current_time);
        funlockfile(ctx->log_file);
        return;
    }
    
    // Get the particles involved
    Particle *p1 = (ev->p1_idx >= 0) ? &(ctx->particles[ev->p1_idx]) : NULL;
    Particle *p2 = (ev->p2_idx >= 0) ? &(ctx->particles[ev->p2_idx]) : NULL;
    
    if (!p1) {
        fprintf(ctx->log_file, "[TIME %.4f] EVENT %s: (invalid particle indices)\n", 
                ctx->current_time, event_type_str(ev->type));
        funlockfile(ctx->log_file);
        return;
    }
    
    if (p2) {
        // Event involving two particles
        fprintf(ctx->log_file,
            "[TIME %.4f] EVENT %s: p%d-p%d @ (%.3f,%.3f)-(%.3f,%.3f)\n",
            ctx->current_time,
            event_type_str(ev->type),
            p1->id, p2->id,
            p1->x, p1->y, p2->x, p2->y);
    } else {
        // Event involving single particle (e.g., cell crossing)
        const char *direction = "";
        switch(ev->type) {
            case EVENT_CELL_CROSS_X_POS: direction = "X+"; break;
            case EVENT_CELL_CROSS_X_NEG: direction = "X-"; break;
            case EVENT_CELL_CROSS_Y_POS: direction = "Y+"; break;
            case EVENT_CELL_CROSS_Y_NEG: direction = "Y-"; break;
            default: direction = "?"; break;
        }
        
        fprintf(ctx->log_file,
            "[TIME %.4f] EVENT %s: p%d @ (%.3f,%.3f) dir=%s\n",
            ctx->current_time,
            event_type_str(ev->type),
            p1->id,
            p1->x, p1->y,
            direction);
    }
    
    funlockfile(ctx->log_file);
}

void log_system_stats(SimContext* ctx) {
    if (!ctx || !ctx->log_file) return;
    int n_shoulder = 0, n_events = 0;
    
    // Count particles in shoulders
    for (int i = 0; i < ctx->num_particles; ++i) {
        if (ctx->particles[i].in_shoulder) ++n_shoulder;
    }
    
    // Count events in BST
    Event* ev = ctx->event_system.event_tree_root;
    while (ev) { 
        ++n_events; 
        ev = ev->right; 
    }
    
    // Count events in Paul lists
    for (int i = 0; i < ctx->event_system.paul_list_size; ++i) {
        Event* e = ctx->event_system.paul_lists[i];
        while (e) { 
            ++n_events; 
            e = e->next_event_in_paul_bin; 
        }
    }
    
    flockfile(ctx->log_file);
    fprintf(ctx->log_file, "[STATS] t=%.4f N=%d shoulder=%d events=%d\n",
        ctx->current_time, ctx->num_particles, n_shoulder, n_events);
    
    // Optional: Add energy tracking information if enabled
    if (ctx->log_verbosity > 1) {
        fprintf(ctx->log_file, "        Energy: tracked=%.6f tolerance=%.6f\n",
            ctx->energy.tracked, ctx->energy.tolerance);
    }
    
    funlockfile(ctx->log_file);
}

void close_logger(SimContext* ctx) {
    if (!ctx || !ctx->log_file) return;
    
    // Log final stats before closing
    log_system_stats(ctx);
    
    // Add a closing message
    flockfile(ctx->log_file);
    fprintf(ctx->log_file, "[LOGGER] Simulation ended at t=%.4f\n", ctx->current_time);
    funlockfile(ctx->log_file);
    
    // Close the log file
    fclose(ctx->log_file);
    ctx->log_file = NULL;
}
