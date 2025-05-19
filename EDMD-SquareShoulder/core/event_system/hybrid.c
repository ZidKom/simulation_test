#include "paul_list.h"
#include "bst_queue.h"
#include "../sim_context.h" // Provides Event, SimContext, Particle
#include "../particle_pool.h" // For release_event, allocate_event_from_pool
#include <stdatomic.h>
#include <math.h> // For fabs or other math functions if needed, though not directly used here
#include <stdio.h>  // For NULL
#include <stdlib.h> // For NULL if not in stdio.h for some compilers

// PAUL_WINDOW might be better defined based on parameters if needed, or removed if not used.
// #define PAUL_WINDOW (ctx->params.sigma * 2.0) // Example, ensure ctx is available or pass as param

// Forward declaration for the helper function if it's to be used by others,
// or keep it static if only used within this file.
static void collect_valid_events_from_bst(Event* node, int particle_idx_to_invalidate,
                                          Event** head, Event** tail, SimContext* context);

void schedule_event(SimContext* ctx, Event* ev) {
    if (!ctx || !ev) return;

    // Nullify linking pointers before scheduling
    ev->left = ev->right = ev->parent = NULL;
    ev->next_event_in_paul_bin = NULL;
    // ev->prev_in_bin = NULL; // Assuming prev_in_bin is managed by add_to_paul
    ev->is_in_bst = 0;
    ev->assigned_paul_bin = -1;


    // Try to schedule in Paul list first
    if (ctx->event_system.paul_list_size > 0 && 
        ev->time >= ctx->event_system.paul_list_window_start_time &&
        ev->time < (ctx->event_system.paul_list_window_start_time + ctx->event_system.paul_list_size * ctx->event_system.paul_dt)) {
        
        double dt = ev->time - ctx->event_system.paul_list_window_start_time;
        // Ensure paul_dt is positive to prevent division by zero or negative dt issues
        if (ctx->event_system.paul_dt <= 0) {
            // Fallback to BST if paul_dt is not configured properly
            bst_insert(&(ctx->event_system.event_tree_root), ev);
            ev->is_in_bst = 1;
            return;
        }
        
        int target_bin_relative = (int)(dt / ctx->event_system.paul_dt);
        
        // Correct modulo handling for target_bin_abs
        int target_bin_abs = (ctx->event_system.current_paul_bin + target_bin_relative);
        target_bin_abs = (target_bin_abs % ctx->event_system.paul_list_size + ctx->event_system.paul_list_size) % ctx->event_system.paul_list_size;

        add_to_paul(ctx, ev, target_bin_abs);
        ev->assigned_paul_bin = target_bin_abs;
        // ev->is_in_paul_list = 1; // This flag is not in the provided Event struct, assigned_paul_bin serves similar role
        return;
    }

    // If not in Paul list's range, schedule in BST
    bst_insert(&(ctx->event_system.event_tree_root), ev);
    ev->is_in_bst = 1;
}


Event* get_next_master_event(SimContext* ctx) {
    if (!ctx) return NULL;
    int max_attempts = ctx->event_system.paul_list_size + 2;
    for (int attempt = 0; attempt < max_attempts; ++attempt) {
        Event* paul_candidate = paul_get_earliest(ctx);
        Event* bst_candidate = ctx->event_system.event_tree_root;
        if (bst_candidate) {
            while (bst_candidate->left) bst_candidate = bst_candidate->left;
        }
        Event* chosen_event = NULL;
        if (paul_candidate && bst_candidate) {
            chosen_event = (paul_candidate->time <= bst_candidate->time) ? paul_candidate : bst_candidate;
        } else {
            chosen_event = paul_candidate ? paul_candidate : bst_candidate;
        }
        if (chosen_event) {
            // Remove from correct structure
            if (chosen_event == paul_candidate && !chosen_event->is_in_bst) {
                remove_from_paul(ctx, chosen_event);
            } else if (chosen_event == bst_candidate && chosen_event->is_in_bst) {
                remove_event_from_bst(ctx, chosen_event);
                chosen_event->is_in_bst = 0;
            } else if (chosen_event) {
                // Ambiguous source, fallback
                if (chosen_event->is_in_bst) remove_event_from_bst(ctx, chosen_event);
                else remove_from_paul(ctx, chosen_event);
            }
            // Warn if event time is before current_time
            if (chosen_event->time < ctx->current_time - 1e-9 && ctx->hooks.log_warning) {
                ctx->hooks.log_warning("Chosen event time %.12f is before current sim_time %.12f", chosen_event->time, ctx->current_time);
            }
            return chosen_event;
        }
        // No event found, advance Paul list and retry
        if (ctx->current_time < ctx->sim_end_time) {
            advance_paul_list(ctx);
        } else {
            break;
        }
    }
    if (ctx->hooks.log_error) ctx->hooks.log_error("No valid event found after %d attempts. SimTime: %.12f", max_attempts, ctx->current_time);
    return NULL;
}


void invalidate_events_for_particle(SimContext* ctx, int p_idx) {
    if (!ctx || p_idx < 0 || p_idx >= ctx->num_particles) return;

    // Invalidate events in Paul lists
    if (ctx->event_system.paul_list_size > 0) {
        for (int i = 0; i < ctx->event_system.paul_list_size; ++i) {
            Event* current = ctx->event_system.paul_lists[i];
            while (current != NULL) {
                Event* next = current->next_event_in_paul_bin;
                if (current->p1_idx == p_idx || (current->p2_idx != -1 && current->p2_idx == p_idx)) {
                    if (current->prev_in_bin) current->prev_in_bin->next_event_in_paul_bin = current->next_event_in_paul_bin;
                    else ctx->event_system.paul_lists[i] = current->next_event_in_paul_bin;
                    if (current->next_event_in_paul_bin) current->next_event_in_paul_bin->prev_in_bin = current->prev_in_bin;
                    current->next_event_in_paul_bin = NULL;
                    current->prev_in_bin = NULL;
                    current->assigned_paul_bin = -1;
                    release_event_to_pool(ctx, current);
                }
                current = next;
            }
        }
    }

    // Invalidate events in BST by collecting valid events and rebuilding
    Event* valid_events_head = NULL;
    Event* valid_events_tail = NULL;

    collect_valid_events_from_bst(ctx->event_system.event_tree_root, p_idx, &valid_events_head, &valid_events_tail, ctx);
    
    // Rebuild BST from the collected valid events
    ctx->event_system.event_tree_root = NULL; // Clear the old BST
    Event* current_valid = valid_events_head;
    while (current_valid != NULL) {
        Event* next_valid = current_valid->right; // Using right pointer temporarily from collect_valid_events_from_bst

        // Important: clear temporary list pointers and set BST flags before rescheduling
        current_valid->left = NULL;
        current_valid->right = NULL;
        current_valid->parent = NULL;
        current_valid->assigned_paul_bin = -1;
        current_valid->is_in_bst = 0; // schedule_event will set this to 1 if placed in BST

        schedule_event(ctx, current_valid); // schedule_event will handle BST insertion logic
        current_valid = next_valid;
    }
}

// Helper function to traverse BST and collect non-invalidated events.
// It builds a singly linked list using the 'right' pointer of events.
static void collect_valid_events_from_bst(Event* node, int particle_idx_to_invalidate,
                                          Event** head, Event** tail, SimContext* context) {
    if (node == NULL) {
        return;
    }
    
    // Store children because node might be released or its pointers modified
    Event* left_child = node->left;
    Event* right_child = node->right;

    collect_valid_events_from_bst(left_child, particle_idx_to_invalidate, head, tail, context);
    
    if (node->p1_idx == particle_idx_to_invalidate || (node->p2_idx != -1 && node->p2_idx == particle_idx_to_invalidate)) {
        release_event_to_pool(context, node); // Release the invalidated event
    } else {
        // Add to the list of valid events
        // Clear BST pointers as they are invalid for the new list structure
        node->left = NULL; 
        node->parent = NULL;
        // Use 'right' pointer to form a temporary linked list
        node->right = NULL; 

        if (*head == NULL) {
            *head = node;
            *tail = node;
        } else {
            (*tail)->right = node; 
            *tail = node;
        }
    }
    // Collect from right subtree
    collect_valid_events_from_bst(right_child, particle_idx_to_invalidate, head, tail, context);
}

// Implement schedule_event, paul_get_earliest, remove_from_paul, advance_paul_list, and get_next_master_event
// according to the robust logic described in the prompt.
