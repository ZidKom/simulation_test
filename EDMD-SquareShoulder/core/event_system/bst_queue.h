#ifndef BST_QUEUE_H
#define BST_QUEUE_H

#include "../sim_context.h"

/**
 * @brief Insert an event into the BST
 * 
 * @param root Pointer to the root of the BST
 * @param ev Event to insert
 */
void bst_insert(Event** root, Event* ev);

/**
 * @brief Get the earliest event in the BST without removing it
 * 
 * @param root Root of the BST
 * @return Event* Earliest event or NULL if tree is empty
 */
Event* bst_get_earliest(Event* root);

/**
 * @brief Remove an event from the BST
 * 
 * @param root Pointer to the root of the BST
 * @param ev Event to remove
 */
void bst_remove(Event** root, Event* ev);

/**
 * @brief Wrapper to remove an event from the BST within the simulation context.
 * 
 * @param ctx The simulation context.
 * @param ev_to_remove The event to remove from the BST.
 */
void remove_event_from_bst(SimContext* ctx, Event* ev_to_remove);

/**
 * @brief Remove all events for a specific particle from the BST
 * 
 * @param root Pointer to the root of the BST
 * @param p_idx Index of the particle
 * @param free_list Pointer to the head of the free event list
 */
void bst_remove_events_for_particle(Event** root, int p_idx, Event** free_list);

/**
 * @brief Helper function for recursively removing events for a particle
 * 
 * @param root Pointer to the root of the BST
 * @param node Current node being examined
 * @param p_idx Index of the particle
 * @param free_list Pointer to the head of the free event list
 */
void bst_remove_events_for_particle_recursive(Event** root, Event* node, int p_idx, Event** free_list);

#endif // BST_QUEUE_H
