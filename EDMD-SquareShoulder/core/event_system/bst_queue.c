#include "bst_queue.h"
#include "../sim_context.h"
#include <stdlib.h>
#include <stdio.h>

/**
 * @brief Insert an event into the BST
 * 
 * @param root Pointer to the root of the BST
 * @param ev Event to insert
 */
void bst_insert(Event** root, Event* ev) {
    if (!root || !ev) return;
    
    // Clear any existing BST pointers
    ev->left = ev->right = ev->parent = NULL;
    
    // If the tree is empty, make the event the root
    if (!*root) {
        *root = ev;
        return;
    }
    
    // Find the correct position in the tree
    Event* current = *root;
    Event* parent = NULL;
    
    while (current) {
        parent = current;
        
        if (ev->time < current->time) {
            current = current->left;
        } else {
            current = current->right;
        }
    }
    
    // Insert the event as a child of the parent
    ev->parent = parent;
    
    if (ev->time < parent->time) {
        parent->left = ev;
    } else {
        parent->right = ev;
    }
}

/**
 * @brief Get the earliest event in the BST without removing it
 * 
 * @param root Root of the BST
 * @return Event* Earliest event or NULL if tree is empty
 */
Event* bst_get_earliest(Event* root) {
    if (!root) return NULL;
    
    // The leftmost node has the earliest time
    Event* current = root;
    while (current->left) {
        current = current->left;
    }
    
    return current;
}

/**
 * @brief Remove an event from the BST
 * 
 * @param root Pointer to the root of the BST
 * @param ev Event to remove
 */
void bst_remove(Event** root, Event* ev) {
    if (!root || !*root || !ev) return;
    
    // Step 1: Find the node to remove if not directly provided
    Event* node_to_remove = ev;
    
    // Step 2: Handle the case where the node has at most one child
    if (!node_to_remove->left || !node_to_remove->right) {
        Event* replacement = node_to_remove->left ? node_to_remove->left : node_to_remove->right;
        
        // Case 2a: Node has no children
        if (!replacement) {
            // If node is the root, set root to NULL
            if (node_to_remove->parent == NULL) {
                *root = NULL;
            } else {
                // Update parent's child pointer
                if (node_to_remove == node_to_remove->parent->left) {
                    node_to_remove->parent->left = NULL;
                } else {
                    node_to_remove->parent->right = NULL;
                }
            }
        } 
        // Case 2b: Node has one child
        else {
            // Update parent's child pointer
            replacement->parent = node_to_remove->parent;
            
            if (node_to_remove->parent == NULL) {
                // Node is the root
                *root = replacement;
            } else if (node_to_remove == node_to_remove->parent->left) {
                node_to_remove->parent->left = replacement;
            } else {
                node_to_remove->parent->right = replacement;
            }
        }
    }
    // Step 3: Handle the case where the node has two children
    else {
        // Find the inorder successor (the leftmost node in the right subtree)
        Event* successor = node_to_remove->right;
        while (successor->left) {
            successor = successor->left;
        }
        
        // If the successor is not the immediate right child
        if (successor->parent != node_to_remove) {
            // Replace successor with its right child in its current position
            if (successor->right) {
                successor->right->parent = successor->parent;
            }
            successor->parent->left = successor->right;
            
            // Update successor's right child
            successor->right = node_to_remove->right;
            successor->right->parent = successor;
        }
        
        // Replace node_to_remove with successor
        successor->parent = node_to_remove->parent;
        successor->left = node_to_remove->left;
        successor->left->parent = successor;
        
        if (node_to_remove->parent == NULL) {
            // Node is the root
            *root = successor;
        } else if (node_to_remove == node_to_remove->parent->left) {
            node_to_remove->parent->left = successor;
        } else {
            node_to_remove->parent->right = successor;
        }
    }
    
    // Clear pointers in the removed node
    node_to_remove->left = node_to_remove->right = node_to_remove->parent = NULL;
}

/**
 * @brief Wrapper to remove an event from the BST within the simulation context.
 * 
 * @param ctx The simulation context.
 * @param ev_to_remove The event to remove from the BST.
 */
void remove_event_from_bst(SimContext* ctx, Event* ev_to_remove) {
    if (!ctx || !ev_to_remove) return;
    bst_remove(&ctx->event_system.event_tree_root, ev_to_remove);
}

/**
 * @brief Remove all events for a specific particle from the BST
 * 
 * @param root Pointer to the root of the BST
 * @param p_idx Index of the particle
 * @param free_list Pointer to the head of the free event list
 */
void bst_remove_events_for_particle(Event** root, int p_idx, Event** free_list) {
    if (!root || !*root || !free_list || p_idx < 0) return;
    
    // We need to use a recursive approach to safely traverse and modify the tree
    bst_remove_events_for_particle_recursive(root, *root, p_idx, free_list);
}

/**
 * @brief Helper function for recursively removing events for a particle
 * 
 * @param root Pointer to the root of the BST
 * @param node Current node being examined
 * @param p_idx Index of the particle
 * @param free_list Pointer to the head of the free event list
 */
void bst_remove_events_for_particle_recursive(Event** root, Event* node, int p_idx, Event** free_list) {
    if (!node) return;
    
    // First, recursively process children
    // We need to save these pointers because we might remove 'node'
    Event* left = node->left;
    Event* right = node->right;
    
    // Check if this node involves the particle
    if (node->p1_idx == p_idx || node->p2_idx == p_idx) {
        // Remove this node
        Event* to_free = node;
        bst_remove(root, node);
        
        // Add to free list
        to_free->next = *free_list;
        *free_list = to_free;
    }
    
    // Continue with children
    bst_remove_events_for_particle_recursive(root, left, p_idx, free_list);
    bst_remove_events_for_particle_recursive(root, right, p_idx, free_list);
}
