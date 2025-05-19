#pragma once
#include <stddef.h> // For size_t
#include <stdio.h>  // For FILE

// Particle structure definition
typedef struct Particle {
    double x, y;           // Position
    double vx, vy;         // Velocity
    double radius;         // Core radius
    double radius_sq;      // Precomputed radius squared for collisions
    double m;              // Mass
    int id;                // Particle ID
    int type;              // Particle type (for LAMMPS compatibility)
    int cellx, celly;      // Current cell indices
    unsigned long coll;    // Collision counter
    int active;            // 1 if active, 0 if in pool
    int in_shoulder;       // Flag if particle is in a shoulder potential with another
    int shoulder_partner_idx; // Index of the particle it's in a shoulder interaction with, -1 if none
    struct Particle* next_in_cell; // For cell list linked list
} Particle;

// Event types enumeration
typedef enum {
    EVENT_NONE, // Should be 0 for default/uninitialized
    EVENT_COLLISION,
    EVENT_SHOULDER_ENTRY, // Renamed from EVENT_SHOULDER for clarity
    EVENT_SHOULDER_EXIT,  // Added for explicit shoulder exit
    EVENT_CELL_CROSS_X_POS, // Cell crossing in positive X direction
    EVENT_CELL_CROSS_X_NEG, // Cell crossing in negative X direction
    EVENT_CELL_CROSS_Y_POS, // Cell crossing in positive Y direction
    EVENT_CELL_CROSS_Y_NEG, // Cell crossing in negative Y direction
    EVENT_INVALID // To mark events that should be ignored/removed
} EventType;

// Event structure definition
typedef struct Event {
    double time;           // Absolute simulation time of the event
    EventType type;        // Type of the event
    int p1_idx;            // Index of first particle
    int p2_idx;            // Index of second particle (if applicable)
    int cross_dir;         // 0: X+, 1: X-, 2: Y+, 3: Y- (for cell crossing events)

    // BST pointers
    struct Event *left, *right, *parent;
    
    // Paul list / free list pointers
    struct Event *next;               // Used for free list linkage
    struct Event *next_event_in_paul_bin; // Pointer to next event in same Paul bin
    struct Event *prev_in_bin;        // Pointer to previous event in Paul bin
    int assigned_paul_bin;            // Index of Paul bin this event belongs to (-1 if not in Paul)
    int is_in_bst;                   // Flag indicating if the event is in the BST
} Event;

// Cell structure for cell lists
typedef struct Cell {
    struct Particle *head;   // Linked list of particles in cell
} Cell;

// Main simulation context structure
typedef struct SimContext {
    // Particle system
    Particle *particles;
    int num_particles;
    double xsize, ysize;   // Simulation box dimensions

    // Event system
    struct {
        Event** paul_lists;          // Array of Paul list heads
        int paul_list_size;          // Number of Paul list bins
        int current_paul_bin;        // Current Paul bin index
        double paul_dt;              // Time width per Paul bin
        double paul_list_window_start_time; // Start time for current Paul window
        Event* event_tree_root;      // BST root
    } event_system;
    
    // Memory pools
    Event *event_pool_storage;       // Preallocated event pool
    int event_pool_size;
    Event *free_event_list;          // Linked list of free events
    
    // Physics parameters
    double current_time;             // Current simulation time
    struct { 
        double U;                    // Square shoulder height
        double sigma;                // Core diameter
        double lambda_shoulder;      // Shoulder width factor
        double deltaE;               // Energy injection per collision
        double default_particle_radius; // Default radius for new particles
        double default_particle_mass;   // Default mass for new particles
    } params;

    // Cell lists
    Cell **cells;                    // 2D array of cells
    int n_cells_x, n_cells_y;        // Number of cells in each dimension
    double cell_size;                // Width/height of cells (assuming square)
    
    // Energy tracking
    struct {
        double tracked;              // Current tracked energy
        double tolerance;            // Tolerance for energy conservation checks
    } energy;
    
    // Validation and I/O
    FILE *log_file;
    int log_verbosity;
    
    // Validation hooks
    struct {
        void (*validate_collision)(struct SimContext*, Event*);
        void (*log_error)(const char* fmt, ...);
        void (*log_warning)(const char* fmt, ...);
        void (*log_info)(const char* fmt, ...);
        void (*on_cell_crossing)(struct SimContext*, Event*);
    } hooks;
    // Simulation control
    double sim_end_time;             // When to end simulation
} SimContext;

// Validation hooks structure
typedef struct ValidationHooks {
    void (*validate_collision)(struct SimContext*, Event*);
    void (*log_error)(const char* fmt, ...);
    void (*log_warning)(const char* fmt, ...);
} ValidationHooks;

// Context initialization/cleanup
void init_sim_context(SimContext *ctx, int num_particles_val, double x_size_val, double y_size_val, 
                     double cell_size_val, int paul_list_size_val, int event_pool_sz_val);
void free_sim_context(SimContext *ctx);

// Helper for aligned allocation
void* aligned_alloc(size_t alignment, size_t size);
void aligned_free(void *ptr);
