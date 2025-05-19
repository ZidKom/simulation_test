/**
 * @file main.c
 * @brief Main simulation driver for Event-Driven Molecular Dynamics with Square Shoulder potential
 * 
 * This implementation provides:
 * 1. Efficient event-driven molecular dynamics in 2D with square shoulder potential
 * 2. Cell-list based neighbor search for O(N) scaling
 * 3. Hybrid event queue (Paul list + BST) for efficient event management
 * 4. Optimized event prediction and processing
 * 5. Rigorous energy conservation checks
 * 6. Support for periodic boundary conditions
 * 
 * Major optimizations:
 * - Only advances positions of particles involved in events
 * - Uses numerically stable quadratic solver for event prediction
 * - Implements full event invalidation for all partner particles
 * - Robustly tracks shoulder interaction states
 * - Properly manages cell list crossings
 * - Uses impulse-based energy injection in collisions
 * 
 * @author EDMD-SquareShoulder Development Team
 */

#include "core/sim_context.h"
#include "core/event_system/hybrid.h" // For get_next_master_event, schedule_event
#include "core/event_system/events.h" // For EventType
#include "core/event_system/paul_list.h" // For advance_paul_list
#include "core/particle_pool.h" // For release_event_to_pool
#include "physics/physics.h"
#include "physics/collisions.h"
#include "physics/shoulders.h"
#include "physics/event_prediction.h"
#include "physics/particle_motion.h" // For advance_particle_position
#include "physics/particle_interactions.h" // For partner invalidation
#include "physics/cell_interactions.h" // For handle_cell_crossing
#include "physics/event_processing.h" // For process_event
#include "utils/cell_list.h"
#include "utils/pbc.h" // For periodic boundary conditions
#include "io/logger.h"
#include "io/snapshot.h"
#include "io/particle_loader.h" // For load_particles
#include "io/lammps_loader.h" // For LAMMPS input format
#include "validation/hooks.h" // For validate_initial_conditions
#include "validation/energy_audit.h"
#include "validation/collision_validation.h"
#include "validation/trajectory_validation.h"
#include "benchmark/benchmark.h"
#include "simulation_hooks.h"
#include "utils/velocity_init.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <getopt.h>

// Default simulation parameters - can be overridden by config file or command line
#define DEFAULT_NUM_PARTICLES 100
#define DEFAULT_SIM_END_TIME 10.0
#define DEFAULT_BOX_SIZE 10.0
#define DEFAULT_PAUL_LIST_SIZE 100
#define DEFAULT_PAUL_DT 0.01
#define DEFAULT_SHOULDER_U 0.5
#define DEFAULT_SHOULDER_SIG 0.2
#define DEFAULT_PARTICLE_RADIUS 0.5
#define DEFAULT_PARTICLE_MASS 1.0
#define DEFAULT_SNAPSHOT_FREQ 1.0
#define DEFAULT_LOG_FREQ 0.1
#define DEFAULT_VALIDATION_FREQ 1000 // Number of events

void print_usage(char *prog_name) {
    fprintf(stderr, "Usage: %s [options]\n", prog_name);
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "  -h, --help             Show this help message\n");
    fprintf(stderr, "  -n, --num_particles    Number of particles (default: %d)\n", DEFAULT_NUM_PARTICLES);
    fprintf(stderr, "  -t, --sim_end_time     Simulation end time (default: %.2f)\n", DEFAULT_SIM_END_TIME);
    fprintf(stderr, "  -L, --box_size         Box size (default: %.2f)\n", DEFAULT_BOX_SIZE);
    fprintf(stderr, "  -s, --snapshot_file    Base name for snapshot files\n");
    fprintf(stderr, "  -l, --log_file         Base name for log files\n");
    fprintf(stderr, "  -i, --input            Input file for initial particle positions\n");
    fprintf(stderr, "  -b, --benchmark        Run in benchmark mode\n");
    fprintf(stderr, "  -v, --validate         Run in validation mode\n");
    fprintf(stderr, "  -r, --reference        Reference trajectory file for validation\n");
    fprintf(stderr, "  --lammps               Use LAMMPS data format for input file\n");
    fprintf(stderr, "  --seed <val>           Random seed for reproducibility\n");
    fprintf(stderr, "  --tolerance <val>      Tolerance for validation (default: 1e-9)\n");
    fprintf(stderr, "  --paul_size <val>      Paul list size (default: %d)\n", DEFAULT_PAUL_LIST_SIZE);
    fprintf(stderr, "  --paul_dt <val>        Paul list dt (default: %.2f)\n", DEFAULT_PAUL_DT);
    fprintf(stderr, "  --shoulder_U <val>     Shoulder potential height U (default: %.2f)\n", DEFAULT_SHOULDER_U);
    fprintf(stderr, "  --shoulder_sig <val>   Shoulder potential width sigma (default: %.2f)\n", DEFAULT_SHOULDER_SIG);
    fprintf(stderr, "  --particle_radius <val> Particle radius (default: %.2f)\n", DEFAULT_PARTICLE_RADIUS);
    fprintf(stderr, "  --particle_mass <val>   Particle mass (default: %.2f)\n", DEFAULT_PARTICLE_MASS);
    fprintf(stderr, "  --snapshot_freq <val>  Snapshot frequency (sim time) (default: %.2f)\n", DEFAULT_SNAPSHOT_FREQ);
    fprintf(stderr, "  --log_freq <val>       Log frequency (sim time) (default: %.2f)\n", DEFAULT_LOG_FREQ);
    fprintf(stderr, "  --validation_freq <val> Validation frequency (events) (default: %d)\n", DEFAULT_VALIDATION_FREQ);
    // TODO: Add options for config file, initial particle configuration file
}

// Placeholder for initializing particles with random positions and velocities
void random_init_particles(SimContext *ctx) {
    srand(time(NULL)); // Seed random number generator
    for (int i = 0; i < ctx->num_particles; ++i) {
        ctx->particles[i].id = i;
        ctx->particles[i].x = ((double)rand() / RAND_MAX) * ctx->xsize;
        ctx->particles[i].y = ((double)rand() / RAND_MAX) * ctx->ysize;
        // No z coordinate in 2D simulation
        ctx->particles[i].vx = ((double)rand() / RAND_MAX - 0.5) * 2.0; // Random velocity between -1 and 1
        ctx->particles[i].vy = ((double)rand() / RAND_MAX - 0.5) * 2.0;
        // No vz in 2D simulation
        ctx->particles[i].radius = ctx->params.default_particle_radius; 
        // mass is m in our Particle structure
        ctx->particles[i].m = ctx->params.default_particle_mass;
        // cellx and celly instead of cell_idx_x and cell_idx_y
        ctx->particles[i].cellx = -1; // Will be updated by init_cell_lists or update_particle_cell
        ctx->particles[i].celly = -1;
        ctx->particles[i].active = 1; // Mark as active
        // Initialize other particle fields if necessary
    }
}

int main(int argc, char **argv) {
    SimContext ctx;
    // Initialize with default parameters first
    ctx.num_particles = DEFAULT_NUM_PARTICLES;
    ctx.sim_end_time = DEFAULT_SIM_END_TIME;
    ctx.xsize = DEFAULT_BOX_SIZE;
    ctx.ysize = DEFAULT_BOX_SIZE;
    // Set event system parameters
    ctx.event_system.paul_list_size = DEFAULT_PAUL_LIST_SIZE;
    ctx.event_system.paul_dt = DEFAULT_PAUL_DT;
    ctx.params.U = DEFAULT_SHOULDER_U;
    ctx.params.sigma = DEFAULT_SHOULDER_SIG; // sigma, not sig
    ctx.params.default_particle_radius = DEFAULT_PARTICLE_RADIUS;
    ctx.params.default_particle_mass = DEFAULT_PARTICLE_MASS;
    double snapshot_freq = DEFAULT_SNAPSHOT_FREQ;
    double log_freq = DEFAULT_LOG_FREQ;
    long validation_freq = DEFAULT_VALIDATION_FREQ;
    char *snapshot_file_base = NULL;
    char *log_file_base = "simulation_log"; // Default log file base name

    // Command line argument parsing
    int opt;
    int benchmark_mode = 0;
    int validate_mode = 0;
    int use_lammps_input = 0;
    char *input_file = NULL;
    char *ref_file = NULL;
    int random_seed = (int)time(NULL);
    double tolerance = 1e-9;
    
    const struct option long_options[] = {
        {"help", no_argument, NULL, 'h'},
        {"num_particles", required_argument, NULL, 'n'},
        {"sim_end_time", required_argument, NULL, 't'},
        {"box_size", required_argument, NULL, 'L'},
        {"snapshot_file", required_argument, NULL, 's'},
        {"log_file", required_argument, NULL, 'l'},
        {"paul_size", required_argument, NULL, 1},
        {"paul_dt", required_argument, NULL, 2},
        {"shoulder_U", required_argument, NULL, 3},
        {"shoulder_sig", required_argument, NULL, 4},
        {"particle_radius", required_argument, NULL, 5},
        {"particle_mass", required_argument, NULL, 6},
        {"snapshot_freq", required_argument, NULL, 7},
        {"log_freq", required_argument, NULL, 8},
        {"validation_freq", required_argument, NULL, 9},
        {"input", required_argument, NULL, 'i'},
        {"lammps", no_argument, NULL, 10},
        {"benchmark", no_argument, NULL, 'b'},
        {"validate", no_argument, NULL, 'v'},
        {"reference", required_argument, NULL, 'r'},
        {"seed", required_argument, NULL, 11},
        {"tolerance", required_argument, NULL, 12},
        {0, 0, 0, 0}
    };

    while ((opt = getopt_long(argc, argv, "hn:t:L:s:l:i:bvr:", long_options, NULL)) != -1) {
        switch (opt) {
            case 'h': print_usage(argv[0]); return 0;
            case 'n': ctx.num_particles = atoi(optarg); break;
            case 't': ctx.sim_end_time = atof(optarg); break;
            case 'L': 
                ctx.xsize = atof(optarg);
                ctx.ysize = atof(optarg); // Assuming square box from single arg
                break;
            case 's': snapshot_file_base = optarg; break;
            case 'l': log_file_base = optarg; break;
            case 'i': input_file = optarg; break;
            case 'b': benchmark_mode = 1; break;
            case 'v': validate_mode = 1; break;
            case 'r': ref_file = optarg; break;
            case 1: ctx.event_system.paul_list_size = atoi(optarg); break;
            case 2: ctx.event_system.paul_dt = atof(optarg); break;
            case 3: ctx.params.U = atof(optarg); break;
            case 4: ctx.params.sigma = atof(optarg); break;
            case 5: ctx.params.default_particle_radius = atof(optarg); break;
            case 6: ctx.params.default_particle_mass = atof(optarg); break;
            case 7: snapshot_freq = atof(optarg); break;
            case 8: log_freq = atof(optarg); break;
            case 9: validation_freq = atol(optarg); break;
            case 10: use_lammps_input = 1; break;
            case 11: random_seed = atoi(optarg); break;
            case 12: tolerance = atof(optarg); break;
            default: print_usage(argv[0]); return 1;
        }
    }

    // Initialize simulation context (allocates memory for particles, event systems, etc.)
    // Heuristic for event pool size
    // int event_pool_heuristic_size = (ctx.num_particles > 1) ? (ctx.num_particles * (ctx.num_particles - 1)) : 1000;
    // if (event_pool_heuristic_size < 10000 && ctx.num_particles > 50) event_pool_heuristic_size = 10000 + ctx.num_particles * 50;
    
    long long event_pool_sz_val; // Use long long to avoid overflow with large num_particles
    if (ctx.num_particles <= 1) {
        event_pool_sz_val = 50000; // Minimum for very few particles
    } else {
        // Heuristic: num_particles * (num_particles + 5) * 4
        // Ensure this doesn't overflow standard int if intermediate products are large
        long long np = ctx.num_particles;
        event_pool_sz_val = np * (np + 5) * 4;
        if (event_pool_sz_val < (50000 + np * 200)) {
            event_pool_sz_val = 50000 + np * 200;
        }
    }
    // Cap the event pool size to a reasonable maximum if necessary, e.g., 100 million, to prevent excessive memory allocation.
    long long max_event_pool_sz = 100 * 1000 * 1000;
    if (event_pool_sz_val > max_event_pool_sz) {
        event_pool_sz_val = max_event_pool_sz;
    }
    if (event_pool_sz_val < 0) { // Check for overflow if long long was not enough (highly unlikely)
        event_pool_sz_val = max_event_pool_sz;
    }


    // Calculate cell_size_val_calc based on particle radius and shoulder parameters
    // Ensure params.default_particle_radius and params.sigma are set (they are by default or cmd line)
    double cell_size_val_calc = 2.0 * ctx.params.default_particle_radius * (1.0 + ctx.params.sigma);
    if (cell_size_val_calc <= 0) { // Fallback if parameters are zero or negative
        cell_size_val_calc = ctx.xsize / 10.0; // Default to 1/10th of box size
        if (cell_size_val_calc <= 0) cell_size_val_calc = 1.0; // Absolute fallback
        if (ctx.hooks.log_warning) {
            ctx.hooks.log_warning("Cell size calculation resulted in non-positive value. Using fallback: %.2f", cell_size_val_calc);
        }
    }

    init_sim_context(&ctx, ctx.num_particles, ctx.xsize, ctx.ysize, 
                    cell_size_val_calc, // Correctly calculated cell_size
                    ctx.event_system.paul_list_size, 
                    (int)event_pool_sz_val); // event_pool_size based on new heuristic
    
    // Initialize the Paul list timing window
    init_paul_event_system_component(&ctx);
    
    // Use consistent random seed for reproducibility
    srand(random_seed);
    printf("Using random seed: %d\n", random_seed);
                    
    // Initialize particles (randomly or from file)
    if (input_file) {
        int load_status = -1;
        if (use_lammps_input) {
            printf("Loading particles from LAMMPS file: %s\n", input_file);
            load_status = load_lammps_data(&ctx, input_file);
        } else {
            printf("Loading particles from regular file: %s\n", input_file);
            load_status = load_particles(&ctx, input_file);
        }
        if (load_status != 0) {
            fprintf(stderr, "Critical Error: Failed to load particles from %s\n", input_file);
            free_sim_context(&ctx);
            exit(EXIT_FAILURE);
        }
        // Initialize velocities for loaded particles (uniform, zero COM)
        initialize_particle_velocities(&ctx, 1.0); // 1.0 is max_abs_velocity_component, can be parameterized
    } else {
        printf("Initializing %d particles randomly\n", ctx.num_particles);
        random_init_particles(&ctx);
    }

    // After loading particles, validate and initialize cell lists
    validate_initial_conditions(&ctx);
    init_cell_lists(&ctx, ctx.n_cells_x, ctx.n_cells_y);
    
    // Make sure all particles are in the right cells
    for(int i=0; i < ctx.num_particles; i++) {
        update_particle_cell(&ctx, &ctx.particles[i]);
    }
    
    // Predict initial events
    predict_cell_crossing_events(&ctx); 
    for(int i=0; i < ctx.num_particles; i++) {
        predict_all_events_for_particle(&ctx, &ctx.particles[i]);
    }
    // Initialize logger
    init_logger(&ctx, log_file_base, 0); // rank 0 for serial

    // Setup validation hooks
    SimulationHooks hooks = create_default_validation_hooks();
    ctx.hooks.log_error = hooks.log_error;
    ctx.hooks.log_warning = hooks.log_warning;
    ctx.hooks.validate_collision = hooks.validate_collision;
    
    // Log initial state
    printf("Starting simulation for %d particles, box size (%.2f, %.2f), end time %.2f\n",
           ctx.num_particles, ctx.xsize, ctx.ysize, ctx.sim_end_time);
    log_system_stats(&ctx);
    
    // Save initial snapshot if requested
    if (snapshot_file_base) {
        char snapshot_filename[256];
        snprintf(snapshot_filename, sizeof(snapshot_filename), "%s_initial.snap", snapshot_file_base);
        save_snapshot(&ctx, snapshot_filename);
    }

    // If in benchmark mode, run the benchmark and exit
    if (benchmark_mode) {
        printf("Running in benchmark mode...\n");
        BenchmarkStats stats = run_benchmark(&ctx);
        print_benchmark_stats(&stats);
        
        // Cleanup and exit
        close_logger(&ctx);
        free_sim_context(&ctx);
        return 0;
    }

    // Main simulation loop
    long event_count = 0;
    double last_snapshot_time = ctx.current_time;
    double last_log_time = ctx.current_time;
    clock_t sim_start_cpu_time = clock();

    while (ctx.current_time < ctx.sim_end_time) {
        Event *next_event = get_next_master_event(&ctx);
        if (!next_event) {
            printf("No more events to process. Simulation time: %.4f\n", ctx.current_time);
            break;
        }

        // Log the event before processing
        log_event(&ctx, next_event);
        
        // Process the event (this will update simulation time and particle positions)
        process_event(&ctx, next_event);
        
        event_count++;

        // Periodic actions
        if (ctx.current_time >= last_log_time + log_freq) {
            log_system_stats(&ctx);
            last_log_time = ctx.current_time;
        }
        
        if (snapshot_file_base && ctx.current_time >= last_snapshot_time + snapshot_freq) {
            char snapshot_filename[256];
            snprintf(snapshot_filename, sizeof(snapshot_filename), "%s_t%.2f.snap", snapshot_file_base, ctx.current_time);
            save_snapshot(&ctx, snapshot_filename);
            last_snapshot_time = ctx.current_time;
        }
        
        if (event_count % validation_freq == 0) {
            validate_energy(&ctx);
            // Additional validation could be added here
        }

        // Check if Paul list needs advancing
        if (ctx.event_system.paul_list_size > 0 && ctx.event_system.paul_dt > 0) {
            advance_paul_list(&ctx);
        }
    }

    clock_t sim_end_cpu_time = clock();
    double cpu_time_used = ((double) (sim_end_cpu_time - sim_start_cpu_time)) / CLOCKS_PER_SEC;

    printf("Simulation finished at time %.4f after %ld events.\n", ctx.current_time, event_count);
    printf("Total CPU time: %.4f seconds.\n", cpu_time_used);
    printf("Average speed: %.2f events/second\n", event_count / (cpu_time_used > 0 ? cpu_time_used : 1));
    
    log_system_stats(&ctx); // Log final state
    
    // Save final snapshot if requested
    if (snapshot_file_base) {
        char snapshot_filename[256];
        snprintf(snapshot_filename, sizeof(snapshot_filename), "%s_final.snap", snapshot_file_base);
        save_snapshot(&ctx, snapshot_filename);
    }

    // If validation mode and a reference file was provided, compare trajectories
    if (validate_mode && ref_file && snapshot_file_base) {
        char final_snap[256];
        snprintf(final_snap, sizeof(final_snap), "%s_final.snap", snapshot_file_base);
        printf("Validating trajectory against reference file...\n");
        if (compare_trajectories(ref_file, final_snap, tolerance)) {
            printf("Validation successful: Trajectory matches reference within tolerance.\n");
        } else {
            printf("Validation failed: Trajectory differs from reference beyond tolerance.\n");
        }
    }

    // Cleanup
    close_logger(&ctx);
    free_sim_context(&ctx);

    return 0;
}
