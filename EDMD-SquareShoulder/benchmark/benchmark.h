#ifndef BENCHMARK_H
#define BENCHMARK_H

#include "../core/sim_context.h"

/**
 * @brief Structure to store benchmark statistics
 */
typedef struct {
    long event_count;          // Total number of events processed
    long collision_count;      // Number of collision events
    long shoulder_entry_count; // Number of shoulder entry events
    long shoulder_exit_count;  // Number of shoulder exit events
    long cell_crossing_count;  // Number of cell crossing events
    long paul_hits;            // Number of events found in Paul list
    long bst_hits;             // Number of events found in BST
    double simulation_time;    // Total simulation time (in simulation units)
    double wall_time;          // Total wall clock time (in seconds)
    double cpu_time;           // Total CPU time (in seconds)
    long max_memory;           // Maximum memory usage (in kilobytes)
} BenchmarkStats;

/**
 * @brief Initializes the benchmark statistics
 * 
 * @param stats Pointer to the benchmark statistics structure
 */
void init_benchmark_stats(BenchmarkStats* stats);

/**
 * @brief Updates statistics when an event is processed
 * 
 * @param ctx Pointer to the simulation context
 * @param stats Pointer to the benchmark statistics structure
 * @param event Pointer to the event being processed
 * @param from_paul_list Flag indicating if the event came from the Paul list
 */
void update_event_stats(SimContext* ctx, BenchmarkStats* stats, Event* event, int from_paul_list);

/**
 * @brief Runs a benchmark simulation and collects performance metrics
 * 
 * @param ctx Pointer to the simulation context
 * @return BenchmarkStats structure with the benchmark results
 */
BenchmarkStats run_benchmark(SimContext* ctx);

/**
 * @brief Prints the benchmark statistics
 * 
 * @param stats Pointer to the benchmark statistics structure
 */
void print_benchmark_stats(BenchmarkStats* stats);

#endif // BENCHMARK_H
