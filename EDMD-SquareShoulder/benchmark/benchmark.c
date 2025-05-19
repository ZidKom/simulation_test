#include "benchmark.h"
#include "../core/event_system/hybrid.h" // For get_next_master_event
#include "../physics/event_processing.h" // For process_event
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#ifdef __linux__
#include <sys/resource.h>
#endif

void init_benchmark_stats(BenchmarkStats* stats) {
    if (!stats) return;
    
    stats->event_count = 0;
    stats->collision_count = 0;
    stats->shoulder_entry_count = 0;
    stats->shoulder_exit_count = 0;
    stats->cell_crossing_count = 0;
    stats->paul_hits = 0;
    stats->bst_hits = 0;
    stats->simulation_time = 0.0;
    stats->wall_time = 0.0;
    stats->cpu_time = 0.0;
    stats->max_memory = 0;
}

void update_event_stats(SimContext* ctx, BenchmarkStats* stats, Event* event, int from_paul_list) {
    if (!stats || !event) return;
    
    stats->event_count++;
    
    // Update event type counters
    switch (event->type) {
        case EVENT_COLLISION:
            stats->collision_count++;
            break;
        case EVENT_SHOULDER_ENTRY:
            stats->shoulder_entry_count++;
            break;
        case EVENT_SHOULDER_EXIT:
            stats->shoulder_exit_count++;
            break;
        case EVENT_CELL_CROSS_X_POS: // Changed from EVENT_CELL_CROSSING
        case EVENT_CELL_CROSS_X_NEG:
        case EVENT_CELL_CROSS_Y_POS:
        case EVENT_CELL_CROSS_Y_NEG:
            stats->cell_crossing_count++;
            break;
    }
    
    // Update event source counters
    if (from_paul_list) {
        stats->paul_hits++;
    } else {
        stats->bst_hits++;
    }
    
    // Update simulation time
    stats->simulation_time = ctx->current_time;
}

BenchmarkStats run_benchmark(SimContext* ctx) {
    BenchmarkStats stats;
    init_benchmark_stats(&stats);
    
    if (!ctx) return stats;
    
    #ifdef __linux__
    struct rusage start_usage, end_usage;
    getrusage(RUSAGE_SELF, &start_usage);
    #endif
    
    clock_t start = clock();
    double start_time = ctx->current_time;
    
    // Process events until end time is reached
    while (ctx->current_time < ctx->sim_end_time) {
        Event* next_event = get_next_master_event(ctx);
        if (!next_event) {
            // No more events in the queue
            if (ctx->hooks.log_warning) {
                ctx->hooks.log_warning("No more events in queue at time %.6f", ctx->current_time);
            }
            break;
        }
        
        // Update event statistics
        int from_paul_list = !next_event->is_in_bst;
        update_event_stats(ctx, &stats, next_event, from_paul_list);
        
        // Process the event
        process_event(ctx, next_event);
    }
    
    // Calculate benchmark metrics
    clock_t end = clock();
    stats.cpu_time = ((double)(end - start)) / CLOCKS_PER_SEC;
    stats.wall_time = stats.cpu_time; // For single-threaded apps, CPU time ~= wall time
    
    #ifdef __linux__
    // Get peak memory usage
    getrusage(RUSAGE_SELF, &end_usage);
    stats.max_memory = end_usage.ru_maxrss;
    #endif
    
    return stats;
}

void print_benchmark_stats(BenchmarkStats* stats) {
    if (!stats) return;
    
    printf("\n===== Simulation Performance Metrics =====\n");
    printf("Events processed: %ld\n", stats->event_count);
    if (stats->cpu_time > 0) {
        printf("Events/sec: %.2f\n", stats->event_count / stats->cpu_time);
    }
    printf("Simulation time: %.2f\n", stats->simulation_time);
    printf("Wall clock time: %.2f seconds\n", stats->wall_time);
    printf("CPU time: %.2f seconds\n", stats->cpu_time);
    
    printf("\nEvent Breakdown:\n");
    printf("- Collisions: %ld (%.1f%%)\n", 
           stats->collision_count, 
           100.0 * stats->collision_count / (stats->event_count ? stats->event_count : 1));
    printf("- Shoulder entries: %ld (%.1f%%)\n", 
           stats->shoulder_entry_count,
           100.0 * stats->shoulder_entry_count / (stats->event_count ? stats->event_count : 1));
    printf("- Shoulder exits: %ld (%.1f%%)\n", 
           stats->shoulder_exit_count,
           100.0 * stats->shoulder_exit_count / (stats->event_count ? stats->event_count : 1));
    printf("- Cell crossings: %ld (%.1f%%)\n", 
           stats->cell_crossing_count,
           100.0 * stats->cell_crossing_count / (stats->event_count ? stats->event_count : 1));
    
    printf("\nEvent System Performance:\n");
    long total_lookups = stats->paul_hits + stats->bst_hits;
    printf("- Paul list hits: %ld (%.1f%%)\n", 
           stats->paul_hits, 
           100.0 * stats->paul_hits / (total_lookups ? total_lookups : 1));
    printf("- BST hits: %ld (%.1f%%)\n", 
           stats->bst_hits,
           100.0 * stats->bst_hits / (total_lookups ? total_lookups : 1));
    
    #ifdef __linux__
    printf("\nMemory Usage:\n");
    printf("- Peak memory: %.2f MB\n", stats->max_memory / 1024.0);
    #endif
    
    printf("=========================================\n");
}
