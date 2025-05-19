#ifndef LOGGER_H
#define LOGGER_H

#include "../core/sim_context.h" // For SimContext, Event

// Initializes the logger. May open files or set up logging parameters.
// rank is for potential MPI integration, pass 0 for serial.
void init_logger(SimContext *ctx, const char *log_filename_base, int rank);

// Logs a specific event.
void log_event(SimContext *ctx, Event *event);

// Logs overall system statistics (e.g., energy, momentum, event counts).
void log_system_stats(SimContext *ctx);

// Closes any open log files and cleans up logger resources.
void close_logger(SimContext *ctx);

#endif // LOGGER_H
