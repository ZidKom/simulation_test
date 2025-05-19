#include "simulation_hooks.h"
#include <stdio.h>
#include <stdarg.h>

// Error logging function
static void default_log_error(const char* fmt, ...) {
    va_list args;
    va_start(args, fmt);
    fprintf(stderr, "[ERROR] ");
    vfprintf(stderr, fmt, args);
    fprintf(stderr, "\n");
    va_end(args);
}

// Warning logging function
static void default_log_warning(const char* fmt, ...) {
    va_list args;
    va_start(args, fmt);
    fprintf(stderr, "[WARNING] ");
    vfprintf(stderr, fmt, args);
    fprintf(stderr, "\n");
    va_end(args);
}

// Remove the definition of create_default_validation_hooks from this file.
// See validation/hooks.c for implementation.
