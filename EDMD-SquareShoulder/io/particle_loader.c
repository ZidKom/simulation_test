#include "particle_loader.h"
#include "../core/sim_context.h"
#include "../core/particle_pool.h" 
#include "../utils/pbc.h"
#include "snapshot.h" // For SnapshotFormat, detect_file_format, and SNAPSHOT_FORMAT_* enums
#include <stdio.h>
#include <string.h>
#include <stdlib.h> // For malloc
#include <ctype.h>

// Define MAX_LINE_LEN if not defined elsewhere
#ifndef MAX_LINE_LEN
#define MAX_LINE_LEN 256
#endif

// Define ParticleData if not defined elsewhere
typedef struct {
    int type;
    double x;
    double y;
    double radius;
    double mass;
} ParticleData;

int load_particles(SimContext* ctx, const char* filename) {
    FILE* fp = fopen(filename, "r");
    if (!fp) {
        if (ctx && ctx->hooks.log_error) ctx->hooks.log_error("Could not open file: %s", filename);
        return -1;
    }

    char line[256];
    SnapshotFormat format = detect_file_format(filename);

    switch(format) {
        case SNAPSHOT_FORMAT_XYZ: {
            // Try to detect XYZ format: first line is number of particles, second line is comment
            if (!fgets(line, sizeof(line), fp)) { fclose(fp); return -1; }
            int num_particles = atoi(line);
            if (num_particles <= 0) { fclose(fp); return -1; }
            ctx->num_particles = num_particles;
            if (!fgets(line, sizeof(line), fp)) { fclose(fp); return -1; } // skip comment

            // Allocate particles
            if (ctx->particles) free(ctx->particles);
            ctx->particles = (Particle*)aligned_alloc(64, ctx->num_particles * sizeof(Particle));
            if (!ctx->particles) { fclose(fp); return -1; }
            memset(ctx->particles, 0, ctx->num_particles * sizeof(Particle));

            // Read particle lines: label x y [z]
            for (int i = 0; i < num_particles; ++i) {
                if (!fgets(line, sizeof(line), fp)) { fclose(fp); return -1; }
                char label[16];
                double x, y;
                int n = sscanf(line, "%15s %lf %lf", label, &x, &y);
                if (n < 3) { fclose(fp); return -1; }
                Particle* p = &ctx->particles[i];
                p->id = i;
                p->x = x;
                p->y = y;
                // Apply PBC wrapping to initial coordinates
                if (ctx->xsize > 0 && ctx->ysize > 0) { // Ensure box dimensions are valid
                    pbc_wrap_position(&p->x, &p->y, ctx->xsize, ctx->ysize);
                }
                p->vx = 0.0;
                p->vy = 0.0;
                p->type = 0;
                p->radius = ctx->params.default_particle_radius;
                p->radius_sq = p->radius * p->radius;
                p->m = ctx->params.default_particle_mass;
            }
            fclose(fp);
            return 0;
        }
        // ...existing code for other formats...
        default:
            fclose(fp);
            return -1;
    }
}
