#include "snapshot.h"
#include "../core/sim_context.h"
#include "../core/event_system/events.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <stdint.h>
#include <endian.h>

#define SNAPSHOT_MAGIC 0xED2025   // Fixed the macro to be a valid hex constant
#define SNAPSHOT_VERSION 1

// State saving/loading stub

// Mark as unused to prevent warnings
static uint32_t calc_crc32(const void* data, size_t len) __attribute__((unused));
static uint32_t calc_crc32(const void* data, size_t len) {
    return crc32(0L, (const Bytef*)data, len);
}

void save_snapshot(SimContext* ctx, const char* path) {
    if (!ctx) return;
    gzFile f = gzopen(path, "wb");
    if (!f) { fprintf(stderr, "[SNAPSHOT] Failed to open %s\n", path); exit(EXIT_FAILURE); }
    uint32_t magic = htole32(SNAPSHOT_MAGIC);
    uint32_t version = htole32(SNAPSHOT_VERSION);
    gzwrite(f, &magic, sizeof(magic));
    gzwrite(f, &version, sizeof(version));
    gzwrite(f, &ctx->current_time, sizeof(ctx->current_time));
    gzwrite(f, &ctx->params, sizeof(ctx->params));
    gzwrite(f, &ctx->num_particles, sizeof(ctx->num_particles));
    for (int i = 0; i < ctx->num_particles; ++i) {
        gzwrite(f, &ctx->particles[i], sizeof(Particle));
    }
    // TODO: Serialize events (BST + Paul lists) if needed
    gzclose(f);
}

void load_snapshot(SimContext* ctx, const char* path) {
    if (!ctx) return;
    gzFile f = gzopen(path, "rb");
    if (!f) { fprintf(stderr, "[SNAPSHOT] Failed to open %s\n", path); exit(EXIT_FAILURE); }
    uint32_t magic, version;
    gzread(f, &magic, sizeof(magic));
    gzread(f, &version, sizeof(version));
    if (le32toh(magic) != SNAPSHOT_MAGIC || le32toh(version) != SNAPSHOT_VERSION) {
        fprintf(stderr, "[SNAPSHOT] Incompatible snapshot file\n");
        gzclose(f); exit(EXIT_FAILURE);
    }
    gzread(f, &ctx->current_time, sizeof(ctx->current_time));
    gzread(f, &ctx->params, sizeof(ctx->params));
    gzread(f, &ctx->num_particles, sizeof(ctx->num_particles));
    if (ctx->particles) free(ctx->particles);
    ctx->particles = calloc(ctx->num_particles, sizeof(Particle));
    for (int i = 0; i < ctx->num_particles; ++i) {
        gzread(f, &ctx->particles[i], sizeof(Particle));
    }
    // TODO: Deserialize events (BST + Paul lists) if needed
    gzclose(f);
}

void export_xyz_frame(SimContext* ctx, const char* path) {
    if (!ctx) return;
    FILE* f = fopen(path, "w");
    if (!f) { fprintf(stderr, "[XYZ] Failed to open %s\n", path); return; }
    fprintf(f, "%d\nComment: t=%.4f\n", ctx->num_particles, ctx->current_time);
    for (int i = 0; i < ctx->num_particles; ++i) {
        Particle* p = &ctx->particles[i];
        fprintf(f, "H %.4f %.4f 0 0.5 0.5\n", p->x, p->y);
    }
    fclose(f);
}

int detect_file_format(const char* filename) {
    if (strstr(filename, "pos.txt")) 
        return SNAPSHOT_FORMAT_POSSPECIAL;
    return SNAPSHOT_FORMAT_XYZ;
}
