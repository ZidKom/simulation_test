#include "lammps_loader.h"
#include "../core/sim_context.h" // Added for aligned_alloc
#include "../utils/pbc.h" // Added for pbc_wrap_position
#include <string.h> // Added for memset
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define LAMMPS_HEADER_LINES 8

int load_lammps_data(SimContext* ctx, const char* path) {
    FILE* fp = fopen(path, "r");
    if (!fp) {
        if (ctx->hooks.log_error) ctx->hooks.log_error("Couldn't open LAMMPS file: %s", path);
        return -1;
    }

    char line[256];
    int atom_style = 0;  // 1=atomic, 2=full
    double xlo = 0, xhi = 0, ylo = 0, yhi = 0;
    int num_particles = 0;
    
    // Parse header
    while (fgets(line, sizeof(line), fp)) {
        if (strstr(line, "atoms")) {
            sscanf(line, "%d", &num_particles);
            ctx->num_particles = num_particles;
        }
        if (strstr(line, "xlo xhi")) sscanf(line, "%lf %lf", &xlo, &xhi);
        if (strstr(line, "ylo yhi")) sscanf(line, "%lf %lf", &ylo, &yhi);
        if (strstr(line, "Atoms #")) {
            atom_style = (strstr(line, "full")) ? 2 : 1;
            break;
        }
    }

    // Set box dimensions
    ctx->xsize = xhi - xlo;
    ctx->ysize = yhi - ylo;
    
    // Allocate particles
    if (ctx->particles) {
        free(ctx->particles);  // Free existing particles if any
    }
    ctx->particles = (Particle*)aligned_alloc(64, ctx->num_particles * sizeof(Particle));
    if (!ctx->particles) {
        if (ctx->hooks.log_error) ctx->hooks.log_error("Memory allocation failed");
        fclose(fp);
        return -1;
    }
    memset(ctx->particles, 0, ctx->num_particles * sizeof(Particle)); // Zero out the allocated memory
    
    // Skip header lines after "Atoms" section
    fgets(line, sizeof(line), fp); // Skip empty line after "Atoms"
    
    // Parse atoms section
    for (int i = 0; i < ctx->num_particles; i++) {
        if (!fgets(line, sizeof(line), fp)) {
            if (ctx->hooks.log_error) ctx->hooks.log_error("Unexpected end of file while reading particles");
            fclose(fp);
            return -1;
        }
        
        Particle* p = &ctx->particles[i];
        p->id = i;
        
        if (atom_style == 2) { // "full" style: id mol type q x y z
            int type;
            sscanf(line, "%*d %*d %d %*f %lf %lf", 
                  &type, &p->x, &p->y);
            p->type = type;
        } else { // "atomic" style: id type x y z
            int type;
            sscanf(line, "%*d %d %lf %lf", &type, &p->x, &p->y);
            p->type = type;
        }
        
        // Convert from LAMMPS units to EDMD units
        p->x = (p->x - xlo) * ctx->params.sigma;
        p->y = (p->y - ylo) * ctx->params.sigma;
        
        // Set default physical properties based on type
        p->radius = ctx->params.default_particle_radius;
        p->radius_sq = p->radius * p->radius;
        p->m = ctx->params.default_particle_mass;
        
        // Initialize velocity (could be read from LAMMPS if available)
        p->vx = 0.0;
        p->vy = 0.0;
        
        // Ensure particle is within simulation box
        pbc_wrap_position(&p->x, &p->y, ctx->xsize, ctx->ysize);
    }
    
    fclose(fp);
    
    // Log successful loading
    if (ctx->hooks.log_error) { // Use log_error as a temporary replacement for log_info
        ctx->hooks.log_error("Loaded %d particles from LAMMPS file: %s", ctx->num_particles, path);
        ctx->hooks.log_error("Box dimensions: %.4f x %.4f", ctx->xsize, ctx->ysize);
    }
    return 0;
}
