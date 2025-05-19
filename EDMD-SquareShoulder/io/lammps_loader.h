#ifndef LAMMPS_LOADER_H
#define LAMMPS_LOADER_H

#include "../core/sim_context.h"

// Change the return type to int for error reporting
int load_lammps_data(SimContext* ctx, const char* path);

#endif // LAMMPS_LOADER_H
