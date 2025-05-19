#pragma once

#include "../core/sim_context.h"

/**
 * @brief Saves the current simulation state to a compressed file
 * 
 * @param ctx Simulation context to save
 * @param path Path to the output file
 */
void save_snapshot(SimContext* ctx, const char* path);

/**
 * @brief Loads a simulation state from a compressed file
 * 
 * @param ctx Simulation context to load into
 * @param path Path to the input file
 */
void load_snapshot(SimContext* ctx, const char* path);

/**
 * @brief Exports the current particle positions to an XYZ format file
 * 
 * @param ctx Simulation context to export
 * @param path Path to the output XYZ file
 */
void export_xyz_frame(SimContext* ctx, const char* path);

/**
 * @brief Detects the file format of a given file
 * 
 * @param filename Name of the file whose format is to be detected
 * @return int Integer representing the file format
 */
int detect_file_format(const char* filename);

/**
 * @brief Enum for different snapshot formats
 * 
 */
typedef enum {
    SNAPSHOT_FORMAT_XYZ,
    SNAPSHOT_FORMAT_EDMD,
    SNAPSHOT_FORMAT_LAMMPS,
    SNAPSHOT_FORMAT_POSSPECIAL // For pos.txt
} SnapshotFormat;
