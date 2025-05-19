#ifndef TRAJECTORY_VALIDATION_H
#define TRAJECTORY_VALIDATION_H

#include <stdio.h>

/**
 * @brief Compares two trajectory files to check for differences
 * 
 * This function reads trajectory files with positions of particles at different
 * time steps and compares them within a specified tolerance. This is useful for
 * regression testing and validation against reference implementations.
 * 
 * @param ref_path Path to the reference trajectory file
 * @param test_path Path to the test trajectory file
 * @param tolerance Maximum allowed difference in positions (defaults to 1e-9)
 * @return 1 if trajectories match within tolerance, 0 otherwise
 */
int compare_trajectories(const char* ref_path, const char* test_path, double tolerance);

#endif // TRAJECTORY_VALIDATION_H
