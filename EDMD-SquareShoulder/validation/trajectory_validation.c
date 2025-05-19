#include "trajectory_validation.h"
#include <stdlib.h>
#include <math.h>

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
int compare_trajectories(const char* ref_path, const char* test_path, double tolerance) {
    if (!ref_path || !test_path) return 0;
    
    // Use default tolerance if not specified
    if (tolerance <= 0.0) tolerance = 1e-9;
    
    FILE* ref_file = fopen(ref_path, "r");
    if (!ref_file) {
        fprintf(stderr, "Error: Could not open reference trajectory file: %s\n", ref_path);
        return 0;
    }
    
    FILE* test_file = fopen(test_path, "r");
    if (!test_file) {
        fprintf(stderr, "Error: Could not open test trajectory file: %s\n", test_path);
        fclose(ref_file);
        return 0;
    }
    
    int match = 1;  // Assume trajectories match until proven otherwise
    char ref_line[1024], test_line[1024];
    int line_count = 0;
    
    // Read header lines and check they match
    if (fgets(ref_line, sizeof(ref_line), ref_file) && 
        fgets(test_line, sizeof(test_line), test_file)) {
        int ref_count, test_count;
        sscanf(ref_line, "%d", &ref_count);
        sscanf(test_line, "%d", &test_count);
        
        if (ref_count != test_count) {
            fprintf(stderr, "Error: Particle count mismatch - ref: %d, test: %d\n", 
                    ref_count, test_count);
            match = 0;
        }
    }
    
    // Reset file pointers to start of file
    rewind(ref_file);
    rewind(test_file);
    
    // Compare timesteps
    double t_ref, t_test;
    int frame = 0;
    
    while (match && fscanf(ref_file, "%lf", &t_ref) == 1 && 
           fscanf(test_file, "%lf", &t_test) == 1) {
        frame++;
        line_count++;
        
        // Check if timesteps match
        if (fabs(t_ref - t_test) > tolerance) {
            fprintf(stderr, "Error: Time mismatch at frame %d - ref: %.12f, test: %.12f\n", 
                    frame, t_ref, t_test);
            match = 0;
            break;
        }
        
        // Read the number of particles
        int num_particles;
        if (fscanf(ref_file, "%d", &num_particles) != 1) {
            fprintf(stderr, "Error: Could not read particle count from reference file at frame %d\n", 
                    frame);
            match = 0;
            break;
        }
        
        int test_num_particles;
        if (fscanf(test_file, "%d", &test_num_particles) != 1) {
            fprintf(stderr, "Error: Could not read particle count from test file at frame %d\n", 
                    frame);
            match = 0;
            break;
        }
        
        if (num_particles != test_num_particles) {
            fprintf(stderr, "Error: Particle count mismatch at frame %d - ref: %d, test: %d\n", 
                    frame, num_particles, test_num_particles);
            match = 0;
            break;
        }
        
        // Compare particle positions
        for (int i = 0; i < num_particles; i++) {
            double x_ref, y_ref, x_test, y_test;
            
            if (fscanf(ref_file, "%lf %lf", &x_ref, &y_ref) != 2) {
                fprintf(stderr, "Error: Could not read particle %d position from reference file at frame %d\n", 
                        i, frame);
                match = 0;
                break;
            }
            
            if (fscanf(test_file, "%lf %lf", &x_test, &y_test) != 2) {
                fprintf(stderr, "Error: Could not read particle %d position from test file at frame %d\n", 
                        i, frame);
                match = 0;
                break;
            }
            
            // Check if positions match within tolerance
            if (fabs(x_ref - x_test) > tolerance || fabs(y_ref - y_test) > tolerance) {
                fprintf(stderr, "Error: Position mismatch for particle %d at frame %d\n", i, frame);
                fprintf(stderr, "  Ref:  (%.12f, %.12f)\n", x_ref, y_ref);
                fprintf(stderr, "  Test: (%.12f, %.12f)\n", x_test, y_test);
                fprintf(stderr, "  Diff: (%.12e, %.12e)\n", 
                        fabs(x_ref - x_test), fabs(y_ref - y_test));
                match = 0;
                break;
            }
            
            line_count++;
        }
    }
    
    // Check if both files have the same number of frames
    int ref_eof = feof(ref_file);
    int test_eof = feof(test_file);
    
    if (match && (ref_eof != test_eof)) {
        if (!ref_eof) {
            fprintf(stderr, "Error: Reference file has more frames than test file\n");
        } else {
            fprintf(stderr, "Error: Test file has more frames than reference file\n");
        }
        match = 0;
    }
    
    // Clean up
    fclose(ref_file);
    fclose(test_file);
    
    if (match) {
        printf("Trajectories match within tolerance of %.12e\n", tolerance);
        printf("Compared %d frames with %d total lines\n", frame, line_count);
    }
    
    return match;
}
