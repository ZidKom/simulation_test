# EDMD-SquareShoulder Simulation

An Event-Driven Molecular Dynamics (EDMD) simulation for particles interacting via square shoulder potentials.

## Overview

This simulation implements an event-driven approach to molecular dynamics, optimized for systems with particles interacting through square shoulder potentials. The event-driven method advances the system from one event to the next rather than using fixed time steps, which is particularly efficient for discontinuous potentials.

## Key Features

- **Event-Driven Architecture**: Efficiently simulates discrete events (collisions, potential changes) without fixed time steps.
- **Hybrid Event Queue**: Combines a binned "Paul list" for near-future events with a binary search tree (BST) for far-future events.
- **Square Shoulder Potentials**: Models particles interacting through a potential with a hard core and square shoulder region.
- **Cell List Optimization**: Spatial partitioning to accelerate collision detection and neighbor finding.
- **Energy Conservation Validation**: Built-in hooks for monitoring energy conservation.

## Project Structure

- **core/**: Core simulation infrastructure
  - **sim_context.h/c**: Main simulation context and data structures
  - **particle_pool.h/c**: Memory management for particles
  - **event_system/**: Event queue implementation
    - **hybrid.h/c**: Hybrid event queue (Paul list + BST)
    - **paul_list.h/c**: Implementation of Paul binning system
    - **bst_queue.h/c**: Binary search tree for far-future events
    
- **physics/**: Physical interaction implementations
  - **collisions.h/c**: Hard-sphere collision handling
  - **shoulders.h/c**: Square shoulder potential interactions
  - **cell_interactions.h/c**: Cell boundary crossing events
  - **event_prediction.c**: Prediction logic for all event types
  
- **utils/**: Utility functions
  - **cell_list.c**: Cell list management
  - **pbc.c**: Periodic boundary condition functions
  - **vec_math.c**: Vector math operations
  
- **validation/**: Simulation validation
  - **energy_audit.c**: Energy conservation checks
  - **hooks.c**: Validation hook system
  
- **tests/**: Test cases
  - **collision_test.c**: Collision physics tests
  - **event_test.c**: Event system tests

## Compilation

To build the simulation:

```bash
gcc -O3 -mavx2 -DVALIDATION_MODE=1 \
    core/*.c physics/*.c validation/*.c tests/*.c \
    -lm -fopenmp -o edmd_sim
```

## Usage

Basic usage:

```bash
./edmd_sim --num_particles 1000 --box_size 50.0 --sim_end_time 100.0
```

For more options:

```bash
./edmd_sim --help
```

## Testing

Run the collision test with:

```bash
./edmd_sim --test collisions
```

## Technical Details

### Event Handling

The simulation handles several types of events:

1. **Collisions**: Hard-core collisions between particles
2. **Shoulder Entry/Exit**: Particles entering or exiting the shoulder potential region
3. **Cell Crossings**: Particles crossing cell boundaries in the spatial grid

### Cell List System

Particles are organized in a grid of cells to optimize neighbor finding. Events are predicted when particles cross cell boundaries, allowing efficient updates to neighbor lists.

### Energy Conservation

The system tracks kinetic energy and validates conservation within a specified tolerance. Energy can be intentionally injected during collisions to model active systems.

## License

This software is provided as-is for research and educational purposes.
