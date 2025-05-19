CC = gcc
CFLAGS = -O3 -mavx2 -Wall -Wextra -fopenmp
LDFLAGS = -lm -lz

# Debug flags
DEBUG_FLAGS = -g -DDEBUG

# Validation flags
VALIDATION_FLAGS = -DVALIDATION_MODE=1

# Source directories
CORE_DIR = core
PHYSICS_DIR = physics
UTILS_DIR = utils
IO_DIR = io
VALIDATION_DIR = validation
TESTS_DIR = tests
BENCHMARK_DIR = benchmark

# Source files
CORE_SRC = $(wildcard $(CORE_DIR)/*.c) $(wildcard $(CORE_DIR)/event_system/*.c)
PHYSICS_SRC = $(wildcard $(PHYSICS_DIR)/*.c)
UTILS_SRC = $(filter-out $(UTILS_DIR)/pbc_new.c, $(wildcard $(UTILS_DIR)/*.c))
IO_SRC = $(wildcard $(IO_DIR)/*.c)
IO_OBJ = $(IO_SRC:.c=.o)
VALIDATION_SRC = $(filter-out $(VALIDATION_DIR)/energy_checks.c, $(filter-out $(VALIDATION_DIR)/energy_validation.c, $(wildcard $(VALIDATION_DIR)/*.c))) # Simplified to include all .c files in validation dir
# EVENT_SRC = $(wildcard event/*.c) # Define if there's an 'event' module and SRC_EVENT was used

SOLVE_PHYSICS_SRC = solve_physics.c
MAIN_SRC = main.c
BENCHMARK_SRC = benchmark/benchmark.c # Ensure this is defined

# Standardize Object File Definitions (OBJ_MODULE -> MODULE_OBJ)
# Add SOLVE_PHYSICS_OBJ
# Remove HOOKS_OBJ

CORE_OBJ = $(CORE_SRC:.c=.o)
PHYSICS_OBJ = $(PHYSICS_SRC:.c=.o)
UTILS_OBJ = $(UTILS_SRC:.c=.o)
IO_OBJ = $(IO_SRC:.c=.o)
VALIDATION_OBJ = $(VALIDATION_SRC:.c=.o) # This will include validation/hooks.o
# EVENT_OBJ = $(EVENT_SRC:.c=.o) # Define if EVENT_SRC is used

# Remove HOOKS_OBJ if it was defined, e.g.:
# HOOKS_OBJ = $(HOOKS_SRC:.c=.o) # This line should be removed

SOLVE_PHYSICS_OBJ = $(SOLVE_PHYSICS_SRC:.c=.o)
MAIN_OBJ = $(MAIN_SRC:.c=.o)
BENCHMARK_OBJ = $(BENCHMARK_SRC:.c=.o) # Ensure this is defined

# Executable
TARGET = edmd_sim

# Test targets - Standardized definitions
PHYSICS_UNIT_TEST_EXE = tests/physics_unit_tests
COLLISION_TEST_EXE = tests/collision_test
EVENT_TEST_EXE = tests/event_test
SHOULDER_COLLISION_TEST_EXE = tests/shoulder_collision_test
FULL_VALIDATION_TEST_EXE = tests/full_validation

# Aggregate of test executables
TEST_TARGETS = $(PHYSICS_UNIT_TEST_EXE) $(COLLISION_TEST_EXE) $(EVENT_TEST_EXE) $(SHOULDER_COLLISION_TEST_EXE) $(FULL_VALIDATION_TEST_EXE)

# Object files for test main sources
PHYSICS_UNIT_TEST_MAIN_OBJ = tests/physics_unit_tests.o
COLLISION_TEST_MAIN_OBJ = tests/collision_test.o
EVENT_TEST_MAIN_OBJ = tests/event_test.o
SHOULDER_COLLISION_TEST_MAIN_OBJ = tests/shoulder_collision_test.o
FULL_VALIDATION_TEST_MAIN_OBJ = tests/full_validation.o

TEST_MAIN_OBJS = $(PHYSICS_UNIT_TEST_MAIN_OBJ) $(COLLISION_TEST_MAIN_OBJ) $(EVENT_TEST_MAIN_OBJ) $(SHOULDER_COLLISION_TEST_MAIN_OBJ) $(FULL_VALIDATION_TEST_MAIN_OBJ)

# Common Objects for linking main target and tests
COMMON_OBJS = $(CORE_OBJ) $(PHYSICS_OBJ) $(UTILS_OBJ) $(IO_OBJ) $(VALIDATION_OBJ) $(SOLVE_PHYSICS_OBJ) # Add $(EVENT_OBJ) if used

# Default target
all: $(TARGET)

# Main executable
$(TARGET): $(COMMON_OBJS) $(BENCHMARK_OBJ) $(MAIN_OBJ)
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

# New rule for physics_unit_tests (becomes standard form)
$(PHYSICS_UNIT_TEST_EXE): $(PHYSICS_UNIT_TEST_MAIN_OBJ) $(COMMON_OBJS)
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

# Standardized rules for other tests
$(COLLISION_TEST_EXE): $(COLLISION_TEST_MAIN_OBJ) $(COMMON_OBJS)
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

$(EVENT_TEST_EXE): $(EVENT_TEST_MAIN_OBJ) $(COMMON_OBJS)
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

$(SHOULDER_COLLISION_TEST_EXE): $(SHOULDER_COLLISION_TEST_MAIN_OBJ) $(COMMON_OBJS)
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

$(FULL_VALIDATION_TEST_EXE): $(FULL_VALIDATION_TEST_MAIN_OBJ) $(COMMON_OBJS)
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)


# Run all tests
test: $(TEST_TARGETS)
	./$(PHYSICS_UNIT_TEST_EXE)
	./$(COLLISION_TEST_EXE)
	./$(EVENT_TEST_EXE)
	./$(SHOULDER_COLLISION_TEST_EXE)
# Note: full_validation might be run separately via `make full_test`

.PHONY: full_test
full_test: $(FULL_VALIDATION_TEST_EXE)
	./$<

# Debug build
debug: CFLAGS += $(DEBUG_FLAGS)
debug: all

# Validation build
validate: CFLAGS += $(VALIDATION_FLAGS)
validate: all

# Clean
clean:
	rm -f $(CORE_OBJ) $(PHYSICS_OBJ) $(UTILS_OBJ) $(IO_OBJ) $(VALIDATION_OBJ) $(MAIN_OBJ) $(SOLVE_PHYSICS_OBJ) $(BENCHMARK_OBJ)
	rm -f $(TARGET)
	rm -f $(TEST_MAIN_OBJS) $(TEST_TARGETS)
	rm -f tests/*.o # Clean any other .o files in tests, if any
	rm -f *.o # Clean root .o files if any (e.g. simulation_hooks.o if not part of a module)
	rm -rf $(DEPDIR)

# Generic rule for compiling .c to .o
DEPDIR = .deps

%.o: %.c
	@mkdir -p $(dir $@) # Ensure directory exists for .o file, including subdirs
	$(CC) $(CFLAGS) -MMD -MP -c $< -o $@

# Specific rule for test objects to ensure they go into tests/ dir if not handled by generic rule with $(dir $@)
tests/%.o: tests/%.c
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) -MMD -MP -c $< -o $@

# Dependencies
# The following directory prerequisites are good for `mkdir -p` but not as direct dependencies for .o files.
# $(CORE_OBJ): | $(CORE_DIR) $(CORE_DIR)/event_system
# $(PHYSICS_OBJ): | $(PHYSICS_DIR)
# $(UTILS_OBJ): | $(UTILS_DIR)
# $(IO_OBJ): | $(IO_DIR)
# $(VALIDATION_OBJ): | $(VALIDATION_DIR)
# $(TEST_OBJ): | $(TESTS_DIR) # TEST_OBJ is not defined, likely meant TEST_MAIN_OBJS
# $(BENCHMARK_OBJ): | $(BENCHMARK_DIR)

# Ensure directories exist (can be removed if `mkdir -p $(dir $@)` is used in rules)
# $(CORE_DIR) $(CORE_DIR)/event_system $(PHYSICS_DIR) $(UTILS_DIR) $(IO_DIR) $(VALIDATION_DIR) $(TESTS_DIR) $(BENCHMARK_DIR):
# 	mkdir -p $@

io/particle_loader.o: io/particle_loader.c io/particle_loader.h core/sim_context.h
	$(CC) $(CFLAGS) -c $< -o $@

# Phony target for running all tests
.PHONY: tests
tests: $(TEST_TARGETS)

# Clean rule - ensure it covers new object files and test executables
# clean:
#	rm -f $(TARGET) $(TEST_TARGETS)
#	rm -f $(MAIN_OBJ) $(COMMON_OBJS) $(BENCHMARK_OBJ_IF_USED) $(TEST_MAIN_OBJS)
#	rm -f core/*.o core/event_system/*.o physics/*.o utils/*.o io/*.o validation/*.o # etc. for all source dirs
#	rm -f tests/*.o
#	rm -f *.o # For any root .o files
.PHONY: all test debug validate clean

# Include dependency files
-include $(wildcard $(DEPDIR)/*.d)
