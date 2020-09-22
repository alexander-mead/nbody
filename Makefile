# Makefile to compile nbody code

# Standard flags
FFLAGS = \
	-fmax-errors=4 \
	-fmax-errors=4 \
	-ffpe-trap=invalid,zero,overflow \
	-fimplicit-none \
	-std=f2008 \
	-ffree-line-length-none \
	-O3 \
	-fdefault-real-16

# Extra debugging flags
DEBUG_FLAGS = \
	-Wall \
	-fcheck=all \
	-fbounds-check \
	-fbacktrace \
	-Og

# Fortran compiler
FC = gfortran 
all: bin

# Source-code directory
SRC_DIR = src

# Build directory
BUILD_DIR = build

MOD_DIR = library/src

# Debug build directory
DEBUG_BUILD_DIR = debug_build

# Library directory
LIB_DIR = lib

# Executable directory
BIN_DIR = bin

# Objects
_OBJ = \
	precision.o \
	constants.o \
	string_operations.o \
	basic_operations.o \
	array_operations.o \
	file_info.o \
	vectors.o

# Add prefixes of build directory to objects
OBJ = $(addprefix $(BUILD_DIR)/,$(_OBJ))
DEBUG_OBJ = $(addprefix $(DEBUG_BUILD_DIR)/,$(_OBJ))

# ?
make_dirs = @mkdir -p $(@D)

# Standard rules
bin: $(BIN_DIR)/nbody

# Debugging rules
debug: FFLAGS += $(DEBUG_FLAGS)
debug: $(BIN_DIR)/nbody_debug

# Rule to make object files
$(BUILD_DIR)/%.o: $(MOD_DIR)/%.f90
	$(make_dirs)
	$(FC) -c -o $@ $< -J$(BUILD_DIR) $(LDFLAGS) $(FFLAGS)

# Rule to make executable
$(BIN_DIR)/nbody: $(OBJ) $(SRC_DIR)/nbody.f90
	@echo "\nBuilding executable.\n"
	$(make_dirs)
	$(FC) -o $@ $^ -J$(BUILD_DIR) $(LDFLAGS) $(FFLAGS)

# Rule to make debugging objects
$(DEBUG_BUILD_DIR)/%.o: $(MOD_DIR)/%.f90
	$(make_dirs)
	$(FC) -c -o $@ $< -J$(DEBUG_BUILD_DIR) $(LDFLAGS) $(FFLAGS)

# Rule to make debugging executable
$(BIN_DIR)/nbody_debug: $(DEBUG_OBJ) $(SRC_DIR)/nbody.f90
	@echo "\nBuilding debugging executable.\n"
	$(FC) -o $@ $^ -J$(DEBUG_BUILD_DIR) $(LDFLAGS) $(FFLAGS)

# Clean up
.PHONY: clean
clean:
	rm -f $(BIN_DIR)/nbody
	rm -f $(BUILD_DIR)/*.o
	rm -f $(BUILD_DIR)/*.mod
	rm -f $(SRC_DIR)/*.mod
	rm -f $(DEBUG_BUILD_DIR)/*.o
	rm -f $(DEBUG_BUILD_DIR)/*.mod
