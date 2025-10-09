# ===========================
# Makefile for spatpg
# ===========================

# Compiler: prefer system h5c++ if available; fallback to g++
H5CXX  := $(shell command -v /usr/bin/h5c++ 2>/dev/null)
CXX    := $(if $(H5CXX),$(H5CXX),g++)

# Compiler flags
CXXFLAGS := -O3 -std=c++11 -Wall

# Linker flags (math + GSL + HDF5)
LDLIBS := -lm -lgsl -lgslcblas

# If h5c++ is not available, add HDF5 include/lib paths manually
ifeq ($(CXX),g++)
    HDF5_INC := -I/usr/include/hdf5/serial
    HDF5_LIB := -L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5_cpp -lhdf5
else
    HDF5_INC :=
    HDF5_LIB :=
endif

# Target binary
TARGET := spatpg

# Source and object files
SRC := main.C func.C mcmc.C
OBJ := $(SRC:.C=.o)

# Default rule
all: $(TARGET)

# Linking rule
$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDLIBS) $(HDF5_INC) $(HDF5_LIB)

# Compilation rule
%.o: %.C
	$(CXX) $(CXXFLAGS) $(HDF5_INC) -c $< -o $@

# Clean up
clean:
	rm -f $(OBJ) $(TARGET)

# Print detected compiler
info:
	@echo "Compiler: $(CXX)"
	@echo "Flags: $(CXXFLAGS)"
	@echo "Libraries: $(LDLIBS) $(HDF5_LIB)"
