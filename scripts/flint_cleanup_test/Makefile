# =============================================================================
# Makefile for a C++ program using OpenMP, FLINT, and GMP
# =============================================================================

# 1. Compiler and Flags
# -----------------------------------------------------------------------------
# Use g++ as the C++ compiler
CXX = g++

# C++ standard and compiler flags:
# -std=c++17: Use the C++17 standard (for thread_local and other modern features)
# -Wall -Wextra: Enable all major and extra warnings for better code quality
# -O2:          Enable level 2 optimizations for performance
# -fopenmp:     Enable OpenMP support for parallelization
CXXFLAGS = -std=c++17 -Wall -Wextra -O2 -fopenmp -pthread

# Linker flags:
# -fopenmp is also needed at the linking stage
LDFLAGS = -fopenmp

# Libraries to link against:
# -lflint: The FLINT library
# -lgmp:   The GMP library (a dependency for FLINT)
LDLIBS = -lflint -lgmp


# 2. Project Files
# -----------------------------------------------------------------------------
# The name of the final executable program
TARGET = flint_program

# List of your C++ source files
SOURCES = main.cpp

# Automatically generate object file names from source file names (e.g., main.cpp -> main.o)
OBJECTS = $(SOURCES:.cpp=.o)


# 3. Build Rules
# -----------------------------------------------------------------------------
# The default rule, executed when you just type "make"
# It depends on the final executable being built.
all: $(TARGET)

# Rule to link the final executable from the object files
$(TARGET): $(OBJECTS)
	@echo "==> Linking executable: $@"
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)
	@echo "==> Build complete: ./$@"

# Generic rule to compile a .cpp source file into a .o object file
%.o: %.cpp
	@echo "==> Compiling source: $<"
	$(CXX) $(CXXFLAGS) -c -o $@ $<


# 4. Utility Rules
# -----------------------------------------------------------------------------
# Rule to run the program
run: $(TARGET)
	@echo "==> Running program..."
	./$(TARGET)

# Rule to clean up the directory by removing generated files
clean:
	@echo "==> Cleaning up..."
	rm -f $(TARGET) $(OBJECTS)

# Declare targets that are not actual files to avoid conflicts
.PHONY: all clean run