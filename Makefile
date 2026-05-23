# ==========================================
# GRACE-DG Solver Makefile
# ==========================================

# 1. Specify the C++ compiler
# Default is usually g++ or clang++ on macOS/Linux
CXX = g++

# 2. Specify compiler flags
# -std=c++17 : Use C++17 standard
# -Wall -Wextra : Enable standard compiler warnings to catch potential issues
# -O3 : Enable maximum optimization for better numerical simulation performance
# -Iinclude : Tell the compiler to look in the 'include' directory for header files (.hpp)
CXXFLAGS = -std=c++17 -Wall -Wextra -O3 -Iinclude

# 3. Specify the target executable name
TARGET = GRACE_DG

# 4. Automatically find all source files
# The wildcard function finds all .cpp files in the src/ directory
SRCS = $(wildcard src/*.cpp)

# 5. Automatically deduce object files
# Replace the .cpp extension in the SRCS list with .o to generate the object file list
OBJS = $(SRCS:.cpp=.o)

# ==========================================
# Build Rules
# ==========================================

# Default target executed when running `make` in the terminal without arguments
all: $(TARGET)

# Linking rule: Link all .o object files to create the final executable
# $@ represents the target name ($(TARGET)), $^ represents all prerequisites ($(OBJS))
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^
	@echo "Build successful! Executable generated: $(TARGET)"

# Compilation rule: Define how to compile a single .cpp file into a .o object file
# $< represents the first prerequisite (.cpp), $@ represents the target (.o)
src/%.o: src/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean rule: Remove all generated .o intermediate files and the final executable
# Execute by running `make clean` in the terminal
clean:
	rm -f $(OBJS) $(TARGET)
	@echo "Clean complete! Intermediate files and executable removed."

# Declare phony targets to prevent conflicts if files named 'all' or 'clean' exist in the directory
.PHONY: all clean
