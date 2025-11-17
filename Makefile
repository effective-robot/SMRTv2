# Makefile for RNA-seq Mapper (refactored modular version)
#
# Build production mapper with high optimization

CXX := g++
CXXFLAGS := -O3 -std=c++17 -march=native -Wall -Wextra -pedantic
LDFLAGS := -lz

# Output binary
TARGET := mapper_v17_modular

# Source directories
SRC_DIR := src
BUILD_DIR := build

# Find all source files
SOURCES := $(shell find $(SRC_DIR) -name '*.cpp')
OBJECTS := $(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(SOURCES))
DEPS := $(OBJECTS:.o=.d)

# Default target
.PHONY: all
all: $(TARGET)

# Link final binary
$(TARGET): $(OBJECTS)
	@echo "Linking $(TARGET)..."
	@$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)
	@echo "Build complete: $(TARGET)"

# Compile source files to object files
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(dir $@)
	@echo "Compiling $<..."
	@$(CXX) $(CXXFLAGS) -MMD -MP -c $< -o $@

# Include dependency files
-include $(DEPS)

# Clean build artifacts
.PHONY: clean
clean:
	@echo "Cleaning build artifacts..."
	@rm -rf $(BUILD_DIR) $(TARGET)
	@echo "Clean complete"

# Test compilation (quick check without full optimization)
.PHONY: test-compile
test-compile: CXXFLAGS := -std=c++17 -Wall -Wextra -pedantic
test-compile: clean all
	@echo "Test compilation successful"

# Show build information
.PHONY: info
info:
	@echo "Compiler: $(CXX)"
	@echo "Flags: $(CXXFLAGS)"
	@echo "Linker flags: $(LDFLAGS)"
	@echo "Sources: $(words $(SOURCES)) files"
	@echo "Target: $(TARGET)"

# Help target
.PHONY: help
help:
	@echo "RNA-seq Mapper - Build System"
	@echo ""
	@echo "Targets:"
	@echo "  all (default)    Build optimized binary"
	@echo "  clean            Remove build artifacts"
	@echo "  test-compile     Quick compilation test"
	@echo "  info             Show build configuration"
	@echo "  help             Show this help message"
	@echo ""
	@echo "Usage:"
	@echo "  make              # Build production binary"
	@echo "  make clean        # Clean build files"
	@echo "  make test-compile # Test compilation"
