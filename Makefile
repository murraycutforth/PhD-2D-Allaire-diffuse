#
# Makefile for the 2D-Allaire-diffuse project
#
# External dependencies
#	- The header-only Eigen library. Download from eigen.tuxfamily.org
#	- My own 2D-HCLframework headers. Download from https://github.com/murraycutforth/2D-HCLsolver-framework
#
# Based on the file at http://hiltmon.com/
#

# Target
TARGET := 2D-Allaire-diffuse

# External header files
EIGEN := -I/home/raid/mcc74/Libraries/eigen-v3.3.4
HCLFRAMEWORK := -I/home/raid/mcc74/Documents/phd_17/2D-HCLsolver-framework

# Folders
BUILDDIR := objectfiles
SRCDIR := sourcecode
INCDIRS := headers
INCLIST := -I./headers

# Files
SOURCES := $(shell find $(SRCDIR) -type f -name *.cpp)
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.cpp=.o))

# Folder Lists
BUILDLIST := $(BUILDDIR)
EXTINCLIST := $(EIGEN) $(HCLFRAMEWORK) $(INCLIST)

# Shared Compiler Flags
OPLEVEL := -O3
CFLAGS := -Wall -c -pg -g -std=c++11 $(OPLEVEL)
LINKFLAGS := -pg -fopenmp $(OPLEVEL)

# Linking Step
$(TARGET): $(OBJECTS)
	@mkdir -p $(BUILDLIST)
	@echo "Linking..."
	@echo "Linking $(TARGET) using options: $(LINKFLAGS)"; g++ $^ $(LINKFLAGS) -o $(TARGET)
	@echo "Success!"

# Compilation Step
$(BUILDDIR)/%.o: $(SRCDIR)/%.cpp
	@mkdir -p $(BUILDLIST)
	@echo "Compiling $< using options: $(CFLAGS)"; g++ $(CFLAGS) $(EXTINCLIST) -o $@ $<

clean:
	@echo "Cleaning $(TARGET)..."; rm $(BUILDDIR)/* $(TARGET)
