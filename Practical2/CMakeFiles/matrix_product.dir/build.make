# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list

# Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /user/2/lijia/Documents/M1/HPC/Practical2

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /user/2/lijia/Documents/M1/HPC/Practical2

# Include any dependencies generated for this target.
include CMakeFiles/matrix_product.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/matrix_product.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/matrix_product.dir/flags.make

CMakeFiles/matrix_product.dir/matrix_product.cxx.o: CMakeFiles/matrix_product.dir/flags.make
CMakeFiles/matrix_product.dir/matrix_product.cxx.o: matrix_product.cxx
	$(CMAKE_COMMAND) -E cmake_progress_report /user/2/lijia/Documents/M1/HPC/Practical2/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/matrix_product.dir/matrix_product.cxx.o"
	/opt/rh/devtoolset-8/root/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/matrix_product.dir/matrix_product.cxx.o -c /user/2/lijia/Documents/M1/HPC/Practical2/matrix_product.cxx

CMakeFiles/matrix_product.dir/matrix_product.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/matrix_product.dir/matrix_product.cxx.i"
	/opt/rh/devtoolset-8/root/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /user/2/lijia/Documents/M1/HPC/Practical2/matrix_product.cxx > CMakeFiles/matrix_product.dir/matrix_product.cxx.i

CMakeFiles/matrix_product.dir/matrix_product.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/matrix_product.dir/matrix_product.cxx.s"
	/opt/rh/devtoolset-8/root/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /user/2/lijia/Documents/M1/HPC/Practical2/matrix_product.cxx -o CMakeFiles/matrix_product.dir/matrix_product.cxx.s

CMakeFiles/matrix_product.dir/matrix_product.cxx.o.requires:
.PHONY : CMakeFiles/matrix_product.dir/matrix_product.cxx.o.requires

CMakeFiles/matrix_product.dir/matrix_product.cxx.o.provides: CMakeFiles/matrix_product.dir/matrix_product.cxx.o.requires
	$(MAKE) -f CMakeFiles/matrix_product.dir/build.make CMakeFiles/matrix_product.dir/matrix_product.cxx.o.provides.build
.PHONY : CMakeFiles/matrix_product.dir/matrix_product.cxx.o.provides

CMakeFiles/matrix_product.dir/matrix_product.cxx.o.provides.build: CMakeFiles/matrix_product.dir/matrix_product.cxx.o

# Object files for target matrix_product
matrix_product_OBJECTS = \
"CMakeFiles/matrix_product.dir/matrix_product.cxx.o"

# External object files for target matrix_product
matrix_product_EXTERNAL_OBJECTS =

matrix_product: CMakeFiles/matrix_product.dir/matrix_product.cxx.o
matrix_product: CMakeFiles/matrix_product.dir/build.make
matrix_product: CMakeFiles/matrix_product.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable matrix_product"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/matrix_product.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/matrix_product.dir/build: matrix_product
.PHONY : CMakeFiles/matrix_product.dir/build

CMakeFiles/matrix_product.dir/requires: CMakeFiles/matrix_product.dir/matrix_product.cxx.o.requires
.PHONY : CMakeFiles/matrix_product.dir/requires

CMakeFiles/matrix_product.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/matrix_product.dir/cmake_clean.cmake
.PHONY : CMakeFiles/matrix_product.dir/clean

CMakeFiles/matrix_product.dir/depend:
	cd /user/2/lijia/Documents/M1/HPC/Practical2 && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /user/2/lijia/Documents/M1/HPC/Practical2 /user/2/lijia/Documents/M1/HPC/Practical2 /user/2/lijia/Documents/M1/HPC/Practical2 /user/2/lijia/Documents/M1/HPC/Practical2 /user/2/lijia/Documents/M1/HPC/Practical2/CMakeFiles/matrix_product.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/matrix_product.dir/depend
