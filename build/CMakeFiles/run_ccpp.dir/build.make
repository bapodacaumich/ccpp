# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.30

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
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
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/brandon/cplusplusport/ccpp

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/brandon/cplusplusport/ccpp/build

# Include any dependencies generated for this target.
include CMakeFiles/run_ccpp.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/run_ccpp.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/run_ccpp.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/run_ccpp.dir/flags.make

CMakeFiles/run_ccpp.dir/src/cone_camera.cpp.o: CMakeFiles/run_ccpp.dir/flags.make
CMakeFiles/run_ccpp.dir/src/cone_camera.cpp.o: /home/brandon/cplusplusport/ccpp/src/cone_camera.cpp
CMakeFiles/run_ccpp.dir/src/cone_camera.cpp.o: CMakeFiles/run_ccpp.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/brandon/cplusplusport/ccpp/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/run_ccpp.dir/src/cone_camera.cpp.o"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/run_ccpp.dir/src/cone_camera.cpp.o -MF CMakeFiles/run_ccpp.dir/src/cone_camera.cpp.o.d -o CMakeFiles/run_ccpp.dir/src/cone_camera.cpp.o -c /home/brandon/cplusplusport/ccpp/src/cone_camera.cpp

CMakeFiles/run_ccpp.dir/src/cone_camera.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/run_ccpp.dir/src/cone_camera.cpp.i"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/brandon/cplusplusport/ccpp/src/cone_camera.cpp > CMakeFiles/run_ccpp.dir/src/cone_camera.cpp.i

CMakeFiles/run_ccpp.dir/src/cone_camera.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/run_ccpp.dir/src/cone_camera.cpp.s"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/brandon/cplusplusport/ccpp/src/cone_camera.cpp -o CMakeFiles/run_ccpp.dir/src/cone_camera.cpp.s

CMakeFiles/run_ccpp.dir/src/cuda_kernels.cu.o: CMakeFiles/run_ccpp.dir/flags.make
CMakeFiles/run_ccpp.dir/src/cuda_kernels.cu.o: CMakeFiles/run_ccpp.dir/includes_CUDA.rsp
CMakeFiles/run_ccpp.dir/src/cuda_kernels.cu.o: /home/brandon/cplusplusport/ccpp/src/cuda_kernels.cu
CMakeFiles/run_ccpp.dir/src/cuda_kernels.cu.o: CMakeFiles/run_ccpp.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/brandon/cplusplusport/ccpp/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CUDA object CMakeFiles/run_ccpp.dir/src/cuda_kernels.cu.o"
	/usr/local/cuda/bin/nvcc -forward-unknown-to-host-compiler $(CUDA_DEFINES) $(CUDA_INCLUDES) $(CUDA_FLAGS) -MD -MT CMakeFiles/run_ccpp.dir/src/cuda_kernels.cu.o -MF CMakeFiles/run_ccpp.dir/src/cuda_kernels.cu.o.d -x cu -c /home/brandon/cplusplusport/ccpp/src/cuda_kernels.cu -o CMakeFiles/run_ccpp.dir/src/cuda_kernels.cu.o

CMakeFiles/run_ccpp.dir/src/cuda_kernels.cu.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CUDA source to CMakeFiles/run_ccpp.dir/src/cuda_kernels.cu.i"
	$(CMAKE_COMMAND) -E cmake_unimplemented_variable CMAKE_CUDA_CREATE_PREPROCESSED_SOURCE

CMakeFiles/run_ccpp.dir/src/cuda_kernels.cu.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CUDA source to assembly CMakeFiles/run_ccpp.dir/src/cuda_kernels.cu.s"
	$(CMAKE_COMMAND) -E cmake_unimplemented_variable CMAKE_CUDA_CREATE_ASSEMBLY_SOURCE

CMakeFiles/run_ccpp.dir/src/obs.cpp.o: CMakeFiles/run_ccpp.dir/flags.make
CMakeFiles/run_ccpp.dir/src/obs.cpp.o: /home/brandon/cplusplusport/ccpp/src/obs.cpp
CMakeFiles/run_ccpp.dir/src/obs.cpp.o: CMakeFiles/run_ccpp.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/brandon/cplusplusport/ccpp/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/run_ccpp.dir/src/obs.cpp.o"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/run_ccpp.dir/src/obs.cpp.o -MF CMakeFiles/run_ccpp.dir/src/obs.cpp.o.d -o CMakeFiles/run_ccpp.dir/src/obs.cpp.o -c /home/brandon/cplusplusport/ccpp/src/obs.cpp

CMakeFiles/run_ccpp.dir/src/obs.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/run_ccpp.dir/src/obs.cpp.i"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/brandon/cplusplusport/ccpp/src/obs.cpp > CMakeFiles/run_ccpp.dir/src/obs.cpp.i

CMakeFiles/run_ccpp.dir/src/obs.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/run_ccpp.dir/src/obs.cpp.s"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/brandon/cplusplusport/ccpp/src/obs.cpp -o CMakeFiles/run_ccpp.dir/src/obs.cpp.s

CMakeFiles/run_ccpp.dir/src/rrtz.cpp.o: CMakeFiles/run_ccpp.dir/flags.make
CMakeFiles/run_ccpp.dir/src/rrtz.cpp.o: /home/brandon/cplusplusport/ccpp/src/rrtz.cpp
CMakeFiles/run_ccpp.dir/src/rrtz.cpp.o: CMakeFiles/run_ccpp.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/brandon/cplusplusport/ccpp/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/run_ccpp.dir/src/rrtz.cpp.o"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/run_ccpp.dir/src/rrtz.cpp.o -MF CMakeFiles/run_ccpp.dir/src/rrtz.cpp.o.d -o CMakeFiles/run_ccpp.dir/src/rrtz.cpp.o -c /home/brandon/cplusplusport/ccpp/src/rrtz.cpp

CMakeFiles/run_ccpp.dir/src/rrtz.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/run_ccpp.dir/src/rrtz.cpp.i"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/brandon/cplusplusport/ccpp/src/rrtz.cpp > CMakeFiles/run_ccpp.dir/src/rrtz.cpp.i

CMakeFiles/run_ccpp.dir/src/rrtz.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/run_ccpp.dir/src/rrtz.cpp.s"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/brandon/cplusplusport/ccpp/src/rrtz.cpp -o CMakeFiles/run_ccpp.dir/src/rrtz.cpp.s

CMakeFiles/run_ccpp.dir/src/station.cpp.o: CMakeFiles/run_ccpp.dir/flags.make
CMakeFiles/run_ccpp.dir/src/station.cpp.o: /home/brandon/cplusplusport/ccpp/src/station.cpp
CMakeFiles/run_ccpp.dir/src/station.cpp.o: CMakeFiles/run_ccpp.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/brandon/cplusplusport/ccpp/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/run_ccpp.dir/src/station.cpp.o"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/run_ccpp.dir/src/station.cpp.o -MF CMakeFiles/run_ccpp.dir/src/station.cpp.o.d -o CMakeFiles/run_ccpp.dir/src/station.cpp.o -c /home/brandon/cplusplusport/ccpp/src/station.cpp

CMakeFiles/run_ccpp.dir/src/station.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/run_ccpp.dir/src/station.cpp.i"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/brandon/cplusplusport/ccpp/src/station.cpp > CMakeFiles/run_ccpp.dir/src/station.cpp.i

CMakeFiles/run_ccpp.dir/src/station.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/run_ccpp.dir/src/station.cpp.s"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/brandon/cplusplusport/ccpp/src/station.cpp -o CMakeFiles/run_ccpp.dir/src/station.cpp.s

CMakeFiles/run_ccpp.dir/src/utils.cpp.o: CMakeFiles/run_ccpp.dir/flags.make
CMakeFiles/run_ccpp.dir/src/utils.cpp.o: /home/brandon/cplusplusport/ccpp/src/utils.cpp
CMakeFiles/run_ccpp.dir/src/utils.cpp.o: CMakeFiles/run_ccpp.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/brandon/cplusplusport/ccpp/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object CMakeFiles/run_ccpp.dir/src/utils.cpp.o"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/run_ccpp.dir/src/utils.cpp.o -MF CMakeFiles/run_ccpp.dir/src/utils.cpp.o.d -o CMakeFiles/run_ccpp.dir/src/utils.cpp.o -c /home/brandon/cplusplusport/ccpp/src/utils.cpp

CMakeFiles/run_ccpp.dir/src/utils.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/run_ccpp.dir/src/utils.cpp.i"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/brandon/cplusplusport/ccpp/src/utils.cpp > CMakeFiles/run_ccpp.dir/src/utils.cpp.i

CMakeFiles/run_ccpp.dir/src/utils.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/run_ccpp.dir/src/utils.cpp.s"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/brandon/cplusplusport/ccpp/src/utils.cpp -o CMakeFiles/run_ccpp.dir/src/utils.cpp.s

CMakeFiles/run_ccpp.dir/src/viewpoint_generator.cpp.o: CMakeFiles/run_ccpp.dir/flags.make
CMakeFiles/run_ccpp.dir/src/viewpoint_generator.cpp.o: /home/brandon/cplusplusport/ccpp/src/viewpoint_generator.cpp
CMakeFiles/run_ccpp.dir/src/viewpoint_generator.cpp.o: CMakeFiles/run_ccpp.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/brandon/cplusplusport/ccpp/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object CMakeFiles/run_ccpp.dir/src/viewpoint_generator.cpp.o"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/run_ccpp.dir/src/viewpoint_generator.cpp.o -MF CMakeFiles/run_ccpp.dir/src/viewpoint_generator.cpp.o.d -o CMakeFiles/run_ccpp.dir/src/viewpoint_generator.cpp.o -c /home/brandon/cplusplusport/ccpp/src/viewpoint_generator.cpp

CMakeFiles/run_ccpp.dir/src/viewpoint_generator.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/run_ccpp.dir/src/viewpoint_generator.cpp.i"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/brandon/cplusplusport/ccpp/src/viewpoint_generator.cpp > CMakeFiles/run_ccpp.dir/src/viewpoint_generator.cpp.i

CMakeFiles/run_ccpp.dir/src/viewpoint_generator.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/run_ccpp.dir/src/viewpoint_generator.cpp.s"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/brandon/cplusplusport/ccpp/src/viewpoint_generator.cpp -o CMakeFiles/run_ccpp.dir/src/viewpoint_generator.cpp.s

CMakeFiles/run_ccpp.dir/src/run_ccpp.cpp.o: CMakeFiles/run_ccpp.dir/flags.make
CMakeFiles/run_ccpp.dir/src/run_ccpp.cpp.o: /home/brandon/cplusplusport/ccpp/src/run_ccpp.cpp
CMakeFiles/run_ccpp.dir/src/run_ccpp.cpp.o: CMakeFiles/run_ccpp.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/brandon/cplusplusport/ccpp/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building CXX object CMakeFiles/run_ccpp.dir/src/run_ccpp.cpp.o"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/run_ccpp.dir/src/run_ccpp.cpp.o -MF CMakeFiles/run_ccpp.dir/src/run_ccpp.cpp.o.d -o CMakeFiles/run_ccpp.dir/src/run_ccpp.cpp.o -c /home/brandon/cplusplusport/ccpp/src/run_ccpp.cpp

CMakeFiles/run_ccpp.dir/src/run_ccpp.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/run_ccpp.dir/src/run_ccpp.cpp.i"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/brandon/cplusplusport/ccpp/src/run_ccpp.cpp > CMakeFiles/run_ccpp.dir/src/run_ccpp.cpp.i

CMakeFiles/run_ccpp.dir/src/run_ccpp.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/run_ccpp.dir/src/run_ccpp.cpp.s"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/brandon/cplusplusport/ccpp/src/run_ccpp.cpp -o CMakeFiles/run_ccpp.dir/src/run_ccpp.cpp.s

# Object files for target run_ccpp
run_ccpp_OBJECTS = \
"CMakeFiles/run_ccpp.dir/src/cone_camera.cpp.o" \
"CMakeFiles/run_ccpp.dir/src/cuda_kernels.cu.o" \
"CMakeFiles/run_ccpp.dir/src/obs.cpp.o" \
"CMakeFiles/run_ccpp.dir/src/rrtz.cpp.o" \
"CMakeFiles/run_ccpp.dir/src/station.cpp.o" \
"CMakeFiles/run_ccpp.dir/src/utils.cpp.o" \
"CMakeFiles/run_ccpp.dir/src/viewpoint_generator.cpp.o" \
"CMakeFiles/run_ccpp.dir/src/run_ccpp.cpp.o"

# External object files for target run_ccpp
run_ccpp_EXTERNAL_OBJECTS =

run_ccpp: CMakeFiles/run_ccpp.dir/src/cone_camera.cpp.o
run_ccpp: CMakeFiles/run_ccpp.dir/src/cuda_kernels.cu.o
run_ccpp: CMakeFiles/run_ccpp.dir/src/obs.cpp.o
run_ccpp: CMakeFiles/run_ccpp.dir/src/rrtz.cpp.o
run_ccpp: CMakeFiles/run_ccpp.dir/src/station.cpp.o
run_ccpp: CMakeFiles/run_ccpp.dir/src/utils.cpp.o
run_ccpp: CMakeFiles/run_ccpp.dir/src/viewpoint_generator.cpp.o
run_ccpp: CMakeFiles/run_ccpp.dir/src/run_ccpp.cpp.o
run_ccpp: CMakeFiles/run_ccpp.dir/build.make
run_ccpp: CMakeFiles/run_ccpp.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/home/brandon/cplusplusport/ccpp/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Linking CXX executable run_ccpp"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/run_ccpp.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/run_ccpp.dir/build: run_ccpp
.PHONY : CMakeFiles/run_ccpp.dir/build

CMakeFiles/run_ccpp.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/run_ccpp.dir/cmake_clean.cmake
.PHONY : CMakeFiles/run_ccpp.dir/clean

CMakeFiles/run_ccpp.dir/depend:
	cd /home/brandon/cplusplusport/ccpp/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/brandon/cplusplusport/ccpp /home/brandon/cplusplusport/ccpp /home/brandon/cplusplusport/ccpp/build /home/brandon/cplusplusport/ccpp/build /home/brandon/cplusplusport/ccpp/build/CMakeFiles/run_ccpp.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : CMakeFiles/run_ccpp.dir/depend

