# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.19

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
CMAKE_COMMAND = /opt/local/bin/cmake

# The command to remove a file.
RM = /opt/local/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/zehaozhou/Desktop/project/toyfssh_cpp/src

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/zehaozhou/Desktop/project/toyfssh_cpp/build

# Include any dependencies generated for this target.
include CMakeFiles/toyFSSH.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/toyFSSH.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/toyFSSH.dir/flags.make

CMakeFiles/toyFSSH.dir/main.cpp.o: CMakeFiles/toyFSSH.dir/flags.make
CMakeFiles/toyFSSH.dir/main.cpp.o: /Users/zehaozhou/Desktop/project/toyfssh_cpp/src/main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/zehaozhou/Desktop/project/toyfssh_cpp/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/toyFSSH.dir/main.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/toyFSSH.dir/main.cpp.o -c /Users/zehaozhou/Desktop/project/toyfssh_cpp/src/main.cpp

CMakeFiles/toyFSSH.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/toyFSSH.dir/main.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/zehaozhou/Desktop/project/toyfssh_cpp/src/main.cpp > CMakeFiles/toyFSSH.dir/main.cpp.i

CMakeFiles/toyFSSH.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/toyFSSH.dir/main.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/zehaozhou/Desktop/project/toyfssh_cpp/src/main.cpp -o CMakeFiles/toyFSSH.dir/main.cpp.s

# Object files for target toyFSSH
toyFSSH_OBJECTS = \
"CMakeFiles/toyFSSH.dir/main.cpp.o"

# External object files for target toyFSSH
toyFSSH_EXTERNAL_OBJECTS =

toyFSSH: CMakeFiles/toyFSSH.dir/main.cpp.o
toyFSSH: CMakeFiles/toyFSSH.dir/build.make
toyFSSH: CMakeFiles/toyFSSH.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/zehaozhou/Desktop/project/toyfssh_cpp/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable toyFSSH"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/toyFSSH.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/toyFSSH.dir/build: toyFSSH

.PHONY : CMakeFiles/toyFSSH.dir/build

CMakeFiles/toyFSSH.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/toyFSSH.dir/cmake_clean.cmake
.PHONY : CMakeFiles/toyFSSH.dir/clean

CMakeFiles/toyFSSH.dir/depend:
	cd /Users/zehaozhou/Desktop/project/toyfssh_cpp/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/zehaozhou/Desktop/project/toyfssh_cpp/src /Users/zehaozhou/Desktop/project/toyfssh_cpp/src /Users/zehaozhou/Desktop/project/toyfssh_cpp/build /Users/zehaozhou/Desktop/project/toyfssh_cpp/build /Users/zehaozhou/Desktop/project/toyfssh_cpp/build/CMakeFiles/toyFSSH.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/toyFSSH.dir/depend

