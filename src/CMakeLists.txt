cmake_minimum_required(VERSION 3.10)

# specify the compiler to be used
# the specific file path is gotten from running `which gcc` in the terminal
set(CMAKE_C_COMPILER )
set(CMAKE_CXX_COMPILER )

# set the project name and version
project(toyFSSH VERSION 1.0)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# find the eigen3 math library
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# add the executable
add_executable(toyFSSH main.cpp)

# link to the eigen3 library
target_link_libraries(toyFSSH Eigen3::Eigen)
