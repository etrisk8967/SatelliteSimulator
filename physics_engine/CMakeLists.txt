cmake_minimum_required(VERSION 3.10)
project(PhysicsEngine)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_TOOLCHAIN_FILE "/home/netrisk/SatelliteSimulator/vcpkg/scripts/buildsystems/vcpkg.cmake" CACHE STRING "Vcpkg toolchain file")

# Find VCPKG-installed packages
find_package(Eigen3 CONFIG REQUIRED)
find_package(Boost REQUIRED)
find_package(spdlog CONFIG REQUIRED)

# Include directories
include_directories(include)

# Gather source files
file(GLOB SOURCES "src/*.cpp")

# Create executable
add_executable(PhysicsEngine ${SOURCES})

# Link libraries
target_link_libraries(PhysicsEngine PRIVATE Eigen3::Eigen)
target_link_libraries(PhysicsEngine PRIVATE Boost::boost)
target_link_libraries(PhysicsEngine PRIVATE spdlog::spdlog)