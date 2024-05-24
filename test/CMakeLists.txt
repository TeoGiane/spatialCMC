cmake_minimum_required(VERSION 3.13.0)
project(test_spatialcmc)
enable_testing()

find_package(PkgConfig REQUIRED)
set(CMAKE_CXX_FLAGS -pthread)

include(FetchContent)
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG        9c332145b71c36a5bad9688312c79184f98601ff # release-1.13
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

add_executable(test_spatialcmc $<TARGET_OBJECTS:spatialCMC>
 hierarchies.cc
)

target_include_directories(test_spatialcmc PUBLIC ${INCLUDE_PATHS})
target_link_libraries(test_spatialcmc PUBLIC ${LINK_LIBRARIES} gtest_main gtest)
target_compile_options(test_spatialcmc PUBLIC ${COMPILE_OPTIONS})

add_test(NAME runtest COMMAND test_spatialcmc)
