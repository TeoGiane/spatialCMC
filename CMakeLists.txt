cmake_minimum_required(VERSION 3.14.0)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}")
message("CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}")

find_program(CCACHE_PROGRAM ccache)
if (CCACHE_PROGRAM)
	set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CCACHE_PROGRAM}")
endif()

project(spatialCMC)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_BUILD_TYPE Release)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -funroll-loops -fopenmp -ftree-vectorize -Wno-deprecated")
set(CMAKE_CXX_FLAGS_DEBUG "-Og -fopenmp -Wno-deprecated")
set(CMAKE_FIND_PACKAGE_PREFER_CONFIG TRUE)
set(CMAKE_FIND_PACKAGE_TARGETS_GLOBAL TRUE)

# Require OpenMP
find_package(OpenMP REQUIRED)

# Include FetchContent
include(FetchContent)

# Set up FetchContent
set(FETCHCONTENT_QUIET OFF)
set(FETCHCONTENT_UPDATES_DISCONNECTED ON)

# Fetch geos library and check if it is found
# message(CHECK_START "Fetching geos library: ")
# list(APPEND CMAKE_MESSAGE_INDENT "  ")
# FetchContent_Declare(
#   geos
#   URL https://github.com/libgeos/geos/archive/refs/tags/3.11.2.zip
#   DOWNLOAD_EXTRACT_TIMESTAMP TRUE
# )
# FetchContent_MakeAvailable(geos)
# message("--- geos library found at: " ${geos_SOURCE_DIR})

# Fetch bayesmix library and check if it is found (SWITCH TO USUAL REPO WHEN UPDATED)
message(CHECK_START "Fetching bayesmix")
list(APPEND CMAKE_MESSAGE_INDENT "  ")
set(DISABLE_TESTS ON)
set(DISABLE_PLOTS ON)
set(DISABLE_EXAMPLES ON)
set(BUILD_RUN OFF)
FetchContent_Declare(
  bayesmix
  GIT_REPOSITORY "https://github.com/bayesmix-dev/bayesmix.git"
  GIT_TAG "master"
)
FetchContent_MakeAvailable(bayesmix)
message("--- bayesmix found at: " ${bayesmix_SOURCE_DIR})

# Compile new protos
file(GLOB ProtoFiles "${CMAKE_CURRENT_LIST_DIR}/src/proto/*.proto")
set(PROTO_DIR proto)
foreach(PROTO_FILE IN LISTS ProtoFiles)
  message(STATUS "protoc proto(cc): ${PROTO_FILE}")
  get_filename_component(PROTO_DIR ${PROTO_FILE} DIRECTORY)
  get_filename_component(PROTO_NAME ${PROTO_FILE} NAME_WE)
  set(PROTO_HDR ${CMAKE_CURRENT_BINARY_DIR}/${PROTO_NAME}.pb.h)
  set(PROTO_SRC ${CMAKE_CURRENT_BINARY_DIR}/${PROTO_NAME}.pb.cc)
  message(STATUS "protoc hdr: ${PROTO_HDR}")
  message(STATUS "protoc src: ${PROTO_SRC}")
  add_custom_command(
    OUTPUT ${PROTO_SRC} ${PROTO_HDR}
    COMMAND protobuf::protoc
    "--proto_path=${BAYESMIX_PROTO_PATH}"
    "--proto_path=${CMAKE_CURRENT_LIST_DIR}/src/proto"
    ${PROTO_DIRS} "--cpp_out=${PROJECT_BINARY_DIR}" ${PROTO_FILE}
    DEPENDS ${PROTO_FILE} protobuf::protoc
    COMMENT "Generate C++ protocol buffer for ${PROTO_FILE}"
    VERBATIM)
  list(APPEND PROTO_HDRS ${PROTO_HDR})
  list(APPEND PROTO_SRCS ${PROTO_SRC})
endforeach()

# Set include paths
set(INCLUDE_PATHS
  ${CMAKE_CURRENT_LIST_DIR}/src
  ${CMAKE_CURRENT_BINARY_DIR}
  ${BAYESMIX_INCLUDE_PATHS}
)

# Set link libraries
set(LINK_LIBRARIES
	pthread
  OpenMP::OpenMP_CXX
  ${BAYESMIX_LINK_LIBRARIES}
  bayesmixlib
)

# Set compile options
set(COMPILE_OPTIONS -D_REENTRANT -fPIC)

# Define spatialCMC library
add_library(spatialCMC OBJECT)

# Add sources for spatialCMC library
target_sources(spatialCMC PUBLIC ${PROTO_SRCS} ${PROTO_HDRS})
add_subdirectory(src)

# Set includes, like and compile options for spatialCMC
target_include_directories(spatialCMC PUBLIC ${INCLUDE_PATHS})
target_link_libraries(spatialCMC PUBLIC ${LINK_LIBRARIES})
target_compile_options(spatialCMC PUBLIC ${COMPILE_OPTIONS})

# Add run_mcmc executable
add_executable(run_mcmc $<TARGET_OBJECTS:spatialCMC> executables/run_mcmc.cc)
target_include_directories(run_mcmc PUBLIC ${INCLUDE_PATHS})
target_link_libraries(run_mcmc PUBLIC ${LINK_LIBRARIES})
target_compile_options(run_mcmc PUBLIC ${COMPILE_OPTIONS})

# Add run_cmc executable
add_executable(run_cmc $<TARGET_OBJECTS:spatialCMC> executables/run_cmc.cc)
target_include_directories(run_cmc PUBLIC ${INCLUDE_PATHS})
target_link_libraries(run_cmc PUBLIC ${LINK_LIBRARIES})
target_compile_options(run_cmc PUBLIC ${COMPILE_OPTIONS})