cmake_minimum_required(VERSION 2.8)
project(parallel)

#find_package(glfw3 3.2 REQUIRED)
find_package(GLEW REQUIRED)
find_package(OpenGL REQUIRED)
find_package( PkgConfig )
find_package(CUDA REQUIRED)

pkg_search_module(GLFW REQUIRED glfw3)

set( CMAKE_EXPORT_COMPILE_COMMANDS ON )
SET(CUDA_PROPAGATE_HOST_FLAGS OFF)

if(APPLE)
    INCLUDE_DIRECTORIES ( /System/Library/Frameworks )
    FIND_LIBRARY(COCOA_LIBRARY Cocoa)
    MARK_AS_ADVANCED (COCOA_LIBRARY)
    SET(EXTRA_LIBS ${COCOA_LIBRARY})
else(APPLE)
    include_directories(/afs/andrew.cmu.edu/usr8/efang/public/glm/)
endif (APPLE)

include_directories(${GLFW_INCLUDE_DIRS})
include_directories(${GLEW_INCLUDE_DIRS})
include_directories(${OPENGL_INCLUDE_DIR})

set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS};-gencode arch=compute_52,code=sm_52 -gencode arch=compute_61,code=sm_61 -std=c++11)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES main.cpp parallel.h)

cuda_add_executable(parallel ${SOURCE_FILES} parallel.cu)

target_link_libraries(parallel ${CUDA_LIBRARIES} ${CUDA_LIBRARY} ${GLFW_LIBRARIES} ${GLEW_LIBRARIES} ${OPENGL_LIBRARIES})
