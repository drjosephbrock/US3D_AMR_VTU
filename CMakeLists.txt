cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

project(writeAMRVTU)

# Set C++ Standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -pg -g -D _GLIBCXX_DEBUG -D _GLIBCXX_DEBUG_PEDANTIC")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
set(VTK_DIR "/h/brock/Work/Software/VTK-9.3.0.rc1/build/")


find_package(VTK COMPONENTS 
  CommonColor
  CommonCore
  CommonDataModel
  IOXML
  InteractionStyle
  RenderingContextOpenGL2
  RenderingCore
  RenderingFreeType
  RenderingGL2PSOpenGL2
  RenderingOpenGL2
)


# Add HDF5
include_directories($ENV{HDF5_HOME}include)
link_directories($ENV{HDF5_HOME}/lib)

# add_library(my_hdf5_module MODULE read_US3D_hdf5.cpp)

# target_link_options(my_hdf5_module PRIVATE "-undefined dynamic_lookup") # requried on MacOS
# target_link_libraries(my_hdf5_module PRIVATE -undefined dynamic_lookup hdf5_cpp hdf5) # undefined dynamic_lookup required on MacOS
# target_link_libraries(my_hdf5_module PRIVATE hdf5)


if (NOT VTK_FOUND)
  message(FATAL_ERROR "writeAMRVTU: Unable to find the VTK build folder.")
endif()

# Prevent a "command line is too long" failure in Windows.
set(CMAKE_NINJA_FORCE_RESPONSE_FILE "ON" CACHE BOOL "Force Ninja to use response files.")

set(SOURCES
source/write_us3d_vtu.cpp
source/read_US3D_hdf5.cpp
  source/main.cpp
)
add_executable(writeAMRVTU MACOSX_BUNDLE ${SOURCES})

target_link_libraries(writeAMRVTU PRIVATE ${VTK_LIBRARIES} hdf5)
# vtk_module_autoinit is needed
vtk_module_autoinit(
  TARGETS writeAMRVTU
  MODULES ${VTK_LIBRARIES}
)
