cmake_minimum_required(VERSION 3.0.0)

project(gamil LANGUAGES Fortran)

if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none -ffpe-trap=invalid,zero,overflow")
  if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O0 -fcheck=all -g -fbacktrace")
  else ()
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Ofast")
  endif ()
  # FIXME: OpenMP causes sigmentation faults.
  # set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fopenmp")
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -no-wrap-margin -O0 -check bounds -g -traceback -mkl=sequential")
  else ()
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -no-wrap-margin -O3 -align array64byte -mkl=sequential")
  endif ()
endif ()

find_package(MPI)

add_subdirectory(lib/container)
add_subdirectory(lib/string)
add_subdirectory(lib/flogger)
add_subdirectory(lib/fiona)

set(sources
  src/kinds_mod.F90
  src/mesh_const_mod.F90
  src/mesh_math_mod.F90
  src/gauss_quad_mod.F90
  src/latlon_mesh_mod.F90
  src/gamil_params_mod.F90
  src/gamil_mod.F90
)

add_library(gamil ${sources})
target_link_libraries(gamil fiona)

add_executable(gamil_driver.exe src/gamil_driver.F90)
target_link_libraries(gamil_driver.exe gamil)
