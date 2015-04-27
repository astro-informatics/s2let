#------------------

cmake_minimum_required(VERSION 2.6)

INCLUDE (CheckLibraryExists)

# http://www.cmake.org/Wiki/CMake:How_To_Write_Platform_Checks
# http://www.cmake.org/cmake/help/cmake2.6docs.html

project(s2let)
enable_language(C)
enable_language(Fortran)

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_LIST_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${CMAKE_CURRENT_LIST_DIR}/lib)

#------------------
# Dependencies
#------------------

# FFTW
link_directories($ENV{FFTW}/lib)

# SSHT
include_directories($ENV{SSHT}/include/c)
link_directories($ENV{SSHT}/lib/c)


# so3
#include_directories($ENV{so3}/include/c)
#link_directories($ENV{so3}/lib/c)

# CFITSIO 
find_path(CFITSIO_INCLUDE_DIR fitsio.h
  PATH_SUFFIXES libcfitsio3 libcfitsio0 cfitsio
  PATHS
  $ENV{CFITSIO}/include  
)
find_library(CFITSIO_LIB NAMES cfitsio
  PATHS
  $ENV{CFITSIO}/lib  
)
if(CFITSIO_INCLUDE_DIR AND CFITSIO_LIB)
  set(CFITSIO_FOUND TRUE)
else (CFITSIO_INCLUDE_DIR AND CFITSIO_LIBRARIES)
  set(CFITSIO_FOUND FALSE)
endif(CFITSIO_INCLUDE_DIR AND CFITSIO_LIB)
# Link if found
if (CFITSIO_FOUND)
   include_directories($ENV{CFITSIO}/include)
   link_directories($ENV{CFITSIO}/lib)
endif (CFITSIO_FOUND)

# HEALPIX 
find_library(HEALPIX_LIB NAMES healpix
  PATHS
  $ENV{HEALPIX}/lib  
)
if(HEALPIX_LIB)
  set(HEALPIX_FOUND TRUE)
else (HEALPIX_LIB)
  set(HEALPIX_FOUND FALSE)
endif(HEALPIX_LIB)
# Link if found
if (HEALPIX_FOUND)
   include_directories($ENV{HEALPIX}/include)
   link_directories($ENV{HEALPIX}/lib)
endif (HEALPIX_FOUND)

#------------------
# Compiler flags : custom if needed
#------------------

# get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
# Wall -O3 -g -fopenmp")
# if (Fortran_COMPILER_NAME STREQUAL "gfortran")
#  #elseif (Fortran_COMPILER_NAME STREQUAL "ifort")
#  set(CMAKE_Fortran_FLAGS "-TODO")
#elseif (Fortran_COMPILER_NAME STREQUAL "g77") 
#  set(CMAKE_Fortran_FLAGS "-TODO")
#else (Fortran_COMPILER_NAME STREQUAL "gfortran")
#  set(CMAKE_Fortran_FLAGS "-TODO")
#endif (Fortran_COMPILER_NAME STREQUAL "gfortran")

#------------------
# Standard configuration, main library
#------------------

set(LIBS_C "-lfftw3 -lssht")
set(LIBTAG_C "-L${LIBRARY_OUTPUT_PATH} -ls2let -lm -lc")
message("Library output :" ${LIBRARY_OUTPUT_PATH})

# Standard library
set(C_SRC "src/main/c/s2let_axisym_lm.c" "src/main/c/s2let_axisym_mw.c" "src/main/c/s2let_idl_mw.c" "src/main/c/s2let_math.c" "src/main/c/s2let_mw.c" "src/main/c/s2let_tiling.c")

# Test executables
set(C_TESTEXE "src/test/c/s2let_test.c")  

# Main executables
set(C_EXE "src/main/c/s2let_denoising_demo.c")  

#------------------
# CFITS and HEALPIX options
#------------------

# Additional files if CFITSIO is present
if (CFITSIO_FOUND)
   message("CFITSIO library was found")
   set(C_SRC ${C_SRC} "src/main/c/s2let_fits.c")
   set(LIBS_C ${LIBS_C} "-lcfitsio")
   set(C_EXE ${C_EXE} "src/main/c/s2let_axisym_mw_analysis_real.c")  
   set(C_EXE ${C_EXE} "src/main/c/s2let_axisym_mw_synthesis_real.c") 
endif (CFITSIO_FOUND)

# Additional files if HEALPIX is present
if (HEALPIX_FOUND)
   message("HEALPIX library was found")
   set(CMAKE_Fortran_COMPILER_INIT gfortran)
   set(CMAKE_Fortran_FLAGS "-O3 -ffree-form -lgfortran -DGFORTRAN -fno-second-underscore -fopenmp")
   message("Setting F90 compilation flags as: ${CMAKE_Fortran_FLAGS}")
   set(C_SRC ${C_SRC} "src/main/f90/s2let_hpx.f90" "src/main/c/s2let_hpx.c" "src/main/c/s2let_axisym_hpx.c" "src/main/c/s2let_idl_hpx.c")
   set(LIBS_C ${LIBS_C} "-lgfortran -fopenmp -lhealpix")
   set(C_EXE ${C_EXE} "src/main/c/s2let_hpx_demo.c")  
   set(C_EXE ${C_EXE} "src/main/c/s2let_axisym_hpx_analysis_real.c") 
   set(C_EXE ${C_EXE} "src/main/c/s2let_axisym_hpx_synthesis_real.c") 
   set(C_TESTEXE ${C_TESTEXE} "src/test/c/s2let_hpx_test.c")  
endif (HEALPIX_FOUND)

# Create library
include_directories(${CMAKE_CURRENT_LIST_DIR}/include)
add_library(s2let ${C_SRC})

#------------------
# Config test executables
#------------------

set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${EXECUTABLE_OUTPUT_PATH})
foreach (p ${C_TESTEXE})
  get_filename_component(q ${p} NAME_WE)
  message( "Executable : ${q} " )
  add_executable (${q} ${p})
  add_dependencies(${q} s2let)
  target_link_libraries(${q} ${LIBTAG_C} ${LIBS_C})
  add_test (${q} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${q})
endforeach (p)

enable_testing()

#------------------
# Config main executables
#------------------
foreach (p ${C_EXE})
  get_filename_component(q ${p} NAME_WE)
  message( "Executable : ${q} " )
  add_executable (${q} ${p})
  add_dependencies(${q} s2let)
  target_link_libraries(${q} ${LIBTAG_C} ${LIBS_C})
endforeach (p)

#------------------
