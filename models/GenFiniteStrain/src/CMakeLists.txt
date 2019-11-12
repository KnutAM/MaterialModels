cmake_minimum_required (VERSION 3.0)

set(model gfs_bc CACHE STRING "Which model to build, allowed values are gfs_bc or gfs_ob")

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING
        "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif()

if (NOT umat_utils_path)
    set(umat_utils_path "${CMAKE_CURRENT_SOURCE_DIR}/../../umat_utils" CACHE STRING "Folder containing umat utility routines")
endif()
    
        
project(GFS)
enable_language (Fortran)

## To build for 64 bit and visual studio, call (from build folder): cmake -G "Visual Studio 11 2012 Win64"
## if you have this particular version of Visual studio installed.

# Get system information
if("${CMAKE_SIZEOF_VOID_P}" EQUAL "8")
    set(A64 TRUE)   # 64 bit compilation system - needs to be same as for mpsim
else()
    set(A64 FALSE)  # 32 bit compilation system - needs to be same as for mpsim
endif()
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
    set(GFORT TRUE)
    set(IFORT FALSE)
    if (A64)
        add_compile_options(-m64 -ffree-line-length-none -ffree-form)
    else ()
        add_compile_options(-ffree-line-length-none -ffree-form)
    endif()
elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
    set(GFORT FALSE)
    set(IFORT TRUE)
	add_compile_options("$<$<CONFIG:DEBUG>:/warn>")
    add_compile_options("$<$<CONFIG:DEBUG>:/check>")
endif()


# Add files
add_library(GFS SHARED "GeneralFiniteStrain.for")
target_sources(GFS PRIVATE
    "${umat_utils_path}/smsutility.f90"
    "${model}.f90"
    "${umat_utils_path}/SolveMatrixEquation.f90"
    "${umat_utils_path}/tensors_module.f90"
    "GFS_module.f90"
    )
    
set_target_properties(GFS PROPERTIES OUTPUT_NAME "${model}")
set_target_properties(GFS PROPERTIES PREFIX "")

# Link to MKL library
if (MSVC)   #If Microsoft Visual compiler
    if (A64)
        target_link_libraries(GFS PRIVATE "mkl_intel_lp64.lib mkl_sequential.lib mkl_core.lib")
    else()
        target_link_libraries(GFS PRIVATE "mkl_intel_c.lib mkl_sequential.lib mkl_core.lib")
    endif()
elseif (A64)
    if (UNIX AND NOT APPLE) # I.e. linux
        set (MKLROOT $ENV{MKLROOT} CACHE STRING "Please specify your MKLROOT path")
        target_include_directories(GFS PRIVATE "${MKLROOT}/include")
        if (GFORT)
            target_link_libraries(GFS PRIVATE "-Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_gf_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl")
        elseif (IFORT)
            target_link_libraries(GFS PRIVATE "-Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl")
        endif ()
    endif()
else()
    message("System configuration not supported")
endif()
