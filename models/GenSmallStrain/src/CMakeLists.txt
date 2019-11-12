## CMake for general small strain plasticity models for matmodfit

## To build for 64 bit and visual studio, call (from build folder): cmake -G "Visual Studio 11 2012 Win64"
## if you have this particular version of Visual studio installed.

## For building of different methods, specify the basemod and ratemod as CMake input:
## E.g. cmake -D basemod=ohnowang -D ratemod=norton ../src
## Or modify using either ccmake or cmake-gui

cmake_minimum_required (VERSION 3.0)

set(basemod chaboche CACHE STRING "Which base model to build")
set_property(CACHE basemod PROPERTY STRINGS "chaboche" "delobelle" "ohnowang")

set(ratemod norate CACHE STRING "Which overstress model to build")
set_property(CACHE ratemod PROPERTY STRINGS "norate" "norton" "cowsym" "delobelle")

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING
        "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif()

if (NOT umat_utils_path)
    set(umat_utils_path "${CMAKE_CURRENT_SOURCE_DIR}/../../umat_utils" CACHE STRING "Folder containing umat utility routines")
endif()
        
project(GeneralSmallStrainPlasticity)
enable_language (Fortran)

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
add_library(GSS SHARED "GeneralSmallStrain.f90")
target_sources(GSS PRIVATE 
    "gss_module.f90"
    "${umat_utils_path}/SolveMatrixEquation.f90"
    "${umat_utils_path}/smsutility.f90"
    )
if ("${ratemod}" STREQUAL "norate")
    target_sources(GSS PRIVATE
    "${basemod}.f90"   # Special version for the rate independent models
    "${basemod}_acegen_mod.f90"
    )
    set_target_properties(GSS PROPERTIES OUTPUT_NAME "${basemod}")
else()
    target_sources(GSS PRIVATE
    "${basemod}_rdep.f90"
    "${basemod}_${ratemod}_acegen_mod.f90"
    )
    set_target_properties(GSS PROPERTIES OUTPUT_NAME "${basemod}_${ratemod}")
endif()

set_target_properties(GSS PROPERTIES PREFIX "") # Avoid the lib prefix on linux

# Link to MKL library
if (MSVC)   #If Microsoft Visual compiler
    if (A64)
        target_link_libraries(GSS PRIVATE "mkl_intel_lp64.lib mkl_sequential.lib mkl_core.lib")
    else()
        target_link_libraries(GSS PRIVATE "mkl_intel_c.lib mkl_sequential.lib mkl_core.lib")
    endif()
elseif (A64)
    if (UNIX AND NOT APPLE) # I.e. linux
        set (MKLROOT $ENV{MKLROOT} CACHE STRING "Please specify your MKLROOT path")
        target_include_directories(GSS PRIVATE "${MKLROOT}/include")
        if (GFORT)
            target_link_libraries(GSS PRIVATE "-Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_gf_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl")
        elseif (IFORT)
            target_link_libraries(GSS PRIVATE "-Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl")
        endif ()
    endif()
else()
    message("System configuration not supported")
endif()
