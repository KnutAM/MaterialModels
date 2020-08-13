# Setup with Windows (example): <cmake -A x64 ../src>
# Or if need of another generator e.g.: <cmake -G "Visual Studio 12 2013 Win64" ../src>
# Build from command line: <cmake --build . --config RELEASE>
cmake_minimum_required (VERSION 3.0)
message("${CMAKE_FORTRAN_COMPILER_ID}")
set(ENV{FC} "ifort")
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING
        "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif()

set(model Cruzado2015 CACHE STRING "Which base model to build")
set_property(CACHE model PROPERTY STRINGS 
	Cruzado2015 
	CP_FCC_Nort_Liniso 
	Hutchinson1976_nlkin_bcc 
	Lindfeldt2012_nlkin_bcc 
	Hutchinson1976_Cowsym_AF_bcc 
	Hutchinson1976_Norton_AF_bcc
	Hutchinson1976exp_Norton_AF0_bcc
	Hutchinson1976exp_Norton_AF_bcc
	Hutchinson1976exp_Lindfeldt_AF_bcc
    Hutchinson1976exp_Lindfeldt_AF_iso_bcc)
    
set(orientation uniform_mod CACHE STRING "Which orientation to use")
set_property(CACHE orientation PROPERTY STRINGS
    uniform_mod
    random_mod
    rand512mod
    rand16384mod)

if (NOT umat_utils_path)
    set(umat_utils_path "${CMAKE_CURRENT_SOURCE_DIR}/../../umat_utils" CACHE STRING "Folder containing umat utility routines")
endif()

enable_language (Fortran)

project(CrystalPlasticity)


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
	if (MSVC)
		add_compile_options("/fpe:0")
		add_compile_options("$<$<CONFIG:DEBUG>:/warn>")
		add_compile_options("$<$<CONFIG:DEBUG>:/check>")
	endif()
endif()


# Add files
add_library(CrystalPlasticity SHARED "taylor_model.f90")
target_sources(CrystalPlasticity PRIVATE
	"crystal_orientations_mod.f90"
    "${orientation}.f90"
	"umat.f90"
	"${umat_utils_path}/tensors_module.f90"
    "${umat_utils_path}/smsutility.f90"
    "${umat_utils_path}/SolveMatrixEquation.f90"
    "umat_module.f90"
	"${model}.f90"
    )
    
set_target_properties(CrystalPlasticity PROPERTIES PREFIX "")
set_target_properties(CrystalPlasticity PROPERTIES OUTPUT_NAME "${model}_${orientation}")
# Link to MKL library
if (MSVC)   #If Microsoft Visual compiler
    if (A64)
        target_link_libraries(CrystalPlasticity PRIVATE "mkl_intel_lp64.lib mkl_sequential.lib mkl_core.lib")
    else()
        target_link_libraries(CrystalPlasticity PRIVATE "mkl_intel_c.lib mkl_sequential.lib mkl_core.lib")
    endif()
elseif (A64)
    if (UNIX AND NOT APPLE) # I.e. linux
        set (MKLROOT $ENV{MKLROOT} CACHE STRING "Please specify your MKLROOT path")
        target_include_directories(CrystalPlasticity PRIVATE "${MKLROOT}/include")
        if (GFORT)
            target_link_libraries(CrystalPlasticity PRIVATE "-Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_gf_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl")
        elseif (IFORT)
            target_link_libraries(CrystalPlasticity PRIVATE "-Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl")
        endif ()
    endif()
else()
    message("System configuration not supported")
endif()
