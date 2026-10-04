cmake_minimum_required(VERSION 3.21)
# CXX e necessario so porque o hip-config.cmake faz find_dependency(Threads)
project(suitesparse_spmv_benchmark LANGUAGES Fortran CXX)

# Backend HIP: "amd" (ROCm) ou "nvidia" (CUDA). Normalmente vem do preset.
set(HIP_PLATFORM "amd" CACHE STRING "Backend HIP: amd | nvidia")
set_property(CACHE HIP_PLATFORM PROPERTY STRINGS amd nvidia)
if(HIP_PLATFORM STREQUAL "amd")
  set(HIPFORT_BACKEND amdgcn)
elseif(HIP_PLATFORM STREQUAL "nvidia")
  set(HIPFORT_BACKEND nvptx)
else()
  message(FATAL_ERROR "HIP_PLATFORM deve ser 'amd' ou 'nvidia' (recebido: '${HIP_PLATFORM}')")
endif()
set(ENV{HIP_PLATFORM} "${HIP_PLATFORM}")

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release)
endif()
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)

# ---- hipfort: .mod e libhipfort-<backend>.a --------------------------------
# Procurados diretamente (os targets exportados do hipfort variam entre
# instalacoes). Layouts aceitos: include/fortran/<compilador>/hipfort/<backend>
# e include/hipfort/<backend>. Para sobrescrever: -DHIPFORT_MOD_DIR / -DHIPFORT_LIB,
# ou aponte o prefixo com HIPFORT_ROOT / CMAKE_PREFIX_PATH.
get_filename_component(_fc_name "${CMAKE_Fortran_COMPILER}" NAME)
set(_fc_dirs "${_fc_name}")
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  list(APPEND _fc_dirs gfortran f95)
endif()
list(REMOVE_DUPLICATES _fc_dirs)

set(_mod_suffixes)
set(_lib_suffixes)
foreach(d IN LISTS _fc_dirs)
  list(APPEND _mod_suffixes include/fortran/${d}/hipfort/${HIPFORT_BACKEND})
  list(APPEND _lib_suffixes lib/fortran/${d})
endforeach()
list(APPEND _mod_suffixes include/hipfort/${HIPFORT_BACKEND})
list(APPEND _lib_suffixes lib lib64)

find_path(HIPFORT_MOD_DIR hipfort_hipsparse.mod
  HINTS ENV HIPFORT_ROOT ${HIPFORT_ROOT}
  PATH_SUFFIXES ${_mod_suffixes})
find_library(HIPFORT_LIB NAMES hipfort-${HIPFORT_BACKEND}
  HINTS ENV HIPFORT_ROOT ${HIPFORT_ROOT}
  PATH_SUFFIXES ${_lib_suffixes})
if(NOT HIPFORT_MOD_DIR OR NOT HIPFORT_LIB)
  message(FATAL_ERROR
    "hipfort (backend ${HIPFORT_BACKEND}) nao encontrado.\n"
    "  HIPFORT_MOD_DIR=${HIPFORT_MOD_DIR}\n  HIPFORT_LIB=${HIPFORT_LIB}\n"
    "Use -DHIPFORT_ROOT=<prefixo> ou -DHIPFORT_MOD_DIR=... -DHIPFORT_LIB=...")
endif()

add_library(hipfort_backend INTERFACE)
target_include_directories(hipfort_backend INTERFACE ${HIPFORT_MOD_DIR})
target_link_libraries(hipfort_backend INTERFACE ${HIPFORT_LIB})

# ---- runtime + sparse do backend -------------------------------------------
if(HIP_PLATFORM STREQUAL "amd")
  find_package(hip REQUIRED)
  find_package(hipsparse REQUIRED)
  target_link_libraries(hipfort_backend INTERFACE hip::host roc::hipsparse)
else()
  find_package(CUDAToolkit REQUIRED)  # use -DCUDAToolkit_ROOT=/usr/local/cuda se nao achar
  # Os simbolos hipsparse* vem de uma libhipsparse compilada para CUDA (wrapper do cuSPARSE)
  find_package(hipsparse QUIET)
  if(TARGET roc::hipsparse)
    target_link_libraries(hipfort_backend INTERFACE roc::hipsparse)
  else()
    find_library(HIPSPARSE_LIB hipsparse REQUIRED
      HINTS ENV HIPFORT_ROOT ${HIPFORT_ROOT} PATH_SUFFIXES lib lib64)
    target_link_libraries(hipfort_backend INTERFACE ${HIPSPARSE_LIB})
  endif()
  target_link_libraries(hipfort_backend INTERFACE CUDA::cudart CUDA::cusparse)
endif()

# ---- alvos -----------------------------------------------------------------
add_library(mmio STATIC mmio.f90)
add_library(hsgen STATIC hipsparse_generic.f90)
target_link_libraries(hsgen PUBLIC hipfort_backend)  # usa os enums do hipfort

add_executable(benchmark_suitesparse benchmark_suitesparse.f90)
target_link_libraries(benchmark_suitesparse PRIVATE mmio hsgen hipfort_backend)

# Teste minimo da API generica (3x4, bases 0 e 1)
add_executable(test_generic test_generic.f90)
target_link_libraries(test_generic PRIVATE hsgen hipfort_backend)