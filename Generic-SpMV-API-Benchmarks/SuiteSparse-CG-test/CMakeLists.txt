cmake_minimum_required(VERSION 3.21)
# CXX e necessario so porque o hip-config.cmake faz find_dependency(Threads)
project(suitesparse_cg_benchmark LANGUAGES Fortran CXX)

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

# ---- runtime + sparse + blas do backend ------------------------------------
if(HIP_PLATFORM STREQUAL "amd")
  find_package(hip REQUIRED)
  find_package(hipsparse REQUIRED)
  find_package(hipblas REQUIRED)
  target_link_libraries(hipfort_backend INTERFACE hip::host roc::hipsparse roc::hipblas)
else()
  find_package(CUDAToolkit REQUIRED)  # use -DCUDAToolkit_ROOT=/usr/local/cuda se nao achar
  # Os simbolos hipsparse*/hipblas* vem de libs compiladas para CUDA (wrappers de cuSPARSE/cuBLAS)
  foreach(lib hipsparse hipblas)
    find_package(${lib} QUIET)
    if(TARGET roc::${lib})
      target_link_libraries(hipfort_backend INTERFACE roc::${lib})
    else()
      find_library(${lib}_LIB NAMES ${lib} REQUIRED
        HINTS ENV HIPFORT_ROOT ${HIPFORT_ROOT} PATH_SUFFIXES lib lib64)
      target_link_libraries(hipfort_backend INTERFACE ${${lib}_LIB})
    endif()
  endforeach()
  target_link_libraries(hipfort_backend INTERFACE CUDA::cudart CUDA::cusparse CUDA::cublas)
endif()

# ---- kernels escalares (HIP/C++), compilados com o hipcc do backend --------
find_program(HIPCC NAMES hipcc REQUIRED HINTS /opt/rocm PATH_SUFFIXES bin)
set(KERNEL_FLAGS "" CACHE STRING "flags extras do hipcc (ex.: --offload-arch=gfx1201)")
separate_arguments(_kernel_flags NATIVE_COMMAND "${KERNEL_FLAGS}")
set(KERNEL_OBJ ${CMAKE_BINARY_DIR}/scalar_kernels.o)
add_custom_command(OUTPUT ${KERNEL_OBJ}
  COMMAND ${CMAKE_COMMAND} -E env HIP_PLATFORM=${HIP_PLATFORM}
          ${HIPCC} -O3 ${_kernel_flags} -c ${CMAKE_CURRENT_SOURCE_DIR}/scalar_kernels.cpp -o ${KERNEL_OBJ}
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/scalar_kernels.cpp
  VERBATIM)
set_source_files_properties(${KERNEL_OBJ} PROPERTIES EXTERNAL_OBJECT TRUE GENERATED TRUE)
add_library(kernels STATIC ${KERNEL_OBJ})
set_target_properties(kernels PROPERTIES LINKER_LANGUAGE Fortran)
target_link_libraries(kernels PUBLIC hipfort_backend stdc++)

# ---- alvos -----------------------------------------------------------------
add_library(mmio STATIC mmio.f90)

add_library(hsgen STATIC hipsparse_generic.f90)
target_link_libraries(hsgen PUBLIC hipfort_backend)  # usa os enums do hipfort

add_library(cgcore STATIC rocm_context.f90 benchmark_cg_mod.f90)
target_link_libraries(cgcore PUBLIC hsgen hipfort_backend)

add_executable(benchmark_cg benchmark_cg.f90)
target_link_libraries(benchmark_cg PRIVATE cgcore kernels mmio hsgen hipfort_backend)