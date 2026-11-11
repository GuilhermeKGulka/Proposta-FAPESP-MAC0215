cmake_minimum_required(VERSION 3.21)
project(gencan_cg_standalone LANGUAGES Fortran CXX)

# Backend HIP
set(HIP_PLATFORM "amd" CACHE STRING "Backend HIP: amd | nvidia")
if(HIP_PLATFORM STREQUAL "amd")
  set(HIPFORT_BACKEND amdgcn)
elseif(HIP_PLATFORM STREQUAL "nvidia")
  set(HIPFORT_BACKEND nvptx)
else()
  message(FATAL_ERROR "HIP_PLATFORM deve ser 'amd' ou 'nvidia'")
endif()
set(ENV{HIP_PLATFORM} "${HIP_PLATFORM}")

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:-ffree-line-length-none>)
endif()

# ---- hipfort ----
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

find_path(HIPFORT_MOD_DIR hipfort_hipsparse.mod HINTS ENV HIPFORT_ROOT ${HIPFORT_ROOT} PATH_SUFFIXES ${_mod_suffixes})
find_library(HIPFORT_LIB NAMES hipfort-${HIPFORT_BACKEND} HINTS ENV HIPFORT_ROOT ${HIPFORT_ROOT} PATH_SUFFIXES ${_lib_suffixes})
if(NOT HIPFORT_MOD_DIR OR NOT HIPFORT_LIB)
  message(FATAL_ERROR "hipfort nao encontrado.")
endif()

add_library(hipfort_backend INTERFACE)
target_include_directories(hipfort_backend INTERFACE ${HIPFORT_MOD_DIR})
target_link_libraries(hipfort_backend INTERFACE ${HIPFORT_LIB})

if(HIP_PLATFORM STREQUAL "amd")
  find_package(hip REQUIRED)
  find_package(hipsparse REQUIRED)
  find_package(hipblas REQUIRED)
  target_link_libraries(hipfort_backend INTERFACE hip::host roc::hipsparse roc::hipblas)
else()
  find_package(CUDAToolkit REQUIRED)
  foreach(lib hipsparse hipblas)
    find_package(${lib} QUIET)
    if(TARGET roc::${lib})
      target_link_libraries(hipfort_backend INTERFACE roc::${lib})
    else()
      find_library(${lib}_LIB NAMES ${lib} REQUIRED HINTS ENV HIPFORT_ROOT ${HIPFORT_ROOT} PATH_SUFFIXES lib lib64)
      target_link_libraries(hipfort_backend INTERFACE ${${lib}_LIB})
    endif()
  endforeach()
  target_link_libraries(hipfort_backend INTERFACE CUDA::cudart CUDA::cusparse CUDA::cublas)
endif()

# ---- Kernels C++ ----
find_program(HIPCC NAMES hipcc REQUIRED HINTS /opt/rocm PATH_SUFFIXES bin)
set(KERNEL_OBJ ${CMAKE_BINARY_DIR}/scalar_kernels.o)
add_custom_command(OUTPUT ${KERNEL_OBJ}
  COMMAND ${CMAKE_COMMAND} -E env HIP_PLATFORM=${HIP_PLATFORM}
          ${HIPCC} -O3 -c ${CMAKE_CURRENT_SOURCE_DIR}/sources/scalar_kernels.cpp -o ${KERNEL_OBJ}
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/sources/scalar_kernels.cpp)
set_source_files_properties(${KERNEL_OBJ} PROPERTIES EXTERNAL_OBJECT TRUE GENERATED TRUE)
add_library(kernels STATIC ${KERNEL_OBJ})
set_target_properties(kernels PROPERTIES LINKER_LANGUAGE Fortran)
target_link_libraries(kernels PUBLIC hipfort_backend stdc++)

# ---- Bibliotecas Fortran ----
# 1. Módulo base (lssit)
add_library(lssit_lib STATIC sources/lssit.f90)
target_link_libraries(lssit_lib PUBLIC hipfort_backend)

# 2. O Gencan modificado para usar GPU
add_library(gencan_gpu STATIC 
  sources/hipsparse_generic.f90 
  sources/gencan_gpu_context.f90 
  sources/cgm_gpu.f90 
  sources/gencan.f90
)
target_link_libraries(gencan_gpu PUBLIC lssit_lib)

# ---- Executáveis ----
# 1. Toyprob (Exemplo de Rosenbrock)
add_executable(toyprob sources/toyprob.f90)
target_link_libraries(toyprob PRIVATE gencan_gpu kernels)

# 2. CUTEst
foreach(v CUTEST SIFDECODE MASTSIF MYARCH)
  if(NOT DEFINED ENV{${v}})
    message(FATAL_ERROR "Variavel de ambiente ${v} nao definida. Rode: source env_cutest.sh")
  endif()
endforeach()

set(PROBLEM "ROSENBR" CACHE STRING "Problema CUTEst")
set(SIF_FILE "$ENV{MASTSIF}/${PROBLEM}.SIF")
find_library(CUTEST_LIB NAMES cutest REQUIRED PATHS "$ENV{CUTEST}/objects/$ENV{MYARCH}/double" NO_DEFAULT_PATH)
set(CUTEST_MOD_DIR "$ENV{CUTEST}/modules/$ENV{MYARCH}/double")

if(NOT "${CUTEST_DECODED}" STREQUAL "${PROBLEM}" OR NOT EXISTS "${CMAKE_BINARY_DIR}/ELFUN.f")
  file(REMOVE "${CMAKE_BINARY_DIR}/ELFUN.f" "${CMAKE_BINARY_DIR}/GROUP.f" "${CMAKE_BINARY_DIR}/RANGE.f" "${CMAKE_BINARY_DIR}/EXTER.f" "${CMAKE_BINARY_DIR}/OUTSDIF.d" "${CMAKE_BINARY_DIR}/AUTOMAT.d")
  execute_process(COMMAND "$ENV{SIFDECODE}/bin/sifdecoder" "${SIF_FILE}" WORKING_DIRECTORY "${CMAKE_BINARY_DIR}")
  if(NOT EXISTS "${CMAKE_BINARY_DIR}/ELFUN.f")
    message(FATAL_ERROR "sifdecoder falhou para ${PROBLEM}")
  endif()
  set(CUTEST_DECODED "${PROBLEM}" CACHE INTERNAL "problema decodificado")
endif()

set(DECODED_SRCS "${CMAKE_BINARY_DIR}/ELFUN.f" "${CMAKE_BINARY_DIR}/GROUP.f" "${CMAKE_BINARY_DIR}/RANGE.f")
if(EXISTS "${CMAKE_BINARY_DIR}/EXTER.f")
  list(APPEND DECODED_SRCS "${CMAKE_BINARY_DIR}/EXTER.f")
endif()

add_executable(cutestma sources/cutestma.f90 sources/cutest_iface.f90 ${DECODED_SRCS})
target_include_directories(cutestma PRIVATE "${CUTEST_MOD_DIR}")
target_link_libraries(cutestma PRIVATE gencan_gpu kernels ${CUTEST_LIB} hipfort_backend)