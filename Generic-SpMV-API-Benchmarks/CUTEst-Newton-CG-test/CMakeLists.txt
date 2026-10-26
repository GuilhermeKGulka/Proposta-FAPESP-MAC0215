cmake_minimum_required(VERSION 3.21)
# CXX e necessario so porque o hip-config.cmake faz find_dependency(Threads)
project(cutest_newton_cg_benchmark LANGUAGES Fortran CXX)

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
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  # sem limite de 132 colunas (o gfortran trata truncamento como erro)
  add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:-ffree-line-length-none>)
endif()

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

# ---- CUTEst: decodifica o problema (sifdecoder) e linka libcutest ----------
# Requer o ambiente do CUTEst: source <raiz do repo>/env_cutest.sh
foreach(v CUTEST SIFDECODE MASTSIF MYARCH)
  if(NOT DEFINED ENV{${v}})
    message(FATAL_ERROR "Variavel de ambiente ${v} nao definida. "
                        "Rode: source <raiz do repositorio>/env_cutest.sh")
  endif()
endforeach()

# Problema a compilar (nome do .SIF em MASTSIF). Trocar -DPROBLEM re-decodifica.
set(PROBLEM "ROSENBR" CACHE STRING "Problema CUTEst")
set(SIF_FILE "$ENV{MASTSIF}/${PROBLEM}.SIF")
if(NOT EXISTS "${SIF_FILE}")
  message(FATAL_ERROR "Arquivo SIF nao encontrado: ${SIF_FILE}")
endif()

find_library(CUTEST_LIB NAMES cutest REQUIRED
  PATHS "$ENV{CUTEST}/objects/$ENV{MYARCH}/double" NO_DEFAULT_PATH)
set(CUTEST_MOD_DIR "$ENV{CUTEST}/modules/$ENV{MYARCH}/double")

# O sifdecoder gera ELFUN.f GROUP.f RANGE.f (e EXTER.f em alguns problemas) e OUTSDIF.d
# no diretorio de build; so roda de novo se o problema mudou.
if(NOT "${CUTEST_DECODED}" STREQUAL "${PROBLEM}" OR NOT EXISTS "${CMAKE_BINARY_DIR}/ELFUN.f")
  file(REMOVE "${CMAKE_BINARY_DIR}/ELFUN.f" "${CMAKE_BINARY_DIR}/GROUP.f" "${CMAKE_BINARY_DIR}/RANGE.f"
              "${CMAKE_BINARY_DIR}/EXTER.f" "${CMAKE_BINARY_DIR}/OUTSDIF.d" "${CMAKE_BINARY_DIR}/AUTOMAT.d")
  message(STATUS "sifdecoder ${PROBLEM}")
  execute_process(COMMAND "$ENV{SIFDECODE}/bin/sifdecoder" "${SIF_FILE}"
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
    RESULT_VARIABLE _rc OUTPUT_VARIABLE _out ERROR_VARIABLE _err)
  if(NOT _rc EQUAL 0 OR NOT EXISTS "${CMAKE_BINARY_DIR}/ELFUN.f")
    set(CUTEST_DECODED "" CACHE INTERNAL "")
    message(FATAL_ERROR "sifdecoder falhou para ${PROBLEM} (codigo ${_rc}):\n${_out}\n${_err}")
  endif()
  set(CUTEST_DECODED "${PROBLEM}" CACHE INTERNAL "problema decodificado")
endif()

set(DECODED_SRCS "${CMAKE_BINARY_DIR}/ELFUN.f" "${CMAKE_BINARY_DIR}/GROUP.f" "${CMAKE_BINARY_DIR}/RANGE.f")
if(EXISTS "${CMAKE_BINARY_DIR}/EXTER.f")
  list(APPEND DECODED_SRCS "${CMAKE_BINARY_DIR}/EXTER.f")
endif()

# ---- alvos -----------------------------------------------------------------
# Parte FIXA (compilada uma vez, independe do problema)
add_library(newtoncore STATIC
  hipsparse_generic.f90 rocm_context.f90 benchmark_cg_mod.f90 cutest_iface.f90 newton_cg_mod.f90)
target_link_libraries(newtoncore PUBLIC hipfort_backend)
target_include_directories(newtoncore PUBLIC "${CUTEST_MOD_DIR}")

# Parte POR PROBLEMA: o programa + as rotinas decodificadas
add_executable(benchmark_newton_cg benchmark_newton_cg.f90 ${DECODED_SRCS})
target_link_libraries(benchmark_newton_cg PRIVATE newtoncore kernels ${CUTEST_LIB} hipfort_backend)
