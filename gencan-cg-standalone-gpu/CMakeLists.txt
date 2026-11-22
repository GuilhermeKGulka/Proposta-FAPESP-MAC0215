cmake_minimum_required(VERSION 3.21)
# So Fortran e obrigatorio. CXX e habilitado apenas com GPU (o hip-config.cmake
# faz find_dependency(Threads)).
project(gencan LANGUAGES Fortran)

# ---- opcoes -----------------------------------------------------------------
# GENCAN_GPU: AUTO liga a GPU se o backend for encontrado (senao, so CPU, com aviso);
#             ON exige o backend (erro se faltar algo); OFF nem procura.
set(GENCAN_GPU "AUTO" CACHE STRING "Suporte a GPU: AUTO | ON | OFF")
set_property(CACHE GENCAN_GPU PROPERTY STRINGS AUTO ON OFF)
if(NOT GENCAN_GPU MATCHES "^(AUTO|ON|OFF)$")
  message(FATAL_ERROR "GENCAN_GPU deve ser AUTO, ON ou OFF (recebido: '${GENCAN_GPU}')")
endif()

# Backend HIP: "amd" (ROCm) ou "nvidia" (CUDA). Normalmente vem do preset.
set(HIP_PLATFORM "amd" CACHE STRING "Backend HIP: amd | nvidia")
set_property(CACHE HIP_PLATFORM PROPERTY STRINGS amd nvidia)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release)
endif()
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  # sem limite de 132 colunas (o gfortran trata truncamento como erro)
  add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:-ffree-line-length-none>)
endif()
# modulos da parte GPU (hipfort, hess_csr, cgm_gpu, gencan_gpu); os do Gencan
# CPU ficam em mod/cpu (ver gencan_cpu abaixo)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod/gpu)

# ---- Gencan CPU (sempre) ----------------------------------------------------
# Mesmos fontes da referencia gencan-cg-standalone/. Pre-processados (-cpp) para
# receber os blocos #ifdef GENCAN_GPU da fase 6; sem a macro, o codigo e o original.
set(GENCAN_SRCS sources/lssit.f90 sources/cgm.f90 sources/gencan.f90)
set_source_files_properties(${GENCAN_SRCS} PROPERTIES Fortran_PREPROCESS ON)
add_library(gencan_cpu STATIC ${GENCAN_SRCS})
set_target_properties(gencan_cpu PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod/cpu)
target_include_directories(gencan_cpu PUBLIC ${CMAKE_BINARY_DIR}/mod/cpu)

add_executable(toyprob_cpu sources/toyprob.f90)
target_link_libraries(toyprob_cpu PRIVATE gencan_cpu)

# ---- deteccao do backend GPU ------------------------------------------------
# A logica de busca e a dos Generic-SpMV-API-Benchmarks (testada em ROCm e CUDA);
# aqui ela so deixa de ser fatal no modo AUTO. O que faltar vai para _gpu_missing.
set(GPU_ENABLED FALSE)
set(_gpu_missing)
if(NOT GENCAN_GPU STREQUAL "OFF")
  if(HIP_PLATFORM STREQUAL "amd")
    set(HIPFORT_BACKEND amdgcn)
  elseif(HIP_PLATFORM STREQUAL "nvidia")
    set(HIPFORT_BACKEND nvptx)
  else()
    message(FATAL_ERROR "HIP_PLATFORM deve ser 'amd' ou 'nvidia' (recebido: '${HIP_PLATFORM}')")
  endif()
  set(ENV{HIP_PLATFORM} "${HIP_PLATFORM}")

  include(CheckLanguage)
  check_language(CXX)
  if(CMAKE_CXX_COMPILER)
    enable_language(CXX)
  else()
    list(APPEND _gpu_missing "compilador C++")
  endif()

  # hipfort: .mod e libhipfort-<backend>.a, procurados diretamente (os targets
  # exportados do hipfort variam entre instalacoes). Layouts aceitos:
  # include/fortran/<compilador>/hipfort/<backend> e include/hipfort/<backend>.
  # Para sobrescrever: -DHIPFORT_MOD_DIR / -DHIPFORT_LIB, ou HIPFORT_ROOT / CMAKE_PREFIX_PATH.
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
    list(APPEND _gpu_missing "hipfort (backend ${HIPFORT_BACKEND}, use -DHIPFORT_ROOT=... ou -DHIPFORT_MOD_DIR/-DHIPFORT_LIB)")
  endif()

  # runtime + sparse + blas do backend
  set(_backend_libs)
  if(HIP_PLATFORM STREQUAL "amd")
    if(CMAKE_CXX_COMPILER)
      foreach(pkg hip hipsparse hipblas)
        find_package(${pkg} QUIET)
        if(NOT ${pkg}_FOUND)
          list(APPEND _gpu_missing "${pkg} (ROCm, use -DCMAKE_PREFIX_PATH=/opt/rocm)")
        endif()
      endforeach()
      set(_backend_libs hip::host roc::hipsparse roc::hipblas)
    endif()
  else()
    find_package(CUDAToolkit QUIET)   # use -DCUDAToolkit_ROOT=/usr/local/cuda se nao achar
    if(NOT CUDAToolkit_FOUND)
      list(APPEND _gpu_missing "CUDA Toolkit (use -DCUDAToolkit_ROOT=...)")
    endif()
    # Os simbolos hipsparse*/hipblas* vem de libs compiladas para CUDA (wrappers de cuSPARSE/cuBLAS)
    foreach(lib hipsparse hipblas)
      if(CMAKE_CXX_COMPILER)
        find_package(${lib} QUIET)
      endif()
      if(TARGET roc::${lib})
        list(APPEND _backend_libs roc::${lib})
      else()
        find_library(${lib}_LIB NAMES ${lib}
          HINTS ENV HIPFORT_ROOT ${HIPFORT_ROOT} PATH_SUFFIXES lib lib64)
        if(${lib}_LIB)
          list(APPEND _backend_libs ${${lib}_LIB})
        else()
          list(APPEND _gpu_missing "${lib} compilado para CUDA")
        endif()
      endif()
    endforeach()
    list(APPEND _backend_libs CUDA::cudart CUDA::cusparse CUDA::cublas)
  endif()

  find_program(HIPCC NAMES hipcc HINTS /opt/rocm PATH_SUFFIXES bin)
  if(NOT HIPCC)
    list(APPEND _gpu_missing "hipcc")
  endif()

  if(_gpu_missing)
    list(JOIN _gpu_missing "\n  - " _msg)
    if(GENCAN_GPU STREQUAL "ON")
      message(FATAL_ERROR "GENCAN_GPU=ON, mas falta (backend ${HIP_PLATFORM}):\n  - ${_msg}")
    endif()
    message(WARNING "GPU desligada (GENCAN_GPU=AUTO, backend ${HIP_PLATFORM}); falta:\n  - ${_msg}\n"
                    "Compilando so a versao CPU.")
  else()
    set(GPU_ENABLED TRUE)
  endif()
endif()

# ---- parte GPU ----------------------------------------------------------------
if(GPU_ENABLED)
  add_library(hipfort_backend INTERFACE)
  target_include_directories(hipfort_backend INTERFACE ${HIPFORT_MOD_DIR})
  target_link_libraries(hipfort_backend INTERFACE ${HIPFORT_LIB} ${_backend_libs})

  # kernels HIP (so os do cgm), compilados com o hipcc do backend
  set(KERNEL_FLAGS "" CACHE STRING "flags extras do hipcc (ex.: --offload-arch=gfx1201)")
  separate_arguments(_kernel_flags NATIVE_COMMAND "${KERNEL_FLAGS}")
  set(KERNEL_OBJ ${CMAKE_BINARY_DIR}/cgm_kernels.o)
  add_custom_command(OUTPUT ${KERNEL_OBJ}
    COMMAND ${CMAKE_COMMAND} -E env HIP_PLATFORM=${HIP_PLATFORM}
            ${HIPCC} -O3 ${_kernel_flags} -c ${CMAKE_CURRENT_SOURCE_DIR}/sources/cgm_kernels.cpp -o ${KERNEL_OBJ}
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/sources/cgm_kernels.cpp VERBATIM)
  set_source_files_properties(${KERNEL_OBJ} PROPERTIES EXTERNAL_OBJECT TRUE GENERATED TRUE)
  add_library(kernels STATIC ${KERNEL_OBJ})
  set_target_properties(kernels PROPERTIES LINKER_LANGUAGE Fortran)
  target_link_libraries(kernels PUBLIC hipfort_backend stdc++)

  add_library(hsgen STATIC sources/hipsparse_generic.f90)
  target_link_libraries(hsgen PUBLIC hipfort_backend)
  target_include_directories(hsgen PUBLIC ${CMAKE_BINARY_DIR}/mod/gpu)

  # Gencan GPU da fase 6: os MESMOS fontes de gencan_cpu, com -DGENCAN_GPU, mais a
  # Hessiana na GPU (lssit_gpu.f90). Modulos em mod/gpu6 (lssit.mod etc. tem o
  # mesmo nome dos da CPU). Nome provisorio: vira gencan_gpu no passo 5.
  add_library(gencan_gpu6 STATIC ${GENCAN_SRCS} sources/lssit_gpu.f90)
  target_compile_definitions(gencan_gpu6 PRIVATE GENCAN_GPU)
  set_target_properties(gencan_gpu6 PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod/gpu6)
  target_include_directories(gencan_gpu6 PUBLIC ${CMAKE_BINARY_DIR}/mod/gpu6)
  target_link_libraries(gencan_gpu6 PUBLIC hsgen kernels hipfort_backend)   # kernels: scatter

  # caminho da fase 5 (sai nos passos 4-5): hess_csr nao usa lssit; cgm_gpu usa o da CPU
  add_library(hesscsr STATIC sources/hess_csr.f90)
  target_link_libraries(hesscsr PUBLIC hsgen hipfort_backend)

  add_library(gpucore STATIC sources/cgm_gpu.f90)
  target_link_libraries(gpucore PUBLIC hesscsr hsgen hipfort_backend gencan_cpu)

  add_library(gencan_gpu STATIC sources/gencan_gpu.f90)
  target_link_libraries(gencan_gpu PUBLIC gpucore)

  # test_hess_csr: lssit GPU (usegpu) contra lssit CPU (usegpu = .false.), ambos de
  # gencan_gpu6; hesscsr so para medir contra a fase 5
  add_executable(test_hess_csr tests/test_hess_csr.f90)
  target_link_libraries(test_hess_csr PRIVATE gencan_gpu6 hesscsr kernels hsgen hipfort_backend)

  add_executable(test_cgm_gpu tests/test_cgm_gpu.f90)
  target_link_libraries(test_cgm_gpu PRIVATE gpucore gencan_cpu kernels hsgen hipfort_backend)
endif()

# ---- CUTEst (opcional): decodifica o problema (sifdecoder) e linka libcutest ----
# Requer o ambiente do CUTEst: source <raiz do repo>/env_cutest.sh. Sem ele, os
# alvos cutestma_* nao sao criados (o resto compila normalmente).
set(CUTEST_ENABLED TRUE)
set(_cutest_missing)
foreach(v CUTEST SIFDECODE MASTSIF MYARCH)
  if(NOT DEFINED ENV{${v}})
    set(CUTEST_ENABLED FALSE)
    list(APPEND _cutest_missing ${v})
  endif()
endforeach()

if(CUTEST_ENABLED)
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

  add_executable(cutestma_cpu sources/cutestma.f90 ${DECODED_SRCS})
  target_include_directories(cutestma_cpu PRIVATE "${CUTEST_MOD_DIR}")
  target_link_libraries(cutestma_cpu PRIVATE gencan_cpu ${CUTEST_LIB})

  if(GPU_ENABLED)
    add_executable(cutestma_gpu sources/cutestma_gpu.f90 ${DECODED_SRCS})
    target_include_directories(cutestma_gpu PRIVATE "${CUTEST_MOD_DIR}")
    target_link_libraries(cutestma_gpu PRIVATE gencan_gpu gpucore kernels hsgen ${CUTEST_LIB} hipfort_backend)
  endif()
endif()

# ---- resumo -----------------------------------------------------------------
if(GPU_ENABLED)
  message(STATUS "Gencan: GPU ligada (backend ${HIP_PLATFORM})")
else()
  message(STATUS "Gencan: GPU desligada (GENCAN_GPU=${GENCAN_GPU}); so CPU")
endif()
if(CUTEST_ENABLED)
  message(STATUS "Gencan: CUTEst ligado (problema ${PROBLEM})")
else()
  list(JOIN _cutest_missing " " _msg)
  message(STATUS "Gencan: CUTEst desligado (faltam ${_msg}, rode source env_cutest.sh); sem cutestma_*")
endif()
