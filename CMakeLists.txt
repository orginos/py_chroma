cmake_minimum_required(VERSION 3.18)

project(py_chroma LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(MPI REQUIRED)
find_package(HDF5 REQUIRED)

set(CHROMA_PREFIX "" CACHE PATH "Chroma install prefix")
set(QDP_PREFIX "" CACHE PATH "QDP++ install prefix")
set(PRIMME_PREFIX "" CACHE PATH "PRIMME install prefix")

if (CHROMA_PREFIX STREQUAL "")
  get_filename_component(_default_chroma "${CMAKE_CURRENT_LIST_DIR}/../install/chroma-qdpxx-double-nd4" ABSOLUTE)
  set(CHROMA_PREFIX "${_default_chroma}" CACHE PATH "Chroma install prefix" FORCE)
endif()

if (QDP_PREFIX STREQUAL "")
  get_filename_component(_default_qdp "${CMAKE_CURRENT_LIST_DIR}/../install/qdpxx-double-nd4" ABSOLUTE)
  set(QDP_PREFIX "${_default_qdp}" CACHE PATH "QDP++ install prefix" FORCE)
endif()

if (PRIMME_PREFIX STREQUAL "")
  get_filename_component(_default_primme "${CMAKE_CURRENT_LIST_DIR}/../install/primme-cpu" ABSOLUTE)
  set(PRIMME_PREFIX "${_default_primme}" CACHE PATH "PRIMME install prefix" FORCE)
endif()

set(QDP_CONFIG "${QDP_PREFIX}/bin/qdp++-config")
if (NOT EXISTS "${QDP_CONFIG}")
  message(FATAL_ERROR "qdp++-config not found at ${QDP_CONFIG}. Set QDP_PREFIX.")
endif()

execute_process(
  COMMAND "${QDP_CONFIG}" --cxxflags
  OUTPUT_VARIABLE QDP_CXXFLAGS
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
execute_process(
  COMMAND "${QDP_CONFIG}" --ldflags
  OUTPUT_VARIABLE QDP_LDFLAGS
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
execute_process(
  COMMAND "${QDP_CONFIG}" --libs
  OUTPUT_VARIABLE QDP_LIBS
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

separate_arguments(QDP_CXXFLAGS_LIST NATIVE_COMMAND "${QDP_CXXFLAGS}")
separate_arguments(QDP_LDFLAGS_LIST NATIVE_COMMAND "${QDP_LDFLAGS}")
separate_arguments(QDP_LIBS_LIST NATIVE_COMMAND "${QDP_LIBS}")

add_library(pychroma SHARED
  src/pychroma_core.cpp
  src/pychroma_hmc.cpp
  src/pychroma_inline.cpp
  src/pychroma_rect_wloops.cpp
  src/pychroma_smd.cpp
  src/pychroma_state.cpp
  src/pychroma_util.cpp
)

target_include_directories(pychroma PRIVATE
  "${CHROMA_PREFIX}/include"
  "${CMAKE_CURRENT_LIST_DIR}/../chroma/lib"
  ${HDF5_INCLUDE_DIRS}
)

target_compile_options(pychroma PRIVATE
  ${QDP_CXXFLAGS_LIST}
)

target_link_options(pychroma PRIVATE
  ${QDP_LDFLAGS_LIST}
)

target_link_libraries(pychroma PRIVATE
  "${CHROMA_PREFIX}/lib/libchroma.a"
  "${CHROMA_PREFIX}/lib/libqdp-lapack.a"
  "${PRIMME_PREFIX}/lib/libprimme.a"
  MPI::MPI_CXX
  ${QDP_LIBS_LIST}
  ${HDF5_LIBRARIES}
)

set_target_properties(pychroma PROPERTIES
  OUTPUT_NAME "pychroma"
)

install(TARGETS pychroma
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  RUNTIME DESTINATION bin
)
install(DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/py_chroma"
  DESTINATION python
)
install(PROGRAMS "${CMAKE_CURRENT_LIST_DIR}/scripts/pychroma_env.sh"
  DESTINATION bin
)
