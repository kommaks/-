# Install script for directory: C:/IDA/ida-6.5.1/src/sundials

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "C:/Program Files (x86)/SUNDIALS")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  MESSAGE("
Install shared components
")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Dd][Ee][Bb][Uu][Gg])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "C:/IDA/BUILDDIR/bin/Debug/sundials_generic.lib")
  elseif(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Rr][Ee][Ll][Ee][Aa][Ss][Ee])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "C:/IDA/BUILDDIR/bin/Release/sundials_generic.lib")
  elseif(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Mm][Ii][Nn][Ss][Ii][Zz][Ee][Rr][Ee][Ll])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "C:/IDA/BUILDDIR/bin/MinSizeRel/sundials_generic.lib")
  elseif(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Rr][Ee][Ll][Ww][Ii][Tt][Hh][Dd][Ee][Bb][Ii][Nn][Ff][Oo])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "C:/IDA/BUILDDIR/bin/RelWithDebInfo/sundials_generic.lib")
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Dd][Ee][Bb][Uu][Gg])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY OPTIONAL FILES "C:/IDA/BUILDDIR/bin/Debug/sundials_generic.lib")
  elseif(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Rr][Ee][Ll][Ee][Aa][Ss][Ee])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY OPTIONAL FILES "C:/IDA/BUILDDIR/bin/Release/sundials_generic.lib")
  elseif(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Mm][Ii][Nn][Ss][Ii][Zz][Ee][Rr][Ee][Ll])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY OPTIONAL FILES "C:/IDA/BUILDDIR/bin/MinSizeRel/sundials_generic.lib")
  elseif(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Rr][Ee][Ll][Ww][Ii][Tt][Hh][Dd][Ee][Bb][Ii][Nn][Ff][Oo])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY OPTIONAL FILES "C:/IDA/BUILDDIR/bin/RelWithDebInfo/sundials_generic.lib")
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Dd][Ee][Bb][Uu][Gg])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "C:/IDA/BUILDDIR/bin/Debug/sundials_generic.dll")
  elseif(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Rr][Ee][Ll][Ee][Aa][Ss][Ee])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "C:/IDA/BUILDDIR/bin/Release/sundials_generic.dll")
  elseif(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Mm][Ii][Nn][Ss][Ii][Zz][Ee][Rr][Ee][Ll])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "C:/IDA/BUILDDIR/bin/MinSizeRel/sundials_generic.dll")
  elseif(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Rr][Ee][Ll][Ww][Ii][Tt][Hh][Dd][Ee][Bb][Ii][Nn][Ff][Oo])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "C:/IDA/BUILDDIR/bin/RelWithDebInfo/sundials_generic.dll")
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/sundials" TYPE FILE FILES
    "C:/IDA/ida-6.5.1/include/sundials/sundials_base.hpp"
    "C:/IDA/ida-6.5.1/include/sundials/sundials_band.h"
    "C:/IDA/ida-6.5.1/include/sundials/sundials_context.h"
    "C:/IDA/ida-6.5.1/include/sundials/sundials_context.hpp"
    "C:/IDA/ida-6.5.1/include/sundials/sundials_convertibleto.hpp"
    "C:/IDA/ida-6.5.1/include/sundials/sundials_dense.h"
    "C:/IDA/ida-6.5.1/include/sundials/sundials_direct.h"
    "C:/IDA/ida-6.5.1/include/sundials/sundials_iterative.h"
    "C:/IDA/ida-6.5.1/include/sundials/sundials_linearsolver.h"
    "C:/IDA/ida-6.5.1/include/sundials/sundials_linearsolver.hpp"
    "C:/IDA/ida-6.5.1/include/sundials/sundials_math.h"
    "C:/IDA/ida-6.5.1/include/sundials/sundials_matrix.h"
    "C:/IDA/ida-6.5.1/include/sundials/sundials_matrix.hpp"
    "C:/IDA/ida-6.5.1/include/sundials/sundials_memory.h"
    "C:/IDA/ida-6.5.1/include/sundials/sundials_mpi_types.h"
    "C:/IDA/ida-6.5.1/include/sundials/sundials_nonlinearsolver.h"
    "C:/IDA/ida-6.5.1/include/sundials/sundials_nonlinearsolver.hpp"
    "C:/IDA/ida-6.5.1/include/sundials/sundials_nvector.h"
    "C:/IDA/ida-6.5.1/include/sundials/sundials_nvector.hpp"
    "C:/IDA/ida-6.5.1/include/sundials/sundials_profiler.h"
    "C:/IDA/ida-6.5.1/include/sundials/sundials_logger.h"
    "C:/IDA/ida-6.5.1/include/sundials/sundials_types.h"
    "C:/IDA/ida-6.5.1/include/sundials/sundials_version.h"
    )
endif()
