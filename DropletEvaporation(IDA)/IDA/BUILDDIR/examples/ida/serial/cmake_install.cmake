# Install script for directory: C:/IDA/ida-6.5.1/examples/ida/serial

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
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "C:/Program Files (x86)/SUNDIALS/examples/ida/serial/idaAnalytic_mels.c;C:/Program Files (x86)/SUNDIALS/examples/ida/serial/idaAnalytic_mels.out")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "C:/Program Files (x86)/SUNDIALS/examples/ida/serial" TYPE FILE FILES
    "C:/IDA/ida-6.5.1/examples/ida/serial/idaAnalytic_mels.c"
    "C:/IDA/ida-6.5.1/examples/ida/serial/idaAnalytic_mels.out"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "C:/Program Files (x86)/SUNDIALS/examples/ida/serial/idaFoodWeb_bnd.c;C:/Program Files (x86)/SUNDIALS/examples/ida/serial/idaFoodWeb_bnd.out")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "C:/Program Files (x86)/SUNDIALS/examples/ida/serial" TYPE FILE FILES
    "C:/IDA/ida-6.5.1/examples/ida/serial/idaFoodWeb_bnd.c"
    "C:/IDA/ida-6.5.1/examples/ida/serial/idaFoodWeb_bnd.out"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "C:/Program Files (x86)/SUNDIALS/examples/ida/serial/idaFoodWeb_kry.c;C:/Program Files (x86)/SUNDIALS/examples/ida/serial/idaFoodWeb_kry.out")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "C:/Program Files (x86)/SUNDIALS/examples/ida/serial" TYPE FILE FILES
    "C:/IDA/ida-6.5.1/examples/ida/serial/idaFoodWeb_kry.c"
    "C:/IDA/ida-6.5.1/examples/ida/serial/idaFoodWeb_kry.out"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "C:/Program Files (x86)/SUNDIALS/examples/ida/serial/idaHeat2D_bnd.c;C:/Program Files (x86)/SUNDIALS/examples/ida/serial/idaHeat2D_bnd.out")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "C:/Program Files (x86)/SUNDIALS/examples/ida/serial" TYPE FILE FILES
    "C:/IDA/ida-6.5.1/examples/ida/serial/idaHeat2D_bnd.c"
    "C:/IDA/ida-6.5.1/examples/ida/serial/idaHeat2D_bnd.out"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "C:/Program Files (x86)/SUNDIALS/examples/ida/serial/idaHeat2D_kry.c;C:/Program Files (x86)/SUNDIALS/examples/ida/serial/idaHeat2D_kry.out")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "C:/Program Files (x86)/SUNDIALS/examples/ida/serial" TYPE FILE FILES
    "C:/IDA/ida-6.5.1/examples/ida/serial/idaHeat2D_kry.c"
    "C:/IDA/ida-6.5.1/examples/ida/serial/idaHeat2D_kry.out"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "C:/Program Files (x86)/SUNDIALS/examples/ida/serial/idaKrylovDemo_ls.c;C:/Program Files (x86)/SUNDIALS/examples/ida/serial/idaKrylovDemo_ls.out;C:/Program Files (x86)/SUNDIALS/examples/ida/serial/idaKrylovDemo_ls_1.out;C:/Program Files (x86)/SUNDIALS/examples/ida/serial/idaKrylovDemo_ls_2.out")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "C:/Program Files (x86)/SUNDIALS/examples/ida/serial" TYPE FILE FILES
    "C:/IDA/ida-6.5.1/examples/ida/serial/idaKrylovDemo_ls.c"
    "C:/IDA/ida-6.5.1/examples/ida/serial/idaKrylovDemo_ls.out"
    "C:/IDA/ida-6.5.1/examples/ida/serial/idaKrylovDemo_ls_1.out"
    "C:/IDA/ida-6.5.1/examples/ida/serial/idaKrylovDemo_ls_2.out"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "C:/Program Files (x86)/SUNDIALS/examples/ida/serial/idaKrylovDemo_ls.c;C:/Program Files (x86)/SUNDIALS/examples/ida/serial/idaKrylovDemo_ls.out;C:/Program Files (x86)/SUNDIALS/examples/ida/serial/idaKrylovDemo_ls_1.out;C:/Program Files (x86)/SUNDIALS/examples/ida/serial/idaKrylovDemo_ls_2.out")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "C:/Program Files (x86)/SUNDIALS/examples/ida/serial" TYPE FILE FILES
    "C:/IDA/ida-6.5.1/examples/ida/serial/idaKrylovDemo_ls.c"
    "C:/IDA/ida-6.5.1/examples/ida/serial/idaKrylovDemo_ls.out"
    "C:/IDA/ida-6.5.1/examples/ida/serial/idaKrylovDemo_ls_1.out"
    "C:/IDA/ida-6.5.1/examples/ida/serial/idaKrylovDemo_ls_2.out"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "C:/Program Files (x86)/SUNDIALS/examples/ida/serial/idaKrylovDemo_ls.c;C:/Program Files (x86)/SUNDIALS/examples/ida/serial/idaKrylovDemo_ls.out;C:/Program Files (x86)/SUNDIALS/examples/ida/serial/idaKrylovDemo_ls_1.out;C:/Program Files (x86)/SUNDIALS/examples/ida/serial/idaKrylovDemo_ls_2.out")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "C:/Program Files (x86)/SUNDIALS/examples/ida/serial" TYPE FILE FILES
    "C:/IDA/ida-6.5.1/examples/ida/serial/idaKrylovDemo_ls.c"
    "C:/IDA/ida-6.5.1/examples/ida/serial/idaKrylovDemo_ls.out"
    "C:/IDA/ida-6.5.1/examples/ida/serial/idaKrylovDemo_ls_1.out"
    "C:/IDA/ida-6.5.1/examples/ida/serial/idaKrylovDemo_ls_2.out"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "C:/Program Files (x86)/SUNDIALS/examples/ida/serial/idaRoberts_dns.c;C:/Program Files (x86)/SUNDIALS/examples/ida/serial/idaRoberts_dns.out")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "C:/Program Files (x86)/SUNDIALS/examples/ida/serial" TYPE FILE FILES
    "C:/IDA/ida-6.5.1/examples/ida/serial/idaRoberts_dns.c"
    "C:/IDA/ida-6.5.1/examples/ida/serial/idaRoberts_dns.out"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "C:/Program Files (x86)/SUNDIALS/examples/ida/serial/idaSlCrank_dns.c;C:/Program Files (x86)/SUNDIALS/examples/ida/serial/idaSlCrank_dns.out")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "C:/Program Files (x86)/SUNDIALS/examples/ida/serial" TYPE FILE FILES
    "C:/IDA/ida-6.5.1/examples/ida/serial/idaSlCrank_dns.c"
    "C:/IDA/ida-6.5.1/examples/ida/serial/idaSlCrank_dns.out"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "C:/Program Files (x86)/SUNDIALS/examples/ida/serial/README")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "C:/Program Files (x86)/SUNDIALS/examples/ida/serial" TYPE FILE FILES "C:/IDA/ida-6.5.1/examples/ida/serial/README")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "C:/Program Files (x86)/SUNDIALS/examples/ida/serial/idaRoberts_dns_stats.csv")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "C:/Program Files (x86)/SUNDIALS/examples/ida/serial" TYPE FILE FILES "C:/IDA/ida-6.5.1/examples/ida/serial/idaRoberts_dns_stats.csv")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "C:/Program Files (x86)/SUNDIALS/examples/ida/serial/CMakeLists.txt")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "C:/Program Files (x86)/SUNDIALS/examples/ida/serial" TYPE FILE FILES "C:/IDA/BUILDDIR/examples/ida/serial/CMakeLists.txt")
endif()

