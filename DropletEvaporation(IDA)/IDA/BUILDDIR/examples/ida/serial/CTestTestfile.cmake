# CMake generated Testfile for 
# Source directory: C:/IDA/ida-6.5.1/examples/ida/serial
# Build directory: C:/IDA/BUILDDIR/examples/ida/serial
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
if(CTEST_CONFIGURATION_TYPE MATCHES "^([Dd][Ee][Bb][Uu][Gg])$")
  add_test(idaRoberts_dns "C:/IDA/BUILDDIR/bin/Debug/idaRoberts_dns.exe")
  set_tests_properties(idaRoberts_dns PROPERTIES  _BACKTRACE_TRIPLES "C:/IDA/ida-6.5.1/cmake/macros/SundialsAddTest.cmake;183;add_test;C:/IDA/ida-6.5.1/examples/ida/serial/CMakeLists.txt;92;sundials_add_test;C:/IDA/ida-6.5.1/examples/ida/serial/CMakeLists.txt;0;")
elseif(CTEST_CONFIGURATION_TYPE MATCHES "^([Rr][Ee][Ll][Ee][Aa][Ss][Ee])$")
  add_test(idaRoberts_dns "C:/IDA/BUILDDIR/bin/Release/idaRoberts_dns.exe")
  set_tests_properties(idaRoberts_dns PROPERTIES  _BACKTRACE_TRIPLES "C:/IDA/ida-6.5.1/cmake/macros/SundialsAddTest.cmake;183;add_test;C:/IDA/ida-6.5.1/examples/ida/serial/CMakeLists.txt;92;sundials_add_test;C:/IDA/ida-6.5.1/examples/ida/serial/CMakeLists.txt;0;")
elseif(CTEST_CONFIGURATION_TYPE MATCHES "^([Mm][Ii][Nn][Ss][Ii][Zz][Ee][Rr][Ee][Ll])$")
  add_test(idaRoberts_dns "C:/IDA/BUILDDIR/bin/MinSizeRel/idaRoberts_dns.exe")
  set_tests_properties(idaRoberts_dns PROPERTIES  _BACKTRACE_TRIPLES "C:/IDA/ida-6.5.1/cmake/macros/SundialsAddTest.cmake;183;add_test;C:/IDA/ida-6.5.1/examples/ida/serial/CMakeLists.txt;92;sundials_add_test;C:/IDA/ida-6.5.1/examples/ida/serial/CMakeLists.txt;0;")
elseif(CTEST_CONFIGURATION_TYPE MATCHES "^([Rr][Ee][Ll][Ww][Ii][Tt][Hh][Dd][Ee][Bb][Ii][Nn][Ff][Oo])$")
  add_test(idaRoberts_dns "C:/IDA/BUILDDIR/bin/RelWithDebInfo/idaRoberts_dns.exe")
  set_tests_properties(idaRoberts_dns PROPERTIES  _BACKTRACE_TRIPLES "C:/IDA/ida-6.5.1/cmake/macros/SundialsAddTest.cmake;183;add_test;C:/IDA/ida-6.5.1/examples/ida/serial/CMakeLists.txt;92;sundials_add_test;C:/IDA/ida-6.5.1/examples/ida/serial/CMakeLists.txt;0;")
else()
  add_test(idaRoberts_dns NOT_AVAILABLE)
endif()
