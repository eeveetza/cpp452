cmake_minimum_required(VERSION 3.24)
project(itu_r_p452)

if(ENABLE_CODE_COVERAGE)
  include(${CMAKE_CURRENT_SOURCE_DIR}/CodeCoverage.cmake)
  option(ENABLE_CODE_COVERAGE "Enable code coverage calculation" OFF)
  message(STATUS "Enabled code coverage analysis")
  append_coverage_compiler_flags()
  set(COVERAGE_EXCLUDES
    'build/*'
    '/usr/include/*'
)
endif()

if(MSVC)
  add_compile_options(/W4)
else()
  add_compile_options(-Wall -Wextra -Wpedantic)
endif(MSVC)

enable_testing()

add_subdirectory(src)
add_subdirectory(tests)
