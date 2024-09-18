include(FetchContent)
FetchContent_Declare(
  googletest
  # Specify the commit you depend on and update it regularly.
  URL https://github.com/google/googletest/archive/5376968f6948923e2411081fd9372e71a59d8e77.zip
  FIND_PACKAGE_ARGS NAMES GTest
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

add_executable(
  p452_test
  test1.cpp
  test2.cpp
  test3.cpp
  test4.cpp
  test5.cpp
  test6.cpp
  test7.cpp
  test8.cpp
  test9.cpp
  test10.cpp)

target_compile_features(p452_test PUBLIC cxx_std_17)

target_link_libraries(p452_test PRIVATE GTest::gtest GTest::gtest_main
                                        p452_library)

add_test(All452Tests p452_test)
