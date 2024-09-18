add_library(p452_library STATIC P452.cpp)
add_library(p452::library ALIAS p452_library)

target_include_directories(p452_library PUBLIC ${PROJECT_SOURCE_DIR}/include)
