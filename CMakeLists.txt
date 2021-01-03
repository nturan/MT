# CMakeList.txt : CMake project for MT, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)

enable_language(Fortran)

set(MT "ttbar_NLO")

project (${MT})

# Add source to this project's executable.
add_subdirectory("src")
link_directories("${PROJECT_SOURCE_DIR}/libs/amps/")
link_directories("/home/turan/src/helas204/")
link_directories("${PROJECT_SOURCE_DIR}/libs/virtual/")
link_directories("/home/turan/src/QCDLoop/ql")
link_directories("/home/turan/src/QCDLoop/ff")
link_directories("/users/local/lhapdf/lib64")

add_library(BSYppttModule "libs/virtual/bsyspinors.f"
                          "libs/virtual/bsyspinit.f"
                          "libs/virtual/ampggint.f"
                          "libs/virtual/ampgg.f"
                          "libs/virtual/ampqqint.f"
                          "libs/virtual/ampqq.f"
                          "libs/virtual/bsypptt.f")

target_link_libraries(BSYppttModule "-lqcdloop" "-lff")
add_executable (${MT} "src/main.cc" 
                      "src/SysInfo.cc"
                      "src/cross_sections.cc" 
                      "src/phase_space_generator.cc"
                      "src/parameters.cc"
                      "src/n_dimensional_integral.cc" 
                      "libs/ranlux-3.3/ranlxd.c" 
                      "src/helper.cc")


#target_include_directories(${MT} PUBLIC "libs/virtual")
target_link_libraries(${MT} BSYppttModule)

target_link_libraries(${MT} "-lMatrix -ldhelas")
target_link_libraries(${MT} "-lgsl -lgslcblas")
target_link_libraries(${MT} "-lLHAPDF")

target_include_directories(${MT} PUBLIC "include")
target_include_directories(${MT} PUBLIC "libs/amps")
target_include_directories(${MT} PUBLIC "/users/local/lhapdf/include")

target_compile_options(${MT} PUBLIC -Wall)
target_compile_options(${MT} PUBLIC -O2)

