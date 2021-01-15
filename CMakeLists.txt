# CMakeList.txt : CMake project for MT, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)

enable_language(Fortran)
set(CMAKE_CXX_STANDARD 17)

set(MT "ttbar_NLO")

project (${MT})

# Add source to this project's executable.
add_subdirectory("src")
add_subdirectory("libs/virtual")
add_subdirectory("libs/amps")
add_subdirectory("libs/ranlux-3.3")

link_directories("/home/turan/src/helas204/")
link_directories("/home/turan/src/QCDLoop/ql")
link_directories("/home/turan/src/QCDLoop/ff")
link_directories("/users/local/lhapdf/lib64")


add_executable (${MT} "src/main.cc" 
                      "src/SysInfo.cc"
                      "src/cross_sections.cc" 
                      "src/phase_space_generator.cc"
                      "src/parameters.cc"
                      "src/n_dimensional_integral.cc"
                      "src/helper.cc" 
                      "src/histogram.cc" 
                      "include/cxxopts.hpp" 
                      "src/event_generator.cc" 
                      "src/matrix_element_method.cc")

target_link_libraries(${MT} virtualAmps)
target_link_libraries(${MT} amps)
target_link_libraries(${MT} ranlux)

target_link_libraries(${MT} "-lgsl -lgslcblas")
target_link_libraries(${MT} "-lLHAPDF")

target_include_directories(${MT} PUBLIC "include")
target_include_directories(${MT} PUBLIC "/users/local/lhapdf/include")

target_compile_options(${MT} PUBLIC -Wall)
target_compile_options(${MT} PUBLIC -O2)

