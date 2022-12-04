# CMakeList.txt : CMake project for MT, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)
enable_language(Fortran)
set(CMAKE_CXX_STANDARD 17)
set(MT "ttbar_NLO")
include(ExternalProject)

project (${MT})

# Add source to this project's executable.
add_subdirectory("src")
add_subdirectory("libs/virtual")
add_subdirectory("libs/amps")
add_subdirectory("libs/ranlux-3.3")
add_subdirectory("libs/sysinfo")
add_subdirectory("libs/vegas")
add_subdirectory("libs/cxxopts")

link_directories("libs/helas204/")
link_directories("libs/QCDLoop/ql")
link_directories("libs/QCDLoop/ff")
link_directories("/users/local/lhapdf/lib64")

list(APPEND PROGRAM_SOURCES "src/cross_sections.cc" 
                      "src/phase_space_generator.cc"
                      "src/parameters.cc"
                      "src/n_dimensional_integral.cc"
                      "src/helper.cc" 
                      "src/histogram.cc" 
                      "src/event_generator.cc" 
                      "src/matrix_element_method.cc")

add_executable (${MT} "src/main.cc" ${PROGRAM_SOURCES})



include(libs/helas204/helas.cmake)
include(libs/QCDLoop/ff/ff.cmake)
include(libs/QCDLoop/ql/qcdloop.cmake)
target_link_libraries(${MT} virtualAmps)
target_link_libraries(${MT} amps)
target_link_libraries(${MT} ranlux)
target_link_libraries(${MT} vegas)
target_link_libraries(${MT} sysinfo)
target_link_libraries(${MT} cxxopts)

target_link_libraries(${MT} "-lgsl -lgslcblas")
target_link_libraries(${MT} "-lLHAPDF")
target_link_libraries(${MT} "-lstdc++")

target_include_directories(${MT} PUBLIC "include")
target_include_directories(${MT} PUBLIC "/usr/include")
target_include_directories(${MT} PUBLIC "libs/cxxopts")
target_include_directories(${MT} PUBLIC "libs/vegas")
target_include_directories(${MT} PUBLIC "libs/sysinfo")
target_include_directories(${MT} PUBLIC "/users/local/lhapdf/include")

target_compile_options(${MT} PUBLIC -Wall)
target_compile_options(${MT} PUBLIC -O2)

