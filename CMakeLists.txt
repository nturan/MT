# CMakeList.txt : CMake project for MT, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)
enable_language(Fortran)
set(CMAKE_Fortran_COMPILER gfortran)
set(CMAKE_CXX_STANDARD 17)
set(MT "ttbar_NLO")

project (${MT})

# Add source to this project's executable.
add_subdirectory("src")
link_directories("${PROJECT_SOURCE_DIR}/libs/amps/")
link_directories("~/src/helas204/")
link_directories("${PROJECT_SOURCE_DIR}/libs/virtual/")
link_directories("~/src/QCDLoop/ql")
link_directories("~/src/QCDLoop/ff")
link_directories("/users/local/lhapdf/lib64")
add_executable (${MT} "src/main.cc" 
                      "src/SysInfo.cc"
                      "src/cross_sections.cc" 
                      "src/phase_space_generator.cc"
                      "src/parameters.cc"
                      "src/n_dimensional_integral.cc" 
                      "libs/ranlux-3.3/ranlxd.c" 
                      "src/helper.cc")

target_link_libraries(${MT} "-lMatrix -ldhelas -lgfortran")
target_link_libraries(${MT} "-lBSYPPTT -lqcdloop -lff")
target_link_libraries(${MT} "-lgsl -lgslcblas")
target_link_libraries(${MT} "-lLHAPDF")

target_include_directories(${MT} PUBLIC "include")
target_include_directories(${MT} PUBLIC "libs/virtual")
target_include_directories(${MT} PUBLIC "libs/amps")
target_include_directories(${MT} PUBLIC "/users/local/lhapdf/include")

target_compile_options(${MT} PUBLIC -Wall)
target_compile_options(${MT} PUBLIC -O2)

