cmake_minimum_required(VERSION 3.15)
project(portillo_pcat_dnest)

set(CMAKE_CXX_STANDARD 14)

set(CMAKE_VERBOSE_MAKEFILE ON)
set(env{VERBOSE} 3)

set(LIB_DEPENDS)

# ===== BOOST =====
find_package (Boost COMPONENTS thread system REQUIRED)
include_directories (SYSTEM ${Boost_INCLUDE_DIRS})
list (APPEND LIB_DEPENDS ${Boost_LIBRARIES})

include_directories(.)
include_directories(../eggplantbren_dnest3/lib/lib)
include_directories(../eggplantbren_dnest3/lib/include/dnest3)
include_directories(../eggplantbren_rjobject/lib/lib)
include_directories(../eggplantbren_rjobject/lib/include/rjobject)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -v" )
add_executable(main
        main.cpp
        FermiData.cpp
        FermiModel.cpp
        FermiModelOptions.cpp
        Makefile
        MyDistribution.cpp
        MyModel.cpp
        MyOptions.cpp
        SloanData.cpp
        SloanModel.cpp
        SloanModelOptions.cpp)
target_link_libraries(main ${LIB_DEPENDS})
