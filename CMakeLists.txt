cmake_minimum_required(VERSION 3.10)
project(portillo_pcat_dnest)

set(CMAKE_CXX_STANDARD 14)

set(CMAKE_VERBOSE_MAKEFILE ON)
set(env{VERBOSE} 3)
set (CMAKE_CXX_FLAGS "-Wall -Wextra -pedantic -Wno-unused-parameter")
#set (CMAKE_CXX_FLAGS_RELEASE "-O2 -DNDEBUG -flto -march=native -g -fno-omit-frame-pointer")
# set (CMAKE_CXX_FLAGS_RELEASE "-O2 -DNDEBUG -flto -march=native -funroll-loops")
set (CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -g -funroll-loops -fno-omit-frame-pointer")

set (CMAKE_CXX_FLAGS_DEBUG "-O0 -DNDEBUG -g -fno-omit-frame-pointer")
#set (CMAKE_CXX_FLAGS_DEBUG "-O2 -DNDEBUG -flto -march=native -g -fno-omit-frame-pointer")


message(STATUS "INFO: Build is type ${CMAKE_BUILD_TYPE}")
message(STATUS "INFO: DNEST lib dir is ${DNEST3_LIBRARY}")
message(STATUS "INFO: RJOBJ lib dir is ${RJOBJ_LIBRARY}")


set(LIB_DEPENDS)
# ===== RJOJECT and DNEST3 =====
find_library(LIB_DNEST3 libdnest3.a HINTS ${DNEST3_LIBRARY})
list (APPEND LIB_DEPENDS ${LIB_DNEST3})
find_library(LIB_RJOBJ librjobject.a HINTS ${RJOBJ_LIBRARY})
list (APPEND LIB_DEPENDS ${LIB_RJOBJ})

if("${DNEST3_INCLUDES}" STREQUAL "")
  message(STATUS "ERR: DNEST3_INCLUDES, the path to dnest3 header includes is required!! Add the path, or rewrite this cmake file so it can find it automatically!")
endif("${DNEST3_INCLUDES}" STREQUAL "")
if("${RJOBJ_INCLUDES}" STREQUAL "")
  message(STATUS "ERR: RJOBJ_INCLUDES, the path to RJObjects header includes is required!! Add the path, or rewrite this cmake file so it can find it automatically!")
endif("${RJOBJ_INCLUDES}" STREQUAL "")
include_directories(${DNEST3_INCLUDES})
include_directories(${RJOBJ_INCLUDES})

# ===== BOOST =====
find_package (Boost COMPONENTS thread system REQUIRED)
include_directories (SYSTEM ${Boost_INCLUDE_DIRS})
list (APPEND LIB_DEPENDS ${Boost_LIBRARIES})

# ===== GSL =====
find_package (GSL REQUIRED)
include_directories (SYSTEM ${GSL_INCLUDES})
list (APPEND LIB_DEPENDS ${GSL_LIBRARIES})

include_directories(.)


#find_package(DNest3 REQUIRED)
#list(APPEND LIB_DEPENDS ${DNEST3_LIBRARY})
#include_directories (SYSTEM ${DNEST3_INCLUDES})

# If attempting to recompile everything from this project:
#file(GLOB DNEST3_SRCS ../eggplantbren_dnest3/include/*.h)
#file(GLOB DNEST3_SRCS ../eggplantbren_dnest3/src/*.cpp)
#file(GLOB RJOBJ_SRCS ../eggplantbren_rjobject/*.cpp)
#file(GLOB RJOBJ_SRCS ../eggplantbren_rjobject/*.h)
#file(GLOB RJOBJ_SRCS ../eggplantbren_rjobject/Distributions/*.cpp)
#file(GLOB RJOBJ_SRCS ../eggplantbren_rjobject/Distributions/*.h)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/buildflags.in
               ${CMAKE_BINARY_DIR}/buildflags)
file(READ ${CMAKE_BINARY_DIR}/buildflags CONTENTS)
file(WRITE ${CMAKE_BINARY_DIR}/buildflags "#### PCAT-DNEST FLAGS:\n${CONTENTS}")
file(READ ${DNEST3_LIBRARY}/../../buildflags CONTENTS)
file(APPEND ${CMAKE_BINARY_DIR}/buildflags "\n#### LIB DNEST3 FLAGS:\n${CONTENTS}")
file(READ ${RJOBJ_LIBRARY}/../../buildflags CONTENTS)
file(APPEND ${CMAKE_BINARY_DIR}/buildflags "\n#### LIB RJObject FLAGS:\n${CONTENTS}######\n")
file(READ ${CMAKE_BINARY_DIR}/buildflags CONTENTS)
file(WRITE ${CMAKE_BINARY_DIR}/buildinfo.h "constexpr char BUILDINFO[] =  R\"(${CONTENTS})\";\n")
execute_process(COMMAND git --no-pager log --graph "--pretty=format:%h -%d %s (%cD) <%an>%n" --abbrev-commit --date=relative --branches -n1
# execute_process(COMMAND git --no-pager log --graph --abbrev-commit --date=relative --branches -n1
                WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                OUTPUT_VARIABLE CONTENTS
                ERROR_VARIABLE CONTENTS)
file(APPEND ${CMAKE_BINARY_DIR}/buildinfo.h "constexpr char SRC_LAST_COMMIT[] =  R\"(${CONTENTS})\";\n")
include_directories(${CMAKE_BINARY_DIR})

add_executable(main
        main.cpp
        FermiData.cpp
        FermiData.h
        FermiModel.cpp
        FermiModel.h
        FermiModelOptions.cpp
        FermiModelOptions.h
        MyDistribution.cpp
        MyDistribution.h
        MyModel.cpp
        MyModel.h
        MyOptions.cpp
        MyOptions.h
        MyRJObject.h
        SloanData.cpp
        SloanData.h
        SloanModel.cpp
        SloanModel.h
        SloanModelOptions.cpp
        SloanModelOptions.h
        ${CMAKE_BINARY_DIR}/buildinfo.h
#        ${DNEST3_SRCS}
#        ${RJOBJ_SRCS}
        )


target_link_libraries(main ${LIB_DEPENDS})
