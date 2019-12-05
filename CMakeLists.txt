cmake_minimum_required(VERSION 3.10)
project(portillo_pcat_dnest_zmod)

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

add_subdirectory(${PROJECT_SOURCE_DIR}/lib_dnest3)
add_subdirectory(${PROJECT_SOURCE_DIR}/lib_rjobject)

set(LIB_DEPENDS)
# ===== BOOST =====
find_package (Boost COMPONENTS thread system REQUIRED)
include_directories (SYSTEM ${Boost_INCLUDE_DIRS})
list (APPEND LIB_DEPENDS ${Boost_LIBRARIES})

# ===== GSL =====
find_package (GSL REQUIRED)
include_directories (SYSTEM ${GSL_INCLUDES})
list (APPEND LIB_DEPENDS ${GSL_LIBRARIES})

include_directories(.)
include_directories(lib_rjobject)
include_directories(lib_rjobject/Distributions)
include_directories(lib_dnest3/include)


configure_file(${CMAKE_CURRENT_SOURCE_DIR}/buildflags.in
               ${CMAKE_BINARY_DIR}/buildflags)
file(READ ${CMAKE_BINARY_DIR}/buildflags CONTENTS)
file(WRITE ${CMAKE_BINARY_DIR}/buildflags "#### PCAT-DNEST FLAGS:\n${CONTENTS}")
file(READ ${CMAKE_BINARY_DIR}/buildflags_dnest3 CONTENTS)
file(APPEND ${CMAKE_BINARY_DIR}/buildflags "\n#### LIB DNEST3 FLAGS:\n${CONTENTS}")
file(READ ${CMAKE_BINARY_DIR}/buildflags_rjobject CONTENTS)
file(APPEND ${CMAKE_BINARY_DIR}/buildflags "\n#### LIB RJObject FLAGS:\n${CONTENTS}######\n")
file(READ ${CMAKE_BINARY_DIR}/buildflags CONTENTS)
file(WRITE ${CMAKE_BINARY_DIR}/buildinfo.h "constexpr char BUILDINFO[] =  R\"(${CONTENTS})\";\n")
execute_process(COMMAND git --no-pager log --graph "--pretty=format:%h -%d %s (%cD) <%an>%n" --abbrev-commit -n1
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

# ===== DNEST3 =====
target_link_libraries(main lib::dnest3)
target_link_libraries(main lib::rjobject)

target_link_libraries(main ${LIB_DEPENDS})
