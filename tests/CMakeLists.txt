# executable CMakeLists.txt template

# minimum cmake version requirement
cmake_minimum_required(VERSION 3.6)

# project name
project(lossywave_test CXX)

# add source files
set(PROJECT_SOURCE_FILES
    lossywave_test.cpp
)

# add header files
set(PROJECT_HEADER_FILES
)

# add executable target
add_executable(${PROJECT_NAME}
    ${PROJECT_SOURCE_FILES}
    ${PROJECT_HEADER_FILES}
)
set(LW_PROGRAMS_BUILT ${PROJECT_NAME})

# set c++ standard
set_target_properties(${PROJECT_NAME} PROPERTIES
    CXX_STANDARD 11
    CXX_STANDARD_REQUIRED ON
    CXX_EXTENSIONS OFF
)

# add include directories
target_include_directories(${PROJECT_NAME}
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# add library dependencies
target_link_libraries(${PROJECT_NAME}
    PRIVATE lossywave
)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()
    
#Add additional sample projects
# add source files
set(DPROJECT_SOURCE_FILES
    LWDecomp.cpp
)
add_executable(LWDecomp ${DPROJECT_SOURCE_FILES})
list(APPEND LW_PROGRAMS_BUILT LWDecomp)
target_link_libraries(LWDecomp PRIVATE lossywave)
target_include_directories(LWDecomp PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
set_target_properties(LWDecomp PROPERTIES
    CXX_STANDARD 11
    CXX_STANDARD_REQUIRED ON
    CXX_EXTENSIONS OFF
)

option(TEST_ENABLE_HDF5 "Enable HDF5 support" OFF)
if(TEST_ENABLE_HDF5)
    find_package(HDF5 COMPONENTS C CXX HL NO_MODULE)
    if(HDF5_FOUND)
        include_directories(${HDF5_INCLUDE_DIR})
        set(HDF5_LIB_DIR ${HDF5_LIB_DIR} CACHE PATH "Location of HDF5 libraries")
        link_directories(${HDF5_LIB_DIR})
        set(HDF5_LIBS ${HDF5_LIBRARIES} libhdf5 libhdf5_cpp libszip libzlib )
    endif()
endif()

if(HDF5_FOUND)
    target_link_libraries(LWDecomp ${HDF5_LIBS})
    add_definitions(-DUSE_HDF5)
else()
    #target_link_libraries(${PROJECT_NAME})
endif()


# install instructions
install(TARGETS ${LW_PROGRAMS_BUILT}
    DESTINATION tests
)