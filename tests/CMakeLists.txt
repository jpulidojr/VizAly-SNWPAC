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

# add executable target
add_executable(${PROJECT_NAME}
    ${PROJECT_SOURCE_FILES}
    ${PROJECT_HEADER_FILES}
)

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

if(HDF5_FOUND)
    target_link_libraries(${PROJECT_NAME} ${HDF5_LIBS})
    add_definitions(-DUSE_HDF5)
else()
    target_link_libraries(${PROJECT_NAME})
endif()

# install instructions
install(TARGETS ${PROJECT_NAME}
    DESTINATION tests
)