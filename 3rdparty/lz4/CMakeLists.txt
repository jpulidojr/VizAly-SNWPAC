# external project CMakeLists.txt template

# minimum cmake version requirement
cmake_minimum_required(VERSION 3.6)

if(NOT lz4)

# project name
project(lz4)

# enable external project cmake module
include(ExternalProject)

if(WIN32)
    # download external project binary
    ExternalProject_Add( ${PROJECT_NAME}-external
        URL https://github.com/lz4/lz4/releases/download/v1.8.2/lz4_v1_8_2_win64.zip
        PREFIX ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-external
        DEPENDS ""
        PATCH_COMMAND ""
        CONFIGURE_COMMAND ""
        # disable build
        BUILD_COMMAND ""
        # disable install
        INSTALL_COMMAND ""
        # disable anoying constant checks during build
        UPDATE_DISCONNECTED TRUE
    )
else()
    # Build external project from source
    ExternalProject_Add( ${PROJECT_NAME}-external
        GIT_REPOSITORY https://github.com/lz4/lz4.git
        GIT_TAG v1.8.2
        PREFIX ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-external
        DEPENDS ""
        PATCH_COMMAND ""
        CONFIGURE_COMMAND ""
        BUILD_COMMAND "make" "CFLAGS=-fPIC" 
        BUILD_IN_SOURCE 1 
        # disable install step
        INSTALL_COMMAND ""
        # disable anoying constant checks during build
        UPDATE_DISCONNECTED TRUE
        # Enable fPIC for static->dynamic builds
        #CMAKE_ARGS "-DCMAKE_POSITION_INDEPENDENT_CODE:BOOL=true"
    )
endif()

message ("cmake_current_bin_dir ${CMAKE_CURRENT_BINARY_DIR} ...")

# get external project source and binary directories
ExternalProject_Get_Property(${PROJECT_NAME}-external
    SOURCE_DIR
    BINARY_DIR
)

message ("binary_dir ${BINARY_DIR} ...")
message ("source_dir ${SOURCE_DIR} ...")

# create source directory with include path.
# this is necessary to configure interface include directories.
# otherwise this doesnt get set in the first pass
execute_process(
    COMMAND ${CMAKE_COMMAND} -E make_directory ${SOURCE_DIR}/lib
)

# create imported target to be used as a dependency by programs
add_library(${PROJECT_NAME} STATIC IMPORTED GLOBAL)

# set project properties
if (UNIX)
set_target_properties(${PROJECT_NAME} PROPERTIES
		IMPORTED_LOCATION ${BINARY_DIR}/lib/liblz4.a
)
else() #WIN32
set_target_properties(${PROJECT_NAME} PROPERTIES
		IMPORTED_LOCATION ${SOURCE_DIR}/dll/liblz4.lib
        #IMPORTED_LOCATION ${SOURCE_DIR}/static/liblz4_static.lib #win32 doesnt work with prebuilt so
)
#rename .dll
if(EXISTS "${SOURCE_DIR}/dll/liblz4.so.1.8.2.dll")
	file(RENAME ${SOURCE_DIR}/dll/liblz4.so.1.8.2.dll ${SOURCE_DIR}/dll/liblz4.dll)
endif()

endif()

# add external project as dependency
add_dependencies(${PROJECT_NAME}
    ${PROJECT_NAME}-external
)

# add external project include directories
# ISSUE: https://gitlab.kitware.com/cmake/cmake/issues/15689
# FIXED: https://gitlab.kitware.com/cmake/cmake/merge_requests/1264
# CMake 3.11+
# target_include_directories(${PROJECT_NAME}
#    INTERFACE ${SOURCE_DIR}/lib
#)

# Should work with older CMakes. Reqs. ln 46
if (UNIX)
set_property(TARGET ${PROJECT_NAME} APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${SOURCE_DIR}/lib)
else() #WIN32
set_property(TARGET ${PROJECT_NAME} APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${SOURCE_DIR}/include)
endif()

endif (NOT lz4)