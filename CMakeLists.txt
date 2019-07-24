# minimum cmake version requirement
cmake_minimum_required(VERSION 3.6)

# TODO: Make sure source files are added into the msvc project to compile correctly

project(snwpac LANGUAGES C CXX)

string(TOUPPER ${PROJECT_NAME} NAME_ID)

option(${NAME_ID}_BUILD_STATIC
	"Build a static version of this library."
	ON
)

# TODO: Automate this
option(${NAME_ID}_INSTALL
	"Install the include headers and the static target (if build) of ${PROJECT_NAME}."
	ON
)

 option(${NAME_ID}_TESTS
 	"Build tests for ${PROJECT_NAME}."
 	ON
 )

set(PROJECT_SOURCE_FILES
	src/snwpac.cpp
	src/writer.cpp
	src/reader.cpp
	src/quick_sort.cpp
	src/bitonic_sort.cpp
)

set(PROJECT_HEADER_FILES
	include/snwpac.hpp
	include/writer.h
	include/writer.inl
	include/reader.h
	include/reader.inl
	include/quick_sort.h
	include/bitonic_sort.h
	include/wavelet.h
	include/wavelet.inl
    include/gsl/gsl_wavelet.h
    include/gsl/gsl_types.h
    include/gsl/gsl_types.h
    include/gsl/bspline.c
    include/gsl/wavelet.c
)

# interface header files for install
set(PROJECT_INTERFACE_HEADER_FILES
	include/snwpac.hpp
)

# External Dependencies
set (EXTRA_LIBS "")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")

OPTION (${NAME_ID}_USE_OPENMP "Use OpenMP" ON)
if(${NAME_ID}_USE_OPENMP)
    find_package(OpenMP)
    if (OPENMP_FOUND)
        set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
  endif()
endif()

# 3rdparty libraries
add_subdirectory(3rdparty/lz4)
set (EXTRA_LIBS ${EXTRA_LIBS} lz4 )

# Testing files
if(${NAME_ID}_TESTS)
	add_subdirectory(tests)
endif()

# add library targets
if(${NAME_ID}_BUILD_STATIC)
	add_library(${PROJECT_NAME} STATIC
		${PROJECT_SOURCE_FILES}
		${PROJECT_HEADER_FILES}
		${PROJECT_INTERFACE_HEADER_FILES}
	)
else()
	add_library(${PROJECT_NAME} SHARED
		${PROJECT_SOURCE_FILES}
		${PROJECT_HEADER_FILES}
		${PROJECT_INTERFACE_HEADER_FILES}
	)
endif()

# set c++ standard
set_target_properties(${PROJECT_NAME} PROPERTIES
	CXX_STANDARD 11
	CXX_STANDARD_REQUIRED ON
	CXX_EXTENSIONS OFF
)

# Enable fPIC for linux
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# add include directories
target_include_directories(${PROJECT_NAME}
	PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# add library dependencies
target_link_libraries(${PROJECT_NAME} ${EXTRA_LIBS} )

# install instructions
if(${NAME_ID}_BUILD_STATIC AND ${NAME_ID}_INSTALL)
	install(TARGETS ${PROJECT_NAME}
		DESTINATION lib
	)
elseif(NOT ${NAME_ID}_BUILD_STATIC)
	if(WIN32)
		install(TARGETS ${PROJECT_NAME}
			DESTINATION bin
		)
	else()
		install(TARGETS ${PROJECT_NAME}
			DESTINATION lib
		)
	endif()
endif()

if(${NAME_ID}_INSTALL)
	install(FILES
		${PROJECT_INTERFACE_HEADER_FILES}
		DESTINATION include
	)
endif()
