cmake_minimum_required(VERSION 3.5)
project(ttk)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

add_subdirectory(core)
option(TTK_BUILD_STANDALONE_APPS "Build the TTK Standalone Applications" OFF)
if (TTK_BUILD_STANDALONE_APPS)
	add_subdirectory(standalone)
endif()

