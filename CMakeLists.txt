cmake_minimum_required(VERSION 3.2)
project(ttk VERSION 0.9.2 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)

add_subdirectory(core)
option(TTK_BUILD_STANDALONE_APPS "Build the TTK Standalone Applications" ON)
if (TTK_BUILD_STANDALONE_APPS)
	add_subdirectory(standalone)
endif()

