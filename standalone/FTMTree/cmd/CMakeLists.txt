cmake_minimum_required(VERSION 3.2)

project(ttkFTMTreeCmd)

if(NOT DEFINED TTK_BUILD_FTMTREE_FILTER)
  # Allows to disable each filter
  option(TTK_BUILD_FTMTREE_FILTER "Build the filter FTMTree" ON)
  mark_as_advanced(TTK_BUILD_FTMTREE_FILTER)
endif()

if(${TTK_BUILD_FTMTREE_FILTER})
  add_executable(${PROJECT_NAME} main.cpp)
  target_link_libraries(${PROJECT_NAME}
    PRIVATE
      ttkFTMTree
      ttkProgramBase
      )
  set_target_properties(${PROJECT_NAME}
    PROPERTIES
      INSTALL_RPATH
        "${CMAKE_INSTALL_PREFIX}/lib64/TopologyToolKit;${CMAKE_INSTALL_PREFIX}/lib/ttk")
  install(
    TARGETS
      ${PROJECT_NAME}
    RUNTIME DESTINATION
      ${TTK_INSTALL_BINARY_DIR}
      )
endif()


