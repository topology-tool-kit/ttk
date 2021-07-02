ttk_add_base_library(webSocketIO
  HEADERS
    WebSocketIO.h
    WebSocketIO.cpp
  DEPENDS
    triangulation
)

if (TTK_ENABLE_WEBSOCKETPP)
  target_compile_definitions(webSocketIO PUBLIC TTK_ENABLE_WEBSOCKETPP)
  target_include_directories(webSocketIO PUBLIC ${WEBSOCKETPP_INCLUDE_DIR})
endif()
