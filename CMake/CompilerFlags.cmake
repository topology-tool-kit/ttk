# Compiler flags to be shared by all base, VTK and standalone targets

if (NOT MSVC) # GCC and Clang

  # warning flags
  list(APPEND TTK_COMPILER_FLAGS
    -Wall
    -Wextra
    -Wtype-limits
    -Wshadow
    )

  # performance and debug flags
  if(TTK_ENABLE_CPU_OPTIMIZATION AND CMAKE_BUILD_TYPE MATCHES Release)
    # -O3 already enabled by CMake's Release configuration
    list(APPEND TTK_COMPILER_FLAGS -march=native -Wfatal-errors)
  endif()

  if(CMAKE_BUILD_TYPE MATCHES Debug)
    list(APPEND TTK_COMPILER_FLAGS -O0 -g -pg)
  endif()

  # hardened linker flags for Clang and GCC
  if(${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
    list(APPEND TTK_LINKER_FLAGS
      -Wl,--as-needed
      -Wl,--no-undefined
      )
  endif()

elseif(${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang" AND "x${CMAKE_CXX_SIMULATE_ID}" STREQUAL "xMSVC")
  # clang-cl on Windows
  # c.f. https://stackoverflow.com/questions/50857779/cmake-detects-clang-cl-as-clang

  # warning flags
  list(APPEND TTK_COMPILER_FLAGS -Wno-unused-parameter)

else() # MSVC

  # warning flags
  list(APPEND TTK_COMPILER_FLAGS /W4)

  # disabled warnings
  list(APPEND TTK_COMPILER_FLAGS /bigobj /wd4005 /wd4061 /wd4100
  /wd4146 /wd4221 /wd4242 /wd4244 /wd4245 /wd4263 /wd4264 /wd4267
  /wd4273 /wd4275 /wd4296 /wd4305 /wd4365 /wd4371 /wd4435 /wd4456
  /wd4457 /wd4514 /wd4619 /wd4625 /wd4626 /wd4628 /wd4668 /wd4701
  /wd4702 /wd4710 /wd4800 /wd4820 /wd4996 /wd5027 /wd5029 /wd5031)

endif()
