ttk_add_base_library(mergeTreeNeuralNetwork
  SOURCES
    MergeTreeNeuralBase.cpp
    MergeTreeNeuralLayer.cpp
    MergeTreeNeuralNetwork.cpp
    MergeTreeTorchUtils.cpp
  HEADERS
    MergeTreeNeuralBase.h
    MergeTreeNeuralLayer.h
    MergeTreeNeuralNetwork.h
    MergeTreeTorchUtils.h
  DEPENDS
    mergeTreePrincipalGeodesics
    geometry
)

if(TTK_ENABLE_TORCH)
  target_include_directories(mergeTreeNeuralNetwork PUBLIC ${TORCH_INCLUDE_DIRS})
  target_compile_options(mergeTreeNeuralNetwork PUBLIC "${TORCH_CXX_FLAGS}")
  target_link_libraries(mergeTreeNeuralNetwork PUBLIC "${TORCH_LIBRARIES}")
  target_compile_definitions(mergeTreeNeuralNetwork PUBLIC TTK_ENABLE_TORCH)
endif()
