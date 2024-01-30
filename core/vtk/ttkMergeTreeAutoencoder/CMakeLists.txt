ttk_add_vtk_module()

if(TTK_ENABLE_TORCH)
  target_link_libraries(ttkMergeTreeAutoencoder PRIVATE "${TORCH_LIBRARIES}")
  target_compile_definitions(ttkMergeTreeAutoencoder PUBLIC TTK_ENABLE_TORCH)
endif()
