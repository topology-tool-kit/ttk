#pragma once

#include <vtkSmartPointer.h>

#include <DataTypes.h>

/// Vertex / Node / Arc data inherit from this
/// master structure.
struct ObjectData {
   template <typename vtkArrayType>
   inline vtkSmartPointer<vtkArrayType> alloc(const char* fieldName, size_t nbElmnt)
   {
      vtkSmartPointer<vtkArrayType> arr = vtkSmartPointer<vtkArrayType>::New();
      arr->SetName(fieldName);
      arr->SetNumberOfComponents(1);
      arr->SetNumberOfTuples(nbElmnt);

#ifndef TTK_ENABLE_KAMIKAZE
      if (!arr) {
         cerr << "[ttkFTMTree] Error, unable to allocate " << fieldName
              << " the program will likely crash" << endl;
      }
#endif
      return arr;
   }
};

struct NodeData : public ObjectData {
};

struct ArcData : public ObjectData {
};

struct VertData : public ObjectData {
};
