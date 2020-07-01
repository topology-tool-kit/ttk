#pragma once

#include <vtkIntArray.h>

#define TTK_COMMA ,

#ifdef TTK_ENABLE_64BIT_IDS
using ttkSimplexIdTypeArray = vtkIdTypeArray;
#else
using ttkSimplexIdTypeArray = vtkIntArray;
#endif

#define ttkVtkTemplateMacroCase(                         \
  dataType, triangulationType, triangulationClass, call) \
  case triangulationType: {                              \
    typedef triangulationClass TTK_TT;                   \
    switch(dataType) { vtkTemplateMacro((call)); };      \
  }; break;

#define ttkVtkTemplateMacro(dataType, triangulationType, call)            \
  switch(triangulationType) {                                             \
    ttkVtkTemplateMacroCase(dataType, ttk::Triangulation::Type::EXPLICIT, \
                            ttk::ExplicitTriangulation, call);            \
    ttkVtkTemplateMacroCase(dataType, ttk::Triangulation::Type::IMPLICIT, \
                            ttk::ImplicitTriangulation, call);            \
    ttkVtkTemplateMacroCase(dataType, ttk::Triangulation::Type::PERIODIC, \
                            ttk::PeriodicImplicitTriangulation, call);    \
  }

#define ttkTemplate2IdMacro(call)                                           \
  vtkTemplate2MacroCase1(VTK_LONG_LONG, long long, call);                   \
  vtkTemplate2MacroCase1(VTK_UNSIGNED_LONG_LONG, unsigned long long, call); \
  vtkTemplate2MacroCase1(VTK_ID_TYPE, vtkIdType, call);                     \
  vtkTemplate2MacroCase1(VTK_LONG, long, call);                             \
  vtkTemplate2MacroCase1(VTK_UNSIGNED_LONG, unsigned long, call);           \
  vtkTemplate2MacroCase1(VTK_INT, int, call);                               \
  vtkTemplate2MacroCase1(VTK_UNSIGNED_INT, unsigned int, call);

#define ttkVtkTemplate2MacroCase(                                    \
  dataType0, dataType1, triangulationType, triangulationClass, call) \
  case triangulationType: {                                          \
    typedef triangulationClass TTK_TT;                               \
    switch(vtkTemplate2PackMacro(dataType0, dataType1)) {            \
      vtkTemplate2Macro((call));                                     \
    };                                                               \
  }; break;

#define ttkVtkTemplate2Macro(dataType0, dataType1, triangulationType, call) \
  switch(triangulationType) {                                               \
    ttkVtkTemplate2MacroCase(dataType0, dataType1,                          \
                             ttk::Triangulation::Type::EXPLICIT,            \
                             ttk::ExplicitTriangulation, call);             \
    ttkVtkTemplate2MacroCase(dataType0, dataType1,                          \
                             ttk::Triangulation::Type::IMPLICIT,            \
                             ttk::ImplicitTriangulation, call);             \
    ttkVtkTemplate2MacroCase(dataType0, dataType1,                          \
                             ttk::Triangulation::Type::PERIODIC,            \
                             ttk::PeriodicImplicitTriangulation, call);     \
  }
