#pragma once

#include <vtkIdTypeArray.h>
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

#ifndef vtkTemplate2MacroCase1
#define vtkTemplate2MacroCase1(type1N, type1, call)                            \
  vtkTemplate2MacroCase2(type1N, type1, VTK_DOUBLE, double, call);             \
  vtkTemplate2MacroCase2(type1N, type1, VTK_FLOAT, float, call);               \
  vtkTemplate2MacroCase2(type1N, type1, VTK_LONG_LONG, long long, call);       \
  vtkTemplate2MacroCase2(                                                      \
    type1N, type1, VTK_UNSIGNED_LONG_LONG, unsigned long long, call);          \
  vtkTemplate2MacroCase2(type1N, type1, VTK_ID_TYPE, vtkIdType, call);         \
  vtkTemplate2MacroCase2(type1N, type1, VTK_LONG, long, call);                 \
  vtkTemplate2MacroCase2(                                                      \
    type1N, type1, VTK_UNSIGNED_LONG, unsigned long, call);                    \
  vtkTemplate2MacroCase2(type1N, type1, VTK_INT, int, call);                   \
  vtkTemplate2MacroCase2(type1N, type1, VTK_UNSIGNED_INT, unsigned int, call); \
  vtkTemplate2MacroCase2(type1N, type1, VTK_SHORT, short, call);               \
  vtkTemplate2MacroCase2(                                                      \
    type1N, type1, VTK_UNSIGNED_SHORT, unsigned short, call);                  \
  vtkTemplate2MacroCase2(type1N, type1, VTK_CHAR, char, call);                 \
  vtkTemplate2MacroCase2(type1N, type1, VTK_SIGNED_CHAR, signed char, call);   \
  vtkTemplate2MacroCase2(type1N, type1, VTK_UNSIGNED_CHAR, unsigned char, call)
#endif

#ifndef vtkTemplate2MacroCase2
#define vtkTemplate2MacroCase2(type1N, type1, type2N, type2, call) \
  case vtkTemplate2PackMacro(type1N, type2N): {                    \
    typedef type1 VTK_T1;                                          \
    typedef type2 VTK_T2;                                          \
    call;                                                          \
  }; break
#endif
