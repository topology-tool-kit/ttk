#pragma once

class vtkIdTypeArray;
class vtkIntArray;

#define TTK_COMMA ,

#ifdef TTK_ENABLE_64BIT_IDS
using ttkSimplexIdTypeArray = vtkIdTypeArray;
#else
using ttkSimplexIdTypeArray = vtkIntArray;
#endif

#ifndef vtkSetEnumMacro
#define vtkSetEnumMacro(name, enumType)                                        \
  virtual void Set##name(enumType _arg) {                                      \
    vtkDebugMacro(<< this->GetClassName() << " (" << this                      \
                  << "): setting " #name " to "                                \
                  << static_cast<std::underlying_type<enumType>::type>(_arg)); \
    if(this->name != _arg) {                                                   \
      this->name = _arg;                                                       \
      this->Modified();                                                        \
    }                                                                          \
  }
#endif

#ifndef vtkGetEnumMacro
#define vtkGetEnumMacro(name, enumType)                                      \
  virtual enumType Get##name() const {                                       \
    vtkDebugMacro(<< this->GetClassName() << " (" << this << "): returning " \
                  << #name " of "                                            \
                  << static_cast<std::underlying_type<enumType>::type>(      \
                       this->name));                                         \
    return this->name;                                                       \
  }
#endif

#define ttkSetEnumMacro(name, enumType)                   \
  virtual void Set##name(int _arg) {                      \
    vtkDebugMacro(<< this->GetClassName() << " (" << this \
                  << "): setting " #name " to " << _arg); \
    if(this->name != static_cast<enumType>(_arg)) {       \
      this->name = static_cast<enumType>(_arg);           \
      this->Modified();                                   \
    }                                                     \
  }                                                       \
  vtkSetEnumMacro(name, enumType);

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
    ttkVtkTemplateMacroCase(dataType, ttk::Triangulation::Type::COMPACT,  \
                            ttk::CompactTriangulation, call);             \
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

// -----------------------------------------------------------------------------

#define ttkTypeMacroErrorCase(idx, type)                          \
  default: {                                                      \
    this->printErr("Unsupported " #idx "-th Template Data Type: " \
                   + std::to_string(static_cast<int>(type)));     \
  } break;

#define ttkTypeMacroCase(enum, type, number, call) \
  case enum: {                                     \
    typedef type T##number;                        \
    call;                                          \
  } break;

#define ttkTypeMacroT(group, call)                                 \
  switch(group) {                                                  \
    ttkTypeMacroCase(ttk::Triangulation::Type::EXPLICIT,           \
                     ttk::ExplicitTriangulation, 0, call);         \
    ttkTypeMacroCase(ttk::Triangulation::Type::IMPLICIT,           \
                     ttk::ImplicitTriangulation, 0, call);         \
    ttkTypeMacroCase(ttk::Triangulation::Type::PERIODIC,           \
                     ttk::PeriodicImplicitTriangulation, 0, call); \
    ttkTypeMacroErrorCase(0, group);                               \
  }

#define ttkTypeMacroR(group, call)                 \
  switch(group) {                                  \
    ttkTypeMacroCase(VTK_FLOAT, float, 0, call);   \
    ttkTypeMacroCase(VTK_DOUBLE, double, 0, call); \
    ttkTypeMacroErrorCase(0, group);               \
  }

#define ttkTypeMacroI(group, call)                                         \
  switch(group) {                                                          \
    ttkTypeMacroCase(VTK_INT, int, 0, call);                               \
    ttkTypeMacroCase(VTK_UNSIGNED_INT, unsigned int, 0, call);             \
    ttkTypeMacroCase(VTK_CHAR, char, 0, call);                             \
    ttkTypeMacroCase(VTK_SIGNED_CHAR, signed char, 0, call);               \
    ttkTypeMacroCase(VTK_UNSIGNED_CHAR, unsigned char, 0, call);           \
    ttkTypeMacroCase(VTK_LONG, long, 0, call);                             \
    ttkTypeMacroCase(VTK_LONG_LONG, long long, 0, call);                   \
    ttkTypeMacroCase(VTK_UNSIGNED_LONG, unsigned long, 0, call);           \
    ttkTypeMacroCase(VTK_UNSIGNED_LONG_LONG, unsigned long long, 0, call); \
    ttkTypeMacroCase(VTK_ID_TYPE, vtkIdType, 0, call);                     \
    ttkTypeMacroErrorCase(0, group);                                       \
  }

#define ttkTypeMacroA(group, call)                                         \
  switch(group) {                                                          \
    ttkTypeMacroCase(VTK_FLOAT, float, 0, call);                           \
    ttkTypeMacroCase(VTK_DOUBLE, double, 0, call);                         \
    ttkTypeMacroCase(VTK_INT, int, 0, call);                               \
    ttkTypeMacroCase(VTK_UNSIGNED_INT, unsigned int, 0, call);             \
    ttkTypeMacroCase(VTK_CHAR, char, 0, call);                             \
    ttkTypeMacroCase(VTK_SIGNED_CHAR, signed char, 0, call);               \
    ttkTypeMacroCase(VTK_UNSIGNED_CHAR, unsigned char, 0, call);           \
    ttkTypeMacroCase(VTK_LONG, long, 0, call);                             \
    ttkTypeMacroCase(VTK_LONG_LONG, long long, 0, call);                   \
    ttkTypeMacroCase(VTK_UNSIGNED_LONG, unsigned long, 0, call);           \
    ttkTypeMacroCase(VTK_UNSIGNED_LONG_LONG, unsigned long long, 0, call); \
    ttkTypeMacroCase(VTK_ID_TYPE, vtkIdType, 0, call);                     \
    ttkTypeMacroErrorCase(0, group);                                       \
  }

#define ttkTypeMacroAT(group0, group1, call)                \
  switch(group1) {                                          \
    ttkTypeMacroCase(ttk::Triangulation::Type::EXPLICIT,    \
                     ttk::ExplicitTriangulation, 1,         \
                     ttkTypeMacroA(group0, call));          \
    ttkTypeMacroCase(ttk::Triangulation::Type::IMPLICIT,    \
                     ttk::ImplicitTriangulation, 1,         \
                     ttkTypeMacroA(group0, call));          \
    ttkTypeMacroCase(ttk::Triangulation::Type::PERIODIC,    \
                     ttk::PeriodicImplicitTriangulation, 1, \
                     ttkTypeMacroA(group0, call));          \
    ttkTypeMacroErrorCase(1, group1);                       \
  }

#define ttkTypeMacroRT(group0, group1, call)                \
  switch(group1) {                                          \
    ttkTypeMacroCase(ttk::Triangulation::Type::EXPLICIT,    \
                     ttk::ExplicitTriangulation, 1,         \
                     ttkTypeMacroR(group0, call));          \
    ttkTypeMacroCase(ttk::Triangulation::Type::IMPLICIT,    \
                     ttk::ImplicitTriangulation, 1,         \
                     ttkTypeMacroR(group0, call));          \
    ttkTypeMacroCase(ttk::Triangulation::Type::PERIODIC,    \
                     ttk::PeriodicImplicitTriangulation, 1, \
                     ttkTypeMacroR(group0, call));          \
    ttkTypeMacroErrorCase(1, group1);                       \
  }

#define ttkTypeMacroIT(group0, group1, call)                \
  switch(group1) {                                          \
    ttkTypeMacroCase(ttk::Triangulation::Type::EXPLICIT,    \
                     ttk::ExplicitTriangulation, 1,         \
                     ttkTypeMacroI(group0, call));          \
    ttkTypeMacroCase(ttk::Triangulation::Type::IMPLICIT,    \
                     ttk::ImplicitTriangulation, 1,         \
                     ttkTypeMacroI(group0, call));          \
    ttkTypeMacroCase(ttk::Triangulation::Type::PERIODIC,    \
                     ttk::PeriodicImplicitTriangulation, 1, \
                     ttkTypeMacroI(group0, call));          \
    ttkTypeMacroErrorCase(1, group1);                       \
  }

#define ttkTypeMacroAI(group0, group1, call)                                  \
  switch(group1) {                                                            \
    ttkTypeMacroCase(VTK_INT, int, 1, ttkTypeMacroA(group0, call));           \
    ttkTypeMacroCase(                                                         \
      VTK_LONG_LONG, long long, 1, ttkTypeMacroA(group0, call));              \
    ttkTypeMacroCase(VTK_ID_TYPE, vtkIdType, 1, ttkTypeMacroA(group0, call)); \
    ttkTypeMacroErrorCase(1, group1);                                         \
  }

#define ttkTypeMacroRR(group0, group1, call)                              \
  switch(group1) {                                                        \
    ttkTypeMacroCase(VTK_FLOAT, float, 1, ttkTypeMacroR(group0, call));   \
    ttkTypeMacroCase(VTK_DOUBLE, double, 1, ttkTypeMacroR(group0, call)); \
    ttkTypeMacroErrorCase(1, group1);                                     \
  }

#define ttkTypeMacroAA(group0, group1, call)                                  \
  switch(group1) {                                                            \
    ttkTypeMacroCase(VTK_FLOAT, float, 1, ttkTypeMacroA(group0, call));       \
    ttkTypeMacroCase(VTK_DOUBLE, double, 1, ttkTypeMacroA(group0, call));     \
    ttkTypeMacroCase(VTK_INT, int, 1, ttkTypeMacroA(group0, call));           \
    ttkTypeMacroCase(                                                         \
      VTK_UNSIGNED_INT, unsigned int, 1, ttkTypeMacroA(group0, call));        \
    ttkTypeMacroCase(VTK_CHAR, char, 1, ttkTypeMacroA(group0, call));         \
    ttkTypeMacroCase(                                                         \
      VTK_SIGNED_CHAR, signed char, 1, ttkTypeMacroA(group0, call));          \
    ttkTypeMacroCase(                                                         \
      VTK_UNSIGNED_CHAR, unsigned char, 1, ttkTypeMacroA(group0, call));      \
    ttkTypeMacroCase(VTK_LONG, long, 1, ttkTypeMacroA(group0, call));         \
    ttkTypeMacroCase(                                                         \
      VTK_LONG_LONG, long long, 1, ttkTypeMacroA(group0, call));              \
    ttkTypeMacroCase(                                                         \
      VTK_UNSIGNED_LONG, unsigned long, 1, ttkTypeMacroA(group0, call));      \
    ttkTypeMacroCase(VTK_UNSIGNED_LONG_LONG, unsigned long long, 1,           \
                     ttkTypeMacroA(group0, call));                            \
    ttkTypeMacroCase(VTK_ID_TYPE, vtkIdType, 1, ttkTypeMacroA(group0, call)); \
    ttkTypeMacroErrorCase(1, group1);                                         \
  }

#define ttkTypeMacroAAA(group0, group1, group2, call)                          \
  switch(group2) {                                                             \
    ttkTypeMacroCase(                                                          \
      VTK_FLOAT, float, 2, ttkTypeMacroAA(group0, group1, call));              \
    ttkTypeMacroCase(                                                          \
      VTK_DOUBLE, double, 2, ttkTypeMacroAA(group0, group1, call));            \
    ttkTypeMacroCase(VTK_INT, int, 2, ttkTypeMacroAA(group0, group1, call));   \
    ttkTypeMacroCase(VTK_UNSIGNED_INT, unsigned int, 2,                        \
                     ttkTypeMacroAA(group0, group1, call));                    \
    ttkTypeMacroCase(VTK_CHAR, char, 2, ttkTypeMacroAA(group0, group1, call)); \
    ttkTypeMacroCase(                                                          \
      VTK_SIGNED_CHAR, signed char, 2, ttkTypeMacroAA(group0, group1, call));  \
    ttkTypeMacroCase(VTK_UNSIGNED_CHAR, unsigned char, 2,                      \
                     ttkTypeMacroAA(group0, group1, call));                    \
    ttkTypeMacroCase(VTK_LONG, long, 2, ttkTypeMacroAA(group0, group1, call)); \
    ttkTypeMacroCase(                                                          \
      VTK_LONG_LONG, long long, 2, ttkTypeMacroAA(group0, group1, call));      \
    ttkTypeMacroCase(VTK_UNSIGNED_LONG, unsigned long, 2,                      \
                     ttkTypeMacroAA(group0, group1, call));                    \
    ttkTypeMacroCase(VTK_UNSIGNED_LONG_LONG, unsigned long long, 2,            \
                     ttkTypeMacroAA(group0, group1, call));                    \
    ttkTypeMacroCase(                                                          \
      VTK_ID_TYPE, vtkIdType, 2, ttkTypeMacroAA(group0, group1, call));        \
    ttkTypeMacroErrorCase(2, group2);                                          \
  }

#define ttkTypeMacroAII(group0, group1, group2, call)                        \
  switch(group2) {                                                           \
    ttkTypeMacroCase(VTK_INT, int, 2, ttkTypeMacroAI(group0, group1, call)); \
    ttkTypeMacroCase(                                                        \
      VTK_LONG_LONG, long long, 2, ttkTypeMacroAI(group0, group1, call));    \
    ttkTypeMacroCase(                                                        \
      VTK_ID_TYPE, vtkIdType, 2, ttkTypeMacroAI(group0, group1, call));      \
    ttkTypeMacroErrorCase(2, group2);                                        \
  }

#define ttkTypeMacroRRR(group0, group1, group2, call)               \
  switch(group2) {                                                  \
    ttkTypeMacroCase(                                               \
      VTK_FLOAT, float, 2, ttkTypeMacroRR(group0, group1, call));   \
    ttkTypeMacroCase(                                               \
      VTK_DOUBLE, double, 2, ttkTypeMacroRR(group0, group1, call)); \
    ttkTypeMacroErrorCase(2, group2);                               \
  }

#define ttkTypeMacroRRI(group0, group1, group2, call)                        \
  switch(group2) {                                                           \
    ttkTypeMacroCase(VTK_INT, int, 2, ttkTypeMacroRR(group0, group1, call)); \
    ttkTypeMacroCase(                                                        \
      VTK_LONG_LONG, long long, 2, ttkTypeMacroRR(group0, group1, call));    \
    ttkTypeMacroCase(                                                        \
      VTK_ID_TYPE, vtkIdType, 2, ttkTypeMacroRR(group0, group1, call));      \
    ttkTypeMacroErrorCase(2, group2);                                        \
  }
