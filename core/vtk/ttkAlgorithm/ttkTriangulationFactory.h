#pragma once

#include <ttkAlgorithmModule.h>

#include <Debug.h>
#include <unordered_map>
#include <vtkType.h>

class vtkDataSet;
class vtkImageData;
class vtkPointSet;
class vtkPoints;
class vtkCellArray;
namespace ttk {
  class Triangulation;
}

using RegistryTriangulation = std::unique_ptr<ttk::Triangulation>;

struct RegistryValue {
  RegistryTriangulation triangulation;
  vtkDataSet *owner;

  vtkMTimeType cellModTime{0};

  int extent[6];
  double origin[3];
  double spacing[3];
  int dimensions[3];

  RegistryValue(vtkDataSet *dataSet, ttk::Triangulation *triangulation_);
  bool isValid(vtkDataSet *dataSet) const;
};

using RegistryKey = long long;
using Registry = std::unordered_map<RegistryKey, RegistryValue>;

class TTKALGORITHM_EXPORT ttkTriangulationFactory : public ttk::Debug {
public:
  static ttk::Triangulation *
    GetTriangulation(int debugLevel, float cacheRatio, vtkDataSet *object);

  static ttkTriangulationFactory Instance;
  static RegistryKey GetKey(vtkDataSet *dataSet);

#ifdef _WIN32
  // to fix a weird MSVC warning about unique_ptr inside
  // unordered_map, this dummy class member should be declared before
  // the Registry
  RegistryTriangulation dummy{};
#endif // _WIN32
  Registry registry;

private:
  RegistryTriangulation CreateImplicitTriangulation(vtkImageData *image);
  RegistryTriangulation CreateExplicitTriangulation(vtkPointSet *pointSet);
  RegistryTriangulation CreateTriangulation(vtkDataSet *dataSet);
  int FindImplicitTriangulation(ttk::Triangulation *&triangulation,
                                vtkImageData *image);

  ttkTriangulationFactory();
};
