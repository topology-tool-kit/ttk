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

struct RegistryValue {
  std::unique_ptr<ttk::Triangulation> triangulation;
  vtkDataSet *owner;

  vtkMTimeType cellModTime{0};

  int extent[6];
  double origin[3];
  double spacing[3];
  int dimensions[3];

  RegistryValue(vtkDataSet *dataSet, ttk::Triangulation *triangulation_);
  bool isValid(vtkDataSet *dataSet) const;
};

typedef long long RegistryKey;
typedef std::unordered_map<RegistryKey, RegistryValue> Registry;

class TTKALGORITHM_EXPORT ttkTriangulationFactory : public ttk::Debug {
public:
  static ttk::Triangulation *GetTriangulation(int debugLevel,
                                              vtkDataSet *object);

  static ttkTriangulationFactory Instance;
  static RegistryKey GetKey(vtkDataSet *dataSet);

  Registry registry;

private:
  std::unique_ptr<ttk::Triangulation>
    CreateImplicitTriangulation(vtkImageData *image);
  std::unique_ptr<ttk::Triangulation>
    CreateExplicitTriangulation(vtkPointSet *pointSet);
  std::unique_ptr<ttk::Triangulation> CreateTriangulation(vtkDataSet *dataSet);
  int FindImplicitTriangulation(ttk::Triangulation *&triangulation,
                                vtkImageData *image);

  ttkTriangulationFactory();
};
