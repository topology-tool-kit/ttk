/// \ingroup vtk
/// \class ttkTopoMap
/// \author Alexandre Talon <alexandre.talon@lip6.fr>
/// \date June 2023
///
/// \brief TTK VTK-filter that wraps the ttk::TopoMap module.
///
/// This class generates a mapper from a data-set. The mapper is computed using
/// some minimum spanning tree and convex hull arguments.
///
/// \param Input vtkDataSet.
/// \param Output vtkDataSet.
///
///
/// \b Related \b publication: \n
/// "Topomap: A 0-dimensional homology preserving projection of high-dimensional
/// data"\n Harish Doraiswamy, Julien Tierny, Paulo J. S. Silva, Luis Gustavo
/// Nonato, and Claudio Silva\n Proc. of IEEE VIS 2020.\n IEEE Transactions on
/// Visualization and Computer Graphics 27(2): 561-571, 2020.

#pragma once

// VTK Module
#include <ttkTopoMapModule.h>

// ttk code includes
#include <TopoMap.h>
#include <ttkAlgorithm.h>
#include <ttkMacros.h>

// vtkTypeMacro(ttkTopoMap, ttkAlgorithm);

class TTKTOPOMAP_EXPORT ttkTopoMap : public ttkAlgorithm,
                                     protected ttk::TopoMap {
private:
  /**
   * TODO 5: Add all filter parameters only as private member variables and
   *         initialize them here.
   */
  std::string OutputArrayName{"AveragedScalarField"};

public:
  vtkSetMacro(SelectMatrixWithRegexp, bool);
  vtkGetMacro(SelectMatrixWithRegexp, bool);

  vtkSetMacro(DistanceMatrixRegexp, const std::string &);
  vtkGetMacro(DistanceMatrixRegexp, std::string);

  void SetDistMat(const std::string &s) {
    this->DistanceMatrixNames.emplace_back(s);
    this->Modified();
  }
  void ClearDistMat() {
    this->DistanceMatrixNames.clear();
    this->Modified();
  }

  vtkSetMacro(AngularSampleNb, unsigned long);
  vtkGetMacro(AngularSampleNb, unsigned long);

  vtkSetMacro(CheckMST, bool);
  vtkGetMacro(CheckMST, bool);

  static ttkTopoMap *New();
  vtkTypeMacro(ttkTopoMap, ttkAlgorithm);

protected:
  ttkTopoMap();
  ~ttkTopoMap() override = default;

  int FillInputPortInformation(int port, vtkInformation *info) override;

  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

private:
  // domain vertices coordinates in high dimension
  bool SelectMatrixWithRegexp{false};
  std::string DistanceMatrixRegexp{".*"};
  std::vector<std::string> DistanceMatrixNames{};
};
