/// \ingroup vtk
/// \class ttkTopologicalSimplificationByPersistence
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 03.06.2021
///
/// \brief TTK VTK-filter that computes a persistence-based simplification of a
/// scalar field.
///
/// Given an input scalar field and a persistence threshold (either as an
/// absolute value or a fraction of the scalar range), this filter modifies the
/// scalar field such that it no longer exhibits persistence pairs below the
/// given threshold. All other pairs are unaffected. To this end the filter uses
/// the persistence-sensitive specialization of localized topological
/// simplification (PLTS). Note that this filter will also compute an
/// unambiguous global vertex order that can be used in subsequent topological
/// data analysis.

/// \b Related \b publications \n
///   "Generalized Topological Simplification of Scalar Fields on Surfaces"\n
///   Julien Tierny, Valerio Pascucci\n
///   Proc. of IEEE VIS 2012.\n
///   IEEE Transactions on Visualization and Computer Graphics, 2012.

///   "Localized Topological Simplification of Scalar Data"\n
///   Jonas Lukasczyk, Christoph Garth, Ross Maciejewski, Julien Tierny\n
///   Proc. of IEEE VIS 2020.\n
///   IEEE Transactions on Visualization and Computer Graphics
///
/// \param Input vtkDataSet.
/// \param Output vtkDataSet.
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutputDataObject()).
///
/// The input data array needs to be specified via the standard VTK call
/// vtkAlgorithm::SetInputArrayToProcess() with the following parameters:
/// \param idx 0 (FIXED: the first array the algorithm requires)
/// \param port 0 (FIXED: first port)
/// \param connection 0 (FIXED: first connection)
/// \param fieldAssociation 0 (FIXED: point data)
/// \param arrayName (DYNAMIC: string identifier of the input array)
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \sa ttkTopologicalSimplification
/// \sa ttk::LocalizedTopologicalSimplification
/// \sa ttkAlgorithm

#pragma once
#include <ttkTopologicalSimplificationByPersistenceModule.h>

// VTK Includes
#include <ttkAlgorithm.h>

// TTK Includes
#include <LocalizedTopologicalSimplification.h>
#include <ttkMacros.h>

class TTKTOPOLOGICALSIMPLIFICATIONBYPERSISTENCE_EXPORT
  ttkTopologicalSimplificationByPersistence
  : public ttkAlgorithm,
    protected ttk::lts::LocalizedTopologicalSimplification {
private:
  double PersistenceThreshold{0.0};
  bool ThresholdIsAbsolute{true};
  bool ComputePerturbation{false};
  PAIR_TYPE PairType{PAIR_TYPE::EXTREMUM_SADDLE};

public:
  vtkSetMacro(PersistenceThreshold, double);
  vtkGetMacro(PersistenceThreshold, double);
  vtkSetMacro(ThresholdIsAbsolute, bool);
  vtkGetMacro(ThresholdIsAbsolute, bool);
  vtkSetMacro(ComputePerturbation, bool);
  vtkGetMacro(ComputePerturbation, bool);
  ttkSetEnumMacro(PairType, PAIR_TYPE);
  vtkGetEnumMacro(PairType, PAIR_TYPE);

  static ttkTopologicalSimplificationByPersistence *New();
  vtkTypeMacro(ttkTopologicalSimplificationByPersistence, ttkAlgorithm);

protected:
  ttkTopologicalSimplificationByPersistence();
  ~ttkTopologicalSimplificationByPersistence() override;

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};
