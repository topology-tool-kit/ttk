#include <ttkPersistenceDiagramUtils.h>

#include <vtkCellData.h>
#include <vtkIntArray.h>
#include <vtkUnstructuredGrid.h>

// used by ttkBottleneckDistance & ttkTrackingFromPersistenceDiagrams
inline int augmentDiagrams(const std::vector<ttk::MatchingType> &matchings,
                           vtkUnstructuredGrid *const vtu0,
                           vtkUnstructuredGrid *const vtu1) {

  if(matchings.empty()) {
    return 0;
  }

  vtkNew<vtkIntArray> matchingIds0{};
  matchingIds0->SetName("MatchingIdentifier");
  matchingIds0->SetNumberOfComponents(1);
  matchingIds0->SetNumberOfTuples(vtu0->GetNumberOfCells());
  vtu0->GetCellData()->AddArray(matchingIds0);

  vtkNew<vtkIntArray> matchingIds1{};
  matchingIds1->SetName("MatchingIdentifier");
  matchingIds1->SetNumberOfComponents(1);
  matchingIds1->SetNumberOfTuples(vtu1->GetNumberOfCells());
  vtu1->GetCellData()->AddArray(matchingIds1);

  // Unaffected by default
  matchingIds0->Fill(-1);
  matchingIds1->Fill(-1);

  // Affect bottleneck matchings
  for(size_t i = 0; i < matchings.size(); ++i) {
    const auto &t = matchings[i];
    matchingIds0->SetTuple1(std::get<0>(t), i);
    matchingIds1->SetTuple1(std::get<1>(t), i);
  }

  return 1;
}
