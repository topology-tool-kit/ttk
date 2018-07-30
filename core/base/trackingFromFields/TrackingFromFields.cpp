#include <TrackingFromFields.h>

int ttk::TrackingFromFields::performDiagramComputation(
  int fieldNumber,
  std::vector<std::vector<diagramTuple>>& persistenceDiagrams,
  const ttk::Wrapper *wrapper)
{

  for (int i = 0; i < fieldNumber; ++i)
  {
    ttk::PersistenceDiagram persistenceDiagram_;
    persistenceDiagram_.setWrapper(wrapper);
    persistenceDiagram_.setupTriangulation(triangulation_);
      // should have been done before

    std::vector<std::tuple<ttk::dcg::Cell,ttk::dcg::Cell>> dmt_pairs;
    persistenceDiagram_.setDMTPairs(&dmt_pairs);
    persistenceDiagram_.setInputScalars(inputData_[i]);
    persistenceDiagram_.setInputOffsets(inputOffsets_);
    persistenceDiagram_.setComputeSaddleConnectors(false);
    std::vector<std::tuple<int, ftm::NodeType, int, ftm::NodeType, dataType, int>> CTDiagram;

    persistenceDiagram_.setOutputCTDiagram(&CTDiagram);
    persistenceDiagram_.execute<dataType, int>();

    // Copy diagram into augmented diagram.
    persistenceDiagrams[i] = std::vector<diagramTuple>(CTDiagram.size());

    for (int j = 0; j < (int) CTDiagram.size(); ++j) {
      float p[3];
      float q[3];
      auto currentTuple = CTDiagram[j];
      const int a = std::get<0>(currentTuple);
      const int b = std::get<2>(currentTuple);
      triangulation_->getVertexPoint(a, p[0], p[1], p[2]);
      triangulation_->getVertexPoint(b, q[0], q[1], q[2]);
      const double sa = ((double*) inputData_[i])[a];
      const double sb = ((double*) inputData_[i])[b];
      diagramTuple dt = std::make_tuple(
        std::get<0>(currentTuple),
        std::get<1>(currentTuple),
        std::get<2>(currentTuple),
        std::get<3>(currentTuple),
        std::get<4>(currentTuple),
        std::get<5>(currentTuple),
        sa, p[0], p[1], p[2],
        sb, q[0], q[1], q[2]);

      persistenceDiagrams[i][j] = dt;
    }
  }

  return 0;
}
