#include <AbstractTriangulation.h>
#include <ttkTriangulationWriter.h>

#include <array>

#include <vtkDataArray.h>
#include <vtkExecutive.h>
#include <vtkImageData.h>
#include <vtkInformation.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkUnstructuredGrid.h>

vtkStandardNewMacro(ttkTriangulationWriter);

ttkTriangulationWriter::ttkTriangulationWriter() {
  this->SetNumberOfInputPorts(1);
  this->setDebugMsgPrefix("TriangulationWriter");
}

int ttkTriangulationWriter::FillInputPortInformation(int port,
                                                     vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  }
  return 0;
}

int ttkTriangulationWriter::OpenFile() {

  std::ofstream f(Filename);

  if(!f.fail()) {
    Stream = std::move(f);
  } else {
    return -1;
  }

  return 0;
}

template <typename T>
void writeBin(std::ofstream &stream, const T var) {
  stream.write(reinterpret_cast<const char *>(&var), sizeof(var));
}

int ttkTriangulationWriter::Write() {

  const auto dataset = vtkDataSet::SafeDownCast(this->GetInput());
  const auto triangulation = ttkAlgorithm::GetTriangulation(dataset);
  if(triangulation == nullptr) {
    this->printErr("Invalid triangulation");
    return 1;
  }

  const auto nVerts = triangulation->getNumberOfVertices();
  const auto nEdges = triangulation->getNumberOfEdges();
  const auto nTriangles = triangulation->getNumberOfTriangles();
  const auto nTetras = triangulation->getNumberOfCells();
  const auto dim = triangulation->getDimensionality();

  std::vector<std::vector<std::string>> rows{
    {"Dimension", std::to_string(dim)},
    {"#Vertices", std::to_string(nVerts)},
    {"#Edges", std::to_string(nEdges)},
    {"#Triangles", std::to_string(nTriangles)},
  };
  if(dim == 3) {
    rows.emplace_back(std::vector<std::string>{
      "#Tetrahedrons", std::to_string(dim == 3 ? nTetras : 0)});
  }
  this->printMsg(rows);

#define TRI_DUMP(FUNC, N_ITER)                          \
  Stream << #FUNC << '\n';                              \
  for(ttk::SimplexId j = 0; j < (N_ITER); ++j) {        \
    Stream << #FUNC << " " << j << " ";                 \
    const auto n = triangulation->get##FUNC##Number(j); \
    if(n < 1) {                                         \
      break;                                            \
    }                                                   \
    std::vector<ttk::SimplexId> vec(n);                 \
    for(ttk::SimplexId k = 0; k < n; ++k) {             \
      triangulation->get##FUNC(j, k, vec[k]);           \
    }                                                   \
    std::sort(vec.begin(), vec.end());                  \
    for(const auto v : vec) {                           \
      Stream << v << ", ";                              \
    }                                                   \
    Stream << '\n';                                     \
  }

#define TRI_DUMP_FLAT(FUNC, N_ITER)              \
  Stream << #FUNC << '\n';                       \
  for(ttk::SimplexId j = 0; j < (N_ITER); ++j) { \
    Stream << triangulation->FUNC(j) << " ";     \
  }                                              \
  Stream << '\n';

  // vertices
  if(triangulation->hasPreconditionedVertexNeighbors()) {
    TRI_DUMP(VertexNeighbor, nVerts);
  }
  if(triangulation->hasPreconditionedVertexEdges()) {
    TRI_DUMP(VertexEdge, nVerts);
  }
  if(triangulation->hasPreconditionedVertexTriangles()) {
    TRI_DUMP(VertexTriangle, nVerts);
  }
  if(triangulation->hasPreconditionedVertexStars()) {
    TRI_DUMP(VertexStar, nVerts);
  }
  if(triangulation->hasPreconditionedVertexLinks()) {
    TRI_DUMP(VertexLink, nVerts);
  }

  // edges
  if(triangulation->hasPreconditionedEdges()) {
    TRI_DUMP(EdgeVertex, nEdges);
  }
  if(triangulation->hasPreconditionedEdgeTriangles()) {
    TRI_DUMP(EdgeTriangle, nEdges);
  }
  if(triangulation->hasPreconditionedEdgeStars()) {
    TRI_DUMP(EdgeStar, nEdges);
  }
  if(triangulation->hasPreconditionedEdgeLinks()) {
    TRI_DUMP(EdgeLink, nEdges);
  }

  // triangles
  if(triangulation->hasPreconditionedTriangles()) {
    TRI_DUMP(TriangleVertex, nTriangles);
  }
  if(triangulation->hasPreconditionedTriangleEdges()) {
    TRI_DUMP(TriangleEdge, nTriangles);
  }
  if(triangulation->hasPreconditionedTriangleStars()) {
    TRI_DUMP(TriangleStar, nTriangles);
  }
  if(triangulation->hasPreconditionedTriangleLinks()) {
    TRI_DUMP(TriangleLink, nTriangles);
  }

  // cells
  TRI_DUMP(CellVertex, nTetras);
  if(triangulation->hasPreconditionedCellEdges()) {
    TRI_DUMP(CellEdge, nTetras);
  }
  if(triangulation->hasPreconditionedCellTriangles()) {
    TRI_DUMP(CellTriangle, nTetras);
  }
  if(triangulation->hasPreconditionedCellNeighbors()) {
    TRI_DUMP(CellNeighbor, nTetras);
  }

  // boundary
  if(triangulation->hasPreconditionedBoundaryVertices()) {
    TRI_DUMP_FLAT(isVertexOnBoundary, nVerts);
  }
  if(triangulation->hasPreconditionedBoundaryEdges()) {
    TRI_DUMP_FLAT(isEdgeOnBoundary, nEdges);
  }
  if(triangulation->hasPreconditionedBoundaryTriangles()) {
    TRI_DUMP_FLAT(isTriangleOnBoundary, nTriangles);
  }

  return 0;
}

vtkDataObject *ttkTriangulationWriter::GetInput() {
  // copied from ParaView's vtkWriter::GetInput()
  if(this->GetNumberOfInputConnections(0) < 1) {
    return nullptr;
  }
  return this->GetExecutive()->GetInputData(0, 0);
}

void ttkTriangulationWriter::SetInputData(vtkDataObject *input) {
  // copied from ParaView's vtkWriter::SetInputData()
  this->SetInputDataInternal(0, input);
}
