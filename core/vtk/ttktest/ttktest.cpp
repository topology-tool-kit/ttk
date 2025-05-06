#include <ttktest.h>

#include <vtkInformation.h>

#include <vtkUnstructuredGrid.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkSmartPointer.h>


#include <ttkMacros.h>
#include <ttkUtils.h>


vtkStandardNewMacro(ttktest);


ttktest::ttktest() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}


int ttktest::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkUnstructuredGrid");
    return 1;
  }
  return 0;
}


int ttktest::FillOutputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
}


int ttktest::RequestData(vtkInformation *ttkNotUsed(request),
                               vtkInformationVector **inputVector,
                               vtkInformationVector *outputVector) {


// Récupérer l'entrée VTU
  vtkUnstructuredGrid* inputUG = vtkUnstructuredGrid::SafeDownCast(
    this->GetInputDataObject(0, 0)
  );
  if(!inputUG) {
    vtkErrorMacro("Entrée invalide: attendu vtkUnstructuredGrid.");
    return 0;
  }

  // Récupérer l'array ComponentLength dans CellData
  vtkDataArray* compLenArr = inputUG->GetCellData()->GetArray("ComponentLength");
  if(!compLenArr) {
    this->printErr("Array 'ComponentLength' introuvable dans CellData.");
    return 0;
  }

  vtkIdType numCells = inputUG->GetNumberOfCells();
  this->printMsg("---- Cellules avec ComponentLength == 299 ----");
  for(vtkIdType cid = 0; cid < numCells; ++cid) {
    if(static_cast<long>(compLenArr->GetTuple1(cid)) != 299) continue;

    // 1) Informations de la cellule
    unsigned char cellType = inputUG->GetCellType(cid);
    this->printMsg("Cellule [" + std::to_string(cid) + "] type=" + std::to_string(cellType)
                   + " ComponentLength=299");

    // 2) Afficher toutes les CellData arrays pour cette cellule
    int nCDArrays = inputUG->GetCellData()->GetNumberOfArrays();
    for(int ai = 0; ai < nCDArrays; ++ai) {
      vtkDataArray* cArr = inputUG->GetCellData()->GetArray(ai);
      std::string name = cArr->GetName();
      double value = cArr->GetTuple1(cid);
      this->printMsg("  CellData[" + name + "] = " + std::to_string(value));
    }

    // 3) Points associés et leurs informations
    vtkSmartPointer<vtkIdList> pts = vtkSmartPointer<vtkIdList>::New();
    inputUG->GetCellPoints(cid, pts);
    vtkIdType nPts = pts->GetNumberOfIds();
    this->printMsg("  Sommets de la cellule (" + std::to_string(nPts) + "):");
    for(vtkIdType k = 0; k < nPts; ++k) {
      vtkIdType pid = pts->GetId(k);
      double coords[3];
      inputUG->GetPoint(pid, coords);
      this->printMsg("    Point " + std::to_string(pid)
                     + " coords=("
                     + std::to_string(coords[0]) + ","
                     + std::to_string(coords[1]) + ","
                     + std::to_string(coords[2]) + ")");

      // Afficher toutes les PointData arrays pour ce sommet
      int nPDArrays = inputUG->GetPointData()->GetNumberOfArrays();
      for(int pj = 0; pj < nPDArrays; ++pj) {
        vtkDataArray* pArr = inputUG->GetPointData()->GetArray(pj);
        std::string pname = pArr->GetName();
        double pval = pArr->GetTuple1(pid);
        this->printMsg("      PointData[" + pname + "] = " + std::to_string(pval));
      }
    }
  }

  // Préparer la sortie: copie superficielle de l'entrée
  vtkUnstructuredGrid* outputUG = vtkUnstructuredGrid::SafeDownCast(
    this->GetOutputDataObject(0)
  );
  outputUG->ShallowCopy(inputUG);

  return 1;}
