#include "Compare.h"
#include "ComparableVector.h"

#include <iterator>

using namespace std;
using namespace ttk;

Compare::Compare() : mesh1_{nullptr}, mesh2_{nullptr}
{
}

Compare::~Compare()
{
}

int Compare::computeMeshDiff(unsigned char* const vertArr, unsigned char* const cellArr)
{
   if (!vertArr) {
      cerr << "[Compare]: needs and array of nb verts char using setVertsArray" << endl;
      return -1;
   }

   if (!cellArr) {
      cerr << "[Compare]: needs and array of nb cells char using setCellsArray" << endl;
      return -1;
   }

   int returnCode = 0;

   const bool hasDiffVerts = computeMeshVertsDiff(vertArr);
   const bool hasDiffCells = computeMeshCellsDiff(cellArr);

   const idVertex nbVerts1 = mesh1_->getNumberOfVertices();
   const idVertex nbVerts2 = mesh2_->getNumberOfVertices();
   if (nbVerts1 != nbVerts2 || hasDiffVerts)
      returnCode = returnCode | 1 << 0;

   const idCell nbCells1 = mesh1_->getNumberOfCells();
   const idCell nbCells2 = mesh2_->getNumberOfCells();
   if (nbCells1 != nbCells2 || hasDiffCells)
      returnCode = returnCode | 1 << 1;

   return returnCode;
}

// Private

bool Compare::computeMeshVertsDiff(unsigned char* const vertArr)
{
   bool           hasDiffVerts = false;
   const idVertex nbVerts1     = mesh1_->getNumberOfVertices();
   const idVertex nbVerts2     = mesh2_->getNumberOfVertices();

   // init all arrays
   vertMapperM1toM2_.resize(nbVerts1);
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
   for (idVertex i = 0; i < nbVerts1; ++i) {
      vertMapperM1toM2_[i] = -1;
      vertArr[i]           = 1;
   }

   // Retains all vertices positions for mesh2
   std::map<ComparableVector<float>, idVertex> posVertsM2;
   for (idVertex v2 = 0; v2 < nbVerts2; ++v2) {
      ComparableVector<float> posV(3);
      mesh2_->getVertexPoint(v2, posV[0], posV[1], posV[2]);
      posVertsM2.emplace(posV, v2);
   }

   // match vertices from mesh1 (fill vertMapperM1toM2_)
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
   for (idVertex v1 = 0; v1 < nbVerts1; ++v1) {
      ComparableVector<float> posV(3);
      mesh1_->getVertexPoint(v1, posV[0], posV[1], posV[2]);
      if (posVertsM2.count(posV)) {
         vertMapperM1toM2_[v1] = posVertsM2[posV];
         vertArr[v1]           = 0;
      } else {
         hasDiffVerts = true;
      }
   }

   return hasDiffVerts;
}

bool Compare::computeMeshCellsDiff(unsigned char* const cellArr)
{
   if (vertMapperM1toM2_.empty()) {
      std::cerr << "[Compare] comuteCellDiff needs vertMapper to be filled" << std::endl;
      return true;
   }

   bool         hasDiffCells = false;
   const idCell nbCells1     = mesh1_->getNumberOfCells();
   const idCell nbCells2     = mesh2_->getNumberOfCells();

   // init all arrays
   cellMapperM1toM2_.resize(nbCells1);
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
   for (idCell i = 0; i < nbCells1; ++i) {
      cellMapperM1toM2_[i] = -1;
      cellArr[i]           = 1;
   }

   // for cells of mesh1, retains their list of vertices using mesh1 vertices
   std::map<ComparableVector<idVertex>, idCell> vertListCellsM2;
   for (idCell c2 = 0; c2 < nbCells2; ++c2) {
      const idVertex             nbVerts = mesh2_->getCellVertexNumber(c2);
      ComparableVector<idVertex> vertsC;
      vertsC.reserve(nbVerts);
      for (idVertex v = 0; v < nbVerts; ++v) {
         idVertex v2Id;
         mesh2_->getCellVertex(c2, v, v2Id);
         vertsC.emplace_back(v2Id);
      }
      std::sort(vertsC.begin(), vertsC.end());
      vertListCellsM2.emplace(vertsC, c2);
   }
   // match cells from mesh1 (fill cellMapperM1toM2_)
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
   for (idCell c1 = 0; c1 < nbCells1; ++c1) {
      bool                       nextCell = false;
      const idVertex             nbVerts  = mesh1_->getCellVertexNumber(c1);
      ComparableVector<idVertex> vertsC;
      vertsC.reserve(nbVerts);
      for (idVertex v = 0; v < nbVerts; ++v) {
         idVertex v1Id;
         mesh1_->getCellVertex(c1, v, v1Id);
         if (vertMapperM1toM2_[v1Id] == -1) {
            // this vertex of mesh1 doen not exists in mesh2
            nextCell = true;
            break;
         }
         vertsC.emplace_back(vertMapperM1toM2_[v1Id]);
      }
      if (nextCell) {
         hasDiffCells = true;
         continue;
      }
      std::sort(vertsC.begin(), vertsC.end());
      if (vertListCellsM2.count(vertsC)) {
         cellMapperM1toM2_[c1] = vertListCellsM2[vertsC];
         cellArr[c1]           = 0;
      } else {
         hasDiffCells = true;
      }
   }

   return hasDiffCells;
}

