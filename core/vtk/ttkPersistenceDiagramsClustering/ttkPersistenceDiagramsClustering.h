#ifndef _TTK_PERSISTENCEDIAGRAMSCLUSTERING_H
#define _TTK_PERSISTENCEDIAGRAMSCLUSTERING_H


#define BLocalMax ttk::CriticalType::Local_maximum
#define BLocalMin ttk::CriticalType::Local_minimum
#define BSaddle1  ttk::CriticalType::Saddle1
#define BSaddle2  ttk::CriticalType::Saddle2

// VTK includes -- to adapt
#include                  <vtkCharArray.h>
#include                  <vtkDataArray.h>
#include                  <vtkDataSet.h>
#include                  <vtkDataSetAlgorithm.h>
#include                  <vtkDoubleArray.h>
#include                  <vtkFiltersCoreModule.h>
#include                  <vtkFloatArray.h>
#include                  <vtkInformation.h>
#include                  <vtkInformationVector.h>
#include                  <vtkIntArray.h>
#include                  <vtkObjectFactory.h>
#include                  <vtkPointData.h>
#include                  <vtkCellData.h>
#include                  <vtkSmartPointer.h>
#include                  <vtkUnstructuredGrid.h>
#include                  <vtkTable.h>

// ttk code includes
#include                  <PersistenceDiagramsClustering.h>
#include                  <Wrapper.h>
#include                  <ttkWrapper.h>

// in this example, this wrapper takes a data-set on the input and produces a
// data-set on the output - to adapt.
// see the documentation of the vtkAlgorithm class to decide from which VTK
// class your wrapper should inherit.
#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkPersistenceDiagramsClustering
#else
class ttkPersistenceDiagramsClustering
#endif
  : public vtkDataSetAlgorithm, public ttk::Wrapper{

  public:
     void setNumberOfInputsFromCommandLine(int number){
        numberOfInputsFromCommandLine = number;
        SetNumberOfInputPorts(number);
        }
    static ttkPersistenceDiagramsClustering* New();

    vtkTypeMacro(ttkPersistenceDiagramsClustering, vtkDataSetAlgorithm);

    // default ttk setters
    vtkSetMacro(debugLevel_, int);

    void SetThreads(){
      if(!UseAllCores)
        threadNumber_ = ThreadNumber;
      else{
        threadNumber_ = ttk::OsCall::getNumberOfCores();
      }
      Modified();
    }

    /*void SetThreadNumber(int threadNumber){
      ThreadNumber = threadNumber;
      SetThreads();
    }*/

    void SetUseAllCores(bool onOff){
      UseAllCores = onOff;
      SetThreads();
    }
    // end of default ttk setters

    // set-getters macros to define from each variable you want to access from
    // the outside (in particular from paraview) - to adapt.

    vtkSetMacro(ScalarField, std::string);
    vtkGetMacro(ScalarField, std::string);

	vtkSetMacro(WassersteinMetric, std::string);
    vtkGetMacro(WassersteinMetric, std::string);

	vtkSetMacro(UseProgressive, int);
    vtkGetMacro(UseProgressive, int);

    vtkSetMacro(TimeLimit, double);
    vtkGetMacro(TimeLimit, double);

	vtkSetMacro(UseOutputMatching, int);
    vtkGetMacro(UseOutputMatching, int);

	vtkSetMacro(ThreadNumber, int);
    vtkGetMacro(ThreadNumber, int);

	vtkSetMacro(Alpha, double);
    vtkGetMacro(Alpha, double);

  vtkSetMacro(Lambda, double);
    vtkGetMacro(Lambda, double);

	vtkSetMacro(NumberOfClusters, int);
	vtkGetMacro(NumberOfClusters, int);

	vtkSetMacro(UseAccelerated, bool);
	vtkGetMacro(UseAccelerated, bool);

	vtkSetMacro(UseKmeansppInit, bool);
	vtkGetMacro(UseKmeansppInit, bool);

  vtkSetMacro(Deterministic, bool);
  vtkGetMacro(Deterministic, bool);

  vtkSetMacro(PairTypeClustering, int);
  vtkGetMacro(PairTypeClustering, int);
  protected:

    ttkPersistenceDiagramsClustering();

    ~ttkPersistenceDiagramsClustering();

	template <typename dataType>
    double getPersistenceDiagram(
      std::vector<diagramTuple>* diagram,
      vtkUnstructuredGrid *CTPersistenceDiagram_,
      const double spacing,
      const int diagramNumber);

	int FillInputPortInformation(int port, vtkInformation *info);
    int FillOutputPortInformation(int port, vtkInformation *info);
    
	template<typename dataType>
    vtkSmartPointer<vtkUnstructuredGrid> createOutputClusteredDiagrams(
        std::vector<std::vector<diagramTuple>>& all_CTDiagrams, std::vector<int> inv_clustering, double max_dimension);

	template<typename dataType>
    vtkSmartPointer<vtkUnstructuredGrid> createOutputCentroids(
        std::vector<std::vector<diagramTuple>>* final_centroids, std::vector<int> inv_clustering, double max_dimension);

    int RequestData(vtkInformation *request,
      vtkInformationVector **inputVector, vtkInformationVector *outputVector);


  private:
    int                 numberOfInputsFromCommandLine;
    int                   PairTypeClustering;
    bool                  Deterministic;
    bool                  UseAllCores;
    int                   ThreadNumber;
	bool                  UseOutputMatching;
	double                Alpha;
  double                Lambda;

	int 				  NumberOfClusters;
	bool 				  UseAccelerated;
	bool 				  UseKmeansppInit;

    std::string                ScalarField;
	std::string                WassersteinMetric;

	bool	              UseProgressive;
	double	              TimeLimit;

    // base code features
    int doIt(vtkDataSet **input,
            vtkUnstructuredGrid *outputClusters,
            vtkUnstructuredGrid *outputCentroids,
             int numInputs);

    bool needsToAbort();

    int updateProgress(const float &progress);

};

template <typename dataType>
double ttkPersistenceDiagramsClustering::getPersistenceDiagram(
  std::vector<diagramTuple>* diagram,
  vtkUnstructuredGrid *CTPersistenceDiagram_,
  const double spacing,
  const int diagramNumber)
{
  vtkIntArray* vertexIdentifierScalars =
    vtkIntArray::SafeDownCast(CTPersistenceDiagram_->
      GetPointData()->GetArray(ttk::VertexScalarFieldName));

  vtkIntArray* nodeTypeScalars =
    vtkIntArray::SafeDownCast(CTPersistenceDiagram_->
      GetPointData()->GetArray("CriticalType"));

  vtkIntArray* pairIdentifierScalars =
    vtkIntArray::SafeDownCast(CTPersistenceDiagram_->
      GetCellData()->GetArray("PairIdentifier"));

  vtkIntArray* extremumIndexScalars =
    vtkIntArray::SafeDownCast(CTPersistenceDiagram_->
      GetCellData()->GetArray("PairType"));

  vtkDoubleArray* persistenceScalars =
    vtkDoubleArray::SafeDownCast(CTPersistenceDiagram_->
      GetCellData()->GetArray("Persistence"));

  vtkDoubleArray* birthScalars =
    vtkDoubleArray::SafeDownCast(CTPersistenceDiagram_->
      GetPointData()->GetArray("Birth"));
      vtkFloatArray* critCoordinates =
        vtkFloatArray::SafeDownCast(CTPersistenceDiagram_->
          GetPointData()->GetArray("Coordinates"));
  vtkDoubleArray* deathScalars =
    vtkDoubleArray::SafeDownCast(CTPersistenceDiagram_->
      GetPointData()->GetArray("Death"));

  vtkPoints* points = (CTPersistenceDiagram_->GetPoints());
  int pairingsSize = (int) pairIdentifierScalars->GetNumberOfTuples();
  // FIX : no more missed pairs
   for(int pair_index = 0; pair_index<pairingsSize; pair_index++){
    const float index_of_pair = pair_index;
    if(*pairIdentifierScalars->GetTuple(pair_index)!=-1)
        pairIdentifierScalars->SetTuple(pair_index, &index_of_pair);
  }
  //auto s = (float) 0.0;

  if (!deathScalars != !birthScalars) return -2;
  //bool Is3D = !(!deathScalars && !birthScalars);
  //if (!Is3D && diagramNumber == 1) s = (float) spacing;

  if (pairingsSize < 1 || !vertexIdentifierScalars
      || !pairIdentifierScalars || !nodeTypeScalars
      || !persistenceScalars || !extremumIndexScalars || !points)
    return -2;

  diagram->resize(pairingsSize+1);
  int nbNonCompact = 0;
  double max_dimension=0;

  for (int i = 0; i < pairingsSize; ++i) {

    int vertexId1 = vertexIdentifierScalars->GetValue(2*i);
    int vertexId2 = vertexIdentifierScalars->GetValue(2*i+1);
    int nodeType1 = nodeTypeScalars->GetValue(2*i);
    int nodeType2 = nodeTypeScalars->GetValue(2*i+1);

    int pairIdentifier = pairIdentifierScalars->GetValue(i);
    int pairType = extremumIndexScalars->GetValue(i);
    double persistence = persistenceScalars->GetValue(i);

    double* critCoords1 = critCoordinates->GetTuple3(2*i);

    auto coordX1 = (float) critCoords1[0];
    auto coordY1 = (float) critCoords1[1];
    auto coordZ1 = (float) critCoords1[2];

    double* critCoords2 = critCoordinates->GetTuple3(2*i+1);
    auto coordX2 = (float) critCoords2[0];
    auto coordY2 = (float) critCoords2[1];
    auto coordZ2 = (float) critCoords2[2];
    int index1 = 2*i;
    double* coords1 = points->GetPoint(index1);
    auto x1 = (float) coords1[0];
    //auto y1 = (float) coords1[1];
    //auto z1 = (float) coords1[2];

    int index2 = index1 + 1;
    double* coords2 = points->GetPoint(index2);
    //auto x2 = (float) coords2[0];
    auto y2 = (float) coords2[1];
    //auto z2 = (float) coords2[2];

    dataType value1 = (!birthScalars) ? (dataType) x1 :
                      (dataType) birthScalars->GetValue(2*i);
    dataType value2 = (!deathScalars) ? (dataType) y2 :
                      (dataType) deathScalars->GetValue(2*i+1);
    
    if(value1 > max_dimension)  max_dimension = value1;
    if(value2 > max_dimension)  max_dimension = value2;

    if (pairIdentifier != -1 && pairIdentifier < pairingsSize){
        if(pairIdentifier ==0){
            diagram->at(0) = std::make_tuple(
                vertexId1, (BNodeType) 1,
                vertexId2, (BNodeType) 3,
                (dataType) persistence,
                pairType,
                value1, coordX1, coordY1, coordZ1,
                value2, coordX2, coordY2, coordZ2
            );

            diagram->at(pairingsSize) = std::make_tuple(
                vertexId1, (BNodeType) 0,
                vertexId2, (BNodeType) 1,
                (dataType) persistence,
                pairType,
                value1, coordX1, coordY1, coordZ1,
                value2, coordX2, coordY2, coordZ2
            );

        }
        else
        {
            diagram->at(pairIdentifier) = std::make_tuple(
                vertexId1, (BNodeType) nodeType1,
                vertexId2, (BNodeType) nodeType2,
                (dataType) persistence,
                pairType,
                value1, coordX1, coordY1, coordZ1,
                value2, coordX2, coordY2, coordZ2
            );
        }
    }
    if (pairIdentifier >= pairingsSize) {
      nbNonCompact++;
      if (nbNonCompact == 0) {
        std::stringstream msg;
        msg << "[TTKPersistenceDiagramClustering] Diagram pair identifiers "
            << "must be compact (not exceed the diagram size). "
            << std::endl;
        dMsg(std::cout, msg.str(), timeMsg);
      }
    }
  }

  if (nbNonCompact > 0) {
    {
      std::stringstream msg;
      msg << "[TTKPersistenceDiagramClustering] Missed "
          << nbNonCompact << " pairs due to non-compactness."
          << std::endl;
      dMsg(std::cout, msg.str(), timeMsg);
    }
  }


  return max_dimension;
}
    
    
// {
//   vtkIntArray* vertexIdentifierScalars =
//     vtkIntArray::SafeDownCast(CTPersistenceDiagram_->
//       GetPointData()->GetArray(ttk::VertexScalarFieldName));

//   vtkIntArray* nodeTypeScalars =
//     vtkIntArray::SafeDownCast(CTPersistenceDiagram_->
//       GetPointData()->GetArray("CriticalType"));

//   vtkIntArray* pairIdentifierScalars =
//     vtkIntArray::SafeDownCast(CTPersistenceDiagram_->
//       GetCellData()->GetArray("PairIdentifier"));

//   vtkIntArray* extremumIndexScalars =
//     vtkIntArray::SafeDownCast(CTPersistenceDiagram_->
//       GetCellData()->GetArray("PairType"));

//   vtkDoubleArray* persistenceScalars =
//     vtkDoubleArray::SafeDownCast(CTPersistenceDiagram_->
//       GetCellData()->GetArray("Persistence"));

//   vtkDoubleArray* birthScalars =
//     vtkDoubleArray::SafeDownCast(CTPersistenceDiagram_->
//       GetPointData()->GetArray("Birth"));

//   vtkDoubleArray* deathScalars =
//     vtkDoubleArray::SafeDownCast(CTPersistenceDiagram_->
//       GetPointData()->GetArray("Death"));

//   vtkPoints* points = (CTPersistenceDiagram_->GetPoints());
//   auto pairingsSize = (int) pairIdentifierScalars->GetNumberOfTuples();

//   // FIX : no more missed pairs
//   for(int pair_index = 0; pair_index<pairingsSize; pair_index++){
//     const float index_of_pair = pair_index;
//     if(*pairIdentifierScalars->GetTuple(pair_index)!=-1)
//         pairIdentifierScalars->SetTuple(pair_index, &index_of_pair);
//   }

//   auto s = (float) 0.0;

//   if (!deathScalars != !birthScalars) return -2;
//   bool Is3D = !(!deathScalars && !birthScalars);
//   if (!Is3D && diagramNumber == 1) s = (float) spacing;

//   if (pairingsSize < 1 || !vertexIdentifierScalars
//       || !pairIdentifierScalars || !nodeTypeScalars
//       || !persistenceScalars || !extremumIndexScalars || !points)
//     return -2;

//   diagram->resize(pairingsSize);
//   int nbNonCompact = 0;

//   for (int i = 0; i < pairingsSize; ++i) {

//     int vertexId1 = vertexIdentifierScalars->GetValue(2*i);
//     int vertexId2 = vertexIdentifierScalars->GetValue(2*i+1);
//     int nodeType1 = nodeTypeScalars->GetValue(2*i);
//     int nodeType2 = nodeTypeScalars->GetValue(2*i+1);

//     int pairIdentifier = pairIdentifierScalars->GetValue(i);
//     int pairType = extremumIndexScalars->GetValue(i);
//     double persistence = persistenceScalars->GetValue(i);

//     int index1 = 2*i;
//     double* coords1 = points->GetPoint(index1);
//     auto x1 = (float) coords1[0];
//     auto y1 = (float) coords1[1];
//     auto z1 = (float) coords1[2];

//     int index2 = index1 + 1;
//     double* coords2 = points->GetPoint(index2);
//     auto x2 = (float) coords2[0];
//     auto y2 = (float) coords2[1];
//     auto z2 = (float) coords2[2];

//     dataType value1 = (!birthScalars) ? (dataType) x1 :
//                       (dataType) birthScalars->GetValue(2*i);
//     dataType value2 = (!deathScalars) ? (dataType) y2 :
//                       (dataType) deathScalars->GetValue(2*i+1);

//     if (pairIdentifier != -1 && pairIdentifier < pairingsSize)
//       diagram->at(pairIdentifier) = std::make_tuple(
//         vertexId1, (BNodeType) nodeType1,
//         vertexId2, (BNodeType) nodeType2,
//         (dataType) persistence,
//         pairType,
//         value1, x1, y1, z1 + s,
//         value2, x2, y2, z2 + s
//       );

//     if (pairIdentifier >= pairingsSize) {
//       nbNonCompact++;
//       if (nbNonCompact == 0) {
//         std::stringstream msg;
//         msg << "[TTKBottleneckDistance] Diagram pair identifiers "
//             << "must be compact (not exceed the diagram size). "
//             << std::endl;
//         dMsg(std::cout, msg.str(), timeMsg);
//       }
//     }
//   }

//   if (nbNonCompact > 0) {
//     {
//       std::stringstream msg;
//       msg << "[TTKBottleneckDistance] Missed "
//           << nbNonCompact << " pairs due to non-compactness."
//           << std::endl;
//       dMsg(std::cout, msg.str(), timeMsg);
//     }
//   }


//   return 0;
// }

template <typename dataType>
vtkSmartPointer<vtkUnstructuredGrid> ttkPersistenceDiagramsClustering::createOutputCentroids(
         std::vector<std::vector<diagramTuple>>* final_centroids,
        std::vector<int> inv_clustering,
         double max_dimension)
{
    if(debugLevel_>3){
        std::cout<<"Creating vtk diagrams"<<std::endl;
    }
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

	vtkSmartPointer<vtkUnstructuredGrid> persistenceDiagram =
		vtkSmartPointer<vtkUnstructuredGrid>::New();

	vtkSmartPointer<vtkIntArray> nodeType =
		vtkSmartPointer<vtkIntArray>::New();
	nodeType->SetName("CriticalType");

	vtkSmartPointer<vtkDoubleArray> persistenceScalars =
		vtkSmartPointer<vtkDoubleArray>::New();
	persistenceScalars->SetName("Persistence");

	vtkSmartPointer<vtkIntArray> idOfPair =
		vtkSmartPointer<vtkIntArray>::New();
	idOfPair->SetName("ID of Pair");

	vtkSmartPointer<vtkDoubleArray> persistenceScalarsPoint =
		vtkSmartPointer<vtkDoubleArray>::New();
	persistenceScalarsPoint->SetName("Persistence");

	vtkSmartPointer<vtkIntArray> idOfDiagramPoint =
		vtkSmartPointer<vtkIntArray>::New();
	idOfDiagramPoint->SetName("ID of Cluster");

	vtkSmartPointer<vtkIntArray> pairType =
		vtkSmartPointer<vtkIntArray>::New();
	pairType->SetName("PairType");

	vtkSmartPointer<vtkFloatArray> coordsScalars=
		vtkSmartPointer<vtkFloatArray>::New();
	coordsScalars->SetNumberOfComponents(3);
	coordsScalars->SetName("Coordinates");


	int count = 0;
	for(unsigned int j = 0; j < final_centroids->size(); ++j){
		std::vector<diagramTuple> *diagram = &((*final_centroids)[j]);

		// First, add diagram points to the global input diagram
		for (unsigned int i = 0; i < diagram->size(); ++i) {
			vtkIdType ids[2];
			diagramTuple t = diagram->at(i);
			double x1 = std::get<6>(t); 
			double y1 = x1;
			x1 += 1.1*max_dimension*j; 
			double z1 = 1;  // Change 1 to j if you want to isolate the diagrams

			float coords1[3];
			coords1[0] = std::get<7>(t);
			coords1[1] = std::get<8>(t);
			coords1[2] = std::get<9>(t);

			double x2 = std::get<6>(t) ;
			double y2 = std::get<10>(t);
			double z2 = 1;  // Change 1 to j if you want to isolate the

			float coords2[3];
			coords2[0] = std::get<11>(t);
			coords2[1] = std::get<12>(t);
			coords2[2] = std::get<13>(t);

			idOfPair->InsertTuple1(count, i);

			points->InsertNextPoint(x1, y1, z1);
			coordsScalars->InsertTuple3(2*count, coords1[0], coords1[1], coords1[2]);
			idOfDiagramPoint->InsertTuple1(2*count, j);
			const ttk::CriticalType n1Type = std::get<1>(t);
			switch (n1Type) {
				case BLocalMin:
					nodeType->InsertTuple1(2*count, 0);
					break;

				case BSaddle1:
					nodeType->InsertTuple1(2*count, 1);
					break;

				case BSaddle2:
					nodeType->InsertTuple1(2*count, 2);
					break;

				case BLocalMax:
					nodeType->InsertTuple1(2*count, 3);
					break;
				default:
					nodeType->InsertTuple1(2*count, 0);
			}

			points->InsertNextPoint(x2+ 1.1*max_dimension*j, y2, z2);
			coordsScalars->InsertTuple3(2*count+1, coords2[0], coords2[1], coords2[2]);
			idOfDiagramPoint->InsertTuple1(2*count+1, j);
			const ttk::CriticalType n2Type = std::get<3>(t);
			switch (n2Type) {
				case BLocalMin:
					nodeType->InsertTuple1(2*count+1, 0);
					break;

				case BSaddle1:
					nodeType->InsertTuple1(2*count+1, 1);
					break;

				case BSaddle2:
					nodeType->InsertTuple1(2*count+1, 2);
					break;

				case BLocalMax:
					nodeType->InsertTuple1(2*count+1, 3);
					break;
				default:
					nodeType->InsertTuple1(2*count+1, 0);
			}

			ids[0] = 2*count;
			ids[1] = 2*count+1;

			persistenceDiagram->InsertNextCell(VTK_LINE, 2, ids);
			persistenceScalars->InsertTuple1(count, y2-x2);
			persistenceScalarsPoint->InsertTuple1(2*count, y2-x2);
			persistenceScalarsPoint->InsertTuple1(2*count+1, y2-x2);
			const ttk::SimplexId type = std::get<5>(t);
			switch (type) {
				case 0:
					pairType->InsertTuple1(count, 0);
					break;

				case 1:
					pairType->InsertTuple1(count, 1);
					break;

				case 2:
					pairType->InsertTuple1(count, 2);
					break;
				default:
					pairType->InsertTuple1(count, 0);
			}
			count++;
		}
	}

	persistenceDiagram->SetPoints(points);
	persistenceDiagram->GetCellData()->AddArray(persistenceScalars);
	persistenceDiagram->GetCellData()->AddArray(pairType);
	persistenceDiagram->GetCellData()->AddArray(idOfPair);
	persistenceDiagram->GetPointData()->AddArray(nodeType);
	persistenceDiagram->GetPointData()->AddArray(coordsScalars);
	persistenceDiagram->GetPointData()->AddArray(idOfDiagramPoint);
	persistenceDiagram->GetPointData()->AddArray(persistenceScalarsPoint);


  return persistenceDiagram;
    

}


template <typename dataType>
vtkSmartPointer<vtkUnstructuredGrid>
  ttkPersistenceDiagramsClustering::createOutputClusteredDiagrams(
    std::vector<std::vector<diagramTuple> > &all_CTDiagrams,
    std::vector<int> inv_clustering,
    double max_dimension)
{
    if(debugLevel_>3){
    std::cout<<"Creating vtk Outputs"<<std::endl;
    }
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

	vtkSmartPointer<vtkUnstructuredGrid> persistenceDiagram =
		vtkSmartPointer<vtkUnstructuredGrid>::New();

	vtkSmartPointer<vtkIntArray> nodeType =
		vtkSmartPointer<vtkIntArray>::New();
	nodeType->SetName("CriticalType");

	vtkSmartPointer<vtkDoubleArray> persistenceScalars =
		vtkSmartPointer<vtkDoubleArray>::New();
	persistenceScalars->SetName("Persistence");

	vtkSmartPointer<vtkIntArray> idOfPair =
		vtkSmartPointer<vtkIntArray>::New();
	idOfPair->SetName("ID of Pair");

	vtkSmartPointer<vtkDoubleArray> persistenceScalarsPoint =
		vtkSmartPointer<vtkDoubleArray>::New();
	persistenceScalarsPoint->SetName("Persistence");

	vtkSmartPointer<vtkIntArray> idOfDiagramPoint =
		vtkSmartPointer<vtkIntArray>::New();
	idOfDiagramPoint->SetName("ID of Diagram");

	vtkSmartPointer<vtkIntArray> idOfCluster =
		vtkSmartPointer<vtkIntArray>::New();
	idOfCluster->SetName("ID of Cluster");

	vtkSmartPointer<vtkIntArray> pairType =
		vtkSmartPointer<vtkIntArray>::New();
	pairType->SetName("PairType");

	vtkSmartPointer<vtkFloatArray> coordsScalars=
		vtkSmartPointer<vtkFloatArray>::New();
	coordsScalars->SetNumberOfComponents(3);
	coordsScalars->SetName("Coordinates");


    std::vector<int> cluster_size;

	int count = 0;
	for(unsigned int j = 0; j < all_CTDiagrams.size(); ++j){
		std::vector<diagramTuple> *diagram = &(all_CTDiagrams[j]);
        
        unsigned int c = inv_clustering[j];
        if(c+1>cluster_size.size()){
            cluster_size.resize(c+1);
            cluster_size[c]=1;
        }
        else{
            cluster_size[c]++;
        }
		// First, add diagram points to the global input diagram
		for (unsigned int i = 0; i < diagram->size(); ++i) {
			vtkIdType ids[2];
			diagramTuple t = diagram->at(i);
			double x1 = std::get<6>(t);
			double y1 = x1;
			x1 +=  1.1*max_dimension*c;
			double z1 = 1-cluster_size[c];  // Change 1 to j if you want to isolate the diagrams

			float coords1[3];
			coords1[0] = std::get<7>(t);
			coords1[1] = std::get<8>(t);
			coords1[2] = std::get<9>(t);

			double x2 = std::get<6>(t) + (1.1*max_dimension*c);
			double y2 = std::get<10>(t);
			double z2 = 1-cluster_size[c];  // Change 1 to j if you want to isolate the

			float coords2[3];
			coords2[0] = std::get<11>(t);
			coords2[1] = std::get<12>(t);
			coords2[2] = std::get<13>(t);

			idOfPair->InsertTuple1(count, i);

			points->InsertNextPoint(x1, y1, z1);
			coordsScalars->InsertTuple3(2*count, coords1[0], coords1[1], coords1[2]);
			idOfDiagramPoint->InsertTuple1(2*count, j);
            // std::cout<<"\nMAX DIM \n"<<max_dimension<<std::endl;
			idOfCluster->InsertTuple1(2*count, c);
			const ttk::CriticalType n1Type = std::get<1>(t);
			switch (n1Type) {
				case BLocalMin:
					nodeType->InsertTuple1(2*count, 0);
					break;

				case BSaddle1:
					nodeType->InsertTuple1(2*count, 1);
					break;

				case BSaddle2:
					nodeType->InsertTuple1(2*count, 2);
					break;

				case BLocalMax:
					nodeType->InsertTuple1(2*count, 3);
					break;
				default:
					nodeType->InsertTuple1(2*count, 0);
			}

			points->InsertNextPoint(x2, y2, z2);
			coordsScalars->InsertTuple3(2*count+1, coords2[0], coords2[1], coords2[2]);
			idOfDiagramPoint->InsertTuple1(2*count+1, j);
			idOfCluster->InsertTuple1(2*count+1, c);
			const ttk::CriticalType n2Type = std::get<3>(t);
			switch (n2Type) {
				case BLocalMin:
					nodeType->InsertTuple1(2*count+1, 0);
					break;

				case BSaddle1:
					nodeType->InsertTuple1(2*count+1, 1);
					break;

				case BSaddle2:
					nodeType->InsertTuple1(2*count+1, 2);
					break;

				case BLocalMax:
					nodeType->InsertTuple1(2*count+1, 3);
					break;
				default:
					nodeType->InsertTuple1(2*count+1, 0);
			}

			ids[0] = 2*count;
			ids[1] = 2*count+1;

			persistenceDiagram->InsertNextCell(VTK_LINE, 2, ids);
			persistenceScalars->InsertTuple1(count, y2-x2);
			persistenceScalarsPoint->InsertTuple1(2*count, y2-x2);
			persistenceScalarsPoint->InsertTuple1(2*count+1, y2-x2);
			const ttk::SimplexId type = std::get<5>(t);
			switch (type) {
				case 0:
					pairType->InsertTuple1(count, 0);
					break;

				case 1:
					pairType->InsertTuple1(count, 1);
					break;

				case 2:
					pairType->InsertTuple1(count, 2);
					break;
				default:
					pairType->InsertTuple1(count, 0);
			}
			count++;
		}
	}

	persistenceDiagram->SetPoints(points);
	persistenceDiagram->GetCellData()->AddArray(persistenceScalars);
	persistenceDiagram->GetCellData()->AddArray(pairType);
	persistenceDiagram->GetCellData()->AddArray(idOfPair);
	persistenceDiagram->GetPointData()->AddArray(nodeType);
	persistenceDiagram->GetPointData()->AddArray(coordsScalars);
	persistenceDiagram->GetPointData()->AddArray(idOfDiagramPoint);
	persistenceDiagram->GetPointData()->AddArray(idOfCluster);
	persistenceDiagram->GetPointData()->AddArray(persistenceScalarsPoint);


  return persistenceDiagram;
}

#endif // _TTK_PERSISTENCEDIAGRAMSCLUSTERING_H
