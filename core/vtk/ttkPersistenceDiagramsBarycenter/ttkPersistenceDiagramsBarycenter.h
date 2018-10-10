/// \ingroup vtk
/// \class ttkPersistenceDiagramsBarycenter
/// \author Michael Michaux <michauxmichael89@gmail.com>
/// \date August 2016.
///
/// \brief TTK VTK-filter that takes an input ensemble data set
/// (represented by a list of scalar fields) and which computes various
/// vertexwise statistics (PDF estimation, bounds, moments, etc.)
///
/// \param Input0 Input ensemble scalar field #0 (vtkDataSet)
/// \param Input1 Input ensemble scalar field #1 (vtkDataSet)\n
/// ...\n
/// \param InputN Input ensemble scalar field #N (vtkDataSet)
/// \param Output0 Lower and upper bound fields (vtkDataSet)
/// \param Output1 Histogram estimations of the vertex probability density
/// functions (vtkDataSet)
/// \param Output2 Mean field (vtkDataSet)
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the corresponding ParaView state file example for a usage example
/// within a VTK pipeline.
///
/// \sa vtkMandatoryCriticalPoints
/// \sa ttk::PersistenceDiagramsBarycenter
#ifndef _TTK_PERSISTENCEDIAGRAMSBARYCENTER_H
#define _TTK_PERSISTENCEDIAGRAMSBARYCENTER_H


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
#include                  <PersistenceDiagramsBarycenter.h>
#include                  <Wrapper.h>
#include                  <ttkWrapper.h>

// in this example, this wrapper takes a data-set on the input and produces a
// data-set on the output - to adapt.
// see the documentation of the vtkAlgorithm class to decide from which VTK
// class your wrapper should inherit.
#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkPersistenceDiagramsBarycenter
#else
class ttkPersistenceDiagramsBarycenter
#endif
  : public vtkDataSetAlgorithm, public ttk::Wrapper{

  public:

    static ttkPersistenceDiagramsBarycenter* New();

    vtkTypeMacro(ttkPersistenceDiagramsBarycenter, vtkDataSetAlgorithm);

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

	vtkSetMacro(ReinitPrices, int);
    vtkGetMacro(ReinitPrices, int);

	vtkSetMacro(EpsilonDecreases, int);
    vtkGetMacro(EpsilonDecreases, int);

	vtkSetMacro(EarlyStoppage, int);
    vtkGetMacro(EarlyStoppage, int);

	vtkSetMacro(ThreadNumber, int);
    vtkGetMacro(ThreadNumber, int);

	vtkSetMacro(Alpha, double);
    vtkGetMacro(Alpha, double);

    vtkSetMacro(Method, int);
    vtkGetMacro(Method, int);
    vtkSetMacro(Deterministic, bool);
    vtkGetMacro(Deterministic, bool);



  protected:

    ttkPersistenceDiagramsBarycenter();

    ~ttkPersistenceDiagramsBarycenter();

	template <typename dataType>
    int getPersistenceDiagram(
      std::vector<diagramTuple>* diagram,
      vtkUnstructuredGrid *CTPersistenceDiagram_,
      const double spacing,
      const int diagramNumber);

	template<typename dataType>
	vtkSmartPointer<vtkUnstructuredGrid> createPersistenceDiagram(
    const std::vector<diagramTuple>* diagram);

	template<typename dataType>
	vtkSmartPointer<vtkUnstructuredGrid> createOutputDiagrams(
    std::vector<std::vector<diagramTuple>> &all_CTDiagrams);

	template<typename dataType>
	vtkSmartPointer<vtkUnstructuredGrid> createMatchings(
		const std::vector<std::vector<matchingTuple>> *matchings,
    const std::vector<diagramTuple>* barycenter,
    std::vector<std::vector<diagramTuple>> &all_CTDiagrams);

	int FillInputPortInformation(int port, vtkInformation *info);
    int FillOutputPortInformation(int port, vtkInformation *info);

    int RequestData(vtkInformation *request,
      vtkInformationVector **inputVector, vtkInformationVector *outputVector);


  private:
    bool                  Deterministic;
    bool                  UseAllCores;
    int                   ThreadNumber;
    int                   Method;
	bool                  UseOutputMatching;
	double                Alpha;
    std::string                ScalarField;
	std::string                WassersteinMetric;

	bool	              UseProgressive;
	double	              TimeLimit;

	bool ReinitPrices;
	bool EarlyStoppage;
	bool EpsilonDecreases;

    // base code features
    int doIt(vtkDataSet **input,
			 vtkUnstructuredGrid *outputBarycenter,
			 vtkUnstructuredGrid *matchings,
			 vtkUnstructuredGrid *outputDiagrams,
             int numInputs);

    bool needsToAbort();

    int updateProgress(const float &progress);

};

template <typename dataType>
int ttkPersistenceDiagramsBarycenter::getPersistenceDiagram(
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

  vtkDoubleArray* deathScalars =
    vtkDoubleArray::SafeDownCast(CTPersistenceDiagram_->
      GetPointData()->GetArray("Death"));

  vtkPoints* points = (CTPersistenceDiagram_->GetPoints());
  int pairingsSize = (int) pairIdentifierScalars->GetNumberOfTuples();
  auto s = (float) 0.0;

  if (!deathScalars != !birthScalars) return -2;
  bool Is3D = !(!deathScalars && !birthScalars);
  if (!Is3D && diagramNumber == 1) s = (float) spacing;

  if (pairingsSize < 1 || !vertexIdentifierScalars
      || !pairIdentifierScalars || !nodeTypeScalars
      || !persistenceScalars || !extremumIndexScalars || !points)
    return -2;

  diagram->resize(pairingsSize);
  int nbNonCompact = 0;

  for (int i = 0; i < pairingsSize; ++i) {

    int vertexId1 = vertexIdentifierScalars->GetValue(2*i);
    int vertexId2 = vertexIdentifierScalars->GetValue(2*i+1);
    int nodeType1 = nodeTypeScalars->GetValue(2*i);
    int nodeType2 = nodeTypeScalars->GetValue(2*i+1);

    int pairIdentifier = pairIdentifierScalars->GetValue(i);
    int pairType = extremumIndexScalars->GetValue(i);
    double persistence = persistenceScalars->GetValue(i);

    int index1 = 2*i;
    double* coords1 = points->GetPoint(index1);
    auto x1 = (float) coords1[0];
    auto y1 = (float) coords1[1];
    auto z1 = (float) coords1[2];

    int index2 = index1 + 1;
    double* coords2 = points->GetPoint(index2);
    auto x2 = (float) coords2[0];
    auto y2 = (float) coords2[1];
    auto z2 = (float) coords2[2];

    dataType value1 = (!birthScalars) ? (dataType) x1 :
                      (dataType) birthScalars->GetValue(2*i);
    dataType value2 = (!deathScalars) ? (dataType) y2 :
                      (dataType) deathScalars->GetValue(2*i+1);

    if (pairIdentifier != -1 && pairIdentifier < pairingsSize)
      diagram->at(pairIdentifier) = std::make_tuple(
        vertexId1, (BNodeType) nodeType1,
        vertexId2, (BNodeType) nodeType2,
        (dataType) persistence,
        pairType,
        value1, x1, y1, z1 + s,
        value2, x2, y2, z2 + s
      );

    if (pairIdentifier >= pairingsSize) {
      nbNonCompact++;
      if (nbNonCompact == 0) {
        std::stringstream msg;
        msg << "[TTKPersistenceDiagramBarycenter] Diagram pair identifiers "
            << "must be compact (not exceed the diagram size). "
            << std::endl;
        dMsg(std::cout, msg.str(), timeMsg);
      }
    }
  }

  if (nbNonCompact > 0) {
    {
      std::stringstream msg;
      msg << "[TTKPersistenceDiagramBarycenter] Missed "
          << nbNonCompact << " pairs due to non-compactness."
          << std::endl;
      dMsg(std::cout, msg.str(), timeMsg);
    }
  }


  return 0;
}


template <typename dataType>
vtkSmartPointer<vtkUnstructuredGrid> ttkPersistenceDiagramsBarycenter::createPersistenceDiagram(
    const std::vector<diagramTuple>* diagram)
{
	printf("Creating vtk Diagram");
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

	vtkSmartPointer<vtkUnstructuredGrid> persistenceDiagram =
		vtkSmartPointer<vtkUnstructuredGrid>::New();

	vtkSmartPointer<vtkIntArray> nodeType =
		vtkSmartPointer<vtkIntArray>::New();
	nodeType->SetName("nodeType");

	vtkSmartPointer<vtkDoubleArray> persistenceScalars =
		vtkSmartPointer<vtkDoubleArray>::New();
	persistenceScalars->SetName("Persistence");

	vtkSmartPointer<vtkIntArray> pairType =
		vtkSmartPointer<vtkIntArray>::New();
	pairType->SetName("pairType");

	vtkSmartPointer<vtkFloatArray> coordsScalars=
		vtkSmartPointer<vtkFloatArray>::New();
	coordsScalars->SetNumberOfComponents(3);
	coordsScalars->SetName("Coordinates");

	for (unsigned int i = 0; i < diagram->size(); ++i) {
		vtkIdType ids[2];
		diagramTuple t = diagram->at(i);
		double x1 = std::get<6>(t);
		double y1 = x1;
		double z1 = 0;

		float coords1[3];
		coords1[0] = std::get<7>(t);
		coords1[1] = std::get<8>(t);
		coords1[2] = std::get<9>(t);

		double x2 = std::get<6>(t);
		double y2 = std::get<10>(t);
		double z2 = 0;

		float coords2[3];
		coords2[0] = std::get<11>(t);
		coords2[1] = std::get<12>(t);
		coords2[2] = std::get<13>(t);


		points->InsertNextPoint(x1, y1, z1);
		coordsScalars->InsertTuple3(2*i, coords1[0], coords1[1], coords1[2]);
		const ttk::CriticalType n1Type = std::get<1>(t);
		switch (n1Type) {
			case BLocalMin:
				nodeType->InsertTuple1(2*i, 0);
				break;

			case BSaddle1:
				nodeType->InsertTuple1(2*i, 1);
				break;

			case BSaddle2:
				nodeType->InsertTuple1(2*i, 2);
				break;

			case BLocalMax:
				nodeType->InsertTuple1(2*i, 3);
				break;
			default:
				nodeType->InsertTuple1(2*i, 0);
		}
		points->InsertNextPoint(x2, y2, z2);
		coordsScalars->InsertTuple3(2*i+1, coords2[0], coords2[1], coords2[2]);

		const ttk::CriticalType n2Type = std::get<3>(t);
		switch (n2Type) {
			case BLocalMin:
				nodeType->InsertTuple1(2*i+1, 0);
				break;

			case BSaddle1:
				nodeType->InsertTuple1(2*i+1, 1);
				break;

			case BSaddle2:
				nodeType->InsertTuple1(2*i+1, 2);
				break;

			case BLocalMax:
				nodeType->InsertTuple1(2*i+1, 3);
				break;
		default :
			nodeType->InsertTuple1(2*i+1, 0);
		}

		ids[0] = 2*i;
		ids[1] = 2*i+1;

		persistenceDiagram->InsertNextCell(VTK_LINE, 2, ids);
		persistenceScalars->InsertTuple1(i, y2-x2);
		const ttk::SimplexId type = std::get<5>(t);
		switch (type) {
			case 0:
				pairType->InsertTuple1(i, 0);
				break;

			case 1:
				pairType->InsertTuple1(i, 1);
				break;

			case 2:
				pairType->InsertTuple1(i, 2);
				break;
			default:
				pairType->InsertTuple1(i, 0);
		}
	}

	persistenceDiagram->SetPoints(points);
	persistenceDiagram->GetCellData()->AddArray(persistenceScalars);
	persistenceDiagram->GetCellData()->AddArray(pairType);
	persistenceDiagram->GetPointData()->AddArray(nodeType);
	persistenceDiagram->GetPointData()->AddArray(coordsScalars);

	return persistenceDiagram;
}


template <typename dataType>
vtkSmartPointer<vtkUnstructuredGrid>
  ttkPersistenceDiagramsBarycenter::createOutputDiagrams(
    std::vector<std::vector<diagramTuple> > &all_CTDiagrams)
{
	printf("Creating vtk Outputs");
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

	vtkSmartPointer<vtkUnstructuredGrid> persistenceDiagram =
		vtkSmartPointer<vtkUnstructuredGrid>::New();

	vtkSmartPointer<vtkIntArray> nodeType =
		vtkSmartPointer<vtkIntArray>::New();
	nodeType->SetName("nodeType");

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

	vtkSmartPointer<vtkIntArray> pairType =
		vtkSmartPointer<vtkIntArray>::New();
	pairType->SetName("pairType");

	vtkSmartPointer<vtkFloatArray> coordsScalars=
		vtkSmartPointer<vtkFloatArray>::New();
	coordsScalars->SetNumberOfComponents(3);
	coordsScalars->SetName("Coordinates");


	int count = 0;
	for(unsigned int j = 0; j < all_CTDiagrams.size(); ++j){
		std::vector<diagramTuple> *diagram = &(all_CTDiagrams[j]);

		// First, add diagram points to the global input diagram
		for (unsigned int i = 0; i < diagram->size(); ++i) {
			vtkIdType ids[2];
			diagramTuple t = diagram->at(i);
			double x1 = std::get<6>(t);
			double y1 = x1;
			double z1 = j;  // Change 1 to j if you want to isolate the diagrams

			float coords1[3];
			coords1[0] = std::get<7>(t);
			coords1[1] = std::get<8>(t);
			coords1[2] = std::get<9>(t);

			double x2 = std::get<6>(t);
			double y2 = std::get<10>(t);
			double z2 = j;  // Change 1 to j if you want to isolate the

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

			points->InsertNextPoint(x2, y2, z2);
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
  ttkPersistenceDiagramsBarycenter::createMatchings(
    const std::vector<std::vector<matchingTuple>> *matchings,
    const std::vector<diagramTuple>* barycenter,
    std::vector<std::vector<diagramTuple>> &all_CTDiagrams)
{
	printf("Creating vtk Matching");
	vtkSmartPointer<vtkPoints> matchingPoints = vtkSmartPointer<vtkPoints>::New();

	vtkSmartPointer<vtkUnstructuredGrid> matchingMesh =
		vtkSmartPointer<vtkUnstructuredGrid>::New();

	vtkSmartPointer<vtkIntArray> idOfDiagramMatchingPoint =
		vtkSmartPointer<vtkIntArray>::New();
	idOfDiagramMatchingPoint->SetName("ID of Diagram");

	vtkSmartPointer<vtkIntArray> idOfPoint =
		vtkSmartPointer<vtkIntArray>::New();
	idOfPoint->SetName("ID of Point");

	vtkSmartPointer<vtkIntArray> idOfDiagramMatching =
		vtkSmartPointer<vtkIntArray>::New();
	idOfDiagramMatching->SetName("ID of Diagram");

	vtkSmartPointer<vtkDoubleArray> cost =
		vtkSmartPointer<vtkDoubleArray>::New();
	cost->SetName("Cost");

	int count=0;
	for(unsigned int j = 0; j < all_CTDiagrams.size(); ++j){
		std::vector<diagramTuple> *diagram = &(all_CTDiagrams[j]);
		std::vector<matchingTuple> matchings_j = matchings->at(j);


		for(unsigned int i = 0; i < matchings_j.size(); ++i) {
			vtkIdType ids[2];
			ids[0] = 2*count;
			ids[1] = 2*count+1;

			matchingTuple m = matchings_j[i];
			int bidder_id = std::get<0>(m);
			int good_id = std::get<1>(m);

			diagramTuple t1 = barycenter->at(good_id);
			double x1 = std::get<6>(t1);
			double y1 = std::get<10>(t1);
			double z1 = 0;

			diagramTuple t2 = diagram->at(bidder_id);
			double x2 = std::get<6>(t2);
			double y2 = std::get<10>(t2);
			double z2 = 1;  // Change 1 to j if you want to isolate the diagrams

			if(y2>x2){
				matchingPoints->InsertNextPoint(x1, y1, z1);
				matchingPoints->InsertNextPoint(x2, y2, z2);
				matchingMesh->InsertNextCell(VTK_LINE, 2, ids);
				idOfDiagramMatching->InsertTuple1(count, j);
				cost->InsertTuple1(count, std::get<2>(m));
				idOfDiagramMatchingPoint->InsertTuple1(2*count, j);
				idOfDiagramMatchingPoint->InsertTuple1(2*count+1, j);
				idOfPoint->InsertTuple1(2*count+1, good_id);
				idOfPoint->InsertTuple1(2*count+1, bidder_id);

				count++;
			}
		}

	}

	matchingMesh->SetPoints(matchingPoints);
	matchingMesh->GetPointData()->AddArray(idOfDiagramMatchingPoint);
	matchingMesh->GetPointData()->AddArray(idOfPoint);
	matchingMesh->GetCellData()->AddArray(idOfDiagramMatching);
	matchingMesh->GetCellData()->AddArray(cost);

  return matchingMesh;
}


#endif // _TTK_PERSISTENCEDIAGRAMSBARYCENTER_H
