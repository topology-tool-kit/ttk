/// \defgroup examples examples
/// \author Yoann Coudert--Osmont <ycoudert@ens-paris-saclay.fr>
/// 		with inspiration from a piece of Adrien's code to compute B1
/// \date October 2019.
///
/// \brief A Program to compute the betti numbers and the euler caracteristic
/// of the VTU file given in argument

#include <CommandLineParser.h>
#include <UnionFind.h>

#include <ttkTriangulation.h>
// #include <ttkBettiNumbers.h>

#include <vtkXMLUnstructuredGridReader.h>
#include <vtkSmartPointer.h>

using namespace ttk;

int main(int argc, char **argv) {

	ttk::globalDebugLevel_ = 3;

	CommandLineParser parser;
	std::string inputFilePath;
	parser.setArgument("i", &inputFilePath, "Path to input VTU file");
	parser.parse(argc, argv);

	vtkSmartPointer<vtkXMLUnstructuredGridReader> reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
	reader->SetFileName(inputFilePath.data());
	reader->Update();

	// vtkSmartPointer<ttkBettiNumbers> betti = vtkSmartPointer<ttkBettiNumbers>::New();
	// betti->SetInputConnection(reader->GetOutputPort());
	// betti->Update();

	vtkUnstructuredGrid* ug = reader->GetOutput();
	Triangulation triangulation;
	vtkIdType np = ug->GetNumberOfPoints();
	std::vector<float> points;
	for(vtkIdType i = 0; i < np; i++) {
		double* p = ug->GetPoint(i);
		for(int j = 0; j < 3; j++) points.push_back(p[j]);
	}
	triangulation.setInputPoints(np, points.data());
	vtkIdType nc = ug->GetNumberOfCells();
	std::vector<LongSimplexId> cells;
	for(vtkIdType i = 0; i < nc; i++) {
		vtkCell* c = ug->GetCell(i);
		cells.push_back(4);
		for(int j = 0; j < 4; j++) cells.push_back(c->GetPointId(j));
	}
	triangulation.setInputPoints(np, points.data());
	triangulation.setInputCells(nc, cells.data());

	triangulation.preprocessCellTriangles();
	triangulation.preprocessTriangles();
	triangulation.preprocessTriangleEdges();
	triangulation.preprocessEdges();

	Timer ti;

	/**** B0 computation ****/
	SimplexId vertexNumber = triangulation.getNumberOfVertices();
	std::vector<UnionFind> uV(vertexNumber);
	SimplexId edgeNumber = triangulation.getNumberOfEdges();
	for(SimplexId e = 0; e < edgeNumber; e++) {
		SimplexId a, b;
		triangulation.getEdgeVertex(e, 0, a);
		triangulation.getEdgeVertex(e, 1, b);
		makeUnion(&uV[a], &uV[b]);
	}
	int b0 = 0;
	for(UnionFind &u : uV)
		if(&u == u.find()) b0 ++;
	/************************/

	/**** B2 and B1 computation ****/

	SimplexId triangleNumber = triangulation.getNumberOfTriangles();
	SimplexId cellNumber = triangulation.getNumberOfCells();

	// Computation of the number of cells to which the triangles belongs
	std::vector<int> ncells(triangleNumber, 0);
	for(SimplexId c = 0; c < cellNumber; c++) {
		SimplexId nfaces = triangulation.getCellTriangleNumber(c);
		for(SimplexId i = 0; i < nfaces; i++) {
			SimplexId t;
			triangulation.getCellTriangle(c, i, t);
			ncells[t] ++;
		}
	}

	// Computation of connected components on surfaces in the UnionFind u
	// And computations of numbers of vertices, edges and triangles on surfaces
	for(UnionFind &u : uV) u.setParent(&u);
	std::vector<bool> surfaceVertex(vertexNumber, false);
	std::vector<bool> surfaceEdge(edgeNumber, false);
	int surfaceVertexNumber = 0;
	int surfaceEdgeNumber = 0;
	int surfaceTriangleNumber = 0;
	for(SimplexId t = 0; t < triangleNumber; t++) {
		if(ncells[t] > 1) continue;
		// Then the triangle belongs to a surface
		surfaceTriangleNumber ++;
		// Treatment of vertices
		SimplexId v, v2=0;
		for(SimplexId i = 0; i < 3; i++) {
			triangulation.getTriangleVertex(t, i, v);
			if(!surfaceVertex[v]) {
				surfaceVertex[v] = true;
				surfaceVertexNumber ++;
			}
			if(i > 0) makeUnion(&uV[v], &uV[v2]);
			v2 = v;
		}
		// Treatment of edges
		SimplexId nedges = triangulation.getTriangleEdgeNumber(t);
		for(SimplexId i = 0; i < nedges; i++) {
			SimplexId e;
			triangulation.getTriangleEdge(t, i, e);
			if(!surfaceEdge[e]) {
				surfaceEdge[e] = true;
				surfaceEdgeNumber ++;
			}
		}
	}

	// B2 computation
	int b2 = -b0;
	for(SimplexId v = 0; v < vertexNumber; v++)
		if(surfaceVertex[v] && &uV[v] == uV[v].find()) b2 ++;
	
	// B1 computation
	int b1 = b0+b2 - (surfaceVertexNumber - surfaceEdgeNumber + surfaceTriangleNumber) / 2;
	/************************/

	/**** X computation ****/
	int X = vertexNumber - edgeNumber + triangleNumber - cellNumber;

	/**** B1 computation thanks to Euler ****/
	int b1Bis = b0 + b2 - X;

	cout << "[BettiNumbers] Data-set (" << vertexNumber << " points, " << edgeNumber << " edges, "
			<< triangleNumber << " triangles and " << cellNumber << " cells) processed in "
			<< ti.getElapsedTime() << " s. (" << 1 << " thread(s))."
			<< std::endl;
	std::cout << std::endl << "Hello World!" << std::endl;
	std::cout << "B0 = " << b0 << std::endl;
	std::cout << "B1 = " << b1 << "  (" << b1Bis << ")" << std::endl;
	std::cout << "B2 = " << b2 << std::endl;
	std::cout << "X = " << X << std::endl;
	std::cout << std::endl;

	return 0;
}

/// @}
