/// \ingroup base
/// \class ttk::BettiNumbers
/// \author Yoann Coudert--Osmont <ycoudert@ens-paris-saclay.fr>
/// \date October 2019
///
/// \brief TTK bettiNumbers processing package.

#pragma once

// base code includes
#include <Triangulation.h>
#include <Wrapper.h>
#include <UnionFind.h>

namespace ttk {

	namespace bettiNumbers {

		class BettiNumbers : public Debug {

		public:
			BettiNumbers();

			~BettiNumbers();

			/// Execute the package.
			/// \pre If this TTK package uses ttk::Triangulation for fast mesh
			/// traversals, the function setupTriangulation() must be called on this
			/// object prior to this function, in a clearly distinct pre-processing
			/// steps. An error will be returned otherwise.
			/// \note In such a case, it is recommended to exclude
			/// setupTriangulation() from any time performance measurement.
			/// \param argment Dummy integer argument.
			/// \return Returns 0 upon success, negative values otherwise.
			template <class dataType>
			int execute(const int &argument) const;

			/// Pass a pointer to an input array representing a scalarfield.
			/// The expected format for the array is the following:
			/// <vertex0-component0> <vertex0-component1> ... <vertex0-componentN>
			/// <vertex1-component0> <vertex1-component1> ... <vertex1-componentN>
			/// <vertexM-component0> <vertexM-component1> ... <vertexM-componentN>.
			/// The array is expected to be correctly allocated.
			/// \param data Pointer to the data array.
			/// \return Returns 0 upon success, negative values otherwise.
			/// \sa setVertexNumber() and setDimensionNumber().
			inline int setInputDataPointer(void *data) {
				inputData_ = data;
				return 0;
			}

			/// Pass a pointer to an output array representing a scalar field.
			/// The expected format for the array is the following:
			/// <vertex0-component0> <vertex0-component1> ... <vertex0-componentN>
			/// <vertex1-component0> <vertex1-component1> ... <vertex1-componentN>
			/// <vertexM-component0> <vertexM-component1> ... <vertexM-componentN>.
			/// The array is expected to be correctly allocated.
			/// \param data Pointer to the data array.
			/// \return Returns 0 upon success, negative values otherwise.
			/// \sa setVertexNumber() and setDimensionNumber().
			inline int setOutputDataPointer(void *data) {
				outputData_ = data;
				return 0;
			}

			// General documentation info:
			//
			/// Setup a (valid) triangulation object for this TTK base object.
			///
			/// \pre This function should be called prior to any usage of this TTK
			/// object, in a clearly distinct pre-processing step that involves no
			/// traversal or computation at all. An error will be returned otherwise.
			///
			/// \note It is recommended to exclude this pre-processing function from
			/// any time performance measurement. Therefore, it is recommended to
			/// call this function ONLY in the pre-processing steps of your program.
			/// Note however, that your triangulation object must be valid when
			/// calling this function (i.e. you should have filled it at this point,
			/// see the setInput*() functions of ttk::Triangulation). See ttkBettiNumbers
			/// for further examples.
			///
			/// \param triangulation Pointer to a valid triangulation.
			/// \return Returns 0 upon success, negative values otherwise.
			/// \sa ttk::Triangulation
			inline int setupTriangulation(Triangulation *triangulation) {
				triangulation_ = triangulation;
				if(triangulation_) {
					triangulation_->preprocessEdges();
					triangulation_->preprocessTriangles();
					triangulation_->preprocessCellTriangles();
				}
				return 0;
			}

		protected:
			void *inputData_, *outputData_;
			Triangulation *triangulation_;
		};
	} // namespace bettiNumbers
} // namespace ttk

// template functions
template <class dataType>
int ttk::bettiNumbers::BettiNumbers::execute(const int &argument) const {

	Timer ti;
	#ifndef TTK_ENABLE_KAMIKAZE
	if(!triangulation_) return -1;
	if(!inputData_) return -2;
	if(!outputData_) return -3;
	#endif

	SimplexId vertexNumber = triangulation_->getNumberOfVertices();
	dataType *outputData = (dataType *)outputData_;
	dataType *inputData = (dataType *)inputData_;
	for(SimplexId i = 0; i < vertexNumber; i++)
		outputData[i] = inputData[i];

	/**** B0 computation ****/
	std::vector<UnionFind> uV(vertexNumber);
	SimplexId edgeNumber = triangulation_->getNumberOfEdges();
	for(SimplexId e = 0; e < edgeNumber; e++) {
		SimplexId a, b;
		triangulation_->getEdgeVertex(e, 0, a);
		triangulation_->getEdgeVertex(e, 1, b);
		makeUnion(&uV[a], &uV[b]);
	}
	int b0 = 0;
	for(UnionFind &u : uV)
		if(&u == u.find()) b0 ++;
	/************************/

	/**** B2 computation ****/
	SimplexId triangleNumber = triangulation_->getNumberOfTriangles();
	SimplexId cellNumber = triangulation_->getNumberOfCells();
	std::vector<int> ncells(triangleNumber, 0);
	for(SimplexId c = 0; c < cellNumber; c++) {
		SimplexId nfaces = triangulation_->getCellTriangleNumber(c);
		for(SimplexId i = 0; i < nfaces; i++) {
			SimplexId t;
			triangulation_->getCellTriangle(c, i, t);
			ncells[t] ++;
		}
	}
	std::vector<bool> surfaceVertex(vertexNumber, false);
	for(UnionFind &u : uV) u.setParent(&u);
	for(SimplexId t = 0; t < triangleNumber; t++) {
		if(ncells[t] > 1) continue;
		SimplexId v, v2=0;
		for(SimplexId i = 0; i < 3; i++) {
			triangulation_->getTriangleVertex(t, i, v);
			surfaceVertex[v] = true;
			if(i > 0) makeUnion(&uV[v], &uV[v2]);
			v2 = v;
		}
	}
	int b2 = -b0;
	for(SimplexId v = 0; v < vertexNumber; v++)
		if(surfaceVertex[v] && &uV[v] == uV[v].find()) b2 ++;
	/************************/

	/**** X computation ****/
	int X = vertexNumber - edgeNumber + triangleNumber - cellNumber;

	/**** B1 computation ****/
	int b1 = b0 + b2 - X;

	{
		std::stringstream msg;
		msg << "[BettiNumbers] Data-set (" << vertexNumber << " points, " << edgeNumber << " edges, "
				<< triangleNumber << " triangles and " << cellNumber << " cells) processed in "
				<< ti.getElapsedTime() << " s. (" << threadNumber_ << " thread(s))."
				<< std::endl;
		dMsg(std::cout, msg.str(), timeMsg);
	}
	std::cout << std::endl << "Hello World!" << std::endl;
	std::cout << "B0 = " << b0 << std::endl;
	std::cout << "B1 = " << b1 << std::endl;
	std::cout << "B2 = " << b2 << std::endl;
	std::cout << "X = " << X << std::endl;
	std::cout << std::endl;

	return 0;
}
