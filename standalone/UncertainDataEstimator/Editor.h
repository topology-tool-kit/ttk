/// \author Michael Michaux <michaelmichaux89@gmail.com>.
/// \date May 2016.
///
/// \brief Data-structures and processing.

#pragma once

#include <cmath>
#include <stdexcept>

// VTK includes
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
// Data types
#include <vtkDataSet.h>
#include <vtkImageData.h>
#include <vtkPolyData.h>
#include <vtkRectilinearGrid.h>
#include <vtkStructuredGrid.h>
#include <vtkUnstructuredGrid.h>
// Readers
#include <vtkDataSetReader.h>
#include <vtkImageReader2.h>
#include <vtkXMLImageDataReader.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLRectilinearGridReader.h>
#include <vtkXMLStructuredGridReader.h>
#include <vtkXMLUnstructuredGridReader.h>
// Writers
#include <vtkDataSetWriter.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLRectilinearGridWriter.h>
#include <vtkXMLStructuredGridWriter.h>
#include <vtkXMLUnstructuredGridWriter.h>
// Array types
#include <vtkCharArray.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkIntArray.h>
#include <vtkShortArray.h>
#include <vtkSignedCharArray.h>
#include <vtkUnsignedCharArray.h>
#include <vtkUnsignedIntArray.h>
#include <vtkUnsignedShortArray.h>

// vtk wrappers
#include <ttkUncertainDataEstimator.h>

class Editor : virtual public ttk::Debug {

public:
  Editor();

  int execute();

  vtkDataSet *getData() const {
    return input_;
  };

  int init(int &argc, char **argv);
  int initRawReader();

  int saveData() const;

  int loadData(const std::string &fileName);

  template <class TReader>
  vtkDataSet *readXMLFile(const std::string &fileName) const {
    vtkSmartPointer<TReader> reader = vtkSmartPointer<TReader>::New();
    reader->SetFileName(fileName.c_str());
    reader->Update();
    reader->GetOutput()->Register(reader);
    return vtkDataSet::SafeDownCast(reader->GetOutput());
  }

  template <class dataType>
  int computeBounds() {
    if(inputFileName_.empty()) {
      return -1;
    }

    // First file in the list as reference
    loadData(inputFileName_[0]);

    // Base code object
    ttk::PDFBounds<dataType> pdfBounds;
    int numberOfVertices = input_->GetNumberOfPoints();
    pdfBounds.setNumberOfVertices(numberOfVertices);
    pdfBounds.evaluateRealization(
      input_->GetPointData()->GetArray(0)->GetVoidPointer(0));

    // Create upper and lower bound arrays
    vtkSmartPointer<vtkDataArray> lowerBound;
    vtkSmartPointer<vtkDataArray> upperBound;
    if(scalarType_ == "float") {
      lowerBound = vtkSmartPointer<vtkFloatArray>::New();
      upperBound = vtkSmartPointer<vtkFloatArray>::New();
    } else if(scalarType_ == "double") {
      lowerBound = vtkSmartPointer<vtkDoubleArray>::New();
      upperBound = vtkSmartPointer<vtkDoubleArray>::New();
    } else if(scalarType_ == "int") {
      lowerBound = vtkSmartPointer<vtkIntArray>::New();
      upperBound = vtkSmartPointer<vtkIntArray>::New();
    } else if(scalarType_ == "unsignedInt") {
      lowerBound = vtkSmartPointer<vtkUnsignedIntArray>::New();
      upperBound = vtkSmartPointer<vtkUnsignedIntArray>::New();
    } else if(scalarType_ == "short") {
      lowerBound = vtkSmartPointer<vtkShortArray>::New();
      upperBound = vtkSmartPointer<vtkShortArray>::New();
    } else if(scalarType_ == "unsignedShort") {
      lowerBound = vtkSmartPointer<vtkUnsignedShortArray>::New();
      upperBound = vtkSmartPointer<vtkUnsignedShortArray>::New();
    } else if(scalarType_ == "char") {
      lowerBound = vtkSmartPointer<vtkCharArray>::New();
      upperBound = vtkSmartPointer<vtkCharArray>::New();
    } else if(scalarType_ == "signedChar") {
      lowerBound = vtkSmartPointer<vtkSignedCharArray>::New();
      upperBound = vtkSmartPointer<vtkSignedCharArray>::New();
    } else if(scalarType_ == "unsignedChar") {
      lowerBound = vtkSmartPointer<vtkUnsignedCharArray>::New();
      upperBound = vtkSmartPointer<vtkUnsignedCharArray>::New();
    }

    lowerBound->SetName("lowerBound");
    upperBound->SetName("upperBound");

    lowerBound->SetNumberOfTuples(input_->GetNumberOfPoints());
    upperBound->SetNumberOfTuples(input_->GetNumberOfPoints());

    // Remove all point arrays from last input before deep copy
    for(int i = 0; i < input_->GetPointData()->GetNumberOfArrays(); i++) {
      input_->GetPointData()->RemoveArray(i);
    }

    // Deep copy
    outputBounds_->DeepCopy(input_);

    // Add bounds arrays
    outputBounds_->GetPointData()->AddArray(lowerBound);
    outputBounds_->GetPointData()->AddArray(upperBound);

    outputBounds_->Register(lowerBound);
    outputBounds_->Register(upperBound);

    // Execute for all files
    for(size_t i = 1; i < inputFileName_.size(); i++) {
      loadData(inputFileName_[i]);
      pdfBounds.evaluateRealization(
        input_->GetPointData()->GetArray(0)->GetVoidPointer(0));
    }

    // Copy the results in the arrays
    for(int i = 0; i < numberOfVertices; i++) {
      lowerBound->SetTuple1(
        i, static_cast<double>((pdfBounds.getLowerBoundPointer())[i]));
      upperBound->SetTuple1(
        i, static_cast<double>((pdfBounds.getUpperBoundPointer())[i]));
    }

    return 0;
  }

protected:
  // Input
  std::string inputDirectory_{};
  std::string inputFormat_{};
  std::vector<std::string> inputFileName_{};
  int numberOfFractions_{};
  int fraction_{};

  // Histogram parameters
  int numberOfBins_{};

  // Parameters for raw files
  bool isBigEndian_{};
  bool isLowerLeft_{};
  std::string scalarType_{};
  int dimension_{};
  int nx_{}, ny_{}, nz_{}; // Dimensions
  double dx_{}, dy_{}, dz_{}; // Spacing
  double ox_{}, oy_{}, oz_{}; // Origin
  // Reader for raw files
  vtkNew<vtkImageReader2> rawReader_{};

  // Input data set
  vtkDataSet *input_{};

  // Output
  std::string outputDirectory_{};
  std::string outputPrefix_{};

  vtkSmartPointer<vtkDataSet> outputBounds_{};
  vtkSmartPointer<vtkDataSet> outputHistograms_{};

  // vtkXMLImageDataReader         *reader_;
  // vtkUncertainDataEstimator     *uncertainDataEstimator_;
};
