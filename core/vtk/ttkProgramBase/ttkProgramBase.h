/// \ingroup vtk
/// \class ttkProgramBase
/// \author Julien Tierny <julien.tierny@lip6.fr>.
/// \date February 2017.
///
/// \brief Base VTK editor class for standalone programs. This class parses the
/// the comamnd line, execute the TTK module and takes care of the IO.

#ifndef _TTK_EDITOR_BASE_H
#define _TTK_EDITOR_BASE_H

// VTK IO
#include <vtkDataSet.h>
#include <vtkDataSetAlgorithm.h>
#include <vtkImageData.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLImageDataReader.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridWriter.h>

// base code includes
#include <ProgramBase.h>
#include <ttkWrapper.h>

class VTKFILTERSCORE_EXPORT ttkProgramBase : public ttk::ProgramBase {

public:
  ttkProgramBase() {

    vtkWrapper_ = NULL;
  };

  virtual ~ttkProgramBase(){};

  /// Set the arguments of your ttk module and execute it here.
  int execute();

  vtkDataSet *getInput(const int &inputId) {
    if((inputId < 0) || (inputId >= (int)inputs_.size()))
      return NULL;
    return inputs_[inputId];
  }

  int getNumberOfInputs() {
    return inputs_.size();
  };

  virtual int run() {

    if(!vtkWrapper_) {
      return -1;
    }

    return execute();
  }

  /// Save the output(s) of the TTK module.
  virtual int save() const;

  virtual int setTTKmodule(vtkDataSetAlgorithm *ttkModule) {

    vtkWrapper_ = ttkModule;
    ttkModule_ = (Debug *)ttkModule;

    return 0;
  }

protected:
  std::vector<vtkDataSet *> inputs_;
  std::vector<vtkSmartPointer<vtkXMLImageDataReader>> imageDataReaders_;
  std::vector<vtkSmartPointer<vtkXMLPolyDataReader>> polyDataReaders_;
  std::vector<vtkSmartPointer<vtkXMLUnstructuredGridReader>>
    unstructuredGridReaders_;
  vtkDataSetAlgorithm *vtkWrapper_;

  template <class vtkReaderClass>
  int load(const std::string &fileName,
           std::vector<vtkSmartPointer<vtkReaderClass>> &readerList);

  /// Load a sequence of input data-sets.
  virtual int load(const std::vector<std::string> &inputPaths);

  template <class vtkWriterClass>
  int save(const int &outputPortId) const;
};

template <class ttkModule>
class vtkProgram : public ttkProgramBase {

public:
  vtkProgram() {
    ttkObject_ = vtkSmartPointer<ttkModule>::New();
    vtkWrapper_ = (vtkDataSetAlgorithm *)ttkObject_.GetPointer();
    ttkModule_ = (Debug *)ttkObject_.GetPointer();
  }

  virtual int run() {

    ttkObject_->setDebugLevel(ttk::globalDebugLevel_);
    ttkObject_->setThreadNumber(ttk::globalThreadNumber_);

    return ttkProgramBase::run();
  }

  vtkSmartPointer<ttkModule> ttkObject_;
};

template <class vtkWriterClass>
int ttkProgramBase::save(const int &outputPortId) const {

  if(!vtkWrapper_)
    return -1;

  std::string extension;

  if((vtkWrapper_->GetOutput(outputPortId)->GetDataObjectType()
      == VTK_IMAGE_DATA)) {
    //     ||(vtkWrapper_->GetOutput(outputPortId)->GetDataObjectType()
    //       == TTK_IMAGE_DATA)){
    extension = "vti";
  }

  if((vtkWrapper_->GetOutput(outputPortId)->GetDataObjectType()
      == VTK_POLY_DATA)) {
    //     ||(vtkWrapper_->GetOutput(outputPortId)->GetDataObjectType()
    //       == TTK_POLY_DATA)){
    extension = "vtp";
  }

  if((vtkWrapper_->GetOutput(outputPortId)->GetDataObjectType()
      == VTK_UNSTRUCTURED_GRID)) {
    //     ||(vtkWrapper_->GetOutput(outputPortId)->GetDataObjectType()
    //       == TTK_UNSTRUCTURED_GRID)){
    extension = "vtu";
  }

  std::stringstream fileName;
  fileName << outputPath_ << "_port#" << outputPortId << "." << extension;

  vtkSmartPointer<vtkWriterClass> writer
    = vtkSmartPointer<vtkWriterClass>::New();
  writer->SetFileName(fileName.str().data());
  writer->SetInputData(vtkWrapper_->GetOutput(outputPortId));
  std::stringstream msg;
  msg << "[ttkProgramBase] Saving output file `" << fileName.str() << "'..."
      << std::endl;
  dMsg(std::cout, msg.str(), Debug::infoMsg);

  writer->Write();

  return 0;
}

template <class vtkReaderClass>
int ttkProgramBase::load(
  const std::string &fileName,
  std::vector<vtkSmartPointer<vtkReaderClass>> &readerList) {

  readerList.resize(readerList.size() + 1);
  readerList.back() = vtkSmartPointer<vtkReaderClass>::New();

  readerList.back()->SetFileName(fileName.data());

  // handle debug messages
  {
    std::stringstream msg;
    msg << "[ttkProgramBase] Reading input data..." << std::endl;
    // choose where to display this message (std::cout, std::cerr, a file)
    // choose the priority of this message (1, nearly always displayed,
    // higher values mean lower priorities)
    dMsg(std::cout, msg.str(), 1);
  }

  readerList.back()->Update();
  inputs_.push_back(readerList.back()->GetOutput());

  if(!inputs_.back())
    return -1;

  if(!inputs_.back()->GetNumberOfPoints())
    return -2;

  if(!inputs_.back()->GetNumberOfCells())
    return -3;

  {
    std::stringstream msg;
    msg << "[ttkProgramBase]   done! (read "
        << inputs_.back()->GetNumberOfPoints() << " vertices, "
        << inputs_.back()->GetNumberOfCells() << " cells)" << std::endl;
    dMsg(std::cout, msg.str(), Debug::infoMsg);
  }

  return 0;
}

#endif // VTK_EDITOR_BASE_H
