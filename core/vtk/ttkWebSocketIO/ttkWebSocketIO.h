/// \ingroup vtk
/// \class ttkWebSocketIO
/// \author Jonas Lukasczyk <jl@jluk.der>
/// \date 01.09.2019.
///
/// TODO

#pragma once

// VTK Module
#include <ttkWebSocketIOModule.h>

// VTK Includes
#include <ttkAlgorithm.h>
#include <vtkSmartPointer.h>

class vtkDataObject;
class vtkMultiBlockDataSet;

// TTK Base Includes
#include <WebSocketIO.h>

class TTKWEBSOCKETIO_EXPORT ttkWebSocketIO : public ttkAlgorithm,
                                             protected ttk::WebSocketIO {
private:
  int PortNumber{9285};
  bool NeedsUpdate{false};

  vtkSmartPointer<vtkDataObject> LastInput;
  vtkSmartPointer<vtkMultiBlockDataSet> LastOutput;

public:
  static ttkWebSocketIO *New();
  vtkTypeMacro(ttkWebSocketIO, ttkAlgorithm);

  vtkSetMacro(PortNumber, int);
  vtkGetMacro(PortNumber, int);

  vtkSetMacro(NeedsUpdate, bool);
  vtkGetMacro(NeedsUpdate, bool);

  int processEvent(const std::string &eventName,
                   const std::string &eventData = "") override;

protected:
  ttkWebSocketIO();
  ~ttkWebSocketIO();

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

  int SendVtkDataObject(vtkDataObject *object);
  int ParseVtkDataObjectFromJSON(const std::string &json,
                                 vtkDataObject *object);
};