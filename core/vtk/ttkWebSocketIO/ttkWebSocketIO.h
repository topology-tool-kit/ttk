/// \ingroup vtk
/// \class ttkWebSocketIO
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 13.04.2021
///
/// This module provides functions for bi-directional communication between two
/// websockets, where the heavy lifting is done by the websocketpp library. The
/// primary use case for this module is to send data that was computed in a c++
/// environment to a browser client, which can then freely process/render the
/// data via JavaScript and HTML. To this end, this module needs to run a
/// websocket server, to which a client can connect to via the ttkWebSocketIO.js
/// library. The server can also receive data from the client, which is then fed
/// into the processEvent function. In order to add custom processing of events,
/// one needs to inherit from the WebSocketIO class and override the
/// processEvent method. The vtk wrapper ttkWebSocketIO provides an example on
/// how this abstract base module can be specialized to the application of
/// sending/receiving vtkDataObjects. Specifically, the ttkWebSocketIO filter
/// will send its input vtkDataObject as a serialized JSON object to all
/// connected clients every time the filter is called with a new input. When the
/// server receives a serialized JSON object form the client, then the filter
/// will instantiate a vtkDataObject and pass it as the filter output.

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
  int ParseVtkDataObjectFromJSON(const std::string &json);
};
