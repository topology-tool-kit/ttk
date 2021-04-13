'use strict';
(function() {

  // https://vtk.org/doc/nightly/html/vtkType_8h.html

  const CONSTS = {
    VTK_VOID: 0,
    VTK_BIT: 1,
    VTK_CHAR: 2,
    VTK_SIGNED_CHAR: 15,
    VTK_UNSIGNED_CHAR: 3,
    VTK_SHORT: 4,
    VTK_UNSIGNED_SHORT: 5,
    VTK_INT: 6,
    VTK_UNSIGNED_INT: 7,
    VTK_LONG: 8,
    VTK_UNSIGNED_LONG: 9,
    VTK_FLOAT: 10,
    VTK_DOUBLE: 11,
    VTK_ID_TYPE: 12,
    VTK_STRING: 13
  };

  // Base Class
  class Base {
    constructor() {
      this.listenersPerEvent = new Map();
    }
    on(event, callback) {
      if (!this.listenersPerEvent.has(event))
        this.listenersPerEvent.set(event, []);
      let listeners = this.listenersPerEvent.get(event);
      listeners.push(callback);
    }
    off(event, callback) {
      if (!this.listenersPerEvent.has(event))
        return;

      let listeners = this.listenersPerEvent.get(event);
      for (let i = listeners.length-1; i >= 0; i--)
        if (listeners[i] === callback)
        listeners.splice(i, 1);
    }
    trigger(event, data) {
      if (!this.listenersPerEvent.has(event))
        return;

      let listeners = this.listenersPerEvent.get(event);
      for (let i = 0; i < listeners.length; i++)
        listeners[i](data);
    }
  }

  class vtkDataSet {
    constructor(className, options) {
      this.className = className || "vtkUnstructuredGrid";

      this.points = {
        coordinates: new vtkArray("coordinates", 0, 3, CONSTS.VTK_FLOAT, [])
      };

      this.cells = {
        offsetsArray: new vtkArray("offsetsArray", 1, 1, CONSTS.VTK_INT, [0]),
        connectivityArray: new vtkArray("connectivityArray", 0, 1, CONSTS.VTK_INT, [])
      };

      this.dimension = null;
      this.origin = null;
      this.spacing = null;

      this.pointData = {};
      this.cellData = {};
      this.fieldData = {};

      options = options || {};
      for (let key of Object.keys(options))
        this[key] = options[key];
    }

    toJSON() {
      let json = {};
      json.className = this.className;

      for (let association of ['points', 'cells', 'pointData', 'cellData', 'fieldData']) {
        json[association] = {};
        for (let name of Object.keys(this[association])) {
          json[association][name] = this[association][name].toJSON();
        }
      }

      return json;
    }

    toString() {
      return JSON.stringify(this.toJSON());
    }

    static async createFromMessageSequence(msgs) {
      const objects = [];

      let object = null;
      for (let i = 0; i < msgs.length; i++) {
        const msg = msgs[i];

        if (typeof msg.data === 'string') {
          const msgJSON = JSON.parse(msg.data);
          if (msgJSON.hasOwnProperty('className')) {
            object = new vtkDataSet();
            object.className = msgJSON.className;
            objects.push(object);
          } else if (msgJSON.hasOwnProperty('dimension')) {
            object.dimension = msgJSON.dimension;
            object.origin = msgJSON.origin;
            object.spacing = msgJSON.spacing;
          } else if (msgJSON.hasOwnProperty('data')) {
            const header = msgJSON;
            object[header.target][header.name] = new vtkArray(
              header.name,
              header.nTuples,
              header.nComponents,
              header.dataType,
              header.data.split(',')
            );
          }
        } else if (msg.data instanceof Blob) {
          const header = JSON.parse(msgs[i-1].data); // header is the previous message
          const array = await vtkArray.createFromBinaryData(header, msg.data);
          object[header.target][header.name] = array;
        }
      }

      return objects;
    }
  }

  class vtkArray {
    constructor(name, nTuples, nComponents, dataType, data) {
      this.name = name;
      this.nComponents = nComponents;
      this.nTuples = nTuples;
      this.dataType = dataType;
      this.data = data;
    }

    toJSON() {
      const json = {};
      json.name = this.name;
      json.nComponents = this.nComponents;
      json.nTuples = this.nTuples;
      json.dataType = this.dataType;
      json.data = this.data.join(',');
      return json;
    }

    toString() {
      return JSON.stringify(this.toJSON());
    }

    static async createFromBinaryData(header, binaryData) {
      const array = new vtkArray(
        header.name,
        header.nTuples,
        header.nComponents,
        header.dataType
      );

      const arrayBuffer = await binaryData.arrayBuffer();
      switch (header.dataType) {
        case CONSTS.VTK_FLOAT:
          array.data = new Float32Array(arrayBuffer, 0, binaryData.size/4);
          break;
        case CONSTS.VTK_DOUBLE:
          array.data = new Float64Array(arrayBuffer, 0, binaryData.size/8);
          break;
        case CONSTS.VTK_LONG:
        case CONSTS.VTK_ID_TYPE:
          array.data = new BigInt64Array(arrayBuffer, 0, binaryData.size/8);
          break;
        case CONSTS.VTK_CHAR:
        case CONSTS.VTK_SIGNED_CHAR:
          array.data = new Int8Array(arrayBuffer, 0, binaryData.size);
          break;
        case CONSTS.VTK_UNSIGNED_CHAR:
          array.data = new Uint8Array(arrayBuffer, 0, binaryData.size);
          break;
        case CONSTS.VTK_INT:
          array.data = new Int32Array(arrayBuffer, 0, binaryData.size/4);
          break;
        case CONSTS.VTK_UNSIGNED_INT:
          array.data = new Uint32Array(arrayBuffer, 0, binaryData.size/4);
          break;
        default:
          console.error('Unsupported Data Type');
          return null;
      }

      return array;
    }
  }

  class ttkWebSocketIO extends Base {

    constructor() {
      super();
      this.socket = null;
      this.ip = null;
      this.port = null;

      this.messageSequence = null;
    }

    connect(ip, port) {
      this.ip = ip || 'localhost';
      this.port = port || 9285;

      // close existing connection
      if (this.socket) {
        this.socket.addEventListener('close', ()=>this.connect(ip, port));
        this.socket.close();
        return;
      }

      // create new connection
      const socket = new WebSocket('ws://'+this.ip+':' + this.port);

      // attach internal listeners
      socket.addEventListener('open', e=> {
        this.socket = socket;
        this.trigger('open', e);
      });
      socket.addEventListener('error', e=> {
        this.trigger('error', e);
      });
      socket.addEventListener('close', e=> {
        this.socket = null;
        this.trigger('close', e);
      });
      socket.addEventListener('message', msg=> {
        // check if msg is internal message
        if (typeof msg.data === 'string' && msg.data.startsWith('ttk_WSIO_')) {
          switch (msg.data) {
            case 'ttk_WSIO_BeginMessageSequence':
              this.messageSequence = [];
              this.sendString('ttk_WSIO_RequestNextMessage');
              break;
            case 'ttk_WSIO_EndMessageSequence':
              const messageSequence = this.messageSequence;
              this.messageSequence = null;
              this.trigger('messageSequence', messageSequence);
              break;
          }
        } else {
          if (this.messageSequence !== null) {
            this.messageSequence.push(msg);
            this.sendString('ttk_WSIO_RequestNextMessage');
          } else {
            this.trigger('message', msg);
          }
        }
      });
    }

    sendString(msg) {
      if (this.socket)
        this.socket.send(msg);
      else {
        console.error('No Socket Connection Established');
        this.trigger('error', 'No Socket Connection Established');
      }
    }

    sendVTKDataSet(object) {
      if (!(object instanceof vtkDataSet)) {
        console.error("Input is not a vtkDataSet");
        return 0;
      }

      this.sendString(JSON.stringify({
        "vtkDataSet": object.toJSON()}));
    }
  }

  if (!window.hasOwnProperty('TTK'))
    window.TTK = {};

  window.TTK.CONSTS = CONSTS;
  window.TTK.Base = Base;
  window.TTK.vtkArray = vtkArray;
  window.TTK.vtkDataSet = vtkDataSet;
  window.TTK.ttkWebSocketIO = ttkWebSocketIO;

})();