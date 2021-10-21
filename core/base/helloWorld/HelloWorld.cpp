#include <HelloWorld.h>

ttk::HelloWorld::HelloWorld() {
  // inherited from Debug: prefix will be printed at the beginning of every msg
  this->setDebugMsgPrefix("HelloWorld");
}
