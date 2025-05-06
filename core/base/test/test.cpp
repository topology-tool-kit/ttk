#include <test.h>

ttk::test::test() {
  // inherited from Debug: prefix will be printed at the beginning of every msg
  this->setDebugMsgPrefix("test");
}
