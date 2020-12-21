/// \author Michael Michaux <michaelmichaux89@gmail.com>.
/// \date May 2016.
///
/// \brief Imports multiple data-sets into one uncertain data-set.

// include the local headers
#include "Editor.h"

int main(int argc, char **argv) {

  // init editor
  Editor editor;
  editor.init(argc, argv);
  // editor.test();

  // execute data processing
  editor.execute();

  // save the output
  editor.saveData();

  return 0;
}
