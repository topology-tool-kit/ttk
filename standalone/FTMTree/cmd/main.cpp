/*
 * file:                  main.cpp
 * description:           dummy program example.
 * author:                Your Name Here <Your Email Address Here>.
 * date:                  The Date Here.
 */

// include the local headers
#include                  "Editor.h"

int main(int argc, char **argv) {

  // init editor
  Editor editor;
  editor.init(argc, argv);

  // execute data processing
  editor.execute();

  // save the output
  //editor.saveData("output.vti");

  return 0;
}
