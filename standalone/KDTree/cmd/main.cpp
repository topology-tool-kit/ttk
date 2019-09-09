

/// \author Your Name Here <Your Email Address Here>.
/// \date The Date Here.
///
/// \brief dummy program example.

#include <KDTree.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <ttkBlank.h>
#include <ttkProgramBase.h>
#include <vector>

using namespace std;
using namespace ttk;

int main(int argc, char **argv) {

  std::cout << "Hello world" << std::endl;

  std::vector<float> coordinates;
  // A
  coordinates.push_back(0);
  coordinates.push_back(0);
  // B
  coordinates.push_back(1);
  coordinates.push_back(1);
  // C
  coordinates.push_back(2);
  coordinates.push_back(-1);
  // D
  coordinates.push_back(2.5);
  coordinates.push_back(-1.5);
  // E
  coordinates.push_back(-2);
  coordinates.push_back(1);
  // F
  coordinates.push_back(-1);
  coordinates.push_back(-1);
  // G
  coordinates.push_back(-0.5);
  coordinates.push_back(-3);

  KDTree<float> kdt;
  kdt.build(coordinates.data(), 7, 2);

  std::cout << "			" << kdt.id_ << std::endl;
  std::cout << "		" << kdt.left_->id_ << "		" << kdt.right_->id_
            << std::endl;
  std::cout << "	   " << kdt.left_->left_->id_ << "	   "
            << kdt.left_->right_->id_ << "	   " << kdt.right_->left_->id_
            << "	   " << kdt.right_->right_->id_ << std::endl;

  KDTree<float> A, B, C, D, E, F, G;
  A = kdt;
  F = *(A.left_);
  C = *(A.right_);
  E = *(F.right_);
  G = *(F.left_);
  B = *(C.right_);
  D = *(C.left_);

  std::cout << A.coords_min_[0] << ", " << A.coordinates_[0] << ", "
            << A.coords_max_[0] << std::endl;
  std::cout << B.coords_min_[0] << ", " << B.coordinates_[0] << ", "
            << B.coords_max_[0] << std::endl;
  std::cout << C.coords_min_[0] << ", " << C.coordinates_[0] << ", "
            << C.coords_max_[0] << std::endl;
  std::cout << D.coords_min_[0] << ", " << D.coordinates_[0] << ", "
            << D.coords_max_[0] << std::endl;
  std::cout << E.coords_min_[0] << ", " << E.coordinates_[0] << ", "
            << E.coords_max_[0] << std::endl;
  std::cout << F.coords_min_[0] << ", " << F.coordinates_[0] << ", "
            << F.coords_max_[0] << std::endl;
  std::cout << G.coords_min_[0] << ", " << G.coordinates_[0] << ", "
            << G.coords_max_[0] << std::endl;

  std::cout << "Bye Bye world" << std::endl;

  /*vtkProgram<ttkBlank> program;

  // TODO-1:
  // specify local parameters to the TTK module with default values.
  bool someOption = false;
  int someIntegerArgument = -1;
  double someDoubleArgument = -1.0;
  // end of TODO-1

  // TODO-2:
  // register these arguments to the command line parser
  program.parser_.setArgument("D", &someDoubleArgument,
    "Some optional double argument", true);
  program.parser_.setArgument("I", &someIntegerArgument,
    "Some optional integer argument", true);
  program.parser_.setOption("O", &someOption,
    "Some option to enable or disable");
  // end of TODO-2

  int ret = 0;
  ret = program.init(argc, argv);

  if(ret != 0)
    return ret;

  // TODO-3:
  // change here the arguments of the vtkWrapper that you want to update prior
  // to execution.
  program.ttkObject_->SetSomeIntegerArgument(someIntegerArgument);
  program.ttkObject_->SetSomeDoubleArgument(someDoubleArgument);
  program.ttkObject_->SetSomeOption(someOption);
  // end of TODO-3

  // execute data processing
  ret = program.run();

  if(ret != 0)
    return ret;

  // save the output
  // optional TODO-4:
  // if you want a different kind of output, re-implement the function save().
  ret = program.save();
  /// end of optional TODO-4

  return ret;*/
}
