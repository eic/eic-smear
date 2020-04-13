R__LOAD_LIBRARY(libeicsmear);
#include "eicsmear/erhic/DisKinematics.h"
#include <eicsmear/functions.h>

void make_tree(std::string filstr){
  erhic::DisKinematics::BoundaryWarning=false;
 
  std::string dirstr = ".";
  std::string inputstr = dirstr + "/" + filstr;
  BuildTree(inputstr,dirstr);
}
