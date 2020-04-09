R__LOAD_LIBRARY(../build/libeicsmear.so);

void make_tree(std::string filstr){

  erhic::DisKinematics::BoundaryWarning=false;
 
  std::string dirstr = ".";
  std::string inputstr = dirstr + "/" + filstr;
  BuildTree(inputstr,dirstr);
}
