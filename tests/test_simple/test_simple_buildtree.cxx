#include <TSystem.h>
#include <eicsmear/functions.h>
int main(){
  // gSystem->Load("libeicsmear");
  BuildTree("./tests/test_simple/beagle_eD.txt", ".", -1);
  return 0;
}



