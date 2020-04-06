#include <TSystem.h>
#include "eicsmear/functions.h"
int main(){
  // gSystem->Load("libeicsmear");
  BuildTree("./tests/beagle_eD.txt", ".", -1);
  return 0;
}



