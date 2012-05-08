#include "eicsmear/erhic/EventDpmjet.h"

#include <sstream>

namespace erhic {

   bool EventDpmjet::Parse(const std::string& line) {
      static std::stringstream ss;
      ss.str( "" );
      ss.clear();
      ss << line;
      ss >> I >> ievent >> process1 >> process2 >> IP >> dtrueW2 
      >> dtrueNu >> dtrueQ2 >> dtrueX >> dtrueY >> theta_Evt >>photonFlux 
      >> tgtparton >> prjparton >> xtgtparton >> xprjparton >> nucleon
      >> nTracks;
      number=ievent;
      return not ss.fail();
   }
} // namespace erhic
