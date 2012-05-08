/**
 hadronicLinkDef.h
 \author Thomas Burton 
 \date 5/3/12
 \copyright 2012 BNL. All rights reserved.
*/

#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclasses;

#pragma link C++ namespace erhic::hadronic;

// Extensions for hadronic events
#pragma link C++ class erhic::hadronic::EventMC+;
#pragma link C++ class erhic::hadronic::EventPythiaPP+;
#pragma link C++ class erhic::hadronic::EventSmear+;
#pragma link C++ class erhic::hadronic::ParticleMC+;
#pragma link C++ class erhic::hadronic::Pythia6EventFactory+;
//#pragma link C++ defined_in "../eic-smear/include/eicsmear/hadronic/EventMC.h"

#endif // __CINT__
