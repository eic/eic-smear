#ifdef __CINT__

// General preamble

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ nestedclass;
#pragma link C++ nestedtypedef;

// namespaces

#pragma link C++ namespace erhic;
#pragma link C++ namespace Smear;

// Functions

#pragma link C++ function BuildTree;
#pragma link C++ function SmearTree;
#pragma link C++ function Smear::ParseInputFunction;

// Particle classes

#pragma link C++ class erhic::VirtualParticle+;
#pragma link C++ class erhic::ParticleMC+;
#pragma link C++ class Smear::ParticleMCS+;

// Base event classes

#pragma link C++ class erhic::VirtualEvent<erhic::ParticleMC>+;
#pragma link C++ class erhic::EventMC+;
#pragma link C++ class erhic::VirtualEvent<Smear::ParticleMCS>+;
#pragma link C++ class Smear::Event+;

// Event classes for individual generators

#pragma link C++ class erhic::EventPythia+;
#pragma link C++ class EventRapgap+;
#pragma link C++ class EventPepsi+;
#pragma link C++ class EventDjangoh+;
#pragma link C++ class EventDpmjet+;
#pragma link C++ class EventMilou+;

// Event building

#pragma link C++ class erhic::VirtualEventFactory+;
#pragma link C++ class erhic::EventFromAsciiFactory<EventPepsi>+;
#pragma link C++ class erhic::EventFromAsciiFactory<EventMilou>+;
#pragma link C++ class erhic::EventFromAsciiFactory<EventRapgap>+;
#pragma link C++ class erhic::EventFromAsciiFactory<EventDjangoh>+;
#pragma link C++ class erhic::EventFromAsciiFactory<EventDpmjet>+;
#pragma link C++ class erhic::EventFromAsciiFactory<erhic::EventPythia>+;
#pragma link C++ class erhic::Pythia6EventBuilder+;
#pragma link C++ class erhic::EventMCFilterABC+;
#pragma link C++ class erhic::Pythia6+;

// Tree and file building

#pragma link C++ class erhic::Forester+;
#pragma link C++ class erhic::Forester::Status+;

// Core smearing components

#pragma link C++ class Smear::Acceptance+;
#pragma link C++ class Smear::Smearer+;
#pragma link C++ class Smear::Device+;
#pragma link C++ class Smear::ParamSimple+;
#pragma link C++ class Smear::Detector+;
#pragma link C++ class Smear::Distributor+;
#pragma link C++ class Smear::EventKinematicsComputer+;
#pragma link C++ class Smear::ParticleID+;

// Specialized smearing devices

#pragma link C++ class Smear::Devious+;
#pragma link C++ class Smear::Tracker+;
#pragma link C++ Class Smear::EMCalorimeter+;
#pragma link C++ class Smear::HCalorimeter+;
#pragma link C++ class Smear::Bremsstrahlung+;

// Monte carlo log file processing

#pragma link C++ class erhic::LogReader+;
#pragma link C++ class erhic::LogReaderPythia+;
#pragma link C++ class erhic::LogReaderPepsi+;
#pragma link C++ class erhic::LogReaderDjangoh+;
#pragma link C++ class erhic::LogReaderMilou+;
#pragma link C++ class erhic::LogReaderFactory;

// Monte carlo file type information

#pragma link C++ class erhic::FileType+;
#pragma link C++ class erhic::File<erhic::EventPythia>+;
#pragma link C++ class erhic::File<EventMilou>+;
#pragma link C++ class erhic::File<EventPepsi>+;
#pragma link C++ class erhic::File<EventRapgap>+;
#pragma link C++ class erhic::File<EventDjangoh>+;
#pragma link C++ class erhic::File<EventDpmjet>+;
#pragma link C++ class erhic::FileFactory;

// Specialised stl templates

#pragma link C++ class std::vector<Particle>+;
#pragma link C++ class std::vector<Particle*>+;
#pragma link C++ class std::vector<const Particle*>+;
#pragma link C++ class std::auto_ptr<erhic::Pid>;

// typedefs

#pragma link C++ class EventS;
#pragma link C++ class ParticleS;
#pragma link C++ typedef erhic::VirtualEvent<erhic::ParticleMC>::TrackType;
#pragma link C++ typedef erhic::VirtualEvent<Smear::ParticleMCS>::TrackType;

// Miscellaneous utilities and helper functions/functors

#pragma link C++ class BeamParticles+;
#pragma link C++ class DoubleAngle+;
#pragma link C++ class EventToDot;
#pragma link C++ class erhic::Pid+;
#pragma link C++ class erhic::Reader+;
#pragma link C++ class KinematicsFromHadrons+;
#pragma link C++ class JacquetBlondel+;
#pragma link C++ class ParticleIdentifier+;

#endif
