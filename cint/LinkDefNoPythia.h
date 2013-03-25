#ifdef __CINT__

// General preamble

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ nestedclass;
#pragma link C++ nestedtypedef;

// namespaces

#pragma link C++ namespace erhic;

// Functions

#pragma link C++ function BuildTree;

// Particle classes

#pragma link C++ class erhic::VirtualParticle+;
#pragma link C++ class erhic::ParticleMC+;

// Base event classes

#pragma link C++ class erhic::EventMC+;
#pragma link C++ class erhic::VirtualEvent+;
#pragma link C++ class erhic::EventDis+;

// Event classes for individual generators

#pragma link C++ class erhic::EventPythia+;
#pragma link C++ class erhic::EventRapgap+;
#pragma link C++ class erhic::EventPepsi+;
#pragma link C++ class erhic::EventDjangoh+;
#pragma link C++ class erhic::EventDpmjet+;
#pragma link C++ class erhic::EventMilou+;
#pragma link C++ class erhic::EventGmcTrans+;

// Event building

#pragma link C++ class erhic::VirtualEventFactory+;
#pragma link C++ class erhic::EventFromAsciiFactory<erhic::EventPepsi>+;
#pragma link C++ class erhic::EventFromAsciiFactory<erhic::EventMilou>+;
#pragma link C++ class erhic::EventFromAsciiFactory<erhic::EventRapgap>+;
#pragma link C++ class erhic::EventFromAsciiFactory<erhic::EventDjangoh>+;
#pragma link C++ class erhic::EventFromAsciiFactory<erhic::EventDpmjet>+;
#pragma link C++ class erhic::EventFromAsciiFactory<erhic::EventPythia>+;
#pragma link C++ class erhic::EventFromAsciiFactory<erhic::EventGmcTrans>+;
#pragma link C++ class erhic::EventMCFilterABC+;

// Tree and file building

#pragma link C++ class erhic::Forester+;
#pragma link C++ class erhic::Forester::Status+;

// Monte carlo log file processing

#pragma link C++ class erhic::LogReader+;
#pragma link C++ class erhic::LogReaderPythia+;
#pragma link C++ class erhic::LogReaderPepsi+;
#pragma link C++ class erhic::LogReaderDjangoh+;
#pragma link C++ class erhic::LogReaderMilou+;
#pragma link C++ class erhic::LogReaderGmcTrans+;
#pragma link C++ class erhic::LogReaderFactory;

// Monte carlo file type information

#pragma link C++ class erhic::FileType+;
#pragma link C++ class erhic::File<erhic::EventPythia>+;
#pragma link C++ class erhic::File<erhic::EventMilou>+;
#pragma link C++ class erhic::File<erhic::EventPepsi>+;
#pragma link C++ class erhic::File<erhic::EventRapgap>+;
#pragma link C++ class erhic::File<erhic::EventDjangoh>+;
#pragma link C++ class erhic::File<erhic::EventDpmjet>+;
#pragma link C++ class erhic::File<erhic::EventGmcTrans>+;
#pragma link C++ class erhic::FileFactory;

// Specialised stl templates

#pragma link C++ class std::vector<Particle>+;
#pragma link C++ class std::vector<Particle*>+;
#pragma link C++ class std::vector<const Particle*>+;
#pragma link C++ class std::vector<const erhic::VirtualParticle*>+;
#pragma link C++ class std::auto_ptr<erhic::Pid>;

// Miscellaneous utilities and helper functions/functors

#pragma link C++ class erhic::DisKinematics+;
#pragma link C++ class BeamParticles+;
#pragma link C++ class EventToDot;
#pragma link C++ class erhic::Pid+;
#pragma link C++ class erhic::Reader+;
#pragma link C++ class erhic::KinematicsComputer+;
#pragma link C++ class erhic::LeptonKinematicsComputer+;
#pragma link C++ class erhic::JacquetBlondelComputer+;
#pragma link C++ class erhic::DoubleAngleComputer+;
#pragma link C++ class ParticleIdentifier+;

#endif
