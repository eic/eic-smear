#ifdef __CINT__

#pragma link C++ nestedclass;
#pragma link C++ nestedtypedef;

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ function BuildTree;
#pragma link C++ function SmearTree;

//#pragma link C++ class Particle+;

//#pragma link C++ class EventBase+;
#pragma link C++ class EventPythia+;
#pragma link C++ class EventRapgap+;
#pragma link C++ class EventPepsi+;
#pragma link C++ class EventDjangoh+;
#pragma link C++ class EventMilou+;

#pragma link C++ class KinematicsFromHadrons+;
#pragma link C++ class JacquetBlondel+;
#pragma link C++ class DoubleAngle+;
//#pragma link C++


#pragma link C++ class ParticleIdentifier+;

#pragma link C++ namespace Smear;
#pragma link C++ class Smear::Device+;
#pragma link C++ class Smear::Detector+;
#pragma link C++ class Smear::ParticleID+;
#pragma link C++ class Smear::Acceptance+;
#pragma link C++ class Smear::Distributor+;
#pragma link C++ class Smear::EventKinematicsComputer+;

//specialized devices
#pragma link C++ class Smear::Devious+;
#pragma link C++ class Smear::Tracker+;
#pragma link C++ Class Smear::EMCalorimeter+;
#pragma link C++ class Smear::HCalorimeter+;

#pragma link C++ namespace erhic;
//#pragma link C++ class erhic::Smearer+;
#pragma link C++ class erhic::Forester+;
#pragma link C++ class erhic::Forester::Status+;
#pragma link C++ class BeamParticles+;

#pragma link C++ namespace erhic;

//#pragma link C++ class erhic::File<int>+;
//#pragma link C++ class erhic::File<double>+;
//#pragma link C++ class erhic::GeneratorType+;
#pragma link C++ class erhic::FileType+;
#pragma link C++ class erhic::File<EventPythia>+;
#pragma link C++ class erhic::File<EventMilou>+;
#pragma link C++ class erhic::File<EventPepsi>+;
#pragma link C++ class erhic::File<EventRapgap>+;
#pragma link C++ class erhic::File<EventDjangoh>+;
#pragma link C++ class erhic::FileFactory;
#pragma link C++ class erhic::LogReader+;
#pragma link C++ class erhic::LogReaderPythia+;
#pragma link C++ class erhic::LogReaderMilou+;
#pragma link C++ class erhic::LogReaderFactory;

#pragma link C++ class std::vector<Particle>+;
#pragma link C++ class std::vector<Particle*>+;
//#pragma link C++ class std::vector<erhic::ParticleMCS*>+;
#pragma link C++ class std::vector<const Particle*>+;
//#pragma link C++ class std::vector<ParticleS*>+;

//#pragma link C++ class ParticleRef;

//#pragma link C++ class erhic::KinematicsCalculation;
//#pragma link C++ class erhic::LeptonKinematics+;
//#pragma link C++ class erhic::JacquetBlondel+;
//#pragma link C++ class erhic::DoubleAngle+;

#pragma link C++ class Smear::Bremsstrahlung+;
#pragma link C++ class Smear::Device::Function+;
#pragma link C++ class Smear::Copier+;
#pragma link C++ class Smear::Gauss+;
//#pragma link C++ class EventBase::ParticleRef;

// Add additional classes that you want to use interactively in ROOT here:
//#pragma link C++ class MyClassName+;

//#pragma link C++ class Smearing+;
//#pragma link C++ class SmearEnergy+;
//#pragma link C++ class SmearPTotal+;
//#pragma link C++ class SmearPTrans+;
//#pragma link C++ class SmearPhi+;
//#pragma link C++ class SmearTheta+;
//#pragma link C++ class SmearSet+;
//#pragma link C++ typedef Smearing::Variable;

//#pragma link C++ class PRef;

// Testing of new base classes

#pragma link C++ class erhic::VirtualParticle+;
#pragma link C++ class erhic::ParticleMC+;
#pragma link C++ class Smear::ParticleMCS+;
//#pragma link C++ class erhic::Event<Particle>+;
#pragma link C++ class erhic::VirtualEvent<erhic::ParticleMC>+;
#pragma link C++ typedef erhic::VirtualEvent<erhic::ParticleMC>::TrackType;
#pragma link C++ class erhic::EventMC+;
#pragma link C++ class erhic::VirtualEvent<Smear::ParticleMCS>+;
#pragma link C++ typedef erhic::VirtualEvent<Smear::ParticleMCS>::TrackType;
#pragma link C++ class Smear::Event+;

#pragma link C++ class EventS;
#pragma link C++ class ParticleS;

#pragma link C++ class EventToDot;

#pragma link C++ class erhic::VirtualEventFactory+;
#pragma link C++ class erhic::EventFromAsciiFactory<EventPepsi>+;
#pragma link C++ class erhic::EventFromAsciiFactory<EventMilou>+;
#pragma link C++ class erhic::EventFromAsciiFactory<EventRapgap>+;
#pragma link C++ class erhic::EventFromAsciiFactory<EventDjangoh>+;
#pragma link C++ class erhic::EventFromAsciiFactory<EventPythia>+;

/*
 #pragma link C++ class Particle+;
 #pragma link C++ class Particle::Filter+;
 #pragma link C++ class ParticleMC+;
 #pragma link C++ class ParticleMCS+;
 #pragma link C++ class Event<ParticleMC>+;
 #pragma link C++ class Event<Particle>+;
 #pragma link C++ class EventMC+;
 #pragma link C++ class Event<ParticleMCS>+;
 #pragma link C++ class Event+;
 */

#pragma link C++ class erhic::Pid+;
#pragma link C++ class std::auto_ptr<erhic::Pid>;
//#pragma link C++ class erhic::EventMCFinalState<erhic::ParticleMC>+;
//#pragma link C++ class erhic::EventMCFinalState<Smear::ParticleMCS>+;

#pragma link C++ class erhic::SizeOf+;

#pragma link C++ class erhic::Reader+;

#pragma link C++ class erhic::Pythia6EventBuilder+;
//#pragma link C++ class std::auto_ptr<erhic:EventPythia>+;

//#pragma link C++ class Header+;
#pragma link C++ class erhic::EventMCFilterABC+;
//#pragma link C++ class erhic::EventMCNullFilter+;
//#pragma link C++ class erhic::RandomFractionFilter+;

#endif
