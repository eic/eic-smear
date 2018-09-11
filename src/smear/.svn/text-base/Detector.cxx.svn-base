/**
 \file
 Implementation of class Smear::Detector.
 
 \author    Michael Savastio
 \date      2011-08-19
 \copyright 2011 Brookhaven National Lab
 */

#include "eicsmear/smear/Detector.h"

#include <algorithm>
#include <functional>
#include <list>
#include <memory>
#include <vector>

#include "eicsmear/erhic/EventDis.h"
#include "eicsmear/smear/EventSmear.h"
#include "eicsmear/erhic/Kinematics.h"
#include "eicsmear/smear/ParticleMCS.h"
#include "eicsmear/smear/Smearer.h"
#include "eicsmear/erhic/VirtualParticle.h"

namespace Smear {

Detector::Detector()
: useNM(false)
, useJB(false)
, useDA(false) {
}

Detector::Detector(const Detector& other)
: TObject(other) {
  useNM = other.useNM;
  useJB = other.useJB;
  useDA = other.useDA;
  Devices = other.CopyDevices();
}

Detector& Detector::operator=(const Detector& that) {
  if (this != &that) {
    useNM = that.useNM;
    useJB = that.useJB;
    useDA = that.useDA;
    Devices = that.CopyDevices();
  }  // if
  return *this;
}

Detector::~Detector() {
  DeleteAllDevices();
}

void Detector::DeleteAllDevices() {
  for (unsigned i(0); i < GetNDevices(); i++) {
    delete Devices.at(i);
    Devices.at(i) = NULL;
  }  // for
  Devices.clear();
}

void Detector::AddDevice(Smearer& dev) {
  Devices.push_back(dev.Clone());
}

void Detector::SetEventKinematicsCalculator(TString s) {
  s.ToLower();
  useNM = s.Contains("nm") || s.Contains("null");
  useJB = s.Contains("jb") || s.Contains("jacquet");
  useDA = s.Contains("da") || s.Contains("double");
}

Smearer* Detector::GetDevice(int n) {
  Smearer* smearer(NULL);
  if (unsigned(n) < Devices.size()) {
    smearer = Devices.at(n);
  }  // if
  return smearer;
}

void Detector::FillEventKinematics(Event* eventS) {
  if (!(useNM || useJB || useDA)) {
    return;
  }  // if
     // Need a bit of jiggery-pokery here, as the incident beam info isn't
     // associated with the smeared event.
     // So, get the beam info from the MC event, but replace the scattered
     // electron with the smeared version.
     // Then we can use the standard JB/DA algorithms on the smeared event.
  const ParticleMCS* scattered = eventS->ScatteredLepton();
  typedef std::auto_ptr<erhic::DisKinematics> KinPtr;
  if (useNM && scattered) {
    KinPtr kin(erhic::LeptonKinematicsComputer(*eventS).Calculate());
    if (kin.get()) {
      eventS->SetLeptonKinematics(*kin);
    }  // if
  } else {
    eventS->SetLeptonKinematics(
                                erhic::DisKinematics(-1., -1., -1., -1., -1.));
  }  // if
  if (useJB) {
    KinPtr kin(erhic::JacquetBlondelComputer(*eventS).Calculate());
    if (kin.get()) {
      eventS->SetJacquetBlondelKinematics(*kin);
    }  // if
  }  // if
  if (useDA && scattered) {
    KinPtr kin(erhic::DoubleAngleComputer(*eventS).Calculate());
    if (kin.get()) {
      eventS->SetDoubleAngleKinematics(*kin);
    }  // if
  }  // if
}

std::list<Smearer*> Detector::Accept(const erhic::VirtualParticle& p) const {
  std::list<Smearer*> devices;
  // Only accept final-state particles, so skip the check against each
  // devices for non-final-state particles.
  if (p.GetStatus() == 1) {
    std::vector<Smearer*>::const_iterator iter;
    for (iter = Devices.begin(); iter != Devices.end(); ++iter) {
      // Store each device that accepts the particle.
      if ((*iter)->Accept.Is(p)) {
        devices.push_back(*iter);
      }  // if
    }  // for
  }  // if
  return devices;
}

ParticleMCS* Detector::Smear(const erhic::VirtualParticle& prt) const {
  // Does the particle fall in the acceptance of any device?
  // If so, we smear it, if not, we skip it (store a NULL pointer).
  std::list<Smearer*> devices = Accept(prt);
  ParticleMCS* prtOut(NULL);
  if (!devices.empty()) {
    // It passes through at least one device, so smear it.
    // Devices in which it doesn't pass won't smear it.
    prtOut = new ParticleMCS();
    std::list<Smearer*>::iterator iter;
    for (iter = devices.begin(); iter != devices.end(); ++iter) {
      (*iter)->Smear(prt, *prtOut);
    }  // for
       // Compute derived momentum components.
    prtOut->px = prtOut->p * sin(prtOut->theta) * cos(prtOut->phi);
    prtOut->py = prtOut->p * sin(prtOut->theta) * sin(prtOut->phi);
    prtOut->pt = sqrt(pow(prtOut->px, 2.) + pow(prtOut->py, 2.));
    prtOut->pz = prtOut->p * cos(prtOut->theta);
  }  // if
  return prtOut;
}

std::vector<Smearer*> Detector::CopyDevices() const {
  std::vector<Smearer*> copies;
  std::transform(Devices.begin(), Devices.end(),
                 std::back_inserter(copies),
                 std::bind2nd(std::mem_fun(&Smearer::Clone), ""));
  return copies;
}

void Detector::Print(Option_t* o) const {
  for (unsigned i(0); i < GetNDevices(); ++i) {
    Devices.at(i)->Print(o);
  }  // for
}

}  // namespace Smear
