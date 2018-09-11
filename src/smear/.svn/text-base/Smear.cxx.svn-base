/**
 \file
 Imolementation of various Smear namespace functions.
 
 \author    Michael Savastio
 \date      2011-08-19
 \copyright 2011 Brookhaven National Lab
 */

#include "eicsmear/smear/Smear.h"

#include <iostream>

#include <TString.h>

namespace Smear {

int ParseInputFunction(TString &s, KinType &kin1, KinType &kin2) {
  int d = 0;
  if (s.Contains("E")) {
    s.ReplaceAll("E", "x");
    d++;
    kin1 = kE;
  }  // if
  if (s.Contains("P")) {
    d++;
    if (d == 2) {
      s.ReplaceAll("P", "y");
      kin2 = kP;
    } else {
      s.ReplaceAll("P", "x");
      kin1 = kP;
    }  // if
  }  // if
  if (s.Contains("theta") || s.Contains("Theta")) {
    d++;
    if (d == 2) {
      s.ReplaceAll("theta", "y"); s.ReplaceAll("Theta", "y");
      kin2 = kTheta;
    } else {
      s.ReplaceAll("theta", "x"); s.ReplaceAll("Theta", "x");
      kin1 = kTheta;
    }  // if
  }  // if
  if (s.Contains("phi")) {
    d++;
    if (d == 2) {
      s.ReplaceAll("phi", "y");
      kin2 = kPhi;
    } else {
      s.ReplaceAll("phi", "x");
      kin1 = kPhi;
    }  // if
  }  // if
  if (s.Contains("pT")) {
    d++;
    if (d == 2) {
      s.ReplaceAll("pT", "y");
      kin2 = kPt;
    } else {
      s.ReplaceAll("pT", "x");
      kin1 = kPt;
    }  // if
  }  // if
  if (s.Contains("pZ")) {
    d++;
    if (d == 2) {
      s.ReplaceAll("pZ", "y");
      kin2 = kPz;
    } else {
      s.ReplaceAll("pZ", "x");
      kin1 = kPz;
    }  // if
  }  // if
  if (d > 2 || d < 0) {
    std::cout << "!WARNING! Not enough, or too many parameters "
    << "detected in parametrization (d=";
    std::cout << d << ").\n";
  }  // if
  if (d == 0) {
    d++;
  }  // if
  return d;
}

}  // namespace Smear
