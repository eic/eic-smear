/**
 \file
 Contains the declaration of class Bremsstrahlung, a specialized Device
 class for modelling radiative losses.
 
 \author Michael Savastio
 */

#ifndef _ERHIC_BUILDTREE_BREMSSTRAHLUNG_
#define _ERHIC_BUILDTREE_BREMSSTRAHLUNG_

#include <memory>

#include <TF1.h>

#include "Device.h"
#include "EventS.h" // For ParticleS
#include "Particle.h" // For Particle

namespace Smear {
   
   struct Bremsstrahlung: public Device {
		
      /**
       Constructor.
       Photon energies are randomly generated in the range
       [epsilon, E - epsilon] for particle energy E.
       traversed is the distance of material through with the particle passes
       and radLength is the radiation length of that material (both in cm).
       */
      Bremsstrahlung(double epsilon = 0.01,
                     double traversed = 10.,
                     double radLength = 47.1);
		
      /**
       Returns dSigmga/dK at k = x[0].
       The arguments have this form to interface with ROOT::TF1.
       The second argument is unused.
       */
		double dSigmadK(double* x, double*);
		
      /** Compute the number of photons emitted */
		int NGamma();
		
		void FixParticleKinematics(ParticleS&);
		
      /** Returns a pointer to a duplicate of this object */
		virtual Bremsstrahlung* Clone();
		
      /** Smear the properties of a Particle and assign them to a ParticleS */
      virtual void DevSmear(const Particle&, ParticleS&);
      
//   protected:
      
      /**
       Set the radiating particle type and configure the dSigma/dK
       function.
       */
		void SetParticle(const Particle&);
      
      /**
       Configure the dSigma/dK function, setting the energy range over
       which to generate photons.
       The energy range is computed from the energy of the current mParticle.
       If the resultant energy range is invalid (e.g. max < min) the range
       is not set and the function returns false.
       */
		bool SetupPDF();
		
      std::auto_ptr<Particle> mParticle; // Copy of the current particle
		
		double mKMin;
		double mKMax;
		double mEpsilon;
		double mTraversed;
		double mRadLength;
		
		TF1* mPdf; // dSigma/dK function
      
		ClassDef(Bremsstrahlung, 1)
	};
   
}

#endif