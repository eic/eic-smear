#ifndef _ERHIC_BUILDTREE_SMEAREVENT_
#define _ERHIC_BUILDTREE_SMEAREVENT_

#include <TROOT.h>
#include <TSystem.h>
#include <TString.h>
#include <TRandom3.h>
#include <TF1.h>
#include <fstream>
#include <TF2.h>
#include <TLorentzVector.h>
#include <iostream>
#include <sstream>
#include <cmath>

#include "EventBase.h"
#include "VirtualParticle.h"
#include "Kinematics.h"
#include "ParticleIdentifier.h"
#include "ParticleMCSmeared.h"

#include "Smear.h"

namespace Smear {
	
	/**
	 This structure is used to generate event-wise kinematic data from smeared values.  Currently three methods
	 are available.  Null momentum (NM) assumes the momentum of all particles is null, and uses the scattered
	 lepton for all calculations.  Jacquet-Blondel (JB) and Double Angle (DA) methods both use the final state
	 hadronic system only, i.e. they do not need the scattered electron.
	 */
	struct EventKinematicsComputer {
		
		EventKinematicsComputer() {
			
			fault = -999.;
			ErrorValues();
			
			needBackground = true;
			bGoodJohnnieLept = false;
			bGoodJohnnieHad = false;
			bTolerated = false;
			bTolerant = true;
			bMute = false;
			
			HadTol = 100.;
			
			LeptonType = 11;
			PartIdent = ParticleIdentifier(LeptonType);
			
		}
		
		JacquetBlondel JBComp;
//      erhic::JacquetBlondel mJBComp;
		DoubleAngle DAComp;
		
		double HadTol;
		bool needBackground;
		bool bGoodJohnnieLept;
		bool bGoodJohnnieHad;
		bool bTolerated;
		bool bTolerant;
		bool bMute;
		double fault;
		int LeptonType;
		ParticleIdentifier PartIdent;
		
		struct InputKin { 
			
			double Ee;
			double EPrime;
			double ENucleon;
			double MNucleon;
			double thetae;
			
		} InKin;
		
		struct OutputKin {
			double x; double y;
			double QSquared; double WSquared;
		} OutKin;
		
	
		/**
		 Returns x being stored by the computer.
		 */
		double x() {
			return OutKin.x;
		}
		
		/**
		 Returns y being stored by the computer.
		 */
		double y() {
			return OutKin.y;
		}
		
		/** 
		 Returns Q^2 being stored by the computer.
		 */
		double QSquared() {
			return OutKin.QSquared;
		}
		
		/**
		 Returns W^2 being stored by the computer.
		 */
		double WSquared() {
			return OutKin.WSquared;
		}
		
		/**
		 Set the beam lepton type.  Must use PDG particle codes.  Default is e- (11) of course.
		 */
		void SetPDGLeptonCode(int n) {
			LeptonType = n;
			PartIdent = ParticleIdentifier(LeptonType);
		}
		
		/**
		 If the argument is false, then ReadEvent will not read in any of the hadronic system.
		 Obviously this would invalidate the JB and DA methods.  This is useful if you are reading in
		 an enormous number of events and only need NM.
		 */
		void SetNeedBackground(bool d) {
			needBackground = d;
		}
		
		void SetSupressWarnings(bool b) {
			bMute = b;
		}
		
		/**
		 Set a default value which the EventKinematicsComputer will store if it is told to make a computation
		 based on incomplete data.  By default, this is -999.
		 */
		void SetDefaultErrorValue(double v) {
			fault = v;
		}
		
		/**
		 Determines if the ParticleS has had all of E, p, theta and phi smeared to sensible values.  
		 */
		bool isAccepted(ParticleS *prt) {
			return prt != NULL && prt->E >= 0. && prt->p >= 0. && prt->theta >= 0. && prt->phi >= 0.;
		}
		
		/** 
		 Determines of ParticleS has had p and theta smeared to sensible values.
		 */
		bool isBeamLeptonAccepted(ParticleS *prt) {
			return prt != NULL && prt->p >= 0. && prt->theta >= 0.;
		}
		
		/**
		 Store the default error value in all (x,y,QSquared,WSquared) output slots.
		 */
		void ErrorValues() {
			OutKin.x = fault;
			OutKin.y = fault;
			OutKin.QSquared = fault;
			OutKin.WSquared = fault;
		}
		
		/**
		 Compute WSqaured using the stored values for x, y and QSquared.  This is automatically called whenver you use
		 ComputeUsingNM, ComputeUsingJB or ComputeUsingDA, so you shouldn't ever need to use it.
		 */
		void ComputeWSquaredFromStoredValues() {
			double &MN = InKin.MNucleon;
			OutKin.WSquared = MN*MN+OutKin.QSquared*(1./OutKin.x-1.);
		}
		
		void SetMissingEnergyTolerance(double E) {
			HadTol = E;
			TolerateBadEvents(false);
		}
		
		void TolerateBadEvents(double b) {
			bTolerant = b;
		}
    
		/**
		 Reads an event into event the EventKinematicsComputer.  event and eventS must both represent
		 events with the same initial conditions for this to make sense.  Currently this does not read
		 in or otherwise use pz or pt values.
       TODO Allow switching between electron momentum and energy for kinematic calculations.
		 */
		void ReadEvent(const EventBase *event, EventS *eventS) {
			
			JBComp.clearParticles();
			DAComp.clearParticles();
			
			bGoodJohnnieLept = false;
			bGoodJohnnieHad  = false;
			bTolerated = true;
			
			int nParticles = event->GetNTracks();
			
			bool foundN = false;
			bool foundL = false;
			bool foundS = false;
			
			bool accept; bool acceptL;
			
			TLorentzVector total;
			
			for (int j=0; j<nParticles; j++) {
				
				accept = isAccepted(eventS->GetTrack(j));
				acceptL = isBeamLeptonAccepted(eventS->GetTrack(j));
				
           const Particle* part = event->GetTrack(j);
			  if (not foundN and PartIdent.isBeamNucleon(*part)) {
				  InKin.ENucleon = event->GetTrack(j)->GetE();
              InKin.MNucleon = event->GetTrack(j)->GetM();
//				  InKin.ENucleon = event->GetTrack(j)->E;
//					InKin.MNucleon = event->GetTrack(j)->m;
					foundN = true;
					JBComp.setBeamHadron(event->GetTrack(j)->Get4Vector());
					DAComp.setBeamHadron(event->GetTrack(j)->Get4Vector());
              // Store the incident beams even in the smeared event, to make
              // use of certain utility functions.
              //              eventS->GetTrack(1) = new Particle(event->particles.at(1));
//              std::cout << "Found beam nucleon"<<std::endl;
//              eventS->GetTrack(1)->Print();
					//              event->GetTrack(j)->Dump();
				}
				else if (not foundL and PartIdent.isBeamLepton(*part)) {
					InKin.Ee = event->GetTrack(j)->GetE();
					foundL = true;
					JBComp.setBeamLepton(event->GetTrack(j)->Get4Vector());
					DAComp.setBeamLepton(event->GetTrack(j)->Get4Vector());
//               std::cout << "Found beam lepton"<<std::endl;
					//               event->GetTrack(j)->Dump();
               //               eventS->GetTrack(0) = new Particle(event->particles.at(0));
               //               eventS->GetTrack(0)->Print();
				}
				else if (!foundS && PartIdent.isScatteredLepton(*part) && acceptL) {
					InKin.EPrime = eventS->GetTrack(j)->p;
					InKin.thetae = eventS->GetTrack(j)->theta;
					foundS = true;
					//               std::cout << "Found scattered lepton"<<std::endl;
					//               event->GetTrack(j)->Dump();
				}
				/*
				 else if (!foundS && PartIdent.isScatteredLepton(event->particles.at(j))) {
				 std::cout << accept << " p=" << eventS->GetTrack(j)->p << " theta=" << eventS->GetTrack(j)->theta << std::endl;
				 }*/
				else if (needBackground && accept) {
					
					TLorentzVector v = eventS->GetTrack(j)->Get4Vector();
//					TLorentzVector lv = event->GetTrack(j)->Get4Vector();
					//               std::cout << "4-vector:"<<std::endl;
					//               std::cout << v.E() << " " << v.GetPx() << " " << v.GetPy() << " " << v.GetPz() << std::endl;
					//               std::cout << tmpE << " " << tmpP << " " << tmpTheta << " " << tmpPhi << std::endl;
					//               std::cout << lv.E() << " " << lv.GetPx() << " " << lv.GetPy() << " " << lv.GetPz() << std::endl;
//               std::cout <<"Adding " << v.E() << " " << v.GetPx() << " " << v.GetPy() << " " << v.GetPz() << std::endl;
					JBComp.addParticle(v);
					DAComp.addParticle(v);
					//               std::cout <<"Found final-state hadron" << std::endl;
					//               event->GetTrack(j)->Dump();
					
					total += v;
					
				}
			}//for particles
			
			if(!bMute && !bTolerant && needBackground && abs(event->HadronicFinalStateMomentum().E() - total.E())>HadTol){
				bTolerated = false;
				std::cout<<"\nERROR in hadronic kinematics!!!\n"<<std::endl;
				std::cout<<"Warning! Smeared hadronic final state energy exceeds that of MC data by "
				<< event->HadronicFinalStateMomentum().E() - total.E() << " while tolerance is " << HadTol << ".\n";
				//std::cout << "From smear:"<<std::endl;
				//total.Dump();
				//std::cout << "From Event:"<<std::endl;
				//event->HadronicFinalStateMomentum().Dump();
			} 
			
			//			JBComp.setLeptonEnergy(InKin.Ee);
			//			DAComp.setLeptonEnergy(InKin.Ee);
			//			JBComp.setMandelstamS(4.*InKin.Ee*InKin.ENucleon);
			//			DAComp.setMandelstamS(4.*InKin.Ee*InKin.ENucleon);
			//         std::cout << "Ee = " << InKin.Ee << " s = " << 4.*InKin.Ee*InKin.ENucleon << std::endl;
			DAComp.setLeptonAngle(InKin.thetae);
			//			DAComp.setProtonEnergy(InKin.ENucleon);
			
			//         myJb.setLeptonEnergy(InKin.Ee);
			//         myJb.setMandelstamS(4.*InKin.Ee*InKin.ENucleon);
			//         std::cout << "My yJB = " << myJb.computeY() << std::endl;
			//         std::cout << "My Q2JB = " << myJb.computeQSquared() << std::endl;
			//         std::cout << "My xJB = " << myJb.computeX() << std::endl;
			
			bGoodJohnnieLept = foundN && foundL && foundS;
			bGoodJohnnieHad  = foundN && foundL;
			
			if (!(foundL && foundN && foundS) && !bMute) {
				std::cerr <<
				"Failed to find all required particles for smeared event kinematics:" <<
				"\nfound nucleon:  " << foundN <<
				"\nfound lepton:   " << foundL <<
				"\nfound scattered:" << foundS << std::endl;
			} 
			
		}
		
		/**
		 Compute event-wise kinematic variables using the null momentum approximation and the scattered
		 electron.  Computation uses event read in via ReadEvent and calculated values are accessible 
		 via x, y, QSquared and WSquared.  If no event was read in, or the last event read in was bad in some
		 way (id est, has missing initial state particles, or missing scattered electron) this will store default
		 bogus values.
		 */
		void ComputeUsingNM() {
			
			if (bGoodJohnnieLept) {
				
				double &E = InKin.Ee;
				double &EPrime = InKin.EPrime;
				double &theta = InKin.thetae;
				double &ENucleon = InKin.ENucleon;
//				double &MN = InKin.MNucleon;
				
				OutKin.QSquared = 2.*E*EPrime*(1.+cos(theta));
				OutKin.y = 1.-(EPrime/(2.*E))*(1.-cos(theta));
				OutKin.x = (EPrime/(2.*OutKin.y*ENucleon))*(1.+cos(theta));
				
				ComputeWSquaredFromStoredValues();
				
			} else ErrorValues();
			
		}
		
		/**
		 Compute event-wise kinematic variables using the Jacquet-Blondel method.  Computation uses event
		 read in via ReadEvent and calculated values are accessible via x, y, QSquared and WSquared.
		 If the event read in was bad, this will return default bogus values, but it doesn't care how many final
		 state particles it picked up.
		 */
		void ComputeUsingJB(EventBase* = NULL) {
//         std::cout << "Computing using JB" << std::endl;
			if (bGoodJohnnieHad && (bTolerated || bTolerant)) {
            /*
            if(e) {
               OutKin.y = mJBComp.ComputeY(*e);
               OutKin.QSquared = mJBComp.ComputeQSquared(*e);
               OutKin.x = mJBComp.ComputeX(*e);
               OutKin.WSquared = mJBComp.ComputeWSquared(*e);
            } // if
            else {
             */
               OutKin.y = JBComp.computeY();
               OutKin.QSquared = JBComp.computeQSquared();
               OutKin.x = JBComp.computeX();
               ComputeWSquaredFromStoredValues();
               /*
               OutKin.y = 0.;
               OutKin.QSquared = 1.;
               OutKin.x = 2.;
               OutKin.WSquared = 3.;
           }
                */
			} else ErrorValues();
		}
		
		/**
		 Compute event-wise kinematic variables using the double angle method.  Coputation uses event
		 read in via ReadEvent and calculated values are accessible via x, y, QSquared and WSquared.
		 If the event read in was bad, this will return default bogus values, but it doesn't care how many final
		 state particles it picked up.
		 */
		void ComputeUsingDA() {
			if (bGoodJohnnieHad && (bTolerated || bTolerant)) {
				OutKin.QSquared = DAComp.computeQSquared();
				OutKin.x = DAComp.computeX();
				OutKin.y = DAComp.computeY();
				ComputeWSquaredFromStoredValues();
			} else ErrorValues();
		}
		
		ClassDef(EventKinematicsComputer,1)
	};
	
	
}

#endif
