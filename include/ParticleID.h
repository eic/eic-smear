#ifndef _ERHIC_BUILDTREE_PARTICLEID_
#define _ERHIC_BUILDTREE_PARTICLEID_

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

#include "VirtualParticle.h"
#include "Kinematics.h"
#include "ParticleIdentifier.h"

#include "Smear.h"
#include "Device.h"
#include "Particle.h"

// TODO make this inherit from Device

namespace Smear {
	
	/**
	 This structure is used to generate particle ID.  It remains somewhat iflexible, but is much
	 better than before.  
	 
	 The input file containing the P matrix must begin with format lines beginning with "!T", "!F"
	 and "!P".  For example
	 
	 !T 211 321 2212
	 !F 211 321 2212 0
	 !P 15
	 
	 The first line tells the PID generator which particles false ID will be generated for.  The second line tells
	 it which particles the particles in the first line can be identified as.  The third line tells it how many
	 momentum bins there are.
	 
	 Lines of data should appear as
	 
	 1 pbinNumber  pmin pmax  FalseID  P1  P2  P3
	 
	 The line must begin with a 1 to be read.  pbinNumber is the number of the momentum bin.  pmin and pmax
	 are the bounds of the momentum bin.  FalseID is the ID that the particle will be mis-identified as,
	 this should be a PDG particle code (unless it is an old Hermes file).  P1 is the probability that the
	 first particle appearing in the "!T" line will be misidentified as FalseID.  Likewise for P2 and P3.
	 
	 Detectors come with their own ParticleID instance.
	 */
	struct ParticleID {
		
		/**
       \todo This needs to be settable by the user - provide static method
       that can be called in a ROOT logon macro
		 Default constructor.  Default path for probability matrix is "pmat_EVT_all.asc".  Also, by
		 default, the particle ID uses smeared kinematic variables as input, rather than Monte Carlo 
		 values.
		 */
		ParticleID() {
			PMatPath = "/afs/rhic.bnl.gov/eic/MACROS/BuildTree/PIDMatrix.dat";
			
			ReadP(PMatPath);
			
			Ran.SetSeed(0);
			
			bUseMC = false;
		}
		
		ParticleID(TString matrix) {
			PMatPath = matrix;
			
			ReadP(PMatPath);
			
			Ran.SetSeed(0);
			
			bUseMC = false;
		}
	
		TRandom3 Ran;
		
		TString PMatPath;
		
		std::vector<int> TrueIdent;
		std::vector<int> FalseIdent;
		std::vector<double> PMin;
		std::vector<double> PMax;
		//indices are momentum bin, true, false
		std::vector< std::vector<std::vector<double> > > Range;
		std::vector< std::vector<std::vector<double> > > PMatrix;
		
		Acceptance Accept;
	
		bool bUseMC;
    
		/** 
		 Set the path to the file containing the particle ID probability matrix, and read it in.
		 */
		void SetPMatrixPath(TString str) {
		  PMatPath = str;
			ReadP(PMatPath);
		}
	
		/** 
		 If true, the ParticleID will use the original Monte Carlo values to generate ID's, rather 
		 than smeared values.  Default is false.
		 */
		void SetPIDUseMC(bool d) {
			bUseMC = d;
		}
		
		/**
		 Return the Root TRandom3 random number generator instance used by the ParticleID.
		 */
		TRandom3 GetRandomGenerator() {
			return Ran;
		}
		
		/**
		 Set the seed of the Root TRandom3 random number generator instance used by the 
		 ParticleID.  See Root TRandom3 documentation for more information.
		 */
		void SetRanSeed(int n) {
			Ran.SetSeed(n);
		}
		
		/**
		 Set the ParticleID instance to have the same acceptance as the device dev.  This includes all acceptance zones as
		 well as the list of specific particles to be detected (if there is one).
		 */
		void GetAcceptanceFromDevice(Device dev) {
			Accept = dev.Accept;
		}
		
		void operator&= (Device dev) {
			GetAcceptanceFromDevice(dev);
		}
		
		//resize the PMatrix, (# of pbins, # of true, # of false)
		void SetPMatrixSize() {
			PMatrix.resize(PMin.size());
			for (unsigned i=0; i<PMatrix.size(); i++) {
				PMatrix.at(i).resize(TrueIdent.size());
				for (unsigned j=0; j<PMatrix.at(i).size(); j++) {
					PMatrix.at(i).at(j).resize(FalseIdent.size());
				}
			}
		}
	
		//setup the range array used by wild. range is used to set up zones in [0,1],
		//if a random number falls in one of the zones, the associated ID is used
		void SetupProbabilityArray() {
			double t=0.;
			Range.resize(PMatrix.size());
			for(unsigned i=0; i<Range.size(); i++) {
				Range.at(i).resize(PMatrix.at(i).size());
				for(unsigned j=0; j<Range.at(i).size(); j++) {
					Range.at(i).at(j).resize(PMatrix.at(i).at(j).size());
					for(unsigned k=0; k<Range.at(i).at(j).size(); k++) {
						t += PMatrix.at(i).at(j).at(k);
						Range.at(i).at(j).at(k) = t;
					}
					t = 0.;
				}
			}
		}
	
		virtual ParticleID* Clone() {
			ParticleID *dev = new ParticleID();
			*dev = *this;
			return dev;
		}
		
		//returns the false id, input P matrix, k (mom bin), j (id of true particle)
		int Wild(int pbin, int trueID) {
			
			double r = Ran.Rndm();
			
			int k = inListofTrue(trueID);
			
			for(unsigned i=0; i<Range.at(pbin).at(k).size(); i++) {
				if (i==0) {
					if (r < Range.at(pbin).at(k).at(i)) {
						return FalseIdent.at(i);
					}
				}
				else if (r > Range.at(pbin).at(k).at(i-1) && r < Range.at(pbin).at(k).at(i)) {
					return FalseIdent.at(i); 
				}
			}
         
         return -999999999;
		}
		
		int inListofTrue(int ID) {
			for(unsigned i=0; i<TrueIdent.size(); i++) {
				if (TrueIdent.at(i)==abs(ID)) return i;
			}
			return -1;
		}
		
		int inListofFalse(int ID) {
			for(unsigned i=0; i<FalseIdent.size(); i++) {
				if (FalseIdent.at(i)==abs(ID)) return i;
			}
			//this is to accomodate the old HERMES format
			if (ID < 5 && ID > 0) {
				return ID-1;
			}
			return -1;
		}

		
		/**
		 Read in a P matrix and set up the ParticleID instance to be ready to generate.  See the documentation for the detector
		 function ReadPIDMatrix.
		 */
		void ReadP(TString filename) {
			
			TrueIdent.clear();
			FalseIdent.clear();
			PMin.clear(); PMax.clear();
			
			char dumb[2];
			int tmpInt1; int tmpInt2; int tmpInt3; int tmpInt4;
			
			float tmpFlt1; float tmpFlt2; float tmpFlt3;// float tmpFlt4;
			float pmin; float pmax;
			
			bool firstTime=true;
			bool gotTrue=false; bool gotFalse=false; bool gotBins=false;
			
			std::ifstream Qfile;
			char oneline[256];
			
			Qfile.open(filename);
			if(!Qfile.is_open()){
				std::cout << "PID ERROR! unable to open file " << filename << std::endl;
				return;
			}
			
			//Loop over each line in the input file
			while (!Qfile.eof() ) {   
				Qfile.getline(oneline,256);
				if (oneline[0]=='!' && oneline[1]=='T') {
					sscanf(oneline,"%s %d %d %d",dumb,&tmpInt1,&tmpInt2,&tmpInt3);
					gotTrue = true;
					TrueIdent.push_back(tmpInt1);
					TrueIdent.push_back(tmpInt2);
					TrueIdent.push_back(tmpInt3);
				}
				if (oneline[0]=='!' && oneline[1]=='F') {
					sscanf(oneline,"%s %d %d %d %d",dumb,&tmpInt1,&tmpInt2,&tmpInt3,&tmpInt4);
					gotFalse = true;
					FalseIdent.push_back(tmpInt1);
					FalseIdent.push_back(tmpInt2);
					FalseIdent.push_back(tmpInt3);
					FalseIdent.push_back(tmpInt4);
				}
				if (oneline[0]=='!' && oneline[1]=='P') {
					sscanf(oneline,"%s %d",dumb,&tmpInt1);
					gotBins = true;
					PMin.resize(tmpInt1);
					PMax.resize(tmpInt1);
				}
				if (oneline[0]=='1' || oneline[0]=='2' || oneline[0]=='3') {
					
					if (!(gotTrue && gotFalse && gotBins)) {
						std::cerr << "PID ERROR! P matrix input file has bad or missing format lines.\n";
						return;
					}
					
					sscanf(oneline,"%d %d %g %g %d %g %g %g",&tmpInt1,&tmpInt2,&pmin,&pmax,&tmpInt3,&tmpFlt1,&tmpFlt2,&tmpFlt3);
					
					if ((unsigned)tmpInt2 > PMin.size()) {
						std::cerr << "PID ERROR! Out of bounds momentum bin listing.\n";
						return;
					}

					tmpInt2 = tmpInt2-1;
					
					PMin.at(tmpInt2) = pmin;
					PMax.at(tmpInt2) = pmax;
					
					tmpInt3 = inListofFalse(tmpInt3);
					
					if ((unsigned)tmpInt3 > FalseIdent.size()) {
						std::cerr << "PID ERROR! P matrix has bad particle listing.\n";
						return;
					}
					
					if (firstTime) {
						SetPMatrixSize();
						firstTime = false;
					}
						
					if (tmpInt1==1) {
						PMatrix.at(tmpInt2).at(0).at(tmpInt3) = tmpFlt1;
						PMatrix.at(tmpInt2).at(1).at(tmpInt3) = tmpFlt2;
						PMatrix.at(tmpInt2).at(2).at(tmpInt3) = tmpFlt3;
					}
				} 
			}
			
			Qfile.close();
			SetupProbabilityArray();
		}
		
		/**
		 Generates particle ID if the particle is in the list of particles to be identified, falls within
		 the momentum range of the PMatrix and falls within acceptence. New id is stored in prtOut.id.  By default, this
		 will use the momentum prtOut.p to make its determination, but you can set it to use
		 the values stored in prt instead using SetPIDUseMC(true).
		 */
		void DevSmear(Particle prt, ParticleS &prtOut) {
			double up;
			if (bUseMC) {
				up = prt.GetP();
			} else {
				up = prtOut.p;
			}
			if (inListofTrue(prt.Id()) not_eq -1 && Accept.Is(prt)) {
				for(unsigned i=0; i<PMin.size(); i++) {
					if (up > PMin.at(i) && up < PMax.at(i)) {
						prtOut.id = sgn(prt.Id())*Wild(i,prt.Id());
					}//if 
				}//for 
			}//if 
		}
			
		ClassDef(ParticleID, 1 )
	}; //end of ParticleID structure
	

	
}

#endif