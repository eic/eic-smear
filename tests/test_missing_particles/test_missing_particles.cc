#include <iomanip>

#include <TDatabasePDG.h>

#include "BeASTDetector.h"
#include "ePHENIXDetector.h"
#include "ZeusDetector.h"

const double deg_to_rad = 0.01745329251; // pi/180

using std::setw;

enum class EicSmearResults {
    null_particle,
    zero_e_smear_p,
    smear_e_zero_p,
    smear_e_smear_p,
    zero_e_zero_p,
};

struct EicSmearStep {
    TLorentzVector not_smeared_lorentz;
    EicSmearResults result;
};

// Some statistics about the smearing process
struct EicSmearStatistics {
    uint64_t total_particles = 0;
    uint64_t null_particles = 0;
    uint64_t zero_e_smear_p = 0;
    uint64_t smear_e_zero_p = 0;
    uint64_t smear_e_smear_p = 0;
    uint64_t zero_e_zero_p = 0;

    std::vector<EicSmearStep> steps;
};

/// This function does smearing itself. by calling detector.Smear(not_smeared_prt);
EicSmearStep DoSmearStep(int pdg, TLorentzVector& input_vect, Smear::Detector& detector) {
    static int particle_index = 0;
    EicSmearStep step;

    step.not_smeared_lorentz = input_vect;

    // Create not smeared particle
    erhic::ParticleMC not_smeared_prt;
    not_smeared_prt.SetId(pdg);                   // PDG particle code
    not_smeared_prt.SetIndex(particle_index++);   // Particle index in event
    not_smeared_prt.SetStatus(1);               // Particle status code: like in PYTHIA, Beagle, etc
    not_smeared_prt.Set4Vector(input_vect);

    // Smear the particle
    Smear::ParticleMCS * smeared_prt;
    smeared_prt = detector.Smear(not_smeared_prt);

    if(!smeared_prt) {
        step.result = EicSmearResults::null_particle;
        return step;
    }

    // Check odds of eic-smear smearing
    double in_p = input_vect.P();
    double sm_p = smeared_prt->GetP();
    double in_e = input_vect.E();
    double sm_e = smeared_prt->GetE();

    // Eic smear return non zero p
    bool zero_p = TMath::Abs(in_p)>0.001 && TMath::Abs(sm_p)<0.00001;
    bool zero_e = TMath::Abs(in_e)>0.001 && TMath::Abs(sm_e)<0.00001;

    if(zero_p && zero_e) {
        step.result = EicSmearResults::zero_e_zero_p;   // we don't take such particle
    }
    else if (zero_p) {
        step.result = EicSmearResults::smear_e_zero_p;  // we don't take such particle
    }
    else if (zero_e) {
        step.result = EicSmearResults::zero_e_smear_p;  // we don't take such particle
    }
    else {
        step.result = EicSmearResults::smear_e_smear_p;
        // that! is the particle we take...

        // In EJana we do:
        // Copy back the particle
        // dest_particle->px = smeared_prt->GetPx();
        // dest_particle->py = smeared_prt->GetPy();
        // dest_particle->pz = smeared_prt->GetPz();
        // dest_particle->tot_e = smeared_prt->GetE();

    }
    return step;
}


// This function sends
EicSmearStatistics Process(int pdg, Smear::Detector& detector) {

    using namespace std;

    // Get the inputs needed for this factory.
    auto db = TDatabasePDG::Instance();

    auto pdg_particle = db->GetParticle(pdg);
    EicSmearStatistics stat;
    
    for(int mom=1; mom < 20; mom+=2) {
        for(int angle_deg=0; angle_deg < 360; angle_deg+=1) {  // theta
            // 4 vector
	    TLorentzVector input_vect(0, 0, mom, sqrt ( mom*mom + pow(pdg_particle->Mass(),2)) );
	    input_vect.RotateX(angle_deg*deg_to_rad);
	    // input_vect.Print();
  
            auto step = DoSmearStep(pdg, input_vect, detector);

            // Update statistics
            switch(step.result) {
                case EicSmearResults::null_particle:
                    stat.null_particles++;
                    break;
                case EicSmearResults::zero_e_smear_p:
                    stat.zero_e_smear_p++;
                    break;
                case EicSmearResults::smear_e_zero_p:
                    stat.smear_e_zero_p++;
                    break;
                case EicSmearResults::zero_e_zero_p:
                    stat.zero_e_zero_p++;
                    break;
                case EicSmearResults::smear_e_smear_p:
                    stat.smear_e_smear_p++;
                    break;
            }

            stat.steps.push_back(step);
        }
    }
    return stat;
}


void PrintSmearStats(const EicSmearStatistics& stat) {
    using namespace std;
    cout<<"SmearingFactory statistics:\n";
    cout << "   total_particles = " << setw(8) << stat.total_particles << " \n";
    cout << "   null_particles  = " << setw(8) << stat.null_particles << " (not smeared)\n";
    cout << "   zero_e_zero_p   = " << setw(8) << stat.zero_e_zero_p << " (not smeared)\n";
    cout << "   zero_e_smear_p  = " << setw(8) << stat.zero_e_smear_p << " (partly smeared)\n";
    cout << "   smear_e_zero_p  = " << setw(8) << stat.smear_e_zero_p << " (partly smeared)\n";
    cout << "   smear_e_smear_p = " << setw(8) << stat.smear_e_smear_p << " (smeared)\n";
}



int main() {

    Smear::Detector beast_detector = BuildBeAST();
    Smear::Detector zeus_detector = BuildZeus();
    Smear::Detector ephoenix_detector = BuildEphoenix();
    
    // auto stat = Process(11, ephoenix_detector); // 11 - electron
    // auto stat = Process(22, ephoenix_detector); // 22 - gamma
    // auto stat = Process(2212, ephoenix_detector); // 2212 - proton
    auto stat = Process(211, beast_detector); // 211 - pi+
    PrintSmearStats(stat);
    return 0;
}
