#include <iomanip>

#include <TDatabasePDG.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>

#include "ePHENIXDetector.h"

#include "eicsmear/smear/Detector.h"
Smear::Detector BuildZeus();
Smear::Detector BuildBeAST();

const double deg_to_rad = 0.01745329251; // pi/180

using std::setw;

enum EicSmearResults {
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
    long int total_particles;
    long int null_particles;
    long int zero_e_smear_p;
    long int smear_e_zero_p;
    long int smear_e_smear_p;
    long int zero_e_zero_p;

    std::vector<EicSmearStep> steps;

    TH1D* null_particles_eta;
    TH1D* zero_e_smear_p_eta;
    TH1D* smear_e_zero_p_eta;
    TH1D* zero_e_zero_p_eta;
    TH1D* smear_e_smear_p_eta;
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
        step.result = null_particle;
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
        step.result = zero_e_zero_p;   // we don't take such particle
    }
    else if (zero_p) {
        step.result = smear_e_zero_p;  // we don't take such particle
    }
    else if (zero_e) {
        step.result = zero_e_smear_p;  // we don't take such particle
    }
    else {
        step.result = smear_e_smear_p;
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
    TDatabasePDG* db = TDatabasePDG::Instance();

    TParticlePDG* pdg_particle = db->GetParticle(pdg);
    EicSmearStatistics stat;
    stat.total_particles=0;
    stat.null_particles=0;
    stat.zero_e_smear_p=0;
    stat.smear_e_zero_p=0;
    stat.smear_e_smear_p=0;
    stat.zero_e_zero_p=0;

    stat.null_particles_eta  = new TH1D("null_particles_eta",  "Unsmeared particles;#eta;counts", 200, -5, 5 );
    stat.zero_e_smear_p_eta  = new TH1D("zero_e_smear_p_eta",  "only p smeared;#eta;counts", 100, -5, 5 );
    stat.smear_e_zero_p_eta  = new TH1D("smear_e_zero_p_eta",  "only e smeared;#eta;counts", 100, -5, 5 );
    stat.zero_e_zero_p_eta   = new TH1D("zero_e_zero_p_eta",   "Unsmeared particles;#eta;counts", 100, -5, 5 );
    stat.smear_e_smear_p_eta = new TH1D("smear_e_smear_p_eta", "smeared e, smeared p;#eta;counts", 100, -5, 5 );

    int angleinc = 1;
    for(int mom=1; mom < 20; mom+=2) {
        for(int angle_deg=angleinc; angle_deg < 180; angle_deg+=angleinc) {  // theta
            // 4 vector
	    TLorentzVector input_vect(0, 0, mom, sqrt ( mom*mom + pow(pdg_particle->Mass(),2)) );
	    input_vect.RotateX(angle_deg*deg_to_rad);
	    // input_vect.Print();
  
            EicSmearStep step = DoSmearStep(pdg, input_vect, detector);

            // Update statistics
            switch(step.result) {
                case null_particle:
                    stat.null_particles++;
		    stat.null_particles_eta->Fill(input_vect.Eta());
                    break;
                case zero_e_smear_p:
                    stat.zero_e_smear_p++;
                    stat.zero_e_smear_p_eta->Fill(input_vect.Eta());
                    break;
                case smear_e_zero_p:
                    stat.smear_e_zero_p++;
                    stat.smear_e_zero_p_eta->Fill(input_vect.Eta());
                    break;
                case zero_e_zero_p:
                    stat.zero_e_zero_p++;
                    stat.zero_e_zero_p_eta->Fill(input_vect.Eta());
                    break;
                case smear_e_smear_p:
                    stat.smear_e_smear_p++;
                    stat.smear_e_smear_p_eta->Fill(input_vect.Eta());
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
    int pid = 211; // pi+
    
    TString detstring = "BeAST";
    // TString detstring = "ePhenix";
    // TString detstring = "ZEUS";
    
    Smear::Detector detector;
    if ( detstring=="BeAST" ) detector = BuildBeAST();
    if ( detstring=="ZEUS" ) detector = BuildZeus();
    if ( detstring=="ePhenix" ) detector = BuildEphoenix();
    
    EicSmearStatistics stat = Process( pid, detector);
    PrintSmearStats(stat);

    gStyle->SetOptStat(0);
    gStyle->SetHistLineWidth(2);
    float lMargin = 0.12;
    float bMargin = 0.12;
    gStyle->SetLabelSize(.05, "XY");
    gStyle->SetTitleSize(.05, "XY");
    gStyle->SetTitleOffset(1.1,"x");	//X-axis title offset from axis
    gStyle->SetTitleOffset(1.1,"y");	//Y-axis title offset from axis
 
    new TCanvas;
    TString title = detstring;
    title += ", pid="; title +=pid;
    TLegend * leg = new TLegend( 0.65, 0.65, 0.95, 0.95, title);
    stat.smear_e_smear_p_eta->Draw("AXIS");
	
    stat.smear_e_smear_p_eta->SetLineColor(kGreen+1);
    stat.smear_e_smear_p_eta->Draw("same");

    stat.zero_e_smear_p_eta->SetLineColor(kBlue);
    stat.zero_e_smear_p_eta->Draw("same");
    
    stat.smear_e_zero_p_eta->SetLineColor(kBlack);
    stat.smear_e_zero_p_eta->Draw("same");
        
    stat.null_particles_eta->SetLineColor(kRed);
    stat.null_particles_eta->Draw("same");
    
    leg->AddEntry( stat.zero_e_smear_p_eta, "P smeared", "l");
    leg->AddEntry( stat.smear_e_zero_p_eta, "E smeared", "l");
    leg->AddEntry( stat.smear_e_smear_p_eta, "Both smeared", "l");
    leg->AddEntry( stat.null_particles_eta, "NONE smeared", "l");
    leg->Draw("same");

    TString outname=detstring;
    outname += "_"; outname+=pid;
    outname += "_eta.png";
    gPad->SaveAs(outname);
    gPad->SaveAs("etaplot.png"); // for quick looks at the most recent one
    
    
    return 0;
}
