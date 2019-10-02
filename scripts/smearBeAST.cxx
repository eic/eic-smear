Smear::Detector BuildBeAST() {

  gSystem->Load("libeicsmear");

   // Calorimeter resolution usually given as sigma_E/E = const% + stocastic%/Sqrt{E}
   // EIC Smear needs absolute sigma: sigma_E = Sqrt{const*const*E*E + stoc*stoc*E}

   // Create the EM Calorimeter
   Smear::Device emcalBck(Smear::kE, "sqrt(0.01*0.01*E*E + 0.015*0.015*E)");
   Smear::Device emcalMidBck(Smear::kE, "sqrt(0.01*0.01*E*E + 0.07*0.07*E)");
   Smear::Device emcalMid(Smear::kE, "sqrt(0.01*0.01*E*E + 0.10*0.10*E)");
   Smear::Device emcalFwd(Smear::kE, "sqrt(0.01*0.01*E*E + 0.07*0.07*E)");


   // Create the Forward/Backward Hadron Calorimeter
   Smear::Device hcalFwd(Smear::kE, "sqrt(0.015*0.015*E*E + 0.50*0.50*E)");
   Smear::Device hcalBck(Smear::kE, "sqrt(0.015*0.015*E*E + 0.50*0.50*E)");

   // Create the Hypothetical Mid Rap Calorimeter
   //Smear::Device hcalMid(Smear::kE, "sqrt(0.07*0.07*E*E + 0.85*0.85*E)"); // ~CMS
   Smear::Device hcalMid(Smear::kE, "sqrt(0.02*0.02*E*E + 0.35*0.35*E)"); // ~Zeus

   // Create Forward/Backward Hadron Calorimeter to Measure Charged Hadrons after the Tracker
   Smear::Device hcalTrkFwd(Smear::kE, "sqrt(0.015*0.015*E*E + 0.50*0.50*E)");
   Smear::Device hcalTrkBck(Smear::kE, "sqrt(0.015*0.015*E*E + 0.50*0.50*E)");


   // Create our tracking capabilities, by a combination of mometum, theta and phi Devices.
   // The momentum parametrization (a*p + b) gives sigma_P/P in percent. 
   // So Multiply through by P and divide by 100 to get absolute sigma_P
   // Theta and Phi parametrizations give absolute sigma in miliradians

   // Track Momentum
   //Smear::Device momentum(Smear::kP, "(P*P*(0.0182031 + 0.00921047*pow((-log(tan(theta/2.0))), 2) - 0.00291243*pow((-log(tan(theta/2.0))), 4) + 0.000264353*pow((-log(tan(theta/2.0))), 6)) + (0.209681 + 0.275144*pow((-log(tan(theta/2.0))), 2) - 0.0436536*pow((-log(tan(theta/2.0))), 4) + 0.00367412*pow((-log(tan(theta/2.0))), 6)))*0.01");
   Smear::Device momentum(Smear::kP, "(P*P*(0.0182031 + 0.00921047*pow((-log(tan(theta/2.0))), 2) - 0.00291243*pow((-log(tan(theta/2.0))), 4) + 0.000264353*pow((-log(tan(theta/2.0))), 6)) + P*(0.209681 + 0.275144*pow((-log(tan(theta/2.0))), 2) - 0.0436536*pow((-log(tan(theta/2.0))), 4) + 0.00367412*pow((-log(tan(theta/2.0))), 6)))*0.01");
   Smear::Device trackTheta(Smear::kTheta, "((1.0/(1.0*P))*(0.752935 + 0.280370*pow((-log(tan(theta/2.0))), 2) - 0.0359713*pow((-log(tan(theta/2.0))), 4) + 0.00200623*pow((-log(tan(theta/2.0))), 6)) + 0.0282315 - 0.00998623*pow((-log(tan(theta/2.0))), 2) + 0.00117487*pow((-log(tan(theta/2.0))), 4) - 0.0000443918*pow((-log(tan(theta/2.0))), 6))*0.001");
   Smear::Device trackPhi(Smear::kPhi, "((1.0/(1.0*P))*(0.743977 + 0.753393*pow((-log(tan(theta/2.0))), 2) + 0.0634184*pow((-log(tan(theta/2.0))), 4) + 0.0128001*pow((-log(tan(theta/2.0))), 6)) + 0.0308753 + 0.0480770*pow((-log(tan(theta/2.0))), 2) - 0.0129859*pow((-log(tan(theta/2.0))), 4) + 0.00109374*pow((-log(tan(theta/2.0))), 6))*0.001");


   // Need these to keep the component not smeared 

   // Momentum for EM
   Smear::Device momentumEM(Smear::kP, "0");
   Smear::Device trackThetaEM(Smear::kTheta, "0");
   Smear::Device trackPhiEM(Smear::kPhi, "0");

   // Momentum for Neutral Hadrons
   Smear::Device momentumHad(Smear::kP, "0");
   Smear::Device trackThetaHad(Smear::kTheta, "0");
   Smear::Device trackPhiHad(Smear::kPhi, "0");

   // Momentum for Charged Hadrons After the Tracker
   Smear::Device momentumHadTrkBck(Smear::kP, "0");
   Smear::Device trackThetaHadTrkBck(Smear::kTheta, "0");
   Smear::Device trackPhiHadTrkBck(Smear::kPhi, "0");

   Smear::Device momentumHadTrkFwd(Smear::kP, "0");
   Smear::Device trackThetaHadTrkFwd(Smear::kTheta, "0");
   Smear::Device trackPhiHadTrkFwd(Smear::kPhi, "0");

   // Energy for Tracks
   Smear::Device trackEnergy(Smear::kE, "0");

   // Create a smearer for momentum
   //Smear::Device momentum(Smear::kP, "0.01 * P");
   //Smear::Device theta(Smear::kTheta, "0.05 * P");
   //Smear::Device phi(Smear::kPhi, "0"); // "0" indicates perf585ect performance i.e. sigma(phi) = 0

   // Create a based on Hermes RICH.
   //Smear::ParticleID rich("PIDMatrix.dat");


   // Set Spatial Extent of Devices (in radians)
   // eta = -4.5 -> theta = 3.1194
   // eta = -4.0 -> theta = 3.1050
   // eta = -3.5 -> theta = 3.0812
   // eta = -3.0 -> theta = 3.0421
   // eta = -2.0 -> theta = 2.8726
   // eta = -1.0 -> theta = 2.4366
   // eta = +1.0 -> theta = 0.7050
   // eta = +2.0 -> theta = 0.2690
   // eta = +3.0 -> theta = 0.0995
   // eta = +3.5 -> theta = 0.0604
   // eta = +4.0 -> theta = 0.0366
   // eta = +4.5 -> theta = 0.0222

   // The BeAST design has calorimetry extending to +-4, here extend to +-4.5 to match the extent of my particle level jets
   // Can place a cut on jet thrust axis in the analysis if need be

   // Set Up EMCal Zones
   Smear::Acceptance::Zone emBck(2.8726,3.1194,0.,TMath::TwoPi(),0.,TMath::Infinity(),0.,TMath::Infinity(),0.,TMath::Infinity(),-TMath::Infinity(),TMath::Infinity());
   Smear::Acceptance::Zone emMidBck(2.4366,2.8726,0.,TMath::TwoPi(),0.,TMath::Infinity(),0.,TMath::Infinity(),0.,TMath::Infinity(),-TMath::Infinity(),TMath::Infinity());
   Smear::Acceptance::Zone emMid(0.7050,2.4366,0.,TMath::TwoPi(),0.,TMath::Infinity(),0.,TMath::Infinity(),0.,TMath::Infinity(),-TMath::Infinity(),TMath::Infinity());
   Smear::Acceptance::Zone emFwd(0.0222,0.7050,0.,TMath::TwoPi(),0.,TMath::Infinity(),0.,TMath::Infinity(),0.,TMath::Infinity(),-TMath::Infinity(),TMath::Infinity());
   Smear::Acceptance::Zone emTot(0.0222,3.1194,0.,TMath::TwoPi(),0.,TMath::Infinity(),0.,TMath::Infinity(),0.,TMath::Infinity(),-TMath::Infinity(),TMath::Infinity());

   // Set Up HCal Zone
   Smear::Acceptance::Zone hBck(2.4366,3.1194,0.,TMath::TwoPi(),0.,TMath::Infinity(),0.,TMath::Infinity(),0.,TMath::Infinity(),-TMath::Infinity(),TMath::Infinity());
   Smear::Acceptance::Zone hFwd(0.0222,0.7050,0.,TMath::TwoPi(),0.,TMath::Infinity(),0.,TMath::Infinity(),0.,TMath::Infinity(),-TMath::Infinity(),TMath::Infinity());
   Smear::Acceptance::Zone hMid(0.7050,2.4366,0.,TMath::TwoPi(),0.,TMath::Infinity(),0.,TMath::Infinity(),0.,TMath::Infinity(),-TMath::Infinity(),TMath::Infinity());
   Smear::Acceptance::Zone hTot(0.0222,3.1194,0.,TMath::TwoPi(),0.,TMath::Infinity(),0.,TMath::Infinity(),0.,TMath::Infinity(),-TMath::Infinity(),TMath::Infinity());
   // Do HCal Over Full Acceptance (-4,4) and just kill neutrals at midrapidity in the jet finder for standard BeAST

   // Set Up HCal After Tracker Zone
   Smear::Acceptance::Zone hTrkBck(3.0812,3.1194,0.,TMath::TwoPi(),0.,TMath::Infinity(),0.,TMath::Infinity(),0.,TMath::Infinity(),-TMath::Infinity(),TMath::Infinity());
   Smear::Acceptance::Zone hTrkFwd(0.0222,0.0604,0.,TMath::TwoPi(),0.,TMath::Infinity(),0.,TMath::Infinity(),0.,TMath::Infinity(),-TMath::Infinity(),TMath::Infinity());

   // Set Up Tracking Zone
   Smear::Acceptance::Zone trk(0.0604,3.0812,0.,TMath::TwoPi(),0.,TMath::Infinity(),0.,TMath::Infinity(),0.,TMath::Infinity(),-TMath::Infinity(),TMath::Infinity());


   // Assign acceptance to calorimeters
   emcalBck.Accept.SetGenre(Smear::kElectromagnetic);
   emcalMidBck.Accept.SetGenre(Smear::kElectromagnetic);
   emcalMid.Accept.SetGenre(Smear::kElectromagnetic);
   emcalFwd.Accept.SetGenre(Smear::kElectromagnetic);

   emcalBck.Accept.AddZone(emBck);
   emcalMidBck.Accept.AddZone(emMidBck);
   emcalMid.Accept.AddZone(emMid);
   emcalFwd.Accept.AddZone(emFwd);

   hcalBck.Accept.AddParticle(2112);
   hcalBck.Accept.AddParticle(-2112);
   hcalBck.Accept.AddParticle(130);
   hcalFwd.Accept.AddParticle(2112);
   hcalFwd.Accept.AddParticle(-2112);
   hcalFwd.Accept.AddParticle(130);
   hcalMid.Accept.AddParticle(2112);
   hcalMid.Accept.AddParticle(-2112);
   hcalMid.Accept.AddParticle(130);

   hcalBck.Accept.AddZone(hBck);
   hcalFwd.Accept.AddZone(hFwd);
   hcalMid.Accept.AddZone(hMid);

   // Assign acceptance to calorimeters for charged hadrons past tracker
   hcalTrkBck.Accept.SetGenre(Smear::kHadronic);
   hcalTrkFwd.Accept.SetGenre(Smear::kHadronic);

   hcalTrkBck.Accept.SetCharge(Smear::kCharged);
   hcalTrkFwd.Accept.SetCharge(Smear::kCharged);
   
   hcalTrkBck.Accept.AddZone(hTrkBck);
   hcalTrkFwd.Accept.AddZone(hTrkFwd);

   // Assign acceptance to tracker
   momentum.Accept.SetGenre(Smear::kHadronic);
   trackTheta.Accept.SetGenre(Smear::kHadronic);
   trackPhi.Accept.SetGenre(Smear::kHadronic);

   momentum.Accept.SetCharge(Smear::kCharged);
   trackTheta.Accept.SetCharge(Smear::kCharged);
   trackPhi.Accept.SetCharge(Smear::kCharged);

   momentum.Accept.AddZone(trk);
   trackTheta.Accept.AddZone(trk);
   trackPhi.Accept.AddZone(trk);

   // Assign acceptance for calorimeter momentum
   momentumEM.Accept.SetGenre(Smear::kElectromagnetic);
   trackThetaEM.Accept.SetGenre(Smear::kElectromagnetic);
   trackPhiEM.Accept.SetGenre(Smear::kElectromagnetic);

   momentumEM.Accept.AddZone(emTot);
   trackThetaEM.Accept.AddZone(emTot);
   trackPhiEM.Accept.AddZone(emTot);

   momentumHad.Accept.AddParticle(2112);
   trackThetaHad.Accept.AddParticle(2112);
   trackPhiHad.Accept.AddParticle(2112);

   momentumHad.Accept.AddParticle(-2112);
   trackThetaHad.Accept.AddParticle(-2112);
   trackPhiHad.Accept.AddParticle(-2112);

   momentumHad.Accept.AddParticle(130);
   trackThetaHad.Accept.AddParticle(130);
   trackPhiHad.Accept.AddParticle(130);

   momentumHad.Accept.AddZone(hTot);
   trackThetaHad.Accept.AddZone(hTot);
   trackPhiHad.Accept.AddZone(hTot);

   // Assign acceptance for calorimeter momentum for charged hadrons past tracker
   momentumHadTrkBck.Accept.SetGenre(Smear::kHadronic);
   trackThetaHadTrkBck.Accept.SetGenre(Smear::kHadronic);
   trackPhiHadTrkBck.Accept.SetGenre(Smear::kHadronic);

   momentumHadTrkBck.Accept.SetCharge(Smear::kCharged);
   trackThetaHadTrkBck.Accept.SetCharge(Smear::kCharged);
   trackPhiHadTrkBck.Accept.SetCharge(Smear::kCharged);

   momentumHadTrkFwd.Accept.SetGenre(Smear::kHadronic);
   trackThetaHadTrkFwd.Accept.SetGenre(Smear::kHadronic);
   trackPhiHadTrkFwd.Accept.SetGenre(Smear::kHadronic);

   momentumHadTrkFwd.Accept.SetCharge(Smear::kCharged);
   trackThetaHadTrkFwd.Accept.SetCharge(Smear::kCharged);
   trackPhiHadTrkFwd.Accept.SetCharge(Smear::kCharged);

   momentumHadTrkBck.Accept.AddZone(hTrkBck);
   trackThetaHadTrkBck.Accept.AddZone(hTrkBck);
   trackPhiHadTrkBck.Accept.AddZone(hTrkBck);

   momentumHadTrkFwd.Accept.AddZone(hTrkFwd);
   trackThetaHadTrkFwd.Accept.AddZone(hTrkFwd);
   trackPhiHadTrkFwd.Accept.AddZone(hTrkFwd);

   // Assign acceptance for track energy
   trackEnergy.Accept.SetGenre(Smear::kHadronic);

   trackEnergy.Accept.SetCharge(Smear::kCharged);

   trackEnergy.Accept.AddZone(trk);


   //emcal.Accept.AddZone(central);
   //momentum.Accept.AddZone(central);
   //theta.Accept.AddZone(central);
   //phi.Accept.AddZone(central);
   ////rich.Accept.AddZone(central);

   // Create a detector and add the devices
   Smear::Detector det;
   det.AddDevice(emcalBck);
   det.AddDevice(emcalMidBck);
   det.AddDevice(emcalMid);
   det.AddDevice(emcalFwd);
   det.AddDevice(hcalBck);
   det.AddDevice(hcalFwd);
   det.AddDevice(hcalMid);
   det.AddDevice(hcalTrkBck);
   det.AddDevice(hcalTrkFwd);
   det.AddDevice(momentum);
   det.AddDevice(trackTheta);
   det.AddDevice(trackPhi);
   det.AddDevice(momentumEM);
   det.AddDevice(trackThetaEM);
   det.AddDevice(trackPhiEM);
   det.AddDevice(momentumHad);
   det.AddDevice(trackThetaHad);
   det.AddDevice(trackPhiHad);
   det.AddDevice(momentumHadTrkBck);
   det.AddDevice(trackThetaHadTrkBck);
   det.AddDevice(trackPhiHadTrkBck);
   det.AddDevice(momentumHadTrkFwd);
   det.AddDevice(trackThetaHadTrkFwd);
   det.AddDevice(trackPhiHadTrkFwd);
   det.AddDevice(trackEnergy);
   det.SetEventKinematicsCalculator("NM JB DA"); // The detector will calculate event kinematics from smeared values

   //det.AddDevice(emcal);
   //det.AddDevice(momentum);
   //det.AddDevice(theta);
   //det.AddDevice(phi);
   ////det.AddDevice(rich);


   

   return det;
}
