#define MCRecoXS_cxx
#include "MCRecoXS.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1.h>

// Primary Particle Start Positions
TH1D *hMCPrimaryStartX = new TH1D("hMCPrimaryStartX", "Primary Particle X_{0}", 200, -50 , 50);
TH1D *hMCPrimaryStartY = new TH1D("hMCPrimaryStartY", "Primary Particle Y_{0}", 200, -50 , 50);
TH1D *hMCPrimaryStartZ = new TH1D("hMCPrimaryStartZ", "Primary Particle Z_{0}", 600, -120 , 180);

// Primary Projected Particle Start Positions
TH1D *hMCPrimaryProjectedStartX = new TH1D("hMCPrimaryProjectedStartX", "Primary Particle X_{0}", 200, -50 , 50);
TH1D *hMCPrimaryProjectedStartY = new TH1D("hMCPrimaryProjectedStartY", "Primary Particle Y_{0}", 200, -50 , 50);
TH1D *hMCPrimaryProjectedStartZ = new TH1D("hMCPrimaryProjectedStartZ", "Primary Particle Z_{0}", 400, -50 , 150);

// Primary Particle End Positions 
TH1D *hMCPrimaryEndX = new TH1D("hMCPrimaryEndX", "Primary Particle X_{f}", 400, -200 , 200);
TH1D *hMCPrimaryEndY = new TH1D("hMCPrimaryEndY", "Primary Particle Y_{f}", 400, -200 , 200);
TH1D *hMCPrimaryEndZ = new TH1D("hMCPrimaryEndZ", "Primary Particle Z_{f}", 600, -120 , 480);

// Primary Particle P
TH1D *hMCPrimaryPx = new TH1D("hMCPrimaryPx", "Primary Particle P_{x}", 300, -150 , 150);
TH1D *hMCPrimaryPy = new TH1D("hMCPrimaryPy", "Primary Particle P_{y}", 300, -159 , 150);
TH1D *hMCPrimaryPz = new TH1D("hMCPrimaryPz", "Primary Particle P_{z}", 3000, -500 , 2500);

// Primary Particle Process
TH1D *hMCPrimaryProcess = new TH1D("hMCPrimaryProcess", "Primary Particle Process", 100, 0 , 13);

// Primary End X vs Z and Y vs Z Positions
TH2D *hMCPrimaryEndXvsZ = new TH2D("hMCPrimaryEndXvsZ", "X_{f} vs Z_{f}", 600, -150, 450, 400, -200, 200);
TH2D *hMCPrimaryEndYvsZ = new TH2D("hMCPrimaryEndYvsZ", "Y_{f} vs Z_{f}", 600, -150, 450, 200, -200, 200);

// Energy Loss in Upstream part of beamline
TH1D *hMCELossUpstream = new TH1D("hMCELossUpstream", "Energy loss prior to entering the TPC", 1000, 0, 1000);

// True Length
TH1D *hTrueLength = new TH1D("hTrueLength", "#True Length of the Primary Particle inside the TPC", 400, 0 , 200);

TH1D *MPV = new TH1D("MPV", "MPV", 200, 0., 10.);
TH2D *MPVvP = new TH2D("MPVvP", "MPVvp", 200, 0., 10., 100, 0., 2000.);
TH2D *dEdxvP = new TH2D("dEdxvP", "dEdxvP", 200, 0., 10., 600., 0., 1000.);

TH1D *PIDA = new TH1D("PIDA", "PIDA", 100, 0. , 25.);
TH1D *PIDA2 = new TH1D("PIDA2", "PIDA2", 100, 0. , 25.);

TH1F *fMagPCross = new TH1F("fMagPCross", "fMagPCross", 1000, 0, 1.0); 

TH1F *fDiffDist = new TH1F("fDiffDist", "fDiffDist", 200, -100.0, 100.0); 

// Reconstructed Information 

// Most upstream Z point of tracks
TH1D *hdataUpstreamZPos = new TH1D("hdataUpstreamZPos", "Most upstream spacepoint of all TPC Tracks", 20, 0, 10);

// Num of tracks in the TPC versus distance
TH2D *hdataNTracksvsZpos = new TH2D("hdataNTracksvsZpos", "Number of TPC tracks vs Z ", 30, 0, 30, 20, 0, 10);

// Alpha between WC and TPC tracks
TH1D *hAlpha = new TH1D("hAlpha", "#alpha between MC Particle and TPC Track", 90, 0, 90);

// Delta Start Positions
TH1D *hDeltaX = new TH1D("hDeltaX", "#Delta X_{0} of the most upstream Reco Track and the Projected Primary Particle X_{0}", 200, -50 , 50);
TH1D *hDeltaY = new TH1D("hDeltaY", "#Delta Y_{0} of the most upstream Reco Track and the Projected Primary Particle Y_{0}", 200, -50 , 50);
TH1D *hDeltaZ = new TH1D("hDeltaZ", "#Delta Z_{0} of the most upstream Reco Track and the Projected Primary Particle Z_{0}", 200, -50 , 50);

// Reconstructed length
TH1D *hRecoLength = new TH1D("hRecoLength", "#Reconstructed Length of the Primary Particle inside the TPC", 200, 0 , 100);

// delta length position
TH1D *hDeltaLength = new TH1D("hDeltaLength", "#Delta Length of the most upstream Reco Track and the Primary Particle Z_{0}", 200, -100 , 100);


// Cross section and Calorimetry information

/////////////////////////////////// "Kaon" initial Kinetic Energy (MeV) //////////////////////////////////////////
TH1D *hMCInitalKE = new TH1D("hMCInitalKE", "Kaon Initial Kinetic Energy (MeV)", 500, 0, 2500);
/////////////////////////////////// "Kaon" dE/dX //////////////////////////////////////////
TH1D *hdataKaondEdX = new TH1D("hdataKaondEdX", "Kaon dE/dX", 200, 0, 50);
/////////////////////////////////// "Kaon" Residual Range //////////////////////////////////////////
TH1D *hdataKaonRR = new TH1D("hdataKaonRR", "Kaon Residual Range", 240, -10, 110);
/////////////////////////////////// "Kaon" Track Pitch //////////////////////////////////////////
TH1D *hdataKaonTrkPitch = new TH1D("hdataKaonTrkPitch", "Track Pitch", 1000, 0, 5);
///////////////////////////////// "Kaon dE/dX vs RR ///////////////////////////////////////////
TH2D *hdataKaondEdXvsRR = new TH2D("", "dE/dX vs Residual Range",200, 0, 100, 200, 0, 50);
/////////////////////////////////// "Kaon" Incident to the slab Kinetic Energy (MeV) //////////////////////////////////////////
TH1D *hdataKaonIncidentKE = new TH1D("hdataKaonIncidentKE", "Kaon Incident Kinetic Energy (MeV)", 40, 0, 2000);
/////////////////////////////////// "Kaon" Exiting the slab Kinetic Energy (MeV) //////////////////////////////////////////
TH1D *fKaonInteractions = new TH1D("fKaonInteractions", "Kaon Out Kinetic Energy (MeV)", 40, 0, 2000);

TH1D *fEndKE = new TH1D("fEndKE", "EndKE", 200, -500, 500);

/////////////////////////////////// Cross-Section //////////////////////////////////////////////////////////////////////
TH1F *fCrossSection = new TH1F("fCrossSection", "Cross-Section", 20, 0, 1000); 

// Kaon Track Start and End Positions
TH1D *hdataKaonTrackEndX = new TH1D("hdataKaonTrackEndX", "Kaon Track End X Position", 50, 0, 50);
TH1D *hdataKaonTrackEndY = new TH1D("hdataKaonTrackEndY", "Kaon Track End Y Position", 40, -20, 20);
TH1D *hdataKaonTrackEndZ = new TH1D("hdataKaonTrackEndZ", "Kaon Track End Z Position", 100, 0, 100);

TH1D *hdataKaonTrackStartX = new TH1D("hdataKaonTrackStartX", "Kaon Track Start X Position", 50, 0, 50);
TH1D *hdataKaonTrackStartY = new TH1D("hdataKaonTrackStartY", "Kaon Track Start Y Position", 40, -20, 20);
TH1D *hdataKaonTrackStartZ = new TH1D("hdataKaonTrackStartZ", "Kaon Track Start Z Position", 100, 0, 100);


// Delta End Positions
TH1D *hDeltaEndX = new TH1D("hDeltaEndX", "#Delta X_{f} of the most upstream Reco Track and the Projected Primary Particle X_{f}", 200, -50 , 50);
TH1D *hDeltaEndY = new TH1D("hDeltaEndY", "#Delta Y_{f} of the most upstream Reco Track and the Projected Primary Particle Y_{f}", 200, -50 , 50);
TH1D *hDeltaEndZ = new TH1D("hDeltaEndZ", "#Delta Z_{f} of the most upstream Reco Track and the Projected Primary Particle Z_{f}", 200, -50 , 50);


// Stacked histograms. Delta End positions.

// Inelastic
TH1D *hDeltaEndZInElastic = new TH1D("hDeltaEndZInElastic", "#Delta Z_{f} InElastic", 200, -50 , 50);
TH1D *hDeltaEndYInElastic = new TH1D("hDeltaEndYInElastic", "#Delta Y_{f} InElastic", 200, -50 , 50);
TH1D *hDeltaEndXInElastic = new TH1D("hDeltaEndXInElastic", "#Delta X_{f} InElastic", 200, -50 , 50);
// Neutron inelastic
TH1D *hDeltaEndZNeutronInElastic = new TH1D("hDeltaEndZNeutronInElastic", "#Delta Z_{f} Neutron InElastic", 200, -50 , 50);
TH1D *hDeltaEndYNeutronInElastic = new TH1D("hDeltaEndYNeutronInElastic", "#Delta Y_{f} Neutron InElastic", 200, -50 , 50);
TH1D *hDeltaEndXNeutronInElastic = new TH1D("hDeltaEndXNeutronInElastic", "#Delta X_{f} Neutron InElastic", 200, -50 , 50);
// HadElastic
TH1D *hDeltaEndZHadElastic = new TH1D("hDeltaEndZHadElastic", "#Delta Z_{f} Hadronic Elastic", 200, -50 , 50);
TH1D *hDeltaEndYHadElastic = new TH1D("hDeltaEndYHadElastic", "#Delta Y_{f} Hadronic Elastic", 200, -50 , 50);
TH1D *hDeltaEndXHadElastic = new TH1D("hDeltaEndXHadElastic", "#Delta X_{f} Hadronic Elastic", 200, -50 , 50);
// Neutron Capture
TH1D *hDeltaEndZnCap = new TH1D("hDeltaEndZnCap", "#Delta Z_{f} Neutron Capture", 200, -50 , 50);
TH1D *hDeltaEndYnCap = new TH1D("hDeltaEndYnCap", "#Delta Y_{f} Neutron Capture", 200, -50 , 50);
TH1D *hDeltaEndXnCap = new TH1D("hDeltaEndXnCap", "#Delta X_{f} Neutron Capture", 200, -50 , 50);
// Nuclear capture at rest
TH1D *hDeltaEndZnuclearCapatureAtRest = new TH1D("hDeltaEndZnuclearCapatureAtRest", "#Delta Z_{f} Nuclear capture at rest ", 200, -50 , 50);
TH1D *hDeltaEndYnuclearCapatureAtRest = new TH1D("hDeltaEndYnuclearCapatureAtRest", "#Delta Y_{f} Nuclear capture at rest ", 200, -50 , 50);
TH1D *hDeltaEndXnuclearCapatureAtRest = new TH1D("hDeltaEndXnuclearCapatureAtRest", "#Delta X_{f} Nuclear capture at rest ", 200, -50 , 50);
// Decay
TH1D *hDeltaEndZDecay = new TH1D("hDeltaEndZDecay", "#Delta Z_{f} Decay ", 200, -50 , 50);
TH1D *hDeltaEndYDecay = new TH1D("hDeltaEndYDecay", "#Delta Y_{f} Decay ", 200, -50 , 50);
TH1D *hDeltaEndXDecay = new TH1D("hDeltaEndXDecay", "#Delta X_{f} Decay ", 200, -50 , 50);
// Kaon zero inelastic
TH1D *hDeltaEndZKaonZeroInElastic = new TH1D("hDeltaEndZKaonZeroInElastic", "#Delta Z_{f} Kaon Zero InElastic ", 200, -50 , 50);
TH1D *hDeltaEndYKaonZeroInElastic = new TH1D("hDeltaEndYKaonZeroInElastic", "#Delta Y_{f} Kaon Zero InElastic ", 200, -50 , 50);
TH1D *hDeltaEndXKaonZeroInElastic = new TH1D("hDeltaEndXKaonZeroInElastic", "#Delta X_{f} Kaon Zero InElastic ", 200, -50 , 50);
// CoulombScat
TH1D *hDeltaEndZCoulombScat = new TH1D("hDeltaEndZCoulombScat", "#Delta Z_{f} Coulomb Scattering ", 200, -50 , 50);
TH1D *hDeltaEndYCoulombScat = new TH1D("hDeltaEndYCoulombScat", "#Delta Y_{f} Coulomb Scattering ", 200, -50 , 50);
TH1D *hDeltaEndXCoulombScat = new TH1D("hDeltaEndXCoulombScat", "#Delta X_{f} Coulomb Scattering ", 200, -50 , 50);
// Mu- Capture
TH1D *hDeltaEndZMuMinusCapture = new TH1D("hDeltaEndZMuMinusCapture", "#Delta Z_{f} MuMinus Capture ", 200, -50 , 50);
TH1D *hDeltaEndYMuMinusCapture = new TH1D("hDeltaEndYMuMinusCapture", "#Delta Y_{f} MuMinus Capture ", 200, -50 , 50);
TH1D *hDeltaEndXMuMinusCapture = new TH1D("hDeltaEndXMuMinusCapture", "#Delta X_{f} MuMinus Capture ", 200, -50 , 50);
// Proton Inelastic
TH1D *hDeltaEndZProtonInelastic = new TH1D("hDeltaEndZProtonInelastic", "#Delta Z_{f} Proton Inelastic ", 200, -50 , 50);
TH1D *hDeltaEndYProtonInelastic = new TH1D("hDeltaEndYProtonInelastic", "#Delta Y_{f} Proton Inelastic ", 200, -50 , 50);
TH1D *hDeltaEndXProtonInelastic = new TH1D("hDeltaEndXProtonInelastic", "#Delta X_{f} Proton Inelastic ", 200, -50 , 50);
// K+ Inelastic
TH1D *hDeltaEndZKaonPlusInelastic = new TH1D("hDeltaEndZKaonPlusInelastic", "#Delta Z_{f} K+ Inelastic ", 200, -50 , 50);
TH1D *hDeltaEndYKaonPlusInelastic = new TH1D("hDeltaEndYKaonPlusInelastic", "#Delta Y_{f} K+ Inelastic ", 200, -50 , 50);
TH1D *hDeltaEndXKaonPlusInelastic = new TH1D("hDeltaEndXKaonPlusInelastic", "#Delta X_{f} K+ Inelastic ", 200, -50 , 50);

//---------------------------------------------------------------------------------------------------------------------------------------

float particle_mass = 493.67;

float rho = 1400; //kg/m^3
//  float cm_per_m = 100;
float molar_mass = 39.9; //g/mol
float g_per_kg = 1000; 
float avogadro = 6.02e+23; //number/mol
float number_density = rho*g_per_kg/molar_mass*avogadro;
float slab_width = 0.0047;//in m

double FirstSpacePointZPos = 2.0;

double DeltaXLowerBound = -1.5;
double DeltaXUpperBound = 1.5;

double DeltaYLowerBound = -1.8;
double DeltaYUpperBound = 1.8;

double XLowerFid = 0;
double XUpperFid = 47;

double YLowerFid = -20;
double YUpperFid = 20;

double ZLowerFid = 0;
double ZUpperFid = 87;

int UpperPartOfTPC = 14.0;
int nLowZTracksAllowed = 4;

double alphaCut = 10;   
float EventWeight = 1.0;

// #################################################   
// ### True  = Use the momentum based weighting  ###
// ### False = Don't weight events               ###
// #################################################
bool UseEventWeight = false;

int agree = 0;
int tot = 0;

int nTotalEvents = 0, nEvtsGoodMC = 0, nEvtsTrackZPos = 0, nEvtsMCTrackMatch = 0;
int nEventsPassingAlpha = 0, nLowZTrkEvents = 0;

int counter = 0;


//std::vector< std::vector<double> > pos;
std::vector<double> posx;
std::vector<double> posy;
std::vector<double> posz;


void MCRecoXS::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L MCRecoXS.C
//      Root > MCRecoXS t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch

  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();

  Long64_t nbytes = 0, nb = 0;
  nentries = 5000;
  for (Long64_t jentry=0; jentry < nentries; jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;


    // Super idiot test, make sure we have Kaons only 
    if(pdg[0] != 321) { printf(" What? \n"); break; }

    // 
    nTotalEvents++;
    if(nTotalEvents%1000 == 0.) { printf("Ping %i) Rn:%i Sb:%i \n", nTotalEvents, run, event); }

    // GEANT 4 Information

    // ### Defining some useful variables we will use ###
    int nG4Primary = 0;
    int nG4TrajPoints = 0;
   
    float g4Primary_X0[100] = {0.}, g4Primary_Y0[100] = {0.}, g4Primary_Z0[100] = {0.};
    float g4Primary_ProjX0[100] = {0.}, g4Primary_ProjY0[100] = {0.}, g4Primary_ProjZ0[100] = {0.};
   
    float g4Primary_Xf[100] = {0.}, g4Primary_Yf[100] = {0.}, g4Primary_Zf[100] = {0.};
    float g4Primary_Px[100] = {0.}, g4Primary_Py[100] = {0.}, g4Primary_Pz[100] = {0.};
   
    int g4Primary_TrkID[100] = {999}, g4PrimaryProcess[100] = {0};
   
    int nG4PriTrj = 0;
    float g4Primary_TrueTrjX[10][10000] = {0.}, g4Primary_TrueTrjY[10][10000] = {0.}, g4Primary_TrueTrjZ[10][10000] = {0.};
    float g4Primary_TrueTrjPx[10][10000] = {0.}, g4Primary_TrueTrjPy[10][10000] = {0.}, g4Primary_TrueTrjPz[10][10000] = {0.};
   
   
    float TrueLength = 0;
    float RecoLength = 0;
    // ##########################################
    // ### Loop over all the GEANT4 particles ###
    // ##########################################
    for (int iG4 = 0; iG4 < geant_list_size; iG4++) {
      // #####################################################
      // ### If this is a primary particle then look at it ###
      // #####################################################
      if(process_primary[iG4] == 1) {

	 
	// ### Recording information for use later ###
	g4Primary_X0[nG4Primary] = StartPointx[iG4];
	g4Primary_Y0[nG4Primary] = StartPointy[iG4];
	g4Primary_Z0[nG4Primary] = StartPointz[iG4];
	 
	//	  printf("lol wut: %i \n\n", pdg[iG4]);   //[geant_list_size]

	g4Primary_Xf[nG4Primary] = EndPointx[iG4];
	g4Primary_Yf[nG4Primary] = EndPointy[iG4];
	g4Primary_Zf[nG4Primary] = EndPointz[iG4];
	 
	g4Primary_Px[nG4Primary] = Px[iG4] * 1000; //<---Converting to MeV
	g4Primary_Py[nG4Primary] = Py[iG4] * 1000; //<---Converting to MeV
	g4Primary_Pz[nG4Primary] = Pz[iG4] * 1000; //<---Converting to MeV
	 
	TrueLength = sqrt( ((EndPointz[iG4]-StartPointz[iG4])*(EndPointz[iG4]-StartPointz[iG4])) + 
			   ((EndPointy[iG4]-StartPointy[iG4])*(EndPointy[iG4]-StartPointy[iG4])) + 
			   ((EndPointx[iG4]-StartPointx[iG4])*(EndPointx[iG4]-StartPointx[iG4])) );

	hTrueLength->Fill(TrueLength);
	 
	// ------------------------------------------------------------------------------------
	// ---------------        Extrapolate the X, Y, Z position of the primary         -----
	// ---------------     if it started upstream of the front face of the TPC        -----
	// ------------------------------------------------------------------------------------
	 
	double b1 = StartPointz[iG4] - StartPointx[iG4]*Pz[iG4]/Px[iG4];
	double b2 = StartPointz[iG4] - StartPointy[iG4]*Pz[iG4]/Py[iG4];
	 
	g4Primary_ProjX0[nG4Primary] = -b1*Px[iG4]/Pz[iG4];
	g4Primary_ProjY0[nG4Primary] = -b2*Py[iG4]/Pz[iG4];
	g4Primary_ProjZ0[nG4Primary] = 0.0;
	 
	// ### Setting the primary particles Track ID ###
	g4Primary_TrkID[nG4Primary] = TrackId[iG4];

	hMCPrimaryEndXvsZ->Fill(EndPointz[iG4], EndPointx[iG4]);
	hMCPrimaryEndYvsZ->Fill(EndPointz[iG4], EndPointy[iG4]);
	 
	hMCPrimaryPx->Fill(g4Primary_Px[nG4Primary], EventWeight);
	hMCPrimaryPy->Fill(g4Primary_Py[nG4Primary], EventWeight);
	hMCPrimaryPz->Fill(g4Primary_Pz[nG4Primary], EventWeight);

	hMCPrimaryStartX->Fill(StartPointx[iG4]);
	hMCPrimaryStartY->Fill(StartPointy[iG4]);
	hMCPrimaryStartZ->Fill(StartPointz[iG4]);
	 
	hMCPrimaryEndX->Fill(EndPointx[iG4]);
	hMCPrimaryEndY->Fill(EndPointy[iG4]);
	hMCPrimaryEndZ->Fill(EndPointz[iG4]);

	hMCPrimaryProjectedStartX->Fill(g4Primary_ProjX0[nG4Primary]);
	hMCPrimaryProjectedStartY->Fill(g4Primary_ProjY0[nG4Primary]);
	hMCPrimaryProjectedStartZ->Fill(g4Primary_ProjZ0[nG4Primary]);
		 
	nG4TrajPoints = 0;
	 
	// ### Recording the primary particles trajectory points ###
	for(int iG4Tr = 0; iG4Tr < NTrTrajPts[iG4]; iG4Tr++) {
	  g4Primary_TrueTrjX[nG4Primary][nG4PriTrj] = MidPosX[iG4][iG4Tr];
	  g4Primary_TrueTrjY[nG4Primary][nG4PriTrj] = MidPosY[iG4][iG4Tr];
	  g4Primary_TrueTrjZ[nG4Primary][nG4PriTrj] = MidPosZ[iG4][iG4Tr];
	    
	  //std::cout<<"g4Primary_TrueTrjZ[nG4Primary][nG4PriTrj] = "<<g4Primary_TrueTrjZ[nG4Primary][nG4PriTrj]<<std::endl;
	    
	  g4Primary_TrueTrjPx[nG4Primary][nG4PriTrj] = MidPx[iG4][iG4Tr]*1000;//<---Converting to MeV
	  g4Primary_TrueTrjPy[nG4Primary][nG4PriTrj] = MidPy[iG4][iG4Tr]*1000;//<---Converting to MeV
	  g4Primary_TrueTrjPz[nG4Primary][nG4PriTrj] = MidPz[iG4][iG4Tr]*1000;//<---Converting to MeV
	    
	  nG4PriTrj++;
	}//<---end looping over this primary particles true trajectory points
	 
	 
	// ### Bumping the counters ###
	nG4Primary++;
	 
      }//<---End Looking only at primaries
    }// <---End iG4 Loop

    hMCPrimaryProcess->Fill(g4PrimaryProcess[nG4Primary - 1]);
   
    // ################################################
    // ### Loop over all the GEANT4 particles again ###
    // ###  to get the process from the daughters   ###
    // ################################################
    for (int iG4 = 0; iG4 < geant_list_size; iG4++) {      
      // ### Looking for the Daughters of the primary ###
      if(Mother[iG4] == g4Primary_TrkID[nG4Primary - 1]) {
	g4PrimaryProcess[nG4Primary - 1] = Process[iG4];
      }//<---End matching daughters      
    }//<---end iG4 loop
        
    //====================================================
    // Only looking at events where the primary particle enters the TPC
    //====================================================
   
    bool GoodMCEventInTPC = true;
   
    // ##############################################
    // ### Looping over all the primary particles ###
    // ##############################################
    for(int npri = 0; npri < nG4Primary; npri++) {
      if(g4Primary_Zf[npri] < 0) { GoodMCEventInTPC = false; }
      
      // ####################################################################
      // ### Calculating the energy loss for particles that enter the TPC ###
      // ####################################################################
      if(GoodMCEventInTPC) {
	float DifferenceInEnergy = 0;
	// ### Loop over the true trajectory points ###
	for(int ntrj = 0; ntrj < nG4PriTrj; ntrj++) {
	  // ### Only looking at point which are upstream of the TPC ###
	  if(g4Primary_TrueTrjZ[npri][ntrj] < 0) {
	       
	    float Momentum_Point1 = sqrt((g4Primary_TrueTrjPx[npri][ntrj]*g4Primary_TrueTrjPx[npri][ntrj]) + 
					 (g4Primary_TrueTrjPy[npri][ntrj]*g4Primary_TrueTrjPy[npri][ntrj]) +
					 (g4Primary_TrueTrjPz[npri][ntrj]*g4Primary_TrueTrjPz[npri][ntrj]));
				       
	    float Momentum_Point2 = sqrt((g4Primary_TrueTrjPx[npri][ntrj+1]*g4Primary_TrueTrjPx[npri][ntrj+1]) + 
					 (g4Primary_TrueTrjPy[npri][ntrj+1]*g4Primary_TrueTrjPy[npri][ntrj+1]) +
					 (g4Primary_TrueTrjPz[npri][ntrj+1]*g4Primary_TrueTrjPz[npri][ntrj+1]));
				       
	    float Energy_Point1 = sqrt( (Momentum_Point1*Momentum_Point1) + (particle_mass*particle_mass)  );
	       
	    float Energy_Point2 = sqrt( (Momentum_Point2*Momentum_Point2) + (particle_mass*particle_mass)  );
	       
	    DifferenceInEnergy +=  Energy_Point1 - Energy_Point2;
	       
	    //std::cout<<"z = "<<g4Primary_TrueTrjZ[npri][ntrj]<<", DifferenceInEnergy = "<<DifferenceInEnergy<<std::endl;	       	       
	  }//<---End only look at points which are upstream of the TPC	    	    
	}//<---End ntrj for loop	 	 
	hMCELossUpstream->Fill(DifferenceInEnergy);	  
      }//<---Only looking at events that actually make it into the TPC            
    }//<---End npri loop
   
    //if(!GoodMCEventInTPC){continue;}
   
    if(EndPointz[0] < ZLowerFid) { continue; }
    nEvtsGoodMC++;
   
   
    //====================================================
    //						Low Z Spacepoint Track Cut
    //====================================================
   
    // ### Boolian for events w/ track which ###
    // ###     starts at the front face      ###
    bool TrackTrjPtsZCut = false;
   
    // ### Recording the index of the track which ###
    // ###   starts at the front face of the TPC  ###
    bool PreLimTrackIndex[500] = {false};
   
    // ##################################################
    // ### Defining a dummy variable used for sorting ###
    // ##################################################
    double dummyXpoint = 999, dummyYpoint = 999, dummyZpoint = 100, dummyPointTrkInd = -1;
    double dummypoint_TempTrjX = 999, dummypoint_TempTrjY = 999, dummypoint_TempTrjZ = 999;
    double dummyTrkX[200] = {0.}, dummyTrkY[200] = {0.}, dummyTrkZ[200] = {100.};
    double dummyTrk_pHat0X[200] = {0.}, dummyTrk_pHat0Y[200] = {0.}, dummyTrk_pHat0Z[200] = {0.};
    double dummyTrk_Theta[200] = {0.}, dummyTrk_Phi[200] = {0.};
    double dummyTrk_Index[200] = {0.};
   
    int nUpStreamTrk = 0;
   
    // ###########################
    // ### Looping over tracks ###
    // ###########################
    for(int iTrk = 0; iTrk < ntracks_reco; iTrk++) {
      
      float tempZpoint = 100;
      // ########################################################
      // ### Looping over the trajectory points for the track ###
      // ########################################################
      for(int iTrjPt = 0; iTrjPt < nTrajPoint[iTrk]; iTrjPt++) {
	 
	// ################################################################################ 
	// ### Tracking the lowest Z point that is inside fiducial boundries of the TPC ###
	// ################################################################################
	// ### Resetting the variables for each track ###
	dummyXpoint = 999, dummyYpoint = 999, dummyZpoint = 999;
	 
	// ###########################################################################
	// ### Setting our dummypoints if this is the lowest Z point on this track ###
	// ###           and still within the active volume of the TPC             ###
	// ###########################################################################
	if(trjPt_Z[iTrk][iTrjPt] < tempZpoint && trjPt_Z[iTrk][iTrjPt] > ZLowerFid && //<---Note: the variable "trjPt_Z[iTrk][iTrjPt]"
	   trjPt_Y[iTrk][iTrjPt] > YLowerFid && trjPt_Y[iTrk][iTrjPt] < YUpperFid &&       // is the z position of the 3d point for the iTrk
	   trjPt_X[iTrk][iTrjPt] > XLowerFid && trjPt_X[iTrk][iTrjPt] < XUpperFid ) {	     
	  dummyXpoint = trjPt_X[iTrk][iTrjPt];
	  dummyYpoint = trjPt_Y[iTrk][iTrjPt];
	  dummyZpoint = trjPt_Z[iTrk][iTrjPt];
	    
	  dummypoint_TempTrjX = pHat0_X[iTrk][iTrjPt];
	  dummypoint_TempTrjY = pHat0_Y[iTrk][iTrjPt];
	  dummypoint_TempTrjZ = pHat0_Z[iTrk][iTrjPt];
	    
	  dummyPointTrkInd = iTrk;
	    
	  tempZpoint = trjPt_Z[iTrk][iTrjPt];
	}//<---End looking for the most upstream point
        
	// ### Only passing events with a track that has ###
	// ###  a spacepoint within the first N cm in Z  ### 
	// ###    And requiring it to be inside the TPC  ###
	if(dummyZpoint < FirstSpacePointZPos) {
	    
	  dummyTrkX[nUpStreamTrk] = dummyXpoint;
	  dummyTrkY[nUpStreamTrk] = dummyYpoint;
	  dummyTrkZ[nUpStreamTrk] = dummyZpoint;
	  
	  dummyTrk_pHat0X[nUpStreamTrk] = dummypoint_TempTrjX;
	  dummyTrk_pHat0Y[nUpStreamTrk] = dummypoint_TempTrjY;
	  dummyTrk_pHat0Z[nUpStreamTrk] = dummypoint_TempTrjZ;
	  dummyTrk_Index[nUpStreamTrk] = dummyPointTrkInd;
	    
	  nUpStreamTrk++;
	    
	  TrackTrjPtsZCut = true;
	}
      }//<---End looping over nspts	
	 
      hdataUpstreamZPos->Fill(tempZpoint);

      
      // ### Recording that this track is a "good Track if ###
      // ###  it has a space point in the first N cm in Z  ###
      if(TrackTrjPtsZCut){ PreLimTrackIndex[iTrk] = true;}
      	 
    }//<---End nTPCtrk loop
      
    // ###############################################
    // ### Skipping events that don't have a track ###
    // ###   in the front of the TPC (Z) Position  ###
    // ###############################################
    if(!TrackTrjPtsZCut) { continue; }
    // ### Counting Events w/ front face TPC Track ###
    nEvtsTrackZPos++; 
      
    //=======================================================================================================================
    //					Cutting on the number of tracks in the upstream TPC
    //=======================================================================================================================
   
    int nLowZTracksInTPC = 0;
    // ################################################################
    // ### Initializing variables for study of low Z track location ###
    // ################################################################
    bool LowZTrackInTPC = false;
   
    float templowz1 = 0;
    float templowz2 = 0;
      
    // #################################################################
    // ### Only keeping events if there is less than N tracks in the ###
    // ###    first ## cm of the TPC (to help cut out EM Showers     ###
    // #################################################################
    for(size_t nTPCtrk = 0; nTPCtrk < ntracks_reco; nTPCtrk++) {
      // ### Start by assuming this track is not in the ###
      // ###          low Z part of the TPC             ###
      LowZTrackInTPC = false;
           

      // ### Looping over the spacepoints for the track ###
      for(size_t nspts = 0; nspts < nTrajPoint[nTPCtrk]; nspts++) {
	 
	// ##################################################
	// ### Count this track if it has a spacepoint in ###
	// ###       the low Z region of the TPC          ###
	// ##################################################
	if(trkz[nTPCtrk][nspts] < UpperPartOfTPC) {
	  if(trky[nTPCtrk][nspts] > YLowerFid && trky[nTPCtrk][nspts] < YUpperFid && 
	     trkx[nTPCtrk][nspts] > XLowerFid && trkx[nTPCtrk][nspts] < XUpperFid)
	    {LowZTrackInTPC = true; templowz1 = trkz[nTPCtrk][nspts];}
		
	}//<---End counting if 
	  
      }//<---End nspts loop
      
      // ##################################################################
      // ### If the track was in the "UpperPartOfTPC", bump the counter ###
      // ##################################################################
      if(LowZTrackInTPC) { nLowZTracksInTPC++; }//<---End counting track in the Upstream part	 	     	    

    }//<---End nTPCtrk
    
        
    // ### Skipping the event if there are too many ###
    // ###       low Z tracks in the event          ###
    if(nLowZTracksInTPC > nLowZTracksAllowed || nLowZTracksInTPC == 0) { continue; }
    
    // ### Counting the event if it passes ###
    nLowZTrkEvents++;

    // Matching the MC Particle to the Reco Track (similar to the WC Track match)      
   
    // ################################################
    // ### Calculating the angles for the Geant4 MC ###
    // ################################################
    TVector3 z_hat_MC(0,0,1);
    TVector3 p_hat_0_MC;
   
    // ### Setting the vector for the MC using the ###
    // ###  extrapolated Momentum vector   ###
    p_hat_0_MC.SetX(g4Primary_Px[0]);
    p_hat_0_MC.SetY(g4Primary_Py[0]);
    p_hat_0_MC.SetZ(g4Primary_Pz[0]); 
   
    // ### Getting everything in the same convention ###
    float mcPhi = 0;
    float mcTheta = 0;
   
    // === Calculating Theta for MC ===
    mcTheta = acos(z_hat_MC.Dot(p_hat_0_MC)/p_hat_0_MC.Mag());
   
    // === Calculating Phi for MC ===
    //---------------------------------------------------------------------------------------------------------------------
    if( p_hat_0_MC.Y() > 0 && p_hat_0_MC.X() > 0 ){ mcPhi = atan(p_hat_0_MC.Y()/p_hat_0_MC.X()); }
    else if( p_hat_0_MC.Y() > 0 && p_hat_0_MC.X() < 0 ){ mcPhi = atan(p_hat_0_MC.Y()/p_hat_0_MC.X())+3.141592654; }
    else if( p_hat_0_MC.Y() < 0 && p_hat_0_MC.X() < 0 ){ mcPhi = atan(p_hat_0_MC.Y()/p_hat_0_MC.X())+3.141592654; }
    else if( p_hat_0_MC.Y() < 0 && p_hat_0_MC.X() > 0 ){ mcPhi = atan(p_hat_0_MC.Y()/p_hat_0_MC.X())+6.28318; }
    else if( p_hat_0_MC.Y() == 0 && p_hat_0_MC.X() == 0 ){ mcPhi = 0; }//defined by convention
    else if( p_hat_0_MC.Y() == 0 ) {
      if( p_hat_0_MC.X() > 0 ){ mcPhi = 0; }
      else{ mcPhi = 3.141592654; }
    } 
    else if( p_hat_0_MC.X() == 0 ) {
      if( p_hat_0_MC.Y() > 0 ){ mcPhi = 3.141592654/2; }
      else{ mcPhi = 3.141592654*3/2; }
    }
    //---------------------------------------------------------------------------------------------------------------------
      
    // ######################################################
    // ### Calculating the angles for the Upstream Tracks ###
    // ######################################################
   
    TVector3 z_hat_TPC(0,0,1);
    TVector3 p_hat_0_TPC;
    for(int aa = 0; aa < nUpStreamTrk; aa++) {
      // ### Setting the TVector ###
      p_hat_0_TPC.SetX(dummyTrk_pHat0X[aa]);
      p_hat_0_TPC.SetY(dummyTrk_pHat0Y[aa]);
      p_hat_0_TPC.SetZ(dummyTrk_pHat0Z[aa]);
           
      // ### Calculating TPC track theta ###
      dummyTrk_Theta[aa] = acos(z_hat_TPC.Dot(p_hat_0_TPC)/p_hat_0_TPC.Mag());
      
      // ### Calculating TPC track phi ###
      //---------------------------------------------------------------------------------------------------------------------
      if( p_hat_0_TPC.Y() > 0 && p_hat_0_TPC.X() > 0 ){ dummyTrk_Phi[aa] = atan(p_hat_0_TPC.Y()/p_hat_0_TPC.X()); }
      else if( p_hat_0_TPC.Y() > 0 && p_hat_0_TPC.X() < 0 ){ dummyTrk_Phi[aa] = atan(p_hat_0_TPC.Y()/p_hat_0_TPC.X())+3.141592654; }
      else if( p_hat_0_TPC.Y() < 0 && p_hat_0_TPC.X() < 0 ){ dummyTrk_Phi[aa] = atan(p_hat_0_TPC.Y()/p_hat_0_TPC.X())+3.141592654; }
      else if( p_hat_0_TPC.Y() < 0 && p_hat_0_TPC.X() > 0 ){ dummyTrk_Phi[aa] = atan(p_hat_0_TPC.Y()/p_hat_0_TPC.X())+6.28318; }
      else if( p_hat_0_TPC.Y() == 0 && p_hat_0_TPC.X() == 0 ){ dummyTrk_Phi[aa] = 0; }//defined by convention
      else if( p_hat_0_TPC.Y() == 0 ) {
	if( p_hat_0_TPC.X() > 0 ){ dummyTrk_Phi[aa] = 0; }
	else{ dummyTrk_Phi[aa] = 3.141592654; }	  
      }
      else if( p_hat_0_TPC.X() == 0 ) {	
	if( p_hat_0_TPC.Y() > 0 ){ dummyTrk_Phi[aa] = 3.141592654/2; }
	else{ dummyTrk_Phi[aa] = 3.141592654*3/2; }
      }
      //---------------------------------------------------------------------------------------------------------------------      
    }//<---End aa loop
   
      
    //Cutting on the DeltaX and Delta Y between MC and Reco Track        
    int nMatches = 0;
    int RecoTrackIndex;
    // ###################################
    // ### Loop over all the US Tracks ###
    // ###################################
    for(int bb = 0; bb < nUpStreamTrk; bb++) {
      float DeltaX = dummyTrkX[bb] - g4Primary_ProjX0[0];
      float DeltaY = dummyTrkY[bb] - g4Primary_ProjY0[0];
      float DeltaZ = dummyTrkZ[bb] - g4Primary_ProjZ0[0];     
      
      hDeltaX->Fill(DeltaX);
      hDeltaY->Fill(DeltaY);
      hDeltaZ->Fill(DeltaZ);

      // ### Matching in Delta X and Delta Y and alpha ###
      if(DeltaX < DeltaXLowerBound || DeltaX > DeltaXUpperBound ||
	 DeltaY < DeltaYLowerBound || DeltaY > DeltaYUpperBound) { continue; }

      // ### Define the unit vectors for the MC and TPC Tracks ###
      TVector3 theUnitVector_MC;
      TVector3 theUnitVector_TPCTrack;
   
      theUnitVector_MC.SetX(sin(mcTheta)*cos(mcPhi));
      theUnitVector_MC.SetY(sin(mcTheta)*sin(mcPhi));
      theUnitVector_MC.SetZ(cos(mcTheta));
   
      theUnitVector_TPCTrack.SetX(sin(dummyTrk_Theta[bb])*cos(dummyTrk_Phi[bb]));
      theUnitVector_TPCTrack.SetY(sin(dummyTrk_Theta[bb])*sin(dummyTrk_Phi[bb]));
      theUnitVector_TPCTrack.SetZ(cos(dummyTrk_Theta[bb]));
      
      // ### Calculating the angle between WCTrack and TPC Track ###
      float alpha = (acos(theUnitVector_MC.Dot(theUnitVector_TPCTrack)) )* (180.0/3.141592654);
      hAlpha->Fill(alpha);

      // ### Setting the boolian for the true match ###
      // Cutting on the Alpha Angle between MC and Reco Track     
      if(alpha > alphaCut) { continue; }
      
      nMatches++;
      RecoTrackIndex = dummyTrk_Index[bb];   
      
    } //<---End bb loop
     
    // If don't have one unique match, get outta here
    if(nMatches != 1) { continue; }
   
    // Tag the Throughgoings
    bool ExitingTrack = false;
    if(trkendx[RecoTrackIndex] < XLowerFid || trkendx[RecoTrackIndex] > XUpperFid || 
       trkendy[RecoTrackIndex] < YLowerFid || trkendy[RecoTrackIndex] > YUpperFid  || 
       trkendz[RecoTrackIndex] > ZUpperFid) { ExitingTrack = true; }

    PIDA->Fill(trkpida[RecoTrackIndex][0]);    
    if(!ExitingTrack) { PIDA2->Fill(trkpida[RecoTrackIndex][0]); }
      
    // Tag the decays
    bool Decay = false;
    if(trkpida[RecoTrackIndex][0] > 12.5) { Decay = true; }

    // So in this context, this means that it is a bad reconstruction (?)
    if(trkpida[RecoTrackIndex][0] > 15.1 && !ExitingTrack) { continue; }


    // ===========================================================================================================================================
    // 						Calculating the cross-section 
    // ===========================================================================================================================================   


    // ### The assumed energy loss between the cryostat and the TPC ###
    float entryTPCEnergyLoss = 48; //MeV

    // ### The assumed mass of the incident particle (here we assume a pion) ###
    float mass = 493.67;

    // #############################################################
    // ### Calculating the momentum from the MC Primary Particle ###
    // #############################################################
    float momentum = sqrt( (g4Primary_Px[0]*g4Primary_Px[0]) + (g4Primary_Py[0]*g4Primary_Py[0]) + (g4Primary_Pz[0]*g4Primary_Pz[0]) );
   
    // ###   Calculating the initial Kinetic Energy    ###
    // ### KE = Energy - mass = (p^2 + m^2)^1/2 - mass ###   
    float kineticEnergy = pow( (momentum*momentum) + (mass*mass) ,0.5) - mass;  
    kineticEnergy -= entryTPCEnergyLoss;
   
    // ### Filling the initial kinetic energy plot ###
    hMCInitalKE->Fill(kineticEnergy);   
   
    // ########################################################################
    // ### Variables for the track we are calculating the cross-section for ###
    // ########################################################################
    double Kaondedx[1000]={0.};
    double Kaonresrange[1000]={0.};
    double Kaonpitchhit[1000]={0.};
    int nKaonSpts = 0;
    double KaonSumEnergy = 0;
   
    // ### Variables for determining the matching of the end point ###
    float TrackEndX = 999, TrackEndY = 999, TrackEndZ = 999;
   
    // ################################
    // ### Loop over all TPC Tracks ###
    // ################################
    for(size_t nTPCtrk = 0; nTPCtrk < ntracks_reco; nTPCtrk++) {
      if(nTPCtrk != RecoTrackIndex){continue;}
      
      // ### Recording the end-point of this track ###
      TrackEndX = trkendx[nTPCtrk];
      TrackEndY = trkendy[nTPCtrk];
      TrackEndZ = trkendz[nTPCtrk];     

      hdataKaonTrackEndX->Fill(TrackEndX);
      hdataKaonTrackEndY->Fill(TrackEndY);
      hdataKaonTrackEndZ->Fill(TrackEndZ);
      
      // ### Recording the start-point of this track ###
      
      hdataKaonTrackStartX->Fill(trkvtxx[nTPCtrk]);
      hdataKaonTrackStartY->Fill(trkvtxy[nTPCtrk]);
      hdataKaonTrackStartZ->Fill(trkvtxz[nTPCtrk]);
      
      RecoLength = trklength[nTPCtrk];
      
      hRecoLength->Fill(RecoLength);
            
      nKaonSpts = 0;
      // ###################################################
      // ### Looping over the spacepoints for this track ###
      // ###################################################
	

      for(size_t nspts = 0; nspts < ntrkhits[nTPCtrk]; nspts++) {

	// ###                 Note: Format for this variable is:             ###
	// ### [trk number][plane 0 = induction, 1 = collection][spts number] ###
	Kaondedx[nKaonSpts]     = trkdedx[nTPCtrk][1][nspts];// * 0.90;//<----Scaling dEdX
	  
	// ### Putting in a fix in the case that the dE/dX is negative in this step ### 
	// ###  then take the point before and the point after and average them
	if(Kaondedx[nKaonSpts] < 0 && nspts < ntrkhits[nTPCtrk] && nspts > 0)
	  {Kaondedx[nKaonSpts] = ( (trkdedx[nTPCtrk][1][nspts - 1] + trkdedx[nTPCtrk][1][nspts + 1]) / 2);}
	 
	// ### If this didn't fix it, then just put in a flat 2.4 MeV / cm fix ###
	if(Kaondedx[nKaonSpts] < 0) {		
	  Kaondedx[nKaonSpts] = 2.4;
	  continue;
	}

	Kaonresrange[nKaonSpts] = trkrr[nTPCtrk][1][nspts];
	Kaonpitchhit[nKaonSpts] = trkpitchhit[nTPCtrk][1][nspts];
         
	KaonSumEnergy = (Kaondedx[nKaonSpts] * Kaonpitchhit[nKaonSpts]) + KaonSumEnergy;	 	  

	// ### Recording the dE/dX ###
	hdataKaondEdX->Fill(Kaondedx[nKaonSpts], EventWeight);
	// ### Recording the residual range ###
	hdataKaonRR->Fill(Kaonresrange[nKaonSpts]);
	// ### Recording the Pitch ###
	hdataKaonTrkPitch->Fill(Kaonpitchhit[nKaonSpts], EventWeight);

	// printf(" %lf \n ", Kaondedx[nKaonSpts] * Kaonpitchhit[nKaonSpts]);
	 
	// ### Filling 2d dE/dX vs RR ###
	hdataKaondEdXvsRR->Fill(Kaonresrange[nKaonSpts], Kaondedx[nKaonSpts]);
	 
	nKaonSpts++;
      } //<---End nspts loop
    } //<---End nTPCtrk loop
   
    bool revCalo = false;
    double ordereddEdx[1000]={0.};
    double orderedRR[1000]={0.};
    double orderedPitch[1000]={0.};

    int dimCalo = nKaonSpts;

    for(size_t caloPoints = 0; caloPoints < nKaonSpts; caloPoints++) {
      //cout <<  Kaondedx[caloPoints] << " " << Kaonresrange[caloPoints] << endl;
      if(caloPoints < nKaonSpts-1) {
	//cout << "pippo" << endl;
	if(Kaonresrange[caloPoints] > Kaonresrange[caloPoints+1]) {
	  //If the previous caloHit RR is higher than the next, that's the starting point of the trac
	  ordereddEdx[caloPoints]=Kaondedx[caloPoints];
	  orderedRR[caloPoints]= Kaonresrange[caloPoints];
	  orderedPitch[caloPoints]=Kaonpitchhit[caloPoints];
	} else {
	
	  revCalo = true;
	  //cout << "pluto" << endl;
	  ordereddEdx[dimCalo-caloPoints]= Kaondedx[caloPoints];
	  orderedRR[dimCalo-caloPoints]= Kaonresrange[caloPoints];
	  orderedPitch[dimCalo-caloPoints]=Kaonpitchhit[caloPoints];
	}
      }
    }//<---End calo points

    for(size_t i = 0; i < nKaonSpts; i++) {
      Kaondedx[i] = ordereddEdx[i];
      Kaonresrange[i] = orderedRR[i];
      Kaonpitchhit[i] = orderedPitch[i];
    }

    // Fixing or vetoing the weird ones
    bool good = true;

    for(size_t i = 0; i < nKaonSpts; i++) {

      if(Kaondedx[i] > 40) { 
	if(i != nKaonSpts-1 && i != 0) { 
	  if(Kaondedx[i+1] < 40 && Kaondedx[i-1] < 40) {     
	    printf(" Before : %lf ", Kaondedx[i]);
	    Kaondedx[i] = (Kaondedx[i-1] + Kaondedx[i+1]) / 2;		
	    printf(" After : %lf \n", Kaondedx[i]);
	  } else { good = false; } // Couldn't be averaged out
	} else { good = false; } // At an end point or start point
      } // Check if the dEdx is too damn high

      if(good == false) { break; }

    } 

    // Veto the event if poorly reconstructed dEdx and can't fix
    if(!good) { continue; }


         
    // ###################################
    // ### Filling the Delta End Point ###
    // ###################################
    float DeltaEndX = g4Primary_Xf[0] - TrackEndX;
    float DeltaEndY = g4Primary_Yf[0] - TrackEndY;
    float DeltaEndZ = g4Primary_Zf[0] - TrackEndZ;
  
    // THROUGH GOINGS MAKE THIS WEIRD

    hDeltaEndX->Fill(DeltaEndX);
    hDeltaEndY->Fill(DeltaEndY);
    hDeltaEndZ->Fill(DeltaEndZ);

    if(DeltaEndZ > 3 || DeltaEndZ < 3){counter++;}
   
    if(g4PrimaryProcess[0] == 1){hDeltaEndZInElastic->Fill(DeltaEndZ);}
    if(g4PrimaryProcess[0] == 2){hDeltaEndZNeutronInElastic->Fill(DeltaEndZ);}
    if(g4PrimaryProcess[0] == 3){hDeltaEndZHadElastic->Fill(DeltaEndZ);}
    if(g4PrimaryProcess[0] == 4){hDeltaEndZnCap->Fill(DeltaEndZ);}
    if(g4PrimaryProcess[0] == 5){hDeltaEndZnuclearCapatureAtRest->Fill(DeltaEndZ);}
    if(g4PrimaryProcess[0] == 6){hDeltaEndZDecay->Fill(DeltaEndZ);}
    if(g4PrimaryProcess[0] == 7){hDeltaEndZKaonZeroInElastic->Fill(DeltaEndZ);}
    if(g4PrimaryProcess[0] == 8){hDeltaEndZCoulombScat->Fill(DeltaEndZ);}
    if(g4PrimaryProcess[0] == 9){hDeltaEndZMuMinusCapture->Fill(DeltaEndZ);}
    if(g4PrimaryProcess[0] == 10){hDeltaEndZProtonInelastic->Fill(DeltaEndZ);}
    if(g4PrimaryProcess[0] == 11){hDeltaEndZKaonPlusInelastic->Fill(DeltaEndZ);}
   
    if(g4PrimaryProcess[0] == 1){hDeltaEndYInElastic->Fill(DeltaEndY);}
    if(g4PrimaryProcess[0] == 2){hDeltaEndYNeutronInElastic->Fill(DeltaEndY);}
    if(g4PrimaryProcess[0] == 3){hDeltaEndYHadElastic->Fill(DeltaEndY);}
    if(g4PrimaryProcess[0] == 4){hDeltaEndYnCap->Fill(DeltaEndY);}
    if(g4PrimaryProcess[0] == 5){hDeltaEndYnuclearCapatureAtRest->Fill(DeltaEndY);}
    if(g4PrimaryProcess[0] == 6){hDeltaEndYDecay->Fill(DeltaEndY);}
    if(g4PrimaryProcess[0] == 7){hDeltaEndYKaonZeroInElastic->Fill(DeltaEndY);}
    if(g4PrimaryProcess[0] == 8){hDeltaEndYCoulombScat->Fill(DeltaEndY);}
    if(g4PrimaryProcess[0] == 9){hDeltaEndYMuMinusCapture->Fill(DeltaEndY);}
    if(g4PrimaryProcess[0] == 10){hDeltaEndYProtonInelastic->Fill(DeltaEndY);}
    if(g4PrimaryProcess[0] == 11){hDeltaEndYKaonPlusInelastic->Fill(DeltaEndY);}
   
    if(g4PrimaryProcess[0] == 1){hDeltaEndXInElastic->Fill(DeltaEndX);}
    if(g4PrimaryProcess[0] == 2){hDeltaEndXNeutronInElastic->Fill(DeltaEndX);}
    if(g4PrimaryProcess[0] == 3){hDeltaEndXHadElastic->Fill(DeltaEndX);}
    if(g4PrimaryProcess[0] == 4){hDeltaEndXnCap->Fill(DeltaEndX);}
    if(g4PrimaryProcess[0] == 5){hDeltaEndXnuclearCapatureAtRest->Fill(DeltaEndX);}
    if(g4PrimaryProcess[0] == 6){hDeltaEndXDecay->Fill(DeltaEndX);}
    if(g4PrimaryProcess[0] == 7){hDeltaEndXKaonZeroInElastic->Fill(DeltaEndX);}
    if(g4PrimaryProcess[0] == 8){hDeltaEndXCoulombScat->Fill(DeltaEndX);}
    if(g4PrimaryProcess[0] == 9){hDeltaEndXMuMinusCapture->Fill(DeltaEndX);}
    if(g4PrimaryProcess[0] == 10){hDeltaEndXProtonInelastic->Fill(DeltaEndX);}
    if(g4PrimaryProcess[0] == 11){hDeltaEndXKaonPlusInelastic->Fill(DeltaEndX);}

    // ###########################################################
    // ### Looping over the spacepoints to fill the histograms ###
    // ###########################################################
    for(size_t npoints = 0; npoints < nKaonSpts; npoints++) {
      // ### Filling the incidient histogram ###

      //if(kineticEnergy < 0) { continue; }

      hdataKaonIncidentKE->Fill(kineticEnergy, EventWeight);
	
      dEdxvP->Fill(Kaondedx[npoints], kineticEnergy);

      // ### Filling the interaction histogram for the last spt ###
      if(npoints == nKaonSpts -1 && !ExitingTrack && !Decay) {
	// (trkpida[RecoTrackIndex][0] < 13.0 || (g4PrimaryProcess[0] != 6 && trkpida[RecoTrackIndex][0] >= 13.0))) {
	fKaonInteractions->Fill(kineticEnergy, EventWeight); // AH haha Alright
	fEndKE->Fill(kineticEnergy);
      }

      //float energyLossInStep = Kaondedx[npoints] * Kaonresrange[npoints] * RecombinationFactor;
      float energyLossInStep = Kaondedx[npoints] * Kaonpitchhit[npoints];
      kineticEnergy -= energyLossInStep;

      /*
      bool fixed = false;
      if(energyLossInStep > 50) { 
	if(npoints != nKaonSpts-1 && npoints != 0) { 
	  if(Kaondedx[npoints+1] * Kaonpitchhit[npoints+1] < 50 && 
	     Kaondedx[npoints-1] * Kaonpitchhit[npoints-1] < 50) {
	    printf(" Before : %lf ", energyLossInStep);
	    energyLossInStep = (Kaondedx[npoints-1] * Kaonpitchhit[npoints-1] + 
				Kaondedx[npoints+1] * Kaonpitchhit[npoints+1]) / 2;
	    printf(" After : %lf \n", energyLossInStep);
	    fixed = true;
	  } else {
	    
	    for(size_t i = npoints; i < nKaonSpts; i++) {
	      printf(" dEdx %lf \n", Kaondedx[i] * Kaonpitchhit[i]);
	    }
	    printf(" \n "); 
	  }
	} 
      } else { fixed = true; }

      if(fixed) { kineticEnergy -= energyLossInStep; }
      else { kineticEnergy -= 40; }
*/
    } // End npoints loop
  } // End of Loop();

  // ###################################################################
  // #### Looping over the exiting bins to extract the cross-section ###
  // ###################################################################
  for( int iBin = 1; iBin <= fKaonInteractions->GetNbinsX(); ++iBin ) {
	
    // ### If an incident bin is equal to zero then skip that bin ###
    if( hdataKaonIncidentKE->GetBinContent(iBin) == 0 ) continue; // Temporary fix to ensure that no Infinities are propagated to pad
	
    // ### Cross-section = (Exit Bins / Incident Bins) * (1/Density) * (1/slab width) ###
    float TempCrossSection = (fKaonInteractions->GetBinContent(iBin)/hdataKaonIncidentKE->GetBinContent(iBin)) * (1/number_density) * (1/slab_width);  
	
    float crossSection = TempCrossSection*(1/1e-28); //To put this into barns
   
    fCrossSection->SetBinContent(iBin,crossSection);
   
    float numError = pow(fKaonInteractions->GetBinContent(iBin),0.5);
    float num = fKaonInteractions->GetBinContent(iBin);
   
    if(num == 0) { num = 1; }
    float term1 = numError/num;
   
    float denomError = pow(hdataKaonIncidentKE->GetBinContent(iBin),0.5);
    float denom = hdataKaonIncidentKE->GetBinContent(iBin);
    if(denom == 0){denom = 1;}
    float term2 = denomError/denom;

    float totalError = (TempCrossSection) * (pow( ( (term1*term1) + (term2*term2) ),0.5)) * (1/number_density) * (1/slab_width) * (1/1e-28) *(1e26);

    fCrossSection->SetBinError(iBin,totalError);
       
  }//<---End iBin loop	 


  printf("\n"); 
  printf("nTotalEvents = %i, nEvtsGoodMC = %i, nEvtsTrackZPos = %i, \n", nTotalEvents, nEvtsGoodMC, nEvtsTrackZPos);
  printf("nEvtsMCTrackMatch = %i, nEventsPassingAlpha = %i, nLowZTrkEvents = %i \n", nEvtsMCTrackMatch, nEventsPassingAlpha, nLowZTrkEvents);

  printf("\n\n agree %i tot %i ratio %lf \n\n", agree, tot, (double)agree/(double)tot);


  //PIDA->Draw();
   
  fCrossSection->Draw();
  //fDiffDist->Draw();

  TFile myfile("histos_mcreco.root", "RECREATE");
  fDiffDist->Write();
  hdataUpstreamZPos->Write();
  hdataNTracksvsZpos->Write();

  // ### MC Info ###
  hMCPrimaryStartX->Write();
  hMCPrimaryStartY->Write();
  hMCPrimaryStartZ->Write();
  hMCPrimaryProjectedStartX->Write();
  hMCPrimaryProjectedStartY->Write();
  hMCPrimaryProjectedStartZ->Write();
  hMCPrimaryEndX->Write();
  hMCPrimaryEndY->Write();
  hMCPrimaryEndZ->Write();
  hMCPrimaryEndXvsZ->Write();
  hMCPrimaryEndYvsZ->Write();
  hMCPrimaryProcess->Write();
  hMCPrimaryPx->Write();
  hMCPrimaryPy->Write();
  hMCPrimaryPz->Write();
  hAlpha->Write();
  hDeltaX->Write();
  hDeltaY->Write();
  hDeltaZ->Write();  
  hMCInitalKE->Write(); 
  hdataKaondEdX->Write();
  hdataKaonRR->Write();
  hdataKaonTrkPitch->Write();
  hdataKaonIncidentKE->Write();
  fKaonInteractions->Write();
  fEndKE->Write();
  fCrossSection->Write();
  fMagPCross->Write();

  PIDA->Write();
  PIDA2->Write();

  MPV->Write();
  MPVvP->Write();
  dEdxvP->Write();

  hdataKaondEdXvsRR->Write();
  hDeltaEndX->Write();
  hDeltaEndY->Write();
  hDeltaEndZ->Write();

  hTrueLength->Write();
  hRecoLength->Write();

  hDeltaEndZInElastic->Write();
  hDeltaEndZNeutronInElastic->Write();
  hDeltaEndZHadElastic->Write();
  hDeltaEndZnCap->Write();
  hDeltaEndZnuclearCapatureAtRest->Write();
  hDeltaEndZDecay->Write();
  hDeltaEndZKaonZeroInElastic->Write();
  hDeltaEndZCoulombScat->Write();
  hDeltaEndZMuMinusCapture->Write();
  hDeltaEndZProtonInelastic->Write();
  hDeltaEndZKaonPlusInelastic->Write();

  hDeltaEndYInElastic->Write();
  hDeltaEndYNeutronInElastic->Write();
  hDeltaEndYHadElastic->Write();
  hDeltaEndYnCap->Write();
  hDeltaEndYnuclearCapatureAtRest->Write();
  hDeltaEndYDecay->Write();
  hDeltaEndYKaonZeroInElastic->Write();
  hDeltaEndYCoulombScat->Write();
  hDeltaEndYMuMinusCapture->Write();
  hDeltaEndYProtonInelastic->Write();
  hDeltaEndYKaonPlusInelastic->Write();

  hDeltaEndXInElastic->Write();
  hDeltaEndXNeutronInElastic->Write();
  hDeltaEndXHadElastic->Write();
  hDeltaEndXnCap->Write();
  hDeltaEndXnuclearCapatureAtRest->Write();
  hDeltaEndXDecay->Write();
  hDeltaEndXKaonZeroInElastic->Write();
  hDeltaEndXCoulombScat->Write();
  hDeltaEndXMuMinusCapture->Write();
  hDeltaEndXProtonInelastic->Write();
  hDeltaEndXKaonPlusInelastic->Write();

  hMCELossUpstream->Write();

  hdataKaonTrackEndX->Write();
  hdataKaonTrackEndY->Write();
  hdataKaonTrackEndZ->Write();
  hdataKaonTrackStartX->Write();
  hdataKaonTrackStartY->Write();
  hdataKaonTrackStartZ->Write();
   
}
