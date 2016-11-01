#define BasicChecks_cxx
#include "BasicChecks.h"
#include <iostream>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>

// Beamline studies: WC, TOF, mass
TH1D *hWCMomentum   = new TH1D("hWCMomentum"   , "All WC Momentum"   , 180, 0, 1800);
TH1D *hTOF          = new TH1D("hTOF"          , "All TOF"           , 180, 15, 60);
TH1D *hBeamLineMass = new TH1D("hBeamLineMass" , "Mass From Beamline", 200, 0, 2000);
TH2D *hWCMomentumTotVsTOF = new TH2D("hWCMomentumTotVsTOF"    , "WC Momentum vs TOF"  , 180, 0, 1800,  180, 15, 60);
TH1D *hWCXFF = new TH1D("hWCXFF","X position of WC track on TPC Front Face",  620 , 100.,  500.);
TH1D *hWCYFF = new TH1D("hWCYFF","Y position of WC track on TPC Front Face",  600 , -300., 300.);

//Calorimetry Study
TH1D *hCollectionPitch    = new TH1D("hCollectionPitch"   , "Pitch On Collection"  , 140, 0.3, 0.9);
TH1D *hInductionPitch     = new TH1D("hInductionPitch"    , "Pitch On Induction "  , 140, 0.3, 0.9);
TH1D *hCollectionDEDX     = new TH1D("hCollectionDEDX"    , "DEDX On Collection"   , 100, 0, 10);
TH1D *hInductionDEDX      = new TH1D("hInductionDEDX"     , "DEDX On Induction "   , 100, 0, 10);
TH1D *hCollectionDXDEDX   = new TH1D("hCollectionDXDEDX"  , "DXDEDX On Collection" , 100, 0, 10);
TH1D *hInductionDXDEDX    = new TH1D("hInductionDXDEDX"   , "DXDEDX On Induction " , 100, 0, 10);
TH1D *hCollectionPIDA     = new TH1D("hCollectionPIDA"    , "PIDA On Collection "  , 200, 0, 50);
TH1D *hInductionPIDA      = new TH1D("hInductionPIDA"     , "PIDA On Induction "   , 200, 0, 50);

//Occupancy Study
TH1D *hNTPCTracksFirst30cm = new TH1D("hNTPCTracksFirst30cm" , "N TPC Tracks in First 30 cm"  , 30, -0.5, 29.5);
TH1D *hNTPCTracksFirst2cm  = new TH1D("hNTPCTracksFirst2cm"  , "N TPC Tracks in First  2 cm"  , 30, -0.5, 29.5);

//Position Study
TH1D *hXStart = new TH1D("hXStart","X position track start",  620 , -10., 52.);
TH1D *hYStart = new TH1D("hYStart","Y position track start",  600 , -30., 30.);
TH1D *hZStart = new TH1D("hZStart","Z position track start", 1100 , -10., 100.);

TH1D *hXEnd = new TH1D("hXEnd","X position track start",  620 , -10., 52.);
TH1D *hYEnd = new TH1D("hYEnd","Y position track start",  600 , -30., 30.);
TH1D *hZEnd = new TH1D("hZEnd","Z position track start", 1100 , -10., 100.);

TH1D *hLinearLenght = new TH1D("hLinearLenght","Linear Lenght", 1100 , -10., 100.);



void BasicChecks::Loop()
{
//   In a ROOT session, you can do:
//      root> .L BasicChecks.C
//      root> BasicChecks t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
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
   std::cout<<"NEntries = "<<nentries<<"\n";
   //nentries = 100;
   std::cout<<"But we'll look at = "<<nentries<<"\n";
   
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      

      if (jentry % 1000 == 0) std::cout<<"Event number  "<<jentry<<std::endl;

      //Fill histograms
      for (int iWC = 0; iWC < nwctrks; iWC++) 
	{
	  //std::cout<<wctrk_XFaceCoor[iWC]<<"\n";
	  hWCMomentum  ->Fill(wctrk_momentum[iWC]);     
	  hWCXFF  ->Fill(wctrk_XFaceCoor[iWC]);     
	  hWCYFF  ->Fill(wctrk_YFaceCoor[iWC]);     
	}
      for (int iTOF = 0; iTOF < ntof; iTOF++) hTOF         ->Fill(tofObject[iTOF]);     
      //      hBeamLineMass       ->Fill();     
      hWCMomentumTotVsTOF ->Fill(wctrk_momentum[0],tofObject[0]);          

      
      int TPCTracksBefore30cm = 0;
      int TPCTracksBefore2cm  = 0;
      
      // Loop on tracks
      for(int iTPCtrk = 0; iTPCtrk < ntracks_reco; iTPCtrk++)
        {
	  // Check which direction the track was reconstructed
          bool   invertedOrder     = false;
          double correctFirstXPoint = trkvtxx[iTPCtrk];
          double correctFirstYPoint = trkvtxy[iTPCtrk];
          double correctFirstZPoint = trkvtxz[iTPCtrk];
	  double correctLastXPoint  = trkendx[iTPCtrk];
          double correctLastYPoint  = trkendy[iTPCtrk];
          double correctLastZPoint  = trkendz[iTPCtrk];

	  
	  double linearLenght= TMath::Sqrt((correctFirstXPoint-correctLastXPoint)*(correctFirstXPoint-correctLastXPoint) +
					   (correctFirstYPoint-correctLastYPoint)*(correctFirstYPoint-correctLastYPoint) +
					   (correctFirstZPoint-correctLastZPoint)*(correctFirstZPoint-correctLastZPoint)
					   );
	  hLinearLenght ->Fill(linearLenght); 
	  
          if ( trkvtxz[iTPCtrk] > trkendz[iTPCtrk] )
            {
              //std::cout<<"Wrong direction, baby! "<<std::endl;
              invertedOrder  = true;
	      correctFirstXPoint = trkendx[iTPCtrk];
	      correctFirstYPoint = trkendy[iTPCtrk];
	      correctFirstZPoint = trkendz[iTPCtrk];
	      correctLastXPoint  = trkvtxx[iTPCtrk];
	      correctLastYPoint  = trkvtxy[iTPCtrk];
	      correctLastZPoint  = trkvtxz[iTPCtrk];
            }
                        
	  //Loop on Collection Plane Hits
	  for (int iHitC = 0; iHitC < trkhits[iTPCtrk][0]; iHitC++ )
	    {
	      hCollectionPitch->Fill(trkpitchhit[iTPCtrk][0][iHitC]);
	      hCollectionDEDX ->Fill(trkdedx[iTPCtrk][0][iHitC]);
	      hCollectionDXDEDX->Fill(trkdedx[iTPCtrk][0][iHitC]*trkpitchhit[iTPCtrk][0][iHitC]);
	    }//End of loop on Collection Plane Hits
	  
          
	  //Loop on Induction Plane Hits
	  for (int iHitI = 0; iHitI < trkhits[iTPCtrk][1]; iHitI++ )
	    {
	      hInductionPitch ->Fill(trkpitchhit[iTPCtrk][1][iHitI]);
	      hInductionDEDX  ->Fill(trkdedx[iTPCtrk][1][iHitI]);
	      hInductionDXDEDX->Fill(trkdedx[iTPCtrk][1][iHitI]*trkpitchhit[iTPCtrk][1][iHitI]);
	    }//End of loop on Induction Plane Hits
	  
	  hXStart ->Fill(correctFirstXPoint);
	  hYStart ->Fill(correctFirstYPoint);
	  hZStart ->Fill(correctFirstZPoint);
	  hXEnd   ->Fill(correctLastXPoint ); 
	  hYEnd   ->Fill(correctLastYPoint ); 
	  hZEnd   ->Fill(correctLastZPoint ); 
	  
	  if (correctFirstZPoint < 30) TPCTracksBefore30cm++;
	  if (correctFirstZPoint <  2)  TPCTracksBefore2cm++;

	  hCollectionPIDA     ->Fill(trkpida[iTPCtrk][0]); 
	  hInductionPIDA      ->Fill(trkpida[iTPCtrk][1]); 
	}
      
      hNTPCTracksFirst30cm ->Fill(TPCTracksBefore30cm);
      hNTPCTracksFirst2cm  ->Fill(TPCTracksBefore2cm);
   }

 
  
   TFile myfile("BasicStudies.root","RECREATE"); 
   hWCMomentum         ->Write();     
   hTOF                ->Write();     
   hBeamLineMass       ->Write();     
   hWCMomentumTotVsTOF ->Write();     
   hWCXFF  ->Write();     
   hWCYFF  ->Write();     
   
   hCollectionPitch    ->Write(); 
   hInductionPitch     ->Write(); 
   hCollectionDEDX     ->Write(); 
   hInductionDEDX      ->Write(); 
   hCollectionDXDEDX   ->Write(); 
   hInductionDXDEDX    ->Write(); 
   hCollectionPIDA     ->Write(); 
   hInductionPIDA      ->Write(); 
   
   hNTPCTracksFirst30cm ->Write();
   hNTPCTracksFirst2cm  ->Write();
   
   hXStart ->Write();
   hYStart ->Write();
   hZStart ->Write();
   
   hXEnd ->Write(); 
   hYEnd ->Write(); 
   hZEnd ->Write(); 
   hLinearLenght->Write(); 
}
