//
// Using Pythia8 with ROOT
//
// Author: Héctor Bello Martínez <hbelloma@cern.ch>
// Update: 2018-10-12
//
// Modification by: 
// Date:
//
//-------------------------------------------------------------------
// To make an event sample (of size 1000) do
//
//    shell> root
//    root [0] .L py8Exam.C
//    root [1] gSystem->Load("libpythia8");
//    root [2] gSystem->Load("libqpythia.so");   
//    root [3] makeEventSample(1000)
//



#ifndef __CINT__
#include "TApplication.h"
//#include "TPythia8.h"
#include "Pythia8/Pythia.h"
#include "TFile.h"
#include "TError.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TH1.h"
#include "TF1.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "Riostream.h"
#include "TDatabasePDG.h"
#include "TDecayChannel.h"
#include "TParticlePDG.h"
#include <cstdlib>
#include "TMCParticle.h"
#include "Event.h"

using namespace std;
using namespace Pythia8;
#endif

#define FILENAME    "py8Exam.root"
#define OUTDIREC     "py8Exam"
#define OUTFILEPOSTPRO "py8Exampospro.root"
#define IMAGEF     "pdf"
#define ECM        13000
#define TREENAME   "treetest"
#define BRANCHNAME "particles"
#define INTEVT 10051


//---------------------------------------------------------------------------
void CreateDir(const Char_t* dirName)
{
   TString pwd(gSystem->pwd());
   gSystem->cd(pwd.Data());

   if(gSystem->cd(dirName)) {
   gSystem->cd(pwd.Data());
   } else {
   gSystem->mkdir(dirName, kTRUE); // kTRUE means recursive
   }
}
//______________________________
TString CompleteDirName( TString nameOut3 ){
   TString nameOut = nameOut3;
   nameOut += "histos";
   return nameOut;
}



int makeEventSample(Int_t nEvents)
{
  Bool_t isfast=kFALSE; //kTRUE;// kFALSE; 
  if(nEvents<101) isfast=kTRUE;

  Pythia pythia;
  pythia.readString("Beams:idA = 2212");
  pythia.readString("Beams:idB = 2212");
  pythia.readString("Beams:eCM = 7000");

  pythia.readString("SoftQCD:nonDiffractive = on");
  pythia.readString("SoftQCD:singleDiffractive = off");
  pythia.readString("SoftQCD:doubleDiffractive = off");
  pythia.readString("SoftQCD:centralDiffractive = off");
  pythia.readString("HardQCD:all = off");
  pythia.readString("PartonLevel:all = on");
  pythia.readString("PartonLevel:ISR = on");
  pythia.readString("PartonLevel:FSR = on");
  pythia.readString("HadronLevel:all = on");

  pythia.readString("Tune:pp = 14");
  pythia.readString("Tune:ee = 7");


  pythia.readString("130:onMode = off");
  pythia.readString("321:onMode = off");
  pythia.readString("313:onMode = off");
  pythia.readString("323:onMode = off");
  pythia.readString("311:onMode = off");
  pythia.readString("310:onMode = off");
  pythia.readString("3122:onMode = off");
  pythia.readString("3112:onMode = off");
  pythia.readString("3222:onMode = off");
  pythia.readString("3212:onMode = off");
  pythia.readString("3312:onMode = off");
  pythia.readString("3322:onMode = off");
  pythia.readString("3334:onMode = off");


  pythia.init(); 

  // Open an output file
  TFile* file = TFile::Open(FILENAME, "RECREATE");
  if (!file || !file->IsOpen()) {
    return 1;
  }

  TH1D* hNch = new TH1D("hNch","",300,0,300);
  TH1F* hetach = new TH1F("hetach", "Pseudorapidity", 240, -12., 12.);
  TH1D* hevents=new TH1D("hevents","Eventos",10,0,10);

 Float_t fy=-999999;
 Double_t ptch=-99999;
 Double_t pxch =-99999;
 Double_t pych =-99999;
 Double_t pzch =-99999;
 Double_t ech = -99999;


     // Make a tree in that file ...
  // read https://root.cern.ch/trees-five-steps
  TTree* myTree = new TTree(TREENAME, "Pythia 8 tree");
//  TTree aliTree("aliTree", "/aliroot");

//  Event *event = new  Event();
//  myTree->Branch("EventBranch", "Event", &event);
//Adding branches with a list of variables is done with the Branch() method which has the following parameters:
//   TBranch *b = tree->Branch ("Ev_Branch",&event, "ntrack/I:nseg:nvtex:flag/i:temp/F");
//Adding branches with a TClonesArray:

  // myTree->Branch( "Track_B", &Track, 64000);
   myTree->Branch( "pT", &ptch, "pT/F");


//Array for TMCParticle objects (only pythia6). 
//  TClonesArray* particles = (TClonesArray*)pythia->GetListOfParticles();
//  tree->Branch(BRANCHNAME, &particles);
  //rama de rapidity
//   Float_t fy=-999999;
   // Need to fix-----------------------ME
   // TLeaf *fyl=particles("fy",&fy);
//   particles->Leaf("fy",&fy);
 //  TLeaf(TBranch* BRANCHNAME, "fy", "f/y")
//   TBranch *newBranch = tree->Branch("fy", &fy, "fy/F");


//Listasdepy

    pythia.settings.listAll();

  // Now we make some events
  for (Int_t ievt = 0; ievt < nEvents; ievt++) {
    // Show how many events we have.
    if(isfast) cout << "Event # " << ievt << endl;
    else{if(ievt % 100 == 0) cout << "Event # " << ievt << endl;}

    // Make event.
    pythia.next();
    //----For 1 event --------
    hevents->Fill(0);

    Int_t np = pythia.event.size();
 //   cout<<"np= "<<np<<endl;
    Int_t nmultch=0;
    Double_t etach=-999.;
    Double_t ych=-999.;
    Double_t phichmv=-999.;
    Double_t phich=-999.;
    

//    cout<<"starting loop of particles np= "<<np<<endl;
    // loop over particles
    for (Int_t ip=0; ip<np; ip++){

      Int_t idp = pythia.event[ip].id();
      Int_t iMother = pythia.event[ip].mother1();
      Int_t iMotherstat = pythia.event[iMother].status();
      Int_t idMother = pythia.event[iMother].id();
      Int_t iMother2 = pythia.event[ip].mother2();
      Int_t iMother2stat = pythia.event[iMother2].status();
      Int_t idMother2 = pythia.event[iMother2].id();

      if(TMath::Abs(idp)>4000012) continue; // protection to not known pdg codes in root  
      Int_t pstat = pythia.event[ip].status();
       pxch = pythia.event[ip].px();
       pych = pythia.event[ip].py();
       pzch = pythia.event[ip].pz();
       ech = pythia.event[ip].e();
       ptch=TMath::Sqrt(pxch*pxch+pych*pych);
       fy=ptch;
       Double_t costhetach=pzch/(TMath::Sqrt(ptch*ptch+pzch*pzch));
       ych= TMath::ATan(pzch/ech);
       etach=TMath::ATanH(costhetach);
       if(ptch!=0){
                  phichmv=TMath::ACos(pxch/ptch);
                  if(pych<0) phich=2*TMath::Pi()-phichmv;
                  else phich=phichmv;
                  }

       Double_t papx= pythia.event[iMother].px();
       Double_t papy= pythia.event[iMother].py();
       Double_t papz= pythia.event[iMother].pz();
       Double_t pae = pythia.event[iMother].e();

   if(pythia.event[ip].isCharged()){
       if(isfast) cout<<"is charged pass"<<endl;
       if(ptch<0.15) continue; //vary cut to 5 GeV see dep on So //newfix
       if(isfast) cout<<"ptch cut pass"<<endl;
       if(pstat<1)continue;
       if(isfast) cout<<"pstat passed"<<endl;
       hetach->Fill(etach);
       if(TMath::Abs(etach)>0.8) continue;//0.8 //ALICE paper eta 0.9
       if(isfast) cout<<"eta cut passed"<<endl;
       nmultch+=1;


   } //end loop of charged within kin cuts
 }// end particleloop

    hNch->Fill(nmultch);
    myTree->Fill();
  }//end loop of events

   myTree->Print();
  file->Write();
  file->Close();
/*  
*/
  return 0;
}





