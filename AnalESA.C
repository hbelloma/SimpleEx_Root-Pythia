/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//*************** Calculador de Event Shape Analysis  ******************** //
//                                                                         //
//  Este programa hace el analisis de la forma de eventos de jets tales    //
//  como el thrust, esfericidad, thrustminor, recoil, de eventos montecarlo//
//  ya almacenados en galice.root  y su respectivo kinematics, asi mismo   //
//  realiza las graficas respectivas.                                      //
//quejas y/o sugerencias con el                                            //
//Author: Héctor Bello Martínez                                            //
//FCFM-BUAP   for ALICE collaboration                                      //
/////////////////////////////////////////////////////////////////////////////

#include <TProfile.h>
#include<math.h>
#include "TMath.h"

#include<stdio.h>
#include<stdlib.h>
#include<float.h>
#include "TArrayD.h"
#include "TMatrixD.h"
#include "TMatrixDBase.h"
#include "TInterpreter.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TFile.h"
#include "TH2.h"
#include "TNtuple.h"
#include "TPaveLabel.h"
#include "TPaveText.h"
#include "TFrame.h"


const Int_t nbines=16;
const char* dirpart = "/home/hector/tareasALICE";
// "/storage/alice/aortizve/simulaciones7TeVMayo2010/PERUGIA0";

TParticle* leading(TObjArray* list);
void AnalESA()
{ 
  gSystem->Load("libTree.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libGeom.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libSTEERBase.so");
  gSystem->Load("libAOD.so");
  gSystem->Load("libESD.so");
  gSystem->Load("libANALYSISalice.so");
  gSystem->Load("libJETAN.so");

  fn = new TH1F("fn", "", 13,-1,12);

  TProfile *fNchTaug[3];
  TProfile *fNchStg[3];
  TProfile *fNchRg[3];
  TProfile *fNchTmg[3];

  for (Int_t i=0; i<3; ++i)
    {
      fNchTaug[i] = new TProfile(Form("NchTau_g_%d", i), Form("NchTau_g_%d", i), 20,0,63,0,1000);
      fNchTaug[i]->GetXaxis()->SetTitle("Multiplicity");
      fNchTaug[i]->GetYaxis()->SetTitle("Thrust (tau)");
      fNchStg[i] = new TProfile(Form("NchSt_g_%d", i), Form("NchSt_g_%d", i), 20,0,63,0,1000);
      fNchStg[i]->GetXaxis()->SetTitle("Multiplicity");
      fNchStg[i]->GetYaxis()->SetTitle("Trans Sphericity");
      fNchRg[i] = new TProfile(Form("NchR_g_%d", i), Form("NchR_g_%d", i), 20,0,63,0,1000);
      fNchRg[i]->GetXaxis()->SetTitle("Multiplicity");
      fNchRg[i]->GetYaxis()->SetTitle("Recoil");
      fNchTmg[i] = new TProfile(Form("NchTm_g_%d", i), Form("NchTm_g_%d", i), 20,0,63,0,1000);   
      fNchTmg[i]->GetXaxis()->SetTitle("Multiplicity");
      fNchTmg[i]->GetYaxis()->SetTitle("Thrust minor");   
    }

  Int_t finalBB=2000;  //2000; //100000;
  Int_t numevents=2000;//2000;
  for(Int_t run=2000; run <finalBB+1; ++run){
    char path[256];
    sprintf(path, "%s/galice.root",dirpart,run);//"%s/%d/galice.root", dirpart,run);
    cout << path << endl;
    AliRunLoader* rl = AliRunLoader::Open(path);
    rl->LoadKinematics();
    rl->LoadHeader(); 
    TFile *outesaA = new TFile("ESAHpt.root","recreate");
    for (Int_t ie=0; ie<numevents; ie++) {  	//mitad de finalBB porque?
      rl->GetEvent(ie);
      AliStack* stack = rl->Stack();
      ///////////////////////////////////////Bulk analysis//////////////////////////////////////////
      Int_t NchCut=30;
      Int_t nPrim  = stack->GetNprimary();
      Double_t etacutoffa=0.8;
      Double_t ptcutoffa=0.5;
      Double_t ptcutoffhard=0.5;
      Double_t etacutoffhard=0.8;
      Bool_t neutral=kFALSE;
      Int_t nana=3;
      ////////////////////////////////////////////////////////////////////////////////
      //AliStack* stack = 0;
      TObjArray* acceptedtracks = new TObjArray(); 
      Int_t  nmctracks2=0; //inic numero de tracks que pasan Trigger
      Int_t  nmctracks1=0; //inic numero de tracks primarios neutros o cargados
      for (Int_t iMCTracks = 0; iMCTracks < nPrim; iMCTracks++) {    
           TParticle* trackmc = stack->Particle(iMCTracks);
           if (!trackmc) continue;
           Double_t etamc =trackmc ->Eta();
           Double_t ptmc=trackmc->Pt();
           Int_t pdgCode = TMath::Abs(trackmc->GetPdgCode());
           if (TMath::Abs(etamc) > etacutoffa) continue; //only particles in |eta|<=etacutoff
           if(ptmc < ptcutoffa) continue;  // PT cut
           Bool_t isprimary = stack->IsPhysicalPrimary(iMCTracks);    // Check if particle is charged, and primary
           if(isprimary == 0) continue;  // only primary particles
           TParticlePDG* pdgPart =trackmc ->GetPDG();
	   if(neutral == 1){//include neutral particles
	        // skip photons and neutrinos
	        if (pdgCode == 22 || pdgCode == 12 || pdgCode == 14 || pdgCode == 16) continue;
           }
           else{ //only charged particles
	       if (pdgPart->Charge() == 0)continue;
	   }    
           nmctracks1++;   //numero de mctracks de part.//se puede quitar
           if (TMath::Abs(etamc) > etacutoffa) continue;     
           if(ptmc < ptcutoffa) continue;                    
	   nmctracks2++; //numero de tracks que pasan cutts   
	   acceptedtracks->Add(trackmc); //son los iMCTracks  que estan en el stack =trackmc que pasan las pruebas de cortes
      }
      //ClassImp(AnalESA)
      //TArrayD* eventshapes = 0;//inicializa eventshapes
      //eventshapes = AnalysisEvShapeH100000::GetThrustParamMC(acceptedtracks, 	, ptcutoffa, etacutoffa, neutral);
      //TArrayD * AnalESA::GetThrustParamMC(acceptedtracks, Int_t  nana, Double_t ptcutoffa, Double_t etacutoffa, Bool_t neutral) 
      //{
      //}
      float ptsuma = 0;
      Double_t pxsuma = 0;
      Double_t pysuma = 0;
      TArrayD* evsh = new TArrayD(4);
      
      
      Int_t j=0;
      pxT = new Double_t[nmctracks2];
      pyT = new Double_t[nmctracks2];
      ptT = new Double_t[nmctracks2];
      for (Int_t i = 0; i < nmctracks2; i++)
        {
          pxT[i] = 0;
          pyT[i] = 0;
          ptT[i] = 0;
        }
      for (Int_t iMCTracks = 0; iMCTracks < nPrim; ++iMCTracks) {    
          TParticle* trackmc = stack->Particle(iMCTracks);
          if (!trackmc) continue;
          Double_t etamc = trackmc ->Eta();
          Double_t pxmc  = trackmc->Px();
          Double_t pymc  = trackmc->Py();
          Double_t ptmc  = trackmc->Pt();
          Int_t pdgCode  = TMath::Abs(trackmc->GetPdgCode());
          if (TMath::Abs(etamc) > etacutoffa) continue;
          if(ptmc < ptcutoffa) continue;
          Bool_t isprimary = stack->IsPhysicalPrimary(iMCTracks); 
          if(isprimary==0) continue;
          TParticlePDG* pdgPart =trackmc ->GetPDG();

          if(neutral==1){  
              if (pdgCode == 22 || pdgCode == 12 || pdgCode == 14 || pdgCode == 16)continue;
          } else {
	      if (pdgPart->Charge() == 0) continue;
          }
     
          ptT[j] = ptmc;  
          pxT[j] = pxmc;
          pyT[j] = pymc;
          ptsuma += ptmc; //SUMA DE Pt
          //printf("\n ptsuma = %f\n ",ptsuma);
          pxsuma+=pxmc;   //SUMA DE pX
          pysuma+=pymc;   //SUMA DE Py
          j++;    
      }
      float numerador = 0;
      float numerador2 = 0;
      Double_t phimax = -1;  
      Double_t pFull = -1;
      Double_t pMax = 0;
      Double_t phi = 0;
      Double_t Thrustg = 80;  // Inicializa Thrust
      Double_t Thrustming = 80;//Inicializa Thrustmin
      Double_t nx = 0;
      Double_t ny = 0;
      Double_t phiparam = 0;
///////////////////////////////GETTING THRUST (obtiene thrust) ////////////////////////////////////////////////////
      for(Int_t i = 0; i < 360; ++i){  //360 //phi
          numerador = 0;
          phiparam  = 0;
          nx = 0;
          ny = 0;
          phiparam=((TMath::Pi()) * i) / 180; // parametrization of the angle
          nx = TMath::Cos(phiparam);            // x component of an unitary vector n
          ny = TMath::Sin(phiparam);            // y component of an unitary vector n
          for(Int_t i1 = 0; i1 < nmctracks2; ++i1){
	      numerador += TMath::Abs(nx * pxT[i1] + ny * pyT[i1]);//product between momentum proyection in XY plane and the unitari vector.
          }
          if(TMath::Abs(ptsuma)==0) 
          Double_t pFull=1;
          else  Double_t pFull=(numerador/ptsuma);
          printf("\n numerador=%f, ptsuma= %f, pFull = %f\n ",numerador, ptsuma, pFull);
          if(pFull > pMax)//maximization of pFull
          {
              pMax = pFull;
              phi = phiparam;
          } 
      }
      phimax=(phi * 180) / TMath::Pi();//angular parameter of the unitary vector which maximiza thrust
      //if n vector and beam axis form a plane, then we can calculate a second unitary vector perpendicular to that plane
      Double_t nx1 = TMath::Cos(phi);
      Double_t ny1 = TMath::Sin(phi);
      for(Int_t i2 =0; i2 < nmctracks2; ++i2){
          numerador2 += TMath::Abs(pxT[i2] * ny1 - nx1 * pyT[i2]);//NORM OF cross product: P_{i} X n, P_{i}=(px_{i},py_{i})
      }
      Thrustg = 1 - pMax;//this is the value of thrust (tau definition)
      if(TMath::Abs(ptsuma)==0){ 
      Double_t Thrustming = 0;
      Double_t recoilg = 0;}
      else {
      Thrustming=(numerador2/ptsuma);//this is the value of thrust minor
      Double_t recoilg = TMath::Abs(TMath::Sqrt(pxsuma * pxsuma + pysuma * pysuma)) / (ptsuma);} //recoil /checar
      printf("\n Evento=%f, Thrust(tau) = %f, Thrustmin =%f, recoil=%f\n ", ie, Thrustg,Thrustming, recoilg); 
///////////////////////////  GETTING SPHERICITY (obtiene esfericidad transversa linearizada) //////////////////////////////////
/////////////////////////// Adaptado del programa GeneraSph2.C del mismo autor///////////////////////////////////////////////////////////////
      float spx2=0; 
      float spxpy=0;
      float spy2=0;
      float invpt=0;
      for(Int_t i3 = 0; i3 < nmctracks2; ++i3){
              if (ptT[i3]==0) { 
                 float invpt=1; }
                 else  {
                 float invpt=(1/ptT[i3]);
                 }
              float px2=(pxT[i3]**2)*invpt;  //esto se han modificado por altos pt y la linearizacion
              float py2=(pyT[i3]**2)*invpt;  //este también
              float pxpy=pxT[i3]*pyT[i3]*invpt;   // y finalmente; este también
              float spx2 = spx2+px2;
              float spy2 = spy2+py2;
              float spxpy=spxpy+pxpy;
      }
      // Sumas totales
      float suma[4];
      suma[0]=spx2; //suma px2_i*invspt
      suma[1]=spxpy;//suma px_ipy_i*invspt
      suma[2]=spxpy;//suma py_ipx_i*invspt
      suma[3]=spy2; //suma py2_i*invspt

      // Ahora escribimos la matriz de momentos transversos linealizada para cada partícula
      TMatrixD ptmatrix(2,2);
      TArrayD ptmatrix_ij(4);

      for (Int_t k= 0; k < 4; k++){ //aquí le decimos como van las entradas
          if(TMath::Abs(ptsuma)==0){ 
          ptmatrix_ij[k]=1;}  // la esfericidad de la matrix {{1,1},{1,1}} con egenvals=0,2 es cero no contribuirian
          else { 
          ptmatrix_ij[k]=suma[k]*(1/ptsuma);}  //
      }
      ptmatrix.SetMatrixArray(ptmatrix_ij.GetArray()); //Esto carga las entradas

      // Esto Encuentra los eigenvalores
      TMatrixDEigen eigen(ptmatrix);
      TVectorD eigenvals = eigen.GetEigenValues();
      //vector<float> eigenvaluess(2);
      float eigenvalue[2];
      eigenvalue[0] = eigenvals[0];
      eigenvalue[1] = eigenvals[1];

      // calculamos la Transversal Sphericity
      if (eigenvalue[0]>eigenvalue[1])
      float mineigenval=eigenvalue[1]; // minimo eigenvalor
      else mineigenval=eigenvalue[0];

      float sphericityg = (2*mineigenval)/(eigenvalue[0]+eigenvalue[1]);
      printf("\n sphericity=%f\n ", sphericityg);
/////////////////////////////////////////////////////////////////////////////////////
      
      evsh->AddAt(Thrustg, 0);//Thrust entrada 0 en evsh array
      evsh->AddAt(Thrustming, 1);//Thrustmin entrada 1 en evsh array 
      evsh->AddAt(recoilg, 2);//recoil entrada 2 en evsh array
      evsh->AddAt(sphericityg, 3); //sphericity entrada 3 en evsh array
      delete [] ptT;
      delete [] pxT;
      delete [] pyT;

      //return evsh;  
      //Double_t Thrustg=eventshapes->GetAt(0);
      //Double_t Thrustming=eventshapes->GetAt(1);
      //Double_t recoilg=eventshapes->GetAt(2);  
      if((Thrustg<0)&&(recoilg<0))continue;

      fn->Fill(0);
      fNchTaug[0]->Fill(nmctracks2,Thrustg);
      fNchTmg[0]->Fill(nmctracks2,Thrustming);
      fNchRg[0]->Fill(nmctracks2,recoilg);
      fNchStg[0]->Fill(nmctracks2,sphericityg);
      
      TParticle * TriggerMC=leading(acceptedtracks);
      if(TriggerMC == NULL)return;
      Double_t pttmc=TriggerMC->Pt();
      //printf("\n TriggerMC=%f\n", TriggerMC); 
      //printf("\n pttmc=%f\n", pttmc); 

      //ALL
      //Int_t caso2=0;
      //if(pttmc<2)caso2=1;//soft scattering
      //if(pttmc>=2)caso2=2;//hard scattering

      //switch(caso2)
      	//{
	//case 1:
	  //{
	    //fn->Fill(1);
	    //fNchTaug[1]->Fill(nmctracks2,Thrustg);
	    //fNchTmg[1]->Fill(nmctracks2,Thrustming);
	    //fNchRg[1]->Fill(nmctracks2,recoilg);
	    //fNchStg[1]->Fill(nmctracks2,sphericityg);
	  //}
	  //break;	  
	//case 2:
	  //{
	    //fn->Fill(2);
	    //fNchTaug[2]->Fill(nmctracks2,Thrustg);
	    //fNchTmg[2]->Fill(nmctracks2,Thrustming);
	    //fNchRg[2]->Fill(nmctracks2,recoilg);
	    ////fNchStg[2]->Fill(nmctracks2,sphericityg);
	  //}
	  //break;  
	//}
            
      delete acceptedtracks;



    }   //loop eventos  
    delete rl;  
  }//loop files
  
  outesaA->cd();
  fn->Write();

  for (Int_t i=0; i<3; ++i)
    {
      fNchTaug[i] -> Write();
      fNchStg[i]  -> Write();
      fNchRg[i]   -> Write();
      fNchTmg[i]  -> Write();
    }
 
  outesaA->Close();
}
/////////////// otro programita que hace el Trigger //////////////////////7
TParticle *leading(TObjArray* list)
{
  Int_t nPriml  = list->GetEntries();
  TParticle* leadingch=NULL;
  Double_t maxipt=0;
  for (Int_t iTracks = 0; iTracks < nPriml; ++iTracks) 
    {
      TParticle* track = dynamic_cast<TParticle*> (list->At(iTracks));
      Double_t pt = track->Pt();
      if(pt>maxipt){
	maxipt=pt;
	leadingch=track;
      }
      
    }
  return leadingch;
}
