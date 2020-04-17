void readHistos(){

  TFile *f1=TFile::Open("neutral_INEL_ppPerugia2011.root");
  TH1D *hLtoK0=(TH1D *)f1->Get("hLambdaToK0s_Sphericity_0");
  hLtoK0->Draw();


}
