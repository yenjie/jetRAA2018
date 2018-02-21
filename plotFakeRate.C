#include <TH1D.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TCut.h>

// for PYTHIA+HYDJET, use the following pthat weights: PTHATWEIGHTS=1., 0.00555693, 0.000213668, .0000107919
//for PYTHIA, use the following pthat weights: PTHATWEIGHTS=1.,0.0130479,0.00032629,0.000049913

// loose jet ID 
// (jtPfNHF<0.99 && jtPfNEF<0.99 && (jtPfCHM+jtPfNHM+jtPfCEM+jtPfNEM+jtPfMUM)>1) && ((abs(jteta)<=2.4 && jtPfCHF>0 && jtPfCHM>0 && jtPfCEF<0.99) || abs(jteta)>2.4) && abs(jteta)<=2.7

// TAA value
// 0-10%: 23.2210-30%: 11.5130-50%: 3.8150-90%: .5435

class ForestTrees
{
   public:
   TFile *inf;
   TTree *Jet;
   TTree *Evt;
   ForestTrees(const char *infname, int cone){
      inf = new TFile(infname);
      char *treeName = Form("akCs%dPU3PFFlowJetAnalyzer/t",cone);
      Jet = (TTree*) inf->Get(treeName);
      Evt = (TTree*) inf->Get("hiEvtAnalyzer/HiTree");
      Jet->AddFriend(Evt);
   }
};

TH1D *myHistogram(const char *hname){

  const int nbins = 7;
  double bin[nbins+1] = {100,120,140,160,200,300,500,1000}
  TH1D *h = new TH1D(hname,nbins,bin);
  h->Sumw2();
  return h;
}

void plotFakeRate(int cone=10, int centBin=1)
{
   ForestTrees fMB("LargeConeRAA_20180218/HiForestAOD_PbPb_MB_20180218_ExcludeTop4_ExcludeToFrac_Frac0p7_Full_5Sigma_MERGED.root",cone);
   ForestTrees f30("LargeConeRAA_20180218/HiForestAOD_PbPb_MCDijet30_20180218_ExcludeTop4_ExcludeToFrac_Frac0p7_Full_5Sigma_MERGED.root",cone);
   ForestTrees f80("LargeConeRAA_20180218/HiForestAOD_PbPb_MCDijet80_20180218_ExcludeTop4_ExcludeToFrac_Frac0p7_Full_5Sigma_MERGED.root",cone);
   ForestTrees f170("LargeConeRAA_20180218/HiForestAOD_PbPb_MCDijet170_20180218_ExcludeTop4_ExcludeToFrac_Frac0p7_Full_5Sigma_MERGED.root",cone);
   ForestTrees f280("LargeConeRAA_20180218/HiForestAOD_PbPb_MCDijet280_20180218_ExcludeTop4_ExcludeToFrac_Frac0p7_Full_5Sigma_MERGED.root",cone);

     
   TCanvas *c = new TCanvas("c","",600,600);
   c->SetLogy();
   TH1D *h30a = myHistogram("h30a");
   TH1D *h30  = myHistogram("h30");
   f30.Jet->Draw("jtpt>>h30a","abs(jteta)<2&&refpt>0");
   f30.Jet->Draw("jtpt>>h30","abs(jteta)<2&&refpt>0&&pthat<80");
   TH1D *h80  = myHistogram("h80");
   f80.Jet->Draw("jtpt>>h80","0.00555693*(abs(jteta)<2&&refpt>0&&pthat<170)");
   TH1D *h170 = myHistogram("h170");
   f170.Jet->Draw("jtpt>>h170","0.0002313668*(abs(jteta)<2&&refpt>0&&pthat<280)");
   TH1D *h280 = myHistogram("h280");
   f280.Jet->Draw("jtpt>>h280","0.0000107919*(abs(jteta)<2&&refpt>0)");
   h170->SetLineColor(2);
   h280->SetLineColor(4);
   h170->SetMarkerColor(2);
   h280->SetMarkerColor(4);
   h30->Draw();
   h30->SetAxisRange(1e-5,1e4,"Y");
   h80->Draw("same");
   h170->Draw("same");
   h280->Draw("same");
   h30->Add(h80);
   h30->Add(h170);
   h30->Add(h280);
   h30->Scale(1./5.3/f30.Jet->GetEntries());
   
   TCut looseJetID= "(jtPfNHF<0.99 && jtPfNEF<0.99 && (jtPfCHM+jtPfNHM+jtPfCEM+jtPfNEM+jtPfMUM)>1) && ((abs(jteta)<=2.4 && jtPfCHF>0&& jtPfCHM>0 && jtPfCEF<0.99) || abs(jteta)>2.4) && abs(jteta)<=2.7";
   TH1D *hFake = new TH1D("hFake","",45,100,1000);
   hFake->Sumw2();
   fMB.Jet->Project("hFake","jtpt",looseJetID&&"abs(jteta)<2&&refpt<0");
//   f30.Jet->Project("hFake","jtpt",looseJetID&&"abs(jteta)<2&&refpt<0");
//   f80.Jet->Project("hFake","jtpt",looseJetID&&"abs(jteta)<2&&refpt<0");
//   f170.Jet->Project("hFake","jtpt",looseJetID&&"abs(jteta)<2&&refpt<0");
//   f280.Jet->Project("hFake","jtpt",looseJetID&&"abs(jteta)<2&&refpt<0");
//   hFake->Scale(1./(f30.Jet->GetEntries()+f80.Jet->GetEntries()+f170.Jet->GetEntries()+f280.Jet->GetEntries()));
   hFake->Scale(1./(fMB.Jet->GetEntries()));
   
   
   h30->Draw();
//   hFake->Scale(70./3.378e-2);
   
   h30a->SetMarkerStyle(24);
   h30->SetXTitle("Jet p_{T} (GeV/c)");
   h30a->Draw("same");
   hFake->SetLineColor(2);
   hFake->SetMarkerColor(2);
   hFake->Draw("same");
 
   TCanvas *c2 = new TCanvas("c2","Fake Rate",600,600);
   TH1D *hFakeRate = (TH1D*)hFake->Clone("hFakeRate");
   TH1D *hAll = (TH1D*)h30->Clone("hAll");
   hAll->Add(hFake);
   hFakeRate->Divide(hAll);  
   hFakeRate->SetAxisRange(0,2,"Y");
   hFakeRate->Draw();
}
