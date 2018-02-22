#include <TH1D.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TCut.h>
#include <TLatex.h>
#include "ncollFunction.h"

// for PYTHIA+HYDJET, use the following pthat weights: PTHATWEIGHTS=1., 0.00555693, 0.000213668, .0000107919
//for PYTHIA, use the following pthat weights: PTHATWEIGHTS=1.,0.0130479,0.00032629,0.000049913

// loose jet ID 
// (jtPfNHF<0.99 && jtPfNEF<0.99 && (jtPfCHM+jtPfNHM+jtPfCEM+jtPfNEM+jtPfMUM)>1) && ((abs(jteta)<=2.4 && jtPfCHF>0 && jtPfCHM>0 && jtPfCEF<0.99) || abs(jteta)>2.4) && abs(jteta)<=2.7

// TAA value
// 0-10%: 23.2210-30%: 11.51 30-50%: 3.81 50-90%: .5435
// pthat 30 0.03378 mb
// 

// Scale facor: sigma_hard = TAA * centFraction * 7745 * sigma_hard
// 
// centBin=1 0-10%  23.221*0.1*7745*0.03378 / (7745*0.1) = 0.784
// centBin=2 10-30% 11.51*0.2*7745*0.03378 / (7745*0.2) = 0.389
// centBin=3 30-50% 3.81*0.2*7745*0.03378 / (7745*0.2) = 0.116
// centBin=4 50-90% 0.543*0.4*7745*0.03378 / (7745*0.4) = 0.0183
// centBin=5 0-100% 208*208*0.3378 / (7745) = 0.189

class ForestTrees
{
   public:
   TFile *inf;
   TTree *Jet;
   TTree *match;
   TTree *Evt;
   ForestTrees(const char *infname, int cone, bool isMB=0){
      inf = new TFile(infname);
      char *treeName = Form("akCs%dPU3PFFlowJetAnalyzer/t",cone);
      Jet = (TTree*) inf->Get(treeName);
      Evt = (TTree*) inf->Get("hiEvtAnalyzer/HiTree");
      Jet->AddFriend(Evt);
      if (isMB==1) {
          TFile *inf2 = new TFile(Form("/data/yjlee/jetRAA2018/outfile_%d.root",cone));
	  TTree *matchingTree = (TTree*) inf2->Get(Form("t%d",cone));
	  Jet->AddFriend(matchingTree);
      }
   }
};

TH1D *myHistogram(const char *hname, int centBin=1)
{
  TH1D *h;
  
  if (centBin==3 || centBin==4) {
     const int nbins = 6;
     double bin[nbins+1] = {160,180,200,250,300,400,1000};
     h = new TH1D(hname,"",nbins,bin);
  } else {
     const int nbins = 7;
     double bin[nbins+1] = {160,180,200,250,300,400,500,1000};
     h = new TH1D(hname,"",nbins,bin);
  }
  h->Sumw2();

  return h;
}

void plotFakeRate(int cone=10, int centBin=1)
{
   ForestTrees fMB("LargeConeRAA_20180218/HiForestAOD_PbPb_MB_20180218_ExcludeTop4_ExcludeToFrac_Frac0p7_Full_5Sigma_MERGED.root",cone,1);
   ForestTrees f30("LargeConeRAA_20180218/HiForestAOD_PbPb_MCDijet30_20180218_ExcludeTop4_ExcludeToFrac_Frac0p7_Full_5Sigma_MERGED.root",cone);
   ForestTrees f80("LargeConeRAA_20180218/HiForestAOD_PbPb_MCDijet80_20180218_ExcludeTop4_ExcludeToFrac_Frac0p7_Full_5Sigma_MERGED.root",cone);
   ForestTrees f170("LargeConeRAA_20180218/HiForestAOD_PbPb_MCDijet170_20180218_ExcludeTop4_ExcludeToFrac_Frac0p7_Full_5Sigma_MERGED.root",cone);
   ForestTrees f280("LargeConeRAA_20180218/HiForestAOD_PbPb_MCDijet280_20180218_ExcludeTop4_ExcludeToFrac_Frac0p7_Full_5Sigma_MERGED.root",cone);
        
   TCanvas *c = new TCanvas("c","",600,600);
   c->SetLogy();
   TH1D *h30a = myHistogram("h30a");
   TH1D *h30  = myHistogram("h30");
   TH1D *h80  = myHistogram("h80");
   TH1D *h170 = myHistogram("h170");
   TH1D *h280 = myHistogram("h280");
   
   TCut centWeight = "1";
   int centL=0,centH=100;
   float scaleFactor=0.189;
   if (centBin==1){
      centL=0;centH=10;
      scaleFactor = 0.784;
   } else if (centBin==2) {
      centL=10;centH=30;
      scaleFactor = 0.389;
   } else if (centBin==3) {
      centL=30;centH=50;
      scaleFactor = 0.116;
   } else if (centBin==4) {
      centL=50;centH=90;
      scaleFactor = 0.0183;
   } else if (centBin==5) {
      centL=0;centH=100;
      scaleFactor = 0.189;
   }
   
   centWeight = TCut(getWeight(centL,centH).c_str());
      
   f30.Jet->Draw("jtpt>>h30a",centWeight*("abs(jteta)<2&&refpt>0"));
   f30.Jet->Draw("jtpt>>h30",centWeight*("abs(jteta)<2&&refpt>0&&pthat<80"));
   f80.Jet->Draw("jtpt>>h80",centWeight*("0.00555693*(abs(jteta)<2&&refpt>0&&pthat<170)"));
   f170.Jet->Draw("jtpt>>h170",centWeight*("0.0002313668*(abs(jteta)<2&&refpt>0&&pthat<280)"));
   f280.Jet->Draw("jtpt>>h280",centWeight*("0.0000107919*(abs(jteta)<2&&refpt>0)"));
   h170->SetLineColor(2);
   h280->SetLineColor(4);
   h170->SetMarkerColor(2);
   h280->SetMarkerColor(4);

   h30->SetAxisRange(1e-5,1e4,"Y");
   h30->Draw();
   h80->Draw("same");
   h170->Draw("same");
   h280->Draw("same");
   h30->Add(h80);
   h30->Add(h170);
   h30->Add(h280);
   h30->Scale(scaleFactor/f30.Jet->GetEntries(centWeight));
   
   TCut looseJetID= "(jtPfNHF<0.99 && jtPfNEF<0.99 && (jtPfCHM+jtPfNHM+jtPfCEM+jtPfNEM+jtPfMUM)>1) && ((abs(jteta)<=2.4 && jtPfCHF>0&& jtPfCHM>0 && jtPfCEF<0.99) || abs(jteta)>2.4) && abs(jteta)<=2.7";
   TH1D *hFake = myHistogram("hFake");
   TH1D *hGen = myHistogram("hGen");
   fMB.Jet->Project("hFake","jtpt",centWeight*(looseJetID&&"jtgenref==-1&&abs(jteta)<2"));
   //fMB.Jet->Project("hGen","genpt",centWeight*(looseJetID&&"abs(geneta)<2"));
   
   /*
   // don't want to double count the statistical error from GenJet subtraction
   
   for (int i=1;i<=hFake->GetNbinsX();i++){
      hFake->SetBinContent(i,hFake->GetBinContent(i)-hGen->GetBinContent(i));
   }
   */

   hFake->Scale(1./(fMB.Jet->GetEntries(centWeight)));

   
   
   h30->Draw();
   
   h30a->SetMarkerStyle(24);
   h30->SetXTitle("Jet p_{T} (GeV/c)");
   h30a->Draw("same");
   hFake->SetLineColor(2);
   hFake->SetMarkerColor(2);
   hFake->Draw("same");
 
   TCanvas *c2 = new TCanvas("c2","Fake Rate",600,600);
   c2->SetLogx();
   TH1D *hFakeRate = (TH1D*)hFake->Clone("hFakeRate");
   TH1D *hAll = (TH1D*)h30->Clone("hAll");
   hAll->Add(hFake);
   hFakeRate->Divide(hAll);  


   for (int i=1;i<=hFakeRate->GetNbinsX();i++){
      hFakeRate->SetBinContent(i,hFakeRate->GetBinContent(i));
      if (hFakeRate->GetBinContent(i)==0) {
         hFakeRate->SetBinError(i,0);}
	 
   }

   hFakeRate->SetAxisRange(-0.01,0.5,"Y");
   hFakeRate->Draw();
   hFakeRate->Draw("hist same");
   hFakeRate->SetXTitle(Form("R=%.1f Jet p_{T} (GeV/c)",cone*0.1));
   hFakeRate->SetYTitle("Fake Rate");
   hFakeRate->SetTitle(Form("%d-%d%%",centL,centH));

   TLatex* latex4 = new TLatex(180,-0.048,"2#times 10^{2}");
   latex4->SetTextFont(42);
   latex4->SetTextSize(0.045);
   latex4->Draw();
   c2->SaveAs(Form("results/fakeRate_%d_bin%d.pdf",cone,centBin));
   c2->SaveAs(Form("results/fakeRate_%d_bin%d.root",cone,centBin));
}

void doAll()
{
   for (int i=1;i<=4;i++) {
      plotFakeRate(3,i);
      plotFakeRate(4,i);
      plotFakeRate(6,i);
      plotFakeRate(8,i);
      plotFakeRate(10,i);
   }
}
