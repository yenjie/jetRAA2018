#include <TTree.h>
#include <TH1D.h>
#include <TFile.h>
#include <iostream>

#include "commonTool.h"

class ForestTreesMatch
{
   public:
   TFile *inf;
   TTree *Jet;
   TTree *Evt;
   ForestTreesMatch(const char *infname, int cone){
      inf = new TFile(infname);
      char *treeName = Form("akCs%dPU3PFFlowJetAnalyzer/t",cone);
      Jet = (TTree*) inf->Get(treeName);
      Evt = (TTree*) inf->Get("hiEvtAnalyzer/HiTree");
      Jet->AddFriend(Evt);
   }
};

void genMatching(Int_t cone=3)
{
   ForestTreesMatch fMB("/data/cmcginn/Forests/Pythia6Dijet/LargeConeRAA_20180218/HiForestAOD_PbPb_MB_20180218_ExcludeTop4_ExcludeToFrac_Frac0p7_Full_5Sigma_MERGED.root",cone);
      
   Int_t nref, ngen;
   Float_t jtpt[100];
   Float_t jteta[100];
   Float_t jtphi[100];
   Int_t jtgenref[100];
   Float_t genpt[100];
   Float_t genphi[100];
   Float_t geneta[100];
   Float_t jtDR[100];
   
   fMB.Jet->SetBranchAddress("nref", &nref);
   fMB.Jet->SetBranchAddress("jtpt", jtpt);
   fMB.Jet->SetBranchAddress("jteta", jteta);
   fMB.Jet->SetBranchAddress("jtphi", jtphi);
   
   fMB.Jet->SetBranchAddress("ngen", &ngen);
   fMB.Jet->SetBranchAddress("genpt", genpt);
   fMB.Jet->SetBranchAddress("geneta", geneta);
   fMB.Jet->SetBranchAddress("genphi", genphi);
   
   TFile *outf = new TFile(Form("outfile_%d.root",cone),"recreate");
   
   TTree *t = new TTree(Form("t%d",cone),"");
   
   t->Branch("nref",&nref,"nref/I");
   t->Branch("ngen",&ngen,"ngen/I");
   t->Branch("jtgenref",jtgenref,"jtgenref[nref]/I");
   t->Branch("jtpt",jtpt,"jtpt[nref]/F");
   t->Branch("jteta",jteta,"jteta[nref]/F");
   t->Branch("jtphi",jtphi,"jtphi[nref]/F");
   t->Branch("jtDR",jtDR,"jtDR[nref]/F");
   t->Branch("genpt",genpt,"genpt[ngen]/F");
   t->Branch("geneta",geneta,"geneta[ngen]/F");
   t->Branch("genphi",genphi,"genphi[ngen]/F");
   
   double dRcut=cone;
   dRcut/=10;
   for (int i=0;i<fMB.Jet->GetEntries();i++) {
      fMB.Jet->GetEntry(i);
      if (i%100000==0) cout <<i<<" / "<<fMB.Jet->GetEntries()<<endl;
      for (int j=0;j<nref;j++) {
         jtgenref[j]=-1;
	 jtDR[j]=-1;
         for (int k=0;k<ngen;k++) {
	    double dR=deltaR(jteta[j],jtphi[j],geneta[k],genphi[k]);
	    if (dR<dRcut&&genpt[k]>50) {
	       jtgenref[j]=k;
	       jtDR[j]=dR;
	       break;   
	    }
	 }
      }
      t->Fill();
   }
   t->Write();
   outf->Close();   
}
