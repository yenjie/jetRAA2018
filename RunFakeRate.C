#include "plotFakeRate.C"
#include "genMatching.C"
#include "/afs/cern.ch/user/s/skanaski/jet_stuff/SKDetector/TMVASlides.C"

void RunFakeRate()
{
  CustomCanvas* c2=new CustomCanvas("c2","Fake Rate",600,600);
  
  std::cout<<"genMatching:"<<std::endl<<"Cone size 0.3:"<<std::endl;
  genMatching(3);
  std::cout<<"Cone size 0.4:"<<std::endl;
  genMatching(4);
  std::cout<<"Cone size 0.6:"<<std::endl;
  genMatching(6);
  std::cout<<"Cone size 0.8:"<<std::endl;
  genMatching(8);
  std::cout<<"Cone size 1.0:"<<std::endl;
  genMatching(10);

  gSystem->cd("results");
  std::cout<<"Creating fake rate histograms"<<std::endl;
  doAll(c2);

  std::cout<<"Generating slides"<<std::endl;
  TMVASlides(new std::vector<std::vector<std::string>*>{c2->GetPointer()},"FakeRateSlides.tex",1);
  delete c2;
  gSystem->cd("..");
}
