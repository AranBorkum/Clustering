
#include "ArbitraryAnaInputManager.hh"
#include "ClusterEngine.hh"
#include "Trigger.hh"

#include "Helper.h"

#include "TCanvas.h"
#include "TEfficiency.h"
#include "TFile.h"
#include "TPad.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TProfile.h"
#include "TTree.h"
#include "WireHit.hh"

#include <iostream>
#include <unistd.h>

int main (int argc, char** argv) {
  int opt;
  // Shut GetOpt error messages down (return '?'):
  extern char *optarg;
  extern int  optopt;

  std::string InFileName ="";
  std::string OutFileName="OutputEfficiency.pdf";
  std::string InTreeName ="";
  int nEvent=-1;
  while ( (opt = getopt(argc, argv, "i:o:n:t:")) != -1 ) {  // for each option...
    switch ( opt ) {
    case 'i':
      InFileName = optarg;
      break;
    case 'o':
      OutFileName = optarg;
      break;
    case 't':
        InTreeName = optarg;
        break;
    case 'n':
      nEvent = atoi(optarg);
      break;
    case '?':  // unknown option...
      std::cerr << "Unknown option: '" << char(optopt) << "'!" << std::endl;
      break;
    }
  }
  
  TFile* InFile = new TFile(InFileName.c_str(),"READ");
  TTree* Tree = (TTree*)InFile->Get(InTreeName.c_str());
  
  // DECLARATION OF ALL THE STORAGE VARIABLES
  int NTotHits;
  
  std::vector<int>   * HitView = NULL;
  std::vector<int>   * HitChan = NULL;
  std::vector<float> * HitTime = NULL;
  std::vector<float> * HitSADC = NULL;
  std::vector<float> * HitRMS  = NULL;
  
  
  // SETTING ALL OF THE BRANCH ADDRESSES
  Tree->SetBranchAddress("NTotHits"          , &NTotHits          );
  Tree->SetBranchAddress("HitView"           , &HitView           );
  Tree->SetBranchAddress("HitChan"           , &HitChan           );
  Tree->SetBranchAddress("HitTime"           , &HitTime           );
  Tree->SetBranchAddress("HitSADC"           , &HitSADC           );
  Tree->SetBranchAddress("HitRMS"            , &HitRMS            );

  
  // THESE ARE THE OUTPUTS I WANT TO GET FROM THE CODE
  TProfile* tprof_nClusterVSnNeutron;
  TH2D*     th2d_nClusterVSnNeutron;
  double efficiency;
  ClusterEngine clusteng;
  SimpleWireTrigger wiretrigger;

  // NOT SURE WHAT THESE DO, NOT GOING TO WORRY ABOUT IT TOO MUCH
  clusteng.SetTimeWindow   (20);
  clusteng.SetChannelWidth (2);
  clusteng.SetTimeWindowOpt(0.2);
  clusteng.SetPositionOpt  (300);
  clusteng.SetBucketSize   (1);
  
  wiretrigger.SetNHitsMin    (6);
  wiretrigger.SetNChanMin    (2);
  wiretrigger.SetChanWidthMin(0);
  wiretrigger.SetSADCMin     (0);

  int fNEvent = nEvent;
  if (fNEvent!=-1) {fNEvent = std::min(fNEvent,(int)Tree->GetEntries());}
  else             {fNEvent = Tree->GetEntries();}
  
  // CONSTRUCTRORS FOR THE VISUAL OUTPUT OF THIS PROGRAM
  tprof_nClusterVSnNeutron = new TProfile("nSignalCluster", ";electron KE [MeV];nSignalCluster", 10, -0.5, 9.5);
  th2d_nClusterVSnNeutron  = new TH2D("nSignalCluste_th2" , "Threshold ADC;electron KE [MeV];nSignalCluster;Event", 10, -0.5, 9.5, 10, -0.5, 9.5);
  th2d_nClusterVSnNeutron->SetStats(0);
  
  // CALCULATING OUTPUT VALUES (I THINK. THIS CODE IS F**KING COMPLICATED)
  for (int CurrentEvent=0; CurrentEvent<fNEvent; ++CurrentEvent) {
    Tree->GetEntry(CurrentEvent);
    std::vector<WireHit*> vec_WireHit;

    for (int j=0; j<HitView->size(); ++j) {
      WireHit* hit = new WireHit((*HitView)[j], (*HitChan)[j] , (*HitTime)[j],
                                 (*HitSADC)[j], (*HitRMS)[j]);
      vec_WireHit.push_back(hit);
    }
    
      bool selected = false;
      int ncluster = 0;
      int nnoisecluster = 0;
      std::vector<WireCluster*> vec_WireCluster;
      clusteng.ClusterHits2(vec_WireHit, vec_WireCluster);
      wiretrigger.SetIsSelected(vec_WireCluster);

//    for (int c=0; c<vec_WireCluster.size(); ++c) {
//      // Calcuklate the efficiency
//      // anbd fill the hist and the tProfile
//      // th2d_nClusterVSnNeutron->Fill(nNeutronGenerated, nClusters);
//      // tprof_nClusterVSnNeutron->Fill(nNeutronGenerated, nClusters);
//      }
    }

  std::cout << "The efficiency of clustering 1 neutron is:" << efficiency << std::endl;

  TCanvas c;
  c.Print((OutFileName+"[").c_str());
  th2d_nClusterVSnNeutron->Draw("COLZ");
  tprof_nClusterVSnNeutron->Draw("SAME");
  c.Print(OutFileName.c_str());

  c.Print((OutFileName+"]").c_str());
  
}


