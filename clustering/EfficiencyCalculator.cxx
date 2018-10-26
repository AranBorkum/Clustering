
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

  std::string InFileName ="/Users/aranborkum/Desktop/SimplePGunAna_hist.root";
  std::string OutFileName="OutputEfficiency.pdf";
  std::string InTreeName ="simplepgunana/PGun";
  
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
  
  std::vector<int>   * HitView        = NULL;
  std::vector<int>   * HitChan        = NULL;
  std::vector<int>   * GenType        = NULL;
  std::vector<int>   * GenParticleID  = NULL;
  std::vector<int>   * BackTrackedID  = NULL;
  std::vector<int>   * GenParticlePDG = NULL;
  std::vector<float> * HitTime        = NULL;
  std::vector<float> * HitSADC        = NULL;
  std::vector<float> * HitRMS         = NULL;
  
  // SETTING ALL OF THE BRANCH ADDRESSES
  Tree->SetBranchAddress("NTotHits"          , &NTotHits          );
  Tree->SetBranchAddress("HitView"           , &HitView           );
  Tree->SetBranchAddress("HitChan"           , &HitChan           );
  Tree->SetBranchAddress("HitTime"           , &HitTime           );
  Tree->SetBranchAddress("HitSADC"           , &HitSADC           );
  Tree->SetBranchAddress("HitRMS"            , &HitRMS            );
  Tree->SetBranchAddress("GenType"           , &GenType           );
  Tree->SetBranchAddress("GenParticleID"     , &GenParticleID     );
  Tree->SetBranchAddress("BackTrackedID"     , &BackTrackedID     );
  Tree->SetBranchAddress("GenParticlePDG"    , &GenParticlePDG    );
 
  // THESE ARE THE OUTPUTS I WANT TO GET FROM THE CODE
  TProfile*         tprof_nClusterVSnNeutron;
  TH2D*             th2d_nClusterVSnNeutron ;
  double            efficiency              ;
  ClusterEngine     clusteng                ;
  SimpleWireTrigger wiretrigger             ;

  // NOT SURE WHAT THESE DO, NOT GOING TO WORRY ABOUT IT TOO MUCH
  clusteng.SetTimeWindow   (20) ;
  clusteng.SetChannelWidth (2)  ;
  clusteng.SetTimeWindowOpt(0.2);
  clusteng.SetPositionOpt  (300);
  clusteng.SetBucketSize   (1)  ;
  
  wiretrigger.SetNHitsMin    (6);
  wiretrigger.SetNChanMin    (2);
  wiretrigger.SetChanWidthMin(0);
  wiretrigger.SetSADCMin     (0);

  // THIS SETS UP THE NUMBER OF EVENTS
  int fNEvent = nEvent;
  if (fNEvent!=-1) {fNEvent = std::min(fNEvent,(int)Tree->GetEntries());}
  else             {fNEvent = Tree->GetEntries();}
  
  // CONSTRUCTRORS FOR THE VISUAL OUTPUT OF THIS PROGRAM
  tprof_nClusterVSnNeutron = new TProfile("", "", 10, -0.5, 9.5);
  th2d_nClusterVSnNeutron  = new TH2D("", "", 10, -0.5, 9.5, 10, -0.5, 9.5);
  th2d_nClusterVSnNeutron->SetStats(0);
  
  // CALCULATING OUTPUT VALUES (I THINK. THIS CODE IS F**KING COMPLICATED)
  for (int CurrentEvent=0; CurrentEvent<fNEvent; ++CurrentEvent) {
    Tree->GetEntry(CurrentEvent);

    std::vector<WireHit*> vec_WireHit;

    for (int j=0; j<HitView->size(); ++j) {
      
      // HERE WE NEED SOMETHING LIKE: IF THE GENPARTICLEID MATCHES THE
      // BACKTRACKED ID FOR THAT PARTICLE ON THAT ROUND MAKE THE SECOND
      // ARGUMENT 1, OTHERWISE, IT WILL BE ZERO...
      
      WireHit* hit = new WireHit((*HitView)[j], (*GenType)[j], (*HitChan)[j],
                                 (*HitTime)[j], (*HitSADC)[j], (*HitRMS)[j] ,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
      vec_WireHit.push_back(hit);
    }
    
      bool selected      = false;
      int  ncluster      = 0    ;
      int  nnoisecluster = 0    ;
    
      std::vector<WireCluster*> vec_WireCluster;
    
      clusteng    .ClusterHits2 (vec_WireHit, vec_WireCluster);
      wiretrigger .SetIsSelected(vec_WireCluster)             ;

    for (int c=0; c<vec_WireCluster.size(); ++c) {
      WireCluster* clust = vec_WireCluster[c];
      if (clust->GetIsSelected()) {
        if (clust->GetType()){
          selected = true;
          ++ncluster;
        }
        else {
          ++nnoisecluster;
        }
        std::cout << "Cluster: " << ncluster << "\tNoise: " << nnoisecluster << std::endl;
      }
      
      
    }

      // Calcuklate the efficiency
      // anbd fill the hist and the tProfile
       th2d_nClusterVSnNeutron ->Fill(3, nnoisecluster);
       tprof_nClusterVSnNeutron->Fill(3, nnoisecluster);
    
  }
  
  std::cout << "\nThe efficiency of clustering 1 neutron is:" << efficiency << std::endl;

  std::cout << std::endl;
  
  TCanvas c;
  c.Print((OutFileName+"[").c_str());
  th2d_nClusterVSnNeutron->Draw("COLZ");
  tprof_nClusterVSnNeutron->Draw("SAME");
  c.Print(OutFileName.c_str());

  c.Print((OutFileName+"]").c_str());
  
}


