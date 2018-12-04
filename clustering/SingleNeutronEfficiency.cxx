
#include "ArbitraryAnaInputManager.hh"
#include "ClusterEngine.hh"
#include "Cluster.hh"
#include "Trigger.hh"

#include "Helper.h"

#include "TCanvas.h"
#include "TEfficiency.h"
#include "TFile.h"
#include "TPad.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TProfile.h"
#include "TGraph.h"
#include "TLine.h"
#include "TTree.h"
#include "WireHit.hh"

#include <iostream>
#include <unistd.h>

int main (int argc, char** argv) {
  int opt;
  // Shut GetOpt error messages down (return '?'):
  extern char *optarg;
  extern int  optopt;
  
  std::string InFileName ="/Users/aranborkum/Docterate/DataSets/100_Neutron_background.root";
  std::string OutFileName="Output_1.pdf";
  std::string InTreeName ="trineutronbackground/TriNeutron";
  
  int nEvent=100;
  int nHitTrials=1;
  while ( (opt = getopt(argc, argv, "i:o:n:t:l:")) != -1 ) {  // for each option...
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
    case 'l':
      nHitTrials = atoi(optarg);
      break;
    case '?':  // unknown option...
      std::cerr << "Unknown option: '" << char(optopt) << "'!" << std::endl;
      break;
    }
  }
  
  TFile* InFile = new TFile(InFileName.c_str(),"READ");
  TTree* Tree = (TTree*)InFile->Get(InTreeName.c_str());
  
  double SingleNeutronEfficiencyNumerator   = 0;
  double SingleNeutronEfficiencyDenominator = 0;
  double NumberOfNeutronsDetected           = 0;
  
  std::vector<int> GoodNeutrons;
  
  std::vector<int>      * ThisNeutronID     = NULL;
  std::vector<int>      * HitView           = NULL;
  std::vector<int>      * HitChan           = NULL;
  std::vector<int>      * GenType           = NULL;
  std::vector<int>      * GenParticleID     = NULL;
  std::vector<int>      * BackTrackedID     = NULL;
  std::vector<int>      * GenParticlePDG    = NULL;
  std::vector<int>      * EventIndex        = NULL;
  std::vector<int>      * NeutronIndex      = NULL;
  std::vector<int>      * ProcessEnding       = NULL;

  std::vector<bool>     * ThisNeutronStaysInTheDetector = NULL;
  std::vector<bool>     * ThisNeutronCapture = NULL;
  std::vector<float>    * HitTime           = NULL;
  std::vector<float>    * HitSADC           = NULL;
  std::vector<float>    * HitRMS            = NULL;
  std::vector<float>    * IonizaitonEnergy  = NULL;
  std::vector<float>    * NumberOfElectrons = NULL;
  std::vector<float>    * Hit_X             = NULL;
  std::vector<float>    * Hit_Y             = NULL;
  std::vector<float>    * Hit_Z             = NULL;
  std::vector<double>   * GenParticleEnergy = NULL;
  std::vector<double>   * GenParticleStartX = NULL;
  std::vector<double>   * GenParticleStartY = NULL;
  std::vector<double>   * GenParticleStartZ = NULL;
  std::vector<double>   * GenParticleEndX   = NULL;
  std::vector<double>   * GenParticleEndY   = NULL;
  std::vector<double>   * GenParticleEndZ   = NULL;
  

  Tree->SetBranchAddress("HitView"          , &HitView            );
  Tree->SetBranchAddress("HitChan"          , &HitChan            );
  Tree->SetBranchAddress("HitTime"          , &HitTime            );
  Tree->SetBranchAddress("HitSADC"          , &HitSADC            );
  Tree->SetBranchAddress("HitRMS"           , &HitRMS             );
  Tree->SetBranchAddress("GenType"          , &GenType            );
  Tree->SetBranchAddress("GenParticleID"    , &GenParticleID      );
  Tree->SetBranchAddress("BackTrackedID"    , &BackTrackedID      );
  Tree->SetBranchAddress("GenParticlePDG"   , &GenParticlePDG     );
  Tree->SetBranchAddress("IonizaitonEnergy" , &IonizaitonEnergy   );
  Tree->SetBranchAddress("NumberOfElectrons", &NumberOfElectrons  );
  Tree->SetBranchAddress("GenParticleEnergy", &GenParticleEnergy  );
  Tree->SetBranchAddress("Hit_X"            , &Hit_X              );
  Tree->SetBranchAddress("Hit_Y"            , &Hit_Y              );
  Tree->SetBranchAddress("Hit_Z"            , &Hit_Z              );
  Tree->SetBranchAddress("GenParticleStartX", &GenParticleStartX  );
  Tree->SetBranchAddress("GenParticleStartY", &GenParticleStartY  );
  Tree->SetBranchAddress("GenParticleStartZ", &GenParticleStartZ  );
  Tree->SetBranchAddress("GenParticleEndX"  , &GenParticleEndX    );
  Tree->SetBranchAddress("GenParticleEndY"  , &GenParticleEndY    );
  Tree->SetBranchAddress("GenParticleEndZ"  , &GenParticleEndZ    );
  Tree->SetBranchAddress("EventIndex"       , &EventIndex         );
  Tree->SetBranchAddress("NeutronIndex"     , &NeutronIndex       );
  Tree->SetBranchAddress("ThisNeutronStaysInTheDetector", &ThisNeutronStaysInTheDetector);
  Tree->SetBranchAddress("ThisNeutronCapture", &ThisNeutronCapture);
  Tree->SetBranchAddress("ThisNeutronID"     , &ThisNeutronID     );
  Tree->SetBranchAddress("ProcessEnding"      , &ProcessEnding      );
  
  ClusterEngine clusteng;
  SimpleWireTrigger wiretrigger;
  
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
  
  
  
  for (int i=0; i<Tree->GetEntries(); ++i){
    
    Tree->GetEntry(i);
    
    if ((*ThisNeutronStaysInTheDetector)[0] == 0 || (*ThisNeutronCapture)[0] == 0)
      continue;
    
    std::vector<WireHit*> vec_WireHit;
    for (int j=0; j<HitView->size(); ++j) {
        
      // HERE WE NEED SOMETHING LIKE: IF THE GENPARTICLEID MATCHES THE
      // BACKTRACKED ID FOR THAT PARTICLE ON THAT ROUND MAKE THE SECOND
      // ARGUMENT 1, OTHERWISE, IT WILL BE ZERO...
      WireHit* hit = new WireHit((*HitView)[j]     , 1            , (*HitChan)[j],
				 (*HitTime)[j]     , (*HitSADC)[j], (*HitRMS)[j] ,
				 0                 , 0            , 0            ,
				 0                 , 0            , 0            ,
				 0                 , 0            , 0            ,
				 (int)(*NeutronIndex)[j], 0           );
      vec_WireHit.push_back(hit);
    }
      
    bool selected     = false;
    int ncluster      = 0    ;
    int nnoisecluster = 0    ;
      
    std::vector<WireCluster*> vec_WireCluster;
    clusteng.ClusterHits2    (vec_WireHit, vec_WireCluster);
   
    for (auto& it: vec_WireCluster) {
      std::map<int, int> number_of_hits;
      for (auto const& hit: it->GetHit()) {number_of_hits[hit->GetMarleyIndex()]++;}
      auto max = std::max_element(number_of_hits.begin(), number_of_hits.end());
      it->SetMarleyIndex(max->first);
    }
    
    wiretrigger.SetIsSelected(vec_WireCluster);
    
    std::map<int, bool> selected_cluster;
    for (int c=0; c<NeutronIndex->size(); ++c) {
      selected_cluster[c] = false;
    }
    
    for (int c=0; c<vec_WireCluster.size(); ++c) {
      WireCluster* clust = vec_WireCluster[c];
      if (clust->GetIsSelected()) {
        if (clust->GetType()){
          selected_cluster[clust->GetMarleyIndex()] = true;
          selected = true;
          ++ncluster;
        } else {
          ++nnoisecluster;
        }
      }
    }

    int NumberOfNeutronsDetected = 0;
    std::map<int, bool> DetectedNeutron;
    DetectedNeutron.clear();
    for (auto const& it: (*ThisNeutronID)) {
      DetectedNeutron[it] = false;
    }
    
    for (auto const& it: vec_WireCluster) {
      if (it->GetIsSelected() && DetectedNeutron.find(it->GetMarleyIndex()) != DetectedNeutron.end()) {
        for (auto const& hit: it->GetHit()){
          DetectedNeutron[hit->GetMarleyIndex()] = true;
        }
      }
    }
    
    for (auto const& it: DetectedNeutron) {
      if (it.second) {
        SingleNeutronEfficiencyNumerator++;
        NumberOfNeutronsDetected++;
      }
      SingleNeutronEfficiencyDenominator++;
    }
    
    
    
    
    
    
    
//    RAISE YOUR HANDS FOR MEMORY MANAGEMENT
    
    for (auto& it: vec_WireHit) {
      delete it;
      it = NULL;
    }
    vec_WireHit.clear();
  }
  
  std::cout << SingleNeutronEfficiencyNumerator << "/" << SingleNeutronEfficiencyDenominator << " = " << SingleNeutronEfficiencyNumerator/SingleNeutronEfficiencyDenominator << std::endl;
  
  
}

  
  

