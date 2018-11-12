
#include "ArbitraryAnaInputManager.hh"
#include "ClusterEngine.hh"
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

  std::string InFileName ="/Users/aranborkum/Docterate/DataSets/1000_tri_Neutron_events.root";
  std::string OutFileName="Output_1.pdf";
  std::string InTreeName ="trineutronbackground/TriNeutron";
  
  int nEvent=1000;
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
  
  // DECLARATION OF ALL THE STORAGE VARIABLES
  int NTotHits;
  double binning[8];
  double clusters[8];
  double efficiency[8];
  
  std::vector<int>      * HitView           = NULL;
  std::vector<int>      * HitChan           = NULL;
  std::vector<int>      * GenType           = NULL;
  std::vector<int>      * GenParticleID     = NULL;
  std::vector<int>      * BackTrackedID     = NULL;
  std::vector<int>      * GenParticlePDG    = NULL;
  std::vector<int>      * EventIndex        = NULL;
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
  std::vector<double>   * NeutOneX          = NULL;

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
  Tree->SetBranchAddress("IonizaitonEnergy"  , &IonizaitonEnergy  );
  Tree->SetBranchAddress("NumberOfElectrons" , &NumberOfElectrons );
  Tree->SetBranchAddress("GenParticleEnergy" , &GenParticleEnergy );
  Tree->SetBranchAddress("Hit_X"             , &Hit_X             );
  Tree->SetBranchAddress("Hit_Y"             , &Hit_Y             );
  Tree->SetBranchAddress("Hit_Z"             , &Hit_Z             );
  Tree->SetBranchAddress("GenParticleStartX", &GenParticleStartX  );
  Tree->SetBranchAddress("GenParticleStartY", &GenParticleStartY  );
  Tree->SetBranchAddress("GenParticleStartZ", &GenParticleStartZ  );
  Tree->SetBranchAddress("GenParticleEndX"  , &GenParticleEndX    );
  Tree->SetBranchAddress("GenParticleEndY"  , &GenParticleEndY    );
  Tree->SetBranchAddress("GenParticleEndZ"  , &GenParticleEndZ    );
  Tree->SetBranchAddress("EventIndex"       , &EventIndex         );
  Tree->SetBranchAddress("NeutOneX"         , &NeutOneX           );
//  Tree->SetBranchAddress("Index"             , &GenType           ); //

  // THESE ARE THE OUTPUTS I WANT TO GET FROM THE CODE
//  TProfile *tprof_nClusterVSnNeutron = new TProfile("", "", nHitTrials, -0.5, 9.5);
//  TH2D *th2d_nClusterVSnNeutron      = new TH2D("", "", 10, -0.5, 9.5, 10, -0.5, 9.5);
//  th2d_nClusterVSnNeutron->SetStats(0);

  for (int i=0; i<8; ++i){
    ClusterEngine                        clusteng;
    SimpleWireTrigger                    wiretrigger;
    
    int numberOfEvents = 0;
    int totalNumberOfClusters = 0;
    int totalNoise = 0;
    
    // NOT SURE WHAT THESE DO, NOT GOING TO WORRY ABOUT IT TOO MUCH
    clusteng.SetTimeWindow   (20) ;
    clusteng.SetChannelWidth (2)  ;
    clusteng.SetTimeWindowOpt(0.2);
    clusteng.SetPositionOpt  (300);
    clusteng.SetBucketSize   (1)  ;

    wiretrigger.SetNHitsMin    (i); //efficiency vs this plot, 6 is the nominal
    wiretrigger.SetNChanMin    (2);
    wiretrigger.SetChanWidthMin(0);
    wiretrigger.SetSADCMin     (0);

    // THIS SETS UP THE NUMBER OF EVENTS
    int fNEvent = nEvent;
    if (fNEvent!=-1) {fNEvent = std::min(fNEvent,(int)Tree->GetEntries());}
    else             {fNEvent = Tree->GetEntries();}

    // CALCULATING OUTPUT VALUES (I THINK. THIS CODE IS F**KING COMPLICATED)
    for (int CurrentEvent=18; CurrentEvent<fNEvent; ++CurrentEvent) {
      Tree->GetEntry(CurrentEvent);
      
      if ((NeutOneX->size() != 0)) {
      

        std::vector<WireHit*> vec_WireHit;
        numberOfEvents++;
        for (int j=0; j<HitView->size(); ++j) {

          // HERE WE NEED SOMETHING LIKE: IF THE GENPARTICLEID MATCHES THE
          // BACKTRACKED ID FOR THAT PARTICLE ON THAT ROUND MAKE THE SECOND
          // ARGUMENT 1, OTHERWISE, IT WILL BE ZERO...
          WireHit* hit = new WireHit((*HitView)[j], 1            , (*HitChan)[j],
                                     (*HitTime)[j], (*HitSADC)[j], (*HitRMS)[j] ,
                                     0            , 0            , 0            ,
                                     0            , 0            , 0            ,
                                     0            , 0            , 0            ,
                                     0/* index here*/            , 0           );
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
          }
        }

        totalNumberOfClusters+=ncluster;
        totalNoise+=nnoisecluster;
    // Calcuklate the efficiency
    // and fill the hist and the tProfile
    // th2d_nClusterVSnNeutron ->Fill(3, ncluster);
    // tprof_nClusterVSnNeutron->Fill(3, ncluster);
      }
    }
    
    binning[i] = static_cast<double>(i);
    clusters[i] = static_cast<double>(totalNumberOfClusters);
    efficiency[i] = static_cast<double>(totalNumberOfClusters)/static_cast<double>(numberOfEvents*3);
    
//    std::cout << "#clusters: " << totalNumberOfClusters
//              << "\t#Events: " << numberOfEvents
//              << "\t#efficiency: "  << efficiency[i]
//    << "\tMin Hits: "  << binning[i] << std::endl;
  }
  
  for (int i=0; i<8; ++i){
    std::cout << binning[i] << "\t" << efficiency[i] << std::endl;
  }
  
  

  TCanvas *c2 = new TCanvas();
  TGraph  *g2 = new TGraph (8, binning, efficiency);
//  TLine   *l2 = new TLine  (6, 0, 6, 0.9              );
  TLegend *legend2 = new TLegend(0.52, 0.7, 0.9, 0.9  );

  legend2->AddEntry(g2, "efficiency vs minimum number of hits");
//  legend2->AddEntry(l2, "nominal number of hits");

  c2->Print((OutFileName+"[").c_str());
  c2->SetLogy();
  g2->SetTitle("");
  g2->GetXaxis()->SetTitle("nHitsMin [#]");
  g2->GetYaxis()->SetTitle("Efficiency");
  g2->Draw("AL*");
//  l2->SetLineColor(kRed);
//  l2->Draw();
  legend2->Draw();
  c2->Print((OutFileName).c_str());
  c2->Print((OutFileName+"]").c_str());

}


