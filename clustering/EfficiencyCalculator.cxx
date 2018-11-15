
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
  double binning[1000];
  double clusters[1000];
  double numEvents[1000];
  double efficiency[1000];
  
  double SingleNeutronEfficiencyNumerator=0;
  double SingleNeutronEfficiencyDenominator=0;
  double MultipleNeutronEfficiencyNumerator=0;
  double MultipleNeutronEfficiencyDenominator=0;
  
  std::map<int, std::vector<int>> Indices;
 
  std::vector<int>      * HitView           = NULL;
  std::vector<int>      * HitChan           = NULL;
  std::vector<int>      * GenType           = NULL;
  std::vector<int>      * GenParticleID     = NULL;
  std::vector<int>      * BackTrackedID     = NULL;
  std::vector<int>      * GenParticlePDG    = NULL;
  std::vector<int>      * EventIndex        = NULL;
  std::vector<int>      * NeutOneIds        = NULL;
  std::vector<int>      * NeutTwoIds        = NULL;
  std::vector<int>      * NeutThreeIds      = NULL;
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
  std::vector<double>   * NeutOneY          = NULL;
  std::vector<double>   * NeutOneZ          = NULL;


  // SETTING ALL OF THE BRANCH ADDRESSES
  Tree->SetBranchAddress("NTotHits"         , &NTotHits           );
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
  Tree->SetBranchAddress("NeutOneX"         , &NeutOneX           );
  Tree->SetBranchAddress("NeutOneY"         , &NeutOneY           );
  Tree->SetBranchAddress("NeutOneZ"         , &NeutOneZ           );
  Tree->SetBranchAddress("EventIndex"       , &EventIndex         );
  Tree->SetBranchAddress("NeutOneIds"       , &NeutOneIds         );
  Tree->SetBranchAddress("NeutTwoIds"       , &NeutTwoIds         );
  Tree->SetBranchAddress("NeutThreeIds"     , &NeutThreeIds       );
    
  
//  Tree->SetBranchAddress("Index"             , &GenType           ); //

  // THESE ARE THE OUTPUTS I WANT TO GET FROM THE CODE
//  TProfile *tprof_nClusterVSnNeutron = new TProfile("", "", nHitTrials, -0.5, 9.5);
//  TH2D *th2d_nClusterVSnNeutron      = new TH2D("", "", 10, -0.5, 9.5, 10, -0.5, 9.5);
//  th2d_nClusterVSnNeutron->SetStats(0);
  

  for (int i=0; i<1; ++i){
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

    wiretrigger.SetNHitsMin    (6); //efficiency vs this plot, 6 is the nominal
    wiretrigger.SetNChanMin    (2);
    wiretrigger.SetChanWidthMin(0);
    wiretrigger.SetSADCMin     (0);

    // THIS SETS UP THE NUMBER OF EVENTS
    int fNEvent = nEvent;
    if (fNEvent!=-1) {fNEvent = std::min(fNEvent,(int)Tree->GetEntries());}
    else             {fNEvent = Tree->GetEntries();}

    // CALCULATING OUTPUT VALUES (I THINK. THIS CODE IS F**KING COMPLICATED)
//    for (int CurrentEvent=185; CurrentEvent<186; ++CurrentEvent) {
    for (int CurrentEvent=0; CurrentEvent<fNEvent; ++CurrentEvent) {
      Tree->GetEntry(CurrentEvent);
      std::vector<WireHit*> vec_WireHit;
      numberOfEvents += (*EventIndex)[0];
      
      Indices[0] = (*NeutOneIds)  ;
      Indices[1] = (*NeutTwoIds)  ;
      Indices[2] = (*NeutThreeIds);
      
//      for (auto it : Indices) {
//        for (auto c : it.second) {
//          std::cout << it.first << "\t" << c << std::endl;
//        }
//      }
      
      for (int j=0; j<HitView->size(); ++j) {

        // HERE WE NEED SOMETHING LIKE: IF THE GENPARTICLEID MATCHES THE
        // BACKTRACKED ID FOR THAT PARTICLE ON THAT ROUND MAKE THE SECOND
        // ARGUMENT 1, OTHERWISE, IT WILL BE ZERO...
        WireHit* hit = new WireHit((*HitView)[j]   , 1            , (*HitChan)[j],
                                   (*HitTime)[j]   , (*HitSADC)[j], (*HitRMS)[j] ,
                                   0               , 0            , 0            ,
                                   0               , 0            , 0            ,
                                   0               , 0            , 0            ,
                                   (*EventIndex)[j], 0           );
        vec_WireHit.push_back(hit);
      }

      bool selected      = false;
      int  ncluster      = 0    ;
      int  nnoisecluster = 0    ;


      std::vector<WireCluster*> vec_WireCluster;
      clusteng    .ClusterHits2 (vec_WireHit, vec_WireCluster);
      wiretrigger .SetIsSelected(vec_WireCluster)             ;

      std::vector<WireCluster*> temp_VecWireCluster;
      double first, second, third;
      
      if (vec_WireCluster.size() > (*EventIndex)[0] && (*EventIndex)[0] != 0){
        
        if ((*EventIndex)[0] >= 1) {
          for (int it=0; it<vec_WireCluster.size(); ++it) {
            if (vec_WireCluster[it]->GetSumPeak() > first) {
              first = vec_WireCluster[it]->GetSumPeak();
            }
          }
        }
        if ((*EventIndex)[0] >= 2) {
          for (int it=0; it<vec_WireCluster.size(); ++it) {
            if (vec_WireCluster[it]->GetSumPeak() > second && vec_WireCluster[it]->GetSumPeak() < first) {
              second = vec_WireCluster[it]->GetSumPeak();
            }
          }
        }
        if ((*EventIndex)[0] >= 3) {
          for (int it=0; it<vec_WireCluster.size(); ++it) {
            if (vec_WireCluster[it]->GetSumPeak() > third && vec_WireCluster[it]->GetSumPeak() < second && vec_WireCluster[it]->GetSumPeak() < first){
              third = vec_WireCluster[it]->GetSumPeak();
            }
          }
        }
        for (int i=0; i<vec_WireCluster.size(); ++i) {
          if (vec_WireCluster[i]->GetSumPeak() == first) { temp_VecWireCluster.push_back(vec_WireCluster[i]); }
        }
        for (int i=0; i<vec_WireCluster.size(); ++i) {
          if (vec_WireCluster[i]->GetSumPeak() == second) { temp_VecWireCluster.push_back(vec_WireCluster[i]); }
        }
        for (int i=0; i<vec_WireCluster.size(); ++i) {
          if (vec_WireCluster[i]->GetSumPeak() == third) { temp_VecWireCluster.push_back(vec_WireCluster[i]); }
        }
        
      }
      vec_WireCluster = temp_VecWireCluster;
      temp_VecWireCluster.clear();
      first  = 0; second = 0; third  = 0;
      
      std::map<int, bool> selected_cluster;
      for (int c=0; c<(*EventIndex)[0]; ++c) {
        selected_cluster[c] = false;
      }

      for (int c=0; c<vec_WireCluster.size(); ++c) {
        WireCluster* clust = vec_WireCluster[c];
        if (clust->GetIsSelected()) {
          if (clust->GetType()){
            selected_cluster[clust->GetMarleyIndex()] = true;
            selected = true;
            ++ncluster;
          }
          else {
            ++nnoisecluster;
          }
        }
      }

      //////////////////////////////////////
      //      CALCULATE ME SOME SHIT      //
      //////////////////////////////////////
      
      int NumberOfNeutronsDetected = 0;
      for (auto const& it: selected_cluster) {
        if (it.second) {
          SingleNeutronEfficiencyNumerator++;
          NumberOfNeutronsDetected++;
        }
        SingleNeutronEfficiencyDenominator++;
      }

      MultipleNeutronEfficiencyDenominator++;
      if (NumberOfNeutronsDetected == selected_cluster.size() && selected_cluster.size() != 0) {
        MultipleNeutronEfficiencyNumerator++;
      }
 
      
      totalNumberOfClusters+=ncluster;
      totalNoise+=nnoisecluster;
//     Calcuklate the efficiency
//     and fill the hist and the tProfile
//     th2d_nClusterVSnNeutron ->Fill(3, ncluster);
//     tprof_nClusterVSnNeutron->Fill(3, ncluster);
    
    binning[CurrentEvent]    = static_cast<double>(CurrentEvent);
    clusters[CurrentEvent]   = static_cast<double>(totalNumberOfClusters);
    numEvents[CurrentEvent]  = static_cast<double>(numberOfEvents);
    efficiency[CurrentEvent] = static_cast<double>(100*totalNumberOfClusters)/static_cast<double>(numberOfEvents);

      
      

    }
  }
  
  std::cout << "\nSingle Neutron Efficiency: \tMultiple Neutron Efficiency: " << std::endl;
  std::cout << "-----------------------------------------------------------------------" << std::endl;
  std::cout << SingleNeutronEfficiencyNumerator << "\t\t\t\t" << MultipleNeutronEfficiencyNumerator << std::endl;
  std::cout << "--\t\t\t\t--" << std::endl;
  std::cout << SingleNeutronEfficiencyDenominator << "\t\t\t\t" << MultipleNeutronEfficiencyDenominator << std::endl;
  std::cout << "-----------------------------------------------------------------------" << std::endl;
  std::cout << SingleNeutronEfficiencyNumerator/SingleNeutronEfficiencyDenominator
            << "\t\t\t" << MultipleNeutronEfficiencyNumerator/MultipleNeutronEfficiencyDenominator << std::endl;
  
  
  
  
  
  
  std::cout << "Event\tnClustFound\tnEventsTrue\teff" << std::endl;
  std::cout << "----------------------------------------------------" << std::endl;
  for (int i=0; i<1000; i++){
    std::cout << binning[i] << "\t" << clusters[i] << "\t\t" << numEvents[i] << "\t\t" << efficiency[i] << std::endl;
  }
  
  
//  TCanvas *c1 = new TCanvas("c1");
//  TLegend *l1 = new TLegend(0.4, 0.78, 0.7, 0.9);
//  TGraph  *g1 = new TGraph (1000, binning, clusters);
//  TGraph  *g2 = new TGraph (1000, binning, efficiency);
//  TGraph  *g3 = new TGraph (1000, binning, numEvents);
//  c1->Print((OutFileName+"[").c_str());
//  g1->SetTitle("");
//  g1->GetXaxis()->SetTitle("Sample number [#]");
//  g1->GetYaxis()->SetTitle("nClusters (black); efficiency (red)");
//  g1->Draw();
//  g2->SetLineColor(kRed);
//  g2->Draw("SAME");
//  l1->SetHeader("nHitsMin = 6", "L");
//  l1->AddEntry(g1, "Number of clusters vs number of samples tested");
//  l1->AddEntry(g2, "Efficiency of clusteringvs number of samples tested (%)");
//  l1->Draw();
//  c1->Print((OutFileName).c_str());
//  c1->Print((OutFileName+"]").c_str());
  

}


