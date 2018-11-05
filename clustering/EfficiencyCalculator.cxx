
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

  std::string InFileName ="/Users/aranborkum/Desktop/SimplePGunAna_hist.root";
  std::string OutFileName="Output_1.pdf";
  std::string InTreeName ="simplepgunana/PGun";
  
  int nEvent=-1;
  int nHitTrials;
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
  
  std::vector<int>      * HitView           = NULL;
  std::vector<int>      * HitChan           = NULL;
  std::vector<int>      * GenType           = NULL;
  std::vector<int>      * GenParticleID     = NULL;
  std::vector<int>      * BackTrackedID     = NULL;
  std::vector<int>      * GenParticlePDG    = NULL;
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
//  Tree->SetBranchAddress("Index"             , &GenType           ); //

  // THESE ARE THE OUTPUTS I WANT TO GET FROM THE CODE
  TProfile*         tprof_nClusterVSnNeutron;
  TH2D*             th2d_nClusterVSnNeutron ;
  ClusterEngine     clusteng                ;
  SimpleWireTrigger wiretrigger             ;
  
  int x[nHitTrials], y[nHitTrials];
  double efficiency[nHitTrials]   ;

  for (int i=0; i<nHitTrials; ++i){
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
  
    // CONSTRUCTRORS FOR THE VISUAL OUTPUT OF THIS PROGRAM
    tprof_nClusterVSnNeutron = new TProfile("", "", 10, -0.5, 9.5);
    th2d_nClusterVSnNeutron  = new TH2D("", "", 10, -0.5, 9.5, 10, -0.5, 9.5);
    int totalNumberOfClusters = 0;
    th2d_nClusterVSnNeutron->SetStats(0);
  
    // CALCULATING OUTPUT VALUES (I THINK. THIS CODE IS F**KING COMPLICATED)
    for (int CurrentEvent=0; CurrentEvent<fNEvent; ++CurrentEvent) {
      Tree->GetEntry(CurrentEvent);
      std::vector<WireHit*> vec_WireHit;

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
    // Calcuklate the efficiency
    // anbd fill the hist and the tProfile
    th2d_nClusterVSnNeutron ->Fill(3, ncluster);
    tprof_nClusterVSnNeutron->Fill(3, ncluster);
      
    }

    x[i] = i;
    y[i] = totalNumberOfClusters;
    efficiency[i] = static_cast<double>(totalNumberOfClusters)/static_cast<double>(3*nEvent);
    
  }
  
  TCanvas *c = new TCanvas();
  TGraph  *g = new TGraph (nHitTrials, x, y   );
  TLine   *l = new TLine  (3, 0, 3, 2500      );
  TLegend *L = new TLegend(0.52, 0.7, 0.9, 0.9);
  
  L->SetHeader("This is the Legend", "C");
  L->AddEntry (l, "number of neutrons produced per event");
  L->AddEntry (g, "number of clusters vs minimum number of hits");
  
  l->SetLineColor(kRed);
  c->Print((OutFileName+"[").c_str());
  g            ->SetTitle(""             );
  g->GetXaxis()->SetTitle("nHitsMin [#]" );
  g->GetYaxis()->SetTitle("nClusters [#]");
  g->Draw("AL*");
  l->Draw();
  L->Draw();
  c->Print((OutFileName).c_str());
  c->Print((OutFileName+"]").c_str());
//  TCanvas c;
//  c.Print((OutFileName+"[").c_str());
//  th2d_nClusterVSnNeutron->Draw("COLZ");
//  tprof_nClusterVSnNeutron->Draw("SAME");
//  c.Print(OutFileName.c_str());
//  c.Print((OutFileName+"]").c_str());
  
}


