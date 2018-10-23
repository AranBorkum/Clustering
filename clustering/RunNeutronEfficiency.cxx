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

#include <iostream>
#include <unistd.h>

int main (int argc, char** argv) {
  int opt;
  // Shut GetOpt error messages down (return '?'):
  extern char *optarg;
  extern int  optopt;
  
  std::string InFileName ="";
  std::string OutFileName="";
  int nEvent=-1;
  while ( (opt = getopt(argc, argv, "i:o:n:")) != -1 ) {  // for each option...
    switch ( opt ) {
      case 'i':
        InFileName = optarg;
        break;
      case 'o':
        OutFileName = optarg;
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
  std::map<int,ArbitraryAnaInputManager*> aaim;
  aaim[1600] = new ArbitraryAnaInputManager();
  aaim[2000] = new ArbitraryAnaInputManager();
  aaim[2400] = new ArbitraryAnaInputManager();
  aaim[2800] = new ArbitraryAnaInputManager();
  aaim[3200] = new ArbitraryAnaInputManager();
  aaim[3600] = new ArbitraryAnaInputManager();
  aaim[4000] = new ArbitraryAnaInputManager();
  
  aaim[1600]->SetInputTree("arbitrary1600/SimTree");
  aaim[2000]->SetInputTree("arbitrary2000/SimTree");
  aaim[2400]->SetInputTree("arbitrary2400/SimTree");
  aaim[2800]->SetInputTree("arbitrary2800/SimTree");
  aaim[3200]->SetInputTree("arbitrary3200/SimTree");
  aaim[3600]->SetInputTree("arbitrary3600/SimTree");
  aaim[4000]->SetInputTree("arbitrary4000/SimTree");
  
  for (auto& it: aaim) {
    it.second->SetInputFile(InFileName);
    it.second->LoadTree();
  }
  std::map<int, TEfficiency*> effEnergy;
  std::map<int, TProfile*>    nNoiseHit;
  std::map<int, TProfile*>    nSignalHit;
  std::map<int, TProfile*>    nNoiseCluster;
  std::map<int, TProfile*>    nSignalCluster;
  
  std::map<int, TH2D*>    th2d_nNoiseHit     ;
  std::map<int, TH2D*>    th2d_nSignalHit    ;
  std::map<int, TH2D*>    th2d_nNoiseCluster ;
  std::map<int, TH2D*>    th2d_nSignalCluster;
  
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
  
  for (auto const& it: aaim) {
    int fNEvent = nEvent;
    int adc_threshold = it.first;
    ArbitraryAnaInputManager* im = it.second;
    
    if (fNEvent!=-1) {
      fNEvent = std::min(fNEvent,(int)im->GetEntries());
    } else {
      fNEvent = im->GetEntries();
    }
    
    effEnergy     [adc_threshold] = new TEfficiency(Form("Eff_Th%i", adc_threshold),
                                                    ";electron KE [MeV];Efficiency", 25, 0, 25);
    nNoiseHit     [adc_threshold] = new TProfile(Form("nNoiseHit_Th%i", adc_threshold),
                                                 ";electron KE [MeV];nNoiseHit", 25, 0, 25);
    nSignalHit    [adc_threshold] = new TProfile(Form("nSignalHit_Th%i", adc_threshold),
                                                 ";electron KE [MeV];nSignalHit", 25, 0, 25);
    nNoiseCluster [adc_threshold] = new TProfile(Form("nNoiseCluster_Th%i", adc_threshold),
                                                 ";electron KE [MeV];nNoiseCluster", 25, 0, 25);
    nSignalCluster[adc_threshold] = new TProfile(Form("nSignalCluster_Th%i", adc_threshold),
                                                 ";electron KE [MeV];nSignalCluster", 25, 0, 25);
    th2d_nNoiseHit     [adc_threshold] = new TH2D(Form("nNoiseHit_Th%i_th2", adc_threshold),
                                                  Form("Threshold %i ADC;electron KE [MeV];nNoiseHit;Event",adc_threshold/100),
                                                  25, 0, 25, 20, 0, 20);
    th2d_nSignalHit    [adc_threshold] = new TH2D(Form("nSignalHit_Th%i_th2", adc_threshold),
                                                  Form("Threshold %i ADC;electron KE [MeV];nSignalHit;Event",adc_threshold/100),
                                                  25, 0, 25, 20, 0, 20);
    th2d_nNoiseCluster [adc_threshold] = new TH2D(Form("nNoiseCluster_Th%i_th2", adc_threshold),
                                                  Form("Threshold %i ADC;electron KE [MeV];nNoiseCluster;Event",adc_threshold/100),
                                                  25, 0, 25, 20, 0, 20);
    th2d_nSignalCluster[adc_threshold] = new TH2D(Form("nSignalCluster_Th%i_th2", adc_threshold),
                                                  Form("Threshold %i ADC;electron KE [MeV];nSignalCluster;Event",adc_threshold/100),
                                                  25, 0, 25, 20, 0, 20);
    
    for (int CurrentEvent=0; CurrentEvent<fNEvent; ++CurrentEvent) {
      im->GetEntry(CurrentEvent);
      std::vector<WireHit*> vec_WireHit;
      double nhit[2] = {0,0};
      for (int j=0; j<im->Hit_View->size(); ++j) {
        WireHit* hit = new WireHit((*im->Hit_View)[j],        (*im->Hit_True_GenType)[j],  (*im->Hit_Chan)[j],
                                   (*im->Hit_Time)[j],        (*im->Hit_SADC)[j],          (*im->Hit_RMS)[j],
                                   (*im->Hit_True_Energy)[j], (*im->Hit_True_EvEnergy)[j], (*im->Hit_True_MainTrID)[j],
                                   0.5, 0.5, 0.5,
                                   (*im->Hit_True_X)[j],      (*im->Hit_True_Y)[j],        (*im->Hit_True_Z)[j],
                                   0, (*im->Hit_True_nElec)[j]);
        nhit[(*im->Hit_True_GenType)[j]]++;
        vec_WireHit.push_back(hit);
      }
      bool selected = false;
      int ncluster = 0;
      int nnoisecluster = 0;
      std::vector<WireCluster*> vec_WireCluster;
      clusteng.ClusterHits2(vec_WireHit, vec_WireCluster);
      wiretrigger.SetIsSelected(vec_WireCluster);
      
      for (int c=0; c<vec_WireCluster.size(); ++c) {
        WireCluster* clust = vec_WireCluster[c];
        if (clust->GetIsSelected()) {
          if (clust->GetType()) {
            selected = true;
            ++ncluster;
          } else {
            ++nnoisecluster;
          }
        }
      }
      
      double VertStartX = im->True_Prim_VertStartX->at(0);
      double VertStartY = im->True_Prim_VertStartY->at(0);
      double VertStartZ = im->True_Prim_VertStartZ->at(0);
      double VertEndX   = im->True_Prim_VertEndX  ->at(0);
      double VertEndY   = im->True_Prim_VertEndY  ->at(0);
      double VertEndZ   = im->True_Prim_VertEndZ  ->at(0);
      double range = sqrt(( VertStartX - VertEndX ) * ( VertStartX - VertEndX ) +
                          ( VertStartY - VertEndY ) * ( VertStartY - VertEndY ) +
                          ( VertStartZ - VertEndZ ) * ( VertStartZ - VertEndZ ));
      
      double KE = im->True_Prim_Energy->at(0) * 1000. - 0.511;
      effEnergy     [adc_threshold]->Fill(selected     , KE);
      nNoiseHit     [adc_threshold]->Fill(KE, nhit[1]      );
      nSignalHit    [adc_threshold]->Fill(KE, nhit[0]      );
      nNoiseCluster [adc_threshold]->Fill(KE, nnoisecluster);
      nSignalCluster[adc_threshold]->Fill(KE, ncluster     );
      th2d_nNoiseHit     [adc_threshold]->Fill(KE, nhit[1]      );
      th2d_nSignalHit    [adc_threshold]->Fill(KE, nhit[0]      );
      th2d_nNoiseCluster [adc_threshold]->Fill(KE, nnoisecluster);
      th2d_nSignalCluster[adc_threshold]->Fill(KE, ncluster     );
      
      
      for (auto& it: vec_WireHit) {
        delete it;
        it = NULL;
      }
      vec_WireHit.clear();
    }
  }
  
  const std::vector<int> color = getColors();
  int i=0;
  TLegend* legend = new TLegend(0.1,0.1,0.9,0.9);
  for (auto const& it: effEnergy) {
    nNoiseHit     [it.first]->SetStats(0);
    nSignalHit    [it.first]->SetStats(0);
    nNoiseCluster [it.first]->SetStats(0);
    nSignalCluster[it.first]->SetStats(0);
    effEnergy     [it.first]->SetLineColor  (color.at(i));
    nNoiseHit     [it.first]->SetLineColor  (color.at(i));
    nSignalHit    [it.first]->SetLineColor  (color.at(i));
    nNoiseCluster [it.first]->SetLineColor  (color.at(i));
    nSignalCluster[it.first]->SetLineColor  (color.at(i));
    effEnergy     [it.first]->SetMarkerColor(color.at(i));
    nNoiseHit     [it.first]->SetMarkerColor(color.at(i));
    nSignalHit    [it.first]->SetMarkerColor(color.at(i));
    nNoiseCluster [it.first]->SetMarkerColor(color.at(i));
    nSignalCluster[it.first]->SetMarkerColor(color.at(i));
    th2d_nNoiseHit     [it.first]->SetStats(0);
    th2d_nSignalHit    [it.first]->SetStats(0);
    th2d_nNoiseCluster [it.first]->SetStats(0);
    th2d_nSignalCluster[it.first]->SetStats(0);
    legend->AddEntry(nNoiseCluster[it.first], Form("ADC Threshold = %i",it.first/100), "L");
    i++;
  }
  
  TCanvas c;
  c.Print((OutFileName+"[").c_str());
  PlotAll(effEnergy);
  c.Print(OutFileName.c_str());
  
  PlotAll(nNoiseHit);
  c.Print(OutFileName.c_str());
  
  PlotAll(nSignalHit);
  c.Print(OutFileName.c_str());
  
  PlotAll(nNoiseCluster);
  c.Print(OutFileName.c_str());
  
  PlotAll(nSignalCluster);
  c.Print(OutFileName.c_str());
  
  c.Clear();
  legend->Draw();
  c.Print(OutFileName.c_str());
  
  gPad->SetRightMargin(1.5*gPad->GetRightMargin());
  for (auto const& it: th2d_nNoiseHit     ) { it.second->Draw("COLZ"); c.Print(OutFileName.c_str()); }
  for (auto const& it: th2d_nSignalHit    ) { it.second->Draw("COLZ"); c.Print(OutFileName.c_str()); }
  for (auto const& it: th2d_nNoiseCluster ) { it.second->Draw("COLZ"); c.Print(OutFileName.c_str()); }
  for (auto const& it: th2d_nSignalCluster) { it.second->Draw("COLZ"); c.Print(OutFileName.c_str()); }
  
  c.Print((OutFileName+"]").c_str());
  
}


