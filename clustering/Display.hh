#ifndef DISPLAY_HH
#define DISPLAY_HH


#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"

#include "Helper.h"
#include "InputManager.hh"


class Display {
public:
  Display(const std::string F, const std::string T):FileName(F), TreeName(T){
    c = new TCanvas();
    im.SetInputFile(FileName.c_str());
    im.SetInputTree(TreeName.c_str());
    im.LoadTree();

  };
  virtual void DisplayEvent(const int, const int)=0;
  ~Display(){
    delete c;
    c = NULL;
  };
protected:
  TCanvas* c;
  SNAnaInputManager im;
  std::string FileName;
  std::string TreeName;
};


#endif
