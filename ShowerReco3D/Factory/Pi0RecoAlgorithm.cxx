#ifndef SHOWERRECO_PI0RECOALGORITHM_CXX
#define SHOWERRECO_PI0RECOALGORITHM_CXX

#include "Pi0RecoAlgorithm.h"

namespace showerreco {
  
  Pi0RecoAlgorithm::Pi0RecoAlgorithm()
  { 

    this->SetVerbose(true);

    ::showerreco::FilterPFPart* filterpfpart = new ::showerreco::FilterPFPart();
    filterpfpart->setMinNHitsAbsolute(5);
    filterpfpart->setMinNHitsLargest(10);
    filterpfpart->setVerbosity(false);
    this->AddAlgo(filterpfpart);
    
    ::showerreco::Angle3DFromVtxQweighted* angle3d = new ::showerreco::Angle3DFromVtxQweighted();
    this->AddAlgo(angle3d);

    ::showerreco::YPlaneStartPoint3D* start3d = new ::showerreco::YPlaneStartPoint3D();
    start3d->setVerbosity(false);
    this->AddAlgo(start3d);

    ::showerreco::LinearEnergy* energy = new ::showerreco::LinearEnergy();
    energy->setVerbosity(false);
    this->AddAlgo(energy);

    ::showerreco::FillLength* length = new ::showerreco::FillLength();
    length->setVerbosity(false);
    this->AddAlgo(length);

    this->PrintModuleList();
  }
  
}

#endif
