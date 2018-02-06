#ifndef FILLLENGTH_CXX
#define FILLLENGTH_CXX

#include "FillLength.h"

namespace showerreco {

  FillLength::FillLength()
  {
    _name        = "FillLength";
  }

void FillLength::do_reconstruction(const ::protoshower::ProtoShower & proto_shower,
                                  Shower_t& resultShower) {
  
  //if the module does not have 2D cluster info -> fail the reconstruction
  if (!proto_shower.hasCluster2D()){
    std::stringstream ss;
      ss << "Fail @ algo " << this->name() << " due to missing 2D cluster";
      throw ShowerRecoException(ss.str());
  }
  
  auto & clusters = proto_shower.clusters();

  // use collection-plane cluster to compute length

  bool filled = false;
  
  for (auto const& clus : clusters) {

    if (clus._plane != 2) continue;
    
    auto sW  = clus._start.w;
    auto eW  = clus._end.w;
    
    // length on Y-plane projection
    auto dLY = fabs(eW - sW);
    
    resultShower.fLength = dLY / fabs(resultShower.fDCosStart[2]) ;
    
    resultShower.fOpeningAngle = clus._opening_angle;
    
    filled = true;
    
  }// for all clusters    
  
  // if collection-plane info not found -> skip
  if (filled == false) {
    std::stringstream ss;
    ss << "Fail @ algo " << this->name() << " due to missing collection-plane cluster";
    throw ShowerRecoException(ss.str());
  }
  
  return;
}

}
  
#endif
  
