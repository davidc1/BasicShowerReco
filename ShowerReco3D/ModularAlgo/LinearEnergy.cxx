#ifndef LINEARENERGY_CXX
#define LINEARENERGY_CXX

#include <iomanip>

#include "LinearEnergy.h"

namespace showerreco {
  
  LinearEnergy::LinearEnergy()
  {
    
    _name = "LinearEnergy";

    return;
  }
  
  void LinearEnergy::initialize()
  {
    return;
  }
  
  void LinearEnergy::do_reconstruction(const ::protoshower::ProtoShower & proto_shower,
				       Shower_t& resultShower) {
    
    //if the module does not have 2D cluster info -> fail the reconstruction
    if (!proto_shower.hasCluster2D()){
      std::stringstream ss;
      ss << "Fail @ algo " << this->name() << " due to missing 2D cluster";
      throw ShowerRecoException(ss.str());
    }
    
    auto & clusters = proto_shower.clusters();
    
    // This function takes the shower cluster set and calculates an energy in MeV for each plane
    
    // check if plane 2 has been used.
    // if so, we will fill the global energy with that from plane 2
    // otherwise, average the other two planes
    bool hasPl2 = false;
    
    // we want an energy for each plane
    for (size_t n = 0; n < clusters.size(); n++) {

      auto const& clus = clusters.at(n);
      
      // get the hits associated with this cluster
      auto const& hits = clus._hits;
      
      // get the plane associated with this cluster
      auto const& pl = clus._plane;
      
      if (pl == 2)
	hasPl2 = true;
      
      // store calculated energy
      double E  = 0.;
      
      // loop over hits
      for (auto const &h : hits) 
	E += h.charge;
      
      if (_verbose)
	std::cout << "energy on plane " << pl << " is : " << E << std::endl;
      
      // set the energy for this plane
      resultShower.fTotalEnergy_v[pl] = E;

      
    }// for all input clusters
    
    if (hasPl2)
      resultShower.fTotalEnergy = resultShower.fTotalEnergy_v[2];
    else
      resultShower.fTotalEnergy = ( resultShower.fTotalEnergy_v[0] + resultShower.fTotalEnergy_v[1] ) / 2.;
    
    return;
  }
  
} //showerreco

#endif
