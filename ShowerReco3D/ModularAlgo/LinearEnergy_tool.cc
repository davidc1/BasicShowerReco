#ifndef LINEARENERGY_CXX
#define LINEARENERGY_CXX

#include <iomanip>
#include <iostream>

#include "uboone/BasicShowerReco/ShowerReco3D/Base/ShowerRecoModuleBase.h"
#include "uboone/BasicShowerReco/ShowerReco3D/Base/Calorimetry.h"

#include "TTree.h"
#include "art/Framework/Services/Optional/TFileService.h"

//#include "uboone/Database/TPCEnergyCalib/TPCEnergyCalibService.h"
//#include "uboone/Database/TPCEnergyCalib/TPCEnergyCalibProvider.h"

/**
   \class ShowerRecoModuleBase
   User defined class ShowerRecoModuleBase ... these comments are used to generate
   doxygen documentation!
 */

namespace showerreco {

  class LinearEnergy : public ShowerRecoModuleBase {
    
  public:
    
    /// Default constructor
    LinearEnergy(const fhicl::ParameterSet& pset);
    
    /// Default destructor
    ~LinearEnergy(){};

    void configure(const fhicl::ParameterSet& pset);

    
    void do_reconstruction(const ::protoshower::ProtoShower &, Shower_t &);
    
    void initialize();
    
  private:

    TTree* _energy_tree;
    double _e0, _e1, _e2;
    int    _nhit0, _nhit1, _nhit2;

  };
  
  LinearEnergy::LinearEnergy(const fhicl::ParameterSet& pset)
  {
    
    configure(pset);
    _name = "LinearEnergy";

    art::ServiceHandle<art::TFileService> tfs;
    _energy_tree = tfs->make<TTree>("_energy_tree","Energy TTree");
    _energy_tree->Branch("_e0",&_e0,"e0/D");
    _energy_tree->Branch("_e1",&_e1,"e1/D");
    _energy_tree->Branch("_e2",&_e2,"e2/D");
    _energy_tree->Branch("_nhit0",&_nhit0,"nhit0/I");
    _energy_tree->Branch("_nhit1",&_nhit1,"nhit1/I");
    _energy_tree->Branch("_nhit2",&_nhit2,"nhit2/I");

    return;
  }

  void LinearEnergy::configure(const fhicl::ParameterSet& pset)
  {
    _verbose   = pset.get<bool>("verbose",false);

    return;
  }
  
  void LinearEnergy::initialize()
  {
    return;
  }
  
  void LinearEnergy::do_reconstruction(const ::protoshower::ProtoShower & proto_shower,
				       Shower_t& resultShower) {
    

    _e0 = 0.;
    _e1 = 0.;
    _e2 = 0.;
    _nhit0 = 0;
    _nhit1 = 0;
    _nhit2 = 0;

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
      
      // store calculated energy [MeV]
      double E  = 0.;
      
      // loop over hits
      for (auto const &h : hits) 
	E += h.charge;// * _ADC_to_e * _e_to_MeV / _recomb;
      
      if (_verbose)
	std::cout << "energy on plane " << pl << " is : " << E << std::endl;
      
      // set the energy for this plane
      resultShower.fTotalEnergy_v[pl] = E;

      if (pl==0) {
	_nhit0 = hits.size();
	_e0     = E;
      }
      if (pl==1) {
	_nhit1 = hits.size();
	_e1     = E;
      }
      if (pl==2) {
	_nhit2 = hits.size();
	_e2     = E;
      }
      
    }// for all input clusters

    _energy_tree->Fill();
    
    if (hasPl2)
      resultShower.fTotalEnergy = resultShower.fTotalEnergy_v[2];
    else
      resultShower.fTotalEnergy = ( resultShower.fTotalEnergy_v[0] + resultShower.fTotalEnergy_v[1] ) / 2.;
    
    return;
  }

  DEFINE_ART_CLASS_TOOL(LinearEnergy)
} //showerreco

#endif
