#ifndef LINEARENERGY_CXX
#define LINEARENERGY_CXX

#include <iomanip>
#include <iostream>

#include "uboone/BasicShowerReco/ShowerReco3D/Base/ShowerRecoModuleBase.h"
#include "uboone/BasicShowerReco/ShowerReco3D/Base/Calorimetry.h"

#include "uboone/Database/TPCEnergyCalib/TPCEnergyCalibService.h"
#include "uboone/Database/TPCEnergyCalib/TPCEnergyCalibProvider.h"

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

    float ChargeCorrection(const double& q, const double& w, const double& t, const TVector3& dir, const TVector3& vtx,
			   const int& pl, const lariov::TPCEnergyCalibProvider& energyCalibProvider);

    //double _recomb, _ADC_to_e, _e_to_MeV;

  };
  
  LinearEnergy::LinearEnergy(const fhicl::ParameterSet& pset)
  {
    
    configure(pset);
    _name = "LinearEnergy";
    return;
  }

  void LinearEnergy::configure(const fhicl::ParameterSet& pset)
  {
    //_recomb    = pset.get<float>("recomb");
    //_ADC_to_e  = pset.get<float>("ADC_to_e");
    _verbose   = pset.get<bool>("verbose",false);
    //_e_to_MeV  = 23.6 * (1e-6);
    return;
  }
  
  void LinearEnergy::initialize()
  {
    return;
  }
  
  void LinearEnergy::do_reconstruction(const ::protoshower::ProtoShower & proto_shower,
				       Shower_t& resultShower) {
    

    //handle to tpc energy calibration provider
    const lariov::TPCEnergyCalibProvider& energyCalibProvider  = art::ServiceHandle<lariov::TPCEnergyCalibService>()->GetProvider();
    

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
      for (auto const &h : hits) {
	auto relcalib = ChargeCorrection(h.charge, h.w, h.t, resultShower.fDCosStart, resultShower.fXYZStart, pl, energyCalibProvider);
	E += h.charge * relcalib;
      }
      
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

  float LinearEnergy::ChargeCorrection(const double& q, const double& w, const double& t, const TVector3& dir, const TVector3& vtx,
				       const int& pl, const lariov::TPCEnergyCalibProvider& energyCalibProvider){

    // find 3D position of hit                                                                                                        
    double z = w;
    double x = t;

    // get 2D distance of hit to vtx                                                                                                  
    double r2D = sqrt( ( (z-vtx.Z()) * (z-vtx.Z()) ) + ( (x-vtx.X()) * (x-vtx.X()) ) );
    double r3D = r2D/dir[1];

    auto xyz = vtx + dir * r3D;

    float yzcorrection = energyCalibProvider.YZdqdxCorrection(pl, xyz.Y(), xyz.Z());
    float xcorrection  = energyCalibProvider.XdqdxCorrection(pl,  xyz.X());

    if (!yzcorrection) yzcorrection = 1.;
    if (!xcorrection ) xcorrection  = 1.;

    return yzcorrection * xcorrection;
  }

  DEFINE_ART_CLASS_TOOL(LinearEnergy)
} //showerreco

#endif
