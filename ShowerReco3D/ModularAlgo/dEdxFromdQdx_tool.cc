#ifndef DEDXFROMDQDX_CXX
#define DEDXFROMDQDX_CXX

#include <iostream>
#include "uboone/BasicShowerReco/ShowerReco3D/Base/ShowerRecoModuleBase.h"
/**
   \class dEdxFromdQdx : ShowerRecoModuleBase
   This is meant to compute the 2D dQdx along the start of the shower.
 */

namespace showerreco {

class dEdxFromdQdx : public ShowerRecoModuleBase {

public:

  /// Default constructor
  dEdxFromdQdx(const fhicl::ParameterSet& pset);

  /// Default destructor
  ~dEdxFromdQdx(){};

  void configure(const fhicl::ParameterSet& pset);

  void do_reconstruction(const ::protoshower::ProtoShower &, Shower_t &);

  void initialize();

private:

  bool _use_pitch;
  double _dEdx;
  int _pl_best;
  int _pl;
};

  dEdxFromdQdx::dEdxFromdQdx(const fhicl::ParameterSet& pset)
  {
    _name = "dEdxFromdQdx";
    configure(pset);
  }

  void dEdxFromdQdx::configure(const fhicl::ParameterSet& pset)
  {
    _use_pitch = pset.get<bool>("use_pitch");
    _verbose   = pset.get<bool>("verbose",false);
    return;
  }

void dEdxFromdQdx::initialize()
{
  if (_tree) delete _tree;
  _tree = new TTree(_name.c_str(), "dQdx Info Tree");
  _tree->Branch("_dEdx", &_dEdx, "dEdx/D");
  _tree->Branch("_pl", &_pl, "pl/I");
  _tree->Branch("_pl_best", &_pl_best, "pl_best/I");
  return;
}

void dEdxFromdQdx::do_reconstruction(const ::protoshower::ProtoShower & proto_shower,
                                     Shower_t& resultShower) {

  auto & clusters = proto_shower.clusters();


  // // get the 3D direction reconstructed hopefully in a previous step
  // auto const& dir3D = resultShower.fDCosStart;

  // get the vector of dQdx filled by someone else
  auto const& dQdx_v = resultShower.fdQdx_v;

  // loop over all input cluster -> calculate a dQdx per plane
  for (size_t n = 0; n < clusters.size(); n++) {

    // get the plane associated with this cluster
    auto const& pl = clusters.at(n)._plane;

    // get the best plane
    auto const& pl_best = resultShower.fBestdQdxPlane;

    double dedx = 0.;

    double dqdx = dQdx_v[pl];

    if (_verbose) { std::cout << "dQdx on plane : " << pl << " -> " << dqdx << std::endl; }

    dedx = dqdx;//larutil::LArProperties::GetME()->ModBoxCorrection(dqdx);

    if (_verbose) { std::cout << "dEdx on plane : " << pl << " -> " << dedx << std::endl; }

    // take the dQdx measured on each plane and convert
    // to dEdx using a recombination model

    _dEdx = dedx;
    _pl = pl;
    _pl_best = pl_best;
    _tree->Fill();
    resultShower.fdEdx_v[pl] = dedx;

  }

  return;
}

  DEFINE_ART_CLASS_TOOL(dEdxFromdQdx)
} //showerreco

#endif
