////////////////////////////////////////////////////////////////////////
// Class:       ShowerReco3D
// Plugin Type: producer (art v2_09_06)
// File:        ShowerReco3D_module.cc
//
// Generated at Fri Feb  9 16:38:52 2018 by David Caratelli using cetskelgen
// from cetlib version v3_01_03.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "ProtoShower/ProtoShowerAlgBase.h"

#include <memory>

class ShowerReco3D;


class ShowerReco3D : public art::EDProducer {
public:
  explicit ShowerReco3D(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  ShowerReco3D(ShowerReco3D const &) = delete;
  ShowerReco3D(ShowerReco3D &&) = delete;
  ShowerReco3D & operator = (ShowerReco3D const &) = delete;
  ShowerReco3D & operator = (ShowerReco3D &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.

};


ShowerReco3D::ShowerReco3D(fhicl::ParameterSet const & p)
// :
// Initialize member data here.
{
  // Call appropriate produces<>() functions here.
}

void ShowerReco3D::produce(art::Event & e)
{
  // Implementation of required member function here.
}

void ShowerReco3D::beginJob()
{
  // Implementation of optional member function here.
}

void ShowerReco3D::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(ShowerReco3D)
