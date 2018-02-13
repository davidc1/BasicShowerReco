////////////////////////////////////////////////////////////////////////
// Class:       ClusterMatcher
// Plugin Type: producer (art v2_09_06)
// File:        ClusterMatcher_module.cc
//
// Generated at Mon Feb 12 21:28:35 2018 by David Caratelli using cetskelgen
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

#include <memory>

class ClusterMatcher;


class ClusterMatcher : public art::EDProducer {
public:
  explicit ClusterMatcher(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  ClusterMatcher(ClusterMatcher const &) = delete;
  ClusterMatcher(ClusterMatcher &&) = delete;
  ClusterMatcher & operator = (ClusterMatcher const &) = delete;
  ClusterMatcher & operator = (ClusterMatcher &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.

};


ClusterMatcher::ClusterMatcher(fhicl::ParameterSet const & p)
// :
// Initialize member data here.
{
  // Call appropriate produces<>() functions here.
}

void ClusterMatcher::produce(art::Event & e)
{
  // Implementation of required member function here.
}

void ClusterMatcher::beginJob()
{
  // Implementation of optional member function here.
}

void ClusterMatcher::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(ClusterMatcher)
