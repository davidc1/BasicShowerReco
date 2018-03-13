////////////////////////////////////////////////////////////////////////
// Class:       TruthVertex
// Plugin Type: filter (art v2_09_06)
// File:        TruthVertex_module.cc
//
// Generated at Mon Mar  5 14:14:25 2018 by David Caratelli using cetskelgen
// from cetlib version v3_01_03.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <memory>

#include "art/Framework/Services/Optional/TFileService.h"

#include "lardataobj/RecoBase/Vertex.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardata/Utilities/AssociationUtil.h"

class TruthVertex;


class TruthVertex : public art::EDFilter {
public:
  explicit TruthVertex(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  TruthVertex(TruthVertex const &) = delete;
  TruthVertex(TruthVertex &&) = delete;
  TruthVertex & operator = (TruthVertex const &) = delete;
  TruthVertex & operator = (TruthVertex &&) = delete;

  // Required functions.
  bool filter(art::Event & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.

  // producers
  std::string fMCTruthProducer;

};


TruthVertex::TruthVertex(fhicl::ParameterSet const & p)
// :
// Initialize member data here.
{
  produces< std::vector< recob::Vertex > >();
  fMCTruthProducer = p.get<std::string>("MCTruthProducer");
}

bool TruthVertex::filter(art::Event & e)
{

  // produce vertex
  std::unique_ptr< std::vector<recob::Vertex> > Vtx_v(new std::vector<recob::Vertex>);

  Double_t xyz[3] = {};

  auto const& mct_h = e.getValidHandle<std::vector<simb::MCTruth> >(fMCTruthProducer);
  auto mct = mct_h->at(0);
  size_t npart = mct.NParticles();

  bool foundPi0 = false;

  for (size_t i=0; i < npart; i++){
    auto const& part = mct.GetParticle(i);
    if ( (part.PdgCode() == 111) and (part.StatusCode() == 1) ){
      xyz[0] = part.Trajectory().X(0);
      xyz[1] = part.Trajectory().Y(0);
      xyz[2] = part.Trajectory().Z(0);
      foundPi0 = true;
      break;
    }
  }
  
  
  if (foundPi0 == false) {
    e.put(std::move(Vtx_v));
    return false;
  }
  
  recob::Vertex vtx(xyz);
  Vtx_v->emplace_back(vtx);
  
  e.put(std::move(Vtx_v));

  return true;
}

void TruthVertex::beginJob()
{
  // Implementation of optional member function here.
}

void TruthVertex::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(TruthVertex)
