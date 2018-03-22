////////////////////////////////////////////////////////////////////////
// Class:       PhotonMerge
// Plugin Type: producer (art v2_09_06)
// File:        PhotonMerge_module.cc
//
// Generated at Wed Mar 21 07:50:29 2018 by David Caratelli using cetskelgen
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

// larsoft data-products
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardata/Utilities/AssociationUtil.h"

#include "lardata/Utilities/FindManyInChainP.h"

class PhotonMerge;


class PhotonMerge : public art::EDProducer {
public:
  explicit PhotonMerge(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PhotonMerge(PhotonMerge const &) = delete;
  PhotonMerge(PhotonMerge &&) = delete;
  PhotonMerge & operator = (PhotonMerge const &) = delete;
  PhotonMerge & operator = (PhotonMerge &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.
  std::string fShrProducer, fVtxProducer, fPhotonProducer;

};


PhotonMerge::PhotonMerge(fhicl::ParameterSet const & p)
// :
// Initialize member data here.
{
  produces<std::vector<recob::Shower> >();
  produces<std::vector<recob::PFParticle> >();
  produces<std::vector<recob::Cluster> >();
  produces<art::Assns <recob::Shower, recob::PFParticle> >();
  produces<art::Assns <recob::Shower, recob::Cluster>    >();
  produces<art::Assns <recob::Shower, recob::Hit>        >();
}

void PhotonMerge::produce(art::Event & e)
{

  // produce recob::Showers
  std::unique_ptr< std::vector<recob::Shower>     > Shower_v    (new std::vector<recob::Shower>     );
  std::unique_ptr< std::vector<recob::PFParticle> > PFParticle_v(new std::vector<recob::PFParticle> );
  std::unique_ptr< std::vector<recob::Cluster>    > Cluster_v   (new std::vector<recob::Cluster>    );
  std::unique_ptr< art::Assns <recob::Shower, recob::PFParticle> > Shower_PFP_assn_v    ( new art::Assns<recob::Shower, recob::PFParticle>);
  std::unique_ptr< art::Assns <recob::Shower, recob::Cluster>    > Shower_Cluster_assn_v( new art::Assns<recob::Shower, recob::Cluster>   );
  std::unique_ptr< art::Assns <recob::Shower, recob::Hit>        > Shower_Hit_assn_v    ( new art::Assns<recob::Shower, recob::Hit>       );


  // load input showers
  auto const& shr_h = e.getValidHandle<std::vector<recob::Shower>>(fShrProducer);
  // load input vertices
  auto const& vtx_h = e.getValidHandle<std::vector<recob::Vertex>>(fVtxProducer);
  // load input photon clusters
  auto const& photon_h = e.getValidHandle<std::vector<recob::Cluster>>(fPhotonProducer);

  // grab clusters associated with shower
  art::FindManyP<recob::Cluster> shr_clus_assn_v(shr_h, e, fShrProducer);
  // grab the hits associated to the showers
  auto shr_hit_assn_v = lar::FindManyInChainP<recob::Hit, recob::Cluster>::find(shr_h, e, fShrProducer);

  if (vtx_h->size() != 1)
    std::cout << "no vertex" << std::endl;
  if (photon_h->size() == 0)
    std::cout << "no photons" << std::endl;

  e.put(std::move(Shower_v));
  e.put(std::move(PFParticle_v));
  e.put(std::move(Cluster_v));
  e.put(std::move(Shower_PFP_assn_v));
  e.put(std::move(Shower_Cluster_assn_v));
  e.put(std::move(Shower_Hit_assn_v));

}

void PhotonMerge::beginJob()
{
  // Implementation of optional member function here.
}

void PhotonMerge::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(PhotonMerge)
