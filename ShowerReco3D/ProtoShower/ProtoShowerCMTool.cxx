#ifndef PROTOSHOWER_PROTOSHOWERALGBASE_CXX
#define PROTOSHOWER_PROTOSHOWERALGBASE_CXX

#include "ProtoShowerCMTool.h"

namespace protoshower {

  ProtoShowerCMTool::GenerateProtoShower(::art::Event & e,
					 const art::Handle<std::vector<recob::PFParticle> >& pfp_v,
					 const size_t proto_shower_pfpart,
					 protoshower::ProtoShower & proto_shower) {

    // grab clusters associated with PFParticles
    art::FindMany<recob::Cluster> pfp_clus_assn_v(pfp_h, e, fClusProducer);

    // load event vertex associated to tagged neutrino interaction
    art::Handle<std::vector<recob::Vertex> > vertex_h;
    e.getByLabel(fVertexProducer,vertex_h);

    std::cout << "there are " << pfp_clus_assn_v.size() << " associations" << std::endl;
    std::cout << "there are " << vertex_h->size() << " vertices" << std::endl;


  }// GenerateProtoShower end

}// namespace

#endif
