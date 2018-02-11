#ifndef PROTOSHOWER_PROTOSHOWERALGBASE_CXX
#define PROTOSHOWER_PROTOSHOWERALGBASE_CXX

#include "ProtoShowerCMTool.h"

namespace protoshower {

  void ProtoShowerCMTool::GenerateProtoShowers(::art::Event & e,
					       const std::string& fPFPproducer,
					       std::vector<protoshower::ProtoShower> & proto_shower_v) {

    std::cout << "grabbing PFParticle" << std::endl;

    // grab PFParticles in event
    art::Handle<std::vector<recob::PFParticle> > pfp_h;
    e.getByLabel(fPFPproducer,pfp_h);    

    std::cout << "grabbing associated Clusters" << std::endl;
    
    // grab clusters associated with PFParticles
    art::FindManyP<recob::Cluster> pfp_clus_assn_v(pfp_h, e, fPFPproducer);

    std::cout << "grabbing associated Hits" << std::endl;

    // grab the hits associated to the PFParticles
    auto pfp_hit_assn_v = lar::FindManyInChainP<recob::Hit, recob::Cluster>::find(pfp_h, e, fPFPproducer);

    std::cout << "grabbing Vertex" << std::endl;

    // load event vertex associated to tagged neutrino interaction
    art::Handle<std::vector<recob::Vertex> > vertex_h;
    e.getByLabel(fPFPproducer,vertex_h);

    // loop through PFParticles
    for (size_t p=0; p < pfp_h->size(); p++) {

      ::protoshower::ProtoShower proto_shower;
      proto_shower.Reset();

      const recob::PFParticle pfp = pfp_h->at(p);
      
      // associated clusters
      const std::vector< art::Ptr<recob::Cluster> >& clus_v = pfp_clus_assn_v.at(p);
      
      // associated hits
      const std::vector< art::Ptr<recob::Hit> >& hit_v = pfp_hit_assn_v.at(p);

      std::cout << "there are " << hit_v.size() << " hits associated to the shower" << std::endl;

      // loop through clusters
      for (size_t c=0; c < clus_v.size(); c++) {

	auto const clus = clus_v.at(c);

	// find the hits associated to this cluster
	std::vector< art::Ptr<recob::Hit> > clusterhits;
	for (auto const& hit : hit_v) {
	  // if hit in same plane as cluster -> add to vector
	  if (hit->WireID().Plane == clus->Plane().Plane )
	    clusterhits.push_back( hit );
	}//for all hits associated to PFParticle

	proto_shower._clusters.at(c) = MakeCluster2D( clus, clusterhits );

	proto_shower.hasCluster2D(true);
	
      }// for all clusters
      
      // require a single vertex!
      if (vertex_h->size() == 1) {
	
	auto const vtx = vertex_h->at(0);
	Double_t xyz[3] = {};
	vtx.XYZ(xyz);
	proto_shower._vertex = xyz;//TVector3(xyz[0],xyz[1],xyz[2]);
	proto_shower.hasVertex(true);
	
      }// if there is only one vertex
      
      proto_shower_v.push_back( proto_shower );
      
    }// for all PFParticles
    
    std::cout << "there are " << pfp_clus_assn_v.size() << " associations" << std::endl;
    std::cout << "there are " << vertex_h->size() << " vertices" << std::endl;
    
    
  }// GenerateProtoShower end

}// namespace

#endif
