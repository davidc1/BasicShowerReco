/**
 * \file ProtoShowerAlgBase.h
 *
 * \ingroup ProtoShower
 *
 * \brief Class def header for a class ProtoShowerAlgBase
 *
 * @author david caratelli
 */

/** \addtogroup ProtoShower

    @{*/
#ifndef PROTOSHOWERALGBASE_H
#define PROTOSHOWERALGBASE_H

#include <iostream>

#include "art/Framework/Principal/Event.h"

#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/Utilities/FindManyInChainP.h"

#include "ProtoShower.h"

/**
   \class ProtoShowerAlgBase
   User defined class ProtoShowerAlgBase ... these comments are used to generate
   doxygen documentation!
 */

namespace protoshower {

class ProtoShowerAlgBase {

public:

  /// Default constructor
  ProtoShowerAlgBase();

  /// Default destructor
  virtual ~ProtoShowerAlgBase() {}

  virtual void GenerateProtoShowers(::art::Event & e,
				    const std::string& fPFPproducer,
				    const size_t proto_shower_pfpart,
				    std::vector<protoshower::ProtoShower> & proto_shower_v) = 0;

  /**
     @brief function which takes recob::Cluster and vector of recob::Hits to create cluster2d::Cluster2D object
     @input art::Ptr to cluster
     @input vector of art::Ptr to hits associated to the cluster
  */
  cluster2d::Cluster2D MakeCluster2D( const art::Ptr<recob::Cluster>& clus, const std::vector< art::Ptr<recob::Hit> >& hit_v);


  std::string name() { return _name; }

protected:

  std::string _name;

  double _wire2cm, _time2cm;

};

}// namespace

#endif
/** @} */ // end of doxygen group

