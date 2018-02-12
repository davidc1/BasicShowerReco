/**
 * \file ClusterMaker.h
 *
 * \ingroup CMToolBase
 * 
 * \brief Class def header for a class Cluster
 *
 * @author david caratelli
 */

/** \addtogroup CMToolBase

    @{*/
#ifndef CLUSTER_CLUSTERMAKER_H
#define CLUSTER_CLUSTERMAKER_H

#include <iostream>
#include <vector>

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardata/Utilities/GeometryUtilities.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "art/Framework/Principal/Event.h"

// Data Products
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardata/Utilities/AssociationUtil.h"

#include "Cluster.h"

/**
   \class Cluster
   User defined class Cluster ... these comments are used to generate
   doxygen documentation!
 */

namespace cluster {


  
  class ClusterMaker {
    
  public:
    
    /// Default constructor
    ClusterMaker();
    
    /// Default destructor
    ~ClusterMaker(){}

    void MakeClusters(::art::Event & e,
		      const std::string& fClusterProducer,
		      const std::string& fVertexProducer,
		      std::vector<::cluster::Cluster>& cluster);

  private:

    bool loadVertex(const art::Handle<std::vector<::recob::Vertex> > vtx_h);

    void GetClusterPts(const std::vector<art::Ptr<recob::Hit> >& hit_v,
		       std::vector<::cluster::pt>& pt_v);

    /// vertex coordinates
    std::vector<double> _vtx_w_cm;
    std::vector<double> _vtx_t_cm;
    
    /// conversion factors for hits
    double _wire2cm, _time2cm;

  };

}
#endif
/** @} */ // end of doxygen group 

