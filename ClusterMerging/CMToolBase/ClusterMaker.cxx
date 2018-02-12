#ifndef CLUSTER_CLUSTERMAKER_CXX
#define CLUSTER_CLUSTERMAKER_CXX

#include "ClusterMaker.h"

namespace cluster {

  ClusterMaker::ClusterMaker() {

    _vtx_w_cm = {0.,0.,0.};
    _vtx_t_cm = {0.,0.,0.};

    // get detector specific properties
    auto const* geom = ::lar::providerFrom<geo::Geometry>();
    auto const* detp = lar::providerFrom<detinfo::DetectorPropertiesService>();
    _wire2cm = geom->WirePitch(0,1,0);
    _time2cm = detp->SamplingRate() / 1000.0 * detp->DriftVelocity( detp->Efield(), detp->Temperature() );

  }

  void ClusterMaker::MakeClusters(::art::Event & e,
				  const std::string& fClusterProducer,
				  const std::string& fVertexProducer,
				  std::vector<::cluster::Cluster>& cluster) {

    // load input clusters
    art::Handle<std::vector<recob::Cluster> > clus_h;
    e.getByLabel(fClusterProducer,clus_h);

    // load associated hits
    art::FindManyP<recob::Hit> clus_hit_assn_v(clus_h, e, fClusterProducer);

    // load vertices
    art::Handle<std::vector<recob::Vertex> > vtx_h;
    e.getByLabel(fVertexProducer,vtx_h);

    if (loadVertex(vtx_h) == false) {
      std::cout << "NO VERTEX" << std::endl;
      return;
    }

    cluster.clear();
    cluster.reserve(clus_hit_assn_v.size());

    // loop trhough hits in each cluster
    for (size_t c=0; c < clus_hit_assn_v.size(); c++) {
      
      const std::vector<art::Ptr<recob::Hit> > hit_v = clus_hit_assn_v.at(c);
      
      // create new cluster
      ::cluster::Cluster clus;
      // and associated hits
      std::vector<::cluster::pt> pts;
      // fill hit information
      GetClusterPts(hit_v, pts);
      // assign hits to cluster
      clus.SetHits(pts);
      // and add to list of clusters
      cluster.push_back(clus);
      
    }// for all input clusters

    return;
  }

  void ClusterMaker::GetClusterPts(const std::vector<art::Ptr<recob::Hit> >& hit_v,
				   std::vector<::cluster::pt>& pt_v) {

    pt_v.clear();
    pt_v.reserve(hit_v.size());
    
    for (art::Ptr<recob::Hit> hit : hit_v) {
      
      auto pl = hit->WireID().Plane;

      float  q = hit->Integral();

      float hw = hit->WireID().Wire * _wire2cm;
      float ht = hit->PeakTime()    * _time2cm;
      
      float r = sqrt( (ht - _vtx_t_cm[pl]) * (ht - _vtx_t_cm[pl]) +
		      (hw - _vtx_w_cm[pl]) * (hw - _vtx_w_cm[pl]) );

      float cos_theta = (hw - _vtx_w_cm[pl]) / r;
      float sin_theta = (ht - _vtx_t_cm[pl]) / r;
      float theta_rad = atan2(sin_theta,cos_theta);
      float theta_deg = theta_rad * 180. / 3.1415;

      theta_deg += 180.;

      ::cluster::pt thispt(r,theta_deg,hw,ht,q,pl);

      thispt._idx = hit->key();
      
      pt_v.push_back(thispt);
      
    }// for all his in cluster

    return;
  }// create polar cluster


  

  bool ClusterMaker::loadVertex(const art::Handle<std::vector<::recob::Vertex> > vtx_h) {

    // load required services to obtain offsets
    auto const* detp = lar::providerFrom<detinfo::DetectorPropertiesService>();
    auto const& geomH = ::util::GeometryUtilities::GetME();
    
    if (vtx_h->size() != 1) return false;

    auto const& vtx = vtx_h->at(0);
    
    Double_t xyz[3] = {};
    vtx.XYZ(xyz);

    for (size_t pl = 0; pl < 3; pl++) {

      auto const& pt = geomH->Get2DPointProjection(xyz,pl);
      _vtx_w_cm[pl] = pt.w;
      _vtx_t_cm[pl] = pt.t + (detp->TriggerOffset() * _time2cm);

      }

    return true;
  }
  

}

#endif
