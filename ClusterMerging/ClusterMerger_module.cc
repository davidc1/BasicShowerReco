////////////////////////////////////////////////////////////////////////
// Class:       ClusterMerger
// Plugin Type: producer (art v2_09_06)
// File:        ClusterMerger_module.cc
//
// Generated at Mon Feb 12 21:28:38 2018 by David Caratelli using cetskelgen
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

#include "uboone/BasicShowerReco/ClusterMerging/CMToolApp/CMergeHelper.h"
#include "uboone/BasicShowerReco/ClusterMerging/CMToolBase/ClusterMaker.h"

#include <memory>

class ClusterMerger;


class ClusterMerger : public art::EDProducer {
public:
  explicit ClusterMerger(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  ClusterMerger(ClusterMerger const &) = delete;
  ClusterMerger(ClusterMerger &&) = delete;
  ClusterMerger & operator = (ClusterMerger const &) = delete;
  ClusterMerger & operator = (ClusterMerger &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  void FillClusterProperties(::recob::Cluster& cluster,
			     const ::cluster::Cluster& CMCluster);

  // Declare member data here.

  float _wire2cm, _time2cm;

  // ClusterMaker class
  ::cluster::ClusterMaker* _CMaker;
  
  // cluster helper
  ::cmtool::CMergeHelper* _merge_helper;

  std::string fClusterProducer, fVertexProducer;

};


ClusterMerger::ClusterMerger(fhicl::ParameterSet const & p)
// :
// Initialize member data here.
{

  // get detector specific properties
  auto const* geom = ::lar::providerFrom<geo::Geometry>();
  auto const* detp = lar::providerFrom<detinfo::DetectorPropertiesService>();
  _wire2cm = geom->WirePitch(0,1,0);
  _time2cm = detp->SamplingRate() / 1000.0 * detp->DriftVelocity( detp->Efield(), detp->Temperature() );

  fClusterProducer = p.get<std::string>("ClusterProducer");
  fVertexProducer  = p.get<std::string>("VertexProducer" );

  produces<std::vector<recob::Cluster> >();
  produces<art::Assns <recob::Cluster, recob::Hit> >();

}

void ClusterMerger::produce(art::Event & e)
{
  // Implementation of required member function here.

  std::unique_ptr< std::vector<recob::Cluster> > Cluster_v(new std::vector<recob::Cluster>);
  std::unique_ptr< art::Assns <recob::Cluster, recob::Hit> > Cluster_Hit_assn_v(new art::Assns<recob::Cluster,recob::Hit>);

  // load data products needed

  // load input clusters
  art::Handle<std::vector<recob::Cluster> > clus_h;
  e.getByLabel(fClusterProducer,clus_h);
  
  // load associated hits
  art::FindManyP<recob::Hit> clus_hit_assn_v(clus_h, e, fClusterProducer);

  // load vertices
  art::Handle<std::vector<recob::Vertex> > vtx_h;
  e.getByLabel(fVertexProducer,vtx_h);
  

  // create cluster::Clusters
  std::vector<::cluster::Cluster> event_clusters;
  _CMaker->MakeClusters(clus_h, clus_hit_assn_v, vtx_h, event_clusters);

  _merge_helper->Process(event_clusters);

  // get output clusters
  auto out_clusters = _merge_helper->GetClusters();

  std::cout << "there are " << out_clusters.size() << " outputs" << std::endl;

  // loop through output clusters
  for (auto const& out_clus : out_clusters) {

    ::recob::Cluster outclus;
    FillClusterProperties(outclus,out_clus);
    
    Cluster_v->emplace_back(outclus);
    
  }// for all output clsters
  
  e.put(std::move(Cluster_v));
  e.put(std::move(Cluster_Hit_assn_v));

}

void ClusterMerger::beginJob()
{
  // Implementation of optional member function here.
}

void ClusterMerger::endJob()
{
  // Implementation of optional member function here.
}

void ClusterMerger::FillClusterProperties(::recob::Cluster& cluster,
					  const ::cluster::Cluster& CMCluster) {
  
  float startW = CMCluster._start_pt._w;
  float startT = CMCluster._start_pt._t;
  
  float endW   = CMCluster._end_pt._w;
  float endT   = CMCluster._end_pt._t;
  
  cluster.set_planeID(geo::PlaneID(0,0,CMCluster._plane));
  
  cluster.set_start_wire( int(startW / _wire2cm), 0.5);
  cluster.set_start_tick( int(startT / _time2cm), 0.5);
  
  cluster.set_end_wire( int(endW / _wire2cm), 0.5);
  cluster.set_end_tick( int(endT / _time2cm), 0.5);
  
  cluster.set_start_angle( CMCluster._angle );
  cluster.set_summedADC( CMCluster._sum_charge, 0., 0.);
  cluster.set_n_hits( CMCluster.size() );

  recob::Cluster clus(startW, 0., startT, 0., 0., CMCluster._angle, 0., 
		      endW,   0., endT,   0., 0., 0., 0., 
		      CMCluster._sum_charge, 0., CMCluster._sum_charge, 0., 
		      CMCluster.size(), 0., 0., n,
		      hits->at(clus_idx_v[0]).View(),
		      geo::PlaneID(0,0,CMCluster._plane));
  
}


DEFINE_ART_MODULE(ClusterMerger)
