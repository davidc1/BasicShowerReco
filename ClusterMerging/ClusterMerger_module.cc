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

#include "art/Persistency/Common/PtrMaker.h"

#include "art/Utilities/make_tool.h"

#include "uboone/BasicShowerReco/ClusterMerging/CMToolApp/CMergeHelper.h"
#include "uboone/BasicShowerReco/ClusterMerging/CMToolBase/ClusterMaker.h"

#include <memory>

class ClusterMerger;


class ClusterMerger : public art::EDProducer {
public:
  explicit ClusterMerger(fhicl::ParameterSet const & pset);
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

  /**
     @brief given a cluster::Cluster fill a recob::Cluster
     @input CMCluster -> where cluster variables are currently stored
     @input n         -> index of current cluster
     @return          -> recob::Cluster
   */
  const recob::Cluster FillClusterProperties(const ::cluster::Cluster& CMCluster,
					     const size_t& n);

  // Declare member data here.

  float _wire2cm, _time2cm;

  // ClusterMaker class
  ::cluster::ClusterMaker* _CMaker;
  
  // cluster helper
  ::cmtool::CMergeHelper* _merge_helper;

  std::string fClusterProducer, fVertexProducer;

  /**
     Given the FindManyP cluster -> hit association
     grab the ProductID and EDProductGetter 
   */
  void GetHitPointer(const art::FindManyP<recob::Hit>& associations);
  art::Ptr<recob::Hit> _hitptr;

};


ClusterMerger::ClusterMerger(fhicl::ParameterSet const & pset)
// :
// Initialize member data here.
{

  std::cout << "DD calling constructor" << std::endl;

  _merge_helper = new ::cmtool::CMergeHelper();

  std::cout << "DD reset Merge manager" << std::endl;
  _merge_helper->GetManager().Reset();
  std::cout << "DD done with reset" << std::endl;

  _CMaker = new ::cluster::ClusterMaker();
  
  // get detector specific properties
  auto const* geom = ::lar::providerFrom<geo::Geometry>();
  auto const* detp = lar::providerFrom<detinfo::DetectorPropertiesService>();
  _wire2cm = geom->WirePitch(0,1,0);
  _time2cm = detp->SamplingRate() / 1000.0 * detp->DriftVelocity( detp->Efield(), detp->Temperature() );

  fClusterProducer = pset.get<std::string>("ClusterProducer");
  fVertexProducer  = pset.get<std::string>("VertexProducer" );

  std::cout << "DD setting up algos" << std::endl;
  
  // grab algorithms for merging
  const fhicl::ParameterSet& mergeTools = pset.get<fhicl::ParameterSet>("MergeTools");
  std::cout << "DD got parameter set. start loop..." << std::endl;
  for (const std::string& mergeTool : mergeTools.get_pset_names()) {
    std::cout << "DD \t in loop..." << std::endl;
    const fhicl::ParameterSet& merge_pset = mergeTools.get<fhicl::ParameterSet>(mergeTool);
    std::cout << "DD \t add merge algor..." << std::endl;
    _merge_helper->GetManager().AddMergeAlgo(art::make_tool<cmtool::CBoolAlgoBase>(merge_pset));
    std::cout << "DD \t done adding algo" << std::endl;
  }// for all algorithms to be added

  _merge_helper->GetManager().ReportAlgoChain();

  produces<std::vector<recob::Cluster> >();
  produces<art::Assns <recob::Cluster, recob::Hit> >();

  std::cout << "DD done with constructor" << std::endl;

}

void ClusterMerger::produce(art::Event & e)
{
  // Implementation of required member function here.

  std::cout << "DD entering producer" << std::endl;

  std::unique_ptr< std::vector<recob::Cluster> > Cluster_v(new std::vector<recob::Cluster>);
  std::unique_ptr< art::Assns <recob::Cluster, recob::Hit> > Cluster_Hit_assn_v(new art::Assns<recob::Cluster,recob::Hit>);

  // cluster pointer maker for later to create associations
  art::PtrMaker<recob::Cluster> ClusPtrMaker(e, *this);

  std::cout << "DD created empty container storer for product objects" << std::endl;

  // load data products needed

  // load input clusters
  auto const& clus_h = e.getValidHandle<std::vector<recob::Cluster>>(fClusterProducer);

   std::cout << "DD loaded " << clus_h->size() << " clusters..." << std::endl;
  
  // load associated hits
  art::FindManyP<recob::Hit> clus_hit_assn_v(clus_h, e, fClusterProducer);

  std::cout << "DD and associations..." << std::endl;

  // get generic hit art::Ptr
  // from it get id() and productGetter()
  // see http://nusoft.fnal.gov/larsoft/doxsvn/html/classart_1_1Ptr.html
  // and then can call Ptr (ProductID const productID, key_type itemKey, EDProductGetter const *prodGetter)
  // to make art::Ptr for the hit I want to associate.
  GetHitPointer(clus_hit_assn_v);

  std::cout << "DD got hit pointer to create output associations..." << std::endl;

  // load vertices
  auto const& vtx_h = e.getValidHandle<std::vector<recob::Vertex>>(fVertexProducer);

  std::cout << "DD vertex loaded" << std::endl;

  // create cluster::Clusters
  std::vector<::cluster::Cluster> event_clusters;
  _CMaker->MakeClusters(clus_h, clus_hit_assn_v, vtx_h, event_clusters);

  std::cout << "DD Make clusters called" << std::endl;

  _merge_helper->Process(event_clusters);

  // get output clusters
  auto out_clusters = _merge_helper->GetClusters();

  std::cout << "there are " << out_clusters.size() << " outputs" << std::endl;

  // loop through output clusters
  for (auto const& out_clus : out_clusters) {

    Cluster_v->emplace_back( FillClusterProperties(out_clus,Cluster_v->size()) );
    art::Ptr<recob::Cluster> const ClusPtr = ClusPtrMaker(Cluster_v->size()-1);

    // create association to hits
    for (auto const& hit : out_clus.GetHits()) {
      auto key = hit._idx;
      art::Ptr<recob::Hit> HitPtr(_hitptr.id(),key,_hitptr.productGetter());
      Cluster_Hit_assn_v->addSingle(ClusPtr,HitPtr);
    }// for all hits associated to the cluster
    
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

const recob::Cluster ClusterMerger::FillClusterProperties(const ::cluster::Cluster& CMCluster,
							  const size_t& n) {
  
  auto const* geom = ::lar::providerFrom<geo::Geometry>();
  
  float startW = CMCluster._start_pt._w;
  float startT = CMCluster._start_pt._t;
  
  float endW   = CMCluster._end_pt._w;
  float endT   = CMCluster._end_pt._t;
  
  auto planeid = geo::PlaneID(0,0,CMCluster._plane);

  recob::Cluster clus(startW, 0., startT, 0., 0., CMCluster._angle, 0., 
		      endW,   0., endT,   0., 0., 0., 0., 
		      CMCluster._sum_charge, 0., CMCluster._sum_charge, 0., 
		      CMCluster.GetHits().size(), 0., 0., n,
		      geom->View(planeid),
		      planeid);

  return clus;
}

/**
   This function returns a pointer to a hit in the event and is used to be
   able to create associations via addSingle once the modules have finished running
 */
void ClusterMerger::GetHitPointer(const art::FindManyP<recob::Hit>& associations)
{

  // find a valid hit pointer!
  if (associations.size()) {
    if (associations.at(0).size()) {
      
      _hitptr = associations.at(0).at(0);

    }// if hit exists
  }// if cluster exists


  return;
}


DEFINE_ART_MODULE(ClusterMerger)
