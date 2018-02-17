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

#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/PFParticle.h"

#include "art/Persistency/Common/PtrMaker.h"

#include "art/Utilities/make_tool.h"

#include "uboone/BasicShowerReco/ClusterMerging/CMToolBase/ClusterMaker.h"
#include "uboone/BasicShowerReco/ClusterMerging/CMToolBase/CMatchManager.h"

#include <memory>

class ClusterMatcher;


class ClusterMatcher : public art::EDProducer {
public:
  explicit ClusterMatcher(fhicl::ParameterSet const & pset);
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

  /**
     @brief given a cluster::Cluster fill a recob::Cluster
     @input CMCluster -> where cluster variables are currently stored
     @input n         -> index of current cluster
     @return          -> recob::Cluster
   */
  const recob::Cluster FillClusterProperties(const ::cluster::Cluster& CMCluster,
					     const size_t& n);

  float _wire2cm, _time2cm;

  // cluster maker
  ::cluster::ClusterMaker* _CMaker;

  std::string fClusterProducer;

  // cluster match manager
  ::cmtool::CMatchManager* _mgr;

};


ClusterMatcher::ClusterMatcher(fhicl::ParameterSet const & pset)
// :
// Initialize member data here.
{

  // get detector specific properties
  auto const* geom = ::lar::providerFrom<geo::Geometry>();
  auto const* detp = lar::providerFrom<detinfo::DetectorPropertiesService>();
  _wire2cm = geom->WirePitch(0,1,0);
  _time2cm = detp->SamplingRate() / 1000.0 * detp->DriftVelocity( detp->Efield(), detp->Temperature() );

  // Call appropriate produces<>() functions here.
  std::cout << "DD calling constructor" << std::endl;

  _mgr = new ::cmtool::CMatchManager(geom->Nplanes());

  std::cout << "DD reset Merge manager" << std::endl;
  _mgr->Reset();
  std::cout << "DD done with reset" << std::endl;

  _CMaker = new ::cluster::ClusterMaker();
  
  fClusterProducer = pset.get<std::string>("ClusterProducer");

  std::cout << "DD setting up algos" << std::endl;
  
  // grab algorithms for merging
  const fhicl::ParameterSet& matchTool = pset.get<fhicl::ParameterSet>("MatchTool");
  _mgr->AddMatchAlgo(art::make_tool<cmtool::CFloatAlgoBase>(matchTool));
  std::cout << "DD \t done adding algo" << std::endl;

  _mgr->ReportAlgoChain();

  produces<std::vector<recob::PFParticle> >();
  produces<art::Assns <recob::PFParticle, recob::Cluster> >();

  std::cout << "DD done with constructor" << std::endl;

}

void ClusterMatcher::produce(art::Event & e)
{
  // Implementation of required member function here.

  std::unique_ptr< std::vector<recob::PFParticle> > PFP_v(new std::vector<recob::PFParticle>);
  std::unique_ptr< art::Assns <recob::PFParticle, recob::Cluster> > PFP_Hit_assn_v(new art::Assns<recob::PFParticle,recob::Cluster>);

  // cluster pointer maker for later to create associations
  art::PtrMaker<recob::PFParticle> PFPPtrMaker(e, *this);

  // load input clusters
  auto const& clus_h = e.getValidHandle<std::vector<recob::Cluster>>(fClusterProducer);

  // load associated hits
  art::FindManyP<recob::Hit> clus_hit_assn_v(clus_h, e, fClusterProducer);

  _mgr->Reset();
  //_mgr->SetClusters();//local_clusters);
  _mgr->Process();


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
