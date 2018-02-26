////////////////////////////////////////////////////////////////////////
// Class:       Showerreco
// Plugin Type: producer (art v2_09_06)
// File:        ShrReco3D_module.cc
//
// Generated at Fri Feb  9 16:38:52 2018 by David Caratelli using cetskelgen
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

// shower-reco classes and utilities
#include "ProtoShower/ProtoShowerAlgBase.h"
#include "Base/ShowerRecoManager.h"
// include specific protoshower and recomanager instances
#include "ProtoShower/ProtoShowerCMTool.h"

// larsoft data-products
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardata/Utilities/AssociationUtil.h"

#include "art/Utilities/make_tool.h"

#include "art/Persistency/Common/PtrMaker.h"

#include <memory>

class ShrReco3D;


class ShrReco3D : public art::EDProducer {
public:
  explicit ShrReco3D(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.
  
  // Plugins should not be copied or assigned.
  ShrReco3D(ShrReco3D const &) = delete;
  ShrReco3D(ShrReco3D &&) = delete;
  ShrReco3D & operator = (ShrReco3D const &) = delete;
  ShrReco3D & operator = (ShrReco3D &&) = delete;
  
  // Required functions.
  void produce(art::Event & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;
  
private:
  
  /// Input producer name
  std::string fPFPproducer;
  std::string fClusproducer;
  std::string fVtxproducer;
  
  /// Shower reco core class instance
  ::showerreco::ShowerRecoManager* _manager;
  
  // ProtoShowerAlgBase to make protoshowers
  ::protoshower::ProtoShowerAlgBase* _psalg;
  
  /**
     @brief Save output showers produced by reconstruction algorithms
     @input s : index in output shower vector, to associate back to PFP index
     @input shower : reconstructed shower to be translated to recob:: object
     @input Shower_v : vector of recob::Showers. This will be placed in the art::Event
     @input pfp_h : handle to input pfparticles
     @input pfp_clus_assn_v : input pfp -> clus assn vector
     @input pfp_hit_assn_v  : input pfp -> hit assn vector
     @input Shower_PFP_assn_v     : output shower -> pfp assn vector
     @input Shower_Cluster_assn_v : output shower -> cluster assn vector
     @input Shower_Hit_assn_v     : output shower -> hit assn vector
  */
  void SaveShower(const size_t idx,
		  const showerreco::Shower_t& shower,
		  std::unique_ptr< std::vector<recob::Shower> >& Shower_v,
		  const art::PtrMaker<recob::Shower> ShowerPtrMaker,
		  const art::ValidHandle<std::vector<recob::PFParticle> > pfp_h,
		  const art::FindManyP<recob::Cluster> pfp_clus_assn_v,
		  const std::vector<std::vector<art::Ptr<recob::Hit> > > pfp_hit_assn_v,
		  std::unique_ptr< art::Assns <recob::Shower, recob::PFParticle> >& Shower_PFP_assn_v,
		  std::unique_ptr< art::Assns <recob::Shower, recob::Cluster>    >& Shower_Cluster_assn_v,
		  std::unique_ptr< art::Assns <recob::Shower, recob::Hit>        >& Shower_Hit_assn_v );
  
  // Declare member data here.
  
};


ShrReco3D::ShrReco3D(fhicl::ParameterSet const & p)
// :
// Initialize member data here.
{

  fPFPproducer  = p.get<std::string>("PFPproducer" );
  fClusproducer = p.get<std::string>("Clusproducer");
  fVtxproducer  = p.get<std::string>("Vtxproducer" );

  // grab algorithms for merging
  _manager = new showerreco::ShowerRecoManager();
  const fhicl::ParameterSet& showerrecoTools = p.get<fhicl::ParameterSet>("ShowerRecoTools");
  std::cout << "DD got parameter set. start loop..." << std::endl;
  for (const std::string& showerrecoTool : showerrecoTools.get_pset_names()) {
    std::cout << "DD \t in loop..." << std::endl;
    const fhicl::ParameterSet& showerreco_pset = showerrecoTools.get<fhicl::ParameterSet>(showerrecoTool);
    std::cout << "DD \t add merge algo..." << std::endl;
    _manager->AddAlgo(art::make_tool<showerreco::ShowerRecoModuleBase>(showerreco_pset));
    std::cout << "DD \t done adding algo" << std::endl;
  }// for all algorithms to be added

  _manager->SetDebug(false);

  //_manager = new ::showerreco::Pi0RecoAlgorithm();
  _psalg   = new ::protoshower::ProtoShowerCMTool();

  produces<std::vector<recob::Shower> >();
  produces<art::Assns <recob::Shower, recob::PFParticle> >();
  produces<art::Assns <recob::Shower, recob::Cluster>    >();
  produces<art::Assns <recob::Shower, recob::Hit>        >();

  _manager->Initialize();

}

void ShrReco3D::produce(art::Event & e)
{

  // produce recob::Showers
  std::unique_ptr< std::vector<recob::Shower> > Shower_v(new std::vector<recob::Shower> );
  std::unique_ptr< art::Assns <recob::Shower, recob::PFParticle> > Shower_PFP_assn_v    ( new art::Assns<recob::Shower, recob::PFParticle>);
  std::unique_ptr< art::Assns <recob::Shower, recob::Cluster>    > Shower_Cluster_assn_v( new art::Assns<recob::Shower, recob::Cluster>   );
  std::unique_ptr< art::Assns <recob::Shower, recob::Hit>        > Shower_Hit_assn_v    ( new art::Assns<recob::Shower, recob::Hit>       );

  // shower pointer maker for later to create associations
  art::PtrMaker<recob::Shower> ShowerPtrMaker(e, *this);

  // pass event to ProtoShowerAlgBase to create ProtoShower objects
  // which will then be fed to shower reco algorithm chain
  std::vector<protoshower::ProtoShower> event_protoshower_v;
  _psalg->GenerateProtoShowers(e, fPFPproducer, fClusproducer, fVtxproducer, event_protoshower_v);
  
  // set protoshowers for algorithms
  _manager->SetProtoShowers(event_protoshower_v);

  // output showers to be saved to event
  std::vector< ::showerreco::Shower_t> output_shower_v;
  _manager->Reconstruct(output_shower_v);


  // load PFP, clus, hit so that associations to showers can be stored
  // grab PFParticles in event
  auto const& pfp_h = e.getValidHandle<std::vector<recob::PFParticle> >(fPFPproducer);
  // grab clusters associated with PFParticles
  art::FindManyP<recob::Cluster> pfp_clus_assn_v(pfp_h, e, fPFPproducer);
  // ADDITION FROM PETRILLO
  e.getValidHandle<std::vector<recob::Cluster>>(fClusproducer);
  // grab the hits associated to the PFParticles
  auto pfp_hit_assn_v = lar::FindManyInChainP<recob::Hit, recob::Cluster>::find(pfp_h, e, fPFPproducer);

  // save output showers
  for (size_t s=0; s < output_shower_v.size(); s++) {
    SaveShower(s, output_shower_v.at(s), Shower_v, ShowerPtrMaker,
	       pfp_h, pfp_clus_assn_v, pfp_hit_assn_v,
	       Shower_PFP_assn_v, Shower_Cluster_assn_v, Shower_Hit_assn_v);
  }// for all output reconstructed showers

  e.put(std::move(Shower_v));
  e.put(std::move(Shower_PFP_assn_v));
  e.put(std::move(Shower_Cluster_assn_v));
  e.put(std::move(Shower_Hit_assn_v));

}

void ShrReco3D::beginJob()
{
  // Implementation of optional member function here.
}

void ShrReco3D::endJob()
{
  _manager->Finalize();
}


void ShrReco3D::SaveShower(const size_t idx,
			   const showerreco::Shower_t& shower,
			   std::unique_ptr< std::vector<recob::Shower> >& Shower_v,
			   const art::PtrMaker<recob::Shower> ShowerPtrMaker,
			   const art::ValidHandle<std::vector<recob::PFParticle> > pfp_h,
			   const art::FindManyP<recob::Cluster> pfp_clus_assn_v,
			   const std::vector<std::vector<art::Ptr<recob::Hit> > > pfp_hit_assn_v,
			   std::unique_ptr< art::Assns <recob::Shower, recob::PFParticle> >& Shower_PFP_assn_v,
			   std::unique_ptr< art::Assns <recob::Shower, recob::Cluster> >& Shower_Cluster_assn_v,
			   std::unique_ptr< art::Assns <recob::Shower, recob::Hit> >& Shower_Hit_assn_v)
  
{
  
  if (shower.fPassedReconstruction == false) {
    return;
  }
  
  // filter out showers with garbage values
  if (shower.fXYZStart.Mag2()  == 0) {
    return;
  }
  if (shower.fDCosStart.Mag2() == 0) {
    return;
  }

  recob::Shower s;
  s.set_id ( Shower_v->size() );
  s.set_total_energy          ( shower.fTotalEnergy_v         );
  s.set_total_energy_err      ( shower.fSigmaTotalEnergy_v    );
  s.set_total_best_plane      ( shower.fBestPlane.Plane       );
  s.set_direction             ( shower.fDCosStart             );
  s.set_direction_err         ( shower.fSigmaDCosStart        );
  s.set_start_point           ( shower.fXYZStart              );
  s.set_start_point_err       ( shower.fSigmaXYZStart         );
  s.set_dedx                  ( shower.fdEdx_v                );
  s.set_dedx_err              ( shower.fSigmadEdx_v           );
  s.set_length                ( shower.fLength                );
  s.set_open_angle            ( shower.fOpeningAngle          );
  
  Shower_v->emplace_back(s);

  art::Ptr<recob::Shower> const ShrPtr = ShowerPtrMaker(Shower_v->size()-1);

  // now take care of associations
  
  // step 1 : pfp
  const art::Ptr<recob::PFParticle> PFPPtr(pfp_h, idx);
  Shower_PFP_assn_v->addSingle( ShrPtr, PFPPtr );

  // step 2 : clusters
  std::vector<art::Ptr<recob::Cluster> > clus_v = pfp_clus_assn_v.at(idx);
  for (size_t c=0; c < clus_v.size(); c++)
    Shower_Cluster_assn_v->addSingle( ShrPtr, clus_v.at(c) );

  // step 3 : hits
  std::vector<art::Ptr<recob::Hit> > hit_v = pfp_hit_assn_v.at(idx);
  for (size_t h=0; h < hit_v.size(); h++)
    Shower_Hit_assn_v->addSingle( ShrPtr, hit_v.at(h) );
  
  return;
}

DEFINE_ART_MODULE(ShrReco3D)
