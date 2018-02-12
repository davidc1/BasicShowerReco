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
#include "Factory/Pi0RecoAlgorithm.h"

// larsoft data-products
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardata/Utilities/AssociationUtil.h"


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
    std::string fInputProducer;

    /// Shower reco core class instance
    ::showerreco::ShowerRecoManager* _manager;

    // ProtoShowerAlgBase to make protoshowers
    ::protoshower::ProtoShowerAlgBase* _psalg;

  /**
     Save output showers produced by reconstruction algorithms
   */
  void SaveShower(const showerreco::Shower_t& shower,
		  std::unique_ptr< std::vector<recob::Shower> >& Shower_v);

  // Declare member data here.

};


ShrReco3D::ShrReco3D(fhicl::ParameterSet const & p)
// :
// Initialize member data here.
{

  fInputProducer  = p.get<std::string>("InputProducer" );

  produces<std::vector<recob::Shower> >();
  produces<art::Assns <recob::Shower, recob::PFParticle> >();
  produces<art::Assns <recob::Shower, recob::Cluster>    >();
  produces<art::Assns <recob::Shower, recob::Hit>        >();

  _manager = new ::showerreco::Pi0RecoAlgorithm();
  _psalg   = new ::protoshower::ProtoShowerCMTool();

}

void ShrReco3D::produce(art::Event & e)
{

  // produce recob::Showers
  std::unique_ptr< std::vector<recob::Shower> > Shower_v(new std::vector<recob::Shower> );
  std::unique_ptr< art::Assns <recob::Shower, recob::PFParticle> > Shower_PFP_assn_v    ( new art::Assns<recob::Shower, recob::PFParticle>);
  std::unique_ptr< art::Assns <recob::Shower, recob::Cluster>    > Shower_Cluster_assn_v( new art::Assns<recob::Shower, recob::Cluster>   );
  std::unique_ptr< art::Assns <recob::Shower, recob::Hit>        > Shower_Hit_assn_v    ( new art::Assns<recob::Shower, recob::Hit>       );

  // pass event to ProtoShowerAlgBase to create ProtoShower objects
  // which will then be fed to shower reco algorithm chain
  std::vector<protoshower::ProtoShower> event_protoshower_v;
  _psalg->GenerateProtoShowers(e, fInputProducer, event_protoshower_v);
  
  // set protoshowers for algorithms
  _manager->SetProtoShowers(event_protoshower_v);

  // output showers to be saved to event
  std::vector< ::showerreco::Shower_t> output_shower_v;
  _manager->Reconstruct(output_shower_v);

  // save output showers
  for (size_t s=0; s < output_shower_v.size(); s++) {
    SaveShower(output_shower_v.at(s), Shower_v);
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
  // Implementation of optional member function here.
}


void ShrReco3D::SaveShower(const showerreco::Shower_t& shower,
			      std::unique_ptr< std::vector<recob::Shower> >& Shower_v) 
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

    return;
}

DEFINE_ART_MODULE(ShrReco3D)