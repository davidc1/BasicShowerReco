////////////////////////////////////////////////////////////////////////
// Class:       CCTrackAna
// Plugin Type: analyzer (art v2_10_03)
// File:        CCTrackAna_module.cc
//
// Generated at Sat Jul 28 14:54:52 2018 by David Caratelli using cetskelgen
// from cetlib version v3_02_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art/Framework/Services/Optional/TFileService.h"

#include "TTree.h"
#include "TVector3.h"

#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardata/Utilities/AssociationUtil.h"

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardata/Utilities/GeometryUtilities.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

class CCTrackAna;


class CCTrackAna : public art::EDAnalyzer {
public:
  explicit CCTrackAna(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CCTrackAna(CCTrackAna const &) = delete;
  CCTrackAna(CCTrackAna &&) = delete;
  CCTrackAna & operator = (CCTrackAna const &) = delete;
  CCTrackAna & operator = (CCTrackAna &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

  void SetTTree();

private:

  // Declare member data here.

  TTree* _trk_tree;

  // variables common to both ttrees
  int _run, _sub, _evt;

  int _n_reco_showers;


  double _rc_vtx_x, _rc_vtx_y, _rc_vtx_z;

  std::vector<float> _dqdx_v, _dedx_v, _rr_v;
  double _dvtx, _length;
  double _xs, _ys, _zs, _xe, _ye, _ze;
  int _pl;

  std::string fTrkProducer;
  std::string fCaloProducer;
  std::string fAssnProducer;


};


CCTrackAna::CCTrackAna(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
 // More initializers here.
{
  fTrkProducer  = p.get<std::string>("TrkProducer" );
  fCaloProducer = p.get<std::string>("CaloProducer");
  fAssnProducer = p.get<std::string>("AssnProducer");
  SetTTree();
}

void CCTrackAna::analyze(art::Event const & e)
{
  // Implementation of required member function here.

  art::Handle< art::Assns<recob::Vertex,recob::Track,void> > numuCCassn_h;
  e.getByLabel(fAssnProducer,numuCCassn_h);
  if (numuCCassn_h->size() != 1){
    std::cout << "Number of vertices != 1 -> ERROR ERROR ERROR" << std::endl;
    return;
  }
  
  Double_t xyz[3] = {};  
  numuCCassn_h->at(0).first->XYZ(xyz);
  _rc_vtx_x = xyz[0];
  _rc_vtx_y = xyz[1];
  _rc_vtx_z = xyz[2];

  // load input tracks
  auto const& trk_h = e.getValidHandle<std::vector<recob::Track>>(fTrkProducer);
  // load calorimetry associated to tracks
  art::FindManyP<anab::Calorimetry> trk_calo_assn_v(trk_h, e, fCaloProducer);

  // track portion of this analysis
  TVector3 nuvtx(_rc_vtx_x,_rc_vtx_y,_rc_vtx_z);
  for (size_t t=0; t < trk_h->size(); t++) {
    
    auto const& trk = trk_h->at(t);
    auto const& beg = trk.Vertex();
    auto const& end = trk.End();
    
    _xs = beg.X();
    _ys = beg.Y();
    _zs = beg.Z();
    _xe = end.X();
    _ye = end.Y();
    _ze = end.Z();

    double dvtx = 1000.;
    // which one is closest?
    if ( (nuvtx-beg).Mag() < (nuvtx-end).Mag() ) 
      dvtx = (nuvtx-beg).Mag();
    else 
      dvtx = (nuvtx-end).Mag();

    if (dvtx < 5.0) {

      _length = trk.Length();
      _dvtx   = dvtx;
      
      // fill calorimetry info for this track
      // grab the associated calorimetry object
      //const std::vector<const anab::Calorimetry*>& Calo_v = trk_calo_assn_v.at(t);
      auto Calo_v = trk_calo_assn_v.at(t);
      
      for (size_t pl=0; pl < Calo_v.size(); pl++){
	
	auto const& calo = Calo_v.at(pl);
	
	_pl = calo->PlaneID().Plane;
	
	// grab point-by-point information
	auto const& dqdx_v = calo->dQdx();
	auto const& dedx_v = calo->dEdx();
	auto const& rr_v   = calo->ResidualRange();
	_dqdx_v.clear();
	_dedx_v.clear();
	_rr_v.clear();
	for (auto const& dqdx : dqdx_v)
	  _dqdx_v.push_back((float)dqdx);
	for (auto const& dedx : dedx_v)
	  _dedx_v.push_back((float)dedx);
	for (auto const& rr : rr_v)
	  _rr_v.push_back((float)rr);

	_trk_tree->Fill();
	
      }// for all planes
    }// if within 5 cm of vertex
  }// for all tracks

  return;
}

void CCTrackAna::beginJob()
{
  // Implementation of optional member function here.
}

void CCTrackAna::endJob()
{
  // Implementation of optional member function here.
}

void CCTrackAna::SetTTree() {

  art::ServiceHandle<art::TFileService> tfs;

  // track tree
  _trk_tree = tfs->make<TTree>("_trk_tree","track tree");
  _trk_tree->Branch("_rc_vtx_x",&_rc_vtx_x,"rc_vtx_x/D");
  _trk_tree->Branch("_rc_vtx_y",&_rc_vtx_y,"rc_vtx_y/D");
  _trk_tree->Branch("_rc_vtx_z",&_rc_vtx_z,"rc_vtx_z/D");
  _trk_tree->Branch("_run",&_run,"run/I");
  _trk_tree->Branch("_sub",&_sub,"sub/I");
  _trk_tree->Branch("_evt",&_evt,"evt/I");
  _trk_tree->Branch("_dqdx_v","std::vector<float>",&_dqdx_v);
  _trk_tree->Branch("_dedx_v","std::vector<float>",&_dedx_v);
  _trk_tree->Branch("_rr_v"  ,"std::vector<float>",&_rr_v  );
  _trk_tree->Branch("_pl",&_pl,"pl/I");
  _trk_tree->Branch("_dvtx",&_dvtx,"dvtx/D");
  _trk_tree->Branch("_length",&_length,"length/D");
  _trk_tree->Branch("_xs",&_xs,"xs/D");
  _trk_tree->Branch("_ys",&_ys,"ys/D");
  _trk_tree->Branch("_zs",&_zs,"zs/D");
  _trk_tree->Branch("_xe",&_xe,"xe/D");
  _trk_tree->Branch("_ye",&_ye,"ye/D");
  _trk_tree->Branch("_ze",&_ze,"ze/D");
  
  return;
}


DEFINE_ART_MODULE(CCTrackAna)
