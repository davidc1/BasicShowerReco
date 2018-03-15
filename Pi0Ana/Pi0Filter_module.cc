////////////////////////////////////////////////////////////////////////
// Class:       Pi0Filter
// Plugin Type: filter (art v2_09_06)
// File:        Pi0Filter_module.cc
//
// Generated at Sat Feb 24 08:47:57 2018 by David Caratelli using cetskelgen
// from cetlib version v3_01_03.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art/Framework/Services/Optional/TFileService.h"

#include <memory>

#include "TTree.h"
#include "TVector3.h"

#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"

class Pi0Filter;

class Pi0Filter : public art::EDFilter {
public:
  explicit Pi0Filter(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  Pi0Filter(Pi0Filter const &) = delete;
  Pi0Filter(Pi0Filter &&) = delete;
  Pi0Filter & operator = (Pi0Filter const &) = delete;
  Pi0Filter & operator = (Pi0Filter &&) = delete;

  // Required functions.
  bool filter(art::Event & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.

  TTree* _tree;
  float  _e1;
  float  _e2;
  float  _dedx1;
  float  _dedx2;
  float  _angle;
  float  _mass;
  int    _nshr;
  int _run,_sub,_evt;

  TTree* _trkangle_tree;
  float  _trkangle;

};


Pi0Filter::Pi0Filter(fhicl::ParameterSet const & p)
// :
// Initialize member data here.
{
  
}

bool Pi0Filter::filter(art::Event & e)
{

  _run = e.run();
  _sub = e.subRun();
  _evt = e.event();

  // load input showers
  auto const& shr_h = e.getValidHandle<std::vector<recob::Shower>>("showerreco3d");
  // load input tracks
  auto const& trk_h = e.getValidHandle<std::vector<recob::Track>>("pandoraCosmic");
  // load input vertices
  auto const& vtx_h = e.getValidHandle<std::vector<recob::Vertex>>("ccvertex");

  // check whether vertex in the middle of a cosmic track
  if (vtx_h->size() != 1) return false;

  auto const& vtx = vtx_h->at(0);
  Double_t xyz[3] = {};
  vtx.XYZ(xyz);
  TVector3 vtxpt(xyz[0],xyz[1],xyz[2]);
  _trkangle = 0.;

  // loop trhough tracks. is the angle ~180?
  for (size_t t=0; t < trk_h->size(); t++) {

    auto const& trk = trk_h->at(t);

    if (trk.Length() < 50) continue;

    auto const& start = trk.Vertex();
    auto const& end   = trk.End();

    if ( ( (start-vtxpt).Mag() < 10) && ((end-vtxpt).Mag() < 10) ) continue;

    float a = (start-vtxpt).Angle(end-vtxpt);
    if (a > _trkangle)
      _trkangle = a;
  }// for all tracks

  _trkangle_tree->Fill();



  _nshr = shr_h->size();

  // if less then two showers -> exit
  if (shr_h->size() < 2){
    _tree->Fill();
    return false;
  }

  // otherwise, check for pi0
  // grab two showers with highest energy
  recob::Shower shr1, shr2;
  _e1 = 0.;
  _e2 = 0.;
  
  // find largest energy shower
  for (size_t s=0; s < shr_h->size(); s++) {
    auto const& shr = shr_h->at(s);
    if (shr.Energy()[2] > _e1) {
      _e1  = shr.Energy()[2];
      shr1 = shr;
    }
  }// for all showers
  // find second largest energy shower
  for (size_t s=0; s < shr_h->size(); s++) {
    auto const& shr = shr_h->at(s);
    if ( (shr.Energy()[2] > _e2) && ((float)shr.Energy()[2] != _e1) ) {
      _e2  = shr.Energy()[2];
      shr2 = shr;
    }
  }// for all showers

  _dedx1 = shr1.dEdx()[2];
  _dedx2 = shr2.dEdx()[2];

  // get opening angle
  auto dir1 = shr1.Direction();
  auto dir2 = shr2.Direction();
  _angle = dir1.Angle(dir2);

  std::cout << "\t opening angle : " << _angle << std::endl;

  _mass = sqrt( 2 * _e1 * _e2 * ( 1 - cos(_angle) ) );

  std::cout << "\t mass          : " << _mass << std::endl << std::endl;

  _tree->Fill();


  if (_trkangle > 2.9) return false;

  if (_e2 < 35.) return false;
  
  std::cout << "Showers of energy " << _e1 << " and " << _e2 << std::endl;
  
  return true;
}

void Pi0Filter::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  _tree = tfs->make<TTree>("_tree","Pi0 Tree TTree");
  _tree->Branch("_run",&_run,"run/I");
  _tree->Branch("_sub",&_sub,"sub/I");
  _tree->Branch("_evt",&_evt,"evt/I");
  _tree->Branch("_nshr",&_nshr,"nshr/I");
  _tree->Branch("_e1",&_e1,"e1/F");  
  _tree->Branch("_e2",&_e2,"e2/F");
  _tree->Branch("_dedx1",&_dedx1,"dedx1/F");
  _tree->Branch("_dedx2",&_dedx2,"dedx2/F");
  _tree->Branch("_angle",&_angle,"angle/F");
  _tree->Branch("_mass" ,&_mass ,"mass/F" );

  _trkangle_tree = tfs->make<TTree>("_trkangle_tree","Track Angle TTree");
  _trkangle_tree->Branch("_trkangle",&_trkangle,"trkangle/F");

}

void Pi0Filter::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(Pi0Filter)
