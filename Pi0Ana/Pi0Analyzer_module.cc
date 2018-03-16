////////////////////////////////////////////////////////////////////////
// Class:       Pi0Analyzer
// Plugin Type: analyzer (art v2_09_06)
// File:        Pi0Analyzer_module.cc
//
// Generated at Wed Mar 14 09:23:33 2018 by David Caratelli using cetskelgen
// from cetlib version v3_01_03.
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

#include <memory>

#include "TTree.h"
#include "TVector3.h"

#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardata/Utilities/AssociationUtil.h"

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardata/Utilities/GeometryUtilities.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

class Pi0Analyzer;


class Pi0Analyzer : public art::EDAnalyzer {
public:
  explicit Pi0Analyzer(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  Pi0Analyzer(Pi0Analyzer const &) = delete;
  Pi0Analyzer(Pi0Analyzer &&) = delete;
  Pi0Analyzer & operator = (Pi0Analyzer const &) = delete;
  Pi0Analyzer & operator = (Pi0Analyzer &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  void Reset();

  //std::vector<int> Match(const sim::MCShower& mcs1, const sim::MCShower& mcs2,
  //			 const std::vector<recob::Shower>&   shr_v);
  std::vector<int> Match(const std::vector<sim::MCShower>& pi0_shower_v,
			 const std::vector<recob::Shower>&   shr_v);
  
  void SetTTree();

  TTree* _tree;

  int _run, _sub, _evt;

  double _w2cm, _t2cm;

  int _n_reco_showers;

  int _event;

  double _nu_e, _pi0_e;

  double _mc_vtx_x, _mc_vtx_y, _mc_vtx_z;
  double _rc_vtx_x, _rc_vtx_y, _rc_vtx_z;

  double _mc_shr1_x,  _mc_shr1_y,  _mc_shr1_z;
  double _mc_shr1_px, _mc_shr1_py, _mc_shr1_pz;
  double _mc_shr1_e;
  double _mc_shr1_cont;
  double _mc_shr2_x,  _mc_shr2_y,  _mc_shr2_z;
  double _mc_shr2_px, _mc_shr2_py, _mc_shr2_pz;
  double _mc_shr2_e;
  double _mc_shr2_cont;
  double _mcradlen1, _mcradlen2;

  double _rc_shr_x,  _rc_shr_y,  _rc_shr_z;
  double _rc_shr_px, _rc_shr_py, _rc_shr_pz;
  double _rc_shr_e;
  double _rcradlen;

  double _mc_oangle;

  double _mc_mass;

  // MC -> RC shower comparisons                                                                                                                                                                                         
  double _dot;
  double _strt;
  double _emc;

  // cluster metrics                                                                                                                                                                                                     
  double _ip;
  double _lin;
  double _ssv;
  double _slope;
  double _slopedirangle;
  double _hitshowerangle;
  
  double _dwallmin;

  // Declare member data here.

};


Pi0Analyzer::Pi0Analyzer(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
 // More initializers here.
{
  SetTTree();
}

void Pi0Analyzer::analyze(art::Event const & e)
{

  _run = e.run();
  _sub = e.subRun();
  _evt = e.event();

  // load input showers
  auto const& shr_h = e.getValidHandle<std::vector<recob::Shower>>("showerreco3d");
  // load input tracks
  //auto const& trk_h = e.getValidHandle<std::vector<recob::Track>>("pandoraCosmic");
  // load input vertices
  auto const& vtx_h = e.getValidHandle<std::vector<recob::Vertex>>("ccvertex");
  // load mcshowers
  auto const& mcs_h = e.getValidHandle<std::vector<sim::MCShower>>("mcreco");

  Double_t xyz[3] = {};

  auto const& mct_h = e.getValidHandle<std::vector<simb::MCTruth> >("generator");

  //std::cout << "got data-products" << std::endl;

  auto mct = mct_h->at(0);
  size_t npart = mct.NParticles();

  bool foundPi0 = false;

  for (size_t i=0; i < npart; i++){
    auto const& part = mct.GetParticle(i);
    if ( (part.PdgCode() == 111) and (part.StatusCode() == 1) ){
      xyz[0] = part.Trajectory().X(0);
      xyz[1] = part.Trajectory().Y(0);
      xyz[2] = part.Trajectory().Z(0);
      foundPi0 = true;
      break;
    }
  }

  //std::cout << " ***** 1 ******" << std::endl;

  if (foundPi0 == false)
    return;

  size_t idx_1 = 0;
  size_t idx_2 = 0;
  size_t n_found = 0;
  for (size_t i=0; i < mcs_h->size(); i++){
    auto const& mcs = mcs_h->at(i);
    // distance from vertex                                                                
    double x = mcs.Start().X();
    double y = mcs.Start().Y();
    double z = mcs.Start().Z();
    double d = sqrt( ( (xyz[0] - x) * (xyz[0] - x) ) +
		     ( (xyz[1] - y) * (xyz[1] - y) ) +
		     ( (xyz[2] - z) * (xyz[2] - z) ) );
    if ( d < 0.01 ){
      if (n_found == 0){
	idx_1 = i;
	n_found += 1;
      }
      else if (n_found == 1){
	idx_2 = i;
	n_found += 1;
      }
      else
	n_found += 1;
    }// if mother is a Pi0   
  }// for all MCShowers                                                               



  //std::cout << " ***** 2 ******" << std::endl;


  size_t idxLARGE = idx_1;
  size_t idxSMALL = idx_2;

  if (mcs_h->at(idx_1).Start().E() < mcs_h->at(idx_2).Start().E() )
    { idxLARGE = idx_2; idxSMALL = idx_1; }

  auto const& mcshr1 = mcs_h->at(idxLARGE);
  auto const& mcshr2 = mcs_h->at(idxSMALL);

  _pi0_e  = mcshr1.MotherEnd().E();

  _mc_shr1_e  = mcshr1.Start().E();
  _mc_shr1_x  = mcshr1.Start().X();
  _mc_shr1_y  = mcshr1.Start().Y();
  _mc_shr1_z  = mcshr1.Start().Z();
  double mom1 = sqrt ( ( mcshr1.Start().Px() * mcshr1.Start().Px() ) +
		       ( mcshr1.Start().Py() * mcshr1.Start().Py() ) +
		       ( mcshr1.Start().Pz() * mcshr1.Start().Pz() ) );

  _mc_shr1_px = mcshr1.Start().Px() / mom1;
  _mc_shr1_py = mcshr1.Start().Py() / mom1;
  _mc_shr1_pz = mcshr1.Start().Pz() / mom1;

  _mc_shr2_e  = mcshr2.Start().E();
  _mc_shr2_x  = mcshr2.Start().X();
  _mc_shr2_y  = mcshr2.Start().Y();
  _mc_shr2_z  = mcshr2.Start().Z();
  double mom2 = sqrt ( ( mcshr2.Start().Px() * mcshr2.Start().Px() ) +
		       ( mcshr2.Start().Py() * mcshr2.Start().Py() ) +
		       ( mcshr2.Start().Pz() * mcshr2.Start().Pz() ) );

  _mc_shr2_px = mcshr2.Start().Px() / mom2;
  _mc_shr2_py = mcshr2.Start().Py() / mom2;
  _mc_shr2_pz = mcshr2.Start().Pz() / mom2;

  _mc_oangle  = mcshr1.Start().Momentum().Vect().Dot( mcshr2.Start().Momentum().Vect() );
  _mc_oangle /= mcshr1.Start().Momentum().Vect().Mag();
  _mc_oangle /= mcshr2.Start().Momentum().Vect().Mag();

  _mcradlen1 = sqrt( ( (_mc_shr1_x - _mc_vtx_x) * (_mc_shr1_x - _mc_vtx_x) ) +
		     ( (_mc_shr1_y - _mc_vtx_y) * (_mc_shr1_y - _mc_vtx_y) ) +
		     ( (_mc_shr1_z - _mc_vtx_z) * (_mc_shr1_z - _mc_vtx_z) ) );

  _mcradlen2 = sqrt( ( (_mc_shr2_x - _mc_vtx_x) * (_mc_shr2_x - _mc_vtx_x) ) +
		     ( (_mc_shr2_y - _mc_vtx_y) * (_mc_shr2_y - _mc_vtx_y) ) +
		     ( (_mc_shr2_z - _mc_vtx_z) * (_mc_shr2_z - _mc_vtx_z) ) );

  _mc_mass = sqrt( 2 * _mc_shr1_e * _mc_shr2_e * ( 1 - _mc_oangle ) );

  _n_reco_showers = shr_h->size();

  //std::cout << " ***** 3 ******" << std::endl;

  // moving on to reconstruction                                                                                                                                                                                         
  if (vtx_h->size() == 1){
    Double_t rcxyz[3] = {};
    auto const& vtx = vtx_h->at(0);
    vtx.XYZ(rcxyz);
    _rc_vtx_x = rcxyz[0];
    _rc_vtx_y = rcxyz[1];
    _rc_vtx_z = rcxyz[2];
  }
  
  std::vector<recob::Shower> reco_shower_v;
  for (size_t i=0; i < shr_h->size(); i++)
    reco_shower_v.push_back( shr_h->at(i) );
  

  std::vector<sim::MCShower> pi0_shower_v = {mcshr1,mcshr2};
  // MC <-> RC matching                                                                                                                                                                                                  
  //auto MCRCmatch_v = Match(mcshr1, mcshr2, reco_shower_v);
  auto MCRCmatch_v = Match(pi0_shower_v, reco_shower_v);
  
  for (size_t mcidx = 0; mcidx < MCRCmatch_v.size(); mcidx++) {

    auto mcshr = pi0_shower_v.at(mcidx);
    //auto mcshr = mcshr1;
    //if (mcidx == 1)
    //mcshr = mcshr2;
    
    if (MCRCmatch_v.at(mcidx) == -1) {
      _rc_shr_e = 0;
      _tree->Fill();
      continue;
    }
    
    //std::cout << " ***** 4 ******" << std::endl;
  
  auto const& rcshr = shr_h->at( MCRCmatch_v.at(mcidx) );
  
  _rc_shr_e = rcshr.Energy()[2];
  _rc_shr_x = rcshr.ShowerStart().X();
  _rc_shr_y = rcshr.ShowerStart().Y();
  _rc_shr_z = rcshr.ShowerStart().Z();
  
  double mom = sqrt( ( rcshr.Direction().X() * rcshr.Direction().X() ) +
		     ( rcshr.Direction().Y() * rcshr.Direction().Y() ) +
		     ( rcshr.Direction().Z() * rcshr.Direction().Z() ) );
  
  _rc_shr_px = rcshr.Direction().X() / mom;
  _rc_shr_py = rcshr.Direction().Y() / mom;
  _rc_shr_pz = rcshr.Direction().Z() / mom;
  
  _rcradlen = sqrt( ( (_rc_shr_x - _rc_vtx_x) * (_rc_shr_x - _rc_vtx_x) ) +
		    ( (_rc_shr_y - _rc_vtx_y) * (_rc_shr_y - _rc_vtx_y) ) +
		    ( (_rc_shr_z - _rc_vtx_z) * (_rc_shr_z - _rc_vtx_z) ) );
  
  
  _dot  = rcshr.Direction().Dot( mcshr.Start().Momentum().Vect() );
  _dot /= mcshr.Start().Momentum().Vect().Mag();
  _dot /= rcshr.Direction().Mag();
  
  _emc = mcshr.Start().E();
  
  _strt = (rcshr.ShowerStart() - mcshr.Start().DetProfile().Vect()).Mag();
  
  _tree->Fill();
  
  //std::cout << " ***** 5 ******" << std::endl;
  
  }// loop through RC showers        
  
  return;
}

void Pi0Analyzer::beginJob()
{
  // Implementation of optional member function here.
}

void Pi0Analyzer::endJob()
{
  // Implementation of optional member function here.
}

void Pi0Analyzer::Reset() {

  _n_reco_showers = 0;
  _nu_e = 0;
  _pi0_e = 0;

  _mc_vtx_x= _mc_vtx_y= _mc_vtx_z= 0;
  _rc_vtx_x= _rc_vtx_y= _rc_vtx_z= 0;

  _mc_shr1_x=  _mc_shr1_y=  _mc_shr1_z= 0;
  _mc_shr1_px= _mc_shr1_py= _mc_shr1_pz= 0;
  _mc_shr1_e= 0;

  _mc_shr2_x=  _mc_shr2_y=  _mc_shr2_z= 0;
  _mc_shr2_px= _mc_shr2_py= _mc_shr2_pz= 0;
  _mc_shr2_e= 0;

  _rc_shr_x=  _rc_shr_y=  _rc_shr_z= 0;
  _rc_shr_px= _rc_shr_py= _rc_shr_pz= 0;
  _rc_shr_e= 0;

  return;
}



//std::vector<int> Pi0Analyzer::Match(const sim::MCShower& mcs1, const sim::MCShower& mcs2,
//				    const std::vector<recob::Shower>&   shr_v) {
std::vector<int> Pi0Analyzer::Match(const std::vector<sim::MCShower>& pi0_shower_v,
				    const std::vector<recob::Shower>&   shr_v) {


  // now match to true showers                                                                                                                                                                                           
  std::vector<int> matched_indices;

  // find best matching RC shower for each MC shower               
  for (auto const& mcshr : pi0_shower_v){

    double dotmax = -1.;
    int idxmax = -1;
    
    // loop through reco showers
    for (size_t i=0; i < shr_v.size(); i++){
      auto const& rcshr = shr_v.at(i);
      
      double dot = rcshr.Direction().Dot( mcshr.Start().Momentum().Vect() );
      dot /= mcshr.Start().Momentum().Vect().Mag();
      dot /= rcshr.Direction().Mag();
      
      if (dot > dotmax) { dotmax = dot; idxmax = i; }
      
    }// for all reco indices                                                                                                                                                                                           
    
    matched_indices.push_back( idxmax );
  }
  
    /*
  dotmax = -1.;
  idxmax = 0.;
  
  // loop through reco shower
  for (size_t i=0; i < shr_v.size(); i++){
    auto const& rcshr = shr_v.at(i);
    
    double dot = rcshr.Direction().Dot( mcs2.Start().Momentum().Vect() );
    dot /= mcs2.Start().Momentum().Vect().Mag();
    dot /= rcshr.Direction().Mag();
    
    if (dot > dotmax) { dotmax = dot; idxmax = i; }
    
  }// for all reco indices                                                                                                                                                                                           
  
  matched_indices.push_back( idxmax );
  */
    
  return matched_indices;
  
}// and of function   

void Pi0Analyzer::SetTTree() {

  art::ServiceHandle<art::TFileService> tfs;
  _tree = tfs->make<TTree>("_tree","Pi0 Tree TTree");

  _tree->Branch("_run",&_run,"run/I");
  _tree->Branch("_sub",&_sub,"sub/I");
  _tree->Branch("_evt",&_evt,"evt/I");

  _tree->Branch("_n_reco_showers",&_n_reco_showers,"n_reco_showers/I");
  _tree->Branch("_event",&_event,"event/I");
  // vertex info                                                                                           
  _tree->Branch("_mc_vtx_x",&_mc_vtx_x,"mc_vtx_x/D");
  _tree->Branch("_mc_vtx_y",&_mc_vtx_y,"mc_vtx_y/D");
  _tree->Branch("_mc_vtx_z",&_mc_vtx_z,"mc_vtx_z/D");
  _tree->Branch("_rc_vtx_x",&_rc_vtx_x,"rc_vtx_x/D");
  _tree->Branch("_rc_vtx_y",&_rc_vtx_y,"rc_vtx_y/D");
  _tree->Branch("_rc_vtx_z",&_rc_vtx_z,"rc_vtx_z/D");
  // mc shower info                                          
  _tree->Branch("_mc_shr1_x",&_mc_shr1_x,"mc_shr1_x/D");
  _tree->Branch("_mc_shr1_y",&_mc_shr1_y,"mc_shr1_y/D");
  _tree->Branch("_mc_shr1_z",&_mc_shr1_z,"mc_shr1_z/D");
  _tree->Branch("_mc_shr1_e",&_mc_shr1_e,"mc_shr1_e/D");
  _tree->Branch("_mc_shr1_cont",&_mc_shr1_cont,"mc_shr1_cont/D");
  _tree->Branch("_mc_shr1_px",&_mc_shr1_px,"mc_shr1_px/D");
  _tree->Branch("_mc_shr1_py",&_mc_shr1_py,"mc_shr1_py/D");
  _tree->Branch("_mc_shr1_pz",&_mc_shr1_pz,"mc_shr1_pz/D");
  _tree->Branch("_mc_shr2_x",&_mc_shr2_x,"mc_shr2_x/D");
  _tree->Branch("_mc_shr2_y",&_mc_shr2_y,"mc_shr2_y/D");
  _tree->Branch("_mc_shr2_z",&_mc_shr2_z,"mc_shr2_z/D");
  _tree->Branch("_mc_shr2_e",&_mc_shr2_e,"mc_shr2_e/D");
  _tree->Branch("_mc_shr2_cont",&_mc_shr2_cont,"mc_shr2_cont/D");
  _tree->Branch("_mc_shr2_px",&_mc_shr2_px,"mc_shr2_px/D");
  _tree->Branch("_mc_shr2_py",&_mc_shr2_py,"mc_shr2_py/D");
  _tree->Branch("_mc_shr2_pz",&_mc_shr2_pz,"mc_shr2_pz/D");

  // reco shower info      
  _tree->Branch("_rc_shr_x",&_rc_shr_x,"rc_shr_x/D");
  _tree->Branch("_rc_shr_y",&_rc_shr_y,"rc_shr_y/D");
  _tree->Branch("_rc_shr_z",&_rc_shr_z,"rc_shr_z/D");
  _tree->Branch("_rc_shr_e",&_rc_shr_e,"rc_shr_e/D");
  _tree->Branch("_rc_shr_px",&_rc_shr_px,"rc_shr_px/D");
  _tree->Branch("_rc_shr_py",&_rc_shr_py,"rc_shr_py/D");
  _tree->Branch("_rc_shr_pz",&_rc_shr_pz,"rc_shr_pz/D");

  _tree->Branch("_rcradlen",&_rcradlen,"rcradlen/D");
  _tree->Branch("_mcradlen1",&_mcradlen1,"mcradlen1/D");
  _tree->Branch("_mcradlen2",&_mcradlen2,"mcradlen2/D");

  // reco cluster info                                             
  _tree->Branch("_ip",&_ip,"ip/D");
  _tree->Branch("_lin",&_lin,"lin/D");
  _tree->Branch("_ssv",&_ssv,"ssv/D");
  _tree->Branch("_slope",&_slope,"slope/D");
  _tree->Branch("_slopedirangle",&_slopedirangle,"slopedirangle/D");
  _tree->Branch("_hitshowerangle",&_hitshowerangle,"hitshowerangle/D");

  // MC -> RC shower comparisons                                          
  _tree->Branch("_dot",&_dot,"dot/D");
  _tree->Branch("_strt",&_strt,"strt/D");
  _tree->Branch("_emc",&_emc,"emc/D");

  // pi0 related MC information                                  
  _tree->Branch("_nu_e",&_nu_e,"nu_e/D");
  _tree->Branch("_pi_e",&_pi0_e,"pi0_e/D");
  _tree->Branch("_mc_oangle",&_mc_oangle,"mc_oangle/D");
  _tree->Branch("_mc_mass"  ,&_mc_mass  ,"_mc_mass/D" );
  _tree->Branch("_dwallmin",&_dwallmin,"dwallmin/D");

  return;
}

DEFINE_ART_MODULE(Pi0Analyzer)
