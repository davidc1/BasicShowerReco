#ifndef DQDXMODULE_CXX
#define DEDXMODULE_CXX

#include <iostream>
#include "uboone/BasicShowerReco/ShowerReco3D/Base/ShowerRecoModuleBase.h"

/**
   \class dedxModule : ShowerRecoModuleBase
   This is meant to compute the 2D dedx along the start of the shower.
*/

#include "math.h"
#include <algorithm>
#include <functional>

#include "art/Framework/Services/Optional/TFileService.h"

namespace showerreco {

  class dEdxModule : public ShowerRecoModuleBase {
    
  public:
    
    /// Default constructor
    dEdxModule(const fhicl::ParameterSet& pset);
    
    /// Default destructor
    ~dEdxModule(){};

    void configure(const fhicl::ParameterSet& pset);
    
    void do_reconstruction(const ::protoshower::ProtoShower &, Shower_t &);
    
    void initialize();
    
  protected:

    void ResetTTree();
    
    double _timetick; // sampling size in usec
    
    // distance along which to calculate dEdx
    double _dtrunk;
    
    // debugging tree
    TTree* _dedx_tree;
    double _dedx0, _dedx1, _dedx2;
    std::vector<double> _dedx0_v, _dedx1_v, _dedx2_v;
    std::vector<double> _dist0_v, _dist1_v, _dist2_v;
    int    _pl0, _pl1, _pl2;
    double _pitch0, _pitch1, _pitch2;
    int    _nhits0, _nhits1, _nhits2;
    int    _ntot0, _ntot1, _ntot2;

    // position-dependent response map
    std::vector< std::vector< std::vector< double > > >_responseMap;
    double _responseStep;
    
  };
  
  dEdxModule::dEdxModule(const fhicl::ParameterSet& pset)
  {
    _name = "dEdxModule";
    configure(pset);

    art::ServiceHandle<art::TFileService> tfs;
    _dedx_tree = tfs->make<TTree>("_dedx_tree","dE/dx TTree");

    _dedx_tree->Branch("_pl0",&_pl0,"pl0/I");
    _dedx_tree->Branch("_pitch0",&_pitch0,"pitch0/D");
    _dedx_tree->Branch("_ntot0",&_ntot0,"ntot0/I");
    _dedx_tree->Branch("_nhits0",&_nhits0,"nhits0/I");
    _dedx_tree->Branch("_dedx0",&_dedx0,"dedx0/D");
    _dedx_tree->Branch("_dedx0_v","std::vector<double>",&_dedx0_v);
    _dedx_tree->Branch("_dist0_v","std::vector<double>",&_dist0_v);

    _dedx_tree->Branch("_pl1",&_pl1,"pl1/I");
    _dedx_tree->Branch("_pitch1",&_pitch1,"pitch1/D");
    _dedx_tree->Branch("_ntot1",&_ntot1,"ntot1/I");
    _dedx_tree->Branch("_nhits1",&_nhits1,"nhits1/I");
    _dedx_tree->Branch("_dedx1",&_dedx1,"dedx1/D");
    _dedx_tree->Branch("_dedx1_v","std::vector<double>",&_dedx1_v);
    _dedx_tree->Branch("_dist1_v","std::vector<double>",&_dist1_v);

    _dedx_tree->Branch("_pl2",&_pl2,"pl2/I");
    _dedx_tree->Branch("_pitch2",&_pitch2,"pitch2/D");
    _dedx_tree->Branch("_ntot2",&_ntot2,"ntot2/I");
    _dedx_tree->Branch("_nhits2",&_nhits2,"nhits2/I");
    _dedx_tree->Branch("_dedx2",&_dedx2,"dedx2/D");
    _dedx_tree->Branch("_dedx2_v","std::vector<double>",&_dedx2_v);
    _dedx_tree->Branch("_dist2_v","std::vector<double>",&_dist2_v);

  }

  void dEdxModule::configure(const fhicl::ParameterSet& pset)
  {
    _dtrunk = pset.get<double>("dtrunk");
    _verbose   = pset.get<bool>("verbose",false);
  }
  
  void dEdxModule::initialize()
  {
    
    return;
  }
  
  void dEdxModule::do_reconstruction(const ::protoshower::ProtoShower & proto_shower, Shower_t & resultShower) {
    
    //if the module does not have 2D cluster info -> fail the reconstruction
    if (!proto_shower.hasCluster2D()){
      std::stringstream ss;
      ss << "Fail @ algo " << this->name() << " due to missing 2D cluster";
      throw ShowerRecoException(ss.str());
    }
    
    auto& clusters = proto_shower.clusters();
    
    // grab shower direction
    auto const& dir3D = resultShower.fDCosStart;

    std::cout << "3D shower direction : " << dir3D[0] << ", " << dir3D[1] << ", " << dir3D[2] << std::endl;

    auto const& geomH = ::util::GeometryUtilities();

    ResetTTree();
    
    // loop through planes
    for (size_t n = 0; n < clusters.size(); n++) {
      
      auto const& clus = clusters.at(n);
      
      // get the hits associated with this cluster
      auto const& hits = clus._hits;
      
      // get the plane associated with this cluster
      auto const& pl = clus._plane;
      
      // get start point on pllane
      auto& start2D = clus._start;
      
      std::cout << std::endl << "PLANE : " << pl << std::endl;

      art::ServiceHandle<geo::Geometry> geom;
      const geo::WireGeo& wire = geom->TPC().Plane(pl).MiddleWire();

      TVector3 wireunitperp = wire.Direction();//(wire.GetStart()-wire.GetEnd()).Unit();
      // rotate by 90 degrees around x
      TVector3 wireunit = {wireunitperp[0], -wireunitperp[2], wireunitperp[1]}; 
      std::cout << "wire unit on plane : " << pl << " is " << wireunit[0] << ", " << wireunit[1] << ", " << wireunit[2] << std::endl;
      double cosPlane = cos(wireunit.Angle(dir3D));

      std::vector<double> dedx_v;
      double dedx;
      int nhits = 0;
      double pitch = 0.3 / fabs(cosPlane);
      dedx_v = std::vector<double>(3 * _dtrunk, 0.);      
      std::cout << " dEdx Module : pitch = " << pitch << " from function" <<  std::endl;      

      // loop through hits and find those within some radial distance of the start point
      // loop over hits
      for (auto const &h : hits) {
	
	double d2D = sqrt( pow(h.w - start2D.w, 2) + pow(h.t - start2D.t, 2) );
	
	double d3D = d2D / cosPlane;
	
	if (d3D >= _dtrunk) continue;
	
	double qcorr = h.charge;

	dedx_v[ d3D * 3 ] += qcorr;
	nhits += 1;
	
	
      }// loop over all hits

      std::vector<double> _dedx_nonzero_v;
      for (auto const& dedx : dedx_v) {
	if (dedx != 0) {
	  _dedx_nonzero_v.push_back(dedx); 
	  std::cout << "dedx Module : \t dedx = " << dedx / pitch << std::endl;
	}
      }// for all dedx values

      if (_dedx_nonzero_v.size() == 0)
	dedx = 0.;

      else {
	std::nth_element(_dedx_nonzero_v.begin(), _dedx_nonzero_v.end(), _dedx_nonzero_v.end() );
	dedx = _dedx_nonzero_v[ _dedx_nonzero_v.size()/2.] / pitch;
      }

      std::cout << "dedx Module : Final dEdx = " << dedx << std::endl;


      if (pl == 0) {
	_pitch0 = pitch;
	_nhits0 = nhits;
	_dedx0_v = dedx_v;
	_dedx0 = dedx;
	_ntot0 = hits.size();
      }
      if (pl == 1) {
	_pitch1 = pitch;
	_nhits1 = nhits;
	_dedx1_v = dedx_v;
	_dedx1 = dedx;
	_ntot1 = hits.size();
      }
      if (pl == 2) {
	_pitch2 = pitch;
	_nhits2 = nhits;
	_dedx2_v = dedx_v;
	_dedx2 = dedx;
	_ntot2 = hits.size();
      }

      resultShower.fBestdEdxPlane = pl;
      resultShower.fdEdx_v[pl] = dedx;
      if (pl == 2)
	resultShower.fBestdEdx   = dedx;
      
    }// for all clusters (planes)

    _dedx_tree->Fill();
    
    return;
  }

  void dEdxModule::ResetTTree() {

    _dedx0 = 0;
    _dedx1 = 0;
    _dedx2 = 0;
    _dedx0_v.clear();
    _dedx1_v.clear();
    _dedx2_v.clear();
    _dist0_v.clear();
    _dist1_v.clear();
    _dist2_v.clear();
    _pl0 = 0;
    _pl1 = 0;
    _pl2 = 0;
    _pitch0 = 0;
    _pitch1 = 0;
    _pitch2 = 0;
    _nhits0 = 0;
    _nhits1 = 0;
    _nhits2 = 0;
    _ntot0 = 0;
    _ntot1 = 0;
    _ntot2 = 0;

    return;
  }

  DEFINE_ART_CLASS_TOOL(dEdxModule)
} //showerreco

#endif
