#ifndef DQDXMODULE_CXX
#define DQDXMODULE_CXX

#include "dQdxModule.h"

#include "math.h"
#include <algorithm>
#include <functional>

namespace showerreco {
  
  dQdxModule::dQdxModule()
  {
    _name = "dQdxModule";
    _tree = nullptr;
    _dtrunk = 0.;
  }
  
  void dQdxModule::initialize()
  {
    
    // degbugging tree
    if (_tree) delete _tree;
    _tree = new TTree(_name.c_str(),"dQdx tree");
    _tree->Branch("_pitch",&_pitch,"pitch/D");
    _tree->Branch("_dqdx",&_dqdx,"dqdx/D");
    _tree->Branch("_dqdx_v","std::vector<double>",&_dqdx_v);
    _tree->Branch("_dmax",&_dmax,"dmax/D");
    _tree->Branch("_nhits",&_nhits,"nhits/I");
    _tree->Branch("_ntot",&_ntot,"ntot/I");
    return;
  }
  
  void dQdxModule::do_reconstruction(const ::protoshower::ProtoShower & proto_shower, Shower_t & resultShower) {
    
    //if the module does not have 2D cluster info -> fail the reconstruction
    if (!proto_shower.hasCluster2D()){
      std::stringstream ss;
      ss << "Fail @ algo " << this->name() << " due to missing 2D cluster";
      throw ShowerRecoException(ss.str());
    }
    
    auto& clusters = proto_shower.clusters();
    
    // grab shower direction
    auto const& dir3D = resultShower.fDCosStart;

    //auto const& geomH = ::util::GeometryUtilities::GetME();
    //auto const* geomH = ::lar::providerFrom<util::GeometryUtilities>();
    
    // loop through planes
    for (size_t n = 0; n < clusters.size(); n++) {
      
      auto const& clus = clusters.at(n);
      
      // get the hits associated with this cluster
      auto const& hits = clus._hits;
      
      // get the plane associated with this cluster
      auto const& pl = clus._plane;
      
      if (pl != 2) continue;
      
      // grab the 2D start point of the cluster
      auto& start2D = clus._start;
      
      double f = (1 - dir3D[1]*dir3D[1] );

      _pitch = 0.001;//geomH->GetPitch(dir3D, (int)pl);
      
      _dmax = 0.;

      _nhits = 0;
      
      _dqdx_v.clear();

      _dqdx_v = std::vector<double>(3 * _dtrunk, 0.);
      
      // loop through hits and find those within some radial distance of the start point
      
      // loop over hits
      for (auto const &h : hits) {
	
	double d2D = sqrt( pow(h.w - start2D.w, 2) + pow(h.t - start2D.t, 2) );

	double d3D = d2D / f;

	if (d3D >= _dtrunk) continue;

	double qcorr = h.charge;

	
	_dqdx_v[ d3D * 3 ] += qcorr;
	
	_nhits += 1;
	
      }// loop over all hits

      std::vector<double> _dqdx_nonzero_v;
      for (auto const& dqdx : _dqdx_v)
	if (dqdx != 0) { _dqdx_nonzero_v.push_back(dqdx); }
      

      if (_dqdx_nonzero_v.size() == 0)
	_dqdx = 0.;

      else {
	std::nth_element(_dqdx_nonzero_v.begin(), _dqdx_nonzero_v.end(), _dqdx_nonzero_v.end() );
	_dqdx = _dqdx_nonzero_v[ _dqdx_nonzero_v.size()/2.] / _pitch;
      }

      _ntot = hits.size();

      resultShower.fBestdQdxPlane = pl;
      resultShower.fdQdx_v[pl] = _dqdx;
      resultShower.fBestdQdx   = _dqdx;
      
    }// for all clusters (planes)
    
    _tree->Fill();
    
    return;
  }

} //showerreco

#endif
