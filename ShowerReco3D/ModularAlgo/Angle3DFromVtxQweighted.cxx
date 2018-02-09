#ifndef ANGLE3DFROMVERTEXQWEIGHTED_CXX
#define ANGLE3DFROMVERTEXQWEIGHTED_CXX

#include "Angle3DFromVtxQweighted.h"

#include <math.h>
#include <sstream>
#include <algorithm>

namespace showerreco {
  
  void Angle3DFromVtxQweighted::do_reconstruction(const ::protoshower::ProtoShower & proto_shower,
					 Shower_t& resultShower) {

    //if the module does not have 2D cluster info -> fail the reconstruction
    if (!proto_shower.hasVertex()){
      std::stringstream ss;
      ss << "Fail @ algo " << this->name() << " due to missing Vertex";
      throw ShowerRecoException(ss.str());
    }

    //if the module does not have 2D cluster info -> fail the reconstruction
    if (!proto_shower.hasCluster2D()){
      std::stringstream ss;
      ss << "Fail @ algo " << this->name() << " due to missing 2D cluster";
      throw ShowerRecoException(ss.str());
    }

    auto const& geomH = ::util::GeometryUtilities::GetME();
    
    if (proto_shower.hasVertex() == false){
      std::cout << "Number of vertices is not one!" << std::endl;
      return;
    }
    
    // get the proto-shower 3D vertex
    auto const& vtx = proto_shower.vertex();

    std::vector<double> dir3D = {0,0,0};

    auto & clusters = proto_shower.clusters();

    // planes with largest number of hits used to get 3D direction
    std::vector<int> planeHits(3,0);
    std::vector<::util::PxPoint> planeDir(3);
    
    // we want an energy for each plane
    for (size_t n = 0; n < clusters.size(); n++) {
      
      // get the hits associated with this cluster
      auto const& hits = clusters.at(n)._hits;
      
      // get the plane associated with this cluster
      auto const& pl = clusters.at(n)._plane;

      // project vertex onto this plane
      auto const& vtx2D = geomH->Get2DPointProjection(vtx,pl);

      // get the charge-averaged 2D vector pointing from vtx in shower direction
      ::util::PxPoint weightedDir;
      weightedDir.w = 0;
      weightedDir.t = 0;
      double Qtot = 0;
      for (auto const& hit : hits){
	weightedDir.w += (hit.w - vtx2D.w) * hit.charge;
	weightedDir.t += (hit.t - vtx2D.t) * hit.charge;
	Qtot += hit.charge;
      }

      weightedDir.w /= Qtot;
      weightedDir.t /= Qtot;

      planeHits[pl] = (int)hits.size();
      planeDir[pl]  = weightedDir;

    }// for all planes

    int pl_max = 0;
    int n_max  = 0;
    int pl_mid = 0;
    int pl_min = 0;
    int n_min  = 1000;
    for (size_t pl=0; pl < planeHits.size(); pl++){
      if (planeHits[pl] > n_max){
	pl_max = pl;
	n_max  = planeHits[pl];
      }
      if (planeHits[pl] < n_min){
	pl_min = pl;
	n_min  = planeHits[pl];
      }
    }
    // find the medium plane
    std::vector<int> planes = {0,1,2};
    planes.erase( std::find( planes.begin(), planes.end(), pl_max) );
    planes.erase( std::find( planes.begin(), planes.end(), pl_min) );
    pl_mid = planes[0];

    double slope_max, slope_mid;
    double angle_max, angle_mid;
    slope_max = planeDir[pl_max].t / planeDir[pl_max].w;
    angle_max = atan(slope_max);
    angle_max = atan2( planeDir[pl_max].t , planeDir[pl_max].w );
    slope_mid = planeDir[pl_mid].t / planeDir[pl_mid].w;
    angle_mid = atan(slope_mid);
    angle_mid = atan2( planeDir[pl_mid].t , planeDir[pl_mid].w );
    
    double theta, phi;
    geomH->Get3DaxisN(pl_max, pl_mid,
		      angle_max, angle_mid,
		      phi, theta);
    
    resultShower.fDCosStart[1] = sin(theta);
    resultShower.fDCosStart[0] = cos(theta) * sin(phi);
    resultShower.fDCosStart[2] = cos(theta) * cos(phi);
    
    return;
  }
  
  
} //showerreco

#endif
