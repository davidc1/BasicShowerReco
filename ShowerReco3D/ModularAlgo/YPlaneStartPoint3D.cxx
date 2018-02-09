#ifndef YPLANESTARTPOINT3D_CXX
#define YPLANESTARTPOINT3D_CXX

#include "YPlaneStartPoint3D.h"


namespace showerreco {

  void YPlaneStartPoint3D::do_reconstruction( const ::protoshower::ProtoShower & proto_shower,
					      Shower_t& resultShower)
{


    auto const& geomH = ::util::GeometryUtilities::GetME();
    auto const* geom  = ::lar::providerFrom<geo::Geometry>();
  
    //if the module does not have 2D cluster info -> fail the reconstruction
    if (!proto_shower.hasCluster2D()) {
        std::stringstream ss;
        ss << "Fail @ algo " << this->name() << " due to missing 2D cluster";
        throw ShowerRecoException(ss.str());
    }

    if (proto_shower.hasVertex() == false){
      std::cout << "Number of vertices is not one!" << std::endl;
      return;
    }
    
    // get the proto-shower 3D vertex
    auto const& vtx = proto_shower.vertex();
    
    auto & clusters = proto_shower.clusters();

    // STEP 1:
    // identify the 2D start point on the collection plane.
    // given the shower's reconstructed 3D direction vector
    // we will find where this point lies in 3D.

    ::util::PxHit strtpt0;
    int    pl0 = 0;
    
    for (size_t i=0; i < clusters.size(); i++) {

      auto const& clus = clusters.at(i);

      // plane
      auto const& pl = clus._plane;

      if (pl != 2) continue;
      
      // project vertex onto this plane
      auto const& vtx2D = geomH->Get2DPointProjection(vtx,pl);

      auto const& start = clus._start;

      pl0 = pl;
      strtpt0 = start;
      
    }// for all clusters

    if (pl0 != 2) {
      std::stringstream ss;
      ss << "Fail @ algo " << this->name() << " due to missing Collection-Plane 2D cluster";
      throw ShowerRecoException(ss.str());
    }

    ::util::PxHit strtpt1;
    int pl1;
    bool found = false; // did we find start point coord. on 2nd plane?

    // find a plane that is not the one with the min start-vtx distance
    for (int p = 0; p < 3; p++) {

      if (p != pl0) {

	// project reconstructed 3D direction and vertex to this plane
	auto const& vtx2D   = geomH->Get2DPointProjection(vtx,p);
	auto const& dir3D   = resultShower.fDCosStart;
	auto const& slope2D = geomH->Get2DangleFrom3D(p,dir3D);

	// the 2D vertex and slope give us the line
	// on which to search for the wire coordinate
	// for the 3D start point on this plane
	// to do so we know the time coordinate
	// from the previously found point "strtpt"
	// find where this time coordinate intersects the vtx -> slope line

	// wire coordinate of the start point on this plane
	double sw = ( strtpt0.t - vtx2D.t ) / slope2D + vtx2D.w;
	double st = strtpt0.t; // start point time is the same! duh...

	strtpt1.w = sw;
	strtpt1.t = st;

	if (sw / geom->WirePitch(0) >= geom->Nwires(p) ) continue;
	if (sw / geom->WirePitch(0) < 0                ) continue;

	pl1 = p;

	found = true;
	break;
	
      }// if "another" plane
      
    }// for all planes

    if (found == false) {
      std::stringstream ss;
      ss << "Fail @ algo " << this->name() << " did not find start point on 2nd plane";
      throw ShowerRecoException(ss.str());
    }

    // now take the 2 2D start points and reconstruct the 3D one

    double sY, sZ;

    unsigned int w0 = strtpt0.w / geom->WirePitch(0);
    unsigned int w1 = strtpt1.w / geom->WirePitch(0);

    // check start/end point range
    if ( (w0 < 0) or (w0 > geom->Nwires(pl0)) or
	 (w1 < 0) or (w1 > geom->Nwires(pl1)) ) {
      std::stringstream ss;
      ss << "Fail @ algo " << this->name() << " due to wires out of range";
      throw ShowerRecoException(ss.str());
    }

    if (_verbose)
      std::cout << "pl " << pl0 << " wire : " << w0
		<< "pl " << pl1 << " wire : " << w1 << std::endl;
    
    geom->IntersectionPoint( w0, w1, pl0, pl1, 0, 0, sY, sZ);

    // check if reconstructed start point is outside of TPC volume    
    /*
    if (geomH->ContainedYZ(sY, sZ) == false){
      std::stringstream ss;
      ss << "Fail @ algo " << this->name() << " due to start point reconstructed out of TPC bouns : [Y,Z] -> ["
	 << sY << ", " << sZ << " ]" ; 
      throw ShowerRecoException(ss.str());
    }
    */

    resultShower.fXYZStart = { strtpt0.t, sY, sZ} ;

    //std::cout << "DONE " << std::endl << std::endl;
    
}

} //showerreco

#endif
