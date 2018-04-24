#ifndef YPLANESTARTPOINT3D_CXX
#define YPLANESTARTPOINT3D_CXX

#include <iostream>
#include "uboone/BasicShowerReco/ShowerReco3D/Base/ShowerRecoModuleBase.h"
/**
   \class ShowerRecoModuleBase
   User defined class ShowerRecoModuleBase ... these comments are used to generate
   doxygen documentation!
 */

namespace showerreco {

  class StartPointfromY2D : public ShowerRecoModuleBase {
    
  public:
    
    /// Default constructor
    StartPointfromY2D(const fhicl::ParameterSet& pset);

    /// Default destructor
    ~StartPointfromY2D() {}
    
    /// Inherited/overloaded function from ShowerRecoModuleBase
    void do_reconstruction(const ::protoshower::ProtoShower &, Shower_t &);

  private:

    double _wire2cm, _time2cm;
    
  };
  
  StartPointfromY2D::StartPointfromY2D(const fhicl::ParameterSet& pset)
  {
    _name = "StartPointfromY2D"; 
    auto const* geom = ::lar::providerFrom<geo::Geometry>();
    auto const* detp = lar::providerFrom<detinfo::DetectorPropertiesService>();
    _wire2cm = geom->WirePitch(0,1,0);
    _time2cm = detp->SamplingRate() / 1000.0 * detp->DriftVelocity( detp->Efield(), detp->Temperature() );

  }

  void StartPointfromY2D::do_reconstruction( const ::protoshower::ProtoShower & proto_shower,
					      Shower_t& resultShower)
{


  
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
    //auto const& vtx2D = util::PxPoint(pl,0,0);//geomH->Get2DPointProjection(vtx,pl);
    auto const* geom = ::lar::providerFrom<geo::Geometry>();
    auto wire = geom->WireCoordinate(vtx[1],vtx[2],geo::PlaneID(0,0,pl)) * _wire2cm;
    auto time = vtx[0];
    util::PxPoint vtx2D(pl,wire,time);
    
    auto const& start = clus._start;
    
    // 2d distance on this plane
    double d2D = sqrt( (start.w - vtx2D.w) * (start.w - vtx2D.w) + (start.t - vtx2D.t) * (start.t - vtx2D.t) );
    
    // convert to 3D distance by accounting for 3D shower direction
    auto const& dir3D   = resultShower.fDCosStart;
    double f = (1 - dir3D[1]*dir3D[1] );
    double d3D = d2D / f;
    
    // now move a distance d3D away from 3D vertex to get 3D start point
    auto strt3D = vtx + (d3D * dir3D);
    
    resultShower.fXYZStart = { strt3D[0], strt3D[1], strt3D[2]} ;
    
    pl0 = pl;
    strtpt0 = start;
    
  }// for all clusters
  
    if (pl0 != 2) {
      std::stringstream ss;
      ss << "Fail @ algo " << this->name() << " due to missing Collection-Plane 2D cluster";
      throw ShowerRecoException(ss.str());
    }
    
}

  DEFINE_ART_CLASS_TOOL(StartPointfromY2D)
} //showerreco

#endif
