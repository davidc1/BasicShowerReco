#ifndef PROTOSHOWERALGBASE_CXX
#define PROTOSHOWERALGBASE_CXX

#include "ProtoShowerAlgBase.h"

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardata/Utilities/GeometryUtilities.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

namespace protoshower {
  
  ProtoShowerAlgBase::ProtoShowerAlgBase()
  {
    
    _name = "ProtoShowerAlgBase";

    // get detector specific properties
    auto const* geom = ::lar::providerFrom<geo::Geometry>();
    auto const* detp = lar::providerFrom<detinfo::DetectorPropertiesService>();
    _wire2cm = geom->WirePitch(0,1,0);
    _time2cm = detp->SamplingRate() / 1000.0 * detp->DriftVelocity( detp->Efield(), detp->Temperature() );
    
  }
  
  ::cluster2d::Cluster2D ProtoShowerAlgBase::MakeCluster2D( const art::Ptr<recob::Cluster>& clus, 
							  const std::vector< art::Ptr<recob::Hit> >& hit_v) 
  {
    
    ::cluster2d::Cluster2D clus2d;
    clus2d.Reset();
    
    for (auto const hit : hit_v) {
      // create PxHit
      ::util::PxHit hit2d( clus->Plane().Plane,
			   hit->WireID().Wire * _wire2cm, hit->PeakTime() * _time2cm,
			   hit->Integral(), hit->Integral(), hit->PeakAmplitude() );
      
      clus2d._hits.push_back( hit2d );
    }// for all hits
    
    clus2d._plane = clus->Plane().Plane;

    auto const& sw = clus->StartWire() * _wire2cm;
    auto const& ew = clus->EndWire()   * _wire2cm;
    auto const& st = clus->StartTick() * _time2cm;
    auto const& et = clus->EndTick()   * _time2cm;
    
    clus2d._start = ::util::PxHit(clus->Plane().Plane, sw, st, 0., 0., 0.);
    clus2d._end   = ::util::PxHit(clus->Plane().Plane, ew, et, 0., 0., 0.);
    
    clus2d._angle_2d = clus->StartAngle();
    
    return clus2d;
    
  }

}

#endif
