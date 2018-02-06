#ifndef PROTOSHOWER_H
#define PROTOSHOWER_H

#include <TVector3.h>
#include <vector>

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "Cluster2D.h"

namespace protoshower {


class ProtoShower
{
  friend class ProtoShowerHelper;

public:
  ProtoShower() {};
  ~ProtoShower() {};

  const std::vector<Double_t*> & vertexes() const {return _vertexes;}
  const std::vector<::cluster2d::Cluster2D> & clusters() const { return _clusters; }

  // getters
  bool hasCluster2D() const {return _hasCluster2D;}
  bool hasCluster3D() const {return _hasCluster3D;}
  bool hasVertex()    const {return _hasVertex;}
  // setters
  void hasCluster2D(bool on) { _hasCluster2D = on; }
  void hasCluster3D(bool on) { _hasCluster3D = on; }
  void hasVertex   (bool on) { _hasVertex    = on; }

  // list of 3D vertices associated to this ProtoShower
  std::vector<Double_t*>  _vertexes;
  
  // list 2D clusters
  std::vector<::cluster2d::Cluster2D> _clusters;

protected:

  bool _hasCluster2D;
  bool _hasCluster3D;
  bool _hasVertex;



  // Not sure what to do with vertexes yet

};

}

#endif
