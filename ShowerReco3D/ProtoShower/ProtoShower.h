#ifndef PROTOSHOWER_H
#define PROTOSHOWER_H

#include <TVector3.h>
#include <vector>

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"

namespace protoshower {


class ProtoShower
{
  friend class ProtoShowerHelper;

public:
  ProtoShower() {};
  ~ProtoShower() {};

  const std::vector<TVector3> & vertexes() const {return _vertexes;}

  // getters
  bool hasCluster2D() const {return _hasCluster2D;}
  bool hasCluster3D() const {return _hasCluster3D;}
  bool hasVertex()    const {return _hasVertex;}
  // setters
  void hasCluster2D(bool on) { _hasCluster2D = on; }
  void hasCluster3D(bool on) { _hasCluster3D = on; }
  void hasVertex   (bool on) { _hasVertex    = on; }

  // list of 3D vertices associated to this ProtoShower
  std::vector<TVector3>  _vertexes;

protected:

  bool _hasCluster2D;
  bool _hasCluster3D;
  bool _hasVertex;



  // Not sure what to do with vertexes yet

};

}

#endif
