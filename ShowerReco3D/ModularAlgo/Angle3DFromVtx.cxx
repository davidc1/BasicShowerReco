#ifndef ANGLE3DFROMVERTEX_CXX
#define ANGLE3DFROMVERTEX_CXX

#include "Angle3DFromVtx.h"
#include <math.h>
#include <sstream>

namespace showerreco {
  
  void Angle3DFromVtx::do_reconstruction(const ::protoshower::ProtoShower & proto_shower,
					 Shower_t& resultShower) {

    //if the module does not have 2D cluster info -> fail the reconstruction
    if (!proto_shower.hasVertex()){
      std::stringstream ss;
      ss << "Fail @ algo " << this->name() << " due to missing Vertex";
      throw ShowerRecoException(ss.str());
    }
    
    // get the 3D start point of the shower
    auto const& start3D = resultShower.fXYZStart;


    if (proto_shower.hasVertex() == false){
      std::cout << "Number of vertices is not one!" << std::endl;
      return;
    }
    // get the proto-shower 3D vertex
    auto const& vtx = proto_shower.vertex();

    std::vector<double> dir3D = {0,0,0};

    dir3D[0] = start3D[0] - vtx[0];
    dir3D[1] = start3D[1] - vtx[1];
    dir3D[2] = start3D[2] - vtx[2];

    // normalize
    double mag = sqrt( dir3D[0]*dir3D[0] + dir3D[1]*dir3D[1] + dir3D[2]*dir3D[2] );
    dir3D[0] /= mag;
    dir3D[1] /= mag;
    dir3D[2] /= mag;
    
    // projet to fDCosStart values
    resultShower.fDCosStart[0] = dir3D[0];
    resultShower.fDCosStart[1] = dir3D[1];
    resultShower.fDCosStart[2] = dir3D[2];


    return;
  }
  
  
} //showerreco

#endif
