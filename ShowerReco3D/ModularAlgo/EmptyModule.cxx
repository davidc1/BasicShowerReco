#ifndef EMPTYMODULE_CXX
#define EMPTYMODULE_CXX

#include "EmptyModule.h"

namespace showerreco {

void EmptyModule::do_reconstruction(const ::protoshower::ProtoShower & proto_shower,
                                    Shower_t& resultShower) {

  // This function takes the shower cluster set and computes the best fit 3D axis
  // and then assigns it to the shower

}


} //showerreco

#endif
