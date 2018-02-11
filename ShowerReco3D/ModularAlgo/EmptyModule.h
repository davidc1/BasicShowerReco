/**
 * \file EmptyModule.h
 *
 * \ingroup ModularAlgo
 *
 * \brief Class def header for a class EmptyModule
 *
 * @author cadams
 */

/** \addtogroup ModularAlgo

    @{*/
#ifndef EMPTYMODULE_H
#define EMPTYMODULE_H

#include <iostream>
#include "uboone/BasicShowerReco/ShowerReco3D/Base/ShowerRecoModuleBase.h"
/**
   \class ShowerRecoModuleBase
   User defined class ShowerRecoModuleBase ... these comments are used to generate
   doxygen documentation!
 */
namespace showerreco {

class EmptyModule : public ShowerRecoModuleBase {

public:

  /// Default constructor
  EmptyModule() {_name = "EmptyModule";}

  /// Default destructor
  ~EmptyModule() {}


  void do_reconstruction(const ::protoshower::ProtoShower &, Shower_t &);

private:

};

} // showerreco

#endif
/** @} */ // end of doxygen group

