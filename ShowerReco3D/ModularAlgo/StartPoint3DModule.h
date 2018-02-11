/**
 * \file StartPoint3DModule.h
 *
 * \ingroup ModularAlgo
 *
 * \brief Class def header for a class StartPoint3DModule
 *
 * @author ariana Hackenburg
 */

/** \addtogroup ModularAlgo

    @{*/
#ifndef STARTPOINT3DMODULE_H
#define STARTPOINT3DMODULE_H

#include <iostream>
#include "uboone/BasicShowerReco/ShowerReco3D/Base/ShowerRecoModuleBase.h"
/**
   \class ShowerRecoModuleBase
   User defined class ShowerRecoModuleBase ... these comments are used to generate
   doxygen documentation!
 */
namespace showerreco {

class StartPoint3DModule : public ShowerRecoModuleBase {

public:

  /// Default constructor
  StartPoint3DModule() {_name = "StartPoint3DModule"; }

  /// Default destructor
  ~StartPoint3DModule() {}

  /// Inherited/overloaded function from ShowerRecoModuleBase
  void do_reconstruction(const ::protoshower::ProtoShower &, Shower_t &);

};

} // showerreco

#endif
/** @} */ // end of doxygen group

