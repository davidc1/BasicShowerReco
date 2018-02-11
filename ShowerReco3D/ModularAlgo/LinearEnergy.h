/**
 * \file LinearEnergy.h
 *
 * \ingroup ModularAlgo
 *
 * \brief Class def header for a class LinearEnergy
 *
 * @author david caratelli
 */

/** \addtogroup ModularAlgo

    @{*/
#ifndef LINEARENERGY_H
#define LINEARENERGY_H

#include <iostream>

#include "uboone/BasicShowerReco/ShowerReco3D/Base/ShowerRecoModuleBase.h"


/**
   \class ShowerRecoModuleBase
   User defined class ShowerRecoModuleBase ... these comments are used to generate
   doxygen documentation!
 */
namespace showerreco {

class LinearEnergy : public ShowerRecoModuleBase {

public:

    /// Default constructor
    LinearEnergy();

    /// Default destructor
    ~LinearEnergy() {}

    void do_reconstruction(const ::protoshower::ProtoShower &, Shower_t &);

    void initialize();

 private:
    
};

} // showerreco

#endif
/** @} */ // end of doxygen group
