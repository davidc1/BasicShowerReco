/**
 * \file ShowerRecoManager.h
 *
 * \ingroup ShowerReco3D
 *
 * \brief Class def header for a class ShowerRecoManager
 *
 * @author kazuhiro
 */

/** \addtogroup ShowerReco3D

    @{*/
#ifndef SHOWERRECO_PI0RECOALGORITHM_H
#define SHOWERRECO_PI0RECOALGORITHM_H

#include "uboone/BasicShowerReco/ShowerReco3D/Base/ShowerRecoManager.h"

// include all individual algorithms
#include "uboone/BasicShowerReco/ShowerReco3D/ModularAlgo/Angle3DFormula.h"
#include "uboone/BasicShowerReco/ShowerReco3D/ModularAlgo/Angle3DFromVtx.h"
#include "uboone/BasicShowerReco/ShowerReco3D/ModularAlgo/Angle3DFromVtxQweighted.h"
#include "uboone/BasicShowerReco/ShowerReco3D/ModularAlgo/EmptyModule.h"
#include "uboone/BasicShowerReco/ShowerReco3D/ModularAlgo/FillLength.h"
#include "uboone/BasicShowerReco/ShowerReco3D/ModularAlgo/FilterPFPart.h"
#include "uboone/BasicShowerReco/ShowerReco3D/ModularAlgo/FilterShowers.h"
#include "uboone/BasicShowerReco/ShowerReco3D/ModularAlgo/LinearEnergy.h"
#include "uboone/BasicShowerReco/ShowerReco3D/ModularAlgo/StartPoint3DModule.h"
#include "uboone/BasicShowerReco/ShowerReco3D/ModularAlgo/YPlaneStartPoint3D.h"
#include "uboone/BasicShowerReco/ShowerReco3D/ModularAlgo/dEdxFromdQdx.h"
#include "uboone/BasicShowerReco/ShowerReco3D/ModularAlgo/dQdxModule.h"

namespace showerreco {

  /**
     Describe class here
  */
  class Pi0RecoAlgorithm : public ShowerRecoManager {
    
  public:

    /// Default constructor
    Pi0RecoAlgorithm();
    
    /// Default destructor
    ~Pi0RecoAlgorithm() {}

  };
}

#endif
/** @} */ // end of doxygen group

