/**
 * \file ProtoShowerAlgBase.h
 *
 * \ingroup ProtoShower
 *
 * \brief Class def header for a class ProtoShowerAlgBase
 *
 * @author david caratelli
 */

/** \addtogroup ProtoShower

    @{*/
#ifndef PROTOSHOWERALGBASE_H
#define PROTOSHOWERALGBASE_H

#include <iostream>

#include "art/Framework/Principal/Event.h"

#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardata/Utilities/AssociationUtil.h"

#include "ProtoShower.h"

/**
   \class ProtoShowerAlgBase
   User defined class ProtoShowerAlgBase ... these comments are used to generate
   doxygen documentation!
 */

namespace protoshower {

class ProtoShowerAlgBase {

public:

  /// Default constructor
  ProtoShowerAlgBase() { _name = "ProtoShowerAlgBase"; }

  /// Default destructor
  virtual ~ProtoShowerAlgBase() {}

  virtual void GenerateProtoShower(::art::Event & e,
				   const std::unique_ptr<std::vector<recob::PFParticle> >& pfp_v,
                                   const size_t proto_shower_pfpart,
                                   protoshower::ProtoShower & proto_shower) = 0;


  std::string name() { return _name; }

protected:

  std::string _name;

};

}// namespace

#endif
/** @} */ // end of doxygen group

