/**
 * \file ProtoShowerCMTool.h
 *
 * \ingroup ProtoShower
 *
 * \brief Class def header for a class ProtoShowerCMTool
 *
 * @author david caratelli
 */

/** \addtogroup ProtoShower

    @{*/
#ifndef PROTOSHOWER_PROTOSHOWERALGBASE_H
#define PROTOSHOWER_PROTOSHOWERALGBASE_H

#include <iostream>

#include "ProtoShowerAlgBase.h"

/**
   \class ProtoShowerCMTool
   User defined class ProtoShowerCMTool ... these comments are used to generate
   doxygen documentation!
 */

namespace protoshower {

  class ProtoShowerCMTool : ProtoShowerAlgBase {

public:

  /// Default constructor
  ProtoShowerCMTool() { _name = "ProtoShowerCMTool"; }

  /// Default destructor
  virtual ~ProtoShowerCMTool() {}

  void GenerateProtoShowers(::art::Event & e,
			    const std::string& fPFPproducer,
			    std::vector<protoshower::ProtoShower> & proto_shower_v);
  

  std::string name() { return _name; }

protected:

  std::string _name;

  // cluster producer name
  std::string fClusProducer;
  // vertex producer name
  std::string fVertexProducer;
  // hit producer name
  std::string fHitProducer;

};

}// namespace

#endif
/** @} */ // end of doxygen group

