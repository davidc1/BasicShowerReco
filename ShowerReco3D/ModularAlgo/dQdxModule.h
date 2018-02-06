/**
 * \file dQdxModule.h
 *
 * \ingroup ModularAlgo
 *
 * \brief Class def header for a class dQdx2DModule
 *
 * @authors David Caratelli, Nevis, dcaratelli@nevis.columbia.edu
 */

/** \addtogroup ModularAlgo

    @{*/
#ifndef DQDXMODULE_H
#define DQDXMODULE_H
#include <iostream>
#include "ShowerRecoModuleBase.h"

//#include "AnalysisAlg/CalorimetryAlg.h"

/**
   \class dQdxModule : ShowerRecoModuleBase
   This is meant to compute the 2D dQdx along the start of the shower.
*/
namespace showerreco {
  
  class dQdxModule : ShowerRecoModuleBase {
    
  public:
    
    /// Default constructor
    dQdxModule();
    
    /// Default destructor
    ~dQdxModule() {}
    
    void do_reconstruction(const ::protoshower::ProtoShower &, Shower_t &);
    
    void initialize();
    
    // set distance along which to calculate dE/dx
    void setTrunkLength(double d) { _dtrunk = d; }

  protected:
    
    double _timetick; // sampling size in usec
    
    // distance along which to calculate dEdx
    double _dtrunk;
    
    // debugging tree
    double _dqdx;
    std::vector<double> _dqdx_v;
    std::vector<double> _dist_v;
    double _pitch;
    double _dmax;
    int    _nhits;
    int    _ntot;
    double _start_w, _start_t;

    // position-dependent response map
    std::vector< std::vector< std::vector< double > > >_responseMap;
    double _responseStep;
    
  };
  
} // showerreco

#endif
/** @} */ // end of doxygen group

