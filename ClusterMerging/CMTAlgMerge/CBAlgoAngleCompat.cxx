#ifndef RECOTOOL_CBALGOANGLECOMPAT_CXX
#define RECOTOOL_CBALGOANGLECOMPAT_CXX

/**
 * \file CBAlgoAngleCompat.h
 *
 * \ingroup CMTool
 * 
 * \brief Class def header for a class CBAlgoAngleCompat
 *
 * @author david caratelli
 */

#include <iostream>
#include "uboone/BasicShowerReco/ClusterMerging/CMToolBase/CBoolAlgoBase.h"

namespace cmtool {

  /**
     \class CBAlgoAngleCompat
     An abstract fake class for merging algorithm. Having this fake class helps
     to have a better overall design of various merging for iterative approach.
     The algorithms are run through CMergeManager.
  */
  class CBAlgoAngleCompat : public CBoolAlgoBase {
    
  public:
    
    /// Default constructor
    CBAlgoAngleCompat();
    
    /// Default destructor
    virtual ~CBAlgoAngleCompat(){};
 
    /**
       Core function: given the ClusterParamsAlg input, return whether a cluster should be
       merged or not.
    */
    virtual bool Bool(const ::cluster::Cluster &cluster1,
                      const ::cluster::Cluster &cluster2);


    /// Function to reset the algorithm instance ... maybe implemented via child class
    virtual void Reset(){}

    void SetMinLargeNHits(size_t n) { _min_size = n; }
    void SetMaxAngleDiff(float a) { _max_angle_diff = a; }

  protected:

    size_t _min_size;
    float  _max_angle_diff;

    bool _flip;
    int _ctr;
  };

  //----------------------------------------
  CBAlgoAngleCompat::CBAlgoAngleCompat() : CBoolAlgoBase()
  //----------------------------------------
  {
    _flip = false;
    _ctr  = 0;
    _min_size = 0.;
    _max_angle_diff = 0.;
    // Nothing to be done in the base class
  }

  //--------------------------------------------------------
  bool CBAlgoAngleCompat::Bool(const ::cluster::Cluster &cluster1,
			       const ::cluster::Cluster &cluster2)
  //--------------------------------------------------------
  {

    // match clusters based on angle compatibility

    // skip if larger cluster (cluster 1) has less then minimum amount of hits
    if (cluster1.size() < _min_size) return false;

    if ( fabs(cluster2._angle - cluster1._angle) < _max_angle_diff ) return true;

    return false;
  }

}

#endif
