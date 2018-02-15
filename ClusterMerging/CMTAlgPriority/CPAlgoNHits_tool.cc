#include <iostream>
#include "uboone/BasicShowerReco/ClusterMerging/CMToolBase/CPriorityAlgoBase.h"

namespace cmtool {

  /**
     \class CPAlgoNHits
     Simple algorithm to determine priority based on # of hits.
     If # hits < set cut value by a user, returns -1.
  */
  class CPAlgoNHits : public CPriorityAlgoBase {
    
  public:
    
    /// Default constructor
    explicit CPAlgoNHits(const fhicl::ParameterSet& pset);
    
    /// Default destructor
    ~CPAlgoNHits(){};

    /**
       Core function: given the CPAN input, return a float which indicates 
       the user-defined priority for analysis.
    */
    float Priority(const ::cluster::Cluster &cluster);

    /// Setter for minimum # hits
    void SetMinHits(size_t n) { _min_hits = n; }

  protected:

    size_t _min_hits;

  };

  //----------------------------------------------
  CPAlgoNHits::CPAlgoNHits(const fhicl::ParameterSet& pset)
  //----------------------------------------------
  {
    _min_hits = 0;
  }

  //------------------------------------------------------------------------
  float CPAlgoNHits::Priority(const ::cluster::Cluster &cluster)
  //------------------------------------------------------------------------
  {
    auto nhit = cluster.size();

    return ( nhit < _min_hits ? -1 : (float)nhit );
  }
  
}