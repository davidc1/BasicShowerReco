#include <iostream>
#include "uboone/BasicShowerReco/ClusterMerging/CMToolBase/CPriorityAlgoBase.h"

namespace cmtool {

  /**
     \class CPAlgoQSum
     Simple algorithm to determine priority based on charge sum
     If charge sum < set cut value by a user, returns -1.     
  */
  class CPAlgoQSum : public CPriorityAlgoBase {
    
  public:
    
    /// Default constructor
    explicit CPAlgoQSum(const fhicl::ParameterSet& pset);
    
    /// Default destructor
    ~CPAlgoQSum(){};

    /**
       Core function: given the CPAN input, return a float which indicates 
       the user-defined priority for analysis.
    */
    float Priority(const ::cluster::Cluster &cluster);

    /// Setter for minimum charge
    void SetMinQ(double v) { _qsum_cut = v; }

  protected:

    double _qsum_cut;

  };

  //----------------------------------------------------------
  CPAlgoQSum::CPAlgoQSum(const fhicl::ParameterSet& pset)
  //----------------------------------------------------------
  {
    _qsum_cut = 0;
  }

  //------------------------------------------------------------------------------
  float CPAlgoQSum::Priority(const ::cluster::Cluster &cluster)
  //------------------------------------------------------------------------------
  {
    if(cluster._sum_charge < _qsum_cut) return -1;

    return cluster._sum_charge;
  }
  
    
}
