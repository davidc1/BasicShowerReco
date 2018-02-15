/**
 * \file CMAlgoBase.h
 *
 * \ingroup CMTool
 * 
 * \brief Class def header for a class CMAlgoBase
 *
 * @author kazuhiro
 */

// art TOOLS
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"

/** \addtogroup CMTool

    @{*/
#ifndef RECOTOOL_CMALGOBASE_H
#define RECOTOOL_CMALGOBASE_H

#include <iostream>
#include "TFile.h"
#include "Cluster.h"

#include "fhiclcpp/ParameterSet.h"

namespace cmtool {

  /**
     \class CMAlgoBase
     An abstract base class for merging/mathcing algorithm. Having this base class helps
     to have a better overall design of various merging for iterative approach.
     The algorithms are run through CMergeManager.
  */
  class CMAlgoBase {
    
  public:
    
    virtual ~CMAlgoBase() noexcept = default;

    void configure(const fhicl::ParameterSet& pset){};

    /// Function to reset the algorithm instance called within CMergeManager/CMatchManager's Reset() ... maybe implemented via child class
    virtual void Reset(){};

    /**
       Optional function: called after Bool() function is called for all possible cluster
       pairs by CMergeManager/CMatchManager IFF run with verbosity level kPerIteration. 
       Maybe useful for debugging.       
     */
    virtual void Report()
    {return;}

    /// Setter function for an output plot TFile pointer
    void SetAnaFile(TFile* fout) { _fout = fout; }

    /// Setter function for verbosity
    void SetVerbose(bool doit=true) { _verbose = doit; }

  protected:

    /// TFile pointer to an output file
    TFile* _fout;

    /// Boolean to choose verbose mode. Turned on if CMergeManager/CMatchManager's verbosity level is >= kPerMerging
    bool _verbose;

  };

}

#endif
/** @} */ // end of doxygen group 

