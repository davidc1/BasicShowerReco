#ifndef SHOWERRECO_SHOWERRECOMANAGER_CXX
#define SHOWERRECO_SHOWERRECOMANAGER_CXX

#include "ShowerRecoManager.h"
#include <iomanip>

namespace showerreco {

ShowerRecoManager::ShowerRecoManager()
{

}

void ShowerRecoManager::Initialize()
{
  for (auto alg : _alg_v)
    alg->Initialize();

  return;
}

void ShowerRecoManager::Reset()
{
  for (auto alg : _alg_v) alg->Reset();

  _proto_showers.clear();

}

void ShowerRecoManager::Reconstruct (std::vector<showerreco::Shower_t>& showers)
{

  showers.clear();
  showers.reserve(_proto_showers.size());

  for (auto const& proto_shower : _proto_showers) {


    for (auto& shower_alg : _alg_v) {
      auto this_shower = shower_alg->RecoOneShower(proto_shower);
      showers.emplace_back(this_shower);
    }// for all shower reco algorithms? not modular algos...

    // why is this a for-loop? DC very confused here...
    for (auto& shower_ana : _ana_v)
      shower_ana->Analyze(showers.back(), proto_shower);

  }// for all pfparticle proto-showers

  // Check that the showers reconstructed are the same length as the proto_showers vector
  if (showers.size() != _proto_showers.size()) {
    throw ShowerRecoException("ERROR: number of reconstructed showers doesn't match input list!!");
  }

  return;
}

void ShowerRecoManager::Finalize(TFile* fout) {
  if (!fout)
    return;
  fout->cd();

  // loop through algorithms
  for (auto& alg : _alg_v)
    alg->Finalize(fout);

  // loop through ana modules
  for (auto& ana : _ana_v) {
    auto tree = ana->GetTree();
    if (tree) tree->Write();
  }

  return;
}

}

#endif
