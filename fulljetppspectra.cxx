// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// fulljet spectra task with Run3 pp dataset

/// \author Archita Rani Dash (archita.rani.dash@cern.ch), University of Muenster

#include <bitset>

#include "TH1F.h"
#include "TTree.h"

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/HistogramRegistry.h"

#include "Common/DataModel/EventSelection.h"

#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/EMCALClusters.h"
#include "PWGJE/Core/JetFinder.h"
#include "PWGJE/Core/FastJetUtilities.h"

#include "EventFiltering/filterTables.h"

using namespace o2;
using namespace o2::analysis;
using namespace o2::framework;
using namespace o2::framework::expressions;
//using selectedClusters = o2::soa::Filtered<o2::aod::EMCALClusters>;
//using filteredJets = o2::soa::Filtered<o2::aod::FullJets>;

struct fulljetppspectra {
  using selectedClusters = o2::soa::Filtered<o2::aod::EMCALClusters>;
  using fulljetsonly = soa::Join<aod::FullJets, aod::FullJetConstituents>;
  using neutraljetsonly = soa::Join<aod::NeutralJets, aod::NeutralJetConstituents>;
  using collisionWithTrigger = soa::Join<aod::Collisions, aod::EvSels, aod::FullJetFilters>::iterator;

  //Define histogram registries here
  HistogramRegistry registry{"registry",
                            {{"h_full_jet_pt", "jet pT;#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 50.}}}},
                             {"h_full_jet_eta", "jet #eta;#eta_{jet};entries", {HistType::kTH1F, {{30, -1.5, 1.5}}}},
                             {"h_full_jet_phi", "jet #phi;#phi{jet};entries", {HistType::kTH1F, {{140, -7.0, 7.}}}},
                             {"h_full_jet_ntracks", "jet N tracks;N_{jet tracks};entries", {HistType::kTH1F, {{40, -0.5, 39.5}}}},
                             {"h_full_jet_nclusters", "jet N clusters;N_{jet clusters};entries", {HistType::kTH1F, {{40, -0.5, 39.5}}}},
                            }};

 //Define Configurables here
  Configurable<float> f_jetPtMin{"f_jetPtMin", 0.0, "minimum jet pT cut"};
  Configurable<float> f_jetR{"f_jetR", 0.4, "jet resolution parameter that you have triggered on"};
  Configurable<float> f_PhiEmcalOrDcal{"f_PhiEmcalOrDcal", 4, "if cluster phi is less than this value, count it to be EMCAL"};
  Configurable<float> f_minClusterTime{"f_minClusterTime", -999, "Min. cluster time for gamma trigger (ns)"};
  Configurable<float> f_maxClusterTime{"f_maxClusterTime", 999, "Max. cluster time for gamma trigger (ns)"};

  Configurable<bool> b_JetsInEmcalOnly{"b_JetsInEmcalOnly", true, "fill histograms only for jets inside the EMCAL"};
  Configurable<bool> b_IgnoreEmcalFlag{"b_IgnoreEmcalFlag", false, "ignore the EMCAL live flag check"};
  Configurable<bool> b_DoFiducialCut{"b_DoFiducialCut", false, "do a fiducial cut on jets to check if they are in the emcal"};
  Configurable<bool> b_RejectExoticClusters{"b_RejectExoticClusters", true, "Reject exotic clusters"};
  Configurable<int> f_GammaObservable{"f_gammaObservable", 0, "Observable for gamma trigger (0 - energy, 1 - pt)"};

  Configurable<float> cfgVertexCut{"cfgVertexCut", 10.0f, "Accepted z-vertex range"};
  Configurable<std::string> mClusterDefinition{"clusterDefinition", "kV3Default", "cluster definition to be selected, e.g. V3Default"};

  Filter jetCuts = aod::jet::pt > f_jetPtMin&& aod::jet::r == nround(f_jetR.node() * 100.0f);

  void init(InitContext const&) {}

  void processFullJetSpectra(collisionWithTrigger const& collision, fulljetsonly const& jets , aod::Tracks const& tracks, selectedClusters const& clusters)
  { //loop over the collisions with trigger once

    registry.fill(HIST("h_full_jet_pt"), jet.pt());
    registry.fill(HIST("h_full_jet_eta"), jet.eta());
    registry.fill(HIST("h_full_jet_phi"), jet.phi());
    registry.fill(HIST("h_full_jet_ntracks"), jet.tracks().size());
    registry.fill(HIST("h_full_jet_nclusters"), jet.clusters().size());

  }

  PROCESS_SWITCH(fulljetppspectra, processDataFullSubstructure, "jet substructure full jets", false);

  void processDummy(aod::Tracks const& track) {}
  PROCESS_SWITCH(fulljetppspectra, processDummy, "Dummy process function turned on by default", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc){
  return WorkflowSpec{
    adaptAnalysisTask<fulljetppspectra>(cfgc, TaskName{"full-jet-pp-spectra"} )
  };
}
