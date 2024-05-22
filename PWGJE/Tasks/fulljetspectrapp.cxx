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

// FullJet Spectra in pp
//
/// \author Archita Rani Dash <archita.rani.dash@cern.ch>
#include <vector>
#include <iostream>
#include <utility>

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "Framework/RunningWorkflowInfo.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "PWGHF/Core/HfHelper.h"

#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/EMCALClusters.h"
#include "PWGJE/Core/JetFinder.h"
// #include "PWGJE/Core/FastJetUtilitites.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/Core/JetFindingUtilities.h"

#include "EventFiltering/filterTables.h"

using namespace std;
using namespace o2;
using namespace o2::analysis;
using namespace o2::framework;
using namespace o2::framework::expressions;

// template <typename JetTableData, typename JetConstituentTableData, typename JetTableMCD, typename JetConstituentTableMCD, typename JetMatchingTableMCDMCP, typename JetTableMCP, typename JetConstituentTableMCP, typename JetMatchingTableMCPMCD>
//haven't included JetTableMCDWeighted and JetTableMCPWeighted yet. TO be included once I check without the weights
struct FullJetSpectrapp {

  HistogramRegistry registry;

  // Event configurables
  Configurable<float> VertexZCut{"VertexZCut", 10.0f, "Accepted z-vertex range"};
  Configurable<float> centralityMin{"centralityMin", -999.0, "minimum centrality"};
  Configurable<float> centralityMax{"centralityMax", 999.0, "maximum centrality"};

  // Jet configurables
  Configurable<double> JetRadii{"JetRadii", 0.4, "jet resolution parameters"};
  Configurable<float> jetpTMin{"jetpTMin", 0., "minimum jet pT"};
  Configurable<float> jetpTMax{"jetpTMax", 350., "maximum jet pT"};
  Configurable<float> jetEtaMin{"jetEtaMin", -1.0, "minimum jet eta"};
  Configurable<float> jetEtaMax{"jetEtaMax", 1.0, "maximum jet eta"};
  Configurable<float> jetPhiMin{"jetPhiMin", 0., "minimum jet phi"};
  Configurable<float> jetPhiMax{"jetPhiMax", 7., "maximum jet phi"};
  Configurable<float> jetAreaFractionMin{"jetAreaFractionMin", -99.0, "used to make a cut on the jet areas"};
  Configurable<float> leadingConstituentPtMin{"leadingConstituentPtMin", -99.0, "minimum pT selection on jet constituent"};


  // Track configurables
  Configurable<float> trackpTMin{"trackpTMin", 0.15, "minimum track pT"};
  Configurable<float> trackpTMax{"trackpTMax", 350., "maximum track pT"};
  Configurable<float> trackEtaMin{"trackEtaMin", -1.0, "minimum track eta"};
  Configurable<float> trackEtaMax{"trackEtaMax", 1.0, "maximum track eta"};
  Configurable<float> trackPhiMin{"trackPhiMin", 0., "minimum track phi"};
  Configurable<float> trackPhiMax{"trackPhiMax", 7., "maximum track phi"};
  Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "set track selections"};
  Configurable<std::string> eventSelections{"eventSelections", "sel8", "choose event selection"};
  Configurable<std::string> particleSelections{"particleSelections", "PhysicalPrimary", "set particle selections"};

  // Cluster configurables
  Configurable<std::string> clusterDefinitionS{"clusterDefinition", "kV3Default", "cluster definition to be selected, e.g. V3Default"};
  Configurable<float> clusterEtaMin{"clusterEtaMin", -0.7, "minimum cluster eta"};
  Configurable<float> clusterEtaMax{"clusterEtaMax", 0.7, "maximum cluster eta"};
  Configurable<float> clusterPhiMin{"clusterPhiMin", 1.39, "minimum cluster phi"};
  Configurable<float> clusterPhiMax{"clusterPhiMax", 3.27, "maximum cluster phi"};
  Configurable<float> clusterEnergyMin{"clusterEnergyMin", 0.5, "minimum cluster energy in EMCAL (GeV)"};
  Configurable<float> clusterTimeMin{"clusterTimeMin", -25., "minimum cluster time (ns)"};
  Configurable<float> clusterTimeMax{"clusterTimeMax", 25., "maximum cluster time (ns)"};
  Configurable<bool> clusterRejectExotics{"clusterRejectExotics", true, "Reject exotic clusters"};

  // Configurable<float> pTHatMaxMCD{"pTHatMaxMCD", 999.0, "maximum fraction of hard scattering for jet acceptance in detector MC"};
  // Configurable<float> pTHatMaxMCP{"pTHatMaxMCP", 999.0, "maximum fraction of hard scattering for jet acceptance in particle MC"};
  // Configurable<float> pTHatExponent{"pTHatExponent", 6.0, "exponent of the event weight for the calculation of pTHat"};

  int trackSelection = -1;
  int eventSelection = -1;
  std::string particleSelection;

  Service<o2::framework::O2DatabasePDG> pdgDatabase;

  void init(o2::framework::InitContext&)
  {
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));
    eventSelection = jetderiveddatautilities::initialiseEventSelection(static_cast<std::string>(eventSelections));
    particleSelection = static_cast<std::string>(particleSelections);

    // JetTrack QA histograms
    if (doprocessTracks) {
      registry.add("h_collisions", "event status; event status;entries", {HistType::kTH1F, {{4, 0., 4.0}}});
      registry.add("h_track_pt", "track pT;#it{p}_{T,track} (GeV/#it{c});entries", {HistType::kTH1F, {{350, 0., 350.}}});
      registry.add("h_track_eta", "track #eta;#eta_{track};entries", {HistType::kTH1F, {{100, -1., 1.}}});
      registry.add("h_track_phi", "track #varphi;#varphi_{track};entries", {HistType::kTH1F, {{160, 0., 7.}}});

      // Cluster QA histograms
      registry.add("h_cluster_pt", "cluster pT;#it{p}_{T_cluster} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}});
      registry.add("h_cluster_eta", "cluster #eta;#eta_{cluster};entries", {HistType::kTH1F, {{100, -1., 1.}}});
      registry.add("h_cluster_phi", "cluster #varphi;#varphi_{cluster};entries", {HistType::kTH1F, {{160, 0., 7.}}});
      registry.add("h_cluster_energy", "cluster #varphi;#varphi_{cluster};entries", {HistType::kTH1F, {{160, 0., 7.}}});
    }

    // Jet QA histograms
    if (doprocessJetsData) {
      registry.add("h_full_jet_pt", "jet pT;#it{p}_{T_jet} (GeV/#it{c});entries", {HistType::kTH1F, {{350, 0., 350.}}});
      registry.add("h_full_jet_eta", "jet #eta;#eta_{jet};entries", {HistType::kTH1F, {{100, -1., 1.}}});
      registry.add("h_full_jet_phi", "jet #varphi;#varphi_{jet};entries", {HistType::kTH1F, {{160, 0., 7.}}});

    }
    // if (doprocessJetsMCP) {
    //   registry.add("h_full_jet_pt_part", "jet pT;#it{p}_{T_jet} (GeV/#it{c});entries", {HistType::kTH1F, {{350, 0., 350.}}});
    //   registry.add("h_full_jet_eta_part", "jet #eta;#eta_{jet};entries", {HistType::kTH1F, {{100, -1., 1.}}});
    //   registry.add("h_full_jet_phi_part", "jet #varphi;#varphi_{jet};entries", {HistType::kTH1F, {{160, 0., 7.}}});
    // }
    // if (doprocessJetsMCPMCDMatched) {
    //   registry.add("h_full_jet_energyscaleDet", "Jet Energy Scale (det); p_{t,det} (GeV/c); (p_{t,det} - p_{t,part})/p_{t,part}", {HistType::kTH2F,{{400, 0., 400.}, {200, -1.,1.}}});
    //   // registry.add("h_full_jet_energyscaleDetCharged", "Jet Energy Scale (det, charged part); p_{t,det} (GeV/c); (p_{t,det} - p_{t,part})/p_{t,part}", {HistType::kTH2F,{{400, 0., 400., 200, -1.,1.}}});
    //   // registry.add("h_full_jet_energyscaleDetNeutral", "Jet Energy Scale (det, neutral part); p_{t,det} (GeV/c); (p_{t,det} - p_{t,part})/p_{t,part}", {HistType::kTH2F,{{400, 0., 400., 200, -1.,1.}}});
    //   // registry.add("h_full_jet_energyscaleDetChargedVsFull", "Jet Energy Scale (det, charged part, vs. full jet pt); p_{t,det} (GeV/c); (p_{t,det} - p_{t,part})/p_{t,part}", {HistType::kTH2F,{{400, 0., 400., 200, -1.,1.}}});
    //   // registry.add("h_full_jet_energyscaleDetNeutralVsFull", "Jet Energy Scale (det, neutral part, vs. full jet pt); p_{t,det} (GeV/c); (p_{t,det} - p_{t,part})/p_{t,part}", {HistType::kTH2F,{{400, 0., 400., 200, -1.,1.}}});
    //   registry.add("h_full_jet_energyscalePart", "Jet Energy Scale (part); p_{t,part} (GeV/c); (p_{t,det} - p_{t,part})/p_{t,part}", {HistType::kTH2F,{{400, 0., 400.}, {200, -1.,1.}}});
    //   // registry.add("h_full_jet_energyscaleCharged", "Jet Energy Scale (charged part); p_{t,part} (GeV/c); (p_{t,det} - p_{t,part})/p_{t,part}", {HistType::kTH2F,{{400, 0., 400., 200, -1.,1.}}});
    //   // registry.add("h_full_jet_energyscaleNeutral", "Jet Energy Scale (neutral part); p_{t,part} (GeV/c); (p_{t,det} - p_{t,part})/p_{t,part}", {HistType::kTH2F,{{400, 0., 400., 200, -1.,1.}}});
    //   // registry.add("h_full_jet_energyscaleChargedVsFull", "Jet Energy Scale (charged part, vs. full jet pt); p_{t,part} (GeV/c); (p_{t,det} - p_{t,part})/p_{t,part}", {HistType::kTH2F,{{400, 0., 400., 200, -1.,1.}}});
    //   // registry.add("h_full_jet_energyscaleNeutralVsFull", "Jet Energy Scale (neutral part, vs. full jet pt); p_{t,part} (GeV/c); (p_{t,det} - p_{t,part})/p_{t,part}", {HistType::kTH2F,{{400, 0., 400., 200, -1.,1.}}});
    // }

  } // init

  using FullJetTableDataJoined = soa::Join<aod::FullJets, aod::FullJetConstituents>;
  // using JetTableMCDJoined = soa::Join<aod::FullMCDetectorLevelJets, aod::FullMCDetectorLevelJetConstituents>;
  // using JetTableMCPJoined = soa::Join<aod::FullMCParticleLevelJets, aod::FullMCParticleLevelJetConstituents>;
  // // using JetTableMCDMatchedJoined = soa::Join<JetTableMCD, JetConstituentTableMCD, JetMatchingTableMCDMCP>;
  // using JetTableMCPMatchedJoined = soa::Join<JetTableMCP, JetConstituentTableMCP, JetMatchingTableMCPMCD>;

  // Applying some cuts(filters) on collisions, tracks, clusters
  // Filter eventCuts = (nabs(aod::jcollision::posZ) < VertexZCut);
  Filter eventCuts = (nabs(aod::jcollision::posZ) < VertexZCut && aod::jcollision::centrality >= centralityMin && aod::jcollision::centrality < centralityMax);
  Filter trackCuts = (aod::jtrack::pt >= trackpTMin && aod::jtrack::pt < trackpTMax && aod::jtrack::eta > trackEtaMin && aod::jtrack::eta < trackEtaMax && aod::jtrack::phi >= trackPhiMin && aod::jtrack::phi <= trackPhiMax);
  aod::EMCALClusterDefinition clusterDefinition = aod::emcalcluster::getClusterDefinitionFromString(clusterDefinitionS.value);
  Filter clusterFilter = (aod::jcluster::definition == static_cast<int>(clusterDefinition) && aod::jcluster::eta > clusterEtaMin && aod::jcluster::eta < clusterEtaMax && aod::jcluster::phi >= clusterPhiMin && aod::jcluster::phi <= clusterPhiMax && aod::jcluster::energy >= clusterEnergyMin && aod::jcluster::time > clusterTimeMin && aod::jcluster::time < clusterTimeMax && (clusterRejectExotics && aod::jcluster::isExotic != true));

  template <typename T, typename U>
  bool isAcceptedJet(U const& jet)
  {

    if (jetAreaFractionMin > -98.0) {
      if (jet.area() < jetAreaFractionMin * M_PI * (jet.r() / 100.0) * (jet.r() / 100.0)) {
        return false;
      }
    }
    if (leadingConstituentPtMin > -98.0) {
      bool isMinleadingConstituent = false;
      for (auto& constituent : jet.template tracks_as<T>()) {
        if (constituent.pt() >= leadingConstituentPtMin) {
          isMinleadingConstituent = true;
          break;
        }
      }
      // could add clusters later if needed
      if (!isMinleadingConstituent) {
        return false;
      }
    }
    return true;
  }

  template <typename T>
  void fillJetHistograms(T const& jet, float weight = 1.0)
  {
      if (jet.r() == round(JetRadii * 100.0f)) {
        registry.fill(HIST("h_full_jet_pt"), jet.pt(), weight);
        registry.fill(HIST("h_full_jet_eta"), jet.eta(), weight);
        registry.fill(HIST("h_full_jet_phi"), jet.phi(), weight);
        // registry.fill(HIST("h_full_jet_energyscaleDet"), jet.phi(), weight);
    }
  }

  // template <typename T>
  // void fillMCPHistograms(T const& jet, float weight = 1.0)
  // {
  //   if (jet.r() == round(JetRadii * 100.0f)) {
  //     registry.fill(HIST("h_full_jet_pt_part"), jet.pt(), weight);
  //     registry.fill(HIST("h_full_jet_eta_part"), jet.eta(), weight);
  //     registry.fill(HIST("h_full_jet_phi_part"), jet.phi(), weight);
  //     // registry.fill(HIST("h_full_jet_ntracks_part"), jet.tracksIds().size(), weight);
  //   }
  // }

  template <typename T, typename U>
  void fillTrackHistograms(T const& tracks, U const& clusters, float weight = 1.0)
  {
    for (auto const& track : tracks) {
      if (!jetderiveddatautilities::selectTrack(track, trackSelection)) {
        continue;
      }
      registry.fill(HIST("h_track_pt"), track.pt(), weight);
      registry.fill(HIST("h_track_eta"), track.eta(), weight);
      registry.fill(HIST("h_track_phi"), track.phi(), weight);
    }
    for (auto const& cluster : clusters) {
      double clusterpt = cluster.energy() / std::cosh(cluster.eta());
      registry.fill(HIST("h_cluster_pt"), clusterpt, weight);
      registry.fill(HIST("h_cluster_eta"), cluster.eta(), weight);
      registry.fill(HIST("h_cluster_phi"), cluster.phi(), weight);
      registry.fill(HIST("h_cluster_energy"), cluster.energy(), weight);
    }
  }

  // template <typename T, typename U>
  // void fillMatchedHistograms(T const& jetBase, float weight = 1.0)
  // {
  //
  //   float pTHat = 10. / (std::pow(weight, 1.0 / pTHatExponent));
  //   if (jetBase.pt() > pTHatMaxMCD * pTHat) { //Here, jetBase = mcd jets and jetTag = mcp jets
  //     return;
  //   }
  //
  //   if (jetBase.has_matchedJetGeo()) {  //geometrical jet matching only needed for pp
  //     for (auto& jetTag : jetBase.template matchedJetGeo_as<std::decay_t<U>>()) {
  //       if (jetTag.pt() > pTHatMaxMCP * pTHat) {
  //         continue;
  //       }
  //       registry.fill(HIST("h_full_jet_energyscaleDet"), jetBase.pt(), (jetBase.pt() - jetTag.pt())/ jetTag.pt() , weight);
  //       registry.fill(HIST("h_full_jet_energyscalePart"), jetTag.pt(), (jetBase.pt() - jetTag.pt())/ jetTag.pt(), weight);
  //
  //       // registry.fill(HIST("h3_jet_r_jet_eta_tag_jet_eta_base_matchedgeo"), jetBase.r() / 100.0, jetTag.eta(), jetBase.eta(), weight);
  //       // registry.fill(HIST("h3_jet_r_jet_phi_tag_jet_phi_base_matchedgeo"), jetBase.r() / 100.0, jetTag.phi(), jetBase.phi(), weight);
  //       // registry.fill(HIST("h3_jet_r_jet_ntracks_tag_jet_ntracks_base_matchedgeo"), jetBase.r() / 100.0, jetTag.tracksIds().size(), jetBase.tracksIds().size(), weight);
  //       // registry.fill(HIST("h3_jet_r_jet_pt_tag_jet_pt_base_diff_matchedgeo"), jetBase.r() / 100.0, jetTag.pt(), (jetTag.pt() - jetBase.pt()) / jetTag.pt(), weight);
  //       // registry.fill(HIST("h3_jet_r_jet_pt_tag_jet_eta_base_diff_matchedgeo"), jetBase.r() / 100.0, jetTag.pt(), (jetTag.eta() - jetBase.eta()) / jetTag.eta(), weight);
  //       // registry.fill(HIST("h3_jet_r_jet_pt_tag_jet_phi_base_diff_matchedgeo"), jetBase.r() / 100.0, jetTag.pt(), (jetTag.phi() - jetBase.phi()) / jetTag.phi(), weight);
  //
  //       // if (jetBase.r() == round(selectedJetsRadius * 100.0f)) {
  //       //   registry.fill(HIST("h3_jet_pt_tag_jet_eta_tag_jet_eta_matchedgeo"), jetTag.pt(), jetTag.eta(), jetBase.eta(), weight);
  //       //   registry.fill(HIST("h3_jet_pt_tag_jet_phi_tag_jet_phi_matchedgeo"), jetTag.pt(), jetTag.phi(), jetBase.phi(), weight);
  //       //   registry.fill(HIST("h3_jet_pt_tag_jet_ntracks_tag_jet_ntracks_matchedgeo"), jetTag.pt(), jetTag.tracksIds().size(), jetBase.tracksIds().size(), weight);
  //       // }
  //     }
  //   }
  // }

  void processDummy(JetCollisions const& collisions)
  {
  }
  PROCESS_SWITCH(FullJetSpectrapp, processDummy, "dummy task", true);

  void processJetsData(soa::Filtered<JetCollisions>::iterator const& collision, FullJetTableDataJoined const& jets, JetTracks const&, JetClusters const&)
  {
    for (auto const& jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (!isAcceptedJet<JetTracks>(jet)) {
        continue;
      }
      fillJetHistograms(jet);
    }
  }
  PROCESS_SWITCH(FullJetSpectrapp, processJetsData, "Full Jets Data", false);

  // void processJetsMCD(soa::Filtered<JetCollisions>::iterator const& collision, JetTableMCDJoined const& jets, JetTracks const& tracks, JetClusters const& clusters)
  // {
  //   for (auto const& jet : jets) {
  //     if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
  //       continue;
  //     }
  //     // fillJetHistograms(jet, collision.centrality());
  //   }
  // }
  // PROCESS_SWITCH(FullJetSpectrapp, processJetsMCD, "Full Jets MCD", false);
  //
  // void processJetsMCP(typename JetTableMCPJoined::iterator const& jet, JetParticles const& particles)
  // {
  //   if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
  //     return;
  //   }
  //   // fillMCPHistograms(jet);
  // }
  // PROCESS_SWITCH(FullJetSpectrapp, processJetsMCP, "Full Jets MCP", false);

  void processTracks(JetCollision const& collision, soa::Filtered<JetTracks> const& tracks, soa::Filtered<JetClusters> const& clusters)
  {
    registry.fill(HIST("h_collisions"), 0.5);
    if (!jetderiveddatautilities::eventEMCAL(collision)) {
      return;
    }
    registry.fill(HIST("h_collisions"), 1.5);
    fillTrackHistograms(tracks, clusters, 1.0);
  }
  PROCESS_SWITCH(FullJetSpectrapp, processTracks, "QA for fulljet tracks", false);

 //  void doprocessJetsMCPMCDMatched(JetCollision const&, JetTableMCDMatchedJoined const& mcdjets, JetTableMCPMatchedJoined const&, JetTracks const&, JetClusters const&, JetParticles const&)
 //  {
 //    for (const auto& mcdjet: mcdjets) {
 //      if (!jetfindingutilities::isInEtaAcceptance(mcdjet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
 //        continue;
 //    }
 //    fillMatchedHistograms<typename JetTableMCDMatchedJoined::iterator, JetTableMCPMatchedJoined>(mcdjet);
 //   }
 // }
 // PROCESS_SWITCH(FullJetSpectrapp, processJetsMCPMCDMatched, "full jet finder MCP matched to MCD", false);

}; // struct

// using FullJetsSpectraTask = FullJetSpectraTask<aod::FullJets, aod::FullJetConstituents, aod::FullMCDetectorLevelJets, aod::FullMCDetectorLevelJetConstituents, aod::FullMCDetectorLevelJetsMatchedToFullMCParticleLevelJets, aod::FullMCParticleLevelJets, aod::FullMCParticleLevelJetConstituents, aod::FullMCParticleLevelJetsMatchedToFullMCDetectorLevelJets>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
  adaptAnalysisTask<FullJetSpectrapp>(cfgc, TaskName{"full-jet-spectra-pp"})};
}
