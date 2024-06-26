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

// **Hadronic Correction in the EMCAL framework: to avoid the double counting of the charged particles' contribution in jets**
/// \author Archita Rani Dash <archita.rani.dash@cern.ch>

#include <algorithm>
#include <iostream>
#include <memory>
#include <unordered_map>
#include <cmath>
#include <TF1.h>
#include <vector>

#include "CCDB/BasicCCDBManager.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/HistogramRegistry.h"

#include "DetectorsBase/GeometryManager.h"

#include "PWGJE/DataModel/EMCALClusters.h"
#include "PWGJE/DataModel/EMCALMatchedTracks.h"
#include "PWGJE/DataModel/EMCALMatchedCollisions.h"
#include "PWGJE/DataModel/emcalCorrectionClusterHadronicCorrectionTask.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/PIDResponse.h"
#include "DataFormatsEMCAL/Cell.h"
#include "DataFormatsEMCAL/Constants.h"
#include "DataFormatsEMCAL/AnalysisCluster.h"
#include "EMCALBase/Geometry.h"
#include "EMCALBase/ClusterFactory.h"
#include "EMCALBase/NonlinearityHandler.h"
#include "EMCALReconstruction/Clusterizer.h"
#include "PWGJE/Core/JetUtilities.h"
#include "TVector2.h"

#include "CommonDataFormat/InteractionRecord.h"


using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using collisionEvSelIt = o2::soa::Join<o2::aod::Collisions, o2::aod::EvSels>::iterator;
using selectedClusters = o2::soa::Filtered<o2::aod::EMCALClusters>;
using myTracks  =  o2::soa::Filtered<o2::soa::Join<o2::aod::pidTPCFullEl, o2::aod::pidTPCFullPi, o2::aod::FullTracks, o2::aod::TrackSelection>>;


struct EmcalCorrectionClusterHadronicCorrectionTask {
  Produces<o2::aod::EmcalHCs> hadroniccorrectedclusters;

  HistogramRegistry registry;

  Preslice<o2::aod::EMCALMatchedTracks> perClusterMatchedTracks = o2::aod::emcalclustercell::emcalclusterId;  // looking at clusterID column in the EMCALMatchedTracks for every cluster

  //define configurables here
  Configurable<int> mClusterDefinition{"clusterDefinition", 10, "cluster definition to be selected, e.g. 10=kV3Default"};
  Configurable<float> minTime{"minTime", -25., "Minimum cluster time for time cut"};
  Configurable<float> maxTime{"maxTime", +20., "Maximum cluster time for time cut"};
  Configurable<float> minM02{"minM02", 0.1, "Minimum M02 for M02 cut"};
  Configurable<float> maxM02{"maxM02", 0.9, "Maximum M02 for M02 cut"};
  Configurable<float> minTrackPt{"minTrackPt", 0.15, "Minimum pT for tracks"};
  Configurable<double> fHadCorr1{"HadCorr1", 1., "hadronic correction fraction for complete cluster energy subtraction for one matched track" };    //100% - default
  Configurable<double> fHadCorr2{"HadCorr2", 0.7, "hadronic correction fraction for systematic studies for one matched track"};                     //70%
  Configurable<double> fHadCorralltrks1{"HadCorralltrks1", 1., "hadronic correction fraction for complete cluster energy subtraction for all matched tracks" };    //100% - all tracks
  Configurable<double> fHadCorralltrks2{"HadCorralltrks2", 0.7, "hadronic correction fraction for systematic studies for all matched tracks"};                     //70%
  Configurable<float> minDEta{"minDEta", 0.01, "Minimum dEta between track and cluster"};
  Configurable<float> minDPhi{"minDPhi", 0.01, "Minimum dPhi between track and cluster"};
  Configurable<double> fEexcl{"Eexcl", 0., "Cell energy that cannot be subtracted"};

  //pT-dependent track-matching configurables
  Configurable<float> Eta0{"eta0", 0.04, "Param 0 in eta for pt-dependent matching"};
  Configurable<float> Eta1{"eta1", 0.010, "Param 1 in eta for pt-dependent matching"};
  Configurable<float> Eta2{"eta2", 2.5, "Param 2 in eta for pt-dependent matching"};
  Configurable<float> Phi0{"phi0", 0.09, "Param 0 in phi for pt-dependent matching"};
  Configurable<float> Phi1{"phi1", 0.015, "Param 1 in phi for pt-dependent matching"};
  Configurable<float> Phi2{"phi2", 2.0, "Param 2 in phi for pt-dependent matching"};
  Configurable<double> PhiMatch{"phiMatch", 0.050, "phi match value in pp"};
  Configurable<double> EtaMatch{"etaMatch", 0.025, "eta match value in pp"};

  Configurable<bool> doHadCorrSyst{"doHadCorrSyst", false, "Do hadronic correction for systematic studies"};
  Configurable<bool> doMomDepMatching{"doMomDepMatching", true, "Do momentum dependent track matching"}; // to be always set to true in Run 3

  // Configurable<bool> doHadCorrOneTrack{"doHadCorrOneTrack", false, "Do hadronic correction with one track only"};  //for clusters with only one matched track
  // Configurable<bool> doHadCorrAllTracks{"doHadCorrAllTracks", true, "Do hadronic correction with all tracks"};    //for clusters with more than one matched tracks

  void init(o2::framework::InitContext&)  {

  //Event histograms
  registry.add("h_allcollisions", "Total events; event status;entries", {HistType::kTH1F, {{1, 0.5, 1.5}}});
  registry.add("h_acceptedcollisions", "Accepted events", {HistType::kTH1F, {{1, 0.5, 1.5}}});

  //Matched-Cluster histograms
  registry.add("h_matchedclusters", "Total matched clusters; cluster status;entries", {HistType::kTH1F, {{1, 0.5, 1.5}}});

  //Matched-Track histograms
  registry.add("h_matchedtracks", "Total matched tracks; track status;entries", {HistType::kTH1F, {{1, 0.5, 1.5}}});

  }

  Filter clusterDefinitionSelection = (o2::aod::emcalcluster::definition == mClusterDefinition) && (o2::aod::emcalcluster::time >= minTime) && (o2::aod::emcalcluster::time <= maxTime) && (o2::aod::emcalcluster::m02 > minM02) && (o2::aod::emcalcluster::m02 < maxM02);
  Filter trackSelection = (o2::aod::track::pt >= minTrackPt);
  //The matching of clusters and tracks is already centralised in the EMCAL framework.
  //One only needs to apply a filter on matched clusters
  //Here looping over all collisions matched to EMCAL clusters
  void processMatchedCollisions(collisionEvSelIt const& collision, selectedClusters const& clusters, o2::aod::EMCALMatchedTracks const& matchedtracks, myTracks const&)
  {
    registry.fill(HIST("h_allcollisions"), 1);

      // skip events with no clusters
      if (clusters.size() == 0) {
        return;
      }

      //Looping over all clusters matched to the collision
      for (const auto& cluster : clusters) {
      registry.fill(HIST("h_matchedclusters"), 1);
      // double pT, mom;
      double Ecluster1; double Ecluster2; double EclusterAll1; double EclusterAll2;
      Ecluster1 = Ecluster2 = EclusterAll1 = EclusterAll2 = cluster.energy();
      // double pT = cluster.energy() / cosh(cluster.eta());

      //selecting ALL MATCHED TRACKS after slicing all entries in perClusterMatchedTracks by the cluster globalIndex
      auto tracksofcluster = matchedtracks.sliceBy(perClusterMatchedTracks, cluster.globalIndex());
      int cmt = 0;             // counter for closest matched track
      double totalTrkP = 0.0;  // counter for total track momentum

      // pT-dependent track-matching instead of PID based track-matching to be adapted from Run 2 - suggested by Markus

      TF1 funcPtDepEta("func", "[1] + 1 / pow(x + pow(1 / ([0] - [1]), 1 / [2]), [2])");
      funcPtDepEta.SetParameters(Eta0, Eta1, Eta2);
      TF1 funcPtDepPhi("func", "[1] + 1 / pow(x + pow(1 / ([0] - [1]), 1 / [2]), [2])");
      funcPtDepEta.SetParameters(Phi0, Phi1, Phi2);

      //Looping over all matched tracks for the cluster
      //Total number of matched tracks = 20 (hard-coded)
      for (const auto& match : tracksofcluster) {
        // bool doHadCorrOneTrack = false;
        double etadiff = 999.;
        double phidiff = 999.;

        double mom = abs(match.track_as<myTracks>().p());
        registry.fill(HIST("h_matchedtracks"), 1);

        //CASE 1: no matched tracks or tracks with a very low pT -> skip clusters and don't subtract any cluster energy
        //Should I also write this non-subracted energy value to the table???
        if ((tracksofcluster.size() == 0) || (mom < 1e-6)) {
          // return Ecluster;
          continue;
        } // end CASE 1

        //CASE 2:
        // - If one matched track-> check if it's the closest one. If yes, subtract 100%.
        // - If more than one matched track -> do the above for the closest one (100%) and then loop over the rest for 100% energy subtraction
        // - If you want to do systematic studies -> perform the above two checks and then subtract 70% energy instead of 100%

        if (tracksofcluster.size() != 0) {
          auto trackEta = match.track_as<myTracks>().eta();
          auto trackPhi = match.track_as<myTracks>().phi();
          double dPhi = trackPhi - cluster.phi();
          double dEta = trackEta - cluster.eta();

          if (fabs(dEta) >= minDEta || fabs(dPhi) >= minDPhi) { // dEta and dPhi cut : ensures that the matched track is within the desired eta/phi window
            continue;
          }

          //Do pT-dependent track matching
          if(doMomDepMatching)
          {
            trackEta     = funcPtDepEta.Eval(mom);
            auto trackPhiHigh = +funcPtDepPhi.Eval(mom);
            auto trackPhiLow  = -funcPtDepPhi.Eval(mom);

          if ((phidiff < trackPhiHigh && phidiff > trackPhiLow) && fabs(etadiff) < trackEta){
            totalTrkP += mom;
          }

          if (((minDPhi < trackPhiHigh && minDPhi > trackPhiLow) && fabs(minDEta) < trackEta) && fHadCorr1 > 1 && cmt == 0)  {         // 100% energy subtraction for only the one closest matched track
            Ecluster1 -= fHadCorr1 * mom;
            //write this corrected energy to the table
            hadroniccorrectedclusters(cluster.energy(), 0.f, 0.f, 0.f);
          }
          if (Ecluster1 < 0) Ecluster1 = 0;

          if (((minDPhi < trackPhiHigh && minDPhi > trackPhiLow) && fabs(minDEta) < trackEta) && fHadCorralltrks1 > 0) {              // 100% energy subtraction for all tracks;
            auto Esub = fHadCorralltrks1 * totalTrkP;
            if (Esub > EclusterAll1) Esub = EclusterAll1;

            //applying Peter's algo from Run 2 : to not subtract the full energy of the cluster
            auto clusEexcl = fEexcl * cluster.nCells();

            if (EclusterAll1 < clusEexcl) clusEexcl = EclusterAll1;
            if ((EclusterAll1 - Esub) < clusEexcl) Esub = EclusterAll1 - clusEexcl;

            // To DO: Implement M02SubtractionScheme
            //
            //

            EclusterAll1 -= Esub;
            hadroniccorrectedclusters(0.f, 0.f, cluster.energy(), 0.f);
          }
          if (EclusterAll1 < 0) EclusterAll1 = 0;

          if (doHadCorrSyst)  {                    // if you want to subtract 70% energy (as was in Run 2) for systematic studies
            if (((minDPhi < trackPhiHigh && minDPhi > trackPhiLow) && fabs(minDEta) < trackEta) && fHadCorr2 > 1 && cmt == 0)  {     // 70% energy subtraction for only the one closest track
              Ecluster2 -= fHadCorr2 * mom;
              hadroniccorrectedclusters(0.f, cluster.energy(), 0.f, 0.f);
            }
            if (Ecluster2 < 0) Ecluster2 = 0;

            if (((minDPhi < trackPhiHigh && minDPhi > trackPhiLow) && fabs(minDEta) < trackEta) && fHadCorralltrks2 > 0) {         // 70% energy subtraction for all tracks
              EclusterAll2 -= fHadCorralltrks2 * mom;
              hadroniccorrectedclusters(0.f, 0.f, 0.f, cluster.energy());
            }
            if (EclusterAll2 < 0) EclusterAll2 = 0;
          }
        } //doMomDepMatching ends
        } // end of CASE 2
      } //end of track loop
    } //end of cluster loop
  } //end of process function
  PROCESS_SWITCH(EmcalCorrectionClusterHadronicCorrectionTask, processMatchedCollisions, "Process matched clusters from collision", true);
}; //end of struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<EmcalCorrectionClusterHadronicCorrectionTask>(cfgc, TaskName{"emcal-correction-cluster-hadronic-correction-task"})}; }
