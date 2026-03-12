
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
#include "PWGDQ/DataModel/ReducedInfoTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/AggregatedRunInfo.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

namespace o2::aod {
DECLARE_SOA_TABLE(RctRawDQ, "AOD", "RCTRAWDQ", evsel::Rct);
};

using namespace o2;

using namespace o2::framework;
using namespace o2::aod::evsel;

static const int32_t nBCsPerOrbit = o2::constants::lhc::LHCMaxBunches;

// Converts Stra Event selections from 004 to 005
struct rctRawProducer {
  Produces<aod::RctRawDQ> straEvSels_005;

  Service<o2::ccdb::BasicCCDBManager> ccdb;

  int lastRun = -1;
  int64_t lastTF = -1;
  uint32_t lastRCT = 0;
  uint64_t sorTimestamp = 0; // default SOR timestamp
  uint64_t eorTimestamp = 1; // default EOR timestamp
  int64_t bcSOR = -1;        // global bc of the start of run
  int64_t nBCsPerTF = -1;    // duration of TF in bcs, should be 128*3564 or 3
  std::map<uint64_t, uint32_t>* mapRCT = nullptr;

  uint32_t getRctRaw(int run, uint64_t timestamp, uint64_t globalBC)
  {
    if (run != lastRun) {
      lastRun = run;
      auto runInfo = o2::parameters::AggregatedRunInfo::buildAggregatedRunInfo(o2::ccdb::BasicCCDBManager::instance(), run);
      // first bc of the first orbit
      bcSOR = runInfo.orbitSOR * nBCsPerOrbit;
      // duration of TF in bcs
      nBCsPerTF = runInfo.orbitsPerTF * nBCsPerOrbit;
      // SOR and EOR timestamps
      sorTimestamp = runInfo.sor;
      eorTimestamp = runInfo.eor;
      // timestamp of the middle of the run used to access run-wise CCDB entries
      int64_t ts = runInfo.sor / 2 + runInfo.eor / 2;
      // QC info
      std::map<std::string, std::string> metadata;
      metadata["run"] = Form("%d", run);
      ccdb->setFatalWhenNull(0);
      mapRCT = ccdb->getSpecific<std::map<uint64_t, uint32_t>>("Users/j/jian/RCT", ts, metadata);
      ccdb->setFatalWhenNull(1);
      if (mapRCT == nullptr) {
        LOGP(info, "rct object missing... inserting dummy rct flags");
        mapRCT = new std::map<uint64_t, uint32_t>;
        mapRCT->insert(std::pair<uint64_t, uint32_t>(sorTimestamp, 0));
      }
    }

    // store rct flags
    uint32_t rct = lastRCT;
    int64_t thisTF = (globalBC - bcSOR) / nBCsPerTF;
    if (mapRCT != nullptr && thisTF != lastTF) { // skip for unanchored runs; do it once per TF
      auto itrct = mapRCT->upper_bound(timestamp);
      if (itrct != mapRCT->begin())
        itrct--;
      rct = itrct->second;
      LOGP(debug, "sor={} eor={} ts={} rct={}", sorTimestamp, eorTimestamp, timestamp, rct);
      lastRCT = rct;
      lastTF = thisTF;
    }
    return rct;
  }

  // void process(soa::Join<aod::StraEvSels_004, aod::StraStamps> const& straEvSels_004)
  void process(soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended> const& straEvSels_004)
  {
    for (auto& values : straEvSels_004) {
      straEvSels_005(getRctRaw(values.runNumber(), values.timestamp(), values.globalBC()) /* Rct value*/);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<straevselsconverter5>(cfgc)};
}
