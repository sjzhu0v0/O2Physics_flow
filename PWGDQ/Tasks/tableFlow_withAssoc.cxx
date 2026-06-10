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
//
// Contact: iarsene@cern.ch, i.c.arsene@fys.uio.no
//
#include "PWGDQ/Core/AnalysisCompositeCut.h"
#include "PWGDQ/Core/AnalysisCut.h"
#include "PWGDQ/Core/CutsLibrary.h"
#include "PWGDQ/Core/HistogramManager.h"
#include "PWGDQ/Core/HistogramsLibrary.h"
#include "PWGDQ/Core/MixingHandler.h"
#include "PWGDQ/Core/MixingLibrary.h"
#include "PWGDQ/Core/VarManager.h"
#include "PWGDQ/DataModel/ReducedInfoTables.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/RCTSelectionFlags.h"
#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/Propagator.h"
#include "Field/MagneticField.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "ITSMFTBase/DPLAlpideParam.h"

#include "TGeoGlobalMagField.h"
#include <TH1F.h>
#include <TH3F.h>
#include <THashList.h>
#include <TList.h>
#include <TString.h>

#include <algorithm>
#include <array>
#include <bitset>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

using std::cout;
using std::endl;
using std::string;

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod;
using namespace o2::aod::rctsel;

// Some definitions
namespace o2::aod
{

namespace dqanalysisflags
{
// TODO: the barrel amd muon selection columns are bit maps so unsigned types should be used, however, for now this is not supported in Filter expressions
// TODO: For now in the tasks we just statically convert from unsigned int to int, which should be fine as long as we do
//      not use a large number of bits (>=30)
// Bcandidate columns for ML analysis of B->Jpsi+K
DECLARE_SOA_COLUMN(MixingHash, mixingHash, int);
DECLARE_SOA_BITMAP_COLUMN(IsEventSelected, isEventSelected, 8);
DECLARE_SOA_BITMAP_COLUMN(IsBarrelSelected, isBarrelSelected, 32);
DECLARE_SOA_BITMAP_COLUMN(IsMuonSelected, isMuonSelected, 32);
DECLARE_SOA_BITMAP_COLUMN(IsBarrelSelectedPrefilter, isBarrelSelectedPrefilter, 32);
DECLARE_SOA_COLUMN(IsPrefilterVetoed, isPrefilterVetoed, int);
DECLARE_SOA_COLUMN(CacheTrackDaughter, cacheTrackDaughter, std::vector<int>);
DECLARE_SOA_COLUMN(HadronicRate, hadronicRate, float);
} // namespace dqanalysisflags

DECLARE_SOA_TABLE(EventCuts, "AOD", "DQANAEVCUTSA", dqanalysisflags::IsEventSelected);
DECLARE_SOA_TABLE(MixingHashes, "AOD", "DQANAMIXHASHA", dqanalysisflags::MixingHash);
DECLARE_SOA_TABLE(BarrelTrackCuts, "AOD", "DQANATRKCUTSA", dqanalysisflags::IsBarrelSelected);
DECLARE_SOA_TABLE(MuonTrackCuts, "AOD", "DQANAMUONCUTSA", dqanalysisflags::IsMuonSelected);
DECLARE_SOA_TABLE(Prefilter, "AOD", "DQPREFILTERA", dqanalysisflags::IsBarrelSelectedPrefilter);
DECLARE_SOA_TABLE(TrackInfo, "AOD", "TRACKINFO", collision::PosX, collision::PosY, collision::PosZ, collision::NumContrib, reducedtrack::Pt, reducedtrack::Eta, reducedtrack::Phi, reducedtrack::Sign, reducedtrack::DcaXY, reducedtrack::DcaZ, track::ITSClusterMap);
DECLARE_SOA_TABLE(HadronicRates, "AOD", "DQHADRATE", dqanalysisflags::HadronicRate);

namespace flowVec
{
DECLARE_SOA_COLUMN(MultFT0C, multFT0C, float);
DECLARE_SOA_COLUMN(MultVtxContri, multVtxContri, float);
DECLARE_SOA_COLUMN(VtxZ, vtxZ, float);
DECLARE_SOA_COLUMN(PT, pT, std::vector<float>);
DECLARE_SOA_COLUMN(Eta, eta, std::vector<float>);
DECLARE_SOA_COLUMN(Phi, phi, std::vector<float>);
DECLARE_SOA_COLUMN(Mass, mass, std::vector<float>);
DECLARE_SOA_COLUMN(Sign, sign, std::vector<float>);
DECLARE_SOA_COLUMN(QvecRe, qvecRe, std::vector<float>);
DECLARE_SOA_COLUMN(QvecIm, qvecIm, std::vector<float>);
DECLARE_SOA_COLUMN(QvecAmp, qvecAmp, std::vector<float>);
DECLARE_SOA_COLUMN(PTREF, pTRef, std::vector<float>);
DECLARE_SOA_COLUMN(EtaREF, etaRef, std::vector<float>);
DECLARE_SOA_COLUMN(PhiREF, phiRef, std::vector<float>);
DECLARE_SOA_COLUMN(ITSChi2NCl, iTSChi2NCl, std::vector<float>);
DECLARE_SOA_COLUMN(TPCNClsCR, tPCNClsCR, std::vector<float>);
DECLARE_SOA_COLUMN(TPCNClsFound, tPCNClsFound, std::vector<float>);
DECLARE_SOA_COLUMN(TPCChi2NCl, tPCChi2NCl, std::vector<float>);
DECLARE_SOA_COLUMN(TPCSignal, tPCSignal, std::vector<float>);
DECLARE_SOA_COLUMN(TPCNSigmaEl, tPCNSigmaEl, std::vector<float>);
DECLARE_SOA_COLUMN(TPCNSigmaPi, tPCNSigmaPi, std::vector<float>);
DECLARE_SOA_COLUMN(TPCNSigmaPr, tPCNSigmaPr, std::vector<float>);
DECLARE_SOA_COLUMN(DcaXY, dcaXY, std::vector<float>);
DECLARE_SOA_COLUMN(DcaZ, dcaZ, std::vector<float>);
DECLARE_SOA_COLUMN(Pt1, pt1, std::vector<float>);
DECLARE_SOA_COLUMN(Eta1, eta1, std::vector<float>);
DECLARE_SOA_COLUMN(Phi1, phi1, std::vector<float>);
DECLARE_SOA_COLUMN(Sign1, sign1, std::vector<int>);
DECLARE_SOA_COLUMN(ITSChi2NCl1, iTSChi2NCl1, std::vector<float>);
DECLARE_SOA_COLUMN(TPCNClsCR1, tPCNClsCR1, std::vector<float>);
DECLARE_SOA_COLUMN(TPCNClsFound1, tPCNClsFound1, std::vector<float>);
DECLARE_SOA_COLUMN(TPCChi2NCl1, tPCChi2NCl1, std::vector<float>);
DECLARE_SOA_COLUMN(TPCSignal1, tPCSignal1, std::vector<float>);
DECLARE_SOA_COLUMN(TPCNSigmaEl1, tPCNSigmaEl1, std::vector<float>);
DECLARE_SOA_COLUMN(TPCNSigmaPi1, tPCNSigmaPi1, std::vector<float>);
DECLARE_SOA_COLUMN(TPCNSigmaPr1, tPCNSigmaPr1, std::vector<float>);
DECLARE_SOA_COLUMN(Pt2, pt2, std::vector<float>);
DECLARE_SOA_COLUMN(Eta2, eta2, std::vector<float>);
DECLARE_SOA_COLUMN(Phi2, phi2, std::vector<float>);
DECLARE_SOA_COLUMN(Sign2, sign2, std::vector<int>);
DECLARE_SOA_COLUMN(ITSChi2NCl2, iTSChi2NCl2, std::vector<float>);
DECLARE_SOA_COLUMN(TPCNClsCR2, tPCNClsCR2, std::vector<float>);
DECLARE_SOA_COLUMN(TPCNClsFound2, tPCNClsFound2, std::vector<float>);
DECLARE_SOA_COLUMN(TPCChi2NCl2, tPCChi2NCl2, std::vector<float>);
DECLARE_SOA_COLUMN(TPCSignal2, tPCSignal2, std::vector<float>);
DECLARE_SOA_COLUMN(TPCNSigmaEl2, tPCNSigmaEl2, std::vector<float>);
DECLARE_SOA_COLUMN(TPCNSigmaPi2, tPCNSigmaPi2, std::vector<float>);
DECLARE_SOA_COLUMN(TPCNSigmaPr2, tPCNSigmaPr2, std::vector<float>);
DECLARE_SOA_COLUMN(ITSClusterMap, itsClusterMap, std::vector<uint8_t>);
} // namespace flowVec

namespace flowPair
{
DECLARE_SOA_COLUMN(PT, pT, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(Mass, mass, float);
DECLARE_SOA_COLUMN(Sign, sign, float);
DECLARE_SOA_COLUMN(PT1, pT1, float);
DECLARE_SOA_COLUMN(Eta1, eta1, float);
DECLARE_SOA_COLUMN(Phi1, phi1, float);
DECLARE_SOA_COLUMN(PT2, pT2, float);
DECLARE_SOA_COLUMN(Eta2, eta2, float);
DECLARE_SOA_COLUMN(Phi2, phi2, float);
} // namespace flowPair

DECLARE_SOA_TABLE(FlowVecs, "AOD", "DQFLOWVECSA", evsel::Selection, flowVec::MultFT0C, flowVec::MultVtxContri, flowVec::VtxZ, flowVec::PT, flowVec::Eta, flowVec::Phi, flowVec::Mass, flowVec::Sign, flowVec::QvecRe, flowVec::QvecIm, flowVec::QvecAmp);

DECLARE_SOA_TABLE(FlowVecD, "AOD", "DQFLOWVECDA", mult::MultTPC, mult::MultTracklets, mult::MultNTracksPV, mult::MultFT0C, collision::NumContrib, collision::PosX, collision::PosY, collision::PosZ, evsel::Selection, dqanalysisflags::HadronicRate, flowVec::PT, flowVec::Eta, flowVec::Phi, flowVec::Mass, flowVec::Sign, flowVec::PTREF, flowVec::EtaREF, flowVec::PhiREF, flowVec::ITSChi2NCl, flowVec::TPCNClsCR, flowVec::TPCNClsFound, flowVec::TPCChi2NCl, flowVec::ITSClusterMap, flowVec::TPCSignal, flowVec::TPCNSigmaEl, flowVec::TPCNSigmaPi, flowVec::TPCNSigmaPr, flowVec::DcaXY, flowVec::DcaZ, flowVec::Pt1, flowVec::Eta1, flowVec::Phi1, flowVec::Sign1, flowVec::ITSChi2NCl1, flowVec::TPCNClsCR1, flowVec::TPCNClsFound1, flowVec::TPCChi2NCl1, flowVec::TPCSignal1, flowVec::TPCNSigmaEl1, flowVec::TPCNSigmaPi1, flowVec::TPCNSigmaPr1, flowVec::Pt2, flowVec::Eta2, flowVec::Phi2, flowVec::Sign2, flowVec::ITSChi2NCl2, flowVec::TPCNClsCR2, flowVec::TPCNClsFound2, flowVec::TPCChi2NCl2, flowVec::TPCSignal2, flowVec::TPCNSigmaEl2, flowVec::TPCNSigmaPi2, flowVec::TPCNSigmaPr2);

DECLARE_SOA_TABLE(FlowPairRR, "AOD", "DQFLOWPAIRRRA", mult::MultTPC, mult::MultTracklets, mult::MultNTracksPV, evsel::Selection, flowVec::MultFT0C, flowVec::MultVtxContri, flowVec::VtxZ, flowPair::PT1, flowPair::Eta1, flowPair::Phi1, flowPair::PT2, flowPair::Eta2, flowPair::Phi2);
DECLARE_SOA_TABLE(FlowPairPR, "AOD", "DQFLOWPAIRPRA", mult::MultTPC, mult::MultTracklets, mult::MultNTracksPV, evsel::Selection, flowVec::MultFT0C, flowVec::MultVtxContri, flowVec::VtxZ, flowPair::PT, flowPair::Eta, flowPair::Phi, flowPair::Mass, flowPair::Sign, flowPair::PT1, flowPair::Eta1, flowPair::Phi1);
DECLARE_SOA_TABLE(RctRawDQ, "AOD", "RCTRAWDQ", evsel::Rct);
} // namespace o2::aod

// Declarations of various short names
using MyEvents = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended>;
using MyEventsSelected = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::EventCuts>;
using MyEventsHashSelected = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::EventCuts, aod::MixingHashes>;
using MyEventsVtxCov = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsVtxCov>;
using MyEventsVtxCovSelected = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsVtxCov, aod::EventCuts, aod::RctRawDQ>;
using MyEventsHashVtxCovSelected = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsVtxCov, aod::EventCuts, aod::MixingHashes>;
using MyEventsVtxCovSelectedQvector = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsVtxCov, aod::EventCuts, aod::ReducedEventsQvector>;
using MyEventsVtxCovSelectedQvectorExtraWithRefFlow = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsVtxCov, aod::EventCuts, aod::ReducedEventsQvector, aod::ReducedEventsQvectorExtra, aod::ReducedEventsRefFlow>;
using MyEventsVtxCovSelectedQvectorCentr = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsVtxCov, aod::EventCuts, aod::ReducedEventsQvectorCentr>;
using MyEventsQvector = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsQvector>;
using MyEventsQvectorMultExtra = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsQvector, aod::ReducedEventsMultPV, aod::ReducedEventsMultAll>;
using MyEventsQvectorCentr = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsQvectorCentr>;
using MyEventsQvectorCentrMultExtra = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsQvectorCentr, aod::ReducedEventsMultPV, aod::ReducedEventsMultAll>;
using MyEventsQvectorExtra = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsQvector, aod::ReducedEventsQvectorExtra>;
using MyEventsHashSelectedQvector = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::EventCuts, aod::MixingHashes, aod::ReducedEventsQvector>;
using MyEventsHashSelectedQvectorExtra = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::EventCuts, aod::MixingHashes, aod::ReducedEventsQvector, aod::ReducedEventsQvectorExtra>;
using MyEventsHashSelectedQvectorCentr = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::EventCuts, aod::MixingHashes, aod::ReducedEventsQvectorCentr>;

using MyBarrelTracks = soa::Join<aod::ReducedTracks, aod::ReducedTracksBarrel, aod::ReducedTracksBarrelPID>;
using MyBarrelTracksWithCov = soa::Join<aod::ReducedTracks, aod::ReducedTracksBarrel, aod::ReducedTracksBarrelCov, aod::ReducedTracksBarrelPID>;
using MyBarrelTracksSelected = soa::Join<aod::ReducedTracks, aod::ReducedTracksBarrel, aod::ReducedTracksBarrelPID, aod::BarrelTrackCuts>;
using MyBarrelTracksSelectedWithPrefilter = soa::Join<aod::ReducedTracks, aod::ReducedTracksBarrel, aod::ReducedTracksBarrelPID, aod::BarrelTrackCuts, aod::Prefilter>;
using MyBarrelTracksSelectedWithCov = soa::Join<aod::ReducedTracks, aod::ReducedTracksBarrel, aod::ReducedTracksBarrelCov, aod::ReducedTracksBarrelPID, aod::BarrelTrackCuts>;
using MyBarrelTracksSelectedWithColl = soa::Join<aod::ReducedTracks, aod::ReducedTracksBarrel, aod::ReducedTracksBarrelPID, aod::BarrelTrackCuts, aod::ReducedTracksBarrelInfo>;
using MyBarrelTrackAssocsSelected = soa::Join<aod::ReducedTracksAssoc, aod::BarrelTrackCuts>;
using MyDielectronCandidates = soa::Join<aod::Dielectrons, aod::DielectronsExtra>;
using MyMuonTracks = soa::Join<aod::ReducedMuons, aod::ReducedMuonsExtra>;
using MyMuonTracksSelected = soa::Join<aod::ReducedMuons, aod::ReducedMuonsExtra, aod::MuonTrackCuts>;
using MyMuonTracksWithCov = soa::Join<aod::ReducedMuons, aod::ReducedMuonsExtra, aod::ReducedMuonsCov>;
using MyMuonTracksSelectedWithCov = soa::Join<aod::ReducedMuons, aod::ReducedMuonsExtra, aod::ReducedMuonsCov, aod::MuonTrackCuts>;
using MyMuonTracksSelectedWithColl = soa::Join<aod::ReducedMuons, aod::ReducedMuonsExtra, aod::ReducedMuonsInfo, aod::MuonTrackCuts>;
using MyMftTracks = soa::Join<aod::ReducedMFTs, aod::ReducedMFTsExtra>;

// bit maps used for the Fill functions of the VarManager
constexpr static uint32_t gkEventFillMap = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventExtended;
constexpr static uint32_t gkEventFillMapWithCov = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventExtended | VarManager::ObjTypes::ReducedEventVtxCov;
constexpr static uint32_t gkEventFillMapWithQvector = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventExtended | VarManager::ObjTypes::ReducedEventQvector;
constexpr static uint32_t gkEventFillMapWithQvectorMultExtra = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventExtended | VarManager::ObjTypes::ReducedEventQvector | VarManager::ObjTypes::ReducedEventMultExtra;
constexpr static uint32_t gkEventFillMapWithQvectorCentr = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventExtended | VarManager::ObjTypes::CollisionQvect;
constexpr static uint32_t gkEventFillMapWithQvectorCentrMultExtra = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventExtended | VarManager::ObjTypes::CollisionQvect | VarManager::ObjTypes::ReducedEventMultExtra;
constexpr static uint32_t gkEventFillMapWithCovQvector = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventExtended | VarManager::ObjTypes::ReducedEventVtxCov | VarManager::ObjTypes::ReducedEventQvector;
constexpr static uint32_t gkEventFillMapWithCovQvectorExtraWithRefFlow = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventExtended | VarManager::ObjTypes::ReducedEventVtxCov | VarManager::ObjTypes::ReducedEventQvector | VarManager::ObjTypes::ReducedEventQvectorExtra | VarManager::ObjTypes::ReducedEventRefFlow;
constexpr static uint32_t gkEventFillMapWithCovQvectorCentr = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventExtended | VarManager::ObjTypes::ReducedEventVtxCov | VarManager::ObjTypes::CollisionQvect;
constexpr static uint32_t gkTrackFillMap = VarManager::ObjTypes::ReducedTrack | VarManager::ObjTypes::ReducedTrackBarrel | VarManager::ObjTypes::ReducedTrackBarrelPID;
constexpr static uint32_t gkTrackFillMapWithCov = VarManager::ObjTypes::ReducedTrack | VarManager::ObjTypes::ReducedTrackBarrel | VarManager::ObjTypes::ReducedTrackBarrelCov | VarManager::ObjTypes::ReducedTrackBarrelPID;
constexpr static uint32_t gkTrackFillMapWithColl = VarManager::ObjTypes::ReducedTrack | VarManager::ObjTypes::ReducedTrackBarrel | VarManager::ObjTypes::ReducedTrackBarrelPID | VarManager::ObjTypes::ReducedTrackCollInfo;

constexpr static uint32_t gkMuonFillMap = VarManager::ObjTypes::ReducedMuon | VarManager::ObjTypes::ReducedMuonExtra;
constexpr static uint32_t gkMuonFillMapWithCov = VarManager::ObjTypes::ReducedMuon | VarManager::ObjTypes::ReducedMuonExtra | VarManager::ObjTypes::ReducedMuonCov;
constexpr static uint32_t gkMuonFillMapWithColl = VarManager::ObjTypes::ReducedMuon | VarManager::ObjTypes::ReducedMuonExtra | VarManager::ObjTypes::ReducedMuonCollInfo;

constexpr static int pairTypeEE = VarManager::kDecayToEE;
constexpr static int pairTypeMuMu = VarManager::kDecayToMuMu;
constexpr static int pairTypeEMu = VarManager::kElectronMuon;

std::vector<int> findCommonBits(int num1, int num2);
std::vector<TString> splitTString(TString str, char c);
// Global function used to define needed histogram classes
void DefineHistograms(HistogramManager* histMan, TString histClasses, Configurable<std::string> configVar); // defines histograms for all tasks

struct AnalysisFlow {
  OutputObj<THashList> fOutputList{"output"};
  Configurable<float> fConfigDileptonLowMass{"cfgDileptonLowMass", 0.0, "Low mass cut for the dileptons used in analysis"};
  Configurable<float> fConfigDileptonHighMass{"cfgDileptonHighMass", 5.0, "High mass cut for the dileptons used in analysis"};
  Configurable<std::string> fConfigQVectorsBarrel2Cal{"cfgQVectorsBarrel", "4,1;2,1:0,1", "Semicolon separated list of q vector to be calculated"};
  Configurable<int> fConfigBarrelTrackCutBitMask{"cfgBarrelTrackCutBitMask", 1, "bit mask of the track cuts"};
  Configurable<int> fConfigMixingDepthRR{"cfgMixingDepthRR", 5, "Event mixing pool depth rr"};
  Configurable<int> fConfigMixingDepthPR{"cfgMixingDepthPR", 5, "Event mixing pool depth pr"};
  Configurable<std::string> fConfigRCTLabel{"cfgRCTLabel", "CBT", "RCT flag labels: CBT, CBT_hadronPID, CBT_electronPID, CBT_calo, CBT_muon, CBT_muon_glo"};

  Produces<aod::FlowVecs> qVectors;
  Produces<aod::FlowVecD> flowVectorsDetailed;
  Produces<aod::FlowPairRR> flowPairsRR;
  Produces<aod::FlowPairPR> flowPairsPR;

  Filter eventFilter = (aod::dqanalysisflags::isEventSelected & static_cast<uint8_t>(1)) > static_cast<uint8_t>(0);
  Filter dileptonFilter = aod::reducedpair::mass > fConfigDileptonLowMass&& aod::reducedpair::mass < fConfigDileptonHighMass;
  Filter filterBarrelTrackSelected = (aod::dqanalysisflags::isBarrelSelected & static_cast<uint32_t>(0xffffffff)) > static_cast<uint32_t>(0);

  std::vector<float> qvecRe;
  std::vector<float> qvecIm;
  std::vector<float> qvecAmp;
  std::vector<std::vector<std::array<int, 2>>> cfgsQVecBarrel2Cal;
  std::vector<int> cfgsIndexCutBitMask;
  int nCfgQvectorBarrel2Cal = 0;

  NoBinningPolicy<aod::dqanalysisflags::MixingHash> hashBin;
  HistogramManager* fHistMan = nullptr;
  bool isWithDilepton = false;
  RCTFlagsChecker rctChecker{"CBT"};

  Preslice<aod::ReducedTracksAssoc> perEventTrackAssocs = aod::reducedtrack_association::reducedeventId;
  Preslice<soa::Filtered<MyDielectronCandidates>> perEventPairs = aod::reducedpair::reducedeventId;

  void init(o2::framework::InitContext& context)
  {
    rctChecker.init(fConfigRCTLabel.value);
    VarManager::SetDefaultVarNames();
    fHistMan = new HistogramManager("analysisHistos", "aa", VarManager::kNVars);
    fHistMan->SetUseDefaultVariableNames(kTRUE);
    fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);

    if (context.mOptions.get<bool>("processSkimmed")) {
      isWithDilepton = true;
    }

    VarManager::SetUseVars(fHistMan->GetUsedVars());
    fOutputList.setObject(fHistMan->GetMainHistogramList());

    TString strConfigurationQVector2Cal = fConfigQVectorsBarrel2Cal.value;
    std::vector<TString> vecStrConfigurationQVector2Cal = splitTString(strConfigurationQVector2Cal, ':');

    for (auto& strCfg : vecStrConfigurationQVector2Cal) {
      std::vector<TString> vecStrCfg = splitTString(strCfg, ';');
      std::vector<std::array<int, 2>> vecCfg;
      for (auto& str : vecStrCfg) {
        std::vector<TString> vecStr = splitTString(str, ',');
        if (vecStr.size() != 2) {
          LOG(fatal) << "Invalid configuration for q vector: " << strCfg;
          exit(1);
        }
        vecCfg.push_back({vecStr[0].Atoi(), vecStr[1].Atoi()});
        nCfgQvectorBarrel2Cal++;
      }
      cfgsQVecBarrel2Cal.push_back(vecCfg);
    }

    LOG(info) << "Number of q vectors to calculate: " << nCfgQvectorBarrel2Cal;

    int numTrackCuts = 0;
    for (int i = 0; i < 32; i++) {
      if (fConfigBarrelTrackCutBitMask.value & (1 << i)) {
        numTrackCuts++;
        cfgsIndexCutBitMask.push_back(i);
      }
    }

    if (numTrackCuts != static_cast<int>(cfgsQVecBarrel2Cal.size())) {
      LOG(fatal) << "Number of track cuts and number of q vector sets to calculate do not match";
      exit(1);
    }
  }

  template <typename TEvent, typename TAssocs, typename TTracks, typename TDileptons>
  void fillQVectors(TEvent const& event, TAssocs const& assocs, TTracks const& tracks, TDileptons const& dileptons, bool withDilepton)
  {
    std::vector<int> daugDileptonGlobalIndexes;

    qvecRe.assign(nCfgQvectorBarrel2Cal, 0.0f);
    qvecIm.assign(nCfgQvectorBarrel2Cal, 0.0f);
    qvecAmp.assign(nCfgQvectorBarrel2Cal, 0.0f);

    if (withDilepton) {
      for (auto dilepton : dileptons) {
        daugDileptonGlobalIndexes.push_back(dilepton.index0Id());
        daugDileptonGlobalIndexes.push_back(dilepton.index1Id());
      }
    }

    for (auto& assoc : assocs) {
      if (!(assoc.isBarrelSelected_raw() & fConfigBarrelTrackCutBitMask.value)) {
        continue;
      }
      auto hadron = assoc.template reducedtrack_as<TTracks>();
      bool isDileptonDaug = false;
      for (int daugDileptonGlobalIndex : daugDileptonGlobalIndexes) {
        if (hadron.globalIndex() == daugDileptonGlobalIndex) {
          isDileptonDaug = true;
          break;
        }
      }
      if (isDileptonDaug) {
        continue;
      }

      float weff = 1.0, wacc = 1.0;
      std::vector<int> vecBits = findCommonBits(assoc.isBarrelSelected_raw(), fConfigBarrelTrackCutBitMask.value);
      for (auto bit : vecBits) {
        int index = 0;
        int indexBit = 0;
        for (int i = 0; i < static_cast<int>(cfgsIndexCutBitMask.size()); ++i) {
          if (cfgsIndexCutBitMask[i] == bit) {
            indexBit = i;
            break;
          }
        }
        for (int i = 0; i < indexBit; ++i) {
          index += cfgsQVecBarrel2Cal[i].size();
        }

        for (int iQvector = 0; iQvector < static_cast<int>(cfgsQVecBarrel2Cal[indexBit].size()); ++iQvector) {
          int npowPhi = cfgsQVecBarrel2Cal[indexBit][iQvector][0];
          int npowWeight = cfgsQVecBarrel2Cal[indexBit][iQvector][1];
          double weights = 1;
          for (int i = 0; i < npowWeight; ++i) {
            weights *= 1. / weff / wacc;
          }
          double phi = hadron.phi();
          qvecRe[index + iQvector] += weights * std::cos(npowPhi * phi);
          qvecIm[index + iQvector] += weights * std::sin(npowPhi * phi);
          qvecAmp[index + iQvector] += weights * weights;
        }
      }
    }

    std::vector<float> vecPT;
    std::vector<float> vecEta;
    std::vector<float> vecPhi;
    std::vector<float> vecMass;
    std::vector<float> vecSign;
    if (withDilepton) {
      for (auto dilepton : dileptons) {
        vecPT.push_back(dilepton.pt());
        vecEta.push_back(dilepton.eta());
        vecPhi.push_back(dilepton.phi());
        vecMass.push_back(dilepton.mass());
        vecSign.push_back(dilepton.sign());
      }
    }

    qVectors(event.selection_raw(), event.multFT0C(), event.numContrib(), event.posZ(), vecPT, vecEta, vecPhi, vecMass, vecSign, qvecRe, qvecIm, qvecAmp);
  }

  template <typename TEvent, typename TAssocs, typename TTracks, typename TDileptons>
  void fillSEFlowPairPoiRef(TEvent const& event, TAssocs const& assocs, TTracks const& tracks, TDileptons const& dileptons)
  {
    std::vector<int> daugDileptonGlobalIndexes;
    std::vector<float> vecPT, vecEta, vecPhi, vecMass, vecSign;
    std::vector<float> vecPTRef, vecEtaRef, vecPhiRef, vecITSChi2NCl, vecTPCNClsCR, vecTPCNClsFound, vecTPCChi2NCl;
    std::vector<uint8_t> vecITSClusterMap;
    std::vector<float> vecTPCSignal, vecTPCNSigmaEl, vecTPCNSigmaPi, vecTPCNSigmaPr, vecDcaXY, vecDcaZ;
    std::vector<float> vecPt1, vecEta1, vecPhi1, vecITSChi2NCl1, vecTPCNClsCR1, vecTPCNClsFound1, vecTPCChi2NCl1, vecTPCSignal1, vecTPCNSigmaEl1, vecTPCNSigmaPi1, vecTPCNSigmaPr1;
    std::vector<int> vecSign1;
    std::vector<float> vecPt2, vecEta2, vecPhi2, vecITSChi2NCl2, vecTPCNClsCR2, vecTPCNClsFound2, vecTPCChi2NCl2, vecTPCSignal2, vecTPCNSigmaEl2, vecTPCNSigmaPi2, vecTPCNSigmaPr2;
    std::vector<int> vecSign2;

    for (auto& dilepton : dileptons) {
      vecPT.push_back(dilepton.pt());
      vecEta.push_back(dilepton.eta());
      vecPhi.push_back(dilepton.phi());
      vecMass.push_back(dilepton.mass());
      vecSign.push_back(dilepton.sign());
      int indexLepton1 = dilepton.index0Id();
      int indexLepton2 = dilepton.index1Id();
      daugDileptonGlobalIndexes.push_back(indexLepton1);
      daugDileptonGlobalIndexes.push_back(indexLepton2);
      auto lepton1 = tracks.rawIteratorAt(indexLepton1);
      auto lepton2 = tracks.rawIteratorAt(indexLepton2);
      vecPt1.push_back(lepton1.pt());
      vecEta1.push_back(lepton1.eta());
      vecPhi1.push_back(lepton1.phi());
      vecSign1.push_back(lepton1.sign());
      vecITSChi2NCl1.push_back(lepton1.itsChi2NCl());
      vecTPCNClsCR1.push_back(lepton1.tpcNClsCrossedRows());
      vecTPCNClsFound1.push_back(lepton1.tpcNClsFound());
      vecTPCChi2NCl1.push_back(lepton1.tpcChi2NCl());
      vecTPCSignal1.push_back(lepton1.tpcSignal());
      vecTPCNSigmaEl1.push_back(lepton1.tpcNSigmaEl());
      vecTPCNSigmaPi1.push_back(lepton1.tpcNSigmaPi());
      vecTPCNSigmaPr1.push_back(lepton1.tpcNSigmaPr());
      vecPt2.push_back(lepton2.pt());
      vecEta2.push_back(lepton2.eta());
      vecPhi2.push_back(lepton2.phi());
      vecSign2.push_back(lepton2.sign());
      vecITSChi2NCl2.push_back(lepton2.itsChi2NCl());
      vecTPCNClsCR2.push_back(lepton2.tpcNClsCrossedRows());
      vecTPCNClsFound2.push_back(lepton2.tpcNClsFound());
      vecTPCChi2NCl2.push_back(lepton2.tpcChi2NCl());
      vecTPCSignal2.push_back(lepton2.tpcSignal());
      vecTPCNSigmaEl2.push_back(lepton2.tpcNSigmaEl());
      vecTPCNSigmaPi2.push_back(lepton2.tpcNSigmaPi());
      vecTPCNSigmaPr2.push_back(lepton2.tpcNSigmaPr());
    }

    for (auto& assoc : assocs) {
      if (!(assoc.isBarrelSelected_raw() & fConfigBarrelTrackCutBitMask.value)) {
        continue;
      }
      auto track = assoc.template reducedtrack_as<TTracks>();
      bool isDileptonDaug = false;
      for (int daugDileptonGlobalIndex : daugDileptonGlobalIndexes) {
        if (track.globalIndex() == daugDileptonGlobalIndex) {
          isDileptonDaug = true;
          break;
        }
      }
      if (isDileptonDaug) {
        continue;
      }
      vecPTRef.push_back(track.pt());
      vecEtaRef.push_back(track.eta());
      vecPhiRef.push_back(track.phi());
      vecITSChi2NCl.push_back(track.itsChi2NCl());
      vecTPCNClsCR.push_back(track.tpcNClsCrossedRows());
      vecTPCNClsFound.push_back(track.tpcNClsFound());
      vecTPCChi2NCl.push_back(track.tpcChi2NCl());
      vecITSClusterMap.push_back(track.itsClusterMap());
      vecTPCSignal.push_back(track.tpcSignal());
      vecTPCNSigmaEl.push_back(track.tpcNSigmaEl());
      vecTPCNSigmaPi.push_back(track.tpcNSigmaPi());
      vecTPCNSigmaPr.push_back(track.tpcNSigmaPr());
      vecDcaXY.push_back(track.dcaXY());
      vecDcaZ.push_back(track.dcaZ());
    }

    flowVectorsDetailed(event.multTPC(), event.multTracklets(), event.multNTracksPV(), event.multFT0C(), event.numContrib(), event.posX(), event.posY(), event.posZ(), event.selection_raw(), -999.f, vecPT, vecEta, vecPhi, vecMass, vecSign, vecPTRef, vecEtaRef, vecPhiRef, vecITSChi2NCl, vecTPCNClsCR, vecTPCNClsFound, vecTPCChi2NCl, vecITSClusterMap, vecTPCSignal, vecTPCNSigmaEl, vecTPCNSigmaPi, vecTPCNSigmaPr, vecDcaXY, vecDcaZ, vecPt1, vecEta1, vecPhi1, vecSign1, vecITSChi2NCl1, vecTPCNClsCR1, vecTPCNClsFound1, vecTPCChi2NCl1, vecTPCSignal1, vecTPCNSigmaEl1, vecTPCNSigmaPi1, vecTPCNSigmaPr1, vecPt2, vecEta2, vecPhi2, vecSign2, vecITSChi2NCl2, vecTPCNClsCR2, vecTPCNClsFound2, vecTPCChi2NCl2, vecTPCSignal2, vecTPCNSigmaEl2, vecTPCNSigmaPi2, vecTPCNSigmaPr2);
  }

  void processSkimmed(soa::Filtered<MyEventsVtxCovSelected> const& events, soa::Filtered<MyBarrelTrackAssocsSelected> const& assocs, MyBarrelTracksWithCov const& tracks, soa::Filtered<MyDielectronCandidates> const& dileptons)
  {
    for (auto& event : events) {
      auto groupedAssocs = assocs.sliceBy(perEventTrackAssocs, event.globalIndex());
      auto groupedDileptons = dileptons.sliceBy(perEventPairs, event.globalIndex());
      fillQVectors(event, groupedAssocs, tracks, groupedDileptons, true);
    }
  }

  void processSkimmedWithoutDilepton(soa::Filtered<MyEventsVtxCovSelected> const& events, soa::Filtered<MyBarrelTrackAssocsSelected> const& assocs, MyBarrelTracksWithCov const& tracks)
  {
    struct EmptyDilepton {
      int index0Id() const { return -1; }
      int index1Id() const { return -1; }
      float pt() const { return 0.f; }
      float eta() const { return 0.f; }
      float phi() const { return 0.f; }
      float mass() const { return 0.f; }
      float sign() const { return 0.f; }
    };
    std::vector<EmptyDilepton> emptyDileptons;
    for (auto& event : events) {
      auto groupedAssocs = assocs.sliceBy(perEventTrackAssocs, event.globalIndex());
      fillQVectors(event, groupedAssocs, tracks, emptyDileptons, false);
    }
  }

  void processSEFlowPairPoiRef(soa::Filtered<MyEventsVtxCovSelected> const& events, soa::Filtered<MyBarrelTrackAssocsSelected> const& assocs, MyBarrelTracksWithCov const& tracks, soa::Filtered<MyDielectronCandidates> const& dileptons)
  {
    for (auto& event : events) {
      if (!rctChecker(event)) {
        continue;
      }
      auto groupedAssocs = assocs.sliceBy(perEventTrackAssocs, event.globalIndex());
      auto groupedDileptons = dileptons.sliceBy(perEventPairs, event.globalIndex());
      fillSEFlowPairPoiRef(event, groupedAssocs, tracks, groupedDileptons);
    }
  }

  void processMEFlowPairRefOnly(soa::Filtered<MyEventsHashVtxCovSelected>& events, soa::Filtered<MyBarrelTrackAssocsSelected> const& assocs, MyBarrelTracksWithCov const& tracks)
  {
    events.bindExternalIndices(&assocs);
    for (auto& [event1, event2] : selfCombinations(hashBin, fConfigMixingDepthRR.value, -1, events, events)) {
      auto evAssocs1 = assocs.sliceBy(perEventTrackAssocs, event1.globalIndex());
      auto evAssocs2 = assocs.sliceBy(perEventTrackAssocs, event2.globalIndex());

      for (auto& assoc1 : evAssocs1) {
        if (!(assoc1.isBarrelSelected_raw() & fConfigBarrelTrackCutBitMask.value)) {
          continue;
        }
        auto track1 = assoc1.template reducedtrack_as<MyBarrelTracksWithCov>();
        for (auto& assoc2 : evAssocs2) {
          if (!(assoc2.isBarrelSelected_raw() & fConfigBarrelTrackCutBitMask.value)) {
            continue;
          }
          auto track2 = assoc2.template reducedtrack_as<MyBarrelTracksWithCov>();
          flowPairsRR(event1.multTPC(), event1.multTracklets(), event1.multNTracksPV(), event1.selection_raw(), event1.multFT0C(), event1.numContrib(), event1.posZ(), track1.pt(), track1.eta(), track1.phi(), track2.pt(), track2.eta(), track2.phi());
        }
      }
    }
  }

  void processMEFlowPairPoiRef(soa::Filtered<MyEventsHashVtxCovSelected>& events, soa::Filtered<MyBarrelTrackAssocsSelected> const& assocs, MyBarrelTracksWithCov const& tracks, soa::Filtered<MyDielectronCandidates> const& dileptons)
  {
    events.bindExternalIndices(&dileptons);
    events.bindExternalIndices(&assocs);
    for (auto& [event1, event2] : selfCombinations(hashBin, fConfigMixingDepthPR.value, -1, events, events)) {
      auto evDileptons = dileptons.sliceBy(perEventPairs, event1.globalIndex());
      auto evAssocs = assocs.sliceBy(perEventTrackAssocs, event2.globalIndex());

      for (auto dilepton : evDileptons) {
        for (auto& assoc : evAssocs) {
          if (!(assoc.isBarrelSelected_raw() & fConfigBarrelTrackCutBitMask.value)) {
            continue;
          }
          auto track = assoc.template reducedtrack_as<MyBarrelTracksWithCov>();
          flowPairsPR(event1.multTPC(), event1.multTracklets(), event1.multNTracksPV(), event1.selection_raw(), event1.multFT0C(), event1.numContrib(), event1.posZ(), dilepton.pt(), dilepton.eta(), dilepton.phi(), dilepton.mass(), dilepton.sign(), track.pt(), track.eta(), track.phi());
        }
      }
    }
  }

  void processDummy(MyEvents&)
  {
  }

  PROCESS_SWITCH(AnalysisFlow, processSkimmed, "Run associated table flow with dileptons", false);
  PROCESS_SWITCH(AnalysisFlow, processSkimmedWithoutDilepton, "Run associated table flow without dileptons", false);
  PROCESS_SWITCH(AnalysisFlow, processSEFlowPairPoiRef, "Run associated same-event POI-reference flow table", false);
  PROCESS_SWITCH(AnalysisFlow, processMEFlowPairRefOnly, "Run associated mixed-event reference-reference flow pairs", false);
  PROCESS_SWITCH(AnalysisFlow, processMEFlowPairPoiRef, "Run associated mixed-event POI-reference flow pairs", false);
  PROCESS_SWITCH(AnalysisFlow, processDummy, "Dummy function", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<AnalysisFlow>(cfgc)};
}

std::vector<int> findCommonBits(int num1, int num2)
{
  std::vector<int> commonBits;
  std::bitset<32> bitset1(num1);
  std::bitset<32> bitset2(num2);
  for (int i = 0; i < 32; ++i) {
    if (bitset1[i] && bitset2[i]) {
      commonBits.push_back(i);
    }
  }
  return commonBits;
}

std::vector<TString> splitTString(TString str, char c)
{
  std::vector<TString> res;
  TString tmp = "";
  for (int i = 0; i < str.Length(); i++) {
    if (str[i] == c) {
      res.push_back(tmp);
      tmp = "";
    } else {
      tmp += str[i];
    }
  }
  if (tmp != "") {
    res.push_back(tmp);
  }
  return res;
}
