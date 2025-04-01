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
#include <iostream>
#include <vector>
#include <algorithm>
#include <TH1F.h>
#include <TH3F.h>
#include <THashList.h>
#include <TList.h>
#include <TString.h>
#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "PWGDQ/DataModel/ReducedInfoTables.h"
#include "PWGDQ/Core/VarManager.h"
#include "PWGDQ/Core/HistogramManager.h"
#include "PWGDQ/Core/MixingHandler.h"
#include "PWGDQ/Core/AnalysisCut.h"
#include "PWGDQ/Core/AnalysisCompositeCut.h"
#include "PWGDQ/Core/HistogramsLibrary.h"
#include "PWGDQ/Core/CutsLibrary.h"
#include "PWGDQ/Core/MixingLibrary.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "Field/MagneticField.h"
#include "TGeoGlobalMagField.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "ITSMFTBase/DPLAlpideParam.h"
#include "Common/CCDB/EventSelectionParams.h"

using std::cout;
using std::endl;
using std::string;

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod;

// Some definitions
namespace o2::aod
{

namespace dqanalysisflags
{
// TODO: the barrel amd muon selection columns are bit maps so unsigned types should be used, however, for now this is not supported in Filter expressions
// TODO: For now in the tasks we just statically convert from unsigned int to int, which should be fine as long as we do
//      not use a large number of bits (>=30)
DECLARE_SOA_COLUMN(MixingHash, mixingHash, int);
DECLARE_SOA_COLUMN(IsEventSelected, isEventSelected, int);
DECLARE_SOA_COLUMN(IsBarrelSelected, isBarrelSelected, int);
DECLARE_SOA_COLUMN(IsMuonSelected, isMuonSelected, int);
DECLARE_SOA_COLUMN(IsBarrelSelectedPrefilter, isBarrelSelectedPrefilter, int);
DECLARE_SOA_COLUMN(IsPrefilterVetoed, isPrefilterVetoed, int);
} // namespace dqanalysisflags

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
} // namespace flowVec

DECLARE_SOA_TABLE(EventCuts, "AOD", "DQANAEVCUTS", dqanalysisflags::IsEventSelected);
DECLARE_SOA_TABLE(MixingHashes, "AOD", "DQANAMIXHASH", dqanalysisflags::MixingHash);
DECLARE_SOA_TABLE(BarrelTrackCuts, "AOD", "DQANATRKCUTS", dqanalysisflags::IsBarrelSelected, dqanalysisflags::IsBarrelSelectedPrefilter);
DECLARE_SOA_TABLE(MuonTrackCuts, "AOD", "DQANAMUONCUTS", dqanalysisflags::IsMuonSelected);
DECLARE_SOA_TABLE(Prefilter, "AOD", "DQPREFILTER", dqanalysisflags::IsPrefilterVetoed);
DECLARE_SOA_TABLE(FlowVecs, "AOD", "DQFLOWVECS", evsel::Selection, flowVec::MultFT0C, flowVec::MultVtxContri, flowVec::VtxZ, flowVec::PT, flowVec::Eta, flowVec::Phi, flowVec::Mass, flowVec::Sign, flowVec::QvecRe, flowVec::QvecIm, flowVec::QvecAmp);
} // namespace o2::aod

// Declarations of various short names
using MyEvents = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended>;
using MyEventsSelected = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::EventCuts>;
using MyEventsHashSelected = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::EventCuts, aod::MixingHashes>;
using MyEventsVtxCov = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsVtxCov>;
using MyEventsVtxCovSelected = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsVtxCov, aod::EventCuts>;
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
  int fCurrentRun; // needed to detect if the run changed and trigger update of calibrations etc.
  Configurable<float> fConfigDileptonLowMass{"cfgDileptonLowMass", 0.0, "Low mass cut for the dileptons used in analysis"};
  Configurable<float> fConfigDileptonHighMass{"cfgDileptonHighMass", 5.0, "High mass cut for the dileptons used in analysis"};
  Configurable<std::string> fConfigQVectorsBarrel2Cal{"cfgQVectorsBarrel", "4,1;2,1:0,1", "Semicolon separated list of q vector to be calculated"};
  Configurable<int> fConfigBarrelTrackCutBitMask{"cfgBarrelTrackCutBitMask", 1, "bit mask of the track cuts"};
  Produces<aod::FlowVecs> qVectors;

  Filter eventFilter = aod::dqanalysisflags::isEventSelected == 1;
  Filter dileptonFilter = aod::reducedpair::mass > fConfigDileptonLowMass&& aod::reducedpair::mass < fConfigDileptonHighMass;
  Filter filterBarrelTrackSelected = aod::dqanalysisflags::isBarrelSelected > 0;

  std::vector<float> qvecRe;
  std::vector<float> qvecIm;
  std::vector<float> qvecAmp;
  std::vector<std::vector<std::array<int, 2>>> cfgsQVecBarrel2Cal;
  std::vector<int> cfgsIndexCutBitMask;

  int nCfgQvectorBarrel2Cal = 0;

  constexpr static uint32_t fgDileptonFillMap = VarManager::ObjTypes::ReducedTrack | VarManager::ObjTypes::Pair; // fill map
  constexpr static uint32_t candidateTypeDilepton = VarManager::PairCandidateType::kDecayToEE;
  constexpr static uint32_t eventFillMapDilepton = VarManager::ObjTypes::ReducedEventVtxCov;
  constexpr static uint32_t fillMapDilepton = VarManager::ObjTypes::ReducedTrackBarrelCov;

  // float* fValuesDilepton;
  // float* fValuesHadron;
  HistogramManager* fHistMan;
  bool isWithDilepton = false;

  void init(o2::framework::InitContext& context)
  {
    fCurrentRun = 0;
    // fValuesDilepton = new float[VarManager::kNVars];
    // fValuesHadron = new float[VarManager::kNVars];
    VarManager::SetDefaultVarNames();
    fHistMan = new HistogramManager("analysisHistos", "aa", VarManager::kNVars);
    fHistMan->SetUseDefaultVariableNames(kTRUE);
    fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);

    // if (context.mOptions.get<bool>("processSkimmed")) {
    //   DefineHistograms(fHistMan, "DileptonsSelected;DileptonHadronInvMass;DileptonHadronCorrelationSE", fConfigAddDileptonHadHistogram); // define all histograms
    // }
    // if (context.mOptions.get<bool>("processMixedEvent")) {
    //   DefineHistograms(fHistMan, "DileptonHadronInvMassME;DileptonHadronCorrelationME", fConfigAddDileptonHadHistogram); // define all histograms
    // }

    if (context.mOptions.get<bool>("processSkimmed")) {
      isWithDilepton = true;
    }

    VarManager::SetUseVars(fHistMan->GetUsedVars());
    fOutputList.setObject(fHistMan->GetMainHistogramList());

    TString str_configuration_qVector2cal = fConfigQVectorsBarrel2Cal.value;
    std::vector<TString> vec_str_configuration_qVector2cal = splitTString(str_configuration_qVector2cal, ':');

    for (auto& str_cfg : vec_str_configuration_qVector2cal) {
      std::vector<TString> vec_str_cfg = splitTString(str_cfg, ';');
      std::vector<std::array<int, 2>> vec_cfg;
      for (auto& str : vec_str_cfg) {
        std::vector<TString> vec_str = splitTString(str, ',');
        if (vec_str.size() != 2) {
          LOG(fatal) << "Invalid configuration for q vector: " << str_cfg;
          exit(1);
        }
        vec_cfg.push_back({vec_str[0].Atoi(), vec_str[1].Atoi()});
        nCfgQvectorBarrel2Cal++;
      }
      cfgsQVecBarrel2Cal.push_back(vec_cfg);
    }

    LOG(info) << "Number of q vectors to calculate: " << nCfgQvectorBarrel2Cal;

    int num_TrackCuts = 0;
    for (int i = 0; i < 32; i++) {
      if (fConfigBarrelTrackCutBitMask.value & (1 << i)) {
        num_TrackCuts++;
      }
    }

    // push all bit in the bit mask to a vector cfgsIndexCutBitMask
    for (int i = 0; i < 32; i++) {
      if (fConfigBarrelTrackCutBitMask.value & (1 << i)) {
        cfgsIndexCutBitMask.push_back(i);
      }
    }

    if (num_TrackCuts != cfgsQVecBarrel2Cal.size()) {
      LOG(fatal) << "Number of track cuts and number of q vector sets to calculate do not match";
      exit(1);
    }
  }

  template <int TCandidateType, uint32_t TEventFillMap, uint32_t TTrackFillMap, typename TEvent, typename TTracks>
  void runFlow(TEvent const& event, TTracks const& tracks, soa::Filtered<MyDielectronCandidates> const& dileptons)
  {
    std::vector<int> trackGlobalIndexes;
    std::vector<int> daugDileptonGlobalIndexes;

    qvecRe.resize(nCfgQvectorBarrel2Cal);
    qvecIm.resize(nCfgQvectorBarrel2Cal);
    qvecAmp.resize(nCfgQvectorBarrel2Cal);

    for (int i = 0; i < nCfgQvectorBarrel2Cal; ++i) {
      qvecRe[i] = 0;
      qvecIm[i] = 0;
      qvecAmp[i] = 0;
    }

    if (isWithDilepton)
      for (auto dilepton : dileptons) {
        int indexLepton1 = dilepton.index0Id();
        int indexLepton2 = dilepton.index1Id();
        daugDileptonGlobalIndexes.push_back(indexLepton1);
        daugDileptonGlobalIndexes.push_back(indexLepton2);
      }

    for (auto& hadron : tracks) {
      if (!(uint32_t(hadron.isBarrelSelected()) & fConfigBarrelTrackCutBitMask.value))
        continue;
      bool isDileptonDaug = false;
      for (int daugDileptonGlobalIndexe : daugDileptonGlobalIndexes) {
        if (hadron.globalIndex() == daugDileptonGlobalIndexe) {
          isDileptonDaug = true;
          break;
        }
      }
      if (isDileptonDaug)
        continue;

      float weff = 1.0, wacc = 1.0;
      std::vector<int> vec_Bits;
      vec_Bits = findCommonBits(hadron.isBarrelSelected(), fConfigBarrelTrackCutBitMask.value);
      for (auto Bit : vec_Bits) {
        int index = 0;
        std::vector<std::vector<std::array<int, 2>>> cfgsQVec2Cal;
        int indexBit = 0;
        // cfgsIndexCutBitMask
        for (int i = 0; i < Bit; ++i) {
          if (cfgsIndexCutBitMask[i] == Bit) {
            indexBit = i;
            break;
          }
        }

        for (int i = 0; i < indexBit; ++i)
          index += cfgsQVecBarrel2Cal[i].size();
        cfgsQVec2Cal = cfgsQVecBarrel2Cal;

        for (int iQvector = 0; iQvector < cfgsQVec2Cal[indexBit].size(); ++iQvector) {
          int npow_phi = cfgsQVec2Cal[indexBit][iQvector][0];
          int npow_weight = cfgsQVec2Cal[indexBit][iQvector][1];
          double weights = 1;
          for (int i = 0; i < npow_weight; ++i) {
            weights *= 1. / weff / wacc;
          }
          double phi = hadron.phi();
          qvecRe[index + iQvector] += weights * std::cos(npow_phi * phi);
          qvecIm[index + iQvector] += weights * std::sin(npow_phi * phi);
          qvecAmp[index + iQvector] += weights * weights;
        }
      }
    }

    std::vector<float> vecPT;
    std::vector<float> vecEta;
    std::vector<float> vecPhi;
    std::vector<float> vecMass;
    std::vector<float> vecSign;
    if (isWithDilepton)
      for (auto dilepton : dileptons) {
        vecPT.push_back(dilepton.pt());
        vecEta.push_back(dilepton.eta());
        vecPhi.push_back(dilepton.phi());
        vecMass.push_back(dilepton.mass());
        vecSign.push_back(dilepton.sign());
      }

    qVectors(event.selection_raw(), event.multFT0C(), event.numContrib(), event.posZ(), vecPT, vecEta, vecPhi, vecMass, vecSign, qvecRe, qvecIm, qvecAmp);
  }

  template <int TCandidateType, uint32_t TEventFillMap, uint32_t TTrackFillMap, typename TEvent, typename TTracks>
  void runFlow2(TEvent const& event, TTracks const& tracks)
  {
    qvecRe.resize(nCfgQvectorBarrel2Cal);
    qvecIm.resize(nCfgQvectorBarrel2Cal);
    qvecAmp.resize(nCfgQvectorBarrel2Cal);

    for (int i = 0; i < nCfgQvectorBarrel2Cal; ++i) {
      qvecRe[i] = 0;
      qvecIm[i] = 0;
      qvecAmp[i] = 0;
    }

    for (auto& hadron : tracks) {
      if (!(uint32_t(hadron.isBarrelSelected()) & fConfigBarrelTrackCutBitMask.value))
        continue;

      float weff = 1.0, wacc = 1.0;
      std::vector<int> vec_Bits;
      vec_Bits = findCommonBits(hadron.isBarrelSelected(), fConfigBarrelTrackCutBitMask.value);
      for (auto Bit : vec_Bits) {
        int index = 0;
        std::vector<std::vector<std::array<int, 2>>> cfgsQVec2Cal;
        int indexBit = 0;
        // cfgsIndexCutBitMask
        for (int i = 0; i < Bit; ++i) {
          if (cfgsIndexCutBitMask[i] == Bit) {
            indexBit = i;
            break;
          }
        }

        for (int i = 0; i < indexBit; ++i)
          index += cfgsQVecBarrel2Cal[i].size();
        cfgsQVec2Cal = cfgsQVecBarrel2Cal;

        for (int iQvector = 0; iQvector < cfgsQVec2Cal[indexBit].size(); ++iQvector) {
          int npow_phi = cfgsQVec2Cal[indexBit][iQvector][0];
          int npow_weight = cfgsQVec2Cal[indexBit][iQvector][1];
          double weights = 1;
          for (int i = 0; i < npow_weight; ++i) {
            weights *= 1. / weff / wacc;
          }
          double phi = hadron.phi();
          qvecRe[index + iQvector] += weights * std::cos(npow_phi * phi);
          qvecIm[index + iQvector] += weights * std::sin(npow_phi * phi);
          qvecAmp[index + iQvector] += weights * weights;
        }
      }
    }

    std::vector<float> vecPT;
    std::vector<float> vecEta;
    std::vector<float> vecPhi;
    std::vector<float> vecMass;
    std::vector<float> vecSign;
    qVectors(event.selection_raw(), event.multFT0C(), event.numContrib(), event.posZ(), vecPT, vecEta, vecPhi, vecMass, vecSign, qvecRe, qvecIm, qvecAmp);
  }

  void processSkimmed(soa::Filtered<MyEventsVtxCovSelected>::iterator const& event, MyBarrelTracksSelectedWithCov const& tracks, soa::Filtered<MyDielectronCandidates> const& dileptons)
  {
    runFlow<VarManager::kBtoJpsiEEK, gkEventFillMapWithCov, gkTrackFillMapWithCov>(event, tracks, dileptons);
  }

  void processSkimmedWithoutDilepton(soa::Filtered<MyEventsVtxCovSelected>::iterator const& event, MyBarrelTracksSelectedWithCov const& tracks)
  {
    runFlow2<VarManager::kBtoJpsiEEK, gkEventFillMapWithCov, gkTrackFillMapWithCov>(event, tracks);
  }

  void processDummy(MyEvents&)
  {
  }

  PROCESS_SWITCH(AnalysisFlow, processSkimmed, "2024.9.30", false);
  PROCESS_SWITCH(AnalysisFlow, processSkimmedWithoutDilepton, "processSkimmedWithoutDilepton", false);
  PROCESS_SWITCH(AnalysisFlow, processDummy, "Dummy function", false);
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
