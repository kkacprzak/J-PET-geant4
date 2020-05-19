/**
 *  @copyright Copyright 2020 The J-PET Monte Carlo Authors. All rights reserved.
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may find a copy of the License in the LICENCE file.
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *  @file HistoManager.h
 */

#ifndef HISTOMANAGER_H
#define HISTOMANAGER_H 1

#include "../Objects/Framework/JPetGeantEventInformation.h"
#include "../Objects/Framework/JPetGeantEventPack.h"
#include "../Objects/Framework/JPetGeantDecayTree.h"
#include "../Objects/Framework/JPetGeantScinHits.h"
#include "../Objects/Geant4/DetectorHit.h"
#include "../Info/VtxInformation.h"
#include "../Info/EventMessenger.h"
#include <G4PrimaryParticle.hh>
#include <G4Event.hh>
#include <globals.hh>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>

class TFile;
class TTree;
const int MaxHisto = 7;
const int MaxHisto2D = 9;

/**
 * @class HistoManager
 * @brief class reach for informations stored during sumulations
 * and saves them into ROOT tree structures
 */
class HistoManager
{
public:
  HistoManager();
  ~HistoManager();

  //! call once; book (create) all trees and histograms
  void Book();
  //! call once; save all trees and histograms
  void Save();

  void SaveEvtPack()
  {
    fTree->Fill();
  };

  void Clear()
  {
    fEventPack->Clear();
  };

  void AddGenInfo(VtxInformation* info);
  void AddGenInfoParticles(G4PrimaryParticle* particle);
  void AddNewHit(DetectorHit*);

  void SetEventNumber(int x)
  {
    fEventPack->SetEventNumber(x);
  };

  void SetHistogramCreation(bool tf)
  {
    fMakeControlHisto = tf;
  };

  bool MakeControlHisto()
  {
    return fMakeControlHisto;
  };

  void FillHistoGenInfo(const G4Event* anEvent);

  void FillCosmicInfo(G4double theta, G4ThreeVector init, G4ThreeVector vtx);

  const JPetGeantEventInformation* GetGeantInfo()
  {
    return fGeantInfo;
  }

private:
  bool fMakeControlHisto = false;

  /**
   * 1D histograms array:
   * 0 - generated gamma
   * 1 - time
   * 2 - energy
   * 3 - z position
   * 4 - lifetime
   * 5 - prompt lifetime
   */
  TH1F* fHisto[MaxHisto];

  /**
  * 2D histograms array:
  * 0 - XY of hits
  * 1 - XY of annihilation
  * 2 - XZ annihilation
  * 3 - YZ annihilation
  * 4 - XY prompt
  * 5 - XZ prompt
  * 6 -YZ prompt
  */
  TH2F* fHisto2D[MaxHisto2D];

  // cosmic generator histograms
  TH1F *h_theta;
  TH2F *h_init_xy, *h_init_xz, *h_init_yz, *h_orig_xy, *h_orig_xz, *h_orig_yz;

  TFile* fRootFile;
  TTree* fTree;
  TBranch* fBranchTrk;
  TBranch* fBranchScin;
  TBranch* fBranchEventPack;

  JPetGeantEventPack* fEventPack;
  JPetGeantEventInformation* fGeantInfo;
  EventMessenger* fEvtMessenger = EventMessenger::GetEventMessenger();

  void BookHistograms();
};

#endif
