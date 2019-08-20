/**
 *  @copyright Copyright 2019 The J-PET Monte Carlo Authors. All rights reserved.
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
 *  @file SteppingAction.cpp
 */

#include "../Info/PrimaryParticleInformation.h"
#include <G4TransportationManager.hh>
#include <G4SteppingManager.hh>
#include <G4PrimaryParticle.hh>
#include "SteppingAction.h"
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "../Info/EventMessenger.h"

SteppingAction::SteppingAction()
{
  G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking()->SetPushVerbosity(0);
}

SteppingAction::~SteppingAction() {}

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{

  // kill event if primary particle escapes the world volume
  if (EventMessenger::GetEventMessenger()->KillEventsEscapingWorld()) {
    G4StepPoint* point = aStep->GetPostStepPoint();
    if (point->GetStepStatus() == G4StepStatus::fWorldBoundary) {
      if (aStep->GetTrack()->GetParentID() == 0) {
        PrimaryParticleInformation* info  = dynamic_cast<PrimaryParticleInformation*> (aStep->GetTrack()->GetDynamicParticle()->GetPrimaryParticle()->GetUserInformation());
        if (info != nullptr ) {
          G4int multiplicity = info->GetGammaMultiplicity();
          const G4int minBoundMultiplicity = EventMessenger::GetEventMessenger()->GetMinRegMultiplicity();
          const G4int maxBoundMultiplicity = EventMessenger::GetEventMessenger()->GetMaxRegMultiplicity();
          const G4int excludedMultiplicity = EventMessenger::GetEventMessenger()->GetExcludedMultiplicity();
          if(  (multiplicity >= minBoundMultiplicity) && (multiplicity <= maxBoundMultiplicity) && (multiplicity != excludedMultiplicity) ) {
            G4RunManager::GetRunManager()->AbortEvent();        
          }
        }
      }
    } 
  }


  //! execute code for
  //! - generated by user gamma quanta
  //! - physical effects that does not occur in Sensitive Detector

  if ("Transportation" == aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()) return;
  if (aStep->GetTrack()->GetParentID() != 0 ) return;
  if (nullptr != aStep->GetPostStepPoint()->GetPhysicalVolume()->GetLogicalVolume()->GetSensitiveDetector()) return;

  PrimaryParticleInformation* info  = static_cast<PrimaryParticleInformation*> (aStep->GetTrack()->GetDynamicParticle()->GetPrimaryParticle()->GetUserInformation());
  if (info != 0 ) {
    //! particle quanta interact in phantom or frame (but not SD!)
    info->SetGammaMultiplicity(PrimaryParticleInformation::kBackground);
  }
}
