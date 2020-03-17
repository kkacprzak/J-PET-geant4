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
 *  @file JPetMC.cpp
 */

#include "Actions/ActionInitialization.h"
#include "Core/DetectorConstruction.h"
#include "Core/PhysicsList.h"
#include <G4VisExecutive.hh>
#include "Core/RunManager.h"
#include <G4UIExecutive.hh>
#include <G4INCLRandom.hh>
#include <G4UImanager.hh>
#include "Info/EventMessenger.h"
#include <fstream> 
#include <random>

void setRandomSeed()
{
  long seeds[2];
  
  std::uniform_int_distribution<long> d(0, LONG_MAX);
  std::random_device rd1;

  seeds[0] = d(rd1);
  seeds[1] = d(rd1);
  G4Random::setTheSeeds(seeds);

  if (EventMessenger::GetEventMessenger()->SaveSeed()) {
    std::ofstream file;
    file.open ("seed", std::ofstream::out | std::ofstream::app);
    file << seeds[0] << "\n" << seeds[1] << "\n";
    file.close();
  }
}


int main (int argc, char** argv)
{

  G4UIExecutive* ui = 0;
  if (argc == 1) {
    ui = new G4UIExecutive(argc, argv);
  }

  RunManager* runManager = new RunManager;
  runManager->SetUserInitialization(DetectorConstruction::GetInstance());
  runManager->SetUserInitialization(new PhysicsList);
  runManager->SetUserInitialization(new ActionInitialization);

  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();

  if (!ui) {
    //! batch mode
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command + fileName);
  } else {
    //! interactive mode
    UImanager->ApplyCommand("/control/execute init_vis.mac");
    //UImanager->ApplyCommand("/control/execute /Users/sushilsharma/Geant4-xcode/project/Trento_Geant4/build/bin/init_vis.mac");
      ui->SessionStart();
    delete ui;
  }

  if (EventMessenger::GetEventMessenger()->SetRandomSeed()) {
    setRandomSeed();
  }


  delete visManager;
  delete runManager;
  return 0;
}


