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
 *  @file DetectorConstruction.h
 */

#ifndef DETECTORCONSTRUCTION_H
#define DETECTORCONSTRUCTION_H 1

#include <G4VUserDetectorConstruction.hh>
#include <G4VPhysicalVolume.hh>
#include <G4GeometryManager.hh>
#include "MaterialExtension.h"
#include <G4SystemOfUnits.hh>
#include <G4MaterialTable.hh>
#include <G4LogicalVolume.hh>
#include <G4VisAttributes.hh>
#include <G4NistManager.hh>
#include <G4PVPlacement.hh>
#include <G4SDManager.hh>
#include <G4Material.hh>
#include <G4Element.hh>
#include "DetectorSD.h"
#include <G4Colour.hh>
#include <CADMesh.hh>
#include <G4Cache.hh>
#include <globals.hh>
#include <G4Box.hh>
#include <string>

//#include "G4GDMLParser.hh"
class DetectorConstructionMessenger;

//! Flag for debugging purposes
const  G4bool checkOverlaps = false;

/**
* @class DetectorConstruction
* @brief creating detector; can read the CAD geometry
*/
class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  //! only single instance can exist
  static DetectorConstruction* GetInstance();
  virtual G4VPhysicalVolume* Construct();
  virtual void ConstructSDandField();
  void LoadGeometryForRun(G4int nr);
  G4int ReturnNumberOfScintillators();
  void UpdateGeometry();
  void ReloadMaterials();

  void LoadFrame(G4bool tf)
  {
    fLoadCADFrame = tf;
  };

  //! Modular layer (known as 4th layer); 24 modules filled with scintillators
  void ConstructModularLayer(bool tf)
  {
    fLoadModularLayer = tf;
  }

  G4int GetRunNumber()
  {
    return fRunNumber;
  };

private:
  static G4ThreadLocal G4bool fConstructedSDandField;
  static DetectorConstruction* fInstance;

  DetectorConstruction();
  virtual ~DetectorConstruction();
  DetectorConstructionMessenger* fMessenger;

  //! Load materials from NIST database
  void InitializeMaterials();
  //! Load detector elements from CAD files
  void ConstructFrameCAD();
  //! Create scintillators only; dimensions are right now fixed in code
  void ConstructScintillators();
  //! Construct modular layer inserted into detector (refered as 4th layer)
  void ConstructScintillatorsModularLayer();
  //! Create target used in run3 - big chamber no XAD material inside
  void ConstructTargetRun3();
  //! Create target used in run5 (and run4) - small chamber + XAD material inside
  void ConstructTargetRun5();
  //! Create target used in run6 - big chamber + XAD material inside tube
  void ConstructTargetRun6();
  //! Create target used in run7
  void ConstructTargetRun7();

  //! Corresponds to JPET measurements; run 0 = user setup
  G4int fRunNumber;
  //! Flag for loading frame from CAD file
  G4bool fLoadCADFrame;
  //! Flag for loading wrapping tf the scintillators
  G4bool fLoadWrapping;
  //! Flag for loading modular (4th) layer
  G4bool fLoadModularLayer;

  G4Box* fWorldSolid;
  G4LogicalVolume* fWorldLogical;
  G4VPhysicalVolume* fWorldPhysical;


  MaterialExtension* fAir;
  MaterialExtension* fKapton;
  MaterialExtension* fVacuum;
  MaterialExtension* fPlexiglass;
  MaterialExtension* fXADMaterial;
  MaterialExtension* fScinMaterial;
  MaterialExtension* fAluminiumMaterial;
  MaterialExtension* fSmallChamberMaterial;
  MaterialExtension* fSmallChamberRun7Material;

  G4LogicalVolume* fScinLog;
  G4LogicalVolume* fScinLogInModule;
  G4Cache<DetectorSD*> fDetectorSD;

//GDML output
    
  //  G4String fWriteFile;
  //  G4GDMLParser* fParser;

};

#endif
