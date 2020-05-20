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
 *  @file BeamParams.h
 */

#ifndef BEAMPARAMS_H
#define BEAMPARAMS_H 1

#include <G4PhysicalConstants.hh>
#include <G4SystemOfUnits.hh>
#include <G4ThreeVector.hh>
#include <globals.hh>

class BeamParams
{
public:
  BeamParams();
  /// \param x is the energy expressed in keV
  void SetEnergy(G4double x) {energy = x * keV;};
  void SetVtxPosition(G4double x, G4double y, G4double z);
  void SetPolarization(G4double x, G4double y, G4double z);
  void SetMomentum(G4double x, G4double y, G4double z);
  void SetVtxPosition(G4ThreeVector v);
  void SetPolarization(G4ThreeVector v);
  void SetMomentum(G4ThreeVector v);

  G4ThreeVector GetVtx() const {return vtxPosition;};
  G4double GetEnergy() const {return energy;};
  G4ThreeVector GetMomentum() const {return momentum;};
  G4ThreeVector GetPolarization() const {return polarization;};

private:
  G4double energy;
  G4ThreeVector vtxPosition;
  G4ThreeVector polarization;
  G4ThreeVector momentum;
};

#endif
