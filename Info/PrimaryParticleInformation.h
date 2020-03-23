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
 *  @file PrimaryParticleInformation.h
 */

#ifndef PRIMARY_PARTICLE_INFORMATION_H
#define PRIMARY_PARTICLE_INFORMATION_H 1

#include <G4VUserPrimaryParticleInformation.hh>
#include <G4ThreeVector.hh>
#include <globals.hh>

class PrimaryParticleInformation : public G4VUserPrimaryParticleInformation
{
public:
  PrimaryParticleInformation();
  virtual ~PrimaryParticleInformation();
  void Clear();
  virtual void Print() const;

  G4int GetIndex() {return fIndex;}
  G4int GetGammaMultiplicity() {return fDecayMultiplicity;}
  G4int GetGeneratedGammaMultiplicity() {return fGeneratedMultiplicity;}
  void SetGeneratedGammaMultiplicity(G4int i) {fGeneratedMultiplicity = i;}
  void SetIndex(G4int i) {fIndex = i;}
  void SetGammaMultiplicity(G4int i) { fDecayMultiplicity = i; }
  void SetGenMomentum(G4double x, G4double y, G4double z);
  G4ThreeVector GenGenMomentum() {return fGenMomentum;}

  //! Multiplicity flags
  static const G4int kBackground = 0;
  static const G4int kPromptGamma = 1;
  static const G4int kBackToBackGamma = 2;
  static const G4int koPsGamma = 3;

private:
  G4int fIndex;
  G4int fDecayMultiplicity;
  G4int fGeneratedMultiplicity;
  G4ThreeVector fGenMomentum;
};

#endif
