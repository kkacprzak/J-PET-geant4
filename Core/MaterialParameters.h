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
 *  @file MaterialParameters.h
 */

#ifndef MATERIALPARAMETERS_H
#define MATERIALPARAMETERS_H 1

#include <G4PhysicalConstants.hh>
#include <G4SystemOfUnits.hh>
#include <G4ThreeVector.hh>
#include <globals.hh>
#include <vector>

struct MaterialConstants
{
  MaterialConstants(const std::vector<G4double> & oPsLF, const std::vector<G4double> & oPsProb, 
                    G4double pPsLF, G4double pPsFrac, 
                    const std::vector<G4double> & DirLF, const std::vector<G4double> & DirProb);
  std::vector<G4double> oPsLifetimes;
  std::vector<G4double> oPsProbabilities;
  G4double pPsLifetime;
  G4double pPsFraction;
  std::vector<G4double> directLifetimes;  
  std::vector<G4double> directProbabilities;
};

class MaterialParameters
{
public:
  enum MaterialID { mUnknown, mXAD4, mAl, mKapton, mPlexiglass, mScin, mPA6, mAir};
    
  static const MaterialConstants XAD4;
  static const MaterialConstants Al;
  static const MaterialConstants Kapton;
  static const MaterialConstants Plexiglass; // PMMA
  static const MaterialConstants Scin;
  static const MaterialConstants PA6;
  static const MaterialConstants Air;
  static MaterialConstants Temp;
  
  MaterialParameters();
  MaterialParameters(const std::vector<G4double> & oPsLF, const std::vector<G4double> & oPsProb, 
                      G4double pPsLF, G4double pPsFrac, 
                      const std::vector<G4double> & DirectLF, const std::vector<G4double> & DirectProb);
  static const G4double oPsTauVaccum;
  static const G4double directTauMax;
  static const G4double pPsTauVaccum;
  static const G4double pPsTauMax;
  static const G4double tauBulk;
  static const G4double direct3Gfraction;
//! doi 10.2478/umcschem-2013-0012          1/372
  static const G4double sodiumGammaEnergy;
  static const G4double sodiumGammaTau;
  static G4String AnnihlationMode;
  
  static void ClearTemp();
  static void SetAnnihilationMode(G4String mode);
  
  void SetMaterial(MaterialConstants mat);
  void SetMaterialByName(MaterialID materialID);
  void AddoPsComponent(G4double Lifetime, G4double Probability);
  void SetoPsComponents(const std::vector<G4double> & oPsLF, const std::vector<G4double> & oPsProb);
  void SetpPsComponent(G4double Lifetime, G4double Fraction);
  void AddDirectComponent(G4double Lifetime, G4double Probability);
  void SetDirectComponents(const std::vector<G4double> & DirectLF, const std::vector<G4double> & DirectProb);

  G4double GetoPsIntensity2G(G4double oPsLifetime, G4double oPsProbability);
  G4double GetoPsIntensity3G(G4double oPsLifetime, G4double oPsProbability);
  void SetComponentsIntensities();
  void RenormalizeIntensities(G4double TotalIntensity);
  
  G4double GetoPs2GTotalIntensity();
  G4double GetoPs3GTotalIntensity();
  G4double GetDirect2GTotalIntensity();
  G4double GetDirect3GTotalIntensity();
  G4double GetpPs2GTotalIntensity();
  G4double GetoPs2GLifetimeFromVector(double randNumber);
  G4double GetoPs3GLifetimeFromVector(double randNumber);
  G4double GetDirectLifetimeFromVector(double randNumber);
  G4double GetpPsLifetime();
  G4double GetLifetimeVector(std::vector<G4double> vectorToCheck);

private:
//There can be more than 1 component of oPs so that is why vector is used
//Different types of pores/free volumes can give different mean lifetime of o-Ps in them
  std::vector<G4double> foPsLifetimes;
//Probability of positron decay with given lifetime. 
//From Probabilities intensities of decay into 2/3 for each component, can be calculated
  std::vector<G4double> foPsProbabilities;
  std::vector<G4double> foPs2GIntensities;
  std::vector<G4double> foPs3GIntensities;
  
  G4double fpPsLifetime;
//pPs Intensity defined as a fraction of pPs to oPs
  G4double fpPsFraction;
  G4double fpPsIntensity;
  
//There can be more than 1 component of direct annihilation of positron so that is why vector is used
//Different defects in material can give different mean lifetime of positron in it
  std::vector<G4double> fdirectLifetimes;
  std::vector<G4double> fdirectProbabilities;
  std::vector<G4double> fdirectIntensities;
};

#endif
