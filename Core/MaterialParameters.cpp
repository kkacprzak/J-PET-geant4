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
 *  @file MaterialParameters.cpp
 */

#include "MaterialParameters.h"
#include "MaterialExtension.h"

MaterialConstants::MaterialConstants(const std::vector<G4double> & oPsLF, const std::vector<G4double> & oPsProb, 
                                        G4double pPsLF, G4double pPsFrac, 
                                        const std::vector<G4double> & DirLF, const std::vector<G4double> & DirProb)
{
  oPsLifetimes = oPsLF;
  oPsProbabilities = oPsProb;
  pPsLifetime = pPsLF;
  pPsFraction = pPsFrac;
  directLifetimes = DirLF;
  directProbabilities = DirProb;
}

const G4double MaterialParameters::oPsTauVaccum =  142.0 * ns;
const G4double MaterialParameters::directTauMax = 0.6 * ns;
const G4double MaterialParameters::pPsTauMax = 0.2 * ns;
const G4double MaterialParameters::tauBulk = 0.3;
const G4double MaterialParameters::direct3Gfraction =  0.002688;
const G4double MaterialParameters::sodiumGammaEnergy = 1.2770 * MeV;
const G4double MaterialParameters::sodiumGammaTau = 3.7 * ps;

G4String MaterialParameters::AnnihlationMode = "";
MaterialConstants MaterialParameters::Temp = MaterialConstants( {}, {}, 0., 0., {}, {} );
//type to probability

const MaterialConstants MaterialParameters::XAD4 = MaterialConstants( {2.5 * ns, 10 * ns, 90 * ns}, 
                                                                      {3.3, 2.8, 44.8}, 0.125 * ns, 1.0/3.0, {0.45 * ns}, {100} );
const MaterialConstants MaterialParameters::Al = MaterialConstants( {142 * ns}, {0}, 0.125 * ns, 0, {0.32 * ns}, {100} );
const MaterialConstants MaterialParameters::Kapton = MaterialConstants( {142 * ns}, {0}, 0.125 * ns, 0, {0.374 * ns}, {100} );
const MaterialConstants MaterialParameters::Plexiglass = MaterialConstants( {1.8 * ns}, {28}, 0.125 * ns, 1.0/3.0, {0.4 * ns}, {100} );
// DOI 10.1063/1.475876 for temperature 293 K
const MaterialConstants MaterialParameters::Scin = MaterialConstants( {1.95 * ns}, {35}, 0.125 * ns, 1.0/3.0, {0.45 * ns}, {100} );
// DOI 10.1515/nuka-2015-0140 for temperature 293 K
const MaterialConstants MaterialParameters::PA6 = MaterialConstants( {1.7 * ns}, {19.5}, 0.125 * ns, 1.0/3.0, {0.356 * ns}, {100} );
// DOI 10.1002/app.10319 for temperature 293 K
const MaterialConstants MaterialParameters::Air = MaterialConstants( {142 * ns}, {75}, 0.125 * ns, 1.0/3.0, {0.5 * ns}, {0} );

MaterialParameters::MaterialParameters()
{
  foPsLifetimes.clear();
  foPsProbabilities.clear();
  foPs2GIntensities.clear();
  foPs3GIntensities.clear();
  fpPsLifetime = 0;
  fpPsFraction = 0;
  fpPsIntensity = 0;
  fdirectLifetimes.clear();
  fdirectIntensities.clear();
  AnnihlationMode = "";
}

MaterialParameters::MaterialParameters(const std::vector<G4double> & oPsLF, const std::vector<G4double> & oPsProb, 
                                        G4double pPsLF, G4double pPsFrac, 
                                        const std::vector<G4double> & DirectLF, const std::vector<G4double> & DirectProb)
{
  foPsLifetimes = oPsLF;
  foPsProbabilities = oPsProb; 
  fpPsLifetime = pPsLF;
  fpPsFraction = pPsFrac;
  fdirectLifetimes = DirectLF;
  fdirectProbabilities = DirectProb;
  SetComponentsIntensities();
  AnnihlationMode = "";
}

void MaterialParameters::ClearTemp()
{
  Temp.oPsLifetimes.clear();
  Temp.oPsProbabilities.clear();
  Temp.pPsLifetime = 0.;
  Temp.pPsFraction = 0.;
  Temp.directLifetimes.clear();
  Temp.directProbabilities.clear();
}

void MaterialParameters::SetMaterial(MaterialConstants mat)
{
  foPsLifetimes = mat.oPsLifetimes;
  foPsProbabilities = mat.oPsProbabilities; 
  fpPsLifetime = mat.pPsLifetime;
  fpPsFraction = mat.pPsFraction;
  fdirectLifetimes = mat.directLifetimes;
  fdirectProbabilities = mat.directProbabilities;
  SetComponentsIntensities();
}

void MaterialParameters::SetMaterialByName(MaterialID materialID)
{
  if(materialID == MaterialID::mXAD4)
    SetMaterial(XAD4);
  else if(materialID == MaterialID::mAl)
    SetMaterial(Al);
  else if(materialID == MaterialID::mKapton)
    SetMaterial(Kapton);
  else if(materialID == MaterialID::mPlexiglass)
    SetMaterial(Plexiglass);
  else if(materialID == MaterialID::mScin)
    SetMaterial(Scin);
  else if(materialID == MaterialID::mPA6)
    SetMaterial(PA6);
  else if(materialID == MaterialID::mAir)
    SetMaterial(Air);
  else
    SetMaterial(Temp);
}

void MaterialParameters::SetAnnihilationMode(G4String mode)
{
  AnnihlationMode = mode;
}

void MaterialParameters::AddoPsComponent(G4double lifetime, G4double probability)
{
  foPsLifetimes.push_back(lifetime);
  foPsProbabilities.push_back(probability);
}

void MaterialParameters::SetoPsComponents(const std::vector<G4double> & oPsLF, const std::vector<G4double> & oPsProb)
{
  foPsLifetimes = oPsLF;
  foPsProbabilities = oPsProb;    
}

void MaterialParameters::SetpPsComponent(G4double lifetime, G4double fraction)
{
  fpPsLifetime = lifetime;
  fpPsFraction = fraction;    
}

void MaterialParameters::AddDirectComponent(G4double lifetime, G4double probability)
{
  fdirectLifetimes.push_back(lifetime);
  fdirectProbabilities.push_back(probability);
}

void MaterialParameters::SetDirectComponents(const std::vector<G4double> & directLF, const std::vector<G4double> & directProb)
{
  fdirectLifetimes = directLF;
  fdirectProbabilities = directProb;
}

void MaterialParameters::SetComponentsIntensities()
{
  foPs2GIntensities.clear();
  foPs3GIntensities.clear();
  fdirectIntensities.clear();
  G4double oPsTotalProbability = 0., TotalIntensity = 0.;
  for(unsigned i=0; i<foPsProbabilities.size(); i++)
  {
    oPsTotalProbability += foPsProbabilities[i];
    foPs2GIntensities.push_back(GetoPsIntensity2G(foPsLifetimes[i], foPsProbabilities[i]));
    foPs3GIntensities.push_back(GetoPsIntensity3G(foPsLifetimes[i], foPsProbabilities[i]));
  }
  TotalIntensity += oPsTotalProbability/100.;
  fpPsIntensity = oPsTotalProbability*fpPsFraction/100; // Fraction from 0 to whatever;
  TotalIntensity += fpPsIntensity;
  G4double directTotalProbability = 0., temp = 0.;
  for(unsigned j=0; j<fdirectProbabilities.size(); j++)
  {
    directTotalProbability = fdirectProbabilities[j];
  }
  for(unsigned j=0; j<fdirectProbabilities.size(); j++)
  {
    temp = (100 - oPsTotalProbability - 100*fpPsIntensity)*fdirectProbabilities[j]/(100*directTotalProbability);
    if(temp > 0. && temp <= 1.)
    {
      fdirectIntensities.push_back(temp);
      TotalIntensity += temp;
    }
    else
      fdirectIntensities.push_back(0);
  }
  if(foPsProbabilities.size() == 0)
  {
    foPs2GIntensities.push_back(0.);
    foPs3GIntensities.push_back(0.);
  }
  if(fdirectProbabilities.size() == 0)
    fdirectIntensities.push_back(0.);
  if(TotalIntensity > 1)
    RenormalizeIntensities(TotalIntensity);
}

void MaterialParameters::RenormalizeIntensities(G4double TotalIntensity)
{
  for(unsigned i=0; i<foPs2GIntensities.size(); i++)
  {
    foPs2GIntensities[i] = foPs2GIntensities[i]/TotalIntensity;
    foPs3GIntensities[i] = foPs3GIntensities[i]/TotalIntensity;
  }
  fpPsIntensity = fpPsIntensity/TotalIntensity;
  for(unsigned i=0; i<fdirectIntensities.size(); i++)
  {
    fdirectIntensities[i] = fdirectIntensities[i]/TotalIntensity;
  }
}

G4double MaterialParameters::GetoPs2GTotalIntensity()
{
  G4double TotalIntensity = 0.;
  for(unsigned i=0; i<foPs2GIntensities.size(); i++)
    TotalIntensity += foPs2GIntensities[i];
  return TotalIntensity;
}

G4double MaterialParameters::GetoPs3GTotalIntensity()
{
  G4double TotalIntensity = 0.;
  for(unsigned i=0; i<foPs3GIntensities.size(); i++)
    TotalIntensity += foPs3GIntensities[i];
  return TotalIntensity;    
}

G4double MaterialParameters::GetDirect2GTotalIntensity()
{
  G4double TotalIntensity = 0.;
  for(unsigned i=0; i<fdirectIntensities.size(); i++)
    TotalIntensity += fdirectIntensities[i];
  return TotalIntensity*(1.-direct3Gfraction);    
}

G4double MaterialParameters::GetDirect3GTotalIntensity()
{
  G4double TotalIntensity = 0.;
  for(unsigned i=0; i<fdirectIntensities.size(); i++)
    TotalIntensity += fdirectIntensities[i];
  return TotalIntensity*direct3Gfraction;      
}

G4double MaterialParameters::GetpPs2GTotalIntensity()
{
  return fpPsIntensity;
}

G4double MaterialParameters::GetoPsIntensity2G(G4double oPsLifetime, G4double oPsProbability)
{
  return ( (oPsTauVaccum - oPsLifetime) / oPsTauVaccum) * oPsProbability / 100.;
}

G4double MaterialParameters::GetoPsIntensity3G(G4double oPsLifetime, G4double oPsProbability)
{
  return (oPsLifetime / oPsTauVaccum) * oPsProbability / 100.;
}

G4double MaterialParameters::GetoPs2GLifetimeFromVector(double randNumber)
{
  double randNumberNorm = randNumber;
  if(AnnihlationMode != "")
    randNumberNorm = randNumber * GetoPs2GTotalIntensity();
  for(unsigned i=0; i<foPsLifetimes.size(); i++)
  {
    if(randNumberNorm > foPs2GIntensities[i])
      return foPsLifetimes[i];
  }
  return GetLifetimeVector(foPsLifetimes);
}

G4double MaterialParameters::GetoPs3GLifetimeFromVector(double randNumber)
{
  double randNumberNorm = randNumber;
  if(AnnihlationMode != "")
    randNumberNorm = randNumber * GetoPs3GTotalIntensity();
  for(unsigned i=0; i<foPsLifetimes.size(); i++)
  {
    if(randNumberNorm > foPs3GIntensities[i])
        return foPsLifetimes[i];
  }
  return GetLifetimeVector(foPsLifetimes);
}

G4double MaterialParameters::GetDirectLifetimeFromVector(double randNumber)
{
  double randNumberNorm = randNumber;
  if(AnnihlationMode != "")
    randNumberNorm = randNumber * GetDirect2GTotalIntensity();
  for(unsigned i=0; i<fdirectLifetimes.size(); i++)
  {
    if(randNumberNorm > fdirectIntensities[i])
      return fdirectLifetimes[i];
  }
  return GetLifetimeVector(fdirectLifetimes);
}

G4double MaterialParameters::GetpPsLifetime()
{ 
  return fpPsLifetime; 
}

G4double MaterialParameters::GetLifetimeVector(std::vector<G4double> vectorToCheck)
{
  if(vectorToCheck.size())
    return vectorToCheck[vectorToCheck.size() - 1];
  else
    return oPsTauVaccum;
}
