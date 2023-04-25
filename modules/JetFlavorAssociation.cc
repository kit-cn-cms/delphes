/*
 *  Delphes: a framework for fast simulation of a generic collider experiment
 *  Copyright (C) 2012-2014  Universite catholique de Louvain (UCL), Belgium
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY || FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/** \class JetFlavorAssociation
 *
 *  Find origin of jet && evaluate jet flavor
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "modules/JetFlavorAssociation.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesFormula.h"

#include "ExRootAnalysis/ExRootClassifier.h"
#include "ExRootAnalysis/ExRootFilter.h"
#include "ExRootAnalysis/ExRootResult.h"

#include "TDatabasePDG.h"
#include "TFormula.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TRandom3.h"
#include "TString.h"

#include <algorithm>
#include <iostream>
#include <sstream>
#include <stdexcept>

using namespace std;

//------------------------------------------------------------------------------

class PartonClassifier : public ExRootClassifier
{
public:
  PartonClassifier() {}
  Int_t GetCategory(TObject *object);
  Double_t fEtaMax, fPTMin;
};

//------------------------------------------------------------------------------
// https://cmssdt.cern.ch/SDT/lxr/source/PhysicsTools/JetMCAlgos/plugins/PartonSelector.cc

Int_t PartonClassifier::GetCategory(TObject *object)
{
  // select parton in the parton list

  Candidate *parton = static_cast<Candidate *>(object);
  const TLorentzVector &momentum = parton->Momentum;
  Int_t pdgCode;

  // inside the eta && momentum range (be a little bit larger that the tracking coverage
  if(momentum.Pt() <= fPTMin || TMath::Abs(momentum.Eta()) > fEtaMax) return -1;

  pdgCode = TMath::Abs(parton->PID);

  if(parton->Status == -1) return -1;
  if(pdgCode != 21 && pdgCode > 5) return -1; // not a parton, skip
  if(parton->Status == 3 || parton->Status == 2) return 0; // if status 3 return

  return 0;
}

//------------------------------------------------------------------------------

class ParticleLHEFClassifier : public ExRootClassifier
{
public:
  ParticleLHEFClassifier() {}
  Int_t GetCategory(TObject *object);
  Double_t fEtaMax, fPTMin;
};

Int_t ParticleLHEFClassifier::GetCategory(TObject *object)
{
  // select parton in the parton list

  Candidate *particleLHEF = static_cast<Candidate *>(object);
  const TLorentzVector &momentum = particleLHEF->Momentum;
  Int_t pdgCode;

  // inside the eta && momentum range (be a little bit larger that the tracking coverage
  if(momentum.Pt() <= fPTMin || TMath::Abs(momentum.Eta()) > fEtaMax) return -1;

  pdgCode = TMath::Abs(particleLHEF->PID);
  if(particleLHEF->Status == -1) return -1;
  if(pdgCode != 21 && pdgCode > 5) return -1; // not a parton, skip
  if(particleLHEF->Status != 1) return -1; // if status 3 return

  return 0;
}

//------------------------------------------------------------------------------

JetFlavorAssociation::JetFlavorAssociation() :
  fPartonClassifier(0), fPartonFilter(0), fParticleLHEFFilter(0),
  fItPartonInputArray(0), fItParticleInputArray(0),
  fItParticleLHEFInputArray(0), fItJetInputArray(0)
{
  fPartonClassifier = new PartonClassifier;
  fParticleLHEFClassifier = new ParticleLHEFClassifier;
}

//------------------------------------------------------------------------------

JetFlavorAssociation::~JetFlavorAssociation()
{
  if(fPartonClassifier) delete fPartonClassifier;
  if(fParticleLHEFClassifier) delete fParticleLHEFClassifier;
}

//------------------------------------------------------------------------------

void JetFlavorAssociation::Init()
{
  ExRootConfParam param;

  fDeltaR = GetDouble("DeltaR", 0.5);

  fPartonClassifier->fPTMin = GetDouble("PartonPTMin", 0.0);
  fPartonClassifier->fEtaMax = GetDouble("PartonEtaMax", 2.5);

  fParticleLHEFClassifier->fPTMin = GetDouble("PartonPTMin", 0.0);
  fParticleLHEFClassifier->fEtaMax = GetDouble("PartonEtaMax", 2.5);

  // import input array(s)
  fPartonInputArray = ImportArray(GetString("PartonInputArray", "Delphes/partons"));
  fItPartonInputArray = fPartonInputArray->MakeIterator();
  fPartonFilter = new ExRootFilter(fPartonInputArray);

  fParticleInputArray = ImportArray(GetString("ParticleInputArray", "Delphes/allParticles"));
  fItParticleInputArray = fParticleInputArray->MakeIterator();

  try
  {
    fParticleLHEFInputArray = ImportArray(GetString("ParticleLHEFInputArray", "Delphes/allParticlesLHEF"));
  }
  catch(runtime_error &e)
  {
    fParticleLHEFInputArray = 0;
  }

  if(fParticleLHEFInputArray)
  {
    fItParticleLHEFInputArray = fParticleLHEFInputArray->MakeIterator();
    fParticleLHEFFilter = new ExRootFilter(fParticleLHEFInputArray);
  }

  fJetInputArray = ImportArray(GetString("JetInputArray", "FastJetFinder/jets"));
  fItJetInputArray = fJetInputArray->MakeIterator();
}

//------------------------------------------------------------------------------

void JetFlavorAssociation::Finish()
{
  if(fPartonFilter) delete fPartonFilter;
  if(fParticleLHEFFilter) delete fParticleLHEFFilter;

  if(fItJetInputArray) delete fItJetInputArray;
  if(fItParticleLHEFInputArray) delete fItParticleLHEFInputArray;
  if(fItParticleInputArray) delete fItParticleInputArray;
  if(fItPartonInputArray) delete fItPartonInputArray;
}

//------------------------------------------------------------------------------

void JetFlavorAssociation::Process()
{

  Candidate *jet;
  TObjArray *partonArray = 0;
  TObjArray *partonLHEFArray = 0;

  // select quark and gluons
  fPartonFilter->Reset();
  partonArray = fPartonFilter->GetSubArray(fPartonClassifier, 0); // get the filtered parton array
  if(partonArray == 0) return;

  if(fParticleLHEFInputArray)
  {
    fParticleLHEFFilter->Reset();
    partonLHEFArray = fParticleLHEFFilter->GetSubArray(fParticleLHEFClassifier, 0); // get the filtered parton array
  }
  // loop over all input jets
  fItJetInputArray->Reset();
  while((jet = static_cast<Candidate *>(fItJetInputArray->Next())))
  {
    // get standard flavor
    GetAlgoFlavor(jet, partonArray, partonLHEFArray);
    if(fParticleLHEFInputArray) GetPhysicsFlavor(jet, partonArray, partonLHEFArray);
  }
}

//------------------------------------------------------------------------------
// Standard definition of jet flavor in
// https://cmssdt.cern.ch/SDT/lxr/source/PhysicsTools/JetMCAlgos/plugins/JetPartonMatcher.cc?v=CMSSW_7_3_0_pre1

void JetFlavorAssociation::GetAlgoFlavor(Candidate *jet, TObjArray *partonArray, TObjArray *partonLHEFArray)
{
  float maxPt = 0;
  int daughterCounter = 0;
  Candidate *parton, *partonLHEF;
  Candidate *tempParton = 0, *tempPartonHighestPt = 0;
  int pdgCode, pdgCodeMax = -1;

  TIter itPartonArray(partonArray);
  TIter itPartonLHEFArray(partonLHEFArray);

  // Pointer to the matched parton
  Candidate* parton_matched = 0;
  // Keep track of the smallest DeltaR
  float deltaRMin = fDeltaR+0.1;
  // Origin of the jet flavor i.e. the relevant mother particle of the matched parton
  int flavorOrigin = 0;

  itPartonArray.Reset();
  while((parton = static_cast<Candidate *>(itPartonArray.Next())))
  {
    // default delphes method
    pdgCode = TMath::Abs(parton->PID);
    if(TMath::Abs(parton->PID) == 21) pdgCode = 0;
    if(jet->Momentum.DeltaR(parton->Momentum) <= fDeltaR)
    {
      if(pdgCodeMax == pdgCode) {
        //pdgCodeMax = pdgCode;
        if(jet->Momentum.DeltaR(parton->Momentum)<deltaRMin)
        {
          parton_matched = parton;
          deltaRMin = jet->Momentum.DeltaR(parton->Momentum);
        }
      }
      else if(pdgCodeMax < pdgCode) {
        pdgCodeMax = pdgCode;
        parton_matched = parton;
      }
    }

    if(!fParticleLHEFInputArray) continue;

    itPartonLHEFArray.Reset();
    while((partonLHEF = static_cast<Candidate *>(itPartonLHEFArray.Next())))
    {
      if(parton->Momentum.DeltaR(partonLHEF->Momentum) < 0.001 && parton->PID == partonLHEF->PID && partonLHEF->Charge == parton->Charge)
      {
        break;
      }

      // check the daughter
      daughterCounter = 0;
      if(parton->D1 != -1 || parton->D2 != -1)
      {
        // partons are only quarks || gluons
        int daughterFlavor1 = -1;
        int daughterFlavor2 = -1;
        if(parton->D1 != -1) daughterFlavor1 = TMath::Abs(static_cast<Candidate *>(fParticleInputArray->At(parton->D1))->PID);
        if(parton->D2 != -1) daughterFlavor2 = TMath::Abs(static_cast<Candidate *>(fParticleInputArray->At(parton->D2))->PID);
        if((daughterFlavor1 == 1 || daughterFlavor1 == 2 || daughterFlavor1 == 3 || daughterFlavor1 == 4 || daughterFlavor1 == 5 || daughterFlavor1 == 21)) daughterCounter++;
        if((daughterFlavor2 == 1 || daughterFlavor2 == 2 || daughterFlavor2 == 3 || daughterFlavor2 == 4 || daughterFlavor2 == 5 || daughterFlavor2 == 21)) daughterCounter++;
      }
      if(daughterCounter > 0) continue;
      if(jet->Momentum.DeltaR(parton->Momentum) <= fDeltaR)
      {
        // if not yet found && pdgId is a c, take as c
        if(TMath::Abs(parton->PID) == 4) tempParton = parton;
        if(TMath::Abs(parton->PID) == 5) tempParton = parton;
        if(parton->Momentum.Pt() > maxPt)
        {
          maxPt = parton->Momentum.Pt();
          tempPartonHighestPt = parton;
        }
      }
    }
  }

  if(!tempParton) tempParton = tempPartonHighestPt;
  jet->FlavorAlgo = tempParton ? TMath::Abs(tempParton->PID) : 0;

  if(pdgCodeMax == 0) pdgCodeMax = 21;
  if(pdgCodeMax == -1) pdgCodeMax = 0;

  // If we have a matched parton, let's get some more information
  if(parton_matched)
  {
    int parton_matched_pid = TMath::Abs(parton_matched->PID);
    // std::cout << "Parton matched " << parton_matched_pid << std::endl;
    // Check if the matched parton is a sensible parton in the sense of flavor matching meaning udscb quark or gluon
    if((parton_matched_pid>=1 && parton_matched_pid<=5) || (parton_matched_pid==21))
    {
      // First check whether the found parton is actually the last in the history or not
      // If not, go down the history and check for further similar partons
      Candidate* parton_daughter = parton_matched;
      while(parton_daughter && TMath::Abs(parton_daughter->PID)==parton_matched_pid)
      {
        // std::cout << "Going down in history" << std::endl;
        if(parton_daughter->D1>=0 && TMath::Abs(static_cast<Candidate *>(fPartonInputArray->At(parton_daughter->D1))->PID)==parton_matched_pid)
        {
          parton_daughter = static_cast<Candidate *>(fPartonInputArray->At(parton_daughter->D1));
        }
        else if(parton_daughter->D2>=0 && TMath::Abs(static_cast<Candidate *>(fPartonInputArray->At(parton_daughter->D2))->PID)==parton_matched_pid)
        {
          parton_daughter = static_cast<Candidate *>(fPartonInputArray->At(parton_daughter->D2));
        }
        else
        {
          break;
        }
      }
      // Set the matched parton to the last in the history
      parton_matched = parton_daughter;
      // Now let's go up the parton history until we find something other than the parton
      // While going up the history, we remove all the intermediate states from the considered parton list
      Candidate* parton_ancestor = parton_matched;
      std::vector<int> previous_mother_indices;
      while(parton_ancestor->M1>=0 && (parton_ancestor->M1<fPartonInputArray->GetEntries()) && std::find(previous_mother_indices.begin(),previous_mother_indices.end(),parton_ancestor->M1)==previous_mother_indices.end() && TMath::Abs(parton_ancestor->PID)==parton_matched_pid)
      {
        //std::cout << "Going up in history" << std::endl;
        previous_mother_indices.push_back(parton_ancestor->M1);
        parton_ancestor = static_cast<Candidate *>(fPartonInputArray->At(parton_ancestor->M1));
        partonArray->Remove(parton_ancestor);
      }
      // Now we are at the parton right before the parton chain
      // Check which particle this is and set a corresponding number
      // std::cout << "Final ancestor " << TMath::Abs(parton_ancestor->PID) << std::endl;
      if(TMath::Abs(parton_ancestor->PID)==6)
      {
        // found top quark
        flavorOrigin = 6;
        // search for the W from the top decay
        Candidate* w_candidate = 0;
        if(parton_ancestor->D1>=0 && TMath::Abs(static_cast<Candidate *>(fPartonInputArray->At(parton_ancestor->D1))->PID)==24)
        {
          w_candidate = static_cast<Candidate *>(fPartonInputArray->At(parton_ancestor->D1));
        }
        else if(parton_ancestor->D2>=0 && TMath::Abs(static_cast<Candidate *>(fPartonInputArray->At(parton_ancestor->D2))->PID)==24)
        {
          w_candidate = static_cast<Candidate *>(fPartonInputArray->At(parton_ancestor->D2));
        }
        // std::cout << "W candidate " << w_candidate << std::endl;
        if(w_candidate)
        {
          // now find the last W in the history
          while(w_candidate && TMath::Abs(w_candidate->PID)==24)
          {
            //std::cout << "Going down in history" << std::endl;
            if(w_candidate->D1>=0 && TMath::Abs(static_cast<Candidate *>(fPartonInputArray->At(w_candidate->D1))->PID)==24)
            {
              w_candidate = static_cast<Candidate *>(fPartonInputArray->At(w_candidate->D1));
            }
            else if(w_candidate->D2>=0 && TMath::Abs(static_cast<Candidate *>(fPartonInputArray->At(w_candidate->D2))->PID)==24)
            {
              w_candidate = static_cast<Candidate *>(fPartonInputArray->At(w_candidate->D2));
            }
            else
            {
              break;
            }
          }
          // after finding the last W, check whether it decays hadronically or leptonically
          if(w_candidate->D1>=0 && TMath::Abs(static_cast<Candidate *>(fPartonInputArray->At(w_candidate->D1))->PID)<=6 && w_candidate->D2>=0 && TMath::Abs(static_cast<Candidate *>(fPartonInputArray->At(w_candidate->D2))->PID)<=6)
          {
            // hadronic decay
            flavorOrigin = 60;
          }
          else if(w_candidate->D1>=0 && TMath::Abs(static_cast<Candidate *>(fPartonInputArray->At(w_candidate->D1))->PID)<=18 && w_candidate->D2>=0 && TMath::Abs(static_cast<Candidate *>(fPartonInputArray->At(w_candidate->D2))->PID)<=18)
          {
            // if not hadronic decay, then leptonic decay
            flavorOrigin = 61;
          }
        }
      }
      else if(TMath::Abs(parton_ancestor->PID)==23)
      {
        // found z boson
        flavorOrigin = 23;
      }
      else if(TMath::Abs(parton_ancestor->PID)==25)
      {
        // found higgs boson
        flavorOrigin = 25;
      }
      else if(TMath::Abs(parton_ancestor->PID)==21)
      {
        // found gluon
        flavorOrigin = 21;
      }
      else if(TMath::Abs(parton_ancestor->PID)==22)
      {
        // found photon
        flavorOrigin = 22;
      }
      else if(TMath::Abs(parton_ancestor->PID)==24)
      {
        // found w boson
        flavorOrigin = 24;
      }
      else if(TMath::Abs(parton_ancestor->PID)==2212)
      {
        // found nothing, history went back to initial proton
        flavorOrigin = 2212;
      }
    }
    // Also remove the matched parton from the considered partons
    partonArray->Remove(parton_matched);
  }
  // std::cout << "Flavor " << pdgCodeMax << std::endl;
  // std::cout << "FlavorOrigin " << flavorOrigin << std::endl;
  jet->Flavor = pdgCodeMax;
  jet->FlavorOrigin = flavorOrigin;
}

//------------------------------------------------------------------------------

void JetFlavorAssociation::GetPhysicsFlavor(Candidate *jet, TObjArray *partonArray, TObjArray *partonLHEFArray)
{
  int partonCounter = 0;
  float biggerConeSize = 0.7;
  float dist;
  bool isGoodCandidate;
  int contaminatingFlavor = 0;
  int motherCounter = 0;
  Candidate *parton, *partonLHEF, *mother1, *mother2;
  Candidate *tempParton = 0;
  vector<Candidate *> contaminations;
  vector<Candidate *>::iterator itContaminations;

  TIter itPartonArray(partonArray);
  TIter itPartonLHEFArray(partonLHEFArray);

  contaminations.clear();

  itPartonLHEFArray.Reset();
  while((partonLHEF = static_cast<Candidate *>(itPartonLHEFArray.Next())))
  {
    dist = jet->Momentum.DeltaR(partonLHEF->Momentum); // take the DR

    if(partonLHEF->Status == 1 && dist <= fDeltaR)
    {
      tempParton = partonLHEF;
      partonCounter++;
    }
  }

  itPartonArray.Reset();
  itPartonLHEFArray.Reset();
  while((parton = static_cast<Candidate *>(itPartonArray.Next())))
  {
    dist = jet->Momentum.DeltaR(parton->Momentum); // take the DR
    isGoodCandidate = true;
    while((partonLHEF = static_cast<Candidate *>(itPartonLHEFArray.Next())))
    {
      if(parton->Momentum.DeltaR(partonLHEF->Momentum) < 0.01 && parton->PID == partonLHEF->PID && partonLHEF->Charge == parton->Charge)
      {
        isGoodCandidate = false;
        break;
      }
    }

    if(!isGoodCandidate) continue;

    if(parton->D1 != -1 || parton->D2 != -1)
    {
      if((TMath::Abs(parton->PID) < 4 || TMath::Abs(parton->PID) == 21)) continue;
      if(dist < biggerConeSize) contaminations.push_back(parton);
    }
  }

  if(partonCounter != 1)
  {
    jet->FlavorPhys = 0;
  }
  else if(contaminations.size() == 0)
  {
    jet->FlavorPhys = TMath::Abs(tempParton->PID);
  }
  else if(contaminations.size() > 0)
  {
    jet->FlavorPhys = TMath::Abs(tempParton->PID);

    for(itContaminations = contaminations.begin(); itContaminations != contaminations.end(); ++itContaminations)
    {
      parton = *itContaminations;
      contaminatingFlavor = TMath::Abs(parton->PID);
      motherCounter = 0;
      if(parton->M1 != -1) motherCounter++;
      if(parton->M2 != -1) motherCounter++;

      if(parton->M1 != -1)
      {
        mother1 = static_cast<Candidate *>(fParticleInputArray->At(parton->M1));
        if(mother1 && motherCounter > 0 && mother1->Momentum.DeltaR(tempParton->Momentum) < 0.001) continue;
      }
      if(parton->M2 != -1)
      {
        mother2 = static_cast<Candidate *>(fParticleInputArray->At(parton->M2));
        if(mother2 && motherCounter > 0 && mother2->Momentum.DeltaR(tempParton->Momentum) < 0.001) continue;
      }
      // mother is the initialParton --> OK
      if(TMath::Abs(tempParton->PID) == 4)
      {
        // keep association --> the initialParton is a c --> the contaminated parton is a c
        if(contaminatingFlavor == 4) continue;
        jet->FlavorPhys = 0; // all the other cases reject!
        break;
      }
    }
  }
}
