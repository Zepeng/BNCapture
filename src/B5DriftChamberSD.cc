//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: B5DriftChamberSD.cc 101036 2016-11-04 09:00:23Z gcosmo $
//
/// \file B5DriftChamberSD.cc
/// \brief Implementation of the B5DriftChamber class

#include "B5DriftChamberSD.hh"
#include "B5DriftChamberHit.hh"

#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include "G4HadronicProcess.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B5DriftChamberSD::B5DriftChamberSD(G4String name)
: G4VSensitiveDetector(name), 
  fHitsCollection(nullptr), fHCID(-1)
{
  collectionName.insert("driftChamberColl");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B5DriftChamberSD::~B5DriftChamberSD()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B5DriftChamberSD::Initialize(G4HCofThisEvent* hce)
{
  fHitsCollection 
    = new B5DriftChamberHitsCollection(SensitiveDetectorName,collectionName[0]);

  if (fHCID<0) { 
     fHCID = G4SDManager::GetSDMpointer()->GetCollectionID(fHitsCollection); 
  }
  hce->AddHitsCollection(fHCID,fHitsCollection);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool B5DriftChamberSD::ProcessHits(G4Step* step, G4TouchableHistory*)
{
  auto charge = step->GetTrack()->GetDefinition()->GetPDGCharge();
  auto preStepPoint = step->GetPreStepPoint();
  //if (charge==0.) return true;
  auto touchable = step->GetPreStepPoint()->GetTouchable();
  auto motherPhysical = touchable->GetVolume(1); // mother
  auto copyNo = motherPhysical->GetCopyNo();

  auto worldPos = preStepPoint->GetPosition();
  auto localPos 
    = touchable->GetHistory()->GetTopTransform().TransformPoint(worldPos);
  
  auto hit = new B5DriftChamberHit(copyNo);
  hit->SetWorldPos(worldPos);
  hit->SetLocalPos(localPos);
  hit->SetTime(preStepPoint->GetGlobalTime());
  
  fHitsCollection->insert(hit);
  
  // count processes
  // 
  const G4StepPoint* endPoint = step->GetPostStepPoint();
  G4VProcess* process   = 
                   const_cast<G4VProcess*>(endPoint->GetProcessDefinedStep());
  
  // check that an real interaction occured (eg. not a transportation)
  G4StepStatus stepStatus = endPoint->GetStepStatus();
  G4bool transmit = (stepStatus==fGeomBoundary || stepStatus==fWorldBoundary);
  if (transmit) return true;
                      
  //real processes : sum track length
  //
  G4double stepLength = step->GetStepLength();
  
  //energy-momentum balance initialisation
  //
  const G4StepPoint* prePoint = step->GetPreStepPoint();
  G4double Q             = - prePoint->GetKineticEnergy();
  G4ThreeVector Pbalance = - prePoint->GetMomentum();
  
  //initialisation of the nuclear channel identification
  //
  G4ParticleDefinition* particle = step->GetTrack()->GetDefinition();
  G4String partName = particle->GetParticleName();
  G4String nuclearChannel = partName;
  G4HadronicProcess* hproc = dynamic_cast<G4HadronicProcess*>(process);
  const G4Isotope* target = NULL;
  if (hproc) target = hproc->GetTargetIsotope();
  G4String targetName = "XXXX";  
  if (target) targetName = target->GetName();
  nuclearChannel += " + " + targetName + " --> ";
    
  //scattered primary particle (if any)
  //
  G4int ih = 1;
  if (step->GetTrack()->GetTrackStatus() == fAlive) {
    G4double energy = endPoint->GetKineticEnergy();      
    //
    G4ThreeVector momentum = endPoint->GetMomentum();
    Q        += energy;
    Pbalance += momentum;
    //
    nuclearChannel += partName + " + ";
    std::cout << "primaty scatter" << energy << std::endl;
  }  
  
  //secondaries
  //
  const std::vector<const G4Track*>* secondary 
                                    = step->GetSecondaryInCurrentStep();  
  for (size_t lp=0; lp<(*secondary).size(); lp++) {
    particle = (*secondary)[lp]->GetDefinition(); 
    G4String name   = particle->GetParticleName();
    G4String type   = particle->GetParticleType();      
    G4double energy = (*secondary)[lp]->GetKineticEnergy();
    //energy-momentum balance
    G4ThreeVector momentum = (*secondary)[lp]->GetMomentum();
    Q        += energy;
    Pbalance += momentum;
    //count e- from internal conversion together with gamma
    if (particle == G4Electron::Electron()) particle = G4Gamma::Gamma();
    std::cout << "secondary" << energy << std::endl;
    //particle flag
    //fParticleFlag[particle]++;
  }
  
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
