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
//
/// \file SteppingAction.cc
/// \brief Implementation of the B1::SteppingAction class

#include "SteppingAction.hh"
#include "EventAction.hh"
#include "DetectorConstruction.hh"

#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"
#include "G4AnalysisManager.hh"

namespace B1
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(EventAction* eventAction)
: fEventAction(eventAction)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* step)
{
  if (!fScoringVolume) {
    const DetectorConstruction* detConstruction
      = static_cast<const DetectorConstruction*>
        (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    fScoringVolume = detConstruction->GetScoringVolume();
  }

  // get volume of the current step
  G4LogicalVolume* volume
    = step->GetPreStepPoint()->GetTouchableHandle()
      ->GetVolume()->GetLogicalVolume();
  G4String volumename = volume->GetName();

  // check if we are in scoring volume
  // if (volume != fScoringVolume) return;

  if(volume == fScoringVolume)//tumor
  {
    // collect energy deposited in this step
    G4double edepStep = step->GetTotalEnergyDeposit();
    G4double stepPostx = step->GetPostStepPoint()->GetPosition().x();
    G4double stepPosty = step->GetPostStepPoint()->GetPosition().y();
    G4double stepPostz = step->GetPostStepPoint()->GetPosition().z();
    G4double stepPrex = step->GetPreStepPoint()->GetPosition().x();
    G4double stepPrey = step->GetPreStepPoint()->GetPosition().y();
    G4double stepPrez = step->GetPreStepPoint()->GetPosition().z();
    G4double stepx = 0.5*(stepPostx+stepPrex);
    G4double stepy = 0.5*(stepPosty+stepPrey);
    G4double stepz = 0.5*(stepPostz+stepPrez);
    
    fEventAction->AddEdep(edepStep);

    auto analysisManager = G4AnalysisManager::Instance();
    analysisManager->FillNtupleIColumn(0, 1);
    analysisManager->FillNtupleDColumn(1, stepx);
    analysisManager->FillNtupleDColumn(2, stepy);
    analysisManager->FillNtupleDColumn(3, stepz);
    analysisManager->FillNtupleDColumn(4, edepStep);
    analysisManager->AddNtupleRow();   
  }
  else if( (volumename=="Chest") || (volumename=="Neck") || (volumename=="Head")
            || (volumename=="Lleg") || (volumename=="Rleg") )
  {
    // collect energy deposited in this step
    G4double edepStep = step->GetTotalEnergyDeposit();
    G4double stepPostx = step->GetPostStepPoint()->GetPosition().x();
    G4double stepPosty = step->GetPostStepPoint()->GetPosition().y();
    G4double stepPostz = step->GetPostStepPoint()->GetPosition().z();
    G4double stepPrex = step->GetPreStepPoint()->GetPosition().x();
    G4double stepPrey = step->GetPreStepPoint()->GetPosition().y();
    G4double stepPrez = step->GetPreStepPoint()->GetPosition().z();
    G4double stepx = 0.5*(stepPostx+stepPrex);
    G4double stepy = 0.5*(stepPosty+stepPrey);
    G4double stepz = 0.5*(stepPostz+stepPrez);
    
    fEventAction->AddEdep(edepStep);

    auto analysisManager = G4AnalysisManager::Instance();
    analysisManager->FillNtupleIColumn(0, 2);
    analysisManager->FillNtupleDColumn(1, stepx);
    analysisManager->FillNtupleDColumn(2, stepy);
    analysisManager->FillNtupleDColumn(3, stepz);
    analysisManager->FillNtupleDColumn(4, edepStep);
    analysisManager->AddNtupleRow();  
  }
  else
    return;
  

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
