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
// $Id: B5DetectorConstruction.cc 103284 2017-03-24 08:27:19Z gcosmo $
//
/// \file B5DetectorConstruction.cc
/// \brief Implementation of the B5DetectorConstruction class

#include "B5DetectorConstruction.hh"
#include "B5Constants.hh"
#include "B5DriftChamberSD.hh"

#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4Mag_UsualEqRhs.hh"

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4MaterialTable.hh"
#include "G4NistManager.hh"

#include "G4VSolid.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4PVReplica.hh"
#include "G4UserLimits.hh"

#include "G4SDManager.hh"
#include "G4VSensitiveDetector.hh"
#include "G4RunManager.hh"
#include "G4GenericMessenger.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4ios.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal G4FieldManager* B5DetectorConstruction::fFieldMgr = 0;
    
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B5DetectorConstruction::B5DetectorConstruction()
: G4VUserDetectorConstruction(), 
  fMessenger(nullptr),
  //fWirePlane1Logical(nullptr), 
  fVisAttributes()
{
  // define commands for this class
  DefineCommands();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B5DetectorConstruction::~B5DetectorConstruction()
{
  delete fMessenger;
  
  for (auto visAttributes: fVisAttributes) {
    delete visAttributes;
  }  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B5DetectorConstruction::Construct()
{
  // Construct materials
  ConstructMaterials();
  auto air = G4Material::GetMaterial("G4_AIR");
  
  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;

  // geometries --------------------------------------------------------------
  // experimental hall (world volume)
  auto worldSolid 
    = new G4Box("worldBox",1.*m,1.*m,1.*m);
  auto worldLogical
    = new G4LogicalVolume(worldSolid,air,"worldLogical");
  auto worldPhysical
    = new G4PVPlacement(0,G4ThreeVector(),worldLogical,"worldPhysical",0,
                        false,0,checkOverlaps);
  
  // visualization attributes ------------------------------------------------
  // Steel Tube
  auto nist = G4NistManager::Instance();
  G4Material* shape1_mat = nist->FindOrBuildMaterial("G4_Fe");
  G4ThreeVector pos1 = G4ThreeVector(0, 0*cm, -6*cm);
  G4double tuberadius = 3*cm;
  G4double tubethickness = 0.2*cm;
  G4double tubelength = 6*cm;
  G4VSolid* tubeSolid = new G4Tubs("SteelTube", tuberadius, tuberadius+tubethickness, tubelength, 0*deg, 360*deg);
  G4LogicalVolume* logicTube = new G4LogicalVolume(tubeSolid, shape1_mat, "SteelTube");
  new G4PVPlacement(0, pos1, logicTube, "SteelTube", worldLogical, false, 0, checkOverlaps);

  //PE cylinder
  G4Material* pe_mat = nist->FindOrBuildMaterial("G4_POLYETHYLENE");
  G4ThreeVector pos2 = G4ThreeVector(0, 0, -2*cm);
  G4double hcylinder = 2*cm;
  G4VSolid* cylinderSolid = new G4Tubs("PE", 0, tuberadius, hcylinder, 0*deg, 360*deg);
  G4LogicalVolume* logicCylinder = new G4LogicalVolume(cylinderSolid, pe_mat, "PE");
  new G4PVPlacement(0, pos2, logicCylinder, "PE", worldLogical, false, 0, checkOverlaps);

  //Shape 3
  G4Material* sio2 = nist->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
  G4Isotope* b11 = new G4Isotope("B11", 5, 11, 11*g/mole);
  G4Isotope* b10 = new G4Isotope("B10", 5, 10, 11*g/mole);
  G4Element* boron = new G4Element("Boron", "Boron", 2);
  boron->AddIsotope(b10, 70*perCent);
  boron->AddIsotope(b11, 30*perCent);
  G4Material* boron_mat = new G4Material("boron_mat", 2.37*g/cm3, 1);
  boron_mat->AddElement(boron, 1);
  G4Material* sphere_mat = new G4Material("sphere_mat", 2.6*g/cm3, 2);
  sphere_mat->AddMaterial(sio2, 90*perCent);
  sphere_mat->AddMaterial(boron_mat, 10*perCent);
  G4double radius = 0.5*cm;
  G4Sphere* solidShape3 =
    new G4Sphere("Shape3",
                  0, radius,
                  0*deg, 360*deg,
                  0*deg, 180*deg);
  logicShape3 =
    new G4LogicalVolume(solidShape3,
                        sphere_mat,
                        "boron_sphere");
  G4ThreeVector pos3 = G4ThreeVector(0, 0*cm, 1*cm);
  new G4PVPlacement(0,
                        pos3,
                        logicShape3,
                        "Shape3",
                        worldLogical,
                        false,
                        0,
                        checkOverlaps);
  //
  
  auto visAttributes = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  visAttributes->SetVisibility(false);
  worldLogical->SetVisAttributes(visAttributes);
  fVisAttributes.push_back(visAttributes);
  
  visAttributes = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  visAttributes->SetVisibility(false);
  fVisAttributes.push_back(visAttributes);
  
  // return the world physical volume ----------------------------------------
  
  return worldPhysical;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B5DetectorConstruction::ConstructSDandField()
{
  // sensitive detectors -----------------------------------------------------
  auto sdManager = G4SDManager::GetSDMpointer();
  G4String SDname;
  auto chamber1 = new B5DriftChamberSD(SDname="/chamber1");
  sdManager->AddNewDetector(chamber1);
  //fWirePlane1Logical->SetSensitiveDetector(chamber1);
  logicShape3->SetSensitiveDetector(chamber1);
}    

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B5DetectorConstruction::ConstructMaterials()
{
  auto nistManager = G4NistManager::Instance();

  // Air 
  nistManager->FindOrBuildMaterial("G4_AIR");
  
  // Argon gas
  nistManager->FindOrBuildMaterial("G4_Ar");

  // Scintillator
  // (PolyVinylToluene, C_9H_10)
  nistManager->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
  
  // CsI
  nistManager->FindOrBuildMaterial("G4_CESIUM_IODIDE");
  
  // Lead
  nistManager->FindOrBuildMaterial("G4_Pb");
  
  G4cout << G4endl << "The materials defined are : " << G4endl << G4endl;
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B5DetectorConstruction::DefineCommands()
{
  // Define /B5/detector command directory using generic messenger class
  fMessenger = new G4GenericMessenger(this, 
                                      "/B5/detector/", 
                                      "Detector control");

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
