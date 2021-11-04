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

#include "DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4SubtractionSolid.hh"

#include "globals.hh"

#include "math.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
: G4VUserDetectorConstruction(),
  fScoringVolume(0)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct(){
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();

  //bool overlapsChecking = false;
  bool overlapsChecking = true;
  //bool buildSingleSiPMArray = true;
  bool buildSingleSiPMArray = false;

  G4double projection_surface_R = 208.0*mm;
  
  //     
  // World
  //
  G4double world_sizeXY = 100*cm;
  G4double world_sizeZ  = 200*cm;
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");
  
  G4Box* solidWorld = new G4Box("World", 0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ);
  G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld,world_mat,"World");                                   
  G4VPhysicalVolume* physWorld = new G4PVPlacement(0,                     //no rotation
						   G4ThreeVector(),       //at (0,0,0)
						   logicWorld,            //its logical volume
						   "World",               //its name
						   0,                     //its mother  volume
						   false,                 //no boolean operation
						   0,                     //copy number
						   overlapsChecking);     //overlaps checking

  G4double siPMEnvelope1_sizeX = 270*mm;
  G4double siPMEnvelope1_sizeY = 140*mm;
  G4double siPMEnvelope1_sizeZ = (2*projection_surface_R + 20*mm);  
  G4Box* solidSiPMEnvelope1 = new G4Box("solidSiPMEnvelope1", siPMEnvelope1_sizeX/2.0, siPMEnvelope1_sizeY/2.0, siPMEnvelope1_sizeZ/2.0);
  G4LogicalVolume* logicSiPMEnvelope1 = new G4LogicalVolume(solidSiPMEnvelope1,world_mat,"logicSiPMEnvelope1");                                   
  /*
  new G4PVPlacement(0,                     //rotation
		    G4ThreeVector(),       //position
		    logicSiPMEnvelope1,     //its logical volume
		    "logicSiPMEnvelope1",               //its name
		    logicWorld,            //its mother  volume
		    false,                 //no boolean operation
		    0,                     //copy number
		    overlapsChecking);     //overlaps checking
  */
  
  G4double siPMEnvelope2_sizeX = siPMEnvelope1_sizeX + 10*mm;
  G4double siPMEnvelope2_sizeY = siPMEnvelope1_sizeY + 10*mm;
  G4double siPMEnvelope2_sizeZ = siPMEnvelope1_sizeZ - 50*mm;  
  G4Box* solidSiPMEnvelope2 = new G4Box("solidSiPMEnvelope2", siPMEnvelope2_sizeX/2.0, siPMEnvelope2_sizeY/2.0, siPMEnvelope2_sizeZ/2.0);
  G4LogicalVolume* logicSiPMEnvelope2 = new G4LogicalVolume(solidSiPMEnvelope2,world_mat,"logicSiPMEnvelope2");
  /*
  new G4PVPlacement(0,                     //rotation
		    G4ThreeVector(),       //position
		    logicSiPMEnvelope2,     //its logical volume
		    "logicSiPMEnvelope2",               //its name
		    logicWorld,            //its mother  volume
		    false,                 //no boolean operation
		    0,                     //copy number
		    overlapsChecking);     //overlaps checking
  */
  G4RotationMatrix* rotMatrix = new G4RotationMatrix();
  G4ThreeVector transVector(0.0,0.0,(siPMEnvelope1_sizeZ - siPMEnvelope2_sizeZ)/2.0+10*mm);
  G4SubtractionSolid *solid_SiPMEnvelope1_m_SiPMEnvelope2 = new G4SubtractionSolid("solid_SiPMEnvelope1_m_SiPMEnvelope2",solidSiPMEnvelope1,solidSiPMEnvelope2,rotMatrix,transVector);
  G4LogicalVolume* logic_SiPMEnvelope1_m_SiPMEnvelope2 = new G4LogicalVolume(solid_SiPMEnvelope1_m_SiPMEnvelope2,world_mat,"logicSiPMEnvelope2");
  new G4PVPlacement(0,                     //rotation
		    G4ThreeVector(0.0,0.0,100*mm),       //position
		    logic_SiPMEnvelope1_m_SiPMEnvelope2,     //its logical volume
		    "logic_SiPMEnvelope1_m_SiPMEnvelope2",               //its name
		    logicWorld,            //its mother  volume
		    false,                 //no boolean operation
		    0,                     //copy number
		    overlapsChecking);     //overlaps checking

  
  //
  // Small box for orientation 
  //
  G4VSolid *boxsmall_solid = new G4Box("boxsmall_solid", 1.0*cm, 3.0*cm, 10.0*cm);
  G4LogicalVolume *boxsmall_logical = new G4LogicalVolume(boxsmall_solid,world_mat,"boxsmall_solid");
  new G4PVPlacement(0,                      //no rotation
		    G4ThreeVector(90*cm/2.0,90*cm/2.0,35*cm),       //at (0,0,0)
		    boxsmall_logical,       //its logical volume
		    "World",                //its name
		    logicWorld,             //its mother  volume
		    false,                  //no boolean operation
		    0,                      //copy number
		    overlapsChecking);      //overlaps checking

  //
  // Small box incenter 
  //
  G4VSolid *boxsmall_centre_solid = new G4Box("boxsmall_centre_solid", 1.0*mm, 1.0*mm, 1.0*mm);
  G4LogicalVolume *boxsmall_centre_logical = new G4LogicalVolume(boxsmall_centre_solid,world_mat,"boxsmall_centre_solid");
  if(!buildSingleSiPMArray)
    new G4PVPlacement(0,                       //no rotation
		      G4ThreeVector(),         //at (0,0,0)
		      boxsmall_centre_logical, //its logical volume
		      "World",                 //its name
		      logicWorld,              //its mother  volume
		      false,                   //no boolean operation
		      0,                       //copy number
		      overlapsChecking);       //overlaps checking
  
  return physWorld;
}
