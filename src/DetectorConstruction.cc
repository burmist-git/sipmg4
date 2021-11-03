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

#include "globals.hh"

#include "math.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
: G4VUserDetectorConstruction(),
  fScoringVolume(0)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{  
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();

  //     
  // World
  //
  G4double world_sizeXY = 100*cm;
  G4double world_sizeZ  = 100*cm;
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
						   false);                //overlaps checking

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
		    false);                 //overlaps checking

  //
  // Small box incenter 
  //
  G4VSolid *boxsmall_centre_solid = new G4Box("boxsmall_centre_solid", 1.0*mm, 1.0*mm, 1.0*mm);
  G4LogicalVolume *boxsmall_centre_logical = new G4LogicalVolume(boxsmall_centre_solid,world_mat,"boxsmall_centre_solid");
  new G4PVPlacement(0,                       //no rotation
		    G4ThreeVector(),         //at (0,0,0)
		    boxsmall_centre_logical, //its logical volume
		    "World",                 //its name
		    logicWorld,              //its mother  volume
		    false,                   //no boolean operation
		    0,                       //copy number
		    false);                  //overlaps checking

  bool overlapsChecking = false;

  G4double R_eff = 65;
  
  G4double L = 14.05*mm;
  G4double L_eff = 12.0*mm;
  G4double R = 208.0*mm;
  G4double d = 0.5*mm;
  unsigned int nOfRings = 10;
  G4double sipm_sizeZ = 1.5*mm;
  G4double alpha_y_0 = -20.3/180.0*CLHEP::pi;

  G4double alpha = 2.0*acos(get_cos_alpha_par_2(L, R, d));
  G4double alpha_new = 0;
  G4double R_new = 0;

  G4double sipm_sizeX = L;
  G4double sipm_sizeY = L;

  G4double sipm_x = 0;
  G4double sipm_y = 0;
  G4double sipm_z = 0;
  
  G4VSolid *projection_sphere_sipm_solid = new G4Box("projection_sphere_sipm_solid", sipm_sizeX/2.0, sipm_sizeY/2.0, sipm_sizeZ/2.0);
  G4LogicalVolume *projection_sphere_sipm_logical = new G4LogicalVolume(projection_sphere_sipm_solid,world_mat,"projection_sphere_sipm_logical");

  G4RotationMatrix *Ra;
  G4ThreeVector Ta;

  unsigned int nn = 0;
  unsigned int ringID = 0;

  //pixel
  G4double sipm_pixel_sizeX = 1.45*mm;
  G4double sipm_pixel_sizeY = 1.45*mm;
  G4double sipm_pixel_sizeZ = 0.5*mm;
  G4double sipm_pixel_pitch = 1.5*mm;
  unsigned int nx_sipm_pixel = 8;
  unsigned int ny_sipm_pixel = 8;
  //
  G4double sipm_pixel_x0 = -L/2.0 + (L - L_eff)/2.0 + sipm_pixel_sizeX/2.0;
  G4double sipm_pixel_y0 = -L/2.0 + (L - L_eff)/2.0 + sipm_pixel_sizeY/2.0;
  G4double sipm_pixel_z0 = 0.0;
  //  
  G4VSolid *sipm_pixel_solid = new G4Box("Sensitive", sipm_pixel_sizeX/2.0, sipm_pixel_sizeY/2.0, sipm_pixel_sizeZ/2.0);
  G4LogicalVolume *sensitive_logical = new G4LogicalVolume(sipm_pixel_solid,world_mat,"Sensitive");
  //
  G4RotationMatrix *Ra_pixel = new G4RotationMatrix();
  for(unsigned i = 0;i<ny_sipm_pixel;i++){
    for(unsigned j = 0;j<nx_sipm_pixel;j++){
      G4ThreeVector Ta_pixel;
      Ta_pixel.setX(sipm_pixel_x0 + sipm_pixel_pitch*j);
      Ta_pixel.setY(sipm_pixel_y0 + sipm_pixel_pitch*i);
      Ta_pixel.setZ(sipm_pixel_z0);
      new G4PVPlacement(Ra_pixel,           //rotation
			Ta_pixel,           //translation
			sensitive_logical,  //its logical volume
			"Sensitive",        //its name
			projection_sphere_sipm_logical, //its mother  volume
			false,              //no boolean operation
			0,                  //copy number
			overlapsChecking);  //overlaps checking
    }
  }
  
  //
  for(unsigned int j = 0;j<nOfRings;j++){
    ringID = j;
    R_new = R*cos(alpha*ringID);
    alpha_new = 2.0*acos(get_cos_alpha_par_2(L, R_new, d));
    nn = (unsigned int)(2*CLHEP::pi/alpha_new);
    for(unsigned int i = 0;i<nn;i++){
      Ra = new G4RotationMatrix();
      Ta = G4ThreeVector();
      sipm_x = 0;
      sipm_y = R*sin(alpha*ringID);
      sipm_z = R_new;
      Ta.setX(sipm_x);
      Ta.setY(sipm_y);
      Ta.setZ(sipm_z);
      Ta.rotateY(alpha_new*i + alpha_y_0);
      Ra->rotateY(-alpha_new*i - alpha_y_0);
      Ra->rotateX(alpha*ringID);
      if(sqrt(Ta.x()*Ta.x() + Ta.y()*Ta.y())<R_eff  && Ta.z()<0.0 ){
	new G4PVPlacement(Ra,                 //rotation
			  Ta,                 //translation
			  projection_sphere_sipm_logical,    //its logical volume
			  "projection_sphere_sipm_physical", //its name
			  logicWorld,         //its mother  volume
			  false,              //no boolean operation
			  0,                  //copy number
			  overlapsChecking);  //overlaps checking
      }
    }
  }

  for(unsigned int j = 1;j<nOfRings;j++){
    ringID = j;
    R_new = R*cos(alpha*ringID);
    alpha_new = 2.0*acos(get_cos_alpha_par_2(L, R_new, d));
    nn = (unsigned int)(2*CLHEP::pi/alpha_new);
    for(unsigned int i = 0;i<nn;i++){
      Ra = new G4RotationMatrix();
      Ta = G4ThreeVector();
      sipm_x = 0;
      sipm_y = -R*sin(alpha*ringID);
      sipm_z = R_new;
      Ta.setX(sipm_x);
      Ta.setY(sipm_y);
      Ta.setZ(sipm_z);
      Ta.rotateY(alpha_new*i + alpha_y_0);
      Ra->rotateY(-alpha_new*i - alpha_y_0);
      Ra->rotateX(-alpha*ringID);
      if(sqrt(Ta.x()*Ta.x() + Ta.y()*Ta.y())<R_eff && Ta.z()<0.0){
	new G4PVPlacement(Ra,                 //rotation
			  Ta,                 //translation
			  projection_sphere_sipm_logical,    //its logical volume
			  "projection_sphere_sipm_physical", //its name
			  logicWorld,         //its mother  volume
			  false,              //no boolean operation
			  0,                  //copy number
			  overlapsChecking);  //overlaps checking
      }
    }
  }

  /*
  //Satellite visualisation attributes
  G4VisAttributes* satelliteVisAtt = new G4VisAttributes();
  G4Color SatelliteColor = G4Color(138.0/255.0,43.0/255.0,226.0/255.0);
  satelliteVisAtt->SetColor(SatelliteColor);
  //satelliteVisAtt->SetVisibility(true);
  satelliteVisAtt->SetVisibility(false);
  terzina_satellite_logical->SetVisAttributes(satelliteVisAtt);
  terzina_satellite_back_wall_logical->SetVisAttributes(satelliteVisAtt);
  terzina_satellite_front_wall_logical->SetVisAttributes(satelliteVisAtt);

  //Absorber visualisation attributes
  G4VisAttributes* absorberVisAtt = new G4VisAttributes();
  G4Color red = G4Color(240.0/255.0,128.0/255.0,128.0/255.0, 0.1);
  absorberVisAtt->SetColor(red);
  //absorberVisAtt->SetVisibility(true);
  absorberVisAtt->SetVisibility(false);
  small_mirror_absorber_logical->SetVisAttributes(absorberVisAtt);
  projection_sphere_absorber_logical->SetVisAttributes(absorberVisAtt);
  */
  
  //Sensitive volume
  G4VisAttributes* sensitiveVisAtt = new G4VisAttributes();
  G4Color darkOrange = G4Color(255/255.0,140/255.0,0/255.0);
  sensitiveVisAtt->SetColor(darkOrange);
  sensitiveVisAtt->SetVisibility(true);
  sensitive_logical->SetVisAttributes(sensitiveVisAtt);

  //SiPM array
  G4VisAttributes* sipmArrayVisAtt = new G4VisAttributes();
  G4Color navyBlue = G4Color(0.0,0.0,128.0/255.0);
  sipmArrayVisAtt->SetColor(navyBlue);
  sipmArrayVisAtt->SetVisibility(true);
  projection_sphere_sipm_logical->SetVisAttributes(sipmArrayVisAtt);
  
  /*  
  for(unsigned j = 0;j<1;j++){
    alpha_i = 2.0*acos(get_cos_alpha_par_2(L, R*cos(alpha*j), d));
    Ra = new G4RotationMatrix();
    Ta = G4ThreeVector();
    sipm_x = 0;
    //sipm_y = R*sin(alpha*j);
    sipm_y = 0;
    //sipm_z = R*cos(alpha*j);
    sipm_z = R;
    Ta.setX(sipm_x);
    Ta.setY(sipm_y);
    Ta.setZ(sipm_z);    
    nn = (int)(2.0*CLHEP::pi/alpha_i);
    //G4cout<<nn<<G4endl;
    for(unsigned i = 0;i<100;i++){
      //Ta.rotateX((alpha_x_0 + alpha*j));
      Ta.rotateY((alpha_y_0 + alpha*i));
      //Ra->rotateX((-alpha_x_0 - alpha*j));
      Ra->rotateY((-alpha_y_0 - alpha*i));
      new G4PVPlacement(Ra,                 //rotation
			Ta,                 //translation
			sensitive_logical,  //its logical volume
			"Sensitive",        //its name
			logicWorld,         //its mother  volume
			false,              //no boolean operation
			0,                  //copy number
			false);             //overlaps checking
    }
  }
  */
  
  
  
  //}
  
  /*
  //  
  G4double sipm_sizeX = 1.0*mm;
  G4double sipm_sizeY = 1.0*mm;
  G4double sipm_sizeZ = 1.0*mm;
  G4double delta      = 2*mm;
  G4double rr         = 45.0*cm;
  G4double dphi       = 10.0/180.0*CLHEP::pi;
  G4double phi0       = 90.0/180.0*CLHEP::pi;
  G4double theta0     =  0.0/180.0*CLHEP::pi;
  G4double dtheta     = 10.0/180.0*CLHEP::pi;
  G4int nSiPM_x = 5;
  G4int nSiPM_y = 5;
  //G4double delta_x0 = -(sipm_sizeX + delta_mm)*(nSiPM_x-1)/2.0;
  //G4double delta_y0 = -(sipm_sizeY + delta_mm)*(nSiPM_y-1)/2.0;
  G4VSolid *projection_sphere_sipm_solid = new G4Box("Sensitive", sipm_sizeX/2.0, sipm_sizeY/2.0, sipm_sizeZ/2.0);
  G4LogicalVolume *sensitive_logical = new G4LogicalVolume(projection_sphere_sipm_solid,world_mat,"Sensitive");
  
  G4double sipm_x = 0;
  G4double sipm_y = 0;
  G4double sipm_z = 0;

  G4double phiy = 45.0/180.0*CLHEP::pi;
  G4double phix = 45.0/180.0*CLHEP::pi;

  G4RotationMatrix *Ra;
  G4ThreeVector Ta;
  G4RotationMatrix Raxtr;
  G4RotationMatrix Raytr;
  
  Ra = new G4RotationMatrix();
  Ta = G4ThreeVector();
  Raxtr = G4RotationMatrix();
  Raytr = G4RotationMatrix();
  sipm_x = 0;
  sipm_y = 0;
  sipm_z = 20*cm;
  Ta.setX(sipm_x);
  Ta.setY(sipm_y);
  Ta.setZ(sipm_z);
  new G4PVPlacement(Ra,                 //rotation
		    Ta,                 //translation
		    sensitive_logical,  //its logical volume
		    "Sensitive",        //its name
		    logicWorld,         //its mother  volume
		    false,              //no boolean operation
		    0,                  //copy number
		    false);             //overlaps checking

  Ra = new G4RotationMatrix();
  Ta = G4ThreeVector();
  Raxtr = G4RotationMatrix();
  Raytr = G4RotationMatrix();
  sipm_x = 0;
  sipm_y = 0;
  sipm_z = 20*cm;
  Ta.setX(sipm_x);
  Ta.setY(sipm_y);
  Ta.setZ(sipm_z);
  new G4PVPlacement(Ra,                //rotation
		    Ta,                //translation
		    sensitive_logical, //its logical volume
		    "Sensitive",       //its name
		    logicWorld,        //its mother  volume
		    false,             //no boolean operation
		    0,                 //copy number
		    false);            //overlaps checking
  */
  
  /*
  for(unsigned i = 0;i<nSiPM_x;i++){
    for(unsigned j = 0;j<nSiPM_y;j++){
      sipm_x = rr*sin(dtheta*j)*cos(dphi*i);
      sipm_y = rr*sin(dtheta*j)*sin(dphi*i);
      sipm_z = rr*cos(dtheta*j);
      G4ThreeVector Ta;
      Ta.setX(sipm_x);
      Ta.setY(sipm_y);
      Ta.setZ(sipm_z);
      new G4PVPlacement(0,                 //rotation
			Ta,                //translation
			sensitive_logical, //its logical volume
			"Sensitive",       //its name
			logicWorld,        //its mother  volume
			false,             //no boolean operation
			0,                 //copy number
			false);            //overlaps checking
    }
  }
  */

  /*
  G4RotationMatrix *Ra;
  G4ThreeVector Ta;
  G4RotationMatrix Raxtr;
  G4RotationMatrix Raytr;
  
  Ra = new G4RotationMatrix();
  Ta = G4ThreeVector();
  Raxtr = G4RotationMatrix();
  Raytr = G4RotationMatrix();
  sipm_x = 0;
  sipm_y = 0;
  sipm_z = -20*cm;
  Ta.setX(sipm_x);
  Ta.setY(sipm_y);
  Ta.setZ(sipm_z);
  new G4PVPlacement(Ra,                  //rotation
		    Ta, //translation
		    sensitive_logical,   //its logical volume
		    "Sensitive",         //its name
		    logicWorld,          //its mother  volume
		    false,               //no boolean operation
		    0,                   //copy number
		    false);              //overlaps checking

  Ra = new G4RotationMatrix();
  Ta = G4ThreeVector();
  Raxtr = G4RotationMatrix();
  Raytr = G4RotationMatrix();
  sipm_x = 0;
  sipm_y = 0;
  sipm_z = -20*cm;
  Raytr.rotateY(phiy);
  Ta.setX(sipm_x);
  Ta.setY(sipm_y);
  Ta.setZ(sipm_z);
  Ra->rotateY(-phiy);
  new G4PVPlacement(Ra,                 //rotation
		    Ta.transform(Raytr), //translation
		    sensitive_logical,  //its logical volume
		    "Sensitive",        //its name
		    logicWorld,         //its mother  volume
		    false,              //no boolean operation
		    0,                  //copy number
		    false);             //overlaps checking
  
  Ra = new G4RotationMatrix();
  Ta = G4ThreeVector();
  Raxtr = G4RotationMatrix();
  Raytr = G4RotationMatrix();
  sipm_x = 0;
  sipm_y = 0;
  sipm_z = -20*cm;
  Raytr.rotateY(-phiy);
  Ta.setX(sipm_x);
  Ta.setY(sipm_y);
  Ta.setZ(sipm_z);
  Ra->rotateY(phiy);
  new G4PVPlacement(Ra,                 //rotation
		    Ta.transform(Raytr), //translation
		    sensitive_logical,  //its logical volume
		    "Sensitive",        //its name
		    logicWorld,         //its mother  volume
		    false,              //no boolean operation
		    0,                  //copy number
		    false);             //overlaps checking

  Ra = new G4RotationMatrix();
  Ta = G4ThreeVector();
  Raxtr = G4RotationMatrix();
  Raytr = G4RotationMatrix();
  sipm_x = 0;
  sipm_y = 0;
  sipm_z = -20*cm;
  Raytr.rotateY(phiy);
  Raxtr.rotateX(phix);
  Ta.setX(sipm_x);
  Ta.setY(sipm_y);
  Ta.setZ(sipm_z);
  Ra->rotateY(-phiy);
  new G4PVPlacement(Ra,                  //rotation
		    Ta.transform(Raytr).transform(Raxtr), //translation
		    sensitive_logical,  //its logical volume
		    "Sensitive",        //its name
		    logicWorld,         //its mother  volume
		    false,              //no boolean operation
		    0,                  //copy number
		    false);             //overlaps checking
  */
  
  /*
  Ra = new G4RotationMatrix();
  Ta = G4ThreeVector();
  Ratr = G4RotationMatrix();
  sipm_x = 0;
  sipm_y = 0;
  sipm_z = -20*cm;
  Ratr.rotateY(-phiy);
  Ta.setX(sipm_x);
  Ta.setY(sipm_y);
  Ta.setZ(sipm_z);
  Ra->rotateY(phiy);
  new G4PVPlacement(Ra,                 //rotation
		    Ta.transform(Ratr), //translation
		    sensitive_logical,  //its logical volume
		    "Sensitive",        //its name
		    logicWorld,         //its mother  volume
		    false,              //no boolean operation
		    0,                  //copy number
		    false);             //overlaps checking
  */

  
  /*

  G4RotationMatrix Ra;
  G4ThreeVector Ta;
  G4Transform3D Tr;

  sipm_x = 0;
  sipm_y = 0;
  sipm_z = -20*cm;
  Ta.setX(sipm_x);
  Ta.setY(sipm_y);
  Ta.setZ(sipm_z);
  Ra.rotateY(phiy);
  Tr = G4Transform3D(Ra, Ta.transform(Ra));
  new G4PVPlacement(Tr,                //Transformation
		    sensitive_logical, //its logical volume				 
		    "Sensitive",       //its name
		    logicWorld,        //its mother  volume
		    false,             //no boolean operation
		    0);	               //copy number

  Ra=G4RotationMatrix();
  Ta=G4ThreeVector();
  Tr=G4Transform3D();
  sipm_x = 0;
  sipm_y = 0;
  sipm_z = -20*cm;
  Ta.setX(sipm_x);
  Ta.setY(sipm_y);
  Ta.setZ(sipm_z);
  Ra.rotateY(-phiy);
  Tr = G4Transform3D(Ra, Ta.transform(Ra));
  new G4PVPlacement(Tr,                //Transformation
		    sensitive_logical, //its logical volume				 
		    "Sensitive",       //its name
		    logicWorld,        //its mother  volume
		    false,             //no boolean operation
		    0);	               //copy number

  Ra=G4RotationMatrix();
  Ta=G4ThreeVector();
  Tr=G4Transform3D();
  sipm_x = 0;
  sipm_y = 0;
  sipm_z = -20*cm;
  Ta.setX(sipm_x);
  Ta.setY(sipm_y);
  Ta.setZ(sipm_z);
  Ra.rotateY(phiy);
  Ra.rotateX(phix);
  Tr = G4Transform3D(Ra, Ta.transform(Ra));
  new G4PVPlacement(Tr,                //Transformation
		    sensitive_logical, //its logical volume				 
		    "Sensitive",       //its name
		    logicWorld,        //its mother  volume
		    false,             //no boolean operation
		    0);	               //copy number
  */
  
  
  /*

  G4RotationMatrix Ra;
  G4ThreeVector Ta;
  G4Transform3D Tr;

  sipm_x = 0;
  sipm_y = 0;
  sipm_z = -20*cm;
  Ta.setX(sipm_x);
  Ta.setY(sipm_y);
  Ta.setZ(sipm_z);
  Ra.rotateY(phi);
  Tr = G4Transform3D(Ra, Ta);
  new G4PVPlacement(Tr,                //Transformation
		    sensitive_logical, //its logical volume				 
		    "Sensitive",       //its name
		    logicWorld,        //its mother  volume
		    false,             //no boolean operation
		    0);	               //copy number
  Ra.rotateY(-phi);
  */

  
  

  return physWorld;
}

G4double DetectorConstruction::get_cos_alpha_par_2(G4double L, G4double R, G4double d){
  //
  G4double A;
  G4double B;
  G4double C;
  G4double D;
  //
  G4double cos1;
  //G4double cos2;
  //
  A = (1.0 + L*L/4.0/R/R);
  B = (d*L/2/R/R);
  C = -(1.0 - d*d/4.0/R/R);
  D = B*B - 4.0*A*C;
  //
  if(D >= 0){
    cos1 = (-B + sqrt(D))/2.0/A;
    //cos2 = (-B - sqrt(D))/2.0/A;  
  }
  else{
    cos1 = -999.0;
    //cos2 = -999.0;
  }
  //  
  //cout<<"cos1 = "<<cos1<<endl
  //    <<"cos2 = "<<cos2<<endl;
  //cout<<"alpha1 = "<<2.0*TMath::ACos(cos1)*180.0/TMath::Pi()<<endl
  //    <<"alpha2 = "<<2.0*TMath::ACos(cos2)*180.0/TMath::Pi()<<endl;
  //
  return cos1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
