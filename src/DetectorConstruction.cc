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
/// \file DetectorConstruction.cc
/// \brief Implementation of the B1::DetectorConstruction class

#include "DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4Ellipsoid.hh"//椭球
#include "G4Tubs.hh"//圆柱

#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

//布尔操作
#include "G4UnionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4SubtractionSolid.hh"

// #include "G4CutTubs.hh"

// #include "G4Para.hh"

// #include "G4Torus.hh"
// #include "G4Trap.hh"





namespace B1
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();

  // Envelope parameters
  //
  G4double env_sizeXY = 100*cm, env_sizeZ = 220*cm;
  G4Material* env_mat = nist->FindOrBuildMaterial("G4_AIR");
  G4Material* body_mat = nist->FindOrBuildMaterial("G4_A-150_TISSUE");
  

  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;

  //
  // World
  //
  G4double world_sizeXY = 1.2*env_sizeXY;
  G4double world_sizeZ  = 1.2*env_sizeZ;
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");

  G4Box* solidWorld =
    new G4Box("World",                       //its name
       0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ);     //its size

  G4LogicalVolume* logicWorld =
    new G4LogicalVolume(solidWorld,          //its solid
                        world_mat,           //its material
                        "World");            //its name

  G4VPhysicalVolume* physWorld =
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking

  //
  // Envelope
  //
  G4Box* solidEnv =
    new G4Box("Envelope",                    //its name
        0.5*env_sizeXY, 0.5*env_sizeXY, 0.5*env_sizeZ); //its size

  G4LogicalVolume* logicEnv =
    new G4LogicalVolume(solidEnv,            //its solid
                        env_mat,             //its material
                        "Envelope");         //its name

  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(),         //at (0,0,0)
                    logicEnv,                //its logical volume
                    "Envelope",              //its name
                    logicWorld,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking

  


  //tumor-----------------------------------------------------------------
  G4double TumorX = 40.*mm;
  G4double TumorY = 20.*mm;
  G4double TumorZ = 60.*mm;
  G4Element* Boron10 = new G4Element("Boron10","B10",1);
  G4Isotope* Boron0 = new G4Isotope("Boron0",5,10,10.01*g/mole);
  Boron10->AddIsotope(Boron0,1.0);
  G4Material* tumor_mat = new G4Material("tumor_mat",1.128*g/cm3,2);
  tumor_mat->AddElement(Boron10,0.0001);
  tumor_mat->AddMaterial(body_mat,0.9999);//质量占比

  // G4Material* materialG4_B = nist->FindOrBuildMaterial("G4_B");

  G4Ellipsoid* solidTumor = 
    new G4Ellipsoid("Tumor",
                    0.5*TumorX,    //x半轴
                    0.5*TumorY,    //y半轴
                    0.5*TumorZ);    //z半轴


  G4LogicalVolume* logicTumor =
      new G4LogicalVolume(solidTumor,            //its solid
                          tumor_mat,             //its material
                          "Tumor");         //its name

  new G4PVPlacement(0,                       //no rotation set 0
                    G4ThreeVector(0.,0.,0.), //at (0,0,0)以肿瘤中心为中心
                    logicTumor,               //its logical volume
                    "Tumor",              //its name
                    logicEnv,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking


  //chest--------------------------------------------------------------  
  G4double ChestX = 260.*mm;
  G4double ChestY = 120.*mm;
  G4double ChestZ = 500.*mm;

  G4double TumorDisx = 50.*mm;
  G4double TumorDisz = 250.*mm;

  G4Ellipsoid* ChestEllip = 
    new G4Ellipsoid("ChestEllip",
                    0.5*TumorX,    //x半轴
                    0.5*TumorY,    //y半轴
                    0.5*TumorZ);    //z半轴

  G4Box* ChestBox =    
      new G4Box("ChestBox",    //its name
                0.5*ChestX, 
                0.5*ChestY, 
                0.5*ChestZ); //its size

  G4RotationMatrix* ChestxRot = new G4RotationMatrix;  // Rotates X and Z axes only
  ChestxRot->rotateX(0);      // Rotates 0
  G4ThreeVector ChestzTrans(-0.5*ChestX+TumorDisx, 0, 0.5*ChestZ-TumorDisz);//平移

  G4SubtractionSolid* solidChest =
    new G4SubtractionSolid("Chest",       //name
                            ChestBox,     //solid a
                            ChestEllip,   //solid b
                            ChestxRot,    //rotation
                            ChestzTrans); //trans b相对于a  以a的中心为原点

  G4LogicalVolume* logicChest =
      new G4LogicalVolume(solidChest,            //its solid
                          body_mat,             //its material
                          "Chest");         //its name

  new G4PVPlacement(0,                       //no rotation set 0
            G4ThreeVector(0.5*ChestX-TumorDisx,0.*mm,0.*mm), //at (80,0,0)mm 以a的中心为整体的中心
            logicChest,               //its logical volume
            "Chest",              //its name
            logicEnv,              //its mother  volume
            false,                   //no boolean operation
            0,                       //copy number
            checkOverlaps);          //overlaps checking


  //neck-----------------------------------------------------------------
  G4double NeckR = 50.*mm;
  G4double NeckZ = 90.*mm;

  G4Tubs* solidNeck =
    new G4Tubs("Neck",0,NeckR, 0.5*NeckZ,0,twopi);   // r:     0 mm -> 50 mm

  G4LogicalVolume* logicNeck =
      new G4LogicalVolume(solidNeck,            //its solid
                          body_mat,             //its material
                          "Neck");         //its name

  new G4PVPlacement(0,                       //no rotation set 0
                    G4ThreeVector(0.5*ChestX-TumorDisx, 0., 0.5*(ChestZ+NeckZ)), //at (80,0,295)
                    logicNeck,               //its logical volume
                    "Neck",              //its name
                    logicEnv,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking


  //head-----------------------------------------------------------------
  G4double HeadR = 90.*mm;
  G4double HeadDisz = 74.833*mm;

  G4Ellipsoid* solidHead = 
    new G4Ellipsoid("Head",
                    HeadR,    //x半轴
                    HeadR,    //y半轴
                    HeadR,     //z半轴
                    -HeadDisz,
                    HeadR);    

  G4LogicalVolume* logicHead =
      new G4LogicalVolume(solidHead,            //its solid
                          body_mat,             //its material
                          "Head");         //its name

  new G4PVPlacement(0,                       //no rotation set 0
                    G4ThreeVector(0.5*ChestX-TumorDisx, 0., 0.5*ChestZ+NeckZ+HeadDisz), //at (80,0,414.833)
                    logicHead,               //its logical volume
                    "Head",              //its name
                    logicEnv,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking



  //left_leg-------------------------------------------------------------- 
  G4double LegR = 55.*mm;
  G4double LegZ = 820.*mm;
  G4Tubs* solidLleg =
    new G4Tubs("Lleg",
                0,//内半径
                LegR,//外半径
                0.5*LegZ,//Z轴方向的半长度
                0*degree,//圆周起始位置弧度值
                360*degree);//该实体的圆心角弧度值

  G4LogicalVolume* logicLleg =
    new G4LogicalVolume(solidLleg,            //its solid
                        body_mat,             //its material
                        "Lleg");         //its name

  new G4PVPlacement(0,                       //no rotation set 0
                    G4ThreeVector(-TumorDisx+LegR, 0., -0.5*(ChestZ+LegZ)), //at (5,0,-660)
                    logicLleg,               //its logical volume
                    "Lleg",              //its name
                    logicEnv,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking


  //right_leg-------------------------------------------------------------- 
  G4Tubs* solidRleg =
    new G4Tubs("Rleg",
                0,//内半径
                LegR,//外半径
                0.5*LegZ,//Z轴方向的半长度
                0*degree,//圆周起始位置弧度值
                360*degree);//该实体的圆心角弧度值

  G4LogicalVolume* logicRleg =
    new G4LogicalVolume(solidRleg,            //its solid
                        body_mat,             //its material
                        "Rleg");         //its name

  new G4PVPlacement(0,                       //no rotation set 0
                    G4ThreeVector(ChestX-TumorDisx-LegR, 0., -0.5*(ChestZ+LegZ)), //at (155,0,-660)
                    logicRleg,               //its logical volume
                    "Rleg",              //its name
                    logicEnv,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking


  // Set Tumor as scoring volume
  //
  fScoringVolume = logicTumor;

  //
  // print parameters
  //
  G4cout
    << G4endl
    << "------------------------------------------------------------" << G4endl
    << "tumor_mat=" << tumor_mat->GetName() << G4endl
    << "------------------------------------------------------------" << G4endl;
  G4cout<<*(G4Material::GetMaterialTable())<<G4endl;
  //
  //always return the physical World
  //
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
