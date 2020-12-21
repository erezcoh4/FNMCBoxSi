/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4RunManager.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

#include "G4RotationMatrix.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction( int _fdebug_ )
:G4VUserDetectorConstruction(),
fWorldMater(0), fPhysiWorld(0),
fDetectorMessenger(0),fdebug(_fdebug_){
    DefineMaterials();
    fDetectorMessenger = new DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
DetectorConstruction::~DetectorConstruction(){ delete fDetectorMessenger;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4VPhysicalVolume* DetectorConstruction::Construct(){
    return ConstructVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::DefineMaterials() {
    // build materials
    G4Element* N  = new G4Element("Nitrogen", "N", 7, 14.01*g/mole);
    G4Element* O  = new G4Element("Oxygen",   "O", 8, 16.00*g/mole);
    G4int ncomponents; G4double fractionmass;
    G4Material* Air20 = new G4Material("Air", 1.205*mg/cm3, ncomponents=2, kStateGas, 293.*kelvin, 1.*atmosphere);
    Air20->AddElement(N, fractionmass=0.7);
    Air20->AddElement(O, fractionmass=0.3);
    fWorldMater = Air20;
    
    G4NistManager* man = G4NistManager::Instance();
    // source holder material - plastic
    fPlasticMater = man->FindOrBuildMaterial("G4_POLYETHYLENE");
    
    // scintillation materials - do not delete
    // scintillation materials from [http://www.sixiangguo.net/code/geant4/AppDevelop/apas06.html]
    //    fDetectorMater = man->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
    fDetectorMater = man->FindOrBuildMaterial("G4_STILBENE");
    // we should implememt AMCRYS UPS-113 here?
    
    
    // SiPMs
    fSiPMMater = man->FindOrBuildMaterial("G4_Si");;
    
    
    // electronics: PCB / Si
    // glass fiber reinforced
    fElectronicsMaterial = man->FindOrBuildMaterial("G4_Si");;
    
    
    if (fdebug>0){
        std::cout << "detector material: " << fDetectorMater -> GetName() << std::endl;
        std::cout << "material table:" << std::endl;
        std::cout << *(G4Material::GetMaterialTable()) << std::endl;
    }
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4VPhysicalVolume* DetectorConstruction::ConstructVolumes(){
    G4bool DoWorld              = TRUE;
    G4bool DoSourceHolder       = TRUE;
    G4bool DoSourcePlaceHolder  = FALSE;
    G4bool DoScintillators      = TRUE;
    G4bool DoSiPMs              = FALSE;
    G4bool DoElectronicPlates   = FALSE;
    
    // Cleanup old geometry
    G4GeometryManager::GetInstance()->OpenGeometry();
    G4PhysicalVolumeStore::GetInstance()->Clean();
    G4LogicalVolumeStore::GetInstance()->Clean();
    G4SolidStore::GetInstance()->Clean();
    
    
    // World
    if (DoWorld){
        G4Box * sWorld = new G4Box("World",                                 //name
                                   0.5*fWorldLength, 0.5*fWorldLength, 0.5*fWorldLength); //dimensions
        
        logicWorld = new G4LogicalVolume(sWorld,                    //shape
                                         fWorldMater,                   //material
                                         "World");                      //name
        
        fPhysiWorld = new G4PVPlacement(0,                          //no rotation
                                        G4ThreeVector(),            //at (0,0,0)
                                        logicWorld,                 //logical volume
                                        "World",                    //name
                                        0,                          //mother volume
                                        false,                      //no boolean operation
                                        0);                         //copy number
    }
    
    // (plastic) source holder
    if (DoSourceHolder) {
        posSourceHolder = G4ThreeVector(0 , 0 ,0 );
        solidSourceHolder = new G4Box("solidSourceHolder",
                                      1.*cm/2., 2.*cm/2. , sourceHolder_dz/2.);
        logicSourceHolder = new G4LogicalVolume(solidSourceHolder,
                                                fPlasticMater ,
                                                "logicSourceHolder",
                                                0,0,0);
        SourceHolderColor = G4Colour(0.8, 0.6, 0.7);
        SourceHolderVisAttributes = new G4VisAttributes(SourceHolderColor);
        logicSourceHolder->SetVisAttributes(SourceHolderVisAttributes);
        new G4PVPlacement(0,                                    // no rotation
                          G4ThreeVector(0,0,0),                 // at (x,y,z)
                          logicSourceHolder,                    // its logical volume
                          "placementSourceHolder",              // its name
                          logicWorld,                           // its mother  volume
                          false,                                // no boolean operations
                          0);                                   // copy number
        // <--- end source holder
    }
    
    // place-holder for the source position
    if (DoSourcePlaceHolder) {
        sourcePlaceHolder = new G4Sphere("sourcePlaceHolder",
                                         0, 0.1 * mm,
                                         0.*deg,360.*deg,0.*deg,180.*deg);
        logicSourcePlaceHolder = new G4LogicalVolume(sourcePlaceHolder,
                                                     fWorldMater ,
                                                     "logicSourcePlaceHolder",
                                                     0,0,0);
        SourceRed = G4Colour(1.0, 0.1, 0.2);
        sourceVisAttributes = new G4VisAttributes(SourceRed);
        logicSourcePlaceHolder->SetVisAttributes(sourceVisAttributes);
        new G4PVPlacement(0,                                  // no rotation
                          G4ThreeVector(0,Source_y,0),        // at (x,y,z)
                          logicSourcePlaceHolder,             // its logical volume
                          "placementSourcePlaceHolder",       // its name
                          logicSourceHolder,                  // its mother  volume
                          false,                              // no boolean operations
                          0);                                 // copy number
        // <--- end place-holder for the source position
    }
    
    
    // Scintillators
    if (DoScintillators){
        ScintBlue = G4Colour(0.1, 0.2, 1.0);
        ScintVisAttributes = new G4VisAttributes(ScintBlue);
        
        for (int facetIdx=0; facetIdx < NFacets; facetIdx++){
            
            // determine z according to facet idx
            G4ThreeVector FacetCentroid = GetFacetCentroid(facetIdx);
            if (fdebug>1){
                std::cout << "building facet " << FacetName(facetIdx) << ", "
                << "centroid at ("
                << FacetCentroid.x() << "," << FacetCentroid.y() << ","  << FacetCentroid.z()
                << ")" << std::endl;
            }
            
            // roatate to facet axes frame
            // G4RotationMatrix * FacetRotation = GetFacetRotation(facetIdx);
            G4Rotate3D FacetRot3D = GetFacetRot3D(facetIdx);
            // facet - a mother volume for multiple scintillators and SiPMs, made of air
            G4Colour FacetGrey(0.6, 0.6, 0.6);
            FacetVisAttributes = new G4VisAttributes(FacetGrey);
            
            solidFacet[facetIdx] = new G4Box( FacetName(facetIdx) + " facet",
                                             FacetSide/2.,
                                             FacetSide/2.,
                                             FacetThickness/2.);
            logicFacet[facetIdx] = new G4LogicalVolume(solidFacet[facetIdx],
                                                       fWorldMater ,
                                                       "logic "+FacetName(facetIdx) + " facet",
                                                       0,0,0);
            
            logicFacet[facetIdx] -> SetVisAttributes(FacetVisAttributes);
            
            // facet translation, followed by rotation, followed by another translation inside facet
            G4Transform3D FacetTransform = G4Translate3D( FacetCentroid ) * FacetRot3D;
            new G4PVPlacement(FacetTransform,
                              logicFacet[facetIdx],    // its logical volume
                              FacetName(facetIdx),       // its name
                              logicWorld,      // its mother  volume
                              false,           // no boolean operations
                              0);              // copy number
            
            for (int cellIdx_i=0; cellIdx_i < NCells_i; cellIdx_i++){
                for (int cellIdx_j=0; cellIdx_j < NCells_j; cellIdx_j++){
                    
                    
                    
                    
                    // scintillator position
                    // ----------------------
                    // 1. translation of facet centre
                    // 2. rotate axes frame to facet
                    // 3. position scintillator in centroid + dl(i) x i + dl(j) x j in i and j directions
                    //                int dl_i = Get_dl_i_facet(cellIdx);
                    //                int dl_j = Get_dl_j_facet(cellIdx);
                    
                    G4ThreeVector posInFacet = GetScintillatorPositionInFacet(cellIdx_i,cellIdx_j);
                    
                    solidScintillator[facetIdx][cellIdx_i][cellIdx_j] = new G4Box("solid"+ScintillatorLabel(facetIdx,cellIdx_i,cellIdx_j),
                                                                                  ScintillatorSide/2.,
                                                                                  ScintillatorSide/2.,
                                                                                  ScintillatorThickness/2.);
                    logicScintillator[facetIdx][cellIdx_i][cellIdx_j] = new G4LogicalVolume(solidScintillator[facetIdx][cellIdx_i][cellIdx_j],
                                                            fDetectorMater ,
                                                            "logic"+ScintillatorLabel(facetIdx,cellIdx_i,cellIdx_j),
                                                            0,0,0);
                    
                    logicScintillator[facetIdx][cellIdx_i][cellIdx_j] -> SetVisAttributes(ScintVisAttributes);
                    
                    // facet translation, followed by rotation, followed by another translation inside facet
                    new G4PVPlacement(0,                // no relative rotation. Rotation is defined by facet
                                      posInFacet,       // position at (x,y,z)
                                      logicScintillator[facetIdx][cellIdx_i][cellIdx_j],    // its logical volume
                                      ScintillatorLabel(facetIdx,cellIdx_i,cellIdx_j),       // its name
                                      logicFacet[facetIdx],      // its mother  volume
                                      false,           // no boolean operations
                                      0);              // copy number
                                        
                    
                    if (fdebug>1) {
                        std::cout << "constructing " << ScintillatorLabel(facetIdx,cellIdx_i,cellIdx_j)
                        << ", position in facet: ("
                        << posInFacet.x() << "," << posInFacet.y() << "," << posInFacet.z() << ")"
                        << std::endl;
                    }

                }  // end for cellIdx_j
            } // end for cellIdx_i
        }// end for facetIdx
        if (fdebug>1) {
            std::cout << "end scintillators" << std::endl;
        }
        // <--- end scintillators
    }
    
    
    // SiPMs
//    if (DoSiPMs){
//        SiPMGreen = G4Colour(0.1, 0.9, 0.2);
//        SiPMVisAttributes = new G4VisAttributes(SiPMGreen);
//
//        positionSiPM_1 = G4ThreeVector(0 , 0 ,
//                                       (sourceHolder_dz/2. + Scint_dz
//                                        + thickness_wrapping/2.
//                                        + Scint_source_dz/2.
//                                        + SiPM_dz/2.) );
//        solidSiPM_1 = new G4Box("solidSiPM_1", SiPM_dx/2., SiPM_dy/2., SiPM_dz/2.);
//        logicSiPM_1 = new G4LogicalVolume(solidSiPM_1,
//                                          fSiPMMater ,
//                                          "logicSiPM_1",
//                                          0,0,0);
//        logicSiPM_1->SetVisAttributes(SiPMVisAttributes);
//
//        new G4PVPlacement(0,              // no rotation
//                          positionSiPM_1, // at (x,y,z)
//                          logicSiPM_1,    // its logical volume
//                          "SiPM_1",       // its name
//                          logicWorld,      // its mother  volume
//                          false,           // no boolean operations
//                          0);              // copy number
//
//        positionSiPM_2 = G4ThreeVector(0 , 0 ,
//                                       -(sourceHolder_dz/2. + Scint_dz
//                                         + thickness_wrapping/2.
//                                         + Scint_source_dz/2.
//                                         + SiPM_dz/2.) );
//        solidSiPM_2 = new G4Box("solidSiPM_2", SiPM_dx/2., SiPM_dy/2., SiPM_dz/2.);
//        logicSiPM_2 = new G4LogicalVolume(solidSiPM_2,
//                                          fSiPMMater ,
//                                          "logicSiPM_2",
//                                          0,0,0);
//
//        logicSiPM_2->SetVisAttributes(SiPMVisAttributes);
//        new G4PVPlacement(0,              // no rotation
//                          positionSiPM_2, // at (x,y,z)
//                          logicSiPM_2,    // its logical volume
//                          "SiPM_2",       // its name
//                          logicWorld,      // its mother  volume
//                          false,           // no boolean operations
//                          0);              // copy number
//
//    }
    
    
//    if (DoElectronicPlates){
//        ElectronicsColor = G4Colour(0.6, 0.8, 0.6);
//        ElectronicsVisAttributes = new G4VisAttributes(ElectronicsColor);
//
//        positionElectronics_1 = G4ThreeVector(0 , 0 ,
//                                              (sourceHolder_dz/2. + Scint_dz
//                                               + thickness_wrapping/2.
//                                               + Scint_source_dz/2.
//                                               + SiPM_dz/2.
//                                               + Electronics_dz/2.) );
//        solidElectronics_1 = new G4Box("solidSiPM_1", Electronics_dx/2., Electronics_dy/2., Electronics_dz/2.);
//        logicElectronics_1 = new G4LogicalVolume(solidElectronics_1,
//                                                 fElectronicsMaterial ,
//                                                 "logicElectronics_1",
//                                                 0,0,0);
//        logicElectronics_1->SetVisAttributes(ElectronicsVisAttributes);
//
//        new G4PVPlacement(0,              // no rotation
//                          positionElectronics_1, // at (x,y,z)
//                          logicElectronics_1,    // its logical volume
//                          "Electronics_1",       // its name
//                          logicWorld,      // its mother  volume
//                          false,           // no boolean operations
//                          0);              // copy number
//
//        positionElectronics_2 = G4ThreeVector(0 , 0 ,
//                                              -(sourceHolder_dz/2. + Scint_dz
//                                                + thickness_wrapping/2.
//                                                + Scint_source_dz/2.
//                                                + SiPM_dz/2.
//                                                + Electronics_dz/2.) );
//        solidElectronics_2 = new G4Box("solidElectronics_2", Electronics_dx/2., Electronics_dy/2., Electronics_dz/2.);
//        logicElectronics_2 = new G4LogicalVolume(solidElectronics_2,
//                                                 fElectronicsMaterial ,
//                                                 "logicElectronics_2",
//                                                 0,0,0);
//
//        logicElectronics_2->SetVisAttributes(ElectronicsVisAttributes);
//        new G4PVPlacement(0,              // no rotation
//                          positionElectronics_2, // at (x,y,z)
//                          logicElectronics_2,    // its logical volume
//                          "Electronics_2",       // its name
//                          logicWorld,      // its mother  volume
//                          false,           // no boolean operations
//                          0);              // copy number
//    }
    
    
    
    
    G4cout << "done G4VPhysicalVolume* DetectorConstruction::ConstructVolumes()." << G4endl;
    // return the root volume
    return fPhysiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::PrintParameters(){
    //    G4cout << "\n Target : Length = " << G4BestUnit(fTargetLength,"Length")
    //    << " Radius = " << G4BestUnit(fTargetRadius,"Length")
    //    << " Material = " << fTargetMater->GetName();
    //    G4cout << "\n Detector : Length = " << G4BestUnit(fDetectorLength,"Length")
    //    << " Tickness = " << G4BestUnit(fDetectorThickness,"Length")
    //    << " Material = " << fDetectorMater->GetName() << G4endl;
    //    G4cout << "\n" << fTargetMater << "\n" << fDetectorMater << G4endl;
}

////....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//void DetectorConstruction::SetTargetRadius(G4double value){
//    fTargetRadius = value;
//    G4RunManager::GetRunManager()->ReinitializeGeometry();
//}

////....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//void DetectorConstruction::SetTargetLength(G4double value){
//    fTargetLength = value;
//    G4RunManager::GetRunManager()->ReinitializeGeometry();
//}

////....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//void DetectorConstruction::SetDetectorThickness(G4double value){
//    fDetectorThickness = value;
//    G4RunManager::GetRunManager()->ReinitializeGeometry();
//}
