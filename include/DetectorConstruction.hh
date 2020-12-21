/// \file DetectorConstruction.hh
/// \brief Definition of the DetectorConstruction class
//
#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "math.h"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4RotationMatrix.hh"

// general characteristics: number of faces, and numebr of Scintillator/SiPM cells per facet
#define NFacets 3
// number of cells
// NCellsPerFacet = NCells_i x NCells_j
#define NCells_i 3
#define NCells_j 3


class G4LogicalVolume;
class G4Material;
class DetectorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction : public G4VUserDetectorConstruction {
  public:
    DetectorConstruction( int _fdebug_=0 );
   ~DetectorConstruction();

    // ----------------------------------------------------------------------------------------------------
    // dimensions
    // ----------------------------------------------------------------------------------------------------
    // scintillators
    G4double ScintillatorThickness  = 40.   *mm; // thickness can genrally differ from side length
    G4double ScintillatorSide       = 50.   *mm; // scintillator transverse cross-section is rectangular
    G4double ScintillatorPitch      = 10    *mm; // separation distance between scintillators
    G4double thickness_wrapping     = 0.1   *mm;
    // source and source holder
    G4double Source_y               = 0.    *mm;
    G4double sourceHolder_dz        = 0.5   *mm;
    // SiPMs [https://www.ketek.net/wp-content/uploads/2018/12/KETEK-PM3325-WB-D0-Datasheet.pdf]
    G4double SiPM_dx                = 3.    *mm;
    G4double SiPM_dy                = 3.    *mm;
    G4double SiPM_dz                = 0.6   *mm;
    // PCB and electronics - rough estimation
    G4double Electronics_dx         = 20.   *mm;
    G4double Electronics_dy         = 20.   *mm;
    G4double Electronics_dz         = 2     *mm;
    // sampling box
    G4double SamplingBoxSide        = 100.  *mm;
    // world
    G4double fWorldLength           = 1.2*SamplingBoxSide + 10*(ScintillatorThickness + SiPM_dz + Electronics_dz) ;
    // facets
    G4double FacetSide              = ScintillatorSide * NCells_i * 1.01;
    G4double FacetThickness         = ScintillatorThickness * 2;
    // ----------------------------------------------------------------------------------------------------
    
    
    
    
    
    
  public:
  
    virtual G4VPhysicalVolume* Construct();
    
    void SetDebug(int _fdebug_=0 ){fdebug=_fdebug_;};
                   
    void PrintParameters();
    
    
       
    
    
    
  public:
      
    
    int fdebug;
    
    
    

    
    
    

    
    G4double                GetTargetRadius();
    G4double                     GetSourceY(){ return Source_y;};
    G4double           GetSamplingBoxSide(){ return SamplingBoxSide;};
    
    G4LogicalVolume * logicScintillator[NFacets][NCells_i][NCells_j];
    G4LogicalVolume *         GetLogicWorld(){ return logicWorld; };
    G4LogicalVolume *  GetLogicSourceHolder(){ return logicSourceHolder; };
    // obselete - delete
    G4LogicalVolume * logicScint_1;
    G4LogicalVolume * logicScint_2;
    G4LogicalVolume *       GetLogicScint_1(){ return logicScint_1; };
    G4LogicalVolume *       GetLogicScint_2(){ return logicScint_2; };
    G4LogicalVolume *        GetLogicSiPM_1(){ return logicSiPM_1; };
    G4LogicalVolume *        GetLogicSiPM_2(){ return logicSiPM_2; };
    G4LogicalVolume * GetLogicElectronics_1(){ return logicElectronics_1; };
    G4LogicalVolume * GetLogicElectronics_2(){ return logicElectronics_2; };
    // ----------------
    
    
    // VolumeName
    std::string VolumeName(int iVol=0){
        std::string name = "unknown-volume";
        switch (iVol) {
            case 0:
                name = "source-holder";
                break;
            case 1:
                name = "scint-1";
                break;
            case 2:
                name = "scint-2";
                break;
            case 3:
                name = "world";
                break;
            case 4:
                name = "SiPM-1";
                break;
            case 5:
                name = "SiPM-2";
                break;
            case 6:
                name = "Electronics-1";
                break;
            case 7:
                name = "Electronics-2";
                break;

            default:
                name = "unknown-volume";
                break;
        }
        return name;
    }
    
    // output file-label
    char * GetOutputFileLabel(){
        char filelabel[100];
        sprintf(filelabel, "SamplingBoxSide%.1fmm_NFacets%d_NCellsPerFacet%d",
                            SamplingBoxSide, NFacets, NCells_i*NCells_j );
        return filelabel;
    };
    
    // scintillator label
    G4String ScintillatorLabel(int facetIdx, int cellIdx_i, int cellIdx_j){
        char scintillatorlabel[50];
        sprintf(scintillatorlabel,"Scintillator[%s][%d][%d]",FacetName(facetIdx).c_str(),cellIdx_i,cellIdx_j);
        return G4String(scintillatorlabel);
    }

    // facet name
    std::string FacetName(int facetIdx=0){
        std::string name = "unknown-facet";
        switch (facetIdx) {
            case 0:
                name = "bottom";
                break;
            case 1:
                name = "top";
                break;
            case 2:
                name = "face";
                break;
            case 3:
                name = "back";
                break;
            case 4:
                name = "left";
                break;
            case 5:
                name = "right";
                break;
            default:
                name = "unknown-volume";
                break;
        }
        return name;
    }

    // get logic scintillator
    G4LogicalVolume *  GetLogicScintillator(int facetIdx, int cellIdx_i, int cellIdx_j){
        return logicScintillator[facetIdx][cellIdx_i][cellIdx_j];
    };
    
    // facet centroid
    G4ThreeVector GetFacetCentroid(int facetIdx=0){
        G4double x = -9999,y = -9999,z = -9999;
        G4double dl = SamplingBoxSide/2 + FacetSide/2 ;
        switch (facetIdx) {
            case 0: // bottom
                x = 0;
                y = -dl;
                z = 0;
                break;
            case 1: // top
                x = 0;
                y = dl;
                z = 0;
                break;
            case 2: // front
                x = dl;
                y = 0;
                z = 0;
                break;
            case 3: // back
                x = -dl;
                y = 0;
                z = 0;
                break;
            case 4: // left
                x = 0;
                y = 0;
                z = -dl;
                break;
            case 5: // right
                x = 0;
                y = 0;
                z = dl;
                break;
        }
        return G4ThreeVector(x,y,z);
    }
    
    // facet rotation
    G4Rotate3D GetFacetRot3D(int facetIdx){
        G4Rotate3D FacetRot3D;
        switch (facetIdx) {
            case 0: // bottom
                FacetRot3D = G4Rotate3D( 90.*deg , G4ThreeVector(1,0,0) );
                break;
            case 1: // top
                FacetRot3D = G4Rotate3D( -90.*deg , G4ThreeVector(1,0,0) );
                break;
            case 2: // front
                FacetRot3D = G4Rotate3D( -90.*deg , G4ThreeVector(0,1,0) );
                break;
            case 3: // back
                FacetRot3D = G4Rotate3D( 90.*deg , G4ThreeVector(0,1,0) );
                break;
            case 4: // left
                break;
            case 5: // right
                FacetRot3D = G4Rotate3D( 180.*deg , G4ThreeVector(1,0,0) );
                break;
        }
    return FacetRot3D;
    }
        
    // compute scintillator position in facet coordinates
    G4ThreeVector GetScintillatorPositionInFacet( int cellIdx_i, int cellIdx_j ){
        // in the facet axes frame (i,j,k)
        // the scintillator is located in k=0
        // and i,j are determined by the cell indices
        G4double dl_i = (ScintillatorSide + ScintillatorPitch) * (cellIdx_i-NCells_j/2) ;//* (NCells_i-1)/2;
        G4double dl_j = (ScintillatorSide + ScintillatorPitch) * (cellIdx_j-NCells_j/2); //(float(cellIdx_j)/NCells_j - 1./2);
        G4double dl_k = 0;
        return G4ThreeVector( dl_i , dl_j , dl_k );
    }

    
    
    
    
//    1
//    0
//    40*0/2
//    2
//    -20 20
//    40/2*-1 40/2*1
//    3
//    -40 0 +40
//    40/2*-2 40/2*0 40/2*2
//    4
//    -60 -20 20 60
//    40/2*-3 40/2*-1 40/2*1 40/2*3
//    5
//    -80 -40 0 40 80
//    40/2*-4 40/2*-2 40/2*0 40/2*2 40/2*4

    
    
    
    
    
    
    
    
    
    
    
  private:
    
    G4double           fTargetLength;
    G4double           fTargetRadius;
    G4double           fDetectorLength;
    G4double           fDetectorThickness;
    
    
    
    
    G4VPhysicalVolume* fPhysiWorld;
    DetectorMessenger* fDetectorMessenger;
    
            
    
    G4Box * solidSourceHolder;
    G4Box * solidFacet[NFacets];
    G4Box * solidScintillator[NFacets][NCells_i][NCells_j];
    G4Box * solidSiPM[NFacets][NCells_i][NCells_j];
    G4Box * solidElectronics[NFacets][NCells_i][NCells_j];
    // obselete - delete
    G4Box * solidScint_1, * solidScint_2;
    G4Box * solidSiPM_1, * solidSiPM_2;
    G4Box * solidElectronics_1, * solidElectronics_2;
    // ----------------
    
    G4Sphere * sourcePlaceHolder;
    
    
    
    G4Material*        fTargetMater;
    G4Material*        fPlasticMater;
    G4Material*        fSiPMMater;
    G4Material*        fElectronicsMaterial;
    G4Material*        fDetectorMater;
    G4Material*        fWorldMater;

    G4ThreeVector posSourceHolder;
    // obselete - delete
    G4ThreeVector positionScint_1, positionScint_2;
    G4ThreeVector positionSiPM_1, positionSiPM_2;
    G4ThreeVector positionElectronics_1, positionElectronics_2;
    // ----------------
    G4ThreeVector positionFacet[NFacets];
    G4ThreeVector positionScintillator[NFacets][NCells_i][NCells_j];
    G4ThreeVector positionSiPM[NFacets][NCells_i][NCells_j];
    G4ThreeVector positionElectronics[NFacets][NCells_i][NCells_j];

    // obselete - delete
    G4LogicalVolume * logicSiPM_1, * logicSiPM_2;
    G4LogicalVolume * logicElectronics_1, * logicElectronics_2;
    // ----------------
    G4LogicalVolume * logicFacet[NFacets];
    G4LogicalVolume * logicSiPM[NFacets][NCells_i][NCells_j];
    G4LogicalVolume * logicElectronics[NFacets][NCells_i][NCells_j];
    G4LogicalVolume * logicWorld;
    G4LogicalVolume * logicSourceHolder;
    G4LogicalVolume * logicSourcePlaceHolder;

    G4Colour SourceHolderColor;
    G4Colour SourceRed;
    G4Colour ScintBlue;
    G4Colour SiPMGreen;
    G4Colour ElectronicsColor;

    G4VisAttributes * SourceHolderVisAttributes;
    G4VisAttributes * sourceVisAttributes;
    G4VisAttributes * FacetVisAttributes;
    G4VisAttributes * ScintVisAttributes;
    G4VisAttributes * SiPMVisAttributes;
    G4VisAttributes * ElectronicsVisAttributes;

    
  private:
    
    void               DefineMaterials();
    G4VPhysicalVolume* ConstructVolumes();     
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#endif

