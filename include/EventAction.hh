/// \file EventAction.hh
/// \brief Definition of the EventAction class
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#ifndef EventAction_h
#define EventAction_h 1
#define NMAXtracks 1000

//#define NVOLUMES 8
#include "G4UserEventAction.hh"
#include "DetectorConstruction.hh"

#include "globals.hh"
#include "G4ThreeVector.hh"
#include <fstream>
#include <vector>
#include <string>
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class EventAction : public G4UserEventAction
{
public:
    
    G4double ThresholdEdepScintillator = 200.*keV;
    
    EventAction(DetectorConstruction * ,int _fdebug_=0);
    ~EventAction();
    
    
    
    
public:
    std::ofstream particlescsvfile, primariescsvfile, eventsscsvfile;
    char particlesfilename[100];
    char primariesfilename[100];
    char eventsfilename[100];
    G4int NScintillatorsFiredPerPrimary;
    
    // prints
    int fdebug;
    G4int eventId;
    void SetDebug(int _fdebug_=0){fdebug=_fdebug_;};
    void Debug (int verobosity_level, G4String text) { if ( fdebug > verobosity_level ) std::cout << text << std::endl; }
    
    virtual void BeginOfEventAction(const G4Event*);
    virtual void   EndOfEventAction(const G4Event*);
    void InitialiseArrays();
    
    void AddEdepScintillator (G4int facetIdx, G4int cellIdx_i, G4int cellIdx_j,
                              G4double edep,
                              G4double time,
                              G4int trackId,
                              G4int PDGCode);
    
    
    void SetFirstPointInScintillator( G4int facetIdx, G4int cellIdx_i, G4int cellIdx_j,
                                     G4int trackId, G4int parentId,
                                     G4ThreeVector pos,
                                     G4double time,
                                     G4double Ek,
                                     G4String fProcessName){
        if (fdebug>1){
            std::cout
            << "pos: (" << pos.x() << "," << pos.y() << ","<< pos.z() << ")"
            << ", time:" << time
            << "ns, Ek:" << Ek <<  "MeV, fProcessName:" << fProcessName
            << std::endl;
        }
        
        //            FirstPointInVolume[iVol][trackId]       = pos;
        //            FirstPointInVolumeTime[iVol][trackId]   = time;
        //            FirstPointInVolumeEk[iVol][trackId]     = Ek;
        
        ProcessName.at(trackId) = fProcessName;
        
    };
    //    void SetLastPointInVolume( G4int iVol, G4int trackId, G4ThreeVector pos, G4double time ){
    //        LastPointInVolume[iVol][trackId] = pos;
    //        LastPointInVolumeTime[iVol][trackId] = time;
    //    };
    
    
    void PrintAction(G4String prefix){
        if (fdebug <= 1) return;
        std::cout
        << "\033[;32m"
        << "------------------------ " << prefix << " event " << eventId << " --------------------------------"
        << "\033[;30m"
        << std::endl;
    }
    
private:
    //    G4double fEdep1,   fEdep2;
    //    G4double fWeight1, fWeight2;
    // G4double fTime0;
    
    // multidimensional array of energy depostion per volume per trak
    // second dimension is track Id
    //    G4double fEdep[NVOLUMES][NMAXtracks];
    G4double fEdepScintillator[NFacets][NCells_i][NCells_j][NMAXtracks];
    
    // total energy deposition in volume per event
    //    G4double fEdepTotScintillator[NFacets][NCells_i][NCells_j];
    
    //    G4double FirstPointInVolumeTime[NVOLUMES][NMAXtracks];
    //    G4double LastPointInVolumeTime[NVOLUMES][NMAXtracks];
    //    G4double FirstPointInVolumeEk[NVOLUMES][NMAXtracks];
    
    
    //    // first and last point of each track in each volume
    //    G4ThreeVector FirstPointInVolume[NVOLUMES][NMAXtracks];
    //    G4ThreeVector LastPointInVolume[NVOLUMES][NMAXtracks];
    
    
    std::vector<G4String> ProcessName;
    DetectorConstruction* fDetector;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


