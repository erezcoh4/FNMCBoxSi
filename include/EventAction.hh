/// \file EventAction.hh
/// \brief Definition of the EventAction class
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#ifndef EventAction_h
#define EventAction_h 1
#define NMAXtracks 1000
#define NVOLUMES 8
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
    EventAction(DetectorConstruction * ,int _fdebug_=0);
   ~EventAction();

  public:
    virtual void BeginOfEventAction(const G4Event*);
    virtual void   EndOfEventAction(const G4Event*);
    void InitialiseArrays();
    
    void AddEdep (G4int iVol,
                  G4double edep,
                  G4double time,
                  G4int trackId,
                  G4int PDGCode,
                  int fdebug = 0);
           
    
    void SetFirstPointInVolume( G4int iVol, G4int trackId,
                               G4ThreeVector pos,
                               G4double time,
                               G4double Ek,
                               G4String fProcessName){
        int fdebug=0;
        if (fdebug>1){
            std::cout
            << "pos: (" << pos.x() << "," << pos.y() << ","<< pos.z() << ")"
            << ", time:" << time
            << "ns, Ek:" << Ek <<  "MeV, fProcessName:" << fProcessName
            << std::endl;
        }
        
        FirstPointInVolume[iVol][trackId]       = pos;
        FirstPointInVolumeTime[iVol][trackId]   = time;
        FirstPointInVolumeEk[iVol][trackId]     = Ek;
        
        ProcessName.at(trackId) = fProcessName;
        
    };
    void SetLastPointInVolume( G4int iVol, G4int trackId, G4ThreeVector pos, G4double time ){
        LastPointInVolume[iVol][trackId] = pos;
        LastPointInVolumeTime[iVol][trackId] = time;
    };
    
    std::ofstream particlescsvfile, eventsscsvfile;
    int fdebug;
    void SetDebug(int _fdebug_=0){fdebug=_fdebug_;};
    
    char particlesfilename[50];
    char eventsfilename[50];

    
  private:
    G4double fEdep1,   fEdep2;
    G4double fWeight1, fWeight2;
    G4double fTime0;
    
    // 2-D array of energy depostion per volume per trak
    // first dimension is volume number
    // (0 = source holder, 1=scintllator 1, 2 = scintillator 2, 3 = world)
    // second dimension is track Id
    G4double fEdep[NVOLUMES][NMAXtracks];
    
    // total energy deposition in volume per event
    G4double fEdepTotVol[NVOLUMES];
    G4double fEdepTotVol_e_g[NVOLUMES]; // energy deposit by e+,e-,gamma
    
    G4double FirstPointInVolumeTime[NVOLUMES][NMAXtracks];
    G4double LastPointInVolumeTime[NVOLUMES][NMAXtracks];
    G4double FirstPointInVolumeEk[NVOLUMES][NMAXtracks];

    
    // first and last point of each track in each volume
    G4ThreeVector FirstPointInVolume[NVOLUMES][NMAXtracks];
    G4ThreeVector LastPointInVolume[NVOLUMES][NMAXtracks];

    
    std::vector<G4String> ProcessName;
    DetectorConstruction* fDetector;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
