/// \file RunAction.hh
/// \brief Definition of the RunAction class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"
#include <fstream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction;
class Run;
class PrimaryGeneratorAction;
class HistoManager;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class RunAction : public G4UserRunAction
{
  public:
    int fdebug;
    RunAction(DetectorConstruction*, PrimaryGeneratorAction*, int _fdebug_ = 0);
   ~RunAction();

  public:
    virtual G4Run* GenerateRun();  
    virtual void BeginOfRunAction(const G4Run*);
    virtual void   EndOfRunAction(const G4Run*);
    std::ofstream particlescsvfile, eventscsvfile, primariescsvfile;
    
    
  private:
    DetectorConstruction*      fDetector;
    PrimaryGeneratorAction*    fPrimary;
    Run*                       fRun;    
    HistoManager*              fHistoManager;
        
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

