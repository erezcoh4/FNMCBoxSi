/// \file SteppingAction.hh
/// \brief Definition of the SteppingAction class
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#ifndef SteppingAction_h
#define SteppingAction_h 1
#include "G4UserSteppingAction.hh"
#include "globals.hh"

class DetectorConstruction;
class EventAction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
class SteppingAction : public G4UserSteppingAction
{
  public:
    SteppingAction(DetectorConstruction*, EventAction*,  int _fdebug_=0);
   ~SteppingAction();

    virtual void UserSteppingAction(const G4Step*);
    int fdebug;
    // prints
    void Debug (int verobosity_level, G4String text) { if ( fdebug > verobosity_level ) std::cout << text << std::endl; }

    
    void PrintAction(G4String prefix){
        if (fdebug <= 1) return;
        std::cout
            << "\033[;31m"
            << "------------------------ " << prefix << " step --------------------------------"
            << "\033[;30m"
            << std::endl;
    }

  private:
    DetectorConstruction* fDetector;  
    EventAction* fEventAction;    
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
