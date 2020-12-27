/// \file SteppingAction.cc
/// \brief Implementation of the SteppingAction class
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#include "SteppingAction.hh"
#include "DetectorConstruction.hh"
#include "Run.hh"
#include "EventAction.hh"
#include "HistoManager.hh"

#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
SteppingAction::SteppingAction(DetectorConstruction* det, EventAction* event, int _fdebug_)
: G4UserSteppingAction(), fDetector(det),  fdebug(_fdebug_), fEventAction(event) { }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
SteppingAction::~SteppingAction(){ }



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void SteppingAction::UserSteppingAction(const G4Step* aStep){
    PrintAction("begin");
    Run* run = static_cast<Run*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
    
    //which volume ?
    G4LogicalVolume* lVolume = aStep -> GetPreStepPoint() -> GetTouchableHandle() -> GetVolume() -> GetLogicalVolume();
    // world: iVol = 0
    // scintillators: iVol = 1
    G4int iVol = 0;
    
    
    // count processes
    const G4StepPoint* startPoint = aStep->GetPreStepPoint();
    const G4StepPoint* endPoint = aStep->GetPostStepPoint();
    
    
    G4double stepLength = aStep -> GetStepLength();
    if (fabs(stepLength / CLHEP::nm)<1) {PrintAction("step length < 1 nm, end"); return;};
    
    const G4VProcess* process   = endPoint->GetProcessDefinedStep();
    //    run->CountProcesses(process, iVol);
        
    G4double trackId = aStep -> GetTrack() -> GetTrackID();
    G4double parentId = aStep -> GetTrack() -> GetParentID();
    G4double edepStep = aStep -> GetTotalEnergyDeposit();
    G4double time = aStep -> GetPreStepPoint() -> GetGlobalTime();
    G4double Ek = aStep -> GetPreStepPoint() -> GetKineticEnergy ();
    G4bool   isFirstStepInVolume = aStep->IsFirstStepInVolume ();
    //    G4bool   isLastStepInVolume = aStep->IsLastStepInVolume ();
    auto creator = aStep -> GetTrack() -> GetCreatorProcess ();
    G4String ProcessName;
    if (creator==nullptr){
        ProcessName = "Primary";
    } else {
        ProcessName = creator -> GetProcessName ();
    }
    //    G4String ProcessName = aStep -> GetTrack() -> GetCreatorProcess () -> GetProcessName ();
        
    if (fdebug>1){
        std::cout << "process: " << process->GetProcessName() << std::endl;
        std::cout << "stepLength: " << stepLength / CLHEP::mm << " mm" << std::endl;
        
        std::cout << "("
        << startPoint->GetPosition().x()/cm << ","
        << startPoint->GetPosition().y()/cm << ","
        << startPoint->GetPosition().z()/cm << ")"
        << " -> ";
        
        std::cout << "("
        << endPoint->GetPosition().x()/cm << ","
        << endPoint->GetPosition().y()/cm << ","
        << endPoint->GetPosition().z()/cm << ") cm"
        << std::endl;
        
        std::cout
        << "iVol: " << iVol << ", "
        << "time: " << time/ns << " ns, "
        << std::endl
        << "track id: "<< aStep->GetTrack()->GetTrackID() << ", "
        << "parentId: "<< aStep->GetTrack()->GetParentID() << ", "
        << "particle: "<< aStep->GetTrack()->GetParticleDefinition ()->GetParticleName() << ", "
        << std::endl
        << "total energy deposit: "<< aStep->GetTotalEnergyDeposit ()/MeV << " MeV, "
        << std::endl
        << "IsFirstStepInVolume?: "<< aStep->IsFirstStepInVolume () << ", "
        << "IsLastStepInVolume?: " << aStep->IsLastStepInVolume () << ", "
        << std::endl;
    }
    
    for (int facetIdx=0; facetIdx < NFacets; facetIdx++){
        for (int cellIdx_i=0; cellIdx_i < NCells_i; cellIdx_i++){
            for (int cellIdx_j=0; cellIdx_j < NCells_j; cellIdx_j++){
                if (lVolume == fDetector -> GetLogicScintillator(facetIdx,cellIdx_i,cellIdx_j) ){
                    iVol = 1;
                    // stream this information into a dedicated buffer at the "event-action"
                    fEventAction -> AddEdepScintillator( facetIdx, cellIdx_i, cellIdx_j,
                                                        edepStep, time, trackId ,
                                                        aStep->GetTrack()->GetParticleDefinition()->GetPDGEncoding());
                    
                    // check if its the first or last point in each volume,
                    // and if so, stream them into the event-action buffer
                    if (isFirstStepInVolume){
                        fEventAction -> SetFirstPointInScintillator( facetIdx, cellIdx_i, cellIdx_j,
                                                              trackId, parentId,
                                                              aStep->GetPreStepPoint()->GetPosition() ,
                                                              time,
                                                              Ek,
                                                              ProcessName );
                    }
                }
            }
        }
    }
    
    
    
    //    if (isLastStepInVolume){
    //        fEventAction -> SetLastPointInVolume( iVol,
    //                                             trackId,
    //                                             aStep->GetPostStepPoint()->GetPosition() ,
    //                                             time );
    //    }
    PrintAction("end");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

