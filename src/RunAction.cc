/// \file RunAction.cc
/// \brief Implementation of the RunAction class
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#include "RunAction.hh"
#include "Run.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "HistoManager.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "Randomize.hh"
#include <iomanip>
#include <stdio.h>
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(DetectorConstruction* det, PrimaryGeneratorAction* prim, int _fdebug_)
: G4UserRunAction(),
fDetector(det), fPrimary(prim), fRun(0), fHistoManager(0), fdebug(_fdebug_)
{
//    // Book predefined histograms
//    fHistoManager = new HistoManager();
    if (fdebug>1) std::cout << "RunAction::RunAction()" << std::endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{
//    delete fHistoManager;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Run* RunAction::GenerateRun()
{ 
    fRun = new Run(fDetector);
    return fRun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RunAction::BeginOfRunAction(const G4Run*) {
    
    bool do_particles_file = false;
    bool do_primaries_file = true;
    bool do_events_file = false;

    
    // ToDo: initialise global time somehow?
    if (fdebug>1) std::cout << "RunAction::BeginOfRunAction(const G4Run*)" << std::endl;
    // save Rndm status
    G4RunManager::GetRunManager()->SetRandomNumberStore(false);
    if (isMaster && fdebug>1) G4Random::showEngineStatus();
    
    if (fdebug>1) std::cout << "keep run condition" << std::endl;
    // keep run condition
    if (fPrimary) {
        G4ParticleDefinition* particle
        = fPrimary->GetParticleGun()->GetParticleDefinition();
        G4double energy = fPrimary->GetParticleGun()->GetParticleEnergy();
        fRun->SetPrimary(particle, energy);
    }
    
    
    if (fdebug>1) std::cout << "output files" << std::endl;
    // output files
    // output particles csv file
    // copy from EventAction::EndOfEventAction(const G4Event*evt)
    char particlesfilename[100];
    sprintf(particlesfilename, "%s_particles.csv",fDetector->GetOutputFileLabel().c_str());
    char primariesfilename[100];
    sprintf(primariesfilename, "%s_primaries.csv",fDetector->GetOutputFileLabel().c_str());
    char eventsfilename[100];
    sprintf(eventsfilename, "%s_events.csv",fDetector->GetOutputFileLabel().c_str());
    
    if (fdebug>1){
        std::cout << "file label: " << fDetector->GetOutputFileLabel()
        << std::endl
        << "particlesfilename: " << particlesfilename
        << std::endl
        << "eventsfilename: " << eventsfilename
        << std::endl;
    }
    
    if (do_particles_file){
    particlescsvfile.open(particlesfilename, std::ios_base::out);
    
    particlescsvfile
    << "eventId"                << ","
    << "trackId"                << ","
    << "parentId"               << ","
    << "PrtclName"              << ","
    << "|p_init|/MeV"           << ","
    << "theta/deg"              << ","
    << "phi/deg"                << ","
    << "CreatorProcessName"     << ","
    << "NScintillatorsFiredPerPrimary"    << ",";
    
    // track hit-position and energy deposition in scintillators
    for (int facetIdx=0; facetIdx < NFacets; facetIdx++){
        for (int cellIdx_i=0; cellIdx_i < NCells_i; cellIdx_i++){
            for (int cellIdx_j=0; cellIdx_j < NCells_j; cellIdx_j++){
                particlescsvfile << "Edep("  << fDetector->ScintillatorLabel(facetIdx,cellIdx_i,cellIdx_j)<<")/MeV"      << ",";
            }
        }
    }
    
    particlescsvfile
    << "PDGcode"            << ","
    << "p_init_x/MeV"       << ","
    << "p_init_y/MeV"       << ","
    << "p_init_z/MeV"       << ",";
    
    // end line
    particlescsvfile << std::endl;
    // close file
    particlescsvfile.close();
    }
    
    if (do_primaries_file){
        primariescsvfile.open(primariesfilename, std::ios_base::out);
        
        primariescsvfile
        << "eventId"                << ","
        << "trackId"                << ","
        << "parentId"               << ","
        << "PrtclName"              << ","
        << "|p_init|/MeV"           << ","
        << "CreatorProcessName"     << ","
        << "NScintillatorsFiredPerPrimary"    << ","
        << "PDGcode"                << ",";
        // end line
        primariescsvfile << std::endl;
        // close file
        primariescsvfile.close();
    }
    
    if (do_events_file){
//    // output events csv file
//    // copy from EventAction::EndOfEventAction(const G4Event*evt)
//    eventscsvfile.open(eventsfilename, std::ios_base::out);
//
//    eventscsvfile << "eventId"            << ",";
//    for (int iVol=0; iVol<4; iVol++){
//        eventscsvfile
//        << "EdepTot("      <<fDetector->VolumeName(iVol)<<") e+e-gamma/MeV"      << ","
//        << "EdepTot("      <<fDetector->VolumeName(iVol)<<")/MeV"      << ",";
//    }
//
//    // end line
//    eventscsvfile << std::endl;
//
//    // close file
//    eventscsvfile.close();

    }
    
    
    if (fdebug>1) std::cout << "done RunAction::BeginOfRunAction(const G4Run*)" << std::endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run*)
{
    if (fdebug>1) std::cout << "RunAction::EndOfRunAction(const G4Run*)" << std::endl;
    if (isMaster) fRun->EndOfRun();
    
    //    //save histograms
    //    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    //    if ( analysisManager->IsActive() ) {
    //        analysisManager->Write();
    //        analysisManager->CloseFile();
    //    }
    
    // show Rndm status
    if (isMaster && fdebug>1) G4Random::showEngineStatus();
    if (fdebug>1) std::cout << "done RunAction::EndOfRunAction(const G4Run*)" << std::endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
