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
/// \file RunAction.cc
/// \brief Implementation of the RunAction class
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
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

void RunAction::BeginOfRunAction(const G4Run*)
{
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
    
    //    //histograms
    //    //
    //    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    //    if ( analysisManager->IsActive() ) {
    //        analysisManager->OpenFile();
    //    }
    
    if (fdebug>1) std::cout << "output files" << std::endl;
    // output files
    // output particles csv file
    // copy from EventAction::EndOfEventAction(const G4Event*evt)
    char particlesfilename[50];
    char eventsfilename[50];
    sprintf(particlesfilename, "%s_particles.csv",fDetector->GetOutputFileLabel());
    sprintf(eventsfilename, "%s_events.csv",fDetector->GetOutputFileLabel());
    if (fdebug>1){
        std::cout << "particlesfilename: " << particlesfilename << std::endl;
        std::cout << "eventsfilename: " << eventsfilename << std::endl;
    }
    
    
    particlescsvfile.open(particlesfilename, std::ios_base::out);
    
    particlescsvfile
    << "eventId"            << ","
    << "trackId"            << ","
    << "parentId"           << ","
    << "PrtclName"          << ","
    << "|p_init|/MeV"       << ","
    << "CreatorProcessName" << ",";

    // track hit-position and energy deposition in scintillators
    for (int iVol=0; iVol<4; iVol++){
        particlescsvfile
        << "Edep("          <<fDetector->VolumeName(iVol)<<")/MeV"      << ",";
    }
     
    particlescsvfile
    << "PDGcode"            << ","
    << "p_init_x/MeV"       << ","
    << "p_init_y/MeV"       << ","
    << "p_init_z/MeV"       << ",";
    
    for (int iVol=0; iVol<4; iVol++){
        std::string VolName = fDetector->VolumeName(iVol);
        particlescsvfile
        << "FirstPt(" <<VolName<<")_Ek/MeV"   << ","
        << "FirstPt(" <<VolName<<")_time/ns"  << ","
        << "FirstPt(" <<VolName<<"v_x/mm"     << ","
        << "FirstPt(" <<VolName<<"v_y/mm"     << ","
        << "FirstPt(" <<VolName<<")_z/mm"     << ","
        << "LastPt("  <<VolName<<")_x/mm"     << ","
        << "LastPt("  <<VolName<<")_y/mm"     << ","
        << "LastPt("  <<VolName<<")_z/mm"     << ",";
    }
    
    // end line
    particlescsvfile << std::endl;
    
    // close file
    particlescsvfile.close();
    
    
    
    
    
    
    // output events csv file
    // copy from EventAction::EndOfEventAction(const G4Event*evt)
    eventscsvfile.open(eventsfilename, std::ios_base::out);
    
    eventscsvfile << "eventId"            << ",";
    for (int iVol=0; iVol<4; iVol++){
        eventscsvfile
        << "EdepTot("      <<fDetector->VolumeName(iVol)<<") e+e-gamma/MeV"      << ","
        << "EdepTot("      <<fDetector->VolumeName(iVol)<<")/MeV"      << ",";
    }
    
    // end line
    eventscsvfile << std::endl;
    
    // close file
    eventscsvfile.close();


    
    
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
