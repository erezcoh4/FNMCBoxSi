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
/// \file EventAction.cc
/// \brief Implementation of the EventAction class
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "EventAction.hh"

#include "Run.hh"
#include "HistoManager.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include <stdio.h>
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction(DetectorConstruction* det,int _fdebug_)
:G4UserEventAction(),
fDetector(det), fEdep1(0.), fEdep2(0.), fWeight1(0.), fWeight2(0.),
fTime0(-1*s)
{
    InitialiseArrays();
    SetDebug(_fdebug_);
            
    sprintf(particlesfilename, "%s_particles.csv",fDetector->GetOutputFileLabel());
    sprintf(eventsfilename, "%s_events.csv",fDetector->GetOutputFileLabel());
    if (fdebug>1){
        std::cout << "particlesfilename: " << particlesfilename << std::endl;
        std::cout << "eventsfilename: " << eventsfilename << std::endl;
    }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
EventAction::~EventAction()
{ }


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void EventAction::InitialiseArrays(){
    
    if (fdebug>1) std::cout << "EventAction::InitialiseArrays()" << std::endl;
    // initialize fEdep
    // fEdep is a 2-D array
    // first dimension is volume number
    // (0 = source holder, 1=scintllator 1, 2 = scintillator 2, 3 = world)
    // second dimension is track Id
    for (int iVol=0; iVol<NVOLUMES; iVol++){
        fEdepTotVol[iVol] = 0; // total energy deposition in volume in this event...
        fEdepTotVol_e_g[iVol] = 0; // total energy deposition in volume by e+,e-,gamma
    }
    for (int trackId=0; trackId<NMAXtracks; trackId++){
        ProcessName.push_back("");
        for (int iVol=0; iVol<NVOLUMES; iVol++){
            fEdep[iVol][trackId] = 0;
            FirstPointInVolume[iVol][trackId] = G4ThreeVector(-999,-999,-999);
            LastPointInVolume[iVol][trackId] = G4ThreeVector(-999,-999,-999);
            FirstPointInVolumeTime[iVol][trackId] = -999.*microsecond;
            LastPointInVolumeTime[iVol][trackId] = -999.*microsecond;
            FirstPointInVolumeEk[iVol][trackId] = -999.*MeV;
            
        }
    }
    if (fdebug>1) std::cout << "done EventAction::InitialiseArrays()" << std::endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event*)
{
    if (fdebug>1) std::cout << "EventAction::BeginOfEventAction(const G4Event*)" << std::endl;
    fEdep1 = fEdep2 = fWeight1 = fWeight2 = 0.;
    fTime0 = -1*s;
    if (fdebug>1) std::cout << "done EventAction::BeginOfEventAction(const G4Event*)" << std::endl;
    InitialiseArrays();
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event*evt) {
    G4int eventId = evt -> GetEventID();
    // --------------------------------
    // write a csv file for each particle that was emitted from the source:
    //
    // event ID
    // partice type (PDG code)
    // track ID
    // parent ID
    // initial momentum [MeV/c]
    // final momentum [MeV/c]
    //
    // hit time scintillator 1 [ps] (-999 if did not hit scintillator)
    // start point in scintillator 1 [cm] ((-999,-999,-999) if did not hit scintillator)
    // end point in scintillator 1 [cm] ((-999,-999,-999) if did not hit scintillator)
    // energy deposition in scintilaltor 1 [MeV] (0 if did not hit scintillator)
    // underwent compton in scintillator 1 ?
    //
    // hit time scintillator 2 [ps] (-999 if did not hit scintillator)
    // start point in scintillator 2 [cm] ((-999,-999,-999) if did not hit scintillator)
    // end point in scintillator 2 [cm] ((-999,-999,-999) if did not hit scintillator)
    // energy deposition in scintilaltor 2 [MeV] (0 if did not hit scintillator)
    // underwent compton in scintillator 2 ?
    //
    bool do_particles_file = true;
    bool do_events_file = true;
    
    
    if (fdebug>0) std::cout << "EventAction::EndOfEventAction(const G4Event*evt)" << std::endl;
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    
    
    // extract event data and write to
    // output files
    if (fdebug>1) std::cout << "event " << evt -> GetEventID() << ",opening output csv and write data" << std::endl;
    if (fdebug>1) std::cout << "trajectory container: " << std::endl;
    
    // trajectory container
    G4TrajectoryContainer * trajCont = evt -> GetTrajectoryContainer();
    if (fdebug>1) std::cout << "got trajCont, of size " << sizeof(trajCont)/sizeof(trajCont[0]) << std::endl;
    if(trajCont==0) {
        if (fdebug>1) std::cout << "trajCont=0, returning" << std::endl;
        return;
    }
    
//        std::vector< G4VTrajectory * > * trajectories = trajCont -> GetVector (); // delete by Dec-30,2020
    TrajectoryVector * trajectories = trajCont -> GetVector ();
    
    if (fdebug>1) std::cout << "for (auto traj:*trajectories)" << std::endl;
    if (fdebug>1) std::cout << "output particles csv file " << particlesfilename << std::endl;
    
    if (do_particles_file){
        // output particles csv file
        // copy from EventAction::EndOfEventAction(const G4Event*evt)
        particlescsvfile.open(particlesfilename, std::ios_base::app);
        for (auto traj:*trajectories){
            
            G4String ParticleName = traj->GetParticleName() ;
            G4int PDGcode = traj->GetPDGEncoding();
            G4int trackId = traj->GetTrackID ();
            G4int parentId = traj->GetParentID ();
            G4ThreeVector p_init = traj->GetInitialMomentum ();
            
            // for data saving, omit neutrinos and heavy isotopes that promptly decay
            if (
                ParticleName=="nu_e" ||
                ParticleName=="Na22" ||
                ParticleName=="Ne22" ||
                ParticleName=="Ne22[1274.577]"
                )
                continue;
            
            // omit tracks that did not deposit energy at all - neither in scintillator 1 nor in scintillator 2
            if ((fEdep[1][trackId]==0) && (fEdep[2][trackId])==0) {
                continue;
                // ToDo: count these tracks for efficiency ?
            }
            
            //  trajectory points can not help, as they do not provide access for the energy deposition etc.
            // for each particle, we dedicate a line in the csv file
            particlescsvfile
            << eventId          << ","
            << trackId          << ","
            << parentId         << ","
            << ParticleName     << ","
            << p_init.mag()/MeV             << ","
            << ProcessName.at(trackId)      << ",";
            
            
            // track hit-position and energy deposition in scintillators
            for (int iVol=0; iVol<NVOLUMES; iVol++){
                particlescsvfile
                << fEdep[iVol][trackId]/MeV << ",";
            }
            
            particlescsvfile
            << PDGcode                      << ","
            << p_init.x()/MeV               << ","
            << p_init.y()/MeV               << ","
            << p_init.z()/MeV               << ",";
            
            for (int iVol=0; iVol<NVOLUMES; iVol++){
                particlescsvfile
                << FirstPointInVolumeEk[iVol][trackId]/MeV << ","
                << FirstPointInVolumeTime[iVol][trackId]/ns << ","
                << FirstPointInVolume[iVol][trackId].x()/mm << ","
                << FirstPointInVolume[iVol][trackId].y()/mm << ","
                << FirstPointInVolume[iVol][trackId].z()/mm << ","
                << LastPointInVolume[iVol][trackId].x()/mm << ","
                << LastPointInVolume[iVol][trackId].y()/mm << ","
                << LastPointInVolume[iVol][trackId].z()/mm << ",";
            }
            
            // end line
            particlescsvfile << std::endl;
            
        }
        particlescsvfile.close();
        if (fdebug>1) std::cout << "wrote to particlescsvfile " << particlesfilename << std::endl;
    }
    
    if (do_events_file) {
        // output events csv file
        // copy from EventAction::EndOfEventAction(const G4Event*evt)
        eventsscsvfile.open(eventsfilename, std::ios_base::app);
        
        eventsscsvfile
        << eventId          << ",";
        for (int iVol=0; iVol<4; iVol++){
            eventsscsvfile
            << fEdepTotVol_e_g[iVol]/MeV    << "," // energy deposit by e+,e-,gamma
            << fEdepTotVol[iVol]/MeV        << ","; // total energy deposition in event (over all tracks)
        }
        eventsscsvfile << std::endl;
        // close file
        eventsscsvfile.close();
        if (fdebug>1) std::cout << "wrote to eventsscsvfile " << eventsfilename << std::endl;
    }
    
    if (fdebug>1) std::cout << "Done EventAction::EndOfEventAction(const G4Event*evt) " << std::endl;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......




//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// AddEdep modified by Erez, Nov-25, 2020
void EventAction::AddEdep(G4int iVol,
                          G4double edep,
                          G4double time,
                          G4int trackId,
                          G4int PDGcode,
                          int fdebug ){
    
    if (fdebug>1) std::cout << "EventAction::AddEdep()" << std::endl;
    
    // omit tracks of id > NMAXtracks
    if (trackId >= NMAXtracks){
        std::cout
        << "(trackId = " << trackId << ") >= (NMAXtracks = " << NMAXtracks << ")"
        <<
        "returning from EventAction::AddEdep() " << std::endl;
        return;
    }
    // omit volume of iVol > NVOLUMES
    if (iVol >= NVOLUMES){
        std::cout
        << "(iVol = " << iVol << ") >= (NVOLUMES = " << NVOLUMES << ")"
        <<
        "returning from EventAction::AddEdep() " << std::endl;
        return;
    }
    // fEdep is a 2-D array
    // first dimension is volume number
    // (0 = source holder, 1=scintllator 1, 2 = scintillator 2, 3 = world)
    // second dimension is track Id
    fEdep[iVol][trackId] += edep;
    fEdepTotVol[iVol] += edep;
    
    if ((PDGcode==22) || (PDGcode==11)  || (PDGcode==-11) ){
        fEdepTotVol_e_g[iVol] += edep;
    }
    if (fdebug>1) std::cout << "fEdep["<<iVol<<"]["<<trackId<<"] = " << fEdep[iVol][trackId]/keV << " keV" << std::endl;
    if (fdebug>1) std::cout << "done EventAction::AddEdep()" << std::endl;
    
    // original
    //    // initialize t0
    //    if (fTime0 < 0.) fTime0 = time;
    //
    //    // out of time window ?
    //    const G4double TimeWindow (1*microsecond);
    //    if (std::fabs(time - fTime0) > TimeWindow) return;
    //    if (iVol == 1) { fEdep1 += edep; fWeight1 += edep*weight;}
    //    if (iVol == 2) { fEdep2 += edep; fWeight2 += edep*weight;}
}


