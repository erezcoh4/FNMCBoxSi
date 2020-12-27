/// \file EventAction.cc
/// \brief Implementation of the EventAction class
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
fDetector(det){
    InitialiseArrays();
    SetDebug(_fdebug_);
    
    sprintf(particlesfilename, "%s_particles.csv",fDetector->GetOutputFileLabel().c_str());
    sprintf(primariesfilename, "%s_primaries.csv",fDetector->GetOutputFileLabel().c_str());
    sprintf(eventsfilename, "%s_events.csv",fDetector->GetOutputFileLabel().c_str());
    if (fdebug>1){
        std::cout << "particlesfilename: " << particlesfilename << std::endl;
        std::cout << "primariesfilename: " << primariesfilename << std::endl;
        std::cout << "eventsfilename: " << eventsfilename << std::endl;
    }
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
EventAction::~EventAction(){ }


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void EventAction::InitialiseArrays(){
    
    if (fdebug>1) std::cout << "EventAction::InitialiseArrays()" << std::endl;
    // initialize fEdep
    // fEdepScintillator is a N-D array
    for (int trackId=0; trackId<NMAXtracks; trackId++){
        ProcessName.push_back("");
        for (int facetIdx=0; facetIdx < NFacets; facetIdx++){
            for (int cellIdx_i=0; cellIdx_i < NCells_i; cellIdx_i++){
                for (int cellIdx_j=0; cellIdx_j < NCells_j; cellIdx_j++){
                    fEdepScintillator[facetIdx][cellIdx_i][cellIdx_j][trackId] = 0;
                }
            }
        }
    }
    if (fdebug>1) std::cout << "done EventAction::InitialiseArrays()" << std::endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void EventAction::BeginOfEventAction(const G4Event*){
    if (fdebug>1) std::cout << "EventAction::BeginOfEventAction(const G4Event*)" << std::endl;
    InitialiseArrays();
    if (fdebug>1) std::cout << "done EventAction::BeginOfEventAction(const G4Event*)" << std::endl;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void EventAction::EndOfEventAction(const G4Event*evt) {
    eventId = evt -> GetEventID();
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
    PrintAction("begin processing");
    bool do_particles_file = false;
    bool do_primaries_file = true;
    bool do_events_file = false;
    
    
    // trajectory container
    G4TrajectoryContainer * trajCont = evt -> GetTrajectoryContainer();
    
    // extract event data and write to output files
    if (fdebug>1) {
        std::cout << "event " << evt -> GetEventID() << ",opening output csv and write data" << std::endl;
        std::cout << "trajectory container: " << std::endl;
        std::cout << "got trajCont, of size " << sizeof(trajCont)/sizeof(trajCont[0]) << std::endl;
        std::cout << "trajCont: " << trajCont << std::endl;
    }
    
    TrajectoryVector * trajectories = trajCont -> GetVector ();
    if(trajCont==0) {
        PrintAction("trajCont=0, not processing");
        return;
    }
    
    // process data
    // count in how many scintillators have
    // the primary neutron knocked out a proton.
    // namely, find all protons that were knocked out by the primary neutron,
    // and if their energy deposition in a given scintillator is above the threshold, determine that the scintillator "fired"
    // carbon knockout to first order does not contribute to neutron detection efficiency
    NScintillatorsFiredPerPrimary = 0;
    // (1) step over all scintillators,
    // (2) for each one, determine if it fired or did not fire,
    // (3) according to energy deposition by protons knocked out by the primary neutron
    for (int facetIdx=0; facetIdx < NFacets; facetIdx++){
        for (int cellIdx_i=0; cellIdx_i < NCells_i; cellIdx_i++){
            for (int cellIdx_j=0; cellIdx_j < NCells_j; cellIdx_j++){
                
                bool ScintillatorFired = false;
                for (auto traj:*trajectories){
                    
                    // record only protons knocked out by primaries
                    G4int parentId = traj->GetParentID ();
                    G4int PDGcode = traj->GetPDGEncoding();
                    if (parentId!=1 || PDGcode!=2212) continue;
                    
                    
                    G4int trackId = traj->GetTrackID ();
                    G4double edep = fEdepScintillator[facetIdx][cellIdx_i][cellIdx_j][trackId];
                    
                    if (edep > ThresholdEdepScintillator) {
                        ScintillatorFired = true;
                    }
                } // end trajectories
                if (ScintillatorFired){
                    NScintillatorsFiredPerPrimary ++ ;
                }
            } // end cells j
        } // end cells i
    } // end facets
        
    
    
    if (fdebug>1) {
        std::cout << "for (auto traj:*trajectories)" << std::endl;
        std::cout << "output particles csv file " << particlesfilename << std::endl;
        std::cout << "output particles summary csv file " << primariesfilename << std::endl;
    }
    
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
            
            // for data saving, omit neutrinos
            if (
                ParticleName=="nu_e"
                )
                continue;
            
            //  trajectory points can not help, as they do not provide access for the energy deposition etc.
            // for each particle, we dedicate a line in the csv file
            particlescsvfile
            << eventId          << ","
            << trackId          << ","
            << parentId         << ","
            << ParticleName     << ","
            << p_init.mag()/MeV             << ","
            << ProcessName.at(trackId)      << ","
            << NScintillatorsFiredPerPrimary << ",";
            
            
            // track hit-position and energy deposition in scintillators
            for (int facetIdx=0; facetIdx < NFacets; facetIdx++){
                for (int cellIdx_i=0; cellIdx_i < NCells_i; cellIdx_i++){
                    for (int cellIdx_j=0; cellIdx_j < NCells_j; cellIdx_j++){
                        particlescsvfile << fEdepScintillator[facetIdx][cellIdx_i][cellIdx_j][trackId]/MeV      << ",";
                    }
                }
            }
            
            
            particlescsvfile
            << PDGcode                      << ","
            << p_init.x()/MeV               << ","
            << p_init.y()/MeV               << ","
            << p_init.z()/MeV               << ",";
            
            // end line
            particlescsvfile << std::endl;
            
        }
        particlescsvfile.close();
        if (fdebug>1) std::cout << "wrote to particlescsvfile " << particlesfilename << std::endl;
    }
    
    if (do_primaries_file){
        // output brief particles csv file
        G4int Nlines = 0;
        primariescsvfile.open(primariesfilename, std::ios_base::app);
        for (auto traj:*trajectories){
            // record only primaries
            G4int parentId = traj->GetParentID ();
            if (parentId!=0) continue;
            
            G4String ParticleName = traj->GetParticleName() ;
            G4int PDGcode = traj->GetPDGEncoding();
            G4int trackId = traj->GetTrackID ();
            G4ThreeVector p_init = traj->GetInitialMomentum ();
            
            
            
            //  trajectory points can not help, as they do not provide access for the energy deposition etc.
            // for each particle, we dedicate a line in the csv file
            primariescsvfile
            << eventId          << ","
            << trackId          << ","
            << ParticleName     << ","
            << p_init.mag()/MeV             << ","
            << p_init.theta()/deg           << ","
            << p_init.phi()/deg             << ","
            << ProcessName.at(trackId)      << ","
            << NScintillatorsFiredPerPrimary<< ","
            << PDGcode                      << ",";
            // end line
            primariescsvfile << std::endl;
            Nlines++;
        }
        primariescsvfile.close();
        if (fdebug>1) std::cout << "wrote " << Nlines
            << " lines to primariescsvfile " << std::endl
            << primariesfilename << std::endl;
    }
    
    if (do_events_file) {
//                // output events csv file
//                // copy from EventAction::EndOfEventAction(const G4Event*evt)
//                eventsscsvfile.open(eventsfilename, std::ios_base::app);
//
//                eventsscsvfile
//                << eventId          << ",";
//                for (int iVol=0; iVol<4; iVol++){
//                    eventsscsvfile
//                    << fEdepTotVol_e_g[iVol]/MeV    << "," // energy deposit by e+,e-,gamma
//                    << fEdepTotVol[iVol]/MeV        << ","; // total energy deposition in event (over all tracks)
//                }
//                eventsscsvfile << std::endl;
//                // close file
//                eventsscsvfile.close();
//                if (fdebug>1) std::cout << "wrote to eventsscsvfile " << eventsfilename << std::endl;
    }
    
    
    PrintAction("done processing");
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......




//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void EventAction::AddEdepScintillator(G4int facetIdx, G4int cellIdx_i, G4int cellIdx_j,
                                      G4double edep,
                                      G4double time,
                                      G4int trackId,
                                      G4int PDGcode){
    
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
    if (facetIdx >= NFacets || cellIdx_i>NCells_i || cellIdx_j>NCells_j){
        std::cout
        << "(facetIdx = " << facetIdx << ", cellIdx_i = " << cellIdx_i << ", cellIdx_j=" << cellIdx_j << ")"
        <<
        "returning from EventAction::AddEdep() " << std::endl;
        return;
    }
    // fEdep is a 2-D array
    // first dimension is volume number
    // (0 = source holder, 1=scintllator 1, 2 = scintillator 2, 3 = world)
    // second dimension is track Id
    fEdepScintillator[facetIdx][cellIdx_i][cellIdx_j][trackId] += edep;
//    fEdepTotVol[iVol] += edep;
    
//    if ((PDGcode==22) || (PDGcode==11)  || (PDGcode==-11) ){
//        fEdepTotVol_e_g[iVol] += edep;
//    }
    if (fdebug>1) {
        std::cout
        << "fEdepScintillator["<<facetIdx<<"]["<<cellIdx_i<<"]["<<cellIdx_j<<"]["<<trackId<<"] = "
        << fEdepScintillator[facetIdx][cellIdx_i][cellIdx_j][trackId]/keV << " keV"
        << std::endl;
        
        
        
        std::cout << "done EventAction::AddEdep()" << std::endl;
    }
    
}


