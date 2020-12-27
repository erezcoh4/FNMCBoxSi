// ********************************************************************
/// \file PrimaryGeneratorAction.cc
/// \brief Implementation of the PrimaryGeneratorAction class
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4Geantino.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction(DetectorConstruction* det, int _fdebug_)
: G4VUserPrimaryGeneratorAction(),fParticleGun(0),fDetector(det), fdebug(_fdebug_){
    G4int n_particle = 1;
    fParticleGun  = new G4ParticleGun(n_particle);
    fParticleGun->SetParticleEnergy(2.*MeV);
    fParticleGun->SetParticlePosition(G4ThreeVector(0.*cm,0.*cm,0.*cm));
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction() {
    delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent){
    //    if (fParticleGun->GetParticleDefinition() == G4Geantino::Geantino()) {
    //
    //        //        // 22Na
    //        //        G4int Z = 11, A = 22;
    //        // 252Cf
    //        //        G4int Z = 98, A = 252;
    ////        // ions
    ////        G4double ionCharge   = 0.*eplus;
    ////        G4double excitEnergy = 0.*keV;
    //
    //        G4ParticleDefinition* ion
    //        = G4IonTable::GetIonTable()->GetIon(Z,A,excitEnergy);
    //        fParticleGun->SetParticleDefinition(ion);
    //        fParticleGun->SetParticleCharge(ionCharge);
    //    }
    
    // isotropic emission of neutrons with 2 MeV energy
    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* particle = particleTable->FindParticle("neutron");
    fParticleGun->SetParticleDefinition(particle);
    
    // isotropic emission
    G4double theta = G4UniformRand()*CLHEP::pi;
    G4double phi = G4UniformRand()*CLHEP::twopi;
    G4double px = std::sin(theta)*std::cos(phi);
    G4double py = std::sin(theta)*std::sin(phi);
    G4double pz = std::cos(theta);
    
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(px,py,pz));
    // 2 MeV energy
    fParticleGun->SetParticleEnergy(2.*MeV);

    
    if (fdebug>1) std::cout << "generating primaries at PrimaryGeneratorAction::GeneratePrimaries()" << std::endl;
    //create vertex
    if (fdebug>1) G4Random::showEngineStatus();
    fParticleGun->GeneratePrimaryVertex(anEvent);
    if (fdebug>1) std::cout << "done PrimaryGeneratorAction::GeneratePrimaries()" << std::endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

