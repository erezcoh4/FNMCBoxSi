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
: G4VUserPrimaryGeneratorAction(),fParticleGun(0),fDetector(det), fdebug(_fdebug_)
{
    G4int n_particle = 1;
    fParticleGun  = new G4ParticleGun(n_particle);
    
    G4double Source_y = fDetector->GetSourceY();
    fParticleGun->SetParticleEnergy(0*eV);
    fParticleGun->SetParticlePosition(G4ThreeVector(0.*cm,Source_y,0.*cm));
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
    delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
    if (fParticleGun->GetParticleDefinition() == G4Geantino::Geantino()) {
        // 22Na
        G4int Z = 11, A = 22;
        G4double ionCharge   = 0.*eplus;
        G4double excitEnergy = 0.*keV;
        
        G4ParticleDefinition* ion
        = G4IonTable::GetIonTable()->GetIon(Z,A,excitEnergy);
        fParticleGun->SetParticleDefinition(ion);
        fParticleGun->SetParticleCharge(ionCharge);
    }
    // standard particle
    //    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    //    G4ParticleDefinition* particle
    //    = particleTable->FindParticle("e+");
    //    fParticleGun->SetParticleDefinition(particle);
    //    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
    //    fParticleGun->SetParticleEnergy(1.*keV);

    if (fdebug>1) std::cout << "generating primaries at PrimaryGeneratorAction::GeneratePrimaries()" << std::endl;
    //create vertex
    if (fdebug>1) G4Random::showEngineStatus();
    fParticleGun->GeneratePrimaryVertex(anEvent);
    if (fdebug>1) std::cout << "done PrimaryGeneratorAction::GeneratePrimaries()" << std::endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

