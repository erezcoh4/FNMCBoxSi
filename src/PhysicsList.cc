/// \file PhysicsList.cc
/// \brief Implementation of the PhysicsList class

#include "PhysicsList.hh"

#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4Radioactivation.hh"
#include "G4UAtomicDeexcitation.hh"
#include "G4LossTableManager.hh"
#include "G4NuclearLevelData.hh"
#include "G4DeexPrecoParameters.hh"
#include "G4NuclideTable.hh"
#include "G4IonConstructor.hh"
#include "G4PhysicsListHelper.hh"
#include "QBBC.hh"
#include "FTFP_BERT.hh"
//#include "HadronElasticPhysicsHP.hh"
#include "G4HadronElasticPhysics.hh"

#include "G4ParticleHPElastic.hh"
#include "G4ParticleHPElasticData.hh"
#include "G4ParticleHPThermalScattering.hh"
#include "G4ParticleHPThermalScatteringData.hh"

#include "G4EmStandardPhysics.hh"
#include "G4EmExtraPhysics.hh"
#include "G4EmParameters.hh"
#include "G4DecayPhysics.hh"
#include "G4NuclideTable.hh"
#include "BiasedRDPhysics.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4HadronPhysicsFTFP_BERT.hh"
#include "G4HadronInelasticQBBC.hh"
#include "G4HadronPhysicsINCLXX.hh"
#include "G4IonElasticPhysics.hh"
#include "G4IonPhysics.hh"
#include "G4IonINCLXXPhysics.hh"

// particles

#include "G4BosonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BosonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

#include "G4EmPenelopePhysics.hh"


PhysicsList::PhysicsList()
:G4VModularPhysicsList()
{
    G4int verb = 0;
    SetVerboseLevel(verb);
    
    //add new units for radioActive decays
    //
    new G4UnitDefinition( "millielectronVolt", "meV", "Energy", 1.e-3*eV);
    //
    const G4double minute = 60*second;
    const G4double hour   = 60*minute;
    const G4double day    = 24*hour;
    const G4double year   = 365*day;
    new G4UnitDefinition("minute", "min", "Time", minute);
    new G4UnitDefinition("hour",   "h",   "Time", hour);
    new G4UnitDefinition("day",    "d",   "Time", day);
    new G4UnitDefinition("year",   "y",   "Time", year);
    
    // Mandatory for G4NuclideTable
    // Half-life threshold must be set small or many short-lived isomers
    // will not be assigned life times (default to 0)
    G4NuclideTable::GetInstance()->SetThresholdOfHalfLife(0.1*picosecond);
    G4NuclideTable::GetInstance()->SetLevelTolerance(1.0*eV);
    
    //read new PhotonEvaporation data set
    G4DeexPrecoParameters* deex =
    G4NuclearLevelData::GetInstance()->GetParameters();
    deex->SetCorrelatedGamma(true);
    deex->SetStoreAllLevels(true);
    deex->SetMaxLifeTime(G4NuclideTable::GetInstance()->GetThresholdOfHalfLife()
                         /std::log(2.));
    
    
    
    // EM physics
    RegisterPhysics(new G4EmStandardPhysics());
    G4EmParameters* param = G4EmParameters::Instance();
    param->SetAugerCascade(true);
    param->SetStepFunction(1., 1*CLHEP::mm);
    param->SetStepFunctionMuHad(1., 1*CLHEP::mm);
    
    // Decay
    RegisterPhysics(new G4DecayPhysics());
    
    // Radioactive decay
    RegisterPhysics(new BiasedRDPhysics());
    
    // Hadron Elastic scattering
    RegisterPhysics( new G4HadronElasticPhysics(verb) );
    
    // Hadron Inelastic physics
    RegisterPhysics( new G4HadronPhysicsFTFP_BERT(verb));
    RegisterPhysics( new G4HadronInelasticQBBC(verb));
    RegisterPhysics( new G4HadronPhysicsINCLXX(verb));
    
    // Ion Elastic scattering
    RegisterPhysics( new G4IonElasticPhysics(verb));
    
    // Ion Inelastic physics
    RegisterPhysics( new G4IonPhysics(verb));
    RegisterPhysics( new G4IonINCLXXPhysics(verb));
    
    // for positron annihilation
    RegisterPhysics( new G4EmPenelopePhysics() );
        
    // Gamma-Nuclear Physics
    G4EmExtraPhysics* gnuc = new G4EmExtraPhysics(verb);
    gnuc->ElectroNuclear(false);
    gnuc->MuonNuclear(false);
    RegisterPhysics(gnuc);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::~PhysicsList()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructParticle()
{
    G4BosonConstructor  pBosonConstructor;
    pBosonConstructor.ConstructParticle();
    
    G4LeptonConstructor pLeptonConstructor;
    pLeptonConstructor.ConstructParticle();
    
    G4MesonConstructor pMesonConstructor;
    pMesonConstructor.ConstructParticle();
    
    G4BaryonConstructor pBaryonConstructor;
    pBaryonConstructor.ConstructParticle();
    
    G4IonConstructor pIonConstructor;
    pIonConstructor.ConstructParticle();
    
    G4ShortLivedConstructor pShortLivedConstructor;
    pShortLivedConstructor.ConstructParticle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetCuts()
{
    SetCutValue(0*mm, "proton");
    SetCutValue(0*mm, "neutron");
    SetCutValue(10*km, "e-");
    SetCutValue(10*km, "e+");
    SetCutValue(10*km, "gamma");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructProcess()
{
    AddTransportation();
    
    G4Radioactivation* radioactiveDecay = new G4Radioactivation();
    
    G4bool ARMflag = false;
    radioactiveDecay->SetARM(ARMflag);        //Atomic Rearangement
    
    // need to initialize atomic deexcitation
    //
    G4LossTableManager* man = G4LossTableManager::Instance();
    G4VAtomDeexcitation* deex = man->AtomDeexcitation();
    if (!deex) {
        ///G4EmParameters::Instance()->SetFluo(true);
        G4EmParameters::Instance()->SetAugerCascade(ARMflag);
        deex = new G4UAtomicDeexcitation();
        deex->InitialiseAtomicDeexcitation();
        man->SetAtomDeexcitation(deex);
    }
    
    // register radioactiveDecay
    //
    G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
    ph->RegisterProcess(radioactiveDecay, G4GenericIon::GenericIon());
    
    
    // penelope annihilation
    G4VPhysicsConstructor* emPhysicsList = new G4EmPenelopePhysics();
    emPhysicsList->ConstructProcess();
    


    
    // neutron interaction with plastics etc., from [examples/extended/hadronic/NeutronSource/src/HadronElasticPhysicsHP.cc]
    g4HadronElasticPhysics = new G4HadronElasticPhysics();
    g4HadronElasticPhysics -> ConstructProcess();
    g4HadronElasticPhysics -> GetNeutronModel() -> SetMinEnergy(100*keV);
    G4HadronicProcess* process = g4HadronElasticPhysics -> GetNeutronProcess();
    G4ParticleHPElastic* model1 = new G4ParticleHPElastic();
    process->RegisterMe(model1);
    process->AddDataSet(new G4ParticleHPElasticData());


    
    //printout
    //
    G4cout << "\n  Set atomic relaxation mode " << ARMflag << G4endl;
}
