#include "Options.h"
#include "changelog.h"

Options* Options::m_Instance = nullptr;

namespace po  = boost::program_options;
namespace cls = po::command_line_style;
namespace fs  = boost::filesystem;


// this is required to set default value for boost program options of type vector<std::string>
namespace std
{
  std::ostream& operator<<(std::ostream &os, const std::vector<std::string> &vec) {    
    for (auto item : vec) os << item << " ";
    return os; 
  }
} 


Options* Options::Instance() {
    if (!m_Instance) {
        m_Instance = new Options();
    }
    return m_Instance;
}

// for convenience

#define COMPLAIN(complainStr)           { std::stringstream _ss; _ss << complainStr; throw _ss.str(); }
#define COMPLAIN_IF(cond, complainStr)  { if (cond) COMPLAIN(complainStr) }


/*
 * Initialise option values
 * 
 * This function initialises the values of the options - this is where the defaults
 * are set.  If a user does not specify an option, either on the commandline or in
 * a grid file, the default values are taken from here.
 * 
 * Note this is a class OptionValues function.
 * 
 * 
 * void Options::OptionValues::Initialise()
 * 
 */
void Options::OptionValues::Initialise() {

    m_Populated = false;        

    // set all options to their default values

    // flags

    m_AllowRLOFAtBirth                                              = false;
    m_AllowTouchingAtBirth                                          = false;

    m_DebugToFile                                                   = false;
    m_ErrorsToFile                                                  = false;

    m_EnableWarnings                                                = false;

	m_BeBinaries                                                    = false;
    m_EvolvePulsars                                                 = false;
	m_EvolveUnboundSystems                                          = false;

    m_DetailedOutput                                                = false;
    m_PopulationDataPrinting                                        = false;
    m_PrintBoolAsString                                             = false;
    m_Quiet                                                         = false;
    m_RlofPrinting                                                  = false;
    m_SwitchLog                                                     = false;

    m_nBatchesUsed                                                  = -1;


    // Evolution mode: SSE or BSE
    m_EvolutionMode                                                 = EVOLUTION_MODE::BSE;
    m_EvolutionModeString                                           = EVOLUTION_MODE_LABEL.at(m_EvolutionMode);

    // Population synthesis variables
    m_ObjectsToEvolve                                               = 10;

    m_FixedRandomSeed                                               = false;                                                // TRUE if --random-seed is passed on command line
    m_RandomSeed                                                    = 0;

    // Specify how long to evolve for
    m_MaxEvolutionTime                                              = 13700.0;
    m_MaxNumberOfTimestepIterations                                 = 99999;
    m_TimestepMultiplier                                            = 1.0;

    // Initial mass options
    m_InitialMass                                                   = 5.0;
    m_InitialMass1                                                  = 5.0;
    m_InitialMass2                                                  = 5.0;

    m_InitialMassFunction                                           = INITIAL_MASS_FUNCTION::KROUPA;
    m_InitialMassFunctionString                                     = INITIAL_MASS_FUNCTION_LABEL.at(m_InitialMassFunction);
    m_InitialMassFunctionMin                                        = 8.0;
    m_InitialMassFunctionMax                                        = 100.0; 
    m_InitialMassFunctionPower                                      = -2.3;


    // Initial mass ratios
    m_MassRatioDistribution                                         = MASS_RATIO_DISTRIBUTION::FLAT;                        // Most likely want FLAT or SANA2012
    m_MassRatioDistributionString                                   = MASS_RATIO_DISTRIBUTION_LABEL.at(m_MassRatioDistribution);
    m_MassRatioDistributionMin                                      = 0.0;
    m_MassRatioDistributionMax                                      = 1.0;

    m_MinimumMassSecondary                                          = 0.0;


    // Initial orbit options
    m_SemiMajorAxis                                                 = 0.1;
    m_SemiMajorAxisDistribution                                     = SEMI_MAJOR_AXIS_DISTRIBUTION::FLATINLOG;              // Most likely want FLATINLOG or SANA2012
    m_SemiMajorAxisDistributionString                               = SEMI_MAJOR_AXIS_DISTRIBUTION_LABEL.at(SEMI_MAJOR_AXIS_DISTRIBUTION::FLATINLOG);
    m_SemiMajorAxisDistributionMin                                  = 0.1;
    m_SemiMajorAxisDistributionMax                                  = 1000.0;
    m_SemiMajorAxisDistributionPower                                = -1.0;

    // Initial orbital period
    m_PeriodDistributionMin                                         = 1.1;
    m_PeriodDistributionMax                                         = 1000.0;

    // Eccentricity
    m_Eccentricity                                                  = 0.0;
    m_EccentricityDistribution                                      = ECCENTRICITY_DISTRIBUTION::ZERO; 
    m_EccentricityDistributionString                                = ECCENTRICITY_DISTRIBUTION_LABEL.at(m_EccentricityDistribution);
    m_EccentricityDistributionMin                                   = 0.0;
    m_EccentricityDistributionMax                                   = 1.0;

    // Kick options
    m_KickMagnitudeDistribution                                     = KICK_MAGNITUDE_DISTRIBUTION::MAXWELLIAN;
    m_KickMagnitudeDistributionString                               = KICK_MAGNITUDE_DISTRIBUTION_LABEL.at(m_KickMagnitudeDistribution);
    m_KickMagnitudeDistributionSigmaCCSN_NS                         = 250;
    m_KickMagnitudeDistributionSigmaCCSN_BH                         = 250;
    m_KickMagnitudeDistributionMaximum                              = -1.0; 
    m_KickMagnitudeDistributionSigmaForECSN                         = 30.0;
    m_KickMagnitudeDistributionSigmaForUSSN   	                    = 30.0;
	m_KickScalingFactor						                        = 1.0;

    // Kick direction option
    m_KickDirectionDistribution                                     = KICK_DIRECTION_DISTRIBUTION::ISOTROPIC;
    m_KickDirectionDistributionString                               = KICK_DIRECTION_DISTRIBUTION_LABEL.at(m_KickDirectionDistribution);
    m_KickDirectionPower                                            = 0.0;

    // Kick magnitude
    m_KickMagnitude                                                 = 0.0;
    m_KickMagnitude1                                                = 0.0;
    m_KickMagnitude2                                                = 0.0;                               

    // Kick magnitude random number (used to draw kick magnitude if necessary)
    m_KickMagnitudeRandom                                           = RAND->Random();
    m_KickMagnitudeRandom1                                          = RAND->Random();
    m_KickMagnitudeRandom2                                          = RAND->Random();

    // Mean anomaly
    m_KickMeanAnomaly1                                              = RAND->Random(0.0, _2_PI);
    m_KickMeanAnomaly2                                              = RAND->Random(0.0, _2_PI);

    // Phi
    m_KickPhi1                                                      = 0.0;                                              // actual value set later
    m_KickPhi2                                                      = 0.0;                                              // actual value set later

    // Theta
    m_KickTheta1                                                    = 0.0;                                              // actual value set later 
    m_KickTheta2                                                    = 0.0;                                              // actual value set later

    // Black hole kicks
    m_BlackHoleKicksOption                                          = BLACK_HOLE_KICK_OPTION::FALLBACK;
    m_BlackHoleKicksOptionString                                    = BLACK_HOLE_KICK_OPTION_LABEL.at(m_BlackHoleKicksOption);


    // Chemically Homogeneous Evolution
    m_CheOption                                                     = CHE_OPTION::NONE;
    m_CheString                                                     = CHE_OPTION_LABEL.at(m_CheOption);


    // Supernova remnant mass prescription options
    m_RemnantMassPrescription                                       = REMNANT_MASS_PRESCRIPTION::FRYER2012;
    m_RemnantMassPrescriptionString                                 = REMNANT_MASS_PRESCRIPTION_LABEL.at(m_RemnantMassPrescription);

    m_FryerSupernovaEngine                                          = SN_ENGINE::DELAYED;
    m_FryerSupernovaEngineString                                    = SN_ENGINE_LABEL.at(m_FryerSupernovaEngine);

    m_NeutrinoMassLossAssumptionBH                                  = NEUTRINO_MASS_LOSS_PRESCRIPTION::FIXED_FRACTION;
    m_NeutrinoMassLossAssumptionBHString                            = NEUTRINO_MASS_LOSS_PRESCRIPTION_LABEL.at(m_NeutrinoMassLossAssumptionBH);
    m_NeutrinoMassLossValueBH                                       = 0.1;


    // Fixed uk options
    m_UseFixedUK                                                    = false;
    m_FixedUK                                                       = -1.0;


    // Pair instability and pulsational pair instability mass loss
    m_UsePairInstabilitySupernovae                                  = false;
    m_PairInstabilityLowerLimit                                     = 60.0;                                                 // Belczynski+ 2016 is 65 Msol
    m_PairInstabilityUpperLimit                                     = 135.0;                                                // Belczynski+ 2016 is 135 Msol

    m_UsePulsationalPairInstability                                 = false;
    m_PulsationalPairInstabilityLowerLimit                          = 35.0;                                                 // Belczynski+ 2016 is 45 Msol
    m_PulsationalPairInstabilityUpperLimit                          = 60.0;                                                 // Belczynski+ 2016 is 65 Msol

    m_PulsationalPairInstabilityPrescription                        = PPI_PRESCRIPTION::COMPAS;
    m_PulsationalPairInstabilityPrescriptionString                  = PPI_PRESCRIPTION_LABEL.at(m_PulsationalPairInstabilityPrescription);

	m_MaximumNeutronStarMass                                        = 3.0;                                                  // StarTrack is 3.0
    
    m_mCBUR1                                                        = MCBUR1HURLEY;                                         // MHurley value, Fryer+ and Belczynski+ use 1.83


    // Output path
    m_OutputPathString                                              = ".";
    m_DefaultOutputPath                                             = boost::filesystem::current_path();
    m_OutputPath                                                    = m_DefaultOutputPath;
    m_OutputContainerName                                           = DEFAULT_OUTPUT_CONTAINER_NAME;
    

    // Mass loss options
    m_UseMassLoss                                                   = false;

    m_MassLossPrescription                                          = MASS_LOSS_PRESCRIPTION::VINK;
    m_MassLossPrescriptionString                                    = MASS_LOSS_PRESCRIPTION_LABEL.at(m_MassLossPrescription);


    // Wind mass loss multiplicitive constants
    m_LuminousBlueVariableFactor                                    = 1.5;
    m_WolfRayetFactor                                               = 1.0;


    // Mass transfer options
    m_UseMassTransfer                                               = true;
	m_CirculariseBinaryDuringMassTransfer         	                = false;
	m_AngularMomentumConservationDuringCircularisation              = false;

    // Case BB/BC mass transfer stability prescription
    m_CaseBBStabilityPrescription                                   = CASE_BB_STABILITY_PRESCRIPTION::ALWAYS_STABLE;
    m_CaseBBStabilityPrescriptionString                             = CASE_BB_STABILITY_PRESCRIPTION_LABEL.at(m_CaseBBStabilityPrescription);

    // Options for mass transfer accretion efficiency
    m_MassTransferAccretionEfficiencyPrescription                   = MT_ACCRETION_EFFICIENCY_PRESCRIPTION::THERMALLY_LIMITED;
    m_MassTransferAccretionEfficiencyPrescriptionString             = MT_ACCRETION_EFFICIENCY_PRESCRIPTION_LABEL.at(m_MassTransferAccretionEfficiencyPrescription);

    m_MassTransferFractionAccreted                                  = 1.0;
    m_MassTransferCParameter                                        = 10.0;
    m_EddingtonAccretionFactor                                      = 1;                                                    // >1 is super-eddington, 0 is no accretion

    // Mass transfer thermally limited options
	m_MassTransferThermallyLimitedVariation                         = MT_THERMALLY_LIMITED_VARIATION::C_FACTOR;
	m_MassTransferThermallyLimitedVariationString                   = MT_THERMALLY_LIMITED_VARIATION_LABEL.at(m_MassTransferThermallyLimitedVariation);

    // Mass transfer angular momentum loss prescription options
    m_MassTransferJloss                                             = 1.0;
    m_MassTransferAngularMomentumLossPrescription                   = MT_ANGULAR_MOMENTUM_LOSS_PRESCRIPTION::ISOTROPIC_RE_EMISSION;
    m_MassTransferAngularMomentumLossPrescriptionString             = MT_ANGULAR_MOMENTUM_LOSS_PRESCRIPTION_LABEL.at(m_MassTransferAngularMomentumLossPrescription);

    // Mass transfer rejuvenation prescriptions
    m_MassTransferRejuvenationPrescription                          = MT_REJUVENATION_PRESCRIPTION::NONE;
    m_MassTransferRejuvenationPrescriptionString                    = MT_REJUVENATION_PRESCRIPTION_LABEL.at(m_MassTransferRejuvenationPrescription);

    // Mass transfer critical mass ratios
    m_MassTransferCriticalMassRatioMSLowMass                        = false;
    m_MassTransferCriticalMassRatioMSLowMassNonDegenerateAccretor   = 1.44;                                                 // Claeys+ 2014 = 1.44
    m_MassTransferCriticalMassRatioMSLowMassDegenerateAccretor      = 1.0;                                                  // Claeys+ 2014 = 1.0

    // AVG
    /*
    m_MassTransferCriticalMassRatioMSHighMass                       = false;
    m_MassTransferCriticalMassRatioMSHighMassNonDegenerateAccretor  = 0.625;                                                // Claeys+ 2014 = 0.625
    m_MassTransferCriticalMassRatioMSHighMassDegenerateAccretor     = 0.0;

    m_MassTransferCriticalMassRatioHG                               = false;
    m_MassTransferCriticalMassRatioHGNonDegenerateAccretor          = 0.40;                                                 // Claeys+ 2014 = 0.25
    m_MassTransferCriticalMassRatioHGDegenerateAccretor             = 0.21;                                                 // Claeys+ 2014 = 0.21

    m_MassTransferCriticalMassRatioGiant                            = false;
    m_MassTransferCriticalMassRatioGiantNonDegenerateAccretor       = 0.0;
    m_MassTransferCriticalMassRatioGiantDegenerateAccretor          = 0.87;                                                 // Claeys+ 2014 = 0.81

    m_MassTransferCriticalMassRatioHeliumMS                         = false;
    m_MassTransferCriticalMassRatioHeliumMSNonDegenerateAccretor    = 0.625;
    m_MassTransferCriticalMassRatioHeliumMSDegenerateAccretor       = 0.0;

    m_MassTransferCriticalMassRatioHeliumHG                         = false;
    m_MassTransferCriticalMassRatioHeliumHGNonDegenerateAccretor    = 0.25;                                                 // Claeys+ 2014 = 0.25
    m_MassTransferCriticalMassRatioHeliumHGDegenerateAccretor       = 0.21;                                                 // Claeys+ 2014 = 0.21

    m_MassTransferCriticalMassRatioHeliumGiant                      = false;
    m_MassTransferCriticalMassRatioHeliumGiantNonDegenerateAccretor = 1.28;                                                 // Claeys+ 2014 = 0.25
    m_MassTransferCriticalMassRatioHeliumGiantDegenerateAccretor    = 0.87;

    m_MassTransferCriticalMassRatioWhiteDwarf                       = false;
	m_MassTransferCriticalMassRatioWhiteDwarfNonDegenerateAccretor  = 0.0;
    m_MassTransferCriticalMassRatioWhiteDwarfDegenerateAccretor     = 1.6;                                                  // Claeys+ 2014 = 1.6
    */

    // Common Envelope options
    m_CommonEnvelopeAlpha                                           = 1.0;
    m_CommonEnvelopeLambda                                          = 0.1;
	m_CommonEnvelopeSlopeKruckow                                    = -4.0 / 5.0;
	m_CommonEnvelopeAlphaThermal                                    = 1.0;
    m_CommonEnvelopeLambdaMultiplier                                = 1.0;
    m_AllowMainSequenceStarToSurviveCommonEnvelope                  = false;

    // Prescription for envelope state (radiative or convective)
    m_EnvelopeStatePrescription                                     = ENVELOPE_STATE_PRESCRIPTION::LEGACY;
    m_EnvelopeStatePrescriptionString                               = ENVELOPE_STATE_PRESCRIPTION_LABEL.at(m_EnvelopeStatePrescription);

    // Accretion during common envelope
    m_CommonEnvelopeMassAccretionPrescription                       = CE_ACCRETION_PRESCRIPTION::ZERO;
    m_CommonEnvelopeMassAccretionPrescriptionString                 = CE_ACCRETION_PRESCRIPTION_LABEL.at(m_CommonEnvelopeMassAccretionPrescription);
    
    m_CommonEnvelopeMassAccretionMin                                = 0.04;
    m_CommonEnvelopeMassAccretionMax                                = 0.1;
    m_CommonEnvelopeMassAccretionConstant                           = 0.0;

	// Common envelope lambda prescription
	m_CommonEnvelopeLambdaPrescription                              = CE_LAMBDA_PRESCRIPTION::NANJING;
	m_CommonEnvelopeLambdaPrescriptionString                        = CE_LAMBDA_PRESCRIPTION_LABEL.at(m_CommonEnvelopeLambdaPrescription);

	// Common envelope Nandez and Ivanova energy formalism
	m_RevisedEnergyFormalismNandezIvanova	                        = false;
	m_MaximumMassDonorNandezIvanova                                 = 2.0;
	m_CommonEnvelopeRecombinationEnergyDensity                      = 1.5E13;


    // Adaptive Importance Sampling options
    m_AISexploratoryPhase                                           = false;
    m_AISDCOtype                                                    = AIS_DCO::ALL;
    m_AISDCOtypeString                                              = AIS_DCO_LABEL.at(AIS_DCO::ALL);
    m_AIShubble                                                     = false;
    m_AISpessimistic                                                = false;
    m_AISrefinementPhase                                            = false;
    m_AISrlof                                                       = false;
    m_KappaGaussians                                                = 2;


	// Zetas
	m_StellarZetaPrescription                                       = ZETA_PRESCRIPTION::SOBERMAN;
	m_StellarZetaPrescriptionString                                 = ZETA_PRESCRIPTION_LABEL.at(m_StellarZetaPrescription);
	m_ZetaAdiabaticArbitrary                                        = 10000.0;                                              // large value favours stable MT
    m_ZetaMainSequence 	                                            = 2.0;
	m_ZetaRadiativeEnvelopeGiant	                                = 6.5;


    // Metallicity options
    m_Metallicity                                                   = ZSOL;


    // Neutron star equation of state
    m_NeutronStarEquationOfState                                    = NS_EOS::SSE;
    m_NeutronStarEquationOfStateString                              = NS_EOSLabel.at(NS_EOS::SSE);


    // Pulsar birth magnetic field distribution
    m_PulsarBirthMagneticFieldDistribution                          = PULSAR_BIRTH_MAGNETIC_FIELD_DISTRIBUTION::ZERO;
    m_PulsarBirthMagneticFieldDistributionString                    = PULSAR_BIRTH_MAGNETIC_FIELD_DISTRIBUTION_LABEL.at(m_PulsarBirthMagneticFieldDistribution);
    m_PulsarBirthMagneticFieldDistributionMin                       = 11.0;
    m_PulsarBirthMagneticFieldDistributionMax                       = 13.0;


    // Pulsar birth spin period distribution string
    m_PulsarBirthSpinPeriodDistribution                             = PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION::ZERO;
    m_PulsarBirthSpinPeriodDistributionString                       = PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION_LABEL.at(m_PulsarBirthSpinPeriodDistribution);
    m_PulsarBirthSpinPeriodDistributionMin                          = 0.0;
    m_PulsarBirthSpinPeriodDistributionMax                          = 100.0;

    m_PulsarMagneticFieldDecayTimescale                             = 1000.0;
    m_PulsarMagneticFieldDecayMassscale                             = 0.025;
    m_PulsarLog10MinimumMagneticField                               = 8.0;


    // Rotational velocity distribution options
    m_RotationalVelocityDistribution                                = ROTATIONAL_VELOCITY_DISTRIBUTION::ZERO;
    m_RotationalVelocityDistributionString                          = ROTATIONAL_VELOCITY_DISTRIBUTION_LABEL.at(m_RotationalVelocityDistribution);


	// grids

	m_GridFilename                                                  = "";


    // debug and logging options

    m_DebugLevel                                                    = 0;
    m_DebugClasses.clear();

    m_LogLevel                                                      = 0;
    m_LogClasses.clear();


    // Logfiles    
    m_LogfileDefinitionsFilename                                    = "";
    m_LogfileDelimiter                                              = DELIMITER::TAB;
    m_LogfileDelimiterString                                        = DELIMITERLabel.at(m_LogfileDelimiter);
    m_LogfileNamePrefix                                             = "";

    m_LogfileSystemParameters                                       = get<0>(LOGFILE_DESCRIPTOR.at(LOGFILE::SYSTEM_PARAMETERS));
    m_LogfileDetailedOutput                                         = (m_EvolutionMode == EVOLUTION_MODE::SSE) ? get<0>(LOGFILE_DESCRIPTOR.at(LOGFILE::SSE_DETAILED_OUTPUT)) : get<0>(LOGFILE_DESCRIPTOR.at(LOGFILE::BSE_DETAILED_OUTPUT));
    m_LogfileDoubleCompactObjects                                   = get<0>(LOGFILE_DESCRIPTOR.at(LOGFILE::DOUBLE_COMPACT_OBJECTS));
    m_LogfileSupernovae                                             = (m_EvolutionMode == EVOLUTION_MODE::SSE) ? get<0>(LOGFILE_DESCRIPTOR.at(LOGFILE::SSE_SUPERNOVAE)) : get<0>(LOGFILE_DESCRIPTOR.at(LOGFILE::BSE_SUPERNOVAE));
    m_LogfileCommonEnvelopes                                        = get<0>(LOGFILE_DESCRIPTOR.at(LOGFILE::COMMON_ENVELOPES));
    m_LogfileRLOFParameters                                         = get<0>(LOGFILE_DESCRIPTOR.at(LOGFILE::RLOF_PARAMETERS));
    m_LogfileBeBinaries                                             = get<0>(LOGFILE_DESCRIPTOR.at(LOGFILE::BE_BINARIES));
    m_LogfilePulsarEvolution                                        = get<0>(LOGFILE_DESCRIPTOR.at(LOGFILE::BSE_PULSAR_EVOLUTION)); // only BSE for now
    m_LogfileSwitchLog                                              = (m_EvolutionMode == EVOLUTION_MODE::SSE) ? get<0>(LOGFILE_DESCRIPTOR.at(LOGFILE::SSE_SWITCH_LOG)) : get<0>(LOGFILE_DESCRIPTOR.at(LOGFILE::BSE_SWITCH_LOG));
    
}
    

/*
 * Sanity check options and option values
 * 
 * We can currently sample mass, metallicity, separation, eccentricity etc. within COMPAS,
 * so those options don't need to be mandatory - but when we move all sampling out of 
 * COMPAS we will need to enforce those as mandatory options - unless we decide to leave
 * some minimal sampling inside COMPAS to allow for missing options.  We would only need
 * to leave a single distribution for each - we wouldn't want to give the user the option
 * of choosing a distribution - the functionality would only be for convenience if an
 * option was missing.
 * 
 * The boost variable map from the parsed options should already have been set before calling
 * this function.  This records, for each option, whether the user specified a value and, if 
 * so, the value specified by the user.  This function sanity checks the user specified values, 
 * sets the values if all pass the sanity checks, then sets the values of the options not 
 * specified by the user the the defaults specifed here.
 * 
 * Note this is a class OptionValues function.
 * 
 * 
 * std::string Options::OptionValues::CheckAndSetOptions(const po::variables_map p_VM)
 * 
 * @return                                      String containing an error string
 *                                              If no error occurred the return string will be the empty string 
 */
std::string Options::OptionValues::CheckAndSetOptions() {
#define DEFAULTED(opt) m_VM[opt].defaulted()    // for convenience and readability - undefined at end of function

    std::string errStr = "";                                                                                // error string

    // check & set prescriptions, distributions, assumptions etc. options - alphabetically

    try {

        bool found;

        m_FixedRandomSeed  = !DEFAULTED("random-seed");                                                     // use random seed if it is provided by the user
        m_UseFixedUK       = !DEFAULTED("fix-dimensionless-kick-magnitude") && (m_FixedUK >= 0.0);          // determine if user supplied a valid kick magnitude


        // Floor
        /*
        if (!DEFAULTED("ais-dcotype")) {                                                                    // Adaptive Importance Sampling DCO type
            std::tie(found, p_OptionValues->m_AISDCOtype) = utils::GetMapKey(p_OptionValues->m_AISDCOtypeString, AIS_DCO_LABEL, p_OptionValues->m_AISDCOtype);
            return "Unknown AIS DCO Type";
        }
        */

        if (!DEFAULTED("black-hole-kicks")) {                                                               // black hole kicks option
            std::tie(found, m_BlackHoleKicksOption) = utils::GetMapKey(m_BlackHoleKicksOptionString, BLACK_HOLE_KICK_OPTION_LABEL, m_BlackHoleKicksOption);
            COMPLAIN_IF(!found, "Unknown Black Hole Kicks Option");
        }

        if (!DEFAULTED("case-bb-stability-prescription")) {                                                 //case BB/BC mass transfer stability prescription
            std::tie(found, m_CaseBBStabilityPrescription) = utils::GetMapKey(m_CaseBBStabilityPrescriptionString, CASE_BB_STABILITY_PRESCRIPTION_LABEL, m_CaseBBStabilityPrescription);
            COMPLAIN_IF(!found, "Unknown Case BB/BC Mass Transfer Stability Prescription");
        }
           
        if (!DEFAULTED("chemically-homogeneous-evolution")) {                                               // Chemically Homogeneous Evolution
            std::tie(found, m_CheOption) = utils::GetMapKey(m_CheString, CHE_OPTION_LABEL, m_CheOption);
            COMPLAIN_IF(!found, "Unknown Chemically Homogeneous Evolution Option");
        }

        if (!DEFAULTED("common-envelope-lambda-prescription")) {                                            // common envelope lambda prescription
            std::tie(found, m_CommonEnvelopeLambdaPrescription) = utils::GetMapKey(m_CommonEnvelopeLambdaPrescriptionString, CE_LAMBDA_PRESCRIPTION_LABEL, m_CommonEnvelopeLambdaPrescription);
            COMPLAIN_IF(!found, "Unknown CE Lambda Prescription");
        }

        if (!DEFAULTED("common-envelope-mass-accretion-prescription")) {                                    // common envelope mass accretion prescription
            std::tie(found, m_CommonEnvelopeMassAccretionPrescription) = utils::GetMapKey(m_CommonEnvelopeMassAccretionPrescriptionString, CE_ACCRETION_PRESCRIPTION_LABEL, m_CommonEnvelopeMassAccretionPrescription);
            COMPLAIN_IF(!found, "Unknown CE Mass Accretion Prescription");
        }
            
        if (!DEFAULTED("envelope-state-prescription")) {                                                    // envelope state prescription
            std::tie(found, m_EnvelopeStatePrescription) = utils::GetMapKey(m_EnvelopeStatePrescriptionString, ENVELOPE_STATE_PRESCRIPTION_LABEL, m_EnvelopeStatePrescription);
            COMPLAIN_IF(!found, "Unknown Envelope State Prescription");
        }

        if (!DEFAULTED("eccentricity-distribution")) {                                                      // eccentricity distribution
            std::tie(found, m_EccentricityDistribution) = utils::GetMapKey(m_EccentricityDistributionString, ECCENTRICITY_DISTRIBUTION_LABEL, m_EccentricityDistribution);
            COMPLAIN_IF(!found, "Unknown Eccentricity Distribution");
        }

        if (!DEFAULTED("fryer-supernova-engine")) {                                                         // Fryer et al. 2012 supernova engine
            std::tie(found, m_FryerSupernovaEngine) = utils::GetMapKey(m_FryerSupernovaEngineString, SN_ENGINE_LABEL, m_FryerSupernovaEngine);
            COMPLAIN_IF(!found, "Unknown Fryer et al. Supernova Engine");
        }

        if (!DEFAULTED("initial-mass-function")) {                                                          // initial mass function
            std::tie(found, m_InitialMassFunction) = utils::GetMapKey(m_InitialMassFunctionString, INITIAL_MASS_FUNCTION_LABEL, m_InitialMassFunction);
            COMPLAIN_IF(!found, "Unknown Initial Mass Function");
        }

        if (!DEFAULTED("kick-direction")) {                                                                 // kick direction
            std::tie(found, m_KickDirectionDistribution) = utils::GetMapKey(m_KickDirectionDistributionString, KICK_DIRECTION_DISTRIBUTION_LABEL, m_KickDirectionDistribution);
            COMPLAIN_IF(!found, "Unknown Kick Direction Distribution");
        }

        if (!DEFAULTED("kick-magnitude-distribution")) {                                                    // kick magnitude
            std::tie(found, m_KickMagnitudeDistribution) = utils::GetMapKey(m_KickMagnitudeDistributionString, KICK_MAGNITUDE_DISTRIBUTION_LABEL, m_KickMagnitudeDistribution);
            COMPLAIN_IF(!found, "Unknown Kick Magnitude Distribution");
        }

        // set values for m_KickPhi[1/2] and m_KickTheta[1/2] here
        // we now have the kick direction distribution and kick direction power (exponent) required by the user (either default or specified)

        bool phi1Defaulted   = DEFAULTED("kick-phi-1");
        bool theta1Defaulted = DEFAULTED("kick-theta-1");

        if (phi1Defaulted || theta1Defaulted) {
            double phi1, theta1;
            std::tie(phi1, theta1) = utils::DrawKickDirection(m_KickDirectionDistribution, m_KickDirectionPower);
            if (phi1Defaulted  ) m_KickPhi1   = phi1;
            if (theta1Defaulted) m_KickTheta1 = theta1;
        }

        bool phi2Defaulted   = DEFAULTED("kick-phi-2");
        bool theta2Defaulted = DEFAULTED("kick-theta-2");

        if (phi2Defaulted || theta2Defaulted) {
            double phi2, theta2;
            std::tie(phi2, theta2) = utils::DrawKickDirection(m_KickDirectionDistribution, m_KickDirectionPower);
            if (phi2Defaulted  ) m_KickPhi2   = phi2;
            if (theta2Defaulted) m_KickTheta2 = theta2;
        }

        if (!DEFAULTED("logfile-delimiter")) {                                                              // logfile field delimiter
            std::tie(found, m_LogfileDelimiter) = utils::GetMapKey(m_LogfileDelimiterString, DELIMITERLabel, m_LogfileDelimiter);
            COMPLAIN_IF(!found, "Unknown Logfile Delimiter");
        }

        if (!DEFAULTED("mass-loss-prescription")) {                                                         // mass loss prescription
            std::tie(found, m_MassLossPrescription) = utils::GetMapKey(m_MassLossPrescriptionString, MASS_LOSS_PRESCRIPTION_LABEL, m_MassLossPrescription);
            COMPLAIN_IF(!found, "Unknown Mass Loss Prescription");
        }

        if (!DEFAULTED("mass-ratio-distribution")) {                                                        // mass ratio distribution
            std::tie(found, m_MassRatioDistribution) = utils::GetMapKey(m_MassRatioDistributionString, MASS_RATIO_DISTRIBUTION_LABEL, m_MassRatioDistribution);
            COMPLAIN_IF(!found, "Unknown Mass Ratio Distribution");
        }

        if (m_UseMassTransfer && !DEFAULTED("mass-transfer-accretion-efficiency-prescription")) {           // mass transfer accretion efficiency prescription
            std::tie(found, m_MassTransferAccretionEfficiencyPrescription) = utils::GetMapKey(m_MassTransferAccretionEfficiencyPrescriptionString, MT_ACCRETION_EFFICIENCY_PRESCRIPTION_LABEL, m_MassTransferAccretionEfficiencyPrescription);
            COMPLAIN_IF(!found, "Unknown Mass Transfer Angular Momentum Loss Prescription");
        }

        if (m_UseMassTransfer && !DEFAULTED("mass-transfer-angular-momentum-loss-prescription")) {          // mass transfer angular momentum loss prescription
            std::tie(found, m_MassTransferAngularMomentumLossPrescription) = utils::GetMapKey(m_MassTransferAngularMomentumLossPrescriptionString, MT_ANGULAR_MOMENTUM_LOSS_PRESCRIPTION_LABEL, m_MassTransferAngularMomentumLossPrescription);
            COMPLAIN_IF(!found, "Unknown Mass Transfer Angular Momentum Loss Prescription");
        }

        if (m_UseMassTransfer && !DEFAULTED("mass-transfer-rejuvenation-prescription")) {                   // mass transfer rejuvenation prescription
            std::tie(found, m_MassTransferRejuvenationPrescription) = utils::GetMapKey(m_MassTransferRejuvenationPrescriptionString, MT_REJUVENATION_PRESCRIPTION_LABEL, m_MassTransferRejuvenationPrescription);
            COMPLAIN_IF(!found, "Unknown Mass Transfer Rejuvenation Prescription");
        }

        if (m_UseMassTransfer && !DEFAULTED("mass-transfer-thermal-limit-accretor")) {                      // mass transfer accretor thermal limit
            std::tie(found, m_MassTransferThermallyLimitedVariation) = utils::GetMapKey(m_MassTransferThermallyLimitedVariationString, MT_THERMALLY_LIMITED_VARIATION_LABEL, m_MassTransferThermallyLimitedVariation);
            COMPLAIN_IF(!found, "Unknown Mass Transfer Accretor Thermal Limit");

            if (m_MassTransferThermallyLimitedVariation == MT_THERMALLY_LIMITED_VARIATION::C_FACTOR) {
                m_MassTransferCParameter = DEFAULTED("mass-transfer-thermal-limit-C") ? 10.0 : m_MassTransferCParameter;
            }

            if (m_MassTransferThermallyLimitedVariation == MT_THERMALLY_LIMITED_VARIATION::RADIUS_TO_ROCHELOBE) {
                m_MassTransferCParameter = DEFAULTED("mass-transfer-thermal-limit-C") ? 1.0 : m_MassTransferCParameter;
            }
        }

        if (!DEFAULTED("mode")) {                                                                           // mode
            std::tie(found, m_EvolutionMode) = utils::GetMapKey(m_EvolutionModeString, EVOLUTION_MODE_LABEL, m_EvolutionMode);
            COMPLAIN_IF(!found, "Unknown Mode");
        }

        if (!DEFAULTED("neutrino-mass-loss-bh-formation")) {                                                // neutrino mass loss assumption
            std::tie(found, m_NeutrinoMassLossAssumptionBH) = utils::GetMapKey(m_NeutrinoMassLossAssumptionBHString, NEUTRINO_MASS_LOSS_PRESCRIPTION_LABEL, m_NeutrinoMassLossAssumptionBH);
            COMPLAIN_IF(!found, "Unknown Neutrino Mass Loss Assumption");
        }

        if (!DEFAULTED("neutron-star-equation-of-state")) {                                                 // neutron star equation of state
            std::tie(found, m_NeutronStarEquationOfState) = utils::GetMapKey(m_NeutronStarEquationOfStateString, NS_EOSLabel, m_NeutronStarEquationOfState);
            COMPLAIN_IF(!found, "Unknown Neutron Star Equation of State");
        }

        if (!DEFAULTED("pulsar-birth-magnetic-field-distribution")) {                                       // pulsar birth magnetic field distribution
            std::tie(found, m_PulsarBirthMagneticFieldDistribution) = utils::GetMapKey(m_PulsarBirthMagneticFieldDistributionString, PULSAR_BIRTH_MAGNETIC_FIELD_DISTRIBUTION_LABEL, m_PulsarBirthMagneticFieldDistribution);
            COMPLAIN_IF(!found, "Unknown Pulsar Birth Magnetic Field Distribution");
        }

        if (!DEFAULTED("pulsar-birth-spin-period-distribution")) {                                          // pulsar birth spin period distribution
            std::tie(found, m_PulsarBirthSpinPeriodDistribution) = utils::GetMapKey(m_PulsarBirthSpinPeriodDistributionString, PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION_LABEL, m_PulsarBirthSpinPeriodDistribution);
            COMPLAIN_IF(!found, "Unknown Pulsar Birth Spin Period Distribution");
        }

        if (!DEFAULTED("pulsational-pair-instability-prescription")) {                                      // pulsational pair instability prescription
            std::tie(found, m_PulsationalPairInstabilityPrescription) = utils::GetMapKey(m_PulsationalPairInstabilityPrescriptionString, PPI_PRESCRIPTION_LABEL, m_PulsationalPairInstabilityPrescription);
            COMPLAIN_IF(!found, "Unknown Pulsational Pair Instability Prescription");
        }

        if (!DEFAULTED("remnant-mass-prescription")) {                                                      // remnant mass prescription
            std::tie(found, m_RemnantMassPrescription) = utils::GetMapKey(m_RemnantMassPrescriptionString, REMNANT_MASS_PRESCRIPTION_LABEL, m_RemnantMassPrescription);
            COMPLAIN_IF(!found, "Unknown Remnant Mass Prescription");
        }

        if (!DEFAULTED("rotational-velocity-distribution")) {                                               // rotational velocity distribution
            std::tie(found, m_RotationalVelocityDistribution) = utils::GetMapKey(m_RotationalVelocityDistributionString, ROTATIONAL_VELOCITY_DISTRIBUTION_LABEL, m_RotationalVelocityDistribution);
            COMPLAIN_IF(!found, "Unknown Rotational Velocity Distribution");
        }

        if (!DEFAULTED("semi-major-axis-distribution")) {                                                   // semi-major axis distribution
            std::tie(found, m_SemiMajorAxisDistribution) = utils::GetMapKey(m_SemiMajorAxisDistributionString, SEMI_MAJOR_AXIS_DISTRIBUTION_LABEL, m_SemiMajorAxisDistribution);
            COMPLAIN_IF(!found, "Unknown Semi-Major Axis Distribution");
        }

        if (!DEFAULTED("stellar-zeta-prescription")) {                                                      // common envelope zeta prescription
            std::tie(found, m_StellarZetaPrescription) = utils::GetMapKey(m_StellarZetaPrescriptionString, ZETA_PRESCRIPTION_LABEL, m_StellarZetaPrescription);
            COMPLAIN_IF(!found, "Unknown stellar Zeta Prescription");
        }

        // constraint/value/range checks - alphabetically (where possible)

        COMPLAIN_IF(!DEFAULTED("common-envelope-alpha") && m_CommonEnvelopeAlpha < 0.0, "CE alpha (--common-envelope-alpha) < 0");
        COMPLAIN_IF(!DEFAULTED("common-envelope-alpha-thermal") && (m_CommonEnvelopeAlphaThermal < 0.0 || m_CommonEnvelopeAlphaThermal > 1.0), "CE alpha thermal (--common-envelope-alpha-thermal) must be between 0 and 1");
        COMPLAIN_IF(!DEFAULTED("common-envelope-lambda-multiplier") && m_CommonEnvelopeLambdaMultiplier < 0.0, "CE lambda multiplie (--common-envelope-lambda-multiplier < 0");
        COMPLAIN_IF(!DEFAULTED("common-envelope-mass-accretion-constant") && m_CommonEnvelopeMassAccretionConstant < 0.0, "CE mass accretion constant (--common-envelope-mass-accretion-constant) < 0");
        COMPLAIN_IF(!DEFAULTED("common-envelope-mass-accretion-max") && m_CommonEnvelopeMassAccretionMax < 0.0, "Maximum accreted mass (--common-envelope-mass-accretion-max) < 0");
        COMPLAIN_IF(!DEFAULTED("common-envelope-mass-accretion-min") && m_CommonEnvelopeMassAccretionMin < 0.0, "Minimum accreted mass (--common-envelope-mass-accretion-min) < 0");

        COMPLAIN_IF(m_DebugLevel < 0, "Debug level (--debug-level) < 0");

        COMPLAIN_IF(m_Eccentricity < 0.0 || m_Eccentricity > 1.0, "Eccentricity (--eccentricity) must be between 0 and 1");
        COMPLAIN_IF(m_EccentricityDistributionMin < 0.0 || m_EccentricityDistributionMin > 1.0, "Minimum eccentricity (--eccentricity-min) must be between 0 and 1");
        COMPLAIN_IF(m_EccentricityDistributionMax < 0.0 || m_EccentricityDistributionMax > 1.0, "Maximum eccentricity (--eccentricity-max) must be between 0 and 1");
        COMPLAIN_IF(m_EccentricityDistributionMax <= m_EccentricityDistributionMin, "Maximum eccentricity (--eccentricity-max) must be > Minimum eccentricity (--eccentricity-min)");

        COMPLAIN_IF(m_InitialMassFunctionMin < 0.0, "Minimum initial mass (--initial-mass-min) < 0");
        COMPLAIN_IF(m_InitialMassFunctionMax < 0.0, "Maximum initial mass (--initial-mass-max) < 0");
        COMPLAIN_IF(m_InitialMassFunctionMax <= m_InitialMassFunctionMin, "Maximum initial mass (--initial-mass-max) must be > Minimum initial mass (--initial-mass-min)");

        if (m_KickMagnitudeDistribution == KICK_MAGNITUDE_DISTRIBUTION::FLAT) {
            COMPLAIN_IF(m_KickMagnitudeDistributionMaximum <= 0.0, "User specified --kick-magnitude-distribution = FLAT with Maximum kick magnitude (--kick-magnitude-max) <= 0.0");
        }

        COMPLAIN_IF(m_LogLevel < 0, "Logging level (--log-level) < 0");
 
        COMPLAIN_IF(m_LuminousBlueVariableFactor < 0.0, "LBV multiplier (--luminous-blue-variable-multiplier) < 0");

        COMPLAIN_IF(m_MassRatioDistributionMin < 0.0 || m_MassRatioDistributionMin > 1.0, "Minimum mass ratio (--mass-ratio-min) must be between 0 and 1");
        COMPLAIN_IF(m_MassRatioDistributionMax < 0.0 || m_MassRatioDistributionMax > 1.0, "Maximum mass ratio (--mass-ratio-max) must be between 0 and 1");
        COMPLAIN_IF(m_MassRatioDistributionMax <= m_MassRatioDistributionMin, "Maximum mass ratio (--mass-ratio-max) must be > Minimum mass ratio (--mass-ratio-min)");

        COMPLAIN_IF(m_MaxEvolutionTime <= 0.0, "Maximum evolution time in Myr (--maxEvolutionTime) must be > 0");

        COMPLAIN_IF(m_Metallicity < 0.0 || m_Metallicity > 1.0, "Metallicity (--metallicity) should be absolute metallicity and must be between 0 and 1");

        COMPLAIN_IF(m_MinimumMassSecondary < 0.0, "Seconday minimum mass (--minimum-secondary-mass) must be >= 0");
        COMPLAIN_IF(m_MinimumMassSecondary > m_InitialMassFunctionMax, "Seconday minimum mass (--minimum-secondary-mass) must be <= Maximum initial mass (--initial-mass-max)");

        if (m_NeutrinoMassLossAssumptionBH == NEUTRINO_MASS_LOSS_PRESCRIPTION::FIXED_MASS) {
            COMPLAIN_IF(m_NeutrinoMassLossValueBH < 0.0, "Neutrino mass loss value < 0");
        }

        COMPLAIN_IF(m_ObjectsToEvolve <= 0, (m_EvolutionMode == EVOLUTION_MODE::SSE ? "Number of stars requested <= 0" : "Number of binaries requested <= 0"));
    
        if (m_NeutrinoMassLossAssumptionBH == NEUTRINO_MASS_LOSS_PRESCRIPTION::FIXED_FRACTION) {
            COMPLAIN_IF(m_NeutrinoMassLossValueBH < 0.0 || m_NeutrinoMassLossValueBH > 1.0, "Neutrino mass loss must be between 0 and 1");
        }

        if (!DEFAULTED("outputPath")) {                                                                     // user specified output path?
                                                                                                            // yes
            fs::path userPath = m_OutputPathString;                                                         // user-specifed path
            if (fs::is_directory(userPath)) {                                                               // valid directory?
                m_OutputPath = userPath;                                                                    // yes - set outputPath to user-specified path
            }
            else {                                                                                          // not a valid directory
                m_OutputPath = m_DefaultOutputPath;                                                         // use default path = CWD
            }
        }

        COMPLAIN_IF(m_PeriodDistributionMin < 0.0, "Minimum orbital period (--orbital-period-min) < 0");
        COMPLAIN_IF(m_PeriodDistributionMax < 0.0, "Maximum orbital period (--orbital-period-max) < 0");

        COMPLAIN_IF(!DEFAULTED("pulsar-magnetic-field-decay-timescale") && m_PulsarMagneticFieldDecayTimescale <= 0.0, "Pulsar magnetic field decay timescale (--pulsar-magnetic-field-decay-timescale) <= 0");
        COMPLAIN_IF(!DEFAULTED("pulsar-magnetic-field-decay-massscale") && m_PulsarMagneticFieldDecayMassscale <= 0.0, "Pulsar Magnetic field decay massscale (--pulsar-magnetic-field-decay-massscale) <= 0");

        COMPLAIN_IF(m_SemiMajorAxisDistributionMin < 0.0, "Minimum semi-major Axis (--semi-major-axis-min) < 0");
        COMPLAIN_IF(m_SemiMajorAxisDistributionMax < 0.0, "Maximum semi-major Axis (--semi-major-axis-max) < 0");

        COMPLAIN_IF(m_TimestepMultiplier <= 0.0, "Timestep multiplier (--timestep-multiplier) <= 0");

        COMPLAIN_IF(m_WolfRayetFactor < 0.0, "WR multiplier (--wolf-rayet-multiplier) < 0");

    }
    catch (po::error& e) {                                                                                  // program options exception
        errStr = e.what();                                                                                  // set the error string
    }
    catch (const std::string eStr) {                                                                        // custom exception
        errStr = eStr;                                                                                      // set the error string
    }
    catch (...) {                                                                                           // unhandled exception
        errStr = ERR_MSG(ERROR::UNHANDLED_EXCEPTION);                                                       // set the error string
    }

    return errStr;
#undef DEFAULTED
}


/*
 * Determine if the user specified a value for the option
 * 
 * 
 * int OptionSpecified(std::string p_OptionString) 
 * 
 * 
 * @param   [IN]    p_OptionString              String containing option name
 * @return                                      Int result:
 *                                                  -1: invalid/unknown option name
 *                                                   0: option was not specified by user - default value used
 *                                                   1: option specified by user - user specified value used
 */
int Options::OptionSpecified(std::string p_OptionString) {

    int  result = -1;                                                                                                               // default = invalid/unknown option
    
    try {
        if (m_GridLine.optionValues.m_VM.count(p_OptionString) > 0) {                                                               // option exists at grid line (evolving object) level?
            result = (!m_GridLine.optionValues.m_Populated || m_GridLine.optionValues.m_VM[p_OptionString].defaulted()) ? 0 : 1;    // yes - set result
        }

        if (result == 0) {                                                                                                          // do we already have a result?
            if (m_CmdLine.optionValues.m_VM.count(p_OptionString) > 0) {                                                            // no - option exists at the commandline (program) level?
                result = (!m_CmdLine.optionValues.m_Populated || m_CmdLine.optionValues.m_VM[p_OptionString].defaulted()) ? 0 : 1;  // yes - set result
            }
            else {                                                                                                                  // option does not exist at the commandline (program-level)
                result = -1;                                                                                                        // set result
            }
        }    
    }
    catch (po::error& e) {                                                                                                          // program options exception
        result = -1;
    }
    catch (...) {                                                                                                                   // unhandled exception
        result = -1;                                                                                                                // set return value - invalid/unknown option
    }
    
    return result;
}


/*
 * This function populates the Boost options_description onject with all defined
 * options.  The options strings are asscoated with the variable to be populated
 * for each option, and the default values are defined.
 * 
 * JR: todo: One day we should construct the list of options shown in the help
 *           string from the maps in constants.h rather than have them here as 
 *           literal strings - too much opportunity for then to get out of sync 
 *           doing it this way (and more work to update the strings here every 
 *           time something changes)
 * 
 * Note that both parameters are modified here.
 * 
 * 
 * bool AddOptions(OptionValues *p_Options, po::options_description *p_OptionsDescription)
 * 
 * @param   [IN/OUT]    p_Options                   Object containing option values
 * @param   [IN/OUT]    p_OptionsDescription        options_description onject
 * @return                                          Boolean result:
 *                                                      true  if options added ok
 *                                                      false if an error occurred
 */
bool Options::AddOptions(OptionValues *p_Options, po::options_description *p_OptionsDescription) {

    bool ok = true;                             // status - unless a problem occurs

    // create default strings for vector<std::string> types (too hard to do inline)

    std::ostringstream ss;

    // debug classes
    std::string defaultDebugClasses;
    ss << "";
    for (auto debugClass = p_Options->m_DebugClasses.begin(); debugClass != p_Options->m_DebugClasses.end(); ++debugClass) ss << *debugClass << ",";
    defaultDebugClasses = ss.str();
    if (defaultDebugClasses.size() > 0) defaultDebugClasses.erase(defaultDebugClasses.size() - 1);

    // log classes
    std::string defaultLogClasses;
    ss << "";
    for (auto logClass = p_Options->m_LogClasses.begin(); logClass != p_Options->m_LogClasses.end(); ++logClass) ss << *logClass << ",";
    defaultLogClasses = ss.str();
    if (defaultLogClasses.size() > 0) defaultLogClasses.erase(defaultLogClasses.size() - 1);


    // add options

    try {

        p_OptionsDescription->add_options()     // begin the list of options to be added - boost syntactic sugar

        // there is no good way of formatting these - the boost syntax doesn't help that much
        // there is just a boatload of options, so this function is just going to be long...


        // switches

        (
            "help,h",                                                      
            po::bool_switch(), "Print this help message"
        )
        (
            "version,v",                                                   
            po::bool_switch(), "Print COMPAS version string"
        )


        // boolean options - alphabetically

        // Floor
        /*
        (
            "ais-exploratory-phase",                                       
            po::value<bool>(&p_Options->m_AISexploratoryPhase)->default_value(p_Options->m_AISexploratoryPhase)->implicit_value(true),                                                            
            ("Run exploratory phase of STROOPWAFEL (default = " + std::string(p_Options->m_AISexploratoryPhase ? "TRUE" : "FALSE") + ")").c_str()
        )
        (
            "ais-hubble",                                                  
            po::value<bool>(&p_Options->m_AIShubble)->default_value(p_Options->m_AIShubble)->implicit_value(true),                                                                                
            ("Excluding not in Hubble time mergers selection in exploratory phase of STROOPWAFEL (default = " + std::string(p_Options->m_AIShubble ? "TRUE" : "FALSE") + ")").c_str()
        )
        (
            "ais-pessimistic",                                             
            po::value<bool>(&p_Options->m_AISpessimistic)->default_value(p_Options->m_AISpessimistic)->implicit_value(true),                                                                      
            ("Optimistic or Pessimistic selection in exploratory phase of STROOPWAFEL (default = " + std::string(p_Options->m_AISpessimistic ? "TRUE" : "FALSE") + ")").c_str()
        )
        (
            "ais-refinement-phase",                                        
            po::value<bool>(&p_Options->m_AISrefinementPhase)->default_value(p_Options->m_AISrefinementPhase)->implicit_value(true),                                                              
            ("Run main sampling phase (step2) of STROOPWAFEL (default = " + std::string(p_Options->m_AISrefinementPhase ? "TRUE" : "FALSE") + ")").c_str()
        )
        (
            "ais-rlof",                                                    
            po::value<bool>(&p_Options->m_AISrlof)->default_value(p_Options->m_AISrlof)->implicit_value(true),                                                                                    
            ("RLOFSecondaryZAMS selection in exploratory phase of STROOPWAFEL (default = " + std::string(p_Options->m_AISrlof ? "TRUE" : "FALSE") + ")").c_str()
       )
        */

        (
            "allow-rlof-at-birth",                                         
            po::value<bool>(&p_Options->m_AllowRLOFAtBirth)->default_value(p_Options->m_AllowRLOFAtBirth)->implicit_value(true),                                                                  
            ("Allow binaries that have one or both stars in RLOF at birth to evolve (default = " + std::string(p_Options->m_AllowRLOFAtBirth ? "TRUE" : "FALSE") + ")").c_str()
        )
        (
            "allow-touching-at-birth",                                     
            po::value<bool>(&p_Options->m_AllowTouchingAtBirth)->default_value(p_Options->m_AllowTouchingAtBirth)->implicit_value(true),                                                          
            ("Allow binaries that are touching at birth to evolve (default = " + std::string(p_Options->m_AllowTouchingAtBirth ? "TRUE" : "FALSE") + ")").c_str()
        )
        (
            "angular-momentum-conservation-during-circularisation",            
            po::value<bool>(&p_Options->m_AngularMomentumConservationDuringCircularisation)->default_value(p_Options->m_AngularMomentumConservationDuringCircularisation)->implicit_value(true),  
            ("Conserve angular momentum when binary is circularised when entering a Mass Transfer episode (default = " + std::string(p_Options->m_AngularMomentumConservationDuringCircularisation ? "TRUE" : "FALSE") + ")").c_str()
        )

        // Serena
        /* 
        (
            "be-binaries",                                                  
            po::value<bool>(&p_Options->m_BeBinaries)->default_value(p_Options->m_BeBinaries)->implicit_value(true),                                                                              
            ("Enable Be Binaries study (default = " + std::string(p_Options->m_BeBinaries ? "TRUE" : "FALSE") + ")").c_str()
        )
        */

        (
            "circularise-binary-during-mass-transfer",                         
            po::value<bool>(&p_Options->m_CirculariseBinaryDuringMassTransfer)->default_value(p_Options->m_CirculariseBinaryDuringMassTransfer)->implicit_value(true),                            
            ("Circularise binary when it enters a Mass Transfer episode (default = " + std::string(p_Options->m_CirculariseBinaryDuringMassTransfer ? "TRUE" : "FALSE") + ")").c_str()
        )
        (
            "common-envelope-allow-main-sequence-survive",                 
            po::value<bool>(&p_Options->m_AllowMainSequenceStarToSurviveCommonEnvelope)->default_value(p_Options->m_AllowMainSequenceStarToSurviveCommonEnvelope)->implicit_value(true),          
            ("Allow main sequence stars to survive common envelope evolution (default = " + std::string(p_Options->m_AllowMainSequenceStarToSurviveCommonEnvelope ? "TRUE" : "FALSE") + ")").c_str()
        )
        (
            "debug-to-file",                                               
            po::value<bool>(&p_Options->m_DebugToFile)->default_value(p_Options->m_DebugToFile)->implicit_value(true),                                                                            
            ("Write debug statements to file (default = " + std::string(p_Options->m_DebugToFile ? "TRUE" : "FALSE") + ")").c_str()
        )
        (
            "detailed-output",                                              
            po::value<bool>(&p_Options->m_DetailedOutput)->default_value(p_Options->m_DetailedOutput)->implicit_value(true),                                                                      
            ("Print detailed output to file (default = " + std::string(p_Options->m_DetailedOutput ? "TRUE" : "FALSE") + ")").c_str()
        )
        (
            "enable-warnings",                                             
            po::value<bool>(&p_Options->m_EnableWarnings)->default_value(p_Options->m_EnableWarnings)->implicit_value(true),                                                                      
            ("Display warning messages to stdout (default = " + std::string(p_Options->m_EnableWarnings ? "TRUE" : "FALSE") + ")").c_str()
        )
        (
            "errors-to-file",                                              
            po::value<bool>(&p_Options->m_ErrorsToFile)->default_value(p_Options->m_ErrorsToFile)->implicit_value(true),                                                                          
            ("Write error messages to file (default = " + std::string(p_Options->m_ErrorsToFile ? "TRUE" : "FALSE") + ")").c_str()
        )
        (
            "evolve-pulsars",                                              
            po::value<bool>(&p_Options->m_EvolvePulsars)->default_value(p_Options->m_EvolvePulsars)->implicit_value(true),                                                                        
            ("Evolve pulsars (default = " + std::string(p_Options->m_EvolvePulsars ? "TRUE" : "FALSE") + ")").c_str()
        )
        (
            "evolve-unbound-systems",                                      
            po::value<bool>(&p_Options->m_EvolveUnboundSystems)->default_value(p_Options->m_EvolveUnboundSystems)->implicit_value(true),                                                          
            ("Continue evolving stars even if the binary is disrupted (default = " + std::string(p_Options->m_EvolveUnboundSystems ? "TRUE" : "FALSE") + ")").c_str()
        )
        (
            "mass-transfer",                                                
            po::value<bool>(&p_Options->m_UseMassTransfer)->default_value(p_Options->m_UseMassTransfer)->implicit_value(true),                                                                    
            ("Enable mass transfer (default = " + std::string(p_Options->m_UseMassTransfer ? "TRUE" : "FALSE") + ")").c_str()
        )
        (
            "pair-instability-supernovae",                                 
            po::value<bool>(&p_Options->m_UsePairInstabilitySupernovae)->default_value(p_Options->m_UsePairInstabilitySupernovae)->implicit_value(true),                                          
            ("Enable pair instability supernovae (PISN) (default = " + std::string(p_Options->m_UsePairInstabilitySupernovae ? "TRUE" : "FALSE") + ")").c_str()
        )
        (
            "population-data-printing",                                      
            po::value<bool>(&p_Options->m_PopulationDataPrinting)->default_value(p_Options->m_PopulationDataPrinting)->implicit_value(true),                                                      
            ("Print details of population (default = " + std::string(p_Options->m_PopulationDataPrinting ? "TRUE" : "FALSE") + ")").c_str()
        )
        (
            "print-bool-as-string",                                        
            po::value<bool>(&p_Options->m_PrintBoolAsString)->default_value(p_Options->m_PrintBoolAsString)->implicit_value(true),                                                                
            ("Print boolean properties as 'TRUE' or 'FALSE' (default = " + std::string(p_Options->m_PrintBoolAsString ? "TRUE" : "FALSE") + ")").c_str()
        )
        (
            "pulsational-pair-instability",                                
            po::value<bool>(&p_Options->m_UsePulsationalPairInstability)->default_value(p_Options->m_UsePulsationalPairInstability)->implicit_value(true),                                        
            ("Enable mass loss due to pulsational-pair-instability (PPI) (default = " + std::string(p_Options->m_UsePulsationalPairInstability ? "TRUE" : "FALSE") + ")").c_str()
        )
        (
            "quiet",                                                       
            po::value<bool>(&p_Options->m_Quiet)->default_value(p_Options->m_Quiet)->implicit_value(true),                                                                                        
            ("Suppress printing (default = " + std::string(p_Options->m_Quiet ? "TRUE" : "FALSE") + ")").c_str()
        )
        (
            "revised-energy-formalism-nandez-ivanova",                     
            po::value<bool>(&p_Options->m_RevisedEnergyFormalismNandezIvanova)->default_value(p_Options->m_RevisedEnergyFormalismNandezIvanova)->implicit_value(true),                            
            ("Enable revised energy formalism (default = " + std::string(p_Options->m_RevisedEnergyFormalismNandezIvanova ? "TRUE" : "FALSE") + ")").c_str()
        )
        (
            "rlof-printing",                                                
            po::value<bool>(&p_Options->m_RlofPrinting)->default_value(p_Options->m_RlofPrinting)->implicit_value(true),                                                                          
            ("Enable output parameters before/after RLOF (default = " + std::string(p_Options->m_RlofPrinting ? "TRUE" : "FALSE") + ")").c_str()
        )
        (
            "switchlog",                                                
            po::value<bool>(&p_Options->m_SwitchLog)->default_value(p_Options->m_SwitchLog)->implicit_value(true),                                                                          
            ("Print switch log to file (default = " + std::string(p_Options->m_SwitchLog ? "TRUE" : "FALSE") + ")").c_str()
        )
        (
            "use-mass-loss",                                               
            po::value<bool>(&p_Options->m_UseMassLoss)->default_value(p_Options->m_UseMassLoss)->implicit_value(true),                                                                            
            ("Enable mass loss (default = " + std::string(p_Options->m_UseMassLoss ? "TRUE" : "FALSE") + ")").c_str()
        )


        // numerical options - alphabetically grouped by type 


        // unsigned long

        (
            "random-seed",                                                 
            po::value<unsigned long>(&p_Options->m_RandomSeed)->default_value(p_Options->m_RandomSeed),                                                                                           
            ("Random seed to use (default = " + std::to_string(p_Options->m_RandomSeed) + ")").c_str()
        )


        // int

        (
            "debug-level",                                                 
            po::value<int>(&p_Options->m_DebugLevel)->default_value(p_Options->m_DebugLevel),                                                                                                     
            ("Determines which print statements are displayed for debugging (default = " + std::to_string(p_Options->m_DebugLevel) + ")").c_str()
        )
        (
            "log-level",                                                   
            po::value<int>(&p_Options->m_LogLevel)->default_value(p_Options->m_LogLevel),                                                                                                         
            ("Determines which print statements are included in the logfile (default = " + std::to_string(p_Options->m_LogLevel) + ")").c_str()
        )
        (
            "maximum-number-timestep-iterations",                          
            po::value<int>(&p_Options->m_MaxNumberOfTimestepIterations)->default_value(p_Options->m_MaxNumberOfTimestepIterations),                                                               
            ("Maximum number of timesteps to evolve binary before giving up (default = " + std::to_string(p_Options->m_MaxNumberOfTimestepIterations) + ")").c_str()
        )

        // Floor
        /*
        (
            "nbatches-used",                                               
            po::value<int>(&p_Options->m_nBatchesUsed)->default_value(p_Options->m_nBatchesUsed),                                                                                                 
            ("Number of batches used, for STROOPWAFEL (AIS), -1 = not required (default = " + std::to_string(p_Options->m_nBatchesUsed) + ")").c_str()
        )
        */

        (
            "number-of-stars",                                        
            po::value<int>(&p_Options->m_ObjectsToEvolve)->default_value(p_Options->m_ObjectsToEvolve),                                                                                                       
            ("Specify the number of stars to simulate (SSE) (default = " + std::to_string(p_Options->m_ObjectsToEvolve) + ")").c_str()
        )


        // double

        (
            "common-envelope-alpha",                                       
            po::value<double>(&p_Options->m_CommonEnvelopeAlpha)->default_value(p_Options->m_CommonEnvelopeAlpha),                                                                                
            ("Common Envelope efficiency alpha (default = " + std::to_string(p_Options->m_CommonEnvelopeAlpha) + ")").c_str()
        )
        (
            "common-envelope-alpha-thermal",                               
            po::value<double>(&p_Options->m_CommonEnvelopeAlphaThermal)->default_value(p_Options->m_CommonEnvelopeAlphaThermal),                                                                  
            ("Defined such that lambda = alpha_th * lambda_b + (1.0 - alpha_th) * lambda_g (default = " + std::to_string(p_Options->m_CommonEnvelopeAlphaThermal) + ")").c_str()
        )
        (
            "common-envelope-lambda",                                      
            po::value<double>(&p_Options->m_CommonEnvelopeLambda)->default_value(p_Options->m_CommonEnvelopeLambda),                                                                              
            ("Common Envelope lambda (default = " + std::to_string(p_Options->m_CommonEnvelopeLambda) + ")").c_str()
        )
        (
            "common-envelope-lambda-multiplier",                           
            po::value<double>(&p_Options->m_CommonEnvelopeLambdaMultiplier)->default_value(p_Options->m_CommonEnvelopeLambdaMultiplier),                                                          
            ("Multiply lambda by some constant (default = " + std::to_string(p_Options->m_CommonEnvelopeLambdaMultiplier) + ")").c_str()
        )
        (
            "common-envelope-mass-accretion-constant",                     
            po::value<double>(&p_Options->m_CommonEnvelopeMassAccretionConstant)->default_value(p_Options->m_CommonEnvelopeMassAccretionConstant),                                                
            ("Value of mass accreted by NS/BH during common envelope evolution if assuming all NS/BH accrete same amount of mass (common-envelope-mass-accretion-prescription CONSTANT). Ignored otherwise (default = " + std::to_string(p_Options->m_CommonEnvelopeMassAccretionConstant) + ")").c_str()
        )
        (
            "common-envelope-mass-accretion-max",                          
            po::value<double>(&p_Options->m_CommonEnvelopeMassAccretionMax)->default_value(p_Options->m_CommonEnvelopeMassAccretionMax),                                                          
            ("Maximum amount of mass accreted by NS/BHs during common envelope evolution in solar masses (default = " + std::to_string(p_Options->m_CommonEnvelopeMassAccretionMax) + ")").c_str()
        )
        (
            "common-envelope-mass-accretion-min",                          
            po::value<double>(&p_Options->m_CommonEnvelopeMassAccretionMin)->default_value(p_Options->m_CommonEnvelopeMassAccretionMin),                                                          
            ("Minimum amount of mass accreted by NS/BHs during common envelope evolution in solar masses (default = " + std::to_string(p_Options->m_CommonEnvelopeMassAccretionMin) + ")").c_str()
        )
        (
            "common-envelope-recombination-energy-density",                
            po::value<double>(&p_Options->m_CommonEnvelopeRecombinationEnergyDensity)->default_value(p_Options->m_CommonEnvelopeRecombinationEnergyDensity),                                      
            ("Recombination energy density in erg/g (default = " + std::to_string(p_Options->m_CommonEnvelopeRecombinationEnergyDensity) + ")").c_str()
        )
        (
            "common-envelope-slope-kruckow",                               
            po::value<double>(&p_Options->m_CommonEnvelopeSlopeKruckow)->default_value(p_Options->m_CommonEnvelopeSlopeKruckow),                                                                  
            ("Common Envelope slope for Kruckow lambda (default = " + std::to_string(p_Options->m_CommonEnvelopeSlopeKruckow) + ")").c_str()
        )

        // AVG - 17/03/2020 - Uncomment mass-ratio options when fully implemented
        /*
        (
            "critical-mass-ratio-giant-degenerate-accretor",
            po::value<double>(&m_MassTransferCriticalMassRatioGiantDegenerateAccretor)->default_value(m_MassTransferCriticalMassRatioGiantDegenerateAccretor),
            ("Critical mass ratio for MT from a giant star (default = " + std::to_string(m_MassTransferCriticalMassRatioGiantDegenerateAccretor) + ") Specify both giant flags to use. 0 is always stable, <0 is disabled").c_str()
        )
        (
            "critical-mass-ratio-giant-non-degenerate-accretor",
            po::value<double>(&m_MassTransferCriticalMassRatioGiantNonDegenerateAccretor)->default_value(m_MassTransferCriticalMassRatioGiantNonDegenerateAccretor),
            ("Critical mass ratio for MT from a giant star (default = " + std::to_string(m_MassTransferCriticalMassRatioGiantNonDegenerateAccretor) + ") Specify both giant flags to use. 0 is always stable, <0 is disabled").c_str()
        )
        (
            "critical-mass-ratio-helium-giant-degenerate-accretor",
            po::value<double>(&m_MassTransferCriticalMassRatioHeliumGiantDegenerateAccretor)->default_value(m_MassTransferCriticalMassRatioHeliumGiantDegenerateAccretor),
            ("Critical mass ratio for MT from a helium giant star (default = " + std::to_string(m_MassTransferCriticalMassRatioHeliumGiantDegenerateAccretor) + ") Specify both helium giant flags to use. 0 is always stable, <0 is disabled").c_str()
        )
        (
            "critical-mass-ratio-helium-giant-non-degenerate-accretor",
            po::value<double>(&m_MassTransferCriticalMassRatioHeliumGiantNonDegenerateAccretor)->default_value(m_MassTransferCriticalMassRatioHeliumGiantNonDegenerateAccretor),
            ("Critical mass ratio for MT from a helium giant star (default = " + std::to_string(m_MassTransferCriticalMassRatioHeliumGiantNonDegenerateAccretor) + ") Specify both helium giant flags to use. 0 is always stable, <0 is disabled").c_str()
        )
        (
            "critical-mass-ratio-helium-hg-degenerate-accretor",
            po::value<double>(&m_MassTransferCriticalMassRatioHeliumHGDegenerateAccretor)->default_value(m_MassTransferCriticalMassRatioHeliumHGDegenerateAccretor),
            ("Critical mass ratio for MT from a helium HG star (default = " + std::to_string(m_MassTransferCriticalMassRatioHeliumHGDegenerateAccretor) + ") Specify both helium HG flags to use. 0 is always stable, <0 is disabled").c_str()
        )
        (
            "critical-mass-ratio-helium-hg-non-degenerate-accretor",
            po::value<double>(&m_MassTransferCriticalMassRatioHeliumHGNonDegenerateAccretor)->default_value(m_MassTransferCriticalMassRatioHeliumHGNonDegenerateAccretor),
            ("Critical mass ratio for MT from a helium HG star (default = " + std::to_string(m_MassTransferCriticalMassRatioHeliumHGNonDegenerateAccretor) + ") Specify both helium HG flags to use. 0 is always stable, <0 is disabled").c_str()
        )
        (
            "critical-mass-ratio-helium-ms-degenerate-accretor",
            po::value<double>(&m_MassTransferCriticalMassRatioHeliumMSDegenerateAccretor)->default_value(m_MassTransferCriticalMassRatioHeliumMSDegenerateAccretor),
            ("Critical mass ratio for MT from a helium MS star (default = " + std::to_string(m_MassTransferCriticalMassRatioHeliumMSDegenerateAccretor) + ") Specify both helium MS flags to use. 0 is always stable, <0 is disabled").c_str()
        )
        (
            "critical-mass-ratio-helium-ms-non-degenerate-accretor",
            po::value<double>(&m_MassTransferCriticalMassRatioHeliumMSNonDegenerateAccretor)->default_value(m_MassTransferCriticalMassRatioHeliumMSNonDegenerateAccretor),
            ("Critical mass ratio for MT from a helium MS star (default = " + std::to_string(m_MassTransferCriticalMassRatioHeliumMSNonDegenerateAccretor) + ") Specify both helium MS flags to use. 0 is always stable, <0 is disabled").c_str()
        )
        (
            "critical-mass-ratio-hg-degenerate-accretor",
            po::value<double>(&m_MassTransferCriticalMassRatioHGDegenerateAccretor)->default_value(m_MassTransferCriticalMassRatioHGDegenerateAccretor),
            ("Critical mass ratio for MT from a HG star (default = " + std::to_string(m_MassTransferCriticalMassRatioHGDegenerateAccretor) + ") Specify both HG flags to use. 0 is always stable, <0 is disabled").c_str()
        )
        (
            "critical-mass-ratio-hg-non-degenerate-accretor",
            po::value<double>(&m_MassTransferCriticalMassRatioHGNonDegenerateAccretor)->default_value(m_MassTransferCriticalMassRatioHGNonDegenerateAccretor),
            ("Critical mass ratio for MT from a HG star (default = " + std::to_string(m_MassTransferCriticalMassRatioHGNonDegenerateAccretor) + ") Specify both HG flags to use. 0 is always stable, <0 is disabled").c_str()
        )
        (
            "critical-mass-ratio-ms-high-mass-degenerate-accretor",
            po::value<double>(&m_MassTransferCriticalMassRatioMSHighMassDegenerateAccretor)->default_value(m_MassTransferCriticalMassRatioMSHighMassDegenerateAccretor),
            ("Critical mass ratio for MT from a MS star to a degenerate accretor (default = " + std::to_string(m_MassTransferCriticalMassRatioMSHighMassDegenerateAccretor) + " Specify both MS high mass flags to use. 0 is always stable, <0 is disabled").c_str()
        )
        (
            "critical-mass-ratio-ms-high-mass-non-degenerate-accretor",
            po::value<double>(&m_MassTransferCriticalMassRatioMSHighMassNonDegenerateAccretor)->default_value(m_MassTransferCriticalMassRatioMSHighMassNonDegenerateAccretor),
            ("Critical mass ratio for MT from a MS star (default = " + std::to_string(m_MassTransferCriticalMassRatioMSHighMassNonDegenerateAccretor) + ") Specify both MS high mass flags to use. 0 is always stable, <0 is disabled").c_str()
        )
        (
            "critical-mass-ratio-ms-low-mass-degenerate-accretor",
            po::value<double>(&m_MassTransferCriticalMassRatioMSLowMassDegenerateAccretor)->default_value(m_MassTransferCriticalMassRatioMSLowMassDegenerateAccretor),
            ("Critical mass ratio for MT from a MS star to a degenerate accretor (default = " + std::to_string(m_MassTransferCriticalMassRatioMSLowMassDegenerateAccretor) + " Specify both MS low mass flags to use. 0 is always stable, <0 is disabled").c_str()
        )
        (
            "critical-mass-ratio-ms-low-mass-non-degenerate-accretor",
            po::value<double>(&m_MassTransferCriticalMassRatioMSLowMassNonDegenerateAccretor)->default_value(m_MassTransferCriticalMassRatioMSLowMassNonDegenerateAccretor),
            ("Critical mass ratio for MT from a MS star (default = " + std::to_string(m_MassTransferCriticalMassRatioMSLowMassNonDegenerateAccretor) + ") Specify both MS low mass flags to use. 0 is always stable, <0 is disabled").c_str()
        )
        (
            "critical-mass-ratio-white-dwarf-degenerate-accretor",
            po::value<double>(&m_MassTransferCriticalMassRatioWhiteDwarfDegenerateAccretor)->default_value(m_MassTransferCriticalMassRatioWhiteDwarfDegenerateAccretor),
            ("Critical mass ratio for MT from a white dwarf (default = " + std::to_string(m_MassTransferCriticalMassRatioWhiteDwarfDegenerateAccretor) + ") Specify both white dwarf flags to use. 0 is always stable, <0 is disabled").c_str()
        )
        (
            "critical-mass-ratio-white-dwarf-non-degenerate-accretor",
            po::value<double>(&m_MassTransferCriticalMassRatioWhiteDwarfNonDegenerateAccretor)->default_value(m_MassTransferCriticalMassRatioWhiteDwarfNonDegenerateAccretor),
            ("Critical mass ratio for MT from a white dwarf (default = " + std::to_string(m_MassTransferCriticalMassRatioWhiteDwarfNonDegenerateAccretor) + ") Specify both white dwarf flags to use. 0 is always stable, <0 is disabled").c_str()
        )
        */

        (
            "eccentricity,e",                                            
            po::value<double>(&p_Options->m_Eccentricity)->default_value(p_Options->m_Eccentricity),                                                                
            ("Eccentricity, e (default = " + std::to_string(p_Options->m_Eccentricity) + ")").c_str()
        )
        (
            "eccentricity-max",                                            
            po::value<double>(&p_Options->m_EccentricityDistributionMax)->default_value(p_Options->m_EccentricityDistributionMax),                                                                
            ("Maximum eccentricity to generate (default = " + std::to_string(p_Options->m_EccentricityDistributionMax) + ")").c_str()
        )
        (
            "eccentricity-min",                                            
            po::value<double>(&p_Options->m_EccentricityDistributionMin)->default_value(p_Options->m_EccentricityDistributionMin),                                                                
            ("Minimum eccentricity to generate (default = " + std::to_string(p_Options->m_EccentricityDistributionMin) + ")").c_str()
        )
        (
            "eddington-accretion-factor",                                  
            po::value<double>(&p_Options->m_EddingtonAccretionFactor)->default_value(p_Options->m_EddingtonAccretionFactor),                                                                      
            ("Multiplication factor for eddington accretion for NS & BH, i.e. >1 is super-eddington and 0. is no accretion (default = " + std::to_string(p_Options->m_EddingtonAccretionFactor) + ")").c_str()
        )

        (
            "fix-dimensionless-kick-magnitude",                            
            po::value<double>(&p_Options->m_FixedUK)->default_value(p_Options->m_FixedUK),                                                                                                        
            ("Fix dimensionless kick magnitude uk to this value (default = " + std::to_string(p_Options->m_FixedUK) + ", -ve values false, +ve values true)").c_str()
        )

        (
            "initial-mass",                                            
            po::value<double>(&p_Options->m_InitialMass)->default_value(p_Options->m_InitialMass),                                                                          
            ("Initial mass (in Msol) for the star (SSE) (default = " + std::to_string(p_Options->m_InitialMass) + ")").c_str()
        )
        (
            "initial-mass-1",                                            
            po::value<double>(&p_Options->m_InitialMass1)->default_value(p_Options->m_InitialMass1),                                                                          
            ("Initial mass (in Msol) for the primary star (BSE) (default = " + std::to_string(p_Options->m_InitialMass1) + ")").c_str()
        )
        (
            "initial-mass-2",                                            
            po::value<double>(&p_Options->m_InitialMass2)->default_value(p_Options->m_InitialMass2),
            ("Initial mass (in Msol) for the secondary star (BSE) (default = " + std::to_string(p_Options->m_InitialMass2) + ")").c_str()
        )
        (
            "initial-mass-max",                                            
            po::value<double>(&p_Options->m_InitialMassFunctionMax)->default_value(p_Options->m_InitialMassFunctionMax),                                                                          
            ("Maximum mass (in Msol) to generate using given IMF (default = " + std::to_string(p_Options->m_InitialMassFunctionMax) + ")").c_str()
        )
        (
            "initial-mass-min",                                            
            po::value<double>(&p_Options->m_InitialMassFunctionMin)->default_value(p_Options->m_InitialMassFunctionMin),                                                                          
            ("Minimum mass (in Msol) to generate using given IMF (default = " + std::to_string(p_Options->m_InitialMassFunctionMin) + ")").c_str()
        )
        (
            "initial-mass-power",                                          
            po::value<double>(&p_Options->m_InitialMassFunctionPower)->default_value(p_Options->m_InitialMassFunctionPower),                                                                      
            ("Single power law power to generate primary mass using given IMF (default = " + std::to_string(p_Options->m_InitialMassFunctionPower) + ")").c_str()
        )

        // Floor
        /*
        (
            "kappa-gaussians",                                             
            po::value<double>(&p_Options->m_KappaGaussians)->default_value(p_Options->m_KappaGaussians),                                                                                          
            ("Scaling factor for the width of the Gaussian distributions in STROOPWAFEL main sampling phase (default = " + std::to_string(p_Options->m_KappaGaussians) + ")").c_str()
        )
        */

        (
            "kick-direction-power",                                        
            po::value<double>(&p_Options->m_KickDirectionPower)->default_value(p_Options->m_KickDirectionPower),                                                                                  
            ("Power for power law kick direction distribution (default = " + std::to_string(p_Options->m_KickDirectionPower) + " = isotropic, +ve = polar, -ve = in plane)").c_str()
        )
        (
            "kick-magnitude-max",                                          
            po::value<double>(&p_Options->m_KickMagnitudeDistributionMaximum)->default_value(p_Options->m_KickMagnitudeDistributionMaximum),                                                      
            ("Maximum drawn kick magnitude in km s^-1. Ignored if < 0. Must be > 0 if using kick-magnitude-distribution=FLAT (default = " + std::to_string(p_Options->m_KickMagnitudeDistributionMaximum) + ")").c_str()
        )
        (
            "kick-magnitude",                                          
            po::value<double>(&p_Options->m_KickMagnitude)->default_value(p_Options->m_KickMagnitude),                                                      
            ("The magnitude of the kick velocity the star receives during the a supernova (default = " + std::to_string(p_Options->m_KickMagnitude) + " km s^-1 )").c_str()
        )
        (
            "kick-magnitude-1",                                          
            po::value<double>(&p_Options->m_KickMagnitude1)->default_value(p_Options->m_KickMagnitude1),                                                      
            ("The magnitude of the kick velocity the primary star receives during the a supernova (default = " + std::to_string(p_Options->m_KickMagnitude1) + " km s^-1 )").c_str()
        )
        (
            "kick-magnitude-2",                                          
            po::value<double>(&p_Options->m_KickMagnitude2)->default_value(p_Options->m_KickMagnitude2),                                                      
            ("The magnitude of the kick velocity the secondary star receives during the a supernova (default = " + std::to_string(p_Options->m_KickMagnitude2) + " km s^-1 )").c_str()
        )
        (
            "kick-magnitude-random",                                          
            po::value<double>(&p_Options->m_KickMagnitudeRandom)->default_value(p_Options->m_KickMagnitudeRandom),                                                      
            "Number used to choose the kick velocity magnitude for the star during the a supernova (default = uniform random number [0.0, 1.0))"
        )
        (
            "kick-magnitude-random-1",                                          
            po::value<double>(&p_Options->m_KickMagnitudeRandom1)->default_value(p_Options->m_KickMagnitudeRandom1),                                                      
            "Number used to choose the kick velocity magnitude for the primary star during the a supernova (default = uniform random number [0.0, 1.0))"
        )
        (
            "kick-magnitude-random-2",                                          
            po::value<double>(&p_Options->m_KickMagnitudeRandom2)->default_value(p_Options->m_KickMagnitudeRandom2),                                                      
            "Number used to choose the kick velocity magnitude for the secondary during the a supernova (default = uniform random number [0.0, 1.0))"
        )
        (
            "kick-magnitude-sigma-ccsn-bh",                                
            po::value<double>(&p_Options->m_KickMagnitudeDistributionSigmaCCSN_BH)->default_value(p_Options->m_KickMagnitudeDistributionSigmaCCSN_BH),                                            
            ("Sigma for chosen kick magnitude distribution for black holes (default = " + std::to_string(p_Options->m_KickMagnitudeDistributionSigmaCCSN_BH) + " km s^-1 )").c_str()
        )
        (
            "kick-magnitude-sigma-ccsn-ns",                                
            po::value<double>(&p_Options->m_KickMagnitudeDistributionSigmaCCSN_NS)->default_value(p_Options->m_KickMagnitudeDistributionSigmaCCSN_NS),                                            
            ("Sigma for chosen kick magnitude distribution for neutron stars (default = " + std::to_string(p_Options->m_KickMagnitudeDistributionSigmaCCSN_NS) + " km s^-1 )").c_str()
        )
        (
            "kick-magnitude-sigma-ecsn",                                   
            po::value<double>(&p_Options->m_KickMagnitudeDistributionSigmaForECSN)->default_value(p_Options->m_KickMagnitudeDistributionSigmaForECSN),                                            
            ("Sigma for chosen kick magnitude distribution for ECSN (default = " + std::to_string(p_Options->m_KickMagnitudeDistributionSigmaForECSN) + " km s^-1 )").c_str()
        )
        (
            "kick-magnitude-sigma-ussn",                                   
            po::value<double>(&p_Options->m_KickMagnitudeDistributionSigmaForUSSN)->default_value(p_Options->m_KickMagnitudeDistributionSigmaForUSSN),                                            
            ("Sigma for chosen kick magnitude distribution for USSN (default = " + std::to_string(p_Options->m_KickMagnitudeDistributionSigmaForUSSN) + " km s^-1 )").c_str()
        )
        (
            "kick-mean-anomaly-1",
            po::value<double>(&p_Options->m_KickMeanAnomaly1)->default_value(p_Options->m_KickMeanAnomaly1),                                                                                  
            "Mean anomaly for the primary star at instantaneous time of the supernova (default = uniform random number [0.0, 2pi))"
        )
        (
            "kick-mean-anomaly-2",
            po::value<double>(&p_Options->m_KickMeanAnomaly2)->default_value(p_Options->m_KickMeanAnomaly2),                                                                                  
            "Mean anomaly for the secondary star at instantaneous time of the supernova (default = uniform random number [0.0, 2pi))"
        )
        (
            "kick-phi-1",
            po::value<double>(&p_Options->m_KickPhi1)->default_value(p_Options->m_KickPhi1),                                                                                  
            "Angle between 'x' and 'y', both in the orbital plane of the supernovae vector, for the primary star (default = drawn from kick direction distribution)"
        )
        (
            "kick-phi-2",
            po::value<double>(&p_Options->m_KickPhi2)->default_value(p_Options->m_KickPhi2),                                                                                  
            "Angle between 'x' and 'y', both in the orbital plane of the supernovae vector, for the secondary star (default = drawn from kick direction distribution)"
        )
        (
            "kick-scaling-factor",                                         
            po::value<double>(&p_Options->m_KickScalingFactor)->default_value(p_Options->m_KickScalingFactor),                                                                                    
            ("Arbitrary factor used to scale kicks (default = " + std::to_string(p_Options->m_KickScalingFactor) + ")").c_str()
        )
        (
            "kick-theta-1",                                        
            po::value<double>(&p_Options->m_KickTheta1)->default_value(p_Options->m_KickTheta1),                                                                                  
            "Angle between the orbital plane and the 'z' axis of the supernovae vector, for the primary star (default = drawn from kick direction distribution)"
        )
        (
            "kick-theta-2",                                        
            po::value<double>(&p_Options->m_KickTheta2)->default_value(p_Options->m_KickTheta2),                                                                                  
            "Angle between the orbital plane and the 'z' axis of the supernovae vector, for the secondary star (default = drawn from kick direction distribution)"
        )

        (
            "luminous-blue-variable-multiplier",                           
            po::value<double>(&p_Options->m_LuminousBlueVariableFactor)->default_value(p_Options->m_LuminousBlueVariableFactor),                                                                  
            ("Multiplicitive constant for LBV mass loss (default = " + std::to_string(p_Options->m_LuminousBlueVariableFactor) + ", use 10 for Mennekens & Vanbeveren 2014)").c_str()
        )

        (
            "mass-ratio-max",                                              
            po::value<double>(&p_Options->m_MassRatioDistributionMax)->default_value(p_Options->m_MassRatioDistributionMax),                                                                      
            ("Maximum mass ratio m2/m1 to generate (default = " + std::to_string(p_Options->m_MassRatioDistributionMax) + ")").c_str()
        )
        (
            "mass-ratio-min",                                              
            po::value<double>(&p_Options->m_MassRatioDistributionMin)->default_value(p_Options->m_MassRatioDistributionMin),                                                                      
            ("Minimum mass ratio m2/m1 to generate (default = " + std::to_string(p_Options->m_MassRatioDistributionMin) + ")").c_str()
        )
        (
            "mass-transfer-fa",                                            
            po::value<double>(&p_Options->m_MassTransferFractionAccreted)->default_value(p_Options->m_MassTransferFractionAccreted),                                                              
            ("Mass Transfer fraction accreted in FIXED prescription (default = " + std::to_string(p_Options->m_MassTransferFractionAccreted) + ", fully conservative)").c_str()
        )
        (
            "mass-transfer-jloss",                                         
            po::value<double>(&p_Options->m_MassTransferJloss)->default_value(p_Options->m_MassTransferJloss),                                                                                    
            ("Specific angular momentum with which the non-accreted system leaves the system (default = " + std::to_string(p_Options->m_MassTransferJloss) + ")").c_str()
        )
        (
            "mass-transfer-thermal-limit-c",                               
            po::value<double>(&p_Options->m_MassTransferCParameter)->default_value(p_Options->m_MassTransferCParameter),                                                                          
            ("Mass Transfer Thermal rate factor fo the accretor (default = " + std::to_string(p_Options->m_MassTransferCParameter) + ")").c_str()
        )
        (
            "maximum-evolution-time",                                      
            po::value<double>(&p_Options->m_MaxEvolutionTime)->default_value(p_Options->m_MaxEvolutionTime),                                                                                      
            ("Maximum time to evolve binaries in Myrs (default = " + std::to_string(p_Options->m_MaxEvolutionTime) + ")").c_str()
        )
        (
            "maximum-mass-donor-nandez-ivanova",                           
            po::value<double>(&p_Options->m_MaximumMassDonorNandezIvanova)->default_value(p_Options->m_MaximumMassDonorNandezIvanova),                                                            
            ("Maximum donor mass allowed for the revised common envelope formalism in Msol (default = " + std::to_string(p_Options->m_MaximumMassDonorNandezIvanova) + ")").c_str()
        )
        (
            "maximum-neutron-star-mass",                                   
            po::value<double>(&p_Options->m_MaximumNeutronStarMass)->default_value(p_Options->m_MaximumNeutronStarMass),                                                                          
            ("Maximum mass of a neutron star (default = " + std::to_string(p_Options->m_MaximumNeutronStarMass) + ")").c_str()
        )
        (
            "mcbur1",                                                      
            po::value<double>(&p_Options->m_mCBUR1)->default_value(p_Options->m_mCBUR1),                                                                                                          
            ("MCBUR1: Min core mass at BAGB to avoid fully degenerate CO core  (default = " + std::to_string(p_Options->m_mCBUR1) + ")").c_str()
        )
        (
            "metallicity,z",                                               
            po::value<double>(&p_Options->m_Metallicity)->default_value(p_Options->m_Metallicity),                                                                                                
            ("Metallicity to use (default " + std::to_string(p_Options->m_Metallicity) + " Zsol)").c_str()
        )
        (
            "minimum-secondary-mass",                                      
            po::value<double>(&p_Options->m_MinimumMassSecondary)->default_value(p_Options->m_MinimumMassSecondary),                                                                              
            ("Minimum mass of secondary to generate in Msol (default = " + std::to_string(p_Options->m_MinimumMassSecondary) + ")").c_str()
        )

        (
            "neutrino-mass-loss-bh-formation-value",                       
            po::value<double>(&p_Options->m_NeutrinoMassLossValueBH)->default_value(p_Options->m_NeutrinoMassLossValueBH),                                                                        
            ("Value corresponding to neutrino mass loss assumption (default = " + std::to_string(p_Options->m_NeutrinoMassLossValueBH) + ")").c_str()
        )

        (
            "orbital-period-max",                                          
            po::value<double>(&p_Options->m_PeriodDistributionMax)->default_value(p_Options->m_PeriodDistributionMax),                                                                            
            ("Maximum period in days to generate (default = " + std::to_string(p_Options->m_PeriodDistributionMax) + ")").c_str()
        )
        (
            "orbital-period-min",                                          
            po::value<double>(&p_Options->m_PeriodDistributionMin)->default_value(p_Options->m_PeriodDistributionMin),                                                                            
            ("Minimum period in days to generate (default = " + std::to_string(p_Options->m_PeriodDistributionMin) + ")").c_str()
        )

        (
            "pisn-lower-limit",                                            
            po::value<double>(&p_Options->m_PairInstabilityLowerLimit)->default_value(p_Options->m_PairInstabilityLowerLimit),                                                                    
            ("Minimum core mass for PISN (default = " + std::to_string(p_Options->m_PairInstabilityLowerLimit) + ")").c_str()
        )
        (
            "pisn-upper-limit",                                            
            po::value<double>(&p_Options->m_PairInstabilityUpperLimit)->default_value(p_Options->m_PairInstabilityUpperLimit),                                                                    
            ("Maximum core mass for PISN (default = " + std::to_string(p_Options->m_PairInstabilityUpperLimit) + ")").c_str()
        )
        (
            "ppi-lower-limit",                                             
            po::value<double>(&p_Options->m_PulsationalPairInstabilityLowerLimit)->default_value(p_Options->m_PulsationalPairInstabilityLowerLimit),                                              
            ("Minimum core mass for PPI (default = " + std::to_string(p_Options->m_PulsationalPairInstabilityLowerLimit) + ")").c_str()
        )
        (
            "ppi-upper-limit",                                             
            po::value<double>(&p_Options->m_PulsationalPairInstabilityUpperLimit)->default_value(p_Options->m_PulsationalPairInstabilityUpperLimit),                                              
            ("Maximum core mass for PPI (default = " + std::to_string(p_Options->m_PulsationalPairInstabilityUpperLimit) + ")").c_str()
        )
        (
            "pulsar-birth-magnetic-field-distribution-max",                
            po::value<double>(&p_Options->m_PulsarBirthMagneticFieldDistributionMax)->default_value(p_Options->m_PulsarBirthMagneticFieldDistributionMax),                                        
            ("Maximum (log10) pulsar birth magnetic field (default = " + std::to_string(p_Options->m_PulsarBirthMagneticFieldDistributionMax) + ")").c_str()
        )
        (
            "pulsar-birth-magnetic-field-distribution-min",                
            po::value<double>(&p_Options->m_PulsarBirthMagneticFieldDistributionMin)->default_value(p_Options->m_PulsarBirthMagneticFieldDistributionMin),                                        
            ("Minimum (log10) pulsar birth magnetic field) (default = " + std::to_string(p_Options->m_PulsarBirthMagneticFieldDistributionMin) + ")").c_str()
        )
        (
            "pulsar-birth-spin-period-distribution-max",                   
            po::value<double>(&p_Options->m_PulsarBirthSpinPeriodDistributionMax)->default_value(p_Options->m_PulsarBirthSpinPeriodDistributionMax),                                              
            ("Maximum pulsar birth spin period in ms (default = " + std::to_string(p_Options->m_PulsarBirthSpinPeriodDistributionMax) + ")").c_str()
        )
        (
            "pulsar-birth-spin-period-distribution-min",                   
            po::value<double>(&p_Options->m_PulsarBirthSpinPeriodDistributionMin)->default_value(p_Options->m_PulsarBirthSpinPeriodDistributionMin),                                              
            ("Minimum pulsar birth spin period in ms (default = " + std::to_string(p_Options->m_PulsarBirthSpinPeriodDistributionMin) + ")").c_str()
        )
        (
            "pulsar-magnetic-field-decay-massscale",                       
            po::value<double>(&p_Options->m_PulsarMagneticFieldDecayMassscale)->default_value(p_Options->m_PulsarMagneticFieldDecayMassscale),                                                    
            ("Mass scale on which magnetic field decays during accretion in solar masses (default = " + std::to_string(p_Options->m_PulsarMagneticFieldDecayMassscale) + ")").c_str()
        )
        (
            "pulsar-magnetic-field-decay-timescale",                       
            po::value<double>(&p_Options->m_PulsarMagneticFieldDecayTimescale)->default_value(p_Options->m_PulsarMagneticFieldDecayTimescale),                                                    
            ("Timescale on which magnetic field decays in Myrs (default = " + std::to_string(p_Options->m_PulsarMagneticFieldDecayTimescale) + ")").c_str()
        )
        (
            "pulsar-minimum-magnetic-field",                               
            po::value<double>(&p_Options->m_PulsarLog10MinimumMagneticField)->default_value(p_Options->m_PulsarLog10MinimumMagneticField),                                                        
            ("log10 of the minimum pulsar magnetic field in Gauss (default = " + std::to_string(p_Options->m_PulsarLog10MinimumMagneticField) + ")").c_str()
        )

        (
            "separation,a",                              
            po::value<double>(&p_Options->m_SemiMajorAxis)->default_value(p_Options->m_SemiMajorAxis),                                                        
            ("Initial semi-major axis, a (default = " + std::to_string(p_Options->m_SemiMajorAxis) + ")").c_str()
        )        
        (
            "semi-major-axis-max",                                         
            po::value<double>(&p_Options->m_SemiMajorAxisDistributionMax)->default_value(p_Options->m_SemiMajorAxisDistributionMax),                                                              
            ("Maximum semi major axis in AU to generate (default = " + std::to_string(p_Options->m_SemiMajorAxisDistributionMax) + ")").c_str()
        )
        (
            "semi-major-axis-min",                                         
            po::value<double>(&p_Options->m_SemiMajorAxisDistributionMin)->default_value(p_Options->m_SemiMajorAxisDistributionMin),                                                              
            ("Minimum semi major axis in AU to generate (default = " + std::to_string(p_Options->m_SemiMajorAxisDistributionMin) + ")").c_str()
        )

        (
            "timestep-multiplier",
            po::value<double>(&p_Options->m_TimestepMultiplier)->default_value(p_Options->m_TimestepMultiplier),
            ("Timestep multiplier for SSE and BSE (default = " + std::to_string(p_Options->m_TimestepMultiplier) + ")").c_str()
        )

        (
            "wolf-rayet-multiplier",                                       
            po::value<double>(&p_Options->m_WolfRayetFactor)->default_value(p_Options->m_WolfRayetFactor),                                                                                        
            ("Multiplicitive constant for WR winds (default = " + std::to_string(p_Options->m_WolfRayetFactor) + ")").c_str()
        )

        (
            "zeta-adiabatic-arbitrary",                                    
            po::value<double>(&p_Options->m_ZetaAdiabaticArbitrary)->default_value(p_Options->m_ZetaAdiabaticArbitrary),                                                                          
            ("Value of mass-radius exponent zeta adiabatic (default = " + std::to_string(p_Options->m_ZetaAdiabaticArbitrary) + ")").c_str()
        )
        (
            "zeta-main-sequence",                                          
            po::value<double>(&p_Options->m_ZetaMainSequence)->default_value(p_Options->m_ZetaMainSequence),                                                                                      
            ("Value of mass-radius exponent zeta on the main sequence (default = " + std::to_string(p_Options->m_ZetaMainSequence) + ")").c_str()
        )
        (
            "zeta-radiative-envelope-giant",                               
            po::value<double>(&p_Options->m_ZetaRadiativeEnvelopeGiant)->default_value(p_Options->m_ZetaRadiativeEnvelopeGiant),                                                                  
            ("Value of mass-radius exponent zeta for radiative envelope giants (default = " + std::to_string(p_Options->m_ZetaRadiativeEnvelopeGiant) + ")").c_str()
        )


        // string options - alphabetically

        // Floor
        /*
        (
            "ais-dcotype",                                                 
            po::value<std::string>(&p_Options->m_AISDCOtypeString)->default_value(p_Options->m_AISDCOtypeString),                                                                                      
            ("DCO type selection in exploratory phase of STROOPWAFEL, (options: [ALL, BBH, BNS, BHNS], default = " + p_Options->m_AISDCOtypeString + ")").c_str()
        )
        */

        (
            "black-hole-kicks",                                            
            po::value<std::string>(&p_Options->m_BlackHoleKicksOptionString)->default_value(p_Options->m_BlackHoleKicksOptionString),                                                                              
            ("Black hole kicks relative to NS kicks (options: [FULL, REDUCED, ZERO, FALLBACK], default = " + p_Options->m_BlackHoleKicksOptionString + ")").c_str()
        )

        (
            "case-bb-stability-prescription",                              
            po::value<std::string>(&p_Options->m_CaseBBStabilityPrescriptionString)->default_value(p_Options->m_CaseBBStabilityPrescriptionString),                                                    
            ("Case BB/BC mass transfer stability prescription (options: [ALWAYS_STABLE, ALWAYS_STABLE_ONTO_NSBH, TREAT_AS_OTHER_MT, ALWAYS_UNSTABLE], default = " + p_Options->m_CaseBBStabilityPrescriptionString + ")").c_str()
        )
        (
            "chemically-homogeneous-evolution",                            
            po::value<std::string>(&p_Options->m_CheString)->default_value(p_Options->m_CheString),                                                                                                    
            ("Chemically Homogeneous Evolution (options: [NONE, OPTIMISTIC, PESSIMISTIC], default = " + p_Options->m_CheString + ")").c_str()
        )
        (
            "common-envelope-lambda-prescription",                         
            po::value<std::string>(&p_Options->m_CommonEnvelopeLambdaPrescriptionString)->default_value(p_Options->m_CommonEnvelopeLambdaPrescriptionString),                                          
            ("CE lambda prescription (options: [LAMBDA_FIXED, LAMBDA_LOVERIDGE, LAMBDA_NANJING, LAMBDA_KRUCKOW, LAMBDA_DEWI], default = " + p_Options->m_CommonEnvelopeLambdaPrescriptionString + ")").c_str()
        )
        (
            "common-envelope-mass-accretion-prescription",                 
            po::value<std::string>(&p_Options->m_CommonEnvelopeMassAccretionPrescriptionString)->default_value(p_Options->m_CommonEnvelopeMassAccretionPrescriptionString),                            
            ("Assumption about whether NS/BHs can accrete mass during common envelope evolution (options: [ZERO, CONSTANT, UNIFORM, MACLEOD], default = " + p_Options->m_CommonEnvelopeMassAccretionPrescriptionString + ")").c_str()
        )
        
        (
            "eccentricity-distribution",                                 
            po::value<std::string>(&p_Options->m_EccentricityDistributionString)->default_value(p_Options->m_EccentricityDistributionString),                                                          
            ("Initial eccentricity distribution (options: [ZERO, FIXED, FLAT, THERMALISED, GELLER+2013], default = " + p_Options->m_EccentricityDistributionString + ")").c_str()
        )
        (
            "envelope-state-prescription",                                 
            po::value<std::string>(&p_Options->m_EnvelopeStatePrescriptionString)->default_value(p_Options->m_EnvelopeStatePrescriptionString),                                                        
            ("Prescription for whether the envelope is radiative or convective (options: [LEGACY, HURLEY, FIXED_TEMPERATURE], default = " + p_Options->m_EnvelopeStatePrescriptionString + ")").c_str()
        )

        (
            "fryer-supernova-engine",                                      
            po::value<std::string>(&p_Options->m_FryerSupernovaEngineString)->default_value(p_Options->m_FryerSupernovaEngineString),                                                                  
            ("If using Fryer et al 2012 fallback prescription, select between 'delayed' and 'rapid' engines (default = " + p_Options->m_FryerSupernovaEngineString + ")").c_str()
        )

        (
            "grid",                                                        
            po::value<std::string>(&p_Options->m_GridFilename)->default_value(p_Options->m_GridFilename)->implicit_value(""),                                                                      
            ("Grid filename (default = " + p_Options->m_GridFilename + ")").c_str()
        )

        (
            "initial-mass-function,i",                                     
            po::value<std::string>(&p_Options->m_InitialMassFunctionString)->default_value(p_Options->m_InitialMassFunctionString),                                                                    
            ("Initial mass function (options: [SALPETER, POWERLAW, UNIFORM, KROUPA], default = " + p_Options->m_InitialMassFunctionString + ")").c_str()
        )

        (
            "kick-direction",                                              
            po::value<std::string>(&p_Options->m_KickDirectionDistributionString)->default_value(p_Options->m_KickDirectionDistributionString),                                                        
            ("Natal kick direction distribution (options: [ISOTROPIC, INPLANE, PERPENDICULAR, POWERLAW, WEDGE, POLES], default = " + p_Options->m_KickDirectionDistributionString + ")").c_str()
        )
        (
            "kick-magnitude-distribution",                                 
            po::value<std::string>(&p_Options->m_KickMagnitudeDistributionString)->default_value(p_Options->m_KickMagnitudeDistributionString),                                                        
            ("Natal kick magnitude distribution (options: [ZERO, FIXED, FLAT, MAXWELLIAN, BRAYELDRIDGE, MULLER2016, MULLER2016MAXWELLIAN, MULLERMANDEL], default = " + p_Options->m_KickMagnitudeDistributionString + ")").c_str()
        )

        // Serena
        /*
        (
            "logfile-be-binaries",                                     
            po::value<std::string>(&p_Options->m_LogfileBeBinaries)->default_value(p_Options->m_LogfileBeBinaries),                                                                              
            ("Filename for Be Binaries logfile (default = " + p_Options->m_LogfileBeBinaries + ")").c_str()
        )
        */

        (
            "logfile-rlof-parameters",                                 
            po::value<std::string>(&p_Options->m_LogfileRLOFParameters)->default_value(p_Options->m_LogfileRLOFParameters),                                                                      
            ("Filename for RLOF Parameters logfile ( default = " + p_Options->m_LogfileRLOFParameters + ")").c_str()
        )
        (
            "logfile-common-envelopes",                                
            po::value<std::string>(&p_Options->m_LogfileCommonEnvelopes)->default_value(p_Options->m_LogfileCommonEnvelopes),                                                                    
            ("Filename for Common Envelopes logfile (default = " + p_Options->m_LogfileCommonEnvelopes + ")").c_str()
        )
        (
            "logfile-detailed-output",                                 
            po::value<std::string>(&p_Options->m_LogfileDetailedOutput)->default_value(p_Options->m_LogfileDetailedOutput),                                                                      
            ("Filename for Detailed Output logfile (default = " + p_Options->m_LogfileDetailedOutput + ")").c_str()
        )
        (
            "logfile-double-compact-objects",                          
            po::value<std::string>(&p_Options->m_LogfileDoubleCompactObjects)->default_value(p_Options->m_LogfileDoubleCompactObjects),                                                          
            ("Filename for Double Compact Objects logfile (default = " + p_Options->m_LogfileDoubleCompactObjects + ")").c_str()
        )
        (
            "logfile-pulsar-evolution",                                
            po::value<std::string>(&p_Options->m_LogfilePulsarEvolution)->default_value(p_Options->m_LogfilePulsarEvolution),                                                                    
            ("Filename for Pulsar Evolution logfile (default = " + p_Options->m_LogfilePulsarEvolution + ")").c_str()
        )
        (
            "logfile-supernovae",                                      
            po::value<std::string>(&p_Options->m_LogfileSupernovae)->default_value(p_Options->m_LogfileSupernovae),                                                                              
            ("Filename for Supernovae logfile (default = " + p_Options->m_LogfileSupernovae + ")").c_str()
        )
        (
            "logfile-system-parameters",                               
            po::value<std::string>(&p_Options->m_LogfileSystemParameters)->default_value(p_Options->m_LogfileSystemParameters),                                                                  
            ("Filename for System Parameters logfile (default = " + p_Options->m_LogfileSystemParameters + ")").c_str()
        )
        (
            "logfile-definitions",                                         
            po::value<std::string>(&p_Options->m_LogfileDefinitionsFilename)->default_value(p_Options->m_LogfileDefinitionsFilename)->implicit_value(""),                                              
            ("Filename for logfile record definitions (default = " + p_Options->m_LogfileDefinitionsFilename + ")").c_str()
        )
        (
            "logfile-delimiter",                                           
            po::value<std::string>(&p_Options->m_LogfileDelimiterString)->default_value(p_Options->m_LogfileDelimiterString),                                                                          
            ("Field delimiter for logfile records (default = " + p_Options->m_LogfileDelimiterString + ")").c_str()
        )
        (
            "logfile-name-prefix",                                         
            po::value<std::string>(&p_Options->m_LogfileNamePrefix)->default_value(p_Options->m_LogfileNamePrefix)->implicit_value(""),                                                                
            ("Prefix for logfile names (default = " + p_Options->m_LogfileNamePrefix + ")").c_str()
        )
        (
            "logfile-switch-log",                                      
            po::value<std::string>(&p_Options->m_LogfileSwitchLog)->default_value(p_Options->m_LogfileSwitchLog),                                                                                
            ("Filename for Switch Log logfile (default = " + p_Options->m_LogfileSwitchLog + ")").c_str()
        )

        (
            "mass-loss-prescription",                                      
            po::value<std::string>(&p_Options->m_MassLossPrescriptionString)->default_value(p_Options->m_MassLossPrescriptionString),                                                                  
            ("Mass loss prescription (options: [NONE, HURLEY, VINK], default = " + p_Options->m_MassLossPrescriptionString + ")").c_str()
        )
        (
            "mass-ratio-distribution,q",                                   
            po::value<std::string>(&p_Options->m_MassRatioDistributionString)->default_value(p_Options->m_MassRatioDistributionString),                                                                
            ("Initial mass ratio distribution for q=m2/m1 (options: [FLAT, DUQUENNOYMAYOR1991, SANA2012], default = " + p_Options->m_MassRatioDistributionString + ")").c_str()
        )
        (
            "mass-transfer-accretion-efficiency-prescription",             
            po::value<std::string>(&p_Options->m_MassTransferAccretionEfficiencyPrescriptionString)->default_value(p_Options->m_MassTransferAccretionEfficiencyPrescriptionString),                    
            ("Mass Transfer Accretion Efficiency prescription (options: [THERMAL, FIXED], default = " + p_Options->m_MassTransferAccretionEfficiencyPrescriptionString + ")").c_str()
        )
        (
            "mass-transfer-angular-momentum-loss-prescription",            
            po::value<std::string>(&p_Options->m_MassTransferAngularMomentumLossPrescriptionString)->default_value(p_Options->m_MassTransferAngularMomentumLossPrescriptionString),                    
            ("Mass Transfer Angular Momentum Loss prescription (options: [JEANS, ISOTROPIC, CIRCUMBINARY, ARBITRARY], default = " + p_Options->m_MassTransferAngularMomentumLossPrescriptionString + ")").c_str()
        )
        (
            "mass-transfer-rejuvenation-prescription",                     
            po::value<std::string>(&p_Options->m_MassTransferRejuvenationPrescriptionString)->default_value(p_Options->m_MassTransferRejuvenationPrescriptionString),                                  
            ("Mass Transfer Rejuvenation prescription (options: [NONE, STARTRACK], default = " + p_Options->m_MassTransferRejuvenationPrescriptionString + ")").c_str()
        )
        (
            "mass-transfer-thermal-limit-accretor",                        
            po::value<std::string>(&p_Options->m_MassTransferThermallyLimitedVariationString)->default_value(p_Options->m_MassTransferThermallyLimitedVariationString),                                
            ("Mass Transfer Thermal Accretion limit (default = " + p_Options->m_MassTransferThermallyLimitedVariationString + ")").c_str()
        )

        (
            "mode",                                                 
            po::value<std::string>(&p_Options->m_EvolutionModeString)->default_value(p_Options->m_EvolutionModeString),                                                                              
            ("Evolution mode (options: [SSE, BSE], default = " + p_Options->m_EvolutionModeString + ")").c_str()
        )

        (
            "neutrino-mass-loss-bh-formation",                             
            po::value<std::string>(&p_Options->m_NeutrinoMassLossAssumptionBHString)->default_value(p_Options->m_NeutrinoMassLossAssumptionBHString),                                                  
            ("Assumption about neutrino mass loss during BH formation (options: [FIXED_FRACTION, FIXED_MASS], default = " + p_Options->m_NeutrinoMassLossAssumptionBHString + ")").c_str()
        )
        (
            "neutron-star-equation-of-state",                              
            po::value<std::string>(&p_Options->m_NeutronStarEquationOfStateString)->default_value(p_Options->m_NeutronStarEquationOfStateString),                                                      
            ("Neutron star equation of state to use (options: [SSE, ARP3], default = " + p_Options->m_NeutronStarEquationOfStateString + ")").c_str()
        )

        (
            "output-container,c",                                          
            po::value<std::string>(&p_Options->m_OutputContainerName)->default_value(p_Options->m_OutputContainerName)->implicit_value(""),                                                            
            ("Container (directory) name for output files (default = " + p_Options->m_OutputContainerName + ")").c_str()
        )
        (
            "outputPath,o",                                                
            po::value<std::string>(&p_Options->m_OutputPathString)->default_value(p_Options->m_OutputPathString)->implicit_value(""),                                                                  
            ("Directory for output (default = " + p_Options->m_OutputPathString + ")").c_str()
        )

        (
            "pulsar-birth-magnetic-field-distribution",                    
            po::value<std::string>(&p_Options->m_PulsarBirthMagneticFieldDistributionString)->default_value(p_Options->m_PulsarBirthMagneticFieldDistributionString),                                  
            ("Pulsar Birth Magnetic Field distribution (options: [ZERO, FIXED, FLATINLOG, UNIFORM, LOGNORMAL], default = " + p_Options->m_PulsarBirthMagneticFieldDistributionString + ")").c_str()
        )
        (
            "pulsar-birth-spin-period-distribution",                       
            po::value<std::string>(&p_Options->m_PulsarBirthSpinPeriodDistributionString)->default_value(p_Options->m_PulsarBirthSpinPeriodDistributionString),                                        
            ("Pulsar Birth Spin Period distribution (options: [ZERO, FIXED, UNIFORM, NORMAL], default = " + p_Options->m_PulsarBirthSpinPeriodDistributionString + ")").c_str()
        )
        (
            "pulsational-pair-instability-prescription",                   
            po::value<std::string>(&p_Options->m_PulsationalPairInstabilityPrescriptionString)->default_value(p_Options->m_PulsationalPairInstabilityPrescriptionString),                              
            ("Pulsational Pair Instability prescription (options: [COMPAS, STARTRACK, MARCHANT], default = " + p_Options->m_PulsationalPairInstabilityPrescriptionString + ")").c_str()
        )

        (
            "remnant-mass-prescription",                                   
            po::value<std::string>(&p_Options->m_RemnantMassPrescriptionString)->default_value(p_Options->m_RemnantMassPrescriptionString),                                                            
            ("Choose remnant mass prescription (options: [HURLEY2000, BELCZYNSKI2002, FRYER2012, MULLER2016, MULLERMANDEL], default = " + p_Options->m_RemnantMassPrescriptionString + ")").c_str()
        )
        (
            "rotational-velocity-distribution",                            
            po::value<std::string>(&p_Options->m_RotationalVelocityDistributionString)->default_value(p_Options->m_RotationalVelocityDistributionString),                                              
            ("Initial rotational velocity distribution (options: [ZERO, HURLEY, VLTFLAMES], default = " + p_Options->m_RotationalVelocityDistributionString + ")").c_str()
        )

        (
            "semi-major-axis-distribution",                              
            po::value<std::string>(&p_Options->m_SemiMajorAxisDistributionString)->default_value(p_Options->m_SemiMajorAxisDistributionString),                                                        
            ("Initial semi-major axis distribution (options: [FLATINLOG, CUSTOM, DUQUENNOYMAYOR1991, SANA2012], default = " + p_Options->m_SemiMajorAxisDistributionString + ")").c_str()
        )        
        (
            "stellar-zeta-prescription",                                   
            po::value<std::string>(&p_Options->m_StellarZetaPrescriptionString)->default_value(p_Options->m_StellarZetaPrescriptionString),                                                            
            ("Prescription for stellar zeta (default = " + p_Options->m_StellarZetaPrescriptionString + ")").c_str()
        )


        // vector (list) options - alphabetically

        (
            "debug-classes",                                               
            po::value<vector<std::string>>(&p_Options->m_DebugClasses)->multitoken()->default_value(p_Options->m_DebugClasses),                                                                        
            ("Debug classes enabled (default = " + defaultDebugClasses + ")").c_str()
        )
        (
            "log-classes",                                                 
            po::value<vector<std::string>>(&p_Options->m_LogClasses)->multitoken()->default_value(p_Options->m_LogClasses),                                                                            
            ("Logging classes enabled (default = " + defaultLogClasses + ")").c_str()
        )
    
        ;   // end the list of options to be added

    }
    catch (po::error& e) {      // program options exception
        ok = false;             // set status
    }
    catch (...) {               // unhandled exception
        ok = false;             // set status
    }

    return ok;
}


/*
 * Retrieve the attributes of an option
 * The option for which the attributes are to be retreived is passed as an iterator
 * pointing at the option in the boost variables map
 * 
 * The attributes are returned as a tuple, described by typedef ATTR, containing
 * 
 *     - dataType       TYPENAME (high-level) data type of the attribute.  Will be one of {NONE, BOOL, INT, FLOAT, STRING}
 *     - defaulted      BOOL     flag to indicate if the option was specified by the user or defaulted to the defaulty value
 *     - typeStr        STRING   detailed data type returned as a string (e.g. "UNSIGNED LONG INT" etc.)
 *     - valueStr       STRING   the value of the option returned as a string (e.g. "2.3", "BSE" etc.)
 * 
 * 
 * Options::ATTR OptionAttributes(const po::variables_map p_VM, const po::variables_map::const_iterator p_IT)
 * 
 *
 * @param   [IN]    p_VM                        The boost variables map
 * @param   [IN]    p_IT                        Iterator for the boost variables map pointing to the option
 *                                              for which the attributed are required to be retrieved
 * @return                                      Tuple (type ATTR) containing the option attributes
 */
Options::ATTR Options::OptionAttributes(const po::variables_map p_VM, const po::variables_map::const_iterator p_IT) {
            
    TYPENAME    dataType  = TYPENAME::NONE;
    std::string typeStr   = "";
    bool        defaulted = false;
    std::string valueStr  = "";

    if (((boost::any)p_IT->second.value()).empty()) return std::make_tuple(TYPENAME::NONE, true, "", "");   // empty option 

    // determine if option values was supplied, or whether the default was used

    defaulted = (p_VM[p_IT->first].defaulted() || p_IT->second.defaulted());

    // find data type and format the option value into a string
    // handles most data types - add others if they cause problems

    bool isCharPtr = false;
    bool isStr     = false;

    // (pre)check for data type = charPtr
    try {
        boost::any_cast<const char *>(p_IT->second.value());
        isCharPtr = true;
    }
    catch (const boost::bad_any_cast &) {
        isCharPtr = false;
    }

    if (!isCharPtr) {
        // (pre)check for data type = string
        try {
            boost::any_cast<std::string>(p_IT->second.value());
            isStr = true;
        }
        catch (const boost::bad_any_cast &) {
            isStr = false;
        }
    }

    // find other data types
    // it's not pretty, but it works

    if (isCharPtr) { 
        dataType = TYPENAME::NONE;                                          // not supported by COMPAS as an option data type
        typeStr  = "CONST_CHAR_*";                                          // ... but we know what type it is, and
        valueStr = p_VM[p_IT->first].as<const char *>();                    // ... we can still format the value
    }

    else if (isStr) {
        dataType = TYPENAME::STRING;
        typeStr  = "STRING";
        std::string tmp = p_VM[p_IT->first].as<std::string>();
        if (tmp.size()) valueStr = "'" + tmp + "'";
        else            valueStr = "''";
    }

    else if (((boost::any)p_IT->second.value()).type() == typeid(signed                )) { dataType = TYPENAME::INT;   typeStr = "SIGNED";                 valueStr = std::to_string(p_VM[p_IT->first].as<signed                >()); }
    else if (((boost::any)p_IT->second.value()).type() == typeid(unsigned              )) { dataType = TYPENAME::INT;   typeStr = "UNSIGNED";               valueStr = std::to_string(p_VM[p_IT->first].as<unsigned              >()); }

    else if (((boost::any)p_IT->second.value()).type() == typeid(short                 )) { dataType = TYPENAME::INT;   typeStr = "SHORT";                  valueStr = std::to_string(p_VM[p_IT->first].as<short                 >()); }
    else if (((boost::any)p_IT->second.value()).type() == typeid(signed short          )) { dataType = TYPENAME::INT;   typeStr = "SIGNED_SHORT";           valueStr = std::to_string(p_VM[p_IT->first].as<signed short          >()); }
    else if (((boost::any)p_IT->second.value()).type() == typeid(unsigned short        )) { dataType = TYPENAME::INT;   typeStr = "UNSIGNED_SHORT";         valueStr = std::to_string(p_VM[p_IT->first].as<unsigned short        >()); }

    else if (((boost::any)p_IT->second.value()).type() == typeid(short int             )) { dataType = TYPENAME::INT;   typeStr = "SHORT_INT";              valueStr = std::to_string(p_VM[p_IT->first].as<short int             >()); }
    else if (((boost::any)p_IT->second.value()).type() == typeid(signed short int      )) { dataType = TYPENAME::INT;   typeStr = "SIGNED_SHORT_INT";       valueStr = std::to_string(p_VM[p_IT->first].as<signed short int      >()); }
    else if (((boost::any)p_IT->second.value()).type() == typeid(unsigned short int    )) { dataType = TYPENAME::INT;   typeStr = "UNSIGNED_SHORT_INT";     valueStr = std::to_string(p_VM[p_IT->first].as<unsigned short int    >()); }

    else if (((boost::any)p_IT->second.value()).type() == typeid(int                   )) { dataType = TYPENAME::INT;   typeStr = "INT";                    valueStr = std::to_string(p_VM[p_IT->first].as<int                   >()); }
    else if (((boost::any)p_IT->second.value()).type() == typeid(signed int            )) { dataType = TYPENAME::INT;   typeStr = "SIGNED_INT";             valueStr = std::to_string(p_VM[p_IT->first].as<signed int            >()); }
    else if (((boost::any)p_IT->second.value()).type() == typeid(unsigned int          )) { dataType = TYPENAME::INT;   typeStr = "UNSIGNED_INT";           valueStr = std::to_string(p_VM[p_IT->first].as<unsigned int          >()); }

    else if (((boost::any)p_IT->second.value()).type() == typeid(long                  )) { dataType = TYPENAME::INT;   typeStr = "LONG";                   valueStr = std::to_string(p_VM[p_IT->first].as<long                  >()); }
    else if (((boost::any)p_IT->second.value()).type() == typeid(signed long           )) { dataType = TYPENAME::INT;   typeStr = "SIGNED_LONG";            valueStr = std::to_string(p_VM[p_IT->first].as<signed long           >()); }
    else if (((boost::any)p_IT->second.value()).type() == typeid(unsigned long         )) { dataType = TYPENAME::INT;   typeStr = "UNSIGNED_LONG";          valueStr = std::to_string(p_VM[p_IT->first].as<unsigned long         >()); }

    else if (((boost::any)p_IT->second.value()).type() == typeid(long int              )) { dataType = TYPENAME::INT;   typeStr = "LONG_INT";               valueStr = std::to_string(p_VM[p_IT->first].as<long int              >()); }
    else if (((boost::any)p_IT->second.value()).type() == typeid(signed long int       )) { dataType = TYPENAME::INT;   typeStr = "SIGNED_LONG_INT";        valueStr = std::to_string(p_VM[p_IT->first].as<signed long int       >()); }
    else if (((boost::any)p_IT->second.value()).type() == typeid(unsigned long int     )) { dataType = TYPENAME::INT;   typeStr = "UNSIGNED_LONG_INT";      valueStr = std::to_string(p_VM[p_IT->first].as<unsigned long int     >()); }

    else if (((boost::any)p_IT->second.value()).type() == typeid(long long             )) { dataType = TYPENAME::INT;   typeStr = "LONG_LONG";              valueStr = std::to_string(p_VM[p_IT->first].as<long long             >()); }
    else if (((boost::any)p_IT->second.value()).type() == typeid(signed long long      )) { dataType = TYPENAME::INT;   typeStr = "SIGNED_LONG_LONG";       valueStr = std::to_string(p_VM[p_IT->first].as<signed long long      >()); }
    else if (((boost::any)p_IT->second.value()).type() == typeid(unsigned long long    )) { dataType = TYPENAME::INT;   typeStr = "UNSIGNED_LONG_LONG";     valueStr = std::to_string(p_VM[p_IT->first].as<unsigned long long    >()); }

    else if (((boost::any)p_IT->second.value()).type() == typeid(long long int         )) { dataType = TYPENAME::INT;   typeStr = "LONG_LONG_INT";          valueStr = std::to_string(p_VM[p_IT->first].as<long long int         >()); }
    else if (((boost::any)p_IT->second.value()).type() == typeid(signed long long int  )) { dataType = TYPENAME::INT;   typeStr = "SIGNED_LONG_LONG_INT";   valueStr = std::to_string(p_VM[p_IT->first].as<signed long long int  >()); }
    else if (((boost::any)p_IT->second.value()).type() == typeid(unsigned long long int)) { dataType = TYPENAME::INT;   typeStr = "UNSIGNED_LONG_LONG_INT"; valueStr = std::to_string(p_VM[p_IT->first].as<unsigned long long int>()); }

    else if (((boost::any)p_IT->second.value()).type() == typeid(float                 )) { dataType = TYPENAME::FLOAT; typeStr = "FLOAT";                  valueStr = std::to_string(p_VM[p_IT->first].as<float                 >()); }
    else if (((boost::any)p_IT->second.value()).type() == typeid(double                )) { dataType = TYPENAME::FLOAT; typeStr = "DOUBLE";                 valueStr = std::to_string(p_VM[p_IT->first].as<double                >()); }
    else if (((boost::any)p_IT->second.value()).type() == typeid(long double           )) { dataType = TYPENAME::FLOAT; typeStr = "LONG_DOUBLE";            valueStr = std::to_string(p_VM[p_IT->first].as<long double           >()); }

    else if (((boost::any)p_IT->second.value()).type() == typeid(char                  )) { dataType = TYPENAME::INT;   typeStr = "CHAR";                   valueStr = std::to_string(p_VM[p_IT->first].as<char                  >()); }
    else if (((boost::any)p_IT->second.value()).type() == typeid(signed char           )) { dataType = TYPENAME::INT;   typeStr = "SIGNED_CHAR";            valueStr = std::to_string(p_VM[p_IT->first].as<signed char           >()); }
    else if (((boost::any)p_IT->second.value()).type() == typeid(unsigned char         )) { dataType = TYPENAME::INT;   typeStr = "UNSIGNED_CHAR";          valueStr = std::to_string(p_VM[p_IT->first].as<unsigned char         >()); }

    else if (((boost::any)p_IT->second.value()).type() == typeid(bool)) {
        dataType = TYPENAME::BOOL;
        typeStr  = "BOOL";
        valueStr = p_VM[p_IT->first].as<bool>() ? "TRUE" : "FALSE";
    } 

    else {  // Assume vector<string>
        try {
            std::ostringstream elemsSS;
            elemsSS << "{ ";
            vector<std::string> tmp = p_VM[p_IT->first].as<vector<std::string>>();
            for (std::vector<std::string>::iterator elem=tmp.begin(); elem != tmp.end(); elem++) {
                elemsSS << "'" << (*elem) << "', ";
            }
            std::string elems = elemsSS.str();
            if (elems.size() > 2) elems.erase(elems.size() - 2);
            else if (elems.size() == 2) elems.erase(elems.size() - 1);
            elems += " }";

            dataType = TYPENAME::NONE;                                                  // not supported by COMPAS as an option data type            
            typeStr  = "VECTOR<STRING>";                                                // ... but we know what type it is, and
            valueStr = elems;                                                           // ... we can still format the value
        }
        catch (const boost::bad_any_cast &) {
            dataType = TYPENAME::NONE;                                                  // unknown data type               
            typeStr  = "<UNKNOWN_DATA_TYPE>";
            valueStr = "<UNKNOWN_DATA_TYPE>";
        }
    }

    return std::make_tuple(dataType, defaulted, typeStr, valueStr);
}


/*
 * Build the output string for the Run_Details file
 * The parameter passed is the opstions descriptor - the grid line options, or
 * the commandline options.  Ordinarily we would build the Run_Details contents
 * from the commandline options, but the flexibility exists to use a set of
 * grid line options (maybe one day we will want to (optionally) produce a
 * per star/binary Run_Details file)
 * 
 * 
 * std::string OptionDetails(const OptionsDescriptorT &p_Options)
 * 
 * @param   [IN]    p_Options                   The options descriptor to use to build the output string
 * @return                                      String containing the Run_Details file contents
 */
std::string Options::OptionDetails(const OptionsDescriptorT &p_Options) {
            
    TYPENAME    dataType  = TYPENAME::NONE;
    std::string typeStr   = "";
    bool        defaulted = false;
    std::string valueStr  = "";

    std::ostringstream ss;                                                                                                              // output string

    ss << "COMMAND LINE OPTIONS\n-------------------\n\n";

    for (po::variables_map::const_iterator it = p_Options.optionValues.m_VM.begin(); it != p_Options.optionValues.m_VM.end(); it++) {   // for all options in the variable map
  
        ss << it->first << " = ";                                                                                                       // add option name to output string

        std::tie(dataType, defaulted, typeStr, valueStr) = OptionAttributes(p_Options.optionValues.m_VM, it);                           // get option attributes

        if (valueStr == "")                                                                                                             // empty option?
            ss << "<EMPTY_OPTION>\n";                                                                                                   // yes - say so
        else                                                                                                                            // no
            ss << valueStr + ", " << (defaulted ? "DEFAULT_USED, " : "USER_SUPPLIED, ") << typeStr << "\n";                             // add option details to output string
    }
  
    ss << "\n\nOTHER PARAMETERS\n----------------\n\n";

    ss << "useFixedUK         = " << (p_Options.optionValues.m_UseFixedUK ? "TRUE" : "FALSE") << ", CALCULATED, BOOL\n";                // useFixedUK
    ss << "outputPath         = " <<  p_Options.optionValues.m_OutputPath.string() << ", CALCULATED, STRING\n";                         // outputPath (fully qualified)
    ss << "fixedRandomSeed    = " << (p_Options.optionValues.m_FixedRandomSeed ? "TRUE" : "FALSE") << ", CALCULATED, BOOL\n";           // fixedRandomSeed

    return ss.str();
}


/*
 * Parse the options provided by the user
 * 
 * Before we give the options to boost we need to determine if the user passed
 * any ranges or sets and, if they did, handle those - boost doesn't know anything
 * about them.
 *
 * A range is allowed only for numeric options (i.e. INT or FLOAT types), but is not 
 * allowed for all numeric options (e.g. --log-level)
 * A set is allowed for numeric, string, and bool options - but not all of them (e.g. --quiet)
 * 
 * We define a vector of options excluded from the range and set constructs (one vector each).
 * We don't need to exclude non-numeric options from range here - that is done later - here we 
 * just exclude options for which range/set makes no sense.  We have defined vectors of
 * option names that are excluded from ranges (m_RangeExluded) and sets (m_SetExcluded).
 * 
 * 
 * std::string ParseOptionValues(int p_ArgCount, char *p_ArgStrings[], OptionsDescriptorT &p_OptionsDescriptor)
 * 
 * 
 * @param   [IN]    p_ArgCount                  The number of argument strings. (note below for p_ArgStrings)
 * @param   [IN]    p_ArgStrings                The argument strings.   The first argument string is expected
 *                                              (by boost) to be the executable name (boost expects the arguments
 *                                              to be commandline arguments passed to main())
 * @param   [IN]    p_OptionsDescriptor         Struct containing options descriptions.  This struct holds the
 *                                              boost options_description object, the option valued, and a struct
 *                                              containing the complex option values (the ranges and sets)
 * @return                                      String containing an error string
 *                                              If no error occurred the return string will be the empty string 
 */
std::string Options::ParseOptionValues(int p_ArgCount, char *p_ArgStrings[], OptionsDescriptorT &p_OptionsDescriptor) {

    bool error         = false;                                                                                             // for now...
    std::string errStr = "";                                                                                                // also for now...

    try {

        p_OptionsDescriptor.complexOptionValues = {};                                                                       // initially

        std::string  optionName        = "";                                                                                // option name
        COMPLEX_TYPE type              = COMPLEX_TYPE::NONE;                                                                // complex arg type (range, set, neither/none)
        std::vector<std::string> parms = {};                                                                                // the range or set parameters

        for (int iArg = 0; iArg < p_ArgCount; iArg++) {                                                                     // for each arg string

            if (iArg == 0) continue;                                                                                        // ignore the executable name

            type = COMPLEX_TYPE::NONE;                                                                                      // initially

            optionName = p_ArgStrings[iArg - 1];                                                                            // get the option name for the argument we're processing
            if (optionName[0] == '-') optionName.erase(0, optionName.find_first_not_of("-"));                               // remove the "-" or "--"

            if (p_ArgStrings[iArg] != nullptr) {                                                                            // null arg?
                                                                                                                            // no
                std::string str(p_ArgStrings[iArg]);                                                                        // convert char* to std::string
                str = utils::ToLower(utils::trim(str));                                                                     // downshift - for comparisons

                // check for RANGE or SET
                // range is indicated by 'range[start,count,inc]', 'r[start,count,inc]', or just '[start,count,inc]'

                if ((str[0] == '[') || (str.rfind("r[", 0) == 0) || (str.rfind("range[", 0) == 0)) {                        // starts with '[', 'r[' or 'range[', so...
                    error = true;                                                                                           // unless set otherwise
                    if (str[str.length()-1] == ']') {                                                                       // ... needs to end with ']' to be a valid RANGE
                        type = COMPLEX_TYPE::RANGE;                                                                         // it did - so RANGE

                        // check for RANGE requested for option in range excluded list
                        if (iArg > 1) {                                                                                     // range not valid for arg[1]
                            if (std::find(m_RangeExcluded.begin(), m_RangeExcluded.end(), optionName) != m_RangeExcluded.end())
                                errStr = ERR_MSG(ERROR::ARGUMENT_RANGE_NOT_SUPPORTED) + std::string("'") + optionName + std::string("'");
                            else
                                error = false;                                                                              // we're good
                        }
                    }
                }

                if (!error) {                                                                                               // still ok?
                                                                                                                            // yes
                    // set is indicated by 'set[elem1,elem2,...,elemN]', or 's[elem1,elem2,...,elemN]'

                    if ((str.rfind("s[", 0) == 0) || (str.rfind("set[", 0) == 0)) {                                         // starts with 's[' or 'set[', so ...
                        error = true;                                                                                       // unless set otherwise
                        if (str[str.length()-1] == ']') {                                                                   // ... needs to end with ']' to be a valid SET
                            type = COMPLEX_TYPE::SET;                                                                       // it did - so SET

                            // check for SET requested for option in set excluded list
                            if (iArg > 1) {                                                                                 // set not valid for arg[1]
                                if (std::find(m_SetExcluded.begin(), m_SetExcluded.end(), optionName) != m_SetExcluded.end())
                                    errStr = ERR_MSG(ERROR::ARGUMENT_SET_NOT_SUPPORTED) + std::string("'") + optionName + std::string("'");
                                else
                                    error = false;                                                                          // we're good
                            }
                        }
                    }
                }

                if (!error && type != COMPLEX_TYPE::NONE) {                                                                 // range or set?
                                                                                                                            // yes
                    // we have what looks like a 'range' or 'set' argument
                    // for now, just stash the details away and substitute the
                    // first value for the argument so we can check parsing

                    // look for comma separated values - there will be no 
                    // spaces - the OS/shell would have complained...

                    if (str.rfind("range", 0) == 0) str.erase(0, 5);                                                        // strip 'range' (range indicator) if present
                    if (str.rfind("set", 0) == 0) str.erase(0, 3);                                                          // strip 'set' (set indicator) if present
                    if (str[0] == 'r' || str[0] == 's') str.erase(0, 1);                                                    // strip 'r' or 's' (range or set indicator) if present
                    str = str.substr(1, str.size() - 2);                                                                    // strip enclosing brackets (must be present)

                    if (str.length() == 0 || str[str.length() - 1] == ',') error = true;                                    // no values, or trailing comma is an error
                    else {

                        parms.clear();                                                                                      // start empty

                        size_t start = 0;                                                                                   // start position
                        size_t pos   = 0;                                                                                   // current position
                        while (!error && start < str.length() && pos != std::string::npos) {                                // comma found before the end of the string?

                            std::string value = "";                                                                         // value

                            pos = str.find(",", start);                                                                     // next comma
                                                                                                        
                            if ((pos - start) > 0) {                                                                        // non-zero length string?
                                value = str.substr(start, pos - start);                                                     // yes - grab it
                                parms.push_back(value);                                                                     // store value

                                start = pos + 1;                                                                            // next start
                            }
                            else error = true;                                                                              // empty value - stop
                        }

                        if (!error) {                                                                                       // still ok?
                                                                                                                            // yes
                            if (type == COMPLEX_TYPE::RANGE && parms.size() != 3) {                                         // ranges require exactly 3 parameters
                                error  = true;                                                                              // error
                                errStr = ERR_MSG(ERROR::ARGUMENT_RANGE_NUM_PARMS); 
                            }
                            else {
                                // if range, then we have 3 parameters (checked above)
                                // if set, we have at least one parameter (check earlier), so we're good

                                RangeOrSetDescriptorT details = {type, TYPENAME::NONE, parms, {}, 0};                       // dummy values for datatype and numerical parms
                                p_OptionsDescriptor.complexOptionValues.push_back(std::make_tuple(optionName, details));    // store the range/set

                                strncpy(p_ArgStrings[iArg], parms[0].c_str(), parms[0].length());                           // replace arg value (temporarily)
                                p_ArgStrings[iArg][parms[0].length()] = '\0';
                            }
                        }
                    }
                }          
            }
            if (error) break;                                                                                               // stop parsing if error encountered
        }

        // boost parse_command_line() expects the first arg to be the program name
        // (it thinks it is getting the values that were passed to main() from the 
        // OS/shell), so for options from a grid file we insert a dummy argument as 
        // arg[0] and set the argument count appropriately.
        //
        // if valid ranges or sets were specified by the user they've been temporariliy
        // replaced for the boost parse, but if they were not valid they've been left
        // in the argument strings that will be passed to boost - so boost will fail
        // and complain about the offending parameter (which is what we want)

        po::parsed_options const parsedOptions = po::parse_command_line(p_ArgCount, p_ArgStrings, p_OptionsDescriptor.optionDescriptions);  // parse user-supplied options
        po::store(parsedOptions, p_OptionsDescriptor.optionValues.m_VM);                                                    // store parsed options into variable map
        po::notify(p_OptionsDescriptor.optionValues.m_VM);                                                                  // populate the variables with option values

        // If we've made it this far then boost parsed the commandline arguments ok.
        //
        // If there were any ranges or sets specified by the user we can now work
        // out the data types of the options for which they (the ranges/sets) were
        // specified and sanity check them.
        //
        // But... because of the way Boost works (maybe one day we'll ditch Boost
        // and parse the options ourselves...) we also have to replace the values
        // for any options changed here manually.  Let me explain...  The Boost
        // options_description object holds the descriptions of all the options we
        // want to allow users to specify - either on the commandline or in the
        // grid file.  We populate that object by calling AddOptions() - and we do
        // that once in Initialise().  We don't really want to do it for every
        // grid file record - that would be unnecessary overhead.  Except that it
        // is kind-of necessary because, as I said earlier, of the way Boost works.
        // So, once we've populated the options_description object, we can then
        // have Boost parse (what it thinks are) our commandline options - and when
        // it does that it updates the options_description object to indicate which
        // options were specified on the commandline and have been parsed and values
        // set.  When it does that it marks the options with values specified as
        // "final" and, no matter how many times we parse a new grid line with
        // options specified, if the option values were previously marked as "final"
        // by Boost, it steadfastly refuses to change the value.  Fair enough -
        // that's how it works - but it makes doing what we do here a little more
        // complicated.  It turns out I can't get in and change the "final" flag
        // inside the Boost construct (not easily anyway, that I know of).  What I
        // can do, though, is overwrite the value that is considered to be "final"!
        // Go figure... But I have to get in and do it manually - Boost won't do it
        // for me (it considers the value "final").  That's ok, I can do that - and
        // it's still quicker than repopulating the options_description object for
        // each grid line (although, having now written the code I'm starting to doubt
        // the veracity of that statement - given how much work I have to do to
        // replace the value it might actually be faster just to repopulate - maybe
        // we should benchmark it one day).
        //
        // So... here we not only sanity check the ranges and sets, we also manually
        // set the values for any options that were specified on the commandline or
        // grid line we just parsed.  We set the values for options that have ranges
        // or sets specified to the first value in the range or set.
        //
        //
        // First iterate through the specified ranges and sets to sanity check - and
        // manually set the value for the option to the first value in the range or set.
        // Need to check:
        //     - ranges have not been specified for non-numeric options
        //     - range values are all numeric
        //     - values for ranges and sets match the data type of the option
        //       (i.e. INTs for integer options, FLOATs for fp options, STRINGs for string options)
        //
        // It would be preferable to range check values against valid values for specific
        // options here but, for now at least, too problematic - we'll just let the
        // evolution fail if an option value is bad (we don't want to read through every
        // record of a grid file and check...)
        //
        // Next, iterate through all other options that were specified on the commandline
        // or grid line and manually set the value as required.


        if (errStr.empty()) {                                                                                               // no need if we've already flagged an error

            RangeOrSetDescriptorT details = {};

            size_t count = p_OptionsDescriptor.complexOptionValues.size();                                                  // count of complex values (ranges or sets)
            for (size_t idx = 0; idx < count; idx++) {                                                                      // for each range or set specified

                error = false;                                                                                              // for now...

                optionName = get<0>(p_OptionsDescriptor.complexOptionValues[idx]);                                          // the option name
                details    = get<1>(p_OptionsDescriptor.complexOptionValues[idx]);                                          // range/set details for this optionName
                type       = details.type;                                                                                  // range or set
                parms      = details.parameters;                                                                            // range/set parameter values (as strings)

                po::variables_map::const_iterator it = p_OptionsDescriptor.optionValues.m_VM.find(optionName);              // yes - find the option in the boost variables map
                if (it != p_OptionsDescriptor.optionValues.m_VM.end()) {                                                    // found?
                    TYPENAME dataType = TYPENAME::NONE;                                                                     // yes
                    std::tie(dataType, std::ignore, std::ignore, std::ignore) = OptionAttributes(p_OptionsDescriptor.optionValues.m_VM, it); // data type
                    details.dataType = dataType;                                                                            // set data type

                    if (idx == (count - 1)) details.currPos = 0;                                                            // initial position for inner iterator

                    if (type == COMPLEX_TYPE::RANGE) {                                                                      // RANGE?
                        if (dataType != TYPENAME::INT && dataType != TYPENAME::FLOAT) {                                     // yes - numeric?
                            error  = true;                                                                                  // no - that's not ok
                            errStr = ERR_MSG(ERROR::ARGUMENT_RANGE_NOT_SUPPORTED) + std::string("'") + optionName + std::string("'");
                        }
                        else {                                                                                              // yes - numeric
                                                                                                                            // yes - determine numerical range parameters
                            if (dataType == TYPENAME::INT) {                                                                // option data type is INT?
                                                                                                                            // yes - expecting integer values for start, count, and inc
                                try {
                                    RangeParameterT tmp = {0.0};                                                            // dummy value
                                    details.rangeParms = {tmp, tmp, tmp};                                                   // create the vector

                                    details.rangeParms[0].iVal = std::stoi(details.parameters[0]);                          // integer start
                                    details.rangeParms[1].iVal = std::stoi(details.parameters[1]);                          // integer count
                                    details.rangeParms[2].iVal = std::stoi(details.parameters[2]);                          // integer inc

                                    p_OptionsDescriptor.complexOptionValues[idx] = std::make_tuple(optionName, details); // reset values
                                }
                                catch (const std::out_of_range& e) {
                                    errStr = ERR_MSG(ERROR::ARGUMENT_RANGE_PARMS_EXPECTED_INT) + std::string("'") + optionName  + std::string("'");
                                }
                            }
                            else {                                                                                          // no - expecting fp values for start and inc, and integer for count
                                try {
                                    RangeParameterT tmp = {0.0};                                                            // dummy value
                                    details.rangeParms = {tmp, tmp, tmp};                                                   // create the vector

                                    details.rangeParms[0].fVal = std::stod(details.parameters[0]);                          // floating point start
                                    details.rangeParms[2].fVal = std::stod(details.parameters[2]);                          // floating point inc

                                    try {
                                        details.rangeParms[1].iVal = std::stoi(details.parameters[1]);                      // integer count

                                        p_OptionsDescriptor.complexOptionValues[idx] = std::make_tuple(optionName, details); // reset values
                                    }
                                    catch (const std::out_of_range& e) {                                                    // not a valid integer
                                        errStr = ERR_MSG(ERROR::ARGUMENT_RANGE_COUNT_EXPECTED_INT) + std::string("'") + optionName  + std::string("'");
                                    }
                                }
                                catch (const std::out_of_range& e) {                                                        // not a valid floating point number
                                    errStr = ERR_MSG(ERROR::ARGUMENT_RANGE_PARMS_EXPECTED_FP) + std::string("'") + optionName  + std::string("'");
                                }
                            }
                        }
                    }
                    else {                                                                                                  // SET
                        // check for numeric/bool data types only that all set parameters are numeric/bool
                        // can't check for string data types 

                        if (dataType == TYPENAME::INT || dataType == TYPENAME::FLOAT) {                                     // numeric?
                            for (size_t ip = 0; ip < parms.size(); ip++) {                                                  // yes - for each set parameter specified
                                if ((dataType == TYPENAME::INT && !utils::IsINT(parms[ip])) ||                              // INT?
                                    (dataType == TYPENAME::FLOAT && !utils::IsFLOAT(parms[ip]))) {                          // FLOAT?
                                    error  = true;                                                                          // no - that's not ok
                                    errStr = ERR_MSG(ERROR::ARGUMENT_SET_EXPECTED_NUMERIC) + std::string("'") + optionName  + std::string("'");
                                    break;
                                }
                            }
                        }
                        else if (dataType == TYPENAME::BOOL) {                                                              // bool?
                            // we allow boolean set values to be specified as 1|0, T|F, TRUE|FALSE, Y|N, YES|NO, ON|OFF - but not mixed
                            // i.e. all must be 0|1, or all must be T|F, or all must be TRUE|FALSE etc. - case is not significant (downshifted already)

                            size_t check1, check2, check3, check4, check5, check6 = 0;
                            for (size_t ip = 0; ip < parms.size(); ip++) {                                                  // yes - for each set parameter specified
                                int check = std::abs(utils::IsBOOL(parms[ip]));                                             // check parm for valid BOOL
                                     if (check == 1) check1++;                                                              // valid: 0|1
                                else if (check == 2) check2++;                                                              // valid: t|f
                                else if (check == 3) check3++;                                                              // valid: true|false
                                else if (check == 4) check4++;                                                              // valid: y|n
                                else if (check == 5) check5++;                                                              // valid: yes|no
                                else if (check == 6) check6++;                                                              // valid: on|off
                                else                 break;                                                                 // not a valid BOOL
                            }

                            if (check1 != parms.size() && check2 != parms.size() && check3 != parms.size() && 
                                check4 != parms.size() && check5 != parms.size() && check6 != parms.size()) {               // all boolean, and consistent?
                                error  = true;                                                                              // no - that's not ok
                                errStr = ERR_MSG(ERROR::ARGUMENT_SET_EXPECTED_BOOL) + std::string("'") + optionName  + std::string("'");
                            }
                        }

                        if (!error) p_OptionsDescriptor.complexOptionValues[idx] = std::make_tuple(optionName, details);    // reset values

                    }
                }
                else {                                                                                                      // option not found in boost variables map
                    error  = true;                                                                                          // that can't be good...
                    errStr = ERR_MSG(ERROR::BOOST_OPTION_INTERNAL_ERROR) + std::string("'") + optionName + std::string("'");
                }

                if (error) break;                                                                                           // stop now
            }
    
            // ok, ranges and sets done - now manually set values for other options
            // we don't have to do so many checks here - by the time we get here much
            // has already been checked and validated

            if (errStr.empty()) {                                                                                           // no need if we've already flagged an error

                int complexCount = p_OptionsDescriptor.complexOptionValues.size();                                          // count of complex values (ranges or sets)

                if ((p_ArgCount - 1 - complexCount) > 0) {                                                                  // options other than those with ranges or sets specified?
                                                                                                                            // yes - process them
                    for (int iArg = 1; iArg < p_ArgCount; iArg ++) {                                                        // for each arg string (except the executable/placeholder)
                        std::string optionName(p_ArgStrings[iArg]);                                                         // get the option name

                        bool complex = false;                                                                               // initially
                        int iComplexArg = 0;                                                                                // first complex option
                        while (!complex && iComplexArg < complexCount) {                                                    // check all complex options
                            if (optionName == get<0>(p_OptionsDescriptor.complexOptionValues[iComplexArg])) {               // this option is complex option?
                                complex = true;                                                                             // yes
                            }
                            iComplexArg++;                                                                                  // next complex option
                        }
                        
                        if (!complex) {                                                                                     // this option complex?

                            if (optionName[0] == '-') optionName.erase(0, optionName.find_first_not_of("-"));               // remove the "-" or "--"
                                                                                                                            // no - replace the value
                            std::string optionValue = "";
                            if (iArg + 1 < p_ArgCount) optionValue = std::string(p_ArgStrings[iArg + 1]);                   // get the (potential) option value string

                            po::variables_map::const_iterator it = p_OptionsDescriptor.optionValues.m_VM.find(optionName);  // find the option in the boost variables map
                            if (it != p_OptionsDescriptor.optionValues.m_VM.end()) {                                        // found?

                                TYPENAME dataType = TYPENAME::NONE;                                                         // yes
                                std::tie(dataType, std::ignore, std::ignore, std::ignore) = OptionAttributes(p_OptionsDescriptor.optionValues.m_VM, it); // data type

                                // don't need to wrap the conversions in try/catch here - already validated
                                switch (dataType) {                                                                         // which data type?

                                    case TYPENAME::INT: {                                                                   // INT
                                        int thisVal = std::stoi(optionValue);                                               // convert to int
                                        p_OptionsDescriptor.optionValues.ModifyVariableMap(p_OptionsDescriptor.optionValues.m_VM, optionName, thisVal); // update Boost variable map
                                        po::notify(p_OptionsDescriptor.optionValues.m_VM);                                  // propagate the change
                                    } break;

                                    case TYPENAME::FLOAT: {                                                                 // FLOAT
                                        double thisVal = std::stod(optionValue);                                            // convert to double                   
                                        p_OptionsDescriptor.optionValues.ModifyVariableMap(p_OptionsDescriptor.optionValues.m_VM, optionName, thisVal); // update Boost variable map
                                        po::notify(p_OptionsDescriptor.optionValues.m_VM);                                  // propagate the change
                                    } break;

                                    case TYPENAME::BOOL: {                                                                  // BOOL
                                        // convert option value to bool
                                        // boolean options can be specified with no value, in which case the value is
                                        // assumed to be true (specifying the option with no value means enable it)
                                        // first check whether value starts with '--', then '-*', where * is any alpha
                                        // if either of those is true assume the value is an option string, so no value
                                        // was specified.
                                        // if neither is true, try to convert the value to a boolean.

                                        bool haveValue = true;                                                              // whether the user specified a value - default true
                                        if (optionValue.length() > 1 && optionValue[0] == '-') {                            // check for '--' or '-*' (* is alpga character)
                                            if (optionValue[1] == '-') haveValue = false;                                   // '--' - no value
                                            else if (isalpha(optionValue[1])) haveValue = false;                            // '-*', where * is alpha character - no value
                                        }

                                        bool thisVal = true;                                                                // default
                                        if (haveValue) {                                                                    // user specified a value?
                                            int boolVal = utils::IsBOOL(optionValue);                                       // yes - try to convert it to a boolean
                                            if (boolVal == 0) COMPLAIN(ERR_MSG(ERROR::INVALID_VALUE_FOR_BOOLEAN_OPTION));   // not a boolean - complain
                                            thisVal = boolVal < 0 ? false : true;                                           // is a boolean - set value
                                        }
                                        p_OptionsDescriptor.optionValues.ModifyVariableMap(p_OptionsDescriptor.optionValues.m_VM, optionName, thisVal); // update Boost variable map
                                        po::notify(p_OptionsDescriptor.optionValues.m_VM);                                  // propagate the change
                                    } break;

                                    case TYPENAME::STRING: {                                                                // STRING
                                        p_OptionsDescriptor.optionValues.ModifyVariableMap(p_OptionsDescriptor.optionValues.m_VM, optionName, optionValue); // update Boost variable map
                                        po::notify(p_OptionsDescriptor.optionValues.m_VM);                                  // propagate the change
                                    } break;
                                
                                    default:                                                                                // that's a problem...
                                        COMPLAIN(ERR_MSG(ERROR::INVALID_DATA_TYPE));                                        // complain
                                }            
                            }
                        }
                    }
                }
            }
        }
    }
    catch (po::error& e) {                                                                                                  // program options exception
        errStr = e.what();                                                                                                  // set error string
    }
    catch (const std::string eStr) {                                                                                        // custom exception
        errStr = eStr;                                                                                                      // set error string
    }
    catch (...) {                                                                                                           // unhandled exception
        errStr = ERR_MSG(ERROR::UNHANDLED_EXCEPTION);                                                                       // set error string
    }

    return errStr;
}


/*
 * Advances the command line options to the next variation.  A "variation" is
 * a combination of options defined by the option values, ranges, and sets
 * the user specified (in this case) on the command line.
 * 
 * When the command line is parsed initially (in Initialise()), the option values 
 * are set to the initial variation of range and set values, then once that is 
 * processed this function is called to advance to the next variation - the ranges 
 * and sets are played out in order.
 * 
 * 
 * int AdvanceCmdLineOptionValues()
 * 
 * @return                                      Int result:
 *                                                  -1: en error occurred
 *                                                   0: no more variations - all done
 *                                                   1: new variation applied - option values are set 
 */
int Options::AdvanceCmdLineOptionValues() {

    int retVal = 0;

    if (m_CmdLine.complexOptionValues.size() == 0) return retVal;

    // upon entry iterators for ranges and sets need to be advanced in order
    // to pick up the correct values to be loaded into the options.  Really
    // all we need to do is pick up the inner (fastest-counting) iterator (the
    // right-most in terms of placement on the commandline).  If we have to 
    // wrap-around to get that value, then we have to increment the next outer
    // (immediately left) iterator.
    
    // traverse each of the complex option values and gather the option values
    // we only need values for the options that have changing values
    // values have already been sanity checked by the time we get here
    // we traverse the complex options values in reverse order because the
    // fastest change is to the right...

    int idx   = m_CmdLine.complexOptionValues.size() - 1;                                                       // start at innermost iterator
    bool stop = false;                                                                                          // stop once we have all the values we need
    while (!stop && idx >= 0) {

        stop = true;                                                                                            // assume we have what we need

        std::string optionName        = get<0>(m_CmdLine.complexOptionValues[idx]);        
        RangeOrSetDescriptorT details = get<1>(m_CmdLine.complexOptionValues[idx]);

        details.currPos++;                                                                                      // advance iterator

        if (details.type == COMPLEX_TYPE::SET) {                                                                // SET
            if (details.currPos >= int(details.parameters.size())) {                                            // currPos is set position - wrap?
                if (idx == 0) {                                                                                 // outermost iterator?
                    retVal = 0;                                                                                 // yes - we're done
                    break;
                }
                else {                                                                                          // no - wrap
                    details.currPos = 0;                                                                        // ... back to zero
                    stop = false;                                                                               // ... and don't stop at this iterator
                }
            }
            std::string optionValue = details.parameters[details.currPos];                                      // option value (as string)

            switch (details.dataType) {                                                                         // which data type?

                case TYPENAME::INT: {                                                                           // INT
                    int thisVal = std::stoi(optionValue);
                    m_CmdLine.optionValues.ModifyVariableMap(m_CmdLine.optionValues.m_VM, optionName, thisVal);
                    po::notify(m_CmdLine.optionValues.m_VM);
                }   break;

                case TYPENAME::FLOAT: {                                                                         // FLOAT
                    double thisVal = std::stod(optionValue);                    
                    m_CmdLine.optionValues.ModifyVariableMap(m_CmdLine.optionValues.m_VM, optionName, thisVal);
                    po::notify(m_CmdLine.optionValues.m_VM);
                }   break;

                case TYPENAME::BOOL: {                                                                          // BOOL
                    bool thisVal = (optionValue == "1" || optionValue == "t" || optionValue == "true") ? true : false;
                    m_CmdLine.optionValues.ModifyVariableMap(m_CmdLine.optionValues.m_VM, optionName, thisVal);
                    po::notify(m_CmdLine.optionValues.m_VM);
                }   break;

                case TYPENAME::STRING: {                                                                        // STRING
                    m_CmdLine.optionValues.ModifyVariableMap(m_CmdLine.optionValues.m_VM, optionName, optionValue);
                    po::notify(m_CmdLine.optionValues.m_VM);
                }   break;
                                
                default: break;                                                                                 // already checked before we get here
            }            
        }
        else {                                                                                                  // RANGE
            if (details.currPos >= details.rangeParms[1].iVal) {                                                // currPos is range count - wrap?
                if (idx == 0) {                                                                                 // outermost iterator?
                    retVal = 0;                                                                                 // yes - we're done
                    break;
                }
                else {                                                                                          // no - wrap
                    details.currPos = 0;                                                                        // ... back to zero
                    stop = false;                                                                               // ... and don't stop at this iterator
                }
            }

            switch (details.dataType) {                                                                         // which data type?

                case TYPENAME::INT: {                                                                           // INT
                    int start = details.rangeParms[0].iVal;
                    int inc   = details.rangeParms[2].iVal;
                    int thisVal = start + (details.currPos * inc);
                    
                    m_CmdLine.optionValues.ModifyVariableMap(m_CmdLine.optionValues.m_VM, optionName, thisVal);
                    po::notify(m_CmdLine.optionValues.m_VM);
                }   break;

                case TYPENAME::FLOAT: {                                                                         // FLOAT
                    double start   = details.rangeParms[0].fVal;
                    double inc     = details.rangeParms[2].fVal;
                    double thisVal = start + (details.currPos * inc);
                    
                    m_CmdLine.optionValues.ModifyVariableMap(m_CmdLine.optionValues.m_VM, optionName, thisVal);
                    po::notify(m_CmdLine.optionValues.m_VM);
                }   break;

                default: break;                                                                                 // already checked before we get here
            }
        }

        m_CmdLine.complexOptionValues[idx] = std::make_tuple(optionName, details);                              // reset values

        bool skip = false;

        if (stop) {                                                                                             // do we think we're done?
                                                                                                                // yes, but should we skip this variation?
            if (m_CmdLine.optionValues.m_EvolutionMode == EVOLUTION_MODE::BSE) {                                // BSE?
                if (std::find(m_SSEOnly.begin(), m_SSEOnly.end(), optionName) != m_SSEOnly.end()) skip = true;  // skip SSE only option
            }
            else {                                                                                              // SSE
                if (std::find(m_BSEOnly.begin(), m_BSEOnly.end(), optionName) != m_BSEOnly.end()) skip = true;  // skip BSE only option
            }
        }

        if (!skip) {                                                                                            // skip this variation?
            retVal = 1;                                                                                         // no - set return value
            idx--;                                                                                              // next (outer) iterator
        }                                                                                            
    }

    return retVal;
}


/*
 * Initialise options service
 * 
 * Intitialises the options service.  Constructs options objects for the program options
 * (options that are specified only on the commandline and that cannot be specified in a 
 * grid file on a per object (star/binary) basis), and the grid file options (options that
 * can be specified in a grid file on a per object (star/binary) basis).
 * 
 * Initialises the commandline (program-level) options object, and the grid line (evolving 
 * object-level) options object.  Populates the commandline (program-level) options object
 * from the commandline arguments passed to main() - this object stays static throughout the
 * life of the program.
 * 
 * 
 * bool Options::Initialise(int p_ArgCount, char *p_ArgStrings[])
 * 
 * @param   [IN]    p_ArgCount                  Integer number of args passed in p_ArgStrings
 * @param   [IN]    p_ArgStrings                Arg strings - 1 per option
 *                                              Note that the first arg string is ignored (expected to be program name)
 * @return                                      Boolean status (true = ok, false = problem) (this function displays any error string)
 */
bool Options::Initialise(int p_ArgCount, char *p_ArgStrings[]) {

    bool ok = true;                                                                                                 // status - unless something changes

    try {

        m_CmdLine.optionValues.Initialise();                                                                        // initialise option variables for program-level options
        m_GridLine.optionValues.Initialise();                                                                       // initialise option variables for evolving object-level options

        po::options_description programLevelOptions("Program Options");                                             // boost options descriptions object for program-level options
        ok = AddOptions(&m_CmdLine.optionValues, &programLevelOptions);                                             // ... add
        if (!ok) {                                                                                                  // ok?
            COMPLAIN(ERR_MSG(ERROR::BOOST_OPTION_CMDLINE));                                                         // no, complain - this throws an exception
        }
        else {                                                                                                      // yes, ok

            m_CmdLine.optionDescriptions.add(programLevelOptions);                                                  // commandline options - stays static throughout the life of the program
    
            // we parse the option values before handing them over to boost
            // boost knows nothing about ranges and sets, so we have to handlde
            // them ourselves first
            m_CmdLine.complexOptionValues = {};                                                                     // no ranges or sets - unless we find them in the parse
            std::string errStr = ParseOptionValues(p_ArgCount, p_ArgStrings, m_CmdLine);                            // parse & populate the option values - specifically for ranges and sets
            if (!errStr.empty()) {                                                                                  // parsed ok?
                COMPLAIN(errStr);                                                                                   // no, complain - this throws an exception
            }
            else {
                errStr = m_CmdLine.optionValues.CheckAndSetOptions();                                               // yes - sanity check, and set, program-level values
                if (!errStr.empty()) {                                                                              // check ok?
                    COMPLAIN(errStr);                                                                               // no, complain - this throws an exception
                }
                else {

                    m_CmdLineOptionsDetails = OptionDetails(m_CmdLine);                                             // yes - get Run_Details contents

                    // initialise evolving object-level options.  The values of options specified in a 
                    // grid file take precedence over the values of the same options specified on the 
                    // commandline, but only for the object (star/binary) corresponding to the grid file
                    // record.

                    po::options_description objectLevelOptions("Program Options");                                  // boost options descriptions object for per object (star/binary) options
                    ok = AddOptions(&m_GridLine.optionValues, &objectLevelOptions);                                 // ... add
                    if (!ok) {                                                                                      // ok?
                        COMPLAIN(ERR_MSG(ERROR::BOOST_OPTION_GRIDLINE));                                            // no, complain - this throws an exception
                    }
                    else {                                                                                          // yes, ok
                        m_GridLine.optionDescriptions.add(objectLevelOptions);                                      // grid line options - stays static throughout the life of the program
                    }

                    // We now have the options the user entered at the commandline, including any ranges and/or 
                    // sets, so this is where we stop the initialisation - from here we just play out the options 
                    // that are specified by any ranges and sets via the AdvanceCmdLineOptionValues() function.
                    //
                    // If the user has specified any ranges or sets we set the options to the first value in each
                    // range or set (already done by the time we get here).  Calls to AdvanceCmdLineOptionValues()
                    // will then advance the option values through the ranges and sets as required - *however*, the
                    // very first call to AdvanceCmdLineOptionValues() just returns the options as they are set here,
                    // so that a call to AdvanceCmdLineOptionValues() can be put in a loop, and the first evaluation
                    // of the loop will be the first star/binary - and if only 1 star/binary is required then no loop,
                    // and no call to AdvanceCmdLineOptionValues() is required (because the initial values have already
                    // been set).
                    //
                    // Note that there are analogous functions for object (star/binary) initialisation and retrievel
                    // option values: InitialiseEvolvingObject() and AdvanceGridLineOptionValues().  These functions 
                    // intitialise and retrieve options specified in grid file records.
                }
            }
        }  
    }
    catch (po::error& e) {                                                                                          // program options exception
        std::cerr << ERR_MSG(ERROR::PROGRAM_OPTIONS_ERROR) << ": " << e.what() << std::endl;                        // show the problem
        std::cerr << m_CmdLine.optionDescriptions << std::endl;                                                     // show help
        ok = false;                                                                                                 // set status
    } 
    catch (const std::string eStr) {                                                                                // custom exception
        std::cerr << ERR_MSG(ERROR::PROGRAM_OPTIONS_ERROR) << ": " << eStr << std::endl;                            // show the problem
        std::cerr << m_CmdLine.optionDescriptions << std::endl;                                                     // show help
        ok = false;                                                                                                 // set status
    }
    catch (...) {                                                                                                   // unhandled exception
        std::cerr << ERR_MSG(ERROR::PROGRAM_OPTIONS_ERROR) << ": " << ERR_MSG(ERROR::UNHANDLED_EXCEPTION) << std::endl; // show the problem
        std::cerr << m_CmdLine.optionDescriptions << std::endl;                                                     // show help
        ok = false;                                                                                                 // set status
    }

    m_CmdLine.optionValues.m_Populated = ok;                                                                        // flag use

    return ok;
}


/*
 * Advances the grid line options to the next variation.  A "variation" is
 * a combination of options defined by the option values, ranges, and sets
 * the user specified (in this case) on the grid file line.
 * 
 * When the grid file line is read and parsed, the option values are set to
 * the initial variation of range and set values, then once that is processed
 * this function is called to advance to the next variation - the ranges and
 * sets are played out in order.
 * 
 * 
 * int AdvanceGridLineOptionValues()
 * 
 * @return                                      Int result:
 *                                                  -1: en error occurred
 *                                                   0: no more variations - all done
 *                                                   1: new variation applied - option values are set
 */
int Options::AdvanceGridLineOptionValues() {

    int retVal = 0;

    if (m_GridLine.complexOptionValues.size() == 0) return retVal;

    // upon entry iterators for ranges and sets need to be advanced in order
    // to pick up the correct values to be loaded into the options.  Really
    // all we need to do is pick up the inner (fastest-counting) iterator (the
    // right-most in terms of placement on the commandline).  If we have to 
    // wrap-around to get that value, then we have to increment the next outer
    // (immediately left) iterator.
    
    // traverse each of the complex option values and gather the option values
    // we only need values for the options that have changing values
    // values have already been sanity checked by the time we get here
    // we traverse the complex options values in reverse order because the
    // fastest change is to the right...

    bool stop  = false;                                                 // stop once we have all the values we need
    size_t idx = m_GridLine.complexOptionValues.size() - 1;
    while (!stop && idx >= 0) {

        stop = true;                                                    // assume we have what we need

        std::string optionName        = get<0>(m_GridLine.complexOptionValues[idx]);        
        RangeOrSetDescriptorT details = get<1>(m_GridLine.complexOptionValues[idx]);

        details.currPos++;                                              // advance iterator

        if (details.type == COMPLEX_TYPE::SET) {                        // SET
            if (details.currPos >= int(details.parameters.size())) {    // currPos is set position - wrap?
                if (idx == 0) {                                         // outermost iterator?
                    retVal = 0;                                         // yes - we're done
                    break;
                }
                else {                                                  // no - wrap
                    details.currPos = 0;                                // ... back to zero
                    stop = false;                                       // ... and don't stop at this iterator
                }
            }
            std::string optionValue = details.parameters[details.currPos];  // option value (as string)

            switch (details.dataType) {                                 // which data type?

                case TYPENAME::INT: {                                   // INT

                    int thisVal = std::stoi(optionValue);
                    m_GridLine.optionValues.ModifyVariableMap(m_GridLine.optionValues.m_VM, optionName, thisVal);
                    po::notify(m_GridLine.optionValues.m_VM);
                }   break;

                case TYPENAME::FLOAT: {                                 // FLOAT
                    double thisVal = std::stod(optionValue);                    
                    m_GridLine.optionValues.ModifyVariableMap(m_GridLine.optionValues.m_VM, optionName, thisVal);
                    po::notify(m_GridLine.optionValues.m_VM);
                }   break;

                case TYPENAME::BOOL: {                                  // BOOL
                    bool thisVal = (optionValue == "1" || optionValue == "t" || optionValue == "true") ? true : false;
                    m_CmdLine.optionValues.ModifyVariableMap(m_CmdLine.optionValues.m_VM, optionName, thisVal);
                    po::notify(m_CmdLine.optionValues.m_VM);
                }   break;

                case TYPENAME::STRING: {                                // STRING
                    m_CmdLine.optionValues.ModifyVariableMap(m_CmdLine.optionValues.m_VM, optionName, optionValue);
                    po::notify(m_CmdLine.optionValues.m_VM);
                }   break;
                                
                default: break;                                         // already checked before we get here
            }            
        }
        else {                                                          // RANGE

            if (details.currPos >= details.rangeParms[1].iVal) {        // currPos is range count - wrap?
                if (idx == 0) {                                         // outermost iterator?
                    retVal = 0;                                         // yes - we're done
                    break;
                }
                else {                                                  // no - wrap
                    details.currPos = 0;                                // ... back to zero
                    stop = false;                                       // ... and don't stop at this iterator
                }
            }

            switch (details.dataType) {                                 // which data type?

                case TYPENAME::INT: {                                   // INT

                    int start = details.rangeParms[0].iVal;
                    int inc   = details.rangeParms[2].iVal;
                    int thisVal = start + (details.currPos * inc);
                    
                    m_GridLine.optionValues.ModifyVariableMap(m_GridLine.optionValues.m_VM, optionName, thisVal);
                    po::notify(m_CmdLine.optionValues.m_VM);
                }   break;

                case TYPENAME::FLOAT: {                                 // FLOAT

                    double start   = details.rangeParms[0].fVal;
                    double inc     = details.rangeParms[2].fVal;
                    double thisVal = start + (details.currPos * inc);
                    
                    m_GridLine.optionValues.ModifyVariableMap(m_GridLine.optionValues.m_VM, optionName, thisVal);
                    po::notify(m_GridLine.optionValues.m_VM);
                }   break;

                default: break;                                         // already checked before we get here
            }

        }

        m_GridLine.complexOptionValues[idx] = std::make_tuple(optionName, details);  // reset values

        retVal = 1;
        idx--;                                                          // next (outer) iterator
    }

    return retVal;
}


/*
 * Initialise grid file options
 * 
 * Initialises and populates the grid line (evolving object-level) options object, using the
 * options the user specified in the grid file record - this object is updated for each grid
 * line read.
 *
 * 
 * bool Options::InitialiseEvolvingObject(int p_ArgCount, char *p_ArgStrings[])
 * 
 * @param   [IN]    p_OptionsString             String containing all options - the grid file record
 * @return                                      Boolean value indicating status: true = ok, false = an error occurred
 */
bool Options::InitialiseEvolvingObject(const std::string p_OptionsString) {

    bool ok = true;                                                                                                 // status - unless something changes

    try {

        // parse the option string (just as the OS/shell would do)

        std::vector<std::string> parsedStrings;                                                                     // parsed option strings

        size_t start      = 0;                                                                                      // start position of parsed option string
        size_t end        = 0;                                                                                      // end position of parsed option strinf
        std::string delim = " ";                                                                                    // delimiter
        while (end != std::string::npos) {                                                                          // iterate over input string
            end = p_OptionsString.find(delim, start);                                                               // find delimiter
            std::string s = p_OptionsString.substr(start, end - start);                                             // grab option/argument string
            parsedStrings.push_back(utils::trim(s));                                                                // trim whitespace and store
            start = end + delim.length();                                                                           // new start position
        }
    
        std::vector<char const *> args {"placeHolder"};                                                             // place-holder - boost expects command name as argv[0]
        for (auto& arg : parsedStrings)                                                                             // iterate over the parsed strings
            args.push_back(arg.c_str());                                                                            // and grab the c_str  

        // check here for excluded grid file options
        // if any are found, issue a warning and remove the option (and any
        // associated value) from the option strings
        // this should work here - I should be able to figure out what are
        // option names and what are option values - if it doesn't pan out
        // then we may need to move it to after Boost has parsed the options

        if (args.size() > 1) {                                                                                                  // args present (other than the executable/placeholder)
            std::vector<int> removeArgs = {};                                                                                   // vector of argument indices to be removed
            size_t iArg = 1;                                                                                                       // start after the executable name/placeholder
            while (iArg < args.size()) {                                                                   // for each arg string (except the executable/placeholder)

                std::string optionName(args[iArg]);                                                                             // get the string (we'll call it the option name for now)

                // check whether the string really is an option name
                // we assume any string starting with "--" is an option name, and
                // any string starting with "-*", where '*' is an alphabetic
                // character, is an option name
                bool haveOptionName = false;                                                              // is the string really an option name - default false
                if (optionName.length() > 1 && optionName[0] == '-') {                            // check for '--' or '-*' (* is alpha character)
                    if (optionName[1] == '-') haveOptionName = true;                                   // '--' - option name
                    else if (isalpha(optionName[1])) haveOptionName = true;                            // '-*', where * is alpha character - option name
                }
                iArg++;                                                                                 // next argument string

                if (haveOptionName) {                                                                   // do we think we have an option name?
                                                                                                        // yes
                    if (optionName[0] == '-') optionName.erase(0, optionName.find_first_not_of("-"));                               // remove the "-" or "--"

                    if (std::find(m_GridLineExcluded.begin(), m_GridLineExcluded.end(), optionName) != m_GridLineExcluded.end()) {  // on excluded list?
                        
                        removeArgs.push_back(iArg - 1);                                                                                 // yes - we need to remove it and any associated value

                        std::cerr << "WARNING: " << ERR_MSG(ERROR::OPTION_NOT_SUPPORTED_IN_GRID_FILE) << ": '" << optionName << "'\n";   // show warning

                        // we need to determine if the next argument string is a value for the
                        // option name we have, or whether it is the next option name - not all
                        // options need to specify a value (e.g. boolean switches)

                        std::string optionValue = "";
                        if (iArg < args.size()) optionValue = std::string(args[iArg]);                   // get the (potential) option value string

                        // check whether the string really is an option value
                        // as noted above, we assume any string starting with "--" is an option name,
                        // and any string starting with "-*", where '*' is an alphabetic character, 
                        // is an option name
                        bool haveOptionValue = false;                                                              // is the string really an option value - default false
                        if (!optionValue.empty()) {                                                             // empty string?
                            haveOptionValue = true;                                                              // no - we'll assume an option value, unless...
                            if (optionValue.length() > 1 && optionValue[0] == '-') {                            // check for '--' or '-*' (* is alpha character)
                                if (optionValue[1] == '-') haveOptionValue = false;                                   // '--' - option name, not value
                                else if (isalpha(optionValue[1])) haveOptionValue = false;                            // '-*', where * is alpha character - option name, not value
                            }
                        }

                        if (haveOptionValue) {                                                                   // do we think we have an option value?
                            removeArgs.push_back(iArg);                                                               // yes - we need to remove it
                            iArg++;                                                                                // next argument string
                        }
                    }
                }
            }

            // remove any argument strings identified as being excluded options
            // and the value associated with excluded options
            // this only works because the removeArgs vector is populated in the
            // same order as the args vector - that allows me to loop through args
            // in reverse order to delete elements (otherwise the index numbers
            // would change from under me as I was deleting them)

            int removeCount(removeArgs.size());
            if (removeCount > 0) {                                                            // anything to remove?
                for (int iArg = removeCount - 1; iArg >= 0; iArg--) {                           // loop through all args identified to be removed
                    args.erase(args.begin() + removeArgs[iArg]);                                                // erase it
                }
            }
        }

        // parse the option values before handing them over to boost
        // boost knows nothing about ranges and sets, so we have to handlde
        // them ourselves first
        m_GridLine.complexOptionValues = {};                                                                        // no ranges or sets - unless we find them in the parse
                                                    /*****  << do *not* try this at home >> *****/
        std::string errStr = ParseOptionValues(args.size(), const_cast<char**>(args.data()), m_GridLine);           // parse the option values - specifically for ranges and sets
        if (!errStr.empty()) {                                                                                      // parsed ok?
            COMPLAIN(errStr);                                                                                       // no, complain - this throws an exception
        }
        else {
            errStr = m_GridLine.optionValues.CheckAndSetOptions();                                                  // yes - sanity check, and set, evolving object-level values
            if (!errStr.empty()) {                                                                                  // check ok?
                COMPLAIN(errStr);                                                                                   // no, complain - this throws an exception
            }
        }
    }
    catch (po::error& e) {                                                                                          // program options exception
        std::cerr << ERR_MSG(ERROR::PROGRAM_OPTIONS_ERROR) << ": " << e.what() << std::endl;                        // show the problem
        std::cerr << m_GridLine.optionDescriptions << std::endl;                                                    // show help
        ok = false;                                                                                                 // set status
    } 
    catch (const std::string eStr) {                                                                                // custom exception
        std::cerr << ERR_MSG(ERROR::PROGRAM_OPTIONS_ERROR) << ": " << eStr << std::endl;                            // show the problem
        std::cerr << m_GridLine.optionDescriptions << std::endl;                                                    // show help
        ok = false;                                                                                                 // set status
    }
    catch (...) {                                                                                                   // unhandled exception
        std::cerr << ERR_MSG(ERROR::PROGRAM_OPTIONS_ERROR) << ": " << ERR_MSG(ERROR::UNHANDLED_EXCEPTION) << std::endl; // show the problem
        std::cerr << m_GridLine.optionDescriptions << std::endl;                                                    // show help
        ok = false;                                                                                                 // set status
    }

    m_GridLine.optionValues.m_Populated = ok;                                                                       // flag use
    
    return ok;
}


/*
 * Read and apply the next record in the grid file
 * 
 * The record from the grid file is read as one string, then passed to
 * InitialiseEvolvingObject() for processing.
 * 
 * In InitialiseEvolvingObject() the record is parsed into separate tokens
 * ready to be passed to the boost program option parser.  Once that is 
 * done the options are handed over to the boost functions for parsing and
 * detting of values.
 * 
 * pon return from this function the option values will be set to the values 
 * specified by the user, or their default values - either way ready for the 
 * star/binary to be evolved.
 * 
 * The grid file struct (m_Gridfile) will be used and updated by this function.
 * The grid file name, error status, and file handle are stored in the struct.  
 * 
 * 
 * int ApplyNextGridLine()
 * 
 * @return                                      Int result:
 *                                                  -1: Error reading grid file record (error value in grid file struct)
 *                                                   0: No record to read - end of file
 *                                                   1: Grid file record read and applied ok
 */
int Options::ApplyNextGridLine() {

    int status = -1;                                                // default status is failure

    if (m_Gridfile.handle.is_open()) {                              // file open?
                                                                    // yes
        std::string record;                                         // the record read
        std::getline(m_Gridfile.handle, record);                    // read the next record
        if (m_Gridfile.handle.fail()) {                             // read ok?
            if (m_Gridfile.handle.eof()) status = 0;                // eof?
            else {                                                  // no
                m_Gridfile.error = ERROR::FILE_READ_ERROR;          // record error
                status = -1;                                        // set status
            }
        }
        else {                                                      // read ok
            status = InitialiseEvolvingObject(record) ? 1 : -1;     // apply record and set status
        }
    }

    return status;
}


/*
 * Open the grid file
 *
 * The grid file is opened at the start of the simulation and stays open until the
 * simulation is complete - we just pick a record off and process the record, and
 * when we hit the end of the file the file is closed and the simulation complete.
 *
 * 
 * ERROR OpenGridFile(const std::string p_GridFilename)
 *
 * @param   [IN]        p_Filename              The filename of the Grid file
 * @return                                      ERROR indicator - will be ERROR::NONE if file opened sccessfully
 */
ERROR Options::OpenGridFile(const std::string p_GridFilename) {

    m_Gridfile.filename = p_GridFilename;                       // record filename

    if (!m_Gridfile.filename.empty()) {                         // have grid filename?
        m_Gridfile.handle.open(m_Gridfile.filename);            // yes - open the file
        if (m_Gridfile.handle.fail()) {                         // open ok?
            m_Gridfile.error = ERROR::FILE_OPEN_ERROR;          // no - record error
        }
        else m_Gridfile.error = ERROR::NONE;                    // open ok - no error
    }
    else m_Gridfile.error = ERROR::EMPTY_FILENAME;              // empty filename

    return m_Gridfile.error;
}


/*
 * Determine the value of the requested program option
 *
 * The property is a boost variant variable, and is one of the following types:
 *
 *      STAR_PROPERTY           - any individual star property
 *      STAR_1_PROPERTY         - property of the primary (m_Star1)
 *      STAR_2_PROPERTY         - property of the secondary (m_Star2)
 *      SUPERNOVA_PROPERTY      - property of the star that has gone supernova
 *      COMPANION_PROPERTY      - property of the companion to the supernova
 *      BINARY_PROPERTY         - property of the binary
 *      PROGRAM_OPTION          - program option
 *
 * This function handles properties of type PROGRAM_OPTION
 *
 * This is the function used to retrieve values for properties required to be printed.
 * This allows the composition of the log records to be dynamically modified - this is
 * how we allow users to specify what properties they want recorded in log files.
 *
 * The functional return is the value of the property requested.  The type of the
 * functional return is a tuple: std::tuple<bool, COMPAS_VARIABLE_TYPE>.  This type
 * is COMPAS_VARIABLE by typedef.
 *
 * The bool returned indicates whether the property value was retrieved ok: true = yes, fales = no
 * The COMPAS_VARIABLE_TYPE variable returned is a boost variant variable, the value of which is the
 * value of the underlying primitive variable.
 *
 *
 * COMPAS_VARIABLE OptionValue(const T_ANY_PROPERTY p_Property) const
 *
 * @param   [IN]    p_Property                  The property for which the value is required
 * @return                                      The value of the requested property
 */
COMPAS_VARIABLE Options::OptionValue(const T_ANY_PROPERTY p_Property) const {

    bool ok = true;                                                                                                     // status - unless a problem occurs

    COMPAS_VARIABLE_TYPE value;                                                                                         // default property value

    PROGRAM_OPTION property = boost::get<PROGRAM_OPTION>(p_Property);                                                   // get property
                                                                                                                        // get property value
    switch (property) {

        // Floor
        /*
        case PROGRAM_OPTION::AIS_DCO_TYPE                                   : value = static_cast<int>(AIS_DCOType());                                      break;
        case PROGRAM_OPTION::AIS_EXPLORATORY_PHASE                          : value = AIS_ExploratoryPhase();                                               break;
        case PROGRAM_OPTION::AIS_HUBBLE                                     : value = AIS_Hubble();                                                         break;
        case PROGRAM_OPTION::AIS_PESSIMISTIC                                : value = AIS_Pessimistic();                                                    break;
        case PROGRAM_OPTION::AIS_REFINEMENT_PHASE                           : value = AIS_RefinementPhase();                                                break;
        case PROGRAM_OPTION::AIS_RLOF                                       : value = AIS_RLOF();                                                           break;
        */

        case PROGRAM_OPTION::ALLOW_MS_STAR_TO_SURVIVE_COMMON_ENVELOPE       : value = AllowMainSequenceStarToSurviveCommonEnvelope();                       break;
        case PROGRAM_OPTION::ALLOW_RLOF_AT_BIRTH                            : value = AllowRLOFAtBirth();                                                   break;
        case PROGRAM_OPTION::ALLOW_TOUCHING_AT_BIRTH                        : value = AllowTouchingAtBirth();                                               break;
        case PROGRAM_OPTION::ANG_MOM_CONSERVATION_DURING_CIRCULARISATION    : value = AngularMomentumConservationDuringCircularisation();                   break;

        // Serena
        //case PROGRAM_OPTION::BE_BINARIES                                    : value = BeBinaries();                                                         break;

        case PROGRAM_OPTION::BLACK_HOLE_KICKS_OPTION                        : value = static_cast<int>(BlackHoleKicksOption());                             break;

        case PROGRAM_OPTION::EVOLUTION_MODE                                 : value = static_cast<int>(EvolutionMode());                                    break;
    
        case PROGRAM_OPTION::CASE_BB_STABILITY_PRESCRIPTION                 : value = static_cast<int>(CaseBBStabilityPrescription());                      break;
    
        case PROGRAM_OPTION::CHE_OPTION                                     : value = static_cast<int>(CHE_Option());                                       break;

        case PROGRAM_OPTION::CIRCULARISE_BINARY_DURING_MT                   : value = CirculariseBinaryDuringMassTransfer();                                break;

        case PROGRAM_OPTION::COMMON_ENVELOPE_ALPHA                          : value = CommonEnvelopeAlpha();                                                break;
        case PROGRAM_OPTION::COMMON_ENVELOPE_ALPHA_THERMAL                  : value = CommonEnvelopeAlphaThermal();                                         break;
        case PROGRAM_OPTION::COMMON_ENVELOPE_LAMBDA                         : value = CommonEnvelopeLambda();                                               break;
        case PROGRAM_OPTION::COMMON_ENVELOPE_LAMBDA_MULTIPLIER              : value = CommonEnvelopeLambdaMultiplier();                                     break;
        case PROGRAM_OPTION::COMMON_ENVELOPE_LAMBDA_PRESCRIPTION            : value = static_cast<int>(CommonEnvelopeLambdaPrescription());                 break;
        case PROGRAM_OPTION::COMMON_ENVELOPE_MASS_ACCRETION_CONSTANT        : value = CommonEnvelopeMassAccretionConstant();                                break;
        case PROGRAM_OPTION::COMMON_ENVELOPE_MASS_ACCRETION_MAX             : value = CommonEnvelopeMassAccretionMax();                                     break;
        case PROGRAM_OPTION::COMMON_ENVELOPE_MASS_ACCRETION_MIN             : value = CommonEnvelopeMassAccretionMin();                                     break;
        case PROGRAM_OPTION::COMMON_ENVELOPE_MASS_ACCRETION_PRESCRIPTION    : value = static_cast<int>(CommonEnvelopeMassAccretionPrescription());          break;
        case PROGRAM_OPTION::COMMON_ENVELOPE_RECOMBINATION_ENERGY_DENSITY   : value = CommonEnvelopeRecombinationEnergyDensity();                           break;
        case PROGRAM_OPTION::COMMON_ENVELOPE_SLOPE_KRUCKOW                  : value = CommonEnvelopeSlopeKruckow();                                         break;

        case PROGRAM_OPTION::ECCENTRICITY                                   : value = Eccentricity();                                                       break;
        case PROGRAM_OPTION::ECCENTRICITY_DISTRIBUTION                      : value = static_cast<int>(EccentricityDistribution());                         break;
        case PROGRAM_OPTION::ECCENTRICITY_DISTRIBUTION_MAX                  : value = EccentricityDistributionMax();                                        break;
        case PROGRAM_OPTION::ECCENTRICITY_DISTRIBUTION_MIN                  : value = EccentricityDistributionMin();                                        break;
        case PROGRAM_OPTION::EDDINGTON_ACCRETION_FACTOR                     : value = EddingtonAccretionFactor();                                           break;
        case PROGRAM_OPTION::ENVELOPE_STATE_PRESCRIPTION                    : value = static_cast<int>(EnvelopeStatePrescription());                        break;

        case PROGRAM_OPTION::FRYER_SUPERNOVA_ENGINE                         : value = static_cast<int>(FryerSupernovaEngine());                             break;

        case PROGRAM_OPTION::INITIAL_MASS                                   : value = InitialMass();                                                        break;
        case PROGRAM_OPTION::INITIAL_MASS_1                                 : value = InitialMass1();                                                       break;
        case PROGRAM_OPTION::INITIAL_MASS_2                                 : value = InitialMass2();                                                       break;

        case PROGRAM_OPTION::INITIAL_MASS_FUNCTION                          : value = static_cast<int>(InitialMassFunction());                              break;
        case PROGRAM_OPTION::INITIAL_MASS_FUNCTION_MAX                      : value = InitialMassFunctionMax();                                             break;
        case PROGRAM_OPTION::INITIAL_MASS_FUNCTION_MIN                      : value = InitialMassFunctionMin();                                             break;
        case PROGRAM_OPTION::INITIAL_MASS_FUNCTIONPOWER                     : value = InitialMassFunctionPower();                                           break;

        case PROGRAM_OPTION::KICK_DIRECTION_DISTRIBUTION                    : value = static_cast<int>(KickDirectionDistribution());                        break;
        case PROGRAM_OPTION::KICK_DIRECTION_POWER                           : value = KickDirectionPower();                                                 break;
        case PROGRAM_OPTION::KICK_SCALING_FACTOR                            : value = KickScalingFactor();                                                  break;
        case PROGRAM_OPTION::KICK_MAGNITUDE_DISTRIBUTION                    : value = static_cast<int>(KickMagnitudeDistribution());                        break;
        case PROGRAM_OPTION::KICK_MAGNITUDE_DISTRIBUTION_MAXIMUM            : value = KickMagnitudeDistributionMaximum();                                   break;

        case PROGRAM_OPTION::KICK_MAGNITUDE_DISTRIBUTION_SIGMA_CCSN_BH      : value = KickMagnitudeDistributionSigmaCCSN_BH();                              break;
        case PROGRAM_OPTION::KICK_MAGNITUDE_DISTRIBUTION_SIGMA_CCSN_NS      : value = KickMagnitudeDistributionSigmaCCSN_NS();                              break;
        case PROGRAM_OPTION::KICK_MAGNITUDE_DISTRIBUTION_SIGMA_FOR_ECSN     : value = KickMagnitudeDistributionSigmaForECSN();                              break;
        case PROGRAM_OPTION::KICK_MAGNITUDE_DISTRIBUTION_SIGMA_FOR_USSN     : value = KickMagnitudeDistributionSigmaForUSSN();                              break;

        case PROGRAM_OPTION::KICK_MAGNITUDE                                 : value = KickMagnitude();                                                      break;
        case PROGRAM_OPTION::KICK_MAGNITUDE_1                               : value = KickMagnitude1();                                                     break;
        case PROGRAM_OPTION::KICK_MAGNITUDE_2                               : value = KickMagnitude2();                                                     break;

        case PROGRAM_OPTION::KICK_MAGNITUDE_RANDOM                          : value = KickMagnitudeRandom();                                                break;
        case PROGRAM_OPTION::KICK_MAGNITUDE_RANDOM_1                        : value = KickMagnitudeRandom1();                                               break;
        case PROGRAM_OPTION::KICK_MAGNITUDE_RANDOM_2                        : value = KickMagnitudeRandom2();                                               break;

        case PROGRAM_OPTION::KICK_MEAN_ANOMALY_1                            : value = SN_MeanAnomaly1();                                                    break;
        case PROGRAM_OPTION::KICK_MEAN_ANOMALY_2                            : value = SN_MeanAnomaly2();                                                    break;
        case PROGRAM_OPTION::KICK_PHI_1                                     : value = SN_Phi1();                                                            break;
        case PROGRAM_OPTION::KICK_PHI_2                                     : value = SN_Phi2();                                                            break;
        case PROGRAM_OPTION::KICK_THETA_1                                   : value = SN_Theta1();                                                          break;
        case PROGRAM_OPTION::KICK_THETA_2                                   : value = SN_Theta2();                                                          break;

        case PROGRAM_OPTION::LBV_FACTOR                                     : value = LuminousBlueVariableFactor();                                         break;

        case PROGRAM_OPTION::MASS_LOSS_PRESCRIPTION                         : value = static_cast<int>(MassLossPrescription());                             break;

        case PROGRAM_OPTION::MASS_RATIO_DISTRIBUTION                        : value = static_cast<int>(MassRatioDistribution());                            break;
        case PROGRAM_OPTION::MASS_RATIO_DISTRIBUTION_MAX                    : value = MassRatioDistributionMax();                                           break;
        case PROGRAM_OPTION::MASS_RATIO_DISTRIBUTION_MIN                    : value = MassRatioDistributionMin();                                           break;

        case PROGRAM_OPTION::MT_ACCRETION_EFFICIENCY_PRESCRIPTION           : value = static_cast<int>(MassTransferAccretionEfficiencyPrescription());      break;
        case PROGRAM_OPTION::MT_ANG_MOM_LOSS_PRESCRIPTION                   : value = static_cast<int>(MassTransferAngularMomentumLossPrescription());      break;
        case PROGRAM_OPTION::MT_THERMAL_LIMIT_C                             : value = MassTransferCParameter();                                             break;

        // AVG
        /*
        case PROGRAM_OPTION::MT_CRIT_MR_MS_LOW_MASS                         : value = MassTransferCriticalMassRatioMSLowMass();                             break;
        case PROGRAM_OPTION::MT_CRIT_MR_MS_LOW_MASS_DEGENERATE_ACCRETOR     : value = MassTransferCriticalMassRatioMSLowMassDegenerateAccretor();           break;
        case PROGRAM_OPTION::MT_CRIT_MR_MS_LOW_MASS_NON_DEGENERATE_ACCRETOR : value = MassTransferCriticalMassRatioMSLowMassNonDegenerateAccretor();        break;
        case PROGRAM_OPTION::MT_CRIT_MR_MS_HIGH_MASS                        : value = MassTransferCriticalMassRatioMSHighMass();                            break;
        case PROGRAM_OPTION::MT_CRIT_MR_MS_HIGH_MASS_DEGENERATE_ACCRETOR    : value = MassTransferCriticalMassRatioMSHighMassDegenerateAccretor();          break;
        case PROGRAM_OPTION::MT_CRIT_MR_MS_HIGH_MASS_NON_DEGENERATE_ACCRETOR: value = MassTransferCriticalMassRatioMSHighMassNonDegenerateAccretor();       break;
        case PROGRAM_OPTION::MT_CRIT_MR_GIANT                               : value = MassTransferCriticalMassRatioGiant();                                 break;
        case PROGRAM_OPTION::MT_CRIT_MR_GIANT_DEGENERATE_ACCRETOR           : value = MassTransferCriticalMassRatioGiantDegenerateAccretor();               break;
        case PROGRAM_OPTION::MT_CRIT_MR_GIANT_NON_DEGENERATE_ACCRETOR       : value = MassTransferCriticalMassRatioGiantNonDegenerateAccretor();            break;
        case PROGRAM_OPTION::MT_CRIT_MR_HG                                  : value = MassTransferCriticalMassRatioHG();                                    break;
        case PROGRAM_OPTION::MT_CRIT_MR_HG_DEGENERATE_ACCRETOR              : value = MassTransferCriticalMassRatioHGDegenerateAccretor();                  break;
        case PROGRAM_OPTION::MT_CRIT_MR_HG_NON_DEGENERATE_ACCRETOR          : value = MassTransferCriticalMassRatioHGNonDegenerateAccretor();               break;
        case PROGRAM_OPTION::MT_CRIT_MR_HE_GIANT                            : value = MassTransferCriticalMassRatioHeliumGiant();                           break;
        case PROGRAM_OPTION::MT_CRIT_MR_HE_GIANT_DEGENERATE_ACCRETOR        : value = MassTransferCriticalMassRatioHeliumGiantDegenerateAccretor();         break;
        case PROGRAM_OPTION::MT_CRIT_MR_HE_GIANT_NON_DEGENERATE_ACCRETOR    : value = MassTransferCriticalMassRatioHeliumGiantNonDegenerateAccretor();      break;
        case PROGRAM_OPTION::MT_CRIT_MR_HE_HG                               : value = MassTransferCriticalMassRatioHeliumHG();                              break;
        case PROGRAM_OPTION::MT_CRIT_MR_HE_HG_DEGENERATE_ACCRETOR           : value = MassTransferCriticalMassRatioHeliumHGDegenerateAccretor();            break;
        case PROGRAM_OPTION::MT_CRIT_MR_HE_HG_NON_DEGENERATE_ACCRETOR       : value = MassTransferCriticalMassRatioHeliumHGNonDegenerateAccretor();         break;
        case PROGRAM_OPTION::MT_CRIT_MR_HE_MS                               : value = MassTransferCriticalMassRatioHeliumMS();                              break;
        case PROGRAM_OPTION::MT_CRIT_MR_HE_MS_DEGENERATE_ACCRETOR           : value = MassTransferCriticalMassRatioHeliumMSDegenerateAccretor();            break;
        case PROGRAM_OPTION::MT_CRIT_MR_HE_MS_NON_DEGENERATE_ACCRETOR       : value = MassTransferCriticalMassRatioHeliumMSNonDegenerateAccretor();         break;
        case PROGRAM_OPTION::MT_CRIT_MR_WD                                  : value = MassTransferCriticalMassRatioWhiteDwarf();                            break;
        case PROGRAM_OPTION::MT_CRIT_MR_WD_DEGENERATE_ACCRETOR              : value = MassTransferCriticalMassRatioWhiteDwarfDegenerateAccretor();          break;
        case PROGRAM_OPTION::MT_CRIT_MR_WD_NONDEGENERATE_ACCRETOR           : value = MassTransferCriticalMassRatioWhiteDwarfNonDegenerateAccretor();       break;
        */

        case PROGRAM_OPTION::MT_FRACTION_ACCRETED                           : value = MassTransferFractionAccreted();                                       break;
        case PROGRAM_OPTION::MT_JLOSS                                       : value = MassTransferJloss();                                                  break;
        case PROGRAM_OPTION::MT_REJUVENATION_PRESCRIPTION                   : value = static_cast<int>(MassTransferRejuvenationPrescription());             break;
        case PROGRAM_OPTION::MT_THERMALLY_LIMITED_VARIATION                 : value = static_cast<int>(MassTransferThermallyLimitedVariation());            break;

        case PROGRAM_OPTION::MCBUR1                                         : value = MCBUR1();                                                             break;

        case PROGRAM_OPTION::METALLICITY                                    : value = Metallicity();                                                        break;

        case PROGRAM_OPTION::MINIMUM_MASS_SECONDARY                         : value = MinimumMassSecondary();                                               break;
        case PROGRAM_OPTION::MAXIMUM_NEUTRON_STAR_MASS                      : value = MaximumNeutronStarMass();                                             break;

        case PROGRAM_OPTION::NEUTRINO_MASS_LOSS_ASSUMPTION_BH               : value = static_cast<int>(NeutrinoMassLossAssumptionBH());                     break;
        case PROGRAM_OPTION::NEUTRINO_MASS_LOSS_VALUE_BH                    : value = NeutrinoMassLossValueBH();                                            break;

        case PROGRAM_OPTION::NS_EOS                                         : value = static_cast<int>(NeutronStarEquationOfState());                       break;

        case PROGRAM_OPTION::PISN_LOWER_LIMIT                               : value = PairInstabilityLowerLimit();                                          break;
        case PROGRAM_OPTION::PISN_UPPER_LIMIT                               : value = PairInstabilityUpperLimit();                                          break;

        case PROGRAM_OPTION::PERIOD_DISTRIBUTION_MAX                        : value = PeriodDistributionMax();                                              break;
        case PROGRAM_OPTION::PERIOD_DISTRIBUTION_MIN                        : value = PeriodDistributionMin();                                              break;

        case PROGRAM_OPTION::PULSAR_MAGNETIC_FIELD_DISTRIBUTION             : value = static_cast<int>(PulsarBirthMagneticFieldDistribution());             break;
        case PROGRAM_OPTION::PULSAR_MAGNETIC_FIELD_DISTRIBUTION_MAX         : value = PulsarBirthMagneticFieldDistributionMax();                            break;
        case PROGRAM_OPTION::PULSAR_MAGNETIC_FIELD_DISTRIBUTION_MIN         : value = PulsarBirthMagneticFieldDistributionMin();                            break;

        case PROGRAM_OPTION::PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION          : value = static_cast<int>(PulsarBirthSpinPeriodDistribution());                break;
        case PROGRAM_OPTION::PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION_MAX      : value = PulsarBirthSpinPeriodDistributionMax();                               break;
        case PROGRAM_OPTION::PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION_MIN      : value = PulsarBirthSpinPeriodDistributionMin();                               break;

        case PROGRAM_OPTION::PULSAR_LOG10_MINIMUM_MAGNETIC_FIELD            : value = PulsarLog10MinimumMagneticField();                                    break;

        case PROGRAM_OPTION::PULSAR_MAGNETIC_FIELD_DECAY_MASS_SCALE         : value = PulsarMagneticFieldDecayMassscale();                                  break;
        case PROGRAM_OPTION::PULSAR_MAGNETIC_FIELD_DECAY_TIME_SCALE         : value = PulsarMagneticFieldDecayTimescale();                                  break;

        case PROGRAM_OPTION::PPI_PRESCRIPTION                               : value = static_cast<int>(PulsationalPairInstabilityPrescription());           break;
        case PROGRAM_OPTION::PPI_LOWER_LIMIT                                : value = PulsationalPairInstabilityLowerLimit();                               break;
        case PROGRAM_OPTION::PPI_UPPER_LIMIT                                : value = PulsationalPairInstabilityUpperLimit();                               break;

        case PROGRAM_OPTION::RANDOM_SEED                                    : value = RandomSeed();                                                         break;
        case PROGRAM_OPTION::RANDOM_SEED_CMDLINE                            : value = RandomSeedCmdLine();                                                  break;

        case PROGRAM_OPTION::REMNANT_MASS_PRESCRIPTION                      : value = static_cast<int>(RemnantMassPrescription());                          break;

        case PROGRAM_OPTION::ROTATIONAL_VELOCITY_DISTRIBUTION               : value = static_cast<int>(RotationalVelocityDistribution());                   break;
   
        case PROGRAM_OPTION::SEMI_MAJOR_AXIS                                : value = SemiMajorAxis();                                                      break;
        case PROGRAM_OPTION::SEMI_MAJOR_AXIS_DISTRIBUTION                   : value = static_cast<int>(SemiMajorAxisDistribution());                        break;
        case PROGRAM_OPTION::SEMI_MAJOR_AXIS_DISTRIBUTION_MAX               : value = SemiMajorAxisDistributionMax();                                       break;
        case PROGRAM_OPTION::SEMI_MAJOR_AXIS_DISTRIBUTION_MIN               : value = SemiMajorAxisDistributionMin();                                       break;
        case PROGRAM_OPTION::SEMI_MAJOR_AXIS_DISTRIBUTION_POWER             : value = SemiMajorAxisDistributionPower();                                     break;

        case PROGRAM_OPTION::STELLAR_ZETA_PRESCRIPTION                      : value = static_cast<int>(StellarZetaPrescription());                          break;

        case PROGRAM_OPTION::WR_FACTOR                                      : value = WolfRayetFactor();                                                    break;

        case PROGRAM_OPTION::ZETA_RADIATIVE_ENVELOPE_GIANT                  : value = ZetaRadiativeEnvelopeGiant();                                         break;
        case PROGRAM_OPTION::ZETA_MS                                        : value = ZetaMainSequence();                                                   break;
        case PROGRAM_OPTION::ZETA_ADIABATIC_ARBITRARY                       : value = ZetaAdiabaticArbitrary();                                             break;

        default:                                                                                                        // unknown property
            ok    = false;                                                                                              // that's not ok...
            value = "UNKNOWN";                                                                                          // default value
            std::cerr << ERR_MSG(ERROR::UNKNOWN_PROGRAM_OPTION) << std::endl;                                           // show warning (don't have logging or errors here...)
    }

    return std::make_tuple(ok, value);
}
