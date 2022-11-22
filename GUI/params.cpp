#include <qstring.h>
#include "params.h"
#include "QDebug"

Params::Params()
{
    static infoStruct label_info[] = {
        {"PARENT_0", "Diffusion coefficient within the blob."},
        {"PARENT_1", "Diffusion coefficient in the medium."},
        {"PARENT_2", "Cell influx parameter Kin.  The rate of mass transport into the cell is Kin.Cex - Kout.Cin (currently no dependence on cell surface area)."},
        {"PARENT_3", "Cell efflux parameter Kout.  The rate of mass transport into the cell is Kin.Cex - Kout.Cin (currently no dependence on cell surface area)."},
        {"PARENT_4", "Half-life of the compound, used to calculate the decay rate.  This is the same in the cell and in the medium."},
        {"PARENT_CT1_0", "Kmet0 is the maximum rate of metabolism.  The actual rate is the product of drug concentration Cdrug, Kmet0 and a sigmoid function of O2 concentration C_O2, with parameters C2 and KO2:\n\
metabolism rate = dMdt = Cdrug.(1 - C2 + C2.KO2^n_O2/(KO2^n_O2 + C_O2^n_O2)).Kmet0   If Vmax > 0, Kmet0 is replaced by (Kmet0 + Vmax*Cdrug/(Km + Cdrug)) \n\
This the rate of transformation of parent drug to metabolite 1, or of metabolite 1 to metabolite 2, or removal of metabolite 2"},
        {"PARENT_CT1_1", "C2 is one of the three parameters of the basic sigmoid function of C_O2 that determines metabolism rate.  When C_O2 = 0, the function = 1, when C_O2 >> KO2, the function = 1 - C2: \n\
metabolism rate = dMdt = Cdrug.(1 - C2 + C2.KO2^n_O2/(KO2^n_O2 + C_O2^n_O2)).Kmet0   If Vmax > 0, Kmet0 is replaced by (Kmet0 + Vmax*Cdrug/(Km + Cdrug)) \n\
This the rate of transformation of parent drug to metabolite 1, or of metabolite 1 to metabolite 2, or removal of metabolite 2"},
        {"PARENT_CT1_2",  "KO2 is one of the three parameters of the basic sigmoid function of C_O2 that determines metabolism rate: \n\
metabolism rate = dMdt = Cdrug.(1 - C2 + C2.KO2^n_O2/(KO2^n_O2 + C_O2^n_O2)).Kmet0   If Vmax > 0, Kmet0 is replaced by (Kmet0 + Vmax*Cdrug/(Km + Cdrug)) \n\
This the rate of transformation of parent drug to metabolite 1, or of metabolite 1 to metabolite 2, or removal of metabolite 2"},
        {"PARENT_CT1_3", "Vmax and Km are parameters that determine the dependence of the maximum rate of metabolism on drug concentration: \n\
metabolism rate = dMdt = Cdrug.(1 - C2 + C2.KO2^n_O2/(KO2^n_O2 + C_O2^n_O2)).Kmet0   If Vmax > 0, Kmet0 is replaced by (Kmet0 + Vmax*Cdrug/(Km + Cdrug))"},
         {"PARENT_CT1_4", "Vmax and Km are parameters that determine the dependence of the maximum rate of metabolism on drug concentration: \n\
metabolism rate = dMdt = Cdrug.(1 - C2 + C2.KO2^n_O2/(KO2^n_O2 + C_O2^n_O2)).Kmet0   If Vmax > 0, Kmet0 is replaced by (Kmet0 + Vmax*Cdrug/(Km + Cdrug))"},
        {"PARENT_CT1_5", "Klesion is currently unused."},
        {"PARENT_CT1_6", "The O2 concentration in the kill experiment."},
        {"PARENT_CT1_7", "The drug concentration in the kill experiment."},
        {"PARENT_CT1_8", "The duration the kill experiment."},
        {"PARENT_CT1_9", "The kill fraction achieved in the kill experiment (1 - SF)."},
        {"PARENT_CT1_10", "Sensitisation of the cells to radiation is determined by three parameters.  The usual radiation kill parameters OER_alpha and OER_beta are multiplied by the sensitisation enhancement ratio SER: \n\
         SER = (C_O2 + SER_KO2*(Cdrug*SER_max + SER_Km)/(Cdrug + SER_Km))/(C_O2 + SER_KO2)\n\
         If a DNA-damage inhibiter, this is a_inhibit"},
        {"PARENT_CT1_11", "Sensitisation of the cells to radiation is determined by three parameters.  The usual radiation kill parameters OER_alpha and OER_beta are multiplied by the sensitisation enhancement ratio SER: \n\
         SER = (C_O2 + SER_KO2*(Cdrug*SER_max + SER_Km)/(Cdrug + SER_Km))/(C_O2 + SER_KO2)\n\
         If a DNA-damage inhibiter, this is b_inhibit"},
        {"PARENT_CT1_12", "Sensitisation of the cells to radiation is determined by three parameters.  The usual radiation kill parameters OER_alpha and OER_beta are multiplied by the sensitisation enhancement ratio SER: \n\
         SER = (C_O2 + SER_KO2*(Cdrug*SER_max + SER_Km)/(Cdrug + SER_Km))/(C_O2 + SER_KO2)\n\
         A negative value indicates that the drug is a DNA-damage inhibiter, making the two SER parameters into inhibition model parameters"},
        {"PARENT_CT1_13",  "n_O2 is one of the three parameters of the basic sigmoid function of C_O2 that determines metabolism rate: \n\
 metabolism rate = dMdt = Cdrug.(1 - C2 + C2.KO2^n_O2/(KO2^n_O2 + C_O2^n_O2)).Kmet0   If Vmax > 0, Kmet0 is replaced by (Kmet0 + Vmax*Cdrug/(Km + Cdrug)) \n\
 This the rate of transformation of parent drug to metabolite 1, or of metabolite 1 to metabolite 2, or removal of metabolite 2"},
        {"PARENT_CT1_14", "The death probability of a drug-tagged cell at time of division."},
        {"PARENT_CT1_15", "The kill probability rate parameter.\n\
 E.g. for kill model 1, kill probability rate r = Kd*dM/dt = Kd*kmet*Cdrug, and for duration t, SF = exp(-rt)"},
        {"PARENT_CT1_16", "This box is ticked if the drug is cytotoxic and kill parameters are provided."},
        {"PARENT_CT1_17", "Using dMdt = Cdrug*(1 - C2 + C2*KO2/(KO2 + C_O2))*Kmet0, the kill probability Pkill in time dt for each model is: \n\
1. Kd*dMdt*dt  2. Kd*Cdrug*dMdt*dt  3. Kd*dMdt^2*dt  4. Kd*Cdrug*dt  5. Kd*Cdrug^2*dt"},
        {"PARENT_CT1_18", "This box is ticked if the drug sensitises the cells to radiation."},
    };

    PARAM_SET params[] = {

{"GUI_VERSION_NAME", 0, 0, 0,
 "GUI0.00",
 "GUI program version number."},

{"DLL_VERSION_NAME", 0, 0, 0,
 "DLL0.00",
 "DLL version number."},

{"INITIAL_COUNT", 10000, 0, 0,
"Initial number of tumour cells",
"Initial number of tumour cells"},

{"USE_LOGNORMAL_DIST", 1, 0, 1,
"Use lognormal distribution",
"The divide time will be a random variate from a log-normal distribution. \n\
 Otherwise checkpoint times are exponentially distributed and base phase times are fixed"},

{"DIVIDE_TIME_1_MEDIAN", 18.62, 0, 0,
"Median (h)",
"The time taken for tumour cell division has a lognormal distribution, described by the median and shape parameters. \n\
[hours]"},

{"DIVIDE_TIME_1_SHAPE", 1.1, 0, 0,
"Shape parameter",
"The time taken for tumour cell division has a lognormal distribution, described by the median and shape parameters."},

{"DIVIDE_TIME_2_MEDIAN", 19, 0, 0,
"Division time median parameter",
"The time taken for tumour cell division has a lognormal distribution, described by the median and shape parameters. \n\
[hours]"},

{"DIVIDE_TIME_2_SHAPE", 1.1, 0, 0,
"Division time shape parameter",
"The time taken for tumour cell division has a lognormal distribution, described by the median and shape parameters."},

{"V_DEPENDENT_GROWTH_RATE", 0, 0, 1,
"V-dependent growth rate",
"The growth rate of a cell is proportional to the volume."},

{"RANDOMISE_INITIAL_V", 1, 0, 1,
"Randomise initial cell volumes",
"The volumes of the initial cell population are randomised."},

{"NDAYS", 8.0, 0.0, 30.0,
"Number of days",
"Length of the simulation.\n\
[days]"},

{"DELTA_T", 600, 0, 0,
"Time step",
"Length of main time step, for cell death, division, etc.  Should be a divisor of 3600. \n\
[mins]"},

{"NT_CONC", 1, 0, 0,
"Number of ODE solver sub-steps.",
"The number of subdivisions of the major time step, for the ODE diffusion-reaction solver."},

{"VCELL_PL", 1.0, 0, 0,
"Cell volume",
"Nominal cell volume."},

{"WELL_AREA", 0.33, 0, 0,
"Well area",
"Cross-sectional area of the well."},

{"MEDIUM_VOLUME", 0.2, 0, 0,
"Medium volume",
"Volume of the medium in which the spheroid is growing."},

{"FULLY_MIXED", 0, 0, 1,
"Medium is fully mixed?",
"The medium is fully mixed"},

{"VDIVIDE0", 1.6, 0, 0,
"Nominal divide volume",
"Nominal multiple of normal cell volume at which division occurs."},

{"DVDIVIDE", 0.3, 0, 0,
"Divide volume variation",
"Variation (+/-) about nominal divide volume multiple."},

{"MM_THRESHOLD", 0.1, 0, 0,
"Michaelis-Menten O2 threshold",
"O2 concentration at which the 'soft-landing' adjustment to the Michaelis-Menten function kicks in.\n\
[uM]"},

{"ANOXIA_THRESHOLD", 0.15, 0, 0,
"Tag threshold",
"A cell begins to experience starvation of oxygen (anoxia) or glucose (aglucosia) leading to cell death at the oxygen/glucose concentration given by this threshold value."},

{"ANOXIA_TAG_TIME", 3.0, 0, 0,
"Tag time limit",
"Length of time under anoxia (O2 < anoxia threshold) or aglucosia (glucose < aglucosia threshold) after which a cell is tagged to die of anoxia or aglucosia."},

{"ANOXIA_DEATH_TIME", 3.0, 0, 0,
"Death delay time",
"Time taken for a cell to die after it is tagged to die of anoxia or aglucosia."},

{"AGLUCOSIA_THRESHOLD", 0.15, 0, 0,
"Aglucosia threshold",
"A cell begins to experience aglucosia leading to cell death at the glucose concentration given by this threshold value."},

{"AGLUCOSIA_TAG_TIME", 3.0, 0, 0,
"Aglucosia time limit",
"Length of time under aglucosia (glucose < aglucosia threshold) after which a cell is tagged to die of aglucosia.\n\
[h]"},

{"AGLUCOSIA_DEATH_TIME", 3.0, 0, 0,
"Aglucosia death delay time",
"Time taken for a cell to die after it is tagged to die of aglucosia.\n\
[h]"},

{"TEST_CASE", 0, 0, 0,
"Test case #",
"Number of the test case to run.  The default value of 0 is for a normal run"},

{"SEED1", 1234, 0, 0,
"First RNG seed",
"The random number generator is seeded by a pair of integers.  Changing the seed generates a different Monte Carlo realization."},

{"SEED2", 5678, 0, 0,
"Second RNG seed",
"The random number generator is seeded by a pair of integers.  Changing the seed generates a different Monte Carlo realization."},

{"NCPU", 1, 1, 8,
"Number of CPUs",
"Number of CPUs to use for the simulation."},

{"NCELLTYPES", 2, 0, 0,
"Number of cell types",
"Maximum number of cell types in the spheroid.  The initial percentage of each type must be specified"},

{"CELLPERCENT_1", 100, 0, 100,
"Percentage of cell type 1",
"Percentage of cell type 1"},

{"CELLPERCENT_2", 0, 0, 100,
"Percentage of cell type 2",
"Percentage of cell type 2"},

{"NT_ANIMATION", 1, 0, 0,
 "Animation interval (timesteps)",
 "Interval between animation screen updates (timesteps).  One timestep = 15 sec."},

{"SHOW_PROGENY", 0, 0, 0,
 "Show descendants of cell #",
 "All the descendants of cell with the specified ID are highlighted.  (0 = no selection)"},

{"USE_OXYGEN", 1, 0, 1,
"Use Oxygen?",
"Oxygen is simulated"},

//{"OXYGEN_GROWTH", 1, 0, 1,
//"Oxygen growth?",
//"The rate of growth of a cell is the maximum rate multiplied by the fractional rates of metabolism of both O2 and glucose"},

//{"OXYGEN_DEATH", 1, 0, 1,
//"Anoxia death?",
//"Oxygen controls death by anoxia"},

{"OXYGEN_DIFF_COEF", 2.0e-5, 0, 0,
 "Spheroid diffusion coeff",
 "Constituent diffusion coefficient in the spheroid"},

{"OXYGEN_MEDIUM_DIFF", 5.0e-5, 0, 0,
 "Medium diffusion coeff",
 "Constituent diffusion coefficient in the medium"},

{"OXYGEN_CELL_DIFF_IN", 600, 0, 0,
 "Cell influx parameter Kin",
 "Cell membrane diffusion constant Kin"},

{"OXYGEN_CELL_DIFF_OUT", 600, 0, 0,
 "Cell efflux parameter Kout",
 "Cell membrane diffusion constant Kout"},

{"OXYGEN_BDRY_CONC", 0.18, 0, 0,
 "Boundary concentration",
 "Constituent concentration in the medium"},

{"OXYGEN_CONSTANT", 0, 0, 1,
 "Constant concentration",
 "Extracellular concentration to be held constant everywhere at the specified boundary value"},

{"OXYGEN_CONSUMPTION", 6.25e-17, 0, 0,
 "Max consumption rate",
 "Maximum rate of consumption of the constituent"},

{"OXYGEN_MM_KM", 1.33, 0, 0,
 "Michaelis-Menten Km",
 "Michaelis-Menten Km (uM)"},

{"OXYGEN_HILL_N", 2, 1, 2,
 "Hill function N",
 "Oxygen uptake rate Hill function N"},

{"USE_GLUCOSE", 1, 0, 1,
"Use Glucose?",
"Glucose is simulated"},

    {"C_G_BASE", 1, 0, 1,
    "C_G_base",
    "If GLUCOSE_BDRY_CONC < 0, it is used as a factor to multiply C_G_BASE to get the actual GLUCOSE_BDRY_CONC (for PEST runs)"},

//{"GLUCOSE_GROWTH", 1, 0, 1,
//"Glucose growth?",
//"The rate of growth of a cell is the maximum rate multiplied by the fractional rates of metabolism of both O2 and glucose"},

//{"GLUCOSE_DEATH", 1, 0, 1,
//"Aglucosia death?",
//"Glucose controls death by aglucosia"},

{"GLUCOSE_DIFF_COEF", 3.0e-7, 0, 0,
 "Spheroid diffusion coeff",
 "GLUCOSE diffusion coefficient"},

{"GLUCOSE_MEDIUM_DIFF", 8.0e-6, 0, 0,
 "Medium diffusion coeff",
 "Constituent diffusion coefficient in the medium"},

{"GLUCOSE_CELL_DIFF_IN", 100, 0, 0,
 "Membrane diff constant",
 "Cell membrane diffusion coefficient Kin"},

{"GLUCOSE_CELL_DIFF_OUT", 100, 0, 0,
 "Membrane diff constant",
 "Cell membrane diffusion coefficient Kout"},

{"GLUCOSE_BDRY_CONC", 5.5, 0, 0,
 "Boundary concentration",
 "GLUCOSE boundary concentration"},

{"GLUCOSE_CONSTANT", 0, 0, 1,
 "Constant concentration",
 "Extracellular concentration to be held constant everywhere at the specified boundary value"},

{"GLUCOSE_CONSUMPTION", 1.3e-16, 0, 0,
 "Max consumption rate",
 "GLUCOSE consumption rate"},

{"GLUCOSE_MM_KM", 300, 0, 0,
 "Michaelis-Menten Km",
 "Michaelis-Menten Km (uM)"},

{"GLUCOSE_HILL_N", 2, 1, 2,
 "Hill function N",
 "Glucose uptake rate Hill function N"},


//==========================
// Radiotherapy parameters
//==========================

{"RADIATION_ALPHA_H_1", 0.0738, 0, 0,
"Alpha (hypoxia)",
"alpha for irradiation of cells under anoxia (zero oxygen)"},

{"RADIATION_BETA_H_1", 0.00725, 0, 0,
"Beta (hypoxia)",
"beta for irradiation of cells under anoxia (zero oxygen)"},

{"RADIATION_OER_ALPHA_1", 2.5, 0, 0,
"OER alpha",
"Maximum oxygen enhancement ratio for alpha component of radiosensitivity "},

{"RADIATION_OER_BETA_1", 2.5, 0, 0,
"OER beta",
"Maximum oxygen enhancement ratio for beta component of radiosensitivity"},

{"RADIATION_KM_1", 4.3e-3, 0, 0,
"Km for radiosensitivity",
"Oxygen concentration for half maximal radiosensitivity relative to hypoxic cell exposure"},

{"RADIATION_ALPHA_H_2", 0.0473, 0, 0,
"Alpha (hypoxia)",
"alpha for irradiation of cells under anoxia (zero oxygen)"},

{"RADIATION_BETA_H_2", 0.0017, 0, 0,
"Beta (hypoxia)",
"beta for irradiation of cells under anoxia (zero oxygen)"},

{"RADIATION_OER_ALPHA_2", 2.5, 0, 0,
"OER alpha",
"Maximum oxygen enhancement ratio for alpha component of radiosensitivity "},

{"RADIATION_OER_BETA_2", 3.0, 0, 0,
"OER beta",
"Maximum oxygen enhancement ratio for beta component of radiosensitivity"},

{"RADIATION_KM_2", 4.3e-3, 0, 0,
"Km for radiosensitivity",
"Oxygen concentration for half maximal radiosensitivity relative to hypoxic cell exposure"},

    {"USE_CELL_CYCLE", 1,0,1,
     "Use cell cycle with G1, S, G2, M phases",
     "Cell cycle parameters determine the time spent in each phase.\n\
     In the case of G1 and G2, an exponentially distributed random delay is added"},

     {"USE_SYNCHRONISE", 0,0,1,
     "Synchronise cell cycles?",
     "Synchronise initial cell phases to start of M phase"},

     {"F_G1_1", 0.36, 0, 0,
     "G1 phase base fraction",
     "Fraction of cycle spent in phase G1"},

     {"F_S_1", 0.49, 0, 0,
     "S phase base fraction",
     "Fraction of cycle spent in phase S"},

     {"F_G2_1", 0.13, 0, 0,
     "G2 phase base fraction",
     "Fraction of cycle spent in phase G2"},

     {"F_M_1", 0.02, 0, 0,
     "M phase base fraction",
     "Fraction of cycle spent in phase M"},

     {"APOPTOSIS_MEDIAN_1", 4, 0, 0,
     "Apoptosis median (hr)",
     "Apoptosis delay = 18 + log-normal variate defined by median, shape"},

     {"APOPTOSIS_SHAPE_1", 1.5, 0, 0,
     "Apoptosis shape",
     "Shape parameter for log-normal apoptosis delay variate"},

      {"F_G1_2", 0.36, 0, 0,
      "G1 phase base fraction)",
      "Fraction of cycle spent in phase G1"},

      {"F_S_2", 0.49, 0, 0,
      "S phase base fraction",
      "Fraction of cycle spent in phase S"},

      {"F_G2_2",0.13, 0, 0,
      "G2 phase base fraction",
      "Fraction of cycle spent in phase G2"},

      {"F_M_2", 0.02, 0, 0,
      "M phase base fraction",
      "Fraction of cycle spent in phase M"},


      {"APOPTOSIS_MEDIAN_2", 4, 0, 0,
      "Apoptosis median (hr)",
      "Apoptosis delay = 18 + log-normal variate defined by median, shape"},

      {"APOPTOSIS_SHAPE_2", 1.5, 0, 0,
      "Apoptosis shape",
      "Shape parameter for log-normal apoptosis delay variate"},

     {"PHASE_HOURS", 1, 0, 0,
     "PEST run case",
     "PEST run case"},

     {"BASERATE", 0.0007, 0, 0,
     "Base apoptosis rate parameter",
     "Base apoptosis rate parameter"},

     {"MITRATE", 0.000067, 0, 0,
     "Mitosis rate parameter",
     "Mitosis rate parameter"},

     {"MSURVIVAL", 0.0016, 0, 0,
     "Mitosis survival probability",
     "Mitosis survival probability"},

     {"KLETHAL", 0.0944, 0, 0,
     "MisrepRate scaling factor",
     "MisrepRate scaling factor"},

//     {"NCPPARAMS", 3, 0, 0,
//     "Fixed",
//     "Fixed"},

     {"KATM1G1", 0, 0, 0,
     "KATM1G1 parameter",
     "G1 pATM production rate parameter: rate = KATM1G1*ATM_DSB"},

     {"KATM1S", 0.01, 0, 0,
     "KATM1S parameter",
     "S pATM production rate parameter: rate = KATM1S*ATM_DSB"},

     {"KATM1G2", 0.01, 0, 0,
     "KATM1G2 parameter",
     "G2 pATM production rate parameter: rate = KATM1G2*ATM_DSB"},

     {"KATM2G1", 0, 0, 0,
     "KATM2G1 parameter",
     "G1 pATM decay rate constant"},

     {"KATM2S", 7.1, 0, 0,
     "KATM2S parameter",
     "S pATM decay rate constant"},

     {"KATM2G2", 0.3465, 0, 0,
     "KATM2G2 parameter",
     "G2 pATM decay rate constant"},

     {"KATM3G1", 0, 0, 0,
     "KATM3G1 parameter",
     "With k1=KATM3, k2=KATM4, x = pATM, ATM slowdown factor = 1 - k1*x/(k2+x)"},

     {"KATM3S", 2.86, 0, 0,
     "KATM3S parameter",
     "With k1=KATM3, k2=KATM4, x = pATM, ATM slowdown factor = 1 - k1*x/(k2+x)"},

     {"KATM3G2", 3.14, 0, 0,
     "KATM3G2 parameter",
     "With k1=KATM3, k2=KATM4, x = pATM, ATM slowdown factor = 1 - k1*x/(k2+x)"},

     {"KATM4G1", 0, 0, 0,
     "KATM4G1 parameter",
     "With k1=KATM3, k2=KATM4, x = pATM, ATM slowdown factor = 1 - k1*x/(k2+x)"},

     {"KATM4S", 3.0, 0, 0,
     "KATM4S parameter",
     "With k1=KATM3, k2=KATM4, x = pATM, ATM slowdown factor = 1 - k1*x/(k2+x)"},

     {"KATM4G2", 59.1, 0, 0,
     "KATM4G2 parameter",
     "With k1=KATM3, k2=KATM4, x = pATM, ATM slowdown factor = 1 - k1*x/(k2+x)"},

     {"KATR1S", 0, 0, 0,
     "KATR1S parameter",
     "pATR production rate parameter: rate = KATR1*ATR_DSB"},

     {"KATR1G2", 0, 0, 0,
     "KATR1G2 parameter",
     "pATR production rate parameter: rate = KATR1*ATR_DSB"},

     {"KATR2S", 0.167, 0, 0,
     "KATR2S parameter",
     "pATR decay rate constant"},

     {"KATR2G2", 1.53, 0, 0,
     "KATR2G2 parameter",
     "pATR decay rate constant"},

     {"KATR3S", 0.0001, 0, 0,
     "KATR3S parameter",
     "With k1=KATR3, k2=KATR4, x = pATR, ATR slowdown factor = 1 - k1*x/(k2+x)"},

     {"KATR3G2", 33.6, 0, 0,
     "KATR3G2 parameter",
     "With k1=KATR3, k2=KATR4, x = pATR, ATR slowdown factor = 1 - k1*x/(k2+x)"},

     {"KATR4S", 0.0002, 0, 0,
     "KATR4S parameter",
     "With k1=KATR3, k2=KATR4, x = pATR, ATR slowdown factor = 1 - k1*x/(k2+x)"},

     {"KATR4G2", 0.0296, 0, 0,
     "KATR4G2 parameter",
     "With k1=KATR3, k2=KATR4, x = pATR, ATR slowdown factor = 1 - k1*x/(k2+x)"},

     {"PCOMPLEX", 0.85, 0, 0,
     "PCOMPLEX parameter",
     "PCOMPLEX parameter"},

     {"PHRSIMPLE", 0.033, 0, 0,
     "PHRSIMPLE parameter",
     "PHRSIMPLE parameter"},

     {"CHALF", 0.2, 0, 0,
     "CHALF parameter",
     "CHALF parameter"},

     {"PREASS", 0.066, 0, 0,
     "PREASS parameter",
     "PREASS parameter"},

     {"MDRREP", 0.55, 0, 0,
     "MDRREP parameter",
     "MDRREP parameter"},

     {"MDRFID", 0.5, 0, 0,
     "MDRFID parameter",
     "MDRFID parameter"},


     {"KCC2A", 4.0, 0, 0,
     "KCC2A parameter",
     "KCC2A parameter"},

     {"KCC2E", 0, 0, 0,
     "KCC2E parameter",
     "KCC2E parameter"},

     {"KD2E", 0, 0, 0,
     "KD2E parameter",
     "KD2E parameter"},

     {"KD2T", 0.785, 0, 0,
     "KD2T parameter",
     "KD2T parameter"},

     {"KE2CC", 0.0005, 0, 0,
     "KE2CC parameter",
     "KE2CC parameter"},

     {"KM1", 108, 0, 0,
     "KM1 parameter",
     "KM1 parameter"},

     {"KM10", 7.83, 0, 0,
     "KM10 parameter",
     "KM10 parameter"},

     {"KT2CC", 4.47, 0, 0,
     "KT2CC parameter",
     "KT2CC parameter"},

     {"KTI2T", 9.61, 0, 0,
     "KTI2T parameter",
     "KTI2T parameter"},

     {"KM10T", 0.002, 0, 0,
     "KM10T parameter",
     "KM10T parameter"},

     {"CC_TOT", 5, 0, 0,
     "CC_TOT parameter",
     "CC_TOT parameter"},

     {"ATR_TOT", 10, 0, 0,
     "ATR_TOT parameter",
     "ATR_TOT parameter"},

     {"ATM_TOT", 10, 0, 0,
     "ATM_TOT parameter",
     "ATM_TOT parameter"},

     {"CC_ACT0", 0, 0, 0,
     "CC_ACT0 parameter",
     "CC_ACT0 parameter"},

     {"CC_THRESHOLD", 0.9, 0, 0,
     "CC_THRESHOLD parameter",
     "CC_THRESHOLD parameter"},

     {"NORM_FACTOR", 0.01, 0, 0,
     "NORM_FACTOR parameter",
     "NORM_FACTOR parameter"},


{"HYPOXIA_1", 0.1, 0, 0,
"Hypoxia threshold 1",
"Hypoxia threshold 1"},

{"HYPOXIA_2", 1.0, 0, 0,
"Hypoxia threshold 2",
"Hypoxia threshold 2"},

{"HYPOXIA_3", 4.0, 0, 0,
"Hypoxia threshold 3",
"Hypoxia threshold 3"},

{"HYPOXIA_THRESHOLD", 4.0, 0, 0,
"Hypoxia threshold",
"Hypoxia threshold"},

{"GROWTH_FRACTION_1", 0.25, 0, 0,
"Growth fraction threshold 1",
"Growth fraction threshold 1"},

{"GROWTH_FRACTION_2", 0.1, 0, 0,
"Growth fraction threshold 2",
"Growth fraction threshold 2"},

{"GROWTH_FRACTION_3", 0.01, 0, 0,
"Growth fraction threshold 3",
"Growth fraction threshold 3"},

{"DRUG_THRESHOLD", 1.0e-6, 0, 0,
 "Drug Threshold",
 "Threshold drug concentration - when all intracellular and extracellular concentrations fall below this level, the drug concentrations everywhere are set to zero"},

{"DRUG_LABEL_THRESHOLD", 0, 0, 0,
"Label Threshold",
"Threshold label-drug concentration - when a labelling drug (e.g. EDU) is used, this is the threshold for a cell to be considered as labelled"},

{"SPCRAD", 200.0, 0, 0,
"Spectral radius",
"Spectral radius value used by RKC solver"},

// {"SAVE_FACS_DATA",0,0,1,
//  "Save FACS data",
//  "Save data for FACS at a specified interval"},

// {"SAVE_FACS_DATA_FILE_NAME",0,0,0,
//  "facs_data",
//  "Base file name for saving FACS data"},

// {"SAVE_FACS_DATA_TSTART",0,0,0,
//  "Start time",
//  "Start time for saving FACS data"},

// {"SAVE_FACS_DATA_INTERVAL",0,0,0,
//  "Interval",
//  "Time interval for saving FACS data"},

// {"SAVE_FACS_DATA_NUMBER",1,0,0,
//  "Number",
//  "Number of times to save FACS data"},


// This is the end of the parameters that are actually read by the DLL
// Entries after this point are QMyLabel dummies, to enable display of explanatory info  - no input data is transmitted,
// followed by the list of time-series and profile plots selected for this run.

{"DUMMY_HYPOXIA_THRESHOLD", 0, 0, 0,
"Hypoxia threshold",
"Select the intracellular O2 level below which the cell is counted as hypoxic"},

{"DUMMY_GROWTH_FRACTION", 0, 0, 0,
"Growth fraction",
"Select the threshold fraction of average growth rate (i.e. with no nutrient limits) used to count slow-growing cells"},

// Time-series plots
    {"nlive",                     1, 0,1,"","Number of live cells"},
    {"nviable",                   1, 0,1,"","Number of viable cells"},
    {"nonviable",                 1, 0,1,"","Total number of non-viable cells"},
    {"ndrugAdead",                0, 0,1,"","Total number of cells that have been killed by drugA"},
//    {"ndrugBdead",                0, 0,1,"","Total number of cells that have been killed by drugB"},
    {"nradiationdead",            0, 0,1,"","Total number of cells that have been killed by radiation"},
    {"ndead",                     1, 0,1,"","Total number of cellls that have died"},
//    {"nATPtagged",                0, 0,1,"","Current number of cells tagged to die by low ATP"},
//    {"nGLNtagged",                1, 0,1,"","Current number of cells tagged to die by low GLN"},
    {"ndrugAtagged",              0, 0,1,"","Current number of cells tagged to die by drugA"},
//    {"ndrugBtagged",              0, 0,1,"","Current number of cells tagged to die by drugB"},
    {"nradiationtagged",          0, 0,1,"","Current number of cells tagged to die by radiation"},
    {"viablefraction",            0, 0,1,"","Fraction of cells that are viable"},
    {"hypoxicfraction",           0, 0,1,"","Fraction of cells with oxygen level below the specified threshold for hypoxia"},
    {"clonohypoxicfraction",      0, 0,1,"","Fraction of clonogenic cells with oxygen level below the specified threshold for hypoxia"},
    {"growthfraction",            0, 0,1,"","Percentage of cells that are growing at a rate less than the specified fraction of the mean growth rate with no nutrient limits"},
    {"nogrowfraction",            0, 0,1,"","Percentage of cells that are not growing (insufficient ATP rate for growth)"},
    {"clonofraction",             0, 0,1,"","Percentage of cells that are clonogenic (will give a colony >= 50)"},
    {"platingefficiency",         0, 0,1,"","Percentage of live cells that are viable"},
    {"ECoxygen",                  1, 0,1,"","EC concentration of oxygen in the medium (bottom)"},
    {"ECglucose",                 1, 0,1,"","EC concentration of glucose in the medium (bottom)"},
//    {"EClactate",                 1, 0,1,"","EC concentration of lactate in the medium (bottom)"},
//    {"ECglutamine",               1, 0,1,"","EC concentration of glutamine in the medium (bottom)"},
//    {"ECother",                   0, 0,1,"","EC concentration of other nutrient in the medium (bottom)"},
    {"ECdrugA",                   1, 0,1,"","EC concentration of drug A in the medium (bottom)"},
    {"ECdrugAmet1",               1, 0,1,"","EC concentration of drug A metabolite 1 in the medium (bottom)"},
//    {"ECdrugAmet2",               0, 0,1,"","EC concentration of drug A metabolite 2 in the medium (bottom)"},
//    {"ECdrugB",                   0, 0,1,"","EC concentration of drug B in the medium (bottom)"},
//    {"ECdrugBmet1",               0, 0,1,"","EC concentration of drug B metabolite 1 in the medium (bottom)"},
//    {"ECdrugBmet2",               0, 0,1,"","EC concentration of drug B metabolite 2 in the medium (bottom)"},
    {"ICoxygen",                  1, 0,1,"","IC concentration of oxygen"},
    {"ICglucose",                 1, 0,1,"","IC concentration of glucose"},
//    {"IClactate",                 1, 0,1,"","IC concentration of lactate "},
//    {"ICglutamine",               1, 0,1,"","IC concentration of glutamine"},
//    {"ICother",                   0, 0,1,"","IC concentration of other nutrient"},
//    {"ICpyruvate",                1, 0,1,"","IC concentration of pyruvate"},
//    {"ICATP",                     1, 0,1,"","IC concentration of ATP"},
    {"ICdrugA",                   1, 0,1,"","IC concentration of drug A"},
    {"ICdrugAmet1",               1, 0,1,"","IC concentration of drug A metabolite 1"},
//    {"ICdrugAmet2",               0, 0,1,"","IC concentration of drug A metabolite 2"},
//    {"ICdrugB",                   0, 0,1,"","IC concentration of drug B"},
//    {"ICdrugBmet1",               0, 0,1,"","IC concentration of drug B metabolite 1"},
//    {"ICdrugBmet2",               0, 0,1,"","IC concentration of drug B metabolite 2"},
    {"Medoxygen",                 1, 0,1,"","Average medium concentration of oxygen"},
    {"Medglucose",                1, 0,1,"","Average medium concentration of glucose"},
//    {"Medlactate",                0, 0,1,"","Average medium concentration of lactate"},
//    {"Medglutamine",              0, 0,1,"","Average medium concentration of glutamine"},
//    {"Medother",                  0, 0,1,"","Average medium concentration of other nutrient"},
    {"MeddrugA",                  1, 0,1,"","Average medium concentration of drug A"},
    {"MeddrugAmet1",              1, 0,1,"","Average medium concentration of drug A metabolite 1"},
//    {"MeddrugAmet2",              0, 0,1,"","Average medium concentration of drug A metabolite 2"},
//    {"MeddrugB",                  0, 0,1,"","Average medium concentration of drug B"},
//    {"MeddrugBmet1",              0, 0,1,"","Average medium concentration of drug B metabolite 1"},
//    {"MeddrugBmet2",              0, 0,1,"","Average medium concentration of drug B metabolite 2"},
    {"doublingtime",              0, 0,1,"","Average doubling time"},
    {"Orate",                     0, 0,1,"","Normalised oxygen consumption rate"},
    {"Grate",                     0, 0,1,"","Normalised glycolysis rate"},
//    {"Prate",                     0, 0,1,"","Normalised pyruvate utilisation rate"},
//    {"Glnrate",                   1, 0,1,"","Normalised glutamine utilisation rate"},
//    {"ONrate",                    1, 0,1,"","Normalised other nutrient utilisation rate"},
//    {"Arate",                     1, 0,1,"","Normalised ATP production rate"},
//    {"Irate",                     1, 0,1,"","Normalised rate of production of anabolic intermediates"},
//    {"f_G",                       0, 0,1,"","f_G"},
//    {"f_P",                       0, 0,1,"","f_P"},
//    {"HIF-1",                     1, 0,1,"","HIF-1"},
//    {"PDK1",                      1, 0,1,"","PDK1"},
    {"dividerate",                0, 0,1,"","# divided/hour"},
    {"G1_phase",                  0, 0,1,"","G1_phase"},
    {"G1_cp",                     0, 0,1,"","G1_cp"},
    {"S_phase",                   0, 0,1,"","S_phase"},
    {"S_cp",                      0, 0,1,"","S_cp"},
    {"G2_phase",                  0, 0,1,"","G2_phase"},
    {"G2_cp",                     0, 0,1,"","G2_cp"},
    {"M_phase",                   0, 0,1,"","M_phase"},
//    {"S_phase_nonarrest",         0, 0,1,"","S_phase_nonarrest"},
//    {"nmutations",                0, 0,1,"","nmutations"},


// Profile plots
    {"MULTI",                     0, 0,1,"","Selected constituent on a line through the blob centre"},
//    {"CFSE",                      0, 0,1,"","Extracellular CFSE concentration on a line through the blob centre"},
    {"Oxygen",                    1, 0,1,"","Extracellular oxygen concentration on a line through the blob centre"},
    {"Glucose",                   1, 0,1,"","Extracellular glucose concentration on a line through the blob centre"},
//    {"Glutamine",                 0, 0,1,"","Extracellular glutamine concentration on a line through the blob centre"},
    {"Drug_A",                    1, 0,1,"","Extracellular drug A concentration on a line through the blob centre"},
    {"Drug_A_metab1",             1, 0,1,"","Extracellular drug A metabolite 1 concentration on a line through the blob centre"},
//    {"Drug_A_metab2",             0, 0,1,"","Extracellular drug A metabolite 2 concentration on a line through the blob centre"},
//    {"Drug_B",                    0, 0,1,"","Extracellular drug Bconcentration on a line through the blob centre"},
//    {"Drug_B_metab1",             0, 0,1,"","Extracellular drug B metabolite 1 concentration on a line through the blob centre"},
//    {"Drug_B_metab2",             0, 0,1,"","Extracellular drug B metabolite 2 concentration on a line through the blob centre"},

/*
    {"IC_MULTI",                  1, 0,1,"","Selected constituent on a line through the blob centre"},
    {"IC_Oxygen",                 0, 0,1,"","Intracellular oxygen concentration on a line through the blob centre"},
    {"IC_Glucose",                0, 0,1,"","Intracellular glucose concentration on a line through the blob centre"},
    {"IC_Glutamine",              0, 0,1,"","Intracellular glutamine concentration on a line through the blob centre"},
    {"IC_Drug_A",                 0, 0,1,"","Intracellular drug A concentration on a line through the blob centre"},
    {"IC_Drug_A_metab1",          0, 0,1,"","Intracellular drug A metabolite 1 concentration on a line through the blob centre"},
    {"IC_Drug_A_metab2",          0, 0,1,"","Intracellular drug A metabolite 2 concentration on a line through the blob centre"},
    {"IC_Drug_B",                 0, 0,1,"","Intracellular drug Bconcentration on a line through the blob centre"},
    {"IC_Drug_B_metab1",          0, 0,1,"","Intracellular drug B metabolite 1 concentration on a line through the blob centre"},
    {"IC_Drug_B_metab2",          0, 0,1,"","Intracellular drug B metabolite 2 concentration on a line through the blob centre"},
    {"IC_CFSE",                   0, 0,1,"","CFSE concentration on a line through the blob centre"},
    {"IC_growthrate",             0, 0,1,"","Cell growth rate on a line through the blob centre"},
    {"IC_cellvolume",             0, 0,1,"","Cell volume fraction on a line through the blob centre"},
    {"IC_O2byvolume",             0, 0,1,"","Cell volume fraction on a line through the blob centre"},
// Distribution plots
//    {"Oxygen",                    0, 0,1,"","Probability distribution of extracellular oxygen concentration"},
//    {"cellvolume",                0, 0,1,"","Probability distribution of cell volume fraction"}
*/

};
	nParams = sizeof(params)/sizeof(PARAM_SET);
	workingParameterList = new PARAM_SET[nParams];
	for (int i=0; i<nParams; i++) {
		workingParameterList[i] = params[i];
	}

    nInfolabel = sizeof(label_info)/sizeof(INFOSTRUCT);
    workingInfolabelList = new INFOSTRUCT[nInfolabel];
    for (int i=0; i<nInfolabel; i++) {
        workingInfolabelList[i] = label_info[i];
    }
    /*
    nInfocheckbox = sizeof(checkbox_info)/sizeof(INFOSTRUCT);
    workingInfocheckboxList = new INFOSTRUCT[nInfocheckbox];
    for (int i=0; i<nInfocheckbox; i++) {
        workingInfocheckboxList[i] = checkbox_info[i];
    }
    */
}


PARAM_SET Params::get_param(int k)
{
	return workingParameterList[k];
}

void Params::set_value(int k, double v)
{
	workingParameterList[k].value = v;
}

void Params::set_label(int k, QString str)
{
	workingParameterList[k].label = str;
}

void Params::get_labeltag(int i, QString *tag)
{
    *tag = workingInfolabelList[i].tag;
}

void Params::infoLabelInfo(QString tag, QString *info)
{
    for (int i=0; i<nInfolabel; i++) {
        if (tag == workingInfolabelList[i].tag) {
            *info = workingInfolabelList[i].info;
            return;
        } else {
            *info = "";
        }
    }
}

