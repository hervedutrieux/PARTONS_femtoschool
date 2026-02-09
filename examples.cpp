#include <ElementaryUtils/logger/LoggerManager.h>
#include <ElementaryUtils/parameters/Parameter.h>
#include <ElementaryUtils/parameters/Parameters.h>
#include <NumA/integration/one_dimension/IntegratorType1D.h>
#include <NumA/integration/one_dimension/QuadratureIntegrator1D.h>
#include <partons/beans/convol_coeff_function/ConvolCoeffFunctionResult.h>
#include <partons/beans/convol_coeff_function/DVCS/DVCSConvolCoeffFunctionKinematic.h>
#include <partons/beans/convol_coeff_function/DVMP/DVMPConvolCoeffFunctionKinematic.h>
#include <partons/beans/gpd/GPDKinematic.h>
#include <partons/beans/gpd/GPDType.h>
#include <partons/beans/KinematicUtils.h>
#include <partons/beans/List.h>
#include <partons/beans/PolarizationType.h>
#include <partons/beans/MesonType.h>
#include <partons/beans/observable/DVCS/DVCSObservableKinematic.h>
#include <partons/beans/observable/DVMP/DVMPObservableKinematic.h>
#include <partons/beans/observable/ObservableResult.h>
#include <partons/beans/PerturbativeQCDOrderType.h>
#include <partons/modules/active_flavors_thresholds/ActiveFlavorsThresholdsConstant.h>
#include <partons/modules/active_flavors_thresholds/ActiveFlavorsThresholdsQuarkMasses.h>
#include <partons/modules/convol_coeff_function/ConvolCoeffFunctionModule.h>
#include <partons/modules/convol_coeff_function/DVCS/DVCSCFFStandard.h>
#include <partons/modules/convol_coeff_function/DVCS/DVCSCFFNN.h>
#include <partons/modules/convol_coeff_function/DVMP/DVMPCFFGK06.h>
#include <partons/modules/evolution/gpd/GPDEvolutionApfel.h>
#include <partons/modules/gpd/GPDGK16Numerical.h>
#include <partons/modules/gpd/GPDGK19.h>
#include <partons/modules/observable/DVCS/asymmetry/DVCSAllMinus.h>
#include <partons/modules/observable/DVMP/cross_section/DVMPCrossSectionUUUMinus.h>
#include <partons/modules/process/DVCS/DVCSProcessGV08.h>
#include <partons/modules/process/DVMP/DVMPProcessGK06.h>
#include <partons/modules/running_alpha_strong/RunningAlphaStrongApfel.h>
#include <partons/modules/scales/DVCS/DVCSScalesQ2Multiplier.h>
#include <partons/modules/scales/DVMP/DVMPScalesQ2Multiplier.h>
#include <partons/modules/xi_converter/DVCS/DVCSXiConverterXBToXi.h>
#include <partons/modules/xi_converter/DVMP/DVMPXiConverterXBToXi.h>
#include <partons/ModuleObjectFactory.h>
#include <partons/Partons.h>
#include <partons/services/ConvolCoeffFunctionService.h>
#include <partons/services/DVCSConvolCoeffFunctionService.h>
#include <partons/services/DVCSObservableService.h>
#include <partons/services/DVMPConvolCoeffFunctionService.h>
#include <partons/services/DVMPObservableService.h>
#include <partons/services/GPDService.h>
#include <partons/services/ObservableService.h>
#include <partons/ServiceObjectRegistry.h>
#include <partons/utils/type/PhysicalType.h>
#include <partons/utils/type/PhysicalUnit.h>
#include <partons/modules/collinear_distribution/CollinearDistributionLHAPDF.h>
#include <partons/services/CollinearDistributionService.h>

// New examples developed for the First International School of Hadron Femtography
// Not official or endorsed by the PARTONS collaboration
// Herve Dutrieux - September 19, 2024

void extract_xdep_GPD(double xi, double t, double mu2, double min_x = 1e-4, int nbre_x = 100) {
    // Outputs the tabulated values of the GPD on a logarithmic grid in x at a given (xi, t, mu2)

    // Retrieve GPD service
    PARTONS::GPDService* pGPDService =
            PARTONS::Partons::getInstance()->getServiceObjectRegistry()->getGPDService();

    // Create GPD module with the BaseModuleFactory
    PARTONS::GPDModule* pGPDModel =
            PARTONS::Partons::getInstance()->getModuleObjectFactory()->newGPDModule(
                    PARTONS::GPDGK16::classId);

    PARTONS::List<PARTONS::GPDKinematic> gpdKinematicList;
    for (int  i = 0; i < nbre_x; i++){
            // Logarithmic grid in x with nbre_x points from min_x to 1
            gpdKinematicList.add( PARTONS::GPDKinematic(exp(-log(min_x) / (nbre_x - 1) * i + log(min_x)), xi, t, mu2, mu2) );
    }

    // Run computation
    PARTONS::Partons::getInstance()->getLoggerManager()->info("main", __func__, "Begin run");
    PARTONS::List<PARTONS::GPDResult> gpdResultList =
            pGPDService->computeManyKinematic(gpdKinematicList, pGPDModel);
    PARTONS::Partons::getInstance()->getLoggerManager()->info("main", __func__, "End run");

    // We form the output as a list of values of x and H^u directly readable by Python
    std::string Hu = "[", grid_x = "[";
    for (int  i = 0; i < nbre_x; i++){
            Hu += gpdResultList[i].getPartonDistribution(PARTONS::GPDType::H).
                        getQuarkDistribution(PARTONS::QuarkFlavor::UP).toString() + ", ";
            grid_x += gpdResultList[i].getKinematic().getX().toString() + ", ";
    }
    Hu += "]";
    grid_x += "]";

    /*     
    // Load list of kinematics from file
    PARTONS::List<PARTONS::GPDKinematic> gpdKinematicList;
    int nbre_x(2), nbre_t(2), nbre_Q2(2);
    for(int k=0; k < nbre_t; k++){
        double t( 8 * k*k / (nbre_t-1) / (nbre_t-1) );
    for(int j=0; j < nbre_Q2; j++){
        double q2( pow(2, 1+0.416666667*j) );
    for(int i=0; i < nbre_x; i++){
        double xi (exp(-log(0.0001) / (nbre_x-1) * i + log(0.0001)) );
        gpdKinematicList.add( PARTONS::GPDKinematic(xi, xi, t, q2, q2) );
    }}}
    

    // Run computation
    PARTONS::List<PARTONS::GPDResult> gpdResultList =
            pGPDService->computeManyKinematic(gpdKinematicList, pGPDModel);

    // Print results
    std::string Hu = "[", kins = "[";
    for(int i=0; i < nbre_x*nbre_t*nbre_Q2; i++){
        Hu += std::to_string( gpdResultList[i].getPartonDistribution(PARTONS::GPDType::H).getQuarkDistribution(PARTONS::QuarkFlavor::UP).getQuarkDistributionPlus() )+",";
        kins += "["+std::to_string(gpdKinematicList[i].getX().getValue())+","+std::tostring(gpdKinematicList[i].getT().getValue())+","+std::tostring(gpdKinematicList[i].getMuF2().getValue())+"],";
    }

    PARTONS::Partons::getInstance()->getLoggerManager()->info("main", __func__, kins);
    PARTONS::Partons::getInstance()->getLoggerManager()->info("main", __func__, Hu);
    */

    PARTONS::Partons::getInstance()->getLoggerManager()->info("main", __func__, grid_x);
    PARTONS::Partons::getInstance()->getLoggerManager()->info("main", __func__, Hu);

    // Remove pointer references
    // Module pointers are managed by PARTONS
    PARTONS::Partons::getInstance()->getModuleObjectFactory()->updateModulePointerReference(
            pGPDModel, 0);
    pGPDModel = 0;
}






std::map<int, double> GPDModel(double const& x, double const&)
{
  // This technical functions maps a GPD Model from the PARTONS format to APFEL for use in high performance evolution

  // Retrieve GPD service
  PARTONS::GPDService* pGPDService =
          PARTONS::Partons::getInstance()->getServiceObjectRegistry()->getGPDService();

  // Create GPD module with the BaseModuleFactory
  PARTONS::GPDModule* pGPDModel =
          PARTONS::Partons::getInstance()->getModuleObjectFactory()->newGPDModule(
                  PARTONS::GPDGK16::classId);

  // Create a GPDKinematic(x, xi, t, MuF2, MuR2) to compute
  const double xi = 1e-2;
  PARTONS::GPDKinematic gpdKinematic(x, xi, -0.1, 4., 4.);
  // Run computation
  PARTONS::GPDResult gpdResult = pGPDService->computeSingleKinematic(
          gpdKinematic, pGPDModel);

  const double uplus = x * gpdResult.getPartonDistribution(PARTONS::GPDType::H).
              getQuarkDistribution(PARTONS::QuarkFlavor::UP).getQuarkDistributionPlus();
  const double uminus = x * gpdResult.getPartonDistribution(PARTONS::GPDType::H).
              getQuarkDistribution(PARTONS::QuarkFlavor::UP).getQuarkDistributionMinus();
  const double dplus = x * gpdResult.getPartonDistribution(PARTONS::GPDType::H).
              getQuarkDistribution(PARTONS::QuarkFlavor::DOWN).getQuarkDistributionPlus();
  const double dminus = x * gpdResult.getPartonDistribution(PARTONS::GPDType::H).
              getQuarkDistribution(PARTONS::QuarkFlavor::DOWN).getQuarkDistributionMinus();
  const double splus = x * gpdResult.getPartonDistribution(PARTONS::GPDType::H).
              getQuarkDistribution(PARTONS::QuarkFlavor::STRANGE).getQuarkDistributionPlus();
  const double sminus = x * gpdResult.getPartonDistribution(PARTONS::GPDType::H).
              getQuarkDistribution(PARTONS::QuarkFlavor::STRANGE).getQuarkDistributionMinus();

  const double upv  = uminus;
  const double dnv  = dminus;
  const double glu  = gpdResult.getPartonDistribution(PARTONS::GPDType::H).
              getGluonDistribution().getGluonDistribution();
  const double dbar = 0.5 * (dplus - dminus);
  const double ubar = 0.5 * (uplus - uminus);
  const double sbar = 0.5 * (splus - sminus);

  // Construct QCD evolution basis conbinations.
  double const Gluon   = glu;
  double const Singlet = dnv + 2 * dbar + upv + 2 * ubar + 2 * sbar;
  double const T3      = upv + 2 * ubar - dnv - 2 * dbar;
  double const T8      = upv + 2 * ubar + dnv + 2 * dbar - 4 * sbar;
  double const Valence = upv + dnv;
  double const V3      = upv - dnv;

  // Fill in map in the QCD evolution basis.
  std::map<int, double> QCDEvMap;
  QCDEvMap[0]  = Gluon;
  QCDEvMap[1]  = Singlet;
  QCDEvMap[2]  = Valence;
  QCDEvMap[3]  = T3;
  QCDEvMap[4]  = V3;
  QCDEvMap[5]  = T8;
  QCDEvMap[6]  = Valence;
  QCDEvMap[7]  = Singlet;
  QCDEvMap[8]  = Valence;
  QCDEvMap[9]  = Singlet;
  QCDEvMap[10] = Valence;
  QCDEvMap[11] = Singlet;
  QCDEvMap[12] = Valence;

  return QCDEvMap;
}






void using_GPDEvolution_high_perf(double mu2F, double mu2I, double min_x = 1e-4, int nbre_x = 100) {
    // Outputs the x dependence of a GPD model at a higher scale on a logarithmic grid in x
    const double xi = 1e-2;

    // Choice of grid
    const apfel::Grid g{{apfel::SubGrid{100, 1e-7, 3}, apfel::SubGrid{50, 1e-1, 5}, apfel::SubGrid{40, 8e-1, 5}}};

    // Running of the coupling
    const std::vector<double> Thresholds = {0, 0, 0, 1.4, 4.75, 175};
    const int PerturbativeOrder = 0;
    const double AlphaQCDRef = 0.118;
    const double MuAlphaQCDRef = 91.1876;
    apfel::AlphaQCD a{AlphaQCDRef, MuAlphaQCDRef, Thresholds, PerturbativeOrder};
    const apfel::TabulateObject<double> Alphas{a, 100, 0.9, 1001, 3};
    const auto as = [&] (double const& mu) -> double{ return Alphas.Evaluate(mu); };

    const auto GpdObj   = InitializeGpdObjects(g, Thresholds, xi);
    const auto EvolvedGPDs = BuildDglap(GpdObj, GPDModel, sqrt(mu2I), PerturbativeOrder, as);
    const apfel::TabulateObject<apfel::Set<apfel::Distribution>> TabulatedGPDs{*EvolvedGPDs, 30, 1, 100, 3};

    // Evaluates the GPD at mu2F (can evaluate at many values of mu simultaneously)
    const std::map<int, apfel::Distribution> gpds = apfel::QCDEvToPhys(TabulatedGPDs.Evaluate(sqrt(mu2F)).GetObjects());
    // Constructs the result of the GPD
    std::string Hu = "[", grid_x = "[";
    for (int  i = 0; i < nbre_x; i++){
            // Logarithmic grid in x with nbre_x points from min_x to 1
            double x = exp(-log(min_x) / (nbre_x - 1) * i + log(min_x));
            grid_x += std::to_string(x) + " [none], ";
            double u = gpds.at(2).Evaluate(x) / x;
            double uMinus = (gpds.at(2).Evaluate(x) - gpds.at(-2).Evaluate(x)) / x;
            double uPlus = (gpds.at(2).Evaluate(x) + gpds.at(-2).Evaluate(x)) / x;
            Hu += "u: " + std::to_string(u) + " u(+): " + std::to_string(uPlus) + " u(-): " + std::to_string(uMinus) + ", ";
    }
    Hu += "]";
    grid_x += "]";

    PARTONS::Partons::getInstance()->getLoggerManager()->info("main", __func__, grid_x);
    PARTONS::Partons::getInstance()->getLoggerManager()->info("main", __func__, Hu);
}






void CFFNeuralNetwork() {

    // Retrieve service
    PARTONS::DVCSConvolCoeffFunctionService* pDVCSConvolCoeffFunctionService =
            PARTONS::Partons::getInstance()->getServiceObjectRegistry()->getDVCSConvolCoeffFunctionService();

    // Create CFF module with the BaseModuleFactory
    PARTONS::DVCSConvolCoeffFunctionModule* pDVCSCFFModule =
            PARTONS::Partons::getInstance()->getModuleObjectFactory()->newDVCSConvolCoeffFunctionModule(
                    PARTONS::DVCSCFFNN::classId);

    std::string res = "[";
    for(int i = 1; i < 101; i++){
        // Select the replica no. i
        ElemUtils::Parameters parameters(
                PARTONS::DVCSCFFNN::PARAMETER_NAME_REPLICA, i);
        pDVCSCFFModule->configure(parameters);

        // Create the list of kinematics
        PARTONS::List<PARTONS::DVCSConvolCoeffFunctionKinematic> cffKinematicList;
        int nbre_xi = 20;
        double min_xi = 1e-3;
        for (int  j = 0; j < nbre_xi; j++){
                // Logarithmic grid in x with 20 points from 1e-3 to 1
                cffKinematicList.add( PARTONS::DVCSConvolCoeffFunctionKinematic
                        (exp(-log(min_xi) / (nbre_xi - 1) * j + log(min_xi)), -0.3, 2., 2., 2.) );
        }
        // Run computation
        PARTONS::List<PARTONS::DVCSConvolCoeffFunctionResult> cffResult =
                pDVCSConvolCoeffFunctionService->computeManyKinematic(
                        cffKinematicList, pDVCSCFFModule);

        res += "[";
        for (int  j = 0; j < nbre_xi; j++){
                res += "(" + std::to_string(cffResult[j].getKinematic().getXi().getValue()) + ", " +
                    std::to_string(cffResult[j].getResult(PARTONS::GPDType::H).imag())  + "), ";
        }
        res += "], ";
    }
    res += "]";

    // Print results for DVCSCFFModule
    PARTONS::Partons::getInstance()->getLoggerManager()->info("main", __func__,
            res);

    // Remove pointer references
    // Module pointers are managed by PARTONS
    PARTONS::Partons::getInstance()->getModuleObjectFactory()->updateModulePointerReference(
            pDVCSCFFModule, 0);
    pDVCSCFFModule = 0;
}
