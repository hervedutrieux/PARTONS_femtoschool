#include <ElementaryUtils/logger/CustomException.h>
#include <ElementaryUtils/logger/LoggerManager.h>
#include <partons/Partons.h>
#include <partons/services/automation/AutomationService.h>
#include <partons/ServiceObjectRegistry.h>
#include <QtCore/qcoreapplication.h>
#include <string>
#include <vector>

// New examples developed for the First International School of Hadron Femtography
// Not official or endorsed by the PARTONS collaboration
// Herve Dutrieux - September 19, 2024

void extract_xdep_GPD(double xi, double t, double mu2, double min_x = 1e-4, int nbre_x = 100) ;
std::map<int, double> GPDModel(double const& x, double const&) ;
void using_GPDEvolution_high_perf(double mu2F, double mu2I, double min_x = 1e-4, int nbre_x = 100) ;
void CFFNeuralNetwork();

// XML Parser
std::vector<std::string> parseArguments(int argc, char** argv) {
    std::vector<std::string> xmlScenarios(argc - 1);

    for (unsigned int i = 1; i < argc; i++) {
        xmlScenarios[i - 1] = argv[i];
    }

    return xmlScenarios;
}


int main(int argc, char** argv) {

    // Init Qt4
    QCoreApplication a(argc, argv);
    PARTONS::Partons* pPartons = 0;

    try {

        // Init PARTONS application
        pPartons = PARTONS::Partons::getInstance();
        pPartons->init(argc, argv);

        // XML scenario
        if (argc <= 1) {

            throw ElemUtils::CustomException("main", __func__,
                    "Missing argument, please provide one or more than one XML scenario file.");
        }

        // Parse arguments to retrieve xml file path list.
        std::vector<std::string> xmlScenarioFilePathList = parseArguments(argc,
                argv);

        // Retrieve automation service parse scenario xml file and play it.
        PARTONS::AutomationService* pAutomationService =
                pPartons->getServiceObjectRegistry()->getAutomationService();

        for (unsigned int i = 0; i < xmlScenarioFilePathList.size(); i++) {
            PARTONS::Scenario* pScenario = pAutomationService->parseXMLFile(
                    xmlScenarioFilePathList[i]);
            pAutomationService->playScenario(pScenario);
        }

    }
    // Exceptions
    catch (const ElemUtils::CustomException &e) {
        pPartons->getLoggerManager()->error(e);
        if (pPartons) { pPartons->close(); }
    }
    catch (const std::exception &e) {
        pPartons->getLoggerManager()->error("main", __func__, e.what());
        if (pPartons) { pPartons->close(); }
    }

    // Close PARTONS application properly
    if (pPartons) {
        pPartons->close();
    }

    return 0;
}
