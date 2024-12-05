// Copyright (C) 2024 Lisa Biermann, Margarete Mühlleitner, Rui Santos, João
// Viana SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and
// Jonas Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

/**
 * @file
 * This program calculates properties of graviational waves sourced by phase
 * transitions
 *
 */

#include "BSMPT/minimum_tracer/minimum_tracer.h"
#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/transition_tracer/transition_tracer.h>
#include <BSMPT/utility/Logger.h>
#include <BSMPT/utility/parser.h>
#include <BSMPT/utility/utility.h>
#include <Eigen/Dense>
#include <iostream>
#include <stdlib.h> // for atoi, EXIT_FAILURE
#include <string>   // for string, operator<<
#include <vector>   // for vector

using namespace std;
using namespace BSMPT;


// extern "C" {double getNuclTemp(double D, double A, double lambda, double T0)
extern "C" {double getNuclTempFullPot(double gDS, double lambda, double vev, double Tmax)
try
{
  // Parameters:
  const auto SMConstants = GetSMConstants();
  // BSMPT::ModelID::ModelIDs Model{ModelID::ModelIDs::DARKPHOTON};
  BSMPT::ModelID::ModelIDs Model{ModelID::ModelIDs::DARKPHOTON};
  // int firstline{0}, lastline{0};
  double templow{0}, temphigh{Tmax};
  double UserDefined_vwall   = 0.95;
  double UserDefined_epsturb = 0.1;
  int MaxPathIntegrations    = 7;
  // std::string inputfile, outputfile;
  bool UseGSL{Minimizer::UseGSLDefault};
  bool UseCMAES{Minimizer::UseLibCMAESDefault};
  bool UseNLopt{Minimizer::UseNLoptDefault};
  int WhichMinimizer{Minimizer::WhichMinimizerDefault};
  bool UseMultithreading{false};
  int UseMultiStepPTMode{-1};
  int CheckEWSymmetryRestoration{1};
  double perc_prbl{.71};
  double compl_prbl{.01};
  int num_check_pts{10};
  int CheckNLOStability{1};
  int WhichTransitionTemperature{
      3}; // 1 = nucl_approx, 2 = nucl, 3 = perc, 4 = compl
  

  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(Model, SMConstants);

  Logger::Write(LoggingLevel::ProgDetailed, "Created modelpointer ");

  std::string linestr, linestr_store;
  int linecounter   = 1;
  std::size_t count = 0;

  std::vector<std::string> transition_history;
  std::vector<std::string> legend;

  // Prepare input for the model
  // linestr = std::to_string(D) + "\t" + std::to_string(A) + "\t" + std::to_string(lambda) + "\t" + std::to_string(T0);
  linestr = std::to_string(gDS) + "\t" + std::to_string(lambda) + "\t" + std::to_string(vev);
  
  linestr_store = linestr;
  modelPointer->setUseIndexCol(linestr_store);

  std::pair<std::vector<double>, std::vector<double>> parameters =
    modelPointer->initModel(linestr);

  auto start = std::chrono::high_resolution_clock::now();

  CheckEWSymmetryRestoration = 0;
  CheckNLOStability = 0;

  user_input input{modelPointer,
		   templow,
		   temphigh,
		   UserDefined_vwall,
		   perc_prbl,
		   compl_prbl,
		   UserDefined_epsturb,
		   MaxPathIntegrations,
		   UseMultiStepPTMode,
		   num_check_pts,
		   CheckEWSymmetryRestoration,
		   CheckNLOStability,
		   WhichMinimizer,
		   UseMultithreading,
		   true,
		   WhichTransitionTemperature};

  // Do the calculation of the phase transition
  TransitionTracer trans(input);

  auto time = std::chrono::duration_cast<std::chrono::milliseconds>(
		 std::chrono::high_resolution_clock::now() - start) .count() / 1000.;

  BSMPT::Logger::Write(BSMPT::LoggingLevel::ProgDetailed,
		       "Took\t" + std::to_string(time) + " seconds.\n");

  auto output = trans.output_store;

  // Check if found a phase transition!
  if (output.vec_trans_data.size() == 0) {
    // Verbose output, maybe better to comment out
    std::cout << "Did not find phase transition for:" << std::endl;
    std::cout << "gDS    = " << gDS << std::endl;
    std::cout << "lambda = " << lambda << std::endl;
    std::cout << "vev    = " << vev << std::endl;
    std::cout << "Tmax   = " << Tmax << std::endl;
    return 0.0;
  }
  
  // Return the nucleation temperature
  double res = output.vec_trans_data.at(0).nucl_temp.value_or(EmptyValue);
  return res;
}
catch (int)
{
  double res = 0;
  return res;
  // return EXIT_SUCCESS;
}
catch (exception &e)
{
  Logger::Write(LoggingLevel::Default, e.what());
  std::cout << "Encountered exception for: " << std::endl;
  std::cout << "Did not find phase transition for:" << std::endl;
  std::cout << "gDS    = " << gDS << std::endl;
  std::cout << "lambda = " << lambda << std::endl;
  std::cout << "vev    = " << vev << std::endl;
  std::cout << "Tmax   = " << Tmax << std::endl;
  double res = 0;
  return res;
}
}

extern "C" {double getNuclTempParametrisedPot(double D, double A, double lambda, double T0, double Tmax)
try
{
  // Parameters:
  const auto SMConstants = GetSMConstants();
  BSMPT::ModelID::ModelIDs Model{ModelID::ModelIDs::DSPARAMETRISED};
  // int firstline{0}, lastline{0};
  double templow{0}, temphigh{Tmax};
  double UserDefined_vwall   = 0.95;
  double UserDefined_epsturb = 0.1;
  int MaxPathIntegrations    = 7;
  // std::string inputfile, outputfile;
  bool UseGSL{Minimizer::UseGSLDefault};
  bool UseCMAES{Minimizer::UseLibCMAESDefault};
  bool UseNLopt{Minimizer::UseNLoptDefault};
  int WhichMinimizer{Minimizer::WhichMinimizerDefault};
  bool UseMultithreading{false};
  int UseMultiStepPTMode{-1};
  int CheckEWSymmetryRestoration{1};
  double perc_prbl{.71};
  double compl_prbl{.01};
  int num_check_pts{10};
  int CheckNLOStability{1};
  int WhichTransitionTemperature{
      3}; // 1 = nucl_approx, 2 = nucl, 3 = perc, 4 = compl
  

  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(Model, SMConstants);

  Logger::Write(LoggingLevel::ProgDetailed, "Created modelpointer ");

  std::string linestr, linestr_store;
  int linecounter   = 1;
  std::size_t count = 0;

  std::vector<std::string> transition_history;
  std::vector<std::string> legend;

  // Prepare input for the model
  linestr = std::to_string(D) + "\t" + std::to_string(A) + "\t" + std::to_string(lambda) + "\t" + std::to_string(T0);
  
  linestr_store = linestr;
  modelPointer->setUseIndexCol(linestr_store);

  std::pair<std::vector<double>, std::vector<double>> parameters =
    modelPointer->initModel(linestr);

  auto start = std::chrono::high_resolution_clock::now();

  CheckEWSymmetryRestoration = 0;
  CheckNLOStability = 0;

  user_input input{modelPointer,
		   templow,
		   temphigh,
		   UserDefined_vwall,
		   perc_prbl,
		   compl_prbl,
		   UserDefined_epsturb,
		   MaxPathIntegrations,
		   UseMultiStepPTMode,
		   num_check_pts,
		   CheckEWSymmetryRestoration,
		   CheckNLOStability,
		   WhichMinimizer,
		   UseMultithreading,
		   true,
		   WhichTransitionTemperature};

  // Do the calculation of the phase transition
  TransitionTracer trans(input);

  auto time = std::chrono::duration_cast<std::chrono::milliseconds>(
		 std::chrono::high_resolution_clock::now() - start) .count() / 1000.;

  BSMPT::Logger::Write(BSMPT::LoggingLevel::ProgDetailed,
		       "Took\t" + std::to_string(time) + " seconds.\n");

  auto output = trans.output_store;

  // Check if found a phase transition 
  if (output.vec_trans_data.size() == 0) {
    std::cout << "Did not find phase transition for:" << std::endl;
    std::cout << "D      = " << D << std::endl;
    std::cout << "A      = " << A << std::endl;
    std::cout << "lambda = " << lambda << std::endl;
    std::cout << "T0     = " << T0 << std::endl;
    std::cout << "Tmax   = " << Tmax << std::endl;
    return 0.0;
  }

  // Return the nucleation temperature
  double res = output.vec_trans_data.at(0).nucl_temp.value_or(EmptyValue);
  return res;
}
catch (int)
{
  double res = 0;
  return res;
  // return EXIT_SUCCESS;
}
catch (exception &e)
{
  Logger::Write(LoggingLevel::Default, e.what());
  std::cout << "Encountered exception for:" << std::endl;
  std::cout << "D      = " << D << std::endl;
  std::cout << "A      = " << A << std::endl;
  std::cout << "lambda = " << lambda << std::endl;
  std::cout << "T0     = " << T0 << std::endl;
  std::cout << "Tmax   = " << Tmax << std::endl;
  double res = 0;
  return res;
}
}
