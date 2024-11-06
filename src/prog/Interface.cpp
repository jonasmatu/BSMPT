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


// vector<double> getGWParams(double gDM, double yDM, double vev)
extern "C" {double getNuclTemp(double A, double lambda3, double lambda4, double T0)
try
{

  // auto argparser         = prepare_parser();
  // argparser.add_input(convert_input(argc, argv));
  // const CLIOptions args(argparser);
  // if (not args.good())
  // {
  //   return EXIT_FAILURE;
  // }

  // std::ifstream infile(args.inputfile);
  // if (!infile.good())
  // {
  //   Logger::Write(LoggingLevel::Default,
  //                 "Input file " + args.inputfile + " not found ");
  //   return EXIT_FAILURE;
  // }

  // Logger::Write(LoggingLevel::ProgDetailed, "Found file");

  // Parameters:
  const auto SMConstants = GetSMConstants();
  // BSMPT::ModelID::ModelIDs Model{ModelID::ModelIDs::CONFORMALDM};
  BSMPT::ModelID::ModelIDs Model{ModelID::ModelIDs::DSPARAMETRISED};
  // int firstline{0}, lastline{0};
  double templow{0}, temphigh{50};
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
  // int num_points    = args.lastline - args.firstline + 1;

  // output contents storage
  // std::vector<std::stringstream> output_contents;
  // output_contents.resize(num_points); // reserve one row per point
  std::vector<std::string> transition_history;
  std::vector<std::string> legend;

  // double gDM = 0.9;
  // double yDM = 0.5;
  // double vev = 2000.0;
  // linestr = std::to_string(gDM) + "\t" + std::to_string(yDM) + "\t" + std::to_string(vev);
  linestr = std::to_string(A) + "\t" + std::to_string(lambda3) + "\t" + std::to_string(lambda4) + "\t" + std::to_string(T0);
  linestr_store = linestr;
  modelPointer->setUseIndexCol(linestr_store);

  // convert parameters to string:
  

  std::pair<std::vector<double>, std::vector<double>> parameters =
    modelPointer->initModel(linestr);

  auto start = std::chrono::high_resolution_clock::now();

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

  TransitionTracer trans(input);

  auto time = std::chrono::duration_cast<std::chrono::milliseconds>(
		 std::chrono::high_resolution_clock::now() - start) .count() / 1000.;

  BSMPT::Logger::Write(BSMPT::LoggingLevel::ProgDetailed,
		       "Took\t" + std::to_string(time) + " seconds.\n");

  auto output = trans.output_store;
  std::cout << "Finished calculation!!" << std::endl;


  // prepare return array
  int i = 0; // asume only one coexisting phase
  // This contains: alpha, betaH, Tperc (or tranistion temp), 
  // std::vector<double> res = {output.vec_gw_data.at(i).alpha.value_or(EmptyValue),
  // 			     output.vec_gw_data.at(i).beta_over_H.value_or(EmptyValue),
  // 			     output.vec_gw_data.at(i).trans_temp.value_or(EmptyValue)};
  
  double res = output.vec_trans_data.at(i).nucl_temp.value_or(EmptyValue);
  return res;
  // return EXIT_SUCCESS;
}
catch (int)
{
  // std::vector<double> res = {0,0,0};
  double res = 0;
  return res;
  // return EXIT_SUCCESS;
}
catch (exception &e)
{
  Logger::Write(LoggingLevel::Default, e.what());
  // std::vector<double> res = {0,0,0};
  double res = 0;
  return res;
  // return EXIT_FAILURE;
}
}
