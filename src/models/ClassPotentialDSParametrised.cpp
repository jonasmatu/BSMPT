// Copyright (C) 2018  Philipp Basler and Margarete Mühlleitner
// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include "Eigen/IterativeLinearSolvers"
#include <BSMPT/models/SMparam.h> // for C_vev0, C_MassTop, C_g
#include <algorithm>              // for max, copy
#include <cmath>
#include <iostream>               // for operator<<, endl, basic_o...
#include <memory>                 // for allocator_traits<>::value...
#include <stddef.h>               // for std::size_t

#include <BSMPT/models/ClassPotentialDSParametrised.h>
#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/utility/Logger.h>
#include <BSMPT/utility/utility.h>
using namespace Eigen;

/**
 * @file
 * DS for adding a new model class
 */

namespace BSMPT
{
namespace Models
{
/**
 * Here you have to adjust NNeutralHiggs, NChargedHiggs, nPar (number of
 * Lagrangian parameters AFTER using the tadpole conditions), nParCT (number of
 * counterterms) as well as nVEV (number of VEVs for minimization)
 */
Class_DSParametrised::Class_DSParametrised(const ISMConstants &smConstants)
    : Class_Potential_Origin(smConstants)
{
  Model =
      ModelID::ModelIDs::DSPARAMETRISED; // global int constant which will be used to
                                   // tell the program which model is called
  NNeutralHiggs = 1;               // number of neutral Higgs bosons at T = 0
  NChargedHiggs = 0; // number of charged Higgs bosons  at T = 0 (all d.o.f.)

  nPar   = 3; // number of parameters in the tree-Level Lagrangian
  nParCT = 0; // number of parameters in the counterterm potential

  nVEV = 1; // number of VEVs to minimize the potential

  NHiggs = NNeutralHiggs + NChargedHiggs;

  VevOrder.resize(nVEV);
  // Here you have to tell which scalar field gets which VEV.
  VevOrder[0] = 0;

  // Set UseVTreeSimplified to use the tree-level potential defined in
  // VTreeSimplified
  UseVTreeSimplified = true;

  // Set UseVCounterSimplified to use the counterterm potential defined in
  // VCounterSimplified
  UseVCounterSimplified = false;
}

Class_DSParametrised::~Class_DSParametrised()
{
}

/**
 * returns a string which tells the user the chronological order of the
 * counterterms. Use this to complement the legend of the given input file
 */
std::vector<std::string> Class_DSParametrised::addLegendCT() const
{
  std::vector<std::string> labels;
  // labels.push_back("dT");
  // labels.push_back("dlambda");
  // labels.push_back("dmsquared");
  return labels;
}

/**
 * returns a string which tells the user the chronological order of the VEVs and
 * the critical temperature. Use this to complement the legend of the given
 * input file
 */
std::vector<std::string> Class_DSParametrised::addLegendTemp() const
{
  std::vector<std::string> labels;
  labels.push_back("T_c"); // Label for the critical temperature
  labels.push_back("v_c"); // Label for the critical vev
  labels.push_back(
      "v_c/T_c"); // Label for v_c/T_c, you could use xi_c also for example
  // out += "Your VEV order"; // Now you have to put the label for your vevs
  labels.push_back("omega");
  return labels;
}

/**
 * returns a string which tells the user the chronological order of the Triple
 * Higgs couplings. Use this to complement the legend of the given input file
 *
 */
std::vector<std::string> Class_DSParametrised::addLegendTripleCouplings() const
{
  std::vector<std::string> labels;
  std::vector<std::string> particles;
  particles.resize(NHiggs);
  // here you have to define the particle names in the vector particles

  particles[0] = "H";

  std::string out = "Tree_";
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    for (std::size_t j = i; j < NHiggs; j++)
    {
      for (std::size_t k = j; k < NHiggs; k++)
      {
        labels.push_back("Tree_" + particles.at(i) + particles.at(j) +
                         particles.at(k));
        labels.push_back("CT_" + particles.at(i) + particles.at(j) +
                         particles.at(k));
        labels.push_back("CW_" + particles.at(i) + particles.at(j) +
                         particles.at(k));
      }
    }
  }

  return labels;
}

/**
 * returns a string which tells the user the chronological order of the VEVs.
 * Use this to complement the legend of the given input file
 */
std::vector<std::string> Class_DSParametrised::addLegendVEV() const
{
  std::vector<std::string> labels;
  // out = "Your VEV order";
  labels.push_back("omega");
  return labels;
}

/**
 * Reads the string linestr and sets the parameter point
 */
void Class_DSParametrised::ReadAndSet(const std::string &linestr,
                                std::vector<double> &par)
{
  std::stringstream ss(linestr);
  double tmp;

  double lA{0}, llambda3{0}, llambda4{0}, lT0{0};

  if (UseIndexCol)
  {
    ss >> tmp;
  }

  for (int k = 1; k <= 3; k++)
  {
    ss >> tmp;
    if (k == 1)
      lA = tmp;
    else if (k == 2)
      llambda3 = tmp;
    else if (k == 3)
      llambda4 = tmp;
    else if (k == 4)
      lT0 = tmp;
  }
  par[0] = lA;
  par[1] = llambda3;
  par[2] = llambda4;
  par[4] = lT0;

  set_gen(par); // This you have to call so that everything will be set
  return;
}

/**
 * Set Class Object as well as the VEV configuration
 */
void Class_DSParametrised::set_gen(const std::vector<double> &par)
{
  A     = par[0]; // Class member is set accordingly to the input parameters
  lambda3 = par[1]; // Class member is set accordingly to the input parameters
  lambda4 = par[2];
  T0      = par[3];
  vev = std::sqrt(4.0 * A * T0*T0 / lambda4);
  scale = vev;
  vevTreeMin.resize(nVEV);
  vevTree.resize(NHiggs);
  // Here you have to set the vector vevTreeMin. The vector vevTree will then be
  // set by the function MinimizeOrderVEV
  vevTreeMin[0] = vev;
  vevTree       = MinimizeOrderVEV(vevTreeMin);
  if (!SetCurvatureDone) SetCurvatureArrays();
}

/**
 * set your counterterm parameters from the entries of par as well as the
 * entries of Curvature_Higgs_CT_L1 to Curvature_Higgs_CT_L4.
 */
void Class_DSParametrised::set_CT_Pot_Par(const std::vector<double> &par)
{
  // No counterterms
  // they hide in the parametrisation
  Curvature_Higgs_CT_L1[0]          = 0.0;
  Curvature_Higgs_CT_L2[0][0]       = 0.0;
  Curvature_Higgs_CT_L4[0][0][0][0] = 0.0;
}

/**
 * console output of all Parameters
 */
void Class_DSParametrised::write() const
{

  std::stringstream ss;
  ss << "Model = " << Model << std::endl;

  ss << "The parameters are : " << std::endl;
  ss << "lambda3 = " << lambda3 <<  std::endl << "lambda4 = " << lambda4 << std::endl;
  ss << "A = " << A << std::endl << "vev = " << vev << std::endl;
  ss << "T0 = " << T0 << std::endl;

  ss << "The scale is given by mu = " << scale << " GeV " << std::endl;
  Logger::Write(LoggingLevel::Default, ss.str());
}

/**
 * Calculates the counterterms. Here you need to work out the scheme and
 * implement the formulas.
 */
std::vector<double> Class_DSParametrised::calc_CT() const
{

  std::vector<double> parCT;

  if (!SetCurvatureDone)
  {
    std::string retmes = __func__;
    retmes += " was called before SetCurvatureArrays()!\n";
    throw std::runtime_error(retmes);
  }
  if (!CalcCouplingsdone)
  {
    std::string retmes = __func__;
    retmes += " was called before CalculatePhysicalCouplings()!\n";
    throw std::runtime_error(retmes);
  }

  std::vector<double> WeinbergNabla, WeinbergHesse;
  WeinbergNabla = WeinbergFirstDerivative();
  WeinbergHesse = WeinbergSecondDerivative();

  VectorXd NablaWeinberg(NHiggs);
  MatrixXd HesseWeinberg(NHiggs, NHiggs), HiggsRot(NHiggs, NHiggs);
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    NablaWeinberg[i] = WeinbergNabla[i];
    for (std::size_t j = 0; j < NHiggs; j++)
      HesseWeinberg(i, j) = WeinbergHesse.at(j * NHiggs + i);
  }

  // Here you have to use your formulae for the counterterm scheme
  return parCT;
}

void Class_DSParametrised::TripleHiggsCouplings()
{
  if (!SetCurvatureDone) SetCurvatureArrays();
  if (!CalcCouplingsdone) CalculatePhysicalCouplings();

  std::vector<double> HiggsOrder(NHiggs);
  // Here you have to set the vector HiggsOrder. By telling e.g. HiggsOrder[0] =
  // 5 you always want your 6th lightest particle to be the first particle in
  // the vector (which has the index 5 because they are sorted by mass)

  // example for keeping the mass order
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    HiggsOrder[i] = i;
  }

  std::vector<double> TripleDeriv;
  TripleDeriv = WeinbergThirdDerivative();
  std::vector<std::vector<std::vector<double>>> GaugeBasis(
      NHiggs,
      std::vector<std::vector<double>>(NHiggs, std::vector<double>(NHiggs)));
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      for (std::size_t k = 0; k < NHiggs; k++)
      {
        GaugeBasis[i][j][k] =
            TripleDeriv.at(i + j * NHiggs + k * NHiggs * NHiggs);
      }
    }
  }

  MatrixXd HiggsRot(NHiggs, NHiggs);
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      HiggsRot(i, j) = HiggsRotationMatrix[i][j];
    }
  }

  MatrixXd HiggsRotSort(NHiggs, NHiggs);

  for (std::size_t i = 0; i < NHiggs; i++)
  {
    HiggsRotSort.row(i) = HiggsRot.row(HiggsOrder[i]);
  }

  TripleHiggsCorrectionsCWPhysical.resize(NHiggs);
  TripleHiggsCorrectionsTreePhysical.resize(NHiggs);
  TripleHiggsCorrectionsCTPhysical.resize(NHiggs);
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    TripleHiggsCorrectionsTreePhysical[i].resize(NHiggs);
    TripleHiggsCorrectionsCWPhysical[i].resize(NHiggs);
    TripleHiggsCorrectionsCTPhysical[i].resize(NHiggs);
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      TripleHiggsCorrectionsCWPhysical[i][j].resize(NHiggs);
      TripleHiggsCorrectionsTreePhysical[i][j].resize(NHiggs);
      TripleHiggsCorrectionsCTPhysical[i][j].resize(NHiggs);
    }
  }

  for (std::size_t i = 0; i < NHiggs; i++)
  {
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      for (std::size_t k = 0; k < NHiggs; k++)
      {
        TripleHiggsCorrectionsCWPhysical[i][j][k]   = 0;
        TripleHiggsCorrectionsTreePhysical[i][j][k] = 0;
        TripleHiggsCorrectionsCTPhysical[i][j][k]   = 0;
        for (std::size_t l = 0; l < NHiggs; l++)
        {
          for (std::size_t m = 0; m < NHiggs; m++)
          {
            for (std::size_t n = 0; n < NHiggs; n++)
            {
              double RotFac =
                  HiggsRotSort(i, l) * HiggsRotSort(j, m) * HiggsRotSort(k, n);
              TripleHiggsCorrectionsCWPhysical[i][j][k] +=
                  RotFac * GaugeBasis[l][m][n];
              TripleHiggsCorrectionsTreePhysical[i][j][k] +=
                  RotFac * LambdaHiggs_3[l][m][n];
              TripleHiggsCorrectionsCTPhysical[i][j][k] +=
                  RotFac * LambdaHiggs_3_CT[l][m][n];
            }
          }
        }
      }
    }
  }
}

void Class_DSParametrised::SetCurvatureArrays()
{
  /*
   *  Here you have to set the vectors
   *  Curvature_Higgs_L1,Curvature_Higgs_L2,Curvature_Higgs_L3,Curvature_Higgs_L4
   *  Curvature_Gauge_G2H2
   *  Curvature_Quark_F2H1, Curvature_Lepton_F2H1
   *  as described in the potential in the paper.
   */

  initVectors();
  SetCurvatureDone = true;
  for (std::size_t i = 0; i < NHiggs; i++)
    HiggsVev[i] = vevTree[i];

  // We set every couplings to zero in order to get only the
  // tree level potential.
  Curvature_Higgs_L2[0][0]       = 0.0;
  Curvature_Higgs_L3[0][0][0]    = 0.0;
  Curvature_Higgs_L4[0][0][0][0] = 0.0;
  // Curvature_Higgs_L2[0][0]       = msq;
  // Curvature_Higgs_L3[0][0][0]    = lambda3;
  // Curvature_Higgs_L4[0][0][0][0] = lambda4;

  Curvature_Gauge_G2H2[0][0][0][0] = 0;

  Curvature_Quark_F2H1[1][0][0] = 0;
  Curvature_Quark_F2H1[0][1][0] = 0;
}

bool Class_DSParametrised::CalculateDebyeSimplified()
{
  /*
   * Use this function if you calculated the Debye corrections to the Higgs mass
   * matrix and implement your formula here and return true. The vector is given
   * by DebyeHiggs[NHiggs][NHiggs]
   */
  // No thermal masses
  return true;
}

bool Class_DSParametrised::CalculateDebyeGaugeSimplified()
{
  /*
   * Use this function if you calculated the Debye corrections to the gauge mass
   * matrix and implement your formula here and return true. The vector is given
   * by DebyeGauge[NGauge][NGauge]
   */

  // We don't want any debey masses
  return true;
}
double Class_DSParametrised::VTreeSimplified(const std::vector<double> &v, double T) const
{
  if (not UseVTreeSimplified) return 0;
  double res = 0;

  double vIn = v[0];
  res = A * (T*T - T0*T0)*vIn*vIn  - lambda3 * T * vIn + lambda4 / 4.0 * std::pow(vIn,4);
  return res;
}
double Class_DSParametrised::VTreeSimplified(const std::vector<double> &v) const
{
  if (not UseVTreeSimplified) return 0;
  double res = 0;
  return res;
}
// Overwrite the default implementation of the effective potential
// This allows for a temperature dependent tree level potential
double Class_DSParametrised::VEff(const std::vector<double> &v,
                                    double Temp,
                                    int diff,
                                    int Order) const
{
  if (v.size() != nVEV and v.size() != NHiggs)
  {
    std::string ErrorString =
        std::string("You have called ") + std::string(__func__) +
        std::string(
            " with an invalid vev configuration. Your vev is of dimension ") +
        std::to_string(v.size()) + std::string(" and it should be ") +
        std::to_string(NHiggs) + std::string(".");
    throw std::runtime_error(ErrorString);
  }
  if (v.size() == nVEV and nVEV != NHiggs)
  {
    std::stringstream ss;
    ss << __func__
       << " is being called with a wrong sized vev configuration. It "
          "has the dimension of "
       << nVEV << " while it should have " << NHiggs
       << ". For now this is transformed but please fix this to reduce "
          "the runtime."
       << std::endl;
    Logger::Write(LoggingLevel::Default, ss.str());
    std::vector<double> Transformedv;
    Transformedv = MinimizeOrderVEV(v);
    return VEff(Transformedv, Temp, diff);
  }

  double resOut = 0;
  resOut        = VTree(v, diff);
  // if (Order != 0 and not UseTreeLevel)
  // {
  //   resOut += CounterTerm(v, diff);
  //   resOut += V1Loop(v, Temp, diff);
  // }
  // for(std::size_t i=0;i<NHiggs;i++) resOut +=
  // DebyeHiggs[i][i]*0.5*std::pow(v.at(i),2)*std::pow(Temp,2);
  return resOut;
}
double Class_DSParametrised::VCounterSimplified(const std::vector<double> &v) const
{
  if (not UseVCounterSimplified) return 0;
  double res = 0.0;
  return res;
}

void Class_DSParametrised::Debugging(const std::vector<double> &input,
                               std::vector<double> &output) const
{
  (void)input;
  (void)output;
}

} // namespace Models
} // namespace BSMPT
