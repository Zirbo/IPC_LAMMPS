#include "printPotential.hpp"

#include <sys/stat.h>

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include <Eigen/Dense>


// #define DEBUG_ORIENT
// #define DEBUG_MAIN

void
PotentialForLammps::initFromEpsilons(std::string const& inputFileName)
{
  // read input.in file
  std::ifstream inputFile(inputFileName);
  if (inputFile.fail()) {
    std::runtime_error("Input file could not be opened. Aborting.");
  }
  std::string modelName;
  inputFile >> modelName;
  inputFile >> delta;
  inputFile >> eccentricity_p1;
  inputFile >> radius_p1;
  inputFile >> e_BB >> e_Bs1 >> e_s1s1;
  inputFile >> e_min;
  if (symmetry == Symmetry::JANUS) {
    e_Bs2 = 0;
    e_s2s2 = 0;
    e_s1s2 = 0;
  } else if (symmetry == Symmetry::SYMMETRIC) {
    e_Bs2 = e_Bs1;
    e_s1s2 = e_s1s1;
    e_s2s2 = e_s1s1;
  } else if (symmetry == Symmetry::ASYMMETRIC) {
    inputFile >> eccentricity_p2;
    inputFile >> radius_p2;
    inputFile >> e_Bs2 >> e_s1s2 >> e_s2s2;
  }

  inputFile.close();
}

void
PotentialForLammps::readContactValues(std::string const& inputFileName)
{
  // read input.in file
  std::ifstream inputFile(inputFileName);
  if (inputFile.fail()) {
    std::runtime_error("Input file could not be opened. Aborting.");
  }
  std::string modelName;
  inputFile >> modelName;
  inputFile >> delta;
  inputFile >> eccentricity_p1;
  inputFile >> radius_p1;
  inputFile >> vEE >> vEP1 >> vP1P1;
  if (symmetry == Symmetry::ASYMMETRIC) {
    inputFile >> eccentricity_p2;
    inputFile >> radius_p2;
    inputFile >> vEP2 >> vP1P2 >> vP2P2;
  } else {
    vEP2 = 0;
    vP1P2 = 0;
    vP2P2 = 0;
  }
  inputFile.close();
}

PotentialForLammps::PotentialForLammps(const std::string& inputFileName,
                                       const Symmetry symmetry,
                                       const Colloid colloid,
                                       bool startFromContactValues)
  : symmetry{ symmetry }
{
  if (symmetry == Symmetry::NONE) {
    throw std::runtime_error("Symmetry and/or colloid type not specified");
  }

  if (startFromContactValues) {
    readContactValues(inputFileName);
  } else {
    initFromEpsilons(inputFileName);
  }

  // define geometry
  HSdiameter = 1.0;
  fakeHScoefficient = 500;
  fakeHSexponent = 15;
  samplingStep = 1.0e-05;
  cutoffValue = 1.0e+06;

  colloidRadius = 0.5 * (HSdiameter + delta);

  if (symmetry != Symmetry::ASYMMETRIC) {
    eccentricity_p2 = eccentricity_p1;
    radius_p2 = radius_p1;
  }

  // interactionRange -> furthest point from center with nonzero interaction
  double range_p1 = radius_p1 + eccentricity_p1;
  double range_p2 = radius_p2 + eccentricity_p2;
  interactionRange = 2 * std::max(std::max(range_p1, range_p2), colloidRadius);
#ifdef DEBUG_MAIN
  std::cout << "interaction range = 2 x max(" << colloidRadius << ", "
            << range_p1 << ", " << range_p2 << ") = " << interactionRange
            << "\n\n";
#endif

  // sanity checks
  if (eccentricity_p1 > 0.5 * HSdiameter || radius_p1 <= 0. ||
      eccentricity_p1 <= 0.) {
    throw std::runtime_error("First patch: invalid geometry");
  }
  if (eccentricity_p2 > 0.5 * HSdiameter || radius_p2 <= 0. ||
      eccentricity_p2 <= 0.) {
    throw std::runtime_error("Second patch: invalid geometry");
  }
  if (delta < 0.) {
    throw std::runtime_error("delta cannot be negative");
  }

  if (colloid == Colloid::IPC) {
    // check IPC geometry requirements
    const double diff_1 = (range_p1 / range_p2) - 1.;
    const double diff_2 = (range_p1 / colloidRadius) - 1.;

    if (std::fabs(diff_1) > 1e-4 || std::fabs(diff_2) > 1e-4) {
      throw std::runtime_error("Invalid IPC geometry");
    }
  }

  // define potential coefficients
  if (startFromContactValues) {
    computeEpsilonsFromContactValues();
  }

  computeSiteSitePotentials();
}

static double
computeOmega(double Ra, double Rb, double rab)
{
  // BKL paper, formula 18
  if (rab > Ra + Rb)
    return 0.;
  else if (rab <= std::fabs(Ra - Rb))
    return 8. * std::pow(std::min(Ra, Rb), 3);
  else {
    const double tempSum = (Ra * Ra - Rb * Rb) / (2. * rab);
    return 2. *
           ((2. * Ra + tempSum + rab / 2.) * pow(Ra - tempSum - rab / 2., 2) +
            (2. * Rb - tempSum + rab / 2.) * pow(Rb + tempSum - rab / 2., 2));
  }
}

void
PotentialForLammps::computeEpsilonsFromContactValues()
{
  bool noPPinEE =
    (radius_p1 < .5*HSdiameter) && (radius_p2 < .5*HSdiameter);
  bool noCPinEE =
    pow(colloidRadius + radius_p1, 2) < pow(HSdiameter, 2) + pow(eccentricity_p2, 2)
 && pow(colloidRadius + radius_p2, 2) < pow(HSdiameter, 2) + pow(eccentricity_p1, 2) ;
  bool noPPinEP =
    pow(radius_p1 + radius_p1, 2) < pow(HSdiameter - eccentricity_p1, 2) + pow(eccentricity_p1, 2)
 && pow(radius_p1 + radius_p2, 2) < pow(HSdiameter - eccentricity_p1, 2) + pow(eccentricity_p2, 2)
 && pow(radius_p2 + radius_p1, 2) < pow(HSdiameter - eccentricity_p2, 2) + pow(eccentricity_p1, 2)
 && pow(radius_p2 + radius_p2, 2) < pow(HSdiameter - eccentricity_p2, 2) + pow(eccentricity_p2, 2);
  reducedMode = noPPinEE && noCPinEE && noPPinEP;
  reducedMode &= false;
  if (reducedMode) {
    computeEpsilonsFromContactValuesReduced();
  } else {
    computeEpsilonsFromContactValuesGeneral();
  }
}

void
PotentialForLammps::computeEpsilonsFromContactValuesReduced()
{
  e_min = 1.0;

  // overlap volumes at contact
  double fBB = computeOmega(colloidRadius, colloidRadius, HSdiameter);
  double fBs1 =
    computeOmega(colloidRadius, radius_p1, HSdiameter - eccentricity_p1);
  double fs1s1 =
    computeOmega(radius_p1, radius_p1, HSdiameter - 2. * eccentricity_p1);
  double fBs2 =
    computeOmega(colloidRadius, radius_p2, HSdiameter - eccentricity_p2);
  double fs1s2 = computeOmega(
    radius_p1, radius_p2, HSdiameter - eccentricity_p1 - eccentricity_p2);
  double fs2s2 =
    computeOmega(radius_p2, radius_p2, HSdiameter - 2. * eccentricity_p2);

  // compute coefficients
  if (fBB == 0.) {
    std::cout << "WARNING: fBB is zero. Did you set delta to zero?\n"
              << "Overriding e_BB and vEE to also be zero\n\n";
    e_BB = 0;
    vEE = 0;
  } else {
    e_BB = vEE / fBB;
  }
  e_Bs1 = (vEP1 - vEE) / fBs1;
  e_s1s1 = (vP1P1 - vEE - 2. * fBs1 * e_Bs1) / fs1s1;

  if (symmetry == Symmetry::JANUS) {
    e_Bs2 = 0;
    e_s2s2 = 0;
    e_s1s2 = 0;
  } else if (symmetry == Symmetry::SYMMETRIC) {
    e_Bs2 = e_Bs1;
    e_s1s2 = e_s1s1;
    e_s2s2 = e_s1s1;
  } else if (symmetry == Symmetry::ASYMMETRIC) {
    e_Bs2 = (vEP2 - vEE) / fBs2;
    e_s1s2 = (vP1P2 - vEE - fBs1 * e_Bs1 - fBs2 * e_Bs2) / fs1s2;
    e_s2s2 = (vP2P2 - vEE - 2. * fBs2 * e_Bs2) / fs2s2;
  }
}

static double real_dist(double x, double y)
{
  return std::sqrt(x * x + y * y);
}

void
PotentialForLammps::computeEpsilonsFromContactValuesGeneral()
{
  e_min = 1.0;
  if (computeOmega(colloidRadius, colloidRadius, HSdiameter) == 0.) {
    std::cout << "WARNING: fBB is zero. Did you set delta to zero?\n"
              << "Overriding e_BB and vEE to also be zero\n\n";
  }

  if (symmetry == Symmetry::JANUS) {
    throw std::runtime_error("General case not implemented for Janus");
  } else if (symmetry == Symmetry::SYMMETRIC) {
    vEP2 = vEP1;
    vP1P2 = vP1P1;
    vP2P2 = vP1P1;
  }

  Eigen::VectorXd V(6);
  V << vEE, vEP1, vEP2, vP1P1, vP1P2, vP2P2;
  // sequence: BB Bs1 Bs2 s1s1 s1s2 s2s2
  Eigen::MatrixXd f(6,6);
  // = vEE
  f(0,0) = computeOmega(colloidRadius, colloidRadius, real_dist(HSdiameter, 0));
  f(0,1) = computeOmega(colloidRadius, radius_p1, real_dist(HSdiameter, eccentricity_p1)) * 2.;
  f(0,2) = computeOmega(colloidRadius, radius_p2, real_dist(HSdiameter, eccentricity_p2)) * 2.;
  f(0,3) = computeOmega(radius_p1, radius_p1, real_dist(HSdiameter, 0));
  f(0,4) = computeOmega(radius_p1, radius_p2, real_dist(HSdiameter, eccentricity_p1 + eccentricity_p2)) * 2.;
  f(0,5) = computeOmega(radius_p2, radius_p2, real_dist(HSdiameter, 0));
  // = VEP1
  f(1,0) = computeOmega(colloidRadius, colloidRadius, real_dist(HSdiameter, 0));
  f(1,1) = computeOmega(colloidRadius, radius_p1, real_dist(HSdiameter - eccentricity_p1, 0))
         + computeOmega(colloidRadius, radius_p1, real_dist(HSdiameter, eccentricity_p1));
  f(1,2) = computeOmega(colloidRadius, radius_p2, real_dist(HSdiameter + eccentricity_p2, 0))
         + computeOmega(colloidRadius, radius_p2, real_dist(HSdiameter, -eccentricity_p2));;
  f(1,3) = computeOmega(radius_p1, radius_p1, real_dist(HSdiameter - eccentricity_p1, eccentricity_p1));
  f(1,4) = computeOmega(radius_p1, radius_p2, real_dist(HSdiameter + eccentricity_p2, eccentricity_p1))
         + computeOmega(radius_p1, radius_p2, real_dist(HSdiameter - eccentricity_p1, -eccentricity_p2));
  f(1,5) = computeOmega(radius_p2, radius_p2, real_dist(HSdiameter + eccentricity_p2, -eccentricity_p2));
  // = VEP2
  f(2,0) = computeOmega(colloidRadius, colloidRadius, real_dist(HSdiameter, 0));
  f(2,1) = computeOmega(colloidRadius, radius_p1, real_dist(HSdiameter + eccentricity_p1, 0))
         + computeOmega(colloidRadius, radius_p1, real_dist(HSdiameter, eccentricity_p1));
  f(2,2) = computeOmega(colloidRadius, radius_p2, real_dist(HSdiameter - eccentricity_p2, 0))
         + computeOmega(colloidRadius, radius_p2, real_dist(HSdiameter, eccentricity_p2));;
  f(2,3) = computeOmega(radius_p1, radius_p1, real_dist(HSdiameter + eccentricity_p1, eccentricity_p1));
  f(2,4) = computeOmega(radius_p1, radius_p2, real_dist(HSdiameter - eccentricity_p2, eccentricity_p1))
         + computeOmega(radius_p1, radius_p2, real_dist(HSdiameter + eccentricity_p1, -eccentricity_p2));
  f(2,5) = computeOmega(radius_p2, radius_p2, real_dist(HSdiameter - eccentricity_p2, -eccentricity_p2));
  // = VP1P1
  f(3,0) = computeOmega(colloidRadius, colloidRadius, real_dist(HSdiameter, 0));
  f(3,1) = computeOmega(colloidRadius, radius_p1, real_dist(HSdiameter - eccentricity_p1, 0)) * 2.;
  f(3,2) = computeOmega(colloidRadius, radius_p2, real_dist(HSdiameter + eccentricity_p2, 0)) * 2.;
  f(3,3) = computeOmega(radius_p1, radius_p1, real_dist(HSdiameter -2*eccentricity_p1, 0));
  f(3,4) = computeOmega(radius_p1, radius_p2, real_dist(eccentricity_p2 + HSdiameter - eccentricity_p1, 0)) * 2.;
  f(3,5) = computeOmega(radius_p2, radius_p2, real_dist(HSdiameter +2*eccentricity_p2, 0));
  // = VP1P2
  f(4,0) = computeOmega(colloidRadius, colloidRadius, real_dist(HSdiameter, 0));
  f(4,1) = computeOmega(colloidRadius, radius_p1, real_dist(HSdiameter + eccentricity_p1, 0))
         + computeOmega(colloidRadius, radius_p1, real_dist(HSdiameter - eccentricity_p1, 0));
  f(4,2) = computeOmega(colloidRadius, radius_p2, real_dist(HSdiameter + eccentricity_p2, 0))
         + computeOmega(colloidRadius, radius_p2, real_dist(HSdiameter - eccentricity_p2, 0));
  f(4,3) = computeOmega(radius_p1, radius_p1, real_dist(HSdiameter, 0));
  f(4,4) = computeOmega(radius_p1, radius_p2, real_dist(eccentricity_p1 + HSdiameter + eccentricity_p2, 0))
         + computeOmega(radius_p1, radius_p2, real_dist(HSdiameter - eccentricity_p1 - eccentricity_p2, 0));
  f(4,5) = computeOmega(radius_p2, radius_p2, real_dist(HSdiameter, 0));
  // = VP2P2
  f(5,0) = computeOmega(colloidRadius, colloidRadius, real_dist(HSdiameter, 0));
  f(5,1) = computeOmega(colloidRadius, radius_p1, real_dist(HSdiameter + eccentricity_p1, 0)) * 2.;
  f(5,2) = computeOmega(colloidRadius, radius_p2, real_dist(HSdiameter - eccentricity_p2, 0)) * 2.;
  f(5,3) = computeOmega(radius_p1, radius_p1, real_dist(HSdiameter +2*eccentricity_p1, 0));
  f(5,4) = computeOmega(radius_p1, radius_p2, real_dist(eccentricity_p1 + HSdiameter - eccentricity_p2, 0)) * 2.;
  f(5,5) = computeOmega(radius_p2, radius_p2, real_dist(HSdiameter -2*eccentricity_p2, 0));

  // solve fe = V
  Eigen::VectorXd e(6);
  //e = f.colPivHouseholderQr().solve(V);
  e = f.completeOrthogonalDecomposition().solve(V);

  // port solutions
  e_BB   = e[0];
  e_Bs1  = e[1];
  e_Bs2  = e[2];
  e_s1s1 = e[3];
  e_s1s2 = e[4];
  e_s2s2 = e[5];
}

static double
computeOmegaRadialDerivative(double Ra, double Rb, double rab)
{
  // BKL paper, derivative of formula 18
  if (rab >= Ra + Rb || rab <= fabs(Ra - Rb))
    return 0.;
  else {
    const double tempSum = (Ra * Ra - Rb * Rb) / (2. * rab);
    const double tempSumMinus = tempSum - rab / 2.;
    const double tempSumPlus = tempSum + rab / 2.;
    return (6. / rab) *
           (tempSumMinus * (Ra - tempSumPlus) * (Ra + tempSumPlus) -
            tempSumPlus * (Rb - tempSumMinus) * (Rb + tempSumMinus));
  }
}

void
PotentialForLammps::computeSiteSitePotentials()
{
  const size_t potentialSteps = size_t(interactionRange / samplingStep) + 2;

  uHS.resize(potentialSteps);
  uBB.resize(potentialSteps);
  uBs1.resize(potentialSteps);
  uBs2.resize(potentialSteps);
  us1s2.resize(potentialSteps);
  us1s1.resize(potentialSteps);
  us2s2.resize(potentialSteps);

  fHS.resize(potentialSteps);
  fBB.resize(potentialSteps);
  fBs1.resize(potentialSteps);
  fBs2.resize(potentialSteps);
  fs1s2.resize(potentialSteps);
  fs1s1.resize(potentialSteps);
  fs2s2.resize(potentialSteps);

  for (size_t i = 0; i < potentialSteps; ++i) {
    const double r = i * samplingStep;
    uBB[i] = (e_BB / e_min) * computeOmega(colloidRadius, colloidRadius, r);
    uBs1[i] = (e_Bs1 / e_min) * computeOmega(colloidRadius, radius_p1, r);
    uBs2[i] = (e_Bs2 / e_min) * computeOmega(colloidRadius, radius_p2, r);
    us1s2[i] = (e_s1s2 / e_min) * computeOmega(radius_p1, radius_p2, r);
    us2s2[i] = (e_s2s2 / e_min) * computeOmega(radius_p2, radius_p2, r);
    us1s1[i] = (e_s1s1 / e_min) * computeOmega(radius_p1, radius_p1, r);

    fBB[i] = (e_BB / e_min) *
             computeOmegaRadialDerivative(colloidRadius, colloidRadius, r);
    fBs1[i] = (e_Bs1 / e_min) *
              computeOmegaRadialDerivative(colloidRadius, radius_p1, r);
    fBs2[i] = (e_Bs2 / e_min) *
              computeOmegaRadialDerivative(colloidRadius, radius_p2, r);
    fs1s2[i] =
      (e_s1s2 / e_min) * computeOmegaRadialDerivative(radius_p1, radius_p2, r);
    fs2s2[i] =
      (e_s2s2 / e_min) * computeOmegaRadialDerivative(radius_p2, radius_p2, r);
    fs1s1[i] =
      (e_s1s1 / e_min) * computeOmegaRadialDerivative(radius_p1, radius_p1, r);

    if (r <= HSdiameter) {
      // setting up a Fake Hard Sphere Core
      double rm = pow(r, -fakeHSexponent);
      uHS[i] += fakeHScoefficient * ((rm - 2.) * rm + 1.);
      fHS[i] += -2. * fakeHSexponent * fakeHScoefficient * (rm - 1.) * rm;
    }
  }
}

void
PotentialForLammps::printLAMMPSpotentialsToFile(
  const std::string& outputDirName)
{
  // create output directory
  const std::string dirName = outputDirName + "/lammpspot";
  if (mkdir(dirName.c_str(), 0777) != 0)
    std::runtime_error("Problem while creating the directory.");

  // prepare strings for defining file names
  std::string interactionType[6];
  interactionType[0] = "BB";
  interactionType[1] = "Bs1";
  interactionType[2] = "Bs2";
  interactionType[3] = "s1s2";
  interactionType[4] = "s1s1";
  interactionType[5] = "s2s2";

  const size_t potentialSteps = uHS.size();

  for (int type = 0; type < 6; ++type) {
    // create the output file
    std::string fileName = dirName + "/" + interactionType[type] + ".table";
    std::ofstream potentialOutputFile(fileName);
    potentialOutputFile << "# potentials for lammps\n\n"
                        << interactionType[type] << "\nN " << potentialSteps - 1
                        << "\n\n";
    // -1 because we don't print r=0
    potentialOutputFile << std::scientific << std::setprecision(6);

    for (size_t i = 1; i < potentialSteps; ++i) {
      const double r = i * samplingStep;
      potentialOutputFile << i << '\t' << r << '\t';

      // compute force and potential depending on type and cutoff
      double printPotential, printForce;
      if (type == 0) {
        printPotential =
          (uHS[i] + uBB[i] > cutoffValue) ? cutoffValue : uHS[i] + uBB[i];
        printForce = ((fHS[i] + fBB[i]) * r < -cutoffValue)
                       ? -cutoffValue
                       : -(fHS[i] + fBB[i]);
      } else if (type == 1) {
        printPotential = (uBs1[i] > cutoffValue) ? cutoffValue : uBs1[i];
        printForce = (fBs1[i] * r < -cutoffValue) ? -cutoffValue : -fBs1[i];
      } else if (type == 2) {
        printPotential = (uBs2[i] > cutoffValue) ? cutoffValue : uBs2[i];
        printForce = (fBs2[i] * r < -cutoffValue) ? -cutoffValue : -fBs2[i];
      } else if (type == 3) {
        printPotential = (us1s2[i] > cutoffValue) ? cutoffValue : us1s2[i];
        printForce = (fs1s2[i] * r < -cutoffValue) ? -cutoffValue : -fs1s2[i];
      } else if (type == 4) {
        printPotential = (us1s1[i] > cutoffValue) ? cutoffValue : us1s1[i];
        printForce = (fs1s1[i] * r < -cutoffValue) ? -cutoffValue : -fs1s1[i];
      } else if (type == 5) {
        printPotential = (us2s2[i] > cutoffValue) ? cutoffValue : us2s2[i];
        printForce = (fs2s2[i] * r < -cutoffValue) ? -cutoffValue : -fs2s2[i];
      }
      // finally, you can print
      potentialOutputFile << printPotential << '\t' << printForce << '\n';
    }
    potentialOutputFile.close();
  }
}

size_t
PotentialForLammps::dist(double x, double y)
{
  double distance = std::sqrt(x * x + y * y);
  // std::cout << "dist = (" << x << ", " << y << ")\n";
  size_t distanceTab = size_t(distance / samplingStep);
  if (distanceTab > uHS.size()) {
    return 0;
  }
  return distanceTab;
}

static double ppot(const std::vector<double>& pot, size_t dist) {
  if (dist >= pot.size() || dist == 0) {
    return 0.;
  }
  return pot[dist];
}

void
PotentialForLammps::printRadialPotentialsToFile(
  std::string const& outputDirName)
{
  // create output directory
  const std::string dirName = outputDirName + "/radial_plots";
  if (mkdir(dirName.c_str(), 0777) != 0)
    std::runtime_error(
      "Problem while creating the directory for radial potentials.");

  std::vector<std::string> plotOrientations;
  if (symmetry == Symmetry::JANUS) {
    // E->quator, P->atch, B->ackside with no patch
    plotOrientations.push_back("EE");
    plotOrientations.push_back("EP");
    plotOrientations.push_back("EB");
    plotOrientations.push_back("PB");
    plotOrientations.push_back("PP");
    plotOrientations.push_back("BB");
  } else {
    plotOrientations.push_back("EE");
    plotOrientations.push_back("Ep1");
    plotOrientations.push_back("Ep2");
    plotOrientations.push_back("p1p2");
    plotOrientations.push_back("p1p1");
    plotOrientations.push_back("p2p2");
  }

  for (int type = 0; type < plotOrientations.size(); ++type) {
    bool saved = false;
    // create the output file
    std::string fileName = dirName + '/' + plotOrientations[type] + ".dat";
    std::ofstream potentialOutputFile(fileName);
    potentialOutputFile << std::scientific << std::setprecision(6);

    for (double r = HSdiameter; r < interactionRange; r += 0.001) {
      potentialOutputFile << r << '\t';
      size_t iCC, iCp1, iCp2, ip1p1, ip1C, ip1p2, ip2p1, ip2C, ip2p2;
      iCC = dist(r, 0);
      if (type == 0) {
        iCp1  = dist(r, eccentricity_p1);
        iCp2  = dist(r, eccentricity_p2);
        ip1p1 = dist(r, 0);
        ip1C  = dist(r, eccentricity_p1);
        ip1p2 = dist(r, eccentricity_p1 + eccentricity_p2);
        ip2p1 = dist(r, eccentricity_p1 + eccentricity_p2);
        ip2C  = dist(r, eccentricity_p2);
        ip2p2 = dist(r, 0);
      } else if (type == 1) {
        iCp1  = dist(r - eccentricity_p1, 0);
        iCp2  = dist(r + eccentricity_p2, 0);
        ip1p1 = dist(r - eccentricity_p1, eccentricity_p1);
        ip1C  = dist(r, eccentricity_p1);
        ip1p2 = dist(r + eccentricity_p2, eccentricity_p1);
        ip2p1 = dist(r - eccentricity_p1, -eccentricity_p2);
        ip2C  = dist(r, -eccentricity_p2);
        ip2p2 = dist(r + eccentricity_p2, -eccentricity_p2);
      } else if (type == 2) {
        iCp1  = dist(r + eccentricity_p1, 0);
        iCp2  = dist(r - eccentricity_p2, 0);
        ip1p1 = dist(r + eccentricity_p1, eccentricity_p1);
        ip1C  = dist(r, eccentricity_p1);
        ip1p2 = dist(r - eccentricity_p2, eccentricity_p1);
        ip2p1 = dist(r + eccentricity_p1, -eccentricity_p2);
        ip2C  = dist(r, -eccentricity_p2);
        ip2p2 = dist(r - eccentricity_p2, -eccentricity_p2);
      } else if (type == 3) {
        iCp1  = dist(r + eccentricity_p1, 0);
        iCp2  = dist(r - eccentricity_p2, 0);
        ip1p1 = dist(r, 0);
        ip1C  = dist(r - eccentricity_p1, 0);
        ip1p2 = dist(r - eccentricity_p1 - eccentricity_p2, 0);
        ip2p1 = dist(r + eccentricity_p1 + eccentricity_p2, 0);
        ip2C  = dist(r + eccentricity_p2, 0);
        ip2p2 = dist(r, 0);
      } else if (type == 4) {
        iCp1  = dist(r - eccentricity_p1, 0);
        iCp2  = dist(r + eccentricity_p2, 0);
        ip1p1 = dist(r -2*eccentricity_p1, 0);
        ip1C  = dist(r - eccentricity_p1, 0);
        ip1p2 = dist(r - eccentricity_p1 + eccentricity_p2, 0);
        ip2p1 = dist(r + eccentricity_p2 - eccentricity_p1, 0);
        ip2C  = dist(r + eccentricity_p2, 0);
        ip2p2 = dist(r + 2*eccentricity_p2, 0);
      } else if (type == 5) {
        iCp1  = dist(r + eccentricity_p1, 0);
        iCp2  = dist(r - eccentricity_p2, 0);
        ip1p1 = dist(r + 2*eccentricity_p1, 0);
        ip1C  = dist(r + eccentricity_p1, 0);
        ip1p2 = dist(r + eccentricity_p1 - eccentricity_p2, 0);
        ip2p1 = dist(r - eccentricity_p2 + eccentricity_p1, 0);
        ip2C  = dist(r - eccentricity_p2, 0);
        ip2p2 = dist(r - 2*eccentricity_p2, 0);
      }
      double printPotential =
        ppot(uHS, iCC) + ppot(uBB, iCC) + ppot(uBs1, iCp1) + ppot(uBs2, iCp2)
        + ppot(us1s1, ip1p1) + ppot(uBs1, ip1C) + ppot(us1s2, ip1p2)
        + ppot(us1s2, ip2p1) + ppot(uBs2, ip2C) + ppot(us2s2, ip2p2);
      // finally, you can print
      potentialOutputFile << printPotential << '\n';
      if (!saved) {
        saved = true;
        if (type == 0) {
          rvEE = printPotential;
        } else if (type == 1) {
          rvEP1 = printPotential;
        } else if (type == 2) {
          rvEP2 = printPotential;
        } else if (type == 3) {
          rvP1P1 = printPotential;
        } else if (type == 4) {
          rvP1P2 = printPotential;
        } else if (type == 5) {
          rvP2P2 = printPotential;
        }
      }
    }
    potentialOutputFile.close();
  }
}

static void
log(std::string orient, size_t dist, double pot)
{
#ifdef DEBUG_ORIENT
  std::cout << orient << " dist: " << dist * 1e-5 << ", pot: " << pot << "\n";
#endif
}

static std::string symmetryToString(Symmetry symmetry)
{
  switch(symmetry) {
    case Symmetry::JANUS:
      return "Janus";
    case Symmetry::SYMMETRIC:
      return "Symmetric";
    case Symmetry::ASYMMETRIC:
      return "Asymmetric";
    default:
      return "ERROR!";
  }
}

void
PotentialForLammps::printRecapFile(std::string const& outputDirName)
{
  // print out all parameters in a tidy tab
  std::string fileName = outputDirName + "/recap.dat";
  std::ofstream recapFile(fileName);

  recapFile
    << "INPUT VALUES GIVEN TO THE PROGRAM:\n"
    << "symmetry: " << symmetryToString(symmetry)
    << "\ndelta: " << delta
    << "\necc1: " << eccentricity_p1
    << "\nrad1: " << radius_p1
    << "\nvEE: " << vEE
    << "\nvEP1: " << vEP1
    << "\nvP1P1: " << vP1P1;
  if (symmetry == Symmetry::ASYMMETRIC) {
    recapFile
      << "\necc2: " << eccentricity_p2
      << "\nrad2: " << radius_p2
      << "\nvEP2: " << vEP2
      << "\nvP1P2: " << vP1P2
      << "\nvP2P2: " << vP2P2;
  }

  if (reducedMode) {
    std::cout << "Solved using the reduced solution\n";
  } else {
    std::cout << "Solved using the general solution\n";
  }

  // overlap volumes at contact
  //double fBB = computeOmega(colloidRadius, colloidRadius, HSdiameter);
  //double fBs1 =
  //  computeOmega(colloidRadius, radius_p1, HSdiameter - eccentricity_p1);
  //double fs1s1 =
  //  computeOmega(radius_p1, radius_p1, HSdiameter - 2. * eccentricity_p1);
  //double fBs2 =
  //  computeOmega(colloidRadius, radius_p2, HSdiameter - eccentricity_p2);
  //double fs1s2 = computeOmega(
  //  radius_p1, radius_p2, HSdiameter - eccentricity_p1 - eccentricity_p2);
  //double fs2s2 =
  //  computeOmega(radius_p2, radius_p2, HSdiameter - 2. * eccentricity_p2);

  recapFile << "\n\n\nOUTPUT VALUES:\n\n"
#ifdef DEBUG_MAIN
    << "overlap volumes:"
    << "\nfBB = " << fBB << "\nfBs1 = " << fBs1 << "\nfBs2 = " << fBs2
    << "\nfs1s1 = " << fs1s1 << "\nfs1s2 = " << fs1s2 << "\nfs2s2 = " << fs2s2
#endif
    << "EPSILON COEFFICIENTS:"
    << "\neps_BB   = " << e_BB << "\neps_Bs1  = " << e_Bs1
    << "\neps_Bs2  = " << e_Bs2 << "\neps_s1s1 = " << e_s1s1
    << "\neps_s1s2 = " << e_s1s2 << "\neps_s2s2 = " << e_s2s2
    << "\neps_min = " << e_min
    << "\n\nNORMALIZED EPSILON COEFFICIENTS:"
    << "\neps_BB   = " << e_BB / e_min << "\neps_Bs1  = " << e_Bs1 / e_min
    << "\neps_Bs2  = " << e_Bs2 / e_min << "\neps_s1s1 = " << e_s1s1 / e_min
    << "\neps_s1s2 = " << e_s1s2 / e_min << "\neps_s2s2 = " << e_s2s2 / e_min
    << "\n\nRESULTING CONTACT VALUES:";
   if (symmetry == Symmetry::JANUS) {
     recapFile << "\nback-back = " << rvP2P2
       << "\nback-patch = " << rvP1P2
       << "\npatch-patch = " << rvP1P1;
     //recapFile << "\nback-back = " << (e_BB * fBB) / e_min
     //  << "\nback-patch = " << (e_BB * fBB + fBs1 * e_Bs1) / e_min
     //  << "\npatch-patch = " << (e_BB * fBB + 2*fBs1 * e_Bs1 + fs1s1 * e_s1s1) / e_min;
   } else {
     recapFile << "\nvEE = " << rvEE
       << "\nvEP1 = " <<  rvEP1
       << "\nvEP2 = " <<  rvEP2
       << "\nvs1s1 = " << rvP1P1
       << "\nvs1s2 = " << rvP1P2
       << "\nvs2s2 = " << rvP2P2;
     //recapFile << "\nvEE = " << (e_BB * fBB) / e_min
     //  << "\nvEP1 = " << (e_BB * fBB + fBs1 * e_Bs1) / e_min
     //  << "\nvEP2 = " << (e_BB * fBB + fBs2 * e_Bs2) / e_min
     //  << "\nvs1s1 = " << (e_BB * fBB + fBs1 * e_Bs1 + fBs1 * e_Bs1 + fs1s1 * e_s1s1) / e_min
     //  << "\nvs1s2 = " << (e_BB * fBB + fBs1 * e_Bs1 + fBs2 * e_Bs2 + fs1s2 * e_s1s2) / e_min
     //  << "\nvs2s2 = " << (e_BB * fBB + fBs2 * e_Bs2 + fBs2 * e_Bs2 + fs2s2 * e_s2s2) / e_min;
   }

   recapFile << "\n\nPlease double check that the INPUT VALUES match the OUTPUT VALUES.";
}

void
PotentialForLammps::printAngularPotentialsToFile(
  std::string const& outputDirName)
{
  // create output directory
  const std::string dirName = outputDirName + "/angular_plots";
  if (mkdir(dirName.c_str(), 0777) != 0)
    std::runtime_error("Problem while creating the directory.");

  struct Orientation
  {
    std::string name;
    int theta_1;
    int theta_2;
  };
  std::vector<Orientation> plotOrientations;
  if (symmetry == Symmetry::JANUS) {
    plotOrientations.push_back({ "EE_EB_EE_EP_EE", 0, 180 });
    plotOrientations.push_back({ "PP_PE_PB_PE_PP", 180, 0 });
    plotOrientations.push_back({ "BB_BE_BP_BE_BB", 0, 0 });
  } else {
    plotOrientations.push_back({ "EE_EP2_EE_EP1_EE", 90, 90 });
    plotOrientations.push_back({ "P1P1_P1E_P1P2_P1E_P1P1", 180, 0 });
    plotOrientations.push_back({ "P2P2_P2E_P2P1_P2E_P2P2", 0, 180 });
  }

  for (auto& type : plotOrientations) {
    // create the output file
    std::string fileName = dirName + '/' + type.name + ".dat";
    std::ofstream potentialOutputFile(fileName);
    potentialOutputFile << std::scientific << std::setprecision(6);

    // set up the startingpoints
    double theta_1 = type.theta_1 * M_PI / 180.;
    double starting_theta_2 = type.theta_2;

    for (int angle = 0; angle < 360; angle += 5) {
      // compute potential depending on type and cutoff
      double theta_2 = (starting_theta_2 + angle) * M_PI / 180.;
      double potential = 0;
      double dx = 0.;
      double dy = 0.;
      size_t distTab = 0;

      // CC
      distTab = dist(HSdiameter, 0.0);
      if (distTab != 0) {
        double pot = uHS[distTab] + uBB[distTab];
        log("CC", distTab, pot);
        potential += pot;
      }
      // Cp1
      dx = HSdiameter - eccentricity_p1 * cos(theta_2);
      dy = eccentricity_p1 * sin(theta_2);
      distTab = dist(dx, dy);
      if (distTab != 0) {
        double pot = uBs1[distTab];
        log("Cp1", distTab, pot);
        potential += pot;
      }
      // Cp2
      dx = HSdiameter + eccentricity_p2 * cos(theta_2);
      dy = -eccentricity_p2 * sin(theta_2);
      distTab = dist(dx, dy);
      if (distTab != 0) {
        double pot = uBs2[distTab];
        log("Cp2", distTab, pot);
        potential += pot;
      }

      // p1C
      dx = eccentricity_p1 * cos(theta_1) + HSdiameter;
      dy = eccentricity_p1 * sin(theta_1);
      distTab = dist(dx, dy);
      if (distTab != 0) {
        double pot = uBs1[distTab];
        log("p1C", distTab, pot);
        potential += pot;
      }
      // p1p1
      dx = eccentricity_p1 * cos(theta_1) + HSdiameter -
           eccentricity_p1 * cos(theta_2);
      dy = eccentricity_p1 * sin(theta_1) - eccentricity_p1 * sin(theta_2);
      distTab = dist(dx, dy);
      if (distTab != 0) {
        double pot = us1s1[distTab];
        log("p1p1", distTab, pot);
        potential += pot;
      }
      // p1p2
      dx = eccentricity_p1 * cos(theta_1) + HSdiameter +
           eccentricity_p2 * cos(theta_2);
      dy = eccentricity_p1 * sin(theta_1) + eccentricity_p2 * sin(theta_2);
      distTab = dist(dx, dy);
      if (distTab != 0) {
        double pot = us1s2[distTab];
        log("p1p2", distTab, pot);
        potential += pot;
      }

      // p2C
      dx = HSdiameter - eccentricity_p2 * cos(theta_1);
      dy = eccentricity_p2 * sin(theta_1);
      distTab = dist(dx, dy);
      if (distTab != 0) {
        double pot = uBs2[distTab];
        log("p2C", distTab, pot);
        potential += pot;
      }
      // p2p1
      dx = HSdiameter - eccentricity_p2 * cos(theta_1) -
           eccentricity_p1 * cos(theta_2);
      dy = eccentricity_p2 * sin(theta_1) + eccentricity_p1 * sin(theta_2);
      distTab = dist(dx, dy);
      if (distTab != 0) {
        double pot = us1s2[distTab];
        log("p2p1", distTab, pot);
        potential += pot;
      }
      // p2p2
      dx = HSdiameter - eccentricity_p2 * cos(theta_1) +
           eccentricity_p2 * cos(theta_2);
      dy = eccentricity_p2 * sin(theta_1) - eccentricity_p2 * sin(theta_2);
      distTab = dist(dx, dy);
      if (distTab != 0) {
        double pot = us2s2[distTab];
        log("p2p2", distTab, pot);
        potential += pot;
      }

      // finally print it
      log("total", 0, potential);
      // std::cout << "\n\n\n";
      potentialOutputFile << angle << '\t' << potential << '\n';
    }
  }
}

void
PotentialForLammps::printPotentialAlongPathToFile(std::string const& outputDirName)
{
  // create output directory
  const std::string dirName = outputDirName + "/potential_on_path";
  if (mkdir(dirName.c_str(), 0777) != 0)
    std::runtime_error("Problem while creating the directory.");
  // create the output file
  std::string fileName = dirName + "/pot.dat";
  std::ofstream potentialPathOutputFile(fileName);
  potentialPathOutputFile << std::scientific << std::setprecision(6);

  /*
   * one IPC is kept constant in the origin, with p1 point up.
   * the second stays at fixed distance, and rotates.
   * four angles are necessary: phi, theta, alpha, beta
   * theta is the longitude of the CM of the second IPC relative to the first
   * phi is the latitude
   * alpha is the latitude of the first patch of the second IPC relative to its CM
   * beta is the longitude
   *
   * theta is 0 on the x axis and goes clockwise
   * phi is 0 on the z axis and grows when going down
   * alpha is like phi
   */

  /*
   * the path is defined as such:
   *        |    phi    |   theta   |   alpha  |   beta
   * start  |     90    |     0     |    90    |   180
   *   I    |     90    |     0     |  90 -> 0 |   180
   *   II   |  90 -> 0  |     0     |     0    |   180
   *  III   |     0     |     0     |  0 -> 90 |   180
   *   IV   |  0 -> 90  |     0     |     90   |   180
   *   V    |     90    |  0 -> 90  |     90   |   180
   *   VI   |     90    |    90     |  90 -> 0 |   180
  */

  double phi = 90;
  double theta = 0.;
  double alpha = 90.;
  double beta = 180.;
  for (alpha = 90.; alpha > 0.; alpha -= 5.) {
    potentialPathOutputFile << 90. - alpha << '\t'
        << computePotRot(phi, theta, alpha, beta) << '\n';
  }
  alpha = 0.;
  for (phi = 90.; phi > 0.; phi -= 5.) {
    potentialPathOutputFile << 90. + (90. - phi) << '\t'
        << computePotRot(phi, theta, alpha, beta) << '\n';
  }
  phi = 0.;
  for (alpha = 0.; alpha < 90.; alpha += 5.) {
    potentialPathOutputFile << 180. + alpha << '\t'
        << computePotRot(phi, theta, alpha, beta) << '\n';
  }
  alpha = 270.;
  for (phi = 0.; phi < 90.; phi += 5.) {
    potentialPathOutputFile << 270. + phi << '\t'
        << computePotRot(phi, theta, alpha, beta) << '\n';
  }
  phi = 90.;
  for (theta = 0.; theta < 90.; theta += 5) {
    potentialPathOutputFile << 360. + theta << '\t'
        << computePotRot(phi, theta, alpha, beta) << '\n';
  }
  theta = 90.;
  for (alpha = 90.; alpha > 0.; alpha -= 5.) {
    potentialPathOutputFile << 450. + (90. - alpha) << '\t'
        << computePotRot(phi, theta, alpha, beta) << '\n';
  }
  alpha = 360.;
}

size_t
PotentialForLammps::dist(const double* xa, const double* xb)
{
  double distance = std::sqrt( std::pow(xa[0] - xb[0], 2)
    + std::pow(xa[1] - xb[1], 2) + std::pow(xa[2] - xb[2], 2) );
  size_t distanceTab = size_t(distance / samplingStep);
  if (distanceTab > uHS.size()) {
    return 0;
  }
  return distanceTab;
}

double PotentialForLammps::computePotRot(double phi, double theta, double alpha, double beta) {
  /*
   * one IPC is kept constant in the origin, with p1 point up.
   * the second stays at fixed distance, and rotates.
   * four angles are necessary: phi, theta, alpha, beta
   * theta is the longitude of the CM of the second IPC relative to the first
   * phi is the latitude
   * alpha is the latitude of the first patch of the second IPC relative to its CM
   * beta would be the longitude, but we don't care about it and we keep it constant
   *
   * theta is 0 on the x axis and goes clockwise
   * phi is 0 on the z axis and grows when going down
   * alpha is like phi
   * beta is, at the moment, ignored.
   */
  double potential = 0;
  size_t distTab = 0;
  phi   *=  M_PI / 180.;
  theta *=  M_PI / 180.;
  alpha *=  M_PI / 180.;
  beta  *=  M_PI / 180.;

  // ipc 1
  double ipc1cb[3] = {0., 0., 0.};
  double ipc1p1[3] = {0., 0.,  eccentricity_p1};
  double ipc1p2[3] = {0., 0., -eccentricity_p2};
  // ipc 2
  double ipc2cb[3] = { std::cos(theta)*std::sin(phi),  std::sin(theta)*std::sin(phi),  std::cos(phi) };
  double ipc2temp[3] = { std::cos(beta)*std::sin(alpha), std::sin(beta)*std::sin(alpha), std::cos(alpha) };
  double ipc2p1[3] = { 0., 0., 0. };
  double ipc2p2[3] = { 0., 0., 0. };
  for (int i = 0; i < 3; ++i) {
    ipc2cb[i] *= HSdiameter;
    ipc2p1[i] = ipc2cb[i] + eccentricity_p1 * ipc2temp[i];
    ipc2p2[i] = ipc2cb[i] - eccentricity_p2 * ipc2temp[i];
  }

  // CC
  distTab = dist(ipc1cb, ipc2cb);
  if (distTab != 0) {
    potential += uHS[distTab] + uBB[distTab];
  }
  // Cp1
  distTab = dist(ipc1cb, ipc2p1);
  if (distTab != 0) {
    potential += uBs1[distTab];
  }
  // Cp2
  distTab = dist(ipc1cb, ipc2p2);
  if (distTab != 0) {
    potential += uBs2[distTab];
  }
  // p1C
  distTab = dist(ipc1p1, ipc2cb);
  if (distTab != 0) {
    potential += uBs1[distTab];
  }
  // p1p1
  distTab = dist(ipc1p1, ipc2p1);
  if (distTab != 0) {
    potential += us1s1[distTab];
  }
  // p1p2
  distTab = dist(ipc1p1, ipc2p2);
  if (distTab != 0) {
    potential += us1s2[distTab];
  }
  // p2C
  distTab = dist(ipc1p2, ipc2cb);
  if (distTab != 0) {
    potential += uBs2[distTab];
  }
  // p2p1
  distTab = dist(ipc1p2, ipc2p1);
  if (distTab != 0) {
    potential += us1s2[distTab];
  }
  // p2p2
  distTab = dist(ipc1p2, ipc2p2);
  if (distTab != 0) {
    potential += us2s2[distTab];
  }
  return potential;
}
