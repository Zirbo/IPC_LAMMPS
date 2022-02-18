#include "printPotential.hpp"

#include <sys/stat.h>

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

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
  if (symmetry == Symmetry::ASYMMETRIC) {
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
  std::cout << "interaction range = 2 x max(" << colloidRadius << ", "
            << range_p1 << ", " << range_p2 << ") = " << interactionRange
            << "\n\n";

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
    std::cout << "WARNING: fBB is zero. That means delta is zero.\n"
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

  std::cout
    << "overlap volumes:\n"
    << "\nfBB = " << fBB << "\nfBs1 = " << fBs1 << "\nfBs2 = " << fBs2
    << "\nfs1s1 = " << fs1s1 << "\nfs1s2 = " << fs1s2 << "\nfs2s2 = " << fs2s2
    << "\n\nNORMALIZED COEFFICIENTS:"
    << "\neps_BB   = " << e_BB / e_min << "\neps_Bs1  = " << e_Bs1 / e_min
    << "\neps_Bs2  = " << e_Bs2 / e_min << "\neps_s1s1 = " << e_s1s1 / e_min
    << "\neps_s1s2 = " << e_s1s2 / e_min << "\neps_s2s2 = " << e_s2s2 / e_min
    << "\n\nResulting contact values:\nvEE = " << e_BB * fBB
    << "\nvEP1 = " << e_BB * fBB + fBs1 * e_Bs1
    << "\nvEP2 = " << e_BB * fBB + fBs2 * e_Bs2
    << "\nvs1s2 = " << e_BB * fBB + fBs1 * e_Bs1 + fBs2 * e_Bs2 + fs1s1 * e_s1s1
    << "\nvs1s1 = " << e_BB * fBB + fBs1 * e_Bs1 + fBs1 * e_Bs1 + fs1s2 * e_s1s2
    << "\nvs2s2 = " << e_BB * fBB + fBs2 * e_Bs2 + fBs2 * e_Bs2 + fs2s2 * e_s2s2
    << "\nDo they match the ones you inserted?\n";
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

void
PotentialForLammps::printRadialPotentialsToFile(
  std::string const& outputDirName)
{
  // create output directory
  const std::string dirName = outputDirName + "/lammpspot_radial_plots";
  if (mkdir(dirName.c_str(), 0777) != 0)
    std::runtime_error(
      "Problem while creating the directory for radial potentials.");

  std::vector<std::string> plotOrientations;
  if (symmetry == Symmetry::JANUS) {
    // E->quator, p->atch, C->backside with no patch
    plotOrientations.push_back("EE");
    plotOrientations.push_back("EP");
    plotOrientations.push_back("EC");
    plotOrientations.push_back("PC");
    plotOrientations.push_back("PP");
    plotOrientations.push_back("CC");
  } else {
    plotOrientations.push_back("EE");
    plotOrientations.push_back("Ep1");
    plotOrientations.push_back("Ep2");
    plotOrientations.push_back("p1p2");
    plotOrientations.push_back("p1p1");
    plotOrientations.push_back("p2p2");
  }

  for (int type = 0; type < plotOrientations.size(); ++type) {
    // create the output file
    std::string fileName = dirName + '/' + plotOrientations[type] + ".dat";
    std::ofstream potentialOutputFile(fileName);
    potentialOutputFile << std::scientific << std::setprecision(6);

    for (double r = HSdiameter; r < interactionRange; r += samplingStep) {
      potentialOutputFile << r << '\t';

      size_t iBB = size_t(r / samplingStep);
      size_t iBs1 = size_t((r - eccentricity_p1) / samplingStep);
      size_t iBs2 = size_t((r - eccentricity_p2) / samplingStep);
      // compute potential depending on type and cutoff
      double printPotential;
      if (type == 0) {
        printPotential = uHS[iBB] + uBB[iBB];
      } else if (type == 1) {
        printPotential = uHS[iBB] + uBB[iBB] + uBs1[iBs1];
      } else if (type == 2) {
        printPotential = uHS[iBB] + uBB[iBB] + uBs2[iBs2];
      } else if (type == 3) {
        size_t is1s2 =
          size_t((r - eccentricity_p1 - eccentricity_p2) / samplingStep);
        printPotential =
          uHS[iBB] + uBB[iBB] + uBs1[iBs1] + uBs2[iBs2] + us1s2[is1s2];
      } else if (type == 4) {
        size_t is1s1 = size_t((r - 2 * eccentricity_p1) / samplingStep);
        printPotential =
          uHS[iBB] + uBB[iBB] + uBs1[iBs1] + uBs1[iBs1] + us1s1[is1s1];
      } else if (type == 5) {
        size_t is2s2 = size_t((r - 2 * eccentricity_p2) / samplingStep);
        printPotential =
          uHS[iBB] + uBB[iBB] + uBs2[iBs2] + uBs2[iBs2] + us2s2[is2s2];
      }
      // finally, you can print
      potentialOutputFile << printPotential << '\n';
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

static void
log(std::string orient, size_t dist, double pot)
{
// #define DEBUG
#ifdef DEBUG
  std::cout << orient << " dist: " << dist * 1e-5 << ", pot: " << pot << "\n";
#endif
}

void
PotentialForLammps::printAngularPotentialsToFile(
  std::string const& outputDirName)
{
  // create output directory
  const std::string dirName = outputDirName + "/lammpspot_angular_plots";
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
    plotOrientations.push_back({ "CC", 0, 180 });
    plotOrientations.push_back({ "PP", 180, 0 });
    plotOrientations.push_back({ "PC", 0, 0 });
  } else {
    plotOrientations.push_back({ "EE", 90, 90 });
    plotOrientations.push_back({ "P1P1", 180, 0 });
    plotOrientations.push_back({ "P2P2", 0, 180 });
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
