#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>

enum class IpcType {
  NONE,
  JANUS,
  IPC,
  ASYM_IPC,
};

class PotentialForLammps {
 public:
  PotentialForLammps(std::string const &inputFileName, const IpcType ipcType,
                     bool startFromContactValues);
  void printLAMMPSpotentialsToFile(std::string const &outputDirName);
  void printRadialPotentialsToFile(std::string const &outputDirName);
  void printAngularPotentialsToFile(std::string const &outputDirName);

 private:
  IpcType ipcType;
  double e_BB, e_Bs1, e_Bs2, e_s1s1, e_s1s2, e_s2s2, e_min;
  double vEE, vEP1, vEp2, vP1P1, vP1P2, vP2P2;

  double delta, ipcRadius, interactionRange;
  double eccentricity_p1, radius_p1;
  double eccentricity_p2, radius_p2;
  double fakeHSdiameter, fakeHScoefficient, fakeHSexponent;

  double samplingStep, cutoffValue;

  std::vector<double> uHS, uBB, uBs1, uBs2, us1s2, us1s1, us2s2;
  std::vector<double> fHS, fBB, fBs1, fBs2, fs1s2, fs1s1, fs2s2;

  void computeSiteSitePotentials();
  void initFromEpsilons(std::string const &inputFileName);
  void readContactValues(std::string const &inputFileName);
  void computeEpsilonsFromContactValues();
  size_t dist(double x, double y);
};
