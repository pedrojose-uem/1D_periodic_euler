#include <chrono>
#include <iomanip>
#include <iostream>
#include <string>

#include "DataStructs.h"
#include "Euler1D.h"
#include "rk4.h"

#ifdef _DOUBLE_
#define FLOATTYPE double
#else
#define FLOATTYPE float
#endif

int main(int narg, char **argv)
{
  int numPoints = 100000;
  std::string layoutArg = "soa";

  if(narg < 2 || narg > 3)
  {
    std::cout << "Usage: ./euler.p NumPoints [soa|aos]" << std::endl;
    return 1;
  }

  numPoints = std::stoi(argv[1]);
  if(narg == 3) layoutArg = argv[2];

  EulerLayout layout;
  if(layoutArg == "aos")
  {
    layout = EulerLayout::AoS;
  }
  else if(layoutArg == "soa")
  {
    layout = EulerLayout::SoA;
  }
  else
  {
    std::cout << "Unknown layout. Use 'soa' or 'aos'." << std::endl;
    return 1;
  }

  const FLOATTYPE gamma = 1.4;

  DataStruct<FLOATTYPE> U(3*numPoints), xj(numPoints);
  initializeSodLikePeriodic(U, xj, numPoints, layout, gamma);

  DataStruct<FLOATTYPE> Uinit;
  Uinit = U;

  RungeKutta4<FLOATTYPE> rk(U);
  EulerCentral1D<FLOATTYPE> rhs(xj, numPoints, gamma, layout);

  const FLOATTYPE dx = xj.getData()[1] - xj.getData()[0];
  const FLOATTYPE CFL = 0.05;
  FLOATTYPE dt = CFL * dx;
  const FLOATTYPE t_final = 0.05;
  FLOATTYPE time = 0.;

  DataStruct<FLOATTYPE> Ui(U.getSize());

  writeEulerToFile(xj, U, numPoints, layout, "euler_initial_" + layoutArg + ".csv");

  const auto start = std::chrono::high_resolution_clock::now();

  while(time < t_final)
  {
    if(time + dt >= t_final) dt = t_final - time;

    rk.initRK();
    for(int s = 0; s < rk.getNumSteps(); s++)
    {
      rk.stepUi(dt);
      Ui = *rk.currentU();
      rhs.eval(Ui);
      rk.setFi(rhs.ref2RHS());
    }
    rk.finalizeRK(dt);
    time += dt;
  }

  const auto end = std::chrono::high_resolution_clock::now();
  const std::chrono::duration<double> elapsed = end - start;

  writeEulerToFile(xj, U, numPoints, layout, "euler_final_" + layoutArg + ".csv");

  const FLOATTYPE variation = calcEulerVariation(U, Uinit);

  std::cout << std::setprecision(6)
            << "Layout: " << layoutArg
            << " NumPoints: " << numPoints
            << " Comp. time: " << elapsed.count() << " sec"
            << " Variation L2: " << variation
            << std::endl;

  return 0;
}
