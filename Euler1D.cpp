#include "Euler1D.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>

template<class T>
int EulerCentral1D<T>::idRho(int j) const
{
  return (layout == EulerLayout::AoS) ? 3*j : j;
}

template<class T>
int EulerCentral1D<T>::idRhoU(int j) const
{
  return (layout == EulerLayout::AoS) ? 3*j + 1 : nPoints + j;
}

template<class T>
int EulerCentral1D<T>::idRhoE(int j) const
{
  return (layout == EulerLayout::AoS) ? 3*j + 2 : 2*nPoints + j;
}

template<class T>
EulerCentral1D<T>::EulerCentral1D(DataStruct<T> &_mesh, int _nPoints, T _gamma, EulerLayout _layout):
mesh(_mesh), nPoints(_nPoints), gamma(_gamma), layout(_layout)
{
  RHS.setSize(3*nPoints);
}

template<class T>
void EulerCentral1D<T>::fluxAt(const T *U, int j, T &f0, T &f1, T &f2) const
{
  const T rho  = U[idRho(j)];
  const T rhou = U[idRhoU(j)];
  const T rhoE = U[idRhoE(j)];

  const T invRho = 1.0 / rho;
  const T u = rhou * invRho;
  const T p = (gamma - 1.0) * (rhoE - 0.5 * rhou * rhou * invRho);

  f0 = rhou;
  f1 = rhou * u + p;
  f2 = u * (rhoE + p);
}

template<class T>
void EulerCentral1D<T>::eval(DataStruct<T> &Uin)
{
  T *dataRHS = RHS.getData();
  const T *dataU = Uin.getData();
  const T *dataMesh = mesh.getData();

  for(int j = 0; j < nPoints; j++)
  {
    const int jm = (j == 0) ? nPoints - 2 : j - 1;
    const int jp = (j == nPoints - 1) ? 1 : j + 1;

    T fm0, fm1, fm2;
    T fp0, fp1, fp2;
    fluxAt(dataU, jm, fm0, fm1, fm2);
    fluxAt(dataU, jp, fp0, fp1, fp2);

    T dx;
    if(j == 0)
    {
      dx = dataMesh[nPoints-1] - dataMesh[nPoints-2];
      dx += dataMesh[1] - dataMesh[0];
    }
    else if(j == nPoints - 1)
    {
      dx = dataMesh[nPoints-1] - dataMesh[nPoints-2];
      dx += dataMesh[1] - dataMesh[0];
    }
    else
    {
      dx = dataMesh[j+1] - dataMesh[j-1];
    }

    dataRHS[idRho(j)]  = -(fp0 - fm0) / dx;
    dataRHS[idRhoU(j)] = -(fp1 - fm1) / dx;
    dataRHS[idRhoE(j)] = -(fp2 - fm2) / dx;
  }

  // Enforce periodic equality between first and last node.
  dataRHS[idRho(nPoints-1)]  = dataRHS[idRho(0)];
  dataRHS[idRhoU(nPoints-1)] = dataRHS[idRhoU(0)];
  dataRHS[idRhoE(nPoints-1)] = dataRHS[idRhoE(0)];
}

template<class T>
DataStruct<T>& EulerCentral1D<T>::ref2RHS()
{
  return RHS;
}

template<class T>
static int rhoId(int j, int nPoints, EulerLayout layout)
{
  return (layout == EulerLayout::AoS) ? 3*j : j;
}

template<class T>
static int rhoUId(int j, int nPoints, EulerLayout layout)
{
  return (layout == EulerLayout::AoS) ? 3*j + 1 : nPoints + j;
}

template<class T>
static int rhoEId(int j, int nPoints, EulerLayout layout)
{
  return (layout == EulerLayout::AoS) ? 3*j + 2 : 2*nPoints + j;
}

template<class T>
void initializeSodLikePeriodic(DataStruct<T> &U, DataStruct<T> &x, int nPoints, EulerLayout layout, T gamma)
{
  T *dataU = U.getData();
  T *dataX = x.getData();

  for(int j = 0; j < nPoints; j++)
  {
    dataX[j] = T(j) / T(nPoints - 1);

    // Smooth periodic perturbation around a uniform state. It avoids the
    // discontinuity of a shock tube and keeps the central scheme stable enough
    // for this educational profiling exercise.
    const T rho = 1.0 + 0.2 * std::sin(2.0 * M_PI * dataX[j]);
    const T u = 0.3;
    const T p = 1.0 + 0.1 * std::sin(2.0 * M_PI * dataX[j]);

    dataU[rhoId<T>(j,nPoints,layout)]  = rho;
    dataU[rhoUId<T>(j,nPoints,layout)] = rho * u;
    dataU[rhoEId<T>(j,nPoints,layout)] = p / (gamma - 1.0) + 0.5 * rho * u * u;
  }

  // Periodic endpoint.
  dataU[rhoId<T>(nPoints-1,nPoints,layout)]  = dataU[rhoId<T>(0,nPoints,layout)];
  dataU[rhoUId<T>(nPoints-1,nPoints,layout)] = dataU[rhoUId<T>(0,nPoints,layout)];
  dataU[rhoEId<T>(nPoints-1,nPoints,layout)] = dataU[rhoEId<T>(0,nPoints,layout)];
}

template<class T>
void writeEulerToFile(DataStruct<T> &x, DataStruct<T> &U, int nPoints, EulerLayout layout, std::string name)
{
  std::ofstream file(name, std::ios_base::trunc);
  if(!file.is_open())
  {
    std::cout << "Couldn't open file " << name << std::endl;
    return;
  }

  file << "x,rho,rhou,rhoE,u,p\n";
  const T gamma = 1.4;
  for(int j = 0; j < nPoints; j++)
  {
    const T rho = U.getData()[rhoId<T>(j,nPoints,layout)];
    const T rhou = U.getData()[rhoUId<T>(j,nPoints,layout)];
    const T rhoE = U.getData()[rhoEId<T>(j,nPoints,layout)];
    const T u = rhou / rho;
    const T p = (gamma - 1.0) * (rhoE - 0.5 * rhou * rhou / rho);
    file << x.getData()[j] << "," << rho << "," << rhou << "," << rhoE << "," << u << "," << p << "\n";
  }
}

template<class T>
T calcEulerVariation(DataStruct<T> &U, DataStruct<T> &Uinit)
{
  T err = 0.;
  const T *dataU = U.getData();
  const T *dataInit = Uinit.getData();

  for(int n = 0; n < U.getSize(); n++)
  {
    const T diff = dataU[n] - dataInit[n];
    err += diff * diff;
  }

  return std::sqrt(err);
}

template class EulerCentral1D<float>;
template class EulerCentral1D<double>;

template void initializeSodLikePeriodic<float>(DataStruct<float>&, DataStruct<float>&, int, EulerLayout, float);
template void initializeSodLikePeriodic<double>(DataStruct<double>&, DataStruct<double>&, int, EulerLayout, double);

template void writeEulerToFile<float>(DataStruct<float>&, DataStruct<float>&, int, EulerLayout, std::string);
template void writeEulerToFile<double>(DataStruct<double>&, DataStruct<double>&, int, EulerLayout, std::string);

template float calcEulerVariation<float>(DataStruct<float>&, DataStruct<float>&);
template double calcEulerVariation<double>(DataStruct<double>&, DataStruct<double>&);
