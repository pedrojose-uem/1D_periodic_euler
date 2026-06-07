#ifndef _EULER_1D
#define _EULER_1D

#include <string>
#include "DataStructs.h"

enum class EulerLayout
{
  AoS, // [rho_1, rhou_1, rhoE_1, rho_2, rhou_2, rhoE_2, ...]
  SoA  // [rho_1..rho_N, rhou_1..rhou_N, rhoE_1..rhoE_N]
};

template<class T>
class EulerCentral1D
{
private:
  DataStruct<T> RHS;
  DataStruct<T> &mesh;
  int nPoints;
  T gamma;
  EulerLayout layout;

  inline int idRho(int j) const;
  inline int idRhoU(int j) const;
  inline int idRhoE(int j) const;

  inline void fluxAt(const T *U, int j, T &f0, T &f1, T &f2) const;

public:
  EulerCentral1D(DataStruct<T> &_mesh, int _nPoints, T _gamma, EulerLayout _layout);

  void eval(DataStruct<T> &Uin);
  DataStruct<T>& ref2RHS();
};

template<class T>
void initializeSodLikePeriodic(DataStruct<T> &U, DataStruct<T> &x, int nPoints, EulerLayout layout, T gamma);

template<class T>
void writeEulerToFile(DataStruct<T> &x, DataStruct<T> &U, int nPoints, EulerLayout layout, std::string name);

template<class T>
T calcEulerVariation(DataStruct<T> &U, DataStruct<T> &Uinit);

#endif // _EULER_1D
