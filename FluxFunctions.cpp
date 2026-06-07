#include "FluxFunctions.h"

template<class T>
FluxFunction<T>::FluxFunction()
{
  
};

template<class T>
LinearFlux<T>::LinearFlux()
{
  c = 1.;
};

template<class T>
void LinearFlux<T>::computeFlux(DataStruct<T> &U, DataStruct<T> &F)
{
  T *dataU = U.getData();
  T *dataF = F.getData();

  for(int n = 0; n < U.getSize(); n++)
  {
    dataF[n] = c * dataU[n];
  };
};

template<class T>
T LinearFlux<T>::computeFlux(const T &Ui)
{
  return c * Ui;
};

template class FluxFunction<float>;
template class FluxFunction<double>;

template class LinearFlux<float>;
template class LinearFlux<double>;