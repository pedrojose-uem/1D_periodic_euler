#ifndef _FLUX_FUNCTIONS
#define _FLUX_FUNCTIONS

#include "DataStructs.h"

// base class for flux functions
template<class T>
class FluxFunction
{
  public:
    FluxFunction();

    // in U; out F
    virtual void computeFlux(DataStruct<T> &U, DataStruct<T> &F) = 0;

    // same as above but at the node level
    virtual T computeFlux(const T &Ui) = 0;
};

template<class T>
class LinearFlux : public FluxFunction<T>
{
  private:
    T c;

  public:
    LinearFlux();

    // in U; out F
    virtual void computeFlux(DataStruct<T> &U, DataStruct<T> &F);

    // same as above but at the node level
    virtual T computeFlux(const T &Ui);
};

#endif // _FLUX_FUNCTIONS