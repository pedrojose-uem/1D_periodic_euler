#ifndef _RUNGE_KUTTA
#define _RUNGE_KUTTA

/*
RUNGE KUTTA 4 CLASS

This class allows time stepping. For flexibility each RK step is taken as follows:
  initRK();
  for(s = 0; s<rk.getNumSteps(); s++)
    1. stepUi(dt)
    2. impose BCs over the current Ui which can be obtained through currentU()
    3. compute current F (i.e., Fi = F(Ui) )
    5. setFi(Fi)
  finalizeRK(dt);

*/

#include "DataStructs.h"
#include "FluxFunctions.h"

template<class T>
class RungeKutta4 
{
  private:
    int nSteps;
    int currentStep;

    T *coeffsA, *coeffsB;

    // reference to solution (Un)
    DataStruct<T> &Un;

    // intermediate solution
    DataStruct<T> Ui;

    // RHS 
    DataStruct<T> *fi;

  public:

    // default constructor
    RungeKutta4(DataStruct<T> &_Un);

    // default destructor
    ~RungeKutta4();

    int getNumSteps();

    // initialize the RK
    void initRK();

    // finalizes the RK (updates Un)
    void finalizeRK(const T dt);

    /*
    For the step to work properly, the user must provide the appropriate F for the current Ui.
    This is done this way becase the user might want to modify the Ui or Fi so that 
    Boundary conditions can be imposed
    */
   void stepUi(T dt);
   void setFi(DataStruct<T> &_F);

   // current Ui
   DataStruct<T>* currentU();
};

#endif // _RUNGE_KUTTA 