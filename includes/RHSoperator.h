#ifndef _RHS_OPERATOR
#define _RHS_OPERATOR

#include "DataStructs.h"
#include "FluxFunctions.h"

template<class T>
class RHSOperator
{
public:
  RHSOperator();
  ~RHSOperator();

  virtual void eval() = 0;
  virtual void eval(DataStruct<T> &Uin) = 0;

  virtual DataStruct<T>& ref2RHS() = 0;
};

template<class T>
class Central1D : public RHSOperator<T>
{
private:

  // structure containing the RHS values
  DataStruct<T> RHS;
  
  // reference to current solution
  DataStruct<T> &U;

  // reference to mesh 
  // TODO: change to a mesh structure
  DataStruct<T> &mesh;

  // reference to flux function
  FluxFunction<T> &F;

  void evalRHS(DataStruct<T> &Uin);

public:
  Central1D(DataStruct<T> &_U, DataStruct<T> &_mesh, FluxFunction<T> &_F);
  ~Central1D();

  virtual void eval();
  virtual void eval(DataStruct<T> &Uin);

  virtual DataStruct<T>& ref2RHS();

};

#endif // _RHS_OPERATOR