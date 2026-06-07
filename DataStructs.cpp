
#include "DataStructs.h"
#include <iostream>

template<class T>
void DataStruct<T>::setSize(int _size)
{
  if(initialized)
  {
    // check if it's the same size
    if(_size != size)
    {
      // this case isn't handled. For the time being
      // simple abort program!
      std::runtime_error("Trying to allocate an already allocated DataStruct. Aborting...");
    }
  }else
  {
    data = new T[_size];
    size = _size;
    initialized = true;
  }
};

template class DataStruct<float>;
template class DataStruct<double>;