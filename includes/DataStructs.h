#ifndef _DATA_STRUCTS
#define _DATA_STRUCTS

#include <cstddef>

template<class T>
class DataStruct 
{
  private:
    int size;
    T *data;
    bool initialized;

  public:

    // constructors & destructor
    inline DataStruct()
    {
      size = 0;
      data = NULL;
      initialized = false;
    };

    inline DataStruct(int _size)
    {
      data = new T[_size];
      size = _size;
      initialized = true;
    };

    inline ~DataStruct() 
    {
      if(initialized) delete[] data;
    };

    // accessors
    inline int getSize()
    {
      return size;
    }

    inline T* getData()
    {
      return data;
    }

    inline T getData(int i)
    {
      return data[i];
    }

    void setSize(int _size);

    inline DataStruct<T>& operator=(DataStruct<T> &rhs)
    {
      this->setSize(rhs.getSize());
      T *dataRHS = rhs.getData();

      // copy data
      for(int n = 0; n < size; n++)
      {
        data[n] = dataRHS[n];
      }

      return *this;
    }

};

#endif // _DATA_STRUCTS