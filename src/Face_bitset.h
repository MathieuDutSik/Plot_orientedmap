#ifndef INCLUDE_Face_bitset_h
#define INCLUDE_Face_bitset_h


#include "Temp_Matrix.h"


#include <bitset>
#include <boost/dynamic_bitset.hpp>


typedef boost::dynamic_bitset<> Face;


bool operator<(Face const& x, Face const& y)
{
  int len=x.size();
  for (int i=0; i<len; i++) {
    if (x[i] == 0 && y[i] == 1)
      return true;
    if (x[i] == 1 && y[i] == 0)
      return false;
  }
  return false;
}


#endif
