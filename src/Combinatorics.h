#ifndef INCLUDE_Combinatorics_h
#define INCLUDE_Combinatorics_h


#include "Temp_Matrix.h"
#include "Face_bitset.h"


template<typename T>
std::vector<T> VectorAsSet(std::vector<T> const& V)
{
  std::set<T> eSet;
  for (auto & eVal : V)
    eSet.insert(eVal);
  std::vector<T> eV;
  for (auto & eVal : eSet)
    eV.push_back(eVal);
  return eV;
}


#endif
