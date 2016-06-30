#ifndef INCLUDE_Temp_common_h
#define INCLUDE_Temp_common_h




#include <ctype.h>
#include <malloc.h>
#include <unistd.h>
#include <getopt.h>
#include <chrono>
#include <ctime>


#include <math.h>


#include <string>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <stdio.h>
#include <stdlib.h>


#include <exception>
#include <vector>
#include <list>
#include <set>
#include <map>


#include <functional>
#include <algorithm>


#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>


typedef unsigned long ulong;
typedef unsigned int uint;


typedef std::vector<std::vector<int> > VectVectInt;


struct TerminalException {
  int eVal;
};


template<typename T>
T VectorMin(std::vector<T> const& eVect)
{
  T eMin=eVect[0];
  for (T eVal : eVect)
    if (eVal < eMin)
      eMin=eVal;
  return eMin;
}


template<typename T>
T VectorMax(std::vector<T> const& eVect)
{
  T eMax=eVect[0];
  for (T eVal : eVect)
    if (eVal > eMax)
      eMax=eVal;
  return eMax;
}


template<typename T>
T T_abs(T const& eVal)
{
  if (eVal > 0)
    return eVal;
  T fVal= - eVal;
  return fVal;
}


template<typename T>
void WriteStdVectorGAP(std::ostream& os, std::vector<T> const& V)
{
  os << "[";
  bool IsFirst=true;
  for (auto & eVal : V) {
    if (IsFirst == false)
      os << ",";
    IsFirst=false;
    os << eVal;
  }
  os << "]";
}


template<typename T>
std::istream& operator>>(std::istream& is, std::vector<T>& obj)
{
  int n;
  is >> n;
  obj.resize(n);
  for (int i=0; i<n; i++) {
    T eVal;
    is >> eVal;
    obj[i]=eVal;
  }
  return is;
}
template<typename T>
std::ostream& operator<<(std::ostream& os, std::vector<T> const& ListVal)
{
  int n=ListVal.size();
  os << " " << n;
  for (int i=0; i<n; i++)
    os << " " << ListVal[i];
  return os;
}


template<typename T>
struct CollectedResult {
  std::vector<T> LVal;
  std::vector<int> LMult;
};


template<typename T>
CollectedResult<T> Collected(std::vector<T> const& eVect)
{
  std::set<T> SetVal;
  for (auto & eVal : eVect)
    SetVal.insert(eVal);
  std::vector<T> LVal;
  for (auto & eVal : SetVal)
    LVal.push_back(eVal);
  int eSize=LVal.size();
  std::vector<int> LMult(eSize,0);
  auto UpPosition=[&](T const& eVal) -> void {
    for (int i=0; i<eSize; i++)
      if (LVal[i] == eVal) {
	LMult[i] += 1;
	return;
      }
    std::cerr << "Should never reach that stage\n";
    throw TerminalException{1};
  };
  for (auto & eVal : eVect)
    UpPosition(eVal);
  return {LVal, LMult};
}


int NextIdx(int const& len,int const& i)
{
  if (i == len-1)
    return 0;
  return i+1;
}


int PrevIdx(int const& len,int const& i)
{
  if (i == 0)
    return len-1;
  return i-1;
}


template<typename T>
int PositionVect(std::vector<T> const& V, T const& eVal)
{
  int len=V.size();
  for (int i=0; i<len; i++)
    if (V[i] == eVal)
      return i;
  return -1;
}


#endif
