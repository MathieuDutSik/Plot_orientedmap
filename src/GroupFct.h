#ifndef INCLUDE_GroupFct_h
#define INCLUDE_GroupFct_h

#include "Temp_common.h"
#include "Basic_file.h"
#include "Basic_string.h"
#include "Boost_bitset.h"

template<typename Telt>
std::vector<int> PermutationOrbit(Telt const& ePerm)
{
  int siz=ePerm.size();
  std::vector<int> StatusOrbit(siz,-1);
  int idxOrbit=0;
  auto GetUnsetPoint=[&](void) -> int {
    for (int i=0; i<siz; i++)
      if (StatusOrbit[i] == -1)
	return i;
    return -1;
  };
  while(true) {
    int iPoint=GetUnsetPoint();
    if (iPoint == -1)
      break;
    int iPointWork=iPoint;
    while(true) {
      StatusOrbit[iPointWork]=idxOrbit;
      iPointWork = ePerm.at(iPointWork);
      if (iPointWork == iPoint) {
	break;
      }
    }
    idxOrbit++;
  }
  return StatusOrbit;
}


#endif
