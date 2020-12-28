#ifndef INCLUDE_GroupFct_h
#define INCLUDE_GroupFct_h


#include "Temp_common.h"
#include "Basic_file.h"
#include "Basic_string.h"
#include "NumberTheory.h"
#include "Boost_bitset.h"


#include <permlib/permlib_api.h>


typedef std::shared_ptr<permlib::PermutationGroup> PermutationGroupPtr;
typedef boost::dynamic_bitset<> DsetList;


std::vector<int> PermutationOrbit(permlib::Permutation const& ePerm)
{
  //  std::cerr << "  Beginning of PermutationOrbit\n";
  int siz=ePerm.size();
  //  for (int i=0; i<siz; i++) {
  //    std::cerr << "i=" << i << " img=" << ePerm.at(i) << "\n";
  //  }
  std::vector<int> StatusOrbit(siz,-1);
  int idxOrbit=0;
  auto GetUnsetPoint=[&](void) -> int {
    for (int i=0; i<siz; i++)
      if (StatusOrbit[i] == -1)
	return i;
    return -1;
  };
  while(1) {
    int iPoint=GetUnsetPoint();
    //    std::cerr << "iPoint=" << iPoint << " idxOrbit=" << idxOrbit << "\n";
    if (iPoint == -1)
      break;
    int iPointWork=iPoint;
    while(1) {
      StatusOrbit[iPointWork]=idxOrbit;
      iPointWork=ePerm.at(iPointWork);
      //      std::cerr << "  iPointWork=" << iPointWork << "\n";
      if (iPointWork == iPoint)
	break;
    }
    idxOrbit++;
  }
  //  std::cerr << "  End of PermutationOrbit\n";
  return StatusOrbit;
}


#endif
