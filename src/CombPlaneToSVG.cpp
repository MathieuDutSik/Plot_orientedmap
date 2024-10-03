#include "PlaneGraphDrawing.h"
#include "Permutation.h"

int main(int argc, char *argv[])
{
  using Tidx = uint32_t;
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  if (argc != 2) {
    std::cerr << "Number of argument is = " << argc << "\n";
    std::cerr << "This program is used as\n";
    std::cerr << "CombPlaneToSVG [FileOption.nml]\n";
    std::cerr << "\n";
    std::cerr << "FileOption.nml : The filename containing the list of options\n";
    return -1;
  }
  //
  try {
    std::cerr << "Reading input\n";
    FullNamelist eFull=NAMELIST_GetStandardCombPlaneGraph();
    std::string eFileName=argv[1];
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    //
    CreateSVGfileOfGraph<Telt>(eFull);
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
  std::cerr << "Completion of the program\n";
}
