#ifndef INCLUDE_Basic_file_h
#define INCLUDE_Basic_file_h


#include "Temp_common.h"


#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>


//{
//  boost::filesystem::create_directories(eDir.c_str() );


bool IsExistingFile(std::string const& eFile)
{
  std::ifstream f(eFile.c_str());
  if (f.good()) {
    f.close();
    return true;
  } else {
    f.close();
    return false;
  }   
}


std::string FILE_RemoveEndingExtension(std::string const& FileName, std::string const& TheExtension)
{
  int len=FileName.size();
  int iCharLast=-1;
  for (int iChar=0; iChar<len; iChar++) {
    std::string eChar=FileName.substr(iChar,1);
    if (eChar == ".")
      iCharLast=iChar;
  }
  if (iCharLast == -1)
    return FileName;
  std::string FileNameRed=FileName.substr(0,iCharLast);
  std::string eExtension=FileName.substr(iCharLast+1,len-1-iCharLast);
  if (eExtension == TheExtension)
    return FileNameRed;
  return FileName;
}


#endif
