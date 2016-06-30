#ifndef INCLUDE_Basic_string_h
#define INCLUDE_Basic_string_h


#include "Temp_common.h"


bool STRING_IsStringReduceToSpace(std::string const& eStr)
{
  int len=eStr.length();
  std::string eChar=" ";
  for (int i=0; i<len; i++) {
    std::string eSubChar=eStr.substr(i,1);
    if (eSubChar != eChar)
      return false;
  }
  return true;
}


int STRING_GetCharPositionInString(std::string const& eStr, std::string const& eChar)
{
  int len=eStr.length();
  for (int i=0; i<len; i++) {
    std::string eSubChar=eStr.substr(i,1);
    if (eSubChar == eChar)
      return i;
  }
  return -1;
}


std::string DoubleTo4dot2f(double const& x)
{
  char buffer[50];
  int n=sprintf(buffer, "%4.2f", x);
  if (n == 0) {
    std::cerr << "Clear error in DoubleTo4dot2f\n";
    throw TerminalException{1};
  }
  return std::string(buffer);
}


std::string DoubleToString(double const& x)
{
  std::stringstream s;
  s << x;
  std::string converted(s.str());
  return converted;
}


std::string IntToString(int const & x)
{
  std::stringstream s;
  s << x;
  std::string converted(s.str());
  return converted;
}


std::string STRING_RemoveSpacesBeginningEnd(std::string const& eStr)
{
  int len=eStr.size();
  std::vector<int> ListIsSpace(len,0);
  std::string eSpace=" ";
  for (int i=0; i<len; i++) {
    std::string eChar=eStr.substr(i, 1);
    if (eChar == eSpace)
      ListIsSpace[i]=1;
  }
  int PosLow=-1;
  for (int i=0; i<len; i++)
    if (PosLow == -1)
      if (ListIsSpace[i] == 0)
	PosLow=i;
  int PosUpp=-1;
  for (int i=0; i<len; i++) {
    int j=len-1-i;
    if (PosUpp == -1)
      if (ListIsSpace[j] == 0)
	PosUpp=j;
  }
  std::string RetStr;
  if (PosLow == -1) {
    return RetStr;
  }
  for (int iPos=PosLow; iPos<PosUpp+1; iPos++) {
    RetStr += eStr.at(iPos);
  }
  return RetStr;
}


std::vector<std::string> STRING_Split(std::string const& eStrA, std::string const& eStrB)
{
  int lenA=eStrA.length();
  int lenB=eStrB.length();
  std::vector<int> ListStatus(lenA,1);
  for (int iA=0; iA<lenA - lenB; iA++)
    if (ListStatus[iA] == 1) {
      bool IsMatch=true;
      for (int iB=0; iB<lenB; iB++) {
	std::string eCharA=eStrA.substr(iA+iB,1);
	std::string eCharB=eStrB.substr(iB,1);
	if (eCharA != eCharB)
	  IsMatch=false;
      }
      if (IsMatch)
	for (int iB=0; iB<lenB; iB++)
	  ListStatus[iA + iB]=0;
    }
  std::vector<std::string> RetList;
  std::string eFound;
  for (int iA=0; iA<lenA; iA++) {
    std::string eChar=eStrA.substr(iA, 1);
    if (ListStatus[iA] == 1)
      eFound += eChar;
    if (ListStatus[iA] == 0) {
      int siz=eFound.length();
      if (siz > 0)
	RetList.push_back(eFound);
      eFound="";
    }
  }
  int siz=eFound.size();
  if (siz > 0)
    RetList.push_back(eFound);
  return RetList;
}


#endif
