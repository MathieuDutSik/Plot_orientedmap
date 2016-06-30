#ifndef INCLUDE_Temp_Matrix_h
#define INCLUDE_Temp_Matrix_h


#include "Temp_common.h"


#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/LU>
#include <Eigen/Eigenvalues>


template <typename T>
using MyVector = Eigen::Matrix<T,Eigen::Dynamic,1>;


template <typename T>
using MyMatrix = Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>;


template<typename T>
MyMatrix<T> ZeroMatrix(int const& nbRow, int const& nbCol)
{
  MyMatrix<T> retMat(nbRow, nbCol);
  T eZero;
  eZero=0;
  for (int iRow=0; iRow<nbRow; iRow++)
    for (int iCol=0; iCol<nbCol; iCol++)
      retMat(iRow, iCol)=eZero;
  return retMat;
}


template<typename T>
void SwapValues(T& val1, T& val2)
{
  T prov;
  prov=val1;
  val1=val2;
  val2=prov;
}


template<typename T>
struct Inverse_exception {
  std::string errmsg;
  T pivot;
};


template<typename T>
void TMat_Inverse_destroy(MyMatrix<T> &Input, MyMatrix<T> &Output)
{
  int iCol, iRow;
  int iColB;
  int nbRow=Input.rows();
  int nbCol=Input.cols();
  T prov1, prov2, eVal;
  if (nbRow != nbCol) {
    std::cerr << "Error on nbRow, nbCol in TMat_Inverse_destroy";
    throw TerminalException{1};
  }
  for (iRow=0; iRow<nbRow; iRow++)
    for (iCol=0; iCol<nbRow; iCol++) {
      if (iRow == iCol)
	prov1=1;
      else
	prov1=0;
      Output(iRow,iCol)=prov1;
    }
  int iRowFound=-400;
  for (iCol=0; iCol<nbCol; iCol++) {
    prov1=0;
    for (iRow=iCol; iRow<nbRow; iRow++)
      if (prov1 == 0) {
	eVal=Input(iRow,iCol);
	if (eVal != 0) {
	  iRowFound=iRow;
	  prov1=1/eVal;
	}
      }
    if (prov1 == 0) {
      Inverse_exception<T> eExcept;
      eExcept.errmsg="Error in matrix inversion";
      eExcept.pivot=0;
      throw eExcept;
    }
    for (iColB=0; iColB<nbCol; iColB++) {
      Input(iRowFound,iColB) *= prov1;
      Output(iRowFound,iColB) *= prov1;
    }
    for (iRow=0; iRow<nbRow; iRow++)
      if (iRow != iRowFound) {
	prov2=Input(iRow, iCol);
	for (iColB=0; iColB<nbCol; iColB++) {
	  Input(iRow, iColB) -= prov2*Input(iRowFound,iColB);
	  //
	  Output(iRow,iColB) -= prov2*Output(iRowFound,iColB);
	}
      }
    if (iRowFound != iCol)
      for (iColB=0; iColB<nbCol; iColB++) {
	prov1=Input(iRowFound, iColB);
	prov2=Input(iCol, iColB);
	SwapValues(prov1, prov2);
	Input(iRowFound, iColB)=prov1;
	Input(iCol     , iColB)=prov2;
	//
	prov1=Output(iRowFound, iColB);
	prov2=Output(iCol, iColB);
	SwapValues(prov1, prov2);
	Output(iRowFound, iColB)=prov1;
	Output(iCol     , iColB)=prov2;
      }
  }
}


template<typename T>
MyMatrix<T> Inverse(MyMatrix<T> const&Input)
{
  int nbRow=Input.rows();
  MyMatrix<T> provMat=Input;
  MyMatrix<T> Output(nbRow, nbRow);
  TMat_Inverse_destroy(provMat, Output);
  return Output;
}


template<typename T>
std::vector<T> ReadStdVector(std::istream &is)
{
  T eVal;
  size_t nbRow;
  is >> nbRow;
  //  std::cerr << "nbRow=" << nbRow << "\n";
  std::vector<T> eVect(nbRow);
  for (size_t iRow=0; iRow<nbRow; iRow++) {
    is >> eVal;
    eVect[iRow]=eVal;
  }
  return eVect;
}


template<typename T>
void WriteMatrix(std::ostream &os, MyMatrix<T> const&TheMat)
{
  int nbRow=TheMat.rows();
  int nbCol=TheMat.cols();
  os << nbRow << " " << nbCol << "\n";
  for (int iRow=0; iRow<nbRow; iRow++) {
    for (int iCol=0; iCol<nbCol; iCol++) {
      T eVal=TheMat(iRow, iCol);
      os << " " << eVal;
    }
    os << "\n";
  }
}


template<typename T>
bool operator==(MyVector<T> const& V1, MyVector<T> const& V2)
{
  int n=V1.size();
  if (V2.size() != n) {
    std::cerr << "We should not have different sizes\n";
    throw TerminalException{1};
  }
  for (int i=0; i<n; i++) {
    if (V1(i) != V2(i))
      return false;
  }
  return true;
}


template<typename T>
bool operator<(MyMatrix<T> const& M1, MyMatrix<T> const& M2)
{
  int nbRow=M1.rows();
  int nbCol=M1.cols();
  for (int iRow=0; iRow<nbRow; iRow++)
    for (int iCol=0; iCol<nbCol; iCol++) {
      if (M1(iRow, iCol) < M2(iRow,iCol))
	return true;
      if (M1(iRow, iCol) > M2(iRow,iCol))
	return false;
    }
  return false;
}


template<typename T>
bool operator<(MyVector<T> const& V1, MyVector<T> const& V2)
{
  int siz=V1.size();
  for (int i=0; i<siz; i++) {
    if (V1(i) < V2(i))
      return true;
    if (V1(i) > V2(i))
      return false;
  }
  return false;
}


namespace std {
  template<typename T>
  struct less<MyVector<T>> {
    bool operator()(MyVector<T> const& V1, MyVector<T> const& V2) const
    {
      int siz=V1.size();
      for (int i=0; i<siz; i++) {
	if (V1(i) < V2(i))
	  return true;
	if (V2(i) < V1(i))
	  return false;
      }
      return false;
    }
  };
}


template<typename T>
T L1_norm(MyVector<T> const& V)
{
  int siz=V.size();
  T norm=0;
  for (int i=0; i<siz; i++)
    norm += T_abs(V(i));
  return norm;
}


#endif
