#ifndef INCLUDE_MapCoordinateFinder_h
#define INCLUDE_MapCoordinateFinder_h


#include "PlaneGraph.h"
#include "SVGfunctions.h"


/* Adapted from the code schlegel.c of cage */


template<typename Tgr>
std::vector<double> TD_EvaluatePosition(std::vector<double> const& ListRadius, Tgr const& eG)
{
  int nbVert=ListRadius.size();
  std::vector<double> ListDefect(nbVert);
  for (int iVert=0; iVert<nbVert; iVert++) {
    double H=M_PI;
    std::vector<size_t> ListAdj=eG.Adjacency(iVert);
    for (auto & jVert : ListAdj)
      H -= atan(ListRadius[jVert] / ListRadius[iVert]);
    ListDefect[iVert]=H;
  }
  return ListDefect;
}


double TD_SquareDefect(std::vector<double> const& ListDefect)
{
  double eDefect=0;
  int nbVert=ListDefect.size();
  for (int iVert=0; iVert<nbVert; iVert++) {
    double eVal=ListDefect[iVert];
    eDefect += eVal*eVal;
  }
  return eDefect;
}


struct TD_result {
  bool test;
  std::vector<double> ListRad;
};


template<typename Tgr>
double TD_ErrorOfListRadius(std::vector<double> const& ListRadius, Tgr const& eG)
{
  return TD_SquareDefect(TD_EvaluatePosition(ListRadius, eG));
}


template<typename Tgr>
TD_result TD_Solv_DecreaseIncrease(std::vector<double> const& ListRadius, Tgr const& eG)
{
  std::vector<double> ListDefect=TD_EvaluatePosition(ListRadius, eG);
  int nbVert=ListRadius.size();
  std::vector<double> dirEvol(nbVert);
  for (int iVert=0; iVert<nbVert; iVert++) {
    double eVal=0;
    if (ListDefect[iVert] > 0)
      eVal=-1;
    if (ListDefect[iVert] < 0)
      eVal= 1;
    dirEvol[iVert]=eVal;
  }
  double delta=0.5; // value so that all terms are well defined
  int maxite=100;
  std::vector<double> NewVect(nbVert);
  std::vector<double> Vect=ListRadius;
  for (int ite=0; ite<maxite; ite++) {
    std::vector<double> EVOL(nbVert);
    for (int iVert=0; iVert<nbVert; iVert++)
      EVOL[iVert] = 1 + delta*dirEvol[iVert];
    //    std::cerr << " |EVOL|=" << EVOL.size() << "\n";
    int nbmove=0;
    while(1) {
      //      std::cerr << " |Vect|=" << Vect.size() << "\n";
      double Norm=TD_ErrorOfListRadius(Vect, eG);
      //      std::cerr << " Norm=" << Norm << "\n";
      for (int iVert=0; iVert<nbVert; iVert++)
	NewVect[iVert]=Vect[iVert] * EVOL[iVert];
      //      std::cerr << " Norm=" << Norm << "\n";
      double NewNorm=TD_ErrorOfListRadius(NewVect, eG);
      if (NewNorm > Norm)
	break;
      Vect=NewVect;
      nbmove++;
    }
    if (nbmove>0)
      return {true, Vect};
    delta /= double(2);
  }
  return {false, ListRadius};
}


struct Point {
  bool fixed;
  coor disp;
  coor v;
};


struct SystemSolution {
  std::vector<int> FixedCycle;
  std::vector<Point> ListVert;
  std::vector<int> ListTriangles;
};


void PLANE_FindCoordinates(SystemSolution & eSol, int const& MAX_ITERATIONS, std::function<double(double const&)> const&AreaToKoef)
{
  int nvertices=eSol.ListVert.size();
  int nbTriangle=eSol.ListTriangles.size()/3;
  coor center;
  //
  std::cerr << "nvertices=" << nvertices << " nbTriangle=" << nbTriangle << "\n";
  for (int i=0; i<nvertices; i++) {
    eSol.ListVert[i].fixed=false;
    eSol.ListVert[i].v.x = 0;
    eSol.ListVert[i].v.y = 0;
  }
  int lenCycle=eSol.FixedCycle.size();
  std::cerr << "lenCycle=" << lenCycle << "\n";
  for (int i=0; i<lenCycle; i++) {
    int iVert=eSol.FixedCycle[i];
    eSol.ListVert[iVert].fixed=true;
    eSol.ListVert[iVert].v.x = 1000*cos(2 * M_PI * i/double(lenCycle));
    eSol.ListVert[iVert].v.y = 1000*sin(2 * M_PI * i/double(lenCycle));
  }
  std::vector<double> ListArea(nbTriangle);
  for (int step = 1; step <= MAX_ITERATIONS; step++) {
    for (int i = 0; i < nvertices; i++)
      eSol.ListVert[i].disp = {0, 0};
    //    std::cerr << "Begin triangle loop, nbTriangle=" << nbTriangle << "\n";
    for (int iTrig = 0; iTrig < nbTriangle; iTrig++) {
      int a=eSol.ListTriangles[3*iTrig];
      int b=eSol.ListTriangles[3*iTrig + 1];
      int c=eSol.ListTriangles[3*iTrig + 2];
      //      std::cerr << "a=" << a << " b=" << b << " c=" << c << "\n";
      coor v_a = eSol.ListVert[a].v;
      coor v_b = eSol.ListVert[b].v;
      coor v_c = eSol.ListVert[c].v;
      //      std::cerr << "After coor v_a, v_b and v_c assignations\n";

      double area = 0.5 * (   ( v_b.y - v_a.y ) * ( v_c.x - v_a.x )
			      - ( v_b.x - v_a.x ) * ( v_c.y - v_a.y ) );
      ListArea[iTrig]=area;
      double koef = AreaToKoef(area);

      center.x = ( v_a.x + v_b.x + v_c.x ) / 3.0;
      center.y = ( v_a.y + v_b.y + v_c.y ) / 3.0;

      eSol.ListVert[a].disp.x += koef * (center.x - v_a.x);
      eSol.ListVert[a].disp.y += koef * (center.y - v_a.y);

      eSol.ListVert[b].disp.x += koef * (center.x - v_b.x);
      eSol.ListVert[b].disp.y += koef * (center.y - v_b.y);

      eSol.ListVert[c].disp.x += koef * (center.x - v_c.x);
      eSol.ListVert[c].disp.y += koef * (center.y - v_c.y);
    }
    //    std::cerr << "After triangle loop\n";

    double temp = 40.0 / exp( 4.0 * (double)step / ( (step<250) ? 250 : (step+1) ) );
    //    std::cerr << "temp=" << temp << "\n";
    for (int i = 0; i < nvertices; i++) {
      if (eSol.ListVert[i].fixed == false) {
	double d = norm(eSol.ListVert[i].disp);
	double koef;
	if (d > temp)
	  koef = temp / d;
	else
	  koef = 1.0;
	eSol.ListVert[i].v.x += eSol.ListVert[i].disp.x * koef;
	eSol.ListVert[i].v.y += eSol.ListVert[i].disp.y * koef;
      }
    }
  }
  std::vector<double> ListX(nvertices), ListY(nvertices);
  for (int i = 0; i < nvertices; i++) {
    eSol.ListVert[i].v.x = eSol.ListVert[i].v.x /1000;
    eSol.ListVert[i].v.y = eSol.ListVert[i].v.y /1000;
    ListX[i]=eSol.ListVert[i].v.x;
    ListY[i]=eSol.ListVert[i].v.y;
  }
  int nbPlus=0;
  int nbMinus=0;
  for (int iTrig = 0; iTrig < nbTriangle; iTrig++) {
    if (ListArea[iTrig] > 0)
      nbPlus++;
    else
      nbMinus++;
  }
  std::cerr << "nbPlus=" << nbPlus << " nbMinus=" << nbMinus << "\n";
  std::cerr << "X(min/max)=" << VectorMin(ListX) << " / " << VectorMax(ListX) << "\n";
  std::cerr << "Y(min/max)=" << VectorMin(ListY) << " / " << VectorMax(ListY) << "\n";
}


struct PLANE_plot_infos {
  std::vector<coor> ListCoord;
  std::vector<int> ListVertStatus;
};


PLANE_plot_infos PLANE_RotationSphericalStructure(PLANE_plot_infos const& eStr, double const& ang)
{
  double eCos=cos(ang);
  double eSin=sin(ang);
  int nbVert=eStr.ListCoord.size();
  std::vector<coor> RetListCoord(nbVert);
  for (int iVert=0; iVert<nbVert; iVert++) {
    double x1=eStr.ListCoord[iVert].x;
    double y1=eStr.ListCoord[iVert].y;
    double x2= eCos*x1 + eSin*y1;
    double y2=-eSin*x1 + eCos*y1;
    RetListCoord[iVert]={x2,y2};
  }
  return {RetListCoord,eStr.ListVertStatus};
}


PLANE_plot_infos PLANE_ComputeCoordinateVertices(VEForiented const& VEFori, std::vector<int> const& ListKillVert, std::vector<int> const& TheCycle, int const& MAX_ITERATIONS, std::function<double(double const&)> const& AreaToKoef)
{
  int nbVert=VEFori.nbVert;
  std::vector<int> ListStatus(nbVert,1);
  for (auto & eVert : ListKillVert)
    ListStatus[eVert]=0;
  int nbActive=0;
  for (int i=0; i<nbVert; i++)
    nbActive += ListStatus[i];
  std::vector<int> ListMap(nbActive);
  std::vector<int> ListRevMap(nbVert,-1);
  int idx=0;
  for (int iVert=0; iVert<nbVert; iVert++)
    if (ListStatus[iVert] == 1) {
      ListRevMap[iVert]=idx;
      ListMap[idx]=iVert;
      idx++;
    }
  int nbFace=VEFori.nbFace;
  std::vector<int> ListStatusFace(nbFace);
  std::vector<int> ListTriangles;
  for (int iFace=0; iFace<nbFace; iFace++) {
    int eCorrVal=1;
    if (VEFori.FaceSet[iFace].size() != 3) {
      std::cerr << "We need a triangulation to make things work\n";
      throw TerminalException{1};
    }
    for (int i=0; i<3; i++) {
      int iDE=VEFori.FaceSet[iFace][i];
      int iVert=VEFori.ListOriginVert[iDE];
      eCorrVal *= ListStatus[iVert];
    }
    ListStatusFace[iFace]=eCorrVal;
    if (eCorrVal == 1) {
      for (int i=0; i<3; i++) {
	int iDE=VEFori.FaceSet[iFace][i];
	int iVert=VEFori.ListOriginVert[iDE];
	int iVertRed=ListRevMap[iVert];
	if (iVertRed == -1) {
	  std::cerr << "iVertRed should never equal -1\n";
	  throw TerminalException{1};
	}
	ListTriangles.push_back(iVertRed);
      }
    }
  }
  std::vector<int> FixedCycle;
  for (auto & eVert : TheCycle) {
    int eVertRed=ListRevMap[eVert];
    FixedCycle.push_back(eVertRed);
  }
  std::vector<Point> ListVert(nbActive);
  SystemSolution eSol{FixedCycle, ListVert, ListTriangles};
  PLANE_FindCoordinates(eSol, MAX_ITERATIONS, AreaToKoef);
  std::vector<coor> ListVertRet(nbVert);
  for (int iAct=0; iAct<nbActive; iAct++) {
    int eVert=ListMap[iAct];
    ListVertRet[eVert] = eSol.ListVert[iAct].v;
  }
  return {ListVertRet, ListStatus};
}


template<typename Tgr>
TD_result TD_NewtonMethod(Tgr const& eG, std::vector<double> const& ListRadius)
{
  int nbVert=ListRadius.size();
  MyMatrix<double> eMat=ZeroMatrix<double>(nbVert, nbVert);
  MyVector<double> delta(nbVert);
  //  std::cerr << "TD_NewtonMethod, step 1\n";
  for (int uVert=0; uVert<nbVert; uVert++) {
    double rad1=ListRadius[uVert];
    std::vector<size_t> ListAdj=eG.Adjacency(uVert);
    double eDelta=M_PI;
    for (auto & vVert : ListAdj) {
      double rad2=ListRadius[vVert];
      double SQV=rad1*rad1 + rad2*rad2;
      eMat(uVert,uVert) += -rad2 / SQV;
      eMat(uVert,vVert) +=  rad1 / SQV;
      eDelta -= atan(rad2 / rad1);
    }
    delta(uVert)=eDelta;
  }
  //  std::cerr << "TD_NewtonMethod, step 2\n";
  MyMatrix<double> PMat=ZeroMatrix<double>(nbVert,nbVert-1);
  for (int iVert=0; iVert<nbVert-1; iVert++) {
    PMat(iVert,iVert)=1;
    PMat(iVert+1,iVert)=-1;
  }
  //  std::cerr << "TD_NewtonMethod, step 3\n";
  MyMatrix<double> eProd=eMat*PMat;
  MyMatrix<double> eProdRed(nbVert-1,nbVert-1);
  MyVector<double> deltaRed(nbVert-1);
  for (int iRow=0; iRow<nbVert-1; iRow++) {
    deltaRed(iRow)=delta(iRow);
    for (int iCol=0; iCol<nbVert-1; iCol++)
      eProdRed(iRow,iCol)=eProd(iRow,iCol);
  }
  //  std::cerr << "TD_NewtonMethod, step 4\n";
  Eigen::FullPivLU<MyMatrix<double>> solver;
  solver.compute(eProdRed);
  MyVector<double> eSolRed=solver.solve(deltaRed);
  MyVector<double> eSol=PMat*eSolRed;
  double pow=1.2;
  double fact=1;
  for (int i=1; i<10; i++) {
    std::vector<double> NewListRadius(nbVert);
    for (int iV=0; iV<nbVert; iV++)
      NewListRadius[iV]=ListRadius[iV] + eSol(iV) / fact;
    //  std::cerr << "TD_NewtonMethod, step 5\n";
    double TheMin=VectorMin(NewListRadius);
    std::cerr << "TheMin=" << TheMin << " fact=" << fact << "\n";
    if (TheMin > 0) {
      double SQR1=TD_ErrorOfListRadius(ListRadius, eG);
      double SQR2=TD_ErrorOfListRadius(NewListRadius, eG);
      std::cerr << "SQR1=" << SQR1 << " SQR2=" << SQR2 << "\n";
      if (SQR2 < SQR1)
	return {true, NewListRadius};
    }
    fact *= pow;
  }
  return {false, {}};
}


template<typename Tgr>
TD_result TD_SeveralMinimization(Tgr const& eG, double const& minimal, int const& MAX_ITERATIONS)
{
  int nbVert=eG.GetNbVert();
  std::vector<double> ListRadius(nbVert,1);
  int iter=0;
  while(1) {
    double SQR=TD_ErrorOfListRadius(ListRadius, eG);
    std::cerr << "iter=" << iter << " nbVert=" << nbVert << " SQR=" << SQR << "\n";
    if (SQR < minimal)
      return {true, ListRadius};
    TD_result eRes1=TD_Solv_DecreaseIncrease(ListRadius, eG);
    if (eRes1.test == false) {
      std::cerr << "Our optimization procedure seems to have failed\n";
      throw TerminalException{1};
    }
    ListRadius=eRes1.ListRad;
    TD_result eRes2=TD_NewtonMethod(eG, ListRadius);
    if (eRes2.test)
      ListRadius=eRes2.ListRad;
    double TheMin=VectorMin(ListRadius);
    for (int iVert=0; iVert<nbVert; iVert++)
      ListRadius[iVert] /= TheMin;
    iter++;
    if (MAX_ITERATIONS != -1) {
      if (iter > MAX_ITERATIONS)
	return {false, ListRadius};
    }
  }
}


struct TD_description {
  std::vector<coor> ListCoord;
  VEForiented VEFori;
  MyMatrix<double> TransLatt;
  MyMatrix<int> ShiftMatrix;
};


MyMatrix<double> TD_RotateBasis(MyMatrix<double> const& TransM, double const& AngRad)
{
  double cosAng=cos(AngRad);
  double sinAng=sin(AngRad);
  MyMatrix<double> TransLatt(2,2);
  for (int i=0; i<2; i++) {
    double x=TransM(0,i);
    double y=TransM(1,i);
    double newCor1=cosAng*x - sinAng*y;
    double newCor2=sinAng*x + cosAng*y;
    TransLatt(0,i)=newCor1;
    TransLatt(1,i)=newCor2;
  }
  return TransLatt;
}


TD_description TD_SimilitudeOperation(TD_description const& eDesc, double const& AngDeg, double const& eScal, double const& shiftX, double const& shiftY)
{
  double AngRad=M_PI*(AngDeg/double(180));
  double cosAng=cos(AngRad);
  double sinAng=sin(AngRad);
  int nbVert=eDesc.ListCoord.size();
  std::vector<coor> ListCoordTot(nbVert);
  for (int iVert=0; iVert<nbVert; iVert++) {
    coor pt=eDesc.ListCoord[iVert];
    double newCor1=eScal * ( cosAng*pt.x - sinAng*pt.y) + shiftX;
    double newCor2=eScal * ( sinAng*pt.x + cosAng*pt.y) + shiftY;
    coor Npt{newCor1, newCor2};
    //    std::cerr << "iVert=" << iVert << " pt.x=" << Npt.x << " pt.y=" << Npt.y << "\n";
    ListCoordTot[iVert]=Npt;
  }
  MyMatrix<double> TransLatt=eScal * TD_RotateBasis(eDesc.TransLatt, AngRad);
  return {ListCoordTot, eDesc.VEFori, TransLatt, eDesc.ShiftMatrix};
}


TD_description TD_CanonicalizationOperation(TD_description const& eDesc)
{
  int iVector=-1;
  double eNormFind=0;
  for (int i=0; i<2; i++) {
    double x=eDesc.TransLatt(0,i);
    double y=eDesc.TransLatt(1,i);
    double eNorm=sqrt(x*x + y*y);
    if (eNorm > eNormFind) {
      iVector=i;
      eNormFind=eNorm;
    }
  }
  double xV=eDesc.TransLatt(0,iVector);
  double yV=eDesc.TransLatt(1,iVector);
  double AngRad= - atan2(yV, xV);
  double AngDeg=AngRad*(double(180) / M_PI);
  MyMatrix<double> TransLatt=TD_RotateBasis(eDesc.TransLatt, AngRad);
  std::cerr << "TransLatt=\n";
  WriteMatrix(std::cerr, TransLatt);
  std::vector<double> ListX{double(0), TransLatt(0,0), TransLatt(0,0) + TransLatt(0,1), TransLatt(0,1)};
  std::vector<double> ListY{double(0), TransLatt(1,0), TransLatt(1,0) + TransLatt(1,1), TransLatt(1,1)};
  double eXmin=VectorMin(ListX);
  double eXmax=VectorMax(ListX);
  double eYmin=VectorMin(ListY);
  double eYmax=VectorMax(ListY);
  // Now equation to solve
  double diffX=eXmax - eXmin;
  double diffY=eYmax - eYmin;
  double scal;
  if (diffX < diffY)
    scal=1/diffY;
  else
    scal=1/diffX;
  std::cerr << "AngDeg=" << AngDeg << " scal=" << scal << "\n";
  return TD_SimilitudeOperation(eDesc, AngDeg, scal, double(0), double(0));
}


std::vector<double> TD_GetListTriangleArea(TD_description const& eDesc)
{
  int nbFace=eDesc.VEFori.nbFace;
  std::vector<double> ListArea(nbFace);
  for (int iFace=0; iFace<nbFace; iFace++) {
    std::vector<int> eFace=eDesc.VEFori.FaceSet[iFace];
    int len=eFace.size();
    std::vector<coor> ListCoor(len);
    std::vector<int> ListIVert(len);
    int j0=0, j1=0;
    for (int i=0; i<len; i++) {
      int iDE=eFace[i];
      int aVert=eDesc.VEFori.ListOriginVert[iDE];
      coor eCoord=eDesc.ListCoord[aVert];
      eCoord.x += j0 * eDesc.TransLatt(0,0) + j1 * eDesc.TransLatt(0,1);
      eCoord.y += j0 * eDesc.TransLatt(1,0) + j1 * eDesc.TransLatt(1,1);
      ListIVert[i]=aVert;
      ListCoor[i]=eCoord;
      j0 += eDesc.ShiftMatrix(iDE,0);
      j1 += eDesc.ShiftMatrix(iDE,1);
    }
    coor center=IsobarycenterPoint(ListCoor);
    double area=0;
    for (int i=0; i<len; i++) {
      int j=NextIdx(len,i);
      coor v_a = center;
      coor v_b = ListCoor[i];
      coor v_c = ListCoor[j];
      double areaTri=0.5 * ( ( v_b.y - v_a.y ) * ( v_c.x - v_a.x )
			     - ( v_b.x - v_a.x ) * ( v_c.y - v_a.y ) );
      area += areaTri;
    }
    ListArea[iFace] = T_abs(area);
  }
  return ListArea;
}


double TD_GetMaximalAreaFraction(TD_description const& eDesc)
{
  std::vector<double> ListArea=TD_GetListTriangleArea(eDesc);
  double maxArea=VectorMax(ListArea);
  double minArea=VectorMin(ListArea);
  double fracArea=maxArea/minArea;
  return fracArea;
}


TD_description TD_CaGeProcessOptimization(TD_description const& eDesc, int const& MAX_ITERATIONS, std::function<double(double const&)> const& AreaToKoef)
{
  //  int nbP=eDesc.VEFori.ListOriginVert.size();
  std::cerr << "Begin TD_CaGeProcessOptimization\n";
  int nbVert=eDesc.VEFori.nbVert;
  int nbFace=eDesc.VEFori.nbFace;
  std::vector<int> ListSiz(nbFace);
  for (int iFace=0; iFace<nbFace; iFace++) {
    int len=eDesc.VEFori.FaceSet[iFace].size();
    ListSiz[iFace]=len;
  }
  int lenMax=VectorMax(ListSiz);
  std::vector<int> ListAtt(lenMax+1,0);
  for (auto & eSize : ListSiz)
    ListAtt[eSize]++;
  for (int i=0; i<=lenMax; i++) {
    int eNB=ListAtt[i];
    if (eNB > 0)
      std::cerr << " f" << i << "=" << eNB;
  }
  std::cerr << "\n";

  std::vector<coor> ListCoord=eDesc.ListCoord;
  std::vector<double> ListArea(nbFace);
  MyMatrix<double> DispMat(nbVert,2);
  for (int step = 1; step <= MAX_ITERATIONS; step++) {
    for (int i = 0; i < nbVert; i++)
      for (int iCol=0; iCol<2; iCol++)
	DispMat(i,iCol)=0;
    for (int iFace=0; iFace<nbFace; iFace++) {
      std::vector<int> eFace=eDesc.VEFori.FaceSet[iFace];
      int len=eFace.size();
      std::vector<coor> ListCoor(len);
      std::vector<int> ListIVert(len);
      int j0=0, j1=0;
      for (int i=0; i<len; i++) {
	int iDE=eFace[i];
	int aVert=eDesc.VEFori.ListOriginVert[iDE];
	coor eCoord=ListCoord[aVert];
	eCoord.x += j0 * eDesc.TransLatt(0,0) + j1 * eDesc.TransLatt(0,1);
	eCoord.y += j0 * eDesc.TransLatt(1,0) + j1 * eDesc.TransLatt(1,1);
	ListIVert[i]=aVert;
	ListCoor[i]=eCoord;
	j0 += eDesc.ShiftMatrix(iDE,0);
	j1 += eDesc.ShiftMatrix(iDE,1);
      }
      coor center=IsobarycenterPoint(ListCoor);
      double area=0;
      for (int i=0; i<len; i++) {
	int j=NextIdx(len,i);
	coor v_a = center;
	coor v_b = ListCoor[i];
	coor v_c = ListCoor[j];
	double areaTri=0.5 * ( ( v_b.y - v_a.y ) * ( v_c.x - v_a.x )
			       - ( v_b.x - v_a.x ) * ( v_c.y - v_a.y ) );
	area += areaTri;
      }
      ListArea[iFace] = area;
      double koef=AreaToKoef(area);
      for (int i=0; i<len; i++) {
	int a=ListIVert[i];
	coor v_a=ListCoor[i];
	DispMat(a,0) += koef * (center.x - v_a.x);
	DispMat(a,1) += koef * (center.y - v_a.y);
      }
    }

    double temp = 40.0 / exp( 4.0 * (double)step / ( (step<250) ? 250 : (step+1) ) );
    for (int i = 0; i < nbVert; i++) {
      coor vec{DispMat(i,0),DispMat(i,1)};
      double d = norm(vec);
      double koef;
      if (d > temp)
	koef = temp / d;
      else
	koef = 1.0;
      ListCoord[i].x += DispMat(i,0) * koef;
      ListCoord[i].y += DispMat(i,1) * koef;
    }
  }
  int nbPlus=0;
  int nbMinus=0;
  for (int iFace = 0; iFace < nbFace; iFace++) {
    if (ListArea[iFace] > 0)
      nbPlus++;
    else
      nbMinus++;
  }
  std::cerr << "nbPlus=" << nbPlus << " nbMinus=" << nbMinus << "\n";
  return {ListCoord, eDesc.VEFori, eDesc.TransLatt, eDesc.ShiftMatrix};
}


template<typename T>
T L1_norm(MyVector<T> const& V)
{
  int len=V.size();
  T eNorm = 0;
  for (int i=0; i<len; i++)
    eNorm += T_abs(V(i));
  return eNorm;
}


template<typename Telt>
TD_description TD_FindDescription(PlanGraphOriented<Telt> const& PL, double const& minimal, int const& MAX_ITERATIONS, double const& tol)
{
  Telt eInv = PL.invers;
  Telt eNext = PL.next;
  Telt ePrev = Inverse(eNext);
  int nbP = PL.nbP;
  auto VEFori = PlanGraphOrientedToVEF(PL);
  int nbVert=VEFori.nbVert;
  GraphSparseImmutable eG=PlanGraphOrientedToGSI(VEFori);
  TD_result eRes=TD_SeveralMinimization(eG, minimal, MAX_ITERATIONS);
  if (eRes.test == false) {
    std::cerr << "No convergence has been achieved in the coordinate determination\n";
    throw TerminalException{1};
  }
  MyMatrix<double> ListCoordDE(nbP,5);
  std::vector<int> ListStatus(nbVert,0);
  auto SetOneVertex=[&](int const& iDE, double const& eX, double const& eY, double const& eA) -> void {
    int iVert=VEFori.ListOriginVert[iDE];
    ListStatus[iVert]=1;
    ListCoordDE(iDE,0)=eX;
    ListCoordDE(iDE,1)=eY;
    ListCoordDE(iDE,2)=eA;
    int iDE1=iDE;
    while(1) {
      int iDE2=eNext.at(iDE1);
      if (iDE2 == iDE)
	break;
      int iDE1rev=PL.invers.at(iDE1);
      int iDE2rev=PL.invers.at(iDE2);
      int iVert1=VEFori.ListOriginVert[iDE1rev];
      int iVert2=VEFori.ListOriginVert[iDE2rev];
      double ang1=atan(eRes.ListRad[iVert1] / eRes.ListRad[iVert]);
      double ang2=atan(eRes.ListRad[iVert2] / eRes.ListRad[iVert]);
      ListCoordDE(iDE2,0)=eX;
      ListCoordDE(iDE2,1)=eY;
      ListCoordDE(iDE2,2)=ListCoordDE(iDE1,2) + ang1 + ang2;
      iDE1=iDE2;
    }
  };
  SetOneVertex(0, 0, 0, 0);
  std::vector<MyVector<double> > ListVect;
  double TotalErrorAngle=0;
  std::cerr << "Before determination of all the positions\n";
  auto AngNormalization=[](double const& x) -> double {
    double xRet = x;
    while(true) {
      bool DoOper=false;
      if (xRet < -M_PI) {
	DoOper=true;
	xRet += 2*M_PI;
      }
      if (xRet > M_PI) {
	DoOper=true;
	xRet -= 2*M_PI;
      }
      if (DoOper == false) {
	return xRet;
      }
    }
  };
  while(true) {
    int nbOper=0;
    for (int iVert=0; iVert<nbVert; iVert++)
      if (ListStatus[iVert] == 1) {
	nbOper++;
	ListStatus[iVert]=2;
	double rad1=eRes.ListRad[iVert];
	for (auto & iDE : VEFori.VertSet[iVert]) {
	  double eX=ListCoordDE(iDE,0);
	  double eY=ListCoordDE(iDE,1);
	  double eA=ListCoordDE(iDE,2);
	  int rDE=eInv.at(iDE);
	  int jVert=VEFori.ListOriginVert[rDE];
	  double rad2=eRes.ListRad[jVert];
	  double dist=sqrt(rad1*rad1 + rad2*rad2);
	  double NewX=eX + dist*cos(eA);
	  double NewY=eY + dist*sin(eA);
	  double NewA=eA + M_PI;
	  if (ListStatus[jVert] == 0)
	    SetOneVertex(rDE, NewX, NewY, NewA);
	  double deltaX = NewX - ListCoordDE(rDE,0);
	  double deltaY = NewY - ListCoordDE(rDE,1);
	  double deltaA = NewA - ListCoordDE(rDE,2);
	  double AngRenormA=AngNormalization(deltaA);
	  //	  std::cerr << "deltaA=" << deltaA << " AngRenorm=" << AngRenormA << "\n";
	  TotalErrorAngle += fabs(AngRenormA);
	  ListCoordDE(iDE,3) = deltaX;
	  ListCoordDE(iDE,4) = deltaY;
	}
      }
    //    std::cerr << "nbOper=" << nbOper << "\n";
    if (nbOper == 0) {
      break;
    }
  }
  std::cerr << "TotalErrorAngle = " << TotalErrorAngle << "\n";
  std::vector<coor> ListCoordTot(nbVert);
  for (int iVert=0; iVert<nbVert; iVert++) {
    int iDE=VEFori.VertSet[iVert][0];
    ListCoordTot[iVert]={ListCoordDE(iDE,0), ListCoordDE(iDE,1)};
  }
  //
  // Now constructing the lattice
  //
  std::vector<MyVector<double> > ListV;
  auto fInsert=[&](MyVector<double> const& V) -> void {
    if (L1_norm(V) < tol)
      return;
    for (auto & eV : ListV) {
      MyVector<double> eDiff = eV - V;
      if (L1_norm(eDiff) < tol)
	return;
    }
    ListV.push_back(V);
  };
  MyVector<double> eDiff(2);
  for (int iDE=0; iDE<nbP; iDE++) {
    for (int i=0; i<2; i++)
      eDiff(i) = ListCoordDE(iDE,i+3);
    fInsert(eDiff);
  }
  double eScalChosen=double(5);
  MyMatrix<double> TheBasis;
  bool FindOneBasis=false;
  int nbVect=ListV.size();
  for (int iVect=0; iVect<nbVect-1; iVect++)
    for (int jVect=iVect+1; jVect<nbVect; jVect++) {
      MyVector<double> U1=ListV[iVect];
      MyVector<double> U2=ListV[jVect];
      double eN1=U1(0)*U1(0) + U1(1)*U1(1);
      double eN2=U2(0)*U2(0) + U2(1)*U2(1);
      double eSc=U1(0)*U2(0) + U1(1)*U2(1);
      double eDefect=eSc*eSc / (eN1*eN2);
      MyMatrix<double> TheMat(2,2);
      for (int i=0; i<2; i++) {
	TheMat(i,0) = U1(i);
	TheMat(i,1) = U2(i);
      }
      MyMatrix<double> Hinv=Inverse(TheMat);
      bool test=true;
      double sumDeltaErr=0;
      for (int k=0; k<nbVect; k++)
	if (test) {
	  MyVector<double> eCoord=Hinv*ListV[k];
	  double delta=fabs(eCoord(0) - round(eCoord(0))) + fabs(eCoord(1) - round(eCoord(1)));
	  sumDeltaErr += delta;
	  if (delta > tol) {
	    test = false;
          }
	}
      std::cerr << "Vect(i/j)=" << iVect << " / " << jVect << " sumDeltaErr=" << sumDeltaErr << "\n";
      if (test && eDefect < eScalChosen) {
	FindOneBasis=true;
	TheBasis = TheMat;
	eScalChosen = eDefect;
      }
    }
  if (FindOneBasis == false) {
    std::cerr << "Failed to find an approximate basis\n";
    std::cerr << "nbVect=" << nbVect << "\n";
    throw TerminalException{1};
  }
  std::cerr << "Ending double loop over possible pairs of vectors\n";
  MyMatrix<double> HinvPerm=Inverse(TheBasis);
  MyMatrix<int> ShiftMatrix(nbP,2);
  std::cerr << "nbP=" << nbP << "\n";
  for (int iDE=0; iDE<nbP; iDE++) {
    MyVector<double> eVect(2);
    for (int i=0; i<2; i++)
      eVect(i) = ListCoordDE(iDE,i+3);
    MyVector<double> eCoord=HinvPerm*eVect;
    for (int i=0; i<2; i++)
      ShiftMatrix(iDE,i)=int(round(eCoord(i)));
  }
  std::cerr << "ShiftMatrix rows=" << ShiftMatrix.rows() << " cols=" << ShiftMatrix.cols() << "\n";
  std::cerr << "Before end of TD_FindDescription\n";
  return {ListCoordTot, VEFori, TheBasis, ShiftMatrix};
}


TD_description TD_FindDescriptionDirect(PlanGraphOriented const& PL, double const& minimal, int const& MAX_ITERATIONS, double const& tol)
{
  int nbP=PL.nbP;
  auto VEFori=PlanGraphOrientedToVEF(PL);
  InfoDualMedialGraph DualMedInfo=DualMedialGraphOriented(PL);
  TD_description eDesc=TD_FindDescription(DualMedInfo.PL, minimal, MAX_ITERATIONS, tol);
  int nbVert=VEFori.nbVert;
  std::vector<coor> ListCoordTot(nbVert);
  for (int iVert=0; iVert<nbVert; iVert++) {
    int iDE=VEFori.VertSet[iVert][0];
    int iDEdualmed=DualMedInfo.MappingDE(iDE,0);
    int iVertDualMed=eDesc.VEFori.ListOriginVert[iDEdualmed];
    coor pt=eDesc.ListCoord[iVertDualMed];
    //    std::cerr << "iVert=" << iVert << " x=" << pt.x << " y=" << pt.y << "\n";
    ListCoordTot[iVert]=pt;
  }
  MyMatrix<int> ShiftMatrix;
  std::cerr << "TD_FindDescriptionDirect nbP=" << nbP << "\n";
  ShiftMatrix.setZero(nbP,2);
  for (int iDE=0; iDE<nbP; iDE++)
    for (int i=0; i<2; i++) {
      int iDEnew=DualMedInfo.DecomposeOldDE(iDE,i);
      for (int j=0; j<2; j++)
	ShiftMatrix(iDE,j) += eDesc.ShiftMatrix(iDEnew,j);
    }
  std::cerr << "TD_FindDescriptionDirect rows=" << ShiftMatrix.rows() << " cols=" << ShiftMatrix.cols() << "\n";
  return {ListCoordTot, VEFori, eDesc.TransLatt, ShiftMatrix};
}


#endif
