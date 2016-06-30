#ifndef INCLUDE_PlaneGraph_h
#define INCLUDE_PlaneGraph_h


#include "GroupFct.h"
#include "Combinatorics.h"
#include "GraphicalBasic.h"


struct PlanGraphOriented {
  int nbP;
  permlib::Permutation invers;
  permlib::Permutation next;
};


PlanGraphOriented ReadPlanGraphOrientedStream(std::istream & is)
{
  int nbP;
  is >> nbP;
  std::vector<permlib::dom_int> eListInv(nbP);
  std::vector<permlib::dom_int> eListNext(nbP);
  for (int i=0; i<nbP; i++) {
    int eVal;
    is >> eVal;
    eListInv[i]=eVal;
  }
  for (int i=0; i<nbP; i++) {
    int eVal;
    is >> eVal;
    eListNext[i]=eVal;
  }
  return {nbP, permlib::Permutation(eListInv), permlib::Permutation(eListNext)};
}


PlanGraphOriented ReadPlanGraphOrientedFile(std::string const& eFile)
{
  std::ifstream is(eFile);
  return ReadPlanGraphOrientedStream(is);
}


struct VEForiented {
  PlanGraphOriented PLori;
  int nbVert;
  int nbEdge;
  int nbFace;
  std::vector<std::vector<int>> VertSet;
  std::vector<std::vector<int>> EdgeSet;
  std::vector<std::vector<int>> FaceSet;
  std::vector<int> ListOriginVert;
  std::vector<int> ListOriginEdge;
  std::vector<int> ListOriginFace;
};


std::vector<std::vector<int> > ConvertListOriginToListList(std::vector<int> const& ListStat)
{
  int len=ListStat.size();
  int nbOrbit = 1 + VectorMax(ListStat);
  std::vector<std::vector<int> > eListList(nbOrbit);
  for (int i=0; i<len; i++)
    eListList[ListStat[i]].push_back(i);
  return eListList;
}


std::vector<std::vector<int> > GetListOrbit(std::vector<int> const& ListOrig, permlib::Permutation const& ePerm)
{
  int nbDE=ListOrig.size();
  int nbOrbit = 1 + VectorMax(ListOrig);
  std::vector<std::vector<int> > eListList(nbOrbit);
  auto TreatIorbit=[&](int const& eOrbit) -> void {
    for (int iDE=0; iDE<nbDE; iDE++)
      if (ListOrig[iDE] == eOrbit) {
	int iDEwork=iDE;
	while(1) {
	  eListList[eOrbit].push_back(iDEwork);
	  iDEwork=ePerm.at(iDEwork);
	  if (iDEwork == iDE)
	    return;
	}
      }
    std::cerr << "Failed to find matching\n";
    throw TerminalException{1};
  };
  for (int iOrbit=0; iOrbit<nbOrbit; iOrbit++)
    TreatIorbit(iOrbit);
  return eListList;
}


struct DEorder {
  int iDE;
  int side;
  int dir;
};


bool operator==(DEorder const& DE1, DEorder const& DE2)
{
  if (DE1.iDE != DE2.iDE)
    return false;
  if (DE1.side != DE2.side)
    return false;
  if (DE1.dir != DE2.dir)
    return false;
  return true;
}


struct DupliEdge {
  int iDE;
  int side;
};


bool operator==(DupliEdge const& DE1, DupliEdge const& DE2)
{
  if (DE1.iDE != DE2.iDE)
    return false;
  if (DE1.side != DE2.side)
    return false;
  return true;
}


struct DualMedialDE {
  int eDE1;
  int eDE2;
  int stat;
};


bool operator==(DualMedialDE const& DE1, DualMedialDE const& DE2)
{
  if (DE1.eDE1 != DE2.eDE1)
    return false;
  if (DE1.eDE2 != DE2.eDE2)
    return false;
  if (DE1.stat != DE2.stat)
    return false;
  return true;
}


struct DEvertface {
  int iDE;
  int nat;
};


VEForiented PlanGraphOrientedToVEF(PlanGraphOriented const& PL)
{
  permlib::Permutation eInv=PL.invers;
  permlib::Permutation eNext=PL.next;
  permlib::Permutation ePrev=~eNext;
  permlib::Permutation eFacePerm=eInv*ePrev;
  std::vector<int> ListOriginVert=PermutationOrbit(eNext);
  std::vector<int> ListOriginEdge=PermutationOrbit(eInv);
  std::vector<int> ListOriginFace=PermutationOrbit(eFacePerm);
  //  std::vector<std::vector<int> > VertSet=ConvertListOriginToListList(ListOriginVert);
  //  std::vector<std::vector<int> > EdgeSet=ConvertListOriginToListList(ListOriginEdge);
  //  std::vector<std::vector<int> > FaceSet=ConvertListOriginToListList(ListOriginFace);
  std::vector<std::vector<int> > VertSet=GetListOrbit(ListOriginVert, eNext);
  std::vector<std::vector<int> > EdgeSet=GetListOrbit(ListOriginEdge, eInv);
  std::vector<std::vector<int> > FaceSet=GetListOrbit(ListOriginFace, eFacePerm);
  int nbVert=VertSet.size();
  int nbEdge=EdgeSet.size();
  int nbFace=FaceSet.size();
  return {PL,
      nbVert, nbEdge, nbFace,
      VertSet, EdgeSet, FaceSet,
      ListOriginVert, ListOriginEdge, ListOriginFace};
}


int EulerPoincareCharacteristic(PlanGraphOriented const& PL)
{
  auto VEFori=PlanGraphOrientedToVEF(PL);
  int nbVertPL=VEFori.nbVert;
  int nbEdgePL=VEFori.nbEdge;
  int nbFacePL=VEFori.nbFace;
  int Charac=nbVertPL - nbEdgePL + nbFacePL;
  return Charac;
}


void PrintInformationOfMap(std::ostream & os, PlanGraphOriented const& PL)
{
  VEForiented VEFori=PlanGraphOrientedToVEF(PL);
  std::vector<int> VertLen;
  for (auto & eVert : VEFori.VertSet)
    VertLen.push_back(int(eVert.size()));
  os << "   nbVert=" << VEFori.VertSet.size() << " Coll(Vert)=";
  CollectedResult<int> eCollVert=Collected(VertLen);
  int lenVert=eCollVert.LVal.size();
  for (int i=0; i<lenVert; i++) {
    if (i>0)
      os << ", ";
    os << "[" << eCollVert.LVal[i] << "," << eCollVert.LMult[i] << "]";
  }
  os << "\n";
  os << "   nbEdge=" << VEFori.EdgeSet.size() << "\n";
  //
  std::vector<int> FaceLen;
  for (auto & eFace : VEFori.FaceSet)
    FaceLen.push_back(int(eFace.size()));
  os << "   nbFace=" << VEFori.FaceSet.size() << " Coll(Face)=";
  CollectedResult<int> eCollFace=Collected(FaceLen);
  int lenFace=eCollFace.LVal.size();
  for (int i=0; i<lenFace; i++) {
    if (i>0)
      os << ", ";
    os << "[" << eCollFace.LVal[i] << "," << eCollFace.LMult[i] << "]";
  }
  os << "\n";
}


struct PairWythoff123 {
  PlanGraphOriented PL;
  std::vector<DEorder> ListDE;
  MyMatrix<int> OldToNew;
  std::vector<int> ListDir;
  std::vector<int> ListSide;
  MyMatrix<int> SmallMat;
};


PairWythoff123 GetWythoff123(PlanGraphOriented const& PL)
{
  int nbP=PL.nbP;
  permlib::Permutation eNext=PL.next;
  permlib::Permutation ePrev=~eNext;
  int nbPtotal=6*nbP;
  std::vector<DEorder> ListDE(nbPtotal);
  std::vector<int> ListSide(nbPtotal);
  std::vector<int> ListDir(nbPtotal);
  MyMatrix<int> OldToNew(nbP,6);
  MyMatrix<int> SmallMat(2,3);
  int posB=0;
  for (int side=0; side<2; side++)
    for (int dir=0; dir<3; dir++) {
      SmallMat(side,dir)=posB;
      posB++;
    }
  int iDEnew=0;
  for (int iDE=0; iDE<nbP; iDE++) {
    int pos=0;
    for (int side=0; side<2; side++)
      for (int dir=0; dir<3; dir++) {
	ListDE[iDEnew]={iDE, side, dir};
	ListSide[iDEnew]=side;
	ListDir[iDEnew]=dir;
	OldToNew(iDE,pos)=iDEnew;
	iDEnew++;
	pos++;
      }
  }
  int nbDEnew=ListDE.size();
  auto GetPosition=[&](DEorder const& eDE) -> int {
    for (int iPos=0; iPos<nbDEnew; iPos++)
      if (ListDE[iPos] == eDE)
	return iPos;
    std::cerr << "Did not find the edge\n";
    throw TerminalException{1};
    //    exit(1);
  };
  std::vector<permlib::dom_int> NewNext(nbDEnew);
  std::vector<permlib::dom_int> NewInv (nbDEnew);
  for (int iDEnew=0; iDEnew<nbDEnew; iDEnew++) {
    DEorder eDE=ListDE[iDEnew];
    DEorder eDEnext = eDE;
    if (eDE.side == 0)
      eDEnext.dir=NextIdx(3,eDE.dir);
    else
      eDEnext.dir=PrevIdx(3,eDE.dir);
    DEorder eDEinv;
    eDEinv.iDE = -400; // just to avoid warning.
    eDEinv.dir = eDE.dir;
    eDEinv.side=1 - eDE.side;
    if (eDE.dir == 0)
      eDEinv.iDE = eDE.iDE;
    if (eDE.dir == 1)
      eDEinv.iDE = PL.invers.at(eDE.iDE);
    if (eDE.dir == 2) {
      if (eDE.side == 0)
	eDEinv.iDE = eNext.at(eDE.iDE);
      else
	eDEinv.iDE = ePrev.at(eDE.iDE);
    }
    int posNext=GetPosition(eDEnext);
    int posInv=GetPosition(eDEinv);
    NewNext[iDEnew]=posNext;
    NewInv [iDEnew]=posInv;
  }
  permlib::Permutation NewPermNext(NewNext);
  permlib::Permutation NewPermInv (NewInv);
  PlanGraphOriented PLnew{nbDEnew, NewPermInv, NewPermNext};
  return {PLnew, ListDE, OldToNew, ListDir, ListSide, SmallMat};
}


struct PairVertexInsert {
  PlanGraphOriented PL;
  std::vector<DupliEdge> ListDE;
  MyMatrix<int> OldToNew;
  std::vector<int> ListSide;
};


PairVertexInsert InsertVertexAllEdges(PlanGraphOriented const& PL)
{
  int nbP=PL.nbP;
  int nbPtotal=2*nbP;
  std::vector<DupliEdge> ListDE(nbPtotal);
  std::vector<int> ListSide(nbPtotal);
  int iDEnew=0;
  MyMatrix<int> OldToNew(nbP,2);
  for (int iDE=0; iDE<nbP; iDE++)
    for (int side=0; side<2; side++) {
      ListDE[iDEnew]={iDE, side};
      ListSide[iDEnew]=side;
      OldToNew(iDE,side)=iDEnew;
      iDEnew++;
    }
  int nbDEnew=ListDE.size();
  auto GetPosition=[&](DupliEdge const& eDE) -> int {
    for (int iPos=0; iPos<nbDEnew; iPos++)
      if (ListDE[iPos] == eDE)
	return iPos;
    std::cerr << "Did not find the edge\n";
    throw TerminalException{1};
    //    exit(1);
  };
  std::vector<permlib::dom_int> NewNext(nbDEnew);
  std::vector<permlib::dom_int> NewInv (nbDEnew);
  for (int iDEnew=0; iDEnew<nbDEnew; iDEnew++) {
    DupliEdge eDE=ListDE[iDEnew];
    DupliEdge eDEnext, eDEinv;
    eDEinv.iDE=eDE.iDE;
    eDEinv.side=1 - eDE.side;
    //
    eDEnext.side=eDE.side;
    if (eDE.side == 0)
      eDEnext.iDE=PL.next.at(eDE.iDE);
    else
      eDEnext.iDE=PL.invers.at(eDE.iDE);
    int posNext=GetPosition(eDEnext);
    int posInv=GetPosition(eDEinv);
    NewNext[iDEnew]=posNext;
    NewInv [iDEnew]=posInv;
  }
  permlib::Permutation NewPermNext(NewNext);
  permlib::Permutation NewPermInv (NewInv);
  PlanGraphOriented PLnew{nbDEnew, NewPermInv, NewPermNext};
  return {PLnew, ListDE, OldToNew, ListSide};
}


struct InfoDualMedialGraph {
  PlanGraphOriented PL;
  std::vector<DualMedialDE> ListDE;
  MyMatrix<int> MappingDE;
  MyMatrix<int> DecomposeOldDE;
};


InfoDualMedialGraph DualMedialGraphOriented(PlanGraphOriented const& PL)
{
  permlib::Permutation eNext=PL.next;
  permlib::Permutation ePrev=~eNext;
  permlib::Permutation eInv =PL.invers;
  int nbP=PL.nbP;
  auto DualMedialDE_to_int=[&](DualMedialDE const& eDE) -> int {
    if (eDE.stat == 1)
      return eDE.eDE1;
    if (eDE.stat == -1)
      return nbP + eDE.eDE1;
    std::cerr << "Shoud not reach that stage\n";
    throw TerminalException{1};
  };
  auto int_to_DualMediaDE=[&](int const& eVal) -> DualMedialDE {
    if (eVal < nbP) {
      int eDE2=eNext.at(eVal);
      return {eVal, eDE2, 1};
    }
    int eDE1=eVal - nbP;
    int eDE2=eNext.at(eDE1);
    return {eDE1, eDE2, -1};
  };
  int NewNbP=2*nbP;
  MyMatrix<int> MappingDE(nbP,2);
  for (int iP=0; iP<nbP; iP++) {
    int iPnew=iP;
    for (int i=0; i<2; i++) {
      if (i == 1)
	iPnew += nbP;
      MappingDE(iP,i)=iPnew;
    }
  }
  std::vector<DualMedialDE> NewListDE(NewNbP);
  for (int i=0; i<NewNbP; i++)
    NewListDE[i]=int_to_DualMediaDE(i);
  std::vector<permlib::dom_int> eListNext(NewNbP);
  std::vector<permlib::dom_int> eListInv (NewNbP);
  for (int iP=0; iP<NewNbP; iP++) {
    DualMedialDE eDE=NewListDE[iP];
    int eDE1=eDE.eDE1;
    int eDE2=eDE.eDE2;
    int stat=eDE.stat;
    DualMedialDE nextDE, invDE;
    if (stat == 1) {
      int eDE2next=eNext.at(eDE2);
      nextDE={eDE2, eDE2next, 1};
      invDE={eDE1, eDE2, -1};
    }
    else {
      int nextDE2=eInv.at(eDE1);
      int nextDE1=ePrev.at(nextDE2);
      nextDE={nextDE1, nextDE2, -1};
      invDE={eDE1, eDE2, 1};
    }
    eListNext[iP]=permlib::dom_int(DualMedialDE_to_int(nextDE));
    eListInv [iP]=permlib::dom_int(DualMedialDE_to_int(invDE));
  }
  permlib::Permutation NewPermNext(eListNext);
  permlib::Permutation NewPermInv (eListInv);
  PlanGraphOriented PLnew{NewNbP, NewPermInv, NewPermNext};
  MyMatrix<int> DecomposeOldDE(nbP,2);
  for (int iDE=0; iDE<nbP; iDE++) {
    int iDE2=eNext.at(iDE);
    DualMedialDE DE1{iDE,iDE2,1};
    //
    int rDE=eInv.at(iDE);
    int rDEprev=ePrev.at(rDE);
    DualMedialDE DE2{rDEprev, rDE,-1};
    DecomposeOldDE(iDE,0)=DualMedialDE_to_int(DE1);
    DecomposeOldDE(iDE,1)=DualMedialDE_to_int(DE2);
  }
  return {PLnew, NewListDE, MappingDE, DecomposeOldDE};
}


struct InfoVertFaceGraph {
  PlanGraphOriented PL;
  std::vector<DEvertface> ListDE;
  MyMatrix<int> OldToNew;
};


InfoVertFaceGraph VertFaceGraphOriented(PlanGraphOriented const& PL)
{
  permlib::Permutation eNext=PL.next;
  permlib::Permutation ePrev=~eNext;
  permlib::Permutation eInv =PL.invers;
  int nbP=PL.nbP;
  auto DEvertface_to_int=[&](DEvertface const& eDE) -> int {
    return eDE.iDE + nbP*eDE.nat;
  };
  auto int_to_DEvertface=[&](int const& eVal) -> DEvertface {
    int iDE=eVal % nbP;
    int nat=(eVal - iDE)/nbP;
    return {iDE, nat};
  };
  int NewNbP=3*nbP;
  MyMatrix<int> OldToNew(nbP,3);
  for (int iP=0; iP<nbP; iP++)
    for (int i=0; i<3; i++) {
      int iPnew =iP + nbP*i;
      OldToNew(iP,i)=iPnew;
    }
  std::vector<DEvertface> NewListDE(NewNbP);
  for (int i=0; i<NewNbP; i++)
    NewListDE[i]=int_to_DEvertface(i);
  std::vector<permlib::dom_int> eListNext(NewNbP);
  std::vector<permlib::dom_int> eListInv (NewNbP);
  for (int iP=0; iP<NewNbP; iP++) {
    DEvertface eDE=NewListDE[iP];
    int iDE=eDE.iDE;
    int nat=eDE.nat;
    DEvertface nextDE{-1,-1}, invDE{-1,-1};
    if (nat == 0) {
      nextDE={iDE,1};
      int rDE=eInv.at(iDE);
      invDE={rDE,0};
    }
    if (nat == 1) {
      int nDE=eNext.at(iDE);
      nextDE={nDE,0};
      invDE={iDE,2};
    }
    if (nat == 2) {
      int iDE1=eInv.at(iDE);
      int iDE2=ePrev.at(iDE1);
      nextDE={iDE2,2};
      invDE={iDE,1};
    }
    eListNext[iP]=permlib::dom_int(DEvertface_to_int(nextDE));
    eListInv [iP]=permlib::dom_int(DEvertface_to_int(invDE));
  }
  permlib::Permutation NewPermNext(eListNext);
  permlib::Permutation NewPermInv (eListInv);
  PlanGraphOriented PLnew{NewNbP, NewPermInv, NewPermNext};
  return {PLnew, NewListDE, OldToNew};
}


GraphSparseImmutable PlanGraphOrientedToGSI(VEForiented const& VEFori)
{
  int nbVert=VEFori.nbVert;
  std::vector<int> ListDeg(nbVert);
  int nbAdj=0;
  for (int iVert=0; iVert<nbVert; iVert++) {
    int siz=VEFori.VertSet[iVert].size();
    nbAdj += siz;
    ListDeg[iVert]=siz;
  }
  std::vector<int> ListStart(nbVert+1);
  ListStart[0]=0;
  for (int iVert=0; iVert<nbVert; iVert++)
    ListStart[iVert+1]=ListStart[iVert] + ListDeg[iVert];
  std::vector<int> ListListAdj(nbAdj);
  int idx=0;
  for (int iVert=0; iVert<nbVert; iVert++) {
    int siz=VEFori.VertSet[iVert].size();
    for (int i=0; i<siz; i++) {
      int iDE=VEFori.VertSet[iVert][i];
      int iDErev=VEFori.PLori.invers.at(iDE);
      int eVert=VEFori.ListOriginVert[iDErev];
      ListListAdj[idx]=eVert;
      idx++;
    }
  }
  return GraphSparseImmutable(nbVert, ListStart, ListListAdj);
}


std::vector<int> GetCycleAdjacent(VEForiented const& VEFori, std::vector<int> const& ListVertRemove)
{
  permlib::Permutation eNext=VEFori.PLori.next;
  permlib::Permutation ePrev=~eNext;
  permlib::Permutation eInv=VEFori.PLori.invers;
  permlib::Permutation FaceNext=eInv*ePrev;
  GraphSparseImmutable eG=PlanGraphOrientedToGSI(VEFori);
  std::vector<int> TotalListDEoriginating;
  for (int const& eVert : ListVertRemove) {
    for (int const& iDE : VEFori.VertSet[eVert]) {
      int rDE=eInv.at(iDE);
      int adjVert=VEFori.ListOriginVert[rDE];
      if (PositionVect(ListVertRemove, adjVert) == -1)
	TotalListDEoriginating.push_back(iDE);
    }
  }
  int siz=TotalListDEoriginating.size();
  if (siz == 0) {
    std::cerr << "We should never reach that stage\n";
    throw TerminalException{1};
  };
  int iDEinitial=TotalListDEoriginating[0];
  int iDEwork=iDEinitial;
  std::vector<int> ListDEinitial;
  std::vector<int> TheCycleVert;
  auto InsertSequence=[&](int const& iDEstart) -> int {
    ListDEinitial.push_back(iDEstart);
    int iDE=iDEstart;
    while(1) {
      int iDEnew=FaceNext.at(iDE);
      int rDE=eInv.at(iDEnew);
      int eVert=VEFori.ListOriginVert[iDEnew];
      TheCycleVert.push_back(eVert);
      int rVert=VEFori.ListOriginVert[rDE];
      if (PositionVect(ListVertRemove, rVert) != -1)
	return rDE;
      iDE=iDEnew;
    }
  };
  while(1) {
    iDEwork=InsertSequence(iDEwork);
    if (iDEwork == iDEinitial)
      break;
  }
  std::vector<int> SetInt1=VectorAsSet(TotalListDEoriginating);
  std::vector<int> SetInt2=VectorAsSet(ListDEinitial);
  if (SetInt1 != SetInt2) {
    std::cerr << "--------------------------------------\n";
    std::cerr << " ERROR in ComputGetCycleAdjacent\n";
    std::cerr << "Compute GetCycleAdjacent |ListVertRemove|=" << ListVertRemove.size() << "\n";
    for (auto & eVert : ListVertRemove) {
      int eDeg=VEFori.VertSet[eVert].size();
      std::cerr << " [" << eVert << "," << eDeg << "]";
    }
    std::cerr << "\n";
    //
    PrintInformationOfMap(std::cerr, VEFori.PLori);
    //
    std::cerr << "ListVertRemove=";
    WriteStdVectorGAP(std::cerr, ListVertRemove);
    std::cerr << "\n";
    //
    std::cerr << "SetInt1=";
    WriteStdVectorGAP(std::cerr, SetInt1);
    std::cerr << "\n";
    //
    std::cerr << "SetInt2=";
    WriteStdVectorGAP(std::cerr, SetInt2);
    std::cerr << "\n";
    //
    std::cerr << "|SetInt1|=" << SetInt1.size() << " |SetInt2|=" << SetInt2.size() << "\n";
    //
    std::cerr << "SetInt1 != SetInt2 which is contrary to what we expect\n";
    throw TerminalException{1};
  }
  return TheCycleVert;
}


#endif
