#ifndef INCLUDE_PlaneGraphDrawing_h
#define INCLUDE_PlaneGraphDrawing_h


#include "Namelist.h"
#include "MapCoordinateFinder.h"


struct ExternalFaceInfo {
  std::string eNature;
  std::vector<int> ListDE;
};


struct ArrowStructuralInfo {
  std::string TypeDir;
  int factormult;
  std::vector<int> ListDE;
};


FullNamelist NAMELIST_GetStandardCombPlaneGraph()
{
  std::map<std::string, SingleBlock> ListBlock;
  // PROC
  std::map<std::string, int> ListIntValues1;
  std::map<std::string, bool> ListBoolValues1;
  std::map<std::string, double> ListDoubleValues1;
  std::map<std::string, std::string> ListStringValues1;
  std::map<std::string, std::vector<std::string> > ListListStringValues1;
  std::map<std::string, std::vector<int> > ListListIntValues1;

  ListStringValues1["PlaneFile"]="PLori";
  ListStringValues1["OutFile"]="unset.svg";
  ListStringValues1["ViewFile"]="unset_view";
  ListIntValues1["MAX_ITER_PrimalDual"]=1000;
  ListIntValues1["MAX_ITER_CaGe"]=1000;
  ListIntValues1["CaGeProcessPolicy"]=0;
  ListDoubleValues1["MaximalFractionArea"]=100;
  ListIntValues1["RoundMethod"]=2;
  ListDoubleValues1["width"]=600;
  ListDoubleValues1["height"]=600;
  ListIntValues1["MethodInsert"]=2;
  ListDoubleValues1["exponent"]=double(2);
  ListDoubleValues1["SphereRotationRad"]=0;
  ListListStringValues1["ListExportFormat"]={};
  SingleBlock BlockPLOT;
  BlockPLOT.ListIntValues=ListIntValues1;
  BlockPLOT.ListBoolValues=ListBoolValues1;
  BlockPLOT.ListDoubleValues=ListDoubleValues1;
  BlockPLOT.ListStringValues=ListStringValues1;
  BlockPLOT.ListListStringValues=ListListStringValues1;
  BlockPLOT.ListListIntValues=ListListIntValues1;
  ListBlock["PLOT"]=BlockPLOT;
  // EDGE information
  std::map<std::string, int> ListIntValues2;
  std::map<std::string, bool> ListBoolValues2;
  std::map<std::string, double> ListDoubleValues2;
  std::map<std::string, std::string> ListStringValues2;
  std::map<std::string, std::vector<std::string> > ListListStringValues2;
  std::map<std::string, std::vector<int> > ListListIntValues2;
  ListBoolValues2["DoMethod1"]=false;
  ListBoolValues2["DoMethod2"]=false;
  ListBoolValues2["DoMethod3"]=false;
  ListDoubleValues2["MultTangent"]=0.5;
  ListIntValues2["NormalTraitSize"]=1;
  ListListIntValues2["ListTraitIDE"]={};
  ListListIntValues2["ListTraitGroup"]={};
  ListListIntValues2["ListTraitSize"]={};
  ListListIntValues2["DefaultRGB"]={0,0,0};
  ListListIntValues2["SpecificRGB_iDE"]={};
  ListListIntValues2["SpecificRGB_Group"]={};
  ListListIntValues2["SpecificRGB_R"]={};
  ListListIntValues2["SpecificRGB_G"]={};
  ListListIntValues2["SpecificRGB_B"]={};
  ListBoolValues2["DoExternalArrow"]=true;
  ListDoubleValues2["ScalExtensionArrow"]=2;
  ListDoubleValues2["LengthExtensionArrow"]=0.5;
  SingleBlock BlockEDGE;
  BlockEDGE.ListIntValues=ListIntValues2;
  BlockEDGE.ListBoolValues=ListBoolValues2;
  BlockEDGE.ListDoubleValues=ListDoubleValues2;
  BlockEDGE.ListStringValues=ListStringValues2;
  BlockEDGE.ListListStringValues=ListListStringValues2;
  BlockEDGE.ListListIntValues=ListListIntValues2;
  ListBlock["EDGE"]=BlockEDGE;
  // VERT information
  std::map<std::string, int> ListIntValues3;
  std::map<std::string, bool> ListBoolValues3;
  std::map<std::string, double> ListDoubleValues3;
  std::map<std::string, std::string> ListStringValues3;
  std::map<std::string, std::vector<std::string> > ListListStringValues3;
  std::map<std::string, std::vector<int> > ListListIntValues3;
  std::map<std::string, std::vector<double> > ListListDoubleValues3;
  ListDoubleValues3["NormalRadius"]=1;
  ListListIntValues3["ListRadiusIDE"]={};
  ListListIntValues3["ListRadiusGroup"]={};
  ListListDoubleValues3["ListRadius"]={};
  ListListIntValues3["DefaultRGB"]={0,0,0};
  ListListIntValues3["SpecificRGB_iDE"]={};
  ListListIntValues3["SpecificRGB_Group"]={};
  ListListIntValues3["SpecificRGB_R"]={};
  ListListIntValues3["SpecificRGB_G"]={};
  ListListIntValues3["SpecificRGB_B"]={};
  SingleBlock BlockVERT;
  BlockVERT.ListIntValues=ListIntValues3;
  BlockVERT.ListBoolValues=ListBoolValues3;
  BlockVERT.ListDoubleValues=ListDoubleValues3;
  BlockVERT.ListStringValues=ListStringValues3;
  BlockVERT.ListListStringValues=ListListStringValues3;
  BlockVERT.ListListIntValues=ListListIntValues3;
  BlockVERT.ListListDoubleValues=ListListDoubleValues3;
  ListBlock["VERT"]=BlockVERT;
  // TORUS information
  std::map<std::string, int> ListIntValues4;
  std::map<std::string, bool> ListBoolValues4;
  std::map<std::string, double> ListDoubleValues4;
  std::map<std::string, std::string> ListStringValues4;
  std::map<std::string, std::vector<std::string> > ListListStringValues4;
  std::map<std::string, std::vector<int> > ListListIntValues4;
  ListDoubleValues4["minimal"]=10^(-8);
  ListDoubleValues4["tol"]=10^(-5);
  ListDoubleValues4["AngDeg"]=0;
  ListDoubleValues4["scal"]=1;
  ListDoubleValues4["shiftX"]=0;
  ListDoubleValues4["shiftY"]=0;
  ListListIntValues4["FundamentalRGB"]={255,0,0};
  ListIntValues4["FundamentalTraitSize"]=1;
  ListBoolValues4["DrawFundamentalDomain"]=false;
  SingleBlock BlockTORUS;
  BlockTORUS.ListIntValues=ListIntValues4;
  BlockTORUS.ListBoolValues=ListBoolValues4;
  BlockTORUS.ListDoubleValues=ListDoubleValues4;
  BlockTORUS.ListStringValues=ListStringValues4;
  BlockTORUS.ListListStringValues=ListListStringValues4;
  BlockTORUS.ListListIntValues=ListListIntValues4;
  ListBlock["TORUS"]=BlockTORUS;
  return {ListBlock, "undefined"};
}


ExternalFaceInfo PLANE_ReadViewFile(std::string const& FileName)
{
  std::ifstream isView(FileName);
  std::string eNature;
  isView >> eNature;
  std::vector<int> ListDE=ReadStdVector<int>(isView);
  return {eNature, ListDE};
}

template<typename Telt>
struct CombinatorialToolPL {
  VEForiented<Telt> VEFori;
  std::function<std::vector<int>(int const&)> MapNativeDE;
  std::function<std::vector<int>(int const&)> GetListVert;
  std::function<std::vector<int>(ExternalFaceInfo const&)> GetListVertRemove;
  std::function<ArrowStructuralInfo(int const&)> GetArrowInfo;
  std::function<int(int const&)> MapSingleDE;
};


template<typename Telt>
CombinatorialToolPL<Telt> COMBI_PL_GeneralCase(PlanGraphOriented<Telt> const& PL)
{
  auto ePair1=GetWythoff123(PL);
  auto ePair2=InsertVertexAllEdges(ePair1.PL);
  auto ePair3=VertFaceGraphOriented(ePair2.PL);
  VEForiented VEFori=PlanGraphOrientedToVEF(ePair3.PL);
  GraphSparseImmutable GR=PlanGraphOrientedToGSI(VEFori);
  std::cerr << "PL=\n";
  PrintInformationOfMap(std::cerr, PL);
  std::cerr << "ePair1.PL=\n";
  PrintInformationOfMap(std::cerr, ePair1.PL);
  std::cerr << "ePair2.PL=\n";
  PrintInformationOfMap(std::cerr, ePair2.PL);
  std::cerr << "ePair3.PL=\n";
  PrintInformationOfMap(std::cerr, ePair3.PL);

  int pos1_00=ePair1.SmallMat(0,0);
  int pos1_10=ePair1.SmallMat(1,0);
  int pos1_01=ePair1.SmallMat(0,1);
  int pos1_02=ePair1.SmallMat(0,2);
  int pos1_11=ePair1.SmallMat(1,1);

  std::function<std::vector<int>(int const&)> MapNativeDE=[=](int const& iDE0) -> std::vector<int> {
    int rDE0=PL.invers.at(iDE0);
    int iDE1_a=ePair1.OldToNew(iDE0, pos1_00);
    int iDE1_b=ePair1.OldToNew(iDE0, pos1_10);
    int iDE1_c=ePair1.OldToNew(rDE0, pos1_10);
    int iDE1_d=ePair1.OldToNew(rDE0, pos1_00);
    int iDE2_a=ePair2.OldToNew(iDE1_a, 1);
    int iDE2_b=ePair2.OldToNew(iDE1_b, 1);
    int iDE2_c=ePair2.OldToNew(iDE1_c, 1);
    int iDE2_d=ePair2.OldToNew(iDE1_d, 1);
    int iDE3_a=ePair3.OldToNew(iDE2_a, 2);
    int iDE3_b=ePair3.OldToNew(iDE2_b, 1);
    int iDE3_c=ePair3.OldToNew(iDE2_c, 2);
    int iDE3_d=ePair3.OldToNew(iDE2_d, 1);
    return {iDE3_a, iDE3_b, iDE3_c, iDE3_d};
  };
  std::function<std::vector<int>(int const&)> GetListVert=[=](int const& iDE) -> std::vector<int> {
    std::vector<int> ListVert;
    std::vector<int> ListDE=MapNativeDE(iDE);
    for (auto & iDEmap : ListDE) {
      int iVert=VEFori.ListOriginVert[iDEmap];
      ListVert.push_back(iVert);
    }
    int siz=ListDE.size();
    int rDE=ePair3.PL.invers.at(ListDE[siz-1]);
    int jVert=VEFori.ListOriginVert[rDE];
    ListVert.push_back(jVert);
    /*
    for (int i=0; i<=siz; i++) {
      int iVert=ListVert[i];
      int nbAdj=VEFori.VertSet[iVert].size();
      std::cerr << "iVert=" << iVert << " nbAdj=" << nbAdj << "\n";
    }
    */
    for (int i=1; i<=siz; i++) {
      int iVert=ListVert[i-1];
      int jVert=ListVert[i];
      if (GR.IsAdjacent(iVert, jVert) == false) {
	std::cerr << "i=" << i << "\n";
	std::cerr << "The points are not adjacent\n";
	std::cerr << "Big error by any mean\n";
	throw TerminalException{1};
      }
    }
    return ListVert;
  };
  std::function<std::vector<int>(ExternalFaceInfo const&)> GetListVertRemove=[=](ExternalFaceInfo const& eExt) -> std::vector<int> {
    std::set<int> SetVertRemove;
    if (eExt.eNature == "vert")
      for (auto & iDE0 : eExt.ListDE) {
	int iDE1_a=ePair1.OldToNew(iDE0,pos1_00);
	int iDE1_b=ePair1.OldToNew(iDE0,pos1_01);
	int iDE1_c=ePair1.OldToNew(iDE0,pos1_11);
	int iDE1_d=ePair1.OldToNew(iDE0,pos1_02);
	int iDE2_a1=ePair2.OldToNew(iDE1_a,0);
	int iDE2_a2=ePair2.OldToNew(iDE1_a,1);
	int iDE2_b1=ePair2.OldToNew(iDE1_b,0);
	int iDE2_b2=ePair2.OldToNew(iDE1_b,1);
	int iDE2_c1=ePair2.OldToNew(iDE1_c,0);
	int iDE2_c2=ePair2.OldToNew(iDE1_c,1);
	int iDE2_d1=ePair2.OldToNew(iDE1_d,1);
	int iDE3_1=ePair3.OldToNew(iDE2_a2,2); // the central vertex corresponding to the removed vertex v
	int iDE3_2=ePair3.OldToNew(iDE2_a2,1); // vertex just opposite to v (class 2)
	int iDE3_3=ePair3.OldToNew(iDE2_a1,0); // vertices between the class 2 vertices
	int iDE3_4=ePair3.OldToNew(iDE2_b1,2); // the vertex corresponding to the face contained in v
	int iDE3_5=ePair3.OldToNew(iDE2_b2,0); // the vertices between two class 3 vertices
	int iDE3_6=ePair3.OldToNew(iDE2_b2,2); // at the center of the face corresponding to the edge
	int iDE3_7=ePair3.OldToNew(iDE2_c2,0); // On the other side of the class 5 vertices
	int iDE3_8=ePair3.OldToNew(iDE2_c1,0); // Similar to C3 but on the other side of v [2v-C3]
	int iDE3_9=ePair3.OldToNew(iDE2_d1,0); // Similar to class 2 but for edges.
	int iVert1=VEFori.ListOriginVert[iDE3_1];
	int iVert2=VEFori.ListOriginVert[iDE3_2];
	int iVert3=VEFori.ListOriginVert[iDE3_3];
	int iVert4=VEFori.ListOriginVert[iDE3_4];
	int iVert5=VEFori.ListOriginVert[iDE3_5];
	int iVert6=VEFori.ListOriginVert[iDE3_6];
	int iVert7=VEFori.ListOriginVert[iDE3_7];
	int iVert8=VEFori.ListOriginVert[iDE3_8];
	int iVert9=VEFori.ListOriginVert[iDE3_9];
	/*
	std::cerr << "------------------------------------------------------------------\n";
	std::cerr << " iVert1=" << iVert1 << " deg=" << VEFori.VertSet[iVert1].size() << "\n";
	std::cerr << " iVert2=" << iVert2 << " deg=" << VEFori.VertSet[iVert2].size() << "\n";
	std::cerr << " iVert3=" << iVert3 << " deg=" << VEFori.VertSet[iVert3].size() << "\n";
	std::cerr << " iVert4=" << iVert4 << " deg=" << VEFori.VertSet[iVert4].size() << "\n";
	std::cerr << " iVert5=" << iVert5 << " deg=" << VEFori.VertSet[iVert5].size() << "\n";
	std::cerr << " iVert6=" << iVert6 << " deg=" << VEFori.VertSet[iVert6].size() << "\n";
	std::cerr << " iVert7=" << iVert7 << " deg=" << VEFori.VertSet[iVert7].size() << "\n";
	std::cerr << " iVert8=" << iVert8 << " deg=" << VEFori.VertSet[iVert8].size() << "\n";
	std::cerr << " iVert9=" << iVert9 << " deg=" << VEFori.VertSet[iVert9].size() << "\n";*/
	SetVertRemove.insert(iVert1);
	SetVertRemove.insert(iVert2);
	SetVertRemove.insert(iVert3);
	SetVertRemove.insert(iVert4);
	SetVertRemove.insert(iVert5);
	SetVertRemove.insert(iVert6);
	SetVertRemove.insert(iVert7);
	SetVertRemove.insert(iVert8);
	SetVertRemove.insert(iVert9);
      }
    if (eExt.eNature == "edge")
      for (auto & iDE0 : eExt.ListDE) {
	int iDE1_a=ePair1.OldToNew(iDE0, pos1_01);
	int iDE2_a1=ePair2.OldToNew(iDE1_a,0);
	int iDE2_a2=ePair2.OldToNew(iDE1_a,1);
	int iDE3_1=ePair3.OldToNew(iDE2_a1,2); // center of the face (adjacent to class 2)
	int iDE3_2=ePair3.OldToNew(iDE2_a2,0); // middle of the edge of Wythoff123(G)
	int iDE3_3=ePair3.OldToNew(iDE2_a2,2); // middle of the edge (only 1 vertex, adjacent to class 2)
	int iVert1=VEFori.ListOriginVert[iDE3_1];
	int iVert2=VEFori.ListOriginVert[iDE3_2];
	int iVert3=VEFori.ListOriginVert[iDE3_3];
	/*	std::cerr << "--------------------------------------------\n";
	std::cerr << " iVert1=" << iVert1 << " deg=" << VEFori.VertSet[iVert1].size() << "\n";
	std::cerr << " iVert2=" << iVert2 << " deg=" << VEFori.VertSet[iVert2].size() << "\n";
	std::cerr << " iVert3=" << iVert3 << " deg=" << VEFori.VertSet[iVert3].size() << "\n";*/
	SetVertRemove.insert(iVert1);
	SetVertRemove.insert(iVert2);
	SetVertRemove.insert(iVert3);
      }
    if (eExt.eNature == "face")
      for (auto & iDE : eExt.ListDE) {
	int iDE1=ePair1.OldToNew(iDE, pos1_01);
	int iDE2=ePair2.OldToNew(iDE1, 0);
	int iDE3=ePair3.OldToNew(iDE2,2);
	int iVert1=VEFori.ListOriginVert[iDE3];
	SetVertRemove.insert(iVert1);
      }
    int siz=SetVertRemove.size();
    //    std::cerr << "siz=" << siz << "\n";
    if (siz == 0) {
      std::cerr << "We should have at least one vertex to remove for this combinatorial options\n";
      throw TerminalException{1};
    }
    std::vector<int> ListVertRemove;
    for (auto & eVert : SetVertRemove)
      ListVertRemove.push_back(eVert);
    return ListVertRemove;
  };
  std::function<ArrowStructuralInfo(int const&)> GetArrowInfo=[=](int const& iDE0) -> ArrowStructuralInfo {
    int iDE1=ePair1.OldToNew(iDE0, pos1_00);
    int iDE2=ePair2.OldToNew(iDE1, 1);
    int iDE3=ePair3.OldToNew(iDE2, 2);
    return {"direct", 1, {iDE3}};
  };
  std::function<int(int const&)> MapSingleDE=[=](int const& iDE0) -> int {
    int iDE1=ePair1.OldToNew(iDE0, pos1_00);
    int iDE2=ePair2.OldToNew(iDE1, 1);
    int iDE3=ePair3.OldToNew(iDE2, 2);
    return iDE3;
  };
  return {VEFori, MapNativeDE, GetListVert, GetListVertRemove, GetArrowInfo, MapSingleDE};
}


template<typename Telt>
void CreateSVGfileOfGraph(FullNamelist const& eFull)
{
  //
  // General : Preparation work of data reading
  //
  SingleBlock BlPLOT=eFull.ListBlock.at("PLOT");
  int MAX_ITER_PrimalDual=BlPLOT.ListIntValues.at("MAX_ITER_PrimalDual");
  int MAX_ITER_CaGe=BlPLOT.ListIntValues.at("MAX_ITER_CaGe");
  int CaGeProcessPolicy=BlPLOT.ListIntValues.at("CaGeProcessPolicy");
  double MaximalFractionArea=BlPLOT.ListDoubleValues.at("MaximalFractionArea");
  std::string OutFile=BlPLOT.ListStringValues.at("OutFile");
  std::string PlaneFile=BlPLOT.ListStringValues.at("PlaneFile");
  int RoundMethod=BlPLOT.ListIntValues.at("RoundMethod");
  int MethodInsert=BlPLOT.ListIntValues.at("MethodInsert");
  double height=BlPLOT.ListDoubleValues.at("height");
  double width=BlPLOT.ListDoubleValues.at("width");
  double exponent=BlPLOT.ListDoubleValues.at("exponent");
  double SphereRotationRad=BlPLOT.ListDoubleValues.at("SphereRotationRad");
  PlanGraphOriented PL=ReadPlanGraphOrientedFile<Telt>(PlaneFile);
  std::vector<std::string> ListExportFormat=BlPLOT.ListListStringValues.at("ListExportFormat");
  // EDGE
  SingleBlock BlEDGE=eFull.ListBlock.at("EDGE");
  bool DoMethod1=BlEDGE.ListBoolValues.at("DoMethod1");
  bool DoMethod2=BlEDGE.ListBoolValues.at("DoMethod2");
  bool DoMethod3=BlEDGE.ListBoolValues.at("DoMethod3");
  double MultTangent=BlEDGE.ListDoubleValues.at("MultTangent");
  bool DoExternalArrow=BlEDGE.ListBoolValues.at("DoExternalArrow");
  double ScalExtensionArrow=BlEDGE.ListDoubleValues.at("ScalExtensionArrow");
  double LengthExtensionArrow=BlEDGE.ListDoubleValues.at("LengthExtensionArrow");
  std::vector<int> EDGE_ListTraitIDE=BlEDGE.ListListIntValues.at("ListTraitIDE");
  std::vector<int> EDGE_ListTraitGroup=BlEDGE.ListListIntValues.at("ListTraitGroup");
  std::vector<int> EDGE_ListTraitSize=BlEDGE.ListListIntValues.at("ListTraitSize");
  std::vector<int> EDGE_DefaultRGB=BlEDGE.ListListIntValues.at("DefaultRGB");
  std::vector<int> EDGE_SpecificRGB_iDE=BlEDGE.ListListIntValues.at("SpecificRGB_iDE");
  std::vector<int> EDGE_SpecificRGB_Group=BlEDGE.ListListIntValues.at("SpecificRGB_Group");
  std::vector<int> EDGE_SpecificRGB_R=BlEDGE.ListListIntValues.at("SpecificRGB_R");
  std::vector<int> EDGE_SpecificRGB_G=BlEDGE.ListListIntValues.at("SpecificRGB_G");
  std::vector<int> EDGE_SpecificRGB_B=BlEDGE.ListListIntValues.at("SpecificRGB_B");
  int EDGE_NormalTraitSize=BlEDGE.ListIntValues.at("NormalTraitSize");
  // VERT
  SingleBlock BlVERT=eFull.ListBlock.at("VERT");
  std::vector<int> VERT_ListRadiusIDE=BlVERT.ListListIntValues.at("ListRadiusIDE");
  std::vector<int> VERT_ListRadiusGroup=BlVERT.ListListIntValues.at("ListRadiusGroup");
  std::vector<double> VERT_ListRadius=BlVERT.ListListDoubleValues.at("ListRadius");
  std::vector<int> VERT_DefaultRGB=BlVERT.ListListIntValues.at("DefaultRGB");
  std::vector<int> VERT_SpecificRGB_iDE=BlVERT.ListListIntValues.at("SpecificRGB_iDE");
  std::vector<int> VERT_SpecificRGB_Group=BlVERT.ListListIntValues.at("SpecificRGB_Group");
  std::vector<int> VERT_SpecificRGB_R=BlVERT.ListListIntValues.at("SpecificRGB_R");
  std::vector<int> VERT_SpecificRGB_G=BlVERT.ListListIntValues.at("SpecificRGB_G");
  std::vector<int> VERT_SpecificRGB_B=BlVERT.ListListIntValues.at("SpecificRGB_B");
  double VERT_NormalRadius=BlVERT.ListDoubleValues.at("NormalRadius");
  //
  // General : combinatorial work of creating the directed edges graphs
  //
  int Charac=EulerPoincareCharacteristic(PL);
  VEForiented VEFori=PlanGraphOrientedToVEF(PL);
  if (Charac != 0 && Charac != 2) {
    std::cerr << "The map is neither a plane nor a torus\n";
    std::cerr << "Maybe buggy map or a map of genus >1 that has not been programmed yet\n";
    throw TerminalException{1};
  }
  CombinatorialToolPL<Telt> eCombi=COMBI_PL_GeneralCase(PL);
  //
  // General: reading colors and other properties
  //
  int nbTraitIDE=EDGE_ListTraitIDE.size();
  int nbTraitGroup=EDGE_ListTraitGroup.size();
  if (nbTraitIDE != nbTraitGroup) {
    std::cerr << "We have nbTraitIDE = " << nbTraitIDE << "\n";
    std::cerr << "  and nbTraitGroup = " << nbTraitGroup << "\n";
    throw TerminalException{1};
  }
  auto Kernel_GetTraitSizeDE=[&](int const& iDE) -> int {
    for (int i=0; i<nbTraitIDE; i++)
      if (EDGE_ListTraitIDE[i] == iDE) {
	int iGroup=EDGE_ListTraitGroup[i];
	return EDGE_ListTraitSize[iGroup];
      }
    return EDGE_NormalTraitSize;
  };
  auto GetTraitSizeDE=[&](int const& iDE) -> int {
    int eSize=Kernel_GetTraitSizeDE(iDE);
    int rDE=PL.invers.at(iDE);
    int rSize=Kernel_GetTraitSizeDE(rDE);
    if (eSize != rSize) {
      std::cerr << "Error in Trait size\n";
      std::cerr << "iDE=" << iDE << " eSize=" << eSize << "\n";
      std::cerr << "rDE=" << rDE << " rSize=" << rSize << "\n";
      throw TerminalException{1};
    }
    return eSize;
  };
  auto GetColorDE=[&](int const& iDE) -> std::vector<int> {
    int len=EDGE_SpecificRGB_iDE.size();
    for (int i=0; i<len; i++)
      if (EDGE_SpecificRGB_iDE[i] == iDE) {
	int iGroup=EDGE_SpecificRGB_Group[i];
	return {EDGE_SpecificRGB_R[iGroup], EDGE_SpecificRGB_G[iGroup], EDGE_SpecificRGB_B[iGroup]};
      }
    return EDGE_DefaultRGB;
  };
  auto GetRadiusVert=[&](int const& iDE) -> double {
    int len=VERT_ListRadiusIDE.size();
    for (int i=0; i<len; i++)
      if (VERT_ListRadiusIDE[i] == iDE) {
	int iGroup=VERT_ListRadiusGroup[i];
	return VERT_ListRadius[iGroup];
      }
    //    std::cerr << "VERT_NormalRadius=" << VERT_NormalRadius << "\n";
    return VERT_NormalRadius;
  };
  auto GetColorVert=[&](int const& iDE) -> std::vector<int> {
    int len=VERT_SpecificRGB_iDE.size();
    for (int i=0; i<len; i++)
      if (VERT_SpecificRGB_iDE[i] == iDE) {
	int iGroup=VERT_SpecificRGB_Group[i];
	return {VERT_SpecificRGB_R[iGroup], VERT_SpecificRGB_G[iGroup], VERT_SpecificRGB_B[iGroup]};
      }
    return VERT_DefaultRGB;
  };
  //
  // General: other selections and miscelaneous properties
  //
  std::vector<int> ListSel;
  if (MethodInsert == 0)
    ListSel={0,1,2,3,4};
  if (MethodInsert == 1)
    ListSel={0,1,3,4};
  if (MethodInsert == 2)
    ListSel={0,2,4};
  if (MethodInsert == 3)
    ListSel={0,4};
  int lenSel=ListSel.size();
  std::vector<SVGclippath> ListClip;
  std::vector<SVGpolyline> ListPolyline;
  std::vector<SVGbezier> ListBezier;
  std::vector<SVGline> ListLine;
  std::vector<SVGellipse> ListEllipse;
  std::string clipname="";
  std::cerr << "DoMethod1 / 2 / 3=" << DoMethod1 << " " << DoMethod2 << " " << DoMethod3 << "\n";
  auto DrawLineOfPlot=[&](int const& iDE, std::vector<coor> const& ListPoint) -> void {
    std::vector<coor> ListTangent=GetListTangent(ListPoint, MultTangent);
    int eSize=GetTraitSizeDE(iDE);
    std::vector<int> color=GetColorDE(iDE);
    if (DoMethod1) {
      for (int j=0; j<lenSel-1; j++) {
	int i1=ListSel[j];
	int i2=ListSel[j+1];
	coor coorM=ListPoint[i1];
	coor coorC=ListPoint[i1] + ListTangent[i1];
	coor coor1=ListPoint[i2] - ListTangent[i2];
	coor coor2=ListPoint[i2];
	SVGbezier eBez{coorM, coorC, coor1, coor2, {color, eSize, "", clipname}};
	ListBezier.push_back(eBez);
      }
    }
    if (DoMethod2) {
      std::vector<coor> ListPointSel(lenSel);
      for (int j=0; j<lenSel; j++)
	ListPointSel[j]=ListPoint[ListSel[j]];
      ListPolyline.push_back({ListPointSel, {color, eSize, "", clipname}});
    }
    if (DoMethod3) {
      for (int j=0; j<lenSel-1; j++) {
	int i1=ListSel[j];
	int i2=ListSel[j+1];
	coor ePt=ListPoint[i1];
	coor fPt=ListPoint[i2];
	SVGline eLine{ePt, fPt, {color, eSize, "", clipname}};
	ListLine.push_back(eLine);
      }
    }
  };
  auto DrawVertexOfPlot=[&](int const& iDE, coor const& ePt) -> void {
    double eRadius=GetRadiusVert(iDE);
    if (eRadius > 0) {
      std::vector<int> color=GetColorVert(iDE);
      coor rad{eRadius, eRadius};
      SVGellipse eEll{ePt, rad, {color, -1, "", clipname}};
      ListEllipse.push_back(eEll);
    }
  };
  //
  // End of general part. Now we compute the specific to each map information
  //
  auto IsSegmentAdmissible=[&](std::vector<coor> const& ListPoint) -> bool {
    std::vector<double> ListX, ListY;
    for (auto & ePoint : ListPoint) {
      ListX.push_back(ePoint.x);
      ListY.push_back(ePoint.y);
    }
    double MinX=VectorMin(ListX);
    double MaxX=VectorMax(ListX);
    double MinY=VectorMin(ListY);
    double MaxY=VectorMax(ListY);
    if (MinX > 1)
      return false;
    if (MaxX < 0)
      return false;
    if (MinY > 1)
      return false;
    if (MaxY < 0)
      return false;
    return true;
  };
  auto IsVertexAdmissible=[&](coor const& ePt) -> bool {
    if (ePt.x > 1)
      return false;
    if (ePt.x < 0)
      return false;
    if (ePt.y > 1)
      return false;
    if (ePt.y < 0)
      return false;
    return true;
  };
  SVGplotDescription eSVGplot;
  eSVGplot.width=width;
  eSVGplot.height=height;
  eSVGplot.RoundMethod=RoundMethod;
  std::cerr << "exponent=" << exponent << "\n";
  std::function<double(double const&)> AreaToKoef=[&exponent](double const& area) -> double {
    double area_a=fabs(area);
    return pow(area_a, exponent);
  };
  if (Charac == 2) {
    std::string ViewFile=BlPLOT.ListStringValues.at("ViewFile");
    ExternalFaceInfo eExtInfo=PLANE_ReadViewFile(ViewFile);
    //
    // Determining vertex status and constructing the cycles
    //
    std::vector<int> ListKillVert = eCombi.GetListVertRemove(eExtInfo);
    std::vector<int> TheCycle=GetCycleAdjacent(eCombi.VEFori, ListKillVert);
    //
    // Compute coordinate of vertices by the iteration process
    //
    PLANE_plot_infos eDescPre1=PLANE_ComputeCoordinateVertices(eCombi.VEFori, ListKillVert, TheCycle, MAX_ITER_CaGe, AreaToKoef);
    PLANE_plot_infos eDesc=PLANE_RotationSphericalStructure(eDescPre1, SphereRotationRad);
    //
    // Creation of SVG information for the plotting
    //
    for (int iDE=0; iDE<PL.nbP; iDE++) {
      int iDErev=PL.invers.at(iDE);
      bool IsCorr=true;
      if (iDErev < iDE)
	IsCorr=false;
      if (eExtInfo.eNature == "vert" || eExtInfo.eNature == "edge") {
	int pos1=PositionVect(eExtInfo.ListDE, iDE);
	if (pos1 != -1)
	  IsCorr=false;
	int pos2=PositionVect(eExtInfo.ListDE, iDErev);
	if (pos2 != -1)
	  IsCorr=false;
      }
      if (IsCorr == true) {
	std::vector<int> ListVert=eCombi.GetListVert(iDE);
	std::vector<coor> ListPoint;
	for (auto & iVert : ListVert)
	  ListPoint.push_back(eDesc.ListCoord[iVert]);
	DrawLineOfPlot(iDE, ListPoint);
      }
    }
    int nbVert=VEFori.VertSet.size();
    for (int iVert=0; iVert<nbVert; iVert++) {
      int iDE=VEFori.VertSet[iVert][0];
      int iDEmap=eCombi.MapSingleDE(iDE);
      int iVertMap=eCombi.VEFori.ListOriginVert[iDEmap];
      coor ePt=eDesc.ListCoord[iVertMap];
      DrawVertexOfPlot(iDE, ePt);
    }
    if (DoExternalArrow) {
      int MethodArrow=2;
      if (eExtInfo.eNature == "vert" || eExtInfo.eNature == "edge") {
	for (auto & iDE : eExtInfo.ListDE) {
	  int rDE=PL.invers.at(iDE);
	  ArrowStructuralInfo eArrInfo=eCombi.GetArrowInfo(rDE);
	  coor vect{0,0};
	  for (auto & iDEmap : eArrInfo.ListDE) {
	    int iVert=eCombi.VEFori.ListOriginVert[iDEmap];
	    int rDE=eCombi.VEFori.PLori.invers.at(iDEmap);
	    int jVert=eCombi.VEFori.ListOriginVert[rDE];
	    if (PositionVect(ListKillVert, iVert) != -1 || PositionVect(ListKillVert, jVert) != -1) {
	      std::cerr << "Error, iVert and jVert should not be in the killed vertices\n";
	      throw TerminalException{1};
	    }
	    coor eDiff=eDesc.ListCoord[jVert] - eDesc.ListCoord[iVert];
	    vect = vect + eDiff;
	  }
	  int iDEmap1=eArrInfo.ListDE[0];
	  int iVert1=eCombi.VEFori.ListOriginVert[iDEmap1];
	  coor point1=eDesc.ListCoord[iVert1];
	  double eFact=ScalExtensionArrow*double(eArrInfo.factormult);
	  coor ePointFinal=point1 + MultScal(vect, eFact);
	  int eSize=GetTraitSizeDE(iDE);
	  std::vector<int> color=GetColorDE(iDE);
	  if (MethodArrow == 1) {
	    std::string MarkerEnd="triangle";
	    SVGline eLine{point1, ePointFinal, {color, eSize, MarkerEnd, clipname}};
	    ListLine.push_back(eLine);
	  }
	  if (MethodArrow == 2) {
	    coor eVectRot1=RotateCoor(vect, 135);
	    coor eVectRot2=RotateCoor(vect, -135);
	    double eFactB=LengthExtensionArrow*double(eArrInfo.factormult);
	    coor ePtFin1=ePointFinal + MultScal(eVectRot1, eFactB);
	    coor ePtFin2=ePointFinal + MultScal(eVectRot2, eFactB);
	    SVGline eLine {point1, ePointFinal, {color, eSize, "", clipname}};
	    SVGline eLine1{ePointFinal, ePtFin1, {color, eSize, "", clipname}};
	    SVGline eLine2{ePointFinal, ePtFin2, {color, eSize, "", clipname}};
	    ListLine.push_back(eLine);
	    ListLine.push_back(eLine1);
	    ListLine.push_back(eLine2);
	  }
	}
      }
    }
    eSVGplot.FrameOption=1;
  }
  if (Charac == 0) {
    clipname="torus-clipping";
    SVGclippath eClip{clipname, 0, 0, width, height};
    ListClip.push_back(eClip);
    //
    SingleBlock BlTOR=eFull.ListBlock.at("TORUS");
    double minimal=BlTOR.ListDoubleValues.at("minimal");
    double tol=BlTOR.ListDoubleValues.at("tol");
    double AngDeg=BlTOR.ListDoubleValues.at("AngDeg");
    double eScal=BlTOR.ListDoubleValues.at("scal");
    double shiftX=BlTOR.ListDoubleValues.at("shiftX");
    double shiftY=BlTOR.ListDoubleValues.at("shiftY");
    bool DrawFundamentalDomain=BlTOR.ListBoolValues.at("DrawFundamentalDomain");
    int nbTranslate=10;
    TD_description<Telt> eDesc1 = TD_FindDescriptionDirect(eCombi.VEFori.PLori, minimal, MAX_ITER_PrimalDual, tol);
    bool ApplyCaGeProcess;
    double fracArea=TD_GetMaximalAreaFraction(eDesc1);
    std::cerr << "fracArea=" << fracArea << "\n";
    if (CaGeProcessPolicy == 0)
      ApplyCaGeProcess=false;
    if (CaGeProcessPolicy == 1)
      ApplyCaGeProcess=true;
    if (CaGeProcessPolicy == 2)
      if (fracArea > MaximalFractionArea)
	ApplyCaGeProcess=true;
    TD_description<Telt> eDesc2;
    if (ApplyCaGeProcess)
      eDesc2=TD_CaGeProcessOptimization(eDesc1, MAX_ITER_CaGe, AreaToKoef);
    else
      eDesc2=eDesc1;
    TD_description<Telt> eDesc3=TD_CanonicalizationOperation(eDesc2);
    TD_description<Telt> eDescF=TD_SimilitudeOperation(eDesc3, AngDeg, eScal, shiftX, shiftY);
    auto translatePoint=[&](coor const& ePt, int const& i0, int const& i1) -> coor {
      coor ePointB=ePt;
      ePointB.x += i0*eDescF.TransLatt(0,0) + i1*eDescF.TransLatt(0,1);
      ePointB.y += i0*eDescF.TransLatt(1,0) + i1*eDescF.TransLatt(1,1);
      return ePointB;
    };
    for (int iDE=0; iDE<PL.nbP; iDE++) {
      int iDErev=PL.invers.at(iDE);
      bool IsCorr=true;
      if (iDErev < iDE)
	IsCorr=false;
      if (IsCorr == true) {
	std::vector<int> ListDE=eCombi.MapNativeDE(iDE);
	std::vector<coor> ListPoint;
	int iDE1=ListDE[0];
	int iVert1=eCombi.VEFori.ListOriginVert[iDE1];
	coor ePoint=eDescF.ListCoord[iVert1];
	ListPoint.push_back(ePoint);
	int i0=0, i1=0;
	for (auto & iDE : ListDE) {
	  i0 += eDescF.ShiftMatrix(iDE,0);
	  i1 += eDescF.ShiftMatrix(iDE,1);
	  int rDE=eCombi.VEFori.PLori.invers.at(iDE);
	  int jVert=eCombi.VEFori.ListOriginVert[rDE];
	  coor ePointB=eDescF.ListCoord[jVert];
	  coor ePointC=translatePoint(ePointB, i0, i1);
	  ListPoint.push_back(ePointC);
	}
	for (int j0=-nbTranslate; j0<=nbTranslate; j0++)
	  for (int j1=-nbTranslate; j1<=nbTranslate; j1++) {
	    std::vector<coor> TransListPoint;
	    for (auto & ePointB : ListPoint) {
	      coor transPt=translatePoint(ePointB, j0, j1);
	      TransListPoint.push_back(transPt);
	    }
	    if (IsSegmentAdmissible(TransListPoint))
	      DrawLineOfPlot(iDE, TransListPoint);
	  }
      }
    }
    int nbVert=VEFori.VertSet.size();
    for (int iVert=0; iVert<nbVert; iVert++) {
      int iDE=VEFori.VertSet[iVert][0];
      int iDEmap=eCombi.MapSingleDE(iDE);
      int iVertMap=eCombi.VEFori.ListOriginVert[iDEmap];
      coor ePt=eDescF.ListCoord[iVertMap];
      for (int j0=-nbTranslate; j0<=nbTranslate; j0++)
	for (int j1=-nbTranslate; j1<=nbTranslate; j1++) {
	  coor ePtTrans=translatePoint(ePt, j0, j1);
	  if (IsVertexAdmissible(ePtTrans))
	    DrawVertexOfPlot(iDE, ePtTrans);
	}
    }
    if (DrawFundamentalDomain) {
      MyMatrix<double> M=eDescF.TransLatt;
      std::vector<double> ListX{double(0), M(0,0), M(0,0) + M(0,1), M(0,1)};
      std::vector<double> ListY{double(0), M(1,0), M(1,0) + M(1,1), M(1,1)};
      double eXmin=VectorMin(ListX);
      double eXmax=VectorMax(ListX);
      double eYmin=VectorMin(ListY);
      double eYmax=VectorMax(ListY);
      double eXmid=(eXmin + eXmax)/2;
      double eYmid=(eYmin + eYmax)/2;
      double transX=double(1)/double(2) - eXmid;
      double transY=double(1)/double(2) - eYmid;
      for (int i=0; i<4; i++) {
	ListX[i] += transX;
	ListY[i] += transY;
      }
      std::vector<int> FundRGB=BlTOR.ListListIntValues.at("FundamentalRGB");
      int FundTraitSize=BlTOR.ListIntValues.at("FundamentalTraitSize");
      for (int i=0; i<4; i++) {
	int j=NextIdx(4,i);
	coor pt1{ListX[i], ListY[i]};
	coor pt2{ListX[j], ListY[j]};
	SVGline eLine{pt1, pt2, {FundRGB, FundTraitSize, "", clipname}};
	ListLine.push_back(eLine);
      }
    }
    eSVGplot.FrameOption=0;
    eSVGplot.scale_factor=width;
  }
  eSVGplot.ListClip=ListClip;
  eSVGplot.ListPolyline=ListPolyline;
  eSVGplot.ListBezier=ListBezier;
  eSVGplot.ListLine=ListLine;
  eSVGplot.ListEllipse=ListEllipse;
  GeneralWriteSVGfile(OutFile, eSVGplot);
  //
  std::string OutFileRed=FILE_RemoveEndingExtension(OutFile, "svg");
  for (auto & eType : ListExportFormat) {
    std::string eCommand="inkscape " + OutFile;
    if (eType == "png")
      eCommand += " --export-png=" + OutFileRed + ".png";
    if (eType == "pdf")
      eCommand += " --export-pdf=" + OutFileRed + ".pdf";
    if (eType == "eps")
      eCommand += " --export-eps=" + OutFileRed + ".eps";
    int iret=system(eCommand.c_str() );
    std::cerr << "eCommand=" << eCommand << "  iret=" << iret << "\n";
    if (iret != 0) {
      std::cerr << "Error in inkscape operation\n";
      throw TerminalException{1};
    }
  }
}


#endif
