#ifndef INCLUDE_SVGfunctions_h
#define INCLUDE_SVGfunctions_h




struct coor {
  double x;
  double y;
};


coor operator+(coor const& c1, coor const& c2)
{
  return {c1.x + c2.x, c1.y + c2.y};
}


coor operator-(coor const& c1, coor const& c2)
{
  return {c1.x - c2.x, c1.y - c2.y};
}


coor MultScal(coor const& c, double const& eScal)
{
  return {eScal*c.x, eScal*c.y};
}


coor IsobarycenterPoint(std::vector<coor> const& ListPt)
{
  double eX=0;
  double eY=0;
  for (auto & ePt : ListPt) {
    eX += ePt.x;
    eY += ePt.y;
  }
  double eF=1/double(ListPt.size());
  return {eX*eF, eY*eF};
}


coor RotateCoor(coor const& c, double const& eAngDeg)
{
  double eAngRad=eAngDeg*(M_PI/double(180));
  //  std::cerr << "eAngRad=" << eAngRad << "\n";
  double eCos = cos(eAngRad);
  double eSin = sin(eAngRad);
  double cX = c.x*eCos - c.y*eSin;
  double cY = c.x*eSin + c.y*eCos;
  return {cX, cY};
}


double norm(coor const& c)
{
  return sqrt(c.x*c.x+c.y*c.y);
}


std::vector<coor> GetListTangent(std::vector<coor> const& ListPoint, double const& eScal)
{
  int len=ListPoint.size();
  std::vector<coor> ListTangent(len);
  for (int i=0; i<len; i++) {
    coor eTangent;
    if (i == 0) {
      eTangent=ListPoint[1] - ListPoint[0];
    }
    else {
      if (i == len-1) {
	eTangent = ListPoint[len-1] - ListPoint[len-2];
      }
      else {
	coor coorSum=ListPoint[i+1] - ListPoint[i-1];
	eTangent = MultScal(coorSum, double(1)/double(2));
      }
    }
    ListTangent[i] = MultScal(eTangent, eScal);
  }
  return ListTangent;
}


struct SVGclippath {
  std::string name;
  double x;
  double y;
  double width;
  double height;
};


struct SVGqualInfo {
  std::vector<int> color;
  int Size;
  std::string MarkerEnd;
  std::string clip;
};


struct SVGbezier {
  coor pointM;
  coor pointC;
  coor point1;
  coor point2;
  SVGqualInfo eQual;
};


struct SVGpolyline {
  std::vector<coor> ListCoor;
  SVGqualInfo eQual;
};


struct SVGline {
  coor ePt;
  coor fPt;
  SVGqualInfo eQual;
};


struct SVGellipse {
  coor c;
  coor r;
  SVGqualInfo eQual;
};


struct SVGplotDescription {
  std::vector<SVGclippath> ListClip;
  std::vector<SVGpolyline> ListPolyline;
  std::vector<SVGbezier> ListBezier;
  std::vector<SVGline> ListLine;
  std::vector<SVGellipse> ListEllipse;
  double height;
  double width;
  double scale_factor;
  double add_offsetX;
  double add_offsetY;
  int FrameOption;
  int RoundMethod;
};


void GeneralWriteSVGfile(std::string const& eFile, SVGplotDescription const& eSVGplot)
{
  //
  // First compute the min/max of the figure.
  //
  double MinX=0, MaxX=0, MinY=0, MaxY=0;
  auto UpdateMinMaxXY=[&](coor const& pt) -> void {
    if (pt.x > MaxX)
      MaxX=pt.x;
    if (pt.x < MinX)
      MinX=pt.x;
    if (pt.y > MaxY)
      MaxY=pt.y;
    if (pt.y < MinY)
      MinY=pt.y;
  };
  for (auto & eLine : eSVGplot.ListLine) {
    UpdateMinMaxXY(eLine.ePt);
    UpdateMinMaxXY(eLine.fPt);
  }
  for (auto & ePolyline : eSVGplot.ListPolyline)
    for (auto & eCoor : ePolyline.ListCoor)
      UpdateMinMaxXY(eCoor);
  for (auto& eBez : eSVGplot.ListBezier) {
    UpdateMinMaxXY(eBez.pointM);
    UpdateMinMaxXY(eBez.point2);
  }
  for (auto& eEll : eSVGplot.ListEllipse) {
    UpdateMinMaxXY(eEll.c);
  }
  std::cerr << "SVG: X(min/max)=" << MinX << " / " << MaxX << "\n";
  std::cerr << "SVG: Y(min/max)=" << MinY << " / " << MaxY << "\n";
  //
  // Second compute the FrameOption
  //
  double scale_factor, add_offsetX, add_offsetY;
  double height=0, width=0;
  bool FrameInit=false;
  // We choose the scaling factor directly from the SVGplot
  if (eSVGplot.FrameOption == 0) {
    height=eSVGplot.height;
    width=eSVGplot.width;
    scale_factor=eSVGplot.scale_factor;
    add_offsetX=eSVGplot.add_offsetX;
    add_offsetY=eSVGplot.add_offsetY;
    FrameInit=true;
  }
  // We have a fixed frame in X and Y and we adjust accordingly.
  if (eSVGplot.FrameOption == 1) {
    height=eSVGplot.height;
    width=eSVGplot.width;
    double FrameX=eSVGplot.width;
    double FrameY=eSVGplot.height;
    double scale_factorX=FrameY / (MaxX - MinX);
    double scale_factorY=FrameX / (MaxY - MinY);
    double MidX=(MaxX + MinX) / 2;
    double MidY=(MaxY + MinY) / 2;
    scale_factor=std::min(scale_factorX, scale_factorY);
    // add_offsetX + scale_factor*MidX = FrameX/2;
    add_offsetX=FrameX/2 - scale_factor*MidX;
    add_offsetY=FrameY/2 - scale_factor*MidY;
    FrameInit=true;
  }
  if (FrameInit == false) {
    std::cerr << "FrameOption has not been used\n";
    std::cerr << "Please correct this\n";
    throw TerminalException{1};
  }
  auto GetStringValue=[&](double const& eVal, double const& add_offset) -> std::string {
    double eValM=add_offset + eVal*scale_factor;
    if (eSVGplot.RoundMethod == 1)
      return DoubleTo4dot2f(eValM);
    if (eSVGplot.RoundMethod == 2)
      return DoubleToString(eValM);
    if (eSVGplot.RoundMethod == 3)
      return DoubleToString(eValM);
    std::cerr << "Failed to find relevant function\n";
    throw TerminalException{1};
  };
  auto GetStringValueX=[&](double const& eVal) -> std::string {
    return GetStringValue(eVal, add_offsetX);
  };
  auto GetStringValueY=[&](double const& eVal) -> std::string {
    return GetStringValue(eVal, add_offsetY);
  };
  auto GetStringPair=[&](coor const& pt) -> std::string {
    return GetStringValueX(pt.x) + " " + GetStringValueY(pt.y);
  };
  //
  // The attribute functionalities
  //
  auto StringColor=[&](std::vector<int> const& eV) -> std::string {
    return "rgb(" + IntToString(eV[0]) + "," + IntToString(eV[1]) + "," + IntToString(eV[2]) + ")";
  };
  auto GetQualityString=[&](SVGqualInfo const& eQual) -> std::string {
    std::string eRet="style=\"stroke:" + StringColor(eQual.color) + ";stroke-width:" + IntToString(eQual.Size) + "\"";
    if (eQual.MarkerEnd != "")
      eRet += " marker-end=\"url(#" + eQual.MarkerEnd + ")\"";
    if (eQual.clip != "")
      eRet += " clip-path=\"url(#" + eQual.clip + ")\"";
    return eRet;
  };
  auto GetQualityStringEllipse=[&](SVGqualInfo const& eQual) -> std::string {
    std::string eRet=" style=\"fill:" + StringColor(eQual.color) + "\"";
    if (eQual.clip != "")
      eRet += " clip-path=\"url(#" + eQual.clip + ")\"";
    return eRet;
  };
  //
  // First the preamble
  //
  std::ofstream os(eFile);
  os << "<svg height=\"" << height << "\" width=\"" << width << "\">\n";
  //
  // Printing the clipping path information
  //
  for (auto & eClip : eSVGplot.ListClip) {
    os << "  <defs>\n";
    os << "    <clipPath id=\"" << eClip.name << "\">\n";
    os << "      <rect x=\"" << eClip.x << "\" y=\"" << eClip.y << "\" width=\"" << eClip.width << "\" height=\"" << eClip.height << "\" />\n";
    os << "    </clipPath>\n";
    os << "  </defs>\n";
  }
  //
  // Single line between points
  //
  std::cerr << "|ListLine|=" << eSVGplot.ListLine.size() << "\n";
  for (auto & eLine : eSVGplot.ListLine) {
    coor ePt=eLine.ePt;
    coor fPt=eLine.fPt;
    os << "  <line x1=\"" << GetStringValueX(ePt.x) << "\" y1=\"" << GetStringValueY(ePt.y) << "\" x2=\"" << GetStringValueX(fPt.x) << "\" y2=\"" << GetStringValueY(fPt.y) << "\" " << GetQualityString(eLine.eQual) << " />\n";
  }
  //
  // Polylines
  //
  std::cerr << "|ListPolyline|=" << eSVGplot.ListPolyline.size() << "\n";
  for (auto & ePolyline : eSVGplot.ListPolyline) {
    os << "<polyline points=\"";
    bool IsFirst=true;
    for (auto & ePt : ePolyline.ListCoor) {
      if (IsFirst == false)
	os << " ";
      os << GetStringValueX(ePt.x) << "," << GetStringValueY(ePt.y);
    }
    os << "\" " << GetQualityString(ePolyline.eQual) << " />\n";
  }
  //
  // Bezier curves
  //
  std::cerr << "|ListBezier|=" << eSVGplot.ListBezier.size() << "\n";
  for (auto& eBez : eSVGplot.ListBezier) {
    os << "  <path d=\"M" << GetStringPair(eBez.pointM) << " C " << GetStringPair(eBez.pointC) << ", " << GetStringPair(eBez.point1) << ", " << GetStringPair(eBez.point2) << "\" fill=\"none\" " << GetQualityString(eBez.eQual) << " />\n";
  }
  //
  // Ellipses
  //
  std::cerr << "|ListEllipse|=" << eSVGplot.ListEllipse.size() << "\n";
  for (auto & eEll : eSVGplot.ListEllipse) {
    os << "  <ellipse";
    os << " cx=\"" << GetStringValueX(eEll.c.x) << "\"";
    os << " cy=\"" << GetStringValueY(eEll.c.y) << "\"";
    os << " rx=\"" << GetStringValue(eEll.r.x, 0) << "\"";
    os << " ry=\"" << GetStringValue(eEll.r.y, 0) << "\"";
    os << " " << GetQualityStringEllipse(eEll.eQual);
    os << " />\n";
  }
  //
  // Finishing it
  //
  os << "</svg>\n";
}


#endif
