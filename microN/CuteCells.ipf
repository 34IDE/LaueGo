#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#requiredPackages "xmlMultiIndex;imageDisplayScaling;"
#pragma ModuleName=CuteCells
#pragma version=0.23
#include "imageDisplayScaling"

Menu "Rotations"
	SubMenu("Cute Cells on Graph")
		MenuItemsWaveClassOnGraph("Add Cute Cell at Cursor A","Random3dArrays*","",csr="A"),/Q,AddCuteCellOnPoint(NaN)
		MenuItemsWaveClassOnGraph("Remove Cute Cell at Cursor A","Random3dArrays*","",csr="A"),/Q,RemoveCuteCellOnPoint(NaN)
		"Set Rotation for 2D View...", SetBeamline2DviewRot("")
	End
End

Menu "GraphMarquee", dynamic
	"-"
	MenuItemsWaveClassOnGraph("Add Cute Cell","Random3dArrays*",""),/Q,AddCuteCellOnMarqueeCenter()
	MenuItemsWaveClassOnGraph("Remove Cute Cell","Random3dArrays*",""),/Q,RemoveCuteCellInMarquee()
End



//	
//	For displaying and everything, else including choosing x and y axes, use rotation matrix, not just the switch 
//	statement assigning EBSDx & EBSDy."
//

// ============================================================================================
// =============================== Start of Set 2D View Rotation ==============================

Function SetBeamline2DviewRot(view,[printIt])
	String view
	Variable printIt
	printIt = ParamIsDefault(printIt) || numtype(printIt) ? strlen(GetRTStackInfo(2))==0 : !(!printIt)

	Init_CuteCells(0)									// in case it was not done previously
	String AllViews="H_X;Y_Z;H_F;Z_X;Y_X;"	// list of known views
	if (WhichListItem(view,AllViews) < 0)		// no valid view given, ask
		Prompt view,"type of 2D View", popup, AllViews
		DoPrompt "2D View", view
		if (V_flag)
			return 1
		endif
		printIt = 1
	endif
	if (printIt)
		printf "SetBeamline2DviewRot(\"%s\"",view
		if (!ParamIsDefault(printIt))
			printf ", printIt=%g", printIt
		endif
		printf ")\r"
	endif
	if (WhichListItem(view,AllViews) < 0)		// no valid view given, ask
		return 1
	endif

	Variable iroot2 = 0.707106781186547524		// APMath/V/N=17 str = 1/sqrt(2)
	Make/N=(3,3)/D/FREE rot=NaN		// changes REPRESENTATION, not the actual direction
	strswitch(view)
		case "H_X":							// +X --> +x,   H --> +y,   F --> +z (inward normal)
			rot[0][0] = {1, 0, 0}		// 45¡ about -X
			rot[0][1] = {0, +iroot2, -iroot2}
			rot[0][2] = {0, +iroot2, +iroot2}
			break

		case "Y_Z":							// +Z --> +x,   +Y --> +y,   X --> -z (outward normal)
			rot[0][0] = {0, 0, -1}		// 90¡ about Y
			rot[0][1] = {0, 1, 0}
			rot[0][2] = {1, 0, 0}
			break

		case "Z_X":							// +Z --> +x,   +X --> -y,   Y --> -z (outward normal)
			rot[0][0] = {0, -1, 0}		// 120¡ about {-1, 1, -1}
			rot[0][1] = {0, 0, -1}
			rot[0][2] = {1, 0, 0}
			break

		case "Y_X":							// +X --> +x,   +Y --> -y,   Z --> +z (inward normal)
			rot = p==q						// no rotation
			break

		case "H_F":							// H --> +y,   F --> +x,   +X --> -z (outward normal)
			rot[0][0] = {0, 0, -1}		// 98.4211¡ about {-0.35741,  0.86286,  0.35741}
			rot[0][1] = {-iroot2, iroot2, 0}
			rot[0][2] = { iroot2, iroot2, 0}
			break

		default:
			rot = NaN
	endswitch

	if (!(abs(MatrixDet(rot)-1)<1e-12))
		print "ERROR -- Determinant =",MatrixDet(rot),"   ********************************"
		return 1
	endif

	Wave Beamline2Dview = root:Packages:micro:CuteCells:Beamline2Dview
	if (!WaveExists(Beamline2Dview))
		Init_CuteCells(0)
		Wave Beamline2Dview = root:Packages:micro:CuteCells:Beamline2Dview
	endif
	Beamline2Dview = rot
	Note/K Beamline2Dview, ReplaceStringByKey("view","",view,"=")

	if (printIt)
		printf "Changed the rotation Beamline2Dview to: %s   (%s)\r" view, ReplaceString("_",view," vs ")
		printWave(Beamline2Dview,name="", brief=1)
//		print " "
//		test_Beamline2Dview(rot)
	endif
	return 0
End
//
Static Function test_Beamline2Dview(rot)
	Wave rot

	Variable iroot2 = 0.707106781186547524		// APMath/V/N=17 str = 1/sqrt(2)
	Make/N=3/D/FREE xBL={1,0,0}, yBL={0,1,0}, zBL={0,0,1}, HBL={0,iroot2,iroot2}, FBL={0,-iroot2,iroot2}
	String name=NameOfWave(rot)
	name = SelectString(StringMatch(name,"_free_"),name,"Beamline2Dview")

	Variable det = MatrixDet(rot)
	if (abs(det-1)>1e-12)
		printf "ERROR -- |%s| = %g\r" ,name, MatrixDet(rot)
	endif

	Make/N=3/D/FREE axis
	Variable angle = axisOfMatrix(rot,axis)
	printf "angle(%s) = %g¡ about %s,   |%s| = %g\r",name,angle,vec2str(axis), name,det
//	printf "angle(%s) = %g¡,   |%s| = %g\r",name,rotationAngleOfMat(rot), name,det

	MatrixOp/FREE ans = rot x xBL
	printf "%s x %s  =  %s\t\t\t%s\r",name,vec2str(xBL),vec2str(ans),"rot x X(BL)"
	MatrixOp/FREE ans = rot x yBL
	printf "%s x %s  =  %s\t\t\t%s\r",name,vec2str(yBL),vec2str(ans),"rot x Y(BL)"
	MatrixOp/FREE ans = rot x zBL
	printf "%s x %s  =  %s\t\t\t%s\r",name,vec2str(zBL),vec2str(ans),"rot x Z(BL)"


	MatrixOp/FREE ans = rot x HBL
	printf "%s x %s  =  %s\t\t\t%s\r",name,vec2str(HBL),vec2str(ans),"rot x H(BL)"
	MatrixOp/FREE ans = rot x FBL
	printf "%s x %s  =  %s\t\t\t%s\r",name,vec2str(FBL),vec2str(ans),"rot x F(BL)"
End

// ================================ End of Set 2D View Rotation ===============================
// ============================================================================================




// ============================================================================================
// ==================================== Start of Cute Cells ===================================

Function AddCuteCellOnMarqueeCenter()
	Init_CuteCells(1)
	Variable pnt=pointAtMarqueeCenter()
	DrawAction/L=UserFront getgroup = $("cell"+num2str(pnt))
	if (V_flag)								// already exists, so don't add another
		return NaN
	endif
	AddCuteCellOnPoint(pnt)
	if (pnt>=0 && !NumVarOrDefault("root:Packages:micro:CuteCells:showNumbers",0))
		printf "added cell on point %d\r",pnt
	endif
	return pnt
End

Function RemoveCuteCellInMarquee()
	Init_CuteCells(1)
	STRUCT box b
	xyRangeOfMarquee("",b)
	Variable pnt = FindCuteCellClosestInBox(b)
	pnt = RemoveCuteCellOnPoint(pnt)
	if (pnt>=0 && !NumVarOrDefault("root:Packages:micro:CuteCells:showNumbers",0))
		printf "removed cell on point %d\r",pnt
	endif
	return pnt
End


Static Function FindCuteCellClosestInBox(b)
	STRUCT box &b
	Variable x0,y0

	DrawAction/L=UserFront commands
	String prefix="SetDrawEnv gstart,gname= cell"
	Variable iprefix=strlen(prefix)
	Make/N=(round(strlen(S_recreation)/20))/FREE pnts=NaN
	Variable i=0, Np=0
	do
		i = strsearch(S_recreation,"SetDrawEnv gstart,gname= cell",i,2)
		if (i>=0)
			pnts[Np] = str2num(S_recreation[i+iprefix,i+iprefix+10])
			Np += 1
			i += iprefix+2
		endif
	while (i>=0)
	Redimension/N=(Np) pnts

	String yName = StringFromList(0,TraceNameList("",";",1))
	Wave yWave = TraceNameToWaveRef("",yName)
	Variable N=DimSize(yWave,0)
	Wave xWave = XWaveRefFromTrace("",yName)
	if (!WaveExists(xWave))
		Make/N=(N)/D/FREE xWaveFREE
		CopyScales/P yWave, xWaveFREE
		xWaveFREE = x
		SetScale/P x, 0, 1, "", xWaveFREE
		Wave xWave=xWaveFREE
	endif

	Make/N=(Np)/D/FREE xs, ys
	xs = xWave[pnts[p]]
	ys = yWave[pnts[p]]
	xs = xs[p] < b.xlo ? NaN : xs[p]		// remove points outside of box
	xs = xs[p] > b.xhi ? NaN : xs[p]
	ys = ys[p] < b.ylo ? NaN : ys[p]
	ys = ys[p] > b.yhi ? NaN : ys[p]

	MatrixOP/FREE dist2 = magSqr(ys - RowRepeat(y0,Np)) + magSqr(xs - RowRepeat(x0,Np))
	WaveStats/M=1/Q dist2

	if (V_minRowLoc>=0)
		return pnts[V_minRowLoc]
	endif	
	return NaN
End
//
Static Function pointAtMarqueeCenter()
	STRUCT box b
	xyRangeOfMarquee("",b)
	Variable x0=(b.xlo + b.xhi)/2, y0=(b.ylo + b.yhi)/2		// center of Marquee

	String yName = StringFromList(0,TraceNameList("",";",1))
	Wave yWave = TraceNameToWaveRef("",yName)
	Variable N=DimSize(yWave,0)
	Wave xWave = XWaveRefFromTrace("",yName)
	if (!WaveExists(xWave))
		Make/N=(N)/D/FREE xWaveFREE
		CopyScales/P yWave, xWaveFREE
		xWaveFREE = x
		SetScale/P x, 0, 1, "", xWaveFREE
		Wave xWave=xWaveFREE
	endif
	MatrixOP/FREE dist2 = magSqr(yWave - RowRepeat(y0,N)) + magSqr(xWave - RowRepeat(x0,N))
	WaveStats/M=1/Q dist2
	Variable pnt = V_minRowLoc
	return pnt
End
//
Static Function xyRangeOfMarquee(gName,b)
	String gName
	STRUCT box &b		// xlo,xhi, ylo,yhi

//		AXTYPE:left
//		CWAVE:HH
//		UNITS:µm
//		CWAVEDF:root:white2d_index2:
//		SETAXISFLAGS:/A/E=0/N=0
//		SETAXISCMD:SetAxis/A left

	String xName, yName
	if (WhichListItem("bottom", AxisList(gName))>=0)
		xName = "bottom"
	elseif (WhichListItem("top", AxisList(gName))>=0)
		xName = "top"
	endif
	if (WhichListItem("left", AxisList(gName))>=0)
		yName = "left"
	elseif (WhichListItem("right", AxisList(gName))>=0)
		yName = "right"
	endif

	GetMarquee/Z $yName, $xName
	if (!V_flag)						// box not set, bad marquee
		b.xlo=NaN ; b.xhi=NaN ; b.ylo=NaN ; b.yhi=NaN
		return 1
	endif

	b.xlo = min(V_left,V_right)
	b.xhi = max(V_left,V_right)
	b.ylo = min(V_bottom,V_top)
	b.yhi = max(V_bottom,V_top)
	return 0
End
//
Function AddCuteCellOnPoint(pnt,[win])
	Variable pnt
	String win
	win = SelectString(ParamIsDefault(win),win,"")
	if (numtype(pnt) || pnt<0)
		pnt = NumberByKey("POINT",CsrInfo(A))
	endif
	if (numtype(pnt))
		return NaN
	endif
	Init_CuteCells(1)

	RemoveCuteCellOnPoint(pnt,win=win)		// first delete in case it is there:

	Wave Beamline2Dview = root:Packages:micro:CuteCells:Beamline2Dview
	Wave XX=XX, HH=HH, gm=gm
	Variable xOff=XX[pnt], yOff=HH[pnt]
	Make/N=(3,3)/D/FREE recipBL=gm[p][q][pnt]
	MatrixOP/FREE/O recipSample = Beamline2Dview x recipBL
	MatrixOP/FREE direct = 2*PI * Inv(recipSample^t)
	Make/N=3/FREE darkGreen={0,30000,0}

	Variable showNums=NumVarOrDefault("root:Packages:micro:CuteCells:showNumbers",0)
	String title = SelectString(showNums,"",num2str(pnt))
	DrawUnitCellOnGraph(num2str(pnt),direct,xOff,yOff,SpaceGroupID="194",thick=2,color=darkGreen, arrow=1, title=title)
	return pnt
End
//
Function RemoveCuteCellOnPoint(pnt,[win])
	Variable pnt
	String win
	win = SelectString(ParamIsDefault(win),win,"")
	if (numtype(pnt) || pnt<0)
		pnt = NumberByKey("POINT",CsrInfo(A))
	endif
	if (numtype(pnt))
		return NaN
	endif

	// get start and stop for the delete
	DrawAction/L=UserFront/W=$win getgroup = $("cell"+num2str(pnt))
	if (!V_flag)
		return NaN							// could not find group of name like cell123
	endif
	DrawAction/L=UserFront/W=$win delete=V_startPos,V_endPos
	return pnt
End
//
Static Function/T DrawUnitCellOnGraph(grpName,direct,xOff,yOff,[SpaceGroupID,win,thick,color,arrow,title])
	String grpName				// name of group of draw commands
	Wave direct					// contains direct lattice (rotated and ready to put on graph)
	Variable xOff, yOff
	String SpaceGroupID
	String win					// graph name, use "" for top graph
	Variable thick
	Wave color
	Variable arrow				// arrow can be one of {0,1,2,3}
	String title				// optional title to add
	win = SelectString(ParamIsDefault(win),win,"")
	if (ParamIsDefault(SpaceGroupID) || strlen(SpaceGroupID)<1)		// SpaceGroupID not given, check wave note
		SpaceGroupID = StringByKey("SpaceGroupID",note(direct),"=")
	endif
	if (strlen(SpaceGroupID)<1)
		return ""
	endif
	if (ParamIsDefault(color))
		Make/N=3/FREE color={55000,55000,55000}	// default to a light gray
	endif
	thick = !ParamIsDefault(thick) && thick>0 ? thick : 2
	arrow = ParamIsDefault(arrow) || numtype(arrow) || arrow<0 ? 0 : limit(round(arrow),0,3)
	title = SelectString(ParamIsDefault(title),title,"")
	if (strlen(grpName)<1)
		return ""
	endif
	grpName = SelectString(isdigit(grpName[0]),"","cell")+grpName
	grpName = CleanupName(grpName,0)
	grpName = ReplaceString("__",grpName,"_")


//	win = SelectString(ParamIsDefault(win),win,"kwTopWin")
//	GetWindow/Z $win  psizeDC
//	Variable size = min(V_right-V_left, V_bottom-V_top)
//	printf "Graph area is [left,right]=[%d, %g] Æ=%d,  [top,bottom]=[%d, %d] Æ=%d\r",V_left,V_right,V_right-V_left,V_top,V_bottom,V_bottom-V_top
//	size = max(5,round(size/15))		// size in pixels
//	print "size = ",size
//	win = SelectString(StringMatch(win,"kwTopWin"),win,"")

	GetAxis/W=$win/Q left
	Variable yLo=V_min, yHi=V_max
	GetAxis/W=$win/Q bottom
	Variable xLo=V_min, xHi=V_max

	Variable size = min(abs(yLo-yHi), abs(xHi-xLo))
//	printf "Graph area is [left,right]=[%d, %g] Æ=%d,  [top,bottom]=[%d, %d] Æ=%d\r",xLo,xHi,(xHi-xLo),yHi,yLo,yHi-yLo
//	size /= 15
//	size /= 10
	size /= 8
	size *= NumVarOrDefault("root:Packages:micro:CuteCells:scale",1)

	// using coordinate system:
	//		X is to left
	//		Y is down
	//		Z is out at you
	Make/N=3/D/FREE xhat={1,0,0}, yhat={0,1,0}
	if (xLo<xHi)
		xhat *= -1
	endif
	if (yLo>yHi)
		yhat *= -1
	endif

	//	print "SpaceGroupID =",SpaceGroupID
	Variable system = LatticeSym#latticeSystem(SpaceGroupID,3)
	//	TRICLINIC=0,MONOCLINIC=1,ORTHORHOMBIC=2,TETRAGONAL=3,TRIGONAL=4,HEXAGONAL=5,CUBIC=6
	if (system==4 || system==5)		// trigonal or hexagonal, special, use a hexagonal cell
		Wave lines=MakeAllLinesAroundHexagonalCell(direct,xhat=xhat,yhat=yhat)
	else
		Wave lines = MakeAllLinesAroundCell(direct,xhat=xhat,yhat=yhat)		// returns a (12x4) matrix (each row is x0y0, x1,y1)
	endif
	lines *= size														// change to desired size
	lines[][0] += xOff												// offset so center is at (xOff, yOff)
	lines[][2] += xOff
	lines[][1] += yOff
	lines[][3] += yOff

	Make/N=3/U/W/FREE red={65535,20000,20000}, blue={25000,25000,65535}
	SetDrawLayer/W=$win UserFront
	SetDrawEnv/W=$win gstart, gname=$grpName
	Variable i, N=DimSize(lines,0)
	for (i=0;i<N;i+=1)
		if ((i==6 && N==18) || (i==4 && N==12))				// make first vertical red
			Wave wc = red
		elseif (i==8 && N==18)										// make vertical from b blue
			Wave wc = blue
		else
			Wave wc = color
		endif
		DrawOneLine(win,lines[i][0],lines[i][1], lines[i][2],lines[i][3], color=wc, arrow=arrow)
	endfor

	if (strlen(title))
		SetDrawEnv/W=$win xcoord=bottom, ycoord=left, textxjust=1, textyjust=1
		DrawText/W=$win xOff,yOff, title
	endif
	SetDrawEnv/W=$win gstop
	return grpName
End
//
Static Function/WAVE MakeAllLinesAroundHexagonalCell(cell,[xhat,yhat])		// returns a (18x4) matrix (each row is x0y0, x1,y1)
	Wave cell					// a 3x3 matrix
	Wave xhat					// NOT assumed to be normalized, default is {-100}
	Wave yhat					// NOT assumed to be normalized, default is {010}
	if (ParamIsDefault(xhat))
		Make/N=3/D/FREE xh = {-1,0,0}
	else
		Make/N=3/D/FREE xh = xhat[p]
		normalize(xh)
	endif
	if (ParamIsDefault(yhat))
		Make/N=3/D/FREE yh = {0,1,0}
	else
		Make/N=3/D/FREE yh = yhat[p]
		normalize(yh)
	endif

	Make/N=3/D/FREE a=cell[p][0], b=cell[p][1], c=cell[p][2], mid
	mid = (a+b)/2
	normalize(mid)
	mid *= norm(a)

	Variable maxLen=max(max(norm(a),norm(b)),norm(c))
	a /= maxLen
	b /= maxLen
	c /= maxLen
	mid /= maxLen

	Variable ax=MatrixDot(a,xh), ay=MatrixDot(a,yh)
	Variable bx=MatrixDot(b,xh), by=MatrixDot(b,yh)
	Variable cx=MatrixDot(c,xh), cy=MatrixDot(c,yh)
	Variable midx = MatrixDot(mid,xh), midy = MatrixDot(mid,yh)

	Make/N=(18,4)/FREE lines=NaN			// holds the 18 line segments

	// first go around base (c==0), six lines
	lines[0][0] =  ax					// x0,y0			a --> mid, mid=(a+b)/2
	lines[0][1] =  ay
	lines[0][2] =  midx				// x1,y1
	lines[0][3] =  midy

	lines[1][0] =  midx				// x0,y0			mid --> b
	lines[1][1] =  midy
	lines[1][2] =  bx					// x1,y1
	lines[1][3] =  by	

	lines[2][0] =  bx					// x0,y0			b --> -a
	lines[2][1] =  by
	lines[2][2] = -ax					// x1,y1
	lines[2][3] = -ay	

	lines[3][0] = -ax					// x0,y0			-a --> -mid
	lines[3][1] = -ay
	lines[3][2] = -midx				// x1,y1
	lines[3][3] = -midy	

	lines[4][0] = -midx				// x0,y0			-mid --> -b
	lines[4][1] = -midy
	lines[4][2] = -bx					// x1,y1
	lines[4][3] = -by	

	lines[5][0] = -bx					// x0,y0			-b --> a
	lines[5][1] = -by
	lines[5][2] =  ax					// x1,y1
	lines[5][3] =  ay	


	// second  the six verticals, all with length of c
	lines[6][0] =  ax					// x0,y0			a --> a+c
	lines[6][1] =  ay
	lines[6][2] =  ax + cx			// x1,y1
	lines[6][3] =  ay + cy

	lines[7][0] =  midx				// x0,y0			mid --> mid+c
	lines[7][1] =  midy
	lines[7][2] =  midx	+ cx		// x1,y1
	lines[7][3] =  midy + cy

	lines[8][0] =  bx					// x0,y0			b --> a+c
	lines[8][1] =  by
	lines[8][2] =  bx + cx			// x1,y1
	lines[8][3] =  by + cy

	lines[9][0] = -ax					// x0,y0			-a --> -a+c
	lines[9][1] = -ay
	lines[9][2] = -ax + cx			// x1,y1
	lines[9][3] = -ay + cy

	lines[10][0] = -midx			// x0,y0			-mid --> -mid+c
	lines[10][1] = -midy
	lines[10][2] = -midx + cx		// x1,y1
	lines[10][3] = -midy + cy

	lines[11][0] = -bx				// x0,y0			-b --> -b+c
	lines[11][1] = -by
	lines[11][2] = -bx + cx		// x1,y1
	lines[11][3] = -by + cy


	// third go around top (c==1), six lines
	lines[12][0] =  ax + cx		// x0,y0			a+c --> mid+c, mid=(a+b)/2
	lines[12][1] =  ay + cy
	lines[12][2] =  midx + cx		// x1,y1
	lines[12][3] =  midy + cy

	lines[13][0] =  midx + cx		// x0,y0			mid+c --> b+c
	lines[13][1] =  midy + cy
	lines[13][2] =  bx + cx		// x1,y1
	lines[13][3] =  by + cy

	lines[14][0] =  bx + cx		// x0,y0			b+c --> -a+c
	lines[14][1] =  by + cy
	lines[14][2] = -ax + cx		// x1,y1
	lines[14][3] = -ay + cy

	lines[15][0] = -ax + cx		// x0,y0			-a+c --> -mid+c
	lines[15][1] = -ay + cy
	lines[15][2] = -midx + cx		// x1,y1
	lines[15][3] = -midy + cy

	lines[16][0] = -midx + cx		// x0,y0			-mid+c --> -b+c
	lines[16][1] = -midy + cy
	lines[16][2] = -bx + cx		// x1,y1
	lines[16][3] = -by + cy

	lines[17][0] = -bx + cx		// x0,y0			-b+c --> a+c
	lines[17][1] = -by + cy
	lines[17][2] =  ax + cx		// x1,y1
	lines[17][3] =  ay + cy


	MatrixOp/FREE maxX0 = maxVal(col(lines,0))	// max of lines[][0]
	MatrixOp/FREE maxX1 = maxVal(col(lines,2))	// max of lines[][2]
	MatrixOp/FREE maxY0 = maxVal(col(lines,1))
	MatrixOp/FREE maxY1 = maxVal(col(lines,3))

	MatrixOp/FREE minX0 = minVal(col(lines,0))	// min of lines[][0]
	MatrixOp/FREE minX1 = minVal(col(lines,2))
	MatrixOp/FREE minY0 = minVal(col(lines,1))
	MatrixOp/FREE minY1 = minVal(col(lines,3))

	Variable xc = (max(maxX0[0],maxX1[0]) + min(minX0[0],minX1[0])) / 2
	Variable yc = (max(maxY0[0],maxY1[0]) + min(minY0[0],minY1[0])) / 2
	lines[][0] -= xc					// offset so lines are centered about (0,0)
	lines[][2] -= xc
	lines[][1] -= yc
	lines[][3] -= yc

	return lines
End
//
Static Function/WAVE MakeAllLinesAroundCell(cell,[xhat,yhat])		// returns a (12x4) matrix (each row is x0y0, x1,y1)
	Wave cell					// a 3x3 matrix
	Wave xhat					// NOT assumed to be normalized, default is {-100}
	Wave yhat					// NOT assumed to be normalized, default is {010}
	if (ParamIsDefault(xhat))
		Make/N=3/D/FREE xh = {-1,0,0}
	else
		Make/N=3/D/FREE xh = xhat[p]
		normalize(xh)
	endif
	if (ParamIsDefault(yhat))
		Make/N=3/D/FREE yh = {0,1,0}
	else
		Make/N=3/D/FREE yh = yhat[p]
		normalize(yh)
	endif

	Make/N=3/D/FREE a=cell[p][0], b=cell[p][1], c=cell[p][2]
	Variable maxLen=max(max(norm(a),norm(b)),norm(c))
	a /= maxLen
	b /= maxLen
	c /= maxLen

	Variable ax=MatrixDot(a,xh), ay=MatrixDot(a,yh)
	Variable bx=MatrixDot(b,xh), by=MatrixDot(b,yh)
	Variable cx=MatrixDot(c,xh), cy=MatrixDot(c,yh)

	Make/N=(12,4)/FREE lines=NaN			// holds the 12 line segments

	// first go around base (c==0)
	lines[0][0] = 0					// x0,y0			(0,0,0) --> (0,1,0)
	lines[0][1] = 0
	lines[0][2] = ax					// x1, y1
	lines[0][3] = ay

	lines[1][0] = ax					// x0				(0,1,0) --> (1,1,0)
	lines[1][1] = ay					// y0
	lines[1][2] = ax + bx			// x1, y1
	lines[1][3] = ay + by

	lines[2][0] = ax + bx			// x0				(1,1,0) --> (0,1,0)
	lines[2][1] = ay + by			// y0
	lines[2][2] = bx					// x1, y1
	lines[2][3] = by

	lines[3][0] = bx					// x0				(0,1,0) --> (0,0,1)
	lines[3][1] = by					// y0
	lines[3][2] = 0					// x1, y1
	lines[3][3] = 0

	// second  the four verticals, all with length of c
	lines[4][0] = 0					// x0				(0,0,0) --> (0,0,1)
	lines[4][1] = 0					// y0
	lines[4][2] = cx					// x1, y1
	lines[4][3] = cy

	lines[5][0] = ax					// x0				(1,0,0) --> (1,0,1)
	lines[5][1] = ay					// y0
	lines[5][2] = ax + cx			// x1, y1
	lines[5][3] = ay + cy

	lines[6][0] = ax + bx			// x0				(1,1,0) --> (1,1,1)
	lines[6][1] = ay + by			// y0
	lines[6][2] = ax + bx + cx	// x1, y1
	lines[6][3] = ay + by + cy

	lines[7][0] = bx					// x0				(0,1,0) --> (0,1,1)
	lines[7][1] = by					// y0
	lines[7][2] = bx + cx			// x1, y1
	lines[7][3] = by + cy

	// third go around top (c==1)
	lines[8][0] = cx					// x0,y0			(0,0,1) --> (0,1,1)
	lines[8][1] = cy
	lines[8][2] = ax + cx			// x1, y1
	lines[8][3] = ay + cy

	lines[9][0] = ax + cx			// x0				(0,1,1) --> (1,1,1)
	lines[9][1] = ay + cy			// y0
	lines[9][2] = ax + bx + cx	// x1, y1
	lines[9][3] = ay + by + cy

	lines[10][0] = ax + bx + cx	// x0				(1,1,1) --> (0,1,1)
	lines[10][1] = ay + by + cy	// y0
	lines[10][2] = bx + cx			// x1, y1
	lines[10][3] = by + cy

	lines[11][0] = bx + cx			// x0				(0,1,1) --> (0,0,1)
	lines[11][1] = by + cy			// y0
	lines[11][2] = cx					// x1, y1
	lines[11][3] = cy	


	MatrixOp/FREE maxX0 = maxVal(col(lines,0))	// max of lines[][0]
	MatrixOp/FREE maxX1 = maxVal(col(lines,2))	// max of lines[][2]
	MatrixOp/FREE maxY0 = maxVal(col(lines,1))
	MatrixOp/FREE maxY1 = maxVal(col(lines,3))

	MatrixOp/FREE minX0 = minVal(col(lines,0))	// min of lines[][0]
	MatrixOp/FREE minX1 = minVal(col(lines,2))
	MatrixOp/FREE minY0 = minVal(col(lines,1))
	MatrixOp/FREE minY1 = minVal(col(lines,3))

	Variable xc = (max(maxX0[0],maxX1[0]) + min(minX0[0],minX1[0])) / 2
	Variable yc = (max(maxY0[0],maxY1[0]) + min(minY0[0],minY1[0])) / 2
	lines[][0] -= xc					// offset so lines are centered about (0,0)
	lines[][2] -= xc
	lines[][1] -= yc
	lines[][3] -= yc

	return lines
End
//
Static Function DrawOneLine(win,x0,y0, x1,y1,[color,thick,arrow])
	String win
	Variable x0,y0, x1,y1
	Wave color
	Variable thick
	Variable arrow				// arrow can be one of {0,1,2,3}
	thick = !ParamIsDefault(thick) && thick>0 ? thick : 2
	arrow = ParamIsDefault(arrow) || numtype(arrow) || arrow<0 ? 0 : limit(round(arrow),0,3)
	Make/N=3/U/W/FREE cw={55000,55000,55000}	// default to a light gray
	if (!ParamIsDefault(color))
		cw = color[p]
	endif

	if (numtype(x0+x1+y0+y1))
		return 1					// bad values
	elseif (abs(y1-y0)<2 && abs(x1-x0)<2)
		return 1					// lines is too small to be worth drawing
	endif
//	printf "drawing line (%g, %g) --> (%g, %g) points\r",x0, y0, x1, y1
	SetDrawEnv/W=$win xcoord=bottom, ycoord=left, linethick=thick, linefgc=(cw[0],cw[1],cw[2]), arrow=arrow
	DrawLine/W=$win x0, y0, x1, y1
	return 0
End


Function Init_CuteCells(brief)
	Variable brief						// returns immediately if data folder exists
	if (brief && DataFolderExists("root:Packages:micro:CuteCells"))
		return 0
	endif
	NewDataFolder/O root:Packages
	NewDataFolder/O root:Packages:micro
	NewDataFolder/O root:Packages:micro:CuteCells
	if (exists("root:Packages:micro:CuteCells:showNumbers")!=2)
		Variable/G root:Packages:micro:CuteCells:showNumbers = 0
	endif
	if (exists("root:Packages:micro:CuteCells:scale")!=2)
		Variable/G root:Packages:micro:CuteCells:scale = 1
	endif
	if (exists("root:Packages:micro:CuteCells:Beamline2Dview")!=1)
		Make/N=(3,3)/D root:Packages:micro:CuteCells:Beamline2Dview = (p==q)
	endif
	if (exists("root:Packages:micro:CuteCells:View2EBSD")!=1)
		Make/N=(3,3)/D root:Packages:micro:CuteCells:View2EBSD /WAVE=View2EBSD = 0
		View2EBSD[0][0] = -1
		View2EBSD[1][1] =  1
		View2EBSD[2][2] = -1
	endif
End

// ===================================== End of Cute Cells ====================================
// ============================================================================================




// ============================================================================================
// ============================ Start of Utilities To Go Elsewhere ============================

// find the step size in a presumable scanned vector coordinate by examining the values
Static Function FindScalingFromVecNEW(vec,threshold,first,stepSize,Ndim,[secondTime])
	Wave vec
	Variable threshold					// a change greater than this is intentional (less is jitter)
	Variable &first						// use 	SetScale/P x first,stepSize,"",waveName
	Variable &stepSize
	Variable &Ndim
	Variable secondTime
	secondTime = ParamIsDefault(secondTime) || numtype(secondTime) ? 0 : !(!secondTime)

	if (!WaveExists(vec))
		first = NaN
		stepSize = NaN
		Ndim = 0
		return 1
	endif
	Variable N=DimSize(vec,0)
	if (!(threshold > 0) || numtype(threshold) || N<2 || WaveDims(vec)!=1)	// invalid inputs
		first = NaN
		stepSize = NaN
		Ndim = 0
		return 1
	endif

	Make/N=(N)/D/FREE vecSort = vec[p]
	SetScale/P x 0,1,"", vecSort
	Sort vecSort, vecSort
	WaveStats/M=1/Q vecSort
	N = V_npnts
	Redimension/N=(N) vecSort								// trim off any NaN's

	Make/N=(N)/FREE/D allStepSizes=NaN, subValues=NaN
	Variable i, last=NaN
	Variable midVal=NaN, iNsubsMax=-1, NsubsMax=-Inf
	Variable stepVal, Nsubs=0								// median(stepVal) and number of samples used to compute step value
	for (i=0,Ndim=0; i<N; i+=1)
		stepVal = subValues[floor(Nsubs/2)]			// median step value
		if ( (vecSort[i]-stepVal) < threshold )		// still on step = stepVal
			subValues[Nsubs] = vecSort[i]				// accumumlate another value for median
			Nsubs += 1											// number of values used get <stepVal>
		else														// a NEW step
			if (Ndim==1)
				first = stepVal								// save the first step value
				last = stepVal	
				midVal = stepVal
			elseif (Ndim > 1)
				allStepSizes[Ndim-2] = stepVal - last	// accumulate the step delta to get stepSize
				last = stepVal
				if (Nsubs > NsubsMax)						// save location of largest step (used to get first)
					NsubsMax = Nsubs
					iNsubsMax = i-Nsubs
					midVal = vecSort[floor((i-Nsubs)+(i-1))/2]	// best value of a step (used to get first)
				endif
			endif
			subValues = NaN									// re-set subValues[]
			subValues[0] = vecSort[i]						//   and start accumulating more values
			Nsubs = 1
			Ndim += 1
		endif
	endfor
	stepVal = subValues[floor(Nsubs/2)]				// median step value
	if (Ndim<2)													// just one step found
		first = stepVal
		stepSize = 0
		Ndim = 1
		return 0
	endif

	allStepSizes[Ndim-2] = stepVal - last				// last point to add to allStepSizes
	WaveStats/M=1/Q allStepSizes							// last get the median step size
	if (Ndim>1)
		stepSize = abs(allStepSizes[floor((Ndim-2)/2)])
		// round stepSize to 9 significant figures
		stepSize = (placesOfPrecision(roundSignificant(stepSize,9))<7) ? roundSignificant(stepSize,9) : stepSize
		threshold = max(stepSize/3, threshold)		// re-set threshold based on step size
		if (!secondTime)										// make a second pass with a threshold more suitable
			FindScalingFromVecNEW(vec,threshold,first,stepSize,Ndim,secondTime=1)
		endif
		// recalc first and Ndim using midVal and stepSize, vec in range [vecSort[0], vecSort[N-1]]
		first = midVal - round((midval-first)/stepSize)*stepSize
		Ndim = round((vecSort[N-1]-first)/stepSize) + 1
	else
		stepSize = 0
	endif
	return 0
End
//
//Function test_FindScalingFromVecNEW(vec,[threshold,test123])
//	Wave vec		// usually something like XX, depth, FF, etc.
//	Variable threshold
//	Variable test123
//	threshold = ParamIsDefault(threshold) || numtype(threshold) || threshold<=0 ? 0.2 : threshold
//	test123 = ParamIsDefault(test123) || numtype(test123) ? 0 : !(!test123)
//
//	if (!WaveExists(vec))
//		String name, wlist="aaa;XX;YY;ZZ;HH;FF;depth;dddd"
//		Variable i,N=ItemsInList(wlist)
//		for (i=N-1;i>=0;i-=1)
//			name = StringFromList(i,wlist)
//			Wave vec=$name
//			if (!WaveExists(vec))
//				wlist = RemoveFromList(name,wlist)
//			endif
//		endfor
//		Prompt name, "Wave to check",popup,wlist
//		DoPrompt "pick a vector",name
//		if (V_flag)
//			return 1
//		endif
//		Wave vec=$name
//	endif
//	if (!WaveExists(vec))
//		return 1
//	endif
//
//	Variable first,stepSize,dimN
//
//	if (WaveExists(vec))
//		CuteCells#FindScalingFromVecNEW(vec,threshold,first,stepSize,dimN)
//		printf "'%s',    start = %g,   Æ = %g,   N= %g\r",NameOfWave(vec),first,stepSize,dimN
//	endif
//
//	if (!test123)
//		return 0
//	endif
//	Make/FREE one={1}, two={1,2}, three={1,2,3}
//	Variable OK123=1
//	CuteCells#FindScalingFromVecNEW(one,threshold,first,stepSize,dimN)
//	if (numtype(first)==0)
//		printf "ERROR -- '%s',    start = %g,   Æ = %g,   N= %g\r",vec2str(one),first,stepSize,dimN
//		OK123 = 0
//	endif
//	CuteCells#FindScalingFromVecNEW(two,threshold,first,stepSize,dimN)
//	if (!(first==1 && stepSize==1 && dimN==2))
//		printf "ERROR -- '%s',    start = %g,   Æ = %g,   N= %g\r",vec2str(two),first,stepSize,dimN
//		OK123 = 0
//	endif
//	CuteCells#FindScalingFromVecNEW(three,threshold,first,stepSize,dimN)
//	if (!(first==1 && stepSize==1 && dimN==3))
//		printf "ERROR -- '%s',    start = %g,   Æ = %g,   N= %g\r",vec2str(three),first,stepSize,dimN
//		OK123 = 0
//	endif
//	if (OK123)
//		print "vectors of length 1, 2, or 3 test OK"
//	endif
//End

// ============================= End of Utilities to Go Elsewhere =============================
// ============================================================================================


