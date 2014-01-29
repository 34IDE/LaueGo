#pragma rtGlobals=1		// Use modern global access method.
#pragma ModuleName=multiIndex
#pragma version=1.68
#include "microGeometryN", version>=1.15
#include "LatticeSym", version>=3.41
//#include "DepthResolvedQueryN"
#include "ArrayOf3dOrientsN", version>=2.41
#include "IndexLotsN", version>=2.18
#include "GizmoZoomTranslate"

Static Constant NindexedMIN = 4
Static StrConstant XMLfilters = "XML Files (*.xml,*.txt):.xml,.txt;All Files:.*;"


Menu "Rotations"
	"-"
	"(XML [from cluster]"
	MenuItemIfWaveClassExists("Pole Figure","Random3dArraysGm",""),NewPoleFigure("")
	MenuItemIfWaveClassExists("  Re-plot Pole Figure","poleXYpoints",""),Make_GraphPoles($"")
	"  Color Hexagon",MakColorHexagon()
	"  Color Triangle [CUBIC]",MakeColorTriangle()
	MenuItemIfWaveClassExists("Strain Refine Many Images","Random3dArraysGm",""),DeviatoricStrainRefineXML_ALL("","")
	MenuItemIfWaveClassExists("  Strain Refine One Image","Random3dArraysGm",""),DeviatoricStrainRefineXML(-1,NaN,"",printit=1)
	MenuItemIfWaveClassExists("Histogram a 3D Array","Interpolated3dArrays",""),multiIndex#HistogramOf3dArray($"")
	MenuItemIfWaveClassExists("  Re-Plot Histogram of 3D Array","HistogramFrom3d",""),multiIndex#DisplayHistogramFrom3d($"")
	" (  --------"
//	"Load XML [micro-diffraction file]",Load3dRecipLatticesFileXML("")
	"Load XML [micro-diffraction file]",Load3dRecipLatticesFileXML("") ; print " " ; ProcessLoadedXMLfile(Inf,NaN)
	multiIndex#MenuItemIfValidRawDataExists("Re-Process Loaded XML"), ProcessLoadedXMLfile(Inf,NaN)
	MenuItemIfWaveClassExists("  Load Pixels from one step in XML file","Random3dArrays",""),/Q,LoadPixelsFromXML(-1)
	MenuItemIfWaveClassExists("2D plot of loaded XML data","Random3dArrays",""), Make2Dplot_xmlData("")
	MenuItemIfWaveClassExists("Make Gizmo of 3D xml data","Random3dArrays",""),MakeGizmo_xmlData($"")
	MenuItemIfWaveClassExists("Make another RGBA for a Gizmo","Random3dArraysXYZ",""),MakeGizmoRGBAforXYZ($"",$"","",NaN,NaN)
	MenuItemIfWaveClassExists("Add Hand Indexed Point to an XML","IndexedPeakList",""),AppendCurrenIndexResult2XML($"","")
	MenuItemIfWaveClassExists("Rotation Between Two Points","Random3dArrays",""),RotationBetweenTwoPoints(NaN,NaN)
	"-"
EndMacro
Menu "Data"
	"Load XML [micro-diffraction file]",Load3dRecipLatticesFileXML("")
End
Menu "Load Waves"
	"Load XML [micro-diffraction file]",Load3dRecipLatticesFileXML("")
End

Static Function/T MenuItemIfValidRawDataExists(item)
	String item
	Variable valid=ValidRawXMLdataAvailable()
	return SelectString(valid,"(","")+item
End


// *****************************************************************************************
// ****************************  Start of Fill Movie More Stuff ****************************

//StrConstant XMLfile = "xmdmac HD:Users:tischler:data:Budai June09:areaTemp.xml"
//StrConstant XMLfile = "work1:July09:Cu001_Indent:Raw1_WW.xml"
//Function FillMovieMore(image,[index])
Function FillMovieMore(image)
	Wave image
//	Variable index						// not used
	String noteStr = note(image)
	String imageName = StringByKey("imageFilePath",noteStr,"=")+StringByKey("imageFileName",noteStr,"=")
	Wave waveXY=waveXY
	if (!WaveExists(waveXY))
		return 1
	endif
 
	// get xml file with values for fitted peaks
	SVAR XMLfileName=XMLfileName
	if (!SVAR_Exists(XMLfileName))
		Variable f
		Open/D/R/T=".xml"/M="XML file" f
		print S_fileName
		if (strlen(S_fileName)<1)
			return 1
		endif
		String/G XMLfileName=S_fileName
	endif

	if (xmlPixelinfoForMovies(XMLfileName,imageName,waveXY)<1)
		waveXY = NaN
	endif
	return 0
End
//
// This routine goes to the xml file and gets the peaksearch and indexing information, which it adds to the plots.
Static Function xmlPixelinfoForMovies(FullFileName,imageName,waveXY)		// returns number of points read into waveXY
	String FullFileName
	String imageName
	Wave waveXY
	imageName = ParseFilePath(3, imageName, ":", 0, 0)

	if (!WaveExists(waveXY))
		Abort "waveXY does not exits"
	endif

	String pathSep
	if (stringmatch(IgorInfo(2),"Macintosh"))
		pathSep = ":"
	elseif (stringmatch(IgorInfo(2),"Windows"))
		pathSep = "\\"
	else
		Abort "This is neither Mac nor Win, what is this computer?"
	endif
	Variable f
	Open/R/P=home/Z=1 f as FullFileName
	if (strlen(S_fileName)<1 || V_flag)
		return 0
	endif
	Close f
	FullFileName = S_fileName

	String names = MakeStepIndexForXML(FullFileName)
	Wave indexPos = $StringFromLIst(0,names)		// start of the <step ...>
	Wave/T indexImage = $StringFromLIst(1,names)
	if (!WaveExists(indexPos) || !WaveExists(indexImage))
		return 0
	endif

	// read step containing imageName from file
	Variable Nindex=DimSize(indexImage,0), index
	for (index=0;index<Nindex;index+=1)				// loop over indexImage to find point number for imageName
		if (stringmatch(indexImage[index],imageName))
			break
		endif
	endfor
	if (index>=Nindex)									// image name not found in indexImage[]
		return 0
	endif
	Variable i0=indexPos[index], i1=indexPos[index+1]
	if (numtype(i0+i1) || i1<=i0)						// invalid range in file [i0,i1]
		return 0
	endif

	String step = PadString("",i1-i0+1,0x20)			// a buffer to hold step (of length i1-i0+1 bytes)
	Open/R/Z=1 f as FullFileName						// start reading here
	FStatus f
	Variable fileLen = V_logEOF							// length of file in bytes
	if (i1>=V_logEOF)									// file too short
		Close f
		return 0
	endif
	FSetPos f, i0
	FBinRead f, step										// initial read
	Close f

	String svec = xmlTagContents("Xpixel",step)		// got <step>, extract data
	svec = ReplaceString(" ",svec,";")
	Variable N=ItemsInList(svec)
	if (N<=0)
		return N
	endif

	Redimension/N=(N,2) waveXY
	SetScale d 0,0,"pixel", waveXY
	waveXY[][0] = str2num(StringFromList(p,svec))
	svec = xmlTagContents("Ypixel",step)
	svec = ReplaceString(" ",svec,";")
	waveXY[][1] = str2num(StringFromList(p,svec))

	Make/N=3/U/W/FREE green={0,65535,0},blue={0,0,65535}, lightBlue={16385,49025,65535}
	Wave waveXYcolor=waveXYcolor, waveXYmarker=waveXYmarker
	if (!WaveExists(waveXYcolor))
		return N
	endif
	Redimension/N=(N,3)/U/W waveXYcolor
	waveXYcolor[][] = lightBlue[q]				// set initial color to light-blue (means un-indexed)
	if (WaveExists(waveXYmarker))
		Redimension/N=(N)/U/W waveXYmarker
		waveXYmarker = 5						// set initial marker shape to open-square
	endif

	Variable px,py,i, dnum=detectorNumFromID(xmlTagContents("detectorID",step))
	Variable depth=str2num(xmlTagContents("depth",step))
	depth = numtype(depth) ? 0 : depth
	String sROI=XMLattibutes2KeyList("ROI",step)
	Variable startx=NumberByKey("startx",sROI,"=")
	Variable starty=NumberByKey("starty",sROI,"=")
	Variable groupx=NumberByKey("groupx",sROI,"=")
	Variable groupy=NumberByKey("groupy",sROI,"=")
	startx = numtype(startx) ? 0 : startx
	starty = numtype(starty) ? 0 : starty
	groupx = numtype(groupx) ? 1 : groupx
	groupy = numtype(groupy) ? 1 : groupy
	String sindex=xmlTagContents("indexing",step)
	Variable x1,x2,x3
	Variable/C pz
	String hvec,kvec,lvec
	STRUCT microGeometry g
	FillGeometryStructDefault(g)
	Variable ipattern, Nindexed, marker
	for (ipattern=0;ipattern<10;ipattern+=1)		// loop over all identified patterns
		Nindexed=NumberByKey("Nindexed",XMLattibutes2KeyList("pattern",sindex,occurance=ipattern),"=")
		if (!(Nindexed>2))							// when nothing indexed, it is time to stop looking
			break
		endif
		String spattern=xmlTagContents("pattern",sindex,occurance=ipattern)
		sscanf xmlTagContents("astar",spattern),"%g %g %g",x1,x2,x3
		if (V_flag != 3)
			return N
		endif
		Make/N=(3,3)/D/FREE recip=NaN
		recip[0][0]=x1;	recip[1][0]=x2;	recip[2][0]=x3
		sscanf xmlTagContents("bstar",spattern),"%g %g %g",x1,x2,x3
		recip[0][1]=x1;	recip[1][1]=x2;	recip[2][1]=x3
		sscanf xmlTagContents("cstar",spattern),"%g %g %g",x1,x2,x3
		recip[0][2]=x1;	recip[1][2]=x2;	recip[2][2]=x3

		hvec=xmlTagContents("h",spattern)
		hvec = ReplaceString(" ",hvec,";")
		kvec=xmlTagContents("k",spattern)
		kvec = ReplaceString(" ",kvec,";")
		lvec=xmlTagContents("l",spattern)
		lvec = ReplaceString(" ",lvec,";")
		Make/N=3/D/FREE hkl
		for (i=0;i<Nindexed;i+=1)
			hkl[0] = str2num(StringFromList(i,hvec))
			hkl[1] = str2num(StringFromList(i,kvec))
			hkl[2] = str2num(StringFromList(i,lvec))
			MatrixOP/FREE/O qvec = recip x hkl
			pz = q2pixel(g.d[dnum],qvec,depth=depth)
			px = (real(pz)-startx - (groupx-1)/2)/groupx
			py = (imag(pz)-starty - (groupy-1)/2)/groupy
			// find pixel closest to pixel in waveXY, and set its color to light-green
			Make/N=2/D/FREE pix={px,py}
			MatrixOP/FREE dist = sqrt(sumRows(magSqr(waveXY-RowRepeat(pix,N))))
			WaveStats/M=1/Q dist
			if (V_min<4)							// within 4 pixesl, set it
				if (ipattern==0)					// pattern 0 is green circles
					waveXYcolor[V_minloc][] = green[q]
					marker = 8
				else									// patterns 1-Inf are blue diamonds
					waveXYcolor[V_minloc][] = blue[q]
					marker = 7
				endif
				if (WaveExists(waveXYmarker))
					waveXYmarker[V_minloc] = marker
				endif
			endif
		endfor
	endfor
	return N
End

// *****************************  End of Fill Movie More Stuff *****************************
// *****************************************************************************************



// *****************************************************************************************
// *****************************  Start of Pole Figure Colors ******************************

Function NewPoleFigure(hklstr,[rad,useCursor,useDialog,gm,visible,iwave])
	String hklstr
	Variable rad				// radius of color circle (degree)
	Variable useCursor
	Variable useDialog			// forces the use of a dialog
	Wave gm
	Wave visible				// flag for each point indicating which points to use (1=use, 0=ignore)
	Wave iwave					// used to set intensity of pixels
	if (ParamIsDefault(gm))
		Wave gm=$""
	endif
	String gmList=WaveListClass("Random3dArraysGm","*","MAXROWS:3,MAXCOLS:3,MINCOLS:3,MINROWS:3,MINLAYERS:2"), gmName
	if (ItemsInList(gmList)==1 && !WaveExists(gm))
		Wave gm=$StringFromList(0,gmList)
	endif
	if (ParamIsDefault(visible))
		Wave visible=$""
	endif
	String visList=WaveListClass("GizmoVisible","*","MINROWS:2,MAXCOLS:0"), visName=""
	Variable askVis = ItemsInList(visList)>0 && !WaveExists(visible)
	if (WaveExists(visible))
		visName = NameOfWave(visible)
	endif
	String iwaveName=""
	if (WaveExists(iwave))
		iwaveName = NameOfWave(iwave)
	endif
	rad = ParamIsDefault(rad) ? 22.5 : rad
	useDialog = ParamIsDefault(useDialog) ? 0 : useDialog
	useDialog = numtype(useDialog) ? 1 : useDialog
	useCursor = ParamIsDefault(useCursor)||numtype(useCursor) ? 1 : !(!useCursor)
	Variable cursorAvailable=0
	if (ItemsInLIst(WinList("*",";","WIN:1"))>0)
		if (stringmatch(StringByKey("TNAME", CsrInfo(A)),"pole*XY"))
			cursorAvailable = numtype(hcsr(A)+vcsr(A))==0
		endif
	endif
	useCursor = cursorAvailable ? useCursor : 0
	useCursor = !(!useCursor)
	hklstr = ReplaceString(",",hklstr," ")
	hklstr = ReplaceString(";",hklstr," ")
	hklstr = ReplaceString("(",hklstr," ")
	hklstr = ReplaceString(")",hklstr," ")
	Variable h,k,l
	if (strlen(hklstr)==3)						// try inserting spaces
		sscanf hklstr, "%1d%1d%1d", h,k,l
	else
		sscanf hklstr, "%g %g %g", h,k,l
	endif
	if (V_flag!=3 || !(rad>0) || useDialog || !WaveExists(gm) || askVis)
		String hklLIst = "0 0 1;0 1 1;1 1 1;other"
		if (WhichListItem(hklstr,hklLIst)<0)
			hklLIst = hklstr + ";" + hklLIst
		endif
		Prompt rad,"radius of color circle (0<rad<90¡)"
		Prompt useCursor,"Center of Color Circle",popup,"Origin at Cursor;Origin at Center"
		Prompt hklstr,"hkl of pole",popup,hklLIst
		Prompt gmName,"Reciprocal Lattice List",popup,gmList
		Prompt iwaveName,"Wave for Intensities of Colors",popup,WaveListClass("Random3dArrays","*","DIMS:1,TEXT:0")
		iwaveName = SelectString(strlen(iwaveName),"Nindexed",iwaveName)
		rad = rad>0 ? rad : 22.5
		useCursor = useCursor ? 1 : 2			// convert to {1,2} for prompt
		if (WaveExists(gm))
			if (cursorAvailable)
				DoPrompt "hkl",hklstr,rad,useCursor,iwaveName
			else
				DoPrompt "hkl",hklstr,rad,iwaveName
			endif
		else
			if (cursorAvailable)
				DoPrompt "hkl",hklstr,rad,useCursor,gmName,iwaveName
			else
				DoPrompt "hkl",hklstr,rad,gmName,iwaveName
			endif
			Wave gm=$gmName
		endif
		if (V_flag)
			return 1
		endif
		Wave iwave=$iwaveName
		useCursor = useCursor==2 ? 0 : 1
		if (strlen(hklstr)==3)						// try inserting spaces
			sscanf hklstr, "%1d%1d%1d", h,k,l
		else
			sscanf hklstr, "%g %g %g", h,k,l
		endif
		if (V_flag!=3)
			h = 0
			k = 0
			l = 1
			Prompt h,"H"
			Prompt k,"K"
			Prompt l,"L"
			DoPrompt "hkl",h,k,l
			if (V_flag)
				return 1
			endif
		endif
		if (askVis)
			Prompt visName,"Visible Mask from Gizmo",popup," none ;"+visList
			DoPrompt "Visible Mask",visName
			if (V_flag)
				return 1
			endif
			if (stringmatch(visName," none "))
				visName=""
			else
				Wave visible=$visName
			endif
		endif
	endif
	if (!WaveExists(gm))
		DoAlert 0,"gm wave does not exist"
		return 1
	endif

	String str = SelectString(WhichListItem("PolePopMenuProc",GetRTStackInfo(0))>=0, "","¥"), cmd
	if (WaveExists(iwave))
		iwaveName = NameOfWave(iwave)
	else
		iwaveName = ""
	endif
	sprintf cmd, "%sNewPoleFigure(\"%g %g %g\",rad=%g, gm=%s, iwave=%s",str,h,k,l,rad,NameOfWave(gm),iwaveName
	if (WaveExists(visible))
		cmd += ", visible="+NameOfWave(visible)
	endif
	if (cursorAvailable)
		sprintf str, ", useCursor=%g)\t\t//cursor = %d\r",useCursor,pcsr(A)
		cmd += str
	else
		cmd += ")\r"
	endif
	print cmd
	if (numtype(h+k+l))
		return 1
	endif
	if (!(rad>0) || !(rad<=90))
		DoAlert 0,"rad is not in range (0,90¡]"
		return 1
	endif

	Make/N=3/O/D/FREE pole
	pole={h,k,l}
	Wave poleXY = $MakePolePoints(gm,pole,rad=rad,useCursor=useCursor, visible=visible,iwave=iwave)
	Make_GraphPoles(poleXY)
	return 0
End
//Function NewPoleFigure(hklstr,[rad,useCursor,useDialog])
//	String hklstr
//	Variable rad				// radius of color circle (degree)
//	Variable useCursor
//	Variable useDialog			// forces the use of a dialog
//
//	rad = ParamIsDefault(rad) ? 22.5 : rad
//	useDialog = ParamIsDefault(useDialog) ? 0 : useDialog
//	useDialog = numtype(useDialog) ? 1 : useDialog
//	useCursor = ParamIsDefault(useCursor)||numtype(useCursor) ? 1 : !(!useCursor)
//	Variable cursorAvailable=0
//	if (ItemsInLIst(WinList("*",";","WIN:1"))>0)
//		if (stringmatch(StringByKey("TNAME", CsrInfo(A)),"pole*XY"))
//			cursorAvailable = numtype(hcsr(A)+vcsr(A))==0
//		endif
//	endif
//	useCursor = cursorAvailable ? useCursor : 0
//	useCursor = !(!useCursor)
//	hklstr = ReplaceString(",",hklstr," ")
//	hklstr = ReplaceString(";",hklstr," ")
//	hklstr = ReplaceString("(",hklstr," ")
//	hklstr = ReplaceString(")",hklstr," ")
//	Variable h,k,l
//	if (strlen(hklstr)==3)						// try inserting spaces
//		sscanf hklstr, "%1d%1d%1d", h,k,l
//	else
//		sscanf hklstr, "%g %g %g", h,k,l
//	endif
//	if (V_flag!=3 || !(rad>0) || useDialog)
//		String hklLIst = "0 0 1;0 1 1;1 1 1;other"
//		if (WhichListItem(hklstr,hklLIst)<0)
//			hklLIst = hklstr + ";" + hklLIst
//		endif
//		Prompt rad,"radius of color circle (0<rad<90¡)"
//		Prompt useCursor,"Center of Color Circle",popup,"Origin at Cursor;Origin at Center"
//		Prompt hklstr,"hkl of pole",popup,hklLIst
//		rad = rad>0 ? rad : 22.5
//		useCursor = useCursor ? 1 : 2			// convert to {1,2} for prompt
//		if (cursorAvailable)
//			DoPrompt "hkl",hklstr,rad,useCursor
//		else
//			DoPrompt "hkl",hklstr,rad
//		endif
//		if (V_flag)
//			return 1
//		endif
//		useCursor = useCursor==2 ? 0 : 1
//		if (strlen(hklstr)==3)						// try inserting spaces
//			sscanf hklstr, "%1d%1d%1d", h,k,l
//		else
//			sscanf hklstr, "%g %g %g", h,k,l
//		endif
//		if (V_flag!=3)
//			h = 0
//			k = 0
//			l = 1
//			Prompt h,"H"
//			Prompt k,"K"
//			Prompt l,"L"
//			DoPrompt "hkl",h,k,l
//			if (V_flag)
//				return 1
//			endif
//		endif
//	endif
//
//	String lead = SelectString(WhichListItem("PolePopMenuProc",GetRTStackInfo(0))>=0, "","¥")
//	if (cursorAvailable)
//		printf "%sNewPoleFigure(\"%g %g %g\",rad=%g,useCursor=%g)\t\t//cursor = %d\r",lead,h,k,l,rad,useCursor,pcsr(A)
//	else
//		printf "%sNewPoleFigure(\"%g %g %g\",rad=%g)\r",lead,h,k,l,rad
//	endif
//	if (numtype(h+k+l))
//		return 1
//	endif
//	if (!(rad>0) || !(rad<=90))
//		DoAlert 0,"rad is not in range (0,90¡]"
//		return 1
//	endif
//
//	Make/N=3/O/D/FREE pole
//	pole={h,k,l}
//	Wave poleXY = $MakePolePoints(gm,pole,rad=rad,useCursor=useCursor)
//	Make_GraphPoles(poleXY)
//	return 0
//End



Function/T MakePolePoints(gm,hkl,[rad,useCursor,visible,iwave])
	Wave gm			// list of reciprocal lattices
	Wave hkl			// 3-vector with hkl of pole
	Variable rad		// radius of color circle (degree)
	Variable useCursor
	Wave visible		// flag for each point indicating which points to use (1=use, 0=ignore)
	Wave iwave			// wave to use for setting intensity

	useCursor = ParamIsDefault(useCursor) ? 1 : useCursor
	Variable x0=0,y0=0
	if (ItemsInLIst(WinList("*",";","WIN:1"))>0 && useCursor)
		if (stringmatch(StringByKey("TNAME", CsrInfo(A)),"pole*XY"))
			x0 = hcsr(A)
			y0 = vcsr(A)
		endif
	endif
	rad = ParamIsDefault(rad) ? 22.5 : rad

	Make/N=3/FREE/D poleHat
	poleHat = hkl
	normalize(poleHat)

	Wave vecs = $makeSymEquivVectors(poleHat)	// list of symmetry equivalent poles to examine
	Variable Nsym = DimSize(vecs,0)
	Nsym = !(Nsym>0) ? 1 : Nsym					// always at least one, the identity

	// define the surface
	Variable ir2 = 1/sqrt(2)
	Make/N=3/FREE/D normal={0,ir2,-ir2}		// the usual surface normal
	Make/N=3/FREE/D roll={0,-ir2,-ir2}			//   associated roll direction
	Make/N=3/FREE/D tilt={1,0,0}				//   associated tilt direction

	Variable Ng=DimSize(gm,2)
	Make/N=(3,3)/FREE/D gmi
	Make/N=3/O/D/FREE vec3
	Make/N=2/FREE/D vec2
	String str = hkl2str(round(hkl[0]),round(hkl[1]),round(hkl[2]))
	String wnote = ReplaceStringByKey("waveClass","","poleXYpoints","=")
	wnote = ReplaceStringByKey("hklStr",wnote,str,"=")
	wnote = ReplaceNumberByKey("rad",wnote,rad,"=")
	if (x0!=0 || y0!=0)
		wnote = ReplaceNumberByKey("x0",wnote,x0,"=")
		wnote = ReplaceNumberByKey("y0",wnote,y0,"=")
	endif
	if (WaveExists(visible))
		wnote = ReplaceStringByKey("visible",wnote,GetWavesDataFolder(visible,2),"=")
	endif
	str = ReplaceString(" ",str,"")
	str = ReplaceString("-",str,"")
	str = ReplaceString(",",str,"")
	String XYname = "pole"+str+"XY"
	Make/N=(Nsym*Ng,3)/O $XYname=NaN			// more than enough (will never get more than Nsym )
	Wave poleXY = $XYname

	String XYrgbname = "pole"+str+"XYrgb"
	Make/N=(Nsym*Ng,3)/O/W/U $XYrgbname=0	// more than enough (will never get more than Nsym )
	Wave poleXYrgb = $XYrgbname
	Make/N=3/FREE rgb

	Make/N=(Ng,3)/O/W/U rotRGBpole=0	// more than enough (will never get more than Nsym )
	Variable r = radiusOnPoleFigure(x0,y0,rad)	// radius on polefigure plot for saturated color
	Variable dist2, dist2min,xmin,ymin, iadd
	Variable i,N, dot, dist, m
	for (i=0,N=0;i<Ng;i+=1)						// loop over indexed points
		if (WaveExists(visible))
			if (!visible[i])							// skip all points that are NOT-visible
				continue
			endif
		endif

		gmi = gm[p][q][i]
		dist2min = Inf
		iadd = 0
		for (m=0;m<Nsym;m+=1)					// loop over symmetry equivalent poles
			vec3 = vecs[m][p]						// set vec3 to each of the symmetry equivalent directions
			MatrixOp/O/FREE vec3 = gmi x vec3
			normalize(vec3)
			dot = MatrixDot(vec3,normal)
			if (dot>=0)								// within 90¡
				vec2[0] = MatrixDot(tilt,vec3)
				vec2[1] = MatrixDot(roll,vec3)
				normalize(vec2)
				dist = sqrt(1-dot*dot) / (1+dot)
				vec2 *= dist
				poleXY[N][0] = vec2[0]
				poleXY[N][1] = vec2[1]
				poleXY[N][2] = i

				dist2 = (vec2[0]-x0)^2 + (vec2[1]-y0)^2

				if (dist2<dist2min)
					dist2min = dist2
					xmin = vec2[0]
					ymin = vec2[1]
				endif
				iadd += 1
				N += 1
			endif
		endfor
		poleXY2rgb(xmin,ymin,x0,y0,r,rgb)		// get best color, and add it
		for (m=N-iadd;m<N;m+=1)
			poleXYrgb[m][] = rgb[q]
		endfor
		rotRGBpole[i][] = rgb[q]					// save for other plots
	endfor
	printf "for (%s), N = %d\r",hkl2str(round(hkl[0]),round(hkl[1]),round(hkl[2])),N
	Redimension/N=(N,3) poleXY					// resize to actual points
	Redimension/N=(N,3) poleXYrgb
	wnote = ReplaceNumberByKey("Npnts",wnote,N,"=")
	Note/K poleXY, wnote

	// make rotRGBpoleIntens, which is probably more suitable for plots
	//	Wave iwave=sumAboveThreshold			// wave for intensity scaling
	//	Wave iwave=Nindexed						// wave for intensity scaling
	if (WaveExists(iwave))
		// Variable hi=waveMax(iwave)
		Variable hi=MedianOfWave(iwave,f=0.95)
		String name=GetWavesDataFolder(rotRGBpole,2)+"Intens"
		Duplicate/O rotRGBpole, $name
		Wave rgbi=$name
		// rgbi = limit(rotRGBpole[p][q]*iwave[p]/hi,0,65535)
		rgbi = limit(rotRGBpole[p][q]*sqrt(iwave[p]/hi),0,65535)
		// rgbi = limit(rotRGBpole[p][q]*(abs(iwave[p]/hi))^0.7,0,65535)
	endif

	if (exists("xyz")==1)							// make RGBA for use with a gizmo
		Duplicate/O rotRGBpole, rotRGBApole
		Redimension/S/N=(-1,4) rotRGBApole
		rotRGBApole = rotRGBpole/65535
		rotRGBApole[][3] = 0.5
		if (WaveExists(rgbi))
			Duplicate/O rotRGBpole, rotRGBApoleIntens
			Redimension/S/N=(-1,4) rotRGBApoleIntens
			rotRGBApoleIntens = rgbi/65535
			rotRGBApoleIntens[][3] = 0.5
		endif
	endif

	KillWaves vecs
	return GetWavesDataFolder(poleXY,2)
End
//
Static Function MedianOfWave(w, [x1, x2,f])					// Returns median value of wave w
	Wave w
	Variable x1, x2											// range of interest
	Variable f
	if (!WaveExists(w))
		return NaN
	endif
	Variable N=numpnts(w)
	Variable xlo=DimOffset(w,0),  xhi=DimOffset(w,0)+(N-1)*DimDelta(w,0)
	x1 = ParamIsDefault(x1) ? xlo : x1
	x2 = ParamIsDefault(x2) ? xhi : x2
	x1 = x1<xlo ? xlo : x1
	x2 = x2>xhi ? xhi : x2
	f = ParamIsDefault(f) ? 0.5 : f
	f = f>0 ? f : 0.5

	Variable result
	if (x1<x2 && numtype(x1+x2)==0)
		Duplicate/R=(x1,x2)/FREE w, tempMedianWave	// Make a clone of wave
	else
		Duplicate/FREE w, tempMedianWave					// Make a clone of wave
	endif
	Sort tempMedianWave, tempMedianWave					// Sort clone
	SetScale/P x 0,1,tempMedianWave
	result = tempMedianWave((numpnts(tempMedianWave)-1)*f)
	return result
End


Static Function radiusOnPoleFigure(x0,y0,rad)
	Variable x0,y0				// origin on pole figure
	Variable rad				// radius of circle around (x0,y0) (degree)

//  this simple equation was the old way of doing it.  It looks wrong.
//	Variable r = sqrt(1-(x0^2 + y0^2))*sin(rad*PI/180)	// radius on polefigure plot for saturated color
	Variable r2d = sqrt(x0^2 + y0^2)				// distance from center of (0,0)
	Variable phi=2*atan(1/r2d), dphi=rad*PI/180// azimuthal angle from pole to (x0,y0)
	phi = dphi<phi ? phi-dphi : phi+dphi			// aximuthal angle from (x0,y0) to dphi away
	Variable r=abs(r2d-sin(phi)/(1-cos(phi)))	// radius on pole figure corresponding to rad
	return r
End


Function MakColorHexagon()
	Variable Ni=512, Nj=512
	Make/N=(Ni,Nj,3)/O/U/W colorHexagon=2000
	SetScale/I x -1,1,"", colorHexagon
	SetScale/I y -1,1,"", colorHexagon
	Variable dx=DimDelta(colorHexagon,0), x0 = DimOffset(colorHexagon,0), px
	Variable dy=DimDelta(colorHexagon,1), y0 = DimOffset(colorHexagon,0), py

	Make/N=3/D/FREE rgb
	Variable i,j
	for (j=0;j<Nj;j+=1)
		py = j*dy + y0
		for (i=0;i<Ni;i+=1)
			px = i*dx + x0
			if ((px*px+py*py)>1)
				colorHexagon[i][j][] = 40000
			else
				xy2saturatedColors(px,py,1,65535,rgb)
				// poleXY2rgb(px,py,0,0,1,rgb)
				colorHexagon[i][j][] = rgb[r]
			endif
		endfor
	endfor

	String gName = StringFromList(0,FindGraphsWithWave(colorHexagon))
 	if (strlen(gName))
		DoWindow/F $gName
		return 0
 	else
		Display /W=(1651,44,1902,295)/K=1 
		AppendImage colorHexagon
		ModifyImage colorHexagon ctab= {*,*,Grays,0}
		ModifyGraph margin(left)=1,margin(bottom)=1,margin(top)=1,margin(right)=1,width={Aspect,1}
		ModifyGraph tick=3, mirror=2, noLabel=2, standoff=0, axThick=0
		SetDrawLayer UserFront
		SetDrawEnv xcoord= bottom,ycoord= left,fillpat= 0
		DrawOval -1,1,1,-1
		SetDrawEnv xcoord= bottom,ycoord= left,dash= 1
		DrawLine -1,0,1,0
		SetDrawEnv xcoord= bottom,ycoord= left,dash= 1
		DrawLine 0,1,0,-1
		DoAlert 1,"Add Direction Names?"
		if (V_flag==1)
			SetDrawEnv fsize= 18,textxjust= 1,textyjust= 1
			DrawText 0.92,0.5,"X"
			SetDrawEnv fsize= 18,textxjust= 1,textyjust= 1
			DrawText 0.29,0.14,"Y"
			SetDrawEnv fsize= 18,textrgb= (65535,65535,65535),textxjust= 1,textyjust= 1
			DrawText 0.29,0.86,"Z"
			SetDrawEnv fsize= 18,textxjust= 1,textyjust= 1
			DrawText 0.71,0.14,"-Z"
			SetDrawEnv fsize= 18,textxjust= 1,textyjust= 1
			DrawText 0.71,0.86,"-Y"
			SetDrawEnv fsize= 18,textxjust= 1,textyjust= 1
			DrawText 0.08,0.5,"-X"
			SetDrawEnv textxjust= 1,textyjust= 1
			DrawText 0.5,0.5,"(111)"
		endif
	endif
	return 0
End
//Window GraphColorHex() : Graph
//	PauseUpdate; Silent 1		// building window...
//	Display /W=(1515,48,1924,426)
//	AppendImage colorHexagon
////	AppendImage :area1:colorHexagon
//	ModifyImage colorHexagon ctab= {*,*,Grays,0}
//	ModifyGraph mirror=2
//EndMacro


Function poleXY2rgb(px,py,x0,y0,rad,rgb)
	Variable px,py							// x,y location on the pole figure
	Variable x0,y0							// center and range of color scan (x0,y0 is center point on pole figure)
	Variable rad								// color saturates at distance rad
	Wave rgb
	Variable dx,dy, hue
	dx = px-x0
	dy = py - y0
	if (dx==0 && dy==0)
		rgb = 1
		return 0
	endif
	hue = xy2saturatedColors(dx,dy,rad,65535,rgb)
	return hue
End
//Function poleXY2rgb(px,py,x0,y0,rad,rgb)
//	Variable px,py										// x,y location on the pole figure
//	Variable x0,y0,rad									// center and range of color scan (x0,y0 is center point on pole figure), color saturate at distance rad
//	Wave rgb
//
//	Variable dx,dy
//	dx = px-x0
//	dy = py - y0
//	if (dx==0 && dy==0)
//		rgb = 1
//		return 0
//	endif
//
//	Variable r = sqrt(dx*dx + dy*dy)
//	r = min(1,r/rad)							// r/rad, but does not exceed 1
//
//	Variable hue = atan2(dy,dx)*180/PI
//	hue += hue<0 ? 360 : 0
//	Variable huei = mod(floor(hue/60),6)
//
//	Variable pp,qq,tt, f
//	f = hue/60 - floor(hue/60)
//	pp = (1-r)
//	qq = 1-f*r
//	tt = 1-(1-f)*r
//	if (huei==0)
//		rgb = {1,tt,pp}
//	elseif(huei==1)
//		rgb = {qq,1,pp}
//	elseif(huei==2)
//		rgb = {pp,1,tt}
//	elseif(huei==3)
//		rgb = {pp,qq,1}
//	elseif(huei==4)
//		rgb = {tt,pp,1}
//	elseif(huei==5)
//		rgb = {1,pp,qq}
//	endif
////print f,hue,huei,"                    ",pp,qq,tt
////print hue,huei,"      ",f,"                    ",pp,qq,tt,"      ",rgb[0],rgb[1],rgb[2],"        f=",f,"         s=",r
//	rgb *= 65535
//	rgb = round(rgb)
//	rgb = limit(rgb,0,65535)
//	return hue*180/PI				// return hue angle (degree)
//End



Static Function/WAVE CubicTriangleColors(hkl,[rgbMax])	// returns rgb wave, used for Inverse Pole Figure Colors
	//	Makes rgb, for the hkl on the inverse pole figure triangle
	//	111=Blue, 001=Red, 101=Green,   all colors are returned saturated so middle is white
	//	returns Black for the (000), and NaN for any hkl that contain NaN
	//	Note, this only works for CUBIC!.  This returns nonsense for anything else.
	Wave hkl
	Variable rgbMax								// scale for RGB values, default=1, usually 1 or 65535
	rgbMax = ParamIsDefault(rgbMax) ? NaN : rgbMax
	rgbMax = rgbMax>0 && numtype(rgbMax)==0 ? rgbMax : 1

	Make/N=3/FREE/D rgb=NaN					// set rgb to scale [0,rgbMax]
	String wNote = ReplaceStringByKey("hkl","waveClass=rgb;",vec2str(hkl,bare=1),"=")

	if (numtype(sum(hkl)))
		Note/K rgb,wNote
		return rgb
	endif

	Duplicate/FREE hkl, vec
	Variable len = norm(vec)
	if (len==0)
		Note/K rgb,ReplaceNumberByKey("delta",wNote,0,"=")
		rgb = 0
		return rgb
	endif
	vec /= len										// input vector now normalized

	// reduce hkl into the triangle
	vec = abs(vec)									// so now vec is all positive and in ascending order
	Sort vec,vec									// and vec[0]²vec[1]²vec[2]

	Make/N=(3,3)/FREE/D poles
	poles[0][0] = {0,0,1}
	poles[0][1] = {0,1,1}
	poles[0][2] = {1,1,1}
	poles[][1] /= sqrt(2)						// normalized vectors in (001), (011), & (111) directions
	poles[][2] /= sqrt(3)

	MatrixOP/FREE coefs = Inv(poles) x vec
	coefs = coefs < 1e-12 ? 0 : coefs

	MatrixOP/FREE delta = sum(abs((poles x coefs) - vec))	// the error, may be caused by trimming coefs
	Note/K rgb,ReplaceNumberByKey("delta",wNote,delta[0],"=")

	rgb = coefs * rgbMax/WaveMax(coefs)
	return rgb
End
//Function testMany(N)
//	Variable N					// 1e4 is a good value (1e4 takes about 1 sec)
//
//	Make/N=3/D/FREE hkl
//	String hklMax=""
//	Variable i, h,k,l, delta, dmax=-Inf
//	Variable microSec=stopMSTimer(-2)
//	for (i=0;i<N;i+=1)
//		h = enoise(10)
//		k = enoise(10)
//		l = enoise(10)
//		hkl={h,k,l}
//		hkl = round(hkl)
//		Wave rgb = CubicTriangleColors(hkl,rgbMax=65532)
//		delta = NumberByKey("delta",note(rgb),"=")
//		if (delta>dmax)
//			hklMax = vec2str(hkl)
//			dmax = delta
//		endif
//	endfor
//	printf "delta max = %g,   at hkl = %s   (took %.2f sec)\r",dmax,hklMax, (stopMSTimer(-2)-microSec)*1e-6
//End
////
//Function test1(h,k,l)
//	Variable h,k,l
//
//	Make/N=3/D/FREE hkl={h,k,l}
//	Wave rgb = CubicTriangleColors(hkl,rgbMax=65532)
//
//	print "rgb = ",vec2str(rgb), "delta = ",NumberByKey("delta",note(rgb),"=")
//	return NumberByKey("delta",note(rgb),"=")
//End
//
Function MakeColorTriangle([N,gray])			// show the color triangle that you can compute with CubicTriangleColors()
	Variable N
	Variable gray
	N = ParamIsDefault(N) ? NaN : N
	N = numtype(N) || N<=1 ? 200 : N
	gray = ParamIsDefault(gray) ? NaN : gray
	gray = numtype(gray)==2 || gray<0 ? 50000 : gray
	gray = limit(gray,0,65535)
	Variable Nx=N,Ny=N

	Variable rmax = 1-(2-sqrt(2))
	Variable ymax = 1/(sqrt(3)+1)
	Make/N=(Nx,Ny,3)/FREE hkls
	SetScale/I x 0,rmax,"", hkls
	SetScale/I y 0,ymax,"", hkls

	hkls[][][0] = 2*x / (1 + x^2 + y^2)
	hkls[][][1] = 2*y / (1 + x^2 + y^2)
	hkls[][][2] = (1 - x^2 - y^2) / (1 + x^2 + y^2)
	hkls[][][] = y>x ? NaN : hkls[p][q][r]

	Make/N=3/D/FREE abc
	Make/N=(3,3)/FREE/D poles
	poles[0][0] = {0,0,1}
	poles[0][1] = {1,0,1}
	poles[0][2] = {1,1,1}
	Make/N=(Nx,Ny,3)/O/W/U ColorTriangleStereoCubic=gray
	SetScale/I x 0,rmax,"", ColorTriangleStereoCubic
	SetScale/I y 0,ymax,"", ColorTriangleStereoCubic
	Make/N=3/D/FREE hkl
	Variable i,j
	for (j=0;j<Ny;j+=1)
		for (i=0;i<Nx;i+=1)
			hkl = hkls[i][j][p]
			if (numtype(sum(hkl)))
				continue
			endif
			MatrixOP/FREE/O abc = Inv(poles) x hkl
			if (abc[0]>=0)
				Wave rgb = multiIndex#CubicTriangleColors(hkl,rgbMax=65535)
				ColorTriangleStereoCubic[i][j][] = rgb[r]
			endif
		endfor
	endfor
	GraphColorTriangle(ColorTriangleStereoCubic)
End
//
Function GraphColorTriangle(image3)
	Wave image3
	if (!WaveExists(image3))				// wave does not exist
		return 1
	elseif (WaveDims(image3)!=3 || DimSize(image3,2)!=3)
		return 1							// wave must have 3 layers (r,g,b)
	endif

	String win=StringFromList(0,FindGraphsWithWave(image3))
	if (strlen(win))						// already on a grapgh, bring to front
		DoWinDow/F $win
		return 0
	endif
	Variable aspect = DimDelta(image3,1)*DimSize(image3,1)/(DimDelta(image3,0)*DimSize(image3,0))

	Display /W=(262,149,770,502)/K=1
	AppendImage image3
	ModifyImage $NameOfWave(image3) ctab= {*,*,Grays,0}
	ModifyGraph height={Aspect,aspect}, mirror=2, noLabel=2, tick=3, margin=5
	TextBox/C/N=text001/F=0/S=3/G=(0,65535,65535)/B=1/A=LB/X=6.57/Y=1.03 "\\Z18(001)"
	TextBox/C/N=text101/F=0/S=3/G=(65535,0,65535)/B=1/A=RB/X=1.95/Y=1.03 "\\Z18(101)"
	TextBox/C/N=text111/F=0/S=3/G=(65535,65535,0)/B=1/X=10/Y=10"\\Z18(111)"
	return 0
End




Function MakepoleImage(poleXY)
	Wave poleXY
//	Make/N=(128,128)/U/W/O poleXYimage=0
//	Make/N=(512,512)/O poleXYimage=0
	Make/N=(301,301)/O poleXYimage=0
	SetScale/I x -1,1,"", poleXYimage
	SetScale/I y -1,1,"", poleXYimage
	Variable dx=DimDelta(poleXYimage,0), x0 = DimOffset(poleXYimage,0)
	Variable dy=DimDelta(poleXYimage,1), y0 = DimOffset(poleXYimage,0)
	Variable i,j,m,N=DimSize(poleXY,0)
	for (m=0;m<N;m+=1)
		i = (poleXY[m][0]-x0)/dx
		j = (poleXY[m][1]-y0)/dy
		poleXYimage[i][j] += 1
	endfor
	poleXYimage = sqrt(sqrt(poleXYimage))
	poleXYimage = numtype(poleXYimage) ? NaN : poleXYimage
End


Function Make_GraphPoles(poleXY)				// makes (or fixes up) a pole figure graph
	Wave poleXY

	if (!WaveExists(poleXY))
		String poleName=""
		Prompt poleName,"pole figure to show",popup,WaveListClass("poleXYpoints","*","DIMS:2;MINCOLS:2")
		DoPrompt "Pole Figure",poleName
		if (V_flag)
			return 1
		endif
		Wave poleXY = $poleName
	endif
	if (!WaveExists(poleXY))
		return 1
	endif
	if (!WaveInClass(poleXY,"poleXYpoints"))
		return 1
	endif
	String gName = StringFromList(0,FindGraphsWithWave(poleXY))
	String wnote=note(poleXY)
	String hkl = "("+StringByKey("hklStr",wnote,"=")+")"
 	if (strlen(gName))
		DoWindow/F $gName
	else
	 	Display /W=(257,198,765,706)/K=1 poleXY[*][1] vs poleXY[*][0]
	 	ModifyGraph margin(left)=1,margin(bottom)=1,margin(top)=1,margin(right)=1,gfMult=115
		ModifyGraph width={Aspect,1}, mode=2, lSize=2, tick=3, zero=2, mirror=2, noLabel=2, standoff=0, axThick=0, zeroThick=1
		String name = GetWavesDataFolder(poleXY,2)+"rgb"
		Wave rgbWave = $name
		if (WaveExists(rgbWave))
			ModifyGraph gbRGB=(40000,40000,40000)
			ModifyGraph zColor($NameOfWave(poleXY))={rgbWave,*,*,directRGB}
		endif
		SetAxis left -1,1
		SetAxis bottom -1,1
		TextBox/C/N=textLabel/F=0/B=1/X=1/Y=1 "\\JR"+hkl
		ShowInfo

		PopupMenu polePopup,pos={2,2},size={52,20},proc=multiIndex#PolePopMenuProc,title="pole"
//		PopupMenu polePopup,mode=0,value= #"\"Re-Calc...;Full Size;scale bar (toggle);grain no.\""
		PopupMenu polePopup,mode=0,value= #"\"Re-Calc...;Full Size;scale bar (toggle);Grain no.;Mark a Grain\""
	endif
	TextBox/C/N=textLabel "\\JR"+hkl
	SetDrawLayer/K UserFront
	SetDrawEnv xcoord= bottom,ycoord= left,fillpat= 0
	DrawOval -1,1,1,-1
	Variable x0=NumberByKey("x0",wnote,"=")		// re-position cursor if appropriate
	Variable y0=NumberByKey("y0",wnote,"=")
	Variable rad=NumberByKey("rad",wnote,"=")
	if (numtype(x0+y0)==0)
		MatrixOp/O/FREE dist2 = (col(poleXY,0) - x0)*(col(poleXY,0) - x0) + (col(poleXY,1) - y0)*(col(poleXY,1) - y0)
		WaveStats/M=1/Q dist2
		if (numtype(V_minloc)+rad)
			Cursor/P A $NameOfWave(poleXY) V_minloc
			SetDrawLayer UserFront
			SetDrawEnv xcoord= bottom,ycoord= left,fillpat= 0, dash= 11
			Variable r = radiusOnPoleFigure(x0,y0,rad)	// radius on polefigure plot for saturated color
			DrawOval x0-r,y0-r,x0+r,y0+r
			DrawMarker(x0,y0,r*0.6,r*0.6,"cross gap")
			TextBox/C/N=textLabel "\\JR"+hkl+"\r\\Zr075Color Circle = "+num2str(rad)+"¡\\M"
		endif
	endif
	Variable N=NumberByKey("Npnts",wnote,"=")
	if (numtype(N)==0)
		AppendText/N=textLabel "\\JR\\Zr075used "+num2istr(N)+" points\\M"
	endif
	DoUpdate
	return 0
End
//
Static Function PolePopMenuProc(pa) : PopupMenuControl
	STRUCT WMPopupAction &pa
	if (pa.eventCode != 2)
		return 0
	endif

	Variable grain
	if (stringmatch(pa.popStr,"Re-Calc..."))
		Wave pole = TraceNameToWaveRef(pa.win,StringFromLIst(0,TraceNameList(pa.win,";",3)))
		String wnote = note(pole)
		String hklStr=StringByKey("hklStr",wnote,"=")
		Variable rad=NumberByKey("rad",wnote,"=")
		NewPoleFigure(hklStr,rad=rad,useDialog=1)
	elseif (stringmatch(pa.popStr,"Full Size"))
		SetAxis/W=$pa.win left -1,1
		SetAxis/W=$pa.win bottom -1,1
	elseif (stringmatch(pa.popStr,"scale bar (toggle)"))
		Wave pole = TraceNameToWaveRef(pa.win,StringFromLIst(0,TraceNameList(pa.win,";",3)))
		Variable barOn = str2num(GetUserData(pa.win,"","scaleBarButton"))
		barOn = numtype(barOn) ? 0 : barOn
		SetWindow $pa.win userdata(scaleBarButton)=num2str(multiIndex#scaleBar(pole,!barOn))
	elseif (stringmatch(pa.popStr,"Grain no."))
		Wave pole = TraceNameToWaveRef(pa.win,StringFromLIst(0,TraceNameList(pa.win,";",3)))
		if (stringmatch(StringByKey("TNAME", CsrInfo(A)),"pole*XY"))
			grain = pole[pcsr(A)][2]
			if (grain>0)
				print "cursor A at grain ",grain
			endif
		endif
		if (stringmatch(StringByKey("TNAME", CsrInfo(B)),"pole*XY"))
			grain = pole[pcsr(B)][2]
			if (grain>0)
				print "cursor B at grain ",grain
			endif
		endif
	elseif (stringmatch(pa.popStr,"Mark a Grain"))
		Wave pole = TraceNameToWaveRef(pa.win,StringFromLIst(0,TraceNameList(pa.win,";",3)))
		Prompt grain,"grain no."
		DoPrompt "grain no.",grain
		if (V_flag==0)			// mark each point on pole figure coming from grain
			Variable i,none=1
			SetDrawLayer/W=$pa.win/K UserFront
			barOn = str2num(GetUserData(pa.win,"","scaleBarButton"))
			barOn = numtype(barOn) ? 0 : barOn
			SetWindow $pa.win userdata(scaleBarButton)=num2str(multiIndex#scaleBar(pole,barOn))
			for (i=0;i<DimSize(pole,0);i+=1)
				if (abs(pole[i][2]-grain) < 0.5)
					none = 0
					printf "for grain %d, draw X at pnt %d\r",grain,i
					DrawMarker(pole[i][0],pole[i][1],0.1,0.1,"X",color="white",dash=2,win=pa.win)
				endif
			endfor
			if (none && grain>=0)
				printf "for grain %d, no points found\r",grain
			endif
		endif
	endif

	return 0
End
//
Static Function FullScaleButtonProc(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	if (ba.eventCode!=2)
		return 0
	endif

	if (stringmatch(ba.ctrlName,"fullScaleButton"))
		SetAxis/W=$ba.win left -1,1
		SetAxis/W=$ba.win bottom -1,1
	elseif (stringmatch(ba.ctrlName,"recalcButton"))
		Wave pole = TraceNameToWaveRef(ba.win,StringFromLIst(0,TraceNameList(ba.win,";",3)))
		String wnote = note(pole)
		String hklStr=StringByKey("hklStr",wnote,"=")
		Variable rad=NumberByKey("rad",wnote,"=")
		NewPoleFigure(hklStr,rad=rad,useDialog=1)
	elseif (stringmatch(ba.ctrlName,"scaleBarButton"))
		Wave pole = TraceNameToWaveRef(ba.win,StringFromLIst(0,TraceNameList(ba.win,";",3)))
		Variable barOn = str2num(GetUserData(ba.win,"","scaleBarButton"))
		barOn = numtype(barOn) ? 0 : barOn
		SetWindow $ba.win userdata(scaleBarButton)=num2str(scaleBar(pole,!barOn))
	endif
	return 0
End
//
Static Function scaleBar(poleXY,barOn)
	Wave poleXY
	Variable barOn

	String wnote=note(poleXY)
	Variable rad=NumberByKey("rad",wnote,"=")
	Variable x0=NumberByKey("x0",wnote,"=")
	Variable y0=NumberByKey("y0",wnote,"=")
//	Variable r = sqrt(1-(x0^2 + y0^2))*sin(rad*PI/180)
	Variable r = radiusOnPoleFigure(x0,y0,rad)	// radius on polefigure plot for saturated color

	SetDrawLayer/K UserFront
	SetDrawEnv xcoord= bottom,ycoord= left,fillpat= 0
	DrawOval -1,1,1,-1

	SetDrawLayer UserFront
	SetDrawEnv xcoord= bottom,ycoord= left,fillpat= 0, dash= 11
	DrawOval x0-r,y0-r,x0+r,y0+r

	if (barOn && (r>0))
		SetDrawEnv xcoord=bottom,ycoord=prel,linethick=0.00,linefgc=(65535,65535,65535)
		DrawRect x0,0.9,x0+r,0.9+0.03
		SetDrawEnv xcoord= bottom,ycoord= prel,textxjust= 1,textyjust= 1
		DrawText x0+r/2,0.9+0.015,num2str(rad)+"¡"
		return 1
	endif
	return 0
End

//Static Function FullScaleButtonProc(ba) : ButtonControl
//	STRUCT WMButtonAction &ba
//	if (ba.eventCode!=2)
//		return 0
//	endif
//
//	if (stringmatch(ba.ctrlName,"fullScaleButton"))
//		SetAxis/W=$ba.win left -1,1
//		SetAxis/W=$ba.win bottom -1,1
//	elseif (stringmatch(ba.ctrlName,"recalcButton"))
//		Wave pole = TraceNameToWaveRef(ba.win,StringFromLIst(0,TraceNameList(ba.win,";",3)))
//		String wnote = note(pole)
//		String hklStr=StringByKey("hklStr",wnote,"=")
//		Variable rad=NumberByKey("rad",wnote,"=")
//		NewPoleFigure(hklStr,rad=rad,useDialog=1)
//	endif
//	return 0
//End


// returns all of the unique vectors that are symmetry equivalent to vecIn
Static Function/T makeSymEquivVectors(vecIn)
	Wave vecIn									// input vector of interest

	if (numpnts(vecIn)!=3)
		return ""
	endif
	Variable v2=norm(vecIn)^2
	if (v2<=0)
		return ""
	endif

	Make/N=(3,3)/FREE/D symOp
	Wave symOps=$StrVarOrDefault("root:Packages:Lattices:SymOps:SymmetryOpsPath","")	// wave with symmetry ops
	Variable Nops = 0							// number of symmetyry operations to check
	if (WaveExists(symOps))
		Nops = NumberByKey("Nproper",note(symOps),"=")
		Nops = numtype(Nops) ? DimSize(symOps,0) : Nops
	endif
	String vecName = NameOfWave(vecIn) + "_Syms"
	Make/N=(Nops,3)/O/D $vecName=NaN
	Wave vecs=$vecName						// sym equivalent vectors created here
	Make/N=3/O/D vec_makeSymEquivVectors, vecs_makeSymEquivVectors
	Wave vec=vec_makeSymEquivVectors, vs=vecs_makeSymEquivVectors

	Variable N=0								// number of symmetry equivalents found
	Variable i,m, dot
	for (i=-0;i<Nops;i+=1)
		symOp = symOps[i][p][q]				// one symmetry operation

		MatrixOp/O vec = symOp x vecIn
		for (m=0;m<N;m+=1)					// now compare against list of equivalents
			vs = vecs[m][p]
			dot = MatrixDot(vec,vs)
			if (abs(dot-v2)<(v2*1e-4))		// found a match
				break
			endif
		endfor
		if (m>=N)								// a new one, add it
			vecs[N][] = vec[q]
			N += 1
		endif
	endfor
	Redimension/N=(N,3) vecs
	KillWaves/Z vec_makeSymEquivVectors, vecs_makeSymEquivVectors
	return GetWavesDataFolder(vecs,2)
End


// This is one way to make the mask used in MakePolePoints()
Static Function/WAVE MakeVisibleMaskInGizmo(waveMatchStr,[resolution])
	String waveMatchStr
	Variable resolution				// used to give a bit of extra space when selecting points
	resolution = ParamIsDefault(resolution) ? NaN : resolution
	if (strlen(waveMatchStr)<1)
		waveMatchStr = "*"
	endif

	String objName = TopObjectTypeInGizmo("Scatter",rejectList="cubeCorners")	// returns top Scatter from top Gizmo
	if (strlen(objName)<1)
		return $""
	endif

	Execute "GetGizmo  objectItemExists="+objName
	NVAR flag=V_flag
	if (!flag)
		return $""
	endif

	Execute "GetGizmo objectList"
	SVAR S_gizmoObjectList=S_gizmoObjectList

	String item=""
	Variable i,N=ItemsInList(S_gizmoObjectList)

	for (i=0;i<N;i+=1)
		item = StringFromList(i,S_gizmoObjectList)
		if (strsearch(item,",name="+objName,Inf,3)>=0)
			break
		endif
		item = ""
	endfor
	if (strlen(item)<1)
		return $""
	endif
	item = TrimFrontBackWhiteSpace(StringFromList(0,item,","))
	String wName = StringByKey("AppendToGizmo Scatter", item,"=")
	if (!stringmatch(wName,waveMatchStr))
		DoAlert 0, "Could not find a Gizmo with a scatter wave matching '"+waveMatchStr+"'"
		return $""
	endif
	KillStrings/Z S_gizmoObjectList
	KillWaves/Z TW_gizmoObjectList

	Wave wxyz=$wName
	if (!waveExists(wxyz))
		return $""
	endif
	resolution = numtype(resolution) ? NumberByKey("resolution",note(wxyz),"=") : resolution
	resolution = numtype(resolution) ? 0.5 : resolution

	Execute "GetGizmo  userBoxLimits"		// get displayed box
	NVAR GizmoBoxXmin=GizmoBoxXmin, GizmoBoxXmax=GizmoBoxXmax
	NVAR GizmoBoxYmin=GizmoBoxYmin, GizmoBoxYmax=GizmoBoxYmax
	NVAR GizmoBoxZmin=GizmoBoxZmin, GizmoBoxZmax=GizmoBoxZmax
	Variable xlo=GizmoBoxXmin-resolution, xhi=GizmoBoxXmax+resolution
	Variable ylo=GizmoBoxYmin-resolution, yhi=GizmoBoxYmax+resolution
	Variable zlo=GizmoBoxZmin-resolution, zhi=GizmoBoxZmax+resolution
	KillVariables/Z GizmoBoxXmin,GizmoBoxXmax,GizmoBoxYmin,GizmoBoxYmax,GizmoBoxZmin,GizmoBoxZmax

	String vName = wName+"Visible"
	N=DimSize(wxyz,0)

	MatrixOP/FREE xi = col(wxyz,0)
	MatrixOP/FREE yi = col(wxyz,1)
	MatrixOP/FREE zi = col(wxyz,2)
	MatrixOP/O $vName = equal(clip(xi,xlo,xhi),xi) && equal(clip(yi,ylo,yhi),yi) && equal(clip(zi,zlo,zhi),zi)
	Wave visible=$vName

	xlo += resolution	;	xhi -= resolution
	ylo += resolution	;	yhi -= resolution
	zlo += resolution	;	zhi -= resolution
	String wnote="waveClass=GizmoVisible"
	wnote = ReplaceStringByKey("source",wnote,wName,"=")
	wnote = ReplaceNumberByKey("xlo",wnote,xlo,"=")
	wnote = ReplaceNumberByKey("xhi",wnote,xhi,"=")
	wnote = ReplaceNumberByKey("ylo",wnote,ylo,"=")
	wnote = ReplaceNumberByKey("yhi",wnote,yhi,"=")
	wnote = ReplaceNumberByKey("zlo",wnote,zlo,"=")
	wnote = ReplaceNumberByKey("zhi",wnote,zhi,"=")
	Note/K visible, wnote
	if (strlen(GetRTStackInfo(2))<1)
		Variable in=sum(visible)
		printf "wName = '%s'\r",wName
		printf "Gizmo Box is:   [%g,  %g],  [%g,  %g],  [%g,  %g]µm,   size is %g x %g x %g µm\r",xlo,xhi, ylo,yhi, zlo,zhi,xhi-xlo,yhi-ylo,zhi-zlo
		printf "vName = '%s',  %.0f points inside box, %.0f outside,  total=%.0f\r",vName,in,N-in,N
	endif
	return visible
End

Static Function/T TopObjectTypeInGizmo(objectType,[rejectList])		// returns name of top object type in the top Gizmo
	String objectType							// type of object, e.g. "Scatter"
	String rejectList							// reject any names in this semicoln separated list, usually "cubeCorners;"
	if (ParamIsDefault(rejectList))
		rejectList = ""
	endif
	if (strlen(objectType)<1)
		return ""
	endif
	Execute "GetGizmo objectList"
	SVAR S_gizmoObjectList=S_gizmoObjectList
	String objName
	Variable i0=0, i1
	do
		objName = ""
		i0 = strsearch(S_gizmoObjectList,"AppendToGizmo "+objectType+"=",i0,2)
		if (i0>0)
			i0 = strsearch(S_gizmoObjectList,",name=",i0)
			i1 = strsearch(S_gizmoObjectList,";",i0+6)
			if (i0>0 && i1>0)
				objName = S_gizmoObjectList[i0+6,i1-1]
			endif
			i0 += 10
		endif
	while(WhichListItem(objName,rejectList)>=0 && i0>0)
	KillWaves/Z TW_gizmoObjectList
	KillStrings/Z S_gizmoObjectList
	return objName
End

// ******************************  this section needs work ******************************
//
//Function moveBallHook(s)
//	STRUCT WMGizmoHookStruct &s
//	if (!stringmatch(s.eventName,"mouseMoved"))
//		return 0
//	elseif (GetKeyState(0)<1)
//		return 0
//	endif
//
////	print s.eventName,"  ",s.winName,"   ",GetKeyState(0),"   ",s
////	print " "
//
//	Variable mouseY = s.mouseY
//	NVAR mouseYlast=root:Gizmo_mouseYlast
//	if (!NVAR_Exists(mouseYlast))
//		Variable/G root:Gizmo_mouseYlast=mouseY
//		print "created"
//		print " "
//		return 0
//	endif
//
//	Variable deltaCursor = mouseY-mouseYlast
//	deltaCursor = abs(deltaCursor)>10 ? 0 : deltaCursor
//
//	mouseYlast = mouseY
//	Variable modifier = GetKeyState(0)
//	if (modifier&1)						// command key x 10
//		deltaCursor *= 10
//	endif
//	if (modifier&2)						// option key X 100
//		deltaCursor *= 100
//	endif
//
//
//	NVAR markerApos=root:Gizmo_markerApos
//	if (!NVAR_Exists(markerApos))
//		DoAlert 0,"gizmoScatterMarker pos does not exist, is it on Gizmo?"
//		return 0
//	endif
//	markerApos += deltaCursor
//	markerApos = max(markerApos,0)
//
////	print "deltaCursor=",deltaCursor,"    ",GetKeyState(0),"     cursor ->",markerApos
//
//	PlaceCursors(markerApos)
//	return 0
//End
////
////		  mouseMoved      Gizmo0       4       
////		STRUCT WMGizmoHookStruct
////		 version: 1
////		 winName[32]: Gizmo0
////		 eventName[32]: mouseMoved
////		 width: 687
////		 height: 656
////		 mouseX: 106
////		 mouseY: 276
////		 xmin: -5561.07
////		 xmax: -5356.93
////		 ymin: 6372.52
////		 ymax: 6576.66
////		 zmin: -1296.76
////		 zmax: -1092.62
////		 eulerA: 117.733
////		 eulerB: 17.355
////		 eulerC: 9.92809
////		 wheelDx: 1.04311e-304
////		 wheelDy: 1.27322e-313
//
//Function PlaceCursors(m)
//	Variable m
//
//	Wave gizmoScatterMarker=gizmoScatterMarker, xyz=xyz
//	gizmoScatterMarker[0][] = xyz[m][q]
//
//	Wave poleXY=pole001XY
//	String char
//	Variable N=DimSize(poleXY,0), i,j
//	for (i=0,j=0; i<N;i +=1)
//		if (poleXY[i][2]==m)
//			char = num2char(char2num("A")+j)
//			Cursor/P/S=2 $char $NameOfWave(poleXY) i
//			j += 1
//		endif
//	endfor
//	DoUpdate
//	Execute "ModifyGizmo update = 2"
//End
//
// ******************************************************************************

// ******************************  End of Pole Figure Colors *******************************
// *****************************************************************************************



// *****************************************************************************************
// *********************************  Start of Histograms **********************************

Static Function HistogramOf3dArray(w3d,[,normType,visible])
	Wave w3d
	String normType			// must be "peak" or "area" or "integral" or "none"  (integral and area are the same)
	Wave visible				// flag for each point indicating which points to use (1=use, 0=ignore)
	if (ParamIsDefault(visible))
		Wave visible=$""
	endif
	String visList=WaveListClass("GizmoVisible","*","MINROWS:2,MAXCOLS:0"), visName=""
	Variable askVisible = ItemsInList(visList)>0 && !WaveExists(visible)
	if (WaveExists(visible))
		visName = NameOfWave(visible)
	endif
	normType = SelectString(ParamIsDefault(normType),normType,"")
	normType = LowerStr(normType)
	Variable askNorm = !stringmatch(normType,"peak") && !stringmatch(normType,"area") && !stringmatch(normType,"integral") && !stringmatch(normType,"*none*")
	if (!WaveExists(w3d) || askVisible || askNorm)
		String wList = WaveListClass("Interpolated3dArrays","*",""), wName
		Prompt wName,"3D wave",popup,wList
		Prompt normType,"Histogram Normalization",popup,"Peak;Integral;None"
		if (askVisible)
			Prompt visName,"Visible Mask from Gizmo",popup, visList+" none "
			DoPrompt "3D Wave",wName,normType,visName
			visName = SelectString(stringmatch(visName," none "),visName,"")
			Wave visible=$visName
		else
			DoPrompt "3D Wave",wName,normType
		endif
		if (V_flag)
			return 1
		endif
		normType = LowerStr(normType)
		Wave w3d=$wName
		if (WaveExists(visible))
			printf "HistogramOf3dArray(%s,normType=\"%s\",visible=%s)\r",NameOfWave(w3d),normType,NameOfWave(visible)
		else
			printf "HistogramOf3dArray(%s,normType=\"%s\")\r",NameOfWave(w3d),normType
		endif
	endif
	if (!WaveExists(w3d))
		return 1
	endif
	if (!stringmatch(normType,"peak") && !stringmatch(normType,"area") && !stringmatch(normType,"integral") && !stringmatch(normType,"*none*"))
		print "normType must be 'peak' or 'area' or 'integral'"
		return 1
	endif
	String Vnote=""
	if (WaveExists(visible))
		Vnote = note(visible)
	endif

	Variable x0=DimOffset(w3d,0), dx=DimDelta(w3d,0)
	Variable y0=DimOffset(w3d,1), dy=DimDelta(w3d,1)
	Variable z0=DimOffset(w3d,2), dz=DimDelta(w3d,2)
	Variable xlo=NumberByKey("xlo", Vnote,"="), xhi=NumberByKey("xhi", Vnote,"=")
	Variable ylo=NumberByKey("ylo", Vnote,"="), yhi=NumberByKey("yhi", Vnote,"=")
	Variable zlo=NumberByKey("zlo", Vnote,"="), zhi=NumberByKey("zhi", Vnote,"=")
	xlo = numtype(xlo) ? -Inf : xlo	;	xhi = numtype(xhi) ? Inf : xhi
	ylo = numtype(ylo) ? -Inf : ylo	;	yhi = numtype(yhi) ? Inf : yhi
	zlo = numtype(zlo) ? -Inf : zlo	;	zhi = numtype(zhi) ? Inf : zhi

	WaveStats/M=1/Q w3d
	Variable Nmax=V_npnts							// maximum possible number of points in the histogram
	Make/N=(Nmax)/FREE values=NaN				// will hold values to be histogramed
	Variable Nx=DimSize(w3d,0),Ny=DimSize(w3d,1), Nz=DimSize(w3d,2), xStart,yStart,zStart
	Nx = numtype(xhi) ? Nx : min(Nx,ceil((xhi-x0)/dx))	// sets volume to loop over
	xStart = numtype(xlo) ? 0 : max(0,floor((xlo-x0)/dx))
	Ny = numtype(yhi) ? Ny : min(Ny,ceil((yhi-y0)/dy))
	yStart = numtype(ylo) ? 0 : max(0,floor((ylo-y0)/dy))
	Nz = numtype(zhi) ? Nz : min(Nz,ceil((zhi-z0)/dz))
	zStart = numtype(zlo) ? 0 : max(0,floor((zlo-z0)/dz))
	Variable ix,iy,iz, N, val
	for(iz=zStart,N=0;iz<Nz;iz+=1)				// loop over all points in w3d
		for(iy=yStart;iy<Ny;iy+=1)
			for(ix=xStart;ix<Nx;ix+=1)
				val = w3d[ix][iy][iz]
				if (numtype(val)==0)				// only histogram good values
					values[N] = val
					N += 1
				endif
			endfor
		endfor
	endfor
	Redimension/N=(N) values

	String units=""
	if (stringmatch(NameOfWave(w3d),"R3*"))
		values = 2*atan(values) * 180/PI					// convert Rodriques vectors to degrees
		units = "¡"
	elseif (stringmatch(NameOfWave(w3d),"GND*"))
		units = "dislocations/cm\S2\M"
	endif

	String hName = NameOfWave(w3d)+"_Hist"
	Variable Nhist=max(round(numpnts(values)/350),51)	// number of points in histogram
	Make/N=(Nhist)/O $hName
	Wave Hout = $hName

	if (WaveMin(values)<0)								// when bipolar, reset average to zero
		Variable avg=faverage(values)
		values -= avg
	endif

	WaveStats/M=1/Q values
	Variable maxmax=max(abs(V_min),abs(V_max)), bipolar=(V_min<0)
	SetScale/I x, (bipolar ? -maxmax : 0), maxmax,units, Hout
	Histogram/B=2 values, Hout
	Variable PeakValue=WaveMax(Hout)					// normalize to a peak value of 1
	Variable Harea=area(Hout)
	Variable Hnorm=1
	if (stringmatch(normType,"peak"))
		Hnorm = PeakValue
	elseif (stringmatch(normType,"integral") || stringmatch(normType,"area"))
		Hnorm = Harea
	endif
	Hout /= Hnorm

	String Hnote=ReplaceStringByKey("waveClass",note(w3d),"HistogramFrom3d","=")
	Hnote = SelectString(numtype(xlo), ReplaceNumberByKey("xlo",Hnote,xlo,"="), Hnote)
	Hnote = SelectString(numtype(xhi), ReplaceNumberByKey("xhi",Hnote,xhi,"="), Hnote)
	Hnote = SelectString(numtype(ylo), ReplaceNumberByKey("ylo",Hnote,ylo,"="), Hnote)
	Hnote = SelectString(numtype(yhi), ReplaceNumberByKey("yhi",Hnote,yhi,"="), Hnote)
	Hnote = SelectString(numtype(zlo), ReplaceNumberByKey("zlo",Hnote,zlo,"="), Hnote)
	Hnote = SelectString(numtype(zhi), ReplaceNumberByKey("zhi",Hnote,zhi,"="), Hnote)
	Hnote = ReplaceNumberByKey("PeakValue",Hnote,PeakValue,"=")
	Hnote = ReplaceNumberByKey("Integral",Hnote,Harea,"=")
	Hnote = ReplaceNumberByKey("norm",Hnote,Hnorm,"=")
	Hnote = ReplaceNumberByKey("pointsInHistogram",Hnote,N,"=")
	normType[0,0] = UpperStr(normType[0])
	Hnote = ReplaceStringByKey("normType",Hnote,normType,"=")
	Note/K Hout, Hnote

	printf "Histogramed values from '%s' into '%s'\r",NameOfWave(w3d), hName
	if ((numtype(xlo)+numtype(xhi) + numtype(ylo)+numtype(yhi) + numtype(zlo)+numtype(zhi))<6)
		printf "   Only using points in volume (%g, %g)  (%g, %g)  (%g, %g)%s\r",xlo,xhi, ylo,yhi, zlo,zhi,WaveUnits(w3d,0)
	endif
	DisplayHistogramFrom3d(Hout)
	return 0
End
//
Static Function/WAVE DisplayHistogramFrom3d(hist)
	Wave hist
	if (!WaveExists(hist))
		String hList=WaveListClass("HistogramFrom3d","*","MINROWS:2,MAXCOLS:0"), hName=""
		Prompt hName,"Histogram",popup,hList
		DoPrompt "Histogram to plot",hName
		if (V_flag)
			return $""
		endif
		Wave hist=$hName
	endif
	if (!WaveExists(hist))
		return $""
	endif
	hName = NameOfWave(hist)

	String gName = StringFromList(0,FindGraphsWithWave(hist))
	if (strlen(gName)>0)						// hist alreay plotted, just bring the plot to the front
		DoWindow/F $gName
		return hist
	endif

	// Plot hist
	Display /W=(227,45,788,392) hist
	ModifyGraph tick=2, mirror=1, minor=1, lowTrip=0.001, standoff=0
	if (DimOffset(hist,0)<0)
		ModifyGraph zero(bottom)=2
	endif
	String units = WaveUnits(hist,0)
	if (strlen(units)>0)
		if (stringmatch(hname,"R3d*"))
			Label bottom "R\\B"+hname[3]+"\\M\\U"
		elseif (stringmatch(hname,"GND*"))
			Label bottom "GND  (\\u)"
		elseif (stringmatch(hname,"totalPeakIntensity*"))
			Label bottom "Total Peak Intensity  (\\u)"
		elseif (stringmatch(hname,"TotalIntensity*"))
			Label bottom "Total Intensity  (\\u)"
		else
			Label bottom "(\\U)"
		endif
	endif
	String normType=StringByKey("normType",note(hist),"=")
	if (strlen(normType)>0 && !stringmatch(normType,"*none*"))
		Label left "Histogram/"+normType+" \\u"
	endif
	SetAxis/A/E=1 left

	String hNote=note(hist)
	String title = StringByKey("title",hNote,"=")
	String sampleName = StringByKey("sampleName",hNote,"=")
	String userName = StringByKey("userName",hNote,"=")
	String dateStr = StringByKey("date",hNote,"=")
	Variable scanNum = NumberByKey("scanNum",hNote,"=")
	String text=title
	if (strlen(sampleName)>0)
		text += SelectString(strlen(text),"","\r")+sampleName
	endif
	if (strlen(userName)>0)
		text += SelectString(strlen(text),"","\r")+userName
	endif
	if (numtype(scanNum)==0)
		text += SelectString(strlen(text),"","\r")+"\\Z06scan #"+num2istr(scanNum)+"\\Z]0"
	endif
	if (strlen(dateStr)>0)
		text += SelectString(strlen(text),"","\r")+"\\Z06"+dateStr+"\\Z]0"
	endif
	if (strlen(text))
		TextBox/C/N=text0/F=0/B=1/X=2/Y=3 text
	endif
	return hist
End

// **********************************  End of Histograms ***********************************
// *****************************************************************************************



// *****************************************************************************************
// ******************************  Start of Deviatoric Strain ******************************

Function DeviatoricStrainRefineXML_ALL(range,constrain,[coords,pattern])
	String range					// numeric range
	String constrain				// constraint on optimization, "111111", a 1 is refine, a 0 is keep constant
	Variable coords					// coordinate system to pass in return {Crystal=0, 1=BL, 2=XHF, 3=Sample (outward surface normal)}
	Variable pattern				// only zero make sense (the default)
	pattern = ParamIsDefault(pattern) ? 0 : pattern

	Wave/T imageNames=imageNames
	if (!WaveExists(imageNames))
		return 1
	endif
	Variable N = DimSize(imageNames,0)

	String xmlFileFull = StringBykey("xmlFileFull",note(imageNames),"=")
	GetFileFolderInfo/P=Igor/Q/Z=1 xmlFileFull
	if (!V_isFile)
		Variable f
		Open/D=2/R/F=XMLfilters/M="XML file: "+xmlFileFull/Z=2 f as xmlFileFull
		if (V_flag)
			return 1
		endif
		xmlFileFull = S_fileName
	endif

	if (ItemsInRange(range)<1 || ParamIsDefault(coords))
		coords = ParamIsDefault(coords) ? 1 : coords
		if (ItemsInRange(range)<1)
			range = "0-"+num2istr(N-1)
		endif
		Prompt range,"range of number, must not exceed [0,"+num2str(N-1)+"]"
		Prompt coords,"coordinate system for strain tensor", popup, "Crystal;Beam LIne;XHF;Sample (outward normal)"
		coords += 1
		DoPrompt "Strain Refine",range,coords
		if (V_flag)
			return 1
		endif
		coords -= 1
	endif

	Variable aFit,bFit,cFit,alphaFit,betFit,gamFit=NaN
	sscanf constrain, "%1d%1d%1d%1d%1d%1d",aFit,bFit,cFit,alphaFit,betFit,gamFit
	if (V_flag!=6 || !((aFit+bFit+cFit+alphaFit+betFit+gamFit)>0) || str2num(constrain[0,2])>=111)
		aFit = NaN												// flags bad input
	endif
	if (numtype(aFit))								// need to ask user
		if (indexing#StartStrainPanel(pattern,1,aFit,bFit,cFit,alphaFit,betFit,gamFit)<0)
			return 1
		endif
	endif
	aFit = aFit ? 1 : 0 ;				bFit = bFit ? 1 : 0 ;			cFit = cFit ? 1 : 0		// ensure only 0 or 1
	alphaFit = alphaFit ? 1 : 0 ;		betFit = betFit ? 1 : 0 ;		gamFit = gamFit ? 1 : 0
	Variable Nfit = aFit+bFit+cFit+alphaFit+betFit+gamFit
	if (Nfit<1 || (aFit&&bFit&&cFit))
		return 1											// this is deviatoric, you cannot fit all 3 lengths
	endif
	sprintf constrain "%d%d%d%d%d%d", aFit,bFit,cFit,alphaFit,betFit,gamFit		// constrain is for {a,b,c,alpha,beta,gamma}

	printf "DeviatoricStrainRefineXML_ALL(\"%s\",\"%s\",coords=%g)\r",range,constrain,coords
	if (!(str2num(range)>=0 && lastInRange(range)<N))
		DoAlert 0,"Illeagal numeric range "+range
		return 1
	endif

	N = ItemsInRange(range)												// actual number to do
	Make/N=(3,3,N)/O epsilonN=NaN
	Make/N=(N)/O vonMisesN=NaN
	Note/K epsilonN, "waveClass=strain_tensor"
	Note/K vonMisesN, "waveClass=vonMisesStrain;Random3dArrays"

	Variable m, once=1,j, ok
	String progressWin = ProgressPanelStart("",stop=1,showTime=1)		// display a progress bar (move to here, constrain set by first call) 
	for (m=str2num(range),j=0,ok=0; numtype(m)==0; m=NextInRange(range,m),j+=1)
		if (mod(j,10)==0)
			if (ProgressPanelUpdate(progressWin,j/N*100))				// update progress bar
				break														//   and break out of loop
			endif
		endif
		Wave eWave = $DeviatoricStrainRefineXML(m,pattern,constrain,coords=coords,xmlFileFull=xmlFileFull,printit=0)
		if (WaveExists(eWave))
			ok += 1
			epsilonN[][][j]  = eWave[p][q]
			vonMisesN[j] = NumberByKey("vonMisesStrain", note(eWave),"=")
			if (once)
				Wave PeaksForStrain=PeaksForStrain						// save first good one
				constrain = StringBykey("constrain",note(PeaksForStrain),"=")
				once = 0
			endif
		endif 
	endfor
	Variable sec = SecondsInProgressPanel(progressWin)
	printf "Successfully extracted %d Deviatoric strains from %d points,   total execution time = %s\r",ok,j,Secs2Time(sec,5,0)
	DoWindow/K $progressWin
	if (sec>10*60)																// if execution took more than 10 min, automatically save
		print "This took more than 10min, so save the experiment"
		SaveExperiment
	endif

	String noteStr="waveClass=strain_tensor"									// re-do wave notes with more information
	noteStr = ReplaceStringByKey("constrain",noteStr,constrain,"=")
	noteStr = ReplaceNumberByKey("patternNum",noteStr,pattern,"=")
	Note/K epsilonN, noteStr
	Note/K vonMisesN, ReplaceStringByKey("waveClass",noteStr,"vonMisesStrain;Random3dArrays","=")

	KillWaves/Z epsilon,epsilonBL,epsilonXHF,epsilonSample
	KillWaves/Z FullPeakList,FullPeakIndexed,PeaksForStrain
	return 0
End

Function/T DeviatoricStrainRefineXML(m,pattern,constrain,[coords,xmlFileFull,printIt])
	Variable m						// point number in list from xml file
	Variable pattern				// pattern number, usually 0
	String constrain				// constraint on optimization, "111111", a 1 is refine, a 0 is keep constant
	Variable coords					// coordinate system to pass in return {Crystal, 1=BL, 2=XHF, 3=Sample (outward surface normal)}
	String xmlFileFull				// full path name to xml file
	Variable printIt					// force full print out of results
	xmlFileFull = SelectString(ParamIsDefault(xmlFileFull),xmlFileFull,"")
	pattern = 0									// do not yet handle other patterns

	Wave/T imageNames=imageNames
	if (!WaveExists(imageNames))
		return ""
	endif
	if (strlen(xmlFileFull)<1)
		xmlFileFull = StringBykey("xmlFileFull",note(imageNames),"=")
	endif
	GetFileFolderInfo/P=Igor/Q/Z=1 xmlFileFull
	if (!V_isFile)
		Variable f
		Open/D=2/R/F=XMLfilters/M="file: "+xmlFileFull/Z=2 f as xmlFileFull
		if (V_flag)
			return ""
		endif
		xmlFileFull = S_fileName
	endif
	String names = MakeStepIndexForXML(xmlFileFull)
	Wave indexPos = $StringFromLIst(0,names)
	Wave/T indexImage = $StringFromLIst(1,names)
	if (!WaveExists(indexPos) || !WaveExists(indexImage))
		return ""
	endif

	String imageName = imageNames[m]
	imageName = ParseFilePath(3, imageName, ":", 0, 0)	// if it is Mac
	imageName = ParseFilePath(3,  imageName, "/", 0, 0)	// or Unix
	imageName = ParseFilePath(3,  imageName, "\\", 0, 0)	// or Windows
	Variable N = DimSize(imageNames,0)

	if (numtype(m) || m<0 || m>=N)
		m = 1
		Prompt m,"point index from xml file [0,"+num2str(N)+"]"
		DoPrompt "index",m
		if (V_flag)
			return ""
		endif
	endif
	if (numtype(m) || m<0 || m>=N)
		return ""
	endif

	// get step containing imageName from file
	Variable Nindex=DimSize(indexPos,0), index
	for (index=0;index<Nindex;index+=1)
		if (stringmatch(indexImage[index],imageName))
			break
		endif
	endfor
	if (index>=Nindex)
		return ""
	endif
	Variable i0=indexPos[index], i1=indexPos[index+1]
	if (numtype(i0+i1) || i1<=i0)
		return ""
	endif

	Open/R/Z=1 f as xmlFileFull						// start reading here
	FStatus f
	Variable fileLen = V_logEOF							// length of file in bytes
	if (i1>=V_logEOF)
		Close f
		return ""
	endif
	FSetPos f, i0
	String step = PadString("",i1-i0+1,0x20)			// a buffer to hold step (of length i1-i0+1 bytes)
	FBinRead f, step										// initial read
	Close f

	String svec = ReplaceString(" ",xmlTagContents("Xpixel",step),";")
	Variable Nlen = ItemsInList(svec)
	if (Nlen<4)
		return ""
	endif
	Make/N=(Nlen,11)/O FullPeakList
	FullPeakList = NaN
	FullPeakList[][0] = str2num(StringFromList(p,svec))
	svec = ReplaceString(" ",xmlTagContents("Ypixel",step),";")
	FullPeakList[][1] = str2num(StringFromList(p,svec))
	svec = ReplaceString(" ",xmlTagContents("Intens",step),";")
	FullPeakList[][10] = str2num(StringFromList(p,svec))

	String roi = xmlTagKeyVals("ROI",step)
	roi = RemoveFromList("/",roi)
	String wnote = ReplaceStringByKey("waveClass","","FittedPeakList","=")
	wnote += roi
	wnote = ReplaceStringByKey("detectorID",wnote,xmlTagContents("detectorID",step),"=")
	SetDimLabel 1,0,x0,FullPeakList		;	SetDimLabel 1,1,y0,FullPeakList
	SetDimLabel 1,2,x0Err,FullPeakList	;	SetDimLabel 1,3,y0Err,FullPeakList
	SetDimLabel 1,4,fwx,FullPeakList		;	SetDimLabel 1,5,fwy,FullPeakList
	SetDimLabel 1,6,fwxErr,FullPeakList	;	SetDimLabel 1,7,fwyErr,FullPeakList
	SetDimLabel 1,8,correlation,FullPeakList ;	SetDimLabel 1,9,correlationErr,FullPeakList
	SetDimLabel 1,10,area,FullPeakList
	Note/K FullPeakList,wnote

	svec = ReplaceString(" ",xmlTagContents("h",step),";")
	Variable Ni = ItemsInList(svec)
	if (Ni<4)
		return ""
	endif
	Make/N=(Ni,12,1)/O FullPeakIndexed=NaN
	FullPeakIndexed[][3] = str2num(StringFromList(p,svec))
	svec = ReplaceString(" ",xmlTagContents("k",step),";")
	FullPeakIndexed[][4] = str2num(StringFromList(p,svec))
	svec = ReplaceString(" ",xmlTagContents("l",step),";")
	FullPeakIndexed[][5] = str2num(StringFromList(p,svec))

	Variable SpaceGroup = str2num(ReplaceString(" ",xmlTagContents("SpaceGroup",step),";"))
	if (!(SpaceGroup>=1 && SpaceGroup<=230))
		return ""
	endif

	wnote = ReplaceStringByKey("waveClass","","IndexedPeakList","=")
	wnote += roi
	wnote = ReplaceNumberByKey("SpaceGroup",wnote,SpaceGroup,"=")
	String str = "{{"
	str += ReplaceString(" ",xmlTagContents("astar",step),",") + "}{"
	str += ReplaceString(" ",xmlTagContents("bstar",step),",") + "}{"
	str += ReplaceString(" ",xmlTagContents("cstar",step),",") + "}}"
	wnote = ReplaceStringByKey("recip_lattice0",wnote,str,"=")
	Note/K FullPeakIndexed, wnote
	SetDimLabel 1,0,Qx,FullPeakIndexed		;	SetDimLabel 1,1,Qy,FullPeakIndexed
	SetDimLabel 1,2,Qz,FullPeakIndexed		;	SetDimLabel 1,3,h,FullPeakIndexed
	SetDimLabel 1,4,k,FullPeakIndexed			;	SetDimLabel 1,5,l,FullPeakIndexed
	SetDimLabel 1,6,Intensity,FullPeakIndexed	;	SetDimLabel 1,7,keV,FullPeakIndexed
	SetDimLabel 1,8,angleErr,FullPeakIndexed	;	SetDimLabel 1,9,pixelX,FullPeakIndexed
	SetDimLabel 1,10,pixelY,FullPeakIndexed	;	SetDimLabel 1,11,detNum,FullPeakIndexed

	if (ParamIsDefault(coords))
		return DeviatoricStrainRefine(pattern,constrain,printIt=0)
	else
		return DeviatoricStrainRefine(pattern,constrain,coords=coords,printIt=0)
	endif
End
//
//Function/T DeviatoricStrainRefineXML(m,pattern,constrain,[coords])
//	Variable m						// point number in list from xml file
//	Variable pattern											// pattern number, usually 0
//	String constrain											// constraint on optimization, "111111", a 1 is refine, a 0 is keep constant
//	Variable coords												// coordinate system to pass in return {1=BL, 2=XHF, 3=Sample (outward surface normal)}
//
//	pattern = 0									// do not yet handle other patterns
//
//	Wave/T imageNames=imageNames
//	if (!WaveExists(imageNames))
//		return ""
//	endif
//	String xmlFileFull = StringBykey("xmlFileFull",note(imageNames),"=")
//	Variable f
//	Open/R/Z=1 f as xmlFileFull
//	if (V_flag)
//		return ""
//	endif
//	Close f
//	String imageName = imageNames[m]
//	imageName = ParseFilePath(3, imageName, ":", 0, 0)	// if it is Mac
//	imageName = ParseFilePath(3,  imageName, "/", 0, 0)	// or Unix
//	imageName = ParseFilePath(3,  imageName, "\\", 0, 0)	// or Windows
//	Variable N = DimSize(imageNames,0)
//
//	if (numtype(m) || m<0 || m>=N)
//		m = 1
//		Prompt m,"point index from xml file [0,"+num2str(N)+"]"
//		DoPrompt "index",m
//		if (V_flag)
//			return ""
//		endif
//	endif
//	if (numtype(m) || m<0 || m>=N)
//		return ""
//	endif
//
//	String inputImage
//	String svec
//	String step
//	Variable i0,i1										// mark start and end of current   <step ></step>
//	Variable i0start=0									// where to start searching for "<step "
//	Variable bytesRead									// keeps track of how many bytes read in (should not exceed fileLen)
//	Variable size=50*1024								// size of a typical read
//
//	// start reading here
//	Open/R/Z=1 f as xmlFileFull
//	FStatus f
//	Variable fileLen = V_logEOF							// length of file in bytes
//	size = min(size,fileLen)
//	String buf = PadString("",size,0x20)				// a buffer for working space (of length size bytes)
//	String bufRead										// a buffer for reading more from the file (read up to this much each time)
//	FBinRead f, buf										// initial read
//	bytesRead = strlen(buf)
//	Variable Nlen=0
//	do
//		i0 = strsearch(buf,"<step ",i0start,2)
//		if (i0<0)
//			size = min(size,fileLen-bytesRead)			// number of bytes to read
//			bufRead = PadString("",size,0x20)
//			FBinRead f, bufRead							// could not find start of step, read more
//			bytesRead += strlen(bufRead)
//			buf = buf[i0,Inf] + bufRead
//			i0 = strsearch(buf,"<step ",i0start,2)		// search again for start of step tag
//			if (i0<0)
//				break									// give up, all done
//			endif
//		endif
//		i1 = strsearch(buf,"</step>",i0,2)
//		if (i1<0)
//			size = min(size,fileLen-bytesRead)			// number of bytes to read
//			bufRead = PadString("",size,0x20)
//			FBinRead f, bufRead							// could not find end of step, read more
//			bytesRead += strlen(bufRead)
//			buf = buf[i0,Inf] + bufRead
//			i0 = 0
//			i1 = strsearch(buf,"</step>",i0,2)			// search again for end of step tag
//			if (i1<0)
//				break									// give up, all done
//			endif
//		endif
//		i1 += 6
//
//		// this step is bracketed by [i0,i1], process it:
//		step = buf[i0,i1]
//		i0start = i1 + 1								// where to start searching for next start tag
//		inputImage = ParseFilePath(3,  xmlTagContents("inputImage",step), "/", 0, 0)
//		if (!stringmatch(inputImage,imageName))
//			continue
//		endif
//		break
//	while(1)
//	Close f
//
//	svec = ReplaceString(" ",xmlTagContents("Xpixel",step),";")
//	Nlen = ItemsInList(svec)
//	if (Nlen<4)
//		return ""
//	endif
//	Make/N=(Nlen,11)/O FullPeakList
//	FullPeakList = NaN
//	FullPeakList[][0] = str2num(StringFromList(p,svec))
//	svec = ReplaceString(" ",xmlTagContents("Ypixel",step),";")
//	FullPeakList[][1] = str2num(StringFromList(p,svec))
//	svec = ReplaceString(" ",xmlTagContents("Intens",step),";")
//	FullPeakList[][10] = str2num(StringFromList(p,svec))
//
//	String roi = xmlTagKeyVals("ROI",step)
//	roi = RemoveFromList("/",roi)
//	String wnote = ReplaceStringByKey("waveClass","","FittedPeakList","=")
//	wnote += roi
//	wnote = ReplaceStringByKey("detectorID",wnote,xmlTagContents("detectorID",step),"=")
//	SetDimLabel 1,0,x0,FullPeakList		;	SetDimLabel 1,1,y0,FullPeakList
//	SetDimLabel 1,2,x0Err,FullPeakList	;	SetDimLabel 1,3,y0Err,FullPeakList
//	SetDimLabel 1,4,fwx,FullPeakList		;	SetDimLabel 1,5,fwy,FullPeakList
//	SetDimLabel 1,6,fwxErr,FullPeakList	;	SetDimLabel 1,7,fwyErr,FullPeakList
//	SetDimLabel 1,8,correlation,FullPeakList ;	SetDimLabel 1,9,correlationErr,FullPeakList
//	SetDimLabel 1,10,area,FullPeakList
//	Note/K FullPeakList,wnote
//
//	svec = ReplaceString(" ",xmlTagContents("h",step),";")
//	Variable Ni = ItemsInList(svec)
//	if (Ni<4)
//		return ""
//	endif
//	Make/N=(Ni,12,1)/O FullPeakIndexed=NaN
//	FullPeakIndexed[][3] = str2num(StringFromList(p,svec))
//	svec = ReplaceString(" ",xmlTagContents("k",step),";")
//	FullPeakIndexed[][4] = str2num(StringFromList(p,svec))
//	svec = ReplaceString(" ",xmlTagContents("l",step),";")
//	FullPeakIndexed[][5] = str2num(StringFromList(p,svec))
//
//	Variable SpaceGroup = str2num(ReplaceString(" ",xmlTagContents("SpaceGroup",step),";"))
//	if (!(SpaceGroup>=1 && SpaceGroup<=230))
//		return ""
//	endif
//
//	wnote = ReplaceStringByKey("waveClass","","IndexedPeakList","=")
//	wnote += roi
//	wnote = ReplaceNumberByKey("SpaceGroup",wnote,SpaceGroup,"=")
//	String str = "{{"
//	str += ReplaceString(" ",xmlTagContents("astar",step),",") + "}{"
//	str += ReplaceString(" ",xmlTagContents("bstar",step),",") + "}{"
//	str += ReplaceString(" ",xmlTagContents("cstar",step),",") + "}}"
//	wnote = ReplaceStringByKey("recip_lattice0",wnote,str,"=")
//	Note/K FullPeakIndexed, wnote
//	SetDimLabel 1,0,Qx,FullPeakIndexed		;	SetDimLabel 1,1,Qy,FullPeakIndexed
//	SetDimLabel 1,2,Qz,FullPeakIndexed		;	SetDimLabel 1,3,h,FullPeakIndexed
//	SetDimLabel 1,4,k,FullPeakIndexed			;	SetDimLabel 1,5,l,FullPeakIndexed
//	SetDimLabel 1,6,Intensity,FullPeakIndexed	;	SetDimLabel 1,7,keV,FullPeakIndexed
//	SetDimLabel 1,8,angleErr,FullPeakIndexed	;	SetDimLabel 1,9,pixelX,FullPeakIndexed
//	SetDimLabel 1,10,pixelY,FullPeakIndexed	;	SetDimLabel 1,11,detNum,FullPeakIndexed
//
//	if (ParamIsDefault(coords))
//		return DeviatoricStrainRefine(pattern,constrain)
//	else
//		return DeviatoricStrainRefine(pattern,constrain,coords=coords)
//	endif
//End

Static Function/T MakeStepIndexForXML(xmlFileFull)
	String xmlFileFull
	Variable f
	Open/R/Z=1 f as xmlFileFull							// cannot open file
	if (V_flag)
		return ""
	endif
	FStatus f
	Close f
	Variable fileLen = V_logEOF							// length of file in bytes
	String hashStr = Hash(S_path+S_fileName+num2istr(V_logEOF),1)

	String namePos = CleanupName("IndexFilePos_"+ParseFilePath(3,xmlFileFull,":",0,0),0)
	String nameImage = CleanupName("IndexFileImage_"+ParseFilePath(3,xmlFileFull,":",0,0),0)
	Wave indexPos = $namePos
	Wave/T indexImage = $nameImage
	Variable redoIndex=!WaveExists(indexPos) || !WaveExists(indexImage)
	if (!redoIndex)
		redoIndex = !stringmatch(hashStr,StringByKey("hash", note(indexPos),"="))
	endif
	if (!redoIndex)										// no need to re-calculate index fwaves, just return them
		return GetWavesDataFolder(indexPos,2)+";"+GetWavesDataFolder(indexImage,2)
	endif

	// re-do the index waves from the xml file
	Variable Nsize=1000
	Make/N=(Nsize)/D/O  $namePos=NaN
	Wave indexPos = $namePos
	//	print namePos
	//	print nameImage
	Make/N=(Nsize)/T/O  $nameImage=""
	Wave/T indexImage = $nameImage

	String noteStr=ReplaceStringByKey("waveClass","","XMLindexPositions","=")
	noteStr=ReplaceStringByKey("path",noteStr,S_path,"=")
	noteStr=ReplaceStringByKey("fileName",noteStr,S_fileName,"=")
	noteStr=ReplaceNumberByKey("fileLength",noteStr,fileLen,"=")
	noteStr=ReplaceStringByKey("hash",noteStr,hashStr,"=")
	Note/K indexPos,noteStr
	Note/K indexImage, ReplaceStringByKey("waveClass",noteStr,"XMLindexImageNames","=")

	String inputImage
	String step
	Variable i0,i1										// mark start and end of current   <step ></step>
	Variable i0start=0									// where to start searching for "<step "
	Variable bytesRead									// keeps track of how many bytes read in (should not exceed fileLen)
//	Variable size=50*1024								// size of a typical read
	Variable size=200*1024							// size of a typical read
	Variable PosOfBuf=0								// position of start of buffer
	Variable N											// <step> index (and afterwards the total number of points)

	// start reading here
	Open/R/Z=1 f as xmlFileFull
	size = min(size,fileLen)
	String buf = PadString("",size,0x20)				// a buffer for working space (of length size bytes)
	String bufRead										// a buffer for reading more from the file (read up to this much each time)
	FBinRead f, buf										// initial read
	bytesRead = strlen(buf)
	N = 0
	do
		i0 = strsearch(buf,"<step ",i0start,2)
		if (i0<0)										// read more (extend buf)
			size = min(size,fileLen-bytesRead)			// number of bytes to read
			bufRead = PadString("",size,0x20)
			FBinRead f, bufRead							// could not find start of step, read more
			bytesRead += strlen(bufRead)
			//	buf = buf[i0,Inf] + bufRead
			buf = buf[0,Inf] + bufRead
			i0 = strsearch(buf,"<step ",i0start,2)		// search again for start of step tag
			if (i0<0)
				break									// give up, all done
			endif
		endif

		if (N>=Nsize)									// extend arrays
			Nsize += 200
			Redimension/N=(Nsize) indexPos, indexImage
		endif

		indexPos[N] = PosOfBuf+i0						// position of start of this <step>
		N += 1											// increment here, because some <step>s may not have an image
		i1 = strsearch(buf,"</step>",i0,2)				// position of end of step
		if (i1<0)
			size = min(size,fileLen-bytesRead)			// number of bytes to read
			bufRead = PadString("",size,0x20)
			FBinRead f, bufRead							// could not find end of step, read more
			bytesRead += strlen(bufRead)
			buf = buf[i0,Inf] + bufRead
			PosOfBuf += i0								// start of buf shifted higher by i0
			i0 = 0
			i1 = strsearch(buf,"</step>",i0,2)			// search again for end of step tag
			if (i1<0)
				break									// give up, all done
			endif
		endif
		i1 += 6

		// this step is bracketed by [i0,i1], process it:
		step = buf[i0,i1]
		i0start = i1 + 1								// where to start searching for next start tag
		inputImage = ParseFilePath(3,  xmlTagContents("inputImage",step), "/", 0, 0)
		if (strlen(inputImage)<1)
			continue
		endif
		indexImage[N-1] = inputImage					// name of image for this <step>
	while(1)
	Close f

	Nsize = N
	Redimension/N=(Nsize) indexPos, indexImage

	return GetWavesDataFolder(indexPos,2)+";"+GetWavesDataFolder(indexImage,2)
End
//
//	Function test_MakeStepIndexForXML()
//		Variable f
//		// String xmlFileFull = "Macintosh HD:Users:tischler:Desktop:SampleN_2D2_resulttest.xml"
//		String xmlFileFull = "Macintosh HD:Users:tischler:Desktop:SampleN_2D2_result.xml"
//	
//		Open/R/Z=1 f as xmlFileFull
//		if (V_flag)
//			Open/R/T=".xml" f
//			Close f
//			xmlFileFull = S_filename
//			print "opening ",xmlFileFull
//		endif
//	
//		String names = MakeStepIndexForXML(xmlFileFull)
//		printf "index positions =  %s\r",StringFromLIst(0,names)
//		printf "index images =  %s\r",StringFromLIst(1,names)
//		Wave indexPos = $StringFromLIst(0,names)
//		Wave/T indexImage = $StringFromLIst(1,names)
//		if (!WaveExists(indexPos) || !WaveExists(indexImage))
//			print "at least one of the waves does not exist!"
//		endif
//	End


Function MakeMany_vonMises(epsilonN)					// create wave with von MIses strain from list of epsilons
	Wave epsilonN
	if (!WaveExists(epsilonN))
		return 1
	endif
	Variable i,N=DimSize(epsilonN,2)
	if (N<1)
		return 1
	endif
	Make/N=(N)/O vonMisesN=NaN
	Variable epsilonvM										// von Mises strain
	Make/N=(3,3)/FREE/D epsilon
	for (i=0;i<N;i+=1)
		epsilon = epsilonN[p][q][i]
		epsilonvM = (epsilon[0][0]-epsilon[1][1])^2 + (epsilon[1][1]-epsilon[2][2])^2 + (epsilon[2][2]-epsilon[0][0])^2 
		epsilonvM += 6*( epsilon[0][1]^2 + epsilon[1][2]^2 + epsilon[2][0]^2 )
		epsilonvM = sqrt(epsilonvM/2)
		vonMisesN[i] = epsilonvM
	endfor
End

// *******************************  End of Deviatoric Strain *******************************
// *****************************************************************************************



// *****************************************************************************************
// *********************************  Start of Make Gizmo **********************************

Function MakeGizmo_xmlData(scatt)
	Wave scatt

	// Do nothing if the Gizmo XOP is not available.
	if(exists("NewGizmo")!=4)
		DoAlert 0, "Gizmo XOP must be installed"
		return 1
	endif

	if (!WaveExists(scatt))
		String wStr="",wList = WaveListClass("Random3dArraysXYZ","*","DIMS:2,MAXCOLS:3,MINCOLS:3")
		if (ItemsInList(wList)<1)
			DoAlert 0,"triplet xyz wave not found"
			return 1
		elseif (ItemsInList(wList)==1)
			wStr = StringFromLIst(0,wList)
		else
			Prompt wStr, "xyz triplet wave",popup, wList
			DoPrompt "xyz Triplets",wStr
			if (V_flag)
				return 1
			endif
		endif
		Wave scatt = $wStr
		printf "MakeGizmo_xmlData(%s)\r",NameOfWave(scatt)
	endif
	if (!WaveExists(scatt))
		DoAlert 0,"triplet xyz wave not found"
		return 1
	endif
	if (WaveDims(scatt)!=2 || DimSize(scatt,1)!=3)
		DoAlert 0,"triplet wave "+NameOfWave(scatt)+" is not of dimension [N][3]"
		return 1
	endif
	Variable N=DimSize(scatt,0)

	String options
	sprintf options, "DIMS:2,MINROWS:%d,MAXROWS:%d,MAXCOLS:4,MINCOLS:4",N,N
	String rgbaList=WaveListClass("Random3dArraysRGBA","*",options)
	if (WhichListItem("rgba",rgbaList)<0 && Exists("rgba")==1)		// rgba not in rgbaList, but it exists
		rgbaList += "rgba;"
	endif
	if (ItemsInList(rgbaList)==1)
		Wave rgba=$StringFromList(0,rgbaList)
	elseif (ItemsInList(rgbaList)>1)
		String rgbaName
		Prompt rgbaName, "RGBA wave",popup,rgbaList
		DoPrompt "RGBA" rgbaName
		if (V_flag)
			return 1
		endif
		Wave rgba=$rgbaName
	else
		Wave rgba=$(GetWavesDataFolder(scatt,1)+"rgba")
	endif

	String title=getTitleFromNote(note(scatt)), title2=""
	if (WaveExists(rgba))
		title2= StringByKey("source",note(rgba),"=")
	endif
	if (strlen(title2))
		title = title2+SelectString(strlen(title),"",", ")+title
	endif
	Wave cubeCorners=$(GetWavesDataFolder(scatt,1)+CleanupName(NameOfWave(scatt)+"Corners",0))
	if (!WaveExists(cubeCorners))
		Wave cubeCorners=$(GetWavesDataFolder(scatt,1)+"cubeCorners")			// the old way for compatibility
	endif

	String name=UniqueName("Gizmo",5,0)
	Execute "NewGizmo/N="+name+"/T=\""+name+"\" /W=(253,44,804,566)"
	Execute "ModifyGizmo startRecMacro"

	if (strlen(title))
		Execute "AppendToGizmo string=\""+title+"\",strFont=\"Geneva\",name=Title"
	endif
	Execute "AppendToGizmo attribute blendFunc={770,771},name=blendFunc0"
	Execute "AppendToGizmo Axes=boxAxes,name=axes0"
	setGizmoAxisLabels("X  (µm)","H  (µm)","F  (µm)")

	Execute "AppendToGizmo Scatter="+GetWavesDataFolder(scatt,2)+",name=scatter0"
	Execute "ModifyGizmo ModifyObject=scatter0 property={ scatterColorType,1}"
	Execute "ModifyGizmo ModifyObject=scatter0 property={ markerType,0}"
	Execute "ModifyGizmo ModifyObject=scatter0 property={ sizeType,0}"
	Execute "ModifyGizmo ModifyObject=scatter0 property={ rotationType,0}"
	Execute "ModifyGizmo ModifyObject=scatter0 property={ Shape,5}"
	Execute "ModifyGizmo ModifyObject=scatter0 property={ size,0.01}"
	if (WaveExists(rgba))
		Execute "ModifyGizmo ModifyObject=scatter0 property={ colorWave,"+GetWavesDataFolder(rgba,2)+"}"
	endif
	Execute "ModifyGizmo ModifyObject=scatter0 property={ CTABScaling,96}"
	Execute "ModifyGizmo ModifyObject=scatter0 property={ MinRGBA,0,0.699992,3.0518e-05,3.0518e-05,1}"
	Execute "ModifyGizmo ModifyObject=scatter0 property={ MaxRGBA,0,0,1,0.987213,1}"

	if (WaveExists(cubeCorners))
		Execute "AppendToGizmo Scatter="+GetWavesDataFolder(cubeCorners,2)+",name=cubeCorners"
		Execute "ModifyGizmo ModifyObject=cubeCorners property={ scatterColorType,0}"
		Execute "ModifyGizmo ModifyObject=cubeCorners property={ markerType,0}"
		Execute "ModifyGizmo ModifyObject=cubeCorners property={ sizeType,0}"
		Execute "ModifyGizmo ModifyObject=cubeCorners property={ rotationType,0}"
		Execute "ModifyGizmo ModifyObject=cubeCorners property={ Shape,5}"
		Execute "ModifyGizmo ModifyObject=cubeCorners property={ size,0.01}"
		Execute "ModifyGizmo ModifyObject=cubeCorners property={ color,0,0,0,1}"
	endif

	Execute "AppendToGizmo light=Directional,name=light0"
	Execute "ModifyGizmo light=light0 property={ position,0.426121,-0.439519,0.790724,0.0}"
	Execute "ModifyGizmo light=light0 property={ direction,0.426121,-0.439519,0.790724}"
	Execute "ModifyGizmo light=light0 property={ specular,0.733333,0.733333,0.733333,0.5}"

	if (strlen(title))
		Execute "ModifyGizmo setDisplayList=0, opName=translateTitle, operation=translate, data={-1.9,1.9,0}"
		Execute "ModifyGizmo setDisplayList=1, opName=scaleTitle, operation=scale, data={0.1,0.1,0.1}"
		Execute "ModifyGizmo setDisplayList=2, opName=rotateTitle, operation=rotate, data={180,1,0,0}"
		Execute "ModifyGizmo setDisplayList=3, object=Title"
		Execute "ModifyGizmo setDisplayList=4, opName=MainTransform, operation=mainTransform"
	endif
	Execute "ModifyGizmo setDisplayList=-1, opName=ortho0, operation=ortho, data={-2,2,-2,2,-3,3}"
	Execute "ModifyGizmo setDisplayList=-1, opName=scale0, operation=scale, data={1.25,1.25,1.25}"
	Execute "ModifyGizmo setDisplayList=-1, attribute=blendFunc0"
	Execute "ModifyGizmo setDisplayList=-1, object=axes0"
	Execute "ModifyGizmo setDisplayList=-1, object=scatter0"
	if (WaveExists(cubeCorners))
		Execute "ModifyGizmo setDisplayList=-1, object=cubeCorners"
	endif
	Execute "ModifyGizmo setDisplayList=-1, opName=clearColor0, operation=clearColor, data={0.733,0.733,0.733,1}"
	Execute "ModifyGizmo SETQUATERNION={0.824516,0.417235,-0.157751,-0.348105}"
	Execute "ModifyGizmo autoscaling=1"
	Execute "ModifyGizmo currentGroupObject=\"\""
	Execute "ModifyGizmo compile"
	Execute "ModifyGizmo endRecMacro"
End
//
Static Function/T getTitleFromNote(list)
	String list
	String title = StringByKey("title",list,"=")
	String sampleName = StringByKey("sampleName",list,"=")
	String userName = StringByKey("userName",list,"=")
	String dateStr = StringByKey("date",list,"=")
	String str=""
	sprintf str,"%s, %s, %s"title,sampleName,userName
	Variable yr,month,day,h,m,s
	sscanf dateStr,"%d-%d-%d %d:%d:%d",yr,month,day,h,m,s
	if (V_flag==6)
		Variable seconds = date2secs(yr,month,day) + 3600*h+60*m+s
		str += ", "+Secs2Date(seconds,2)+"  "+Secs2Time(seconds,0)
	endif
	return str
End



Function/WAVE MakeGizmoRGBAforXYZ(wxyz,values,cTab,lo,hi)
	Wave wxyz				// list of xyz values, need this to know how many xyz triplets to search for in values
	Wave values
	String cTab				// name of color table
	Variable lo,hi			// scaling range

	Variable printIt=0
	String options="DIMS:2,MINCOLS:3,MAXCOLS:3"
	String xyzList=WaveListClass("Random3dArraysXYZ","*",options)
	if (!WaveExists(wxyz) && ItemsInLIst(xyzList)<1)			// no xyz waves
		DoAlert 0,"No xyz waves in this data folder"
		return $""
	elseif (!WaveExists(wxyz) && ItemsInLIst(xyzList)==1)	// just one xyz wave, use it
		Wave wxyz=$StringFromList(0,xyzList)
	endif

	if (!WaveExists(wxyz))
		String sxyz
		Prompt sxyz,"xyz wave", popup, xyzList
		DoPrompt "xyz wave", sxyz
		if (V_flag)
			return $""
		endif
		Wave wxyz=$sxyz
		printIt = 1
	endif
	if (!WaveExists(wxyz))
		return $""
	endif
	Variable N=DimSize(wxyz,0)
	if (N<1)
		return $""
	endif

	sprintf options,"DIMS:1,MINROWS:%d,MAXROWS:%d",N,N
	String valueList=WaveListClass("Random3dArrays","*",options)
	if (!WaveExists(values) && ItemsInLIst(valueList)<1)		// no value waves
		DoAlert 0,"No value waves in this data folder"
		return $""
	elseif (!WaveExists(values) && ItemsInLIst(valueList)==1)	// just one value wave, use it
		Wave values=$StringFromList(0,valueList)
	endif

	// add ability to do full RxRyRz too
	if (WhichListItem("RX",valueList)>=0 && WhichListItem("RH",valueList)>=0 && WhichListItem("RF",valueList)>=0)
		valueList = AddListItem("RxRyRz",valueList)
	endif

	Variable Rxyz=0,InvPoleCubic=0
	String svalues=""
	if (!WaveExists(values))
		STRUCT crystalStructure xtal
		if (!FillCrystalStructDefault(xtal))			// get xtal
			if (xtal.SpaceGroup >= 195)			//	CUBIC
				valueList = AddListItem("Inverse Pole Figure (CUBIC)",valueList,";",Inf)
			endif
		endif
		Prompt svalues,"value wave", popup, valueList
		DoPrompt "Values for xyz", svalues
		if (V_flag)
			return $""
		endif
		Rxyz = stringmatch(svalues,"RxRyRz")
		InvPoleCubic = stringmatch(svalues,"Inverse Pole Figure (CUBIC)")
		Wave values=$svalues
		printIt = 1
	endif
	if (Rxyz)
		svalues = "RxRyRz"
	elseif (InvPoleCubic)
		svalues = "InvPoleCubic"
	else
		svalues = NameOfWave(values)
	endif

	if (InvPoleCubic)						// do not need lo or hi for InvPoleCubic, but need a surface normal
		lo = 0
		hi = 0
		String normStr="0,1,-1"
		Prompt normStr, "Inverse Pole Figure Surface Normal (Beam Line Coords)"
		DoPrompt "Surface Normal", normStr
		if (V_flag)
			return $""
		endif
		Wave InvPoleNormal=str2vec(normStr)
		if (!WaveExists(InvPoleNormal))
			return $""
		elseif (DimSize(InvPoleNormal,0)!=3 || WaveDims(InvPoleNormal)!=1 || norm(InvPoleNormal)<=0)
			return $""
		endif
		normStr = hkl2str(InvPoleNormal[0],InvPoleNormal[1],InvPoleNormal[2])
		cTab = "InvPole"
	elseif (numtype(lo+hi) && Rxyz)		// prompt for RxRyRz
		hi = max(WaveMax(Rx),WaveMax(Rh))
		hi = max(hi,WaveMax(Rh))
		lo = min(WaveMin(Rx),WaveMin(Rh))
		lo = min(lo,WaveMin(Rh))
		hi = max(hi,abs(lo))
		hi = !(hi>0) ? 0.1 : 2*atan(hi)*180/PI	// Rodriques = tan(theta/2)
		Prompt hi,"maximum rotation, for saturated color (degree)"
		DoPrompt "Max Angle (degre)",hi
		if (V_flag)
			return $""
		endif
		lo = -hi							// leave hi in degrees
		cTab = "RxRyRz"
		printIt = 1
	elseif (numtype(lo+hi))			// prompt for a single value
		String colorList = CTabList(), colorDefault
		lo = roundSignificant(WaveMin(values),2)
		hi = roundSignificant(WaveMax(values),2)
		if (lo*hi < 0)
			hi = max(abs(hi),abs(lo))
			lo = -hi
			colorDefault = "RedWhiteBlue"
		else
			colorDefault = "SeaLandAndFire"
			if (hi>0)
				lo = 0
			elseif (lo<0)
				lo = -hi
				hi =0
			endif
		endif
		cTab = SelectString(WhichListItem(cTab,colorList)<0,cTab,colorDefault)
		Prompt lo,"low value for color scale"
		Prompt hi,"high value for color scale"
		Prompt cTab,"Color Table",popup,colorList
		DoPrompt "Color Scaling",lo,hi,cTab
		if (V_flag)
			return $""
		endif
		printIt = 1
	endif
	if (printIt)
		sxyz = SelectString(WaveExists(wxyz),"$\"\"",NameOfWave(wxyz))
		String sss = SelectString(WaveExists(values),"$\"\"",NameOfWave(values))
		printf "MakeGizmoRGBAforXYZ(%s, %s, \"%s\", %g, %g)",sxyz,sss,cTab,lo,hi
		if (Rxyz)
			printf "		// using {RH,RH,RF}"
		endif
		if (InvPoleCubic)
			printf "		// Inverse Pole Figure Colors (CUBIC), with surf normal {%s}",normStr
		endif
		printf "\r"
	endif

	if (!Rxyz && !InvPoleCubic)
		if (WhichListItem(cTab,CTabList())<0)
			DoAlert 0, "The requested color table '"+cTab+"' does not exist"
			return $""
		endif
		ColorTab2Wave $cTab
		Wave M_colors=M_colors
		Variable Nc=DimSize(M_colors,0)
		Make/N=(Nc,4)/FREE colors
		colors[][0,2] = M_colors/65535					// colors appropriate for a Gizmo
		colors[][3] = 1
		KillWaves/Z M_colors
	endif
	String rName=svalues+"_RGBA"
	Make/N=(N,4)/O $rName
	Wave rgba=$rName
	if (printIt)
		printf "made RGBA  '%s'  for '%s' using %s over range [%g, %g]%s\r",rName,NameOfWave(wxyz),cTab,lo,hi,SelectString(Rxyz,"","¡")
	endif

	if (InvPoleCubic)
		Wave gm = $(GetWavesDataFolder(wxyz,1)+"gm")
		if (!WaveExists(gm))
			print "ERROR -- Cannot Find gm Wave"
			return $""
		endif
		Make/N=(3,3)/D/FREE recip
		Variable i
		for (i=0;i<N;i+=1)
			recip = gm[p][q][i]
			MatrixOP/FREE/O hkl = Inv(recip) x InvPoleNormal
			Wave rgbai=CubicTriangleColors(hkl)			// returns rgb, hkl is hkl of surface normal
			rgba[i][0,2] = rgbai[q]
		endfor
		rgba = limit(rgba,0,1)
		rgba[][3] = 1
	elseif (Rxyz)
		Wave RX=RX, RH=RH, RF=RF
		Wave rotrgb=$makeRGBJZT(RX,RH,RF,hi)	// hi is in degrees
		rgba = rotrgb[p][q]/65535
		rgba = limit(rgba,0,1)
		rgba[][3] = 1
		KillWaves/Z rotrgb
	else
		Variable m, ic
		for (m=0;m<N;m+=1)
			ic = (values[m]-lo) / (hi-lo)  * (Nc-1)
			ic = limit(round(ic),0,Nc-1)
			rgba[m][] = colors[ic][q]
		endfor
	endif

	String sampleName=StringBykey("sampleName",note(wxyz),"="), userName=StringBykey("userName",note(wxyz),"=")
	String wNote="waveClass=Random3dArraysRGBA"
	wNote=ReplaceStringByKey("fldrName",wNote,GetDataFolder(1),"=")
	if (strlen(sampleName))
		wNote=ReplaceStringByKey("sampleName",wNote,sampleName,"=")
	endif
	if (strlen(userName))
		wNote=ReplaceStringByKey("userName",wNote,userName,"=")
	endif
	wNote=ReplaceStringByKey("source",wNote,svalues,"=")
	Note/K rgba, wNote
	return rgba
End
//Function/WAVE MakeGizmoRGBAforXYZ(wxyz,values,cTab,lo,hi)
//	Wave wxyz				// list of xyz values, need this to know how many xyz triplets to search for in values
//	Wave values
//	String cTab				// name of color table
//	Variable lo,hi			// scaling range
//
//	Variable printIt=0
//	String options="DIMS:2,MINCOLS:3,MAXCOLS:3"
//	String xyzList=WaveListClass("Random3dArraysXYZ","*",options)
//	if (!WaveExists(wxyz) && ItemsInLIst(xyzList)<1)			// no xyz waves
//		DoAlert 0,"No xyz waves in this data folder"
//		return $""
//	elseif (!WaveExists(wxyz) && ItemsInLIst(xyzList)==1)	// just one xyz wave, use it
//		Wave wxyz=$StringFromList(0,xyzList)
//	endif
//
//	if (!WaveExists(wxyz))
//		String sxyz
//		Prompt sxyz,"xyz wave", popup, xyzList
//		DoPrompt "xyz wave", sxyz
//		if (V_flag)
//			return $""
//		endif
//		Wave wxyz=$sxyz
//		printIt = 1
//	endif
//	if (!WaveExists(wxyz))
//		return $""
//	endif
//	Variable N=DimSize(wxyz,0)
//	if (N<1)
//		return $""
//	endif
//
//	sprintf options,"DIMS:1,MINROWS:%d,MAXROWS:%d",N,N
//	String valueList=WaveListClass("Random3dArrays","*",options)
//	if (!WaveExists(values) && ItemsInLIst(valueList)<1)		// no value waves
//		DoAlert 0,"No value waves in this data folder"
//		return $""
//	elseif (!WaveExists(values) && ItemsInLIst(valueList)==1)	// just one value wave, use it
//		Wave values=$StringFromList(0,valueList)
//	endif
//	if (!WaveExists(values))
//		String svalues
//		Prompt svalues,"value wave", popup, valueList
//		DoPrompt "Values for xyz", svalues
//		if (V_flag)
//			return $""
//		endif
//		Wave values=$svalues
//		printIt = 1
//	endif
//
//	String colorList = CTabList(), colorDefault
//	if (numtype(lo+hi))
//		lo = roundSignificant(WaveMin(values),2)
//		hi = roundSignificant(WaveMax(values),2)
//		if (lo*hi < 0)
//			hi = max(abs(hi),abs(lo))
//			lo = -hi
//			colorDefault = "RedWhiteBlue"
//		else
//			colorDefault = "SeaLandAndFire"
//			if (hi>0)
//				lo = 0
//			elseif (lo<0)
//				lo = -hi
//				hi =0
//			endif
//		endif
//		cTab = SelectString(WhichListItem(cTab,colorList)<0,cTab,colorDefault)
//		Prompt lo,"low value for color scale"
//		Prompt hi,"high value for color scale"
//		Prompt cTab,"Color Table",popup,colorList
//		DoPrompt "Color Scaling",lo,hi,cTab
//		if (V_flag)
//			return $""
//		endif
//		printIt = 1
//	endif
//	if (printIt)
//		printf "MakeGizmoRGBAforXYZ(%s, %s, \"%s\", %g, %g)\r",NameOfWave(wxyz),NameOfWave(values),cTab,lo,hi
//	endif
//
//	ColorTab2Wave $cTab
//	Wave M_colors=M_colors
//	Variable Nc=DimSize(M_colors,0)
//	Make/N=(Nc,4)/FREE colors
//	colors[][0,2] = M_colors/65535						// colors appropriate for a Gizmo
//	colors[][3] = 1
//	KillWaves/Z M_colors
//
//	String rName=NameOfWave(values)+"_RGBA"
//	Make/N=(N,4)/O $rName
//	Wave rgba=$rName
//	if (printIt)
//		printf "made RGBA  '%s'  for '%s' using %s over range [%g, %g]\r",rName,NameOfWave(wxyz),cTab,lo,hi
//	endif
//
//	Variable m, ic
//	for (m=0;m<N;m+=1)
//		ic = (values[m]-lo) / (hi-lo)  * (Nc-1)
//		ic = limit(round(ic),0,Nc-1)
//		rgba[m][] = colors[ic][q]
//	endfor
//
//	String sampleName=StringBykey("sampleName",note(wxyz),"="), userName=StringBykey("userName",note(wxyz),"=")
//	String wNote="waveClass=Random3dArraysRGBA"
//	wNote=ReplaceStringByKey("fldrName",wNote,GetDataFolder(1),"=")
//	if (strlen(sampleName))
//		wNote=ReplaceStringByKey("sampleName",wNote,sampleName,"=")
//	endif
//	if (strlen(userName))
//		wNote=ReplaceStringByKey("userName",wNote,userName,"=")
//	endif
//	wNote=ReplaceStringByKey("source",wNote,NameOfWave(values),"=")
//	Note/K rgba, wNote
//	return rgba
//End

// **********************************  End of Make Gizmo ***********************************
// *****************************************************************************************



// *****************************************************************************************
// *****************************  Start of Plot Random Points ******************************

Function/T Make2Dplot_xmlData(fldr,[minRange,printIt,ForceNew])
	String fldr					// folder with data to display
	Variable minRange			// a range must be this big to be real
	Variable printIt
	Variable ForceNew			// forces a new plot to be made
	minRange = ParamIsDefault(minRange) ? 1 : minRange
	minRange = numtype(minRange) || minRange<=0 ? 1 : minRange
	printIt = ParamIsDefault(printIt) ? NaN : printIt
	printIt = numtype(printIt) ? strlen(GetRTStackInfo(2))==0 : !(!printIt)
	ForceNew = ParamIsDefault(ForceNew) ? 0 : ForceNew
	ForceNew = numtype(ForceNew) ? 0 : !(!ForceNew)

	String fldrSav=GetDataFolder(1)
	if (strsearch(fldr,"root:",0)!=0 && char2num(fldr[0])!=char2num(":"))
		fldr = ":" + fldr
	endif
	Variable len=strlen(fldr)
	if (strsearch(fldr,"::",Inf,1)==(len-2) && len>1)			// if it ends in "::" then add a ":"
		fldr += ":"
	elseif (!StringMatch(fldr[len-1],":") && len>1)	// if it does not end in a ":" add one
		fldr += ":"
	endif	
	if (!DataFolderExists(fldr))
		return ""
	endif
	Wave XX=$(fldr+"XX"), YY=$(fldr+"YY"), ZZ=$(fldr+"ZZ")
	Wave HH=$(fldr+"HH"), FF=$(fldr+"FF"), depth=$(fldr+"depth")
	if (!WaveExists(XX) || !WaveExists(YY) || !WaveExists(ZZ) || !WaveExists(HH) || !WaveExists(FF))
		DoAlert 0,"Could not find one of XX,YY,ZZ,HH,FF  in \""+fldr+"\""
		return ""
	endif
	if (WaveDims(XX)!=1 || DimSize(XX,0)<2)
		DoAlert 0,"XX wave is not 1D or is too short"
		return ""
	endif
	Variable N=DimSize(XX,0)

	// determine what to plot
	Variable dX=WaveMax(XX)-WaveMin(XX), dY=WaveMax(YY)-WaveMin(YY), dZ=WaveMax(ZZ)-WaveMin(ZZ)
	Variable dH=WaveMax(HH)-WaveMin(HH), dF=WaveMax(FF)-WaveMin(FF), ddepth=0
	if (WaveExists(depth))
		ddepth=WaveMax(depth)-WaveMin(depth)
	endif
	dX = dX<minRange ? 0 : dX					// set to zero if too small to be real
	dY = dY<minRange ? 0 : dY
	dZ = dZ<minRange ? 0 : dZ
	dH = dH<minRange ? 0 : dH
	dF = dF<minRange ? 0 : dF
	ddepth = ddepth<1 ? 0 : ddepth

	Variable yRev=0, xRev=0					// flags whether to revers axis
	Wave Xw=$"", Yw=$""							// waves to use for X & Y axes
	if (dX==0 && ddepth>0 && dH>0 && dF>0)	// a wire scan with X-constant, and sample scanned in H
		Wave Xw=HH, Yw=FF
		yRev = 1
	elseif (dX==0 && ddepth>0)				// wire scan at constant X, ???
		Wave Xw=ZZ, Yw=YY
	elseif (dY==0 && ddepth>0)				// slice into the sample with moving X, Sample-Y constant
		Wave Xw=XX, Yw=depth
	elseif (ddepth==0)
		Wave Xw=XX, Yw=HH							// no wire scan, just scanned sample surface
	endif
	if (!WaveExists(Xw) || !WaveExists(Yw))
		String str="Could not figure out suitable X and Y axes"
		printf "ranges in sample are: ÆX=%g, ÆY=%g, ÆZ=%g, ÆH=%g, ÆF=%g, Ædepth=%g\r",dX,dY,dZ, dH,dF, ddepth
		print str
		DoAlert 0, str
		return ""
	endif

	String win=StringFromList(0,FindGraphsWithWave(Yw))
	if (!ForceNew && strlen(win))			// plot exists, bring it to front
		DoWindow/F $win
		return win
	endif

	String listInFldr=WaveListClass("Random3dArrays","*","DIMS:1,TEXT:0")
	listInFldr = RemoveFromList("XX;YY;ZZ;HH;FF;depth",listInFldr)

	String listUni0="Nindexed;totalSum;sumAboveThreshold;numAboveThreshold;rmsIndexed;totalAngles;"
	String listBi0="RX;RH;RF;RY;RZ;"
	listInFldr += WaveListClass("Random3dArraysRGB","*","DIMS:2,MINCOLS:3,MAXCOLS:3")
	String Xname=NameOfWave(Xw), Yname=NameOfWave(Yw), Zname
	Prompt Zname,"Wave for Z-axis Color", popup,listInFldr
	sprintf str, "Color for '%s' vs '%s' Graph",Yname,Xname
	DoPrompt str, Zname
	if (V_flag)
		return ""
	endif
	Wave Zw=$Zname
	if (!WaveExists(Zw) || !WaveExists(Yw))
		str="Could not figure out suitable Z-axis"
		print str
		DoAlert 0, str
		return ""
	endif

	Variable isRGB=DimSize(Zw,1)==3, BiPolar=0, UniPolar=0
	if (!isRGB)
		//	BiPolar = (WaveMin(Zw)*WaveMax(Zw)) < 0
		BiPolar = WaveMin(Zw)<0
		UniPolar = !BiPolar
	endif
	if (printIt)
		str = SelectString(BiPolar-isRGB,"3-Column RGB Wave","Uni-Polar","Bi-Polar")
		printf "Graphing '%s' vs '%s',  using Z-color from '%s' (%s)\r",Yname,Xname,Zname,str
	endif

	Variable Zmin=0, Zmax=WaveMax(Zw)
	if (BiPolar)
		Zmax = max(Zmax,-WaveMin(Zw))
		Zmin = -Zmax
	endif

	Display /W=(205,61,754,395) Yw vs Xw
	ModifyGraph mode=3, marker=16, tick=2, mirror=1, minor=1, lowTrip=0.001
	if (isRGB)
		ModifyGraph zColor($Yname)={Zw,0,0,directRGB}
	else
		String colorWave=SelectString(UniPolar,"RedWhiteBlue256","Terrain256")
		ModifyGraph zColor($Yname)={Zw,Zmin,Zmax,$colorWave,UniPolar}
	endif
	Label left Yname[0] + "  (\\U)"
	Label bottom Xname[0] + "  (\\U)"
	if (yRev)
		SetAxis/A/R left
	endif
	if (xRev)
		SetAxis/A/R bottom
	endif
	DoUpdate
	SetAspectToSquarePixels("")

	String wnote=note(Xw), title=""
	title = StringByKey("title",wnote,"=")
	str = StringByKey("sampleName",wnote,"=")
	if (strlen(str)>0)
		title += "\r"+str
	endif
	str = StringByKey("userName",wnote,"=")
	if (strlen(str)>0)
		title += "\r"+str
	endif
	Variable iref=NumberBykey("iref",wnote,"="), mref=NumberBykey("mref",wnote,"=")
	if (iref>=0 || mref>=0)
		title += "\r"
		str = ""
		if (iref>0)
			sprintf str,"iref = %g",iref
			title += str
		endif
		if (mref>0 && mref!=iref)
			title += SelectString(strlen(str),"",",   ")
			sprintf str,"mref = %g",mref
			title += str
		endif
	endif

	str = StringByKey("date",wnote,"=")
	if (strlen(str)>0)
		title += "\r\\Zr075"+ISOtime2niceStr(str)
	endif
	str = StringByKey("xmlFileFull",wnote,"=")
	if (strlen(str))
		title += "\r\\Zr075"+str
	endif
	title = TrimFrontBackWhiteSpace(title)
	TextBox/N=titleText/F=0/B=1/A=LB/X=2/Y=2 title

	if (mref>=0 && mref<DimSize(Xw,0))
		Cursor/P A $Yname mref
		ShowInfo
	endif
	return StringFromList(0,WinList("*",";","WIN:1"))
End



Function RotationBetweenTwoPoints(i1,i2)
	Variable i1,i2
	if (!(i1>=0))
		i1 = strlen(CsrInfo(A)) ? pcsr(A) : NaN
		Wave ww = TraceNameToWaveRef("",StringByKey("TNAME",CsrInfo(A)))
		if (!WaveInClass(ww,"Random3dArrays"))
			return NaN
		endif
	endif
	if (!(i2>=0))
		i2 = strlen(CsrInfo(B)) ? pcsr(B) : NaN
		Wave ww = TraceNameToWaveRef("",StringByKey("TNAME",CsrInfo(B)))
		if (!WaveInClass(ww,"Random3dArrays"))
			return NaN
		endif
	endif
	if (numtype(i1+i2))
		DoAlert 0,"Cursors A & B not on Graph"
		return NaN
	endif
	String fldr=""
	if (WaveExists(ww))
		fldr=GetWavesDataFolder(ww,1)
	endif
	Wave RX=$(fldr+"RX"), RH=$(fldr+"RH"), RF=$(fldr+"RF")
	if (!WaveExists(RX) || !WaveExists(RH) || !WaveExists(RF))
		return NaN
	endif

	Variable printIt = strlen(GetRTStackInfo(2))==0 || stringmatch(GetRTStackInfo(2),"IndexButtonProc")
	Make/N=3/D/FREE R1,R2, axis
	R1 = {RX[i1], RH[i1], RF[i1]}
	R2 = {RX[i2], RH[i2], RF[i2]}

	Make/N=(3,3)/FREE/D rot1, rot2
	rotationMatAboutAxis(R1,NaN,rot1)
	rotationMatAboutAxis(R2,NaN,rot2)
	MatrixOp/O/FREE rot12 = rot2 x rot1^t	// rot12 is rotation from point1 to point2
	Variable angle = axisOfMatrix(rot12,axis)	// total rotation angle between point1 and point2 (degrees)
	if (printIt)
		printf "rotation from point %d to %d is about the axis {XHF} = %s  by  %g¡\r",i1,i2,vec2str(axis),angle
	endif

	// to go from (XYZ) -> (XHF) is a +45¡ rotation about the X-axis
	Make/N=(3,3)/D/FREE rotFrame
	rotFrame[0][0] = 1
	rotFrame[1][1] = cos(45*PI/180)			// this matrix rotates direction of vector by +45¡ (this is changed below)
	rotFrame[2][2] = cos(45*PI/180)
	rotFrame[1][2] = -sin(45*PI/180)
	rotFrame[2][1] = sin(45*PI/180)
	MatrixOP/FREE/O R1 = rotFrame x R1		// rotate from XHF to XYZ
	MatrixOP/FREE/O R2 = rotFrame x R2
	rotationMatAboutAxis(R1,NaN,rot1)
	rotationMatAboutAxis(R2,NaN,rot2)
	MatrixOp/O/FREE rot12 = rot2 x rot1^t	// rot12 is rotation from point1 to point2
	axisOfMatrix(rot12,axis)					// rotation axis in XYZ coordinates

	Make/N=(3,3)/D/FREE ref					// reference reciprocal lattice in beam-line XYZ coordinates
	if (matString2mat(StringByKey("recipRef",note(RX),"="),ref))
		DoAlert 0, "Could not get ref recip lattice from note(RX)"
		print "Could not get ref recip lattice from note(RX)"
		return NaN
	endif

	// find hkl of axis in both grains
//	MatrixOP/FREE gm1 = Inv(rot1) x ref		// gm is the symmetry-reduced recip lattice
//	MatrixOP/FREE gm2 = Inv(rot2) x ref
	MatrixOP/FREE gm1 = rot1 x ref			// gm is the symmetry-reduced recip lattice
	MatrixOP/FREE gm2 = rot2 x ref

	MatrixOP/O/FREE hkl1 = Normalize(Inv(gm1) x axis)
	MatrixOP/O/FREE hkl2 = Normalize(Inv(gm2) x axis)
	Variable h=round(24*hkl1[0]), k=round(24*hkl1[1]), l=round(24*hkl1[2])
	lowestOrderHKL(h,k,l)
	if (printIt)
		Make/FREE hkl0 = {h,k,l}
		printf "in both crystals, the axis is an hkl=%s,    approximately {%s} (off by %.2f¡)\r",vec2str(hkl1),hkl2str(h,k,l),angleVec2Vec(hkl0,hkl1)
	endif
	return angle								// total rotation angle (degrees)
End
//Function RotationBetweenTwoPoints(i1,i2)
//	Variable i1,i2
//	if (!(i1>=0))
//		i1 = strlen(CsrInfo(A)) ? pcsr(A) : NaN
//		Wave ww = TraceNameToWaveRef("",StringByKey("TNAME",CsrInfo(A)))
//		if (!WaveInClass(ww,"Random3dArrays"))
//			return NaN
//		endif
//	endif
//	if (!(i2>=0))
//		i2 = strlen(CsrInfo(B)) ? pcsr(B) : NaN
//		Wave ww = TraceNameToWaveRef("",StringByKey("TNAME",CsrInfo(B)))
//		if (!WaveInClass(ww,"Random3dArrays"))
//			return NaN
//		endif
//	endif
//	if (numtype(i1+i2))
//		DoAlert 0,"Cursors A & B not on Graph"
//		return NaN
//	endif
//	String fldr=""
//	if (WaveExists(ww))
//		fldr=GetWavesDataFolder(ww,1)
//	endif
//	Wave RX=$(fldr+"RX"), RH=$(fldr+"RH"), RF=$(fldr+"RF")
//	if (!WaveExists(RX) || !WaveExists(RH) || !WaveExists(RF))
//		return NaN
//	endif
//
//	Variable printIt = strlen(GetRTStackInfo(2))==0
//	Make/N=3/D/FREE R1,R2, axis
//	R1 = {RX[i1], RH[i1], RF[i1]}
//	R2 = {RX[i2], RH[i2], RF[i2]}
//
//	Make/N=(3,3)/FREE/D rot1, rot2
//	rotationMatAboutAxis(R1,NaN,rot1)
//	rotationMatAboutAxis(R2,NaN,rot2)
//	MatrixOp/O/FREE rot12 = rot2 x rot1^t	// rot12 is rotation from point1 to point2
//	Variable angle = axisOfMatrix(rot12,axis)	// total rotation angle between point1 and point2 (degrees)
//	if (printIt)
//		printf "rotation from point %d to %d is about the axis {XHF} = %s  by  %g¡\r",i1,i2,vec2str(axis),angle
//	endif
//
//	Make/N=(3,3)/D/FREE ref					// reference reciprocal lattice (the RX,RY,RZ are relative to this)
//	if (matString2mat(StringByKey("recipRef",note(RX),"="),ref))
//		DoAlert 0, "Could not get ref recip lattice from note(RX)"
//		print "Could not get ref recip lattice from note(RX)"
//		return NaN
//	endif
//
//	// find hkl of axis in both grains
//	MatrixOP/FREE gm1 = Inv(rot1) x ref		// gm is the symmetry-reduced recip lattice
//	MatrixOP/FREE gm2 = Inv(rot2) x ref
//
//	Variable h,k,l
//	Variable ax=axis[0], ay=(axis[1]-axis[2])/sqrt(2), az=(axis[1]+axis[2])/sqrt(2)	// rotate axis from XHF to XYZ
//	axis = {ax,ay,az}							// axis described in beam-line (XYZ) coordinates
//	MatrixOP/O/FREE hkl1 = Normalize(Inv(gm1) x axis)
//	MatrixOP/O/FREE hkl2 = Normalize(Inv(gm2) x axis)
//	h = round(24*hkl1[0]);	k = round(24*hkl1[1]);	l = round(24*hkl1[2])
//	lowestOrderHKL(h,k,l)
//	String str1 = hkl2str(h,k,l)
//	h = round(24*hkl2[0]);	k = round(24*hkl2[1]);	l = round(24*hkl2[2])
//	lowestOrderHKL(h,k,l)
//	if (printIt)
//		printf "in point 1, axis is an hkl=%s,    approximately (%s)\r",vec2str(hkl1),str1
//		printf "in point 2, axis is an hkl=%s,    approximately (%s)\r",vec2str(hkl2),hkl2str(h,k,l)
//	endif
//	return angle								// total rotation angle (degrees)
//End
//
Static Function matString2mat(str,mat)
	String str
	Wave mat
	Variable a0,a1,a2, b0,b1,b2, c0,c1,c2
	sscanf str, "{{%g,%g,%g}{%g,%g,%g}{%g,%g,%g}}", a0,a1,a2, b0,b1,b2, c0,c1,c2
	if (V_flag==9)
		mat[0][0]=a0;		mat[0][1]=b0;		mat[0][2]=c0
		mat[1][0]=a1;		mat[1][1]=b1;		mat[1][2]=c1
		mat[2][0]=a2;		mat[2][1]=b2;		mat[2][2]=c2
		return 0
	endif
	mat = NaN
	return 1
End

// ******************************  End of Plot Random Points *******************************
// *****************************************************************************************



// *****************************************************************************************
// **********************************  Start of Read XML ***********************************

// after running Load3dRecipLatticesFileXML(), process the loaded waves for viewing
Function/T ProcessLoadedXMLfile(maxAngle,refType,[iref,Xoff,Yoff,Zoff,centerVolume,printIt])
	Variable maxAngle
	Variable refType						// method used for reference orientation, 0=std, 1=average, 2=iref
	Variable iref							// if >= 0, then use this point number as the reference orientation
	Variable Xoff,Yoff,Zoff			// optional offsets
	Variable centerVolume				// if true, center the volume on x,y,z
	Variable printIt
	refType = limit(refType,0,2)==refType ? refType : -1
	iref = ParamIsDefault(iref) ? -1 : iref
	iref = numtype(iref) ? -1 : iref
	refType = iref>=0 && refType<0 ? 2 : refType	// if iref is given, then reftype must be 2 if it is not specified
	Xoff = ParamIsDefault(Xoff) ? 0 : Xoff
	Xoff = numtype(Xoff) ? 0 : Xoff
	Yoff = ParamIsDefault(Yoff) ? 0 : Yoff
	Yoff = numtype(Yoff) ? 0 : Yoff
	Zoff = ParamIsDefault(Zoff) ? 0 : Zoff
	Zoff = numtype(Zoff) ? 0 : Zoff
	centerVolume = ParamIsDefault(centerVolume) ? NaN : centerVolume
	centerVolume = numtype(centerVolume) ? 0 : !(!centerVolume)	// default is False
	printIt = ParamIsDefault(printIt) ? NaN : printIt
	printIt = numtype(printIt) ? strlen(GetRTStackInfo(2))==0 : !(!printIt)

	String rawFldr = GetDataFolder(1)+"raw:"
	if (!ValidRawXMLdataAvailable())
		DoAlert 0, "raw folder \""+rawFldr+"\" not found, first read in the xml file, or it does not have data in it"
		return ""
	endif
	Wave Xsample = $(rawFldr+"Xsample")
	Wave Ysample = $(rawFldr+"Ysample")
	Wave Zsample = $(rawFldr+"Zsample")
	Wave depthRaw = $(rawFldr+"depth")
	Wave totalSumRaw = $(rawFldr+"totalSum")
	Wave sumAboveThresholdRaw = $(rawFldr+"sumAboveThreshold")
	Wave numAboveThresholdRaw = $(rawFldr+"numAboveThreshold")
	Wave NindexedRaw = $(rawFldr+"Nindexed")
	Wave rmsIndexedRaw = $(rawFldr+"rmsIndexed")
	Wave gmRaw = $(rawFldr+"gm")
	Wave std = $(rawFldr+"stdLattice")
	Wave/T imageNamesRaw = $(rawFldr+"imageNames")
	Variable Nraw=DimSize(Xsample,0)					// number of points read in (raw points)
	Variable useSymmetry=NumVarOrDefault(rawFldr+"useSymmetry",0)

	if (!(maxAngle>0) || refType<0 || (refType==2 && !(iref>=0)))
		if (iref<0 && strlen(WinList("*",";","WIN:1"))>0)
			iref = NumberByKey("POINT",CsrInfo(A))	// set iref to cursor A as a default
			iref = iref>=0 ? iref : -1
		endif
		maxAngle = !(maxAngle>0) ? Inf : maxAngle
		Prompt maxAngle, "reject all points with rotation angle greater than this (degree)"
		Prompt refType,"method for picking reference orientation",popup,"Standard;Average;Choose one point (iref)"
		Prompt iref,"point number to use for reference orientation (optional)"
		Prompt Xoff, "X-axis offset, beamline coords (µm)"
		Prompt Yoff, "Y-axis offset, beamline coords (µm)"
		Prompt Zoff, "Z-axis offset, beamline coords (µm)"
		Prompt centerVolume, "Recenter Volume in (X,Y,Z)", popup, "Absolute Positions;Recenter to (XYZ)=(000)"
		centerVolume += 1
		refType +=1
		DoPrompt "max Rotation & reference",maxAngle,Xoff,refType,Yoff,iref,Zoff, centerVolume
		if (V_flag)
			return ""
		endif
		refType -=1
		centerVolume -=1
		printIt = 1
	endif
	if (printIt)
		printf "ProcessLoadedXMLfile(%g, %g",maxAngle,refType
		if (refType==2)
			printf ", iref=%g",iref
		endif
		if (!(Xoff==0 && Yoff==0 && Zoff==0))
			printf ", Xoff=%g, Yoff=%g, Zoff=%g",Xoff,Yoff,Zoff
		endif
		if (centerVolume || !ParamIsDefault(centerVolume))
			printf ", centerVolume=%g",centerVolume
		endif
		printf ")\r"
	endif
	if (!(maxAngle>0))
		DoAlert 0, "ERROR in ProcessLoadedXMLfile(), maxAngle = "+num2str(maxAngle)
		return ""
	elseif (limit(refType,0,2)!=refType)
		DoAlert 0, "ERROR in ProcessLoadedXMLfile(), refType = "+num2str(refType)
		return ""
	elseif (refType==2 && !(iref>=0))
		DoAlert 0, "ERROR in ProcessLoadedXMLfile(), refType==2 and iref = "+num2str(iref)
		return ""
	endif
	if (numtype(Xoff+Yoff+Zoff))
		String str = "ERROR in ProcessLoadedXMLfile()"
		str += SelectString(numtype(Xoff),"",", Xoff="+num2str(Xoff))
		str += SelectString(numtype(Yoff),"",", Yoff="+num2str(Yoff))
		str += SelectString(numtype(Zoff),"",", Zoff="+num2str(Zoff))
		str += " is bad"
		DoAlert 0, str
		return ""
	endif
	iref = refType==2 ? iref : -1						// if not using iref set it to -1

	String noteStr = note(Xsample)
	noteStr = ReplaceNumberByKey("maxAngle",noteStr,maxAngle,"=")
	noteStr = ReplaceNumberByKey("Xoff",noteStr,Xoff,"=")
	noteStr = ReplaceNumberByKey("Yoff",noteStr,Yoff,"=")
	noteStr = ReplaceNumberByKey("Zoff",noteStr,Zoff,"=")
	noteStr = ReplaceNumberByKey("Hoff",noteStr,YZ2H(Yoff,Zoff)	,"=")
	noteStr = ReplaceNumberByKey("Foff",noteStr,YZ2F(Yoff,Zoff),"=")
	noteStr = ReplaceNumberByKey("X0",noteStr,0,"=")
	noteStr = ReplaceNumberByKey("Y0",noteStr,0,"=")
	noteStr = ReplaceNumberByKey("Z0",noteStr,0,"=")
	noteStr = ReplaceNumberByKey("H0",noteStr,0,"=")
	noteStr = ReplaceNumberByKey("F0",noteStr,0,"=")
	noteStr = ReplaceNumberByKey("dX",noteStr,0,"=")
	noteStr = ReplaceNumberByKey("dY",noteStr,-0.707107,"=")	// an inward pointing vector (matches F direction, not -F)
	noteStr = ReplaceNumberByKey("dZ",noteStr,0.707107,"=")
	noteStr = ReplaceNumberByKey("dH",noteStr,0,"=")
	noteStr = ReplaceNumberByKey("dF",noteStr,-1,"=")

	Make/N=(Nraw)/O depth, ZZ,YY, FF, HH, XX
	Make/N=(Nraw)/O rmsIndexed, numAboveThreshold, sumAboveThreshold, totalSum, Nindexed
	Make/N=(3,3,Nraw)/O gm
	Make/N=(Nraw)/T/O imageNames
	Make/N=(Nraw)/O totalAngles, RZ, RY, RF, RH, RX
	Make/N=(Nraw)/O IndexBackTrack = -1
	SetScale d 0,0,"µm", depth, XX, HH, FF, YY, ZZ
	SetScale d 0,0,"¡", totalAngles

	// transform to voxel coordinate in sample, not sample position
	XX = -(Xsample - Xoff)
	YY = -(Ysample - Yoff)
	ZZ = -(Zsample - Zoff) + depth + (numtype(depthRaw) ? 0 : depthRaw)
	depth = depthRaw
	HH = YZ2H(YY,ZZ)
	FF = YZ2F(YY,ZZ)

	if (centerVolume)											// shift so that average center of volume is zero
		Variable X0, Y0, Z0, H0, F0						// average center values
		WaveStats/M=1/Q XX;		X0 = V_avg;		XX -= X0
		WaveStats/M=1/Q YY;		Y0 = V_avg;		HH -= Y0
		WaveStats/M=1/Q ZZ;		Z0 = V_avg;		FF -= Z0
		H0 = YZ2H(Y0,Z0)
		F0 = YZ2F(Y0,Z0)
		if (printIt)
			printf "re-set origin of data to XYZ={%g, %g, %g},  XHF={%g, %g, %g},   so now {0,0,0} is center\r",X0,Y0,Z0,X0,H0,F0
		endif
		noteStr = ReplaceNumberByKey("X0",noteStr,X0,"=")
		noteStr = ReplaceNumberByKey("Y0",noteStr,Y0,"=")
		noteStr = ReplaceNumberByKey("Z0",noteStr,Z0,"=")
		noteStr = ReplaceNumberByKey("H0",noteStr,H0,"=")
		noteStr = ReplaceNumberByKey("F0",noteStr,F0,"=")
		XX -= X0
		YY -= Y0
		ZZ -= Z0
		HH -= H0
		FF -= F0
	endif

	// set rl0 to the reference orientation, how this is done will depend upon 'refType'
	Variable i, angle
	Make/N=(3,3)/D/FREE rl0, gmi, mat3, id33=p==q
	Make/N=3/D/FREE vec3
	String progressWin = ProgressPanelStart("",stop=1,showTime=1,status="starting")	// display a progress bar
	if (refType==2 && (iref>=0 && iref<Nraw))		// a specific point was chosen as the reference
		gmi = gmRaw[p][q][iref]
		if (useSymmetry)										// symmetry reduce this structure
			angle = symReducedRecipLattice(std,id33,gmi,mat3)
		else
			MatrixOp/O/FREE mat3 = gmi x Inv(std)		// gmi = mat3 x std
		endif
		MatrixOp/O/FREE rl0 = mat3 x std				// reference reciprocal lattice

	elseif (refType==1)										// find central orientation by averaging
		Make/N=3/D/FREE vecCenter=0
		MatrixOp/O/FREE stdInv = Inv(std)
		for (i=0;i<Nraw;i+=1)
			if (mod(i,1000) == 0)
				if (ProgressPanelUpdate(progressWin,i/Nraw*100,status="looking for reference, number "+num2istr(i)))	// update progress bar
					break											//   and break out of loop
				endif
			endif
			gmi = gmRaw[p][q][i]
			if (useSymmetry)									// symmetry reduce this structure
				angle = symReducedRecipLattice(std,id33,gmi,mat3)
			else
				MatrixOp/O/FREE mat3 = gmi x stdInv	// gmi = mat3 x std
			endif
			angle = axisOfMatrix(mat3,vec3)				// returned angle (degrees)
			vec3 *= angle
			vecCenter += vec3
		endfor
		if (printIt)
			printf "and finding reference,  execution time = %s\r",Secs2Time(SecondsInProgressPanel(progressWin),5,0)
		endif
		vecCenter /= Nraw										// the central rotation vector
		angle = norm(vecCenter)
		rotationMatAboutAxis(vecCenter,angle,mat3)	// mat3 rotates std to the central orientation (the new reference)
		MatrixOp/O/FREE rl0 = mat3 x std				// reference reciprocal lattice

	elseif (refType==0)										// use standard orientation for reference orientation
		rl0 = std
		mat3 = p==q

	else
		str = "invalid combination refType="+num2str(refType)+",  and iref="+num2str(iref)
		DoAlert 0, str
		print str
		DoWindow/K $progressWin
		return ""
	endif
	sprintf str, "{{%g,%g,%g}{%g,%g,%g}{%g,%g,%g}}",rl0[0][0],rl0[1][0],rl0[2][0],rl0[0][1],rl0[1][1],rl0[2][1],rl0[0][2],rl0[1][2],rl0[2][2]
	noteStr = ReplaceStringByKey("recipRef",noteStr,str,"=")// save reference reicp lattice, as column vectors
	if (printIt)
		print "using a reference reciprocal lattice of: ",StringFromList(refType,"Standard Orientation;Average Orientation;Orientation from point ")+SelectString(refType==2,"",num2str(iref))
		printf "%.2f\t%.2f\t%.2f (1/nm)\r",rl0[0][0],rl0[0][1],rl0[0][2]
		printf "%.2f\t%.2f\t%.2f\r",rl0[1][0],rl0[1][1],rl0[1][2]
		printf "%.2f\t%.2f\t%.2f\r",rl0[2][0],rl0[2][1],rl0[2][2]
		if (refType==2 && (iref>=0 && iref<Nraw))		// a valid reference point
			printf "from iref=%d, at position (%g, %g, %g)\r",iref,XX[iref],HH[iref],FF[iref]
		endif

		Make/FREE/N=3/D hkl = {0,1,-1}
		MatrixOP/O/FREE hkl = Inv(rl0) x hkl
		Variable maxH = max(-WaveMin(hkl),WaveMax(hkl))
		hkl = round(hkl*24/maxH)
		Variable h=hkl[0], k=hkl[1], l=hkl[2]
		lowestOrderHKL(h,k,l)
		hkl = {h,k,l}
		printf "hkl of surface normal (reference) = %s\r",vec2str(hkl)
	endif

	// use current reference orientation (rl0) to compute all Rodriques vectors
	Make/N=(3,3)/D/FREE rot33
	Variable N=0, mref=NaN									// use this for rejecting rotations greater than maxAngle
	for (i=0;i<Nraw;i+=1)									// using rl0, find the Rodriques vectors
		if (mod(i,1000) == 0)
			if (ProgressPanelUpdate(progressWin,i/Nraw*100,status="setting rotation, number "+num2istr(i)))	// update progress bar
				break												//   and break out of loop
			endif
		endif
		gmi = gmRaw[p][q][i]
		if (useSymmetry)										// symmetry reduce this structure
			angle = symReducedRecipLattice(std,mat3,gmi,rot33)	// angle between go and gm,  gmi = rot33 x rl0 = rot33 x (mat3 x std)
		else
			MatrixOp/O/FREE rot33 = gmi x Inv(rl0)	// gmi = rot33 x rl0
		endif
		angle = axisOfMatrix(rot33,vec3)				// returned angle (degrees)
		if (numtype(angle) && i==0)
			str = "Cannot find axis of rotation matrix, Probably the lattice parameters in Igor do not match those in the xml file.\r  Check the 'Xtal' tab"
			DoAlert 0, str
			print str
			DoWindow/K $progressWin
			return ""
		endif
		vec3 *= tan(angle*PI/180/2)						// this is now the Rodriques vector
		if (angle<=maxAngle)								// keep this one
			depth[N] = depth[i]
			XX[N] = XX[i]
			YY[N] = YY[i]
			ZZ[N] = ZZ[i]
			HH[N] = HH[i]
			FF[N] = FF[i]
			totalSum[N] = totalSumRaw[i]
			sumAboveThreshold[N] = sumAboveThresholdRaw[i]
			numAboveThreshold[N] = numAboveThresholdRaw[i]
			Nindexed[N] = NindexedRaw[i]
			rmsIndexed[N] = rmsIndexedRaw[i]
			gm[][][N] = gmRaw[p][q][i]
			imageNames[N] = imageNamesRaw[i]
			RX[N] = vec3[0]
			RY[N] = vec3[1]
			RZ[N] = vec3[2]
			totalAngles[N] = angle
			IndexBackTrack[N] = i
			if (i==iref)
				mref = N
			endif
			N += 1
		endif
	endfor
	if (printIt)
		printf "total  execution time = %s\r",Secs2Time(SecondsInProgressPanel(progressWin),5,0)
	endif
	DoWindow/K $progressWin
	Redimension/N=(N) depth,XX,YY,ZZ,HH,FF,totalSum,sumAboveThreshold,numAboveThreshold,Nindexed,rmsIndexed
	Redimension/N=(N) RX,RY,RZ, RH,RF, totalAngles
	Redimension/N=(N) imageNames
	Redimension/N=(-1,-1,N) gm
	if (refType==2 && (iref>=0 && iref<N))			// a valid reference point
		if (printIt)
			printf "requested iref=%d,  this corresponds to point %d\r",iref,mref
		endif
		noteStr = ReplaceNumberByKey("iref",noteStr,iref,"=")
		noteStr = ReplaceNumberByKey("mref",noteStr,mref,"=")
	endif

	// to go from (XYZ) -> (XHF) is a +45¡ rotation about the X-axis
	Make/N=(3,3)/D/FREE rotFrame=(p==q)				// rotates about X-axis (rotates vector not frame) to get from beam-line coordinates to (XHF)
	Variable cosTheta = NumVarOrDefault("root:Packages:geometry:cosThetaWire",cos(PI/4))
	Variable sinTheta = NumVarOrDefault("root:Packages:geometry:sinThetaWire",sin(PI/4))
	rotFrame[1][1] = cosTheta								// this matrix rotates direction of vector by +45¡ (this is changed below)
	rotFrame[2][2] = cosTheta
	rotFrame[1][2] = -sinTheta
	rotFrame[2][1] = sinTheta
	MatrixOp/O/FREE rotFrame = Inv(rotFrame)		// rotFrame now rotates the description of a vector from (XYZ) frame to (XHF) frame
	Make/N=(N,3)/D/FREE Rxyz
	Rxyz[][0] = RX[p]
	Rxyz[][1] = RY[p]
	Rxyz[][2] = RZ[p]
	MatrixOP/FREE/O Rxhf = (rotFrame x Rxyz^t)^t	// rotate from beam line frame to XHF frame (outward surface normal is Z)
	RH = Rxhf[p][1]
	RF = Rxhf[p][2]
	WaveClear Rxyz, Rxhf

	WaveStats/M=1/Q totalAngles
	Variable maxAngleFound = V_max
	if (printIt)
		printf "rotation angles range from:   %.3g¡ at %d  to  %.3g¡ at %d\r",V_min,V_minloc,V_max,V_maxloc
		printf "remaining in Igor data folder  '%s'\r",GetDataFolder(1)
	endif

	Variable Xlo=WaveMin(XX), Xhi=WaveMax(XX)
	Variable Ylo=WaveMin(YY), Yhi=WaveMax(YY)
	Variable Zlo=WaveMin(ZZ), Zhi=WaveMax(ZZ)
	Variable Hlo=WaveMin(HH), Hhi=WaveMax(HH)
	Variable Flo=WaveMin(FF), Fhi=WaveMax(FF)
	if (printIt)
		print "  Range of sample in Sample Coordinates"
		printf "X = [%g, %g] µm\r",Xlo,Xhi
		printf "Y = [%g, %g] µm\r",Ylo,Yhi
		printf "Z = [%g, %g] µm\r",Zlo,Zhi
		printf "H = [%g, %g] µm\r",Hlo,Hhi
		printf "F = [%g, %g] µm\r",Flo,Fhi
	endif
	noteStr = ReplaceNumberByKey("Xlo",noteStr,Xlo,"=")
	noteStr = ReplaceNumberByKey("Xhi",noteStr,Xhi,"=")
	noteStr = ReplaceNumberByKey("Ylo",noteStr,Ylo,"=")
	noteStr = ReplaceNumberByKey("Yhi",noteStr,Yhi,"=")
	noteStr = ReplaceNumberByKey("Zlo",noteStr,Zlo,"=")
	noteStr = ReplaceNumberByKey("Zhi",noteStr,Zhi,"=")
	noteStr = ReplaceNumberByKey("Hlo",noteStr,Hlo,"=")
	noteStr = ReplaceNumberByKey("Hhi",noteStr,Hhi,"=")
	noteStr = ReplaceNumberByKey("Flo",noteStr,Flo,"=")
	noteStr = ReplaceNumberByKey("Fhi",noteStr,Fhi,"=")

	Note/K XX,noteStr
	Note/K YY,noteStr
	Note/K ZZ,noteStr
	Note/K HH,noteStr
	Note/K FF,noteStr
	Note/K depth,noteStr
	Note/K RX,noteStr
	Note/K RY,noteStr
	Note/K RZ,noteStr
	Note/K RH,noteStr
	Note/K RF,noteStr
	Note/K totalSum,noteStr
	Note/K sumAboveThreshold,noteStr
	Note/K numAboveThreshold,noteStr
	Note/K Nindexed,noteStr
	Note/K rmsIndexed,noteStr
	Note/K totalAngles,noteStr
	Note/K imageNames,noteStr
	Note/K gm,ReplaceStringByKey("waveClass", noteStr, "Random3dArraysGm","=")

	Wave rotrgb=$makeRGBJZT(RX,RH,RF,maxAngleFound)		// rgb suitable for a 2-D plot
	Variable threeD=0
	WaveStats/Q/M=1 XX
	threeD += (V_max-V_min) > 0.1
	WaveStats/Q/M=1 HH
	threeD += (V_max-V_min) > 0.1
	WaveStats/Q/M=1 FF
	threeD += (V_max-V_min) > 0.1
	threeD = (threeD==3)
	if (threeD)
		Make/N=(N,3)/O xyz
		SetScale d 0,0,WaveUnits(XX,-1), xyz
		xyz[][0] = XX[p]
		xyz[][1] = HH[p]
		xyz[][2] = FF[p]
		MakeGizmocubeCorners(xyz)
		Make/N=(N,4)/O rgba
		Make/N=(N)/O/FREE brite
		Variable hi=WaveMax(totalsum)/4
		brite = sqrt(limit(totalsum/hi,0.09,1))
		// brite = sqrt(limit(totalsum/100e3,0.09,1))
		rgba[][0,2] = rotrgb[p][q]/65535
		rgba[][3] = brite[p]
		Note/K xyz,ReplaceStringByKey("waveClass",noteStr,"Random3dArraysXYZ","=")

		String note2= ReplaceStringByKey("waveClass",noteStr,"Random3dArraysRGBA","=")
		note2= ReplaceStringByKey("source",noteStr,"RxRyRz","=")
		String title=StringByKey("title",note(RX),"=")
		title = "RxRyRz"+SelectString(strlen(title),"",", ")+title
		note2= ReplaceStringByKey("title",noteStr,title,"=")
		Note/K rgba,note2
	endif
	return noteStr
End



Function/T Load3dRecipLatticesFileXML(FullFileName,[printIt])
	String FullFileName
	Variable printIt
	printIt = ParamIsDefault(printIt) ? NaN : printIt
	printIt = numtype(printIt) ? strlen(GetRTStackInfo(2))==0 : !(!printIt)

	Variable f
	Open/R/Z=2/M="pick the XML file"/F=XMLfilters f as FullFileName
	printIt = StringMatch(FullFileName,S_fileName) ? printIt : 1
	FullFileName = S_fileName
	if (V_flag)
		return ""
	endif
	if (printIt)
		printf "Load3dRecipLatticesFileXML(\"%s\")\r",FullFileName
	endif

	String pathSep
	if (stringmatch(IgorInfo(2),"Macintosh"))
		pathSep = ":"
	elseif (stringmatch(IgorInfo(2),"Windows"))
		pathSep = "\\"
	else
		Abort "This is neither Mac nor Win, what is this computer?"
	endif

	GetFileFolderInfo/Q/Z=1 FullFileName
	if (!V_isFile)
		return ""
	endif
	String noteStr = ReplaceStringByKey("waveClass","","Random3dArrays","=")
	noteStr = ReplaceStringByKey("xmlFileFull",noteStr,FullFileName,"=")	// add full file name to note string
	noteStr = ReplaceStringByKey("localFile",noteStr,ParseFilePath(0,FullFileName,":",1,0),"=")	// add file name to note string
	String fldrSav= GetDataFolder(1)
	String fldrName = ParseFilePath(3,FullFileName, pathSep, 0, 0)
	fldrName = CleanupName(fldrName, 0)
	if (stringmatch(fldrName,"raw"))
		fldrName = UniqueName(fldrName,11,0)
	endif
	if (CheckName(fldrName,11))
		DoAlert 2, "This data folder already exists, delete the folder and re-load the data?"
		if (V_flag==1)
			KillDataFolder $fldrName
		elseif (V_flag==2)
			fldrName = UniqueName(fldrName,11,0)
		else
			SetDataFolder fldrSav
			return ""
		endif
	endif

	if (StringMatch(GetDataFolder(0),fldrName))
		DoAlert 2,"In folder:\r  \""+fldrName+"\",\rreload and overwrite data from xml file?"
		if (V_flag!=1)
			return ""
		endif
	else
		NewDataFolder/O/S $fldrName
	endif

	fldrName = GetDataFolder(1)
	NewDataFolder/O/S raw
	noteStr = ReplaceStringByKey("fldrName",noteStr,fldrName,"=")

	STRUCT crystalStructure xtal
	if (LatticeParametersFromXML(FullFileName,xtal))				// fill the lattice from first occurance in xml file
		if (FillCrystalStructDefault(xtal))		// if no lattice in xml file, use current values
			DoAlert 0, "no crystal structure found"
			return ""
		endif
	endif
	LatticeSym#ForceLatticeToStructure(xtal)	// mostly calling this to call reMakeAtomXYZs(xtal)
	initSymmetryOperations()							// initialize all symmetry operations
	Variable/G useSymmetry=strlen(LatticeSym#MakeSymmetryOps(xtal))>0	// make a wave with the symmetry operation
	Make/N=(3,3)/O/D stdLattice						// reciprocal lattice in standard orientation (used for sym reduction)
	stdLattice[0][0] = xtal.as0 ;	stdLattice[0][1] = xtal.bs0 ;	stdLattice[0][2] = xtal.cs0
	stdLattice[1][0] = xtal.as1 ;	stdLattice[1][1] = xtal.bs1 ;	stdLattice[1][2] = xtal.cs1
	stdLattice[2][0] = xtal.as2 ;	stdLattice[2][1] = xtal.bs2 ;	stdLattice[2][2] = xtal.cs2

	Variable Nalloc=1000
	Make/N=(Nalloc)/O depth, Zsample,Ysample,Xsample
	Make/N=(Nalloc)/O sumAboveThreshold, totalSum		// total sum of pixels, and sum of pixels over threshold
	Make/N=(Nalloc)/O numAboveThreshold			// number of pixels that exceed threshold (by any amount)
	Make/N=(Nalloc)/O Nindexed						// number of spots indexed in first pattern
	Make/N=(Nalloc)/O rmsIndexed						// rms error indexing each pattern (degree)
	Make/N=(3,3,Nalloc)/D/O gm
	SetScale d 0,0,"µm", Xsample,Ysample,Zsample,depth
	Make/N=(Nalloc)/O/T imageNames

	Variable Yn,Zn
	String pattern, patternKeyVals
	String svec
	Variable as0,as1,as2, bs0,bs1,bs2, cs0,cs1,cs2
	Variable first = 1									// flags first time through
	Variable value, i
	String step
	Variable i0,i1											// mark start and end of current   <step ></step>
	Variable i0start=0									// where to start searching for "<step "
	Variable N=0											// counts number of points read in
	Variable bytesRead									// keeps track of how many bytes read in (should not exceed fileLen)
	Variable size=200*1024								// size of a typical read
	Variable N0												// number of spots indexed in first pattern
	Variable rms											// rms error indexing a pattern

	// start reading here
	String progressWin = ProgressPanelStart("",stop=1,showTime=1,status="starting")	// display a progress bar
	Open/R/Z=1 f as FullFileName
	FStatus f
	Variable fileLen = V_logEOF						// length of file in bytes
	size = min(size,fileLen)
	String buf = PadString("",size,0x20)			// a buffer for working space (of length size bytes)
	String bufRead, str									// a buffer for reading more from the file (read up to this much each time)
	FBinRead f, buf										// initial read
	bytesRead = strlen(buf)
	do
		// read in data from xml file
		if (mod(N,1000) == 0)
			if (ProgressPanelUpdate(progressWin,bytesRead/fileLen*100,status="reading number "+num2istr(N)))	// update progress bar
				break											//   and break out of loop
			endif
		endif
		i0 = strsearch(buf,"<step ",i0start,2)
		if (i0<0)
			size = min(size,fileLen-bytesRead)		// number of bytes to read
			bufRead = PadString("",size,0x20)
			FBinRead f, bufRead							// could not find start of step, read more
			bytesRead += strlen(bufRead)
			buf = buf[i0,Inf] + bufRead
			i0 = strsearch(buf,"<step ",i0start,2)		// search again for start of step tag
			if (i0<0)
				break											// give up, all done
			endif
		endif
		i1 = strsearch(buf,"</step>",i0,2)
		if (i1<0)
			size = min(size,fileLen-bytesRead)		// number of bytes to read
			bufRead = PadString("",size,0x20)
			FBinRead f, bufRead							// could not find end of step, read more
			bytesRead += strlen(bufRead)
			buf = buf[i0,Inf] + bufRead
			i0 = 0
			i1 = strsearch(buf,"</step>",i0,2)		// search again for end of step tag
			if (i1<0)
				break											// give up, all done
			endif
		endif
		i1 += 6

		// this step is bracketed by [i0,i1], process it:
		step = buf[i0,i1]
		i0start = i1 + 1									// where to start searching for next start tag

		if (first)
			first = 0
			noteStr = ReplaceStringByKey("title",noteStr,xmlTagContents("title",step),"=")
			noteStr = ReplaceStringByKey("sampleName",noteStr,xmlTagContents("sampleName",step),"=")
			noteStr = ReplaceStringByKey("userName",noteStr,xmlTagContents("userName",step),"=")
			noteStr = ReplaceNumberByKey("scanNum",noteStr,str2num(xmlTagContents("scanNum",step)),"=")
			noteStr = ReplaceStringByKey("date",noteStr,xmlTagContents("date",step),"=")
			value = str2num(xmlTagContents("hutchTemperature",step))
			if (value != 0)
				noteStr = ReplaceNumberByKey("hutchTemperature",noteStr,value,"=")
			endif
			value = str2num(xmlTagContents("beamBad",step))
			if (value)
				noteStr = ReplaceNumberByKey("beamBad",noteStr,value,"=")
			endif
			str = xmlTagContents("CCDshutter",step)
			if (!stringmatch(str,"out"))
				noteStr = ReplaceStringByKey("CCDshutter",noteStr,"in","=")
			endif
			str = xmlTagContents("monoMode",step)
			noteStr = ReplaceStringByKey("monoMode",noteStr,str,"=")
			if (stringmatch(str,"monochromatic"))
				noteStr = ReplaceNumberByKey("energy",noteStr,str2num(xmlTagContents("energy",step)),"=")
			endif
		endif
		patternKeyVals = multiIndex#xmlTagKeyVals("pattern",step)
		N0 = NumberByKey("Nindexed", patternKeyVals,"=")
		rms = NumberByKey("rms_error", patternKeyVals,"=")
		if (NumberByKey("num", patternKeyVals,"=")!=0 || N0<NindexedMIN)
			//if (strlen(patternKeyVals))
			//	print xmlTagContents("inputImage",step)
			//	print patternKeyVals
			//endif
			continue											// skip this image, it was not adequately indexed
		endif
		pattern = xmlTagContents("pattern",step)
		svec = xmlTagContents("astar",pattern)
		sscanf svec,"%g %g %g",as0,as1,as2
		svec = xmlTagContents("bstar",pattern)
		sscanf svec,"%g %g %g",bs0,bs1,bs2
		svec = xmlTagContents("cstar",pattern)
		sscanf svec,"%g %g %g",cs0,cs1,cs2

		if (N>=Nalloc)										// need longer arrays
			Nalloc += 1000
			Redimension/N=(Nalloc) Xsample,Ysample,Zsample,depth
			Redimension/N=(Nalloc) totalSum, sumAboveThreshold, numAboveThreshold
			Redimension/N=(Nalloc) Nindexed
			Redimension/N=(Nalloc) rmsIndexed
			Redimension/N=(Nalloc) imageNames
			Redimension/N=(3,3,Nalloc) gm
		endif
		depth[N] = str2num(xmlTagContents("depth",step))
		Xsample[N] = str2num(xmlTagContents("Xsample",step))
		Ysample[N] = str2num(xmlTagContents("Ysample",step))
		Zsample[N] = str2num(xmlTagContents("Zsample",step))
		Nindexed[N] = N0
		rmsIndexed[N] = rms
		totalSum[N] = str2num(xmlTagContents("totalSum",step))
		sumAboveThreshold[N] = str2num(xmlTagContents("sumAboveThreshold",step))
		numAboveThreshold[N] = str2num(xmlTagContents("numAboveThreshold",step))
		imageNames[N] = xmlTagContents("inputImage",step)
		gm[0][0][N] = as0;		gm[0][1][N] = bs0;	gm[0][2][N] = cs0	// save measured reciprocal lattice for this point
		gm[1][0][N] = as1;		gm[1][1][N] = bs1;	gm[1][2][N] = cs1
		gm[2][0][N] = as2;		gm[2][1][N] = bs2;	gm[2][2][N] = cs2
		N += 1
	while(1)
	Close f
	if (printIt)
		printf "  execution time for reading = %s\r",Secs2Time(SecondsInProgressPanel(progressWin),5,0)
	endif
	DoWindow/K $progressWin

	Redimension/N=(-1,-1,N) gm						// trim to actual length
	Redimension/N=(N) depth,Xsample,Ysample,Zsample,totalSum,sumAboveThreshold,numAboveThreshold,Nindexed,rmsIndexed
	Redimension/N=(N) imageNames
	SetDataFolder $fldrName
	printf "remaining in Igor data folder  '%s'\r",fldrName
	// done reading

	Duplicate/FREE Ysample, Hsample, Fsample
	Hsample = YZ2H(Ysample[p],Zsample[p])			// H =  Y*sin(angle) + Z*cos(angle))	Yn and Zn are sample system from the file
	Fsample = YZ2F(Ysample[p],Zsample[p])			// F = -Y*cos(angle) + Z*sin(angle)
	Variable Xlo=WaveMin(Xsample), Xhi=WaveMax(Xsample)
	Variable Ylo=WaveMin(Ysample), Yhi=WaveMax(Ysample)
	Variable Zlo=WaveMin(Zsample), Zhi=WaveMax(Zsample)
	Variable Hlo=WaveMin(Hsample), Hhi=WaveMax(Hsample)
	Variable Flo=WaveMin(Fsample), Fhi=WaveMax(Fsample)
	if (printIt)
		printf "read in %d steps\r",N
		print "   Positioner ranges:"
		printf "Xsample = [%g, %g] µm\r",Xlo,Xhi
		printf "Ysample = [%g, %g] µm\r",Ylo,Yhi
		printf "Zsample = [%g, %g] µm\r",Zlo,Zhi
		printf "Hsample = [%g, %g] µm\r",Hlo,Hhi
		printf "Fsample = [%g, %g] µm\r",Flo,Fhi
	endif

	Note/K Xsample,noteStr
	Note/K Ysample,noteStr
	Note/K Zsample,noteStr
	Note/K depth,noteStr
	Note/K totalSum,noteStr
	Note/K sumAboveThreshold,noteStr
	Note/K numAboveThreshold,noteStr
	Note/K Nindexed,noteStr
	Note/K rmsIndexed,noteStr
	Note/K imageNames,noteStr
	Note/K gm,ReplaceStringByKey("waveClass", noteStr, "Random3dArraysGm","=")

	return noteStr
End


Static Function ValidRawXMLdataAvailable()
	Variable valid = DataFolderExists(":raw")
	if (!valid)
		return 0
	endif
	valid = valid && Exists(":raw:stdLattice")==1
	valid = valid && Exists(":raw:Xsample")==1
	valid = valid && Exists(":raw:Ysample")==1
	valid = valid && Exists(":raw:Zsample")==1
	valid = valid && Exists(":raw:depth")==1
	valid = valid && Exists(":raw:totalSum")==1
	valid = valid && Exists(":raw:sumAboveThreshold")==1
	valid = valid && Exists(":raw:numAboveThreshold")==1
	valid = valid && Exists(":raw:Nindexed")==1
	valid = valid && Exists(":raw:rmsIndexed")==1
	valid = valid && Exists(":raw:gm")==1
	valid = valid && Exists(":raw:imageNames")==1
	valid = valid && Exists(":raw:useSymmetry")==2
	if (!valid)
		return 0
	endif

	Wave Xsample = $":raw:Xsample"
	valid = valid && DimSize(Xsample,0)>1
	return valid
End
















Function/T Load3dRecipLatticesFileXML_OLD(FullFileName,maxAngle,refType,[iref,Xoff,Yoff,Zoff])
	String FullFileName
	Variable maxAngle
	Variable refType						// method used for reference orientation, 0=std, 1=average, 2=iref
	Variable iref							// if >= 0, then use this point number as the reference orientation
	Variable Xoff,Yoff,Zoff					// optional offsets
	refType = limit(refType,0,2)==refType ? refType : -1
	iref = ParamIsDefault(iref) ? -1 : iref
	iref = numtype(iref) ? -1 : iref
	refType = iref>=0 && refType<0 ? 2 : refType	// if iref is given, then reftype must be 2 if it is not specified
	Xoff = ParamIsDefault(Xoff) ? 0 : Xoff
	Xoff = numtype(Xoff) ? 0 : Xoff
	Yoff = ParamIsDefault(Yoff) ? 0 : Yoff
	Yoff = numtype(Yoff) ? 0 : Yoff
	Zoff = ParamIsDefault(Zoff) ? 0 : Zoff
	Zoff = numtype(Zoff) ? 0 : Zoff

	Variable f
	Open/R/Z=2/M="pick the XML file"/T=".xml"  f as FullFileName
	FullFileName = S_fileName
	if (V_flag)
		return ""
	endif
	if (!(maxAngle>0) || refType<0 || (refType==2 && !(iref>=0)))
		maxAngle = !(maxAngle>0) ? Inf : maxAngle
		Prompt maxAngle, "reject all points with rotation angle greater than this (degree)"
		Prompt refType,"method for picking reference orientation",popup,"Standard;Average;Choose one point (iref)"
		Prompt iref,"point number to use for reference orientation (optional)"
		Prompt Xoff, "X-axis offset, beamline coords (µm)"
		Prompt Yoff, "Y-axis offset, beamline coords (µm)"
		Prompt Zoff, "Z-axis offset, beamline coords (µm)"
		refType +=1
		DoPrompt "max Rotation & reference",maxAngle,Xoff,refType,Yoff,iref,Zoff
		if (V_flag)
			return ""
		endif
		refType -=1
	endif
	printf "¥Load3dRecipLatticesFileXML(\"%s\",%g,%g",FullFileName,maxAngle,refType
	if (refType==2)
		printf ", iref=%g",iref
	endif
	if (!(Xoff==0 && Yoff==0 && Zoff==0))
		printf ", Xoff=%g, Yoff=%g, Zoff=%g",Xoff,Yoff,Zoff
	endif
	printf ")\r"
	if (!(maxAngle>0))
		DoAlert 0, "ERROR in Load3dRecipLatticesFileXML(), maxAngle = "+num2str(maxAngle)
		return ""
	elseif (limit(refType,0,2)!=refType)
		DoAlert 0, "ERROR in Load3dRecipLatticesFileXML(), refType = "+num2str(refType)
		return ""
	elseif (refType==2 && !(iref>=0))
		DoAlert 0, "ERROR in Load3dRecipLatticesFileXML(), refType==2 and iref = "+num2str(iref)
		return ""
	endif
	if (numtype(Xoff+Yoff+Zoff))
		String str = "ERROR in Load3dRecipLatticesFileXML()"
		str += SelectString(numtype(Xoff),"",", Xoff="+num2str(Xoff))
		str += SelectString(numtype(Yoff),"",", Yoff="+num2str(Yoff))
		str += SelectString(numtype(Zoff),"",", Zoff="+num2str(Zoff))
		str += " is bad"
		DoAlert 0, str
		return ""
	endif

	String pathSep
	if (stringmatch(IgorInfo(2),"Macintosh"))
		pathSep = ":"
	elseif (stringmatch(IgorInfo(2),"Windows"))
		pathSep = "\\"
	else
		Abort "This is neither Mac nor Win, what is this computer?"
	endif
	Variable refNum
	Open/R/Z=1 refNum as FullFileName					// test that file exists
	if (strlen(S_fileName)<1 || V_flag)
		return ""
	endif
	Close refNum

	String noteStr = ReplaceStringByKey("waveClass","","Random3dArrays","=")
	noteStr = ReplaceStringByKey("xmlFileFull",noteStr,FullFileName,"=")	// add full file name to note string
	noteStr = ReplaceStringByKey("localFile",noteStr,ParseFilePath(0,FullFileName,":",1,0),"=")	// add file name to note string
	noteStr = ReplaceNumberByKey("maxAngle",noteStr,maxAngle,"=")
	String fldrSav= GetDataFolder(1)
	String fldrName = ParseFilePath(3,FullFileName, pathSep, 0, 0)
	fldrName = CleanupName(fldrName, 0)
	if (stringmatch(fldrName,"raw"))
		fldrName = UniqueName(fldrName,11,0)
	endif
	if (CheckName(fldrName,11))
		DoAlert 2, "This data folder already exists, delete the folder and re-load the data?"
		if (V_flag==1)
			KillDataFolder $fldrName
		elseif (V_flag==2)
			fldrName = UniqueName(fldrName,11,0)
		else
			SetDataFolder fldrSav
			return ""
		endif
	endif
	NewDataFolder/O/S $fldrName
	noteStr = ReplaceStringByKey("fldrName",noteStr,GetDataFolder(1),"=")

	STRUCT crystalStructure xtal
	if (LatticeParametersFromXML(FullFileName,xtal))	// fill the lattice from first occurance in xml file
		if (FillCrystalStructDefault(xtal))					// if no lattice in xml file, use current values
			DoAlert 0, "no crystal structure found"
			return ""
		endif
	endif
	LatticeSym#ForceLatticeToStructure(xtal)				// mostly calling this to call reMakeAtomXYZs(xtal)
	initSymmetryOperations()								// initialize all symmetry operations
	Variable useSymmetry									// flag, know how to reduce this symmetry
	useSymmetry = strlen(LatticeSym#MakeSymmetryOps(xtal))>0	// make a wave with the symmetry operation

	Make/N=(3,3)/O/D std, id33=p==q						// id33 is 3x3 identity matrix
	std[0][0] = xtal.as0 ;	std[0][1] = xtal.bs0 ;	std[0][2] = xtal.cs0	// reciprocal lattice in standard orientation (needed for sym reduction)
	std[1][0] = xtal.as1 ;	std[1][1] = xtal.bs1 ;	std[1][2] = xtal.cs1
	std[2][0] = xtal.as2 ;	std[2][1] = xtal.bs2 ;	std[2][2] = xtal.cs2

	Make/N=(3,3)/O/D rl0, gmi, mat3
	Make/N=3/O/D vec3

	Variable Nalloc=1000
	Make/N=(Nalloc)/O FF,HH,XX, depth
	Make/N=(Nalloc)/O sumAboveThreshold, totalSum		// total sum of pixels, and sum of pixels over threshold
	Make/N=(Nalloc)/O numAboveThreshold					// number of pixels that exceed threshold (by any amount)
	Make/N=(Nalloc)/O Nindexed							// number of spots indexed in first pattern
	Make/N=(Nalloc)/O rmsIndexed							// rms error indexing each pattern (degree)
	Make/N=(3,3,Nalloc)/D/O gm
	SetScale d 0,0,"µm", XX,HH,FF,depth
	Make/N=(Nalloc)/O/T imageNames

	Variable Yn,Zn
	String pattern, patternKeyVals
	String svec
	Variable as0,as1,as2, bs0,bs1,bs2, cs0,cs1,cs2
	Variable angle
	Variable first = 1										// flags first time through
	Variable value, i
	String step
	Variable i0,i1											// mark start and end of current   <step ></step>
	Variable i0start=0										// where to start searching for "<step "
	Variable N=0											// counts number of points read in
	Variable bytesRead										// keeps track of how many bytes read in (should not exceed fileLen)
	Variable size=200*1024								// size of a typical read
	Variable N0												// number of spots indexed in first pattern
	Variable rms											// rms error indexing a pattern

	// start reading here
	String progressWin = ProgressPanelStart("",stop=1,showTime=1,status="starting")	// display a progress bar
	Open/R/Z=1 refNum as FullFileName
	FStatus refNum
	Variable fileLen = V_logEOF								// length of file in bytes
	size = min(size,fileLen)
	String buf = PadString("",size,0x20)					// a buffer for working space (of length size bytes)
	String bufRead											// a buffer for reading more from the file (read up to this much each time)
	FBinRead refNum, buf									// initial read
	bytesRead = strlen(buf)
	do
		// read in data from xml file
		if (mod(N,400) == 0)
			if (ProgressPanelUpdate(progressWin,bytesRead/fileLen*100,status="reading number "+num2istr(N)))	// update progress bar
				break										//   and break out of loop
			endif
		endif
		i0 = strsearch(buf,"<step ",i0start,2)
		if (i0<0)
			size = min(size,fileLen-bytesRead)				// number of bytes to read
			bufRead = PadString("",size,0x20)
			FBinRead refNum, bufRead						// could not find start of step, read more
			bytesRead += strlen(bufRead)
			buf = buf[i0,Inf] + bufRead
			i0 = strsearch(buf,"<step ",i0start,2)			// search again for start of step tag
			if (i0<0)
				break										// give up, all done
			endif
		endif
		i1 = strsearch(buf,"</step>",i0,2)
		if (i1<0)
			size = min(size,fileLen-bytesRead)				// number of bytes to read
			bufRead = PadString("",size,0x20)
			FBinRead refNum, bufRead						// could not find end of step, read more
			bytesRead += strlen(bufRead)
			buf = buf[i0,Inf] + bufRead
			i0 = 0
			i1 = strsearch(buf,"</step>",i0,2)				// search again for end of step tag
			if (i1<0)
				break										// give up, all done
			endif
		endif
		i1 += 6

		// this step is bracketed by [i0,i1], process it:
		step = buf[i0,i1]
		i0start = i1 + 1									// where to start searching for next start tag

		if (first)
			first = 0
			noteStr = ReplaceNumberByKey("Xoff",noteStr,Xoff,"=")
			noteStr = ReplaceNumberByKey("Yoff",noteStr,Yoff,"=")
			noteStr = ReplaceNumberByKey("Zoff",noteStr,Zoff,"=")
			noteStr = ReplaceNumberByKey("Hoff",noteStr,YZ2H(Yoff,Zoff)	,"=")
			noteStr = ReplaceNumberByKey("Foff",noteStr,YZ2F(Yoff,Zoff),"=")
			noteStr = ReplaceNumberByKey("X0",noteStr,0,"=")
			noteStr = ReplaceNumberByKey("Y0",noteStr,0,"=")
			noteStr = ReplaceNumberByKey("Z0",noteStr,0,"=")
			noteStr = ReplaceNumberByKey("H0",noteStr,0,"=")
			noteStr = ReplaceNumberByKey("F0",noteStr,0,"=")
			noteStr = ReplaceNumberByKey("dX",noteStr,0,"=")
			noteStr = ReplaceNumberByKey("dY",noteStr,-0.707107,"=")		// an inward pointing vector (matches F direction, not -F)
			noteStr = ReplaceNumberByKey("dZ",noteStr,0.707107,"=")
			noteStr = ReplaceNumberByKey("dH",noteStr,0,"=")
			noteStr = ReplaceNumberByKey("dF",noteStr,-1,"=")
			noteStr = ReplaceStringByKey("title",noteStr,xmlTagContents("title",step),"=")
			noteStr = ReplaceStringByKey("sampleName",noteStr,xmlTagContents("sampleName",step),"=")
			noteStr = ReplaceStringByKey("userName",noteStr,xmlTagContents("userName",step),"=")
			noteStr = ReplaceNumberByKey("scanNum",noteStr,str2num(xmlTagContents("scanNum",step)),"=")
			noteStr = ReplaceStringByKey("date",noteStr,xmlTagContents("date",step),"=")
			value = str2num(xmlTagContents("hutchTemperature",step))
			if (value != 0)
				noteStr = ReplaceNumberByKey("hutchTemperature",noteStr,value,"=")
			endif
			value = str2num(xmlTagContents("beamBad",step))
			if (value)
				noteStr = ReplaceNumberByKey("beamBad",noteStr,value,"=")
			endif
			str = xmlTagContents("CCDshutter",step)
			if (!stringmatch(str,"out"))
				noteStr = ReplaceStringByKey("CCDshutter",noteStr,"in","=")
			endif
			str = xmlTagContents("monoMode",step)
			noteStr = ReplaceStringByKey("monoMode",noteStr,str,"=")
			if (stringmatch(str,"monochromatic"))
				noteStr = ReplaceNumberByKey("energy",noteStr,str2num(xmlTagContents("energy",step)),"=")
			endif
		endif
		patternKeyVals = xmlTagKeyVals("pattern",step)
		N0 = NumberByKey("Nindexed", patternKeyVals,"=")
		rms = NumberByKey("rms_error", patternKeyVals,"=")
		if (NumberByKey("num", patternKeyVals,"=")!=0 || N0<NindexedMIN)
			//if (strlen(patternKeyVals))
			//	print xmlTagContents("inputImage",step)
			//	print patternKeyVals
			//endif
			continue
		endif
		pattern = xmlTagContents("pattern",step)
		svec = xmlTagContents("astar",pattern)
		sscanf svec,"%g %g %g",as0,as1,as2
		svec = xmlTagContents("bstar",pattern)
		sscanf svec,"%g %g %g",bs0,bs1,bs2
		svec = xmlTagContents("cstar",pattern)
		sscanf svec,"%g %g %g",cs0,cs1,cs2

		if (N>=Nalloc)										// need longer arrays
			Nalloc += 1000
			Redimension/N=(Nalloc) XX,HH,FF,depth
			Redimension/N=(Nalloc) totalSum, sumAboveThreshold, numAboveThreshold
			Redimension/N=(Nalloc) Nindexed
			Redimension/N=(Nalloc) rmsIndexed
			Redimension/N=(Nalloc) imageNames
			Redimension/N=(3,3,Nalloc) gm
		endif
		depth[N] = str2num(xmlTagContents("depth",step))
		XX[N] = str2num(xmlTagContents("Xsample",step)) - Xoff
		Yn = str2num(xmlTagContents("Ysample",step)) - Yoff
		Zn = str2num(xmlTagContents("Zsample",step)) - Zoff
		HH[N] = YZ2H(Yn,Zn)								// H =  Y*sin(angle) + Z*cos(angle))	Yn and Zn are sample system in the file
		FF[N] = YZ2F(Yn,Zn)								// F = -Y*cos(angle) + Z*sin(angle)
		Nindexed[N] = N0
		rmsIndexed[N] = rms
		totalSum[N] = str2num(xmlTagContents("totalSum",step))
		sumAboveThreshold[N] = str2num(xmlTagContents("sumAboveThreshold",step))
		numAboveThreshold[N] = str2num(xmlTagContents("numAboveThreshold",step))
		imageNames[N] = xmlTagContents("inputImage",step)
		gm[0][0][N] = as0;		gm[0][1][N] = bs0;	gm[0][2][N] = cs0	// save measured reciprocal lattice for this point
		gm[1][0][N] = as1;		gm[1][1][N] = bs1;	gm[1][2][N] = cs1
		gm[2][0][N] = as2;		gm[2][1][N] = bs2;	gm[2][2][N] = cs2
		N += 1
	while(1)
	Close refNum
	printf "for reading,  execution time = %s\r",Secs2Time(SecondsInProgressPanel(progressWin),5,0)
	Redimension/N=(N) XX,HH,FF,depth						// trim to actual length
	Redimension/N=(N) totalSum, sumAboveThreshold, numAboveThreshold
	Redimension/N=(N) Nindexed
	Redimension/N=(N) rmsIndexed
	Redimension/N=(3,3,N) gm
	Redimension/N=(N) imageNames
	Make/N=(N)/O/D totalAngles,RF, RH, RX
	SetScale d 0,0,"¡", totalAngles
	// done reading

	Variable Xlo=WaveMin(XX), Xhi=WaveMax(XX)
	Variable Hlo=WaveMin(HH), Hhi=WaveMax(HH)
	Variable Flo=WaveMin(FF), Fhi=WaveMax(FF)
	printf "read in %d steps\r",N
	printf "X range = [%g, %g]\r",Xlo,Xhi
	printf "H range = [%g, %g]\r",Hlo,Hhi
	printf "F range = [%g, %g]\r",Flo,Fhi
	noteStr = ReplaceNumberByKey("Xlo",noteStr,Xlo,"=")
	noteStr = ReplaceNumberByKey("Xhi",noteStr,Xhi,"=")
	noteStr = ReplaceNumberByKey("Hlo",noteStr,Hlo,"=")
	noteStr = ReplaceNumberByKey("Hhi",noteStr,Hhi,"=")
	noteStr = ReplaceNumberByKey("Flo",noteStr,Flo,"=")
	noteStr = ReplaceNumberByKey("Fhi",noteStr,Fhi,"=")
	Duplicate/FREE XX,yztemp
	yztemp = HF2Y(HH[p],FF[p])
	WaveStats/M=1/Q yztemp
	noteStr = ReplaceNumberByKey("Ylo",noteStr,V_min,"=")
	noteStr = ReplaceNumberByKey("Yhi",noteStr,V_max,"=")
	yztemp = HF2Z(HH[p],FF[p])
	WaveStats/M=1/Q yztemp
	noteStr = ReplaceNumberByKey("Zlo",noteStr,V_min,"=")
	noteStr = ReplaceNumberByKey("Zhi",noteStr,V_max,"=")
	WaveClear yztemp

	// shift so that average center of volume is zero
	if (0)
		Variable X0, H0, F0
		WaveStats/M=1/Q XX;		X0 = V_avg;		XX -= X0	// find average center value
		WaveStats/M=1/Q HH;		H0 = V_avg;	HH -= H0
		WaveStats/M=1/Q FF;		F0 = V_avg;		FF -= F0
		printf "re-set origin of data to XHF={%g, %g, %g},   so now {0,0,0} is center\r",X0,H0,F0
		noteStr = ReplaceNumberByKey("X0",noteStr,X0,"=")
		noteStr = ReplaceNumberByKey("H0",noteStr,H0,"=")
		noteStr = ReplaceNumberByKey("F0",noteStr,F0,"=")
		noteStr = ReplaceNumberByKey("Y0",noteStr,HF2Y(H0,F0),"=")
		noteStr = ReplaceNumberByKey("Z0",noteStr,HF2Z(H0,F0),"=")
	endif

	// set rl0 to the reference orientation, how this is done will depend upon 'refType'
	if (refType==2 && (iref>=0 && iref<N))				// a valid reference point
		gmi = gm[p][q][iref]
		if (useSymmetry)									// symmetry reduce this structure
			angle = symReducedRecipLattice(std,id33,gmi,mat3)
		else
			MatrixOp/O mat3 = gmi x Inv(std)				// gmi = mat3 x std
		endif
		MatrixOp/O rl0 = mat3 x std						// reference reciprocal lattice

	elseif (refType==1)									// find central orientation by averaging
		Make/N=3/O/D vecCenter=0
		MatrixOp/O stdInv = Inv(std)
		for (i=0;i<N;i+=1)
			if (mod(i,400) == 0)
				if (ProgressPanelUpdate(progressWin,i/N*100,status="looking for reference, number "+num2istr(i)))	// update progress bar
					break									//   and break out of loop
				endif
			endif
			gmi = gm[p][q][i]
			if (useSymmetry)								// symmetry reduce this structure
				angle = symReducedRecipLattice(std,id33,gmi,mat3)
			else
				MatrixOp/O mat3 = gmi x stdInv				// gmi = mat3 x std
			endif
			angle = axisOfMatrix(mat3,vec3)				// returned angle (degrees)
			vec3 *= angle
			vecCenter += vec3
		endfor
		printf "and finding reference,  execution time = %s\r",Secs2Time(SecondsInProgressPanel(progressWin),5,0)
		vecCenter /= N										// the central rotation vector
		angle = norm(vecCenter)
		rotationMatAboutAxis(vecCenter,angle,mat3)		// mat3 rotates std to the central orientation (the new reference)
		MatrixOp/O rl0 = mat3 x std						// reference reciprocal lattice
		KillWaves/Z vecCenter, stdInv

	elseif (refType==0)									// find central orientation by averaging
		rl0 = std
		mat3 = p==q

	else
		str = "invalid combination refType="+num2str(refType)+",  and iref="+num2str(iref)
		DoAlert 0, str
		print str
		DoWindow/K $progressWin
		return ""
	endif
	sprintf str, "{{%g,%g,%g}{%g,%g,%g}{%g,%g,%g}}",rl0[0][0],rl0[1][0],rl0[2][0],rl0[0][1],rl0[1][1],rl0[2][1],rl0[0][2],rl0[1][2],rl0[2][2]
	noteStr = ReplaceStringByKey("recipRef",noteStr,str,"=")// save reference reicp lattice, as column vectors
	print "using a reference reciprocal lattice of: ",StringFromList(refType,"Standard Orientation;Average Orientation;Orientation from point ")+SelectString(refType==2,"",num2str(iref))
	printf "%.2f\t%.2f\t%.2f\r",rl0[0][0],rl0[0][1],rl0[0][2]
	printf "%.2f\t%.2f\t%.2f\r",rl0[1][0],rl0[1][1],rl0[1][2]
	printf "%.2f\t%.2f\t%.2f\r",rl0[2][0],rl0[2][1],rl0[2][2]
	if (refType==2 && (iref>=0 && iref<N))				// a valid reference point
		printf "from iref=%d, at position (%g, %g, %g)\r",iref,XX[iref],HH[iref],FF[iref]
	endif
	Make/FREE/N=3/D hkl = {0,1,-1}
	MatrixOP/O/FREE hkl = Inv(rl0) x hkl

	Variable maxH = max(-WaveMin(hkl),WaveMax(hkl))
	hkl = round(hkl*24/maxH)
	Variable h=hkl[0], k=hkl[1], l=hkl[2]
	lowestOrderHKL(h,k,l)
	hkl = {h,k,l}
	printf "hkl of surface normal (reference) = %s\r",vec2str(hkl)

	// use current reference orientation (rl0) to compute all Rodriques vectors
	Make/N=(3,3)/O/D rot33
	Make/N=(3,3)/D/FREE rot_Load3dRecipXML=0			// rotates about X-axis (rotates vector not frame) to get from beam-line coordinates to (XHF)
	Wave rotFrame=rot_Load3dRecipXML

	// to go from (XYZ) -> (XHF) is a +45¡ rotation about the X-axis
	rotFrame[0][0] = 1
	rotFrame[1][1] = cos(45*PI/180)						// this matrix rotates direction of vector by +45¡ (this is changed below)
	rotFrame[2][2] = cos(45*PI/180)
	rotFrame[1][2] = -sin(45*PI/180)
	rotFrame[2][1] = sin(45*PI/180)
	MatrixOp/O rotFrame = Inv(rotFrame)					// rotFrame now rotates the description of a vector from (XYZ) frame to (XHF) frame
	Make/N=(N) IndexBackTrack = -1
	Variable m=0, mref=NaN								// use this for rejecting rotations greater than maxAngle
	for (i=0;i<N;i+=1)										// using rl0, find the Rodriques vectors
		if (mod(i,400) == 0)
			if (ProgressPanelUpdate(progressWin,i/N*100,status="setting rotation, number "+num2istr(i)))	// update progress bar
				break										//   and break out of loop
			endif
		endif
		gmi = gm[p][q][i]
		if (useSymmetry)									// symmetry reduce this structure
			angle = symReducedRecipLattice(std,mat3,gmi,rot33)	// angle between go and gm,  gmi = rot33 x rl0 = rot33 x (mat3 x std)
		else
			MatrixOp/O rot33 = gmi x Inv(rl0)				// gmi = rot33 x rl0
		endif
		angle = axisOfMatrix(rot33,vec3)					// returned angle (degrees)
		if (numtype(angle) && i==0)
			str = "Cannot find axis of rotation matrix, Probably the lattice parameters in Igor do not match those in the xml file.\r  Check the 'Xtal' tab"
			DoAlert 0, str
			print str
			DoWindow/K $progressWin
			return ""
		endif
		vec3 *= tan(angle*PI/180/2)						// this is now the Rodriques vector
		MatrixOp/O vec3 = rotFrame x vec3					// rotate from beam line frame to XHF frame (outward surface normal is Z)
		if (angle<=maxAngle)								// keep this one
			depth[m] = depth[i]
			XX[m] = XX[i]
			HH[m] = HH[i]
			FF[m] = FF[i]
			totalSum[m] = totalSum[i]
			sumAboveThreshold[m] = sumAboveThreshold[i]
			numAboveThreshold[m] = numAboveThreshold[i]
			Nindexed[m] = Nindexed[i]
			rmsIndexed[m] = rmsIndexed[i]
			gm[][][m] = gm[p][q][i]
			imageNames[m] = imageNames[i]
			RX[m] = vec3[0]
			RH[m] = vec3[1]
			RF[m] = vec3[2]
			totalAngles[m] = angle
			IndexBackTrack[m] = i
			if (i==iref)
				mref = m
			endif
			m += 1
		endif
	endfor
	printf "total  execution time = %s\r",Secs2Time(SecondsInProgressPanel(progressWin),5,0)
	DoWindow/K $progressWin
	N = m
	Redimension/N=(N) depth,XX,HH,FF,totalSum,sumAboveThreshold,numAboveThreshold,Nindexed,rmsIndexed,RX,RH,RF,totalAngles
	Redimension/N=(N) imageNames
	Redimension/N=(-1,-1,N) gm
	if (refType==2 && (iref>=0 && iref<N))				// a valid reference point
		printf "requested iref=%d,  this corresponds to point %d\r",iref,mref
		noteStr = ReplaceNumberByKey("iref",noteStr,iref,"=")
		noteStr = ReplaceNumberByKey("mref",noteStr,mref,"=")
	endif

	WaveStats/M=1/Q totalAngles
	printf "rotation angles range from:   %.3g¡ at %d  to  %.3g¡ at %d\r",V_min,V_minloc,V_max,V_maxloc
	printf "remaining in Igor data folder  '%s'\r",GetDataFolder(1)
	Wave rotrgb=$makeRGBJZT(RX,RH,RF,V_max)
	KillWaves/Z id33,std,rl0, gmi,vec3,mat3,rot33, rot_Load3dRecipXML

	Note/K XX,noteStr
	Note/K HH,noteStr
	Note/K FF,noteStr
	Note/K depth,noteStr
	Note/K RX,noteStr
	Note/K RH,noteStr
	Note/K RF,noteStr
	Note/K totalSum,noteStr
	Note/K sumAboveThreshold,noteStr
	Note/K numAboveThreshold,noteStr
	Note/K Nindexed,noteStr
	Note/K rmsIndexed,noteStr
	Note/K totalAngles,noteStr
	Note/K imageNames,noteStr
	Note/K gm,ReplaceStringByKey("waveClass", noteStr, "Random3dArraysGm","=")

	Duplicate/O XX XXorig									// transform to voxel coordinate in sample, not sample position
	Duplicate/O HH HHorig
	Duplicate/O FF FForig
	XX = -XXorig
	HH = -HHorig
	FF = -FForig
	Duplicate/O depth depthAddTemp
	depthAddTemp = numtype(depth) ? 0 : depth/sqrt(2)
	FF += depthAddTemp
	HH += depthAddTemp
	KillWaves/Z depthAddTemp
	Duplicate/O HH, YY,ZZ									// for some samples you want to work in XYZ, not XHF
	YY = HF2Y(HH[p],FF[p])
	ZZ = HF2Z(HH[p],FF[p])
	Duplicate/O RH, RY,RZ
	RY = (Rh-RF)/sqrt(2)
	RZ = (Rh+RF)/sqrt(2)

	Variable threeD=0
	WaveStats/Q/M=1 XX
	threeD += (V_max-V_min) > 0.1
	WaveStats/Q/M=1 HH
	threeD += (V_max-V_min) > 0.1
	WaveStats/Q/M=1 FF
	threeD += (V_max-V_min) > 0.1
	threeD = (threeD==3)
	if (threeD)
		Make/N=(N,3)/O xyz
		SetScale d 0,0,WaveUnits(XX,-1), xyz
		xyz[][0] = XX[p]
		xyz[][1] = HH[p]
		xyz[][2] = FF[p]
		MakeGizmocubeCorners(xyz)
		Make/N=(N,4)/O rgba
		Make/N=(N)/O brite
		Variable hi=WaveMax(totalsum)/4
		brite = sqrt(limit(totalsum/hi,0.09,1))
		// brite = sqrt(limit(totalsum/100e3,0.09,1))
		rgba = rotrgb[p][q]/65535 ;  rgba[][3] = brite[p]
		Note/K xyz,ReplaceStringByKey("waveClass",noteStr,"Random3dArraysXYZ","=")

		String note2= ReplaceStringByKey("waveClass",noteStr,"Random3dArraysRGBA","=")
		note2= ReplaceStringByKey("source",noteStr,"RxRyRz","=")
		String title=StringByKey("title",note(RX),"=")
		title = "RxRyRz"+SelectString(strlen(title),"",", ")+title
		note2= ReplaceStringByKey("title",noteStr,title,"=")
		Note/K rgba,note2
	endif
	return noteStr
End
//Static Function axisOfMatrix(rot,axis)
//	// returns total rotation angle, and sets axis to the axis of the total rotation
//	Wave rot										// rotation matrix
//	Wave axis										// axis of the rotation (angle is returned)
//
//	Variable sumd = notRotationMatrix(rot)	// accept positive sumd that are less than 1e-4
//	if (sumd<0 || sumd>1e-4)
//		DoAlert 0, "'"+NameOfWave(rot)+"' is not a rotation matrix in axisOfMatrix(), did you set the correct lattice?"
//		axis = NaN
//		return NaN
//	elseif (0<sumd)								// close enough to a roation matix, but tidy it up first
//		Make/N=(3,3)/O/D axisMat__
//		axisMat__ = rot
//		if (SquareUpMatrix(axisMat__))
//			DoAlert 0, "cannot square up '"+NameOfWave(rot)+"' in axisOfMatrix(), did you set the correct lattice?"
//			axis = NaN
//			return NaN
//		endif
//		Wave rr = axisMat__
//	else
//		Wave rr = rot
//	endif
//
//	Variable cosine = (MatrixTrace(rr)-1)/2	// trace = 1 + 2*cos(theta)
//	cosine = limit(cosine,-1,1)
//	if (cosine<= -1)								// special for 180¡ rotation,
//		axis[0] = sqrt((rr[0][0]+1)/2)
//		axis[1] = sqrt((rr[1][1]+1)/2)
//		axis[2] = sqrt((rr[2][2]+1)/2)		// always assume z positive
//		axis[0] = (rr[0][2]+rr[2][0])<0 ? -axis[0] : axis[0]
//		axis[1] = (rr[1][2]+rr[2][1])<0 ? -axis[1] : axis[1]
//	else											// rotaion < 180¡, usual formula works
//		axis[0] = rr[2][1] - rr[1][2]
//		axis[1] = rr[0][2] - rr[2][0]
//		axis[2] = rr[1][0] - rr[0][1]
//		axis /= 2
//	endif
//	normalize(axis)
//	KillWaves /Z axisMat__
//	return acos(cosine)*180/PI					// rotation angle in degrees
//End
//
Function/T makeRGBJZT(RX,RH,RF,maxAngle)
	Wave RX,RH,RF
	Variable maxAngle								// angle where colors saturate (degree)
	Variable i, N=numpnts(RX)
	Make/N=(N,3)/U/W/O rotRGB=0

	Make/N=3/O/D vec3
	Variable cx,cy,cz, angle
	Variable r,g,b
	for (i=0;i<N;i+=1)
		vec3[0] = RX[i]
		vec3[1] = RH[i]
		vec3[2] = RF[i]
		angle = 2*atan(normalize(vec3)) * 180/PI
		vec3 *= angle/maxAngle					// length is now rotation normalized to maxAngle
		cx = limit(vec3[0],-1,1)*65535
		cy = limit(vec3[1],-1,1)*65535
		cz = limit(vec3[2],-1,1)*65535
		r=0;	g=0;	b=0
		if (cx>0)
			r += cx
		else
			g += abs(cx)/2
			b += abs(cx)/2
		endif
		if (cy>0)
			g += cy
		else
			r += abs(cy)/2
			b += abs(cy)/2
		endif
		if (cz>0)
			b += cz
		else
			r += abs(cz)/2
			g += abs(cz)/2
		endif
		r =min(r,65535)
		g =min(g,65535)
		b =min(b,65535)
		if (r>65535 || g>65535 || b>65535)
			print "number too big at i=",i
		endif
		rotRGB[i][0] = r
		rotRGB[i][1] = g
		rotRGB[i][2] = b
	endfor

	String noteStr = ReplaceStringByKey("waveClass",note(RX),"Random3dArraysRGB","=")
	noteStr = ReplaceNumberByKey("maxAngle",noteStr,maxAngle,"=")
	Note/K rotRGB, noteStr
	KillWaves/Z vec3
	return GetWavesDataFolder(rotRGB,2)
End

//Function MakeGizmocubeCorners(xyz)			// this routine improved and moved to GizmoZoomTranslate.ipf
//	Wave xyz							// xyz[N][3]
//	Variable N=DimSize(xyz,0)
//
//	Make/N=(N)/O/D cubeCorners_temp_
//	cubeCorners_temp_ = xyz[p][0]
//	WaveStats/Q cubeCorners_temp_
//	Variable Xlo=V_min, Xhi=V_max
//	cubeCorners_temp_ = xyz[p][1]
//	WaveStats/Q cubeCorners_temp_
//	Variable Ylo=V_min, Yhi=V_max
//	cubeCorners_temp_ = xyz[p][2]
//	WaveStats/Q cubeCorners_temp_
//	Variable Zlo=V_min, Zhi=V_max
//	KillWaves/Z cubeCorners_temp_
//
//	Variable dX=Xhi-Xlo, dY=Yhi-Ylo, dZ=Zhi-Zlo
//	Variable d = max(max(dX,dY),dZ)/2, mid
//
//	Variable printIt = strlen(GetRTStackInfo(2))<=0
//	if (printIt)
//		printf "X = [%g, %g], 	Æ=%g\r",Xlo,Xhi,dX
//		printf "Y = [%g, %g], 		Æ=%g\r",Ylo,Yhi,dY
//		printf "Z = [%g, %g], 	Æ=%g\r",Zlo,Zhi,dZ
//	endif
//	mid = (Xlo+Xhi)/2;		Xlo = mid-d;	Xhi = mid+d
//	mid = (Ylo+Yhi)/2;		Ylo = mid-d;	Yhi = mid+d
//	mid = (Zlo+zhi)/2;		Zlo = mid-d;	Zhi = mid+d
//	if (printIt)
//		print " "
//		printf "X = [%g, %g], 	Æ=%g\r",Xlo,Xhi,Xhi-Xlo
//		printf "Y = [%g, %g], 		Æ=%g\r",Ylo,Yhi,Yhi-Ylo
//		printf "Z = [%g, %g], 	Æ=%g\r",Zlo,Zhi,Zhi-Zlo
//	endif
//
//	Make/N=(8,3)/O cubeCorners=NaN
//	cubeCorners[0][0] = Xlo;	cubeCorners[0][1] = Ylo; 	cubeCorners[0][2] = Zlo
//	cubeCorners[1][0] = Xhi;	cubeCorners[1][1] = Ylo; 	cubeCorners[1][2] = Zlo
//	cubeCorners[2][0] = Xlo;	cubeCorners[2][1] = Yhi; 	cubeCorners[2][2] = Zlo
//	cubeCorners[3][0] = Xhi;	cubeCorners[3][1] = Yhi; 	cubeCorners[3][2] = Zlo
//	cubeCorners[4][0] = Xlo;	cubeCorners[4][1] = Ylo; 	cubeCorners[4][2] = Zhi
//	cubeCorners[5][0] = Xhi;	cubeCorners[5][1] = Ylo; 	cubeCorners[5][2] = Zhi
//	cubeCorners[6][0] = Xlo;	cubeCorners[6][1] = Yhi; 	cubeCorners[6][2] = Zhi
//	cubeCorners[7][0] = Xhi;	cubeCorners[7][1] = Yhi; 	cubeCorners[7][2] = Zhi
//End

// ***********************************  End of Read XML ************************************
// *****************************************************************************************



// ****************************  Start of XML Reading Utilities ****************************
// *****************************************************************************************

Function LoadPixelsFromXML(m)	// for the point m in imageNames[m], load in the fitted pixels from the xml file
	Variable m						// point number in list

	Wave/T imageNames=imageNames
	if (!WaveExists(imageNames))
		return 1
	endif
	Variable Ni=DimSize(imageNames,0),N
	if (!(m>=0 && m<Ni))
		Prompt m, "index into list [0,"+num2str(Ni-1)+"]"
		DoPrompt "index",m
		if (V_flag)
			return 1
		endif
	endif
	if (!(m>=0 && m<Ni))
		return 1
	endif
	String fileName = StringByKey("xmlFileFull",note(imageNames),"=")	// xml file name from wave note
	GetFileFolderInfo/Q/Z=1 fileName
	if (!V_isFile && strlen(GetRTStackInfo(2))==0)
		Variable refNum
		Open/D/M="xml file"/R/T=".xml" refNum
		fileName = S_fileName
	endif
	if (strlen(fileName)<1)
		return 1
	endif

	String wName="pxy_"+num2istr(m)
	Variable newWave = (exists(wName)!=1)
	Make/N=(N,2)/O $wName
	Wave pxy = $wName

	Make/N=(1,2)/O getPixelsFromXML_tempXY						// temp wave to hold XY values
	String iname = ReplaceString("/",imageNames[m], ":")

	N = getPixelXY(fileName,iname,pxy)
	if (N<1)
		if (newWave)
			KillWaves/Z pxy
		endif
		return 1
	endif
	return 0
End


Function getPixelXY(FullFileName,imageName,waveXY)		// returns number of points read into waveXY
	String FullFileName
	String imageName
	Wave waveXY
	imageName = ParseFilePath(3, imageName, ":", 0, 0)

	if (!WaveExists(waveXY))
		Abort "waveXY does not exits"
	endif

	String pathSep
	if (stringmatch(IgorInfo(2),"Macintosh"))
		pathSep = ":"
	elseif (stringmatch(IgorInfo(2),"Windows"))
		pathSep = "\\"
	else
		Abort "This is neither Mac nor Win, what is this computer?"
	endif
	Variable f
	Open/R/P=home/Z=1 f as FullFileName
	if (strlen(S_fileName)<1 || V_flag)
		return 0
	endif
	Close f
	FullFileName = S_fileName

	String names = MakeStepIndexForXML(FullFileName)
	Wave indexPos = $StringFromLIst(0,names)
	Wave/T indexImage = $StringFromLIst(1,names)
	if (!WaveExists(indexPos) || !WaveExists(indexImage))
		return 0
	endif

	// read step containing imageName from file
	Variable Nindex=DimSize(indexImage,0), index
	for (index=0;index<Nindex;index+=1)				// loop over indexImage to find point number for imageName
		if (stringmatch(indexImage[index],imageName))
			break
		endif
	endfor
	if (index>=Nindex)									// image name not found in indexImage[]
		return 0
	endif
	Variable i0=indexPos[index], i1=indexPos[index+1]
	if (numtype(i0+i1) || i1<=i0)						// invalid range in file [i0,i1]
		return 0
	endif

	String step = PadString("",i1-i0+1,0x20)			// a buffer to hold step (of length i1-i0+1 bytes)
	Open/R/Z=1 f as FullFileName						// start reading here
	FStatus f
	Variable fileLen = V_logEOF							// length of file in bytes
	if (i1>=V_logEOF)									// file too short
		Close f
		return 0
	endif
	FSetPos f, i0
	FBinRead f, step										// initial read
	Close f

	String svec = xmlTagContents("Xpixel",step)		// got <step>, extract data
	svec = ReplaceString(" ",svec,";")
	Variable N=ItemsInList(svec)
	if (N>=1)
		Redimension/N=(N,2) waveXY
		waveXY[][0] = str2num(StringFromList(p,svec))
		svec = xmlTagContents("Ypixel",step)
		svec = ReplaceString(" ",svec,";")
		waveXY[][1] = str2num(StringFromList(p,svec))
	endif
	return N
End
//	Function test_getPixelXY()
//		Variable f
//		// String xmlFileFull = "Macintosh HD:Users:tischler:Desktop:SampleN_2D2_resulttest.xml"
//		String xmlFileFull = "Macintosh HD:Users:tischler:Desktop:SampleN_2D2_result.xml"
//	
//		Open/R/Z=1 f as xmlFileFull
//		if (V_flag)
//			Open/R/T=".xml" f
//			Close f
//			xmlFileFull = S_filename
//			print "opening ",xmlFileFull
//		endif
//		String imageName = "SampleN_2D2_571"
//		Make/N=(2)/O waveXYtest
//		Variable N = getPixelXY(xmlFileFull,imageName,waveXYtest)		// returns number of points read into waveXY
//		print "N=",N
//		print "created wave  'waveXYtest', check it out"
//	End



Function LatticeParametersFromXML(FullFileName,xtal)
	String FullFileName
	STRUCT crystalStructure &xtal

	Variable f
	Open/R/Z=1 f as FullFileName
	if (V_flag)
		return 1
	endif
	Variable filePos = 0									// current file pos
	Variable fileLen										// length of file in bytes
	Variable size,sizeMax=20*1024					// size of a typical read
	Variable bytesRead									// keeps track of how many bytes read in (should not exceed fileLen)
	String buf = ""										// a buffer for working space (of length size bytes)

	Variable i0=-1, i1=-1
	do
		FStatus f
		filePos = V_filePos								// current file pos
		fileLen = V_logEOF								// length of file in bytes
		size = min(sizeMax,fileLen-filePos)
		if (size<2)
			break
		endif
		buf = PadString("",size,0x20)					// a buffer for working space (of length size bytes)
		FBinRead f, buf									// initial read

		i0 = strsearch(buf,"<xtl>",0)					// search for start of <xtl>
		if (i0>=0)										// found start of <xtl>, so look for end
			i1 = strsearch(buf,"</xtl>",0)				// end of <xtl>
			if (i1>0)
				break
			endif
			FSetPos f, filePos+i0						// set position to start of <xtl>
			FStatus f
			filePos = V_filePos							// current file pos
			size = min(sizeMax,V_logEOF-filePos)
			if (size<2)
				break
			endif
			buf = PadString("",size,0x20)
			FBinRead f, buf
			i1 = strsearch(buf,"</xtl>",0)				// do I also have end?
			break
		endif
	while(1)
	Close f
	if (i0<0 || i1<i0)
		return 1
	endif
	buf = xmlTagContents("xtl",buf[i0,i1+5])
	buf = ReplaceString("\t",buf," ")

	Variable SpaceGroup = str2num(xmlTagContents("SpaceGroup",buf))
	if (!(SpaceGroup>=1 && SpaceGroup<=230))
		return 1
	endif
	xtal.SpaceGroup = SpaceGroup

	String str = xmlTagContents("latticeParameters",buf)
	Variable a,b,c,alpha,bet,gam
	sscanf str, "%g %g %g %g %g %g",a,b,c,alpha,bet,gam
	if (V_flag!=6 || numtype(a+b+c+alpha+bet+gam))
		return 1
	endif
	xtal.a = a
	xtal.b = b
	xtal.c = c
	xtal.alpha = alpha
	xtal.beta = bet
	xtal.gam = gam
	xtal.Unconventional00 = NaN
	xtal.N = 0							// init to no atoms

	// read the atoms
	String list,atomLabel
	Variable n,Zatom,occ,DebyeT, Ntotal=0
	Variable ax,ay,az
	i0 = 0
	do
		i0 = strsearch(buf,"<atom ",i0)
		if (i0<0)
			break
		endif

		str = xmlTagContents("atom",buf[i0,Inf])
		sscanf str,"%g %g %g",ax,ay,az
		if (V_flag!=3 || numtype(ax+ay+az))
			break
		endif

		list = xmlTagKeyVals("atom",buf[i0,Inf])
		n = NumberByKey("n",list,"=")
		Zatom = NumberByKey("Z",list,"=")
		atomLabel = StringByKey("label",list,"=")
		atomLabel = atomLabel[0,59]
		if (numtype(n) || n>=20 || numtype(Zatom) || Zatom<0 || Zatom>100 || strlen(atomLabel)<1)
			break
		endif
		Ntotal = max(Ntotal,n)
		occ = NumberByKey("occ",list,"=")
		occ = numtype(occ) ? 1 : occ
		DebyeT = NumberByKey("DebyeT",list,"=")
		DebyeT = numtype(DebyeT) ? 0 : DebyeT
		xtal.N = Ntotal
		xtal.atom[n-1].x = ax
		xtal.atom[n-1].y = ay
		xtal.atom[n-1].z = az
		xtal.atom[n-1].name = atomLabel
		xtal.atom[n-1].Zatom = Zatom
		xtal.atom[n-1].occ = occ
		xtal.atom[n-1].DebyeT = DebyeT

		i0 = strsearch(buf,"</atom>",i0)
	while(n>0 && n<20)

	LatticeSym#ForceLatticeToStructure(xtal)
	return 0
End
//
//	Function test_LatticeParametersFromXML()
//		Variable f
//		// String xmlFileFull = "Macintosh HD:Users:tischler:Desktop:SampleN_2D2_resulttest.xml"
//		String xmlFileFull = "Macintosh HD:Users:tischler:Desktop:SampleN_2D2_result.xml"
//		Open/R/Z=1 f as xmlFileFull
//		if (V_flag)
//			Open/R/T=".xml" f
//			Close f
//			xmlFileFull = S_filename
//			print "opening ",xmlFileFull
//		endif
//	
//		STRUCT crystalStructure xtal
//		xtal.SpaceGroup = 0
//		LatticeParametersFromXML(xmlFileFull,xtal)
//		print_crystalStructure(xtal)					// prints out the value in a crystalStructure structure
//	End



//Static Function/T xmlTagContents(key,buf)
//	String key
//	String buf
//
//	Variable j0,jend,jmid
//	j0 = strsearch(buf,"<"+key+">",0)
//	if (j0<0)
//		j0 = strsearch(buf,"<"+key+" ",0)
//		if (j0<0)
//			j0 = strsearch(buf,"<"+key+" \t",0)
//		endif
//	endif
//	if (j0<0)
//		return ""
//	endif
//
//	jmid = strsearch(buf,">",j0)
//	if (jmid<0)
//		return ""
//	endif
//	jmid += 1
//
//	jend = strsearch(buf,"</"+key+">",0)
//	if (jend<0)
//		return ""
//	endif
//	jend -= 1
//
//	String contents = buf[jmid,jend]
//	return contents
//End


Static Function/T xmlTagKeyVals(key,buf)
	String key
	String buf

	Variable j0,jend,jmid
	j0 = strsearch(buf,"<"+key+" ",0)
	if (j0<0)
		j0 = strsearch(buf,"<"+key+" \t",0)
	endif
	if (j0<0)
		return ""
	endif
	j0 += 2+strlen(key)
	jmid = strsearch(buf,">",j0)-1
	if (jmid<=j0)
		return ""
	endif

	String str
	str = ReplaceString("\t", buf[j0,jmid]," ")
	str = ReplaceString("  ", str," ")
	str = ReplaceString(" ", str,";")
	str = ReplaceString("\"", str,"")
	return str
End

// *********************************  End of XML Utilities *********************************
// *****************************************************************************************



// *****************************************************************************************
// *************************  Start Add Hand Indexed Point to XML **************************

Static Constant STRUCTURE_N_DETECTORS_MAX=3
Static Constant STRUCTURE_N_PEAKS_MAX=200
Static Constant STRUCTURE_N_PATTERNS_MAX=10

Static Structure step3DinfoStruct		// structure definition for contents of one 3D step in an XML
	char title[100] 					// title
	char sampleName[100] 			// sample ame
	char userName[100] 				// user name
	char beamline[100] 				// name of beam line
	double scanNum						// scan number
	char date[50]						// date & time e.g. 2013-06-16T19:59:54-06:00
	int16 beamBad						// 0 or 1
	char monoMode[50]					// something like "white slitted"
	double Xsample						// sample X position
	double Ysample	
	double Zsample	
	double depth
	double energy						// energy, unit="keV"
	double hutchTemperature		// hutch temperature (C)
	double sampleDistance			// probably 0

	int16 Ndetectors					// number of detector structures defined (<=STRUCTURE_N_DETECTORS_MAX)
	Struct detector3DinfoStruct d[STRUCTURE_N_DETECTORS_MAX]
	Struct indexed3DinfoStruct index
EndStructure
	//
	Static Structure detector3DinfoStruct		// structure definition for contents of one detector in a step3DinfoStruct
		char inputImage[250]			// name of image file
		char	detectorID[100]			// detector ID
		double exposure					// exposure time unit="sec"
		int16 Nx,Ny							// number of pixels in detector
		double totalSum					// sum of all pixels in image
		double sumAboveThreshold		// sum all pixels that pass threshold
		double numAboveThreshold		// number of pixels that pass threshold
		char geoFile[250]					// name of geo file
		Struct roi3DinfoStruct roi	// roi saved from this detector
		Struct peaks3DinfoStruct peaks
	EndStructure
		//
		Static Structure roi3DinfoStruct			// structure definition for contents of one roi on a detector
			int16 startx, endx, groupx
			int16 starty, endy, groupy
		EndStructure
		//
		Static Structure peaks3DinfoStruct		// structure definition for contents of each peak on a detector
			char peakProgram[100]			// name of peak program
			double minwidth
			double threshold
			double thresholdRatio
			double maxRfactor
			double maxwidth
			double maxCentToFit
			double boxsize
			double max_number
			double min_separation
			int16 smooth						// 1 or 0
			char peakShape[50]				// probably "Gaussian" or Lorentzian
			double executionTime
			int16 Npeaks						// number of peaks in this structure
			double Xpixel[STRUCTURE_N_PEAKS_MAX]
			double Ypixel[STRUCTURE_N_PEAKS_MAX]
			double Intens[STRUCTURE_N_PEAKS_MAX]
			double Integral[STRUCTURE_N_PEAKS_MAX]
			double hwhmX[STRUCTURE_N_PEAKS_MAX]
			double hwhmY[STRUCTURE_N_PEAKS_MAX]
			double tilt[STRUCTURE_N_PEAKS_MAX]
			double chisq[STRUCTURE_N_PEAKS_MAX]
			double Qx[STRUCTURE_N_PEAKS_MAX]
			double Qy[STRUCTURE_N_PEAKS_MAX]
			double Qz[STRUCTURE_N_PEAKS_MAX]
		EndStructure
	//
	Static Structure indexed3DinfoStruct		// defines result of indexing
		char indexProgram[100]
		int16 Nindexed
		int16 Npeaks
		double keVmaxCalc
		double keVmaxTest
		double angleTolerance
		double cone	
		char hklPrefer[50]
		double executionTime
	
		char xtl_structureDesc[250]
		int16 xtl_SpaceGroup
		double xtl_a									// units are nm
		double xtl_b
		double xtl_c
		double xtl_alpha
		double xtl_beta
		double xtl_gamma
	
		int16 Npatterns
		Struct pattern3DinfoStruct pattern[STRUCTURE_N_PATTERNS_MAX]
	EndStructure
		//
		Static Structure pattern3DinfoStruct		// defines one pattern in the indexing
			int16 num
			double rms_error
			double goodness
			int16 Nindexed
			double astar[3]			// unit="1/nm"
			double bstar[3]
			double cstar[3]
			double h[STRUCTURE_N_PEAKS_MAX]
			double k[STRUCTURE_N_PEAKS_MAX]
			double l[STRUCTURE_N_PEAKS_MAX]
		EndStructure


Function/T AppendCurrenIndexResult2XML(wIndex,fileName,[overwrite])
	Wave wIndex
	String fileName
	Variable overwrite						// 1=overwrite existing (or new), 0=append to existing
	overwrite = ParamIsDefault(overwrite) ? NaN : overwrite

	if (!WaveExists(wIndex))
		String list = reverseList(WaveListClass("IndexedPeakList*","*",""))
		Variable i = ItemsInLIst(list)
		if (i<1)
			DoAlert 0, "No waves of this type in current folder,\rNothing Done"
			return ""
		elseif (i==1)
			Wave wIndex = $StringFromList(0,list)
		else
			String wName = StringFromList(i-1,list)
			Prompt wName, "Wave with list of indexed peaks",popup,list
			Prompt overwrite,"Overwrite existing or New XML", popup, "Overwrite/Create new;Append to Existing"
			overwrite = !(!overwrite) && numtype(overwrite)==0 ? 1 : 2
			DoPrompt "FullPeakIndexed",wName, overwrite
			if (V_flag)
				return ""
			endif			
			overWrite = overWrite==1
			Wave wIndex = $wName
		endif
		if (!WaveExists(wIndex))
			DoAlert 0, "No waves of name "+wName
			return ""
		endif
	endif
	if (numtype(overwrite))				// have not yet set overwrite
		Prompt overwrite,"Overwrite existing or New XML", popup, "Overwrite/Create new;Append to Existing"
		overwrite = !(!overwrite) && numtype(overwrite)==0 ? 1 : 2
		DoPrompt "OverWrite XML?",overwrite
		if (V_flag)
			return ""
		endif			
		overWrite = overWrite==1
	endif
	if (numtype(overwrite))
		DoAlert 0, "overwrite = "+num2str(overwrite)
		return ""
	endif

	PathInfo home
	String path = SelectString(V_flag,"","home")		// default to home if available
	String last3DxmlFileName = StrVarOrDefault("root:Packages:micro:xml3D:last3DxmlFileName","")

	String fileFilter = "XML Files (*.xml):.xml;All Files:.*;"
	Variable f
	if (overwrite)													// these two Open commands only get file name, do not open a file
		Open/D/C="R*ch"/F=fileFilter/M="New XML, or Existing to overwrite"/P=$path/Z=2 f as last3DxmlFileName
	else
		Open/D/A/F=fileFilter/M="Existing XML for appending"/P=$path/Z=2 f as last3DxmlFileName
	endif

	STRUCT step3DinfoStruct step
	if (FillStepStructure(step,wIndex))
		print "ERROR -- FillStepStructure(), could not fill the 'step' structure"
		return ""
	endif
	last3DxmlFileName = AppendTo3dXML(step,path,S_filename,overwrite=overwrite)	// AppendTo3dXML(step,path,filename,[overwrite])
	if (StringMatch(path,"home"))							// home exists, shorten last3DxmlFileName if it is in home
		PathInfo home
		if (strsearch(last3DxmlFileName,S_path,0,2)==0)
			last3DxmlFileName = last3DxmlFileName[strlen(S_path),Inf]
		endif
	endif

	if (strlen(last3DxmlFileName) && DataFolderExists("root:Packages:micro"))
		if (!DataFolderExists("root:Packages:micro:xml3D"))
			NewDataFolder/O root:Packages:micro:xml3D
		endif
		String/G root:Packages:micro:xml3D:last3DxmlFileName = last3DxmlFileName
		printf "%s Current indexed point to \"%s\"\r",SelectString(overwrite,"Appended","New file, wrote"),last3DxmlFileName
	else
		print "Nothing done"
	endif
	return last3DxmlFileName
End
//
Static Function FillStepStructure(step,FullPeakIndex)
	STRUCT step3DinfoStruct &step			// this structure gets filled here
	Wave FullPeakIndex

	STRUCT microGeometry g						// this will be needed later to comput Q's from pixels
	FillGeometryStructDefault(g)				//fill the geometry structure with current values

	String wnote = note(FullPeakIndex)
	Wave FullPeakList = $StringByKey("peakListWave",wnote,"=")
	if (!WaveExists(FullPeakIndex) || !WaveExists(FullPeakList))
		return 1
	endif

	// load the index part
	Variable Npatterns=NumberByKey("NpatternsFound",wnote,"=")
//	Npatterns = min(1,Npatterns)			// This routine is limited to only 1 pattern, not yet knows how to write multiple
	step.index.Npatterns = Npatterns
	step.index.indexProgram = StrVarOrDefault("root:Packages:micro:indexingExecutableMac","euler")
	step.index.Nindexed = NumberByKey("Nindexed",wnote,"=")
	step.index.Npeaks = NumberByKey("NiData",wnote,"=")
	step.index.keVmaxCalc = NumberByKey("keVmaxCalc",wnote,"=")
	step.index.keVmaxTest = NumberByKey("keVmaxTest",wnote,"=")
	step.index.angleTolerance = NumberByKey("angleTolerance",wnote,"=")
	step.index.cone = NumberByKey("cone",wnote,"=")
	String str = StringByKey("hklPrefer",wnote,"=")
	step.index.hklPrefer = str[1,strlen(str)-2]
	step.index.executionTime = NumberByKey("executionTime",wnote,"=")
	str = StringByKey("latticeParameters",wnote,"=")
	str = str[1,strlen(str)-2]
	step.index.xtl_a = str2num(StringFromList(0,str,","))
	step.index.xtl_b = str2num(StringFromList(1,str,","))
	step.index.xtl_c = str2num(StringFromList(2,str,","))
	step.index.xtl_alpha = str2num(StringFromList(3,str,","))
	step.index.xtl_beta = str2num(StringFromList(4,str,","))
	step.index.xtl_gamma = str2num(StringFromList(5,str,","))
	step.index.xtl_SpaceGroup = NumberByKey("SpaceGroup",wnote,"=")
	step.index.xtl_structureDesc = StringByKey("structureDesc",wnote,"=")

	Make/N=(DimSize(FullPeakIndex,0))/FREE testCount
	Variable i, m, Nindexed
	for (m=0;m<Npatterns;m+=1)
		testCount = FullPeakIndex[p][0][m]		// check first column of FullPeakIndex to count valid points in layer m
		WaveStats/Q/M=1 testCount
		Nindexed = V_npnts							// number of points indexed in pattern m
		step.index.pattern[m].Nindexed = Nindexed
		step.index.pattern[m].num = m
		step.index.pattern[m].rms_error = NumberByKey("rms_error"+num2istr(m),wnote,"=")
		step.index.pattern[m].goodness = NumberByKey("goodness"+num2istr(m),wnote,"=")

		str = StringByKey("recip_lattice"+num2istr(m),wnote,"=")
		Wave recip = str2recip(str)
		step.index.pattern[m].astar[0] = recip[0][0]
		step.index.pattern[m].astar[1] = recip[1][0]
		step.index.pattern[m].astar[2] = recip[2][0]
		step.index.pattern[m].bstar[0] = recip[0][1]
		step.index.pattern[m].bstar[1] = recip[1][1]
		step.index.pattern[m].bstar[2] = recip[2][1]
		step.index.pattern[m].cstar[0] = recip[0][2]
		step.index.pattern[m].cstar[1] = recip[1][2]
		step.index.pattern[m].cstar[2] = recip[2][2]
		for (i=0;i<Nindexed;i+=1)
			step.index.pattern[m].h[i] = FullPeakIndex[i][3][m]
			step.index.pattern[m].k[i] = FullPeakIndex[i][4][m]
			step.index.pattern[m].l[i] = FullPeakIndex[i][5][m]
		endfor
	endfor

	// done with index part, load the peak searching part
	wnote = note(FullPeakList)
	Variable Npeaks = DimSize(FullPeakList,0)
	step.title = StringByKey("title",wnote,"=")
	step.sampleName = StringByKey("sampleName",wnote,"=")
	step.userName = StringByKey("userName",wnote,"=")
	step.beamline = StringByKey("beamline",wnote,"=")
	step.scanNum = NumberByKey("scanNum",wnote,"=")
	step.date = StringByKey("file_time",wnote,"=")
	step.beamBad = NumberByKey("BeamBad",wnote,"=")
	step.beamBad = (step.beamBad==-1) ? 0 : step.beamBad
	step.Xsample = NumberByKey("X1",wnote,"=")
	step.Ysample = NumberByKey("Y1",wnote,"=")
	step.Zsample = NumberByKey("Z1",wnote,"=")
	Variable depth = NumberByKey("depth",wnote,"=")
	step.depth = depth
	depth = numtype(depth) ? 0 : depth			// use this later when calculating Qxyz
	str = StringByKey("monoMode",wnote,"=")
	step.monoMode = str
	if (strsearch(str,"white",0,2)==0)
		step.energy = NaN								// no monochromator energy when in white beam mode
	else
		step.energy = NumberByKey("keV",wnote,"=")
	endif
	step.hutchTemperature = NumberByKey("HutchTemperature",wnote,"=")
	step.sampleDistance = NumberByKey("sampleDistance",wnote,"=")

	Make/N=3/D/FREE qhat							// used to calculate the Qxyz
////////////////////
//			step.Ndetectors = ???
step.Ndetectors = 1			// Code is so far limited to only 1 detector
Variable id = 0				// detector number
////////////////////
////////////////////
	step.d[id].inputImage = StringByKey("file_name",wnote,"=")
	step.d[id].detectorID = StringByKey("detectorID",wnote,"=")
	step.d[id].exposure = NumberByKey("exposure",wnote,"=")
	step.d[id].Nx = NumberByKey("xDimDet",wnote,"=")
	step.d[id].Ny = NumberByKey("yDimDet",wnote,"=")
	step.d[id].totalSum = NumberByKey("totalIntensity",wnote,"=")
	step.d[id].sumAboveThreshold = NumberByKey("totalPeakIntensity",wnote,"=")
	step.d[id].numAboveThreshold = NumberByKey("NumPixelsAboveThreshold",wnote,"=")
	step.d[id].geoFile = StringByKey("geoFile",wnote,"=")
	step.d[id].roi.startx = NumberByKey("startx",wnote,"=")
	step.d[id].roi.endx = NumberByKey("endx",wnote,"=")
	step.d[id].roi.groupx = NumberByKey("groupx",wnote,"=")
	step.d[id].roi.starty = NumberByKey("starty",wnote,"=")
	step.d[id].roi.endy = NumberByKey("endy",wnote,"=")
	step.d[id].roi.groupy = NumberByKey("groupy",wnote,"=")
	str = StringByKey("peakProgram",wnote,"=")
	if (strlen(str)<1)
		str = StrVarOrDefault("root:Packages:micro:PeakFit:fitPeakFuncLast","")
	endif
	step.d[id].peaks.peakProgram = str
	step.d[id].peaks.minwidth = NumberByKey("minPeakWidth",wnote,"=")
	step.d[id].peaks.threshold = NumberByKey("threshAboveAvg",wnote,"=")
	step.d[id].peaks.thresholdRatio = NumberByKey("thresholdRatio",wnote,"=")
	step.d[id].peaks.maxRfactor = NumberByKey("maxRfactor",wnote,"=")
	step.d[id].peaks.maxwidth = NumberByKey("maxPeakWidth",wnote,"=")
	step.d[id].peaks.maxCentToFit = NumberByKey("maxCentToFit",wnote,"=")
	step.d[id].peaks.boxsize = NumberByKey("boxsize",wnote,"=")
	step.d[id].peaks.max_number = NumberByKey("max_number",wnote,"=")
	step.d[id].peaks.min_separation = NumberByKey("minSpotSeparation",wnote,"=")
	step.d[id].peaks.smooth = NumberByKey("smooth",wnote,"=")
	step.d[id].peaks.smooth = (step.d[id].peaks.smooth)==-1 ? 0 : step.d[id].peaks.smooth
	step.d[id].peaks.peakShape = StringByKey("peakShape",wnote,"=")
	step.d[id].peaks.Npeaks = Npeaks
	step.d[id].peaks.executionTime = NumberByKey("executionTime",wnote,"=")
	for (i=0;i<Npeaks;i+=1)
		step.d[id].peaks.Xpixel[i] = FullPeakList[i][0]
		step.d[id].peaks.Ypixel[i] = FullPeakList[i][1]
		step.d[id].peaks.Intens[i] = FullPeakList[i][10]
		step.d[id].peaks.Integral[i] = FullPeakList[i][10]
		step.d[id].peaks.hwhmX[i] = FullPeakList[i][4]/2
		step.d[id].peaks.hwhmY[i] = FullPeakList[i][5]/2
		step.d[id].peaks.tilt[i] = FullPeakList[i][8]
		step.d[id].peaks.chisq[i] = FullPeakList[i][9]
		pixel2q(g.d[id],FullPeakList[i][0],FullPeakList[i][1],qhat,depth=depth)
		step.d[id].peaks.Qx[i] = qhat[0]
		step.d[id].peaks.Qy[i] = qhat[1]
		step.d[id].peaks.Qz[i] = qhat[2]
	endfor
	return 0
End
//
Static Function/T AppendTo3dXML(step,path,filename,[overwrite])	// append data (or optionally start fresh) a 3d XML file
	STRUCT step3DinfoStruct &step
	String path
	String filename			// path+filename is where the new file should go
	Variable overwrite		// defaults to FALSE
	overwrite = ParamIsDefault(overwrite) ? 0 : overwrite
	overwrite = numtype(overwrite) ? 0 : !(!overwrite)

	if (strlen(filename)<1)
		return ""
	endif
	PathInfo $path
	path = SelectString(V_flag,"",path)		// if path does not exist, set to ""

	// now sent out to the file
	GetFileFolderInfo/P=$path/Q/Z filename
	if (V_isFolder)
		printf "ERROR -- in AppendTo3dXML(), the requested file '%s' is a Folder\r",filename
		return ""
	endif
	overwrite = overwrite || !V_isFile		// if it does not exist, then must be starting new

	// make contents to write to file
	String out = SelectString(overwrite,"\n","<?xml version=\"1.0\" encoding=\"UTF-8\" ?>\n\n")
	out += "<step xmlns=\"http://sector34.xray.aps.anl.gov/34ide:indexResult\">\n"
	out += make_str_xml_line("\t","title",step.title)
	out += make_str_xml_line("\t","sampleName",step.sampleName)
	out += make_str_xml_line("\t","userName",step.userName)
	out += make_str_xml_line("\t","beamline",step.beamline)
	out += make_num_xml_line("\t","scanNum",step.scanNum)
	out += make_str_xml_line("\t","date",step.date)
	out += make_num_xml_line("\t","beamBad",step.beamBad)
	out += make_str_xml_line("\t","monoMode",step.monoMode)
	out += make_num_xml_line("\t","Xsample",step.Xsample)
	out += make_num_xml_line("\t","Ysample",step.Ysample)
	out += make_num_xml_line("\t","Zsample",step.Zsample)
	out += make_num_xml_line("\t","depth",step.depth)
	out += make_num_xml_line("\t","energy",step.energy,unit="keV")
	out += make_num_xml_line("\t","hutchTemperature",step.hutchTemperature)
	out += make_num_xml_line("\t","sampleDistance",step.sampleDistance)

	String str
	Variable i
	for (i=0;i<step.Ndetectors;i+=1)
		out += "\t<detector>\n"				// start of each detector

		out += make_str_xml_line("\t\t","inputImage",step.d[i].inputImage)
		out += make_str_xml_line("\t\t","detectorID",step.d[i].detectorID)
		out += make_num_xml_line("\t\t","exposure",step.d[i].exposure,unit="sec")
		out += make_num_xml_line("\t\t","Nx",step.d[i].Nx)
		out += make_num_xml_line("\t\t","Ny",step.d[i].Ny)
		out += make_num_xml_line("\t\t","totalSum",step.d[i].totalSum,int=1)
		out += make_num_xml_line("\t\t","sumAboveThreshold",step.d[i].sumAboveThreshold,int=1)
		out += make_num_xml_line("\t\t","numAboveThreshold",step.d[i].numAboveThreshold,int=1)
		out += make_str_xml_line("\t\t","geoFile",step.d[i].geoFile)
//		sprintf str, "\t\t<ROI startx=\"%d\" endx=\"%d\" groupx=\"%d\" starty=\"%d\" endy=\"%d\" groupy=\"%d\" /ROI>\n",step.d[i].roi.startx,step.d[i].roi.endx,step.d[i].roi.groupx,step.d[i].roi.starty,step.d[i].roi.endy,step.d[i].roi.groupy
		sprintf str, "\t\t<ROI startx=\"%d\" endx=\"%d\" groupx=\"%d\" starty=\"%d\" endy=\"%d\" groupy=\"%d\"></ROI>\n",step.d[i].roi.startx,step.d[i].roi.endx,step.d[i].roi.groupx,step.d[i].roi.starty,step.d[i].roi.endy,step.d[i].roi.groupy
		out += str

		Variable Npeaks = step.d[i].peaks.Npeaks
		str = "\t\t<peaksXY"
		str += singleXMLkeyValStr("peakProgram",step.d[i].peaks.peakProgram)
		str += singleXMLkeyValNum("minwidth",step.d[i].peaks.minwidth)
		str += singleXMLkeyValNum("threshold",step.d[i].peaks.threshold)
		str += singleXMLkeyValNum("thresholdRatio",step.d[i].peaks.thresholdRatio)
		str += singleXMLkeyValNum("maxRfactor",step.d[i].peaks.maxRfactor)
		str += singleXMLkeyValNum("maxwidth",step.d[i].peaks.maxwidth)
		str += singleXMLkeyValNum("maxCentToFit",step.d[i].peaks.maxCentToFit)
		str += singleXMLkeyValNum("boxsize",step.d[i].peaks.boxsize)
		str += singleXMLkeyValNum("max_number",step.d[i].peaks.max_number)
		str += singleXMLkeyValNum("min_separation",step.d[i].peaks.min_separation)
		str += singleXMLkeyValStr("smooth",SelectString(step.d[i].peaks.smooth,"False","True"))
		str += singleXMLkeyValStr("peakShape",step.d[i].peaks.peakShape)
		str += singleXMLkeyValNum("Npeaks",Npeaks)
		str += singleXMLkeyValNum("executionTime",step.d[i].peaks.executionTime)
		str += ">\n"
		out += str

		Make/N=(Npeaks)/FREE/D vec
		vec = step.d[i].peaks.Xpixel[p]
		out += "\t\t\t<Xpixel>"+vec2str(vec,bare=1,sep=" ",places=7,maxPrint=Npeaks+1)+"</Xpixel>\n"
		vec = step.d[i].peaks.Ypixel[p]
		out += "\t\t\t<Ypixel>"+vec2str(vec,bare=1,sep=" ",places=7,maxPrint=Npeaks+1)+"</Ypixel>\n"
		vec = step.d[i].peaks.Intens[p]
		out += "\t\t\t<Intens>"+vec2str(vec,bare=1,sep=" ",places=8,maxPrint=Npeaks+1)+"</Intens>\n"
		vec = step.d[i].peaks.Integral[p]
		out += "\t\t\t<Integral>"+vec2str(vec,bare=1,sep=" ",places=9,maxPrint=Npeaks+1)+"</Integral>\n"
		vec = step.d[i].peaks.hwhmX[p]
		out += "\t\t\t<hwhmX unit=\"pixel\">"+vec2str(vec,bare=1,sep=" ",maxPrint=Npeaks+1)+"</hwhmX>\n"
		vec = step.d[i].peaks.hwhmY[p]
		out += "\t\t\t<hwhmY unit=\"pixel\">"+vec2str(vec,bare=1,sep=" ",maxPrint=Npeaks+1)+"</hwhmY>\n"
		vec = step.d[i].peaks.tilt[p]
		out += "\t\t\t<tilt unit=\"degree\">"+vec2str(vec,bare=1,sep=" ",places=7,maxPrint=Npeaks+1)+"</tilt>\n"
		vec = step.d[i].peaks.chisq[p]
		out += "\t\t\t<chisq>"+vec2str(vec,bare=1,sep=" ",maxPrint=Npeaks+1)+"</chisq>\n"
		vec = step.d[i].peaks.Qx[p]
		out += "\t\t\t<Qx>"+vec2str(vec,bare=1,sep=" ",places=7,maxPrint=Npeaks+1)+"</Qx>\n"
		vec = step.d[i].peaks.Qy[p]
		out += "\t\t\t<Qy>"+vec2str(vec,bare=1,sep=" ",places=7,maxPrint=Npeaks+1)+"</Qy>\n"
		vec = step.d[i].peaks.Qz[p]
		out += "\t\t\t<Qz>"+vec2str(vec,bare=1,sep=" ",places=7,maxPrint=Npeaks+1)+"</Qz>\n"
		out += "\t\t</peaksXY>\n"
		out += "\t</detector>\n"				// end of each detector
	endfor

	Variable Nindexed, Npatterns=step.index.Npatterns
	if (Npatterns>0)
		str = "\t<indexing"
		str += singleXMLkeyValStr("indexProgram",step.index.indexProgram)
		str += singleXMLkeyValNum("Nindexed",step.index.Nindexed)
		str += singleXMLkeyValNum("Npeaks",step.index.Npeaks)
		str += singleXMLkeyValNum("Npatterns",step.index.Npatterns)
		str += singleXMLkeyValNum("keVmaxCalc",step.index.keVmaxCalc)
		str += singleXMLkeyValNum("keVmaxTest",step.index.keVmaxTest)
		str += singleXMLkeyValNum("angleTolerance",step.index.angleTolerance)
		str += singleXMLkeyValNum("cone",step.index.cone)
		str += singleXMLkeyValStr("hklPrefer",step.index.hklPrefer)
		str += singleXMLkeyValNum("executionTime",step.index.executionTime)
		str += ">\n"
		out += str
		out += "\t\t<!-- Result of indexing. -->\n"

		for (i=0;i<Npatterns;i+=1)
			Nindexed = step.index.pattern[i].Nindexed
			Redimension/N=(Nindexed) vec
			sprintf str, "\t\t<pattern num=\"%d\" rms_error=\"%g\" goodness=\"%g\" Nindexed=\"%d\">\n", step.index.pattern[i].num, step.index.pattern[i].rms_error, step.index.pattern[i].goodness, Nindexed
			out += str
			out += "\t\t\t<recip_lattice unit=\"1/nm\">\n"
			sprintf str "\t\t\t\t<astar>%.9g %.9g %.9g</astar>\n", step.index.pattern[i].astar[0], step.index.pattern[i].astar[1], step.index.pattern[i].astar[2]
			out += str
			sprintf str "\t\t\t\t<bstar>%.9g %.9g %.9g</bstar>\n", step.index.pattern[i].bstar[0], step.index.pattern[i].bstar[1], step.index.pattern[i].bstar[2]
			out += str
			sprintf str "\t\t\t\t<cstar>%.9g %.9g %.9g</cstar>\n", step.index.pattern[i].cstar[0], step.index.pattern[i].cstar[1], step.index.pattern[i].cstar[2]
			out += str
			out += "\t\t\t</recip_lattice>\n"

			out += "\t\t\t<hkl_s>\n"
			vec = step.index.pattern[i].h[p]
			out += "\t\t\t\t<h>"+vec2str(vec,bare=1,sep=" ",places=7,maxPrint=Nindexed+1)+"</h>\n"
			vec = step.index.pattern[i].k[p]
			out += "\t\t\t\t<k>"+vec2str(vec,bare=1,sep=" ",places=7,maxPrint=Nindexed+1)+"</k>\n"
			vec = step.index.pattern[i].l[p]
			out += "\t\t\t\t<l>"+vec2str(vec,bare=1,sep=" ",places=7,maxPrint=Nindexed+1)+"</l>\n"
			out += "\t\t\t</hkl_s>\n"
			out += "\t\t</pattern>\n"
		endfor

		out += "\t\t<xtl>\n"
		if (strlen(step.index.xtl_structureDesc))
			out += "\t\t\t<structureDesc>"+step.index.xtl_structureDesc+"</structureDesc>\n"
		endif
		out += "\t\t\t<SpaceGroup>"+num2istr(step.index.xtl_SpaceGroup)+"</SpaceGroup>\n"
		sprintf str, "\t\t\t<latticeParameters unit=\"nm\">%.9g %.9g %.9g %.9g %.9g %.9g</latticeParameters>\n",step.index.xtl_a,step.index.xtl_b,step.index.xtl_c,step.index.xtl_alpha,step.index.xtl_beta,step.index.xtl_gamma
		out += str
		out += "\t\t</xtl>\n"
		out += "\t</indexing>\n"
	endif
	out += "</step>\n"

	// finshed building te output string, write it to the file
	Variable f
	if (overwrite)
		Open/C="R*ch"/P=$path/Z f as filename
	else
		Open/A/P=$path/Z f as filename
	endif
	if (V_flag)
		printf "ERROR -- in AppendTo3dXML(), unable to open '%s' at path '%s'\r",filename,path
	else
		FBinWrite f, out
		Close f
	endif
	return S_fileName
End
//
Static Function/T make_str_xml_line(pre,key,svalue)
	String pre, key, svalue

	String str=""
	if (strlen(key)>0 && strlen(svalue)>0)
		sprintf str "%s<%s>%s</%s>\n",pre,key,svalue,key
	endif
	return str
End
//
Static Function/T make_num_xml_line(pre,key,value,[unit,int])
	String pre, key
	Variable value
	String unit
	Variable int
	unit = SelectString(ParamIsDefault(unit),unit,"")
	int = ParamIsDefault(int) ? 0 : !(!int)

	String str=""
	if (numtype(value)==0)
		unit = SelectString(strlen(unit),""," unit=\""+unit+"\"")
		if (int)
			sprintf str "%s<%s%s>%.13g</%s>\n",pre,key,unit,value,key
		else
			sprintf str "%s<%s%s>%g</%s>\n",pre,key,unit,value,key
		endif
	endif
	return str
End
//
Static Function/T singleXMLkeyValStr(key,svalue)
	String key, svalue

	String str=""
	if (strlen(key)>0 && strlen(svalue)>0)
		sprintf str, " %s=\"%s\"",key,svalue
	endif
	return str
End
//
Static Function/T singleXMLkeyValNum(key,value)
	String key
	Variable value

	String str=""
	if (strlen(key)>0 && numtype(value)==0)
		sprintf str, " %s=\"%g\"",key,value
	endif
	return str
End

// **************************  End Add Hand Indexed Point to XML ***************************
// *****************************************************************************************




//	<step xmlns="http://sector34.xray.aps.anl.gov/34ide:indexResult">
//	<title>Low dislocation Cu</title>
//	<sampleName>Cu X575</sampleName>
//	<userName>Ben Larson</userName>
//	<beamline>34ID-E</beamline>
//	<scanNum>5647</scanNum>
//	<date>2009-06-10 19:18:50-0600</date>
//	<beamBad>0</beamBad>
//	<CCDshutter>out</CCDshutter>
//	<lightOn>1</lightOn>
//	<monoMode>white slitted</monoMode>
//	<Xsample>2805.05</Xsample>
//	<Ysample>-11527.8</Ysample>
//	<Zsample>-8369.85</Zsample>
//	<depth>nan</depth>
//	<energy unit="keV">12.6866</energy>
//	<hutchTemperature>0.0</hutchTemperature>
//	<sampleDistance>0.0</sampleDistance>
//	<detector>
//		<inputImage>/data34a/Budai_June09/Si_PR5_2/2D_Area1/images/area1_3617.h5</inputImage>
//		<detectorID>PE1621 723-3335</detectorID>
//		<exposure unit="sec">0.1332</exposure>
//		<Nx>2048</Nx>
//		<Ny>2048</Ny>
//		<geoFile>/data34a/Budai_June09/Si_PR5_2/geometry.txt</geoFile>
//		<ROI startx="0" endx="2047" groupx="1" starty="0" endy="2047" groupy="1" />
//		<peaksXY minwidth="0.0113" threshold="1396.8" maxRfactor="0.5" maxwidth="18.0" maxCentToFit="18.0" boxsize="18">
//			<Xpixel>30.423 162.912 169.949 251.233 339.521 453.003 460.515 482.599 536.614 636.415 638.533 698.351 715.685 741.311 773.251 773.273 777.688 842.456 1041.986 1042.703 1122.396 1123.328 1212.051 1372.249 1416.477 1430.587 1483.09 1558.439 1610.788 1647.618 1720.685 1722.404 1798.49 1835.772 1837.875 1879.11</Xpixel>
//			<Ypixel>487.987 2013.307 1309.067 1552.414 967.562 1999.665 767.566 1509.188 1211.719 192.894 869.603 1989.615 678.56 148.604 1459.415 556.537 1154.477 1983.64 577.541 754.507 0.001 1279.518 1549.544 821.413 1964.83 1233.251 1016.48 1609.929 1510.546 617.204 1333.938 212.696 742.754 1181.853 1954.333 258.375</Ypixel>
//			<Intens>5709.0 3989.0 19469.0 3644.0 3840.0 58788.0 3590.0 19347.0 6273.0 30154.0 44619.0 4981.0 14469.0 5703.0 6897.0 2211.0 6969.0 26085.0 22541.0 65535.0 40357.0 1791.0 10409.0 4859.0 18940.0 5997.0 44124.0 2448.0 43131.0 4085.0 65535.0 1881.0 20032.0 6161.0 6067.0 11330.0</Intens>
//		</peaksXY>
//	</detector>
//	<indexing Nindexed="36" Npeaks="36" Npatterns="1" keVmaxCalc="17.0" keVmaxTest="35.0" angleTolerance="0.2" cone="40.0" hklPrefer="2 2 0">
//		<!-- Result of indexing. -->
//		<pattern num="0" rms_error="0.00819" goodness="2079.09" Nindexed="36">
//			<recip_lattice unit="1/nm">
//				<astar>3.6395862 10.6046376 -2.8527945</astar>
//				<bstar>-8.1043903 0.5657142 -8.2366369</bstar>
//				<cstar>-7.4105082 4.5896614 7.6067561</cstar>
//			</recip_lattice>
//		</pattern>
//		<xtl>
//			<structureDesc>Si</structureDesc>
//			<SpaceGroup>227</SpaceGroup>
//			<latticeParameters unit="nm">0.54310206 0.54310206 0.54310206 90 90 90</latticeParameters>
//			<citation>physics.nist.gov/cuu/Constants/index.html, uses 2006, CODATA numbers</citation>
//		   <atom n="1" occ="1" symbol="Si" Z="14" label="Si001">0 0 0</atom>
//		</xtl>
//	</indexing>
//</step>
