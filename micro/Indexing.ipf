#pragma rtGlobals=1		// Use modern global access method.
#pragma ModuleName=Indexing
#pragma IgorVersion = 6.0
#pragma version = 2.65
#include "microGeometry", version>=2.53
#include "WinView", version>=1.92
#include "LatticeSym", version>=3.39

//	with version 2.50, I changed so that yc in the file header always superceedes yc in the geo panel

//	If you want to use more than 250 spots for indexing, make the global variable "root:NumSpotsEulerMax" and set it to the new maximum that you want.
//	The default is 250, to get 300 type into the command line:
//	Variable/G root:NumSpotsEulerMax = 300
//	If you want to remove this limit, either delete this global in the DataBrowser, or set the value to zero or less, or execute the following line:
//	To delete the global use:
//	KillVariables/Z root:NumSpotsEulerMax

//	To change the maximum value for hkl that is used when displaying un-indexed points, set the global "root:MAX_HKL_Displayed", This is done with a line like:
//	Variable/G root:MAX_HKL_Displayed=72
//	To undo this, either delete this global in the DataBrowser, or set the value to zero or less, or execute the following line:
//	To delete the global use:
//	KillVariables/Z root:MAX_HKL_Displayed


Menu "micro"
	SubMenu "Crystal Lattice"
		"Set Crystal Lattice...", /Q, MakeLatticeParametersPanel("")
		"Show Crystal Lattice",/Q,showLattice()
		help={"Shows the crystal structure and lattice that is currently defined"}
		"d[hkl]",/Q,get_dhkl(NaN,NaN,NaN)
		help={"show d-spacing for the hkl reflection"}
	End
	"Load WinView File...", LoadWinViewFile("")
	help={"Load a WinView Image from the file, and then display the image.  If you just want to display a loaded file, go to 'WinView' menu"}
	MenuItemIfWaveClassExists("Remove Background from Image...","speImage","DIMS:2"),RemoveBackgroundImage($"",NaN)
	help={"remove background from an image, do not use on depth sotred images"}
	MenuItemIfWaveClassExists("Fit Peaks on Image...","speImage*","DIMS:2"),FitPeaks($"",NaN,NaN)
	help={"Identify peaks in an Image, does not do any background removal"}
	MenuItemIfWaveClassExists("   TESTING, Fit Peaks step-wise on Image...","speImageNoBkg;speImage","DIMS:2"), FitPeaksStepWise($"",NaN,NaN,NaN,NaN)
	help={"Identify peaks in an Image, does not do any background removal, uses a simple threshold to find peaks"}
	MenuItemIfWaveClassExists("   TESTING, Fit Peaks on Image...","speImageNoBkg;speImage","DIMS:2"),FitPeaksNew($"",NaN,NaN)
	help={"Identify peaks in an Image, does not do any background removal, uses a simple threshold to find peaks"}
	     MenuItemIfWaveClassExists("   OLD Remove Bkg & Fit Peaks on Image...","speImage","DIMS:2"),OLD_removeBkg_FitPeaks($"",NaN,NaN)
	help={"removes bkg and then finds and fits peaks, do not use anymore"}
	  pkLIstWinViewMenuItem("   Reset Fitted Peak list from boxes"),setFittedPeaksFromList($"",NaN,pkList)
	  help={"reset list of peaks for fitting to the boxes on an image plot"}
	MenuItemIfWaveClassExists("Index and Display...","FittedPeakList","DIMS:2,MAXCOLS:11,MINCOLS:11"),IndexAndDisplay($"",NaN,NaN,NaN,NaN,NaN,NaN,NaN)
	help={"Index the identified peaks, using the defined lattice, and then show it all on a plot"}
	   MenuItemIfWaveClassExists("  Refine Strain","IndexedPeakList",""),DeviatoricStrainRefine(NaN,"")
	   help={"compute deviatoric strain for the pattern just idexed"}
	   MenuItemIfWaveClassExists("  Re-Draw hkl tags","IndexedPeakList",""),DisplayResultOfIndexing($"",NaN)
	   help={"for an image plot, clear old hkl tags (if any) and redraw new ones using results of an indexing"}
	   MenuItemIfWaveClassExists("  Clear hkl tags off Graph","IndexedPeakList",""),ClearhklOffGraph("")
	   help={"clear the hkl tags off of an image"}
	   MenuItemIfWaveClassExists(SelectString(exists("MakeStereo")==6 , "  Load Stereographic","  Stereographic Projection..."),"IndexedPeakList",""),StereoOfIndexedPattern($"",NaN)
	   help={"Make a Stereographic projection with the same orientation as the indexed pattern"}
	   MenuItemIfWaveClassExists("  Report Sheet of Indexed Image","IndexedPeakList",""),MakeIndexingReportSheet()
	   help={"if the top window is an images showing the hkl's, this makes a layout with reflections for printing"}
	   MenuItemIfWaveClassExists("  calc energy of [hkl]","IndexedPeakList",""),EnergyOfhkl($"",NaN,NaN,NaN,NaN)
	help={"calculate the energy of an hkl using the current indexing"}
	//	SubMenu("Tables")
	MenuItemsWaveClassOnGraph("Report angle between two cursors","speImage*",""),/Q,AngleBetweenTwoPoints($"",nan, nan,nan,nan)
	SubMenu(MenuItemIfWaveClassExists("Tables","FittedPeakList;IndexedPeakList",""))
		MenuItemIfWaveClassExists("Fitted Peaks","FittedPeakList",""),TableFullPeakList($"")
		help={"Show table of identified peak positions"}
		MenuItemIfWaveClassExists("Indexed Peaks","IndexedPeakList",""),TableFullPeakIndexed($"")
		help={"Show table of indexed peaks"}
		MenuItemIfWaveClassExists("Strain Peaks","IndexedPeakList;StrainPeakList",""),TableFullPeakStrain($"")
		help={"Show table of peaks from the strain refinement"}
	End
	"  Choose Custom Indexing Program...",pickIndexingFunction("")
End
Menu "GraphMarquee", dynamic
	"-"
	MenuItemIndexPeak("Indexed Peak Info"),/Q,IndexedPeakInfoFromMarquee()
End
//
Function/S MenuItemIndexPeak(item)
	String item
	GetMarquee/Z
	if (!V_flag)
		return "("+item
	endif
	if (strlen(GetUserData("","","Indexing"))<1)
		return "("+item
	endif
	return item
End




Static Constant hc = 1.239841857			// keV-nm

Function IndexAndDisplay(FullPeakList,keVmaxCalc,keVmaxTest,angleTolerance,hp,kp,lp,cone)
	Wave FullPeakList				// contains the result of a peak fitting
	Variable keVmaxCalc				// 17, maximum energy to calculate (keV)
	Variable keVmaxTest				// 26, maximum energy to test (keV)  [-t]
	Variable angleTolerance			// 0.25, angular tolerance (deg)
	Variable hp,kp,lp					// preferred hkl
	Variable cone						// angle from preferred hkl, (0 < cone < 180¡)

	if (NumVarOrDefault("root:Packages:geometry:PanelValues:dirty",0))	// geo waiting
		Variable/G root:Packages:geometry:PanelValues:dirty=0
		DoAlert 0,"You have made changes in the geometry panel, deal with them first"
		MakeGeometryParametersPanel("")
		return 1
	endif
	Variable badWave = !WaveExists(FullPeakList)
	if (!badWave)
		badWave = (DimSize(FullPeakList,0)<1 || DimSize(FullPeakList,1)!=11)
	endif
	Variable badNums= !(keVmaxCalc>1 && keVmaxCalc<40) || !(keVmaxTest>1 && keVmaxTest<51)
	badNums += !(angleTolerance>=0.01 && angleTolerance<10)
	badNums += numtype(hp+kp+lp)
	badNums += !(cone>1 && cone<180)
	badNums += (abs(hp)+abs(kp)+abs(lp))==0
	if (badWave || badNums)
		keVmaxCalc = (keVmaxCalc>1 && keVmaxCalc<40) ? keVmaxCalc : 17
		keVmaxTest = (keVmaxTest>1 && keVmaxTest<51) ? keVmaxTest : 26
		angleTolerance = (angleTolerance>=0.01 && angleTolerance<10) ? angleTolerance : 0.25
		hp = numtype(hp) ? 0 : hp
		kp = numtype(kp) ? 0 : kp
		lp = numtype(lp) ? 2 : lp
		cone = (cone>1 && cone<180) ? cone : 72
		String hkl
		sprintf hkl,"%d, %d, %d",hp,kp,lp
		Prompt hkl,"preferred hkl of center"
		Prompt cone,"max angle¡ from hkl prefer to spot"
		String peakListStr=SelectString(badWave,NameOfWave(FullPeakList),"FullPeakList")
//		Prompt peakListStr,"wave with fitted peaks",popup,WaveList("*Peak*",";","DIMS:2,MAXCOLS:11,MINCOLS:11")
//		Prompt peakListStr,"wave with fitted peaks",popup,WaveList_Tags("fittedIgorImage","peakFile","*","DIMS:2,MAXCOLS:11,MINCOLS:11")
		Prompt peakListStr,"wave with fitted peaks",popup,reverseList(WaveListClass("FittedPeakList","*","DIMS:2,MAXCOLS:11,MINCOLS:11"))
		Prompt keVmaxCalc,"max energy for searching (keV)"
		Prompt keVmaxTest,"max energy for matching (keV)"
		Prompt angleTolerance "max angle¡ between peak & hkl"
		DoPrompt/Help="3D-Xray Diffraction[Indexing]" "index",peakListStr,keVmaxCalc,hkl,keVmaxTest,cone,angleTolerance
		if (V_flag)
			return 1
		endif
		sscanf hkl,"%d, %d, %d", hp,kp,lp
		printf "IndexAndDisplay(%s,%g,%g,%g, %d,%d,%d,%g)\r",peakListStr,keVmaxCalc,keVmaxTest,angleTolerance,hp,kp,lp,cone
		Wave FullPeakList=$peakListStr
		badWave = !WaveExists(FullPeakList)
		if (!badWave)
			badWave = (DimSize(FullPeakList,0)<1 || DimSize(FullPeakList,1)!=11)
		endif
		badNums= !(keVmaxCalc>1 && keVmaxCalc<40) || !(keVmaxTest>1 && keVmaxTest<51)
		badNums += !(angleTolerance>=0.01 && angleTolerance<10)
		badNums += numtype(hp+kp+lp) + !(cone>1 && cone<180)
		badNums += (abs(hp)+abs(kp)+abs(lp))==0
		if (badWave || badNums)
			DoAlert 0, "Invalid inputs sent to IndexAndDisplay()"
			return 1
		endif
	endif
	Wave FullPeakIndexed=$(runEulerCommand(FullPeakList,keVmaxCalc,keVmaxTest,angleTolerance,hp,kp,lp,cone))
	if (!WaveExists(FullPeakIndexed))
		return 1
	endif
	String wnote = note(FullPeakIndexed)
	Variable NpatternsFound = NumberByKey("NpatternsFound",wnote,"=")
	Variable Nindexed = NumberByKey("Nindexed",wnote,"=")
	Variable NiData = NumberByKey("NiData",wnote,"=")
	Variable executionTime = NumberByKey("executionTime",wnote,"=")
	Variable rms_error = NumberByKey("rms_error0",wnote,"=")

	if (ItemsInList(GetRTStackInfo(0))<2 || stringmatch(StringFromList(0,GetRTStackInfo(0)),"IndexButtonProc"))
		if (executionTime<60)
			printf "from Euler, found %d patterns, indexed %d out of %d spots with rms=%g¡ in %g sec\r",NpatternsFound,Nindexed,NiData,rms_error,executionTime
		else
			printf "from Euler, found %d patterns, indexed %d out of %d spots  with rms=%g¡ in %s (%g sec)\r",NpatternsFound,Nindexed,NiData, rms_error,Secs2Time(executionTime,5,1),executionTime
		endif
		DoAlert 1, "Examine fit on a picture"
		if (V_flag==1)
			DisplayResultOfIndexing(FullPeakIndexed,NaN)
		endif
	endif
	return 0
End
//
Function ClearhklOffGraph(win)
	String win
	String item,list = AnnotationList("")
	Variable i
	for (i=0,item=StringFromList(0,list); strlen(item); i+=1,item=StringFromList(i,list))
		if (strsearch(item, "hkl",0)==0)
			Tag/W=$win/K/N=$item
		endif
	endfor
	list = GetUserData("","","Indexing")
	SetWindow kwTopWin,userdata(Indexing)=RemoveByKey("FullPeakIndexed",list,"=")
End
//
Function DisplayResultOfIndexing(FullPeakIndexed,pattern)
	Wave FullPeakIndexed				// contains the result of indexing peaks
	Variable pattern					// pattern in FullPeakIndexed to draw

	Variable i, badWave = !WaveExists(FullPeakIndexed)
	if (!badWave)
		badWave = (DimSize(FullPeakIndexed,0)<1 || DimSize(FullPeakIndexed,1)!=11)
	endif
	if (badWave)
		String wList = reverseList(WaveListClass("IndexedPeakList","*",""))
		i = ItemsInList(wList)
		if (i==1)						// only one, so choose it
			Wave FullPeakIndexed = $(StringFromList(0,wList))
		elseif(i>0)						// more than 1, so ask
			String wName=""
			Prompt wName,"Indexed Peak List",popup,wList
			DoPrompt "Indexed Peak List",wName
			if (V_flag)
				return 1
			endif
			Wave FullPeakIndexed = $wName
		endif
	endif
	badWave = !WaveExists(FullPeakIndexed)
	if (!badWave)
		badWave = (DimSize(FullPeakIndexed,0)<1 || DimSize(FullPeakIndexed,1)!=11)
	endif
	String wnote = note(FullPeakIndexed)
	Variable NpatternsFound = max(DimSize(FullPeakIndexed,2),1)
	NpatternsFound = min(NpatternsFound,NumberByKey("NpatternsFound",wnote,"="))
	if (NpatternsFound<1)
		return 1
	elseif (!(NpatternsFound!=1))
		pattern = 0
	elseif (!(pattern>=0))
		pattern = 0
		Prompt pattern, "pattern number to fit [0,"+num2istr(NpatternsFound-1)+",]",popup,expandRange("0-"+num2istr(NpatternsFound-1),";")
		DoPrompt "pattern number",pattern
		if (V_flag)
			return 1
		endif
		pattern = min(NpatternsFound-1,pattern-1)
	endif
	if (!(pattern>=0))
		return 1
	endif
	// get a graph of the image at the front, with boxes on to the front
	Wave image = $""					// init to NULL
 	String imageName					// name of image to look at
	String win							// name of window displaying image
	imageName = StringByKey("fittedIgorImage",wnote,"=")

	if (exists(imageName)==1)
		Wave image = $imageName
		win = FindGraphWithImage(image)
	endif
	if (strlen(win)<1)					// found fitted image, but no window, see if raw window is open
		Wave rawimage = $(StringByKey("rawIgorImage",wnote,"="))
		win = FindGraphWithImage(rawimage)
	endif
	if (!WaveExists(image))			// could not find the fitted image, try for raw image
		imageName = StringByKey("rawIgorImage",wnote,"=")
		Wave image = $imageName
		win = FindGraphWithImage(image)
	endif

	if (!WaveExists(image))			// could not find the image
		return 1
	endif
	if (strlen(win))
		DoWindow/F $win				// bring window with image to the top
	else
		Graph_imageMake(image,1)	// no window exists, so make one showing image
  	endif
	NVAR boxesOn=boxesOn
	boxesOn = 0
	ButtonProc_Boxes("")				// make sure that boxes are on
	ClearhklOffGraph("")

	String peakListName = StringByKey("peakListWave", wnote,"=")
	if (strlen(peakListName))
		SetWindow kwTopWin userdata(FitPeaks )="FullPeakList="+peakListName
	endif

	Variable wid = (NumberByKey("maxPeakWidth",wnote,"=")+NumberByKey("minSpotSeparation",wnote,"="))/2
	wid = numtype(wid) ? NumberByKey("minSpotSeparation",wnote,"=") : wid
	//	wid = numtype(wid) ? 30 : wid		// width of the cross
	wid = numtype(wid) ? 20 : wid		// width of the cross
	Variable N=DimSize(FullPeakIndexed,0)
	Variable px,py
	SetDrawLayer UserFront			// draw the X's
	for (i=0;i<N;i+=1)
		px = FullPeakIndexed[i][9][pattern]
		py = FullPeakIndexed[i][10][pattern]
		if (!numtype(px+py))
			DrawMarker(px,py,wid,wid,"X",color="55000,5000,0")
//			DrawX(px,py,wid,wid)
			DrawhklTags(i,FullPeakIndexed[i][3][pattern],FullPeakIndexed[i][4][pattern],FullPeakIndexed[i][5][pattern],px,py,wid,0)
		endif
	endfor
	TextBox/K/N=indexedPeakInfo
	TextBox/C/N=indexedPeakInfo/F=0 ""
	String list = GetUserData("","","Indexing")
	list = ReplaceStringByKey("FullPeakIndexed",list,GetWavesDataFolder(FullPeakIndexed,2),"=")
	list = ReplaceNumberByKey("patternNum",list,pattern,"=")
	SetWindow kwTopWin,userdata(Indexing)=list
	return 0
End
//
Static Function DrawhklTags(i,h,k,l,px,py,dx,dy)
	Variable i
	Variable h,k,l				// hkl value
	Variable px,py				// pixel position
	Variable dx,dy				// off set of label from pixel center
	String str
	str = "\\F'Comic Sans MS'\\Zr075\\[0"		// use for gfMult=100
	str += SelectString(h<0,"","\\S \\f01Ñ\\f00\\M\\X0")+" "+num2istr(abs(h))
	str += SelectString(k<0,"","\\[1\\S\\f01 Ñ\\f00\\M\\X1")+" "+num2istr(abs(k))
	str += SelectString(l<0,"","\\[2\\S \\f01Ñ\\f00\\M\\X2")+" "+num2istr(abs(l))
	Wave image = ImageNameToWaveRef("",StringFromList(0,ImageNameList("",";")))
	px = limit(px,0,DimSize(image,0)-1)
	py = limit(py,0,DimSize(image,1)-1)
	Variable p = round(px + DimSize(image,0)*round(py))
	String name = "hkl"+num2istr(i)
	String wname=NameOfWave(image)
	Tag/C/N=$name/I=1/F=0/G=(55000,5000,0)/B=1/A=MB/L=0/X=1/Y=0.5 $wname,p,str
	return 0
End
//
Function IndexedPeakInfoFromMarquee()	// find the energy of a spot in the marquee
	String imageName = StringFromList(0,ImageNameList("",";"))
	Wave image = ImageNameToWaveRef("",imageName)
	STRUCT box b
	if (pointRangeOfMarquee("",imageName,b)!=0)
		return NaN
	endif

	String list = GetUserData("","","Indexing")
	Variable ip = NumberByKey("patternNum",list,"=")		// pattern number
	Wave FullPeakIndexed=$(StringByKey("FullPeakIndexed",list,"="))
	if (!WaveExists(FullPeakIndexed))
		DoAlert 0,"FullPeakIndexed cannot be found"
		return NaN
	endif
	ip = (ip<0 || numtype(ip)) ? 0 : ip
	ip = limit(ip,0,max(DimSize(FullPeakIndexed,2)-1,0))
	Variable px,py,i,N=DimSize(FullPeakIndexed,0)
	for (i=0;i<N;i+=1)				// first find the peak
		px = FullPeakIndexed[i][9][ip]
		py = FullPeakIndexed[i][10][ip]
		if (b.xlo<px && px<b.xhi && b.ylo<py && py<b.yhi)
			break
		endif
	endfor
	if (i>=N)
		DoAlert 0,"No indexed peak is in the Marquee"
		return NaN
	endif
	Variable h=FullPeakIndexed[i][3][ip], k=FullPeakIndexed[i][4][ip], l=FullPeakIndexed[i][5][ip]
	Variable Intens=FullPeakIndexed[i][6][ip], keV=FullPeakIndexed[i][7][ip],angleErr=FullPeakIndexed[i][8][ip]
	Variable Np = NumberByKey("NpatternsFound",note(FullPeakIndexed),"=")
	if (Np>1)
		printf "Indexed peak %d at pixel(%g, %g),     in pattern %d of %d, hkl=(%d %d %d),   Intensity=%g,  Energy=%.4f keV,  angleErr=%.3f (deg)\r",i,px,py,ip,Np,h,k,l,Intens,keV,angleErr
	else
		printf "Indexed peak %d at pixel(%g, %g),     hkl=(%d %d %d),   Intensity=%g,  Energy=%.4f keV,  angleErr=%.3f (deg)\r",i,px,py,h,k,l,Intens,keV,angleErr
	endif
End
// This puts up a window with info about an indexed peak on the current graph, when you shift-click on a spot
//
//	if (exists("getIndexedPeakInfoHook")==6)
//		SetWindow Graph1 hook(peakInfo)=getIndexedPeakInfoHook
//	endif
//
//
Function getIndexedPeakInfoHook(s)	// Command=fitted peak,  Shift=Indexed peak,  Command-Shift=strained peak,  &  Command-mouseDown=follow mouse
	STRUCT WMWinHookStruct &s

	// use shift to look at the indexing result, use command to look at fitted peaks
	String win = (s.winName)
	if (!(s.eventMod&10 && (s.eventCode==4 || s.eventCode==3)) || s.keycode)		// eventMod==2 is shift key, 8 is command, 1 is mouse, and no regular key held down
		Tag/K/N=indexedPeakInfo/W=$win
		return 0								// shift key not down, ignore
	endif

	Variable useIndex = s.eventMod&2		// eventMod==2 is shift key,  use info from indexed peaks
	Variable useFitted = s.eventMod&8		// eventMod==8 is command key, use info from fitted peaks
	Variable useStrain = useIndex && useFitted
	useFitted = useFitted && !useStrain		// both means use strain
	Variable useMouse = s.eventMod==9		// evendMod==8+1, is command & mouse down, use mouse position
	if (useMouse)
		useIndex = 0
		useFitted = 0
		useStrain = 0
	endif

	String imageName = StringFromList(0,ImageNameList(win,";"))
	Wave image = ImageNameToWaveRef(win,StringFromList(0,ImageNameList(win,";")))
	if (!WaveExists(image))
		return 1
	elseif (WaveDims(image)!=2)
		return 1
	endif
	GetWindow $win psize
	Variable vert = limit((s.mouseLoc.v-V_top)/(V_bottom-V_top),0,1)	// fractional position on graph
	Variable horiz = limit((s.mouseLoc.h-V_left)/(V_right-V_left),0,1)

	Variable mx,my							// pixel position at the mouse-click
	GetAxis/W=$win/Q bottom
	if (V_flag)
		return 1
	endif
	mx = (V_max-V_min)*horiz + V_min
	GetAxis/W=$win/Q left
	if (V_flag)
		return 1
	endif
	my = (V_min-V_max)*vert + V_max

	String tagStr="", wnote, str
	Variable h,k,l, angleErr,keV,SpaceGroup
	Variable dist2, m
	Variable px,py,i,N
	Variable ip = NumberByKey("patternNum",GetUserData("","","Indexing"),"=")	// pattern number
	ip = (ip<0 || numtype(ip)) ? 0 : ip

	Wave FullPeakIndexed=$(StringByKey("FullPeakIndexed",GetUserData("","","Indexing"),"="))
	if (WaveExists(FullPeakIndexed))
		ip = limit(ip,0,max(DimSize(FullPeakIndexed,2)-1,0))
		Wave PeaksForStrain=$( GetWavesDataFolder(FullPeakIndexed,1)+"PeaksForStrain")
	endif
	useStrain = useStrain && (modDate(PeaksForStrain) >= modDate(FullPeakIndexed))	// only use strain if not stale
	useStrain = useStrain && (NumberByKey("patternNum", note(PeaksForStrain),"=")==ip)
	useIndex = useIndex && !useStrain

	if (useFitted || useMouse)				// examining fitted peaks (or mouse position)
		Wave FullPeakList=$(StringByKey("FullPeakList",GetUserData("","","FitPeaks"),"="))
		if (!WaveExists(FullPeakList) && useFitted)
			return 0
		endif

		Variable fwx=NaN, fwy=NaN, xc=NaN
		if (useFitted)
			dist2=Inf
			m=-1								// find the indexed peak closest to the mouse-click
			N=DimSize(FullPeakList,0)
			for (i=0;i<N;i+=1)				// search through all fitted peaks
				px = FullPeakList[i][0]		// FullPeakList contains binned pixels (of fitted peaks)
				py = FullPeakList[i][1]
				if ((px-mx)^2+(py-my)^2 < dist2)
					m = i
					dist2 = (px-mx)^2+(py-my)^2
				endif
			endfor
			if (m<0)
				return 1						// no indexed peak found
			endif
			px = FullPeakList[m][0]			// binned peaks
			py = FullPeakList[m][1]
			fwx=FullPeakList[m][4]
			fwy=FullPeakList[m][5]
			xc=FullPeakList[m][8]
		elseif (useMouse)
			px = mx
			py = my
		else
			Abort "not useFitted or useMouse"
		endif

		STRUCT microGeometry geo
		if (FillGeometryStructDefault(geo))					//fill the geometry structure with test values
			DoAlert 0, "no geometry structure found, did you forget to set it?"
			return 1
		endif
		if (WaveExists(FullPeakList))
			wnote = note(FullPeakList)
		else
			wnote = note(image)		
		endif
		Variable ddLocal = geo.ddOffset + NumberByKey("CCDy",wnote,"=")	// get dd from the image (in wave note) if it is there
		geo.dd = ddLocal>0 ? ddLocal :  geo.dd
		Variable ycLocal = NumberByKey("yc",wnote,"=")					// get ycLocal from the image (in wave note) if it is there
		geo.ycent = ycLocal>0 ? ycLocal :  geo.ycent

		Make/N=3/O/D PeakInfoHook_qhatW,PeakInfoHook_qhatBL
		Wave qhatW = PeakInfoHook_qhatW,qBL=PeakInfoHook_qhatBL
		Variable startx,groupx, starty,groupy					// ROI of the original image
		startx = NumberByKey("startx",wnote,"=")
		groupx = NumberByKey("groupx",wnote,"=")
		starty = NumberByKey("starty",wnote,"=")
		groupy = NumberByKey("groupy",wnote,"=")
		startx = numtype(startx) ? 1 : startx
		groupx = numtype(groupx) ? 1 : groupx
		starty = numtype(starty) ? 1 : starty
		groupy = numtype(groupy) ? 1 : groupy
		Variable pxUnb = (startx-1) + groupx*px + (groupx-1)/2		// change to un-binned pixels
		Variable pyUnb = (starty-1) + groupy*py + (groupy-1)/2		// pixels are still zero based
		Variable theta = pixel2q(geo,pxUnb,pyUnb,qhatW)		// get theta, and q^ in Wenge system
		Variable/C pz = peakcorrection2($("root:Packages:geometry:xymap"),pxUnb,pyUnb)
		wenge2IdealBeamline(geo,qhatW,qBL)					// rotate qhatW from Wenge to Ideal Beam-line
		if (useFitted)
			sprintf tagStr,"\\Zr090Fitted peak position (%.2f, %.2f)\rFWHM: x=%.2f, y=%.2f,  x-corr=%.3f\r\\F'Symbol'q\\F]0 = %.4f\\F'Symbol'°\\F]0",px,py,fwx,fwy,xc,theta*180/PI
		elseif (useMouse)
			sprintf tagStr,"\\Zr090Mouse position (%.2f, %.2f)\r\\F'Symbol'q\\F]0 = %.4f\\F'Symbol'°\\F]0,     distort=%.2f px",px,py,theta*180/PI,cabs(pz)
			sprintf str,"\r\[0\\Zr075q\X0\y+20^\y-20\BBL\\M\\Zr075 = {%.3f, %.3f, %.3f}\M\\Zr090",qBL[0],qBL[1],qBL[2]
			tagStr += str
		endif

		if (WaveExists(FullPeakIndexed))						// an indexation exists, so try to add some hkl info too
			wnote = note(FullPeakIndexed)
			str = StringByKey("recip_lattice"+num2istr(ip),wnote,"=")
			Variable r00,r10,r20,r01,r11,r21,r02,r12,r22
			sscanf str, "{{%g, %g, %g}{ %g, %g, %g}{%g, %g,%g}}",r00,r10,r20,r01,r11,r21,r02,r12,r22
			if (V_flag==9)
				Make/N=3/O/D PeakInfoHook_recip
				Wave recip = PeakInfoHook_recip
				recip[0][0]= {r00,r10,r20}
				recip[0][1]= {r01,r11,r21}
				recip[0][2]= {r02,r12,r22}
				MatrixOp/O PeakInfoHook_hkl = Inv(recip) x qBL	// go from q in Ideal beam line to hkl (but wrong length)
				Wave hkl = PeakInfoHook_hkl
				sprintf str, "hkl = (%.3f, %.3f, %.3f)",hkl[0]*24,hkl[1]*24,hkl[2]*24
				tagStr += "\r"+str
				Make/N=3/O/D PeakInfoHook_hkli=hkl			// try to make a nice integral hkl
				Wave hkli = PeakInfoHook_hkli
				CloseAllowedhkl(hkli)
				MatrixOp/O PeakInfoHook_xyz = recip x hkli
				keV = hc*norm(PeakInfoHook_xyz)/(4*PI*sin(theta))
				sprintf str, "hkl = (%g %g %g),  %.4f keV",hkli[0],hkli[1],hkli[2],keV
				tagStr += "\r"+str
				if (useFitted)
					MatrixOp/O PeakInfoHook_hkl = recip x hkli	// go from integral hkl to Q in BL
					Variable dtheta = angleVec2Vec(hkl,qBL)	// angle between fitted peak and indexed peak
					sprintf str,"\r\\F'Symbol'Dq\\F]0 = %.4f\\F'Symbol'°\\F]0,   distort=%.2f px",dtheta,cabs(pz)
					tagStr += str
				endif
			endif
		endif
		KillWaves/Z PeakInfoHook_qhatW,PeakInfoHook_qhatBL, PeakInfoHook_recip, PeakInfoHook_hkl, PeakInfoHook_hkli, PeakInfoHook_xyz
	elseif (useIndex || useStrain)								// examining indexed peaks
		if (!WaveExists(FullPeakIndexed))
			return 1
		endif
		dist2=Inf
		m=-1														// find the indexed peak closest to the mouse-click
		N=DimSize(FullPeakIndexed,0)
		for (i=0;i<N;i+=1)										// search through all indexed peaks
			px = FullPeakIndexed[i][9][ip][ip]					// binned pixels
			py = FullPeakIndexed[i][10][ip]
			if ((px-mx)^2+(py-my)^2 < dist2)
				m = i
				dist2 = (px-mx)^2+(py-my)^2
			endif
		endfor
		if (m<0)
			return 1												// no indexed peak found
		endif
		h=FullPeakIndexed[m][3][ip]
		k=FullPeakIndexed[m][4][ip]
		l=FullPeakIndexed[m][5][ip]
		if (useIndex)
			px = FullPeakIndexed[m][9][ip]
			py = FullPeakIndexed[m][10][ip]
			keV=FullPeakIndexed[m][7][ip]
			angleErr=FullPeakIndexed[m][8][ip]
			SpaceGroup=NumberByKey("SpaceGroup",note(FullPeakIndexed),"=")
			if (ip>0)
				sprintf str,"hkl=(%d %d %d),   %.4f keV\rpixel(%.2f, %.2f),   #%d, %d\rangleErr=%.4f (deg)",h,k,l,keV,px,py,m,ip,angleErr
			else
				sprintf str,"hkl=(%d %d %d),   %.4f keV\rpixel(%.2f, %.2f),   #%d\rangleErr=%.4f (deg)",h,k,l,keV,px,py,m,angleErr
			endif
			tagStr = "\\Zr090Indexed peak position\r" + str
			tagStr += SelectString(numtype(SpaceGroup),"\r"+getSymString(SpaceGroup)+"    Space Group #"+num2istr(SpaceGroup),"")
		else														// show strain data
			px = PeaksForStrain[m][11]
			py = PeaksForStrain[m][12]
			sprintf str,"pixel(%.2f, %.2f)\r%.4f keV\rhkl=(%d %d %d),   #%d",px,py, PeaksForStrain[m][10],h,k,l,m
			tagStr = "\\Zr090Strained peak position\r" + str
		endif
	endif

	px = limit(px,0,DimSize(image,0)-1)						// needed in case (px,py) is outside the image
	py = limit(py,0,DimSize(image,1)-1)
	Variable index = round(px) + DimSize(image,0)*round(py)	// convert px,py back to index into image
	GetAxis/W=$win/Q bottom
	horiz = (px-V_min)/ (V_max-V_min)
	GetAxis/W=$win/Q left
	vert = (py-V_max)/(V_min-V_max)
	Variable x0 = horiz<0.5 ? 5 : -5,  y0 = vert<0.5 ? -5 : 5
	String anchor = SelectString(horiz<0.5,"R","L")+ SelectString(vert<0.5,"B","T")
	Tag/C/N=indexedPeakInfo/W=$win/A=$anchor/F=2/L=2/X=(x0)/Y=(y0)/P=1 $imageName,index,tagStr
	DoUpdate
	return 1														// 1 means do not send back to Igor for more processing
End


Function MakeIndexingReportSheet()	// make a report of the indexation using the top graph
	String gName = StringFromList(0,WinList("*",";","WIN:1"))
	String list = GetUserData(gName,"","Indexing")
	Wave FullPeakIndexed = $(StringByKey("FullPeakIndexed",list,"="))
	Variable ip = NumberByKey("patternNum",list,"=")
	ip = numtype(ip) || ip<0 ? 0 : ip
	if (!WaveExists(FullPeakIndexed))
		DoAlert 0,"cannot find indexed peaks on the top graph"
		return 1
	endif

	String imageName = StringFromList(0,ImageNameList(gName,";"))
	Wave image = ImageNameToWaveRef(gName,imageName)
	String str, textOut="\\Zr150"+imageName+SelectString(ip,"","  pattern "+num2istr(ip))
	if (WaveExists(image))
		String dateExposed="", imageFileName=""
		dateExposed = StringByKey("dateExposed",note(image),"=")
		imageFileName = StringByKey("imageFileName",note(image),"=")
		sprintf str,"\\Zr075   from file '%s'   taken on  %s",imageFileName,dateExposed
		textOut += str
	endif

	Wave FullPeakIndexed=$(GetWavesDataFolder(image,1)+"FullPeakIndexed")
	str = ParagraphsFromIndexing(FullPeakIndexed,ip)
	String par1 = StringFromLIst(0,str), par2 = StringFromLIst(1,str)
	NewLayout/C=1/K=1/P=Portrait
	AppendLayoutObject/F=0/R=(43,23,582,506) graph $gName
	if (strlen(par1))
		TextBox/N=text0/C/F=0/T={31,80,126,177,216,252,288,324,360}/A=LT/X=6/Y=69.67 par1
	endif
	if (strlen(par2))
		TextBox/N=text1/C/F=0/T={31,80,126,177,216,252,288,324,360}/A=RT/X=6/Y=69.67 par2
	endif
	TextBox/N=text2/F=0/X=4.70/Y=66.80 textOut
	str = StringByKey("structureDesc",note(FullPeakIndexed),"=")
	if (strlen(str))
		TextBox/C/N=textDesc/B=1/F=0/A=RT/X=5/Y=0 "\\Zr150"+str
	endif
	Textbox/C/N=stamp0/F=0/A=RB/X=0.1/Y=0.1 "\\Z06\\{\"%s %s\",date(), time()}"
	Textbox/C/N=stamp1/F=0/A=LB/X=0.1/Y=0.1 "\\Z06\\{\"%s\",CornerStamp1_()}"+":"+WinName(0, 1)
End
//
Static Function/S ParagraphsFromIndexing(fpi,ip)
	Wave fpi									// usually FullPeakIndexed
	Variable ip								// pattern number (usually 0)

	String line, par1="",par2=""
	String topLine = "\\Z09#\t(hkl)\tIntens\tkeV\t err(deg)"

	Variable i, N = DimSize(fpi,0)
	for (i=0;i<min(N,20);i+=1)			// set par1
		if (numtype(fpi[i][0][ip]))
			continue
		endif
		sprintf line,"\r%02d	(%s)	%.g	%.3f	%.3f",i,hkl2str(fpi[i][3][ip],fpi[i][4][ip],fpi[i][5][ip]),fpi[i][6][ip],fpi[i][7][ip],fpi[i][8][ip]
		par1 += line
	endfor
	if (N>i)
		N = min(N,i+20)
		for (;i<N;i+=1)						// set par2
			if (numtype(fpi[i][0][ip]))
				continue
			endif
			sprintf line,"\r%02d	(%s)	%.g	%.3f	%.3f",i,hkl2str(fpi[i][3][ip],fpi[i][4][ip],fpi[i][5][ip]),fpi[i][6][ip],fpi[i][7][ip],fpi[i][8][ip]
			par2 += line
		endfor
	endif
	par1 = SelectString(strlen(par1),"",topLine+par1)
	par2 = SelectString(strlen(par2),"",topLine+par2)
	return par1+";"+par2
End


Function/T pickIndexingFunction(path)
	String path						// path to folder on disk to be searched
	if (strlen(path)<1)
		path = ParseFilePath(1,FunctionPath("runEulerCommand"),":",1,0)// default path to the executables
	endif

	String defaultExe,exe=""
	Variable isMac = stringmatch(igorInfo(2),"Macintosh")
	if (isMac && stringmatch(igorinfo(4),"PowerPC"))
		defaultExe = "Euler_ppc"
		isMac = 1
	elseif (isMac && stringmatch(igorinfo(4),"Intel"))
		defaultExe = "Euler_i386"
		isMac = 2
	else
		defaultExe = "Euler.exe"
	endif
	exe = StrVarOrDefault("root:Packages:micro:indexingExecutableMac","")

	// get list of possibilities
	NewPath/O/Q/Z IndexingSearchPath, path
	if (V_Flag)
		DoAlert 0,"Unable to set path to "+path
		return ""
	endif
	String list=""
	if (isMac)
		list = IndexedFile(IndexingSearchPath,-1,"????")
		Variable i, bad
		String fname
		for (i=ItemsInList(list)-1;i>=0;i-=1)
			fname = StringFromList(i,list)
			GetFileFolderInfo/Q/Z path+fname
			bad = (strlen(S_fileType+S_creator) || V_isInvisible || !V_isFile || V_isAliasShortcut || V_isStationery || V_logEOF<300e3)
			bad = bad || stringmatch(fname,"*.exe")				// on a Mac, exclude .exe files
			if (isMac ==1)											// for ppc, exclude files ending in i386 or intel
				bad = bad || stringmatch(fname,"*i386") || stringmatch(fname,"*intel")
			elseif (isMac==2)										// for intel, exclude files ending with ppc
				bad = bad || stringmatch(fname,"*ppc")	
			endif
			if (bad)
				list = RemoveFromList(fname,list)
			endif
		endfor
	else
		list = IndexedFile(IndexingSearchPath,-1,".exe")
	endif
	if (ItemsInList(list)<1)			// nothing to select from
		DoAlert 0,"No executables found, do nothing."
		return ""
	endif

	// choose from list of possibilities
	Prompt exe,"name of indexing program",popup,"_default_;"+list
	DoPrompt/HELP="Choose a different program for indexing.\r  Euler is fast\r  IndexAll is slow\r  default is usually best." "indexing program",exe
	if (V_flag)
		return ""
	endif

	// set global, or delete it if default
	if (stringmatch(exe,defaultExe) || stringmatch(exe,"_default_")|| strlen(exe)<1)	// set to use default
		KillStrings/Z root:Packages:micro:indexingExecutableMac
 		printf "\r\r ==========  Setting Indexing Executable to use default  ==========\r\r\r"
	else
		String/G root:Packages:micro:indexingExecutableMac = exe
		printf "\r\r ==========  Setting Indexing Executable to '%s'  ==========\r\r\r",exe
	endif
	return exe
End


#if stringmatch(igorInfo(2),"Macintosh")				// Mac version of runEulerCommand()
//
// using the result from FitPeaks() run Euler, and read in the results from the index file
Function/S runEulerCommand(FullPeakList,keVmaxCalc,keVmaxTest,angleTolerance,hp,kp,lp,cone)
	Wave FullPeakList
	Variable keVmaxCalc				// 17, maximum energy to calculate (keV)
	Variable keVmaxTest				// 26, maximum energy to test (keV)  [-t]
	Variable angleTolerance			// 0.25, angular tolerance (deg)
	Variable hp,kp,lp					// preferred hkl
	Variable cone						// angle from preferred hkl, (0 < cone < 180¡)

	Variable badNums= !(keVmaxCalc>1 && keVmaxCalc<40) || !(keVmaxTest>1 && keVmaxTest<51)
	badNums += !(angleTolerance>=0.01 && angleTolerance<10)
	badNums += numtype(hp+kp++lp)
	badNums += !(cone>1 && cone<180)
	if (badNums)
		DoAlert 0, "Invalid inputs sent to runEulerCommand()"
		printf "runEulerCommand(%s,%g,%g,%g)\r",FullPeakList,keVmaxCalc,keVmaxTest,angleTolerance
		return ""
	elseif (!WaveExists(FullPeakList))
		DoAlert 0, "the in put wave does not exist in runEulerCommand()"
		return ""
	elseif (DimSize(FullPeakList,0)<1 || DimSize(FullPeakList,1)!=11)
		DoAlert 0, "Full peak list '"+NameOfWave(FullPeakList)+"' is empty or the wrong size"
		return ""
	endif
	Variable isMac = stringmatch(igorInfo(2),"Macintosh")

	// first write the command file to drive Euler
	String upath=SpecialDirPath("Temporary",0,1,0)	// local path (probably unix path) for the command line
	String mpath=SpecialDirPath("Temporary",0,0,0)	// mac style path for use only within Igor
	PathInfo home
	if (V_flag)									// the path "home" exists, use it
		mpath = S_path
		upath = HFSToPosix("",S_path,1,1)
	endif
	NewPath/O/Q/Z EulerCalcFolder, mpath		// need a new path name since "home" may not exist
	if (V_flag)
		DoAlert 0, "Unable to create path to do Euler calculation, try saving the experiment first."
		return ""
	endif
	String peakFile="generic_Peaks.txt"		// name of file with input peak positions
	if(FullPeakList2Qfile(FullPeakList,peakFile,"EulerCalcFolder"))// convert peaks to a Qlist+intens, and write to a file
		return ""									// nothing happened
	endif

	// find the full path name of the Euler executable
	String name,EulerPath = ParseFilePath(1,FunctionPath("runEulerCommand"),":",1,0)// path to the Euler executable
	if (isMac && stringmatch(igorinfo(4),"PowerPC"))	// find the default name for this architecture
		name = "Euler_ppc"
	elseif (isMac && stringmatch(igorinfo(4),"Intel"))
		name = "Euler_i386"
	else
		name = "Euler"
	endif
	name = StrVarOrDefault("root:Packages:micro:indexingExecutableMac",name)	// over ride default if indexingExecutableMac is set
	GetFileFolderInfo/Q/Z EulerPath+name
	EulerPath += SelectString(V_Flag==0 && V_isFile,"Euler",name)	// use just plane Euler if Euler_ppc or Euler_i386 does not exist
	EulerPath = HFSToPosix("",EulerPath,1,1)							// convert from Mac to unix style paths
	if (strlen(EulerPath)<1)
		DoAlert 0, "cannot find the executable '"+name+"'"
		return ""
	endif

	Variable maxSpots=NumVarOrDefault("root:NumSpotsEulerMax",-1)	//	-n max num. of spots from data file to use, default is 250
	maxSpots = ((maxSpots>2) && numtype(maxSpots)==0) ? maxSpots : -1
	String cmd
	if (maxSpots>2)
		sprintf cmd "do shell script \"cd \\\"%s\\\" ; \\\"%s\\\" -k %g -t %g -a %g -h %d %d %d -c %g -n %d -f %s\"",upath,EulerPath,keVmaxCalc,keVmaxTest,angleTolerance,hp,kp,lp,cone,maxSpots,peakFile
	else
		sprintf cmd "do shell script \"cd \\\"%s\\\" ; \\\"%s\\\" -k %g -t %g -a %g -h %d %d %d -c %g -f %s\"",upath,EulerPath,keVmaxCalc,keVmaxTest,angleTolerance,hp,kp,lp,cone,peakFile
	endif
	ExecuteScriptText cmd
	if (ItemsInList(GetRTStackInfo(0))<2)		// print everything if run from command line
		printf S_value
	endif
	Variable i = strsearch(S_value,"writing output to file '",0)
	if (i<0)
		DoAlert 0, "failure in runEulerCommand()"
		print "\r\r"
		print cmd
		print "\r\r"
		print S_value
		return ""								// there is no index file
	endif
	String line = S_value[i,Inf]
	i = strsearch(line,"'",0)+1
	line = line[i,Inf]
	String indexFile = line[0,strsearch(line,"'",0)-1]
	if (!strlen(indexFile))
		DoAlert 0, "in runEulerCommand(), cannot find the index file"
		print "\r\r"
		print S_value
		return ""								// there is no index file
	endif
	return readIndexFile(indexFile,"EulerCalcFolder")	// read an index file and create an array to hold the results, return array name
End
//
//
#elif stringmatch(igorInfo(2),"Windows")			// Windows version of runEulerCommand()
//
//	This is just a stub until I get a final version
Function/S runEulerCommand(FullPeakList,keVmaxCalc,keVmaxTest,angleTolerance,hp,kp,lp,cone)
	Wave FullPeakList
	Variable keVmaxCalc				// 17, maximum energy to calculate (keV)
	Variable keVmaxTest				// 26, maximum energy to test (keV)  [-t]
	Variable angleTolerance			// 0.25, angular tolerance (deg)
	Variable hp,kp,lp					// preferred hkl
	Variable cone						// angle from preferred hkl, (0 < cone < 180¡)

	DoAlert 0, "runEulerCommand() not yet implemented for Windows"
	printf "runEulerCommand() not yet implemented for Windows\r"
	return ""
End
//
#endif
//
Static Function/S readIndexFile(indexFile,path)
	String indexFile						// name of index file, the output from Euler
	String path							// name of Igor path to go with indexFile

	Variable f = OpenFileOf_ftype(indexFile,"IndexFile",path)// open an index file
	if (f<1)
		return ""						// could not get any info from indexFile 
	endif
	FStatus f
	String buffer=""
	buffer = PadString(buffer,V_logEOF-V_filePos,0x20)
	FBinRead f, buffer					// read entire file into buffer
	Close f
	buffer = ReplaceString("\r\n",buffer,"\n")				// ensure that line term is a new-line only
	buffer = ReplaceString("\n\r",buffer,"\n")
	buffer = ReplaceString("\r",buffer,"\n")
	Variable Nbuf = strlen(buffer)

	STRUCT microGeometry geo
	if (FillGeometryStructDefault(geo))					//fill the geometry structure with test values
		DoAlert 0, "no geometry structure found, did you forget to set it?"
		return GetWavesDataFolder(FullPeakIndexed,2)		// can not compute (px,py), but still a valid result
	endif

	Variable ipattern=0										// number of current patterns
	Variable i1,i0 = strsearch(buffer,"$pattern"+num2istr(ipattern),0)		// put into "key=value" everything up pattern ipattern
	i0= i0<0 ? Inf : i0-1
	String wnote = keyStrFromBuffer(buffer[0,i0])
	Variable Npatterns = NumberByKey("NpatternsFound",wnote,"=")
	if (!(Npatterns>0))
			return ""
	endif
	Variable startx,groupx, starty,groupy, ddLocal
	startx = NumberByKey("startx",wnote,"=")
	groupx = NumberByKey("groupx",wnote,"=")
	starty = NumberByKey("starty",wnote,"=")
	groupy = NumberByKey("groupy",wnote,"=")
	startx = numtype(startx) ? 1 : startx
	groupx = numtype(groupx) ? 1 : groupx
	starty = numtype(starty) ? 1 : starty
	groupy = numtype(groupy) ? 1 : groupy
	ddLocal = geo.ddOffset + NumberByKey("CCDy",wnote,"=")	// get dd from the image (in wave note) if it is there
	geo.dd = ddLocal>0 ? ddLocal :  geo.dd
	Variable ycLocal = NumberByKey("yc",wnote,"=")				// get ycLocal from the image (in wave note) if it is there
	geo.ycent = ycLocal>0 ? ycLocal :  geo.ycent

	Make/N=3/O/D runEulerCommand_qBL, runEulerCommand_qW
	Wave qBL=runEulerCommand_qBL,  qW=runEulerCommand_qW
	Make/N=3/O/D runEulerCommand_qhat
	Wave qvec = runEulerCommand_qhat
	String list1
	Variable/C pz
	Variable qx,qy,qz,hh,kk,ll,int,keV,err
	Variable Ni=-1, cols,i,val, nLines
	do
		i0 +=1												// range of header info for this pattern
		i1 = strsearch(buffer,"\n$array"+num2istr(ipattern),i0)
		list1 = keyStrFromBuffer(buffer[i0,i1])
		list1 = RemoveByKey("pattern"+num2istr(ipattern),list1,"=")
		wnote = MergeKeywordLists(wnote,list1,0,"=","")

		// buffer is now the data (and the leading $arrayN line)
		i0 = i1+1
		i1 = strsearch(buffer,"$pattern"+num2istr(ipattern+1),i0)
		i1 = (i1<1 ? Nbuf : i1)-1
		sscanf buffer[i0,i1],"$array%d %d %d",i,cols,nLines
		if (ipattern==0)									// for first pattern, make the output wave
			Ni = nLines
			Make/N=(Ni,cols+1,Npatterns)/O FullPeakIndexed// make is long enough to also hold the pixel positions (an extra 2 cols)
			FullPeakIndexed = NaN
		elseif (nLines>Ni)
			Ni = nLines
			Redimension/N=(Ni,cols+1,Npatterns) FullPeakIndexed
		endif

		i0 = strsearch(buffer,"\n",i0+1)+1				// start of first data
		for (i=0;i<nLines;i+=1)
			sscanf buffer[i0,i1], "    [ %d]   (%g  %g %g)     ( %d %d  %d)   %g,    %g,    %g",val,qx,qy,qz,hh,kk,ll,int,keV,err
			qBL[0]=qx ; qBL[1]=qy ; qBL[2]=qz
			IdealBeamLine2wenge(geo,qBL,qW)
			FullPeakIndexed[i][0][ipattern] = qW[0]		// transform Ideal Beam-line to Wenge
			FullPeakIndexed[i][1][ipattern] = qW[1]
			FullPeakIndexed[i][2][ipattern] = qW[2]
			FullPeakIndexed[i][3][ipattern] = hh
			FullPeakIndexed[i][4][ipattern] = kk
			FullPeakIndexed[i][5][ipattern] = ll
			FullPeakIndexed[i][6][ipattern] = int
			FullPeakIndexed[i][7][ipattern] = keV
			FullPeakIndexed[i][8][ipattern] = err
			i0 = strsearch(buffer,"\n",i0+1)+1			// start of next line
		endfor
		for (i=0;i<nLines;i+=1)							// convert the indexed qhats to pixel positions
		qvec = FullPeakIndexed[i][p][ipattern]
			pz = q2pixel(geo,qvec)
			FullPeakIndexed[i][9][ipattern] = (real(pz)-(startx-1)-(groupx-1)/2)/groupx		// change to binned pixels
			FullPeakIndexed[i][10][ipattern] = (imag(pz)-(starty-1)-(groupy-1)/2)/groupy		// pixels are still zero based
		endfor
		ipattern += 1
	while (ipattern<Npatterns && i0<(Nbuf-1))
	for (ipattern=0;ipattern<Npatterns;ipattern+=1)		// remove all (000) for all patterns
		for (i=0;i<Ni;i+=1)
			if (FullPeakIndexed[i][0][ipattern]==0 && FullPeakIndexed[i][1][ipattern]==0 && FullPeakIndexed[i][2][ipattern]==0)
				FullPeakIndexed[i][][ipattern] = NaN
			endif
		endfor
	endfor
	wnote = ReplaceStringByKey("waveClass",wnote,"IndexedPeakList","=")

	// a legacy section
	Variable SpaceGroup = NumberByKey("SpaceGroup",wnote,"=")
	if (numtype(SpaceGroup))					// could not find 'SpaceGroup', see if it is still called 'latticeStructure'
		SpaceGroup = NumberByKey("latticeStructure",wnote,"=")
		if (numtype(SpaceGroup)==0)
			wnote = ReplaceNumberByKey("SpaceGroup",wnote,spaceGroup,"=")
		endif
	endif

	Note/K FullPeakIndexed, wnote
	SetDimLabel 1,0,Qx,FullPeakIndexed		;	SetDimLabel 1,1,Qy,FullPeakIndexed
	SetDimLabel 1,2,Qz,FullPeakIndexed		;	SetDimLabel 1,3,h,FullPeakIndexed
	SetDimLabel 1,4,k,FullPeakIndexed			;	SetDimLabel 1,5,l,FullPeakIndexed
	SetDimLabel 1,6,Intensity,FullPeakIndexed	;	SetDimLabel 1,7,keV,FullPeakIndexed
	SetDimLabel 1,8,angleErr,FullPeakIndexed	 ;	SetDimLabel 1,9,pixelX,FullPeakIndexed
	SetDimLabel 1,10,pixelY,FullPeakIndexed	

	KillWaves/Z runEulerCommand_qhat, runEulerCommand_qBL, runEulerCommand_qW
	KillWaves/Z M_Lower,M_Upper,W_LUPermutation,M_x
	return GetWavesDataFolder(FullPeakIndexed,2)
End
//
Static Function OpenFileOf_ftype(fname,ftype,path)// open a file only if it is of type ftype
	String fname									// full path name to file with tagged geometry values
	String ftype										// the required file identifier, included as a tag (ftype is optional)
	String path										// name of Igor path

	Variable OK = 0
	Variable f = 0


	if (strlen(ftype)<1)
		OK = 1
	else
		String list = getListOfTypesInFile(fname,path)
		OK = (WhichListItem(ftype,list)>=0)
	endif

	if (OK)
		Open/M="file with a keyword list"/P=$path/R/Z=2 f as fname
		if (strlen(S_fileName)<1 || !f)
			if (f)
				Close f
			endif
			return 0
		endif
	endif
	return f
End
//
Static Function/S keyStrFromBuffer(buffer)// get key list from the buffer
	String buffer 							// long sting with line feeds, each line is of the form:   "$tag   value   // comment\n"

	Variable N=strlen(buffer)
	Variable i0,i1							// indicies into buffer
	String tagName, value					// the found tag and value (if they exist)
	String line								// one line of buffer
	String list=""							// the result
	buffer = ReplaceString("\r\n",buffer,"\n")	// ensure that line term is a new-line only
	buffer = ReplaceString("\n\r",buffer,"\n")
	buffer = ReplaceString("\r",buffer,"\n")
	i0 = 0
	i1 = strsearch(buffer,"\n",i0)
	do
		line = buffer[i0,i1]						// one line from the buffer
		tagValueFromLine(line,tagName,value)	// get tag and value from this line
		if (strlen(tagName) && !keyInList(tagName,list,"=","")) // valid tag, and not yet in list
			list = ReplaceStringByKey(tagName,list,value,"=")
		endif
		i0 = i1+1
		i1 = strsearch(buffer,"\n",i0)
		i1 = (i1<0) ? N-1 : i1
	while(i0<N)
	return list
End
//
Static Function/S tagValueFromLine(lineIn,tagName,value)// find the tag and value from a line
	String lineIn								// line of text to process
	String &tagName							// the found tag (if it exists)
	String &value								// the found value for that tag (if it exists)

	String line=lineIn
	tagName = ""
	value = ""
	Variable i,dollar = char2num("$")
	if (char2num(line)!=dollar)				// if it starts with a $, probable tag so process it
		return ""
	endif

	i = strsearch(line,"//",0)					// strip off comments
	if (i>=0)
		line = line[0,i-1]
	endif

	for (i=0;char2num(line[i+1])>32;i+=1)	// find end of tag, it ends with a space or lower
	endfor
	tagName = line[1,i]
	if (strlen(tagName)<1)						// no tag
		return ""
	endif

	for (i=i+1;char2num(line[i])<=32;i+=1)	// find first non-white space, start of value
	endfor
	value = line[i,Inf]							// value associated with tagName

	for (i=strlen(value)-1;char2num(value[i])<=32 && i>0;i-=1)	// strip off trailing whitespace
	endfor
	value = ReplaceString(";",value[0,i],":")	// cannot have a semicolon or equals in the value
	return ReplaceStringByKey(tagName,"",value,"=")
End
//Function test()
//	Variable f
//	Open/M="test file of tagged lines"/P=home/R/T="TEXT" f as "Macintosh HD:Users:tischler:data:Sector 34:August 2006:calibration:generic_Index.txt"
//	print S_fileName
//	FStatus f
//	print "V_logEOF=",V_logEOF
//	String buffer=""
//	buffer = PadString(buffer,V_logEOF,0x20)
//	FBinRead f, buffer
//	Close f
//	print "strlen(buffer) = ",strlen(buffer)
//	String list = keyStrFromBuffer(buffer)
//	print list
//End
//
//Function test()
//	String tagName,value,line
//	line = "$EulerAngles {-171.29728810, 149.32171846, 229.41456865}	// Euler angles for this pattern (deg)\n"
//	String output = tagValueFromLine(line,tagName,value)
//	printf "tagName = '%s'\r",tagName
//	printf "value = '%s'\r",value
//	printf "output = '%s'\r",output
//End


// remove bkg from an image, identify the peaks, and fit them all
Function/S RemoveBackgroundImage(rawImage,fractionBkg)
	Wave rawImage				// the raw image from the spe or HDF file
	Variable fractionBkg			// fraction of image that is background, use something like 0.99 or 0.995 (even 0.7 works OK)

	Variable printIt = ItemsInList(GetRTStackInfo(0))<2 || stringmatch(StringFromList(0,GetRTStackInfo(0)),"IndexButtonProc")
	if (!WaveExists(rawImage) || WaveDims(rawImage)!=2 || !(fractionBkg==limit(fractionBkg,1e-3,1)))
		String imageName = SelectString(WaveExists(rawImage),"",NameOfWave(rawImage))
		fractionBkg = fractionBkg==limit(fractionBkg,1e-3,1) ? fractionBkg : 0.8
		Prompt imageName,"name of image with peaks to find",popup,reverseList(WaveListClass("speImage","*","DIMS:2"))
		Prompt fractionBkg,"fraction of image that is background, usually in range [0.7, 0.995]"
		DoPrompt/Help="3D-Xray Diffraction[Fit Peaks]" "peak fitting",imageName,fractionBkg
		if (V_flag)
			return ""
		endif
		Wave rawImage = $imageName
		printf "FitPeaks(%s,%g)\r",imageName,fractionBkg
		printIt = 1
	endif
	if (!WaveExists(rawImage) || WaveDims(rawImage)!=2 || !(fractionBkg==limit(fractionBkg,1e-3,1)))
		return ""
	endif
	Variable timer=startMSTimer
//	Variable timer1 = startMSTimer
	Variable thresh = numpnts(rawImage)*fractionBkg

	WaveStats/Q/M=1 rawImage
	Make/N=200/D/O FitPeaks_Hist
	Variable i
	for (i=-1;i<10 || numtype(i);)							// keep expanding the histogram range while i<10,   I want i to be bigger than 10
		SetScale/I x,V_min,V_max,"",FitPeaks_Hist			// lower high end of histogram to zoom in around zero
		Histogram/B=2 rawImage,FitPeaks_Hist
		Integrate/P FitPeaks_Hist/D=FitPeaks_Hist_INT	// integrate the histogram
		i=BinarySearchInterp(FitPeaks_Hist_INT,thresh)	// find point where FitPeaks_Hist_INT[i] == thresh
		V_max = numtype(i) ? V_max/2 : pnt2x(FitPeaks_Hist_INT,i+4)
	endfor
//	print "first part took ",stopMSTimer(timer1)/1e6,"sec"
//	timer1 = startMSTimer

	thresh = DimDelta(FitPeaks_Hist,0)*i + DimOffset(FitPeaks_Hist,0)
	ImageThreshold/Q/T=(thresh)/M=0 rawImage					// create M_ImageThresh
	Wave M_ImageThresh=M_ImageThresh

	ImageStats/M=1 M_ImageThresh
	Variable fracBlack = V_avg*V_npnts/255/V_npnts
	//	printf "number of points in peak is %d out of %d points (%.2f%%)\r",fracBlack*V_npnts,V_npnts,fracBlack*100
	Make/N=(9,9)/O/B/U RemoveBackgroundImage_black
	RemoveBackgroundImage_black = 255
	for(;fracBlack<0.10;)											// loop until fracBlack is atleast 10%
		Imagemorphology/S=RemoveBackgroundImage_black/O BinaryDilation M_ImageThresh
		ImageThreshold/Q/T=(2)/M=0 M_ImageThresh
		ImageStats/M=1 M_ImageThresh
		fracBlack = V_avg*V_npnts/255/V_npnts
		//	printf "number of points in peak is %d out of %d points (%.2f%%)\r",fracBlack*V_npnts,V_npnts,fracBlack*100
	endfor
	ImageTransform/O invert M_ImageThresh						//	M_ImageThresh = !M_ImageThresh

//	print "second part took ",stopMSTimer(timer1)/1e6,"sec"
//	timer1 = startMSTimer

	ImageInterpolate /f={.125,.125} bilinear rawImage
	Redimension/S M_InterpolatedImage

	KillWaves/Z FitPeaks_rawImage_smaller
	Rename M_InterpolatedImage FitPeaks_rawImage_smaller		// make a source image using fewer pixels
//	print "\t\t1st interpolate part took ",stopMSTimer(timer1)/1e6,"sec"

//	timer1 = startMSTimer
 	ImageInterpolate /f={.125,.125} bilinear M_ImageThresh		// make a mask with fewer pixels
	KillWaves/Z FitPeaks_mask_smaller
	Rename M_InterpolatedImage FitPeaks_mask_smaller
	Redimension/B/U FitPeaks_mask_smaller
	ImageRemoveBackground/F/R=FitPeaks_mask_smaller/P=2/W FitPeaks_rawImage_smaller	// creates M_RemovedBackground, uses 3rd order polynomial
	Wave W_BackgroundCoeff=W_BackgroundCoeff
	//	print W_BackgroundCoeff
//	print "\t\t2nd interpolate part took ",stopMSTimer(timer1)/1e6,"sec"

//	timer1 = startMSTimer
	ImageInterpolate /f={8,8} bilinear M_RemovedBackground		// creates M_InterpolatedImage, a bkg with full number of pixels
	Wave M_InterpolatedImage=M_InterpolatedImage

	String outName = NameOfWave(rawImage)+"NoBkg"
	String fldr = GetWavesDataFolder(rawImage,1)
	if (CheckName(fldr+outName,1))
		outName = CleanupName(outName,0)
	endif
	outName = fldr+outName
	Duplicate/O rawImage $outName
	Wave noBkg = $outName

	String wnote=ReplaceStringByKey("rawIgorImage",note(rawImage),GetWavesDataFolder(rawImage,2),"=")
	wnote = ReplaceNumberByKey("fractionBkg",wnote,fractionBkg,"=")
	wnote = ReplaceStringByKey("waveClass",wnote,"speImageNoBkg","=")
	Note/K noBkg, wnote
	Redimension/I noBkg											// signed int32 version of rawImage
	noBkg = rawImage - M_InterpolatedImage[p][q]					// p,q because M_InterpolatedImage is smaller

	Variable sec0=stopMSTimer(timer)/1e6
	if (printIt)
		printf "removing background  took %g sec\r",sec0
		printf "  result stored in the wave '%s'\r", GetWavesDataFolder(noBkg,2)
	endif

	KillWaves/Z FitPeaks_rawImage_smaller,FitPeaks_mask_smaller, M_ImageThresh, M_RemovedBackground, M_InterpolatedImage, W_BackgroundCoeff
	KillWaves/Z FitPeaks_Hist_INT,FitPeaks_Hist, RemoveBackgroundImage_black
	return GetWavesDataFolder(noBkg,2)
End





// First get a background:
//		1) down sample the image by 2 (could be 4)
//		2) do a median filter to get a background level for each pixel
//		3) gaussian filter the background to smooth it
//		4) upsample to get a background the correct size
//	Subtract this background from the original image
//	and median filter (only 5) and gauss smooth to get rid of specks
//	Thresold this image to generate a mask of regions that represent the peak locations (an ROI)
//	Examine each region in the ROI separately, fitting a gaussian
Function/S FitPeaksWithSeedFill(image,minPeakWidth,maxPeakWidth,minSep,threshAboveAvg,[mask,maxNu,whoami])
	Wave image						// the image from the spe or HDF file
	Variable minPeakWidth			// minimum width for an acceptable peak (pixels)
	Variable maxPeakWidth			// maximum width for an acceptable peak (pixels)
	Variable minSep					// minimum separation between two peaks (pixels)
	Variable threshAboveAvg			// threshold above average value,  use 25
	Wave mask						// starting mask for the fit (this mask is unchanged by this routine)
	Variable maxNu					// maximum number of peaks to find, normally goes to completion
	Variable whoami					// if passed, then just return the name of this routine
	if (!ParamIsDefault(whoami))
		return GetRTStackInfo(1)
	endif
	maxNu = ParamIsDefault(maxNu) ? Inf : maxNu

	String imageName
	Variable printIt = ItemsInList(GetRTStackInfo(0))<2 || stringmatch(StringFromList(0,GetRTStackInfo(0)),"IndexButtonProc")
	if (!WaveExists(image) || WaveDims(image)!=2 || !(minPeakWidth>0) || !(maxPeakWidth>0) || !(minSep>3) || !(threshAboveAvg>0))
		imageName = SelectString(WaveExists(image),"",NameOfWave(image))
		minPeakWidth = minPeakWidth>.5 ? minPeakWidth : 1.2
		maxPeakWidth = maxPeakWidth>.5 ? maxPeakWidth : 70
		minSep = minSep>minPeakWidth ? minSep : 50
		threshAboveAvg = threshAboveAvg>0 ? threshAboveAvg : NaN
		Prompt imageName,"name of image with peaks to find",popup,reverseList(WaveListClass("speImageNoBkg;speImage","*","DIMS:2"))
		Prompt minPeakWidth,"minimum allowed width of a peak"
		Prompt maxPeakWidth,"maximum allowed width of a peak"
		Prompt minSep,"minimum distance between two peaks, reject peaks closer than this"
		Prompt threshAboveAvg,"min peak height, use NaN for auto detect"
		DoPrompt/Help="3D-Xray Diffraction[Fit Peaks]" "peak fitting",imageName,minPeakWidth,maxPeakWidth,minSep,threshAboveAvg
		if (V_flag)
			return ""
		endif
		Wave image = $imageName
		printf "FitPeaksWithSeedFill(%s,%g,%g,%g,%g)\r",imageName,minPeakWidth,maxPeakWidth,minSep,threshAboveAvg
		printIt = 1
	endif
	if (!WaveExists(image) || WaveDims(image)!=2)
		return ""
	endif
	if (!(threshAboveAvg>0))
		if (strlen(StringByKey("depthSi",note(image),"=")))
			// threshAboveAvg = 10
			threshAboveAvg = noiseInImage(image)/2
		else
			threshAboveAvg = StatsMedian(image)/4 
			// !(threshAboveAvg>0) = !(threshAboveAvg>0) ? StatsMedian(image)/4 : threshAboveAvg
		endif
		if (printIt)
			printf "   using a threshAboveAvg = %.3g\r",threshAboveAvg
		endif
	endif
	if (!(minPeakWidth>0) || !(maxPeakWidth>0) || !(minSep>3) || !(threshAboveAvg>0))
		return ""
	endif
	imageName = GetWavesDataFolder(image,2)
	Variable Nx=DimSize(image,0), Ny=DimSize(image,1)
	if (Nx<2 || Ny<2)
		return ""
	endif
	Variable timer=startMSTimer

	Variable Nfilt=20
	Variable downSample = NumVarOrDefault("downSample", 2 )						// could be 4 (don't be greedy),  1 also works but is slow

	Duplicate/O image FitPeaksSeedFill_smth											// image to use when finding next hot pixel
	Wave imageS = FitPeaksSeedFill_smth
	Redimension/D imageS

	ImageInterpolate/PXSZ={(downSample),(downSample)} Pixelate imageS			// down sample the image to get background
	Wave M_PixelatedImage=M_PixelatedImage

//	Variable Ndown = Nfilt/downSample
	Variable Ndown = max(7,round(Nfilt/downSample))
	MatrixFilter/N=(Ndown)/M=(Ndown*Ndown*0.5) median M_PixelatedImage
	MatrixFilter/N=(2*Nfilt) gauss M_PixelatedImage								// this is now the background

	ImageInterpolate/F={(downSample),(downSample)} bilinear M_PixelatedImage	// up sample the background
	Wave M_InterpolatedImage=M_InterpolatedImage
	Variable rowsToAdd = DimSize(image,0)-DimSize(M_InterpolatedImage,0)		// redimension to exactly the right size
	Variable colsToAdd = DimSize(image,1)-DimSize(M_InterpolatedImage,1)
	ImageTransform/N={(rowsToAdd),(colsToAdd)}/O padimage M_InterpolatedImage

	ImageStats M_InterpolatedImage
	Variable bkg = V_rms
	MatrixOp/O imageS = image - M_InterpolatedImage			// imageS = image - M_InterpolatedImage
	MatrixFilter/N=5 median imageS
	MatrixFilter/N=10 gauss imageS
	ImageThreshold/T=(threshAboveAvg)/I imageS
	Wave M_ImageThresh=M_ImageThresh
	Duplicate/O M_ImageThresh FitPeaks_ImageMask 
	if (ParamIsDefault(mask))									// preset FitPeaks_ImageMask with 'mask' if it was passed and is valid
		if (WaveExists(mask) && DimSize(mask,0)==Nx && DimSize(mask,1)==Ny)	// if mask exists and is right size, use it as starting point
			FitPeaks_ImageMask = mask || FitPeaks_ImageMask
		endif
	endif
	Duplicate/O imageS FitPeaksWithSeedFill_weights__		// weights to use in the peak fitting
	Wave weights = FitPeaksWithSeedFill_weights__
	Duplicate/O weights FitPeaksWithSeedFill_ones__			// weights to use in the peak fitting
	Wave ones = FitPeaksWithSeedFill_ones__
	Redimension/B/U ones
	ones = 1

	String str
	Make/O/T JZT_Constraints={"K3 > 0.1","K5 > 0.1","K6 > 0","K6 < 1","","","",""}
	sprintf str,"K3 > %g",minPeakWidth/100
	JZT_Constraints[0] = str
	sprintf str,"K5 > %g",minPeakWidth/100
	JZT_Constraints[1] = str
	sprintf str,"K2 > %g",-maxPeakWidth
	JZT_Constraints[4] = str
	sprintf str,"K2 < %g",DimOffset(image,0) + (DimSize(image,0)-1)*DimDelta(image,0) + maxPeakWidth
	JZT_Constraints[5] = str
	sprintf str,"K4 > %g",-maxPeakWidth
	JZT_Constraints[6] = str
	sprintf str,"K4 < %g",DimOffset(image,1) + (DimSize(image,1)-1)*DimDelta(image,1) + maxPeakWidth
	JZT_Constraints[7] = str
//	sprintf str,"K1 > %g",-100*threshAboveAvg
//	JZT_Constraints[8] = str

	Variable Nu=0, Nlen=50
	Make/N=(Nlen,11)/O FullPeakList
	FullPeakList = NaN
	String wnote = note(image)
	wnote = ReplaceNumberByKey("minPeakWidth",wnote,minPeakWidth,"=")
	wnote = ReplaceNumberByKey("maxPeakWidth",wnote,maxPeakWidth,"=")
	wnote = ReplaceNumberByKey("minSpotSeparation",wnote,minSep,"=")
	wnote = ReplaceNumberByKey("threshAboveAvg",wnote,threshAboveAvg,"=")
	wnote = ReplaceStringByKey("fittedIgorImage",wnote,GetWavesDataFolder(image,2),"=")
	wnote = ReplaceStringByKey("waveClass",wnote,"FittedPeakList","=")
	SetDimLabel 1,0,x0,FullPeakList		;	SetDimLabel 1,1,y0,FullPeakList
	SetDimLabel 1,2,x0Err,FullPeakList	;	SetDimLabel 1,3,y0Err,FullPeakList
	SetDimLabel 1,4,fwx,FullPeakList		;	SetDimLabel 1,5,fwy,FullPeakList
	SetDimLabel 1,6,fwxErr,FullPeakList	;	SetDimLabel 1,7,fwyErr,FullPeakList
	SetDimLabel 1,8,correlation,FullPeakList ;	SetDimLabel 1,9,correlationErr,FullPeakList
	SetDimLabel 1,10,area,FullPeakList
	Note/K FullPeakList,wnote

	Variable fw= ( 2*sqrt(2*ln(2)) )							// this factor corrects for 2-d sigma to FWHM, apply to K3, K5, and the corresponding errors
	Variable xx,yy,fwx,fwy,sigX,sigY, amp						// holds the resutls of gaussian fit
	Variable px,py, left,right,top,bot, err, i
	do
		// DoUpdate
		ImageStats/M=1/R=FitPeaks_ImageMask imageS		// find the highest point
		if (V_max<threshAboveAvg || V_npnts<2)
			break
		endif
		px=V_maxRowLoc ; py=V_maxColLoc					// look for a peak at (px,py)
		// fit a gaussian about px,py
		left = limit(floor(px-maxPeakWidth),0,Nx-1)		// box to use for fitting
		right = limit(ceil(px+maxPeakWidth),0,Nx-1)
		top = limit(floor(py-maxPeakWidth),0,Ny-1)
		bot = limit(ceil(py+maxPeakWidth),0,Ny-1)
		ImageSeedFill seedP=px,seedQ=py,min=0,max=0,target=1,srcwave=FitPeaks_ImageMask
		Wave M_SeedFill=M_SeedFill
		MatrixOp/O weights = equal(M_SeedFill,ones)			// weights = M_SeedFill==1 ? 1 : 0
		err = FitGaussianPeak(imageS,left,right,top,bot,weight=weights)	// returns 0 == OK,  non-zero if error
		ImageSeedFill/O seedP=px,seedQ=py,min=0,max=0,target=255,srcwave=FitPeaks_ImageMask	// fil in this region, already used
		if (err)
			continue
		endif
		Wave W_coef=W_coef, W_sigma=W_sigma
		xx = W_coef[2]
		yy = W_coef[4]
		sigX = W_sigma[2]
		sigY = W_sigma[4]
		fwx = abs(W_coef[3])*fw
		fwy = abs(W_coef[5])*fw
		amp = abs(W_coef[1])

//		if (ItemsInList(GetRTStackInfo(0))<2)
//			printf "%d    Gaussian peak in at (%.2f±%.2f, %.2f±%.2f) with FWHM of (%.2f±%.2f, %.2f±%.2f),   amp=%g,  image[x,y]= %g\r",Nu,xx,sigX,yy,sigY,fwx,W_sigma[3]*fw,fwy,W_sigma[5]*fw,amp,image[xx][yy]
//		endif
		err = ( sigX>maxPeakWidth || sigY>maxPeakWidth || numtype(sigX+sigY) )						// errors bars are huge
		if (!(min(fwx,fwy)>minPeakWidth) || !(max(fwx,fwy)<maxPeakWidth) || !(max(amp,image[xx][yy])>bkg) || err)	// cannot have spots too wide or too narrow
			continue
		endif
		for (i=0;i<Nu;i+=1)
			if (((xx-FullPeakList[i][0])^2 + (yy-FullPeakList[i][1])^2) < minsep^2)
				i = -1
				break
			endif
		endfor
		if (i<0)
			continue
		endif

		if (Nu>=Nlen)												// FullPeakList is too short, extend it
			Nlen += 50
			Redimension/N=(Nlen,11) FullPeakList
		endif
		FullPeakList[Nu][0]=xx				;	FullPeakList[Nu][1]=yy
		FullPeakList[Nu][2]=sigX				;	FullPeakList[Nu][3]=sigY
		FullPeakList[Nu][4]=fwx				;	FullPeakList[Nu][5]=fwy
		FullPeakList[Nu][6]=W_sigma[3]*fw	;	FullPeakList[Nu][7]=W_sigma[5]*fw
		FullPeakList[Nu][8]=W_coef[6]		;	FullPeakList[Nu][9]=W_sigma[6]
		FullPeakList[Nu][10]=W_coef[1]*fwx*fwy
		Nu += 1
	while(Nu<maxNu)
	Redimension/N=(Nu,11) FullPeakList						// set to exact length
	ImageStats/M=1/G={0,Nu-1,10,10}/Q FullPeakList
	wnote = ReplaceNumberByKey("totalPeakIntensity",wnote,V_avg*V_npnts,"=")	// total intensity in all fitted peaks
	wnote = ReplaceNumberByKey("totalIntensity",wnote,sum(image),"=")			// total intensity in the image
	Note/K FullPeakList,wnote

	String/G pkLIst=""
	for (i=0;i<Nu;i+=1)											// re-set pkLIst to the fitted peaks
		sprintf str,"%.0f,%.0f;",FullPeakList[i][0],FullPeakList[i][1]
		pkList += str
	endfor

	Variable sec=stopMSTimer(timer)/1e6
	if (printIt)
		printf "found and fitted %d peaks (using a threshold of %.1f),  this all took %g sec\r",Nu,threshAboveAvg,sec
		printf "  result stored in the wave '%s'\r", GetWavesDataFolder(FullPeakList,2)
	endif
	KillWaves/Z M_ImageThresh, M_SeedFill, FitPeaksSeedFill_smth, FitPeaksWithSeedFill_ones__
	KillWaves/Z M_PixelatedImage,M_InterpolatedImage
	KillWaves/Z FitPeaks_ImageMask, JZT_Constraints, FitPeaksWithSeedFill_weights__
	return GetWavesDataFolder(FullPeakList,2)
End

Function noiseInImage(image)
	Wave image
	if (!WaveExists(image))
		return NaN
	endif
	Duplicate/O/I image, diff_noiseInImage
	Wave diff = diff_noiseInImage
	diff = image[p+1][q] - image[p][q]
	WaveStats/Q diff
	KillWaves/Z diff_noiseInImage

	if (ItemsInList(GetRTStackInfo(0))<2)
		print "V_sdev =",V_sdev
	endif
	return V_sdev
End


// identify the peaks, and fit them all (no background removal), a lot like FitPeaksNew(), but proceeds in a more step wise fashion
Function/S FitPeaksStepWise(image,minPeakWidth,maxPeakWidth,minSep,threshAboveAvg,[mask,maxNu,whoami])
	Wave image						// the image from the spe or HDF file
	Variable minPeakWidth			// minimum width for an acceptable peak (pixels)
	Variable maxPeakWidth			// maximum width for an acceptable peak (pixels)
	Variable minSep					// minimum separation between two peaks (pixels)
	Variable threshAboveAvg			// threshold above average value
	Wave mask						// starting mask for the fit (this mask is unchanged by this routine)
	Variable maxNu					// maximum number of peaks to find, normally goes to completion
	Variable whoami					// if passed, then just return the name of this routine
	if (!ParamIsDefault(whoami))
		return GetRTStackInfo(1)
	endif
	maxNu = ParamIsDefault(maxNu) ? Inf : maxNu

	String imageName
	Variable printIt = ItemsInList(GetRTStackInfo(0))<2 || stringmatch(StringFromList(0,GetRTStackInfo(0)),"IndexButtonProc")
	if (!WaveExists(image) || WaveDims(image)!=2 || !(minPeakWidth>0) || !(maxPeakWidth>0) || !(minSep>3) || !(threshAboveAvg>0))
		imageName = SelectString(WaveExists(image),"",NameOfWave(image))
		minPeakWidth = minPeakWidth>.5 ? minPeakWidth : 1.2
		maxPeakWidth = maxPeakWidth>.5 ? maxPeakWidth : 10
		minSep = minSep>minPeakWidth ? minSep : 40
		threshAboveAvg = threshAboveAvg>0 ? threshAboveAvg : 100
		Prompt imageName,"name of image with peaks to find",popup,reverseList(WaveListClass("speImageNoBkg;speImage","*","DIMS:2"))
		Prompt minPeakWidth,"minimum allowed width of a peak"
		Prompt maxPeakWidth,"maximum allowed width of a peak"
		Prompt minSep,"minimum distance between two peaks"
		Prompt threshAboveAvg,"min peak height"
		DoPrompt/Help="3D-Xray Diffraction[Fit Peaks]" "peak fitting",imageName,minPeakWidth,maxPeakWidth,minSep,threshAboveAvg
		if (V_flag)
			return ""
		endif
		Wave image = $imageName
		printf "FitPeaksStepWise(%s,%g,%g,%g,%g)\r",imageName,minPeakWidth,maxPeakWidth,minSep,threshAboveAvg
		printIt = 1
	endif
	if (!WaveExists(image) || WaveDims(image)!=2 || !(minPeakWidth>0) || !(maxPeakWidth>0) || !(minSep>3) || !(threshAboveAvg>0))
		return ""
	endif
	imageName = GetWavesDataFolder(image,2)
	Variable Nx=DimSize(image,0), Ny=DimSize(image,1)
	if (Nx<2 || Ny<2)
		return ""
	endif

	Variable timer0=startMSTimer

	Duplicate/O image FitPeaksStepWise_smth					// image to use when finding next hot pixel
	Wave imageS = FitPeaksStepWise_smth
	Redimension/D imageS
	imageS = image[p][q]<=0 ? NaN : image[p][q]
	MatrixFilter/N=3 min  imageS

	ImageStats/M=1/Q imageS									// find the average value
	Variable threshold = threshAboveAvg+max(V_avg,0)		// the actual threshold, only consider pixels above this level
	if (printIt)
		printf "average for the image = %g,  so only consider pixels values above %g\r",V_avg,threshold
	endif

	Make/N=(Nx,Ny)/O/B/U FitPeaks_ImageMask				// make a mask for this image
	FitPeaks_ImageMask = 0										// set to 0, this includes all of the image
	if (ParamIsDefault(mask))									// preset FitPeaks_ImageMask with 'mask' if it was passed and is valid
		if (WaveExists(mask) && DimSize(mask,0)==Nx && DimSize(mask,1)==Ny)	// if mask exists and is right size, use it as starting point
			FitPeaks_ImageMask = !(!mask[p][q])
		endif
	endif

	Duplicate/O image image_weights_temp_					// weights to use in the peak fitting
	Wave weights = image_weights_temp_
	weights = FitPeaks_ImageMask<1

	Make/N=(3,3)/O/U/I FitPeaks_3x3
	Make/O/T JZT_Constraints={"K3 > 0","K5 > 0","K6 > 0","K6 < 1"}

	Variable Nu=0, Nlen=50
	Make/N=(Nlen,11)/O FullPeakList
	FullPeakList = NaN
	String wnote = note(image)
	wnote = ReplaceNumberByKey("minPeakWidth",wnote,minPeakWidth,"=")
	wnote = ReplaceNumberByKey("maxPeakWidth",wnote,maxPeakWidth,"=")
	wnote = ReplaceNumberByKey("minSpotSeparation",wnote,minSep,"=")
	wnote = ReplaceNumberByKey("threshAboveAvg",wnote,threshAboveAvg,"=")
	wnote = ReplaceStringByKey("fittedIgorImage",wnote,GetWavesDataFolder(image,2),"=")
	wnote = ReplaceStringByKey("waveClass",wnote,"FittedPeakList","=")
	SetDimLabel 1,0,x0,FullPeakList		;	SetDimLabel 1,1,y0,FullPeakList
	SetDimLabel 1,2,x0Err,FullPeakList	;	SetDimLabel 1,3,y0Err,FullPeakList
	SetDimLabel 1,4,fwx,FullPeakList		;	SetDimLabel 1,5,fwy,FullPeakList
	SetDimLabel 1,6,fwxErr,FullPeakList	;	SetDimLabel 1,7,fwyErr,FullPeakList
	SetDimLabel 1,8,correlation,FullPeakList ;	SetDimLabel 1,9,correlationErr,FullPeakList
	SetDimLabel 1,10,area,FullPeakList

	Variable fw= ( 2*sqrt(2*ln(2)) )							// this factor corrects for 2-d sigma to FWHM, apply to K3, K5, and the corresponding errors
	Variable xx,yy,fwx,fwy,sigX,sigY,amp,sigAmp				// holds the resutls of gaussian fit
	Variable px,py, left,right,top,bot
	do
		// DoUpdate
		ImageStats/M=1/R=FitPeaks_ImageMask imageS		// find the highest point
		if (V_max<threshold || V_npnts<2)
			break
		endif
		px=V_maxRowLoc ; py=V_maxColLoc					// look for a peak at (px,py)
		//	if (px==491 && py==157)
		//		printf "point = [%d, %d]\r",px,py
		//	endif

		// check if this is just a lone pixel
		FitPeaks_3x3 = image[p+px-1][q+py-1]
		FitPeaks_3x3[1][1] = 0
		ImageStats/M=1 FitPeaks_3x3
		if (!(V_max>threshold))									// this is a lone pixel
			image[px][py] = V_avg								// set lone pixels to surrounding average, and continue
			FitPeaks_ImageMask[px][py] = 255
			weights[px][py] = 0
			continue
		endif

		// fit a gaussian about px,py
		left = limit(floor(px-maxPeakWidth/2),0,Nx-1)		// box to use for fitting
		right = limit(ceil(px+maxPeakWidth/2),0,Nx-1)
		top = limit(floor(py-maxPeakWidth/2),0,Ny-1)
		bot = limit(ceil(py+maxPeakWidth/2),0,Ny-1)
		if (FitGaussianPeak(image,left,right,top,bot,weight=weights))			// returns 0 == OK,  non-zero if error
			FitPeaks_ImageMask[left,right][top,bot] = 255	// extend the FitPeaks_ImageMask to reject spots within maxPeakWidth
			weights[left,right][top,bot] = 0
			continue
		endif
		Wave W_coef=W_coef, W_sigma=W_sigma
		xx = W_coef[2]
		yy = W_coef[4]
		fwx = abs(W_coef[3] * fw)
		fwy = abs(W_coef[5] * fw)
		sigX = W_sigma[2]
		sigY = W_sigma[4]
		amp = abs(W_coef[1])
		sigAmp = W_sigma[2]
		//	if (ItemsInList(GetRTStackInfo(0))<2)
		//		printf "%d    Gaussian peak in at (%.2f±%.2f, %.2f±%.2f) with FWHM of (%.2f±%.2f, %.2f±%.2f),   amp=%g\r",Nu,xx,sigX,yy,sigY,fwx,W_sigma[3]*fw,fwy,W_sigma[5]*fw,W_coef[1]*fwx*fwy
		//	endif

		if (min(fwx,fwy)<minPeakWidth || max(fwx,fwy)>maxPeakWidth || sigX>maxPeakWidth || sigY>maxPeakWidth || amp<threshAboveAvg || amp/sigAmp<3)	// cannot have spots too wide or too narrow
//			left = limit(floor(px-min(minSep,fwx)),0,Nx-1)// fitted peak is invalid, delete area, use a smaller region than minSep
//			right = limit(ceil(px+min(minSep,fwx)),0,Nx-1)
//			top = limit(floor(py-min(minSep,fwy)),0,Ny-1)
//			bot = limit(ceil(py+min(minSep,fwy)),0,Ny-1)
			left = limit(floor(px-min(maxPeakWidth/2,fwx)),0,Nx-1)// fitted peak is invalid, delete area
			right = limit(ceil(px+min(maxPeakWidth/2,fwx)),0,Nx-1)
			top = limit(floor(py-min(maxPeakWidth/2,fwy)),0,Ny-1)
			bot = limit(ceil(py+min(maxPeakWidth/2,fwy)),0,Ny-1)
			FitPeaks_ImageMask[left,right][top,bot] = 255	// extend the FitPeaks_ImageMask to reject spots within maxPeakWidth
			weights[left,right][top,bot] = 0
			continue
		endif

		// exclusion area around (px,py) for future peaks when this peak is good
//px = round(xx)
//py = round(yy)
		left = limit(px-minSep,0,Nx-1)						// fitted peak was valid, delete area so it does not get counted twice
		right = limit(px+minSep,0,Nx-1)
		top = limit(py-minSep,0,Ny-1)
		bot = limit(py+minSep,0,Ny-1)
		FitPeaks_ImageMask[left,right][top,bot] = 255		// extend the FitPeaks_ImageMask to reject spots within maxPeakWidth
		weights[left,right][top,bot] = 0

		if (Nu>=Nlen)												// FullPeakList is too short, extend it
			Nlen += 50
			Redimension/N=(Nlen,11) FullPeakList
		endif
		FullPeakList[Nu][0]=xx				;	FullPeakList[Nu][1]=yy
		FullPeakList[Nu][2]=sigX				;	FullPeakList[Nu][3]=sigY
		FullPeakList[Nu][4]=fwx				;	FullPeakList[Nu][5]=fwy
		FullPeakList[Nu][6]=W_sigma[3]*fw	;	FullPeakList[Nu][7]=W_sigma[5]*fw
		FullPeakList[Nu][8]=W_coef[6]		;	FullPeakList[Nu][9]=W_sigma[6]
		FullPeakList[Nu][10]=W_coef[1]*fwx*fwy
		Nu += 1
	while(Nu<maxNu)
	Redimension/N=(Nu,11) FullPeakList						// set to exact length
	ImageStats/M=1/G={0,Nu-1,10,10}/Q FullPeakList
	wnote = ReplaceNumberByKey("totalPeakIntensity",wnote,V_avg*V_npnts,"=")	// total intensity in all fitted peaks
	wnote = ReplaceNumberByKey("totalIntensity",wnote,sum(image),"=")			// total intensity in the image
	Note/K FullPeakList,wnote

	String/G pkLIst=""
	String str
	Variable i
	for (i=0;i<Nu;i+=1)											// re-set pkLIst to the fitted peaks
		sprintf str,"%.0f,%.0f;",FullPeakList[i][0],FullPeakList[i][1]
		pkList += str
	endfor

	Variable sec0=stopMSTimer(timer0)/1e6
	if (printIt)
		printf "found and fitted %d peaks,  this all took %g sec\r",Nu,sec0
		printf "  result stored in the wave '%s'\r", GetWavesDataFolder(FullPeakList,2)
	endif
	KillWaves/Z FitPeaks_ImageMask, FitPeaks_3x3, image_weights_temp_
	KillWaves/Z FitPeaksStepWise_smth, JZT_Constraints
	return GetWavesDataFolder(FullPeakList,2)
End
//
Function/S FitPeaksProto(image,minPeakWidth,maxPeakWidth,minSep,threshAboveAvg,[mask,maxNu,whoami])
	Wave image						// the image from the spe or HDF file
	Variable minPeakWidth			// minimum width for an acceptable peak (pixels)
	Variable maxPeakWidth			// maximum width for an acceptable peak (pixels)
	Variable minSep					// minimum separation between two peaks (pixels)
	Variable threshAboveAvg			// threshold above average value,  use 25
	Wave mask						// starting mask for the fit (this mask is unchanged by this routine)
	Variable maxNu					// maximum number of peaks to find, normally goes to completion
	Variable whoami					// if passed, then just return the name of this routine
	return ""
End
// 
Static Function FitGaussianPeak(image,left,right,bottom,top,[weight])	// returns 1 if error,  0 is OK
	// the result is passed to calling function by W_coef or {K0, K1, ... K6}
	// this fits a 2-d gaussina to a region on an image.  The region is selected with a marquee in an image plot, or specified
	Wave image
	Variable left, right, bottom, top								// range to use for fit (defaults to full image for bad values)
	Wave weight
	Variable noWeight = ParamIsDefault(weight)

	left = left==limit(left,0,DimSize(image,0)-1) ? left : 0	// for bad values, use whole image
	bottom = bottom==limit(bottom,0,DimSize(image,1)-1) ? bottom : 0
	right = right==limit(right,0,DimSize(image,0)-1) ? right : DimSize(image,0)-1
	top = top==limit(top,0,DimSize(image,1)-1) ? top : DimSize(image,1)-1

	// set up for the gaussian fit
	Wave/T JZT_Constraints=JZT_Constraints					// = {"K3 > 0","K5 > 0","K6 > 0","K6 < 1"}
	Variable V_FitOptions=4, V_FitError=0,V_FitQuitReason=0
	V_FitOptions=4

	if (!noWeight)												// limit box to non-zero weights
		Variable l,r,t,b
		maskBoundingBox(weight,l,r,t,b)
		left = max(left,l)
		right = min(right,r)
		top = min(top,b)
		bottom = max(bottom,t)
	endif
	if ((right-left)<3 || (top-bottom)<3)
		return 1													// need bigger region for a 2d fit
	endif

	if (noWeight)
		CurveFit/Q/N Gauss2D image[left,right][bottom,top]/C=JZT_Constraints
	else
		CurveFit/Q/N Gauss2D image[left,right][bottom,top]/C=JZT_Constraints/W=weight
	endif

	if (!V_FitError)												// if no error, but constraints failed, set bit 5 as error anyhow
		Wave W_coef = W_coef
		V_FitError += (W_coef[3]<-1e-7 || W_coef[5]<-1e-7 || W_coef[6]<-1e-7 || W_coef[6]>1.0001) ? 32 : 0
	endif
	return V_FitError
End
//
Static Function maskBoundingBox(mask,left,right,top,bot)
	Variable &left,&right,&top,&bot
	Wave mask
	left = Inf
	right = -Inf
	top = Inf
	bot = -Inf

	Variable i,Nx=DimSize(mask,0), Ny=DimSize(mask,1)
	for (i=0;i<Nx;i+=1)
		ImageStats/M=1/G={i,i,0,Ny-1} mask
		if (V_max)
			left = min(left,V_maxRowLoc)
			right = max(right,V_maxRowLoc)
		endif
	endfor
	for (i=0;i<Ny;i+=1)
		ImageStats/M=1/G={0,Nx-1,i,i} mask
		if (V_max)
			top = min(top,V_maxColLoc)
			bot = max(bot,V_maxColLoc)
		endif
	endfor
End


//Static Function FitGaussianPeak(image,left,right,bottom,top,[weight])	// returns 1 if error,  0 is OK
//	// the result is passed to calling function by W_coef or {K0, K1, ... K6}
//	// this fits a 2-d gaussina to a region on an image.  The region is selected with a marquee in an image plot, or specified
//	Wave image
//	Variable left, right, bottom, top								// range to use for fit
//	Wave weight
//	Variable noWeight = ParamIsDefault(weight)
//
//	// set up for the gaussian fit
//	Wave JZT_Constraints=JZT_Constraints					// = {"K3 > 0","K5 > 0","K6 > 0","K6 < 1"}
//	Variable V_FitOptions=4, V_FitError=0,V_FitQuitReason=0
//	V_FitOptions=4
//	if (noWeight)
//		CurveFit/Q/N Gauss2D image[left,right][bottom,top]/C=JZT_Constraints
//	else
//		CurveFit/Q/N Gauss2D image[left,right][bottom,top]/C=JZT_Constraints/W=weight
//	endif
//
//	if (!V_FitError)												// if no error, but constraints failed, set bit 5 as error anyhow
//		Wave W_coef = W_coef
//		V_FitError += (W_coef[3]<-1e-7 || W_coef[5]<-1e-7 || W_coef[6]<-1e-7 || W_coef[6]>1.0001) ? 32 : 0
//	endif
//	return V_FitError
//End

// identify the peaks, and fit them all (no background removal)
Function/S FitPeaksNew(image,threshAboveAvg,dist)
	Wave image						// the image from the spe or HDF file
	Variable threshAboveAvg		// threshold above average value
	Variable dist					// minimum distance between spots, spots have to be at least this far apart (pixels)

	Variable printIt = ItemsInList(GetRTStackInfo(0))<2 || stringmatch(StringFromList(0,GetRTStackInfo(0)),"IndexButtonProc")
	if (!WaveExists(image) || WaveDims(image)!=2 || !(threshAboveAvg>0) || !(dist>3))
		String imageName = SelectString(WaveExists(image),"",NameOfWave(image))
		threshAboveAvg = threshAboveAvg>0 ? threshAboveAvg : 100
		dist = dist>3 ? dist : 30
		Prompt imageName,"name of image with peaks to find",popup,reverseList(WaveListClass("speImageNoBkg;speImage","*","DIMS:2"))
		Prompt threshAboveAvg,"fraction of image that is background, usually in range [0.7, 0.995], use 0 to skip bkg removal"
		Prompt dist, "minimum distance (pixels) between any two peaks (use ~30 for unbinned Si)"
		DoPrompt/Help="3D-Xray Diffraction[Fit Peaks]" "peak fitting",imageName,threshAboveAvg,dist
		if (V_flag)
			return ""
		endif
		Wave image = $imageName
		printf "FitPeaks(%s,%g,%g)\r",imageName,threshAboveAvg,dist
		printIt = 1
	endif
	if (!WaveExists(image) || WaveDims(image)!=2 || !(threshAboveAvg>0) || !(dist>3))
		return ""
	endif
	Variable Nx=DimSize(image,0), Ny=DimSize(image,1)

	Variable timer0=startMSTimer
	ImageStats/M=1/Q image									// find the average value
	Variable threshold = threshAboveAvg+V_avg					// the actual threshold, only consider pixels above this level
if (printIt)
print "average for the image = ",V_avg
endif
	Duplicate/O image, FitPeaks_ImageMask
	Wave FitPeaks_ImageMask=FitPeaks_ImageMask
	Redimension/B/U FitPeaks_ImageMask						// make a mask for the image
	FitPeaks_ImageMask = 0									// set to 0, this includes all of the image
	Variable localAvg											// average of local spots (not including center)

	Make/N=(3,3)/O/U/I FitPeaks_3x3
	Variable xx,yy,fwx,fwy,sigX,sigY
	Variable fw= ( 2*sqrt(2*ln(2)) )							// this factor corrects for 2-d sigma to FWHM, apply to K3, K5, and the corresponding errors

	Make/N=(50,11)/O FullPeakList
	FullPeakList = NaN
	String wnote = note(image)
	wnote = ReplaceNumberByKey("minSpotSeparation",wnote,dist,"=")
	wnote = ReplaceNumberByKey("threshAboveAvg",wnote,threshAboveAvg,"=")
	wnote = ReplaceStringByKey("fittedIgorImage",wnote,GetWavesDataFolder(image,2),"=")
	wnote = ReplaceStringByKey("waveClass",wnote,"FittedPeakList","=")
	Note/K FullPeakList,wnote
	SetDimLabel 1,0,x0,FullPeakList		;	SetDimLabel 1,1,y0,FullPeakList
	SetDimLabel 1,2,x0Err,FullPeakList	;	SetDimLabel 1,3,y0Err,FullPeakList
	SetDimLabel 1,4,fwx,FullPeakList		;	SetDimLabel 1,5,fwy,FullPeakList
	SetDimLabel 1,6,fwxErr,FullPeakList	;	SetDimLabel 1,7,fwyErr,FullPeakList
	SetDimLabel 1,8,correlation,FullPeakList ;	SetDimLabel 1,9,correlationErr,FullPeakList
	SetDimLabel 1,10,area,FullPeakList

	Variable i
	Variable px,py, left,right,top,bot
	Variable Nu=0, Nlen=50
	do
		ImageStats/M=1/R=FitPeaks_ImageMask image			// find the highest point
		if (V_max<threshold || V_npnts<2)
			break
		endif
		px = V_maxRowLoc
		py = V_maxColLoc

		// check if this is just a lone pixel
		FitPeaks_3x3 = image[p+px][q+py]
		FitPeaks_3x3[1][1] = 0
		ImageStats/M=1 FitPeaks_3x3
		localAvg = V_avg
		if (!(V_max>threshold))								// this is a lone pixel
			image[px][py] = localAvg							// set lone pixels to local average, and continue
			FitPeaks_ImageMask[px][py] = 255
			continue
		endif

		// fit a gaussian about px,py
		left = limit(px-dist/2,0,Nx-1)
		right = limit(px+dist/2,0,Nx-1)
		top = limit(py-dist/2,0,Ny-1)
		bot = limit(py+dist/2,0,Ny-1)
		if (FitOneGaussianPeak(image,left,right,top,bot))	// returns 0 == OK,  non-zero if error
			FitPeaks_ImageMask[left,right][top,bot] = 255		// extend the FitPeaks_ImageMask to include this spot
			continue
		endif
		Wave W_coef=W_coef, W_sigma=W_sigma
		xx = W_coef[2]
		yy = W_coef[4]
		fwx = W_coef[3] * fw
		fwy = W_coef[5] * fw
		sigX = W_sigma[2]
		sigY = W_sigma[4]
		if (fwx>dist || fwy>dist || fwx<1 || fwy<1 || sigX>dist || sigY>dist)	// cannot have spots wider than dist
			FitPeaks_ImageMask[left,right][top,bot] = 255		// extend the FitPeaks_ImageMask to include this spot
			continue
		endif
		// printf "%d  %d  Gaussian peak in at (%.2f±%.2f, %.2f±%.2f) with FWHM of (%.2f±%.2f, %.2f±%.2f),   amp=%g\r",i,Nu,xx,sigX,yy,sigY,fwx,W_sigma[3]*fw,fwy,W_sigma[5]*fw,W_coef[1]*fwx*fwy
		// check here that the found xx,yy is not too close to existing peaks (if if is stamp out with roi again)
		for (i=0;i<Nu;i+=1)									// check that this is not too close to an existing peak
			if (sqrt((FullPeakList[i][0]-xx)^2 + (FullPeakList[i][1]-yy)^2)<dist)
				break
			endif
		endfor
		if (i<Nu)												// there is already a spot too close to here
			FitPeaks_ImageMask[left,right][top,bot] = 255		// extend the FitPeaks_ImageMask to include this spot
			continue
		endif
		if (Nu>=Nlen)											// FullPeakList is too short, extend it
			Nlen += 50
			Redimension/N=(Nlen,11) FullPeakList
		endif
if (printIt)
printf "found Gaussian start @ (%d, %d)=%d   with left=%.1f,  right=%.1f,  top=%.1f,  bot=%.1f      result=(%.3f, %.3f)\r",px,py,image[px][py],left,right,top,bot,xx,yy
endif
		FullPeakList[Nu][0]=xx				;	FullPeakList[Nu][1]=yy
		FullPeakList[Nu][2]=sigX				;	FullPeakList[Nu][3]=sigY
		FullPeakList[Nu][4]=fwx				;	FullPeakList[Nu][5]=fwy
		FullPeakList[Nu][6]=W_sigma[3]*fw	;	FullPeakList[Nu][7]=W_sigma[5]*fw
		FullPeakList[Nu][8]=W_coef[6]		;	FullPeakList[Nu][9]=W_sigma[6]
		FullPeakList[Nu][10]=W_coef[1]*fwx*fwy
		Nu += 1
		FitPeaks_ImageMask[left,right][top,bot] = 255			// extend the FitPeaks_ImageMask to include this spot
	while(1)
	Redimension/N=(Nu,11) FullPeakList

	String/G pkLIst=""
	String str
	for (i=0;i<Nu;i+=1)
		sprintf str,"%.0f,%.0f;",FullPeakList[i][0],FullPeakList[i][1]
		pkList += str
	endfor

	Variable sec0=stopMSTimer(timer0)/1e6
	if (printIt)
		printf "found and fitted %d peaks,  this all took %g sec\r",i,sec0
		printf "  result stored in the wave '%s'\r", GetWavesDataFolder(FullPeakList,2)
	endif
	KillWaves/Z FitPeaks_ImageMask, FitPeaks_3x3
	return GetWavesDataFolder(FullPeakList,2)
End


// remove bkg from an image, identify the peaks, and fit them all
Function/S FitPeaks(image,fractionBkg,dist)
	Wave image					// the image from the spe or HDF file (with or without bkg removed)
	Variable fractionBkg			// fraction of image that is background, use something like 0.99 or 0.995 (even 0.7 works OK)
	Variable dist					// minimum distance between spots, spots have to be at least this far apart (pixels)

	Variable printIt = ItemsInList(GetRTStackInfo(0))<2 || stringmatch(StringFromList(0,GetRTStackInfo(0)),"IndexButtonProc")
	if (!WaveExists(image) || WaveDims(image)!=2 || !(fractionBkg==limit(fractionBkg,1e-3,1) || !(dist>0)))
		String imageName = SelectString(WaveExists(image),"",NameOfWave(image))
		fractionBkg = fractionBkg==limit(fractionBkg,1e-3,1) ? fractionBkg : 0.9
		dist = dist>0 ? dist : 30
		Prompt imageName,"name of image with peaks to find",popup,reverseList(WaveListClass("speImage*","*","DIMS:2"))
		Prompt fractionBkg,"fraction of image that is background, usually in range [0.7, 0.995]"
		Prompt dist, "minimum distance (pixels) between any two peaks (use ~30 for unbinned Si)"
		DoPrompt/Help="3D-Xray Diffraction[Fit Peaks]" "peak fitting",imageName,fractionBkg,dist
		if (V_flag)
			return ""
		endif
		Wave image = $imageName
		printf "FitPeaks(%s,%g,%g)\r",imageName,fractionBkg,dist
		printIt = 1
	endif
	if (!WaveExists(image) || WaveDims(image)!=2 || !(fractionBkg==limit(fractionBkg,1e-3,1) || !(dist>0)))
		return ""
	endif
	Variable timer=startMSTimer

//	timer2 = startMSTimer
	Make/N=100/O hist_
	SetScale/P x 0,2,"", hist_
	Histogram/B=2 image,hist_
//	print "\t\tHistogram part took ",stopMSTimer(timer2)/1e6,"sec"

	Variable i = BinarySearchInterp(hist_, hist_(0)/10)
	i = numtype(i) ? numpnts(hist_)/2 : i
 	Variable cutOff = pnt2x(hist_,i)
// 	Variable cutOff = pnt2x(hist_,BinarySearchInterp(hist_, hist_(0)/10))
	printf "search for peaks above a threshold of %g\r",cutOff
//	print "cutOff=",cutOff											// threshold in image wave for finding peaks
//	timer2 = startMSTimer
	Duplicate/O image,FitPeaks_mask_bigger

	ImageThreshold/O/Q/I/T=(cutOff)/M=0 FitPeaks_mask_bigger
//	print "\t\tmaking mask part took ",stopMSTimer(timer2)/1e6,"sec"

	DoUpdate
//	timer2 = startMSTimer
	Variable id = dist>20 ? 5 : 1
//	Imagemorphology/O/E=5 BinaryDilation FitPeaks_mask_bigger
	Imagemorphology/O/E=(id) BinaryDilation FitPeaks_mask_bigger
//	print "\t\tImagemorphology part took ",stopMSTimer(timer2)/1e6,"sec"

//	timer2 = startMSTimer
	ImageAnalyzeParticles/A=5/M=3/Q stats FitPeaks_mask_bigger
	KillWaves/Z M_RawMoments,M_Particle
	Sort/R W_ImageObjArea,W_ImageObjArea,W_SpotX,W_SpotY,W_circularity,W_rectangularity,W_ImageObjPerimeter,W_xmin,W_xmax,W_ymin,W_ymax
//	print "\t\tImageAnalyzeParticles part took ",stopMSTimer(timer2)/1e6,"sec"

//	timer2 = startMSTimer
	i = reFit_GaussianPkList(image,dist)
//	print "\t\t refitting the peaks took ",stopMSTimer(timer2)/1e6,"sec"

	Variable sec0=stopMSTimer(timer)/1e6
	if (printIt)
		printf "fitted %d peaks,  this all took %g sec\r",i,sec0
		printf "  result stored in the wave '%s'\r", GetWavesDataFolder(image,2)
	endif
	Wave FullPeakList=FullPeakList

	KillWaves/Z hist_, W_sigma,W_coef,W_ParamConfidenceInterval, FitPeaks_mask_bigger
	KillWaves/Z W_xmin,W_xmax,W_ymin,W_ymax,W_circularity,W_rectangularity,W_ImageObjPerimeter,W_ImageObjArea,W_SpotX,W_SpotY
	return GetWavesDataFolder(image,2)
End



// remove bkg from an image, identify the peaks, and fit them all
Function/S OLD_removeBkg_FitPeaks(rawImage,fractionBkg,dist)
	Wave rawImage					// the raw image from the spe or HDF file
	Variable fractionBkg			// fraction of image that is background, use something like 0.99 or 0.995 (even 0.7 works OK)
	Variable dist					// minimum distance between spots, spots have to be at least this far apart (pixels)

	Variable printIt = ItemsInList(GetRTStackInfo(0))<2
	if (!WaveExists(rawImage) || WaveDims(rawImage)!=2 || !(fractionBkg==limit(fractionBkg,1e-3,1) || !(dist>0)))
		String imageName = SelectString(WaveExists(rawImage),"",NameOfWave(rawImage))
		fractionBkg = fractionBkg==limit(fractionBkg,1e-3,1) ? fractionBkg : 0.8
		dist = dist>0 ? dist : 30
//		Prompt imageName,"name of image with peaks to find",popup,WaveList("*",";","DIMS:2")
//		Prompt imageName,"name of image with peaks to find",popup,WaveList_Tags("imageFileName","rawIgorImage","*","DIMS:2")
		Prompt imageName,"name of image with peaks to find",popup,reverseList(WaveListClass("speImage","*","DIMS:2"))
		Prompt fractionBkg,"fraction of image that is background, usually in range [0.7, 0.995]"
		Prompt dist, "minimum distance (pixels) between any two peaks (use ~30 for unbinned Si)"
		DoPrompt/Help="3D-Xray Diffraction[Fit Peaks]" "peak fitting",imageName,fractionBkg,dist
		if (V_flag)
			return ""
		endif
		Wave rawImage = $imageName
		printf "FitPeaks(%s,%g,%g)\r",imageName,fractionBkg,dist
		printIt = 1
	endif
	if (!WaveExists(rawImage) || WaveDims(rawImage)!=2 || !(fractionBkg==limit(fractionBkg,1e-3,1) || !(dist>0)))
		return ""
	endif
//	Variable timer=startMSTimer
	Variable timer0=startMSTimer
	Variable thresh = numpnts(rawImage)*fractionBkg

	WaveStats/Q/M=1 rawImage
	Make/N=200/D/O FitPeaks_Hist
	Variable i
	for (i=-1;i<10 || numtype(i);)							// keep expanding the histogram range while i<10,   I want i to be bigger than 10
		SetScale/I x,V_min,V_max,"",FitPeaks_Hist			// lower high end of histogram to zoom in around zero
		Histogram/B=2 rawImage,FitPeaks_Hist
		Integrate/P FitPeaks_Hist/D=FitPeaks_Hist_INT	// integrate the histogram
		i=BinarySearchInterp(FitPeaks_Hist_INT,thresh)	// find point where FitPeaks_Hist_INT[i] == thresh
		V_max = numtype(i) ? V_max/2 : pnt2x(FitPeaks_Hist_INT,i+4)
	endfor
//	print "first part took ",stopMSTimer(timer)/1e6,"sec"
//	timer = startMSTimer

	thresh = DimDelta(FitPeaks_Hist,0)*i + DimOffset(FitPeaks_Hist,0)
	ImageThreshold/Q/T=(thresh)/M=0 rawImage					// create M_ImageThresh
	Wave M_ImageThresh=M_ImageThresh

	ImageStats/M=1 M_ImageThresh
	Variable fracBlack = V_avg*V_npnts/255/V_npnts
	//	printf "number of points in peak is %d out of %d points (%.2f%%)\r",fracBlack*V_npnts,V_npnts,fracBlack*100
	Make/N=(9,9)/O/B/U black
	black = 255
	for(;fracBlack<0.10;)											// loop until fracBlack is atleast 10%
		Imagemorphology/S=black/O BinaryDilation M_ImageThresh
		ImageThreshold/Q/T=(2)/M=0 M_ImageThresh
		ImageStats/M=1 M_ImageThresh
		fracBlack = V_avg*V_npnts/255/V_npnts
		//	printf "number of points in peak is %d out of %d points (%.2f%%)\r",fracBlack*V_npnts,V_npnts,fracBlack*100
	endfor
	ImageTransform/O invert M_ImageThresh						//	M_ImageThresh = !M_ImageThresh

//	print "second part took ",stopMSTimer(timer)/1e6,"sec"
//	timer = startMSTimer

//	Variable timer2 = startMSTimer
	ImageInterpolate /f={.125,.125} bilinear rawImage
	Redimension/S M_InterpolatedImage

	KillWaves/Z FitPeaks_rawImage_smaller
	Rename M_InterpolatedImage FitPeaks_rawImage_smaller		// make a source image using fewer pixels
//	print "\t\t1st interpolate part took ",stopMSTimer(timer2)/1e6,"sec"

//	timer2 = startMSTimer
 	ImageInterpolate /f={.125,.125} bilinear M_ImageThresh		// make a mask with fewer pixels
	KillWaves/Z FitPeaks_mask_smaller
	Rename M_InterpolatedImage FitPeaks_mask_smaller
	Redimension/B/U FitPeaks_mask_smaller
	ImageRemoveBackground/F/R=FitPeaks_mask_smaller/P=2/W FitPeaks_rawImage_smaller	// creates M_RemovedBackground, uses 3rd order polynomial
	Wave W_BackgroundCoeff=W_BackgroundCoeff
	//	print W_BackgroundCoeff
//	print "\t\t2nd interpolate part took ",stopMSTimer(timer2)/1e6,"sec"

//	timer2 = startMSTimer
	ImageInterpolate /f={8,8} bilinear M_RemovedBackground		// creates M_InterpolatedImage, a bkg with full number of pixels
	Wave M_InterpolatedImage=M_InterpolatedImage

	String outName = NameOfWave(rawImage)+"NoBkg"
	String fldr = GetWavesDataFolder(rawImage,1)
	if (CheckName(fldr+outName,1))
		outName = CleanupName(outName,0)
	endif
	outName = fldr+outName
	Duplicate/O rawImage $outName
	Wave noBkg = $outName

	String wnote=ReplaceStringByKey("rawIgorImage",note(rawImage),GetWavesDataFolder(rawImage,2),"=")
	wnote = ReplaceNumberByKey("fractionBkg",wnote,fractionBkg,"=")
	wnote = ReplaceNumberByKey("minSpotSeparation",wnote,dist,"=")
	wnote = ReplaceStringByKey("waveClass",wnote,"speImageNoBkg","=")
	Note/K noBkg, wnote
	Redimension/I noBkg											// signed int32 version of rawImage
	noBkg = rawImage - M_InterpolatedImage[p][q]					// p,q because M_InterpolatedImage is smaller
//	print "\t\t",DimSize(rawImage,0),DimSize(rawImage,1)
//	print "\t\t",DimSize(noBkg,0),DimSize(noBkg,1)
//	print "\t\t",DimSize(M_InterpolatedImage,0),DimSize(M_InterpolatedImage,1)
//	print "\t\t3rd interpolate part took ",stopMSTimer(timer2)/1e6,"sec"


//	timer2 = startMSTimer
	Make/N=100/O hist_
	SetScale/P x 0,2,"", hist_
	Histogram/B=2 noBkg,hist_
//	print "\t\tHistogram part took ",stopMSTimer(timer2)/1e6,"sec"

	i = BinarySearchInterp(hist_, hist_(0)/10)
	i = numtype(i) ? numpnts(hist_)/2 : i
 	Variable cutOff = pnt2x(hist_,i)
// 	Variable cutOff = pnt2x(hist_,BinarySearchInterp(hist_, hist_(0)/10))
	printf "after background removal search for peaks above a threshold of %g\r",cutOff
//	print "cutOff=",cutOff											// threshold in NoBkg wave for finding peaks
//	timer2 = startMSTimer
	Duplicate/O noBkg,FitPeaks_mask_bigger

	ImageThreshold/O/Q/I/T=(cutOff)/M=0 FitPeaks_mask_bigger
//	print "\t\tmaking mask part took ",stopMSTimer(timer2)/1e6,"sec"

	DoUpdate
//	timer2 = startMSTimer
	Variable id = dist>20 ? 5 : 1
//	Imagemorphology/O/E=5 BinaryDilation FitPeaks_mask_bigger
	Imagemorphology/O/E=(id) BinaryDilation FitPeaks_mask_bigger
//	print "\t\tImagemorphology part took ",stopMSTimer(timer2)/1e6,"sec"

//	timer2 = startMSTimer
	ImageAnalyzeParticles/A=5/M=3/Q stats FitPeaks_mask_bigger
	KillWaves/Z M_RawMoments,M_Particle
	Sort/R W_ImageObjArea,W_ImageObjArea,W_SpotX,W_SpotY,W_circularity,W_rectangularity,W_ImageObjPerimeter,W_xmin,W_xmax,W_ymin,W_ymax
//	print "\t\tImageAnalyzeParticles part took ",stopMSTimer(timer2)/1e6,"sec"

//	timer2 = startMSTimer
	i = reFit_GaussianPkList(noBkg,dist)

//	print "\t\t refitting the peaks took ",stopMSTimer(timer2)/1e6,"sec"

//	print "fancy part took ",stopMSTimer(timer)/1e6,"sec"	
	Variable sec0=stopMSTimer(timer0)/1e6
	if (printIt)
		printf "found and fitted %d peaks,  this all took %g sec\r",i,sec0
//		print "peak finding and fitting all took ",sec0,"sec"	
		printf "  result stored in the wave '%s'\r", GetWavesDataFolder(noBkg,2)
	endif
	Wave FullPeakList=FullPeakList

	KillWaves/Z hist_, W_sigma,W_coef,W_ParamConfidenceInterval
	KillWaves/Z FitPeaks_rawImage_smaller,FitPeaks_mask_smaller, M_ImageThresh, M_RemovedBackground, M_InterpolatedImage, W_BackgroundCoeff
	KillWaves/Z FitPeaks_Hist_INT,FitPeaks_Hist, black, FitPeaks_mask_bigger
	KillWaves/Z W_xmin,W_xmax,W_ymin,W_ymax,W_circularity,W_rectangularity,W_ImageObjPerimeter,W_ImageObjArea,W_SpotX,W_SpotY
	return GetWavesDataFolder(noBkg,2)
End



// Fit every peak from the output of ImageAnalyzeParticles,  it sets the string pkList based on the fits
// it also sets FullPeakList[][11] which contains the exact info about the fit of each peak
Static Function reFit_GaussianPkList(image,dist)
	Wave image
	Variable dist						// minimum distance between spots (pixels)
	Wave W_spotX=W_spotX, W_spotY=W_spotY
	Wave W_xmin=W_xmin, W_xmax=W_xmax, W_ymin=W_ymin, W_ymax=W_ymax
	if (exists("pkLIst")!=2)
		String/G pkLIst=""
	endif
	SVAR pkList = pkList
	pkList = ""
	String str
	Variable Nx=DimSize(image,0), Ny=DimSize(image,1)
	Variable j,i,N=numpnts(W_spotX),err
	Variable fw= ( 2*sqrt(2*ln(2)) )			// this factor corrects for 2-d sigma to FWHM, apply to K3, K5, and the corresponding errors
	Variable fwx,fwy,sigX,sigY
	Variable xx,yy, tooClose, Nu=0
	Variable left,right,top,bot
	Make/N=(N,11)/O FullPeakList
	String wnote = ReplaceStringByKey("fittedIgorImage",note(image),GetWavesDataFolder(image,2),"=")
	wnote = ReplaceStringByKey("waveClass",wnote,"FittedPeakList","=")
	Note/K FullPeakList,wnote
	SetDimLabel 1,0,x0,FullPeakList		;	SetDimLabel 1,1,y0,FullPeakList
	SetDimLabel 1,2,x0Err,FullPeakList	;	SetDimLabel 1,3,y0Err,FullPeakList
	SetDimLabel 1,4,fwx,FullPeakList		;	SetDimLabel 1,5,fwy,FullPeakList
	SetDimLabel 1,6,fwxErr,FullPeakList	;	SetDimLabel 1,7,fwyErr,FullPeakList
	SetDimLabel 1,8,correlation,FullPeakList ;	SetDimLabel 1,9,correlationErr,FullPeakList
	SetDimLabel 1,10,area,FullPeakList
	FullPeakList = NaN
	for (i=0;i<N;i+=1)
		xx = (W_xmin[i]+W_xmax[i])/2
		yy = (W_ymin[i]+W_ymax[i])/2

		tooClose = 0
		for (j=0;j<i;j+=1)						// check only spots already fitted
			if ((xx-W_spotX[j])^2 + (yy-W_spotY[j])^2 <= dist^2)
				tooClose = 1
				break
			endif
		endfor
		if (tooClose)
			continue
		endif
		left = floor(W_xmin[i])
		right = ceil(W_xmax[i])
		top = floor(W_ymin[i])
		bot = ceil(W_ymax[i])
		if (right-left<dist/2)					// fit box cannot be too small
			left = xx - dist/2
			right = xx + dist/2
		endif
		if (bot-top<dist/2)
			top = yy - dist/2
			bot = yy + dist/2
		endif
		left = limit(left,0,Nx-1)
		right = limit(right,0,Nx-1)
		top = limit(top,0,Ny-1)
		bot = limit(bot,0,Ny-1)

		err = FitOneGaussianPeak(image,left,right,top,bot)	// returns 0 == OK,  non-zero if error
		if (err)
			continue
		endif
		Wave W_coef=W_coef, W_sigma=W_sigma
		xx = W_coef[2]
		yy = W_coef[4]
		fwx = W_coef[3] * fw
		fwy = W_coef[5] * fw
		sigX = W_sigma[2]
		sigY = W_sigma[4]
		if (fwx>dist || fwy>dist || fwx<1 || fwy<1 || sigX>dist || sigY>dist)	// cannot have spots wider than dist
			continue
		endif
		// printf "%d  %d  Gaussian peak in at (%.2f±%.2f, %.2f±%.2f) with FWHM of (%.2f±%.2f, %.2f±%.2f),   amp=%g\r",i,Nu,xx,sigX,yy,sigY,fwx,W_sigma[3]*fw,fwy,W_sigma[5]*fw,W_coef[1]*fwx*fwy

		tooClose = 0													// check that this spot is not too close to another already included
		for (j=0;j<Nu-1;j+=1)					// check all spots already fitted
			if ((FullPeakList[j][0]-xx)^2 + (FullPeakList[j][1]-yy)^2 <= dist^2)
				tooClose = 1
				break
			endif
		endfor
		if (tooClose)
			continue
		endif
		FullPeakList[Nu][0]=xx				;	FullPeakList[Nu][1]=yy
		FullPeakList[Nu][2]=sigX				;	FullPeakList[Nu][3]=sigY
		FullPeakList[Nu][4]=fwx				;	FullPeakList[Nu][5]=fwy
		FullPeakList[Nu][6]=W_sigma[3]*fw	;	FullPeakList[Nu][7]=W_sigma[5]*fw
		FullPeakList[Nu][8]=W_coef[6]		;	FullPeakList[Nu][9]=W_sigma[6]
		FullPeakList[Nu][10]=W_coef[1]*fwx*fwy
		sprintf str,"%.0f,%.0f;",xx,yy
//		sprintf str,"%d,%d;",round(xx),round(yy)
		pkList += str
		Nu += 1
	endfor
	Redimension/N=(Nu,11) FullPeakList
	KillWaves/Z W_sigma,W_coef,W_ParamConfidenceInterval
	return Nu
End


// This routine allows the editing of which peaks to fit
// Fit every peak in the list 'peaks',  it sets the string pkList based on the fits
// it also sets FullPeakList[][11] which contains the exact info about the fit of each peak
Function setFittedPeaksFromList(image,dist,peaks)
	Wave image
	Variable dist						// minimum distance between spots (pixels)
	String peaks						// list of peaks, basically this is pkList

	peaks = SelectString(strlen(peaks),StrVarOrDefault(":pkList",""),peaks)
	if (strlen(peaks)<1)
		DoAlert 0,"no list of peaks available in setFittedPeaksFromList(), cannot do anything"
		return 0
	endif
	if (!WaveExists(image) || !(dist>0))
		dist = dist>0 ? dist : 10
		String wList = reverseList(WaveListClass("speImage*","*","DIMS:2")), wName=""
		Prompt dist,"minimum distance between peaks (pixels)"
		Prompt wName,"Indexed Peak List",popup,wList
		if (dist>0 && ItemsInList(wList)==1)	// only need image, and there is only one image
			Wave image = $(StringFromList(0,wList))
		else
			if (dist>0)					// only need image, and need to ask
				DoPrompt "image to use",wName
				Wave image = $wName
			elseif (WaveExists(image))		// only need dist
				DoPrompt "dist between peakst",dist
			else								// need both image and dist
				DoPrompt "image to use",wName,dist
				Wave image = $wName
			endif
			if (V_flag)
				return 0
			endif
		endif
		if (strlen(peaks)>100)
			printf "setFittedPeaksFromList(%s,%g,\"%s ...\")\r",GetWavesDataFolder(image,2),dist,peaks[0,100]
		else
			printf "setFittedPeaksFromList(%s,%g,\"%s\")\r",GetWavesDataFolder(image,2),dist,peaks
		endif
	endif
	if (!WaveExists(image) || !(dist>0))
		return 0
	endif

	if (exists("pkLIst")!=2)
		String/G pkLIst=""
	endif
	SVAR pkList = pkList
	pkList = ""
	String imageName = GetWavesDataFolder(image,2)
	Variable Nx=DimSize(image,0), Ny=DimSize(image,1)
	Variable j,i,N=ItemsInList(peaks),err,xx,yy

	String str
	Make/N=(N)/O peakxx_,peakyy_,peakInt_		// re-order the points by peak intensity
	for (i=0;i<N;i+=1)
		str = StringFromList(i,peaks)
		xx = str2num(StringFromList(0,str,","))
		yy = str2num(StringFromList(1,str,","))
		peakxx_[i] = xx							// fill arrays with poistions
		peakyy_[i] = yy
		peakInt_[i] = image[xx][yy]				// and intensities
	endfor
	Sort/R peakInt_, peakxx_,peakyy_,peakInt_	// sort by intensity
	peaks=""
	for (i=0;i<N;i+=1)								// and reset the peak list
		sprintf str,"%.0f,%.0f;",peakxx_[i],peakyy_[i]
		peaks += str
	endfor
	KillWaves/Z peakxx_,peakyy_,peakInt_

	Variable fw= ( 2*sqrt(2*ln(2)) )			// this factor corrects for 2-d sigma to FWHM, apply to K3, K5, and the corresponding errors
	Variable fwx,fwy,sigX,sigY
	Variable Nu=0, tooClose
	Variable left,right,top,bot
	Make/N=(N,11)/O FullPeakList
	String wnote = ReplaceStringByKey("fittedIgorImage",note(image),GetWavesDataFolder(image,2),"=")
	wnote = ReplaceStringByKey("waveClass",wnote,"FittedPeakList","=")
	Note/K FullPeakList,wnote
	SetDimLabel 1,0,x0,FullPeakList		;	SetDimLabel 1,1,y0,FullPeakList
	SetDimLabel 1,2,x0Err,FullPeakList	;	SetDimLabel 1,3,y0Err,FullPeakList
	SetDimLabel 1,4,fwx,FullPeakList		;	SetDimLabel 1,5,fwy,FullPeakList
	SetDimLabel 1,6,fwxErr,FullPeakList	;	SetDimLabel 1,7,fwyErr,FullPeakList
	SetDimLabel 1,8,correlation,FullPeakList ;	SetDimLabel 1,9,correlationErr,FullPeakList
	SetDimLabel 1,10,area,FullPeakList
	FullPeakList = NaN
	for (i=0;i<N;i+=1)
		str = StringFromList(i,peaks)
		xx = str2num(StringFromList(0,str,","))
		yy = str2num(StringFromList(1,str,","))
		left = limit(xx-dist/2,0,Nx-1)
		right = limit(xx+dist/2,0,Nx-1)
		top = limit(yy-dist/2,0,Ny-1)
		bot = limit(yy+dist/2,0,Ny-1)
		err = FitOneGaussianPeak(image,left,right,top,bot)	// returns 0 == OK,  non-zero if error
		if (err)
			continue
		endif
		Wave W_coef=W_coef, W_sigma=W_sigma
		xx = W_coef[2]
		yy = W_coef[4]
		fwx = W_coef[3] * fw
		fwy = W_coef[5] * fw
		sigX = W_sigma[2]
		sigY = W_sigma[4]
		if (fwx>dist*2 || fwy>dist*2 || fwx<1 || fwy<1 || sigX>dist || sigY>dist)	// cannot have spots wider than 2*dist
			continue
		endif
		for (j=0,tooClose=0;j<Nu;j+=1)							// check that this spot is not too close to an existing spot
			if (abs(FullPeakList[j][0]-xx)<dist && abs(FullPeakList[j][1]-yy)<dist)
				tooClose = 1
			endif
		endfor
		if (tooClose)
			continue
		endif
		FullPeakList[Nu][0]=xx				;	FullPeakList[Nu][1]=yy
		FullPeakList[Nu][2]=sigX				;	FullPeakList[Nu][3]=sigY
		FullPeakList[Nu][4]=fwx				;	FullPeakList[Nu][5]=fwy
		FullPeakList[Nu][6]=W_sigma[3]*fw	;	FullPeakList[Nu][7]=W_sigma[5]*fw
		FullPeakList[Nu][8]=W_coef[6]		;	FullPeakList[Nu][9]=W_sigma[6]
		FullPeakList[Nu][10]=W_coef[1]*fwx*fwy
		sprintf str,"%.0f,%.0f;",xx,yy
		pkList += str
		Nu += 1
	endfor
	Redimension/N=(Nu,11) FullPeakList
	KillWaves/Z W_sigma,W_coef,W_ParamConfidenceInterval
	return Nu
End



Static Function FullPeakList2Qfile(FullPeakList,fname,pathName)	// convert peaks to a Qlist+intens, peaks file for analysis be Euler
	Wave FullPeakList
	String fname
	String pathName

	fname = SelectString(strlen(fname),"generic_Peaks.txt",fname)
	if (!WaveExists(FullPeakList))
		DoAlert 0, "input wave for FullPeakList2File() does not exists"
		return 1
	elseif (DimSize(FullPeakList,1)!=11)
		DoAlert 0, "the passed full peak list '"+NameOfWave(FullPeakList)+"' is not the right size"
		return 1
	endif
	Variable N=DimSize(FullPeakList,0)
	if (N<1)													// nothing to write
		return 1
	endif

	STRUCT microGeometry geo								// note, dd and yc are reset from wave note below if it exists
	if (FillGeometryStructDefault(geo))					//fill the geometry structure with default values
		DoAlert 0, "no geometry structure found, did you forget to set it?"
		return 1
	endif
	STRUCT crystalStructure xtal
	if (FillCrystalStructDefault(xtal))
		DoAlert 0, "no lattice structure found, did you forget to set it?"
		return 1
	endif

	Make/N=(N,4)/O FullPeakList2Qfile_Q					// create temp waves for the q info
	Wave Qs=FullPeakList2Qfile_Q							// holds Qhat and intensity
	Make/N=3/D/O FullPeakList2Qfile_qhat, FullPeakList2Qfile_qBL
	Wave qhat=FullPeakList2Qfile_qhat, qBL=FullPeakList2Qfile_qBL
	String wnote = note(FullPeakList)
	Variable startx,groupx, starty,groupy					// ROI of the original image
	startx = NumberByKey("startx",wnote,"=")
	groupx = NumberByKey("groupx",wnote,"=")
	starty = NumberByKey("starty",wnote,"=")
	groupy = NumberByKey("groupy",wnote,"=")
	startx = numtype(startx) ? 1 : startx
	groupx = numtype(groupx) ? 1 : groupx
	starty = numtype(starty) ? 1 : starty
	groupy = numtype(groupy) ? 1 : groupy
	Variable ddLocal = geo.ddOffset + NumberByKey("CCDy",wnote,"=") // get dd from the image (in wave note) if it is there
	geo.dd = ddLocal>0 ? ddLocal :  geo.dd
	Variable ycLocal = NumberByKey("yc",wnote,"=")						// get ycLocal from the image (in wave note) if it is there
	geo.ycent = ycLocal>0 ? ycLocal :  geo.ycent

	Variable i, m, px,py
	for (i=0,m=0;i<N;i+=1)
		px = (startx-1) + groupx*FullPeakList[i][0] + (groupx-1)/2		// change to un-binned pixels
		py = (starty-1) + groupy*FullPeakList[i][1] + (groupy-1)/2		// pixels are still zero based
		pixel2q(geo,px,py,qhat)						// in Wenge-coord system
		if (norm(qhat)>0)									// check for a valid Q
			wenge2IdealBeamline(geo,qhat,qBL)
			Qs[m][0,2] = qBL[q]
			Qs[m][3] = FullPeakList[i][10]					// the intensity
			m += 1
		endif
	endfor
	KillWaves/Z FullPeakList2Qfile_qhat, FullPeakList2Qfile_qBL
	N = m

	Variable refNum
	Open/C="R*ch"/P=$pathName/T="TEXT" refNum as fname
	if (!strlen(S_fileName))
		DoAlert 0, "Cannot open file '"+fname+"'"
		KillWaves/Z FullPeakList2Qfile_Q
		return 1
	endif

	String str
	Variable val
	fprintf refNum,"$PeaksFile\n"
	str = xtal.desc
	if (strlen(str))
		fprintf refNum,"$structureDesc		%s\n",xtal.desc
	endif
	fprintf refNum,"$latticeParameters	{ %g, %g, %g, %g, %g, %g }// using nm and degrees\n",xtal.a,xtal.b,xtal.c,xtal.alpha,xtal.beta,xtal.gam
	fprintf refNum,"$lengthUnit			nm					// length unit for lattice constants a,b,c\n"
	fprintf refNum,"$SpaceGroup			%d					// Structure number from International Tables\n",xtal.SpaceGroup
	if (stringmatch(StrVarOrDefault("root:Packages:micro:indexingExecutableMac",""),"EulerOrig"))
		fprintf refNum,"$latticeStructure		%d					// Structure number from International Tables\n",xtal.SpaceGroup
	endif
	str = "\t// {element  x y z occupancy}"
	for (i=0;i<xtal.N;i+=1)
		fprintf refNum,"$AtomDesctiption%d	{%s  %g %g %g %g}%s\n",i+1,xtal.atom[i].name,xtal.atom[i].x,xtal.atom[i].y,xtal.atom[i].z,xtal.atom[i].occ,str
		str = ""
	endfor

	str = StringByKey("imageFileName",wnote,"=")
	if (strlen(str))
		fprintf refNum,"$imageFileName			%s		// name of image file read in by Igor\n",str
	endif
	fprintf refNum,"$peakListWave	%s	// Igor wave containng list of peak positions\n",GetWavesDataFolder(FullPeakList,2)
	val = NumberByKey("exposure",wnote,"=")
	if (!numtype(val))
		fprintf refNum,"$exposure			%g					// CCD exposure of original image (sec)\n",val
	endif
	str = StringByKey("dateExposed",wnote,"=")
	if (strlen(str))
		fprintf refNum,"$dateExposed		%s			// date of original CCD image\n",str
	endif
	str = StringByKey("rawIgorImage",wnote,"=")
	if (strlen(str))
		fprintf refNum,"$rawIgorImage		%s		// name of Igor image wave\n",str
	endif
	str = StringByKey("fittedIgorImage",wnote,"=")
	if (strlen(str))
		fprintf refNum,"$fittedIgorImage	%s	// name of flattened Igor image wave used to fit peaks\n",str
	endif
	val = NumberByKey("fractionBkg",wnote,"=")
	if (!numtype(val))
		fprintf refNum,"$fractionBkg		%g					// in peak analysis, the given fraction of image that is bkg\n",val
	endif
	val = NumberByKey("minSpotSeparation",wnote,"=")
	if (!numtype(val))
		fprintf refNum,"$minSpotSeparation	%g					// in peak analysis, min separaion between spots (pixels)\n",val
	endif
	val = NumberByKey("minPeakWidth",wnote,"=")
	if (!numtype(val))
		fprintf refNum,"$minPeakWidth	%g					// in peak analysis, min allowed fwhm for a spot (pixels)\n",val
	endif
	val = NumberByKey("maxPeakWidth",wnote,"=")
	if (!numtype(val))
		fprintf refNum,"$maxPeakWidth	%g					// in peak analysis, max allowed fwhm for a spot (pixels)\n",val
	endif
	val = NumberByKey("totalPeakIntensity",wnote,"=")
	if (!numtype(val))
		fprintf refNum,"$totalPeakIntensity	%g					// sum of the areas of all fitted peaks\n",val
	endif
	val = NumberByKey("totalIntensity",wnote,"=")
	if (!numtype(val))
		fprintf refNum,"$totalIntensity	%g						// sum of all the intensity in the image\n",val
	endif
	val = NumberByKey("startx",wnote,"=")					// write the ROI of the image
	if (!numtype(val))
		fprintf refNum,"$startx	%g								// ROI of the raw image file\n",val
	endif
	val = NumberByKey("groupx",wnote,"=")
	if (!numtype(val))
		fprintf refNum,"$groupx	%g\n",val
	endif
	val = NumberByKey("starty",wnote,"=")
	if (!numtype(val))
		fprintf refNum,"$starty	%g\n",val
	endif
	val = NumberByKey("groupy",wnote,"=")
	if (!numtype(val))
		fprintf refNum,"$groupy	%g\n",val
	endif
	val = NumberByKey("CCDy",wnote,"=")
	if (!numtype(val))
		fprintf refNum,"$CCDy	%g\n",val
	endif
	fprintf refNum,"\n// the following table contains xyz compotnents of G^ and the integral of the peak\n"
	fprintf refNum,"$N_Ghat+Intens 	%d		// number of G^ vectors\n",N
	for (i=0;i<N;i+=1)
		fprintf refNum, "% .7f, % .7f, % .7f,     %g\n",Qs[i][0],Qs[i][1],Qs[i][2],Qs[i][3]
	endfor
	Close refNum
	KillWaves/Z FullPeakList2Qfile_Q
	return 0
End


Function TableFullPeakList(w)
	Wave w
	if (!WaveExists(w))
		String list = reverseList(WaveListClass("FittedPeakList","*",""))
		Variable i = ItemsInLIst(list)
		if (i<1)
			DoAlert 0, "No waves of this type in current folder"
			return 1
		elseif (i==1)
			Wave w = $StringFromList(0,list)
		else
			String wName = StringFromList(i-1,list)
			Prompt wName, "Wave with list of fitted peaks",popup,list
			DoPrompt "FullPeakList",wName
			if (V_flag)
				return 1
			endif			
			Wave w = $wName
		endif
		if (!WaveExists(w))
			DoAlert 0, "No waves of name "+wName
			return 1
		endif
	endif

	String table = FindTableWithWave(w)
	if (strlen(table)>0)
		DoWindow/F $table
		return 0
	endif
	Edit/K=1/W=(5,44,891,684) w.ld
	ModifyTable width(Point)=30,title(w.d)="Fitted Peaks",width(w.l)=20
	ModifyTable format(w.d)=3,width(w.d)=74
End

Function TableFullPeakIndexed(w)
	Wave w
	if (!WaveExists(w))
		String list = reverseList(WaveListClass("IndexedPeakList","*",""))
		Variable i = ItemsInLIst(list)
		if (i<1)
			DoAlert 0, "No waves of this type in current folder"
			return 1
		elseif (i==1)
			Wave w = $StringFromList(0,list)
		else
			String wName = StringFromList(i-1,list)
			Prompt wName, "Wave with list of indexed peaks",popup,list
			DoPrompt "FullPeakIndexed",wName
			if (V_flag)
				return 1
			endif			
			Wave w = $wName
		endif
		if (!WaveExists(w))
			DoAlert 0, "No waves of name "+wName
			return 1
		endif
	endif

	String table = FindTableWithWave(w)
	if (strlen(table)>0)
		DoWindow/F $table
		return 0
	endif
	Edit/K=1/W=(5,44,891,684) w.ld
	ModifyTable width(Point)=30,title(w.d)="Indexed Peaks",width(w.l)=20
	ModifyTable format(w.d)=3,width(w.d)=74
End

Function TableFullPeakStrain(w)
	Wave w
	if (!WaveExists(w))
		String list = reverseList(WaveListClass("StrainPeakList","*",""))
		Variable i = ItemsInLIst(list)
		if (i<1)
			DoAlert 0, "No waves of this type in current folder"
			return 1
		elseif (i==1)
			Wave w = $StringFromList(0,list)
		else
			String wName = StringFromList(i-1,list)
			Prompt wName, "Wave with list of indexed peaks",popup,list
			DoPrompt "FullPeakIndexed",wName
			if (V_flag)
				return 1
			endif			
			Wave w = $wName
		endif
		if (!WaveExists(w))
			DoAlert 0, "No waves of name "+wName
			return 1
		endif
	endif

	String table = FindTableWithWave(w)
	if (strlen(table)>0)
		DoWindow/F $table
		return 0
	endif
	Edit/K=1/W=(5,44,1035,681) w.ld
	ModifyTable width(Point)=30,title(w.d)="Strain Peaks",width(w.l)=20
	ModifyTable format(w.d)=3,width(w.d)=74
End



// find hkl that is paralllel to the input and allowed
Static Function CloseAllowedhkl(hkl)	// in and out can be the same wave
	Wave hkl							// real numbers on intput, set to h,k,l on output

	Variable maxhkl = NumVarOrDefault("root:MAX_HKL_Displayed",25 )
	maxhkl = ((maxhkl>0) && numtype(maxhkl)==0) ? maxhkl : 25
	Variable xin=hkl[0],yin=hkl[1],zin=hkl[2]		// save input values
	hkl = abs(hkl)
	Variable i, maxVal = WaveMax(hkl)
	Variable ibest=0, dotMax=-Inf, dot
	for (i=1;i<=maxhkl;i+=1)
		hkl = {xin,yin,zin}
		hkl = round(hkl*i/maxVal)
		dot = (hkl[0]*xin + hkl[1]*yin + hkl[2]*zin)/norm(hkl)
		if (dot > dotMax)
			ibest = i
			dotMax = dot
		endif
	endfor
	hkl = {xin,yin,zin}
	hkl = round(hkl*ibest/maxVal)
	hkl = hkl[p]==0 ? 0 : hkl[p]		// avoid "-0"
	return 0
End
//Static  Function CloseAllowedhkl(hkl)	// in and out can be the same wave
//	Wave hkl							// real numbers on intput, set to h,k,l on output
//
//	Variable xin=hkl[0],yin=hkl[1],zin=hkl[2]		// save input values
//	hkl = abs(hkl)
//	Variable minNonZeroVal = WaveMax(hkl)/60	// if it is less than this, then it is a zero
//	hkl = (hkl[p]<minNonZeroVal) ? Inf : hkl[p]
//	WaveStats/Q/M=1 hkl							// do this to get V_minloc
//
//	hkl = {xin,yin,zin}
//	Variable minh = abs(hkl[V_minloc])			// location of smallest non-zero value, e.g. for {0,2,7}, minh=2
//	hkl = round(hkl[p]/minh)
//	hkl = hkl[p]==0 ? 0 : hkl[p]					// get rid of "-0"
//
//	xin=hkl[0]; yin=hkl[1]; zin=hkl[2]
//	lowestAllowedHKL(xin,yin,zin)					// convert hkl to allowed reflection
//	hkl = {xin,yin,zin}
//	return 0
//End
//Static  Function CloseAllowedhkl(in,out)	// in and out can be the same wave
//	Wave in,out								// input and output hkl
//	String wName=""
//	if (stringmatch(GetWavesDataFolder(in,2),GetWavesDataFolder(out,2)))
//		wName = UniqueName("CloseAllowedhkl",1,0)
//		Make/N=3/O/D $wName=in
//		Wave inTemp = $wName
//	else
//		Wave inTemp=in
//	endif
//
//	out = abs(inTemp)
//	Variable minNonZeroVal = WaveMax(out)/60	// if it is less than this, then it is a zero
//	out = (out[p]<minNonZeroVal) ? Inf : out[p]
//	WaveStats/Q/M=1 out
//	Variable minh = abs(inTemp[V_minloc])		// location of smallest non-zero value, e.g. for {0,2,7}, minh=2
//	out = round(inTemp[p]/minh)
//	out = out[p]==0 ? 0 : out[p]				// get rid of "-0"
//	Variable h=out[0], k=out[1], l=out[2]
//	lowestAllowedHKL(h,k,l)					// convert hkl to allowed reflection
//	out = {h,k,l}
//	KillWaves/Z $wName
//	return 0
//End
//Static  Function CloseAllowedhkl(in,out)
//	Wave in,out								// input and output hkl
//
//	out = abs(in)
//	Variable minNonZeroVal = WaveMax(out)/60	// if it is less than this, then it is a zero
//	out = (out[p]<minNonZeroVal) ? Inf : out[p]
//	WaveStats/Q/M=1 out
//	Variable minh = abs(in[V_minloc])		// location of smallest non-zero value, e.g. for {0,2,7}, minh=2
//	out = round(in[p]/minh)
//	out = out[p]==0 ? 0 : out[p]				// get rid of "-0"
//	Variable h=out[0], k=out[1], l=out[2]
//	lowestAllowedHKL(h,k,l)					// convert hkl to allowed reflection
//	out = {h,k,l}
//	return 0
//End
//Static Function CloseAllowedhkl(in,out)
//	Wave in,out								// input and output hkl
//	out = abs(round(in*60)/60)
//	out = out[p]>0 ?out[p] : Inf
//	WaveStats/Q/M=1 out
//	Variable minh = abs(in[V_minloc])		// smallest non-zero value, e.g. for {0,2,7}, minh=2
//	out = round(in[p]/minh*12)
//	out = out[p]==0 ? 0 : out[p]				// get rid of "-0"
//	Variable h=out[0], k=out[1], l=out[2]
//	lowestAllowedHKL(h,k,l)					// convert hkl to allowed reflection
//	out = {h,k,l}
//	return 0
//End
//// find hkl that is paralllel to the input and allowed
//Static Function CloseAllowedhkl(SpaceGroup,in,out)
//	Wave in,out								// input and output hkl
//	Variable SpaceGroup					// Space Group number from International tables
//	out = abs(round(in*60)/60)
//	out = out[p]>0 ?out[p] : Inf
//	WaveStats/Q/M=1 out
//	Variable minh = abs(in[V_minloc])		// smallest non-zero value, e.g. for {0,2,7}, minh=2
//	out = round(in[p]/minh*12)
//	out = out[p]==0 ? 0 : out[p]				// get rid of "-0"
//	Variable h=out[0], k=out[1], l=out[2]
//	lowestAllowedHKL(h,k,l,SpaceGroup)	// convert hkl to allowed reflection
//	out = {h,k,l}
//	return 0
//End



Function EnergyOfhkl(FullPeakIndexed,pattern,h,k,l)
	Wave FullPeakIndexed
	Variable pattern
	Variable h,k,l

	Variable i
	if (!WaveExists(FullPeakIndexed))
		String wList = reverseList(WaveListClass("IndexedPeakList","*",""))
		i = ItemsInList(wList)
			if (i==1)						// only one, so choose it
			Wave FullPeakIndexed = $(StringFromList(0,wList))
		elseif(i>0)						// more than 1, so ask
			String wName=""
			Prompt wName,"Indexed Peak List",popup,wList
			DoPrompt "Indexed Peak List",wName
			if (V_flag)
				return NaN
			endif
			Wave FullPeakIndexed = $wName
		endif
	endif
	if (!WaveExists(FullPeakIndexed))
		DoAlert 0, "cannot find wave of indexed peak information 'FullPeakIndexed'"
		return NaN
	endif

	Variable Np = max(DimSize(FullPeakIndexed,2),1)
	pattern = (Np==1) ? 0 : pattern



	// need to fix this so it takes the user data from the graph, not the panel
	if (numtype(pattern))
		pattern = NumberByKey("patternNum",Getuserdata("","","Indexing"),"=")
	endif



	if (Np>1 && (pattern<0 || pattern>=Np || numtype(pattern)))
		pattern = pattern>=0 ? pattern : 0
		Prompt pattern, "pattern to choose [0,"+num2istr(Np-1)+",]",popup,expandRange("0-"+num2istr(Np-1),";")
		DoPrompt "pattern number",pattern
		if (V_flag || !(pattern>=0))
			return NaN
		endif
		pattern -= 1
	endif
	pattern = pattern>=0 ? pattern : 0
	if (numtype(h+k+l))
		Prompt h, "H"
		Prompt k, "K"
		Prompt l, "L"
		DoPrompt "(hkl)",h,k,l
		if (V_flag || numtype(h+k+l))
			return NaN
		endif
	endif
	if (numtype(h+k+l))
		return NaN
	endif

	Variable d									// d-spacing
	Variable as0,as1,as2						// componenets of a*
	Variable bs0,bs1,bs2						// componenets of b*
	Variable cs0,cs1,cs2						// componenets of c*
	Make/N=3/O/D EnergyOfhkl_qvec,EnergyOfhkl_ki={0,0,1}
	Wave qvec=EnergyOfhkl_qvec, ki=EnergyOfhkl_ki
	String wnote = note(FullPeakIndexed)
	String recip_lattice = StringByKey("recip_lattice"+num2istr(pattern),wnote,"=")	// the starting point
	sscanf recip_lattice, "{{%g,%g,%g}{%g,%g,%g}{%g,%g,%g}}",as0,as1,as2,bs0,bs1,bs2,cs0,cs1,cs2
	if (V_flag!=9)
		DoAlert 0, "Unable to read recip_lattice"
		return NaN
	endif
	qvec[0] = h*as0 + k*bs0 + l*cs0
	qvec[1] = h*as1 + k*bs1 + l*cs1
	qvec[2] = h*as2 + k*bs2 + l*cs2
	d = 2*PI/normalize(qvec)

	Variable sineTheta = -MatrixDot(qvec,ki)	// sin(theta) = -ki dot qhat
	Variable keV = hc/(2*d*sineTheta)			// energy

	if (ItemsInList(GetRTStackInfo(0))<2 || stringmatch(StringFromList(0,GetRTStackInfo(0)),"IndexButtonProc"))
		Variable px=NaN, py=NaN
		STRUCT microGeometry geo
		if (!FillGeometryStructDefault(geo))			//fill the geometry structure with test values
			Variable startx,groupx, starty,groupy, ddLocal, ycLocal
			startx = NumberByKey("startx",wnote,"=")
			groupx = NumberByKey("groupx",wnote,"=")
			starty = NumberByKey("starty",wnote,"=")
			groupy = NumberByKey("groupy",wnote,"=")
			startx = numtype(startx) ? 1 : startx
			groupx = numtype(groupx) ? 1 : groupx
			starty = numtype(starty) ? 1 : starty
			groupy = numtype(groupy) ? 1 : groupy
			ddLocal = geo.ddOffset + NumberByKey("CCDy",wnote,"=")// get dd from the image (in wave note) if it is there
			geo.dd = ddLocal>0 ? ddLocal :  geo.dd
			ycLocal = NumberByKey("yc",wnote,"=")					// get ycLocal from the image (in wave note) if it is there
			geo.ycent = ycLocal>0 ? ycLocal :  geo.ycent
			IdealBeamLine2wenge(geo,qvec,qvec)
			Variable/C pz = q2pixel(geo,qvec)
			px = (real(pz)-(startx-1)-(groupx-1)/2)/groupx		// change to binned pixels
			py = (imag(pz)-(starty-1)-(groupy-1)/2)/groupy		// pixels are still zero based
		endif
		printf "for %s pattern #%d,  d[(%s)] = %.9g nm,   E(%s) = %.4f keV   (theta=%.3f¡)",NameOfWave(FullPeakIndexed),pattern,hkl2str(h,k,l),d,hkl2str(h,k,l),keV,asin(sineTheta)*180/PI
		if (numtype(px+py)==0)
			printf "       should be at pixel [%.2f, %.2f]",px,py
		endif
		printf "\r"
	endif
	KillWaves/Z EnergyOfhkl_qvec, EnergyOfhkl_ki
	return keV
End




// find angle between two Qvectors from points on an image
Function AngleBetweenTwoPoints(image,px0,py0,px1,py1)
	Wave image
	Variable px0,py0				// first pixel (in binning of image)
	Variable px1,py1				// second pixel (in binning of image)

	Variable printIt = (ItemsInList(GetRTStackInfo(0))<2)

	if (!WaveExists(image))		// use image on top graph if not passed
		Wave image = ImageNameToWaveRef("",StringFromList(0,ImageNameList("",";")))
	endif
	if (!WaveExists(image))
		if (printIt)
			DoAlert 0, "could not find the image"
		endif
		return 1
	endif
	if (numtype(px0+py0+px1+py1))		// bad pixels, use cursors A and B
		String infoA=CsrInfo(A), infoB=CsrInfo(B)
		if (strlen(infoA)>1 && strlen(infoB)>1)
			px0 = hcsr(A) ;		py0 = vcsr(A)
			px1 = hcsr(B) ;	py1 = vcsr(B)
		endif
	endif
	if (numtype(px0+py0+px1+py1))		// bad pixels, use cursors A and B
		if (printIt)
			DoAlert 0, "could not get the two pixel positions"
		endif
		return 1
	endif
	STRUCT microGeometry geo
	if (FillGeometryStructDefault(geo))					//fill the geometry structure with test values
		if (printIt)
			DoAlert 0, "no geometry structure found, did you forget to set it?"
		endif
		return 1
	endif

	String wnote = note(image)
	Variable startx,groupx, starty,groupy					// ROI of the original image
	startx = NumberByKey("startx",wnote,"=")
	groupx = NumberByKey("groupx",wnote,"=")
	starty = NumberByKey("starty",wnote,"=")
	groupy = NumberByKey("groupy",wnote,"=")
	startx = numtype(startx) ? 1 : startx
	groupx = numtype(groupx) ? 1 : groupx
	starty = numtype(starty) ? 1 : starty
	groupy = numtype(groupy) ? 1 : groupy
	Variable ddLocal = geo.ddOffset + NumberByKey("CCDy",wnote,"=")	// get dd from the image (in wave note) if it is there
	geo.dd = ddLocal>0 ? ddLocal :  geo.dd
	Variable ycLocal = NumberByKey("yc",wnote,"=")						// get ycLocal from the image (in wave note) if it is there
	geo.ycent = ycLocal>0 ? ycLocal :  geo.ycent

	Make/N=3/O/D AngleBetweenTwoPoints_qhat1, AngleBetweenTwoPoints_qhat2
	Wave qhat1=AngleBetweenTwoPoints_qhat1, qhat2=AngleBetweenTwoPoints_qhat2

	Variable angle					// angle between points (degree)
	Variable px,py
	px = (startx-1) + groupx*px0 + (groupx-1)/2		// change to un-binned pixels
	py = (starty-1) + groupy*py0 + (groupy-1)/2		// pixels are still zero based
	pixel2q(geo,px,py,qhat1)										// in Wenge-coord system
	px = (startx-1) + groupx*px1 + (groupx-1)/2		// change to un-binned pixels
	py = (starty-1) + groupy*py1 + (groupy-1)/2		// pixels are still zero based
	pixel2q(geo,px,py,qhat2)										// in Wenge-coord system
	Variable cosine = limit(MatrixDot(qhat1,qhat2),-1,1)
	KillWaves/Z AngleBetweenTwoPoints_qhat1, AngleBetweenTwoPoints_qhat2
	angle = acos(cosine)*180/PI
	if (printIt)
		printf "Between points [%.2f, %.2f],  and [%.2f, %.2f]  the angle between Q vectors is %.3f¡\r",px0,py0,px1,py1,angle
	endif
	return angle
End







// if the image cannot be found in existing graphs, it returns empty string
Function/S FindGraphWithImage(image)	// find the graph window which contains the specified image
	Wave image
	if (!WaveExists(image))
		return ""
	endif
	String name, name0=GetWavesDataFolder(image,2)
	String ilist, win,wlist = WinList("*",";","WIN:1")
	Variable i,Ni, m,Nm=ItemsInList(wlist)
	for (m=0;m<Nm;m+=1)
		win = StringFromList(m,wlist)
		ilist = ImageNameList(win,";")		// list of images in this window
		Ni = ItemsInList(ilist)
		for (i=0;i<ItemsInList(ilist);i+=1)
			name = GetWavesDataFolder(ImageNameToWaveRef(win,StringFromList(i,ilist)),2)
			if (stringmatch(name,name0))
				return win
			endif
		endfor
	endfor
	return ""
End

Static Function/S FindTableWithWave(w0)	// find the table window which contains the specified wave
	Wave w0
	if (!WaveExists(w0))
		return ""
	endif
	String name0=GetWavesDataFolder(w0,2)		// full path name of desired wave
	String table,tlist = WinList("*",";","WIN:2")
	Variable i, m,Nm=ItemsInList(tlist)
	for (m=0;m<Nm;m+=1)						// loop over all displayed tables
		table = StringFromList(m,tlist)				// table name
		i = 0
		do											// check each wave in this table
			Wave wav = WaveRefIndexed(table,i,1)	// 1 means table data
			if (!WaveExists(wav))					// no more waves in this table
				break
			endif
			if (stringmatch(GetWavesDataFolder(wav,2),name0))
				return table						// found our wave, return table name
			endif
			i += 1
		while(1)
	endfor
	return ""										// no exising displayed table found
End


// this routine moved to Utility_JZT.ipf
//// return a list of waves in current folder having all tags in 'require', and none of the tags in 'avoid'
//// This is similar to WaveList(), but with a finer selection
//Function/T WaveListClass(waveClassList,search,options)
//	String waveClassList				// a list of acceptable wave classes (semi-colon separated)
//	String search						// same as first argument in WaveList()
//	String options						// same as last argument in WaveList()
//
//	String in = WaveList(search,";",options), out=""
//	String name,key, wnote
//	Variable m
//	for (m=0, name=StringFromList(0,in); strlen(name); m+=1,name=StringFromList(m,in))
//		wnote = note($name)
//		if (WhichListItem(StringByKey("waveClass",wnote,"="),waveClassList)<0)
//			continue
//		endif
//		out += name+";"
//	endfor
//	return out
//End
//
// return a list of waves in current folder having all tags in 'require', and none of the tags in 'avoid'
// This is similar to WaveList(), but with a finer selection
Function/T WaveList_Tags(require,avoid,search,options)
	String require						// list of required tags in wave note
	String avoid							// list of tags to avoid in wave note
	String search						// same as first argument in WaveList()
	String options						// same as last argument in WaveList()

	String in = WaveList(search,";",options), out=""
	String name,key, wnote
	Variable m,i,Na=ItemsInList(avoid),Nr=ItemsInList(require), bad
	for (m=0, name=StringFromList(0,in); strlen(name); m+=1,name=StringFromList(m,in))
		wnote = note($name)
		bad = 0

		for (i=0;i<Nr;i+=1)				// check for required keys
			key = StringFromList(i,require)
			if (strlen(StringByKey(key,wnote,"="))==0)
				bad = 1						// a required key not found, set bad and break
				break
			endif
		endfor

		for (i=0;i<Na;i+=1)				// check for keys to avoid
			key = StringFromList(i,avoid)
			if (strlen(StringByKey(key,wnote,"=")))
				bad = 1						// found a key to avoid, so stop and set bad flag
				break
			endif
		endfor
		if (!bad)
			out += name+";"
		endif
	endfor
	return out
End



// =========================================================================
// =========================================================================
//	Start of Index strain refinement

Function/T DeviatoricStrainRefine(pattern,constrain,[coords])
	Variable pattern											// pattern number, usually 0
	String constrain											// constraint on optimization, "111111", a 1 is refine, a 0 is keep constant
	Variable coords											// coordinate system to pass in return {1=BL, 2=XHF, 3=Sample (outward surface normal)}
	coords = ParamIsDefault(coords) ? 0 : round(coords)	// Beam Line system is the default
	if (coords<0 || coords>2)									// coords must be 1, 2, or 3
		return ""
	endif
	Wave FullPeakList=FullPeakList, FullPeakIndexed=FullPeakIndexed
	if (!WaveExists(FullPeakList) || !WaveExists(FullPeakIndexed))
		return ""
	endif
	STRUCT microGeometry geo
	if (FillGeometryStructDefault(geo))					//fill the geometry structure with test values
		DoAlert 0, "no geometry structure found, did you forget to set it?"
		return ""
	endif
	Variable printIt=0, Npatterns=DimSize(FullPeakIndexed,2)

	Variable aFit,bFit,cFit,alphaFit,betFit,gamFit=NaN
	sscanf constrain, "%1d%1d%1d%1d%1d%1d",aFit,bFit,cFit,alphaFit,betFit,gamFit
	if (V_flag!=6 || !((aFit+bFit+cFit+alphaFit+betFit+gamFit)>0) || str2num(constrain[0,2])>=111)
		aFit = NaN												// flags bad input
	endif

	if (!(pattern>=0) || numtype(aFit))					// need to ask user
		pattern = limit(pattern,0,Npatterns-1)
		if (StartStrainPanel(pattern,Npatterns,aFit,bFit,cFit,alphaFit,betFit,gamFit)<0)
			return ""
		endif
		printIt = 1
	endif
	aFit = aFit ? 1 : 0 ;				bFit = bFit ? 1 : 0 ;			cFit = cFit ? 1 : 0		// ensure only 0 or 1
	alphaFit = alphaFit ? 1 : 0 ;		betFit = betFit ? 1 : 0 ;		gamFit = gamFit ? 1 : 0
	Variable Nfit = aFit+bFit+cFit+alphaFit+betFit+gamFit
	if (Nfit<1 || (aFit&&bFit&&cFit))
		return ""												// this is deviatoric, you cannot fit all 3 lengths
	endif
	sprintf constrain "%d%d%d%d%d%d", aFit,bFit,cFit,alphaFit,betFit,gamFit		// constrain is for {a,b,c,alpha,beta,gamma}

	Variable N=DimSize(FullPeakIndexed,0)
	printIt = printIt || ((ItemsInList(GetRTStackInfo(0))<2) || stringmatch(StringFromList(0,GetRTStackInfo(0)),"IndexButtonProc"))
	KillWaves/Z epsilon,epsilonBL,epsilonXHF,epsilonSample
	if (N<(Nfit+3))											// with fixed Vc, there are (Nfit-1) free parameters, and 3 rotations
		return ""
	endif
	Wave Qs = $(FullPeakList2Qwave(FullPeakList))
	if (!WaveExists(Qs))
		return ""
	endif
	String funcName = "Indexing#latticeMismatch"+num2istr(Nfit+3)	// need 3 extra for the rotation, which is always there
	//
	// done with I/O, now setup the calculation

	// from initial recip lattice, find lattice constants, and initial rotation (stored as rx,ry,rz)
	Make/N=3/O/D hkl_latticeMismatch, axis_latticeMismatch
	Wave hkl=hkl_latticeMismatch, axis=axis_latticeMismatch
	Variable as0,as1,as2,bs0,bs1,bs2,cs0,cs1,cs2
	Make/N=(3,3)/O/D RLmeas_latticeMismatch
	Wave  RLmeas=RLmeas_latticeMismatch
	sscanf StringByKey("recip_lattice0",note(FullPeakIndexed),"="), "{{%g,%g,%g}{%g,%g,%g}{%g,%g,%g}}",as0,as1,as2,bs0,bs1,bs2,cs0,cs1,cs2
	RLmeas[0][0] = as0 ;		RLmeas[0][1] = bs0 ;		RLmeas[0][2] = cs0		// the original measured RL
	RLmeas[1][0] = as1 ;		RLmeas[1][1] = bs1 ;		RLmeas[1][2] = cs1
	RLmeas[2][0] = as2 ;		RLmeas[2][1] = bs2 ;		RLmeas[2][2] = cs2
	Make/N=6/O/D optimize_LatticeConstantsWave
	Wave LC=optimize_LatticeConstantsWave
	Variable Vc = RL2latticeConstants(RLmeas,LC)			// calc original lattice constants from RL, stored in LC, save original Vc too
	forceLattice(LC,NumberByKey("SpaceGroup",note(FullPeakIndexed),"="),Vc=Vc)	// force lattice constants to match Space Group exactly, preserve Vc

	Make/N=(3,3)/O/D DL_latticeMismatch, RL_latticeMismatch, RL0_latticeMismatch, rho_latticeMismatch
	Wave DL=DL_latticeMismatch, RL=RL_latticeMismatch, RL0=RL0_latticeMismatch, rho=rho_latticeMismatch
	RLfromLatticeConstants(LC,DL,RL0)					// make reference RL from lattice constants (used to compute rho), this RL exactly matches Space Group
	MatrixOp/O rho = RLmeas x Inv(RL0)					// RLmeas = rho x RL0,  the rotation (or almost a perfect rotation matrix)
	Variable angle = axisOfMatrix(rho,axis)				// rho is almost a  perfect rotation matrix, by remaking it, it will be a perfect rotation matrix
	if (numtype(angle))
		return ""
	endif
	axis *= angle												// length of axis is now rotation angle (rad)
	rotationMatAboutAxis(axis,angle,rho)					// forces rho to be a perfect rotation matrix
	MatrixOp/O RLmeas = rho x RL0							// ensure that starting point perfectly agrees with Space Group

	Make/N=(3,3)/O/D Ameas_Deviatoric, Astart_Deviatoric	// Vcartesian = A x Vcell,   used to make the strain tensor
	Wave Ameas=Ameas_Deviatoric, Astart=Astart_Deviatoric
	Astart = DL												// the starting direct lattice
	if (printIt)
		if ( stringmatch(StringFromList(0,GetRTStackInfo(0)),"IndexButtonProc"))
			printf "¥"
		endif
		printf "DeviatoricStrainRefine(%d,\"%d%d%d%d%d%d\"",pattern,aFit,bFit,cFit,alphaFit,betFit,gamFit
		if (!ParamIsDefault(coords))
			printf ",coords=%g",coords
		endif
		printf ")\r"
		if (Npatterns>1)
			printf "For pattern  %d  of  [0, %d]\r",pattern,Npatterns-1
		endif
		printf "starting reciprocal lattice,	a* = {%+.6f, %+.6f, %+.6f} (1/nm)\r",as0,as1,as2
		printf "							b* = {%+.6f, %+.6f, %+.6f}\r",bs0,bs1,bs2
		printf "							c* = {%+.6f, %+.6f, %+.6f}\r",cs0,cs1,cs2
		print " "
		printf "a = %.7g nm,   b = %.7g,   c = %.7g,   alpha = %.8g¡,   beta = %.8g¡,   gamma = %.8g¡,   Vc = %.7g (nm^3)\r",LC[0],LC[1],LC[2],LC[3],LC[4],LC[5],Vc
		printf "a/c = %.7g,   b/c = %.7g,   a/b = %.7g\r"LC[0]/LC[2],LC[1]/LC[2],LC[0]/LC[1]
	endif

	Make/N=3/O/D qhat_PeaksForStrain, qhat_PeaksForStrain2, qhat_latticeMismatch, ghat_latticeMismatch
	Wave qhat = qhat_PeaksForStrain
	Wave qhat2 = qhat_PeaksForStrain2

	Make/N=(N,14)/O/D PeaksForStrain=NaN
	SetDimLabel 1,0,h,PeaksForStrain ;			SetDimLabel 1,1,k,PeaksForStrain ;			SetDimLabel 1,2,l,PeaksForStrain
	SetDimLabel 1,3,fitQx,PeaksForStrain ;		SetDimLabel 1,4,fitQy,PeaksForStrain ;		SetDimLabel 1,5,fitQz,PeaksForStrain
	SetDimLabel 1,6,indexQx,PeaksForStrain ;	SetDimLabel 1,7,indexQy,PeaksForStrain ;	SetDimLabel 1,8,indexQz,PeaksForStrain
	SetDimLabel 1,9,fit_Intens,PeaksForStrain;	SetDimLabel 1,10,keV,PeaksForStrain
	SetDimLabel 1,11,pixelX,PeaksForStrain;	SetDimLabel 1,12,pixelY,PeaksForStrain;	SetDimLabel 1,13,angErr,PeaksForStrain
	PeaksForStrain[][0,2] = FullPeakIndexed[p][q+3][pattern]	// set the hkl
	String wnote = ReplaceNumberByKey("patternNum","",pattern,"=")
	wnote = ReplaceStringByKey("constrain",wnote,constrain,"=")
	wnote = ReplaceStringByKey("waveClass",wnote,"StrainPeakList","=")
	Note/K PeaksForStrain,wnote

	Variable Nj=DimSize(Qs,0), cosMax, jmax, dot			// find the measured peak that goes with each indexed peak (in Q^)
	Variable i,j,m
	for (i=0,m=0;i<N;i+=1)
		hkl[0] = PeaksForStrain[i][0]
		hkl[1] = PeaksForStrain[i][1]
		hkl[2] = PeaksForStrain[i][2]
		MatrixOp/O qhat = RLmeas x hkl						// predicted g^ from each measured hkl
		normalize(qhat)
		PeaksForStrain[i][6,8] = qhat[q-6]
		cosMax = -4
		jmax = -1
		for (j=0;j<Nj;j+=1)									// search fitted peaks to find the closest one
			qhat2 = Qs[j][p]
			dot = MatrixDot(qhat,qhat2)
			if (dot>cosMax)
				cosMax = dot
				jmax = j
			endif
		endfor
		if (jmax>=0)
			qhat2 = Qs[jmax][p]								// closest fitted peak to this hkl
			normalize(qhat2)
			PeaksForStrain[m][3,5] = qhat2[q-3]			// store fitted q^, and its intensity
			PeaksForStrain[m][9] = Qs[jmax][3]
			m += 1
		endif
	endfor
	N = m
	Variable rmsErr = latticeMismatchAll(PeaksForStrain,axis[0],axis[1],axis[2],LC)
	if (printIt)
//		printf "at start, rms error = %.5f¡\r",latticeMismatchAll(PeaksForStrain,axis[0],axis[1],axis[2],LC)		// 3 rotations and all 6 lattice constants
		printf "at start, rms error = %.5f¡\r",rmsErr		// 3 rotations and all 6 lattice constants
	endif

	KillWaves/Z qhat_PeaksForStrain2, RL0_latticeMismatch
	Redimension/N=(N,14) PeaksForStrain
	if (N<(Nfit+3))											// with fixed Vc, there are (Nfit-1) free parameters, and 3 rotations
		return ""
	endif

	Make/N=(Nfit+3)/O/D optimize_typXWave, optimize_xWave
	optimize_typXWave[0,2] = axis[p]
	for (m=3,i=0;i<6 && m<(Nfit+3);i+=1)
		if (str2num(constrain[i]))
			optimize_typXWave[m] = LC[i]
			m += 1
		endif
	endfor
	optimize_xWave = optimize_typXWave

	Variable/G printOnlyFirstStrainNaN=1
//	Optimize/Q/A=0/M={0,0}/R=optimize_typXWave/X=optimize_xWave/I=1000 $funcName, PeaksForStrain
	Optimize/Q/A=0/M={0,0}/R=optimize_typXWave/X=optimize_xWave/I=1000/Y=(rmsErr) $funcName, PeaksForStrain
	if (V_flag==791 && V_OptTermCode==6)				// too close to a critical point, change starting point and try again
		Variable del=1.001
		optimize_xWave[0]*=del
		if (numpnts(optimize_xWave)>3)
			optimize_xWave[3]*=del
			optimize_xWave[6]*=del
		endif
		printOnlyFirstStrainNaN=1
		Optimize/Q/A=0/M={0,0}/R=optimize_typXWave/X=optimize_xWave/I=1000/Y=(rmsErr) $funcName, PeaksForStrain
	elseif (V_flag==789 ||  V_OptTermCode==5)			// optimization got lost, try a different method
		printf "firtst try at optimization using a 'LIne Search' failed with V_flag=%g, V_OptTermCode=%g,   try again using 'More-Hebdon'\r",V_flag,V_OptTermCode
		optimize_xWave = optimize_typXWave
		printOnlyFirstStrainNaN=1
		Optimize/Q/A=0/M={2,0}/R=optimize_typXWave/X=optimize_xWave/I=1000/Y=(rmsErr) $funcName, PeaksForStrain
	endif
	KillVariables/Z printOnlyFirstStrainNaN

	if (V_flag)												// Optimization failed
		if (printIt)
			printf "Optimize failed with V_flag=%d,   V_OptTermCode=%d,   V_OptNumIters=%d,   V_OptNumFunctionCalls=%d\r",V_flag,V_OptTermCode,V_OptNumIters,V_OptNumFunctionCalls
		endif
		KillWaves/Z Ameas_Deviatoric, Astart_Deviatoric, totalTrans
		KillWaves/Z qhat_PeaksForStrain, FullPeakList2Qfile_Q, W_Cross
		KillWaves/Z W_OptGradient, optimize_xWave, optimize_typXWave
		KillWaves/Z ghat_latticeMismatch,qhat_latticeMismatch, axis_latticeMismatch, rho_latticeMismatch
		return ""
	endif

	axis = optimize_xWave[p]								// axis = {rx,ry,rz}
	fillLC(LC,constrain,optimize_xWave[3],optimize_xWave[4],optimize_xWave[5],optimize_xWave[6],optimize_xWave[7],optimize_xWave[8])
	KillWaves/Z W_OptGradient, optimize_xWave, optimize_typXWave
	forceLattice(LC,1,Vc=Vc)								// re-adjusts a,b,c to give correct Vc, (SpaceGroup=1 is triclinic, so no readjusting is done)
	rmsErr = latticeMismatchAll(PeaksForStrain,axis[0],axis[1],axis[2],LC)		// 3 rotations and all 6 lattice constants

	angle= norm(axis)
	rotationMatAboutAxis(axis,angle,rho)					// calculate rho, the rotation matrix from {rx,ry,rz}
	RLfromLatticeConstants(LC,DL,RL,Vc=Vc)				// calculate the reciprocal lattice from {a,b,c,alpha,bet,gam}
	MatrixOp/O RL = rho x RL								// rotate the reciprocal lattice by rho, this is the strained RL
	MatrixOp/O DL = rho x DL
	Ameas = DL												// Vcartesian = A x Vcell,   used to make the strain tensor

	String wName = UniqueName("sqrtMat",1,0)
	Duplicate/O Astart $wName
	Wave F = $wName											// holds the transformation from Astart to Ameas
	wName = UniqueName("sqrtMat",1,0)
	Duplicate/O Astart $wName
	Wave U = $wName										// holds U,  F = R x U

	MatrixOp/O F = Ameas x Inv(Astart)					// transformed un-strained to strained,  Ameas = F x Ao
	MatrixOp/O U = F^t x F									// U is symmetric, so F^t x F = U^t x R^t x R x U = U x R^-1 x R x U = U x U = U^2
	KillWaves/Z Ameas_Deviatoric, Astart_Deviatoric, F
	if (sqrtSymmetricMat(U))								// U starts as F^t x F, and is replaced by sqrt(U^2)
		KillWaves/Z U
		KillWaves/Z ghat_latticeMismatch,qhat_latticeMismatch, axis_latticeMismatch, rho_latticeMismatch
		KillWaves/Z hkl_latticeMismatch, RLmeas_latticeMismatch
		return ""
	endif
	MatrixOp/O epsilon = U - Identity(3)					// U = I + epsilon
	KillWaves/Z U

	Variable trace = MatrixTrace(epsilon)
	epsilon -= (p==q)*trace/3								// only deviatoric part, subtract trace/3 from diagonal
	trace = MatrixTrace(epsilon)							// re-set trace to be just the epsilon part, should be zero
	trace = abs(trace)<1e-15 ? 0 : trace					// set really small numbers to zero
	Wave epsilonAbs = $(microGeo#MakeUnique3x3Mat($""))
	epsilonAbs = abs(epsilon)
	WaveStats/Q epsilonAbs
	Variable SumAbsEpsilon=V_Sum
	KillWaves/Z epsilonAbs
	Variable epsilonvM										// von Mises strain
	epsilonvM = (epsilon[0][0]-epsilon[1][1])^2 + (epsilon[1][1]-epsilon[2][2])^2 + (epsilon[2][2]-epsilon[0][0])^2 
	epsilonvM += 6*( epsilon[0][1]^2 + epsilon[1][2]^2 + epsilon[2][0]^2 )
	epsilonvM = sqrt(epsilonvM/2)

	MatrixOp/O epsilonBL = rho x epsilon x Inv(rho)		// epsilon in beam line coordinates

	Wave Rlocal = $(microGeo#MakeUnique3x3Mat($""))// make a new 3x3 matrix,for the rotation from BL --> XHF coordinates (45¡ about X)
	Rlocal = 0
	Rlocal[0][0] = 1
	Rlocal[1,2][1,2] = 1/sqrt(2)
	Rlocal[2][1] *= -1
	MatrixOp/O epsilonXHF = Rlocal x epsilonBL x Inv(Rlocal)	// get epsilon in the XHF coordinate system
	Rlocal[2][1] *= 1/sqrt(2)								// now Rlocal transforms BL -->  "sample system" (uses outward surface normal for sample at 45¡)
	Rlocal[1][2] *= 1/sqrt(2)
	MatrixOp/O epsilonSample = Rlocal x epsilonBL x Inv(Rlocal)	// get epsilon in the Sample coordinate system
	KillWaves/Z Rlocal

	Variable startx,groupx, starty,groupy
	startx = NumberByKey("startx",note(FullPeakIndexed),"=")
	groupx = NumberByKey("groupx",note(FullPeakIndexed),"=")
	starty = NumberByKey("starty",note(FullPeakIndexed),"=")
	groupy = NumberByKey("groupy",note(FullPeakIndexed),"=")
	startx = numtype(startx) ? 1 : startx
	groupx = numtype(groupx) ? 1 : groupx
	starty = numtype(starty) ? 1 : starty
	groupy = numtype(groupy) ? 1 : groupy
	Variable ddLocal = geo.ddOffset + NumberByKey("CCDy",note(FullPeakIndexed),"=")	// get dd from the image (in wave note) if it is there
	geo.dd = ddLocal>0 ? ddLocal :  geo.dd
	Variable ycLocal = NumberByKey("yc",note(FullPeakIndexed),"=")						// get ycLocal from the image (in wave note) if it is there
	geo.ycent = ycLocal>0 ? ycLocal :  geo.ycent
	Variable/C pz
	Variable sine, d, keV
	Make/N=3/O/D fithat_PeaksForStrain
	Wave fithat = fithat_PeaksForStrain
	for (i=0;i<N;i+=1)
		qhat = PeaksForStrain[i][p+6]						// strained peak directions
		fithat = PeaksForStrain[i][p+3]						// fitted peak directions
		PeaksForStrain[i][13] = acos(limit(MatrixDot(qhat,fithat),-1,1))*180/PI
		sine = -qhat[3]										// sin(Bragg angle)
		d = dSpacingFromLatticeConstants(PeaksForStrain[i][0],PeaksForStrain[i][1],PeaksForStrain[i][2],LC[0],LC[1],LC[2],LC[3],LC[4],LC[5])
		keV = hc/(2*d*sine)
		PeaksForStrain[i][10][pattern] = keV				// energy of strained reflection (at strained position, not fitted position)
		IdealBeamLine2wenge(geo,qhat,qhat)				// transform qhat from Ideal Beam-line to Wenge, to get pixel on detector
		pz = q2pixel(geo,qhat)
		PeaksForStrain[i][11][pattern] = (real(pz)-(startx-1)-(groupx-1)/2)/groupx		// change to binned pixels
		PeaksForStrain[i][12][pattern] = (imag(pz)-(starty-1)-(groupy-1)/2)/groupy	// pixels are still zero based
	endfor

	Variable a=LC[0], b=LC[1], c=LC[2], alpha=LC[3], bet=LC[4], gam=LC[5]
	String str
	if (printIt)
		print " "
		printf "final reciprocal lattice,	a* = {%+.6f, %+.6f, %+.6f} (1/nm)\r",RL[0][0],RL[1][0],RL[2][0]
		printf "							b* = {%+.6f, %+.6f, %+.6f}\r",RL[0][1],RL[1][1],RL[2][1]
		printf "							c* = {%+.6f, %+.6f, %+.6f}\r",RL[0][2],RL[1][2],RL[2][2]
		printf "final lattice,	a = {%+.6f, %+.6f, %+.6f} (nm)\r",DL[0][0],DL[1][0],DL[2][0]
		printf "				b = {%+.6f, %+.6f, %+.6f}\r",DL[0][1],DL[1][1],DL[2][1]
		printf "				c = {%+.6f, %+.6f, %+.6f}\r",DL[0][2],DL[1][2],DL[2][2]
		print " "
		printf "a = %.7g nm,   b = %.7g,   c = %.7g,   alpha = %.8g¡,   beta = %.8g¡,   gamma = %.8g¡,   Vc = %.7g (nm^3)\r",a,b,c,alpha,bet,gam,Vc
		printf "a/c = %.7g,   b/c = %.7g,   a/b = %.7g\r",a/c,b/c,a/b
		printf "at end, rms error    = %.5f¡\r",rmsErr
		print " "
		sprintf str,"  \t*****  Tr(epsilonBL) = %.2g  *****",trace
		printf "epsilonCrystal =	{%+.6f, %+.6f, %+.6f}  \tvon Mises strain = %.3g%s\r",epsilon[0][0],epsilon[0][1],epsilon[0][2], epsilonvM, SelectString(abs(trace)>1e-12,"",str)
		printf "					{%+.6f, %+.6f, %+.6f}  \tmax{ | epsilon(ij) | } = %.3g\r",epsilon[1][0],epsilon[1][1],epsilon[1][2], V_max
		printf "					{%+.6f, %+.6f, %+.6f}  \tSum{ | epsilon | } = %.3g\r\r",epsilon[2][0],epsilon[2][1],epsilon[2][2], SumAbsEpsilon
		WaveStats/Q epsilonBL
		sprintf str,"  \t*****  Tr(epsilonBL) = %.2g  *****",MatrixTrace(epsilonBL)
		printf "epsilonBL =		{%+.6f, %+.6f, %+.6f}  \tmax{ | epsilon(ij) | } = %.3g\r",epsilonBL[0][0],epsilonBL[0][1],epsilonBL[0][2],max(abs(V_max),abs(V_min))
		printf "					{%+.6f, %+.6f, %+.6f}%s\r",epsilonBL[1][0],epsilonBL[1][1],epsilonBL[1][2], SelectString(abs(MatrixTrace(epsilonBL))>1e-12,"",str)
		printf "					{%+.6f, %+.6f, %+.6f}\r\r",epsilonBL[2][0],epsilonBL[2][1],epsilonBL[2][2]
		WaveStats/Q epsilonXHF
		sprintf str,"  \t*****  Tr(epsilonXHF) = %.2g  *****",MatrixTrace(epsilonXHF)
		printf "epsilonXHF =		{%+.6f, %+.6f, %+.6f}  \tmax{ | epsilon(ij) | } = %.3g\r",epsilonXHF[0][0],epsilonXHF[0][1],epsilonXHF[0][2],max(abs(V_max),abs(V_min))
		printf "					{%+.6f, %+.6f, %+.6f}%s\r",epsilonXHF[1][0],epsilonXHF[1][1],epsilonXHF[1][2], SelectString(abs(MatrixTrace(epsilonXHF))>1e-12,"",str)
		printf "					{%+.6f, %+.6f, %+.6f}\r\r",epsilonXHF[2][0],epsilonXHF[2][1],epsilonXHF[2][2]
		WaveStats/Q epsilonSample
		sprintf str,"  \t*****  Tr(epsilonSample) = %.2g  *****",MatrixTrace(epsilonSample)
		printf "epsilonSample =	{%+.6f, %+.6f, %+.6f}  \tmax{ | epsilon(ij) | } = %.3g\r",epsilonSample[0][0],epsilonSample[0][1],epsilonSample[0][2],max(abs(V_max),abs(V_min))
		printf "					{%+.6f, %+.6f, %+.6f}%s\r",epsilonSample[1][0],epsilonSample[1][1],epsilonSample[1][2], SelectString(abs(MatrixTrace(epsilonSample))>1e-12,"",str)
		printf "					{%+.6f, %+.6f, %+.6f}\r\r",epsilonSample[2][0],epsilonSample[2][1],epsilonSample[2][2]
	endif

	wnote = note(FullPeakIndexed)
	sprintf str, "{{%.8g,%.8g,%.8g}{%.8g,%.8g,%.8g}{%.8g,%.8g,%.8g}}",RL[0][0],RL[1][0],RL[2][0],  RL[0][1],RL[1][1],RL[2][1],  RL[0][2],RL[1][2],RL[2][2]
	wnote = ReplaceStringByKey("recip_lattice_refined",wnote,str,"=")
	sprintf str, "{%.8g %.8g %.8g %.8g %.8g %.8g}",a,b,c,alpha,bet,gam
	wnote = ReplaceStringByKey("lattice_constants_refined",wnote,str,"=")
	wnote = ReplaceStringByKey("waveClass",wnote,"strain_tensor","=")
	wnote = ReplaceNumberByKey("SumAbsEpsilon",wnote,SumAbsEpsilon,"=")
	Note/K epsilon, wnote

	wnote = note(PeaksForStrain)
	sprintf str, "{{%.8g,%.8g,%.8g}{%.8g,%.8g,%.8g}{%.8g,%.8g,%.8g}}",RL[0][0],RL[1][0],RL[2][0],  RL[0][1],RL[1][1],RL[2][1],  RL[0][2],RL[1][2],RL[2][2]
	wnote = ReplaceStringByKey("recip_lattice_refined",wnote,str,"=")
	sprintf str, "{%.8g %.8g %.8g %.8g %.8g %.8g}",a,b,c,alpha,bet,gam
	wnote = ReplaceStringByKey("lattice_constants_refined",wnote,str,"=")
	Note/K PeaksForStrain, wnote

	KillWaves/Z ghat_latticeMismatch,qhat_latticeMismatch, axis_latticeMismatch, rho_latticeMismatch
	KillWaves/Z hkl_latticeMismatch, RLmeas_latticeMismatch
	KillWaves/Z fithat_PeaksForStrain, optimize_LatticeConstantsWave
	KillWaves/Z RL_latticeMismatch, DL_latticeMismatch
	KillWaves/Z qhat_PeaksForStrain, FullPeakList2Qfile_Q, W_Cross
	KillWaves/Z M_Lower,M_Upper,W_LUPermutation,M_x
	return SelectString(coords-2,GetWavesDataFolder(epsilon,2),GetWavesDataFolder(epsilonXHF,2),GetWavesDataFolder(epsilonSample,2))
End
//
Static Function latticeMismatch4(PeaksForStrain,rx,ry,rz,f1)	// only fit rx,ry,rz, and one others
	Wave PeaksForStrain
	Variable rx,ry,rz
	Variable f1
	Wave LC=optimize_LatticeConstantsWave
	fillLC(LC,StringByKey("constrain",note(PeaksForStrain),"="),f1,NaN,NaN,NaN,NaN,NaN)
	return latticeMismatchAll(PeaksForStrain,rx,ry,rz,LC)
End
//
Static Function latticeMismatch5(PeaksForStrain,rx,ry,rz,f1,f2)	// only fit rx,ry,rz, and two others
	Wave PeaksForStrain
	Variable rx,ry,rz
	Variable f1,f2
	Wave LC=optimize_LatticeConstantsWave
	fillLC(LC,StringByKey("constrain",note(PeaksForStrain),"="),f1,f2,NaN,NaN,NaN,NaN)
	return latticeMismatchAll(PeaksForStrain,rx,ry,rz,LC)
End
//
Static Function latticeMismatch6(PeaksForStrain,rx,ry,rz,f1,f2,f3)	// only fit rx,ry,rz, and three others
	Wave PeaksForStrain
	Variable rx,ry,rz
	Variable f1,f2,f3
	Wave LC=optimize_LatticeConstantsWave
	fillLC(LC,StringByKey("constrain",note(PeaksForStrain),"="),f1,f2,f3,NaN,NaN,NaN)
	return latticeMismatchAll(PeaksForStrain,rx,ry,rz,LC)
End
//
Static Function latticeMismatch7(PeaksForStrain,rx,ry,rz,f1,f2,f3,f4)	// only fit rx,ry,rz, and four others
	Wave PeaksForStrain
	Variable rx,ry,rz
	Variable f1,f2,f3,f4
	Wave LC=optimize_LatticeConstantsWave
	fillLC(LC,StringByKey("constrain",note(PeaksForStrain),"="),f1,f2,f3,f4,NaN,NaN)
	return latticeMismatchAll(PeaksForStrain,rx,ry,rz,LC)
End
//
Static Function latticeMismatch8(PeaksForStrain,rx,ry,rz,f1,f2,f3,f4,f5)	// only fit rx,ry,rz, and five others
	Wave PeaksForStrain
	Variable rx,ry,rz
	Variable f1,f2,f3,f4,f5
	Wave LC=optimize_LatticeConstantsWave
	fillLC(LC,StringByKey("constrain",note(PeaksForStrain),"="),f1,f2,f3,f4,f5,NaN)
	return latticeMismatchAll(PeaksForStrain,rx,ry,rz,LC)
End
//
// This should never be called when optimizing Deviatoric strain, something has to be fixed
Static Function latticeMismatch9(PeaksForStrain,rx,ry,rz,f1,f2,f3,f4,f5,f6)	// only fit rx,ry,rz, and five others
	Wave PeaksForStrain
	Variable rx,ry,rz
	Variable f1,f2,f3,f4,f5,,f6
	Wave LC=optimize_LatticeConstantsWave
	fillLC(LC,StringByKey("constrain",note(PeaksForStrain),"="),f1,f2,f3,f4,f5,f6)
	return latticeMismatchAll(PeaksForStrain,rx,ry,rz,LC)
End
//
Static Function latticeMismatchAll(PeaksForStrain,rx,ry,rz,LC)		// 3 rotations and all 6 lattice constants
	Wave PeaksForStrain
	Variable rx,ry,rz
	Wave LC														//	lattice parameters:   a,b,c,alpha,bet,gam

	Wave qhat=qhat_latticeMismatch, ghat=ghat_latticeMismatch, axis=axis_latticeMismatch, rho=rho_latticeMismatch
	Wave RL=RL_latticeMismatch, hkl=hkl_latticeMismatch, LC=optimize_LatticeConstantsWave
	if (!WaveExists(qhat) || !WaveExists(ghat) || !WaveExists(axis) || !WaveExists(rho) || !WaveExists(RL) || !WaveExists(hkl) || !WaveExists(LC))
		Abort "some of the waves do not exist in latticeMismatchAll()"
	endif
	Variable Vc = NumberByKey("Vc",note(PeaksForStrain),"=")	// This is to be held constant
//	LC = {a,b,c,alpha,bet,gam}

	axis = {rx,ry,rz}
	Variable angle= norm(axis)
	rotationMatAboutAxis(axis,angle,rho)				// calculate rho, the rotation matrix from {rx,ry,rz}

	RLfromLatticeConstants(LC,$"",RL,Vc=Vc)			// calculate the reciprocal lattice from {a,b,c,alpha,bet,gam}
	MatrixOp/O RL = rho x RL							// rotate the reciprocal lattice by rho

	Variable N = DimSize(PeaksForStrain,0)
	Variable i,j, err=0
	Variable weight, SumWeight
	for (i=0,SumWeight=0; i<N; i+=1)
		qhat = PeaksForStrain[i][p+3]				// measured direction
		hkl[0] = PeaksForStrain[i][0]				// hkl of reflection
		hkl[1] = PeaksForStrain[i][1]
		hkl[2] = PeaksForStrain[i][2]
		MatrixOp/O ghat = RL x hkl
		normalize(ghat)
		PeaksForStrain[i][6,8] = ghat[q-6]		// save for later checking
//		weight = sqrt(sqrt(PeaksForStrain[i][9]))
		weight = 1
		err += (1-MatrixDot(qhat,ghat))*weight	// 1-dot is ~ 0.5*Ætheta^2,  using cos(theta) ~ 1-(theta^2)/2!
		SumWeight += weight
	endfor
	err =sqrt(err*2/SumWeight)					// change to rms
	err = numtype(err) ? Inf : err
	if (numtype(err) && NumVarOrDefault("printOnlyFirstStrainNaN",1))
		printf "err = %g,    LC={%g, %g, %g, %g, %g, %g}\r",err,LC[0],LC[1],LC[2],LC[3],LC[4],LC[5]
		Variable/G printOnlyFirstStrainNaN = 0
	endif
	return err
End
//
Static Function fillLC(LC,constrain,f1,f2,f3,f4,f5,f6)
	Wave LC										// lattice constants
	String constrain
	Variable f1,f2,f3,f4,f5,f6

	Variable i
	i = strsearch(constrain,"1",0)
	LC[i] = f1
	i = strsearch(constrain,"1",i+1)
	if (i<0)
		return 0
	endif
	LC[i] = f2
	i = strsearch(constrain,"1",i+1)
	if (i<0)
		return 0
	endif
	LC[i] = f3
	i = strsearch(constrain,"1",i+1)
	if (i<0)
		return 0
	endif
	LC[i] = f4
	i = strsearch(constrain,"1",i+1)
	if (i<0)
		return 0
	endif
	LC[i] = f5
	i = strsearch(constrain,"1",i+1)
	if (i<0)
		return 0
	endif
	LC[i] = f6
	return 0
End


Static Function forceLattice(LC,SG,[Vc])
	Wave LC
	Variable SG
	Variable Vc										// if present, then preserve Vc

	STRUCT crystalStructure xtal					// this sruct is set in this routine
	xtal.SpaceGroup = SG
	xtal.a = LC[0]
	xtal.b = LC[1]
	xtal.c = LC[2]
	xtal.alpha = LC[3]
	xtal.beta = LC[4]
	xtal.gam = LC[5]
	LatticeSym#ForceLatticeToStructure(xtal)		// also updates xtal.Vc
	LC[0] = xtal.a
	LC[1] = xtal.b
	LC[2] = xtal.c
	LC[3] = xtal.alpha
	LC[4] = xtal.beta
	LC[5] = xtal.gam
	if (!ParamIsDefault(Vc) && (Vc>0))			// passed a valid Vc, preserve it by adjusting a,b,c
		Variable factor = (Vc/xtal.Vc)^(1/3)
		LC[0,2] *= factor
	endif
End


// The choice of DL from lattice constants is same as section Vol B, 3.3 of the International Tables
Static Function RLfromLatticeConstants(LC,DL,RL[,Vc])	// set a*, b*, c* from lattice constants, option to preserve Vc
	Wave LC											// wave with 6 lattice constants
	Wave DL, RL										// direct and recip lattices
	Variable Vc										// rescale so volume of unit cell is this (if valid Vc)
	Vc = ParamIsDefault(Vc) ? NaN : Vc

	Variable a=LC[0], b=LC[1], c=LC[2]
	Variable alpha=LC[3], bet=LC[4], gam=LC[5]

	// calculate the reciprocal lattice from {a,b,c,alpha,bet,gam}
	Variable sa = sin((alpha)*PI/180), ca = cos((alpha)*PI/180)
	Variable cb = cos((bet)*PI/180), cg = cos((gam)*PI/180)
	Variable phi = sqrt(1.0 - ca*ca - cb*cb - cg*cg + 2*ca*cb*cg)	// = Vc/(a*b*c), dimensionless
	if (Vc>0)											// re-scale {a,b,c} to preserve volume of unit cell (deviatoric strain only)
		Variable scale = (Vc/(a*b*c*phi))^(1/3)	// a*b*c*phi = Vc
		a *= scale
		b *= scale
		c *= scale
	else
		Vc = a*b*c * phi								// no Vc passed, use computed Vc
	endif

	Variable pv = (2*PI) / (Vc)						// used for scaling reciprocal lattice
	Variable a0,a1,a2,  b0,b1,b2,  c0,c1,c2		//components of the direct lattice vectors
	Variable rhom = (a==b && b==c && alpha==bet && bet==gam && abs(bet-90)>0.1)	// rhombohedral
	if (rhom)											// special for rhombohedral axes
		Variable pp=a/3*sqrt(1+2*ca), qq=a/3*sqrt(1-ca)
		a0=pp+2*qq	; a1=pp-qq				; a2=pp-qq
		b0=pp-qq		; b1=pp+2*qq				; b2=pp-qq
		c0=pp-qq		; c1=pp-qq				; c2=pp+2*qq
	else												// for all others
		a0=a*phi/sa	; a1=a*(cg-ca*cb)/sa	; a2=a*cb
		b0=0			; b1=b*sa					; b2=b*ca
		c0=0			; c1=0						; c2=c
	endif
	if (WaveExists(DL))								// set direct lattice
		DL[0][0] = a0	;	DL[0][1] = b0	;	DL[0][2] = c0
		DL[1][0] = a1	;	DL[1][1] = b1	;	DL[1][2] = c1
		DL[2][0] = a2	;	DL[2][1] = b2	;	DL[2][2] = c2
	endif
	if (WaveExists(RL))								// set reciprocal lattice
		// a* = (b x c)*2¹/Vc				b* = (c x a)*2¹/Vc					c*=(a x b)*2¹/Vc
		RL[0][0]=(b1*c2-b2*c1)*pv	;	RL[0][1]=(c1*a2-c2*a1)*pv ;	RL[0][2]=(a1*b2-a2*b1)*pv
		RL[1][0]=(b2*c0-b0*c2)*pv ;	RL[1][1]=(c2*a0-c0*a2)*pv ;	RL[1][2]=(a2*b0-a0*b2)*pv
		RL[2][0]=(b0*c1-b1*c0)*pv ;	RL[2][1]=(c0*a1-c1*a0)*pv ;	RL[2][2]=(a0*b1-a1*b0)*pv
	endif
	return Vc
End


Static Function RL2latticeConstants(RL,LC)			// calc lattice constants from RL
	Wave RL											// recip lattice
	Wave LC											// wave to hold 6 lattice constants

	Make/N=3/O/D astar_latticeMismatch, bstar_latticeMismatch, cstar_latticeMismatch
	Wave astar=astar_latticeMismatch, bstar=bstar_latticeMismatch, cstar=cstar_latticeMismatch
	astar = RL[p][0]
	bstar = RL[p][1]
	cstar = RL[p][2]

	Cross bstar,cstar
	Wave W_Cross=W_Cross
	Variable pv = 2*PI/MatrixDot(astar,W_Cross)	// real lattice for un-strained material
	Make/N=3/O/D avec_Deviatoric,bvec_Deviatoric,cvec_Deviatoric
	Wave avec=avec_Deviatoric,bvec=bvec_Deviatoric,cvec=cvec_Deviatoric
	Cross bstar,cstar ;	avec = pv*W_Cross
	Cross cstar,astar ;	bvec = pv*W_Cross
	Cross astar,bstar ;	cvec = pv*W_Cross
	Cross bvec,cvec
	Variable Vc = MatrixDot(avec,W_Cross)			// volume of unstrained real unit cell
	LC[0] = norm(avec)							// the lattice parameters
	LC[1] = norm(bvec)
	LC[2] = norm(cvec)
	LC[3] = angleVec2Vec(bvec,cvec)
	LC[4] = angleVec2Vec(avec,cvec)
	LC[5] = angleVec2Vec(avec,bvec)
	KillWaves/Z astar_latticeMismatch, bstar_latticeMismatch, cstar_latticeMismatch
	KillWaves/Z avec_Deviatoric,bvec_Deviatoric,cvec_Deviatoric
	return Vc
End


Static Function StartStrainPanel(ipat,npat,aIn,bIn,cIn,alphaIn,betIn,gamIn)
	Variable &ipat, npat						// number of patterns, and initial pattern
	Variable &aIn,&bIn,&cIn, &alphaIn,&betIn,&gamIn

	NewDataFolder/O root:Packages:micro
	NewDataFolder/O root:Packages:micro:strainRefine
	String path = "root:Packages:micro:strainRefine:"
	if (exists(path+"aFit")!=2)				// ensure globals exist
		Variable/G $(path+"aFit")=1,$(path+"bFit")=1, $(path+"cFit")=0
		Variable/G $(path+"alphaFit")=1,$(path+"betFit")=1, $(path+"gamFit")=1
		Variable/G $(path+"patFit")=0
		String/G $(path+"patternNumberList")="0"
	endif
	NVAR aFit=$(path+"aFit"), bFit=$(path+"bFit"), cFit=$(path+"cFit")
	NVAR alphaFit=$(path+"alphaFit"), betFit=$(path+"betFit"), gamFit=$(path+"gamFit")
	cFit = (aFit&&bFit&&cFit) ? 0 : cFit	// cannot set all three lengths
	SVAR patternNumberList=$(path+"patternNumberList")

	ipat = !(ipat>=0) ? 0 : ipat
	npat = !(npat>ipat) ? ipat+1 : npat
	patternNumberList = expandRange("0-"+num2istr(max(0,npat-1)),";")

	NewPanel /W=(865,44,1020,177)/K=1
	PopupMenu popupPattern,pos={12,2},size={111,20},title="pattern #",fSize=12, disable=(npat==1 ? 2 : 0)
	PopupMenu popupPattern,mode=4,bodyWidth= 52,popvalue=num2istr(ipat),value=StrVarOrDefault("root:Packages:micro:strainRefine:patternNumberList","0")
	CheckBox check_a,pos={22,33},size={28,23},title="a",fSize=18,variable=$(path+"aFit"), proc=Indexing#UpdateStartStrainCheck
	CheckBox check_b,pos={22,55},size={29,23},title="b",fSize=18,variable=$(path+"bFit"), proc=Indexing#UpdateStartStrainCheck
	CheckBox check_c,pos={22,77},size={28,23},title="c",fSize=18,variable=$(path+"cFit"), proc=Indexing#UpdateStartStrainCheck
	CheckBox check_alpha,pos={92,35},size={30,18},title="a",font="Symbol",fSize=18
	CheckBox check_alpha,variable=$(path+"alphaFit"), proc=Indexing#UpdateStartStrainCheck
	CheckBox check_bet,pos={92,57},size={28,18},title="b",font="Symbol",fSize=18
	CheckBox check_bet,variable=$(path+"betFit"), proc=Indexing#UpdateStartStrainCheck
	CheckBox check_gam,pos={92,79},size={26,18},title="g",font="Symbol",fSize=18
	CheckBox check_gam,variable=$(path+"gamFit"), proc=Indexing#UpdateStartStrainCheck
	Button DoIt_button,pos={11,107},size={50,20},proc=Indexing#StartStrainButtonProc,title="Go"
	Button Cancel_button,pos={86,107},size={50,20},proc=Indexing#StartStrainButtonProc,title="Cancel"
	DoUpdate
	String win = StringFromList(0,WinList("*",";","WIN:64"))
	EnableDisableStartStrainCheck(win)
	PauseForUser $win

	aIn = NumVarOrDefault(path+"aFit",1)
	bIn = NumVarOrDefault(path+"bFit",1)
	cIn = NumVarOrDefault(path+"cFit",1)
	alphaIn = NumVarOrDefault(path+"alphaFit",1)
	betIn = NumVarOrDefault(path+"betFit",1)
	gamIn = NumVarOrDefault(path+"gamFit",1)
	ipat = NumVarOrDefault(path+"patFit",1)
	Variable Nfit = aIn+bIn+cIn+alphaIn+betIn+gamIn
	if (!(Nfit>0 && Nfit<6))
		ipat = -1
	endif
	aIn = aIn ? 1 : 0				// ensure only 0 or 1
	bIn = bIn ? 1 : 0
	cIn = cIn ? 1 : 0
	alphaIn = alphaIn ? 1 : 0
	betIn = betIn ? 1 : 0
	gamIn = gamIn ? 1 : 0
	return ipat					// <0 is error, >=0 is OK
End
//
Static Function StartStrainButtonProc(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	if (ba.eventCode != 2)		// not mouse up
		return 0
	endif
	Variable pattern
	Variable/G root:Packages:micro:strainRefine:patFit
	NVAR patFit= root:Packages:micro:strainRefine:patFit
	if (stringmatch(ba.ctrlName,"DoIt_button"))
		ControlInfo/W=$(ba.win) popupPattern
		patFit = V_flag ? str2num(S_Value): 0				// if no control present, then only pattern 0 allowed
	else
		patFit = -1				// flags a "Cancel"
	endif
	DoWindow/K $(ba.win)
	KillStrings/Z patternNumberList
	return 0
End
//
Static Function UpdateStartStrainCheck(cba) : CheckBoxControl
	STRUCT WMCheckboxAction &cba
	if (cba.eventCode != 2)
		return 0
	endif
	EnableDisableStartStrainCheck(cba.win)
	return 0
End
//
Static Function EnableDisableStartStrainCheck(win)
	String win
	DoUpdate

	String ctrlName, list="check_a;check_b;check_c;check_alpha;check_bet;check_gam"
	Variable i, n, nAll
	for (n=0,i=0; i<3; i+=1)				//	count how many boxes are checked only a,b,c
		ctrlName = StringFromList(i,list)
		ControlInfo/W=$win $ctrlName
		n += V_Value
	endfor
	if (n==3)								// un-check c box, for deviatoric, you cannot fit a,b, & c
		CheckBox check_c,value=0,disable=2
		n = 5
	elseif(n==2)							// disable the only un-checked box
		for (i=0;i<3;i+=1)
			ctrlName = StringFromList(i,list)
			ControlInfo/W=$win $ctrlName
			CheckBox $ctrlName,disable=(V_Value ? 0 : 2)
		endfor
	else	
		for (i=0;i<3;i+=1)					// enable all
			ctrlName = StringFromList(i,list)
			CheckBox $ctrlName,disable=0
		endfor
	endif

	for (nAll=n,i=3; i<6; i+=1)				//	count how many boxes are checked, all lattice constants
		ctrlName = StringFromList(i,list)
		ControlInfo/W=$win $ctrlName
		nAll += V_Value
	endfor
	Button DoIt_button,disable=((nAll<1 || n>2 || nAll>5) ? 2 : 0)
	return 0
End

Static Function/T FullPeakList2Qwave(FullPeakList) // convert fitted peaks to a Qlist+intens, peaks file for analysis be Euler
	Wave FullPeakList
	if (!WaveExists(FullPeakList))
		DoAlert 0, "input wave for FullPeakList2File() does not exists"
		return ""
	elseif (DimSize(FullPeakList,1)!=11)
		DoAlert 0, "the passed full peak list '"+NameOfWave(FullPeakList)+"' is not the right size"
		return ""
	endif
	Variable N=DimSize(FullPeakList,0)
	if (N<1)											// nothing to write
		return ""
	endif
	STRUCT microGeometry geo
	if (FillGeometryStructDefault(geo))			//fill the geometry structure with test values
		DoAlert 0, "no geometry structure found, did you forget to set it?"
		return ""
	endif

	Make/N=(N,4)/O/D FullPeakList2Qfile_Q=NaN // create temp waves for the q + intensity
	Wave Qs=FullPeakList2Qfile_Q					// holds Qhat and intensity
	Make/N=3/O/D FullPeakList2Qfile_qhat, FullPeakList2Qfile_qBL
	Wave qhat=FullPeakList2Qfile_qhat, qBL=FullPeakList2Qfile_qBL
	String wnote = note(FullPeakList)
	Variable startx,groupx, starty,groupy			// ROI of the original image
	startx = NumberByKey("startx",wnote,"=")
	groupx = NumberByKey("groupx",wnote,"=")
	starty = NumberByKey("starty",wnote,"=")
	groupy = NumberByKey("groupy",wnote,"=")
	startx = numtype(startx) ? 1 : startx
	groupx = numtype(groupx) ? 1 : groupx
	starty = numtype(starty) ? 1 : starty
	groupy = numtype(groupy) ? 1 : groupy
	Variable ddLocal = geo.ddOffset + NumberByKey("CCDy",wnote,"=")	// get dd from the image (in wave note) if it is there
	geo.dd = ddLocal>0 ? ddLocal :  geo.dd
	Variable ycLocal = NumberByKey("yc",wnote,"=")						// get ycLocal from the image (in wave note) if it is there
	geo.ycent = ycLocal>0 ? ycLocal :  geo.ycent
	Variable i, m, px,py
	for (i=0,m=0;i<N;i+=1)
		px = (startx-1) + groupx*FullPeakList[i][0] + (groupx-1)/2		// change to un-binned pixels
		py = (starty-1) + groupy*FullPeakList[i][1] + (groupy-1)/2		// pixels are still zero based
		pixel2q(geo,px,py,qhat)				// in Wenge-coord system
		if (norm(qhat)>0)							// check for a valid Q
			wenge2IdealBeamline(geo,qhat,qBL)		// convert to beam line coords, normalize and save
			normalize(qBL)
			Qs[m][0,2] = qBL[q]
			Qs[m][3] = max(FullPeakList[i][10],0) // save intensity too
			m += 1
		endif
	endfor
	KillWaves/Z FullPeakList2Qfile_qhat, FullPeakList2Qfile_qBL
	N = m
	Redimension/N=(N,4) FullPeakList2Qfile_Q
	ImageStats/M=1/G={0,N-1,3,3} FullPeakList2Qfile_Q
	FullPeakList2Qfile_Q[][3] /= V_max
	return GetWavesDataFolder(FullPeakList2Qfile_Q,2)
End
//
Function sqrtSymmetricMat(symmetricMat)			// computes sqrt(symmetricMat),  returns 1=error, 0=OK,   symmetricMat must be symmetric
	Wave symmetricMat									// Note, symmetricMat must be a "positive definite" matrix, or some of the eigen values will be negative

	MatrixEigenV /SYM/EVEC/Z symmetricMat
	if (V_flag)
		return 1
	endif
	Wave eigenValues = W_eigenValues, V = M_eigenVectors

//	// start of check
//	MatrixOp/O diagCheck = Inv(V) x symmetricMat x V
//	diagCheck = abs(diagCheck[p][q]) > 1e-15 ? 1 : 0
//	print diagCheck
//	// end of chek

	String wName = UniqueName("sqrtMat",1,0)
	Duplicate/O symmetricMat $wName
	Wave D = $wName								// holds the diagonal wave
	D = (p==q) ? sqrt(eigenValues[p]) : 0			// put sqrt(eigen values) on the diagonal, eigenValues better be >=0
	MatrixOp/O symmetricMat = V x D x Inv(V)	// transform back, putting root into original matrix
	KillWaves/Z W_eigenValues, M_eigenVectors, D
	return 0
End


Function VolumeOfLattice(a0,a1,a2,b0,b1,b2,c0,c1,c2)		// Volume = a .dot. (b .cross. c), the triple product
	Variable a0,a1,a2,b0,b1,b2,c0,c1,c2
	Variable Vc=0
	Vc += (b1*c2 - b2*c1)*a0
	Vc += (b2*c0 - b0*c2)*a1
	Vc += (b0*c1 - b1*c0)*a2
	return Vc
End

//	End of Index strain refinement
// =========================================================================
// =========================================================================



// =========================================================================
// =========================================================================
//	Start of Stereographic projection

Function StereoOfIndexedPattern(FullPeakIndexed,pattern)
	Wave FullPeakIndexed								// provides the reciprocal lattice
	Variable pattern										// in case more than one, default is 0

	if (exists("MakeStereo")!=6)
		print "Loading Stereographic package"
		Execute/P "INSERTINCLUDE  \"StereographicProjection\", version>=2.7";Execute/P "COMPILEPROCEDURES "
		return 1
	endif
	Variable printIt=0
	if (!WaveExists(FullPeakIndexed))
		String wList=WaveListClass("IndexedPeakList","*","")
		String wName=StringFromList(0,wList)
		if (ItemsInList(wList)!=1)						// only ask if more than one choice
			Prompt wName, "name of Full Peak Indexed list",popup,WaveListClass("IndexedPeakList","*","")
			DoPrompt "Full Peak Indexed",wName
			if (V_flag)
				return 1
			endif
		endif
		Wave FullPeakIndexed=$wName
		printIt = 1
	endif
	if (!WaveExists(FullPeakIndexed))
		return 1
	endif
	Variable Npatterns=DimSize(FullPeakIndexed,2)
	pattern = (Npatterns==1) ? 0 : pattern
	pattern = limit(pattern,0,Npatterns-1)
	if (numtype(pattern))
		Prompt pattern, "pattern number",popup,expandRange("0-"+num2istr(Npatterns-1),";")
		DoPrompt "pattern number",pattern
		if (V_flag)
			return 1
		endif
		printIt = 1
	endif
	printIt = printIt || ((ItemsInList(GetRTStackInfo(0))<2) || stringmatch(StringFromList(0,GetRTStackInfo(0)),"IndexButtonProc"))

	Make/N=3/O/D hkl_StereoOfIndexedPattern, vec_StereoOfIndexedPattern
	Make/N=(3,3)/O/D RL_StereoOfIndexedPattern
	Wave hkl=hkl_StereoOfIndexedPattern, vec=vec_StereoOfIndexedPattern, RL=RL_StereoOfIndexedPattern
	Variable as0,as1,as2,bs0,bs1,bs2,cs0,cs1,cs2
	sscanf StringByKey("recip_lattice"+num2istr(pattern),note(FullPeakIndexed),"="), "{{%g,%g,%g}{%g,%g,%g}{%g,%g,%g}}",as0,as1,as2,bs0,bs1,bs2,cs0,cs1,cs2
	RL[0][0] = as0 ;		RL[0][1] = bs0 ;		RL[0][2] = cs0		// the original measured RL
	RL[1][0] = as1 ;		RL[1][1] = bs1 ;		RL[1][2] = cs1
	RL[2][0] = as2 ;		RL[2][1] = bs2 ;		RL[2][2] = cs2
	vec = {0,1,-1}
	MatrixOp/O hkl = Inv(RL) x vec
	normalize(hkl)
	hkl = round(hkl*10000)/1000
	Variable h=hkl[0], k=hkl[1], l=hkl[2]
	vec = {1,0,0}
	MatrixOp/O hkl = Inv(RL) x vec
	normalize(hkl)
	hkl = round(hkl*10000)/1000
	if (printIt)
		if (stringmatch(StringFromList(0,GetRTStackInfo(0)),"IndexButtonProc"))
			printf "¥"
		endif
		printf "StereoOfIndexedPattern(%s,%d)\r",NameOfWave(FullPeakIndexed),pattern
		printf "Making Stereographic projection with pole at (%g, %g, %g),  and  x toward (%g, %g, %g)\r",h,k,l,round(hkl[0]*1000)/100, round(hkl[1]*1000)/100, round(hkl[2]*1000)/100
	endif
	FUNCREF MakeStereoProto func = $"MakeStereo"
	func(h,k,l,  NaN,9, 0,hklPerp=hkl)		// MakeStereo(h,k,l,  NaN,9, 0,hklPerp=hkl)
	KillWaves/Z hkl_StereoOfIndexedPattern, vec_StereoOfIndexedPattern, RL_StereoOfIndexedPattern
	return 0
End
Function MakeStereoProto(Hz,Kz,Lz,hklmax,fntSize,phi,[Qmax,hklPerp, WulffStepIn,WulffPhiIn])
	Variable Hz,Kz,Lz			// hkl of the pole
	Variable hklmax
	Variable fntSize
	Variable phi				// aximuthal angle (degree)
	Variable Qmax				// maximum Q (1/nm)
	Wave hklPerp				// optional hkl giving the perpendicular direction long=0
	Variable WulffStepIn, WulffPhiIn	// for Wulff Net, use WulffStep=0 for no Wulff Net
End

//	End of Stereographic projection
// =========================================================================
// =========================================================================



// =========================================================================
// =========================================================================
//	Start of Index set panel

Function/T FillIndexParametersPanel(strStruct,hostWin,left,top)
	String strStruct									// optional passed value of xtal structure, this is used if passed
	String hostWin										// name of home window
	Variable left, top									// offsets from the left and top

	NewDataFolder/O root:Packages:micro				// ensure that the needed data folders exist
	NewDataFolder/O root:Packages:micro:Index
	String/G root:Packages:Lattices:PanelValues:desc
	Variable/G root:Packages:micro:Index:h_e
	Variable/G root:Packages:micro:Index:k_e
	Variable/G root:Packages:micro:Index:l_e
	NVAR h=root:Packages:micro:Index:h_e
	NVAR k=root:Packages:micro:Index:k_e
	NVAR l=root:Packages:micro:Index:l_e

	SetWindow kwTopWin,userdata(IndexPanelName)=hostWin+"#IndexPanel"
	NewPanel/K=1/W=(left,top,left+221,top+445+30)/HOST=$hostWin
	ModifyPanel frameStyle=0, frameInset=0
	RenameWindow #,IndexPanel

	Button buttonLoadImage,pos={17,5},size={80,20},proc=IndexButtonProc,title="Load Image"
	Button buttonLoadImage,help={"Load an Image from the file, and then display the image."}
	Button buttonViewImage,pos={111,5},size={80,20},proc=IndexButtonProc,title="View SPE"
	Button buttonViewImage,help={"Display an existing Image that is already loaded"}

	Button buttonBkgRemove,pos={29,40},size={160,20},proc=IndexButtonProc,title="Remove Background"
	Button buttonBkgRemove,help={"remove background from an image, do not use on depth sotred images"}
	Button buttonFitPeaks,pos={29,65},size={160,20},proc=IndexButtonProc,title="Fit Peaks"
	Button buttonFitPeaks,help={"Identify peaks in an Image, does not do any background removal"}
	Button buttonPeaksFromBoxes,pos={29,90},size={160,20},proc=IndexButtonProc,title="set peaks from boxes"
	Button buttonPeaksFromBoxes,help={"reset list of peaks for fitting to the boxes on an image plot"}

	Button buttonIndex,pos={29,125},size={160,20},proc=IndexButtonProc,title="Index & Display"
	Button buttonIndex,help={"Index the identified peaks, using the defined lattice, and then show it all on a plot"}
	PopupMenu popuphklTags,pos={69,150},size={76,20},proc=hklTagsPopMenuProc,title="hkl Tags..."
	PopupMenu popuphklTags,help={"re-draw or clear the hkl tags from a plot"}
	PopupMenu popuphklTags,fSize=14
	PopupMenu popuphklTags,mode=0,value= #"\"Re-Draw hkl tags;Clear hkl tags off Graph\""
	Button buttonStrainRefine,pos={29,175},size={160,20},proc=IndexButtonProc,title="Strain Refine"
	Button buttonStrainRefine,help={"after indexing, use this to refine the lattice and find the deviatoric strain"}
	if (exists("MakeStereo")==6)
		Button buttonStereographic,pos={29,205},size={160,20},proc=IndexButtonProc,title="Stereographic..."
	else
		Button buttonStereographic,pos={29,205},size={160,20},proc=IndexButtonProc,title="Load Stereographic"
	endif
	Button buttonStereographic,help={"Make a Stereographic projection with the same orientation as the indexed pattern"}
	Button buttonReportIndex,pos={29,205+30},size={160,20},proc=IndexButtonProc,title="Make Report Sheet"
	Button buttonReportIndex,help={"if the top window is an images showing the hkl's, this makes a layout with reflections for printing"}

	Button buttonCalcEnergyIndex,pos={29,235+30},size={160,20},proc=IndexButtonProc,title="Calc Energy of (hkl)"
	Button buttonCalcEnergyIndex,help={"calculate the energy of an hkl using the current indexing"}
	SetVariable h_IndexEneVar,pos={27,257+30},size={50,15},proc=CalcE_SetVarProc,title="H"
	SetVariable h_IndexEneVar,value= root:Packages:micro:Index:h_e
	SetVariable k_IndexEneVar,pos={91,257+30},size={50,15},proc=CalcE_SetVarProc,title="K"
	SetVariable k_IndexEneVar,value= root:Packages:micro:Index:k_e
	SetVariable L_IndexEneVar,pos={155,257+30},size={50,15},proc=CalcE_SetVarProc,title="L"
	SetVariable L_IndexEneVar,value= root:Packages:micro:Index:l_e

	PopupMenu popupTables,pos={62,290+30},size={66,20},proc=TablesPopMenuProc,title="Display Tables..."
	PopupMenu popupTables,help={"Show table of identified peak positions indexed peaks, or the peaks for strain refinement"}
	PopupMenu popupTables,fSize=14,mode=0,value= #"\"Fitted Peaks;Indexed Peaks;Strain Peaks\""

	if (exists("IndexLots")==6)
		SetDrawLayer UserBack
		SetDrawEnv linethick= 2
		DrawLine 0,330+30,221,330+30
		Button buttonIndexLots,pos={29,340+30},size={160,20},proc=IndexButtonProc,title="Index Lots"
		Button buttonIndexLots,help={"index lots of images"}
		Button buttonWrite3dFile,pos={39,365+30},size={160,20},proc=IndexButtonProc,title="Write 3d Orients file"
		Button buttonWrite3dFile,help={"Write 3d Orientations file"}
		Button buttonFindClosestPoint,pos={39,390+30},size={160,20},proc=IndexButtonProc,title="find closest point"
		Button buttonFindClosestPoint,help={"find point in original indexing list that is closest to the cursor"}
		Button buttonIndexingListTable,pos={39,415+30},size={160,20},proc=IndexButtonProc,title="Indexing List Table"
		Button buttonIndexingListTable,help={"Make a table of the indexing list"}
	else
		Button buttonInitIndexLots,pos={29,340+30},size={150,50},title="Init IndexLots\rpackage",proc=LoadPackageButtonProc
		Button buttonInitIndexLots,help={"Load Package for Indexing Lots of Images"}
	endif

	EnableDisableIndexControls(hostWin+"#IndexPanel")
	return "#IndexPanel"
End
//
Function EnableDisableIndexControls(win)				// here to enable/disable
	String win												// window (or sub window) to use
	Variable d

	Button buttonLoadImage,win=$win,disable=0		// always OK to load an image

	d = strlen(WaveListClass("spe*","*","DIMS:2"))<1 ? 2 : 0
	Button buttonViewImage,win=$win,disable= d

	d = strlen(WaveListClass("speImage","*","DIMS:2"))<1 ? 2 : 0
	Button buttonBkgRemove,win=$win,disable= d

	d = strlen(WaveListClass("speImage*","*","DIMS:2"))<1 ? 2 : 0
	Button buttonFitPeaks,win=$win,disable= d

	d = strlen(StrVarOrDefault("pkList","")) ? 0 : 2
	Button buttonPeaksFromBoxes,win=$win,disable= d

	d = strlen(WaveListClass("FittedPeakList","*","DIMS:2,MAXCOLS:11,MINCOLS:11"))<1 ? 2 : 0
	Button buttonIndex,win=$win,disable= d

	d = (strlen(WaveListClass("FittedPeakList","*","DIMS:2,MAXCOLS:11,MINCOLS:11")) && strlen(WaveListClass("IndexedPeakList","*",""))) ? 0 : 2
	Button buttonStrainRefine,win=$win,disable= d
	Button buttonStereographic,win=$win,disable= d

	d = strlen(WaveListClass("IndexedPeakList","*",""))<1 ? 2 : 0
	PopupMenu popuphklTags,win=$win,disable= d
	Button buttonReportIndex,win=$win,disable= d
	Button buttonCalcEnergyIndex,win=$win,disable= d
	if (d)
		SetVariable h_IndexEneVar,win=$win,noedit=1,frame=0,disable= d
		SetVariable k_IndexEneVar,win=$win,noedit=1,frame=0,disable= d
		SetVariable L_IndexEneVar,win=$win,noedit=1,frame=0,disable= d
	else
		SetVariable h_IndexEneVar,win=$win,noedit=0,frame=2,disable= d
		SetVariable k_IndexEneVar,win=$win,noedit=0,frame=2,disable= d
		SetVariable L_IndexEneVar,win=$win,noedit=0,frame=2,disable= d
	endif

	d = strlen(WaveListClass("FittedPeakList;IndexedPeakList","*",""))<1 ? 2 : 0
	PopupMenu popupTables,win=$win,disable= d

	if (exists("IndexLots")==6)
		Button buttonIndexLots,win=$win,disable= 0					// always OK to do this

		d = strlen(WaveListClass("indexationResultList","*","TEXT:1"))<1 ? 2 : 0
		Button buttonWrite3dFile,win=$win,disable= d

		d = strlen(WaveListClass("OrientationSliceWave*","*",""))<1 ? 2 : 0
		Button buttonFindClosestPoint,win=$win,disable= d

		d = strlen(WaveListClass("indexationResultList*","*","TEXT:1"))<1 ? 2 : 0
		Button buttonIndexingListTable,win=$win,disable= d
	endif
End


Function IndexButtonProc(B_Struct) : ButtonControl
	STRUCT WMButtonAction &B_Struct
	if (B_Struct.eventCode != 2)
		return 0
	endif
	String ctrlName=B_Struct.ctrlName

	if (stringmatch(ctrlName,"buttonLoadImage"))
		Wave image = $LoadWinViewFile("")
		if (WaveExists(image))
			String wnote=note(image), bkgFile=StringByKey("bkgFile",wnote,"=")
			Variable xdim=NumberByKey("xdim",wnote,"="), ydim=NumberByKey("ydim",wnote,"=")
			SVAR S_fileName=S_fileName
			printf "for file '"+S_fileName+"'"
			if (strlen(bkgFile)>0)
				printf ",             background file = '%s'",  bkgFile
			endif
			printf "\r"
			printf "total length = %d x %d  = %d points\r", xdim,ydim,xdim*ydim
			print "number type is  '"+WinViewFileTypeString(NumberByKey("numType", wnote,"="))+"'"
			print "Created a 2-d wave    '"+GetWavesDataFolder(image,2)+"'"
			DoAlert 1, "Display this image"
			if (V_Flag==1)
				Graph_imageMake(image,1)
			endif
		endif
	elseif (stringmatch(ctrlName,"buttonViewImage") && strlen(WaveListClass("spe*","*","DIMS:2")))
		Graph_imageMake($"",1)
	elseif (stringmatch(ctrlName,"buttonBkgRemove") && strlen(WaveListClass("speImage","*","DIMS:2")))
		RemoveBackgroundImage($"",NaN)
	elseif (stringmatch(ctrlName,"buttonFitPeaks") && strlen(WaveListClass("speImage*","*","DIMS:2")))
		// FitPeaks($"",NaN,NaN)
		//		FunctionList("FitPeak*",";" , "WIN:[Indexing]")
		String fitPeakFunc="FitPeaksWithSeedFill"		// "FitPeaksStepWise"
		prompt fitPeakFunc,"peak fitting method",popup,"FitPeaksWithSeedFill;FitPeaks;FitPeaksNew;FitPeaksStepWise"
		DoPrompt "Peak Fit",fitPeakFunc
		if (!V_flag)
			if (stringmatch(fitPeakFunc,"FitPeaks"))
				FitPeaks($"",NaN,NaN)
			elseif (stringmatch(fitPeakFunc,"FitPeaksNew"))
				FitPeaksNew($"",NaN,NaN)
			elseif (stringmatch(fitPeakFunc,"FitPeaksStepWise"))
				FitPeaksStepWise($"",NaN,NaN,NaN,NaN)
			elseif (stringmatch(fitPeakFunc,"FitPeaksWithSeedFill"))
				FitPeaksWithSeedFill($"",NaN,NaN,NaN,NaN)
			endif
		endif
	elseif (stringmatch(ctrlName,"buttonPeaksFromBoxes") && strlen(StrVarOrDefault("pkList","")))
		SVAR pkList=pkList
		setFittedPeaksFromList($"",NaN,pkList)
	elseif (stringmatch(ctrlName,"buttonIndex") && strlen(WaveListClass("FittedPeakList","*","DIMS:2,MAXCOLS:11,MINCOLS:11")))
		IndexAndDisplay($"",NaN,NaN,NaN,NaN,NaN,NaN,NaN)
	elseif (stringmatch(ctrlName,"buttonStrainRefine"))
		DeviatoricStrainRefine(NaN,"")
	elseif (stringmatch(ctrlName,"buttonStereographic"))
		StereoOfIndexedPattern($"",NaN)
		if (exists("MakeStereo")!=6)
			Button buttonStereographic,win=$(B_Struct.win),title="Stereographic..."
		endif
	elseif (stringmatch(ctrlName,"buttonReportIndex"))
		MakeIndexingReportSheet()
	elseif (stringmatch(ctrlName,"buttonCalcEnergyIndex"))
		NVAR h=root:Packages:micro:Index:h_e
		NVAR k=root:Packages:micro:Index:k_e
		NVAR l=root:Packages:micro:Index:l_e
		EnergyOfhkl($"",NaN,h,k,l)
#if (Exists("IndexLots")==6)
	elseif (stringmatch(ctrlName,"buttonIndexLots") )
		IndexLots("","","","",NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN)
	elseif (stringmatch(ctrlName,"buttonWrite3dFile") && strlen(WaveListClass("indexationResultList","*","TEXT:1")))
		Write3dOrientsFile($"","")
	elseif (stringmatch(ctrlName,"buttonFindClosestPoint") && strlen(WaveListClass("OrientationSliceWave*","*","")))
		findClosestPoint(NaN,NaN)
	elseif (stringmatch(ctrlName,"buttonIndexingListTable") && strlen(WaveListClass("indexationResultList*","*","TEXT:1")))
		TableOfIndexingList($"")
#endif
	endif
	EnableDisableIndexControls(GetUserData("microPanel","","IndexPanelName"))
End
//
Function hklTagsPopMenuProc(ctrlName,popNum,popStr) : PopupMenuControl
	String ctrlName
	Variable popNum
	String popStr
	if (strlen(WaveListClass("IndexedPeakList","*",""))<1)
		return 1
	endif
	if (stringmatch(popStr,"Re-Draw hkl tags"))
		DisplayResultOfIndexing($"",NaN)
	elseif (stringmatch(popStr,"Clear hkl tags off Graph"))
		ClearhklOffGraph("")
	endif
	EnableDisableIndexControls(GetUserData("microPanel","","IndexPanelName"))
End
//
Function CalcE_SetVarProc(ctrlName,varNum,varStr,varName) : SetVariableControl
	String ctrlName
	Variable varNum
	String varStr
	String varName
	if (strlen(WaveListClass("IndexedPeakList","*",""))<1)
		return 1
	endif
	NVAR h=root:Packages:micro:Index:h_e
	NVAR k=root:Packages:micro:Index:k_e
	NVAR l=root:Packages:micro:Index:l_e
	// EnergyOfhkl($"",NaN,h,k,l)
	EnableDisableIndexControls(GetUserData("microPanel","","IndexPanelName"))
End
//
Function TablesPopMenuProc(ctrlName,popNum,popStr) : PopupMenuControl
	String ctrlName
	Variable popNum
	String popStr
	if (strlen(WaveListClass("FittedPeakList;IndexedPeakList","*",""))<1)
		return 1
	endif
	if (stringmatch(popStr,"Fitted Peaks"))
		TableFullPeakList($"")
	elseif (stringmatch(popStr,"Indexed Peaks"))
		TableFullPeakIndexed($"")
	elseif (stringmatch(popStr,"Strain Peaks"))
		TableFullPeakStrain($"")
	endif
End


//	End of Index set panel
// =========================================================================
// =========================================================================


Function InitIndexingPackage()					// used to initialize this package
	 InitLatticeSymPackage()					// init Lattice Sym package
	NewDataFolder/O root:Packages:geometry
End