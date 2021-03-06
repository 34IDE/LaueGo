#pragma rtGlobals=3		// Use modern global access method.
#pragma TextEncoding = "MacRoman"
#pragma ModuleName=JZTutil
#pragma IgorVersion = 6.11
#pragma version = 4.90
// #pragma hide = 1

Menu "Graph"
	"Append Multiple Graphs to a Layout",AppendGraph2LayoutN(NaN,"","")
End
Menu "Layout"
	"Append Multiple Graphs to a Layout",AppendGraph2LayoutN(NaN,"","")
End
Menu "Data"
	"Generic Wave Note Info", GenericWaveNoteInfo($"","")
End


#if (IgorVersion()>7)
	// Igor 7 use UNICODE
	strConstant DEGREESIGN = "\xC2\xB0"		// UTF8, DEGREE SIGN
	strConstant BULLET = "\xE2\x80\xA2"
	strConstant ARING = "\xC3\x85"				// Aring, Angstrom sign
	strConstant BCHAR = "\xE2\x80\x94"		// EM DASH
	strConstant GDELTA = "\xCE\x94"			// UTF8, Greek DELTA
	strConstant HORIZ_ELLIPSIS = "\xE2\x80\xA6"	// UTF8, Horizontal Ellipsis
	strConstant Gmu = "\xCE\xBC"				// UTF8, Greek mu
	strConstant PLUSMINUS = "\xC2\xB1"		// UTF8, plus-minus sign
	strConstant INTEGRAL_SIGN = "\xE2\x88\xAB"		// UTF8, integral sign
	strConstant ANGLE_SIGN = "\xE2\x88\xA0"		// angle sign
#elif StringMatch(IgorInfo(2),"Windows")
	// Igor 6 Windows
	strConstant BULLET = "\7"				// Bullet, MS Alt-7
	strConstant HORIZ_ELLIPSIS = "\205"	// ellipsis, was: "\0133"  MS Alt-0133
	strConstant DEGREESIGN = "\260"		// Degree sign, was: "\370"  MS Alt-248d=o370
	strConstant ARING = "\305"				// Angstrom sign, was: "\217"  MS Alt-143d=217o
	strConstant BCHAR = "\257"				// EM DASH, MS Alt-257, or "\226" or "\227"
	strConstant GDELTA = "\916"				// ???
	strConstant PLUSMINUS = "\261"			// MS plus-minus sign, was: "\361"  MS Alt-241, 241d=361
	strConstant Gmu = "\265"					// Greek mu, was: "\346"  MS Alt-230d=346o
	strConstant INTEGRAL_SIGN = "\222B"	// Integral sign, MS Alt-222B
	strConstant ANGLE_SIGN = "angle"		// angle sign
#else
	// Igor 6 Mac
	strConstant BULLET = "\245"				// Mac option-8
	strConstant HORIZ_ELLIPSIS = "\311"	// Mac option-; xC9, o311, d301
	strConstant DEGREESIGN = "\241"		// option-shift-8
	strConstant ARING = "\201"				// Angstrom sign, option-shift-A
	strConstant BCHAR = "\321"				// EM DASH
	strConstant PLUSMINUS = "\261"			// Mac option-+, plus-minus sign
	strConstant GDELTA = "\306"				// Mac option-j, xC6, o306, d198
	strConstant Gmu = "\265"					// Mac option-m, Greek mu
	strConstant INTEGRAL_SIGN = "\272"	// option-b
	strConstant ANGLE_SIGN = "angle"		// angle sign
#endif

StrConstant ELEMENT_Symbols = "H;He;Li;Be;B;C;N;O;F;Ne;Na;Mg;Al;Si;P;S;Cl;Ar;K;Ca;Sc;Ti;V;Cr;Mn;Fe;Co;Ni;Cu;Zn;Ga;Ge;As;Se;Br;Kr;Rb;Sr;Y;Zr;Nb;Mo;Tc;Ru;Rh;Pd;Ag;Cd;In;Sn;Sb;Te;I;Xe;Cs;Ba;La;Ce;Pr;Nd;Pm;Sm;Eu;Gd;Tb;Dy;Ho;Er;Tm;Yb;Lu;Hf;Ta;W;Re;Os;Ir;Pt;Au;Hg;Tl;Pb;Bi;Po;At;Rn;Fr;Ra;Ac;Th;Pa;U;Np;Pu;Am;Cm;Bk;Cf;Es;Fm;Md;No;Lr;Rf;Db;Sg;Bh;Hs;Mt;Ds;Rg;Cn;Nh;Fl;Mc;Lv;Ts;Og"
StrConstant MonthNamesFull = "January;February;March;April;May;June;July;August;September;October;November;December"
StrConstant MonthNamesShort = "Jan;Feb;Mar;Apr;May;Jun;Jul;Aug;Sep;Oct;Nov;Dec"
StrConstant DayNamesFull = "Sunday;Monday;Tuesday;Wednesday;Thursday;Friday;Saturday"
StrConstant DayNamesShort = "Sun;Mon;Tue;Wed;Thu;Fri;Sat"
Static Constant Smallest32bitFloat = 1.40129846432482e-45			// see DefaultZeroThresh(ww) below for use and finding
Static Constant Smallest64bitFloat = 4.94065645841247e-324
Static Constant maxIgorWaveNameLen = 31
StrConstant WhiteSpaceChars = " \t\r\n"

#if (IgorVersion()<7)
Constant GIZMO_WIN_TYPE = 13			// numbers for Igor 6 and under
Constant GIZMO_WIN_BIT = 4096
#else
Constant GIZMO_WIN_TYPE = 17			// numbers for Igor 7
Constant GIZMO_WIN_BIT = 65536	
#endif
#if StringMatch(IgorInfo(2)"Windows")
strConstant GenevaEquivFont = "Verdana"			// on windows, this is the Geneva equivalent
#else
strConstant GenevaEquivFont = "Geneva"			// Geneva is for Mac
#endif
StrConstant XMLfilters = "XML Files (*.xml,*.txt):.xml,.txt;All Files:.*;"
StrConstant XMLfiltersStrict = "XML Files (*.xml):.xml,;All Files:.*;"


// Sections:
//	1	function for optionally showing menus, e.g. MenuItemIfWaveClassExists(), and others
//	2	macros to put labels on lower left/right of plots to identify where they came from
//	3	tool for putting multiple graph on one page
//	4	String ranges, deals with "1,4,5-20",  for handling a non-consecutive range of integers, This one is good
//	5	Progress panels
//	6	Generic XML support
//			XMLNodeList()					returns a list of node tag at top most level in buf
//			XMLtagContents()				returns contents of an element
//			XMLremoveComments()			remove all xml comments from str, returns string without comments
//			XMLattibutes2KeyList()		return a list with all of the attribute value pairs for xmltag
//			XMLtagContents2List()		reads a tag contensts and converts it to a list, kind of specialized
//	7	WaveListClass() & WaveInClass(), used for handling the "waveClass=ccccc;" in wave notes
//		  also AddClassToWaveNote(), RemoveClassFromWaveNote(), ExcludeWavesInClass(), IncludeOnlyWavesInClass()
//		  IncludeOnlyWavesInClass(), removes waves from the list if they are not of correct class
//		  WindowsWithClass(), 			return list of windows that display a wave of appropriate waveClass
//		  FoldersWithWaveClass(), returns list of sub-folders with wave of a certain class
//		  TraceNamesInClass(), like TraceNameList(), but limit result to a list of wave classes
//	8	Contains lots of utility stuff
//		RecompileAllProcedures(), FORCE ALL procedures to recompile
//		FixUpFolderName(fldr), return foldr name with the needed ":", so fldr+"wavename" is valid
//		AddEndingToWaveName(wName,waveNameEnd), add an ending to wName, without exceeding max permitted name length
//		WavesWithMatchingKeyVals(), further filter a list of waves, look for those with matching key=value pairs
//		keyInList(), MergeKeywordLists(), & keysInList() "key=value" list utilities
//		joinLists(a,b,[sep]) returns lists a & b joined together in one list
//		intersectionOfLists(a,b,[sep]) returns intersection of two lists a & b
//		RemoveDuplicatesFromList(list,[listSepStr,matchCase]), return list with duplicate items removed
//		OnlyWavesThatAreDisplayed(), removes waves that are not displayed from a list of wave
//		NextLineInBuf(), returns next line in buf (which has <nl> separators) each time it is called.
//		AxisLabelFromGraph(), gets the axis label
//		WindowsWithWave(w,flag), return list of all windows with w, flag (a mask): graph=1, table=2, gizmo=4
//		isWaveOnGizmo(gName,ww),  returns True when ww is in object list of Gizmo
//		DisplayTableOfWave(...), display table taking into account any col/row labels
//		getListOfTypesInFile(), returns file type (using the $filetype) in my old standard files (trying to start only using xml)
//		DrawMarker(), draw a marker
//		xy2saturatedColors(), computes saturated colors for and RGB wheel on an xy graph
//		reverseList(), reverses a list, handy for prompt,popups
//		UniqueWaveValues() & UniqueWaveValuesT(), get a new wave with all duplicate values in the source wave removed
//		monotonic(a), checks if a wave is monotonic
//		isWaveConstant(), returns 1 if all elements of wav are the same
//		isdigit(c) & isletter(c), handy utilities
//		countOccurancesOfStr(main,find), count number of times find occurs in main
//		angleVec2Vec(a,b), finds angle between two vectors (degree)
//		perp2Vector(a), returns a free vector that is perpendicular to a, the direction around a is unspecified.
//		rotationAngleOfMat(rot)  finds the total rotation angle of a matrix 'rot'
//		isRotationMat(mat,[tol]), returns true if mat is a rotation matrix
//		axisOfMatrix(mat,axis,[squareUp]), find axis and angle from the rotation matrix mat
//		SquareUpMatrix(rot), turns rot from almost a rotation matrix to a true rotation matrix
//		rotationMatFromAxis() returns the rotation matrix about axis
//		smallestNonZeroValue(vec,[tol]), returns the abs( smallest non-zero element ), e.g. {0, -0.1, 3} returns 0.1
//		MedianOfWave(), returns median (or other percentile) of a wave, useful for picking scaling ranges
//		rangeOfVec(wav,[row,col,layer,chunk]), returns the max and min number in the specified vector
//		FindScalingFromVec(), find the scaling (for a SetScale command) that fit the values of vec[] (position of a step scan)
//		FindStepSizeInVec(),  find the step size in a vec (positions of a step scan) by examining the values
//		RangeOfValuesContainingFraction(),  returns range of VALUES of wave that constitute a given fraction of all values in wave
//		fractionalLevelFind(wwIN,frac), find the val such that frac of the values in wwIN are < val.  Faster than sorting for large waves
//		roundSignificant(val,N), returns val rounded to N places
//		placesOfPrecision(a), returns number of places of precision in a
//		ValErrStr(val,err), returns string  "val � err" formatted correctly
//		num2strFull(n,[fmt,tol]), like num2str, but for integers, gives full precision (needed in the range processing)
//		normalize(a), normalizes a if it is a vector or square matrix
//		isPositiveInts(ww), returns 1 only if ww is positive ints, useful for determining if ww are counts
//		arithmetic(expression), return value of expression as a real number, e.g. "1/3 + 5" --> 5.33333   or "sqrt(2)" --> 1.41421
//		maxValueOfType(w), returns largest number that can be stored in a wave of specified type
//		DefaultZeroThresh(w), returns a smallest number that the wave type can store, useful for printing zeros.
//		OrderValues(x0,x1), x0 and x1 are optionally swapped so that x0<x1 (x0 & x1 are passed by reference)
//		GetCursorRangeFromGraph(), get the cursor range from a graph (returns key=value string)
//		FWHM of fitted peaks: GaussianFWHM(W_coef), LorentzianFWHM(W_coef)
//		Area of fitted peaks: LorentzianIntegral(W_coef), GaussianIntegral(W_coef), Gaussian2DIntegral(W_coef)
//		computeCOM(), compute the Center of Mass
//		PowerIntegerScale(), rescale a waves values by ^n or ^(1/n), preserves sign for non-complex values
//		FitErrorString(FitError,FitQuitReason), return a string representation of the fitting error
//		Posix2HFS, a replacement for PosixToHFS(), (using ParseFilePath() for HFSToPosix()) we no longer need HFSAndPosix.xop
//		pingHost(host), returns ping time in seconds, returns NaN CANNOT ping the host
//		cpuFrequency(), systemUserName(), sytemHostname(), localTimeZoneName(), getEnvironment(), returns system info
//		FindFirstCharInStr(str,chars), finds first occurance of any of the chars in chars that occur in str
//		RemoveLeadingString(str,head,ignoreCase), removes head from start of str
//		RemoveTrailingString(str,tail,ignoreCase), removes tail from end of str
//		TrimBoth(str,[chars,ignoreCase]), TrimFront(), & TrimEnd(),  trim white space or given set of characters
//		use functions in line above:  TrimFrontBackWhiteSpace(str), TrimLeadingWhiteSpace(str), TrimTrailingWhiteSpace(str), trims whitespace
//		countChars(buf, chars), count number of times one of the characters in chars occur in buf
//		PrintLongStrings(buf,[sep]), prints really long strings, returns number of lines printed
//		IgorFileTypeString() gives descriptive string from the NUMTYPE from WaveInfo()
//		GenericWaveNoteInfo(), returns wave note info
//		StopAllTimers(), stops all the Igor timers
//		dateStr2ISO8601Str(), convert a date to an ISO 8601 format
//		ISOtime2niceStr(iso), convert ISO8601 string to a nice format for graph annotations
//		ISOtime2IgorEpoch(iso), convert ISO8601 string to an Igor Epoch (error returns NaN), if TZ present, then returns epoch in UTC
//		epoch2ISOtime(seconds), convert an Igor epoch (in seconds) to an ISO8601 format string
//		AskForUserForDateTime(epoch), puts up a dialog to select a date & time, returns the epoch
//		ElapsedTime2Str(seconds,[showSec,fracDigits]), convert seconds to a nice time interval string
//		WeekDayNumber(month, idate, year), returns week day number [1,7]
//		vec2str(), convert a vector to a string
//		str2vec(), convert a string to a free vector
//		encodeMatAsStr(mat,[places]), convert mat to a string interpretable by decodeMatFromStr(), used for wave notes
//		decodeMatFromStr(str), returns a FREE wave defined by str, inverse of encodeMatAsStr()
//		cmplx2str(zz,[places,mag]), convert complex number to a printable string
//		str2cmplx(str), this is like str2num, but for complex
//		num2fraction(val,tol,[addPlus]), convert val to closest fraction string that is within tol
//		vec2MINstr(vecIN), similar to hkl2str(), convert vecIN to a string of acceptable minimal length
//		minStr2Vec(inStr,Nreq), similar to str2hkl(), convert a string of numbers to a wave of Nreq values, pretty forgiving about format
//		ReplaceCharacters(chars,inStr,replacement)  replace all occurance of a character in chars[] with replacement
//		nice printing of vecs & mats: printWave() and associated functions: {printvec(),printmat(),printmatOneListReal(),printmatOneListComplex(),printmatOneListText()}
//		SetAspectToSquarePixels(), Used to square up a graph window
//		GetGraphAspectRatio(gName), returns aspect ratio of graph in units
//		SquareUpGizmo(gName), Used to square up a graph gizmo
//		SetGizmoZoom(zoom), Interactive Set zoomFactor for a Gizmo (only for Igor 7)
//		Name2SymbolCharacter(name), return character equivalent in Symbol font, "theta" returns "q"
//		Letter2SymbolOrUnicode(letter), for "theta", returns either unicode or "\F'Symbol'q\F]0"
//		added LetterName2Unicode(letter), returns unicode version of letter {alpha,Beta,...} and {Alef, bet, gimel...}
//			for Igor6, returns Mac Roman keyboard characters as well as it can.
//		num2Ordinal(n), converts integer n to ordinal string, e.g. 1 --> "1st"
//		ChangeStrEnding(oldEnd, inStr, newEnd)  if inStr ends in oldEnd, replace oldEnd with newEnd
//		SIprefix2factor(prefix)  get the factor for an SI prefix
//		ConvertUnits2meters(unit,[defaultLen])  returns conversion factor from length unit to meters
//		ConvertTemperatureUnits(Tin,unitIn,unitOut)  returns Temperature in units of (unitOut)
//		ConvertNumWithUnitsLength(in,outUnit,[rnd,space])  can convert "2in" to "0.0508 m"
//		SplitNumAndUnits(in)  take string such as "-9.3e2 mi" --> "-9.3e2;mi",  or "2 nm^-1" --> "2;nm^-1"
//		ConvertUnits2seconds(unit,[defaultLen])  returns conversion factor from time unit to seconds
//		ConvertUnits2kg(unit,[defaultMass])  returns conversion factor from unit to killo-grams
//		ConvertUnits2Joules(unit,[defaultEnergy])  returns conversion factor from unit to Joules
//		ConvertUnits2Pascal(unit,[defaultEnergy])  returns conversion factor from unit to Pascal
//		ConvertAngleUnits(angle,unitIn,unitOut, [defaultU])  returns angle(unitIn) in units of (unitOut)
//		ConvertUnits2Degree(unit,[defaultFactor, isCos])  returns conversion factor from unit to degree for angles
//		RomanNumeral(j) converts a number to a Roman Numeral string, NO upper limit, so watch out for string length
//		RomanNumeral2Int(str) convert a Roman Numeral string to an integer
//	9	Old legacy or deprecated functions



//  ======================================================================================  //
//  =========================== Start of Option Menu Functions ===========================  //

// generic, used for lots of Menus
Function/S MenuItemIfWaveClassExists(item,classes,optionsStr,[invisible,all])
	String item						// string that you want to appear in menu
	String classes					// semi-colon separated list of wave classes, can use "*" in classes
	String optionsStr				// options that are found in WaveList function optionsStr
	Variable invisible			// controls menu item when conditions not met: true -> menu item absent, false or absent -> menu item grayed out
	Variable all					// when all is TRUE, then all of the classes in waveClassList must be present, not just one
	invisible = ParamIsDefault(invisible) ? 0 : invisible
	invisible = numtype(invisible) ? 0 : invisible
	all = ParamIsDefault(all) ? 0 : !(!all)
	all = numtype(all) ? 0 : all
	String list = WaveListClass(classes,"*",optionsStr,all=all)
	if (invisible)
		return SelectString(strlen(list),"",item)
	endif
	return SelectString(strlen(list),"(","")+item
End

// menu active only if top graph is an image
Function/S MenuItemIfTopGraphImage(item)
	String item
	if (strlen(ImageNameList("",""))<1)
		return "("+item					// top graph does not contain an image, so disable menu item
	endif
	return item
End

// enable menu item when a particular wave class is displayed (on a graph, table, or gizmo)
Function/S MenuItemIfWaveClassDisplayed(item,classes,optionsStr,flag)
	String item					// item name to return, or "(item" to disable
	String classes				// list of classes to check for (semi-colon separated)
	String optionsStr			// options that are found in WaveList function optionsStr
	Variable flag				// use: 1=graphs, 2=tables, 4=gizmos, or any sum of 1,2,4.  -1 or 7 gives all
	String wList=WaveListClass(classes,"*",optionsStr)	// a list of waves in requested classes
	Variable i, N=ItemsInList(wList)
	for (i=0;i<N;i+=1)
		WAVE wav=$StringFromList(i,wList)			// for each ODF wave
		if (strlen(WindowsWithWave(wav,flag))>0)// found a window displaying wav
			return item
		endif
	endfor
	return "("+item										// found none
End

// menu active only if a graph is displayed, useful when using a Cursor on the graph
Function/S MenuItemIfAnyGraphExists(item,[visible])
	String item
	Variable visible						// default is to ONLY show visible graphs, use visible=0 to include both
	visible = ParamIsDefault(visible) || numtype(visible) ? 1 : visible
	String optStr = SelectString(visible, "WIN:1", "WIN:1,VISIBLE:1")
	if (strlen(WinList("*",";",optStr))<1)
		return "("+item					// top graph does not contain an image, so disable menu item
	endif
	return item
End

// for menus, enable menu item when a particular wave class is present on the graphName
Function/S MenuItemsWaveClassOnGraph(item,classes,graphName,[csr])
	String item					// item name to return, or "(item" to disable
	String classes					// list of classes to check for (semi-colon separated)
	String graphName				// name of graph, use "" for top graph
	String csr						// requires cursor csr, can be one of "A", "B", ..., "J"
	csr = SelectString(ParamIsDefault(csr),csr,"")	// default to none

	if (strlen(ImageNameList(graphName,";")+TraceNameList(graphName,";",3))<1)
		return "("+item									// no images or graphs displayed
	elseif (strlen(csr) && strlen(CsrInfo($csr,graphName))<1)
		return "("+item									// cursor csr is NOT on graph
	endif

	Variable i
	String list = ImageNameList(graphName,";")	// first check the images
	for (i=0;i<ItemsInList(list);i+=1)
		Wave ww = ImageNameToWaveRef(graphName, StringFromList(i,list))
		if (WaveInClass(ww,classes))
			return item
		endif
	endfor

	list = TraceNameList(graphName,";",3)		// second check ordinary traces and contours
	for (i=0;i<ItemsInList(list);i+=1)
		Wave ww = TraceNameToWaveRef(graphName, StringFromList(i,list))
		if (WaveInClass(ww,classes))
			return item
		endif
	endfor
	return "("+item						// top graph does not contain a wave of desired class
End


// This is for when we want to check if any folder under the current folder contains wave with a certain class
Function/S MenuItemFolderWithClassExists(item,classes,optionsStr,[invisible,all])
	String item						// string that you want to appear in menu
	String classes					// semi-colon separated list of wave classes, can use "*" in classes
	String optionsStr				// options that are found in WaveList function optionsStr
	Variable invisible			// controls menu item when conditions not met: true -> menu item absent, false or absent -> menu item grayed out
	Variable all					// when all is TRUE, then all of the classes in waveClassList must be present, not just one
	invisible = ParamIsDefault(invisible) ? 0 : invisible
	invisible = numtype(invisible) ? 0 : invisible
	all = ParamIsDefault(all) ? 0 : !(!all)
	all = numtype(all) ? 0 : all

	String list = WaveListClass(classes,"*",optionsStr,all=all)	// check current folder
	if (strlen(list)<1)
		list = FoldersWithWaveClass("",classes,"*",optionsStr,all=all)
	endif
	if (invisible)
		return SelectString(strlen(list),"",item)
	endif
	return SelectString(strlen(list),"(","")+item
End


// This is really useful with an  Execute/P ... command
// e.g.  MenuItemIfWindowAbsent("Include ABC Support","ABC.ipf","WIN:128"), Execute/P "INSERTINCLUDE  \"ABC\"";Execute/P "COMPILEPROCEDURES "
Function/S MenuItemIfWindowAbsent(item,win,options)		// Shows menu item if win NOT present
	String item			// the string that appears in the menu
	String win			// may contains wildcard
	String options		// options that are in WinList()
	Variable present = (strlen(WinList(win,"",options)))>0
	return SelectString(present,"","(")+item
End
//
Function/S MenuItemIfWindowHidden(item,win)		// Shows menu item if win NOT present
	String item
	String win
	GetWindow/Z $win, hide
	return SelectString(V_flag,"(","")+item
End

// only show menu item when a single specific wave exists
Function/S MenuItemIfWaveExists(item,wname)
	String item
	String wname

	Variable there = Exists(wname)==1
	return SelectString(there,"(","")+item
End

// useful when just looking for a type of wave, no class information
Function/S MenuItemIfWavesExists(item,matchStr,optionsStr,[invisible])
	String item						// string that you want to appear in menu
	String matchStr
	String optionsStr				// options that are found in WaveList function optionsStr
	Variable invisible				// controls menu item when conditions not met: true -> menu item absent, false or absent -> menu item grayed out
	invisible = ParamIsDefault(invisible) ? 0 : invisible
	invisible = numtype(invisible) ? 0 : invisible
	String list = WaveList(matchStr,";",optionsStr)
	if (invisible)
		return SelectString(strlen(list),"",item)
	endif
	return SelectString(strlen(list),"(","")+item
End

// useful when working with a specific window, particularly with a Panel
Function/S MenuItemIfWinExists(item,win,optionsStr)
	String item								// menu item text
	String win								// win name (can contain a wild card)
	String optionsStr						// optionsStr from WinList()
	if (strlen(WinList(win,"",optionsStr))<1)
		return "("+item					// win not found, so disable menu item
	endif
	return item
End

// only show menu item when a funcName exists with specified values in optStr
//	e.g.  MenuItemIfFunctionExists("Write Detector values to EPICS...","EPICS_put_PV_num")
Function/S MenuItemIfFunctionExists(item,funcName,[optStr])
	String item
	String funcName
	String optStr					// key-values must match result of FunctionInfo() "key:value;"
	optStr = SelectString(ParamIsDefault(optStr),optStr,"")
	// optStr example is:  "PROCWIN:Procedure;RETURNTYPE:3;"

	String keyVals = FunctionInfo(funcName)
	if (strlen(keyVals)<1)		// no function info, add "("
		return "("+item
	endif

	if (strlen(optStr)>1)		// an optStr was passed, check for matching
		String key, value
		Variable i, N=ItemsInList(optStr)
		for (i=0;i<N;i+=1)
			key = StringFromList(0, StringFromList(i,optStr),":")
			value = StringByKey(key,optStr)
			if (!StringMatch(StringByKey(key,keyVals), value))
				return "("+item
			endif
		endfor
	endif

	return item						// either no optStr, or all key:values in optStr match
End

// only show menu item when the scrap (clipboard) contains valid text as decided by validTextFuncName()
Function/S MenuItemIfScrapValid(item, validTextFuncName)
	String item						// string that you want to appear in menu
	String validTextFuncName		// name of function used to test the scrap
	FUNCREF EmptyIsInvalidTextProto func = $validTextFuncName
	Variable show = func(GetScrapText())	// test the scrap
	return SelectString(show,"(","")+item
End
//
Function EmptyIsInvalidTextProto(text)		// returns true if text is NOT empty
	String text
	return strlen(text)>0
End

//  ============================ End of Option Menu Functions ============================  //
//  ======================================================================================  //


//  ======================================================================================  //
//  ============================== Start of Corner Labels ================================  //

// Puts text on the lower right corner showing when the layout was made, and on the lower left showing where it came from
Proc Layout_Corner_Labels_Style_()		// This is needed to make it show up in the "Layout Macros" menu
	Add_Corner_Labels_To_Layout()
End
Proc Add_Corner_Labels_To_Layout() : LayoutStyle
	PauseUpdate; Silent 1		// modifying window...
	Textbox/C/N=stamp0/F=0/A=RB/X=0.1/Y=0.1 "\\Z06\\{\"%s %s\",date(), time()}"
	Textbox/C/N=stamp1/F=0/A=LB/X=0.1/Y=0.1 "\\Z06\\{\"%s\",JZTutil#CornerStampWindow()}"
EndMacro
//
Static Function/S CornerStampWindow()
	String win=WinName(0,127+GIZMO_WIN_BIT,1)// Graph, Table, Layout, Notebook,Panel, or XOP window
	String gwin=""

	switch(WinType(win))
		case 3:						// a layout, find the first graph on this layout
			String list="x"
			Variable i
			for (i=0; strlen(gwin)==0 && strlen(list); i+=1)		// get first graph name
				list = LayoutInfo(win,num2istr(i))
				if (stringmatch(StringByKey("TYPE",list),"Graph"))
					gwin = StringByKey("NAME",list)
				endif
			endfor
			break
		case 1:						// a graph, gwin is win
		case 2:						// a table, gwin is win
		case GIZMO_WIN_TYPE:	// from an XOP, probably a gizmo
			gwin = win
			break
		default:						// only layout or graph
			gwin = ""
	endswitch
	gwin = SelectString(strlen(gwin),"",":")+gwin

	PathInfo home					// creates String s_path
	return S_path+IgorInfo(1)+gwin
End


// Puts text on the lower right corner showing when the graph was made, and on the lower left showing where it came from
Proc Add_Corner_Labels_To_Graph() : GraphStyle
	AddCornerLabelsToGraph()
EndMacro
Function AddCornerLabelsToGraph([outside,size])
	Variable outside
	Variable size					// font size
	outside = ParamIsDefault(outside)||numtype(outside) ? 0 : !(!outside)
	Variable defSize = outside ? 4 : 6
	size = ParamIsDefault(size) ? defSize : size
	size = size==limit(size,4,48) ? size : defSize
	if (strlen(WinList("*",";","WIN:1"))<1)
		return 1
	endif
	Variable bottom=0
	if (outside)
		GetWindow kwTopWin , psize
		Variable pBottom=V_bottom
		Variable pheight = V_bottom-V_top
		GetWindow kwTopWin , gsize
		Variable gBottom=V_bottom
		bottom = -floor(100*(gBottom-pBottom)/pheight)
	endif
	String str
	sprintf str,"%02d",size
	Textbox/C/N=stamp0/F=0/A=RB/X=0.2/Y=(bottom)/E=2 "\\Z"+str+"\\{\"%s %s\",date(), time()}"
	Textbox/C/N=stamp1/F=0/A=LB/X=0.2/Y=(bottom)/E=2 "\\Z"+str+"\\{\"%s\",JZTutil#CornerStampWindow()}"
	return 0
End

//  ================================ End of Corner Labels ================================  //
//  ======================================================================================  //


//  ======================================================================================  //
//  ========================== Start of Multiple Graphs on Page ==========================  //

// Put up multigraph layouts.  Each time you call this, it adds another graph to the layout
Window AppendMultiGraph2LayoutN() : Layout			// this is needed to make it appear in the "Layout Macros" menu
	AppendGraph2LayoutN(NaN,"","")
End
// Append a graph to a multiple graph layout.  This is for making multi-graph layouts
Function/S AppendGraph2LayoutN(Nmax,orientation,gName)
	Variable Nmax									// requested maximum number of graphs in the layout, one of Nlist
	String orientation								// this must be "portrait" or "landscape"
	String gName									// name of graph to append, if empty string, DoPrompt

	if (!stringmatch("portrait",orientation) && !stringmatch("landscape",orientation))
		orientation = ""
	endif

	String Nlist = "1;2;3;4;6;8;9;12;16;20"		// list of allowed values for Nmax
	Make/N=(ItemsInList(Nlist))/FREE Nvals=str2num(StringFromList(p,Nlist))
	Variable m=BinarySearchInterp(Nvals,Nmax)
	Nmax = (m>=0) ? Nvals[ceil(m)] : -1			// Nmax is now one of Nlist, or Nmax is -1

	Variable printIt=0
	String lname
	if (strlen(gName)<1 || Nmax<0)
		if (Nmax<0)
			lname = StringFromList(0,WinList("Layout*",";","WIN:4"))
			Nmax = NumberByKey("Nmax",GetUserData(lname,"","AppendGraph2Layout"),"=")
			Nmax = numtype(Nmax) ? -1 : Nmax
		endif
		Nmax = Nmax<0 ? 12 : Nmax				// 12 is default for absurd inputs
		Nmax = WhichListItem(num2istr(Nmax),Nlist)+1
		Prompt gName, "Graph to add to layout",popup,NOTonLayout(lname,WinList("*",";","WIN:1"))
		Prompt Nmax,"number of graphs/layout",popup,Nlist
		DoPrompt "Pick a Graph",gName,Nmax
		if (V_flag)
			return ""
		endif
		Nmax = str2num(StringFromList(Nmax-1,Nlist))
		printIt = 1
	endif
	if (strlen(gName)<1)
		DoAlert 0, "no graph to add"
		return ""
	endif

	lname = StringFromList(0,WinList("Layout"+num2istr(Nmax)+"_*",";","WIN:4"))		// get layout to append to
	if (strlen(LayoutInfo(lname,gName))>1)
		return lname								// graph is on layout, all done
	endif
	// check if this layout is full
	Variable i,N,NL = NumberByKey("NUMOBJECTS", LayoutInfo(lname,"Layout"))
	for (i=0,N=0;i<NL;i+=1)					// check each object in the layout
		N += stringmatch(StringByKey("TYPE",LayoutInfo("", num2istr(i))),"Graph")	// increlement if obj is graph
	endfor
	if (N>=Nmax)									// this layout is full, force creation of a new one
		lname = ""
	endif
	if (strlen(lname)<1)
		if (strlen(orientation)<1)					// orientation was not passed, ask now
			orientation = SelectString(Nmax==8,"Landscape","Portrait")
			Prompt orientation, "page orientation",popup,"Portrait;Landscape"
			DoPrompt "Pick a Graph",orientation
			if (V_flag)
				return ""
			endif
			printIt = 1
		endif
		lname = UniqueName("Layout"+num2istr(Nmax)+"_",8,0)		// need to make a newlayout
		NewLayout/P=$orientation
		ModifyLayout mag=0.5, units=0, frame=0
		DoWindow/C $lname
		SetWindow kwTopWin ,userdata(AppendGraph2Layout)=ReplaceNumberByKey("Nmax","",Nmax,"=")
	endif

	String page = StringByKey("PAGE", LayoutInfo(lname,"Layout"))	// get page size to determing landscape or portrait
	Variable aspect = (str2num(StringFromList(3,page,","))-str2num(StringFromList(1,page,",")))	// aspect = height/width
	aspect /= (str2num(StringFromList(2,page,","))-str2num(StringFromList(0,page,",")))
	orientation = SelectString(aspect>1,"landscape","portrait")
	String rc = StringByKey(num2istr(Nmax),"1:1,1;2:2,1;3:3,1;4:2,2;6:3,2;8:4,2;9:3,3;12:4,3;16:4,4;20:4,5")
	Variable rows=str2num(StringFromList(0,rc,","))
	Variable columns=str2num(StringFromList(1,rc,","))
	if (stringmatch("landscape",orientation))		// swap rows and columns for Landscape orientation
		Variable swap = rows
		rows = columns
		columns = swap
	endif
	AppendLayoutObject/F=0/W=$lname graph $gName
	String cmd
	sprintf cmd, "Tile /A=(%d,%d)/O=1",rows,columns
	Execute cmd
	RemoveLayoutObjects/W=$lname/Z stamp0,stamp1	// make sure that text boxes are written last, so they show
	Textbox/C/N=stamp0/W=$lname/F=0/A=RB/X=0.1/Y=0.1 "\\Z06\\{\"%s %s\",date(), time()}"
	Textbox/C/N=stamp1/W=$lname/F=0/A=LB/X=0.1/Y=0.1 "\\Z06\\{\"%s\", JZTutil#CornerStampWindow()}"

	if (printIt && (strlen(GetRTStackInfo(2))<1))
		printf "\tAppendGraph2LayoutN(%d,\"%s\",\"%s\")\r",Nmax,orientation,gName
	endif
	return lname
End
//
Static Function/T NOTonLayout(lname,listIN)
	String lname				// name of Layout
	String listIN				// list of things to check, remove items that are on layout

	if (strlen(lname)<1)	// no Layout, none of list are on it
		return listIN
	endif

	String item, listOUT=""
	Variable i, N=ItemsInList(listIN)
	for (i=0;i<N;i+=1)
		item = StringFromList(i,listIN)
		if (strlen(LayoutInfo(lname,item))<1)
			listOUT += item+";"
		endif
	endfor
	return listOUT
End
//
//Function xxx(Nmax)		// a test routine for AppendGraph2LayoutN()
//	Variable Nmax
//	Make/N=100/O y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12
//	SetScale/I x 0,10,"", y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12
//	y0 = 0.2*x -1
//	y1 = sin(x/2) ; y2 = sin(x) ; y3 = sin(2*x) ; y4 = sin(3*x) ; y5 = sin(4*x)
//	y5 = cos(x/2) ; y6 = cos(x) ; y7 = cos(2*x) ; y8 = cos(3*x) ; y9 = cos(4*x)
//	y10 = exp(x/10)/3  ;  y11 = exp(-x/4)  ;  y12 = 2*exp(-x/2)-1
//	Variable i
//	String win
//	for (i=0;i<=12;i+=1)
//		win = "Graph"+num2istr(i)
//		if (strlen(WinList(win,"","WIN:1"))<1)
//			Display/K=1 $("y"+num2istr(i))
//			DoWindow/C $win
//			ModifyGraph tick=2, minor=1, standoff=0, mirror=1
//			ModifyGraph axOffset(left)=-4.7,axOffset(bottom)=-1.5
//		endif
//		AppendGraph2LayoutN(Nmax,"",win)
//	endfor
//End

//  =========================== End of Multiple Graphs on Page ===========================  //
//  ======================================================================================  //


//  ======================================================================================  //
//  =============================== Start of String Ranges ===============================  //

// This section is for dealing with random or usually non contiguous sequence of integers
// i.e.  you took data in scans 1-30, but scans 17, and 25 were no good.  So the valid range is "1-16,18-24,26-30"
//  or perhaps you want to combine scans "3,7,9"  The following routines handle those situations in a simple fashion.

ThreadSafe Function isValidRange(range)	// returns 1 if range is a valid string range
	String range								// note, "1-Inf" is a valid sub range, but "-Inf-6" is not

	range = ReplaceString(" ",range,"")
	String subRange
	Variable i, m, N=ItemsInList(range,","), first,last, step,prev=-Inf
	for (i=0;i<N;i+=1)						// check each item for validity
		subRange = StringFromList(i,range,",")
		first = str2num(subRange)
		if (numtype(first))					// this sub range is not even numeric or is Inf
			return 0
		elseif (first<=prev)				// start of this sub range is not after previous one
			return 0
		endif

		m = strsearch(subRange,"-",1)	// occurance of first "-", after first character
		if (m<0)									// a single number, a simple sub range
			prev = first
			continue
		endif
		last = str2num(subRange[m+1,Inf])
		step = str2num(StringFromList(1,subRange,":"))
		step = numtype(step)==2 ? 1 : step	// no step found, default to 1

		if (!(last>=first) || step<=0)	// must have first<=last and positive non-zero step
			return 0
		endif
		prev = floor((last-first)/step)*step+first	// last number in this sub range
	endfor
	return 1
End

ThreadSafe Function NextInRange(range,last)	// given a string like "2-5,7,9-12,50" get the next number in this compound range
										// the range is assumed to be monotonic, it returns NaN if no more values
	String range					// list defining the range
	Variable last					// last number obtained from this range, use -Inf to get start of range, it returns the next

	// find first item in the list that should use next
	String item,item0,item1
	Variable m,i,j,step,ret
	Variable first						// first value in an item
	do
		item = StringFromList(j,range,",")
		first = str2num(item)
		// do we need to check for NaN in first?

		if (numtype(first))
			return NaN
		elseif (last<first)				// skipping to the next item
			return first
		endif

		// if last>=first, check to see if item is a '2-5' type range
		m=-1							// remove any leading white space from item
		do
			m += 1
		while (char2num(item[m])<=32)
		item = item[m,strlen(item)-1]

		//check if a step size is specified
		item1=StringFromList(1,item,":")
		item0 = StringFromList(0,item,":") // now item0 may be a "2-5" type range or a single number
		step = str2num(item1)
		step = numtype(step)==2 ? 1 : step
		if (!(step>0))								// step must be > 0
			return NaN
		endif

		i = strsearch(item0,"-",1)		// location of first '-' after the first character
		if (i<0)							// only a single number, not a dash type range, keep looking
			j += 1
			continue
		endif
		// check to see if last was in the range of item, but not the last value
		if (last>=str2num(item0) && last<str2num(item0[i+1,Inf]))
			if(step ==1)
				return last+1
			else	// a little more math if step > 1
				if(mod(last - first, step) ==0) // an exact member of the sequence
					ret = last + step
				else
					ret = first + ceil((last - first)/step)*step
				endif
				if(ret <= str2num(item0[i+1,Inf]) ) // return only if the result is no greater than the end of the range, e.g. in the case of NextInRange("2-7:2",6)
					return ret
				endif
			endif
		endif
		j += 1
	while(strlen(item)>0)
	return NaN
End
//Function TestNextInRange(range)
//	String range
//	range = SelectString(strlen(range) ,"-16,-7--3,-1-2,5,50-54,99",range)
//	Variable i = -Inf
//	printf "for range = {%s},   ", range
//	do
//		i = NextInRange(range,i)
//		if (numtype(i))
//			break
//		endif
//		printf "%d  ", i
//	while (!numtype(i))
//	print ""
//End

ThreadSafe Function PointInRange(i,range)	// return the ith number in range, first number is i==0
	Variable i								// index to number in range
	String range							// a string range

	Variable N=ItemsInRange(range)
	i = round(i)
	if (numtype(i) || i<0 || i>=N)		// check for valid index
		return NaN
	endif
	Variable j, m=-Inf
	for (j=0;j<=i;j+=1)					// this is not the fastest, but it is simple
		m = NextInRange(range,m)
	endfor
	return m
End

ThreadSafe Function ItemsInRange(range)	// given a string like "2-5,7,9-12,50" get the total number of values in the range
	String range				// list defining the range
								// the range is assumed to be monotonic, it returns NaN on error
	String item							// each of the comma sepated items
	String item0, item1
	Variable len=0						// the result, number of values represented
	Variable m,i,j,step,N=ItemsInList(range,",")

	for (j=0;j<N;j+=1)				// loop over each item
		item = StringFromList(j,range,",")

		m=-1							// remove any leading white space from item
		do
			m += 1
		while (char2num(item[m])<=32)
		item = item[m,strlen(item)-1]

		item1 = StringFromList(1,item,":")	// stepsize, if specified with ":"
		item0 = StringFromList(0,item,":")	// start-end of range is specified before ":"
		step = str2num(item1)
		step = numtype(step)==2 ? 1 : step
		if (!(step>0))								// step must be > 0
			return NaN
		endif

		i = strsearch(item0,"-",1)		// location of first '-' after the first character
		if (i<0)							// only a single number, not a dash type range, keep looking
			len += 1
		else							// item0 is a dash type 
			len += ceil((str2num(item0[i+1,Inf])-str2num(item0)+1)/step)
		endif
	endfor
	return len
End
//	print ItemsInRange("-16,-7--3,-1-2,5,50-54,99")
//  17
//	print ItemsInRange("1-5")
//  5


//Function firstInRange(range)
//	String range
//	return str2num(range)
//End
//

ThreadSafe Function lastInRange(range)		// returns the last number in the range, lastInRange("3,5,9-20") returns 20
	String range
	Variable i,last,first,step
	String item0, item1

	i = strsearch(range,",",Inf,1)
	if (i>=0)									// remove previous comma separated items
		range = range[i+1,Inf]
	endif
	range = ReplaceString(" ",range,"")		// spaces do not count
	range = ReplaceString("\t",range,"")	// tabs do not count

	item1 = StringFromList(1,range,":")	// stepsize, if specified with ":"
	item0 = StringFromList(0,range,":")	// start-end of range is specified before ":"
	step = str2num(item1)
	step = numtype(step)==2 ? 1 : step
	if (!(step>0))									// step must be > 0
		return NaN
	endif

	i = strsearch(item0,"-",Inf,1)
	if (char2num(item0[i-1])==45)		// a double dash, we have a minus sign
		i -= 1
	endif
	i += i>0 ? 1 : 0							// increment i to skip over dash unless there is no dash in string, i<0 means no dash present
	last = str2num(item0[i,Inf])
	first = str2num(item0)
	return floor((last-first)/step)*step+first
End

ThreadSafe Function isInRange(range,m)// returns TRUE if m is a number in range
	String range							// list defining the range
	Variable m								// number of interest, is m part of range?

	String item							// each of the comma sepated items
	Variable i,j,N=ItemsInList(range,",")
	String item0, item1
	Variable step
	
	for (j=0;j<N;j+=1)					// loop over each comma separated item
		item = StringFromList(j,range,",")
		item = TrimLeadingWhiteSpace(item)
		i = strsearch(item,"-",1)		// location of first '-' after the first character
		if (i<0)							// only a single number, not a dash type range, keep looking
			if (m==str2num(item))
				return 1					// this number is m
			endif
		else								// item is a dash type
			item1 = StringFromList(1,item,":")	// stepsize, if specified with ":"
			item0 = StringFromList(0,item,":")	// everything else before ":"
			step = str2num(item1)
			step = numtype(step)==2 ? 1 : step
			if (!(step>0))								// step must be > 0
				return 0
			endif

			if (m>=str2num(item0) && m<=str2num(item0[i+1,Inf]))  //  m is in this continuous range
				return (mod((m - str2num(item0)), step)== 0)  // true only if (m - first) can be divided by step
			endif
		endif
	endfor
	return 0
End

Function/T subRangeOfRange(fullRange,jlo,jhi)
	// returns new string range containing only values in [jlo,jhi]
	String fullRange							// original string range
	Variable jlo,jhi							// want new range to go from jlo ... jhi
	jlo = round(jlo)
	jhi = round(jhi)

	String item, subRange=""
	Variable i, delta, ifirst, ilast, stride
	for (i=0; 1; i+=1)
		item = StringFromList(i,fullRange,",")
		ifirst = str2num(item)
		ilast = lastInRange(item)
		if (numtype(ifirst+ilast))		// finished all the simple ranges
			break
		elseif (ilast < jlo)				// starting point is further on, try next item
			continue
		elseif (ifirst > jhi)				// any further simple ranges are past jhi, done
			break
		endif

		stride = str2num(StringFromList(1,item,":"))
		stride = stride > 0 ? stride : 1
		if (strlen(subRange)==0)
			delta = max(ceil((jlo-ifirst)/stride),0)
			jlo = ifirst + stride*delta	// first point equal or greater than jlo in item
			subRange = num2istr(jlo)		// start to build the first piece of subRange
			if (jhi >= ilast && ilast>ifirst)	// add all of this simple range
				subRange += "-"+num2istr(ilast)
				subRange += SelectString(stride>1, "", ":"+num2istr(stride) )	// add the :stride
			elseif (ilast>ifirst)			// can only use part of this simple range
				jhi = ifirst + stride*floor((jhi-ifirst)/stride)	// last point less than or equal to jhi in item
				subRange += "-"+num2istr(jhi)
				subRange += SelectString(stride>1, "", ":"+num2istr(stride) )	// add the :stride
			endif
		elseif (ilast <= jhi)				// jhi is past this simple range, just append item to subRange
			subRange += ","+item
		else										// ilast (from item) is in this item (this is the last simple range)
			jhi = ifirst + stride*floor((jhi-ifirst)/stride)	// last point less than or equal to jhi in item
			subRange += ","+num2istr(ifirst)	// start to build the first piece of subRange
			if (jhi > ifirst)					// something to add to subRange
				subRange += "-"+num2istr(jhi)
				subRange += SelectString(stride>1, "", ":"+num2istr(stride) )	// add the :stride
			endif
		endif
	endfor

	return subRange
End
//	Function test()
//		Variable bad
//		bad = bad || test1("1-100,105-140",10,20, "10-20")
//		bad = bad || test1("1-100,105-140",99,107, "99-100,105-107")
//		bad = bad || test1("1-15",10,20, "10-15")
//		bad = bad || test1("-10-100,105-140",99,107, "99-100,105-107")
//		bad = bad || test1("-10-100,105-140",-5,107, "-5-100,105-107")
//		bad = bad || test1("-10-100,105-140",-50,107, "-10-100,105-107")
//		bad = bad || test1("10-100,105-140",10,20, "10-20")
//		bad = bad || test1("10-100,105-140",5,20, "10-20")
//		bad = bad || test1("-10-100,105-140",-50,500, "-10-100,105-140")
//		bad = bad || test1("1-100,105-140",-50,107, "1-100,105-107")
//		bad = bad || test1("1-100,105-140",5,107, "5-100,105-107")
//		bad = bad || test1("1-100,105-140",5,500, "5-100,105-140")
//		bad = bad || test1("1-100,105-140,490-510",5,500, "5-100,105-140,490-500")
//	
//		bad = bad || test1("0-100:2",10,20, "10-20:2")
//		bad = bad || test1("0-100:2",-10,20, "0-20:2")
//		bad = bad || test1("0-100:2",-10,100, "0-100:2")
//		bad = bad || test1("0-100:2",-10,101, "0-100:2")
//		bad = bad || test1("0-100:2,105-140",98,107, "98-100:2,105-107")
//		bad = bad || test1("1-100:2,105-140",95,107, "95-99:2,105-107")
//		print SelectString(bad,"all OK", "ERROR")
//	End
//	Static Function test1(fullRange,jlo,jhi, answer)
//		String fullRange
//		Variable jlo,jhi
//		String answer
//		String subRange = subRangeOfRange(fullRange,jlo,jhi)
//		Variable bad = !StringMatch(subRange,answer)
//		if (bad)
//			printf "%ssubRangeOfRange(\"%s\", %d, %d) --> \"%s\"\r",SelectString(bad, "    ", "*** "), fullRange,jlo,jhi,subRange
//		endif
//		return bad
//	End



// Caution expandRange("1-100000",";") will produce a very long string!  Use NextInRange() to avoid this problem
ThreadSafe Function/S expandRange(range,sep)	// expand a string like "2-5,7,9-12,50" to "2,3,4,5,7,9,10,11,12,50"
	String range					// series of numberseparated by commas and dashes, white is space igaored
	String sep						// separates final list, usually ";"
	if (strsearch(range,"Inf",0,2)>=0)	// cannot put an infinite number of characters into a string
		return ""
	endif
	if (strlen(sep)<1)				// sep defaults to ';'
		sep = ";"
	endif
	Variable i1,i2,i
	String str,out=""
	Variable N=ItemsInList(range,",")
	if (N<1)
		return ""
	endif
	Variable j=0, step
	String item1,item0
	do
		str = StringFromList(j, range, ",")
		Variable m=-1				// remove any leading white space
		do
			m += 1
		while (char2num(str[m])<=32)
		str = str[m,strlen(str)-1]

		// now check str to see if it is a range like "20-23" or "20-30:5"
		item1 = StringFromList(1,str,":")	// stepsize, if specified with ":"
		item0 = StringFromList(0,str,":")	// start-end of range is specified before ":"
		step = str2num(item1)
		step = numtype(step)==2 ? 1 : step
		if (!(step>0))								// step must be > 0
			return ""
		endif
		
		i1 = str2num(item0)                        // i1 is the start of the range
		i = strsearch(item0,"-",strlen(num2strFull(i1)))		// position of "-" after first number
		if (i>0)
			i2 = str2num(item0[i+1,inf])  // i2 is the specified end of the range
			i = i1
			do
				out += num2strFull(i)+sep
				i += step
			while (i<=i2)
		else
			out += num2strFull(i1)+sep
		endif
		j += 1
	while (j<N)

	i = strlen(out)-1
	if (char2num(out[i])==char2num(sep))  // remove the trailing sep character
		out = out[0,i-1]
	endif
	return out
End

// This is the inverse of expandRange()
ThreadSafe Function/S compressRange(range,sep) 	// take a range like "1;2;3;4;5;9;10;11" change it to "1-5,9-11"
	String range
	String sep							// sep is the separator used, will be replaced with commas and dashes

	String comp=""						// the compressed string
	String num
	Variable j,first,last,i=0,aa,bb,cc,step,in_subrange,write_subrange,writeend
	Variable N=ItemsInList(range,sep)
	if (N<1)
		return ""
	endif
	range = SortList(range,sep,2)				// make list monotonic
	
	if(N<3)  	// make no change if list has 1 or 2 elements
		comp = ","+ReplaceString(sep,range,",")			// change all sep to ","
		if (StringMatch(comp[strlen(comp)-1],",")>=0)	// remove a possible trailing ","
			comp = comp[0,strlen(comp)-2]
		endif
	else	// if list has more than 3 elements
		in_subrange = 0
		write_subrange = 0
		writeend = 0
		last = str2num(StringFromList(0,range,sep)) - 2
		for(i=2;i<N;i+=1)
			aa = str2num(StringFromList(i-2,range,sep))
			bb = str2num(StringFromList(i-1,range,sep))
			cc = str2num(StringFromList(i,range,sep))
			//if((cc-bb) == (bb-aa) && aa>last) // if the 3 consecutive numbers are equally spaced, and aa is not part of previous subrange
			if((cc-bb) == (bb-aa)) // if the 3 consecutive numbers are equally spaced
				if(!in_subrange && aa>last)  // no existing subrange, and aa is not part of previous subrange, then start a new one
					first = aa
					step = bb-aa
					last = cc
					in_subrange = 1
				elseif(!in_subrange && aa==last) // no existing subrange, but aa is the end of the previous subrange, then write down the subrange
					write_subrange = 1
				else  // already have an active subrange
					last = cc
				endif
				if(i==(N-1) && in_subrange) // if hit the last number while still in a subrange, then just write it down.
					write_subrange = 1
				elseif(i==(N-1)) // if hit the last number, and there's no active subrange, then write down the last two numbers
					writeend = 2
				endif
			else  // if they are not equally spaced
				if(in_subrange) // already started a subrange, then close it out
					in_subrange = 0
				elseif(aa==last) // no active subrange, but aa is the end of previous subrange, then write down the subrange
					write_subrange = 1
				else  // no active or recent subrange at all, then write down single item aa
					comp += ","+num2strFull(aa) 
				endif
				if(i == (N-1)) // if hitting the last number, then write down single elements cc or bb&cc
					if(bb == last)  // if bb is the end of the recently-closed subrange
						writeend = 1
						write_subrange = 1
					else    // if no active or recent subrange
						writeend = 2
					endif
				endif
			endif
			if(write_subrange)
				if(step ==1)
					comp += ","+num2strFull(first)+"-"+num2strFull(last)
				elseif (step==0)
					comp += ","+num2strFull(first)
				else
					comp += ","+num2strFull(first)+"-"+num2strFull(last)+":"+num2strFull(step)
				endif
				write_subrange = 0
			endif
			switch(writeend)
				case 1:
					comp += ","+num2strFull(cc) 
					break
				case 2:
					comp += ","+num2strFull(bb) + ","+num2strFull(cc) 
			endswitch
		endfor
	endif
	
	//remove the leading ","
	if(strlen(comp)>0)
		comp = comp[1,Inf]
	endif
	return comp
End
//Function test_compressRange()
//	String range
//	range = ""
//	printf "'%s' ---> '%s'\r",range, compressRange(range,";")
//	range = "1;2;3;4;5;9;10;11"
//	printf "'%s' ---> '%s'\r",range, compressRange(range,";")
//	range = "4"
//	printf "'%s' ---> '%s'\r",range, compressRange(range,";")
//	range = "4;7"
//	printf "'%s' ---> '%s'\r",range, compressRange(range,";")
//	range = "-10;-9;-8;-7;-6;-5;-3;-2;-1;0;1;2;7;9;22"
//	printf "'%s' ---> '%s'\r",range, compressRange(range,";")
//	range = "7;6;4;5;2;3"
//	printf "'%s' ---> '%s'\r",range, compressRange(range,";")
//	range = "-1;1;0;-3;2;-2"
//	printf "'%s' ---> '%s'\r",range, compressRange(range,";")
//	range = "1;1;1;10;10;10;20;21;22;"
//	printf "'%s' ---> '%s'\r",range, compressRange(range,";")
//End

////// not used ////////////// not yet compatible with "1-9:3" format /////////////
ThreadSafe Function/T compactRange(in) // take a range and return the most compact way of describing that range, e.g.  "1,2,3,4,6" -> "1-4,6"
	String in							// list defining the starting range, in is assumed to be monotonic
	String out=""						// list defining the minimalized output range

	Variable i,j0,j1					// i's come from in, j's refer to out
	i = str2num(in)
	if (numtype(i)==2)
		return ""
	endif

	String item						// one comma separated item
	Variable pitem					// index of item
	j0 = i
	j1 = i
	for (pitem=0,item=StringFromList(0,in,","); strlen(item); pitem+=1,item=StringFromList(pitem,in,","))	// for each comma separated item
		i = str2num(item)
		if ((i-j1)>1)					// this item is more then 1 from the last j, add to out
			if (j0==j1)				// a single number range
				out += num2istr(j0)+","
			else						// more then one in a dash range
				out += num2istr(j0)+"-"+num2istr(j1)+","
			endif
			j0 = i
		endif
		j1 = lastInRange(item)
	endfor
	if (j0==j1)						// a single number range
		out += num2istr(j0)
	else								// more then one in a dash range
		out += num2istr(j0)+"-"+num2istr(j1)
	endif
	return out
End

//  ================================ End of String Ranges ================================  //
//  ======================================================================================  //


//  ======================================================================================  //
//	 =============================== Start of Progress Panel =============================== //

Function/T ProgressPanelStart(percentName,[stop,showTime,status,wname,part,total])	// display a progress bar
//	if the user does a "DoUpdate/W=$wname", then V_flag will be set to 2 if the stop button was pushed
	String percentName									// name of global value to use, if not a global variable then set value to 0
	Variable stop											// if true, include a stop button
	Variable showTime									// if true, include the ending time on panel
	String status											// if passed, make room for a status line and set it to status
	String wname											// if Progress window named wname exists, exit without doing anything
	Variable part,total									// if entered, currently doing part 'part' of 'total' parts, used with wname, part starts at 0
	stop = ParamIsDefault(stop) ? 0 : stop
	stop = NumberByKey("IGORVERS",IgorInfo(0))>=6.10 ? stop : 0	// stop feature requires version 6.10 or above
	showTime = ParamIsDefault(showTime) ? 0 : showTime
	part = ParamIsDefault(part) ? 0 : part
	total = ParamIsDefault(total) ? 1 : total
	Variable showStatus=1
	if (ParamIsDefault(status))
		status = ""
		showStatus = 0
	endif
	if (ParamIsDefault(wname))
		wname = ""
	endif

	if (WinType(wname)!=7 || strlen(wname)==0)	// passed wname is not a panel, create new one
		Variable top = 87 + 93*mod(CountProgressWindows(),7)
		Variable right = stop ? 695 : 630,  bottom=top+30, new=0
		bottom += showTime ? 20 : 0
		bottom += showStatus ? 20 : 0
		NewPanel/K=1/W=(330,top,right,bottom)/N=Progress
		wname = S_name									// save the name of the window
		new = 1
	endif

	if (!new)												// panel already exists
		DoWindow/F $wname								// do not create any buttons or display fields, just bring to front
		return wname
	endif

	// create a new panel
	ValDisplay progressValdisp, win=$wname, pos={5,5},size={288,25},font="Helvetica",fSize=18,format="%d%%"
	ValDisplay progressValdisp, win=$wname, limits={0,100,0},barmisc={0,50},highColor= (19368,32650,65196)
	if (exists(percentName)==2)						// the global exists, use it
		NVAR percent = $percentName
		percent = limit(percent,0,100)
		ValDisplay progressValdisp,win=$wname, value= #percentName
		SetWindow $wname userdata(percent)=num2str(percent)
	else													// no global, just set to zero
		ValDisplay progressValdisp,win=$wname, value= #"0"
		SetWindow $wname userdata(percent)="0"
	endif

	if (stop)
		Button stopButton,pos={299,4},size={50,20},title="Break"
	endif

	top = 29
	if (showStatus)
		TitleBox statusTitle,pos={1,top},size={213,16},fSize=12,frame=0,title=status
		top += 20
	endif
	if (showTime)
		TitleBox timeTitle,pos={1,top},size={213,16},fSize=12,frame=0,title=" "
	endif

#if NumberByKey("IGORVERS",IgorInfo(0))>=6.10
	DoUpdate/W=$wname/E=1							// mark as a progress window
#else
	DoUpdate
#endif
	String str = ""
	sprintf str,"%.3f",stopMSTimer(-2)*1e-6
	SetWindow $wname userdata(startTimeSec)=str
	SetWindow $wname userdata(lastTimeSec)=str		// last time updated
	sprintf str,"%.3f",stopMSTimer(-2)*1e-6 - DateTime
	SetWindow $wname userdata(epoch)=str			// needed to get correct time of day
	SetWindow $wname userdata(showTime)=num2istr(showTime)
	SetWindow $wname userdata(status)=status
	SetWindow $wname userdata(stop)=num2istr(stop)
	SetWindow $wname userdata(startingFunc)=GetRTStackInfo(2)
	SetWindow $wname userdata(part)=num2str(part)
	SetWindow $wname userdata(total)=num2str(total)
//printf "starting function = '%s'\r",GetRTStackInfo(2)
	return wname
End
//
Function ProgressPanelUpdate(wname,percent,[status,resetClock])	// update a progress bar
	String wname											// window name
	Variable percent										// percent done (% !!)
	String status											// if passed, update the status line
	Variable resetClock									// reset start time to now if TRUE

	Variable newStatus = !ParamIsDefault(status)
	if (newStatus)
		SetWindow $wname userdata(status)=status
	else
		status = GetUserData(wname,"","status")
	endif

	Variable part = str2num(GetUserData(wname,"","part"))
	Variable total = str2num(GetUserData(wname,"","total"))
	percent = 100*part/total + percent/total			// percent ot total done

	String str
	if (resetClock)
		sprintf str,"%.3f",stopMSTimer(-2)*1e-6
		SetWindow $wname userdata(startTimeSec)=str	// reset all times to now
		SetWindow $wname userdata(lastTimeSec)=str
	endif

	ValDisplay progressValdisp,win=$wname, value= #num2istr(floor(limit(percent,0,100)))
	Variable showTime = str2num(GetUserData(wname,"","showTime"))
	if (showTime && percent>0)
		Variable start = str2num(GetUserData(wname,"","startTimeSec"))
		Variable last = str2num(GetUserData(wname,"","lastTimeSec"))
		Variable now = stopMSTimer(-2)*1e-6
		if (now-last > 2)
			Variable epoch = str2num(GetUserData(wname,"","epoch"))
			Variable remain = (now-start) * (100-percent)/percent
			Variable ifmt = (remain < 6*60)
			sprintf str, "should end in %s,   at    %s",Secs2Time(remain,4+ifmt),Secs2Time(now+remain-epoch,ifmt)
			str += SelectString(remain > 24*3600,"","  on  "+Secs2Date(now+remain,2))	// more than 1 day
			TitleBox timeTitle, win=$wname, title=str
			sprintf str,"%.3f",now
			SetWindow $wname userdata(lastTimeSec)=str		// last time updated, is now
		endif
	endif
	if (strlen(status) || newStatus)
		TitleBox statusTitle, win=$wname, title=status
	endif
	if (strlen(wname))
		DoWindow/F $wname								// keep progress window in front
	endif

	Variable done=0
#if NumberByKey("IGORVERS",IgorInfo(0))>=6.10
	DoUpdate/W=$wname									// update only progress window
	done = V_Flag == 2									// stop button on progress window was pushed
#else
	DoUpdate												// update all windows, no stop button
#endif
	return done
End
//
Function SecondsInProgressPanel(wname)				// total time in the Progress Panel
	String wname											// window name
	Variable start = str2num(GetUserData(wname,"","startTimeSec"))
	return (stopMSTimer(-2)*1e-6 - start)
End

Function ProgressPanelKill(wname)					// kills a progress bar, will not kill it if  it was created by a routine higher up in the current stack
	String wname											// window name

	String stackList=GetRTStackInfo(0)
	String startingFunc=GetUserData(wname,"","startingFunc")
	Variable i = WhichListItem(startingFunc,stackList)
	Variable seconds = SecondsInProgressPanel(wname)
	if (i<0 || i >= ItemsInList(stackList)-2)
		DoWindow/K $wname
	endif
	return seconds
End
//
Function ProgressPanelNewPart(wname,part,total)	// update part and/or total in the progress window, only updates vlaues if they are valid
	String wname											// window name
	Variable part,total									// if entered, currently doing part 'part' of 'total' parts, used with wname, part starts at 0

	if (WinType(wname)!=7)								// window not a panel
		return 1
	endif
	if (part>=0)											// part is valid, update it
		SetWindow $wname userdata(part)=num2str(part)
	endif
	part = numtype(part) ? 0 : part
	if (part<total)										// total is valid update it
		SetWindow $wname userdata(total)=num2str(total)
	endif
	return 0
End

Static Function CountProgressWindows()			// returns number of open ProgressWindows
	String str, wlist=WinList("Progress*",";","WIN:64")
	Variable i, N
	for (i=0,N=0; i<ItemsInList(wlist); i+=1)
		str = GetUserData(StringFromList(i,wlist),"","startTimeSec")
		N += strlen(str) ? 1 : 0
	endfor
	return N
End
//Function testProgress()
//	String progressWin = ProgressPanelStart("",stop=1,showTime=1,status="starting")	// display a progress bar
//	Variable i,N=300
//	for (i=0;i<N;i+=1)
//		if (mod(i,100)==0)
//			if (ProgressPanelUpdate(progressWin,i/N*100,status="on number "+num2istr(i)))	// update progress bar
//				break											//   and break out of loop
//			endif
//		else
//			if (ProgressPanelUpdate(progressWin,i/N*100))	// update progress bar
//				break											//   and break out of loop
//			endif
//		endif
//	endfor
//	printf "total execution time = %s\r",Secs2Time(SecondsInProgressPanel(progressWin),5,0)
//	ProgressPanelUpdate(progressWin,i/N*100,status="done")
//	Sleep/S 1
//	DoWindow/K $progressWin
//End
//  =============================== End of Progress Panel ================================  //
//  ======================================================================================  //


//  ======================================================================================  //
//  ================================ Start of Generic XML ================================  //
//
//	XML support	 (occurance optionally allows selecting the the occuranceth instance of xmltag), note vectors usually delimited by a space
//
//	XMLNodeList(buf)											returns a list with all top level nodes in buf
//	XMLtagContents(xmltag,buf,[occurance])			returns the contents of xmltag
//	XMLtagContents2List(xmltag,buf,[occurance,delimiters])	returns the contents of xmltag as a list, useful for vectors in the contents
//	XMLattibutes2KeyList(xmltag,buf)					return a list with all of the attribute value pairs for xmltag
//	XMLremoveComments(str)									remove all xml comments from str
//
//	for XMLtagContents() and XMLattibutes2KeyList()
// when there are MANY occurances of xmltag, do not use occurance, but rather:
//	Variable start=0
//	String feed
//	do
//		feed = XMLtagContents("feed",buf, start=start)
//		other code goes here ...
//	while(strlen(feed))
//
//
//		Element names are case-sensitive
//		Element names must start with a letter or underscore
//		Element names cannot start with the letters xml (or XML, or Xml, etc)
//		Element names can contain letters, digits, hyphens, underscores, and periods
//		Element names cannot contain spaces

ThreadSafe Function/T XMLNodeList(buf)				// returns a list of node names at top most level in buf
	String buf
	String tName, names=""
	Variable i0=0, i1,i2
	do
		i0 = XMLfindNextLeadingTagPos(buf,i0)			// i0 is on the leading "<"
		if (i0<0)
			break
		endif
		tName = XMLfindNextLeadingTagName(buf,i0)
		if (XMLvalidElementName(tName))					// skip prolog and things like that
			names += ReplaceString(";",tName,"_")+";"	// tName cannot contain semi-colons
		endif
		i0 = XMLfindCloser(buf,tName,i0+1)				// position just after closer
	while(i0>0)
	return names
End


ThreadSafe Function/T XMLtagContents(xmltag,buf,[occurance,start,utf8])
	String xmltag
	String buf
	Variable occurance				// use 0 for first occurance, 1 for second, ...
	Variable &start					// offset in buf, start searching at buf[start], new start is returned
	Variable utf8					// flag, if true convert html codes to utf8, e.g. "&#x03B2;" --> beta, only for Igor>=7
										// both occurance and start may be used together, but usually you only want to use one of them
	occurance = ParamIsDefault(occurance) ? 0 : occurance
	Variable startLocal = ParamIsDefault(start) ? 0 : start
	startLocal = numtype(startLocal) || startLocal<1 ? 0 : round(startLocal)
	utf8 = ParamIsDefault(utf8) || numtype(utf8) ? 0 : utf8
	if (!XMLvalidElementName(xmltag))
		return ""
	endif

	Variable i1, i0=startOfxmltag(xmltag,buf,occurance, start=startLocal)
	if (i0<0)
		return ""
	endif
	i0 = strsearch(buf,">",i0)						// character after '>' in intro
	if (i0<0)												// this is an ERROR
		return ""
	endif
	i0 += 1													// start of contents

	i1 = XMLfindCloser(buf,xmltag,i0-2)			// character just after final '>'
	startLocal = i1										// save character just after closing '<tag>'

	if (i1<i0 || i1<0)									// could not find a valid closer
		if (!ParamIsDefault(start))
			start = -1
		endif
		return ""
	endif

	if (!ParamIsDefault(start))						// pass start back if used
		start = startLocal
	endif

	i1 -= strlen(xmltag)+4
	if (i1<i0)
		return ""										// no content
	endif

	buf = buf[i0,i1]
	if (utf8 && IgorVersion()>=7)
		buf = ConvertHTML2UTF8(buf)					// convrert html codes to UTF8 if requested
	endif
	return buf
End


ThreadSafe Function/T XMLtagContents2List(xmltag,buf,[occurance,delimiters]) //reads a tag contents and converts it to a list
	String xmltag
	String buf
	Variable occurance				// use 0 for first occurance, 1 for second, ...
	String delimiters					// characters that might be used for delimiters (NOT semi-colon), default is space, tab, cr, or nl = " \t\r\n"
	occurance = ParamIsDefault(occurance) ? 0 : occurance
	if (ParamIsDefault(delimiters) || strlen(delimiters)==0)
		delimiters = WhiteSpaceChars											// the usual white-space characters
	endif

	String str = XMLtagContents(xmltag,buf,occurance=occurance)	// get the contents
	if (strlen(str)<1)
		return ""
	endif
	return XMLContents2List(str,delimiters=delimiters)				// convert contents to a list
End
//
ThreadSafe Static Function/T XMLContents2List(str,[delimiters,sep])
	// take a tag contents (just text, no xml) and convert it to a semi-colon separated list
	// assumes that list elements can be quoted with single or double quotes
	// list elements are separated by characters in delimiters (CANNOT be semi-colon)
	String str							// contents of an xml tag pair, the contents in: <tag>contents</tag>
	String delimiters					// characters that might be used for delimiters (NOT semi-colon), default is space, tab, cr, or nl = " \t\r\n"
	String sep							// list separator, usually ";"
	if (ParamIsDefault(delimiters) || strlen(delimiters)==0)
		delimiters = WhiteSpaceChars					// the usual white-space characters
	endif
	if (ParamIsDefault(sep) || strlen(sep)==0)
		sep = ";"											// the usual white-space characters
	endif
	String sep2 = sep+sep

	str = ReplaceString(sep,str,"_")				// cannot have any semi-colons in input string, that is the list separator

	String first, list=""
	Variable i0,i1
	do
		first = str[0]
		if (strsearch(delimiters,first,0)>=0)		// starts with a delimiter, skip this character
			str = str[1,Inf]
			continue
		elseif (cmpstr(first,"\"")==0)				// starts with a ", find the matching "
			i0 = 1
			i1 = strsearch(str,"\"", 1)-1
		elseif (cmpstr(first,"'")==0)				// starts with a ', find the matching '
			i0 = 1
			i1 = strsearch(str,"'", 1)-1
		else
			i0 = 0											// no quotes, find next delimiter
			i1 = FindFirstCharInStr(str[1,Inf],delimiters+"\"'")
		endif

		if (i1>=i0)
			list += str[i0,i1]+sep						// add to list
			str = str[i1+2,Inf]
		else
			break
		endif
	while(numtype(i1)==0)

	do
		list = ReplaceString(sep2,list,sep)		// replace all multiple semi-colons with a single semi-colon
	while(strsearch(list,sep2,0)>=0)
	list = RemoveLeadingString(list,sep,1)		// remove any leaing semi-colon

	return list
End


ThreadSafe Function/T XMLattibutes2KeyList(xmltag,buf,[occurance,start])// return a list with all of the attribute value pairs for xmltag
	String xmltag											// name of tag to find
	String buf												// buf containing xml
	Variable occurance									// use 0 for first occurance, 1 for second, ...
	Variable &start										// offset in buf, start searching at buf[start], new start is returned
																// both occurance and start may be used together, but usually you only want to use one of them
	occurance = ParamIsDefault(occurance) ? 0 : occurance
	Variable startLocal = ParamIsDefault(start) ? 0 : start
	startLocal = numtype(startLocal) || startLocal<1 ? 0 : round(startLocal)
	if (!XMLvalidElementName(xmltag))
		return ""
	endif

	Variable i1, i0=startOfxmltag(xmltag,buf,occurance, start=startLocal)
	if (i0<0)
		return ""
	endif
	i0 += strlen(xmltag)+2								// first char of the attributes
	i1 = strsearch(buf,">",i0)-1						// last char of the attributes

	String key, value, keyVals=""
	if (i1 < i0)											// this is an ERROR
		startLocal = -1
	else
		i1 -= (cmpstr(buf[i1],"/")==0) ? 1 : 0			// in case of a short tag, i.e. ends with "/>"
		startLocal = XMLfindCloser(buf,xmltag,i1+1)	// get character just AFTER closing '>'
		// parse buf into key=value pairs
		buf = buf[i0,i1]
		buf = ReplaceString("\t",buf," ")
		buf = ReplaceString("\r",buf," ")
		buf = ReplaceString("\n",buf," ")
		buf = TrimFrontBackWhiteSpace(buf)
		i0 = 0
		do
			i1 = strsearch(buf,"=",i0,0)
			key = TrimFrontBackWhiteSpace(buf[i0,i1-1])
			i0 = strsearch(buf,"\"",i1,0)+1			// character after the first double quote around value
			i1 = strsearch(buf,"\"",i0,0)-1			// character before the second double quote around value
			value = buf[i0,i1]
			if (strlen(key)>0)
				keyVals = ReplaceStringByKey(key,keyVals,value,"=")
			endif
			i0 = strsearch(buf," ",i1,0)				// find space separator, set up for next key="val" pair
		while(i0>0 && strlen(key))
	endif

	if (!ParamIsDefault(start))						// set start if it was passed
		start = startLocal
	endif
	return keyVals
End


ThreadSafe Function/T XMLremoveComments(str)	// remove all xml comments from str
	String str
	Variable i0,i1
	do
		i0 = strsearch(str,"<!--",0)					// start of a comment
		i1 = strsearch(str,"-->",0)					// end of a comment
		if (i0<0 || i1<=i0)
			break
		endif
		str[i0,i1+2] = ""									// snip out comment
	while(1)
	return str
End


ThreadSafe Static Function startOfxmltag(xmltag,buf,occurance,[start])	// returns the index into buf pointing to the start of xmltag
	String xmltag, buf
	Variable occurance									// use 0 for first occurance, 1 for second, ...
	Variable start
	start = ParamIsDefault(start) || start<=0 || numtype(start) ? 0 : round(start)

	String name												// current tag
	Variable i0,i1, iocc, match
	for (iocc=0,i0=start; iocc<=occurance; )
		i0 = XMLfindNextLeadingTagPos(buf,i0)
		if (i0<0)
			return -1
		endif
		name = XMLfindNextLeadingTagName(buf,i0)

		match = cmpstr(xmltag,name,1)==0
		iocc += match
		if (!match || iocc<=occurance)
			i0 = XMLfindCloser(buf,name,i0+1)		// start after the beginning of the opener
		endif
	endfor
	return i0
End
//
ThreadSafe Static Function/T XMLfindNextLeadingTagName(buf,start)
	// return string with the first tag that we find in buf, start searching in buf[start]
	String buf												// a string containing xml
	Variable start											// start searching at buf[start]
	if (start<0)
		return ""
	endif

	Variable i0, i1,i2
	i0 = XMLfindNextLeadingTagPos(buf,start)
	if (i0 < 0)
		return ""
	endif

	i0 += 1													// first char in name of tag

	i1 = strsearch(buf," ",i0) - 1
	i2 = strsearch(buf,">",i0) - 1
	i1 = i1<i0 ? Inf : i1
	i2 = i2<i0 ? Inf : i2
	i1 = min(i1,i2)

	if (numtype(i1)==0)
		return buf[i0,i1]
	endif
	return ""
End
//
ThreadSafe Static Function XMLfindNextLeadingTagPos(buf,start)
	// return position of "<" that is the next xml tag in buf[start]
	String buf												// a string containing xml
	Variable start											// start searching at buf[start]
	if (start<0)
		return -1
	endif

	return strsearch(buf,"<",start)
End
//
ThreadSafe Static Function XMLfindCloser(buf,xmltag,start)
	// returns index into buf that points to first character AFTER "</xmltag>"
	String buf												// a string containing xml
	String xmltag											// the name of the xml tag
	Variable start											// start searching buf at buf[start]

	Variable i
	if (char2num(xmltag[0])==63)							// prolog, tag starts with a "?"==63
		i = strsearch(buf,"?>",start)
		i += i>0 ? 2 : 0
		return i
	endif

	String close1="</"+xmltag+">", close2="/>"
	if (cmpstr(xmltag,"!--")==0)							// special rule for comments
		close1 = "-->"
		close2 = close1
	endif
	Variable j0,j1, i1,i2, iClose						// jn are openers, in are closers

	do
		j0 = strsearch(buf,"<"+xmltag+" ",start)	// next opener of type "<tag ...>"
		j0 = j0<0 ? Inf : j0
		j1 = strsearch(buf,"<"+xmltag+">",start)	// next opener of type "<tag>"
		j1 = j1<0 ? Inf : j1
		j0 = min(j0,j1)										// next opener of tag
		// j0 is position of next tag opener (j1 no longer used)

		i1 = strsearch(buf,close1,start)				// next closer of type "</tag>"
		i1 = i1<0 ? Inf : i1
		i2 = strsearch(buf,close2,start)				// next closing of type "/>"
		i2 = i2<0 ? Inf : i2
		i2 = i2<strsearch(buf,"<",start) ? i2 : Inf	// required for sort tags

		// j0 is position of next tag opener (ignore j1)
		//	i1 is next "</tag>"
		//	i2 is next "/>" (it is inf if a </tag> comes first

		if (j0<i1 && j0<i2)									// found another open before the close, dig deeper
			start = XMLfindCloser(buf,xmltag,j0+strlen(xmltag))
			iClose = NaN										// NaN causes a loop in the while()
		elseif (i1<i2)
			iClose = i1 + strlen(close1)					// done, advance to char after "</xmltag>"
		elseif (i2<i1)
			iClose = i2 + strlen(close2)					// done, advance to char after "/>"
		else
			iClose = -1											// failed to find a closer
		endif
	while (numtype(iClose))
	return iClose
End
//
//	ThreadSafe Static Function XMLfindCloser(buf,xmltag,start)
//		// returns index into buf that points to first character AFTER "</xmltag>"
//		String buf												// a string containing xml
//		String xmltag											// the name of the xml tag
//		Variable start											// start searching buf at buf[start]
//	
//		Variable i, len=strlen(xmltag)+3
//	
//		if (char2num(xmltag[0])==63)						// tag starts with a "?"==63
//			i = strsearch(buf,"?>",start)
//			i += i>0 ? 2 : 0
//			return i
//		endif
//	
//		i = strsearch(buf,"</"+xmltag+">",start)	// find the closer for this tag, search for '</xmltag>'
//		if (i<0)
//			i = strsearch(buf,"/>",start)				// no '</xmltag>', just an empty-element tag, ends with '/>'
//			len = 2
//		endif
//	
//		if (i<0)													// failed to find closer
//			return -1
//		endif
//		return i+len
//	End


ThreadSafe Static Function XMLvalidElementName(tagName)
	// returns 1 if tagName is a valid xml name, 0 if not valid
	String tagName
	if (!isletter(tagName) && char2num(tagName)!=95)	// must start with letter or '_'
		return 0
	endif
	if (strsearch(tagName," ",0) >= 0)						// cannot contain " "
		return 0
	endif
	if (StringMatch(tagName,"xml*"))						// cannot start with xml, XML, Xml, etc.
		return 0
	endif
	return 1	
End


ThreadSafe Function/S ConvertHTML2UTF8(buf)
	// replace all occurances of somethng like "&#x03B2;" or "&#946;" with the UTF8 character
	String buf			// buffer containing html codes
	String html, utf8
	Variable i0=strsearch(buf,"&#",0,2), i1=strsearch(buf,";",i1,2)
	for(; i0>=0 && i1>i0; )
		html = buf[i0,i1]							// the html string, something like "&#x03B2;" or "&#946;"
		utf8 = HTMLcode2UTF8(html)				// a single UTF8 character
		buf = ReplaceString(html, buf, utf8)	// replace all occurances of html --> utf8
		i0 = strsearch(buf,"&#",i0+1,2)		// find the next html code to replace, start searching after the one replaced
		i1 = strsearch(buf,";",i0,2)
	endfor
	return buf
End
//
ThreadSafe Static Function/S HTMLcode2UTF8(code)
	// Converts a string like "&#x03B2;" or "&#946;" to UTF8 character
	String code

	Variable lead = strsearch(code,"&#",0)
	Variable last = strsearch(code,";",0)
	String num = code[2,last-1]

	Variable ichar
	if (strlen(code)>10)						// code cannot be longer than 10
		return code
	elseif (lead != 0 || last-2 != strlen(num))
		return code
	elseif (strsearch(num,"x",0,2)==0)
		sscanf num, "x%x;", ichar			// Hex code
	else
		sscanf num, "%d;", ichar				// Decimal code
	endif
	ichar = V_flag==1 ? ichar : 0

	if (numtype(ichar)==0 && ichar>0 && ichar<=1114112)	// max value for ichar is 17*(2^16) = 1114112
		return num2char(ichar)
	endif
	return code
End

//  ================================= End of Generic XML =================================  //
//  ======================================================================================  //


//  ======================================================================================  //
//  ============================ Start WaveClass in Wave Note ============================  //

// return a list of waves in current folder having a "waveClass" that is a member of the list waveClassList
// The waveClassList, is a semicolon separated list, and the members can have wildcards. e.g. "speImage*"
// The class of a wave is given by a key=value pair in the wavenote:      "key1=val1;waveClass=class1,class2,class3;key5=val5;"
// This is similar to WaveList(), but with a finer selection
//
// returns true if any one of the classes of ww matches one of the classes in waveClassList
// note that the items in waveClassList can have wild cards
Function/T WaveListClass(waveClassList,search,options,[all,win,fldr])
	String waveClassList				// a list of acceptable wave classes (semicolon separated)
	String search						// same as first argument in WaveList()
	String options						// same as last argument in WaveList()
	Variable all						// when all is TRUE, then all of the classes in waveClassList must be present, not just one
	String win							// when present, only provide waves also displayed in win
	String fldr							// name of data folder to check, default is current folder
	all = ParamIsDefault(all) ? 0 : !(!all)
	all = numtype(all) ? 0 : all
	win = SelectString(ParamIsDefault(win),win,"")
	fldr = SelectString(ParamIsDefault(fldr),fldr,"")
	if (strlen(fldr))
		fldr += SelectString(StringMatch(fldr,"*:"),":","")	// ensure fldr ends with ":"
		fldr = SelectString(strsearch(fldr,"root:",0) && strsearch(fldr,":",0),"",":")+fldr	// add leading ":" unless "root:..."
		if (!DataFolderExists(fldr))	// if requested folder does not exist, just return
			return ""
		endif
	endif

	String fldrSav=GetDataFolder(1)
	SetDataFolder fldr
	String in = WaveList(search,";",options), name, out=""
	SetDataFolder fldrSav

	Variable m, displayed=1
	for (m=0, name=StringFromList(0,in); strlen(name); m+=1,name=StringFromList(m,in))
		name = fldr+name
		if (strlen(win))
			if (WinType(win)==GIZMO_WIN_TYPE)		// a gizmo, cannot use CheckDisplayed
				displayed = isWaveOnGizmo(win, $name)
			else
				CheckDisplayed/W=$win $name
				displayed = !(!V_flag)
			endif
		endif
		if (displayed && WaveInClass($name,waveClassList,all=all))
			out += name+";"
		endif
	endfor
	return out
End


// removes waves from inList that are also in includeClassList
Function/T IncludeOnlyWavesInClass(inList,includeClassList)
	String inList						// a list of input waves (semicolon separated)
	String includeClassList				// a list of classes to exclude (semicolon separated)

	if (strlen(inList)<1)				// nothing
		return ""
	endif
	String out=""
	String name
	Variable m, use
	for (m=0, name=StringFromList(0,inList); strlen(name); m+=1,name=StringFromList(m,inList))
		Wave ww=$name
		if (WaveExists(ww))
			if (WaveInClass(ww,includeClassList))	// keep this wave
				out += name+";"
			endif
		endif
	endfor
	return out
End


// returns true if any one of the classes of ww matches one of the classes in waveClassList
// note that the items in waveClassList can have wild cards
ThreadSafe Function WaveInClass(ww,waveClassList,[all])
	Wave ww						// Wave to check
	String waveClassList			// a list of acceptable wave classes (semicolon separated)
	Variable all					// when all is TRUE, then all of the classes in waveClassList must be present, not just one
	if (!WaveExists(ww))
		return 0
	elseif (strlen(waveClassList)==1 && char2num(waveClassList)==42)		// just a "*" matches everything
		return 1
	elseif (strlen(waveClassList)<1)
		return 1
	endif

	waveClassList = ReplaceString(",",waveClassList,";")
	all = ParamIsDefault(all) ? 0 : !(!all)
	all = numtype(all) ? 0 : all
	String class = StringByKey("waveClass",note(ww),"=")	// class list stored in wave note (comma separated)
	String matchClass, wavClass
	Variable m, i,Nc, Nclasses=ItemsInList(waveClassList)
	Variable Nmatch=(all ? Nclasses : 1)						// number of classes I have to match
	for (m=0,Nc=0; m<Nclasses; m+=1)						// check each item in waveClassList
		matchClass = StringFromList(m,waveClassList)
		for (i=0;i<ItemsInList(class,",");i+=1)
			wavClass = StringFromLIst(i,class,",")				// class item from the wave
			Nc += stringmatch(wavClass, matchClass)			// note that matchClass can have wild cards
			if (Nc>=Nmatch)
				return 1
			endif
		endfor
	endfor
	return 0
End


// Adds a class to a waveClass
ThreadSafe Function/T AddClassToWaveNote(wnote,addClass)
	String wnote									// existing wavenote with key=value parirs
	String addClass								// class is a comma separated list of wave classes (or just one)

	String waveClasses=StringByKey("waveClass",wnote,"=")	// get current list of wave classes
	String oneClass								// one of the classes from addClass
	Variable i
	for (i=0;i<ItemsInList(addClass,",");i+=1)
		oneClass = StringFromList(i,addClass)
		if (WhichListItem(oneClass,waveClasses,",")<0)					// addClass not in waveClasses, add it
			waveClasses = AddListItem(oneClass,waveClasses,",",Inf)
		endif
	endfor
	wnote = ReplaceStringByKey("waveClass",wnote,waveClasses,"=")	// fix up wave note
	return wnote
End


// Remove a class(s) from the waveClass in a wavenote
Function/T RemoveClassFromWaveNote(wnote,removeClass)
	String wnote							// existing wavenote with key=value parirs
	String removeClass					// class is a COMMA separated list of wave classes (or just one) to remove

	String waveClasses=StringByKey("waveClass",wnote,"=")	// get current list of wave classes
	String oneClass								// one of the classes from removeClass
	Variable i
	for (i=0;i<ItemsInList(removeClass,",");i+=1)
		oneClass = StringFromList(i,removeClass)	// each class in removeClass
		waveClasses = RemoveFromList(oneClass,waveClasses,",")
	endfor

	for (i=strlen(waveClasses)-1; char2num(waveClasses[i])==44 && i>=0; i-=1)	// find last not comma
	endfor											// trim off trailing commas
	waveClasses = waveClasses[0,i]
	wnote = ReplaceStringByKey("waveClass",wnote,waveClasses,"=")	// fix up wave note
	return wnote
End


// removes waves from inList that are also in excludeClassList
ThreadSafe Function/T ExcludeWavesInClass(inList,excludeClassList)
	String inList						// a list of input waves (semicolon separated)
	String excludeClassList			// a list of classes to exclude (semicolon separated)

	if (strlen(inList)<1)				// nothing
		return ""
	endif
	String out=""
	String name
	Variable m, use
	for (m=0, name=StringFromList(0,inList); strlen(name); m+=1,name=StringFromList(m,inList))
		Wave ww=$name
		if (WaveExists(ww))
			if (!WaveInClass(ww,excludeClassList))	// keep this wave
				out += name+";"
			endif
		endif
	endfor
	return out
End


// return list of windows that display a wave of appropriate waveClass
Function/S WindowsWithClass(classes,search,options,flag)
	String classes						// a list of acceptable wave classes (semicolon separated)
	String search						// same as first argument in WaveList()
	String options						// same as last argument in WaveList()
	Variable flag						// bit flag: 1=graphs, 2=tables, 4=gizmos, or any sum of 1,2,4, 7 or -1 gives all

	Variable bit = flag & 3 + (flag & 4 ? GIZMO_WIN_BIT : 0)
	String allWindows=WinList("*",";","WIN:"+num2istr(bit))

	String win, wList=""
	Variable i, N=ItemsInList(allWindows)
	for (i=0;i<N;i+=1)
		win = StringFromList(i,allWindows)
		if (strlen(WaveListClass(classes,search,options,win=win)))
			wList += win+";"
		endif
	endfor
	return wList
End


// returns list of folders containing waves with waveClassList
Function/T FoldersWithWaveClass(fldrPath,waveClassList,search,options,[all,win,fullPath])
	String fldrPath					// folder in which to look, defaults to current folder
	String waveClassList			// a list of acceptable wave classes (semicolon separated)
	String search						// same as first argument in WaveList()
	String options						// same as last argument in WaveList()
	Variable all						// when all is TRUE, then all of the classes in waveClassList must be present, not just one
	String win							// when present, only provide waves also displayed in win
	Variable fullPath					// flag, when true, returned list has full path
	all = ParamIsDefault(all) ? 0 : !(!all)
	all = numtype(all) ? 0 : all
	win = SelectString(ParamIsDefault(win),win,"")
	fullPath = ParamIsDefault(fullPath) || numtype(fullPath) ? 0 : !(!fullPath)
	fldrPath = SelectString(strlen(fldrPath),":",fldrPath)
	fldrPath = FixUpFolderName(fldrPath)	// add necessary ":"

	DFREF DFR = $fldrPath
	Variable Nlist = CountObjectsDFR(DFR,4)
	String name, flist, list=""
	Variable i
	for (i=0;i<Nlist;i+=1)		// check each folder in fldrPath
		name = GetIndexedObjNameDFR(DFR,4,i)
		if (ItemsInList(WaveListClass(waveClassList,search,options,all=all,win=win,fldr=fldrPath+name)))
			list += SelectString(fullPath,"",fldrPath) + name + ";"
		endif
	endfor
	return list
End


Function/T TraceNamesInClass(waveClassList,win,[optionsFlag,Xwaves])
	// like TraceNameList(), but limit result to a list of wave classes
	// returns a list of acceptable trace names, use TraceNameToWaveRef() to get wave ref.
	String waveClassList			// list of classes, semi-colon separated
	String win
	Variable optionsFlag
	Variable	Xwaves					// flag, if true, also include Xwaves
	optionsFlag = ParamIsDefault(optionsFlag) || numtype(optionsFlag) ? 1 : optionsFlag	// only normal graph traces (exclude contours & hidden)
	Xwaves = ParamIsDefault(Xwaves) || numtype(Xwaves) ? 0 : !(!Xwaves)
	String listAll=TraceNameList(win,";",optionsFlag), list="", trName

	Variable i, N=ItemsInlist(listAll)
	for (i=0;i<N;i+=1)			// for each trace, check if it is in waveClassList
		trName = StringFromList(i,listAll)
		Wave ww = TraceNameToWaveRef(win, trName)
		if (WhichListItem(trName,list)<0)		// not already in list, add it
			list += SelectString(WaveInClass(ww,waveClassList), "", trName+";")
		endif
		if (Xwaves)
			Wave ww = XWaveRefFromTrace(win, trName)
			trName = NameOfWave(ww)
			if (WhichListItem(trName,list)<0)		// not already in list, add it
				list += SelectString(WaveInClass(ww,waveClassList), "", trName+";")
			endif
		endif
	endfor
	return list
End

//  ============================= End WaveClass in Wave Note =============================  //
//  ======================================================================================  //


//  ======================================================================================  //
//  ====================== Start of some general utility functions =======================  //


Proc RecompileAllProcedures()						// FORCE ALL procedures to recompile,  This must be a Proc or Macro (NOT Function)
	SetIgorOption poundDefine=DOESNTMATTER		// mark all procedures as needing compile 
	SetIgorOption poundUnDefine=DOESNTMATTER	// don't leave this defined.
	Execute/P/Q/Z "COMPILEPROCEDURES "				// re-compile (all)
	// print "ran RecompileAllProcedures"
End


Function/T FixUpFolderName(fldr)		// return foldr name with the needed ":", so fldr+"wavename" is valid
	String fldr
	if (strlen(fldr)<1)
		return ""
	endif

	if (strsearch(fldr,"root:",0,2) != 0)				// if it does not start with "root:" then it must start with ":"
		fldr = SelectString(strsearch(fldr,":",0), "", ":") + fldr
	endif
	Variable i = strlen(fldr)-1
	if (i>0 && char2num(fldr[i])!=char2num(":"))	// it must end with a ":"
		fldr += ":"
	endif
	return fldr
End


Function/T AddEndingToWaveName(wName,waveNameEnd)
	// create a new wave name from wav but with a new ending, needed due to maxIgorWaveNameLen
	// also call CleanupName() so the returned name is a valid wave name
	String wName
	String waveNameEnd

	String path=""								// a possible path preceeding the name
	Variable i=strsearch(wName,":",Inf,1)
	if (i>=0)
		path = wName[0,i]						// strip off the path and save it
		wName = wName[i+1,Inf]				// wName, now only the name part (no path)
	endif
	wName = wName[0,maxIgorWaveNameLen-strlen(waveNameEnd)-1]
	wName = path + CleanupName(wName + waveNameEnd,0)	// reassemble the full name
	return wName
End


Function/T WavesWithMatchingKeyVals(inList,keyVals)
	String inList						// a list of input waves (semicolon separated)
	String keyVals						// key value pairs that must match for acceptance (optional)

	if (strlen(inList)<1)				// nothing
		return ""
	elseif (strlen(keyVals)<1)		// nothing to do
		return inList
	endif

	String out=""
	String name, item, key,val
	Variable m, i
	for (m=0, name=StringFromList(0,inList); strlen(name); m+=1,name=StringFromList(m,inList))
		Wave ww=$name
		if (WaveExists(ww))
			for (i=0;i<ItemsInlist(keyVals);i+=1)
				item = StringFromList(i,keyVals)
				key = StringFromList(0,item,"=")
				val = StringFromList(1,item,"=")
				if (StringMatch(StringByKey(key,note(ww),"="), val))
					out += name+";"
				endif
			endfor
		endif
	endfor
	return out
End


Function/T OnlyWavesThatAreDisplayed(inList,[not])
	String inList						// a list of input waves (semicolon separated)
	Variable not						// if TRUE, then only return wave from inList that are NOT displayed
	not = ParamIsDefault(not) ? 0 : !(!not)

	if (strlen(inList)<1)				// nothing
		return ""
	endif
	String out=""
	String name, item
	Variable m, i, use
	for (m=0, name=StringFromList(0,inList); strlen(name); m+=1,name=StringFromList(m,inList))
		Wave ww=$name
		if (WaveExists(ww))
			use = strlen(WindowsWithWave(ww,1)) > 0
			use = not ? !use : use
			out += SelectString(use,"",name+";")
		endif
	endfor
	return out
End


ThreadSafe Function keyInList(key,keyWordList,keySepStr,listSepStr)	// returns true if key=value pair is in the keyWordList
	String key						// string with key
	String keyWordList			// list of keyword=value pairs
	String keySepStr				// separates key and value, defaults to colon
	String listSepStr				// separates key value pairs, defaults to semicolon
	keySepStr = SelectString(strlen(keySepStr),":",keySepStr)	// default to colon
	listSepStr = SelectString(strlen(listSepStr),";",listSepStr)	// default to semicolon

	String find=key+keySepStr								// find this
	if (strsearch(keyWordList,find,0)==0)				// check if it is at the start
		return 1												// found key=value is first pair
	endif
	if ( strsearch(keyWordList,listSepStr+find,0)>0)		// check if key is after first key=value pair
		return 1												// found key=value is a later pair
	endif
	return 0													// no key=value found
End


//	Merges two key=value lists, if priority=0, then list0 has priority, if priority=1 then list1
ThreadSafe Function/S MergeKeywordLists(list0,list1,priority,keySepStr,listSepStr,[keys])
	String list0,list1
	Variable priority				// 0 or 1
	String keySepStr				// separates key and value, defaults to colon
	String listSepStr				// separates key value pairs, defaults to semicolon
	String keys						// an optional list of keys to transfer ("" or "*" means use all)
										// when keys is used, list0 is added to from list1
	keySepStr = SelectString(strlen(keySepStr),":",keySepStr)	// default to colon
	listSepStr = SelectString(strlen(listSepStr),";",listSepStr)	// default to semicolon
	keys = SelectString(ParamIsDefault(keys),keys,"")					// defaults to all
	Variable check_keys = (strlen(keys)>0) && (strsearch(keys,"*",0)!=0)
	String item, key,value
	Variable i,N=ItemsInList(list1)
	for (i=0;i<N;i+=1)				// for each keyword=value pair in list1
		item = StringFromList(i,list1,listSepStr)
		key = StringFromList(0,item,keySepStr)
		value = StringFromList(1,item,keySepStr)
		if (check_keys)
			if (WhichListItem(key,keys)<0)	// checking keys, but key not in list, so skip
				continue
			endif
		endif
		if (keyInList(key,list0,keySepStr,listSepStr) && priority==0)
			continue				// skip because key already in list0, and list0 has priority
		endif
		list0 = ReplaceStringByKey(key,list0,value,keySepStr,listSepStr)
	endfor
	return list0
End


//	returns a list containing only the keys from keyVals, keySepStr defaults to ":", listSepStr to ";"
ThreadSafe Function/S keysInList(keyVals,keySepStr,listSepStr)
	String keyVals					// a key=val; list
	String keySepStr
	String listSepStr
	listSepStr = SelectString(strlen(listSepStr),";",listSepStr)
	keySepStr = SelectString(strlen(keySepStr),":",keySepStr)

	Variable i0,i, N=	ItemsInList(keyVals,listSepStr)
	String key, out=""
	for (i=0;i<N;i+=1)
		key = StringFromList(i,keyVals,listSepStr)
		i0 = strsearch(key,keySepStr,0)
		if (i0==-1)
			out += key+listSepStr
		elseif(i0>0)
			out += key[0,i0-1]+listSepStr
		endif
	endfor
	return out
End


Function/T joinLists(a,b,[sep])
	// join two lists using separator sep.  This takes care of terminated or un-terminated lists
	String a,b			// two list with sep as separator
	String sep			// defaults to ";"
	sep = SelectString(ParamIsDefault(sep), sep, ";")	// not given, use ";"
	sep = SelectString(strlen(sep)<1, sep, ";")			// empty string, use ";"
	if (strlen(a)<1)			// a is empty, just b
		return b
	elseif (strlen(b)<1)	// b is empty, just a
		return a
	else							// join a & b
		return RemoveEnding(a,sep)+ sep + b
	endif
End


Function/T intersectionOfLists(a,b,[sep])
	// return intersection of two lists using separator sep.  This takes care of terminated or un-terminated lists
	String a,b			// two list with sep as separator
	String sep			// defaults to ";"
	sep = SelectString(ParamIsDefault(sep), sep, ";")	// not given, use ";"
	sep = SelectString(strlen(sep)<1, sep, ";")			// empty string, use ";"

	if (strlen(a)<1 || strlen(b)<1)		// if a or b is empty, so is intersection
		return ""
	endif
	String out="", item
	Variable i
	for (i=0;i<ItemsInList(a);i+=1)
		item = StringFromList(i,a,sep)
		if (WhichListItem(item,b,sep)>=0)
			out += item+sep
		endif
	endfor
	return out
End


// remove duplicate items in list
Function/T RemoveDuplicatesFromList(list,[listSepStr,matchCase])
	String list
	String listSepStr				// separates key value pairs, defaults to semicolon
	Variable matchCase			// if true, match the case too (default is MATCH CASE)
	listSepStr = SelectString(ParamIsDefault(listSepStr),listSepStr,";")	// default to semicolon
	listSepStr = SelectString(strlen(listSepStr),";",listSepStr)	// default to semicolon
	matchCase = ParamIsDefault(matchCase) || numtype(matchCase) ? 1 : matchCase

	String item
	Variable i, m
	for (i=0;i<ItemsInList(list);i+=1)
		item = StringFromList(i,list)
		m = WhichListItem(item,list,listSepStr,i+1,matchCase)
		if (m>0)
			list = RemoveListItem(m,list,listSepStr)
			i -= 1					// forces a search again for item, there may be more than one occurance of item
		endif
	endfor
	return list
End


Function/S NextLineInBuf(buf,i,[NL])
	// eturns the next line in buf (which has <nl> separators) each time it is called, see the following test_NextLineInBuf() for usage
	String buf
	Variable &i					// next line starts at buf[i], this gets set for each call in here
	String NL					// a single character, probably "\n", but may be "\r"
	NL = SelectString(ParamIsDefault(NL),NL,"\n")
	NL = SelectString(strlen(NL)<1,NL,"\n")

	Variable len = strlen(buf)
	if (i<0 || i>=len)
		i = -1
		return ""
	endif

	Variable i0=i, i1
	i1 = strsearch(buf,NL,i0)
	if (i1<0)
		i = -1					// flags that we are done
		return buf[i0,Inf]
	endif
	i = i1+1						// first char after the NL
	i = i>=len ? -1 : i		// new starting point

	i1 = limit(i1,0,len-1)
	return buf[i0,i1]			// the next line
End
//	Function test_NextLineInBuf()
//		String buf = "\nabc\nxyz\n\nhij\n", line
//		Variable i=0
//		do
//			line = NextLineInBuf(buf,i)
//			printf "line = --%s--   i = %d\r",line,i
//		while(i>=0)
//	End


Function/T AxisLabelFromGraph(gName,w,axis)	// returns specified axis label for the top graph containing the wave
	String gName										// name of graph, use "" for the top graph
	Wave w											// if w does not exist, use the first wave (image or trace) in the graph gName
	String axis										// use "HORIZONTAL" or "VERTICAL" for the x and y axes

	if (strlen(WinList("*", "", "WIN:1"))<1)
		return ""
	endif
	if (!WaveExists(w))								// find top wave
		Wave w = ImageNameToWaveRef(gName,StringFromList(0,ImageNameList(gName,";")))
		if (!WaveExists(w))
			Wave w = TraceNameToWaveRef(gName,StringFromList(0,TraceNameList(gName,";",1)))
		endif
		if (!WaveExists(w))
			return ""
		endif
	endif
	gName = SelectString(strlen(gName),StringFromList(0,WindowsWithWave(w,1)),gName)
	String infoStr=""
	if (WaveExists(w))
		infoStr=ImageInfo(gName,NameOfWave(w),0)
		if (strlen(infoStr)<1)
			infoStr = TraceInfo(gName,NameOfWave(w),0)
		endif
	endif
	if (cmpstr(axis,"VERTICAL",1)==0)			// return label for vertical axis used by wave w
		axis = AxisLabelFromGraph(gName,w,StringByKey("YAXIS",infoStr))
		return axis
	elseif (cmpstr(axis,"HORIZONTAL",1)==0)	// return label for horzontal axis used by wave w
		axis = AxisLabelFromGraph(gName,w,StringByKey("XAXIS",infoStr))
		return axis
	endif

	String recreation=WinRecreation(gName,0)
	String str = "\tLabel "+axis
	Variable i0,i1
	i0 = strsearch(recreation,str,0,2)
	if (i0<0)
		return ""
	endif
	i0 += strlen(str)+1
	i1 = strsearch(recreation[i0,Inf],"\r",0,2)-1 + i0

	str=recreation[i0+1,i1-1]						// strip off leading and trailing double quotes
	str = ReplaceString("\\\\",str,"\\")
	return str
End



Function/S WindowsWithWave(w,flag)
	// returns a list of all graph, table or gizmo windows that make use of w in any way, x,y,color,...
	Wave w										// wave to look for, checks is wave y or x or image or used to set color...
	Variable flag								// use: 1=graphs, 2=tables, 4=gizmos, or any sum of 1,2,4, 7 gives all
	flag = round(flag) & 7					// flag is a bit mask
	if (!WaveExists(w) || numtype(flag) || flag==0)
		return ""
	endif
	String out=""
	if (flag & 4)								// flag includes Gizmos
		out = FindGizmosWithWave(w)
	endif
	if ((flag&3) == 0)						// if no Graphs or Tables, we are done
		return out
	endif
	String win,wlist = WinList("*",";","WIN:"+num2istr(flag)), clist, cwin
	Variable i,m,Nm=ItemsInList(wlist)
	for (m=0;m<Nm;m+=1)
		win = StringFromList(m,wlist)
		CheckDisplayed/W=$win w
		if (V_flag>0)
			out += win+";"
		else
			clist = ChildWindowList(win)	// check for sub-windows
			for (i=0;i<ItemsInLIst(clist);i+=1)
				cwin = StringFromList(i,clist)
				CheckDisplayed/W=$(win+"#"+cwin) w
				if (V_flag>0)
					out += win+"#"+cwin+";"
					break
				endif
			endfor
		endif
	endfor
	return out
End
//
//		********** This Function is DEPRECATED, just call WindowsWithWave(w,1) directly **********
Function/S FindGraphsWithWave(w)	// find the graph window which contains the specified wave
	Wave w
	return WindowsWithWave(w,1)
End
//
//		********** This Function is DEPRECATED, just call WindowsWithWave(w,2) directly **********
Function/S FindTablesWithWave(w)	// find the table windows which contains the specified wave
	Wave w
	return WindowsWithWave(w,2)
End
//
//		********** Do not directly call this, use: WindowsWithWave(w,4) **********
Function/T FindGizmosWithWave(w)	// find list of Gizmos that contain the specified wave
	// the wave may be ANY wave that shows up in the objectList, this includes RGBA, line size...
	Wave w
	if (!WaveExists(w) || exists("NewGizmo")!=4)
		return ""
	endif

#if (IgorVersion()<7)
	Execute "GetGizmo/Z gizmoNameList"
	String gList=StrVarOrDefault("S_GizmoNames","")
	KillStrings/Z S_GizmoNames
#else
	GetGizmo/Z gizmoNameList
	String gList = SelectString(strlen(S_GizmoNames)>0, "", S_GizmoNames)
#endif
	String gName, objetAll, list=""
	Variable i, m, N=ItemsInList(gList)
	for (i=0;i<N;i+=1)
		gName=StringFromList(i,gList)				// check each gizmo that is displayed
#if (IgorVersion()<7)
		Execute "GetGizmo/Z/N="+gName+" objectList"
		KillStrings/Z S_gizmoObjectList
#else
		GetGizmo/Z/N=$gName objectList
#endif
		Wave/T TW_gizmoObjectList=TW_gizmoObjectList
		objetAll = ""
		for (m=0;m<numpnts(TW_gizmoObjectList);m+=1)
			objetAll += TW_gizmoObjectList[m]+";"
		endfor
		objetAll = ReplaceString("=",objetAll,";")
		objetAll = ReplaceString(",",objetAll,";")
		objetAll = ReplaceString("}",objetAll,";")
		objetAll = ReplaceString("{",objetAll,";")
		objetAll = ReplaceString("\r",objetAll,";")
		objetAll = ReplaceString("\n",objetAll,";")
		if (strsearch(objetAll,";"+GetWavesDataFolder(w,2)+";",0,2)>=0)
			list += gName+";"
		endif
	endfor
	KillWaves/Z TW_gizmoObjectList
	return list
End


Function isWaveOnGizmo(gName, ww)	// returns True when ww is in object list of Gizmo
	String gName			// name of gizmo, use "" for top gizmo
	Wave ww					// a wave in the object list of the gizmo
	if (!WaveExists(ww) || exists("NewGizmo")!=4)
		return 0
	endif

#if (IgorVersion()<7)
	if (strlen(gName))
		Execute "GetGizmo/Z/N="+gName+" objectList"
	else
		Execute "GetGizmo/Z objectList"
	endif
	KillStrings/Z S_gizmoObjectList
#else
	GetGizmo/Z/N=$gName objectList
#endif
	Wave/T TW_gizmoObjectList=TW_gizmoObjectList
	String objetAll="", list=""
	Variable m
	for (m=0;m<numpnts(TW_gizmoObjectList);m+=1)
		objetAll += TW_gizmoObjectList[m]+";"
	endfor
	KillWaves TW_gizmoObjectList

	objetAll = ReplaceString("=",objetAll,";")
	objetAll = ReplaceString(",",objetAll,";")
	objetAll = ReplaceString("}",objetAll,";")
	objetAll = ReplaceString("{",objetAll,";")
	objetAll = ReplaceString("\r",objetAll,";")
	objetAll = ReplaceString("\n",objetAll,";")
	return strsearch(objetAll,";"+GetWavesDataFolder(ww,2)+";",0,2)>=0
End
//
//Function/S FindGraphsWithWave(w)	// find the graph window which contains the specified wave
//	Wave w
//	if (!WaveExists(w))
//		return ""
//	endif
//	String name0=GetWavesDataFolder(w,2), out=""
//	String win,wlist = WinList("*",";","WIN:1"), clist, cwin
//	Variable i,m,Nm=ItemsInList(wlist)
//	for (m=0;m<Nm;m+=1)
//		win = StringFromList(m,wlist)
//		CheckDisplayed/W=$win w
//		if (V_flag>0)
//			out += win+";"
//		else
//			clist = ChildWindowList(win)
//			for (i=0;i<ItemsInLIst(clist);i+=1)
//				cwin = StringFromList(i,clist)
//				CheckDisplayed/W=$(win+"#"+cwin) w
//				if (V_flag>0)
//					out += win+"#"+cwin+";"
//					break
//				endif
//			endfor
//		endif
//	endfor
//	return out
//End
//
//Function/T FindTablesWithWave(w)	// find the table windows which contains the specified wave
//	Wave w
//	if (!WaveExists(w))
//		return ""
//	endif
//	String name0=GetWavesDataFolder(w,2), out=""
//	String win,wlist = WinList("*",";","WIN:2"), clist, cwin
//	Variable i,m,Nm=ItemsInList(wlist)
//	for (m=0;m<Nm;m+=1)
//		win = StringFromList(m,wlist)
//		CheckDisplayed/W=$win w
//		if (V_flag>0)
//			out += win+";"
//		else
//			clist = ChildWindowList(win)
//			for (i=0;i<ItemsInLIst(clist);i+=1)
//				cwin = StringFromList(i,clist)
//				CheckDisplayed/W=$(win+"#"+cwin) w
//				if (V_flag>0)
//					out += win+"#"+cwin+";"
//					break
//				endif
//			endfor
//		endif
//	endfor
//	return out
//End



Function/WAVE DisplayTableOfWave(ww,[classes,promptStr,names,options,colWid,top,left,maxRows,maxCols,align])
	// put up a table of a 2D wave
	Wave ww						// the ZonesWave wave (4 columns)
	String classes				// can use commas or semicolons and *
	String promptStr			// string for use in prompts
	String names				// name wildcards as used in WaveList()
	String options				// WaveList() options (e.g. "DIMS:2")
	Variable colWid			// fixed column width (points)
	Variable top,left			// top left corner of table
	Variable maxRows			// maximum number of rows to show
	Variable maxCols			// maximum number of columns to show
	Variable align				// alignment, 0=left, 1=center, 2=right (default=2)
	classes = SelectString(ParamIsDefault(classes),classes,"*")
	promptStr = SelectString(ParamIsDefault(promptStr),promptStr,"Select a Wave for Table")
	names = SelectString(ParamIsDefault(names),names,"*")
	options = SelectString(ParamIsDefault(options),options,"")
	colWid = ParamIsDefault(colWid) || colWid<=1 || numtype(colWid) ? 80 : colWid
	top = ParamIsDefault(top) || numtype(top) ? 44 : limit(top,44,500)
	left = ParamIsDefault(left) || numtype(left) ? 5 : limit(left,0,784)
	maxRows = ParamIsDefault(maxRows) || numtype(maxRows) ? 90 : limit(maxRows,3,200)
	maxCols = ParamIsDefault(maxCols) || numtype(maxCols) ? 15 : limit(maxCols,3,30)
	align = ParamIsDefault(align) || numtype(align) ? 2 : limit(round(align),0,2)

	if (!WaveExists(ww))
		String list=WaveListClass(classes,names,options)
		if(ItemsInList(list)==1)
			Wave ww = $StringFromList(0,list)
		elseif (ItemsInList(list)>1)
			String name
			Prompt name,promptStr+" Wave",popup,list
			DoPrompt promptStr,name
			if (!V_flag)
				Wave ww= $name
			endif
		endif
	endif
	if (!WaveExists(ww))
		return $""
	endif
	String win=StringFromList(0,WindowsWithWave(ww,2))	// find an existing table windows containing ww
	if (strlen(win))
		DoWindow/F $win								// Table already exists, just bring it to front
		return ww
	endif

	String fontName=GetDefaultFont("")
	Variable fontSize = 12

	Variable i, Nr=DimSize(ww,0), Nc=DimSize(ww,1), hasColLabels, hasRowLabels
	for (i=0,hasColLabels=0; i<Nc; i+=1)
		hasColLabels = hasColLabels || strlen(GetDimLabel(ww,1,i))
	endfor
	for (i=0,hasRowLabels=0; i<Nr; i+=1)
		hasRowLabels = max(hasRowLabels, FontSizeStringWidth(fontName,fontSize,0,GetDimLabel(ww,0,i)))
	endfor
	Nr = min(maxRows,Nr)							// limit display to maxRows rows
	Nc = min(maxCols,Nc)							// limit display to maxCols columns

	String screen1=StringByKey("SCREEN1",IgorInfo(0))
	i = strsearch(screen1,"RECT=",0)
	String rect = StringByKey("RECT",screen1[i,Inf],"=")
	Variable scrWidth=str2num(StringFromList(2,rect,",")), scrHeight=str2num(StringFromList(3,rect,","))
	scrWidth = numtype(scrWidth) || scrWidth<800 ? 1000 : scrWidth
	scrHeight = numtype(scrHeight) || scrHeight<600 ? 600 : scrHeight

	Variable width = 53 + colWid*Nc, lwidth=0
	Variable height = 86 + (FontSizeHeight(fontName,fontSize,0))*Nr
	if (hasColLabels || hasRowLabels)
		lwidth = max(hasRowLabels+3, 20)
		width += lwidth
		height +=  Nr*4
	endif
	width = min(scrWidth-2*left, width)
	height =  min(scrHeight-top-5, height)

	if (lwidth)											// Make the table /W=(left,top,right,bottom)
		Edit/W=(left,top,left+width,top+height)/K=1 ww.ld
		ModifyTable width(ww.l)=lwidth
	else
		Edit/W=(left,top,left+width,top+height)/K=1 ww
	endif
	ModifyTable font=fontName, size=fontSize, format(Point)=1,width(Point)=36, width(ww.d)=colWid, alignment=align
	if (WaveType(ww) & 0x38)				// if ww is an integer (8, 16, or 32 bit) use integer format
		ModifyTable format(ww.d)=1
	endif
	return ww
End



// make a function getFiletype that goes to a file and determines it's filetype, this applies to old $tag type files
Function/T getListOfTypesInFile(fname,path)
	String fname									// full path name to file with tagged geometry values
	String path									// name of Igor path

	Variable refNum
	Open/M="file containing tagged values"/P=$path/R/Z=2 refNum as fname
	if (strlen(S_fileName)<1 || !refNum)
		return ""
	endif
	String line
	FReadLine refNum, line
	Close refNum
	if (!stringmatch(line[0],"$"))				// if it does not start with a '$' then forget it.
		return ""
	endif

	Variable i
	i = strlen(line)
	i -= char2num(line[i-1])<32 ? 1 : 0
	i -= char2num(line[i-2])<32 ? 1 : 0
	line = line[0,i-1]							// trim off any c/r or new lines
	i = strsearch(line,"//",0)					// remove comment
	if (i>0)
		line = line[0,i-1]
	endif
	line = ReplaceString("\t",line," ")			// change all tabs to spaces
	line += " "										// ensure a space terminator

	if (strsearch(line,"$filetype ",0)!=0)		// does not starts with $filetype, pass back the tag
		i = strsearch(line," ",0)
		return line[1,i-1]
	endif
	line = line[9,Inf]								// strip off the $filetype tag

	// now get the value
	for (i=0;i<strlen(line) && char2num(line[i])<=32;i+=1)
	endfor
	line = line[i,Inf]								// trim off leading spaces

	for (i=strlen(line)-1;i>=0 && char2num(line[i])<=32;i-=1)
	endfor
	line = line[0,i]

	line = ReplaceString(",",line,";")			// change all commas, and spaces to semi-colons
	line = ReplaceString(" ",line,";")
	return line
End


Function DrawMarker(x0,y0,dx,dy,style,[color,thick,dash,win,layer,grpName])
	Variable x0,y0,dx,dy
	String style
	String color						// one of "red", "blue", "green", "yellow", "magenta", "cyan", or a triplet like "1,1,0" or "65535,0,0"
	Variable thick						// line thickness, defaults to 0.50 (thin)
	Variable dash
	String win							// optional name of window
	String layer						// Drawing layer to use.  Default is "UserFront"
	String grpName						// name of group of draw commands
	if (numtype(x0+y0+dx+dy) || dx<=0 || dy<=0)
		return 1
	endif
	if (ParamIsDefault(color))
		color = SelectString(stringmatch(style,"BoxWithTicks"),"black","30583,30583,30583")	// special default color for BoxWithTicks
	endif
//	color = SelectString(ParamIsDefault(color),color,"black")	// default color to black
//	color = SelectString(ParamIsDefault(color) && stringmatch(style,"BoxWithTicks"),color,"30583,30583,30583")	// special default color for BoxWithTicks
	thick = ParamIsDefault(thick) ? 0.5 : thick
	dash = ParamIsDefault(dash) ? 0 : dash
	if (ParamIsDefault(win))
		win = ""
	endif
	if (ParamIsDefault(layer))
		layer = "UserFront"
	endif
	String layerList = "ProgBack;UserBack;ProgAxes;UserAxes;ProgFront;UserFront"
	layer = SelectString(strlen(layer),"UserFront",layer)
	if (WhichListItem(layer,layerList)<0)
		printf "ERROR -- layer = '%s' is invlaid, must be one of '%s'\r",layer,layerList
		return 1
	endif

	String rgb = StringByKey(color,"red:1,0,0;blue:0,0,1;green:0,1,0;yellow:1,1,0;magenta:1,0,1;cyan:0,1,1;black:0,0,0;white:1,1,1")	// name -> numbers
	rgb = SelectString(strlen(rgb),color,rgb)
	Variable r=str2num(StringFromList(0,rgb,",")), g=str2num(StringFromList(1,rgb,",")), b=str2num(StringFromList(2,rgb,","))
	r = !(r>=0) ? 0 : r											// default invalid to black
	g = !(g>=0) ? 0 : g
	b = !(b>=0) ? 0 : b
	if (r+g+b<=3)												// if numbers are all small change from [0,1] -> [0,65535]
		r *= 65535
		g *= 65535
		b *= 65535
	endif

	style = SelectString(strlen(style),"cross",style)			// defaults to Cross
	String list = ImageInfo(win,"",0)
	if (strlen(list)<1)											// perhaps no image, only a graph
		list = TraceInfo(win,StringFromList(0,TraceNameList(win,";",1+4)),0)
	endif
	String xaxis=StringByKey("XAXIS",list), yaxis=StringByKey("YAXIS",list)
	String aList = AxisList(win)
	aList = RemoveFromList(yaxis,aList)
	if (strlen(xaxis)<1)										// I am getting desperate, last thing to try
		Variable i
		for (i=0;i<ItemsInList(aList);i+=1)
			xaxis = StringFromList(i,aList)
			if (strsearch(xaxis,"bottom",0,2)>=0)
				break
			elseif (strsearch(xaxis,"top",0,2)>=0)
				break
			endif
		endfor
	endif
	if (strlen(yaxis)<1)										// I am getting desperate, last thing to try
		aList = RemoveFromList(xaxis,aList)
		for (i=0;i<ItemsInList(aList);i+=1)
			yaxis = StringFromList(i,aList)
			if (strsearch(yaxis,"left",0,2)>=0)
				break
			elseif (strsearch(yaxis,"right",0,2)>=0)
				break
			endif
		endfor
	endif

	dy = abs(dy/2)		// change FW to HW
	dx = abs(dx/2)

	if (strlen(grpName))
		grpName = CleanupName(grpName,0)
		grpName = ReplaceString("__",grpName,"_")
		if (char2num(grpName)==char2num("_"))
			grpName = "M"+grpName
		endif
	else
		init_UtiltyJZT()
		NVAR markerIndex = root:Packages:UtilityJZT:markerIndex
		grpName = "M"+num2istr(markerIndex)
		markerIndex += 1
	endif

	SetDrawLayer/W=$win $layer
	SetDrawEnv/W=$win gstart, gname=$grpName

	Variable err = 0
	if (stringmatch(style,"cross gap"))
		SetDrawEnv/W=$win  xcoord= $xaxis,ycoord= $yaxis,linefgc=(r,g,b), linethick= thick, dash=dash
		DrawLine/W=$win  x0,(y0-dy),x0,(y0-0.1*dy)
		SetDrawEnv/W=$win  xcoord= $xaxis,ycoord= $yaxis,linefgc=(r,g,b), linethick= thick, dash=dash
		DrawLine/W=$win  x0,(y0+0.1*dy),x0,(y0+dy)
		SetDrawEnv/W=$win xcoord= $xaxis,ycoord= $yaxis,linefgc=(r,g,b), linethick= thick, dash=dash
		DrawLine/W=$win (x0-dx),y0,(x0-0.1*dx),y0
		SetDrawEnv/W=$win xcoord= $xaxis,ycoord= $yaxis,linefgc=(r,g,b), linethick= thick, dash=dash
		DrawLine/W=$win (x0+0.1*dx),y0,(x0+dx),y0
	elseif (stringmatch(style,"cross"))
		SetDrawEnv/W=$win  xcoord= $xaxis,ycoord= $yaxis,linefgc=(r,g,b), linethick= thick, dash=dash
		DrawLine/W=$win  (x0),(y0-dy),x0,(y0+dy)
		SetDrawEnv/W=$win xcoord= $xaxis,ycoord= $yaxis,linefgc=(r,g,b), linethick= thick, dash=dash
		DrawLine/W=$win (x0-dx),y0,(x0+dx),y0
	elseif (stringmatch(style,"X"))
		SetDrawEnv/W=$win xcoord= $xaxis,ycoord= $yaxis,linefgc=(r,g,b), linethick= thick, dash=dash
		DrawLine/W=$win (x0-dx),(y0-dy),(x0+dx),(y0+dy)
		SetDrawEnv/W=$win xcoord= $xaxis,ycoord= $yaxis,linefgc=(r,g,b), linethick= thick, dash=dash
		DrawLine/W=$win (x0+dx),(y0-dy),(x0-dx),(y0+dy)
	elseif (stringmatch(style,"BoxWithTicks"))
		Variable hw							// half width of box in pixels, usually ~10
		SetDrawEnv/W=$win xcoord= $xaxis,ycoord= $yaxis,fillpat= 0, linethick= thick, dash=dash
		DrawRect x0-dx,y0-dy,x0+dx,y0+dy
		SetDrawEnv/W=$win xcoord= $xaxis,ycoord= $yaxis,linefgc=(r,g,b), linethick= thick, dash=dash
		DrawLine/W=$win x0,y0-dy,x0,y0-dy/2
		SetDrawEnv/W=$win xcoord= $xaxis,ycoord= $yaxis,linefgc=(r,g,b), linethick= thick, dash=dash
		DrawLine/W=$win x0,y0+dy/2,x0,y0+dy
		SetDrawEnv/W=$win xcoord= $xaxis,ycoord= $yaxis,linefgc=(r,g,b), linethick= thick, dash=dash
		DrawLine/W=$win x0-dx,y0,x0-dx/2,y0
		SetDrawEnv/W=$win xcoord= $xaxis,ycoord= $yaxis,linefgc=(r,g,b), linethick= thick, dash=dash
		DrawLine/W=$win x0+dx/2,y0,x0+dx,y0
	else
		err = 1
	endif
	SetDrawEnv/W=$win gstop
	SetDrawLayer/W=$win UserFront
	return err
End



// make saturated color RGB from point (dx,dy) on circle of radius rmax. RGB scaled to rgbMax
ThreadSafe Function xy2saturatedColors(dx,dy,rmax,rgbMax,rgb)
	Variable dx,dy												// x,y location on the pole figure
	Variable rmax												// max radius (default is 1)
	Variable rgbMax											// scale for RGB values, usually 1 or 65535
	Wave rgb													// set rgb to scale [0,rgbMax]
	if (numtype(dx+dy+rgbMax+rmax) || !(rgbMax>0) || !(rmax>0))
		rgb = NaN
		return NaN
	endif

	Variable rad = min(1,sqrt(dx*dx + dy*dy)/rmax)		// r/rmax, but does not exceed 1
	Variable hue = mod(atan2(dy,dx)*180/PI+360,360)	// angle in range [0,360)
	Variable huei = mod(floor(hue/60),6)					// {0, 1, 2, 3, 4, 5, 6}
	Variable pp,qq,tt, f
	f = hue/60 - floor(hue/60)
	pp = (1-rad)*rgbMax
	qq = (1-f*rad)*rgbMax
	tt = (1-(1-f)*rad)*rgbMax
	Variable r, g, b
	if (huei==0)
		r = rgbMax
		g = tt
		b = pp
	elseif(huei==1)
		r = qq
		g = rgbMax
		b = pp
	elseif(huei==2)
		r = pp
		g = rgbMax
		b = tt
	elseif(huei==3)
		r = pp
		g = qq
		b = rgbMax
	elseif(huei==4)
		r = tt
		g = pp
		b = rgbMax
	elseif(huei==5)
		r = rgbMax
		g = pp
		b = qq
	endif
	rgb = {r,g,b}
	rgb = limit(rgb,0,rgbMax)
	return hue						// return hue angle (degree)
End
//Function test_xy2saturatedColors()
//	Variable angle, dx, dy, rgbMax=65535, rmax=1
//	Make/N=3/D/FREE rgbi
//	Make/N=(361,3)/O rgbAll
//	SetScale/P x 0,1,"", rgbAll
//	for (angle=0;angle<361;angle+= 1)
//		dx=cos(angle*PI/180)
//		dy=sin(angle*PI/180)
//		xy2saturatedColors(dx,dy,rmax,rgbMax,rgbi)
//		rgbAll[angle][] = rgbi[q]
//	endfor
//	rgbAll = abs(rgbAll) < 1e-10 ? 0 : rgbAll
//End



ThreadSafe Function/T reverseList(in,[sep])	// reverse the order of items in a list
	String in										// original list
	String sep										// separator between items
	if (ParamIsDefault(sep))
		sep = ";"
	endif
	String out=""
	Variable i, N=ItemsInList(in,sep)

	for (i=0;i<N;i+=1)
		out = AddListItem(StringFromList(i,in,sep),out,sep,0)
	endfor
	return out
End



#if (IgorVersion()<7)
//	For Igor7 use:
//		FindDuplicates/RT=typeNoDups type
//		Duplicate/FREE/T typeNoDups, atomTypes
//		KillWaves/Z typeNoDups
//	For Igor6 use:
//		Wave/T atomTypes = UniqueWaveValuesT(type)
//
//	I have made these static since it is only for Igor6, and it is time to start using Igor7 (or 8)
//
Static Function/WAVE UniqueWaveValues(wwIN, [tol])		// only for Igor 6
	Wave wwIN
	Variable tol
	tol = ParamIsDefault(tol) || numtype(tol) || tol<=0 ? 0: tol
	if (!WaveExists(wwIN) || numtype(tol) || tol<0)
		return $""
	endif
	Duplicate/FREE wwIN, ww, unique
	Variable i, m, val, N=numpnts(ww)

	for (i=0,m=0; i<N; i+=1)
		val = ww[i]
		if (numtype(val)==2)
			continue
		endif		
		MatrixOP/FREE flags = greater(Abs(ww - val), tol)
		flags = !flags
		ww = flags[p] ? NaN : ww[p]
		unique[m] = val
		m += 1
	endfor
	Redimension/N=(m) unique
	return unique
End
//
Static Function/WAVE UniqueWaveValuesT(wwIN)			// only for Igor 6
	Wave/T wwIN
	if (!WaveExists(wwIN))
		return $""
	endif
	Duplicate/T/FREE wwIN, ww, unique
	Variable i, m, N=numpnts(ww)
	Make/N=(N)/B/FREE flags

	String str, invalid="INVALID_JZT_897"
	for (i=0,m=0; i<N; i+=1)
		str = ww[i]
		if (cmpstr(str,invalid)==0)
			continue
		endif
		flags[] = cmpstr(str,ww[p],1)
		ww = SelectString(flags[p],invalid,ww[p])
		unique[m] = str
		m += 1
	endfor
	Redimension/N=(m) unique
	return unique
End
#endif



ThreadSafe Function monotonic(a,[printIt])
	// determines whether values of a particular wave are monotonic increasing (if any NaN, returns false)
	Wave a		// the wave
	Variable printIt
	printIt = ParamIsDefault(printIt) ? 0 : printIt
	printIt = numtype(printIt) ? 0 : !(!printIt)

	if (!WaveExists(a))
		if (printIt)
			print "ERROR -- monotonic(), wave not found"
		endif
		return NaN
	elseif (numpnts(a)<=0)
		if (printIt)
			print "ERROR -- monotonic(), wave is empty"
		endif
		return NaN
	elseif (WaveDims(a)>1)
		if (printIt)
			print "ERROR -- monotonic(), wave is not 1D"
		endif
		return NaN
	endif
	MatrixOP/FREE delta = a - rotateRows(a,1)
	delta[0] = 0
	WaveStats/M=1/Q delta
	return !(V_min<0 || V_numNans)
End


ThreadSafe Function isWaveConstant(wav,[places,checkCase])
	// returns 1 if all elements of wav are the same
	Wave wav
	Variable places								// used to compare floating point waves
	Variable checkCase							// only used for text waves, defaults to 1

	Variable N=numpnts(wav), i,j, x0

	if (WaveType(wav,1)==0)					// bad input wave, null wave
		return 0
	elseif (N<=1)									// empty or only 1 item, must be constant
		return 1

	elseif (WaveType(wav,1)==2)										// text wave
		checkCase = ParamIsDefault(checkCase) || numtype(checkCase) ? 1 : !(!checkCase)
		Wave/T wt = wav
		String str=wt[0]
		for (i=1; i<N; i+=1)
			if (cmpstr(wt[i],str,checkCase))
				return 0
			endif
		endfor
		return 1

	elseif (WaveType(wav,1)==1 && WaveType(wav,0) & 0x78)	// integer wave
		x0 = wav[0]
		MatrixOP/FREE delta = maxVal( Abs(wav - x0) )
		return (delta[0] < 0.05)

	elseif (WaveType(wav,1)==1 && WaveType(wav,0) & 0x07)	// floating point wave
		places = ParamIsDefault(places) || numtype(places) || places<=0 ? 0 : !(!places)
		if (places==0)
			places = WaveType(wav,0) & 0x02 ? 5 : 15			// 32bit-->5 places,  64bit-->15 places
		endif
		x0 = max(WaveMax(wav),-WaveMin(wav))
		MatrixOP/FREE delta = maxVal( Abs(wav - x0) )
		if (delta[0]==0)
			return 1
		endif
		return (delta[0]/abs(x0)) < (10^(-places))

	elseif (WaveType(wav,1)==3)										// wav holds data folder references
		Wave/DF wd = wav
		DFREF DFR = wd[0]
		for (i=1;i<N;i+=1)
			if (!DataFolderRefsEqual(wd[i],DFR))
				return 0
			endif
		endfor
		return 1

	elseif (WaveType(wav,1)==3)										// wav holds wave references
		Wave/WAVE ww = wav
		Wave/WAVE w0 = ww[0]
		for (i=1;i<N;i+=1)
			if (!WaveRefsEqual(ww[i],w0))
				return 0
			endif
		endfor
		return 1
	endif

	return 0					// wav is not known
End


ThreadSafe Function isdigit(c)				// returns 1 if c is a digit, 0-9 (otherwise returns 0)
	String c
	Variable i=char2num(c)
	return (48<=i && i<=57)
End


ThreadSafe Function isletter(c)				// returns 1 if c is an upper or lower case letter (otherwise returns 0)
	String c
	Variable i=char2num(c)
	return (65<=i && i<=90) || (97<=i && i<=122)
End


ThreadSafe Function countOccurancesOfStr(main,find)	// count number of times find occurs in main
	String main				// main string
	String find				// count occurances of find in main
	Variable n=0, i0=0, len=strlen(find)
	do
		i0 = strsearch(main,find,i0+len,2)
		n += 1
	while (i0>=0)
	return (n-1)
End


ThreadSafe Function angleVec2Vec(a,b)		// return the angle between two vectors (degree)
	Wave a,b
	Variable dot = MatrixDot(a,b) / (norm(a)*norm(b))
	dot = limit(dot,-1,1)						// ensure that the acos will exist
	return acos(dot)*180/PI
End


ThreadSafe Function/WAVE perp2Vector(a)
	// returns a new vector that is perpendicular to a, the direction around a is unspecified.
	Wave a
	Variable i,j, imin=abs(a[0])<abs(a[1]) ? 0 : 1
	imin = abs(a[2])<abs(a[imin]) ? 2 : imin
	i = mod(imin+2,3)
	j = mod(imin+1,3)
	Make/N=(3,3)/D/FREE mat=0
	mat[i][j] = 1						//	see:  https://en.wikipedia.org/wiki/Cross_product
	mat[j][i] = -1
	MatrixOP/FREE cr = Normalize(mat x a)
	return cr
End


ThreadSafe Function rotationAngleOfMat(rot)	// returns the total rotation angle of a matrix 'rot'
	Wave rot											// the rotation matrix
	Variable trace = MatrixTrace(rot)		// trace = 1 + 2*cos(theta)
	Variable cosine = (trace-1)/2			// cosine of the rotation angle
	cosine = (cosine>1) ? (2-cosine) : cosine
	return acos(cosine)*180/PI				// rotation angle in degrees
End


ThreadSafe Function isRotationMat(mat,[tol])	// true if mat is a rotation matrix
	Wave mat
	Variable tol							// positive threshold for zero
	if (!WaveExists(mat))
		return 0								// mat does not exist
	elseif (DimSize(mat,0)<2 || (DimSize(mat,0) != DimSize(mat,1)))
		return 0								// mat not a square matrix, at least 2x2
	endif 
	tol = ParamIsDefault(tol) || numtype(tol) || tol<=0 ? NaN : tol
	if (numtype(tol))
		tol = WaveType(mat) & 0x04 ? 1e-13 : 1e-6
	endif
	MatrixOP/FREE err = sum(abs((mat x (mat^t)) - Identity(3)))
	return (err[0]<tol)
End


// compute angle and axis of a rotation matrix
// Aug 2, 2007, this was giving the wrong sign for the rotation, so I reversed the "curl" in defn of axis.  JZT
//		changed "axis[0] = rot[1][2] - rot[2][1]"   -->   "axis[0] = rot[2][1] - rot[1][2]"
// Oct 12, 2017, changed to use a Quaternion based method:
//			http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/
//	*** NOTE: for DeviatoricStrainRefine(), angle MUST be POSITIVE
ThreadSafe Function axisOfMatrix(matIN,axis,[squareUp])
	// returns total rotation angle (deg), and sets axis to the axis of the total rotation
	Wave matIN								// should be a rotation matrix
	Wave axis								// axis of the rotation (angle is returned)
	Variable squareUp						// optionally square up mat (default is NOT square up)
	squareUp = ParamIsDefault(squareUp) ? NaN : squareUp
	squareUp = numtype(squareUp) ? 0 : !(!squareUp)

	Duplicate/FREE matIN, mat
	Variable det = MatrixDet(mat)
	if (abs(det)<1e-12)					// a zero matrix, no rotation
		axis = 0
		return 0
	elseif (det<0)							// contains an inversion
		mat *= -1							// ensure that |mat| is always positive
	endif

	if (squareUp)
		if (SquareUpMatrix(mat))
			axis = NaN						// default for error
			return NaN
		endif
	endif

	Variable qw, SS, tr=mat[0][0]+mat[1][1]+mat[2][2]
	if (tr>0)
		SS = 2*sqrt(tr+1)							// S=4*qw 
		qw = SS/4
		axis[0] = mat[2][1] - mat[1][2]		// any common factors are irrelevant since I next normalize axis
		axis[1] = mat[0][2] - mat[2][0]
		axis[2] = mat[1][0] - mat[0][1]
		axis /= SS
	elseif ((mat[0][0] > mat[1][1]) & (mat[0][0] > mat[2][2]))
		SS = 2*sqrt(1 + mat[0][0] - mat[1][1] - mat[2][2])		// S=4*qx
		qw = (mat[2][1] - mat[1][2]) / SS
		axis[0] = 0.25 * SS
		axis[1] = (mat[0][1] + mat[1][0]) / SS 
		axis[2] = (mat[0][2] + mat[2][0]) / SS 
	elseif (mat[1][1] > mat[2][2])
		SS = 2*sqrt(1 + mat[1][1] - mat[0][0] - mat[2][2])		// S=4*qy
		qw = (mat[0][2] - mat[2][0]) / SS
		axis[0] = (mat[0][1] + mat[1][0]) / SS 
		axis[1] = SS / 4
		axis[2] = (mat[1][2] + mat[2][1]) / SS 
	else
		SS = 2*sqrt(1.0 + mat[2][2] - mat[0][0] - mat[1][1])	// S=4*qz
		qw = (mat[1][0] - mat[0][1]) / SS
		axis[0] = (mat[0][2] + mat[2][0]) / SS
		axis[1] = (mat[1][2] + mat[2][1]) / SS
		axis[2] = SS / 4
	endif
	normalize(axis)

	Variable angle = 2*acos(limit(qw,-1,1))
	if (sum(axis)<0)
		axis *= -1
		angle *= -1
	endif

	angle += angle<0 ? (2*PI) : 0			// always return angle in range [0, 360)
	return mod(angle*180/PI, 360)			// rotation angle in degrees (always a positive number)
End
//ThreadSafe Function axisOfMatrix(mat,axis,[squareUp])
//	// returns total rotation angle (deg), and sets axis to the axis of the total rotation
//	Wave mat									// should be a rotation matrix
//	Wave axis								// axis of the rotation (angle is returned)
//	Variable squareUp						// optionally square up mat (default is NOT square up)
//	squareUp = ParamIsDefault(squareUp) ? NaN : squareUp
//	squareUp = numtype(squareUp) ? 0 : !(!squareUp)
//
//	Make/N=(3,3)/FREE/D rot=mat
//	if (squareUp)
//		if (SquareUpMatrix(rot))
//			axis = NaN								// default for error
//			return NaN
//		endif
//	else
//		MatrixOp/FREE sumd = sum(Abs((mat x mat^t) - Identity(3)))
//		if (sumd[0]<0 || sumd[0]>1e-4)		// not close enough to a rotation mat, an error
//			axis = NaN								// default for error
//			return NaN
//		endif
//	endif
//
//	Variable cosine = (MatrixTrace(rot)-1)/2	// trace = 1 + 2*cos(theta)
//	cosine = limit(cosine,-1,1)
//	if (cosine<= -1)								// special for 180 degree rotation,
//		axis[0] = sqrt((rot[0][0]+1)/2)
//		axis[1] = sqrt((rot[1][1]+1)/2)
//		axis[2] = sqrt((rot[2][2]+1)/2)		// always assume z positive
//		axis[0] = (rot[0][2]+rot[2][0])<0 ? -axis[0] : axis[0]
//		axis[1] = (rot[1][2]+rot[2][1])<0 ? -axis[1] : axis[1]
//		if (numtype(sum(axis)))				// this is for very special cases such as diag = {1,-1,-1}
//			WaveStats/M=1/Q axis
//			axis = p==V_maxloc
//		endif
//	else												// rotaion < 180 degree, usual formula works
//		axis[0] = rot[2][1] - rot[1][2]
//		axis[1] = rot[0][2] - rot[2][0]
//		axis[2] = rot[1][0] - rot[0][1]
//		axis /= 2
//	endif
//	normalize(axis)
//	return acos(cosine)*180/PI				// rotation angle in degrees
//End


// Oct 5, 2014, used new proper method for changing rot to an exact rotation matrix
ThreadSafe Function SquareUpMatrix(rot)
	// see http://en.wikipedia.org/wiki/Kabsch_algorithm
	Wave rot
	MatrixSVD/U=0/V=0/Z rot
	if (V_flag)
		return 1
	endif
	Wave M_V=M_U			// 3x3 column-orthonormal matrix
	Wave M_WT=M_VT			// transpose of NxN orthonormal matrix

	// A = M_V x W_S x M_WT		// from wikipedia
	MatrixOp/FREE deter = Det(M_WT^t x M_V^t)
	Make/N=(3,3)/FREE/D Idet=(p==q)
	Idet[2][2] = sign(deter[0])

	MatrixOp/FREE rr = (M_WT^t x Idet x M_V^t)^t
	rot = rr
	KillWaves/Z W_W, M_U, M_VT
	return 0
End


// set mat to be a rotation matrix about axis with angle
ThreadSafe Function/WAVE rotationMatFromAxis(axis,angle)
	Wave axis			// axis about which to rotate (or possibly rotation vector in radians)
	Variable angle		// angle to rotate (degrees), assumes axis is true Rotation vector if angle invalid

	if (!WaveExists(axis) || numpnts(axis)<3)
		return $""
	endif
	Variable len = norm(axis)
	angle = numtype(angle) ? len : angle*PI/180	// the rotation angle (rad)
	if (numtype(angle))
		return $""
	endif

	Make/N=(3,3)/D/FREE rot								// desired rotation matrix
	if (angle==0)												// zero angle rotation is just the identity matrix
		rot = (p==q)
		return rot
	endif

	Variable nx=axis[0]/len, ny=axis[1]/len, nz=axis[2]/len
	Variable cosa=cos(angle), sina=sin(angle)
	Variable c1 = 1-cosa
	// from		http://mathworld.wolfram.com/RodriguesRotationFormula.html (I double checked this too.)
	rot[0][0] = cosa+nx*nx*c1;			rot[0][1] =nx*ny*c1-nz*sina;			rot[0][2] = nx*nz*c1+ny*sina;
	rot[1][0] = nx*ny*c1+nz*sina;		rot[1][1] = cosa+ny*ny*c1;			rot[1][2] =ny*nz*c1-nx*sina;
	rot[2][0] = nx*nz*c1-ny*sina;		rot[2][1] = ny*nz*c1+nx*sina;		rot[2][2] = cosa+nz*nz*c1;
	return rot
End


ThreadSafe Function smallestNonZeroValue(vec,[tol])
	// returns the abs( smallest non-zero element ), this routine ignores NaN's
	// {0, 0.1, 3} returns 0.1,  {0, -0.1, 3} returns 0.1
	//  This is useful if you want to display {0, 0.1, 3} as integers, then dividing gives {0, 1, 30}
	Wave vec								// a vector or matrix of any dimension
	Variable tol						// positive threshold for zero

	if (!WaveExists(vec))
		return NaN						// vec does not exist
	elseif (numpnts(vec)<1)
		return NaN
	endif
	tol = ParamIsDefault(tol) || numtype(tol) || tol<=0 ? NaN : tol
	if (numtype(tol))
		tol = WaveType(vec) & 0x04 ? 1e-13 : 1e-6
	endif

	if (WaveType(vec) & 0x80)		// special, MatrixOP does not handle 64 bit ints
		Duplicate/FREE, vec, temp
		Redimension/D temp				// MatrixOP does handle double precision
		MatrixOP/FREE temp = abs(temp)
	else
		MatrixOP/FREE temp = abs(vec)
	endif
	temp = temp < tol ? NaN : temp
	return WaveMin(temp)
End

Function/C waveRange(ww,[name,printIt])
	// returns the min and max value in wave, only works for numeric waves, ignores number of dimensions
	Wave ww
	String name							// when printIt is true, use this as the name
	Variable printIt
	printIt = ParamIsDefault(printIt) || numtype(printIt) ? strlen(GetRTStackInfo(2))==0 : printIt

	if (WaveType(ww,1) != 1)
		if (printIt)
			if (WaveType(ww,1)<1)
				printf "The given wave is NULL.\r"
			else
				printf "\"%s\"   is invalid, it is NOT numeric.\r", NameOfWave(ww)
			endif
		endif
		return cmplx(NaN,NaN)
	endif
	WaveStats/M=1/Q ww
	if (printIt)
		name = SelectString(ParamIsDefault(name), name, NameOfWave(ww))
		Variable variation = (V_max-V_min)==0 ? 0 : (V_max-V_min)/V_avg
		printf "\"%s\"   range=[%g, %g], %s=%.3g  --> %.2g%% variation\r", name,V_min,V_max, gDELTA, V_max-V_min, variation*100
	endif
	return cmplx(V_min,V_max)
End


ThreadSafe Function MedianOfWave(wwIN,f,[x1,x2,p1,p2])
	//	return median (or other percentile) of a wave, useful for picking min/max range for a color scale
	Wave wwIN
	Variable f					// fraction (e.g. percentile) in range [0,1], use f=0.5 for usual median
	Variable x1, x2			// range of interest (x-values) (cannot specify both x1 & p1)
	Variable p1, p2			// range of interest (point-values)
	if (!WaveExists(wwIN) || WaveDims(wwIN)!=1 || numtype(f) || f<0 || f>1)
		return NaN
	endif
	Variable setLo = (ParamIsDefault(x1) ? 0 : 2) + (ParamIsDefault(p1) ? 0 : 1)
	Variable setHi = (ParamIsDefault(x2) ? 0 : 2) + (ParamIsDefault(p2) ? 0 : 1)
	if (setLo>2 || setHi>2)
		return NaN								// cannot specify both x and point for a range
	endif

	Variable N=numpnts(wwIN), x0=DimOffset(wwIN,0), dx=DimDelta(wwIN,0)
	Variable pLo=0, pHi=N-1				// init pLo,pHi to full range
	pLo = setLo==1 ? p1 : pLo				// p1 given
	pHi = setHi==1 ? p2 : pHi				// p2 given
	pLo = setLo==2 ? (x1-x0)/dx : pLo	// x1 given, convert to p1
	pHi = setHi==2 ? (x2-x0)/dx : pHi	// x2 given, convert to p2
	if (numtype(pLo+pHi)==2)				// cannot deal with NaN's, probably from bad input
		return NaN
	endif
	pLo = limit(round(pLo),0,N-1)		// trim to allowed point range, and force to integer
	pHi = limit(round(pHi),0,N-1)
	Variable Nrange = (pHi-pLo+1)		// number of points in the sub-range
	if (Nrange < 1)							// nothing in range
		return NaN
	elseif (Nrange==1)
		return wwIN[0]							// trivial result
	endif

	Duplicate/R=[pLo,pHi]/FREE wwIN, temp
	Sort temp, temp							// Sort range of wwIN
	WaveStats/M=1/Q temp
	if (V_npnts<1)
		return NaN
	endif
	return temp[f*(V_npnts-1)]			// note that V_npnts*f need not be an integer
End
//
//Function test_MedianOfWave(bitFlag)
//	Variable bitFlag
//
//	Variable p1, p2
//	printf "\t0.5\t\t0\t\t\t1\t\t\t w\r"
//
//	if (bitFlag & 1)
//		Make/FREE tw={7}
//		printf "%4.2f,   %4.2f,   %4.2f,\t\t%s\r",MedianOfWave(tw,0.5),MedianOfWave(tw,0),MedianOfWave(tw,1),vec2str(tw)
//	endif
//	if (bitFlag & 2)
//		Make/FREE tw={4.1,2.1}
//		printf "%4.2f,   %4.2f,   %4.2f,\t\t%s\r",MedianOfWave(tw,0.5),MedianOfWave(tw,0),MedianOfWave(tw,1),vec2str(tw)
//	endif
//	if (bitFlag & 4)
//		Make/FREE tw={4.1,2.1,4}
//		printf "%4.2f,   %4.2f,   %4.2f,\t\t%s\r",MedianOfWave(tw,0.5),MedianOfWave(tw,0),MedianOfWave(tw,1),vec2str(tw)
//	endif
//
//	Make/FREE tw={1,2,3,4,5,0,6,7,8,9}
//	if (bitFlag & 8)
//		printf "%4.2f,   %4.2f,   %4.2f,\t\t%s\r",MedianOfWave(tw,0.5),MedianOfWave(tw,0),MedianOfWave(tw,1),vec2str(tw)
//	endif
//	if (bitFlag & 16)
//		p1=1; p2=10
//		Make/N=(min(p2,numpnts(tw)-1)-p1+1)/FREE tw12
//		tw12 = tw[p+p1]
//		printf "%4.2f,   %4.2f,   %4.2f,\t\t%s, p1=%d, p2=%d\r",MedianOfWave(tw,0.5,p1=p1,p2=p2),MedianOfWave(tw,0,p1=p1,p2=p2),MedianOfWave(tw,1,p1=p1,p2=p2),vec2str(tw12), p1,p2
//	endif
//	if (bitFlag & 32)
//		p1=2; p2=6
//		Make/N=(min(p2,numpnts(tw)-1)-p1+1)/FREE tw12
//		tw12 = tw[p+p1]
//		printf "%4.2f,   %4.2f,   %4.2f,\t\t%s, p1=%d, p2=%d\r",MedianOfWave(tw,0.5,p1=p1,p2=p2),MedianOfWave(tw,0,p1=p1,p2=p2),MedianOfWave(tw,1,p1=p1,p2=p2),vec2str(tw12), p1,p2
//	endif
//
//	if (bitFlag & 64)
//		p1=12; p2=6
//		printf "%4.2f,   %4.2f,   %4.2f,\t\t%s, p1=%d, p2=%d\r",MedianOfWave(tw,0.5,p1=p1,p2=p2),MedianOfWave(tw,0,p1=p1,p2=p2),MedianOfWave(tw,1,p1=p1,p2=p2),vec2str(tw), p1,p2
//	endif
//End


ThreadSafe Function/C rangeOfVec(wav,[row,col,layer,chunk])		// returns the max and min number in the specified vector
	Wave wav						// a 1D array (if more than 1D, MUST specify appropriate row, col, layer, chunk)
	Variable row,col,layer,chunk	// for selecting out a vector at a specific row, col, or layer, chunk
	if (WaveType(wav,1)!=1)				// only understand numeric
		return cmplx(NaN,NaN)
	endif

	Variable rowGiven=!ParamIsDefault(row), colGiven=!ParamIsDefault(col)
	Variable layerGiven=!ParamIsDefault(layer), chunkGiven=!ParamIsDefault(chunk)
	Variable dims = WaveDims(wav), type=WaveType(wav), N

	Wave vec = $""
	if (dims==1)
		Wave vec = wav
	elseif (dims==2)
		if (colGiven)							// want the row
			N = DimSize(wav,0)
			Make/N=(N)/Y=(type)/FREE vec
			vec = wav[p][col]
		endif
		if (rowGiven)							// want the column
			N = DimSize(wav,1)
			Make/N=(N)/Y=(type)/FREE vec
			vec = wav[row][p]
		endif
	elseif (dims==3)
		if (colGiven && layerGiven)			// want the row
			N = DimSize(wav,0)
			Make/N=(N)/Y=(type)/FREE vec
			vec = wav[p][col][layer]
		elseif (rowGiven && layerGiven)		// want the column
			N = DimSize(wav,1)
			Make/N=(N)/Y=(type)/FREE vec
			vec = wav[row][p][layer]
		elseif (rowGiven && layerGiven)		// want the layer
			N = DimSize(wav,2)
			Make/N=(N)/Y=(type)/FREE vec
			vec = wav[row][col][p]
		endif
	elseif (dims==4)
		if (colGiven && layerGiven && chunkGiven)		// want the row
			N = DimSize(wav,0)
			Make/N=(N)/Y=(type)/FREE vec
			vec = wav[p][col][layer][chunk]
		elseif (rowGiven && layerGiven && chunkGiven)	// want the column
			N = DimSize(wav,1)
			Make/N=(N)/Y=(type)/FREE vec
			vec = wav[row][p][layer][chunk]
		elseif (rowGiven && colGiven && chunkGiven)	// want the layer
			N = DimSize(wav,2)
			Make/N=(N)/Y=(type)/FREE vec
			vec = wav[row][col][p][chunk]
		elseif (rowGiven && colGiven && layerGiven)	// want the chunk
			N = DimSize(wav,3)
			Make/N=(N)/Y=(type)/FREE vec
			vec = wav[row][col][layer][p]
		endif
	endif
	if (!WaveExists(vec))
		return cmplx(NaN,NaN)
	endif

	WaveStats/M=1/Q vec
	WaveClear vec
	return cmplx(V_min,V_max)
End


Function FindScalingFromVec(vecIN,threshold,first,stepSize,dimN)
	// find the scaling suitable for a SetScale command that fits the values of vecIN[]
	// vecIN[] are the positions of a step scan.
	// first & stepSize are rounded to show no significant digits past threshold/10
	// dimN (number of points) does not count steps when vecIN does not change, this occurs in outer dim of a 2D scan
	// note, this should properly handle zig-zag scans too.
	Wave vecIN
	Variable threshold					// a change greater than this is intentional (less is jitter)
	Variable &first						// use 	SetScale/P x first,stepSize,"",waveName
	Variable &stepSize					// returned |step size|
	Variable &dimN							// returned number of points (=steps+1), does not count when standing still

	threshold = abs(threshold)							// no negative thresholds
	first = NaN													// init to bad values
	stepSize = NaN
	dimN = 0
	if (!WaveExists(vecIN) ||numpnts(vecIN)<1)		// failed
		return 1
	endif
	Variable i,m, N=numpnts(vecIN)
	Make/N=(N)/D/FREE vec=vecIN							// ensures that vec is a floating type

	// get the step size
	stepSize = FindStepSizeInVec(vec,threshold,signed=0)	// returns the un-signed |step size| from values in vec
	if (numtype(stepSize))									// failed
		stepSize = NaN
		return 1
	elseif (stepSize==0)									// only one point, no actual scan
		first = vecIN[0]
		stepSize = 0
		dimN = 1
		return 0
	endif

	// find the starting point, the median of all the points within |stepSize/2| of the lowest
	// note that first is always the smallest, regardless of the direction of scan, stepSize is always positive
	Duplicate/FREE vec vecSort
	Sort vecSort, vecSort
	i = floor(BinarySearchInterp(vecSort, vecSort[0]+stepSize/2))
	first = vecSort(i/2)									// this gives median of the points within |stepSize/2| of the lowest
	if (threshold!=0)											// use threshold/10 to round first
		Variable num = threshold/10
		Variable times = num<1 ? 10^ceil(-log(num)) : 10^floor(log(num))
		first = round(first*times)/times
	endif
	WaveClear vecSort

	Make/N=(N-1)/D/FREE dVec								// ensures that vec is a floating type
	dVec = vec[p+1]-vec[p]									// make dVec the differences
	vec = NaN
	vec[0] = vecIN[0]
	for (i=1,m=1; i<N; i+=1)								// set vec to be vecIN, but WITHOUT any constant steps
		if (abs(dVec[i-1]) > threshold)	
			vec[m] = vecIN[i]									// copy vecIN --> vec, but skip constant steps
			m += 1
		endif
	endfor
	N = m
	Redimension/N=(N) vec									// vec now contains NO contant steps, and NO NaN's

	// find number of points (=steps+1) in "one scan"
	//   does not count steps when vec does not change, this occurs in outer dim of a 2D scan
	//   a new scan is a change of more than |2*stepSize|
	//   a step is a change of more than |stepSize/5|
	Duplicate/FREE vec dVec
	Redimension/N=(N-1) dVec
	dVec = vec[p+1]-vec[p]									// make dVec the differences
	dVec = numtype(dVec) ? NaN : dVec					// change Inf --> NaN
	Variable tol = max(abs(stepSize/5),threshold)	// must be greater than this to look like a real step
	dVec = abs(dVec)<tol ? NaN : dVec					// step is a repeat, do not count it

	// test for ZigZag scans
	Duplicate/FREE dVec, ztest
	ztest = abs(dVec)>(2*stepSize) ? NaN : dVec	// remove all steps that look like snap-backs
	WaveStats/Q ztest											// average value of all normal scan steps
	Variable istart, Nsizes, isZigZag
	isZigZag = V_sdev>abs(2*V_avg) && (V_sdev/stepSize)>0.5
	WaveClear ztest
	if (isZigZag)
		// count the number of points in each scan, stored in sizes[]
		Make/N=(N)/U/I/FREE sizes=0						// holds the size of each scan that was found
		Variable plusMinus=sign(dVec[0])
		for (i=0,istart=0,Nsizes=0; i<(N-1); i+=1)
			if (sign(dVec[i]) != plusMinus)				// found end of a scan, process
				WaveStats/M=1/Q/R=[istart,i-1] dVec	// want number of valid (not NaN) points in this range of a single scan
				sizes[Nsizes] = V_npnts+1					// number of valid (not NaN) points in this scan [istart,i]
				plusMinus *= -1
				Nsizes += 1
				istart = i
			elseif (i==(N-2))									// the last point is alwasy the end of a scan, process
				WaveStats/M=1/Q/R=[istart,i] dVec		// want number of valid (not NaN) points in this range of a single scan
				sizes[Nsizes] = V_npnts+1					// number of valid (not NaN) points in this scan [istart,i]
				Nsizes += 1
			endif
		endfor

	else
		// at this point, dVec is non-zero at the end of each scan, NaN at bad points, and 0 at normal steps
		// count the number of points in each scan, stored in sizes[]
		dVec = abs(dVec)<(2*stepSize) ? 0 : dVec	// set all normal steps to 0, step is: delta < (2*step)
		Make/N=(N)/U/I/FREE sizes=0						// holds the size of each scan that was found
		for (i=0,istart=0,Nsizes=0; i<(N-1); i+=1)
			if (dVec[i])										// found a scan end
				WaveStats/M=1/Q/R=[istart,i] dVec		// want number of valid (not NaN) points in this range of a single scan
				sizes[Nsizes] = V_npnts					// number of valid (not NaN) points in this scan [istart,i]
				Nsizes += 1
				istart = i+1
			endif
		endfor

		if ((N-2-istart) > 0)								// add a last point, since there was probably not a big step at end
			WaveStats/M=1/Q/R=[istart,N-2] dVec
			sizes[Nsizes] = V_npnts+1
			Nsizes += 1
		elseif (Nsizes<1)
			dimN = 2
			return 0
		endif
	endif

	Redimension/N=(Nsizes) sizes
	sizes = !sizes ? NaN : sizes							// remove all zeros
	WaveStats/M=1/Q sizes
	Redimension/N=(V_npnts) sizes
	dimN = StatsMedian(sizes)								// the median number of points in one scan
	return 0
End


Function FindStepSizeInVec(vec,threshold,[signed])
	// find the step size in vec (position of a step scan) by examining the values
	// returns the step size from values in vec
	// stepSize is rounded to show no significant digits past threshold/10
	Wave vec
	Variable threshold								// changes greater than this are intentional (less is jitter)
	Variable signed									// if True, return signed stepSize, use signed=0 for zig-zag scans or SetScale commands
	signed = numtype(signed) ? 0 : signed		//		default is UNsigned

	if (!WaveExists(vec))							// failed, ERROR
		return NaN
	elseif (numpnts(vec)<2)						// do not have 2 point, just return step size of 0
		return 0
	endif

	threshold = abs(threshold)					// no negative thresholds
	if (numtype(threshold))						//  and you must pass a valid looking threshold
		return NaN
	endif

	Duplicate/FREE vec dVec
	Variable N = numpnts(vec)-1
	Redimension/N=(N) dVec
	SetScale/P x 0,1,"", dVec

	dVec = vec[p+1]-vec[p]							// make dVec the differences
	if (!signed)
		dVec = abs(dVec)								// abs() will show correct step size for signed, and for SetScale commands
	endif
	dVec = abs(dVec[p])<threshold ? NaN : dVec[p]
	Sort dVec, dVec

	WaveStats/Q dVec
	N = V_npnts
	Redimension/N=(N) dVec							// remove all NaN's (tiny steps) 
	if (N<1)												// do not have 2 point, just return step size of 0
		return 0
	endif

	Variable stepSize=dVec[floor((N-1)/2)]	// take the median value (avoids problems with average)
	if (threshold!=0)	
		Variable num = threshold/10				// use threshold/10 to round
		Variable times = num<1 ? 10^ceil(-log(num)) : 10^floor(log(num))
		stepSize = round(stepSize*times)/times
	endif

	return stepSize
End


Function/C RangeOfValuesContainingFraction(wwIN,frac, [zero])	// this is a 1D calculation
	// returns range of VALUES of wave that constitute a given fraction of all values in wave
	// if zero is True, then values are assumed to be all positive, and returns rad s.t. frac of wwIN are < rad
	Wave wwIN							// a 1D wave
	Variable frac					// want range that contains frac, either [lo,hi], or [0,rad]
	Variable zero					// if True, then want only fraction starting from zero, assumes no negatives
	zero = ParamIsDefault(zero) || numtype(zero) ? 0 : zero

	if (frac<=0 || numtype(frac)==2)
		return cmplx(NaN,NaN)
	endif
	frac = min(frac,1)

	Variable ylo, yhi
	ylo = fractionalLevelFind(wwIN,frac)
	yhi = zero ? NaN : fractionalLevelFind(wwIN,1-frac)
	return cmplx(ylo,yhi)
End
//Function/C RangeOfValuesContainingFraction(wwIN,frac, [zero])	// this is a 1D calculation
//	// returns range of VALUES of wave that constitute a given fraction of all values in wave
//	// if zero is True, then values are assumed to be all positive, and returns rad s.t. frac of wwIN are < rad
//	Wave wwIN							// a 1D wave
//	Variable frac						// want range that contains frac, either [lo,hi], or [0,rad]
//	Variable zero						// if True, then want only fraction starting from zero, assumes no negatives
//	zero = ParamIsDefault(zero) || numtype(zero) ? 0 : zero
//
//	if (frac<=0 || numtype(frac)==2)
//		return cmplx(NaN,NaN)
//	endif
//	frac = min(frac,1)
//
//	Duplicate/FREE wwIN, ww
//	Sort ww, ww
//	WaveStats/M=1/Q ww
//	if (V_npnts<2)
//		return cmplx(NaN,NaN)
//	endif
//	Redimension/N=(V_npnts) ww				// this is needed since wwIN may contain NaN
//
//	Variable N=DimSize(ww,0), xlo=NaN, xhi=NaN
//	Variable x0=(1/N), dx=(1-1/N)/(N-1)	// a scaling for ww
//	if (zero)
//		Variable i = (frac-x0)/dx + 0.5
//		if (i<0)
//			xlo = ww[0]/2	
//		else
//			i = limit(i,0,N-1)
//			xlo = ww[i]									// xhi is left as NaN
//		endif
//	else
//		Variable lo = ((1-frac)/2 - x0)/dx - 0.5
//		Variable hi = ((1+frac)/2 - x0)/dx + 0.5
//		lo = limit(lo,0,N-1)
//		hi = limit(hi,0,N-1)
//		xlo = ww[lo]
//		xhi = ww[hi]
//	endif
//	return cmplx(xlo,xhi)
//End


ThreadSafe Function fractionalLevelFind(wwIN,frac)
	// find the val such that frac of the values in wwIN are < val
	// e.g. if frac is 0.05, then find val such that 5% of all elements in ww have values <= val
	// you can use this to get the [5%, 95%] range by calling it twice
	// Note, this ignores NaN's
	//
	// I use this routine since sorting to find the median can take a long time for large waves.
	// for 8e6 points, this routine take 1 sec, sorting takes 30 sec
	// at 1000 pnts they are about the same, so I Sort for smaller waves
	Wave wwIN					// wave with values
	Variable frac			// fraction to exclude

	Duplicate/FREE wwIN, ww							// don't change the original ww
	Redimension/N=(numpnts(ww)) ww
	SetScale/P x,0,1,"" ww
	Note/K ww
	WaveStats/Q/M=1 ww
	Variable lo=V_min, hi=V_max, Ntotal=V_npnts
	if (numtype(frac) || frac<0 || Ntotal<1)
		return NaN
	elseif (frac==0)
		return lo
	elseif (frac>=1)
		return hi
	endif

	Variable Nfind = round(frac*Ntotal)			// number of points that should have ww[i] > val
	if (Ntotal<1001)									// for wave smaller than 1000, sorting is faster
		Sort ww,ww
		return ww[limit(Nfind-0.5, 0, Ntotal-1)]
	endif
	Variable val, Nval	, minDelta=(hi-lo)*1e-10
	Variable count=0
	do
		count += 1
		val = (lo + hi)/2								// trial value
		MatrixOp/FREE Nval0 = sum(greater(ww,val))
		Nval = Ntotal - Nval0[0]						// number below val
		if (Nval > Nfind)								// number found is too big, this is the new hi
			hi = val
		elseif (Nval < Nfind)							// number found is too small, this is the new lo
			lo = val
		endif
	while (Nval!=Nfind && (hi-lo)>minDelta && count<40)
	return val
End
//	Function test_fractionalLevelFind(N)
//		Variable N
//		Variable val, timer, ss0,ss1
//		Make/N=(N)/D/FREE wav=enoise(1)
//		wav += 1
//	
//		timer = StartMSTimer
//		val = fractionalLevelFind(wav,0.3)
//		ss0 = StopMSTimer(timer )
//		printf "for N=%d,  val=%g,   time=%g\r",N,val,ss0*1e-6
//	
//		timer = StartMSTimer
//		Sort wav, wav
//		ss1 = StopMSTimer(timer )
//		printf "sorting took %g,   %.2g times longer\r",ss1*1e-6, ss1/ss0
//	End


// This routine is much faster than going through an [sprintf str,"%g",val] conversion
ThreadSafe Function roundSignificant(val,N)	// round val to N significant figures
	Variable val			// input value to round
	Variable N			// number of significant figures

	if (val==0 || numtype(val))
		return val
	endif
	Variable is,tens
	is = sign(val) 
	val = abs(val)
	tens = 10^(N-floor(log(val))-1)
	return is*round(val*tens)/tens
End

ThreadSafe Function placesOfPrecision(a)	// number of significant figures in a number (at most 16)
	Variable a
	a = roundSignificant(abs(a),17)
	Variable i
	for (i=1;i<18;i+=1)
		if (abs(a-roundSignificant(a,i))/a<1e-15)
			break
		endif
	endfor
	return i
End


ThreadSafe Function/T ValErrStr(val,err,[sp])	// returns string  "val � err"
	Variable val							// value
	Variable err							// err in value
	Variable sp								// optionally put spaces around the �
	sp = ParamIsDefault(sp) ? 0 : !(!sp)
	sp = numtype(sp) ? 0 : sp
	err = numtype(err) ? 0 : abs(err)

	Variable n = ceil(log(abs(val/err)))+1// number of places of precision to show
	n = numtype(n) ? 6  : n					// default precision is 6
	n = limit(n,1,15)

	String vfmt, evfmt, str
	vfmt = "%."+num2istr(n+1)+"g"				// format for value
	String pm = SelectString(sp,PLUSMINUS," "+PLUSMINUS+" ")
	evfmt = vfmt + pm + SelectString(n>=2,"%.1g","%.2g")

	if (err>0)
		sprintf str,evfmt,val,err				// have a valid err, show it
	else
		sprintf str,vfmt,val					// no valid err, just show value
	endif
	return str
End


ThreadSafe Function/S num2strFull(n,[fmt,tol])// like num2str, but for integers, gives full precision
	Variable n
	String fmt										// format to use for non-integers
	Variable tol									// tolerance for identifying integers
	tol = ParamIsDefault(tol) || tol<=0 || numtype(tol) ? 1e-13 : tol
	if (abs(mod(n,1))<tol)						// this is an integer
		return num2istr(n)						// return integer format
	endif
	String str
	fmt = SelectString(ParamIsDefault(fmt), fmt, "%.15g")
	sprintf str, fmt, n							// real format
	return str
End



ThreadSafe Function normalize(a)	// normalize a and return the initial magnitude
	Wave a
	Variable norm_a
	if (WaveDims(a)==1)											// for a 1-d wave, normalize the vector
		norm_a = norm(a)
	elseif(WaveDims(a)==2 && DimSize(a,0)==DimSize(a,1))	// for an (n x n) wave, divide by the determinant
		norm_a = MatrixDet(a)^(1/DimSize(a,0))
	endif
	if (norm_a==0 || numtype(norm_a))
		return 0
	endif

	if (WaveType(a)&1)											// for a complex wave
		FastOp/C a = (1/norm_a)*a								//	a /= norm_a
	else
		FastOp a = (1/norm_a)*a									//	a /= norm_a
	endif
	return norm_a
End



// returns 1 if ww is only positive integers (or 0), tolerance to an int is tol
//	This is useful for determining if ww are counts (only positive ints)
ThreadSafe Function isPositiveInts(ww,[tol])
	Wave ww
	Variable tol							// tolerance for deciding if a number is an integer
	tol = ParamIsDefault(tol) ? 0.0001 : tol
	if (numtype(tol) || tol <0)
		return 0
	elseif (WaveType(ww,1) != 1)		// check if a numeric wave, NOT counts
		return 0
	elseif (WaveMin(ww)<0)				// found negatives, NOT counts
		return 0
	elseif (WaveType(ww) & 0x38)		// an integer type, IS counts
		return 1
	elseif (WaveType(ww) & 0x01)		// complex, NOT counts
		return 0
	endif

	// so ww is a float type, check if all integers
	MatrixOp/FREE tw = ReplaceNaNs(ww,0)	// ensure no NaN, (Inf is passed)
	MatrixOp/FREE test = sum(greater(abs(round(tw)-tw),tol))
	return test[0]<1					// if test[0] is 0, then all are ints
End


Function arithmetic(expression,[def])
	String expression		// a string that evaluates to a NUMBER, e.g. "1/3 + 5" or "-1/3 + 5"
	Variable def			// default value, usually 1 or 0 or NaN
	def = ParamIsDefault(def) ? NaN : def
	if ((strsearch(expression,"x",0,2)+strsearch(expression,"y",0,2)+strsearch(expression,"z",0,2)) > -3)
		return def			// just x, can give incorrect result??
	endif
	String cmd
	sprintf cmd, "Variable val__temp__ = (%s)",expression
	Execute/Z cmd
	Variable val=def
	if (V_flag==0)
		val = NumVarOrDefault("val__temp__",def)
	endif
	KillVariables/Z val__temp__
	return val
End


ThreadSafe Function maxValueOfType(ww)		// returns largest number that can be stored in a wave of specified type
	Wave ww
	if (WaveType(ww,1) != 1)		// returns NaN for non-numeric waves
		return NaN
	endif

	switch(WaveType(ww))
		case 0x02:					// single float
		case 0x03:					// include complex
			return 3.40282356e+38
		case 0x04:					// double float
		case 0x05:					// include complex
			return 1.79769313486231e+308
		case 0x08:					// 8 bit int, 2^7 - 1
			return 127
		case 0x48:					// 8 bit unsigned int, 2^8 - 1
			return 255
		case 0x10:					// 16 bit int, 2^15 - 1
			return 32767
		case 0x50:					// 16 bit unsigned int, 2^16 - 1
			return 65535
		case 0x20:					// 32 bit int, 2^31 - 1
			return 2147483647
		case 0x60:					// 32 bit unsigned int, 2^32 - 1
			return 4294967295
	endswitch
	return NaN
End
//Function FindMaxValueOfDouble()
//	Variable current, lastOK=1e307, step=lastOK
//	do
//		current = lastOK + step
//		if (numtype(current)==0)
//			lastOK = current
//		else
//			step /= 2
//		endif
//	while (step/lastOK > 1e-15)
//	print/D lastOK, current, step, step/lastOK
//End
//
//Function FindMaxValueOfSingle()
//	Make/N=1/FREE single
//	Variable small=0, big=1e50, mid
//	do
//		mid = (small+big)/2
//		single[0] = mid
//
//		if (numtype(single[0]))
//			big = (small+big)/2
//		else
//			small = mid
//		endif
//	while ((big-small)/small > 1e-15)
//	print/D small,mid,big,single[0]
//	single[0] = small
//	print/D mid,single[0]
//End


ThreadSafe Static Function DefaultZeroThresh(ww)
	Wave ww
	switch(WaveType(ww) & 0x3E)
		case 0x02:						// 32 bit float
			return Smallest32bitFloat
		case 0x04:						// 64 bit float
			return Smallest64bitFloat
		default:							// all integer types
			return 0
	endswitch
End
//Function FindSmallestFloatMachine()
//	Variable last
//	Make/N=1/D/FREE double=1
//	do
//		last = double[0]
//		double[0] /= 2
//	while(double[0]!=0)
//	print/d last
//
//	Variable mult=0.9999999, factor=mult
//	double[0] = last
//	do
//		last = double[0]
//		double[0] *= factor
//		factor *= mult
//	while(double[0]!=0)
//	print/d last
//
//	print " "
//	Make/N=1/FREE single=1
//	do
//		last = single[0]
//		single[0] /= 2
//	while(single[0]!=0)
//	print/d last
//
//	factor=mult
//	single[0] = last
//	do
//		last = single[0]
//		single[0] *= factor
//		factor *= mult
//	while(single[0]!=0)
//	print/d last
//End


ThreadSafe Function OrderValues(x0,x1)
	// on return, x0 and x1 are optionally swapped so that x0<x1
	Variable &x0, &x1
	if (x0 > x1)
		Variable swap
		swap = x0
		x0 = x1
		x1 = swap
	endif
End


// return the cursor range from a graph as a cmplx(lo,hi) (cmplx(NaN,NaN) on failure)
Function/S GetCursorRangeFromGraph(gName)
	String gName				// name of graph, use "" for top graph

	String infoA=CsrInfo(A,gName), infoB=CsrInfo(B,gName)
	if (strlen(infoA)<1 || strlen(infoB)<1)
		return ""
	endif

	String trName=StringByKey("TNAME",infoA)
	Wave yData=TraceNameToWaveRef(gName, trName )
	if (!WaveExists(yData))
		return ""
	endif
	Wave xData=XWaveRefFromTrace(gName,trName)

	Variable pLo=pcsr(A,gName), pHi=pcsr(B,gName) 
	Variable xLo=hcsr(A,gName), xHi=hcsr(B,gName) 
	if (numtype(xlo+xhi+pLo+pHi))
		return ""
	endif
	OrderValues(pLo,pHi)
	OrderValues(xLo,xHi)

	String out=""
	out = ReplaceStringByKey("yData",out,trName,"=")
	if (WaveExists(xData))
		out = ReplaceStringByKey("xData",out,GetWavesDataFolder(xData,2),"=")
	endif
	out = ReplaceNumberByKey("pLo",out,pLo,"=")
	out = ReplaceNumberByKey("pHi",out,pHi,"=")
	out = ReplaceNumberByKey("xLo",out,xLo,"=")
	out = ReplaceNumberByKey("xHi",out,xHi,"=")
	return out
End


// These next three functions convert Igor Peak fitting parameters to useful numbers
ThreadSafe Function GaussianFWHM(W_coef)			// FWHM of a Gaussian Peak
	// y = K0+K1*exp(-((x-K2)/K3)^2)
	Wave W_coef
	if (WaveType(W_coef,1)==1)							// W_coef exists, and it is numeric
		return 2*W_coef[3]*sqrt(ln(2))
	endif
	return NaN
End


ThreadSafe Function LorentzianFWHM(W_coef)			// FWHM of a Lorentzian Peak
	// y = K0+K1/((x-K2)^2+K3)
	Wave W_coef
	if (WaveType(W_coef,1)==1)							// W_coef exists, and it is numeric
		return 2*sqrt(W_coef[3])
	endif
	return NaN
End


ThreadSafe Function LorentzianIntegral(W_coef)	// peak integral from a Lorentzian fit
	// y = K0+K1/((x-K2)^2+K3)
	Wave W_coef
	if (WaveType(W_coef,1)==1)							// W_coef exists, and it is numeric
		return PI*W_coef[1]/sqrt(W_coef[3])
	endif
	return NaN
End

ThreadSafe Function GaussianIntegral(W_coef)		// peak integral from a Gaussian fit
	// y = K0+K1*exp(-((x-K2)/K3)^2)
	Wave W_coef
	if (WaveType(W_coef,1)==1)							// W_coef exists, and it is numeric
		return W_coef[1] * abs(W_coef[3]) * sqrt(PI)
	endif
	return NaN
End

ThreadSafe Function Gaussian2DIntegral(W_coef)	// peak integral from a 2D-Gaussian fit
	// z = K0+K1*exp((-1/(2*(1-K6^2)))*(((x-K2)/K3)^2 + ((y-K4)/K5)^2 - (2*K6*(x-K2)*(y-K4)/(K3*K5))))
	Wave W_coef
	if (WaveType(W_coef,1)==1)							// W_coef exists, and it is numeric
		Variable K1=W_coef[1], K3=W_coef[3], K6=W_coef[6], K5=W_coef[5]
		return 2*PI * K1 * abs(K3*K5) * sqrt(1-K6^2)
	endif
	return NaN
End



ThreadSafe Function computeCOM(ywave,xwave,[pLo,pHi])	// computes center of mass,  xwave is optional (use $"" if no xwave)
	Wave ywave,xwave
	Variable pLo,pHi										// optional range to check
	Variable com=NaN
	if (WaveType(ywave,1)!=1)							// ywave must exist, and must be numeric
		return NaN
	endif
	Variable N=DimSize(ywave,0)
	pLo = ParamIsDefault(pLo) || numtype(pLo) ? 0 : pLo
	pHi = ParamIsDefault(pHi) || numtype(pHi) ? N-1 : pHi

	MatrixOP/FREE yy = ReplaceNaNs(ywave,0)

	if (pLo>0)
		yy[0,pLo-1] = 0
	endif
	if (pHi < (N-1))
		yy[pHi,N-1] = 0
	endif
	if (WaveExists(xwave))
		yy = numtype(xwave) ? 0 : yy
		MatrixOP/FREE comW = sum(ReplaceNaNs(yy*xwave,0)) / sum(yy))
		com = comW[0]
	else
		Variable syy = sum(yy)
		CopyScales ywave,yy
		yy *= x
		com = sum(yy)/syy
	endif
	return com
End



Function/S FitErrorString(FitError,FitQuitReason)
	Variable FitError,FitQuitReason			// values of V_FitError  &  V_FitQuitReason

	if (FitError==0 && FitQuitReason==0)
		return ""
	endif

	String out=""
	out += SelectString(FitError & 2,"","Singular matrix; ")
	out += SelectString(FitError & 4,"","Out of memory; ")
	out += SelectString(FitError & 8,"","Function returned NaN or INF; ")
	out += SelectString(FitError & 16,"","Fit function requested stop; ")
	out += SelectString(FitError & 32,"","Reentrant curve fitting; ")
	out += StringFromList(FitQuitReason,";the iteration limit was reached;the user stopped the fit;the limit of passes without decreasing chi-square was reached;")
	out = SelectString(strlen(out),"Unspecified Fit Error",out)
	return out
End



Function/Wave PowerIntegerScale(image,power)	// scale the values of image by ^(power) AND for real values preserve the sign
	Wave image					// assumed to be an image, but works for waves of dimension of 1,2, & 3
	Variable power				// must be > 0, uses integer power or integer roots

	String wname, scaling
	Variable ipower
	if (numtype(power) || WaveType(image,1)!=1)
		return $""				// really bad input, power must be a number & image a numeric wave
	elseif (power>1)
		ipower = round(power)
		power = ipower
		wname = GetWavesDataFolder(image,2)+"_pow"+num2istr(ipower)
		sprintf scaling, "^%d",ipower
	elseif (power<1 && power>0)
		ipower = round(1/power)
		power = 1/ipower
		wname = GetWavesDataFolder(image,2)+"_root"+num2istr(ipower)
		sprintf scaling, "^(1/%d)",ipower
	else
		return $""				// power must be > 0
	endif
	if (numtype(ipower+power) || WaveType(image,1)!=1)
		return $""
	endif

	Variable type=WaveType(image)
	Variable itype = (type & 0x06) ? type : 0x02	// if a float or complex then the same type, otherwise single precision outout
	Make/N=(2)/Y=(itype)/O $wname/WAVE=imagePow
	if (type & 0x01)			// for complex images
		MatrixOp/O/C imagePow = powC(abs(image),power)	// don't try to fix sign for complex values
	else
		MatrixOp/O imagePow = powR(abs(image),power) * (2.0* greater(image,0.0) - 1.0)
	endif
	Redimension/Y=(itype) imagePow						// otherwise the MatrixOp always retrurn double
	CopyScales image, imagePow
	SetScale d 0,0,WaveUnits(imagePow,-1)+scaling, imagePow
	Note/K imagePow, ReplaceStringByKey("valueScaling",note(image),scaling,"=")
	return imagePow 
End



Function/T Posix2HFS(posixName,[printIt])	// This is a replacement for PosixToHFS()
	// Since we can replace HFSToPosix() with ParseFilePath(), we no longer need HFSAndPosix.xop
	String posixName
	Variable printIt									// normally returns "" on error, but this causes more printout
	printIt = ParamIsDefault(printIt) ? (strlen(GetRTStackInfo(2))==0) : printIt
	printIt = numtype(printIt) ? 0 : !(!printIt)

	String cmd
	sprintf cmd, "(POSIX file \"%s\") as Unicode text",posixName
#if (IgorVersion()<7)
	ExecuteScriptText/Z cmd
	S_value = TrimBoth(S_value,chars="\"")
#else
	ExecuteScriptText/UNQ/Z cmd
#endif
	if (V_flag)
		if (printIt)
			printf "failure in Posix2HFS(\"%s\")\t\tcmd = '%s'\r",posixName,cmd
			print S_value
		endif
		return ""										// there is no valid file path
	endif
	return S_value
End


Function pingHost(host)				// returns ping time in seconds, returns NaN CANNOT ping the host
	String host								// host name to ping
	if (!stringmatch(StringByKey("OS",IgorInfo(3)),"Macintosh OS X"))
		DoAlert 0, "Only know how to ping from a Mac"
		return NaN								// cannot get answer
	endif
	ExecuteScriptText/Z "do shell script \"ping -Q -t 1 -c 1 "+host+"\""

	Variable i = strsearch(S_value,"round-trip min/avg",0,2), Tavg
	if (i<0)
		return NaN
	endif
	i = strsearch(S_value,"=",i)
	if (i<0)
		return NaN
	endif
	Tavg = str2num(S_value[i+1,Inf])
	return Tavg
End

Function cpuFrequency()					// return the cpu frequency (Hz)
	if (!stringmatch(StringByKey("OS",IgorInfo(3)),"Macintosh OS X"))
		print "Only know how to get cpu speed from a Mac"
		// DoAlert 0, "Only know how to get cpu speed from a Mac"
		return NaN								// cannot get answer
	endif
	String cmd
	sprintf cmd "do shell script \"sysctl hw.cpufrequency\""
#if (IgorVersion()<7)
	ExecuteScriptText cmd						//returns something like: 	"hw.cpufrequency: 2000000000"
	S_value = TrimBoth(S_value,chars="\"")
#else
	ExecuteScriptText/UNQ cmd				//returns something like: 	"hw.cpufrequency: 2000000000"
#endif
	Variable freq = NumberByKey("hw.cpufrequency",S_value)
	if (numtype(freq))
		DoAlert 0, "Unable to get cpu frequency, message is '"+S_value+"'"
	endif
	return freq
End

Function/T systemUserName()				// return unix username
#if (IgorVersion()<7)
	if (!stringmatch(StringByKey("OS",IgorInfo(3)),"Macintosh OS X"))
		print "Only know how to get user name from a Mac"
		return ""								// cannot get answer
	endif
	String cmd
	sprintf cmd "do shell script \"whoami\""
	ExecuteScriptText cmd					//returns the user name in double quotes
	Variable i=strlen(S_value)-2
	String userName=S_value[1,i]
	if (i<=0)									// no name,
		DoAlert 0, "Unable to get user name, message is '"+S_value+"'"
	endif
	return userName
#else
	return IgorInfo(7)
#endif
End

#if (IgorVersion()>7 && stringmatch(IgorInfo(2),"Windows"))
Function/T sytemHostname()				// returns the full hostname as a string e.g. bob.xray.aps.anl.gov  (not ip address)
	String fqdn = GetEnvironmentVariable("COMPUTERNAME") + "." + GetEnvironmentVariable("USERDNSDOMAIN")
	return fqdn
End
#elif (stringmatch(IgorInfo(2),"Macintosh"))
Function/T sytemHostname()				// returns the full hostname as a string e.g. bob.xray.aps.anl.gov  (not ip address)
	if (!stringmatch(StringByKey("OS",IgorInfo(3)),"Macintosh OS X"))
		DoAlert 0, "Only know how to get hostname from a Mac"
		return ""								// cannot get answer
	endif
#if (IgorVersion()<7)
	ExecuteScriptText "do shell script \"host $(hostname)\""		//returns something like:  "bob.xray.aps.anl.gov"
	S_value = ReplaceString("\"",S_value,"")
#else
	ExecuteScriptText/UNQ/Z "do shell script \"host $(hostname)\""	//returns something like:  "bob.xray.aps.anl.gov"
	if (V_flag)
		return "ERROR -- "+S_value
	endif
#endif
	S_value = TrimBoth(S_value)
	return StringFromList(0,S_value," ")
End
#else
Function/T sytemHostname()				// returns the full hostname as a string e.g. bob.xray.aps.anl.gov  (not ip address)
	DoAlert 0, "Do not know how to get hostname from Windows in Igor 6 or earlier."
	return ""									// cannot get answer
End
#endif

Function/T localTimeZoneName()
	if (!stringmatch(StringByKey("OS",IgorInfo(3)),"Macintosh OS X"))
		print "Only know how to get Time Zone name from a Mac"
		return ""									// cannot get answer
	endif
#if (IgorVersion()<7)
	ExecuteScriptText "do shell script \"date +'%Z'\""			//returns something like: 	"CST"
	S_value = ReplaceString("\"",S_value,"")
#else
	ExecuteScriptText/UNQ "do shell script \"date +'%Z'\""			//returns something like: 	"CST"
#endif
	if (strlen(S_value)>3)
		DoAlert 0, "Unable to get local Time Zone name, message is '"+S_value+"'"
	endif
	return S_value
End

#if (stringmatch(IgorInfo(2),"Macintosh"))
// Get a list of system environment variables, semi-colon separated
Function/T getEnvironment(name)
	String name										// name of variable, use "=" to get all
	String cmd
	sprintf cmd "do shell script \"source ~/.bash_profile ; env\""
#if (IgorVersion()<7)
	ExecuteScriptText cmd
	S_value = TrimBoth(S_value,chars="\"")
#else
	ExecuteScriptText/UNQ cmd
#endif
	if (strsearch(S_value, ";", 0)>=0)
		DoAlert 0, "env has semi-colon"
		return ""
	endif
	S_value = ReplaceString("\r",S_value,";")

	if (StringMatch(name,"="))
		return S_value
	else
		return StringByKey(name,S_value,"=")
	endif
End
#elif (IgorVersion() >= 7)					// Windows in Igor7
Function/T getEnvironment(name)
	String name										// name of variable, use "=" to get all
	return GetEnvironmentVariable(name)
End
#else													// Windows in Igor6
Function/T getEnvironment(name)
	String name										// name of variable, use "=" to get all
	DoAlert 0, "Do not know how to get envirnment variable on Windows in Igor 6 or earlier"
	return ""										// cannot get answer
End
#endif


ThreadSafe Function FindFirstCharInStr(str,chars)
	// finds fist occurance of any of the chars in chars that occur in str
	// returns position of first one to occur, returns Inf if none are found
	String str						// the string to look in
	String chars					// the characers to search for

	String ch
	Variable i,j,Nc=strlen(chars), first=Inf
	for (i=0;i<Nc;i+=1)			// loop over each of the characters in chars
		ch = chars[i]
		j = strsearch(str,ch,0)
		first = j>=0 ? min(first,j) : first
	endfor
	return first					// will return Inf if no characters found
End


ThreadSafe Function/S RemoveLeadingString(str,head,ignoreCase)
	// returns str with head removed (if it starts with head)
	// if head=="", returns str
	// ignoreCase: 1=upper or lower case,   0=case must match
	String str				// string to search
	String head				// string that may be at start
	Variable ignoreCase	// if true then ignore case of head

	Variable strict = ignoreCase ? 3 : 1
	if (strsearch(str,head,Inf, strict)==0)	// found head at position 0
		return str[strlen(head),Inf]				// remove head
	endif
	return str
End

ThreadSafe Function/S RemoveTrailingString(str,tail,ignoreCase)
	// returns str with tail removed (if it ends with tail)
	// if tail=="", returns str
	// ignoreCase: 1=upper or lower case,   0=case must match
	String str				// string to search
	String tail				// string that may be at end
	Variable ignoreCase	// if true then ignore case of tail

	Variable strict = ignoreCase ? 3 : 1
	Variable i = strsearch(str,tail,Inf, strict)
	if ((i+strlen(tail)) == strlen(str))	// found tail at end
		str = str[0,strlen(str)-strlen(tail)-1]	// trim off tail
	endif
	return str
End

// These three functions are replacements for TrimFrontBackWhiteSpace() functions that are now DEPRECATED (see bottom of file)
ThreadSafe Function/S TrimBoth(str,[chars,ignoreCase])
	String str					// string to process
	String chars				// a set of characters to strip, defaults to white space if not given
	Variable ignoreCase		// True means NOT case sensitive, default is ignore case
	chars = SelectString(ParamIsDefault(chars), chars, "--WHITE SPACE--")
	ignoreCase = ParamIsDefault(ignoreCase) || numtype(ignoreCase) ? 2 : ignoreCase
	ignoreCase = ignoreCase ? 2 : 0
	str = TrimFront(str,chars=chars,ignoreCase=ignoreCase)
	str = TrimEnd(str,chars=chars,ignoreCase=ignoreCase)
	return str
End
//
ThreadSafe Function/S TrimFront(str,[chars,ignoreCase])		// remove specified leading chars
	String str					// string to process
	String chars				// a set of characters to strip, defaults to white space if not given
	Variable ignoreCase		// True means NOT case sensitive, default is ignore case
	chars = SelectString(ParamIsDefault(chars), chars, "--WHITE SPACE--")

	Variable i, N=strlen(str)
	if (StringMatch(chars,"--WHITE SPACE--"))
		for (i=0; char2num(str[i])<=32 && i<N; i+=1)	// find first non-white space
		endfor
	else
		ignoreCase = ParamIsDefault(ignoreCase) || numtype(ignoreCase) ? 2 : ignoreCase
		ignoreCase = ignoreCase ? 2 : 0
		for (i=0; strsearch(chars,str[i],0,ignoreCase)>=0 && i<N; i+=1)	// find first char not in chars
		endfor
	endif
	return str[i,Inf]
End
//
ThreadSafe Function/S TrimEnd(str,[chars,ignoreCase])		// remove specified trailing chars
	String str					// string to process
	String chars				// a set of characters to strip, defaults to white space if not given
	Variable ignoreCase		// True means NOT case sensitive, default is ignore case
	chars = SelectString(ParamIsDefault(chars), chars, "--WHITE SPACE--")

	Variable i, N=strlen(str)
	if (StringMatch(chars,"--WHITE SPACE--"))
		for (i=N-1; char2num(str[i])<=32 && i>=0; i-=1)	// find last non-white space
		endfor
	else
		ignoreCase = ParamIsDefault(ignoreCase) || numtype(ignoreCase) ? 2 : ignoreCase
		ignoreCase = ignoreCase ? 2 : 0
		for (i=N-1; strsearch(chars,str[i],0,ignoreCase)>=0 && i>=0; i-=1)	// find first char not in chars
		endfor
	endif
	return str[0,i]
End
//	DEPRECATED, the old functions: TrimFrontBackWhiteSpace(), TrimLeadingWhiteSpace(), and TrimTrailingWhiteSpace() are DEPRECATED


Function countChars(buf, chars)
	// count number of times one of the characters in chars occur in buf
	String buf				// a string buffer
	String chars			// string with characters to search for
	Variable i, m, N
	String char
	for (m=0, N=0; m<strlen(chars); m+=1)	// loop over all the characters in chars
		char = chars[m]								// a single char in chars
		i = -1
		do 
			i = strsearch(buf, char, i+1)		// count each occurance in buf
			N += i<0 ? 0 : 1
		while(i>=0)
	endfor 
	return N
End


Function PrintLongStrings(buf,[sep])	// prints really long strings, returns number of lines printed
	String buf
	String sep
	String white=" \t\r\n"
	sep = SelectString(ParamIsDefault(sep),sep,white)	// white space {space,tab,return,newline}
	String trim = SelectString(ParamIsDefault(sep), sep+white, sep)

	if (strlen(buf)<1)						// nothing
		return 0
	elseif (strlen(buf)<250)				// not that long
		print TrimBoth(buf,chars=trim)
		return 1
	endif

	String ch
	Variable m,i,imax, lines=0, Nch=strlen(sep)
	do
		buf = TrimBoth(buf,chars=trim)
		if (strlen(buf)<1)
			break
		elseif (strlen(buf)<=200)			// only a little bit left, print and end
			print TrimBoth(buf,chars=trim)
			lines += 1
			break
		endif

		for(m=0,imax=0; m<Nch; m+=1)		// loop over characters in sep to find best break
			i = strsearch(buf,sep[m],200,1)
			i = numtype(i) ? Inf : i
			imax = max(i,imax)
		endfor
		imax = imax<10 ? 200 : imax
		print TrimBoth(buf[0,imax],chars=trim)
		lines += 1
		buf = buf[imax+1,Inf]
	while (strlen(buf)>0)
	return lines
End


ThreadSafe Function/T IgorFileTypeString(itype)		// string associated with NUMTYPE from WaveInfo() function, or a WaveType() function
	Variable itype
	//	1:	Complex, added to one of the following:
	//	2:	32-bit (single precision) floating point
	//	4:	64-bit (double precision) floating point
	//	8:	8-bit signed integer
	//	16:	16-bit signed integer
	//	32:	32-bit signed integer
	//	64:	Unsigned, added to 8, 16 or 32 if wave is unsigned

	if (itype>127)
		return ""
	endif
	Variable s
	String stype=""

	stype += SelectString(itype&1,"","complex ")
	stype += SelectString(itype&2,"","single precision float (32-bit)")
	stype += SelectString(itype&4,"","double precision float (64-bit)")
	if (itype & 6)						// done with floating point type
		s = !!(itype&2) + !!(itype&4)
		stype = SelectString(s>1,stype,"")
		return stype
	endif

	s = !!(itype&8) + !!(itype&16) + !!(itype&32)
	if (strlen(stype) || s>1)
		return ""
	endif

	stype = SelectString(itype&64,"signed ","unsigned ")
	stype += SelectString(itype&8,"","8-bit integer ")
	stype += SelectString(itype&16,"","16-bit integer ")
	stype += SelectString(itype&32,"","32-bit integer ")
	return stype
End



Function/T GenericWaveNoteInfo(ww,key,[class,options,type])
	Wave ww
	String key
	String class							// default is all classes
	String options							// default is no options
	String type								// name of wave type (used in prompt), e.g. "image"
	class = SelectString(ParamIsDefault(class),class,"")
	options = SelectString(ParamIsDefault(options),options,"")
	type = SelectString(ParamIsDefault(type),type,"wave")

	Variable printIt = 0
	String wName, wList = WaveListClass(class,"*",options)
	if (!WaveExists(ww) && ItemsInList(wList)==1)
		Wave ww = $StringFromList(0,wList)
		printIt = 1
	endif
	if (!WaveExists(ww))
		String item=""
		Prompt wName,type+" to use", popup, wList
		DoPrompt "pick "+type, wName
		if(V_flag)
			return ""
		endif
		Wave ww=$wName
		printIt = 1
	endif
	wName = NameOfWave(ww)
	if (!WaveExists(ww))
		return ""
	endif
	String list = note(ww)

	FUNCREF MoreWaveNoteInfoProto more = $"MoreWaveNoteInfo"
	list = more(ww,list)

	if (strlen(key)>1)
		item = StringByKey(key,list,"=")
		if (strlen(item)>1 && numtype(str2num(item)) && printIt)
			printf "from '%s',  %s = '%s'\r",wName,key,item
		endif
	else
		Prompt item,"item",popup,list
		DoPrompt "info about "+wName, item
		if(V_flag)
			return ""
		endif
		key = StringFromList(0,item,"=")
		item = StringFromList(1,item,"=")
		printf "from '%s',  %s = %s\r",wName,key,item
	endif
	return item
End
//
Function/T MoreWaveNoteInfoProto(ww,list)		// additions to the list, only called by GenericWaveNoteInfo() in Utility_JZT.ipf
	Wave ww
	String list
	return list
End


Function StopAllTimers()
	Variable i
	String str=""
	for(i=0;i<=9;i+=1)
		str += SelectString(stopMSTimer(i),"",num2istr(i)+" ")
	endfor
	if (strlen(str))
		printf "timers %s were running\r",str
	endif
End


ThreadSafe Function/T dateStr2ISO8601Str(dateStr,[timeStr,zoneHr,zoneMin])	// convert input date string and converts to ISO8601 format
	String dateStr							// something like "4/2/2010" or "2010-4-2", if dateStr=="now", use now
	String timeStr
	Variable zoneHr, zoneMin				// the time zone in hours or minutes, OPTIONAL, (you can also use hour=6.5)
	timeStr = SelectString(ParamIsDefault(timeStr),timeStr,"")
	if (CmpStr(dateStr,"now")==0)			// special for dateStr = "now"
		return epoch2ISOtime(datetime)		// a shortcut for "now"
	endif

	Variable d1,d2,d3, dateSec
	dateStr = ReplaceString("-",dateStr,"/")
	sscanf dateStr,"%d/%d/%d",d1,d2,d3
	if (V_flag!=3)
		return ""								// failed to find date like string
	endif
	String outStr="", zoneStr=""
	if (d1>1903)									// assume year/month/day
		sprintf outStr,"%04d-%02d-%02d",d1,d2,d3
	else											// assume month/day/year
		sprintf outStr,"%04d-%02d-%02d",d3,d1,d2
	endif

	if (strlen(timeStr))						// time was also passed
		Variable hour,minute,second
		sscanf timeStr,"%d:%d:%d",hour,minute,second
		if (V_flag<2)							// invalid time of day
			return ""
		elseif(V_flag==2)
			second = 0							// no seconds given, assume 0
		endif
		Variable add12=(strsearch(timeStr,"p",0)>0) && (hour<12)
		hour += add12 ? 12 : 0
		outStr += "T"+Secs2Time(hour*3600+minute*60+second,3)
	endif

	hour = ParamIsDefault(zoneHr) ? 0 : zonehr	// decode the the time zone
	minute = round(ParamIsDefault(zoneMin) ? 0 : zoneMin)
	second = hour*3600 + 60*minute
	if (abs(second)>24*3600)					// time zone > 24h????
		return ""
	endif
	if ((!ParamIsDefault(zoneHr) || !ParamIsDefault(zoneMin)) && numtype(second)==0)	// have a time zone
		if (mod(second/3600,1))				// have minutes
			zoneStr = Secs2Time(abs(second),2)
			zoneStr = SelectString(second<0,"+","-")+zoneStr
		else
			sprintf zoneStr,"%+03d",round(second/3600)
		endif
	endif
	outStr += zoneStr

	return outStr
End


ThreadSafe Function ISOtime2IgorEpoch(iso,[NOTZ])	// convert ISO8601 string to an Igor Epoch (error returns NaN)
	// if no timezone info present, then the result is independent of timezone
	// if timezone info present, then returned seconds are in UTC
	// valid timezones are of the form: "+01", "+0130", "-01:30", "Z" (only Z, not other letters)
	// to get the old behavior, use new ISOtime2IgorEpoch(iso, NOTZ=1) instead
	String iso									// format     2013-10-02T01:22:35-06:00  (the seconds are optional, the time zone is optional)
	Variable NOTZ								// if True, then ignore any time zone info, just use local (this is like the old version)
	NOTZ = ParamIsDefault(NOTZ) || numtype(NOTZ) ? 0 : NOTZ		// default is interpret timezone info

	Variable year,month,day, hr,mn,se

	sscanf iso,"%4d-%2d-%2dT%2d:%2d:%2d", year,month,day,hr,mn,se
	Variable N=V_flag, UTC
	if (N<3)
		return NaN
	endif
	UTC = date2secs(year, month, day )
	UTC += N>=4 ? hr*3600 : 0
	UTC += N>=5 ? mn*60 : 0
	UTC += N>=6 ? se : 0

	Variable i0
	i0 = strsearch(iso,"T",i0)+1
	if (N==3 || NOTZ || i0<1)				// no timezone info, done
		return UTC
	endif

	if (isletter(iso[strlen(iso)-1]))
		if (cmpstr(iso[strlen(iso)-1],"Z"))
			return NaN							// the only allowed letter timezone is "Z"
		endif
	endif
	Variable ip,im,i1
	ip = strsearch(iso,"+",i0)
	im = strsearch(iso,"-",i0)
	i1 = max(ip,im)
	if (i1<0)
		return UTC								// not timezone info, done
	endif
	String tzStr=iso[i1,Inf]

	Variable tzH=0,tzM=0, tzLen=strlen(tzStr)
	if (tzLen==6)
		sscanf tzStr,"%3d:%2d", tzH,tzM
	elseif (tzLen==5)
		sscanf tzStr,"%3d%2d", tzH,tzM
	elseif (tzLen==3)
		sscanf tzStr,"%3d", tzH
	endif
	tzM *= im>0 ? -1 : 1
	tzH += (V_flag==2) ? tzM/60 : 0

	UTC -= tzH*3600
	return UTC
End


ThreadSafe Function/T ISOtime2niceStr(iso, [TZ])	// convert ISO8601 string to a nice format for graph annotations
	String iso					// format     2013-10-02T01:22:35  (the seconds are optional), if iso=="now", then use now
	Variable TZ					// flag, if true, include Time Zone
	TZ = ParamIsDefault(TZ) || numtype(TZ) ? 0 : TZ

	Variable N, seconds=NaN
	if (CmpStr(iso,"now")==0)					// special for iso == "now"
		seconds = DateTime
		N = 6
	else
		Variable year=NaN,month=NaN,day=NaN, hr=NaN,mn=NaN,se=NaN
		sscanf iso,"%4d-%2d-%2dT%2d:%2d:%2d", year,month,day,hr,mn,se
		N = V_flag
		if (N<3 || numtype(year+month+day))
			return ""
		endif
		seconds = date2secs(year, month, day)
		seconds += N>=5 ? hr*3600 + mn*60 : 0
		seconds += N>=6 ? se : 0
	endif

	String out = Secs2Date(seconds,2)
	if (N>=5)
		Variable fmt = (N>=6) ? 1 : 0
		out += SelectString(numtype(seconds),"  "+Secs2Time(seconds,fmt),"")
	endif
	if (TZ)
		out += " ("+num2str(date2secs(-1,-1,-1)/3600)+")"
	endif
	return out
End

ThreadSafe Function/T epoch2ISOtime(epoch,[localTZ])	// convert an Igor epoch (in seconds) to an ISO8601 format
	Variable epoch
	Variable localTZ											// flag if true, then add offset for local time zone
	localTZ = ParamIsDefault(localTZ) || numtype(localTZ) ? 0 : localTZ

	Variable s = mod(epoch,60)							// number of seconds
	Variable frac=roundSignificant(mod(epoch,1),5)// fractional seconds
	Variable places=placesOfPrecision(frac)			// places of precision in frac

	String out=Secs2Date(epoch,-2)+"T"					// date part with "T"
	if (s==0 && frac==0)
		out += Secs2Time(epoch,2)							// no seconds, just hr:min
	elseif (frac==0)
		out += Secs2Time(epoch,3)							// seconds, but no fracitonal seconds
	else
		out += Secs2Time(epoch,3,places)				// fractional seconds
	endif

	if (localTZ)												// add time zone info
		Variable tz = date2secs(-1,-1,-1)
		if (tz==0)
			out += "Z"
		else
			String str, ssign=SelectString(tz>0,"-","+")	
			tz = abs(tz)										// sign already taken care of by ssign
			Variable hr = trunc(tz/3600)
			Variable minutes = round(mod(tz/60,60))
			sprintf str,"%02d", tz/3600
			out += ssign+str
			if (minutes)			// minutes too
				sprintf str,":%02d",minutes
				out += str
			endif
		endif
	endif
	return out
End


Function AskForUserForDateTime(epoch)
	// puts up a dialog to select a date & time, returns the epoch (NaN on error or Cancel)
	// return NaN for invalid dates (e.g. Apr 31, or Feb 29 2005)
	Variable epoch
	epoch = (numtype(epoch) || epoch<=0) ? DateTime : epoch

	Variable day,month,year
	sscanf Secs2Date(epoch,-2), "%4d-%2d-%2d",year,month,day

	Variable midnight = epoch-date2secs(year,month,day)		// seconds since midnight(year,month,day)
	Variable hour=floor(midnight/3600), minute=floor(mod(midnight/60,60)), second=floor(mod(midnight,30))

	Prompt day,"Date", popup, expandRange("1-31",";")
	Prompt month,"Month", popup, MonthNamesFull
	Prompt year,"Year", popup, expandRange("1995-2025",";")
	Prompt hour,"Hours [0-23]",popup,expandRange("0-23",";")
	Prompt minute,"Minutes [0-59]",popup,expandRange("0-59",";")
	Prompt second,"Seconds [0-59]",popup,expandRange("0-59",";")
	year -= 1994
	hour += 1
	minute += 1
	second += 1
	DoPrompt "Date/Time",day,hour,month,minute,year,second
	if (V_flag)
		return NaN
	endif
	year += 1994
	hour -= 1
	minute -= 1
	second -= 1

	Variable err=0
	switch(month)
		case 2:			// February is special 28 days except leap year
			err = err || (mod(year,4) && day>28)		// on non-leap year, Feb only to 28
			err = err || (mod(year,4)==0 && day>29)	// on leap year, Feb goes to 29
		case 4:			// April, these four months have 30 days
		case 6:			// June
		case 9:			// September
		case 11:			// November
			err = err || day>30
	endswitch
	epoch = err ? NaN : date2secs(year,month,day) + 3600*hour + 60*minute + second
	return epoch
End


Function/T ElapsedTime2Str(seconds,[showSec,fracDigits])	// convert seconds to a nice time interval string
	Variable seconds			// number of seconds
	Variable showSec			// flag, if True, show seconds
	Variable fracDigits		// number of fractional seconds to show

	if (ParamIsDefault(showSec) && ParamIsDefault(fracDigits))
		if (seconds < 5*60)
			showSec = 1
		endif
		if (seconds < 2)
			fracDigits = 2
		elseif (seconds < 60)
			fracDigits = 1
		else
			fracDigits = 0
		endif
	else
		showSec = ParamIsDefault(showSec) || numtype(showSec) ? 0 : !(!showSec)
		fracDigits = ParamIsDefault(fracDigits) || numtype(fracDigits) ? 0 : fracDigits
		fracDigits = fracDigits<=0 || !showSec ? 0 : round(fracDigits)
	endif

	Variable years=trunc(seconds/(365*24*3600))	// integer number of years (365.0 days/year)
	seconds -= years * (365*24*3600)
	Variable weeks=trunc(seconds/(7*24*3600))		// integer number of weeks (<= 52)
	seconds -= weeks * (7*24*3600)
	Variable days=trunc(seconds/(24*3600))			// integer number of days (<= 7)
	seconds -= days * (24*3600)

	String str=""
	if (years>0)
		str += num2istr(years)+" yr,  "
	endif
	if (weeks>0)
		str += num2istr(weeks)+" wk,  "
	endif
	if (days>0)
		str += num2istr(days)+" d,  "
	endif
	if (seconds>0)
		showSec = showSec ? 5 : 4
		str += Secs2Time(seconds,showSec,fracDigits)
	endif
	str = RemoveEnding(str,",  ")
	return str
End


Function WeekDayNumber(month, idate, year, [printIt])	// from: https://cs.uwaterloo.ca/~alopez-o/math-faq/node73.html
	Variable month			// month in range [1, 12]
	Variable idate			// date in range [1,31], depending upon month
	Variable year			// year, use 2020, not just 20
	Variable printIt
	printIt = ParamIsDefault(printIt) || numtype(printIt) ? strlen(GetRTStackInfo(2))==0 : printIt
	Variable k=idate
	Variable m = month + (month<3 ? 10 : -2)
	Variable C = floor(year/100)
	Variable Y = mod(year,100)
	Y += (month < 3) ? 1 : 0				// move Jan & Feb into next year
	Variable W = (k + floor(2.6*m - 0.2) - 2*C +Y + floor(Y/4) + floor(C/4))
	W = mod(W,7) + 1

	if (printIt)
		printf "%s %d, %d  falls on a %s\r",StringFromList(month-1,MonthNamesFull), idate, year, StringFromList(W-1,DayNamesFull)
	endif
	return W										// weekday in [1,7]
End



ThreadSafe Function/T vec2str(w1,[places,fmt,maxPrint,bare,zeroThresh,sep])		// convert vector to s string suitable for printing, does not include name
	Wave w1										// 1d wave to print
	Variable places							// number of places, for default, use negative or NaN
	String fmt									// optional format for on number, this OVERRIDES places
	Variable maxPrint							// maximum number of elements to print, defaults to 20
	Variable bare								// if bare is TRUE, then suppress the "{}" in the output
	Variable zeroThresh						// |values| < zeroThresh show as a "0", in many vectors 1e-15 is really zero
	String sep									// optional separator, default is ",  "   a comma and 2 spaces

	maxPrint = ParamIsDefault(maxPrint) ? 20 : maxPrint
	maxPrint = maxPrint>0 ? maxPrint : 20
	places = ParamIsDefault(places) ? -1 : places
	fmt = SelectString(ParamIsDefault(fmt),fmt,"")
	bare = ParamIsDefault(bare) ? 0 : !(!bare)
	zeroThresh = ParamIsDefault(zeroThresh) || numtype(zeroThresh) || zeroThresh<=0 ? NaN : zeroThresh
	sep = SelectString(ParamIsDefault(sep),sep,",  ")

	if (!WaveExists(w1) || numpnts(w1)<1)
		return SelectString(bare,"{}","")
	endif

//	Wave/T tw=$GetWavesDataFolder(w1,2)
//	Wave/C cw=$GetWavesDataFolder(w1,2)
	Variable waveIsComplex = WaveType(w1) %& 0x01
	Variable numeric = (WaveType(w1)!=0)

	if (strlen(fmt))
		// fmt was passed, do not reset it
	elseif (waveIsComplex)
		Wave/C cw=w1
		places = places>=0 ? min(20,places) : 5	// default to 5 for unacceptable values
		sprintf fmt,"(%%.%dg, %%.%dg)",places,places
	elseif (numeric)
		places = places>=0 ? min(20,places) : 5	// default to 5 for unacceptable values
		sprintf fmt,"%%.%dg",places
	elseif (places>0)										// must be text, and a maximum length given
		Wave/T tw=w1
		sprintf fmt, "\"%d%%s\"",places
	else														// must be text with no preferred length
		Wave/T tw=w1
		fmt = "\"%%s\""
	endif

	Duplicate/FREE w1, wInternal
	if (numeric)
		if (!zeroThresh || numtype(zeroThresh))
			zeroThresh = DefaultZeroThresh(w1)
		endif
		wInternal = abs(wInternal)<zeroThresh ? 0 : wInternal
	endif

	Variable i=0, n
	n = numpnts(wInternal)
	maxPrint = min(n,maxPrint)
	String str, out=SelectString(bare,"{","")

	do
		if (waveIsComplex)						// a complex wave
			sprintf str,fmt, real(cw[i]),imag(cw[i])
		elseif (numeric && (!waveIsComplex))	// a simple number wave
			sprintf str,fmt, wInternal[i]
		elseif (!numeric)							// a text wave
			sprintf str,"\"%s\"", tw[i]
		endif
		out += str
		if (i<(n-1))
			sprintf str,sep
			out += str
		endif
		i += 1
	while (i<maxPrint)
	if (n>maxPrint)
		sprintf str,"...}\ronly printed %d of %d values\r",maxPrint,n
		out += str
	else
		out += SelectString(bare,"}","")
	endif
	return out
End


ThreadSafe Function/WAVE str2vec(str,[sep])// returns a free vector based on the string
	String str
	String sep											// optional separator, if not included will figure it out
	sep = SelectString(ParamIsDefault(sep),sep,"")

	if (strlen(str)<1)
		return $""
	elseif (strlen(sep)>1)							// sep can only be 1 character long
		return $""
	endif

	Variable i0=strsearch(str,"{",0), i1=strsearch(str,"}",0)
	if (i1>(i0+1) && i0>=0)						// vec is in {...} form, remove leading "{" and trailing "}"
		str = str[i0+1,i1-1]
	endif

	if (strlen(sep)==0)								// determine the separator and set it so semi-colon
		if (strsearch(str,";",0)>=0)
			str = ReplaceString(",",str,";")	// change comma separators to semi-colons, this represents a mixture of separators
		elseif (strsearch(str,",",0)>=0)
			str = ReplaceString(",",str,";")	// change comma separators to semi-colons
		elseif (strsearch(str,"\t",0)>=0)
			str = ReplaceString("\t",str,";")	// change tab separators to semi-colons
		else												// assume space separators
			str = TrimFrontBackWhiteSpace(str)	// remove leading & trailing white space
			do
				str = ReplaceString("  ",str," ")	// change all multi-space runs to single spaces
			while (strsearch(str,"  ",0)>=0)
			str = ReplaceString(" ",str,";")	// change space separators to semi-colons
		endif
		sep = ";"
	endif

	Variable N=ItemsInList(str,sep)				// number of values in str
	Make/N=(N)/FREE/D ww
	ww = str2num(StringFromList(p,str,sep))
	return ww
End
//
//	Function test_str2vec(str,[sep])			// try test_str2vec("1  3 2",sep="  ")
//		String str
//		String sep
//	
//		if (ParamIsDefault(sep))
//			Wave vec = str2vec(str)
//		else
//			Wave vec = str2vec(str,sep=sep)
//		endif
//		printf "'%s' --> %s\r",str,vec2str(vec)
//	End


// The routines encodeMatAsStr() and decodeMatFromStr() are NOT intended for viewing by user, but for a wave note
ThreadSafe Function/T encodeMatAsStr(mat,[places,vsep,t])	// write a string interpretable by decodeMatFromStr
	Wave mat
	Variable places		// defaults to 15 places of precisiont
	String vsep			// optional separator between vectors, use "," for Igor readable strings (the default)
	Variable t			// t=1 write the transpose of mat, only valid if 2D mat, default is t=0
	places = ParamIsDefault(places) || !(places>0) ? 15 : places
	vsep = SelectString(ParamIsDefault(vsep),vsep,",")
	t = ParamIsDefault(t) || numtype(t) ? 0 : t

	Variable i, Nr=DimSize(mat,0), Nc=DimSize(mat,1)
	if (Nc<1)
		return vec2str(mat,places=places,sep=",")
	elseif (t && Nc>1)								// transpose was requested for a 2D mat
		MatrixOP/FREE mfree = mat^t
		Variable swap = Nc							// swap the dimensions
		Nc = Nr
		Nr = swap
	else
		Duplicate/FREE mat, mfree
	endif

	Make/N=(Nr)/D/FREE vec
	String str="{"
	for (i=0;i<Nc;i+=1)
		vec = mfree[p][i]
		str += SelectString(i,"",vsep)			// no separator before first vector
		str += vec2str(vec,places=places,sep=",")
	endfor
	str += "}"
	return str
End


// This is also may be used as a replacement for str2recip() in LatticeSym.ipf
ThreadSafe Function/WAVE decodeMatFromStr(strIn, [t])	// returns a FREE wave defined by str
	// strIn looks something like "{{1.3,0,0},{0,1.3,0},{0,0,1.3},{3.14,2,19.666}}"
	//    or "{1.3,0,0}{0,1.3,0}{0,0,1.3}" is also OK
	String strIn
	Variable t			// t=1 return the transpose of mat, only valid if 2D mat, default is t=0
	t = ParamIsDefault(t) || numtype(t) ? 0 : t
	strIn = ReplaceString(" ",strIn,"")		// remove all spaces
	strIn = ReplaceString("\t",strIn,"")		//   and tabs
	strIn = ReplaceString("},{",strIn,";")	// use ';' to separate the column vectors
	strIn = ReplaceString("}{",strIn,";")	//   accept either "},{" or "}{"
	strIn = ReplaceString("{",strIn,"")		// remove leading "{"
	strIn = ReplaceString("}",strIn,"")		//   and remove and trailing "}"

	String str=StringFromList(0,strIn)
	Variable ic, Nr=ItemsInList(str,","), Nc=ItemsInList(strIn)	// number of rows, columns
	if (!(Nr*Nc>0))
			return $""								// Nr & Nc must both be non-zero
	endif
	Make/N=(Nr,Nc)/D/FREE mat=NaN				// mat is the result, fill it
	for (ic=0;ic<Nc;ic+=1)
		Wave vec = str2vec(StringFromList(ic,strIn))
		if (!(numpnts(vec)==Nr))					// every vector should have length of Nr
			return $""
		endif
		mat[][ic] = vec[p]
	endfor
	if (Nc==1)										// a 1D wave
		Redimension/N=(Nr) mat
	elseif (Nc>1 && t)								// transpose requested for a 2D mat
		MatrixTranspose mat
	endif
	return mat
End


ThreadSafe Function/T cmplx2str(zz,[places,pow])		// convert complex number to a printable string
	Variable/C zz
	Variable places
	Variable pow						// optional parameter for magnitude, adds |zz|^pow
	places = ParamIsDefault(places) ? NaN : places
	places = round(places)
	pow = ParamIsDefault(pow) ? 0 : pow
	pow = numtype(pow) ? 0 : pow

	Variable powValue=NaN
	String str, fmt1="%g", fmt="(%g, %g)", powStr=""
	if (places == limit(places,0,20))
		fmt1 = "%."+num2istr(places)+"g"
	endif
	if (pow)
		powStr = SelectString(pow==1, "^"+num2str(pow), "")
		powValue = sqrt(magsqr(zz))^pow
		sprintf fmt,"|(%s, %s)|%%s = %s",fmt1,fmt1,fmt1
		sprintf str,fmt,real(zz),imag(zz),powStr,powValue
	else
		sprintf fmt,"(%s, %s)",fmt1,fmt1
		sprintf str,fmt,real(zz),imag(zz)
	endif

	return str
End


ThreadSafe Function/C str2cmplx(str)	// this is like str2num, but for complex
	String str

	str = ReplaceCharacters("()[]{}",str,"")
	str = ReplaceCharacters(",;:",str," ")
	Variable rr,ii
	sscanf str,"%g %g",rr,ii
	if (V_flag==2)
		return cmplx(rr,ii)
	endif
	return cmplx(NaN,NaN)					// fail
End


ThreadSafe Function/T num2fraction(val,tol,[addPlus])	// turn val into a fraction string, 0.25 --> "1/4"
	Variable val
	Variable tol					// requested tolerance
	Variable addPlus				// if True, always include the '+' sign
	addPlus = ParamIsDefault(addPlus) || numtype(addPlus) ? 0 : addPlus

	if (numtype(val+tol)==2)
		return "NaN"
	elseif (numtype(val)==1)
		return num2str(val)
	elseif (val==0)
		return "0"
	endif

	RatioFromNumber/MERR=(tol) abs(val)
	Variable numer=V_numerator, denom=V_denominator
	String strSign = SelectString(addPlus && val>0, "", "+")
	strSign = SelectString(val<0,strSign,"-")
	String numerStr = strSign + num2istr(numer)

	if (numtype(numer)==2)
		return "NaN"
	elseif (denom==1)
		return numerStr
	endif
	return numerStr+"/"+num2istr(denom)
End
//
//ThreadSafe Function/C real2rational(val,maxDenom)
//	// for val, find closest rational number real/imag, where real & imag are ints
//	// note, for maxDenom=9e7, I get a memory fault, but 8e7 takes 5 seconds, 1e6 is only 1/16 sec.
//	Variable val					// the real value
//	Variable maxDenom				// maximum denominator to consider, allowed maximum of 1e6
//
//	if (maxDenom>8e7 || numtype(maxDenom) || maxDenom<1 || numtype(val)==2)
//		return cmplx(NaN,NaN)
//	elseif (numtype(val)==1)
//		return cmplx(val,1)
//	endif
//	Make/N=(maxDenom-1)/D/FREE pp=p+1
//	MatrixOP/FREE err = Abs((pp*val) - round(pp*val))
//	WaveStats/Q/M=1 err
//	Variable denom = V_minloc+1
//	Variable numer = round(val*denom)
//	denom = (numer==0) ? 1 : denom
//	return cmplx(numer,denom)	// numer/denom is closest rational to val
//End


ThreadSafe Function/S vec2MINstr(vecIN)		// a replacement for hkl2str()
	// format values in vecIN into a string of acceptable minimal length
	Wave vecIN
	if (numtype(sum(vecIN)))
		return ""
	elseif (numpnts(vecIN)>50)
		return ""
	endif

	Variable tooSmall = 0.5				// 1+tooSmall == 1,  (0.5 works for integer waves)
	tooSmall = (WaveType(vecIN) & 0x02) ? 1e-7 : tooSmall	// 32 bit float
	tooSmall = (WaveType(vecIN) & 0x04) ? 1e-15 : tooSmall	// 64 bit float

	Duplicate/FREE vecIN, vec
	vec = abs(vec)<tooSmall ? 0 : vec
	MatrixOP/FREE deltas = maxVal(Abs(round(vec)-vec))
	Variable isInts = deltas[0] < tooSmall
	MatrixOP/FREE biggest = maxVal(Abs(vec))

	if (isInts && biggest[0]<10)			// just single digit ints, no separator needed
		return vec2str(vec,fmt="%d",bare=1,maxPrint=50,sep="")
	endif
	return vec2str(vec,places=6,bare=1,maxPrint=50,zeroThresh=1e-6,sep=" ")
End
//
//	Function test_vec2MINstr()
//		Make/FREE vec={9,2,3,2e-8}
//		printf "%s --> %s     (hkl2str = %s)\r", vec2str(vec), vec2MINstr(vec), hkl2str(vec[0],vec[1],vec[2])
//	End


ThreadSafe Function/WAVE minStr2Vec(inStr,Nreq)		// a replacement for str2hkl()
	// returns the numeric values from a string, pretty forgiving about format in string
	// moved to here from Dynamical.ipf versions <=1.14
	// moved to here from LatticeSym.ipf, 5.16
	String inStr
	Variable Nreq				// rquired number of numbers to find

	inStr = ReplaceString("+",inStr," ")		// change all '+' to space
	inStr = ReplaceString("-",inStr," -")	// change all '-' to space+'-'
	inStr = ReplaceString(" ",inStr,";")		// change all ' ' to semi-colon
	inStr = ReplaceString(",",inStr,";")		// change all ',' to semi-colon
	inStr = ReplaceString(":",inStr,";")		// change all ':' to semi-colon
	String digits = "0123456789"
	Variable semicolon = 0x3B, period = 0x2E
	String skips = "-.;"
	for (; strsearch(inStr,";;",0)>=0; )		// change all ";;" to a single ";"
		inStr = ReplaceString(";;",inStr,";")
	endfor
	inStr = TrimBoth(inStr,chars=" ;")		// trim off leading white and isolated semi-colons
	Wave wN = str2vec(inStr)						// try to interpret
	Variable NwN=numpnts(wN)

	if (NwN<Nreq)											// too few values, add separators
		Variable i
		i = strlen(inStr)-1
		inStr = SelectString(char2num(inStr[i])==semicolon, inStr, inStr[0,i-1])	// remove trailing ';'
		i = strlen(inStr)-1
		inStr = SelectString(char2num(inStr[i])==period, inStr, inStr[0,i-1])		// remove trailing '.'
		String ch											// one character
		for (i=0;i<(strlen(inStr)-1);i+=1)
			ch = inStr[i]
			if (strsearch(skips,ch,0)>=0)				// a "." or "-"
				continue
			elseif (strsearch(digits,ch,0)>=0)		// a digit type of character
				i += 1
				inStr[i] = ";"
			else
				return $""
			endif
		endfor
		do
			inStr = ReplaceString(";;",inStr,";")	// remove all double ';;' --> ';'
		while(strsearch(inStr,";;",0)>=0)
		Wave wN = str2vec(inStr)							// try to interpret again
	endif

	if (numpnts(wN)==Nreq && numtype(sum(wN))==0)
		return wN
	endif
	return $""
End
//
//	Function test_minStr2Vec(str,N)
//		String str
//		Variable N
//		Wave vint = minStr2Vec(str,N)
//			if (WaveExists(vint))
//			printf "'%s' --> %s\r",str,vec2str(vint)
//		else
//			printf "'%s' --> ****** ERROR *****\r",str
//		endif
//	End



ThreadSafe Function/T ReplaceCharacters(chars,inStr,replacement)
	// replace all occurances of a character in chars[] with replacement
	String chars			// look for each of the characters in chars
	String inStr			// input string
	String replacement	// replace every occurance of chars[i] with this

	Variable i, N=strlen(chars)
	for (i=0;i<N;i+=1)
		inStr = ReplaceString(chars[i],inStr,replacement)	
	endfor
	return inStr
End


//  ====================================================================================  //
//  ============================== Start of Wave Printing ==============================  //

ThreadSafe Function printWave(w1,[name,brief,fmt,zeroThresh])	// print a wave (vector or matrix) to history
	Wave w1
	String name									// optional user supplied name to use
	Variable brief								// print in briefer form
	String fmt									// optional format for on number
	Variable zeroThresh						// |values| < zeroThresh show as a "0", in many vectors 1e-15 is really zero
	if (ParamIsDefault(name))
		name = NameOfWave(w1)
	endif
	brief = ParamIsDefault(brief) ? 0 : !(!brief)
	fmt = SelectString(ParamIsDefault(fmt),fmt,"")
	zeroThresh = ParamIsDefault(zeroThresh) || numtype(zeroThresh) || zeroThresh<=0 ? NaN : zeroThresh
	if (!WaveExists(w1))
		print "in 'printWave', wave does not exist"
		return 1
	endif

	Variable N=numpnts(w1), N0=DimSize(w1,0), N1=DimSize(w1,1), N2=DimSize(w1,2), N3=DimSize(w1,3)
	Variable dim=4
	dim -= N0<=1 ? 1 : 0
	dim -= N1<=1 ? 1 : 0
	dim -= N2<=1 ? 1 : 0
	dim -= N3<=1 ? 1 : 0

	Variable Nx=N0, Ny=N1
	if (N0<=1 && N1<=1)
		Nx = N2	;	Ny = N3
	elseif (N0<=1 && N2<=1)
		Nx = N1	;	Ny = N3
	elseif (N0<=1 && N3<=1)
		Nx = N1	;	Ny = N2
	elseif (N1<=1 && N2<=1)
		Nx = N0	;	Ny = N3
	elseif (N1<=1 && N3<=1)
		Nx = N0	;	Ny = N2
	elseif (N2<=1 && N2<=1)
		Nx = N0	;	Ny = N2
	endif

	String dimStr=""
	if (dim==1)						// a vector either row, column, layer, or beam
		printvec(w1,name=name,fmt=fmt,zeroThresh=zeroThresh)
	elseif (dim==2)					// a 2D matrix
		Duplicate/FREE w1, w2D
		if (N2<=1 && N3<=0)		// for simple 2-d matrix, do nothing
		elseif (N0==1 && N3<=1)	// a 1xN1xN2 matrix
			Redimension/N=(Nx,Ny) w2D
			w2D = w1[0][p][q]
			sprintf dimStr, "[][%d][%d]",Nx,Ny
		elseif (N1==1 && N3<=1)	// a N0x1xN2 matrix
			Redimension/N=(Nx,Ny) w2D
			w2D = w1[p][0][q]
			sprintf dimStr, "[%d][][%d]",Nx,Ny
		elseif (N0==1 && N1==1)	// a 1x1xN2xN3 matrix
			Redimension/N=(Nx,Ny) w2D
			w2D = w1[0][0][p][q]
			sprintf dimStr, "[][][%d][%d]",Nx,Ny
		elseif (N0==1 && N2==1)	// a 1xN1x1xN3 matrix
			Redimension/N=(Nx,Ny) w2D
			w2D = w1[0][p][0][q]
			sprintf dimStr, "[][%d][][%d]",Nx,Ny
		elseif (N1==1 && N2==1)	// a N0x1x1xN3 matrix
			Redimension/N=(Nx,Ny) w2D
			w2D = w1[p][0][0][q]
			sprintf dimStr, "[%d][][][%d]",Nx,Ny
		endif
		return printmat(w2D,name=name+dimStr,brief=brief,fmt=fmt,zeroThresh=zeroThresh)
	else
		print "cannot yet handle dimensions 3 or 4"
	endif
	return 0
End
//
ThreadSafe Static Function printvec(w1,[name,fmt,zeroThresh])	// print a vector to screen
	Wave w1
	String name									// optional user supplied name to use
	String fmt									// optional format for on number
	Variable zeroThresh						// |values| < zeroThresh show as a "0", in many vectors 1e-15 is really zero
	if (ParamIsDefault(name))
		name = NameOfWave(w1)
	endif
	fmt = SelectString(ParamIsDefault(fmt),fmt,"")
	zeroThresh = ParamIsDefault(zeroThresh) || numtype(zeroThresh) || zeroThresh<=0 ? NaN : zeroThresh

	if (strlen(name))
		printf "%s = %s\r", name,vec2str(w1, fmt=fmt, zeroThresh=zeroThresh)
	else
		printf "%s\r", vec2str(w1, fmt=fmt, zeroThresh=zeroThresh)
	endif
End
//
ThreadSafe Static Function printmat(m1,[name,brief,fmt,rowMax,zeroThresh])
	Wave m1
	String name									// optional user supplied name to use
	Variable brief								// print in briefer form
	String fmt									// optional format for on number
	Variable rowMax							// maximum number of rows to print
	Variable zeroThresh						// |values| < zeroThresh show as a "0", in many vectors 1e-15 is really zero
	if (ParamIsDefault(name))
		name = NameOfWave(m1)
	endif
	brief = ParamIsDefault(brief) ? 0 : !(!brief)
	fmt = SelectString(ParamIsDefault(fmt),fmt,"")
	rowMax = ParamIsDefault(rowMax) ? 50 : rowMax
	rowMax = (rowMax>0) ? rowMax : 50
	zeroThresh = ParamIsDefault(zeroThresh) || numtype(zeroThresh) || zeroThresh<=0 ? NaN : zeroThresh
	if (DimSize(m1,1)==0 || DimSize(m1,2)!=0)	// for 2-d matrix only
		print "Can only print 2-d matricies with printmat"
		// DoAlert 0, "Can only print 2-d matricies with printmat"
		return 1
	endif

	if (brief && strlen(name))
		printf "%s:\r",name
	endif
	Variable Nrow=DimSize(m1,0), row
	Nrow = min(Nrow,rowMax)
	for (row=0;row<Nrow;row+=1)
		if (WaveType(m1,1) == 2)				// true for text waves
			print printmatOneListText(m1,row, name=name, brief=brief, fmt=fmt)
		elseif (WaveType(m1) %& 0x01)		// true for complex numbers
			print printmatOneListComplex(m1,row, name=name, brief=brief, fmt=fmt, zeroThresh=zeroThresh)
		else
			print printmatOneListReal(m1,row, name=name, brief=brief, fmt=fmt, zeroThresh=zeroThresh)
		endif
	endfor
	if (DimSize(m1,0)>Nrow)
		print "      ."
		print "      ."
		printf "      .\t\t printed only %d of the %d rows\r",Nrow,DimSize(m1,1)
	endif
//	if (DimSize(mw,0)>Nrow || DimSize(mw,1)>Nrow)
//		printf "Only printed part of the (%d x %d) matrix\r",DimSize(mw,0),DimSize(mw,1)
//	endif
	return 0
End
//
ThreadSafe Static Function/T printmatOneListReal(m1,row,[name,brief,fmt,zeroThresh])// print one line for real (not complex) matricies
	Wave m1
	Variable row								// row number (starts with 0)
	String name									// optional user supplied name to use
	Variable brief								// print in briefer form
	String fmt									// optional format for on number
	Variable zeroThresh						// |values| < zeroThresh show as a "0", in many vectors 1e-15 is really zero
	if (ParamIsDefault(name))
		name = NameOfWave(m1)
	endif
	brief = ParamIsDefault(brief) || numtype(brief) ? 0 : !(!brief)
	fmt = SelectString(ParamIsDefault(fmt)|| strlen(fmt)<1,fmt,"%g")
	zeroThresh = ParamIsDefault(zeroThresh) || numtype(zeroThresh) || zeroThresh<=0 ? NaN : zeroThresh
	if (!zeroThresh || numtype(zeroThresh))
		zeroThresh = DefaultZeroThresh(m1)
	endif

	Duplicate/FREE m1, mInternal
	mInternal = abs(m1)<zeroThresh ? 0 : m1
	if (brief)
		fmt = fmt + "    "
	else
		fmt = "%s[%d][%d] = "+fmt+";    "
	endif

	String line="", str
	Variable j, Ncol=DimSize(m1,1)
	for (j=0;j<Ncol;j+=1)
		if (strlen(line)>100)
			line += "  ..."
			break
		elseif (brief)
			sprintf str, fmt,mInternal[row][j]
		else
			sprintf str, fmt,name,row,j,mInternal[row][j]
		endif
		line += str
	endfor
	line = line[0,strlen(line)-4-1]		// strip off trailing 4 spaces
	return line
End
//
ThreadSafe Static Function/T printmatOneListComplex(m1,row,[name,brief,fmt,zeroThresh])// print one line for complex (not real) matricies
	Wave/C m1
	Variable row								// row number (starts with 0)
	String name									// optional user supplied name to use
	Variable brief								// print in briefer form
	String fmt									// optional format for on number
	Variable zeroThresh						// |values| < zeroThresh show as a "0", in many vectors 1e-15 is really zero
	if (ParamIsDefault(name))
		name = NameOfWave(m1)
	endif
	brief = ParamIsDefault(brief) || numtype(brief) ? 0 : !(!brief)
	fmt = SelectString(ParamIsDefault(fmt),fmt,"%g,%g")
	zeroThresh = ParamIsDefault(zeroThresh) || numtype(zeroThresh) || zeroThresh<=0 ? NaN : zeroThresh

	if (!zeroThresh || numtype(zeroThresh))
		zeroThresh = DefaultZeroThresh(m1)
	endif

	Make/N=(DimSize(m1,0))/D/FREE rInternal,iInternal
	rInternal = real(m1)
	iInternal = imag(m1)
	rInternal = abs(rInternal)<zeroThresh ? 0 : rInternal
	iInternal = abs(iInternal)<zeroThresh ? 0 : iInternal

	if (brief)
		fmt = "("+fmt+")    "
	else
		fmt = "%s[%d][%d] = ("+fmt+");    "
	endif

	String line="", str
	Variable j, Ncol=DimSize(m1,1)
	for (j=0;j<Ncol;j+=1)
		if (strlen(line)>100)
			line += "  ..."
			break
		elseif (brief)
			sprintf str, fmt,rInternal[row][j],iInternal[row][j]
		else
			sprintf str, fmt,name,row,j,rInternal[row][j],iInternal[row][j]
		endif
		line += str
	endfor
	line = line[0,strlen(line)-4-1]		// strip off trailing 4 spaces
	return line
End
//
ThreadSafe Static Function/T printmatOneListText(m1,row,[name,brief,fmt])// print one line for real (not complex) matricies
	Wave/T m1
	Variable row									// row number (starts with 0)
	String name									// optional user supplied name to use
	Variable brief								// print in briefer form
	String fmt									// optional format for on number
	if (ParamIsDefault(name))
		name = NameOfWave(m1)
	endif
	brief = ParamIsDefault(brief) || numtype(brief) ? 0 : !(!brief)
	fmt = SelectString(ParamIsDefault(fmt)|| strlen(fmt)<1,fmt,"%s")
	if (brief)
		fmt = fmt + "    "
	else
		fmt = "%s[%d][%d] = "+fmt+";    "
	endif

	String line="", str
	Variable j, Ncol=DimSize(m1,1)
	for (j=0;j<Ncol;j+=1)
		if (strlen(line)>100)
			line += "  ..."
			break
		elseif (brief)
			sprintf str, fmt,m1[row][j]
		else
			sprintf str, fmt,name,row,j,m1[row][j]
		endif
		line += str
	endfor
	line = line[0,strlen(line)-4-1]		// strip off trailing 4 spaces
	return line
End

//  =============================== End of Wave Printing ===============================  //
//  ====================================================================================  //



//  ====================================================================================  //
//  ============================== Start of Square Pixels ==============================  //
Function SetAspectToSquarePixels(gName, [printIt])
	// Used to square up a graph window
	String gName										// name of the graph, use "" for the top graph
	Variable printIt
	printIt = ParamIsDefault(printIt) || numtype(printIt) ? strlen(GetRTStackInfo(2))==0 : printIt

	Variable aspect = GetGraphAspectRatio(gName, printIt=printIt)
	if (aspect<1)
		ModifyGraph/W=$gName height={Aspect,aspect}, width=0
	elseif (aspect>=1)
		ModifyGraph/W=$gName width={Aspect,1/aspect}, height=0
	endif
	return aspect
End


Function GetGraphAspectRatio(gName, [printIt])
	// returns aspect ratio of graph in units
	String gName										// name of the graph, use "" for the top graph
	Variable printIt
	printIt = ParamIsDefault(printIt) || numtype(printIt) ? strlen(GetRTStackInfo(2))==0 : printIt
	if (strlen(gName)<1)
		gName = StringFromList(0,WinList("*",";","WIN:1"))
	endif
	if (WinType(gName)!=1)
		if (printIt)
			print "ERROR -- GetGraphAspectRatio(), '"+gName+"' is not an graph"
		endif
		return NaN										// if no image on graph, do not try to set aspect ratio
	endif

	GetAxis/W=$gName/Q bottom
	if (V_flag)											// if no bottom, try top
		GetAxis/W=$gName/Q top
	endif
	if (V_flag)
		if (printIt)
			print "ERROR --- GetGraphAspectRatio(), unable to get size of vertical axis"
		endif
		return NaN
	endif
	Variable width = abs(V_max-V_min)

	GetAxis/W=$gName/Q left
	if (V_flag)											// if no left, try right
		GetAxis/W=$gName/Q right
	endif
	if (V_flag)
		if (printIt)
			print "ERROR -- GetGraphAspectRatio(), unable to get size of horizontal axis"
		endif
		return NaN
	endif
	Variable height = abs(V_max-V_min)

	Variable aspect = height / width
	if (numtype(aspect) || aspect<=0)
		return NaN
	endif
	return aspect
End


Function SquareUpGizmo(gName)
	// Used to square up a graph gizmo
	String gName
	if(exists("NewGizmo")!=4)				// Do nothing if the Gizmo XOP is not available.
		DoAlert 0, "Gizmo XOP must be installed"
		return 1
	elseif (strlen(gName)<1)
#if (IgorVersion()<7)
		Execute "GetGizmo gizmoName"
		gName = StrVarOrDefault("S_GizmoName","")
		KillStrings/Z S_GizmoName
#else
		GetGizmo gizmoName
		gName = S_GizmoName
#endif
	endif
	if (strlen(gName)<1)
		return 1
	endif

#if (IgorVersion()<7)
	Execute "GetGizmo winPixels"			// get window position & size
	NVAR V_left=V_left, V_top=V_top, V_bottom=V_bottom
#else
	ModifyGizmo keepPlotSquare=1
	GetGizmo winPixels						// get window position & size
#endif
	Variable height=V_bottom-V_top, top=max(44,V_top)
	MoveWindow/W=$gName V_left, top, V_left+height, V_top+height
	KillVariables/Z V_left, V_right, V_top, V_bottom
End

//  =============================== End of Square Pixels ===============================  //
//  ====================================================================================  //



#if (IgorVersion()>7)
Function SetGizmoZoom(zoom)
	// Interactive Set zoomFactor for a Gizmo
	Variable zoom
	String win = StringFromList(0,WinList("*",";","WIN:"+num2istr(GIZMO_WIN_BIT)))	// top Gizmo
	if (strlen(win)<1)
		return NaN
	endif
	if (numtype(zoom) || zoom<=0)
		zoom = GetGizmoZoom(win)
		Prompt zoom, "zoom must be >0 and finite, (0.8 is nice)"
		DoPrompt "Gizmo Zoom",zoom
		if (V_flag)
			return NaN
		endif
	endif
	if (numtype(zoom) || zoom<=0)
		printf "ERROR -- SetGizmoZoom(), Invalid zoom = %g,  must be positive finite\r",zoom
	endif

	ModifyGizmo/N=$win zoomFactor = zoom

	if (exists("ResetScaleBarLength")==6)				// and update scale bars if GizmoUtility.ipf is loaded
		String keyList=GetUserData(win,"","ScaleBar")	// if there is a scaleBar, update it too.
		if (strlen(keyList)>1)
			keyList = ReplaceNumberByKey("scaleFactor",keyList,zoom,"=")
			SetWindow $win userdata(ScaleBar)=keyList
			FUNCREF ResetScaleBarLengthProto func = $"ResetScaleBarLength"
			func()
		endif
	endif
	return zoom
End
//
Function ResetScaleBarLengthProto([units,scaleFactor,font,fontSize,win])
	String units
	Variable scaleFactor
	String font
	Variable fontSize
	String win
	return NaN
End

Static Function GetGizmoZoom(win)
	String win
	String str = WinRecreation(win,0)
	if (strlen(str)<1)
		return NaN
	endif
	str = ReplaceString(" ",str,"")
	str = ReplaceString("\t",str,"")
	Variable zoom = NumberByKey("ModifyGizmozoomFactor",str,"=","\r")
	zoom = numtype(zoom) || zoom<=0 ? 1 : zoom
	return zoom
End

//Function SetGizmoZoom(zoom)
//	// Interactive Set zoomFactor for a Gizmo
//	Variable zoom
//	String win = StringFromList(0,WinList("*",";","WIN:"+num2istr(GIZMO_WIN_BIT)))	// top Gizmo
//	if (strlen(win)<1)
//		return NaN
//	endif
//	if (numtype(zoom) || zoom<=0)
//		String str = WinRecreation(win,0)
//		str = ReplaceString(" ",str,"")
//		str = ReplaceString("\t",str,"")
//		zoom = NumberByKey("ModifyGizmozoomFactor",str,"=","\r")
//		zoom = numtype(zoom) || zoom<=0 ? 1 : zoom
//		Prompt zoom, "zoom must be >0 and finite, (0.8 is nice)"
//		DoPrompt "Gizmo Zoom",zoom
//		if (V_flag)
//			return NaN
//		endif
//	endif
//	if (numtype(zoom) || zoom<=0)
//		printf "ERROR -- SetGizmoZoom(), Invalid zoom = %g,  must be positive finite\r",zoom
//	endif
//	ModifyGizmo/N=$win zoomFactor = zoom
//	return zoom
//End
#else
Function SetGizmoZoom(zoom)			// There is no equivalent for Igor 6
	Variable zoom
	return NaN
End
#endif



//  ====================================================================================  //
//  ============================= Start of Unicode Letters =============================  //

Function/T Name2SymbolCharacter(name)				// returns the character equivalent for the Symbol font
	String name												// name of a Symbol character, e.g. "theta", or "aleph"

	String symNames="alpha:A;beta:B;gamma:G;delta:D;epsilon:E;zeta:Z;eta:H;theta:Q;iota:I;kappa:K;lambda:L;"
	symNames += "mu:M;nu:N;xi:X;xsi:X;omicron:O;pi:P;rho:R;sigma:S;tau:T;upsilon:U;phi:F;chi:C;psi:Y;omega:W;"
	symNames += "forAll:34;Exists:36;EqualTilde:64;therefore:92;perp:94;inf:165;club:167;diamond:168;heart169;spade:170;"
	symNames += "leftRight:171;left:172;up:173;right:174;down:175;"
	symNames += "ex:180;proportional:181;partialDif:182;bullet:183;div:184;NotEqual:185;"
	symNames += "equivalent:186;approx:187;3dot:188;return:191;"
	symNames += "aleph:192;Imag:193;Real:194;InnerProd:196;InnerSum:197;circleCross:198;"
	symNames += "Intersection:199;Union:200;subset:201;superset:204;"
	symNames += "angle:208;del:209;registered:210;copyright:211;TM:212;Product:213;root:214;"
	symNames += "leftRightDouble:219;leftDouble:220;upDouble:221;rightDouble:222;downDouble:223;"
	symNames += "bra:225;bar:231;ket:241;Sum:229;"

	String ch = StringByKey(name,symNames)
	Variable i = str2num(ch)

	if (i>0)
		ch = num2char(i)
	else
		ch = SelectString(char2num(name) < 97 , LowerStr(ch), UpperStr(ch))
	endif
	return ch
End


#if (IgorVersion()<7)
Function/T Letter2SymbolOrUnicode(letter)		// Igor 6 does not support unicode
	String letter
	String ch = Name2SymbolCharacter(letter)	
	if (strlen(ch))
		return "\\F'Symbol'"+ch+"\\F]0"
	endif
	return ""
End
#else
Function/T Letter2SymbolOrUnicode(letter)
	String letter
	return LetterName2Unicode(letter)
End
#endif


#if (IgorVersion()<7)
Function/T LetterName2Unicode(letter)		// Igor 6 does not support unicode
	String letter

	if (cmpstr(letter,"Delta",0)==0)
		letter = num2char(198)
	elseif (cmpstr(letter,"mu",1)==0)
		letter = num2char(181)
	elseif (cmpstr(letter,"Sigma",1)==0)
		letter = num2char(183)
	elseif (cmpstr(letter,"Pi",1)==0)
		letter = num2char(184)
	elseif (cmpstr(letter,"pi",1)==0)
		letter = num2char(185)
	elseif (cmpstr(letter,"Omega",1)==0)
		letter = num2char(189)
	endif
	return letter
End
#else
Function/T LetterName2Unicode(letter)
	// returns a unicode letter from the name "letter", e.g.
	//		LetterName2Unicode("Q") --> "Q"
	//		LetterName2Unicode("Delta") --> Delta
	//		LetterName2Unicode("alpha") --> alpha
	String letter
	if (strlen(letter)<2)		// just a roman letter, not a  "name"
		return letter
	endif
	String greek = Greek2Unicode(letter)
	if (strlen(greek)>0)
		return greek	
	endif
	String hebrew = Hebrew2Unicode(letter)
	if (strlen(hebrew)>0)
		return hebrew	
	endif
	return ""			// failed to find anything
End
//
Static Function/T Greek2Unicode(letter)
	String letter
	Variable upper=0, lower=0, i=char2num(letter)
	upper = (65<=i && i<=90)
	lower = (97<=i && i<=122)
	letter = LowerStr(letter)
	letter = ReplaceString(" ",letter,"")
	letter = ReplaceString("lamda",letter,"lambda")		// optional spelling
	String allGreekLetters = "alpha;beta;gamma;delta;epsilon;zeta;eta;theta;iota;kappa;lambda;mu;nu;xi;omicron;pi;rho;finalsigma;sigma;tau;upsilo;phi;chi;psi;omega;"
	i = WhichListItem(letter,allGreekLetters)
	if (i<0 || (upper+lower)<1)
		return ""
	endif
	i += lower ? 0x03B1 : 0x0391		// offset to Greek section (lowercase starts at 0x03B1 = 0x0391+32)
	return num2char(i)
End
//
//	for a vav with a dot in middle, use:  LetterName2Unicode("vav")+LetterName2Unicode("Dagesh")
Static Function/T Hebrew2Unicode(letter)
	String letter							//	http://www.i18nguy.com/unicode/hebrew.html
	Variable i=char2num(letter)
	letter = LowerStr(letter)
	letter = ReplaceString(" ",letter,"")

	letter = ReplaceString("Mapiq",letter,"Dagesh")		// optional spellings
	letter = ReplaceString("Shuruq",letter,"Dagesh")
	Variable addSinDot = 0
	if (StringMatch(letter,"sin"))
		letter = "shin"
		addSinDot = 1
	endif

	String allHebrewLetters = "alef;bet;gimel;dalet;he;vav;zayin;het;tet;yod;finalkaf;kaf;lamed;finalmem;mem;finalnun;nun;samekh;ayin;finalpe;pe;finaltsadi;tsadi;qof;resh;shin;tav;"
	i = WhichListItem(letter,allHebrewLetters)
	if (i>=0)
		i += 0x05D0							// offset to start of Hebrew section
	elseif (StringMatch(letter,"Dagesh"))
		i = 0x05BC
	elseif (i<0)
		return ""
	endif

	String out = num2char(i)
	if (addSinDot)
		out += num2char(0x05C2)
//	elseif (StringMatch(letter,"shin"))
//		out += num2char(0x05C1)
	endif

	return out
End
#endif

//  ============================== End of Unicode Letters ==============================  //
//  ====================================================================================  //



ThreadSafe Function/T num2Ordinal(n)
	// converts integer n to ordinal string, e.g. 1 --> "1st"
	// if n is not an integer then just append "th"
	Variable n		// an integer
	Variable i = mod(abs(n),10)
	i = (abs(n) == limit(abs(n),10,19)) ? 0 : i		// because 11, 12, 13 just get "th", but 31 --> "31st"
	String ending=StringFromList(i,"th;st;nd;rd;")
	ending = SelectString(strlen(ending),"th",ending)
	if (numtype(n))				// NaN or Inf
		ending = ""
	elseif (abs(round(n)-n)>1e-5)
		ending = "th"				// not an integer
	else
		n = round(n)
	endif
	return num2str(n)+ending
End


ThreadSafe Function/T ChangeStrEnding(oldEnd, inStr, newEnd)
	// if inStr ends in oldEnd, replace oldEnd with newEnd
	// if oldEnd =="", always add newEnd
	// if newEnd=="", removes oldEnd
	String oldEnd					// existing ending to replace
	String inStr					// given string
	String newEnd					// replace oldEnd with this
	if (StringMatch(inStr,"*"+oldEnd))
		inStr = RemoveEnding(inStr,oldEnd)+newEnd
	endif
	return inStr
End


//  ====================================================================================  //
//  ============================= Start of Unit Conversions ============================  //
Static Constant c = 299792458						// exact speed of light (m/s)
Static Constant e_C = 1.6021766208e-19			// Charge on electron (C)
Static Constant hbar = 1.054571800e-34			// reduced Planck Constant, h/2pi (J s)
Static Constant kB = 8.6173303e-5					// Boltzmann constant (eV/K), updated March 2017
Static Constant GN = 6.67408e-11					// Newton Gravity Constant (m^3 kg^-1 s^-2)
Static Constant inch = 0.0254						// length of inch (m)
Static Constant kgPerPound = 0.45359237			// 1 pound = 0.45359237 kgm [definition of pound]
Static Constant gStd = 9.80665						// std acceleration of gravity (m * s^-2)
Static Constant tropicalYear = 31556925.216	// = 365.24219 * 24*3600, seconds in a tropical year (NOT sidereal), there are 365.24219 days in 1 tropical year
Static Constant julianYear = 31557600				// = 365.25 * 24*3600, seconds in a Julian, there are 365.25 days in 1 Julian year
Static Constant AstronomicalUnit = 149597870700		// IAU 2009,2012
Static Constant stdAtmosphere = 101325			// standard atmospheric pressure (Pascal)
Static Constant inH2O = 249.082						// pressure of 1 inch of water (Pascal)


ThreadSafe Function SIprefix2factor(prefix)
	String prefix
	// return the product of the prefixes, e.g. "mp" returns 1e-15
	// white space (anything <= a space) is ignored.
	// for bad prefix values, (e.g. "v" or "3") returns NaN
	// Note, internally I use "o" instead of "�" since the NumberByKey() routine does not work with a key="�"
	// also note that all prefixes are case sensitive except for "H" and "K", which can be either upper or lower.
	//		deci	= "d" = 0.1
	//		centi	= "c" = 0.01
	//		milli	= "m" = 1e-3
	//		micro	= "�" = 1e-6
	//		nano	= "n" = 1e-9
	//		pico	= "p" = 1e-12
	//		femto	= "f" = 1e-15
	//		atto	= "a" = 1e-18
	//		zepto	= "z" = 1e-21
	//		yocto	= "y" = 1e-24
	//
	//		hecto	= "H" = 100  (or "h")
	//		kilo	= "K" = 1e3  (or "k")
	//		Mega	= "M" = 1e6
	//		Giga	= "G" = 1e9
	//		Tera	= "T" = 1e12
	//		Peta	= "P" = 1e15
	//		Exa	= "E" = 1e18
	//		Zeta	= "Z" = 1e21
	//		Yotta	= "Y" = 1e24

	if (strsearch(prefix,"o",0)>=0)				// need to check for "o", which is also invalid
		return NaN
	endif
	prefix = ReplaceString(Gmu,prefix,"o")	// NumberByKey() routine does not work with a key="�", so use internally use "o" instead

	String keyVals="d:0.1;c:0.01;m:1e-3;o:1e-6;n:1e-9;p:1e-12;f:1e-15;a:1e-18;z:1e-21;y:1e-24;"
	keyVals += "h:100;H:100;k:1e3;K:1e3;M:1e6;G:1e9;T:1e12;P:1e15;E:1e18;Z:1e21;Y:1e24;"

	Variable i, value
	String ch
	for (i=0,value=1; i<strlen(prefix); i+=1)
		ch = prefix[i]
		if (char2num(ch)>32)
			value *= NumberByKey(ch,keyVals,":",";",1)
		endif
	endfor
	return value
End
//	Function test_SIprefix2factor()
//		String chars=" dcm�npfazyHhKkMGTPEZYv"
//		Variable i
//		for (i=0;i<strlen(chars);i+=1)
//			printf " '%s'  %g\r",chars[i],SIprefix2factor(chars[i])
//		endfor
//		print " "
//		print "following should be 1's"
//		printf "%g, %g, %g, %g, %g, ", SIprefix2factor("cH"), SIprefix2factor("mk"), SIprefix2factor("�M"), SIprefix2factor("nG"), SIprefix2factor("pT")
//		printf "%g, %g, %g, %g\r", SIprefix2factor("fP"), SIprefix2factor("aE"), SIprefix2factor("zZ"), SIprefix2factor("yY")
//	End


ThreadSafe Function ConvertUnits2meters(unit,[defaultLen])
	// returns conversion factor from unit to meters
	String unit
	Variable defaultLen
	defaultLen = ParamIsDefault(defaultLen) ? NaN : defaultLen
	//
	//	meter, metre, m		1 m
	//	inch, inches, in		0.0254 m
	//	foot, feet, ft			12 inches
	//	yard, yd					36 inches
	//	mile, mi					5280 feet
	//	nauticalmile			1852 m
	//	Angstrom, Ang,ARING	1e-10 m
	//	micron, micrometer	1e-6 m
	//	RackUnit, Rack, U		1.75 inch
	//	parsec, pc				1 parsec = AU * (180*3600)/pi = 3.085677581e16 (m), IAU definition
	//	lightYear, ly			9460730472580800 m  =  c * julianYear (this is the IAU definition)
	//	astronomicalunit, au	149597870700 m
	//	BohrRadius, ao, a0	0.52917721092e-10 m
	//	fermi, fm				1e-15 m == 1 fm
	//	mil						0.001 inch
	//	CuXunit, CuXU			Cu X-unit, 1.00207697e-13 m
	//	MoXunit, MoXU			Mo X-unit, 1.00209952e-13 m
	//	Xunit, XU				X-unit (just average of Mo & Cu), 1.002088e-13 m
	//	Planck, PlanckLength sqrt(hbar*GN/c^3)  =  1.616199e-35 m
	//	fathom					6 feet
	//	chain						66 feet
	//	link						661/100 feet (0.01 of chain)
	//	rod						16.5 feet
	//	league					3 miles
	//	furlong					660 feet
	//	cubit						0.525 m (a rough number)
	//	point						1/72 inch
	//	pica						1/72 foot
	//	Li							Chinese mile is 500m

	unit = ReplaceString(" ",unit,"")			// no spaces

	// check for powers  "^N"
	Variable ipow=strsearch(unit,"^",0), power=1
	if (ipow>=0)										// found a power
		power = str2num(unit[ipow+1,Inf])		// number after "^"
		unit = unit[0,ipow-1]						// string before "^"
	endif

	// fix spellings
	unit = ChangeStrEnding("metre",unit,"meter")	// British spelling
	unit = ChangeStrEnding("feet",unit,"foot")
	unit = ChangeStrEnding("inches",unit,"inch")
	unit = ChangeStrEnding("fermi",unit,"fm")
	unit = ChangeStrEnding("√Ö",unit,ARING)		// funny encodings of Angstrom symbol
	unit = ChangeStrEnding("\201",unit,ARING)		// Mac opt-shift-A
	unit = ChangeStrEnding("\xC3\x85",unit,ARING)	// UTF8
	unit = ChangeStrEnding("Ang",unit,ARING)		// lots of ways to write Angstrom
	unit = ChangeStrEnding("Angstrom",unit,ARING)
	unit = ChangeStrEnding("micrometer",unit,Gmu+"m")		// Gmu is Greek mu
	unit = ChangeStrEnding("micron",unit,Gmu+"m")
	unit = ChangeStrEnding("micro",unit,Gmu)
	unit = ChangeStrEnding("micrometer",unit,Gmu+"m")
	unit = ChangeStrEnding("micron",unit,Gmu+"m")
	unit = ChangeStrEnding("micro",unit,Gmu)
	unit = ChangeStrEnding("RackUnit",unit,"U")
	unit = ChangeStrEnding("Rack",unit,"U")
	unit = RemoveEnding(unit,"s")				// remove any trailing "s"

	String prefix
	Variable value=NaN, i = max(0,strlen(unit)-1)
	if (strsearch(unit,"m",i)==i)				// ends in 'm', means meters
		value = 1
		prefix = unit[0,strlen(unit)-2]
	elseif(StringMatch(unit,"*"+ARING))		// the Angstrom
		value = 1e-10
		prefix = unit[0,strlen(unit)-2]
	elseif(StringMatch(unit,"*CuXunit") || StringMatch(unit,"*CuXU"))	// the Cu X-unit
		value = 1.00207697e-13 
		i = StringMatch(unit,"*CuXU") ? 5 : 8
		prefix = unit[0,strlen(unit)-8]
	elseif(StringMatch(unit,"*MoXunit") || StringMatch(unit,"*MoXU")) // the Mo X-unit
		value = 1.00209952e-13 
		i = StringMatch(unit,"*MoXU") ? 5 : 8
		prefix = unit[0,strlen(unit)-8]
	elseif(StringMatch(unit,"*Xunit") || StringMatch(unit,"*XU"))		// the X-unit (just average of Mo & Cu)
		value = 1.002088e-13
		i = StringMatch(unit,"*XU") ? 3 : 6
		prefix = unit[0,strlen(unit)-i]
	elseif(StringMatch(unit,"*pc") || StringMatch(unit,"*parsec"))
		value = AstronomicalUnit * 180*3600 / pi			// = 3.085677581e16 (m)
		i = StringMatch(unit,"*parsec") ? 7 : 3
		prefix = unit[0,strlen(unit)-i]
	elseif(StringMatch(unit,"*ly") || StringMatch(unit,"*lightYear"))
		value = c * julianYear				// = 9460730472580800 m
		i = StringMatch(unit,"*lightYear") ? 10 : 3
		prefix = unit[0,strlen(unit)-i]
	elseif(StringMatch(unit,"*au") || StringMatch(unit,"*astronomicalunit"))
		value = AstronomicalUnit
		i = StringMatch(unit,"*au") ? 3 : 17
		prefix = unit[0,strlen(unit)-i]
	elseif(StringMatch(unit,"*a0") || StringMatch(unit,"*ao") || StringMatch(unit,"*BohrRadiu"))
		value = 0.52917721092e-10
		i = StringMatch(unit,"*BohrRadiu") ? 10 : 3		// the "s" got trimmed off!
		prefix = unit[0,strlen(unit)-i]
	elseif(StringMatch(unit,"*Planck") || StringMatch(unit,"*PlanckLength"))
		value = sqrt(hbar*GN/c^3)
		i = StringMatch(unit,"*Planck") ? 6 : 13
		prefix = unit[0,strlen(unit)-i]
	elseif(StringMatch(unit,"*in") || StringMatch(unit,"*inch"))	 	// a few english units just for fun
		value = inch
		i = StringMatch(unit,"*inch") ? 5 : 3
		prefix = unit[0,strlen(unit)-i]
	elseif(StringMatch(unit,"*ft") || StringMatch(unit,"*foot") || StringMatch(unit,"*feet"))
		value = 12*inch
		i = StringMatch(unit,"*ft") ? 3 : 5
		prefix = unit[0,strlen(unit)-i]
	elseif(StringMatch(unit,"*nauticalmile"))
		value = 1852
		prefix = unit[0,strlen(unit)-13]
	elseif(StringMatch(unit,"*mi") || StringMatch(unit,"*mile"))
		value = 5280 * 12*inch
		i = StringMatch(unit,"*mile") ? 5 : 3
		prefix = unit[0,strlen(unit)-i]
	elseif(StringMatch(unit,"*yd") || StringMatch(unit,"*yard"))
		value = 3 * 12*inch
		i = StringMatch(unit,"*yard") ? 5 : 3
		prefix = unit[0,strlen(unit)-i]
	elseif(StringMatch(unit,"*mil"))			// 0.001 inch
		value = 0.001 * inch
		prefix = unit[0,strlen(unit)-4]
	elseif(StringMatch(unit,"*fathom"))		// 6 feet
		value = 6 * 12*inch
		prefix = unit[0,strlen(unit)-7]
	elseif(StringMatch(unit,"*chain"))		// 66 feet
		value = 66 * 12*inch
		prefix = unit[0,strlen(unit)-6]
	elseif(StringMatch(unit,"*link"))			// 66/100 feet
		value = 66 * 12*inch /100
		prefix = unit[0,strlen(unit)-5]
	elseif(StringMatch(unit,"*rod"))			// 16.5 feet
		value = 16.5 * 12*inch
		prefix = unit[0,strlen(unit)-4]
	elseif(StringMatch(unit,"*league"))		// 3 miles
		value = 3 * 5280 * 12*inch
		prefix = unit[0,strlen(unit)-7]
	elseif(StringMatch(unit,"*furlong"))		// 660 feet
		value = 660 * 12*inch
		prefix = unit[0,strlen(unit)-8]
	elseif(StringMatch(unit,"*cubit"))		// a rough number
		value = 0.525
		prefix = unit[0,strlen(unit)-6]
	elseif(StringMatch(unit,"*point"))		// 1/72 inch
		value = inch / 72
		prefix = unit[0,strlen(unit)-6]
	elseif(StringMatch(unit,"*pica"))			// 1/72 foot
		value = 12*inch / 72
		prefix = unit[0,strlen(unit)-5]
	elseif(StringMatch(unit,"*Li"))			// Chinese mile is 500m
		value = 500
		prefix = unit[0,strlen(unit)-3]
	elseif(StringMatch(unit,"*U"))				// 1.75 unch
		value = 1.75 * inch
		prefix = unit[0,strlen(unit)-2]
	else
		return defaultLen								// cannot find base value
	endif

	value *= SIprefix2factor(prefix)
	return (power==1) ? value : (value ^ power)
End


Function/T ConvertNumWithUnitsLength(in,outUnit,[rnd,space])
	// converts string like "2in" to "0.0508 m"
	// ConvertNumWithUnitsLength("2 in","cm") --> "5.08 cm"
	// if you know the units and just want to conver a number, use ConvertUnits2meters()
	String in				// input string e.g. "2 in" or "1nm^-2"...
	String outUnit			// desired final unit, e.g. "cm" or "Angstrom^-2"
	Variable rnd			// if rnd True, use same number of places in final as input (default=False)
	Variable space			// flag, if True put space between number and unit (default=True)
	rnd = ParamIsDefault(rnd) ? NaN : rnd
	rnd = numtype(rnd) ? 0 : !(!rnd)
	space = ParamIsDefault(space) ? 1 : !(!space)

	String list = SplitNumAndUnits(in)
	Variable val = str2num(StringFromList(0,list))
	String inUnit = StringFromList(1,list)
	if (strlen(inUnit)<1 || numtype(val))
		return ""
	endif
	String powIn,powOut				// check that units have same power
	Variable i = strsearch(inUnit,"^",strlen(inUnit)-1,1)
	i = i<0 ? Inf : i
	powIn = inUnit[i,Inf]
	i = strsearch(outUnit,"^",strlen(outUnit)-1,1)
	i = i<0 ? Inf : i
	powOut = outUnit[i,Inf]
	if (cmpstr(powIn, powOut))
		return ""
	endif

	Variable places = rnd ? placesOfPrecision(val)+1 : 16
	places = limit(places,1,16)
	val *= ConvertUnits2meters(inUnit)/ConvertUnits2meters(outUnit)

	String out, sspace=SelectString(space,""," ")
	String fmt="%."+num2istr(places)+"g"+sspace+"%s"
	sprintf out, fmt,val,outUnit
	return out
End


Function/T SplitNumAndUnits(in)
	// take string such as "-9.3e2mi" --> "-9.3e2;mi",  or "2 nm^-1" --> "2;nm^-1"
	// this is useful with ConvertNumWithUnitsLength() above
	String in

	in = ReplaceString(" ",in,"")
	in = ReplaceString("\t",in,"")
	Variable i,N=strlen(in)
	i = strsearch(in,"^",N-1,1)-1
	i = i<0 ? N-1 : i
	for ( ;i>=0;i-=1)
		if (strsearch("01234567890.,",in[i],0)>=0)
			i += 1
			break
		endif
	endfor
	String unit = in[i,N-1]
	String num = in[0,i-1]
	num = RemoveEnding(num,".")
	num = RemoveEnding(num,",")
	return num+";"+unit
End


ThreadSafe Function ConvertTemperatureUnits(Tin,unitIn,unitOut)	// returns Temperature in units of (unitOut)
	// ConvertTemperatureUnits(400,"K","F") --> "260.33"
	// work for Celsius, Kelvin, Fahrenheit, Rankine, & Planck units
	// for Planck Temperature, see: http://en.wikipedia.org/wiki/Planck_temperature
	String unitIn, unitOut			// input and output units, NO defaults allowed
	Variable Tin						// input temperature in units of (unitIn)

	unitIn = ReplaceString("Temperature",unitIn,"")
	unitOut = ReplaceString("Temperature",unitOut,"")
	unitIn = ReplaceString("Temp",unitIn,"")
	unitOut = ReplaceString("Temp",unitOut,"")

	unitIn = ReplaceString("Kelvin",unitIn,"K")	// use single letter abreviations
	unitOut = ReplaceString("Kelvin",unitOut,"K")
	unitIn = ReplaceString("Celsius",unitIn,"C")
	unitOut = ReplaceString("Celsius",unitOut,"C")
	unitIn = ReplaceString("Fahrenheit",unitIn,"F")
	unitOut = ReplaceString("Fahrenheit",unitOut,"F")
	unitIn = ReplaceString("Rankine",unitIn,"R")
	unitOut = ReplaceString("Rankine",unitOut,"R")
	unitIn = ReplaceString("Planck",unitIn,"P")	// Planck Temperature
	unitOut = ReplaceString("Planck",unitOut,"P")

	unitIn = ReplaceString(" ",unitIn,"")				// no spaces
	unitIn = ReplaceString("\241",unitIn,"")		// no degree signs, Mac opt-shift-9
	unitIn = ReplaceString("\xC2\xB0",unitIn,"")	// no degree signs, UTF8
	unitIn = RemoveEnding(unitIn,"s")					// and no trailing 's'

	unitOut = ReplaceString(" ",unitOut,"")			// no spaces
	unitOut = ReplaceString("\241",unitOut,"")		// no degree signs, Mac opt-shift-9
	unitOut = ReplaceString("\xC2\xB0",unitOut,"")// no degree signs, UTF8
	unitOut = RemoveEnding(unitOut,"s")				// and no trailing 's'

	if (strlen(unitIn)<1 || strlen(unitOut)<1)
		return NaN
	endif

	Variable T_Planck = sqrt(hbar * c^5 / (GN * (kB*e_C)^2))	// Planck Temperature = 1.416833e32 (K)
	Variable Tout, Kelvin, n
	if (strlen(unitIn)>1)
		n = strlen(unitIn)
		Tin *= SIprefix2factor(unitIn[0,n-2])
		unitIn = unitIn[n-1]								// the last character
	endif
	strswitch(unitIn)											// first convert Tin to Kelvin
		case "F":
			Kelvin = (Tin - 32) * 5/9 + 273.15
			break
		case "R":
			Kelvin = Tin * 5/9
			break
		case "C":
			Kelvin = Tin + 273.15
			break
		case "K":
			Kelvin = Tin
			break
		case "P":
			Kelvin = Tin / T_Planck
			break
		default:
			return NaN
		endswitch

	Variable factor=1
	if (strlen(unitOut)>1)
		n = strlen(unitOut)
		factor = SIprefix2factor(unitOut[0,n-2])
		unitOut = unitOut[n-1]							// the last character
	endif
	strswitch(unitOut)									// convert Kelvin to unitOut
		case "F":
			Tout = (Kelvin-273.15)*9/5 + 32
			break
		case "R":
			Tout = Kelvin * 9/5
			break
		case "C":
			Tout = Kelvin - 273.15
			break
		case "K":
			Tout = Kelvin
			break
		case "P":
			Tout = Kelvin / T_Planck
			break
		default:
			return NaN
		endswitch
	Tout /= factor

	return Tout
End


ThreadSafe Function ConvertUnits2seconds(unit,[defaultSeconds])
	// returns conversion factor from unit to seconds (NO plural units allowed)
	// note, dy is deci-year, not day
	String unit
	Variable defaultSeconds
	defaultSeconds = ParamIsDefault(defaultSeconds) ? NaN : defaultSeconds
	//
	//	second, sec, s		1
	//	minute, min, m		60 s
	//	hour, hr, h			3600 s
	//	day, d				24 hour
	//	week, wk				7 days
	//	fortnight				14 days
	//	lunar month, moon, lune	29.530588 days
	//	year, yr, y			tropical year = 365.24219 days
	//	olympiad				4 years
	//	lustrum				5 years
	//	indiction				15 years
	//	decade				10 years
	//	century				100 years
	//	millennium			1000 years
	//	jiffy					1 fm/c
	//	shake					1e-8 s
	//	beat					0.001 day = 3.6 s
	//	Planck time			1.616199e-35 / c  =  sqrt(hbar*GN / c^5)
	//	Svedberg				1e-13 s, Abbreviation is "S" which is too close to seconds
	//	galactic year		230e6 tropical years
	//	sidereal day			23.9344699 hour
	//	sidereal year		365.256363004 days
	//	helek					3 + 1/3 seconds
	//	pahar					3 hours

	unit = ReplaceString(" ",unit,"")				// no spaces

	// check for powers  "^N"
	Variable ipow=strsearch(unit,"^",0), power=1
	if (ipow>=0)											// found a power
		power = str2num(unit[ipow+1,Inf])			// number after "^"
		unit = unit[0,ipow-1]							// string before "^"
	endif

	String prefix
	Variable value=NaN, i = max(0,strlen(unit)-1)
	Variable hour = 3600								// seconds in 1 hour
	Variable day = 24*hour								// seconds in 1 day

	if(StringMatch(unit,"*second") || StringMatch(unit,"*sec") || StringMatch(unit,"*s"))				// the second
		i = StringMatch(unit,"*second") ? 7 : 2
		i = StringMatch(unit,"*sec") ? 4 : i
		value = 1
		prefix = unit[0,strlen(unit)-i]
	elseif(StringMatch(unit,"*decade"))			// decade = 10 years
		value = 10 * tropicalYear
		prefix = unit[0,strlen(unit)-7]
	elseif(StringMatch(unit,"*century"))			// century = 100 years
		value = 100 * tropicalYear
		prefix = unit[0,strlen(unit)-8]
	elseif(StringMatch(unit,"*millennium"))		// millennium = 1000 years
		value = 1000 * tropicalYear
		prefix = unit[0,strlen(unit)-11]
	elseif(StringMatch(unit,"*fortnight") )		// fortnight = 2 weeks
		value = 14 * day
		prefix = unit[0,strlen(unit)-10]
	elseif(StringMatch(unit,"*lune") || StringMatch(unit,"*lunarmonth") || StringMatch(unit,"*moon"))	// 1 lunar month
		value = 29.530588 * day
		i = StringMatch(unit,"*lunarmonth") ? 11 : 5
		prefix = unit[0,strlen(unit)-i]
	elseif(StringMatch(unit,"*olympiad"))			// olympiad = 4 years
		value = 4 * tropicalYear
		prefix = unit[0,strlen(unit)-9]
	elseif(StringMatch(unit,"*lustrum"))			// lustrum = 5 years
		value = 5 * tropicalYear
		prefix = unit[0,strlen(unit)-8]
	elseif(StringMatch(unit,"*indiction"))		// indiction = 5 years
		value = 15 * tropicalYear
		prefix = unit[0,strlen(unit)-10]
	elseif(StringMatch(unit,"*jiffy"))				// time to go one fermi
		value = 1e-15 / c	
		prefix = unit[0,strlen(unit)-6]
	elseif(StringMatch(unit,"*shake"))				// 1e-8 s
		value = 1e-8
		prefix = unit[0,strlen(unit)-6]
	elseif(StringMatch(unit,"*beat"))				// 1/1000 of a day (Swatch decimal time)
		value = 3.6
		prefix = unit[0,strlen(unit)-5]
	elseif(StringMatch(unit,"*Planck") || StringMatch(unit,"*Plancktime"))
		value = sqrt(hbar*GN / c^5)
		i = StringMatch(unit,"*Plancktime") ? 10 : 6
		prefix = unit[0,strlen(unit)-i]
	elseif(StringMatch(unit,"*Svedberg"))
		value = 1e-13
		i = StringMatch(unit,"*Svedberg") ? 9 : 3
		prefix = unit[0,strlen(unit)-i]
	elseif(StringMatch(unit,"*galacticyear"))	// galactic year = 230e6 years
		value = 230e6 * tropicalYear
		prefix = unit[0,strlen(unit)-13]
	elseif(StringMatch(unit,"*siderealday"))	// sidereal day = 23.9344699 hour
		value = 23.9344699 * hour
		prefix = unit[0,strlen(unit)-12]
	elseif(StringMatch(unit,"*siderealyear"))	// 1 sidereal year = 365.256363004 day
		value = 365.256363004 * day
		prefix = unit[0,strlen(unit)-13]
	elseif(StringMatch(unit,"*helek"))				// 1 helek = 3 + 1/3 seconds
		value = 10/3
		prefix = unit[0,strlen(unit)-6]
	elseif(StringMatch(unit,"*pahar"))				// 1 pahar = 3 hours
		value = 3 * hour
		prefix = unit[0,strlen(unit)-6]
	elseif(StringMatch(unit,"*minute") || StringMatch(unit,"*min") || StringMatch(unit,"*m"))	// minute = 60 seconds
		value = 60
		i = StringMatch(unit,"*minute") ? 7 : 2
		i = StringMatch(unit,"*min") ? 4 : i
		prefix = unit[0,strlen(unit)-i]
	elseif(StringMatch(unit,"*hour") || StringMatch(unit,"*hr") || StringMatch(unit,"*h"))	// hour
		value = hour
		i = StringMatch(unit,"*hour") ? 5 : 2
		i = StringMatch(unit,"*hr") ? 3 : i
		prefix = unit[0,strlen(unit)-i]
	elseif(StringMatch(unit,"*d") || StringMatch(unit,"*day"))	// day
		value = day
		i = StringMatch(unit,"*day") ? 4 : 2
		prefix = unit[0,strlen(unit)-i]
	elseif(StringMatch(unit,"*wk") || StringMatch(unit,"*week") )	// week
		value = 7*day
		i = StringMatch(unit,"*week") ? 5 : 3
		prefix = unit[0,strlen(unit)-i]
	elseif(StringMatch(unit,"*year") || StringMatch(unit,"*yr") || StringMatch(unit,"*y"))
		value = tropicalYear							// 1 year = 365.24219 days
		i = StringMatch(unit,"*year") ? 5 : 2
		i = StringMatch(unit,"*yr") ? 3 : i
		prefix = unit[0,strlen(unit)-i]
	else
		return defaultSeconds							// cannot find base value
	endif
	prefix = ReplaceString("micro",prefix,Gmu)	// Greek mu

	value *= SIprefix2factor(prefix)
	return (power==1) ? value : (value ^ power)
End
//	Function testAllTimes(SIprefix)
//		String SIprefix
//	
//		if (numtype(SIprefix2factor(SIprefix)))
//			printf "bad SI prefix '%s'\r",SIprefix
//			return 1
//		endif
//	
//		String list = "second;sec;s;minute;min;m;hour;hr;h;day;d;week;wk;fortnight;lunar month;lune;moon;year;yr;y;"
//		list += "olympiad;lustrum;indiction;decade;century;millennium;jiffy;shake;beat;Planck time;"
//		list += "Svedberg;galactic year;sidereal day;sidereal year;helek;pahar;abc"
//		String unit
//		Variable value, i, N=ItemsInList(list)
//		for (i=0;i<N;i+=1)
//			unit = StringFromList(i,list)
//			value = ConvertUnits2seconds(SIprefix+unit)
//			if (numtype(value))
//				printf "bad value for:  '%s' + '%s'\r",SIprefix,unit
//			endif
//		endfor
//		return 0
//	End


ThreadSafe Function ConvertUnits2kg(unit,[defaultMass])
	// returns conversion factor from unit to killo-grams
	String unit
	Variable defaultMass
	defaultMass = ParamIsDefault(defaultMass) ? NaN : defaultMass
	//
	//	kg, kilogram			1 kg
	//	g, gram					0.001 kg
	//	carat, ct				2e-4 kg						metric carat, 200 mg
	//	gr, grain				kgPerPound/7000.0
	//	fir, firkin				90.0*kgPerPound
	//	lb, lbm, #, pound		0.45359237 kg				Avoirdupois pound
	//	oz, ounce				kgPerPound/16.0			Avoirdupois ounce
	//	slug						kgPerPound*gEarth/(12.0*0.0254)
	//	st, stone				14.0*kgPerPound
	//	t, tonne, metricton 1000 kg						metric ton
	//	t.short, shortton	, ton	2000*kgPerPound		short ton (2000 pounds)
	//	t.long, long			2240.0*kgPerPound			long ton (2240 pounds)
	//	tlb, troyPound			5.760*0.06479891 kg		Troy pound
	//	toz, troyOunce			5.760*0.06479891/12.0 kg	Troy ounce
	//	amu, dalton				1.66053904e-27 kg			atomic mass unit
	//	mP, Planck				sqrt(hbar*c/GN)			Planck mass
	//	sun, sol, solar		1.9891E30 kg				solar mass
	//	me							9.10938356e-31 kg			electron mass
	//	mp							1.672621898e-27 kg		proton mass
	//	mn							1.674927471e-27 kg		neutron mass
	//	mmu						1.883531594e-28 kg		muon mass

	unit = ReplaceString(" ",unit,"")			// no spaces
	// check for powers  "^N"
	Variable ipow=strsearch(unit,"^",0), power=1
	if (ipow>=0)										// found a power
		power = str2num(unit[ipow+1,Inf])		// number after "^"
		unit = unit[0,ipow-1]						// string before "^"
	endif
	unit = RemoveEnding(unit,"s")				// none of the units end in 's'

	// alternate spellings:
	unit = ChangeStrEnding("metricton",unit,"t")

	String prefix
	Variable value=NaN, i = max(0,strlen(unit)-1)
	if(StringMatch(unit,"*Planck")) 			// Planck mass
		value = sqrt(hbar*c/GN)
		prefix = unit[0,strlen(unit)-7]
	elseif(StringMatch(unit,"*amu") || StringMatch(unit,"*dalton"))
		value = 1.66053904e-27
		i = StringMatch(unit,"*dalton") ? 7 : 4
		prefix = unit[0,strlen(unit)-i]
	elseif(StringMatch(unit,"*sun") || StringMatch(unit,"*sol") || StringMatch(unit,"*solar"))
		value = 1.9891E30								// solar mass
		i = StringMatch(unit,"*solar") ? 6 : 4
		prefix = unit[0,strlen(unit)-i]
	elseif(StringMatch(unit,"*mmu"))	 		// muon mass
		value = 1.883531594e-28
		prefix = unit[0,strlen(unit)-4]
	elseif(StringMatch(unit,"*stone") || StringMatch(unit,"*st"))
		value = 14 * kgPerPound					// stone = 14 pounds
		i = StringMatch(unit,"*stone") ? 6 : 3
		prefix = unit[0,strlen(unit)-i]
	elseif(StringMatch(unit,"*slug")) 			// slug, acceleration at earth surface is 9.80665 N
		value = kgPerPound*gStd/(12*inch)
		prefix = unit[0,strlen(unit)-5]
	elseif(StringMatch(unit,"*fir") || StringMatch(unit,"*firkin"))
		value = 90 * kgPerPound 					// firkin
		i = StringMatch(unit,"*firkin") ? 7 : 4
		prefix = unit[0,strlen(unit)-i]
	elseif(StringMatch(unit,"*me")) 			// electron mass
		value = 9.10938356e-31
		prefix = unit[0,strlen(unit)-3]
	elseif(StringMatch(unit,"*mp")) 			// proton mass
		value = 1.672621898e-27
		prefix = unit[0,strlen(unit)-3]
	elseif(StringMatch(unit,"*mn")) 			// neutron mass
		value = 1.674927471e-27
		prefix = unit[0,strlen(unit)-3]
	elseif(StringMatch(unit,"*grain") || StringMatch(unit,"*gr"))
		value = kgPerPound/7000					// grain
		i = StringMatch(unit,"*grain") ? 6 : 3
		prefix = unit[0,strlen(unit)-i]
	elseif(StringMatch(unit,"*carat") || StringMatch(unit,"*ct"))
		value = 2e-4									// metric carat
		i = StringMatch(unit,"*carat") ? 6 : 3
		prefix = unit[0,strlen(unit)-i]
	elseif(StringMatch(unit,"*troypound") || StringMatch(unit,"*tlb"))
		value = 5.760*0.06479891					// Troy pound
		i = StringMatch(unit,"*troypound") ? 10 : 4
		prefix = unit[0,strlen(unit)-i]
	elseif(StringMatch(unit,"*troyounce") || StringMatch(unit,"*toz"))
		value = 5.760*0.06479891/12.0			// Troy ounce
		i = StringMatch(unit,"*troyounce") ? 10 : 4
		prefix = unit[0,strlen(unit)-i]
	elseif(StringMatch(unit,"*ounce") || StringMatch(unit,"*oz"))
		value = kgPerPound/16.0					// Avoirdupois ounce
		i = StringMatch(unit,"*ounce") ? 6 : 3
		prefix = unit[0,strlen(unit)-i]
	elseif(StringMatch(unit,"*pound") || StringMatch(unit,"*lbm") || StringMatch(unit,"*lb") || StringMatch(unit,"*#"))
		value = kgPerPound							// Avoirdupois pound
		i = StringMatch(unit,"*#") ? 2 : NaN
		i = StringMatch(unit,"*lb") ? 3 : i
		i = StringMatch(unit,"*lbm") ? 4 : i
		i = StringMatch(unit,"*pound") ? 6 : i
		prefix = unit[0,strlen(unit)-i]
	elseif(StringMatch(unit,"*longton") || StringMatch(unit,"*long") || StringMatch(unit,"*t.long"))
		value = 2240 * kgPerPound					// long ton, (2240 pounds)
		i = StringMatch(unit,"*long") ? 5 : NaN
		i = StringMatch(unit,"*t.long") ? 7 : i
		i = StringMatch(unit,"*longton") ? 8 : i
		prefix = unit[0,strlen(unit)-i]
	elseif(StringMatch(unit,"*shortton") || StringMatch(unit,"*t.short") || StringMatch(unit,"*short") || StringMatch(unit,"*ton"))
		value = 2000 * kgPerPound					// short ton, (2000 pounds)
		i = StringMatch(unit,"*ton") ? 4 : NaN
		i = StringMatch(unit,"*short") ? 6 : i
		i = StringMatch(unit,"*t.short") ? 8 : i
		i = StringMatch(unit,"*shortton") ? 9 : i
		prefix = unit[0,strlen(unit)-i]
	elseif(StringMatch(unit,"*tonne") || StringMatch(unit,"*t"))
		value = 1000									// metric ton, (1000 kg)
		i = StringMatch(unit,"*tonne") ? 6 : 2
		prefix = unit[0,strlen(unit)-i]
	elseif(StringMatch(unit,"*gram") || StringMatch(unit,"*g"))
		value = 0.001									// gram
		i = StringMatch(unit,"*gram") ? 5 : 2
		prefix = unit[0,strlen(unit)-i]
	else
		return defaultMass							// cannot find base value
	endif

	value *= SIprefix2factor(prefix)
	return (power==1) ? value : (value ^ power)
End


ThreadSafe Function ConvertUnits2Joules(unit,[defaultEnergy])
	// returns conversion factor from unit to Joules
	String unit
	Variable defaultEnergy
	defaultEnergy = ParamIsDefault(defaultEnergy) ? NaN : defaultEnergy	// no default value

	// Joule								1 J							(kg * m^2 * s^-2 = N*m)
	// erg								1e-7 J
	// cal, calorie					4.184 J						(Thermochemical calorie)
	// Btu, BTU							1055.06 J
	// Wh, watt-hour					3600 J						usually used as 'kiloWatt-hour'
	// eV, electron-volt			1.6021766208e-19			e_C = 1.6021766208e-19 C
	// Ry, Rydberg						2.179872325e-18 J			R(inf) * hc = (me * e^4) / (8 * epsilon0^2 * h^2) = ( 1/2 * me * c^2 * alpha^2 )
	// Ha, Hartree						2*Ry								
	// ft-lbf, ft-lb, Foot-pound	1.3558179483314 J			(12*inch) *(g*kg/pound)) = ((12*0.0254) * (9.80665*0.45359237))
	// therm								1e5 BTU		
	// quad								1e15 BTU
	// TWyr, terawatt-year			31.556925216e18 J			= 365.24219 * 24*3600 * 1e12
	// 'TNT' or 'kg of TNT'		4.184e6 J					(1 kg of TNT), releases 4.184e6 J of energy = cal*1e6 = Mcal
	// 'ton of TNT', , 'kiloton, ton	4.184e9 J			(1 metric ton of TNT),  '1 kiloton = (1000 metric tons of TNT)
	// Planck, Planck energy		1956113859.56355 J		sqrt(hbar * c^5 / GN)
	// foe, Bethe						1e44 J						= 10^51 erg

	Variable BTU = 1055.06							// one BTU (J) [[definition of BTU]
	Variable Ry_hc = 2.179872325e-18			// R(inf) * hc, Rydberg energy (J)
	Variable cal = 4.184							// 1 Thermochemical calorie = 4.184 J

	// check for powers  "^N"
	Variable ipow=strsearch(unit,"^",0), power=1
	if (ipow>=0)										// found a power
		power = str2num(unit[ipow+1,Inf])		// number after "^"
		unit = unit[0,ipow-1]						// string before "^"
	endif

	unit = RemoveEnding(unit,"s")				// none of the units end in 's'
	unit = ReplaceString(" of ",unit," ")
	unit = ReplaceString("-",unit,"")			// no dashes
	unit = ReplaceString(" ",unit,"")			// no spaces

	unit = ReplaceString("tonTNT",unit,"kTNT")
	unit = ReplaceString("ton",unit,"kTNT")
	unit = ReplaceString("kg",unit,"TNT")

	unit = ReplaceString("kilo",unit,"k")
	unit = ReplaceString("mega",unit,"M")
	unit = ReplaceString("giga",unit,"G")

	String prefix=""
	Variable value=NaN, i = max(0,strlen(unit)-1)

	if(StringMatch(unit,"*Joule") || StringMatch(unit,"*J"))
		value = 1
		i = strsearch(unit,"J",inf,3)
		prefix = unit[0,i-1]
	elseif(StringMatch(unit,"*Planck*"))
		value = sqrt(hbar*c^5/GN)					// Planck Energy = 1956113859.56355
		i = strsearch(unit,"Planck",0,2)
		prefix = unit[0,i-1]
	elseif(StringMatch(unit,"*cal*"))
		value = cal										// (Thermochemical calorie)
		i = strsearch(unit,"cal",0,2)
		prefix = unit[0,i-1]
	elseif(StringMatch(unit,"*BTU"))
		value = BTU
		i = strsearch(unit,"BTU",0,2)
		prefix = unit[0,i-1]
	elseif(StringMatch(unit,"*watthour"))
		value = 3600
		i = strsearch(unit,"watt",0,2)
		prefix = unit[0,i-1]
	elseif(StringMatch(unit,"*Wh"))				// Watt-hour
		value = 3600
		i = strsearch(unit,"Wh",0,2)
		prefix = unit[0,i-1]
	elseif(StringMatch(unit,"*eV") || StringMatch(unit,"*electronVolt"))
		value = e_C										// elementary charge
		i = strsearch(unit,"eV",0,2)
		i = i<0 ? strsearch(unit,"electron",0,2) : i
		prefix = unit[0,i-1]
	elseif(StringMatch(unit,"*Ry") || StringMatch(unit,"*Rydberg"))
		value = Ry_hc									// R(inf) * hc
		i = strsearch(unit,"Ry",inf,3)
		prefix = unit[0,i-1]
	elseif(StringMatch(unit,"*Ha") || StringMatch(unit,"*Hartree"))
		value = 2*Ry_hc								// 2*Ry
		i = strsearch(unit,"Ha",0,2)
		prefix = unit[0,i-1]
	elseif(StringMatch(unit,"*erg"))
		value = 1e-7
		i = strsearch(unit,"erg",0,2)
		prefix = unit[0,i-1]
	elseif(StringMatch(unit,"*therm"))
		value = 1e5 * BTU
		i = strsearch(unit,"therm",0,2)
		prefix = unit[0,i-1]
	elseif(StringMatch(unit,"*quad"))
		value = 1e15 * BTU
		i = strsearch(unit,"quad",0,2)
		prefix = unit[0,i-1]
	elseif(StringMatch(unit,"*foe") || StringMatch(unit,"*Bethe"))
		value = 1e44									// 10^51 ergs
		i = strsearch(unit,"foe",0,2)
		i = i<0 ? strsearch(unit,"Bethe",0,2) : i
		prefix = unit[0,i-1]
	elseif(StringMatch(unit,"*TWyr") || StringMatch(unit,"*terawattYear"))
		value = tropicalYear * 1e12		// (days in year) * (seconds in day) * 1e12 = 31.556925216e18
		i = strsearch(unit,"TWyr",0,2)
		i = i<0 ? strsearch(unit,"terawatt",0,2) : i
		prefix = unit[0,i-1]
	elseif(StringMatch(unit,"*ftlb*") || StringMatch(unit,"*FootPound*"))
		// distance * pound force
		value = (12*inch) * (gStd*kgPerPound)	// (12*inch) * (9.80665 * 0.45359237)
		i = strsearch(unit,"ft",0,2)
		i = i<0 ? strsearch(unit,"Foot",0,2) : i
		prefix = unit[0,i-1]
	elseif(StringMatch(unit,"*TNT"))
		value = cal*1e6								// (1 kg of TNT) = 1 Mcal
		i = strsearch(unit,"TNT",0,2)
		prefix = unit[0,i-1]
	else
		return defaultEnergy						// cannot find base value, use default
	endif

	value *= SIprefix2factor(prefix)
	return (power==1) ? value : (value ^ power)
End
//
//	Function test_ConvertUnits2Joules()
//		String AllValues  = "J:1;Joule:1;erg:1e-7;cal:4.184;calorie:4.184;Btu:1055.06;BTU:1055.06;kWh:3.6e6;kilowatt-hour:3.6e6;"
//		AllValues += "eV:1.6021766208e-19;electron-volt:1.6021766208e-19;Ry:2.179872325e-18;Rydberg:2.179872325e-18;"
//		AllValues += "Ha:4.35974465e-18;Hartree:4.35974465e-18;ft-lbf:1.3558179483314;ft-lb:1.3558179483314;Foot-pound:1.3558179483314;"
//		AllValues += "therm:1055.06e5;quad:1055.06e15;TWyr:31.556925216e18;terawatt-year:31.556925216e18;ton of TNT:4.184e9;ton:4.184e9;"
//		AllValues += "Planck:1956113859.56355;Planck energy:1956113859.56355;foe:1e44;Bethe:1e44"
//		String SI=";d;c;m;"+Gmu+";n;p;f;a;z;y;h;k;M;G;T;P;E;Z;Y"		// first item is empty, SIprefix2factor("")==1
//	
//		String prefix, unitBase
//		Variable ip, Np=ItemsInList(SI), valExpect, valCalc
//		Variable Nerr=0, i, Nv=ItemsInlist(AllValues)
//		for (i=0;i<Nv;i+=1)
//			unitBase = StringFromList(0,StringFromList(i,AllValues),":")
//	
//			for (ip=0;ip<Np;ip+=1)					// check all SI prefixes
//				prefix = StringFromList(ip,SI)
//				valExpect = SIprefix2factor(prefix) * NumberByKey(unitBase,AllValues)
//				valCalc = ConvertUnits2Joules(prefix+unitBase)
//				if (abs(valCalc-valExpect)/valExpect > 1e-15 || numtype(valCalc+valExpect))
//					printf "%s\t\t\tCalc=%g,\t\texpected=%g\r",prefix+unitBase,valCalc,valExpect
//					Nerr += 1
//				endif
//			endfor
//		endfor
//		if (Nerr)
//			printf "\r  completed with %g errors\r",Nerr
//		else
//			print "completed with NO errors"
//		endif
//	End


ThreadSafe Function ConvertUnits2Pascal(unit,[defaultP])
	// returns conversion factor from unit to Pascal
	String unit
	Variable defaultP
	defaultP = ParamIsDefault(defaultP) ? NaN : defaultP
	//
	//	Pascal, Pa				1 Pascal					(= 1  N * m^-2)
	//	bar						1e5
	//	atmosphere, atm		stdAtmosphere			(=101325 = standard atmosphere)
	// psi						gStd*kgPerPound/(inch^2)	(pound/square inch)
	//	Torr						stdAtmosphere/760		(exact defintion is 1atm/760, ~mm of Hg)
	//	mHg						stdAtmosphere/0.760	(exact defintion is 1000*atm/760, ~m of Hg)
	//	inH2O, wc				inH2O						(=249.082, at 4�C, also inch or inches of water)
	// ftH2O,'ft of water'	12*inH2O					(also foot or feet of water)
	// mH2O,mwater				inH2O/inch				(1 m of water)
	//	msw						1e4						(meters of sea water)
	//	fsw						1e4*0.0254*12			(feet of sea water)

	unit = ReplaceString(" ",unit,"")			// no spaces

	// check for powers  "^N"
	Variable ipow=strsearch(unit,"^",0), power=1
	if (ipow>=0)										// found a power
		power = str2num(unit[ipow+1,Inf])		// number after "^"
		unit = unit[0,ipow-1]						// string before "^"
	endif

	// fix spellings
	unit = ChangeStrEnding("Pascal",unit,"Pa")
	unit = ReplaceString("of",unit,"")
	unit = ReplaceString("pounds",unit,"pound")
	unit = ReplaceString("pound/squareinch",unit,"psi")
	unit = ReplaceString("poundpersquareinch",unit,"psi")
	unit = ReplaceString("poundpersquareinch",unit,"psi")
	unit = ReplaceString("poundspersquareinch",unit,"psi")
	unit = ReplaceString("seawater",unit,"sw")
	unit = ReplaceString("water",unit,"H2O")
	unit = ReplaceString("mercury",unit,"Hg")
	unit = ReplaceString("wc",unit,"inH2O")
	unit = ReplaceString("inches",unit,"in")
	unit = ReplaceString("inch",unit,"in")
	unit = ReplaceString("foot",unit,"ft")
	unit = ReplaceString("feet",unit,"ft")
	unit = ReplaceString("metres",unit,"m")	// British spelling
	unit = ReplaceString("meters",unit,"m")	// meters
	unit = ReplaceString("metre",unit,"m")
	unit = ReplaceString("meter",unit,"m")
	unit = ChangeStrEnding("ftsw",unit,"fsw")
	unit = ChangeStrEnding("atmosphere",unit,"atm")
	unit = RemoveEnding(unit,"s")				// remove any trailing "s"

	String prefix=""
	Variable value=NaN, i=max(0,strlen(unit)-1)

	if(StringMatch(unit,"*Pa"))					// ends in 'Pa', means Pascal
		value = 1
		prefix = unit[0,i-2]
	elseif(StringMatch(unit,"*bar"))			// 1 bar = 1e5 Pa
		value = 1e5
		prefix = unit[0,i-3]
	elseif(StringMatch(unit,"*atm"))			// 1 atm = 101325 Pa
		value = stdAtmosphere
		prefix = unit[0,i-3]
	elseif(StringMatch(unit,"*psi"))			// pound per square inch
		value = gStd*kgPerPound/(inch^2)
		prefix = unit[0,i-3]
	elseif(StringMatch(unit,"*Torr"))			// 1 Torr = 1atm/760
		value = stdAtmosphere/760
		prefix = unit[0,i-4]
	elseif(StringMatch(unit,"*mHg"))			// 1 mHg = 1000atm/760
		value = 1000*stdAtmosphere/760
		prefix = unit[0,i-3]
	elseif(StringMatch(unit,"*inH2O"))			// inches of water (at 4�C)
		value = inH2O
		prefix = unit[0,i-5]
	elseif(StringMatch(unit,"*ftH2O"))			// feet of water (at 4�C)
		value = 12*inH2O
		prefix = unit[0,i-5]
	elseif(StringMatch(unit,"*mH2O"))			// meters of water (at 4�C)
		value = inH2O/inch
		prefix = unit[0,i-4]
	elseif(StringMatch(unit,"*msw"))			// meters of sea water
		value = 1e4
		prefix = unit[0,i-3]
	elseif(StringMatch(unit,"*fsw"))			// feet of sea water
		value = 1e4*0.0254*12
		prefix = unit[0,i-3]
	else
		return defaultP								// cannot find base value, use default
	endif

	value *= SIprefix2factor(prefix)
	return (power==1) ? value : (value ^ power)
End
//
//	Function test_ConvertUnits2Pascal()
//		String AllValues  = "Pa:1;Pascal:1;bar:1e5;bar:1e5;atm:101325;atmosphere:101325;Torr:133.322368421053;"
//		AllValues += "psi:6894.75729316836;pounds per square inch:6894.75729316836;pound per square inch:6894.75729316836;"
//		AllValues += "mHg:133322.368421053;meters of mercury:133322.368421053;meter of mercury:133322.368421053;meter of Hg:133322.368421053;"
//		AllValues += "inH2O:249.082;wc:249.082;inches of H2O:249.082;inches of water:249.082;inches water:249.082;in water:249.082;"
//		AllValues += "ftH2O:2988.984;feet of H2O:2988.984;feet of water:2988.984;foot water:2988.984;foot of water:2988.984;"
//		AllValues += "msw:1e4;meters of sea water:1e4;fsw:3048;feet of sea water:3048;foot sea water:3048;"
//		AllValues += "mH2O:9806.37795275591;meters of H2O:9806.37795275591;meters of water:9806.37795275591;meter of water:9806.37795275591;"
//		String SI=";d;c;m;"+Gmu+";n;p;f;a;z;y;h;k;M;G;T;P;E;Z;Y"		// first item is empty, SIprefix2factor("")==1
//	
//		String prefix, unitBase
//		Variable ip, Np=ItemsInList(SI), valExpect, valCalc
//		Variable Nerr=0, i, Nv=ItemsInlist(AllValues)
//		for (i=0;i<Nv;i+=1)
//			unitBase = StringFromList(0,StringFromList(i,AllValues),":")
//			for (ip=0;ip<Np;ip+=1)					// check all SI prefixes
//				prefix = StringFromList(ip,SI)
//				valExpect = SIprefix2factor(prefix) * NumberByKey(unitBase,AllValues)
//				valCalc = ConvertUnits2Pascal(prefix+unitBase)
//				if (abs(valCalc-valExpect)/valExpect > 1e-14 || numtype(valCalc+valExpect))
//					printf "%s\t\t\tCalc=%g,\t\texpected=%g\r",prefix+unitBase,valCalc,valExpect
//					Nerr += 1
//				endif
//			endfor
//		endfor
//		if (Nerr)
//			printf "\r  completed with %g errors\r",Nerr
//		else
//			print "completed with NO errors"
//		endif
//	End


ThreadSafe Function ConvertAngleUnits(angle,unitIn,unitOut, [defaultUnit])	// returns angle(unitIn) in units of (unitOut)
	// ConvertAngleUnits(0.001, "degree", "�rad")	 --> 17.4533,		0.001*PI/180*1e6 = 17.453
	String unitIn, unitOut			// input and output units, both default to defaultUnit if empty
	Variable angle						// input angle  in units of (unitIn)
	String defaultUnit				// default units for unitIn & unitOut
	defaultUnit = SelectString(ParamIsDefault(defaultUnit), defaultUnit, "")
	//	degree, deg, �			1 degree (360 degrees in a circle)
	// radian, rad, r			180/PI
	// circle, cir				360
	// arcmin, min '			1/60  (1 degree = 60 arcminutes)
	// arcsec, sec	 " ''		1/3600  (1 arcmin = 60 arcseconds)
	// grad						360/400 (400 grad = 1 circle)
	// cos, cosine, co		cos(degree) = angle   SPECIAL	
	unitIn = SelectString(strlen(unitIn), defaultUnit, unitIn)
	unitOut = SelectString(strlen(unitOut), defaultUnit, unitOut)
	Variable factor, isCos, out, degree
	degree = angle * ConvertUnits2Degree(unitIn, isCos=isCos)	// conversion factor from unitIn to degrees
	degree = isCos ? acos(degree)*180/PI : degree

	factor = ConvertUnits2Degree(unitOut, isCos=isCos)	// conversion factor from degrees to unitOut
	out = isCos ? cos(degree*PI/180) : degree				// convert degree-->unitOut, also handles cosine
	out /= factor
	return out
End
//
ThreadSafe Static Function ConvertUnits2Degree(unit,[defaultFactor, isCos])
	// returns conversion factor from unit to degree for angles
	String unit
	Variable defaultFactor
	Variable &isCos				// if present, this is set to 1 if cosing, 0 if not

	defaultFactor = ParamIsDefault(defaultFactor) ? NaN : defaultFactor
	//
	//	degree, deg, �			1 degree (360 degrees in a circle)
	// radian, rad, r			180/PI
	// circle, cir				360
	// arcmin, min				1/60  (1 degree = 60 arcminutes)
	// arcsec, sec				1/3600  (1 arcmin = 60 arcseconds)
	// grad						360/400 (400 grad = 1 circle)
	// cos, cosine, co		cos(degree) = angle   SPECIAL

	unit = ReplaceString(" ",unit,"")				// no spaces, so 'arc min' --> 'arcmin'
	// check for powers  "^N"
	Variable ipow=strsearch(unit,"^",0), power=1
	if (ipow>=0)											// found a power
		power = str2num(unit[ipow+1,Inf])			// number after "^"
		unit = unit[0,ipow-1]							// string before "^"
	endif

	// fix spellings
	unit = RemoveEnding(unit,".")					// e.g. 'deg.' --> 'deg'
	unit = ReplaceString("cosine",unit,"Q")		// internally use Q for cosine
	unit = ReplaceString("cos",unit,"Q")
	unit = ChangeStrEnding("co",unit,"Q")	
	unit = ReplaceString("ofarc",unit,"")
	unit = ReplaceString("arc",unit,"")

	unit = ReplaceString("degree",unit,"deg")
	unit = ReplaceString(DEGREESIGN,unit,"deg")
	unit = ReplaceString("\241",unit,"deg")		// no degree signs, Mac opt-shift-9
	unit = ReplaceString("\xC2\xB0",unit,"deg")// no degree signs, UTF8

	unit = ReplaceString("GRAD",unit,"O",1)		// Grad is assumed to be (Giga-radian)
	unit = ReplaceString("grad",unit,"O",1)		// internally use 'O' for grad
	if (strsearch(unit,"circle",0,2)<0)			// circle NOT present
		unit = ReplaceString("cir",unit,"circle")	// change cir --> circle
	endif
	unit = ReplaceString("radian",unit,"R")		// internally use 'R' for radian
	unit = ReplaceString("rad",unit,"R")
	unit = ChangeStrEnding("r",unit,"R")
	unit = ReplaceString("&quot",unit,"sec")
	unit = ReplaceString("\"",unit,"sec")
	unit = ReplaceString("''",unit,"sec")
	unit = ReplaceString("'",unit,"min")
	unit = ReplaceString("minute",unit,"min")
	unit = ReplaceString("second",unit,"sec")
	unit = RemoveEnding(unit,"s")					// remove any trailing "s", e.g. 'degrees' --> 'degree'

	String prefix=""
	Variable value=NaN, i=max(0,strlen(unit)-1), isCosLocal=0

	if (StringMatch(unit,"*deg"))					// ends in 'deg', degrees
		value = 1
		prefix = unit[0,i-3]
	elseif(StringMatch(unit,"*circle"))			// 1 circle = 360 deg
		value = 360
		prefix = unit[0,i-6]
	elseif(StringMatch(unit,"*min"))				// 1 arcmin = 1/60 deg
		value = 1/60
		prefix = unit[0,i-3]
	elseif(StringMatch(unit,"*sec"))				// 1 arcsec = 1/3600 deg
		value = 1/3600
		prefix = unit[0,i-3]
	elseif(StringMatch(unit,"*O"))					// 1 grad = 360/400 deg (400 grad = 1 circle)
		value = 360/400
		prefix = unit[0,i-1]
	elseif(StringMatch(unit,"*Q"))					// angle = cos(degree)   SPECIAL
		isCosLocal = 1
		value = 1
		prefix = unit[0,i-1]
	elseif(StringMatch(unit,"*R"))					// 1 rad = 180/PI deg
		value = 180/PI
		prefix = unit[0,i-1]
	else
		return defaultFactor							// cannot find base value, use default
	endif
	if (!ParamIsDefault(isCos))
		isCos = isCosLocal
	endif

	value *= SIprefix2factor(prefix)
	if (power !=1 )
		value = isCosLocal ? NaN : (value ^ power)
	endif
	return value
End
//
//	Function test_AngleUnits()
//		Variable err = 0
//	
//		String ch, all="�;deg;degree;degree;rad;radian;grad;min;arcmin;';sec;arcsec;\";'';mcos;Cosine"
//		Variable i
//		for (i=0;i<ItemsInList(all);i+=1)
//			ch = StringFromList(i,all)
//			err += test_ConvertAngle(0.5,ch,"degree")
//		endfor
//		print " "
//		all="�;deg;degree;degree;rad;grad;min;sec;';\";'';mcos;Cosine"
//		for (i=0;i<ItemsInList(all);i+=1)
//			ch = StringFromList(i,all)
//			err += test_ConvertAngle(2,"degree",ch)
//		endfor
//		print " "
//		err += test_ConvertAngle(0.5,"mcos","degree", expected=89.9, explanation="intentional")
//		err += test_ConvertAngle(0.0005,"cos", "degree", expected=89.9713521090498)
//		err += test_ConvertAngle(0.5,"mcos","degree", expected=89.9713521090498)
//		err += test_ConvertAngle(0.5,"mco","degree", expected=89.9713521090498)
//		err += test_ConvertAngle(0.5,"mCosine","degree", expected=89.9713521090498)
//	
//		print " "
//		err += test_ConvertAngle(0.01,"deg","cos")
//		err += test_ConvertAngle(89.,"deg","cos")
//		err += test_ConvertAngle(89.,"deg","mcos")
//		err += test_ConvertAngle(180.,"","rad", explanation="no unitIN, empty")
//		err += test_ConvertAngle(180.,"rad","", explanation="no unitOUT, empty")
//		print " "
//		err += test_ConvertAngle(0.05,"cos","deg")
//		err += test_ConvertAngle(50,"mcos","deg")
//		err += test_ConvertAngle(0.01,"mCosine","mdegree")
//		err += test_ConvertAngle(89999.42704220486,"mdegree","mCosine")
//		err += test_ConvertAngle(0,"co","mcir")
//		print " "
//		err += test_ConvertAngle(0.324,"M\"","deg")
//		err += test_ConvertAngle(324,"k''","deg")
//		err += test_ConvertAngle(5.4,"k'","deg")
//		err += test_ConvertAngle(1,"arcsec","�rad")
//		err += test_ConvertAngle(4.85,"�r","arcsec")
//		err += test_ConvertAngle(1,"min","mrad")
//		err += test_ConvertAngle(1,"min","mr")
//		err += test_ConvertAngle(1,"min","r")
//		print " "
//		err += test_ConvertAngle(1, "deg^2", "min^2", expected=3600)
//		err += test_ConvertAngle(1, "min^2", "sec^2", expected=3600)
//		err += test_ConvertAngle(1, "min^-1", "deg^-1", expected=60)
//		err += test_ConvertAngle(1, "sec^-2", "min^-2", expected=3600)
//		err += test_ConvertAngle(.05,"s","deg", explanation="'s' is not an angular unit")
//		if (err)
//			print "\r  	*********** Finished with ERROR ! ***********"
//		else
//			print "	----------- Finished with no errors -----------"
//		endif
//		return err
//	End
//
//	Static Function test_ConvertAngle(valueIN, unitIN, unitOUT, [expected, explanation])
//		Variable valueIN
//		String unitIN, unitOUT
//		Variable expected
//		String explanation
//		expected = ParamIsDefault(expected) ? NaN : expected
//		explanation = SelectString(ParamIsDefault(explanation), explanation, "")
//	
//		Variable tol = 1e-15, err=0
//		Variable out = ConvertAngleUnits(valueIN,unitIN,unitOUT)
//	
//		if (!ParamIsDefault(expected))
//			if (abs((out-expected)/out)>tol)
//				err = 1
//			endif
//		else
//			err = numtype(out)>0 
//		endif
//	
//		if (!err)
//			printf "     %g (%s)  -->  %g (%s)\r", valueIN,unitIN,out,unitOUT
//			return 0
//		else
//			String errStr = SelectString(strlen(explanation), "ERR  ", "     ")
//			printf "%sINVALID -- %g (%s)  -->  %s     should be %g not %g", errStr,valueIN,unitIN, unitOUT, expected, out
//			if (strlen(explanation))
//				printf "\t\texplanation = \"%s\"", explanation
//			endif
//			printf "\r"
//		endif
//		return (strlen(explanation)<1)
//	End
//
//	ThreadSafe Function endswith(str, ending, case)		// returns True if str ends with ending
//		String str				// string to check
//		String ending			// ending to find
//		Variable case			// true=match case,  false, ignore case
//	
//		Variable i, Nstr=strlen(str), Nend=strlen(ending)
//		if (Nend<1)
//			return 1				// everything ends with nothing
//		endif
//		i = strsearch(str,ending,Inf,case?1:3)
//		return i == (Nstr-Nend) && (i>=0)
//	End

//  ============================== End of Unit Conversions =============================  //
//  ====================================================================================  //


ThreadSafe Function/T RomanNumeral(j)	// convert integer j to a Roman Numeral String, NO upper limit, so watch out for string length
	Variable j

	String str = SelectString(j<0,"","-")		// start str with the sign
	j = round(abs(j))

	if (j<1 || numtype(j))					// retrns "" for zero, �Inf, NaN
		str = ""
	elseif (j>=1000)
		str += "M"+RomanNumeral(j-1000)	// add an M & remove 1000
	elseif (j>=900)
		str += "CM"+RomanNumeral(j-900)	// add a CM & remove 900
	elseif (j>=500)
		str += "D"+RomanNumeral(j-500)	// add a D & remove 500
	elseif (j>=400)
		str += "CD"+RomanNumeral(j-400)	// add a CD & remove 400
	elseif (j>=100)
		str += "C"+RomanNumeral(j-100)	// add a C & remove 100
	elseif (j>=90)
		str += "XC"+RomanNumeral(j-90)	// add a XC & remove 90
	elseif (j>=50)
		str += "L"+RomanNumeral(j-50)	// add a L & remove 50
	elseif (j>=40)
		str += "XL"+RomanNumeral(j-40)	// add a XL & remove 40
	elseif (j>=10)
		str += "X"+RomanNumeral(j-10)	// add a X & remove 10
	elseif (j>=9)
		str += "IX"+RomanNumeral(j-9)	// add a IX & remove 9
	elseif (j>=5)
		str += "V"+RomanNumeral(j-5)		// add a V & remove 5
	elseif (j>=4)
		str += "IV"+RomanNumeral(j-4)	// add a IV & remove 4
	elseif (j>=1)
		str += "I"+RomanNumeral(j-1)		// add a I & remove 1
	endif
	return str
End


ThreadSafe Function RomanNumeral2Int(str)	// convert a Roman Numeral string to an integer
	String str
	if (strlen(str)<1)						// empty strings return NaN, there is no zero
		return NaN
	endif
	str = UpperStr(str )

	Variable iM=char2num("M"),iC=char2num("C"),iX=char2num("X"),iI=char2num("I")
	Variable j, isign=1
	if (char2num(str)==char2num("-"))
		str = str[1,Inf]						// strip off leading "-", and set isign to -1
		isign = -1
	endif

	for (;char2num(str[0])==iM;)			// fournd all leading "M", add 1000
		str = str[1,Inf]
		j += 1000
	endfor

	if (strsearch(str,"CM",0)==0)		// fournd a leading "CM", add 900
		str = str[2,Inf]
		j += 900
	endif
	if (strsearch(str,"D",0)==0)			// fournd a leading "D", add 500
		str = str[1,Inf]
		j += 500
	endif
	if (strsearch(str,"CD",0)==0)		// fournd a leading "CD", add 400
		str = str[2,Inf]
		j += 400
	endif
	for (;char2num(str[0])==iC;)			// fournd all leading "C", add 100
		str = str[1,Inf]
		j += 100
	endfor

	if (strsearch(str,"XC",0)==0)		// fournd a leading "XC", add 90
		str = str[2,Inf]
		j += 90
	endif
	if (strsearch(str,"L",0)==0)			// fournd a leading "L", add 50
		str = str[1,Inf]
		j += 50
	endif
	if (strsearch(str,"XL",0)==0)		// fournd a leading "XL", add 40
		str = str[2,Inf]
		j += 40
	endif
	for (;char2num(str[0])==iX;)			// fournd all leading "X", add 10
		str = str[1,Inf]
		j += 10
	endfor

	if (strsearch(str,"IX",0)==0)		// fournd a leading "IX", add 9
		str = str[2,Inf]
		j += 9
	endif
	if (strsearch(str,"V",0)==0)			// fournd a leading "V", add 5
		str = str[1,Inf]
		j += 5
	endif
	if (strsearch(str,"IV",0)==0)		// fournd a leading "IV", add 4
		str = str[2,Inf]
		j += 4
	endif
	for (;char2num(str[0])==iI;)			// fournd all leading "I", add 1
		str = str[1,Inf]
		j += 1
	endfor

	j = strlen(str) ? NaN : j*isign		// check that str is all used up & fix sign
	return j
End

//  ======================= End of some general utility functions ========================  //
//  ======================================================================================  //


Function init_UtiltyJZT()
	if (DataFolderExists("root:Packages:UtilityJZT"))
		return 0
	endif
	NewDataFolder/O root:Packages
	NewDataFolder/O root:Packages:UtilityJZT

	if (exists("root:Packages:UtilityJZT:markerIndex")!=2)
		Variable/G root:Packages:UtilityJZT:markerIndex = 0
	endif
	return 0
End



//  ======================================================================================  //
//  ======================== Start of legacy deprecated functions ========================  //

// use new ISOtime2IgorEpoch(iso, NOTZ=1) instead
//ThreadSafe Function ISOtime2IgorEpoch(iso)	// convert ISO8601 string to an Igor Epoch (error returns NaN)
//	String iso									// format     2013-10-02T01:22:35  (the seconds are optional)
//
//	Variable year=NaN,month=NaN,day=NaN, hr=NaN,mn=NaN,se=NaN
//	sscanf iso,"%4d-%2d-%2dT%2d:%2d:%2d", year,month,day,hr,mn,se
//	Variable N=V_flag
//
//	if (N<3 || numtype(year+month+day))
//		return NaN
//	endif
//	Variable epoch=date2secs(year, month, day )
//	if (N>=5)
//		epoch += hr*3600
//		epoch += mn*60
//		se = (N>=6) ? se : 0
//		epoch += se
//	endif
//	return epoch
//End


Function/S CornerStamp1_()		// ONLY for backwards compatibility, DEPRECATED, use JZTutil#CornerStampWindow()
	PathInfo home
	return S_path+IgorInfo(1)
End

Static Function isWaveOnGraph(w,win)				// DEPRECATED use CheckDisplayed directly
	Wave w
	String win
	CheckDisplayed/W=$win w
	return (V_flag>0)
End
//Function isWaveOnGraph(w,win)
//	Wave w
//	String win
//	if (!WaveExists(w))
//		return 0
//	endif
//	String name, name0=GetWavesDataFolder(w,2)
//	String ilist
//	Variable i,Ni
//		ilist = ImageNameList(win,";")		// list of images in this window
//		Ni = ItemsInList(ilist)
//		if (strsearch(win,"#",0)<0)
//			for (i=0;;i+=1)					// first check for graphs (x or y trace)
//				Wave wi = WaveRefIndexed(win,i,3)
//				if (!WaveExists(wi))
//					break
//				endif
//				if (stringmatch(GetWavesDataFolder(wi,2),name0))
//					return 1
//				endif
//			endfor
//		endif
//		for (i=0;i<ItemsInList(ilist);i+=1)// next check all the images
//			name = GetWavesDataFolder(ImageNameToWaveRef(win,StringFromList(i,ilist)),2)
//			if (stringmatch(name,name0))
//				return 1
//			endif
//		endfor
//	return 0
//End


//	DEPRECATED,  use the above three TrimBoth, TrimEnd, & TrimFront instead
//
//		These three functions have been moved here from from LoadJZTtaggedData.ipf, BeamProcedures, LoadEPICSscans.ipf, and oldJZTdataLoad.ipf.
//			spec.ipf uses its own version of these so it is stand-alone
//
ThreadSafe Function/T TrimFrontBackWhiteSpace(str)
	String str
	str = TrimLeadingWhiteSpace(str)
	str = TrimTrailingWhiteSpace(str)
	return str
End
//
ThreadSafe Function/T TrimLeadingWhiteSpace(str)
	String str
	Variable i, N=strlen(str)
	for (i=0;char2num(str[i])<=32 && i<N;i+=1)	// find first non-white space
	endfor
	return str[i,Inf]
End
//
ThreadSafe Function/T TrimTrailingWhiteSpace(str)
	String str
	Variable i
	for (i=strlen(str)-1; char2num(str[i])<=32 && i>=0; i-=1)	// find last non-white space
	endfor
	return str[0,i]
End


ThreadSafe Function compareWaves(a,b)		// DEPRECATED,  just use EqualWaves(a,b,1)
	WAVE a,b
	return EqualWaves(a,b,1)
End
//ThreadSafe Function compareWaves(a,b)
//	WAVE a,b
//	Variable n=numpnts(a)
//	Variable i
//	if (n!=numpnts(b))
//		return 0
//	endif
//	for (i=0;i<n;i+=1)
//		if (numtype(a[i]))
//			if (numtype(a[i])!=numtype(b[i]))
//				return 0
//			endif
//		elseif (a[i]!=b[i])
//			return 0
//		endif
//	endfor
//	return 1
//End


// do not use num2sexigesmal() anymore, for new code use Secs2Time() directly.
ThreadSafe Function/S num2sexigesmal(seconds,places)	// convert seconds into a hh:mm:ss.sss   string
	Variable seconds
	Variable places
	return Secs2Time(seconds,5,places)
End


// DO NOT use this anymore.  In your code, do the test directly with:
//
//	if (strlen(GetRTStackInfo(2))<1)	// do this if this routine is at top level
//		// code goes here
//	endif
//
//
// returns 1 if the  calling function was invoked from a  menu item or command line (otherwise 0)
Function topOfStack()
	return strlen(GetRTStackInfo(2))<1
End

//Function abc()
//	String cmd, pathStr="\\\"/Users/tischler/data/Sector 34/July 4, 2006/EW5/recon/\\\""
//	print PathStr
//	sprintf cmd, "do shell script \"ls %s\"",pathStr
//	ExecuteScriptText cmd
//	print cmd
//	print S_value[0,72]
//End

//  ========================= End of legacy deprecated functions =========================  //
//  ======================================================================================  //
