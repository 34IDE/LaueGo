#pragma rtGlobals=2		// Use modern global access method.
#pragma ModuleName=JZTutil
#pragma IgorVersion = 6.11
#pragma version = 3.23
#pragma hide = 1

Menu "Graph"
	"Append Multiple Graphs to a Layout",AppendGraph2LayoutN(NaN,"","")
End
Menu "Layout"
	"Append Multiple Graphs to a Layout",AppendGraph2LayoutN(NaN,"","")
End
Menu "Data"
	"Generic Wave Note Info",/Q,GenericWaveNoteInfo($"","")
End

// Sections:
//	1	function for optionally showing menus, e.g. MenuItemIfWaveClassExists(), and others
//	2	macros to put labels on lower left/right of plots to identify where they came from
//	3	tool for putting multiple graph on one page
//	4	String ranges, deals with "1,4,5-20",  for handling a non-consecutive range of integers, This one is good
//	5	Progress panels
//	6	Generic XML support
//	7	Contains lots of utility stuff
//		returns file type (using the $filetype) in my old standard files (trying to start only using xml)
//		keyInList() & MergeKeywordLists(), "key=value" list utilities
//		WaveListClass() & WaveInClass(), used for handling the "waveClass=ccccc;" in wave notes
//		  also AddClassToWaveNote(), ExcludeWavesInClass()
//		WavesWithMatchingKeyVals(), further filter a list of waves, look for those with matching key=value pairs
//		AxisLabelFromGraph(), gets the axis label
//		FindGraphsWithWave() & FindGizmoWithWave(), finds an existing graph or gizmo with a particular wave
//		OnlyWavesThatAreDisplayed(), removes waves that are not displayed from a list of wave
//		IncludeOnlyWavesInClass(), removes waves from the list if they are not of correct class
//		xy2saturatedColors(), computes saturated colors for and RGB wheel on an xy graph
//		reverseList(), reverses a list, handy for prompt,popups
//		monotonic(a), checks if a wave is monotonic
//		isdigit(c) & isletter(c), handy utilities
//		roundSignificant(val,N), returns val rounded to N places
//		placesOfPrecision(a), returns number of places of precision in a
//		ValErrStr(val,err), returns string  "val � err" formatted correctly
//		normalize(a), normalizes a if it is a vector or square matrix
//		FWHM of fitted peaks: GaussianFWHM(W_coef), LorentzianFWHM(W_coef)
//		Area of fitted peaks: LorentzianIntegral(W_coef), GaussianIntegral(W_coef), Gaussian2DIntegral(W_coef)
//		computeCOM(), compute the Center of Mass
//		PowerIntegerScale(), rescale a waves values by ^n or ^(1/n), preserves sign for non-complex values
//		Posix2HFS, a replacement for PosixToHFS(), (using ParseFilePath() for HFSToPosix()) we no longer need HFSAndPosix.xop
//		cpuFrequency(), systemUserName(), getEnvironment(), returns system info
//		TrimFrontBackWhiteSpace(str), TrimLeadingWhiteSpace(str), TrimTrailingWhiteSpace(str), trims whitespace
//		IgorFileTypeString() gives descriptive string from the NUMTYPE from WaveInfo()
//		GenericWaveNoteInfo(), returns wave note info
//		StopAllTimers(), stops all the Igor timers
//		dateStr2ISO8601Str(), convert a date to an ISO 8601 format
//		ISOtime2niceStr(iso), convert ISO8601 string to a nice format for graph annotations
//		ISOtime2IgorEpoch(iso), convert ISO8601 string to an Igor Epoch (error returns NaN)
//		vec2str(), convert a vector to a string
//		str2vec(), convert a string to a free vector
//		RomanNumeral(j) converts a number to a Roman Numeral string
//	7	Old legacy or deprecated functions



//  ============================================================================  //
//  =========================== Start of Option Menu Functions ===========================  //

// generic, used for lots of Menus
//Function/S MenuItemIfWaveClassExists(item,classes,options)
//	String item
//	String classes
//	String options
//	String list = WaveListClass(classes,"*",options)
//	return SelectString(strlen(list),"(","")+item
//End
Function/S MenuItemIfWaveClassExists(item,classes,optionsStr,[invisible])
	String item						// string that you want to appear in menu
	String classes					// semi-colon separated list of wave classes, can use "*" in classes
	String optionsStr				// options that are found in WaveList function optionsStr
	Variable invisible				// controls menu item when conditions not met: true -> menu item absent, false or absent -> menu item grayed out
	invisible = ParamIsDefault(invisible) ? 0 : invisible
	invisible = numtype(invisible) ? 0 : invisible
	String list = WaveListClass(classes,"*",optionsStr)
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

// for menus, enable menu item when a particular wave class is present on the graphName
Function/S MenuItemsWaveClassOnGraph(item,classes,graphName)
	String item					// item name to return, or "(item" to disable
	String classes					// list of classes to check for (semi-colon separated)
	String graphName				// name of graph, use "" for top graph

	Variable i
	String list = imageNameList("",";")	// first check the images
	for (i=0;i<ItemsInList(list);i+=1)
		Wave ww = ImageNameToWaveRef(graphName, StringFromList(i,list))
		if (WaveInClass(ww,classes))
			return item
		endif
	endfor

	list = TraceNameList("",";",3)		// second check ordinary traces and contours
	for (i=0;i<ItemsInList(list);i+=1)
		Wave ww = TraceNameToWaveRef(graphName, StringFromList(i,list))
		if (WaveInClass(ww,classes))
			return item
		endif
	endfor
	return "("+item						// top graph does not contain a wave of desired class
End

Function/S MenuItemIfWaveExists(item,wname)
	String item
	String wname

	Variable there = Exists(wname)==1
	return SelectString(there,"(","")+item
End

// This is really useful with an  Execute/P ... command
// e.g.  MenuItemIfWindowAbsent("Include ABC Support","ABC.ipf"), Execute/P "INSERTINCLUDE  \"ABC\"";Execute/P "COMPILEPROCEDURES "
Function/S MenuItemIfWindowAbsent(item,win)		// Shows menu item if win NOT present
	String item
	String win

	GetWindow/Z $win, hide
	return SelectString(V_flag,"(","")+item
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

//  =========================== End of Option Menu Functions ============================  //
//  ============================================================================  //


//  ============================================================================  //
//  ============================== Start of Corner Labels ==============================  //

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
	String win=WinName(0,127+4096,1)// Graph, Table, Layout, Notebook,Panel, or XOP window
	String gwin=""

	switch(WinType(win))
		case 3:		// a layout, find the first graph on this layout
			String list="x"
			Variable i
			for (i=0; strlen(gwin)==0 && strlen(list); i+=1)		// get first graph name
				list = LayoutInfo(win,num2istr(i))
				if (stringmatch(StringByKey("TYPE",list),"Graph"))
					gwin = StringByKey("NAME",list)
				endif
			endfor
			break
		case 1:		// a graph, gwin is win
		case 2:		// a table, gwin is win
		case 13:	// from an XOP, probably a gizmo
			gwin = win
			break
		default:		// only layout or graph
			gwin = ""
	endswitch
	gwin = SelectString(strlen(gwin),"",":")+gwin

	PathInfo home			// creates String s_path
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

//  =============================== End of Corner Labels ==============================  //
//  ============================================================================  //


//  ============================================================================  //
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
		Prompt gName, "Graph to add to layout",popup,WinList("*",";","WIN:1")
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
	if (strlen(LayoutInfo(lname, gName ))>1)
		return lname								// graph is on layout, all done
	endif
	// check if this layout is full
	Variable i,N,NL = NumberByKey("NUMOBJECTS", LayoutInfo(lname,"Layout"))
	for (i=0,N=0;i<NL;i+=1)						// check each object in the layout
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
//  ============================================================================  //


//  ============================================================================  //
//  ============================== Start of String Ranges ==============================  //

// This section is for dealing with random or usually non contiguous sequence of integers
// i.e.  you took data in scans 1-30, but scans 17, and 25 were no good.  So the valid range is "1-16,18-24,26-30"
//  or perhaps you want to combine scans "3,7,9"  The following routines handle those situations in a simple fashion.

ThreadSafe Function NextInRange(range,last)	// given a string like "2-5,7,9-12,50" get the next number in this compound range
									// the range is assumed to be monotonic, it returns NaN if no more values
	String range					// list defining the range
	Variable last					// last number obtained from this range, use -Inf to get start of range, it returns the next

	// find first item in the list that should use next
	String item
	Variable m,i,j
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

		i = strsearch(item,"-",1)		// location of first '-' after the first character
		if (i<0)							// only a single number, not a dash type range, keep looking
			j += 1
			continue
		endif
		// check to see if last was in the range of item, but not the last value
		if (last>=str2num(item) && last<str2num(item[i+1,Inf]))
			return last+1
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
	Variable len=0						// the result, number of values represented
	Variable m,i,j,N=ItemsInList(range,",")

	for (j=0;j<N;j+=1)				// loop over each item
		item = StringFromList(j,range,",")

		m=-1							// remove any leading white space from item
		do
			m += 1
		while (char2num(item[m])<=32)
		item = item[m,strlen(item)-1]

		i = strsearch(item,"-",1)		// location of first '-' after the first character
		if (i<0)							// only a single number, not a dash type range, keep looking
			len += 1
		else								// item is a dash type
			len += str2num(item[i+1,Inf])-str2num(item)+1
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
	Variable i,last

	i = strsearch(range,",",Inf,1)
	if (i>=0)									// remove previous comma separated items
		range = range[i+1,Inf]
	endif
	range = ReplaceString(" ",range,"")		// spaces do not count
	range = ReplaceString("\t",range,"")	// tabs do not count

	i = strsearch(range,"-",Inf,1)
	if (char2num(range[i-1])==45)		// a double dash, we have a minus sign
		i -= 1
	endif
	i += i>0 ? 1 : 0							// increment i to skip over dash unless there is no dash in string, i<0 means no dash present
	last = str2num(range[i,Inf])
	return last
End

ThreadSafe Function isInRange(range,m)// returns TRUE if m is a number in range
	String range							// list defining the range
	Variable m								// number of interest, is m part of range?

	String item							// each of the comma sepated items
	Variable i,j,N=ItemsInList(range,",")
	for (j=0;j<N;j+=1)					// loop over each comma separated item
		item = StringFromList(j,range,",")
		item = TrimLeadingWhiteSpace(item)
		i = strsearch(item,"-",1)		// location of first '-' after the first character
		if (i<0)							// only a single number, not a dash type range, keep looking
			if (m==str2num(item))
				return 1					// this number is m
			endif
		else								// item is a dash type
			if (m>=str2num(item) && m<=str2num(item[i+1,Inf]))
				return 1					//  m is in this continuous range
			endif
		endif
	endfor
	return 0
End

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
	Variable j=0
	do
		str = StringFromList(j, range, ",")
		Variable m=-1				// remove any leading white space
		do
			m += 1
		while (char2num(str[m])<=32)
		str = str[m,strlen(str)-1]

		// now check str to see if it is a range like "20-23"
		i1 = str2num(str)
		i = strsearch(str,"-",strlen(num2str(i1)))		// position of "-" after first number
		if (i>0)
			i2 = str2num(str[i+1,inf])
			i = i1
			do
				out += num2str(i)+sep
				i += 1
			while (i<=i2)
		else
			out += num2str(i1)+sep
		endif
		j += 1
	while (j<N)

	i = strlen(out)-1
	if (char2num(out[i])==char2num(sep))
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
	Variable j,first,last,i=0
	Variable N=ItemsInList(range,sep)
	if (N<1)
		return ""
	endif
	range = SortList(range,sep,2)				// make list monotonic
	last = str2num(StringFromList(0,range,sep))-2	// ensure that first item is at the start
	for (i=0;i<N;i+=1)
		j = str2num(StringFromList(i,range,sep))
		num = num2str(j)
		if (numtype(j))
			return ""
		elseif ((j-last)==1)					// keep counting
			last = j
		elseif ((j-last)!=1)					// new sub-range
			if (i==0)							// special for first point
				comp = num
			elseif (first==last)					// just add a single number range
				comp += ","+num
			else									// close out previous range, and add single number
				comp += "-"+num2str(last)+","+num
			endif
			last = j
			first = j
		endif
	endfor
	if (first!=last)
		comp += "-"+num2str(last)
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
//End

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

//  =============================== End of String Ranges ==============================  //
//  ============================================================================  //


//	===================================================================================
//	================================= Start of Progress Panel =================================

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
		Variable right = stop ? 695 : 630,  bottom=117, new=0
		bottom += showTime ? 20 : 0
		bottom += showStatus ? 20 : 0
		NewPanel/K=1/W=(330,87,right,bottom)/N=Progress
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

	Variable top=29
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

	if (WinType(wname)!=7)							// window not a panel
		return 1
	endif
	if (part>=0)											// part is valid, update it
		SetWindow $wname userdata(part)=num2str(part)
	endif
	part = numtype(part) ? 0 : part
	if (part<total)											// total is valid update it
		SetWindow $wname userdata(total)=num2str(total)
	endif
	return 0
End


//Function/T ProgressPanelStart(percentName,[stop,showTime])	// display a progress bar
////	if the user does a "DoUpdate/W=$wname", then V_flag will be set to 2 if the stop button was pushed
//	String percentName									// name of global value to use, if not a global variable then set value to 0
//	Variable stop											// if true, include a stop button
//	Variable showTime									// if true, include the ending time on panel
//	stop = ParamIsDefault(stop) ? 0 : stop
//	stop = NumberByKey("IGORVERS",IgorInfo(0))>=6.10 ? stop : 0	// stop feature requires version 6.10 or above
//	showTime = ParamIsDefault(showTime) ? 0 : showTime
//
//	Variable right = stop ? 695 : 630,  bottom = showTime ? 137 : 117
//	NewPanel/K=1/W=(330,87,right,bottom)/N=Progress
//	String wname = S_name								// save the name of the window
//	ValDisplay progressValdisp, win=$wname, pos={5,5},size={288,25},font="Helvetica",fSize=18,format="%d%%"
//	ValDisplay progressValdisp, win=$wname, limits={0,100,0},barmisc={0,50},highColor= (19368,32650,65196)
//	if (exists(percentName)==2)						// the global exists, use it
//		NVAR percent = $percentName
//		percent = limit(percent,0,100)
//		ValDisplay progressValdisp,win=$wname, value= #percentName
//	else													// no global, just set to zero
//		ValDisplay progressValdisp,win=$wname, value= #"0"
//	endif
//
//	if (stop)
//		Button stopButton,pos={299,4},size={50,20},title="Break"
//	endif
//	if (showTime)
//		TitleBox timeTitle,pos={1,29},size={213,16},fSize=12,frame=0,title=" "
//	endif
//#if NumberByKey("IGORVERS",IgorInfo(0))>=6.10
//	DoUpdate/W=$wname/E=1							// mark as a progress window
//#else
//	DoUpdate
//#endif
//
//	String str = num2istr(DateTime)
//	SetWindow $wname userdata(startTimeSec)=str
//	SetWindow $wname userdata(lastTimeSec)=str		// last time updated
//	SetWindow $wname userdata(showTime)=num2istr(showTime)
//	SetWindow $wname userdata(stop)=num2istr(stop)
//	SetWindow $wname userdata(percent)=num2str(percent)
//	return wname
//End
////
//Function ProgressPanelUpdate(wname,percent)		// update a progress bar
//	String wname											// window name
//	Variable percent										// percent done (% !!)
//
//	ValDisplay progressValdisp,win=$wname, value= #num2istr(floor(limit(percent,0,100)))
//	Variable showTime = str2num(GetUserData(wname,"","showTime"))
//	if (showTime && percent>0)
//		Variable start = str2num(GetUserData(wname,"","startTimeSec"))
//		Variable last = str2num(GetUserData(wname,"","lastTimeSec"))
//		Variable now = DateTime
//		if (now-last > 2)
//			Variable remain = (now-start) * (100-percent)/percent
//			String str
//			if (remain < 5*60)							// under 10 min
//				sprintf str, "should end in %s,   at    %s",Secs2Time(remain,5),Secs2Time(now+remain,1)
//			else
//				sprintf str, "should end in %s,   at    %s",Secs2Time(remain,4),Secs2Time(now+remain,0)
//			endif
//			str += SelectString(remain > 24*3600,"","  on  "+Secs2Date(now+remain,2))	// more than 1 day
//			TitleBox timeTitle,title=str
//			SetWindow $wname userdata(lastTimeSec)=num2istr(now)	// last time updated, is now
//		endif
//	endif
//
//	Variable done=0
//#if NumberByKey("IGORVERS",IgorInfo(0))>=6.10
//	DoUpdate/W=$wname									// update only progress window
//	done = V_Flag == 2									// stop button on progress window was pushed
//#else
//	DoUpdate												// update all windows, no stop button
//#endif
//	return done
//End
//

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
//	================================== End of Progress Panel =================================
//	===================================================================================


//  ============================================================================  //
//  ============================== Start of Generic XML ===============================  //
//
//	XML support	 (occurance optionally allows selecting the the occuranceth instance of xmltag), note vectors usually delimited by a space
//
//	XMLNodeList(buf)											returns a list with all top level nodes in buf
//	XMLtagContents(xmltag,buf,[occurance])					returns the contents of xmltag
//	XMLtagContents2List(xmltag,buf,[occurance,delimiters])	returns the contents of xmltag as a list, useful for vectors in the contents
//	XMLattibutes2KeyList(xmltag,buf)							return a list with all of the attribute value pairs for xmltag
//	XMLremoveComments(str)									remove all xml comments from str

ThreadSafe Function/T XMLNodeList(buf)				// returns a list of node names at top most level in buf
	String buf
	String name,nodes=""
	Variable i0=0, i1,i2
	do
		i0 = strsearch(buf,"<",i0)						// find start of a tag
		if (i0<0)
			break
		endif
		i1 = strsearch(buf," ",i0)						// find end of tag name using i1 or i2, end will be in i1
		i1 = i1<0 ? Inf : i1
		i2 = strsearch(buf,">",i0)
		i2 = i2<0 ? Inf : i2
		i1 = min(i1,i2)
		if (numtype(i1) || (i1-i0-1)<1)
			break
		endif
		name = ReplaceString(";",buf[i0+1,i1-1],"_")// name cannot contain semi-colons
		nodes += name+";"

		i2 = strsearch(buf,"</"+name+">",i0)			// find the closer for this tag, check for '</name>'
		if (i2<0)
			i0 = strsearch(buf,">",i1+1)				// no '</name>', just a simple node
		else
			i0 = i2 + strlen(name) + 3					// first character after '</name>'
		endif
	while(i0>0)
	return nodes
End


ThreadSafe Function/T XMLtagContents(xmltag,buf,[occurance])
	String xmltag
	String buf
	Variable occurance									// use 0 for first occurance, 1 for second, ...
	occurance = ParamIsDefault(occurance) ? 0 : occurance

	Variable i0,i1
	i0 = startOfxmltag(xmltag,buf,occurance)
	if (i0<0)
		return ""
	endif
	i0 = strsearch(buf,">",i0)							// character after '>' in intro
	if (i0<0)
		return ""
	endif
	i0 += 1												// start of contents

	i1 = strsearch(buf,"</"+xmltag+">",i0)-1		// character just before closing '<tag>'
	if (i1<i0 || i1<0)
		return ""
	endif

	return buf[i0,i1]
End


ThreadSafe Function/T XMLtagContents2List(xmltag,buf,[occurance,delimiters]) //reads a tag contensts and converts it to a list
	String xmltag
	String buf
	Variable occurance									// use 0 for first occurance, 1 for second, ...
	String delimiters									// characters that might be used for delimiters (NOT semi-colon), default is space, tab, cr, or nl = " \t\r\n"
	occurance = ParamIsDefault(occurance) ? 0 : occurance
	if (ParamIsDefault(delimiters) || strlen(delimiters)==0)
		delimiters = " \t\r\n"							// the usual white-space characters
	endif

	String str = XMLtagContents(xmltag,buf,occurance=occurance)
	str = ReplaceString(";",str,"_")					// cannot have any semi-colons in input string

	Variable i
	for (i=0;i<strlen(delimiters);i+=1)
		str = ReplaceString(delimiters[i],str,";")		// replace every occurance of a character in delimiters with a semi-colon
	endfor

	do
		str = ReplaceString(";;",str,";")				// replace all multiple semi-colons with a single semi-colon
	while(strsearch(str,";;",0)>=0)

	if (char2num(str[0])==char2num(";"))			// remove any leaing semi-colon
		str = str[1,Inf]
	endif
	return str
End


ThreadSafe Function/T XMLattibutes2KeyList(xmltag,buf,[occurance])// return a list with all of the attribute value pairs for xmltag
	String xmltag											// name of tag to find
	String buf												// buf containing xml
	Variable occurance									// use 0 for first occurance, 1 for second, ...
	occurance = ParamIsDefault(occurance) ? 0 : occurance

	Variable i0,i1
	i0 = startOfxmltag(xmltag,buf,occurance)
	if (i0<0)
		return ""
	endif
	i0 += strlen(xmltag)+2								// start of attributes
	i1 = strsearch(buf,">",i0)-1						// end of attributes
	if (i1<i0)
		return ""
	endif

	// parse buf into key=value pairs
	buf = buf[i0,i1]
	buf = ReplaceString("\t",buf," ")
	buf = ReplaceString("\r",buf," ")
	buf = ReplaceString("\n",buf," ")
	buf = TrimFrontBackWhiteSpace(buf)
	String key, value, keyVals=""
	i0 = 0
	do
		i1 = strsearch(buf,"=",i0,0)
		key = TrimFrontBackWhiteSpace(buf[i0,i1-1])
		i0 = strsearch(buf,"\"",i1,0)+1				// character after the first double quote around value
		i1 = strsearch(buf,"\"",i0,0)-1				// character before the second double quote around value
		value = buf[i0,i1]
		if (strlen(key)>0)
			keyVals = ReplaceStringByKey(key,keyVals,value,"=")
		endif
		i0 = strsearch(buf," ",i1,0)						// find space separator, set up for next key="val" pair
	while(i0>0 && strlen(key) && strlen(value))
	return keyVals
End


ThreadSafe Function/T XMLremoveComments(str)		// remove all xml comments from str
	String str
	Variable i0,i1
	do
		i0 = strsearch(str,"<!--",0)					// start of a comment
		i1 = strsearch(str,"-->",0)						// end of a comment
		if (i0<0 || i1<=i0)
			break
		endif
		str[i0,i1+2] = ""									// snip out comment
	while(1)
	return str
End
//
ThreadSafe Static Function startOfxmltag(xmltag,buf,occurance)	// returns the index into buf pointing to the start of xmltag
	String xmltag, buf
	Variable occurance									// use 0 for first occurance, 1 for second, ...

	Variable i0,i1, i, start
	for (i=0,i0=0;i<=occurance;i+=1)
		start = i0
		i0 = strsearch(buf,"<"+xmltag+" ",start)		// find start of a tag with attributes
		i1 = strsearch(buf,"<"+xmltag+">",start)		// find start of a tag without attributes
		i0 = i0<0 ? Inf : i0
		i1 = i1<0 ? Inf : i1
		i0 = min(i0,i1)
		i0 += (i<occurance) ? strlen(xmltag)+2 : 0	// for more, move starting point forward
	endfor
	i0 = numtype(i0) || i0<0 ? -1 : i0
	return i0
End

//  =============================== End of Generic XML ===============================  //
//  ============================================================================  //


//  ============================================================================  //
//  ========================= Start of some general utility functions ========================  //

// make a function getFiletype that goes to a file and determines it's filetype
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
ThreadSafe Function/S MergeKeywordLists(list0,list1,priority,keySepStr,listSepStr)
	String list0,list1
	Variable priority				// 0 or 1
	String keySepStr				// separates key and value, defaults to colon
	String listSepStr				// separates key value pairs, defaults to semicolon
	keySepStr = SelectString(strlen(keySepStr),":",keySepStr)	// default to colon
	listSepStr = SelectString(strlen(listSepStr),";",listSepStr)	// default to semicolon
	String item, key,value
	Variable i,N=ItemsInList(list1)
	for (i=0;i<N;i+=1)				// for each keyword=value pair in list1
		item = StringFromList(i,list1,listSepStr)
		key = StringFromList(0,item,keySepStr)
		value = StringFromList(1,item,keySepStr)
		if (keyInList(key,list0,keySepStr,listSepStr) && priority==0)
			continue				// skip because key already in list0, and list0 has priority
		endif
		list0 = ReplaceStringByKey(key,list0,value,keySepStr,listSepStr)
	endfor
	return list0
End


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
			CheckDisplayed/W=$win $name
			displayed = !(!V_flag)
		endif
		if (displayed && WaveInClass($name,waveClassList,all=all))
			out += name+";"
		endif
	endfor
	return out
End
//Function/T WaveListClass(waveClassList,search,options,[all])
//	String waveClassList				// a list of acceptable wave classes (semicolon separated)
//	String search						// same as first argument in WaveList()
//	String options						// same as last argument in WaveList()
//	Variable all						// when all is TRUE, then all of the classes in waveClassList must be present, not just one
//	all = ParamIsDefault(all) ? 0 : !(!all)
//	all = numtype(all) ? 0 : all
//
//	String in = WaveList(search,";",options), out=""
//	String name
//	Variable m
//	for (m=0, name=StringFromList(0,in); strlen(name); m+=1,name=StringFromList(m,in))
//		if (WaveInClass($name,waveClassList,all=all))
//			out += name+";"
//		endif
//	endfor
//	return out
//End
//Function/T WaveListClass(waveClassList,search,options)
//	String waveClassList				// a list of acceptable wave classes (semicolon separated)
//	String search						// same as first argument in WaveList()
//	String options						// same as last argument in WaveList()
//
//	String in = WaveList(search,";",options), out=""
//	String name
//	Variable m
//	for (m=0, name=StringFromList(0,in); strlen(name); m+=1,name=StringFromList(m,in))
//		if (WaveInClass($name,waveClassList))
//			out += name+";"
//		endif
//	endfor
//	return out
//End


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
			use = strlen(FindGraphsWithWave(ww)) > 0
			use = not ? !use : use
			out += SelectString(use,"",name+";")
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
//// returns true if any one of the classes of ww matches one of the classes in waveClassList
//// note that the items in waveClassList can have wild cards
//ThreadSafe Function WaveInClass(ww,waveClassList)
//	Wave ww						// Wave to check
//	String waveClassList			// a list of acceptable wave classes (semicolon separated)
//	if (!WaveExists(ww) || strlen(waveClassList)<1)
//		return 0
//	endif
//	String class = StringByKey("waveClass",note(ww),"=")	// class list stored in wave note (comma separated)
//	String wavClass, matchClass
//	Variable m, i
//	for (m=0;m<ItemsInList(waveClassList);m+=1)			// check each item in waveClassList
//		matchClass = StringFromList(m,waveClassList)
//		for (i=0;i<ItemsInList(class,",");i+=1)
//			wavClass = StringFromLIst(i,class,",")				// class item from the wave
//			if (stringmatch(wavClass, matchClass))			// note that matchClass can have wild cards
//				return 1
//			endif
//		endfor
//	endfor
//	return 0
//End


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
	gName = SelectString(strlen(gName),StringFromList(0,FindGraphsWithWave(w)),gName)
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



Function/S FindGraphsWithWave(w)	// find the graph window which contains the specified wave
	Wave w
	if (!WaveExists(w))
		return ""
	endif
	String name0=GetWavesDataFolder(w,2), out=""
	String win,wlist = WinList("*",";","WIN:1"), clist, cwin
	Variable i,m,Nm=ItemsInList(wlist)
	for (m=0;m<Nm;m+=1)
		win = StringFromList(m,wlist)
		CheckDisplayed/W=$win w
		if (V_flag>0)
			out += win+";"
		else
			clist = ChildWindowList(win)
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


Function/T FindGizmoWithWave(w)	// find the Gizmo which contains the specified wave
	Wave w
	if (!WaveExists(w) || exists("NewGizmo")!=4)
		return ""
	endif

	Execute "GetGizmo/Z gizmoNameList"
	String gName, objetAll, gList=StrVarOrDefault("S_GizmoNames","")
	Variable i, N=ItemsInList(gList)
	for (i=0;i<N;i+=1)
		gName=StringFromList(i,gList)				// check each gizmo that is displayed
		Execute "GetGizmo/Z/N="+gName+" objectList"
		objetAll=StrVarOrDefault("S_gizmoObjectList","")
		objetAll = ReplaceString("=",objetAll,";")
		objetAll = ReplaceString(",",objetAll,";")
		objetAll = ReplaceString("}",objetAll,";")
		objetAll = ReplaceString("{",objetAll,";")
		objetAll = ReplaceString("\r",objetAll,";")
		objetAll = ReplaceString("\n",objetAll,";")
		if (strsearch(objetAll,";"+GetWavesDataFolder(w,2)+";",0,2)>=0)
			break
		endif
	endfor
	KillStrings/Z S_GizmoNames, S_gizmoObjectList
	KillWaves/Z TW_gizmoObjectList
	return SelectString(i<N,"",gName)
End



Function/T FindTablesWithWave(w)	// find the table windows which contains the specified wave
	Wave w
	if (!WaveExists(w))
		return ""
	endif
	String name0=GetWavesDataFolder(w,2), out=""
	String win,wlist = WinList("*",";","WIN:2"), clist, cwin
	Variable i,m,Nm=ItemsInList(wlist)
	for (m=0;m<Nm;m+=1)
		win = StringFromList(m,wlist)
		CheckDisplayed/W=$win w
		if (V_flag>0)
			out += win+";"
		else
			clist = ChildWindowList(win)
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



Function DrawMarker(x0,y0,dx,dy,style,[color,thick,dash,win,layer])
	Variable x0,y0,dx,dy
	String style
	String color							// one of "red", "blue", "green", "yellow", "magenta", "cyan", or a triplet like "1,1,0" or "65535,0,0"
	Variable thick						// line thickness, defaults to 0.50 (thin)
	Variable dash
	String win							// optional name of window
	String layer						// Drawing layer to use.  Default is "UserFront"
	if (numtype(x0++y0+dx+dy) || dx<=0 || dy<=0)
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
	SetDrawLayer/W=$win $layer

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
	SetDrawLayer/W=$win UserFront
	return err
End



// make saturated color RGB from point (dx,dy) on circle of radius rmax. RGB scaled to rgbMax
Function xy2saturatedColors(dx,dy,rmax,rgbMax,rgb)
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



ThreadSafe Function monotonic(a)
	// determines whether values of a particular wave are monotonic increasing (if any NaN, returns false)
	Wave a		// the wave

	if (!WaveExists(a))
//		Abort "ERROR -- monotonic(), wave not found"
		print "ERROR -- monotonic(), wave not found"
		return NaN
	elseif (numpnts(a)<=0)
//		Abort "ERROR -- monotonic(), wave is empty"
		print "ERROR -- monotonic(), wave is empty"
		return NaN
	elseif (WaveDims(a)>1)
//		Abort "ERROR -- monotonic(), wave is not 1D"
		print "ERROR -- monotonic(), wave is not 1D"
		return NaN
	endif
	MatrixOP/O delta = a - rotateRows(a,1)
	delta[0] = 0
	WaveStats/M=1/Q delta
	return !(V_min<0 || V_numNans)
End
//	Variable n=numpnts(a)
//	Variable i=1
//	do
//		if (a[i-1]>a[i])
//			return 0
//		endif
//		i += 1
//	while (i<n)
//	return 1


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
	Variable val								// value
	Variable err								// err in value
	Variable sp								// optionally put spaces around the �
	sp = ParamIsDefault(sp) ? 0 : !(!sp)
	sp = numtype(sp) ? 0 : sp
	err = numtype(err) ? 0 : abs(err)

	Variable n = ceil(log(abs(val/err)))+1// number of places of precision to show
	n = numtype(n) ? 6  : n					// default precision is 6
	n = limit(n,1,15)

	String vfmt, evfmt, str
	vfmt = "%."+num2istr(n)+"g"			// format for value
	String pm = SelectString(sp,"�"," � ")
	evfmt = vfmt + pm + SelectString(n>=2,"%.1g","%.2g")

	if (err>0)
		sprintf str,evfmt,val,err				// have a valid err, show it
	else
		sprintf str,vfmt,val					// no valid err, just show value
	endif
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
//Function normalize(a)	// normalize a and return the initial magnitude
//	Wave a
//	Variable norm_a = norm(a)
//	a /= norm_a
//	return norm_a
//End



// These next three functions convert Igor Peak fitting parameters to useful numbers
ThreadSafe Function GaussianFWHM(W_coef)			// FWHM of a Gaussian Peak
	// y = K0+K1*exp(-((x-K2)/K3)^2)
	Wave W_coef
	if (WaveType(W_coef,1)==1)						// W_coef exists, and it is numeric
		return 2*W_coef[3]*sqrt(ln(2))
	endif
	return NaN
End


ThreadSafe Function LorentzianFWHM(W_coef)			// FWHM of a Lorentzian Peak
	// y = K0+K1/((x-K2)^2+K3)
	Wave W_coef
	if (WaveType(W_coef,1)==1)						// W_coef exists, and it is numeric
		return 2*sqrt(W_coef[3])
	endif
	return NaN
End


ThreadSafe Function LorentzianIntegral(W_coef)		// peak integral from a Lorentzian fit
	// y = K0+K1/((x-K2)^2+K3)
	Wave W_coef
	if (WaveType(W_coef,1)==1)						// W_coef exists, and it is numeric
		return PI*W_coef[1]/sqrt(W_coef[3])
	endif
	return NaN
End

ThreadSafe Function GaussianIntegral(W_coef)			// peak integral from a Gaussian fit
	// y = K0+K1*exp(-((x-K2)/K3)^2)
	Wave W_coef
	if (WaveType(W_coef,1)==1)						// W_coef exists, and it is numeric
		return W_coef[1] * abs(W_coef[3]) * sqrt(PI)
	endif
	return NaN
End

ThreadSafe Function Gaussian2DIntegral(W_coef)		// peak integral from a 2D-Gaussian fit
	// z = K0+K1*exp((-1/(2*(1-K6^2)))*(((x-K2)/K3)^2 + ((y-K4)/K5)^2 - (2*K6*(x-K2)*(y-K4)/(K3*K5))))
	Wave W_coef
	if (WaveType(W_coef,1)==1)						// W_coef exists, and it is numeric
		Variable K1=W_coef[1], K3=W_coef[3], K6=W_coef[6], K5=W_coef[5]
		return 2*PI * K1 * abs(K3*K5) * sqrt(1-K6^2)
	endif
	return NaN
End



ThreadSafe Function computeCOM(ywave,xwave)		// computes center of mass,  xwave is optional (use $"" if no xwave)
	Wave ywave,xwave
	Variable com=NaN
	if (WaveType(ywave,1)!=1)							// ywave must exist, and must be numeric
		return NaN
	endif
	MatrixOP/FREE yy = ReplaceNaNs(ywave,0)
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
	ExecuteScriptText/Z cmd
	if (V_flag)
		if (printIt)
			printf "failure in Posix2HFS(\"%s\")\t\tcmd = '%s'\r",posixName,cmd
			print S_value
		endif
		return ""										// there is no valid file path
	endif
	return S_value[1,strlen(S_value)-2]		// strip off leading & trailing double-quotes
End


Function cpuFrequency()		// return the cpu frequency (Hz)
	if (!stringmatch(StringByKey("OS",IgorInfo(3)),"Macintosh OS X"))
		print "Only know how to get cpu speed from a Mac"
		// DoAlert 0, "Only know how to get cpu speed from a Mac"
		return NaN								// cannot get answer
	endif
	String cmd
	sprintf cmd "do shell script \"sysctl hw.cpufrequency\""
	ExecuteScriptText cmd						//returns something like: 	"hw.cpufrequency: 2000000000"
	Variable freq = NumberByKey("hw.cpufrequency",ReplaceString("\"",S_value,""))
	if (numtype(freq))
		DoAlert 0, "Unable to get cpu frequency, message is '"+S_value+"'"
	endif
	return freq
End

Function/T systemUserName()		// return unix username
	if (!stringmatch(StringByKey("OS",IgorInfo(3)),"Macintosh OS X"))
		print "Only know how to get user name from a Mac"
		return ""								// cannot get answer
	endif
	String cmd
	sprintf cmd "do shell script \"whoami\""
	ExecuteScriptText cmd						//returns the user name in double quotes
	Variable i=strlen(S_value)-2
	String userName=S_value[1,i]
	if (i<=0)									// no name,
		DoAlert 0, "Unable to get user name, message is '"+S_value+"'"
	endif
	return userName
End
//Function/S ipAddress()			// return the IP address as a string
//	if (!stringmatch(StringByKey("OS",IgorInfo(3)),"Macintosh OS X"))
//		DoAlert 0, "Only know how to get hostname from a Mac"
//		return ""								// cannot get answer
//	endif
//	ExecuteScriptText "do shell script \"hostname\""						//returns something like: 	"tischler.uni.aps.anl.gov"
//	return  ReplaceString("\"",S_value,"")
//End


// Get a list of system environment variables, semi-colon separated
#if (stringmatch(IgorInfo(2),"Macintosh"))
Function/T getEnvironment()
	String cmd
	sprintf cmd "do shell script \"source ~/.bash_profile ; env\""
	ExecuteScriptText cmd
	Variable N=strlen(S_value)
	if (char2num(S_value[0]) == 34)			// strip leading double quote
		S_value = S_value[1,N-1]
		N -= 1
	endif
	if (char2num(S_value[N-1]) == 34)		// strip trailing double quote
		S_value = S_value[0,N-2]
	endif
	if (strsearch(S_value, ";", 0)>=0)
		DoAlert 0, "env has semi-colon"
		return ""
	endif
	S_value = ReplaceString("\r",S_value,";")
	return S_value
End
#else
Function/T getEnvVariables()
	return ""
End
#endif



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


Function/T dateStr2ISO8601Str(dateStr,[timeStr,zoneHr,zoneMin])	// convert input date string and converts to ISO8601 format
	String dateStr								// something like "4/2/2010" or "2010-4-2"
	String timeStr
	Variable zoneHr, zoneMin					// the time zone in hours or minutes, OPTIONAL, (you can also use hour=6.5)
	timeStr = SelectString(ParamIsDefault(timeStr),timeStr,"")

	Variable d1,d2,d3, dateSec
	dateStr = ReplaceString("-",dateStr,"/")
	sscanf dateStr,"%d/%d/%d",d1,d2,d3
	if (V_flag!=3)
		return ""								// failed to find date like string
	endif
	String outStr="", zoneStr=""
	if (d1>1903)								// assume year/month/day
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


Function ISOtime2IgorEpoch(iso)	// convert ISO8601 string to an Igor Epoch (error returns NaN)
	String iso							// format     2013-10-02T01:22:35  (the seconds are optional)

	Variable year=NaN,month=NaN,day=NaN, hr=NaN,mn=NaN,se=NaN
	sscanf iso,"%4d-%2d-%2dT%2d:%2d:%2d", year,month,day,hr,mn,se
	Variable N=V_flag

	if (N<3 || numtype(year+month+day))
		return NaN
	endif
	Variable epoch=date2secs(year, month, day )
	if (N>=5)
		epoch += hr*3600
		epoch += mn*60
		se = (N>=6) ? se : 0
		epoch += se
	endif
	return epoch
End


Function/T ISOtime2niceStr(iso)	// convert ISO8601 string to a nice format for graph annotations
	String iso							// format     2013-10-02T01:22:35  (the seconds are optional)

	Variable year=NaN,month=NaN,day=NaN, hr=NaN,mn=NaN,se=NaN
	sscanf iso,"%4d-%2d-%2dT%2d:%2d:%2d", year,month,day,hr,mn,se
	Variable N = V_flag

	if (N<3 || numtype(year+month+day))
		return ""
	endif
	Variable epoch = date2secs(year, month, day )
	String out = Secs2Date(epoch,2)
	if (N>=5)
		epoch += hr*3600
		epoch += mn*60
		se = (N>=6) ? se : 0
		epoch += se
		Variable fmt = (N>=6) ? 1 : 0
		out += SelectString(numtype(epoch),"  "+Secs2Time(epoch,fmt),"")
	endif
	return out
End


ThreadSafe Function/T vec2str(w,[places,maxPrint,bare,sep])		// convert vector to s string suitable for printing, does not include name
	Wave w										// 1d wave to print
	Variable places								// number of places, for default, use negative or NaN
	Variable maxPrint							// maximum number of elements to print, defaults to 20
	Variable bare								// if bare is TRUE, then suppress the "{}" in the output
	String sep									// optional separator, default is ",  "   a comma and 2 spaces

	maxPrint = ParamIsDefault(maxPrint) ? 20 : maxPrint
	maxPrint = maxPrint>0 ? maxPrint : 20
	places = ParamIsDefault(places) ? -1 : places
	bare = ParamIsDefault(bare) ? 0 : !(!bare)
	sep = SelectString(ParamIsDefault(sep),sep,",  ")

	if (!WaveExists(w))
		return SelectString(bare,"{}","")
	endif

	Wave/T tw=$GetWavesDataFolder(w,2)
	Wave/C cw=$GetWavesDataFolder(w,2)
	Variable waveIsComplex = WaveType(w) %& 0x01
	Variable numeric = (WaveType(w)!=0)

	String fmt
	if (waveIsComplex)
		places = places>=0 ? min(20,places) : 5	// default to 5 for unacceptable values
		sprintf fmt,"(%%.%dg, %%.%dg)",places,places
	elseif (numeric)
		places = places>=0 ? min(20,places) : 5	// default to 5 for unacceptable values
		sprintf fmt,"%%.%dg",places
	elseif (places>0)								// must be string, and a maximum length given
		sprintf fmt, "\"%d%%s\"",places
	else												// string with no preferred length
		fmt = "\"%%s\""
	endif

	Variable i=0, n
	n = numpnts(w)
	maxPrint = min(n,maxPrint)
	String str, out=SelectString(bare,"{","")

	do
		if (waveIsComplex)						// a complex wave
			sprintf str,fmt, real(cw[i]),imag(cw[i])
		elseif (numeric && (!waveIsComplex))	// a simple number wave
			sprintf str,fmt, w[i]
		elseif (!numeric)						// a text wave
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


Function SIprefix2factor(prefix)
	String prefix

	Variable i, value=1
	String ch
	for (i=0;i<strlen(prefix);i+=1)
		ch = prefix[i]

		if (strlen(prefix)==0)								// no prefix
			value *= 1
		elseif (strsearch(ch,"d",0)==0)			// deci
			value *= 1e-1
		elseif (strsearch(ch,"c",0)==0)			// centi
			value *= 1e-2
		elseif (strsearch(ch,"m",0)==0)			// milli
			value *= 1e-3
		elseif (strsearch(ch,"�",0)==0)			// micro
			value *= 1e-6
		elseif (strsearch(ch,"n",0)==0)			// nano
			value *= 1e-9
		elseif (strsearch(ch,"p",0)==0)			// pico
			value *= 1e-12
		elseif (strsearch(ch,"f",0)==0)			// femto
			value *= 1e-15
		elseif (strsearch(ch,"a",0)==0)			// atto
			value *= 1e-18
		elseif (strsearch(ch,"z",0)==0)			// zepto
			value *= 1e-21
		elseif (strsearch(ch,"y",0)==0)			// yocto
			value *= 1e-24

		elseif (strsearch(ch,"h",0,2)==0)		// hecto (h or H) is acceptable)
			value *= 100
		elseif (strsearch(ch,"k",0,2)==0)		// kilo (k or K is acceptable)
			value *= 1e3
		elseif (strsearch(ch,"M",0)==0)			// Mega
			value *= 1e6
		elseif (strsearch(ch,"G",0)==0)			// Giga
			value *= 1e9
		elseif (strsearch(ch,"T",0)==0)			// Tera
			value *= 1e12
		elseif (strsearch(ch,"P",0)==0)			// Peta
			value *= 1e15
		elseif (strsearch(ch,"E",0)==0)			// Exa
			value *= 1e18
		elseif (strsearch(ch,"Z",0)==0)			// Zeta
			value *= 1e21
		elseif (strsearch(ch,"Y",0)==0)			// Yotta
			value *= 1e24

		else												// unknown prefix char
			value = NaN
		endif
	endfor

	return value
End


Function/T RomanNumeral(j)	// convert integer j to a Roman Numeral String
	Variable j
	j = round(j)
	if (j<1 || j>3999 || numtype(j))
		return ""
	elseif (j>=1000)
		return "M"+RomanNumeral(j-1000)
	elseif (j>=900)
		return "CM"+RomanNumeral(j-900)
	elseif (j>=500)
		return "D"+RomanNumeral(j-500)
	elseif (j>=400)
		return "CD"+RomanNumeral(j-400)
	elseif (j>=100)
		return "C"+RomanNumeral(j-100)
	elseif (j>=90)
		return "XC"+RomanNumeral(j-90)
	elseif (j>=50)
		return "L"+RomanNumeral(j-50)
	elseif (j>=40)
		return "XL"+RomanNumeral(j-40)
	elseif (j>=10)
		return "X"+RomanNumeral(j-10)
	elseif (j>=9)
		return "IX"+RomanNumeral(j-9)
	elseif (j>=5)
		return "V"+RomanNumeral(j-5)
	elseif (j>=4)
		return "IV"+RomanNumeral(j-4)
	elseif (j>=1)
		return "I"+RomanNumeral(j-1)
	endif
	return ""
End

Function RomanNumeral2Int(str)	// convert a Roman Numeral string to an integer
	String str
	if (strlen(str)<1)
		return NaN
	endif
	str = UpperStr(str )

	Variable iM=char2num("M"),iC=char2num("C"),iX=char2num("X"),iI=char2num("I")
	Variable j

	for (;char2num(str[0])==iM;)
		str = str[1,Inf]
		j += 1000
	endfor

	if (strsearch(str,"CM",0)==0)
		str = str[2,Inf]
		j += 900
	endif
	if (strsearch(str,"D",0)==0)
		str = str[1,Inf]
		j += 500
	endif
	if (strsearch(str,"CD",0)==0)
		str = str[2,Inf]
		j += 400
	endif
	for (;char2num(str[0])==iC;)
		str = str[1,Inf]
		j += 100
	endfor

	if (strsearch(str,"XC",0)==0)
		str = str[2,Inf]
		j += 90
	endif
	if (strsearch(str,"L",0)==0)
		str = str[1,Inf]
		j += 50
	endif
	if (strsearch(str,"XL",0)==0)
		str = str[2,Inf]
		j += 40
	endif
	for (;char2num(str[0])==iX;)
		str = str[1,Inf]
		j += 10
	endfor

	if (strsearch(str,"IX",0)==0)
		str = str[2,Inf]
		j += 9
	endif
	if (strsearch(str,"V",0)==0)
		str = str[1,Inf]
		j += 5
	endif
	if (strsearch(str,"IV",0)==0)
		str = str[2,Inf]
		j += 4
	endif
	for (;char2num(str[0])==iI;)
		str = str[1,Inf]
		j += 1
	endfor

	j = strlen(str) ? NaN : j
	return j
End

//  ========================= End of some general utility functions =========================  //
//  ============================================================================  //


//  ============================================================================  //
//  ========================== Start of legacy deprecated functions =========================  //

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
//	return (ItemsInList(GetRTStackInfo(0))<3)
End

//Function abc()
//	String cmd, pathStr="\\\"/Users/tischler/data/Sector 34/July 4, 2006/EW5/recon/\\\""
//	print PathStr
//	sprintf cmd, "do shell script \"ls %s\"",pathStr
//	ExecuteScriptText cmd
//	print cmd
//	print S_value[0,72]
//End

//  ========================== End of legacy deprecated functions ==========================  //
//  ============================================================================  //
