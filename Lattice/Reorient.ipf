#pragma rtGlobals=1		// Use modern global access method.
#pragma version = 2.01
#pragma ModuleName=Reorient
#include "microGeometryN", version>=1.28
#include "IndexingN", version>=3.69


Menu "BackReflection"
	"Load Image File...", LoadGenericImageFile("")
	help={"Load an Image from the file, and then display the image"}
	MenuItemIfWaveClassExists("Re-Display image plot with buttons...","tif*;rawImage*","DIMS:2"), NewImageGraph($"",1)
	help={"Re-Display an image and put buttons on the plot for listing spots"}
		MenuItemIfWaveClassExists("   Find Z range of image","*","DIMS:2"), getZrange($"",NaN)
		help={"for an image, find the z range for the specified % range"}
	"-"
	"Reorient Crystal", ReorientCrystal()
End


//Static Function/S MarqueeReorientMenuItem(item)
//	String item
//	if (strlen(WinList("*","","WIN:1"))<1)
//		return "("+item
//	endif
//	Variable V_flag
//	GetMarquee/Z								// These menu items requre a Marquee present
//	if (V_flag==0)
//		return "("+item
//	endif
//	Wave image = ImageNameToWaveRef("",StringFromList(0,ImageNameList("",";")))
//	if (!WaveExists(image))
//		return "("+item
//	endif
//	String class = StringByKey("waveClass", note(image),"=")
//	if (!stringmatch(class,"spe*"))				// is not an editable mas, disable (only works for spe images
//		item = "("+item
//		// item = "("
//	endif
//	return item
//End



Function ReorientPanel() : Panel
	NewPanel /W=(765,206,1549,909)/K=1
	SetDrawLayer UserBack
	DrawText 65,628,"Axis along beam: positive rotates sample CCW when viewed in beam direction"
	DrawText 67,642,"Horizontal axis: positive tilts sample face upward"
	DrawText 66,659,"Vertical axis: positive tilts sample toward the left side of the image"
	GroupBox group0,pos={41,25},size={237,544},title="\\JC Current Orientation"
	GroupBox group1,pos={328,31},size={216,537},title="\\JCNew Orientation"
	SetVariable CurrentBeamH,pos={76,186},size={107,16},disable=2,title="H = "
	SetVariable CurrentBeamH,value= CurrentInvRL[0][2]
	SetVariable CurrentBeamL,pos={74,237},size={107,16},disable=2,title="L = "
	SetVariable CurrentBeamL,value= CurrentInvRL[2][2]
	SetVariable CurrentBeamK,pos={73,212},size={107,16},disable=2,title="K = "
	SetVariable CurrentBeamK,value= CurrentInvRL[1][2]
	GroupBox group2,pos={62,167},size={142,101},title="Beam Direction"
	GroupBox group3,pos={65,283},size={144,114},title="Image Up Direction"
	SetVariable CurrentUpH,pos={80,313},size={107,16},disable=2,title="H = "
	SetVariable CurrentUpH,value= CurrentInvRL[0][1]
	SetVariable CurrentUpK,pos={79,340},size={107,16},disable=2,title="K = "
	SetVariable CurrentUpK,value= CurrentInvRL[1][1]
	SetVariable CurrentUpL,pos={79,369},size={107,16},disable=2,title="L = "
	SetVariable CurrentUpL,value= CurrentInvRL[2][1]
	GroupBox group4,pos={65,413},size={144,114},title="Image Right Direction"
	SetVariable CurrentRightH,pos={80,443},size={107,16},disable=2,title="H = "
	SetVariable CurrentRightH,value= CurrentInvRL[0]
	SetVariable CurrentRightK,pos={79,470},size={107,16},disable=2,title="K = "
	SetVariable CurrentRightK,value= CurrentInvRL[1]
	SetVariable CurrentRightL,pos={80,495},size={107,16},disable=2,title="L = "
	SetVariable CurrentRightL,value= CurrentInvRL[2]
	GroupBox group5,pos={58,46},size={167,94},title="Goniometer Settings"
	SetVariable CurrentGonioOuter,pos={78,73},size={107,16},title="Outer = "
	SetVariable CurrentGonioOuter,value= CurrentGonio[1]
	SetVariable CurrentGonioInner,pos={76,103},size={118,16},title="Inner = "
	SetVariable CurrentGonioInner,value= CurrentGonio[0]
	SetVariable NewBeamH,pos={375,202},size={107,16},title="H = "
	SetVariable NewBeamH,value= NewInvRL[0][2]
	SetVariable NewBeamK,pos={372,228},size={107,16},title="K = "
	SetVariable NewBeamK,value= NewInvRL[1][2]
	SetVariable NewBeamL,pos={373,253},size={107,16},title="L = "
	SetVariable NewBeamL,value= NewInvRL[2][2]
	GroupBox group6,pos={361,183},size={142,101},title="Beam Direction"
	GroupBox group7,pos={364,299},size={144,114},title="Image Up Direction"
	SetVariable NewUpH1,pos={379,329},size={107,16},disable=2,title="H = "
	SetVariable NewUpH1,value= NewInvRL[0][1]
	SetVariable NewUpK1,pos={378,356},size={107,16},disable=2,title="K = "
	SetVariable NewUpK1,value= NewInvRL[1][1]
	SetVariable NewUpL1,pos={379,381},size={107,16},disable=2,title="L = "
	SetVariable NewUpL1,value= NewInvRL[2][1]
	GroupBox group8,pos={364,429},size={144,114},title="Image Right Direction"
	SetVariable NewRightH,pos={379,459},size={107,16},disable=2,title="H = "
	SetVariable NewRightH,value= NewInvRL[0]
	SetVariable NewRight,pos={378,486},size={107,16},disable=2,title="K = "
	SetVariable NewRight,value= NewInvRL[1]
	SetVariable NewRightL,pos={379,511},size={107,16},disable=2,title="L = "
	SetVariable NewRightL,value= NewInvRL[2]
	GroupBox group9,pos={357,62},size={167,94},title="Goniometer Settings"
	SetVariable NewGonioOuter,pos={377,89},size={128,16},title="Outer = "
	SetVariable NewGonioOuter,value= NewGonio[1]
	SetVariable NewGonioInner,pos={375,119},size={131,16},title="Inner = "
	SetVariable NewGonioInner,value= NewGonio[0]
	Button ButtonCalcGonio,pos={581,48},size={114,50},proc=Reorient#ButtonCalcGonio,title="Calculate New\rGoniometer Settings"
	Button ButtonCalcRLV,pos={582,114},size={114,50},proc=Reorient#ButtonCalcRLV,title="Calculate New\rReciprocal Lattice\rOrientation"
	Button ExitButton,pos={612,190},size={50,20},proc=Reorient#ButtonExit,title="Exit"
	GroupBox GonioConv,pos={56,588},size={456,89},title="Goniometer Conventions"
	GroupBox group10,pos={559,225},size={222,322},title="Goniometer Choice"
	CheckBox RadioBlake,pos={587,253},size={116,39},proc=Reorient#GonioRadioButtons,title="Blue Blake\rInner axis: Horizontal\rOuter axis: Vertical"
	CheckBox RadioBlake,value= 1,mode=1
	CheckBox RadioSouthBayHorizontal,pos={588,312},size={124,39},proc=Reorient#GonioRadioButtons,title="Black SouthBay\rInner axis: along beam\rOuter axis: horizontal"
	CheckBox RadioSouthBayHorizontal,value= 0,mode=1
	CheckBox RadioSouthBayVertical,pos={591,370},size={124,39},proc=Reorient#GonioRadioButtons,title="Black SouthBay\rInner axis: along beam\rOuter axis: vertical"
	CheckBox RadioSouthBayVertical,value= 0,mode=1
EndMacro

// Entry function for this procedure
Function ReorientCrystal()
	ReorientCrystalAssignGlobals()

	Wave FullPeakIndexed								// provides the reciprocal lattice
	Wave CurrentInvRL
	Variable pattern									// in case more than one, default is 0
	SVAR GonioType = root:GonioType
	
	String wList=WaveListClass("IndexedPeakList*","*","")
	String wName=StringFromList(0,wList)
	if (ItemsInList(wList)!=1)						// only ask if more than one choice
		Prompt wName, "name of Full Peak Indexed list",popup,WaveListClass("IndexedPeakList*","*","")
		DoPrompt "Full Peak Indexed",wName
		if (V_flag)
			return 1
		endif
	endif
	Wave FullPeakIndexed=$wName
	
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
		pattern -= 1
	endif
	
	// read Reciprocal Lattice matrix
	Make/N=(3,3)/O/D RL_reorient	
	Wave RL=RL_reorient
	Variable as0,as1,as2,bs0,bs1,bs2,cs0,cs1,cs2
	String wnote = note(FullPeakIndexed)
	sscanf StringByKey("recip_lattice"+num2istr(pattern),wnote,"="), "{{%g,%g,%g}{%g,%g,%g}{%g,%g,%g}}",as0,as1,as2,bs0,bs1,bs2,cs0,cs1,cs2
	RL[0][0]= {as0,as1,as2}								// the original measured RL
	RL[0][1]= {bs0,bs1,bs2}
	RL[0][2]= {cs0,cs1,cs2}
	
	MatrixInverse RL
	Wave M_Inverse
	CurrentInvRL = M_Inverse/(MatrixDet(M_Inverse))^(1/3)
	DisplayReciprocalLattice() // transfer control to dialog box
	CheckBox RadioBlake, value= stringmatch(GonioType ,"Blake")
	CheckBox RadioSouthBayHorizontal,value= stringmatch(GonioType, "SouthBayHorizontal")
	CheckBox RadioSouthBayVertical,value= stringmatch(GonioType, "SouthBayVertical")
End

// Load the panel; further action is controlled by the buttons
Function DisplayReciprocalLattice()
	if ( strlen( WinList("ReorientPanel",";","" )) > 0 )
		KillWindow ReorientPanel
	endif
	Execute "ReorientPanel()"
End	

Static Function ReorientCrystalAssignGlobals()

	if( ! WaveExists(CurrentInvRL))
		make /N=(3,3) CurrentInvRL
	endif

	if( ! WaveExists(NewInvRL))
		make /N=(3,3) NewInvRL
	endif
	
	if( ! WaveExists(HKLNew))
		make /N=(3) HKLNew
	endif

	if( ! WaveExists(NewBeam))
		make /N=(3) NewBeam
	endif

	if( ! WaveExists(CurrentGonio))
		make /N=(2) CurrentGonio
	endif

	if( ! WaveExists(NewGonio))
		make /N=(2) NewGonio
	endif

	if( ! WaveExists(CurrentGonRot))
		make /N=(3,3) CurrentGonRot
	endif

	if( ! WaveExists(NewGonRot))
		make /N=(3,3) NewGonRot
	endif

	if( ! WaveExists(GonRotIn))
		make /N=(3,3) GonRotIn
	endif

	if( ! WaveExists(GonRotOut))
		make /N=(3,3) GonRotOut
	endif

	SVAR /Z GonioType = root:GonioType
	if( ! SVAR_Exists(GonioType))
		String /G root:GonioType
	endif

End

Static Function ButtonCalcGonio(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			// click code here
			CalcNewGonioAngles()
			break
	endswitch

	return 0
End

Static Function ButtonCalcRLV(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			// click code here
			CalcNewRLV()
			break
	endswitch

	return 0
End

Static Function ButtonExit(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			// click code here
			KillWindow ReorientPanel
			break
	endswitch

	return 0
End

Static Function GonioRadioButtons(RadioButtonName,checked) : CheckBoxControl
	String RadioButtonName
	Variable checked

	SVAR GonioType = root:GonioType
		
	strswitch (RadioButtonName)
		case "RadioBlake":
			GonioType = "Blake"
			break
		case "RadioSouthBayHorizontal":
			GonioType = "SouthBayHorizontal"
			break
		case "RadioSouthBayVertical":
			GonioType = "SouthBayVertical"
			break
	endswitch
	CheckBox RadioBlake, value= stringmatch(GonioType ,"Blake")
	CheckBox RadioSouthBayHorizontal,value= stringmatch(GonioType, "SouthBayHorizontal")
	CheckBox RadioSouthBayVertical,value= stringmatch(GonioType, "SouthBayVertical")
End

// calculate orientation for specified goniometer angles
Static Function CalcNewRLV()
	Wave CurrentInvRL
	Wave CurrentGonRot
	Wave NewGonRot
	CalcGonRot(CurrentGonRot,CurrentGonio)
	CalcGonRot(NewGonRot,NewGonio)
	MatrixOp /O NewInvRL = CurrentInvRL x Inv(CurrentGonRot) x NewGonRot
End

// calculate rotation matrix GonRot corresponding to goniometer angles GonAng
Static Function CalcGonRot(GonRot,GonAng)
	Wave GonRot
	Wave GonAng
	Wave GonRotIn
	Wave GonRotOut	
	SVAR GonioType = root:GonioType

	strswitch (GonioType)
		case "SouthBayHorizontal":
			ZRot(GonRotIn, GonAng[0])
			XRot(GonRotOut, GonAng[1])
			break
		case "Blake":
			XRot(GonRotIn, GonAng[0])
			YRot(GonRotOut, GonAng[1])
			break
		case "SouthBayVertical":
			ZRot(GonRotIn, GonAng[0])
			YRot(GonRotOut, GonAng[1])
			break
	endswitch
	MatrixOp /O GonRot = GonRotIn x GonRotOut
End

Static Function XRot(RotMat, RotAng)
	Wave RotMat
	Variable RotAng
	
	RotMat = 0.0
	RotMat[0][0] = 1.0
	RotMat[1][1] = cos(RotAng*Pi/180)
	RotMat[2][2] = cos(RotAng*Pi/180)
	RotMat[1][2] = -sin(RotAng*Pi/180)
	RotMat[2][1] = sin(RotAng*Pi/180)
End

Static Function YRot(RotMat, RotAng)
	Wave RotMat
	Variable RotAng
	
	RotMat = 0.0
	RotMat[1][1] = 1.0
	RotMat[2][2] = cos(RotAng*Pi/180)
	RotMat[0][0] = cos(RotAng*Pi/180)
	RotMat[2][0] = -sin(RotAng*Pi/180)
	RotMat[0][2] = sin(RotAng*Pi/180)
End

Static Function ZRot(RotMat, RotAng)
	Wave RotMat
	Variable RotAng
		
	RotMat = 0.0
	RotMat[2][2] = 1.0
	RotMat[0][0] = cos(RotAng*Pi/180)
	RotMat[1][1] = cos(RotAng*Pi/180)
	RotMat[0][1] = -sin(RotAng*Pi/180)
	RotMat[1][0] = sin(RotAng*Pi/180)
End

Static Function CalcNewGonioAngles()
	Wave NewInvRL
	Wave CurrentGonRot
	Wave NewGonio
	Wave HKLNew
	Wave NewBeam
	Wave RL = RL_reorient
	SVAR GonioType = root:GonioType
	
	CalcGonRot(CurrentGonRot,CurrentGonio)
	
	HKLNew[0] = NewInvRL[0][2]
	HKLNew[1] = NewInvRL[1][2]
	HKLNew[2] = NewInvRL[2][2]
	
// NewBeam is the direction which must be rotated into the beam
	MatrixOp /O NewBeam = Normalize(CurrentGonRot x RL x HKLNew)
	if (NewBeam[3] < 0.0)		// take plane normal which is closer to current beam direction
		NewBeam = -NewBeam
	endif
	
// Solve the goniometer equation
	strswitch (GonioType)
		case "SouthBayHorizontal":
			NewGonio[1] = acos(NewBeam[2])	// outer tilt
			if (sin(NewGonio[1]) == 0.0)
				NewGonio[0] = 0.0
			else
				NewGonio[0] = atan2( NewBeam[0]/sin(NewGonio[1]), -NewBeam[1]/sin(NewGonio[1]) )
			endif
			WhichPhi()
			break
		case "Blake":
			NewGonio[1] = asin(NewBeam[0]) // outer tilt
			If(cos(NewGonio[1]) == 0.0)
				NewGonio[0] = 0.0
			else
				NewGonio[0] = atan2( -NewBeam[1]/cos(NewBeam[0]), NewBeam[2]/cos(NewBeam[0]) )
			endif
			break
		case "SouthBayVertical":
			NewGonio[1] = acos(NewBeam[2])	// outer tilt
			if (sin(NewGonio[1]) == 0.0)
				NewGonio[0] = 0.0
			else
				NewGonio[0] = atan2( NewBeam[1]/sin(NewGonio[1]), NewBeam[0]/sin(NewGonio[1]) )
			endif
			WhichPhi()
			break
	endswitch	
	NewGonio = NewGonio * 180/Pi
End

// When the inner rotation is along the beam, there are two setting which give the same orientation.
// WhichPhi picks the one closest to the current settings
Static function WhichPhi()
	Wave CurrentGonio
	Wave NewGonio
	Variable CurrentPhi = CurrentGonio[0]*Pi/180
	Variable PhiDiff = abs(mod(NewGonio[0] - CurrentPhi, 2*Pi))
	if( PhiDiff > Pi/2 && PhiDiff < 3*Pi/2 )
		NewGonio[0] = NewGonio[0] + Pi
		NewGonio[1] = -NewGonio[1]
	endif
	NewGonio[0] = mod(NewGonio[0],2*Pi)
	if(NewGonio[0] < 0)
		NewGonio[0] = NewGonio[0] + 2*Pi
	endif

//	do
//		if( NewGonio[0]-CurrentPhi < Pi )
//			break
//		endif
//		NewGonio[0] = NewGonio[0] - 2*Pi
//	while(1)
//	do
//		if( NewGonio[0]-CurrentPhi> -Pi )
//			break
//		endif
//		NewGonio[0] = NewGonio[0] + 2*Pi
//	while(1)
End