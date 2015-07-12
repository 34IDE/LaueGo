#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#pragma ModuleName=GizmoMovies
#pragma version = 2.05
#include "GizmoUtility", version>=0.16

Static Constant MAX_MOVIE_STEPS = 50		// maximum number of steps in a movie process (not max number of frames)
Static StrConstant MovieStepTypesBase = " ;Begin;Quaternion;Rotate;Clip;ClipReset;Add_Static"
Static StrConstant directionTypes = "+X;-X;+Y;-Y;+Z;-Z"
Static StrConstant GizmoMovieClipPlaneGroupName = "gizmoClipPlaneGroupMovie"



Static Function AfterFileOpenHook(refNum,file,pathName,type,creator,kind)
	Variable refNum, kind
	String file,pathName,type,creator
	if ((kind==1) || (kind==2))		// an experiment (packed or unpacked)
		GizmoMovies#InitGizmoMovies()
	endif
	return 0
End
Static Function IgorStartOrNewHook(IgorApplicationNameStr)
	String IgorApplicationNameStr
	GizmoMovies#InitGizmoMovies()
	return 0
End


#if (IgorVersion()>=7)
Menu "Gizmo"
	"Movie Steps Panel...", MakeMovieStepsPanel()
End
#endif


Menu "Gizmo-Movie"
	"Gizmo Movie Steps Panel...", MakeMovieStepsPanel()
	MenuItemIfWaveClassExists("Fill Panel from a Wave","MovieStepsList*",""), FillMovieStepsPanelFromWave($"")
	MenuItemIfWinExists("Save Panel to a Wave","MovieStepsPanel","WIN:64,VISIBLE:1"), SaveMovieSteps2Wave()
	"<I  Test Data & Gizmo", GizmoMovies#MakeTestGizmoData(); GizmoMovies#MakeGizmoMovieTestGizmo(); MakeMovieStepsPanel()
End


//  ======================================================================================  //
//  =========================== Start of Make a Movie Section ============================  //

Structure gizmoMovieStruct
	char		gizmoName[100]	// name of gizmo being used
	char		imageName[100]	// name of temporary image being used
	char		graphName[100]	// name of graph displaying imageName[]
	char		clipPlaneGroupName[100]	// name of clip plane group, "" means none on Gizmo
	double	quaternion[4]	// current quaternion, starts as quaternion of base orientation, but is continuously updated 
	int32		NframesAll		// total number of frames in movie, count those before Begin too.
	int32		Nframes			// total number of frames in movie, does not count anything before "Begin"
	int32		m					// number of frames so far written to the movie file
EndStructure



Function MakeMovieOfGizmo(gizName,movieSteps,[printIt,testing])
	String gizName
	Wave/T movieSteps
	Variable printIt
	Variable testing
	printIt = ParamIsDefault(printIt) ? NaN : printIt
	printIt = numtype(printIt) ? strlen(GetRTStackInfo(2))==0 : !(!printIt)
	testing = ParamIsDefault(testing) ? 0 : testing

	KillVariables/Z MovieIsOpenGlobal
	Variable MovieIsOpen=0

	Variable Nframes=NaN, NframesAll=0, fromPanel=0
	if (WaveExists(movieSteps))		// first choice is from passed wave
		Nframes = NumberByKey("Nframes",note(movieSteps),"=")
		NframesAll = NframesFromMovieStepListWave(movieSteps,all=1)
	endif
	if (!(NframesAll>0))				// second choice is from panel
		Wave/T movieSteps = GatherMovieStepFromPanel()
		Nframes = NframesFromMovieStepListWave(movieSteps)
		NframesAll = NframesFromMovieStepListWave(movieSteps,all=1)
		fromPanel = 1
	endif

	String GizmoNamesList = WinList("*",";","WIN:"+num2istr(GIZMO_WIN_BIT))
	if (WhichListItem(gizName,GizmoNamesList)<0)
		gizName = ""
		if (ItemsInList(GizmoNamesList)<1)
			return 1
		elseif (ItemsInList(GizmoNamesList)==1)
			gizName = StringFromList(0,GizmoNamesList)
		else
			Prompt gizName,"Gizmo to Use",popup,GizmoNamesList
			DoPrompt "Pick Gizmo to Use",gizName
			if (V_flag)
				return 1
			endif
		endif
	endif
	if (WhichListItem(gizName,GizmoNamesList)<0)
		return 1
	endif

	if (!(NframesAll>0))				// third choice is look for a wave
		String slist=reverseList(WaveListClass("MovieStepsList*","*","TEXT:1"))
		if (ItemsInList(slist)==1)
			Wave/T movieSteps=$StringFromList(0,slist)
		else
			String name
			Prompt name,"Movie Step List",popup,slist
			DoPrompt "Wave with Movie Steps",name
			if (V_flag)
				return 1
			endif
			Wave/T movieSteps=$name
			Nframes = NumberByKey("Nframes",note(movieSteps),"=")
			NframesAll = NframesFromMovieStepListWave(movieSteps)
		endif
		printIt = 1
	endif

	if (printIt)
		name = "$\"\""
		if (WaveType(movieSteps,2)==1)
			name = NameOfWave(movieSteps)
		endif
		printf "\rMakeMovieOfGizmo(\"%s\", %s",gizName,name
		if (testing)
			printf ", testing=%g",testing
		endif
		printf ")\t\t// "
		if (fromPanel)
			printf "  input from Panel,"
		endif
		printf "  Nframes=%g",Nframes
		printf "\r"
	endif
	Variable Nstep=DimSize(movieSteps,0)
	if (!(NframesAll>0) || Nstep<1)				// valid?
		return 1
	endif

	STRUCT gizmoMovieStruct gm
	gm.Nframes = Nframes
	gm.NframesAll = NframesAll
	gm.m = 0
	gm.gizmoName = gizName
	gm.imageName = UniqueName("GizmoImage",1,0)		// name of wave to hold image on graph
	gm.graphName = UniqueName("GizmoMovieGraph",6,0)	// name of graph
	gm.clipPlaneGroupName = ""

	Variable i, istep
	for (istep=0;i<Nstep;istep+=1)				// find first step in movieSteps that does something
		if (strlen(movieSteps[istep]))
			break
		endif
	endfor

	Make/N=(6,3)/FREE axes
	axes[0][0]= {1,-1,0,0,0,0}					// axes for {+X;-X;+Y;-Y;+Z;-Z}
	axes[0][1]= {0,0,1,-1,0,0}
	axes[0][2]= {0,0,0,0,1,-1}

	String list = StringByKey("Quaternion",movieSteps[istep])
	if (strlen(list))									// first step is quaternion
		istep += 1										// continue iwht next step
		gm.quaternion[0] = str2num(StringFromList(0,list,","))	
		gm.quaternion[1] = str2num(StringFromList(1,list,","))
		gm.quaternion[2] = str2num(StringFromList(2,list,","))
		gm.quaternion[3] = str2num(StringFromList(3,list,","))
	else													// first is not a quaternion, so get the current quaternion
		Wave q4=getGizmoQuaternion("")
		gm.quaternion[0] = q4[0]
		gm.quaternion[1] = q4[1]
		gm.quaternion[2] = q4[2]
		gm.quaternion[3] = q4[3]
	endif
print "starting gm.quaternion =",gm.quaternion[0],gm.quaternion[1],gm.quaternion[2],gm.quaternion[3]

	String progressWin = ProgressPanelStart("",showTime=1), cmd
	DoWindow/F $(gm.gizmoName)
#if (IgorVersion()<7)
	Execute "ModifyGizmo/N="+gm.gizmoName+" stopRotation"
	sprintf cmd, "ModifyGizmo/N=%s SETQUATERNION={%g,%g,%g,%g}",gm.gizmoName,gm.quaternion[0],gm.quaternion[1],gm.quaternion[2],gm.quaternion[3]
	Execute cmd
	Execute "ModifyGizmo/N="+gm.gizmoName+" update = 2"

	Execute "GetGizmo/N="+gm.gizmoName+" winPixels"		// get gizmo window position & size
	Variable left=NumVarOrDefault("V_left",NaN), right=NumVarOrDefault("V_right",NaN)
	Variable top=NumVarOrDefault("V_top",NaN), bottom=NumVarOrDefault("V_bottom",NaN)
	KillVariables/Z V_left, V_right, V_top, V_bottom
	Execute "ExportGizmo wave as \""+gm.imageName+"\""	// get a picture of the gizmo
#else
	ModifyGizmo/N=$(gm.gizmoName) stopRotation
	ModifyGizmo/N=$(gm.gizmoName) SETQUATERNION={gm.quaternion[0],gm.quaternion[1],gm.quaternion[2],gm.quaternion[3]}
	ModifyGizmo/N=$(gm.gizmoName) update = 2
	GetGizmo/N=$(gm.gizmoName) winPixels						// get gizmo window position & size
	Variable left=V_left, right=V_right, top=V_top, bottom=V_bottom
	ExportGizmo wave as $(gm.imageName)						// get a picture of the gizmo
#endif
	Wave GizmoImage=$(gm.imageName)

	// make the Graph to hold images of the Gizmo
	Display/W=(left,top,right,bottom )/N=$(gm.graphName)
	AppendImage GizmoImage
	ModifyGraph margin=1,width={Aspect,1}, tick=3, mirror=2, noLabel=2, standoff(left)=0
	SetAxis/A/R left

	Variable clipVal, c0,c1,dc,  angle, total, da, Nstatic
	String plane, axis, key, value, extName
	for (;istep<Nstep;istep+=1)
		key = StringFromList(0,movieSteps[istep],":")			// the key from istep
		value = StringFromList(1,movieSteps[istep],":")		// the value from istep
		if (strsearch(key,"Begin",0)==0)	//
			MovieIsOpen = NumVarOrDefault("MovieIsOpenGlobal",0)
			if (!MovieIsOpen)
				DoWindow/F $(gm.graphName)							// bring graph with new picture to front
				DoUpdate/W=$(gm.graphName)							// update graph to show picture
				NewMovie/P=home
				AddMovieFrame
				MovieIsOpen = 1
				Variable/G MovieIsOpenGlobal=1						// set global flag that movie is Open
				gm.m += 1
			endif

		elseif (strsearch(key,"Quaternion",0)==0)
			list = value
			sprintf cmd, "ModifyGizmo/N=%s SETQUATERNION={%s}",gm.gizmoName,list
			Execute cmd
			AddGizmoFrames2Movie(gm,1)								// add 1 frame to Movie
			ProgressPanelUpdate(progressWin,(gm.m)/(gm.NframesAll)*100)

		elseif (strsearch(key,"ClipReset",0)==0)				// remove clip plane
			RemoveGizmoMovieClipPlane(gm)
			AddGizmoFrames2Movie(gm,1)
			ProgressPanelUpdate(progressWin,(gm.m)/(gm.NframesAll)*100)

		elseif (strsearch(key,"Clip",0)==0)						// move the clip plane
			list = value													//	"Clip:+X, 0, 60, 5"
			plane = StringFromList(0,list,",")
			c0 = str2num(StringFromList(1,list,","))
			c1 = str2num(StringFromList(2,list,","))
			dc = str2num(StringFromList(3,list,","))
			dc = c1>c0 ? dc : -dc
			dc = testing ? 4*dc : dc
			if (strlen(gm.clipPlaneGroupName)<1)					// add clip plane group if it is not there
				gm.clipPlaneGroupName = AddGizmoClipPlaneGroup(GizmoMovieClipPlaneGroupName)
				InsertGizmoMovieClipPlanesGroup(gm.clipPlaneGroupName)
			endif

			if (dc<0)
				for (clipVal=c0; clipVal>c1; clipVal+=dc)
					ModifyGizmoClipPlaneGroup(gm.clipPlaneGroupName,plane,clipVal)
					AddGizmoFrames2Movie(gm,1)
					ProgressPanelUpdate(progressWin,(gm.m)/(gm.NframesAll)*100)
				endfor
			else
				for (clipVal=c0; clipVal<c1; clipVal+=dc)
					ModifyGizmoClipPlaneGroup(gm.clipPlaneGroupName,plane,clipVal)
					AddGizmoFrames2Movie(gm,1)
					ProgressPanelUpdate(progressWin,(gm.m)/(gm.NframesAll)*100)
				endfor
			endif

		elseif (strsearch(key,"Rotate",0)==0)						// rotate
			list = value													//	"Rotate:+Y, -60, 5"
			axis = StringFromList(0,list,",")
			total = str2num(StringFromList(1,list,","))
			da = abs(str2num(StringFromList(2,list,",")))
			da = total<0 ? -da : da
			da = testing ? 4*da : da

			i = WhichListItem(axis,directionTypes)				// directionTypes = "+X;-X;+Y;-Y;+Z;-Z"
			Make/N=3/D/FREE rotAxis = axes[i][p]
			for (angle=0; angle<=abs(total); angle+=abs(da))
				ProgressPanelUpdate(progressWin,(gm.m)/(gm.NframesAll)*100)
				AddRotationMovieFrame(gm,da,rotAxis)
			endfor

		elseif (strsearch(key,"Add_Static",0)==0)				// add some static frames
			Nstatic = NumberByKey("Add_Static",MovieSteps[istep])
			Nstatic = Nstatic==limit(Nstatic,0,500) ? round(Nstatic) : 1
			AddGizmoFrames2Movie(gm,Nstatic)						// add Nstatic frames

		else																	// an external function movie step
			extName = SelectString(strlen(key),"","AddGizmoMovieFrame_"+key)
			list = FunctionInfo(extName)
			if (StringMatch(StringByKey("TYPE",list),"UserDefined") && NumberByKey("N_PARAMS",list)==2 && NumberByKey("PARAM_0_TYPE",list)==4608 && NumberByKey("PARAM_1_TYPE",list)==8192)
				FUNCREF AddGizmoMovieFrame_Proto func=$extName
				func(gm,value)
			endif
		endif
	endfor

	DoWindow/K $(gm.graphName)												// kill the graph, we are done with it
	KillWaves/Z GizmoImage														// kill the image that was on the graph
	if (MovieIsOpen)
		CloseMovie
		MovieIsOpen = 0
	endif
	KillVariables/Z MovieIsOpenGlobal
	printf "total execution time = %s,   %d frames\r",Secs2Time(SecondsInProgressPanel(progressWin),5,0),gm.m
	DoWindow/K $progressWin
	RemoveGizmoMovieClipPlane(gm)										// clean up
	return 0
End


Function AddGizmoMovieFrame_Proto(gm,data)
	STRUCT gizmoMovieStruct &gm
	String data								// string with data
	return 0
End




Static Function AddGizmoFrames2Movie(gm,N)
	// take an image from the gizmo and add N frames
	STRUCT gizmoMovieStruct &gm
	Variable N			// number of frames to add, usually 1

#if (IgorVersion()<7)
	Execute "ModifyGizmo/N="+gm.gizmoName+" update = 1"
	Execute "ExportGizmo wave as \""+gm.imageName+"\""	// get a picture of the gizmo
#else
	ModifyGizmo/N=$(gm.gizmoName) update = 1
	ExportGizmo wave as $(gm.imageName)						// get a picture of the gizmo
#endif
	DoWindow/F $(gm.graphName)									// bring graph with new picture to front
	DoUpdate/W=$(gm.graphName)									// update graph to show picture
	Variable i, MovieIsOpen = NumVarOrDefault("MovieIsOpenGlobal",0)
	for (i=0;i<N;i+=1)					// add N frames
		if (MovieIsOpen)
			AddMovieFrame					// add the graph to the movie N times
		endif
		gm.m += 1							// update gm.m to reflect number of frames added to movie
	endfor
End

//  ============================ End of Make a Movie Section =============================  //
//  ======================================================================================  //



//  ======================================================================================  //
//  ============================ Start of Clip Plane Section =============================  //

Static Function InsertGizmoMovieClipPlanesGroup(group)
	String group
	// returns position in display list where group was placed

	if (strlen(group)<1)
		return NaN								// ERROR, no group name
	endif

#if (IgorVersion()<7)
	Execute "GetGizmo displayNameList"
	String DisplayList=StrVarOrDefault("S_DisplayNames","")
	KillStrings/Z S_DisplayNames
#else
	GetGizmo displayNameList
	String DisplayList=S_DisplayNames
#endif

	Variable m=WhichListItem(group,DisplayList)
	if (m>=0)
		return m									// group is already displayed
	endif

	Variable N=ItemsInList(DisplayList)
	for (m=N-1;m>=0;m-=1)
		if (StringMatch(StringFromList(m,DisplayList),"axes*"))
			m += 1
			break		
		elseif (StringMatch(StringFromList(m,DisplayList),"BeamLineAxesGroup*"))
			m += 1
			break
		endif
	endfor
	m = max(0,m)								// first is at least 0
	m = m>=N ? 0 : m							// if no axes found, put at top (i.e. 0)

#if (IgorVersion()<7)
	Execute "ModifyGizmo insertDisplayList="+num2istr(m)+", object="+group
#else
	ModifyGizmo insertDisplayList=m, object=$group
#endif
	return m
End


Function RemoveGizmoMovieClipPlane(gm)	// remove gizmo movie clip plane group from Gizmo
	STRUCT gizmoMovieStruct &gm
	if (strlen(gm.clipPlaneGroupName)<1)	// if no clip plane group, do nothing
		return 0
	endif

#if (IgorVersion()<7)
	String cmd
	sprintf cmd, "ModifyGizmo/N=%s userString={%s,\"\"}",gm.gizmoName,gm.clipPlaneGroupName
	Execute cmd
	sprintf cmd, "RemoveFromGizmo/N=%s/U=2/Z object=%s",gm.gizmoName,gm.clipPlaneGroupName
	Execute cmd
#else
	ModifyGizmo/N=$(gm.gizmoName) userString={gm.clipPlaneGroupName,""}
	RemoveFromGizmo/N=$(gm.gizmoName)/U=2/Z object=$(gm.clipPlaneGroupName)
#endif
	gm.clipPlaneGroupName = ""			// flags that clip plane group is gone
	return 0
End

//  ============================= End of Clip Plane Section ==============================  //
//  ======================================================================================  //


//  ======================================================================================  //
//  ============================= Start of Rotation Section ==============================  //

Static Function AddRotationMovieFrame(gm,da,axis)
	STRUCT gizmoMovieStruct &gm
	Variable da						// rotate by da about axis from current orientation (deg)
	Wave axis						// axis of rotation

	da *= PI/180 / 2				// convert to half-angle (rad)
	Make/N=4/D/FREE q0 = gm.quaternion[p]	// the base (e.g. current quaternion)

	Make/N=4/D/FREE qr=0		// quaternion representing the roattion
	qr[1,3] = axis[p-1]
	normalize(qr)
	qr *= sin(da)
	qr[0] = cos(da)
	Wave qnew = GizmoMovies#MultiplyQuaternions(qr,q0)

#if (IgorVersion()<7)
	String cmd, Nswitch=SelectString(strlen(gm.gizmoName),"","/N="+(gm.gizmoName))
	sprintf cmd, "ModifyGizmo%s SETQUATERNION={%s}",Nswitch,vec2str(qnew,bare=1)
	Execute cmd
#else
	ModifyGizmo/N=$(gm.gizmoName) SETQUATERNION={qnew[0],qnew[1],qnew[2],qnew[3]}
#endif
	gm.quaternion[0] = qnew[0]	// save as current quaternion
	gm.quaternion[1] = qnew[1]
	gm.quaternion[2] = qnew[2]
	gm.quaternion[3] = qnew[3]
	AddGizmoFrames2Movie(gm,1)
	return 0
End


Function/WAVE getGizmoQuaternion(gwin)
	String gwin

#if (IgorVersion()<7)
	Execute "GetGizmo gizmoNameList"
	String GizmoNames = StrVarOrDefault("S_GizmoNames","")
	KillVariables/Z S_GizmoNames
#else
	GetGizmo gizmoNameList
	String GizmoNames = S_GizmoNames
#endif
	if (strlen(gwin) && WhichListItem(gwin,GizmoNames)<0)
		return $""
	endif

#if (IgorVersion()<7)
	if (strlen(gwin) && WhichListItem(gwin,GizmoNames)<0)
		return $""
	endif
	String Nswitch=SelectString(strlen(gwin),"","/N="+(gwin))
	Execute "GetGizmo"+Nswitch+"/Z curMatrix"
#else
	if (strlen(gwin) && WhichListItem(gwin,GizmoNames)<0)
		return $""
	endif
	GetGizmo/N=$gwin/Z curMatrix
#endif
	Wave M_RotationMatrix=M_RotationMatrix
	if (!WaveExists(M_RotationMatrix))
		return $""
	endif

	Make/N=(3,3)/D/FREE rot = M_RotationMatrix[p][q]
	KillWaves/Z M_RotationMatrix
	Wave quat = RotMatrix2Quaternion(rot)
	return quat
End
//
Static Function/WAVE RotMatrix2Quaternion(rot)	// This is basically an Igor only routine, only with Gizmo
	Wave rot
	Variable q0,q1,q2,q3

	q0 = ( rot[0][0] + rot[1][1] + rot[2][2] + 1) / 4
	q1 = ( rot[0][0] - rot[1][1] - rot[2][2] + 1) / 4
	q2 = (-rot[0][0] + rot[1][1] - rot[2][2] + 1) / 4
	q3 = (-rot[0][0] - rot[1][1] + rot[2][2] + 1) / 4
	q0 = max(q0,0)		// if(q0 < 0) q0 = 0
	q1 = max(q1,0)		// if(q1 < 0) q1 = 0
	q2 = max(q2,0)		// if(q2 < 0) q2 = 0
	q3 = max(q3,0)		// if(q3 < 0) q3 = 0
	q0 = sqrt(q0)
	q1 = sqrt(q1)
	q2 = sqrt(q2)
	q3 = sqrt(q3)

	if(q0 >= q1 && q0 >= q2 && q0 >= q3)
		q0 *= 1
		q1 *= sign(rot[2][1] - rot[1][2])
		q2 *= sign(rot[0][2] - rot[2][0])
		q3 *= sign(rot[1][0] - rot[0][1])
	elseif(q1 >= q0 && q1 >= q2 && q1 >= q3)
		q0 *= sign(rot[2][1] - rot[1][2])
		q1 *= 1
		q2 *= sign(rot[1][0] + rot[0][1])
		q3 *= sign(rot[0][2] + rot[2][0])
	elseif(q2 >= q0 && q2 >= q1 && q2 >= q3)
		q0 *= sign(rot[0][2] - rot[2][0])
		q1 *= sign(rot[1][0] + rot[0][1])
		q2 *= 1
		q3 *= sign(rot[2][1] + rot[1][2])
	elseif(q3 >= q0 && q3 >= q1 && q3 >= q2)
		q0 *= sign(rot[1][0] - rot[0][1])
		q1 *= sign(rot[2][0] + rot[0][2])
		q2 *= sign(rot[2][1] + rot[1][2])
		q3 *= 1
	endif

	Make/N=4/D/FREE quaternion={q1,q2,q3,q0}		// I do not understand the order
	normalize(quaternion)
	return quaternion
End


Static Function/WAVE MultiplyQuaternions(u,v)	// calculate the product: u v
	Wave u,v					// input quaternions (4-vectors)
	//	See:	http://en.wikipedia.org/wiki/Quaternion
	//	and	http://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation

	Make/N=4/FREE/D uv
	uv[0] = u[0]*v[0] - u[1]*v[1] - u[2]*v[2] - u[3]*v[3]
	uv[1] = u[0]*v[1] + u[1]*v[0] + u[2]*v[3] - u[3]*v[2]
	uv[2] = u[0]*v[2] - u[1]*v[3] + u[2]*v[0] + u[3]*v[1]
	uv[3] = u[0]*v[3] + u[1]*v[2] - u[2]*v[1] + u[3]*v[0]
	return uv
End

//  ============================== End of Rotation Section ===============================  //
//  ======================================================================================  //






//  ======================================================================================  //
//  ========================= Start of Movie Steps Panel Section =========================  //

Function MakeMovieStepsPanel([N])
	Variable N
	N = ParamIsDefault(N) ? NaN : N
	N = N==limit(round(N),2,MAX_MOVIE_STEPS) ? N : NaN
	N = numtype(N) ? 15 : N

	if (ItemsInList(WinList("MovieStepsPanel",";","WIN:64")))
		DoWindow/F MovieStepsPanel
		return 0
	endif

	String stepList = reverseList(WaveListClass("MovieStepsList*","*","TEXT:1"))
	Wave/T MovieStepsList= $StringFromList(0,stepList)
	if (!WaveExists(MovieStepsList))
		Wave/T MovieStepsList= DefaultMovieStepListWave()
	endif
	print "    Starting MovieStepsPanel with data in ",NameOfWave(MovieStepsList)
	N = min(N,DimSize(MovieStepsList,0))
	GizmoMovies#ResetMovieStepTypes()			// re-set list of step types for popup menus

	Variable left=200, top=200, startControls=40
	NewPanel/W=(left,top,left+490,top+1+startControls+25*N)/N=MovieStepsPanel/K=1

	Variable y0, i
	SetDrawLayer UserBack
	SetDrawEnv linethick= 2
	DrawLine 0,32,500,32
	for (i=0;i<(N-1);i+=1)
		SetDrawEnv dash= 1
		y0 = 22 +  i*25 + startControls
		DrawLine 0,y0,500,y0
	endfor

	Button LoadWaveButton,pos={10,5},size={110,20},proc=GizmoMovies#MovieStepPanelButtonProc,title="Load from Wave"
	Button LoadWaveButton,fColor=(16385,28398,65535), disable=2
	Button SaveWaveButton,pos={140,5},size={100,20},proc=GizmoMovies#MovieStepPanelButtonProc,title="Save to Wave"
	Button SaveWaveButton,fColor=(16385,28398,65535), disable=2
	Button MakeMovieWaveButton,pos={260,5},size={100,20},proc=GizmoMovies#MovieStepPanelButtonProc,title="Make Movie"
	Button MakeMovieWaveButton,fColor=(32769,65535,32768), disable=2
	ValDisplay NframesValdisp,pos={390,7},size={90,17},title="\\f02frames",fSize=12, format="%d", frame=0
	ValDisplay NframesValdisp,limits={0,0,0},barmisc={0,1000}, value= _NUM:0

	String name
	String MovieTypesStr = "\""+StrVarOrDefault("root:Packages:GizmoMovies:MovieStepTypes",MovieStepTypesBase)+"\""

	String DirectionTypesStr = "\""+directionTypes+"\""
	for (i=0;i<N;i+=1)
		top = i*25 + startControls
		sprintf name, "MoviePop%d", i
		PopupMenu $name,pos={5,top},size={33,20},proc=GizmoMovies#MovieStepPopupProc
		PopupMenu $name,mode=1,popvalue=" ",value= #MovieTypesStr

		// add all input fields, and set them all to be hidden (disable=1)
		sprintf name,"MoviePop%d_plane",i
		PopupMenu $name,pos={110,top},size={46,20},disable=1
		PopupMenu $name,mode=1,popvalue="+X",value= #DirectionTypesStr
		sprintf name,"MoviePop%d_StartClip",i
		SetVariable $name,pos={175,top},size={89,19},title="start",disable=1,proc=GizmoMovies#MovieStepPanelSetVarProc
		SetVariable $name,fSize=12,limits={-Inf,Inf,0},value= _NUM:0
		sprintf name,"MoviePop%d_EndClip",i
		SetVariable $name,pos={282,top},size={89,19},title="end",fSize=12,disable=1,proc=GizmoMovies#MovieStepPanelSetVarProc
		SetVariable $name,limits={-Inf,Inf,0},value= _NUM:10
		sprintf name,"MoviePop%d_DeltaClip",i
		SetVariable $name,pos={388,top},size={61,19},title="Æ",fSize=12,disable=1,proc=GizmoMovies#MovieStepPanelSetVarProc
		SetVariable $name,limits={-Inf,Inf,0},value= _NUM:1

		sprintf name,"MoviePop%d_axis",i
		PopupMenu $name,pos={110,top},size={46,20},disable=1
		PopupMenu $name,mode=1,popvalue="+X",value= #DirectionTypesStr
		sprintf name,"MoviePop%d_TotalAngle",i
		SetVariable $name,pos={170,top},size={89,19},title="total¡",disable=1,proc=GizmoMovies#MovieStepPanelSetVarProc
		SetVariable $name,fSize=12,limits={-720,720,0},value= _NUM:60
		sprintf name,"MoviePop%d_DeltaAngle",i
		SetVariable $name,pos={275,top},size={61,19},title="Æ¡",fSize=12,disable=1,proc=GizmoMovies#MovieStepPanelSetVarProc
		SetVariable $name,limits={-720,720,0},value= _NUM:5

		sprintf name,"MoviePop%d_Qw",i
		SetVariable $name,pos={110,top},size={89,19},title="W",fSize=12,disable=1,proc=GizmoMovies#MovieStepPanelSetVarProc
		SetVariable $name,limits={-1,1,0},value= _NUM:0
		sprintf name,"MoviePop%d_Qx",i
		SetVariable $name,pos={205,top},size={89,19},title="X",fSize=12,disable=1,proc=GizmoMovies#MovieStepPanelSetVarProc
		SetVariable $name,limits={-1,1,0},value= _NUM:0
		sprintf name,"MoviePop%d_Qy",i
		SetVariable $name,pos={300,top},size={89,19},title="Y",fSize=12,disable=1,proc=GizmoMovies#MovieStepPanelSetVarProc
		SetVariable $name,limits={-1,1,0},value= _NUM:0
		sprintf name,"MoviePop%d_Qz",i
		SetVariable $name,pos={395,top},size={89,19},title="Z",fSize=12,disable=1,proc=GizmoMovies#MovieStepPanelSetVarProc
		SetVariable $name,limits={-1,1,0},value= _NUM:0

		sprintf name,"MoviePop%d_Nstatic",i
		SetVariable $name,pos={110,top},size={100,19},proc=GizmoMovies#MovieStepPanelSetVarProc
		SetVariable $name,title="# of frames",disable=1,fSize=12,limits={0,inf,0},value= _NUM:1

		sprintf name,"MoviePop%d_externalData",i
		SetVariable $name,pos={110,top},size={370,19},value= _STR:"",proc=GizmoMovies#MovieStepPanelSetVarProc
		SetVariable $name,disable=1,fSize=12,value= _STR:""
	endfor
	FillMovieStepsPanelFromWave(MovieStepsList)
	SetWindow kwTopWin userdata(MovieSteps)="N:"+num2istr(N)
	UpdateMovieStepPanelButtons()
End
//
Static Function MovieStepPopupProc(pa) : PopupMenuControl
	STRUCT WMPopupAction &pa
	if (pa.eventCode == 2)		// mouse up
		Variable i
		sscanf pa.ctrlName, "MoviePop%d", i
		EnableHideInputFields(i,pa.popStr)
		UpdateMovieStepPanelButtons()
	endif
	return 0
End
//
Static Function EnableHideInputFields(i,type)
	Variable i						// step index, [0,N-1]
	String type						// {Quaternion,Rotate,Clip,ClipReset}

	String name
	Variable disable, foundOne=0

	disable = !StringMatch(type,"Quaternion")
	foundOne += !disable
	sprintf name, "MoviePop%d_Qw", i
	SetVariable $name,disable=disable
	sprintf name, "MoviePop%d_Qx", i
	SetVariable $name,disable=disable
	sprintf name, "MoviePop%d_Qy", i
	SetVariable $name,disable=disable
	sprintf name, "MoviePop%d_Qz", i
	SetVariable $name,disable=disable

	disable = !StringMatch(type,"Rotate")
	foundOne += !disable
	sprintf name, "MoviePop%d_axis", i
	PopupMenu $name,disable=disable
	sprintf name, "MoviePop%d_TotalAngle", i
	SetVariable $name,disable=disable
	sprintf name, "MoviePop%d_DeltaAngle", i
	SetVariable $name,disable=disable

	disable = !StringMatch(type,"Clip")
	foundOne += !disable
	sprintf name, "MoviePop%d_plane", i
	PopupMenu $name,disable=disable
	sprintf name, "MoviePop%d_StartClip", i
	SetVariable $name,disable=disable
	sprintf name, "MoviePop%d_EndClip", i
	SetVariable $name,disable=disable
	sprintf name, "MoviePop%d_DeltaClip", i
	SetVariable $name,disable=disable

	disable = !StringMatch(type,"Add_Static")
	foundOne += !disable
	sprintf name, "MoviePop%d_Nstatic", i
	SetVariable $name,disable=disable

	sprintf name, "MoviePop%d_externalData", i
	disable = WhichListItem(type,MovieStepTypesBase)>=0
	SetVariable $name,disable=disable
End
//
Static Function MovieStepPanelButtonProc(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	if (ba.eventCode == 2)		// mouse up
		if (StringMatch(ba.ctrlName,"LoadWaveButton"))
			FillMovieStepsPanelFromWave($"",printIt=1)
		elseif (StringMatch(ba.ctrlName,"SaveWaveButton"))
			SaveMovieSteps2Wave(printIt=1)
		elseif (StringMatch(ba.ctrlName,"MakeMovieWaveButton"))
			MakeMovieOfGizmo("",$"",printIt=1)
		endif
	endif
	return 0
End
//
Static Function MovieStepPanelSetVarProc(sva) : SetVariableControl
	STRUCT WMSetVariableAction &sva
	switch( sva.eventCode )
		case 1:		// mouse up
		case 2:		// Enter key
			UpdateMovieStepPanelButtons()
	endswitch
	return 0
End
//
Static Function UpdateMovieStepPanelButtons()
	// update button enable/disable for all the buttons on the MovieStepsPanel

	Wave/T steps=GatherMovieStepFromPanel()
	if (!WaveExists(steps))
		return 1
	endif
	Variable Nframes = NframesFromMovieStepListWave(steps)
	ValDisplay NframesValdisp,value= _NUM:Nframes

	Variable disable = ItemsInList(WaveListClass("MovieStepsList*","*","TEXT:1"))==0 ? 2 : 0
	Button LoadWaveButton, disable=disable

	disable = NframesFromMovieStepListWave(steps,all=1)>0 ? 0 : 2
	Button SaveWaveButton, disable=0

	disable = strlen(WinList("*","","WIN:"+num2istr(GIZMO_WIN_BIT)))>0 ? 0 : 2
	Button MakeMovieWaveButton, disable=disable

	return 0
End



Function FillMovieStepsPanelFromWave(ms,[printIt])
	Wave/T ms
	Variable printIt
	printIt = ParamIsDefault(printIt) ? NaN : printIt
	printIt = numtype(printIt) ? strlen(GetRTStackInfo(2))==0 : !(!printIt)

	if (!WaveExists(ms))
		String slist=reverseList(WaveListClass("MovieStepsList*","*","TEXT:1"))
		if (ItemsInList(slist)<1)
			return 1
		elseif (ItemsInList(slist)==1)
			Wave/T ms=$StringFromList(0,slist)
		else
			String name
			Prompt name,"Movie Step List",popup,slist
			DoPrompt "Wave with Movie Steps",name
			if (V_flag)
				return 1
			endif
			Wave/T ms=$name
		endif
		printIt = 1
	endif
	if (!WaveExists(ms))
		return 1
	endif

	if (printIt)
		printf "FillMovieStepsPanelFromWave(%s",NameOfWave(ms)
		//		if (!ParamIsDefault(printIt))
		//			printf ", printIt=1"
		//		endif
		printf ")\r"
	endif

	String win="MovieStepsPanel"
	String pnote=GetUserData(win,"","MovieSteps")
	Variable Npanel=round(NumberByKey("N",pnote)), Nwave=DimSize(ms,0)
	Variable N = numtype(Npanel) ? Nwave : min(Npanel,Nwave), i, mode
	N = round(min(N,MAX_MOVIE_STEPS))
	if (!(N>=1))
		return 1
	endif

	String type, cname, list
	for (i=0;i<N;i+=1)
		type = StringFromList(0,ms[i],":")
		type = SelectString(strlen(type)," ",type)
		sprintf name, "MoviePop%d", i
		mode = WhichListItem(type,StrVarOrDefault("root:Packages:GizmoMovies:MovieStepTypes",MovieStepTypesBase))+1
		PopupMenu $name,mode=mode

		list = StringByKey(type,ms[i])
		if (StringMatch(type,"Quaternion"))
			cname = name+"_Qw"
			SetVariable $cname,value= _NUM:str2num(StringFromList(0,list,","))
			cname = name+"_Qx"
			SetVariable $cname,value= _NUM:str2num(StringFromList(1,list,","))
			cname = name+"_Qy"
			SetVariable $cname,value= _NUM:str2num(StringFromList(2,list,","))
			cname = name+"_Qz"
			SetVariable $cname,value= _NUM:str2num(StringFromList(3,list,","))

		elseif (StringMatch(type,"Clip"))
			cname = name+"_plane"
			mode = WhichListItem(StringFromList(0,list,","),directionTypes)+1
			PopupMenu $cname,mode=mode
			cname = name+"_StartClip"
			SetVariable $cname,value= _NUM:str2num(StringFromList(1,list,","))
			cname = name+"_EndClip"
			SetVariable $cname,value= _NUM:str2num(StringFromList(2,list,","))
			cname = name+"_DeltaClip"
			SetVariable $cname,value= _NUM:abs(str2num(StringFromList(3,list,",")))

		elseif (StringMatch(type,"Rotate"))
			cname = name+"_axis"
			mode = WhichListItem(StringFromList(0,list,","),directionTypes)+1
			PopupMenu $cname,mode=mode
			cname = name+"_TotalAngle"
			SetVariable $cname,value= _NUM:str2num(StringFromList(1,list,","))
			cname = name+"_DeltaAngle"
			SetVariable $cname,value= _NUM:abs(str2num(StringFromList(2,list,",")))

		elseif (StringMatch(type,"Add_Static"))
			cname = name+"_Nstatic"
			SetVariable $cname,value= _NUM:str2num(StringFromList(0,list,","))

		elseif (WhichListItem(type,MovieStepTypesBase)<0)
			cname = name+"_externalData"					// for external, just load whole string
			SetVariable $cname,value= _STR:list
		endif
		EnableHideInputFields(i,type)
	endfor
	UpdateMovieStepPanelButtons()
	return 0
End



Function/WAVE SaveMovieSteps2Wave([printIt])
	Variable printIt
	printIt = ParamIsDefault(printIt) ? NaN : printIt
	printIt = numtype(printIt) ? strlen(GetRTStackInfo(2))==0 : !(!printIt)

	Wave/T stepsPanel=GatherMovieStepFromPanel()
	if (!WaveExists(stepsPanel))
		return $""
	endif

	Variable N=DimSize(stepsPanel,0), Nframes=NframesFromMovieStepListWave(stepsPanel)
	Variable NframesAll = NframesFromMovieStepListWave(stepsPanel,all=1)
	if (NframesAll<1)
		if (printIt)
			printf "SaveMovieSteps2Wave()		// FAILED,  %g frames in movie\r",NframesAll
		endif
		return $""
	endif

	String name=UniqueName("MovieStepsList",1,0)
	Make/N=(N)/T/O $name/WAVE=stepsWave
	stepsWave = stepsPanel

	String wnote="waveClass=MovieStepsList", list
	wnote = ReplaceNumberByKey("Nframes",wnote,Nframes,"=")
	Note/K stepsWave, wnote
	UpdateMovieStepPanelButtons()

	if (printIt)
		printf "SaveMovieSteps2Wave()		// saved to wave \"%s\",  %g frames in movie\r",NameOfWave(stepsWave),Nframes
	endif
	return stepsWave
End
//
Static Function/WAVE GatherMovieStepFromPanel()
	if (strlen(WinList("MovieStepsPanel","","WIN:64"))==0)
		return $""
	endif

	String pnote=GetUserData("MovieStepsPanel","","MovieSteps")
	Variable N=round(NumberByKey("N",pnote))
	if (!(N>=1 && N<MAX_MOVIE_STEPS))
		return $""
	endif

	Make/N=(N)/T/FREE steps=""
	String name, cname
	Variable i
	for (i=0;i<N;i+=1)
		sprintf name, "MoviePop%d", i
		ControlInfo/W=MovieStepsPanel $name
		if (stringMatch(S_Value,"Begin"))			// Begin
			steps[i] = "Begin"
		elseif (stringMatch(S_Value,"Quaternion"))	// Quaternion
			steps[i] = "Quaternion:"
			cname = name+"_Qw"
			ControlInfo/W=MovieStepsPanel $cname
			steps[i] += num2str(V_Value)
			cname = name+"_Qx"
			ControlInfo/W=MovieStepsPanel $cname
			steps[i] += ", "+num2str(V_Value)
			cname = name+"_Qy"
			ControlInfo/W=MovieStepsPanel $cname
			steps[i] += ", "+num2str(V_Value)
			cname = name+"_Qz"
			ControlInfo/W=MovieStepsPanel $cname
			steps[i] += ", "+num2str(V_Value)
		elseif (stringMatch(S_Value,"Rotate"))		// Rotate
			steps[i] = "Rotate:"
			cname = name+"_axis"
			ControlInfo/W=MovieStepsPanel $cname
			steps[i] += S_Value
			cname = name+"_TotalAngle"
			ControlInfo/W=MovieStepsPanel $cname
			steps[i] += ", "+num2str(V_Value)
			cname = name+"_DeltaAngle"
			ControlInfo/W=MovieStepsPanel $cname
			steps[i] += ", "+num2str(abs(V_Value))
		elseif (stringMatch(S_Value,"Clip"))			// Clip
			steps[i] = "Clip:"
			cname = name+"_plane"
			ControlInfo/W=MovieStepsPanel $cname
			steps[i] += S_Value
			cname = name+"_StartClip"
			ControlInfo/W=MovieStepsPanel $cname
			steps[i] += ", "+num2str(V_Value)
			cname = name+"_EndClip"
			ControlInfo/W=MovieStepsPanel $cname
			steps[i] += ", "+num2str(V_Value)
			cname = name+"_DeltaClip"
			ControlInfo/W=MovieStepsPanel $cname
			steps[i] += ", "+num2str(abs(V_Value))
		elseif (stringMatch(S_Value,"ClipReset"))	// ClipReset
			steps[i] = "ClipReset"
		elseif (stringMatch(S_Value,"Add_Static"))	// add N static frames
			steps[i] = "Add_Static:"
			cname = name+"_Nstatic"
			ControlInfo/W=MovieStepsPanel $cname
			steps[i] += num2str(V_Value)
		elseif (exists("AddGizmoMovieFrame_"+S_Value)==6)	// an external function
			steps[i] = S_Value+":"
			cname = name+"_externalData"
			ControlInfo/W=MovieStepsPanel $cname
			steps[i] += S_Value
		else
			steps[i] = ""										// empty, so nothing
		endif
	endfor
	return steps
End



Static Function NframesFromMovieStepListWave(steps,[all])
	Wave/T steps
	Variable all
	all = ParamIsDefault(all) ? 0 : all
	all = numtype(all) ? 0 : !(!all)
	if (!WaveExists(steps))
		return NaN
	endif

	String list
	Variable i, N=DimSize(steps,0), Nframes
	Variable started=all
	for (i=0,Nframes=0; i<N; i+=1)
		if (strsearch(steps[i],"Begin",0)==0)
			started = 1					// nothing counts until the first "Begin"
			Nframes += 1
		elseif (!started)				// if not started, the following do not count
			continue
		elseif (strsearch(steps[i],"Quaternion",0)==0)
			Nframes += 1
		elseif (strsearch(steps[i],"ClipReset",0)==0)
			Nframes += 1
		elseif (strsearch(steps[i],"Clip",0)==0)
			Variable c0,c1,dc
			list = StringByKey("Clip",steps[i])
			c0 = str2num(StringFromList(1,list,","))
			c1 = str2num(StringFromList(2,list,","))
			dc = abs(str2num(StringFromList(3,list,",")))
			dc = c1>c0 ? dc : -dc
			Nframes += max(0,floor((c1-c0)/dc + 1))
		elseif (strsearch(steps[i],"Rotate",0)==0)
			Variable total,da
			list = StringByKey("Rotate",steps[i])
			total = str2num(StringFromList(1,list,","))
			da = abs(str2num(StringFromList(2,list,",")))
			da = total > 0 ? da : -da
			Nframes += max(0,floor(total/da + 1))
		elseif (strsearch(steps[i],"Add_Static",0)==0)
			Variable Nstatic=NumberByKey("Add_Static",steps[i])
			Nstatic = Nstatic==limit(Nstatic,0,500) ? round(Nstatic) : 1
			Nframes += Nstatic
		else								// try an external command
			String data = StringFromList(1,steps[i],":")
			Variable Nextern = NumberByKey("Nframes",data,"=",",")
			Nextern = Nextern>0 && Nextern<500 ? round(Nextern) : 0
			Nframes += Nextern
		endif
	endfor
	return Nframes
End


Static Function/WAVE DefaultMovieStepListWave()
	// a defualt movie step list, used as a starting point
	Make/N=15/T/FREE ms=""
	ms[0] = "Quaternion:-0.70888, -0.33848, 0.29885, 0.54185"
	ms[1] = "Begin"
	ms[2] = "Add_Static:1"
	ms[3] = "Rotate:+Y, 60, 5"
	ms[4] = "Clip:-Z, 0, 1, 0.2"
	ms[5] = "Clip:-Z, 1, 0, 0.2"
	ms[6] = "ClipReset"
	ms[7] = "Rotate:+Y, -60, 5"
	ms[8] = "Add_Static:2"
	Note/K ms, "waveClass=MovieStepsList;Nframes=39;"
	return ms
End

//  ========================== End of Movie Steps Panel Section ==========================  //
//  ======================================================================================  //





//  ======================================================================================  //
//  ========================== Start of TESTING Movie Of Gizmo ===========================  //

Static Function MakeTestGizmoData()			// some test data for a simple Gizmo
	Make/N=(9,3)/O xyzTestGizmoMovie
	xyzTestGizmoMovie[0][0]= {0.2,1.2,0.2,1.2,0.2,1.2,0.2,1.2,0.7}
	xyzTestGizmoMovie[0][1]= {0.2,0.2,1.2,1.2,0.2,0.2,1.2,1.2,0.7}
	xyzTestGizmoMovie[0][2]= {0.2,0.2,0.2,0.2,1.2,1.2,1.2,1.2,0.7}

	Make/N=(9,4)/O rgbaTestGizmoMovie			// rgba that goes with the xyz
	rgbaTestGizmoMovie[0][0]= {0,1,0,0,0,1,1,0.5,0.5}
	rgbaTestGizmoMovie[0][1]= {0,0,1,0,1,0,1,0.5,0}
	rgbaTestGizmoMovie[0][2]= {0,0,0,1,1,1,0,0.5,0.5}
	rgbaTestGizmoMovie[0][3]= {1,1,1,1,1,1,1,1,0}
End

Static Function MakeGizmoMovieTestGizmo()	// show the test data in a Gizmo
	if(exists("NewGizmo")!=4)						// Do nothing if the Gizmo XOP is not available.
		DoAlert 0, "Gizmo XOP must be installed"
		return 1
	endif

#if (IgorVersion()<7)
	Execute "NewGizmo/K=1/N=GizmoMovieTest0/T=\"GizmoMovieTest0\" /W=(151,297,746,892)"
	Execute "ModifyGizmo startRecMacro"
	Execute "AppendToGizmo Scatter=root:xyzTestGizmoMovie,name=scatter0"
	Execute "ModifyGizmo ModifyObject=scatter0 property={ scatterColorType,1}"
	Execute "ModifyGizmo ModifyObject=scatter0 property={ Shape,2}"
	Execute "ModifyGizmo ModifyObject=scatter0 property={ size,3}"
	Execute "ModifyGizmo ModifyObject=scatter0 property={ colorWave,root:rgbaTestGizmoMovie}"
	Execute "AppendToGizmo Axes=boxAxes,name=axes0"
	Execute "ModifyGizmo ModifyObject=axes0,property={-1,axisScalingMode,1}"
	Execute "ModifyGizmo ModifyObject=axes0,property={-1,axisColor,0,0,0,1}"
	Execute "ModifyGizmo ModifyObject=axes0,property={0,ticks,3}"
	Execute "ModifyGizmo ModifyObject=axes0,property={1,ticks,3}"
	Execute "ModifyGizmo ModifyObject=axes0,property={2,ticks,3}"
	Execute "ModifyGizmo ModifyObject=axes0,property={0,axisLabel,1}"
	Execute "ModifyGizmo ModifyObject=axes0,property={1,axisLabel,1}"
	Execute "ModifyGizmo ModifyObject=axes0,property={2,axisLabel,1}"
	Execute "ModifyGizmo ModifyObject=axes0,property={0,axisLabelText,\"X\"}"
	Execute "ModifyGizmo ModifyObject=axes0,property={1,axisLabelText,\"Y\"}"
	Execute "ModifyGizmo ModifyObject=axes0,property={2,axisLabelText,\"Z\"}"
	Execute "ModifyGizmo modifyObject=axes0 property={Clipped,0}"
	Execute "ModifyGizmo setDisplayList=0, object=axes0"
	Execute "ModifyGizmo setDisplayList=1, object=scatter0"
	Execute "ModifyGizmo SETQUATERNION={-0.708885,-0.338485,0.298854,0.541845}"
	Execute "ModifyGizmo autoscaling=1"
	Execute "ModifyGizmo currentGroupObject=\"\""
	Execute "ModifyGizmo compile"
	Execute "ModifyGizmo showAxisCue=1"
	Execute "ModifyGizmo endRecMacro"
#else
	NewGizmo/K=1/N=GizmoMovieTest0/T="GizmoMovieTest0" /W=(151,297,746,892)
	ModifyGizmo startRecMacro
	AppendToGizmo Scatter=root:xyzTestGizmoMovie,name=scatter0
	ModifyGizmo ModifyObject=scatter0 property={ scatterColorType,1}
	ModifyGizmo ModifyObject=scatter0 property={ Shape,2}
	ModifyGizmo ModifyObject=scatter0 property={ size,3}
	ModifyGizmo ModifyObject=scatter0 property={ colorWave,root:rgbaTestGizmoMovie}
	AppendToGizmo Axes=boxAxes,name=axes0
	ModifyGizmo ModifyObject=axes0,property={-1,axisScalingMode,1}
	ModifyGizmo ModifyObject=axes0,property={-1,axisColor,0,0,0,1}
	ModifyGizmo ModifyObject=axes0,property={0,ticks,3}
	ModifyGizmo ModifyObject=axes0,property={1,ticks,3}
	ModifyGizmo ModifyObject=axes0,property={2,ticks,3}
	ModifyGizmo ModifyObject=axes0,property={0,axisLabel,1}
	ModifyGizmo ModifyObject=axes0,property={1,axisLabel,1}
	ModifyGizmo ModifyObject=axes0,property={2,axisLabel,1}
	ModifyGizmo ModifyObject=axes0,property={0,axisLabelText,"X"}
	ModifyGizmo ModifyObject=axes0,property={1,axisLabelText,"Y"}
	ModifyGizmo ModifyObject=axes0,property={2,axisLabelText,"Z"}
	ModifyGizmo modifyObject=axes0 property={Clipped,0}
	ModifyGizmo setDisplayList=0, object=axes0
	ModifyGizmo setDisplayList=1, object=scatter0
	ModifyGizmo SETQUATERNION={-0.708885,-0.338485,0.298854,0.541845}
	ModifyGizmo autoscaling=1
	ModifyGizmo currentGroupObject=""
	ModifyGizmo compile
	ModifyGizmo showAxisCue=1
	ModifyGizmo endRecMacro
#endif
End

//  =========================== End of TESTING Movie Of Gizmo ============================  //
//  ======================================================================================  //


Static Function/T ResetMovieStepTypes()
	// reset the Global string root:Packages:GizmoMovies:MovieStepTypes
	// check for all functions of the type "AddGizmoMovieFrame_*" containing two arguments (struct,string)
	String/G root:Packages:GizmoMovies:MovieStepTypes=MovieStepTypesBase
	SVAR MovieStepTypes=root:Packages:GizmoMovies:MovieStepTypes

	String flist=FunctionList("AddGizmoMovieFrame_*",";","NPARAMS:2,KIND:2")
	flist = RemoveFromList("AddGizmoMovieFrame_Proto",flist)		// remove the proto from the list
	String name, ilist
	Variable i,N=ItemsInList(flist)
	for (i=0;i<N;i+=1)
		name = StringFromList(i,flist)
			ilist = FunctionInfo(name)
			if (StringMatch(StringByKey("TYPE",ilist),"UserDefined") && NumberByKey("N_PARAMS",ilist)==2 && NumberByKey("PARAM_0_TYPE",ilist)==4608 && NumberByKey("PARAM_1_TYPE",ilist)==8192)
				MovieStepTypes += ";"+ReplaceString("AddGizmoMovieFrame_",name,"")
			endif
	endfor
	return MovieStepTypes
End

Static Function InitGizmoMovies()
	GizmoUtil#InitGizmoUtilityGeneral()
#if (IgorVersion()<7)
	Execute/Q/Z "GizmoMenu AppendItem={JZTmov0,\"Movie Steps Panel...\", \"MakeMovieStepsPanel()\"}"
#endif
	NewDataFolder/O root:Packages
	NewDataFolder/O root:Packages:GizmoMovies
	ResetMovieStepTypes()
End
