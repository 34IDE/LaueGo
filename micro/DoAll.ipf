#pragma rtGlobals=1		// Use modern global access method.
#pragma ModuleName=DoAll
#pragma version = 0.02

//Static Constant pxPrimary=NaN, pyPrimary=NaN, ycPrimary = NaN		// needed by DoAll(),  actually DoAllSingleSpotInfo(), and also by FillMovieMore()
//Static Constant hPrimary=NaN, kPrimary=NaN, lPrimary=NaN			// following 3 lines only needed when using two spos
//Static Constant pxSecondary=NaN, pySecondary=NaN
//Static Constant hSecondary=NaN, kSecondary=NaN, lSecondary=NaN


Menu "Macros"
	"Display 3D Rotation Gizmo",MakeGizmoRot($"")
	"Calculate Rotation3D from R3dXYZ",Rotation3DfromR3dXYZ("",1)
	"Trim junk from top of gizmo",TrimTopRgizmo()
	"-"
	"Make Graph RX,RH,RF",Make_Graph_RX_RH_RF()
	"-"
	"Do All Spots [for single spots]",DoAll(-1)
	"Change X0&H0 in ALL wave notes", ChangeAllX0H0(NaN,NaN)
	"Show DoAll procedure",/Q,DisplayProcedure "DoAll"
End



Function TrimTopRgizmo()						// find the top surface, and trim all points above it
	Wave Rgizmo=Rgizmo, total=totalintensity3d
	Variable topValue = Zintensity(total)
	// print topValue

	Variable i,j,Nx=DimSize(Rgizmo,0), Ny=DimSize(Rgizmo,1), Nz=DimSize(Rgizmo,2)
	Make/N=(Nz)/O JZT_columnIntensity
	Wave wcol = JZT_columnIntensity
	for (j=0;j<Ny;j+=1)
		for (i=0;i<Nx;i+=1)
			wcol = total[i][j][p]
			FindLevel/EDGE=1/P/Q wcol, topValue
			if (!V_flag && V_LevelX>1)
				Rgizmo[i][j][0,floor(V_LevelX-1)] = NaN
			endif
		endfor
	endfor
	KillWaves/Z JZT_columnIntensity
	return 0
End
//
Function Zintensity(mat)			// finds top of sample volume (40% above starting point to max)
	Wave mat
	if (!WaveExists(mat))
		Wave mat=totalintensity3d
	endif
	if (!WaveExists(mat))
		DoAlert 0,"totalintensity3d does not exist"
		return NaN
	endif

	Variable i, Nz=DimSize(mat,2)
	Make/N=(Nz)/O JZT_Zintensity
	for (i=0;i<Nz;i+=1)
		ImageTransform/P=(i)/Q getPlane mat
		Wave M_ImagePlane=M_ImagePlane
		WaveStats/Q/M=1 M_ImagePlane
		JZT_Zintensity[i] = V_avg
	endfor
	SetScale/P x,DimOffset(mat,2),DimDelta(mat,2),WaveUnits(mat,2),JZT_Zintensity
	WaveStats/Q/M=1 JZT_Zintensity
	JZT_Zintensity = numtype(JZT_Zintensity) ? V_min : JZT_Zintensity
	Variable level = (V_max-JZT_Zintensity[0])*0.4 + JZT_Zintensity[0]	// 40% above start
	FindLevel/EDGE=1/P/Q JZT_Zintensity, level
	level = V_flag ? NaN : JZT_Zintensity[floor(V_LevelX)]
	KIllWaves/Z JZT_Zintensity, M_ImagePlane
	return level
End



Function trimVolume(lim)
	Variable lim
	if (!(lim>1))
		DoAlert 0,"lim = "+num2str(lim)+", too small"
		return 1
	endif

	Wave XX=XX, HH=HH,FF=FF, RX=RX, RH=RH, RF=RF
	Wave totalAngles=totalAngles, totalIntensity=totalIntensity, totalPeakIntensity=totalPeakIntensity
	Wave/T imageNames=imageNames

	Variable i, N=numpnts(XX)
	Variable timer=startMSTimer
	printf "at start N=%d		",N
	for (i=N-1;i>=0;i-=1)
		if (abs(XX[i])>lim || abs(HH[i])>lim || abs(FF[i])>lim || numtype(RX[i]+RH[i]+RF[i]))	// delete all undesirable points
			DeletePoints i, 1, XX,HH,FF, RX,RH,RF
			DeletePoints i, 1, totalAngles, totalIntensity, totalPeakIntensity
			DeletePoints i, 1, imageNames
		endif
	endfor
	N=numpnts(XX)
	printf "at end N=%d,   took %.2f sec\r",N,stopMSTimer(timer)*1e-6
End




Function/S Rotation3DfromR3dXYZ(rotation,epsilonIsZero)	// calculate a volume for gizmo from Rodriques vectors
	String rotation				// specifys the desired rotation axis, can be one of "X;H;F;Y;Z;R_Cyl;Phi_Cyl;Axis_Cyl;Total"
	Variable epsilonIsZero		// flag, true means ignore epsilon when computing alpha and GND (i.e. assume that epsilon is zero)

	if (strlen(rotation)<1 || numtype(epsilonIsZero))
		epsilonIsZero = numtype(epsilonIsZero) ? 1 : epsilonIsZero
		epsilonIsZero = epsilonIsZero ? 1 : 2		// setup for prompt's popup
		rotation = SelectString(strlen(rotation),"X",rotation)
		Prompt rotation,"rotation type",popup,listOfAllRotations
		Prompt epsilonIsZero,"is epsilon zero",popup,"is zero;non-zero"
		DoPrompt "3D rotation type"rotation,epsilonIsZero
		if (V_flag)
			return ""
		endif
		epsilonIsZero = (epsilonIsZero==1)
	endif
	Wave R3dX=R3dX, R3dH=R3dH, R3dF=R3dF
	Wave Rvol=R3dX0
	if (!WaveExists(Rvol))
		Wave Rvol=R3dX
	endif
	if (!WaveExists(Rvol))
		DoAlert 0, "Cannot find R3dX, did you forget to interpolate?"
		return ""
	endif

	String errstr=""
	if (!WaveExists(R3dX) || !WaveExists(R3dH) || !WaveExists(R3dF))
	elseif (WhichListItem(rotation,listOfAllRotations)<0)
		errstr = "Rotation3DfromR3dXYZ(), rotation = '"+rotation+"'"
	elseif (numtype(epsilonIsZero))
		errstr = "Rotation3DfromR3dXYZ(), epsilonIsZero = "+num2str(epsilonIsZero)
	endif
	if (strlen(errstr))
		DoAlert 0, errstr
		return ""
	endif

	Wave totalPeakIntensity3d=totalPeakIntensity3d, totalIntensity3d=totalIntensity3d
	if (WaveExists(totalPeakIntensity3d))
		Wave intensWave = totalPeakIntensity3d
	elseif (WaveExists(totalIntensity3d))
		Wave intensWave = totalIntensity3d
	else
		Wave intensWave = $""
	endif

	String GizmoWaveUnit
	Variable ctensor=0, strain=0						// flags, need lattice curvature tensor, or strain tensor (epsilon)
	Variable grad=0, calcAxis=0							// flags, neeed gradient of rotation angle, or need to calculate the rotation axis normal (cylindrical coordinates)
	Variable ddtensor=0									// flag, need the dislocation density tensor
	Variable aRotation=0								// flags, just a plain rotation
	Variable i0,j0										// specifies a particular component
	if (stringmatch(rotation,"GND"))
		ctensor = 1										// need the curvature tensor
		strain = 1										// needs alpha, so needs both kappa and epsilon
		ddtensor = 1									// and needs the dislocation density tensor
		i0 = -1
		GizmoWaveUnit = "dislocations / cm\\S2"
	elseif (strsearch(rotation,"alpha(",0)==0)
		ctensor = 1										// need the dislocation density tensor, so need curvature tensor
		strain = 1										// and needs epsilon
		ddtensor = 1									// and needs the dislocation density tensor
		if (ArrayOf3dOrients#xyComponent2ij(rotation,i0,j0))			// get i0,j0 index into alpha[i0][j0]
			return ""
		endif
		GizmoWaveUnit = "radian / µm"
	elseif (strsearch(rotation,"kappa(",0)==0)
		ctensor = 1										// this cut will neede the lattice curvature tensor
		if(ArrayOf3dOrients#xyComponent2ij(rotation,i0,j0))			// get i0,j0 index into kappa[i0][j0]
			return ""
		endif
		GizmoWaveUnit = "radian / µm"
	elseif (strsearch(rotation,"dRd",0)==0 || stringmatch(rotation,"|gradR|"))
		grad = 1										// this cut will need the gradient of the total rotation
		i0 = WhichListItem(rotation,"dRdX;dRdY;dRdZ;|gradR|")
		if (i0<0)
			return ""
		endif
		i0 = i0==3 ? -1 : i0							// use -1 to flag magnitude
		GizmoWaveUnit = "radian / µm"
	elseif(WhichListItem(rotation,"exx;eyy;ezz;exy;exz;eyz")>=0)
		strain = 1										// need the strain tensor
		GizmoWaveUnit = ""								// no unit for epsilon
		Wave e3d=$(rotation+"3d")
		if (!WaveExists(e3d))
			return ""
		endif
	elseif(WhichListItem(rotation,"a;b;c")>=0)
		strain = 1
		GizmoWaveUnit = "nm"							// lattice parameter (units=nm)
		Wave e3d=$(rotation+"3d")
		if (!WaveExists(e3d))
			return ""
		endif
	elseif(WhichListItem(rotation,"alpha;beta;gamma")>=0)
		strain = 1
		GizmoWaveUnit = "¡"							// lattice parameter (units=¡)
		Wave e3d=$(rotation+"3d")
		if (!WaveExists(e3d))
			return ""
		endif
	elseif (!stringmatch(rotation,"Total"))				// need to compute angle component
		aRotation = 1
		Wave rotationAxis = $MakeDirectionVector(rotation)// want rotation component about rotationAxis[]
		calcAxis = numtype(rotationAxis[0])!=0		// axis of rotation is not constant, needs to be re-calculated for each point
		GizmoWaveUnit = "radian"
	else
		aRotation = 1
		GizmoWaveUnit = "radian"
	endif

	Duplicate/O Rvol,Rgizmo							// Rgizmo holds the result

	String noteStr = note(Rvol)
	Variable x0,h0,f0, dXaxis, dHaxis, dFaxis
	x0 = NumVarOrDefault(":X0",NumberByKey("X0",noteStr,"=")) //origin to work from (especially important for cylindrical coordinates)
	h0 = NumVarOrDefault(":H0",NumberByKey("H0",noteStr,"="))
	f0 = NumVarOrDefault(":F0",NumberByKey("F0",noteStr,"="))
	dXaxis = NumVarOrDefault(":dXaxis",NumberByKey("dX",noteStr,"=")) // direction of the axis for cylindrical coordinates
	dHaxis = NumVarOrDefault(":dHaxis",NumberByKey("dH",noteStr,"="))
	dFaxis = NumVarOrDefault(":dFaxis",NumberByKey("dF",noteStr,"="))
	Wave axis=$MakeUnique3Vector($"")				// axis of cylindrical coords, needed to calculate rotationAxis[]
	axis = {dXaxis,dHaxis,dFaxis}						// axis is assumed to go through xhf0[]
	normalize(axis)

	//	We are assuming that xhf0[] in the 3d volume is the origin that we want.
	//
	Wave xhf0=$MakeUnique3Vector($"")				// origin
	Wave xhf=$MakeUnique3Vector($"")				// a point in Rgizmo, XHF space
	Wave Rod=$MakeUnique3Vector($"")				// direction of Rodriques vector
	Wave rvec=$MakeUnique3Vector($"")				// vector from xhf0[] to point in Rgizmo
	Wave rhat=$MakeUnique3Vector($"")				// like xhf, but axial part removed
	xhf0 =  {x0,h0,f0}
	if (ItemsInList(GetRTStackInfo(0))<2)
		printf "Rotation3DfromR3dXYZ(\"%s\",%d)		// rotationAxis=\"%s\", epsilonIsZero=%d\r",rotation,epsilonIsZero,rotation,epsilonIsZero
		if (!calcAxis && aRotation && WaveExists(rotationAxis))
			printf " the rotation axis = (%g, %g, %g)\r",rotationAxis[0],rotationAxis[1],rotationAxis[2]
		endif
		printf "xhf0={%.3g, %.3g, %.3g},     axis={%.3g, %.3g, %.3g}\r",xhf0[0],xhf0[1],xhf0[2],axis[0],axis[1],axis[2]
	endif
	noteStr = ReplaceNumberByKey("epsilonIsZero",noteStr,epsilonIsZero,"=")
	noteStr = ReplaceStringByKey("rotation",noteStr,rotation,"=")
	Note/K Rgizmo, ReplaceStringByKey("waveClass",noteStr,"OrientationGizmoWave","=")
	SetScale d 0,0,GizmoWaveUnit, Rgizmo
//	Note/K sliceWaveIntens, ReplaceStringByKey("waveClass",noteStr,"InensitySliceWave","=")
//	Note/K sliceOrientsD, ReplaceStringByKey("waveClass",noteStr,"RGBfullOrientationMap","=")
//	SetScale d 0,0,"radian", sliceOrientsD
//	SetScale d 0,0,"", sliceWaveIntens
//	sliceWaveIntens = NaN
//	sliceOrientsD = NaN

	// next fill the Rgizmo wave
	Variable i,j,k, xx,yy,zz, angle, GND, magGrad, peakIntensity
	Make/N=(3,3)/D/O SlicePlaneIn3d_alpha, SlicePlaneIn3d_kappa
	Make/N=(3)/O/D SlicePlaneIn3d_grad
	Wave alpha=SlicePlaneIn3d_alpha, kappa=SlicePlaneIn3d_kappa, gradR=SlicePlaneIn3d_grad
	Variable ii,jj,kk													// indicies into R3dX (can be different than R3dX0)
	Variable dot
	Variable Nx=DimSize(Rgizmo,0), Ny=DimSize(Rgizmo,1), Nz=DimSize(Rgizmo,2)
	for (k=0;k<Nz;k+=1)
		zz = DimOffset(Rgizmo,2) + k*DimDelta(Rgizmo,2)
		kk = (zz-DimOffset(R3dX,2))/DimDelta(R3dX,2)
		for (j=0;j<Ny;j+=1)
			yy = DimOffset(Rgizmo,1) + j*DimDelta(Rgizmo,1)
			jj = (yy-DimOffset(R3dX,1))/DimDelta(R3dX,1)
			for (i=0;i<Nx;i+=1)
				xx = DimOffset(Rgizmo,0) + i*DimDelta(Rgizmo,0)
				ii = (xx-DimOffset(R3dX,0))/DimDelta(R3dX,0)
				xhf = {xx,yy,zz}
	
				Rod[0] = R3dX[ii][jj][kk]							// Rodriques vector at xhf
				Rod[1] = R3dH[ii][jj][kk]
				Rod[2] = R3dF[ii][jj][kk]

				peakIntensity = WaveExists(intensWave) ? intensWave[ii][jj][kk] : 1
//				sliceWaveIntens[i][j] = peakIntensity
//				sliceOrientsD[i][j][] = Rod[r]

				if (aRotation)											// looking for an angle, not a dislocation denstity tensor component
					// compute the angle for this point and save it in sliceWave[i][j]
					angle = 2*atan(normalize(Rod))					// Rodriques angle (radian), Rod is now normalized
					if (calcAxis)										// need to recalculate rotationAxis[] each time
						rvec = xhf - xhf0								// rvec is now direction of point relative to the origin of cylinder
						dot = MatrixDot(axis,rvec)
						rhat = rvec - dot*axis							// put rhat in the plane where z along axis is 0, only r & phi
						normalize(rhat)
						if (stringmatch(rotation,"R_Cyl"))				// rotation around axis sticking out radially from axis, r^ direction
							rotationAxis = rhat
						elseif (stringmatch(rotation,"Phi_Cyl"))		// rotation around axis in phi^ direction, probably the biggest for an indent
							Cross axis, rhat								// rotationAxis = axis x rHat - phi^
							Wave W_Cross=W_Cross
							rotationAxis = W_Cross
						elseif (stringmatch(rotation,"Axis_Cyl"))		// rotation around axis of cylinder, circulating direction
							rotationAxis = axis
						else
							Abort "SlicePlaneIn3d(), invalid rotation"	// should not be possible to reach this line
						endif
					endif
	
					dot = stringmatch(rotation,"Total") ? 1 : MatrixDot(rotationAxis,Rod)
					Rgizmo[i][j][k] = angle * dot						// amount of rotation about rotationAxis[]
	
				elseif (ddtensor)										// GND or alpha component, need a dislocation density tensor component, or GND
					curvatureTensorFromRods3d(xhf[0],xhf[1],xhf[2],R3dX,R3dH,R3dF,kappa)
					if (epsilonIsZero)
						GND = dislocationTensorFromRods3d(xhf[0],xhf[1],xhf[2],$"",$"",$"",$"",$"",$"",kappa,alpha)
					else
						GND = dislocationTensorFromRods3d(xhf[0],xhf[1],xhf[2],exx3d,eyy3d,ezz3d,exy3d,exz3d,eyz3d,kappa,alpha)
					endif
					Rgizmo[i][j][k] = (i0<0 || j0<0) ? GND : alpha[i0][j0]	// component of alpha(i0,j0) or GND
	
				elseif (ctensor)													// kappa, a lattice curvature tensor component
					GND = curvatureTensorFromRods3d(xhf[0],xhf[1],xhf[2],R3dX,R3dH,R3dF,kappa) // this GND assumes no strain
					Rgizmo[i][j][k] = (i0<0 || j0<0) ? GND : kappa[i0][j0]	// component of kappa(i0,j0) or GND (GND, without strain part)
	
				elseif (strain)
					Rgizmo[i][j][k] = Interp3d(e3d, xhf[0],xhf[1],xhf[2])
	
				elseif (grad)														// need a part of the gradient of the total rotation angle
					// calculate the gradient of the total rotation angle at the point (xx,hh,ff) from the 3d-Rodriques vectors
					magGrad = GradtRotationAnglesFromRods3d(xhf[0],xhf[1],xhf[2],R3dX,R3dH,R3dF,gradR)
					Rgizmo[i][j][k] = (i0<0) ? magGrad : gradR[i0]			// d|R| / dx,    or dy, or dz   or | gradient (R) |
				endif
			endfor
		endfor
	endfor
//	Rgizmo = numtype(Rgizmo) ? 0 : Rgizmo
//	ArrayOf3dOrients#RodriquesToRGB(sliceOrientsD,"sliceOrients")
//	ArrayOf3dOrients#sliceWave2RGB(sliceWave,sliceWaveIntens)
//	KillWaves/Z sliceOrientsD
	KillWaves/Z SlicePlaneIn3d_alpha, SlicePlaneIn3d_kappa, SlicePlaneIn3d_grad
	KillWaves/Z rotationAxis,axis,Rod,xhf0,xhf,rhat,rvec,W_Cross
	return GetWavesDataFolder(Rgizmo,2)
End












Function MakeScatter(imax,rmax)
	Variable imax				// imax=8e5
	Variable rmax				// 0.01
	Variable Fsurf=-7
	Variable abslen = 12		// absorption length
	Variable zoom=1			// 3,  zoom about the weighted center

	Wave XX=XX, HH=HH, FF=FF, RX=RX, RH=RH, RF=RF
//	Wave intensity=totalPeakIntensity		//total Peak Intensity
	Wave intensity=totalIntensity			//totalIntensity
	Variable N=numpnts(XX)

	String axis="X"
	String colorTable="RedWhiteBlue256"
	Prompt axis,"rotation axis",popup,"X;H;F"
	Prompt colorTable,"color table",popup,CTabList()
	DoPrompt "color table",colorTable,axis
	if (V_flag)
		return 1
	endif
	if (stringmatch(axis,"X"))
		Wave RR = RX
	elseif (stringmatch(axis,"H"))
		Wave RR = RH
	elseif (stringmatch(axis,"F"))
		Wave RR = RF
	else
		DoAlert 0,"rotation axis '"+axis+"' is impossible"
		return 1
	endif
	WaveStats/M=1/Q RR
	Variable RRmin = V_max/50
	printf "max rotation is %g (rad) = %g¡\r",V_max,V_max*180/PI
	printf "min rotation is %g (rad) = %g¡\r",RRmin,RRmin*180/PI

	ColorTab2Wave $colorTable
	Wave M_colors=M_colors
	Variable ic=DimSize(M_colors,0)-1			// highest color value (probably 255 here)

	Variable rmin = -rmax
	printf "colors cover the rotatation range R%s[%g,  %g](rad)\r",axis,rmin,rmax
	Make/N=(N,3)/O xyz=NaN						// position of each point
	Make/N=(N,4)/O xyzRGBA=NaN					// RGBA for each point
	Make/N=(N,3)/O xyzSize=NaN					// size for each point
	Make/N=(N,4)/O xyzRot=NaN					// direction vector for each point

	Variable i,j,m,intens, markerSize
	for (i=0,j=0; i<N; i+=1)
		intens = limit(intensity[i]/imax,0,1)
		intens *= (FF[i] > Fsurf) ? exp((FF[i]-Fsurf)/abslen) : 1	// correct for absorption
		markerSize = intens/6
		if (markerSize<.05 || !(abs(RR[i])>=RRmin))
			continue
		endif

		xyz[j][0] = XX[i]
		xyz[j][1] = HH[i]
		xyz[j][2] = FF[i]
		m = ic*(RR[i]-rmin)/(rmax-rmin)
		m = round(limit(m,0,ic))
		xyzRGBA[j][0,2] = M_colors[m][q]/65536

		xyzRGBA[j][3] = intens/2
		xyzSize[j][] = markerSize

		xyzRot[j][0] = 2*atan(sqrt(RX[i]^2 + RH[i]^2 + RF[i]^2)) * 180/PI	// angle in degrees
		xyzRot[j][1] = RX[i]
		xyzRot[j][2] = RH[i]
		xyzRot[j][3] = RF[i]
		j += 1
	endfor
	KillWaves/Z M_colors

	WaveStats/M=1/Q xyzSize						// trim off all un-used points
	N = V_npnts/3 + 8								// eight exta points to set the box size (want it to be square)
	Redimension/N=(N,3) xyz, xyzSize
	Redimension/N=(N,4) xyzRGBA, xyzRot
	Variable i8=N-8								// first index of the last eight points
	xyzRGBA[i8,N-1][] = 0
	xyzSize[i8,N-1][] = 0
	xyzRot[i8,N-1][] = 0

	// find the range of xyz, and set last 8 points to make the volume have equal length edges
	Make/N=(N-8)/O rangeWaveScatter
	rangeWaveScatter = xyz[p][0]					// find the x-range
	WaveStats/M=1/Q rangeWaveScatter
	Variable xRange = V_max-V_min, xmid=(V_max+V_min)/2
	rangeWaveScatter = xyz[p][1]					// find the y-range
	WaveStats/M=1/Q rangeWaveScatter
	Variable yRange = V_max-V_min, ymid=(V_max+V_min)/2
	rangeWaveScatter = xyz[p][2]					// find the z-range
	WaveStats/M=1/Q rangeWaveScatter
	Variable zRange = V_max-V_min, zmid=(V_max+V_min)/2
	KillWaves/Z rangeWaveScatter
	Variable l2 = max(xRange,yRange)				// half the  longest length
	l2 = max(l2,zRange)/2

	if (zoom!=1)
		// find the weighted center of the sample
		Make/N=(N-8)/O weightedCenterScatter
		Variable xc,yc,zc
		weightedCenterScatter = xyz[p][0]*xyzSize[p][0]
		WaveStats/M=1/Q weightedCenterScatter
		xc = V_npnts*V_avg
		weightedCenterScatter = xyzSize[p][0]
		WaveStats/M=1/Q weightedCenterScatter
		xc /= (V_npnts*V_avg)

		weightedCenterScatter = xyz[p][1]*xyzSize[p][1]
		WaveStats/M=1/Q weightedCenterScatter
		yc = V_npnts*V_avg
		weightedCenterScatter = xyzSize[p][1]
		WaveStats/M=1/Q weightedCenterScatter
		yc /= (V_npnts*V_avg)

		weightedCenterScatter = xyz[p][2]*xyzSize[p][2]
		WaveStats/M=1/Q weightedCenterScatter
		zc = V_npnts*V_avg
		weightedCenterScatter = xyzSize[p][2]
		WaveStats/M=1/Q weightedCenterScatter
		zc /= (V_npnts*V_avg)
		printf "weighted center is [%g,  %g,  %g] (µm)\r",xc,yc,zc
		KillWaves/Z weightedCenterScatter

		l2 /= zoom
		xmid = xc
		ymid = yc
		zmid = zc
		for (i=0;i<N-8;i+=1)			// now remove all points outside of this box
			if (xyz[i][0]<(xc-l2) || xyz[i][0]>(xc+l2) || xyz[i][1]<(yc-l2) || xyz[i][1]>(yc+l2) || xyz[i][2]<(zc-l2) || xyz[i][2]>(zc+l2))
			xyz[i][] = NaN
			xyzRGBA[i] = NaN
			xyzSize[i][] = NaN
			xyzRot[i][] = NaN
			endif
		endfor
		xyzSize *= sqrt(zoom)
	endif

	xyz[i8+0,i8+3][2] = zmid-l2
	xyz[i8+4,i8+7][2] = zmid+l2

	xyz[i8+0,i8+1][1] = ymid-l2
	xyz[i8+2,i8+3][1] = ymid+l2
	xyz[i8+4,i8+5][1] = ymid-l2
	xyz[i8+6,i8+7][1] = ymid+l2

	xyz[i8+0][0] = xmid-l2
	xyz[i8+1][0] = xmid+l2
	xyz[i8+2][0] = xmid-l2
	xyz[i8+3][0] = xmid+l2
	xyz[i8+4][0] = xmid-l2
	xyz[i8+5][0] = xmid+l2
	xyz[i8+6][0] = xmid-l2
	xyz[i8+7][0] = xmid+l2
	print "Z of plane =",(7-zmid)/l2

	return 0
End




Function MakeGizmoRot(Rgizmo,[valueList])	// make a gizmo of a
	Wave Rgizmo
	String valueList		// = "0.015;-0.006;0.008;-0.004;0.03"
	if (ParamIsDefault(valueList ))
		valueList = "0.015;-0.006;0.008;-0.004;0.03"
	endif

	if(exists("NewGizmo")!=4)				// Do nothing if the Gizmo XOP is not available.
		DoAlert 0, "Gizmo XOP must be installed"
		return 1
	endif
	String Rname = ""
	if (!WaveExists(Rgizmo))
		String wList = WaveListClass("OrientationGizmoWave","*","DIMS:3")
		if (ItemsInList(wList)==1)
			Rname = StringFromList(0,wList)	// only one choice, so do not ask
		else
			Prompt Rname,"3D wave for Gizmo",popup,WaveListClass("OrientationGizmoWave","*","DIMS:3")
			DoPrompt "3D wave",Rname
			if (V_flag)
				return 1
			endif
		endif
		Wave Rgizmo = $Rname
	endif
	if (!WaveExists(Rgizmo))
		DoAlert 0,"Rgizmo does not exist"
		return 1
	endif
	Rname = GetWavesDataFolder(Rgizmo,2)
	String wnote = note(Rgizmo)
	String rotation = StringByKey("rotation", wnote,"=")
	String title=IgorInfo(1)+",  "+GetDataFolder(0)+",  "+rotation
	printf "using values of \"%s\"\r",valueList
	String colorList="0.25,0.25,1,  0.5;1,0.5,0.5,  0.5;0,0,0.8,  1;1,0,0,  1;0,1,1,  0.5"

	Execute "NewGizmo/N=GizmoRot/T=\"GizmoRot\" /W=(40,56,800,728)"
	Execute "ModifyGizmo startRecMacro"
	String titleObjectName = addFixedTitleStringToGizmo("TL",title)
	Execute "AppendToGizmo attribute shininess={5,1032},name=shininess0"
	Execute "AppendToGizmo attribute blendFunc={770,771},name=blendFunc0"
	Execute "AppendToGizmo attribute specular={0.1,0.1,0.1,1,1032},name=specular0"
	Execute "AppendToGizmo attribute ambient={0.1,0.1,0.1,1,1032},name=ambient0"
	Execute "AppendToGizmo attribute color={1,1,1,1},name=colorWhite"

	Execute "AppendToGizmo light=Directional,name=light0"
	Execute "ModifyGizmo light=light0 property={ position,0.260909,0,-0.965363,0}"
//	Execute "ModifyGizmo light=light0 property={ direction,0.260909,0,-0.965363}"
//	Execute "ModifyGizmo light=light0 property={ direction,0.939693,-0.000001,-0.342020}"
	Execute "ModifyGizmo light=light0 property={ direction,0.984808,-0.000001,-0.173648}"
	Execute "ModifyGizmo light=light0 property={ ambient,0,0,0,0.4}"
	Execute "ModifyGizmo light=light0 property={ specular,1,1,1,0.4}"
	Execute "ModifyGizmo light=light0 property={ diffuse,1,1,1,0.4}"
	Execute "AppendToGizmo light=Directional,name=light1"
	Execute "ModifyGizmo light=light1 property={ position,0.987580,0,0.157115,0}"
//	Execute "ModifyGizmo light=light1 property={ direction,0.987580,0,0.157115}"
//	Execute "ModifyGizmo light=light1 property={ direction,0.000002,0.939693,0.342020}"
	Execute "ModifyGizmo light=light1 property={ direction,0.000002,0.984808,0.173648}"
	Execute "ModifyGizmo light=light1 property={ ambient,0,0,0,0.4}"
	Execute "ModifyGizmo light=light1 property={ specular,1,1,1,0.4}"
	Execute "ModifyGizmo light=light1 property={ diffuse,1,1,1,0.4}"
	Execute "AppendToGizmo Axes=boxAxes,name=axes0"
	setAxesInGizmo("X  (µm)","H  (µm)","F  (µm)")

	Variable X0=NumberByKey("X0",wnote,"="), H0=NumberByKey("H0",wnote,"=")
	if (numtype(X0+H0)==0)
		String axisName = GetWavesDataFolder(Rgizmo,1)+"axisPathGizmo"
		Make/N=(2,3)/O $axisName
		Wave axisPath=$axisName
		axisPath[][0] = X0
		axisPath[][1] = H0
		axisPath[][2] = {DimOffset(Rgizmo,2),DimOffset(Rgizmo,2)+DimDelta(Rgizmo,0)*(DimSize(Rgizmo,2)-1)}

		Execute "AppendToGizmo Path="+axisName+",name=pathAxis"
		Execute "ModifyGizmo ModifyObject=pathAxis property={ pathColorType,1}"
		Execute "ModifyGizmo ModifyObject=pathAxis property={ lineWidthType,1}"
		Execute "ModifyGizmo ModifyObject=pathAxis property={ lineWidth,3}"
		Execute "ModifyGizmo ModifyObject=pathAxis property={ pathColor,1,0.25,0.85,1}"
	endif

	String boxName = GetWavesDataFolder(Rgizmo,1)+"boxScatterGizmo"		// append a box to make the volume Rgizmo scaled uniformly
	Make/N=(8,3)/O $boxName
	Wave boxScatter=$boxName
	Variable Xstart=DimOffset(Rgizmo,0), Hstart=DimOffset(Rgizmo,1), Fstart=DimOffset(Rgizmo,2)
	Variable Xwid=DimDelta(Rgizmo,0)*(DimSize(Rgizmo,0)-1), Hwid=DimDelta(Rgizmo,1)*(DimSize(Rgizmo,1)-1), Fwid=DimDelta(Rgizmo,2)*(DimSize(Rgizmo,2)-1)
	Variable dw, width=max(max(Xwid,Hwid),Fwid)
	dw = max(width-Xwid,0)								// extend to width, this will make a cube
	Xstart -= dw/2 ;	Xwid += dw
	dw = max(width-Hwid,0)
	Hstart -= dw/2 ;	Hwid += dw
	dw = max(width-Fwid,0)
	Fstart -= dw/2 ;	Fwid += dw
	boxScatter[][0] = {Xstart,Xstart+Xwid,Xstart+Xwid,Xstart,Xstart,Xstart+Xwid,Xstart+Xwid,Xstart}
	boxScatter[][1] = {Hstart,Hstart,Hstart+Hwid,Hstart+Hwid,Hstart,Hstart,Hstart+Hwid,Hstart+Hwid}
	boxScatter[][2] = {Fstart,Fstart,Fstart,Fstart,Fstart+Fwid,Fstart+Fwid,Fstart+Fwid,Fstart+Fwid}
	Execute "AppendToGizmo Scatter="+boxName+",name=scatterBox"
	Execute "ModifyGizmo ModifyObject=scatterBox property={ markerType,0}"
	Execute "ModifyGizmo ModifyObject=scatterBox property={ sizeType,0}"
	Execute "ModifyGizmo ModifyObject=scatterBox property={ Shape,1}"
	Execute "ModifyGizmo ModifyObject=scatterBox property={ size,1}"

	Variable dx=0, dh=0
	if (numtype(X0+H0)==0)						// put axis cue through X0,H0 if they exist
		dx = (X0 - Xstart-Xwid/2) * 2/Xwid
		dh = (H0 - Hstart-Hwid/2) * 2/Hwid
	endif
	Execute "AppendToGizmo freeAxesCue={"+num2str(dx)+","+num2str(dh)+",0,1},name=freeAxesCue0"

	String color,isoList=""
	Variable i, value

	Execute "AppendToGizmo group,name=SurfaceGroup0"
	// ************************* Group Object Start *******************
	Execute "ModifyGizmo currentGroupObject=\"SurfaceGroup0\""
	for (i=0;i<itemsInList(valueList);i+=1)		// add in reverse order (biggest to smallest)
		value = str2num(StringFromList(i,valueList))
		color = StringFromList(i, colorList)
		isoList += Add_isoSurface(i,value,color,Rname)+";"
	endfor
	for (i=itemsInList(isoList)-1;i>=0;i-=1)
		Execute "ModifyGizmo setDisplayList=-1, object="+StringFromList(i,isoList)
	endfor
	Execute "ModifyGizmo currentGroupObject=\"::\""
	// ************************* Group Object End *******************

	Execute "ModifyGizmo setDisplayList=0, object="+titleObjectName
	Execute "ModifyGizmo setDisplayList=1, opName=MainTransform, operation=mainTransform"
	Execute "ModifyGizmo setDisplayList=2, attribute=blendFunc0"
	Execute "ModifyGizmo setDisplayList=3, object=light0"
	Execute "ModifyGizmo setDisplayList=4, object=light1"
	Execute "ModifyGizmo setDisplayList=5, attribute=shininess0"
	Execute "ModifyGizmo setDisplayList=6, attribute=ambient0"
	Execute "ModifyGizmo setDisplayList=7, attribute=specular0"
	Execute "ModifyGizmo setDisplayList=8, opName=rotateView, operation=rotate, data={0,0,0,0}"
	Execute "ModifyGizmo setDisplayList=9, object=scatterBox"
	if (numtype(X0+H0)==0)
		Execute "ModifyGizmo setDisplayList=-1, object=pathAxis"
		Execute "ModifyGizmo setDisplayList=-1, attribute=colorWhite"
	endif
	Execute "ModifyGizmo setDisplayList=-1, object=SurfaceGroup0"
	Execute "ModifyGizmo setDisplayList=-1, object=axes0"
	Execute "ModifyGizmo setDisplayList=-1, object=freeAxesCue0"

	Execute "ModifyGizmo SETQUATERNION={-0.192446,0.565210,-0.720302,0.353137}"
	Execute "ModifyGizmo autoscaling=1"
	Execute "ModifyGizmo currentGroupObject=\"\""
	Execute "ModifyGizmo compile"
	Execute "ModifyGizmo showInfo"
//	Execute "ModifyGizmo infoWindow={793,92,1268,548}"
	Execute "ModifyGizmo infoWindow={794,92,1270,648}"
	Execute "ModifyGizmo bringToFront"
//	Execute "ModifyGizmo showAxisCue=1"
	Execute "ModifyGizmo endRecMacro"

	if (!isMacBookPro() && !stringmatch(systemUserName(),"tischler"))			// all done if not a MacBookPro
		return 0
	endif
	DoAlert 1,"Convert Iso-Surfaces to triangle waves?"
	if (V_flag!=1)
		return 0
	endif

	Execute "ModifyGizmo startRecMacro"
	// make the triangle waves so things display properly on the MacBookPro
	wnote = ReplaceStringByKey("waveClass",wnote,"OrientationGizmoTri","=")
	String item,triList=""
	Execute "AppendToGizmo group,name=SurfaceGroup0_Tri"
	// ************************* Group Object Start *******************
	Execute "ModifyGizmo currentGroupObject=\"SurfaceGroup0_Tri\""
	for (i=itemsInList(isoList)-1;i>=0;i-=1)
		KillWaves/Z MakeGizmoRot_tri
		item = StringFromList(i,isoList)
		Execute "ModifyGizmo currentGroupObject=\"::\""
		Execute "ModifyGizmo currentGroupObject=\"SurfaceGroup0\""
		Execute "ModifyGizmo modifyObject="+item+",property={savetowave,\"MakeGizmoRot_tri\"}"
		Execute "ModifyGizmo currentGroupObject=\"::\""
		Execute "ModifyGizmo currentGroupObject=\"SurfaceGroup0_Tri\""
		if (exists("MakeGizmoRot_tri")!=1)
			continue
		endif
		Duplicate/O MakeGizmoRot_tri $(item+"_tri")
		Wave triW=$(item+"_tri")
		Note/K triW, wnote
		KillWaves/Z MakeGizmoRot_tri
		color = GizmoObjectProperty("",item,"frontColor")// get color of the isoSurface
		triList += Add_triangleSurface(item+"_tri",GetWavesDataFolder(triW,2),color)+";"
	endfor
	for (i=0;i<itemsInList(triList);i+=1)				// add the triangle waves to display list
		Execute "ModifyGizmo insertDisplayList=-1, object="+StringFromList(i,triList)
	endfor
	Execute "ModifyGizmo currentGroupObject=\"::\""
	// ************************* Group Object End *******************

	Execute "GetGizmo/Z displayNameList"					// find the index into display list of the iso-surfaces group
	SVAR S_DisplayNames=S_DisplayNames
	Variable location = WhichListItem("SurfaceGroup0",S_DisplayNames)
	KillStrings/Z S_DisplayNames
	Execute "ModifyGizmo insertDisplayList="+num2istr(location)+", object=SurfaceGroup0_Tri"	// insert triangle wave group
	Execute "RemoveFromGizmo/Z displayItem=SurfaceGroup0"										// delete iso-surface wave group
	Execute "ModifyGizmo endRecMacro"
	return 0
End
//Function MakeGizmoRot(Rgizmo,[valueList])	// make a gizmo of a
//	Wave Rgizmo
//	String valueList		// = "0.015;-0.006;0.008;-0.004;0.03"
//	if (ParamIsDefault(valueList ))
//		valueList = "0.015;-0.006;0.008;-0.004;0.03"
//	endif
//
//	if(exists("NewGizmo")!=4)				// Do nothing if the Gizmo XOP is not available.
//		DoAlert 0, "Gizmo XOP must be installed"
//		return 1
//	endif
//	String Rname = ""
//	if (!WaveExists(Rgizmo))
//		String wList = WaveListClass("OrientationGizmoWave","*","DIMS:3")
//		if (ItemsInList(wList)==1)
//			Rname = StringFromList(0,wList)	// only one choice, so do not ask
//		else
//			Prompt Rname,"3D wave for Gizmo",popup,WaveListClass("OrientationGizmoWave","*","DIMS:3")
//			DoPrompt "3D wave",Rname
//			if (V_flag)
//				return 1
//			endif
//		endif
//		Wave Rgizmo = $Rname
//	endif
//	if (!WaveExists(Rgizmo))
//		DoAlert 0,"Rgizmo does not exist"
//		return 1
//	endif
//	Rname = GetWavesDataFolder(Rgizmo,2)
//	String wnote = note(Rgizmo)
//	String rotation = StringByKey("rotation", wnote,"=")
//	String title=IgorInfo(1)+",  "+GetDataFolder(0)+",  "+rotation
//	printf "using values of \"%s\"\r",valueList
//	String colorList="0.25,0.25,1,  0.5;1,0.5,0.5,  0.5;0,0,0.8,  1;1,0,0,  1;0,1,1,  0.5"
//
//	Execute "NewGizmo/N=GizmoRot/T=\"GizmoRot\" /W=(40,56,800,728)"
//	Execute "ModifyGizmo startRecMacro"
//	String titleObjectName = addFixedTitleStringToGizmo("TL",title)
//	Execute "AppendToGizmo attribute shininess={5,1032},name=shininess0"
//	Execute "AppendToGizmo attribute blendFunc={770,771},name=blendFunc0"
//	Execute "AppendToGizmo attribute specular={0.1,0.1,0.1,1,1032},name=specular0"
//	Execute "AppendToGizmo attribute ambient={0.1,0.1,0.1,1,1032},name=ambient0"
//	Execute "AppendToGizmo attribute color={1,1,1,1},name=colorWhite"
//	Execute "AppendToGizmo light=Directional,name=light0"
//	Execute "ModifyGizmo light=light0 property={ position,0.260909,0,-0.965363,0}"
//	Execute "ModifyGizmo light=light0 property={ direction,0.260909,0,-0.965363}"
//	Execute "ModifyGizmo light=light0 property={ ambient,0,0,0,0.4}"
//	Execute "ModifyGizmo light=light0 property={ specular,1,1,1,0.4}"
//	Execute "ModifyGizmo light=light0 property={ diffuse,1,1,1,0.4}"
//	Execute "AppendToGizmo light=Directional,name=light1"
//	Execute "ModifyGizmo light=light1 property={ position,0.987580,0,0.157115,0}"
//	Execute "ModifyGizmo light=light1 property={ direction,0.987580,0,0.157115}"
//	Execute "ModifyGizmo light=light1 property={ ambient,0,0,0,0.4}"
//	Execute "ModifyGizmo light=light1 property={ specular,1,1,1,0.4}"
//	Execute "ModifyGizmo light=light1 property={ diffuse,1,1,1,0.4}"
//	Execute "AppendToGizmo Axes=boxAxes,name=axes0"
//	setAxesInGizmo("X  (µm)","H  (µm)","F  (µm)")
//
//	Variable X0=NumberByKey("X0",wnote,"="), H0=NumberByKey("H0",wnote,"=")
//	if (numtype(X0+H0)==0)
//		String axisName = GetWavesDataFolder(Rgizmo,1)+"axisPathGizmo"
//		Make/N=(2,3)/O $axisName
//		Wave axisPath=$axisName
//		axisPath[][0] = X0
//		axisPath[][1] = H0
//		axisPath[][2] = {DimOffset(Rgizmo,2),DimOffset(Rgizmo,2)+DimDelta(Rgizmo,0)*(DimSize(Rgizmo,2)-1)}
//
//		Execute "AppendToGizmo Path="+axisName+",name=pathAxis"
//		Execute "ModifyGizmo ModifyObject=pathAxis property={ pathColorType,1}"
//		Execute "ModifyGizmo ModifyObject=pathAxis property={ lineWidthType,1}"
//		Execute "ModifyGizmo ModifyObject=pathAxis property={ lineWidth,3}"
//		Execute "ModifyGizmo ModifyObject=pathAxis property={ pathColor,1,0.25,0.85,1}"
//	endif
//
//	String boxName = GetWavesDataFolder(Rgizmo,1)+"boxScatterGizmo"		// append a box to make the volume Rgizmo scaled uniformly
//	Make/N=(8,3)/O $boxName
//	Wave boxScatter=$boxName
//	Variable Xstart=DimOffset(Rgizmo,0), Hstart=DimOffset(Rgizmo,1), Fstart=DimOffset(Rgizmo,2)
//	Variable Xwid=DimDelta(Rgizmo,0)*(DimSize(Rgizmo,0)-1), Hwid=DimDelta(Rgizmo,1)*(DimSize(Rgizmo,1)-1), Fwid=DimDelta(Rgizmo,2)*(DimSize(Rgizmo,2)-1)
//	Variable dw, width=max(max(Xwid,Hwid),Fwid)
//	dw = max(width-Xwid,0)								// extend to width, this will make a cube
//	Xstart -= dw/2 ;	Xwid += dw
//	dw = max(width-Hwid,0)
//	Hstart -= dw/2 ;	Hwid += dw
//	dw = max(width-Fwid,0)
//	Fstart -= dw/2 ;	Fwid += dw
//	boxScatter[][0] = {Xstart,Xstart+Xwid,Xstart+Xwid,Xstart,Xstart,Xstart+Xwid,Xstart+Xwid,Xstart}
//	boxScatter[][1] = {Hstart,Hstart,Hstart+Hwid,Hstart+Hwid,Hstart,Hstart,Hstart+Hwid,Hstart+Hwid}
//	boxScatter[][2] = {Fstart,Fstart,Fstart,Fstart,Fstart+Fwid,Fstart+Fwid,Fstart+Fwid,Fstart+Fwid}
//	Execute "AppendToGizmo Scatter="+boxName+",name=scatterBox"
//	Execute "ModifyGizmo ModifyObject=scatterBox property={ markerType,0}"
//	Execute "ModifyGizmo ModifyObject=scatterBox property={ sizeType,0}"
//	Execute "ModifyGizmo ModifyObject=scatterBox property={ Shape,1}"
//	Execute "ModifyGizmo ModifyObject=scatterBox property={ size,1}"
//
//	Variable dx=0, dh=0
//	if (numtype(X0+H0)==0)						// put axis cue through X0,H0 if they exist
//		dx = (X0 - Xstart-Xwid/2) * 2/Xwid
//		dh = (H0 - Hstart-Hwid/2) * 2/Hwid
//	endif
//	Execute "AppendToGizmo freeAxesCue={"+num2str(dx)+","+num2str(dh)+",0,1},name=freeAxesCue0"
//
//	String color,isoList=""
//	Variable i, value
//	for (i=0;i<itemsInList(valueList);i+=1)		// add in reverse order (biggest to smallest)
//		value = str2num(StringFromList(i,valueList))
//		color = StringFromList(i, colorList)
//		isoList += Add_isoSurface(i,value,color,Rname)+";"
//	endfor
//
//	Execute "ModifyGizmo setDisplayList=0, object="+titleObjectName
//	Execute "ModifyGizmo setDisplayList=1, opName=MainTransform, operation=mainTransform"
//	Execute "ModifyGizmo setDisplayList=2, attribute=blendFunc0"
//	Execute "ModifyGizmo setDisplayList=3, object=light0"
//	Execute "ModifyGizmo setDisplayList=4, object=light1"
//	Execute "ModifyGizmo setDisplayList=5, attribute=shininess0"
//	Execute "ModifyGizmo setDisplayList=6, attribute=ambient0"
//	Execute "ModifyGizmo setDisplayList=7, attribute=specular0"
//	Execute "ModifyGizmo setDisplayList=8, opName=rotateView, operation=rotate, data={0,0,0,0}"
//	Execute "ModifyGizmo setDisplayList=9, object=scatterBox"
//
//	if (numtype(X0+H0)==0)
//		Execute "ModifyGizmo setDisplayList=-1, object=pathAxis"
//		Execute "ModifyGizmo setDisplayList=-1, attribute=colorWhite"
//	endif
//	for (i=itemsInList(isoList)-1;i>=0;i-=1)
//		Execute "ModifyGizmo setDisplayList=-1, object="+StringFromList(i,isoList)
//	endfor
//	Execute "ModifyGizmo setDisplayList=-1, object=axes0"
//	Execute "ModifyGizmo setDisplayList=-1, object=freeAxesCue0"
//
//	Execute "ModifyGizmo SETQUATERNION={-0.192446,0.565210,-0.720302,0.353137}"
//	Execute "ModifyGizmo autoscaling=1"
//	Execute "ModifyGizmo currentGroupObject=\"\""
//	Execute "ModifyGizmo compile"
//	Execute "ModifyGizmo showInfo"
//	Execute "ModifyGizmo infoWindow={793,92,1268,548}"
//	Execute "ModifyGizmo bringToFront"
////	Execute "ModifyGizmo showAxisCue=1"
//	Execute "ModifyGizmo endRecMacro"
//
//	if (!isMacBookPro() && !stringmatch(systemUserName(),"tischler"))			// all done if not a MacBookPro
//		return 0
//	endif
//	DoAlert 1,"Convert Iso-Surfaces to triangle waves?"
//	if (V_flag!=1)
//		return 0
//	endif
//
//	Execute "ModifyGizmo startRecMacro"
//	// make the triangle waves so things display properly on the MacBookPro
//	String item,triList=""
//	wnote = ReplaceStringByKey("waveClass",wnote,"OrientationGizmoTri","=")
//	for (i=itemsInList(isoList)-1;i>=0;i-=1)
//		KillWaves/Z MakeGizmoRot_tri
//		item = StringFromList(i,isoList)
//		Execute "ModifyGizmo modifyObject="+item+",property={savetowave,\"MakeGizmoRot_tri\"}"
//		if (exists("MakeGizmoRot_tri")!=1)
//			continue
//		endif
//		Duplicate/O MakeGizmoRot_tri $(item+"_tri")
//		Wave triW=$(item+"_tri")
//		Note/K triW, wnote
//		KillWaves/Z MakeGizmoRot_tri
//		color = GizmoObjectProperty("",item,"frontColor")// get color of the isoSurface
//		triList += Add_triangleSurface(item+"_tri",GetWavesDataFolder(triW,2),color)+";"
//	endfor
//
//	Execute "GetGizmo/Z displayNameList"					// find the index into display list of the iso-surfaces
//	SVAR S_DisplayNames=S_DisplayNames
//	Variable location = WhichListItem(StringFromList(ItemsInList(isoList)-1,isoList),S_DisplayNames)
//	KillStrings/Z S_DisplayNames
//	for (i=itemsInList(triList)-1;i>=0;i-=1)				// add the triangle waves to display list
//		Execute "ModifyGizmo insertDisplayList="+num2istr(location)+", object="+StringFromList(i,triList)
//	endfor
//	for (i=itemsInList(isoList)-1;i>=0;i-=1)				// remove isosurfaces from display list
//		Execute "RemoveFromGizmo/Z displayItem="+StringFromList(i,isoList)
//	endfor
//	Execute "ModifyGizmo endRecMacro"
//	return 0
//End
//
Static Function isMacBookPro()		// returns true for a MacBookPro
	if (!stringmatch(StringByKey("OS",IgorInfo(3)),"Macintosh OS X"))
		return 0
	endif
	String cmd
	sprintf cmd "do shell script \"sysctl hw.model\""
	ExecuteScriptText cmd						//returns something like: 	"hw.model:MacBookPro2,2"
	String model = StringByKey("hw.model",ReplaceString("\"",S_value,""))
	model = ReplaceString(" ",model,"")
	return stringmatch(model,"MacBookPro*")
End
//
Static Function/T Add_isoSurface(i,value,color,wName)// adds an iso-surface for MakeGizmoRot
	Variable i
	Variable value
	String color											// a color+alpha, e.g. "0.25,0.25,1,  0.5"	// light blue
	String wName

	String name="isoSurface"+num2istr(i)
	Execute "AppendToGizmo isoSurface="+wName+",name="+name
	Execute "ModifyGizmo ModifyObject="+name+" property={ surfaceColorType,1}"
	Execute "ModifyGizmo ModifyObject="+name+" property={ lineColorType,0}"
	Execute "ModifyGizmo ModifyObject="+name+" property={ lineWidthType,0}"
	Execute "ModifyGizmo ModifyObject="+name+" property={ fillMode,2}"
	Execute "ModifyGizmo ModifyObject="+name+" property={ lineWidth,1}"
	Execute "ModifyGizmo ModifyObject="+name+" property={ isoValue,"+num2str(value)+"}"
	Execute "ModifyGizmo ModifyObject="+name+" property={ frontColor,"+color+"}"
	Execute "ModifyGizmo ModifyObject="+name+" property={ backColor,"+color+"}"
	Execute "ModifyGizmo modifyObject="+name+" property={calcNormals,1}"
	return name
End
//
Static Function/T Add_triangleSurface(objName,wName,color)		// adds an triangle-surface for MakeGizmoRot
	String objName
	String wName
	String color										// color with alpha,  "0.25,0.25,1,  0.5"

	Execute "AppendToGizmo Surface="+wName+",name="+objName
	Execute "ModifyGizmo ModifyObject="+objName+" property={ surfaceColorType,1}"
	Execute "ModifyGizmo ModifyObject="+objName+" property={ srcMode,1}"
	Execute "ModifyGizmo ModifyObject="+objName+" property={ frontColor,"+color+"}"
	Execute "ModifyGizmo ModifyObject="+objName+" property={ backColor,"+color+"}"
	Execute "ModifyGizmo modifyObject="+objName+" property={calcNormals,1}"
	return objName
End
//
Function/T GizmoObjectProperty(gizmoName,objectName,propName)		// get value of a property for some gizmo object
	String gizmoName			// name of gizmo, use "" for top one
	String objectName			// name of object
	String propName			// name of property, e.g. frontColor or lineWidthType

	Execute "GetGizmo"+SelectString(strlen(gizmoName),"","/N="+gizmoName)+" objectList"
	String list = StrVarOrDefault("S_gizmoObjectList","")
	KillWaves/Z TW_gizmoObjectList
	KillStrings/Z S_gizmoObjectList
	Variable I,N=ItemsInLIst(list), j

	String item, match=" ModifyObject="+objectName+" ", matchProp=" property={ "+propName+","
	for (i=0;i<N;i+=1)
		item = StringFromList(i,list)
		j = strsearch(item,match,0)		// right object
		if (j<0)
			continue
		endif
		j = strsearch(item,matchProp,0)	// right property
		if (j<0)
			continue
		endif
		j += strlen(matchProp)
		item = item[j,Inf]
		j =  strsearch(item,"}",0)			// find and remove trailing "}"
		return item[0,j-1]
	endfor

	// did not find it yet, check the AppendToGizmo part
	list = ReplaceString("\t",list," ")
	match = ",name="+objectName+";"
	matchProp=" "+propName+"="

	for (i=0;i<N;i+=1)
		item = StringFromList(i,list)+";"
		j = strsearch(item,"AppendToGizmo",0)		// only search AppendToGizmo commands
		if (j<0 || j>5)
			continue
		endif
		j = strsearch(item,match,0)		// right object
		if (j<0)
			continue
		endif
		j = strsearch(item,matchProp,0)	// right property
		if (j<0)
			continue
		endif
		j += strlen(matchProp)
		item = item[j,Inf]
		j =  strsearch(item,",",0)			// find and remove trailing ","
		return item[0,j-1]
	endfor
	return ""
End
//Function/T GizmoObjectProperty(gizmoName,objectName,propName)		// get value of a property for some gizmo object
//	String gizmoName			// name of gizmo, use "" for top one
//	String objectName			// name of object
//	String propName			// name of property, e.g. frontColor or lineWidthType
//
//	Execute "GetGizmo"+SelectString(strlen(gizmoName),"","/N="+gizmoName)+" objectList"
//	String list = StrVarOrDefault("S_gizmoObjectList","")
//	KillWaves/Z TW_gizmoObjectList
//	KillStrings/Z S_gizmoObjectList
//	Variable I,N=ItemsInLIst(list), j
//
//	String item, match=" ModifyObject="+objectName+" ", match2=" property={ "+propName+","
//	for (i=0;i<N;i+=1)
//		item = StringFromList(i,list)
//		j = strsearch(item,match,0)		// right object
//		if (j<0)
//			continue
//		endif
//		j = strsearch(item,match2,0)		// right property
//		if (j<0)
//			continue
//		endif
//		j += strlen(match2)
//		item = item[j,Inf]
//		j =  strsearch(item,"}",0)			// find and remove trailing "}"
//		return item[0,j-1]
//	endfor
//	return ""
//End
//
//
//	setAxesInGizmo("Qx  (1/nm)","Qy  (1/nm)","Qz  (1/nm)")
Function setAxesInGizmo(xname,yname,zname)
	String xname,yname,zname
	if (ItemsInList(GetRTStackInfo(0))<2)
		Execute "ModifyGizmo startRecMacro"
	endif
	Execute "ModifyGizmo ModifyObject=axes0,property={-1,axisScalingMode,1}"
	Execute "ModifyGizmo ModifyObject=axes0,property={-1,axisColor,0,0,0,1}"
	Execute "ModifyGizmo ModifyObject=axes0,property={0,ticks,3}"
	Execute "ModifyGizmo ModifyObject=axes0,property={1,ticks,3}"
	Execute "ModifyGizmo ModifyObject=axes0,property={2,ticks,3}"
	Execute "ModifyGizmo ModifyObject=axes0,property={0,axisLabel,1}"
	Execute "ModifyGizmo ModifyObject=axes0,property={1,axisLabel,1}"
	Execute "ModifyGizmo ModifyObject=axes0,property={2,axisLabel,1}"
	Execute "ModifyGizmo ModifyObject=axes0,property={0,axisLabelText,\""+xname+"\"}"
	Execute "ModifyGizmo ModifyObject=axes0,property={1,axisLabelText,\""+yname+"\"}"
	Execute "ModifyGizmo ModifyObject=axes0,property={2,axisLabelText,,\""+zname+"\"}"
	Execute "ModifyGizmo ModifyObject=axes0,property={0,axisLabelCenter,-0.2}"
	Execute "ModifyGizmo ModifyObject=axes0,property={1,axisLabelCenter,-0.2}"
	Execute "ModifyGizmo ModifyObject=axes0,property={2,axisLabelCenter,-0.5}"
	Execute "ModifyGizmo ModifyObject=axes0,property={0,axisLabelDistance,0.05}"
	Execute "ModifyGizmo ModifyObject=axes0,property={1,axisLabelDistance,0.05}"
	Execute "ModifyGizmo ModifyObject=axes0,property={2,axisLabelDistance,0.3}"
	Execute "ModifyGizmo ModifyObject=axes0,property={0,axisLabelScale,0.5}"
	Execute "ModifyGizmo ModifyObject=axes0,property={1,axisLabelScale,0.5}"
	Execute "ModifyGizmo ModifyObject=axes0,property={2,axisLabelScale,0.5}"
	Execute "ModifyGizmo ModifyObject=axes0,property={0,axisLabelRGBA,0,0,0,1}"
	Execute "ModifyGizmo ModifyObject=axes0,property={1,axisLabelRGBA,0,0,0,1}"
	Execute "ModifyGizmo ModifyObject=axes0,property={2,axisLabelRGBA,0,0,0,1}"
	if (ItemsInList(GetRTStackInfo(0))<2)
		Execute "ModifyGizmo endRecMacro"
	endif
End



Static Function/T addFixedTitleStringToGizmo(pos,title)
	String pos									// position TL=top left, etc.
	String title									// String to display
	pos = SelectString(strlen(pos),"TL",pos)	// defaults to top left

	NewDataFolder/O root:Packages				// ensure Packages exists
	NewDataFolder/O root:Packages:JZT_Gizmo	// ensure geometry exists
	if (exists("root:Packages:JZT_Gizmo:stringNumber")!=2)
		Variable/G root:Packages:JZT_Gizmo:stringNumber=-1
	endif
	NVAR stringNumber = root:Packages:JZT_Gizmo:stringNumber
	stringNumber += 1
	String name="gizmoStringGroup"+num2istr(stringNumber)

	Execute "AppendToGizmo group,name="+name
	// ************************* Group Object Start *******************
	Execute "ModifyGizmo currentGroupObject=\""+name+"\""
	Execute "AppendToGizmo string=\""+title+"\",strFont=\"Geneva\",name=string0"
	Execute "AppendToGizmo attribute color={0,0,0,1},name=colorBlack"

	Execute "ModifyGizmo setDisplayList=0, opName=pushMatrix0, operation=pushMatrix"
	if (stringmatch(pos,"TL"))
		Execute "ModifyGizmo setDisplayList=1, opName=rotate0, operation=rotate, data={180,1,0,0}"
		Execute "ModifyGizmo setDisplayList=2, opName=translate0, operation=translate, data={-1.95,-1.9,0}"
	elseif (stringmatch(pos,"BL"))
		Execute "ModifyGizmo setDisplayList=1, opName=rotate0, operation=rotate, data={180,1,0,0}"
		Execute "ModifyGizmo setDisplayList=2, opName=translate0, operation=translate, data={-1.95,1.95,0}"
	endif

	Execute "ModifyGizmo setDisplayList=3, opName=scale0, operation=scale, data={0.1,0.1,0.1}"
	Execute "ModifyGizmo setDisplayList=4, attribute=colorBlack"
	Execute "ModifyGizmo setDisplayList=5, object=string0"
	Execute "ModifyGizmo setDisplayList=6, opName=loadIdentity0, operation=loadIdentity"
	Execute "ModifyGizmo setDisplayList=7, opName=popMatrix0, operation=popMatrix"
	Execute "ModifyGizmo currentGroupObject=\"::\""
	// ************************* Group Object End *******************
	return name

	//	Don't forget to put this before a 
	//	ModifyGizmo setDisplayList=-1, opName=MainTransform, operation=mainTransform
End
//
Static Function/T addStringToGizmo(pos,title)
	String pos						// position TL=top left, etc.
	String title						// String to display

	NewDataFolder/O root:Packages						// ensure Packages exists
	NewDataFolder/O root:Packages:JZT_Gizmo			// ensure geometry exists
	if (exists("root:Packages:JZT_Gizmo:stringNumber")!=2)
		Variable/G root:Packages:JZT_Gizmo:stringNumber=-1
	endif
	NVAR stringNumber = root:Packages:JZT_Gizmo:stringNumber
	stringNumber += 1
	String name="gizmoStringGroup"+num2istr(stringNumber)

	Execute "AppendToGizmo group,name="+name
	// ************************* Group Object Start *******************
	Execute "ModifyGizmo currentGroupObject=\""+name+"\""

	Execute "AppendToGizmo string=\""+title+"\",strFont=\"Geneva\",name=string0"
	Execute "AppendToGizmo attribute color={0,0,0,1},name=colorBlack"
	Execute "ModifyGizmo setDisplayList=0, opName=pushMatrix0, operation=pushMatrix"
	Execute "ModifyGizmo setDisplayList=1, opName=rotate0, operation=rotate, data={180,1,0,0}"
	Execute "ModifyGizmo setDisplayList=2, opName=translate0, operation=translate, data={0,0,0}"
	Execute "ModifyGizmo setDisplayList=3, opName=scale0, operation=scale, data={0.1,0.1,0.1}"
//	Execute "ModifyGizmo setDisplayList=-1, attribute=colorBlack"
	Execute "ModifyGizmo setDisplayList=-1, object=string0"
	Execute "ModifyGizmo setDisplayList=-1, opName=popMatrix0, operation=popMatrix"
	Execute "ModifyGizmo currentGroupObject=\"::\""
	// ************************* Group Object End *******************
	return name
End




Function Make_Graph_RX_RH_RF()
	Wave RX=RX, RH=RH, RF=RF
	Display /W=(264,53,844,301) RX,RH,RF
	ModifyGraph gfMult=130
	ModifyGraph rgb(RH)=(0,65535,0),rgb(RF)=(1,12815,52428)
	ModifyGraph tick=2, zero(left)=2, mirror=1, minor=1, lowTrip=0.001, standoff=0
	Legend/C/N=text0/J/A=RB/X=8.22/Y=11.50 "\\s(RX) RX\r\\s(RH) RH\r\\s(RF) RF"
EndMacro





//  ============================================================================  //
//  =============================== Start of Single Spot ===============================  //

Function DoAll(method)
	Variable method			//0=proper from peaks,  1=find peaks,  2=findpeaks+proper

	Wave singlePeakInfo=singlePeakInfo
	Wave/T imageNamesOriginal=imageNamesOriginal
	Variable peaksListFound = WaveExists(singlePeakInfo) && WaveExists(imageNamesOriginal)

	Variable valid = method==0 || method==1 || method==2
	valid = peaksListFound ? valid : (method>=1)&&valid

	if (!valid)
		method = 1
		if (peaksListFound)
			Prompt method, "calc method",popup,"Proper Only;Only Find Peaks;Find Peaks + Proper"
			DoPrompt "DoAlll",method
			method -=1
		else
			Prompt method, "calc method",popup,"Only Find Peaks;Find Peaks + Proper"
			DoPrompt "DoAlll",method
		endif
		if (V_flag)
			return -1
		endif
	endif

	Variable N=-1										// N is number of single spots images processed
	if (method==1 || method==2)						// for method 1 or 2, find the peaks
		printf "calling DoAllSingleSpotInfo()...\r"
		N = DoAllSingleSpotInfo(NaN,NaN)
	endif
	if (method==0 || method==2)						// for method not 1, compute rotations too
		printf "calling DoAllProper()...\r"
		N = DoAllProper()
	endif
	return N
End
//
Function DoAllSingleSpotInfo(width,threshold)	// calculates all positions of the single spot, and saves this in a waves singlePeakInfo & imageNamesOriginal
	Variable width									// width used in FitPeaksWithSeedFill()
	Variable threshold							// threshold above background for FitPeaksWithSeedFill(), and to find peak intensity

	Variable X0=2.8, H0=1.2							// values for wave note
X0 = 0
H0 = 4
	DoAlert 0, "Processing with X0="+num2str(X0)+",  H0="+num2str(H0)+",  you can change these later using ChangeAllX0H0()"
	if (!(width>0) || !(threshold>0))
		width = (width>0) ? width : 30
		threshold = (threshold>0) ? threshold : 20		// 20 is good for reconstructed images, for raw, use bigger value ~200
		Prompt width, "max peak width"
		Prompt threshold,"threshold above background for accepting peaks in FitPeaksWithSeedFill()"
		DoPrompt "peak parameters",width,threshold
		if (V_flag)
			return -1
		endif
	endif
	if (!(width>0) || !(threshold>0))
		String str
		sprintf str,"ERROR, width = %g,  threshold = %g,  both must be positive",width,threshold
		DoAlert 0, str
		return -1
	endif
//	String fileRoot = "WW_", rangePos="58211-62285", rangeDepths="12-60"
	String fileRoot="", rangePos="", rangeDepths=""
	String list = Get2RangesOfSPEfile()
	String path = StringByKey("path",list,"=")
	NewPath/O reconPath path
	PathInfo reconPath
	if (!V_flag)
		Abort "'reconPath' does not exist"
	endif
	fileRoot = StringByKey("name",list,"=")
	rangePos = StringByKey("range",list,"=")
	rangeDepths = StringByKey("depths",list,"=")
	if (numtype(str2num(rangePos)))
		return -1
	endif
	Prompt fileRoot,"first part of file name"
	Prompt rangePos,"range of positions"
	Prompt rangeDepths,"range of depths"
	DoPrompt "inputs",fileRoot,rangePos,rangeDepths
	if (V_flag)
		return -1
	endif
	printf "DoAllSingleSpotInfo(%g,%g), fileRoot='%s',  ranges='%s',  depths='%s'\r",width,threshold,fileRoot,rangePos,rangeDepths

	String progressWin = ProgressPanel("",stop=1)	// display a progress bar
	Variable Ndepths = (lastInRange(rangeDepths)-str2num(rangeDepths))		// number of depths to check (for progress)
	Variable Nmax = (lastInRange(rangePos)-str2num(rangePos)+1)*Ndepths	// guess at total number of images to check (for progress)
	Variable percent=0																// value of progress bar

	Make/N=100/O/T imageNamesOriginal
	Make/N=(100,8)/O/D singlePeakInfo				// to save X,H,F,px,py,totalIntensity,totalPeakIntensity,yc
	String fName
	Variable clip_top = NumVarOrDefault("root:CLIP_TOP",0)
	clip_top = (clip_top>0 && numtype(clip_top)==0) ? clip_top : 0
	Variable X1,Y1,Z1, depthSi, w2=width/2
	Variable idepth,ipos, N=0							// number of images found
	Variable px,py, idim
	Variable groupx,groupy,startx,starty, CCDy
	Variable timer=startMSTimer,timer1percent=startMSTimer
	for (ipos=str2num(rangePos),N=0; !numtype(ipos); ipos=NextInRange(rangePos,ipos))
		sprintf fName, "%s%s%d_%d.SPE",path,fileRoot,ipos,str2num(rangeDepths)
		GetFileFolderInfo/Q/Z fName						// check if ipos is valid
		if (V_Flag || !V_isFile)
			continue										// this ipos not valid, skip it
		endif
		for (idepth=str2num(rangeDepths); !numtype(idepth); idepth=NextInRange(rangeDepths,idepth))
			percent = (ipos-str2num(rangePos))*Ndepths / Nmax * 100	// update progress bar
			ValDisplay progressValdisp,win=$progressWin, value= #num2istr(floor(limit(percent,0,100)))
#if NumberByKey("IGORVERS",IgorInfo(0))>=6.10
			DoUpdate/W=$progressWin					// update only progress window
			if( V_Flag == 2 )								// stop button on progress window was pushed
				ipos = Inf									// this line causes outer loop to also end
				break										//   and break out of inner loop
			endif
#else
			DoUpdate										// update all windows, no stop button
#endif
			if (timer1percent>=0 && percent>1)		// estimate execution time based on first 1%
				Variable remaining = stopMSTimer(timer1percent)*1e-6 * (100-percent)
				timer1percent = -1
				if (remaining > 10)						// tell user if more than 10 sec
					printf "This process should finish in another %s sec,   which will be at %s\r",Secs2Time(remaining,5,1),Secs2Time(DateTime+remaining,1)
				endif
			endif

			sprintf fName, "%s%s%d_%d.SPE",path,fileRoot,ipos,idepth
			Wave image = $(WinViewReadROI(fName,0,-1,0,-1))	// load image
			if (!WaveExists(image))
				continue
			endif

			if (DimSize(singlePeakInfo,0)<=N)			// extend waves
				idim = DimSize(singlePeakInfo,0)+100
				Redimension/N=(idim), imageNamesOriginal
				Redimension/N=(idim,8) singlePeakInfo// to save X,H,F,px,py,totalIntensity,totalPeakIntensity,yc
			endif
			imageNamesOriginal[N] = NameOfWave(image)
			singlePeakInfo[N][] = NaN

			X1 = WinViewInfo(image,"X1")				// things that we always know
			Y1 = WinViewInfo(image,"Y1")
			Z1 = WinViewInfo(image,"Z1")
			//	allow user to optionally override using a Keyence correction, default is to use it.
			if (NumVarOrDefault("UseKeyenceCorrection",1 ))	// look for optional global flag
				X1 = X1correctedKeyence(X1)			// takes PM500 X1 and returns the "real" X1
				Y1 = Y1correctedKeyence(Y1)			// takes PM500 Y1 and returns the "real" Y1
				Z1 = Z1correctedKeyence(Z1)			// takes PM500 Z1 and returns the "real" Z1
			endif

			singlePeakInfo[N][0] = -X1					// we want coordinates in the sample, not coords of the PM500
			singlePeakInfo[N][1] = -YZ2H(Y1,Z1)		// H
			singlePeakInfo[N][2] = -YZ2F(Y1,Z1)		// F

			depthSi = WinViewInfo(image,"depthSi")
			singlePeakInfo[N][1] +=depthSi/sqrt(2)	// correct H and Z for the depth
			singlePeakInfo[N][2] +=depthSi/sqrt(2)
			singlePeakInfo[N][5] = sum(image)
			singlePeakInfo[N][7] = WinViewInfo(image,"yc")


//if(NumVarOrDefault("root:CLIP_TOP",0)>0)
//	clip_top = root:CLIP_TOP
//image[][0,50]=0
//endif
if(clip_top)
image[][0,clip_top]=0
endif
//			Wave FullPeakList=$FitPeaksWithSeedFill(image,1.2,70,50,20,maxNu=1)		// include a box ±width/2 around the peak
			Wave FullPeakList=$FitPeaksWithSeedFill(image,1.2,2*width,2*width,threshold,maxNu=1)		// include a box ±width around the peak
			if (WaveExists(FullPeakList) && DimSize(FullPeakList,0)>0)
				px = FullPeakList[0][0]
				py = FullPeakList[0][1]
				singlePeakInfo[N][3] = px
				singlePeakInfo[N][4] = py
				ImageStats/M=1/G={px-w2,px+w2,py-w2,py+w2} image
				singlePeakInfo[N][6] = V_avg*V_npnts	// totalPeakIntensity
				groupx = WinViewInfo(image,"groupx")	// really only need to do set these once
				groupy = WinViewInfo(image,"groupy")
				startx = WinViewInfo(image,"startx")
				starty = WinViewInfo(image,"starty")
				CCDy = WinViewInfo(image,"CCDy")
			endif
			KillWaves/Z image
			N += 1
		endfor
	endfor
	Redimension/N=(N), imageNamesOriginal
	Redimension/N=(N,8) singlePeakInfo

	String noteStr = ""
	noteStr = ReplaceStringByKey("fldrName",noteStr,GetDataFolder(1),"=")
	noteStr = ReplaceNumberByKey("dX",noteStr,0,"=")
	noteStr = ReplaceNumberByKey("dH",noteStr,0,"=")
	noteStr = ReplaceNumberByKey("dF",noteStr,1,"=")
	noteStr = ReplaceNumberByKey("X0",noteStr,X0,"=")
	noteStr = ReplaceNumberByKey("H0",noteStr,H0,"=")
	noteStr = ReplaceNumberByKey("F0",noteStr,0,"=")
	noteStr = ReplaceStringByKey("waveClass",noteStr,"Random3dSingleSpotInfo","=")
	noteStr = ReplaceNumberByKey("groupx",noteStr,groupx,"=")
	noteStr = ReplaceNumberByKey("groupy",noteStr,groupy,"=")
	noteStr = ReplaceNumberByKey("startx",noteStr,startx,"=")
	noteStr = ReplaceNumberByKey("starty",noteStr,starty,"=")
	noteStr = ReplaceNumberByKey("CCDy",noteStr,CCDy,"=")
	Note/K imageNamesOriginal, noteStr
	Note/K singlePeakInfo, noteStr
	SetDimLabel 1,0,XX,singlePeakInfo			;	SetDimLabel 1,1,HH,singlePeakInfo		;		SetDimLabel 1,2,FF,singlePeakInfo
	SetDimLabel 1,3,px,singlePeakInfo			;	SetDimLabel 1,4,py,singlePeakInfo
	SetDimLabel 1,5,totalIntensity,singlePeakInfo;	SetDimLabel 1,6,totalPeakIntensity,singlePeakInfo
	SetDimLabel 1,7,yc,singlePeakInfo
	beep

	DoWindow/K $progressWin
	Variable elapsed = stopMSTimer(timer)*1e-6
	printf "the whole peak fitting process took %s\r",Secs2Time(elapsed,5,0)
	if (timer1percent>=0 )						// ensure this timer was stopped
		timer1percent = stopMSTimer(timer1percent)
	endif
	if (elapsed>10*60)									// if execution took more than 10 min, automatically save
		print "This took more than 10min, so save the experiment"
		SaveExperiment
	endif
	return N
End
//
//
Static Function DoAllProper()				// calculates rotation angles using proper pixel2qhat routines from precalculated peak positions in singlePeakInfo
	Variable pkX0 = NumVarOrDefault("root:Packages:micro:Arrays3d:pxPrimary",NaN)	// center of an un-tilted spot, uses binned pixels, not re-set to zero either
	Variable pkY0 = NumVarOrDefault("root:Packages:micro:Arrays3d:pyPrimary",NaN)	// used to define the reference (un-tilted) orientation
	Variable yc0 = NumVarOrDefault("root:Packages:micro:Arrays3d:ycPrimary",NaN)
	if (numtype(pkX0+pkY0+yc0))
		DoAlert 0, "bad values for reference peak position"
		printf "pkX0=%g,  pkY0=%g,  yc0=%g\r",pkX0,pkY0,yc0
		return 1
	endif
	Variable pkX1 = NumVarOrDefault("root:Packages:micro:Arrays3d:pxSecondary",NaN)	// center of an un-tilted spot, uses binned pixels, not re-set to zero either
	Variable pkY1 = NumVarOrDefault("root:Packages:micro:Arrays3d:pySecondary",NaN)	// used to define the secondary reference (un-tilted) orientation

	Wave singlePeakInfo=singlePeakInfo
	Wave/T imageNamesOriginal=imageNamesOriginal
	if (!WaveExists(singlePeakInfo) || !WaveExists(imageNamesOriginal))
		DoAlert 0, "singlePeakInfo or imageNamesOriginal does not exist, find peaks first"
		return -1
	endif
	String noteStr = note(singlePeakInfo)
	if (!stringmatch(StringByKey("waveClass",noteStr,"="),"Random3dSingleSpotInfo"))
		DoAlert 0, "singlePeakInfo is of wrong waveClass"
		return -1
	endif
	Variable N=DimSize(singlePeakInfo,0)
	Variable doubleSpot = (DimSize(singlePeakInfo,1)>=11)
	Wave hkl0=root:Packages:micro:Arrays3d:hkl0, hkl1=root:Packages:micro:Arrays3d:hkl1
	if (doubleSpot && (!WaveExists(hkl0) || !WaveExists(hkl1)))
		DoAlert 0,"Double Spot data, cannot find hkl0 or hkl1"
		return 1
	endif

	STRUCT microGeometry geo
	if (FillGeometryStructDefault(geo))							//fill the geometry structure with test values
		DoAlert 0, "no geometry structure found, did you forget to set it?"
		return -1
	endif
	Variable timer=startMSTimer
	Make/N=(N)/O XX=NaN,HH=NaN,FF=NaN, RX=NaN,RH=NaN,RF=NaN, totalPeakIntensity=NaN,totalIntensity=NaN, totalAngles=NaN
	Make/N=(N)/O/T imageNames
	SetScale d 0,0,"µm", XX,HH,FF
	Variable i
	Variable px,py, angle, ddLocal, ycLocal, pxb,pyb
	Variable groupx,groupy, startx, starty, CCDy
	groupx = NumberByKey("groupx",noteStr,"=")
	groupy = NumberByKey("groupy",noteStr,"=")
	startx = NumberByKey("startx",noteStr,"=")
	starty = NumberByKey("starty",noteStr,"=")
	CCDy = NumberByKey("CCDy",noteStr,"=")

	ddLocal = geo.ddOffset + CCDy										// get dd from the wave note of singlePeakInfo
	geo.dd = ddLocal>0 ? ddLocal :  geo.dd
	geo.ycent = (yc0>0) ? yc0 : geo.ycent

	XX = singlePeakInfo[p][0]
	HH = singlePeakInfo[p][1]
	FF = singlePeakInfo[p][2]
	totalIntensity = singlePeakInfo[p][5]
	totalPeakIntensity = singlePeakInfo[p][6]
	imageNames = imageNamesOriginal[p]

	// Reference Orientation
	Make/N=3/O/D DoAllProper_Rodriques, DoAllProper_qhat, DoAllProper_qhat0
	Wave Rodriques=DoAllProper_Rodriques, qhat=DoAllProper_qhat, qhat0=DoAllProper_qhat0
	pixel2q(geo,pkX0*groupx+startx,pkY0*groupy+starty,qhat0)// direction on primary reference spot (only used for single spot analysis)
	wenge2IdealBeamline(geo,qhat0,qhat0)							// qhat0 is used for single spot analysis
	if (doubleSpot)
		Make/N=3/O/D DoAllProper_Qhkl0, DoAllProper_xyz0, DoAllProper_Qhkl1, DoAllProper_xyz1
		Wave Qhkl0=DoAllProper_Qhkl0, xyz0=DoAllProper_xyz0, Qhkl1=DoAllProper_Qhkl1, xyz1=DoAllProper_xyz1
		Make/N=(3,3)/O/D DoAllProper_rhoRef, DoAllProper_rhoMeas
		Wave rhoRef=DoAllProper_rhoRef, rhoMeas=DoAllProper_rhoMeas
		xyz0 = qhat0
		pixel2q(geo,pkX1*groupx+startx,pkY1*groupy+starty,xyz1)
		wenge2IdealBeamline(geo,xyz1,xyz1)						// measured direction of secondary spot
		hkl2Q(hkl0[0],hkl0[1],hkl0[2], Qhkl0,normal=1)		// get Q of reference spots in std orientation Qhkl
		hkl2Q(hkl1[0],hkl1[1],hkl1[2], Qhkl1,normal=1)
		Rfrom2spots(Qhkl0,xyz0,Qhkl1,xyz1,rhoRef)				// compute rhoRef, rotation that takes recp lattice hkl0 & hkl1 to the reference orientation
	endif

	for (i=0;i<N;i+=1)
		px = singlePeakInfo[i][3]									// primary spot
		py = singlePeakInfo[i][4]
		geo.ycent = singlePeakInfo[i][7]
		if (doubleSpot)
			pxb = singlePeakInfo[i][8]								// secondary spot
			pyb = singlePeakInfo[i][9]
		else
			pxb=NaN;	pyb=NaN
		endif
		if (numtype(px+py+pxb+pyb)==0)						// a valid two spot image
			pixel2q(geo,px*groupx+startx,py*groupy+starty,xyz0)// use 0 based un-binned pixels
			wenge2IdealBeamline(geo,xyz0,xyz0)
			pixel2q(geo,pxb*groupx+startx,pyb*groupy+starty,xyz1)// use 0 based un-binned pixels
			wenge2IdealBeamline(geo,xyz1,xyz1)

			Rfrom2spots(Qhkl0,xyz0,Qhkl1,xyz1,rhoMeas)		// compute rhoRef, rotation that takes standard orientation to the reference orientation

			MatrixOp/O rhoMeas = rhoMeas x Inv(rhoRef)
			angle = axisOfMatrix(rhoMeas,Rodriques,squareUp=1)		// construct the Rodriques vector
			totalAngles[i] = angle
			angle = tan(angle/2 * PI/180)							// normalize Rodriques so that length is tan(theta/2)
			Rodriques *= angle
			RX[i] = Rodriques[0]
			RH[i] = -YZ2H(Rodriques[1], Rodriques[2])			// H =  Y*sin(angle) + Z*cos(angle).	This is really sample Y, not H
			RF[i] = -YZ2F(Rodriques[1], Rodriques[2])			// F = -Y*cos(angle) + Z*sin(angle).	This is really sample Z, not F
		elseif (numtype(px+py)==0 && !doubleSpot)				// a valid single spot (only expecting one spot)
			pixel2q(geo,px*groupx+startx,py*groupy+starty,qhat)// use 0 based un-binned pixels
			wenge2IdealBeamline(geo,qhat,qhat)
			Cross qhat0, qhat										// qhat0 x qhat is almost the Rodriques vector (just the wrong length)
			Wave W_Cross=W_Cross
			Rodriques = W_Cross
			angle = asin(normalize(Rodriques))						// length is sin(theta), want tan(theta/2) for Rodriques vector
			totalAngles[i] = angle * 180/PI
			Rodriques *= tan(angle/2)
			RX[i] = Rodriques[0]
			RH[i] = -YZ2H(Rodriques[1], Rodriques[2])			// H =  Y*sin(angle) + Z*cos(angle)
			RF[i] = -YZ2F(Rodriques[1], Rodriques[2])			// F = -Y*cos(angle) + Z*sin(angle)
		endif
	endfor
	KillWaves/Z DoAllProper_qhat, DoAllProper_qhat0, DoAllProper_Rodriques, W_Cross
	KillWaves/Z DoAllProper_Qhkl0, DoAllProper_xyz0, DoAllProper_Qhkl1, DoAllProper_xyz1, DoAllProper_rhoRef, DoAllProper_rhoMeas
	KillWaves/Z Rfrom2spots_rho0, Rfrom2spots_rho1, Rfrom2spots_h1,Rfrom2spots_x1, Rfrom2spots_axis

	// center coordinates at the point {X0,H0,F0}.  This resets the origin to {0,0,0}
	WaveStats/Q/M=1 XX
	XX -= V_avg
	WaveStats/Q/M=1 HH
	HH -= V_avg
	WaveStats/Q/M=1 FF
	FF -= V_avg

	noteStr = ReplaceStringByKey("waveClass",noteStr,"Random3dArrays","=")
	Note/K XX, noteStr ;			Note/K HH, noteStr ;	Note/K FF, noteStr
	Note/K RX, noteStr ;			Note/K RH, noteStr ;	Note/K RF, noteStr
	Note/K totalIntensity, noteStr ;						Note/K totalPeakIntensity, noteStr
	Note/K totalAngles, noteStr ;							Note/K imageNames, noteStr

	beep
//	TrimBadPointsOnOutside()
//	Interpolate3dRodriquesClosest(1)

	beep
	Variable elapsed = stopMSTimer(timer)*1e-6
	if (elapsed>3)
		printf "this process took %s\r",Secs2Time(elapsed,5,0)
	endif
	if (elapsed>10*60)									// if execution took more than 10 min, automatically save
		print "This took more than 10min, so saving the experiment"
		SaveExperiment
	endif
	return N
End
//	compute rotation matrix that will take unrotated reciprocal lattice to alignment with measured xyz0, xyz1
Static Function Rfrom2spots(Qhkl0,xyz0,Qhkl1,xyz1,rmat)
	Wave Qhkl0						// direction of reciprocal lattice vector for hkl of main spot (normalized!!!), direction of from the standard recip-lattice
	Wave xyz0						// unit vector pointing to main spot (obtained from pixel position) (NOT normalized)
	Wave Qhkl1						// reciprocal lattice vector for hkl of secondary spot (normalized!!!), this is actually Q, not hkl
	Wave xyz1						// unit vector pointing to secondary spot (obtained from pixel position) (NOT normalized)
	Wave rmat						// resultant rotation matrix

	Make/N=(3,3)/O/D Rfrom2spots_rho0, Rfrom2spots_rho1
	Wave rho0=Rfrom2spots_rho0						// first rotation, get Qhkl0 aligned with xyz0
	Wave rho1=Rfrom2spots_rho1						// second rotation, aligns Qhkl1 with xyz1
	Make/N=3/O/D Rfrom2spots_h1,Rfrom2spots_x1, Rfrom2spots_axis	// second spot in rotated space (by rho0)
	Wave h1=Rfrom2spots_h1, x1=Rfrom2spots_x1, axis=Rfrom2spots_axis
	Variable dot, angle, xx,yy

	// find a rotation that will align Qhkl0 to xyz0 (assumes that Qhkl0 is normalized)
	Cross Qhkl0, xyz0									// qhat0 x qhat is almost the Rotation vector to align xyz0 and Qhkl0, but it has the wrong length
	Wave W_Cross=W_Cross
	yy = normalize(W_Cross)							// length of cross is sin(theta)
	xx = MatrixDot(Qhkl0,xyz0)						// cos(theta)
	angle = atan2(yy,xx)*180/PI						// total rotation angle
	rotationMatAboutAxis(W_Cross,angle,rho0)			// calculate rotation matrix from Qhkl0 to xyz0

	// rotate about xyz0 to get Qhkl1 and xyz1 into best alignment (assumes that Qhkl1 is normalized)
	axis = xyz0											// second rotation is about direction of primary spot (i.e. "axis")
	normalize(axis)

	dot=MatrixDot(xyz1,axis)
	x1 = xyz1 - dot*axis								// component of xyz1 perpendicular to axis (==xyz0)

	MatrixOp/O h1 = rho0 x Qhkl1						// rotate Qhkl1 into orientation where (xyz0 || Qhkl0)
	dot=MatrixDot(h1,axis)
	h1 -= dot*axis										// component of Qhkl1 perpendicular to axis

	dot = MatrixDot(x1,h1)/(norm(x1)*norm(h1))
	angle = acos(limit(dot,-1,1))*180/PI				// rotate around axis by this angle (but don't know sign yet)
	Cross h1, x1										// used to find sign for angle
	dot = MatrixDot(W_Cross,axis)						// if dot==0, then problem is unsolvable (two spots are colinear)
	angle = dot>0 ? angle : -angle						// correctly signed angle
	rotationMatAboutAxis(axis,angle,rho1)				// calculate rotation matrix from Qhkl1 to xyz1

	MatrixOp/O rmat = rho1 x rho0						// rmat is now the total rotation matrix that best aligns both Qhkl0 and Qhkl1

	// finally do the 180 rotation about X (not sure why this is needed, but it is)
	rho0 = -(p==q)
	rho0[0][0] = 1
	MatrixOp/O rmat = rho0 x rmat						// rmat is now the total rotation matrix that best aligns both Qhkl0 and Qhkl1
End
//
//Function test_Rfrom2spots2()
//	STRUCT microGeometry geo
//	if (FillGeometryStructDefault(geo))							//fill the geometry structure with test values
//		DoAlert 0, "no geometry structure found, did you forget to set it?"
//		return -1
//	endif
//
//	Variable pkX0 = NumVarOrDefault("root:Packages:micro:Arrays3d:pxPrimary",NaN)	// center of an un-tilted spot, uses binned pixels, not re-set to zero either
//	Variable pkY0 = NumVarOrDefault("root:Packages:micro:Arrays3d:pyPrimary",NaN)	// used to define the reference (un-tilted) orientation
//	Variable yc0 = NumVarOrDefault("root:Packages:micro:Arrays3d:ycPrimary",NaN)
//	//pkX0=1824.91;	pkY0=789.54;	pkX1=410.18;	pkY1=111.93
//
//	if (numtype(pkX0+pkY0+yc0))
//		DoAlert 0, "bad values for reference peak position"
//		printf "pkX0=%g,  pkY0=%g,  yc0=%g\r",pkX0,pkY0,yc0
//		return 1
//	endif
//	Variable pkX1 = NumVarOrDefault("root:Packages:micro:Arrays3d:pxSecondary",NaN)	// center of an un-tilted spot, uses binned pixels, not re-set to zero either
//	Variable pkY1 = NumVarOrDefault("root:Packages:micro:Arrays3d:pySecondary",NaN)	// used to define the secondary reference (un-tilted) orientation
//	Wave hkl0=root:Packages:micro:Arrays3d:hkl0, hkl1=root:Packages:micro:Arrays3d:hkl1
//
//	Wave singlePeakInfo=singlePeakInfo
//	String noteStr=note(singlePeakInfo)
//
//	Variable groupx = NumberByKey("groupx",noteStr,"=")
//	Variable groupy = NumberByKey("groupy",noteStr,"=")
//	Variable startx = NumberByKey("startx",noteStr,"=")
//	Variable starty = NumberByKey("starty",noteStr,"=")
//	Variable CCDy = NumberByKey("CCDy",noteStr,"=")
//	Variable ddLocal = geo.ddOffset + CCDy							// get dd from the wave note of singlePeakInfo
//	geo.dd = ddLocal>0 ? ddLocal :  geo.dd
//	geo.ycent = (yc0>0) ? yc0 : geo.ycent
//
//
//	// Two Spot Orientation
//	print "from two point fit:"
//	Make/N=3/O/D DoAllProper_Qhkl0, DoAllProper_xyz0, DoAllProper_Qhkl1, DoAllProper_xyz1
//	Wave Qhkl0=DoAllProper_Qhkl0, xyz0=DoAllProper_xyz0, Qhkl1=DoAllProper_Qhkl1, xyz1=DoAllProper_xyz1
//	Make/N=(3,3)/O/D DoAllProper_rhoTwo
//	Wave rhoTwo=DoAllProper_rhoTwo
//
//	pixel2q(geo,pkX0*groupx+startx,pkY0*groupy+starty,xyz0)	// direction on primary reference spot (only used for single spot analysis)
//	wenge2IdealBeamline(geo,xyz0,xyz0)							// xyz0 is used for single spot analysis
//	printf "pixel(%g, %g) --> %s ideal beam line (%s)\r",pkX0,pkY0,vec2str(xyz0),hkl2str(hkl0[0],hkl0[1],hkl0[2])
//
//	pixel2q(geo,pkX1*groupx+startx,pkY1*groupy+starty,xyz1)
//	wenge2IdealBeamline(geo,xyz1,xyz1)							// measured direction of secondary spot
//	printf "pixel(%g, %g) --> %s ideal beam line (%s)\r",pkX1,pkY1,vec2str(xyz1),hkl2str(hkl1[0],hkl1[1],hkl1[2])
//	hkl2Q(hkl0[0],hkl0[1],hkl0[2], Qhkl0,normal=1)				// get Q of reference spots in std orientation Qhkl
//	hkl2Q(hkl1[0],hkl1[1],hkl1[2], Qhkl1,normal=1)
//	printf "hkl(%s) --> Qhkl0 = {%s}	reference orientation\r",hkl2str(hkl0[0],hkl0[1],hkl0[2]),hkl2str(Qhkl0[0],Qhkl0[1],Qhkl0[2])
//	printf "hkl(%s) --> Qhkl1 = {%s}\r",hkl2str(hkl1[0],hkl1[1],hkl1[2]),hkl2str(Qhkl1[0],Qhkl1[1],Qhkl1[2])
//	printf "Qhkl0 ^ Qhkl1 = %g¡\r",angleBetweenVectors(Qhkl0,Qhkl1)
//	printf "xyz0 ^ xyz1 = %g¡\r",angleBetweenVectors(xyz0,xyz1)
//	Rfrom2spots(Qhkl0,xyz0,Qhkl1,xyz1,rhoTwo)					// compute rhoTwo, rotation that takes reference recip to measured hkl0 & hkl1
////	printWave(rhoTwo)
//
//	Make/N=3/O/D axisTwo,axisIndex
//	Variable angle = axisOfMatrix(rhoTwo,axisTwo)
//	printf "axisTwo = %s,  angle(rhoTwo) = %g¡\r",vec2str(axisTwo),angle
//	axisTwo *= tan(angle/2*PI/180)								// convert to Rodriques vector
//
//	// Indexing result (using full pattern and Euler)
//	Make/N=(3,3)/O/D rLattice,rLattice0,rhoIndex=nan
//	rLattice[0][0] = 17.3617088;		rLattice[0][1] = 0.2654619;		rLattice[0][2] = -0.7758887	// from generic_Index.txt
//	rLattice[1][0] = 0.7165537;		rLattice[1][1] = -12.9080763;	rLattice[1][2] = 11.6176383	// indexing result that has correct orientation
//	rLattice[2][0] = -0.3987782;		rLattice[2][1] = -11.6366877;	rLattice[2][2] = -12.9046457
//	Variable len = MatrixDet(rLattice)^(1/3)
//	rLattice /= len
//	printf "\r  from full frame indexing:\r"
////	printWave(rLattice)
//
//	STRUCT crystalStructure xtal
//	if (FillCrystalStructDefault(xtal))								//fill the lattice structure with current values
//		DoAlert 0, "No Lattice, please set one"
//		return 1
//	endif
//	rLattice0[0][0] = xtal.as0;		rLattice0[0][1] = xtal.bs0;		rLattice0[0][2] = xtal.cs0		// the reference reciprocal lattice (un-rotated)
//	rLattice0[1][0] = xtal.as1;		rLattice0[1][1] = xtal.bs1;		rLattice0[1][2] = xtal.cs1
//	rLattice0[2][0] = xtal.as2;		rLattice0[2][1] = xtal.bs2;		rLattice0[2][2] = xtal.cs2
//	len = MatrixDet(rLattice0)^(1/3)
//	rLattice0 /= len
//
//	MatrixOp/O rhoIndex = rLattice x Inv(rLattice0)
////	printWave(rhoIndex)
//	angle = axisOfMatrix(rhoIndex,axisIndex)
//	printf "axisIndex = %s,  angle(rhoIndex) = %g¡\r",vec2str(axisIndex),angle
//	printf "axisIndex ^ axisTwo = %g¡\r",angleBetweenVectors(axisIndex,axisTwo)
//	axisIndex *= tan(angle/2*PI/180)								// convert to Rodriques vector
//	printf "total angle between axisTwo and axisIndex = %g¡\r",AngleBetweenRotationVectors(axisTwo,axisIndex)
//	KillWaves/Z mat2_JZT,mat1_JZT,temp_angle_11_JZT
//	Make/N=(3,3)/O/D MatrixTemp
//	MatrixTemp = rhoTwo - rhoIndex
//	printf "max difference (rhoTwo-rhoIndex)= %g\r",max(abs( WaveMax(MatrixTemp)),abs(WaveMin(MatrixTemp)))
//	KillWaves/Z MatrixTemp
//
//	// more checking
//	printf "\r  more testing:\r"
//	Make/N=3/O/D xyz
//	Variable px=1825, py=789
//	pixel2XYZ(geo,px,py,xyz)
//	wenge2IdealBeamline(geo,xyz,xyz)
//	xyz /= 1000
//	printf "(%g, %g) --> %s\r",px,py,vec2str(xyz)
//
//	printf "rho x {001} shold be near {0,1,-1} in beam line coords (rotate about X,  approx. 135¡)\r"
//	xyz = {0,0,1}
//	MatrixOp/O axisTwo = rhoTwo x xyz
//	axisTwo *= sqrt(2)
//	printf "rhoTwo x {%s} = {%s}\r",hkl2str(xyz[0],xyz[1],xyz[2]),hkl2str(axisTwo[0],axisTwo[1],axisTwo[2])
//
//	xyz = {2,0,8}
//	MatrixOp/O axisTwo = rhoTwo x xyz
//	normalize(axisTwo)
//	printf "rhoTwo x {%s} = {%s}\r",hkl2str(xyz[0],xyz[1],xyz[2]),hkl2str(axisTwo[0],axisTwo[1],axisTwo[2])
//
//	xyz = {-1,1,9}
//	MatrixOp/O axisTwo = rhoTwo x xyz
//	normalize(axisTwo)
//	printf "rhoTwo x {%s} = {%s}\r",hkl2str(xyz[0],xyz[1],xyz[2]),hkl2str(axisTwo[0],axisTwo[1],axisTwo[2])
//
//
//	xyz = {0,0,1}
//	MatrixOp/O axisIndex = rhoIndex x xyz
//	axisIndex *= sqrt(2)
//	printf "rhoIndex x {%s} = {%s}\r",hkl2str(xyz[0],xyz[1],xyz[2]),hkl2str(axisIndex[0],axisIndex[1],axisIndex[2])
//
//	xyz = {2,0,8}
//	MatrixOp/O axisIndex = rhoIndex x xyz
//	normalize(axisIndex)
//	printf "rhoIndex x {%s} = {%s}\r",hkl2str(xyz[0],xyz[1],xyz[2]),hkl2str(axisIndex[0],axisIndex[1],axisIndex[2])
//
//	xyz = {-1,1,9}
//	MatrixOp/O axisIndex = rhoIndex x xyz
//	normalize(axisIndex)
//	printf "rhoIndex x {%s} = {%s}\r",hkl2str(xyz[0],xyz[1],xyz[2]),hkl2str(axisIndex[0],axisIndex[1],axisIndex[2])
//
//
//	// check the orientations and coordinate systems
//	Make/N=3/O/D xyz0Temp
//	pixel2q(geo,0,937,xyz0Temp)	// spot on left edge of detector
//	printf "\r  pixel(%g, %g) --> %s ideal beam line (wenge)\r",pkX0,pkY0,vec2str(xyz0Temp)
//	wenge2IdealBeamline(geo,xyz0Temp,xyz0Temp)
//	printf "pixel(%g, %g) --> %s ideal beam line (ideal beam line)\r",pkX0,pkY0,vec2str(xyz0Temp)
//	print "I remember that low xPixel is out the door (positive X)"
//	KillWaves/Z xyz0Temp
//
//
//	// check orientations by putting in known offsets and observing the kinds of rotations
//	print " "
//	printf "pixel(%g, %g) --> %s ideal beam line (%s)\r",pkX0,pkY0,vec2str(xyz0),hkl2str(hkl0[0],hkl0[1],hkl0[2])
//	printf "pixel(%g, %g) --> %s ideal beam line (%s)\r",pkX1,pkY1,vec2str(xyz1),hkl2str(hkl1[0],hkl1[1],hkl1[2])
//
//	Make/N=(3,3)/O/D DoAllProper_rhoRef
//	Wave rhoRef=DoAllProper_rhoRef
//	rhoRef = rhoTwo
//
//	print "\r  both dx=0 & dy = 0 pixels, should be 0 rotation"
//	Variable dx=0,dy=0
//	pixel2q(geo,(pkX0+dx)*groupx+startx,(pkY0+dy)*groupy+starty,xyz0)	// direction on primary reference spot (only used for single spot analysis)
//	wenge2IdealBeamline(geo,xyz0,xyz0)							// xyz0 is used for single spot analysis
//	printf "pixel(%g, %g) --> %s ideal beam line (%s)\r",(pkX0+dx),(pkY0+dy),vec2str(xyz0),hkl2str(hkl0[0],hkl0[1],hkl0[2])
//	pixel2q(geo,(pkX1+dx)*groupx+startx,(pkY1+dy)*groupy+starty,xyz1)
//	wenge2IdealBeamline(geo,xyz1,xyz1)							// measured direction of secondary spot
//	printf "pixel(%g, %g) --> %s ideal beam line (%s)\r",(pkX1+dx),(pkY1+dy),vec2str(xyz1),hkl2str(hkl1[0],hkl1[1],hkl1[2])
//	Rfrom2spots(Qhkl0,xyz0,Qhkl1,xyz1,rhoTwo)					// compute rhoTwo, rotation that takes reference recip to measured hkl0 & hkl1
//	Make/N=(3,3)/O/D Rfrom2spots_rho0
//	Wave rho0=Rfrom2spots_rho0
//	MatrixOp/O rho0 = rhoTwo x Inv(rhoRef)
//	rho0 = abs(rho0[p][q])<1e-12 ? 0 : rho0[p][q]
//	// printWave(rho0)
//	angle = axisOfMatrix(rho0,axisTwo)
//	printf "axis(rho0) = %s,  angle(rho0) = %g¡\r",vec2str(axisTwo),angle
//
//
//	print "\r  both dx=0 & dy = 10 pixels, should be a positive X rotation"
//	dx=0;	dy=10
//	pixel2q(geo,(pkX0+dx)*groupx+startx,(pkY0+dy)*groupy+starty,xyz0)	// direction on primary reference spot (only used for single spot analysis)
//	wenge2IdealBeamline(geo,xyz0,xyz0)							// xyz0 is used for single spot analysis
//	printf "pixel(%g, %g) --> %s ideal beam line (%s)\r",(pkX0+dx),(pkY0+dy),vec2str(xyz0),hkl2str(hkl0[0],hkl0[1],hkl0[2])
//	pixel2q(geo,(pkX1+dx)*groupx+startx,(pkY1+dy)*groupy+starty,xyz1)
//	wenge2IdealBeamline(geo,xyz1,xyz1)							// measured direction of secondary spot
//	printf "pixel(%g, %g) --> %s ideal beam line (%s)\r",(pkX1+dx),(pkY1+dy),vec2str(xyz1),hkl2str(hkl1[0],hkl1[1],hkl1[2])
//	Rfrom2spots(Qhkl0,xyz0,Qhkl1,xyz1,rhoTwo)					// compute rhoTwo, rotation that takes reference recip to measured hkl0 & hkl1
//	MatrixOp/O rho0 = rhoTwo x Inv(rhoRef)
//	rho0 = abs(rho0[p][q])<1e-12 ? 0 : rho0[p][q]
//	// printWave(rho0)
//	angle = axisOfMatrix(rho0,axisTwo)
//	printf "axis(rho0) = %s,  angle(rho0) = %g¡\r",vec2str(axisTwo),angle
//
//
//	print "\r  both dx=10 & dy = 0 pixels, should be a Positive Z rotation"
//	dx=10;	dy=0
//	pixel2q(geo,(pkX0+dx)*groupx+startx,(pkY0+dy)*groupy+starty,xyz0)	// direction on primary reference spot (only used for single spot analysis)
//	wenge2IdealBeamline(geo,xyz0,xyz0)							// xyz0 is used for single spot analysis
//	printf "pixel(%g, %g) --> %s ideal beam line (%s)\r",(pkX0+dx),(pkY0+dy),vec2str(xyz0),hkl2str(hkl0[0],hkl0[1],hkl0[2])
//	pixel2q(geo,(pkX1+dx)*groupx+startx,(pkY1+dy)*groupy+starty,xyz1)
//	wenge2IdealBeamline(geo,xyz1,xyz1)							// measured direction of secondary spot
//	printf "pixel(%g, %g) --> %s ideal beam line (%s)\r",(pkX1+dx),(pkY1+dy),vec2str(xyz1),hkl2str(hkl1[0],hkl1[1],hkl1[2])
//	Rfrom2spots(Qhkl0,xyz0,Qhkl1,xyz1,rhoTwo)					// compute rhoTwo, rotation that takes reference recip to measured hkl0 & hkl1
//	MatrixOp/O rho0 = rhoTwo x Inv(rhoRef)
//	rho0 = abs(rho0[p][q])<1e-12 ? 0 : rho0[p][q]
//	// printWave(rho0)
//	angle = axisOfMatrix(rho0,axisTwo)
//	printf "axis(rho0) = %s,  angle(rho0) = %g¡\r",vec2str(axisTwo),angle
//
//
//	print "\r  rotate peaks about detector center 1¡, should be a Negative F rotation"
//	RotatePixelCenter(1, pkX0,pkY0, dx,dy)
//	pixel2q(geo,(pkX0+dx)*groupx+startx,(pkY0+dy)*groupy+starty,xyz0)	// direction on primary reference spot (only used for single spot analysis)
//	wenge2IdealBeamline(geo,xyz0,xyz0)							// xyz0 is used for single spot analysis
//	printf "pixel(%g, %g) --> %s ideal beam line (%s)\r",(pkX0+dx),(pkY0+dy),vec2str(xyz0),hkl2str(hkl0[0],hkl0[1],hkl0[2])
//	RotatePixelCenter(1, pkX1,pkY1, dx,dy)
//	pixel2q(geo,(pkX1+dx)*groupx+startx,(pkY1+dy)*groupy+starty,xyz1)
//	wenge2IdealBeamline(geo,xyz1,xyz1)							// measured direction of secondary spot
//	printf "pixel(%g, %g) --> %s ideal beam line (%s)\r",(pkX1+dx),(pkY1+dy),vec2str(xyz1),hkl2str(hkl1[0],hkl1[1],hkl1[2])
//	Rfrom2spots(Qhkl0,xyz0,Qhkl1,xyz1,rhoTwo)					// compute rhoTwo, rotation that takes reference recip to measured hkl0 & hkl1
//	MatrixOp/O rho0 = rhoTwo x Inv(rhoRef)
//	rho0 = abs(rho0[p][q])<1e-12 ? 0 : rho0[p][q]
//	// printWave(rho0)
//	angle = axisOfMatrix(rho0,axisTwo)
//	printf "axis(rho0) = %s,  angle(rho0) = %g¡\r",vec2str(axisTwo),angle
//
//	KillWaves/Z xyz, DoAllProper_rhoRef
//	KillWaves/Z axisTwo,axisIndex
//	KillWaves/Z rLattice,rLattice0,rhoIndex
//	KillWaves/Z W_Cross
//	KillWaves/Z DoAllProper_Qhkl0, DoAllProper_xyz0, DoAllProper_Qhkl1, DoAllProper_xyz1, DoAllProper_rhoTwo
//	KillWaves/Z Rfrom2spots_rho0, Rfrom2spots_rho1, Rfrom2spots_h1,Rfrom2spots_x1, Rfrom2spots_axis
//End
////
//Static Function RotatePixelCenter(angle, px,py, dx,dy)
//	Variable angle					// angle (degree)
//	Variable px, py
//	Variable &dx,&dy
//	Variable c=cos(angle*PI/180), s=sin(angle*PI/180)
//	Variable xc=1064.6, yc=937.02
//	Variable x1,y1, x0=px-xc,y0=py-yc
//	x1 = c*x0 - s*y0
//	y1 = s*x0 + c*y0
//	dx = x1 - x0
//	dy = y1 - y0
//End
////
//Static Function angleBetweenVectors(a,b)
//	Wave a,b
//	Variable dot = MatrixDot(a,b)/(norm(a)*norm(b))
//	dot = limit(dot,-1,1)
//	Variable angle = acos(dot)*180/PI
//	angle = abs(angle)<1e-6 ? 0 : angle
//	return angle
//End
Function TwoSpotSetup()
	Variable H=NaN, K=NaN, L=NaN
	if (exists("root:Packages:micro:Arrays3d:hkl0")==1)
		Wave hkl0=root:Packages:micro:Arrays3d:hkl0
		H=hkl0[0]; K=hkl0[1]; L=hkl0[2]
	endif
	Variable px = NumVarOrDefault("root:Packages:micro:Arrays3d:pxPrimary",NaN)
	Variable py = NumVarOrDefault("root:Packages:micro:Arrays3d:pyPrimary",NaN)
	Variable yc = NumVarOrDefault("root:Packages:micro:Arrays3d:ycPrimary",NaN)
	Prompt H,"H"
	Prompt K,"K"
	Prompt L,"L"
	Prompt px,"X pixel (approximate)"
	Prompt py,"Y pixel (approximate)"
	Prompt yc,"yc for this reference pixel (µm)"
	DoPrompt "Primary reflection",H,px,K,py,L, yc
	if (V_flag)
		return 1
	endif
	printf "Primary:  (%s) --> {%g, %g}pixel,     yc=%gµm\r",hkl2str(H,K,L),px,py,yc
	if (numtype(H+K+L+px+py+yc))
		DoAlert 0, "Invalid inputs for primary spot, nothing done"
		return 1
	else													// update globals with current vales for primary spot
		NewDataFolder/O root:Packages					// ensure Packages exists
		NewDataFolder/O root:Packages:micro			// ensure micro exists
		NewDataFolder/O root:Packages:micro:Arrays3d	// ensure Arrays3d exists
		Variable/G root:Packages:micro:Arrays3d:pxPrimary=px, root:Packages:micro:Arrays3d:pyPrimary=py, root:Packages:micro:Arrays3d:ycPrimary=yc
		Make/N=3/D/O root:Packages:micro:Arrays3d:hkl0={H,K,L}
	endif


	if (exists("root:Packages:micro:Arrays3d:hkl1")==1)
		Wave hkl1=root:Packages:micro:Arrays3d:hkl1
		H=hkl1[0]; K=hkl1[1]; L=hkl1[2]
	else
		H=NaN; K=NaN; L=NaN
	endif
	px = NumVarOrDefault("root:Packages:micro:Arrays3d:pxSecondary",NaN)
	py = NumVarOrDefault("root:Packages:micro:Arrays3d:pySecondary",NaN)
	DoPrompt "Secondary reflection",H,K,L,px,py
	if (V_flag)
		return 1
	endif
	printf "Secondary:  (%s) --> {%g, %g}pixel\r",hkl2str(H,K,L),px,py
	if (numtype(H+K+L+px+py))
		DoAlert 0, "Invalid inputs for secondary spot, only primary was changed"
		return 1
	else													// update globals with current vales for secondary spot
		Variable/G root:Packages:micro:Arrays3d:pxSecondary=px, root:Packages:micro:Arrays3d:pySecondary=py
		Make/N=3/D/O root:Packages:micro:Arrays3d:hkl1={H,K,L}
	endif
End


Function ChangeAllX0H0(X0,H0)		// resets X0 & H0 in all wave notes where there is an X0= & H0=
	Variable X0,H0
	if (numtype(X0+H0))
		X0 = numtype(X0) ? 0 : X0
		H0 = numtype(H0) ? 0 : H0
		Prompt X0, "X offset to the indent axis (µm)"
		Prompt H0, "H offset to the indent axis (µm)"
		DoPrompt "axis offset",X0,H0
		if (V_flag)
			return 1
		endif
	endif
	if (numtype(X0+H0))
		DoAlert 0, "invalid values for X0 or H0"
		return 1
	endif

	String list = WaveList("*",";",""), wnote
	Variable N=ItemsInList(list),i
	for (i=0;i<N;i+=1)
		Wave w=$StringFromList(i,list)
		wnote = note(w)
		if (keyInList("X0",wnote,"=",";") && keyInList("H0",wnote,"=",";"))
			wnote = ReplaceNumberByKey("X0",wnote,X0,"=")
			wnote = ReplaceNumberByKey("H0",wnote,H0,"=")
			Note/K w, wnote
		endif
	endfor
End


Function/T Get2RangesOfSPEfile()
	Variable f
	PathInfo reconPath
	if (V_flag)
		Open/D/M="pick first of the spe files"/P=reconPath/R/T=".SPE" f
	else
		Open/D/M="pick first of the spe files"/R/T=".SPE" f
	endif

	String out="", fullPathName=S_fileName, namePart=""
	if (strlen(fullPathName)<1)
		return out
	endif
	String pathPart = ParseFilePath(1, fullPathName, ":", 1, 0)
	if (strlen(pathPart)<1)
		return out
	endif
	out = ReplaceStringByKey("path",out,pathPart,"=")
	String firstName = ParseFilePath(3, fullPathName, ":", 0, 0)

	Variable i, range0=NaN,depth0=Nan, N1=strlen(firstName)-1
	i = strsearch(firstName,"_",N1,1)
	if (i<0)
		return out
	endif
	depth0 = str2num(firstName[i+1,N1])
	if (numtype(depth0))
		return out
	endif
	namePart = firstName[0,i]							// possible name part, includes final underscore
	firstName = firstName[0,i-1]						// trim of the final _124
	N1 = i-1

	i = strsearch(firstName,"_",N1,1)
	if (i<0)
		range0 = depth0									// no depths, just a simple range
		depth0 = NaN
	else
		range0 = str2num(firstName[i+1,N1])
		if (numtype(range0))
			range0 = depth0								// no depths, just a simple range
			depth0 = NaN
		else
			namePart = firstName[0,i]					// name part, includes final underscore
		endif
	endif
	out = ReplaceStringByKey("name",out,namePart,"=")
	String range = SelectString(numtype(range0),num2istr(range0),"")
	String depths = SelectString(numtype(depth0),num2istr(depth0),"")

	NewPath/O/Q Get2RangesOfSPEfilePath pathPart
	String list = IndexedFile(Get2RangesOfSPEfilePath,-1,".SPE")	// all SPE files
	KillPath/Z Get2RangesOfSPEfilePath
	String fName
	Variable m,N=ItemsInList(list),rr,dd=NaN
	Variable range1=range0,depth1=depth0

	for (m=0;m<N;m+=1)
		fName = StringFromList(m,list)
		if (strsearch(fName,namePart,0)!=0)			// must start with namePart
			continue
		endif
		if (numtype(depth0))							// no depths, only range
			sscanf fName, namePart+"%d", rr
			if (V_flag!=1)
				rr = NaN
			endif
		else
			sscanf fName, namePart+"%d_%d", rr,dd
			if (V_flag!=2)
				rr = NaN
				dd = NaN
			endif
		endif
		range1 = (rr>range1) ? rr : range1
		depth1 = (dd>depth1) ? dd : depth1
	endfor
	range += SelectString(range1>range0,"","-"+num2istr(range1))
	depths += SelectString(depth1>depth0,"","-"+num2istr(depth1))
	out = ReplaceStringByKey("range",out,range,"=")
	out = ReplaceStringByKey("depths",out,depths,"=")
	return out
End

Static Function/T ProgressPanel(percentName,[stop])	// display a progress bar
//	if the user does a "DoUpdate/W=$wname", then V_flag will be set to 2 if the stop button was pushed
	String percentName									// name of global value to use, if not a global variable then set value to 0
	Variable stop											// if true, include a stop button
	stop = ParamIsDefault(stop) ? 0 : stop
	stop = NumberByKey("IGORVERS",IgorInfo(0))>=6.10 ? stop : 0	// stop feature requires version 6.10 or above

	Variable right = stop ? 695 : 630
	NewPanel/K=1/W=(330,87,right,117)/N=Progress
	String wname = S_name								// save the name of the window
	ValDisplay progressValdisp, win=$wname, pos={5,5},size={288,25},font="Helvetica",fSize=18,format="%d%%"
	ValDisplay progressValdisp, win=$wname, limits={0,100,0},barmisc={0,50},highColor= (19368,32650,65196)
	if (exists(percentName)==2)						// the global exists, use it
		NVAR percent = $percentName
		percent = limit(percent,0,100)
		ValDisplay progressValdisp,win=$wname, value= #percentName
	else													// no global, just set to zero
		ValDisplay progressValdisp,win=$wname, value= #"0"
	endif

	if (stop)
		Button stopButton,pos={299,4},size={50,20},title="Stop"
	endif
#if NumberByKey("IGORVERS",IgorInfo(0))>=6.10
	DoUpdate/W=$wname/E=1							// mark as a progress window
#else
	DoUpdate
#endif
	return wname
End
//	ValDisplay progressValdisp,win=$wname, value= #"30"		// example of how to set value

//  ================================ End of Single Spot ===============================  //
//  ============================================================================  //