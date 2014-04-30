#pragma rtGlobals=1		// Use modern global access method.
#pragma ModuleName=Triangulation
#pragma version = 2.11
#include "microGeometry", version>=2.1


Menu "micro"
	SubMenu "Refinement"
		"Triangulate...",Triangulate("")
	End
	SubMenu "Tables"
		"Pixel Rays for Triang",/Q,EditTriangPixelRays($"")
	End
End


Function plot_all_rays()
	if (exists("Triangulate_rays")==1)
		Wave rays = Triangulate_rays
	elseif (exists("TriPixelRays")==1)
		Wave rays = TriPixelRays
	else
		return 1
	endif
	Variable i, N=DimSize(rays,1)
	display
	for (i=0;i<N;i+=1)
		appendtograph rays[][i][1]
		appendtograph/r rays[][i][2]
	endfor

	String wname, list=TraceNameList("", ";", 1)
	N = ItemsInList(list)
	for (i=0;i<N;i+=1)
		wname = StringFromList(i,list)
		if (mod(i,2))
			ModifyGraph rgb($wname)=(0,0,65280)
		endif
	endfor
	ModifyGraph/Z mode=4,marker=19,msize=1, lowTrip=0.001
	ModifyGraph/Z tick=2, minor=1, standoff=0, mirror(bottom)=1
End



Function Triangulate(range)
	String range
	if (!(ItemsInRange(range)>0))
		Prompt range,"range of file numbers to triangulate"
		DoPrompt "range",range
		if (V_flag)
			return 1
		endif
		printf "Triangulate(\"%s\")\r",range
	endif
	Variable N=ItemsInRange(range)
	if (!(N>0))
		DoAlert 0,"No items in range '"+range+"'"
	endif
	String prefixList = "FullPeakList_tri_"

	String strStruct=StrVarOrDefault(":geoStructStr","")		// get to values in current directory
	if (strlen(strStruct)<1)
		strStruct=StrVarOrDefault("root:Packages:geometry:geoStructStr","")	// try the default values
	endif
	if (strlen(strStruct)<1)
		Abort "you must set the geomery parameters first"
	endif
	STRUCT microGeometry geo
	StructGet/S/B=2 geo, strStruct							// load previous values
	Variable xc=geo.xcent, yc=geo.ycent, ddOffset=geo.ddOffset
	xc -= 1													// change from 1 based to 0 based
	yc -= 1

	Variable m, allHere
	for(m=str2num(range),allHere=1; !numtype(m)&&allHere; m=NextInRange(range,m))	// loop over every image
		allHere = (exists(prefixList+num2istr(m))==1)
	endfor
	if (allHere)
		DoAlert 1,"Use the already loaded data?"
		allHere = (V_flag==1)
	endif
	if (!allHere)											// load and process each image
		Variable fid
		Open/D/R/P=home/T=".SPE"/M="pick one of the spe files" fid
		if (strlen(S_fileName)<1)
			DoAlert 0,"No file chosen"
			return 1
		endif

		String fileRoot = ParseFilePath(1,S_fileName,":",1,0)+ParseFilePath(3,S_fileName,":",1,0)
		Variable i, ichar, zero=char2num("0"), nine=char2num("9")
			for (i=strlen(fileRoot)-1;i>=0 && isdigit(fileRoot[i]);i-=1)
		endfor
		fileRoot = fileRoot[0,i]
		printf "Triangulating with '%s',  range=\"%s\"\r",fileRoot,range
		for(m=str2num(range); !numtype(m); m=NextInRange(range,m))	// loop over every image
			if (ProcessOneTriangImage(fileRoot,m,prefixList))
				DoAlert 0,"Unable to process image "+num2istr(m)
				return 1
			endif
			printf "processed image %s%d.SPE\r",fileRoot,m
		endfor
	endif
	// done loading and processing all images (delted them too), have left the fitted peak lists

	Wave fitted = $(prefixList+num2istr(LastInRange(range)))
	Variable Nrays = round(DimSize(fitted,0)*0.7)
	printf "doing triangulation on %d different rays\r",Nrays

	Make/N=3/O/D Triangulate_xyz
	Wave xyz = Triangulate_xyz
	Make/N=(N,Nrays,5)/O/D Triangulate_rays, TriPixelRays	// only used for looking at, nothing really uses this wave
	Wave rays = Triangulate_rays
	SetDimLabel 0,-1,'height',TriPixelRays
	String str
	for(i=0;i<Nrays;i+=1)
		sprintf str,"spot[][%d][]",i
		SetDimLabel 1,i,$str,TriPixelRays
	endfor
	Note/K TriPixelRays, "waveClass=TriPixelRays"
	rays = NaN
	TriPixelRays = NaN
	Wave fitted = $(prefixList+num2istr(LastInRange(range)))// using the last image (farthest)
	Variable j
	for (j=0;j<Nrays;j+=1)
		rays[N-1][j][0] = NumberByKey("CCDy",note(fitted),"=")+ddOffset	// height of this peak
		rays[N-1][j][1] = fitted[j][0]							// fitted peak positions PIXELS not mm, just temporary
		rays[N-1][j][2] = fitted[j][1]
	endfor

	Variable geo_dd_save = geo.dd
	Variable px,py, height, mx,my, imin
	for (j=0;j<Nrays;j+=1)									// for each of the rays
		height = rays[N-1][j][0]								//use this pixel position to get approx line
		px = rays[N-1][j][1]									//  for finding the other fitted peak positions to use
		py = rays[N-1][j][2]
		// px,py is now the pixel on the last image
		// ray goes from (xc,yc) at height=0   to (px,py) at height=CCDy+ddOffset  (note xc,yc are only approximate here)
		// equations are px = mx*height + xc  and  py = my*height + yc
		mx = (px-xc)/height									// slope for x-pixels, not exact, just used to find right peak
		my = (py-yc)/height									// slope for y-pixels
		for(i=0,m=str2num(range); !numtype(m); m=NextInRange(range,m),i+=1)	// loop over every FullPeakList wave
			Wave fitted = $(prefixList+num2istr(m))			// find closest pixel in fitted peaks wave
			height = NumberByKey("CCDy",note(fitted),"=")+ddOffset
			sprintf str,"CCDy=%.1f",height
			SetDimLabel 0,i,$str,TriPixelRays
			px = mx*height + xc
			py =my*height + yc
			imin = ClosestPixelInPeakWave(fitted,px,py)		// closest pixel in fitted wave
			if (abs(mx*height+xc-px)+abs(my*height+yc-py) < 10)	// if not within 10 pixels, skip it
				geo.dd = NumberByKey("CCDy",note(fitted),"=")// equivalent to ddOffset=0, that is what we need to find
				pixel2XYZ(geo,fitted[imin][0],fitted[imin][1],xyz)	// convert from pixel to xyz (Wenge Coords)
				rays[i][j][0] = xyz[2]							// height (mm)
				rays[i][j][1] = xyz[0]							// x-direction on CCD (mm)
				rays[i][j][2] = xyz[1]							// y-direction on CCD (mm)
				rays[i][j][3] = fitted[imin][2]*geo.dpsx/geo.NxCCD	// rescale error from pixels to mm
				rays[i][j][4] = fitted[imin][3]*geo.dpsy/geo.NyCCD
				TriPixelRays[i][j][0] = geo.dd
				TriPixelRays[i][j][1,4] = fitted[imin][r-1]
			endif
		endfor
	endfor
	WaveStats/Q/M=1 rays
	if (V_numNaNs)
		printf "skipped %d points\r",V_numNaNs
	endif

	Make/N=(N)/O/D Triangulate_height, Triangulate_xy, Triangulate_weight
	Make/N=3/O/D Triangulate_xWave, Triangulate_typXWave
	Wave xWave=Triangulate_xWave, typXWave=Triangulate_typXWave
	xWave = {0,0,ddOffset}										// initial guess
	typXWave = {0.1,0.1,0.01}									// typical values, used for scaling
	Variable ystart = TriangNeckErr(rays,xWave[0],xWave[1],xWave[2])
	printf "before, geo.xc = %.3f,  geo.yc=%.3f,  geo.ddOffset=%.4f\r",geo.xcent, geo.ycent, geo.ddOffset
	printf "at start of optimization, TriangNeckErr(rays,xo=%g, yo=%g, ddOffset=%.4f) = %.9g\r",xWave[0],xWave[1],xWave[2],ystart
	Optimize/Q/M={0,0}/R=typXWave/X=xWave/Y=(ystart) TriangNeckErr,rays
	if (V_flag)
		printf "Optimize terminated improperly, V_flag=%g,   V_OptTermCode=%g\r",V_flag,V_OptTermCode
	endif
	// printf "after optimize, xWave = %s and TriangNeckErr()=%g\r",vec2str(xWave),TriangNeckErr(rays,xWave[0],xWave[1],xWave[2])
	ddOffset = xWave[2]
	xyz[0] = xWave[0]											// xyz of center after optimization
	xyz[1] = xWave[1]
	xyz[2] = geo.dd - ddOffset

	Variable/C pz = xyz2pixelUn(geo,xyz)						// returns 1 based undistorted pixels (the new xc,yc)
	xc = real(pz)
	yc = imag(pz)
	printf "W_OptGradient = %s\r",vec2str(W_OptGradient)
	printf "at end of optimization, TriangNeckErr(rays, xc=%.3f, yc=%.3f, ddOffset=%.4g) = %.9g\r",xc,yc,ddOffset,V_min
	printf "changed  Æxc=%.3f pixel,   Æyc=%.3f pixel,   ÆddOffset=%.4f mm\r",(xc-geo.xcent),(yc-geo.ycent),(ddOffset-geo.ddOffset)
	geo.dd = geo_dd_save										// restore dd
	KillWaves/Z Triangulate_height, Triangulate_xy, Triangulate_weight, Triangulate_xWave, Triangulate_typXWave
	KillWaves/Z W_coef,W_sigma, W_OptGradient
	KillWaves/Z Triangulate_rays, Triangulate_xyz
End
//
// returns rms value of distance from (x0,y0) for each ray projected to height=0
Function TriangNeckErr(rays,x0,y0,ddOffset)	// used by Optimize in Triangulate()
	Wave rays
	Variable x0,y0,ddOffset
	Wave height=Triangulate_height, xy=Triangulate_xy, weight=Triangulate_weight
	Variable j, Nrays=DimSize(rays,1), err
	for (err=0,j=0;j<Nrays;j+=1)
		height = rays[p][j][0]+ddOffset
		xy = rays[p][j][1]				// first the x pixel
		weight = rays[p][j][3]
		CurveFit/Q line xy /X=height/W=weight/I=1
		Wave W_coef=W_coef, W_sigma=W_sigma
		err += ((W_coef[0]-x0)/W_sigma[0])^2
		xy = rays[p][j][2]				// next the y pixel
		weight = rays[p][j][4]
		CurveFit/Q line xy /X=height/W=weight/I=1
		err += ((W_coef[0]-y0)/W_sigma[0])^2
	endfor
	return sqrt(err/(2*Nrays))
End
//
Static Function ClosestPixelInPeakWave(w,px,py)
	Wave w
	Variable &px,&py

	Variable dist, i, imin, N=DimSize(w,0)
	String class = StringByKey("waveClass",note(w),"=")
	if (!stringmatch(class,"FittedPeakList"))
		Abort "ClosestPixelInPeakWave() does not undstand the waveClass = "+class
	endif

	for (i=0,dist=Inf;i<N;i+=1)						// just find the closest pixel
		if ((w[i][0]-px)^2+(w[i][1]-py)^2 < dist)
			dist = (w[i][0]-px)^2+(w[i][1]-py)^2
			imin = i
		endif
	endfor
	px = w[imin][0]
	py = w[imin][1]
	return imin
End
//
Static Function ProcessOneTriangImage(fileRoot,fileNum,prefixList)
	String fileRoot											// full pathname of file, up to "4.SPE" part
	Variable fileNum
	String prefixList										// prefixes for fitted peaks

	Wave image = $LoadWinViewFile(fileRoot+num2istr(fileNum)+".SPE")
	if (!WaveExists(image))
		DoAlert 0, "unable to load image "+num2istr(fileNum)
		return 1
	endif
	Wave imageBkg = $FitPeaks(image,0.8,30)
	if (!WaveExists(imageBkg))
		DoAlert 0, "unable to find background image "+num2istr(fileNum)
		return 1
	endif
	Wave FullPeakList=FullPeakList
	String str = prefixList+num2istr(fileNum)
	KillWaves/Z $str
	Rename FullPeakList $str
	Wave FullPeakList = $str

 	Killwaves/Z image, imageBkg
	return 0
End



//  CAUTION, this is a special routine, and should not ordinarily be used.  See microGeometry.ipf for the usual routines.
// convert pixel location xyz[3](µm) to pixel values,  NOTE the returned pixel is NOT distortion corrected, and 1 based
// since there is no guarantee that xyz lies on the detector, treat xyz as a direction, not a position
Static Function/C xyz2pixelUn(geo,xyz)						// returns pixel values as a complex number cmplx(px,py)
	STRUCT microGeometry &geo
	Wave xyz										// absolute position of pixel in Wenge coords (a 3-vector) (µm), origin at Si
	Variable px,py									// final pixel position, returned (zero based pixels)

	Wave vecB=root:Packages:geometry:q2pixel_vecB
	Wave mat=root:Packages:geometry:q2pixel_mat
	Wave kout=root:Packages:geometry:q2pixel_kout
	kout = xyz
	normalize(kout)

	//	k k^ = a px^ + b py^ + dd {001}				// need to solve for a & b, (know all unit vectors and dd)
	//  -dd {001} = a px^ + b py^ - k k^			// rewrite this way to solve for a, b, k, (I know dd)
	// vecB = mat x {a,b,k}
	mat[0][0] = geo.rded00
	mat[1][0] = geo.rded01
	mat[2][0] = geo.rded02
	mat[0][1] = geo.rded10
	mat[1][1] = geo.rded11
	mat[2][1] = geo.rded12
	mat[][2] = -kout[p]
	vecB = {0,0,-geo.dd}
	MatrixLUD mat
	MatrixLUBkSub M_Lower, M_Upper, W_LUPermutation, vecB 
	Wave M_x = M_x
	M_x[0] = CCD_X_REVERSE ? -M_x[0] : M_x[0]
	px = geo.xcent + M_x[0]*geo.NxCCD/geo.dpsx	// this formula is for 1 based pixels
	py = geo.ycent - M_x[1]*geo.NyCCD/geo.dpsy
	KillWaves/Z M_Lower,M_Upper,W_LUPermutation,M_x
	return cmplx(px,py)
End



Function EditTriangPixelRays(pixelRays)				// table of pixel positions used to make xyz positions for triangulation
	Wave pixelRays

	if (!WaveExists(pixelRays))
		String list =WaveListClass("TriPixelRays","*","")
		if (stringmatch(list,"TriPixelRays;"))			// only one choice, use it
			Wave pixelRays = :TriPixelRays
		else
			String wName=""
			Prompt wName,"TriPixelRays",popup,list
			DoPrompt "TriPixelRays",wName
			if (V_flag)
				return 1
			endif
			Wave pixelRays = $wName
		endif
	endif
	if (!WaveExists(pixelRays))
		DoAlert 0, "TriPixelRays does not exist in this data folder"
		return 1
	endif
	String table = Refinement#FindTableWithWave(pixelRays)
	if (strlen(table))
		DoWindow/F $table
	else
		Edit/K=1/W=(5,44,950,297) pixelRays.ld
		ModifyTable width(Point)=42,format(pixelRays.d)=3,digits(pixelRays.d)=4,width(pixelRays.d)=114
		ModifyTable width(pixelRays.l)=84
	endif
	return 0
End

//Static Function isdigit(c)		// moved to Utility_JZT
//	String c
//	Variable i=char2num(c)
//	return (48<=i && i<=57)
//End


