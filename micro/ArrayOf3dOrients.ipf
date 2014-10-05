#pragma rtGlobals=1		// Use modern global access method.
#pragma version=2.46	// removed all of the old stuff, for 1.6 added the tensor parts, 2.0 changed surfer->gizmo
#pragma ModuleName=ArrayOf3dOrients
#include "DepthResolvedQuery",version>=1.16
#include "ImageDisplayScaling", version>=1.47
#include "microGeometry", version>=2.51
//
//
//	Strconstant micron ="µm"			// use for Mac, I think that this works for Windows too!
//	Strconstant micron ="micron"		// use for Win
//	mu is hex B5 = 181
//	deg is hex A1 = 161
//	Windows deg is hex C1 = 193


// July 10, 2006		version 1.7
//  Added the gradient of the total rotation vector.  So you can now also plot any of the 3 components of the gradient or the magnitude.
//  This uses the new function GradtRotationAnglesFromRods3d().  which calculates the gradient of the total rotaion angle.  Not
//  the components (as used in the dislocation tensor).
//
//
// July 10, 2006		version 1.6
//  Added the dislocation tensor (and GND) capability.  So in addition to rotation, you can now plot any of the 9 
//  componensts of the dislocation tensor, and the GND.   For the GND, the scaling is set by the (optonal) 
//  global variable 'GND_DislocationDensity'.  Which should be in the current datafolder or in the root data folder.
//
//
// Jan 26, 2006
//  Changed SlicePlaneIn3d(resolution,normalW,SliceValue,rotation)
//  The components of SliceAngleAxis were incorrect.  Since the matricies are in "sample" coordinate 
//  system (X,-H,-F), the correct unit vectors for the axes are:
//  Hhat=(0,-1,0);  Fhat=(0,0,-1);  Xhat=(1,0,0);  Yhat=(0,-1,1)/sqrt(2);  Zhat=(0,-1,-1)/sqrt(2)
//
//
// Aug 16, 2006
//  Changed the graph background to light gray so that NaN will show through
//
//
// Jan 24, 2007
//  Changed SlicePlaneIn3d() to also make an RGB image showing orientation relative to 001, 011 and 111.
//  also changed CutSlicePopUpProc() to supprot sliceOrients[][]
//  also changed SetSliceProc() to support sliceOrients[][]
//  added MakeGraphOfSliceOrients() to plot sliceOrients[][]
//
//
// Jan 29, 2007
//  added waveClass to "sliceWave" and "sliceOrients"
//  also changed CutSlicePopUpProc() to supprot sliceOrients[][]
//  also changed SetSliceProc() to support sliceOrients[][]
//  added MakeGraphOfSliceOrients() to plot sliceOrients[][]
//
//
// Feb 15, 2007
//  changed readInNewRawDataFile() so that the input data are recip lattices not rotation matricies.
//  changed readInNewRawDataFile() so that bad points are skipped
//  added "micrGeometry.ipf", so that YZ -> HF uses functions that permit angles other than 45°
//  added MakeGraphOfSliceOrients() to plot sliceOrients[][]
//
// Jun 7, 2007
// Changed SlicePlaneIn3d() & Interpolate3dRodriquesFromXHF()
// made it work with R2dX, R2dH, R2dF for a single X
//	added InterpolateRodriquesFromXHF() to call the 2d or 3d interpolation
//  added maxAngle to readInNewRawDataFile()
//
// Jun 11, 2007
// Changed Surfer to Gizmo
// needed to change:   SetSliceProc(), MakeGraph3dSurface(), Set3dSurfaceZscale()
// also changed from using readInNewRawDataFile() to read3dRecipLatticesFile(), so automaticaly deals with Rodriques or recip
//
// Jun 20, 2007,  version 2.02
// added interpolate3dRodFromXHFSheets(), for interpolating with few X values, changed InterpolateRodriquesFromXHF() too.
// in all three of the interpolation routines, bad points are set to NaN (used to be zeros)
// also fixed read3dRecipLatticesFile(), the definition of rho was wrong
//
// Jun 22, 2007,  version 2.03
// include intensity of the pattern in the information passed around.  Also the slice wave plots can now show the intensity as
// a darkening of the image.
// changed SlicePlaneIn3d(), sliceWave2RGB(), interpolate3dRodFromXHFSheets(), read3dRecipLatticesFile()
//
// Jun 25, 2007,  version 2.05
// when using intensity, make un-indexed pixels green but show the intensity
// changed sliceWave2RGB(), and read3dRecipLatticesFile()
//
// Jun 27, 2007,  version 2.06, 2.07
// added the ability to also view strain (epsilon) information
// changed interpolate3dRodFromXHFSheets(), SurfacePlotStyle()
//
// Jul 3, 2007,  version 2.09
// add kappa (lattice curvature tensor), change calc of alpha to be kappa^t - Tr(kappa)*I - (del x epsilon)
// added curvatureTensorFromRods3d(), changed dislocationTensorFromRods3d(), SlicePlaneIn3d(), SurfacePlotStyle(), MakeDirectionVector()
//
// Jul 5, 2007,  version 2.10
// fixed FixAspectRatio(), forgot graphName in the ModifyGraph
// removed FixAspectRatio(), replace with SetAspectToSquarePixels() in ImageDisplayScaling.ipf
//
// Jul 5, 2007,  version 2.10
// fixed interpolate3dRodFromXHFSheets(), skip sheets with less than 4 points at one X (cannot interpolate)
//
// Jul 15, 2007,  version 2.12
// added ability to ignore epsilon when calculating alpha or GND
// added epsilonZeroProc(), and modified, dislocationTensorFromRods3d(), MakePanelSlicer(), SetSliceProc(), SurfacePlotStyle(), SlicePlaneIn3d()
//
// Jul 15, 2007,  version 2.13
// fixed calculation of the alphas, had a typo for calc of exy0, exz0, eyz0
// fixed small bug in FindScalingFromVec()
//
// Jul 30, 2007,  version 2.14
// added graphing of 9 tensor components in one plot
// added routines GraphTensorcomponents(), ReCalcAllTensorComponents(), & DisplayBoxFor3x3fromAB()
// also added Change_GND_DislocationDensity() as a convienence
//
// Aug 2, 2007,  version 2.15
// changed SlicePlaneIn3d() to prefer totalIntensity over totalPeakIntensity
// changed sliceWave2RGB() , use trim=0.1 for picking intensity range, this brightens the RGB image
//
// Aug 2, 2007,  version 2.16
//	removed R3dX,R3dH,R3dF from Fixed dislocationTensorFromRods3d()
//	fixed a bug in axisOfMatrix(), it was coming out negative, reversed some things and it now looks right
//	changed read3dRecipLatticesFile(), so that RX,RH,RF is the rotation axis from axisOfMatrix()
//	and little fix on InterpolateRodriquesFromXHF(), changed the order of the if(...)
//
// Aug 10, 2007,  version 2.17
//	fixed Set3dSurfaceZscale() to work for Gizmo
//
// Aug 19, 2007,  version 2.18
//	added DuplicateTheXs() so that data sets with only one X value will appear to have two X values,
//  also changed readInXYZorientationsFromFile() to call DuplicateTheXs()
//
// Aug 27, 2007,  version 2.19
//	in "Rotations" menu, enabled z-scaling (I forgot to do this earlier)
//
// Sept 26, 2007,  version 2.21
//	fixed rotationMatAboutAxis(), it now agrees with the C code, and I checked it with axisOfMatrix(), applied sequentially you get back where you started
//	changed rotationMatAboutAxis(), now can use length of axis[] to give angle
//	added AngleBetweenRotationVectors(), finds angle between two Rodriques vectors
//
// Oct 1, 2007,  version 2.22
//	in read3dRecipLatticesFile(), added the array indicies.  It is needed for backtracking, used by findClosestPoint()
//
// Oct 24, 2007,  version 2.23
//	added FillArrays3dParametersPanel(), EnableDisableArrays3dControls(), Arrays3dButtonProc(), testing3dArraysPopMenuProc()
//		all to make the functions here available on the micro panel.
//	changed RodriquesVector2str() to work with the micro panel
// 	ALSO, added waveClasses to the waves that are read in and interpolated
//
// Oct 26, 2007,  version 2.24
//	modified FillArrays3dParametersPanel(), EnableDisableArrays3dControls(), Arrays3dButtonProc(), Trim3dArraysPopMenuProc()
//	added buttons for trimming surface and bottom, the added funcitons are all in IndexLots.ipf
//	Also changed FindScalingFromVec() once more
//
// Oct 29, 2007,  version 2.26
//	added optional method argument to InterpolateRodriquesFromXHF()
//	added TrimBadPointsOnOutside() and removeAllBadPoints()
//	small change to FindScalingFromVec(), added the part with the round(stepSize*25)
//
// Nov 7, 2007,  version 2.27
//	changed TableOfIndexingList()
//
// Nov 10, 2007,  version 2.27
//	modified DisplayBoxFor3x3fromAB() to work for 3 angle plots
//	added Graph3AngleComponents() and ReCalc3AngleComponents()
//
// Jan 13, 2008,  version 2.29
//	modified DisplayBoxFor3x3fromAB() to work for 3 angle plots
//	added Graph3AngleComponents() and ReCalc3AngleComponents()
//
// Jan 16, 2008,  version 2.31
//	modified FillArrays3dParametersPanel(), EnableDisableArrays3dControls(), Arrays3dButtonProc(), to add "mark reference points" button
//
// Feb 7, 2008,  version 2.32
//	modified rotationMatAboutAxis(), so that it check dimensions of mat, but does not re-dimension it
//
// Feb 8, 2008,  version 2.33
//	moved rotationAngleOfMat(), axisOfMatrix(), AngleBetweenRotationVectors(), rotationMatAboutAxis(), notRotationMatrix(), and SquareUpMatrix() to microGeometry.ipf
//
// Feb 12, 2008,  version 2.34
//	changed Set3dSurfaceZscale() and SetSliceProc() so that colors on gizmo plot are right when you zoom z-scale
//	modified MakeGraph3dSurface() so that the grid s a bit translucent, had to add a blend too.
//
// Mar 6, 2008,  version 2.36
//	changed SurfacePlotStyle(), deal with situation when data does not come from a file
//
// Mar 6, 2008,  version 2.37
//	made phi range [-360, 360]
//
// Mar 13, 2008,  version 2.37
//	added MakeOrientationsGizmo() and rawCalcToGizmo()
//
// Mar 27, 2008,  version 2.39
//	changed InterpolateRodriquesFromXHF so that method is handled correctly
//
// May 15, 2008,  version 2.40
//	changed Interpolate3dRodriquesFromXHF so that the intensities are also interpolated
//
// Aug 21, 2008,  version 2.41
//	changed interpolate3dRodFromXHFSheets(), so that it handles a non-uniform x-spacing properly.
//
// Apr 22, 2009,  version 2.42
//	changed step size in PanelSlicer so that the Interp3d did not get NaN
//
// Apr 23, 2009,  version 2.43
//	changed FindScalingFromVec() to use median instesad of average to determine first point
//
// Apr 25, 2009,  version 2.44
//	changed SetSliceProc() so that 3d zscale is sticky
//
// Apr 30, 2014,  version 2.45
//	changed line terminations CR -> LF
//
// Oct 5, 2014,  version 2.46
//	added squareUp parameter to axisOfMatrix()

// modified read3dRecipLatticesFile() to deal with multiple reference matricies
// and added SetMatrixFromString() to help with this task.
//
//	 e.g. Variable/G GND_DislocationDensity=2e10
//
Constant GND_DislocationDensity_Al = 77.92e10	// converts from (rad/µm) to (dislocations/cm^2)
	// I was told (1°/µm) = 1.36e10 (dislocations/cm^2) which should be 1.36e10*180/π = 77.92e10 for (rad/µm)
	// for Al with a=4.05Å, I would guess that GND_DislocationDensity_Al=1e4/4.05e-8 = 24.691e10,  so where does 77.92e10 come from?

Static Constant GRAY = 55000			// used to set color for points that did not index gray = {GRAY,GRAY,GRAY}}


// change rotationAxis[3]

Menu "Rotations", dynamic
	"Plot slice of 3d array",MakeGraphOfSlice()
	"Plot slice of 3d array with intensity",MakeGraphOfSliceWaveRGB()
	"Plot slice of RGB orientations",MakeGraphOfSliceOrients()
	"put up PanelSlicer", /Q,MakePanelSlicer(sliceWave)
	Submenu "3d surface"
		"Show 3d surface",MakeGraph3dSurface()
		"  autoscale z",Set3dSurfaceZscale(Inf)
		"  fix z scale",Set3dSurfaceZscale(NaN)
//		"( autoscale z",ModifySurfer autoscale=1
//		"( fix z scale",Set3dSurfaceZscale()
		BeamDirectionMenuItem(0) ,/Q, SetBeamDirection(0)
		BeamDirectionMenuItem(1) ,/Q, SetBeamDirection(1)
	End
	Submenu "   Read & Interpolate"
		"2nd — Interpolate 3d grid",InterpolateRodriquesFromXHF(NaN)
		"1st — Read In Raw Data",readInXYZorientationsFromFile(NaN)
	End
	"-"
	MenuItemIfWaveClassExists("Graph of all 9 components","TensorOrientationSlice",""),/Q,GraphTensorcomponents("")
	MenuItemIfWaveClassExists("   Re-Calc all 9 tensor components","OrientationSliceWave",""),/Q,ReCalcAllTensorComponents("")
	//	"Graph of all 9 components",/Q,GraphTensorcomponents("")
	//	"   Re-Calc all 9 tensor components",/Q,ReCalcAllTensorComponents("")
	"Re-Set GND density conversion factor",Change_GND_DislocationDensity(NaN)
	"-"
	Submenu "   testing"
		"Table of IndexingLIst",TableOfIndexingList($"")
		"print a Rodriques vector",RodriquesVector2str(NaN,NaN,NaN)
		"Table of Raw Data",TableOfRawData()
		"Make test data",MakeTestData("")
	End
	"-"
	"Cut Off One Side",CutOffOneSide(NaN)
	"Show Cubic Color Triangle",showCubicTriangleColors(NaN,NaN)

	Submenu "Gizmo of Raw Rotations"
		MenuItemIfWaveClassExists("Gizmo of {RX,RH,RF}","Orient3dRawGizmo",""),MakeOrientationsGizmo(NaN)
		MenuItemIfWaveClassExists("     Re-Calc Orients for Gizmo","Random3dArrays",""),rawCalcToGizmo(NaN)
	End



End
Function/S BeamDirectionMenuItem(m)			// note, "!\022(" causes a menu item to be checked and disabled
	Variable m
	Variable/G :reverse3dSurfacePlot
	NVAR flag = :reverse3dSurfacePlot
	if (m==0)
		return SelectString(flag,"!\022(","")+"beam from Right «--"
	else
		return SelectString(flag,"","!\022(")+"beam from Left   --»"
	endif
End
Function SetBeamDirection(flag)
	Variable flag
	NVAR rev = :reverse3dSurfacePlot
	rev = flag
	if (strlen(WinList("PanelSlicer","","WIN:64")))					// update PanleSlicer if showing
		String str = SelectString(flag,"--» beam from Left","«-- beam from Right")
		Execute "PopupMenu popupBeamDirection win=PanelSlicer, value= "+"\""+str+"\""
		str = SelectString(flag,"«-- beam from Right","--» beam from Left")
		PopupMenu popupBeamDirection,title=str
	endif
	SetSliceProc("", NumVarOrDefault("SliceValue",0),"","SliceValue")
End

// ************************************************************************************************************************
// ************************************************************************************************************************
//
// calculates the dislocation tensor (alpha), and GND from alpha
Function dislocationTensorFromRods3d(xx,hh,ff,exx,eyy,ezz,exy,exz,eyz,kappa,alpha)
	Variable xx,hh,ff					// desired xyz position
	Wave exx,eyy,ezz,exy,exz,eyz		// components of the strain tensor (if these are null waves, assume strain is zero)
	Wave kappa							// Lattice Curvature Tensor
	Wave alpha							// Dislocation Density Tensor

	//	alpha = kappa^T - Tr(kappa)*I - (del x epsilon)
	// (del x epsilon) = SUM{ e(ikl)*d[epsilon(lj)]/dk }

	Variable trace=MatrixTrace(kappa)
	if (numtype(trace))												// cannot calculate
		alpha = NaN
		return NaN
	endif
	Variable dexx_dX=0,dexy_dX=0,dexz_dX=0,deyy_dX=0,deyz_dX=0,dezz_dX=0,deyx_dX=0,dezx_dX=0,dezy_dX=0
	Variable dexx_dY=0,dexy_dY=0,dexz_dY=0,deyy_dY=0,deyz_dY=0,dezz_dY=0,deyx_dY=0,dezx_dY=0,dezy_dY=0
	Variable dexx_dZ=0,dexy_dZ=0,dexz_dZ=0,deyy_dZ=0,deyz_dZ=0,dezz_dZ=0,deyx_dZ=0,dezx_dZ=0,dezy_dZ=0
	if (WaveExists(exx) )											// compute curl of epsilon if epsilons exist, otherwise, curl is zero
		Variable dX=DimDelta(exx,0), dY=DimDelta(exx,1), dZ=DimDelta(exx,2)
		if (xx+dX == DimOffset(exx,0) + (DimSize(exx,0)-1)*DimDelta(exx,0))	// if xx at high edge, reduce tiny bit
			dX -= DimDelta(exx,0)*1e-9												// This section needed so Interp3d does not clip off the last plane
		endif																				// and ditto for dimensions 1 and 2
		if (hh+dY == DimOffset(exx,1) + (DimSize(exx,1)-1)*DimDelta(exx,1))
			dY -= DimDelta(exx,1)*1e-9
		endif
		if (ff+dZ == DimOffset(exx,2) + (DimSize(exx,2)-1)*DimDelta(exx,2))
			dZ -= DimDelta(exx,2)*1e-9
		endif
		if (numtype(dX+dY+dZ))									// cannot calculate
			alpha = NaN
			return NaN
		endif
		Variable exx0 = Interp3d(exx, xx,hh,ff)
		Variable eyy0 = Interp3d(eyy, xx,hh,ff)
		Variable ezz0 = Interp3d(ezz, xx,hh,ff)
		Variable exy0 = Interp3d(exy, xx,hh,ff)
		Variable exz0 = Interp3d(exz, xx,hh,ff)
		Variable eyz0 = Interp3d(eyz, xx,hh,ff)
		Variable eyx0=exy0, ezx0=exz0, ezy0=eyz0

//		dexx_dX =(Interp3d(exx, xx+dX,hh,ff) - exx0)/dX	// all three derivatives of the 9 components of epsilon
		dexy_dX =(Interp3d(exy, xx+dX,hh,ff) - exy0)/dX
		dexz_dX =(Interp3d(exz, xx+dX,hh,ff) - exz0)/dX
		deyy_dX =(Interp3d(eyy, xx+dX,hh,ff) - eyy0)/dX
		deyz_dX =(Interp3d(eyz, xx+dX,hh,ff) - eyz0)/dX
		dezz_dX =(Interp3d(ezz, xx+dX,hh,ff) - ezz0)/dX
		deyx_dX = dexy_dX									// epsilon is a symmetric tensor
		dezx_dX = dexz_dX
		dezy_dX = deyz_dX
		dexx_dY =(Interp3d(exx, xx,hh+dY,ff) - exx0)/dY
		dexy_dY =(Interp3d(exy, xx,hh+dY,ff) - exy0)/dY
		dexz_dY =(Interp3d(exz, xx,hh+dY,ff) - exz0)/dY
//		deyy_dY =(Interp3d(eyy, xx,hh+dY,ff) - eyy0)/dY
		deyz_dY =(Interp3d(eyz, xx,hh+dY,ff) - eyz0)/dY
		dezz_dY =(Interp3d(ezz, xx,hh+dY,ff) - ezz0)/dY
		deyx_dY = dexy_dY
		dezx_dY = dexz_dY
		dezy_dY = deyz_dY
		dexx_dZ =(Interp3d(exx, xx,hh,ff+dZ) - exx0)/dZ
		dexy_dZ =(Interp3d(exy, xx,hh,ff+dZ) - exy0)/dZ
		dexz_dZ =(Interp3d(exz, xx,hh,ff+dZ) - exz0)/dZ
		deyy_dZ =(Interp3d(eyy, xx,hh,ff+dZ) - eyy0)/dZ
		deyz_dZ =(Interp3d(eyz, xx,hh,ff+dZ) - eyz0)/dZ
//		dezz_dZ =(Interp3d(ezz, xx,hh,ff+dZ) - ezz0)/dZ
		deyx_dZ = dexy_dZ
		dezx_dZ = dexz_dZ
		dezy_dZ = deyz_dZ
	endif

	alpha[0][0] = kappa[0][0] - trace - (dezx_dy - deyx_dz)	//	alpha = kappa^T - Tr(kappa) - (del x epsilon), for diagonal terms
	alpha[1][1] = kappa[1][1] - trace - (dexy_dz - dezy_dx)
	alpha[2][2] = kappa[2][2] - trace - (deyz_dx - dexz_dy)

	alpha[0][1] = kappa[1][0] - (dezy_dy - deyy_dz)				//	alpha = kappa^T - (del x epsilon), for off diagonal terms
	alpha[0][2] = kappa[2][0] - (dezz_dy - deyz_dz)

	alpha[1][0] = kappa[0][1] - (dexx_dz - dezx_dx)
	alpha[1][2] = kappa[2][1] - (dexz_dz - dezz_dx)

	alpha[2][0] = kappa[0][2] - (deyx_dx - dexx_dy)
	alpha[2][1] = kappa[1][2] - (deyy_dx - dexy_dy)

	Variable GND_Density = NumVarOrDefault("GND_DislocationDensity",GND_DislocationDensity_Al)// convertion factor from dislocation tensor to dislocation density
	if (exists("GND_DislocationDensity")!=2)																// if not in current data folder, check root folder
		GND_Density = NumVarOrDefault("root:GND_DislocationDensity",GND_DislocationDensity_Al)
	endif
	Variable GND		// sum of the absolute values
	GND   = abs(alpha[0][0]) + abs(alpha[0][1]) + abs(alpha[0][2]) 
	GND += abs(alpha[1][0]) + abs(alpha[1][1]) + abs(alpha[1][2])
	GND += abs(alpha[2][0]) + abs(alpha[2][1]) + abs(alpha[2][2])
	GND *= GND_Density
	Note/K alpha, ReplaceNumberByKey("GND_Density",note(alpha),GND_Density,"=")
	SetScale d 0,0,"radian / µm",alpha
	return GND			// return the sum |alpha_ij|, proportional to GND
End


// this calculates the Lattice Curvature Tensor, it returns the GND assuming that all elastic strains (epsilon) are zero!
Function curvatureTensorFromRods3d(xx,hh,ff,R3dX,R3dH,R3dF,kappa)
	Variable xx,hh,ff					// desired xyz position
	Wave R3dX,R3dH,R3dF				// three components of the Roriques vectors in 3-d
	Wave kappa							// the lattice curvature tensor

	// in here I need to calculate the gradients at (x,h,f) and build the dislocation tensor
	//
	// Lattice Curvature Tensor = kappa^t = 
	//		dRx/dx			dRy/dx			dRz/dx		// NOTE do not forget the transpose,  kappa = Transpose(grad Theta)
	//		dRx/dy			dRy/dy			dRz/dy
	//		dRx/dz			dRy/dz			dRz/dz
	//
	// so get the gradients, and then fill kappa

	Variable dX=DimDelta(R3dX,0), dY=DimDelta(R3dX,1), dZ=DimDelta(R3dX,2)
	if (xx+dX == DimOffset(R3dX,0) + (DimSize(R3dX,0)-1)*DimDelta(R3dX,0))// if xx at high edge, reduce tiny bit
		dX -= DimDelta(R3dX,0)*1e-9													// This section needed so Interp3d does not clip off the last plane
	endif																					// and ditto for dimensions 1 and 2
	if (hh+dY == DimOffset(R3dX,1) + (DimSize(R3dX,1)-1)*DimDelta(R3dX,1))
		dY -= DimDelta(R3dX,1)*1e-9
	endif
	if (ff+dZ == DimOffset(R3dX,2) + (DimSize(R3dX,2)-1)*DimDelta(R3dX,2))
		dZ -= DimDelta(R3dX,2)*1e-9
	endif

	Variable Rx = Interp3d(R3dX, xx,hh,ff)
	Variable Ry = Interp3d(R3dH, xx,hh,ff)
	Variable Rz = Interp3d(R3dF, xx,hh,ff)
	kappa[0][0] =(Interp3d(R3dX, xx+dX,hh,ff) - Rx)/dX	// kappa = Transpose(grad Theta)
	kappa[1][0] =(Interp3d(R3dX, xx,hh+dY,ff) - Rx)/dY
	kappa[2][0] =(Interp3d(R3dX, xx,hh,ff+dZ) - Rx)/dZ
	kappa[0][1] =(Interp3d(R3dH, xx+dX,hh,ff) - Ry)/dX
	kappa[1][1] =(Interp3d(R3dH, xx,hh+dY,ff) - Ry)/dY
	kappa[2][1] =(Interp3d(R3dH, xx,hh,ff+dZ) - Ry)/dZ
	kappa[0][2] =(Interp3d(R3dF, xx+dX,hh,ff) - Rz)/dX
	kappa[1][2] =(Interp3d(R3dF, xx,hh+dY,ff) - Rz)/dY
	kappa[2][2] =(Interp3d(R3dF, xx,hh,ff+dZ) - Rz)/dZ

	Variable Rmag = sqrt(Rx^2+Ry^2+Rz^2)				// magnitude of Rodriques vector
	Variable Rod2radian = 2*atan(Rmag)/Rmag			// converts Rodriques length to angle (radian)
	kappa *= Rod2radian										// convert to radian/µm
	SetScale d 0,0,"radian / µm",kappa

	Variable Tr = MatrixTrace(kappa)
	Variable GND_Density = NumVarOrDefault("GND_DislocationDensity",GND_DislocationDensity_Al)	// convertion factor from dislocation tensor to dislocation density
	if (exists("GND_DislocationDensity")!=2)																	// if not in current data folder, check root folder
		GND_Density = NumVarOrDefault("root:GND_DislocationDensity",GND_DislocationDensity_Al)
	endif
	Variable GND		// sum of the absolute values
	GND   = abs(kappa[0][0]-Tr)	+ abs(kappa[0][1])		+ abs(kappa[0][2]) 
	GND += abs(kappa[1][0])		+ abs(kappa[1][1]-Tr)	+ abs(kappa[1][2])
	GND += abs(kappa[2][0])		+ abs(kappa[2][1])		+ abs(kappa[2][2]-Tr)
	GND *= GND_Density
	Note/K kappa, ReplaceNumberByKey("GND_Density",note(kappa),GND_Density,"=")
	return GND			// return the sum |kappa_ij - I*Tr(kappa)|, proportional to GND
End
//
Function GraphTensorcomponents(type,[toPrint])
	String type								// needs to be "alpha"
	Variable toPrint
	toPrint = ParamIsDefault(toPrint) ? 0 : toPrint

	String typeList = "alpha;theta;epsilon"
	if (WhichListItem(type,typeList)<0)
		Prompt type,"type",popup,typeList		
		DoPrompt "type",type
		if (V_flag)
			return 1
		endif
		printf "•GraphTensorcomponents(\"%s\")\r",type
	endif
	String wPrefix="",tName=""

	String size = SelectString(toPrint,"\\Zr175","\\Zr100")
	if (stringmatch(type,"alpha"))
		wPrefix = "SliceAlpha"
		tName=size+"\\F'Symbol'a\\F]0\\Zr050 \\M"+size+"\\B"
	elseif (stringmatch(type,"theta"))
		wPrefix = "SliceTheta"
		tName=size+"\\F'Symbol'q\\F]0\\B"
	elseif (stringmatch(type,"epsilon"))
		wPrefix = "SliceEpsilon"
		tName=size+"\\F'Symbol'e\\F]0\\B"
	else
		DoAlert 0,"do not know about type = '"+type+"'"
		return 1
	endif
	Variable epsilonIsZero=0
	if (stringmatch(type,"alpha"))			// special for alpha with epsilon=0
		Wave sw = SliceAlphaXX
		epsilonIsZero = (NumberByKey("epsilonIsZero",note(sw),"=")) ? 1 : 0
	endif
	String ab, wName=wPrefix+"xx"
	Wave ww = $(wName)
	if (!WaveExists(ww))
		DoAlert 0,"the wave "+wName+" does not exist here, either in wrong data folder, or need run ReCalcAllTensorComponents()"
		return 1
	endif
	String list = "xx;xy;xz;yx;yy;yz;zx;zy;zz"
	Variable zrange, i
	for (i=0,zrange=-Inf;i<9;i+=1)
		wName=wPrefix+StringFromList(i,list)
		Wave ww = $wName
		if (WaveExists(ww))
			zrange = max(zrange, ArrayOf3dOrients#PickSymmetricRange(ww,0.025))
		endif
	endfor
	Prompt zrange,"symmetric color range for all 9 of the graphs (give the positive value)"
	DoPrompt "color range",zrange
	if (V_flag || !(zrange>0))
		return 1
	endif
	print "using zrange=",zrange

	Variable left,top,right,bottom			// /W=(left,top,right,bottom )
	Display /W=(173,236,957,710)

	Variable aspect=((DimSize(ww,1)-1)*DimDelta(ww,1)) / ((DimSize(ww,0)-1)*DimDelta(ww,0))
	ModifyGraph height={Aspect,aspect}
	for (i=0;i<9;i+=1)
		ab = StringFromList(i,list)
		wName=wPrefix+ab
		DisplayBoxFor3x3fromAB(ab,left,top,right,bottom,toPrint=toPrint)
		Display/W=(left,top,right,bottom)/HOST=# 
		RenameWindow #,$ab
		AppendImage $wName
		ModifyImage $wName ctab= {-zrange,zrange,RedWhiteBlue,0}
		ModifyGraph margin=5
		ModifyGraph gfMult=200,gbRGB=(GRAY,GRAY,GRAY), tick=2, mirror=1, minor=1
//		ModifyGraph nticks(left)=3
		if (toPrint)
			// ModifyGraph nticks=6, sep=3, stLen=4, btLen=8
			ModifyGraph manTick={0,20,0,0},manMinor={9,5}
		endif

		if (stringmatch(type,"theta"))
			TextBox/C/N=text0/F=0/S=3/B=1/A=LB/X=9/Y=5 tName+ab[0]+","+ab[1]+"\\M"
		else
			TextBox/C/N=text0/F=0/S=3/B=1/A=LB/X=9/Y=5 tName+ab+"\\M"
		endif
		if (epsilonIsZero && stringmatch(type,"alpha") && stringmatch(ab,"xx"))
			AppendText/N=text0  size+"\\F'Symbol'e\\F]0\\M=0"
		endif
		Label bottom "\\E"
		Label Left "\\E"
		SetAxis/A/R left
		if (toPrint)
			ModifyGraph noLabel = 2
		else
			ModifyGraph noLabel(left) = stringmatch(ab[1],"x") ? 0 : 2
			ModifyGraph noLabel(bottom) = stringmatch(ab[0],"z") ? 0 : 2
		endif
		if (stringmatch(ab,"xz"))
			ColorScale/C/N=colorScale1/B=1/A=RC/X=0.00/Y=0.00 image=$wName, heightPct=60
			ColorScale/C/N=colorScale1 lblMargin=0, lowTrip=0.01, tickUnit=1, ZisZ=1
			ColorScale/C/N=colorScale1/F=0
//			ColorScale/C/N=colorScale1 lblMargin=0, lowTrip=0.001, tickUnit=1, ZisZ=1
//			AppendText WaveUnits($wName,-1)	// AppendText "radian / µm"
			AppendText ""
			ColorScale/C/N=colorScale1 fsize=1
		endif
		SetActiveSubwindow ##
	endfor

	if (!toPrint)
		SetDrawLayer UserFront
		SetDrawEnv fsize= 18
		DrawText 0.52,1.108, "Y (µm)"
		SetDrawEnv fsize= 18,textrot= 90
		DrawText 0.01530612244898,0.546413502109705, "Z (µm)"
		ModifyGraph margin(bottom)=50
	endif
	return 0
End
//
Static Function DisplayBoxFor3x3fromAB(ab,left,top,right,bottom,[toPrint])
	String ab
	Variable &left,&top,&right,&bottom			// /W=(left,top,right,bottom )
	Variable toPrint
	toPrint = ParamIsDefault(toPrint) ? 0 : toPrint
	ab = LowerStr(ab)							// ensure lower case

	Variable delta = toPrint ? 0.33 : 0.31
	Variable start = toPrint ? 0.01 : 0.07

	Variable i,j
	if (strlen(ab)==2)
		i = char2num(ab[1])-120			// 120 == 'x'
		j = char2num(ab[0])-120
		top = j*delta
		bottom = (j+1)*delta
	elseif (strlen(ab)==1)
		i = char2num(ab)-120	
		top = 0
//		bottom = toPrint ? 0.93 : 0.79
		bottom = toPrint ? 0.99 : 0.79
	endif
	left = start + i*delta
	right = start + (i+1)*delta
	return 0
End
//Static Function DisplayBoxFor3x3fromAB(ab,left,top,right,bottom)
//	String ab
//	Variable &left,&top,&right,&bottom			// /W=(left,top,right,bottom )
//	ab = LowerStr(ab)							// ensure lower case
//
//	Variable i,j
//	if (strlen(ab)==2)
//		i = char2num(ab[1])-120			// 120 == 'x'
//		j = char2num(ab[0])-120
//		top = j*0.31
//		bottom = (j+1)*0.31
//	elseif (strlen(ab)==1)
//		i = char2num(ab)-120	
//		top = 0
//		bottom = 0.79
//	endif
//	left = 0.07 + i*0.31
//	right = 0.07 + (i+1)*0.31
//	return 0
//End
//
Function ReCalcAllTensorComponents(type)
	String type								// needs to be "alpha or theta"

	String typeList = "alpha;theta;epsilon"
	if (WhichListItem(type,typeList)<0)
		Prompt type,"type",popup,typeList		
		DoPrompt "type",type
		if (V_flag)
			return 1
		endif
		printf "•ReCalcAllTensorComponents(\"%s\")\r",type
	endif

	String wPrefix="", prefix=""
	Variable transpose=0, epsilon=0
	if (stringmatch(type,"alpha"))
		wPrefix = "SliceAlpha"
		prefix = "alpha"
	elseif (stringmatch(type,"theta"))
		wPrefix = "SliceTheta"
		prefix = "kappa"
		transpose = 1
	elseif (stringmatch(type,"epsilon"))
		wPrefix = "SliceEpsilon"
		prefix = "e"
		epsilon = 1
	endif
	if (strlen(wPrefix)<1)
		DoAlert 0,"do not know about type = '"+type+"'"
		return 1
	endif

	Wave sliceW = ImageNameToWaveRef("GraphSlice", "sliceWave" )	// the 2-d wave displayed as an image and contours
	String fldr
	if (WaveExists(sliceW))
		fldr = GetWavesDataFolder(sliceW,1)
	else
		fldr = GetDataFolder(1)
	endif
	NVAR SliceValue=$(fldr+"SliceValue")
	SVAR SliceNormal=$(fldr+"SliceNormal")					// name of normal to sliceWave
	SVAR SliceAngleAxis=$(fldr+"SliceAngleAxis")				// name of desired rotation axis

	String noteStr=note(SliceWave),waveClass
	waveClass = AddListItem("TensorOrientationSlice", StringByKey("waveClass",noteStr,"="),",")
	noteStr = ReplaceStringByKey("waveClass",noteStr, waveClass,"=")

	String list = "xx;xy;xz;yx;yy;yz;zx;zy;zz"
	String ab, wName
	Variable i
	String SliceAngleAxisOld=SliceAngleAxis						// save current setting
	for (i=0;i<9;i+=1)
		ab = StringFromList(i,list)
		SliceAngleAxis = prefix+"("+ab+")"
		if (epsilon)
			if (char2num(ab[0]) > char2num(ab[1]))
				SliceAngleAxis = prefix+ab[1]+ab[0]
			else
				SliceAngleAxis = prefix+ab
			endif
		endif
		SetSliceProc("",SliceValue,"","SliceValue")				// force recalculation
		wName = wPrefix+UpperStr(ab)
		if (transpose)
			wName = wPrefix+UpperStr(ab[1]+ab[0])			// kappa is transpose of grad theta
		endif
		Duplicate/O SliceWave $wName
		Note/K $wName, noteStr
	endfor
	SliceAngleAxis = SliceAngleAxisOld							// restore original setting
	SetSliceProc("",SliceValue,"","SliceValue")					// force recalculation
End


Function Graph3AngleComponents(type,[toPrint])
	String type								// needs to be "theta"
	Variable toPrint
	toPrint = ParamIsDefault(toPrint) ? 0 : toPrint

	String typeList = "theta"
	if (WhichListItem(type,typeList)<0)
		Prompt type,"type",popup,typeList		
		DoPrompt "type",type
		if (V_flag)
			return 1
		endif
		printf "•Graph3AngleComponents(\"%s\")\r",type
	endif
	String wPrefix="",tName=""
	if (stringmatch(type,"theta"))
		wPrefix = "SliceTheta"
		if (toPrint)
			tName="\\Zr100\\F'Symbol'q\\F]0\\B"
		else
			tName="\\Zr175\\F'Symbol'q\\F]0\\B"
		endif
	else
		DoAlert 0,"do not know about type = '"+type+"'"
		return 1
	endif
	String ab, wName=wPrefix+"X"
	Wave ww = $(wName)
	if (!WaveExists(ww))
		DoAlert 0,"the wave "+wName+" does not exist here, either in wrong data folder, or need run ReCalc3AngleComponents()"
		return 1
	endif
	String list = "x;h;f"
	Variable zrange, i
	for (i=0,zrange=-Inf;i<3;i+=1)
		wName=wPrefix+StringFromList(i,list)
		Wave ww = $wName
		if (WaveExists(ww))
			zrange = max(zrange, PickSymmetricRange(ww,0.025))
		endif
	endfor
	Prompt zrange,"symmetric color range for all 9 of the graphs (give the positive value)"
	DoPrompt "color range",zrange
	if (V_flag || !(zrange>0))
		return 1
	endif
	print "using zrange=",zrange

	Variable left,top,right,bottom			// /W=(left,top,right,bottom )
	Display /W=(173,236,957,710)
	Variable aspect=((DimSize(ww,1)-1)*DimDelta(ww,1)) / ((DimSize(ww,0)-1)*DimDelta(ww,0))
	aspect /= 3
	ModifyGraph height={Aspect,aspect}
	for (i=0;i<3;i+=1)
		ab = StringFromList(i,list)
		wName=wPrefix+ab
		ab = SelectString(stringmatch(ab,"h"),ab,"y")
		ab = SelectString(stringmatch(ab,"f"),ab,"z")
		DisplayBoxFor3x3fromAB(ab,left,top,right,bottom,toPrint=toPrint)
		Display/W=(left,top,right,bottom)/HOST=# 
		RenameWindow #,$ab
		AppendImage $wName
		ModifyImage $wName ctab= {-zrange,zrange,RedWhiteBlue,0}
		ModifyGraph margin=5
		ModifyGraph gfMult=200,gbRGB=(GRAY,GRAY,GRAY), tick=2, mirror=1, minor=1
		if (toPrint)
			//	ModifyGraph nticks=6, sep=3, stLen=4, btLen=8
			ModifyGraph manTick={0,20,0,0},manMinor={9,5}
		endif
		TextBox/C/N=text0/F=0/S=3/B=1/A=LB/X=9/Y=5 tName+ab+"\\M"
		Label bottom "\\E"
		Label Left "\\E"
		SetAxis/A/R left
		if (toPrint)
			ModifyGraph noLabel = 2
		else
			ModifyGraph noLabel(left) = stringmatch(ab,"x") ? 0 : 2
		endif
		if (stringmatch(ab,"z"))
			ColorScale/C/N=colorScale1/B=1/A=RC/X=0.00/Y=0.00 image=$wName, heightPct=60
			ColorScale/C/N=colorScale1 lblMargin=0, lowTrip=0.01, tickUnit=1, ZisZ=1
			ColorScale/C/N=colorScale1/F=0
			AppendText ""
			ColorScale/C/N=colorScale1 fsize=1
		endif
		SetActiveSubwindow ##
	endfor

	if (!toPrint)
		SetDrawLayer UserFront
		SetDrawEnv fsize= 18,textxjust= 1,textyjust= 2
		DrawText 0.55,1.15,"Y (µm)"
		SetDrawEnv fsize= 18,textyjust= 1,textrot= 90
		DrawText 0,0.5, "Z (µm)"
		ModifyGraph margin(bottom)=50
	endif
	return 0
End
//
Function ReCalc3AngleComponents(type)
	String type								// needs to be "theta"

	String typeList = "theta"
	if (WhichListItem(type,typeList)<0)
		Prompt type,"type",popup,typeList		
		DoPrompt "type",type
		if (V_flag)
			return 1
		endif
		printf "•ReCalc3AngleComponents(\"%s\")\r",type
	endif

	String wPrefix=""
	if (stringmatch(type,"theta"))
		wPrefix = "SliceTheta"
	endif
	if (strlen(wPrefix)<1)
		DoAlert 0,"do not know about type = '"+type+"'"
		return 1
	endif

	Wave sliceW = ImageNameToWaveRef("GraphSlice", "sliceWave" )	// the 2-d wave displayed as an image and contours
	String fldr
	if (WaveExists(sliceW))
		fldr = GetWavesDataFolder(sliceW,1)
	else
		fldr = GetDataFolder(1)
	endif
	NVAR SliceValue=$(fldr+"SliceValue")
	SVAR SliceNormal=$(fldr+"SliceNormal")					// name of normal to sliceWave
	SVAR SliceAngleAxis=$(fldr+"SliceAngleAxis")				// name of desired rotation axis

	String noteStr=note(SliceWave),waveClass
	waveClass = AddListItem("3AngleOrientationSlice", StringByKey("waveClass",noteStr,"="),",")
	noteStr = ReplaceStringByKey("waveClass",noteStr, waveClass,"=")

	String list = "X;H;F"
	String ab, wName
	Variable i
	String SliceAngleAxisOld=SliceAngleAxis						// save current setting
	for (i=0;i<3;i+=1)
		ab = StringFromList(i,list)
		SliceAngleAxis = ab
		SetSliceProc("",SliceValue,"","SliceValue")				// force recalculation
		wName = wPrefix+UpperStr(ab)
		Duplicate/O SliceWave $wName
		Note/K $wName, noteStr
	endfor
	SliceAngleAxis = SliceAngleAxisOld							// restore original setting
	SetSliceProc("",SliceValue,"","SliceValue")					// force recalculation
End



Function Change_GND_DislocationDensity(GND_Density)		// used to set the GND dislocation denstiy conversion factor
	Variable GND_Density

	if (!(GND_Density>0))												// invalid input given
		GND_Density = NumVarOrDefault("GND_DislocationDensity",GND_DislocationDensity_Al)// current value
		if (exists("GND_DislocationDensity")!=2)																// if not in current data folder, check root folder
			GND_Density = NumVarOrDefault("root:GND_DislocationDensity",GND_DislocationDensity_Al)
		endif
		Prompt GND_Density,"GND dislocation density (converts rad/µm -> dislocations/cm2)"
		DoPrompt "GND conversion factor",GND_Density
		if (V_flag)
			return NaN
		endif
		if (!(GND_Density>0))
			DoAlert 0,"Invalid input, setting GND_Density to the Aluminum value, 77.92e10"
			GND_Density = GND_DislocationDensity_Al
		endif
	endif

	if (ItemsInList(GetRTStackInfo(0))<2)
		printf "setting the GND Dislocation Density converion factor to %g, (converts rad/µm -> dislocations/cm^2)\r",GND_Density
	endif
	Variable/G GND_DislocationDensity=GND_Density
	return GND_Density
End

// ************************************************************************************************************************
// ************************************************************************************************************************







// calculate the gradient of the total rotation angle at the point (xx,hh,ff) from the 3d-Rodriques vectors
Function GradtRotationAnglesFromRods3d(xx,hh,ff,R3dX,R3dH,R3dF,grad)
	Variable xx,hh,ff					// desired xyz position
	Wave R3dX,R3dH,R3dF				// three components of the Roriques vectors in 3-d
	Wave grad							// 3-vector that will receive the gradient(|R|)

	// in here I need to calculate the gradient at (x,h,f) or the total rotation angle (a scalar field)
	// the gradient is the 3-vector (d|R|/dX, d|R|/dY, d|R|/dZ)
	// for R is the Rodriques vector, then the magnitude of R = |R| = tan(angle/2), where angle is the total rotation angle
	//	angle = 2*atan(|R|)

	Variable dX=DimDelta(R3dX,0), dY=DimDelta(R3dX,1), dZ=DimDelta(R3dX,2)
	if (xx+dX == DimOffset(R3dX,0) + (DimSize(R3dX,0)-1)*DimDelta(R3dX,0))	// if xx at high edge, reduce tiny bit
		dX -= DimDelta(R3dX,0)*1e-9												// This section needed so Interp3d does not clip off the last plane
	endif																			// and ditto for dimensions 1 and 2
	if (hh+dY == DimOffset(R3dX,1) + (DimSize(R3dX,1)-1)*DimDelta(R3dX,1))
		dY -= DimDelta(R3dX,1)*1e-9
	endif
	if (ff+dZ == DimOffset(R3dX,2) + (DimSize(R3dX,2)-1)*DimDelta(R3dX,2))
		dZ -= DimDelta(R3dX,2)*1e-9
	endif
	Variable Rx,Ry,Rz									// individual components of Rodriques vector at some point
	Rx = Interp3d(R3dX, xx,hh,ff)						// Rx, Ry, Rz are the components of the Rodriques vector at a point
	Ry = Interp3d(R3dH, xx,hh,ff)						// (Rx,Ry,Rz) is the Rodriquies vector
	Rz = Interp3d(R3dF, xx,hh,ff)
	Variable Ro = 2*atan(sqrt(Rx*Rx+Ry*Ry+Rz*Rz))// total rotation angle at (xhf) (radian)

	Rx = Interp3d(R3dX, xx+dX,hh,ff)					// recompute R at (x+dX,h,f)
	Ry = Interp3d(R3dH, xx+dX,hh,ff)
	Rz = Interp3d(R3dF, xx+dX,hh,ff)
	Variable Rdx = 2*atan(sqrt(Rx*Rx+Ry*Ry+Rz*Rz))// total rotation angle at (x+dX,h,f) (radian)

	Rx = Interp3d(R3dX, xx,hh+dY,ff)					// recompute R at (x,h+dY,f)
	Ry = Interp3d(R3dH, xx,hh+dY,ff)
	Rz = Interp3d(R3dF, xx,hh+dY,ff)
	Variable Rdy = 2*atan(sqrt(Rx*Rx+Ry*Ry+Rz*Rz))// total rotation angle at (x,h+dY,f) (radian)

	Rx = Interp3d(R3dX, xx,hh,ff+dZ)					// recompute R at (x,h,f+dZ)
	Ry = Interp3d(R3dH, xx,hh,ff+dZ)
	Rz = Interp3d(R3dF, xx,hh,ff+dZ)
	Variable Rdz = 2*atan(sqrt(Rx*Rx+Ry*Ry+Rz*Rz))// total rotation angle at (x,h,f+dZ) (radian)

	grad[0] = (Rdx-Ro)/dX								// d|R| / dX
	grad[1] = (Rdy-Ro)/dY								// d|R| / dY
	grad[2] = (Rdz-Ro)/dZ								// d|R| / dZ
	SetScale d 0,0,"radian / µm",grad
	return norm(grad)
End






// special routine just for Ben, probably temporary that will cut off one side of a constant F surface
Function CutOffOneSide(cutOff)
	Variable cutOff			// 1=-X,   2=+X,   3=-H,   4=+H

	Wave sliceWave=sliceWave, sliceWaveAbs=sliceWaveAbs
	if (!WaveExists(sliceWave) || !WaveExists(sliceWave) || !WaveExists(sliceWave))
		DoAlert 0, "Cannot find the wave named 'sliceWave'"
		return 1
	endif
	String noteStr = note(sliceWave)
	if (!stringmatch(StringByKey("XaxisName",noteStr,"="),"X") || !stringmatch(StringByKey("YaxisName",noteStr,"="),"H"))
		DoAlert 0, "Can only slice off one side for a cut perpendicular to F"
		return 1
	endif
	if (cutOff!=limit(cutOff,1,4) || cutOff!=round(cutOff))
		Prompt cutOff, "side to cut off",popup,"-X;+X;-H;+H"
		DoPrompt "side to remove",cutOff
		if (V_flag)
			return 1
		endif
	endif

	Variable XisX = stringmatch(StringByKey("XaxisName",noteStr,"="),"X")	// true if X is first dimension
	Variable i,j				// indicies corresponding to the (0,0) in (X,H)
	i = -DimOffset(sliceWave,0)/DimDelta(sliceWave,0)
	j = -DimOffset(sliceWave,1)/DimDelta(sliceWave,1)
	Variable k,m
	k = -DimOffset(sliceWaveAbs,0)/DimDelta(sliceWaveAbs,0)
	m = -DimOffset(sliceWaveAbs,1)/DimDelta(sliceWaveAbs,1)
	Variable delta=0.5		// used to force the strip containing zero to be kept
	if (cutOff==1)			// remove -X
		i = floor(i-delta)
		j = floor(j-delta)
		k = floor(k-delta)
		m = floor(m-delta)
		if (XisX)
			sliceWave[0,i][] = NaN
			sliceWaveAbs[0,k][] = NaN
//			sliceWaveLog[0,k][] = NaN
		else
			sliceWave[][0,m] = NaN
			sliceWaveAbs[][0,m] = NaN
//			sliceWaveLog[][0,m] = NaN
		endif
	elseif (cutOff==2)		// remove +X
		i = ceil(i+delta)
		j = ceil(j+delta)
		k = ceil(k+delta)
		m = ceil(m+delta)
		if (XisX)
			sliceWave[i,Inf][] = NaN
			sliceWaveAbs[k,Inf][] = NaN
//			sliceWaveLog[k,Inf][] = NaN
		else
			sliceWave[][j,Inf] = NaN
			sliceWaveAbs[][m,Inf] = NaN
//			sliceWaveLog[][m,Inf] = NaN
		endif
	elseif (cutOff==3)		// remove -H
		i = floor(i-delta)
		j = floor(j-delta)
		k = floor(k-delta)
		m = floor(m-delta)
		if (!XisX)
			sliceWave[0,i][] = NaN
			sliceWaveAbs[0,k][] = NaN
//			sliceWaveLog[0,k][] = NaN
		else
			sliceWave[][0,j] = NaN
			sliceWaveAbs[][0,m] = NaN
//			sliceWaveLog[][0,m] = NaN
		endif
	elseif (cutOff==4)		// remove +H"
		i = ceil(i+delta)
		j = ceil(j+delta)
		k = ceil(k+delta)
		m = ceil(m+delta)
		if (!XisX)
			sliceWave[i,Inf][] = NaN
			sliceWaveAbs[k,Inf][] = NaN
//			sliceWaveLog[k,Inf][] = NaN
		else
			sliceWave[][j,Inf] = NaN
			sliceWaveAbs[][m,Inf] = NaN
//			sliceWaveLog[][m,Inf] = NaN
		endif
	else
		Abort "Invalid number for cutOff"
	endif
End


Function/S RodriquesVector2str(x,h,f)	// print out one Rodriques vector at (X,H,F)
	Variable x,h,f
	if (numtype(x+h+f))
		x = numtype(x) ? 0 : x
		h = numtype(h) ? 0 : h
		f = numtype(f) ? 0 : f
		Prompt x, "X position in sample (µm)"
		Prompt h, "H position in sample (µm)"
		Prompt f, "F position in sample (µm)"
		DoPrompt "position in sample", x,h,f
		if (V_flag)
			return ""
		endif
	endif
	Variable printIt = ItemsInList(GetRTStackInfo(0))<2 || WhichListItem("testing3dArraysPopMenuProc", GetRTStackInfo(0))>=0
	if (numtype(x+h+f))
		if (printIt)
			printf "Invalid inputs XHF={%.3g, %.3g, %.3g} µm\r",x,h,f
		endif
		return ""
	endif

	String str
	Wave rhat=$MakeUnique3Vector($"")		// Rodriques vector
	Wave R3dX=R3dX, R3dH=R3dH, R3dF=R3dF
	rhat[0] = Interp3d(R3dX, x,h,f)		// Rodriques vector at XHF
	rhat[1] = Interp3d(R3dH, x,h,f)
	rhat[2] = Interp3d(R3dF, x,h,f)
	rhat = (rhat[p]==0) ? 0 : rhat[p]
	Variable angle = 2*atan(normalize(rhat))*180/PI
	sprintf str, "at {%.3g, %.3g, %.3g}µm,   Rodriques^={%.3g, %.3g, %.3g},  rotation angle=%.3g°",x,h,f,rhat[0],rhat[1],rhat[2],angle
	if (printIt)
		print str
	endif
	KillWaves/Z rhat
	return str
End


Function TableOfRawData()		// put up a table of the raw data
	if (strlen(WinList("RawDataTable", "", "WIN:2")))
		DoWindow /F RawDataTable
		return 0
	endif
	Edit/W=(5,44,547,708)/K=1 XX,HH,FF,RX,RH,RF
	Wave irefMat=irefMat
	if (WaveExists(irefMat))
		AppendToTable irefMat
	endif
	DoWindow /C RawDataTable
	ModifyTable width(Point)=32,format(XX)=3,sigDigits(XX)=4,format(HH)=3,sigDigits(HH)=4
	ModifyTable format(FF)=3,sigDigits(FF)=4,format(RX)=3,digits(RX)=7,sigDigits(RX)=3
	ModifyTable format(RH)=3,digits(RH)=7,sigDigits(RH)=3,format(RF)=3,digits(RF)=7
	ModifyTable sigDigits(RF)=3
	if (WaveExists(irefMat))
		ModifyTable width(irefMat)=45
	endif
	ModifyTable width(XX)=60, width(HH)=60, width(FF)=60
End


Function TableOfIndexingList(wlist)		// put up a table of the raw data
	Wave wlist
	if (strlen(WinList("IndexingListTable", "", "WIN:2")))
		DoWindow /F IndexingListTable
		return 0
	endif

	String list = WaveListClass("indexationResultList","*","TEXT:1")
	if (!WaveExists(wlist) && ItemsInLIst(list)>0)
		String IndexingListName=StringFromLIst(0,list)
		if ( ItemsInLIst(list)>1)
			Prompt IndexingListName,"Text Wave with Indexations",popup,list
			DoPrompt "pick IndexingList",IndexingListName
			if (V_flag)
				return 1
			endif
		endif
		Wave wlist = $IndexingListName
	endif
	if (!WaveExists(wlist))
		Wave wlist = indexingList
	endif
	if (!WaveExists(wlist))
		Wave wlist = ::IndexingList
	endif
	if (!WaveExists(wlist))
		if (ItemsInList(GetRTStackInfo(0))<2 || stringmatch(GetRTStackInfo(0),"*ButtonProc;findClosestPoint;"))
			DoAlert 0,"cannot find indexingList"
		endif
		return 1
	endif

	Edit/W=(5,44,1228,684)/K=1 wlist
	DoWindow /C IndexingListTable
	ModifyTable format(Point)=1,alignment(wlist)=0,width(wlist)=1138
	return 0
End



Function MakeGraphOfSlice()				// create the standard surface cut through the volume
	Wave R3dX=R3dX, R3dH=R3dH, R3dF=R3dF
	Wave R2dX=R2dX, R2dH=R2dH, R2dF=R2dF
	if (!WaveExists(R3dX) || !WaveExists(R3dH) || !WaveExists(R3dF))
		if (!WaveExists(R2dX) || !WaveExists(R2dH) || !WaveExists(R2dF))
			DoAlert 0, "MakeGraphOfSlice(), cannot find R3dX, R3dH, & R3dF  or  R2dX, R2dH, & R2dF"
			return 1
		endif
	endif
	if (strlen(WinList("GraphSlice", "", "WIN:1")))
		DoWindow /F GraphSlice
		return 1
	endif
	String fldr= GetDataFolder(1)			// this is the folder where everything is created and stored
	if (exists("sliceWave")!=1)			// sliceWave does not exist, yet
		Make/N=(10,10)/O sliceWave
	else
		Wave sliceWave=sliceWave
	endif
	String noteStr = note(sliceWave)
	Display/W=(466,47,1016,660)/K=1 
	DoWindow /C GraphSlice
	AppendMatrixContour sliceWave

	Variable/G $( fldr+"SliceValue")		// ensure existance of SliceValue
	Variable/G $(fldr+"epsilonIsZero")		// ensure existance of epsilonIsZero

	String/G $(fldr+"SliceAngleAxis")		// ensure existance of SliceAngleAxis, the name of the rotation axis
	SVAR SliceAngleAxis=$(fldr+"SliceAngleAxis")
	SliceAngleAxis = SelectString(strlen(SliceAngleAxis),"X",SliceAngleAxis)

	String/G $( fldr+"SliceNormal")		// ensure existance of SliceNormal
	SVAR SliceNormal=$( fldr+"SliceNormal")
	SliceNormal = SelectString(strlen(SliceNormal),"H",SliceNormal)

	Variable/G $(fldr+"phiSlice")			// ensure existance of phiSlice
	NVAR phi=$(fldr+"phiSlice")

	Variable/G $(fldr+"CutResolution")		// ensure existance of resolution
	NVAR CutResolution=$(fldr+"CutResolution")
	CutResolution = (CutResolution>0.0) ? CutResolution : 1.0

	ModifyContour sliceWave autoLevels={*,*,5}, rgbLines=(0,0,0)	// ModifyContour sliceWave autoLevels={0,0.001,11}, rgbLines=(0,0,0)
	ModifyContour sliceWave labelBkg=1
	AppendImage sliceWave
	ModifyImage sliceWave ctab= {*,*,RedWhiteBlue,0}	// red is negative, white=0, blue is postive
	ModifyGraph gfMult=130
	ModifyGraph gbRGB=(GRAY,GRAY,GRAY)				// set graph background to gray so that NaN will show through
	ModifyGraph mirror=1,minor=1
	ModifyGraph/Z lSize('sliceWave=0')=2
	ModifyGraph zero=2
	Label left StringByKey("YaxisName",noteStr,"=")+" (\\U)"
	Label bottom StringByKey("XaxisName",noteStr,"=")+" (\\U)"
//	SetAxis/A/R left
	SetAxis/A left
	SetAxis/A bottom
	Cursor/P/I A sliceWave 1,1
	ShowInfo
	TextBox/N=titleText/F=0/S=3/A=LT/X=3.06/Y=2.06 "lattice tilt (radian)\r"+GetWavesDataFolder(sliceWave,2)

	ColorScale/N=colorScale1/A=RC/X=-4.75/Y=0 image=sliceWave, heightPct=60
	ColorScale/C/N=colorScale1 lowTrip=0.001,lblMargin=0,ZisZ=1
	ColorScale/C/N=colorScale1 "radian"
	ColorScale/C/N=colorScale1 tickUnit=1
	CutSlicePopUpProc("",0,"")					// force an update using current values
	DoUpdate
	MakePanelSlicer(sliceWave)					// put up the panel with controls
End

Function MakeGraphOfSliceWaveRGB()			// create the standard surface cut through the volume
	Wave sliceWaveRGB = sliceWaveRGB
	if (!WaveExists(sliceWaveRGB))
		DoAlert 0, "MakeGraphOfSliceWaveRGB(), cannot find 'sliceWaveRGB'"
		return 1
	endif
	String gName = StringFromList(0,FindGraphsWithWave(sliceWaveRGB))
	if (strlen(gName))
		DoWindow/F $gName
		return 1
	endif
	Display/K=1/W=(35,44,780,405)
	AppendImage sliceWaveRGB
	SurfacePlotStyle("")
	return 0
End

Function MakeGraphOfSliceOrients()			// create the standard surface cut through the volume
	Wave sliceOrients = sliceOrients
	if (!WaveExists(sliceOrients))
		DoAlert 0, "MakeGraphOfSliceOrients(), cannot find 'sliceOrients'"
		return 1
	endif
	if (strlen(WinList("GraphSliceOrients", "", "WIN:1")))
		DoWindow /F GraphSliceOrients
		return 1
	endif

	String noteStr = note(sliceOrients)
	Display/K=1/W=(5,44,700,589)
	DoWindow /C GraphSliceOrients
	AppendImage sliceOrients

	String fldr= GetDataFolder(1)			// this is the folder where everything is created and stored
	Variable/G $(fldr+"SliceValue")		// ensure existance of SliceValue
	String/G $(fldr+"SliceAngleAxis")		// ensure existance of SliceAngleAxis, the name of the rotation axis
	SVAR SliceAngleAxis=$(fldr+"SliceAngleAxis")
	SliceAngleAxis = SelectString(strlen(SliceAngleAxis),"X",SliceAngleAxis)
	String/G $( fldr+"SliceNormal")		// ensure existance of SliceNormal
	SVAR SliceNormal=$( fldr+"SliceNormal")
	SliceNormal = SelectString(strlen(SliceNormal),"H",SliceNormal)
	Variable/G $(fldr+"phiSlice")			// ensure existance of phiSlice
	NVAR phi=$(fldr+"phiSlice")
	Variable/G $(fldr+"CutResolution")		// ensure existance of resolution
	NVAR CutResolution=$(fldr+"CutResolution")
	CutResolution = (CutResolution>0.0) ? CutResolution : 1.0
	Variable/G $(fldr+"epsilonIsZero")		// ensure existance of epsilonIsZero

	ModifyGraph gfMult=130, mirror=1,minor=1, zero=2
	Label left StringByKey("YaxisName",noteStr,"=")+" (\\U)"
	Label bottom StringByKey("XaxisName",noteStr,"=")+" (\\U)"
	SetAxis/A left
	SetAxis/A bottom
	TextBox/N=titleText/F=0/S=3/A=LT/X=3.06/Y=2.06 "Orientation\r"+GetWavesDataFolder(sliceOrients,2)
	DoUpdate
	SetAspectToSquarePixels("GraphSliceOrients")
End

Function showCubicTriangleColors(Np,saturate)
	Variable Np					// dimension of color triangle image
	Variable saturate			// make all colors saturated (max rgb always 65535)

	if (!(Np>10) || numtype(saturate))
		Np = !(Np>10) ? 400 : Np
		saturate = numtype(saturate) ? 1 : saturate
		saturate = saturate ? 1 : 2
		Prompt Np, "dimension of color triangle image"
		Prompt saturate, "make all colors saturated (max rgb always 65535)",popup,"saturate;un-saturated (no white spot)"
		DoPrompt "dimensions",Np,saturate
		saturate = saturate==1
		if (V_flag)
			return 1
		endif
	endif
	if (!(Np>10) || numtype(saturate))
		return 1
	endif

	Make/N=3/O/D xhat,yhat,zhat, hkl
	zhat = {1,2,3}
	FindPerpVector(zhat,xhat)
	Cross zhat, xhat
	Wave W_Cross=W_Cross
	yhat = W_Cross
	normalize(xhat)
	normalize(zhat)
	normalize(yhat)

	Wave n001=$MakeUnique3Vector($"")		// unit vector towards (001)
	Wave n011=$MakeUnique3Vector($"")		// unit vector towards (011)
	Wave n111=$MakeUnique3Vector($"")		// unit vector towards (111)
	Wave n001p=$MakeUnique3Vector($"")	// unit vector of n011 x n111, defines edge oppotiset 001
	Wave n011p=$MakeUnique3Vector($"")	// unit vector of n111 x n001, defines edge oppotiset 011
	Wave n111p=$MakeUnique3Vector($"")	// unit vector of n001 x n011, defines edge oppotiset 111
	n001 = {0,0,1}
	n011 = {0,1,1}
	n111 = {1,1,1}
	normalize(n011)
	normalize(n111)
	// change vectors to be in the local system
	hkl=n111
	n111 = {MatrixDot(hkl,yhat),-MatrixDot(hkl,xhat),MatrixDot(hkl,zhat)}
	hkl=n011
	n011 = {MatrixDot(hkl,yhat),-MatrixDot(hkl,xhat),MatrixDot(hkl,zhat)}
	hkl=n001
	n001 = {MatrixDot(hkl,yhat),-MatrixDot(hkl,xhat),MatrixDot(hkl,zhat)}
	Cross n011,n111	;	n001p = W_Cross	;	normalize(n001p)
	Cross n111,n001	;	n011p = W_Cross	;	normalize(n011p)
	Cross n001,n011	;	n111p = W_Cross	;	normalize(n111p)

	Variable x001,y001, x011,y011, x111,y111
	Variable xmax,xmin,ymax,ymin
	CubicHKL2stereo(0,0,1,x001,y001)
	CubicHKL2stereo(0,1,1,x011,y011)
	CubicHKL2stereo(1,1,1,x111,y111)
	xmin = min(min(x001,x011),x111)		// and get the needed range of box
	xmax = max(max(x001,x011),x111)
	ymin = min(min(y001,y011),y111)
	ymax = max(max(y001,y011),y111)

	Make/N=(Np,Np,3)/O/W/U CubicColors		// make rgb image wave
	SetScale/I x xmin,xmax,"", CubicColors
	SetScale/I y ymin,ymax,"", CubicColors
	CubicColors = 0

	Make/N=3/D/O oneColor,edgeVec, dn
	Variable dx=DimDelta(CubicColors,0), dy=DimDelta(CubicColors,1)
	Variable xx,yy
	Variable i,j, maxc, t
	for (j=0;j<Np;j+=1)
		yy = j*dy + ymin
		for (i=0;i<Np;i+=1)
			xx = i*dx + xmin
			CubicStereo2HKL(xx,yy,hkl)			// hkl is returned normalized
			if (MatrixDot(hkl,n001p)<=0 && MatrixDot(hkl,n011p)<=0 && MatrixDot(hkl,n111p)<=0)
				dn = hkl - n001
				t = -MatrixDot(n001,n001p)/MatrixDot(dn,n001p)
				edgeVec = n001 + t*dn										// vector at end of line from n001 through hkl
				normalize(edgeVec)
				oneColor[0] = 1 - (1-MatrixDot(n001,hkl))/(1-MatrixDot(n001,edgeVec))

				dn = hkl - n011
				t = -MatrixDot(n011,n011p)/MatrixDot(dn,n011p)
				edgeVec = n011 + t*dn										// vector at end of line from n011 through hkl
				normalize(edgeVec)
				oneColor[1] = 1 - (1-MatrixDot(n011,hkl))/(1-MatrixDot(n011,edgeVec))

				dn = hkl - n111
				t = -MatrixDot(n111,n111p)/MatrixDot(dn,n111p)
				edgeVec = n111 + t*dn										// vector at end of line from n111 through hkl
				normalize(edgeVec)
				oneColor[2] = 1 - (1-MatrixDot(n111,hkl))/(1-MatrixDot(n111,edgeVec))

				maxc = saturate ? max(max(oneColor[0],oneColor[1]),oneColor[2]) : 1	// max color at this pixel
				oneColor = limit(oneColor*65535/maxc,0,65535)			// also limit to allowed 16 bit range
				CubicColors[i][j][] = oneColor[r]
			endif
		endfor
	endfor

	if (strlen(FindGraphsWithWave(CubicColors))<1)		// make the plot
		// Display/K=1/W=(1004,56,1424,322)
		Display/K=1/W=(650,57,1151,372)
		AppendImage CubicColors
		ModifyGraph tick=3,noLabel=2
		ModifyGraph axOffset(left)=-5,axOffset(bottom)=-1.41667
		TextBox/N=text001/F=0/S=3/A=LC/X=2.89/Y=7.34 "\\Zr125001"
		TextBox/N=text011/F=0/S=3/X=11.56/Y=2.75 "\\Zr125011"
		TextBox/N=text111/F=0/S=3/A=RB/X=17.92/Y=1.83 "\\Zr125111"
		DoUpdate
		SetAspectToSquarePixels("")
	endif
	if (0)
		SetDrawLayer /K UserFront				// always clear first
		SetDrawLayer UserFront
		SetDrawEnv xcoord= bottom,ycoord= left
		DrawLine x001,y001,x011,y011		// clear above this line
		SetDrawEnv xcoord= bottom,ycoord= left
		DrawLine x011,y011,x111,y111		// clear to right of this line
		SetDrawEnv xcoord= bottom,ycoord= left
		DrawLine x111,y111,x001,y001		// clear below this line
	endif
	KillWaves/Z n001,n011,n111, n001p,n011p,n111p, hkl,zhat,yhat,xhat,W_Cross, oneColor, edgeVec, dn
	return 0
End
//
//Static Function InternalToTriangle(xx,yy,x001,y001, x011,y011, x111,y111)
//	Variable xx,yy
//	Variable x001,y001, x011,y011, x111,y111
//
//	Variable m,b,h
//	m = (y011-y001)/(x011-x001)	// clear above this line
//	b = y011 - m*x011
//	h = m*xx+b
//	if (yy>h)
//		return 0
//	endif
//
//	m = (y111-y011)/(x111-x011)	// clear to right of this line
//	b = y111 - m*x111
//	h =(yy-b)/m
//	if (xx>h)
//		return 0
//	endif
//
//	m = (y001-y111)/(x001-x111)	// clear below this line
//	b = y001 - m*x001
//	h = m*xx+b
//	if (yy<h)
//		return 0
//	endif
//	return 1
//End
//
Static Function CubicStereo2HKL(xx,yy,hkl)
	Variable xx, yy
	Wave hkl

	Wave xhat=xhat, yhat=yhat, zhat=zhat, hkl=hkl
	Variable r, r2
	r2 = xx*xx + yy*yy
	r = sqrt(r2)
	Variable s = 4*r/(4+r2)			// sin(bet)
	Variable t = (r - s)/r
	hkl[0] = xx*(t-1)					// this is now a point on the sphere
	hkl[1] = yy*(t-1)					//  origin of sphere is (0,0,1)
	hkl[2] = 1- 2*t 					// shift to using origin of sphere
	return 0
End
//
Static Function CubicHKL2stereo(h,k,l,xx,yy)
	Variable h,k,l
	Variable &xx, &yy
	Wave xhat=xhat, yhat=yhat, zhat=zhat, hkl=hkl
	Variable xp,yp,zp,r, bet, alpha
	hkl = {h,k,l}
	normalize(hkl)
	bet = acos(min(MatrixDot(zhat,hkl),1))	// angle from pole to intersection with sphere (radians)
	alpha = atan2(MatrixDot(yhat,hkl),MatrixDot(xhat,hkl))+PI/2
	r = sin(bet)
	xp = r*cos(alpha)							// coords of point of sphere containing projection line
	yp = r*sin(alpha)							// the projection line starts at (0,0,2)
	zp = 1-cos(bet)
	xx = -2/(zp-2) * xp
	yy = -2/(zp-2) * yp
End
//
Static Function FindPerpVector(ref,perp)		// find a vector perpendicular to ref
	Wave ref									// reference vector, find a vec perp to this
	Wave perp									// resulting perpendicular vector, not normalized
	Make/N=3/O/D FindPerpVector_hat
	Wave hat=FindPerpVector_hat

	hat = abs(ref)
	WaveStats/M=1/Q hat
	hat = 0
	hat[V_maxLoc] = 1
	Variable dot = MatrixDot(ref,hat)
	perp = ref - dot*hat
	Cross ref, perp
	Wave W_cross=W_cross
	perp = W_cross
	KillWaves/Z FindPerpVector_hat, W_cross
	return 0
End


Function MakeGraph3dSurface()					// make a 3-d surface to go with cut
	if(exists("NewGizmo")!=4)								// Do nothing if the Gizmo XOP is not available.
		DoAlert 0, "Gizmo XOP must be installed"
		return 1
	endif
	if (exists("sliceWaveAbs")!=1 || exists("sliceWaveRGBA")!=1)
		DoAlert 0, "sliceWaveAbs or  sliceWaveRGBA does not exists, try MakeGraphOfSlice() first"
		return 1
	endif
	Wave sliceWaveAbs=sliceWaveAbs, sliceWaveRGBA=sliceWaveRGBA
	String str = WaveUnits(sliceWaveAbs,0)
	String XaxisName = StringByKey("XaxisName",note(sliceWaveAbs),"=")+SelectString(strlen(str),"","  ("+str+")")
	str = WaveUnits(sliceWaveAbs,1)
	String YaxisName = StringByKey("YaxisName",note(sliceWaveAbs),"=")+SelectString(strlen(str),"","  ("+str+")")

	Execute "NewGizmo/N=Gizmo0/T=\"GizmoSurfaceOrients\" /W=(448,185,1065,753)"
	Execute "ModifyGizmo startRecMacro"
	Execute "AppendToGizmo attribute blendFunc={770,771},name=blendFunc0"
	Execute "AppendToGizmo Surface="+ GetWavesDataFolder(sliceWaveAbs,2)+",name=surface0"
	Execute "ModifyGizmo ModifyObject=surface0 property={ surfaceColorType,3}"
	Execute "ModifyGizmo ModifyObject=surface0 property={ lineColorType,1}"
	Execute "ModifyGizmo ModifyObject=surface0 property={ lineWidthType,1}"
	Execute "ModifyGizmo ModifyObject=surface0 property={ fillMode,3}"
	Execute "ModifyGizmo ModifyObject=surface0 property={ srcMode,0}"
	Execute "ModifyGizmo ModifyObject=surface0 property={ lineColor,0,0,0,0.8}"
	Execute "ModifyGizmo ModifyObject=surface0 property={ surfaceColorWave, "+GetWavesDataFolder(sliceWaveRGBA,2)+"}"
	Execute "AppendToGizmo Axes=boxAxes,name=axes0"
	Execute "ModifyGizmo ModifyObject=axes0,property={0,axisRange,-1,-1,-1,1,-1,-1}"
	Execute "ModifyGizmo ModifyObject=axes0,property={1,axisRange,-1,-1,-1,-1,1,-1}"
	Execute "ModifyGizmo ModifyObject=axes0,property={2,axisRange,-1,-1,-1,-1,-1,1}"
	Execute "ModifyGizmo ModifyObject=axes0,property={3,axisRange,-1,1,-1,-1,1,1}"
	Execute "ModifyGizmo ModifyObject=axes0,property={4,axisRange,1,1,-1,1,1,1}"
	Execute "ModifyGizmo ModifyObject=axes0,property={5,axisRange,1,-1,-1,1,-1,1}"
	Execute "ModifyGizmo ModifyObject=axes0,property={6,axisRange,-1,-1,1,-1,1,1}"
	Execute "ModifyGizmo ModifyObject=axes0,property={7,axisRange,1,-1,1,1,1,1}"
	Execute "ModifyGizmo ModifyObject=axes0,property={8,axisRange,1,-1,-1,1,1,-1}"
	Execute "ModifyGizmo ModifyObject=axes0,property={9,axisRange,-1,1,-1,1,1,-1}"
	Execute "ModifyGizmo ModifyObject=axes0,property={10,axisRange,-1,1,1,1,1,1}"
	Execute "ModifyGizmo ModifyObject=axes0,property={11,axisRange,-1,-1,1,1,-1,1}"
	Execute "ModifyGizmo ModifyObject=axes0,property={-1,axisScalingMode,1}"
	Execute "ModifyGizmo ModifyObject=axes0,property={-1,axisColor,0,0,0,1}"

	Execute "ModifyGizmo ModifyObject=axes0,property={9,ticks,3}"
	Execute "ModifyGizmo ModifyObject=axes0,property={8,ticks,3}"
//	Execute "ModifyGizmo ModifyObject=axes0,property={9,ticks,3}"
	Execute "ModifyGizmo ModifyObject=axes0,property={5,ticks,3}"
	Execute "ModifyGizmo ModifyObject=axes0,property={9,axisLabel,1}"
	Execute "ModifyGizmo ModifyObject=axes0,property={8,axisLabel,1}"
//	Execute "ModifyGizmo ModifyObject=axes0,property={0,axisLabel,1}"
	if (strlen(XaxisName))
		Execute "ModifyGizmo ModifyObject=axes0,property={9,axisLabelText,\""+XaxisName+"\"}"
		Execute "ModifyGizmo ModifyObject=axes0,property={9,axisLabelCenter,-0.3}"
		Execute "ModifyGizmo ModifyObject=axes0,property={9,axisLabelDistance,0.1}"
		Execute "ModifyGizmo ModifyObject=axes0,property={9,axisLabelScale,0.5}"
		Execute "ModifyGizmo ModifyObject=axes0,property={9,axisLabelRGBA,0,0,0,1}"
	endif
	if (strlen(YaxisName))
		Execute "ModifyGizmo ModifyObject=axes0,property={8,axisLabelText,\""+YaxisName+"\"}"
		Execute "ModifyGizmo ModifyObject=axes0,property={8,axisLabelCenter,-0.5}"
		Execute "ModifyGizmo ModifyObject=axes0,property={8,axisLabelDistance,0.1}"
		Execute "ModifyGizmo ModifyObject=axes0,property={8,axisLabelScale,0.5}"
		Execute "ModifyGizmo ModifyObject=axes0,property={8,axisLabelRGBA,0,0,0,1}"
	endif
	Execute "ModifyGizmo setDisplayList=0, attribute=blendFunc0"
	Execute "ModifyGizmo setDisplayList=1, object=surface0"
	Execute "ModifyGizmo setDisplayList=2, object=axes0"
	Execute "ModifyGizmo SETQUATERNION={-0.102663,-0.598222,-0.788212,-0.101455}"
	Execute "ModifyGizmo autoscaling=1"
	Execute "ModifyGizmo currentGroupObject=\"\""
	Execute "ModifyGizmo compile"

//	Execute "ModifyGizmo showInfo"
//	Execute "ModifyGizmo infoWindow={2085,657,2558,920}"
	Execute "ModifyGizmo bringToFront"
	Execute "ModifyGizmo endRecMacro"
End
//Function MakeGraph3dSurface()					// make a 3-d surface to go with cut
//	if(exists("NewGizmo")!=4)								// Do nothing if the Gizmo XOP is not available.
//		DoAlert 0, "Gizmo XOP must be installed"
//		return 1
//	endif
//	if (exists("sliceWaveAbs")!=1 || exists("sliceWaveRGBA")!=1)
//		DoAlert 0, "sliceWaveAbs or  sliceWaveRGBA does not exists, try MakeGraphOfSlice() first"
//		return 1
//	endif
//	Wave sliceWaveAbs=sliceWaveAbs, sliceWaveRGBA=sliceWaveRGBA
//	String str = WaveUnits(sliceWaveAbs,0)
//	String XaxisName = StringByKey("XaxisName",note(sliceWaveAbs),"=")+SelectString(strlen(str),"","  ("+str+")")
//	str = WaveUnits(sliceWaveAbs,1)
//	String YaxisName = StringByKey("YaxisName",note(sliceWaveAbs),"=")+SelectString(strlen(str),"","  ("+str+")")
//
//	Execute "NewGizmo/N=Gizmo0/T=\"GizmoSurfaceOrients\" /W=(448,185,1065,753)"
//	Execute "ModifyGizmo startRecMacro"
//	Execute "AppendToGizmo Surface="+ GetWavesDataFolder(sliceWaveAbs,2)+",name=surface0"
//	Execute "ModifyGizmo ModifyObject=surface0 property={ surfaceColorType,3}"
//	Execute "ModifyGizmo ModifyObject=surface0 property={ lineColorType,1}"
//	Execute "ModifyGizmo ModifyObject=surface0 property={ lineWidthType,1}"
//	Execute "ModifyGizmo ModifyObject=surface0 property={ fillMode,3}"
//	Execute "ModifyGizmo ModifyObject=surface0 property={ srcMode,0}"
//	Execute "ModifyGizmo ModifyObject=surface0 property={ lineColor,0,0,0,1}"
//	Execute "ModifyGizmo ModifyObject=surface0 property={ surfaceColorWave, "+GetWavesDataFolder(sliceWaveRGBA,2)+"}"
//	Execute "AppendToGizmo Axes=boxAxes,name=axes0"
//	Execute "ModifyGizmo ModifyObject=axes0,property={0,axisRange,-1,-1,-1,1,-1,-1}"
//	Execute "ModifyGizmo ModifyObject=axes0,property={1,axisRange,-1,-1,-1,-1,1,-1}"
//	Execute "ModifyGizmo ModifyObject=axes0,property={2,axisRange,-1,-1,-1,-1,-1,1}"
//	Execute "ModifyGizmo ModifyObject=axes0,property={3,axisRange,-1,1,-1,-1,1,1}"
//	Execute "ModifyGizmo ModifyObject=axes0,property={4,axisRange,1,1,-1,1,1,1}"
//	Execute "ModifyGizmo ModifyObject=axes0,property={5,axisRange,1,-1,-1,1,-1,1}"
//	Execute "ModifyGizmo ModifyObject=axes0,property={6,axisRange,-1,-1,1,-1,1,1}"
//	Execute "ModifyGizmo ModifyObject=axes0,property={7,axisRange,1,-1,1,1,1,1}"
//	Execute "ModifyGizmo ModifyObject=axes0,property={8,axisRange,1,-1,-1,1,1,-1}"
//	Execute "ModifyGizmo ModifyObject=axes0,property={9,axisRange,-1,1,-1,1,1,-1}"
//	Execute "ModifyGizmo ModifyObject=axes0,property={10,axisRange,-1,1,1,1,1,1}"
//	Execute "ModifyGizmo ModifyObject=axes0,property={11,axisRange,-1,-1,1,1,-1,1}"
//	Execute "ModifyGizmo ModifyObject=axes0,property={-1,axisScalingMode,1}"
//	Execute "ModifyGizmo ModifyObject=axes0,property={-1,axisColor,0,0,0,1}"
//
//	Execute "ModifyGizmo ModifyObject=axes0,property={9,ticks,3}"
//	Execute "ModifyGizmo ModifyObject=axes0,property={8,ticks,3}"
////	Execute "ModifyGizmo ModifyObject=axes0,property={9,ticks,3}"
//	Execute "ModifyGizmo ModifyObject=axes0,property={5,ticks,3}"
//	Execute "ModifyGizmo ModifyObject=axes0,property={9,axisLabel,1}"
//	Execute "ModifyGizmo ModifyObject=axes0,property={8,axisLabel,1}"
////	Execute "ModifyGizmo ModifyObject=axes0,property={0,axisLabel,1}"
//	if (strlen(XaxisName))
//		Execute "ModifyGizmo ModifyObject=axes0,property={9,axisLabelText,\""+XaxisName+"\"}"
//		Execute "ModifyGizmo ModifyObject=axes0,property={9,axisLabelCenter,-0.3}"
//		Execute "ModifyGizmo ModifyObject=axes0,property={9,axisLabelDistance,0.1}"
//		Execute "ModifyGizmo ModifyObject=axes0,property={9,axisLabelScale,0.5}"
//		Execute "ModifyGizmo ModifyObject=axes0,property={9,axisLabelRGBA,0,0,0,1}"
//	endif
//	if (strlen(YaxisName))
//		Execute "ModifyGizmo ModifyObject=axes0,property={8,axisLabelText,\""+YaxisName+"\"}"
//		Execute "ModifyGizmo ModifyObject=axes0,property={8,axisLabelCenter,-0.5}"
//		Execute "ModifyGizmo ModifyObject=axes0,property={8,axisLabelDistance,0.1}"
//		Execute "ModifyGizmo ModifyObject=axes0,property={8,axisLabelScale,0.5}"
//		Execute "ModifyGizmo ModifyObject=axes0,property={8,axisLabelRGBA,0,0,0,1}"
//	endif
//	Execute "ModifyGizmo setDisplayList=0, object=surface0"
//	Execute "ModifyGizmo setDisplayList=1, object=axes0"
//	Execute "ModifyGizmo SETQUATERNION={-0.102663,-0.598222,-0.788212,-0.101455}"
//	Execute "ModifyGizmo autoscaling=1"
//	Execute "ModifyGizmo currentGroupObject=\"\""
//	Execute "ModifyGizmo compile"
//
////	Execute "ModifyGizmo showInfo"
////	Execute "ModifyGizmo infoWindow={2085,657,2558,920}"
//	Execute "ModifyGizmo bringToFront"
//	Execute "ModifyGizmo endRecMacro"
//End
//Function MakeGraph3dSurface()					// make a 3-d surface to go with cut
//	// Mac V. 3.0.0
//	if (exists("CreateSurfer") !=4)		// Do nothing if the Surface Plotter XOP is not available.
//		DoAlert 0, "Surface Plotter XOP must be installed"
//		return 1
//	endif
//	if (exists("sliceWaveAbs")!=1 || exists("sliceWaveLog")!=1)
//		DoAlert 0, "sliceWaveAbs or  sliceWaveLog does not exists, try MakeGraphOfSlice() first"
//		return 1
//	endif
//	Execute "CreateSurfer/K=1"
//	MoveWindow 7,53,675,571 		//	left, top, right, bottom
//	Execute "ModifySurfer  FactoryDefaults, Update=0"
//	Execute "ModifySurfer expandWindow=0"
//	Execute "ModifySurfer srcWave=sliceWaveAbs"
//	Execute "ModifySurfer secondWave=sliceWaveLog"
//	Execute "ModifySurfer surfaceColorWave=sliceWaveLog"
//	Execute "ModifySurfer  srcType=1,plotType=3,  setControlView=3"
//	Execute "ModifySurfer  theta=28.9,  phi=311.8,  zScale=1,  xStep=1,  yStep=1"
//	Execute "ModifySurfer  frame=895,  drawFrame=0, drawBox=1"
//	Execute "ModifySurfer  numXTicks=2,  numYTicks=2,  numZTicks=2, drawTicks=5"
//	Execute "ModifySurfer  noTicksX=0,  noTicksY=0,  noTicksZ=0"
//	Execute "ModifySurfer palette=RedWhiteBlue, reverseColorMap=1"
//	Execute "ModifySurfer  xmin=NaN,  xmax=NaN,  ymin=NaN,  ymax=NaN"
//	Execute "ModifySurfer  surfaceColorFill=1,  grids=1"
//	Execute "ModifySurfer  numContourLevels=12,  contour=1027"
//	Execute "ModifySurfer  marker=16,  markerSize=1,  markerColorType=4"
//	Execute "ModifySurfer  scatterDepthCue=1, imageType=1, BPP=16"
//	Execute "ModifySurfer  Update=1"
//	Execute "ModifySurfer fullUpdate"
//End


Function Set3dSurfaceZscale(zmax)
	Variable zmax											// max z value

	if (strlen(WinList("*", "","WIN:4096"))<1)		// no surface plot is up
		return NaN
	endif
	Wave SliceWaveAbs=SliceWaveAbs
	if (!WaveExists(SliceWaveAbs))
		DoAlert 0, "not in right Igor data folder, SliceWaveAbs does not exist here"
		return NaN
	endif
	zmax = (numtype(zmax)==1 && zmax<0) ? NaN : zmax
	if (numtype(zmax)==2)
		WaveStats/M=1/Q SliceWaveAbs
		zmax = roundSignificant(V_max,3)
		Prompt zmax,"max Z value"
		DoPrompt "max Z value",zmax
		if (V_flag)
			return NaN
		endif
	endif
	if (numtype(zmax)==2)
		return NaN
	endif

	String wName = SquaredVersionOfSurf(sliceWave)	// squared version of sliceWave
	Wave sliceSquare=$wName
	if (NumVarOrDefault(":reverse3dSurfacePlot",0))
		Variable xlo=DimOffset(sliceSquare,0), N=DimSize(sliceSquare,0), dx=DimDelta(sliceSquare,0)
		String xstr = GetDimLabel(sliceSquare,0,-1)
		ImageTransform flipRows sliceSquare
		SetScale/I x ((N-1)*dx+xlo),xlo,xstr, sliceSquare
	endif
	Duplicate/O sliceSquare, sliceWaveAbs				// these waves needed to make the 3d surface plot
	sliceWaveAbs = abs(sliceSquare)
	sliceSquare = sqrt(abs(sliceWaveAbs)) * sign(sliceSquare)
	WaveStats/M=1/Q sliceWaveAbs
	if (numtype(zmax)==1)								// gives auto-scaling
		zmax = V_max
		sliceSquare /= (4*zmax)
	elseif (V_max<zmax)									// need to shrink the surface
		Variable i=0, j=0
		FindUnusedCorner(SliceWaveAbs,i,j)
		SliceWaveAbs[i][j] = zmax
		sliceSquare /= (4*V_max)						// keep the colors bright
	elseif(V_max>zmax)									// need to expand the surface
		SliceWaveAbs = SliceWaveAbs[p][q] < zmax ? SliceWaveAbs[p][q] : NaN
		sliceSquare /= (4*zmax)
	endif
	Wave sliceWaveRGBA=sliceWaveRGBA
	ColorTab2Wave RedWhiteBlue
	Wave M_colors=M_colors
	Variable Ncolors = DimSize(M_colors,0)/2
	sliceWaveRGBA = M_colors[sliceSquare[p][q]*Ncolors + Ncolors][r]/65535.0
	sliceWaveRGBA[][][3] = 1
	sliceWaveRGBA = limit(sliceWaveRGBA,0,1)
	KillWaves/Z M_colors, sliceSquare,sliceSquare
	if (numtype(zmax)==0)								// gives auto-scaling
		Note/K SliceWaveAbs, ReplaceNumberByKey("zmaxGizmo",note(SliceWaveAbs),zmax,"=")
	endif
	return zmax
End
//
Static Function FindUnusedCorner(mat,i,j)		// finds indicies to a corner whose surrounding pixels are NaN (returns 1 is none found)
	Wave mat										// 2-d wave to be checked
	Variable &i, &j								// results passed back

	Variable ni=DimSize(mat,0)-1, nj=DimSize(mat,1)-1		// last index, not size in i,j directions
	if (numtype(mat[0][1]) && numtype(mat[1][0]))			// check the [0][0] corner
		i=0  ;  j=0
		return 0
	elseif (numtype(mat[ni-1][0]) && numtype(mat[ni][1]))		// check the [ni][0] corner
		i=ni  ;  j=0
		return 0
	elseif (numtype(mat[0][nj-1]) && numtype(mat[1][nj]))		// check the [0][nj] corner
		i=0  ;  j=nj
		return 0
	elseif (numtype(mat[ni][nj-1]) && numtype(mat[ni-1][nj]))	// check the [ni][nj] corner
		i=0  ;  j=nj
		return 0
	endif
	return 1																// could not find an un-occupied corner (i and j unchanged)
End


//Strconstant listOfAllRotations = "X;H;F;Y;Z;R_Cyl;Phi_Cyl;Axis_Cyl;Total;alpha(xx);alpha(xy);alpha(xz);alpha(yx);alpha(yy);alpha(yz);alpha(zx);alpha(zy);alpha(zz);GND;dRdX;dRdY;dRdZ;|gradR|;exx;eyy;ezz;exy;exz;eyz;a;b;c;alpha;beta;gamma"
//Strconstant listOfAllRotationsNoEps = "X;H;F;Y;Z;R_Cyl;Phi_Cyl;Axis_Cyl;Total;alpha(xx);alpha(xy);alpha(xz);alpha(yx);alpha(yy);alpha(yz);alpha(zx);alpha(zy);alpha(zz);GND;dRdX;dRdY;dRdZ;|gradR|"
Strconstant listOfAllRotations = "X;H;F;Y;Z;R_Cyl;Phi_Cyl;Axis_Cyl;Total;kappa(xx);kappa(xy);kappa(xz);kappa(yx);kappa(yy);kappa(yz);kappa(zx);kappa(zy);kappa(zz);alpha(xx);alpha(xy);alpha(xz);alpha(yx);alpha(yy);alpha(yz);alpha(zx);alpha(zy);alpha(zz);GND;dRdX;dRdY;dRdZ;|gradR|;exx;eyy;ezz;exy;exz;eyz;a;b;c;alpha;beta;gamma"
Strconstant listOfAllRotationsNoEps = "X;H;F;Y;Z;R_Cyl;Phi_Cyl;Axis_Cyl;Total;kappa(xx);kappa(xy);kappa(xz);kappa(yx);kappa(yy);kappa(yz);kappa(zx);kappa(zy);kappa(zz);GND;dRdX;dRdY;dRdZ;|gradR|"



Function MakePanelSlicer(sliceW) : Panel	// panel to go with GraphSlice
	Wave sliceW
	if (!WaveExists(sliceW))
		DoAlert 0, "MakePanelSlicer(), cannot put up panel by hand, it goes up automatically with MakeGraphOfSlice()"
		return 1
	endif
	if (strlen(WinList("GraphSlice",";","WIN:1"))<=0)
		DoAlert 0, "MakePanelSlicer(), needs to have GraphSlice already up, it goes up automatically with MakeGraphOfSlice()"
		return 1
	endif

	String fldr= GetWavesDataFolder(sliceW,1)		// this is the folder where everything is created and stored
	String noteStr = note(sliceWave)

	String str = fldr+"SliceValue"
	Variable/G $str							// ensure existance of SliceValue
	NVAR SliceValue=$str
	String/G $( fldr+"SliceNormal")		// ensure existance of SliceNormal
	SVAR SliceNormal=$( fldr+"SliceNormal")
	SliceNormal = SelectString(strlen(SliceNormal),"H",SliceNormal)
	String/G $(fldr+"SliceAngleAxis")		// ensure existance of SliceAngleAxis, the name of the rotation axis
	SVAR SliceAngleAxis=$(fldr+"SliceAngleAxis")
	SliceAngleAxis = SelectString(strlen(SliceAngleAxis),"X",SliceAngleAxis)
	Variable/G $(fldr+"phiSlice")			// ensure existance of phiSlice
	NVAR phi=$(fldr+"phiSlice")
	Variable/G $(fldr+"CutResolution")		// ensure existance of resolution
	NVAR CutResolution=$(fldr+"CutResolution")
	CutResolution = (CutResolution>0.0) ? CutResolution : 1.0
	str = fldr+"epsilonIsZero"
	Variable/G $str							// ensure existance of epsilonIsZero
	NVAR epsilonIsZero=$str

	if (strlen(WinList("PanelSlicer","","WIN:64")))
		DoWindow/K PanelSlicer							// kill panel if it exists
	endif
	NewPanel /W=(229,44,449,182)/K=1				// W=(left, top, right, bottom )
	DoWindow /C PanelSlicer
	AutoPositionWindow /M=0/R=GraphSlice PanelSlicer	// position panel by top most gizmo
	Variable mode, lo=ceil(NumberByKey("Hlo",noteStr,"=")), hi=floor(NumberByKey("Hhi",noteStr,"="))
	SetVariable setvarSlice,pos={1,2},size={75,18},proc=SetSliceProc,title=SliceNormal
	SetVariable setvarSlice,fSize=12
//	SetVariable setvarSlice,limits={lo,hi,1},value= SliceValue
	Variable stepSize=min(abs(hi-lo),1)
	stepSize = stepSize==0 ? 1 : stepSize
	if (stepSize<1)
		stepSize *= stepSize<1 ? (1-1e-14) : 1
	endif
	SetVariable setvarSlice,limits={lo,hi,stepSize},value= SliceValue
	SetVariable setvarPhi,pos={141,2},size={71,18},disable=1,proc=SetPhiCutProc,title="phi°"
	SetVariable setvarPhi,fSize=12,limits={-360,360,0},value=phiSlice
	SetVariable setvarResolution,pos={9,25},size={145,18},proc=SetResolutionProc,title="resolution (µm)"
	SetVariable setvarResolution,fSize=12,limits={0,inf,0},value=CutResolution
	mode = WhichListItem(SliceNormal,"X;H;F;Y;Z;phiZ")+1
	PopupMenu CutNormalPopUp,pos={81,1},size={35,20},proc=CutSlicePopUpProc
	PopupMenu CutNormalPopUp,mode=mode,popvalue=SliceNormal,value= #"\"X;H;F;Y;Z;phiZ\""
	SetVariable setvarPhi disable=!stringmatch(SliceNormal,"phiZ")
	PopupMenu CutAnglePopUp,pos={10,48},size={135,20},proc=CutSlicePopUpProc,title="rotation axis"
	if (exists("exx3d")==1)
		mode = WhichListItem(SliceAngleAxis,listOfAllRotations)+1
		PopupMenu CutAnglePopUp,mode=mode,popvalue=SliceAngleAxis,value= listOfAllRotations
	else
		mode = WhichListItem(SliceAngleAxis,listOfAllRotationsNoEps)+1
		PopupMenu CutAnglePopUp,mode=mode,popvalue=SliceAngleAxis,value= listOfAllRotationsNoEps
	endif

	PopupMenu popup3dSurface,pos={14,84},size={89,20},proc=PopMenu3dSurfaceProc,title="3d surface"
	PopupMenu popup3dSurface,fSize=12,mode=0,value= #"\"Show 3d surface;autoscale z;fix z scale\""
	PopupMenu popupBeamDirection,pos={14,110},size={89,20},proc=BeamDirectionProc,fSize=12,mode=0
	CheckBox epsilon0,pos={165,50},size={45,18},title="e=0",font="Symbol",fSize=18
	CheckBox epsilon0,variable= epsilonIsZero,proc=epsilonZeroProc
	CheckBox epsilon0,help={"set epsilon to 0 for computing alphas and GND"}
	NVAR flag = :reverse3dSurfacePlot
	BeamDirectionProc("",flag+1,"")
End
Function PopMenu3dSurfaceProc(ctrlName,popNum,popStr) : PopupMenuControl	//control for 3d-surface view
	String ctrlName
	Variable popNum
	String popStr
	strswitch(popStr)
		case "Show 3d surface":
			MakeGraph3dSurface()
			break
		case "autoscale z":
			Set3dSurfaceZscale(Inf)
			break
		case "fix z scale":
			Set3dSurfaceZscale(NaN)
			break
	endswitch
	return 0
End
Function BeamDirectionProc(ctrlName,popNum,popStr) : PopupMenuControl
	String ctrlName
	Variable popNum
	String popStr
	SetBeamDirection(strsearch(popStr,"left",0,2)>=0)
	return 0
End

Function epsilonZeroProc(cba) : CheckBoxControl
	STRUCT WMCheckboxAction &cba
	switch( cba.eventCode )
		case 2: // mouse up
			Variable checked = cba.checked
			break
		default:
			return 1
	endswitch
	String gName = "GraphSlice"
	Wave sliceW = ImageNameToWaveRef(gName, "sliceWave" )	// the 2-d wave displayed as an image and contours
	String fldr = GetWavesDataFolder(sliceW,1)
	NVAR SliceValue=$(fldr+"SliceValue")
	SetSliceProc("",SliceValue,"","SliceValue")					// force recalculation
	return 0
End


// this updates the slice surface.  It is called each time you want a new level.  To change the kind of slice (x, y, or z) call CutSlicePopUpProc
Function SetSliceProc(ctrlName,value,varStr,varName) : SetVariableControl
	String ctrlName
	Variable value
	String varStr
	String varName
	if (cmpstr("SliceValue",varName))
		Abort "SetSliceProc only callable from the button"
	endif

	String gName = "GraphSlice"
	Wave sliceW = ImageNameToWaveRef(gName, "sliceWave" )
	String fldr = GetWavesDataFolder(sliceW,1)
	NVAR CutResolution=$(fldr+"CutResolution")
	NVAR SliceValue=$(fldr+"SliceValue")
	SVAR SliceNormal=$(fldr+"SliceNormal")						// name of normal to the slice
	SVAR SliceAngleAxis=$(fldr+"SliceAngleAxis")					// name of desired rotation axis
	NVAR epsilonIsZero=$(fldr+"epsilonIsZero")					// flag to ignore epsilon when computing alphas or GND

	Wave surfNormal = $MakeDirectionVector(SliceNormal)

	if (strlen(SlicePlaneIn3d(CutResolution,surfNormal,SliceNormal,SliceValue,SliceAngleAxis,epsilonIsZero))==0)	// creates or recalcs sliceWave
		return 1
	endif
	String noteStr = note(sliceW)
	noteStr = ReplaceNumberByKey("SliceValue",noteStr,SliceValue,"=")
	Note/K sliceW
	Note sliceW, noteStr
	KillWaves/Z surfNormal

	String wName = SquaredVersionOfSurf(sliceW)					// squared version of sliceWave
	Wave sliceSquare=$wName
	if (NumVarOrDefault(":reverse3dSurfacePlot",0))
		Variable xlo=DimOffset(sliceSquare,0), N=DimSize(sliceSquare,0), dx=DimDelta(sliceSquare,0)
		String xstr = GetDimLabel(sliceSquare,0,-1)
		ImageTransform flipRows sliceSquare
		SetScale/I x ((N-1)*dx+xlo),xlo,xstr, sliceSquare
	endif
	Variable zmaxOld = NaN
	if (exists("sliceWaveAbs")==1)
		Wave sliceWaveAbs=sliceWaveAbs
		zmaxOld = NumberByKey("zmaxGizmo",note(sliceWaveAbs),"=")
	endif
	Duplicate/O sliceSquare, sliceWaveAbs, sliceWaveRGBA		// these waves needed to make the 3d surface plot
	sliceWaveAbs = abs(sliceSquare)
	WaveStats/Q/M=1 sliceWaveAbs
	Variable zmax = numtype(zmaxOld) ? V_max : zmaxOld
	SliceWaveAbs = SliceWaveAbs[p][q] < zmax ? SliceWaveAbs[p][q] : NaN

	Redimension/N=(-1,-1,4) sliceWaveRGBA						// make the RGBA wave for the 3d gizmo plot
	ColorTab2Wave RedWhiteBlue
	Wave M_colors=M_colors
	Variable Ncolors = DimSize(M_colors,0)/2
	sliceSquare = sqrt(abs(sliceWaveAbs)) * sign(sliceSquare)
	sliceSquare /= (4*zmax)
	sliceWaveRGBA = M_colors[sliceSquare[p][q]*Ncolors + Ncolors][r]/65535.0
	sliceWaveRGBA[][][3] = 1
	sliceWaveRGBA = limit(sliceWaveRGBA,0,1)
	KillWaves/Z M_colors, sliceSquare

	SurfacePlotStyle("gName")
	SurfacePlotStyle("")
//	KillWaves/Z surfNormal_
	return 0
End
//
Function SurfacePlotStyle(gName)
	String gName
	if (strlen(gName)>0 && strlen(WinList(gName,"","WIN:1"))<1)	// graph not on top, return
		return 1
	endif
	Wave image = ImageNameToWaveRef(gName,StringFromList(0,ImageNameList(gName,";")))
	if (!WaveExists(image))
		return 1
	endif
	if (!WaveInClass(image,"OrientationSliceWave*;RGBfullOrientationMap;epsilonSlice*"))
		return 1
	endif

	String noteStr = note(image)
	Variable SliceValue = NumberByKey("SliceValue",noteStr,"=")
	String SliceNormal = StringByKey("SliceNormal",noteStr,"=")
	String SliceAngleAxis = StringByKey("SliceAngleAxis",noteStr,"=")
	Variable epsilonIsZero = NumberByKey("epsilonIsZero",noteStr,"=")
	String fldr = GetWavesDataFolder(image,1)
	Variable rgbWave = WaveDims(image)==3 && DimSize(image,2)==3

	if (WhichListItem(SliceAngleAxis,"a;b;c;alpha;beta;gamma")>=0)	// display full auto-scale
		ModifyImage $NameOfWave(image) ctab= {*,*,Rainbow,0}
		ColorScale/C/N=colorScale1/A=RC/B=1/X=0/Y=0 image=$NameOfWave(image), heightPct=60, lowTrip=0.001,lblMargin=0,  tickUnit=1,ZisZ=1
		ColorScale/W=$gName/C/N=colorScale1 StringByKey("DUNITS",WaveInfo(image,0))
	elseif (!rgbWave)											// not an RGB wave
		// set range of color to be symmetric if it was not set by the user (this is a little tricky)
		String str = StringByKey("RECREATION",ImageInfo(gName,NameOfWave(image),0))	//ctab= {*,*,RedWhiteBlue,0}
		str = StringByKey("ctab",str,"=")
		str = str[2,strlen(str)-2]
		if (strlen(StringFromList(0,str,",")+StringFromList(1,str,","))>17 || strsearch(str,"*,",0)>=0)	// probably auto-scaled
			WaveStats/Q/M=1 image				// these 3 lines make red negative, blue positive
			Variable colRange = max(abs(V_min),abs(V_max))
			colRange = PickSymmetricRange(image,0.015)
			colRange = (colRange==0||numtype(colRange)) ? sqrt(2)*1e-4 : colRange
			ModifyImage $NameOfWave(image) ctab= {-colRange,colRange,RedWhiteBlue,0}	// red is negative, white=0, blue is postive
		endif
		ColorScale/C/N=colorScale1/A=RC/B=1/X=0/Y=0 image=$NameOfWave(image), heightPct=60, lowTrip=0.001,lblMargin=0,  tickUnit=1,ZisZ=1
		ColorScale/W=$gName/C/N=colorScale1 StringByKey("DUNITS",WaveInfo(image,0))
	elseif (WaveInClass(image,"OrientationSliceWave*"))
		Variable rmin=NumberByKey("rmin",noteStr,"="), rmax=NumberByKey("rmax",noteStr,"=")
		ColorScale/C/N=colorScale1/A=RC/B=1/X=0/Y=0 heightPct=60, lblMargin=0, lowTrip=0.001, tickUnit=1,ZisZ=1
		ColorScale/C/N=colorScale1 ctab={rmin,rmax,RedWhiteBlue,0}
		ColorScale/W=$gName/C/N=colorScale1 StringByKey("DUNITS",WaveInfo(image,0))
	endif

	// update the title for the plot
	NVAR phi=$(fldr+"phiSlice")
	String orientsText
	Variable fepsilon = WaveInClass(image,"epsilonSlice*") || (WhichListItem(SliceAngleAxis,"exx;eyy;ezz;exy;exz;eyz;a;b;c;alpha;beta;gamma")>=0)
	if (WhichListItem(SliceAngleAxis,"a;b;c")>=0)
		sprintf str "\\Zr150%s (nm)\\M",SliceAngleAxis
	elseif (WhichListItem(SliceAngleAxis,"alpha;beta;gamma")>=0)
		sprintf str "\\Zr150\\F'Symbol'%s∞\\F]0\\M",SliceAngleAxis[0,0]
	elseif (WhichListItem(SliceAngleAxis,"exx;eyy;ezz;exy;exz;eyz")>=0)
		sprintf str "\\Zr150\\F'Symbol'e\\F]0\\B%s\\M",SliceAngleAxis[1,2]
	elseif (strsearch(SliceAngleAxis, "kappa(",0)==0)	// using lattice curvature tensor component
		sprintf str "Lattice Curvature Tensor, \\Zr125\\F'Symbol'k\\F]0\\B%s\\M (rad/µm)",SliceAngleAxis[6,7]
	elseif (strsearch(SliceAngleAxis, "alpha(",0)==0)	// using dislocation tensor component
		sprintf str "Dislocation Tensor, \\Zr125\\F'Symbol'a\\F]0\\B%s\\M (rad/µm)",SliceAngleAxis[6,7]
		str += SelectString(epsilonIsZero,"",",  \\Zr150\\F'Symbol'e\\F]0\\M=0")
	elseif (stringmatch(SliceAngleAxis,"GND"))			// using dislocation tensor GND
		str = "Dislocation Tensor,   GND"
		str += SelectString(epsilonIsZero,"",",  \\Zr150\\F'Symbol'e\\F]0\\M=0")
	elseif (strsearch(SliceAngleAxis, "dRd",0)==0)		// using one component of gradient
		str = SliceAngleAxis[0,1]+"/"+SliceAngleAxis[2,3]+"  (rad/µm)"
	elseif (stringmatch(SliceAngleAxis,"|gradR|"))		// mag of gradient of total rotation angle
		str = "| grad(total Rotation) | (rad/µm)"
	else
		String listAxis = "X:X;H:H;F:F;Y:Y;Z:Z;R_Cyl:r\\X0\\S^\\M\\Bcylinder\\M;Phi_Cyl:\\F'Symbol'f\\F]0\\X0\\S^\\M\\B cylinder\\M;Axis_Cyl:cylindrical axis;Total:total angle"
		String axisText = StringByKey(SliceAngleAxis,listAxis)
		axisText += SelectString(strlen(axisText)==1,"","-axis")
		str = "lattice tilt about \\[0"+axisText+" (radian)"							// using a tilt
	endif
	TextBox/W=$gName/C/N=titleText/A=LT/B=1/F=0 str
	orientsText = "Orientation\r"+fldr
	if (stringmatch(SliceNormal,"phiZ"))
		sprintf str "Surface normal toward phi=%g",phi
	else
		sprintf str "Surface normal toward %s",SliceNormal
	endif
	AppendText/W=$gName/N=titleText str
	orientsText += "\r"+str
	sprintf str "cut perp to %s at %g µm from origin",SliceNormal,SliceValue
	AppendText/W=$gName/N=titleText str
	orientsText += "\r"+str
	if (stringmatch(SliceAngleAxis,"Total"))
		AppendText/W=$gName/N=titleText "total rotaion angle"
	endif
	if (strlen(WinList("GraphSliceOrients","","WIN:1" ))>0)
		TextBox/W=GraphSliceOrients/C/N=titleText/F=0/S=3/A=LT/X=2/Y=2 orientsText
	endif
	sprintf str, "%s [%s] _ [%s].SPE",StringByKey("filePrefix", noteStr,"="),StringByKey("positions", noteStr,"="),StringByKey("depths", noteStr,"=")
	str = SelectString(strlen(str)>12,"",str+"\r")								// skip line if no positions, depths, or filePrefix
	AppendText/W=$gName/N=titleText "\\Zr060"+str+IgorInfo(1)+"    "+StringByKey("fldrName", noteStr,"=")
	ModifyGraph gfMult=130,zero=2,minor=1,mirror=1
	Label left StringByKey("YaxisName",noteStr,"=")+" (\\U)"
	Label bottom StringByKey("XaxisName",noteStr,"=")+" (\\U)"
	DoUpdate
	SetAspectToSquarePixels(gName)
	return 0
End
//
//Function SetSliceProc(ctrlName,value,varStr,varName) : SetVariableControl
//	String ctrlName
//	Variable value
//	String varStr
//	String varName
//	if (cmpstr("SliceValue",varName))
//		Abort "SetSliceProc only callable from the button"
//	endif
//
//	String gName = "GraphSlice"
//	Wave sliceW = ImageNameToWaveRef(gName, "sliceWave" )
//	String fldr = GetWavesDataFolder(sliceW,1)
//	NVAR CutResolution=$(fldr+"CutResolution")
//	NVAR SliceValue=$(fldr+"SliceValue")
//	SVAR SliceNormal=$(fldr+"SliceNormal")						// name of normal to the slice
//	SVAR SliceAngleAxis=$(fldr+"SliceAngleAxis")					// name of desired rotation axis
//
//	Wave surfNormal = $MakeDirectionVector(SliceNormal)
//	SlicePlaneIn3d(CutResolution,surfNormal,SliceNormal,SliceValue,SliceAngleAxis)	// creates or recalcs sliceWave
//
//	String noteStr = note(sliceW)
//	noteStr = ReplaceNumberByKey("SliceValue",noteStr,SliceValue,"=")
//	Note/K sliceW
//	Note sliceW, noteStr
//	KillWaves/Z surfNormal
//
//	String wName = SquaredVersionOfSurf(sliceW)					// squared version of sliceWave
//	Wave sliceSquare=$wName
//	if (NumVarOrDefault(":reverse3dSurfacePlot",0))
//		Variable xlo=DimOffset(sliceSquare,0), N=DimSize(sliceSquare,0), dx=DimDelta(sliceSquare,0)
//		String xstr = GetDimLabel(sliceSquare,0,-1)
//		ImageTransform flipRows sliceSquare
//		SetScale/I x ((N-1)*dx+xlo),xlo,xstr, sliceSquare
//	endif
//	Duplicate/O sliceSquare, sliceWaveAbs, sliceWaveRGBA		// these waves needed to make the 3d surface plot
//	sliceWaveAbs = abs(sliceSquare)
//
//	Redimension/N=(-1,-1,4) sliceWaveRGBA						// make the RGBA wave for the 3d gizmo plot
//	ColorTab2Wave RedWhiteBlue
//	Wave M_colors=M_colors
//	Variable Ncolors = DimSize(M_colors,0)/2
//	sliceSquare = sqrt(abs(sliceWaveAbs)) * sign(sliceSquare)
//	WaveStats/Q/M=1 sliceWaveAbs
//	sliceSquare /= (4*V_max)
//	sliceWaveRGBA = M_colors[sliceSquare[p][q]*Ncolors + Ncolors][r]/65535.0
//	sliceWaveRGBA[][][3] = 1
//	sliceWaveRGBA = limit(sliceWaveRGBA,0,1)
//	KillWaves/Z M_colors, sliceSquare
//
//	// set range of color to be symmetric if it was not set by the user (this is a little tricky)
//	String str = StringByKey("RECREATION",ImageInfo("","sliceWave",0))	//ctab= {*,*,RedWhiteBlue,0}
//	str = StringByKey("ctab",str,"=")
//	str = str[2,strlen(str)-2]
//	if (strlen(StringFromList(0,str,",")+StringFromList(1,str,","))>17 || strsearch(str,"*,",0)>=0)	// probably auto-scaled
//		WaveStats/Q/M=1 sliceW				// these 3 lines make red negative, blue positive
//		Variable colRange = max(abs(V_min),abs(V_max))
//		colRange = PickSymmetricRange(sliceW,0.015)
//		colRange = (colRange==0||numtype(colRange)) ? sqrt(2)*1e-4 : colRange
//		ModifyImage sliceWave ctab= {-colRange,colRange,RedWhiteBlue,0}	// red is negative, white=0, blue is postive
//	endif
//
//	// update the title for the plot
//	NVAR phi=$(fldr+"phiSlice")
//	String listAxis = "X:X;H:H;F:F;Y:Y;Z:Z;R_Cyl:r\\X0\\S^\\M;Phi_Cyl:\\F'Symbol'f\\F]0\\X0\\S^\\M;Axis_Cyl:axis;Total:total angle"
//
//	String orientsText
//	Variable tensor = strsearch(SliceAngleAxis, "alpha(",0)==0 || stringmatch(SliceAngleAxis,"GND")
//	String scale = StringByKey("DUNITS",WaveInfo(sliceWave,0))
//	if (strsearch(SliceAngleAxis, "alpha(",0)==0)				// using dislocation tensor component
//		sprintf str "Dislocation Tensor, \\F'Symbol'a\\F]0\\B%s\\M (rad/µm)",SliceAngleAxis[6,7]
//		TextBox/W=$gName/C/N=titleText str
//		ColorScale/W=$gName/C/N=colorScale1 scale
//	elseif (stringmatch(SliceAngleAxis,"GND"))				// using dislocation tensor GND
//		TextBox/W=$gName/C/N=titleText "Dislocation Tensor,   GND"
//		ColorScale/W=$gName/C/N=colorScale1 scale
//	elseif (strsearch(SliceAngleAxis, "dRd",0)==0)				// using one component of gradient
//		TextBox/W=$gName/C/N=titleText SliceAngleAxis+"  (rad/µm)"
//		ColorScale/W=$gName/C/N=colorScale1 scale
//	elseif (stringmatch(SliceAngleAxis,"|gradR|"))				// mag of gradient of total rotation angle
//		TextBox/W=$gName/C/N=titleText "| grad(total Rotation) | (rad/µm)"
//		ColorScale/W=$gName/C/N=colorScale1 scale
//	else
//		TextBox/W=$gName/C/N=titleText "lattice tilt (radian)"		// using a tilt
//		ColorScale/W=$gName/C/N=colorScale1 scale
//	endif
//	AppendText/W=$gName/N=titleText fldr
//	orientsText = "Orientation\r"+fldr
//	if (stringmatch(SliceNormal,"phiZ"))
//		sprintf str "Surface normal toward phi=%g",phi
//	else
//		sprintf str "Surface normal toward %s",SliceNormal
//	endif
//	AppendText/W=$gName/N=titleText str
//	orientsText += "\r"+str
//	sprintf str "cut perp to %s at %g µm from origin",SliceNormal,SliceValue
//	AppendText/W=$gName/N=titleText str
//	orientsText += "\r"+str
//	if (stringmatch(SliceAngleAxis,"Total"))
//		AppendText/W=$gName/N=titleText "total rotaion angle"
//	elseif (!tensor)
//		sprintf str "rotation axis is along \\[0%s direction",StringByKey(SliceAngleAxis,listAxis)
//		AppendText/W=$gName/N=titleText str
//	endif
//	if (strlen(WinList("GraphSliceOrients","","WIN:1" ))>0)
//		TextBox/W=GraphSliceOrients/C/N=titleText/F=0/S=3/A=LT/X=2/Y=2 orientsText
//	endif
//	KillWaves/Z surfNormal_
//	return 0
//End
Static Function PickSymmetricRange(w,frac)
	Wave w
	Variable frac
	Duplicate/O w, PickSymmetricRangeW
	Redimension/N=(numpnts(PickSymmetricRangeW)) PickSymmetricRangeW
	SetScale/P x 0,1,"", PickSymmetricRangeW
	Sort PickSymmetricRangeW PickSymmetricRangeW
	WaveStats/M=1/Q PickSymmetricRangeW
	Variable n = V_npnts-1
	Variable xlo=PickSymmetricRangeW[n*frac], xhi=PickSymmetricRangeW[n*(1-frac)]
	KillWaves/Z PickSymmetricRangeW
	return max(abs(xlo),abs(xhi))
End
Function/S SquaredVersionOfSurf(wav)		// return the name of a new wave that is the squared version of the input
	Wave wav								// a 2d wave that needs to be squared
	Variable iN=DimSize(wav,0), jN=DimSize(wav,1)	// dimenstions of the input wave
	Variable i,j
	String wName, noteStr

	WaveStats/Q/M=1 wav
	if ((V_npnts+V_numINFs)<=0)						// no valid points, just return one the same size
		wName = UniqueName("surface_",1,0)
		i = max(iN,jN)
		Make/N=(i,i) $wName
		Wave surf = $wName
		noteStr = note(wav)
		Note surf, noteStr
		CopyScales/P wav, surf
		return GetWavesDataFolder(surf,2)
	endif

	wName = UniqueName("vec_",1,0)
	Make/N=(iN)/D $wName
	Wave veci = $wName
	wName = UniqueName("vec_",1,0)
	Make/N=(jN)/D $wName
	Wave vecj = $wName

	// find the first i, that has a valid number (is not NaN)
	Variable dx=DimDelta(wav,0), dy=DimDelta(wav,1)
	Variable ilo, jlo, ihi, jhi							// new size of surace
	for (jlo=0;jlo<jN;jlo+=1)
		veci = wav[p][jlo]
		WaveStats/Q/M=1 veci
		if (V_numNans<iN)
			break
		endif
	endfor
	for (ilo=0;ilo<iN;ilo+=1)
		vecj = wav[ilo][p]
		WaveStats/Q/M=1 vecj
		if (V_numNans<jN)
			break
		endif
	endfor
	for (jhi=jN-1;jhi>=0;jhi-=1)
		veci = wav[p][jhi]
		WaveStats/Q/M=1 veci
		if (V_numNans<iN)
			break
		endif
	endfor
	for (ihi=iN-1;ihi>=0;ihi-=1)
		vecj = wav[ihi][p]
		WaveStats/Q/M=1 vecj
		if (V_numNans<jN)
			break
		endif
	endfor

	Variable iNew, jNew, x0,y0
	Variable size = max(abs((jhi-jlo)*dy), abs((ihi-ilo)*dx))
	iNew = ceil(size/dx)+1				// new dimension in i
	jNew = ceil(size/dy)+1				// new dimension in j
	x0 = ilo*dx+DimOffset(wav,0)
	y0 = jlo*dy+DimOffset(wav,1)
	//	printf "i=[%d,%d],   j=[%d,%d]\r",ilo,ihi,jlo,jhi
	//	printf "size=%g,   x0=%g,  y0=%g,  iNew=%d,  jNew=%d\r",size,x0,y0,iNew,jNew
	wName = UniqueName("surface_",1,0)
	Make/N=(iNew,jNew) $wName
	Wave surf = $wName
	surf = NaN
	CopyScales/P wav, surf
	SetScale/P x, x0, dx, surf
	SetScale/P y, y0, dy, surf
	surf[0,ihi-ilo][0,jhi-jlo] = wav[ilo+p][jlo+q]	// set the values
	noteStr = note(wav)
	Note surf, noteStr
	KillWaves/Z veci, vecj
	return GetWavesDataFolder(surf,2)		// return full path for assingment to a wave
End



// this is called by the popup that selects the direction of the cut and the direction of the angle for sliceWave
Function CutSlicePopUpProc(ctrlName,popNum,popStr) : PopupMenuControl
	String ctrlName
	Variable popNum
	String popStr

	String gName = "GraphSlice"
	Wave sliceW = ImageNameToWaveRef(gName, "sliceWave" )	// the 2-d wave displayed as an image and contours
	String fldr = GetWavesDataFolder(sliceW,1)
	NVAR SliceValue=$(fldr+"SliceValue")
	SVAR SliceNormal=$(fldr+"SliceNormal")					// name of normal to sliceWave
	SVAR SliceAngleAxis=$(fldr+"SliceAngleAxis")				// name of desired rotation axis

	if (stringmatch(ctrlName,"CutNormalPopUp"))				// called by a change to the normal to sliceWave
		SliceNormal = popStr									// set the name of the normal to sliceWave
		SetVariable setvarPhi disable=!stringmatch(SliceNormal,"phiZ")
	elseif (stringmatch(ctrlName,"CutAnglePopUp"))			// called by a change to the requested rotaion axis {X,H,F,Y,Z,phi,axis,total}
		SliceAngleAxis = popStr									// name of desired rotation axis
	endif

	if(SetSliceProc("",SliceValue,"","SliceValue"))				// force recalculation
		return 1
	endif
	String noteStr = note(sliceW)
	Label/W=$gName left StringByKey("YaxisName",noteStr,"=")+" (\\U)"
	Label/W=$gName bottom StringByKey("XaxisName",noteStr,"=")+" (\\U)"
	DoUpdate					// this is needed to refresh values for call to SetAspectToSquarePixels()
	SetAspectToSquarePixels(gName)
	if (strlen(WinList("GraphSliceOrients","","WIN:1" ))>0)
		gName = "GraphSliceOrients"
		Label/W=$gName left StringByKey("YaxisName",noteStr,"=")+" (\\U)"
		Label/W=$gName bottom StringByKey("XaxisName",noteStr,"=")+" (\\U)"
		DoUpdate					// this is needed to refresh values for call to SetAspectToSquarePixels()
		SetAspectToSquarePixels(gName)
	endif
	Variable lo=NumberByKey("nlo",noteStr,"="), hi=NumberByKey("nhi",noteStr,"=")
	if (strlen(WinList("PanelSlicer","","WIN:64"))>1)
		Variable stepSize=min(abs(hi-lo),1)
		stepSize = stepSize==0 ? 1 : stepSize
		if (stepSize<1)
			stepSize *= stepSize<1 ? (1-1e-14) : 1
		endif
		SetVariable setvarSlice,win=PanelSlicer,limits={lo,hi,stepSize}
//		SetVariable setvarSlice,win=PanelSlicer,limits={lo,hi,1}
		SetVariable setvarSlice,win=PanelSlicer,title=SelectString(stringmatch(SliceNormal,"phiZ"),SliceNormal,"r")
	endif
	return 0
End

// This function replaced by SetAspectToSquarePixels() in ImageDisplayScaling.ipf
//
//Function FixAspectRatio(graphName)		// set aspect ratio of graph to match the axes
//	String graphName
//	if (strlen(graphName)<1)
//		graphName = StringFromList(0,WinList("*",";","WIN:1"))
//	endif
//	if (WinType(graphName)!=1)
//		DoAlert 0, "FixAspectRatio(), unable to get size of bottom axis"
//		return 1
//	endif
//	GetAxis /W=$graphName/Q bottom
//	if (V_flag)
//		DoAlert 0, "FixAspectRatio(), unable to get size of bottom axis"
//		return 1
//	endif
//	Variable xRange=abs(V_max-V_min)
//	GetAxis /W=$graphName/Q left
//	if (V_flag)
//		DoAlert 0, "FixAspectRatio(), unable to get size of left axis"
//		return 1
//	endif
//	Variable yRange=abs(V_max-V_min)
//	Variable aspect = yRange/xRange
//	//	printf "xRange=%g, yRange=%g,  aspect=%g\r",xRange,yRange,aspect
//	if (aspect>=1)
//		ModifyGraph/W=$graphName width={Aspect,1/aspect},height=0		// tall & thin
//	else
//		ModifyGraph/W=$graphName height={Aspect,aspect},width=0		// short & fat 		
//	endif
//End
Function SetPhiCutProc(ctrlName,varNum,varStr,varName) : SetVariableControl
	String ctrlName
	Variable varNum
	String varStr
	String varName
	CutSlicePopUpProc("",0,"")
//	XHFPopMenuProc("",-1,"phi")
End
Function SetResolutionProc(ctrlName,varNum,varStr,varName) : SetVariableControl // proc for CutResolution
	String ctrlName
	Variable varNum
	String varStr
	String varName
	if (varNum<=0)
		DoAlert 0, "you set the resolution to ≤ 0, what were you thinking?"
	elseif (abs(varNum)<0.1)
		DoAlert 0, "Caution, this fine a resolution will take a long (really long) time to interpolate!"
	endif
End



// Jan 24, 2007
// Changed SlicePlaneIn3d(resolution,normalW,SliceValue,rotation)
// also makes sliceOrients, a 3-color image, red is component || to 001, green is 011 component, blue is 111 component
//
// Jan 26, 2006
// Changed SlicePlaneIn3d(resolution,normalW,SliceValue,rotation)
// The components of SliceAngleAxis were incorrect.  Since the matricies are in "sample" coordinate 
// system (X,-H,-F), the correct unit vectors for the axes are:
// Hhat=(0,-1,0);  Fhat=(0,0,-1);  Xhat=(1,0,0);  Yhat=(0,-1,1)/sqrt(2);  Zhat=(0,-1,-1)/sqrt(2)
//
// Jun 7, 2007
// Changed SlicePlaneIn3d(resolution,normalW,SliceValue,rotation)
// made it work with R2dX, R2dH, R2dF for a single X
//
// Jul 3, 2007
// Changed SlicePlaneIn3d(resolution,normalW,SliceValue,rotation)
// for calculating alpha, use full formula with epsilon and kappa
//
// Jul 15, 2007
// added epsilonIsZero, to control use of epsilon when computing alphas and GND
//
Function/S SlicePlaneIn3d(resolution,normalW,SliceNormal,SliceValue,rotation,epsilonIsZero)	// calculate a cut, result goes into sliceWave
	Variable resolution			// resolution, grid spacing of output wave (µm)
	Wave normalW				// surface normal of cut plane, (need not be of length 1)
	String SliceNormal			// name of normal to the slice
	Variable SliceValue			// intercept displaced along normal
	String rotation				// specifys the desired rotation axis, can be one of "X;H;F;Y;Z;R_Cyl;Phi_Cyl;Axis_Cyl;Total"
	Variable epsilonIsZero		// flag, true means ignore epsilon when computing alpha and GND (i.e. assume that epsilon is zero)

	Wave R3dX=R3dX, R3dH=R3dH, R3dF=R3dF,  R2dX=R2dX, R2dH=R2dH, R2dF=R2dF
	Wave Rvol=R3dX
	String errstr=""
	if (!WaveExists(R3dX) || !WaveExists(R3dH) || !WaveExists(R3dF))
		Wave Rvol=R2dX
		if (!WaveExists(R2dX) || !WaveExists(R2dH) || !WaveExists(R2dF))
			errstr = "Error in SlicePlaneIn3d(), R3dX or R3dH or R3dF absent"
		endif
	elseif (!WaveExists(normalW))
		errstr = "SlicePlaneIn3d(), normal does not exists"
	elseif (numtype(norm(normalW)))
		errstr = "SlicePlaneIn3d(), normal is invalid"
	elseif (norm(normalW)<1e-10)
		errstr = "SlicePlaneIn3d(), normal is too short"
	elseif (!(resolution>0))
		errstr = "SlicePlaneIn3d(), resolution = "+num2str(resolution)
	elseif (numtype(SliceValue))
		errstr = "SlicePlaneIn3d(), SliceValue = "+num2str(SliceValue)
	elseif (WhichListItem(rotation,listOfAllRotations)<0)
		errstr = "SlicePlaneIn3d(), rotation = '"+rotation+"'"
	elseif (numtype(epsilonIsZero))
		errstr = "SlicePlaneIn3d(), epsilonIsZero = "+num2str(epsilonIsZero)
	endif
	if (strlen(errstr))
		DoAlert 0, errstr
		return ""
	endif

	Wave totalPeakIntensity3d=totalPeakIntensity3d, totalIntensity3d=totalIntensity3d
	if (WaveExists(totalIntensity3d))
		Wave intensWave = totalIntensity3d
	elseif (WaveExists(totalPeakIntensity3d))
		Wave intensWave = totalPeakIntensity3d
	else
		Wave intensWave = $""
	endif

	String SliceWaveUnit
	Variable ctensor=0, strain=0						// flags, need lattice curvature tensor, or strain tensor (epsilon)
	Variable grad=0, calcAxis=0							// flags, neeed gradient of rotation angle, or need to calculate the rotation axis normal (cylindrical coordinates)
	Variable ddtensor=0									// flag, need the dislocation density tensor
	Variable aRotation=0									// flags, just a plain rotation
	Variable i0,j0											// specifies a particular component
	if (stringmatch(rotation,"GND"))
		ctensor = 1										// need the curvature tensor
		strain = 1											// needs alpha, so needs both kappa and epsilon
		ddtensor = 1										// and needs the dislocation density tensor
		i0 = -1
		SliceWaveUnit = "dislocations / cm\\S2"
	elseif (strsearch(rotation,"alpha(",0)==0)
		ctensor = 1										// need the dislocation density tensor, so need curvature tensor
		strain = 1											// and needs epsilon
		ddtensor = 1										// and needs the dislocation density tensor
		if (xyComponent2ij(rotation,i0,j0))			// get i0,j0 index into alpha[i0][j0]
			return ""
		endif
		SliceWaveUnit = "radian / µm"
	elseif (strsearch(rotation,"kappa(",0)==0)
		ctensor = 1										// this cut will neede the lattice curvature tensor
		if(xyComponent2ij(rotation,i0,j0))			// get i0,j0 index into kappa[i0][j0]
			return ""
		endif
		SliceWaveUnit = "radian / µm"
	elseif (strsearch(rotation,"dRd",0)==0 || stringmatch(rotation,"|gradR|"))
		grad = 1											// this cut will need the gradient of the total rotation
		i0 = WhichListItem(rotation,"dRdX;dRdY;dRdZ;|gradR|")
		if (i0<0)
			return ""
		endif
		i0 = i0==3 ? -1 : i0								// use -1 to flag magnitude
		SliceWaveUnit = "radian / µm"
	elseif(WhichListItem(rotation,"exx;eyy;ezz;exy;exz;eyz")>=0)
		strain = 1											// need the strain tensor
		SliceWaveUnit = ""								// no unit for epsilon
		Wave e3d=$(rotation+"3d")
		if (!WaveExists(e3d))
			return ""
		endif
	elseif(WhichListItem(rotation,"a;b;c")>=0)
		strain = 1
		SliceWaveUnit = "nm"							// lattice parameter (units=nm)
		Wave e3d=$(rotation+"3d")
		if (!WaveExists(e3d))
			return ""
		endif
	elseif(WhichListItem(rotation,"alpha;beta;gamma")>=0)
		strain = 1
		SliceWaveUnit = "°"								// lattice parameter (units=°)
		Wave e3d=$(rotation+"3d")
		if (!WaveExists(e3d))
			return ""
		endif
	elseif (!stringmatch(rotation,"Total"))				// need to compute angle component
		aRotation = 1
		Wave rotationAxis = $MakeDirectionVector(rotation)// want rotation component about rotationAxis[]
		calcAxis = numtype(rotationAxis[0])!=0		// axis of rotation is not constant, needs to be re-calculated for each point
		SliceWaveUnit = "radian"
	else
		aRotation = 1
		SliceWaveUnit = "radian"
	endif
	Wave normal=$MakeUnique3Vector(normalW)		// make a new 3-vector, and optionally set it to preset[]
	normalize(normal)

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

	Variable xlo,xhi, hlo,hhi, flo, fhi					// extent of box in x,h,f
	xlo = NumberByKey("Xlo",noteStr,"=")
	xhi = NumberByKey("Xhi",noteStr,"=")
	hlo = NumberByKey("Hlo",noteStr,"=")
	hhi = NumberByKey("Hhi",noteStr,"=")
	flo = NumberByKey("Flo",noteStr,"=")
	fhi = NumberByKey("Fhi",noteStr,"=")

	//	Need point representing center of cut surface.  It is:     xhf0[] + normal[]*SliceValue
	//	We are assuming that xhf0[] in the 3d volume is the origin that we want.
	// 	Then need x and y directions in the slice, this depends upon the surfNormal[]
	//	call these directions sx and sy, along the x and y in the slice
	//
	// get the sx and sy axes, the 3vectors giving x, y directions for the slice, and the center of the slice
	Wave xhf0=$MakeUnique3Vector($"")				// origin
	Wave hat=$MakeUnique3Vector($"")				// approximate direction of x-axis
	Wave sx=$MakeUnique3Vector($"")					// direction of x-axis in surface
	Wave sy=$MakeUnique3Vector($"")					// direction of y-axis in surface
	Wave center=$MakeUnique3Vector($"")			// center of the surface
	Wave xhf=$MakeUnique3Vector($"")				// a point on surface in XHF space
	Wave Rod=$MakeUnique3Vector($"")				// direction of Rodriques vector
	Wave rvec=$MakeUnique3Vector($"")				// vector from xhf0[] to point on slice
	Wave rhat=$MakeUnique3Vector($"")				// like xhf, but axial part removed

	xhf0 =  {x0,h0,f0}
	center =xhf0[p] + SliceValue*normal[p]// center[] = {X0,H0,F0} + SliceValue*normal[], center of slice
	if (abs(normal[0])<=abs(normal[1]))
		hat = {1,0,0}										// use xhat to define sx
	else	
		hat = {0,1,0}										// use yhat to define sx
	endif
	Variable dot
	dot = MatrixDot(normal,hat)
	sx = hat[p] - dot*normal[p]							// sx is unit vector perpendicular to normal, toward hat direction
	normalize(sx)
	Cross normal, hat
	Wave W_Cross=W_Cross
	sy = W_Cross
	normalize(sy)

	// next what is the range of the plane, what are its starting and stopping values along sx and sy
	Variable sxlo, sxhi, sylo, syhi,snlo,snhi				// range of plane, used to set the scaling, use cube corners as limits (snlo,snhi refers to distance along normal)
	Variable nsx,nsy										// sliceWave[nsx][nsy]
	FindRangeForBox(xlo,xhi,hlo,hhi,flo,fhi,sx,xhf0,sxlo,sxhi)		// determine max range along sx[] (relative to xhf0[])
	FindRangeForBox(xlo,xhi,hlo,hhi,flo,fhi,sy,xhf0,sylo,syhi)		// determine max range along sx[] (relative to xhf0[])
	FindRangeForBox(xlo,xhi,hlo,hhi,flo,fhi,normal,xhf0,snlo,snhi)// determine max range along normal[] (relative to xhf0[])
	nsx = ceil((sxhi-sxlo)/resolution)+1				// number of points in sliceWave
	nsy = ceil((syhi-sylo)/resolution)+1
	if (ItemsInList(GetRTStackInfo(0))<2)
		printf "•SlicePlaneIn3d(resolution=%g, surfNorm={%.3g, %.3g, %.3g}, SliceValue=%g, rotationAxis='%s',epsilonIsZero=%d)\r",resolution,normalW[0],normalW[1],normalW[2],SliceValue,rotation,epsilonIsZero
		printf " the rotation axis = (%g, %g, %g)\r",rotationAxis[0],rotationAxis[1],rotationAxis[2]
		printf "the range of slice is:     sx = [%g, %g] (%d pts),  sy=[%g, %g] (%d pts),  along normal = [%.4g, %.4g]\r",sxlo,sxhi,nsx,sylo,syhi,nsy,snlo,snhi
		printf "normal={%.3g, %.3g, %.3g},     sx={%.3g, %.3g, %.3g},     sy={%.3g, %.3g, %.3g}\r",normal[0],normal[1],normal[2],sx[0],sx[1],sx[2],sy[0],sy[1],sy[2]
		printf "xhf0={%.3g, %.3g, %.3g},     axis={%.3g, %.3g, %.3g}\r",xhf0[0],xhf0[1],xhf0[2],axis[0],axis[1],axis[2]
	endif
	noteStr = ReplaceNumberByKey("nlo",noteStr,snlo,"=")
	noteStr = ReplaceNumberByKey("nhi",noteStr,snhi,"=")
	Make/N=(nsx,nsy)/O sliceWave,sliceWaveIntens
	Make/N=(nsx,nsy,3)/O/D sliceOrientsD
	SetScale/I x sxlo,sxhi,"µm", sliceWave,sliceWaveIntens,sliceOrientsD
	SetScale/I y sylo,syhi,"µm", sliceWave,sliceWaveIntens,sliceOrientsD
	noteStr = ReplaceNumberByKey("SliceValue",noteStr,SliceValue,"=")
	noteStr = ReplaceStringByKey("SliceNormal",noteStr,SliceNormal,"=")
	noteStr = ReplaceStringByKey("SliceAngleAxis",noteStr,rotation,"=")
	noteStr = ReplaceStringByKey("XaxisName",noteStr,axisName(sx),"=")
	noteStr = ReplaceStringByKey("YaxisName",noteStr,axisName(sy),"=")
	noteStr = ReplaceNumberByKey("epsilonIsZero",noteStr,epsilonIsZero,"=")
	Note/K sliceWave, ReplaceStringByKey("waveClass",noteStr,"OrientationSliceWave","=")
	Note/K sliceWaveIntens, ReplaceStringByKey("waveClass",noteStr,"InensitySliceWave","=")
	Note/K sliceOrientsD, ReplaceStringByKey("waveClass",noteStr,"RGBfullOrientationMap","=")
	SetScale d 0,0,"radian", sliceOrientsD
	SetScale d 0,0,"", sliceWaveIntens
	sliceWaveIntens = NaN
	sliceOrientsD = NaN
	SetScale d 0,0,SliceWaveUnit, sliceWave

	// next fill the sliceWave
	Variable i,j, xx,yy, angle, GND, magGrad, peakIntensity
	Make/N=(3,3)/D/O SlicePlaneIn3d_alpha, SlicePlaneIn3d_kappa
	Make/N=(3)/O/D SlicePlaneIn3d_grad
	Wave alpha=SlicePlaneIn3d_alpha, kappa=SlicePlaneIn3d_kappa, gradR=SlicePlaneIn3d_grad
	for (j=0;j<nsy;j+=1)
		yy = DimOffset(sliceWave,1) + j*DimDelta(sliceWave,1)
		for (i=0;i<nsx;i+=1)
			xx = DimOffset(sliceWave,0) + i*DimDelta(sliceWave,0)
			xhf = center[p] + xx*sx[p] + yy*sy[p]		// point in XHF space corresponding to point (i,j) in slice
			if (xhf[0] == DimOffset(R3dX,0) + (DimSize(R3dX,0)-1)*DimDelta(R3dX,0))	// if xhf[0] at high edge, reduce tiny bit
				xhf[0] -= DimDelta(R3dX,0)*1e-9													// This section needed so Interp3d does not clip off the last plane
			endif																						// and ditto for dimensions 1 and 2
			if (xhf[1] == DimOffset(R3dX,1) + (DimSize(R3dX,1)-1)*DimDelta(R3dX,1))
				xhf[1] -= DimDelta(R3dX,1)*1e-9
			endif
			if (xhf[2] == DimOffset(R3dX,2) + (DimSize(R3dX,2)-1)*DimDelta(R3dX,2))
				xhf[2] -= DimDelta(R3dX,2)*1e-9
			endif

			if (WaveDims(Rvol)==2)
				Rod[0] = Interp2D(R2dX, xhf[1],xhf[2])			// Rodriques vector at xhf
				Rod[1] = Interp2D(R2dH, xhf[1],xhf[2])
				Rod[2] = Interp2D(R2dF, xhf[1],xhf[2])
				peakIntensity = NaN
//				peakIntensity = Interp2D(intensWave, xhf[1],xhf[2])	// need to fix this
			else
				Rod[0] = Interp3d(R3dX, xhf[0],xhf[1],xhf[2])	// Rodriques vector at xhf
				Rod[1] = Interp3d(R3dH, xhf[0],xhf[1],xhf[2])
				Rod[2] = Interp3d(R3dF, xhf[0],xhf[1],xhf[2])
				peakIntensity = Interp3d(intensWave, xhf[0],xhf[1],xhf[2])
			endif

			sliceWaveIntens[i][j] = peakIntensity
			sliceOrientsD[i][j][] = Rod[r]
			if (aRotation)												// looking for an angle, not a dislocation denstity tensor component
				// compute the angle for this point and save it in sliceWave[i][j]
				angle = 2*atan(normalize(Rod))					// Rodriques angle (radian), Rod is now normalized
				if (calcAxis)											// need to recalculate rotationAxis[] each time
					rvec = xhf - xhf0									// rvec is now direction of point relative to the origin of cylinder
					dot = MatrixDot(axis,rvec)
					rhat = rvec - dot*axis							// put rhat in the plane where z along axis is 0, only r & phi
					normalize(rhat)
					if (stringmatch(rotation,"R_Cyl"))				// rotation around axis sticking out radially from axis, r^ direction
						rotationAxis = rhat
					elseif (stringmatch(rotation,"Phi_Cyl"))		// rotation around axis in phi^ direction, probably the biggest for an indent
						Cross axis, rhat								// rotationAxis = axis x rHat - phi^
						rotationAxis = W_Cross
					elseif (stringmatch(rotation,"Axis_Cyl"))		// rotation around axis of cylinder, circulating direction
						rotationAxis = axis
					else
						Abort "SlicePlaneIn3d(), invalid rotation"	// should not be possible to reach this line
					endif
				endif

				if (i==floor(nsx/2) && j==floor(nsy/2) && ItemsInList(GetRTStackInfo(0))<2)
					printf "at center of sliceWave[%d][%d],      angle=%.3g°   Rod={%.3g, %.3g, %.3g},      XHF={%.3g, %.3g, %.3g}",i,j,angle*180/PI,Rod[0],Rod[1],Rod[2],xhf[0],xhf[1],xhf[2]
					if (WaveExists(rotationAxis))
						printf ",      rotAxis={%.3g, %.3g, %.3g}",rotationAxis[0],rotationAxis[1],rotationAxis[2]
					endif
					if (calcAxis)
						printf ",      rvec={%.4g, %.4g, %.4g},       rhat={%.4g, %.4g, %.4g}",rvec[0],rvec[1],rvec[2],rhat[0],rhat[1],rhat[2]
					endif
					print ""
				endif

				dot = stringmatch(rotation,"Total") ? 1 : MatrixDot(rotationAxis,Rod)
				sliceWave[i][j] = angle * dot									// amount of rotation about rotationAxis[]

			elseif (ddtensor)													// GND or alpha component, need a dislocation density tensor component, or GND
				curvatureTensorFromRods3d(xhf[0],xhf[1],xhf[2],R3dX,R3dH,R3dF,kappa)
				if (epsilonIsZero)
					GND = dislocationTensorFromRods3d(xhf[0],xhf[1],xhf[2],$"",$"",$"",$"",$"",$"",kappa,alpha)
				else
					GND = dislocationTensorFromRods3d(xhf[0],xhf[1],xhf[2],exx3d,eyy3d,ezz3d,exy3d,exz3d,eyz3d,kappa,alpha)
				endif
				sliceWave[i][j] = (i0<0 || j0<0) ? GND : alpha[i0][j0]	// component of alpha(i0,j0) or GND

			elseif (ctensor)													// kappa, a lattice curvature tensor component
				GND = curvatureTensorFromRods3d(xhf[0],xhf[1],xhf[2],R3dX,R3dH,R3dF,kappa) // this GND assumes no strain
				sliceWave[i][j] = (i0<0 || j0<0) ? GND : kappa[i0][j0]	// component of kappa(i0,j0) or GND (GND, without strain part)

			elseif (strain)
				sliceWave[i][j] = Interp3d(e3d, xhf[0],xhf[1],xhf[2])

			elseif (grad)														// need a part of the gradient of the total rotation angle
				// calculate the gradient of the total rotation angle at the point (xx,hh,ff) from the 3d-Rodriques vectors
				magGrad = GradtRotationAnglesFromRods3d(xhf[0],xhf[1],xhf[2],R3dX,R3dH,R3dF,gradR)
				sliceWave[i][j] = (i0<0) ? magGrad : gradR[i0]			// d|R| / dx,    or dy, or dz   or | gradient (R) |
			endif
		endfor
	endfor
	RodriquesToRGB(sliceOrientsD,"sliceOrients")
	sliceWave2RGB(sliceWave,sliceWaveIntens)
	KillWaves/Z sliceOrientsD, SlicePlaneIn3d_alpha, SlicePlaneIn3d_kappa, SlicePlaneIn3d_grad
	KillWaves/Z rotationAxis,axis,normal,Rod,xhf0,xhf,center,sx,sy,hat,rhat,rvec,W_Cross
	return GetWavesDataFolder(sliceWave,2)
End
//
Static Function xyComponent2ij(str,i,j)			// get i,j component to match something like alpha(xy)
	String str
	Variable &i, &j
	Variable m=strsearch(str,"(",0)
	if (m<0 || strsearch(str[m+3],")",0)<0)		// check for matching (xy)
		i=-1 ; j=-1									// invalid values
		return 1
	endif
	i = WhichListItem(str[m+1],"x;y;z")
	j = WhichListItem(str[m+2],"x;y;z")
	if (i<0 || j<0)
		i=-1 ; j=-1									// invalid values
	endif
	return 0
End
//
// Make a RGB version of slice wave so that you can also show the effect of intensity
Static Function sliceWave2RGB(sliceWave,SliceWaveIntens)
	Wave sliceWave								// 2d wave of angles
	Wave SliceWaveIntens						// 2d wave of intensities
	if (!WaveExists(sliceWave) || !WaveExists(SliceWaveIntens))
		return 1
	endif
	Variable nsx=DimSize(sliceWave,0), nsy=DimSize(sliceWave,1)

	// get angle range
	Variable trim = 0.015						// trim off bottom 1% and top 1% to make colors brighter (not whiter)
	Variable rmin, rmax							// rotation range to map to color table (i.e. saturate at [rmin,rmax]
	rmax = PickSymmetricRange(sliceWave,trim)
	rmin = -rmax
//	printf "symetric,  rmin=%g,  rmax=%g\r", rmin,rmax
//	print "rmax-rmin=",rmax-rmin

	Make/N=(nsx,nsy,2)/O/D sliceWaveRGBd=NaN							// used to make the RGB wave (plane0=color, plane1=intensity)
	sliceWaveRGBd[][][0] = (limit(sliceWave[p][q],rmin,rmax)-rmin)/(rmax-rmin)// color plane is scaled to [0,1]
	sliceWaveRGBd[][][0] = 2*sliceWaveRGBd[p][q][0]-1					// color plane is scaled to [-1,1]

	// get intensity range
	trim = 0.1
	Variable imax = PickSymmetricRange(SliceWaveIntens,trim)				// maximum intensity to use
//	print "imax = ",imax
	sliceWaveRGBd[][][1] = limit(SliceWaveIntens[p][q],0,imax)/imax		// intensity plane is scaled to [0,1]

	// next calculate the RGB values
	Make/N=(nsx,nsy,3)/O/W/U sliceWaveRGB								// RGB wave is 16 bit unsigned (0-65535)
	CopyScales sliceWave, sliceWaveRGB

	String noteStr = note(sliceWave)
	noteStr = ReplaceNumberByKey("rmin",noteStr,rmin,"=")
	noteStr = ReplaceNumberByKey("rmax",noteStr,rmax,"=")
	noteStr = ReplaceStringByKey("waveClass",noteStr,"OrientationSliceWaveRGB","=")
	Note/K sliceWaveRGB, noteStr
	//	r = 65535*(c<=0 ? 1 : (1-c))
	//	g = 65535*abs(1-c)
	//	b = 65535*(c>=0 ? 1 : (1+c))
	sliceWaveRGB[][][0] = round(65535*(sliceWaveRGBd[p][q][0]<=0 ? 1 : (1-sliceWaveRGBd[p][q][0])))	// set the colors
	sliceWaveRGB[][][1] = round(65535*(1-abs(sliceWaveRGBd[p][q][0])))
	sliceWaveRGB[][][2] = round(65535*(sliceWaveRGBd[p][q][0]>=0 ? 1 : (1+sliceWaveRGBd[p][q][0])))
	sliceWaveRGB[][][] *= sliceWaveRGBd[p][q][1]																	// dim to match intensity

	Make/N=3/O/D GreenVec={0,65535,0}			// use Green for un-indexed patterns,  yellow={1,1,0},  magenta={1,0,1},  cyan={0,1,1}
	sliceWaveRGB[][][] = numtype(sliceWaveRGBd[p][q][0]) ? GreenVec[r]*sliceWaveRGBd[p][q][1] : sliceWaveRGB[p][q][r]
	KillWaves/Z GreenVec
	sliceWaveRGB[][][] = numtype(sliceWaveRGBd[p][q][1]) ? GRAY : sliceWaveRGB[p][q][r]					// use Gray for non-existant points
	KillWaves/Z sliceWaveRGBd
End
//
Static Function RodriquesToRGB(sliceOrientsD,RGBimageName)
	Wave sliceOrientsD							// 2-d array of Rodriques vectors
	String RGBimageName						// name of wave to use for rgb image

	Variable zoom=NumVarOrDefault("zoomRGB",1)	// extra zoom for the sliceOrients magnitude scaling

	if (!WaveExists(sliceOrientsD))
		return 1
	endif
	if (stringmatch(NameOfWave(sliceOrientsD),"sliceOrientsD") && strlen(RGBimageName)<1)
		RGBimageName = "sliceOrients"
	endif
	Wave Rod=$MakeUnique3Vector($"")		// direction of Rodriques vector
	Wave rgb=$MakeUnique3Vector($"")		// vector from xhf0[] to point on slice

	Wave vred   =  $MakeUnique3Vector($"")	// unit vector towards (001)
	Wave vgreen= $MakeUnique3Vector($"")	// unit vector towards (011)
	Wave vblue =  $MakeUnique3Vector($"")	// unit vector towards (111)
	vred = {0,0,1}		;	normalize(vred)
	vgreen = {0,1,1}	;	normalize(vgreen)
	vblue = {1,1,1}		;	normalize(vblue)

	vred = {0,0,1}
	vgreen = {0,1,0}
	vblue = {1,0,0}

	Variable i,j, nsx=DimSize(sliceOrientsD,0), nsy=DimSize(sliceOrientsD,1)
	for (j=0;j<nsy;j+=1)
		for (i=0;i<nsx;i+=1)
			Rod = sliceOrientsD[i][j][p]
			if (norm(Rod)<1e-9)
				sliceOrientsD[i][j][] = NaN
				continue
			endif
			rgb = {MatrixDot(vred,Rod),MatrixDot(vgreen,Rod),MatrixDot(vblue,Rod)}
			sliceOrientsD[i][j][] = rgb[r]
		endfor
	endfor

	WaveStats/M=1/Q sliceOrientsD
	Variable NH = 401							// do this to let the top 1% of pixels saturate, makes picture brighter
	Make/N=(NH)/O sliceOrientsD_Hist
	SetScale/I x V_min,V_max,"radian", sliceOrientsD_Hist
	Histogram/B=2 sliceOrientsD,sliceOrientsD_Hist
	Integrate sliceOrientsD_Hist

	Variable rmin, rmax						// trim to only use this range (out of ±65K)
	Variable trim = 0.01						// trim off bottom 1% and top 1% to make colors brighter (not whiter)
	rmin = pnt2x(sliceOrientsD_Hist,ceil(BinarySearchInterp(sliceOrientsD_Hist,trim*sliceOrientsD_Hist[NH-1])))
	rmax = pnt2x(sliceOrientsD_Hist,floor(BinarySearchInterp(sliceOrientsD_Hist,(1-trim)*sliceOrientsD_Hist[NH-1])))
	rmin = limit(rmin,-1,1)
	rmax = limit(rmax,-1,1)
	Variable mr,mg,mb, br,bg,bb

	// zoom about the position of cursor A (or zero if no cursor)
	Variable ci=NaN, cj=NaN					// cursor position on image (points, not scaled values)
	DoWindow GraphSliceOrients 				// checking if GraphSliceOrients exists
	if( V_Flag )
		//	if (strlen(WinList("GraphSliceOrients","","WIN:1"))>1)
		ci = NumberByKey("POINT",CsrInfo(A,"GraphSliceOrients"))
		cj = NumberByKey("YPOINT",CsrInfo(A,"GraphSliceOrients"))
	endif
	if (numtype(ci+cj)==0 && zoom!=1)		// zoom about the cursor color at point [ci,cj]
		printf "sliceOrientsD at cursor = (%g, %g, %g)\r",sliceOrientsD[ci][cj][0],sliceOrientsD[ci][cj][1],sliceOrientsD[ci][cj][2]
		ImageStats/M=1/P=0 sliceOrientsD		// find range of red, and expanded linear eqn.
		zoomRange(V_min,V_max,sliceOrientsD[ci][cj][0],zoom)
		mr=65535/(V_max-V_min)	;	br=-V_min*65535/(V_max-V_min)
		ImageStats/M=1/P=1 sliceOrientsD		// find range of green, and expanded linear eqn.
		mg=65535/(V_max-V_min)	;	bg=-V_min*65535/(V_max-V_min)
		zoomRange(V_min,V_max,sliceOrientsD[ci][cj][1],zoom)
		ImageStats/M=1/P=2 sliceOrientsD		// find range of blue, and expanded linear eqn.
		zoomRange(V_min,V_max,sliceOrientsD[ci][cj][2],zoom)
		mb=65535/(V_max-V_min)	;	bb=-V_min*65535/(V_max-V_min)
	else											// no cursor, zooom about zero (not the middle)
		rmax /= (rmax>0) ? zoom : 1
		rmin /= (rmin<0) ? zoom : 1
		rmin = (rmin>=rmax) ? rmax/2 : rmin	// in case rmin is positive and rmax is reduced to less than rmin
		rmax = (rmax<=rmin) ? rmin/2 : rmax// in case rmax is negative and rmin is reduced to less than rmax
		Variable m=65535/(rmax-rmin), b=-rmin*65535/(rmax-rmin)
		mr=m  ;  mg=m  ;  mb=m  ;  br=b  ;  bg=b  ;  bb=b
		if (zoom!=1)
			printf "rmin=%g,  rmax=%g\r", rmin,rmax
		endif
	endif
	Duplicate/O sliceOrientsD, $RGBimageName
	Wave sliceOrients = $RGBimageName
	Redimension/W/U sliceOrients
//	sliceOrients = limit(m*sliceOrientsD+b,0,65535)
	sliceOrientsD[][][0] = mr*sliceOrientsD[p][q][0] + br
	sliceOrientsD[][][1] = mg*sliceOrientsD[p][q][1] + bg
	sliceOrientsD[][][2] = mb*sliceOrientsD[p][q][2] + bb
	sliceOrients = limit(sliceOrientsD,0,65535)
	sliceOrients[][][] = numtype(sliceOrientsD[p][q][r]) ? GRAY : sliceOrients[p][q][r]	// use Gray for un-indexed patterns
	KillWaves/Z vred,vgreen,vblue, Rod, rgb, sliceOrientsD_Hist
End
Static Function zoomRange(lo,hi,center,zoom)
	Variable &lo,&hi							// range to zoom,  THESE NUMBER ARE CHANGED ON OUTPUT
	Variable center								// value to zoom about
	Variable zoom								// zoom factor
	if (zoom==1)
		return 0
	endif
	if (hi<=lo)									// hi must be greater than lo
		return 1
	endif
	if (lo>center)
		center = lo								// if center is low, move it up to bottom
	elseif (hi<center)
		center = hi								// if center is high, move it down to top
	endif

	Variable dlo=center-lo, dhi=hi-center		// distance from center to both ends
	lo = center - dlo/zoom
	hi = center + dhi/zoom
	return 0
End



// Find the range [lo,hi] along direction vec[] through point v0[] for each corner of a box.
// So (v0[]+lo*n[]), (v0[]+hi*n[]), just span range of box
Static Function FindRangeForBox(xlo,xhi,ylo,yhi,zlo,zhi,vec,v0,lo,hi)
	Variable xlo,xhi,ylo,yhi,zlo,zhi
	Wave vec						// direction of axis
	Wave v0						// location of origin
	Variable &lo, &hi				// the values retrned

	Variable dot
	Wave xyz=$MakeUnique3Vector($"")			// xyz holds the position of each of the box corners
	Wave n=$MakeUnique3Vector($"")				// n is normalized version of vec
	n = vec
	normalize(n)

	lo = Inf	;	hi = -Inf
	// check each of the 8 corners
	xyz = {xlo,ylo,zlo}	;	xyz -= v0	;	dot = MatrixDot(xyz,n)	;	hi = max(hi,dot)	;	lo = min(lo,dot)
	xyz = {xhi,ylo,zlo}	;	xyz -= v0	;	dot = MatrixDot(xyz,n)	;	hi = max(hi,dot)	;	lo = min(lo,dot)
	xyz = {xlo,yhi,zlo}	;	xyz -= v0	;	dot = MatrixDot(xyz,n)	;	hi = max(hi,dot)	;	lo = min(lo,dot)
	xyz = {xhi,yhi,zlo}	;	xyz -= v0	;	dot = MatrixDot(xyz,n)	;	hi = max(hi,dot)	;	lo = min(lo,dot)
	xyz = {xlo,ylo,zhi}	;	xyz -= v0	;	dot = MatrixDot(xyz,n)	;	hi = max(hi,dot)	;	lo = min(lo,dot)
	xyz = {xhi,ylo,zhi}	;	xyz -= v0	;	dot = MatrixDot(xyz,n)	;	hi = max(hi,dot)	;	lo = min(lo,dot)
	xyz = {xlo,yhi,zhi}	;	xyz -= v0	;	dot = MatrixDot(xyz,n)	;	hi = max(hi,dot)	;	lo = min(lo,dot)
	xyz = {xhi,yhi,zhi}	;	xyz -= v0	;	dot = MatrixDot(xyz,n)	;	hi = max(hi,dot)	;	lo = min(lo,dot)
	KillWaves/Z n, xyz
End


// Used to make temporary internal waves, a typical use would be:    Wave xxx=$MakeUnique3Vector($"")
// You will probably want to KillWaves/Z xxx, at the end of the routine.
Function/S MakeUnique3Vector(preset)		// make a new 3-vector, and optionally set it to preset[]
	Wave preset							// optional wave, if none call with $""
	String wName = UniqueName("vec3_",1,0)
	Make/N=3/D $wName
	Wave vec=$wName
	if (WaveExists(preset))
		vec = preset
	endif
	return GetWavesDataFolder(vec,2)		// return full path for assingment to a wave
End

// Used to make temporary internal waves, a typical use would be:    Wave xxx=$MakeUnique3Vector($"")
// You will probably want to KillWaves/Z xxx, at the end of the routine.
Function/S MakeUnique3x3Mat(preset)		// make a new 3x3 matrix, and optionally set it to preset[]
	Wave preset							// optional wave, if none call with $""
	String wName = UniqueName("mat3x3_",1,0)
	Make/N=(3,3)/D $wName
	Wave mat=$wName
	if (WaveExists(preset))
		mat = preset
	endif
	return GetWavesDataFolder(mat,2)		// return full path for assingment to a wave
End



// The big interpolation from random (XYZ)&Rodriques, to three grids of RX, RH, RF
// It calls Interpolate3dRodriquesFromXHF() or Interpolate2dRodriquesFromXHF() as appropriate
Function InterpolateRodriquesFromXHF(resolution,[method])			// choose between the 2d and 3d interpolation
	Variable resolution	// resolution in final array (µm)
	String method

	Wave XX,HH,FF		// positions (beamline coords) for the (RX,RH,RF)
	Wave RX,RH,RF		// three components of the Rodriques vector
	Wave totalAngles		// total rotation angle (degree)
	if (!WaveExists(XX) || !WaveExists(HH) || !WaveExists(FF) || !WaveExists(RX) || !WaveExists(RH) || !WaveExists(RF))
		DoAlert 0, "the waves XX,HH,FF, RX,RH,RF cannot all be found, are you in the right data folder?"
		return 1
	endif
	Variable N=numpnts(RX)
	if (ParamIsDefault(method))
		method = ""
	endif

	Variable defaultResolution = NumVarOrDefault("CutResolution",1.0)
	defaultResolution = defaultResolution<=0 ? NaN : defaultResolution
	resolution = (resolution<=0) ? NaN : resolution			// force 0 or negative to be invalid
	if (numtype(resolution) || resolution<=0)
		Prompt resolution, "resolution in final array (µm)"
		resolution = numtype(resolution) ? defaultResolution : max(abs(resolution),resolution)
		DoPrompt "set resolution", resolution
		if (V_flag)
			return 1
		endif
		if (WhichListItem("Arrays3dButtonProc",GetRTStackInfo(0))>=0)
			printf "•InterpolateRodriquesFromXHF(%g)\r",resolution
		endif
	endif
	if (numtype(resolution) || resolution<=0)
		DoAlert 0, "a resolution of "+num2str(resolution)+" is not valid"
		return 1
	endif

	Variable threshold = 0.18							// (µm) a motion of the PM500 greater than this is intentional (less is jitter)
	Variable Xoff,dX, dimX=1							// for setting scaling in X-dimension
	FindScalingFromVec(XX,threshold,Xoff,dX,dimX)
	if (strlen(method)<1)
		if (dimX==1)
			method = "2d"
		elseif (dimX<=4 )
			method = "sheets"
		elseif (dimX>4)
			method = "3d"
		endif
	endif

	if (stringmatch(method,"2d"))
		print "interpolating a single X-value individually"
		interpolate2dRodriquesFromXHF(resolution)
	elseif (stringmatch(method,"closest"))
		print "not really interpolating by, just taking the closest measured points"
		Interpolate3dRodriquesClosest(resolution)
	elseif (stringmatch(method,"sheets"))
		print "interpolating each X-sheet individually"
		interpolate3dRodFromXHFSheets(resolution)
	elseif (stringmatch(method,"3d"))
		print "interpolating a full 3d array"
		Interpolate3dRodriquesFromXHF(resolution)
	else
		DoAlert 0, "Unable to interpolate, dimX="+num2str(dimX)
	endif
//	if (dimX==1 || stringmatch(method,"2d"))
//		print "interpolating a single X-value individually"
//		interpolate2dRodriquesFromXHF(resolution)
//	elseif (stringmatch(method,"closest"))
//		print "not really interpolating by, just taking the closest measured points"
//		Interpolate3dRodriquesClosest(resolution)
//	elseif (dimX<=4 || stringmatch(method,"sheets"))
//		print "interpolating each X-sheet individually"
//		interpolate3dRodFromXHFSheets(resolution)
//	elseif (dimX>4 || stringmatch(method,"3d"))
//		print "interpolating a full 3d array"
//		Interpolate3dRodriquesFromXHF(resolution)
//	else
//		DoAlert 0, "Unable to interpolate, dimX="+num2str(dimX)
//	endif
End
//Function InterpolateRodriquesFromXHF(resolution)			// choose between the 2d and 3d interpolation
//	Variable resolution	// resolution in final array (µm)
//
//	Wave XX,HH,FF		// positions (beamline coords) for the (RX,RH,RF)
//	Wave RX,RH,RF		// three components of the Rodriques vector
//	Wave totalAngles
//	if (!WaveExists(XX) || !WaveExists(HH) || !WaveExists(FF) || !WaveExists(RX) || !WaveExists(RH) || !WaveExists(RF))
//		DoAlert 0, "the waves XX,HH,FF, RX,RH,RF cannot all be found, are you in the right data folder?"
//		return 1
//	endif
//	Variable N=numpnts(RX)
//
//	Variable defaultResolution = NumVarOrDefault("CutResolution",1.0)
//	defaultResolution = defaultResolution<=0 ? NaN : defaultResolution
//	resolution = (resolution<=0) ? NaN : resolution			// force 0 or negative to be invalid
//	if (numtype(resolution) || resolution<=0)
//		Prompt resolution, "resolution in final array (µm)"
//		resolution = numtype(resolution) ? defaultResolution : max(abs(resolution),resolution)
//		DoPrompt "set resolution", resolution
//		if (V_flag)
//			return 1
//		endif
//	endif
//	if (numtype(resolution) || resolution<=0)
//		DoAlert 0, "a resolution of "+num2str(resolution)+" is not valid"
//		return 1
//	endif
//
//	Variable threshold = 0.18							// (µm) a motion of the PM500 greater than this is intentional (less is jitter)
//	Variable Xoff,dX, dimX=1							// for setting scaling in X-dimension
//	FindScalingFromVec(XX,threshold,Xoff,dX,dimX)
//	if (dimX==1)
//		interpolate2dRodriquesFromXHF(resolution)
//	else
//		Interpolate3dRodriquesFromXHF(resolution)
//	endif
//End
// The big interpolation from random (XYZ)&Rodriques, to three grids of RX, RH, RF
// This has been changed to also interpolate the intensities.
Static Function Interpolate3dRodriquesFromXHF(resolution)
	Variable resolution	// resolution in final array (µm)

	Wave XX,HH,FF		// positions (beamline coords) for the (RX,RH,RF)
	Wave RX,RH,RF		// three components of the Rodriques vector
	Wave totalAngles		// total rotation angle (degree)
	if (!WaveExists(XX) || !WaveExists(HH) || !WaveExists(FF) || !WaveExists(RX) || !WaveExists(RH) || !WaveExists(RF))
		DoAlert 0, "the waves XX,HH,FF, RX,RH,RF cannot all be found, are you in the right data folder?"
		return 1
	endif
	Variable N=numpnts(RX)
	Wave totalPeakIntensity=totalPeakIntensity, totalIntensity=totalIntensity		// these may not exist

	Variable defaultResolution = NumVarOrDefault("CutResolution",1.0)
	defaultResolution = defaultResolution<=0 ? NaN : defaultResolution
	resolution = (resolution<=0) ? NaN : resolution			// force 0 or negative to be invalid
	if (numtype(resolution) || resolution<=0)
		Prompt resolution, "resolution in final array (µm)"
		resolution = numtype(resolution) ? defaultResolution : max(abs(resolution),resolution)
		DoPrompt "set resolution", resolution
		if (V_flag)
			return 1
		endif
	endif
	if (numtype(resolution) || resolution<=0)
		DoAlert 0, "a resolution of "+num2str(resolution)+" is not valid"
		return 1
	endif

	Variable threshold = 0.18							// (µm) a motion of the PM500 greater than this is intentional (less is jitter)
	Variable Xoff,dX, dimX=1							// for setting scaling in X-dimension
	Variable Hoff,dH, dimH=1
	Variable Foff,dF, dimF=1
	FindScalingFromVec(XX,threshold,Xoff,dX,dimX)
	FindScalingFromVec(HH,threshold,Hoff,dH,dimH)
	FindScalingFromVec(FF,threshold,Foff,dF,dimF)
	//	printf "X=[%g,%g], dimX=%d\r",Xoff,Xoff+dX*(dimX-1),dimX
	//	printf "H=[%g,%g], dimH=%d\r",Hoff,Hoff+dH*(dimH-1),dimH
	//	printf "F=[%g,%g], dimF=%d\r",Foff,Foff+dF*(dimF-1),dimF

	Variable Xlo, Hlo, Flo								// low end of 3d volume
	Variable Xhi, Hhi, Fhi								// high end of 3d volume
	Variable NX, NH, NF									// number of points (not intervals) along each axis

	Xlo = resolution*floor(Xoff/resolution)				// low end of 3d volume
	Hlo = resolution*floor(Hoff/resolution)
	Flo = resolution*floor(Foff/resolution)

	NX=ceil((dimX*dX+Xoff-Xlo)/resolution)			// number of points (not intervals) along each axis
	NH=ceil((dimH*dH+Hoff-Hlo)/resolution)
	NF=ceil((dimF*dF+Foff-Flo)/resolution)
	NX = max(NX,1)									// have to have at least one point (no intervals)
	NH = max(NH,1)
	NF = max(NF,1)

	Xhi = Xlo + (NX-1)*resolution						// high end of 3d volume
	Hhi = Hlo + (NH-1)*resolution
	Fhi = Flo + (NF-1)*resolution

	String noteStr = note(RX)
	noteStr = ReplaceNumberByKey("Xlo",noteStr,Xlo,"=")
	noteStr = ReplaceNumberByKey("Xhi",noteStr,Xhi,"=")
	noteStr = ReplaceNumberByKey("Hlo",noteStr,Hlo,"=")
	noteStr = ReplaceNumberByKey("Hhi",noteStr,Hhi,"=")
	noteStr = ReplaceNumberByKey("Flo",noteStr,Flo,"=")
	noteStr = ReplaceNumberByKey("Fhi",noteStr,Fhi,"=")
	noteStr = ReplaceNumberByKey("resolution",noteStr,resolution,"=")
	noteStr = ReplaceStringByKey("waveClass",noteStr,"Interpolated3dArrays","=")

	Variable N3 = NX*NH*NF							// total number of points
	printf "Interpolating Rodriques vectors, with resolution = %gµm\r",resolution
	printf "interpolating on a grid [%g,%g] [%g,%g] [%g,%g], NX=%d,  NH=%d,  NF=%d\r",Xlo,Xhi,Hlo,Hhi,Flo,Fhi,NX,NH,NF
	printf "total number of points in volume = (NX*NH*NF) = %d\r",N3
	Variable predicted = ceil(0.000836*N^1.19 + 3.56e-6*N3* N^1.07)	// predicted execution time (seconds)
	predicted *= numtype(cpuFrequency()) ? 1 : 2e9/cpuFrequency()		// scale to 2GHz
	DoAlert 1,"Interpolating could take ~"+Secs2Time(ceil(predicted),5,0)+",  do it?"
	if (V_flag>1)
		printf "apparently %s was too long to wait for interpolation\r",Secs2Time(ceil(predicted),5,0)
		return 1
	endif
	Variable now = dateTime
	printf "started at %s,   should finish at %s\r",Secs2Time(now,1),Secs2Time(now+predicted,1)
	DoUpdate

	KillWaves/Z totalPeakIntensity3d, TotalIntensity3d
	Make/N=(N,4)/O triangVerts
	triangVerts[][0] = XX[p]
	triangVerts[][1] = HH[p]
	triangVerts[][2] = FF[p]
	Variable timer1=startMSTimer, timer2=startMSTimer
	triangVerts[][0,2] += enoise(resolution * 0.0001)	// need this to keep 3d interpolation from breaking
	Triangulate3d triangVerts
	Wave M_3DVertexList=M_3DVertexList				// vertex list is the same for all three, since the xyz are the same
	//	printf "to triangulate the Rodriques vectors took %s\r",num2sexigesmal(stopMSTimer(timer1)*1e-6,0)
	printf "to triangulate the Rodriques vectors took %s\r",Secs2Time(stopMSTimer(timer1)*1e-6,5,0)

	triangVerts[][3] = RX[p]							// Rodriques X values plus XHF
	Interpolate3D /RNGX={Xlo,resolution,NX}/RNGY={Hlo,resolution,NH}/RNGZ={Flo,resolution,NF}/DEST=R3dX triangulationWave=M_3DVertexList, srcWave=triangVerts
	triangVerts[][3] = RH[p]							// Rodriques H values plus XHF
	Interpolate3D /RNGX={Xlo,resolution,NX}/RNGY={Hlo,resolution,NH}/RNGZ={Flo,resolution,NF}/DEST=R3dH triangulationWave=M_3DVertexList, srcWave=triangVerts
	triangVerts[][3] = RF[p]							// Rodriques F values plus XHF
	Interpolate3D /RNGX={Xlo,resolution,NX}/RNGY={Hlo,resolution,NH}/RNGZ={Flo,resolution,NF}/DEST=R3dF triangulationWave=M_3DVertexList, srcWave=triangVerts
	if (WaveExists(totalPeakIntensity))
		triangVerts[][3] = totalPeakIntensity[p]		// Rodriques total Peak Intensity values plus XHF
		Interpolate3D /RNGX={Xlo,resolution,NX}/RNGY={Hlo,resolution,NH}/RNGZ={Flo,resolution,NF}/DEST=totalPeakIntensity3d triangulationWave=M_3DVertexList, srcWave=triangVerts
	endif
	if (WaveExists(totalIntensity))
		triangVerts[][3] = totalIntensity[p]			// Rodriques total Intensity values plus XHF
		Interpolate3D /RNGX={Xlo,resolution,NX}/RNGY={Hlo,resolution,NH}/RNGZ={Flo,resolution,NF}/DEST=TotalIntensity3d triangulationWave=M_3DVertexList, srcWave=triangVerts
	endif

	SetScale/I x xlo,xhi, "µm", R3dX,R3dH,R3dF		// components of the interpolated Rodriques vectors interpolated into new volume
	SetScale/I y hlo,hhi, "µm", R3dX,R3dH,R3dF
	SetScale/I z flo,fhi, "µm", R3dX,R3dH,R3dF
	Note/K R3dX, noteStr
	Note/K R3dH, noteStr
	Note/K R3dF, noteStr
	if (WaveExists(totalPeakIntensity))
		CopyScales/P R3dX, totalPeakIntensity3d
		Note/K totalPeakIntensity3d, noteStr
	endif
	if (WaveExists(totalIntensity))
		CopyScales/P R3dX, TotalIntensity3d
		Note/K TotalIntensity3d, noteStr
	endif
	Variable/G CutResolution = resolution
//	R3dX = numtype(R3dX[p][q][r]) ? 0 : R3dX[p][q][r]	// set all bad numbers to zero
//	R3dH = numtype(R3dH[p][q][r]) ? 0 : R3dH[p][q][r]
//	R3dF = numtype(R3dF[p][q][r]) ? 0 : R3dF[p][q][r]

	beep
	Variable elapsed = stopMSTimer(timer2)*1e-6
	printf "the whole interpolation process took %s\r",Secs2Time(elapsed,5,0)
	if (elapsed>10*60)										// if execution took more than 10 min, automatically save
		print "This took more than 10min, so save the experiment"
		SaveExperiment
	endif
	KillWaves/Z triangVerts,M_3DVertexList
	return 0
End
//
// The big interpolation from random (XYZ)&Rodriques, to three grids of RX, RH, RF
// this is only a 2-d inerpolation, it assumes that the X is a constant
Static Function Interpolate2dRodriquesFromXHF(resolution)
	Variable resolution	// resolution in final array (µm)

	Wave XX,HH,FF		// positions (beamline coords) for the (RX,RH,RF)
	Wave RX,RH,RF		// three components of the Rodriques vector
	Wave totalAngles		// total rotation angle (degree)
	if (!WaveExists(XX) || !WaveExists(HH) || !WaveExists(FF) || !WaveExists(RX) || !WaveExists(RH) || !WaveExists(RF))
		DoAlert 0, "the waves XX,HH,FF, RX,RH,RF cannot all be found, are you in the right data folder?"
		return 1
	endif
	Variable N=numpnts(RX)

	Variable defaultResolution = NumVarOrDefault("CutResolution",1.0)
	defaultResolution = defaultResolution<=0 ? NaN : defaultResolution
	resolution = (resolution<=0) ? NaN : resolution			// force 0 or negative to be invalid
	if (numtype(resolution) || resolution<=0)
		Prompt resolution, "resolution in final array (µm)"
		resolution = numtype(resolution) ? defaultResolution : max(abs(resolution),resolution)
		DoPrompt "set resolution", resolution
		if (V_flag)
			return 1
		endif
	endif
	if (numtype(resolution) || resolution<=0)
		DoAlert 0, "a resolution of "+num2str(resolution)+" is not valid"
		return 1
	endif

	Variable threshold = 0.18							// (µm) a motion of the PM500 greater than this is intentional (less is jitter)
	Variable Xoff,dX, dimX=1							// for setting scaling in X-dimension
	Variable Hoff,dH, dimH=1
	Variable Foff,dF, dimF=1
	FindScalingFromVec(XX,threshold,Xoff,dX,dimX)
	if (dimX!=1)
		DoAlert 0, "Interpolate2dRodriquesFromXHF() is only for a constant X"
		return 1
	endif
	FindScalingFromVec(HH,threshold,Hoff,dH,dimH)
	FindScalingFromVec(FF,threshold,Foff,dF,dimF)

	Variable Xlo, Hlo, Flo								// low end of 3d volume
	Variable Hhi, Fhi									// high end of 3d volume
	Variable NH, NF										// number of points (not intervals) along each axis

	Xlo = resolution*floor(Xoff/resolution)				// low end of 3d volume
	Hlo = resolution*floor(Hoff/resolution)
	Flo = resolution*floor(Foff/resolution)

	NH=ceil((dimH*dH+Hoff-Hlo)/resolution)			// number of points (not intervals) along each axis
	NF=ceil((dimF*dF+Foff-Flo)/resolution)
	NH = max(NH,1)									// have to have at least one point (no intervals)
	NF = max(NF,1)

	Hhi = Hlo + (NH-1)*resolution						// high end of 3d volumen
	Fhi = Flo + (NF-1)*resolution

	String noteStr = note(RX)
	noteStr = ReplaceNumberByKey("Xlo",noteStr,Xlo,"=")
	noteStr = ReplaceNumberByKey("Xhi",noteStr,Xlo,"=")
	noteStr = ReplaceNumberByKey("Hlo",noteStr,Hlo,"=")
	noteStr = ReplaceNumberByKey("Hhi",noteStr,Hhi,"=")
	noteStr = ReplaceNumberByKey("Flo",noteStr,Flo,"=")
	noteStr = ReplaceNumberByKey("Fhi",noteStr,Fhi,"=")
	noteStr = ReplaceNumberByKey("resolution",noteStr,resolution,"=")

	Variable N3 = NH*NF								// total number of points
	printf "Interpolating Rodriques vectors, with resolution = %gµm\r",resolution
	printf "interpolating on a grid [%g,%g] [%g,%g], NH=%d,  NF=%d\r",Hlo,Hhi,Flo,Fhi,NH,NF
	printf "total number of points in volume = (NH*NF) = %d\r",N3
	Variable predicted = ceil(0.000836*N^1.19 + 3.56e-6*N3* N^1.07)	// predicted execution time (seconds)
	predicted *= numtype(cpuFrequency()) ? 1 : 2e9/cpuFrequency()		// scale to 2GHz
	DoAlert 1,"Interpolating could take ~"+Secs2Time(ceil(predicted),5,0)+",  do it?"
	if (V_flag>1)
		printf "apparently %s was too long to wait for interpolation\r",Secs2Time(ceil(predicted),5,0)
		return 1
	endif
	Variable now = dateTime
	printf "started at %s,   should finish at %s\r",Secs2Time(now,1),Secs2Time(now+predicted,1)
	DoUpdate

	KillWaves/Z R2dX, R2dH, R2dF
	Make/N=(N,3)/O/D myTripletWave
	myTripletWave[][0] = HH[p]
	myTripletWave[][1] = FF[p]
	myTripletWave[][2] = RX[p]

	Variable timer1=startMSTimer, timer2=startMSTimer
	ImageInterpolate/STW/S={0,1,0,1,1,1} Voronoi myTripletWave
	printf "to triangulate the Rodriques vectors took %s\r",Secs2Time(stopMSTimer(timer1)*1e-6,5,0)

	ImageInterpolate/PTW=W_TriangulationData/S={hlo,resolution,hhi,flo,resolution,fhi} Voronoi RX
	Duplicate/O M_InterpolatedImage R2dX
	ImageInterpolate/PTW=W_TriangulationData/S={hlo,resolution,hhi,flo,resolution,fhi} Voronoi RH
	Duplicate/O M_InterpolatedImage R2dH
	ImageInterpolate/PTW=W_TriangulationData/S={hlo,resolution,hhi,flo,resolution,fhi} Voronoi RF
	Duplicate/O M_InterpolatedImage R2dF

	SetScale/I x hlo,hhi, "µm", R2dX,R2dH,R2dF
	SetScale/I y flo,fhi, "µm", R2dX,R2dH,R2dF
	noteStr = ReplaceStringByKey("waveClass",noteStr,"Interpolated3dArrays","=")
	Note/K R2dX, noteStr
	Note/K R2dH, noteStr
	Note/K R2dF, noteStr
	Variable/G CutResolution = resolution
//	R2dX = numtype(R2dX[p][q][r]) ? 0 : R2dX[p][q][r]	// set all bad numbers to zero
//	R2dH = numtype(R2dH[p][q][r]) ? 0 : R2dH[p][q][r]
//	R2dF = numtype(R2dF[p][q][r]) ? 0 : R2dF[p][q][r]

	beep
	Variable elapsed = stopMSTimer(timer2)*1e-6
	printf "the whole interpolation process took %s\r",Secs2Time(elapsed,5,0)
	if (elapsed>10*60)										// if execution took more than 10 min, automatically save
		print "This took more than 10min, so save the experiment"
		SaveExperiment
	endif
	KillWaves/Z W_TriangulationData, myTripletWave, M_InterpolatedImage
	return 0
End
//
Static Function Interpolate3dRodriquesClosest(resolution)
	Variable resolution	// resolution in final array (µm)

	Wave XX,HH,FF		// positions (beamline coords) for the (RX,RH,RF)
	Wave RX,RH,RF		// three components of the Rodriques vector
	Wave totalAngles
	Wave totalIntensity,totalPeakIntensity
	if (!WaveExists(XX) || !WaveExists(HH) || !WaveExists(FF) || !WaveExists(RX) || !WaveExists(RH) || !WaveExists(RF))
		DoAlert 0, "the waves XX,HH,FF, RX,RH,RF cannot all be found, are you in the right data folder?"
		return 1
	endif
	Variable N=numpnts(RX)

	Variable defaultResolution = NumVarOrDefault("CutResolution",1.0)
	defaultResolution = defaultResolution<=0 ? NaN : defaultResolution
	resolution = (resolution<=0) ? NaN : resolution	// force 0 or negative to be invalid
	if (numtype(resolution) || resolution<=0)
		Prompt resolution, "resolution in final array (µm)"
		resolution = numtype(resolution) ? defaultResolution : max(abs(resolution),resolution)
		DoPrompt "set resolution", resolution
		if (V_flag)
			return 1
		endif
	endif
	if (numtype(resolution) || resolution<=0)
		DoAlert 0, "a resolution of "+num2str(resolution)+" is not valid"
		return 1
	endif

	Variable threshold = 0.18							// (µm) a motion of the PM500 greater than this is intentional (less is jitter)
	Variable Xoff,dX, dimX=1								// for setting scaling in X-dimension
	Variable Hoff,dH, dimH=1
	Variable Foff,dF, dimF=1
	FindScalingFromVec(XX,threshold,Xoff,dX,dimX)
	FindScalingFromVec(HH,threshold,Hoff,dH,dimH)
	FindScalingFromVec(FF,threshold,Foff,dF,dimF)
	//	printf "X=[%g,%g], dimX=%d\r",Xoff,Xoff+dX*(dimX-1),dimX
	//	printf "H=[%g,%g], dimH=%d\r",Hoff,Hoff+dH*(dimH-1),dimH
	//	printf "F=[%g,%g], dimF=%d\r",Foff,Foff+dF*(dimF-1),dimF
	Variable Xlo, Hlo, Flo									// low end of 3d volume
	Variable Xhi, Hhi, Fhi									// high end of 3d volume
	Variable NX, NH, NF									// number of points (not intervals) along each axis
	Xlo = resolution*floor(Xoff/resolution)				// low end of 3d volume
	Hlo = resolution*floor(Hoff/resolution)
	Flo = resolution*floor(Foff/resolution)
	NX=ceil((dimX*dX+Xoff-Xlo)/resolution)			// number of points (not intervals) along each axis
	NH=ceil((dimH*dH+Hoff-Hlo)/resolution)
	NF=ceil((dimF*dF+Foff-Flo)/resolution)
	NX = max(NX,1)										// have to have at least one point (no intervals)
	NH = max(NH,1)
	NF = max(NF,1)
	Xhi = Xlo + (NX-1)*resolution						// high end of 3d volume
	Hhi = Hlo + (NH-1)*resolution
	Fhi = Flo + (NF-1)*resolution

	String noteStr = note(RX)
	noteStr = ReplaceNumberByKey("Xlo",noteStr,Xlo,"=")
	noteStr = ReplaceNumberByKey("Xhi",noteStr,Xhi,"=")
	noteStr = ReplaceNumberByKey("Hlo",noteStr,Hlo,"=")
	noteStr = ReplaceNumberByKey("Hhi",noteStr,Hhi,"=")
	noteStr = ReplaceNumberByKey("Flo",noteStr,Flo,"=")
	noteStr = ReplaceNumberByKey("Fhi",noteStr,Fhi,"=")
	noteStr = ReplaceNumberByKey("resolution",noteStr,resolution,"=")
	noteStr = ReplaceStringByKey("waveClass",noteStr,"Interpolated3dArrays","=")

	Variable N3 = NX*NH*NF								// total number of points
	printf "Interpolating Rodriques vectors, with resolution = %gµm\r",resolution
	printf "interpolating on a grid [%g,%g] [%g,%g] [%g,%g], NX=%d,  NH=%d,  NF=%d\r",Xlo,Xhi,Hlo,Hhi,Flo,Fhi,NX,NH,NF
	printf "total number of points in volume = (NX*NH*NF) = %d\r",N3
	Variable predicted = (2.6e-8)*N3*N				// predicted execution time (seconds)
	predicted *= numtype(cpuFrequency()) ? 1 : 2e9/cpuFrequency()		// scale to 2GHz
	if (ItemsInList(GetRTStackInfo(0))<2)
		DoAlert 1,"Interpolating could take ~"+Secs2Time(ceil(predicted),5,0)+",  do it?"
		if (V_flag>1)
			printf "apparently %s was too long to wait for interpolation\r",Secs2Time(ceil(predicted),5,0)
			return 1
		endif
	endif
	Variable now = dateTime
	printf "started at %s,   should finish at %s\r",Secs2Time(now,1),Secs2Time(now+predicted,1)
	DoUpdate
	Variable timer=startMSTimer

	// to improve efficiency, order by H
	Duplicate/O XX,XXsort,index
	Duplicate/O HH,HHsort
	Duplicate/O FF,FFsort
	Duplicate/O XX,index
	index = p
	Sort FFsort ,XXsort,HHsort,FFsort,index
	SetScale/P x 0,1,"",XXsort,HHsort,FFsort
	Variable Ns = numpnts(XXsort), N1,N2

	Make/N=(NX,NH,NF)/O R3dX=NaN,R3dH=NaN,R3dF=NaN, totalPeakIntensity3d=NaN, TotalIntensity3d=NaN
	Variable x,h,f, i,j,k
	Variable m,mbest,dist,minDist
	Variable mf0,mf1, mh0,mh1, mx0,mx1
	Variable timer2=startMSTimer

	for (k=0;k<NF;k+=1)
		//	printf "k=%d of %d\r",k,NF
		if (k==1)
			Variable seconds = stopMSTimer(timer2)*1e-6
			printf "first point took to triangulate the Rodriques vectors took %s\r",Secs2Time(seconds,5,0)
			printf "   the whole thing will take %s\r",Secs2Time(seconds*NF,5,0)
		endif
		f = Flo + k*resolution
		FindLevel/P/Q FFsort, f-3.01*resolution
		mf0 = V_flag ? 0 : floor(V_LevelX)
		FindLevel/P/Q FFsort, f+3.01*resolution
		mf1 = V_flag ? Ns-1 : ceil(V_LevelX)

		Make/N=(mf1-mf0+1)/O XX1sort,HH1sort,FF1sort,index1
		N1 = numpnts(XX1sort)
		XX1sort = XXsort[p+mf0]
		HH1sort = HHsort[p+mf0]
		FF1sort = FFsort[p+mf0]
		index1 = index[p+mf0]
		Sort HH1sort ,XX1sort,HH1sort,FF1sort,index1

		for (j=0;j<NH;j+=1)
			h = Hlo + j*resolution
			FindLevel/P/Q HH1sort, h-3.01*resolution
			mh0 = V_flag ? 0 : floor(V_LevelX)
			FindLevel/P/Q HH1sort, h+3.01*resolution
			mh1 = V_flag ? N1-1 : ceil(V_LevelX)

			Make/N=(mh1-mh0+1)/O XX2sort,HH2sort,FF2sort,index2
			N2 = numpnts(XX2sort)
			XX2sort = XX1sort[p+mh0]
			HH2sort = HH1sort[p+mh0]
			FF2sort = FF1sort[p+mh0]
			index2 = index1[p+mh0]
			Sort XX2sort ,XX2sort,HH2sort,FF2sort,index2

			for (i=0;i<NX;i+=1)
				x = Xlo + i*resolution
				FindLevel/P/Q XX2sort, x-3.01*resolution
				mx0 = V_flag ? 0 : floor(V_LevelX)
				FindLevel/P/Q XX2sort, x+3.01*resolution
				mx1 = V_flag ? N2 : ceil(V_LevelX)

				for (m=mx0,minDist=Inf,mbest=-1; m<=mx1; m+=1)
					if (numtype(RX[index2[m]]))			// skip invalid points
						continue
					endif
					dist = (XX2sort[m]-x)^2 + (HH2sort[m]-h)^2 + (FF2sort[m]-f)^2
					if (dist < minDist)
						minDist = dist
						mbest = m
					endif
				endfor
				minDist = sqrt(minDist)
				if (minDist<3*resolution && mbest>=0)	// close enough and found valid mbest, so aassign the value
					R3dX[i][j][k] = RX[index2[mbest]]
					R3dH[i][j][k] = RH[index2[mbest]]
					R3dF[i][j][k] = RF[index2[mbest]]
					totalPeakIntensity3d[i][j][k] = totalPeakIntensity[index2[mbest]]
					TotalIntensity3d[i][j][k] = totalIntensity[index2[mbest]]
				endif

			endfor
		endfor
	endfor
	KillWaves/Z index,XXsort,HHsort,FFsort,index1,FF1sort,HH1sort,XX1sort,index2,FF2sort,HH2sort,XX2sort
	WaveStats/M=1/Q R3dX
	printf "\rfor a total of %d interpolated points, there are %d (%.1f%%) valid points,   and %d (%.1f%%) invalid points\r",N3,V_npnts,V_npnts/N3*100,V_numNans,V_numNans/N3*100

	SetScale/I x xlo,xhi, "µm", R3dX,R3dH,R3dF, totalPeakIntensity3d, TotalIntensity3d		// components of the interpolated Rodriques vectors interpolated into new volume
	SetScale/I y hlo,hhi, "µm", R3dX,R3dH,R3dF, totalPeakIntensity3d, TotalIntensity3d
	SetScale/I z flo,fhi, "µm", R3dX,R3dH,R3dF, totalPeakIntensity3d, TotalIntensity3d
	Note/K R3dX, noteStr
	Note/K R3dH, noteStr
	Note/K R3dF, noteStr
	Note/K totalPeakIntensity3d, noteStr
	Note/K TotalIntensity3d, noteStr

	Variable/G CutResolution = resolution
	beep
	Variable elapsed = stopMSTimer(timer)*1e-6
	printf "the whole interpolation process took %s\r",Secs2Time(elapsed,5,0)
	if (elapsed>10*60)										// if execution took more than 10 min, automatically save
		print "This took more than 10min, so save the experiment"
		SaveExperiment
	endif
	return 0
End
//
// The big interpolation from random (XYZ)&Rodriques, to three grids of RX, RH, RF
// this is only a 2-d inerpolation followed by a 1d linear interpolation
// this routine assumes that the X coord is the on the OUTER most loop
Static Function interpolate3dRodFromXHFSheets(resolution)
	Variable resolution	// resolution in final array (µm)

	Wave XX=XX,HH=HH,FF=FF		// positions (beamline coords) for the (RX,RH,RF)
	Wave RX=RX,RH=RH,RF=RF		// three components of the Rodriques vector
	Wave totalAngles=totalAngles		// total rotation angle (degree)
	if (!WaveExists(XX) || !WaveExists(HH) || !WaveExists(FF) || !WaveExists(RX) || !WaveExists(RH) || !WaveExists(RF))
		DoAlert 0, "the waves XX,HH,FF, RX,RH,RF cannot all be found, are you in the right data folder?"
		return 1
	endif
	Variable N=numpnts(RX)

	Variable defaultResolution = NumVarOrDefault("CutResolution",1.0)
	defaultResolution = defaultResolution<=0 ? NaN : defaultResolution
	resolution = (resolution<=0) ? NaN : resolution	// force 0 or negative to be invalid
	if (numtype(resolution) || resolution<=0)
		Prompt resolution, "resolution in final array (µm)"
		resolution = numtype(resolution) ? defaultResolution : max(abs(resolution),resolution)
		DoPrompt "set resolution", resolution
		if (V_flag)
			return 1
		endif
	endif
	if (numtype(resolution) || resolution<=0)
		DoAlert 0, "a resolution of "+num2str(resolution)+" is not valid"
		return 1
	endif

	Variable threshold = 0.18							// (µm) a motion of the PM500 greater than this is intentional (less is jitter)
	Variable Xoff,dX, dimX=1								// for setting scaling in X-dimension
	Variable Hoff,dH, dimH=1
	Variable Foff,dF, dimF=1
	FindScalingFromVec(XX,threshold,Xoff,dX,dimX)
	FindScalingFromVec(HH,threshold,Hoff,dH,dimH)
	FindScalingFromVec(FF,threshold,Foff,dF,dimF)
	printf "from raw data in range:  [%g,%g] [%g,%g] [%g,%g], N = (%d, %d, %d),  ∆=(%g, %g, %g)\r",Xoff,Xoff+dX*(dimX-1),Hoff,Hoff+dH*(dimH-1),Foff,Foff+dF*(dimF-1),dimX,dimH,dimF,dX,dH,dF
	//	printf "X=[%g,%g], dimX=%d,  ∆X%.3f\r",Xoff,Xoff+dX*(dimX-1),dimX,dX
	//	printf "H=[%g,%g], dimH=%d,  ∆H=%.3f\r",Hoff,Hoff+dH*(dimH-1),dimH,dH
	//	printf "F=[%g,%g], dimF=%d,  ∆F=%.3f\r",Foff,Foff+dF*(dimF-1),dimF,dF
	Variable Xlo, Hlo, Flo									// low end of 3d volume
	Variable Xhi, Hhi, Fhi									// high end of 3d volume
	Variable NX, NH, NF									// number of points (not intervals) along each axis
	//	Xlo = resolution*floor(Xoff/resolution)			// low end of 3d volume
	Xlo = abs(1-resolution/dX)<=0.05 ? resolution*floor(Xoff/resolution) : Xoff	// low end of 3d volume
	Hlo = resolution*floor(Hoff/resolution)
	Flo = resolution*floor(Foff/resolution)
	NX=ceil((dimX*dX+Xoff-Xlo)/resolution)			// number of points (not intervals) along each axis
	NH=ceil((dimH*dH+Hoff-Hlo)/resolution)
	NF=ceil((dimF*dF+Foff-Flo)/resolution)
	NX = max(NX,1)										// have to have at least one point (no intervals)
	NH = max(NH,1)
	NF = max(NF,1)
	Xhi = Xlo + (NX-1)*resolution						// high end of 3d volume
	Hhi = Hlo + (NH-1)*resolution
	Fhi = Flo + (NF-1)*resolution

	String noteStr = note(RX)
	noteStr = ReplaceNumberByKey("Xlo",noteStr,Xlo,"=")
	noteStr = ReplaceNumberByKey("Xhi",noteStr,Xhi,"=")
	noteStr = ReplaceNumberByKey("Hlo",noteStr,Hlo,"=")
	noteStr = ReplaceNumberByKey("Hhi",noteStr,Hhi,"=")
	noteStr = ReplaceNumberByKey("Flo",noteStr,Flo,"=")
	noteStr = ReplaceNumberByKey("Fhi",noteStr,Fhi,"=")
	noteStr = ReplaceNumberByKey("resolution",noteStr,resolution,"=")
	noteStr = ReplaceStringByKey("waveClass",noteStr,"Interpolated3dArrays","=")

	Variable N3 = NX*NH*NF								// total number of points
	printf "Interpolating Rodriques vectors, with resolution = %gµm\r",resolution
	printf "interpolating on a grid [%g,%g] [%g,%g] [%g,%g], NX=%d,  NH=%d,  NF=%d\r",Xlo,Xhi,Hlo,Hhi,Flo,Fhi,NX,NH,NF
	printf "total number of points in volume = (NX*NH*NF) = %d\r",N3
	Variable predicted = ceil(N*0.00268)				// predicted execution time (seconds)
	predicted *= numtype(cpuFrequency()) ? 1 : 2e9/cpuFrequency()		// scale to 2GHz
	DoAlert 1,"Interpolating should take ~"+Secs2Time(ceil(predicted),5,0)+",  do it?"
	if (V_flag>1)
		printf "apparently %s was too long to wait for interpolation\r",Secs2Time(ceil(predicted),5,0)
		return 1
	endif
	String progressWin = ProgressPanelStart("",stop=1,showTime=1)		// display a progress bar
	Variable now = dateTime
	printf "started at %s,   should finish at %s\r",Secs2Time(now,1),Secs2Time(now+predicted,1)
	DoUpdate

	Variable timer=startMSTimer
	Make/N=(NX,NH,NF)/O R3dX,R3dH,R3dF
	SetScale/I x xlo,xhi, "µm", R3dX,R3dH,R3dF		// components of the interpolated Rodriques vectors interpolated into new volume
	SetScale/I y hlo,hhi, "µm", R3dX,R3dH,R3dF
	SetScale/I z flo,fhi, "µm", R3dX,R3dH,R3dF
	Note/K R3dX, noteStr
	Note/K R3dH, noteStr
	Note/K R3dF, noteStr
	Wave totalPeakIntensity=totalPeakIntensity, totalIntensity=totalIntensity
	KillWaves/Z totalPeakIntensity3d, TotalIntensity3d
	if (WaveExists(totalPeakIntensity))
		Make/N=(NX,NH,NF)/O totalPeakIntensity3d
		SetScale/I x xlo,xhi, "µm", totalPeakIntensity3d
		SetScale/I y hlo,hhi, "µm", totalPeakIntensity3d
		SetScale/I z flo,fhi, "µm", totalPeakIntensity3d
		Note/K totalPeakIntensity3d, noteStr
	endif
	if (WaveExists(totalIntensity))
		Make/N=(NX,NH,NF)/O TotalIntensity3d
		SetScale/I x xlo,xhi, "µm", TotalIntensity3d
		SetScale/I y hlo,hhi, "µm", TotalIntensity3d
		SetScale/I z flo,fhi, "µm", TotalIntensity3d
		Note/K TotalIntensity3d, noteStr
	endif
	Wave exx=exx, eyy=eyy, ezz=ezz, exy=exy, exz=exz, eyz=eyz
	Wave a_nm=a_nm, b_nm=b_nm, c_nm=c_nm, alpha=alpha, bet=bet, gam=gam
	Variable fepsilon = WaveExists(exx)				// flags whether epsilon is present
	if (fepsilon)
		Make/N=(NX,NH,NF)/O exx3d, eyy3d, ezz3d, exy3d, exz3d, eyz3d
		SetScale/I x xlo,xhi, "µm",exx3d, eyy3d, ezz3d, exy3d, exz3d, eyz3d
		SetScale/I y hlo,hhi, "µm", exx3d, eyy3d, ezz3d, exy3d, exz3d, eyz3d
		SetScale/I z flo,fhi, "µm", exx3d, eyy3d, ezz3d, exy3d, exz3d, eyz3d
		Note/K exx3d, noteStr;	Note/K eyy3d, noteStr;	Note/K ezz3d, noteStr
		Note/K exy3d, noteStr;	Note/K exz3d, noteStr;	Note/K eyz3d, noteStr

		Make/N=(NX,NH,NF)/O a3d, b3d, c3d, alpha3d, beta3d, gamma3d
		SetScale/I x xlo,xhi, "µm",a3d, b3d, c3d, alpha3d, beta3d, gamma3d
		SetScale/I y hlo,hhi, "µm", a3d, b3d, c3d, alpha3d, beta3d, gamma3d
		SetScale/I z flo,fhi, "µm", a3d, b3d, c3d, alpha3d, beta3d, gamma3d
		Note/K a3d, noteStr;		Note/K b3d, noteStr;		Note/K c3d, noteStr
		Note/K alpha3d, noteStr;	Note/K beta3d, noteStr;	Note/K gamma3d, noteStr
	endif

	// find a list of X-values in the raw data
	Make/N=(N)/O/D Xvalues=NaN
	Variable i,m, xval, nXval=0
	for (m=0,nXval=0; m<N; m+=1)
		if (numtype(XX[m]))
			continue
		endif
		xval = XX[m]
		for (i=0;i<nXval;i+=1)							// check if xval is different than existing values in Xvalues[]
			if (abs(Xvalues[i]-xval)<threshold)
				break
			endif
		endfor
		if (i>=nXval)
			Xvalues[nXval] = xval
			nXval += 1
		endif
	endfor
	Sort Xvalues, Xvalues
	Redimension/N=(nXval) Xvalues						// sorted (increasing) list of valid measured X-values

	Make/N=(NH,NF)/O RXlo, RXhi						// planes of RX, interpolate between these two planes
	Make/N=(NH,NF)/O RHlo, RHhi						//  RH, etc.
	Make/N=(NH,NF)/O RFlo, RFhi
	if (WaveExists(totalPeakIntensity))
		Make/N=(NH,NF)/O totalPeakIntensityLo, totalPeakIntensityHi
	endif
	if (WaveExists(totalIntensity))
		Make/N=(NH,NF)/O totalIntensityLo, totalIntensityHi
	endif
	if (fepsilon)
		Make/N=(NH,NF)/O exxLo, exxHi, eyyLo, eyyHi, ezzLo, ezzHi
		Make/N=(NH,NF)/O exyLo, exyHi, exzLo, exzHi, eyzLo, eyzHi
		Make/N=(NH,NF)/O aLo, aHi, bLo, bHi, cLo, cHi
		Make/N=(NH,NF)/O alphaLo, alphaHi, betaLo, betaHi, gammaLo, gammaHi
	endif

	Variable hiX=Xvalues[0]								// current X value of the to measured planes
	Variable ihiX=0										// index into Xvalues[] of current hiX
	Make/N=(N)/O/D indexWave=NaN					// only used in this routine
	Variable ix, mx
	for (i=0,mx=0;i<N;i+=1)							// find all of the raw data points with this X, and save location in indexWave[]
//		if (abs(hiX-XX[i])<threshold)
		if (abs(hiX-XX[i])<threshold && numtype(RX[i])==0)
			indexWave[mx] = i
			mx += 1
		endif
	endfor
	if (mx<4)
		print "\r\t\tnot enough points to define a sheet, need at least 4 at this X (the first one)"
		return 1
	endif

	String indexNote = ""
	indexNote = ReplaceNumberByKey("mx",indexNote,mx,"=")
	indexNote = ReplaceNumberByKey("hlo",indexNote,hlo,"=")
	indexNote = ReplaceNumberByKey("hhi",indexNote,hhi,"=")
	indexNote = ReplaceNumberByKey("flo",indexNote,flo,"=")
	indexNote = ReplaceNumberByKey("fhi",indexNote,fhi,"=")
	indexNote = ReplaceNumberByKey("resolution",indexNote,resolution,"=")
	indexNote = ReplaceNumberByKey("Xval",indexNote,hiX,"=")
	Note/K indexWave, indexNote

	Make/N=(mx,3)/O/D myTripletWave				// interpolate all sheets at this measured X = hiX
	myTripletWave[][2] = 0
	myTripletWave[][0] = HH[indexWave[p]]
	myTripletWave[][1] = FF[indexWave[p]]
	ImageInterpolate/STW/S={0,resolution,0,1,resolution,1} Voronoi myTripletWave	// only calculates W_TriangulationData
	// reversed order for the first one, only, also copy the wave note
	interpOneDataSheetValue(RX,indexWave,RXhi) ;				CopyWave(RXlo,RXhi)
	interpOneDataSheetValue(RH,indexWave,RHhi) ;				CopyWave(RHlo,RHhi)
	interpOneDataSheetValue(RF,indexWave,RFhi) ;				CopyWave(RFlo,RFhi)
	if (WaveExists(totalPeakIntensity))
		interpOneDataSheetValue(totalPeakIntensity,indexWave,totalPeakIntensityHi)
		CopyWave(totalPeakIntensityLo,totalPeakIntensityHI)
	endif
	if (WaveExists(totalIntensity))
		interpOneDataSheetValue(totalIntensity,indexWave,totalIntensityHi)
		CopyWave(totalIntensityLo,totalIntensityHi)
	endif
	if (fepsilon)
		interpOneDataSheetValue(exx,indexWave,exxHi) ;		CopyWave(exxLo,exxHi)
		interpOneDataSheetValue(eyy,indexWave,eyyHi) ;		CopyWave(eyyLo,eyyHi)
		interpOneDataSheetValue(ezz,indexWave,ezzHi) ;		CopyWave(ezzLo,ezzHi)
		interpOneDataSheetValue(exy,indexWave,exyHi) ;		CopyWave(exxLo,exxHi)
		interpOneDataSheetValue(exz,indexWave,exzHi) ; 		CopyWave(exxLo,exxHi)
		interpOneDataSheetValue(eyz,indexWave,eyzHi) ;		CopyWave(exxLo,exxHi)

		interpOneDataSheetValue(a_nm,indexWave,aHi) ;		CopyWave(aLo,aHi)
		interpOneDataSheetValue(b_nm,indexWave,bHi) ;		CopyWave(bLo,bHi)
		interpOneDataSheetValue(c_nm,indexWave,cHi) ;		CopyWave(cLo,cHi)
		interpOneDataSheetValue(alpha,indexWave,alphaHi);	CopyWave(alphaLo,alphaHi)
		interpOneDataSheetValue(bet,indexWave,betaHi)	 ;		CopyWave(betaLo,betaHi)
		interpOneDataSheetValue(gam,indexWave,gammaHi) ;	CopyWave(gammaLo,gammaHi)
	endif

	for (ix=0;ix<NX;ix+=1)								// for each of the X-planes in R3dX,R3dH,R3dF, ...
		if (ProgressPanelUpdate(progressWin,ix/Nx*100	))	// update progress bar
			break											//   and break out of loop
		endif
		xval = xlo + ix*DimDelta(R3dX,0)				// x value for the ix plane
		if (xval > hiX)										// need to increment data sheets (set hi=lo, and compte new hi)
			ihiX += 1										// set new hiX
			hiX = Xvalues[ihiX]
			for (i=0,mx=0;i<N;i+=1)					// find all of the raw data points with this X, and save location in indexWave[]
//				if (abs(hiX-XX[i])<threshold)
				if (abs(hiX-XX[i])<threshold && numtype(RX[i])==0)
					indexWave[mx] = i
					mx += 1
				endif
			endfor
			if (mx<4)
				print "\r\t\tnot enough points to define a sheet, need at least 4 at X = ",hiX
				return 1
			endif
			indexNote = ReplaceNumberByKey("mx",indexNote,mx,"=")
			indexNote = ReplaceNumberByKey("Xval",indexNote,hiX,"=")
			Note/K indexWave, indexNote
			Make/N=(mx,3)/O/D myTripletWave		// interpolate all sheets at this measured X = hiX
			myTripletWave[][2] = 0
			myTripletWave[][0] = HH[indexWave[p]]
			myTripletWave[][1] = FF[indexWave[p]]
			ImageInterpolate/STW/S={0,resolution,0,1,resolution,1} Voronoi myTripletWave	// re-calculates W_TriangulationData
			CopyWave(RXlo,RXhi);	interpOneDataSheetValue(RX,indexWave,RXhi)
			CopyWave(RHlo,RHhi);	interpOneDataSheetValue(RH,indexWave,RHhi)
			CopyWave(RFlo,RFhi);	interpOneDataSheetValue(RF,indexWave,RFhi)
			if (WaveExists(totalPeakIntensity))
				CopyWave(totalPeakIntensityLo,totalPeakIntensityHi)
				interpOneDataSheetValue(totalPeakIntensity,indexWave,totalPeakIntensityHi)
			endif
			if (WaveExists(totalIntensity))
				CopyWave(totalIntensityLo,totalIntensityHi)
				interpOneDataSheetValue(totalIntensity,indexWave,totalIntensityHi)
			endif
			if (fepsilon)
				CopyWave(exxLo,exxHi) ;		CopyWave(eyyLo,eyyHi) ;	CopyWave(ezzLo,ezzHi)
				CopyWave(exyLo,exyHi);			CopyWave(exzLo, exzHi);	CopyWave(eyzLo,eyzHi)
				CopyWave(aLo,aHi) ;				CopyWave(bLo,bHi) ;			CopyWave(cLo,cHi)
				CopyWave(alphaLo,alphaHi) ;	CopyWave(betaLo,betaHi) ;	CopyWave(gammaLo,gammaHi)
				interpOneDataSheetValue(exx,indexWave,exxHi)
				interpOneDataSheetValue( eyy,indexWave,eyyHi)
				interpOneDataSheetValue(ezz,indexWave,ezzHi)
				interpOneDataSheetValue(exy,indexWave,exyHi)
				interpOneDataSheetValue(exz,indexWave,exzHi)
				interpOneDataSheetValue(eyz,indexWave,eyzHi)
				interpOneDataSheetValue(a_nm,indexWave,aHi)
				interpOneDataSheetValue(b_nm,indexWave,bHi)
				interpOneDataSheetValue(c_nm,indexWave,cHi)
				interpOneDataSheetValue(alpha,indexWave,alphaHi)
				interpOneDataSheetValue(bet,indexWave,betaHi)
				interpOneDataSheetValue(gam,indexWave,gammaHi)
			endif
		endif												// have all of the lo and hi sheets

		// now actually interp the 3d arrays
		interpBetweenTwoSheetsInto3d(RXlo,RXhi,xval,R3dX,ix)
		interpBetweenTwoSheetsInto3d(RHlo,RHhi,xval,R3dH,ix)
		interpBetweenTwoSheetsInto3d(RFlo,RFhi,xval,R3dF,ix)
		if (WaveExists(totalPeakIntensity))
			interpBetweenTwoSheetsInto3d(totalPeakIntensityLo,totalPeakIntensityHi,xval,totalPeakIntensity3d,ix)
		endif
		if (WaveExists(totalIntensity))
			interpBetweenTwoSheetsInto3d(totalIntensityLo,totalIntensityHi,xval,TotalIntensity3d,ix)
		endif
		if (fepsilon)
			interpBetweenTwoSheetsInto3d(exxLo,exxHi,xval,exx3d,ix)
			interpBetweenTwoSheetsInto3d(eyyLo,eyyHi,xval,eyy3d,ix)
			interpBetweenTwoSheetsInto3d(ezzLo,ezzHi,xval,ezz3d,ix)
			interpBetweenTwoSheetsInto3d(exyLo,exyHi,xval,exy3d,ix)
			interpBetweenTwoSheetsInto3d(exzLo,exzHi,xval,exz3d,ix)
			interpBetweenTwoSheetsInto3d(eyzLo,eyzHi,xval,eyz3d,ix)
			interpBetweenTwoSheetsInto3d(aLo,aHi,xval,a3d,ix)
			interpBetweenTwoSheetsInto3d(bLo,bHi,xval,b3d,ix)
			interpBetweenTwoSheetsInto3d(cLo,cHi,xval,c3d,ix)
			interpBetweenTwoSheetsInto3d(alphaLo,alphaHi,xval,alpha3d,ix)
			interpBetweenTwoSheetsInto3d(betaLo,betaHi,xval,beta3d,ix)
			interpBetweenTwoSheetsInto3d(gammaLo,gammaHi,xval,gamma3d,ix)
		endif
	endfor

	Variable/G CutResolution = resolution
	Variable elapsed = stopMSTimer(timer)*1e-6
	beep
	printf "the whole interpolation process took %s\r",Secs2Time(elapsed,5,0)
	if (elapsed>10*60)									// if execution took more than 10 min, automatically save
		print "This took more than 10min, so save the experiment"
		SaveExperiment
	endif
	DoWindow/K $progressWin

	KillWaves/Z Xvalues
	KillWaves/Z W_TriangulationData, myTripletWave, M_InterpolatedImage, indexWave
	KillWaves/Z RXlo, RXhi, RHlo, RHhi, RFlo, RFhi
	KillWaves/Z totalPeakIntensityLo, totalPeakIntensityHi
	KillWaves/Z totalIntensityLo, totalIntensityHi
	KillWaves/Z exxLo, exxHi, eyyLo, eyyHi, ezzLo, ezzHi, exyLo, exyHi, exzLo, exzHi, eyzLo, eyzHi
	KillWaves/Z aLo, aHi, bLo, bHi, cLo, cHi, alphaLo, alphaHi, betaLo, betaHi, gammaLo, gammaHi
	return 0
End
//
Static Function interpOneDataSheetValue(valWave,indexWave,surf)
	Wave valWave								// such as RF
	Wave indexWave
	Wave surf									// resultant surface

	String indexNote = note(indexWave)
	Variable resolution, hlo,hhi,flo,fhi, mx
	resolution = NumberByKey("resolution", indexNote,"=")
	hlo = NumberByKey("hlo", indexNote,"=")
	hhi = NumberByKey("hhi", indexNote,"=")
	flo = NumberByKey("flo", indexNote,"=")
	fhi = NumberByKey("fhi", indexNote,"=")
	mx = NumberByKey("mx", indexNote,"=")
	Make/N=(mx)/O/D tempVal

	tempVal = valWave[indexWave[p]]
	Wave W_TriangulationData=W_TriangulationData
	ImageInterpolate/PTW=W_TriangulationData/S={hlo,resolution,hhi+resolution,flo,resolution,fhi+resolution} Voronoi tempVal
	Wave M_InterpolatedImage=M_InterpolatedImage
	surf = M_InterpolatedImage[p][q]
	Note/K surf, indexNote
	CopyScales/P M_InterpolatedImage, surf
	KillWaves/Z tempVal
End
//
Static Function CopyWave(dest,source)		// copy a wave, and its wave note
	Wave dest,source
	dest = source
	Note/K dest, note(source)
	CopyScales/P source, dest
End
//
Static Function interpBetweenTwoSheetsInto3d(sheetLo,sheetHi,xval,r3d,ix)
	Wave sheetLo, sheetHi
	Variable xval
	Wave r3d
	Variable ix								// first index of r3d (loop over others)

	// interpolate at xval, which should be in range (xlo,xhi)
	Variable xlo = NumberByKey("Xval", note(sheetLo),"=")
	Variable xhi = NumberByKey("Xval", note(sheetHi),"=")
	Variable fLo, fHi							// fractions of Lo and Hi sheets
	if (xhi==xlo)
		fLo = 1
		fhi = 0
	else
		fHI = (xval-xlo)/(xhi-xlo)
		fHi = limit(fHi,0,1)
		fLo = 1-fHi
	endif
	r3d[ix][][] = fLo*sheetLo[q][r] + fHi*sheetHi[q][r]
End
//Static Function interpBetweenTwoSheets(sheetLo,sheetHi,xval,result)
//	Wave sheetLo, sheetHi
//	Variable xval
//	Wave result
//
//	// interpolate at xval, which should be in range (xlo,xhi)
//	Variable xlo = NumberByKey("Xval", note(sheetLo),"=")
//	Variable xhi = NumberByKey("Xval", note(sheetHi),"=")
//	Variable fLo, fHi				// fractions of Lo and Hi sheets
//	fLo = (xval-xlo)/(xhi-xlo)
//	fLo = limit(fLo,0,1)
//	fHi = 1-fLo
//	result = fLo*sheetLo + fHi*sheetHi
//End




// makes pretend test data for testing
Function MakeTestData(type)
	String type
	Prompt type,"type of test data",popup,"X;H;F;circulating;toroidal"
	DoPrompt "type",type
	if (V_flag)
		return 1
	endif
	printf "   MakeTestData(\"%s\")\r",type
	NewDataFolder/O/S root:testData
	Make/N=(51,51,51)/O R3dX,R3dH,R3dF

	SetScale/I x -25,25,"µm", R3dF,R3dH,R3dX
	SetScale/I y -25,25,"µm", R3dF,R3dH,R3dX
	SetScale/I z 0,50,"µm", R3dF,R3dH,R3dX
	//	String noteStr = "X0=0;H0=0;F0=0;dX=0;dH=0;dF=-1;Xlo=-25;Xhi=25;Hlo=-25;Hhi=25;Fhi=50;Flo=0;resolution=1;"
	String noteStr = "X0=0;H0=0;F0=0;dX=0;dH=0;dF=1;Xlo=-25;Xhi=25;Hlo=-25;Hhi=25;Fhi=50;Flo=0;resolution=1;"
	noteStr = ReplaceStringByKey("waveClass",noteStr,"Interpolated3dArrays","=")
	Note/K R3dX	;	Note R3dX, noteStr
	Note/K R3dH	;	Note R3dH, noteStr
	Note/K R3dX	;	Note R3dX, noteStr
	Note/K R3dF	;	Note R3dF, noteStr

	Variable angle = tan(.5)
	Duplicate/O R3dX, tempPhi, tempDist
	Variable f0=20,r0=10
	tempPhi = atan2(y,x)
	tempDist = sqrt((x-r0*cos(tempPhi))^2+(y-r0*sin(tempPhi))^2 + (z-f0)^2)

	strswitch(type)
		case "X":
			R3dX = angle
			R3dH = 0
			R3dF = 0
			break
		case "H":
			R3dX = 0
			R3dH = angle
			R3dF = 0
			break
		case "F":
			R3dX = 0
			R3dH = 0
			R3dF = angle
			break
		case "circulating":
			R3dX = 0
			R3dH = 0
			R3dF = angle*exp(-tempDist/7)
			break
		case "toroidal":
			R3dX = -angle*sin(tempPhi)*exp(-tempDist/7)				// direction of R is same as theta^
			R3dH = angle*cos(tempPhi)*exp(-tempDist/7)
			R3dF = 0
			break
		default:
			Abort "illegal type, '"+type+"'"
	endswitch
	KillWaves/Z tempPhi, tempDist
	DoAlert 0, "Current data folder is now 'root:testData:'"
End



// make the smallest box that contains all of the good points, then delete all points outside the box
Function TrimBadPointsOnOutside()
	Wave FF=FF, HH=HH, XX=XX
	Wave RF=RF, RH=RH, RX=RX
	Wave totalAngles=totalAngles, totalIntensity=totalIntensity, totalPeakIntensity=totalPeakIntensity
	if (!WaveExists(XX) || !WaveExists(HH) || !WaveExists(FF) || !WaveExists(RX) || !WaveExists(RH) || !WaveExists(RF) || !WaveExists(totalAngles) || !WaveExists(totalIntensity))
		return 1
	endif

	Duplicate/O XX posTemp
	posTemp = numtype(RX[p]) ? NaN : FF[p]
	WaveStats/M=1/Q posTemp
	Variable Flo=V_min, Fhi=V_max
	posTemp = numtype(RX[p]) ? NaN : HH[p]
	WaveStats/M=1/Q posTemp
	Variable Hlo=V_min, Hhi=V_max
	posTemp = numtype(RX[p]) ? NaN : XX[p]
	WaveStats/M=1/Q posTemp
	Variable Xlo=V_min, Xhi=V_max
	KillWaves/Z posTemp

	Variable i,N=numpnts(XX)
	for (i=N-1;i>=0;i-=1)
		if (numtype(RX[i]) || XX[i]<Xlo || XX[i]>Xhi  ||  HH[i]<Hlo || HH[i]>Hhi  ||  FF[i]<Flo || FF[i]>Fhi)
			DeletePoints i,1, totalAngles,totalIntensity,RF,RH,RX,FF,HH,XX
			if (WaveExists(totalPeakIntensity))
				DeletePoints i,1, totalPeakIntensity
			endif
			if (WaveExists(imageNames))
				DeletePoints i,1, imageNames
			endif
		endif
	endfor
	print "final length is ",numpnts(XX)

	String fldr = GetWavesDataFolder(XX,1)
	Note/K XX,ReplaceStringByKey("fldrName",note(xx),fldr,"=")
	Note/K HH,ReplaceStringByKey("fldrName",note(HH),fldr,"=")
	Note/K FF,ReplaceStringByKey("fldrName",note(FF),fldr,"=")
	Note/K RX,ReplaceStringByKey("fldrName",note(RX),fldr,"=")
	Note/K RH,ReplaceStringByKey("fldrName",note(RH),fldr,"=")
	Note/K RF,ReplaceStringByKey("fldrName",note(RF),fldr,"=")
	Note/K totalAngles,ReplaceStringByKey("fldrName",note(totalAngles),fldr,"=")
	Note/K totalIntensity,ReplaceStringByKey("fldrName",note(totalIntensity),fldr,"=")
	if (WaveExists(totalPeakIntensity))
		Note/K totalPeakIntensity,ReplaceStringByKey("fldrName",note(totalPeakIntensity),fldr,"=")
	endif
End
//
Function removeAllBadPoints()		// just delete all points that have invalid rotation vectors (all bad points)
	Wave FF=FF, HH=HH, XX=XX
	Wave RF=RF, RH=RH, RX=RX
	Wave totalAngles=totalAngles, totalIntensity=totalIntensity, totalPeakIntensity=totalPeakIntensity
	Variable i,N=numpnts(XX)
	for (i=N-1;i>=0;i-=1)
		if (numtype(RX[i]))
			DeletePoints i,1, totalAngles,totalIntensity,totalPeakIntensity,RF,RH,RX,FF,HH,XX
			if (WaveExists(imageNames))
				DeletePoints i,1, imageNames
			endif
		endif
	endfor

	String fldr = GetWavesDataFolder(XX,1)
	Note/K XX,ReplaceStringByKey("fldrName",note(xx),fldr,"=")
	Note/K HH,ReplaceStringByKey("fldrName",note(HH),fldr,"=")
	Note/K FF,ReplaceStringByKey("fldrName",note(FF),fldr,"=")
	Note/K RX,ReplaceStringByKey("fldrName",note(RX),fldr,"=")
	Note/K RH,ReplaceStringByKey("fldrName",note(RH),fldr,"=")
	Note/K RF,ReplaceStringByKey("fldrName",note(RF),fldr,"=")
	Note/K totalAngles,ReplaceStringByKey("fldrName",note(totalAngles),fldr,"=")
	Note/K totalIntensity,ReplaceStringByKey("fldrName",note(totalIntensity),fldr,"=")
	Note/K totalPeakIntensity,ReplaceStringByKey("fldrName",note(totalPeakIntensity),fldr,"=")
End



// start of readin process
Function/T readInXYZorientationsFromFile(maxAngle)
	Variable maxAngle
	if (!(maxAngle>0))
		maxAngle = !(maxAngle>0) ? Inf : maxAngle
		Prompt maxAngle, "reject all points with rotation angle greater than this (degree)"
		DoPrompt "max Rotation",maxAngle
		if (V_flag)
			return ""
		endif
		if (WhichListItem("Arrays3dButtonProc",GetRTStackInfo(0))>=0)
			printf "•readInXYZorientationsFromFile(%g)\r",maxAngle
		endif
	endif
	if (!(maxAngle>0))
		DoAlert 0, "ERROR in readInXYZorientationsFromFile(), maxAngle = "+num2str(maxAngle)
		return ""
	endif

	String/G root:rawDataFullPathName
	SVAR rawDataFullPathName=root:rawDataFullPathName
	Variable refNum

	Open/R/Z=1 refNum rawDataFullPathName	// see if rawDataFullPathName exists?
	if (refNum>0)
		Close refNum
	endif

	if (!V_flag)									// V_flag==0 if file exists
		DoAlert 1, "use: "+rawDataFullPathName
		V_flag = (V_flag==2)
	endif
	if (V_flag)
		rawDataFullPathName = ""
	endif
	Open/R/M="Wirescan peak orientations"/P=home/Z=2 refNum rawDataFullPathName
	if (refNum<=0)
		Abort "unable to re-open the data file"
	endif
	Close refNum
	rawDataFullPathName = S_fileName			// save file name for future use

	String noteStr=""
	String list = getListOfTypesInFile(rawDataFullPathName,"")				// get list of file types

	if (WhichListItem("ArrayOf3dRotationsFile",list)>=0)
		noteStr = readInRodriquesDataFile(rawDataFullPathName,maxAngle)	//	Rodriquez vectors (rotations)
	elseif (WhichListItem("ArrayOf3dOrientsFile",list)>=0)
		noteStr = read3dRecipLatticesFile(rawDataFullPathName,maxAngle)	//	reciprocal lattices
//		noteStr = readInNewRawDataFile(rawDataFullPathName,maxAngle)	//	reciprocal lattices
	else
		Abort "You have chosen an unknown type of file, it is not supported"
	endif

	if (strlen(noteStr)<1)
		DoAlert 0, "Unable to properly read "+rawDataFullPathName
		return ""
	endif

	DuplicateTheXs()								// change data with only one X value to have two X values
	return noteStr
End
//// start of readin process
//Function/T readInXYZorientationsFromFile()
//	String/G root:rawDataFullPathName
//	SVAR rawDataFullPathName=root:rawDataFullPathName
//	Variable refNum
//
//	Open/R/Z=1 refNum rawDataFullPathName	// see if rawDataFullPathName exists?
//	if (refNum>0)
//		Close refNum
//	endif
//
//	if (!V_flag)								// V_flag==0 if file exists
//		DoAlert 1, "use: "+rawDataFullPathName
//		V_flag = (V_flag==2)
//	endif
//	if (V_flag)
//		rawDataFullPathName = ""
//	endif
//	Open/R/M="Wirescan peak positions from Wenge"/P=home/Z=2 refNum rawDataFullPathName
//	if (refNum<=0)
//		Abort "unable to re-open the data file"
//	endif
//	rawDataFullPathName = S_fileName		// save file name for future use
//
//	String line								// determine if this file is a new or old type
//	FReadLine refNum, line
//	FSetPos refNum,0
//	Close refNum
//	String noteStr=""
//	if (char2num(line[0])==char2num("$"))
//		noteStr = readInNewRawDataFile(rawDataFullPathName)	//	new
//	else
//		Abort "You have chosen an unknown type of file, it is not supported"
//	endif
//	if (strlen(noteStr)<1)
//		DoAlert 0, "Unable to properly read "+rawDataFullPathName
//		return ""
//	endif
//	return noteStr
//End

// read in data from file, and create arrays XX,HH,FF of the positions of each voxel, and (RX,RH,RF) the Rodriques vector for each voxel
// the positions are in µm, and the Rodriques vectors are in the XHF system, and  |Rodriques|=tan(angle/2).
Function/T readInRodriquesDataFile(FullFileName,maxAngle)
	String FullFileName
	Variable maxAngle
	if (!(maxAngle>0))
		maxAngle = !(maxAngle>0) ? 10 : maxAngle
		Prompt maxAngle, "reject all points with rotation angle greater than this (degree)"
		DoPrompt "max Rotation",maxAngle
		if (V_flag)
			return ""
		endif
	endif
	if (!(maxAngle>0))
		DoAlert 0, "ERROR in readInRodriquesDataFile(), maxAngle = "+num2str(maxAngle)
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
	Open /R/Z=1 refNum as FullFileName
	if (strlen(S_fileName)<1 || V_flag)
		return ""
	endif
	Close refNum

	String noteStr = microGeo#keyStrFromFile(FullFileName,"ArrayOf3dRotationsFile","")// read in all of the tagged values into a keyword list
	if (strlen(noteStr)<1)
		return ""
	endif
	Variable skipLines = DepthResolve#countLinesToTag(FullFileName,"XYZ_Rot")		// returns the number of lines in file before $XYZ_Rot

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

	// 1st three columns are X Y Z		X is PM500, so change X -> -X, Y and Z are in sample coordinates, so no negative signs
	// matij are in Wenge sample coordinates (X, -H, -F).  So negate the H and F coordinates when computing Rodriques vector
	String columnInfoStr = ""
	columnInfoStr += "N=XX;N=YY;N=ZZ;"
	columnInfoStr += "N=RX;N=RH;N=RF;"
	LoadWave/G/B=columnInfoStr/Q/A/L={0,skipLines,0,0,6} FullFileName	// start load after skipLines lines
	if (V_flag!=6)
		SetDataFolder fldrSav
		return ""
	endif
	Wave XX=XX,YY=YY,ZZ=ZZ					// PM500 positions (µm)
	Wave RX=RX, RH=RH, RF=RF
	Variable N=numpnts(XX)
	Make/N=(N) totalAngles
	SetScale d 0,0,"°", totalAngles
	Make/N=(3,3)/O/D refMat
	printf "read in %d voxels of raw data from   '%s'\r",N,FullFileName

	String list = StringByKey("refMat",noteStr,"=")
	refMat[0][0] = str2num(StringFromList(0,list,","))		// as of version 1.9 refMat is a recip-lattice not a rotation matrix
	refMat[0][1] = str2num(StringFromList(1,list,","))
	refMat[0][2] = str2num(StringFromList(2,list,","))
	refMat[1][0] = str2num(StringFromList(3,list,","))
	refMat[1][1] = str2num(StringFromList(4,list,","))
	refMat[1][2] = str2num(StringFromList(5,list,","))
	refMat[2][0] = str2num(StringFromList(6,list,","))
	refMat[2][1] = str2num(StringFromList(7,list,","))
	refMat[2][2] = str2num(StringFromList(8,list,","))

	Duplicate/O XX,HH,FF							// make XX,HH,FF in the XHF system (see note above)
	Variable i, rlen, degrees
	for (i=0;i<N;i+=1)
		HH[i] = YZ2H(YY[i],ZZ[i])				// H =  Y*sin(angle) + Z*cos(angle))	YY and ZZ are sample system in the file
		FF[i] = YZ2F(YY[i],ZZ[i])				// F = -Y*cos(angle) + Z*sin(angle)
		rlen = sqrt(RX[i]*RX[i] + RH[i]*RH[i] + RF[i]*RF[i])
		totalAngles[i] = 2*atan(rlen) * 180/PI// rotation angle in degrees
//		RX[i] = vec3[0]
//		RH[i] = -vec3[1]							// Wenge reversed definition of F and H for the matricies
//		RF[i] = -vec3[2]
	endfor
	KillWaves/Z YY,ZZ
	WaveStats/M=1/Q totalAngles
	printf "rotation angles range from:   %.3g° at %d  to  %.3g° at %d\r",V_min,V_minloc,V_max,V_maxloc


	if (maxAngle<180)							// may need to reject some points
		Variable rejected = 0
		for (i=N-1;i>=0;i-=1)
			if (totalAngles[i]>maxAngle)			// reject this point
				DeletePoints i, 1, XX,RX,RH,RF,totalAngles,FF,HH
				rejected += 1
			endif
		endfor
		printf "rejected %d points > %g (degree)\r",rejected, maxAngle
	endif

	Variable X0, Y0,Z0, H0, F0,  dX, dY, dZ, dH, dF
	X0 = NumberByKey("X0",noteStr,"=")
	Y0 = NumberByKey("Y0",noteStr,"=")
	Z0 = NumberByKey("Z0",noteStr,"=")
	if (numtype(X0+Y0+Z0))					// X0, Y0, or Z0 not in file, so use center of region as center
		WaveStats/Q XX
		X0 = V_avg
		WaveStats/Q HH
		H0 = V_avg
		WaveStats/Q FF
		F0 = V_avg
	else
		//	X0 = -X0								//	X0 in file is PM500 posiition
		H0 = YZ2H(Y0,Z0)						// Y0 and Z0 are already sample frame
		F0 = YZ2F(Y0,Z0)
	endif
	dX = NumberByKey("dX",noteStr,"=")
	dY = NumberByKey("dY",noteStr,"=")
	dZ = NumberByKey("dZ",noteStr,"=")
	//	dX = -dX									// in file dX is PM500 position, so reverse
	dH = YZ2H(dY,dZ)								// dY and dZ are in sample system in file
	dF = YZ2F(dY,dZ)
	Make/N=3/O/D vec3 = {dX,dH,dF}
	if (numtype(normalize(vec3)))
		dX=0  ;  dH=0  ;  dF=1					// default to the surface normal
	else
		dX = vec3[0]
		dH = vec3[1]
		dF = vec3[2]
	endif
	KillWaves/Z vec3
	noteStr = ReplaceNumberByKey("dX",noteStr,dX,"=")
	noteStr = ReplaceNumberByKey("dH",noteStr,dH,"=")
	noteStr = ReplaceNumberByKey("dF",noteStr,dF,"=")

	// center coordinates at the point {X0,H0,F0}.  This resets the origin to {0,0,0}
	XX -= X0
	HH -= H0
	FF -= F0
	printf "re-set origin of data to XHF={%g, %g, %g},   so now {0,0,0} is center\r",X0,H0,F0
	X0 = 0
	H0 = 0
	F0 = 0
	noteStr = ReplaceNumberByKey("X0",noteStr,X0,"=")	// save new X0, H0, F0
	noteStr = ReplaceNumberByKey("H0",noteStr,H0,"=")
	noteStr = ReplaceNumberByKey("F0",noteStr,F0,"=")
	noteStr = RemoveByKey("Y0", noteStr,"=")// remove the old Y0,Z0, dY, dZ
	noteStr = RemoveByKey("Z0", noteStr,"=")
	noteStr = RemoveByKey("dY", noteStr,"=")
	noteStr = RemoveByKey("dZ", noteStr,"=")
	noteStr = ReplaceStringByKey("waveClass",noteStr,"Random3dArrays","=")
	Note/K XX  ;  Note XX, noteStr
	Note/K HH  ;  Note HH, noteStr
	Note/K FF  ;  Note FF, noteStr
	Note/K RX  ;  Note RX, noteStr
	Note/K RH  ;  Note RH, noteStr
	Note/K RF  ;  Note RF, noteStr
	printf "remaining Igor data folder  '%s'\r",GetDataFolder(1)
//	SetDataFolder fldrSav
	return noteStr
End


// Read in reciprocal lattice from file, and create arrays XX,HH,FF of the positions of each voxel, and (RX,RH,RF) the 
// Rodriques vector for each voxel.  The positions are in µm, and the Rodriques vectors are in the XHF system, 
// and  |Rodriques|=tan(angle/2).
Function/T read3dRecipLatticesFile(FullFileName,maxAngle)
	String FullFileName
	Variable maxAngle
	if (!(maxAngle>0))
		maxAngle = !(maxAngle>0) ? 10 : maxAngle
		Prompt maxAngle, "reject all points with rotation angle greater than this (degree)"
		DoPrompt "max Rotation",maxAngle
		if (V_flag)
			return ""
		endif
	endif
	if (!(maxAngle>0))
		DoAlert 0, "ERROR in read3dRecipLatticesFile(), maxAngle = "+num2str(maxAngle)
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
	Open /R/Z=1 refNum as FullFileName
	if (strlen(S_fileName)<1 || V_flag)
		return ""
	endif
	Close refNum

	String noteStr = ReplaceStringByKey("localFile","",ParseFilePath(0,FullFileName,":",1,0),"=")	// add file name to note string
	noteStr = readSmarterHeader(FullFileName)
	if (strlen(noteStr)<1)
		return ""
	endif

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

	// turn info in XYZ_matrix into the columnInfoStr
	String list = StringByKey("XYZ_matrix",noteStr,"="), inWaves="", item
	list = ReplaceString("\t",list," ")
	do
		list = ReplaceString("  ",list," ")
	while(strsearch(list,"  ",0)>=0)
	list = ReplaceString(" ",list,";")					// list is now a semi-colon separated list
	Variable i, OKflag=0
	for (i=0;i<ItemsInLIst(list);i+=1)
		item = StringFromList(i,list)
		if (stringmatch(item,"X"))
			item = "XX"
			OKflag = OKflag | 1
		elseif (stringmatch(item,"Y"))
			item = "YY"
			OKflag = OKflag | 2
		elseif (stringmatch(item,"Z"))
			item = "ZZ"
			OKflag = OKflag | 4
		elseif (stringmatch(item,"ReciprocalLattice") || stringmatch(item,"RotationMatrix"))
			item = "mat11;mat12;mat13;mat21;mat22;mat23;mat31;mat32;mat33"
			OKflag = OKflag | 8
		else
			item = CleanupName(item,0)
		endif
		inWaves += item+";"							// list of waves read in
	endfor
	if (OKflag < 15)
		printf "read3dRecipLatticesFile,  Invalid columns\r%s\r",list
		return ""
	endif
	String columnInfoStr = ""
	for (i=0;i<ItemsInLIst(inWaves);i+=1)
		columnInfoStr += "N="+StringFromList(i,inWaves)+";"
	endfor
	// 1st three columns are X Y Z		X is PM500, so change X -> -X, Y and Z are in sample coordinates, so no negative signs
	// matij are in Wenge sample coordinates (X, -H, -F).  So negate the H and F coordinates when computing Rodriques vector
//	String columnInfoStr = ""
//	columnInfoStr += "N=XX;N=YY;N=ZZ;"
//	columnInfoStr += "N=mat11;N=mat12;N=mat13;"
//	columnInfoStr += "N=mat21;N=mat22;N=mat23;"
//	columnInfoStr += "N=mat31;N=mat32;N=mat33;"
	Variable iwaves = ItemsInList(columnInfoStr)
	LoadWave/G/B=columnInfoStr/Q/A/L={0,0,0,0,iwaves} FullFileName
	if (V_flag!=iwaves)
		SetDataFolder fldrSav
		return ""
	endif
	Wave XX=XX,YY=YY,ZZ=ZZ							// PM500 positions (µm)
	Wave mat11=mat11,mat12=mat12,mat13=mat13
	Wave mat21=mat21,mat22=mat22,mat23=mat23
	Wave mat31=mat31,mat32=mat32,mat33=mat33
	Wave irefMat=irefMat
	Variable N=numpnts(XX)
	Make/N=(N) RX,RH,RF,totalAngles
	SetScale d 0,0,"°", totalAngles
	Make/N=3/O/D vec3
	Make/N=(3,3)/O/D mat3, refMatTemp
	printf "read in %d voxels of raw data from   '%s'\r",N,FullFileName

	Variable NrefMat=NumberByKey("NrefMat",noteStr,"=")	// number of reference matrices present
	NrefMat = !(NrefMat>0) ? 1 : NrefMat						// defaults to 1
	Make/N=(3,3,NrefMat)/O/D refMat
	printf "setting %d reference matricies\r",NrefMat
	for (i=0;i<NrefMat;i+=1)									// set all of the reference matricies
		list = StringByKey("refMat"+SelectString(i,"",num2istr(i)),noteStr,"=")	// refMat, refMat1, refMat2, ...
		SetMatrixFromString(list,mat3)						// set mat from the values in list
		refMat[][][i] = mat3[p][q]
		MatrixOp/O mat3 = Inv(mat3) x mat3
		if (!(abs(3-MatrixTrace(mat3))<1e-9))				// check if refMat is invertible
			DoAlert 0, "refMat is not invertible, to it does not represent a valid 3D lattice.  This data is bad!"
			SetDataFolder fldrSav
			return ""
		endif
	endfor
	if (NrefMat>1 && !WaveExists(irefMat))			// make irefMat if it does not exist
		Abort "found multiple ref matricies, but no irefMat wave"
	endif

	Duplicate/O XX,HH,FF								// make XX,HH,FF in the XHF system (see note above)
	Make/N=(N)/I/U/O indicies						// contains original line numbers for backtracking, used by findClosestPoint()
	indicies = p
	Wave rho = $MakeUnique3x3Mat($"")				// make a new 3x3 matrix,for the rotations
	Variable j, angle, iref
	for (i=0,j=0;i<N;i+=1)
		mat3[0][0] = mat11[i]							// input matrix
		mat3[0][1] = mat12[i]
		mat3[0][2] = mat13[i]
		mat3[1][0] = mat21[i]
		mat3[1][1] = mat22[i]
		mat3[1][2] = mat23[i]
		mat3[2][0] = mat31[i]
		mat3[2][1] = mat32[i]
		mat3[2][2] = mat33[i]
		if (numtype(sum(mat3)))
			angle = NaN
			vec3 = NaN
		else
			iref = WaveExists(irefMat) ? irefMat[i] : 0
			refMatTemp = refMat[p][q][iref]			// for multiple refMat's, choose the correct refMat from  irefMat[]
			MatrixOp/O rho = mat3 x Inv(refMatTemp)	// rho x g0 = gm, the rotation matrix from ref to measured
			angle = axisOfMatrix(rho,vec3,squareUp=1)		// returned angle (degrees)
			if (numtype(angle))						// unable to find rotation axis
				printf "skip bad point (%d) found at XYZ= (%g, %g, %g)\r",i,XX[i],YY[i],ZZ[i]
				continue
			endif
			vec3 *= tan(angle*PI/180/2)				// this is now the Rodriques vector
		endif
		XX[j] = XX[i]
		HH[j] = YZ2H(YY[i],ZZ[i])						// H =  Y*sin(angle) + Z*cos(angle))	YY and ZZ are sample system in the file
		FF[j] = YZ2F(YY[i],ZZ[i])						// F = -Y*cos(angle) + Z*sin(angle)
		RX[j] = vec3[0]
		RH[j] = vec3[1]
		RF[j] = vec3[2]
		totalAngles[i] = angle
		j += 1											// only increments j for valid points
	endfor
	KillWaves/Z rho, YY,ZZ, mat11,mat12,mat13,mat21,mat22,mat23,mat31,mat32,mat33, mat3, refMatTemp
	N = j												// reset number of points
	inWaves = RemoveFromList("mat11;mat12;mat13;mat21;mat22;mat23;mat31;mat32;mat33",inWaves)	// change list so it is still waves of interest
	inWaves = RemoveFromList("YY;ZZ",inWaves)+"HH;FF;RX;RH;RF;totalAngles;indicies;"

	for (i=0;i<ItemsInList(inWaves);i+=1)
		Wave ww = $(StringFromList(i,inWaves))
		Redimension/N=(N)ww							// trim off end because I skipped invalid points
	endfor
	WaveStats/M=1/Q totalAngles
	printf "rotation angles range from:   %.3g° at %d  to  %.3g° at %d\r",V_min,V_minloc,V_max,V_maxloc

	if (maxAngle<180)									// may need to reject some points
		Variable rejected = 0
		for (i=N-1;i>=0;i-=1)
			if (totalAngles[i]>maxAngle)				// reject this point
				for (j=0;j<ItemsInList(inWaves);j+=1)
					Wave ww = $StringFromList(j,inWaves)
					DeletePoints i, 1, ww
				endfor
				rejected += 1
			endif
		endfor
		printf "rejected %d points > %g (degree)\r",rejected, maxAngle
	endif

	Variable X0, Y0,Z0, H0, F0,  dX, dY, dZ, dH, dF
	X0 = NumberByKey("X0",noteStr,"=")
	Y0 = NumberByKey("Y0",noteStr,"=")
	Z0 = NumberByKey("Z0",noteStr,"=")
	if (numtype(X0+Y0+Z0))							// X0, Y0, or Z0 not in file, so use center of region as center
		WaveStats/Q XX
		X0 = V_avg
		WaveStats/Q HH
		H0 = V_avg
		WaveStats/Q FF
		F0 = V_avg
	else
		//	X0 = -X0									//	X0 in file is PM500 posiition
		H0 = YZ2H(Y0,Z0)								// Y0 and Z0 are already sample frame
		F0 = YZ2F(Y0,Z0)
	endif
	dX = NumberByKey("dX",noteStr,"=")
	dY = NumberByKey("dY",noteStr,"=")
	dZ = NumberByKey("dZ",noteStr,"=")
	//	dX = -dX										// in file dX is PM500 position, so reverse
	dH = YZ2H(dY,dZ)									// dY and dZ are in sample system in file
	dF = YZ2F(dY,dZ)
	vec3 = {dX,dH,dF}
	if (numtype(normalize(vec3)))
		dX=0  ;  dH=0  ;  dF=1							// default to the surface normal
	else
		dX = vec3[0]
		dH = vec3[1]
		dF = vec3[2]
	endif
	noteStr = ReplaceNumberByKey("dX",noteStr,dX,"=")
	noteStr = ReplaceNumberByKey("dH",noteStr,dH,"=")
	noteStr = ReplaceNumberByKey("dF",noteStr,dF,"=")

	// center coordinates at the point {X0,H0,F0}.  This resets the origin to {0,0,0}
	XX -= X0
	HH -= H0
	FF -= F0
	printf "re-set origin of data to XHF={%g, %g, %g},   so now {0,0,0} is center\r",X0,H0,F0
	X0 = 0
	H0 = 0
	F0 = 0
	noteStr = ReplaceNumberByKey("X0",noteStr,X0,"=")	// save new X0, H0, F0
	noteStr = ReplaceNumberByKey("H0",noteStr,H0,"=")
	noteStr = ReplaceNumberByKey("F0",noteStr,F0,"=")
	noteStr = RemoveByKey("Y0", noteStr,"=")			// remove the old Y0,Z0, dY, dZ
	noteStr = RemoveByKey("Z0", noteStr,"=")
	noteStr = RemoveByKey("dY", noteStr,"=")
	noteStr = RemoveByKey("dZ", noteStr,"=")
	noteStr = ReplaceStringByKey("waveClass",noteStr,"Random3dArrays","=")
	for (i=0;i<ItemsInList(inWaves);i+=1)
		Wave ww = $StringFromList(i,inWaves)
		Note/K ww, noteStr
	endfor
	KillWaves/Z vec3
	printf "remaining Igor data folder  '%s'\r",GetDataFolder(1)
//	SetDataFolder fldrSav
	return noteStr
End
//
Static Function SetMatrixFromString(list,mat)			// set mat from the values in list, works for vectors too
	String list											// comma, semi-colon, or space separated list of elements
	Wave mat											// matrix already dimensioned to match the list
	if (DimSize(mat,0)<1 || strlen(list)<1)
		return 1
	endif

	list = ReplaceString(",",list,";")					// change all commas -> semi-colons
	list = ReplaceString("  ",list," ")					// change all double spaces to single
	list = ReplaceString(" ",list,";")					// change all spaces -> semi-colons
	Variable i,j, N=DimSize(mat,1)
	if (N<1)
		mat[] = str2num(StringFromList(p,list))
	else
		for (i=0;i<N;i+=1)
			mat[i][] = str2num(StringFromList(q+N*i,list))
		endfor
	endif
	return 0
	//	mat[0][0] = str2num(StringFromList(0,list,","))
	//	mat[0][1] = str2num(StringFromList(1,list,","))
	//	mat[0][2] = str2num(StringFromList(2,list,","))
	//	mat[1][0] = str2num(StringFromList(3,list,","))
	//	mat[1][1] = str2num(StringFromList(4,list,","))
	//	mat[1][2] = str2num(StringFromList(5,list,","))
	//	mat[2][0] = str2num(StringFromList(6,list,","))
	//	mat[2][1] = str2num(StringFromList(7,list,","))
	//	mat[2][2] = str2num(StringFromList(8,list,","))
End
//
Static Function DuplicateTheXs()
	Wave XX=XX, FF=FF, HH=HH
	Wave totalPeakIntensity=totalPeakIntensity
	Wave totalIntensity=totalIntensity
	Wave totalAngles=totalAngles
	Wave RF=RF, RH=RH, RX=RX

	WaveStats/Q/M=1 XX
	if ((V_max-V_min)>0.15)
		return 0					// nothing to do
	endif
	DoAlert 1,"This data set has only one X value, make it two X values?"
	if (V_flag!=1)
		return 1
	endif

	Variable N=numpnts(XX)
	Redimension/N=(2*N) XX,HH,FF,totalPeakIntensity,totalIntensity,totalAngles,RX,RH,RF
	Variable pLast = 2*N-1
	XX[N,pLast] = XX[p-N]+1
	HH[N,pLast] = HH[p-N]
	FF[N,pLast] = FF[p-N]
	RX[N,pLast] = RX[p-N]
	RH[N,pLast] = RH[p-N]
	RF[N,pLast] = RF[p-N]
	totalAngles[N,pLast] = totalAngles[p-N]
	totalIntensity[N,pLast] = totalIntensity[p-N]
	totalPeakIntensity[N,pLast] = totalPeakIntensity[p-N]
	return 0
End




//// This following routine is OLD and should not be used (only for data from Wenge)
//// read in data from file, and create arrays XX,HH,FF of the positions of each voxel, and (RX,RH,RF) the Rodriques vector for each voxel
//// the positions are in µm, and the Rodriques vectors are in the XHF system, and  |Rodriques|=tan(angle/2).
//Function/T readInNewRawDataFile(FullFileName,maxAngle)
//	String FullFileName
//	Variable maxAngle
//	if (!(maxAngle>0))
//		maxAngle = !(maxAngle>0) ? 10 : maxAngle
//		Prompt maxAngle, "reject all points with rotation angle greater than this (degree)"
//		DoPrompt "max Rotation",maxAngle
//		if (V_flag)
//			return ""
//		endif
//	endif
//	if (!(maxAngle>0))
//		DoAlert 0, "ERROR in readInNewRawDataFile(), maxAngle = "+num2str(maxAngle)
//		return ""
//	endif
//
//	String pathSep
//	if (stringmatch(IgorInfo(2),"Macintosh"))
//		pathSep = ":"
//	elseif (stringmatch(IgorInfo(2),"Windows"))
//		pathSep = "\\"
//	else
//		Abort "This is neither Mac nor Win, what is this computer?"
//	endif
//	Variable refNum
//	Open /R/Z=1 refNum as FullFileName
//	if (strlen(S_fileName)<1 || V_flag)
//		return ""
//	endif
//	Close refNum
//
//	String noteStr = ReplaceStringByKey("localFile","",ParseFilePath(0,FullFileName,":",1,0),"=")	// add file name to note string
//	noteStr = readSmarterHeader(FullFileName)
//	if (strlen(noteStr)<1)
//		return ""
//	endif
//
//	String fldrSav= GetDataFolder(1)
//	String fldrName = ParseFilePath(3,FullFileName, pathSep, 0, 0)
//	fldrName = CleanupName(fldrName, 0)
//	if (stringmatch(fldrName,"raw"))
//		fldrName = UniqueName(fldrName,11,0)
//	endif
//	if (CheckName(fldrName,11))
//		DoAlert 2, "This data folder already exists, delete the folder and re-load the data?"
//		if (V_flag==1)
//			KillDataFolder $fldrName
//		elseif (V_flag==2)
//			fldrName = UniqueName(fldrName,11,0)
//		else
//			return ""
//		endif
//	endif
//	NewDataFolder/O/S $fldrName
//	noteStr = ReplaceStringByKey("fldrName",noteStr,GetDataFolder(1),"=")
//
//	// 1st three columns are X Y Z		X is PM500, so change X -> -X, Y and Z are in sample coordinates, so no negative signs
//	// matij are in Wenge sample coordinates (X, -H, -F).  So negate the H and F coordinates when computing Rodriques vector
//	String columnInfoStr = ""
//	columnInfoStr += "N=XX;N=YY;N=ZZ;"
//	columnInfoStr += "N=mat11;N=mat12;N=mat13;"
//	columnInfoStr += "N=mat21;N=mat22;N=mat23;"
//	columnInfoStr += "N=mat31;N=mat32;N=mat33;"
//	LoadWave/G/B=columnInfoStr/Q/A/L={0,0,0,0,12} FullFileName
//	if (V_flag!=12)
//		SetDataFolder fldrSav
//		return ""
//	endif
//	Wave XX=XX,YY=YY,ZZ=ZZ					// PM500 positions (µm)
//	Wave mat11=mat11,mat12=mat12,mat13=mat13
//	Wave mat21=mat21,mat22=mat22,mat23=mat23
//	Wave mat31=mat31,mat32=mat32,mat33=mat33
//	Variable N=numpnts(XX)
//	Make/N=(N) RX,RH,RF,totalAngles
//	Make/N=3/O/D vec3
//	Make/N=(3,3)/O/D mat3, refMat
//	printf "read in %d voxels of raw data from   '%s'\r",N,FullFileName
//
//	String list = StringByKey("refMat",noteStr,"=")
//	refMat[0][0] = str2num(StringFromList(0,list,","))		// as of version 1.9 refMat is a recip-lattice not a rotation matrix
//	refMat[0][1] = str2num(StringFromList(1,list,","))
//	refMat[0][2] = str2num(StringFromList(2,list,","))
//	refMat[1][0] = str2num(StringFromList(3,list,","))
//	refMat[1][1] = str2num(StringFromList(4,list,","))
//	refMat[1][2] = str2num(StringFromList(5,list,","))
//	refMat[2][0] = str2num(StringFromList(6,list,","))
//	refMat[2][1] = str2num(StringFromList(7,list,","))
//	refMat[2][2] = str2num(StringFromList(8,list,","))
//	MatrixOp/O mat3 = Inv(refMat) x refMat
//	if (abs(3-MatrixTrace(mat3))>1e-9)		// check if refMat is invertible
//		DoAlert 0, "refMat is not invertible, to it does not represent a valid 3D lattice.  This data is bad!"
//		SetDataFolder fldrSav
//		return ""
//	endif
//
//	Duplicate/O XX,HH,FF							// make XX,HH,FF in the XHF system (see note above)
//	Variable i,j, angle
//	for (i=0,j=0;i<N;i+=1)
//		mat3[0][0] = mat11[i]					// input matrix
//		mat3[0][1] = mat12[i]
//		mat3[0][2] = mat13[i]
//		mat3[1][0] = mat21[i]
//		mat3[1][1] = mat22[i]
//		mat3[1][2] = mat23[i]
//		mat3[2][0] = mat31[i]
//		mat3[2][1] = mat32[i]
//		mat3[2][2] = mat33[i]
//		MatrixOp/O rho__ = Inv(mat3) x refMat// rho = gm^-1 x g0, rho should be a rotation matrix
//		if (notRotationMatrix(rho__))			// not a rotation matrix, skip this point
//			printf "skip bad point (%d) found at XYZ= (%g, %g, %g)\r",i,XX[i],YY[i],ZZ[i]
//			continue
//		endif
//		XX[j] = XX[i]
//		HH[j] = YZ2H(YY[i],ZZ[i])				// H =  Y*sin(angle) + Z*cos(angle))	YY and ZZ are sample system in the file
//		FF[j] = YZ2F(YY[i],ZZ[i])				// F = -Y*cos(angle) + Z*sin(angle)
//		angle = axisOfMatrix(rho__,vec3)		// returned angle (degrees)
//		vec3 *= tan(angle*PI/180/2)			// this is now the Rodriques vector
//		RX[j] = vec3[0]
//		RH[j] = -vec3[1]							// Wenge reversed definition of F and H for the matricies
//		RF[j] = -vec3[2]
//		totalAngles[i] = angle
//		j += 1										// only increments j for valid points
//	endfor
//	KillWaves/Z rho__, YY,ZZ
//	N = j											// reset number of points
//	Redimension/N=(N) XX,HH,FF,RX,RH,RF,totalAngles	// trim off end because I skipped invalid points
//	WaveStats/M=1/Q totalAngles
//	printf "rotation angles range from:   %.3g° at %d  to  %.3g° at %d\r",V_min,V_minloc,V_max,V_maxloc
//
//	if (maxAngle<180)							// may need to reject some points
//		Variable rejected = 0
//		for (i=N-1;i>=0;i-=1)
//			if (totalAngles[i]>maxAngle)			// reject this point
//				DeletePoints i, 1, XX,RX,RH,RF,totalAngles,FF,HH
//				rejected += 1
//			endif
//		endfor
//		printf "rejected %d points > %g (degree)\r",rejected, maxAngle
//	endif
//
//	Variable X0, Y0,Z0, H0, F0,  dX, dY, dZ, dH, dF
//	X0 = NumberByKey("X0",noteStr,"=")
//	Y0 = NumberByKey("Y0",noteStr,"=")
//	Z0 = NumberByKey("Z0",noteStr,"=")
//	if (numtype(X0+Y0+Z0))					// X0, Y0, or Z0 not in file, so use center of region as center
//		WaveStats/Q XX
//		X0 = V_avg
//		WaveStats/Q HH
//		H0 = V_avg
//		WaveStats/Q FF
//		F0 = V_avg
//	else
//		//	X0 = -X0								//	X0 in file is PM500 posiition
//		H0 = YZ2H(Y0,Z0)						// Y0 and Z0 are already sample frame
//		F0 = YZ2F(Y0,Z0)
//	endif
//	dX = NumberByKey("dX",noteStr,"=")
//	dY = NumberByKey("dY",noteStr,"=")
//	dZ = NumberByKey("dZ",noteStr,"=")
//	//	dX = -dX									// in file dX is PM500 position, so reverse
//	dH = YZ2H(dY,dZ)								// dY and dZ are in sample system in file
//	dF = YZ2F(dY,dZ)
//	vec3 = {dX,dH,dF}
//	if (numtype(normalize(vec3)))
//		dX=0  ;  dH=0  ;  dF=1					// default to the surface normal
//	else
//		dX = vec3[0]
//		dH = vec3[1]
//		dF = vec3[2]
//	endif
//	noteStr = ReplaceNumberByKey("dX",noteStr,dX,"=")
//	noteStr = ReplaceNumberByKey("dH",noteStr,dH,"=")
//	noteStr = ReplaceNumberByKey("dF",noteStr,dF,"=")
//
//	// center coordinates at the point {X0,H0,F0}.  This resets the origin to {0,0,0}
//	XX -= X0
//	HH -= H0
//	FF -= F0
//	printf "re-set origin of data to XHF={%g, %g, %g},   so now {0,0,0} is center\r",X0,H0,F0
//	X0 = 0
//	H0 = 0
//	F0 = 0
//	noteStr = ReplaceNumberByKey("X0",noteStr,X0,"=")	// save new X0, H0, F0
//	noteStr = ReplaceNumberByKey("H0",noteStr,H0,"=")
//	noteStr = ReplaceNumberByKey("F0",noteStr,F0,"=")
//	noteStr = RemoveByKey("Y0", noteStr,"=")// remove the old Y0,Z0, dY, dZ
//	noteStr = RemoveByKey("Z0", noteStr,"=")
//	noteStr = RemoveByKey("dY", noteStr,"=")
//	noteStr = RemoveByKey("dZ", noteStr,"=")
//	Note/K XX  ;  Note XX, noteStr
//	Note/K HH  ;  Note HH, noteStr
//	Note/K FF  ;  Note FF, noteStr
//	Note/K RX  ;  Note RX, noteStr
//	Note/K RH  ;  Note RH, noteStr
//	Note/K RF  ;  Note RF, noteStr
//
//	KillWaves/Z mat11,mat12,mat13,mat21,mat22,mat23,mat31,mat32,mat33
//	KillWaves/Z mat3, vec3
//	printf "remaining Igor data folder  '%s'\r",GetDataFolder(1)
////	SetDataFolder fldrSav
//	return noteStr
//End
////// read in data from file, and create arrays XX,HH,FF of the positions of each voxel, and (RX,RH,RF) the Rodriques vector for each voxel
////// the positions are in µm, and the Rodriques vectors are in the XHF system, and  |Rodriques|=tan(angle/2).
////Function/T readInNewRawDataFile(FullFileName)
////	String FullFileName
////
////	String pathSep
////	if (stringmatch(IgorInfo(2),"Macintosh"))
////		pathSep = ":"
////	elseif (stringmatch(IgorInfo(2),"Windows"))
////		pathSep = "\\"
////	else
////		Abort "This is neither Mac nor Win, what is this computer?"
////	endif
////	Variable refNum
////	Open /R/Z=1 refNum as FullFileName
////	if (strlen(S_fileName)<1 || V_flag)
////		return ""
////	endif
////	Close refNum
////
////	String noteStr = ReplaceStringByKey("localFile","",ParseFilePath(0,FullFileName,":",1,0),"=")	// add file name to note string
////	noteStr = readSmarterHeader(FullFileName)
////	if (strlen(noteStr)<1)
////		return ""
////	endif
////
////	String fldrSav= GetDataFolder(1)
////	String fldrName = ParseFilePath(3,FullFileName, pathSep, 0, 0)
////	fldrName = CleanupName(fldrName, 0)
////	if (stringmatch(fldrName,"raw"))
////		fldrName = UniqueName(fldrName,11,0)
////	endif
////	if (CheckName(fldrName,11))
////		DoAlert 2, "This data folder already exists, delete the folder and re-load the data?"
////		if (V_flag==1)
////			KillDataFolder $fldrName
////		elseif (V_flag==2)
////			fldrName = UniqueName(fldrName,11,0)
////		else
////			// Close refNum
////			return ""
////		endif
////	endif
////	NewDataFolder/O/S $fldrName
////
////	// 1st three columns are X Y Z		X is PM500, so change X -> -X, Y and Z are in sample coordinates, so no negative signs
////	// matij are in Wenge sample coordinates (X, -H, -F).  So negate the H and F coordinates when computing Rodriques vector
////	String columnInfoStr = ""
////	columnInfoStr += "N=XX;N=YY;N=ZZ;"
////	columnInfoStr += "N=mat11;N=mat12;N=mat13;"
////	columnInfoStr += "N=mat21;N=mat22;N=mat23;"
////	columnInfoStr += "N=mat31;N=mat32;N=mat33;"
////	LoadWave/G/B=columnInfoStr/Q/A/L={0,0,0,0,12} FullFileName
////	if (V_flag!=12)
////		SetDataFolder fldrSav
////		return ""
////	endif
////	Wave XX=XX,YY=YY,ZZ=ZZ					// PM500 positions (µm)
////	Wave mat11=mat11,mat12=mat12,mat13=mat13
////	Wave mat21=mat21,mat22=mat22,mat23=mat23
////	Wave mat31=mat31,mat32=mat32,mat33=mat33
////	Variable N=numpnts(XX)
////	Make/N=(N) RX,RH,RF,totalAngles
////	Make/N=3/O/D vec3
////	Make/N=(3,3)/O/D mat3, refMat
////	printf "read in %d voxels of raw data from   '%s'\r",N,FullFileName
////
////	String list = StringByKey("refMat",noteStr,"=")
////	refMat[0][0] = str2num(StringFromList(0,list,","))		// as of version 1.9 refMat is a recip-lattice not a rotation matrix
////	refMat[0][1] = str2num(StringFromList(1,list,","))
////	refMat[0][2] = str2num(StringFromList(2,list,","))
////	refMat[1][0] = str2num(StringFromList(3,list,","))
////	refMat[1][1] = str2num(StringFromList(4,list,","))
////	refMat[1][2] = str2num(StringFromList(5,list,","))
////	refMat[2][0] = str2num(StringFromList(6,list,","))
////	refMat[2][1] = str2num(StringFromList(7,list,","))
////	refMat[2][2] = str2num(StringFromList(8,list,","))
////	MatrixOp/O mat3 = Inv(refMat) x refMat
////	if (abs(3-MatrixTrace(mat3))>1e-9)		// check if refMat is invertible
////		DoAlert 0, "refMat is not invertible, to it does not represent a valid 3D lattice.  This data is bad!"
////		SetDataFolder fldrSav
////		return ""
////	endif
////
////	Duplicate/O XX,HH,FF							// make XX,HH,FF in the XHF system (see note above)
////	Variable i,j, angle
////	for (i=0,j=0;i<N;i+=1)
////		mat3[0][0] = mat11[i]					// input matrix
////		mat3[0][1] = mat12[i]
////		mat3[0][2] = mat13[i]
////		mat3[1][0] = mat21[i]
////		mat3[1][1] = mat22[i]
////		mat3[1][2] = mat23[i]
////		mat3[2][0] = mat31[i]
////		mat3[2][1] = mat32[i]
////		mat3[2][2] = mat33[i]
////		MatrixOp/O rho__ = Inv(mat3) x refMat// rho = gm^-1 x g0, rho should be a rotation matrix
////		if (notRotationMatrix(rho__))			// not a rotation matrix, skip this point
////			printf "skip bad point (%d) found at XYZ= (%g, %g, %g)\r",i,XX[i],YY[i],ZZ[i]
////			continue
////		endif
////		XX[j] = XX[i]
////		HH[j] = YZ2H(YY[i],ZZ[i])				// H =  Y*sin(angle) + Z*cos(angle))	YY and ZZ are sample system in the file
////		FF[j] = YZ2F(YY[i],ZZ[i])				// F = -Y*cos(angle) + Z*sin(angle)
////		angle = axisOfMatrix(rho__,vec3)		// returned angle (degrees)
////		vec3 *= tan(angle*PI/180/2)			// this is now the Rodriques vector
////		RX[j] = vec3[0]
////		RH[j] = -vec3[1]							// Wenge reversed definition of F and H for the matricies
////		RF[j] = -vec3[2]
////		totalAngles[i] = angle
////		j += 1										// only increments j for valid points
////	endfor
////	KillWaves/Z rho__, YY,ZZ
////	N = j											// reset number of points
////	Redimension/N=(N) XX,HH,FF,RX,RH,RF,totalAngles	// trim off end because I skipped invalid points
////
////	Variable X0, Y0,Z0, H0, F0,  dX, dY, dZ, dH, dF
////	X0 = NumberByKey("X0",noteStr,"=")
////	Y0 = NumberByKey("Y0",noteStr,"=")
////	Z0 = NumberByKey("Z0",noteStr,"=")
////	if (numtype(X0+Y0+Z0))					// X0, Y0, or Z0 not in file, so use center of region as center
////		WaveStats/Q XX
////		X0 = V_avg
////		WaveStats/Q HH
////		H0 = V_avg
////		WaveStats/Q FF
////		F0 = V_avg
////	else
////		//	X0 = -X0								//	X0 in file is PM500 posiition
////		H0 = YZ2H(Y0,Z0)						// Y0 and Z0 are already sample frame
////		F0 = YZ2F(Y0,Z0)
////	endif
////	dX = NumberByKey("dX",noteStr,"=")
////	dY = NumberByKey("dY",noteStr,"=")
////	dZ = NumberByKey("dZ",noteStr,"=")
////	//	dX = -dX									// in file dX is PM500 position, so reverse
////	dH = YZ2H(dY,dZ)								// dY and dZ are in sample system in file
////	dF = YZ2F(dY,dZ)
////	vec3 = {dX,dH,dF}
////	if (numtype(normalize(vec3)))
////		dX=0  ;  dH=0  ;  dF=1					// default to the surface normal
////	else
////		dX = vec3[0]
////		dH = vec3[1]
////		dF = vec3[2]
////	endif
////	noteStr = ReplaceNumberByKey("dX",noteStr,dX,"=")
////	noteStr = ReplaceNumberByKey("dH",noteStr,dH,"=")
////	noteStr = ReplaceNumberByKey("dF",noteStr,dF,"=")
////
////	// center coordinates at the point {X0,H0,F0}.  This resets the origin to {0,0,0}
////	XX -= X0
////	HH -= H0
////	FF -= F0
////	printf "re-set origin of data to XHF={%g, %g, %g},   so now {0,0,0} is center\r",X0,H0,F0
////	X0 = 0
////	H0 = 0
////	F0 = 0
////	noteStr = ReplaceNumberByKey("X0",noteStr,X0,"=")	// save new X0, H0, F0
////	noteStr = ReplaceNumberByKey("H0",noteStr,H0,"=")
////	noteStr = ReplaceNumberByKey("F0",noteStr,F0,"=")
////	noteStr = RemoveByKey("Y0", noteStr,"=")// remove the old Y0,Z0, dY, dZ
////	noteStr = RemoveByKey("Z0", noteStr,"=")
////	noteStr = RemoveByKey("dY", noteStr,"=")
////	noteStr = RemoveByKey("dZ", noteStr,"=")
////	Note/K XX  ;  Note XX, noteStr
////	Note/K HH  ;  Note HH, noteStr
////	Note/K FF  ;  Note FF, noteStr
////	Note/K RX  ;  Note RX, noteStr
////	Note/K RH  ;  Note RH, noteStr
////	Note/K RF  ;  Note RF, noteStr
////
////	KillWaves/Z mat11,mat12,mat13,mat21,mat22,mat23,mat31,mat32,mat33
////	KillWaves/Z mat3, vec3
////	SetDataFolder fldrSav
////	return noteStr
////End


// readin the header part
Function/T readSmarterHeader(fullFileName)
	String fullFileName
	Variable refNum
	if (WhichListItem("ArrayOf3dOrientsFile",getListOfTypesInFile(fullFileName,"home"))<0)
		return ""												// this file is wrong type
	endif

	Open /R/M="List of Reciprocal Lattices"/P=home/Z=1 refNum as fullFileName
	if (strlen(S_fileName)<1 || V_flag)
		return ""
	endif
	FStatus refNum
	String buf = PadString("",min(V_logEOF,3*1024),0)
	FBinRead refNum, buf
	Close refNum
	String noteStr = microGeo#keyStrFromBuffer(buf)

	if (!keyInList("refMat",noteStr,"=",";"))					// set default refMat if none present
		noteStr = ReplaceStringByKey("refMat",noteStr,"1,0,0,0,1,0,0,0,1","=")// default ref matrix
	endif
	return noteStr
End
//Function/T readSmarterHeader(fullFileName)
//	String fullFileName
//	Variable refNum
//	Open /R/M="Wirescan peak positions from Wenge"/P=home/Z=1 refNum as fullFileName
//	if (strlen(S_fileName)<1 || V_flag)
//		return ""
//	endif
//	if (!stringmatch(findTagInFile(refNum,"filetype"),"ArrayOf3dOrientsFile"))
//		Close refNum
//		return ""												// this file is wrong type
//	endif
//
//	String noteStr = ReplaceStringByKey("refMat","","1,0,0,0,1,0,0,0,1","=")// default ref matrix
//	noteStr = fetchAndAddNumber2List(refNum,"X0",noteStr,1)		// X-coord of cylindrical symmetry axis (µm)
//	noteStr = fetchAndAddNumber2List(refNum,"Y0",noteStr,1)		// Y-coord of cylindrical symmetry axis (µm)
//	noteStr = fetchAndAddNumber2List(refNum,"Z0",noteStr,1)		// Z-coord of cylindrical symmetry axis (µm)
//	noteStr = fetchAndAddNumber2List(refNum,"dX",noteStr,1)		// X component of unit vector in direction of cylindrical axis, beam line coordinates
//	noteStr = fetchAndAddNumber2List(refNum,"dY",noteStr,1)		// Y component of unit vector in direction of cylindrical axis
//	noteStr = fetchAndAddNumber2List(refNum,"dZ",noteStr,1)		// Z component of unit vector in direction of cylindrical axis
//	noteStr = fetchAndAddSring2List(refNum,"refMat",noteStr)		// reference matrix identifies zero angle
//	noteStr = fetchAndAddSring2List(refNum,"date",noteStr)		// date file was created
//	Close refNum
//	return noteStr
//End


////	Dduplicate/O xx AA
////	AA = 4
////	test(AA)
////	test(XX)
//Function testFindScalingFromVec(ww)
//	Wave ww
//	String name = NameOfWave(ww)
//	Variable threshold = 0.18							// (µm) a motion of the PM500 greater than this is intentional (less is jitter)
//	Variable off,d, dim, err
//	err = FindScalingFromVec(ww,threshold,off,d,dim)
//	printf "err=%d,   %s=[%g,%g], dim%s=%d,   ∆%s=%g\r",err,name,off,off+d*(dim-1),name,dim,name,d
//End
//
// find the step size in a coordinate by examining the values
Function FindScalingFromVec(vec,threshold,first,stepSize,dimN)
	Wave vec
	Variable threshold			// a change greater than this is intentional (less is jitter)
	Variable &first				// use 	SetScale/P x first,stepSize,"",waveName
	Variable &stepSize
	Variable &dimN

	Duplicate/O/D vec testDim
	Variable N = numpnts(vec)-1
	Redimension/N=(N) testDim
	SetScale/P x 0,1,"", testDim
	testDim = vec[p+1] - vec[p]
	testDim = abs(testDim)<threshold ? NaN : testDim
	Sort testDim, testDim
	WaveStats/M=1/Q testDim
	N = V_npnts
	if (N<1)								// no stepss, all the same value
		WaveStats/M=1/Q vec
		first = V_avg
		stepSize = 0
		dimN=1
		KillWaves/Z testDim
		return 0
	endif
	Redimension/N=(N) testDim			// trim off trailing NaNs
	stepSize = abs(testDim[round(N/2)])
	stepSize = mod(round(stepSize*25)/25,1)==0 ? round(0.98*25)/25 : stepSize

	Redimension/N=(numpnts(vec)) testDim	// next get the first value
	testDim = vec
	Sort testDim, testDim
	WaveStats/M=1/Q testDim
	N = V_npnts
	Redimension/N=(N) testDim
	Variable i = floor(BinarySearchInterp(testDim, testDim[0]+threshold))
//	first = faverage(testDim,0,i)			// average value
	first = testDim[round(i/2)]			// median value
	first = (first==0) ? 0 : first			// I do not like "-0"
	if (ItemsInList(GetRTStackInfo(0))<=1)
		print stepSize, first
	endif
	i = ceil(BinarySearchInterp(testDim, testDim[N-1]-threshold))
	dimN = round((faverage(testDim,i,N-1)-first)/stepSize + 1)
	KillWaves/Z testDim
	return 0
End
//Function FindScalingFromVec(vec,threshold,first,stepSize,dimN)
//	Wave vec
//	Variable threshold			// a change greater than this is intentional (less is jitter)
//	Variable &first				// use 	SetScale/P x first,stepSize,"",waveName
//	Variable &stepSize
//	Variable &dimN
//
//	Variable last,i, Ndim=1
//	Duplicate/O vec testDim
//	Variable N = numpnts(testDim)
//	SetScale/P x 0,1,"", testDim
//	Sort testDim, testDim
//	WaveStats/M=1/Q testDim
//	Redimension/N=(V_npnts) testDim		// trim off trailing NaNs
//
//	Variable delta
//	stepSize=Inf
//	for (i=0;i<(N-1);i+=1)
//		delta = testDim[i+1]-testDim[i]
//		if (delta<stepSize && delta>threshold)
//			stepSize = delta
//		endif
//	endfor
//	if (numtype(stepSize)==1)				// no steps, all are the same
//		WaveStats/M=1/Q testDim
//		first = V_avg
//		stepSize = 0
//		dimN=1
//		return 0
//	endif
//	// stepSize is now approximately right
//
//	Variable Nj, stepSum=0
//	for (i=0;i<(N-1);i+=1)
//		delta = testDim[i+1]-testDim[i]
//		if (round(delta/stepSize)==1)
//			stepSum += delta
//			Nj += 1
//		endif
//	endfor
//	if (Nj>0)
//		stepSize = stepSum/Nj				// this is now a good average step size
//	endif
//
//	i = floor(BinarySearchInterp(testDim, testDim[0]+threshold))
//	first = faverage(testDim,0,i)
//	first = (first==0) ? 0 : first			// I do not like "-0"
//	if (ItemsInList(GetRTStackInfo(0))<=1)
//		print stepSize, first
//	endif
//
//	i = ceil(BinarySearchInterp(testDim, testDim[N-1]-threshold))
//	dimN = round((faverage(testDim,i,N-1)-first)/stepSize + 1)
//	KillWaves/Z testDim
//	return 0
//End



// retruns axis name assuming vec is in XHF coordinates
Function/S axisName(vec)
	Wave vec
	Wave axis=$MakeUnique3Vector(vec)	// normalized version of vec
	normalize(axis)
	Wave hat=$MakeUnique3Vector($"")
	Variable one=cos(0.1*PI/180)	// almost one (0.1 degree off)
	Variable dot

	hat = {1,0,0}
	dot = MatrixDot(axis,hat)
	if (abs(dot)>one)
		if (dot>0)
			KillWaves/Z axis,hat
			return "X"
		else
			KillWaves/Z axis,hat
			return "-X"
		endif
	endif

	hat = {0,1,0}
	dot = MatrixDot(axis,hat)
	if (abs(dot)>one)
		if (dot>0)
			KillWaves/Z axis,hat
			return "H"
		else
			KillWaves/Z axis,hat
			return "-H"
		endif
	endif

	hat = {0,0,1}
	dot = MatrixDot(axis,hat)
	if (abs(dot)>one)
		if (dot>0)
			KillWaves/Z axis,hat
			return "F (depth)"
		else
			KillWaves/Z axis,hat
			return "-F (depth)"
		endif
	endif

	hat = {0,1,-1}			// Y direction,  Y = (H-F)/sqrt(2)
	normalize(hat)
	dot = MatrixDot(axis,hat)
	if (abs(dot)>one)
		if (dot>0)
			KillWaves/Z axis,hat
			return "Y"
		else
			KillWaves/Z axis,hat
			return "-Y"
		endif
	endif

	hat = {0,1,1}			// Z direction,  Z = (H+F)/sqrt(2)
	normalize(hat)
	dot = MatrixDot(axis,hat)
	if (abs(dot)>one)
		if (dot>0)
			KillWaves/Z axis,hat
			return "Z"
		else
			KillWaves/Z axis,hat
			return "-Z"
		endif
	endif

	// not one of the standard directions, so just give direction
	String str
	sprintf str, "along (%.3g, %.3g, %.3g)",axis[0],axis[1],axis[2]
	KillWaves/Z axis,hat
	return str
End


// make a unit vector in the direction indicated by the string 'direction'
// returns the full name of the new unit vector, returns "" for invalid input
Function/S MakeDirectionVector(direction)		// make a unit vector in the specified direction
	String direction

	Wave normal = $MakeUnique3Vector($"")	// normal in the requested direction
	Variable H,F
	strswitch(direction)
		case "X":
			normal={1,0,0}
			break
		case "H":
			normal={0,1,0}
			break
		case "F":
			normal={0,0,1}
			break
		case "Y":
			H = YZ2H(1,0)
			F = YZ2F(1,0)
			normal={0,H,F}
			break
		case "Z":
			H = YZ2H(0,1)
			F = YZ2F(0,1)
			normal={0,H,F}
			break
		case "phiZ":
			Variable phi = NumVarOrDefault("phiSlice",0)*PI/180
			normal[0] = cos(phi)				// for a phi cut around F
			normal[1] = sin(phi)
			normal[2] = 0
			break
		case "R_Cyl":							// for cylindrical coordinates, need to recalc directions for each point
		case "Phi_Cyl":
		case "Axis_Cyl":
		case "kappa(xx)":						// for lattice curvature tensor components, need to calc for each point too.
		case "kappa(xy)":
		case "kappa(xz)":
		case "kappa(yx)":
		case "kappa(yy)":
		case "kappa(yz)":
		case "kappa(zx)":
		case "kappa(zy)":
		case "kappa(zz)":
		case "alpha(xx)":						// for dislocation tensor components, need to calc for each point too.
		case "alpha(xy)":
		case "alpha(xz)":
		case "alpha(yx)":
		case "alpha(yy)":
		case "alpha(yz)":
		case "alpha(zx)":
		case "alpha(zy)":
		case "alpha(zz)":
		case "GND":
		case "dRdX":
		case "dRdY":
		case "dRdZ":
		case "|gradR|":
		case "exx":
		case "eyy":
		case "ezz":
		case "exy":
		case "exz":
		case "eyz":
			normal = NaN
			break
		default:
			DoAlert 0, "Illegal direction vector ''"+direction+"'"
			KillWaves/Z normal
			return ""
	endswitch
	return GetWavesDataFolder(normal,2)
End


// moved to microGeometry.ipf
//
//// fix up a matrix so that it is exactly a rotation matrix, not just close to one
//Function SquareUpMatrix(mat)
//	Wave mat
//
//	Wave vec0=$MakeUnique3Vector($"")
//	Wave vec1=$MakeUnique3Vector($"")
//	Wave vec2=$MakeUnique3Vector($"")
//	Variable err
//	vec0 = mat[p][0]
//	vec1 = mat[p][1]
//	vec2 = mat[p][2]
//	Variable norm0=norm(vec0), norm1=norm(vec1), norm2=norm(vec2)
//
//	// Start with the longest vector, and assume it is correct
//	if (norm0>=norm1 && norm0>=norm2)	// X is longest
//		Cross vec0,vec1						// Z = X x Y
//		Wave W_Cross=W_Cross
//		vec2 = W_Cross
//		Cross vec2,vec0						// Y = Z x X
//		vec1 = W_Cross
//	elseif (norm1>=norm0 && norm1>=norm2)// Y is longest
//		Cross vec1,vec2						// X = Y x Z
//		Wave W_Cross=W_Cross
//		vec0 = W_Cross
//		Cross vec0,vec1						// Z = X x Y
//		vec2 = W_Cross
//	else											// Z is longest
//		Cross vec2,vec0						// Y = Z x X
//		Wave W_Cross=W_Cross
//		vec1 = W_Cross
//		Cross vec1,vec2						// X = Y x Z
//		vec0 = W_Cross
//	endif
//
//	err = normalize(vec0)
//	err += normalize(vec1)
//	err += normalize(vec2)
//	mat[][0] = vec0[p]
//	mat[][1] = vec1[p]
//	mat[][2] = vec2[p]
//
//	err = ( notRotationMatrix(mat) != 0 )		// 0 is a rotation matrix
//	KillWaves/Z vec0,vec1,vec2,W_Cross
//	return numtype(err)
//End
//
//// fix up a matrix so that it is exactly a rotation matrix, not just close to one
//OverRide Function SquareUpMatrix(mat)
//	Wave mat
//
//	Wave vec0=$MakeUnique3Vector($"")
//	Wave vec1=$MakeUnique3Vector($"")
//	Wave vec2=$MakeUnique3Vector($"")
//	Variable err
//
//	vec0 = mat[p][0]
//	vec1 = mat[p][1]
//
//	Cross vec0,vec1
//	Wave W_Cross=W_Cross
//	vec2 = W_Cross
//
//	Cross vec2,vec0
//	vec1 = W_Cross
//
//	err = normalize(vec0)
//	err += normalize(vec1)
//	err += normalize(vec2)
//	mat[][0] = vec0[p]
//	mat[][1] = vec1[p]
//	mat[][2] = vec2[p]
//
//	err = ( notRotationMatrix(mat) != 0 )	// 0 is a rotation matrix
//
//	KillWaves/Z vec0,vec1,vec2,W_Cross
//	return numtype(err)
//End
//
//	// fix up a matrix so that it is exactly a rotation matrix, not just close to one
//Function SquareUpMatrix(mat)
//	Wave mat
//
//	Wave vec=$MakeUnique3Vector($"")
//	Variable len
//	Variable err
//
//	vec[] = mat[0][p]
//	len = norm(vec)
//	mat[0][] = vec[q]/len
//	err = len
//
//	vec[] = mat[1][p]
//	len = norm(vec)
//	mat[1][] = vec[q]/len
//	err += len
//
//	vec[] = mat[2][p]
//	len = norm(vec)
//	mat[2][] = vec[q]/len
//	err += len
//
//	KillWaves/Z vec
//	return numtype(err)
//End




// =========================================================================
// =========================================================================
//	Start of plotting for raw data

Function rawCalcToGizmo(maxDist)
	Variable maxDist

	Wave XX=XX, HH=HH, FF=FF, RX=RX, RH=RH, RF=RF, totalAngles=totalAngles
	if (!WaveExists(XX) || !WaveExists(HH) || !WaveExists(FF) || !WaveExists(RX) || !WaveExists(RH) || !WaveExists(RF) || !WaveExists(totalAngles))
		DoAlert 0,"cannot find XX, HH, & etc.\rProbably in wrong DataFolder."
		return 1
	endif
	Variable N=numpnts(XX)
	if (!(maxDist>0))
		maxDist = N>200 ? 20 : Inf
		Prompt maxDist, "max dist to graph"
		DoPrompt "max distance",maxDist
		if (V_flag)
			return 1
		endif
	endif
	String wnote = ReplaceStringByKey("waveClass",note(XX),"Orient3dRawGizmo","=")

	Make/N=(N,3)/O OrientsXYZ=NaN, OrientsSize=NaN
	Make/N=(N,4)/O OrientsRGBA=NaN, OrientsRot=NaN
	Make/N=3/O/D tempVec_rawCalcToGizmo, zhat_rawCalcToGizmo={0,0,1}
	Wave vec=tempVec_rawCalcToGizmo, zhat=zhat_rawCalcToGizmo

	Variable r,g,b, angle, swap
	Variable i,m
	for (m=0,i=0;i<N;i+=1)
		if (abs(XX[i])>maxDist || abs(HH[i])>maxDist || abs(FF[i])>maxDist)
			continue
		endif

//if (abs(FF[i]) < 12)
//continue
//endif

		OrientsXYZ[m][0] = XX[i]	;	OrientsXYZ[m][1] = HH[i]	;	OrientsXYZ[m][2] = FF[i]
		r = RX[i]					;	g = RH[i]					;	b = RF[i]
		vec = {r,g,b}												// save here for rotation
		getColor(r,g,b)
		OrientsRGBA[m][0] = r		;	OrientsRGBA[m][1] = g		;	OrientsRGBA[m][2] = b
		OrientsSize[m][] = sqrt(totalAngles[i])

		normalize(vec)												// set direction of the arrows
		Cross/Z zhat,vec											// for no rotation, arrows point along +Z
		Wave W_Cross=W_Cross
		angle = asin(normalize(W_Cross))*180/PI
		if (angle<0)
			W_Cross *= -1
			angle = -angle
		endif
		OrientsRot[m][0] = angle
		OrientsRot[m][1,3] = W_Cross[q-1]
		m += 1
	endfor
	N = m
	KillWaves/Z tempVec_rawCalcToGizmo, zhat_rawCalcToGizmo, W_Cross

	N += numtype(maxDist) ? 0 : 8									// add extra points for corners if using maxDist
	Redimension/N=(N,3) OrientsXYZ, OrientsSize					// trim off un-used  (because of maxDist), and add points for corners
	Redimension/N=(N,4) OrientsRGBA, OrientsRot
	if (numtype(maxDist)==0)
		OrientsXYZ[N-8][0] = -maxDist ;		OrientsXYZ[N-8][1] = -maxDist ;		OrientsXYZ[N-8][2] = 0		// add corner points
		OrientsXYZ[N-7][0] = -maxDist ;		OrientsXYZ[N-7][1] = +maxDist ;		OrientsXYZ[N-7][2] = 0
		OrientsXYZ[N-6][0] = +maxDist ;		OrientsXYZ[N-6][1] = -maxDist ;		OrientsXYZ[N-6][2] = 0
		OrientsXYZ[N-5][0] = +maxDist ;		OrientsXYZ[N-5][1] = +maxDist ;		OrientsXYZ[N-5][2] = 0
		OrientsXYZ[N-4][0] = -maxDist ;		OrientsXYZ[N-4][1] = -maxDist ;		OrientsXYZ[N-4][2] = maxDist
		OrientsXYZ[N-3][0] = -maxDist ;		OrientsXYZ[N-3][1] = +maxDist ;		OrientsXYZ[N-3][2] = maxDist
		OrientsXYZ[N-2][0] = +maxDist ;		OrientsXYZ[N-2][1] = -maxDist ;		OrientsXYZ[N-2][2] = maxDist
		OrientsXYZ[N-1][0] = +maxDist ;		OrientsXYZ[N-1][1] = +maxDist ;		OrientsXYZ[N-1][2] = maxDist
		OrientsRGBA[N,N-1][] = 0									// make corners invisible
		OrientsSize[N,N-1][] = 0
	endif

	WaveStats/M=1/Q OrientsRGBA
	OrientsRGBA /= V_max
	//printf "V_max=%g,   or an angle of %g°\r",V_max,atan(V_max)*2*180/PI
//	OrientsRGBA[][3] = 0.4
	OrientsRGBA[][3] = 1-0.8*sqrt(OrientsRGBA[p][0]^2 + OrientsRGBA[p][1]^2 + OrientsRGBA[p][2]^2)
	OrientsRGBA = limit(OrientsRGBA[p][q],0,1)

	WaveStats/M=1/Q OrientsSize
	Variable maxAngle = V_max
	OrientsSize /= maxAngle*1.5
	OrientsSize = limit(OrientsSize[p][q],0,1)
	Note/K OrientsXYZ, wnote
	Note/K OrientsSize, wnote
	Note/K OrientsRGBA, wnote
	Note/K OrientsRot, wnote
	printf "for a maximum coordinate of %+g,  the maximum rotation angle = %g°\r",maxDist,maxAngle
	return 0
End
//
Static Function getColor(r,g,b)
	Variable &r, &g, &b									// x, y, z on input, returns r,g,b

	Variable x=r, y=g, z=b								//initial values
	Variable len = sqrt(x*x + y*y + z*z)
	if (len==0)
		return  0
	endif

	r=0; g=0; b=0
	if (x>0)
		r += x
	else
		r += 0 ;	g += -x ;	b += -x
	endif
	if (y>0)
		g += y
	else
		r += -y ;	g += 0 ;	b += -y
	endif
	if (z>0)
		b += z
	else
		r += -z ;	g += -z ;	b += 0
	endif

	Variable ratio = len/sqrt(r*r + g*g + b*b)		// re-scale so length is the same
	r *= ratio
	g *= ratio
	b *= ratio
	return len
End


Function MakeOrientationsGizmo(arrows)					// : GizmoPlot
	Variable arrows
	if(exists("NewGizmo")!=4)						// Do nothing if the Gizmo XOP is not available.
		DoAlert 0, "Gizmo XOP must be installed"
		return 1
	endif
	String fldr = GetDataFolder(1)
	if (!(arrows>=0))
		arrows = !(arrows>=0) ? 2 : arrows+1
		Prompt arrows,"maker type",popup,"balls;arrows"
		DoPrompt "markers",arrows
		if (V_flag)
			return 1
		endif
		arrows -= 1								// change from {1,2} to {0,1}
	endif

	Execute "NewGizmo/N="+UniqueName("Gizmo",5, 0)+"/T=\""+GetDataFolder(0)+"\" /W=(147,80,1064,926)"
	Execute "ModifyGizmo startRecMacro"

	if (arrows)
		// ************************* Group Object Start *******************
		Execute "AppendToGizmo group,name=group0"
		Execute "ModifyGizmo currentGroupObject=\"group0\""
		Execute "AppendToGizmo cylinder={0.01,0.01,0.75,25,25},name=cylinder0"
		Execute "	ModifyGizmo setObjectAttribute={cylinder0,black}"
		Execute "AppendToGizmo cylinder={0.1,0,0.3,25,25},name=cylinder1"
		Execute "AppendToGizmo attribute color={0.8,0,0,0.3},name=red"
		Execute "AppendToGizmo attribute color={0,0,1,0.3},name=blue"
		Execute "AppendToGizmo attribute color={0,0,0,0.2},name=black"
		Execute "ModifyGizmo setDisplayList=0, opName=enable0, operation=enable, data=2903"
		Execute "ModifyGizmo setDisplayList=1, opName=colorMaterial0, operation=colorMaterial, data={1032,5634}"
		Execute "ModifyGizmo setDisplayList=2, opName=translate1, operation=translate, data={0,0,-0.5}"
		Execute "ModifyGizmo setDisplayList=3, object=cylinder0"
		Execute "ModifyGizmo setDisplayList=4, opName=translate0, operation=translate, data={0,0,0.75}"
		Execute "ModifyGizmo setDisplayList=5, object=cylinder1"
		Execute "ModifyGizmo setDisplayList=6, opName=disable0, operation=disable, data=2903"
		Execute "ModifyGizmo currentGroupObject=\"::\""
		// ************************* Group Object End *******************
	endif

	Execute "AppendToGizmo Scatter="+fldr+"OrientsXYZ,name=scatter0"
	Execute "ModifyGizmo ModifyObject=scatter0 property={ scatterColorType,1}"
	Execute "ModifyGizmo ModifyObject=scatter0 property={ markerType,0}"
	Execute "ModifyGizmo ModifyObject=scatter0 property={ sizeType,1}"
	Execute "ModifyGizmo ModifyObject=scatter0 property={ size,0.2}"
	if (arrows)
		Execute "ModifyGizmo ModifyObject=scatter0 property={ rotationType,1}"
		Execute "ModifyGizmo ModifyObject=scatter0 property={ Shape,7}"
		Execute "ModifyGizmo ModifyObject=scatter0 property={ rotationWave,"+fldr+"OrientsRot}"
		Execute "ModifyGizmo ModifyObject=scatter0 property={ objectName,group0}"
	else
		Execute "ModifyGizmo ModifyObject=scatter0 property={ rotationType,0}"
		Execute "ModifyGizmo ModifyObject=scatter0 property={ Shape,2}"
	endif
	Execute "ModifyGizmo ModifyObject=scatter0 property={ colorWave,"+fldr+"OrientsRGBA}"
	Execute "ModifyGizmo ModifyObject=scatter0 property={ sizeWave,"+fldr+"OrientsSize}"

	// the axes
	Execute "AppendToGizmo Axes=boxAxes,name=axes0"
	Execute "ModifyGizmo ModifyObject=axes0,property={-1,axisScalingMode,1}"
	Execute "ModifyGizmo ModifyObject=axes0,property={-1,axisColor,0,0,0,1}"
	Execute "ModifyGizmo ModifyObject=axes0,property={0,ticks,3}"
	Execute "ModifyGizmo ModifyObject=axes0,property={1,ticks,3}"
	Execute "ModifyGizmo ModifyObject=axes0,property={2,ticks,3}"
	Execute "ModifyGizmo ModifyObject=axes0,property={4,ticks,3}"
	Execute "ModifyGizmo ModifyObject=axes0,property={0,axisLabel,1}"
	Execute "ModifyGizmo ModifyObject=axes0,property={1,axisLabel,1}"
	Execute "ModifyGizmo ModifyObject=axes0,property={2,axisLabel,1}"
	Execute "ModifyGizmo ModifyObject=axes0,property={4,axisLabel,1}"
	Execute "ModifyGizmo ModifyObject=axes0,property={0,axisLabelText,\"X  (µm)\"}"
	Execute "ModifyGizmo ModifyObject=axes0,property={1,axisLabelText,\"H  (µm)\"}"
	Execute "ModifyGizmo ModifyObject=axes0,property={2,axisLabelText,\"F  (µm)\"}"
	Execute "ModifyGizmo ModifyObject=axes0,property={4,axisLabelText,\"F  (µm)\"}"
	Execute "ModifyGizmo ModifyObject=axes0,property={0,axisLabelCenter,0}"
	Execute "ModifyGizmo ModifyObject=axes0,property={1,axisLabelCenter,0}"
	Execute "ModifyGizmo ModifyObject=axes0,property={2,axisLabelCenter,0}"
	Execute "ModifyGizmo ModifyObject=axes0,property={4,axisLabelCenter,0}"
	Execute "ModifyGizmo ModifyObject=axes0,property={0,axisLabelDistance,0.05}"
	Execute "ModifyGizmo ModifyObject=axes0,property={1,axisLabelDistance,0.05}"
	Execute "ModifyGizmo ModifyObject=axes0,property={2,axisLabelDistance,0.2}"
	Execute "ModifyGizmo ModifyObject=axes0,property={4,axisLabelDistance,0.2}"
	Execute "ModifyGizmo ModifyObject=axes0,property={0,axisLabelScale,0.5}"
	Execute "ModifyGizmo ModifyObject=axes0,property={1,axisLabelScale,0.5}"
	Execute "ModifyGizmo ModifyObject=axes0,property={2,axisLabelScale,0.5}"
	Execute "ModifyGizmo ModifyObject=axes0,property={4,axisLabelScale,0.5}"
	Execute "ModifyGizmo ModifyObject=axes0,property={0,axisLabelRGBA,0,0,0,1}"
	Execute "ModifyGizmo ModifyObject=axes0,property={1,axisLabelRGBA,0,0,0,1}"
	Execute "ModifyGizmo ModifyObject=axes0,property={2,axisLabelRGBA,0,0,0,1}"
	Execute "ModifyGizmo ModifyObject=axes0,property={4,axisLabelRGBA,0,0,0,1}"

	// a line going along the axis, along (0,0,z)
	Execute "AppendToGizmo line={0,0,-1,0,0,1}, name=axisLine"
	Execute "AppendToGizmo attribute lineWidth=2, name=lineWidth0"
	Execute "AppendToGizmo attribute color={0,0.5,1,1},name=axisColor"
	Execute "AppendToGizmo attribute blendFunc={770,771},name=blendFunc0"

	Execute "ModifyGizmo setDisplayList=0, attribute=blendFunc0"
	Execute "ModifyGizmo setDisplayList=1, object=scatter0"
	Execute "ModifyGizmo setDisplayList=2, object=axes0"
	Execute "ModifyGizmo setDisplayList=3, attribute=lineWidth0"
	Execute "ModifyGizmo setDisplayList=4, attribute=axisColor"
	Execute "ModifyGizmo setDisplayList=5, object=axisLine"

	Execute "ModifyGizmo SETQUATERNION={-0.974948,0.122056,-0.008176,0.185790}"
	Execute "ModifyGizmo autoscaling=1"
	Execute "ModifyGizmo currentGroupObject=\"\""
	Execute "ModifyGizmo compile"
	Execute "ModifyGizmo endRecMacro"
End

//	End of plotting for raw data
// =========================================================================
// =========================================================================



// =========================================================================
// =========================================================================
//	Start of Arrays3d set panel

Function/T FillArrays3dParametersPanel(strStruct,hostWin,left,top)
	String strStruct									// optional passed value of xtal structure, this is used if passed
	String hostWin										// name of home window
	Variable left, top									// offsets from the left and top

	NewDataFolder/O root:Packages:micro				// ensure that the needed data folders exist
	NewDataFolder/O root:Packages:micro:Arrays3d
	SetWindow kwTopWin,userdata(Arrays3dPanelName)=hostWin+"#Arrays3dPanel"
	NewPanel/K=1/W=(left,top,left+221,top+500)/HOST=$hostWin
	ModifyPanel frameStyle=0, frameInset=0
	RenameWindow #,Arrays3dPanel

	Variable vert=30
	SetDrawLayer UserBack
	SetDrawEnv fname= "Lucida Grande",fsize= 14
	DrawText 30,vert,"Plotting Results"
	Button buttonPlotSlice3dArray,pos={19,vert},size={180,20},proc=Arrays3dButtonProc,title="Plot slice of 3D Array"
	Button buttonPlotSlice3dArray,help={"make a plot of one slice of a 3D Array"}
	Button buttonPlotSliceIntens3dArray,pos={19,vert+25},size={180,20},proc=Arrays3dButtonProc,title="Slice with Intensity, RGB"
	Button buttonPlotSliceIntens3dArray,help={"Plot slice of 3d array with intensity, an RGB image"}
	Button buttonPlotSliceOrients3dArray,pos={19,vert+50},size={180,20},proc=Arrays3dButtonProc,title="Slice of Orientations, RGB"
	Button buttonPlotSliceOrients3dArray,help={"Plot slice of 3d array showing Orientations, an RGB image"}
	Button buttonGraph9Components3dArray,pos={19,vert+75},size={180,20},proc=Arrays3dButtonProc,title="Show 9 Tensor Comps."
	Button buttonGraph9Components3dArray,help={"Graph of all 9 of the tensor components"}
	Button buttonReCalc9Components3dArray,pos={19,vert+100},size={180,20},proc=Arrays3dButtonProc,title="Re-Calc 9 Tensor Comps."
	Button buttonReCalc9Components3dArray,help={"Re-Calculate all 9 of the tensor components"}
	Button buttonReSetGNDconversion3dArray,pos={19,vert+125},size={180,20},proc=Arrays3dButtonProc,title="Re-Set GND conversion"
	Button buttonReSetGNDconversion3dArray,help={"Re-Set the GND conversion factor"}
	Button buttonPanelSlicer3dArray,pos={39,vert+150},size={160,20},proc=Arrays3dButtonProc,title="Show PanelSlicer"
	Button buttonPanelSlicer3dArray,help={"put up PanelSlicer, to control view of  slices"}

	vert=170+75
	SetDrawLayer UserBack
	SetDrawEnv fname= "Lucida Grande",fsize= 14
	DrawText 30,vert,"Read & Interpolate"
	Button buttonRead3dArrays,pos={19,vert},size={180,20},proc=Arrays3dButtonProc,title="Read In 3D data, 1st"
	Button buttonRead3dArrays,help={"Read 3D data from a file"}
	Button buttonInterp3dArrays,pos={19,vert+25},size={180,20},proc=Arrays3dButtonProc,title="Interpolate 3D data, 2nd"
	Button buttonInterp3dArrays,help={"interpolate the data just read in to a 3D volume"}

	vert=265+75
	SetDrawLayer UserBack
	SetDrawEnv fname= "Lucida Grande",fsize= 14
	DrawText 30,vert,"miscelaneous"
	Button buttonFindClosestPoint,pos={19,vert},size={180,20},proc=Arrays3dButtonProc,title="find closest point"
	Button buttonFindClosestPoint,help={"find point in original indexing list that is closest to the cursor"}
	Button buttonMarkRererencePoints,pos={19,vert+25},size={180,20},proc=Arrays3dButtonProc,title="mark reference point(s)"
	Button buttonMarkRererencePoints,help={"put markers on the graph to identify the reference point(s)"}
	PopupMenu popupTrim3dArrays,pos={19,vert+50},size={180,20},proc=Trim3dArraysPopMenuProc,title="Trim Top Surface & Bottom..."
	PopupMenu popupTrim3dArrays,help={"Trim Top Surface & Bottom from 3d arrays"}
	PopupMenu popupTrim3dArrays,mode=0,value= #"\"Both Top & Bottom;  Trim Only Top;  Trim Only Bottom\""
	Button buttonCutSide3dArrays,pos={19,vert+75},size={180,20},proc=Arrays3dButtonProc,title="Cut Off One Side"
	Button buttonCutSide3dArrays,help={"Cut Off One Side"}
	Button buttonColorTriangle3dArrays,pos={19,vert+100},size={180,20},proc=Arrays3dButtonProc,title="Show Cubic Color Triangle"
	Button buttonColorTriangle3dArrays,help={"Show Cubic Color Triangle"}
	PopupMenu popupTesting3dArrays,pos={69,vert+125},size={180,20},proc=testing3dArraysPopMenuProc,title="testing..."
	PopupMenu popupTesting3dArrays,help={"misc debugging & utility items"}
	PopupMenu popupTesting3dArrays,mode=0,value= #"\"Table of IndexingList;print a Rodriques Vector;Table of Raw Data;Make Test Data\""

	EnableDisableArrays3dControls(hostWin+"#Arrays3dPanel")
	return "#Arrays3dPanel"
End
//
Function EnableDisableArrays3dControls(win)			// here to enable/disable
	String win												// window (or sub window) to use
	Variable d

	d = strlen(WaveListClass("Interpolated3dArrays*","*",""))<1 ? 2 : 0
	Button buttonPlotSlice3dArray,win=$win,disable=d
	d = strlen(WaveListClass("Interpolated3dArrays*","*",""))<1 ? 2 : 0
	Button buttonPlotSliceIntens3dArray,win=$win,disable=d
	d = strlen(WaveListClass("Interpolated3dArrays*","*",""))<1 ? 2 : 0
	Button buttonPlotSliceOrients3dArray,win=$win,disable=d
	d = strlen(WaveListClass("TensorOrientationSlice","*",""))<1 ? 2 : 0
	Button buttonGraph9Components3dArray,win=$win,disable=d
	d = strlen(WaveListClass("OrientationSliceWave","*",""))<1 ? 2 : 0
	Button buttonReCalc9Components3dArray,win=$win,disable=d
	Button buttonReSetGNDconversion3dArray,win=$win,disable=0
	d = strlen(WaveListClass("OrientationSliceWave","*",""))<1 ? 2 : 0
	Button buttonPanelSlicer3dArray,win=$win,disable=d

	Button buttonRead3dArrays,win=$win,disable=0			// always OK to read in data
	d = strlen(WaveListClass("Random3dArrays*","*","DIMS:1"))<1 ? 2 : 0
	Button buttonInterp3dArrays,win=$win,disable=d
	d = strlen(WaveListClass("OrientationSliceWave*","*",""))<1 || exists("IndexLots")!=6 ? 2 : 0
	Button buttonFindClosestPoint,win=$win,disable= d
	d = strlen(WaveListClass("OrientationSliceWave*","*",""))<1 || exists("IndexLots")!=6 ? 2 : 0
	Button buttonMarkRererencePoints,win=$win,disable= d
	d = strlen(WaveListClass("Interpolated3dArrays*","*",""))<1 ? 2 : 0
	PopupMenu popupTrim3dArrays,win=$win,disable= d
	d = strlen(WaveListClass("Interpolated3dArrays*","*",""))<1 ? 2 : 0
	Button buttonCutSide3dArrays,win=$win,disable=d
	Button buttonColorTriangle3dArrays,win=$win,disable=0			// always OK to view the triangle
	PopupMenu popupTesting3dArrays,win=$win,disable= 0
End

Function Arrays3dButtonProc(ctrlName) : ButtonControl
	String ctrlName
	if (stringmatch(ctrlName,"buttonPlotSlice3dArray") && strlen(WaveListClass("Interpolated3dArrays*","*","")))
		MakeGraphOfSlice()
	elseif (stringmatch(ctrlName,"buttonPlotSliceIntens3dArray") && strlen(WaveListClass("Interpolated3dArrays*","*","")))
		MakeGraphOfSliceWaveRGB()
	elseif (stringmatch(ctrlName,"buttonPlotSliceOrients3dArray") && strlen(WaveListClass("Interpolated3dArrays*","*","")))
		MakeGraphOfSliceOrients()
	elseif (stringmatch(ctrlName,"buttonGraph9Components3dArray") && strlen(WaveListClass("TensorOrientationSlice","*","")))
		GraphTensorcomponents("")
	elseif (stringmatch(ctrlName,"buttonReCalc9Components3dArray") && strlen(WaveListClass("OrientationSliceWave","*","")))
		ReCalcAllTensorComponents("")
	elseif (stringmatch(ctrlName,"buttonReSetGNDconversion3dArray"))
		Change_GND_DislocationDensity(NaN)
	elseif (stringmatch(ctrlName,"buttonPanelSlicer3dArray") && strlen(WaveListClass("OrientationSliceWave","sliceWave","")))
		Wave sliceWave=sliceWave
		MakePanelSlicer(sliceWave)
#if (Exists("IndexLots")==6)
	elseif (stringmatch(ctrlName,"buttonFindClosestPoint") && strlen(WaveListClass("OrientationSliceWave*","*","")))
		findClosestPoint(NaN,NaN)
	elseif (stringmatch(ctrlName,"buttonMarkRererencePoints") && strlen(WaveListClass("OrientationSliceWave*","*","")))
		putMarkersAtRefPoints("FF")
#endif
	elseif (stringmatch(ctrlName,"buttonRead3dArrays"))
		readInXYZorientationsFromFile(NaN)
	elseif (stringmatch(ctrlName,"buttonInterp3dArrays") &&  strlen(WaveListClass("Random3dArrays*","*","DIMS:1")))
		InterpolateRodriquesFromXHF(NaN)
	elseif (stringmatch(ctrlName,"buttonCutSide3dArrays") &&  strlen(WaveListClass("Random3dArrays*","*","DIMS:1")))
		CutOffOneSide(NaN)
	elseif (stringmatch(ctrlName,"buttonColorTriangle3dArrays"))
		showCubicTriangleColors(NaN,NaN)
	endif
	EnableDisableArrays3dControls(GetUserData("microPanel","","Arrays3dPanelName"))
End

Function Trim3dArraysPopMenuProc(ctrlName,popNum,popStr) : PopupMenuControl
	String ctrlName
	Variable popNum
	String popStr
	Variable top, bot
	top = stringmatch(popStr,"Both Top & Bottom") || stringmatch(popStr,"  Trim Only Top")
	bot = stringmatch(popStr,"Both Top & Bottom") || stringmatch(popStr,"  Trim Only Bottom")
	top = exists("TrimTopSurface3d")==6 ? top : 0
	bot = exists("TrimBotSurface3d")==6 ? bot : 0
	if (top)
#if (Exists("TrimTopSurface3d")==6)
		TrimTopSurface3d()
#endif
	endif
	if (bot)
#if (Exists("TrimTopSurface3d")==6)
		TrimBotSurface3d(NaN)
#endif
	endif
	EnableDisableArrays3dControls(GetUserData("microPanel","","Arrays3dPanelName"))
End

Function testing3dArraysPopMenuProc(ctrlName,popNum,popStr) : PopupMenuControl
	String ctrlName
	Variable popNum
	String popStr

	if (stringmatch(popStr,"Table of IndexingList"))
		TableOfIndexingList($"")
	elseif (stringmatch(popStr,"print a Rodriques Vector") && strlen(WaveListClass("Interpolated3dArrays*","*","")))
		RodriquesVector2str(NaN,NaN,NaN)
	elseif (stringmatch(popStr,"Table of Raw Data") && strlen(WaveListClass("Random3dArrays*","*","")))
		TableOfRawData()
	elseif (stringmatch(popStr,"Make Test Data"))
		MakeTestData("")
	endif
	EnableDisableArrays3dControls(GetUserData("microPanel","","Arrays3dPanelName"))
End

//	End of Arrays3d set panel
// =========================================================================
// =========================================================================




// ==================================================================
// =============== Start of utility section, very general stuff =====================

// find tagged numerical value in file, and return it
Function/T fetchAndAddNumber2List(refNum,tagStr,list,scaleFactor)
	Variable refNum
	String tagStr
	String list
	Variable scaleFactor			// a scale factor to apply to numbers, if NaN, use 1
	scaleFactor = numtype(scaleFactor) ? 1 : scaleFactor

	Variable value
	FSetPos refNum,0
	value = str2num(findTagInFile(refNum,tagStr))*scaleFactor
	if (numtype(value)!=2)
		list = ReplaceNumberByKey(tagStr,list,value,"=")
	endif
	return list
End


// find tagged string in file, and return it
Function/T fetchAndAddSring2List(refNum,tagStr,list)
	Variable refNum
	String tagStr
	String list

	String str
	FSetPos refNum,0
	str = findTagInFile(refNum,tagStr)
	if (strlen(str)>1 || (strlen(str)==1 && char2num(str[0])!=32))
		list = ReplaceStringByKey(tagStr,list,str,"=")
	endif
	return list
End


// return as a string the tagged value in file
// this is used by fetchAndAddNumber2List() and fetchAndAddSring2List()
Function/T findTagInFile(fileRef,tagStr)	// returns the 'value' string associated with the tag
	Variable fileRef				// ref to an opened file
	String tagStr				// the tag to search for

	String line
	tagStr = "$"+tagStr
	Variable i,c, tagLen=strlen(tagStr)
	do
		FReadLine fileRef, line
		i = strsearch(line,tagStr,0,2)
		i = (i==0 && char2num(line[tagLen])>32) ? -1 : i		// in case tag is a substring of line
	while (i && strlen(line)>0)

	if (strlen(line)<1)								// end of file
		return ""
	endif

	// find index to first useful character
	i = strlen(tagStr)-1
	do
		i += 1
		c = char2num(line[i])
	while (c<=32 || c==char2num("="))
	line = line[i,Inf]

	if (char2num(line[0])==char2num("'"))		// a quoted string
		i = strsearch(line,"'",1)
		if (i<0)
			Abort "encountered a partially quoted string in tagged data file"
		endif
		line = line[1,i-1]
		return line
	endif

	// find end of value
	i = strsearch(line,"//",0)						// remove possible comment
	line = SelectString(i,line," ",line[0,i-1])
	// i<0		no comment
	// i=0		only a comment, so use 1 space
	// i>0		comment found, remove it

	// remove trailing white space
	i = strlen(line)
	do
		i -= 1
		c = char2num(line[i])
	while(c<=32 && i>0)
	line = line[0,i]
	line = SelectString(strlen(line)==0,line," ")	// empty string means EOF, so make contents at least one space
	return line
End


Static Function normalize(a)	// normalize a and return the initial magnitude, also in Utility_JZT.ipf
	Wave a
	Variable norm_a = norm(a)
	a /= norm_a
	return norm_a
End

// ================ End of utility section, very general stuff =====================
// ==================================================================





// *************************************************************************************
// *************************************************************************************
// *************************************************************************************
// ********************************* Start of Grain Math Items *********************************

//	moved to microGeometry.ipf
//
//Function rotationAngleOfMat(rot)				// returns the total rotation angle of a matrix 'rot'
//	Wave rot									// the rotation matrix
//	Variable trace = MatrixTrace(rot)			// trace = 1 + 2*cos(theta)
//	Variable cosine = (trace-1)/2				// cosine of the rotation angle
//	cosine = (cosine>1) ? (2-cosine) : cosine
//	return acos(cosine)*180/PI				// rotation angle in degrees
//End


//	moved to microGeometry.ipf
//
//// compute angle and axis of a rotation matrix
//// Aug 2, 2007, this was giving the wrong sign for the rotation, so I reversed the "curl" in defn of axis.  JZT
////		changed "axis[0] = rr[1][2] - rr[2][1]"   -->   "axis[0] = rr[2][1] - rr[1][2]"
//Function axisOfMatrix(rot,axis)
//	// returns total rotation angle, and sets axis to the axis of the total rotation
//	Wave rot										// rotation matrix
//	Wave axis										// axis of the rotation (angle is returned)
//
//	Variable sumd = notRotationMatrix(rot)	// accept positive sumd that are less than 1e-4
//	if (sumd<0 || sumd>1e-4)
//		DoAlert 0, "'"+NameOfWave(rot)+"' is not a rotation matrix in axisOfMatrix()"
//		axis = NaN
//		return NaN
//	elseif (0<sumd)								// close enough to a roation matix, but tidy it up first
//		Make/N=(3,3)/O/D axisMat__
//		axisMat__ = rot
//		if (SquareUpMatrix(axisMat__))
//			DoAlert 0, "cannot square up '"+NameOfWave(rot)+"' in axisOfMatrix()"
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
//	if (cosine<= -1)								// special for 180° rotation,
//		axis[0] = sqrt((rr[0][0]+1)/2)
//		axis[1] = sqrt((rr[1][1]+1)/2)
//		axis[2] = sqrt((rr[2][2]+1)/2)		// always assume z positive
//		axis[0] = (rr[0][2]+rr[2][0])<0 ? -axis[0] : axis[0]
//		axis[1] = (rr[1][2]+rr[2][1])<0 ? -axis[1] : axis[1]
//	else											// rotaion < 180°, usual formula works
//		axis[0] = rr[2][1] - rr[1][2]
//		axis[1] = rr[0][2] - rr[2][0]
//		axis[2] = rr[1][0] - rr[0][1]
//		axis /= 2
//	endif
//	normalize(axis)
//	KillWaves /Z axisMat__
//	return acos(cosine)*180/PI					// rotation angle in degrees
//End
////	else											// rotaion < 180°, usual formula works
////		axis[0] = rr[1][2] - rr[2][1]
////		axis[1] = rr[2][0] - rr[0][2]
////		axis[2] = rr[0][1] - rr[1][0]
////	endif


//	moved to microGeometry.ipf
//
//Function AngleBetweenRotationVectors(R1,R2)	// returns rotation between orientations defined by Rodriques vectors (degrees)
//	Wave R1,R2								// two Rodriques vectors
//
//	Make/N=(3,3)/O/D mat1_JZT, mat2_JZT
//	Wave mat1=mat1_JZT, mat2=mat2_JZT
//	rotationMatAboutAxis(R1,NaN,mat1)
//	rotationMatAboutAxis(R2,NaN,mat2)
//	MatrixOp/O temp_angle_11_JZT = Trace(mat2 x mat1^t)
//	Variable cosine = (temp_angle_11_JZT[0][0]-1) / 2
//	cosine = min(max(-1,cosine),1)			// ensure cosine in range [-1,1]
//	return acos(cosine)*180/PI
////	KillWaves/Z temp_angle_11_JZT, mat1_JZT, mat2_JZT
//End


//	moved to microGeometry.ipf
//
//// set mat to be a rotation matrix about axis with angle
//Function rotationMatAboutAxis(axis,angle,mat)
//	Wave axis				// axis about which to rotate (or possibly Rodriques vector)
//	Variable angle			// angle to rotate (degrees), assumes axis is true Rodriques vector if angle invalid
//	Wave mat				// desired rotation matrix
//
//	Variable len = norm(axis)
//	angle = numtype(angle) ? 2*atan(len) : angle*PI/180	// the rotation angle (rad)
//
//	if (angle==0)			// zero angle rotation is just the identity matrix
//		mat = (p==q)
//		return 0
//	endif
//
//	Variable nx=axis[0]/len, ny=axis[1]/len, nz=axis[2]/len
//	Variable cosa=cos(angle), sina=sin(angle)
//	Variable c1 = 1-cosa
//	if (DimSize(mat,0)!=3 && DimSize(mat,1)!=3)
//		Abort "in rotationMatAboutAxis(), mat must be (3,3)"
//	endif
////	Redimension/N=(3,3) mat
//	// from		http://mathworld.wolfram.com/RodriguesRotationFormula.html (I double checked this too.)
//	mat[0][0] = cosa+nx*nx*c1;			mat[0][1] =nx*ny*c1-nz*sina;			mat[0][2] = nx*nz*c1+ny*sina;
//	mat[1][0] = nx*ny*c1+nz*sina;		mat[1][1] = cosa+ny*ny*c1;			mat[1][2] =ny*nz*c1-nx*sina;
//	mat[2][0] = nx*nz*c1-ny*sina;		mat[2][1] = ny*nz*c1+nx*sina;		mat[2][2] = cosa+nz*nz*c1;
//	return 0
//End


// rotate mat about x, y, or z by angle, not a general axis
Function rotateMat(mat,axis,angle)
	Wave mat				// matrix or vector to rotate
	Variable axis			// axis to rotate about 0=x, 1=y, 2=z
	Variable angle			// angle to rotate (radian)

	Variable cosa = cos(angle)
	Variable sina = sin(angle)

	String wName=UniqueName("mat33_",1,0)
	if (NumberByKey("NUMTYPE",WaveInfo(mat,0))==2)	// 2 means single precision
		Make/N=(3,3)/O $wName
	else
		Make/N=(3,3)/O/D $wName						// if mat not single float, use double
	endif
	Wave rot=$wName
	rot = (p==q)

	switch(axis)
		case 0:			// rotate about x-axis
			rot[1][1] = cosa
			rot[2][2] = cosa
			rot[1][2] = -sina
			rot[2][1] =   sina
			break
		case 1:			// rotate about y-axis
			rot[0][0] = cosa
			rot[2][2] = cosa
			rot[0][2] =   sina
			rot[2][0] = -sina
			break
		case 2:			// rotate about z-axis
			rot[0][0] = cosa
			rot[1][1] = cosa
			rot[0][1] = -sina
			rot[1][0] =   sina
			break
		default:
			rot = NaN
	endswitch

	MatrixOp/O mat = rot x mat								// apply the rotation matrix to mat
	Killwaves/Z rot
End
//
//Function test(axis)
//	Variable axis			// axis to rotate about 0=x, 1=y, 2=z
//
//	Make/N=3/O gg={0,1,0}
//	rotateMat(gg,axis,10*PI/180)
//	printWave(gg)
//End

//	moved to microGeometry.ipf
//
//Function notRotationMatrix(mat)			// false if mat is a rotation matrix, i.e. mat^T = mat^-1
//	Wave mat
//	String wName = UniqueName("rot",1,0)
//	MatrixOp/O $wName = Abs(( mat x mat^t ) - Identity(3))
//	Wave diff = $wName
//	Variable sumd = sum(diff)
//	// printWave(diff)
//	// print "sum(diff) = ",sumd
//	KillWaves/Z $wName
//
////	Variable returnVal, thresh=2e-10
//	Variable returnVal, thresh=2e-5
//	if (sumd>thresh)
//		returnVal = sumd						// not a rotation matrix
//	elseif (MatrixDet(mat)<0)
//		returnVal = -1						// an improper rotation matrix
//	else
//		returnVal = 0							// yes it is a rotation matrix
//	endif
//	if (ItemsInList(GetRTStackInfo(0))<2)
//		printf "'%s'  is %s matrix\r",NameOfWave(mat),SelectString(returnVal ,"an improper rotation","a rotation","NOT a rotation")
//		if (sumd>=thresh && sumd<1e-4)
//			printf "	but, it is almost a rotation matrix, sum(|diff|) = %g\r",sumd
//		endif
//	endif
//	return returnVal
//End


// ********************************** End of Grain Math Items *********************************
// *************************************************************************************
// *************************************************************************************
// *************************************************************************************