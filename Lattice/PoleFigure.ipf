#pragma rtGlobals=1		// Use modern global access method.
#include "vector-math"
#pragma version = 1.01

// add the following to the "resizeHookFunc" so that things appear square even without width={Aspect,1}
//	GetWindow kwTopWin psize
//	Variable width = V_right-V_left,  height=V_bottom-V_top
//	// width and height of plottable area,  I need widht/height to correctly scale the axes


//		tan(delta/2) = x
//		where delta = polar (90�-lattitude)
//		x is position in stereographic projection plane
//		at maximum delta of 90�, the value is 1


//	here is the way tog get rotation axis and angle.
//	From rotation matrices A, B, you get C=A*B^(-1)
//	The angle = acos((C11+C22+C33-1)/2)
			// cos(angle) = { Trace(A*Binverse)-1 } / 2
//	axis direction (x,y,z): x=C23-C32, y=C31-C13, z=C12-C21
//	Make sure to normalize A and B first.
//		Wenge


Menu "PoleFigure"
	"Display 001 Pole Figure", DisplayPoleFigure(0,0,1,"orientations")
	help = {"Display a 001 pole figure stereographic projection based on matricies in 'orientations'"}
	"Display 011 Pole Figure", DisplayPoleFigure(0,1,1,"orientations")
	help = {"Display a 011 pole figure stereographic projection based on matricies in 'orientations'"}
	"Display 111 Pole Figure", DisplayPoleFigure(1,1,1,"orientations")
	help = {"Display a 111 pole figure stereographic projection based on matricies in 'orientations'"}
	"Display hkl Pole Figures", DisplayPoleFigure(NaN,NaN,NaN,"orientations")
	help = {"Display an hkl pole figure stereographic projection based on matricies in 'orientations'"}
	"Show the rotation axis table", ShowAxisTable()
	help = {"Show a table of the rotation axis for the current pole figure'"}
	"-"
	"Make StepWise Rotation Angles", MakeStepWiseRotationMats(orientations,"rmat") ; Make_rotationAngles()
	help = {"calculate grain to grain rotation angles, store in a wave called 'rotAngles'"}
	"  add Signs to Total Rotation Angles", MakeRotationAnglesSigned(NaN,NaN,NaN)
	help = {"makes signed grain to grain rotation angles, sign is referenced to hkl direction of orientations[][][0]."}
	"  make Total Rotation Angles Unsigned", rotAngles = abs(rotAngles)
	help = {"change 'rotAngles' to absolute values"}
	"Display Rotation Angles", Display rotAngles  ;  Graph_rotationAnglesStyle()
	help = {"display the 'rotAngles' wave"}
	"Make Total Rotations from one Ref Matrix", MakeRotationAnglesFromReference(orientations,NaN)
	help = {"Make a wave of the total rotation from one specified reference matrix"}
	"Make Wave of Rotation Angle Components", MakeRotationAnglesComponent(NaN,NaN,NaN,NaN,NaN)
	help = {"make a wave of the component of the rotation angle in the xyz direction"}
	Submenu "only needed once"
		"Read in all the Matricies", ReadInManyMatricies(10,300,"")
		help = {"read in a group of matricies from a file, and store all of them in 'orientations'"}
		"Make and Apply Swap Vectors",MakeSwapVectors(orientations,"swaps") ; ApplySwapsToOrients(orientations,swaps)
		help = {"make swap vectors for rotating orientation matricies, and then apply them"}
		"Reset Matricies to as-read values", orientations=orientationsOrig
		help = {"reset the matricies to values as read in (before swapping)"}
	End
End



// from the email to Budai and Ben
//
//	this is the final one
//
//	In the .IND file, I take the reciprocal lattice vectors from the part:
//
//	matrix h  k  l -> X  Y  Z
//	  -2.28186  -0.54961   0.76759
//	   0.93727  -1.55936   1.66974
//	   0.11308   1.83424   1.64952
//
//	I assume that:
//	a* = (-2.28186,  0.93727, 0.11308)
//	b* = (-0.54961, -1.55936, 1.83424)
//	c* = ( 0.76759,  1.66974, 1.64952)
//
//
//	These get shuffled to:
//	a* = (-2.28186,  0.93727, 0.11308)
//	b* = ( 0.76759,  1.66974, 1.64952)
//	c* = (-0.54961, -1.55936, 1.83424)
//
//
//	From the file that I send "data3D.txt"
//
//	  a*(x)     a*(y)     a*(z)   b*(x)     b*(y)    b*(z)    c*(x)     c*(y)     c*(z)
//	  or00      or10      or20    or01      or11     or21     or02      or12      or22
//	2.28186  -0.93727  -0.11308  0.76759  1.66974  1.64952  -0.54961  -1.55936  1.83424
//
//
//	so that 
//	a* = (or00, or10, or20)
//	b* = (or01, or11, or21)
//	c* = (or02, or12, or22)
//
//
//	Note that there is a swap between c* and b*.  This is to get all of the orientaions into the smallest 
//	region, and is correct.  When using this data, use a* = (or00, or01, or02).







// ============================================================================
// ========  this section is for checking that each matrix is the shortest 90� rotation from the next =======


Function MakeRotationAnglesComponent(type,x,y,z,index)
	// make a wave of the component of the rotation angle in the xyz direction
	// this can be done in four ways depending on 'type'
	Variable type	// type=1	total rotation from indexth grain, about (x*,y*,z*) component
					// type=2	total rotation from indexth grain, about (h,k,l) component
					// type=3	grain to grain rotation about (h,k,l) component
	Variable x,y,z	// components of reference direction
	Variable index	// index into orientations to select reference matrix, index<0 flags use of global ref matrix

	Variable needPrompt=0
	if ( (type==1 || type==2 || numtype(type)) && (numtype(x)||numtype(y)||numtype(z)) )	// need reference direction information
		x = numtype(x) ? 1 : x
		y = numtype(y) ? 0 : y
		z = numtype(z) ? 0 : z
		needPrompt = 1
	endif
	if ( numtype(type) || type<1 || type >3)
		type = numtype(type) ? 1 : type
		type = min(3,max(1,type))
		needPrompt = 1
	endif
	index = numtype(index) ? 0 : index			// zero is a good default for this
	if (needPrompt)
		String reference
		sprintf reference "%d %d %d",x,y,z
		Prompt type, "type of rotation angles to make", popup, "total rotation from indexth grain, about (x*,y*,z*) component;total rotation from indexth grain, about (h,k,l) component;grain to grain rotation about (h,k,l) component"
		DoPrompt "type of rotation angle to calculate", type
		if (V_flag)
			return 1
		endif
		if (type==1)
			Prompt reference, "reference (x* y* z*)"
			Prompt index, "point index of reference orientation matrix (index<0 flags use of global ref matrix)"
			DoPrompt "total rotation from indexth grain, about (x*,y*,z*) component", reference,index
		elseif (type==2)
			Prompt reference, "reference (hkl)"
			Prompt index, "point index of refernece orientaion matrix (index<0 flags use of global ref matrix)"
			DoPrompt "total rotation from indexth grain, about (h,k,l) component", reference,index
		elseif (type==3)
			Prompt reference, "reference (hkl)"
			DoPrompt "grain to grain rotation about (h,k,l) component", reference
		else
			Abort "how did I get here"
		endif
		if (V_flag)
			return 1
		endif
		sscanf reference, "%g %g %g", x,y,z
		printf "  MakeRotationAnglesComponent(%d,%g,%g,%g,%d)\r",type,x,y,z,index
	endif

	String angleName							// name of wave to make and fill
	String sx=num2istr(abs(x)), sy=num2istr(abs(y)), sz=num2istr(abs(z))
	if (sign(x)<0)
		sx = "n"+sx
	endif
	if (sign(y)<0)
		sy = "n"+sy
	endif
	if (sign(z)<0)
		sz = "n"+sz
	endif

	Make/O ref_axis__ = {x,y,z}
	Make/O/N=3 axis___
	Make/N=(3,3)/O mat__, mat_ref__, rot_temp__
	Wave orientations=orientations

	if (type==1)								// total from reference (x*,y*,z*)
		sprintf angleName, "angles_xyz_%s_%s_%s" sx,sy,sz
	elseif (type==2)
		if (index<0)
			Wave tmat=::ref_mat0				// index<0 flags use to global ref matrix
			mat__ = tmat
		else
			mat__ = orientations[p][q][index]
		endif
		ref_axis__ = {x,y,z}					// convert hkl to (x*,y*,z*)
		MatrixMultiply mat__, ref_axis__
		Wave M_product=M_product
		ref_axis__ = M_product
		sprintf angleName, "angles_hkl_%s_%s_%s" sx,sy,sz
	else
		sprintf angleName, "angles_pt2pt_%s_%s_%s" sx,sy,sz
	endif

	if (type==1 || type==2)
		if (index<0)
			Wave tmat=::ref_mat0				// index<0 flags use to global ref matrix
			mat_ref__ = tmat
		else
			mat_ref__ =  orientations[p][q][index]	// reference orientation matrix (not point to point)
		endif
	endif
	normalize(ref_axis__)						// normalize reference axis
	Variable N=DimSize(orientations,2)		// length of wave angleName
	Make/N=(N)/O $angleName
	print "  created wave of angles named '"+angleName+"'"
	Wave angles=$angleName
	Wave rmat=rmat
	Variable i=0
	do
		if (type==3)							// for point to point, hkl is know, xyz is from each orientation
			mat__ = orientations[p][q][i]		// to get reference axis
			ref_axis__ = {x,y,z}				// convert hkl to (x*,y*,z*)
			MatrixMultiply mat__, ref_axis__
			Wave M_product=M_product
			ref_axis__ = M_product
			normalize(ref_axis__)				// normalize reference axis
			mat__ = rmat[p][q][i]				// the rotation matrix
		else										// for types 1 & 2, ref_axis__ is already known
			mat__ = orientations[p][q][i]
			MakeRotationMatFromAB(mat_ref__,mat__,rot_temp__)	// rotation mat between ref and orientations[i]
			mat__ = rot_temp__
		endif
		angles[i] = axisOfMatrix(mat__,axis___)	// compute angle and axis of a rotation matrix
		normalize(axis___)
		angles[i] *= MatrixDot(ref_axis__,axis___)
		i += 1
	while (i<N)
	Killwaves/Z ref_axis__, axis___, mat__, mat_ref__, rot_temp__, M_product
	WaveStats/Q angles
//	printf "range of rotation angles is [%g�,%g�],  with %d bad points\r",V_min,V_max,V_numNaNs
End


Function MakeRotationAnglesComponentOLD(type,x,y,z,index)
	// make a wave of the component of the rotation angle in the xyz direction
	// this can be done in three ways depending on 'type'
	// type=1	rotation from index, (x*,y*,z*) component
	// type=2	rotation from index (hkl) component
	// type=3	point to point rotation (do not need x*,y*,z*)
	Variable type
	Variable x,y,z				// components of reference direction
	Variable index				// index into orientations to select reference orientation

	Wave orientations=orientations

	Variable needPrompt=0
	if ( !(type==3) && (numtype(x)||numtype(y)||numtype(z)) )		// need reference information
		x = numtype(x) ? 1 : x
		y = numtype(y) ? 0 : y
		z = numtype(z) ? 0 : z
		needPrompt = 1
	endif
	if ( numtype(type) || type<1 || type >3)
		type = numtype(type) ? 1 : type
		type = min(3,max(1,type))
		needPrompt = 1
	endif
	index = numtype(index) ? 0 : index			// zero is a good default for this
	if (needPrompt)
		String reference
		sprintf reference "%d %d %d",x,y,z
		Prompt type, "type of rotation angles to make", popup, "rotation from index, (x*,y*,z*) component;rotation from index (hkl) component;point to point rotation (do not need index)"
		Prompt index, "point index (used for hkl and to give starting orientation)"
		Prompt reference, "reference (x* y* z*) or (hkl) if not point to point"
		DoPrompt "rotation angle reference", type,reference,index
		if (V_flag)
			return 1
		endif
		sscanf reference, "%g %g %g", x,y,z
		printf "  MakeRotationAnglesComponent(%d,%g,%g,%g,%d)\r",type,x,y,z,index
	endif

	String angleName							// name of wave to make and fill
	String sx=num2istr(abs(x)), sy=num2istr(abs(y)), sz=num2istr(abs(z))
	if (sign(x)<0)
		sx = "n"+sx
	endif
	if (sign(y)<0)
		sy = "n"+sy
	endif
	if (sign(z)<0)
		sz = "n"+sz
	endif

	Make/O ref_axis__ = {x,y,z}
	Make/O/N=3 axis___
	Make/N=(3,3)/O mat__, mat_ref__, rot_temp__

	if (type==1)								// total from reference (x*,y*,z*)
		sprintf angleName, "angles_xyz_%s_%s_%s" sx,sy,sz
	elseif (type==2)
		mat__ = orientations[p][q][index]
		ref_axis__ = {x,y,z}					// convert hkl to (x*,y*,z*)
		MatrixMultiply mat__, ref_axis__
		Wave M_product=M_product
		ref_axis__ = M_product
		sprintf angleName, "angles_hkl_%s_%s_%s" sx,sy,sz
	else
		sprintf angleName, "angles_pt2pt_%s_%s_%s" sx,sy,sz
	endif

	if (type==1 || type==2)
		mat_ref__ =  orientations[p][q][index]	// reference orientation matrix (not point to point)
	endif
	normalize(ref_axis__)						// normalize reference axis
	Variable N=DimSize(orientations,2)		// length of wave angleName
	Make/N=(N)/O $angleName
	print "  created wave of angles named '"+angleName+"'"
	Wave angles=$angleName
	Wave rmat=rmat
	Variable i=0
	do
		if (type==3)							// for point to point, hkl is know, xyz is from each orientation
			mat__ = orientations[p][q][i]		// to get reference axis
			ref_axis__ = {x,y,z}				// convert hkl to (x*,y*,z*)
			MatrixMultiply mat__, ref_axis__
			Wave M_product=M_product
			ref_axis__ = M_product
			normalize(ref_axis__)				// normalize reference axis
			mat__ = rmat[p][q][i]				// the rotation matrix
		else										// for types 1 & 2, ref_axis__ is already known
			mat__ = orientations[p][q][i]
			MakeRotationMatFromAB(mat_ref__,mat__,rot_temp__)	// rotation mat between ref and orientations[i]
			mat__ = rot_temp__
		endif
		angles[i] = axisOfMatrix(mat__,axis___)	// compute angle and axis of a rotation matrix
		normalize(axis___)
		angles[i] *= MatrixDot(ref_axis__,axis___)
		i += 1
	while (i<N)
	Killwaves/Z ref_axis__, axis___, mat__, mat_ref__, rot_temp__, M_product
	WaveStats/Q angles
	printf "range of rotation angles is [%g�,%g�],  with %d bad points\r",V_min,V_max,V_numNaNs
End


Function MakeRotationAnglesFromReference(refMatrix,index)
	// make a wave of the total rotation from one specified reference matrix
	Wave refMatrix
	Variable index	// index into refMatrix to select reference matrix, only used if refMatrix is 3d

	if ((DimSize(refMatrix,2)>0 && index<0) || DimSize(refMatrix,1)!=3 || numtype(index))	// need a prompt
		String fldrSav=GetDataFolder(1)
		String list = "orientations;"
		if (DataFolderExists(":::"))
			SetDataFolder :::
			list += WaveList("*",";","DIMS:2")
			// String list = "orientations;"+WaveList("*",";","DIMS:2")
			SetDataFolder fldrSav
		endif
		String refName="orientations"
		Prompt refName, "name of reference matrix", popup, list
		DoPrompt "pick reference orientation wave", refName
		if (V_flag)
			return 1
		endif
		if (!stringmatch(refName,"orientations"))
			refName = ":::"+refName
		endif
		Wave refMatrix=$refName
		if (!WaveExists(refMatrix))
			DoAlert 0, "reference wave '"+NameOfWave(refMatrix)+"' does not exist"
			return 1
		endif
		if (DimSize(refMatrix,2)>0)
			index = numtype(index) ? 0 : index
			Prompt index, "point index of reference orientation matrix"
			DoPrompt "index into "+NameOfWave(refMatrix)+" for reference matrix",index
			if (V_flag)
				return 1
			endif
		endif
		if (index<0 || index>=DimSize(refMatrix,2))
			DoAlert 0, "Illegal index, outside of wave "+refName
			return 1
		endif
		printf "  MakeRotationAnglesFromReference(%s,%d)\r",refName,index
	endif
	index = (DimSize(refMatrix,2)==0) ? -1 : index

	String angleName="angles_ref" 					// name of wave to make and fill
	Make/N=(3,3)/O mat__, mat_ref__
	if (index<0)
		mat_ref__ = refMatrix
	else
		mat_ref__ = refMatrix[p][q][index]
	endif

	Wave orientations=orientations
	Variable N=DimSize(orientations,2)		// length of wave angleName
	Make/N=(N)/O $angleName
//	print "  created wave of angles named '"+angleName+"'"
	Wave angles=$angleName
	Wave rmat=rmat
	Variable i=0
	do
		mat__ = orientations[p][q][i]
		angles[i] = AngleBetweenMats(mat__,mat_ref__)
		i += 1
	while (i<N)
	Killwaves/Z mat__, mat_ref__, M_product
//	WaveStats/Q angles
//	printf "range of rotation angles is [%g�,%g�],  with %d bad points\r",V_min,V_max,V_numNaNs
End


Function ApplySwapsToOrients(orients,swaps)
	Wave orients										// matrix of all of the orientation matricies
	wave swaps											// matrixs with the individual swap vectors
	// applys the swaps in 'swaps' to the set of matricies in 'orients'

	Variable N=DimSize(orients,2)
	if (DimSize(orients,0)!=3 || DimSize(orients,1)!=3 || N<1)
		Abort "ApplySwapsToOrients, orients matrix is illegal size"
	endif
	if (DimSize(swaps,0)!=3 || DimSize(swaps,1)!=N)
		Abort "ApplySwapsToOrients, swaps matrix is illegal size"
	endif

	Make/N=(3,3)/O mat__
	Variable i,m
	for(i=0;i<N;i+=1)										// loop through all orientation matricies
		mat__ = orients[p][q][i]
		m = SwapXYZ2xyz(mat__,swaps[0][i],swaps[1][i],swaps[2][i])	// swap around mat__
		if (!m)
			print "bad swap vectors at point "+num2str(i)
			Abort "bad swap vectors at point "+num2str(i)
		endif
		orients[][][i] = mat__[p][q]						// put swapped matrix back into orients
	endfor
	Killwaves/Z mat__
End


Function MakeSwapVectors(orients,Sswaps)
	Wave orients										// matrix of all of the orientation matricies
	String Sswaps										// name of matrixs to get the individual swap vectors
	// From the orientation matricies, make swap vectors that will swap a*, b*, c* so that they are all
	// as close together as possible (less than 90�)

	Variable lastGood=3
	Make/N=(3,3)/O orientRef__
	orientRef__ = orients[p][q][lastGood]				// the orientation matrix to test against
	MakeSwapVectorsReference(orients,Sswaps,orientRef__)
	Killwaves/Z orientRef__
End


Function MakeSwapVectorsReference(orients,Sswaps,orientRef)
	Wave orients										// matrix of all of the orientation matricies
	String Sswaps										// name of matrixs to get the individual swap vectors
	Wave orientRef										// reference orientation matrix
	// From the orientation matricies, make swap vectors that will swap a*, b*, c* so that they are all
	// as close together as possible (less than 90�)

	Variable N=DimSize(orients,2)
	if (DimSize(orients,0)!=3 || DimSize(orients,1)!=3 || N<1)
		Abort "MakeSwapVectorsReference, orients matrix is illegal size"
	endif
	if (!WaveExists(orientRef))
		Abort "no refrence orientation matrix given"
	endif

	Variable i,j,k,eps
	Variable i0=1,j0=2,k0=3
	Variable index, minVal, a

	Make/N=(3,N)/O $Sswaps
	Wave swaps = $Sswaps
	for(i=0;i<N;i+=1)
		swaps[0,2][i]={1,2,3}								// set each swap to default of {1,2,3}
	endfor
	Make/N=(3,3)/O mat__, rmat__
	for(index=0;index<N;index+=1)						// loop through all orientation matricies
		mat__ = orients[p][q][index]
		minVal=1e20
		for(k=-3;k<4;k+=1)								// check all 48 of the 90� type rotations
			for(j=-3;j<4;j+=1)
				for(i=-3;i<4;i+=1)
					eps = epsilon_ijk(i,j,k)				// one possible swapping
					if (eps==1)							// only 24 right handed swappings are checked
						rmat__ = mat__
						SwapXYZ2xyz(rmat__,i,j,k)		// rmat__ is now swapped version of mat__
						a = AngleBetweenMats(rmat__,orientRef)
						if (a<minVal)
							minVal = a
							i0=i  ;  j0=j  ;  k0=k
						endif
					endif
				endfor
			endfor
		endfor
		swaps[0,2][index]={i0,j0,k0}					// save best swap for this index
	endfor
	Killwaves/Z mat__, rmat__, M_product
End



Function MakeStepWiseRotationMats(orientations,Srmat)	// makes rmat based on orientations
	// calculate grain to grain rotation matricies, all rot matricies in a wave named by Srmat
	Wave orientations
	String Srmat						// name of wave to recieve stepwise rotation matricies
	// makes rmat, which is the  rotation matricies between succesive grains
	// rmat[i] is the rotation needed to go from orient[i-1] to orient[i]
	Duplicate/O orientations $Srmat
	Wave rmat=$Srmat					// wave to recieve stepwise rotation matricies

	Make/N=(3,3)/O amat,bmat, lastGood__
	Variable det							// value of determinant
	Variable N=DimSize(orientations,2)
	rmat = NaN

	rmat[0,2][0,2][0] = p==q			// first one is always identity matrix
	Variable i=1
	lastGood__ = orientations[p][q][0]
	do
		if (!numtype(orientations[0][0][i-1]))
			lastGood__ = orientations[p][q][i-1]	// last orientation matrix that was not NaN
			amat = lastGood__
		endif
		amat = lastGood__
		bmat = (p==q)
		MatrixLLS /O/Z  amat bmat		// bmat is now a^-1
		amat = bmat					// amat is now a^-1
		bmat = orientations[p][q][i]
		MatrixMultiply bmat,amat		// rotation matrix is now in M_product
		Wave M_product=M_product

		det = MatrixDet(M_product)
		if (det>1)
			M_product /= det
		endif

		rmat[][][i] = M_product[p][q]
		i += 1
	while (i<N)
	Killwaves/Z amat,bmat,M_product, lastGood__
EndMacro


Function MakeRotationMatFromAB(a,b,rot)	// makes rotation mat between a and b
	Wave a,b							// input matricies
	Wave rot							// rotation matrix

	// makes rmat, which is the  rotation matricies between succesive grains
	// rmat[i] is the rotation needed to go from orient[i-1] to orient[i]

	if (numtype(a[0][0]) || numtype(b[0][0]))
		rot = NaN
		return 1
	endif

	Make/N=(3,3)/O temp__, invA_temp__
	rot = a
	temp__ = (p==q)
	MatrixLLS /O/Z  rot temp__		// temp__ is now a^-1
	invA_temp__ = temp__				// invA_temp__ is now a^-1
	MatrixMultiply b,invA_temp__		// rotation matrix is now in M_product
	Wave M_product=M_product

	Variable det = MatrixDet(M_product)
	if (det>1)
		M_product /= det
	endif
	rot = M_product
	Killwaves/Z temp__,M_product, invA_temp__
	return 0
EndMacro




Function SwapXYZ2xyz(mat,x2,y2,z2)		// permute axes so X->x, Y->y, Z->z
	Wave mat				// matrix to permute (this is NOT a rotation)
	Variable x2,y2,z2		// destinations of x,y,z (use �1, �2, and �3 (cannot use 0, since 0==�0)

	Variable eps = epsilon_ijk(x2,y2,z2)
	if (!eps)
		DoAlert 0, "SwapXYZ2xyz x2,y2,z2, must be a permutation of �1,2, or 3"
//		Abort "SwapXYZ2xyz x2,y2,z2, must be a permutation of �1,2, or 3"
	endif

	Make/N=(3,3)/O axses__
	axses__ = mat
	mat[][0] = axses__[p][abs(x2)-1] * sign(x2)
	mat[][1] = axses__[p][abs(y2)-1] * sign(y2)
	mat[][2] = axses__[p][abs(z2)-1] * sign(z2)

	Killwaves/Z axses__
	return eps
End


//	here is the way to get rotation axis and angle.
//	From rotation matrices A, B, you get C=A*B^(-1)
//	The angle = acos((C11+C22+C33-1)/2)
			// cos(angle) = { Trace(A*Binverse)-1 } / 2
//	axis direction (x,y,z): x=C23-C32, y=C31-C13, z=C12-C21
//	Make sure to normalize A and B first.
//		Wenge

Function AngleBetweenMats(am,bm)
	Wave am,bm										// the two matricies

	Make/N=(3,3)/O am__, bm__
	Make/N=3/O oneaxis__

	oneaxis__ = am[p][0]
	normalize(oneaxis__)
	am__[][0] = oneaxis__[p]
	oneaxis__ = am[p][1]
	normalize(oneaxis__)
	am__[][1] = oneaxis__[p]
	oneaxis__ = am[p][2]
	normalize(oneaxis__)
	am__[][2] = oneaxis__[p]

	oneaxis__ = bm[p][0]
	normalize(oneaxis__)
	bm__[][0] = oneaxis__[p]
	oneaxis__ = bm[p][1]
	normalize(oneaxis__)
	bm__[][1] = oneaxis__[p]
	oneaxis__ = bm[p][2]
	normalize(oneaxis__)
	bm__[][2] = oneaxis__[p]

	MatrixTranspose bm__
	MatrixMultiply am__,bm__		// M_product now contain product
	KillWaves/Z am__, bm__, oneaxis__
	return rotationAngleOfMat(M_product)
//	return acos((MatrixTrace(M_product)-1)/2)*180/PI
End


Function Make_rotationAngles()
	// make total rotation angles from rmat
	Variable N=DimSize(rmat,2)
	Make/N=(N)/O rotAngles
	Make/N=(3,3)/O mat__
	Wave rmat=rmat
	Variable i=0
	do
		mat__ = rmat[p][q][i]
		rotAngles[i] = rotationAngleOfMat(mat__)
		i += 1
	while (i<N)
	Killwaves/Z mat__
//	WaveStats/Q rotAngles
//	printf "range of rotation angles is [%g�,%g�],  with %d bad points\r",V_min,V_max,V_numNaNs
End
Proc Make_rotationAnglesOLD()
	// make total rotation angles from rmat
	Silent 1 ; PauseUpdate
	Variable N=DimSize(rmat,2)
	Make/N=(N)/O rotAngles
	Make/N=(3,3)/O mat__
	Variable i=0
	do
		mat__ = rmat[p][q][i]
		rotAngles[i] = rotationAngleOfMat(mat__)
		i += 1
	while (i<N)
	Killwaves/Z mat__
	WaveStats/Q rotAngles
	printf "range of rotation angles is [%g�,%g�],  with %d bad points\r",V_min,V_max,V_numNaNs
EndMacro
Function MakeRotationAnglesSigned(h,k,l)
	// makes signed grain to grain rotation angles, sign is referenced to hkl direction of orientations[][][0].
	// It does not show only the hkl component of the rotation
	Variable h,k,l							// reference direction for axis (for determinin sign)
	if (numtype(h) || numtype(k) || numtype(l))	// invalid inputs, so prompt user
		h = numtype(h) ? 1 : h
		k = numtype(k) ? 1 : k
		l = numtype(l) ? 1 : l
		Prompt h, "h"
		Prompt k, "k"
		Prompt l, "l"
		DoPrompt "for grain to grain angle, give reference direction for sign", h,k,l
		if (V_flag)
			return 1
		endif
	endif

	Wave orientations=orientations
	Wave rmat=rmat
	Variable N=DimSize(rmat,2)
	Make/N=(N)/O rotAngles
	Make/N=(3,3)/O mat__
	Make/N=3/O temp_axis__, ref_axis__
	mat__ = orientations[p][q][0]
	ref_axis__[0] = h
	ref_axis__[1] = k
	ref_axis__[2] = l
	MatrixMultiply mat__, ref_axis__
	Wave M_product=M_product
	ref_axis__ = M_product					// this is now the reference axis (in hkl direction)
	Variable i=0
	do
		mat__ = rmat[p][q][i]
		rotAngles[i] = axisOfMatrix(mat__,temp_axis__)	// compute total rotation angle, and axis
		rotAngles[i] = rotationAngleOfMat(mat__)			// recompute total rotation angle
		if (MatrixDot(ref_axis__,temp_axis__)<0)		// add sign
			rotAngles[i] *= -1
		endif
		i += 1
	while (i<N)
	Killwaves/Z mat__, temp_axis__, ref_axis__, M_product
	WaveStats/Q rotAngles
	printf "range of rotation angles is [%g�,%g�],  with %d bad points\r",V_min,V_max,V_numNaNs
End
//Function rotationAngleOfMat(rot)				// returns the total rotation angle of a matrix 'rot'
//	Wave rot									// the rotation matrix
//	Variable trace = MatrixTrace(rot)			// trace = 1 + 2*cos(theta)
//	Variable cosine = (trace-1)/2				// cosine of the rotation angle
//	cosine = (cosine>1) ? (2-cosine) : cosine
//	return acos(cosine)*180/PI				// rotation angle in degrees
//End
Function rotationCosineOfMat(rot)				// returns cos(total rotation angle of matrix) for matrix 'rot'
	Wave rot
	Variable trace = MatrixTrace(rot)			// trace = 1 + 2*cos(theta)
	return (trace-1)/2						// cosine of the rotation angle
End

// use version in Utility_JZT.ipf
//Function axisOfMatrix(rot,axis)				// compute angle and axis of a rotation matrix
//	// returns total rotation angle, and sets axis to the axis of the total rotation
//	Wave rot									// rotation matrix
//	Wave axis									// axis of the rotation (angle is returned)
//
//DoAlert 0,"check definition of routine 'rotationMatAboutAxis()'"
//	Make/N=(3,3)/O axisMat__
//	axis = rot[p][0]
//	normalize(axis)
//	axisMat__[][0] = axis[p]
//
//	axis = rot[p][1]
//	normalize(axis)
//	axisMat__[][1] = axis[p]
//
//	axis = rot[p][2]
//	normalize(axis)
//	axisMat__[][2] = axis[p]
//
//	axis[0] = axisMat__[1][2] - axisMat__[2][1]
//	axis[1] = axisMat__[2][0] - axisMat__[0][2]
//	axis[2] = axisMat__[0][1] - axisMat__[1][0]
//	normalize(axis)
//	Variable trace = MatrixTrace(axisMat__)			// trace = 1 + 2*cos(theta)
//	KillWaves /Z axisMat__
//	return acos((trace-1)/2)*180/PI					// rotation angle in degrees
//End


//Function rotationMatAboutAxis(axis,angle,mat)
//	Wave axis				// axis about which to rotate
//	Variable angle			// angle to rotate (passed as degrees, but convert to radians)
//	Wave mat				// desired rotation matrix
//	angle *= PI/180
//
//	Make/N=3/O/D xhat_rotate__, yhat_rotate__, zhat_rotate__
//	Wave xhat=xhat_rotate__
//	Wave zhat=zhat_rotate__
//
//	zhat = axis
//	if (normalize(zhat) <=0)
//		KillWaves/Z xhat_rotate__, W_Cross, zhat_rotate__
//		return 1			// error, cannot rotate about a zero length axis
//	endif
//
//	Variable i
//	i = ( abs(zhat[0])<= abs(zhat[1]) ) ? 0 : 1
//	i = ( abs(zhat[2])< abs(zhat[i]) ) ? 2 : i
//	xhat = zhat
//	xhat[i] = 2							// choose x-z plane
//	normalize(xhat)
//	Variable dotxz = MatrixDot(xhat,zhat)
//	xhat -= dotxz*zhat					// xhat is now perpendicular to zhat
//	normalize(xhat)					// xhat is now normalized vector perp to zhat
//	Cross zhat,xhat						//	was: cross(zhat,xhat,yhat)
//	Wave yhat=W_Cross
//	normalize(yhat)					// yhat is now normalized vector perp to zhat and xhat
////	xhat, yhat, zhat is now an orthonormal system with zhat || to axis
//
////	now rotate around zhat in the xhat -> yhat direction, an angle of angle
//	Variable cosa = cos(angle)
//	Variable sina = sin(angle)
//	Redimension/N=(3,3) mat
//	for(i=0;i<3;i+=1)
//		mat[][i] = (xhat[i]*cosa+yhat[i]*sina)*xhat[p] + (-xhat[i]*sina + yhat[i]*cosa)*yhat[p] + zhat[i]*zhat[p]
//	endfor
//	KillWaves/Z xhat_rotate__, W_Cross, zhat_rotate__
//End


Function rotation(mat,axis,angle)		// rotate mat about x, y, or z by angle
	Wave mat				// matrix to rotate
	Variable axis			// axis to rotate about 0=x, 1=y, 2=z
	Variable angle			// angle to rotate (radians)

	Variable cosa = cos(angle)
	Variable sina = sin(angle)
	Make/N=(3,3)/O mat_rot_
	mat_rot_ = (p==q)

	switch(axis)
		case 0:			// rotate about x-axis
			mat_rot_[1][1] = cosa
			mat_rot_[2][2] = cosa
			mat_rot_[1][2] = -sina
			mat_rot_[2][1] =   sina
			break
		case 1:			// rotate about y-axis
			mat_rot_[0][0] = cosa
			mat_rot_[2][2] = cosa
			mat_rot_[0][2] =   sina
			mat_rot_[2][0] = -sina
			break
		case 2:			// rotate about z-axis
			mat_rot_[0][0] = cosa
			mat_rot_[1][1] = cosa
			mat_rot_[0][1] = -sina
			mat_rot_[1][0] =   sina
			break
		default:
	endswitch
	Killwaves/Z mat_rot_
End


// ============================================================================
// ========  this section is for displaying stereo graphic projections ============================


Function DisplayPoleFigure(h0,k0,l0,Sorients)
	Variable h0,k0,l0
	String Sorients

	if (numtype(h0) || numtype(k0) || numtype(l0) || exists(Sorients)!=1)
		h0 = (numtype(h0)) ? 0 : h0
		k0 = (numtype(k0)) ? 0 : k0
		l0 = (numtype(l0)) ? 1 : l0
		if (strlen(Sorients)<1)
			Sorients = "orientations"
		endif
		Prompt h0, "h"
		Prompt k0, "k"
		Prompt l0, "l"
		Prompt Sorients, "orientations matrix name", popup, WaveList("*orient*",";","")
		DoPrompt "hkl and range", h0,k0,l0,Sorients
		if (V_Flag)
			return 1
		endif
		printf "     DisplayPoleFigure(%d,%d,%d,\"%s\")\r" ,h0,k0,l0,Sorients
	endif
	if (!exists(Sorients))
		Abort "orientations matrix '"+Sorients+"' does not exist"
	endif
	Wave orients=$Sorients

	h0 = round(h0)
	k0 = round(k0)
	l0 = round(l0)
	Variable i, factor = gcf3(h0,k0,l0)					// simplify hkl
	h0 = abs(h0)/factor
	k0 = abs(k0)/factor
	l0 = abs(l0)/factor
	String namstr
	sprintf namstr, "%d%d%d_", abs(h0),abs(k0),abs(l0)

	Variable nf				// number of points in a family
	nf = MakeCubicFamilyOf_hkl(h0,k0,l0)
	String xList="", yList=""
	String xName="", yName=""
	i=0
	do
		xList += "Xstereo"+namstr+num2istr(i)+";"	// build the list of wave names
		yList += "Ystereo"+namstr+num2istr(i)+";"
		i += 1
	while (i<nf)

	i=0
	do
		xName = StringFromList(i,xList)			// null out all of the x & y waves
		yName = StringFromList(i,yList)			// create them if necessary
		if (exists(xName)==1)
			Redimension/N=0 $xName
		else
			Make/N=0 $xName
		endif
		if (exists(yName)==1)
			Redimension/N=0 $yName
		else
			Make/N=0 $yName
		endif
		i += 1
	while (i<nf)

	Make/N=(3,3)/O orient_temp__
	i=0
	do													// loop over each grain
		orient_temp__ = orients[p][q][i]
		AppendFamilyOfhkl2xyStereo(hh,kk,ll,orient_temp__,xList,yList)
		i += 1
	while (i<DimSize(orients,2))
	KillWaves/Z orient_temp__

	ResumeUpdate
	PauseUpdate
	Display /W=(8,56,372,406)
	i=0
	do
		xName = StringFromList(i,xList)			// null out all of the x & y waves
		yName = StringFromList(i,yList)			// create them if necessary
		WaveStats/Q $yName
		if (V_npnts>0)
			AppendToGraph $yName vs $xName
		endif
		i += 1
	while (i<nf)

	String str
	Wave hh=hh
	Wave kk=kk
	Wave ll=ll
	sprintf str, "pole figure\r(%d%d%d)", abs(hh[0]),abs(kk[0]),abs(ll[0])
	TextBox/N=text0/F=0/S=3/A=LT/X=-13/Y=-3.5 str
	TextBox/N=text1/F=0/S=3/A=LB/X=-17/Y=-9 "\\Z10"+note(orients)+"n.IND"
	ResumeUpdate
	Execute "Graph_streo_Style()"

	Variable/G polarIncrement=10
	Variable/G phiIncrement=30
	Variable /G deltaAngle
	Make/N=3/O axis_of_rotation, hkl_of_axis
	String fldr= GetDataFolder(1)

	SetVariable polarControl,pos={280,3},size={80,17},proc=StereoSetVarProc,title="polar"
	SetVariable polarControl,font="Helvetica",fSize=14
	SetVariable polarControl,limits={0,90,1},value= $fldr+"polarIncrement"
	SetVariable phiControl,pos={281,23},size={80,17},proc=StereoSetVarProc,title="phi"
	SetVariable phiControl,font="Helvetica",fSize=14
	SetVariable phiControl,limits={0,360,1},value= $fldr+"phiIncrement"
	Button FullSize0,pos={311,41},size={50,20},proc=FullSizeButtonProc,title="Full"
	Button Clear0,pos={311,65},size={50,20},proc=ClearButtonProc,title="Clear"
	Button Axis0,pos={311,89},size={50,20},proc=ShowAxisButtonProc,title="�axis"
	ValDisplay delta_angle,pos={202,5},size={72,17},title="ơ",fSize=12
	ValDisplay delta_angle,format="%.3f",limits={0,0,0},barmisc={0,1000}
	String command = "ValDisplay delta_angle,value= #"+fldr+"deltaAngle"
	Execute command
	if (NumberByKey("IGORVERS",igorinfo(0)) > 4.04)
		SetWindow kwTopWin, hook=resizeHookFunc , hookcursor=0, hookEvents=0
	else
		SetWindow kwTopWin, hook=resizeHookFunc , hookcursor=0, hookEvents=1
	endif
End

Function resizeHookFunc(infoStr)
	String infoStr
// print "hook called",infoStr
	String event = StringByKey("EVENT",infoStr)
	if (!cmpstr("modified",event) || !cmpstr("resize",event) || !cmpstr("mouseup",event))
		SquarizeGraph()
	endif
	return 0
End


Function CursorMovedHook(infoStr)
	String infoStr
//			 print "CursorMovedHook called  ",infoStr

	String wname=StringByKey("TNAME",infoStr)
	if (strsearch(wname,"Ystereo",0)!=0)				// wrong kind of wave
		return 0
	endif
	String csr=StringByKey("CURSOR",infoStr)
	Variable point=NumberByKey("POINT", infoStr)

	String fldr = GetWavesDataFolder(CsrWaveRef($csr),1)

	if (exists(fldr+"orientations")!=1)
		return 0
	endif
	Wave orientations=$(fldr+"orientations")

	strswitch(csr)
		case "A":
			Make/N=(3,3)/O $(fldr+"cursor_mat_a__")
			Wave mata = $(fldr+"cursor_mat_a__")
			mata = orientations[p][q][point]
			break
		case "B":
			Make/N=(3,3)/O $(fldr+"cursor_mat_b__")
			Wave matb = $(fldr+"cursor_mat_b__")
			matb = orientations[p][q][point]
			break
		default:
			return 0
	endswitch

	String dname = fldr+"deltaAngle"
	if (exists(dname)!=2)
		return 0
	endif
	String NmatA = fldr+"cursor_mat_a__"
	String NmatB = fldr+"cursor_mat_b__"
	String Naxis = fldr+"axis_of_rotation"
	String Nhkl = fldr+"hkl_of_axis"
	if (strlen(CsrWave(A))<1 || strlen(CsrWave(B))<1)
		NVAR deltaAngle=$dname
		deltaAngle = NaN
	elseif (exists(NmatA)==1 && exists(NmatB)==1 && exists(Naxis)==1)
		NVAR deltaAngle=$dname
		Wave matA=$NmatA
		Wave matB=$NmatB
		Wave axis=$Naxis
		Wave hkl=$Nhkl
//		print AngleBetweenMats(matA,matB)
		Make/N=(3,3)/O rot_temp__
		MakeRotationMatFromAB(matA,matB,rot_temp__)	// makes rotation mat between a and b
		deltaAngle = axisOfMatrix(rot_temp__,axis)
//		print "deltaAngle=",deltaAngle
		if (!WaveExists(hkl))
			Duplicate $Naxis $Nhkl
			Wave hkl=$Nhkl
		endif
		rot_temp__ = mata
		hklFromAxisOfRotation(rot_temp__,axis,hkl)
		hkl = round(24*hkl[p])							// picked 24 because it is 2^3 * 3^1
		Variable gcf = gcf3(hkl[0],hkl[1],hkl[2])
		hkl = hkl/gcf
		KillWaves/Z rot_temp__

		Variable xx,yy		// this section updates rotation axis marker on plot and creates it if necessary
		Variable mark		// used to set the marker
		mark = ComputeXYfromDirection(axis,xx,yy)
//		mark = mark ? 41 : 43
		Make/O axisDot={yy}
		SetScale/P x xx,1,"", axisDot
		if (strlen(WaveList("axisDot",";","WIN:"))<1)
			AppendToGraph axisDot
			ModifyGraph mode(axisDot)=3, msize(axisDot)=5, marker(axisDot)=1,mrkThick(axisDot)=1.5
			Tag/F=0/S=3/A=MB/X=0.00/Y=3/L=0 axisDot, 0
			AppendText "\\Z12\\{\"%g %g %g\",root:AlTJ2Ind:hkl_of_axis[0],root:AlTJ2Ind:hkl_of_axis[1],root:AlTJ2Ind:hkl_of_axis[2]}"
		endif
		ModifyGraph marker(axisDot)=mark

	endif
	return 0
End



Function DrawRadialLines(increment)
	Variable increment
	Variable lo, hi
	increment = (increment<=0) ? 1 : increment

	SquarizeGraph()
	Variable xlo,xhi,ylo,yhi
	GetAxis /Q bottom
	xlo = V_min
	xhi = V_max
	GetAxis /Q left
	ylo = V_min
	yhi = V_max

	Variable p1,p2,p3,p4			// phi of the four corners
	p1 = atan2(ylo,xlo)*180/PI
	p2 = atan2(ylo,xhi)*180/PI
	p3 = atan2(yhi,xlo)*180/PI
	p4 = atan2(yhi,xhi)*180/PI

	if (OriginExterior(xlo,xhi,ylo,yhi))
		lo = min(p1,p2)
		lo = min(lo,p3)
		lo = min(lo,p4)
		hi = max(p1,p2)
		hi = max(hi,p3)
		hi = max(hi,p4)
		if ((hi-lo)>180)
			Variable swap
			swap = lo
			lo = hi
			hi = swap+360
		endif
		lo = ceil(lo/increment)*increment
	else
		lo= 0
		hi = 360
	endif

	SetDrawLayer UserFront
	Variable xx,yy,phi=lo
	do
		xx = cos(phi*PI/180)
		yy = sin(phi*PI/180)
		SetDrawEnv xcoord= bottom,ycoord= left,dash= 1
		DrawLine 0,0,xx,yy
		phi += increment
	while (phi<hi)
	return 0
End


Function DrawPhiCircles(increment)
	Variable increment
	increment = (increment<=0) ? 1 : increment

	SquarizeGraph()
	Variable xlo,xhi,ylo,yhi
	GetAxis /Q bottom
	xlo = V_min
	xhi = V_max
	GetAxis /Q left
	ylo = V_min
	yhi = V_max

	Variable c1,c2,c3,c4, lo,hi				// determine range of radii for circles
	c1 = sqrt(xlo^2 + ylo^2)
	c2 = sqrt(xhi^2 + ylo^2)
	c3 = sqrt(xlo^2 + yhi^2)
	c4 = sqrt(xhi^2 + yhi^2)
	hi = max(c1,c2)
	hi = max(hi,c3)
	hi = max(hi,c4)
	hi = min(1-1e-5,hi)
	hi = 360*atan(hi)/PI
	if (OriginExterior(xlo,xhi,ylo,yhi))
		lo = min(c1,c2)
		lo = min(lo,c3)
		lo = min(lo,c4)
		lo = max(0,lo)
		lo = 360*atan(lo)/PI
		lo = ceil(lo/increment)*increment
	else
		lo = 0
	endif

	SetDrawLayer /K UserFront
	SetDrawEnv xcoord= bottom,ycoord= left,linepat= 3,fillpat= 0
	DrawOval -1,1,1,-1
	Variable xx,i=lo+increment
	do
		xx = tan(i/2*PI/180)
		SetDrawEnv xcoord= bottom,ycoord= left,linepat= 3,fillpat= 0
		DrawOval -xx,xx,xx,-xx
		i += increment
	while (i<hi)
	return 0
End


Function OriginExterior(xlo,xhi,ylo,yhi)
	Variable xlo,xhi,ylo,yhi
	if ((xlo<0) && (0<xhi)  &&  (ylo<0) && (0<yhi))
		return 0
	endif
	return 1
End


Function SquarizeGraph()
	Variable xlo,xhi,ylo,yhi
	GetAxis /Q bottom
	xlo = V_min
	xhi = V_max
	GetAxis /Q left
	ylo = V_min
	yhi = V_max
	Variable dax,cen
	if ((yhi-ylo)<(xhi-xlo))						// square up the plot
		dax=(xhi-xlo)/2
		cen = (yhi+ylo)/2
		SetAxis left cen-dax,cen+dax
		GetAxis /Q left
		ylo = V_min
		yhi = V_max
	elseif ((xhi-xlo)<(yhi-ylo))
		dax=(yhi-ylo)/2
		cen = (xhi+xlo)/2
		SetAxis bottom cen-dax,cen+dax
		GetAxis /Q bottom
		xlo = V_min
		xhi = V_max
	endif
End


Function StereoSetVarProc(ctrlName,varNum,varStr,varName) : SetVariableControl
	String ctrlName
	Variable varNum
	String varStr
	String varName
	NVAR polarIncrement=polarIncrement
	NVAR phiIncrement=phiIncrement
	DrawPhiCircles(polarIncrement)
	DrawRadialLines(phiIncrement)
End


Function FullSizeButtonProc(ctrlName) : ButtonControl
	String ctrlName
//	ModifyGraph/Z axOffset(left)=-4.8,axOffset(bottom)=-1.77778
//	ModifyGraph/Z axThick=0
	SetAxis/Z left -1.0,1.0
	SetAxis/Z bottom -1.0,1.0
	ShowInfo
	SetDrawLayer UserFront
	SetDrawEnv xcoord= bottom,ycoord= left,linepat= 3,fillpat= 0
	DrawOval -1,1,1,-1
//	Execute "Graph_streo_Style()"
End


Function ClearButtonProc(ctrlName) : ButtonControl
	String ctrlName
	SetDrawLayer /K UserFront
	SetDrawEnv xcoord= bottom,ycoord= left,linepat= 3,fillpat= 0
	DrawOval -1,1,1,-1
End


Function ShowAxisButtonProc(ctrlName) : ButtonControl
	String ctrlName
	String aname=GetWavesDataFolder(WaveRefIndexed("",0,1),1)+"axis_of_rotation"
	String hname=GetWavesDataFolder(WaveRefIndexed("",0,1),1)+"hkl_of_axis"
	String adotname=GetWavesDataFolder(WaveRefIndexed("",0,1),1)+"axisDot"
	Wave adot=$adotname

	//	ShowAxisTable()
	if (strlen(WaveList("axisDot",";","WIN:")))
		RemoveFromGraph axisDot
	elseif (WaveExists(adot))
		AppendToGraph $adotname
		ModifyGraph mode(axisDot)=3, msize(axisDot)=5, marker(axisDot)=1,mrkThick(axisDot)=1.5
		Tag/F=0/S=3/A=MB/X=0.00/Y=3/L=0 axisDot, 0
		AppendText "\\Z12\\{\"%g %g %g\","+hname+"[0],"+hname+"[1],"+hname+"[2]}"
	endif
End


Function ShowAxisTable()
	String aname=GetWavesDataFolder(WaveRefIndexed("",0,1),1)+"axis_of_rotation"
	String hname=GetWavesDataFolder(WaveRefIndexed("",0,1),1)+"hkl_of_axis"
	if (exists(aname)==1 && exists(hname)==1)
		Edit/W=(381,56,610,174) $aname,$hname
		Execute "ModifyTable alignment(Point)=1,width(Point)=42,width(axis_of_rotation)=90,alignment(hkl_of_axis)=1"
		Execute "ModifyTable width(hkl_of_axis)=68"
	elseif (exists(aname)==1)
		Edit/W=(349,57,541,180) $aname
		Execute "ModifyTable width(Point)=42"
	endif
End
// ============================================================================
// ========  this section is for reading matricies from the files ================================


Function ReadInManyMatricies(i0,i1,base)			// read in a group of matricies from a file, and store all of them in 'orientations'
	Variable i0,i1
	String base

	if (strlen(base)<1)									// get base name of file to use
		Variable fileRef
		Open /D/R/T="????" fileRef
		base = S_fileName[0,strsearch(S_fileName,".",0)-1]
		Prompt base, "base file name (no extension)"
		Prompt i0, "index of first file"
		Prompt i1, "index of last file"
		DoPrompt "file and range", base,i0,i1
		if (V_Flag || (i1<i0))
			return 1
		endif
		print "base file name = ", base
	endif
	String ext=".IND"
	Make/N=(3,3,i1-i0+1)/O orientations
	Make/N=3/O astar,bstar,cstar
	Note /K orientations
	Note orientations, "base="+base
	orientations = NaN
	Variable i=i0
	do													// loop over each grain
		if (readOneMatrix(base+num2istr(i)+ext))
			orientations[0,2][0][i-i0] = astar[p]
			orientations[0,2][1][i-i0] = bstar[p]
			orientations[0,2][2][i-i0] = cstar[p]
		endif
		i += 1
	while (i<=i1)
	KillWaves/Z astar,bstar,cstar

	if (exists("orientationsOrig")==1)
		DoAlert 1, "Overwite backup orientations (orientationsOrig)?"
		if (V_flag!=1)
			return 0
		endif
		print "overwrote orientationsOrig with new values"
	endif
	Duplicate/O orientations orientationsOrig
End


Function readOneMatrix(fullFileName)
	String fullFileName

	Variable fileRef
	Open /R/Z fileRef  fullFileName
	if (V_Flag)
		return 0
	endif

	String line
	do
		FReadLine /N=255 fileRef, line				// find the beginning of the matrix
		if (strlen(line)>250)						// this check is beacuse of the silly files I sometimes get
			printf "\rfor file:  %s\rthe lines in this file are too line!! What is it?", fullFileName
		Abort "the lines in this file are too ling!! What is it?"
		endif
	while (strsearch(line,"matrix h  k  l -> X  Y  Z",0)<0)

	// next read will get the a*, b*, c*
	Wave astar=astar
	Wave bstar=bstar
	Wave cstar=cstar
	Variable aa,bb,cc
	FReadLine /N=255 fileRef, line
	sscanf line, "%g %g %g", aa,bb,cc
	if (V_flag!=3)									// this check is beacuse of the silly files I sometimes get
		printf "\rfor file:  %s\rthe lines in this file are too line!! What is it?", fullFileName
		Abort "the lines in this file are probably too ling!! What is it?"
	endif
	astar[0] = aa ; bstar[0] = bb ; cstar[0] = cc

	FReadLine /N=255 fileRef, line
	sscanf line, "%g %g %g", aa,bb,cc
	if (V_flag!=3)									// this check is beacuse of the silly files I sometimes get
		printf "\rfor file:  %s\rthe lines in this file are too line!! What is it?", fullFileName
		Abort "the lines in this file are probably too ling!! What is it?"
	endif
	astar[1] = aa ; bstar[1] = bb ; cstar[1] = cc

	FReadLine /N=255 fileRef, line
	sscanf line, "%g %g %g", aa,bb,cc
	if (V_flag!=3)									// this check is beacuse of the silly files I sometimes get
		printf "\rfor file:  %s\rthe lines in this file are too line!! What is it?", fullFileName
		Abort "the lines in this file are probably too ling!! What is it?"
	endif
	astar[2] = aa ; bstar[2] = bb ; cstar[2] = cc

	Close fileRef
	return 1
End


// ============================================================================
// ========  this section is for going from matrix to hkl =====================================


Function AppendFamilyOfhkl2xyStereo(hh,kk,ll,orient,xList,yList)
	Wave hh,kk,ll							// family of hkl
	Wave orient								// orientation for this crystallite
	String xList,yList						// list of wave names to recieve the x,y coords

	Variable nf=numpnts(hh)				// number of points in a family
	Variable i=0
	if (nf != ItemsInList(xList))
		Abort "AppendFamilyOfhkl2xyStereo, number of hkl not equal to number in list"
	endif

	// extend the x & y waves by one
	String xName, yName
	xName = StringFromList(nf-1,xList)
	yName = StringFromList(nf-1,yList)
	if (exists(xName)!=1 || exists(yName)!=1)		// if they do not exist, exit
		Abort "AppendFamilyOfhkl2xyStereo, waves to recieve x,y must already exist"
	endif

	Variable g0,g1,g2, ge
	Variable delta											// polar angle (=90-lattitude)
	Variable rad,phi										// coords of stereographic projection in polar coords
	Variable xx,yy											// x,y version of rad, phi
	Variable N=numpnts($xName)
	Variable/C lz
	i = 0
	do														// loop over the family
		// for each hkl create plotable x-y coord for stereographic plot
		g0 = hh[i]*orient[0][0] + kk[i]*orient[0][1] + ll[i]*orient[0][2]	// the point in reciprocal space
		g1 = hh[i]*orient[1][0] + kk[i]*orient[1][1] + ll[i]*orient[1][2]
		g2 = hh[i]*orient[2][0] + kk[i]*orient[2][1] + ll[i]*orient[2][2]
		if (g2>=0)
			ge = sqrt(g0^2+g1^2)							// length of g in equatorial plane
			delta = atan2(ge,g2)							// polar angle of recip lattice vector
			rad = tan(delta/2)								// radius of stereographic projection
			phi = atan2(g1,g0)
			xx = rad*cos(phi)
			yy = rad*sin(phi)
		else
			xx = NaN
			yy = NaN
		endif

		xName = StringFromList(i,xList)
		yName = StringFromList(i,yList)
		InsertPoints N, 1, $xName, $yName
		Wave xw=$xName
		Wave yw = $yName
		xw[N] = xx
		yw[N] = yy	
		N += 1
		i += 1
	while(i<nf)
	return numpnts(xw)
End


Function ComputeXYfromDirection(axis,xx,yy)
	Wave axis								// direction vector in reciprocal space
	Variable &xx, &yy						// x,y version of rad, phi
											// returns 0 for normal, 1 for reversed
	Variable g0,g1,g2, ge
	Variable delta									// polar angle (=90-lattitude)
	Variable rad,phi								// coords of stereographic projection in polar coords
	// create plotable x-y coord for stereographic plot
	g0 = axis[0]									// the point in reciprocal space
	g1 = axis[1]
	g2 = axis[2]
	if (g2==0)
		xx = NaN
		yy = NaN
		return 0
	endif
	if (g2<0)
		g0  *= -1
		g1  *= -1
		g2  *= -1
	endif
	ge = sqrt(g0^2+g1^2)							// length of g in equatorial plane
	delta = atan2(ge,g2)							// polar angle of recip lattice vector
	rad = tan(delta/2)								// radius of stereographic projection
	phi = atan2(g1,g0)
	xx = rad*cos(phi)
	yy = rad*sin(phi)
	return 	(axis[2]<0)
End


Function hklFromAxisOfRotation(orient,G,hkl)		// in reciprocal space, convert (x*,y*,z*) to (hkl)
	Wave G											// vector in reciprocal space (x,y,z)
	Wave orient										// orientation matrix
	Wave hkl										// desired hkl
	//	G = orient*hkl
	//	orientInv * G = orientInv * orient*hkl
	//	hkl = orientInv * G

	Duplicate/O $NameOfWave(orient) orient_Inv__, temp__
	orient_Inv__ = (p==q)
	MatrixLLS /O/Z  temp__ orient_Inv__			// orient_Inv__ is now orient^-1
	MatrixMultiply orient_Inv__, G				// hkl is now in M_product
	Wave M_product=M_product
	hkl = M_product
	KillWaves /Z orient_Inv__, temp__
End


Function MakeCubicFamilyOf_hkl(h1,k1,l1)
	Variable h1,k1,l1						// the root hkl
	//	h1=1  ;  k1=2  ;  l1=3

	Make/N=48/O hh,kk,ll
	Variable i,n=0
	hh=nan
	kk=nan
	ll=nan

	Variable factor = gcf3(h1,k1,l1)	
	h1 /= factor
	k1 /= factor
	l1 /= factor

	i = 0
	do
		hh[n] = h1 * sign(.5-(1&i))
		kk[n] = k1 * sign(.5-(2&i))
		ll[n] = l1 * sign(.5-(4&i))
		n += 1
		i += 1
 	while(i<8)
	i = 0
	do
		hh[n] = k1 * sign(.5-(1&i))
		kk[n] = h1 * sign(.5-(2&i))
		ll[n] = l1 * sign(.5-(4&i))
		n += 1
		i += 1
 	while(i<8)

	i = 0
	do
		hh[n] = l1 * sign(.5-(1&i))
		kk[n] = h1 * sign(.5-(2&i))
		ll[n] = k1 * sign(.5-(4&i))
		n += 1
		i += 1
 	while(i<8)
	i = 0
	do
		hh[n] = l1 * sign(.5-(1&i))
		kk[n] = k1 * sign(.5-(2&i))
		ll[n] = h1 * sign(.5-(4&i))
		n += 1
		i += 1
 	while(i<8)

	i = 0
	do
		hh[n] = h1 * sign(.5-(1&i))
		kk[n] = l1 * sign(.5-(2&i))
		ll[n] = k1 * sign(.5-(4&i))
		n += 1
		i += 1
 	while(i<8)
	i = 0
	do
		hh[n] = k1 * sign(.5-(1&i))
		kk[n] = l1 * sign(.5-(2&i))
		ll[n] = h1 * sign(.5-(4&i))
		n += 1
		i += 1
 	while(i<8)

	RemoveZeroTriplets(hh,kk,ll)
	RemoveDuplicateDirections(hh,kk,ll)
	return numpnts(hh)
End


Function RemoveDuplicateDirections(aa,bb,cc)	// remove duplicate directions in h,k,l triplets
	Wave aa,bb,cc
	if ((numpnts(aa) != numpnts(bb)) || (numpnts(aa) != numpnts(cc)))
		Abort "in RemoveDuplicateDirections, waves must have same number of points"
	endif
	Variable N=numpnts(aa)
	Variable a0,b0,c0					// direction cosines
	Variable ax,bx,cx
	Variable norm
	Variable tol=1e-5
	Variable i0=0						// point to test against
	Variable i							// possible duplicate to check
	do
		norm = sqrt(aa[i0]^2+bb[i0]^2+cc[i0]^2)
		if (norm>0)
			a0 = aa[i0]/norm
			b0 = bb[i0]/norm
			c0 = cc[i0]/norm
		else
			a0=0; b0=0; c0=0
		endif
		i = i0+1
		do
			norm = sqrt(aa[i]^2+bb[i]^2+cc[i]^2)
			if (norm>0)
				ax = aa[i]/norm
				bx = bb[i]/norm
				cx = cc[i]/norm
			else
				ax=0; bx=0; cx=0
			endif
			if ((abs(ax-a0)+abs(bx-b0)+abs(cx-c0)) < tol)
				DeletePoints i,1,aa,bb,cc
			else
				i += 1
			endif
		while (i<numpnts(aa))

		i0 += 1
	while (i0<(numpnts(aa)-1))
	return numpnts(aa)
End


Function RemoveZeroTriplets(aa,bb,cc)
	Wave aa,bb, cc
	if ((numpnts(aa) != numpnts(bb)) || (numpnts(aa) != numpnts(cc)))
		Abort "in RemoveZeroTriplets, waves must have same number of points"
	endif
	Variable i=0
	do							// remove points with NaN
		if ( (aa[i]^2 + bb[i]^2 + cc[i]^2) <= 0 )
			DeletePoints i,1,aa,bb, cc
		else
			i += 1
		endif
	while (i<numpnts(aa))
	return numpnts(aa)
End


Function gcf3(a,b,c)			// find greatest common factor of a,b,c (all integers)
	Variable a,b,c
	a = round(abs(a))
	b = round(abs(b))
	c = round(abs(c))

	Variable fmax				// biggest possible factor to search for
	Variable a0,b0,c0
	fmax = max(a,b)
	fmax = max(fmax,c)
	a0 = a ? a : fmax
	b0 = b ? b : fmax
	c0 = c ? c : fmax
	fmax = min(a0,b0)			// gcf cannot be > smallest of the three (not including zero)
	fmax = min(c0,fmax)
	fmax = floor(fmax)

	Variable factor=1,i=2
	do
		if ((mod(a,i)+mod(b,i)+mod(c,i))<=0)
			factor = i
		endif
		i += 1
	while (i<=fmax)

	return factor
End


Function epsilon_ijk(i,j,k)
	Variable i,j,k					// returns +1 for even permutaion, -1 for odd permutaions

	Variable ps = sign(i)*sign(j)*sign(k)
	i = round(abs(i))
	j = round(abs(j))
	k = round(abs(k))

	if (i*j*k ==0)					// no zeros allowed
		return 0
	endif
	if (i>3 || j>3 || k>3)			// must be in range [1,3]
		return 0
	endif
	if (i==j || j==k || k==i)		// do duplicates
		return 0
	endif
	Variable perm = (i==1) + 2*(j==2)  + 4*(k==3)
	perm = ((perm==0) || (perm==7)) ? 1 : -1
	return ps*perm
End
//Function test48()
//	Variable i,j,k,eps,n=1
//	Variable pos=0,neg=0
//	for(k=-3;k<4;k+=1)
//		for(j=-3;j<4;j+=1)
//			for(i=-3;i<4;i+=1)
//				eps = epsilon_ijk(i,j,k)
//				if (eps)
//					printf "%3d    %+d, %+d, %+d     %+d\r",n,i,j,k,eps
//					pos = (eps>0) ? pos+1 : pos
//					neg = (eps>0) ? neg+1 : neg
//					n += 1
//				endif
//			endfor
//		endfor
//	endfor
//	print "pos = ",pos,"     neg=",neg
//End


// ============================================================================
// ========  display and style macros ==================================================



Proc Graph_streo_Style() : GraphStyle
	PauseUpdate; Silent 1		// modifying window...
	ModifyGraph/Z gfSize=18,width={Aspect,1}
	ModifyGraph/Z mode=3
	ModifyGraph/Z marker=19
	ModifyGraph/Z rgb[0]=(0,0,65280),rgb[2]=(2,39321,1),rgb[3]=(65535,16385,55749),rgb[4]=(0,0,0)
	ModifyGraph/Z rgb[5]=(52428,34958,1), rgb[6]=(1,52428,52428), rgb[7]=(32768,32770,65535)
	ModifyGraph/Z msize=2
	ModifyGraph/Z tick=3
	ModifyGraph/Z noLabel=2
	ModifyGraph/Z lowTrip=0.001
	ModifyGraph/Z standoff=0
	ModifyGraph/Z axOffset(left)=-4.8,axOffset(bottom)=-1.77778
	ModifyGraph/Z axThick=0
	SetAxis/Z left -1.0,1.0
	SetAxis/Z bottom -1.0,1.0
	ShowInfo
	SetDrawLayer UserFront
	SetDrawEnv xcoord= bottom,ycoord= left,linepat= 3,fillpat= 0
	DrawOval -1,1,1,-1
EndMacro

Window Table_hkl() : Table
	PauseUpdate; Silent 1		// building window...
	String fldrSav= GetDataFolder(1)
	SetDataFolder root:test:
	Edit/W=(5,44,197,343) hh,kk,ll
	ModifyTable width=30
	SetDataFolder fldrSav
EndMacro


Proc Graph_rotationAnglesStyle() : GraphStyle
	PauseUpdate; Silent 1		// modifying window...
	ModifyGraph/Z gfSize=18
	ModifyGraph/Z mode=3
	ModifyGraph/Z marker=19
	ModifyGraph/Z tick=2
	ModifyGraph/Z mirror=1
	ModifyGraph/Z minor=1
	ModifyGraph/Z lowTrip=0.001
	ModifyGraph/Z lblMargin(left)=3
	ModifyGraph/Z standoff=0
	ModifyGraph/Z axOffset(left)=-2.4,axOffset(bottom)=-0.777778
	ModifyGraph/Z lblLatPos(left)=-1
	ModifyGraph zero(left)=2
	Label/Z left "rotation angle�"
	Label/Z bottom "grain index"
	ShowInfo
	Cursor/P A $StringFromList(0, TraceNameList("",";",1)) 0
EndMacro



// ======================  END ==================================================
// ============================================================================

