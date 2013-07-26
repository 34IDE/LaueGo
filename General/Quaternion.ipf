#pragma rtGlobals=1		// Use modern global access method.
#pragma ModuleName=Quaternion
#pragma version = 1.0

// For information about Quaternions, see:
//		http://en.wikipedia.org/wiki/Quaternions
//  and	http://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation


Function QuarternionRotate(quat,v)			// rotates vector v by quaternion quat, v is changed!
	Wave quat								// quaternion, 4-vector
	Wave v									// 3-vector

	Make/N=3/D/FREE vp					// rotated vector
	Make/N=3/D/FREE u=quat[p+1]		// vector part of quat
	Variable c2=quat[0]					// scalar part of quat = cos(angle/2)
 	Variable s2 = normalize(u)				// u becomes unit vector, s2 = sin(angle/2)
 	Variable cosa = c2*c2 - s2*s2			// cos(2a) = cos^2(a) = sin^2(a), use half-angle formulas
 	Variable sina = 2*s2*c2				// sin(2a) = 2 sin(a) cos(a)
	Variable dot = MatrixDot(u,v) 
	Cross u,v
	Wave C=W_Cross
	vp = (v-u*dot)*cosa + C*sina + u*dot	// (v - u(u.v))cosa + (u x v) sina + u(u.v)
	KillWaves/Z W_Cross
	v = vp									// do the replacement
End


Function QuaternionMultiply(A,B,C)		// compute  C = (A)(B)
	Wave A, B, C							// all waves are 1d waves with length=4
	Make/N=4/D/FREE Ctemp				// in case A or B is the same as C
	Ctemp[0] = A[0]*B[0] - A[1]*B[1] - A[2]*B[2] - A[3]*B[3]
	Ctemp[1] = A[0]*B[1] + A[1]*B[0] + A[2]*B[3] - A[3]*B[2]
	Ctemp[2] = A[0]*B[2] - A[1]*B[3] + A[2]*B[0] + A[3]*B[1]
	Ctemp[3] = A[0]*B[3] + A[1]*B[2] - A[2]*B[1] + A[3]*B[0]
	C = Ctemp								// now do the final assignment
	//
	// a is scalar, vector = (b,c,d),  so quaternion is {a,b,c,d}
	//
	//	(a1,b1,c1,d1) (a2,b2,c2,d1) = 
	//		  (	a1a2 - b1b2 - c1c2 - d1d2,
	//			a1b2 + b1a2 + c1d2 - d1c2,
	//			a1c2 - b1d2 + c1a2 + d1b2,
	//			a1d2 + b1c2 - c1b2 + d1a2  )
End


Function QuaternionNormalize(quat)		// Normalizes the quaternion, makes a pure rotation
	Wave quat
	Variable c2 = quat[0]					// scalar part = cos(angle/2), use this to define the angle
	Variable s2 = sin(acos(c2))
	Make/N=3/D/FREE u
	u = quat[p+1]
	normalize(u)
	u *= s2
End


Function RotationVector2Quaternion(rvec,quat)// convert rotation vector to a quaternion
	Wave rvec								// rotation vector length is angle radians (3-vector)
	Wave quat								// quaternion to receive the result (4-vector)

	Make/N=3/D/FREE u=rvec
	Variable angle=normalize(u)
	quat[0] = cos(angle/2)
	quat[1,3] = sin(angle/2) * u[p-1]
End


Function Quaternion2RotationVector(quat,rvec)// convert quaternion to a rotation vector
	Wave quat								// input quaternion to receive the result (4-vector)
	Wave rvec								// resultant rotation vector length is angle radians (3-vector)

	Variable angle = 2*acos(quat[0])		// w = cos(angle/2)
	rvec = quat[p+1]						// axis is last three terms in quaternion
	normalize(rvec)
	rvec *= angle							// length of rvec is rotation angle in radians
End


Function Quaternion2RotationMat(quat,rmat)// fill in the rotaion matrix from the quaternion
	Wave quat								// quaternion, the input, 4-vector
	Wave rmat								// 3x3 rotation matrix

	Variable a=quat[0], b=quat[1], c=quat[2], d=quat[3]
	Variable a2=a*a, b2=b*b, c2=c*c, d2=d*d
	rmat[0][0] = { a2 + b2 - c2 - d2, 2*b*c + 2*a*d, 2*b*d - 2*a*c }
	rmat[0][1] = { 2*b*c - 2*a*d, a2 - b2 + c2 - d2, 2*c*d + 2*a*b }
	rmat[0][2] = { 2*b*d + 2*a*c, 2*c*d - 2*a*b, a2 - b2 - c2 + d2 }
	//
	// rmat =
	//	a2+b2-c2-d2		2bc - 2ad			2bd + 2ac
	//	2bc + 2ad			a2-b2+c2-d2		2cd - 2ab
	//	2bd - 2ac			2cd + 2ab			a2-b2-c2+d2
End


Function RotationMat2Quaternion(rmat,quat)// fill in the quaternion from thee rotaion matrix, returns angle (degree)
	Wave rmat								// 3x3 rotation matrix, input
	Wave quat								// quaternion (4-vector), result

	Variable u=abs(rmat[0][0]), v=abs(rmat[1][1]), w=abs(rmat[2][2])
	Variable uu,uv,uw,vu,vv,vw,wu,wv,ww	// 9 matrix elements
	Variable id
	if (u>v && u>w)						// u is the biggest
		id = 0
		uu = rmat[0][0];	uv = rmat[0][1];	uw = rmat[0][2]
		vu = rmat[1][0];	vv = rmat[1][1];	vw = rmat[1][2]
		wu = rmat[2][0];	wv = rmat[2][1];	ww = rmat[2][2]
	elseif (v>u && v>w)					// v is the biggest
		id = 1
		ww = rmat[0][0];	wu = rmat[0][1];	wv = rmat[0][2]
		uw = rmat[1][0];	uu = rmat[1][1];	uv = rmat[1][2]
		vw = rmat[2][0];	vu = rmat[2][1];	vv = rmat[2][2]
	elseif (w>u && w>v)					// w is the biggest
		id = 2
		vv = rmat[0][0];	vw = rmat[0][1];	vu = rmat[0][2]
		wv = rmat[1][0];	ww = rmat[1][1];	wu = rmat[1][2]
		uv = rmat[2][0];	uw = rmat[2][1];	uu = rmat[2][2]
	endif

	Variable qw = sqrt(1 + uu + vv + ww)/2
	qw = min(1,qw)						// in case round off makes qw>1
	qw = numtype(qw) ? 0 : qw				// in case roudoff make sum negative
	if (qw==1)								// identity matrix, retrurn known answer
		quat = {1,0,0,0}
		return 0
	endif
	Variable w4 = 4*qw
	if (id==0)
		quat = { qw, (wv-vw)/w4, (uw-wu)/w4, (vu-uv)/w4 }
	elseif(id==1)
		quat = { qw, (vu-uv)/w4, (wv-vw)/w4, (uw-wu)/w4 }
	elseif(id==2)
		quat = { qw, (uw-wu)/w4, (vu-uv)/w4, (wv-vw)/w4 }
	endif
End
//
//Function RotationMat2Quaternion(rmat,quat)// fill in the quaternion from thee rotaion matrix, returns angle (degree)
//	Wave rmat								// 3x3 rotation matrix, input
//	Wave quat								// quaternion (4-vector), result
//
//	Variable u=abs(rmat[0][0]), v=abs(rmat[1][1]), w=abs(rmat[2][2])
//	Variable uu,uv,uw,vu,vv,vw,wu,wv,ww	// 9 matrix elements
//	if (u>v && u>w)						// u is the biggest
//		uu = rmat[0][0];	uv = rmat[0][1];	uw = rmat[0][2]
//		vu = rmat[1][0];	vv = rmat[1][1];	vw = rmat[1][2]
//		wu = rmat[2][0];	wv = rmat[2][1];	ww = rmat[2][2]
//	elseif (v>u && v>w)					// v is the biggest
//		ww = rmat[0][0];	wu = rmat[0][1];	wv = rmat[0][2]
//		uw = rmat[1][0];	uu = rmat[1][1];	uv = rmat[1][2]
//		vw = rmat[2][0];	vu = rmat[2][1];	vv = rmat[2][2]
//	elseif (w>u && w>v)					// w is the biggest
//		vv = rmat[0][0];	vw = rmat[0][1];	vu = rmat[0][2]
//		wv = rmat[1][0];	ww = rmat[1][1];	wu = rmat[1][2]
//		uv = rmat[2][0];	uw = rmat[2][1];	uu = rmat[2][2]
//	endif
//
//	Variable r = sqrt(1 + uu*uu - vv*vv - ww*ww)	// identity matrix, no rotation
//	if (r==0)
//		quat = {1,0,0,0}
//		return 0
//	endif
//	quat = { (wv-vw)/(2*r), r/2, (uv+vu)/(r*r), (wu+uw)/(r*r) }
////	quat[0] = (wv - vw)/(2*r)
////	quat[1] = r/2
////	quat[2] = (uv + vu)/(r*r)
////	quat[3] = (wu + uw)/(r*r)
//End






//Function testQuaternion()
//	Variable err
//	Make/N=3/D/FREE rvec = {PI/4,0,0}	//45¡ rotation about x-axis
//	printWave(rvec,name="rvec (45¡ about X)",brief=1)
//	if (!testQuaternion1(1,0,0,rvec))
//		print "rotation of 100 is OK"
//	else	
//		return 1
//	endif
//	if (!testQuaternion1(0,1,0,rvec))
//		print "rotation of 010 is OK"
//	else
//		return 1
//	endif
//	if (!testQuaternion1(0,0,1,rvec))
//		print "rotation of 001 is OK"
//	else
//		return 1
//	endif
//
//	rvec = {0,PI/4,0}						// 45¡ about y-axis
//	printWave(rvec,name="rvec (45¡ about Y)",brief=1)
//	if (!testQuaternion1(0,1,0,rvec))
//		print "rotation of 010 is OK"
//	else
//		return 1
//	endif
//	if (!testQuaternion1(1,0,0,rvec))
//		print "rotation of 100 is OK"
//	else
//		return 1
//	endif
//	if (!testQuaternion1(0,0,1,rvec))
//		print "rotation of 001 is OK"
//	else
//		return 1
//	endif
//
//	rvec = {0,0,PI/4}						// 45¡ about z-axis
//	printWave(rvec,name="rvec (45¡ about Z)",brief=1)
//	if (!testQuaternion1(0,0,1,rvec))
//		print "rotation of 001 is OK"
//	else
//		return 1
//	endif
//	if (!testQuaternion1(1,0,0,rvec))
//		print "rotation of 100 is OK"
//	else
//		return 1
//	endif
//	if (!testQuaternion1(0,1,0,rvec))
//		print "rotation of 010 is OK"
//	else
//		return 1
//	endif
//
//	print " "
//	print "ALL OK!"
//	return 0
//End
////
//Static Function testQuaternion1(x0,y0,z0,rvec)
//	Variable x0,y0,z0						// test point
//	Wave rvec
//	Make/N=3/D/FREE in={x0,y0,z0}, out
//	Variable printit = strlen(GetRTStackInfo(2))==0
//
//	Make/N=4/D/FREE quat1,quat2
//	Make/N=(3,3)/D/FREE rmat
//
//	RotationVector2Quaternion(rvec,quat1)
//	Quaternion2RotationMat(quat1,rmat)	// fill in the rotaion matrix from the quaternion
//	rmat = abs(rmat)<1e-14 ? 0 : rmat
//	RotationMat2Quaternion(rmat,quat2)
//
//	MatrixOP/FREE delta = sum(magSqr(quat1-quat2))
//	Variable err = delta[0]>1e-13
//	if (printIt || err)
//		printWave(rvec,name="rvec",brief=1)
//		printWave(quat1,name="quat1",brief=1)
//		printWave(quat2,name="quat2",brief=1)
//		printWave(rmat,name="rmat",brief=1)
//		print " "
//		MatrixOP/FREE out = rmat x in
//		printf "rmat x %s  -->  %s\r",vec2str(in),vec2str(out)
//		out = in
//		QuarternionRotate(quat1,out)	;	out = abs(out)<1e-14 ? 0 : out
//		printf "quat1 x %s  -->  %s\r",vec2str(in),vec2str(out)
//		out = in
//		QuarternionRotate(quat2,out)	;	out = abs(out)<1e-14 ? 0 : out
//		printf "quat2 x %s  -->  %s\r",vec2str(in),vec2str(out)
//	endif
//	return err
//End
