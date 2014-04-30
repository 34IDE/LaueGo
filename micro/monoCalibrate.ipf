#pragma rtGlobals=1		// Use modern global access method.
#pragma version = 1.61


//Function formSignal()
//	Wave keV=keV
//	Wave seconds=seconds
//	Wave after=after
//	Wave before=before
//	Wave monoDial=monoDial
//	Wave darkAfter=darkAfter
//	Wave darkBefore=darkBefore
//	Duplicate/O after signal
//	signal = ln(after - darkAfter*seconds) / (before - darkBefore*seconds)
//End



//	epicsDial = {18.4798,  14.6029,  9.712,  8.8219}
//	¥epicsDial-= 10.9651
//	¥printwave(epicsdial)
//	  epicsDial = {7.51467,  3.6378,  -1.2531,  -2.1432}
//	¥printwave(energy)
//	  energy = {8.9813,  11.568,  18.0028,  20.003}



Menu "micro-mono"
	"Enter and Fit Mono Calib Data", EnterEdgeDataAndFit()
	"Set Starting Parameters", MonoParametersPanel()
	"Show Comparison Table", Make_Table_compare()
	SubMenu "   Edges"
		"Show Fitted Edges Table", ShowFittedEdgesTable()
		"Fit and Fix Edge Plot", FitFixUpEdgePlot()
	End
End



Function EnterEdgeDataAndFit()		// returns 0 if edge data entered correctly
	Variable Ndata=6				// number of data points
	Variable rbv=1
	Variable offset=1.07			// probable values

	Prompt Ndata, "number of data points (at least 1)",popup,"1;2;3;4;5;6;7;8;9;10;11"
	Prompt rbv, "theta drive in RBV or DIAL?", popup, "RBV values;DIAL values"
	Prompt offset, "epics motor offset if RBV values (ignored for dial)"
	DoPrompt "input data", Ndata,rbv,offset
	if (V_flag)
		return -1
	endif
	if (Ndata<1)
		DoAlert 0,"cannot proceed, no data entered"
		return -1
	endif
	rbv = (rbv==1)
	if (rbv)
		print "data entered as readback values (RBV), with motor offset of ",offset
	else
		print "the dial values entered directly, (motor offset not used)"
	endif
	Make/N=(Ndata)/O/D energy, epicsDial
	epicsDial = NaN
	energy = NaN
	if (Ndata==3)
		energy = {8.9813,11.5681,20.016}				// probable values
	elseif (Ndata==4)
		energy = {8.9813,11.5681,18.0028,20.003}		// probable values
	elseif (Ndata==7)
		energy = {20.016, 18.014, 13.078, 11.947, 11.568, 9.669, 8.9813}
	endif

	if (Ndata==4)									// example for testing
		rbv=1
		epicsDial = {18.4798,  14.6029,  9.712,  8.8219}		// with offset=10.9651
		// epicsDial = {7.51467,  3.6378,  -1.2531,  -2.1432}	// with offset=10.9651
	elseif (Ndata==3)								// example for testing
		energy = {8.9813,11.5681,20.003}		// probable values
		rbv=1
		epicsDial = {18.5025,14.6307,8.8512}	// with offset=10.9651

	elseif (Ndata==5 && exists("root:Packages:monoCalibrate:measured") && exists("root:Packages:monoCalibrate:monoDial")) // example for testing
		Redimension/N=5 energy,epicsDial
		Wave trueEdge = root:Packages:monoCalibrate:trueEdge
		Wave wDial = root:Packages:monoCalibrate:monoDial
		epicsDial = wDial
		energy = trueEdge
		rbv = 0
	elseif (Ndata==7)								// example for testing
		epicsDial = {9.3259, 10.2143, 13.5408, 14.6802, 15.1105, 17.7551, 18.9827}
		rbv=1
	elseif (Ndata==10)
		epicsDial = {19.2281,  17.8937,  14.585,  16.6797,  20.7086,  22.3626,  13.6226,  14.0279,  9.14064,  8.24561}
		energy = {8.333,  8.979,  11.103,  9.659,  7.709,  7.112,  11.919,  11.564,  17.998,  20.016}
		rbv=0
	elseif (Ndata==6)
		epicsDial = {20.210828,16.315065,11.789337,10.358793,9.241537,7.609617}
		energy = {7.90816, 9.8852, 13.8393, 15.8163, 17.7934, 21.7474}
		rbv=0
		SetDimLabel 0,0,$"(444)",energy,epicsDial
		SetDimLabel 0,1,$"(555)",energy,epicsDial
		SetDimLabel 0,2,$"(777)",energy,epicsDial
		SetDimLabel 0,3,$"(888)",energy,epicsDial
		SetDimLabel 0,4,$"(999)",energy,epicsDial
		SetDimLabel 0,5,$"(11,11,11)",energy,epicsDial
	endif

	SetScale d 0,0,"mm", epicsDial
	SetScale d 0,0,"keV", energy

	ShowFittedEdgesTable()
	if (Ndata==6)
		Edit/W=(354,219,671,438)/K=1 epicsDial.ld,energy as "Close To Proceed"
		ModifyTable alignment=1,size(Point)=9,width(Point)=24,size(epicsDial.y)=12,width(epicsDial.y)=116
		ModifyTable format(epicsDial.y)=3,digits(epicsDial.y)=6,width(epicsDial.l)=66
		ModifyTable title(epicsDial.y)="theta Drive (mm)",title(epicsDial.l)="Si (hkl)"
		ModifyTable size(energy)=12,format(energy)=3,digits(energy)=4,width(energy)=94,title(energy)="energy (keV)"
	else
		Edit/W=(354,219,635,435)/K=1 epicsDial,energy as "Close To Proceed"
		ModifyTable alignment=1,size(Point)=9,width(Point)=24,size(epicsDial)=12,width(epicsDial)=116
		ModifyTable title(epicsDial)="theta Drive (mm)",size(energy)=12,width(energy)=94
		ModifyTable title(energy)="energy (keV)"
		ModifyTable format(epicsDial)=3,digits(epicsDial)=6
		ModifyTable format(energy)=3,digits(energy)=4
	endif
	DoWindow /C InputTable
	PauseForUser InputTable

	// do some checks
	Ndata = numpnts(energy)
	if (numpnts(energy)!=numpnts(epicsDial))
		Abort "different number of energys and positions"
	endif
	print "data entered as"
	printWave(epicsDial)
	printWave(energy)
	if (rbv)
		epicsDial -= offset
	endif

	Duplicate/O energy,motor_cor,angle,theta,lambda,calc_keV,delta_eV
	SetScale d 0,0,"¡", theta,angle
	SetScale d 0,0,"Angstrom", lambda
	SetScale d 0,0,"keV", calc_keV
	SetScale d 0,0,"eV", delta_eV

	//	set up for the fit
	if (!WaveExists(fit_energy))
		Make/N=200/O fit_energy
		WaveStats/Q epicsDial
		SetScale/I x V_min,V_max,"mm", fit_energy
		SetScale d 0,0,"keV", fit_energy
	endif
	Make_Graph_fit()
	WaveStats/Q epicsDial
	Variable N=V_npnts
	WaveStats/Q energy
	N = (V_npnts<N) ? V_npnts : N		// number of points I can fit
	// w[0] = motor offset
	// w[1] = correction to scale factor
	// w[2] = correction to theta (deg)
	// w[3] = fudgeCurve
	// w[4] = yaw (deg),   yaw of the monochromator (means theta axis not perpendicular to beam)
	String Hstr
	if (N<1)
		DoAlert 0,"cannot proceed, no data entered"
		return -1
	elseif (N==1)
		Hstr="01111"
	elseif (N==2)
		Hstr="00111"
	elseif (N==3)
		Hstr="00111"
	elseif (N==4)
		Hstr="00011"
	else
		Hstr="00001"
	endif
	if (!WaveExists(ww))
		Make/O/D ww={1.07917,0.998282,-.858946,0.000114654,1}
	endif

	if (char2num(Hstr[4,4])>48)
		DoAlert 0, "The yaw parameter will not be fit.  It does no good anyhow"
		print "The yaw parameter will not be fit.  It does no good anyhow"
	endif

ww = {1,  1,  -.5,  0.0001,  1e-12}
ww = {1.07917,0.998282,-.858946,0.000114654,1e-12}
ww = {1.08494, 0.99735,-.855,0.00027,1e-12}
ww = {2.1, 0.99735,-.855,0.00027,1e-12}				// the 2.1 comes from the piece of glass glued in encoder
	if (exists("ww_start")==1 && DimSize(ww_start,0))
		Wave ww_start=ww_start
		ww = ww_start
	endif

	Variable H0,H1,H2,H3,H4
	H0 = str2num(Hstr[0,0])+1
	H1 = str2num(Hstr[1,1])+1
	H2 = str2num(Hstr[2,2])+1
	H3 = str2num(Hstr[3,3])+1
	H4 = str2num(Hstr[4,4])+1
	Prompt H0 "motor offset",popup,"Yes;No"
	Prompt H1 "correction to scale factor",popup,"Yes;No"
	Prompt H2 "correction to theta (deg)",popup,"Yes;No"
	Prompt H3 "fudgeCurve",popup,"Yes;No"
	Prompt H4 "yaw (deg),   yaw of the monochromator (means theta axis not perpendicular to beam)",popup,"Yes;No"
	if (N>4)
		DoPrompt "parameters to fit", H0,H1,H2,H3,H4
	elseif (N>3)
		DoPrompt "parameters to hold", H0,H1,H2,H3
	elseif (N>3)
		DoPrompt "parameters to hold", H0,H1,H2
	elseif (N>3)
		DoPrompt "parameters to hold", H0,H1
	elseif (N>3)
		DoPrompt "parameters to hold", H0
	endif
	sprintf Hstr,"%d%d%d%d%d", H0-1,H1-1,H2-1,H3-1,H4-1
	if (H3-1)
		ww[3] = 1e-12
	endif

	print "Hstr=",Hstr
	printWave(ww)

	if (!WaveExists(sigma))
		Make/N=(Ndata)/D sigma=0.001
	else
		Redimension/N=(Ndata) sigma
	endif
	WaveStats/Q sigma							// make sure that none of the sigma are zero
	if (V_min<=0)
		sigma = .001
		DoAlert 0, "found some sigma=0, set them all to 0.001"
	endif

	FuncFit/M=2/H=Hstr mono ww energy /X=epicsDial /W=sigma /I=1 /D /R 
	// assume that all parameters less than 1e-10 are really zero
	ww = (abs(ww[p])<1e-10) ? 0 : ww[p]
	Wave Res_energy=Res_energy
	Res_energy = energy - (mono(ww,epicsDial))
	ModifyGraph zero(Res_Left)=2
	Label Res_Left "error  (\\U)"
	Edit/K=1/W=(8,51,525,214) M_Covar
	Execute "ModifyTable format(M_Covar)=3,digits(M_Covar)=8,width(M_Covar)=92,width(Point)=32"

	// change Res_energy from keV to eV
	if (exists("Res_energy") && !cmpstr(WaveUnits(Res_energy,1),"keV"))
		SetScale d 0,0,"eV", Res_energy
		Res_energy *= 1000
	endif

	NVAR fudgeCurve=root:Packages:monoCalibrate:fudgeCurve
	if (fudgeCurve)
		Make/N=100/O yy					// for plotting fudge factor
		SetScale/I x 5,15,"mm", yy
		SetScale d 0,0,"mm", yy
		yy = fudgeCurve*(x-6.6)*(x-13)
		Make_Graph_fudge()
	endif
	UpdateTable()
	LayoutResultsFunc(GetDataFolder(1))

	WaveStats/Q delta_eV
	printf "The average deviation = %g   ",V_adev
	maxError(delta_eV)
End



Function FitFixUpEdgePlot()
	String edges="8.333;8.979;11.564;17.998;20"
	String edgeNames="Ni K;Cu K;Pt L3;Zr K;Mo K"
	String stdName="NiK;CuK;PtL3;ZrK;MoK"

	GetAxis/Q bottom
	if (V_min==V_max)
		Abort "no axis"
	endif

	Variable edge,i=0, N=ItemsInList(edges)
	do
		edge = str2num(StringFromList(i,edges))
		if (V_min<=edge && edge<=V_max)
			break
		endif
		i += 1
	while(i<N)
	if (i>=N)
		Abort "x-axis does not bracket an absorption edge"
	endif

	String edgeName = StringFromList(i,edgeNames)
	String sig = "signal"
	Prompt sig, "wave with measured signal",popup,TraceNameList("",";",1)
	Prompt edgeName, "which edge to use",popup,edgeNames
	DoPrompt "signal",edgeName,sig
	if (V_Flag)
		return 1
	endif
	i = WhichListItem(edgeName,edgeNames)
	if (i<0)
		Abort "invalid edge selected"
	endif
	Wave signal = $sig
	if (!WaveExists(signal))
		Abort "the wave '"+sig+"' is not in the current data folder"
	endif

//		ModifyGraph gfSize=18, tick=2, mirror(bottom)=1, minor=1, lowTrip=0.001, standoff=0
//		ModifyGraph mode(signal)=3, marker(signal)=19,rgb(signal)=(0,0,65280), msize(signal)=2

	Make/N=4/O/D W_coef
	WaveStats/Q signal
	W_coef[0] = V_min				// offset
	W_coef[1] = V_max-V_min		// amp
	W_coef[2] = 0.001				// xo, always start with 1 eV
	W_coef[3] = i					// tells which edge to fit
	Wave keV = XWaveRefFromTrace("",sig)
	if (WaveExists(keV))
		FuncFit/H="0001"/Q EdgeFit W_coef signal /X=keV /D 
	else
		FuncFit/H="0001"/Q EdgeFit W_coef signal /D 
	endif
	Variable measured = edge + W_coef[2]	// measured edge value
	printf "edge is at %g keV,  measured at %g keV (delta = %g eV)\r",edge,measured,W_coef[2]*1000
	String text, str
	sprintf text, "%s at %g keV\rmeasured at %.4f keV\r",StringFromList(i,edgeNames),edge,measured
	text += SelectString(W_coef[2],"low by","off by","high by")
	sprintf str, " %.1f eV", 1000*abs(W_coef[2])
	text += str

	if (exists("monoDial") && exists("keV"))
		Variable dial = interp(measured,keV,monoDial)
		print "monoDial = ",dial
		sprintf str, "\redge at dial = %.4f",dial
		text += str
		if (exists("root:Packages:monoCalibrate:measured") && exists("root:Packages:monoCalibrate:monoDial"))
			Wave wmeasured = root:Packages:monoCalibrate:measured
			Wave wDial = root:Packages:monoCalibrate:monoDial
			Wave wTrue = root:Packages:monoCalibrate:trueEdge
			wmeasured[i] = measured
			wdial[i] = dial
			wTrue[i] = edge
		endif
	endif
	TextBox/C/N=text0/F=0/S=3/A=LT/X=5.7/Y=7 text

	SetDrawLayer /K UserFront
	SetDrawLayer UserFront
	SetDrawEnv xcoord= bottom,ycoord= prel,dash= 1
	DrawLine edge,0,edge,1
	SetDrawEnv xcoord= bottom,ycoord= prel,dash= 1,linefgc= (65535,0,0)
	DrawLine measured,0,measured,1
End

Function EdgeFit(w,x) : FitFunc
	Wave w
	Variable x

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(x) = offset + amp * edge(x-x0)
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ x
	//CurveFitDialog/ Coefficients 4
	//CurveFitDialog/ w[0] = offset
	//CurveFitDialog/ w[1] = amp
	//CurveFitDialog/ w[2] = x0
	//CurveFitDialog/ w[3] = edgeNumber

	Variable i = round(w[3])			// needed because Igor changes w[3] even if it is fixed!
	if (i==0)
		Wave edge = root:Packages:monoCalibrate:Nickel:NiK
	elseif (i==1)
		Wave edge = root:Packages:monoCalibrate:Copper:CuK
	elseif (i==2)
		Wave edge = root:Packages:monoCalibrate:Platinum:PtL3
	elseif (i==3)
		Wave edge = root:Packages:monoCalibrate:Zirconium:ZrK
	elseif (i==4)
		Wave edge = root:Packages:monoCalibrate:Moly:MoK
	else
		Abort "i not in the range [0,4] in EdgeFit"
	endif
	return w[0] + w[1] * edge(x-w[2])
End



Function UpdateTable()
	Variable hc=12.39841857					// hc in keV-Angstroms
	Wave epicsDial=epicsDial
	Wave energy=energy
	Wave ww=ww
	Wave motor_cor=motor_cor
	Wave angle=angle
	Wave theta=theta
	Wave lambda=lambda
	Wave calc_keV=calc_keV
	Wave delta_eV=delta_eV
	motor_cor = (epicsDial[p]+ww[0])*ww[1]
	angle = motor2angle(motor_cor[p])+ww[2]
	theta = angle2theta(angle)
	lambda = theta2lambda(theta[p])
	calc_keV = hc/lambda
	delta_eV = 1000*(energy - calc_keV)
	Make_Table_compare()
End


Function maxError(wav)
	Wave wav
	WaveStats/Q wav
	Variable maxE = max(abs(V_max),abs(V_min))
	printf "max error from zero in %s is %g\r", NameOfWave(wav),maxE
	return maxE
End



Function angle2theta(angle)
	Variable angle				// angle of crystal from mechanism (deg)
	NVAR yaw=root:Packages:monoCalibrate:yaw
	return 180/PI*asin(cos(yaw*PI/180)*sin(angle*PI/180))
End



Function mono(w,mot) : FitFunc		// returns energy (keV)
	Wave w
	// w[0] = motor offset
	// w[1] = correction to scale factor
	// w[2] = correction to theta (deg)
	// w[3] = fudgeCurve
	// w[4] = yaw (deg),   yaw of the monochromator (means theta axis not perpendicular to beam)
	Variable mot			// motor position (mm)
	Variable hc=12.39841857					// hc in keV-Angstroms
	NVAR fudgeCurve=root:Packages:monoCalibrate:fudgeCurve
	fudgeCurve = w[3]

	mot = (mot+w[0])*w[1]
	Variable angle = motor2angle(mot)

	angle += w[2]
	if (angle<=0)
		return 1e4
	endif

	// add the yaw correction, w[4]
	NVAR yaw=root:Packages:monoCalibrate:yaw
	yaw = w[4]
	Variable theta = angle2theta(angle)		// corrects for yaw

	Variable lambda = theta2lambda(theta)
	return hc/lambda
End





Function motor2angle(motor)
	Variable motor					// motor position (mm)

	Variable mm = 25.4				// usually 25.4
	Variable h = 0.688*mm
	Variable v = 3.25*mm
	Variable d = 13.763*mm
	Variable del = 0.246*mm
	Variable vdel = v-del
	Variable theta0 = abs(atan(h/v))

	Variable r = sqrt(v*v + h*h)						// |distance from P1 to P2|
	Variable c = sqrt(d*d + del*del)					// |distance from P2 to P3|

	// put in the curvature
	NVAR fudgeCurve=root:Packages:monoCalibrate:fudgeCurve
	motor = motor + fudgeCurve*(motor-6.6)*(motor-13)

	// the above variable is fixed at all angles (set by size of parts)
	// all following variables vary with the motor

	Variable l = sqrt((motor+d+h)^2 + vdel*vdel)		// |distance between P1 and P3|
	Variable eta = asin(vdel/l)
	Variable gamma = acos((r*r + l*l - c*c)/(2*l*r))
	Variable angle = PI/2-eta-gamma-theta0			// Bragg angle (except for yaw correction)

	NVAR V_print=root:Packages:monoCalibrate:V_print
	if (V_print)
		printf "h = %g\",         v = %g\"         d = %g\",         del = %g\",         vdel = %g\"\r",h/mm,v/mm,d/mm,del/mm,vdel/mm
		printf "gamma = %g¡,         eta = %g¡,         theta0 = %g¡,         angle = %g¡\r",gamma*180/PI,eta*180/PI,theta0*180/PI,angle*180/PI
		Variable alpha = acos((r*r+d*d-l*l)/(2*r*d))
		Variable beta = acos((l*l+d*d-r*r)/(2*l*d))
		printf "alpha = %g¡,     beta = %g¡    alpha+beta+gamma=%g¡\r",alpha*180/PI,beta*180/PI,(alpha+beta+gamma)*180/PI
		printf "P1 = (%g\",%g\")\r",0,0
		printf "P2 = (%g\",%g\")\r",r*cos(-PI/2+theta0+angle)/mm,r*sin(-PI/2+theta0+angle)/mm
		printf "P3 = (%g\",%g\")\r",(d+h+motor)/mm,(-v+del)/mm
	endif
	return angle*180/PI								// return in degrees
End


Function angle2motor(angle)
	Variable angle					// Bragg angle, except for yaw (locally radians, returns degrees)
	angle *= PI/180				// radians internally, degrees externally

	Variable mm = 25.4				// usually 25.4
	Variable h = 0.688*mm
	Variable v = 3.25*mm
	Variable d = 13.763*mm
	Variable del = 0.246*mm
	Variable vdel = v-del
	Variable theta0 = atan(h/v)
	Variable motor										// motor position (mm)

	Variable a = -PI/2 +theta0 + angle					// angle to P2
	Variable r = sqrt(v*v + h*h)						// distance from P1 to P2
	Variable c2 = d*d + del*del							// distance from P2 to P3 squared

	motor = sqrt(c2 - (r*sin(a)+vdel)^2) - h - d + r*cos(a)

	// take out the curvature
	NVAR fudgeCurve=root:Packages:monoCalibrate:fudgeCurve
	Variable bb = 1/(2*fudgeCurve) - 9.8				// b/2 in quadratic equation
	motor = -bb +sign(fudgeCurve)* sqrt( bb^2 + motor/fudgeCurve - 85.8)

	return motor
End



Function theta2lambda(angle)		// convert Bragg angle to wavelength
	Variable angle								//Bragg angle (degree)
	Variable re=2.817940285e-5				//Classical electron radius (Angstrom)
	Variable ao=5.43102088					//lattice constant of Si () at 22.5¡ C
	NVAR T_C = root:Packages:monoCalibrate:TemperatureC
	ao = ao*(1 + (T_C-22.5)*2.56e-6)		// dL/L=2.56e-6,   Temperature is 19C
	Variable Z=14								// atomic number of Si (no. of electrons)
	Variable N = 8/ao^3							// no. of Si atoms per cubic Angstrom
	Variable dspacing = ao/sqrt(3)				// compute dspacing ()
	angle *= PI/180
	Variable sintheta = sin(angle)
	Variable lambda = 2. * dspacing * sintheta	// approximate wavelength
	if (sintheta<=0)
		return lambda
	endif
	Variable delta = re*N*lambda^2*Z/(2*PI)	// n = 1-delta
	lambda = 2. * dspacing * (1-delta/(sintheta^2)) * sintheta	// corrected wavelength
	return lambda
End

Function lambda2theta(lambda)		// convert wavelength to Bragg angle (degree)
	Variable lambda								// wavelength ()
	Variable angle
	Variable re=2.817940285e-5				//Classical electron radius (Angstrom)
	Variable ao=5.43102088					//lattice constant of Si () at 22.5¡ C
	NVAR T_C = root:Packages:monoCalibrate:TemperatureC
	ao = ao*(1 + (T_C-22.5)*2.56e-6)		// dL/L=2.56e-6,   Temperature is 19C
	Variable Z=14								// atomic number of Si (no. of electrons)
	Variable N = 8/ao^3							// no. of Si atoms per cubic Angstrom
	Variable dspacing = ao/sqrt(3)				// compute dspacing ()
	Variable delta = re*N*lambda^2*Z/(2*PI)	// n = 1-delta
	Variable l2d = lambda/(2.*dspacing)
	Variable sinth = ( l2d + sqrt(l2d^2+4*delta) ) / 2
	angle = (sinth>0) ? asin(sinth) : 0			// avoid problem with lambda > 2d
	return (angle*180./PI)
End


Function ShowFittedEdgesTable()
	if (strlen(WinList("Table_FittedEdges", ";", "WIN:2"))>0)
		DoWindow /F Table_FittedEdges
		return 0
	endif
	Edit/W=(168,49,468,199)/K=1 root:Packages:monoCalibrate:monoDial,root:Packages:monoCalibrate:measured,root:Packages:monoCalibrate:trueEdge
	Execute "ModifyTable width(Point)=32,format=4,digits=5"
	Execute "ModifyTable width(Point)=32,format(Point)=1"
	DoWindow /C Table_FittedEdges
End


Function Make_Table_compare()
	if (strlen(WinList("Table_compare", ";", "WIN:2"))>0)
		DoWindow /F Table_compare
		return 1
	endif

	Edit/K=1/W=(596,405,1095,576) epicsDial,motor_cor,angle,theta,lambda,calc_keV,energy
	AppendToTable delta_eV
	Execute "ModifyTable font=\"Helvetica\",alignment=1,size=10,width(Point)=24,width(epicsDial)=54"
	Execute "ModifyTable width(motor_cor)=60,width(angle)=54,width(theta)=54,width(lambda)=60"
	Execute "ModifyTable width(calc_keV)=56,width(energy)=48,width(delta_eV)=72,size(Point)=9"
	Execute "ModifyTable format(epicsDial)=3,digits(epicsDial)=6"
	Execute "ModifyTable format(energy)=3,digits(energy)=4"
	Execute "ModifyTable format(delta_eV)=3,digits(delta_eV)=2,width(delta_eV)=50"
	DoWindow /C Table_compare
End



Function Make_Graph_fit()
	if (strlen(WinList("Graph_fit", ";", "WIN:1"))>0)
		DoWindow /F Graph_fit
		return 1
	endif
	Display /W=(425,61,987,411)/K=1 energy vs epicsDial
	AppendToGraph fit_energy
	ModifyGraph mode(energy)=3,marker(energy)=19
	ModifyGraph lSize(fit_energy)=2,rgb(energy)=(0,3,65280), tick=2
	ModifyGraph mirror=1,minor=1,lowTrip=0.001,standoff=0
	ModifyGraph lblPos(left)=35
	Label left "energy  (\\U)"
	Label bottom "motor position  (\\U)"
//	SetAxis left 8,25
//	SetAxis bottom -4,8
	ShowInfo
	DoWindow/C/R Graph_fit
End

Function Make_Graph_fudge()
	if (strlen(WinList("Graph_fudge", ";", "WIN:1"))>0)
		DoWindow /F Graph_fudge
		return -1
	endif
	Display /W=(25,56,564,346)/K=1 yy
	ModifyGraph lSize=2,tick=2,zero(left)=2,mirror=1,minor=1
	ModifyGraph lowTrip=0.001,standoff=0
	Label left "motor fudge (\\U)"
	Label bottom "motor position (\\U)"
	SetAxis left -0.02,0.02
	DoWindow /C Graph_fudge
End


Function LayoutResultsFunc(fldr)
	String fldr					// folder to use for all waves and variables
	PauseUpdate
	NewLayout /C=1/W=(125,45,572,553)/K=1

	TextBox/N=text0/F=0/S=3/A=LB/X=19.86/Y=93.03
	AppendText "\\Z12\\{\"motor(internal mm) = (epics_motor + %g mm) * %g\", "+fldr+"ww[0],"+fldr+"ww[1]}"
	AppendText "\\{\"BraggAngle¡ = motor2theta(epicsDial) + %g¡\", "+fldr+"ww[2]}"
	AppendText "\\{\"fudge curvature = %g\", "+fldr+"ww[3]}"

	AppendLayoutObject table  Table_compare
	ModifyLayout left(Table_compare)=59, top(Table_compare)=75, width(Table_compare)=484, height(Table_compare)=138

	TextBox/N=text1/F=0/S=3/A=LB/X=75.96/Y=65.16 "\\Z12h = +0.688\"\rd = +13.763\"\rv = +3.25\"\rdel = +0.246\""
	TextBox/N=text2/S=1/A=LB/X=29.44/Y=63.11 "\\Z12\\{\"yum:m58:c2:m1.OFF = %g\", "+fldr+"ww[0]}"
	AppendText "\\{\"yum:mono:thetaOff = %g¡\", "+fldr+"ww[2]}"
	AppendText "\\{\"yum:mono:fudgeGain = %g\", "+fldr+"ww[1]}"
	AppendText "\\{\"yum:mono:fudgeCurve = %.6g\", "+fldr+"ww[3]}"
	NVAR yaw=root:Packages:monoCalibrate:yaw
	if (yaw>1e-9)
		AppendText "\\{\"yum:mono:yaw = %.6g\", "+fldr+"ww[4]}"
	endif

	AppendLayoutObject graph  Graph_fit
	ModifyLayout left(Graph_fit)=75,top(Graph_fit)=297,width(Graph_fit)=449,height(Graph_fit)=243
	if (!cmpstr("Graph_fudge",WinList("Graph_fudge","","WIN:1")))
		AppendLayoutObject graph  Graph_fudge
		ModifyLayout left(Graph_fudge)=127,top(Graph_fudge)=555,width(Graph_fudge)=348,height(Graph_fudge)=170
	endif

	TextBox/N=stamp0/F=0/A=RB/X=0.18/Y=0.14 "\\Z06\\{\"%s %s\",date(), time()}"
	TextBox/N=stamp1/F=0/A=LB/X=0.18/Y=0.14 "\\Z06\\{\"%s\",CornerStamp1_()}:Graph_fit"
	ModifyLayout mag=.5
End
Window Layout_results() :Layout
	PauseUpdate; Silent 1		// building window...
	String list = StringByKey("FOLDERS", DataFolderDir(1))
	list = "root,"+RemoveFromList("Packages", list,",")
	String fldr = pickFolder(list)
	if (strlen(fldr)>1 && !stringmatch(fldr,"root"))
		fldr="root:"+fldr
	endif
	if (!stringmatch(fldr,"*:"))
		fldr += ":"
	endif
	LayoutResultsFunc(fldr)
EndMacro
Function/T pickFolder(list)
	String list
	String fldr=""
	Variable i
	do
		i = strsearch(list, ",",0)
		if (i<0)
			break
		endif
		list[i,i]=";"
	while (1)
	Prompt fldr, "folder name", popup,list
	DoPrompt "pick a folder", fldr
	return fldr
End


Function MonoParametersPanel()
	monoCalibrateInitPackage()

	if (strlen(WinList("PanelMonoParameters","","WIN:64")))
		DoWindow/F PanelMonoParameters
		return 0
	endif

	NVAR motorOff = root:Packages:monoCalibrate:motorOff_start
	NVAR thetaOff = root:Packages:monoCalibrate:thetaOff_start
	NVAR gainFactor = root:Packages:monoCalibrate:gainFactor_start
	NVAR fudgeCurve = root:Packages:monoCalibrate:fudgeCurve_start
	NVAR yaw = root:Packages:monoCalibrate:yaw_start
	// w[0] = motor offset
	// w[1] = correction to scale factor
	// w[2] = correction to theta (deg)
	// w[3] = fudgeCurve
	// w[4] = yaw (deg),   yaw of the monochromator (means theta axis not perpendicular to beam)

	if (exists("ww_start")!=1)
		Make/D ww_start = {motorOff,gainFactor,thetaOff,fudgeCurve,yaw}
	else
		Wave ww_start = ww_start
		motorOff = ww_start[0]
		gainFactor = ww_start[1]
		thetaOff = ww_start[2]
		fudgeCurve = ww_start[3]
		yaw = ww_start[4]
	endif

	NewPanel /W=(721,45,972,201)/K=1
	DoWindow/C PanelMonoParameters
	SetVariable setvar_motorOff,pos={14,5},size={191,21},title="motor offset"
	SetVariable setvar_motorOff,help={"motor offset"},fSize=14
	SetVariable setvar_motorOff,limits={-inf,inf,0},value= root:Packages:monoCalibrate:motorOff_start,bodyWidth= 100
	SetVariable setvar_thetaOff,pos={20,30},size={186,21},title="theta offset"
	SetVariable setvar_thetaOff,help={"theta offset"},fSize=14
	SetVariable setvar_thetaOff,limits={-inf,inf,0},value= root:Packages:monoCalibrate:thetaOff_start,bodyWidth= 100
	SetVariable setvar_gainOff,pos={28,55},size={178,21},title="gain factor"
	SetVariable setvar_gainOff,help={"gain factor"},fSize=14
	SetVariable setvar_gainOff,limits={-inf,inf,0},value= root:Packages:monoCalibrate:gainFactor_start,bodyWidth= 100
	SetVariable setvar_fudgeCurve,pos={22,80},size={185,21},title="fudge curve"
	SetVariable setvar_fudgeCurve,help={"fudge curve"},fSize=14
	SetVariable setvar_fudgeCurve,limits={-inf,inf,0},value= root:Packages:monoCalibrate:fudgeCurve_start,bodyWidth= 100
	SetVariable setvar_yaw,pos={76,105},size={131,21},title="yaw",help={"yaw"}
	SetVariable setvar_yaw,fSize=14
	SetVariable setvar_yaw,limits={-inf,inf,0},value= root:Packages:monoCalibrate:yaw_start,bodyWidth= 100
	Button buttonSetWW,pos={51,128},size={123,23},proc=ButtonSetWWProc,title="Set w_start"
End
//
Function ButtonSetWWProc(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			NVAR motorOff = root:Packages:monoCalibrate:motorOff_start
			NVAR thetaOff = root:Packages:monoCalibrate:thetaOff_start
			NVAR gainFactor = root:Packages:monoCalibrate:gainFactor_start
			NVAR fudgeCurve = root:Packages:monoCalibrate:fudgeCurve_start
			NVAR yaw = root:Packages:monoCalibrate:yaw_start
			Make/D/O ww_start = {motorOff,gainFactor,thetaOff,fudgeCurve,yaw}
			break
	endswitch
	return 0
End


Function monoCalibrateInitPackage()
	NewDataFolder /O root:Packages
	NewDataFolder /O root:Packages:monoCalibrate

	if (exists("root:Packages:Elements:edgeStrings")!=1)
		Execute/P "INSERTINCLUDE  \"Elements\", version>=1.0"
		Execute/P "COMPILEPROCEDURES "
		Execute/P "ElementDataInitPackage()"
	endif

	if (exists("root:Packages:monoCalibrate:TemperatureC")!=2)
		Variable /G root:Packages:monoCalibrate:TemperatureC=19
		NVAR T_C = root:Packages:monoCalibrate:TemperatureC
		printf "setting the monochromator temperature to %g¡C\r", T_C
	endif
	if (exists("root:Packages:monoCalibrate:measured")!=1)
		Make/N=5/D root:Packages:monoCalibrate:measured
	endif
	if (exists("root:Packages:monoCalibrate:monoDial")!=1)
		Make/N=5/D root:Packages:monoCalibrate:monoDial
	endif
	if (exists("root:Packages:monoCalibrate:trueEdge")!=1)
		Make/N=5/D root:Packages:monoCalibrate:trueEdge
	endif

	if (exists("root:Packages:monoCalibrate:motorOff_start")!=2)
		Variable /G root:Packages:monoCalibrate:motorOff_start=2
	endif
	if (exists("root:Packages:monoCalibrate:thetaOff_start")!=2)
		Variable /G root:Packages:monoCalibrate:thetaOff_start=-0.9884
	endif
	if (exists("root:Packages:monoCalibrate:gainFactor_start")!=2)
		Variable /G root:Packages:monoCalibrate:gainFactor_start=0.99796
	endif
	if (exists("root:Packages:monoCalibrate:fudgeCurve_start")!=2)
		Variable /G root:Packages:monoCalibrate:fudgeCurve_start=0.99796
	endif
	if (exists("root:Packages:monoCalibrate:yaw_start")!=2)
		Variable /G root:Packages:monoCalibrate:yaw_start=0.99796
	endif


	if (exists("root:Packages:monoCalibrate:fudgeCurve")!=2)
		Variable /G root:Packages:monoCalibrate:fudgeCurve=-0.00009
	endif
	if (exists("root:Packages:monoCalibrate:yaw")!=2)
		Variable /G root:Packages:monoCalibrate:yaw=0.01
	endif
	Variable /G root:Packages:monoCalibrate:V_print=0
End