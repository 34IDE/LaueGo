#pragma rtGlobals=1		// Use modern global access method.
#pragma version = 3.3
#pragma ModuleName=ionChamber

// Note on usage:	Jon Tischler, ORNL
//	If you want to interactively get a result, call IonChamberAsk() or use the "ionChamber" entry in the "Analysis" menu.  
//	If you want to do a more automatic analysis, use Io_ionChamber(...) to the get flux incident on an ion chamber,
//	and optionally, call IonChamber_detector(...) to get the intensity projected all the way back to the detector

// Sept27,2003     Andrei Tkachuk.
//  Added control panel.  All previous functionality is intact (see above)
//  Note on new usage:
//  Select  menu item in toolbar : Analysis->Ion Chamber             to use new interface (Panel)
//  Select menu item in toolbar:   Analysis->Ion Chamber_old      to use interactive prompt input
//
// July18, 2007
//	changed xxx() and xxx1() to static function to remove conflicts
//  also added MoculeName to the top
//
// Apr 30, 2014
//	changed line terminations CR -> LF

Menu "Analysis"
       "Ion Chamber",BuildIonChamberPanel()
       help = {"calculate intensities using ion chamber data using new panel interface"}
	"Ion Chamber_old", IonChamberAsk()
	 help = {"calculate intensities using ion chamber data using old interface"}
End

Proc BuildIonChamberPanel()
       if (wintype("IonChamberPanel") == 0)
	            ionChamberInitPackage()
	            IonChamberPanel()
	        else    
	             DoWindow/F IonChamberPanel
	  endif           	
end

//	computes flux that would be measured at the detector thru an ion chamber
//
//	Program calculates the flux through an ion chamber by the function photon()
//	The flux found is that which would be measured in the detector.
Function IonChamberAsk()		// dives the ion chamber calculation
	Variable energy=9.			// photon energy (KeV)
	Variable VperA=1.e7		// gain of Keithley (V/A)
	Variable cps=1.e5			// #counts/sec output by the VF converter
	Variable aln=60.			// length of active region of ion chamber
	Variable dln1=17.5			// length of upstream dead region of ion chamber
	Variable dln2=17.5			// length of downstream dead region of ion chamber
	Variable Kapton1=.001		// thickness of kapton window before active region
	Variable Kapton2=.001		// thickness of Kapton window on rear of ion chamber
	Variable Kapton3=.0		// thickness of Kapton in beam path after ion chamber (but not including rear window)
	Variable dHe=0.				// He path after ion chamber (was 160.)
	Variable dAir=0.			// air path after ion chamber (was 310.)
	Variable dBe=0				// thickness of all Be windows after the ion chamber (inches)
	Variable dEff=1.			// detector efficiency at this energy in range of [0,1] (if <0, then inches of Ar at 1atm)
	Variable id1=4,id2=0,id3=0,id4		// gas type in ion chamber, 1=He, 2=N, 3=Ar, 4=Air, 7=CO2
	Variable frac1=1,frac2,frac3,frac4		// fraction of id1,id2, id3, and id4
	Variable flux				// incident flux (photons/sec) before ion chamber
	Variable detector			// result, photons/sec measured in detector
	String idList				// gas mix as "id:frac" pairs, id 1=He, 2=N, 3=Ar, 4=Air, 7=CO2, (1atm air = "4:1")

	Prompt cps, "Counts/sec collected in the Io scaler"			// counts = volts*VtoF*cntTime
	Prompt VperA, "Keithley amplifier gain (V/A)"
	Prompt energy, "X-ray energy (KeV) (<0 for wavelength)"
	DoPrompt "Ion Chanber calc",cps,VperA,energy
	if (V_flag)
		return 0
	endif

	// All internal thickness and distnces are in cm
	Prompt Kapton1, "Kapton before active region (inches)"
	Prompt aln, "active length of ion chamber (mm)"
	Prompt dln1, "dead length at front of ion chamber (mm)"
	Prompt dln2, "dead length at rear of ion chamber (mm)"
	Prompt Kapton2, "Kapton on rear of ion chamber (inches)"
	DoPrompt "ionChamber",Kapton1,aln,dln1,dln2,Kapton2
	if (V_flag)
		return 0
	endif

	// ion chamber fill gas mixture
	frac1 = (id1>0) ? frac1 : 0
	frac2 = (id2>0) ? frac2 : 0
	frac3 = (id3>0) ? frac3 : 0
	frac4 = (id4>0) ? frac4 : 0
	Prompt id1,"ion chamber fill gas 1", popup, "Helium;Nitrogen;Argon;Air;Carbon Dioxide"
	Prompt id2,"ion chamber fill gas 2", popup, "none;Helium;Nitrogen;Argon;Air;Carbon Dioxide"
	Prompt id3,"ion chamber fill gas 3", popup, "none;Helium;Nitrogen;Argon;Air;Carbon Dioxide"
	Prompt id4,"ion chamber fill gas 4", popup, "none;Helium;Nitrogen;Argon;Air;Carbon Dioxide"
	Prompt frac1, "fraction of gas1 (atmospheres)"
	Prompt frac2, "fraction of gas2 (atmospheres)"
	Prompt frac3, "fraction of gas3 (atmospheres)"
	Prompt frac4, "fraction of gas4 (atmospheres)"
	DoPrompt "ionChamber",id1,frac1,id2,frac2,id3,frac3,id4,frac4
	if (V_flag)
		return 0
	endif
	id2 -= 1							// allow for "none" in id2, id3, id4
	id3 -= 1
	id4 -= 1
	id1 = (id1==5) ? 7 : id1			// CO2 is 7, not 5
	id2 = (id2==5) ? 7 : id2
	id3 = (id3==5) ? 7 : id3
	id4 = (id4==5) ? 7 : id4
	frac1 = (id1>0) ? frac1 : 0			// if id<1, then force frac to zero
	frac2 = (id2>0) ? frac2 : 0
	frac3 = (id3>0) ? frac3 : 0
	frac4 = (id4>0) ? frac4 : 0
	id1 = (frac1>0) ? id1 : 0			// if frac less than or equal to 0, then set id to zero
	id2 = (frac2>0) ? id2 : 0
	id3 = (frac3>0) ? id3 : 0
	id4 = (frac4>0) ? id4 : 0
	if ((frac1+frac2+frac3+frac4 )<= 0) 
		Abort "No fill gas"
	endif
	idList = MakeIdStringList(id1,frac1,id2,frac2,id3,frac3,id4,frac4)

	// path length after the ion chamber
	Prompt dHe,"He path after ion chamber (mm)"
	Prompt dAir, "Air path after ion chamber (mm)"
	Prompt dBe"Be after ion chamber (inches)"
	Prompt Kapton3, "Kapton in beam path after, but not including rear of ion chamber (inches)"
	DoPrompt "Path after IonChamber",dHe,dAir,dBe,Kapton3
	if (V_flag)
		return 0
	endif

	Variable isProportional=2
	Prompt isProportional, "detector is an Ar proportional counter", popup,"is proportional counter; is not"
	Prompt dEff, "efficiency if NOT proportional counter, or inches of Ar at 1atm in proportional"
	DoPrompt "detector efficiency", isProportional,dEff
	if (V_flag)
		return 0
	endif
	if (isProportional==1)
		dEff = (dEff>0) ? -dEff : dEff	// <0 flags to use proportional counter calculation
		dEff = (dEff==0) ? -4 : dEff		// 2inch dia and 2atm gives 4
	endif


	flux = Io_ionChamber(cps,VperA,energy,idList,Kapton1,aln,dln1)
       print cps,VperA,energy,idList,Kapton1,aln,dln1
//	Calculate absorption from front of ion chamber to detector
	detector = IonChamber_detector(flux,energy,idList,Kapton1,dln1+aln+dln2,Kapton2,dHe,dAir,dBe,Kapton3,dEff)
	printf "%g photons before ion chamber --> %g photons in detector\r",flux,detector
	print energy,idList,Kapton1,dln1+aln+dln2,Kapton2,dHe,dAir,dBe,Kapton3,dEff
	print energy, flux, detector
	return detector
End




Function Io_ionChamber(counts,VperA,energy,idList,Kapton1,aln,dln1)	// calculates intensity before an ion chamber
	Variable counts			// #counts output by the VF converter (=1.e5)
	Variable VperA			// gain of Keithley (V/A) (=1.e8)
	Variable energy			// photon energy (KeV) (=9)
	String idList			// gas mix as "id:frac" pairs, id 1=He, 2=N, 3=Ar, 4=Air, 7=CO2, (1atm air = "4:1")
	Variable Kapton1		// thickness of kapton window before active region (=0.001)
	Variable aln				// length of active region of ion chamber (mm) (=60.)
	Variable dln1			// length of upstream dead region of ion chamber (mm) (=17.5)

	Variable incidentFlux		// result, photons/sec
	Variable mu					// mu of fill gas mix (1/cm)
	
	Variable hc=12.39842435	// KeV-Angstroms
	Variable VtoF=1.e5			// output gain of the voltage-to-frequency converter
	Variable cntTime=1.		// count time in seconds

	incidentFlux = photon(counts,cntTime,VtoF,VperA,aln/10,idList,energy)
	mu = muOfList(idList,energy)							// mu of fill gas (1/cm)
	incidentFlux *= exp(mu*dln1/10)						// absorption of dead length at front
	incidentFlux *= exp(absorb(6,energy)*Kapton1*2.54)	// absorption of Kapton
	return incidentFlux
End



Function IonChamber_detector(flux,energy,idList,Kapton1,ionLength,Kapton2,dHe,dAir,dBe,Kapton3,dEff)
	Variable flux				// incident flux (photons/sec) before ion chamber
	Variable energy				// photon energy (KeV), 8.0546
	String idList				// gas mix as "id:frac" pairs, id 1=He, 2=N, 3=Ar, 4=Air, 7=CO2, (1atm air = "4:1")
	Variable Kapton1			// thickness of kapton window before active region (.001)
	Variable ionLength			// total length of active region of ion chamber (mm) (=95.)
	Variable Kapton2			// thickness of Kapton window on rear of ion chamber (=.001)
	Variable dHe				// He path after ion chamber (mm) (160.)
	Variable dAir				// air path after ion chamber (mm) (310.)
	Variable dBe				// thickness of all Be windows after the ion chamber (inches)
	Variable Kapton3			// thickness of Kapton in beam path after ion chamber (but not including rear window) (=.001)
	Variable dEff				// detector efficiency at this energy in range of [0,1] (if <0, then inches of Ar at 1atm) (usually 1)

	Variable mu				// mu of fill gas mix (1/cm)
	Variable AirAbs			// absorption of dAir mm in air
	Variable HeAbs			// absorption of dHe mm in He
	Variable KapAbs			// absorption of all the Kapton
	Variable BeAbs			// absorption of all the Be windows
	Variable ionAbs			// absorption in the fill gas of the ion chamber
	Variable detector		// result, photons/sec measured in detector

	if (dEff<0)				// for proportional counter, dEff is length of Ar at 1atm (2"dia, 2atm Ar is 4)
		dEff = abs(dEff)
		dEff =  1 - exp(-dEff*2.54*absArPhoto(energy))	// photo-electric part of Argon absorption (4" at 1atm)
		dEff *= exp(-.005*2.54*absorb(5,energy))		// 5 mil Be window
	endif

	// All internal thickness and distnces are in cm
	//	Calculate absorption from front of ion chamber to detector
	mu = muOfList(idList,energy)					// mu of fill gas (1/cm)
	ionAbs = exp(-mu*ionLength/10)				// absorption due to fill gas
	AirAbs = exp(-dAir/10 * absorb(4,energy) )
	HeAbs = exp(-dHe/10 * absorb(1,energy) )
	KapAbs = exp(-(Kapton1+Kapton2+Kapton3)*2.54 * absorb(6,energy) )
	BeAbs = exp(-dBe*2.54 * absorb(5,energy) )
	detector = flux * AirAbs*HeAbs*KapAbs*BeAbs*ionAbs * dEff
	return detector
End





// calculates the number of photons incident from the measured current
Function photon(counts,cntTime,VtoF,VperA,activeL,idList,energy)
	Variable counts			// the number of pulses output by the VtoF converter
	Variable cntTime		// count time in sec
	Variable VtoF			// pulses/sec/volt output by the V-F converter
	Variable VperA			// current amplifier gain in volts/amp
	Variable activeL			// active length of the ion chamber (cm)
	String idList			// gas mix as "id:frac" pairs, id 1=He, 2=N, 3=Ar, 4=Air, 7=CO2, (1atm air = "4:1")
	Variable energy			// enegy  (KeV)

	Variable photo=0		// photo-electron absorption length  (1/cm)
	Variable work			// work function in eV/ion pair for the ion chamber gas
//								suggested values are,  He: 29.6 eV or 27.8 eV
//								N2: 36.3 eV,  O2: 32.2 eV,  CO2: 33.5 eV,  Ne: 27.4 eV
//								Ar: 24.4 eV,  Kr: 22.8 eV,  Xe: 20.8 eV
	Variable id, frac
	Variable convHe	=6.647
	Variable rhoHe	=0.0001785
	Variable convN	=23.26
	Variable rhoN	=0.00125
	Variable convAr	=66.32
	Variable rhoAr	=0.001784
	Variable edgeAr	=3.202

	Make/O root:Packages:ionChamber:aHe= { 6.06488,-3.29055,-0.107256, 0.0144465}
	Make/O root:Packages:ionChamber:aN = { 11.2765,-2.65400,-0.200445, 0.0200765}
	Make/O root:Packages:ionChamber:aAr1={ 13.9491,-1.82276,-0.328827, 0.0274382}
	Make/O root:Packages:ionChamber:aAr2={ 12.2960,-2.63279,-0.073660, 0.0      }
	Make/O root:Packages:ionChamber:aO = { 11.7130,-2.57229,-0.205893, 0.0199244}
	Make/O root:Packages:ionChamber:aC = {10.6879, -2.71400, -0.200530, 0.0207248}
	Wave aHe=root:Packages:ionChamber:aHe
	Wave aN=root:Packages:ionChamber:aN
	Wave aAr1=root:Packages:ionChamber:aAr1
	Wave aAr2=root:Packages:ionChamber:aAr2
	Wave aO=root:Packages:ionChamber:aO
	Wave aC=root:Packages:ionChamber:aC

	Variable E1,E2,E3
	Variable total
	Variable part1,part2

	E1 = ln(energy)
	E2 = E1*E1
	E3 = E2*E1

	Variable photoSum=0				// photo-electron absorption length  (1/cm)
	Variable work_photo=0				// sum of mu/work
	Variable i=0
	String str
	do
		str = StringFromList(i,idList)
		id = str2num(StringFromList(0, str,":"))
		frac = str2num(StringFromList(1, str,":"))
		if (numtype(id) || numtype(frac))
			break
		endif

		//	Calculate the photo-electric cross section for the ion chamber gas

		if (id==1)						// Photo-electric cross section for Helium
			total = exp(aHe[0] + aHe[1]*E1 + aHe[2]*E2 + aHe[3]*E3)
			photo = total*rhoHe/convHe
			work_photo += photo*frac / 29.6
			photoSum += photo*frac
		elseif (id==2)					// Photo-electric cross section for Nitrogen
			total = exp( aN[0] +  aN[1]*E1 +  aN[2]*E2 +  aN[3]*E3)
			photo = total *rhoN /convN
			work_photo += photo*frac / 36.3
			photoSum += photo*frac
		elseif (id==3)					// Photo-electric cross section for Argon
			if (E1>edgeAr)
				total = exp(aAr1[0]+aAr1[1]*E1 +aAr1[2]*E2 +aAr1[3]*E3)
			else
			total = exp(aAr2[0]+aAr2[1]*E1 +aAr2[2]*E2 +aAr2[3]*E3)
			endif
			photo = total*rhoAr/convAr
			work_photo += photo*frac / 24.4
			photoSum += photo*frac
		elseif (id==4)					// Photo-electric cross section for Air .79N2, .20O2, .01Ar
			total = exp( aN[0] +  aN[1]*E1 +  aN[2]*E2 +  aN[3]*E3)
			photo = total * 0.000922 / convN
			work_photo += photo*frac / 36.3
			photoSum += photo*frac

			total = exp( aO[0] +  aO[1]*E1 +  aO[2]*E2 +  aO[3]*E3)
			photo = total * 0.000266 / convN		// should use convO here, but I do not know it
			work_photo += photo*frac / 36.3
			photoSum += photo*frac

			if (E1>edgeAr)
				total = exp(aAr1[0]+aAr1[1]*E1 +aAr1[2]*E2 +aAr1[3]*E3)
			else
				total = exp(aAr2[0]+aAr2[1]*E1 +aAr2[2]*E2 +aAr2[3]*E3)
			endif
			photo = total * 1.66E-5 / convAr
			work_photo += photo*frac / 24.4
			photoSum += photo*frac
		elseif (id==7)					// Photo-electric cross section for CO2
			photo = 0
			total = exp( aC[0] +  aC[1]*E1 +  aC[2]*E2 +  aC[3]*E3)
			photo += total * 0.0005396 / convN		// should use convC here, but I do not know it
			total = exp( aO[0] +  aO[1]*E1 +  aO[2]*E2 +  aO[3]*E3)
			photo += total * 0.0014374 / convN		// should use convO here, but I do not know it
			work_photo += photo*frac / 35.36
			photoSum += photo*frac
		else
			Abort "illegal gas id of "+num2str(id)
		endif
		i += 1
	while(1)

//	Parts 1 and 2 calculate the flux in photons/sec
	part1 = counts/cntTime/(1.602176462e-19*VtoF*VperA*energy*1000.)
	part2 = (1 - exp(-photoSum*activeL)) * work_photo / photoSum
	KillWaves/Z aHe, aN aAr1, aAr2, aO
	return(part1/part2)
End







Function absArPhoto(energy)		// Linear absorption coeff for photo-electric part of Argon (cm^2/gm)
	Variable energy					// enegy  (KeV)
	Variable conv=66.32
	Variable rho=0.001784		// g/cm^3
	Variable E1,E2,E3
	Variable photo
	Variable edge=3.202			// K-edge in KeV
	Make/O root:Packages:ionChamber:a1_abs_ = { 1.39491E+01,-1.82276E+00,-3.28827E-01, 2.74382E-02 }
	Make/O root:Packages:ionChamber:a2_abs_ = { 1.22960E+01,-2.63279E+00,-7.36600E-02, 0.}
	Wave a1=root:Packages:ionChamber:a1_abs_, a2=root:Packages:ionChamber:a2_abs_

	E1 = ln(energy)
	E2 = E1*E1
	E3 = E2*E1

	if (energy>edge)
		photo = exp( a1[0] + a1[1]*E1 + a1[2]*E2 + a1[3]*E3)
	else
		photo = exp( a2[0] + a2[1]*E1 + a2[2]*E2 + a2[3]*E3)
	endif
	return(photo/conv*rho)
End



//***********************************************************************



//double absH(double energy)		// for Hydrogen, calculate mass absorption coefficient
//double absHe(double energy)		//  "  Helium
//double absBe(double energy)		//  "  Beryllium
//double absC(double energy)		//  "  Carbon
//double absN(double energy)		//  "  Nitrogen
//double absO(double energy)		//  "  Oxygen
//double absAr(double energy)		//  "  Argon

Function absorb(id,energy)			// linear absorption coefficient (1/cm)
	Variable id						// integer specifying type of absorber
	Variable energy					// enegy  (KeV)
// This routine calculates the absorbtion coefficient in 1/cm for
// 	a chosen material and a given photon energy.
// 
// 	Variable definitions:
// 
// 	id = integer value which specifies the xray absorber
// 		1 = Helium
// 		2 = Nitrogen
// 		3 = Argon
// 		4 = Air (dry)
// 		5 = Beryllium
// 		6 = Kapton
// 		7 = CO2
//
// energy = photon energy in KeV

	Variable rhoH				// for Hydrogen, density  (gm/cm^3)
	Variable rhoHe				//  "  Helium
	Variable rhoBe				//  "  Beryllium
	Variable rhoC				//  "  Carbon
	Variable rhoN				//  "  Nitrogen
	Variable rhoO				//  "  Oxygen
	Variable rhoAr				//  "  Argon

	switch (id)
	case 1:								// Helium
		rhoHe = 0.0001785
		return(absHe(energy)*rhoHe)
	case 5:								// Beryllium
		rhoBe = 1.848
		return(absBe(energy)*rhoBe)
	case 2:								// Nitrogen
		rhoN  = 0.00125
		return(absN(energy)*rhoN)
	case 4:								// Dry Air
		rhoN  = 0.000922					// 79% N2, 20% O2 AND 1% Ar by volume
		rhoO  = 0.000266					// rho Air = 1.2047E-3 gm/cc at 20C
		rhoAr = 1.66E-5
		return(absN(energy)*rhoN + absO(energy)*rhoO + absAr(energy)*rhoAr)
	case 3:								// Argon
		rhoAr = 0.001784
		return(absAr(energy)*rhoAr)
	case 6:								// Kapton
		rhoC  = 0.981					// C22 H10 O5 N2
		rhoH  = 0.037					// MolWt=382 gm/mole    rho=1.42 gm/cc
		rhoO  = 0.297
		rhoN  = 0.105
		return(absC(energy)*rhoC + absH(energy)*rhoH + absO(energy)*rhoO + absN(energy)*rhoN)
	case 7:								// CO2,   carbon dioxide density = 0.001977 g/cc
		rhoC  = 0.0005396				// 12.01115/44.00995 * 0.001977
		rhoO  = 0.0014374				// 2*15.9994/44.00995 * 0.001977
		return(absC(energy)*rhoC + absO(energy)*rhoO)
	endswitch
	return 0
End


//***********************************************************************


Function absH(energy)			// Mass absorption coefficient for Hydrogen (cm^2/gm)
Variable energy					// enegy  (KeV)
	Variable	conv
	Variable	E1,E2,E3
	Variable	photo			// photo-electron scattering
	Variable	coherent		// coherent scattering
	Variable	compton		// compton scattering

	Make/O root:Packages:ionChamber:a_abs_ = { 2.44964E+00,-3.34953E+00,-4.71370E-02, 7.09962E-03 }
	Make/O root:Packages:ionChamber:b_abs_ = {-1.19075E-01,-9.37086E-01,-2.00538E-01, 1.06587E-02 }
	Make/O root:Packages:ionChamber:c_abs_ = {-2.15772E+00, 1.32685E+00,-3.05620E-01, 1.85025E-02 }
	Wave a=root:Packages:ionChamber:a_abs_, b=root:Packages:ionChamber:b_abs_,  c=root:Packages:ionChamber:c_abs_

//	rho = 8.987E-5
	conv = 1.674

	E1 = ln(energy)
	E2 = E1*E1
	E3 = E2*E1

	photo = exp( a[0] + a[1]*E1 + a[2]*E2 + a[3]*E3)		// Photoelectron
	coherent = exp( b[0] + b[1]*E1 + b[2]*E2 + b[3]*E3)	// Coherent
	compton = exp( c[0] + c[1]*E1 + c[2]*E2 + c[3]*E3)		// Compton
	return( (photo+coherent+compton)/conv)					// Mass absorption coef (cm^2/gm)
End


//***********************************************************************


Function absHe(energy)		// Mass absorption coefficient for Helium (cm^2/gm)
Variable energy				// enegy  (KeV)
	Variable	conv
	Variable	E1,E2,E3
	Variable	photo
	Variable	coherent
	Variable	compton
	Make/O root:Packages:ionChamber:a_abs_ = { 6.06488E+00,-3.29055E+00,-1.07256E-01, 1.44465E-02 }
	Make/O root:Packages:ionChamber:b_abs_ = { 1.04768E+00,-8.51805E-02,-4.03527E-01, 2.69398E-02 }
	Make/O root:Packages:ionChamber:c_abs_ = {-2.56357E+00, 2.02536E+00,-4.48710E-01, 2.79691E-02 }
	Wave a=root:Packages:ionChamber:a_abs_, b=root:Packages:ionChamber:b_abs_,  c=root:Packages:ionChamber:c_abs_

//	rho = 1.785E-04
	conv = 6.647

	E1 = ln(energy)
	E2 = E1*E1
	E3 = E2*E1

	photo = exp( a[0] + a[1]*E1 + a[2]*E2 + a[3]*E3)
	coherent = exp( b[0] + b[1]*E1 + b[2]*E2 + b[3]*E3)
	compton = exp( c[0] + c[1]*E1 + c[2]*E2 + c[3]*E3)
	return((photo+coherent+compton)/conv)
End


//***********************************************************************


Function absBe(energy)		// Mass absorption coefficient for Beryllium (cm^2/gm) 
Variable energy				// enegy  (KeV)
	Variable	conv
	Variable	E1,E2,E3
	Variable	photo
	Variable	coherent
	Variable	compton
	Make/O root:Packages:ionChamber:a_abs_ =  {9.04511E+00,-2.83487E+00,-2.10021E-01, 2.29526E-02 }
	Make/O root:Packages:ionChamber:b_abs_ =  {2.00860E+00,-4.61920E-02,-3.37018E-01, 1.86939E-02 }
	Make/O root:Packages:ionChamber:c_abs_ = {-6.90079E-01, 9.46448E-01,-1.71142E-01, 6.51413E-03 }
	Wave a=root:Packages:ionChamber:a_abs_, b=root:Packages:ionChamber:b_abs_,  c=root:Packages:ionChamber:c_abs_

//	rho = 1.848
	conv = 14.96

	E1 = ln(energy)
	E2 = E1*E1
	E3 = E2*E1

	photo = exp( a[0] + a[1]*E1 + a[2]*E2 + a[3]*E3)
	coherent = exp( b[0] + b[1]*E1 + b[2]*E2 + b[3]*E3)
	compton = exp( c[0] + c[1]*E1 + c[2]*E2 + c[3]*E3)
	return((photo+coherent+compton)/conv)
End


//***********************************************************************


Function absC(energy)		// Mass absorption coefficient for Carbon (cm^2/gm)
Variable energy				// enegy  (KeV)
	Variable	conv
	Variable	E1,E2,E3
	Variable	photo
	Variable	coherent
	Variable	compton
	Make/O root:Packages:ionChamber:a_abs_ = { 1.06879E+01,-2.71400E+00,-2.00530E-01, 2.07248E-02 }
	Make/O root:Packages:ionChamber:b_abs_ = { 3.10861E+00,-2.60580E-01,-2.71974E-01, 1.35181E-02 }
	Make/O root:Packages:ionChamber:c_abs_ = {-9.82878E-01, 1.46693E+00,-2.93743E-01, 1.56005E-02 }
	Wave a=root:Packages:ionChamber:a_abs_, b=root:Packages:ionChamber:b_abs_,  c=root:Packages:ionChamber:c_abs_

//	rho = 1.580
	conv = 19.94

	E1 = ln(energy)
	E2 = E1*E1
	E3 = E2*E1

	photo = exp( a[0] + a[1]*E1 + a[2]*E2 + a[3]*E3)
	coherent = exp( b[0] + b[1]*E1 + b[2]*E2 + b[3]*E3)
	compton = exp( c[0] + c[1]*E1 + c[2]*E2 + c[3]*E3)
	return((photo+coherent+compton)/conv)
End


//***********************************************************************


Function absN(energy)			// Mass absorption coefficient for Nitrogen (cm^2/gm)
	Variable energy				// enegy  (KeV)
	Variable	conv
	Variable	E1,E2,E3
	Variable	photo
	Variable	coherent
	Variable	compton
	Make/O root:Packages:ionChamber:a_abs_ = { 1.12765E+01,-2.65400E+00,-2.00445E-01, 2.00765E-02 }
	Make/O root:Packages:ionChamber:b_abs_ = { 3.47760E+00,-2.15762E-01,-2.88874E-01, 1.51312E-02 }
	Make/O root:Packages:ionChamber:c_abs_ = {-1.23693E+00, 1.74510E+00,-3.54660E-01, 1.98705E-02 }
	Wave a=root:Packages:ionChamber:a_abs_, b=root:Packages:ionChamber:b_abs_,  c=root:Packages:ionChamber:c_abs_

//	rho = 0.001250
	conv = 23.26

	E1 = ln(energy)
	E2 = E1*E1
	E3 = E2*E1

	photo = exp( a[0] + a[1]*E1 + a[2]*E2 + a[3]*E3)
	coherent = exp( b[0] + b[1]*E1 + b[2]*E2 + b[3]*E3)
	compton = exp( c[0] + c[1]*E1 + c[2]*E2 + c[3]*E3)
	return((photo+coherent+compton)/conv)
End


//***********************************************************************


Function absO(energy)		// Mass absorption coefficient for Oxygen (cm^2/gm)
	Variable energy			// enegy  (KeV)
	Variable	conv
	Variable	E1,E2,E3
	Variable	photo
	Variable	coherent
	Variable	compton
	Make/O root:Packages:ionChamber:a_abs_ = { 1.17130E+01,-2.57229E+00,-2.05893E-01, 1.99244E-02 }
	Make/O root:Packages:ionChamber:b_abs_ = { 3.77239E+00,-1.48539E-01,-3.07124E-01, 1.67303E-02 }
	Make/O root:Packages:ionChamber:c_abs_ = {-1.73679E+00, 2.17686E+00,-4.49050E-01, 2.64733E-02 }
	Wave a=root:Packages:ionChamber:a_abs_, b=root:Packages:ionChamber:b_abs_,  c=root:Packages:ionChamber:c_abs_

//	rho = 0.001429
	conv = 26.57

	E1 = ln(energy)
	E2 = E1*E1
	E3 = E2*E1

	photo = exp( a[0] + a[1]*E1 + a[2]*E2 + a[3]*E3)
	coherent = exp( b[0] + b[1]*E1 + b[2]*E2 + b[3]*E3)
	compton = exp( c[0] + c[1]*E1 + c[2]*E2 + c[3]*E3)
	return((photo+coherent+compton)/conv)
End


//***********************************************************************


Function absAr(energy)			// Mass absorption coefficient for Argon (cm^2/gm)
	Variable energy				//enegy  (KeV)

	Variable	conv
	Variable	E1,E2,E3
	Variable	photo
	Variable	coherent
	Variable	compton
	Variable	edge			// K-edge in KeV
	Make/O root:Packages:ionChamber:a1_abs_ = { 1.39491E+01,-1.82276E+00,-3.28827E-01, 2.74382E-02 }
	Make/O root:Packages:ionChamber:a2_abs_ = { 1.22960E+01,-2.63279E+00,-7.36600E-02, 0.}
	Make/O root:Packages:ionChamber:b_abs_ = { 5.21079E+00, 1.35618E-01,-3.47214E-01, 1.84333E-02 }
	Make/O root:Packages:ionChamber:c_abs_ = {-6.82105E-01, 1.74279E+00,-3.17646E-01, 1.56467E-02 }
	Wave a1=root:Packages:ionChamber:a1_abs_, a2=root:Packages:ionChamber:a2_abs_, b=root:Packages:ionChamber:b_abs_,  c=root:Packages:ionChamber:c_abs_

//	rho = 0.001784
	conv = 66.32
	edge = 3.202

	E1 = ln(energy)
	E2 = E1*E1
	E3 = E2*E1

	if (energy>edge)
		photo = exp( a1[0] + a1[1]*E1 + a1[2]*E2 + a1[3]*E3)
	else
		photo = exp( a2[0] + a2[1]*E1 + a2[2]*E2 + a2[3]*E3)
	endif
	coherent = exp( b[0] + b[1]*E1 + b[2]*E2 + b[3]*E3)
	compton = exp( c[0] + c[1]*E1 + c[2]*E2 + c[3]*E3)
	return((photo+coherent+compton)/conv)
End






Function/T MakeIdStringList(id1,frac1,id2,frac2,id3,frac3,id4,frac4)
	Variable id1,frac1,id2,frac2,id3,frac3,id4,frac4
	String out=""
	if (id1>0 && frac1>0)
		out += num2istr(id1)+":"+num2str(frac1)+";"
	endif
	if (id2>0 && frac2>0)
		out += num2istr(id2)+":"+num2str(frac2)+";"
	endif
	if (id3>0 && frac3>0)
		out += num2istr(id3)+":"+num2str(frac3)+";"
	endif
	if (id4>0 && frac4>0)
		out += num2istr(id4)+":"+num2str(frac4)+";"
	endif
	return out
End



Function muOfList(idList,energy)		// returns mu = 1/absorption length (1/cm)
	String idList			// gas mixsure , 1=He, 2=N2, 3=Ar, 4=Air
	Variable energy			// enegy  (KeV)

	Variable id, frac		// individual gas id and fraction of fill (atm)
	Variable mu=0			// 1/absorption length (1/cm)
	Variable i=0
	String str
	do
		str = StringFromList(i,idList)
		id = str2num(StringFromList(0, str,":"))
		frac = str2num(StringFromList(1, str,":"))
		if (numtype(id) || numtype(frac))
			break
		endif
		mu += absorb(id,energy)*frac			// Calculate mu for each component
		i += 1
	while(1)
	return mu
End

//====================================================================================

Function ionChamberInitPackage()
	NewDataFolder /O root:Packages
	NewDataFolder /O root:Packages:ionChamber
	NewDataFolder /O root:Packages:ionChamber:PanelGlobals
	
	if(exists("root:Packages:ionChamber:PanelGlobals:gcounts")!=2)
		Variable/G root:Packages:ionChamber:PanelGlobals:gcounts=1.0e+5        // #counts output by the VF converter (=1.e5)
	endif
	if(exists("root:Packages:ionChamber:PanelGlobals:gVperA")!=2)	
		Variable/G root:Packages:ionChamber:PanelGlobals:gVperA=1.0e+08	      // gain of Keithley (V/A) (=1.e8)
	endif
	if(exists("root:Packages:ionChamber:PanelGlobals:gVtoF")!=2)	
		Variable/G root:Packages:ionChamber:PanelGlobals:gVtoF=1.0e+05	      // output gain of the voltage-to-frequency converter
	endif

	if(exists("root:Packages:ionChamber:PanelGlobals:genergy")!=2)		
		Variable/G  root:Packages:ionChamber:PanelGlobals:genergy=9.			// photon energy (KeV) (=9)
	endif
	if(exists("root:Packages:ionChamber:PanelGlobals:gKapton1")!=2)
		Variable/G root:Packages:ionChamber:PanelGlobals:gKapton1=0.001		// thickness of kapton window before active region (=0.001)
	endif	
	if(exists("root:Packages:ionChamber:PanelGlobals:gKapton2")!=2)
	       Variable/G root:Packages:ionChamber:PanelGlobals:gKapton2=0.001		// thickness of Kapton window on rear of ion chamber
	endif  
	if(exists("root:Packages:ionChamber:PanelGlobals:galn")!=2)     
		Variable/G root:Packages:ionChamber:PanelGlobals:galn=60.	// length of active region of ion chamber (mm) (=60.)
	endif
	if(exists("root:Packages:ionChamber:PanelGlobals:gdln1")!=2)     
		Variable/G root:Packages:ionChamber:PanelGlobals:gdln1=17.5		// length of upstream dead region of ion chamber (mm) (=17.5)
	endif
	if(exists("root:Packages:ionChamber:PanelGlobals:gdln2")!=2)    
	       Variable/G root:Packages:ionChamber:PanelGlobals:gdln2=17.5			// length of downstream dead region of ion chamber
	endif
	if(exists("root:Packages:ionChamber:PanelGlobals:gdHe")!=2)    
      		Variable/G root:Packages:ionChamber:PanelGlobals:gdHe=0		// He path after ion chamber (mm)
      	endif
      	if(exists("root:Packages:ionChamber:PanelGlobals:gdAir")!=2)	
		Variable/G root:Packages:ionChamber:PanelGlobals:gdAir=762.		// air path after ion chamber (mm)
	endif
	if(exists("root:Packages:ionChamber:PanelGlobals:gdBe")!=2)	
		Variable/G root:Packages:ionChamber:PanelGlobals:gdBe	=0.08			// thickness of all Be windows after the ion chamber (inches)
	endif
	if(exists("root:Packages:ionChamber:PanelGlobals:gKapton3")!=2)	
		Variable/G root:Packages:ionChamber:PanelGlobals:gKapton3=0.		// thickness of Kapton in beam path after ion chamber (but not including rear window)
	endif
	if(exists("root:Packages:ionChamber:PanelGlobals:gdEff")!=2)
		Variable/G root:Packages:ionChamber:PanelGlobals:gdEff=1.			// detector efficiency at this energy in range of [0,1] (if <0, then inches of Ar at 1atm)
       endif
	if(exists("root:Packages:ionChamber:PanelGlobals:gID1")!=2)
		Variable/G root:Packages:ionChamber:PanelGlobals:gID1=4			// iD1 of fill gas 1  
       endif
     	if(exists("root:Packages:ionChamber:PanelGlobals:gID2")!=2)
		Variable/G root:Packages:ionChamber:PanelGlobals:gID2=1			// iD1 of fill gas 2  
       endif
	if(exists("root:Packages:ionChamber:PanelGlobals:gID3")!=2)
		Variable/G root:Packages:ionChamber:PanelGlobals:gID3=1			// iD1 of fill gas 3 
       endif
	if(exists("root:Packages:ionChamber:PanelGlobals:gID4")!=2)
		Variable/G root:Packages:ionChamber:PanelGlobals:gID4=1			// iD1 of fill gas 4  
       endif
       if(exists("root:Packages:ionChamber:PanelGlobals:gfrac1")!=2)
		Variable/G root:Packages:ionChamber:PanelGlobals:gfrac1=1			//fraction of fill gas (atm) with gID1 
       endif
       if(exists("root:Packages:ionChamber:PanelGlobals:gfrac2")!=2)
		Variable/G root:Packages:ionChamber:PanelGlobals:gfrac2=0.0			//fraction of fill gas (atm) with gID2 
       endif
       if(exists("root:Packages:ionChamber:PanelGlobals:gfrac3")!=2)
		Variable/G root:Packages:ionChamber:PanelGlobals:gfrac3=0			//fraction of fill gas (atm) with gID3 
       endif
       if(exists("root:Packages:ionChamber:PanelGlobals:gfrac4")!=2)
		Variable/G root:Packages:ionChamber:PanelGlobals:gfrac4=0			//fraction of fill gas (atm) with gID4 
       endif
       if(exists("root:Packages:ionChamber:PanelGlobals:gflux")!=2)
		Variable/G root:Packages:ionChamber:PanelGlobals:gflux=0			//photons before ion chamber 
       endif
       if(exists("root:Packages:ionChamber:PanelGlobals:gdetector")!=2)
		Variable/G root:Packages:ionChamber:PanelGlobals:gdetector=0			// photons in detector
       endif
       if(exists("root:Packages:ionChamber:PanelGlobals:gRadioVal")!=2)
		Variable/G root:Packages:ionChamber:PanelGlobals:gRadioVal=1	// Energy or Wavelengh Panel radio button value (1 for energy)
       endif
       
End


Static Function xxx(ene)

	Variable ene
	Variable flux,detector

	flux =  Io_ionChamber(100000,1e+08,ene,"1:1",0.001,60,17.5)
	detector = IonChamber_detector(flux,ene,"1:1",0.001,95,0.001, 0, 762,  0.08,  0,  1 )
printf "%g keV %g photons before ion chamber --> %g photons in detector\r",ene,flux,detector
	return detector
//    return flux
end

Static Function xxx1(counts,VperA,energy,Kapton1,Kapton2,aln,dln1,dln2,dHe,dAir,dBe,dEff,id1,id2,id3,id4,frac1,frac2,frac3,frac4)
	Variable counts		// #counts output by the VF converter (=1.e5)
	Variable VperA		      // gain of Keithley (V/A) (=1.e8)
	Variable  energy			// photon energy (KeV) (=9)
	Variable Kapton1		// thickness of kapton window before active region (=0.001)
       Variable Kapton2		// thickness of Kapton window on rear of ion chamber
	Variable aln				// length of active region of ion chamber (mm) (=60.)
	Variable dln1			// length of upstream dead region of ion chamber (mm) (=17.5)
       Variable dln2		// length of downstream dead region of ion chamber
      Variable dHe		// He path after ion chamber (mm)
	Variable dAir		// air path after ion chamber (mm)
	Variable dBe			// thickness of all Be windows after the ion chamber (inches)	
	Variable dEff			// detector efficiency at this energy in range of [0,1] (if <0, then inches of Ar at 1atm)
	Variable id1,id2,id3,id4
	Variable frac1,frac2,frac3,frac4  // fraction of id1,id2, id3, and id4
       
	Variable Kapton3=0		// thickness of Kapton in beam path after ion chamber (but not including rear window)
       
       NVAR flux=root:Packages:ionChamber:PanelGlobals:gflux
       NVAR detector=root:Packages:ionChamber:PanelGlobals:gdetector
       
       
       String idList=ProcessGasMixtureInput(id1,frac1,id2,frac2,id3,frac3,id4,frac4)// gas mix as "id:frac" pairs, id 1=He, 2=N, 3=Ar, 4=Air, 7=CO2, (1atm air = "4:1")

//       print idList
	flux =  Io_ionChamber(counts,VperA,energy,idList,Kapton1,aln,dln1)	// calculates intensity before an ion chamber
//	Io_ionChamber(100000,1e+08,ene,"1:1",0.001,60,17.5)
	detector = IonChamber_detector(flux,energy,idList,Kapton1,dln1+aln+dln2,Kapton2, dHe, dAir,  dBe,  Kapton3,  deff )
       printf "%g keV %g photons before ion chamber --> %g photons in detector\r",energy,flux,detector
	return detector
//    return flux
end

Function ButtonProc_CalculatePhotons(ctrlName) : ButtonControl
	String ctrlName
	
	NVAR counts=root:Packages:ionChamber:PanelGlobals:gcounts   // #counts output by the VF converter (=1.e5)
	NVAR VperA=root:Packages:ionChamber:PanelGlobals:gVperA    //// gain of Keithley (V/A) (=1.e8)
	NVAR  energy=root:Packages:ionChamber:PanelGlobals:genergy // photon energy (KeV) (=9)
	NVAR Kapton1=root:Packages:ionChamber:PanelGlobals:gKapton1 // thickness of kapton window before active region (=0.001)
	NVAR Kapton2=root:Packages:ionChamber:PanelGlobals:gKapton2  // thickness of Kapton window on rear of ion chamber
	NVAR aln=root:Packages:ionChamber:PanelGlobals:galn // length of active region of ion chamber (mm) (=60.)
	NVAR dln1=root:Packages:ionChamber:PanelGlobals:gdln1 // length of upstream dead region of ion chamber (mm) (=17.5)
	NVAR dln2=root:Packages:ionChamber:PanelGlobals:gdln2  // length of downstream dead region of ion chamber
	NVAR dHe=root:Packages:ionChamber:PanelGlobals:gdHe // He path after ion chamber (mm)
	NVAR dAir=root:Packages:ionChamber:PanelGlobals:gdAir // air path after ion chamber (mm)
	NVAR dBe=root:Packages:ionChamber:PanelGlobals:gdBe  //// thickness of all Be windows after the ion chamber (inches)
	NVAR Kapton3=root:Packages:ionChamber:PanelGlobals:gKapton3 // thickness of Kapton in beam path after ion chamber (but not including rear window)
	NVAR dEff=root:Packages:ionChamber:PanelGlobals:gdEff   // detector efficiency at this energy in range of [0,1] (if <0, then inches of Ar at 1atm)
	
	NVAR id1=root:Packages:ionChamber:PanelGlobals:gID1
	NVAR id2=root:Packages:ionChamber:PanelGlobals:gID2
	NVAR id3=root:Packages:ionChamber:PanelGlobals:gID3
	NVAR id4=root:Packages:ionChamber:PanelGlobals:gID4
	
	NVAR frac1=root:Packages:ionChamber:PanelGlobals:gfrac1
	NVAR frac2=root:Packages:ionChamber:PanelGlobals:gfrac2
	NVAR frac3=root:Packages:ionChamber:PanelGlobals:gfrac3
	NVAR frac4=root:Packages:ionChamber:PanelGlobals:gfrac4

      Variable energy1   //if wavelength panel radio button is selected instead of energy 
      Variable hc=12.39842435	// KeV-Angstroms
	ControlInfo/W=IonChamberPanel check_wavelength
	if (V_Value==1)
		    energy1=hc/energy
		else
		    energy1=energy    //radio button is set to energy in (keV)
	endif
	
//	dEff = (dEff==0) ? -4 : dEff		// 2inch dia and 2atm gives 4
	
      xxx1(counts,VperA,energy1,Kapton1,Kapton2,aln,dln1,dln2,dHe,dAir,dBe,dEff,id1,id2,id3,id4,frac1,frac2,frac3,frac4)
End

Function SetVarProc_Calculate(ctrlName,varNum,varStr,varName) : SetVariableControl
	String ctrlName
	Variable varNum
	String varStr
	String varName
	
       ButtonProc_CalculatePhotons("button_calculate")    //same as  pushing button on the panel
End


Function/S ProcessGasMixtureInput(id1,frac1,id2,frac2,id3,frac3,id4,frac4)
	Variable id1,frac1, id2,frac2,id3,frac3,id4,frac4		// gas type in ion chamber, 1=He, 2=N, 3=Ar, 4=Air, 7=CO2
              
       String idList       
              	// ion chamber fill gas mixture
	frac1 = (id1>0) ? frac1 : 0
	frac2 = (id2>0) ? frac2 : 0
	frac3 = (id3>0) ? frac3 : 0
	frac4 = (id4>0) ? frac4 : 0

	id2 -= 1							// allow for "none" in id2, id3, id4
	id3 -= 1
	id4 -= 1
	id1 = (id1==5) ? 7 : id1			// CO2 is 7, not 5
	id2 = (id2==5) ? 7 : id2
	id3 = (id3==5) ? 7 : id3
	id4 = (id4==5) ? 7 : id4
	frac1 = (id1>0) ? frac1 : 0			// if id<1, then force frac to zero
	frac2 = (id2>0) ? frac2 : 0
	frac3 = (id3>0) ? frac3 : 0
	frac4 = (id4>0) ? frac4 : 0
	id1 = (frac1>0) ? id1 : 0			// if frac less than or equal to 0, then set id to zero
	id2 = (frac2>0) ? id2 : 0
	id3 = (frac3>0) ? id3 : 0
	id4 = (frac4>0) ? id4 : 0
	if ((frac1+frac2+frac3+frac4 )<= 0) 
		Abort "No fill gas"
	endif
	idList = MakeIdStringList(id1,frac1,id2,frac2,id3,frac3,id4,frac4)

    return idList
end


Function PopMenuProc_FillGasID(ctrlName,popNum,popStr) : PopupMenuControl
	String ctrlName
	Variable popNum
	String popStr
	
	NVAR id1=root:Packages:ionChamber:PanelGlobals:gID1
	NVAR id2=root:Packages:ionChamber:PanelGlobals:gID2
	NVAR id3=root:Packages:ionChamber:PanelGlobals:gID3
	NVAR id4=root:Packages:ionChamber:PanelGlobals:gID4
	
	NVAR frac1=root:Packages:ionChamber:PanelGlobals:gfrac1
	NVAR frac2=root:Packages:ionChamber:PanelGlobals:gfrac2
	NVAR frac3=root:Packages:ionChamber:PanelGlobals:gfrac3
	NVAR frac4=root:Packages:ionChamber:PanelGlobals:gfrac4

	
		strswitch (ctrlName)
		case "popup_gas1":
			id1=popNum
			break
		case "popup_gas2":
			id2=popNum
			break
		case "popup_gas3":
			id3=popNum
			break
		case "popup_gas4":
			id4=popNum
			break
	endswitch
//	print ProcessGasMixtureInput(id1,frac1,id2,frac2,id3,frac3,id4,frac4)

End


Window IonChamberPanel() : Panel
	PauseUpdate; Silent 1		// building window...
       NewPanel /K =1/W=(550,50,980,711) as "Ion Chamber Panel"
       DoWindow/C IonChamberPanel
//	ShowTools
//	SetDrawLayer UserBack
	DrawText 31,29,"Counts collected in the Io scaler"
	DrawText 31,58,"Current amplifier gain (V/A)"
	DrawText 31,85,"Voltage to Frequency gain (Hz/V)"
	SetDrawEnv fillpat= 0
	DrawRect 65,135,380,271
	DrawText 14,313,"Active length of ion chamber"
	DrawText 10,338,"dead length at front of ion chamber "
	DrawText 81,364,"before active region"
	DrawText 10,390,"dead length at rear of ion chamber "
	DrawText 79,417,"after active region"
	DrawText 41,446,"path after ion chamber"
	DrawText 41,476,"path after ion chamber "
	DrawText 41,505,"after ion chamber "
	DrawText 41,505,"after ion chamber "
	SetDrawEnv fsize= 14,fstyle= 1
	DrawText 15,612,"photons before ion chamber"
	SetDrawEnv fsize= 14,fstyle= 1
	DrawText 15,642,"photons in detector"
	DrawText 109,153,"Ion Chamber Fill Gas (atmospheres)"
	SetDrawEnv fsize= 14,fstyle= 1
	DrawText 7,477,"Air"
	SetDrawEnv fsize= 14,fstyle= 1
	DrawText 9,505,"Be"
	SetDrawEnv fsize= 14,fstyle= 1
	DrawText 9,447,"He"
	DrawText 331,419,"(inches)"
	SetDrawEnv fsize= 14,fstyle= 1
	DrawText 10,418,"Kapton"
	SetDrawEnv fsize= 14,fstyle= 1
	DrawText 10,365,"Kapton"
	DrawText 331,501,"(inches)"
	DrawText 331,367,"(inches)"
	DrawText 332,444,"(mm)"
	DrawText 333,472,"(mm)"
	DrawText 331,393,"(mm)"
	DrawText 331,340,"(mm)"
	DrawText 332,312,"(mm)"
	DrawText 9,531,"detector efficiency from 0 to 1"
	DrawText 9,545,"(Note: type -4 for 2\" dia, 2 atm Ar proportional counter)" 
	SetVariable setvar_counts title=" ",pos={255,15},size={129,15},format="%e"
	SetVariable setvar_counts,limits={0,Inf,1000},value= root:Packages:ionChamber:PanelGlobals:gcounts
	SetVariable setvar_VperA title=" ",pos={255,39},size={131,15},format="%e"
	SetVariable setvar_VperA,limits={0,Inf,1000},value= root:Packages:ionChamber:PanelGlobals:gVperA
	SetVariable setvar_VtoF title=" ",pos={255,67},size={130,15},format="%e"
	SetVariable setvar_VtoF,limits={1.0e+5,1.0e+5,0.0e+5},value= root:Packages:ionChamber:PanelGlobals:gVtoF
	SetVariable setvar_energy title=" ",pos={88,112},size={77,15}
	SetVariable setvar_energy,limits={0,Inf,1},value= root:Packages:ionChamber:PanelGlobals:genergy

	CheckBox check_energy,pos={200,112},size={83,14},proc=EnergyOrWavelengthCheckProc,title="energy (keV)",value= 1,mode=1
	CheckBox check_wavelength,pos={303,112},size={104,14},proc=EnergyOrWavelengthCheckProc,title="wavelength (Ang)"
	CheckBox check_wavelength,value= 0,mode=1
       
	PopupMenu popup_gas1,pos={120,160},size={76,20},proc=PopMenuProc_FillGasID
	PopupMenu popup_gas1,mode=root:Packages:ionChamber:PanelGlobals:gid1,value= #"\"Helium;Nitrogen;Argon;Air;Carbon Dioxide\""
	PopupMenu popup_gas2,pos={120,186},size={55,20},proc=PopMenuProc_FillGasID
	PopupMenu popup_gas2,mode=root:Packages:ionChamber:PanelGlobals:gid2,value= #"\"none;Helium;Nitrogen;Argon;Air;Carbon Dioxide\""
	PopupMenu popup_gas3,pos={120,213},size={55,20},proc=PopMenuProc_FillGasID
	PopupMenu popup_gas3,mode=root:Packages:ionChamber:PanelGlobals:gid3,value= #"\"none;Helium;Nitrogen;Argon;Air;Carbon Dioxide\""
	PopupMenu popup_gas4,pos={120,238},size={55,20},proc=PopMenuProc_FillGasID
	PopupMenu popup_gas4,mode=root:Packages:ionChamber:PanelGlobals:gid4,value= #"\"none;Helium;Nitrogen;Argon;Air;Carbon Dioxide\""

	SetVariable setvar_frac1 title=" ",pos={230,161},size={77,15}
	SetVariable setvar_frac1,limits={0,Inf,1},value= root:Packages:ionChamber:PanelGlobals:gfrac1
	SetVariable setvar_frac2 title=" ",pos={229,187},size={77,15}
	SetVariable setvar_frac2,limits={0,Inf,1},value= root:Packages:ionChamber:PanelGlobals:gfrac2
	SetVariable setvar_frac3 title=" ",pos={231,214},size={77,15}
	SetVariable setvar_frac3,limits={0,Inf,1},value= root:Packages:ionChamber:PanelGlobals:gfrac3
	SetVariable setvar_frac4 title=" ",pos={232,241},size={77,15}
	SetVariable setvar_frac4,limits={0,Inf,1},value= root:Packages:ionChamber:PanelGlobals:gfrac4

	SetVariable setvar_aln title=" ",pos={243,298},size={77,15}
	SetVariable setvar_aln,value= root:Packages:ionChamber:PanelGlobals:galn
	SetVariable setvar_dln1 title=" ",pos={243,325},size={77,15}
	SetVariable setvar_dln1,value= root:Packages:ionChamber:PanelGlobals:gdln1
	SetVariable setvar_Kapton1 title=" ",pos={243,352},size={77,15}
	SetVariable setvar_Kapton1,value= root:Packages:ionChamber:PanelGlobals:gKapton1
	SetVariable setvar_dln2 title=" ",pos={243,378},size={77,15}
	SetVariable setvar_dln2,value= root:Packages:ionChamber:PanelGlobals:gdln2
	SetVariable setvar_Kapton2 title=" ",pos={243,404},size={77,15}
	SetVariable setvar_Kapton2,value= root:Packages:ionChamber:PanelGlobals:gKapton2
	SetVariable setvar_dHe title=" ",pos={243,430},size={77,15}
	SetVariable setvar_dHe,value= root:Packages:ionChamber:PanelGlobals:gdHe
	SetVariable setvar_dBe title=" ",pos={243,486},size={77,15}
	SetVariable setvar_dBe,value= root:Packages:ionChamber:PanelGlobals:gdBe
	SetVariable setvar_dAir title=" ",pos={243,456},size={77,15}
	SetVariable setvar_dAir,value= root:Packages:ionChamber:PanelGlobals:gdAir
	Button button_calculate,pos={147,562},size={76,21},proc=ButtonProc_CalculatePhotons,title="calculate"
	SetVariable setvar_PhotonsI0  title=" ",fSize=14,limits={-inf,Inf,0},pos={250,593},size={146,15},value=root:Packages:ionChamber:PanelGlobals:gflux
	SetVariable setvar_PhotonsDet title=" ",fSize=14,limits={-inf,Inf,0},pos={250,623},size={146,15},value= root:Packages:ionChamber:PanelGlobals:gdetector
	SetVariable setvar_dEff,pos={242,515},size={77,15},title=" "
	SetVariable setvar_dEff,limits={-inf,1,0.1},value= root:Packages:ionChamber:PanelGlobals:gdEff
EndMacro


Function EnergyOrWavelengthCheckProc(name,value)
	String name
	Variable value
	
	NVAR gRadioVal= root:Packages:ionChamber:PanelGlobals:gRadioVal    //Energy or Wavelengh Panel radio button value
	
	strswitch (name)
		case "check_energy":
			gRadioVal= 1
			break
		case "check_wavelength":
			gRadioVal= 2
			break
//		case "check2":
//			gRadioVal= 3
//			break
	endswitch
	CheckBox check_energy,value= gRadioVal==1
	CheckBox check_wavelength,value= gRadioVal==2
//	CheckBox check2,value= gRadioVal==3
end