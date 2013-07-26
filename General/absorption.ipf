#pragma rtGlobals=1		// Use modern global access method.
#pragma version = 0.01

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

	Make/O a_abs_ = { 2.44964E+00,-3.34953E+00,-4.71370E-02, 7.09962E-03 }
	Make/O b_abs_ = {-1.19075E-01,-9.37086E-01,-2.00538E-01, 1.06587E-02 }
	Make/O c_abs_ = {-2.15772E+00, 1.32685E+00,-3.05620E-01, 1.85025E-02 }
	Wave a=a_abs_, b=b_abs_,  c=c_abs_

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
	Make/O a_abs_ = { 6.06488E+00,-3.29055E+00,-1.07256E-01, 1.44465E-02 }
	Make/O b_abs_ = { 1.04768E+00,-8.51805E-02,-4.03527E-01, 2.69398E-02 }
	Make/O c_abs_ = {-2.56357E+00, 2.02536E+00,-4.48710E-01, 2.79691E-02 }
	Wave a=a_abs_, b=b_abs_,  c=c_abs_

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
	Make/O a_abs_ =  {9.04511E+00,-2.83487E+00,-2.10021E-01, 2.29526E-02 }
	Make/O b_abs_ =  {2.00860E+00,-4.61920E-02,-3.37018E-01, 1.86939E-02 }
	Make/O c_abs_ = {-6.90079E-01, 9.46448E-01,-1.71142E-01, 6.51413E-03 }
	Wave a=a_abs_, b=b_abs_,  c=c_abs_

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
	Make/O a_abs_ = { 1.06879E+01,-2.71400E+00,-2.00530E-01, 2.07248E-02 }
	Make/O b_abs_ = { 3.10861E+00,-2.60580E-01,-2.71974E-01, 1.35181E-02 }
	Make/O c_abs_ = {-9.82878E-01, 1.46693E+00,-2.93743E-01, 1.56005E-02 }
	Wave a=a_abs_, b=b_abs_,  c=c_abs_

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
	Make/O a_abs_ = { 1.12765E+01,-2.65400E+00,-2.00445E-01, 2.00765E-02 }
	Make/O b_abs_ = { 3.47760E+00,-2.15762E-01,-2.88874E-01, 1.51312E-02 }
	Make/O c_abs_ = {-1.23693E+00, 1.74510E+00,-3.54660E-01, 1.98705E-02 }
	Wave a=a_abs_, b=b_abs_,  c=c_abs_

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
	Make/O a_abs_ = { 1.17130E+01,-2.57229E+00,-2.05893E-01, 1.99244E-02 }
	Make/O b_abs_ = { 3.77239E+00,-1.48539E-01,-3.07124E-01, 1.67303E-02 }
	Make/O c_abs_ = {-1.73679E+00, 2.17686E+00,-4.49050E-01, 2.64733E-02 }
	Wave a=a_abs_, b=b_abs_,  c=c_abs_

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
	Make/O a1_abs_ = { 1.39491E+01,-1.82276E+00,-3.28827E-01, 2.74382E-02 }
	Make/O a2_abs_ = { 1.22960E+01,-2.63279E+00,-7.36600E-02, 0.}
	Make/O b_abs_ = { 5.21079E+00, 1.35618E-01,-3.47214E-01, 1.84333E-02 }
	Make/O c_abs_ = {-6.82105E-01, 1.74279E+00,-3.17646E-01, 1.56467E-02 }
	Wave a1=a1_abs_, a2=a2_abs_, b=b_abs_,  c=c_abs_

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
