#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#pragma ModuleName=BondLattice
#pragma version = 0.03
#include "LatticeSym", version>=7.28

#define USE_BONDS_VERSION_2			// should be 1 or 2, for anything else defaults to the one in AtomView.ipf
// #define USE_BONDS_VERSION_1			// should be 1 or 2, for anything else defaults to the one in AtomView.ipf

Menu "Analysis"
	SubMenu "Lattice"
		"-"
		"Show Bonds Calculated from Xtal", ShowCalculatedBonds_xtal()
		help={"Show Calculatee Bonds for current xtal. Only Prints, does NOT change anything."}
	End
End



#ifdef USE_BONDS_VERSION_2

Static Constant STRUCTURE_ATOMS_MAX=50	// max number of atom types in a material structure (for Si only need 1)

Function ShowCalculatedBonds_xtal()
	STRUCT crystalStructure xtal
	if (FillCrystalStructDefault(xtal))			//fill the lattice structure with test values
		DoAlert 0, "ERROR AnalyzeBonds()\rno lattice structure found"
		return 1
	endif
	AddBondsToXtal(xtal,NaN)

	String type, typeName="match;covalent;ionic"
	Variable i, mult, distCalc
	printf "for '%s'  defined %d bonds (nm):\r",xtal.desc,xtal.Nbonds
	for (i=0;i<xtal.Nbonds;i+=1)
		Make/N=(xtal.bond[i].N)/FREE/D lens
		lens = xtal.bond[i].len[p]
		mult = xtal.bond[i].len[4]		// multiplicity
		distCalc = xtal.bond[i].len[2]	// distCalc
		type = StringFromList(xtal.bond[i].len[3]	,typeName)	// type 0=same atom,  1=covalent,  ±2=ionic
		printf "     %s  <-->  %s:  %s  [mult=%d] '%s'  ∆=%.4f\r",xtal.bond[i].label0,xtal.bond[i].label1,vec2str(lens,bare=1),mult,type,abs(distCalc-lens[0])
	endfor

	if (xtal.Nbonds > 0)
		String list=LatticeSym#UnBondedAtomsList(xtal)
		if (ItemsInList(list)<1)
			print "      All atom types are associated with at least 1 bond"
		else				// print list of those atoms not associated with a bond
			printf "These atom types have NO bonds: {%s}\r",TrimEnd(ReplaceString(";",list,", "),chars=", ")
		endif
	endif





End


Function/WAVE MakeBondList_blen(prefix,xyz,[bLen])	// This is a guess when no bonds are given
	// use with MakeBondList_blen()
	String prefix
	Wave xyz				// list of atom xyz positions
	Variable bLen		// maximum distance that gets a bond
	bLen = ParamIsDefault(bLen) || blen<=0 || numtype(blen) ? LatticeSym_maxBondLen : bLen

	STRUCT crystalStructure xtal
	if (FillCrystalStructDefault(xtal))			//fill the lattice structure with test values
		DoAlert 0, "ERROR AnalyzeBonds()\rno lattice structure found"
		return $""
	endif

	if (!WaveExists(xyz))
		return $""
	elseif (DimSize(xyz,0)<2 || DimSize(xyz,1) < xtal.dim)
		return $""
	endif

	AddBondsToXtal(xtal,blen)
	return AtomView#MakeBondList_Given(prefix,xtal,xyz)	// This makes the bond list from given bond lengths
End


Static Function AddBondsToXtal(xtal,blen)
	// this removes any existing bonds in xtal, and re-computes a new set of bonds
	STRUCT crystalStructure &xtal
	Variable bLen										// maximum distance that gets a bond
	blen = numtype(blen) || blen<=0 ? LatticeSym_maxBondLen : blen
	xtal.Nbonds = 0										// ignore/remove any pre-existing bonds, we are here to find bonds

	Wave enegWave = root:Packages:Elements:electroneg	// use enegWave[Z]
	Wave atomRad = root:Packages:Elements:atomRadius		// use atomRad[Z] (Å)
	Wave covRad = root:Packages:Elements:covRadius		// use covRad[Z] (Å)

	Variable tick0=stopMSTimer(-2)
	Variable i, 	Natoms=xtal.N
	String labelList=""
	for (i=0;i<Natoms;i+=1)
		labelList += xtal.atom[i].name + ";"
	endfor

	Variable MaxBonds = min(Natoms*Natoms,2*STRUCTURE_ATOMS_MAX)
	Make/N=(MaxBonds)/I/FREE atomTypes = -20	// init to unknown
	Variable Nbonds=0
	Make/N=(MaxBonds,8)/D/FREE bondsWave=NaN	// columns are: m0,Z0,m1,Z1,len,mult,type,distCalc
	SetDimLabel 1,0,m0,bondsWave		;	SetDimLabel 1,1,Z0,bondsWave
	SetDimLabel 1,2,m1,bondsWave		;	SetDimLabel 1,3,Z1,bondsWave
	SetDimLabel 1,4,len,bondsWave		;	SetDimLabel 1,5,mult,bondsWave
	SetDimLabel 1,6,itype,bondsWave	;	SetDimLabel 1,7,distCalc,bondsWave

	String label0, label1, type, typeName="match;covalent;ionic"
	Variable Z0,Z1, deltaEneg, distCalc, itype0,itype1, m0,m1, bad, dd
	Variable distCurrent = 1e-4			// 0.001 Å
	for (i=0;i<MaxBonds;i+=1)
		if (NextClosestPair(xtal,distCurrent,blen))	// true if could not find a valid looking bond
			break
		endif
		label0 = xtal.bond[i].label0
		label1 = xtal.bond[i].label1
		Z0 = WhichListItem(label0[0,1], ELEMENT_Symbols)
		Z0 = Z0<0 ? WhichListItem(label0[0], ELEMENT_Symbols)+1 : Z0+1
		Z1 = WhichListItem(label1[0,1], ELEMENT_Symbols)
		Z1 = Z1<0 ? WhichListItem(label1[0], ELEMENT_Symbols)+1 : Z1+1

		deltaEneg = enegWave[Z1] - enegWave[Z0]
		if (deltaEneg==0)
			itype0 = 0						// elements match
			itype1 = 0
		elseif (abs(deltaEneg) < 1.5)	
			itype0 = 1						// covalent
			itype1 = 1
		elseif (deltaEneg > 0)	
			itype0 = 2						// ionic, Z0+
			itype1 = -2						// ionic, Z1-
		else
			itype0 = -2						// ionic, Z0-
			itype1 = 2						// ionic, Z1+
		endif

		m0 = WhichListItem(label0,labelList)	// use pre-existing atomTypes if it was set
		itype0 = atomTypes[m0]<-10 ? itype0 : atomTypes[m0]
		atomTypes[m0] = itype0
		m1 = WhichListItem(label1,labelList)	// use pre-existing atomTypes if it was set
		itype1 = atomTypes[m1]<-10 ? itype1 : atomTypes[m1]
		atomTypes[m1] = itype1
		bad = (abs(atomTypes[m0]) > 1 && abs(atomTypes[m1]) > 1 && atomTypes[m0]*atomTypes[m1] > 0)	// ion type mismatch

		xtal.bond[i].len[3] = abs(itype0)
		distCalc = abs(itype0)<2 ? covRad[Z0] + covRad[Z1] : atomRad[Z0] + atomRad[Z1]
		distCalc /= 10						// convert Å --> nm

		dd = abs(distCalc - xtal.bond[i].len[0])
		bad = dd>0.08 ? 1 : bad			// distances should match withing 0.8Å of calculated (was 0.6Å)
		type = StringFromList(xtal.bond[i].len[3],typeName)
		if (!bad)
			//printf "'%s'[%d] <--> '%s'[%d]  =  %g nm   mult=%d   '%s'   calc=%.3f  ∆=%.4f   %d\r",label0,Z0, label1,Z1, xtal.bond[i].len[0], xtal.bond[i].len[4], type, distCalc, dd,bad
			bondsWave[Nbonds][0] = m0	;	bondsWave[Nbonds][1] = Z0
			bondsWave[Nbonds][2] = m1	;	bondsWave[Nbonds][3] = Z1
			bondsWave[Nbonds][4] = xtal.bond[i].len[0]
			bondsWave[Nbonds][5] = xtal.bond[i].len[4]		// multiplicity
			bondsWave[Nbonds][6] = xtal.bond[i].len[3]		// type 0=same atom,  1=covalent,  ±2=ionic
			bondsWave[Nbonds][7] = distCalc
			Nbonds += 1
		endif
	endfor
	if (Nbonds<1)
		xtal.Nbonds = 0
		printf "NO bonds found  (in %.3g sec):\r", Nbonds, (stopMSTimer(-2)-tick0)*1e-6
		return 1
	endif
	Redimension/N=(Nbonds,-1) bondsWave

	// now check if too many bonds going to one atom
	Variable mult, m
	for (m=0;m<Natoms;m+=1)
		MatrixOP/FREE eq_m = equal(m, col(bondsWave,0))
		MatrixOP/FREE eq_m = eq_m && equal(m, col(bondsWave,2))
		Make/N=(Nbonds)/I/FREE mults = eq_m[p] ? bondsWave[p][5] : 0
		mult = sum(mults)
		if (mult>10)			// too many bonds to this atom, take the <= 10 shortest bonds to atom m
			Variable ilast = BinarySearch(mults, 10)		// last entry to get 10
			if (ilast >= 0)
				mults[0,ilast] = 0		// keeps 0-ilast
			elseif (ilast == -2)		//should never happen
				mults = 0					// keep all
			endif							// if ilast=-2, then delete all

			for (i=0;i<Nbonds;i+=1)
				if (mults[i])			// delete rows with non-zero mult[i]
					DeletePoints/M=0 i, 1, bondsWave
				endif
			endfor
		endif
	endfor
	Nbonds = DimSize(bondsWave,0)

	xtal.Nbonds = Nbonds
	for (i=0;i<Nbonds;i+=1)			// transfer info from bondsWave[i][] --> xtal.bond[i].xxx
		xtal.bond[i].label0 = StringFromList(bondsWave[i][0],labelList)
		xtal.bond[i].label1 = StringFromList(bondsWave[i][2],labelList)
		xtal.bond[i].N = 1
		xtal.bond[i].len[0] = bondsWave[i][4]	// bond length
		xtal.bond[i].len[4] = bondsWave[i][5]	// multiplicity
		xtal.bond[i].len[3] = bondsWave[i][6]	// type
		xtal.bond[i].len[2] = bondsWave[i][7]	// distCalc
	endfor
	// Duplicate/O bondsWave, bondsWaveView		// use:   DisplayTableOfWave(bondsWaveView)
	printf "Computed %d bonds  (in %.3g sec):\r", Nbonds, (stopMSTimer(-2)-tick0)*1e-6
End
//
Static Function NextClosestPair(xtal,distCurrent,blen)
	// find next bond, and add it to xtal.bond[]
	// returns 0 if a bond was added to xtal, 1 if no bond found
	STRUCT crystalStructure &xtal
	Variable distCurrent
	Variable bLen										// maximum distance that gets a bond
	blen = numtype(blen) || blen<=0 ? LatticeSym_maxBondLen : blen

	Wave direct = directFrom_xtal(xtal)
	Variable dim = xtal.dim
	String lj, li, 	ljMin="", liMin=""
	Variable i,j, Natoms=xtal.N
	Make/N=(Natoms)/WAVE/FREE atomRefs
	Make/N=(Natoms)/I/FREE Zs
	Make/N=(Natoms)/D/FREE occupys
	Make/N=(Natoms)/T/FREE labels
	for (i=0;i<Natoms;i+=1)
		Wave atomi = $("root:Packages:Lattices:atom"+num2istr(i))
		Make/N=(DimSize(atomi,0),dim)/D/FREE localWave = atomi
		atomRefs[i] = localWave
		Zs[i] = xtal.atom[i].Zatom
		labels[i] = xtal.atom[i].name
		occupys[i] = xtal.atom[i].occ
	endfor

	Variable multMin, distMin, Ni
	for (j=0,distMin=Inf; j<Natoms; j+=1)
		lj = labels[j]
		Wave xyzj = Extend_XYZ_by1(atomRefs[j])
		MatrixOP/FREE xyzj = (direct x xyzj^t )^t	// xyz of atoms of type atomj in real units
		Wave xyzMid = LatticeSym#FindCentralAtom(xyzj)	// find the centermost atom in xyzj, not average, an actual atom position

		for (i=j;i<Natoms;i+=1)
			li = labels[i]
			if (BondExists(xtal,lj,li))							// skip existing bonds
				continue
			endif

			Wave xyzi = Extend_XYZ_by1(atomRefs[i])
			MatrixOP/FREE xyzi = (direct x xyzi^t )^t		// xyz of atoms of type atomi in real units
		 	Ni = DimSize(xyzi,0)
		
			// find the atom in xyzi that is closest to xyzMid
			MatrixOP/FREE distij = sqrt(sumRows( magSqr(xyzi - rowRepeat(xyzMid,Ni)) ))
			distij = distij > distCurrent ? distij : Inf	// avoid two atoms at less than current distance
			WaveStats/M=1/Q distij
			if (distij[V_minloc] < distMin)					// new smallest distance
				distMin = distij[V_minloc]						// closest distance from an atomj to an atomi
				MatrixOP/FREE sumNeighbors = sum( greater( 1e-2, Abs(distij - rowRepeat(distMin,Ni)) ) )
				multMin = sumNeighbors[0]						// number of neighbors of atomj at distance of distMin
				ljMin = lj
				liMin = li
			endif
		endfor
	endfor

	if (numtype(distMin)==0 && distMin < blen)				// found a valid looking bond, add it
		Variable Nbonds = 	xtal.Nbonds
		xtal.bond[Nbonds].label0 = ljMin
		xtal.bond[Nbonds].label1 = liMin
		xtal.bond[Nbonds].N = 1
		xtal.bond[Nbonds].len[0] = distMin
		xtal.bond[Nbonds].len[4] = multMin					// NOTE, this is a mis-use of N, it must be fixed later
		Nbonds += 1
		xtal.Nbonds = Nbonds
		return 0
	endif
	return 1															// nothing found
End
//
Static Function BondExists(xtal,l0,l1)
	STRUCT crystalStructure &xtal
	String l0, l1

	Variable i
	for (i=0;i<xtal.Nbonds;i+=1)
		if (CmpStr(l0,xtal.bond[i].label0)==0 && CmpStr(l1,xtal.bond[i].label1)==0)
			return 1
		elseif (CmpStr(l0,xtal.bond[i].label1)==0 && CmpStr(l1,xtal.bond[i].label0)==0)
			return 1
		endif
	endfor
	return 0
End



#else

#ifdef USE_BONDS_VERSION_1

Function ShowCalculatedBonds_xtal()
	// Show Calculatee Bonds for current xtal. Only Prints, does NOT change anything.
	STRUCT crystalStructure xtal
	if (FillCrystalStructDefault(xtal))			//fill the lattice structure with test values
		DoAlert 0, "ERROR AnalyzeBonds()\rno lattice structure found"
		return 1
	endif

	Struct ExcessiveBondList FinalBonds
	Variable err = AnalyzeBonds(xtal, FinalBonds)
	if (err)
		print "ERROR -- AnalyzeBonds(xtal, FinalBonds)  failed"
	else
		print "Bonds Calculated from Current xtal (Nothing Changed):"
		PrintExcessiveBondList(FinalBonds)
	endif
End


Function/WAVE MakeBondList_blen(prefix,xyz,[bLen])	// This is a guess when no bonds are given
	// use with MakeBondList_blen()
	String prefix
	Wave xyz				// list of atom xyz positions
	Variable bLen		// maximum distance that gets a bond
	bLen = ParamIsDefault(bLen) || blen<=0 || numtype(blen) ? LatticeSym_maxBondLen : bLen

	STRUCT crystalStructure xtal
	if (FillCrystalStructDefault(xtal))			//fill the lattice structure with test values
		DoAlert 0, "ERROR AnalyzeBonds()\rno lattice structure found"
		return $""
	endif

	if (!WaveExists(xyz))
		return $""
	elseif (DimSize(xyz,0)<2 || DimSize(xyz,1) < xtal.dim)
		return $""
	endif

	Struct ExcessiveBondList FinalBonds
	Variable err = AnalyzeBonds(xtal, FinalBonds)
	if (err)
		print "ERROR -- AnalyzeBonds(FinalBonds)  failed"
		return AtomView#MakeBondList_blenProto(prefix,xyz, bLen=bLen)		// fall back in case of error
	else
		print "FinalBonds are:"
		PrintExcessiveBondList(FinalBonds)
		// transfer FinalBonds into xtal
		Variable i, m, N=FinalBonds.N
		for (i=0,m=0; i<N; i+=1)
			if (FinalBonds.all[i].valid)
				xtal.bond[m].label0 = FinalBonds.all[i].label0
				xtal.bond[m].label1 = FinalBonds.all[i].label1
				xtal.bond[m].N = 1
				xtal.bond[m].len[0] = FinalBonds.all[i].len
				m += 1
			endif
		endfor
		xtal.Nbonds = m					// number of valid bonds
	endif
	return AtomView#MakeBondList_Given(prefix,xtal,xyz)	// This makes the bond list from given bond lengths
End


Static Function AnalyzeBonds(xtal, FinalBonds)
	STRUCT crystalStructure &xtal
	Struct ExcessiveBondList &FinalBonds			// result gets returned in FinalBonds

	Variable Natoms=xtal.N							// number of distinct atom types, both element and wyckoff location
	if (Natoms<1)
		return 1
	endif
	Variable dim = xtal.dim
	Wave direct = directFrom_xtal(xtal)			// a FREE wave with real lattice (a,b,c, are the columns)

	Wave enegWave = root:Packages:Elements:electroneg	// use enegWave[Z]
	Wave atomRad = root:Packages:Elements:atomRadius		// use atomRad[Z] (Å)
	Wave covRad = root:Packages:Elements:covRadius		// use covRad[Z] (Å)

	Variable Nbonds=0
	Struct ExcessiveBondList AllBonds
	initExcessiveBondList(AllBonds)

	String labeli, labelj, str, bondType, labelList=""
	Variable Zj,Zi, i,j
	Variable deltaEneg, ionic, covalent, distCalc, Nj,Ni, distMin, mult, dd
	for (j=0;j<Natoms;j+=1)
		Wave atomj = $("root:Packages:Lattices:atom"+num2istr(j))
		Zj = NumberByKey("Zatom",note(atomj),"=")
		labelj = StringByKey("atomType",note(atomj),"=")
		labelList += labelj + ";"
		Wave xyz0 = Extend_XYZ_by1(atomj)

		MatrixOP/FREE xyz0 = (direct x xyz0^t )^t	// xyz of atoms of type atomj in real units
		Nj = DimSize(xyz0,0)
		// find the centermost atom in xyz0
		MatrixOP/FREE xyzMid = sumCols(xyz0)
		Redimension/N=(dim) xyzMid
		xyzMid /= Nj											// average value
		MatrixOP/FREE dist2 = sumRows( magSqr(xyz0 - rowRepeat(xyzMid,Nj)) )
		WaveStats/M=1/Q dist2
		xyzMid = xyz0[V_minloc][p]						// now xyzMid is the actual atom xyz closest to average
		WaveClear xyz0

		for (i=j;i<Natoms;i+=1)
			Wave atomi = $("root:Packages:Lattices:atom"+num2istr(i))
 			Zi = NumberByKey("Zatom",note(atomi),"=")
			labeli = StringByKey("atomType",note(atomi),"=")

			Wave xyz1 = Extend_XYZ_by1(atomi)
			MatrixOP/FREE xyz1 = (direct x xyz1^t )^t	// xyz of atoms of type atomi in real units
 			Ni = DimSize(xyz1,0)

			// find the atom in xyz1 that is closest to xyzMid
			MatrixOP/FREE distij = sqrt(sumRows( magSqr(xyz1 - rowRepeat(xyzMid,Ni)) ))
			distij = distij > 1e-4 ? distij : Inf	// avoid two atoms at same positione
			WaveStats/M=1/Q distij
			distMin = distij[V_minloc]					// closest distance from an atomj to an atomi
			MatrixOP/FREE neighbors = greater( 1e-2, Abs(distij - rowRepeat(distMin,Ni)) )
			mult = sum(neighbors)							// number of neighbors of atomj at distance of distMin

			deltaEneg = abs(enegWave[Zj] - enegWave[Zi])
			ionic = deltaEneg > 1.5
			covalent = !ionic
			distCalc = covalent ? covRad[Zj] + covRad[Zi] : atomRad[Zj] + atomRad[Zi]
			distCalc /= 10									// convert Å --> nm
			bondType = SelectString(covalent,"ionic","covalent")
			dd = abs(distCalc - distMin)
			if (dd<0.05)										// a bond, within 0.5 Å of predicted
				AllBonds.all[Nbonds].valid = 1
				AllBonds.all[Nbonds].label0 = labelj		;	AllBonds.all[Nbonds].label1 = labeli
				AllBonds.all[Nbonds].Z0 = Zj				;	AllBonds.all[Nbonds].Z1 = Zi
				AllBonds.all[Nbonds].type = bondType		;	AllBonds.all[Nbonds].mult = mult
				AllBonds.all[Nbonds].len = distMin		;	AllBonds.all[Nbonds].calc = distCalc
				Nbonds += 1
				AllBonds.N = Nbonds
			endif
		endfor
	endfor

	// relabel some atom types: covalent --> ionic, an atom cannot be both covalent and ionic
	String l0, ionicList=""
	for(i=0;i<Nbonds;i+=1)
		if (CmpStr(AllBonds.all[i].type,"ionic")==0)
			l0 = AllBonds.all[i].label0
			if (WhichListItem(l0,ionicList)<0)
				ionicList += l0+";"
			endif
			l0 = AllBonds.all[i].label1
			if (WhichListItem(l0,ionicList)<0)
				ionicList += l0+";"
			endif
		endif
	endfor
	for(i=0;i<Nbonds;i+=1)
		if (WhichListItem(AllBonds.all[i].label0,ionicList)>=0 || WhichListItem(AllBonds.all[i].label1,ionicList)>=0)
			AllBonds.all[i].valid = CmpStr(AllBonds.all[i].type,"ionic") ? 0 : 1	// if one of the atoms is ionic, then coavalent is in-valid
		endif
	endfor

	// re-check distCalc since some covalents were changed to ionic
	for(i=0;i<Nbonds;i+=1)
		if ( CmpStr(AllBonds.all[i].type,"ionic") || !AllBonds.all[i].valid )	// only re-check ionic bonds
			continue
		endif
		distCalc = (atomRad[AllBonds.all[i].Z0] + atomRad[AllBonds.all[i].Z1])/10
		dd = abs(distCalc - AllBonds.all[i].len)
		if (dd<0.05)																// a bond, within 0.5 Å of predicted
			AllBonds.all[i].calc = distCalc
		else
			AllBonds.all[i].valid = 0
		endif
	endfor

	initExcessiveBondList(FinalBonds)
	Variable Nfinal=0
	Struct ExcessiveBondList OneLabelBonds
	for (i=0;i<ItemsInList(labelList);i+=1)								// loop over all atom types
		l0 = StringFromList(i,labelList)
		initExcessiveBondList(OneLabelBonds)
		DetermineOneAtomBonds(l0,AllBonds, OneLabelBonds)			// returns OneLabelBonds aound l0
		for (j=0;j<OneLabelBonds.N;j+=1)
			if (WhichBondItem(OneLabelBonds.all[j],FinalBonds)<0)	// bond not present, add it
				FinalBonds.all[Nfinal] = OneLabelBonds.all[j]
				Nfinal += 1
				FinalBonds.N = Nfinal
			endif
		endfor
	endfor

	FinalBonds.unused = ""
	String usedAtoms="", unUsedAtoms=""
	for (i=0;i<FinalBonds.N;i+=1)
		usedAtoms += FinalBonds.all[i].label0 + ";"				// this list may have duplicates
		usedAtoms += FinalBonds.all[i].label1 + ";"
	endfor

	String lab=""
	for (i=0;i<xtal.N;i+=1)
		lab = xtal.atom[i].name
		if (WhichListItem(lab,usedAtoms) < 0)
			unUsedAtoms += lab + ";"
		endif
	endfor
	FinalBonds.unused = unUsedAtoms									// list of un-used atom types
	return 0
End
//
Static Function DetermineOneAtomBonds(l0, AllBonds, OneLabelBonds)
	String l0
	Struct ExcessiveBondList &AllBonds
	Struct ExcessiveBondList &OneLabelBonds
	initExcessiveBondList(OneLabelBonds)

	Variable Nall=AllBonds.N
	Variable i, Ni, mult, isCovalent=1
	for (i=0; i<Nall; i+=1)		// loop over all bonds containing label0
		if (!AllBonds.all[i].valid)
			continue
		endif
		if (CmpStr(AllBonds.all[i].label0,l0) && CmpStr(AllBonds.all[i].label1,l0))	// l0 not persent
			continue
		endif
		isCovalent = CmpStr(AllBonds.all[i].type,"covalent") ? 0 : isCovalent
	endfor

	for (i=0,Ni=0,mult=0; i<Nall; i+=1)		// loop over all bonds containing label0, save the ones I am interested in
		if (!AllBonds.all[i].valid)
			continue
		endif
		if (CmpStr(AllBonds.all[i].label0,l0) && CmpStr(AllBonds.all[i].label1,l0))	// l0 not persent
			continue
		endif
		if (isCovalent || CmpStr(AllBonds.all[i].type,"covalent"))
			OneLabelBonds.all[Ni] = AllBonds.all[i]
			Ni += 1
			OneLabelBonds.N = Ni
			mult += AllBonds.all[i].mult
		endif
		isCovalent = CmpStr(AllBonds.all[i].type,"covalent") ? 0 : isCovalent
	endfor

	// reorder by bond length
	Make/N=(OneLabelBonds.N)/D/FREE lengths = OneLabelBonds.all[p].len, points=p
	Sort lengths, points	// does not sort lengths, only points
	Struct ExcessiveBondList SwapBonds
	SwapBonds = OneLabelBonds
	for (i=0;i<OneLabelBonds.N;i+=1)
		OneLabelBonds.all[points[i]] = SwapBonds.all[i]
	endfor

	if (mult > 10)			// too many bonds to this atom, reduce to <= 10
		Make/N=(OneLabelBonds.N)/I/FREE mults = OneLabelBonds.all[p].mult
		if (numpnts(mults)>1)
		mults[1,OneLabelBonds.N-1] += mults[p-1]
		endif

		Variable imax = BinarySearch(mults,10)		// 10 is the maximum number of bonds for one atom
		imax = imax==-1 ? 0 : imax
		imax = imax==-2 ? OneLabelBonds.N-1 : imax
		OneLabelBonds.N = imax + 1
		for (i=imax+1; i<100; i+=1)
			OneLabelBonds.all[i].valid = 0
		endfor
		mult = mults[imax]
		if (mult>10)
			OneLabelBonds.N = 0
			OneLabelBonds.all[0].valid = 0
			printf "For '%s' <--> '%s', mult = %d, Skipping\r",OneLabelBonds.all[0].label0,OneLabelBonds.all[0].label1,OneLabelBonds.all[0].mult
		endif
	endif
End


Static Structure ExcessiveBondList
	int16 N						// number of bonds in len (often just 1)
	Struct tempBondStructure all[100]
	char unused[400]			// list of unused labels
EndStructure
//
Static Structure tempBondStructure	// defines the type of bond between two atom types
	int16 valid
	char label0[60]			// label for first atom, usually starts with atomic symbol
	char label1[60]			// label for second atom, usually starts with atomic symbol
	int16 Z0, Z1
	char type[10]
	int16 mult
	double len
	double calc
EndStructure
//
Static Function PrintExcessiveBondList(AllBonds)
	Struct ExcessiveBondList &AllBonds
	Struct tempBondStructure bs
	Variable i, dd, Ntrue=0
	for (i=0,Ntrue=0; i<AllBonds.N; i+=1)
		if (AllBonds.all[i].valid)
			Ntrue += 1
		endif
	endfor
	for (i=0; i<AllBonds.N; i+=1)
		PrintTempBondStructure(AllBonds.all[i])
	endfor

	if (ItemsInList(AllBonds.unused)>0)
		printf "The atoms: {%s} are NOT associated with any bonds.\r",TrimEnd(ReplaceString(";",AllBonds.unused,", "),chars=", ")
	else				// print list of those atoms not associated with a bond
		print "    All atom types are associated with at least 1 bond"
	endif
	return 0
End
//
Static Function PrintTempBondStructure(bs)
	Struct tempBondStructure &bs
	if (bs.valid)
		Variable dd = abs(bs.len - bs.calc)
		printf "   [%3s %2d] <--> [%3s %2d] = %.4f nm,  mult=%2d,  '%s'  (Calc=%.3f nm --> ∆=%.4f nm)\r",bs.label0,bs.Z0,bs.label1,bs.Z1,bs.len,bs.mult, bs.type, bs.calc, dd
	else
		printf " empty"
	endif
	return 0
End
//
Static Function initExcessiveBondList(AllBonds)
	Struct ExcessiveBondList &AllBonds
	Struct tempBondStructure bs
	Variable i
	AllBonds.N = 0
	for (i=0;i<100;i+=1)
		AllBonds.all[i].valid = 0
	endfor
	AllBonds.unused = ""
End
//
Static Function WhichBondItem(bs,AllBonds)
	// sort of like WhichListItem(), returns index into AllBonds, or -1 is bs not in AllBonds
	Struct tempBondStructure &bs
	Struct ExcessiveBondList &AllBonds
	Variable i
	for (i=0;i<AllBonds.N;i+=1)
		if (BSequal(bs,AllBonds.all[i]))
			return i
		endif
	endfor
	return -1
End
//
Static Function BSequal(b0,b1)
	Struct tempBondStructure &b0
	Struct tempBondStructure &b1
	if (!b0.valid || !b1.valid)
		return 0						// both must be valid
	endif
	Variable match0 = CmpStr(b0.label0,b1.label0)==0 && CmpStr(b0.label1,b1.label1)==0
	Variable match1 = CmpStr(b0.label0,b1.label1)==0 && CmpStr(b0.label1,b1.label0)==0
	if (!match0 && !match1)
		return 0						// label mis-match
	endif
	if (match0)
		if (b0.Z0 != b1.Z0 || b0.Z1 != b1.Z1)
			return 0					// Z mis-match
		endif
	else		
		if (b0.Z0 != b1.Z1 || b0.Z1 != b1.Z0)
			return 0					// Z mis-match
		endif
	endif
	if (CmpStr(b0.type,b1.type))
		return 0						// type mis-match
	endif
	if (b0.mult != b1.mult)
		return 0						// mult mis-match
	endif
	if (abs(b0.len-b1.len)>1e-6)
		return 0						// length mis-match
	endif
	if (abs(b0.calc-b1.calc)>1e-6)
		return 0						// calculated length mis-match
	endif
	return 1							// everything matches
End

#else
	// neither USE_BONDS_VERSION_1 nor USE_BONDS_VERSION_2 are defined,  defaults to the one in AtomView.ipf
#endif
#endif



Static Function/WAVE Extend_XYZ_by1(xyz0)
	// extend xyz in all directions by ±1
	Wave xyz0					// wave with RELATIVE coordinates in base cell, NOT altered

	Variable dim = (DimSize(xyz0,1) == 2) ? 2 : 3	// dimension (2 or 3)
	if (dim == 3)								// 3D
		Make/D/FREE offsets = { {-1,0,1,-1,0,1,-1,0,1,-1,0,1,-1,0,1,-1,0,1,-1,0,1,-1,0,1,-1,0,1},{-1,-1,-1,0,0,0,1,1,1,-1,-1,-1,0,0,0,1,1,1,-1,-1,-1,0,0,0,1,1,1},{-1,-1,-1,-1,-1,-1,-1,-1,-1,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1} }
	elseif (dim == 2)							// 2D
		Make/D/FREE offsets = { {-1,0,1,-1,0,1,-1,0,1}, {-1,-1,-1,0,0,0,1,1,1} }
	else
		return $""
	endif
	Variable Noff = DimSize(offsets,0)	// number of offsets

	Variable N0=DimSize(xyz0,0), Nmax=Noff*N0, i, j, i1
	Make/N=(Nmax,dim)/D/FREE xyz
	for (i=0;i<Nmax;i+=N0)
		i1 = i + (N0-1)
		j = floor(i/N0)
		xyz[i,i1][] = xyz0[p-i][q] + offsets[j][q]
	endfor
	return xyz
End

