#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#pragma version = 3.2

Static Function IgorStartOrNewHook(IgorApplicationNameStr)
	String IgorApplicationNameStr

	Variable g = NumVarOrDefault("root:Packages:MICRO_GEOMETRY_VERSION",NaN)
	if (numtype(g) && !DataFolderExists("root:Packages:"))
		Execute/P/Q "LaueGoFirst#setMICRO_GEOMETRY_VERSION_PATH()"
		Execute/P "INSERTINCLUDE  \"mdaFiles\", version>=1.08";Execute/P "COMPILEPROCEDURES "
	endif
End
