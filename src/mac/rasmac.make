#   File:       RasMol.symantec.make
#   Target:     RasMol.symantec
#   Sources:
#               abstree.c applemac.c command.c molecule.c outfile.c pixutils.c  
#		rasmac.c render.c script.c transfor.c 
#   Created:    Tue Oct 10 15:08:59 1995


PPC_C = SMrC
PPC_OPTIONS = 
OPTIONS = -ansi relaxed -align power -opt all
C_OPTIONS = {OPTIONS} 



COBJECTS = abstree.c.o applemac.c.o command.c.o infile.c.o �
           molecule.c.o outfile.c.o pixutils.c.o rasmac.c.o �
		   render.c.o repres.c.o script.c.o transfor.c.o 

RasMol.symantec � RasMol.symantec.make RasMac.rsrc {COBJECTS} 
	PPCLink   -warn  {COBJECTS} �
		"{SharedLibraries}"InterfaceLib �
		"{SharedLibraries}"MathLib  �
		"{SharedLibraries}"StdCLib  �
		"{PPCLibraries}"StdCRuntime.o  �
		"{PPCLibraries}"PPCCRuntime.o �
		-o RasMol.symantec
	Echo "include �"RasMac.rsrc�" ;" �
		| Rez -a  -o RasMol.symantec
	SetFile -c RSML RasMol.symantec


abstree.c.o � abstree.c
	{PPC_C} abstree.c {PPC_OPTIONS} {C_OPTIONS}

applemac.c.o � applemac.c
	{PPC_C} applemac.c {PPC_OPTIONS} {C_OPTIONS}

command.c.o � command.c
	{PPC_C} command.c {PPC_OPTIONS} {C_OPTIONS}

infile.c.o � infile.c
	{PPC_C} infile.c {PPC_OPTIONS} {C_OPTIONS}

molecule.c.o � molecule.c
	{PPC_C} molecule.c {PPC_OPTIONS} {C_OPTIONS}

outfile.c.o �  outfile.c
	{PPC_C} outfile.c {PPC_OPTIONS} {C_OPTIONS}

pixutils.c.o �  pixutils.c
	{PPC_C} pixutils.c {PPC_OPTIONS} {C_OPTIONS}

rasmac.c.o �  rasmac.c
	{PPC_C} rasmac.c {PPC_OPTIONS} {C_OPTIONS}

render.c.o �  render.c
	{PPC_C} render.c {PPC_OPTIONS} {C_OPTIONS}

repres.c.o �  repres.c
	{PPC_C} repres.c {PPC_OPTIONS} {C_OPTIONS}

script.c.o �  script.c
	{PPC_C} script.c {PPC_OPTIONS} {C_OPTIONS}

transfor.c.o �  transfor.c
	{PPC_C} transfor.c {PPC_OPTIONS} {C_OPTIONS}

