###########################################################################
#                            RasMol 2.7.1.1                               # 
#                                                                         #
#                                 RasMol                                  #
#                 Molecular Graphics Visualisation Tool                   #
#                            17 January 2001                              #
#                                                                         #
#                   Based on RasMol 2.6 by Roger Sayle                    #
# Biomolecular Structures Group, Glaxo Wellcome Research & Development,   #
#                      Stevenage, Hertfordshire, UK                       #
#         Version 2.6, August 1995, Version 2.6.4, December 1998          #
#                   Copyright (C) Roger Sayle 1992-1999                   #
#                                                                         #
#               Version 2.7.0 Mods by Herbert J. Bernstein                #
#           Bernstein + Sons, P.O. Box 177, Bellport, NY, USA             #
#                      yaya@bernstein-plus-sons.com                       #
#                               March 1999                                #
#              Copyright (C) Herbert J. Bernstein 1998-2001               #
#                                                                         #
# Please read the file NOTICE for important notices which apply to this   #
# package. If you are not going to make changes to RasMol, you are not    # 
# only permitted to freely make copies and distribute them, you are       # 
# encouraged to do so, provided you do the following:                     #
#   * 1. Either include the complete documentation, especially the file   #
#     NOTICE, with what you distribute or provide a clear indication      #
#     where people can get a copy of the documentation; and               #
#   * 2. Please give credit where credit is due citing the version and    # 
#     original authors properly; and                                      # 
#   * 3. Please do not give anyone the impression that the original       # 
#     authors are providing a warranty of any kind.                       #
#                                                                         #
# If you would like to use major pieces of RasMol in some other program,  #
# make modifications to RasMol, or in some other way make what a lawyer   #
# would call a "derived work", you are not only permitted to do so, you   #
# are encouraged to do so. In addition to the things we discussed above,  #
# please do the following:                                                #
#   * 4. Please explain in your documentation how what you did differs    #
#     from this version of RasMol; and                                    #
#   * 5. Please make your modified source code available.                 #
#                                                                         #
# This version of RasMol is not in the public domain, but it is given     #
# freely to the community in the hopes of advancing science. If you make  #
# changes, please make them in a responsible manner, and please offer us  #
# the opportunity to include those changes in future versions of RasMol.  #
###########################################################################
#   File:       RasMol.symantec.make
#   Target:     RasMol.symantec
#   Sources:
#      abstree.c cmndline.c command.c infile.c molecule.c outfile.c
#      pixutils.c render.c repres.c script.c transfor.c tokens.c
#      cif_fract.c cif.c cif_ctonum.c cif_stx.c
#   Created:    Monday March 29, 1999


PPC_C = SMrC
PPC_OPTIONS = 
OPTIONS = -ansi relaxed -align power -opt all
C_OPTIONS = {OPTIONS} 



COBJECTS = abstree.c.o applemac.c.o cmndline.c.o command.c.o infile.c.o ¶
           molecule.c.o outfile.c.o pixutils.c.o rasmac.c.o ¶
		   render.c.o repres.c.o script.c.o transfor.c.o tokens.c.o¶
		   cif_fract.c.o cif.c.o cif_ctonum.c.o cif_stx.c.o langsel.c.o

RasMol.symantec Ä RasMol.symantec.make RasMac.rsrc {COBJECTS} 
	PPCLink   -warn  {COBJECTS} ¶
		"{SharedLibraries}"InterfaceLib ¶
		"{SharedLibraries}"MathLib  ¶
		"{SharedLibraries}"StdCLib  ¶
		"{PPCLibraries}"StdCRuntime.o  ¶
		"{PPCLibraries}"PPCCRuntime.o ¶
		-o RasMol.symantec
	Echo "include ¶"RasMac.rsrc¶" ;" ¶
		| Rez -a  -o RasMol.symantec
	SetFile -c RSML RasMol.symantec


abstree.c.o Ä abstree.c
	{PPC_C} abstree.c {PPC_OPTIONS} {C_OPTIONS}

applemac.c.o Ä applemac.c
	{PPC_C} applemac.c {PPC_OPTIONS} {C_OPTIONS}

cmndline.c.o Ä cmndline.c
	{PPC_C} cmndline.c {PPC_OPTIONS} {C_OPTIONS}

command.c.o Ä command.c
	{PPC_C} command.c {PPC_OPTIONS} {C_OPTIONS}

infile.c.o Ä infile.c
	{PPC_C} infile.c {PPC_OPTIONS} {C_OPTIONS}

molecule.c.o Ä molecule.c
	{PPC_C} molecule.c {PPC_OPTIONS} {C_OPTIONS}

outfile.c.o Ä  outfile.c
	{PPC_C} outfile.c {PPC_OPTIONS} {C_OPTIONS}

pixutils.c.o Ä  pixutils.c
	{PPC_C} pixutils.c {PPC_OPTIONS} {C_OPTIONS}

rasmac.c.o Ä  rasmac.c
	{PPC_C} rasmac.c {PPC_OPTIONS} {C_OPTIONS}

render.c.o Ä  render.c
	{PPC_C} render.c {PPC_OPTIONS} {C_OPTIONS}

repres.c.o Ä  repres.c
	{PPC_C} repres.c {PPC_OPTIONS} {C_OPTIONS}

script.c.o Ä  script.c
	{PPC_C} script.c {PPC_OPTIONS} {C_OPTIONS}

transfor.c.o Ä  transfor.c
	{PPC_C} transfor.c {PPC_OPTIONS} {C_OPTIONS}

trokens.c.o Ä  tokens.c
	{PPC_C} tokens.c {PPC_OPTIONS} {C_OPTIONS}

cif_fract.c.o Ä  cif_fract.c
	{PPC_C} cif_fract.c {PPC_OPTIONS} {C_OPTIONS}

cif.c.o Ä  cif.c
	{PPC_C} cif.c {PPC_OPTIONS} {C_OPTIONS}

cif_ctonum.c.o Ä  cif_ctonum.c
	{PPC_C} cif_ctonum.c {PPC_OPTIONS} {C_OPTIONS}

cif_stx.c.o Ä  cif_stx.c
	{PPC_C} cif_stx.c {PPC_OPTIONS} {C_OPTIONS}

langsel.c.o Ä  langsel.c
	{PPC_C} langsel.c {PPC_OPTIONS} {C_OPTIONS}

