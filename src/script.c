/* script.c
 * RasMol2 Molecular Graphics
 * Roger Sayle, October 1994
 * Version 2.5
 */
#define SCRIPT
#include "rasmol.h"

#ifdef IBMPC
#include <windows.h>
#include <malloc.h>
#endif
#ifdef APPLEMAC
#ifdef __CONDITIONALMACROS__
#include <Printing.h>
#else
#include <PrintTraps.h>
#endif
#include <Types.h>
#endif
#ifndef sun386
#include <stdlib.h>
#endif
 

#include <stdio.h>
#include <ctype.h>
#include <math.h>

#include "script.h"
#include "molecule.h"
#include "command.h"
#include "abstree.h"
#include "transfor.h"
#include "render.h"
#include "graphics.h"
#include "pixutils.h"

#ifdef INVERT
#define InvertY(y) (y)
#else
#define InvertY(y) (-(y))
#endif

#define Round(x)       ((int)(x))
#define DatWirFlag  (Long)0x10000
#define DatCylFlag  (Long)0x20000

/* Special thanks to Profs Jane and David Richardson
 * for the following Kinemage colour lookup table */

#define MAXMAGECOL 20
static struct {
        int r, g, b;
        char *name;
        } KinemageCol[MAXMAGECOL] = {
    { 255,   0,   0, "red"        },  /*  1 */
    {  23, 256,   0, "green"      },  /*  2 */
    {  62,  62, 255, "blue"       },  /*  3 */
    {   0, 242, 226, "cyan"       },  /*  4 */
    { 255, 246,   0, "yellow"     },  /*  5 */
    { 255,   0, 234, "magenta"    },  /*  6 */
    { 255, 255, 255, "white"      },  /*  7 */
    { 255, 101, 117, "pink"       },  /*  8 */
    { 255,  93,   0, "orange"     },  /*  9 */
    { 140,  54, 255, "purple"     },  /* 10 */
    {  58, 144, 255, "skyblue"    },  /* 11 */
    { 175, 117,  89, "brown"      },  /* 12 */
    { 125, 125, 125, "gray"       },  /* 13 */
    { 255, 156,   0, "gold"       },  /* 15 */
    { 246, 246, 117, "yellowtint" },  /* 16 */
    {   0, 250, 109, "seagreen"   },  /* 17 */
    { 255, 171, 187, "pinktint"   },  /* 18 */
    { 175, 214, 255, "bluetint"   },  /* 19 */
    { 152, 255, 179, "greentint"  },  /* 20 */
    { 255,   0, 101, "hotpink"    }   /* 21 */
        };


typedef struct {
        Long datum;
        Long count;
    } FreqEntry;

#define FREQSIZE  8
static FreqEntry Freq[FREQSIZE];


static Atom __far *MagePrev;
static char *MageCol;
static FILE *OutFile;
static int SelectAll;

/* Macros for commonly used loops */
#define ForEachAtom  for(chain=Database->clist;chain;chain=chain->cnext) \
                     for(group=chain->glist;group;group=group->gnext)    \
                     for(aptr=group->alist;aptr;aptr=aptr->anext)
#define ForEachBond  for(bptr=Database->blist;bptr;bptr=bptr->bnext)
#define ForEachBack  for(chain=Database->clist;chain;chain=chain->cnext) \
                     for(bptr=chain->blist;bptr;bptr=bptr->bnext)


static void FatalScriptError( ptr )
    char *ptr;
{
    if( CommandActive ) WriteChar('\n');
    WriteString("Script Error: Unable to create file `");
    WriteString( ptr );  WriteString("'!\n");
    CommandActive = False;
}


#ifdef FUNCPROTO
static void IncFreqTable( Long );
static Long GetBondDatum( Bond __far* );
static Long GetHBondDatum( HBond __far* );
#endif

static void ResetFreqTable()
{
    register int i;

    for( i=0; i<FREQSIZE; i++ )
        Freq[i].count = 0;
}

static void IncFreqTable( datum )
    Long datum;
{
    register Long count;
    register int i;

    for( i=0; i<FREQSIZE; i++ )
        if( !Freq[i].count )
        {   Freq[i].datum = datum;
            Freq[i].count = 1;
            return;
        } else if( Freq[i].datum == datum )
        {   count = Freq[i].count+1;
            while( i && (Freq[i-1].count<=count) )
            {   Freq[i] = Freq[i-1];  
                i--;
            }
            Freq[i].datum = datum;
            Freq[i].count = count;
            return;
        }

    /* Replace Singletons! */
    if( Freq[FREQSIZE-1].count == 1 )
        Freq[FREQSIZE-1].datum = datum;
}


static Long GetBondDatum( bptr )
    Bond __far *bptr;
{
    if( bptr->flag & CylinderFlag )
    {   return( DatCylFlag | bptr->radius );
    } else if( bptr->flag & WireFlag )
    {   return( DatWirFlag );
    } else return( (Long)0 );
}

static Long GetHBondDatum( bptr )
    HBond __far *bptr;
{
    if( bptr->flag & CylinderFlag )
    {   return( DatCylFlag | bptr->radius );
    } else if( bptr->flag & WireFlag )
    {   return( DatWirFlag );
    } else return( (Long)0 );
}

 
#ifndef SCRIPT
#ifdef FUNCPROTO
static void WriteMolScriptAtomSel( Chain __far*, Group __far*, Atom __far* );
#endif

static void WriteMolScriptColour( r, g, b )
    int r, g, b;
{
    fprintf(OutFile," rgb %#g %#g %#g",r/255.0,g/255.0,b/255.0);
}

static void WriteMolScriptAtomSel( chain, group, aptr )
    Chain __far *chain;
    Group __far *group;
    Atom __far *aptr;
{
    register char *ptr;
    register int i;

    fputs("require atom ",OutFile);
    ptr = ElemDesc[aptr->refno];
    for( i=0; i<4; i++ )
        if( ptr[i]=='*' )
        {   fputc('\'',OutFile);
        } else if( ptr[i]!=' ' )
            fputc(ptr[i],OutFile);

    fputs(" and in residue ",OutFile);
    if( chain->ident!=' ' && !isdigit(chain->ident) )
        fputc(chain->ident,OutFile);
    fprintf(OutFile,"%d",group->serno);
}

static void WriteMolScriptAtoms()
{
    register Chain __far *chain;
    register Group __far *group;
    register Atom __far *aptr;
    register ShadeDesc *shade;
    register Label *label;
    register char *ptr;

    char buffer[80];


    ForEachAtom
    {   if( aptr->flag & SphereFlag )
        {   /* Atom Colour */
            fputs("set atomcolour ",OutFile);
            WriteMolScriptAtomSel(chain,group,aptr);
            shade = Shade + Colour2Shade(aptr->col);
            WriteMolScriptColour(shade->r,shade->g,shade->b);
            fputs(";\n",OutFile);

            /* CPK Sphere */
            fputs("set atomradius ",OutFile);
            WriteMolScriptAtomSel(chain,group,aptr);
            fprintf(OutFile," %#g;\ncpk ",aptr->radius/250.0);
            WriteMolScriptAtomSel(chain,group,aptr);
            fputs(";\n",OutFile);
        }

        if( aptr->label )
        {   /* Atom Label */
            label = (Label*)aptr->label;
            FormatLabel(chain,group,aptr,label->label,buffer);

            fputs("label ",OutFile);
            WriteMolScriptAtomSel(chain,group,aptr);
            fputs(" \"",OutFile);
            for( ptr=buffer; *ptr; ptr++ )
                if( *ptr!='%' ) fputc(*ptr,OutFile);
            fputs("\";\n",OutFile);
        }
    }
}
#endif

static void MolScriptSegment( ptr, src, dst, chain )
    char *ptr;  int src, dst;  char chain;
{   
    if( (chain!=' ') && !isdigit(chain) ) 
    {   fprintf(OutFile,"  %s from %c%d to %c%d;\n",ptr,chain,src,chain,dst);
    } else fprintf(OutFile,"  %s from %d to %d;\n",ptr,src,dst);
}


int WriteMolScriptFile( name )
    char *name;
{
    register Real temp;
    register Real psi, phi, theta;
    register Chain __far *chain;
    register Group __far *group;
    register Group __far *next;
    register Group __far *prev;
    register int flag,len;
    register char *ptr;

    if( !Database )
        return(False);

    OutFile = fopen(name,"w");
    if( !OutFile )
    {   FatalScriptError(name);
        return(False);
    }
    fprintf(OutFile,"! File: %s\n",name);
    fputs("! Creator: RasMol Version 2.5\n",OutFile);
    fputs("! Version: MolScript v1.3\n\n",OutFile);

    fputs("plot\n",OutFile);
    if( BackR || BackG || BackB )
    {   fputs("  background rgb ",OutFile);
        fprintf(OutFile,"%#g ",  BackR/255.0);
        fprintf(OutFile,"%#g ",  BackG/255.0);
        fprintf(OutFile,"%#g;\n",BackB/255.0);
    }
    temp = 0.004/Scale;
    fprintf(OutFile,"  window %g;\n",temp*Range);
    if( UseSlabPlane )
        fprintf(OutFile,"  slab %g;\n",SideLen/250.0);
    fputc('\n',OutFile);

    fprintf(OutFile,"  read mol \"%s\";\n",InfoFileName);
    fputs("  transform atom *\n",OutFile);
    fputs("    by centre position atom *\n",OutFile);
    fputs("    by rotation x 180.0",OutFile);

    phi = Rad2Deg*asin(RotX[2]);
    if( (int)phi == 90 )
    {   theta = -Rad2Deg*atan2(RotY[0],RotY[1]);
        psi = 0;
    } else if( (int)phi == -90 )
    {   theta = Rad2Deg*atan2(RotY[0],RotY[1]);
        psi = 0;
    } else /* General Case! */
    {   theta = Rad2Deg*atan2(RotY[2],RotZ[2]);
        psi =  -Rad2Deg*atan2(RotX[1],RotX[0]);
    }

    if( (int)psi )   fprintf(OutFile,"\n    by rotation z %#g",InvertY(psi));
    if( (int)phi )   fprintf(OutFile,"\n    by rotation y %#g",phi);
    if( (int)theta ) fprintf(OutFile,"\n    by rotation x %#g",InvertY(-theta));

    if( UseSlabPlane || (XOffset!=WRange) || (YOffset!=HRange) )
    {   fputs("\n    by translation ",OutFile);
        fprintf(OutFile,"%#g ",(XOffset-WRange)*temp);
        fprintf(OutFile,"%#g ",-(YOffset-HRange)*temp);
        if( UseSlabPlane )
        {   temp = (1.0-DialValue[7])/500.0;
            fprintf(OutFile,"%#g",SideLen*temp);
        } else fputs("0.0",OutFile);
    }
    fputs(";\n\n",OutFile);

    /* fputs("  trace amino-acids;\n",OutFile); */

    if( Database->clist )
    {   if( InfoHelixCount<0 )
            DetermineStructure();

        for( chain=Database->clist; chain; chain=chain->cnext )
        {   prev = (Group __far*)0;
            for( group=chain->glist; group && group->gnext; group=next )
            {   next = group->gnext;
                if( next->serno < group->serno )
                {   if( prev && prev!=group )
                        MolScriptSegment("coil",prev->serno,group->serno,
                                                chain->ident);
                    prev = (Group __far*)0;
                    continue;
                }
                flag = group->struc & next->struc;

                if( flag&HelixFlag )
                {   flag = HelixFlag;
                    ptr = "helix";
                } else if( flag&SheetFlag )
                {   flag = SheetFlag;
                    ptr = "strand";
                } else 
                {   if( flag&TurnFlag )
                    {   fputs("  turn residue ",OutFile);
                        if( chain->ident != ' ' )
                            fputc(chain->ident,OutFile);
                        fprintf(OutFile,"%d;\n",group->serno);
                    }
                    if( !prev ) prev = group;
                    continue;
                }

                len = 2;  /* Determine Structure Length */
                while( next->gnext && (next->gnext->struc&flag)
                           && (next->serno<=next->gnext->serno) )
                {   next = next->gnext;
                    len++;
                }

                if( len>2 )
                {   if( prev ) /* MolScript coil or turn? */
                       MolScriptSegment("coil",prev->serno,group->serno,
                                              chain->ident);
                    MolScriptSegment(ptr,group->serno,next->serno,
                                         chain->ident);
	            prev = next;
                } 
            }

            if( prev && prev!=group )  /* C-terminal coil/turn */
                MolScriptSegment("coil",prev->serno,group->serno,chain->ident);
        }
    }

#ifndef SCRIPT
    WriteMolScriptAtoms();
#endif

    fputs("end_plot\n",OutFile);
    fclose(OutFile);
#ifdef APPLEMAC
    SetFileInfo(name,'ttxt','TEXT',133);
#endif
    return( True );
}




#ifdef FUNCPROTO
static void WriteScriptDatum( char*, Long );
static void WriteScriptHBonds( char*, HBond __far* );
#endif

static void WriteScriptAll()
{
    if( !SelectAll )
    {   fputs("select all\n",OutFile);
        SelectAll = True;
    }
}

static void WriteScriptColour( ptr, col )
    char *ptr;  int col;
{
    register ShadeDesc *shade;
    
    if( col )
    {   shade = Shade + Colour2Shade(col);
        fprintf(OutFile,"colour %s [%d,%d,%d]\n",ptr,
                shade->r,shade->g,shade->b);
    } else fprintf(OutFile,"colour %s none\n",ptr);
}


static void WriteScriptBetween( lo, hi )
    int lo, hi;
{
    if( lo != hi )
    {   fprintf(OutFile,"select (atomno>=%d) and (atomno<=%d)\n",lo,hi);
    } else fprintf(OutFile,"select atomno=%d\n",lo);
    SelectAll = False;
}

static void WriteScriptAtoms()
{
    register Chain __far *chain;
    register Group __far *group;
    register Atom __far *aptr;
    register int first,last;
    register int same,init;
    register int cpk,vdw;
    register int col,rad;

    fputs("\n# Atoms\n",OutFile);

    same = True;
    init = False;
    ForEachAtom
        if( !init )
        {   first = last = aptr->serno;
            cpk = IsCPKColour( aptr );
            col = aptr->col;
            init = True;
        } else if( cpk && IsCPKColour(aptr) )
        {   last = aptr->serno;
            if( aptr->col != col )
                col = 0;
        } else if( aptr->col == col )
        {   last = aptr->serno;
            cpk = False;
        } else if( aptr->col != col )
        {   WriteScriptBetween( first, last );
            if( !col )
            {   fputs("colour atoms cpk\n",OutFile);
            } else WriteScriptColour("atoms",col);
                
            first = last = aptr->serno;
            cpk = IsCPKColour( aptr );
            col = aptr->col;
            same = False;
        } else last = aptr->serno; 
        
    if( init )
    {   if( !same )
        {   WriteScriptBetween(first,last);
        } else WriteScriptAll();

        if( !col )
        {   fputs("colour atoms cpk\n",OutFile);
        } else WriteScriptColour("atoms",col);
    }

    if( DrawAtoms )
    {   same = True;
        init = False;
        ForEachAtom
            if( !init )
            {   rad = aptr->flag&SphereFlag? aptr->radius : 0;
                first = last = aptr->serno;
                vdw = IsVDWRadius( aptr );
                init = True;
            } else if( rad == ((aptr->flag&SphereFlag)? aptr->radius : 0) )
            {   if( vdw ) vdw = IsVDWRadius( aptr );
                last = aptr->serno;
            } else if( vdw && IsVDWRadius(aptr) )
            {   last = aptr->serno;
                rad = -1;
            } else 
            {   WriteScriptBetween(first,last);
                if( rad == -1 )
                {   fputs("spacefill on\n",OutFile);
                } else if( rad )
                {   fprintf(OutFile,"spacefill %d\n",rad);
                } else fputs("spacefill off\n",OutFile); 

                rad = aptr->flag&SphereFlag? aptr->radius : 0;
                first = last = aptr->serno;
                vdw = IsVDWRadius( aptr );
                same = False;
            }

        if( !same )
        {   WriteScriptBetween(first,last);
        } else WriteScriptAll();

        if( rad == -1 )
        {   fputs("spacefill on\n",OutFile);
        } else if( rad )
        {   fprintf(OutFile,"spacefill %d\n",rad);
        } else fputs("spacefill off\n",OutFile); 

        if( UseShadow )
        {   fputs("set shadow on\n",OutFile);
        } else fputs("set shadow off\n",OutFile);

    } else
    {   WriteScriptAll();
        fputs("spacefill off\n",OutFile);
    }
        
}


static void WriteScriptDatum( ptr, datum )
    char *ptr;  Long datum;
{
    if( datum & DatCylFlag )
    {   fprintf(OutFile,"%s %d\n",ptr,(int)(datum-DatCylFlag));
    } else if( datum & DatWirFlag )
    {   fprintf(OutFile,"%s on\n",ptr);
    } else fprintf(OutFile,"%s off\n",ptr);
}


static void WriteScriptBonds()
{
    register Bond __far *bptr;
    register Long defdat;
    register Long datum;
    register int col;

    fputs("\n# Bonds\n",OutFile);

    ResetFreqTable();
    ForEachBond
        IncFreqTable(GetBondDatum(bptr));

    WriteScriptAll();
    defdat = Freq[0].datum;
    WriteScriptDatum("wireframe",defdat);

    if( Freq[1].count )
    {   ForEachBond
        {   datum = GetBondDatum(bptr);
            if( datum != defdat )
            {    fprintf(OutFile,"select (atomno=%d) or (atomno=%d)\n",
                         bptr->srcatom->serno, bptr->dstatom->serno );
                 WriteScriptDatum("wireframe",datum);
                 SelectAll = False;
            }
        }
    } else if( !defdat )
        return;

    ResetFreqTable();
    ForEachBond
        IncFreqTable(bptr->col);

    WriteScriptAll();
    if( (col=(int)Freq[0].datum) )
        WriteScriptColour("bonds",col);

    if( Freq[1].count )
        ForEachBond
            if( bptr->col != col )
            {   fprintf(OutFile,"select (atomno=%d) or (atomno=%d)\n",
                    bptr->srcatom->serno, bptr->dstatom->serno );
                WriteScriptColour("bonds",bptr->col);
                SelectAll = False;
            }
}

static void WriteScriptBackbone()
{
    register Chain __far *chain;
    register Bond __far *bptr;

    register Long defdat;
    register Long datum;
    register int col;

    fputs("\n# Backbone\n",OutFile);

    ResetFreqTable();
    ForEachBack
        IncFreqTable(GetBondDatum(bptr));

    WriteScriptAll();
    defdat = Freq[0].datum;
    WriteScriptDatum("backbone",defdat);

    if( Freq[1].count )
    {   ForEachBack
        {   datum = GetBondDatum(bptr);
            if( datum != defdat )
            {    fprintf(OutFile,"select (atomno=%d) or (atomno=%d)\n",
                         bptr->srcatom->serno, bptr->dstatom->serno );
                 WriteScriptDatum("backbone",datum);
                 SelectAll = False;
            }
        }
    } else if( !defdat )
        return;

    ResetFreqTable();
    ForEachBack
        IncFreqTable(bptr->col);

    WriteScriptAll();
    if( (col=(int)Freq[0].datum) )
        WriteScriptColour("backbone",col);

    if( Freq[1].count )
        ForEachBack
            if( bptr->col != col )
            {   fprintf(OutFile,"select (atomno=%d) or (atomno=%d)\n",
                    bptr->srcatom->serno, bptr->dstatom->serno );
                WriteScriptColour("backbone",bptr->col);
                SelectAll = False;
            }
}

static void WriteScriptRibbons()
{
    register Chain __far *chain;
    register Group __far *group;
    register Atom __far *aptr;

    fputs("\n# Ribbons\n",OutFile);

    if( DrawRibbon )
    {   fprintf(OutFile,"set strands %d\n",SplineCount);

        for( chain=Database->clist; chain; chain=chain->cnext )
            for( group=chain->glist; group; group=group->gnext )
            {    if( IsAmino(group->refno) )
                 {   aptr = FindGroupAtom(group,1);
                 } else aptr = FindGroupAtom(group,7);
                 if( !aptr ) continue;

                 fprintf(OutFile,"select atomno=%d\n",aptr->serno);
                 SelectAll = False;

                 if( group->flag & RibbonFlag )
                 {   fprintf(OutFile,"ribbons %d\n",group->width);
                 } else if( group->flag & StrandFlag )
                 {   fprintf(OutFile,"strands %d\n",group->width);
                 } else fputs("ribbons off\n",OutFile);

                 if( group->col1 != group->col2 )
                 {   WriteScriptColour("ribbon1",group->col1);
                     WriteScriptColour("ribbon2",group->col2);
                 } else WriteScriptColour("ribbons",group->col1);
            }
    } else
    {   WriteScriptAll();
        fputs("ribbons off\n",OutFile);
    }
}


static void WriteScriptHBonds( obj, list )
    char *obj;  HBond __far *list;
{
    register HBond __far *ptr;
    register Long defdat;
    register Long datum;
    register int col;

    ResetFreqTable();
    for( ptr=list; ptr; ptr=ptr->hnext )
        IncFreqTable(GetHBondDatum(ptr));

    WriteScriptAll();
    defdat = Freq[0].datum;
    WriteScriptDatum(obj,defdat);

    if( Freq[1].count )
    {   for( ptr=list; ptr; ptr=ptr->hnext )
        {   datum = GetHBondDatum(ptr);
            if( datum != defdat )
            {    fprintf(OutFile,"select (atomno=%d) or (atomno=%d)\n",
                         ptr->src->serno, ptr->dst->serno );
                 WriteScriptDatum(obj,datum);
                 SelectAll = False;
            }
        }
    } else if( !defdat )
        return;

    ResetFreqTable();
    for( ptr=list; ptr; ptr=ptr->hnext )
        IncFreqTable(ptr->col);

    WriteScriptAll();
    if( (col=(int)Freq[0].datum) )
        WriteScriptColour(obj,col);

    if( Freq[1].count )
        for( ptr=list; ptr; ptr=ptr->hnext )
            if( ptr->col != col )
            {   fprintf(OutFile,"select (atomno=%d) or (atomno=%d)\n",
                    ptr->src->serno, ptr->dst->serno );
                WriteScriptColour(obj,ptr->col);
                SelectAll = False;
            }
}

static void WriteScriptLabels()
{
    register Chain __far *chain;
    register Group __far *group;
    register Atom __far *aptr;
    register int first,last;
    register Label *label;

    fputs("\n# Labels\n",OutFile);
    WriteScriptAll();
    fputs("labels off\n",OutFile);
    if( !DrawLabels ) return;

    if( UseLabelCol )
    {   fprintf(OutFile,"colour labels [%d,%d,%d]\n",LabR,LabG,LabB);
    } else fputs("colour labels none\n",OutFile);
    fprintf(OutFile,"set fontsize %d\n",FontSize);

    label = (Label*)0;
    ForEachAtom
        if( aptr->label != label )
        {   if( label )
            {   WriteScriptBetween(first,last);
                fprintf(OutFile,"label \"%s\"\n",label->label);
            }
            label = (Label*)aptr->label;
            first = last = aptr->serno;
        } else last = aptr->serno;

    if( label )
    {   WriteScriptBetween(first,last);
        fprintf(OutFile,"label \"%s\"",label->label);
    }
}


int WriteScriptFile( name )
    char *name;
{
    register int theta,phi,psi;
    register char *ptr;
    register int temp;

    OutFile = fopen(name,"w");
    if( !OutFile )
    {   FatalScriptError(name);
        return(False);
    }

    fprintf(OutFile,"#!rasmol -script\n# File: %s\n",name);
    fputs("# Creator: RasMol Version 2.5\n\n",OutFile);
    fprintf(OutFile,"background [%d,%d,%d]\n",BackR,BackG,BackB);
    if( !Database )
    {   /* No Molecule! */
        fputs("zap\n",OutFile);
        fclose(OutFile);
#ifdef APPLEMAC
        SetFileInfo(name,'RSML','RSML',133);
#endif
        return(True);
    }

    /* Molecule File Name */
    switch( DataFileFormat )
    {   default:
        case(FormatPDB):      ptr = "pdb";      break;
        case(FormatAlchemy):  ptr = "alchemy";  break;
        case(FormatCharmm):   ptr = "charmm";   break;
        case(FormatMol2):     ptr = "mol2";     break;
        case(FormatMDL):      ptr = "mdl";      break;
        case(FormatXYZ):      ptr = "xyz";      break;
    }
    fprintf(OutFile,"load %s \"%s\"\n",ptr,InfoFileName);

    /* Colour Details */
    fprintf(OutFile,"set ambient %d\n", (int)(100*Ambient) );
    fputs("set specular ",OutFile);
    if( FakeSpecular )
    {   fprintf(OutFile,"on\nset specpower %d\n",SpecPower);
    } else fputs("off\n",OutFile);
    putc('\n',OutFile);

    /* Transformation */
    fputs("reset\n",OutFile);
    if( UseSlabPlane )
    {   if( (temp = (int)(50.0*DialValue[7])) )
        {   fprintf(OutFile,"slab %d\n",temp+50);
        } else fputs("slab on\n",OutFile);

        fputs("set slabmode ",OutFile);
        switch( SlabMode )
        {   default:            
            case(SlabClose):    ptr = "solid";    break;
            case(SlabReject):   ptr = "reject";   break;
            case(SlabHalf):     ptr = "half";     break;
            case(SlabHollow):   ptr = "hollow";   break;
            case(SlabSection):  ptr = "section";
        }
        fputs(ptr,OutFile);
        putc('\n',OutFile);
    } else fputs("slab off\n",OutFile);

    phi = Round(Rad2Deg*asin(RotX[2]));
    if( phi == 90 )
    {   theta = -Round(Rad2Deg*atan2(RotY[0],RotY[1]));
        psi = 0;
    } else if( phi == -90 )
    {   theta = Round(Rad2Deg*atan2(RotY[0],RotY[1]));
        psi = 0;
    } else /* General Case! */
    {   theta = Round(Rad2Deg*atan2(RotY[2],RotZ[2]));
        psi =  Round(-Rad2Deg*atan2(RotX[1],RotX[0]));
    }

    if( psi )   fprintf(OutFile,"rotate z %d\n",InvertY(-psi));
    if( phi )   fprintf(OutFile,"rotate y %d\n",phi);
    if( theta ) fprintf(OutFile,"rotate x %d\n",InvertY(-theta));

    if( (temp = (int)(100.0*DialValue[4])) )
        fprintf(OutFile,"translate x %d\n",temp);
    if( (temp = (int)(100.0*DialValue[5])) )
        fprintf(OutFile,"translate y %d\n",InvertY(-temp));

    if( DialValue[3] != 0.0 )
    {   if( DialValue[3]<0.0 )
        {   temp = (int)(100*DialValue[3]);
        } else temp = (int)(100*MaxZoom*DialValue[3]);
        fprintf(OutFile,"zoom %d\n",temp+100);
    }
    putc('\n',OutFile);

    /* Rendering */
    if( DrawAxes || DrawBoundBox || DrawUnitCell )
        fprintf(OutFile,"colour axes [%d %d %d]\n",BoxR,BoxG,BoxB);
    fprintf(OutFile,"set axes %s\n", DrawAxes? "on":"off" );
    fprintf(OutFile,"set boundingbox %s\n", DrawBoundBox? "on":"off" );
    fprintf(OutFile,"set unitcell %s\n", DrawUnitCell? "on":"off" );

    if( Database->hlist )
    {   fputs("set hbond ",OutFile);
        fputs(HBondMode?"backbone":"sidechain",OutFile);
        fputc('\n',OutFile);
    }

    if( Database->slist )
    {   fputs("set ssbond ",OutFile);
        fputs(SSBondMode?"backbone":"sidechain",OutFile);
        fputc('\n',OutFile);
    }

    fputs("set bondmode and\ndots off\n\n",OutFile); 
    fputs("\n# Avoid Colour Problems!\nselect all\n",OutFile);
    fputs("colour bonds none\ncolour backbone none\n",OutFile);
    fputs("colour hbonds none\ncolour ssbonds none\n",OutFile);
    fputs("colour ribbons none\ncolour white\n",OutFile);
    SelectAll = True;

    WriteScriptAtoms();
    if( UseSlabPlane && (SlabMode==SlabSection) )
    {   /* Section Mode Slabbing! */
        fclose(OutFile);
#ifdef APPLEMAC
        SetFileInfo(name,'RSML','RSML',133);
#endif
        return(True);
    }

    WriteScriptBonds();
    WriteScriptRibbons();
    WriteScriptBackbone();
    WriteScriptLabels();
    fputc('\n',OutFile);
    
    WriteScriptHBonds("ssbonds",Database->slist);
    WriteScriptHBonds("hbonds",Database->hlist);
    WriteScriptAll();

    fclose(OutFile);
#ifdef APPLEMAC
    SetFileInfo(name,'RSML','RSML',133);
#endif
    return( True );
}



#ifdef FUNCPROTO
static void OutputKinemageVector(Atom __far*,Atom __far*,Chain __far*,int);
static void WriteKinemageBonds( Chain __far* );
static void WriteKinemageSpheres( Chain __far* );
static void WriteKinemageRibbons( Chain __far* );
static void WriteKinemageLabels( Chain __far* );
#endif


static char *FindKinemageCol( r, g, b )
    int r, g, b;
{
    register Long dist,best;
    register int dr,dg,db;
    register int i,res;

    res = 0;
    dr = KinemageCol[0].r - r;
    dg = KinemageCol[0].g - g;
    db = KinemageCol[0].b - b;
    best = (Long)dr*dr + (Long)dg*dg + (Long)db*db;

    for( i=1; i<MAXMAGECOL; i++ )
    {   dr = KinemageCol[i].r - r;
        dg = KinemageCol[i].g - g;
        db = KinemageCol[i].b - b;
        dist = (Long)dr*dr + (Long)dg*dg + (Long)db*db;  

        if( dist < best ) 
        {   best = dist;
            res = i;
        }
    }
    return( KinemageCol[res].name );
}


static char *GetKinemageCol( col )
    int col;
{
    register ShadeDesc *ptr;

    ptr = Shade + Colour2Shade(col);
    return( FindKinemageCol(ptr->r,ptr->g,ptr->b) );
}




static void OutputKinemageVector( src, dst, chain, col )
    Atom __far *src, __far *dst;  
    Chain __far *chain;
    int col;
{
    register Group __far *group;
    register Group __far *sgrp;
    register Group __far *dgrp;
    register Atom __far *aptr;
    register Real x, y, z;
    register char *col1;
    register char *col2;


    /* Determine Chain and Groups */
    sgrp = dgrp = (Group __far*)0;

    if( chain )
    {   for( group=chain->glist; group; group=group->gnext )
        {   for( aptr=group->alist; aptr; aptr=aptr->anext )
            {   if( aptr == src ) sgrp = group;
                if( aptr == dst ) dgrp = group;
            }
            if( sgrp && dgrp ) break;
        }
    } else for( chain=Database->clist; chain; chain=chain->cnext )
    {   for( group=chain->glist; group; group=group->gnext )
        {   for( aptr=group->alist; aptr; aptr=aptr->anext )
            {   if( aptr == src ) sgrp = group;
                if( aptr == dst ) dgrp = group;
            }
            if( sgrp && dgrp ) break;
        }
        if( group ) break;
    }


    if( !sgrp || !dgrp ) 
        return;

    if( !col )
    {   col1 = GetKinemageCol(src->col);
        col2 = GetKinemageCol(dst->col);
    } else col1 = col2 = GetKinemageCol(col);

    if( (col1!=MageCol) && (col2==MageCol) )
    {   aptr = src;  src = dst;  dst = aptr;
        col2 = col1;  col1 = MageCol;
    }


    if( col1 != MageCol )
    {   fprintf(OutFile,"@vectorlist {} color= %s\n",col1);
        MagePrev = (Atom __far*)0;
    }

    if( src != MagePrev )
    {   if( MainGroupCount>1 )
        {   fprintf(OutFile,"{%.4s %.3s %d}", ElemDesc[src->refno], 
                    Residue[sgrp->refno], sgrp->serno );
        } else fprintf(OutFile,"{%.4s %d}",ElemDesc[src->refno],src->serno);
        fprintf(OutFile," P %g %g %g\n", src->xorg/250.0, 
                InvertY(src->yorg)/250.0, src->zorg/250.0 );
    }

    if( col1 != col2 )
    {   x = (src->xorg+dst->xorg)/500.0;
        y = (src->yorg+dst->yorg)/500.0;
        z = (src->zorg+dst->zorg)/500.0;

        fprintf(OutFile,"{} L %g %g %g\n", x, InvertY(y), z );
        fprintf(OutFile,"@vectorlist {} color= %s\n",col2);
        fprintf(OutFile,"{} P %g %g %g\n", x, InvertY(y), z );
    }

    if( MainGroupCount>1 )
    {   fprintf(OutFile,"{%.4s %.3s %d}", ElemDesc[dst->refno],
                Residue[dgrp->refno], dgrp->serno );
    } else fprintf(OutFile,"{%.4s %d}",ElemDesc[dst->refno],dst->serno);
    fprintf(OutFile," L %g %g %g\n", dst->xorg/250.0,
            InvertY(dst->yorg)/250.0, dst->zorg/250.0 );

    MagePrev = dst;
    MageCol = col2;
}


static void WriteKinemageBonds( chain )
    Chain __far *chain;
{
    register Bond __far *bptr;
    register Bond __far *flag;

    MagePrev = (Atom __far*)0;  
    MageCol = (char*)0;

    ForEachBond
        if( KinemageFlag || (bptr->flag&DrawBondFlag) )
        {   if( !MagePrev ) 
                fputs("@subgroup {wireframe} dominant\n",OutFile);
            OutputKinemageVector(bptr->srcatom, bptr->dstatom, 
                                 chain, bptr->col);
        }

    if( !chain->blist ) 
        return;

    /* Test for displayed backbone */
    for( flag=chain->blist; flag; flag=flag->bnext )
        if( flag->flag & DrawBondFlag ) break;
    if( !KinemageFlag && !flag ) return;

    MagePrev = (Atom __far*)0;  
    MageCol = (char*)0;

    for( bptr=chain->blist; bptr; bptr=bptr->bnext )
        if( KinemageFlag || (bptr->flag&DrawBondFlag) )
        {   if( !MagePrev )
            {   fputs("@subgroup {alpha trace} dominant",OutFile);
                fputs( (KinemageFlag && !flag)? " off\n":"\n",OutFile);
            }
            OutputKinemageVector(bptr->srcatom,bptr->dstatom,
                                 chain,bptr->col);
        }
}


static void WriteKinemageSpheres( chain )
    Chain __far *chain;
{
    register Group __far *group;
    register Atom __far *aptr;
    register char *col;
    register int radius;

    MageCol = (char*)0;
    for( group=chain->glist; group; group=group->gnext )
        for( aptr=group->alist; aptr; aptr=aptr->anext )
            if( aptr->flag & SphereFlag )
            {   if( !MageCol )
                    fputs("@subgroup {CPK spheres} dominant\n",OutFile);

                col = GetKinemageCol(aptr->col);
                if( (col!=MageCol) || (aptr->radius!=radius) )
                {   fprintf(OutFile,"@balllist {} color= %s radius= %g\n",
                                    col, aptr->radius/250.0);
                    radius = aptr->radius;
                    MageCol = col;
                }

                if( MainGroupCount>1 )
                {   fprintf(OutFile,"{%.4s %.3s %d}", ElemDesc[aptr->refno],
                                    Residue[group->refno], group->serno );
                } else fprintf(OutFile,"{%.4s %d}",ElemDesc[aptr->refno],
                                    aptr->serno);
                fprintf(OutFile," %g %g %g\n", aptr->xorg/250.0,
                       InvertY(aptr->yorg)/250.0, aptr->zorg/250.0 );
            }
}


static void WriteKinemageRibbons( chain )
    Chain __far *chain;
{
}



static void WriteKinemageLabels( chain )
    Chain __far *chain;
{
    register Group __far *group;
    register Atom __far *aptr;
    register Label *label;
    register char *col;

    auto char buffer[256];

    MageCol = (char*)0;
    for( group=chain->glist; group; group=group->gnext )
        for( aptr=group->alist; aptr; aptr=aptr->anext )
            if( aptr->label )
            {   if( !MageCol )
                    fputs("@subgroup {labels} dominant\n",OutFile);
                if( UseLabelCol )
                {   col = FindKinemageCol(LabR,LabG,LabB);
                } else col = GetKinemageCol(aptr->col);

                if( col != MageCol )
                    fprintf(OutFile,"@labellist {} color= %s\n",col);
                label = (Label*)aptr->label;
                FormatLabel(chain,group,aptr,label->label,buffer);
                fprintf(OutFile,"{%s} %g %g %g\n", buffer, aptr->xorg/250.0, 
                        InvertY(aptr->yorg)/250.0, aptr->zorg/250.0);
                MageCol = col;
            }
}


static void OutputKinemageUnitCell()
{
}


static void WriteKinemageData()
{
    register DotStruct __far *ptr;
    register HBond __far *hptr;
    register Real dx, dy, dz;
    register char *col;
    register int i;

    /* Hydrogen Bonds */
    for( hptr=Database->hlist; hptr; hptr=hptr->hnext )
        if( hptr->flag & DrawBondFlag ) break;

    if( KinemageFlag && Database->hlist )
    {   fputs("@group {h-bonds}",OutFile);
        fputs( hptr? "\n" : "off", OutFile );

        fputs("@subgroup {sidechain} dominant",OutFile);
        fputs( HBondMode? " off\n" : "\n", OutFile);
        for( hptr=Database->hlist; hptr; hptr=hptr->hnext )
            OutputKinemageVector(hptr->src,hptr->dst,
                                 (Chain __far*)0, hptr->col);

        fputs("@subgroup {mainchain} dominant",OutFile);
        fputs( HBondMode? "\n" : " off\n", OutFile);
        for( hptr=Database->hlist; hptr; hptr=hptr->hnext )
            OutputKinemageVector(hptr->srcCA,hptr->dstCA,
                                 (Chain __far*)0, hptr->col);
    } else if( hptr )
    {   fputs("@group {h-bonds} dominant\n",OutFile);
        for( hptr=Database->hlist; hptr; hptr=hptr->hnext )
            if( hptr->flag&DrawBondFlag )
                OutputKinemageVector( HBondMode?hptr->srcCA:hptr->src,
                                      HBondMode?hptr->dstCA:hptr->dst,
                                     (Chain __far*)0, hptr->col);
    }


    /* Disulphide Bridges */
    for( hptr=Database->slist; hptr; hptr=hptr->hnext )
        if( hptr->flag & DrawBondFlag ) break;

    if( KinemageFlag && Database->slist )
    {   fputs("@group {S-S bridges}",OutFile);
        fputs( hptr? "\n" : "off", OutFile );

        fputs("@subgroup {sidechain} dominant",OutFile);
        fputs( SSBondMode? " off\n" : "\n", OutFile);
        for( hptr=Database->slist; hptr; hptr=hptr->hnext )
            OutputKinemageVector(hptr->src,hptr->dst,
                                 (Chain __far*)0, hptr->col);

        fputs("@subgroup {mainchain} dominant",OutFile);
        fputs( SSBondMode? "\n" : " off\n", OutFile);
        for( hptr=Database->slist; hptr; hptr=hptr->hnext )
            OutputKinemageVector(hptr->srcCA,hptr->dstCA,
                                 (Chain __far*)0, hptr->col);
    } else if( hptr )
    {   fputs("@group {S-S bridges} dominant\n",OutFile);
        for( hptr=Database->slist; hptr; hptr=hptr->hnext )
            if( hptr->flag&DrawBondFlag )
                OutputKinemageVector( SSBondMode?hptr->srcCA:hptr->src,
                                      SSBondMode?hptr->dstCA:hptr->dst,
                                     (Chain __far*)0, hptr->col);
    }


    /* Dot Surfaces */
    if( DrawDots )
    {   fputs("@group {dot surface} dominant\n",OutFile);
        MageCol = (char*)0;

        for( ptr=DotPtr; ptr; ptr=ptr->next )
            for( i=0; i<ptr->count; i++ )
            {   col = GetKinemageCol(ptr->col[i]);
                if( col != MageCol )
                {   fprintf(OutFile,"@dotlist {} color= %s\n",col);
                    MageCol = col;
                }
                fprintf(OutFile, "{} %g %g %g\n", ptr->xpos[i]/250.0,
                        InvertY(ptr->ypos[i])/250.0, ptr->zpos[i]/250.0 );
            }
    }

    /* Draw `Background' Objects */
    if( !KinemageFlag && !DrawAxes && 
        !DrawBoundBox && !DrawUnitCell )
        return;

    dx = MaxX/250.0;  dy = MaxY/250.0;  dz = MaxZ/250.0;
    MageCol = FindKinemageCol( BoxR, BoxG, BoxB );

    if( DrawAxes || KinemageFlag )
    {   fputs("@group {coord axes} dominant",OutFile);
        fputs( (DrawAxes?"\n":" off\n"), OutFile );
        fprintf(OutFile,"@vectorlist {} color= %s\n",MageCol);

        fprintf(OutFile,"{} P %g 0 0\n{} L %g 0 0\n",-dx,dx);
        fprintf(OutFile,"{} P 0 %g 0\n{} L 0 %g 0\n",-dy,dy);
        fprintf(OutFile,"{} P 0 0 %g\n{} L 0 0 %g\n",-dz,dz);
    }

    if( DrawBoundBox || KinemageFlag )
    {   fputs("@group {bound box} dominant",OutFile);
        fputs( (DrawAxes?"\n":" off\n"), OutFile );
        fprintf(OutFile,"@vectorlist {} color= %s\n",MageCol);

        fprintf(OutFile,"{} P %g %g %g\n",-dx,-dy,-dz);
        fprintf(OutFile,"{} L %g %g %g\n", dx,-dy,-dz);
        fprintf(OutFile,"{} L %g %g %g\n", dx, dy,-dz);
        fprintf(OutFile,"{} L %g %g %g\n",-dx, dy,-dz);
        fprintf(OutFile,"{} L %g %g %g\n",-dx,-dy,-dz);
        fprintf(OutFile,"{} L %g %g %g\n",-dx,-dy, dz);
        fprintf(OutFile,"{} L %g %g %g\n", dx,-dy, dz);
        fprintf(OutFile,"{} L %g %g %g\n", dx, dy, dz);
        fprintf(OutFile,"{} L %g %g %g\n",-dx, dy, dz);
        fprintf(OutFile,"{} L %g %g %g\n",-dx,-dy, dz);

        fprintf(OutFile,"{} P %g %g %g\n", dx,-dy,-dz);
        fprintf(OutFile,"{} L %g %g %g\n", dx,-dy, dz);
        fprintf(OutFile,"{} P %g %g %g\n", dx, dy,-dz);
        fprintf(OutFile,"{} L %g %g %g\n", dx, dy, dz);
        fprintf(OutFile,"{} P %g %g %g\n",-dx, dy,-dz);
        fprintf(OutFile,"{} L %g %g %g\n",-dx, dy, dz);
    }

    if( DrawUnitCell || KinemageFlag )
    {   fputs("@group {unit cell} dominant",OutFile);
        fputs( (DrawAxes?"\n":" off\n"), OutFile );
        fprintf(OutFile,"@vectorlist {} color= %s\n",MageCol);
        OutputKinemageUnitCell();
    }
}


int WriteKinemageFile( name )
    char *name;
{
    register Chain __far *chain;
    register Real zoom;

    if( !Database )
        return(True);

    OutFile = fopen(name,"w");
    if( !OutFile )
    {   FatalScriptError(name);
        return(False);
    }

    fputs("@text\nRasMol v2.5 generated Kinemage\n \n",OutFile);
    if( *InfoMoleculeName )
        fprintf(OutFile,"Molecule name ....... %s\n",InfoMoleculeName);
    if( *InfoClassification )
        fprintf(OutFile,"Classification ...... %s\n",InfoClassification);
    if( *InfoIdentCode )
        fprintf(OutFile,"Brookhaven Code ..... %s\n",InfoIdentCode);

    fputs("@kinemage 1\n@caption RasMol v2.5 generated Kinemage\n",OutFile);
    fputs("@onewidth\n",OutFile);

    if( DialValue[3] != 0.0 )
    {   if( DialValue[3]<0.0 )
        {   zoom = DialValue[3];
        } else zoom = MaxZoom*DialValue[3];
        fprintf(OutFile,"zoom %g\n",zoom+1.0);
    }

    if( InfoChainCount>1 )
    {   for( chain=Database->clist; chain; chain=chain->cnext )
        {   fprintf(OutFile,"@group {chain %c}\n",chain->ident);
            WriteKinemageSpheres( chain );
            WriteKinemageBonds( chain );
            WriteKinemageRibbons( chain );
            WriteKinemageLabels( chain );
        }
    } else
    {   fputs("@group {molecule}\n",OutFile);
        chain = Database->clist;
        WriteKinemageSpheres( chain );
        WriteKinemageBonds( chain );
        WriteKinemageRibbons( chain );
        WriteKinemageLabels( chain );
    }

    WriteKinemageData();
    fclose(OutFile);
#ifdef APPLEMAC
    SetFileInfo(name,'MAGE','TEXT',135);
#endif
    return( True );
}

 int WritePOVRayFile( name )
    char *name;
{
    return( True );
}

