/* molecule.c
 * RasMol2 Molecular Graphics
 * Roger Sayle, October 1994
 * Version 2.5
 */
#include "rasmol.h"

#ifdef IBMPC
#include <windows.h>
#include <malloc.h>
#endif
#ifdef APPLEMAC
#include <Types.h>
#endif
#ifndef sun386
#include <stdlib.h>
#endif

#include <string.h>
#include <ctype.h>
#include <stdio.h>
#include <math.h>

#define MOLECULE
#include "molecule.h"
#include "command.h"
#include "abstree.h"
#include "transfor.h"
#include "render.h"


#define GroupPool    8
#define HBondPool   32
#define BondPool    32
#define AtomPool    32

#define NoLadder     0x00
#define ParaLadder   0x01
#define AntiLadder   0x02

#define Cos70Deg     0.34202014332567

#define FeatHelix    1
#define FeatSheet    2
#define FeatTurn     3

#define SourceNone   0
#define SourcePDB    1
#define SourceCalc   2

#define MaxHBondDist   ((Long)300*300)
#define MaxBondDist    ((Long)475*475)
#define MinBondDist    ((Long)100*100)
#define AbsMaxBondDist 500

typedef struct {
	int init, term;
	char chain;
	char type;
	} FeatEntry;

#ifdef APPLEMAC
#define AllocSize   256
typedef struct _AllocRef {
	struct _AllocRef *next;
	void *data[AllocSize];
	int count;
	} AllocRef;
static AllocRef *AllocList;  
#endif


#define FeatSize    32
typedef struct _Feature {
	    struct _Feature __far *fnext;
	FeatEntry data[FeatSize];
	int count;
	} Feature;

typedef struct {
	  char src[4];
	  char dst[4];
	  } ConvTable;

#define MAXALCATOM   5
static ConvTable AlcAtomTable[MAXALCATOM] = {
    { { 'S', 'O', '2', ' ' }, { ' ', 'S', '2', ' ' } },  /*  1 */
    { { 'C', 'A', 'R', ' ' }, { ' ', 'C', ' ', ' ' } },  /*  2 */
    { { 'N', 'A', 'R', ' ' }, { ' ', 'N', ' ', ' ' } },  /*  3 */
    { { 'N', 'A', 'M', ' ' }, { ' ', 'N', ' ', ' ' } },  /*  4 */
    { { 'N', 'P', 'L', '3' }, { ' ', 'N', '3', ' ' } },  /*  5 */
				 };

static Molecule __far *FreeMolecule;
static HBond __far *FreeHBond;
static Chain __far *FreeChain;
static Group __far *FreeGroup;
static Atom __far *FreeAtom;
static Bond __far *FreeBond;

static char PDBInsert;
static Feature __far *FeatList;
static HBond __far * __far *CurHBond;
static Molecule __far *CurMolecule;
static Chain __far *CurChain;
static Group __far *CurGroup;
static Atom __far *CurAtom;
static Atom __far *ConnectAtom;
static int InfoStrucSource;
static int HMinMaxFlag;
static int MMinMaxFlag;
static int MemSize;

static char Record[82];
static FILE *DataFile;

/* Macros for commonly used loops */
#define ForEachAtom  for(chain=Database->clist;chain;chain=chain->cnext) \
		     for(group=chain->glist;group;group=group->gnext)    \
		     for(aptr=group->alist;aptr;aptr=aptr->anext)
#define ForEachBond  for(bptr=Database->blist;bptr;bptr=bptr->bnext)


#ifdef APPLEMAC
/* External RasMac Function Declaration! */
void SetFileInfo( char*, OSType, OSType, short );
#endif


static void FatalDataError(ptr)
    char *ptr;
{
    char buffer[80];

    sprintf(buffer,"Database Error: %s!",ptr);
    RasMolFatalExit(buffer);
}


static void FetchRecord()
{
    register char *ptr;
    register int ch;

    ptr = Record + 1;
    if( !feof(DataFile) )
    {   do {
	    ch = getc(DataFile);
	    if( (ch=='\n') || (ch==EOF) )
	    {   *ptr = 0;
		return;
            } else if( (ch=='\r') )
            {   ch = getc(DataFile);
                if( ch != '\n' )
                    ungetc(ch,DataFile);
                *ptr = 0;
                return;
	    } else *ptr++ = ch;
	} while( ptr < Record+80 );

	/* skip to the end of the line! */
	do { ch = getc(DataFile);
	} while( (ch!='\n') && (ch!='\r') && (ch!=EOF) );

        if( ch=='\r' )
        {   ch = getc(DataFile);
            if( ch != '\n' )
                ungetc(ch,DataFile);
        }
    }
    *ptr = 0;
}


static void ExtractString( len, src, dst )
    int len;  char *src, *dst;
{
    register char *ptr;
    register char ch;
    register int i;

    ptr = dst;
    for( i=0; i<len; i++ )
    {   if( (ch = *src++) )
	{   if( (*dst++ = ch) != ' ' ) 
		ptr = dst;
	} else break;
    }
    *ptr = 0;
}


static Long ReadValue( pos, len )
    int pos, len;
{
    register Long result;
    register char *ptr;
    register char ch;
    register int neg;

    result = 0;
    neg = False;
    ptr = Record+pos;
    while( len-- )
    {   ch = *ptr++;
	if( (ch>='0') && (ch<='9') )
	{   result = (10*result)+(ch-'0');
	} else if( ch=='-' )
	    neg = True;
    }
    return( neg? -result : result );
}


void DescribeMolecule()
{
    char buffer[40];

    if( CommandActive )
	WriteChar('\n');
    CommandActive=False;

    if( *InfoMoleculeName )
    {   WriteString("Molecule name ....... ");
	WriteString(InfoMoleculeName);
	WriteChar('\n');
    }

    if( *InfoClassification )
    {   WriteString("Classification ...... ");
	WriteString(InfoClassification);
	WriteChar('\n');
    }

    if( Database && (MainGroupCount>1) )
    {   WriteString("Secondary Structure . ");
	if( InfoStrucSource==SourceNone )
	{   WriteString("No Assignment\n");
	} else if( InfoStrucSource==SourcePDB )
	{   WriteString("PDB Data Records\n");
	} else WriteString("Calculated\n");
    }


    if( *InfoIdentCode )
    {   WriteString("Brookhaven Code ..... ");
	WriteString(InfoIdentCode);
	WriteChar('\n');
    }

    if( InfoChainCount>1 )
    {   sprintf(buffer,"Number of Chains .... %d\n",InfoChainCount);
	WriteString(buffer);
    }

    sprintf(buffer,"Number of Groups .... %d",MainGroupCount);
    WriteString(buffer);
    if( HetaAtomCount )
    {   sprintf(buffer," (%d)\n",HetaGroupCount);
	WriteString(buffer);
    } else WriteChar('\n');

    sprintf(buffer,"Number of Atoms ..... %ld",MainAtomCount);
    WriteString(buffer);
    if( HetaAtomCount )
    {   sprintf(buffer," (%d)\n",HetaAtomCount);
	WriteString(buffer);
    } else WriteChar('\n');

    if( InfoBondCount )
    {   sprintf(buffer,"Number of Bonds ..... %ld\n",InfoBondCount);
	WriteString(buffer);
    }

    if( InfoSSBondCount != -1 )
    {   WriteString("Number of Bridges ... ");
	sprintf(buffer,"%d\n\n",InfoSSBondCount);
	WriteString(buffer);
    }

    if( InfoHBondCount != -1 )
    {   WriteString("Number of H-Bonds ... ");
	sprintf(buffer,"%d\n",InfoHBondCount);
	WriteString(buffer);
    }

    if( InfoHelixCount != -1 )
    {   WriteString("Number of Helices ... ");
	sprintf(buffer,"%d\n",InfoHelixCount);
	WriteString(buffer);

	WriteString("Number of Strands ... ");
	sprintf(buffer,"%d\n",InfoLadderCount);
	WriteString(buffer);

	WriteString("Number of Turns ..... ");
	sprintf(buffer,"%d\n",InfoTurnCount);
	WriteString(buffer);
    }
}


#ifdef APPLEMAC
/* Avoid System Memory Leaks! */
static void RegisterAlloc( data )
    void *data;
{
    register AllocRef *ptr;
    
    if( !AllocList || (AllocList->count==AllocSize) )
    {   ptr = (AllocRef *)_fmalloc( sizeof(AllocRef) );
	if( !ptr ) FatalDataError("Memory allocation failed");
	
	ptr->next = AllocList;
	ptr->data[0] = data;
	ptr->count = 1;
	AllocList = ptr;
    } else AllocList->data[AllocList->count++] = data;
}
#else
#define RegisterAlloc(x)
#endif


static void CreateChain( ident )
    char ident;
{
    register Chain __far *prev;

    if( !CurMolecule )
    {   if( !(CurMolecule = FreeMolecule) )
	{   MemSize += sizeof(Molecule);
	    CurMolecule = (Molecule __far *)_fmalloc(sizeof(Molecule));
	    if( !CurMolecule ) FatalDataError("Memory allocation failed");
	    RegisterAlloc( CurMolecule );
	} else FreeMolecule = (void __far*)0;

	CurChain = (void __far*)0;
	CurMolecule->slist = (void __far*)0;
	CurMolecule->hlist = (void __far*)0;
	CurMolecule->blist = (void __far*)0;
	CurMolecule->clist = (void __far*)0;
	Database = CurMolecule;
    }

    /* Handle chain breaks! */
    if( !(prev=CurChain) )
	if( (prev=CurMolecule->clist) )
	    while( prev->cnext )
		prev = prev->cnext;

    if( !(CurChain = FreeChain) )
    {   MemSize += sizeof(Chain);
	CurChain = (Chain __far *)_fmalloc(sizeof(Chain));
	if( !CurChain ) FatalDataError("Memory allocation failed");
	RegisterAlloc( CurChain );
    } else FreeChain = FreeChain->cnext;

    if( prev )
    {   prev->cnext = CurChain;
    } else CurMolecule->clist = CurChain;
    CurChain->cnext = (void __far*)0;
     
    CurChain->ident = ident;
    CurChain->glist = (void __far*)0;
    CurChain->blist = (void __far*)0;
    ConnectAtom = (void __far*)0;
    CurGroup = (void __far*)0;
    InfoChainCount++;
}


static void CreateGroup( pool )
    int pool;
{
    register Group __far *ptr;
    register int i;

    if( !(ptr = FreeGroup) )
    {   MemSize += pool*sizeof(Group);
	ptr = (Group __far *)_fmalloc( pool*sizeof(Group) );
	if( !ptr ) FatalDataError("Memory allocation failed");
	RegisterAlloc( ptr );
	for( i=1; i<pool; i++ )
	{   ptr->gnext = FreeGroup;
	    FreeGroup = ptr++;
	} 
    } else FreeGroup = ptr->gnext;
    
    if( CurGroup )
    {   CurGroup->gnext = ptr;
    } else CurChain->glist = ptr;
    CurGroup = ptr;

    CurAtom = (void __far*)0;
    ptr->gnext = (void __far*)0;
    ptr->alist = (void __far*)0;
    ptr->struc = 0;
    ptr->flag = 0;
    ptr->col1 = 0;
    ptr->col2 = 0;
}


static Atom __far *CreateAtom()
{
    register Atom __far *ptr;
    register int i;

    if( !(ptr = FreeAtom) )
    {   MemSize += AtomPool*sizeof(Atom);
	ptr = (Atom __far *)_fmalloc( AtomPool*sizeof(Atom) );
	if( !ptr ) FatalDataError("Memory allocation failed");
	RegisterAlloc( ptr );
	for( i=1; i<AtomPool; i++ )
	{   ptr->anext = FreeAtom;
	    FreeAtom = ptr++;
	} 
    } else FreeAtom = ptr->anext;

    if( CurAtom )
    {   CurAtom->anext = ptr;
    } else CurGroup->alist = ptr;
    ptr->anext = (void __far*)0;
    CurAtom = ptr;

    SelectCount++;
    ptr->flag = SelectFlag | NonBondFlag;
    ptr->label = (void*)0;
    ptr->radius = 375;
    ptr->altl = ' ';
    ptr->mbox = 0;
    ptr->col = 0;

    return( ptr );
}


static void ProcessAtom( ptr )
    Atom __far *ptr;
{
    if( GetElemNumber(ptr) == 1 )
    {   ptr->flag |= HydrogenFlag;
	HasHydrogen = True;
    }

    if( !(ptr->flag&(HydrogenFlag|HeteroFlag)) )
	ptr->flag |= NormAtomFlag;

#ifdef INVERT
    ptr->yorg = -ptr->yorg;
#endif

    if( HMinMaxFlag || MMinMaxFlag )
    {   if( ptr->xorg < MinX ) 
	{   MinX = ptr->xorg;
	} else if( ptr->xorg > MaxX ) 
	    MaxX = ptr->xorg;

	if( ptr->yorg < MinY ) 
	{   MinY = ptr->yorg;
	} else if( ptr->yorg > MaxY ) 
	    MaxY = ptr->yorg;

	if( ptr->zorg < MinZ ) 
	{   MinZ = ptr->zorg;
	} else if( ptr->zorg > MaxZ ) 
	    MaxZ = ptr->zorg;
    } else 
    {   MinX = MaxX = ptr->xorg;
	MinY = MaxY = ptr->yorg;
	MinZ = MaxZ = ptr->zorg;
    }
	    
    if( ptr->flag & HeteroFlag )
    {   if( HMinMaxFlag )
	{   if( ptr->temp < MinHetaTemp ) 
	    {   MinHetaTemp = ptr->temp;
	    } else if( ptr->temp > MaxHetaTemp ) 
		MaxHetaTemp = ptr->temp;
	} else MinHetaTemp = MaxHetaTemp = ptr->temp;
	HMinMaxFlag = True;
	HetaAtomCount++;
    } else
    {   if( MMinMaxFlag )
	{   if( ptr->temp < MinMainTemp ) 
	    {   MinMainTemp = ptr->temp;
	    } else if( ptr->temp > MaxMainTemp ) 
		MaxMainTemp = ptr->temp;
	} else MinMainTemp = MaxMainTemp = ptr->temp;
	MMinMaxFlag = True;
	MainAtomCount++;
    }
}


static int FindResNo( ptr )
    char *ptr;
{
    register int refno;

    for( refno=0; refno<ResNo; refno++ )
	if( !strncmp(Residue[refno],ptr,3) )
	    return( refno );

    if( !strncmp(ptr,"CSH",3) ||
	!strncmp(ptr,"CYH",3) ||
	!strncmp(ptr,"CSM",3) )
    {   return(17);  /* cystine */
    } else if( !strncmp(ptr,"WAT",3) ||
	       !strncmp(ptr,"H2O",3) ||
	       !strncmp(ptr,"SOL",3) ||
	       !strncmp(ptr,"TIP",3) )
    /* check HHO, OHH, 0H2? */
    {   return(31);  /* Water (HOH) */
    } else if( !strncmp(ptr,"D2O",3) )
    {   return(32);  /* Heavy water (DOD) */
    } else if( !strncmp(ptr,"SUL",3) )
    {   return(33);  /* Sulphate (SO4) */
    } else if( !strncmp(ptr,"CPR",3) )
    {   return(11);  /* cis-proline */
    } else if( !strncmp(ptr,"TRY",3) )
	return(15);  /* tryptophan */
    

    if( ResNo++ == MAXRES )
	FatalDataError("Too many new residues");
    Residue[refno][0] = *ptr++;
    Residue[refno][1] = *ptr++;
    Residue[refno][2] = *ptr;
    return( refno );
}


static Long ReadPDBCoord( offset )
    int offset;
{
    register int len,neg;
    register Long result;
    register char *ptr;
    register char ch;

    result = 0;
    neg = False;
    len = 8;

    ptr = Record+offset;
    while( len-- )
    {   ch = *ptr++;
	if( (ch>='0') && (ch<='9') )
	{   result = (10*result)+(ch-'0');
	} else if( ch=='-' )
	    neg = True;
    }

    /* Handle Chem3D PDB Files! */
    if( Record[offset+3]=='.' )
	result /= 10;
    return( neg? -result : result );
}


static void ProcessPDBColourMask()
{
    register MaskDesc *ptr;
    register char *mask;
    register int i;

    if( MaskCount==MAXMASK )
	FatalDataError("Too many COLOR records in file");
    ptr = &UserMask[MaskCount];
    mask = ptr->mask;


    ptr->flags = 0;
    for( i=7; i<12; i++ )
	if( (*mask++ = Record[i]) != '#' )
	    ptr->flags |= SerNoFlag;

    for( i=13; i<21; i++ )
	*mask++ = Record[i];
    *mask++ = Record[22];

    for( i=23; i<27; i++ )
	if( (*mask++ = Record[i]) != '#' )
	    ptr->flags |= ResNoFlag;
    *mask++ = Record[27];

    ptr->r = (int)(ReadPDBCoord(31)>>2) + 5;
    ptr->g = (int)(ReadPDBCoord(39)>>2) + 5;
    ptr->b = (int)(ReadPDBCoord(47)>>2) + 5;
    ptr->radius = (short)(5*ReadValue(55,6))>>1;
    MaskCount++;
}


static Bond __far *ProcessBond( src, dst, flag )
    Atom __far *src, __far *dst;
    int flag;
{
    register Bond __far *ptr;
    register int i;

    if( !(ptr = FreeBond) )
    {   MemSize += BondPool*sizeof(Bond);
	ptr = (Bond __far *)_fmalloc( BondPool*sizeof(Bond) );
	if( !ptr ) FatalDataError("Memory allocation failed");
	RegisterAlloc( ptr );
	for( i=1; i<BondPool; i++ )
	{   ptr->bnext = FreeBond;
	    FreeBond = ptr++;
	} 
    } else FreeBond = ptr->bnext;
    
    ptr->flag = flag | SelectFlag;
    ptr->srcatom = src;
    ptr->dstatom = dst;
    ptr->radius = 0;
    ptr->col = 0;

    return( ptr );
}


static void ConnectAtoms( src, dst, flag )
    int src, dst, flag;
{
    register Chain __far *chain;
    register Group __far *group;
    register Atom __far *aptr;
    register Atom __far *sptr;
    register Atom __far *dptr;
    register Bond __far *bptr;
    register int done;

    if( src == dst )
	return;

    sptr = (void __far*)0;
    dptr = (void __far*)0;

    done = False;
    for( chain=Database->clist; chain && !done; chain=chain->cnext )
	for( group=chain->glist; group && !done; group=group->gnext )
	    for( aptr=group->alist; aptr; aptr=aptr->anext )
	    {   if( aptr->serno == src )
		{   sptr = aptr;
		    if( dptr )
		    {   done = True;
			break;
		    }
		} else if( aptr->serno == dst )
		{   dptr = aptr;
		    if( sptr )
		    {   done = True;
			break;
		    }
		}
	    }


    /* Both found! */
    if( done ) 
    {   /* Reset Non-bonded flags! */
	sptr->flag &= ~NonBondFlag;
	dptr->flag &= ~NonBondFlag;

	bptr = ProcessBond( sptr, dptr, flag );
	bptr->bnext = CurMolecule->blist;
	CurMolecule->blist = bptr;
	InfoBondCount++;
    }
}


static void PDBConnectAtoms( src, dst )
    int src, dst;
{
    register Bond __far *bptr;
    register int bs,bd;

    ForEachBond
    {   bs = bptr->srcatom->serno;
	bd = bptr->dstatom->serno;

	if( ((bs==src)&&(bd==dst)) || ((bs==dst)&&(bd==src)) )
	{   if( bptr->flag & NormBondFlag )
	    {  /* Convert Single to Double */
	       bptr->flag &= ~(NormBondFlag);
	       bptr->flag |= DoubBondFlag;
	    } else if( bptr->flag & DoubBondFlag )
	    {  /* Convert Double to Triple */
	       bptr->flag &= ~(DoubBondFlag);
	       bptr->flag |= TripBondFlag;
	    }
	    return;
	}
    }

    ConnectAtoms( src, dst, NormBondFlag );
}



static void ProcessPDBGroup( heta, serno )
    int heta, serno;
{
    PDBInsert = Record[27];
    if( !CurChain || (CurChain->ident!=Record[22]) )
	CreateChain( Record[22] );
    CreateGroup( GroupPool );

    CurGroup->refno = FindResNo( &Record[18] );
    CurGroup->serno = serno;

    /* Solvents should be hetero! */
    if( IsSolvent(CurGroup->refno) )
	heta = True;

    if( heta )
    {   HetaGroupCount++;
	if( HMinMaxFlag )
	{   if( serno > MaxHetaRes ) 
	    {   MaxHetaRes = serno;
	    } else if( serno < MinHetaRes )
		MinHetaRes = serno;
	} else MinHetaRes = MaxHetaRes = serno;
    } else 
    {   MainGroupCount++;
	if( MMinMaxFlag )
	{   if( serno > MaxMainRes )
	    {   MaxMainRes = serno;
	    } else if( serno < MinMainRes )
		MinMainRes = serno;
	} else MinMainRes = MaxMainRes = serno; 
    }
}


static void ProcessPDBAtom( heta )
    int heta;
{
    auto char name[4];
    register Bond __far *bptr;
    register Atom __far *ptr;
    register Long dx,dy,dz;
    register int serno;
    register int refno;
    register int i;

    /* Ignore Pseudo Atoms!! */
    if( (Record[13]==' ') && (Record[14]=='Q') )
	return; 

    dx = ReadPDBCoord(31);
    dy = ReadPDBCoord(39);
    dz = ReadPDBCoord(47);

    /* Ignore XPLOR Pseudo Atoms!! */
    if( (dx==9999000L) && (dy==9999000L) && (dz==9999000L) )
	return;

    serno = (int)ReadValue(23,4);
    if( !CurGroup || (CurGroup->serno!=serno) 
	|| (CurChain->ident!=Record[22]) 
	|| (PDBInsert!=Record[27]) )
	ProcessPDBGroup( heta, serno );

    /* Solvents should be hetero! */
    if( IsSolvent(CurGroup->refno) )
	heta = True;


    ptr = CreateAtom();
    if( isdigit(Record[14]) )
    {   name[1] = ToUpper(Record[13]);
	name[2] = name[3] = ' ';
	name[0] = ' ';

    } else for( i=0; i<4; i++ )
	name[i] = ToUpper(Record[i+13]);

    for( refno=0; refno<ElemNo; refno++ )
	if( !strncmp(ElemDesc[refno],name,4) )
	    break;

    if( refno==ElemNo )
    {   if( ElemNo++ == MAXELEM )
	    FatalDataError("Too many new atom types");

	for( i=0; i<4; i++ )
	    ElemDesc[refno][i] = name[i];
    }

    ptr->refno = refno;
    ptr->altl = Record[17];
    ptr->serno = (int)ReadValue(7,5);
    ptr->temp = (int)ReadValue(61,6);

    ptr->xorg =  dx/4;
    ptr->yorg =  dy/4;
    ptr->zorg = -dz/4;

    if( heta ) 
	ptr->flag |= HeteroFlag;

    ProcessAtom( ptr );
    if( IsAlphaCarbon(refno) && IsProtein(CurGroup->refno) )
    {   if( ConnectAtom )
	{   dx = ConnectAtom->xorg - ptr->xorg;
	    dy = ConnectAtom->yorg - ptr->yorg;
	    dz = ConnectAtom->zorg - ptr->zorg;

	    /* Break backbone if CA-CA > 7.00A */
	    if( dx*dx+dy*dy+dz*dz < (Long)1750*1750 )
	    {   bptr = ProcessBond(ptr,ConnectAtom,NormBondFlag);
		bptr->bnext = CurChain->blist;
		CurChain->blist = bptr;
	    } else ptr->flag |= BreakFlag;
	}
	ConnectAtom = ptr;
    } else if( IsSugarPhosphate(refno) && IsNucleo(CurGroup->refno) )
    {   if( ConnectAtom )
	{   bptr = ProcessBond(ConnectAtom,ptr,NormBondFlag);
	    bptr->bnext = CurChain->blist;
	    CurChain->blist = bptr;
	}
	ConnectAtom = ptr;
    }
}


static FeatEntry __far *AllocFeature()
{
    register Feature __far *ptr;

    if( !FeatList || (FeatList->count==FeatSize) )
    {   ptr = (Feature __far*)_fmalloc(sizeof(Feature));
	if( !ptr ) FatalDataError("Memory allocation failed");
	/* Features are always deallocated! */
	
	ptr->fnext = FeatList;
	ptr->count = 0;
	FeatList = ptr;
    } else ptr = FeatList;

    return( &(ptr->data[ptr->count++]) );
}


#ifdef FUNCPROTO
static void UpdateFeature( FeatEntry __far*, int );
#endif

static void UpdateFeature( ptr, mask )
    FeatEntry __far *ptr;  int mask;
{
    register Chain __far *chain;
    register Group __far *group;

    for( chain=Database->clist; chain; chain=chain->cnext )
	if( chain->ident == ptr->chain )
	{   group=chain->glist;
	    while( group && (group->serno<ptr->init) )
		group = group->gnext;

	    while( group && (group->serno<=ptr->term) )
	    {   group->struc |= mask;
		group = group->gnext;
	    }
	    return;
	}
}


static void ProcessPDBFeatures()
{
    register Feature __far *next;
    register Feature __far *ptr;
    register int i;

    InfoTurnCount = 0;
    InfoHelixCount = 0;
    InfoLadderCount = 0;
    InfoStrucSource = SourcePDB;

    for( ptr=FeatList; ptr; ptr=next )
    {    if( Database )
	     for( i=0; i<ptr->count; i++ )
		 if( ptr->data[i].type==FeatHelix )
		 {   UpdateFeature( &ptr->data[i], HelixFlag );
		     InfoHelixCount++;
		 } else if( ptr->data[i].type==FeatSheet )
		 {   UpdateFeature( &ptr->data[i], SheetFlag );
		     InfoLadderCount++;
		 } else /* FeatTurn */
		 {   UpdateFeature( &ptr->data[i], TurnFlag );
		     InfoTurnCount++;
		 }

	 /* Deallocate Memory */
	 next = ptr->fnext;
	 _ffree( ptr );
    }
}


int LoadPDBMolecule( fp )
    FILE *fp;
{
    register FeatEntry __far *ptr;
    register int srcatm, dstatm;
    register char *src, *dst;
    register int model;

    model = True;
    FeatList = (void __far*)0;
    DataFile = fp;


    do {   
	FetchRecord();

	if( !strncmp("ATOM",Record+1,4) )
	{   if( model ) ProcessPDBAtom(False);
	} else if( !strncmp("HETA",Record+1,4) )
	{   if( model ) ProcessPDBAtom(True);
	} else if( !strncmp("SHEE",Record+1,4) )
	{   if( !model ) continue;
	
	    /* Remaining SHEET record fields   */
	    /* 39-40 .... Strand Parallelism   */
	    /* 33 ....... Same Chain as 22?    */
	    ptr = AllocFeature();
	    ptr->type = FeatSheet;
	    ptr->chain = Record[22];
	    ptr->init = (int)ReadValue(23,4);
	    ptr->term = (int)ReadValue(34,4);
       
	} else if( !strncmp("HELI",Record+1,4) )
	{   if( !model ) continue;
	
	    /* Remaining HELIX record fields   */
	    /* 39-40 .... Helix Classification */
	    /* 32 ....... Same Chain as 20?    */
	    ptr = AllocFeature();
	    ptr->type = FeatHelix;
	    ptr->chain = Record[20];
	    ptr->init = (int)ReadValue(22,4);
	    ptr->term = (int)ReadValue(34,4);

	} else if( !strncmp("TURN",Record+1,4) )
	{   if( !model ) continue;
	
	    ptr = AllocFeature();
	    ptr->type = FeatTurn;
	    ptr->chain = Record[20];
	    ptr->init = (int)ReadValue(21,4);
	    ptr->term = (int)ReadValue(32,4);

	} else if( !strncmp("CONE",Record+1,4) )
	{   if( !model ) continue;

	    if( (srcatm = (int)ReadValue(7,5)) )
	    {   if( (dstatm = (int)ReadValue(12,5)) )
		    if( dstatm > srcatm )
			PDBConnectAtoms( srcatm, dstatm );

		if( !Record[17] ) continue;
		if( (dstatm = (int)ReadValue(17,5)) )
		    if( dstatm > srcatm )
			PDBConnectAtoms( srcatm, dstatm );

		if( !Record[22] ) continue;
		if( (dstatm = (int)ReadValue(22,5)) )
		    if( dstatm > srcatm )
			PDBConnectAtoms( srcatm, dstatm );

		if( !Record[27] ) continue;
		if( (dstatm = (int)ReadValue(27,5)) )
		    if( dstatm > srcatm )
			PDBConnectAtoms( srcatm, dstatm );
	    }

	} else if( !strncmp("COLO",Record+1,4) )
	{   ProcessPDBColourMask();
	} else if( !strncmp("TER",Record+1,3) )
	{   if( !Record[4] || (Record[4]==' ') )
	    {   ConnectAtom = (void __far*)0;
		CurGroup = (void __far*)0;
		CurChain = (void __far*)0;
	    }
	} else if( !strncmp("HEAD",Record+1,4) )
	{   ExtractString(40,Record+11,InfoClassification);
	    ExtractString( 4,Record+63,InfoIdentCode);
	    
	} else if( !strncmp("COMP",Record+1,4) )
	{   if( Record[10]==' ' )  /* First COMPND record */
		ExtractString(60,Record+11,InfoMoleculeName);

	} else if( !strncmp("CRYS",Record+1,4) )
	{   dst = InfoSpaceGroup;
	    for( src=Record+56; *src && src<Record+67; src++ )
		if( *src!=' ' ) *dst++ = *src;
	    *dst = 0;

	    InfoCellA = ReadValue( 7,9)/1000.0;
	    InfoCellB = ReadValue(16,9)/1000.0;
	    InfoCellC = ReadValue(25,9)/1000.0;

	    InfoCellAlpha = Deg2Rad*(ReadValue(34,7)/100.0);
	    InfoCellBeta =  Deg2Rad*(ReadValue(41,7)/100.0);
	    InfoCellGamma = Deg2Rad*(ReadValue(48,7)/100.0);

	} else if( !strncmp("ENDM",Record+1,4) )
	{   model = False;
	} else if( !strncmp("END",Record+1,3) )
	    if( !Record[4] || (Record[4]==' ') )
	    {   /* Treat END same as TER! */
		ConnectAtom = (void __far*)0;
		CurGroup = (void __far*)0;
		CurChain = (void __far*)0;
	    }
    } while( !feof(DataFile) );

    if( Database )
	strcpy(InfoFileName,DataFileName);
    if( FeatList ) ProcessPDBFeatures();
    return( True );
}


static int SimpleAtomType( type )
    char *type;
{
    register int i, refno;
    auto char name[4];

    name[2] = name[3] = ' ';
    if( type[1] && (type[1]!=' ') )
    {   name[0] = ToUpper(type[0]);
	name[1] = ToUpper(type[1]);
    } else
    {   name[1] = ToUpper(type[0]);
	name[0] = ' ';
    }

    for( refno=0; refno<ElemNo; refno++ )
	if( !strncmp(ElemDesc[refno],name,4) )
	    return( refno );

    if( ElemNo++ == MAXELEM )
	FatalDataError("Too many new atom types");

    for( i=0; i<4; i++ )
	ElemDesc[refno][i] = name[i];
    return( refno );
}


static void CreateMolGroup()
{
    strcpy(InfoFileName,DataFileName);

    CreateChain( ' ' );
    CreateGroup( 1 );

    CurGroup->refno = FindResNo( "MOL" );
    CurGroup->serno = 1;
	
    MinMainRes = MaxMainRes = 1;
    MinHetaRes = MaxHetaRes = 0;
    MainGroupCount = 1;
}


int LoadXYZMolecule( fp )
    FILE *fp;
{
    auto char type[12];
    auto double xpos, ypos, zpos;
    auto double charge, u, v, w;
    auto int atoms;

    register Atom __far *ptr;
    register char *src,*dst;
    register int i,count;


    DataFile = fp;

    /* Number of Atoms */
    FetchRecord();
    sscanf(Record+1,"%d",&atoms);

    /* Molecule (step) Description */
    FetchRecord();
    src = Record+1;
    while( *src == ' ' )
	src++;

    dst = InfoMoleculeName;
    for( i=0; i<78; i++ )
	if( *src ) *dst++ = *src++;
    *dst = '\0';

    if( atoms )
    {   CreateMolGroup();
	for( i=0; i<atoms; i++ )
	{   FetchRecord();
	    ptr = CreateAtom();
	    ptr->serno = i;

	    xpos = ypos = zpos = 0.0;
	    count = sscanf(Record+1,"%s %lg %lg %lg %lg %lg %lg %lg",
			   type, &xpos, &ypos, &zpos, &charge, &u, &v, &w );

	    ptr->refno = SimpleAtomType(type);
	    ptr->xorg =  (Long)(250.0*xpos);
	    ptr->yorg =  (Long)(250.0*ypos);
	    ptr->zorg = -(Long)(250.0*zpos);

	    if( (count==5) || (count==8) )
	    {   ptr->temp = (short)(100.0*charge);
	    } else ptr->temp = 0;
	    ProcessAtom( ptr );
	}
    }
    return( True );
}


static int FindSybylResNo( ptr )
    char *ptr;
{
    register char *src,*dst;
    register int i, j;
    auto char name[4];

    src = ptr;
    dst = name;
    if( ptr[1] && (ptr[1]!='.') )
    {   *dst++ = ToUpper(*src);  src++;
	*dst++ = ToUpper(*src);  src++;
    } else
    {   *dst++ = ToUpper(*src);  src++;
	*dst++ = ' ';
    }

    if( *src )
    {   src++;

	if( *src == 'a' )
	{   *dst++ = ' ';
	    *dst = ' ';
	} else if( *src == 'p' )
	{   *dst++ = '3';
	    *dst = ' ';
	} else
	{   *dst++ = *src++;
	    if( *src && (*src!='+') )
	    {   *dst = *src;
	    } else *dst = ' ';
	}
    } else
    {   *dst++ = ' ';
	*dst = ' ';
    }

    for( j=0; j<ElemNo; j++ )
	if( !strncmp(ElemDesc[j],name,4) )
	    return( j );

    if( ElemNo++ == MAXELEM )
	FatalDataError("Too many new atom types");
    for( i=0; i<4; i++ ) ElemDesc[j][i] = name[i];
    return( j );
}


int LoadMol2Molecule( fp )
    FILE *fp;
{
    auto int srcatm, dstatm;
    auto int features, sets, serno;
    auto int atoms, bonds, structs;
    auto double xpos, ypos, zpos;

    auto char name[8];
    auto char type[4];

    register Atom __far *ptr;
    register char *src, *dst;
    register int i;


    DataFile = fp;

    while( !feof(DataFile) )
    {   FetchRecord();
	if( !Record[1] || Record[1]=='#' )
	    continue;

	if( !strncmp("@<TRIPOS>MOLECULE",Record+1,17) )
	{   FetchRecord();  /* Molecule Name */
	    src = Record+1;
	    while( *src==' ' )
		src++;

	    dst = InfoMoleculeName;
	    while( (*dst++ = *src++) );

	    FetchRecord();
	    atoms = bonds = structs = features = sets = 0;
	    sscanf(Record+1,"%d %d %d %d %d", &atoms, &bonds, &structs,
					      &features, &sets );

	    FetchRecord();  /* Molecule Type  */
	    FetchRecord();  /* Charge Type    */

	} else if( !strncmp("@<TRIPOS>ATOM",Record+1,13) )
	{   if( !atoms ) continue;

	    CreateMolGroup();
	    for( i=0; i<atoms; i++ )
	    {    FetchRecord();
		 ptr = CreateAtom();

		 sscanf(Record+1,"%d %s %lg %lg %lg %s", &serno, name,
				  &xpos, &ypos, &zpos, type );

		 ptr->refno = FindSybylResNo( type );
		 ptr->serno = serno;
		 /* ptr->serno = i; */

		 ptr->xorg =  (Long)(250.0*xpos);
		 ptr->yorg =  (Long)(250.0*ypos);
		 ptr->zorg = -(Long)(250.0*zpos);
		 ProcessAtom( ptr );
	    }

	} else if( !strncmp("@<TRIPOS>BOND",Record+1,13) )
	    for( i=0; i<bonds; i++ )
	    {   FetchRecord();
		sscanf(Record+1,"%d %d %d %.2s",&serno,&srcatm,&dstatm,type);
		if( !strncmp(type,"ar",2) )
		{   ConnectAtoms(srcatm,dstatm,AromBondFlag);
		} else if( *type == '2' )
		{   ConnectAtoms(srcatm,dstatm,DoubBondFlag);
		} else /* *type == '1' */
		    ConnectAtoms(srcatm,dstatm,NormBondFlag);
	    }
    }
    return( True );
}



static int FindAlchemyResNo()
{
    register char *ptr;
    register int i, j;
    char name[4];

    ptr = Record+7;
    if( !isalpha(ptr[1]) )
    {   name[0] = ' ';
	for( i=0; i<3; i++ )
	    name[i+1] = ToUpper(ptr[i]);
	ptr = name;
    } else
    {   for( i=0; i<4; i++ )
	    ptr[i] = ToUpper(ptr[i]);

	for( i=0; i<MAXALCATOM; i++ )
	    if( !strncmp(AlcAtomTable[i].src,ptr,4) )
	    {   ptr = AlcAtomTable[i].dst;
		break;
	    }
    }

    for( j=0; j<ElemNo; j++ )
	if( !strncmp(ElemDesc[j],ptr,4) )
	    return( j );

    if( ElemNo++ == MAXELEM )
	FatalDataError("Too many new atom types");
    for( i=0; i<4; i++ ) ElemDesc[j][i] = ptr[i];
    return( j );
}


int LoadAlchemyMolecule( fp )
    FILE *fp;
{
    auto char type[12];
    auto int serno,srcatm,dstatm;
    register Atom __far *ptr;
    register int atoms, bonds;
    register int i;

    DataFile = fp;
    FetchRecord();
    atoms = (int)ReadValue(1,5);
    bonds = (int)ReadValue(14,5);
    ExtractString(38,Record+42,InfoMoleculeName);

    if( !atoms )
	return( False );

    CreateMolGroup();
    for( i=0; i<atoms; i++ )
    {   FetchRecord();
	ptr = CreateAtom();

	ptr->refno = FindAlchemyResNo();
	ptr->temp = (int)ReadValue(41,8);
	ptr->serno = (int)ReadValue(1,5);
	/* ptr->serno = i+1; */

	ptr->xorg =  ReadValue(13,7)/4;
	ptr->yorg =  ReadValue(22,7)/4;
	ptr->zorg = -ReadValue(31,7)/4;
	ProcessAtom( ptr );
    }

    for( i=0; i<bonds; i++ )
    {   FetchRecord();
	sscanf(Record+1,"%d %d %d %.10s",&serno,&srcatm,&dstatm,type);

	if( *type =='A' )             /* AROMATIC */
	{   ConnectAtoms(srcatm,dstatm,AromBondFlag);
	} else if( *type == 'D' )     /* DOUBLE */
	{   ConnectAtoms(srcatm,dstatm,DoubBondFlag);
	} else if( *type == 'T' )     /* TRIPLE */
	{   ConnectAtoms(srcatm,dstatm,TripBondFlag);
	} else /* (*type == 'S') */   /* SINGLE */
	    ConnectAtoms(srcatm,dstatm,NormBondFlag);
    }
    return( True );
}


int LoadCharmmMolecule( fp )
    FILE *fp;
{
    auto char buffer[9];
    register Atom __far *ptr;
    register int refno,serno;
    register int resno,atoms;
    register int init,chain;
    register int i;

    DataFile = fp;

    do { 
	FetchRecord();
    } while( Record[1]=='*' );
    atoms = (int)ReadValue(1,5);
    if( !atoms ) return False;

    MinHetaRes = MaxHetaRes = 0;
    strcpy(InfoFileName,DataFileName);
    MainGroupCount = 0;

    chain = 0;
    init = False;
    CurChain = NULL;
    for( serno=0; serno<atoms; serno++ )
    {   FetchRecord();

	if( !CurChain || strncmp(Record+52,buffer,9) )
	{   for( i=0; i<9; i++ )
		buffer[i] = Record[52+i];
	    CreateChain(chain+49);
	    chain++;
	}

	resno = (int)ReadValue(6,5);
	if( !CurGroup || (CurGroup->serno!=resno) )
	{   CreateGroup( GroupPool );
	    CurGroup->refno = FindResNo(Record+12);
	    CurGroup->serno = resno;

	    if( !init )
	    {   MinMainRes = resno;
		MaxMainRes = resno;
		init = True;
	    } else if( resno > MaxMainRes )
	    {   MaxMainRes = resno;
	    } else if( resno < MinMainRes )
		MinMainRes = resno;
	    MainGroupCount++;
	}

	ptr = CreateAtom();
	for( refno=0; refno<ElemNo; refno++ )
	    if( !strncmp(ElemDesc[refno],Record+16,4) )
		break;

	if( refno == ElemNo )
	{   if( ElemNo++ == MAXELEM )
		FatalDataError("Too many new atom types");
	    for( i=0; i<4; i++ )
		ElemDesc[refno][i] = Record[16+i];
	}

	ptr->refno = refno;
	ptr->temp = (int)ReadValue(61,9);
	ptr->serno = (int)ReadValue(1,5);
	/* ptr->serno = serno+1; */

	ptr->xorg =  ReadValue(21,8)/4;
	ptr->yorg =  ReadValue(31,8)/4;
	ptr->zorg = -ReadValue(41,8)/4;
	ProcessAtom( ptr );
    }
    return( True );
}


int LoadMDLMolecule( fp )
    FILE *fp;
{
    register Bond __far *bptr;
    register Atom __far *src;
    register Atom __far *dst;
    register Atom __far *ptr;

    register int atoms, bonds;
    register int srcatm,dstatm;
    register int i,type,done;
    register Long dx, dy, dz;
    register Card dist2;
    register Real scale;

    DataFile = fp;

    FetchRecord(); /* Molecule Name */
    ExtractString(78,Record+1,InfoMoleculeName);

    FetchRecord(); /* Program Details */
    FetchRecord(); /* Comments */

    FetchRecord();
    atoms = (int)ReadValue(1,3);
    bonds = (int)ReadValue(4,3);

    if( !atoms )
	return( False );

    CreateMolGroup();
    for( i=1; i<=atoms; i++ )
    {   FetchRecord();
	ptr = CreateAtom();
	ptr->refno = SimpleAtomType(Record+32);

	switch( (int)ReadValue(37,3) )
	{   case(1):  ptr->temp =  300;  break;
	    case(2):  ptr->temp =  200;  break;
	    case(3):  ptr->temp =  100;  break;
	    case(5):  ptr->temp = -100;  break;
	    case(6):  ptr->temp = -200;  break;
	    case(7):  ptr->temp = -300;  break;
	    default:  ptr->temp = 0;
	}
	ptr->serno = i;

	ptr->xorg =  ReadValue(1,10)/40;
	ptr->yorg =  ReadValue(11,10)/40;
	ptr->zorg = -ReadValue(21,10)/40;
	ProcessAtom( ptr );
    }

    for( i=0; i<bonds; i++ )
    {   FetchRecord();
	srcatm = (int)ReadValue(1,3);
	dstatm = (int)ReadValue(4,3);
	type = (int)ReadValue(7,3);

	if( type==2 )                 /* DOUBLE */
	{   ConnectAtoms(srcatm,dstatm,DoubBondFlag);
	} else if( type==3 )          /* TRIPLE */
	{   ConnectAtoms(srcatm,dstatm,TripBondFlag);
	} else if( type==4 )          /* AROMATIC */
	{   ConnectAtoms(srcatm,dstatm,AromBondFlag);
	} else                        /* SINGLE */
	    ConnectAtoms(srcatm,dstatm,NormBondFlag);
    }

    done = False;
    for( bptr=CurMolecule->blist; bptr; bptr=bptr->bnext )
	if( bptr->flag & NormBondFlag )
	{   src = bptr->srcatom;
	    dst = bptr->dstatom;
	    if( (src->refno==2) && (dst->refno==2) )
	    {   dx = dst->xorg - src->xorg;
		dy = dst->yorg - src->yorg;
		dz = dst->zorg - src->zorg;
		if( dx || dy || dz )
		{   dist2 = dx*dx + dy*dy + dz*dz;
		    scale = 385.0/sqrt(dist2);
		    break;
		}
	    }
	}

    if( bptr )
    {   for( ptr=CurGroup->alist; ptr; ptr=ptr->anext )
	{   ptr->xorg = (Long)(ptr->xorg*scale);
	    ptr->yorg = (Long)(ptr->yorg*scale);
	    ptr->zorg = (Long)(ptr->zorg*scale);
	}
	MinX = (Long)(MinX*scale);  MaxX = (Long)(MaxX*scale);
	MinY = (Long)(MinY*scale);  MaxY = (Long)(MaxY*scale);
	MinZ = (Long)(MinZ*scale);  MaxZ = (Long)(MaxZ*scale);
    }
    return( True );
}



int SavePDBMolecule( filename )
    char *filename;
{
    register double x, y, z;
    register Group __far *prev;
    register Chain __far *chain;
    register Group __far *group;
    register Atom __far *aptr;
    register char *ptr;
    register int count;
    register char ch;
    register int i;

    if( !Database )
	return( False );

    DataFile = fopen( filename, "w" );
    if( !DataFile )
    {   if( CommandActive )
	    WriteChar('\n');
	WriteString("Error: Unable to create file!\n\n");
	CommandActive=False;
	return( False );
    }

    if( *InfoClassification || *InfoIdentCode )
    {   fputs("HEADER    ",DataFile);

	ptr = InfoClassification;
	for( i=11; i<=50; i++ )
	    putc( (*ptr ? *ptr++ : ' '), DataFile );
	fprintf(DataFile,"13-JUL-93   %.4s\n",InfoIdentCode);
    }

    if( *InfoMoleculeName )
	fprintf(DataFile,"COMPND    %.60s\n",InfoMoleculeName);

    prev = (void __far*)0;

    count = 1;
    ForEachAtom
	if( aptr->flag&SelectFlag )
	{   if( prev && (chain->ident!=ch) )
		fprintf( DataFile, "TER   %5d      %.3s %c%4d \n", 
			 count++, Residue[prev->refno], ch, prev->serno);

	    if( aptr->flag&HeteroFlag )
	    {      fputs("HETATM",DataFile);
	    } else fputs("ATOM  ",DataFile);
	    fprintf( DataFile, "%5d %.4s %.3s %c%4d    ",
		     count++, ElemDesc[aptr->refno], Residue[group->refno],
		     chain->ident, group->serno );

	    x = (double)aptr->xorg/250.0;
	    y = (double)aptr->yorg/250.0;
	    z = (double)aptr->zorg/250.0;

#ifdef INVERT
	    fprintf(DataFile,"%8.3f%8.3f%8.3f",x,-y,-z);
#else
	    fprintf(DataFile,"%8.3f%8.3f%8.3f",x,y,-z);
#endif
	    fprintf(DataFile,"  1.00%6.2f\n",aptr->temp/100.0);
	    
	    ch = chain->ident;
	    prev = group;
	}

    if( prev )
	fprintf( DataFile, "TER   %5d      %.3s %c%4d \n", 
		 count, Residue[prev->refno], ch, prev->serno);

    fputs("END   \n",DataFile);
    fclose( DataFile );
#ifdef APPLEMAC
    SetFileInfo(filename,'RSML','TEXT',131);
#endif
    return( True );
}

int SaveXYZMolecule( filename )
    char *filename;
{
    return( True );
}


int SaveCIFMolecule( filename )
    char *filename;
{
    if( !filename )
	    return( False );
    return( True );
}


int SaveAlchemyMolecule( filename )
    char *filename;
{
    register Real x, y, z;
    register float xpos, ypos, zpos;
    register Chain __far *chain;
    register Group __far *group;
    register Atom __far *aptr;
    register Bond __far *bptr;
    register char *ptr;
    register int atomno;
    register int bondno;
    register int num;

    if( !Database )
	return( False );

    DataFile = fopen( filename, "w" );
    if( !DataFile )
    {   if( CommandActive )
	    WriteChar('\n');
	WriteString("Error: Unable to create file!\n\n");
	CommandActive=False;
	return( False );
    }

    atomno = 0;
    ForEachAtom
	if( aptr->flag & SelectFlag )
	{   aptr->mbox = 0;
	    atomno++;
	}

    bondno = 0;
    ForEachBond
	if( (bptr->srcatom->flag&bptr->dstatom->flag) & SelectFlag )
	{   if( bptr->flag&AromBondFlag )
	    {   bptr->srcatom->mbox = -1;
		bptr->dstatom->mbox = -1;
	    } else if( !(bptr->flag&HydrBondFlag) )
	    {   num = (bptr->flag&DoubBondFlag)? 2 : 1;
		if( bptr->srcatom->mbox>0 ) 
		    bptr->srcatom->mbox += num;
		if( bptr->dstatom->mbox>0 ) 
		    bptr->dstatom->mbox += num;
	    }
	    bondno++;
	}

    fprintf(DataFile,"%5d ATOMS, %5d BONDS, ",atomno,bondno);
    fprintf(DataFile,"    0 CHARGES, %s\n", InfoMoleculeName );

    atomno = 1;
    ForEachAtom
	if( aptr->flag & SelectFlag )
	{   aptr->mbox = atomno;
	    fprintf(DataFile,"%5d ",atomno++);

	    switch( GetElemNumber(aptr) )
	    {   case( 6 ):  if( aptr->mbox == -1 )
			    {   ptr = "CAR ";
			    } else if( aptr->mbox == 1 )
			    {   ptr = "C3  ";
			    } else ptr = "C2  ";
			    fputs( ptr, DataFile );
			    break;

		case( 7 ):  if( aptr->mbox == -1 )
			    {   ptr = "NAR ";
			    } else ptr = "N2  ";
			    fputs( ptr, DataFile );
			    break;

		case( 8 ):  if( aptr->mbox == 2 )
			    {   ptr = "O2  ";
			    } else ptr = "O3  ";
			    fputs( ptr, DataFile );
			    break;

		case( 1 ):  fputs( "H   ", DataFile );  break;

		default:    ptr = ElemDesc[aptr->refno];
			    if( *ptr==' ' )
			    {   fprintf(DataFile,"%.3s ",ptr+1);
			    } else fprintf(DataFile,"%.4s",ptr);
	    }

	    x = aptr->xorg/250.0;
	    y = aptr->yorg/250.0;
	    z = aptr->zorg/250.0;

	    /* Apply Current Viewpoint Rotation Matrix */
	    xpos = (float)(x*RotX[0] + y*RotX[1] + z*RotX[2]);
	    ypos = (float)(x*RotY[0] + y*RotY[1] + z*RotY[2]);
	    zpos = (float)(x*RotZ[0] + y*RotZ[1] + z*RotZ[2]);

#ifdef INVERT
	    fprintf(DataFile,"  %8.4f %8.4f %8.4f",xpos,-ypos,-zpos);
#else
	    fprintf(DataFile,"  %8.4f %8.4f %8.4f",xpos,ypos,-zpos);
#endif
	    fprintf(DataFile,"    %7.4f\n",aptr->temp/1000.0);
	}

    bondno = 1;
    ForEachBond
	if( (bptr->srcatom->flag&bptr->dstatom->flag) & SelectFlag )
	{   fprintf(DataFile,"%5d %5d %5d  ", bondno++,
			bptr->srcatom->mbox, bptr->dstatom->mbox );
	    if( bptr->flag & AromBondFlag )
	    {   ptr = "AROMATIC\n";
	    } else if( bptr->flag & TripBondFlag )
	    {   ptr = "TRIPLE\n";
	    } else if( bptr->flag & DoubBondFlag )
	    {   ptr = "DOUBLE\n";
	    } else ptr = "SINGLE\n";
	    fputs( ptr, DataFile );
	}

    ForEachAtom
	    if( aptr->flag & SelectFlag )
	        aptr->mbox = 0;
    fclose( DataFile );
#ifdef APPLEMAC
    SetFileInfo(filename,'RSML','TEXT',131);
#endif
    return( True );
}


static void TestBonded( sptr, dptr, flag )
    Atom __far *sptr, __far *dptr; 
    int flag;
{
    register Bond __far *bptr;
    register Long dx, dy, dz;
    register Long max, dist;

    if( flag )
    {    /* Sum of covalent radii with 0.5A tolerance */
         dist = Element[GetElemNumber(sptr)].covalrad + 
                Element[GetElemNumber(dptr)].covalrad + 125;
         max = dist*dist;  
    } else 
    {    /* Fast Bio-Macromolecule Bonding Calculation */
         if( (sptr->flag|dptr->flag) & HydrogenFlag )
	 {      max = MaxHBondDist;
         } else max = MaxBondDist;
    }

    dx = sptr->xorg-dptr->xorg;   if( (dist=dx*dx)>max ) return;
    dy = sptr->yorg-dptr->yorg;   if( (dist+=dy*dy)>max ) return;
    dz = sptr->zorg-dptr->zorg;   if( (dist+=dz*dz)>max ) return;

    if( dist > MinBondDist )
    {   /* Reset Non-bonded flags! */
	sptr->flag &= ~NonBondFlag;
	dptr->flag &= ~NonBondFlag;

        if( (sptr->flag|dptr->flag) & HydrogenFlag )
        {      flag = HydrBondFlag;
        } else flag = NormBondFlag;
	bptr = ProcessBond(sptr,dptr,flag);
	bptr->bnext = CurMolecule->blist;
	CurMolecule->blist = bptr;
	InfoBondCount++;
    }
}


static void ReclaimHBonds( ptr )
    HBond __far *ptr;
{
    register HBond __far *temp;

    if( (temp = ptr) )
    {   while( temp->hnext )
	    temp=temp->hnext;
	temp->hnext = FreeHBond;
	FreeHBond = ptr;
    }
}


static void ReclaimBonds( ptr )
    Bond __far *ptr;
{
    register Bond __far *temp;

    if( (temp = ptr) )
    {   while( temp->bnext )
	    temp=temp->bnext;
	temp->bnext = FreeBond;
	FreeBond = ptr;
    }
}


void CreateMoleculeBonds( info, flag )
    int info, flag;
{
    register int i, x, y, z;
    register Long tx, ty, tz;
    register Long mx, my, mz; 
    register Long dx, dy, dz;
    register int lx, ly, lz, hx, hy, hz;
    register Atom __far *aptr, __far *dptr;
    register Chain __far *chain;
    register Group __far *group;
    char buffer[40];


    if( !Database ) 
	return;

    dx = (MaxX-MinX)+1;
    dy = (MaxY-MinY)+1;
    dz = (MaxZ-MinZ)+1;

    InfoBondCount = 0;
    ReclaimBonds( CurMolecule->blist );
    CurMolecule->blist = (void __far*)0;
    ResetVoxelData();

    for( chain=Database->clist; chain; chain=chain->cnext )
    {   ResetVoxelData();
	for( group=chain->glist; group; group=group->gnext )
	    for( aptr=group->alist; aptr; aptr=aptr->anext )
	    {   /* Initially non-bonded! */
		aptr->flag |= NonBondFlag;

		mx = aptr->xorg-MinX;
		my = aptr->yorg-MinY;
		mz = aptr->zorg-MinZ;

		tx = mx-AbsMaxBondDist;  
		ty = my-AbsMaxBondDist;  
		tz = mz-AbsMaxBondDist;  

		lx = (tx>0)? (int)((VOXORDER*tx)/dx) : 0;
		ly = (ty>0)? (int)((VOXORDER*ty)/dy) : 0;
		lz = (tz>0)? (int)((VOXORDER*tz)/dz) : 0;

		tx = mx+AbsMaxBondDist;  
		ty = my+AbsMaxBondDist;  
		tz = mz+AbsMaxBondDist;

		hx = (tx<dx)? (int)((VOXORDER*tx)/dx) : VOXORDER-1;
		hy = (ty<dy)? (int)((VOXORDER*ty)/dy) : VOXORDER-1;
		hz = (tz<dz)? (int)((VOXORDER*tz)/dz) : VOXORDER-1;
	
		for( x=lx; x<=hx; x++ )
		{   i = VOXORDER2*x + VOXORDER*ly;
		    for( y=ly; y<=hy; y++ )
		    {   for( z=lz; z<=hz; z++ )
			    if( (dptr = (Atom __far*)HashTable[i+z]) )
				do { TestBonded(aptr,dptr,flag);
				} while( (dptr = dptr->next) );
			i += VOXORDER;
		    }
		}
		
		x = (int)((VOXORDER*mx)/dx);
		y = (int)((VOXORDER*my)/dy);
		z = (int)((VOXORDER*mz)/dz);

		i = VOXORDER2*x + VOXORDER*y + z;
		aptr->next = (Atom __far*)HashTable[i];
		HashTable[i] = (void __far*)aptr;
	    }
	VoxelsClean = False;
    }

    if( info )
    {   if( CommandActive )
	    WriteChar('\n');
	CommandActive=False;
	sprintf(buffer,"Number of Bonds ..... %ld\n\n",InfoBondCount);
	WriteString(buffer);
    }
}


Atom __far *FindGroupAtom( group, n )
    Group __far *group;  Byte n;
{
    register Atom __far *ptr;

    ptr = group->alist;
    while( ptr )
    {   if( ptr->refno == n )
	    return( ptr );
	ptr = ptr->anext;
    }
    return( (Atom __far*)0 );
}


void TestDisulphideBridge( group1, group2, cys1 )
    Group __far *group1, __far *group2;  
    Atom __far *cys1;
{
    register HBond __far *ptr;
    register Atom __far *cys2;
    register int dx, dy, dz;
    register Long max,dist;

    if( !(cys2=FindGroupAtom(group2,20)) )
	return;

    max = (Long)750*750;
    dx = (int)(cys1->xorg-cys2->xorg);   if( (dist=(Long)dx*dx)>max ) return;
    dy = (int)(cys1->yorg-cys2->yorg);   if( (dist+=(Long)dy*dy)>max ) return;
    dz = (int)(cys1->zorg-cys2->zorg);   if( (dist+=(Long)dz*dz)>max ) return;

    if( !(ptr = FreeHBond) )
    {   MemSize += sizeof(HBond);
	ptr = (HBond __far*)_fmalloc(sizeof(HBond));
	if( !ptr ) FatalDataError("Memory allocation failed");
	RegisterAlloc( ptr );
    } else FreeHBond = ptr->hnext;

    ptr->dst = cys1;
    if( !(ptr->dstCA=FindGroupAtom(group1,1)) )
	ptr->dstCA = cys1;

    ptr->src = cys2;
    if( !(ptr->srcCA=FindGroupAtom(group2,1)) )
	ptr->srcCA = cys2;

    ptr->offset = 0;
    ptr->energy = 0;
    ptr->flag = 0;
    ptr->col = 0;

    ptr->hnext = CurMolecule->slist;
    CurMolecule->slist = ptr;

    group1->flag |= CystineFlag;
    group2->flag |= CystineFlag;
    InfoSSBondCount++;
}


void FindDisulphideBridges()
{
    register Chain __far *chn1;
    register Chain __far *chn2;
    register Group __far *group1;
    register Group __far *group2;
    register Atom __far *cys;
    char buffer[40];

    if( !Database ) return;
    ReclaimHBonds( CurMolecule->slist );
    InfoSSBondCount = 0;

    for(chn1=Database->clist;chn1;chn1=chn1->cnext)
	for(group1=chn1->glist;group1;group1=group1->gnext)
	    if( IsCysteine(group1->refno) && (cys=FindGroupAtom(group1,20)) )
	    {   for(group2=group1->gnext;group2;group2=group2->gnext)
		    if( IsCysteine(group2->refno) )
			TestDisulphideBridge(group1,group2,cys);

		for(chn2=chn1->cnext;chn2;chn2=chn2->cnext)
		    for(group2=chn2->glist;group2;group2=group2->gnext)
			if( IsCysteine(group2->refno) )
			    TestDisulphideBridge(group1,group2,cys);
	    }

    if( CommandActive )
	WriteChar('\n');
    CommandActive=False;
    
    sprintf(buffer,"Number of Bridges ... %d\n\n",InfoSSBondCount);
    WriteString(buffer);
}


#ifdef FUNCPROTO
static void CreateHydrogenBond( Atom __far*, Atom __far*,
				Atom __far*, Atom __far*, int, int );
static int IsHBonded( Atom __far*, Atom __far*, HBond __far* );
static void TestLadder( Chain __far* );
#endif


static void CreateHydrogenBond( srcCA, dstCA, src, dst, energy, offset )
    Atom __far *srcCA, __far *dstCA;
    Atom __far *src, __far *dst;
    int energy, offset;
{
    register HBond __far *ptr;
    register int i,flag;

    if( !(ptr = FreeHBond) )
    {   MemSize += HBondPool*sizeof(HBond);
	ptr = (HBond __far *)_fmalloc( HBondPool*sizeof(HBond) );
	if( !ptr ) FatalDataError("Memory allocation failed");
	RegisterAlloc( ptr );
	for( i=1; i<HBondPool; i++ )
	{   ptr->hnext = FreeHBond;
	    FreeHBond = ptr++;
	} 
    } else FreeHBond = ptr->hnext;

    if( (offset>=-128) && (offset<127) )
    {   ptr->offset = (Char)offset;
    } else ptr->offset = 0;

    flag = ZoneBoth? src->flag&dst->flag : src->flag|dst->flag;
    ptr->flag = flag & SelectFlag;

    ptr->src = src;
    ptr->dst = dst;
    ptr->srcCA = srcCA;
    ptr->dstCA = dstCA;
    ptr->energy = energy;
    ptr->col = 0;

    *CurHBond = ptr;
    ptr->hnext = (void __far*)0;
    CurHBond = &ptr->hnext;
    InfoHBondCount++;
}

/* Coupling constant for Electrostatic Energy   */
/* QConst = -332 * 0.42 * 0.2 * 1000.0 [*250.0] */
#define QConst (-6972000.0)
#define MaxHDist ((Long)2250*2250)
#define MinHDist ((Long)125*125)


/* Protein Donor Atom Coordinates */
static int hxorg,hyorg,hzorg;
static int nxorg,nyorg,nzorg;
static Atom __far *best1CA;
static Atom __far *best2CA;
static Atom __far *best1;
static Atom __far *best2;
static Atom __far *optr;
static int res1,res2;
static int off1,off2;


static int CalculateBondEnergy( group )
    Group __far *group;
{
    register double dho,dhc;
    register double dnc,dno;

    register Atom __far *cptr;
    register Long dx,dy,dz;
    register Long dist;
    register int result;

    if( !(cptr=FindGroupAtom(group,2)) )  return(0);
    if( !(optr=FindGroupAtom(group,3)) )  return(0);

    dx = hxorg-optr->xorg;  
    dy = hyorg-optr->yorg;  
    dz = hzorg-optr->zorg;
    dist = dx*dx+dy*dy+dz*dz;
    if( dist < MinHDist ) 
	return( -9900 );
    dho = sqrt((double)dist);

    dx = hxorg-cptr->xorg;  
    dy = hyorg-cptr->yorg;  
    dz = hzorg-cptr->zorg;
    dist = dx*dx+dy*dy+dz*dz;
    if( dist < MinHDist ) 
	return( -9900 );
    dhc = sqrt((double)dist);

    dx = nxorg-cptr->xorg;  
    dy = nyorg-cptr->yorg;  
    dz = nzorg-cptr->zorg;
    dist = dx*dx+dy*dy+dz*dz;
    if( dist < MinHDist ) 
	return( -9900 );
    dnc = sqrt((double)dist);

    dx = nxorg-optr->xorg;  
    dy = nyorg-optr->yorg;  
    dz = nzorg-optr->zorg;
    dist = dx*dx+dy*dy+dz*dz;
    if( dist < MinHDist ) 
	return( -9900 );
    dno = sqrt((double)dist);

    result = (int)(QConst/dho - QConst/dhc + QConst/dnc - QConst/dno);

    if( result<-9900 ) 
    {   return( -9900 );
    } else if( result>-500 ) 
	return( 0 );

    return( result );

}


void CalcHydrogenBonds()
{
    register int energy, offset, refno;
    register Chain __far *chn1, __far *chn2;
    register Group __far *group1;
    register Group __far *group2;
    register Group __far *best;
    register Atom __far *ca1;
    register Atom __far *ca2;
    register Atom __far *pc1;
    register Atom __far *po1;
    register Atom __far *n1;
    register Long max,dist;
    register int pos1,pos2;
    register int dx,dy,dz;
    register double dco;
    char buffer[40];

    if( !Database ) return;
    ReclaimHBonds( CurMolecule->hlist );
    CurMolecule->hlist = (void __far*)0;
    CurHBond = &CurMolecule->hlist;
    InfoHBondCount = 0;

    if( MainAtomCount > 10000 )
    {   if( CommandActive )
	    WriteChar('\n');
	WriteString("Please wait... ");
	CommandActive=True;
    }

    for(chn1=Database->clist;chn1;chn1=chn1->cnext)
    {   if( !chn1->glist ) continue;

	if( IsProtein(chn1->glist->refno) )
	{   pc1 = po1 = (void __far*)0;
	    pos1 = 0;
	    for(group1=chn1->glist;group1;group1=group1->gnext)
	    {   pos1++;
		if( pc1 && po1 )
		{   dx = (int)(pc1->xorg - po1->xorg);
		    dy = (int)(pc1->yorg - po1->yorg);
		    dz = (int)(pc1->zorg - po1->zorg);
		} else
		{   pc1 = FindGroupAtom(group1,2);
		    po1 = FindGroupAtom(group1,3);
		    continue;
		}

		pc1 = FindGroupAtom(group1,2);
		po1 = FindGroupAtom(group1,3);

		if( !IsAmino(group1->refno) || IsProline(group1->refno) )
		    continue;

		if( !(ca1=FindGroupAtom(group1,1)) ) continue;
		if( !(n1=FindGroupAtom(group1,0)) )  continue;

		dist = (Long)dx*dx + (Long)dy*dy + (Long)dz*dz;
		dco = sqrt( (double)dist )/250.0;

		nxorg = (int)n1->xorg;   hxorg = nxorg + (int)(dx/dco);
		nyorg = (int)n1->yorg;   hyorg = nyorg + (int)(dy/dco);
		nzorg = (int)n1->zorg;   hzorg = nzorg + (int)(dz/dco);
		res1 = res2 = 0;

		/* Only Hydrogen Bond within a single chain!       */
		/* for(chn2=Database->clist;chn2;chn2=chn2->cnext) */

		chn2 = chn1;
		{   /* Only consider non-empty peptide chains! */
		    if( !chn2->glist || !IsProtein(chn2->glist->refno) )
			continue;

		    pos2 = 0;
		    for(group2=chn2->glist;group2;group2=group2->gnext)
		    {   pos2++;
			if( (group2==group1) || (group2->gnext==group1) )
			    continue;

			if( !IsAmino(group2->refno) ) 
			    continue;
			if( !(ca2=FindGroupAtom(group2,1)) ) 
			    continue;

			dx = (int)(ca1->xorg-ca2->xorg);
			if( (dist=(Long)dx*dx) > MaxHDist )
			    continue;

			dy = (int)(ca1->yorg-ca2->yorg);
			if( (dist+=(Long)dy*dy) > MaxHDist )
			    continue;

			dz = (int)(ca1->zorg-ca2->zorg);
			if( (dist+=(Long)dz*dz) > MaxHDist )
			    continue;

			if( (energy = CalculateBondEnergy(group2)) )
			{   if( chn1 == chn2 )
			    {   offset = pos1 - pos2;
			    } else offset = 0;

			    if( energy<res1 )
			    {   best2CA = best1CA;  best1CA = ca2;
				best2 = best1;      best1 = optr;
				res2 = res1;        res1 = energy;
				off2 = off1;        off1 = offset;
			    } else if( energy<res2 )
			    {   best2CA = ca2;
				best2 = optr;
				res2 = energy;
				off2 = offset;
			    }
			}
		    }  /* group2 */
		}      /* chn2 */

		if( res1 ) 
		{   if( res2 ) 
			CreateHydrogenBond(ca1,best2CA,n1,best2,res2,off2);
		    CreateHydrogenBond(ca1,best1CA,n1,best1,res1,off1);
		}
	    }
	} else if( IsNucleo(chn1->glist->refno) )
	    for(group1=chn1->glist;group1;group1=group1->gnext)
	    {   if( !IsPurine(group1->refno) ) continue;
		/* Find N1 of Purine Group */
		if( !(n1=FindGroupAtom(group1,21)) )
		    continue;

		/* Maximum N1-N3 distance 5A */
		refno = NucleicCompl(group1->refno);
		max = (Long)1250*1250;
		best = (void __far*)0;

		for(chn2=Database->clist;chn2;chn2=chn2->cnext)
		{   /* Only consider non-empty peptide chains! */
		    if( (chn1==chn2) || !chn2->glist || 
			!IsNucleo(chn2->glist->refno) )
			continue;

		    for(group2=chn2->glist;group2;group2=group2->gnext)
			if( group2->refno == (Byte)refno )
			{   /* Find N3 of Pyramidine Group */
			    if( !(ca1=FindGroupAtom(group2,23)) )
				continue;

			    dx = (int)(ca1->xorg - n1->xorg);
			    if( (dist=(Long)dx*dx) >= max ) 
				continue;

			    dy = (int)(ca1->yorg - n1->yorg);
			    if( (dist+=(Long)dy*dy) >= max ) 
				continue;

			    dz = (int)(ca1->zorg - n1->zorg);
			    if( (dist+=(Long)dz*dz) >= max )
				continue;

			    best1 = ca1;
			    best = group2;
			    max = dist;
			}
		}

		if( best )
		{   /* Find the sugar phosphorous atoms */
		    ca1 = FindGroupAtom( group1, 7 );
		    ca2 = FindGroupAtom( best, 7 );

		    CreateHydrogenBond( ca1, ca2, n1, best1, 0, 0 );
		    if( IsGuanine(group1->refno) )
		    {   /* Guanine-Cytosine */
			if( (ca1=FindGroupAtom(group1,22)) &&  /* G.N2 */
			    (ca2=FindGroupAtom(best,26)) )     /* C.O2 */
			    CreateHydrogenBond( (void __far*)0, (void __far*)0,
						ca1, ca2, 0, 0 );

			if( (ca1=FindGroupAtom(group1,28)) &&  /* G.O6 */
			    (ca2=FindGroupAtom(best,24)) )     /* C.N4 */
			    CreateHydrogenBond( (void __far*)0, (void __far*)0,
						ca1, ca2, 0, 0 );

		    } else /* Adenine-Thymine */
			if( (ca1=FindGroupAtom(group1,25)) &&  /* A.N6 */
			    (ca2=FindGroupAtom(best,27)) )     /* T.O4 */
			    CreateHydrogenBond( (void __far*)0, (void __far*)0,
						ca1, ca2, 0, 0 );
		}
	    }
    }

    if( CommandActive )
	WriteChar('\n');
    CommandActive=False;
    
    sprintf(buffer,"Number of H-Bonds ... %d\n",InfoHBondCount);
    WriteString(buffer);
}


static int IsHBonded( src, dst, ptr )
    Atom __far *src, __far *dst;
    HBond __far *ptr;
{
    while( ptr && (ptr->srcCA==src) )
	if( ptr->dstCA == dst )
	{   return( True );
	} else ptr=ptr->hnext;
    return( False );
}


static void FindAlphaHelix( pitch, flag )
    int pitch, flag;
{
    register HBond __far *hbond;
    register Chain __far *chain;
    register Group __far *group;
    register Group __far *first;
    register Group __far *ptr;
    register Atom __far *srcCA;
    register Atom __far *dstCA;
    register int res,dist,prev;

    /* Protein chains only! */
    hbond = Database->hlist;
    for( chain=Database->clist; chain; chain=chain->cnext )
    if( (first=chain->glist) && IsProtein(first->refno) )
    {   prev = False; dist = 0;
	for( group=chain->glist; hbond && group; group=group->gnext )
	{   if( IsAmino(group->refno) && (srcCA=FindGroupAtom(group,1)) )
	    {   if( dist==pitch )
		{   res = False;
		    dstCA=FindGroupAtom(first,1);

		    while( hbond && hbond->srcCA == srcCA )
		    {   if( hbond->dstCA==dstCA ) res=True;
			hbond = hbond->hnext;
		    }

		    if( res )
		    {   if( prev )
			{   if( !(first->struc&HelixFlag) ) 
				InfoHelixCount++;

			    ptr = first;
			    do {
				ptr->struc |= flag;
				ptr = ptr->gnext;
			    } while( ptr != group );
			} else prev = True;
		    } else prev = False;
		} else while( hbond && hbond->srcCA==srcCA )
		    hbond = hbond->hnext;
	    } else prev = False;

	    if( group->struc&HelixFlag )
	    {   first = group; prev = False; dist = 1;
	    } else if( dist==pitch )
	    {   first = first->gnext;
	    } else dist++;
	}
    } else if( first && IsNucleo(first->refno) )
	while( hbond && !IsAminoBackbone(hbond->src->refno) )
	    hbond = hbond->hnext;
}


static Atom __far *cprevi, __far *ccurri, __far *cnexti;
static HBond __far *hcurri, __far *hnexti;
static Group __far *curri, __far *nexti;



static void TestLadder( chain )
    Chain __far *chain;
{
    register Atom __far *cprevj, __far *ccurrj, __far *cnextj;
    register HBond __far *hcurrj, __far *hnextj;
    register Group __far *currj, __far *nextj;
    register int count, result, found;

    /* Already part of atleast one ladder */
    found = curri->flag & SheetFlag;
    nextj = nexti->gnext;

    hnextj = hnexti;
    while( hnextj && hnextj->srcCA==cnexti )
	hnextj = hnextj->hnext;

    while( True )
    {   if( nextj )
	    if( IsProtein(chain->glist->refno) )
	    {   count = 1;
		do {
		    cnextj = FindGroupAtom(nextj,1);
		    if( count == 3 )
		    {   if( IsHBonded(cnexti,ccurrj,hnexti) &&
			    IsHBonded(ccurrj,cprevi,hcurrj) )
			{   result = ParaLadder;
			} else if( IsHBonded(cnextj,ccurri,hnextj) &&
				   IsHBonded(ccurri,cprevj,hcurri) )
			{   result = ParaLadder;
			} else if( IsHBonded(cnexti,cprevj,hnexti) &&
				   IsHBonded(cnextj,cprevi,hnextj) )
			{   result = AntiLadder;
			} else if( IsHBonded(ccurrj,ccurri,hcurrj) &&
				   IsHBonded(ccurri,ccurrj,hcurri) )
			{   result = AntiLadder;
			} else result = NoLadder;

			if( result )
			{   curri->struc |= SheetFlag;
			    currj->struc |= SheetFlag;
			    if( found ) return;
			    found = True;
			}
		    } else count++;

		    cprevj = ccurrj; ccurrj = cnextj; 
		    currj = nextj;   hcurrj = hnextj;

		    while( hnextj && hnextj->srcCA==cnextj )
			hnextj = hnextj->hnext;
		} while( (nextj = nextj->gnext) );

	    } else if( IsNucleo(chain->glist->refno) )
		while( hnextj && !IsAminoBackbone(hnextj->src->refno) )
		    hnextj = hnextj->hnext;

	if( (chain = chain->cnext) ) 
	{   nextj = chain->glist;
	} else return;
    }
}


static void FindBetaSheets()
{
    register Chain __far *chain;
    register int ladder;
    register int count;

    hnexti = Database->hlist;
    for( chain=Database->clist; chain; chain=chain->cnext )
	if( (nexti = chain->glist) )
	    if( IsProtein(nexti->refno) )
	    {   count = 1;
		ladder = False;
		do {
		    cnexti = FindGroupAtom(nexti,1);

		    if( count == 3 )
		    {   TestLadder( chain );
			if( curri->struc & SheetFlag )
			{   if( !ladder )
			    {   InfoLadderCount++;
				ladder = True;
			    }
			} else ladder = False;
		    } else count++;

		    cprevi = ccurri; ccurri = cnexti; 
		    curri = nexti;   hcurri = hnexti;
		    while( hnexti && hnexti->srcCA==cnexti )
			hnexti = hnexti->hnext;
		} while( (nexti = nexti->gnext) );

	    } else if( IsNucleo(nexti->refno) )
		while( hnexti && !IsAminoBackbone(hnexti->src->refno) )
		    hnexti = hnexti->hnext;
}


static void FindTurnStructure()
{
    static Atom __far *aptr[5];
    register Chain __far *chain;
    register Group __far *group;
    register Group __far *prev;
    register Atom __far *ptr;
    register Long ux,uy,uz,mu;
    register Long vx,vy,vz,mv;
    register int i,found,len;
    register Real CosKappa;

    for( chain=Database->clist; chain; chain=chain->cnext )
	if( chain->glist && IsProtein(chain->glist->refno) )
	{   len = 0;  found = False;
	    for( group=chain->glist; group; group=group->gnext )
	    {    ptr = FindGroupAtom(group,1);
		 if( ptr && (ptr->flag&BreakFlag) )
		 {   found = False;
		     len = 0;
		 } else if( len==5 )
		 {   for( i=0; i<4; i++ )
			 aptr[i] = aptr[i+1];
		     len = 4;
		 } else if( len==2 )
		     prev = group;

		 aptr[len++] = ptr;
		 if( len==5 ) 
		 {   if( !(prev->struc&(HelixFlag|SheetFlag)) &&
			 aptr[0] && aptr[2] && aptr[4] )
		     {   ux = aptr[2]->xorg - aptr[0]->xorg;
			 uy = aptr[2]->yorg - aptr[0]->yorg;
			 uz = aptr[2]->zorg - aptr[0]->zorg;

			 vx = aptr[4]->xorg - aptr[2]->xorg;
			 vy = aptr[4]->yorg - aptr[2]->yorg;
			 vz = aptr[4]->zorg - aptr[2]->zorg;

			 mu = ux*ux + uy*uy + uz*uz;
			 mv = vx*vx + vz*vz + vy*vy;
			 if( mu && mv )
			 {   CosKappa = (Real)(ux*vx + uy*vy + uz*vz);
			     CosKappa /= sqrt( (Real)mu*mv );
			     if( CosKappa<Cos70Deg )
			     {   if( !found )
				     InfoTurnCount++;
				 prev->struc |= TurnFlag;
			     }
			 }
		     }
		     found = prev->struc&TurnFlag;
		     prev = prev->gnext;
		 } /* len==5 */
	    }
	}
}

void DetermineStructure()
{
    register Chain __far *chain;
    register Group __far *group;
    char buffer[40];

    if( !Database )
	return;

    if( InfoHBondCount<0 )
	CalcHydrogenBonds();

    if( InfoHelixCount>=0 )
	for( chain=Database->clist; chain; chain=chain->cnext )
	    for( group=chain->glist; group; group=group->gnext )
		group->struc = 0;

    InfoStrucSource = SourceCalc;
    InfoLadderCount = 0;
    InfoHelixCount = 0;
    InfoTurnCount = 0;

    if( InfoHBondCount )
    {   FindAlphaHelix(4,Helix4Flag);
	FindBetaSheets();
	FindAlphaHelix(3,Helix3Flag);
	FindAlphaHelix(5,Helix5Flag);
	FindTurnStructure();
    }

    if( CommandActive )
	WriteChar('\n');
    CommandActive=False;

    sprintf(buffer,"Number of Helices ... %d\n",InfoHelixCount);
    WriteString(buffer);
    sprintf(buffer,"Number of Strands ... %d\n",InfoLadderCount);
    WriteString(buffer);
    sprintf(buffer,"Number of Turns ..... %d\n",InfoTurnCount);
    WriteString(buffer);
}


void RenumberMolecule( start )
    int start;
{
    register Chain __far *chain;
    register Group __far *group;
    register int hinit, minit;
    register int resno;

    if( !Database )
	return;

    hinit = minit = False;
    for( chain=Database->clist; chain; chain=chain->cnext )
    {   resno = start;
	for( group=chain->glist; group; group=group->gnext )
	{   if( group->alist->flag & HeteroFlag )
	    {   if( hinit )
		{   if( resno > MaxHetaRes )
		    {   MaxHetaRes = resno;
		    } else if( resno < MinHetaRes )
			MinHetaRes = resno;
		} else MinHetaRes = MaxHetaRes = resno;
		hinit = True;
	    } else
	    {   if( minit )
		{   if( resno > MaxMainRes )
		    {   MaxMainRes = resno;
		    } else if( resno < MinMainRes )
			MinMainRes = resno;
		} else MinMainRes = MaxMainRes = resno;
		minit = True;
	    }
	    group->serno = resno++;
	}
    }
}


static void ReclaimAtoms( ptr )
    Atom __far *ptr;
{
    register Atom __far *temp;

    if( (temp = ptr) )
    {   while( temp->anext )
	    temp=temp->anext;
	temp->anext = FreeAtom;
	FreeAtom = ptr;
    }
}

static void ResetDatabase()
{
    
    Database = CurMolecule = (void __far*)0;
    MainGroupCount = HetaGroupCount = 0;
    InfoChainCount = HetaAtomCount = 0;
    MainAtomCount = InfoBondCount = 0;  
    SelectCount = 0;

    InfoStrucSource = SourceNone;
    InfoSSBondCount = InfoHBondCount = -1;
    InfoHelixCount = InfoLadderCount = -1;
    InfoTurnCount = -1;

    CurGroup = (void __far*)0;
    CurChain = (void __far*)0;
    CurAtom = (void __far*)0;

    MinX = MinY = MinZ = 0;
    MaxX = MaxY = MaxZ = 0;

    MinMainTemp = MaxMainTemp = 0;
    MinHetaTemp = MaxHetaTemp = 0;
    MinMainRes = MaxMainRes = 0;
    MinHetaRes = MaxHetaRes = 0;

    *InfoMoleculeName = 0;
    *InfoClassification = 0;
    *InfoIdentCode = 0;
    *InfoSpaceGroup = 0;
    *InfoFileName = 0;

    VoxelsClean = False;
    HMinMaxFlag = False;
    MMinMaxFlag = False;
    HasHydrogen = False;
    ElemNo = MINELEM;
    ResNo = MINRES;
    MaskCount = 0;
}


void DestroyDatabase()
{
    register void __far *temp;
    register Group __far *gptr;

    if( Database )
    {   ReclaimHBonds( Database->slist );
	ReclaimHBonds( Database->hlist );
	ReclaimBonds( Database->blist );

	while( Database->clist )
	{   ReclaimBonds(Database->clist->blist);
	    if( (gptr = Database->clist->glist) )
	    {   ReclaimAtoms(gptr->alist);
		while( gptr->gnext )
		{   gptr = gptr->gnext;
		    ReclaimAtoms(gptr->alist);
		}
		gptr->gnext = FreeGroup;
		FreeGroup = Database->clist->glist;
	    }
	    temp = Database->clist->cnext;
	    Database->clist->cnext = FreeChain;
	    FreeChain = Database->clist;
	    Database->clist = temp;
	}

	FreeMolecule = Database;
	Database = (void __far*)0;
    }
    ResetDatabase();
}


void PurgeDatabase()
{
#ifdef APPLEMAC
    register AllocRef *ptr;
    register AllocRef *tmp;
    register int i;
    
    /* Avoid Memory Leaks */
    for( ptr=AllocList; ptr; ptr=tmp )
    {   for( i=0; i<ptr->count; i++ )
	    _ffree( ptr->data[i] );
	tmp = ptr->next;
	_ffree( ptr );
    }
#endif
}


void InitialiseDatabase()
{
    FreeMolecule = (void __far*)0;
    FreeHBond = (void __far*)0;
    FreeChain = (void __far*)0;
    FreeGroup = (void __far*)0;
    FreeAtom = (void __far*)0;
    FreeBond = (void __far*)0;

#ifdef APPLEMAC
    AllocList = (void*)0;
#endif

    ResetDatabase();
}

