/* molecule.c
 * RasMol2 Molecular Graphics
 * Roger Sayle, February 1994
 * Version 2.3
 */
#include "rasmol.h"

#ifdef IBMPC
#include <windows.h>
#include <malloc.h>
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

typedef struct {
        int init, term;
        char chain;
        char type;
        } FeatEntry;

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

#define MAXALCATOM   20
static ConvTable AlcAtomTable[MAXALCATOM] = {
    { { 'C', '3', ' ', ' ' }, { ' ', 'C', '3', ' ' } },  /*  0 */
    { { 'O', '3', ' ', ' ' }, { ' ', 'O', '3', ' ' } },  /*  1 */
    { { 'N', '3', ' ', ' ' }, { ' ', 'N', '3', ' ' } },  /*  2 */
    { { 'S', 'O', '2', ' ' }, { ' ', 'S', '2', ' ' } },  /*  3 */
    { { 'C', 'A', 'R', ' ' }, { ' ', 'C', ' ', ' ' } },  /*  4 */
    { { 'C', '2', ' ', ' ' }, { ' ', 'C', '2', ' ' } },  /*  5 */
    { { 'C', '1', ' ', ' ' }, { ' ', 'C', '1', ' ' } },  /*  6 */
    { { 'O', '2', ' ', ' ' }, { ' ', 'O', '2', ' ' } },  /*  7 */
    { { 'N', '2', ' ', ' ' }, { ' ', 'N', '2', ' ' } },  /*  8 */
    { { 'N', '1', ' ', ' ' }, { ' ', 'N', '1', ' ' } },  /*  9 */
    { { 'N', 'A', 'R', ' ' }, { ' ', 'N', ' ', ' ' } },  /* 10 */
    { { 'N', 'A', 'M', ' ' }, { ' ', 'N', ' ', ' ' } },  /* 11 */
    { { 'N', 'P', 'L', '3' }, { ' ', 'N', '3', ' ' } },  /* 12 */
    { { 'N', '3', '+', ' ' }, { ' ', 'N', '3', ' ' } },  /* 13 */
    { { 'P', '3', ' ', ' ' }, { ' ', 'P', '3', ' ' } },  /* 14 */
    { { 'P', '3', 'D', ' ' }, { ' ', 'P', '3', ' ' } },  /* 15 */
    { { 'P', '4', ' ', ' ' }, { ' ', 'P', '4', ' ' } },  /* 16 */
    { { 'S', '3', ' ', ' ' }, { ' ', 'S', '3', ' ' } },  /* 17 */
    { { 'S', '2', ' ', ' ' }, { ' ', 'S', '2', ' ' } },  /* 18 */
    { { 'S', 'O', ' ', ' ' }, { ' ', 'S', ' ', ' ' } },  /* 19 */
                                 };

static Molecule __far *FreeMolecule;
static HBond __far *FreeHBond;
static Chain __far *FreeChain;
static Group __far *FreeGroup;
static Atom __far *FreeAtom;
static Bond __far *FreeBond;

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
            if( (ch=='\n') || (ch=='\r') )
            {   if( ptr != Record+1 )
                {   *ptr = 0;
                    return;
                }
            } else if( ch==EOF )
            {   *ptr = 0;
                return;
            } else *ptr++ = ch;
        } while( ptr < Record+80 );

        /* skip to the end of the line! */
        do { ch = getc(DataFile);
        } while( (ch!='\n') && (ch!='\r') && (ch!=EOF) );
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

    if( Database )
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

        WriteString("Number of Ladders ... ");
        sprintf(buffer,"%d\n",InfoLadderCount);
        WriteString(buffer);

        WriteString("Number of Turns ..... ");
        sprintf(buffer,"%d\n",InfoTurnCount);
        WriteString(buffer);
    }
}


static void CreateChain( ident )
    char ident;
{
    register Chain __far *prev;

    if( !CurMolecule )
    {   if( !(CurMolecule = FreeMolecule) )
        {   MemSize += sizeof(Molecule);
            CurMolecule = (Molecule __far *)_fmalloc(sizeof(Molecule));
            if( !CurMolecule ) FatalDataError("Memory allocation failed");
        } else FreeMolecule = (void __far*)0;

        CurChain = (void __far*)0;
        CurMolecule->slist = (void __far*)0;
        CurMolecule->hlist = (void __far*)0;
        CurMolecule->blist = (void __far*)0;
        Database = CurMolecule;
    }

    prev = CurChain;
    if( !(CurChain = FreeChain) )
    {   MemSize += sizeof(Chain);
        CurChain = (Chain __far *)_fmalloc(sizeof(Chain));
        if( !CurChain ) FatalDataError("Memory allocation failed");
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


static Atom __far *CreateAtom()
{
    register Atom __far *ptr;
    register int i;

    if( !(ptr = FreeAtom) )
    {   MemSize += AtomPool*sizeof(Atom);
        ptr = (Atom __far *)_fmalloc( AtomPool*sizeof(Atom) );
        if( !ptr ) FatalDataError("Memory allocation failed");
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
    ptr->flag = SelectFlag;
    ptr->radius = 375;
    ptr->altl = ' ';
    ptr->mbox = 0;
    ptr->col = 0;

    return( ptr );
}


static void ProcessAtom( ptr )
    Atom __far *ptr;
{
    if( GetElemIdent(ptr) == HandH )
        ptr->flag |= HydrogenFlag;
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
               !strncmp(ptr,"H20",3) )
    /* check HHO, OHH, 0H2 and SOL? */
    {   return(31);
    } else if( !strncmp(ptr,"D20",3) )
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

    ptr->r = (int)(ReadValue(31,8)>>2) + 5;
    ptr->g = (int)(ReadValue(39,8)>>2) + 5;
    ptr->b = (int)(ReadValue(47,8)>>2) + 5;
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



static void ProcessPDBGroup( heta, serno )
    int heta, serno;
{
    register Group __far *ptr;
    register int i;

    if( !CurChain || (CurChain->ident!=Record[22]) )
        CreateChain( Record[22] );

    if( !(ptr = FreeGroup) )
    {   MemSize += GroupPool*sizeof(Group);
        ptr = (Group __far *)_fmalloc( GroupPool*sizeof(Group) );
        if( !ptr ) FatalDataError("Memory allocation failed");
        for( i=1; i<GroupPool; i++ )
        {   ptr->gnext = FreeGroup;
            FreeGroup = ptr++;
        } 
    } else FreeGroup = ptr->gnext;
    
    if( CurGroup )
    {   CurGroup->gnext = ptr;
    } else CurChain->glist = ptr;
    ptr->gnext = (void __far*)0;
    CurGroup = ptr;

    CurAtom = (void __far*)0;
    ptr->refno = FindResNo( &Record[18] );
    ptr->alist = (void __far*)0;
    ptr->serno = serno;
    ptr->flag = 0;
    ptr->col = 0;

    /* Solvents should be hetero! */
    if( IsSolvent(ptr->refno) )
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
            } else if( serno < MinHetaRes )
                MinHetaRes = serno;
        } else MinMainRes = MaxMainRes = serno; 
    }
}


static void ProcessPDBAtom( heta )
    int heta;
{
    register Bond __far *bptr;
    register Atom __far *ptr;
    register Long dx,dy,dz;
    register int serno;
    register int refno;
    register int i;

    /* Ignore Pseudo Atoms!! */
    if( (Record[13]==' ') && (Record[14]=='Q') )
        return; 

    serno = (int)ReadValue(23,4);
    if( !CurGroup || (CurGroup->serno!=serno) 
        || (CurChain->ident!=Record[22]) )
        ProcessPDBGroup( heta, serno );

    /* Solvents should be hetero! */
    if( IsSolvent(CurGroup->refno) )
        heta = True;


    ptr = CreateAtom();
    for( refno=0; refno<ElemNo; refno++ )
        if( !strncmp(ElemDesc[refno],Record+13,4) )
            break;

    if( refno==ElemNo )
    {   if( ElemNo++ == MAXELEM )
            FatalDataError("Too many new atom types");

        for( i=0; i<4; i++ )
            ElemDesc[refno][i] = Record[i+13];
    }

    ptr->refno = refno;
    ptr->altl = Record[17];
    ptr->serno = (int)ReadValue(7,5);
    ptr->temp = (int)ReadValue(61,6);

    ptr->xorg =  ReadValue(31,8)/4;
    ptr->yorg =  ReadValue(39,8)/4;
    ptr->zorg = -ReadValue(47,8)/4;

    if( heta ) 
        ptr->flag |= HeteroFlag;

    ProcessAtom( ptr );
    if( IsAlphaCarbon(refno) && IsAmino(CurGroup->refno) )
    {   if( ConnectAtom )
        {   dx = ConnectAtom->xorg - ptr->xorg;
            dy = ConnectAtom->yorg - ptr->yorg;
            dz = ConnectAtom->zorg - ptr->zorg;

            /* Break backbone if CA-CA > 7.00A */
            if( dx*dx+dy*dy+dz*dz < (Long)1750*1750 )
            {   bptr = ProcessBond(ConnectAtom,ptr,NormBondFlag);
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
        ptr->fnext = FeatList;
        ptr->count = 0;
        FeatList = ptr;
    } else ptr = FeatList;

    return( &(ptr->data[ptr->count++]) );
}


#if defined(__STDC__) || defined(IBMPC)
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
            {   group->flag |= mask;
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
            
        } else if( !strncmp("COLO",Record+1,4) )
        {   ProcessPDBColourMask();
        } else if( !strncmp("TER ",Record+1,4) )
        {   ConnectAtom = (void __far*)0;
        
        } else if( !strncmp("HEAD",Record+1,4) )
        {   ExtractString(40,Record+11,InfoClassification);
            ExtractString( 4,Record+63,InfoIdentCode);
            AdviseUpdate(AdvClass);
            AdviseUpdate(AdvIdent);
            
        } else if( !strncmp("COMP",Record+1,4) )
        {   if( Record[10]==' ' )  /* First COMPND record */
            {   ExtractString(60,Record+11,InfoMoleculeName);
                AdviseUpdate(AdvName);
            }
        } else if( !strncmp("CRYS",Record+1,4) )
        {   dst = InfoSpaceGroup;
            for( src=Record+56; *src && src<Record+67; src++ )
                if( *src!=' ' ) *dst++ = *src;
            *dst = 0;
        } else if( !strncmp("ENDM",Record+1,4) )
        {   model = False;
        } else if( !strncmp("END ",Record+1,4) ) 
            break;
    } while( !feof(DataFile) );

    if( Database )
        strcpy(InfoFileName,DataFileName);
    if( FeatList ) ProcessPDBFeatures();
    return( True );
}


int LoadXYZMolecule( fp )
    FILE *fp;
{
    return( False );
}


static int FindAlchemyResNo()
{
    register char *ptr;
    register int i, j;
    char name[4];

    if( Record[8]==' ' )
    {   name[0] = name[2] = name[3] = ' ';
        if( islower(Record[7]) )
        {   name[1] = toupper(Record[7]); 
        } else name[1] = Record[7];
        ptr = name;
    } else
    {   ptr = Record+7;
        for( i=0; i<4; i++ )
	    if( islower(ptr[i]) )
                ptr[i] = toupper(ptr[i]);

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
    register Atom __far *ptr;
    register int atoms, bonds;
    register int i;

    DataFile = fp;
    FetchRecord();
    atoms = (int)ReadValue(1,5);
    bonds = (int)ReadValue(14,5);
    ExtractString(38,Record+42,InfoMoleculeName);

    if( atoms )
    {   CreateChain( ' ' );
        if( !(CurGroup = FreeGroup) )
        {   MemSize += sizeof(Group);
            CurGroup = (Group __far *)_fmalloc(sizeof(Group));
            if( !CurGroup ) FatalDataError("Memory allocation failed");
        } else FreeGroup = CurGroup->gnext;
        strcpy(InfoFileName,DataFileName);

        CurAtom = (void __far*)0;
        CurChain->glist = CurGroup;
        CurGroup->refno = FindResNo( "MOL" );
        CurGroup->gnext = (void __far*)0;
        CurGroup->alist = (void __far*)0;
        CurGroup->serno = 1;
        CurGroup->flag = 0;
        CurGroup->col = 0;
        
        MinMainRes = MaxMainRes = 1;
        MinHetaRes = MaxHetaRes = 0;
        MainGroupCount = 1;

        for( i=0; i<atoms; i++ )
        {   FetchRecord();
            ptr = CreateAtom();

            ptr->refno = FindAlchemyResNo();
            ptr->temp = (int)ReadValue(41,8);
            /* ptr->serno = (int)ReadValue(1,5); */
            ptr->serno = i;

            ptr->xorg = ReadValue(13,7)/4;
            ptr->yorg = ReadValue(22,7)/4;
            ptr->zorg = ReadValue(31,7)/4;
            ProcessAtom( ptr );
        }

        for( i=0; i<bonds; i++ )
        {   /* InfoBondCount++; */
            FetchRecord();
        }
    }
    return( True );
}



int SavePDBMolecule( filename )
    char *filename;
{
    register Real x, y, z;
    register float xpos, ypos, zpos;
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

            x = aptr->xorg/250.0;
            y = aptr->yorg/250.0;
            z = aptr->zorg/250.0;

            /* Apply Current Viewpoint Rotation Matrix */
            xpos = (float)(x*RotX[0] + y*RotX[1] + z*RotX[2]);
            ypos = (float)(x*RotY[0] + y*RotY[1] + z*RotY[2]);
            zpos = (float)(x*RotZ[0] + y*RotZ[1] + z*RotZ[2]);

#ifdef INVERT
            fprintf(DataFile,"%8.3f%8.3f%8.3f",xpos,-ypos,-zpos);
#else
            fprintf(DataFile,"%8.3f%8.3f%8.3f",xpos,ypos,-zpos);
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
    return( True );
}


int SaveXYZMolecule( filename )
    char *filename;
{
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

            switch( GetElemIdent(aptr) )
            {   case( HandC ):   if( aptr->mbox == -1 )
                                 {   ptr = "CAR ";
                                 } else if( aptr->mbox == 1 )
                                 {   ptr = "C3  ";
                                 } else ptr = "C2  ";
                                 fputs( ptr, DataFile );
                                 break;

                case( HandN ):   if( aptr->mbox == -1 )
                                 {   ptr = "NAR ";
                                 } else ptr = "N2  ";
                                 fputs( ptr, DataFile );
                                 break;

                case( HandO ):   if( aptr->mbox == 2 )
                                 {   ptr = "O2  ";
                                 } else ptr = "O3  ";
                                 fputs( ptr, DataFile );
                                 break;

                case( HandH ):   fputs( "H   ", DataFile );  break;

                default:  ptr = ElemDesc[aptr->refno];
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
            } else if( bptr->flag & DoubBondFlag )
            {   ptr = "DOUBLE\n";
            } else ptr = "SINGLE\n";
            fputs( ptr, DataFile );
        }

    ForEachAtom
        if( aptr->flag & SelectFlag )
            aptr->mbox = 0;
    fclose( DataFile );
    return( True );
}


static void TestBonded( sptr, dptr )
    Atom __far *sptr, __far *dptr;
{
    register Bond __far *bptr;
    register Long dx, dy, dz;
    register Long max, min;
    register Long dist;
    register int flag;


    if( ((sptr->flag)&HydrogenFlag) || ((dptr->flag)&HydrogenFlag) )
    {    flag = HydrBondFlag;
         max = MaxHBondDist;
         min = MinHBondDist;
    } else
    {    flag = NormBondFlag;
         max = MaxBondDist;
         min = MinBondDist;
    }

    dx = sptr->xorg-dptr->xorg;   if( (dist=dx*dx)>max ) return;
    dy = sptr->yorg-dptr->yorg;   if( (dist+=dy*dy)>max ) return;
    dz = sptr->zorg-dptr->zorg;   if( (dist+=dz*dz)>max ) return;
    if( dist>=min )  
    {   /* Reset Non-bonded flags! */
        sptr->flag &= ~NonBondFlag;
        dptr->flag &= ~NonBondFlag;

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


void CreateMoleculeBonds( info )
    int info;
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

    ReclaimBonds( CurMolecule->blist );
    ResetVoxelData();

    ForEachAtom
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
                        do { TestBonded(aptr,dptr);
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


void TestDisulphideBridge( grp1, grp2, cys1 )
    Group __far *grp1, __far *grp2;  
    Atom __far *cys1;
{
    register HBond __far *ptr;
    register Atom __far *cys2;
    register int dx, dy, dz;
    register Long max,dist;

    if( !(cys2=FindGroupAtom(grp2,20)) )
        return;

    max = (Long)750*750;
    dx = (int)(cys1->xorg-cys2->xorg);   if( (dist=(Long)dx*dx)>max ) return;
    dy = (int)(cys1->yorg-cys2->yorg);   if( (dist+=(Long)dy*dy)>max ) return;
    dz = (int)(cys1->zorg-cys2->zorg);   if( (dist+=(Long)dz*dz)>max ) return;

    if( !(ptr = FreeHBond) )
    {   MemSize += sizeof(HBond);
        ptr = (HBond __far*)_fmalloc(sizeof(HBond));
        if( !ptr ) FatalDataError("Memory allocation failed");
    } else FreeHBond = ptr->hnext;

    ptr->dst = cys1;
    if( !(ptr->dstCA=FindGroupAtom(grp1,1)) )
        ptr->dstCA = cys1;

    ptr->src = cys2;
    if( !(ptr->srcCA=FindGroupAtom(grp2,1)) )
        ptr->srcCA = cys2;

    ptr->offset = 0;
    ptr->energy = 0;
    ptr->flag = 0;
    ptr->col = 0;

    ptr->hnext = CurMolecule->slist;
    CurMolecule->slist = ptr;

    grp1->flag |= CystineFlag;
    grp2->flag |= CystineFlag;
    InfoSSBondCount++;
}


void FindDisulphideBridges()
{
    register Chain __far *chn1;
    register Chain __far *chn2;
    register Group __far *grp1;
    register Group __far *grp2;
    register Atom __far *cys;
    char buffer[40];

    if( !Database ) return;
    ReclaimHBonds( CurMolecule->slist );
    InfoSSBondCount = 0;

    for(chn1=Database->clist;chn1;chn1=chn1->cnext)
        for(grp1=chn1->glist;grp1;grp1=grp1->gnext)
            if( IsCysteine(grp1->refno) && (cys=FindGroupAtom(grp1,20)) )
            {   for(grp2=grp1->gnext;grp2;grp2=grp2->gnext)
                    if( IsCysteine(grp2->refno) )
                        TestDisulphideBridge(grp1,grp2,cys);

                for(chn2=chn1->cnext;chn2;chn2=chn2->cnext)
                    for(grp2=chn2->glist;grp2;grp2=grp2->gnext)
                        if( IsCysteine(grp2->refno) )
                            TestDisulphideBridge(grp1,grp2,cys);
            }

    if( CommandActive )
        WriteChar('\n');
    CommandActive=False;
    
    sprintf(buffer,"Number of Bridges ... %d\n\n",InfoSSBondCount);
    WriteString(buffer);
}

#if defined(__STDC__) || defined(IBMPC)
static void CreateHydrogenBond( Atom __far*, Atom __far*,
                                Atom __far*, Atom __far*, int, int );
#endif

static void CreateHydrogenBond( srcCA, dstCA, src, dst, energy, offset )
    Atom __far *srcCA, __far *dstCA;
    Atom __far *src, __far *dst;
    int energy, offset;
{
    register HBond __far *ptr;
    register int i;

    if( !(ptr = FreeHBond) )
    {   MemSize += HBondPool*sizeof(HBond);
        ptr = (HBond __far *)_fmalloc( HBondPool*sizeof(HBond) );
        if( !ptr ) FatalDataError("Memory allocation failed");
        for( i=1; i<HBondPool; i++ )
        {   ptr->hnext = FreeHBond;
            FreeHBond = ptr++;
        } 
    } else FreeHBond = ptr->hnext;

    if( (offset>=-128) && (offset<127) )
    {   ptr->offset = (Char)offset;
    } else ptr->offset = 0;

    ptr->src = src;
    ptr->dst = dst;
    ptr->srcCA = srcCA;
    ptr->dstCA = dstCA;
    ptr->energy = energy;
    ptr->flag = 0;
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


static int CalculateBondEnergy( grp )
    Group __far *grp;
{
    register double dho,dhc;
    register double dnc,dno;

    register Atom __far *cptr;
    register Long dx,dy,dz;
    register Long dist;
    register int result;

    if( !(cptr=FindGroupAtom(grp,2)) )  return(0);
    if( !(optr=FindGroupAtom(grp,3)) )  return(0);

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
    register Group __far *grp1, __far *grp2;
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
            for(grp1=chn1->glist;grp1;grp1=grp1->gnext)
            {   pos1++;
                if( pc1 && po1 )
                {   dx = (int)(pc1->xorg - po1->xorg);
                    dy = (int)(pc1->yorg - po1->yorg);
                    dz = (int)(pc1->zorg - po1->zorg);
                } else
                {   pc1 = FindGroupAtom(grp1,2);
                    po1 = FindGroupAtom(grp1,3);
                    continue;
                }

                pc1 = FindGroupAtom(grp1,2);
                po1 = FindGroupAtom(grp1,3);

                if( !IsAmino(grp1->refno) || IsProline(grp1->refno) )
                    continue;

                if( !(ca1=FindGroupAtom(grp1,1)) ) continue;
                if( !(n1=FindGroupAtom(grp1,0)) )  continue;

                dist = (Long)dx*dx + (Long)dy*dy + (Long)dz*dz;
                dco = sqrt( (double)dist )/250.0;

                nxorg = (int)n1->xorg;   hxorg = nxorg + (int)(dx/dco);
                nyorg = (int)n1->yorg;   hyorg = nyorg + (int)(dy/dco);
                nzorg = (int)n1->zorg;   hzorg = nzorg + (int)(dz/dco);
                res1 = res2 = 0;

                for(chn2=Database->clist;chn2;chn2=chn2->cnext)
                {   /* Only consider non-empty peptide chains! */
                    if( !chn2->glist || !IsProtein(chn2->glist->refno) )
                        continue;

                    pos2 = 0;
                    for(grp2=chn2->glist;grp2;grp2=grp2->gnext)
                    {   pos2++;
                        if( (grp2==grp1) || (grp2->gnext==grp1) )
                            continue;

                        if( !IsAmino(grp2->refno) ) 
                            continue;
                        if( !(ca2=FindGroupAtom(grp2,1)) ) 
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

                        if( (energy = CalculateBondEnergy(grp2)) )
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
                    }  /* grp2 */
                }      /* chn2 */

                if( res1 ) 
                {   if( res2 ) 
                        CreateHydrogenBond(ca1,best2CA,n1,best2,res2,off2);
                    CreateHydrogenBond(ca1,best1CA,n1,best1,res1,off1);
                }
            }
        } else if( IsNucleo(chn1->glist->refno) )
            for(grp1=chn1->glist;grp1;grp1=grp1->gnext)
            {   if( !IsPurine(grp1->refno) ) continue;
                /* Find N1 of Purine Group */
                if( !(n1=FindGroupAtom(grp1,21)) )
                    continue;

                /* Maximum N1-N3 distance 5A */
                refno = NucleicCompl(grp1->refno);
                max = (Long)1250*1250;
                best = (void __far*)0;

                for(chn2=Database->clist;chn2;chn2=chn2->cnext)
                {   /* Only consider non-empty peptide chains! */
                    if( (chn1==chn2) || !chn2->glist || 
                        !IsNucleo(chn2->glist->refno) )
                        continue;

                    for(grp2=chn2->glist;grp2;grp2=grp2->gnext)
                        if( grp2->refno == (Byte)refno )
                        {   /* Find N3 of Pyramidine Group */
                            if( !(ca1=FindGroupAtom(grp2,23)) )
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
                            best = grp2;
                            max = dist;
                        }
                }

                if( best )
                {   /* Find the sugar phosphorous atoms */
                    ca1 = FindGroupAtom( grp1, 7 );
                    ca2 = FindGroupAtom( best, 7 );

                    CreateHydrogenBond( ca1, ca2, n1, best1, 0, 0 );
                    if( IsGuanine(grp1->refno) )
                    {   /* Guanine-Cytosine */
                        if( (ca1=FindGroupAtom(grp1,22)) &&    /* G.N2 */
                            (ca2=FindGroupAtom(best,26)) )     /* C.O2 */
                            CreateHydrogenBond( (void __far*)0, (void __far*)0,
                                                ca1, ca2, 0, 0 );

                        if( (ca1=FindGroupAtom(grp1,28)) &&    /* G.O6 */
                            (ca2=FindGroupAtom(best,24)) )     /* C.N4 */
                            CreateHydrogenBond( (void __far*)0, (void __far*)0,
                                                ca1, ca2, 0, 0 );

                    } else /* Adenine-Thymine */
                        if( (ca1=FindGroupAtom(grp1,25)) &&    /* A.N6 */
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

#if defined(__STDC__) || defined(IBMPC)
static int IsHBonded( Atom __far*, Atom __far*, HBond __far* );
static void TestLadder( Chain __far* );
#endif

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
                        {   if( !(first->flag&HelixFlag) ) 
                                InfoHelixCount++;

                            ptr = first;
                            do {
                                ptr->flag |= flag;
                                ptr = ptr->gnext;
                            } while( ptr != group );
                        } else prev = True;
                    } else prev = False;
                } else while( hbond && hbond->srcCA==srcCA )
                    hbond = hbond->hnext;
            } else prev = False;

            if( group->flag&HelixFlag )
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
                        {   curri->flag |= SheetFlag;
                            currj->flag |= SheetFlag;
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
                        if( curri->flag & SheetFlag )
                        {   if( ladder )
                            {   InfoLadderCount++;
                            } else ladder = True;
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
    register Long ux,uy,uz,mu;
    register Long vx,vy,vz,mv;
    register int i,found,len;
    register Real CosKappa;

    for( chain=Database->clist; chain; chain=chain->cnext )
        if( chain->glist && IsProtein(chain->glist->refno) )
        {   len = 0;  found = False;
            for( group=chain->glist; group; group=group->gnext )
            {    if( group->flag & BreakFlag )
                 {   found = False;
                     len = 0;
                 } else if( len==5 )
                 {   for( i=0; i<4; i++ )
                         aptr[i] = aptr[i+1];
                     len = 4;
                 } else if( len==2 )
                     prev = group;

                 aptr[len++] = FindGroupAtom(group,1);
                 if( len==5 ) 
                 {   if( !(prev->flag&(HelixFlag|SheetFlag)) &&
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
                                 prev->flag |= TurnFlag;
                             }
                         }
                     }
                     found = prev->flag&TurnFlag;
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
                group->flag &= ~(HelixFlag|SheetFlag|TurnFlag);

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
    sprintf(buffer,"Number of Ladders ... %d\n",InfoLadderCount);
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

    *InfoMoleculeName = 0;	AdviseUpdate(AdvName);
    *InfoClassification = 0;	AdviseUpdate(AdvClass);
    *InfoIdentCode = 0;		AdviseUpdate(AdvIdent);
    *InfoSpaceGroup = 0;
    *InfoFileName = 0;

    VoxelsClean = False;
    HMinMaxFlag = False;
    MMinMaxFlag = False;
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


void InitialiseDatabase()
{
    MinHBondDist = (Long)175*175;  MaxHBondDist = (Long)300*300;
    MinBondDist  = (Long)250*250;  MaxBondDist = (Long)475*475;
    AbsMaxBondDist = 475;

    FreeMolecule = (void __far*)0;
    FreeHBond = (void __far*)0;
    FreeChain = (void __far*)0;
    FreeGroup = (void __far*)0;
    FreeAtom = (void __far*)0;
    FreeBond = (void __far*)0;

    ResetDatabase();
}

