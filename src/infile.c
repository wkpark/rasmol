/* infile.c
 * RasMol2 Molecular Graphics
 * Roger Sayle, August 1995
 * Version 2.6
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

#define INFILE
#include "infile.h"
#include "molecule.h"
#include "abstree.h"
#include "command.h"
#include "transfor.h"

#ifndef APPLEMAC
#ifndef IBMPC
#ifndef VMS
#include <sys/types.h>
#endif
#include <sys/time.h>
#endif
#include <time.h>
#endif
 

#ifdef MMIO
#include "mmio.h"
#endif

#define GroupPool    8

#define FeatHelix    1
#define FeatSheet    2
#define FeatTurn     3


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

#define MAXALCATOM   5
static ConvTable AlcAtomTable[MAXALCATOM] = {
    { { 'S', 'O', '2', ' ' }, { ' ', 'S', '2', ' ' } },  /*  1 */
    { { 'C', 'A', 'R', ' ' }, { ' ', 'C', ' ', ' ' } },  /*  2 */
    { { 'N', 'A', 'R', ' ' }, { ' ', 'N', ' ', ' ' } },  /*  3 */
    { { 'N', 'A', 'M', ' ' }, { ' ', 'N', ' ', ' ' } },  /*  4 */
    { { 'N', 'P', 'L', '3' }, { ' ', 'N', '3', ' ' } },  /*  5 */
				 };

static char PDBInsert;
static Feature __far *FeatList;
static Atom __far *ConnectAtom;
static char Record[202];
static FILE *DataFile;

/* Macros for commonly used loops */
#define ForEachAtom  for(chain=Database->clist;chain;chain=chain->cnext) \
		     for(group=chain->glist;group;group=group->gnext)    \
		     for(aptr=group->alist;aptr;aptr=aptr->anext)
#define ForEachBond  for(bptr=Database->blist;bptr;bptr=bptr->bnext)


/* Forward Reference */
void DestroyDatabase();


#ifdef APPLEMAC
/* External RasMac Function Declaration! */
void SetFileInfo( char*, OSType, OSType, short );
#endif


static void FatalInFileError(ptr)
    char *ptr;
{
    char buffer[80];

    sprintf(buffer,"InFile Error: %s!",ptr);
    RasMolFatalExit(buffer);
}


/*================================*/
/* File/String Handling Functions */
/*================================*/

static int FetchRecord()
{
    register char *ptr;
    register int ch;

    if( feof(DataFile) )
    {   *Record = '\0';
        return( False );
    }

    ptr = Record;
    do {
        ch = getc(DataFile);
        if( ch == '\n' )
        {   *ptr = 0;
            return( True );
        } else if( ch == '\r' )
        {   ch = getc(DataFile);
            if( ch != '\n' )
                ungetc(ch,DataFile);
            *ptr = 0;
            return( True );
        } else if( ch == EOF )
        {   *ptr = 0;
            return( ptr != Record+1 );
        } else *ptr++ = ch;
    } while( ptr < Record+200 );

    /* skip to the end of the line! */
    do { ch = getc(DataFile);
    } while( (ch!='\n') && (ch!='\r') && (ch!=EOF) );

    if( ch == '\r' )
    {   ch = getc(DataFile);
        if( ch != '\n' )
            ungetc(ch,DataFile);
    }
    *ptr = 0;
    return( True );
}


static void ExtractString( len, src, dst )
    int len;  char *src, *dst;
{
    register char *ptr;
    register char ch;
    register int i;

    ptr = dst;
    for( i=0; i<len; i++ )
    {   if( *src )
	{   ch = *src++;
            *dst++ = ch;
            if( ch != ' ' ) 
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


/*===================================*/
/* File Format Independent Functions */
/*===================================*/

#ifdef FUNCPROTO
static void UpdateFeature( FeatEntry __far*, int );
#endif

static FeatEntry __far *AllocFeature()
{
    register Feature __far *ptr;
 
    if( !FeatList || (FeatList->count==FeatSize) )
    {   ptr = (Feature __far*)_fmalloc(sizeof(Feature));
        if( !ptr ) FatalInFileError("Memory allocation failed");
        /* Features are always deallocated! */
 
        ptr->fnext = FeatList;
        ptr->count = 0;
        FeatList = ptr;
    } else ptr = FeatList;
 
    return( &(ptr->data[ptr->count++]) );
}


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
 
            if( NMRModel )
            {  continue;
            } else return;
        }
}
 
 
static void ProcessFeatures()
{
    register Feature __far *next;
    register Feature __far *ptr;
    register int i;
 
    Info.turncount = 0;
    Info.helixcount = 0;
    Info.laddercount = 0;
    Info.structsource = SourcePDB;
 
    for( ptr=FeatList; ptr; ptr=next )
    {    if( Database )
             for( i=0; i<ptr->count; i++ )
                 if( ptr->data[i].type==FeatHelix )
                 {   UpdateFeature( &ptr->data[i], HelixFlag );
                     Info.helixcount++;
                 } else if( ptr->data[i].type==FeatSheet )
                 {   UpdateFeature( &ptr->data[i], SheetFlag );
                     Info.laddercount++;
                 } else /* FeatTurn */
                 {   UpdateFeature( &ptr->data[i], TurnFlag );
                     Info.turncount++;
                 }
 
         /* Deallocate Memory */
         next = ptr->fnext;
         _ffree( ptr );
    }
}
 
 


/*==============================*/
/* Molecule File Format Parsing */
/*==============================*/

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


static void ProcessPDBGroup( heta, serno )
    int heta, serno;
{
    PDBInsert = Record[26];
    if( !CurChain || (CurChain->ident!=Record[21]) )
    {   ConnectAtom = (Atom __far*)0;
        CreateChain( Record[21] );
    }
    CreateGroup( GroupPool );
 
    CurGroup->refno = FindResNo( Record+17 );
    CurGroup->serno = serno;
    ProcessGroup( heta );
}
 

static void ProcessPDBAtom( heta )
    int heta;
{
    register Bond __far *bptr;
    register Atom __far *ptr;
    register Long dx,dy,dz;
    register int temp,serno;
 
    dx = ReadPDBCoord(30);
    dy = ReadPDBCoord(38);
    dz = ReadPDBCoord(46);
 
    /* Process Pseudo Atoms Limits!! */
    if( (Record[13]=='Q') && (Record[12]==' ') )
    {   temp = (int)ReadValue(60,6);
        if( MMinMaxFlag )
        {   if( temp < MinMainTemp )
            {   MinMainTemp = temp;
            } else if( temp > MaxMainTemp )
                MaxMainTemp = temp;
        }
 
        /* Dummy co-ordinates! */
        if( (dx==dy) && (dx==dz) )
        {   if( !dx || (dx == 9999000L) )
                return;
        }
 
        if( HMinMaxFlag || MMinMaxFlag )
        {   if( dx < MinX )
            {   MinX = dx;
            } else if( dx > MaxX )
                MaxX = dx;
 
            if( dy < MinY )
            {   MinY = dy;
            } else if( dy > MaxY )
                MaxY = dy;
 
            if( dz < MinZ )
            {   MinZ = dz;
            } else if( dz > MaxZ )
                MaxZ = dz;
        }
        return;
    }
 
 
    /* Ignore XPLOR Pseudo Atoms!! */
    if( (dx==9999000L) && (dy==9999000L) && (dz==9999000L) )
        return;
 
    serno = (int)ReadValue(22,4);
    if( !CurGroup || (CurGroup->serno!=serno)
        || (CurChain->ident!=Record[21])
        || (PDBInsert!=Record[26]) )
        ProcessPDBGroup( heta, serno );
 
 
    ptr = CreateAtom();
    ptr->refno = ComplexAtomType(Record+12);
    ptr->serno = ReadValue(6,5);
    ptr->temp = (int)ReadValue(60,6);
    ptr->altl = Record[16];
 
    /* Handle CONCORD PDB Files */
    if( (Record[12]==' ') && islower(Record[14])
        && !strncmp(Record+15,"       ",7) ) 
        ptr->refno = SimpleAtomType(Record+13);

    ptr->xorg =  dx/4;
    ptr->yorg =  dy/4;
    ptr->zorg = -dz/4;
 
    if( heta ) ptr->flag |= HeteroFlag;
    ProcessAtom( ptr );
 
    /* Create biopolymer Backbone */
    if( IsAlphaCarbon(ptr->refno) && IsProtein(CurGroup->refno) )
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
    } else if( IsSugarPhosphate(ptr->refno) && IsNucleo(CurGroup->refno) )
    {   if( ConnectAtom )
        {   bptr = ProcessBond(ConnectAtom,ptr,NormBondFlag);
            bptr->bnext = CurChain->blist;
            CurChain->blist = bptr;
        }
        ConnectAtom = ptr;
    }
}
 

static void ProcessPDBColourMask()
{
    register MaskDesc *ptr;
    register char *mask;
    register int i;
 
    if( MaskCount==MAXMASK )
        FatalInFileError("Too many COLOR records in file");
    ptr = &UserMask[MaskCount];
    mask = ptr->mask;
 
 
    ptr->flags = 0;
    for( i=6; i<11; i++ )
        if( (*mask++ = Record[i]) != '#' )
            ptr->flags |= SerNoFlag;
 
    for( i=12; i<20; i++ )
        *mask++ = Record[i];
    *mask++ = Record[21];
 
    for( i=22; i<26; i++ )
        if( (*mask++ = Record[i]) != '#' )
            ptr->flags |= ResNoFlag;
    *mask++ = Record[26];
 
    ptr->r = (int)(ReadPDBCoord(30)>>2) + 5;
    ptr->g = (int)(ReadPDBCoord(38)>>2) + 5;
    ptr->b = (int)(ReadPDBCoord(46)>>2) + 5;
    ptr->radius = (short)(5*ReadValue(54,6))>>1;
    MaskCount++;
}
 

int LoadPDBMolecule( fp, flag )
    FILE *fp;  int flag;
{
    register FeatEntry __far *ptr;
    register int srcatm, dstatm;
    register char *src, *dst;
    register int i,ignore;
 
    ignore = False;
    FeatList = (void __far*)0;
    DataFile = fp;
    NMRModel = 0;
 
    while( FetchRecord() )
    {   if( *Record == 'A' )
        {   if( !ignore && !strncmp("ATOM",Record,4) )
                ProcessPDBAtom( False );

        } else switch(*Record)
        {   case('C'):    if( !strncmp("CONE",Record,4) )
                          {   if( ignore || flag ) continue;
 
                              srcatm = (int)ReadValue(6,5);
                              if( srcatm )
                                  for( i=11; i<=36 && Record[i]; i+=5 )
                                  {   dstatm = (int)ReadValue(i,5);
                                      if( dstatm && (dstatm>srcatm) )
                                          CreateBondOrder(srcatm,dstatm);
                                  }
                               
                          } else if( !strncmp("COMP",Record,4) )
                          {   /* First or MOLECULE: COMPND record */
                              if( (Record[9]==' ') && 
                                  strncmp(Record+10,"MOL_ID:",7) )  
                              {   ExtractString(60,Record+10,
                                                Info.moleculename);
                              } else if( !Info.moleculename[0] &&
                                         !strncmp(Record+11,"MOLECULE: ",10) )
                                  ExtractString(49,Record+21,
                                                Info.moleculename);
                          } else if( !strncmp("CRYS",Record,4) )
                          {   dst = Info.spacegroup;  src=Record+55;
                              while( *src && src<Record+66 )
                                  if( *src!=' ' ) 
                                  {   *dst++ = *src++;
                                  } else src++;
                              *dst = 0;
 
                              Info.cella = ReadValue( 6,9)/1000.0;
                              Info.cellb = ReadValue(15,9)/1000.0;
                              Info.cellc = ReadValue(24,9)/1000.0;
 
                              Info.cellalpha = Deg2Rad*(ReadValue(33,7)/100.0);
                              Info.cellbeta =  Deg2Rad*(ReadValue(40,7)/100.0);
                              Info.cellgamma = Deg2Rad*(ReadValue(47,7)/100.0);
 
                          } else if( !strncmp("COLO",Record,4) )
                              ProcessPDBColourMask();
                          break;

            case('E'):    if( !strncmp("ENDM",Record,4) )
                          {   /* break after single model??? */
                              if( flag )
                              {   ConnectAtom = (void __far*)0;
                                  CurGroup = (void __far*)0;
                                  CurChain = (void __far*)0;
                              } else ignore = True;
 
                          } else if( !strncmp("END",Record,3) )
                              if( !Record[4] || (Record[4]==' ') )
                              {   /* Treat END same as TER! */
                                  ConnectAtom = (void __far*)0;
                                  CurGroup = (void __far*)0;
                                  CurChain = (void __far*)0;
                              }
                          break;

            case('H'):    if( !strncmp("HETA",Record,4) )
                          {   if( !ignore ) ProcessPDBAtom(True);
                          } else if( !strncmp("HELI",Record,4) )
                          {   if( ignore ) continue;
 
                              /* Remaining HELIX record fields   */
                              /* 38-39 .... Helix Classification */
                              /* 31 ....... Same Chain as 19?    */
                              ptr = AllocFeature();
                              ptr->type = FeatHelix;
                              ptr->chain = Record[19];
                              ptr->init = (int)ReadValue(21,4);
                              ptr->term = (int)ReadValue(33,4);
                              
                          } else if( !strncmp("HEAD",Record,4) )
                          {   ExtractString(40,Record+10,Info.classification);
                              ExtractString( 4,Record+62,Info.identcode);
                          }
                          break;

            case('M'):    if( !strncmp("MODE",Record,4) )
                              if( flag ) NMRModel++;
                          break;
 
            case('S'):    if( !strncmp("SHEE",Record,4) )
                          {   if( ignore ) break;
                              /* Remaining SHEET record fields   */
                              /* 38-39 .... Strand Parallelism   */
                              /* 32 ....... Same Chain as 21?    */
                              ptr = AllocFeature();
                              ptr->type = FeatSheet;
                              ptr->chain = Record[21];
                              ptr->init = (int)ReadValue(22,4);
                              ptr->term = (int)ReadValue(33,4);
                          }
                          break;

            case('T'):    if( !strncmp("TURN",Record,4) )
                          {   if( ignore ) continue;
 
                              ptr = AllocFeature();
                              ptr->type = FeatTurn;
                              ptr->chain = Record[19];
                              ptr->init = (int)ReadValue(20,4);
                              ptr->term = (int)ReadValue(31,4);
                          } else if( !strncmp("TER",Record,3) )
                          {   if( !Record[3] || (Record[3]==' ') )
                              {   ConnectAtom = (void __far*)0;
                                  CurGroup = (void __far*)0;
                                  CurChain = (void __far*)0;
                              }
                          }
                          break;
        }
    }
 
    if( Database )
        strcpy(Info.filename,DataFileName);
    if( FeatList ) ProcessFeatures();
    return( True );
}


int LoadMDLMolecule( fp )
    FILE *fp;
{
    register Bond __far *bptr;
    register Atom __far *src;
    register Atom __far *dst;
    register Atom __far *ptr;
 
    register int i,type;
    register int atoms, bonds;
    register int srcatm,dstatm;
    register Long dx, dy, dz;
    register Card dist2;
    register Real scale;
    register char *cptr;
 
    DataFile = fp;
 
    FetchRecord(); /* Molecule Name */
    ExtractString(78,Record,Info.moleculename);
 
    FetchRecord(); /* Program Details */
    FetchRecord(); /* Comments */
 
    FetchRecord();
    atoms = (int)ReadValue(0,3);
    bonds = (int)ReadValue(3,3);
 
    if( !atoms )
        return( False );
 
    CreateMolGroup();
    for( i=1; i<=atoms; i++ )
    {   FetchRecord();
        ptr = CreateAtom();
 
        cptr = Record+31;
        while( *cptr == ' ' ) cptr++;
        ptr->refno = SimpleAtomType(cptr);
 
        switch( (int)ReadValue(36,3) )
        {   case(1):  ptr->temp =  300;  break;
            case(2):  ptr->temp =  200;  break;
            case(3):  ptr->temp =  100;  break;
            case(5):  ptr->temp = -100;  break;
            case(6):  ptr->temp = -200;  break;
            case(7):  ptr->temp = -300;  break;
            default:  ptr->temp = 0;
        }
        ptr->serno = i;
 
        ptr->xorg =  ReadValue( 0,10)/40;
        ptr->yorg =  ReadValue(10,10)/40;
        ptr->zorg = -ReadValue(20,10)/40;
        ProcessAtom( ptr );
    }
 
    for( i=0; i<bonds; i++ )
    {   FetchRecord();
        srcatm = (int)ReadValue(0,3);
        dstatm = (int)ReadValue(3,3);
        type =   (int)ReadValue(6,3);
 
        if( type==2 )                 /* DOUBLE */
        {   CreateBond(srcatm,dstatm,DoubBondFlag);
        } else if( type==3 )          /* TRIPLE */
        {   CreateBond(srcatm,dstatm,TripBondFlag);
        } else if( type==4 )          /* AROMATIC */
        {   CreateBond(srcatm,dstatm,AromBondFlag);
        } else                        /* SINGLE */
            CreateBond(srcatm,dstatm,NormBondFlag);
    }
 
    for( bptr=Database->blist; bptr; bptr=bptr->bnext )
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
 

int LoadXYZMolecule( fp )
    FILE *fp;
{
    auto char type[12];
    auto double xpos, ypos, zpos;
    auto double charge, u, v, w;
    auto long atoms;
 
    register Atom __far *ptr;
    register char *src,*dst;
    register int count;
    register Long i;
 
 
    DataFile = fp;
    /* Number of Atoms */
    FetchRecord();
    sscanf(Record,"%ld",&atoms);
 
    /* Molecule (step) Description */
    FetchRecord();
    src = Record;
    while( *src == ' ' )
        src++;
 
    dst = Info.moleculename;
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
            count = sscanf(Record,"%s %lf %lf %lf %lf %lf %lf %lf",
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
 

static int FindSybylRefNo( ptr )
    char *ptr;
{
    register char *src,*dst;
    auto char name[4];
 
    src = ptr;
    dst = name;
    if( ptr[1] && (ptr[1]!='.') )
    {   *dst++ = ToUpper(*src);  src++;
        *dst++ = ToUpper(*src);  src++;
    } else
    {   *dst++ = ' ';
        *dst++ = ToUpper(*src);
        src++;
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
    return( NewAtomType(name) );
}
 
 
int LoadMol2Molecule( fp )
    FILE *fp;
{
    double xpos, ypos, zpos;
    long features, sets, serno;
    long atoms, bonds, structs;
    long srcatm, dstatm;
 
    char name[8];
    char type[8];
 
    register Atom __far *ptr;
    register char *src, *dst;
    register Long i;
 
 
    DataFile = fp;
    while( FetchRecord() )
    {   if( !*Record || *Record=='#' )
            continue;
 
        if( !strncmp("@<TRIPOS>MOLECULE",Record,17) ||
            !strncmp("@MOLECULE",Record,9) )
        {   FetchRecord();  /* Molecule Name */
            src = Record;
            while( *src==' ' )
                src++;
 
            dst = Info.moleculename;
            while( (*dst++ = *src++) );
 
            FetchRecord();
            atoms = bonds = structs = features = sets = 0;
            sscanf(Record,"%ld %ld %ld %ld %ld", &atoms, &bonds, 
                          &structs, &features, &sets );
 
            FetchRecord();  /* Molecule Type  */
            FetchRecord();  /* Charge Type    */
 
        } else if( !strncmp("@<TRIPOS>ATOM",Record,13) ||
                   !strncmp("@ATOM",Record,5) )
        {   if( !atoms ) continue;
 
            CreateMolGroup();
            for( i=0; i<atoms; i++ )
            {    FetchRecord();
                 ptr = CreateAtom();
 
                 sscanf(Record,"%ld %s %lf %lf %lf %s", &serno, name,
                                &xpos, &ypos, &zpos, type );
 
                 ptr->refno = FindSybylRefNo( type );
                 ptr->serno = serno;
                 /* ptr->serno = i; */
 
                 ptr->xorg =  (Long)(250.0*xpos);
                 ptr->yorg =  (Long)(250.0*ypos);
                 ptr->zorg = -(Long)(250.0*zpos);
                 ProcessAtom( ptr );
            }
 
        } else if( !strncmp("@<TRIPOS>BOND",Record,13) ||
                   !strncmp("@BOND",Record,5) )
            for( i=0; i<bonds; i++ )
            {   FetchRecord();
                sscanf(Record,"%ld %ld %ld %s",
                              &serno,&srcatm,&dstatm,type);
                if( !strncmp(type,"ar",2) )
                {   CreateBond(srcatm,dstatm,AromBondFlag);
                } else if( *type == '2' )
                {   CreateBond(srcatm,dstatm,DoubBondFlag);
                } else /* *type == '1' */
                    CreateBond(srcatm,dstatm,NormBondFlag);
            }
    }
    return( True );
}


static int FindAlchemyRefNo()
{
    register char *ptr;
    register int i;
    char name[4];
 
    ptr = Record+6;
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
    return( NewAtomType(ptr) );
}
 

int LoadAlchemyMolecule( fp )
    FILE *fp;
{
    auto long serno,srcatm,dstatm;
    register Atom __far *ptr;
    register Long atoms, bonds;
    register char *chptr;
    register Long i;
 
    DataFile = fp;
    FetchRecord();
    atoms = ReadValue(0,5);
    bonds = ReadValue(13,5);
    ExtractString(38,Record+41,Info.moleculename);
 
    if( !atoms )
        return( False );
 
    CreateMolGroup();
    for( i=0; i<atoms; i++ )
    {   FetchRecord();
        ptr = CreateAtom();
 
        ptr->refno = FindAlchemyRefNo();
        ptr->temp = (int)ReadValue(40,8);
        ptr->serno = (int)ReadValue(0,5);
        /* ptr->serno = i+1; */
 
        ptr->xorg =  ReadValue(12,7)/4;
        ptr->yorg =  ReadValue(21,7)/4;
        ptr->zorg = -ReadValue(30,7)/4;
        ProcessAtom( ptr );
    }
 
    for( i=0; i<bonds; i++ )
    {   FetchRecord();
        sscanf(Record,"%ld %ld %ld",&serno,&srcatm,&dstatm);
 
        chptr = Record;
        while( *chptr && !isalpha(*chptr) )
            chptr++;
 
        if( *chptr =='A' )             /* AROMATIC */
        {   CreateBond(srcatm,dstatm,AromBondFlag);
        } else if( *chptr == 'D' )     /* DOUBLE */
        {   CreateBond(srcatm,dstatm,DoubBondFlag);
        } else if( *chptr == 'T' )     /* TRIPLE */
        {   CreateBond(srcatm,dstatm,TripBondFlag);
        } else if( *chptr == 'H' )     /* HYDROGEN */
        {   CreateBond(srcatm,dstatm,HydrBondFlag);
        } else /* (*chptr == 'S') */   /* SINGLE */
            CreateBond(srcatm,dstatm,NormBondFlag);
    }
    return( True );
}
 
 
int LoadCharmmMolecule( fp )
    FILE *fp;
{
    auto char buffer[4];
    register Atom __far *ptr;
    register Long atoms,serno;
    register int chain,resno;
    register int i;

    DataFile = fp;
 
    do {
        FetchRecord();
    } while( *Record=='*' );
    atoms = ReadValue(0,5);
    if( !atoms ) return False;
 
    MinHetaRes = MaxHetaRes = 0;
    strcpy(Info.filename,DataFileName);
    MainGroupCount = 0;
 
    chain = 0;
    CurChain = NULL;
    for( serno=0; serno<atoms; serno++ )
    {   FetchRecord();
 
        if( !CurChain || strncmp(Record+51,buffer,4) )
        {   for( i=0; i<4; i++ )
                buffer[i] = Record[51+i];
            ConnectAtom = (Atom __far*)0;
            CreateChain(chain+49);
            chain++;
        }

        /* Non-standard Sequential Residue Numbering */
        /* resno = (int)ReadValue(5,5); */

        resno = (int)ReadValue(56,5);
        if( !CurGroup || (CurGroup->serno!=resno) )
        {   CreateGroup( GroupPool );
            CurGroup->refno = FindResNo(Record+11);
            CurGroup->serno = resno;
            ProcessGroup( False );
        }
 
        ptr = CreateAtom();
        ptr->refno = ComplexAtomType(Record+15);
        ptr->temp = (int)ReadValue(60,9);
        ptr->serno = ReadValue(0,5);
        /* ptr->serno = serno+1; */
 
        ptr->xorg =  ReadValue(20,8)/4;
        ptr->yorg =  ReadValue(30,8)/4;
        ptr->zorg = -ReadValue(40,8)/4;
        ProcessAtom( ptr );
    }
    return( True );
}
 

static int MOPACAtomType( type )
    char *type;
{
    auto char name[4];
    register char ch1,ch2;
    register int i;
 
    if( *type == ' ' )
        type++;
 
    name[2] = name[3] = ' ';
    if( isdigit(type[0]) )
    {   i = *type++ - '0';
        while( isdigit(*type) )
            i = (10*i) + (*type++ - '0');
 
        /* Last Atom in File! */
        if( i == 0 )
        {   return( -1 );
        } else if( i >= 99 )
            return( 1 );
 
        /* Construct Name */
        ch1 = Element[i].symbol[0];
        ch2 = Element[i].symbol[1];
        if( ch2 == ' ' )
        {   name[1] = ch1;
            name[0] = ' ';
        } else
        {   name[1] = ToUpper(ch2);
            name[0] = ch1;
        }
 
    } else
    {   ch1 = ToUpper(type[0]);
        ch2 = ToUpper(type[1]);
        if( (ch1=='X') || (ch1=='+') || (ch1=='-') )
        {   return( 1 );
        } else if( (ch1=='T') && (ch2=='V') )
            return( 1 );
 
        if( ch2 && (ch2!=' ') && (ch2!='(') && !isdigit(ch2) )
        {   name[0] = ch1;
            name[1] = ch2;
        } else
        {   name[1] = ch1;
            name[0] = ' ';
        }
    }
    return( NewAtomType(name) );
}
 
 
static int ReadMOPACOutputFile()
{
    register Atom __far *atm;
    register int i,init;
    register char *ptr;
 
    init = False;
    while( FetchRecord() )
    {   ptr = Record;
        while( *ptr == ' ' )
            ptr++;
 
        if( !strncmp(ptr,"CARTESIAN COORDINATES",21) )
        {   for( i=0; i<3; i++ )
                FetchRecord();
 
            if( Database )
            {   atm = CurGroup->alist;
                MMinMaxFlag = False;
                HasHydrogen = False;
                MainAtomCount = 0;
            }
 
            while( FetchRecord() && *Record && isdigit(Record[5]) )
            {   if( !Database )
                {   atm = (Atom __far*)0;
                    CreateMolGroup();
                    init = True;
                }
 
                if( !atm )
                {   atm = CreateAtom();
                    atm->serno = (int)ReadValue(0,6);
                    atm->refno = MOPACAtomType(Record+14);
                    atm->temp = 0;
 
                    atm->xorg = ReadValue(20,10)/40;
                    atm->yorg = ReadValue(30,10)/40;
                    atm->zorg = ReadValue(40,10)/40;
                } else
                {   atm->xorg = ReadValue(30,10)/40;
                    atm->yorg = ReadValue(40,10)/40;
                    atm->zorg = ReadValue(50,10)/40;
                }
                ProcessAtom(atm);
                atm = atm->anext;
            }
             
        } else if( !strncmp(ptr,"NET ATOMIC CHARGES",18) )
        {   FetchRecord();
            FetchRecord();
 
            if( Database )
            {   atm = CurGroup->alist;
                MMinMaxFlag = False;
                HasHydrogen = False;
                MainAtomCount = 0;
            }
 
            while( FetchRecord() && strncmp(Record," DIPOLE",7) )
            {   if( !Database )
                {   atm = (Atom __far*)0;
                    CreateMolGroup();
                }
 
                if( !atm )
                {   atm = CreateAtom();
                    atm->serno = (int)ReadValue(0,12);
                    atm->refno = MOPACAtomType(Record+21);
                    atm->temp = (int)(ReadValue(27,13)/100);
                    atm->xorg = atm->yorg = atm->zorg = 0;
                } else
                    atm->temp = (int)(ReadValue(27,13)/100);
                ProcessAtom(atm);
                atm = atm->anext;
            }
        }
    }
 
    if( !init )
    {   if( Database )
            DestroyDatabase();
        return( False );
    } else return( True );
}
 
 
static int MoreMOPACKeywords()
{
    register char *ptr;
 
    ptr = Record;
    while( *ptr )
    {   if( *ptr == '+' )
            if( !ptr[1] || (ptr[1]==' ') )
                return(True);
 
        /* Skip Next Keyword */
        while( *ptr && *ptr!=' ' ) ptr++;
        while( *ptr == ' ' ) ptr++;
    }
    return( False );
}
 
 
int LoadMOPACMolecule( fp )
    FILE *fp;
{
    static int na,nb,nc,lopt;
    static double dist,angle,dihed;
    register IntCoord __far *coord;
    register Atom __far *aptr;
    register int count,refno;
    register int cartflag;
    register char *ptr;
 
    DataFile = fp;
    FetchRecord();
 
    /* Test for MOPAC output file */
    if( !strncmp(Record," ***",4) )
        return( ReadMOPACOutputFile() );
 
    /* MOPAC Keywords */
    while( MoreMOPACKeywords() )
        FetchRecord();
 
    FetchRecord(); /* Molecule Name */
    ExtractString(78,Record,Info.moleculename);
    FetchRecord(); /* Comments */
 
    count = 0;
    cartflag = False;
    InitInternalCoords();
    while( FetchRecord() )
    {   /* Process Record! */
        for( ptr=Record; *ptr; ptr++ )
            if( (*ptr==',') || (*ptr==0x09) )
                *ptr = ' ';
 
        ptr = Record;
        while( *ptr == ' ' )
           ptr++;
 
        if( !*ptr ) break;
        refno = MOPACAtomType(ptr);
        if( refno == -1 ) break;
 
        while( *ptr && (*ptr!=' ') )
            if( *ptr == '(' )
            {   /* Skip Atom Label */
                while( *ptr && (*ptr!=')') )
                    ptr++;
            } else ptr++;
 
        na = nb = nc = 0;
        dist = angle = dihed = 0.0;
        sscanf(ptr,"%lf %*d %lf %*d %lf %d %d %d %d",
               &dist, &angle, &dihed, &lopt, &na, &nb, &nc );
        count++;
 
        if( count == 3 )
        {   /* Handle missing dihedral */
            if( lopt == 2 )
            {   na = 1;  nb = 2;
                dihed = 0.0;
            } else if( lopt == 1)
            {   /* Safe FP comparison for Cartesian */
                if( (dihed>=1.9999) && (dihed<=2.0001) )
                {   na = 2;  nb = 1;
                    dihed = 0.0;
                }
            }
        } else if( count == 4 )
            cartflag = (na == 0);
 
        coord = AllocInternalCoord();
        coord->na = na; coord->nb = nb; coord->nc = nc;
        coord->refno = refno;
        coord->angle = angle;
        coord->dihed = dihed;
        coord->dist = dist;
    }
 
    if( !count )
        return( False );
 
    /* Co-ordinate conversion! */
    if( !cartflag )
        if( !ConvertInternal2Cartesian() )
        {   FreeInternalCoords();
            if( CommandActive )
                WriteChar('\n');
            WriteString("Error: Invalid MOPAC z-matrix file!\n\n");
            CommandActive=False;
            return( False );
        }
 
    count = 0;
    for( coord=IntList; coord; coord=coord->inext )
        if( coord->refno != 1 )
        {   if( !Database )
                CreateMolGroup();
 
            aptr = CreateAtom();
            aptr->refno = (Byte)coord->refno;
            aptr->serno = ++count;
            aptr->temp = 0;
 
            aptr->xorg =  (Long)(250.0*coord->dist);
            aptr->yorg =  (Long)(250.0*coord->angle);
            aptr->zorg = -(Long)(250.0*coord->dihed);
            ProcessAtom(aptr);
        } else count++;
 
    FreeInternalCoords();
    return( True );
}
 

int LoadMacroModelMolecule( fp )
    FILE *fp;
{
#ifdef MMIO
    return( True );
#else
    register char *src,*dst;
    register int i;
    auto int atoms;
 
    DataFile = fp;
 
    /* Number of Atoms & Description */
    FetchRecord();
    sscanf(Record,"%d",&atoms);
 
    src = Record;
    dst = Info.moleculename;
    for( i=0; i<78; i++ )
        if( *src ) *dst++ = *src++;
    *dst = '\0';
 
    for( i=0; i<atoms; i++ )
    {    FetchRecord();
    }
    return( True );
#endif
}


int LoadBiosymMolecule( fp )
    FILE *fp;
{
    return( True );
}
 
 
int LoadSHELXMolecule( fp )
    FILE *fp;
{
    return( True );
}
 

int LoadFDATMolecule( fp )
    FILE *fp;
{
    return( True );
}


/*=================================*/
/* Molecule File Format Generation */
/*=================================*/

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
 
    if( *Info.classification || *Info.identcode )
    {   fputs("HEADER    ",DataFile);
 
        ptr = Info.classification;
        for( i=11; i<=50; i++ )
            putc( (*ptr ? *ptr++ : ' '), DataFile );
        fprintf(DataFile,"13-JUL-93   %.4s\n",Info.identcode);
    }
 
    if( *Info.moleculename )
        fprintf(DataFile,"COMPND    %.60s\n",Info.moleculename);
 
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
 
 
int SaveMDLMolecule( filename )
    char *filename;
{
    register Chain __far *chain;
    register Group __far *group;
    register Atom __far *aptr;
    register Bond __far *bptr;
    register int atoms,bonds;
    register double x,y,z;
    register int ch,temp;
    register int atomno;

#ifndef APPLEMAC
    register struct tm *ptr;
    static time_t curtime;
#endif
 
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

    /* Write Mol file header */
    fprintf(DataFile,"%.80s\n  RasMol  ",Info.moleculename);

#ifndef APPLEMAC
    curtime = time((time_t*)0);
    ptr = localtime(&curtime);
    fprintf(DataFile,"%02d%02d%02d%02d%02d",ptr->tm_mon+1,ptr->tm_mday,
                             ptr->tm_year%100,ptr->tm_hour,ptr->tm_min);
#else
    fputs("0101951200",DataFile);
#endif

    atoms = 0;
    ForEachAtom
        if( aptr->flag & SelectFlag )
            atoms++;

    bonds = 0;
    ForEachBond
        if( bptr->srcatom->flag & bptr->dstatom->flag & SelectFlag )
             bonds++;

    fprintf(DataFile,"%cD\n\n", (MinZ||MaxZ)?'3':'2' );
    fprintf(DataFile,"%3d%3d",atoms,bonds);
    fputs("  0        0              1 V2000\n",DataFile);

    atomno = 1;
    ForEachAtom
        if( aptr->flag & SelectFlag )
        {   x = (double)aptr->xorg/250.0;
            y = (double)aptr->yorg/250.0;
            z = (double)aptr->zorg/250.0;
#ifdef INVERT
#ifdef __STDC__
            fprintf(DataFile,"%10.4f%10.4f%10.4f ",x,-y,-z);
#else
            fprintf(DataFile,"%10.4lf%10.4lf%10.4lf ",x,-y,-z);
#endif
#else
#ifdef __STDC__
            fprintf(DataFile,"%10.4f%10.4f%10.4f ",x,y,-z);
#else
            fprintf(DataFile,"%10.4lf%10.4lf%10.4lf ",x,y,-z);
#endif
#endif

            fputc(Element[aptr->elemno].symbol[0],DataFile);
            fputc(Element[aptr->elemno].symbol[1],DataFile);
            fputs("  0  ",DataFile);

            ch = '0';
            /* Test for charges or b-factors */
            if( (MinMainTemp<0) || (MinHetaTemp<0) )
            {   temp = aptr->temp;
                if( temp > 50 )
                {   if( temp > 250 )
                    {   ch = '1';
                    } else if( temp > 150 )
                    {   ch = '2';
                    } else ch = '3';
                } else if( temp < -50 )
                {   if( temp < -250 )
                    {   ch = '7';
                    } else if( temp < -150 )
                    {   ch = '6';
                    } else ch = '5';
                } else ch = '0';
            } 
            fputc(ch,DataFile);

            fputs("  0  0  0  0  0  0  0  0  0  0\n",DataFile);
            aptr->mbox = atomno++;
        }

    ForEachBond
        if( bptr->srcatom->flag & bptr->dstatom->flag & SelectFlag )
        {   fprintf(DataFile,"%3d%3d  ",bptr->srcatom->mbox,
                                        bptr->dstatom->mbox);
            if( bptr->flag & AromBondFlag )
            {      fputc('4',DataFile);
            } else if( bptr->flag & TripBondFlag )
            {      fputc('3',DataFile);
            } else if( bptr->flag & DoubBondFlag )
            {      fputc('2',DataFile);
            } else fputc('1',DataFile);
            fputs("  0\n",DataFile);
        }

    fputs("M  END\n",DataFile);
    fclose( DataFile );
#ifdef APPLEMAC
    SetFileInfo(filename,'RSML','mMOL',131);
#endif
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
            aptr->mbox = 0;
 
    ForEachBond
        if( ((bptr->srcatom->flag&bptr->dstatom->flag)&SelectFlag) &&
           !((bptr->srcatom->flag|bptr->dstatom->flag)&HydrogenFlag) )
        {   if( bptr->flag&AromBondFlag )
            {   bptr->srcatom->mbox = -1;
                bptr->dstatom->mbox = -1;
            } else
            {   num = (bptr->flag&DoubBondFlag)? 2 : 1;
                if( bptr->srcatom->mbox>0 )
                    bptr->srcatom->mbox += num;
                if( bptr->dstatom->mbox>0 )
                    bptr->dstatom->mbox += num;
            }
        }
 
    fprintf(DataFile,"%5ld ATOMS, ",(long)(MainAtomCount+HetaAtomCount));
    fprintf(DataFile,"%5ld BONDS, ",(long)Info.bondcount);
    fprintf(DataFile,"    0 CHARGES, %s\n", Info.moleculename );
 
    atomno = 1;
    ForEachAtom
        if( aptr->flag & SelectFlag )
        {   aptr->mbox = atomno;
            fprintf(DataFile,"%5d ",atomno++);
 
            switch( aptr->elemno )
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


int SaveXYZMolecule( filename )
    char *filename;
{
    if( !Database )
        return( False );
    return( True );
}
 
 
int SaveCIFMolecule( filename )
    char *filename;
{
    if( !Database )
        return( False );
    return( True );
}

