/* abstree.c
 * RasMol2 Molecular Graphics
 * Roger Sayle, October 1994
 * Version 2.5
 */
#include "rasmol.h"

#ifdef IBMPC
#include <malloc.h>
#endif
#ifndef sun386
#include <stdlib.h>
#endif
#include <string.h>
#include <ctype.h>
#include <stdio.h>
#include <math.h>

#define ABSTREE
#include "molecule.h"
#include "abstree.h"


#define ExprPool    16
#define SetPool     4

typedef struct _SymEntry {
                struct _SymEntry __far *lft;
                struct _SymEntry __far *rgt;
                AtomSet __far *defn;
                char *ident;
                } SymEntry;

static SymEntry __far *SymbolTable;
static SymEntry __far *FreeEntry;
static AtomSet __far *FreeSet;
static Expr *FreeExpr;
static Expr FalseExpr;
static Expr TrueExpr;

/* Macros for commonly used loops */
#define ForEachAtom  for(chain=Database->clist;chain;chain=chain->cnext) \
                     for(group=chain->glist;group;group=group->gnext)    \
                     for(aptr=group->alist;aptr;aptr=aptr->anext)


#define BitAcidic       0x001
#define BitAliphatic    0x002
#define BitAromatic     0x004
#define BitBasic        0x008
#define BitBuried       0x010
#define BitCyclic       0x020
#define BitHydrophobic  0x040
#define BitMedium       0x080
#define BitNeutral      0x100
#define BitSmall        0x200

#define BitCharged      0x009
#define BitNotLarge     0x280

/* Acyclic = !Cyclic         */
/* Large = !Medium && !Small */
/* Polar = !Hydrophobic      */
/* Surface = !Buried         */


static int AminoProp[] = {
        /*ALA*/  BitAliphatic | BitBuried | BitHydrophobic | BitNeutral |
                 BitSmall,
        /*GLY*/  BitAliphatic | BitHydrophobic | BitNeutral | BitSmall,
        /*LEU*/  BitAliphatic | BitBuried | BitHydrophobic | BitNeutral,
        /*SER*/  BitNeutral | BitSmall,
        /*VAL*/  BitAliphatic | BitBuried | BitHydrophobic | BitMedium |
                 BitNeutral,
        /*THR*/  BitMedium | BitNeutral,
        /*LYS*/  BitBasic,
        /*ASP*/  BitAcidic | BitMedium,
        /*ILE*/  BitAliphatic | BitBuried | BitHydrophobic | BitNeutral,
        /*ASN*/  BitMedium | BitNeutral,
        /*GLU*/  BitAcidic,
        /*PRO*/  BitCyclic | BitHydrophobic | BitMedium | BitNeutral,
        /*ARG*/  BitBasic,
        /*PHE*/  BitAromatic | BitBuried | BitCyclic | BitHydrophobic |
                 BitNeutral,
        /*GLN*/  BitNeutral,
        /*TYR*/  BitAromatic | BitCyclic | BitNeutral,
        /*HIS*/  BitAromatic | BitBasic | BitCyclic | BitNeutral,
        /*CYS*/  BitBuried | BitMedium | BitNeutral,
        /*MET*/  BitBuried | BitHydrophobic | BitNeutral,
        /*TRP*/  BitAromatic | BitBuried | BitCyclic | BitHydrophobic |
                 BitNeutral,

        /*ASX*/  BitMedium | BitNeutral,
        /*GLX*/  BitNeutral,
        /*PCA*/  BitCyclic | BitHydrophobic | BitMedium | BitNeutral,
        /*HYP*/  BitCyclic | BitMedium | BitNeutral
        };


static void FatalExprError(ptr)
    char *ptr;
{
    char buffer[80];

    sprintf(buffer,"Expression Error: %s!",ptr);
    RasMolFatalExit(buffer);
}


#ifdef FUNCPROTO
/* Function Prototypes */
static AtomSet __far *SetInsert( AtomSet __far*, Atom __far* );
static int IsWithinRadius( AtomSet __far*, Long );
static int IsSetMember( AtomSet __far* );
#endif


static AtomSet __far *SetInsert( ptr, item )
    AtomSet __far *ptr;  Atom __far *item;
{
    register AtomSet __far *temp;
    register int i;

    if( ptr && (ptr->count<SetSize) )
    {   ptr->data[ptr->count] = item;
        ptr->count++;
        return( ptr );
    }

    if( !FreeSet )
    {   temp = (AtomSet __far*)_fmalloc( SetPool*sizeof(AtomSet) );
        if( !temp ) FatalExprError("Memory allocation failed");
        for( i=1; i<SetPool; i++ )
        {    temp->next = FreeSet;
             FreeSet = temp++;
        }
    } else
    {   temp = FreeSet;
        FreeSet = temp->next;
    }

    temp->data[0] = item;
    temp->next = ptr;
    temp->count = 1;
    return( temp );
}


static int IsWithinRadius( ptr, limit )
    AtomSet __far *ptr;
    Long limit;
{
    register Atom __far *aptr;
    register Long dx,dy,dz;
    register Long dist;
    register int i;

    while( ptr )
    {   for( i=0; i<ptr->count; i++ )
        {    aptr = ptr->data[i];
             dx = QAtom->xorg-aptr->xorg; if( (dist=dx*dx)>limit ) continue;
             dy = QAtom->yorg-aptr->yorg; if( (dist+=dy*dy)>limit ) continue;
             dz = QAtom->zorg-aptr->zorg; if( (dist+=dz*dz)>limit ) continue;
             return( True );
        }
        ptr = ptr->next;
    }
    return( False );
}


static int IsSetMember( ptr )
    AtomSet __far *ptr;
{
    register int i;

    while( ptr )
    {   for( i=0; i<ptr->count; i++ )
            if( ptr->data[i] == QAtom )
                return( True );
        ptr = ptr->next;
    }
    return( False );
}


void DeleteAtomSet( ptr )
    AtomSet __far *ptr;
{
    register AtomSet __far *temp;

    if( (temp = ptr) )
    {   while( temp->next )
            temp = temp->next;
        temp->next = FreeSet;
        FreeSet = ptr;
    }
}


Expr *AllocateNode()
{
    register Expr *ptr;
    register int i;

    if( !FreeExpr )
    {   ptr = (Expr*)malloc( ExprPool*sizeof(Expr) );
        if( !ptr ) FatalExprError("Memory allocation failed");
        for( i=1; i<ExprPool; i++ )
        {   ptr->rgt.ptr = FreeExpr;
            FreeExpr = ptr++;
        }
    } else
    {   ptr = FreeExpr;
        FreeExpr = ptr->rgt.ptr;
    }
    ptr->rgt.ptr = NULL;
    ptr->lft.ptr = NULL;
    return( ptr );
}


void DeAllocateExpr( expr )
    Expr *expr;
{
    if( !expr ) return;
    if( (expr == &FalseExpr) ||
        (expr == &TrueExpr) ) 
        return;

    if( expr->type!=OpWithin )
    {   if( !(expr->type&(OpLftProp|OpLftVal)) )
            DeAllocateExpr( expr->lft.ptr );
        if( !(expr->type&(OpRgtProp|OpRgtVal)) )
            DeAllocateExpr( expr->rgt.ptr );
    } else DeleteAtomSet( expr->rgt.set );
        
    expr->rgt.ptr = FreeExpr;
    FreeExpr = expr;
}


int GetElemNumber( aptr )
    Atom __far *aptr;
{
    register char *ptr;
    register char ch;

    ptr = ElemDesc[aptr->refno];
    ch = ptr[1];

    switch( *ptr )
    {   case(' '):  switch( ch )
                    {   case('B'):  return(  5 );
                        case('C'):  return(  6 );
                        case('D'):  return(  1 );
                        case('F'):  return(  9 );
                        case('H'):  return(  1 );
                        case('I'):  return( 53 );
                        case('K'):  return( 19 );
                        case('L'):  return(  1 );
                        case('N'):  return(  7 );
                        case('O'):  return(  8 );
                        case('P'):  return( 15 );
                        case('S'):  return( 16 );
                        case('U'):  return( 92 );
                        case('V'):  return( 23 );
                        case('W'):  return( 74 );
                        case('Y'):  return( 39 );
                    }
                    break;

        case('A'):  switch( ch )
                    {   case('C'):  return( 89 );
                        case('G'):  return( 47 );
                        case('L'):  return( 13 );
                        case('M'):  return( 95 );
                        case('R'):  return( 18 );
                        case('S'):  return( 33 );
                        case('T'):  return( 85 );
                        case('U'):  return( 79 );
                    }
                    break;

        case('B'):  switch( ch )
                    {   case('A'):  return( 56 );
                        case('E'):  return(  4 );
                        case('I'):  return( 83 );
                        case('K'):  return( 97 );
                        case('R'):  return( 35 );
                    }
                    break;

        case('C'):  switch( ch )
                    {   case('A'):  return( 20 );
                        case('D'):  return( 48 );
                        case('E'):  return( 58 );
                        case('F'):  return( 98 );
                        case('L'):  return( 17 );
                        case('M'):  return( 96 );
                        case('O'):  return( 27 );
                        case('R'):  return( 24 );
                        case('S'):  return( 55 );
                        case('U'):  return( 29 );
                    }
                    break;

        case('D'):  if( ch=='Y' )
                        return( 66 );
                    break;

        case('E'):  if( ch=='R' )
                    {   return( 68 );
                    } else if( ch=='S' )
                    {   return( 99 );
                    } else if( ch=='U' )
                        return( 63 );
                    break;

        case('F'):  if( ch=='E' )
                    {   return(  26 );
                    } else if( ch=='M' )
                    {   return( 100 );
                    } else if( ch=='R' )
                        return(  87 );
                    break;

        case('G'):  if( ch=='A' )
                    {   return( 31 );
                    } else if( ch=='D' )
                    {   return( 64 );
                    } else if( ch=='E' )
                        return( 32 );
                    break;

        case('H'):  if( ch=='E' )
                    {   return(  2 );
                    } else if( ch=='F' )
                    {   return( 72 );
                    } else if( ch=='G' )
                    {   return( 80 );
                    } else if( ch=='O' )
                        return( 67 );
                    break;

        case('I'):  if( ch=='N' )
                    {   return( 49 );
                    } else if( ch=='R' )
                        return( 77 );
                    break;

        case('K'):  if( ch=='R' )
                        return( 36 );
                    break;

        case('L'):  if( ch=='A' )
                    {   return(  57 );
                    } else if( ch=='I' )
                    {   return(   3 );
                    } else if( (ch=='R') || (ch=='W') )
                    {   return( 103 );
                    } else if( ch=='U' )
                        return(  71 );
                    break;

        case('M'):  if( ch=='D' )
                    {   return( 101 );
                    } else if( ch=='G' )
                    {   return(  12 );
                    } else if( ch=='N' )
                    {   return(  25 );
                    } else if( ch=='O' )
                        return(  42 );
                    break;

        case('N'):  switch( ch )
                    {   case('A'):  return(  11 );
                        case('B'):  return(  41 );
                        case('D'):  return(  60 );
                        case('E'):  return(  10 );
                        case('I'):  return(  28 );
                        case('O'):  return( 102 );
                        case('P'):  return(  93 );
                    }
                    break;

        case('O'):  if( ch=='S' )
                        return( 76 );
                    break;

        case('P'):  switch( ch )
                    {   case('A'):  return( 91 );
                        case('B'):  return( 82 );
                        case('D'):  return( 46 );
                        case('M'):  return( 61 );
                        case('O'):  return( 84 );
                        case('R'):  return( 59 );
                        case('T'):  return( 78 );
                        case('U'):  return( 94 );
                    }
                    break;

        case('R'):  switch( ch )
                    {   case('A'):  return( 88 );
                        case('B'):  return( 37 );
                        case('E'):  return( 75 );
                        case('H'):  return( 45 );
                        case('N'):  return( 86 );
                        case('U'):  return( 44 );
                    }
                    break;
                    break;

        case('S'):  switch( ch )
                    {   case('B'):  return( 51 );
                        case('C'):  return( 21 );
                        case('E'):  return( 34 );
                        case('I'):  return( 14 );
                        case('M'):  return( 62 );
                        case('N'):  return( 50 );
                        case('R'):  return( 38 );
                    }
                    break;

        case('T'):  switch( ch )
                    {   case('A'):  return( 73 );
                        case('B'):  return( 65 );
                        case('C'):  return( 43 );
                        case('E'):  return( 52 );
                        case('H'):  return( 90 );
                        case('I'):  return( 22 );
                        case('L'):  return( 81 );
                        case('M'):  return( 69 );
                    }
                    break;

        case('X'):  if( ch=='E' )
                        return( 54 );
                    break;

        case('Y'):  if( ch=='B' )
                        return( 70 );
                    break;

        case('Z'):  if( ch=='N' )
                    {   return( 30 );
                    } else if( ch=='R' )
                        return( 40 );
                    break;
    }

    if( (*ptr>='0') && (*ptr<='9') )
        if( (ch=='H') || (ch=='D') )
            return( 1 ); /* Hydrogen */

    return( 0 );
}


int ElemVDWRadius( elem )
    int elem;
{
    if( !HasHydrogen )
        switch( elem )
        {   case(  6 ):  return( VDWCarbon );
            case(  7 ):  return( VDWNitrogen );
            case(  8 ):  return( VDWOxygen );    
            case( 16 ):  return( VDWSulphur );
        }
    return( Element[elem].vdwrad );
}


static int EvaluateProperty( prop )
    int prop;
{
    switch( prop )
    {   case( PropIdent ):    return( QAtom->serno );
        case( PropXCord ):    return( (int)QAtom->xorg );
        case( PropYCord ):    return( (int)QAtom->yorg );
        case( PropZCord ):    return( (int)QAtom->zorg );
        case( PropTemp ):     return( QAtom->temp );
        case( PropName ):     return( QAtom->refno );
        case( PropResId ):    return( QGroup->serno );
        case( PropResName ):  return( QGroup->refno );
        case( PropChain ):    return( QChain->ident );
        case( PropSelect ):   return( QAtom->flag&SelectFlag );
        case( PropElemNo ):   return( GetElemNumber(QAtom) );
        case( PropRad ):      if( QAtom->flag&SphereFlag )
                              {   return( QAtom->radius );
                              } else return( 0 );

        
        /* Predicates stored in flags */
        case( PredBonded ):       return( !(QAtom->flag&NonBondFlag) );
        case( PredHydrogen ):     return( QAtom->flag&HydrogenFlag );
        case( PredHetero ):       return( QAtom->flag&HeteroFlag );
        case( PredCystine ):      return( QGroup->flag&CystineFlag );
        case( PredHelix ):        return( QGroup->struc&HelixFlag );
        case( PredSheet ):        return( QGroup->struc&SheetFlag );
        case( PredTurn ):         return( QGroup->struc&TurnFlag );

        /* Residue type predicates */
        case( PredDNA ):
        case( PredRNA ):
        case( PredNucleic ):      return( IsNucleo(QGroup->refno) );
        case( PredProtein ):      return( IsProtein(QGroup->refno) );
        case( PredAmino ):        return( IsAmino(QGroup->refno) );
        case( PredWater ):        return( IsWater(QGroup->refno) );
        case( PredSolvent ):      return( IsSolvent(QGroup->refno) );
        case( PredIon ):          return( IsIon(QGroup->refno) );

        /* General Predicates */
        case( PredAlpha ):        return( IsAmino(QGroup->refno) &&
                                          IsAlphaCarbon(QAtom->refno) );
        case( PredMainChain ):    return( (IsAmino(QGroup->refno) && 
                                           IsAminoBackbone(QAtom->refno)) ||
                                          (IsNucleo(QGroup->refno) &&
                                           IsNucleicBackbone(QAtom->refno)) );
        case( PredSidechain ):    return( IsAmino(QGroup->refno) &&
                                          !IsAminoBackbone(QAtom->refno) );
        case( PredLigand ):       return( (QAtom->flag&HeteroFlag) &&
                                          !IsSolvent(QGroup->refno) );

        /* Nucleic Acid Classifications */
        case( PredAT ):           return( IsAdenine(QGroup->refno) ||
                                          IsThymine(QGroup->refno) );
        case( PredCG ):           return( IsCytosine(QGroup->refno) ||
                                          IsGuanine(QGroup->refno) );
        case( PredPyrimidine ):   return( IsPyrimidine(QGroup->refno) );
        case( PredPurine ):       return( IsPurine(QGroup->refno) );


        /* Amino Acid Classifications */
        case( PredAcidic ):       return( IsAmino(QGroup->refno) &&
                                  AminoProp[QGroup->refno]&BitAcidic );

        case( PredAcyclic ):      return( IsAmino(QGroup->refno) &&
                                  !(AminoProp[QGroup->refno]&BitCyclic) );

        case( PredAliphatic ):    return( IsAmino(QGroup->refno) &&
                                  AminoProp[QGroup->refno]&BitAliphatic );

        case( PredAromatic ):     return( IsAmino(QGroup->refno) &&
                                  AminoProp[QGroup->refno]&BitAromatic );

        case( PredBasic ):        return( IsAmino(QGroup->refno) &&
                                  AminoProp[QGroup->refno]&BitBasic );

        case( PredBuried ):       return( IsAmino(QGroup->refno) &&
                                  AminoProp[QGroup->refno]&BitBuried );

        case( PredCharged ):      return( IsAmino(QGroup->refno) &&
                                  AminoProp[QGroup->refno]&BitCharged );

        case( PredCyclic ):       return( IsAmino(QGroup->refno) &&
                                  AminoProp[QGroup->refno]&BitCyclic );

        case( PredHydrophobic ):  return( IsAmino(QGroup->refno) &&
                                  AminoProp[QGroup->refno]&BitHydrophobic );
                 
        case( PredLarge ):        return( IsAmino(QGroup->refno) &&
                                  !(AminoProp[QGroup->refno]&BitNotLarge) );

        case( PredMedium ):       return( IsAmino(QGroup->refno) &&
                                  AminoProp[QGroup->refno]&BitMedium );

        case( PredNeutral ):      return( IsAmino(QGroup->refno) &&
                                  AminoProp[QGroup->refno]&BitNeutral );

        case( PredPolar ):        return( IsAmino(QGroup->refno) &&
                                  !(AminoProp[QGroup->refno]&BitHydrophobic) );

        case( PredSmall ):        return( IsAmino(QGroup->refno) &&
                                  AminoProp[QGroup->refno]&BitSmall );

        case( PredSurface ):      return( IsAmino(QGroup->refno) &&
                                  !(AminoProp[QGroup->refno]&BitBuried) );

    }
    return( True );
}


int EvaluateExpr( expr )
    Expr *expr;
{
    register int lft, rgt;

    if( !expr )
        return( True );

    if( expr->type==OpWithin )
    {   if( expr->lft.limit )
        {   return( IsWithinRadius(expr->rgt.set,expr->lft.limit) );
        } else return( IsSetMember(expr->rgt.set) );
    } else if( expr->type==OpMember )
        return( IsSetMember(expr->rgt.set) );

    if( expr->type & OpLftVal )
    {   lft = expr->lft.val;
    } else if( expr->type & OpLftProp )
    {   lft = EvaluateProperty( expr->lft.val );
    } else lft = EvaluateExpr( expr->lft.ptr );

    if( OpCode(expr)==OpConst ) return( lft );
    if( (OpCode(expr)==OpAnd) && !lft ) return( False );
    if( (OpCode(expr)==OpOr) && lft ) return( True );
    if( OpCode(expr)==OpNot ) return( !lft );

    if( expr->type & OpRgtVal )
    {   rgt = expr->rgt.val;
    } else if( expr->type & OpRgtProp )
    {   rgt = EvaluateProperty( expr->rgt.val );
    } else rgt = EvaluateExpr( expr->rgt.ptr );

    switch( OpCode(expr) )
    {   case(OpOr):
        case(OpAnd):     return( rgt );
        case(OpLess):    return( lft<rgt );
        case(OpMore):    return( lft>rgt );
        case(OpEqual):   return( lft==rgt );
        case(OpNotEq):   return( lft!=rgt );
        case(OpLessEq):  return( lft<=rgt );
        case(OpMoreEq):  return( lft>=rgt );
    }
    return( True );
}


AtomSet __far *BuildAtomSet( expr )
    Expr *expr;
{
    register AtomSet __far *ptr;

    ptr = (AtomSet __far*)0;

    if( Database )
        for( QChain=Database->clist; QChain; QChain=QChain->cnext )
            for( QGroup=QChain->glist; QGroup; QGroup=QGroup->gnext )
                for( QAtom=QGroup->alist; QAtom; QAtom=QAtom->anext )
                    if( EvaluateExpr(expr) )
                        ptr = SetInsert( ptr, QAtom );

    return( ptr );
}



int DefineSetExpr( ident, expr )
    char *ident;  Expr *expr;
{
    register SymEntry __far * __far *prev;
    register SymEntry __far *ptr;
    register AtomSet __far *set;
    register int result;

    result = True;
    prev = &SymbolTable;
    while( (ptr = *prev) )
    {   result = strcmp(ident,ptr->ident);
        if( !result ) break;  /* Entry Exists! */
        prev = (result<0)? &(ptr->lft) : &(ptr->rgt);
    }

    if( result )
    {   if( FreeEntry )
        {   ptr = FreeEntry;
            FreeEntry = ptr->rgt;
        } else /* Allocate SymEntry! */
            if( !(ptr = (SymEntry __far*)_fmalloc(sizeof(SymEntry))) )
                return( False );

        *prev = ptr;
        ptr->ident = ident;
        ptr->defn = (void __far*)0;
        ptr->lft = (void __far*)0;
        ptr->rgt = (void __far*)0;
    } else free(ident);

    if( expr )
    {   set = BuildAtomSet(expr);
        if( ptr->defn )
            DeleteAtomSet(ptr->defn);
        DeAllocateExpr(expr);
        ptr->defn = set;
    } else ptr->defn = (void __far*)0;
    return( True );
}


Expr *LookUpSetExpr( ident )
    char *ident;
{
    register SymEntry __far *ptr;
    register Expr *expr;
    register int result;

    result = True;
    ptr = SymbolTable;
    while( ptr )
    {   result = strcmp(ident,ptr->ident);
        if( !result ) break;  /* Entry Exists! */
        ptr = (result<0)? ptr->lft : ptr->rgt;
    }

    if( !result )
    {   expr = AllocateNode();
        expr->type = OpMember;
        expr->rgt.set = ptr->defn;
    } else expr = (Expr*)0;
    return( expr );
}


static int ElemCompare( ident, elem )
    char *ident, *elem;
{
    while( *elem )
        if( *elem++ != *ident++ )
            return( False );

    /* Handle Plurals */
    if( (ident[0]=='S') && !ident[1] )
        return( (elem[-1]!='S') && (elem[-1]!='Y') );
    return( !*ident );
}


Expr *LookUpElement( ident )
    char *ident;
{
    register Expr *expr;
    register int elem;

    for( elem=1; elem<MAXELEMNO; elem++ )
        if( ElemCompare(ident,Element[elem].name) )
            break;

    /* Handle Difficult Plurals & US Spelling! */
    if( elem == MAXELEMNO )
    {   if( *ident=='A' )
        {   if( ElemCompare(ident,"ALUMINUM") )
            {   elem = 13;
            } else if( !strcmp(ident,"ANTIMONIES") )
                elem = 51;
        } else if( *ident=='C' )
        {   if( ElemCompare(ident,"CESIUM") )
                elem = 55;
        } else if( *ident=='M' )
        {   if( !strcmp(ident,"MERCURIES") )
                elem = 80;
        } else if( *ident=='P' )
        {   if( !strcmp(ident,"PHOSPHORUSES") )
                elem = 8;
        } else if( *ident=='S' )
        {   if( ElemCompare(ident,"SULFUR") )
                elem = 16;
        }
    }

    if( elem<MAXELEMNO )
    {   expr = AllocateNode();
        expr->type = OpEqual|OpLftProp|OpRgtVal;
        expr->lft.val = PropElemNo;
        expr->rgt.val = elem;
    } else expr = (Expr*)0;
    return( expr );
}



static int MatchWildName( src, dst, size, len )
    char *src, *dst; int size, len;
{
    register int i, left;

    left = size;
    while( *dst==' ' )
    {   dst++; left--;
    }

    for( i=0; i<len; i++ )
    {   if( left )
        {   if( (*dst==*src) || (*src=='?') )
            {   dst++;  src++;  left--;
            } else return( False );
        } else if( *src++ != '?' )
            return( False );
    }

    while( left )
         if( *dst++!=' ' )
         {   return( False );
         } else left--;
    return( True );
}



int ParsePrimitiveExpr( orig )
    char **orig;
{
    static char NameBuf[4];
    register Expr *tmp1,*tmp2;
    register Expr *wild;
    register char *ptr;
    register int i, j;
    register int neg;
    register int ch;
    
    QueryExpr = &TrueExpr;
    ptr = *orig;
    ch = *ptr++;
    i = 0;

    if( ch != '*' )
    {   if( ch == '[' )
        {   i = 0;
            while( (ch = *ptr++) != ']' )
                if( ch && (i<3) )
                {   NameBuf[i++] = ToUpper(ch);
                } else return( False );
            ch = *ptr++;
        } else
            for( i=0; i<3; i++ )
                if( isalpha(ch) )
                {   NameBuf[i] = ToUpper(ch);
                    ch = *ptr++;
                } else if( (ch=='?') || (ch=='%') )
                {   NameBuf[i] = '?';
                    ch = *ptr++;
                } else break;
        if( !i ) return( False );

        wild = &FalseExpr;
        for( j=0; j<ResNo; j++ )
            if( MatchWildName(NameBuf,Residue[j],3,i) )
            {   tmp1 = AllocateNode();
                tmp1->type = OpEqual | OpLftProp | OpRgtVal;
                tmp1->lft.val = PropResName;
                tmp1->rgt.val = j;

                tmp2 = AllocateNode();
                tmp2->type = OpOr;
                tmp2->lft.ptr = tmp1;
                tmp2->rgt.ptr = wild;
                wild = tmp2;
            }
        QueryExpr = wild;
    } else ch = *ptr++;

    if( ch != '*' )
    {   if( ch == '-' )
        {   ch = *ptr++;
            neg = True;
        } else neg = False;

        if( isdigit(ch) )
        {   i = ch-'0';
            while( isdigit(*ptr) )
                i = 10*i + (*ptr++)-'0';

            tmp1 = AllocateNode();
            tmp1->type = OpEqual | OpLftProp | OpRgtVal;
            tmp1->rgt.val = neg? -i : i;
            tmp1->lft.val = PropResId;
            if( QueryExpr != &TrueExpr )
            {   tmp2 = AllocateNode();
                tmp2->type = OpAnd;
                tmp2->rgt.ptr = QueryExpr;
                tmp2->lft.ptr = tmp1;
                QueryExpr = tmp2;
            } else QueryExpr = tmp1;
            ch = *ptr++;
        } else if( neg )
            return( False );
    } else ch = *ptr++;

    if( ch==':' )
        ch = *ptr++;

    if( isalnum(ch) )
    {   ch = ToUpper(ch);

        tmp1 = AllocateNode();
        tmp1->type = OpEqual | OpLftProp | OpRgtVal;
        tmp1->lft.val = PropChain;
        tmp1->rgt.val = ch;
        if( QueryExpr != &TrueExpr )
        {   tmp2 = AllocateNode();
            tmp2->type = OpAnd;
            tmp2->rgt.ptr = QueryExpr;
            tmp2->lft.ptr = tmp1;
            QueryExpr = tmp2;
        } else QueryExpr = tmp1;
        ch = *ptr++;
    } else if( (ch=='?') || (ch=='%') || (ch=='*') )
        ch = *ptr++;

    if( ch == '.' )
    {   ch = *ptr++;
        if( ch!='*' )
        {   for( i=0; i<4; i++ )
                if( isalnum(ch) || ch=='\'' || ch=='*' )
                {   NameBuf[i] = ToUpper(ch);
                    ch = *ptr++;
                } else if( (ch=='?') || (ch=='%') || (ch=='#') )
                {   NameBuf[i] = '?';
                    ch = *ptr++;
                } else break;
            if( !i ) return( False );


            wild = &FalseExpr;
            for( j=0; j<ElemNo; j++ )
                if( MatchWildName(NameBuf,ElemDesc[j],4,i) )
                {   tmp1 = AllocateNode();
                    tmp1->type = OpEqual | OpLftProp | OpRgtVal;
                    tmp1->lft.val = PropName;
                    tmp1->rgt.val = j;

                    tmp2 = AllocateNode();
                    tmp2->type = OpOr;
                    tmp2->lft.ptr = tmp1;
                    tmp2->rgt.ptr = wild;
                    wild = tmp2;
                }

            if( (QueryExpr == &TrueExpr) || (wild == &FalseExpr) )
            {   DeAllocateExpr(QueryExpr);
                QueryExpr=wild;
            } else
            {   tmp1 = AllocateNode();
                tmp1->type = OpAnd;
                tmp1->lft.ptr = QueryExpr;
                tmp1->rgt.ptr = wild;
                QueryExpr = tmp1;
            }
        } else ch = *ptr++;
    }
    *orig = --ptr;
    return( !ch || isspace(ch) || ispunct(ch) );
}


static char *FormatInteger( ptr, value )
    char *ptr; int value;
{
    auto char buffer[10];
    register char *tmp;

    if( value<0 )
    {   value = -value;
        *ptr++ = '-';
    }

    if( value>9 )
    {   tmp = buffer;
        while( value>9 )
        {   *tmp++ = (value%10) + '0';
            value /= 10;
        }

        *ptr++ = value + '0';
        do { tmp--; 
            *ptr++ = *tmp;
        } while( tmp != buffer );
    } else *ptr++ = value + '0';
    return( ptr );
}


void FormatLabel( chain, group, aptr, label, ptr )
    Chain __far *chain;
    Group __far *group;
    Atom __far *aptr;
    char *label, *ptr;
{
    register char ch;
    register int i,j;

    while( ch = *label++ )
       if( ch=='%' )
       {   ch = *label++;
           if( isupper(ch) )
             ch = tolower(ch);

           switch( ch )
           {   case('a'):  /* Atom Name */
                           i = aptr->refno;
                           for( j=0; j<4; j++ )
                               if( ElemDesc[i][j]!=' ' )
                                   *ptr++ = ElemDesc[i][j];
                           break;

               case('b'):  /* Temperature/B-factor */
               case('t'):  ptr = FormatInteger(ptr,aptr->temp);
                           break;

               case('c'):  /* Chain Identifier */
               case('s'):  *ptr++ = chain->ident;
                           break;

               case('e'):  /* Element Type */
                           i = GetElemNumber(aptr);
                           *ptr++ = Element[i].symbol[0];
                           if( Element[i].symbol[1]!=' ' )
                               *ptr++ = Element[i].symbol[1];
                           break;

               case('i'):  /* Atom Number */
                           ptr = FormatInteger(ptr,aptr->serno);
                           break;

               case('n'):  /* Residue Name   */
                           i = group->refno;
                           for( j=0; j<3; j++ )
                               if( Residue[i][j]!=' ' )
                                   *ptr++ = Residue[i][j];
                           break;

               case('r'):  /* Residue Number */
                           ptr = FormatInteger(ptr,group->serno);
                           break;

               case('%'):  *ptr++ = '%';
                           break;
           }
       } else if( (ch>=' ') && (ch<='~') )
           *ptr++ = ch;
    *ptr = '\0';
}


#ifdef FUNCPROTO
/* Function Prototypes */
static void DeleteSymEntry( SymEntry __far* );
#endif


static void DeleteSymEntry( ptr )
    SymEntry __far *ptr;
{
    if( ptr->lft )
        DeleteSymEntry( ptr->lft );
    if( ptr->rgt )
        DeleteSymEntry( ptr->rgt );

    if( ptr->defn )
        DeleteAtomSet( ptr->defn );
    free( ptr->ident );

    ptr->rgt = FreeEntry;
    FreeEntry = ptr;
}


void ResetSymbolTable()
{
    if( SymbolTable )
    {   DeleteSymEntry(SymbolTable);
        SymbolTable = (void __far*)0;
    }
}


void InitialiseAbstree()
{
    FalseExpr.type = OpConst | OpLftVal | OpRgtVal;
    FalseExpr.rgt.val = FalseExpr.lft.val = 0;

    TrueExpr.type = OpConst | OpLftVal | OpRgtVal;
    TrueExpr.rgt.val = TrueExpr.lft.val = 1;

    QChain = (void __far*)0;
    QGroup = (void __far*)0;
    QAtom = (void __far*)0;

    SymbolTable = (void __far*)0;

    FreeEntry = (void __far*)0;
    FreeSet = (void __far*)0;
    FreeExpr = NULL;
}
