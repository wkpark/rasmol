/* abstree.c
 * RasMol2 Molecular Graphics
 * Roger Sayle, February 1994
 * Version 2.3
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


#if defined(__STDC__) || defined(IBMPC)
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


int GetElemIdent( aptr )
    Atom __far *aptr;
{
    register char *ptr;

    ptr = ElemDesc[aptr->refno];
    switch( *ptr )
    {   case(' '):  switch( ptr[1] )
                    {   case('B'): return( HandB );
                        case('C'): return( HandC );
                        case('D'): return( HandH );
                        case('F'): return( HandF );
                        case('H'): return( HandH );
                        case('I'): return( HandI );
                        case('K'): return( HandK );
                        case('N'): return( HandN );
                        case('O'): return( HandO );
                        case('P'): return( HandP );
                        case('S'): return( HandS );
                    }
                    break;

        case('A'):  if( ptr[1]=='G' )
                    {   return( HandAg );
                    } else if( ptr[1]=='L' )
                    {   return( HandAl );
                    } else if( ptr[2]=='U' )
                        return( HandAu );
                    break;

        case('B'):  if( ptr[1]=='R' )
                        return( HandBr );
                    break;

        case('C'):  switch( ptr[1] )
                    {   case('A'): return( HandCa );
                        case('L'): return( HandCl );
                        case('R'): return( HandCr );
                        case('U'): return( HandCu );
                    }
                    break;

        case('F'):  if( ptr[1]=='E' )
                        return( HandFe );
                    break;

        case('H'):  if( ptr[1]=='E' )
                        return( HandHe );
                    break;

        case('L'):  if( ptr[1]=='I' )
                        return( HandLi );
                    break;

        case('M'):  if( ptr[1]=='G' )
                    {   return( HandMg );
                    } else if( ptr[1]=='N' )
                        return( HandMn );
                    break;

        case('N'):  if( ptr[1]=='A' )
                    {   return( HandNa );
                    } else if( ptr[1]=='I' )
                        return( HandNi );
                    break;

        case('S'):  if( ptr[1]=='I' )
                        return( HandSi );
                    break;

        case('Z'):  if( ptr[1]=='N' )
                        return( HandZn );
                    break;
    }

    if( (*ptr>='0') && (*ptr<='9') )
        if( (ptr[1]=='H') || (ptr[1]=='D') )
            return( HandH );

    return( 0 );
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
        case( PropRad ):      if( QAtom->flag&SphereFlag )
                              {   return( QAtom->radius );
                              } else return( 0 );

        
        /* Predicates stored in flags */
        case( PredBonded ):       return( !(QAtom->flag&NonBondFlag) );
        case( PredHydrogen ):     return( QAtom->flag&HydrogenFlag );
        case( PredHetero ):       return( QAtom->flag&HeteroFlag );
        case( PredCystine ):      return( QGroup->flag&CystineFlag );
        case( PredHelix ):        return( QGroup->flag&HelixFlag );
        case( PredSheet ):        return( QGroup->flag&SheetFlag );
        case( PredTurn ):         return( QGroup->flag&TurnFlag );

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
        case( PredBackbone ):     return( (IsAmino(QGroup->refno) && 
                                           IsAminoBackbone(QAtom->refno)) ||
                                          (IsNucleo(QGroup->refno) &&
                                           IsNucleicBackbone(QAtom->refno)) );
        case( PredLigand ):       return( (QAtom->flag&HeteroFlag) &&
                                          !IsSolvent(QGroup->refno) );
        case( PredSidechain ):    return( IsAmino(QGroup->refno) &&
                                          !IsAminoBackbone(QAtom->refno) );

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
                {   NameBuf[i++] = islower(ch)? toupper(ch) : ch;
                } else return( False );
            ch = *ptr++;
        } else
            for( i=0; i<3; i++ )
                if( isalpha(ch) )
                {   NameBuf[i] = islower(ch)? toupper(ch) : ch;
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
    {   if( islower(ch) )
            ch = toupper(ch);

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
                {   NameBuf[i] = islower(ch)? toupper(ch) : ch;
                    ch = *ptr++;
                } else if( (ch=='?') || (ch=='%') )
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


#if defined(__STDC__) || defined(IBMPC)
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
