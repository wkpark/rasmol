/* abstree.h
 * RasMol2 Molecular Graphics
 * Roger Sayle, February 1994
 * Version 2.3
 */

#define OpCode(x) (((x)->type)&0x0f)

/* Operator Types */
#define OpAnd            0x01
#define OpOr             0x02
#define OpNot            0x03
#define OpEqual          0x04
#define OpNotEq          0x05
#define OpLess           0x06
#define OpMore           0x07
#define OpLessEq         0x08
#define OpMoreEq         0x09
#define OpConst          0x0a
#define OpWithin         0x0b
#define OpMember         0xac

#define OpLftProp        0x10
#define OpLftVal         0x20
#define OpRgtProp        0x40
#define OpRgtVal         0x80

/* Property fields */
#define PropIdent        1
#define PropXCord        2
#define PropYCord        3
#define PropZCord        4
#define PropTemp         5
#define PropRad          6
#define PropResId        7
#define PropName         8
#define PropChain        9
#define PropResName      10
#define PropSelect       11

#define PredAbsOrd(x)    ((x)-20)
#define PredAbsChr(x)    ((x)+20)

#define PredAlpha        20
#define PredAmino        21
#define PredAT           22
#define PredBackbone     23
#define PredBonded       24
#define PredCG           25
#define PredCystine      26
#define PredDNA          27
#define PredHelix        28
#define PredHetero       29
#define PredHydrogen     30
#define PredIon          31
#define PredLigand       32
#define PredNucleic      33
#define PredProtein      34
#define PredPurine       35
#define PredPyrimidine   36
#define PredRNA          37
#define PredSelected     38 /* Unused! */
#define PredSheet        39
#define PredSidechain    40
#define PredSolvent      41
#define PredTurn         42
#define PredWater	 43

#define PredAcidic       44
#define PredAcyclic      45
#define PredAliphatic    46
#define PredAromatic     47
#define PredBasic        48
#define PredBuried       49
#define PredCharged      50
#define PredCyclic       51
#define PredHydrophobic  52
#define PredLarge        53
#define PredMedium       54
#define PredNeutral      55
#define PredPolar        56
#define PredSmall        57
#define PredSurface      58


#define MAXHANDLE   28
#define HandAg       1
#define HandAl       2
#define HandAu       3
#define HandB        4
#define HandBr       5
#define HandC        6
#define HandCa       7
#define HandCl       8
#define HandCr       9
#define HandCu      10
#define HandF       11
#define HandFe      12
#define HandH       13
#define HandHe      14
#define HandI       15
#define HandK       16
#define HandLi      17
#define HandMg      18
#define HandMn      19
#define HandN       20
#define HandNa      21
#define HandNi      22
#define HandO       23
#define HandP       24
#define HandS       25
#define HandSi      26
#define HandZn      27


#define SetSize     10
typedef struct _AtomSet {
	struct _AtomSet __far *next;
	Atom __far *data[SetSize];
        int count;
        } AtomSet;
        
typedef union {
	AtomSet __far *set;
	struct _Expr *ptr;
        Long limit;
	int val;
	} Branch;

typedef struct _Expr {
	int type;
        Branch rgt;
        Branch lft;
	} Expr;


#ifdef ABSTREE
int VanWaalRadius[MAXHANDLE] = { 360, /* Unknown */
     318, /*Ag*/  425, /*Al*/  342, /*Au*/  250, /* B*/
     488, /*Br*/  468, /* C*/  425, /*Ca*/  450, /*Cl*/
     342, /*Cr*/  180, /*Cu*/  350, /* F*/  425, /*Fe*/
     250, /* H*/  300, /*He*/  538, /* I*/  360, /* K*/
     465, /*Li*/  218, /*Mg*/  342, /*Mn*/  375, /* N*/
     280, /*Na*/  172, /*Ni*/  350, /* O*/  425, /* P*/
     462, /* S*/  482, /*Si*/  185  /*Zn*/ };

Expr *QueryExpr;
Chain __far *QChain;
Group __far *QGroup;
Atom __far *QAtom;

#else
extern int VanWaalRadius[MAXHANDLE];

extern Expr *QueryExpr;
extern Chain __far *QChain;
extern Group __far *QGroup;
extern Atom __far *QAtom;

#ifdef __STDC__
Expr *AllocateNode();
void DeAllocateExpr( Expr* );
int GetElemIdent( Atom __far* );
int EvaluateExpr( Expr* );
int DefineSetExpr( char*, Expr* );
Expr *LookUpSetExpr( char* );
AtomSet __far *BuildAtomSet( Expr* );
void DeleteAtomSet( AtomSet __far* );
int ParsePrimitiveExpr( char** );
void InitialiseAbstree();
void ResetSymbolTable();

#else /* non-ANSI C compiler */
Expr *AllocateNode();
void DeAllocateExpr();
int GetElemIdent();
int EvaluateExpr();
int DefineSetExpr();
Expr *LookUpSetExpr();
AtomSet __far *BuildAtomSet();
void DeleteAtomSet();
int ParsePrimitiveExpr();
void InitialiseAbstree();
void ResetSymbolTable();

#endif
#endif

