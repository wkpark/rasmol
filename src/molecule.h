/***************************************************************************
 *                              RasMol 2.7.1                               *
 *                                                                         *
 *                                 RasMol                                  *
 *                 Molecular Graphics Visualisation Tool                   *
 *                              22 June 1999                               *
 *                                                                         *
 *                   Based on RasMol 2.6 by Roger Sayle                    *
 * Biomolecular Structures Group, Glaxo Wellcome Research & Development,   *
 *                      Stevenage, Hertfordshire, UK                       *
 *         Version 2.6, August 1995, Version 2.6.4, December 1998          *
 *                   Copyright (C) Roger Sayle 1992-1999                   *
 *                                                                         *
 *                  and Based on Mods by Arne Mueller                      *
 *                      Version 2.6x1, May 1998                            *
 *                   Copyright (C) Arne Mueller 1998                       *
 *                                                                         *
 *           Version 2.7.0, 2.7.1 Mods by Herbert J. Bernstein             *
 *           Bernstein + Sons, P.O. Box 177, Bellport, NY, USA             *
 *                      yaya@bernstein-plus-sons.com                       *
 *                    2.7.0 March 1999, 2.7.1 June 1999                    *
 *              Copyright (C) Herbert J. Bernstein 1998-1999               *
 *                                                                         *
 * Please read the file NOTICE for important notices which apply to this   *
 * package. If you are not going to make changes to RasMol, you are not    *
 * only permitted to freely make copies and distribute them, you are       *
 * encouraged to do so, provided you do the following:                     *
 *   * 1. Either include the complete documentation, especially the file   *
 *     NOTICE, with what you distribute or provide a clear indication      *
 *     where people can get a copy of the documentation; and               *
 *   * 2. Please give credit where credit is due citing the version and    *
 *     original authors properly; and                                      *
 *   * 3. Please do not give anyone the impression that the original       *
 *     authors are providing a warranty of any kind.                       *
 *                                                                         *
 * If you would like to use major pieces of RasMol in some other program,  *
 * make modifications to RasMol, or in some other way make what a lawyer   *
 * would call a "derived work", you are not only permitted to do so, you   *
 * are encouraged to do so. In addition to the things we discussed above,  *
 * please do the following:                                                *
 *   * 4. Please explain in your documentation how what you did differs    *
 *     from this version of RasMol; and                                    *
 *   * 5. Please make your modified source code available.                 *
 *                                                                         *
 * This version of RasMol is not in the public domain, but it is given     *
 * freely to the community in the hopes of advancing science. If you make  *
 * changes, please make them in a responsible manner, and please offer us  *
 * the opportunity to include those changes in future versions of RasMol.  *
 ***************************************************************************/

/* molecule.h
 */
#define MAXMASK 40
#define MAXELEM 1024
#define MINELEM 29
#define MAXRES  100
#define MINRES  54
#define CIS     90  /* max. omega-angle to form a cis-peptidbond */


#define IsAmino(x)       ((x)<=23)
#define IsAminoNucleo(x) ((x)<=42)
#define IsNucleo(x)      (((x)>=24) && ((x)<=42))
#define IsProtein(x)     (((x)<=23) || (((x)>=43) && ((x)<=45)))
#define IsDNA(x)         (((x)>=24) && ((x)<=27))
#define IsSolvent(x)     (((x)>=46) && ((x)<=49))
#define IsWater(x)       (((x)==46) || ((x)==47))
#define IsIon(x)         (((x)==48) || ((x)==49))

#define IsPyrimidine(x)  (IsCytosine(x) || IsThymine(x))
#define IsPurine(x)      (IsAdenine(x) || IsGuanine(x))
#define IsRNA(x)         (IsNucleo(x) && !IsThymine(x))
#define NucleicCompl(x)  ((x)^3)


#define IsProline(x)     ((x)==11)
#define IsCysteine(x)    ((x)==17)
#define IsAdenine(x)     ((x)==24)
#define IsCytosine(x)    ((x)==25)
#define IsGuanine(x)     ((x)==26)
#define IsThymine(x)     ((x)==27)



#define IsAlphaCarbon(x)     ((x)==1)
#define IsSugarPhosphate(x)  ((x)==7)
#define IsAminoBackbone(x)   ((x)<=3)
#define IsShapelyBackbone(x) ((x)<=7)
#define IsNucleicBackbone(x) (((x)>=7) && ((x)<=18))
#define IsShapelySpecial(x)  ((x)==19)
#define IsCysteineSulphur(x) ((x)==20)
#define IsCoenzyme(x)        (((x)>=50) && ((x)<=53))


/*=================*/
/*  Database Flags */
/*=================*/

#define SelectFlag      0x01
#define DrawBondFlag    0x0e
#define AllAtomFlag     0x1c
#define HelixFlag       0x03
#define DrawKnotFlag    0x7e
#define WideKnotFlag    0x0e

/* Atom Flags */
#define SphereFlag      0x02     /* Sphere representation */
#define HeteroFlag      0x04     /* HETATM record         */
#define HydrogenFlag    0x08     /* Hydrogen atom         */
#define NormAtomFlag    0x10
#define NonBondFlag     0x20
#define BreakFlag       0x40     /* Break in backbone     */
#define StarFlag        0x80     /* Star representation   */

/* Bond Flags */
#define WireFlag        0x02     /* Depth-cued wireframe         */
#define DashFlag        0x04     /* Dashed Depth-cued wireframe  */
#define CylinderFlag    0x08     /* Line/Cylinder representation */

#define HydrBondFlag    0x00     /* Hydrogen bond [place keeper] */
#define NormBondFlag    0x10
#define DoubBondFlag    0x20
#define TripBondFlag    0x40
#define AromBondFlag    0x80


/* Group Flags */
#define CystineFlag     0x01     /* Disulphide bonded cysteine  */
#define StrandFlag      0x02     /* Strands representation      */
#define DashStrandFlag  0x04     /* Dash Strands representation */
#define RibbonFlag      0x08     /* Solid Ribbon representation */
#define TraceFlag       0x10     /* Smooth trace representation */
#define CartoonFlag     0x20     /* Richardson protein cartoon  */
#define DotsFlag        0x40     /* Dotted trace representation */
#define CisBondFlag     0x80     /* Cis bonded residue          */

/* Structure Flags */
#define Helix3Flag      0x01     /* 3,10-Helix structure       */
#define Helix4Flag      0x02     /* Alpha Helix structure      */
#define Helix5Flag      0x03     /* 5-Helix structure          */
#define SheetFlag       0x04     /* Beta Sheet structure       */
#define TurnFlag        0x08     /* Turn Secondary structure   */


/*=====================*/
/*  Molecule Database  */
/*=====================*/

typedef struct _Atom {
        struct _Atom __far *anext;        /* Linked list of atoms  */
        struct _Atom __far *bucket;       /* Sphere Y-Bucket       */
        struct _Atom __far *next;         /* Active Object List    */
        Long   xorg, yorg, zorg;          /* World Co-ordinates    */
        short  xtrl, ytrl, ztrl;          /* Trailing Bits         */
        short  x, y, z;                   /* Image Co-ordinates    */
        short  radius;                    /* World Radius          */
        short  temp;                      /* Temperature Factor    */
        short  col;                       /* Atom Colour           */
        Long   serno;                     /* Atom Serial Number    */
        void   *label;                    /* Atom Label Structure  */
        Byte   elemno;                    /* Atomic Number         */
        Byte   refno;                     /* ElemDesc index number */
        Byte   flag;                      /* Database flags        */
        char   altl;                      /* Alternate Location    */
        short  irad;                      /* Image Radius          */
        short  mbox;                      /* Shadow Casting NOnce  */
        short  model;                     /* Atom Model Number     */
    } Atom;


typedef struct _Bond {
        struct _Bond __far *bnext;       /* Linked list of bonds  */
        Atom __far *srcatom;             /* Source Atom Ptr       */
        Atom __far *dstatom;             /* Destination Atom Ptr  */
        short radius;                    /* World Radius          */
        short irad;                      /* Image Radius          */
        short col;                       /* Bond Colour           */
        Byte  flag;                      /* Database flags        */
        char  altl;                      /* Bond Alternate Loc    */
    } Bond;

typedef struct _Group {
        struct _Group __far *gnext;       /* Linked list of groups */
        Atom __far *alist;                /* Linked list of atoms  */
        short serno;                      /* Group serial number   */
        short sserno;                     /* Secondary serial no.  */ 
        short width;                      /* Ribbon Width          */
        short col1;                       /* Ribbon Colour #1      */
        short col2;                       /* Ribbon Colour #2      */
        char  insert;                     /* PDB insertion code    */
        char  sinsert;                    /* Secondary insert code */
        Byte  refno;                      /* Residue index number  */
        Byte  struc;                      /* Secondary Structure   */
        Byte  flag;                       /* Database flags        */
        short model;                      /* Group Model Number    */
    } Group;
 
#ifdef APPLEMAC
/* Avoid Name Clash! */
#define Chain ChainSeg
#endif

typedef struct _ChainSeg {
        struct _ChainSeg __far *cnext;    /* Linked list of chains     */
        Group __far *glist;               /* Linked list of groups     */
        Bond __far *blist;                /* Linked list of back bonds */
        char ident;                       /* Chain identifier          */
        Byte model;                       /* NMR Model / Symmetry      */
    } Chain;

typedef struct _AtomRef {
        Chain __far *chn;
        Group __far *grp;
        Atom  __far *atm;
    } AtomRef;

typedef struct _HBond {
        struct _HBond __far *hnext;       /* Ordered list of hbonds   */
        Atom __far *srcCA;                /* Source Alpha Carbon      */
        Atom __far *dstCA;                /* Destination Alpha Carbon */
        Atom __far *dst;                  /* Acceptor [=CO] Atom Ptr  */
        Atom __far *src;                  /* Donor [=NH] Atom Ptr     */
        short energy;                     /* Hydrogen bond energy     */
        short radius;                     /* World Radius             */
        short irad;                       /* Image Radius             */
        Char offset;                      /* Signed Offset            */
        Byte flag;                        /* Database flags           */
        Byte col;                         /* Hydrogen bond colour     */
        char altl;                        /* Bond Alternate Loc       */
    } HBond;

typedef struct _Molecule {
        HBond __far *slist;               /* Linked list of SS bonds  */
        HBond __far *hlist;               /* Linked list of hbonds    */
        Chain __far *clist;               /* Linked list of chains    */
        Bond __far *blist;                /* Linked list of bonds     */
    } Molecule;



/*========================*/
/* Other Consts & Structs */
/*========================*/

#define SourceNone   0
#define SourcePDB    1
#define SourceCalc   2
 
#define SerNoFlag 0x01
#define ResNoFlag 0x02

typedef struct {
        short radius;
        char  mask[19];
        Byte  flags;
        Byte  r;
        Byte  g;
        Byte  b;
        } MaskDesc;

typedef struct _IntCoord {
        struct _IntCoord __far *inext;
        short na,nb,nc;
        short refno;
        Real dihed;
        Real angle;
        Real dist;
    } IntCoord;

typedef struct _InfoStruct {
        char filename[1024];
        char moleculename[80];
        char classification[42];
        char date[12];
        char technique[80];
        char identcode[80];

        char spacegroup[12];
        Real cellalpha, cellbeta, cellgamma;
        Real cella, cellb, cellc;

        double vecf2o[3], veco2f[3], matf2o[3][3], mato2f[3][3];
        double cell[6];

        Long bondcount;
        int chaincount;
        int ssbondcount;
        int hbondcount;
        int cisbondcount;

        int structsource;
        int laddercount;
        int helixcount;
        int turncount;
    } InfoStruct;

#ifdef APPLEMAC
void RegisterAlloc(void __far * );
#else
#define RegisterAlloc(x)
#endif
void FreeAlloc(void __far * );

/* used to describe an defined part of the selected molecule */
typedef enum{NO, ATM, CRD, GRP, CHN} Selection;


#ifdef MOLECULE
/* Avoid SGI Compiler Warnings! */
char Residue[MAXRES][4] = {
    /*===============*/
    /*  Amino Acids  */
    /*===============*/

/* Ordered by Cumulative Frequency in Brookhaven *
 * Protein Databank, December 1991               */

          "ALA", /* 8.4% */     "GLY", /* 8.3% */
          "LEU", /* 8.0% */     "SER", /* 7.5% */
          "VAL", /* 7.1% */     "THR", /* 6.4% */
          "LYS", /* 5.8% */     "ASP", /* 5.5% */
          "ILE", /* 5.2% */     "ASN", /* 4.9% */
          "GLU", /* 4.9% */     "PRO", /* 4.4% */
          "ARG", /* 3.8% */     "PHE", /* 3.7% */
          "GLN", /* 3.5% */     "TYR", /* 3.5% */
          "HIS", /* 2.3% */     "CYS", /* 2.0% */
          "MET", /* 1.8% */     "TRP", /* 1.4% */

          "ASX", "GLX", "PCA", "HYP",

    /*===================*/
    /*  DNA Nucleotides  */
    /*===================*/
          "  A", "  C", "  G", "  T",

    /*===================*/
    /*  RNA Nucleotides  */
    /*===================*/
          "  U", " +U", "  I", "1MA", 
          "5MC", "OMC", "1MG", "2MG", 
          "M2G", "7MG", "OMG", " YG", 
          "H2U", "5MU", "PSU",

    /*=================*/
    /*  Miscellaneous  */ 
    /*=================*/
          "UNK", "ACE", "FOR", "HOH",
          "DOD", "SO4", "PO4", "NAD",
          "COA", "NAP", "NDP"  };


/* Avoid SGI Compiler Warnings! */
char ElemDesc[MAXELEM][4] = {
    { ' ', 'N', ' ', ' ' },  /* 0*/
    { ' ', 'C', 'A', ' ' },  /* 1*/
    { ' ', 'C', ' ', ' ' },  /* 2*/
    { ' ', 'O', ' ', ' ' },  /* 3*/   /* 0-3   Amino Acid Backbone    */
    { ' ', 'C', '\'', ' ' }, /* 4*/
    { ' ', 'O', 'T', ' ' },  /* 5*/
    { ' ', 'S', ' ', ' ' },  /* 6*/
    { ' ', 'P', ' ', ' ' },  /* 7*/   /* 4-7   Shapely Amino Backbone */
    { ' ', 'O', '1', 'P' },  /* 8*/
    { ' ', 'O', '2', 'P' },  /* 9*/
    { ' ', 'O', '5', '*' },  /*10*/
    { ' ', 'C', '5', '*' },  /*11*/
    { ' ', 'C', '4', '*' },  /*12*/
    { ' ', 'O', '4', '*' },  /*13*/
    { ' ', 'C', '3', '*' },  /*14*/
    { ' ', 'O', '3', '*' },  /*15*/
    { ' ', 'C', '2', '*' },  /*16*/
    { ' ', 'O', '2', '*' },  /*17*/
    { ' ', 'C', '1', '*' },  /*18*/   /* 7-18  Nucleic Acid Backbone  */
    { ' ', 'C', 'A', '2' },  /*19*/   /* 19    Shapely Special        */
    { ' ', 'S', 'G', ' ' },  /*20*/   /* 20    Cysteine Sulphur       */
    { ' ', 'N', '1', ' ' },  /*21*/
    { ' ', 'N', '2', ' ' },  /*22*/
    { ' ', 'N', '3', ' ' },  /*23*/
    { ' ', 'N', '4', ' ' },  /*24*/
    { ' ', 'N', '6', ' ' },  /*25*/
    { ' ', 'O', '2', ' ' },  /*26*/
    { ' ', 'O', '4', ' ' },  /*27*/
    { ' ', 'O', '6', ' ' }   /*28*/   /* 21-28 Nucleic Acid H-Bonding */
    };


InfoStruct Info;
int MainGroupCount;
int HetaGroupCount;
Long MainAtomCount; 
Long HetaAtomCount;
int CisBondCutOff;

Long MinX, MinY, MinZ;
Long MaxX, MaxY, MaxZ;

int HMinMaxFlag, MMinMaxFlag;
int MinMainTemp, MaxMainTemp;
int MinHetaTemp, MaxHetaTemp;
int MinMainRes,  MaxMainRes;
int MinHetaRes,  MaxHetaRes;
int MinAltl,     MaxAltl;

short MinModel, MaxModel;

Molecule __far *CurMolecule;
Chain __far *CurChain;
Group __far *CurGroup;
Atom __far *CurAtom;

IntCoord __far *IntList;
Molecule __far *Database;
MaskDesc UserMask[MAXMASK];
Long MinHBondDist, MaxHBondDist;
Long MinBondDist,  MaxBondDist;
int ElemNo,ResNo;
int HasHydrogen;
int MaskCount;
int NMRModel;
int NullBonds;
int MarkAtoms;

int HBondChainsFlag;

#else
extern char Residue[MAXRES][4];
extern char ElemDesc[MAXELEM][4];
extern InfoStruct Info;

extern int MainGroupCount;
extern int HetaGroupCount;
extern Long MainAtomCount;
extern Long HetaAtomCount;
extern int CisBondCutOff;

extern Long MinX, MinY, MinZ;
extern Long MaxX, MaxY, MaxZ;

extern int HMinMaxFlag, MMinMaxFlag;
extern int MinMainTemp, MaxMainTemp;
extern int MinHetaTemp, MaxHetaTemp;
extern int MinMainRes,  MaxMainRes;
extern int MinHetaRes,  MaxHetaRes;
extern int MinAltl,     MaxAltl;

extern short MinModel, MaxModel;

extern Molecule __far *CurMolecule;
extern Chain __far *CurChain;
extern Group __far *CurGroup;
extern Atom __far *CurAtom;

extern IntCoord __far *IntList;
extern Molecule __far *Database;
extern MaskDesc UserMask[MAXMASK];
extern Long MinHBondDist, MaxHBondDist;
extern Long MinBondDist,  MaxBondDist;
extern int ElemNo,ResNo;
extern int HasHydrogen;
extern int MaskCount;
extern int NMRModel;
extern int NullBonds;
extern int MarkAtoms;

extern int HBondChainsFlag;
#ifndef APPLEMAC
#define RegisterAlloc(x)
#endif
#endif

void CreateChain( int );
void CreateGroup( int );
void ProcessGroup( int );
void CreateMolGroup( void );
void CreateNextMolGroup( void );
int FindResNo( char* );

Atom __far *CreateAtom( void );
Atom __far *FindGroupAtom( Group __far*, int );
void ProcessAtom( Atom __far* );

int NewAtomType( char* );
int SimpleAtomType( char* );
int ComplexAtomType( char* );

Bond __far *ProcessBond( Atom __far*, Atom __far*, int );
void CreateBond( Long, Long, int );
void CreateBondOrder( Long, Long );
void CreateNewBond( Long, Long );
void CreateMoleculeBonds( int, int, int );
Atom __far *FindCysSulphur( Group __far *group );
void FindDisulphideBridges( void );
void FindCisBonds( void );
void CalcHydrogenBonds( void );

void InitInternalCoords( void );
IntCoord __far* AllocInternalCoord( void );
int ConvertInternal2Cartesian( void );
void FreeInternalCoords( void );

void DetermineStructure( int );
void RenumberMolecule( int );

void InitialiseDatabase( void );
void ReviseTitle( void );
void DescribeMolecule( void );
void DestroyDatabase( void );
void PurgeDatabase( void );
#ifdef APPLEMAC
void RegisterAlloc( void *);
#endif

