/* molecule.h
 * RasMol2 Molecular Graphics
 * Roger Sayle, October 1994
 * Version 2.5
 */
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

#define MAXMASK 40
#define MAXELEM 256
#define MINELEM 29
#define MAXRES  80
#define MINRES  35


#define IsAmino(x)       ((x)<=23)
#define IsAminoNucleo(x) ((x)<=27)
#define IsNucleo(x)      (((x)>=24) && ((x)<=27))
#define IsProtein(x)     (((x)<=23) || (((x)>=28) && ((x)<=30)))
#define IsSolvent(x)     (((x)>=31) && ((x)<=34))
#define IsWater(x)       (((x)==31) || ((x)==32))
#define IsIon(x)         (((x)==33) || ((x)==34))

#define IsPyrimidine(x)  (IsCytosine(x) || IsThymine(x))
#define IsPurine(x)      (IsAdenine(x) || IsGuanine(x))
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


/*=================*/
/*  Database Flags */
/*=================*/

#define SelectFlag      0x01
#define DrawBondFlag    0x06
#define AllAtomFlag     0x1c
#define HelixFlag       0x03

/* Atom Flags */
#define SphereFlag      0x02     /* Sphere representation */
#define HeteroFlag      0x04     /* HETATM record         */
#define HydrogenFlag    0x08     /* Hydrogen atom         */
#define NormAtomFlag    0x10
#define NonBondFlag     0x20
#define BreakFlag       0x40     /* Break in backbone     */

/* Bond Flags */
#define WireFlag        0x02     /* Depth-cued wireframe         */
#define CylinderFlag    0x04     /* Line/Cylinder representation */
#define HydrBondFlag    0x08     /* Hydrogen *-H bond            */
#define NormBondFlag    0x10
#define DoubBondFlag    0x20
#define TripBondFlag    0x40
#define AromBondFlag    0x80

/* Group Flags */
#define CystineFlag     0x01     /* Disulphide bonded cysteine  */
#define RibbonFlag      0x02     /* Solid Ribbon representation */
#define StrandFlag      0x04     /* Strands representation      */

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
        short  x, y, z;                   /* Image Co-ordinates    */
        short  radius;                    /* World Radius          */
        short  serno;                     /* Atom Serial Number    */
        short  temp;                      /* Temperature Factor    */
        short  col;                       /* Atom Colour           */
        void   *label;                    /* Atom Label Structure  */
        Byte   refno;                     /* ElemDesc index number */
        Byte   flag;                      /* Database flags        */
        char   altl;                      /* Alternate Location    */
        short  irad;                      /* Image Radius          */
        short  mbox;                      /* Shadow Casting NOnce  */
	} Atom;


typedef struct _Bond {
        struct _Bond __far *bnext;       /* Linked list of bonds  */
        Atom __far *srcatom;             /* Source Atom Ptr       */
        Atom __far *dstatom;             /* Destination Atom Ptr  */
        short radius;                    /* World Radius          */
        short irad;                      /* Image Radius          */
        short col;                       /* Bond Colour           */
        Byte  flag;                      /* Database flags        */
	} Bond;

typedef struct _Group {
        struct _Group __far *gnext;       /* Linked list of groups */
        Atom __far *alist;                /* Linked list of atoms  */
        short serno;                      /* Group serial number   */
        short width;                      /* Ribbon Width          */
        short col1;			  /* Ribbon Colour #1      */
        short col2;			  /* Ribbon Colour #2      */
        Byte  refno;                      /* Residue index number  */
        Byte  struc;                      /* Secondary Structure   */
        Byte  flag;                       /* Database flags        */
	} Group;
 
#ifdef APPLEMAC
/* Avoid Name Clash! */
#define Chain ChainSeg
#endif

typedef struct _ChainSeg {
        struct _ChainSeg __far *cnext;       /* Linked list of chains     */
        Group __far *glist;               /* Linked list of groups     */
        Bond __far *blist;                /* Linked list of back bonds */
        char  ident;                      /* Chain identifier          */
	} Chain;

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
        } HBond;

typedef struct _Molecule {
        HBond __far *slist;               /* Linked list of SS bonds  */
        HBond __far *hlist;               /* Linked list of hbonds    */
        Chain __far *clist;               /* Linked list of chains    */
        Bond __far *blist;                /* Linked list of bonds     */
	} Molecule;



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

    /*===============*/
    /*  Nucleotides  */
    /*===============*/
          "  A", "  C", "  G", "  T",

    /*=================*/
    /*  Miscellaneous  */ 
    /*=================*/
          "UNK", "ACE", "FOR", "HOH",
          "DOD", "SO4", "PO4"  };


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



char InfoFileName[256];
char InfoClassification[42];
char InfoMoleculeName[80];
char InfoSpaceGroup[11];
char InfoIdentCode[6];

Real InfoCellAlpha, InfoCellBeta, InfoCellGamma;
Real InfoCellA, InfoCellB, InfoCellC;

int InfoSSBondCount;
int InfoLadderCount;
int InfoChainCount;
int InfoHBondCount;
int InfoHelixCount;
int InfoTurnCount;
Long InfoBondCount;

int MainGroupCount,HetaGroupCount;
Long MainAtomCount; 
int HetaAtomCount;

Long MinX, MinY, MinZ;
Long MaxX, MaxY, MaxZ;

int MinMainTemp, MaxMainTemp;
int MinHetaTemp, MaxHetaTemp;
int MinMainRes,  MaxMainRes;
int MinHetaRes,  MaxHetaRes;

Molecule __far *Database;
MaskDesc UserMask[MAXMASK];
Long MinHBondDist, MaxHBondDist;
Long MinBondDist,  MaxBondDist;
int AbsMaxBondDist;
int ElemNo,ResNo;
int HasHydrogen;
int MaskCount;

#else
extern char InfoFileName[256];
extern char Residue[MAXRES][4];
extern char ElemDesc[MAXELEM][4];
extern char InfoClassification[42];
extern char InfoMoleculeName[80];
extern char InfoSpaceGroup[11];
extern char InfoIdentCode[6];

extern Real InfoCellAlpha, InfoCellBeta, InfoCellGamma;
extern Real InfoCellA, InfoCellB, InfoCellC;

extern int InfoSSBondCount;
extern int InfoLadderCount;
extern int InfoChainCount;
extern int InfoHBondCount;
extern int InfoHelixCount;
extern int InfoTurnCount;
extern Long InfoBondCount;

extern int MainGroupCount,HetaGroupCount;
extern Long MainAtomCount;
extern int HetaAtomCount;

extern Long MinX, MinY, MinZ;
extern Long MaxX, MaxY, MaxZ;

extern int MinMainTemp, MaxMainTemp;
extern int MinHetaTemp, MaxHetaTemp;
extern int MinMainRes,  MaxMainRes;
extern int MinHetaRes,  MaxHetaRes;

extern Molecule __far *Database;
extern MaskDesc UserMask[MAXMASK];
extern Long MinHBondDist, MaxHBondDist;
extern Long MinBondDist,  MaxBondDist;
extern int AbsMaxBondDist;
extern int ElemNo,ResNo;
extern int HasHydrogen;
extern int MaskCount;

#ifdef FUNCPROTO
int LoadAlchemyMolecule( FILE* );
int LoadCharmmMolecule( FILE* );
int LoadMol2Molecule( FILE* );
int LoadPDBMolecule( FILE* );
int LoadXYZMolecule( FILE* );
int LoadMDLMolecule( FILE* );

int SaveAlchemyMolecule( char* );
int SaveCIFMolecule( char* );
int SavePDBMolecule( char* );
int SaveXYZMolecule( char* );

void CreateMoleculeBonds( int, int );
void FindDisulphideBridges();
void CalcHydrogenBonds();
void DetermineStructure();
void RenumberMolecule( int );
void InitialiseDatabase();
void DescribeMolecule();
void DestroyDatabase();
void PurgeDatabase();

Atom __far *FindGroupAtom( Group __far*, int );

#else /* non-ANSI C compiler */
int LoadAlchemyMolecule();
int LoadCharmmMolecule();
int LoadMol2Molecule();
int LoadPDBMolecule();
int LoadXYZMolecule();
int LoadMDLMolecule();

int SaveAlchemyMolecule();
int SaveMol2Molecule();
int SaveCIFMolecule();
int SavePDBMolecule();
int SaveXYZMolecule();

void CreateMoleculeBonds();
void FindDisulphideBridges();
void CalcHydrogenBonds();
void DetermineStructure();
void RenumberMolecule();
void InitialiseDatabase();
void DescribeMolecule();
void DestroyDatabase();
void PurgeDatabase();

Atom __far *FindGroupAtom();

#endif
#endif

