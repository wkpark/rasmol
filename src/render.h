/* render.h
 * RasMol2 Molecular Graphics
 * Roger Sayle, October 1994
 * Version 2.5
 */

#define SlabReject       0x00
#define SlabHalf         0x01
#define SlabHollow       0x02
#define SlabFinal        0x03
#define SlabClose        0x04
#define SlabSection      0x05

#define ViewLeft         0
#define ViewRight        1

#define ColBits          24

#define VOXORDER       21
#define VOXORDER2      (VOXORDER*VOXORDER)
#define VOXSIZE        (VOXORDER2*VOXORDER)


#define DotMax    100
typedef struct _DotStruct {
        struct _DotStruct __far *next;
        short col[DotMax];
        Long xpos[DotMax];
        Long ypos[DotMax];
        Long zpos[DotMax];
        int count;
    } DotStruct;


typedef struct _Label {
        struct _Label *next;
        Long  refcount;
        char *label;
    } Label;



#ifdef RENDER
DotStruct __far *DotPtr;
int UseStereo,StereoView;
int UseShadow,DisplayMode;
int UseClipping,UseSlabPlane;
int SolventDots,ProbeRadius;
int SlabMode,SlabValue;
int SlabInten,SliceValue;
int ImageRadius,ImageSize;
int SSBondMode,HBondMode;

int DrawDots,DrawLabels;
int DrawBoundBox,DrawAxes;
int DrawUnitCell;

Pixel __huge *FBuffer;
short __huge *DBuffer;
int VoxelsClean;
int BucketFlag;
int FBClear;


#ifdef IBMPC
HGLOBAL FBufHandle;
HGLOBAL DBufHandle;
#endif

#ifdef APPLEMAC
Handle FBufHandle;
Handle DBufHandle;
#endif

#if defined(IBMPC) || defined(APPLEMAC)
void __far * __far *HashTable;
Byte __far * __far *LookUp;
Card __far *ColConst;
Byte __far *Array;

#else /* UNIX */
void *HashTable[VOXSIZE];
Byte *LookUp[120];
Card ColConst[120];
Byte Array[7261];
#endif

#else
extern DotStruct __far *DotPtr;
extern int UseShadow, DisplayMode;
extern int UseClipping,UseSlabPlane;
extern int ProbeRadius,SolventDots;
extern int SlabMode,SlabValue;
extern int SlabInten,SliceValue;
extern int ImageRadius,ImageSize;
extern int SSBondMode, HBondMode;

extern int DrawDots,DrawLabels;
extern int DrawBoundBox,DrawAxes;
extern int DrawUnitCell;

extern Pixel __huge *FBuffer;
extern short __huge *DBuffer;
extern int VoxelsClean;
extern int BucketFlag;
extern int FBClear;

#ifdef IBMPC
extern HGLOBAL FBufHandle;
extern HGLOBAL DBufHandle;
#endif

#ifdef APPLEMAC
extern Handle FBufHandle;
extern Handle DBufHandle;
#endif

#if defined(IBMPC) || defined(APPLEMAC)
extern void __far * __far *HashTable;
extern Byte __far * __far *LookUp;
extern Card __far *ColConst;
extern Byte __far *Array;

#else /* UNIX or VMS */
extern void *HashTable[VOXSIZE];
extern Byte *LookUp[120];
extern Card ColConst[120];
extern Byte Array[7261];
#endif

#ifdef FUNCPROTO
void ClearBuffers();
void ReSizeScreen();
void ReAllocBuffers();
void ResetVoxelData();
void ShadowTransform();

int DeleteLabels();
void DefineLabels( char* );
void CalculateSurface( int );
void DeleteSurface();

void DrawFrame();
void ResetRenderer();
void InitialiseRenderer();
void SetStereoMode( int );
void IdentifyAtom( int, int );
int isqrt( Card );

#else /* non-ANSI C compiler */
void ClearBuffers();
void ReSizeScreen();
void ReAllocBuffers();
void ResetVoxelData();
void ShadowTransform();

void DeleteLabels();
void DefineLabels();
void CalculateSurface();
void DeleteSurface();

void DrawFrame();
void ResetRenderer();
void InitialiseRenderer();
void SetStereoMode();
void IdentifyAtom();
int isqrt();

#endif
#endif

