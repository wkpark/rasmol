/* render.h
 * RasMol2 Molecular Graphics
 * Roger Sayle, February 1994
 * Version 2.3
 */

#define SlabReject       0x00
#define SlabHalf         0x01
#define SlabHollow       0x02
#define SlabFinal        0x03
#define SlabClose        0x04
#define SlabSection      0x05

#define ColBits          24

#define VOXORDER       21
#define VOXORDER2      (VOXORDER*VOXORDER)
#define VOXSIZE        (VOXORDER2*VOXORDER)


#ifdef RENDER
int UseShadow, DisplayMode;
int UseClipping,UseSlabPlane;
int SlabMode, SlabValue;
int SlabInten, SliceValue;
int ImageRadius,ImageSize;
int SSBondMode, HBondMode;
int RibbonMode;
int DrawBoundBox;
int DrawUnitCell;
int DrawAxes;

int VoxelsClean;
int BucketFlag;
int FBClear;

Pixel __huge *FBuffer;
short __huge *DBuffer;

#ifdef IBMPC
void __far * __far *HashTable;
Byte __far * __far *LookUp;
Card __far *ColConst;
Byte __far *Array;

HGLOBAL FBufHandle;
HGLOBAL DBufHandle;
#else /* UNIX */
void *HashTable[VOXSIZE];
char *LookUp[120];
int ColConst[120];
char Array[7261];
#endif

#else
extern int UseShadow, DisplayMode;
extern int UseClipping,UseSlabPlane;
extern int SlabMode,SlabValue;
extern int SlabInten,SliceValue;
extern int ImageRadius,ImageSize;
extern int SSBondMode, HBondMode;
extern int RibbonMode;
extern int DrawBoundBox;
extern int DrawUnitCell;
extern int DrawAxes;

extern int VoxelsClean;
extern int BucketFlag;
extern int FBClear;

extern Pixel __huge *FBuffer;
extern short __huge *DBuffer;

#ifdef IBMPC
extern void __far * __far *HashTable;
extern Byte __far * __far *LookUp;
extern Card __far *ColConst;
extern Byte __far *Array;

extern HGLOBAL FBufHandle;
extern HGLOBAL DBufHandle;
#else /* UNIX */
extern void *HashTable[VOXSIZE];
extern char *LookUp[120];
extern int ColConst[120];
extern char Array[7261];
#endif

#ifdef __STDC__
void InitialiseRenderer();
void ResetRenderer();

void ClearBuffers();
void ReSizeScreen();
void ReAllocBuffers();
void ResetVoxelData();
void CreateVoxelData();
void ShadowTransform();
void DrawFrame();
void IdentifyAtom( int, int );
int isqrt( Card );

#else /* non-ANSI C compiler */
void InitialiseRenderer();
void ResetRenderer();

void ClearBuffers();
void ReSizeScreen();
void ReAllocBuffers();
void ResetVoxelData();
void CreateVoxelData();
void ShadowTransform();
void DrawFrame();
void IdentifyAtom();
int isqrt();

#endif
#endif

