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

/* render.h
*/

/* These values set the sizes of the sphere rendering
 * tables. The first value, maxrad, is the maximum
 * sphere radius and the second value is the table
 * size = (maxrad*(maxrad+1))/2 + 1
 */
/* #define MAXRAD    120   256   */
/* #define MAXTABLE  7261  32897 */
#define MAXRAD    255
#define MAXTABLE  32641


#define SlabReject       0x00
#define SlabHalf         0x01
#define SlabHollow       0x02
#define SlabFinal        0x03
#define SlabClose        0x04
#define SlabSection      0x05

#define PickNone         0x00
#define PickIdent        0x01
#define PickDist         0x02
#define PickAngle        0x03
#define PickTorsn        0x04
#define PickLabel        0x05
#define PickMonit        0x06
#define PickCentr        0x07
#define PickOrign        0x08
#define PickCoord        0x09


#define ViewLeft         0
#define ViewRight        1

#define ColBits          24

#define VOXORDER       21
#define VOXORDER2      (VOXORDER*VOXORDER)
#define VOXSIZE        (VOXORDER2*VOXORDER)

typedef struct _Item {
        struct _Item __far *list;
        Atom  __far *data;
    } Item;
 


#ifdef RENDER
int UseDepthCue;
int UseStereo,StereoView;
int UseShadow,DisplayMode;
int UseClipping,UseSlabPlane;
int SlabMode,SlabValue;
int SlabInten,SliceValue;
int ImageRadius,ImageSize;
int SSBondMode,HBondMode;

double StereoAngle;
int PickMode;

int DrawBoundBox,DrawAxes;
int DrawDoubleBonds;
int DrawUnitCell;

Real IVoxRatio;
int VoxelsClean;
int BucketFlag;
int FBClear;


Card __far *ColConst;
#if defined(IBMPC) || defined(APPLEMAC)
void __far * __far *HashTable;
Byte __far * __far *LookUp;
Byte __far *Array;

#else /* UNIX or VMS */
void *HashTable[VOXSIZE];
Byte *LookUp[MAXRAD];
Byte Array[MAXTABLE];
#endif

#else
extern int UseDepthCue;
extern int UseStereo,StereoView;
extern int UseShadow, DisplayMode;
extern int UseClipping,UseSlabPlane;
extern int SlabMode,SlabValue;
extern int SlabInten,SliceValue;
extern int ImageRadius,ImageSize;
extern int SSBondMode, HBondMode;

extern double StereoAngle;
extern int PickMode;

extern int DrawBoundBox,DrawAxes;
extern int DrawDoubleBonds;
extern int DrawUnitCell;

extern Real IVoxRatio;
extern int VoxelsClean;
extern int BucketFlag;
extern int FBClear;


extern Card __far *ColConst;
#if defined(IBMPC) || defined(APPLEMAC)
extern void __far * __far *HashTable;
extern Byte __far * __far *LookUp;
extern Byte __far *Array;

#else /* UNIX or VMS */
extern void *HashTable[VOXSIZE];
extern Byte *LookUp[MAXRAD];
extern Byte Array[MAXTABLE];
#endif
#endif


void ClearBuffers( void );
void ReSizeScreen( void );
void ReAllocBuffers( void );
void ShadowTransform( void );

void ResetVoxelData( void );
void CreateVoxelData( int );

void DrawFrame( void );
void ResetRenderer( void );
void InitialiseRenderer( void );
void SetStereoMode( int );
void SetPickMode( int );
int PickAtom( int, int, int );
unsigned int isqrt( Card );

