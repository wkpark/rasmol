/***************************************************************************
 *                             RasMol 2.7.2.1                              *
 *                                                                         *
 *                                 RasMol                                  *
 *                 Molecular Graphics Visualisation Tool                   *
 *                              14 April 2001                              *
 *                                                                         *
 *                   Based on RasMol 2.6 by Roger Sayle                    *
 * Biomolecular Structures Group, Glaxo Wellcome Research & Development,   *
 *                      Stevenage, Hertfordshire, UK                       *
 *         Version 2.6, August 1995, Version 2.6.4, December 1998          *
 *                   Copyright (C) Roger Sayle 1992-1999                   *
 *                                                                         *
 *                          and Based on Mods by                           *
 *Author             Version, Date             Copyright                   *
 *Arne Mueller       RasMol 2.6x1   May 98     (C) Arne Mueller 1998       *
 *Gary Grossman and  RasMol 2.5-ucb Nov 95     (C) UC Regents/ModularCHEM  *
 *Marco Molinaro     RasMol 2.6-ucb Nov 96         Consortium 1995, 1996   *
 *                                                                         *
 *Philippe Valadon   RasTop 1.3     Aug 00     (C) Philippe Valadon 2000   *
 *                                                                         *
 *Herbert J.         RasMol 2.7.0   Mar 99     (C) Herbert J. Bernstein    * 
 *Bernstein          RasMol 2.7.1   Jun 99         1998-2001               *
 *                   RasMol 2.7.1.1 Jan 01                                 *
 *                   RasMol 2.7.2   Aug 00                                 *
 *                   RasMol 2.7.2.1 Apr 01                                 *
 *                                                                         *
 *                    and Incorporating Translations by                    *
 *  Author                               Item                      Language*
 *  Isabel Serván Martínez,                                                *
 *  José Miguel Fernández Fernández      2.6   Manual              Spanish *
 *  José Miguel Fernández Fernández      2.7.1 Manual              Spanish *
 *  Fernando Gabriel Ranea               2.7.1 menus and messages  Spanish *
 *  Jean-Pierre Demailly                 2.7.1 menus and messages  French  *
 *  Giuseppe Martini, Giovanni Paolella, 2.7.1 menus and messages          *
 *  A. Davassi, M. Masullo, C. Liotto    2.7.1 help file           Italian *
 *                                                                         *
 *                             This Release by                             *
 * Herbert J. Bernstein, Bernstein + Sons, P.O. Box 177, Bellport, NY, USA *
 *                       yaya@bernstein-plus-sons.com                      *
 *               Copyright(C) Herbert J. Bernstein 1998-2001               *
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
 $Log: render.h,v $
 Revision 1.1  2001/01/31 02:13:45  yaya
 Initial revision

 Revision 1.4  2000/08/26 18:12:59  yaya
 Updates to header comments in all files

 Revision 1.3  2000/08/13 20:56:43  yaya
 Conversion from toolbar to menus

 Revision 1.2  2000/08/09 01:18:37  yaya
 Rough cut with ucb

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
#define PickAtom         0x0A
#define PickGroup        0x0B
#define PickChain        0x0C
#define PickBond         0x0D

#define RotBond          0x01
#define RotMol           0x02
#define RotAll           0x04

#define ViewLeft         0
#define ViewRight        1

#define ColBits          24

#define VOXORDER       21
#define VOXORDER2      (VOXORDER*VOXORDER)
#define VOXSIZE        (VOXORDER2*VOXORDER)

typedef struct _Item {
        struct _Item __far *list;
        RAtom  __far *data;
    } Item;
 


#ifdef RENDER
int UseDepthCue;
int UseStereo,StereoView;
int UseShadow,DisplayMode;
int UseClipping,UseSlabPlane;
int UseAutoDepthCue,UseDepthPlane;
int SlabMode,SlabValue,DepthValue;
int SlabInten,SliceValue;
int ImageRadius,ImageSize;
int SSBondMode,HBondMode;
int PickCount;
int LabelOptFlag;

double StereoAngle;
int PickMode;
int RotMode;
int DrawArea;
int AreaX1, AreaX2, AreaY1, AreaY2;

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
extern int UseAutoDepthCue,UseDepthPlane;
extern int SlabMode,SlabValue,DepthValue;
extern int SlabInten,SliceValue;
extern int ImageRadius,ImageSize;
extern int SSBondMode, HBondMode;
extern int PickCount;
extern int LabelOptFlag;

extern double StereoAngle;
extern int PickMode;
extern int RotMode;
extern int DrawArea;
extern int AreaX1, AreaX2, AreaY1, AreaY2;

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
int PickAtoms( int, int, int );
unsigned int isqrt( Card );

