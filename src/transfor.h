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

/* transfor.h
 $Log: transfor.h,v $
 Revision 1.1  2001/01/31 02:13:45  yaya
 Initial revision

 Revision 1.7  2000/08/26 18:13:02  yaya
 Updates to header comments in all files

 Revision 1.6  2000/08/26 17:31:24  yaya
 Fix for world rot, remove refs to toolbar

 Revision 1.5  2000/08/26 03:14:29  yaya
 Mods for mac compilations

 Revision 1.4  2000/08/21 21:08:09  yaya
 semi-final ucb mods

 Revision 1.3  2000/08/09 01:18:42  yaya
 Rough cut with ucb

 Revision 1.2  2000/08/03 18:32:43  yaya
 Parametrization for alt conformer bond radius

 */

#define GroupAttr       0x00
#define ChainAttr       0x01
#define TempAttr	    0x02
#define ChargeAttr      0x03
#define AltAttr         0x04
#define ModelAttr       0x05

#define MaskColourFlag  0x01
#define MaskRadiusFlag  0x02
#define MaskBothFlag    0x03

#define RibColInside    0x01
#define RibColOutside   0x02
#define RibColBoth      0x03

#ifdef EIGHTBIT
#define ColourDepth       16
#define ColourMask        15
#define AltlDepth         4
#ifdef APPLEMAC
#define LastShade         14
#else
#define LastShade         15
#endif
#else
#define ColourDepth       32
#define ColourMask        31
#define LastShade         31
#define AltlDepth         8
#endif


#ifdef __esv
/* Evans & Sutherland Gamma! */
#define DefaultAmbient    0.05
#else
#define DefaultAmbient    0.4
#endif


typedef struct { 
        Long refcount;
        unsigned char r;
        unsigned char g;
        unsigned char b;
    } ShadeDesc;

typedef struct {
        short col;
        short shade;
        unsigned char r;
        unsigned char g;
        unsigned char b;
    } ShadeRef;


#ifdef IBMPC
#define BackCol    0
#define BoxCol     1
#define LabelCol   2
#define FirstCol   3
#endif

#ifdef APPLEMAC
#define BackCol    1
#define BoxCol     2
#define LabelCol   3
#define FirstCol   4
#endif

#if !defined(IBMPC) && !defined(APPLEMAC)
#define BackCol    5
#define BoxCol     6
#define LabelCol   7
#define FirstCol   8
#endif

#define Colour2Shade(x)  ((int)((x)-FirstCol)/ColourDepth)
#define Shade2Colour(x)  ((x)*ColourDepth+FirstCol)


#ifdef TRANSFORM
ShadeDesc Shade[LastShade];
Real RotX[3],RotY[3],RotZ[3];
Real LRotX[3],LRotY[3],LRotZ[3];
Real LOffset[3];
Real MatX[3],MatY[3],MatZ[3];
Real InvX[3],InvY[3],InvZ[3];
Long OrigCX,OrigCY,OrigCZ;
Long CenX, CenY, CenZ;
Long ShiftS;
int XlateCen;

int FakeSpecular,SpecPower;
int ShadePower;
int BackR,BackG,BackB;
int DotR,DotG,DotB;
int LabR,LabG,LabB;
int BoxR,BoxG,BoxB;
int UseLabelCol;
int UseBackFade;
Real Ambient;
int UseDotColPot;

Real Scale,MaxZoom;
Real DScale,IScale;
Long SideLen,Offset;
Card WorldRadius,WorldSize,LocalRadius;
int XOffset,YOffset,ZOffset;
int UseScreenClip;
int ZoomRange;

int Hydrogens,HetaGroups;
int DrawAtoms,MaxAtomRadius;
int DrawBonds,MaxBondRadius;
int DrawStars;
int DrawRibbon;
int ZoneBoth;
int ModelInclude;

int ScaleCount;
ShadeRef ScaleRef[LastShade];
int AltlColours[AltlDepth];

Real LastRX,LastRY,LastRZ;
Real LastTX,LastTY,LastTZ;
Real Zoom;


#else
extern ShadeDesc Shade[LastShade];
extern Real RotX[3],RotY[3],RotZ[3];
extern Real LRotX[3],LRotY[3],LRotZ[3];
extern Real LOffset[3];
extern Real MatX[3],MatY[3],MatZ[3];
extern Real InvX[3],InvY[3],InvZ[3];
extern Long OrigCX, OrigCY, OrigCZ;
extern Long CenX, CenY, CenZ;
extern Long ShiftS;
extern int XlateCen;


extern int FakeSpecular,SpecPower;
extern int ShadePower;
extern int BackR,BackG,BackB;
extern int DotR,DotG,DotB;
extern int LabR,LabG,LabB;
extern int BoxR,BoxG,BoxB;
extern int UseLabelCol;
extern int UseBackFade;
extern Real Ambient;
extern int UseDotColPot;

extern Real Scale,MaxZoom;
extern Real DScale,IScale;
extern Long SideLen,Offset;
extern Card WorldRadius,WorldSize,LocalRadius;
extern int XOffset,YOffset,ZOffset;
extern int UseScreenClip;
extern int ZoomRange;

extern int Hydrogens,HetaGroups;
extern int DrawAtoms,MaxAtomRadius;
extern int DrawBonds,MaxBondRadius;
extern int DrawStars;
extern int DrawRibbon;
extern int ZoneBoth;
extern int ModelInclude;

extern int ScaleCount;
extern ShadeRef ScaleRef[LastShade];
extern int AltlColours[AltlDepth];

extern ShadeDesc Shade[LastShade];
extern Real RotX[3],RotY[3],RotZ[3];

extern Real LastRX,LastRY,LastRZ;
extern Real LastTX,LastTY,LastTZ;
extern Real Zoom;

#endif


void SetRadiusValue( int, int  );
void SetRadiusTemperature( int );
void SetVanWaalRadius( int );
void DisableSpacefill( void );
void SetHBondStatus( int, int, int, int );
void SetRibbonStatus( int, int, int );
void SetRibbonCartoons( void );
void SetTraceTemperature( void );
void EnableWireframe( int, int, int );
void EnableBackbone( int, int, int );
void DisableWireframe( void );
void DisableBackbone( void );

void DisplaySelectCount( void );
void SelectZoneExpr( Expr* );
void RestrictZoneExpr( Expr* );
void RestrictZone( int );
void SelectArea( int, int, int, int, int, int );
void SelectZone( int );
void SelectAtom( int, RAtom __far *, Group __far * );
void SelectGroup( int, Group __far * );
void SelectChain( int, Chain __far * );

int IsCPKColour( RAtom __far * );
int IsVDWRadius( RAtom __far * );

void DefineColourMap( void );
void ResetColourMap( void );

void ColourBackNone( void );
void ColourBondNone( void );
void ColourHBondType( void );
void ColourHBondNone( int );
void ColourRibbonNone( int );
void ColourMonitNone( void );
void ColourBackAttrib( int, int, int );
void ColourBondAttrib( int, int, int );
void ColourHBondAttrib( int, int, int, int );
void ColourRibbonAttrib( int, int, int, int );
void ColourMonitAttrib( int, int, int );
void ColourDotsAttrib( int, int, int );
void ColourDotsPotential( void );
void MonoColourAttrib( int, int, int );
void ScaleColourAttrib( int );
void CPKColourAttrib( void );
void AminoColourAttrib( void );
void ShapelyColourAttrib( void );
void StructColourAttrib( void );
void UserMaskAttrib( int );

void DefaultRepresentation( void );

void DetermineClipping( void );
void InitialiseTransform( void );
void InitialTransform( void );
void RMatVec( Real [3], Real [3], Real [3], Real [3], Real [3] );
void RMatRMat( Real [3][3], Real [3], Real [3], Real [3], Real [3][3] );
void RMatInv( Real [3], Real [3], Real [3], Real [3], Real [3], Real [3] );
void PrepareTransform( void );
void ReviseInvMatrix( void );
void ApplyTransform( void );
void CentreTransform( int, int, int, int );
void ResetTransform( void );

void SetLutEntry( int, int, int, int );
