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

/* transfor.h
 */

#define GroupAttr       0x00
#define ChainAttr       0x01
#define TempAttr	0x02
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
Real MatX[3],MatY[3],MatZ[3];
Real InvX[3],InvY[3],InvZ[3];
Long OrigCX,OrigCY,OrigCZ;
Long CenX, CenY, CenZ;

int FakeSpecular,SpecPower;
int BackR,BackG,BackB;
int LabR,LabG,LabB;
int BoxR,BoxG,BoxB;
int UseLabelCol;
int UseBackFade;
Real Ambient;

Real Scale,MaxZoom;
Real DScale,IScale;
Long SideLen,Offset;
Card WorldRadius,WorldSize;
int XOffset,YOffset,ZOffset;
int UseScreenClip;
int ZoomRange;

int Hydrogens,HetaGroups;
int DrawAtoms,MaxAtomRadius;
int DrawBonds,MaxBondRadius;
int DrawStars;
int DrawRibbon;
int ZoneBoth;

int ScaleCount;
ShadeRef ScaleRef[LastShade];
int AltlColours[AltlDepth];

#else
extern ShadeDesc Shade[LastShade];
extern Real RotX[3],RotY[3],RotZ[3];
extern Real MatX[3],MatY[3],MatZ[3];
extern Real InvX[3],InvY[3],InvZ[3];
extern Long OrigCX, OrigCY, OrigCZ;
extern Long CenX, CenY, CenZ;

extern int FakeSpecular,SpecPower;
extern int BackR,BackG,BackB;
extern int LabR,LabG,LabB;
extern int BoxR,BoxG,BoxB;
extern int UseLabelCol;
extern int UseBackFade;
extern Real Ambient;

extern Real Scale,MaxZoom;
extern Real DScale,IScale;
extern Long SideLen,Offset;
extern Card WorldRadius,WorldSize;
extern int XOffset,YOffset,ZOffset;
extern int UseScreenClip;
extern int ZoomRange;

extern int Hydrogens,HetaGroups;
extern int DrawAtoms,MaxAtomRadius;
extern int DrawBonds,MaxBondRadius;
extern int DrawStars;
extern int DrawRibbon;
extern int ZoneBoth;

extern int ScaleCount;
extern ShadeRef ScaleRef[LastShade];
extern int AltlColours[AltlDepth];

extern ShadeDesc Shade[LastShade];
extern Real RotX[3],RotY[3],RotZ[3];
extern Real MatX[3],MatY[3],MatZ[3];
extern Real InvX[3],InvY[3],InvZ[3];
extern Long OrigCX, OrigCY, OrigCZ;
extern Long CenX, CenY, CenZ;
#endif


void SetRadiusValue( int, int  );
void SetRadiusTemperature( int );
void SetVanWaalRadius( int );
void DisableSpacefill();
void SetHBondStatus( int, int, int );
void SetRibbonStatus( int, int, int );
void SetRibbonCartoons();
void SetTraceTemperature();
void EnableWireframe( int, int );
void EnableBackbone( int, int );
void DisableWireframe();
void DisableBackbone();

void SelectZoneExpr( Expr* );
void RestrictZoneExpr( Expr* );
void RestrictZone( int );
void SelectZone( int );

int IsCPKColour( Atom __far * );
int IsVDWRadius( Atom __far * );

void DefineColourMap();
void ResetColourMap();

void ColourBackNone();
void ColourBondNone();
void ColourHBondType();
void ColourHBondNone( int );
void ColourRibbonNone( int );
void ColourMonitNone();
void ColourBackAttrib( int, int, int );
void ColourBondAttrib( int, int, int );
void ColourHBondAttrib( int, int, int, int );
void ColourRibbonAttrib( int, int, int, int );
void ColourMonitAttrib( int, int, int );
void ColourDotsAttrib( int, int, int );
void ColourDotsPotential();
void MonoColourAttrib( int, int, int );
void ScaleColourAttrib( int );
void CPKColourAttrib();
void AminoColourAttrib();
void ShapelyColourAttrib();
void StructColourAttrib();
void UserMaskAttrib( int );

void DefaultRepresentation();

void DetermineClipping();
void InitialiseTransform();
void InitialTransform();
void PrepareTransform();
void ReviseInvMatrix();
void ApplyTransform();
void ResetTransform();
