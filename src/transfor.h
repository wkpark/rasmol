/* transfor.h
 * RasMol2 Molecular Graphics
 * Roger Sayle, February 1994
 * Version 2.3
 */

#define GroupAttr       0x00
#define ChainAttr       0x01
#define TempAttr	0x02

#define MaskColourFlag  0x01
#define MaskRadiusFlag  0x02
#define MaskBothFlag    0x03


#ifdef EIGHTBIT
#define DefaultAmbient    0.6
#define DefaultColDepth   16
#else
#define DefaultAmbient    0.05
#define DefaultColDepth   32
#endif


#ifdef IBMPC
#define Colour2Shade(x)  ((int)((x)-1)/ColourDepth)
#define Shade2Colour(x)  ((x)*ColourDepth+1)
#else
#define Colour2Shade(x)  ((int)((x)-6)/ColourDepth)
#define Shade2Colour(x)  ((x)*ColourDepth+6)
#endif


#ifdef TRANSFORM
Real RotX[3],RotY[3],RotZ[3];
Real MatX[3],MatY[3],MatZ[3];
Real InvX[3],InvY[3],InvZ[3];
Real DirX[3],DirY[3],DirZ[3];
Long CenX, CenY, CenZ;

Real Ambient;
Real Scale, MaxZoom;
Real DScale, IScale;
Long SideLen, Offset;
Card WorldRadius, WorldSize;
int FakeSpecular, SpecPower;
int ColourDepth, ColourMask;
int BackR, BackG, BackB;
int XOffset, YOffset;
int UseScreenClip;
int ZoomRange;

int Hydrogens, HetaGroups;
int DrawAtoms, MaxAtomRadius;
int DrawBonds, MaxBondRadius;
int DrawRibbon;
int ZoneBoth;

#else
extern Real RotX[3],RotY[3],RotZ[3];
extern Real MatX[3],MatY[3],MatZ[3];
extern Real InvX[3],InvY[3],InvZ[3];
extern Real DirX[3],DirY[3],DirZ[3];
extern Long CenX, CenY, CenZ;

extern Real Ambient;
extern Real Scale,MaxZoom;
extern Real DScale,IScale;
extern Long SideLen, Offset;
extern Card WorldRadius,WorldSize;
extern int ColourDepth,ColourMask;
extern int FakeSpecular,SpecPower;
extern int BackR,BackG,BackB;
extern int XOffset,YOffset;
extern int UseScreenClip;
extern int ZoomRange;

extern int Hydrogens,HetaGroups;
extern int DrawAtoms,MaxAtomRadius;
extern int DrawBonds,MaxBondRadius;
extern int DrawRibbon;
extern int ZoneBoth;

#ifdef __STDC__
void SetRadiusValue( int );
void SetRadiusTemperature();
void SetVanWaalRadius();
void DisableSpacefill();
void EnableWireFrame( int, int );
void EnableBackBone( int, int );
void SetHBondStatus( int, int, int );
void SetRibbonStatus( int, int );
void DisableWireFrame();
void DisableBackBone();
void RestrictZone( int );
void SelectZone( int );


void DefineColourMap();
void ResetColourMap();
void ColourBackNone();
void ColourBondNone();
void ColourHBondNone( int );
void ColourHBondType();
void ColourRibbonNone();
void ColourBackAttrib( int, int, int );
void ColourBondAttrib( int, int, int );
void ColourHBondAttrib( int, int, int, int );
void ColourRibbonAttrib( int, int, int );
void MonoColourAttrib( int, int, int );
void ScaleColourAttrib( int );
void CPKColourAttrib();
void AminoColourAttrib();
void ShapelyColourAttrib();
void StructColourAttrib();
void UserMaskAttrib( int );

void DetermineClipping();
void InitialiseTransform();
void InitialTransform();
void PrepareTransform();
void ReviseInvMatrix();
void ApplyTransform();
void ResetTransform();

#else /* non-ANSI C compiler */
void SetRadiusValue();
void SetRadiusTemperature();
void SetVanWaalRadius();
void DisableSpacefill();
void EnableWireFrame();
void EnableBackBone();
void SetHBondStatus();
void SetRibbonStatus();
void DisableWireFrame();
void DisableBackBone();
void RestrictZone();
void SelectZone();

void DefineColourMap();
void ResetColourMap();
void ColourBackNone();
void ColourBondNone();
void ColourHBondNone();
void ColourHBondType();
void ColourRibbonNone();
void ColourBackAttrib();
void ColourBondAttrib();
void ColourHBondAttrib();
void ColourRibbonAttrib();
void MonoColourAttrib();
void ScaleColourAttrib();
void CPKColourAttrib();
void AminoColourAttrib();
void ShapelyColourAttrib();
void StructColourAttrib();
void UserMaskAttrib();

void DetermineClipping();
void InitialiseTransform();
void InitialTransform();
void PrepareTransform();
void ReviseInvMatrix();
void ApplyTransform();
void ResetTransform();

#endif
#endif

