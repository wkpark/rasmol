/* transfor.h
 * RasMol2 Molecular Graphics
 * Roger Sayle, October 1994
 * Version 2.5
 */

#define GroupAttr       0x00
#define ChainAttr       0x01
#define TempAttr	0x02
#define ChargeAttr      0x03

#define MaskColourFlag  0x01
#define MaskRadiusFlag  0x02
#define MaskBothFlag    0x03

#define RibColInside    0x01
#define RibColOutside   0x02
#define RibColBoth      0x03

#ifdef EIGHTBIT
#define DefaultAmbient    0.6
#define DefaultColDepth   16
#else
#define DefaultAmbient    0.05
#define DefaultColDepth   32
#endif


#define MAXSHADE 32
typedef struct { 
        Long refcount;
        unsigned char r;
        unsigned char g;
        unsigned char b;
    } ShadeDesc;


#ifdef IBMPC
#define Colour2Shade(x)  ((int)((x)-3)/ColourDepth)
#define Shade2Colour(x)  ((x)*ColourDepth+3)
#define BackCol    0
#define BoxCol     1
#define LabelCol   2
#endif

#ifdef APPLEMAC
#define Colour2Shade(x)  ((int)((x)-4)/ColourDepth)
#define Shade2Colour(x)  ((x)*ColourDepth+4)
#define BackCol    1
#define BoxCol     2
#define LabelCol   3
#endif

#if !defined(IBMPC) && !defined(APPLEMAC)
#define Colour2Shade(x)  ((int)((x)-8)/ColourDepth)
#define Shade2Colour(x)  ((x)*ColourDepth+8)
#define BackCol    5
#define BoxCol     6
#define LabelCol   7
#endif


#ifdef TRANSFORM
ShadeDesc Shade[MAXSHADE];
Real RotX[3],RotY[3],RotZ[3];
Real MatX[3],MatY[3],MatZ[3];
Real InvX[3],InvY[3],InvZ[3];
Long OrigCX,OrigCY,OrigCZ;
Long CenX, CenY, CenZ;

Real Ambient;
Real Scale,MaxZoom;
Real DScale,IScale;
Long SideLen,Offset;
Card WorldRadius,WorldSize;
int XOffset,YOffset,ZOffset;
int FakeSpecular,SpecPower;
int ColourDepth,ColourMask;
int BackR,BackG,BackB;
int LabR,LabG,LabB;
int BoxR,BoxG,BoxB;
int UseLabelCol;
int UseScreenClip;
int ZoomRange;

int Hydrogens,HetaGroups;
int DrawAtoms,MaxAtomRadius;
int DrawBonds,MaxBondRadius;
int DrawRibbon;
int ZoneBoth;

#else
extern ShadeDesc Shade[MAXSHADE];
extern Real RotX[3],RotY[3],RotZ[3];
extern Real MatX[3],MatY[3],MatZ[3];
extern Real InvX[3],InvY[3],InvZ[3];
extern Long OrigCX, OrigCY, OrigCZ;
extern Long CenX, CenY, CenZ;

extern Real Ambient;
extern Real Scale,MaxZoom;
extern Real DScale,IScale;
extern Long SideLen,Offset;
extern Card WorldRadius,WorldSize;
extern int XOffset,YOffset,ZOffset;
extern int ColourDepth,ColourMask;
extern int FakeSpecular,SpecPower;
extern int BackR,BackG,BackB;
extern int LabR,LabG,LabB;
extern int BoxR,BoxG,BoxB;
extern int UseLabelCol;
extern int UseScreenClip;
extern int ZoomRange;

extern int Hydrogens,HetaGroups;
extern int DrawAtoms,MaxAtomRadius;
extern int DrawBonds,MaxBondRadius;
extern int DrawRibbon;
extern int ZoneBoth;

#ifdef FUNCPROTO
void SetRadiusValue( int );
void SetRadiusTemperature();
void SetVanWaalRadius();
void DisableSpacefill();
void EnableWireFrame( int, int );
void EnableBackBone( int, int );
void SetHBondStatus( int, int, int );
void SetRibbonStatus( int, int, int );
void DisableWireFrame();
void DisableBackBone();
void RestrictZone( int );
void SelectZone( int );

int IsCPKColour( Atom __far * );
int IsVDWRadius( Atom __far * );

void DefineColourMap();
void ResetColourMap();
void ColourBackNone();
void ColourBondNone();
void ColourHBondNone( int );
void ColourHBondType();
void ColourRibbonNone( int );
void ColourBackAttrib( int, int, int );
void ColourBondAttrib( int, int, int );
void ColourHBondAttrib( int, int, int, int );
void ColourRibbonAttrib( int, int, int, int );
void ColourDotsAttrib( int, int, int );
void ColourDotsPotential();
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

int IsCPKColour();
int IsVDWRadius();

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
void ColourDotsAttrib();
void ColourDotsPotential();
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

