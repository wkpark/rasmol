/* pixutils.h
 * RasMol2 Molecular Graphics
 * Roger Sayle, October 1994
 * Version 2.5
 */

typedef struct {
        int px, py, pz;
        int wx, wy, wz;
        int tx, ty, tz;
        int cx, cy, cz;
        int inten;
	} Knot;

#define ZValid(z)     ((!UseSlabPlane) || ((z)<SlabValue))
#define XValid(x)     (((x)>=0)&&((x)<XRange))
#define YValid(y)     (((y)>=0)&&((y)<YRange))


#ifdef PIXUTILS
int SplineCount;
int FontSize;

#else
extern int SplineCount;
extern int FontSize;

#ifdef FUNCPROTO
void PlotDeepPoint( int, int, int, int );
void ClipDeepPoint( int, int, int, int );
void DrawTwinLine( int, int, int, int, int, int, int, int );
void ClipTwinLine( int, int, int, int, int, int, int, int );
void DrawTwinVector( int, int, int, int, int, int, int, int );
void ClipTwinVector( int, int, int, int, int, int, int, int );
void ClipDashVector( int, int, int, int, int, int, int, int );

void DrawCylinder( int, int, int, int, int, int, int, int, int );
void ClipCylinder( int, int, int, int, int, int, int, int, int );
void StrandRibbon( Knot __far*, Knot __far*, int, int );
void SolidRibbon( Knot __far*, Knot __far*, int );
void DrawSphere( int, int, int, int, int );
void ClipSphere( int, int, int, int, int );

void SetFontSize( int );
void DisplayString( int, int, int, char*, int );
void InitialisePixUtils();

#else /* non-ANSI C compiler */
void PlotDeepPoint();
void ClipDeepPoint();
void DrawTwinLine();
void ClipTwinLine();
void DrawTwinVector();
void ClipTwinVector();
void ClipDashVector();

void DrawCylinder();
void ClipCylinder();
void StrandRibbon();
void SolidRibbon();
void DrawSphere();
void ClipSphere();

void SetFontSize();
void DisplayString();
void InitialisePixUtils();

#endif
#endif

