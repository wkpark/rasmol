/* pixutils.h
 * RasMol2 Molecular Graphics
 * Roger Sayle, February 1994
 * Version 2.3
 */

typedef struct {
        int px, py, pz;
        int wx, wy, wz;
        int tx, ty, tz;
        int cx, cy, cz;
        int inten;
	} Knot;

#ifdef PIXUTILS
int SplineCount;

#else
extern int SplineCount;

#ifdef __STDC__
void PlotPoint( int, int, int, int );
void PlotDeepPoint( int, int, int, int );
void DrawTwinLine( int, int, int, int, int, int, int, int );
void ClipTwinLine( int, int, int, int, int, int, int, int );
void DrawTwinVector( int, int, int, int, int, int, int, int );
void ClipTwinVector( int, int, int, int, int, int, int, int );
void ClipDashVector( int, int, int, int, int, int, int, int );

void DrawCylinder( int, int, int, int, int, int, int, int, int );
void ClipCylinder( int, int, int, int, int, int, int, int, int );
void StrandRibbon( Knot __far*, Knot __far*, int );
void SolidRibbon( Knot __far*, Knot __far*, int );
void DrawSphere( int, int, int, int, int );
void ClipSphere( int, int, int, int, int );

void InitialisePixUtils();

#else /* non-ANSI C compiler */
void PlotPoint();
void PlotDeepPoint();
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
void InitialisePixUtils();

#endif
#endif

