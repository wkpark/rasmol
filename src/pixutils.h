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

/* pixutils.h
 $Log: pixutils.h,v $
 Revision 1.1  2001/01/31 02:13:45  yaya
 Initial revision

 Revision 1.4  2000/08/26 18:12:58  yaya
 Updates to header comments in all files

 Revision 1.3  2000/08/13 20:56:42  yaya
 Conversion from toolbar to menus

 Revision 1.2  2000/08/03 18:32:43  yaya
 Parametrization for alt conformer bond radius

 */

typedef struct {
        int px, py, pz;      /* Spline Control Co-ordinate */
        int tx, ty, tz;      /* Spline Direction Vector    */
        int hnx, hny, hnz;   /* Horizontal Normal Vector   */
        int vnx, vny, vnz;   /* Vertical Normal Vector     */
        int dx, dy, dz;      /* Ribbon Height Vector       */
        int wx, wy, wz;      /* Ribbon Width Vector        */
        char hinten;         /* Horizontal Intensity       */
        char vinten;         /* Vertical Intensity         */
        short hsize;         /* Horizontal Vector Length   */
        short vsize;         /* Vertical Vector Length     */
        short wide;          /* Ribbon Width               */
    } Knot;

typedef struct {
        Pixel __huge *fbuf;
        short __huge *dbuf;
        int xmax, ymax;
        int yskip;
    } ViewStruct;


#define MAXVERT  6
typedef struct {
        int x, y, z;
        int inten;
    } Vert;

typedef struct {
        Vert v[MAXVERT];
        int count;
    } Poly;

#define ZValid(z)     ((!UseSlabPlane) || ((z)<SlabValue))
#define ZBack(z)      ((!UseDepthPlane) || ((z)>DepthValue))
#define XValid(x)     (((x)>=0)&&((x)<View.xmax))
#define YValid(y)     (((y)>=0)&&((y)<View.ymax))


#ifdef PIXUTILS
ViewStruct View;
int SplineCount;
int FontSize;
int FontPS;
int FontStroke;

#else
extern ViewStruct View;
extern int SplineCount;
extern int FontSize;
extern int FontPS;
extern int FontStroke;
#endif

void PlotDeepPoint( int, int, int, int );
void ClipDeepPoint( int, int, int, int );

void DrawCylinder( int, int, int, int, int, int, int, int, int, char, int );
void ClipCylinder( int, int, int, int, int, int, int, int, int, char, int );

/* void OutLinePolygon( Poly* ); */
/* void DrawPolygon( Poly* );    */
/* void ClipPolygon( Poly* );    */

void DashRibbon( Knot __far*, Knot __far*, int, int );
void StrandRibbon( Knot __far*, Knot __far*, int, int );
void SolidRibbon2( Knot __far*, Knot __far*, int, int );
void SolidRibbon( Knot __far*, Knot __far*, int );
void RectRibbon( Knot __far*, Knot __far*, int );
void DrawSphere( int, int, int, int, int );
void ClipSphere( int, int, int, int, int );
void DrawStar( int, int, int, int, int );
void ClipStar( int, int, int, int, int );

void ClipTwinLine( int, int, int, int, int, int, int, int, char );

void SetFontSize( int );
void DisplayRasString( int, int, int, unsigned char *, int );
void InitialisePixUtils( void );

