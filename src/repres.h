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

/* repres.h
 $Log: repres.h,v $
 Revision 1.1  2001/01/31 02:13:45  yaya
 Initial revision

 Revision 1.4  2000/08/27 16:09:43  yaya
 monitor dynamics extensions

 Revision 1.3  2000/08/26 18:13:00  yaya
 Updates to header comments in all files

 Revision 1.2  2000/08/09 01:18:38  yaya
 Rough cut with ucb

*/

#define DotMax    100
typedef struct _DotStruct {
        struct _DotStruct __far *next;
        short col[DotMax];
        Long xpos[DotMax];
        Long ypos[DotMax];
        Long zpos[DotMax];
        int count;
    } DotStruct;


typedef struct _Monitor {
        struct _Monitor *next;
        RAtom __far *src;
        RAtom __far *mid1;
        RAtom __far *mid2;
        RAtom __far *dst;
        int monmode;
        int dist;
        short col;
        unsigned char units;
    } Monitor;


typedef struct _Label {
        struct _Label *next;
        Long  refcount;
        char *label;
    } Label;



#ifdef REPRES
DotStruct __far *DotPtr;
Monitor *MonitList;
Label *LabelList;

int CartoonHeight;
int SolventDots;
int ProbeRadius;

int SurfaceChainsFlag;
int DotDensity;
int DotSize;
int DrawMonitDistance;
int DrawBetaArrows;

char LabelFormat[128];

#else
extern DotStruct __far *DotPtr;
extern Monitor *MonitList;
extern Label *LabelList;

extern int CartoonHeight;
extern int ProbeRadius;
extern int SolventDots;

extern int SurfaceChainsFlag;
extern int DotDensity;
extern int DotSize;
extern int DrawMonitDistance;
extern int DrawBetaArrows;

extern char LabelFormat[128];
#endif


int DeleteLabels( void );
void DeleteLabel( Label* );
Label *CreateLabel( char*, int );
void LabelTerminii( int );
void DefaultLabels( int );
void DefineLabels( char* );
void DisplayLabels( void );

void DeleteMonitors( void );
void AddMonitors2( RAtom __far*, RAtom __far*,
  RAtom __far*, RAtom __far*, 
  unsigned short, unsigned char, int );
void AddMonitors( RAtom __far*, RAtom __far* );
void CreateMonitor( Long, Long );
void DisplayMonitors( void );

void DeleteSurface( void );
void CalculateSurface( int );
void DisplaySurface( void );

/* Ribbons & Cartoons */
void DisplayRibbon( Chain __far* );

void ResetRepres( void );
void InitialiseRepres( void );

