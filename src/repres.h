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

/* repres.h
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
        Atom __far *src;
        Atom __far *dst;
        unsigned short dist;
        short col;
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
int DrawMonitDistance;
int DrawBetaArrows;

#else
extern DotStruct __far *DotPtr;
extern Monitor *MonitList;
extern Label *LabelList;

extern int CartoonHeight;
extern int ProbeRadius;
extern int SolventDots;

extern int SurfaceChainsFlag;
extern int DrawMonitDistance;
extern int DrawBetaArrows;
#endif


int DeleteLabels( void );
void DeleteLabel( Label* );
Label *CreateLabel( char*, int );
void LabelTerminii( int );
void DefaultLabels( int );
void DefineLabels( char* );
void DisplayLabels( void );

void DeleteMonitors( void );
void AddMonitors( Atom __far*, Atom __far* );
void CreateMonitor( Long, Long );
void DisplayMonitors( void );

void DeleteSurface( void );
void CalculateSurface( int );
void DisplaySurface( void );

/* Ribbons & Cartoons */
void DisplayRibbon( Chain __far* );

void ResetRepres();
void InitialiseRepres();

