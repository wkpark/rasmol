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

/* graphics.h
 */

#ifdef APPLEMAC      
#define DefaultWide  400
#define DefaultHigh  400
#endif

#ifdef IBMPC
#define DefaultWide  480
#define DefaultHigh  480
#endif

#ifndef DefaultWide
#define DefaultWide  576
#define DefaultHigh  576
#endif


#ifdef EIGHTBIT
#define LutSize  256
#else
#define LutSize  1024
#endif


#define RFRotateX  0x0001
#define RFRotateY  0x0002
#define RFRotateZ  0x0004
#define RFZoom     0x0008
#define RFTransX   0x0010
#define RFTransY   0x0020
#define RFTransZ   0x0040
#define RFSlab     0x0080
#define RFReSize   0x0100
#define RFColour   0x0200
#define RFRefresh  0x0400
#define RFRotBond  0x0800

#define RFTrans    0x0070
#define RFRotate   0x0807
#define RFApply    0x097F
#define RFDials    0x00FF
#define RFMagnify  0x0108
#define RFInitial  0x09FF


#define DialRX      0
#define DialRY      1
#define DialRZ      2
#define DialZoom    3
#define DialTX      4
#define DialTY      5
#define DialTZ      6
#define DialSlab    7
#define DialBClip   8
#define DialBRot    9



#ifdef GRAPHICS
double DialValue[10];
int XRange, WRange;
int YRange, HRange;
int ZRange;
int UseHourGlass;
int DisableMenu;
int ReDrawFlag;
int Range;

int MouseCaptureStatus;
int MouseUpdateStatus;

Pixel __huge *FBuffer;
short __huge *DBuffer;

Pixel Lut[LutSize];
Byte RLut[LutSize];
Byte GLut[LutSize];
Byte BLut[LutSize];
Byte ULut[LutSize];


#ifdef MSWIN
LOGPALETTE __far *Palette;
HPALETTE ColourMap;
HGLOBAL FBufHandle;
HGLOBAL DBufHandle;
HBITMAP PixMap;
HWND CanvWin;
#endif /* MSWIN */

#ifdef APPLEMAC
ControlHandle HScroll;
ControlHandle VScroll;
CursHandle CanvCursor;
CursHandle CmndCursor;
CursHandle WaitCursor;
WindowPtr CanvWin;
WindowPtr CmndWin;
THPrint PrintHand;
Handle FBufHandle;
Handle DBufHandle;
#endif /* APPLEMAC */

#ifdef X11WIN
Display *dpy;
Pixmap icon;
Visual *vis;
int Monochrome;
int PixDepth;
int LocalMap;
Colormap cmap;
Window PopUpWin;
Window MainWin;
Window CanvWin;
Window MenuWin;
Window RootWin;

#endif

#else /* GRAPHICS */
extern double DialValue[10];
extern int XRange, WRange;
extern int YRange, HRange;
extern int ZRange;
extern int UseHourGlass;
extern int DisableMenu;
extern int ReDrawFlag;
extern int Range;

extern int MouseCaptureStatus;
extern int MouseUpdateStatus;

extern Pixel __huge *FBuffer;
extern short __huge *DBuffer;

extern Pixel Lut[LutSize];
extern Byte RLut[LutSize];
extern Byte GLut[LutSize];
extern Byte BLut[LutSize];
extern Byte ULut[LutSize];


#ifdef MSWIN
extern LOGPALETTE __far *Palette;
extern HPALETTE ColourMap;
extern HGLOBAL FBufHandle;
extern HGLOBAL DBufHandle;
extern HBITMAP PixMap;
extern HWND CanvWin;
#endif /* MSWIN */

#ifdef APPLEMAC
extern ControlHandle HScroll;
extern ControlHandle VScroll;
extern CursHandle CanvCursor;
extern CursHandle CmndCursor;
extern CursHandle WaitCursor;
extern WindowPtr CanvWin;
extern WindowPtr CmndWin;
extern THPrint PrintHand;
extern Handle FBufHandle;
extern Handle DBufHandle;
#endif /* APPLEMAC */
#endif

int CreateImage( void );
void TransferImage( void );
int ClipboardImage( void );
void ClearImage( void );
int PrintImage( void );

void AllocateColourMap( void );
void UpdateScrollBars( void );
int LookUpColour( char*, int*, int*, int* );
void SetMouseUpdateStatus( int );
void SetMouseCaptureStatus( int );
void SetCanvasTitle( char* );
void ReDrawWindow( void );
void EnableMenus( int );
void CloseDisplay( void );
void BeginWait( void );
void EndWait( void );

#ifdef MSWIN
int OpenDisplay( HANDLE, int );
#else
int OpenDisplay( int, int );
#endif

#if !defined(IBMPC) && !defined(APPLEMAC)
int FetchEvent( int );
#endif

#ifdef X11WIN
void FatalGraphicsError( char *ptr );
#endif

