/* graphics.h
 * RasMol2 Molecular Graphics
 * Roger Sayle, October 1994
 * Version 2.5
 */

#ifdef APPLEMAC      
#define InitialWide  400
#define InitialHigh  400
#endif

#ifdef IBMPC
#define InitialWide  480
#define InitialHigh  480
#endif

#ifndef InitialWide
#define InitialWide  576
#define InitialHigh  576
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
#define RFPoint    0x0200
#define RFColour   0x0400
#define RFRefresh  0x0800

#define RFTrans    0x0070
#define RFRotate   0x0007
#define RFApply    0x017F
#define RFDials    0x00FF
#define RFMagnify  0x0108
#define RFInitial  0x01FF


#define MMRasMol   0x00
#define MMInsight  0x01
#define MMQuanta   0x02


#define ButMax   8


#ifdef GRAPHICS
double DialValue[8];
int WinHigh, WinWide;
int PointX, PointY;
int XRange, WRange;
int YRange, HRange;
int UseHourGlass;
int DisableMenu;
int ReDrawFlag;
int MouseMode;
int Range;

Pixel Lut[LutSize];
Byte RLut[LutSize];
Byte GLut[LutSize];
Byte BLut[LutSize];
Byte ULut[LutSize];


#ifdef IBMPC
HWND CanvWin;
HBITMAP PixMap;
HPALETTE ColourMap;
LOGPALETTE __far *Palette;
#endif /* IBMPC */

#ifdef APPLEMAC
ControlHandle HScroll;
ControlHandle VScroll;
CursHandle CanvCursor;
CursHandle CmndCursor;
CursHandle WaitCursor;
WindowPtr CanvWin;
WindowPtr CmndWin;
THPrint PrintHand;
#endif /* APPLEMAC */

#else /* GRAPHICS */
extern double DialValue[8];
extern int WinHigh, WinWide;
extern int PointX, PointY;
extern int XRange, WRange;
extern int YRange, HRange;
extern int UseHourGlass;
extern int DisableMenu;
extern int ReDrawFlag;
extern int MouseMode;
extern int Range;

extern Pixel Lut[LutSize];
extern Byte RLut[LutSize];
extern Byte GLut[LutSize];
extern Byte BLut[LutSize];
extern Byte ULut[LutSize];


#ifdef IBMPC
extern HWND CanvWin;
extern HBITMAP PixMap;
extern HPALETTE ColourMap;
extern LOGPALETTE __far *Palette;
#endif /* IBMPC */

#ifdef APPLEMAC
extern ControlHandle HScroll;
extern ControlHandle VScroll;
extern CursHandle CanvCursor;
extern CursHandle CmndCursor;
extern CursHandle WaitCursor;
extern WindowPtr CanvWin;
extern WindowPtr CmndWin;
extern THPrint PrintHand;
#endif /* APPLEMAC */


#ifdef FUNCPROTO
int CreateImage();
void TransferImage();
int ClipboardImage();
void ClearImage();
int PrintImage();

void AllocateColourMap();
void UpdateScrollBars();
int LookUpColour( char*, int*, int*, int* );
void SetMouseMode( int );
void EnableMenus( int );
void CloseDisplay();
void BeginWait();
void EndWait();

#ifdef IBMPC
int OpenDisplay( HANDLE, int );
#else
int OpenDisplay( int, int );
#endif

#if !defined(IBMPC) && !defined(APPLEMAC)
int FetchEvent( int );
#endif

#else /* non-ANSI C compiler */
int CreateImage();
void TransferImage();
int ClipboardImage();
void ClearImage();
int PrintImage();

int OpenDisplay();
void AllocateColourMap();
void UpdateScrollBars();
int LookUpColour();
void SetMouseMode();
void EnableMenus();
void CloseDisplay();
void BeginWait();
void EndWait();

#if !defined(IBMPC) && !defined(APPLEMAC)
int FetchEvent();
#endif

#endif
#endif /* GRAPHICS */

