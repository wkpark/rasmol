/* graphics.h
 * RasMol2 Molecular Graphics
 * Roger Sayle, February 1994
 * Version 2.3
 */

#ifdef IBMPC
#define CanvWidth  480
#define CanvHeight 480
#else
#define CanvWidth  576
#define CanvHeight 576
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
Pixel Lut[LutSize];
int PointX, PointY;
int XRange, WRange;
int YRange, HRange;
int UseHourGlass;
int MenuDisable;
int ReDrawFlag;
int Range;

#ifdef IBMPC
HWND CanvWin;
HBITMAP PixMap;
HPALETTE ColourMap;
LOGPALETTE __far *Palette;
#endif /* IBMPC */

Byte RLut[256];
Byte GLut[256];
Byte BLut[256];
Byte ULut[256];

#else /* GRAPHICS */
extern double DialValue[8];
extern Pixel Lut[LutSize];
extern int PointX, PointY;
extern int XRange, WRange;
extern int YRange, HRange;
extern int UseHourGlass;
extern int MenuDisable;
extern int ReDrawFlag;
extern int Range;

#ifdef IBMPC
extern HWND CanvWin;
extern HBITMAP PixMap;
extern HPALETTE ColourMap;
extern LOGPALETTE __far *Palette;
#endif /* IBMPC */

extern Byte RLut[256];
extern Byte GLut[256];
extern Byte BLut[256];
extern Byte ULut[256];

#ifdef __STDC__
void AllocateColourMap();
void ClearImage();
void TransferImage();
void UpdateScrollBars();
void SetMouseMode( int );
int LookUpColour( char*, int*, int*, int* );
void CloseDisplay();
void BeginWait();
void EndWait();

#ifdef IBMPC
int OpenDisplay( HANDLE, int );
int PrintImage();
#else
int OpenDisplay( int, int );
void NewMenu( int, char** );
int FetchEvent( int );
Pixel *CreateImage();
#endif /* IBMPC */

#else /* non-ANSI C compiler */
int OpenDisplay();
void AllocateColourMap();
void ClearImage();
void TransferImage();
void UpdateScrollBars();
void SetMouseMode();
int LookUpColour();
void CloseDisplay();
void BeginWait();
void EndWait();

#ifdef IBMPC
int PrintImage();
#else
void NewMenu();
int FetchEvent();
Pixel *CreateImage();
#endif /* IBMPC */

#endif
#endif /* GRAPHICS */

