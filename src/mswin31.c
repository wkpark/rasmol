/***************************************************************************
 *                            RasMol 2.7.1.1                               *
 *                                                                         *
 *                                RasMol                                   *
 *                 Molecular Graphics Visualisation Tool                   *
 *                            21 January 2001                              *
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
 *       Version 2.7.0, 2.7.1, 2.7.1.1 Mods by Herbert J. Bernstein        *
 *           Bernstein + Sons, P.O. Box 177, Bellport, NY, USA             *
 *                      yaya@bernstein-plus-sons.com                       *
 *           2.7.0 March 1999, 2.7.1 June 1999, 2.7.1.1 Jan 2001           *
 *              Copyright (C) Herbert J. Bernstein 1998-2001               *
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

/* mswin31.c
 */

#include <windows.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>
#include <stdio.h>


#define GRAPHICS
#include "rasmol.h"
#include "raswin.idm"
#include "graphics.h"
#include "langsel.h"


static int ColCount;
static BITMAPINFO __far *BitInfo;
static HCURSOR WaitCursor;
static HCURSOR OldCursor;
static HMENU hMenu;


void AllocateColourMap( void )
{
    register COLORREF ref;      
    register int i;
    
    if( ColourMap )
        DeleteObject(ColourMap);
        

    ColCount = 0;
    for( i=0; i<256; i++ )
        if( ULut[i] ) 
        {  Palette->palPalEntry[ColCount].peFlags = 0;
           Palette->palPalEntry[ColCount].peRed   = RLut[i];
           Palette->palPalEntry[ColCount].peGreen = GLut[i];
           Palette->palPalEntry[ColCount].peBlue  = BLut[i];

           BitInfo->bmiColors[ColCount].rgbBlue     = BLut[i];
           BitInfo->bmiColors[ColCount].rgbGreen    = GLut[i];
           BitInfo->bmiColors[ColCount].rgbRed      = RLut[i];
           BitInfo->bmiColors[ColCount].rgbReserved = 0;
           ColCount++;
        }   
    Palette->palNumEntries = ColCount;
    BitInfo->bmiHeader.biClrUsed = ColCount;   
    ColourMap = CreatePalette(Palette);

    for( i=0; i<256; i++ )
       if( ULut[i] )
       {   ref = RGB(RLut[i],GLut[i],BLut[i]);
           Lut[i] = GetNearestPaletteIndex(ColourMap,ref);
       }    
}



int CreateImage( void )
{
    register Long size;

    if( FBufHandle ) GlobalFree(FBufHandle);
    size = (Long)XRange*YRange*sizeof(Pixel)+16;
    FBufHandle = GlobalAlloc(GMEM_MOVEABLE,size);
    return (int)FBufHandle;
}


void TransferImage( void )
{
    HPALETTE OldCMap;
    HDC hDC;
        
    if( PixMap )
        DeleteObject(PixMap);

    BitInfo->bmiHeader.biWidth = XRange;
    BitInfo->bmiHeader.biHeight = YRange;
        
    hDC = GetDC(NULL);
    FBuffer = (Pixel  __huge*)GlobalLock(FBufHandle);
    /* CreateBitMap(XRange,YRange,1,8,FBuffer); */

    if( ColourMap )
    {   OldCMap = SelectPalette(hDC,ColourMap,FALSE);
        RealizePalette(hDC);  /* GDI Bug?? */
    }
        
    PixMap = CreateDIBitmap( hDC, (BITMAPINFOHEADER __far *)BitInfo, 
                             CBM_INIT, FBuffer, BitInfo, DIB_RGB_COLORS);
        
    if( ColourMap && OldCMap )                         
        SelectPalette(hDC,OldCMap,False);

    GlobalUnlock(FBufHandle);
    ReleaseDC(NULL,hDC);
    
    InvalidateRect(CanvWin,NULL,FALSE);
    UpdateWindow(CanvWin);
}


void ClearImage( void )
{
    HBRUSH hand;
    RECT rect;
    HDC hDC;
    
    hDC = GetDC(CanvWin);
    hand = CreateSolidBrush(RGB(RLut[0],GLut[0],BLut[0]));
    GetClientRect(CanvWin,&rect);
    FillRect(hDC,&rect,hand);
    ReleaseDC(CanvWin,hDC);
    DeleteObject(hand);

    if( PixMap )
    {   DeleteObject(PixMap);
        PixMap = NULL;
    }
}


int PrintImage( void )
{
    register char *device, *driver, *output;
    register int xsize, xres, yres;
    register int dx, dy, caps;
    char printer[255];
    PRINTDLG mypdlg;
    

    DOCINFO info;
    RECT rect;
    HDC hDC;

    GetProfileString("windows","device",",,,", printer, 255 );
    if( !(device = strtok(printer,",")) ||
        !(driver = strtok((char*)NULL,", ")) ||
        !(output = strtok((char*)NULL,", ")) ) {
        memset ((void *) &mypdlg, 0, sizeof(PRINTDLG));
        mypdlg.lStructSize = sizeof(PRINTDLG);
        mypdlg.hwndOwner = NULL;
        mypdlg.Flags = PD_RETURNDC;
        mypdlg.hInstance = NULL;
        
        if (PrintDlg(&mypdlg) == True) {
          if (mypdlg.hDevMode) GlobalFree (mypdlg.hDevMode);
          if (!mypdlg.hDC) return False;
          hDC = mypdlg.hDC;
        } else {
          WriteString("Unable to locate printer\n");
          return False;
        }
    } else {
      hDC = CreateDC(driver,device,output,NULL);
      if( !hDC ) {
        WriteString("Unable to locate printer\n");
        return False;
      }
    }

    caps = GetDeviceCaps( hDC, RASTERCAPS );
    if( !(caps & RC_STRETCHDIB) ) {
       WriteString("Unable to get necessary caps\n");
       return False;
    }
    
    xres = GetDeviceCaps( hDC, LOGPIXELSX );
    yres = GetDeviceCaps( hDC, LOGPIXELSY );
    xsize = GetDeviceCaps( hDC, HORZRES );

    dx = xsize - xres;
    dy = (int)(((long)dx*YRange)/XRange);

    /* Should set printer abort procedure */
    /* Position Image on Printed Page */
    rect.top = yres;        rect.bottom = rect.top + dy;
    rect.left = xres>>1;    rect.right = rect.left + dx;
    Escape( hDC, SET_BOUNDS, sizeof(RECT), (char __far*)&rect, NULL );

    /* Start RasWin Document */
    info.cbSize = sizeof(DOCINFO);
    info.lpszDocName = "RasWin";
    info.lpszOutput = NULL;
    StartDoc( hDC, &info );
    StartPage( hDC );
    


    BitInfo->bmiHeader.biWidth = XRange;
    BitInfo->bmiHeader.biHeight = YRange;
    FBuffer = (Pixel  __huge*)GlobalLock(FBufHandle);

    StretchDIBits( hDC, xres>>1, yres, dx, dy, 
                        0, 0, XRange, YRange, 
                        FBuffer, BitInfo, DIB_RGB_COLORS, SRCCOPY );

    GlobalUnlock(FBufHandle);

    EndPage( hDC );
    EndDoc( hDC );

    DeleteDC( hDC );
    return True;
}


int ClipboardImage( void )
{
    register BITMAPINFO __far *bitmap;
    register Pixel __huge *src;
    register Pixel __huge *dst;
    register long size,len;
    register HANDLE hand;
    register int i;


    if( OpenClipboard(CanvWin) )
    {   EmptyClipboard();

        /* SetClipboardData(CF_DIB,NULL);     */
        /* SetClipboardData(CF_PALETTE,NULL); */

        if( PixMap )
        {   len = (long)XRange*YRange*sizeof(Pixel);
            size = sizeof(BITMAPINFOHEADER) + 256*sizeof(RGBQUAD);
            if( (hand=GlobalAlloc(GHND,size+len)) )
            {   bitmap = (BITMAPINFO __far *)GlobalLock(hand);
                bitmap->bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
                bitmap->bmiHeader.biWidth = XRange;
                bitmap->bmiHeader.biHeight = YRange;
                bitmap->bmiHeader.biPlanes = 1;
                bitmap->bmiHeader.biBitCount = 8;
                bitmap->bmiHeader.biCompression = BI_RGB;
                bitmap->bmiHeader.biSizeImage = len;
                bitmap->bmiHeader.biXPelsPerMeter = 0;
                bitmap->bmiHeader.biYPelsPerMeter = 0;
                bitmap->bmiHeader.biClrImportant = 0;
                bitmap->bmiHeader.biClrUsed = 0;

                for( i=0; i<256; i++ )
                    if( ULut[i] )
                    {   bitmap->bmiColors[Lut[i]].rgbBlue  = BLut[i];
                        bitmap->bmiColors[Lut[i]].rgbGreen = GLut[i];
                        bitmap->bmiColors[Lut[i]].rgbRed   = RLut[i];
                    }


                src = (Pixel __huge*)GlobalLock(FBufHandle);
                dst = ((Pixel __huge*)bitmap)+size;

                /* Transfer the frame buffer */
                while( len-- ) *dst++ = *src++;

                GlobalUnlock(FBufHandle);
                GlobalUnlock(hand);
                SetClipboardData(CF_DIB,hand);
            }
        }

        if( ColourMap )
        {   if( (hand = CreatePalette(Palette)) )
                SetClipboardData(CF_PALETTE,hand);
        }
        CloseClipboard();
        return True;
    } else return False;
}


void SetCanvasTitle( char *ptr )
{
    SetWindowText(CanvWin,ptr);
}


void UpdateScrollBars( void )
{
    register int pos;
    
    pos = 50-(int)(50.0*DialValue[0]);
    SetScrollPos(CanvWin,SB_VERT,pos,TRUE);
    
    pos = (int)(50.0*DialValue[1])+50;
    SetScrollPos(CanvWin,SB_HORZ,pos,TRUE);
}

void ReDrawWindow( void )
{
 
  HMENU mentop[7];
  int ii;

  ModifyMenu(hMenu,IDM_OPEN,     MF_STRING, IDM_OPEN,    MsgStrs[StrMOpen]);
  ModifyMenu(hMenu,IDM_INFO,     MF_STRING, IDM_INFO,    MsgStrs[StrMInfo]);
  ModifyMenu(hMenu,IDM_CLOSE,    MF_STRING, IDM_CLOSE,   MsgStrs[StrMClose]);
  ModifyMenu(hMenu,IDM_PRINT,    MF_STRING, IDM_PRINT,   MsgStrs[StrMPrint]);
  ModifyMenu(hMenu,IDM_SETUP,    MF_STRING, IDM_SETUP,   MsgStrs[StrMPSetup]);
  ModifyMenu(hMenu,IDM_EXIT,     MF_STRING, IDM_EXIT,    MsgStrs[StrMExit]);

  ModifyMenu(hMenu,IDM_SELECT,   MF_STRING, IDM_SELECT,  MsgStrs[StrMSelAll]);
  ModifyMenu(hMenu,IDM_CUT,      MF_STRING, IDM_CUT,     MsgStrs[StrMCut]);
  ModifyMenu(hMenu,IDM_COPY,     MF_STRING, IDM_COPY,    MsgStrs[StrMCopy]);
  ModifyMenu(hMenu,IDM_PASTE,    MF_STRING, IDM_PASTE,   MsgStrs[StrMPaste]);
  ModifyMenu(hMenu,IDM_DELETE,   MF_STRING, IDM_DELETE,  MsgStrs[StrMDelete]);

  ModifyMenu(hMenu,IDM_WIREFRAME,MF_STRING,IDM_WIREFRAME,MsgStrs[StrMWirefr]);
  ModifyMenu(hMenu,IDM_BACKBONE, MF_STRING,IDM_BACKBONE, MsgStrs[StrMBackbn]);
  ModifyMenu(hMenu,IDM_STICKS,   MF_STRING,IDM_STICKS,   MsgStrs[StrMSticks]);
  ModifyMenu(hMenu,IDM_SPHERES,  MF_STRING,IDM_SPHERES,  MsgStrs[StrMSpacefl]);
  ModifyMenu(hMenu,IDM_BALLSTICK,MF_STRING,IDM_BALLSTICK,MsgStrs[StrMBallStk]);
  ModifyMenu(hMenu,IDM_RIBBONS,  MF_STRING,IDM_RIBBONS,  MsgStrs[StrMRibbons]);
  ModifyMenu(hMenu,IDM_STRANDS,  MF_STRING,IDM_STRANDS,  MsgStrs[StrMStrands]);
  ModifyMenu(hMenu,IDM_CARTOONS, MF_STRING,IDM_CARTOONS, MsgStrs[StrMCartoon]);

  ModifyMenu(hMenu,IDM_MONO,     MF_STRING,IDM_MONO,     MsgStrs[StrMMonochr]);
  ModifyMenu(hMenu,IDM_CPK,      MF_STRING,IDM_CPK,      MsgStrs[StrMCPK]);
  ModifyMenu(hMenu,IDM_SHAPELY,  MF_STRING,IDM_SHAPELY,  MsgStrs[StrMShapely]);
  ModifyMenu(hMenu,IDM_GROUP,    MF_STRING,IDM_GROUP,    MsgStrs[StrMGroup]);
  ModifyMenu(hMenu,IDM_CHAIN,    MF_STRING,IDM_CHAIN,    MsgStrs[StrMChain]);
  ModifyMenu(hMenu,IDM_TEMPER,   MF_STRING,IDM_TEMPER,   MsgStrs[StrMTemp]);
  ModifyMenu(hMenu,IDM_STRUCT,   MF_STRING,IDM_STRUCT,   MsgStrs[StrMStruct]);
  ModifyMenu(hMenu,IDM_USER,     MF_STRING,IDM_USER,     MsgStrs[StrMUser]);
  ModifyMenu(hMenu,IDM_MODEL,    MF_STRING,IDM_MODEL,    MsgStrs[StrMModel]);
  ModifyMenu(hMenu,IDM_ALT,      MF_STRING,IDM_ALT,      MsgStrs[StrMAlt]);

  ModifyMenu(hMenu,IDM_SLAB,     MF_STRING,IDM_SLAB,     MsgStrs[StrMSlab]);
  ModifyMenu(hMenu,IDM_HYDROGEN, MF_STRING,IDM_HYDROGEN, MsgStrs[StrMHydr]);
  ModifyMenu(hMenu,IDM_HETERO,   MF_STRING,IDM_HETERO,   MsgStrs[StrMHet]);
  ModifyMenu(hMenu,IDM_SPECULAR, MF_STRING,IDM_SPECULAR, MsgStrs[StrMSpec]);
  ModifyMenu(hMenu,IDM_SHADOW,   MF_STRING,IDM_SHADOW,   MsgStrs[StrMShad]);
  ModifyMenu(hMenu,IDM_STEREO,   MF_STRING,IDM_STEREO,   MsgStrs[StrMStereo]);
  ModifyMenu(hMenu,IDM_LABELS,   MF_STRING,IDM_LABELS,   MsgStrs[StrMLabel]);

  ModifyMenu(hMenu,IDM_BMP,      MF_STRING,IDM_BMP,      MsgStrs[StrMBMP]);
  ModifyMenu(hMenu,IDM_GIF,      MF_STRING,IDM_GIF,      MsgStrs[StrMGIF]);
  ModifyMenu(hMenu,IDM_EPSF,     MF_STRING,IDM_EPSF,     MsgStrs[StrMPostscr]);
  ModifyMenu(hMenu,IDM_PPM,      MF_STRING,IDM_PPM,      MsgStrs[StrMPPM]);
  ModifyMenu(hMenu,IDM_RAST,     MF_STRING,IDM_RAST,     MsgStrs[StrMSRast]);

  ModifyMenu(hMenu,IDM_ABOUT,    MF_STRING,IDM_ABOUT,    MsgStrs[StrMAbout]);
  ModifyMenu(hMenu,IDM_HELP,     MF_STRING,IDM_HELP,     MsgStrs[StrMUserM]);
  
  for (ii = 7; ii > 0; ii--) {
    mentop[ii-1] = GetSubMenu(hMenu,ii-1);
    RemoveMenu(hMenu,ii-1,MF_BYPOSITION);
  }
  for (ii = 0; ii < 7; ii++) {
    AppendMenu(hMenu, MF_POPUP | MF_STRING , (UINT) mentop[ii], MsgStrs[StrMFile+ii]);
  }
  
  DrawMenuBar(CanvWin);

}

void SetMouseUpdateStatus( int bool )
{
    MouseUpdateStatus = bool;
}
                         
                         
void SetMouseCaptureStatus( int bool )
{
    if( bool )
    {   if( !MouseCaptureStatus )
            SetCapture(CanvWin);
    } else
        if( MouseCaptureStatus )
            ReleaseCapture();
    MouseCaptureStatus = bool;
}


int LookUpColour( char *name, int *r, int *g, int *b )
{
    UnusedArgument(name);
    UnusedArgument(r);
    UnusedArgument(g);
    UnusedArgument(b);

    return False;
}    


void EnableMenus( int flag )
{
    if( flag )
    {   SetMenu(CanvWin,hMenu);
    } else SetMenu(CanvWin,0);
    DisableMenu = !flag;
}


int OpenDisplay( HANDLE instance, int mode )
{
    register int i,size;
    long style;
    RECT rect;
    static char VersionStr[50];

    sprintf (VersionStr,"RasMol Version %s", VERSION);

    PixMap = NULL;
    ColourMap = NULL;

    MouseCaptureStatus = False;
    MouseUpdateStatus = False;   
    UseHourGlass = True;
    DisableMenu = False;

    for( i=0; i<8; i++ )
         DialValue[i] = 0.0;

    ULut[0] = True;
    RLut[0] = GLut[0] = BLut[0] = 0;
    XRange = DefaultWide;   WRange = XRange>>1;
    YRange = DefaultHigh;   HRange = YRange>>1;
    Range = MinFun(XRange,YRange);
    
    rect.top  = 0;   rect.bottom = YRange;
    rect.left = 0;   rect.right  = XRange;

    
    style = WS_OVERLAPPEDWINDOW | WS_HSCROLL | WS_VSCROLL;

    hMenu = LoadMenu(instance,"RasWinMenu");
    AdjustWindowRect(&rect,style,TRUE);
    CanvWin = CreateWindow("RasWinClass",VersionStr, style,
                            CW_USEDEFAULT, CW_USEDEFAULT,
                            rect.right-rect.left, 
                            rect.bottom-rect.top,
                            NULL,hMenu,instance,NULL);
                            
    if( !CanvWin) return False;

    size = sizeof(LOGPALETTE) + 256*sizeof(PALETTEENTRY);
    Palette = (LOGPALETTE __far*)_fmalloc( size );
    size = sizeof(BITMAPINFOHEADER) + 256*sizeof(RGBQUAD);
    BitInfo = (BITMAPINFO __far*)_fmalloc( size );

    if( !Palette || !BitInfo )
        return False;
   
    WaitCursor = LoadCursor(NULL,IDC_WAIT);
        
    Palette->palVersion = 0x300;   
    
    BitInfo->bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
    BitInfo->bmiHeader.biCompression = BI_RGB;
    BitInfo->bmiHeader.biXPelsPerMeter = 0;
    BitInfo->bmiHeader.biYPelsPerMeter = 0;
    BitInfo->bmiHeader.biClrImportant = 0;
    BitInfo->bmiHeader.biSizeImage = 0;
    BitInfo->bmiHeader.biBitCount = 8;
    BitInfo->bmiHeader.biPlanes = 1;

    /* Initialise Palette! */
    for( i=1; i<256; i++ )
        ULut[i] = False;
    AllocateColourMap();

    ShowWindow(CanvWin,mode);
    UpdateScrollBars();
    UpdateWindow(CanvWin);
    return True;                       
}

    
void BeginWait( void )
{
    if( UseHourGlass )
        OldCursor = SetCursor(WaitCursor);
}


void EndWait( void )
{
    if( UseHourGlass )
        SetCursor(OldCursor);
}


void CloseDisplay( void )
{
    if( ColourMap )
        DeleteObject(ColourMap);
    if( PixMap )
        DeleteObject(PixMap);
}

