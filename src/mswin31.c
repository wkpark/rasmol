/* mswin31.c
 * RasMol2 Molecular Graphics
 * Roger Sayle, February 1994
 * Version 2.3
 */
#include <windows.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>
#include <stdio.h>


#define GRAPHICS
#include "rasmol.h"
#include "graphics.h"
#include "render.h"


static BITMAPINFO __far *BitInfo;
static HCURSOR WaitCursor;
static HCURSOR OldCursor;
static int ColCount;


void AllocateColourMap()
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


void ClearImage()
{
    RECT rect;
    HDC hDC;
    
    hDC = GetDC(CanvWin);
    GetClientRect(CanvWin,&rect);
    FillRect(hDC,&rect,GetStockObject(BLACK_BRUSH));
    ReleaseDC(CanvWin,hDC);

    if( PixMap )
    {   DeleteObject(PixMap);
        PixMap = NULL;
    }
}


void TransferImage()
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
        
    if( ColourMap )                         
        SelectPalette(hDC,OldCMap,False);

    GlobalUnlock(FBufHandle);
    ReleaseDC(NULL,hDC);
    
    InvalidateRect(CanvWin,NULL,FALSE);
    UpdateWindow(CanvWin);
}


int PrintImage()
{
    register char *device, *driver, *output;
    register int xsize, ysize, xres, yres;
    register int dx, dy, caps;
    char printer[80];

    DOCINFO info;
    RECT rect;
    HDC hDC;

    GetProfileString("windows","device", "", printer, 80 );
    if( !(device = strtok(printer,",")) ) return( False );
    if( !(driver = strtok((char*)NULL,", ")) ) return( False );
    if( !(output = strtok((char*)NULL,", ")) ) return( False );

    hDC = CreateDC(driver,device,output,NULL);
    if( !hDC ) return( False );

    caps = GetDeviceCaps( hDC, RASTERCAPS );
    if( !(caps & RC_STRETCHDIB) ) return( False );
    
    xres = GetDeviceCaps( hDC, LOGPIXELSX );
    yres = GetDeviceCaps( hDC, LOGPIXELSY );
    xsize = GetDeviceCaps( hDC, HORZRES );
    ysize = GetDeviceCaps( hDC, VERTRES );

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
    return( True );
}


void UpdateScrollBars()
{
    register int pos;
    
    pos = 50-(int)(50.0*DialValue[0]);
    SetScrollPos(CanvWin,SB_VERT,pos,TRUE);
    
    pos = (int)(50.0*DialValue[1])+50;
    SetScrollPos(CanvWin,SB_HORZ,pos,TRUE);
}


int LookUpColour( name, r, g, b )
    char *name; int *r, *g, *b;
{
    return( False );
}    


int OpenDisplay( instance, mode )
    HANDLE instance; int mode;
{
    register int i,size;
    long style;
    RECT rect;

    PixMap = NULL;
    ColourMap = NULL;
    UseHourGlass = True;

    for( i=0; i<8; i++ )
         DialValue[i] = 0.0;


    RLut[0] = GLut[0] = BLut[0] = 0;
    XRange = CanvWidth;   WRange = XRange>>1;
    YRange = CanvHeight;  HRange = YRange>>1;
    Range = MinFun(XRange,YRange);
    
    rect.top  = 0;   rect.bottom = YRange;
    rect.left = 0;   rect.right  = XRange;

    
    style = WS_OVERLAPPEDWINDOW | WS_HSCROLL | WS_VSCROLL;
    
    
    AdjustWindowRect(&rect,style,TRUE);
    CanvWin = CreateWindow("RasWinClass","RasMol Version 2.3", style,
			    CW_USEDEFAULT, CW_USEDEFAULT,
                            rect.right-rect.left, 
                            rect.bottom-rect.top,
                            NULL,NULL,instance,NULL);
                            

    WaitCursor = LoadCursor(NULL,IDC_WAIT);

    size = sizeof(LOGPALETTE) + 256*sizeof(PALETTEENTRY);
    Palette = (LOGPALETTE __far*)_fmalloc( size );
    size = sizeof(BITMAPINFOHEADER) + 256*sizeof(RGBQUAD);
    BitInfo = (BITMAPINFO __far*)_fmalloc( size );


    if( !CanvWin || !Palette || !BitInfo )
        return( False );
        
    Palette->palVersion = 0x300;   
    
    BitInfo->bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
    BitInfo->bmiHeader.biCompression = BI_RGB;
    BitInfo->bmiHeader.biXPelsPerMeter = 0;
    BitInfo->bmiHeader.biYPelsPerMeter = 0;
    BitInfo->bmiHeader.biClrImportant = 0;
    BitInfo->bmiHeader.biSizeImage = 0;
    BitInfo->bmiHeader.biBitCount = 8;
    BitInfo->bmiHeader.biPlanes = 1;
    
  
    ShowWindow(CanvWin,mode);
    UpdateScrollBars();
    UpdateWindow(CanvWin);
    return(True);                       
}

    
void BeginWait()
{
    if( UseHourGlass )
        OldCursor = SetCursor(WaitCursor);
}


void EndWait()
{
    if( UseHourGlass )
        SetCursor(OldCursor);
}


void CloseDisplay()
{
    if( ColourMap )
        DeleteObject(ColourMap);
    if( PixMap )
        DeleteObject(PixMap);
}

