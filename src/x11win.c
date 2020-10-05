/***************************************************************************
 *                            RasMol 2.7.1.1                               *
 *                                                                         *
 *                                RasMol                                   *
 *                 Molecular Graphics Visualisation Tool                   *
 *                            17 January 2001                              *
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

/* x11win.c
 */

#ifndef sun386
#include <stdlib.h>
#endif
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>

#ifdef VMS
#include <in.h>
#endif

#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xatom.h>
#include <X11/keysym.h>
#include <X11/cursorfont.h>

#define GRAPHICS
#include "rasmol.h"
#include "graphics.h"
#include "bitmaps.h"
#include "command.h"
#include "cmndline.h"
#include "langsel.h"


/* Menu Definitions */
#define mbEnable    0x01
#define mbOption    0x02
#define mbCheck     0x04
#define mbSepBar    0x08
#define mbAccel     0x10


typedef struct _MenuItem {
            char **text;
            int flags;
            int *pos;
            int *len;
        } MenuItem;

static MenuItem FilMenu[5] = {
    { &MsgStrs[StrMOpen],   0x11,&MsgAuxl[StrMOpen],   &MsgLens[StrMOpen]    },
    { &MsgStrs[StrMSaveAs], 0x11,&MsgAuxl[StrMSaveAs], &MsgLens[StrMSaveAs]  },
    { &MsgStrs[StrMClose],  0x11,&MsgAuxl[StrMClose],  &MsgLens[StrMClose]   },
    { &MsgStrs[StrMEmpty],  0x08,&MsgAuxl[StrMEmpty],  &MsgLens[StrMEmpty]   },
    { &MsgStrs[StrMExit],   0x11,&MsgAuxl[StrMExit],   &MsgLens[StrMExit]    }};

static MenuItem DisMenu[8] = {
    { &MsgStrs[StrMWirefr], 0x11,&MsgAuxl[StrMWirefr], &MsgLens[StrMWirefr]  },
    { &MsgStrs[StrMBackbn], 0x11,&MsgAuxl[StrMBackbn], &MsgLens[StrMBackbn]  },
    { &MsgStrs[StrMSticks], 0x11,&MsgAuxl[StrMSticks], &MsgLens[StrMSticks]  },
    { &MsgStrs[StrMSpacefl],0x11,&MsgAuxl[StrMSpacefl],&MsgLens[StrMSpacefl] },
    { &MsgStrs[StrMBallStk],0x11,&MsgAuxl[StrMBallStk],&MsgLens[StrMBallStk] },
    { &MsgStrs[StrMRibbons],0x11,&MsgAuxl[StrMRibbons],&MsgLens[StrMRibbons] },
    { &MsgStrs[StrMStrands],0x11,&MsgAuxl[StrMStrands],&MsgLens[StrMStrands] },
    { &MsgStrs[StrMCartoon],0x11,&MsgAuxl[StrMCartoon],&MsgLens[StrMCartoon] }};

static MenuItem ColMenu[10] = {
    { &MsgStrs[StrMMonochr],0x11,&MsgAuxl[StrMMonochr],&MsgLens[StrMMonochr] },
    { &MsgStrs[StrMCPK],    0x11,&MsgAuxl[StrMCPK],    &MsgLens[StrMCPK]     },
    { &MsgStrs[StrMShapely],0x11,&MsgAuxl[StrMShapely],&MsgLens[StrMShapely] },
    { &MsgStrs[StrMGroup],  0x11,&MsgAuxl[StrMGroup],  &MsgLens[StrMGroup]   },
    { &MsgStrs[StrMChain],  0x11,&MsgAuxl[StrMChain],  &MsgLens[StrMChain]   },
    { &MsgStrs[StrMTemp],   0x11,&MsgAuxl[StrMTemp],   &MsgLens[StrMTemp]    },
    { &MsgStrs[StrMStruct], 0x11,&MsgAuxl[StrMStruct], &MsgLens[StrMStruct]  },
    { &MsgStrs[StrMUser],   0x11,&MsgAuxl[StrMUser],   &MsgLens[StrMUser]    }, 
    { &MsgStrs[StrMModel],  0x11,&MsgAuxl[StrMModel],  &MsgLens[StrMModel]   },
    { &MsgStrs[StrMAlt],    0x11,&MsgAuxl[StrMAlt],    &MsgLens[StrMAlt]     }};

static MenuItem OptMenu[7] = {
    { &MsgStrs[StrMSlab],   0x13,&MsgAuxl[StrMSlab],   &MsgLens[StrMSlab]    },
    { &MsgStrs[StrMHydr],   0x17,&MsgAuxl[StrMHydr],   &MsgLens[StrMHydr]    },
    { &MsgStrs[StrMHet],    0x17,&MsgAuxl[StrMHet],    &MsgLens[StrMHet]     },
    { &MsgStrs[StrMSpec],   0x13,&MsgAuxl[StrMSpec],   &MsgLens[StrMSpec]    },
    { &MsgStrs[StrMShad],   0x13,&MsgAuxl[StrMShad],   &MsgLens[StrMShad]    },
    { &MsgStrs[StrMStereo], 0x13,&MsgAuxl[StrMStereo], &MsgLens[StrMStereo]  },
    { &MsgStrs[StrMLabel],  0x13,&MsgAuxl[StrMLabel],  &MsgLens[StrMLabel]   }};

static MenuItem ExpMenu[7] = {
    { &MsgStrs[StrMGIF],    0x11,&MsgAuxl[StrMGIF],    &MsgLens[StrMGIF]     },
    { &MsgStrs[StrMPostscr],0x11,&MsgAuxl[StrMPostscr],&MsgLens[StrMPostscr] },
    { &MsgStrs[StrMPPM],    0x11,&MsgAuxl[StrMPPM],    &MsgLens[StrMPPM]     },
    { &MsgStrs[StrMIRGB],   0x11,&MsgAuxl[StrMIRGB],   &MsgLens[StrMIRGB]    },
    { &MsgStrs[StrMSRast],  0x11,&MsgAuxl[StrMSRast],  &MsgLens[StrMSRast]   },
    { &MsgStrs[StrMBMP],    0x11,&MsgAuxl[StrMBMP],    &MsgLens[StrMBMP]     },
    { &MsgStrs[StrMPICT],   0x11,&MsgAuxl[StrMPICT],   &MsgLens[StrMPICT]    }};

static MenuItem HelMenu[2] = {
    { &MsgStrs[StrMAbout],  0x10,&MsgAuxl[StrMAbout],  &MsgLens[StrMAbout]   },
    { &MsgStrs[StrMUserM],  0x10,&MsgAuxl[StrMUserM],  &MsgLens[StrMUserM]   }};


typedef struct _BarItem {
            MenuItem *menu;
            char **text;
            int count;
            int flags;
            int *pos;
            int *len;
        } BarItem;

#define MenuBarMax 6
static BarItem MenuBar[MenuBarMax] = { 
    { FilMenu,  &MsgStrs[StrMFile],     5, 0x01, &MsgAuxl[StrMFile],
                                                 &MsgLens[StrMFile] },
    { DisMenu,  &MsgStrs[StrMDisplay],  8, 0x01, &MsgAuxl[StrMDisplay],
                                                 &MsgLens[StrMDisplay] },
    { ColMenu,  &MsgStrs[StrMColour],  10, 0x01, &MsgAuxl[StrMColour],
                                                 &MsgLens[StrMColour] },
    { OptMenu,  &MsgStrs[StrMOpt],      7, 0x01, &MsgAuxl[StrMOpt],
                                                 &MsgLens[StrMOpt] },
    { ExpMenu,  &MsgStrs[StrMExport],   7, 0x01, &MsgAuxl[StrMExport],
                                                 &MsgLens[StrMExport] },
    { HelMenu,  &MsgStrs[StrMHelp],     2, 0x01, &MsgAuxl[StrMHelp],
                                                 &MsgLens[StrMHelp] } };

static int MenuFocus;
static int ItemFocus;
static int MenuItemSelect;
static int MenuBarSelect;
static int MenuBarCount;
static int PopUpWide;
static int PopUpHigh;
static int PopUpFlag;
static int ItemFlag;


#ifdef DIALBOX
#include <X11/extensions/XInput.h>

static char *DialLabel[] = { "ROTATE X", "ROTATE Y", "ROTATE Z", "  ZOOM  ",
                             "TRANS X ", "TRANS Y ", "TRANS Z ", "  SLAB  " };
static int *DialMap;
static int ESVDialMap[8] = { 0, 1, 2, 3, 4, 5, 6, 7 };
static int SGIDialMap[8] = { 3, 7, 2, 6, 1, 5, 0, 4 };

static Real DialRes[8];
static int DialPrev[8];
static int DialMode;

static int UseDialLEDs;
static XDevice *Dials;
static int DialEvent;
static int UseDials;
#endif


#ifdef MITSHM
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <X11/extensions/XShm.h>

XShmSegmentInfo xshminfo;
int SharedMemOption;
int SharedMemFlag;
#endif

#define XScrlDial  1 /*1*/
#define YScrlDial  0 /*0*/
#define XScrlSkip  8
#define YScrlSkip  8

typedef union {
    Long longword;
    Byte bytes[4];
    } ByteTest;


static int MenuHigh;
static int FontHigh;

static Cursor cross;
static Cursor arrow;
static Cursor hglass;
static Pixmap Scrl;
static Pixmap tilepix;
static Pixmap uppix, dnpix;
static Pixmap lfpix, rgpix;
static XFontStruct *MenuFont;
static XSetWindowAttributes attr;
static Window XScrlWin, YScrlWin;
static Window PopUpWin;
static Window MainWin;
static Window CanvWin;
static Window MenuWin;
static Window RootWin;
static XWMHints hints;
static Colormap cmap;
static Colormap lmap;
static XImage *image;
static Display *dpy;
static Visual *vis;
static GC gcon;

#ifdef EIGHTBIT
static unsigned long Ident[256];
static int IdentCount;
#endif

static int HeldButton;
static int HeldStep;

static Byte Intensity[LutSize];
static Pixel WhiteCol;
static Pixel BlackCol;
static int Monochrome;

#ifdef THIRTYTWOBIT
static int SwapBytes;
#endif

static int MaxWidth, MaxHeight;
static int MinWidth, MinHeight;
static int MainWide, MainHigh;
static int ScrlX,NewScrlX;
static int ScrlY,NewScrlY;
static int PixDepth;
static int LocalMap;


/* WM_PROTOCOLS */
static char TkInterp[10];
static Atom AppNameAtom;
static Atom DelWinXAtom;
static Atom ProtoXAtom;
static Atom InterpAtom;
static Atom CommAtom;


/*=======================*/
/*  Function Prototypes  */
/*=======================*/
  
extern int ProcessCommand( void );
static int HandleMenuLoop( void );



static void FatalGraphicsError( char *ptr )
{
    char buffer[80];

    sprintf(buffer,"Graphics Error: %s!",ptr);
    RasMolFatalExit(buffer);
}


void AllocateColourMap( void )
{
#ifdef EIGHTBIT
    static XColor Col;
    register int i,j;

    if( Monochrome )
    {   for( i=0; i<LutSize; i++ )
            Intensity[i] = (Byte)((int)(20*RLut[i]+32*GLut[i]+12*BLut[i])>>6);
        return;
    }

    if( LocalMap )
    {   XSetWindowColormap(dpy,MainWin,cmap);
        XSetWindowColormap(dpy,CanvWin,cmap);
        XUninstallColormap(dpy,lmap);
        XFreeColormap(dpy,lmap);
        LocalMap = False;
    } else if( IdentCount )
        XFreeColors(dpy,cmap,Ident,IdentCount,(long)0);
    IdentCount = 0;


    for( i=0; i<LutSize; i++ )
        if( ULut[i] )
        {   Col.red   = RLut[i]<<8 | RLut[i];
            Col.green = GLut[i]<<8 | GLut[i];
            Col.blue  = BLut[i]<<8 | BLut[i];
            Col.flags = DoRed | DoGreen | DoBlue;
            if( !XAllocColor(dpy,cmap,&Col) )
                break;
            Ident[IdentCount++] = Col.pixel;
            Lut[i] = (Pixel)Col.pixel;
        }

    if( i<LutSize )
    {   lmap = XCopyColormapAndFree(dpy,cmap);
        LocalMap = True;

        for( j=0; j<5; j++ )
        {   Col.red   = RLut[j]<<8 | RLut[j];
            Col.green = GLut[j]<<8 | GLut[j];
            Col.blue  = BLut[j]<<8 | BLut[j];
            XAllocColor(dpy,cmap,&Col);
            Lut[i] = (Pixel)Col.pixel;
        }

        for( j=i; j<LutSize; j++ )
            if( ULut[j] )
            {   Col.red   = RLut[j]<<8 | RLut[j];
                Col.green = GLut[j]<<8 | GLut[j];
                Col.blue  = BLut[j]<<8 | BLut[j];
                XAllocColor(dpy,lmap,&Col);
                Lut[j] = (Pixel)Col.pixel;
            }
        XSetWindowColormap(dpy,MainWin,lmap);
        XSetWindowColormap(dpy,CanvWin,lmap);
        XInstallColormap(dpy,lmap);
    }
#else /* EIGHTBIT */
#ifdef THIRTYTWOBIT
    static XColor Col;
    static ByteTest buf;
    register Byte temp;
    register int i;

    for( i=0; i<LutSize; i++ )
        if( ULut[i] )
        {   Col.red   = RLut[i]<<8 | RLut[i];
            Col.green = GLut[i]<<8 | GLut[i];
            Col.blue  = BLut[i]<<8 | BLut[i];
            XAllocColor(dpy,cmap,&Col);
            if( SwapBytes )
            {   buf.longword = (Long)Col.pixel;
                temp = buf.bytes[0];
                buf.bytes[0] = buf.bytes[3];
                buf.bytes[3] = temp;

                temp = buf.bytes[1];
                buf.bytes[1] = buf.bytes[2];
                buf.bytes[2] = temp;
                Lut[i] = buf.longword;
            } else Lut[i] = (Long)Col.pixel;
       }
#else /* THIRTYTWOBIT */
    static XColor Col;
    register int i;

    for( i=0; i<LutSize; i++ )
        if( ULut[i] )
        {   Col.red   = RLut[i]<<8 | RLut[i];
            Col.green = GLut[i]<<8 | GLut[i];
            Col.blue  = BLut[i]<<8 | BLut[i];
            XAllocColor(dpy,cmap,&Col);
            Lut[i] = (Pixel)Col.pixel;
       }
#endif /* THIRTYTWOBIT */
#endif /* EIGHTBIT */
    XSetWindowBackground(dpy,CanvWin,(unsigned long)Lut[5]);
}


static void OpenCanvas( int x, int y )
{
    register unsigned long mask;

    mask = CWEventMask;
    attr.event_mask = ExposureMask | ButtonPressMask | ButtonMotionMask 
                    | ButtonReleaseMask;
    attr.cursor = cross;                           mask |= CWCursor;
    attr.background_pixel = Lut[0];                mask |= CWBackPixel;

    CanvWin = XCreateWindow(dpy, MainWin, 14, MenuHigh+14, x, y, 0, 
                            CopyFromParent, InputOutput, vis, mask, &attr );
}


static void OpenFonts( void )
{
    static char *fontname[] = { "-*-helvetica-bold-o-normal-*-14-*",
                                     "-*-serf-bold-o-normal-*-14-*",
                                        "-*-*-bold-o-normal-*-14-*" };
    register int i;

    cross = XCreateFontCursor(dpy,XC_tcross);
    arrow = XCreateFontCursor(dpy,XC_top_left_arrow);

    for( i=0; i<3; i++ )
        if( (MenuFont=XLoadQueryFont(dpy,fontname[i])) ) 
            break;

    if( !MenuFont )
        FatalGraphicsError("Unable to find suitable font");
    FontHigh = MenuFont->max_bounds.descent +
               MenuFont->max_bounds.ascent + 1;
    MenuHigh = FontHigh+6;
}


static void OpenCursors( void )
{
    Pixmap source,mask;
    XColor black,white;

    white.red = 65535;     black.red = 0;      
    white.green = 65535;   black.green = 0;
    white.blue = 65535;    black.blue = 0;
     
    white.flags = DoRed | DoGreen | DoBlue;
    black.flags = DoRed | DoGreen | DoBlue;

    source = XCreateBitmapFromData(dpy,MainWin,(char*)HGlassData,16,16);
    mask   = XCreateBitmapFromData(dpy,MainWin,(char*)HGlassMask,16,16);
    hglass = XCreatePixmapCursor(dpy,source,mask,&black,&white,7,7);
}


static void OpenColourMap( void )
{
    static XColor Col;
    register int i;

#ifdef EIGHTBIT
    if( !Monochrome )
    {   Col.flags = DoRed | DoGreen | DoBlue;

        for( i=0; i<5; i++ )
        {   Col.red   = RLut[i]<<8 | RLut[i];
            Col.green = GLut[i]<<8 | GLut[i];
            Col.blue  = BLut[i]<<8 | BLut[i];
            if( !XAllocColor(dpy,cmap,&Col) )
            {   cmap = XCopyColormapAndFree(dpy,cmap);
                XAllocColor(dpy,cmap,&Col);
            } 
            Lut[i] = (Pixel)Col.pixel;
        }
        Lut[5] = Lut[0];
    } else /* Black & White */
    {   Lut[0] = BlackCol;
        Lut[1] = BlackCol;
        Lut[2] = WhiteCol;
        Lut[3] = BlackCol;
        Lut[4] = WhiteCol;

        Intensity[5] = 0;
        Lut[5] = 5;         
    }

    LocalMap = False;
    IdentCount = 0;
#else
    Col.flags = DoRed | DoGreen | DoBlue;

    for( i=0; i<5; i++ )
    {   Col.red   = RLut[i]<<8 | RLut[i];
        Col.green = GLut[i]<<8 | GLut[i];
        Col.blue  = BLut[i]<<8 | BLut[i];
        XAllocColor(dpy,cmap,&Col);
        Lut[i] = (Pixel)Col.pixel;
    }
    Lut[5] = Lut[0];
#endif
}


static int RegisterInterpName( char *name )
{
    static unsigned char *registry;
    static unsigned long len,left;
    static char buffer[32];
    static int format;
    static Atom type;

    register int result;
    register char *ptr;

    registry = NULL;
    result = XGetWindowProperty(dpy, RootWindow(dpy,0), InterpAtom,
                                0, 100000, False, XA_STRING, &type,
                                &format, &len, &left, &registry );

    if( (result!=Success) || (format!=8) || (type!=XA_STRING) )
    {   if( (type!=None) && registry ) XFree( (char*)registry );

        sprintf(buffer,"%x %s",(int)MainWin,name);
        XChangeProperty( dpy, RootWindow(dpy,0), InterpAtom, XA_STRING, 
                         8, PropModeReplace, (unsigned char*)buffer, 
                         strlen(buffer)+1 );
        return( True );
    }

    ptr = (char*)registry;
    while( *ptr )
    {   /* Skip Window ID */
        while( *ptr++ != ' ' )
            if( !*ptr ) break;

        /* Compare Interp Name */
        if( !strcmp(ptr,name) )
        {   XFree( (char*)registry );
            return False;
        }

        while( *ptr++ );
    }

    XFree( (char*)registry );
    sprintf(buffer,"%x %s",(int)MainWin,name);
    XChangeProperty( dpy, RootWindow(dpy,0), InterpAtom, XA_STRING, 
                     8, PropModeAppend, (unsigned char*)buffer, 
                     strlen(buffer)+1 );
    return( True );
}


static void DeRegisterInterpName( char *name )
{
    static unsigned char *registry;
    static unsigned long len,left;
    static int format;
    static Atom type;

    register char *src, *dst;
    register int result;

    registry = NULL;
    result = XGetWindowProperty(dpy, RootWindow(dpy,0), InterpAtom,
                                0, 100000, False, XA_STRING, &type,
                                &format, &len, &left, &registry );
    if( type==None )
        return;

    if( (result!=Success) || (format!=8) || (type!=XA_STRING) )
    {   XDeleteProperty( dpy, RootWindow(dpy,0), InterpAtom );
        if( registry ) XFree( (char*)registry );
        return;
    }

    dst = (char*)registry;
    while( *dst )
    {   /* Skip Window ID */
        src = dst;
        while( *src++ != ' ' )
            if( !*src ) break;

        /* Compare Interp Name */
        if( strcmp(src,name) )
        {   while( *dst++ );
        } else break;
    }

    if( *dst )
    {   /* Skip Interp Name */
        while( *src++ );
        
        /* Shuffle Registry */
        while( *src )
            while( (*dst++ = *src++) );
        *dst = 0;

        XChangeProperty( dpy, RootWindow(dpy,0), InterpAtom, XA_STRING,
                         8, PropModeReplace, registry, dst-(char*)registry );
    }
    XFree( (char*)registry );
}


static void OpenIPCComms( void )
{
    auto char buffer[16];
    register int i;

    CommAtom = XInternAtom( dpy, "Comm", False );
    InterpAtom = XInternAtom( dpy, "InterpRegistry", False );
    AppNameAtom = XInternAtom(dpy, "TK_APPLICATION", False );
    DelWinXAtom = XInternAtom(dpy, "WM_DELETE_WINDOW", False);
    /* XSetWMProtocols(dpy,MainWin,&DelWinXAtom,True); */
    if( (ProtoXAtom = XInternAtom(dpy,"WM_PROTOCOLS",False)) )
        XChangeProperty( dpy, MainWin, ProtoXAtom, XA_ATOM, 32, 
                        PropModeReplace, (Byte*)&DelWinXAtom, True );

    i = 0;
    XGrabServer( dpy );
    if( !RegisterInterpName("rasmol") )
    {   strcpy(TkInterp,"rasmol #0");
        for( i=1; i<10; i++ )
        {    TkInterp[8] = i+'0';
             if( RegisterInterpName(TkInterp) )
                 break;
        }

        if( i < 10 ) 
        {   /* Tk4.0 and later! */
            strcpy(buffer,"{rasmol #0}");  buffer[9] = i+'0';
            XChangeProperty( dpy, MainWin, AppNameAtom, XA_STRING, 
                             8, PropModeReplace, (Byte*)buffer, 12 );
        } else *TkInterp = 0;
    } else 
    {   XChangeProperty( dpy, MainWin, AppNameAtom, XA_STRING,
                         8, PropModeReplace, (Byte*)"rasmol", 7 );
        strcpy(TkInterp,"rasmol");
    }
    XUngrabServer( dpy );
}


static void DrawUpBox( Drawable wdw, int x1, int y1, int x2, int y2 )
{
    register int lx,ly,ux,uy;

    lx = x1+1;  ly = y1+1;
    ux = x2-1;  uy = y2-1;

    XSetForeground(dpy,gcon,(unsigned long)Lut[3]);
    XDrawLine(dpy,wdw,gcon,x1,y1,x2,y1);
    XDrawLine(dpy,wdw,gcon,x1,y1,x1,y2);
    XDrawLine(dpy,wdw,gcon,lx,ly,ux,ly);
    XDrawLine(dpy,wdw,gcon,lx,ly,lx,uy);

    XSetForeground(dpy,gcon,(unsigned long)Lut[1]);
    XDrawLine(dpy,wdw,gcon,x2,y1,x2,y2);
    XDrawLine(dpy,wdw,gcon,x1,y2,x2,y2);
    XDrawLine(dpy,wdw,gcon,ux,ly,ux,uy);
    XDrawLine(dpy,wdw,gcon,lx,uy,ux,uy);
}


static void DrawDnBox(  Drawable wdw, int x1, int y1, int x2, int y2 )
{
    register int lx,ly,ux,uy;

    lx = x1+1;  ly = y1+1;
    ux = x2-1;  uy = y2-1;

    XSetForeground(dpy,gcon,(unsigned long)Lut[1]);
    XDrawLine(dpy,wdw,gcon,x1,y1,x2,y1);
    XDrawLine(dpy,wdw,gcon,x1,y1,x1,y2);
    XDrawLine(dpy,wdw,gcon,lx,ly,ux,ly);
    XDrawLine(dpy,wdw,gcon,lx,ly,lx,uy);

    XSetForeground(dpy,gcon,(unsigned long)Lut[3]);
    XDrawLine(dpy,wdw,gcon,x2,y1,x2,y2);
    XDrawLine(dpy,wdw,gcon,x1,y2,x2,y2);
    XDrawLine(dpy,wdw,gcon,ux,ly,ux,uy);
    XDrawLine(dpy,wdw,gcon,lx,uy,ux,uy);
}


static void DrawNoBox( Drawable wdw, int x1, int y1, int x2, int y2 )
{
    register int lx,ly,ux,uy;

    lx = x1+1;  ly = y1+1;
    ux = x2-1;  uy = y2-1;

    XSetForeground(dpy,gcon,(unsigned long)Lut[2]);

    XDrawLine(dpy,wdw,gcon,x1,y1,x2,y1);
    XDrawLine(dpy,wdw,gcon,x2,y1,x2,y2);
    XDrawLine(dpy,wdw,gcon,x2,y2,x1,y2);
    XDrawLine(dpy,wdw,gcon,x1,y2,x1,y1);

    XDrawLine(dpy,wdw,gcon,lx,ly,ux,ly);
    XDrawLine(dpy,wdw,gcon,ux,ly,ux,uy);
    XDrawLine(dpy,wdw,gcon,ux,uy,lx,uy);
    XDrawLine(dpy,wdw,gcon,lx,uy,lx,ly);
}


static void OpenMenuBar( void )
{
    register unsigned long mask;

    mask = CWEventMask;
    attr.event_mask = ExposureMask | ButtonPressMask | ButtonReleaseMask;
    MenuWin = XCreateWindow( dpy, MainWin, 2, 2, XRange+49, FontHigh+5, 0,
                             CopyFromParent, InputOnly, vis, mask, &attr );

    /* Create Unmapped PopUp Window! */
    mask = CWEventMask;
    attr.event_mask = ExposureMask | ButtonPressMask | ButtonReleaseMask | 
                      KeyPressMask;
    attr.background_pixel = Lut[2];     mask |= CWBackPixel;
    attr.border_pixel = Lut[2];         mask |= CWBorderPixel;
    attr.override_redirect = True;      mask |= CWOverrideRedirect;
    attr.save_under = True;             mask |= CWSaveUnder;
    attr.colormap = cmap;               mask |= CWColormap;

    PopUpWin = XCreateWindow(dpy, RootWin, 0, 0, 100, 100, 0, 
                             PixDepth, InputOutput, vis,
                             mask, &attr );
    MenuFocus = False;
    PopUpFlag = False;
}


static void OpenScrollBars( void )
{
    register unsigned long mask;

    Scrl = XCreatePixmap( dpy, MainWin, 16, 16, PixDepth );
    XSetForeground(dpy,gcon,(unsigned long)Lut[2]); 
    XFillRectangle(dpy,Scrl,gcon,0,0,15,15);
    XSetForeground(dpy,gcon,(unsigned long)Lut[0]); 
    XDrawRectangle(dpy,Scrl,gcon,0,0,15,15);
    DrawUpBox( Scrl, 1, 1, 14, 14 );

    tilepix = XCreatePixmapFromBitmapData( dpy, MainWin, (char*)ScrlTile, 8, 8,
                                           (unsigned long)Lut[0], 
                                           (unsigned long)Lut[2], PixDepth );

    mask = CWEventMask;
    attr.event_mask = ExposureMask | ButtonPressMask | ButtonMotionMask 
                    | ButtonReleaseMask;
    attr.background_pixmap = tilepix;              mask |= CWBackPixmap;

    XScrlWin = XCreateWindow(dpy,MainWin,14,YRange+MenuHigh+24,XRange,16, 
                             0,CopyFromParent,InputOutput,vis,mask,&attr);
    lfpix = XCreatePixmapFromBitmapData( dpy, MainWin, (char*)LfArrow, 16, 16,
                                         (unsigned long)Lut[0], 
                                         (unsigned long)Lut[2], PixDepth );
    rgpix = XCreatePixmapFromBitmapData( dpy, MainWin, (char*)RgArrow, 16, 16, 
                                         (unsigned long)Lut[0], 
                                         (unsigned long)Lut[2], PixDepth );

    YScrlWin = XCreateWindow(dpy,MainWin,XRange+24,MenuHigh+14,16,YRange, 
                             0,CopyFromParent,InputOutput,vis,mask,&attr);
    uppix = XCreatePixmapFromBitmapData( dpy, MainWin, (char*)UpArrow, 16, 16,
                                         (unsigned long)Lut[0], 
                                         (unsigned long)Lut[2], PixDepth );
    dnpix = XCreatePixmapFromBitmapData( dpy, MainWin, (char*)DnArrow, 16, 16,
                                         (unsigned long)Lut[0], 
                                         (unsigned long)Lut[2], PixDepth );

    ScrlX = (XRange/2)-8;
    ScrlY = (YRange/2)-8;
}


static void DrawXScroll( void )
{
    XCopyArea(dpy,rgpix,XScrlWin,gcon,0,0,16,16,XRange-16,0);
    XCopyArea(dpy,Scrl ,XScrlWin,gcon,0,0,16,16,ScrlX,0);
    XCopyArea(dpy,lfpix,XScrlWin,gcon,0,0,16,16,0,0);
}


static void DrawYScroll( void )
{
    XCopyArea(dpy,dnpix,YScrlWin,gcon,0,0,16,16,0,YRange-16);
    XCopyArea(dpy,Scrl ,YScrlWin,gcon,0,0,16,16,0,ScrlY);
    XCopyArea(dpy,uppix,YScrlWin,gcon,0,0,16,16,0,0);
}


void UpdateScrollBars( void )
{
    register int temp;

    temp = (DialValue[YScrlDial]+1.0)*(YRange-48);  
    NewScrlY = (temp>>1)+16;

    if( NewScrlY != ScrlY )
    {   XClearArea(dpy,YScrlWin,0,ScrlY,16,16,False);
        XCopyArea(dpy,Scrl,YScrlWin,gcon,0,0,16,16,0,NewScrlY);
        ReDrawFlag |= (1<<YScrlDial);
        ScrlY = NewScrlY; 
    }

    temp = (DialValue[XScrlDial]+1.0)*(XRange-48);  
    NewScrlX = (temp>>1)+16;

    if( NewScrlX != ScrlX )
    {   XClearArea(dpy,XScrlWin,ScrlX,0,16,16,False);
        XCopyArea(dpy,Scrl,XScrlWin,gcon,0,0,16,16,NewScrlX,0);
        ReDrawFlag |= (1<<XScrlDial);
        ScrlX = NewScrlX;
    }
    XFlush(dpy);
}



#ifdef DIALBOX
static void SetDialLabel( int num, char *ptr )
{
    static XStringFeedbackControl ctrl;
    static KeySym text[8];
    register int length;

    length = 0;
    while( *ptr )
       text[length++] = *ptr++;

    ctrl.id = num;
    ctrl.num_keysyms = length;
    ctrl.class = ValuatorClass;
    ctrl.syms_to_display = text;
    XChangeFeedbackControl(dpy,Dials,DvString,
                           (XFeedbackControl*)&ctrl);
}


static void GetDialState( void )
{
    register XValuatorState *ptr;
    register XDeviceState *stat;
    register int i,j,max;

    stat = XQueryDeviceState(dpy,Dials);
    ptr = (XValuatorState*)stat->data;
    for( i=0; i<stat->num_classes; i++ )
    {   if( ptr->class == ValuatorClass )
        {   if( ptr->mode & 0x01 )
            {   DialMode = Absolute;
                max = MinFun(ptr->num_valuators,8);
                for( j=0; j<max; j++ )
                    DialPrev[j] = ptr->valuators[j];
            } else DialMode = Relative;
            break;
        } else ptr = (XValuatorState*)(((char*)ptr) + 
                                       ptr->length);
    }
    XFreeDeviceState(stat);
}


static void OpenDialsBox( void )
{
    register XValuatorInfo *valptr;
    register XFeedbackState *list;
    register XFeedbackState *feed;
    register XDeviceInfo *devlist;
    register XDeviceInfo *ptr;
    register Atom devtype;
    register int i,j,max;

    static XEventClass dclass;
    static int count;

    UseDials = False;
    /* Avoid X Server's without the extension */
    if( !XQueryExtension(dpy,"XInputExtension",
                         &count,&count,&count) )
        return;
    
    devlist = XListInputDevices(dpy,&count);
    devtype = XInternAtom(dpy,XI_KNOB_BOX,True );
    if( (devtype==None) || !devlist ) return;

    ptr = devlist;
    for( i=0; i<count; i++ )
        if( (ptr->use==IsXExtensionDevice) && (ptr->type==devtype) )
        {   valptr = (XValuatorInfo*)ptr->inputclassinfo;
            for( j=0; j<ptr->num_classes; j++ )
            {   if( valptr->class == ValuatorClass )
                    if( (Dials=XOpenDevice(dpy,ptr->id)) )
                    {   UseDials = True;
                        break;
                    }
                valptr = (XValuatorInfo*)(((char*)valptr) +
                                          valptr->length);
            }
            if( UseDials ) break;
        } else ptr++;
    /* XFreeDeviceList(devlist); */

    if( UseDials ) 
    {   /* Determine Dial Mapping! */
        if( !strcmp(ServerVendor(dpy),"Silicon Graphics") )
        {      DialMap = SGIDialMap;
        } else DialMap = ESVDialMap;

        DialMode = valptr->mode;
        max = MinFun(valptr->num_axes,8);
        for( i=0; i<max; i++ )
            DialRes[i] = (Real)valptr->axes[i].resolution;
        GetDialState();
    } else return;

    UseDialLEDs = 0;
    feed = list = XGetFeedbackControl( dpy, Dials, &count );
    for( i=0; i<count; i++ )
    {   if( feed->class == StringFeedbackClass ) UseDialLEDs++;
        feed = (XFeedbackState*)(((char*)feed) + feed->length);
    }
    XFreeFeedbackList( list );

    if( UseDialLEDs >= 8 )
    {   for( i=0; i<8; i++ )
            SetDialLabel(i,DialLabel[DialMap[i]]);
    } else UseDialLEDs = False;

    DeviceMotionNotify( Dials, DialEvent, dclass );
    XSelectExtensionEvent( dpy, MainWin, &dclass, 1 );
    XSelectExtensionEvent( dpy, MenuWin, &dclass, 1 );
    XSelectExtensionEvent( dpy, CanvWin, &dclass, 1 );
    XSelectExtensionEvent( dpy, XScrlWin, &dclass, 1 );
    XSelectExtensionEvent( dpy, YScrlWin, &dclass, 1 );
}


static void HandleDialEvent( XDeviceMotionEvent *ptr )
{
    register double temp;
    register int count;
    register int value;
    register int index;
    register int num;

    /* Limit Number of Dials */
    count = 8 - ptr->first_axis;
    if( count > (int)ptr->axes_count )
        count = (int)ptr->axes_count;

    for( index=0; index<count; index++ )
    {   num = ptr->first_axis+index;
        if( DialMode == Absolute )
        {   value = ptr->axis_data[index] - DialPrev[num];
            DialPrev[num] = ptr->axis_data[index];
        } else value = ptr->axis_data[index];

        if( value )
        {   temp = (Real)value/DialRes[num];
            num = DialMap[num];
            temp += DialValue[num];
            ReDrawFlag |= (1<<num);

            if( num<3 )
            {   while( temp<-1.0 ) temp += 2.0;
                while( temp>1.0 )  temp -= 2.0;
            } else
            {   if( temp<-1.0 ) temp = -1.0;
                if( temp>1.0 )  temp = 1.0;
            }
            DialValue[num] = temp;

            if( num==YScrlDial )
            {   value = (temp+1.0)*(YRange-48);
                NewScrlY = (value>>1)+16;
            }

            if( num==XScrlDial )
            {   value = (temp+1.0)*(XRange-48);
                NewScrlX = (value>>1)+16;
            }
        }
    }
}
#endif


static void DrawMainWin( void )
{
    register int temp;

    DrawUpBox(MainWin,0,0,MainWide,MainHigh);
    DrawUpBox(MainWin,0,0,MainWide-2,FontHigh+7);

    temp = YRange+MenuHigh;
    DrawDnBox(MainWin,12,MenuHigh+12,XRange+16,temp+16);
    DrawDnBox(MainWin,XRange+22,MenuHigh+12,XRange+41,temp+16);
    DrawDnBox(MainWin,12,temp+22,XRange+16,temp+41);
}


/*====================*/
/*  Menu Bar Display  */
/*====================*/

static void DisplayMenuBarText( BarItem *ptr, int x, int y )
{
    register unsigned long col;
    register int under, pos, i, index, wide;

    if( ptr->flags&mbEnable && !DisableMenu )
    {      col = Lut[0];
    } else col = Lut[1];
    XSetForeground( dpy, gcon, col );

    XDrawString( dpy, MainWin, gcon, x, y, *(ptr->text), *(ptr->len) );

    under = y + MenuFont->descent;

    pos = x;
    for( i=0; i<*(ptr->pos); i++)
    {   index = (*(ptr->text))[i] - MenuFont->min_char_or_byte2;
        pos += MenuFont->per_char[index].width;
    }

    index = (*(ptr->text))[*(ptr->pos)] - MenuFont->min_char_or_byte2;
    wide = pos+MenuFont->per_char[index].rbearing;
    pos += MenuFont->per_char[index].lbearing;

    XDrawLine( dpy, MainWin, gcon, pos, under, wide, under );
}


static void DrawMenuBar( void )
{
    register BarItem *ptr;
    register int wide;
    register int x,y;
    register int i;

    x = 6; y = MenuFont->ascent+4;
    XSetFont( dpy, gcon, MenuFont->fid );

    for( i=0; i<MenuBarMax; i++ )
    {   ptr = MenuBar+i;
        wide = XTextWidth( MenuFont, *(ptr->text), *(ptr->len) );
        if( x+wide+24 > MainWide ) break;

        /* Right Justify "Help" */
        if( i == MenuBarMax-1 )
            x = MainWide - (wide+24);

        DisplayMenuBarText( ptr, x+8, y );

        if( MenuFocus && (i==MenuBarSelect) )
        {      DrawUpBox( MainWin, x, 2, x+wide+16, FontHigh+5 );
        } else DrawNoBox( MainWin, x, 2, x+wide+16, FontHigh+5 );
        x += wide+24;
    }
    MenuBarCount = i;
    /* XSync(dpy,False); */
    XFlush(dpy);
}


/*=======================*/
/*  Pop-up Menu Display  */
/*=======================*/

static void DisplayPopUpText( MenuItem *ptr, int x, int y )
{
    register unsigned long col;
    register int pos, wide;
    register int i,under;
    register int index;

    col = (ptr->flags&mbEnable)? Lut[0] : Lut[1];
    XSetForeground( dpy, gcon, col );

    XDrawString( dpy, PopUpWin, gcon, x, y, *(ptr->text), *(ptr->len) );

    if( ptr->flags & mbAccel )
    {   under = y + MenuFont->descent;

        pos = x;
        for( i=0; i<*(ptr->pos); i++ )
        {   index = (*(ptr->text))[i] - MenuFont->min_char_or_byte2;
            pos += MenuFont->per_char[index].width;
        }

        index = (*(ptr->text))[*(ptr->pos)] - MenuFont->min_char_or_byte2;
        wide = pos+MenuFont->per_char[index].rbearing;
        pos += MenuFont->per_char[index].lbearing;

        XDrawLine( dpy, PopUpWin, gcon, pos, under, wide, under );
    }
}


static void DrawPopUpMenu( void )
{
    register MenuItem *ptr;
    register int count;
    register int x,y;
    register int i;

    DrawUpBox(PopUpWin,0,0,PopUpWide,PopUpHigh);

    ptr = MenuBar[MenuBarSelect].menu;
    count = MenuBar[MenuBarSelect].count;

    y = 2;  x = 2;
    for( i=0; i<count; i++ )
    {   if( !(ptr->flags&mbSepBar) )
        {   DisplayPopUpText( ptr, x+8, y+MenuFont->ascent+2 );

            if( ItemFlag && (i==MenuItemSelect) )
            {      DrawUpBox(PopUpWin,2,y,PopUpWide-2,y+FontHigh+3);
            } else DrawNoBox(PopUpWin,2,y,PopUpWide-2,y+FontHigh+3);
            y += FontHigh+4;
        } else
        {   XSetForeground( dpy, gcon, (unsigned long)Lut[1] );
            XDrawLine(dpy,PopUpWin,gcon,2,y,PopUpWide-2,y);
            XSetForeground( dpy, gcon, (unsigned long)Lut[3] );
            XDrawLine(dpy,PopUpWin,gcon,2,y+1,PopUpWide-2,y+1);
            y += 2;
        }
        ptr++;
    }
    /* XSync(dpy,False); */
    XFlush(dpy);
}


static void DisplayPopUpMenu( int i, int x )
{
    register int wide, count;
    register MenuItem *ptr;
    static int xpos, ypos;
    static Window win;


    MenuBarSelect = i;
    DrawMenuBar();

    ptr = MenuBar[i].menu;
    count = MenuBar[i].count;

    PopUpHigh = 4;
    PopUpWide = 4;
    for( i=0; i<count; i++ )
    {   if( !(ptr->flags&mbSepBar) )
        {   wide = XTextWidth(MenuFont,*(ptr->text),*(ptr->len));
            if( wide+28 > PopUpWide ) PopUpWide = wide+28;
            PopUpHigh += FontHigh+4;
        } else PopUpHigh += 2;
        ptr++;
    }

    /* Determine pop-up menu position! */
    XTranslateCoordinates(dpy,MainWin,RootWin,x,FontHigh+6,
                          &xpos, &ypos, &win );

    if( ypos+PopUpHigh > MaxHeight )
        ypos -= (PopUpHigh+FontHigh+6);
    if( xpos+PopUpWide > MaxWidth )
        xpos = MaxWidth-PopUpWide;
    if( xpos < 0 ) xpos = 0;

    XUnmapWindow(dpy,PopUpWin);
    XMoveResizeWindow(dpy,PopUpWin,xpos,ypos,PopUpWide+1,PopUpHigh+1);
    XRaiseWindow(dpy,PopUpWin);
    XMapWindow(dpy,PopUpWin);
    PopUpFlag = True;
    DrawPopUpMenu();
}



/*==============================*/
/*  Pop-Up Menu Event Handling  */
/*==============================*/

static void HandleItemClick( int x, int y )
{
    register MenuItem *ptr;
    register int count,i;

    static int xpos, ypos;
    static Window win;

    XTranslateCoordinates(dpy,MenuWin,PopUpWin,x,y,
                          &xpos,&ypos,&win);

    /* Ignore by not setting ItemFocus! */
    if( (xpos<0) || (xpos>PopUpWide) ) return;
    if( (ypos<0) || (ypos>PopUpHigh) ) return;
    ItemFocus = True;

    ptr = MenuBar[MenuBarSelect].menu;
    count = MenuBar[MenuBarSelect].count;

    y = 2;
    for( i=0; i<count; i++ )
    {   if( !(ptr->flags&mbSepBar) )
        {   if( (ypos>=y) && (ypos<=y+FontHigh+3) )
            {   if( ptr->flags & mbEnable )
                {   if( !ItemFlag || (MenuItemSelect!=i) )
                    {   /* Avoid Flickering */
                        MenuItemSelect = i;
                        ItemFlag = True;
                        DrawPopUpMenu();
                    }
                    return;
                } else break;
            }
            y += FontHigh+4;
        } else y += 2;
        ptr++;
    }

    if( ItemFlag )
    {   ItemFlag = False;
        DrawPopUpMenu();
    }
}


static void HandleItemMove( int x, int y )
{
    register MenuItem *ptr;
    register int count,i;

    static int xpos, ypos;
    static Window win;

    XTranslateCoordinates(dpy,MenuWin,PopUpWin,x,y,
                          &xpos,&ypos,&win);

    if( (xpos>=0) && (xpos<=PopUpWide) )
    {   ptr = MenuBar[MenuBarSelect].menu;
        count = MenuBar[MenuBarSelect].count;

        y = 2;
        for( i=0; i<count; i++ )
        {   if( !(ptr->flags&mbSepBar) )
            {   if( (ypos>=y) && (ypos<=y+FontHigh+3) )
                {   if( !ItemFlag || (MenuItemSelect!=i) )
                    {   /* Avoid Flicker! */
                        MenuItemSelect = i;
                        ItemFlag = True;
                        DrawPopUpMenu();
                    }
                    ItemFocus = True;
                    return;
                }
                y += FontHigh+4;
            } else y += 2;
            ptr++;
        }
    }

    if( ItemFlag )
    {   /* Avoid Flicker! */
        ItemFlag = False;
        DrawPopUpMenu();
    }
}


static int HandleItemKey( int key )
{
    register MenuItem *ptr;
    register int count;
    register int item;
    register int ch;
    register int i;

    key = ToUpper( key );
    item = MenuItemSelect;
    ptr = &MenuBar[MenuBarSelect].menu[item];
    count = MenuBar[MenuBarSelect].count;
    for( i=0; i<count; i++ )
    {   if( (ptr->flags&(mbEnable|mbAccel)) && 
           !(ptr->flags&mbSepBar) )
        {   ch = (*(ptr->text))[*(ptr->pos)];
            if( ToUpper(ch) == key )
                return( (MenuBarSelect<<8)+item+1 );
        }

        /* Advance to next item! */
        if( item == count-1 )
        {   ptr = MenuBar[MenuBarSelect].menu;
            item = 0;
        } else 
        {   item++;
            ptr++;
        }
    }
    return 0;
}


static void SelectFirstItem( int menu )
{
    register MenuItem *ptr;
    register int count;
    register int i;

    count = MenuBar[menu].count;
    ptr = MenuBar[menu].menu;

    ItemFlag = False;
    for( i=0; i<count; i++ )
        if( (ptr->flags&mbEnable) &&
           !(ptr->flags&mbSepBar) )
        {   MenuItemSelect = i;
            ItemFlag = True;
            break;
        } else ptr++;
}


static void SelectPrevItem( void )
{
    register BarItem *ptr;
    register int flags;
    register int item;
    register int i;

    if( !ItemFlag )
        return;

    item = MenuItemSelect;
    ptr = MenuBar + MenuBarSelect;
    for( i=0; i<ptr->count; i++ )
    {   if( !item )
        {   item = ptr->count-1;
        } else item--;

        flags = ptr->menu[item].flags;
        if( (flags&mbEnable) && !(flags&mbSepBar) )
            break;
    }

    if( item != MenuItemSelect )
    {   MenuItemSelect = item;
        DrawPopUpMenu();
    }
}


static void SelectNextItem( void )
{
    register BarItem *ptr;
    register int flags;
    register int item;
    register int i;

    if( !ItemFlag )
        return;

    item = MenuItemSelect;
    ptr = MenuBar + MenuBarSelect;
    for( i=0; i<ptr->count; i++ )
    {   if( item == ptr->count-1 )
        {   item = 0;
        } else item++;

        flags = ptr->menu[item].flags;
        if( (flags&mbEnable) && !(flags&mbSepBar) )
            break;
    }

    if( item != MenuItemSelect )
    {   MenuItemSelect = item;
        DrawPopUpMenu();
    }
}



/*===========================*/
/*  Menu Bar Event Handling  */
/*===========================*/

static void SelectMenu( int menu )
{
    register BarItem *ptr;
    register int wide;
    register int i,x;

    if( !PopUpFlag )
    {   MenuBarSelect = menu;
        DrawMenuBar();
        return;
    }

    if( menu != MenuBarMax-1 )
    {   x = 6;
        for( i=0; i<menu; i++ )
        {   ptr = MenuBar+i;
            wide = XTextWidth(MenuFont,*(ptr->text),*(ptr->len));
            x += wide+24;
        }
    } else 
    {   ptr = MenuBar+menu;
        wide = XTextWidth(MenuFont,*(ptr->text),*(ptr->len));
        x = MainWide - (wide+24);
    }

    SelectFirstItem( menu );
    DisplayPopUpMenu( menu, x );
    ItemFocus = False;
}


static int HandleMenuClick( int pos )
{
    register BarItem *ptr;
    register int wide;
    register int x,i;

    x = 6;
    for( i=0; i<MenuBarCount; i++ )
    {   ptr = MenuBar+i;
        wide = XTextWidth( MenuFont, *(ptr->text), *(ptr->len) );
        if( i == MenuBarMax-1 ) x = MainWide - (wide+24);

        if( (pos>=x) && (pos<=x+wide+16) )
        {   if( !PopUpFlag || (MenuBarSelect!=i) )
            {   ItemFlag = False;
                DisplayPopUpMenu(i,x);
            } else if( ItemFlag )
            {   ItemFlag = False;
                DrawPopUpMenu();
            }
            ItemFocus = True;
            return True;
        } else x += wide+24;
    }
    return False;
}


static int HandleMenuKey( char key )
{
    register int i;

    key = ToUpper(key);
    for( i=0; i<MenuBarCount; i++ )
        if( ToUpper((*(MenuBar[i].text))[*(MenuBar[i].pos)]) == key )
        {   if( !PopUpFlag || (MenuBarSelect!=i) )
            {   PopUpFlag = True;
                SelectMenu( i );
            }
            return True;
        }
    return False;
}


void EnableMenus( int flag )
{
    DisableMenu = !flag;
    if( Interactive )
        DrawMenuBar();
}


static void ReSizeWindow( int wide, int high )
{
    register Real xpos;
    register Real ypos;
    register int dx;

    xpos = (XRange>48)? (Real)(ScrlX-16)/(XRange-48) : 0.0;
    ypos = (YRange>48)? (Real)(ScrlY-16)/(YRange-48) : 0.0;

    YRange = high-(MenuHigh+53);
    XRange = wide-53;

    if( (dx = XRange%4) )
        XRange += 4-dx;

    MainHigh = YRange+(MenuHigh+53);  HRange = YRange>>1;
    MainWide = XRange+53;             WRange = XRange>>1;
    Range = MinFun(XRange,YRange);

    XResizeWindow( dpy, CanvWin, XRange, YRange);
    XResizeWindow( dpy, MenuWin, XRange+49, FontHigh+5 );
    XMoveResizeWindow( dpy, XScrlWin, 14, YRange+MenuHigh+24, XRange, 16 );
    XMoveResizeWindow( dpy, YScrlWin, XRange+24, MenuHigh+14, 16, YRange );

    NewScrlX = ScrlX = (xpos*(XRange-48))+16;
    NewScrlY = ScrlY = (ypos*(YRange-48))+16;

    XClearWindow( dpy, MainWin );
    XClearWindow( dpy, CanvWin );

    DrawXScroll();
    DrawYScroll();
    DrawMainWin();
    DrawMenuBar();

    ReDrawFlag |= RFReSize;
    XSync(dpy,True);
}

void ReDrawWindow( void )
{
    if( Interactive )
        ReSizeWindow( MainWide, MainHigh );
}



int FatalXError( Display *ptr )
{
     /* Avoid Compiler Warnings! */
    UnusedArgument(ptr);

    dpy = (Display*)NULL;
    RasMolFatalExit("*** Fatal X11 I/O Error! ***");
    return 0;
}


int OpenDisplay( int x, int y )
{
#ifdef THIRTYTWOBIT
    static ByteTest test;
#endif
    register unsigned long mask;
    register int i,num;
    register char *ptr;
 
    static XVisualInfo visinfo;
    static XClassHint xclass;
    static XSizeHints size;
    static Pixmap icon;
    static int temp;
    static char VersionStr[50];

    sprintf (VersionStr,"RasMol Version %s", VERSION);

    image = (XImage*)NULL;

    MouseCaptureStatus = False;
    MouseUpdateStatus = False;
    UseHourGlass = True;
    DisableMenu = False;
    Monochrome = False;
    HeldButton = -1;

    for( i=0; i<8; i++ )
         DialValue[i] = 0.0;

    RLut[0]=0;   GLut[0]=0;   BLut[0]=0;    ULut[0]=True;
    RLut[1]=100; GLut[1]=100; BLut[1]=100;  ULut[1]=True;
    RLut[2]=150; GLut[2]=150; BLut[2]=150;  ULut[2]=True;
    RLut[3]=200; GLut[3]=200; BLut[3]=200;  ULut[3]=True;
    RLut[4]=255; GLut[4]=255; BLut[4]=255;  ULut[4]=True;

    XRange = x;  WRange = XRange>>1;
    YRange = y;  HRange = YRange>>1;
    Range = MinFun(XRange,YRange);

    if( !Interactive ) return( False );
    if( (dpy=XOpenDisplay(NULL)) == NULL )
        return 0;

    num = DefaultScreen(dpy);
    RootWin = RootWindow(dpy,num);
    XSetIOErrorHandler( FatalXError );

#ifdef EIGHTBIT
    if( !(XMatchVisualInfo(dpy,num,8,PseudoColor,&visinfo) ||
          XMatchVisualInfo(dpy,num,8,GrayScale,&visinfo)) )
    {   /* Attempt to use Monochrome Mode! */
        if( !(XMatchVisualInfo(dpy,num,1,StaticColor,&visinfo) ||
              XMatchVisualInfo(dpy,num,1,StaticGray,&visinfo)) )
        {   XCloseDisplay(dpy);
            return 0;
        }
        Monochrome = True;
        PixDepth = 1;
    } else PixDepth = 8;
#else
#ifdef THIRTYTWOBIT
    if( XMatchVisualInfo(dpy,num,32,TrueColor,&visinfo) ||
        XMatchVisualInfo(dpy,num,32,DirectColor,&visinfo) )
    {   PixDepth = 32;
    } else if( XMatchVisualInfo(dpy,num,24,TrueColor,&visinfo) ||
               XMatchVisualInfo(dpy,num,24,DirectColor,&visinfo) )
    {   PixDepth = 24;
    } else /* No suitable display! */
    {   XCloseDisplay(dpy);
        return(0);
    }
#else /* SIXTEENBIT */
    if( XMatchVisualInfo(dpy,num,16,TrueColor,&visinfo) ||
        XMatchVisualInfo(dpy,num,16,DirectColor,&visinfo) )
    {   PixDepth = 16;
    } else if( XMatchVisualInfo(dpy,num,15,TrueColor,&visinfo) ||
               XMatchVisualInfo(dpy,num,15,DirectColor,&visinfo) )
    {   PixDepth = 15;
    } else /* No suitable display! */
    {   XCloseDisplay(dpy);
        return 0;
    }
#endif
#endif

    if( !Monochrome )
    {   vis = visinfo.visual;
        if( vis != DefaultVisual(dpy,num) )
        {   cmap = XCreateColormap(dpy,RootWin,vis,AllocNone);
        } else cmap = DefaultColormap(dpy,num);
    } else /* Black & White */
    {   /* PixDepth = DefaultDepth(dpy,num); */
        cmap = DefaultColormap(dpy,num);
        vis = visinfo.visual;

        BlackCol = (Pixel)(BlackPixel(dpy,num)&1);
        WhiteCol = (Pixel)(WhitePixel(dpy,num)&1);
    }

    OpenFonts();
    OpenColourMap();

    MaxHeight = DisplayHeight(dpy,num);  MinHeight = MenuHigh+101;
    MaxWidth = DisplayWidth(dpy,num);    MinWidth = 101;

    MainHigh = YRange+MenuHigh+53;
    MainWide = XRange+53;

    mask = CWEventMask;
    attr.event_mask = ExposureMask | KeyPressMask | StructureNotifyMask
                    | EnterWindowMask | LeaveWindowMask | PropertyChangeMask;
    attr.background_pixel = Lut[2];     mask |= CWBackPixel;
    attr.border_pixel = Lut[2];         mask |= CWBorderPixel;
    attr.colormap = cmap;               mask |= CWColormap;
    attr.cursor = arrow;                mask |= CWCursor;

    MainWin = XCreateWindow(dpy, RootWin, 0, 0, MainWide, MainHigh, 2,
			    PixDepth, InputOutput, vis, mask, &attr );

    gcon = XCreateGC(dpy,MainWin,0L,NULL);
    /* DefaultGC(dpy,num) */

    XSetGraphicsExposures(dpy,gcon,False);
    icon = XCreateBitmapFromData(dpy,MainWin,(char*)icon_bits,
                                 icon_width,icon_height );

    size.flags = PMinSize | PMaxSize;
    size.min_width = MinWidth;    size.max_width = MaxWidth;
    size.min_height = MinHeight;  size.max_height = MaxHeight;
    XSetStandardProperties(dpy, MainWin, VersionStr,
                           "RasMol", icon, NULL, 0, &size );

    xclass.res_name = "rasmol";
    xclass.res_class = "RasMol";
    XSetClassHint(dpy,MainWin,&xclass);

    hints.icon_pixmap = icon;       
    hints.flags = IconPixmapHint;
    XSetWMHints(dpy,MainWin,&hints);

    OpenCanvas( XRange, YRange );
    OpenScrollBars();
    OpenMenuBar();
    OpenCursors();
    OpenIPCComms();

#ifdef DIALBOX
    OpenDialsBox();
#endif

#ifdef MITSHM
    ptr = DisplayString(dpy);
    if( !ptr || (*ptr==':') || !strncmp(ptr,"localhost:",10) || 
        !strncmp(ptr,"unix:",5) || !strncmp(ptr,"local:",6) )
    {   SharedMemOption = XQueryExtension(dpy,"MIT-SHM",&temp,&temp,&temp);
        if( Monochrome && (PixDepth!=1) ) SharedMemOption = False;
    } else SharedMemOption = False;
    SharedMemFlag = False;
#endif

#ifdef THIRTYTWOBIT
    /* Determine Byte Ordering */
    test.longword = (Long)0x000000ff;
    if( ImageByteOrder(dpy) == MSBFirst )
    {      SwapBytes = test.bytes[0];
    } else SwapBytes = test.bytes[3];
#endif

    XMapSubwindows(dpy,MainWin);
    XMapWindow(dpy,MainWin);

    DrawXScroll();
    DrawYScroll();
    DrawMainWin();
    DrawMenuBar();

    XClearWindow( dpy, CanvWin );
    XSync(dpy,False);

    num = 1<<ConnectionNumber(dpy);
    return( num );
}


int CreateImage( void )
{
    register Long size, temp;
    register int format;
    register Pixel *ptr;

    if( !Interactive )
    {   if( FBuffer ) free(FBuffer);
        size = (Long)XRange*YRange*sizeof(Pixel);
        FBuffer = (Pixel*)malloc( size+32 );
        return( (int)FBuffer );
    }

    format = Monochrome? XYPixmap : ZPixmap;

    if( image ) 
    {   /* Monochrome Mode Frame Buffer! */
        if( FBuffer && (FBuffer!=(Pixel*)image->data) )
            free(FBuffer);
#ifdef MITSHM
        if( SharedMemFlag )
        {   XShmDetach( dpy, &xshminfo );
            image->data = (char*)NULL;
            shmdt( xshminfo.shmaddr );
        }
#endif
        XDestroyImage( image );
        image = (XImage*)NULL;
    }

    if( Monochrome )
    {   /* Monochrome Mode Frame Buffer! */
        size = (Long)XRange*YRange*sizeof(Pixel);
        FBuffer = (Pixel*)malloc( size+32 );
        if( !FBuffer ) return False;

        /* Bit per Pixel ScanLines! */
        temp = ((XRange+31)>>5)<<2;
        size = (Long)temp*YRange + 32;
    } else 
        size = (Long)XRange*YRange*sizeof(Pixel) + 32;

#ifdef MITSHM
    if( SharedMemOption )
    {   SharedMemFlag = False;
        image = XShmCreateImage( dpy, vis, PixDepth, format,
                                 NULL, &xshminfo, XRange, YRange );

        if( image )
        {   temp = (Long)image->bytes_per_line * image->height;
            if( temp > size ) size = temp;
            xshminfo.shmid = shmget( IPC_PRIVATE, size, IPC_CREAT|0777 );
            if( xshminfo.shmid != -1 ) 
            {   xshminfo.shmaddr = (char*)shmat(xshminfo.shmid,0,0);
                if( xshminfo.shmaddr != (char*)-1 )
                {   image->data = xshminfo.shmaddr;
                    if( !Monochrome )
                        FBuffer = (Pixel*)image->data;
                    xshminfo.readOnly = True;

                    SharedMemFlag = XShmAttach( dpy, &xshminfo );
                    XSync(dpy,False);
                }
                /* Always Destroy Shared Memory Ident */
                shmctl( xshminfo.shmid, IPC_RMID, 0 );
            }

            if( SharedMemFlag )
            {   if( Monochrome )
                {   if( BlackCol )
                    {      memset((void*)image->data,255,size);
                    } else memset((void*)image->data,255,size);
                }
                return True;
            } else 
            {   XDestroyImage( image );
                image = (XImage*)NULL;
            }
        }
    }
#endif

    /* Allocate Frame Buffer! */
    ptr = (Pixel*)malloc( size );
    if( !ptr ) return False;

    if( !Monochrome ) FBuffer = ptr;
    image = XCreateImage( dpy, vis, PixDepth, format, 0, (char*)ptr, 
                          XRange, YRange, sizeof(Pixel)<<3, 0 );
    return (int)image;
}


static void DitherImage( void )
{
    register Card bits;
    register Card *dst;
    register Pixel *src;
    register int xmax,ymax;
    register int count,x,y;
    register int error;

    register Card bmask,wmask;
    register Card bhigh,whigh;
    register int wlen;
    register int blen;
    register int len;

    src = (Pixel*)FBuffer;
    dst = (Card*)image->data;

    wlen = XRange>>5;
    blen = XRange&31;
    if( blen )
    {   wmask = WhiteCol << (blen-1);
        bmask = BlackCol << (blen-1);
        len = wlen+1;
    } else len = wlen;

    whigh = WhiteCol << 31;
    bhigh = BlackCol << 31;

    /* Allow Compiler Optimisation */
    xmax = XRange;  ymax = YRange;

    error = 0;
    for( y=0; y<ymax; y++ )
    {    for( x=0; x<wlen; x++ )
         {   for( count=0; count<32; count++ )
             {   error += Intensity[*src++];
                 bits <<= 1;
                 if( error >= 128  )
                 {   error -= 255;
                        bits |= WhiteCol;
                 } else bits |= BlackCol;
             }
             *dst++ = bits;
         }

         if( blen )
         {   for( count=0; count<blen; count++ )
             {   error += Intensity[*src++];
                 bits <<= 1;
                 if( error >= 128  )
                 {   error -= 255;
                        bits |= WhiteCol;
                 } else bits |= BlackCol;
             }
             *dst++ = bits;
           
             /* Asymmetric Loop Unrolling! */
             if( ++y == ymax ) break;
             src += xmax;
             dst += wlen;

             bits = 0;
             for( count=0; count<blen; count++ )
             {   error += Intensity[*(--src)];
                 bits >>= 1;
                 if( error >= 128  )
                 {   error -= 255;
                        bits |= wmask;
                 } else bits |= bmask;
             }
             *(--dst) = bits;
         } else
         {   /* Asymmetric Loop Unrolling! */
             if( ++y == ymax ) break;
             src += xmax;
             dst += len;
         }

         for( x=wlen-1; x>=0; x-- )
         {   for( count=0; count<32; count++ )
             {   error += Intensity[*(--src)];
                 bits >>= 1;
                 if( error >= 128  )
                 {   error -= 255;
                        bits |= whigh;
                 } else bits |= bhigh;
             }
             *(--dst) = bits;
         }
         src += xmax;
         dst += len;
    }
}


void TransferImage( void )
{
    if( Monochrome )
        DitherImage();

#ifdef MITSHM
    if( SharedMemFlag )
    {   XShmPutImage(dpy,CanvWin,gcon,image,0,0,0,0,XRange,YRange,False);
        XSync(dpy,False);
    } else
    {   XPutImage( dpy, CanvWin, gcon,image,0,0,0,0,XRange,YRange);
        XFlush(dpy);
    }
#else
    XPutImage( dpy, CanvWin, gcon, image, 0, 0, 0, 0, XRange, YRange );
    XFlush(dpy);
#endif
}


void ClearImage( void )
{
    XClearWindow( dpy, CanvWin );
    XFlush(dpy);
}


int PrintImage( void )
{
    return False;
}


int ClipboardImage( void )
{
    return False;
}


void SetCanvasTitle( char *ptr )
{
    if( Interactive ) {
      XStoreName(dpy,MainWin,ptr);
    }
}



static int HandleIPCError( Display *dpy, XErrorEvent *ptr )
{
    /* Avoid Compiler Warnings! */
    UnusedArgument(dpy);
    UnusedArgument(ptr);
    return 0;
}


static void HandleIPCCommand( void )
{
    static unsigned long len,left;
    static unsigned char *command;
    static Window source;
    static int serial;
    static int format;
    static Atom type;
    char buffer[32];

    register int rlen;
    register int result;
    register int (*handler)();
    register char *cmnd;
    register char *ptr;

    command = NULL;
    result = XGetWindowProperty( dpy, MainWin, CommAtom, 0, 1024, True, 
                                 XA_STRING, &type, &format, &len, &left,
                                 &command );
    if( (result!=Success) || (type!=XA_STRING) || (format!=8) )
    {   if( command ) XFree( (char*)command );
        return;
    }

    result = 0;
    ptr = (char*)command;
    if( !*ptr )
    {   /* Tcl/Tk4.0 and later */

        ptr++;
        while( ptr < (char*)command+len )
        {    if( (ptr[0]=='c') && (ptr[1]=='\0') )
             {   ptr += 2;
                 cmnd = (char*)NULL;
                 source = serial = 0;
                 while( (ptr<(char*)command+len) && (*ptr=='-') )
                 {   if( (ptr[1]=='r') && (ptr[2]==' ') )
                     {   sscanf(ptr+3,"%x %d\n",(int*)&source,&serial);
                     } else if( (ptr[1]=='s') && (ptr[2]==' ') )
                         cmnd = ptr+3;
                     while( *ptr ) ptr++;
                     ptr++;
                 }

                 if( !cmnd ) continue;
                 result = ExecuteIPCCommand(cmnd);
                 if( !source || !serial ) continue;

                 buffer[0]='\0';
                 buffer[1]='r';
                 buffer[2]='\0';
                 buffer[3]='-';
                 buffer[4]='r';
                 buffer[5]=' ';
                 buffer[6]= result? '1' : '0';
                 buffer[7]='\0';
                 sprintf(buffer+8,"-s %d",serial);
                 rlen = strlen(buffer+8)+9;

                 /* Return Tcl/Tk v4.0 result! */
                 handler = XSetErrorHandler( HandleIPCError );
                 XChangeProperty(dpy,source,CommAtom, XA_STRING, 8,
                                 PropModeAppend,(unsigned char*)buffer,rlen);
                 XSync(dpy,False);
                 XSetErrorHandler(handler);
             } else /* Unrecognised command! */
             {   while( *ptr ) ptr++;
                 ptr++;
             }
        }

    } else while( *ptr )
    {   /* Tcl/Tk3.0 and later */
        if( *ptr=='C' )
        {   sscanf(ptr+1,"%x %x\n",(int*)&source,&serial);
            while( *ptr && (*ptr!='|') ) ptr++;
            if( *ptr=='|' )
            {   result = ExecuteIPCCommand(ptr+1);
            } else result = 0;

            sprintf(buffer,"R %x 0 %d",serial,result);
            handler = XSetErrorHandler( HandleIPCError );
            XChangeProperty( dpy, source, CommAtom, XA_STRING, 8,
                             PropModeAppend, (unsigned char*)buffer, 
                             strlen(buffer)+1 );
            XSync(dpy,False);
            XSetErrorHandler(handler);
        } 

        /* Next Command! */
        while( *ptr++ );
    }
    XFree( (char*)command );

    if( (result==IPC_Quit) || (result==IPC_Exit) )
        RasMolExit();
}


void SetMouseUpdateStatus( int bool )
{
    if( MouseUpdateStatus != bool )
    {   /* Enable/Disable Pointer Motion Events! */
        attr.event_mask = ExposureMask | ButtonPressMask | ButtonMotionMask 
                        | ButtonReleaseMask;
        if( bool ) attr.event_mask |= PointerMotionMask;
        XChangeWindowAttributes( dpy, CanvWin, CWEventMask, &attr );
    }
    MouseUpdateStatus = bool;
}


void SetMouseCaptureStatus( int bool )
{
    MouseCaptureStatus = bool;
}
                         

static int GetStatus( int mask )
{
    register int status;
    
    status = 0;                             
    if( mask & Button1Mask ) status |= MMLft;
    if( mask & Button2Mask ) status |= MMMid;
    if( mask & Button3Mask ) status |= MMRgt;
    if( mask & ControlMask ) status |= MMCtl;          
    if( mask & ShiftMask )   status |= MMSft;
    return status;
}
  

static int CropRange( int val, int min, int  max )
{
    if( val<min ) return min;
    if( val>max ) return max;
    return val;
}


static void DoneEvents( void )
{
    register Real temp;
    register int index;

    if( HeldButton == YScrlDial )
    {   index = NewScrlY+HeldStep;
#ifdef ORIG
        if( YScrlDial < 3 )
        {   if( index<16 )             
            {   index += YRange-48;
            } else if( index > YRange-32 ) 
                index -= YRange-48;
            NewScrlY = index;
        } else NewScrlY = CropRange(index,16,YRange-32);
#else
        if( index < 16 )
        {   index += YRange-48;
        } else if( index > YRange-32 )
            index -= YRange-48;
        NewScrlY = index;
#endif
    }

    if( NewScrlY != ScrlY )
    {   XClearArea(dpy,YScrlWin,0,ScrlY,16,16,False);
        XCopyArea(dpy,Scrl,YScrlWin,gcon,0,0,16,16,0,NewScrlY);

        temp = ((Real)(NewScrlY-16))/(YRange-48);
        DialValue[YScrlDial] = 2.0*temp - 1.0;
        ReDrawFlag |= (1<<YScrlDial);
        ScrlY = NewScrlY;
    }

    if( HeldButton == XScrlDial )
    {   index = NewScrlX+HeldStep;
#ifdef ORIG
        if( XScrlDial<3 )
        {   if( index < 16 ) 
            {   index += XRange-48;
            } else if( index > XRange-32 ) 
                index -= XRange-48;
            NewScrlX = index;
        } else NewScrlX = CropRange(index,16,XRange-32);
#else
        if( index < 16 )
        {   index += XRange-48;
        } else if( index > XRange-32 )
            index -= XRange-48;
        NewScrlX = index;
#endif
    }

    if( NewScrlX != ScrlX )
    {   XClearArea(dpy,XScrlWin,ScrlX,0,16,16,False);
        XCopyArea(dpy,Scrl,XScrlWin,gcon,0,0,16,16,NewScrlX,0);

        temp = ((Real)(NewScrlX-16))/(XRange-48);
        DialValue[XScrlDial] = 2.0*temp - 1.0;
        ReDrawFlag |= (1<<XScrlDial);
        ScrlX = NewScrlX;
    }
    /* XSync(dpy,False); */
    XFlush(dpy);
}


static int ProcessEvent(  XEvent *event )
{
    register int result;
    register int index;
    register int stat;

    result = 0;
    switch( event->type )
    {   case(ButtonPress):
            {   XButtonPressedEvent *ptr;

                HeldButton = -1;
                ptr = (XButtonPressedEvent*)event;

                if( ptr->window==CanvWin )
                {   stat = GetStatus(ptr->state);
                    ProcessMouseDown(ptr->x,ptr->y,stat);
                } else if( ptr->window==MenuWin )
                {   if( !DisableMenu )
                        if( HandleMenuClick(ptr->x) )
                            result = HandleMenuLoop();
                } else if( ptr->window==XScrlWin )
                {   ReDrawFlag |= RFRotateY;
                    if( ptr->x<16 )
                    {   HeldButton = XScrlDial;
                        HeldStep = -XScrlSkip;
                    } else if( ptr->x>=XRange-16 )
                    {   HeldButton = XScrlDial;
                        HeldStep = XScrlSkip;
                    } else
                    {   index = ptr->x-8;
#ifdef ORIG
                        if( XScrlDial<3 )
                        {   if( index>XRange-32 )
                            {   index -= XRange-48;
                            } else if( index<16 )
                                index += XRange-48;
                            NewScrlX = index;
                        } else NewScrlX = CropRange(index,16,XRange-32);
#else
                        if( index > XRange-32 )
                        {   index -= XRange-48;
                        } else if( index < 16 )
                            index += XRange-48;
                        NewScrlX = index;
#endif
                    }

                } else if( ptr->window==YScrlWin )
                {   ReDrawFlag |= RFRotateX;
                    if( ptr->y<16 )
                    {   HeldButton = YScrlDial;
                        HeldStep = -YScrlSkip;
                    } else if( ptr->y>=YRange-16 )
                    {   HeldButton = YScrlDial;
                        HeldStep = YScrlSkip;
                    } else
                    {   index = ptr->y-8;
#ifdef ORIG
                        if( YScrlDial<3 )
                        {   if( index > YRange-32 )
                            {   index -= YRange-48;
                            } else if( index < 16 )
                                index += YRange-48;
                            NewScrlY = index;
                        } else NewScrlY = CropRange(index,16,YRange-32);
#else
                        if( index > YRange-32 )
                        {   index -= YRange-48;
                        } else if( index < 16 )
                            index += YRange-48;
                        NewScrlY = index;
#endif
                    } 

	        }
            } break;

        case(MotionNotify):
            {   XMotionEvent *ptr;

                ptr = (XMotionEvent*)event;
                if( ptr->window==CanvWin )
                {   stat = GetStatus(ptr->state);
                    ProcessMouseMove(ptr->x,ptr->y,stat);
                } else if( HeldButton == -1 )
                {   if( ptr->window==XScrlWin )
                    {   index = ptr->x-8;
                        NewScrlX = CropRange(index,16,XRange-32);
                    } else /* if( ptr->window==YScrlWin ) */
                    {   index = ptr->y-8;
                        NewScrlY = CropRange(index,16,YRange-32);
                    }
                }
            } break;
             
        case(ButtonRelease):
            {   XButtonReleasedEvent *ptr;

                if( HeldButton != -1 )
                {   /* Three button emulation fix! */
                    DoneEvents();  HeldButton = -1;
                }

                ptr = (XButtonReleasedEvent*)event;
                if( ptr->window==CanvWin )
                {   stat = GetStatus(ptr->state);
                    ProcessMouseUp(ptr->x,ptr->y,stat);
                }
            } break;

        case(KeyPress):
            {   XKeyPressedEvent *ptr;
                static KeySym symbol;
                static char keychar;

                keychar = '\0';
                ptr = (XKeyPressedEvent*)event;
                index = XLookupString(ptr,&keychar,1,&symbol,NULL);
                switch( symbol )
                {   case(XK_Begin):
                    case(XK_Home):  ProcessCharacter(0x01);  break;
                    case(XK_Right): ProcessCharacter(0x06);  break;
                    case(XK_Left):  ProcessCharacter(0x02);  break;
                    case(XK_End):   ProcessCharacter(0x05);  break;
                    case(XK_Up):
                    case(XK_Prior): ProcessCharacter(0x10);  break;
                    case(XK_Down):
                    case(XK_Next):  ProcessCharacter(0x0e);  break;

                    case(XK_F10):   if( !DisableMenu )
                                    {   SelectMenu(0);
                                        result = HandleMenuLoop();
                                    }
                                    break;

                    default:        if( index == 1 )
                                        if( !(ptr->state&Mod1Mask) )
                                        {   if( ProcessCharacter(keychar) )
                                            {   if( ProcessCommand() )
                                                    RasMolExit();

                                                if( !CommandActive )
                                                    ResetCommandLine(0);
                                            }
                                        } else if( !DisableMenu )
                                            if( HandleMenuKey(keychar) )
                                                result = HandleMenuLoop();
                }
            } break;


        case(Expose):
            {   XExposeEvent *ptr;

                ptr = (XExposeEvent*)event;
                if( ptr->window==CanvWin )
                {   if( image ) {
#ifdef MITSHM
                        if( SharedMemFlag )
                        {   XShmPutImage( dpy, CanvWin, gcon, image,
                                          ptr->x, ptr->y, ptr->x, ptr->y,
                                          ptr->width, ptr->height, False);
                            XSync(dpy,False);
                        } else
#endif 
                        XPutImage( dpy, CanvWin, gcon, image,
                                   ptr->x, ptr->y, ptr->x, ptr->y,
                                   ptr->width, ptr->height );
                    } else XClearWindow( dpy, CanvWin );
		    
                } else if( ptr->window==MainWin )
                {   DrawMainWin();
                    DrawMenuBar();
                } else if( ptr->window==XScrlWin )
                {   DrawXScroll();
                } else if( ptr->window==YScrlWin )
                    DrawYScroll();
                XFlush(dpy);
            } break;

        case(EnterNotify):
            {   XCrossingEvent *ptr;

                ptr = (XCrossingEvent*)event;
                if( ptr->detail != NotifyInferior )
                {   if( LocalMap )
                        XInstallColormap(dpy,lmap);
                }
#ifdef DIALBOX
                if( UseDials )
                    GetDialState();
#endif
            }
            break;

        case(LeaveNotify):
            if( LocalMap )
            {   XCrossingEvent *ptr;

                ptr = (XCrossingEvent*)event;
                if( ptr->detail != NotifyInferior )
                    XUninstallColormap(dpy,lmap);
            }
            break;

        case(ConfigureNotify):
            {   XConfigureEvent *ptr;
                register int wide,high;

                ptr = (XConfigureEvent*)event;
                high = CropRange(ptr->height,MinHeight,MaxHeight);
                wide = CropRange(ptr->width, MinWidth, MaxWidth );

                if( (wide!=MainWide) || (high!=MainHigh) )
                    ReSizeWindow(wide,high);
            } break;

        case(ClientMessage):
            {   XClientMessageEvent *ptr;

                ptr = (XClientMessageEvent*)event;
                if( (ptr->message_type==ProtoXAtom) && 
                    (ptr->data.l[0]==DelWinXAtom) )
                    RasMolExit();
            } break;

        case(PropertyNotify):
            {   XPropertyEvent *ptr;

                ptr = (XPropertyEvent*)event;
                if( (ptr->atom==CommAtom) &&
                    (ptr->state==PropertyNewValue) )
                    HandleIPCCommand();
            } break;

        case(MapNotify):
            DrawXScroll();
            DrawYScroll();
            DrawMainWin();
            DrawMenuBar();
            break;

        default:  
#ifdef DIALBOX
            if( event->type == DialEvent )
                HandleDialEvent(  (XDeviceMotionEvent*)event );
#endif
            break;
    }
    return result;
}



/*=========================*/
/*  Modal Dialog Handling  */
/*=========================*/

static int HandleMenuLoop( void )
{
    register unsigned int mask;
    register int result;
    register int done;
    auto XEvent event;

    /* Passive Pointer Grab */
    mask = ButtonPressMask | ButtonReleaseMask | ButtonMotionMask;
    XGrabPointer(dpy,MenuWin,False,mask,
                 GrabModeAsync,GrabModeAsync,
                 None,None,CurrentTime);

    HeldButton = -1;
    MenuFocus = True;
    DrawMenuBar();

    result = 0;
    done = False;
    while( !done )
    {   XNextEvent( dpy, &event );
        switch( event.type )
        {   case(Expose):
                {   XExposeEvent *ptr;

                    ptr = (XExposeEvent*)&event;
                    if( ptr->window==PopUpWin )
                    {   DrawPopUpMenu();
                    } else ProcessEvent(&event);
                } break;

            case(ButtonPress): 
                {   XButtonPressedEvent *ptr;

                    ptr = (XButtonPressedEvent*)&event;
                    /* All Events Relative to MenuWin */
                    if( (ptr->y>=0) && (ptr->y<=FontHigh+5) )
                    {   HandleMenuClick(ptr->x);
                    } else if( PopUpFlag )
                    {   HandleItemClick(ptr->x,ptr->y);
                    } else done = True;
                } break;

            case(MotionNotify):
                    if( ItemFocus )
                    {   XMotionEvent *ptr;

                        ptr = (XMotionEvent*)&event;
                        /* All Events Relative to MenuWin */
                        if( (ptr->y>=0) && (ptr->y<=FontHigh+5) )
                        {   HandleMenuClick( ptr->x );
                        } else if( PopUpFlag )
                            HandleItemMove(ptr->x,ptr->y);
                    } break;

            case(ButtonRelease):
                    {   XButtonReleasedEvent *ptr;

                        ptr = (XButtonReleasedEvent*)&event;
                        /* All Events Relative to MenuWin */
                        if( (ptr->y>=0) && (ptr->y<=FontHigh+5) )
                        {   if( HandleMenuClick( ptr->x ) )
                            {   SelectFirstItem(MenuBarSelect);
                                DrawPopUpMenu();
                            } else done = True;
                        } else if( PopUpFlag )
                        {   if( ItemFocus )
                                HandleItemClick(ptr->x,ptr->y);
                            if( ItemFlag )
                                result = (MenuBarSelect<<8) +
                                         MenuItemSelect+1;
                            done = True;
                        } else done = False;
                        ItemFocus = False;
                    }
                    break;
     
            case(KeyPress):
                if( !ItemFocus )
                {   XKeyPressedEvent *ptr;
                    static KeySym symbol;
                    static char keychar;
                    register int index;

                    ptr = (XKeyPressedEvent*)&event;
                    index = XLookupString(ptr,&keychar,1,&symbol,NULL);
                    switch( symbol )
                    {   case(XK_Right): index = MenuBarSelect+1;
                                        if( index != MenuBarCount )
                                        {   SelectMenu( index );
                                        } else SelectMenu( 0 );
                                        break;

                        case(XK_Left):  if( MenuBarSelect )
                                        {   SelectMenu( MenuBarSelect-1 );
                                        } else SelectMenu( MenuBarCount-1 );
                                        break;

                        case(XK_Up):    if( !PopUpFlag )
                                        {   PopUpFlag = True;
                                            SelectMenu(MenuBarSelect);
                                        } else SelectPrevItem();
                                        break;

                        case(XK_Down):  if( !PopUpFlag )
                                        {   PopUpFlag = True;
                                            SelectMenu(MenuBarSelect);
                                        } else SelectNextItem();
                                        break;

                        case(XK_KP_Enter):
                        case(XK_Linefeed):
                        case(XK_Return):   if( PopUpFlag && ItemFlag )
                                               result = (MenuBarSelect<<8) +
                                                        MenuItemSelect+1;
                                           done = True;
                                           break;

                        default:    if( (index==1) && (keychar>=' ') ) 
                                    {   if( !(ptr->state&Mod1Mask) )
                                        {   if( PopUpFlag )
                                            {   result = HandleItemKey(keychar);
                                                if( result ) done = True;
                                            } else HandleMenuKey(keychar);
                                        } else HandleMenuKey(keychar);
                                    }
                    }
                } break;


            case(ConfigureNotify):  /* done = True; */
            default:                ProcessEvent(&event);
        }
    }

    /* Passive Grab Release */
    XUngrabPointer(dpy,CurrentTime);

    XUnmapWindow(dpy,PopUpWin);
    PopUpFlag = False;
    MenuFocus = False;
    DrawMenuBar();
    return result;
}


int FetchEvent( int wait )
{
    register int result;
    auto XEvent event;

    NewScrlX = ScrlX;
    NewScrlY = ScrlY;

    if( HeldButton != -1 ) wait = False;
    while( XPending(dpy) || (wait && !ReDrawFlag) )
    {   XNextEvent( dpy, &event );
        result = ProcessEvent(&event);
        if( result ) return result;
    }
    DoneEvents();
    return 0;
}


int LookUpColour( char *name, int *red, int *grn, int *blu )
{
    static XColor exact, close;
    register Colormap map;

    map = (LocalMap)? lmap : cmap;
    if( XLookupColor(dpy,map,name,&exact,&close) )
    {   *red = exact.red>>8;
        *grn = exact.green>>8;
        *blu = exact.blue>>8;
        return True;
    }
    return False;
}


void BeginWait( void )
{
    if( UseHourGlass )
    {   XDefineCursor(dpy,CanvWin,hglass);
        XDefineCursor(dpy,MainWin,hglass);
        XFlush(dpy);
    }
}


void EndWait( void )
{
    if( UseHourGlass )
    {   XDefineCursor(dpy,CanvWin,cross);
        XDefineCursor(dpy,MainWin,arrow);
        XFlush(dpy);
    }
}


void CloseDisplay( void )
{
#ifdef DIALBOX
    register int num;
#endif

    /* FatalXError! */
    if( !dpy ) return;

    if( image ) 
    {
#ifdef MITSHM
        if( SharedMemFlag )
        {   XShmDetach( dpy, &xshminfo );
            image->data = (char*)NULL;
            shmdt( xshminfo.shmaddr );
        }
#endif
        XDestroyImage( image );
    }

    if( *TkInterp )
    {   XGrabServer( dpy );
        DeRegisterInterpName(TkInterp);
        XUngrabServer( dpy );
    }

#ifdef DIALBOX
    if( UseDials )
    {   if( UseDialLEDs )
            for( num=0; num<8; num++ )
                SetDialLabel(num,"");
        XCloseDevice(dpy,Dials);
    }
#endif
    XCloseDisplay( dpy );
}

