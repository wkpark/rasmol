/* x11win.c
 * RasMol2 Molecular Graphics
 * Roger Sayle, February 1994
 * Version 2.3
 */
#ifndef sun386
#include <stdlib.h>
#endif
#include <string.h>
#include <stdio.h>
#include <math.h>

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


#ifdef DIALBOX
#include <X11/extensions/XInput.h>

static int UseDialLEDs;
static XEventClass DialClass;
static XDevice *Dials;
static XID DialIdent;
static int DialEvent;
static int UseDials;
#endif


#ifdef MITSHM
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <X11/extensions/XShm.h>

#ifdef __STDC__
Bool XShmQueryExtension( Display* );
#else /* non-ANSI C compiler */
Bool XShmQueryExtension();
#endif

XShmSegmentInfo xshminfo;
int SharedMemOption;
int SharedMemFlag;
#endif


#define ButWide    96
#define ButHigh    32
#define XScrlDial  1 /*1*/
#define YScrlDial  0 /*0*/
#define XScrlSkip  8
#define YScrlSkip  8


/* Determine Mouse Sensitivity! */
#define IsClose(u,v) (((u)>=(v)-1) && ((u)<=(v)+1))
#define MinHeight    (ButMax*(ButHigh+18)+18)
#define MinWidth     (ButWide+127+32)


/* Increased point size for NCD X-terminals! */
#define NrmFontMax 4
static char *NrmFont[] = {
        "-*-lucida-bold-r-*-*-12-*",
        "-*-helvetica-bold-r-*-*-12-*",
        "-*-lucida-bold-r-*-*-14-*",
        "-*-helvetica-bold-r-*-*-14-*" };

static Cursor cross;
static Cursor arrow;
static Cursor hglass;
static Pixmap Scrl;
static Pixmap tilepix;
static Pixmap uppix, dnpix;
static Pixmap lfpix, rgpix;
static Pixmap ButUp, ButDn;
static XFontStruct *NrmInfo;
static XSetWindowAttributes attr;
static Window XScrlWin, YScrlWin;
static Window MainWin;
static Window CanvWin;
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


static int MouseMode;
static int InitX, InitY;
static int HeldButton;
static int HeldStep;

static Window OptWin[ButMax];
static char *OptPtr[ButMax];
static int OptLen[ButMax];
static int OptCount;

static int MaxWidth;
static int MaxHeight;
static int MainWide, MainHigh;
static int ScrlX, NewScrlX;
static int ScrlY, NewScrlY;
static int PixDepth;
static int LocalMap;


/* WM_PROTOCOLS */
static char TkInterp[10];
static Atom DelWinXAtom;
static Atom ProtoXAtom;
static Atom InterpAtom;
static Atom CommAtom;



static void FatalGraphicsError(ptr)
    char *ptr;
{
    char buffer[80];

    sprintf(buffer,"Graphics Error: %s!",ptr);
    RasMolFatalExit(buffer);
}


void AllocateColourMap()
{
#ifdef EIGHTBIT
    static XColor Col;
    register int i,j;

    if( LocalMap )
    {   XSetWindowColormap(dpy,MainWin,cmap);
        XSetWindowColormap(dpy,CanvWin,cmap);
        XUninstallColormap(dpy,lmap);
        XFreeColormap(dpy,lmap);
        LocalMap = False;
    } else if( IdentCount )
        XFreeColors(dpy,cmap,Ident,IdentCount,(long)0);
    IdentCount = 0;


    for( i=0; i<256; i++ )
        if( ULut[i] )
        {   Col.red   = RLut[i]<<8;
            Col.green = GLut[i]<<8;
            Col.blue  = BLut[i]<<8;
            Col.flags = DoRed | DoGreen | DoBlue;
            if( !XAllocColor(dpy,cmap,&Col) )
                break;
            Ident[IdentCount++] = Col.pixel;
            Lut[i] = Col.pixel;
        }

    if( i<256 )
    {   lmap = XCopyColormapAndFree(dpy,cmap);
        LocalMap = True;

        for( j=0; j<5; j++ )
        {   Col.red   = RLut[j]<<8;
            Col.green = GLut[j]<<8;
            Col.blue  = BLut[j]<<8;
            XAllocColor(dpy,cmap,&Col);
            Lut[i] = Col.pixel;
        }

        for( j=i; j<256; j++ )
            if( ULut[j] )
            {   Col.red   = RLut[j]<<8;
                Col.green = GLut[j]<<8;
                Col.blue  = BLut[j]<<8;
                XAllocColor(dpy,lmap,&Col);
                Lut[j] = Col.pixel;
            }
        XSetWindowColormap(dpy,MainWin,lmap);
        XSetWindowColormap(dpy,CanvWin,lmap);
        XInstallColormap(dpy,lmap);
    }
#endif

    XSetWindowBackground(dpy,CanvWin,(long)Lut[5]);
}



static void OpenCanvas( x, y )
    int x, y;
{
    register unsigned long mask;

    mask = CWEventMask;
    attr.event_mask = ExposureMask | ButtonPressMask | ButtonMotionMask 
                    | ButtonReleaseMask;
    attr.cursor = cross;                           mask |= CWCursor;
    attr.background_pixel = Lut[0];                mask |= CWBackPixel;

    CanvWin = XCreateWindow(dpy, MainWin, 18, 18, x, y, 0, CopyFromParent,
                            InputOutput, vis, mask, &attr );
}


static void OpenFonts()
{
    register int i;
    cross = XCreateFontCursor(dpy,XC_tcross);
    arrow = XCreateFontCursor(dpy,XC_top_left_arrow);

    for( i=0; i<NrmFontMax; i++ )
        if( (NrmInfo=XLoadQueryFont(dpy,NrmFont[i])) ) 
            return;
    FatalGraphicsError("Unable to find suitable font");
}


static void OpenCursors()
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


static void OpenColourMap()
{
#ifdef EIGHTBIT
    static XColor Col;
    register int i;

    Col.flags = DoRed | DoGreen | DoBlue;

    for( i=0; i<5; i++ )
    {   Col.red   = RLut[i]<<8;
        Col.green = GLut[i]<<8;
        Col.blue  = BLut[i]<<8;
        if( !XAllocColor(dpy,cmap,&Col) )
        {   cmap = XCopyColormapAndFree(dpy,cmap);
            XAllocColor(dpy,cmap,&Col);
        } 
        Lut[i] = Col.pixel;
    }

    LocalMap = False;
    IdentCount = 0;
    Lut[5]=Lut[0];
#endif
}


static int RegisterInterpName( name )
    char *name;
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
    {   if( (type!=None) && registry ) XFree(registry);

        sprintf(buffer,"%x %s",MainWin,name);
        XChangeProperty( dpy, RootWindow(dpy,0), InterpAtom, XA_STRING, 
                         8, PropModeReplace, buffer, strlen(buffer)+1 );
        return( True );
    }

    ptr = (char*)registry;
    while( *ptr )
    {   /* Skip Window ID */
        while( *ptr++ != ' ' )
            if( !*ptr ) break;

        /* Compare Interp Name */
        if( !strcmp(ptr,name) )
        {   XFree(registry);
            return(False);
        }

        while( *ptr++ );
    }

    XFree(registry);
    sprintf(buffer,"%x %s",MainWin,name);
    XChangeProperty( dpy, RootWindow(dpy,0), InterpAtom, XA_STRING, 
                     8, PropModeAppend, buffer, strlen(buffer)+1 );
    return( True );
}


static void DeRegisterInterpName( name )
    char *name;
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
        if( registry ) XFree(registry);
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
    XFree( registry );
}


static void OpenIPCComms()
{
    register int i;

    CommAtom = XInternAtom( dpy, "Comm", False );
    DelWinXAtom = XInternAtom(dpy, "WM_DELETE_WINDOW", False);
    /* XSetWMProtocols(dpy,MainWin,&DelWinXAtom,True); */
    if( (ProtoXAtom = XInternAtom(dpy,"WM_PROTOCOLS",False)) )
        XChangeProperty( dpy, MainWin, ProtoXAtom, XA_ATOM, 32, 
                        PropModeReplace, (Byte*)&DelWinXAtom, True );

    InterpAtom = XInternAtom( dpy, "InterpRegistry", False );

    XGrabServer( dpy );
    if( !RegisterInterpName("rasmol") )
    {   strcpy(TkInterp,"rasmol #0");
        for( i=1; i<10; i++ )
        {    TkInterp[8] = i+'0';
             if( RegisterInterpName(TkInterp) )
                 break;
        }

        if( i==10 ) *TkInterp = 0;
    } else strcpy(TkInterp,"rasmol");
    XUngrabServer( dpy );
}


static void DrawBox(wdw,pos,x1,y1,x2,y2)
    Drawable wdw; int pos,x1,y1,x2,y2;
{
    register unsigned long colour;
    register int ux,uy,lx,ly;
    register int index;

    lx=x1; ly=y1; ux=x2; uy=y2;
    colour = Lut[pos? 3 : 1 ];
    XSetForeground(dpy,gcon,colour);
    for( index=0; index<3; index++ )
    {   XDrawLine(dpy,wdw,gcon,lx,ly,ux,ly);
        XDrawLine(dpy,wdw,gcon,lx,ly,lx,uy);
        lx++; ly++; ux--; uy--;
    }

    
    lx=x1; ly=y1; ux=x2; uy=y2;
    colour = Lut[pos? 1 : 3 ];
    XSetForeground(dpy,gcon,colour);
    for( index=0; index<3; index++ )
    {   XDrawLine(dpy,wdw,gcon,ux,ly,ux,uy);
        XDrawLine(dpy,wdw,gcon,lx,uy,ux,uy);
        lx++; ly++; ux--; uy--;
    }
}


static void OpenButtons()
{
    register unsigned long mask;
    register int index;

    ButUp = XCreatePixmap( dpy, MainWin, ButWide+8, ButHigh+8, PixDepth );
    ButDn = XCreatePixmap( dpy, MainWin, ButWide+8, ButHigh+8, PixDepth );

    XSetForeground( dpy, gcon, (unsigned long)Lut[2] );
    XFillRectangle( dpy, ButUp, gcon, 0, 0, ButWide+7, ButHigh+7 );
    XFillRectangle( dpy, ButDn, gcon, 0, 0, ButWide+7, ButHigh+7 );

    XSetForeground( dpy, gcon, (unsigned long)Lut[0] );
    XDrawRectangle( dpy, ButUp, gcon, 0, 0, ButWide+7, ButHigh+7 );
    XDrawRectangle( dpy, ButDn, gcon, 0, 0, ButWide+7, ButHigh+7 );
    
    DrawBox( ButUp, True,  1, 1, ButWide+6, ButHigh+6 );
    DrawBox( ButDn, False, 1, 1, ButWide+6, ButHigh+6 );


    mask = CWEventMask;
    attr.event_mask = ButtonPressMask | ButtonReleaseMask;

    for( index=0; index<ButMax; index++ )
        OptWin[index] = XCreateWindow( dpy, MainWin,
                                       XRange+58, index*(ButHigh+18)+15,
                                       ButWide+6, ButHigh+6, 0,
                                       CopyFromParent, InputOnly, vis,
                                       mask, &attr );
    OptCount = 0;
}


static void ReDrawButton( num, pos )
    int num, pos;
{
    register int deltaX, deltaY;
    register int xpos, ypos;

    xpos = XRange+61;
    ypos = num*(ButHigh+18)+18;

    deltaY = (ButHigh-(NrmInfo->ascent+NrmInfo->descent))/2+NrmInfo->ascent;
    deltaX = (ButWide-XTextWidth(NrmInfo,OptPtr[num],OptLen[num]))/2;

    XSetFont( dpy, gcon, NrmInfo->fid );
    XSetForeground( dpy, gcon, (unsigned long)Lut[0] );
    XCopyArea( dpy, pos? ButUp : ButDn, MainWin, gcon,
               0, 0, ButWide+8, ButHigh+8, xpos-4, ypos-4 );
    XDrawString( dpy, MainWin, gcon, xpos+deltaX, ypos+deltaY, 
                 OptPtr[num], OptLen[num] );
    XFlush(dpy);
}


void NewMenu( num, option )
    int num; char **option;
{
    register char *ptr;
    register int start,stop;
    register int index;
    register int len;

    for( index=0; index<num; index++ )
    {   len = 0;
        OptPtr[index] = ptr = *option++;
        while( *ptr++ ) len++;
        OptLen[index] = len;

        ReDrawButton( index, True );
    }
    MenuDisable = False;
    OptCount = num;

    if( num<ButMax )
    {    start = num*(ButHigh+18)+13;
         stop = ButMax*(ButHigh+18)+5;
         XSetForeground(dpy,gcon,(unsigned long)Lut[2]);
         XFillRectangle( dpy, MainWin, gcon, 
                         XRange+56, start, ButWide+10, stop-start );
    }
    XSync(dpy,False);
}


static void OpenScrollBars()
{
    register unsigned long mask;

    Scrl = XCreatePixmap( dpy, MainWin, 16, 16, PixDepth );
    XSetForeground(dpy,gcon,(unsigned long)Lut[2]); 
    XFillRectangle(dpy,Scrl,gcon,0,0,15,15);
    XSetForeground(dpy,gcon,(unsigned long)Lut[0]); 
    XDrawRectangle(dpy,Scrl,gcon,0,0,15,15);
    DrawBox( Scrl, True, 1, 1, 14, 14 );

    tilepix = XCreatePixmapFromBitmapData( dpy, MainWin, (char*)ScrlTile,
                                           8, 8, Lut[0], Lut[2], PixDepth );

    mask = CWEventMask;
    attr.event_mask = ExposureMask | ButtonPressMask | ButtonMotionMask 
                    | ButtonReleaseMask;
    attr.background_pixmap = tilepix;              mask |= CWBackPixmap;

    XScrlWin = XCreateWindow( dpy, MainWin, 18, YRange+27, XRange, 16, 0,
                              CopyFromParent, InputOutput, vis, mask, &attr );
    lfpix = XCreatePixmapFromBitmapData( dpy, MainWin, (char*)LfArrow,
                                         16, 16, Lut[0], Lut[2], PixDepth );
    rgpix = XCreatePixmapFromBitmapData( dpy, MainWin, (char*)RgArrow,
                                         16, 16, Lut[0], Lut[2], PixDepth );


    YScrlWin = XCreateWindow( dpy, MainWin, XRange+27, 18, 16, YRange, 0,
                              CopyFromParent, InputOutput, vis, mask, &attr );
    uppix = XCreatePixmapFromBitmapData( dpy, MainWin, (char*)UpArrow,
                                         16, 16, Lut[0], Lut[2], PixDepth );
    dnpix = XCreatePixmapFromBitmapData( dpy, MainWin, (char*)DnArrow,
                                         16, 16, Lut[0], Lut[2], PixDepth );

    ScrlX = (XRange/2)-8;
    ScrlY = (YRange/2)-8;
}

static void ReDrawXScroll()
{
    XCopyArea(dpy,rgpix,XScrlWin,gcon,0,0,16,16,XRange-16,0);
    XCopyArea(dpy,Scrl ,XScrlWin,gcon,0,0,16,16,ScrlX,0);
    XCopyArea(dpy,lfpix,XScrlWin,gcon,0,0,16,16,0,0);
}

static void ReDrawYScroll()
{
    XCopyArea(dpy,dnpix,YScrlWin,gcon,0,0,16,16,0,YRange-16);
    XCopyArea(dpy,Scrl ,YScrlWin,gcon,0,0,16,16,0,ScrlY);
    XCopyArea(dpy,uppix,YScrlWin,gcon,0,0,16,16,0,0);
}


void UpdateScrollBars()
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
static void SetDialLabel( num, ptr )
    int num; char *ptr;
{
    static XStringFeedbackControl ctrl;
    static KeySym text[8];
    register int length;

    length = 0;
    while( *ptr )
       text[length++] = *ptr++;

    ctrl.id = num;
    ctrl.class = DialIdent;
    ctrl.num_keysyms = length;
    ctrl.syms_to_display = text;
    XChangeFeedbackControl(dpy,Dials,DvString,&ctrl);
}


static void OpenDialsBox()
{
    register XFeedbackState *list;
    register XFeedbackState *feed;
    register XDeviceInfo *devlist;
    register XDeviceInfo *ptr;
    register Atom devtype;
    register int index;
    static int count;

    UseDials = False;
    devlist = XListInputDevices(dpy,&count);
    devtype = XInternAtom(dpy,XI_KNOB_BOX,True );
    if( (devtype==None) || !devlist ) return;

    ptr = devlist;
    for( index=0; index<count; index++ )
        if( (ptr->use==IsXExtensionDevice) && (ptr->type==devtype) )
        {   Dials = XOpenDevice(dpy,ptr->id);
            DialIdent = ptr->inputclassinfo->class;
            UseDials = True;
            break;
        } else ptr++;
    /* XFreeDeviceList(devlist); */
    if( !UseDials ) return;

    UseDialLEDs = 0;
    feed = list = XGetFeedbackControl( dpy, Dials, &count );
    for( index=0; index<count; index++ )
    {   if( feed->class == StringFeedbackClass )
            UseDialLEDs++;
        feed = (XFeedbackState*)(((char*)feed) + feed->length);
    }
    XFreeFeedbackList( list );

    if( UseDialLEDs >= 8 )
    {   SetDialLabel(0,"ROTATE X");
        SetDialLabel(1,"ROTATE Y");
        SetDialLabel(2,"ROTATE Z");
        SetDialLabel(3,"  ZOOM  ");

        SetDialLabel(4,"TRANS X ");
        SetDialLabel(5,"TRANS Y ");
        SetDialLabel(6,"TRANS Z ");
        SetDialLabel(7,"  SLAB  ");
    } else UseDialLEDs = False;

    DeviceMotionNotify( Dials, DialEvent, DialClass );
    XSelectExtensionEvent( dpy, MainWin, &DialClass, 1 );
    XSelectExtensionEvent( dpy, CanvWin, &DialClass, 1 );
    XSelectExtensionEvent( dpy, XScrlWin, &DialClass, 1 );
    XSelectExtensionEvent( dpy, YScrlWin, &DialClass, 1 );

    for( index=0; index<ButMax; index++ )
       XSelectExtensionEvent( dpy, OptWin[index], &DialClass, 1 );
}


static void HandleDialEvent( ptr )
    XDeviceMotionEvent *ptr;
{
    register double temp;
    register int count;
    register int value;
    register int index;
    register int num;

    /* Limit Number of Dials */
    count = 8 - ptr->first_axis;
    if( count > ptr->axes_count )
        count = ptr->axes_count;

    for( index=0; index<count; index++ )
        if( value = ptr->axis_data[index] )
        {   num = ptr->first_axis+index;
            temp = DialValue[num]+(value/1024.0);
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
#endif


static void ReDrawMain()
{
    register int index;

    DrawBox(MainWin,True,0,0,MainWide,MainHigh);
    DrawBox(MainWin,False,15,15,XRange+20,YRange+20);
    DrawBox(MainWin,False,XRange+24,15,XRange+45,YRange+20);
    DrawBox(MainWin,False,15,YRange+24,XRange+20,YRange+45);

    for( index=0; index<OptCount; index++ )
        ReDrawButton( index, True );
}


static void ReSizeWindow( wide, high )
    int wide, high;
{
    register Real xpos;
    register Real ypos;
    register int index;

    xpos = (XRange>48)? (Real)(ScrlX-16)/(XRange-48) : 0.0;
    ypos = (YRange>48)? (Real)(ScrlY-16)/(YRange-48) : 0.0;

    wide = (wide & ~3) | ((ButWide+79) & 3);
    MainWide = wide;  XRange = wide-(ButWide+79);  WRange = XRange>>1;
    MainHigh = high;  YRange = high-61;            HRange = YRange>>1;
    Range = MinFun(XRange,YRange);

    XResizeWindow( dpy, CanvWin, XRange, YRange);
    XMoveResizeWindow( dpy, XScrlWin, 18, YRange+27, XRange, 16 );
    XMoveResizeWindow( dpy, YScrlWin, XRange+27, 18, 16, YRange );

    for( index=0; index<ButMax; index++ )
        XMoveWindow( dpy, OptWin[index], XRange+58, index*(ButHigh+18)+15 );

    NewScrlX = ScrlX = (xpos*(XRange-48))+16;
    NewScrlY = ScrlY = (ypos*(YRange-48))+16;

    XClearWindow( dpy, MainWin );
    XClearWindow( dpy, CanvWin );

    ReDrawXScroll();
    ReDrawYScroll();
    ReDrawMain();

    ReDrawFlag |= RFReSize;
    XSync(dpy,True);
}


int FatalXError( ptr )
    Display *ptr;
{
    dpy = (Display*)NULL;
    RasMolFatalExit("*** Fatal X11 I/O Error! ***");
    /* Avoid Compilation Warnings! */
    return( (int)ptr );
}


int OpenDisplay( x, y )
    int x, y;
{
    register int i,num;
    register char *ptr;
    register unsigned long mask;
    static XVisualInfo visinfo;
    static XSizeHints size;
    static Window rootwin;
    static Pixmap icon;


    image = (XImage*)NULL;

    MouseMode = MMRasMol;
    UseHourGlass = True;
    MenuDisable = False;
    HeldButton = -1;

    for( i=0; i<8; i++ )
         DialValue[i] = 0.0;

#ifdef EIGHTBIT
    RLut[0]=0;   GLut[0]=0;   BLut[0]=0;    ULut[0]=True;
    RLut[1]=100; GLut[1]=100; BLut[1]=100;  ULut[1]=True;
    RLut[2]=150; GLut[2]=150; BLut[2]=150;  ULut[2]=True;
    RLut[3]=200; GLut[3]=200; BLut[3]=200;  ULut[3]=True;
    RLut[4]=255; GLut[4]=255; BLut[4]=255;  ULut[4]=True;
#else
    Lut[0] = 65793*0;
    Lut[1] = 65793*64;
    Lut[2] = 65793*128;
    Lut[3] = 65793*196;
    Lut[4] = 65793*255;
#endif

    XRange = x;  WRange = XRange>>1;
    YRange = y;  HRange = YRange>>1;
    Range = MinFun(XRange,YRange);

    if( (dpy=XOpenDisplay(NULL)) == NULL )
        return( 0 );

    num = DefaultScreen(dpy);
    rootwin = RootWindow(dpy,num);
    XSetIOErrorHandler( FatalXError );

#ifdef EIGHTBIT
    if( !XMatchVisualInfo(dpy,num,8,PseudoColor,&visinfo) )
    {   XCloseDisplay(dpy);
        return(0);
    } else PixDepth = 8;
#else
    if( XMatchVisualInfo(dpy,num,32,TrueColor,&visinfo) )
    {   PixDepth = 32;
    } else if( XMatchVisualInfo(dpy,num,24,TrueColor,&visinfo) )
    {   PixDepth = 24;
    } else /* No suitable display! */
    {   XCloseDisplay(dpy);
        return(0);
    }
#endif

    vis = visinfo.visual;
    if( vis != DefaultVisual(dpy,num) )
    {   cmap = XCreateColormap(dpy,rootwin,vis,AllocNone);
    } else cmap = DefaultColormap(dpy,num);

    OpenFonts();
    OpenColourMap();

    MaxWidth = DisplayWidth(dpy,num);
    MaxHeight = DisplayHeight(dpy,num);
    MainWide = x+ButWide+79;
    MainHigh = y+61;

    mask = CWEventMask;
    attr.event_mask = ExposureMask | KeyPressMask | StructureNotifyMask
                    | EnterWindowMask | LeaveWindowMask | PropertyChangeMask;
    attr.background_pixel = Lut[2];     mask |= CWBackPixel;
    attr.border_pixel = Lut[2];         mask |= CWBorderPixel;
    attr.colormap = cmap;               mask |= CWColormap;
    attr.cursor = arrow;                mask |= CWCursor;

    MainWin = XCreateWindow(dpy, rootwin, 0, 0, MainWide, MainHigh, 2,
			    PixDepth, InputOutput, vis, mask, &attr );

    gcon = XCreateGC(dpy,MainWin,0L,NULL);
    /* DefaultGC(dpy,num) */

    XSetGraphicsExposures(dpy,gcon,False);
    icon = XCreateBitmapFromData(dpy,MainWin,(char*)icon_bits,
                                 icon_width,icon_height );

    size.flags = PMinSize | PMaxSize;
    size.min_width = MinWidth;    size.max_width = MaxWidth;
    size.min_height = MinHeight;  size.max_height = MaxHeight;
    XSetStandardProperties(dpy, MainWin, "RasMol Version 2.3",
                           "RasMol", icon, NULL, 0, &size );

    hints.icon_pixmap = icon;
    hints.flags = IconPixmapHint;
    XSetWMHints(dpy,MainWin,&hints);

    OpenCanvas( x, y );
    OpenScrollBars();
    OpenCursors();
    OpenButtons();
    OpenIPCComms();

#ifdef DIALBOX
    OpenDialsBox();
#endif

#ifdef MITSHM
    ptr = DisplayString(dpy);
    if( !ptr || (*ptr==':') || !strncmp(ptr,"localhost:",10) || 
        !strncmp(ptr,"unix:",5) || !strncmp(ptr,"local:",6) )
    {   SharedMemOption = XShmQueryExtension(dpy);
    } else SharedMemOption = False;
    SharedMemFlag = False;
#endif

    XMapSubwindows(dpy,MainWin);
    XMapWindow(dpy,MainWin);

    ReDrawMain();
    ReDrawXScroll();
    ReDrawYScroll();
    XSync(dpy,False);

    return( ConnectionNumber(dpy) );
}


#ifdef MITSHM
Pixel *AllocSharedMem( size )
    int size;
{
    register Pixel *ptr;

    SharedMemFlag = False;
    if( SharedMemOption )

    if( !SharedMemFlag ) 
        ptr = (Pixel*)malloc(size);
    return( ptr );
}


void FreeSharedMem( ptr )
    Pixel *ptr;
{
    if( SharedMemFlag )
    {   XShmDetach(dpy,&xshminfo);
        shmdt(xshminfo.shmaddr);
    } else free( ptr );
}
#endif



Pixel *CreateImage()
{
    register Long size, temp;
    register Pixel *ptr;

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

    size = (Long)XRange*YRange*sizeof(Pixel) + 32;

#ifdef MITSHM
    if( SharedMemOption )
    {   SharedMemFlag = False;
        image = XShmCreateImage( dpy, vis, PixDepth, ZPixmap, 
                                 NULL, &xshminfo, XRange, YRange );
        if( image )
        {   temp = (Long)image->bytes_per_line * image->height;
            if( temp > size ) size = temp;
            xshminfo.shmid = shmget( IPC_PRIVATE, size, IPC_CREAT|0777 );
            if( xshminfo.shmid != -1 ) 
            {   xshminfo.shmaddr = (char*)shmat(xshminfo.shmid,0,0);
                if( xshminfo.shmaddr != (char*)-1 )
                {   image->data = xshminfo.shmaddr;
                    ptr = (Pixel*)xshminfo.shmaddr;
                    xshminfo.readOnly = True;

                    SharedMemFlag = XShmAttach( dpy, &xshminfo );
                    XSync(dpy,False);
                }
                /* Always Destroy Shared Memory Ident */
                shmctl( xshminfo.shmid, IPC_RMID, 0 );
            }

            if( !SharedMemFlag )
            {   XDestroyImage( image );
            } else return( ptr );
        }
    }

#endif
    if( (ptr = (Pixel*)malloc( size )) )
    {   image = XCreateImage( dpy, vis, PixDepth, ZPixmap, 0, (char*)ptr, 
                              XRange, YRange, ((PixDepth>8)?32: 8) , 0 );
        if( !image ) return( (Pixel*)NULL );
    } else image = (XImage*)NULL;
    return( ptr );
}


void TransferImage()
{   
#ifdef MITSHM
    if( SharedMemFlag )
    {   XShmPutImage(dpy,CanvWin,gcon,image,0,0,0,0,XRange,YRange,False);
    } else
#endif
    XPutImage( dpy, CanvWin, gcon, image, 0, 0, 0, 0, XRange, YRange );
    XFlush(dpy);
}


void ClearImage()
{
    XClearWindow( dpy, CanvWin );
    XFlush(dpy);
}

static int HandleIPCError( disp, ptr )
    Display *disp;  XErrorEvent *ptr;
{
    return( 0 );
}


static void HandleIPCCommand()
{
    static unsigned long len,left;
    static unsigned char *command;
    static Window source;
    static int serial;
    static int format;
    static Atom type;
    char buffer[32];
    int (*handler)();

    register int result;
    register char *ptr;

    command = NULL;
    result = XGetWindowProperty( dpy, MainWin, CommAtom, 0, 1024, True, 
                                 XA_STRING, &type, &format, &len, &left,
                                 &command );
    if( (result!=Success) || (type!=XA_STRING) || (format!=8) )
    {   if( command ) XFree(command);
        return;
    }

    result = 0;
    ptr = (char*)command;
    while( *ptr )
    {   if( *ptr=='C' )
        {   sscanf(ptr+1,"%x %x\n",&source,&serial);
            while( *ptr && (*ptr!='|') ) ptr++;
            if( *ptr=='|' )
            {   result = ExecuteIPCCommand(ptr+1);
            } else result = 0;

            sprintf(buffer,"R %x 0 %d",serial,result);
            handler = XSetErrorHandler( HandleIPCError );
            XChangeProperty( dpy, source, CommAtom, XA_STRING, 8,
                             PropModeAppend, buffer, strlen(buffer)+1 );
            XSync(dpy,False);
            XSetErrorHandler(handler);
        } 

        /* Next Command! */
        while( *ptr++ );
    }
    XFree(command);

    if( result==2 )
        RasMolExit();
}


static int CropRange( val, min, max )
    int val, min, max;
{
    if( val<min ) return( min );
    if( val>max ) return( max );
    return( val );
}


static void ClampDial( dial, value )
    int dial;  Real value;
{
    register Real temp;

    temp = DialValue[dial] + value;

    if( temp > 1.0 )
    {   DialValue[dial] = 1.0;
    } else if( temp < -1.0 )
    {   DialValue[dial] = -1.0;
    } else DialValue[dial] = temp;
}


static void WrapDial( dial, value )
    int dial;  Real value;
{
    register Real temp;

    temp = DialValue[dial] + value;
    while( temp < -1.0 )  temp += 2.0;
    while( temp > 1.0 )   temp -= 2.0;
    DialValue[dial] = temp;
}


void SetMouseMode( mode )
{
    if( mode==MouseMode )
        return;

    if( (mode==MMQuanta) || (MouseMode==MMQuanta) )
    {   /* Enable/Disable Pointer Motion Events! */
        attr.event_mask = ExposureMask | ButtonPressMask | ButtonMotionMask 
                        | ButtonReleaseMask;
        if( mode==MMQuanta ) attr.event_mask |= PointerMotionMask;
        XChangeWindowAttributes( dpy, CanvWin, CWEventMask, &attr );
    }
    MouseMode = mode;
}


static void MouseMove( status, dx, dy )
    int status, dx, dy;
{
    register int index;

    if( MouseMode == MMRasMol )
    {   if( status & ShiftMask )
        {   if( status & Button1Mask ) 
            {   if( dy ) /* Zoom Vertical */
                {   ClampDial( 3, (Real)dy/HRange );
                    ReDrawFlag |= RFZoom;
                }
            } else if( status & (Button2Mask|Button3Mask) )
                if( dx ) /* Z Rotation Horizontal */
                {   WrapDial( 2, (Real)dx/WRange );
                    ReDrawFlag |= RFRotateZ;
                }
        } else if( status & ControlMask )
        {   if( status & Button1Mask )
            {   if( dy ) /* Slab Vertical */
                {   ClampDial( 7, (Real)dy/YRange );
                    ReDrawFlag |= RFSlab;
                }
            }

        } else /* Unmodified! */
            if( status & Button1Mask )
            {   if( dx ) /* Rotate Y Horizontal */
                {   WrapDial( 1, (Real)dx/WRange );
                    index = (DialValue[1]+1.0)*(XRange-48);
                    NewScrlX = (index>>1)+16;
                    ReDrawFlag |= RFRotateY;
                }

                if( dy ) /* Rotate X Vertical */
                {   WrapDial( 0, (Real)dy/HRange );
                    index = (DialValue[0]+1.0)*(YRange-48);
                    NewScrlY = (index>>1)+16;
                    ReDrawFlag |= RFRotateX;
                }
            } else if( status & (Button2Mask|Button3Mask) )
            {   if( dx ) /* Translate X Horizontal */
                {   ClampDial( 4, (Real)dx/XRange );
                    ReDrawFlag |= RFTransX;
                }

                if( dy ) /* Translate Y Vertical */
                {   ClampDial( 5, (Real)dy/YRange );
                    ReDrawFlag |= RFTransY;
                }
            }
    } else if( MouseMode==MMQuanta )
    {   if( status & ShiftMask )
        {   if( status & Button1Mask )
            {   if( dy ) /* Slab Vertical */
                {   ClampDial( 7, (Real)dy/YRange );
                    ReDrawFlag |= RFSlab;
                }
            } else if( status & Button2Mask )
            {   if( dx ) /* Translate X Horizontal */
                {   ClampDial( 4, (Real)dx/XRange );
                    ReDrawFlag |= RFTransX;
                }

                if( dy ) /* Translate Y Vertical */
                {   ClampDial( 5, (Real)dy/YRange );
                    ReDrawFlag |= RFTransY;
                }
            } else if( !(status & Button3Mask) )
                if( dy ) /* Zoom Vertical */
                {   ClampDial( 3, (Real)dy/HRange );
                    ReDrawFlag |= RFZoom;
                }
        } else if( status & Button2Mask )
        {   if( dx ) /* Rotate Y Horizontal */
            {   WrapDial( 1, (Real)dx/WRange );
                index = (DialValue[1]+1.0)*(XRange-48);
                NewScrlX = (index>>1)+16;
                ReDrawFlag |= RFRotateY;
            }

            if( dy ) /* Rotate X Vertical */
            {   WrapDial( 0, (Real)dy/HRange );
                index = (DialValue[0]+1.0)*(YRange-48);
                NewScrlY = (index>>1)+16;
                ReDrawFlag |= RFRotateX;
            }
        } else if( status & Button3Mask )
            if( dx ) /* Z Rotation Horizontal */
            {   WrapDial( 2, (Real)dx/WRange );
                ReDrawFlag |= RFRotateZ;
            }
    } else /* MMInsight */
        switch( status & (Button1Mask|Button2Mask|Button3Mask) )
        {   case( Button1Mask ):
                    if( dx ) /* Rotate Y Horizontal */
                    {   WrapDial( 1, (Real)dx/WRange );
                        index = (DialValue[1]+1.0)*(XRange-48);
                        NewScrlX = (index>>1)+16;
                        ReDrawFlag |= RFRotateY;
                    }

                    if( dy ) /* Rotate X Vertical */
                    {   WrapDial( 0, (Real)dy/HRange );
                        index = (DialValue[0]+1.0)*(YRange-48);
                        NewScrlY = (index>>1)+16;
                        ReDrawFlag |= RFRotateX;
                    }
                    break;

            case( Button2Mask ):
                    if( dx ) /* Translate X Horizontal */
                    {   ClampDial( 4, (Real)dx/XRange );
                        ReDrawFlag |= RFTransX;
                    }

                    if( dy ) /* Translate Y Vertical */
                    {   ClampDial( 5, (Real)dy/YRange );
                        ReDrawFlag |= RFTransY;
                    }
                    break;

            case( Button1Mask|Button2Mask ):
                    ClampDial( 3, (Real)dx/WRange - (Real)dy/HRange );
                    ReDrawFlag |= RFZoom;
                    break;

            case( Button1Mask|Button3Mask ):
                    WrapDial( 2, (Real)dx/WRange - (Real)dy/HRange );
                    ReDrawFlag |= RFRotateZ;
                    break;

            case( Button1Mask|Button2Mask|Button3Mask ):
                    ClampDial( 7, (Real)dx/XRange - (Real)dy/YRange );
                    ReDrawFlag |= RFSlab;
                    break;
        }
}


int FetchEvent( wait )
    int wait;
{ 
    static XEvent event;
    register Real temp;
    register int index;


    NewScrlX = ScrlX;
    NewScrlY = ScrlY;

    if( HeldButton != -1 ) 
        wait = False;

    while( XPending(dpy) || (wait && !ReDrawFlag) )
    {   XNextEvent(dpy,&event);
        switch( event.type )
        {   case(ButtonPress):
                {   XButtonPressedEvent *ptr;

                    HeldButton = -1;
                    ptr = (XButtonPressedEvent*)&event;
                    if( !MenuDisable )
                        for( index=0; index<OptCount; index++ )
                            if( ptr->window==OptWin[index] )
                            {   ReDrawButton( index, False );
                                break;
                            }

                    if( ptr->window==XScrlWin )
                    {   wait = False;
                        if( ptr->x<16 )
                        {   HeldButton = XScrlDial;
                            HeldStep = -XScrlSkip;
                        } else if( ptr->x>=XRange-16 )
                        {   HeldButton = XScrlDial;
                            HeldStep = XScrlSkip;
                        } else
                        {   index = ptr->x-8;
                            if( XScrlDial<3 )
                            {   if( index>XRange-32 ) index -= XRange-48;
                                else if( index<16 ) index += XRange-48;
                                NewScrlX = index;
                            } else NewScrlX = CropRange(index,16,XRange-32);
                        }

                    } else if( ptr->window==YScrlWin )
                    {   wait = False;
                        if( ptr->y<16 )
                        {   HeldButton = YScrlDial;
                            HeldStep = -YScrlSkip;
                        } else if( ptr->y>=YRange-16 )
                        {   HeldButton = YScrlDial;
                            HeldStep = YScrlSkip;
                        } else
                        {   index = ptr->y-8;
                            if( YScrlDial<3 )
                            {   if( index>YRange-32 ) index -= YRange-48;
                                else if( index<16 ) index += YRange-48;
                                NewScrlY = index;
                            } else NewScrlY = CropRange(index,16,YRange-32);
                        }

                    } else if( ptr->window==CanvWin )
                    {   InitX = PointX = ptr->x;
                        InitY = PointY = ptr->y;
                    }
                } break;

            case(MotionNotify):
                {   XMotionEvent *ptr;
                    int dx, dy;

                    ptr = (XMotionEvent*)&event;
                    if( ptr->window==CanvWin )
                    {   if( !IsClose(ptr->x,InitX) || !IsClose(ptr->y,InitY) )
                        {   dx = ptr->x-PointX;  dy = ptr->y-PointY;
                            MouseMove( ptr->state, dx, dy );

                            PointX = ptr->x;
                            PointY = ptr->y;
                        }
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

		    HeldButton = -1;
                    ptr = (XButtonReleasedEvent*)&event;
                    if( ptr->window==CanvWin )
                    {   PointX = ptr->x;  PointY = ptr->y;
                        if( IsClose(PointX,InitX) && IsClose(PointY,InitY) )
                            ReDrawFlag |= RFPoint;
                    } else if( !MenuDisable )
                        for( index=0; index<OptCount; index++ )
                            if( ptr->window==OptWin[index] )
                                return(index-ButMax);

                } break;

            case(KeyPress):
                {   XKeyPressedEvent *ptr;
                    static KeySym symbol;
                    static char keychar;

                    ptr = (XKeyPressedEvent*)&event;
                    index = XLookupString(ptr,&keychar,1,&symbol,NULL);
                    switch( symbol )
                    {   case(XK_Begin):
                        case(XK_Home):  return( 0x01 );  break;
                        case(XK_Right): return( 0x06 );  break;
                        case(XK_Left):  return( 0x02 );  break;
                        case(XK_End):   return( 0x05 );  break;
                        case(XK_Up):
                        case(XK_Prior): return( 0x10 );  break;
                        case(XK_Down):
                        case(XK_Next):  return( 0x0e );  break;

                        default:        if( index==1 ) 
                                            return( keychar );
                    }
                } break;


            case(Expose):
                {   XExposeEvent *ptr;

                    ptr = (XExposeEvent*)&event;
                    if( ptr->window==CanvWin )
                    {   if( image ) {
#ifdef MITSHM
                            if( SharedMemFlag )
                            {   XShmPutImage( dpy, CanvWin, gcon, image,
                                              ptr->x, ptr->y, ptr->x, ptr->y,
                                              ptr->width, ptr->height, False);
                            } else
#endif 
                            XPutImage( dpy, CanvWin, gcon, image,
                                       ptr->x, ptr->y, ptr->x, ptr->y,
                                       ptr->width, ptr->height );
                        }
                    } else if( ptr->window==MainWin )
                    {   ReDrawMain();
                    } else if( ptr->window==XScrlWin )
                    {   ReDrawXScroll();
                    } else if( ptr->window==YScrlWin )
                        ReDrawYScroll();
                    XSync(dpy,False);
                } break;

            case(EnterNotify):
                if( LocalMap )
                    XInstallColormap(dpy,lmap);
                break;

            case(LeaveNotify):
                if( LocalMap )
                    XUninstallColormap(dpy,lmap);
                break;

            case(ConfigureNotify):
                {   XConfigureEvent *ptr;
                    register int wide,high;

                    ptr = (XConfigureEvent*)&event;
                    wide = CropRange(ptr->width, MinWidth, MaxWidth );
                    high = CropRange(ptr->height,MinHeight,MaxHeight);
                    if( (wide!=MainWide) || (high!=MainHigh) )
                        ReSizeWindow(wide,high);
                     
                } break;

            case(ClientMessage):
                {   XClientMessageEvent *ptr;

                    ptr = (XClientMessageEvent*)&event;
                    if( (ptr->message_type==ProtoXAtom) && 
                        (ptr->data.l[0]==DelWinXAtom) )
                        RasMolExit();
                } break;

            case(PropertyNotify):
                {   XPropertyEvent *ptr;

                    ptr = (XPropertyEvent*)&event;
                    if( (ptr->atom==CommAtom) &&
                        (ptr->state==PropertyNewValue) )
                        HandleIPCCommand();
                } break;

            case(MapNotify):
                ReDrawXScroll();
                ReDrawYScroll();
                ReDrawMain();
		break;

            default:  
#ifdef DIALBOX
                if( event.type==DialEvent )
		    HandleDialEvent( &event );
#endif
                break;
        }
    }

    if( HeldButton == YScrlDial )
    {   index = NewScrlY+HeldStep;
        if( YScrlDial<3 )
        {   if( index<16 )             
            {   index += YRange-48;
            } else if( index>YRange-32 ) 
                index -= YRange-48;
            NewScrlY = index;
        } else NewScrlY = CropRange(index,16,YRange-32);
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
        if( XScrlDial<3 )
        {   if( index<16 ) 
            {   index += XRange-48;
            } else if( index>XRange-32 ) 
                index -= XRange-48;
            NewScrlX = index;
        } else NewScrlX = CropRange(index,16,XRange-32);
    }

    if( NewScrlX != ScrlX )
    {   XClearArea(dpy,XScrlWin,ScrlX,0,16,16,False);
        XCopyArea(dpy,Scrl,XScrlWin,gcon,0,0,16,16,NewScrlX,0);

        temp = ((Real)(NewScrlX-16))/(XRange-48);
        DialValue[XScrlDial] = 2.0*temp - 1.0;
        ReDrawFlag |= (1<<XScrlDial);
        ScrlX = NewScrlX;
    }
    XSync(dpy,False);
    return( 0 );
}


int LookUpColour( name, red, grn, blu )
    char *name; int *red, *grn, *blu;
{
    static XColor exact, close;
    register Colormap map;

    map = (LocalMap)? lmap : cmap;
    if( XLookupColor(dpy,map,name,&exact,&close) )
    {   *red = exact.red>>8;
        *grn = exact.green>>8;
        *blu = exact.blue>>8;
        return(True);
    } else 
        return(False);
}


void BeginWait()
{
    if( UseHourGlass )
    {   XDefineCursor(dpy,CanvWin,hglass);
        XDefineCursor(dpy,MainWin,hglass);
        XFlush(dpy);
    }
}


void EndWait()
{
    if( UseHourGlass )
    {   XDefineCursor(dpy,CanvWin,cross);
        XDefineCursor(dpy,MainWin,arrow);
        XFlush(dpy);
    }
}


void CloseDisplay()
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
    if( UseDials && UseDialLEDs )
        for( num=0; num<8; num++ )
            SetDialLabel(num,"");
#endif
    XCloseDisplay( dpy );
}

