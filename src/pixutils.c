/* pixutils.c
 * RasMol2 Molecular Graphics
 * Roger Sayle, February 1994
 * Version 2.3
 */
#include "rasmol.h"

#ifdef IBMPC
#include <windows.h>
#include <malloc.h>
#endif
#ifdef sun386
#include <stdlib.h>
#endif

#include <stdio.h>
#include <math.h>

#define PIXUTILS
#include "pixutils.h"
#include "graphics.h"
#include "transfor.h"
#include "render.h"

/* Sutherland-Cohen Line Clipping Macros */
#define BitAbove    0x01
#define BitBelow    0x02
#define BitRight    0x04
#define BitLeft     0x08
#define BitFront    0x10

#define Reject(x,y)   ((x)&(y))
#define Accept(x,y)   (!((x)|(y)))


#define ZValid(z)     ((!UseSlabPlane) || ((z)<SlabValue))
#define XValid(x)     (((x)>=0)&&((x)<XRange))
#define YValid(y)     (((y)>=0)&&((y)<YRange))
#define RootSix       2.44948974278

#define SETPIXEL(dptr,fptr,d,c)    if( (d) > *(dptr) )              \
                                   {   *(dptr) = (d);               \
                                       *(fptr) = (c);               \
                                   }


typedef struct {
                int x, y, z;
                int inten;
               } Vertex;

typedef struct {
                short dx,dy,dz;
                short inten;
		Long offset;
	       } ArcEntry;


/* Note: DrawCylinderCaps currently employs an
 *       extremely crude hack to avoid stripes
 *       appearing along cylinders.
 */
#define ARCSIZE  2048

static ArcEntry __far *ArcAcPtr;
static ArcEntry __far *ArcDnPtr;
#ifdef IBMPC
static ArcEntry __far *ArcAc;
static ArcEntry __far *ArcDn;
#else
static ArcEntry ArcAc[ARCSIZE];
static ArcEntry ArcDn[ARCSIZE];
#endif


static int ClipStatus;


static int OutCode(x,y,z)
    register int x,y,z;
{
    register int result;

    result=0;
    if( y<0 )
    {   result |= BitAbove;
    } else if( y>=YRange )
        result |= BitBelow;

    if( x<0 )
    {   result |= BitLeft;
    } else if( x>=XRange )
        result |= BitRight;

    if( !ZValid(z) )
        result |= BitFront;
    return result;
}


#ifndef PIXUTILS  /* Unused Function */
void PlotPoint(x,y,z,col)
    int x,y,z,col;
{
    register Pixel __huge *fptr;
    register short __huge *dptr;
    register Long offset;

    /* SETPIXEL(dptr,fptr,z,Lut[col+ColourMask]); */

    offset = (Long)y*XRange+x;
    dptr = DBuffer+offset;
    if( z > *dptr )
    {   fptr = FBuffer+offset;
        *fptr = Lut[col+ColourMask];
        *dptr = z;
    }
}
#endif

#ifndef PIXUTILS  /* Unused Function */
void ClipPoint(x,y,z,col)
    int x,y,z,col;
{
    register Pixel __huge *fptr;
    register short __huge *dptr;
    register Long offset;

    if( XValid(x) && YValid(y) && ZValid(z) )
    {   /* PlotPoint(x,y,z,col); */
        offset = (Long)y*XRange+x;
        dptr = DBuffer+offset;
        if( z > *dptr )
        {   fptr = FBuffer+offset;
            *fptr = Lut[col+ColourMask];
            *dptr = z;
        }
    }
}
#endif

#ifndef PIXUTILS  /* Unused Function */
void PlotDeepPoint(x,y,z,col)
    int x,y,z,col;
{
    register Long offset;
    register Pixel __huge *fptr;
    register short __huge *dptr;
    register int inten;

    offset = (Long)y*XRange+x;
    dptr = DBuffer+offset;

    if( z > *dptr )
    {  fptr = FBuffer+offset;
       inten = (int)((Long)(ColourDepth*z)/ImageSize);
       *fptr = Lut[col+inten];
       *dptr = z;
    }
}
#endif

#ifndef PIXUTILS  /* Unused Function */
void ClipDeepPoint(x,y,z,col)
    int x,y,z,col;
{
    register Long offset;
    register Pixel __huge *fptr;
    register short __huge *dptr;
    register int inten;

    if( XValid(x) && YValid(y) && ZValid(z) )
    {   /* PlotDeepPoint(x,y,z,col); */
        offset = (Long)y*XRange+x;
        dptr = DBuffer+offset;

        if( z > *dptr )
        {  fptr = FBuffer+offset;
           inten = (int)((Long)(ColourDepth*z)/ImageSize);
           *fptr = Lut[col+inten];
           *dptr = z;
        }
    }
}
#endif


/* Macros for Bresenhams Line Drawing Algorithm */
#define CommonStep(s)  z1 += zrate; SETPIXEL(dptr,fptr,z1,c);     \
                       if( (zerr+=dz)>0 ) { zerr-=(s); z1+=iz; }

#define XStep  { if((err+=dy)>0) { fptr+=ystep; dptr+=ystep; err-=dx; } \
                 fptr+=ix; dptr+=ix; x1+=ix; CommonStep(dx); }

#define YStep  { if((err+=dx)>0) { fptr+=ix; dptr+=ix; err-=dy; } \
                 fptr+=ystep; dptr+=ystep; y1+=iy; CommonStep(dy); }
                     

void DrawTwinLine(x1,y1,z1,x2,y2,z2,col1,col2)
    int x1,y1,z1,x2,y2,z2,col1,col2;
{
    register Long offset;
    register Pixel __huge *fptr;
    register short __huge *dptr;
    register int zrate, zerr;
    register int ystep,err;
    register int ix,iy,iz;
    register int dx,dy,dz;
    register int mid;
    register Pixel c;

    c = Lut[col1+ColourMask];

    offset = (Long)y1*XRange + x1;
    fptr = FBuffer+offset;
    dptr = DBuffer+offset;

    SETPIXEL(dptr,fptr,z1,c);

    dx = x2-x1;  dy = y2-y1; dz = z2-z1;
    if( !dx && !dy ) return;

    ystep = XRange;
    ix = iy = iz = 1;
    if( dy<0 ) { dy = -dy; iy = -1; ystep = -ystep; }
    if( dx<0 ) { dx = -dx; ix = -1; }
    if( dz<0 ) { dz = -dz; iz = -1; }

    if( dx>dy )
    {   zrate = dz/dx;  
        dz -= dx*zrate;
        err = zerr = -(dx>>1);

        if( col1 != col2 )
        {   mid = (x1+x2)>>1;
            while( x1!=mid ) XStep;
            c = Lut[col2+ColourMask];
        }
        while( x1!=x2 ) XStep;
    } else
    {   zrate = dz/dy;
        dz -= dy*zrate;
        err = zerr = -(dy>>1);

        if( col1 != col2 )
        {   mid = (y1+y2)>>1;
            while( y1!=mid ) YStep;
            c = Lut[col2+ColourMask];
        }
        while( y1!=y2 ) YStep;
    }
}


static void ClipLine(x1,y1,z1,x2,y2,z2,col)
    int x1,y1,z1,x2,y2,z2,col;
{
    register int code1,code2;
    register int delta,rest;
    register int temp;

    while( True )
    {   code1 = OutCode(x1,y1,z2);
        code2 = OutCode(x2,y2,z2);
        if( Reject(code1,code2) ) return;
        if( Accept(code1,code2) ) break;

        if( !code1 )
        {   temp=x1; x1=x2; x2=temp;
            temp=y1; y1=y2; y2=temp;
            temp=z1; z1=z2; z2=temp;
            code1 = code2;
        }

        if( code1 & BitAbove )
        {   delta = y2-y1;
            x1 += (y1*(x1-x2))/delta;  
            z1 += (y1*(z1-z2))/delta;
            y1 = 0;
        } else if( code1 & BitLeft )
        {   delta = x2-x1;
            y1 += (x1*(y1-y2))/delta;
            z1 += (x1*(z1-z2))/delta;
            x1 = 0;
        } else if( code1 & BitRight )
        {   delta = x2-x1;
            temp=XRange-1; rest=temp-x1;
            y1 += (rest*(y2-y1))/delta;
            z1 += (rest*(z2-z1))/delta;
            x1 = temp;
        } else if( code1 & BitBelow )
        {   delta = y2-y1;
            temp=YRange-1; rest=temp-y1;
            y1 += (rest*(x2-x1))/delta;
            z1 += (rest*(z2-z1))/delta;
            y1 = temp;
        } else /* SLAB */
        {   delta = z2-z1;
            rest = SlabValue-z1;
            x1+=(rest*(x2-x1))/delta;
            y1+=(rest*(y2-y1))/delta;
            z1 = SlabValue;
        }
    }
    DrawTwinLine(x1,y1,z1,x2,y2,z2,col,col);
}


void ClipTwinLine(x1,y1,z1,x2,y2,z2,col1,col2)
    int x1,y1,z1,x2,y2,z2,col1,col2;
{
    register int xmid,ymid,zmid;
    register int code1,code2;


    if( col1!=col2 )
    {   code1 = OutCode(x1,y1,z1);
        code2 = OutCode(x2,y2,z2);
        if( !Reject(code1,code2) )
        {   if( !Accept(code1,code2) )
            {  xmid = (x1+x2)/2;
               ymid = (y1+y2)/2;
               zmid = (z1+z2)/2;
               ClipLine(x1,y1,z1,xmid,ymid,zmid,col1);
               ClipLine(xmid,ymid,zmid,x2,y2,z2,col2);
            } else
               DrawTwinLine(x1,y1,z1,x2,y2,z2,col1,col2);
        }
    } else
        ClipLine(x1,y1,z1,x2,y2,z2,col1);
}



/* Macros for 3D Bresenhams Vector Algorithm */
#define CommonVectStep(s)  z1 += zrate;   c1 += crate;                    \
                           SETPIXEL(dptr,fptr,z1,Lut[col+c1]);            \
                           if( (zerr+=dz)>0 ) { zerr -= (s); z1 += iz; }  \
                           if( (cerr+=dc)>0 ) { cerr -= (s); c1 += ic; }

#define XVectStep  { if((err+=dy)>0) { fptr+=ystep; dptr+=ystep; err-=dx; } \
                     fptr+=ix; dptr+=ix; x1+=ix; CommonVectStep(dx); }

#define YVectStep  { if((err+=dx)>0) { fptr+=ix; dptr+=ix; err-=dy; } \
                     fptr+=ystep; dptr+=ystep; y1+=iy; CommonVectStep(dy); }


void DrawTwinVector(x1,y1,z1,x2,y2,z2,col1,col2)
    int x1,y1,z1,x2,y2,z2,col1,col2;
{
    register Long offset;
    register Pixel __huge *fptr;
    register short __huge *dptr;
    register int ix,iy,iz,ic;
    register int dx,dy,dz,dc;
    register int crate, cerr;
    register int zrate, zerr;
    register int ystep,err;
    register int col, mid;
    register int c1, c2;

    c1 = (ColourDepth*z1)/ImageSize;
    c2 = (ColourDepth*z2)/ImageSize;

    offset = (Long)y1*XRange + x1;
    fptr = FBuffer+offset;
    dptr = DBuffer+offset;

    SETPIXEL(dptr,fptr,z1,Lut[col1+c1]);

    dx = x2 - x1;  dy = y2 - y1;
    dz = z2 - z1;  dc = c2 - c1;
    if( !dx && !dy ) return;

    ystep = XRange;
    ix = iy = iz = ic = 1;
    if( dy<0 ) { dy = -dy; iy = -1; ystep = -ystep; }
    if( dx<0 ) { dx = -dx; ix = -1; }
    if( dz<0 ) { dz = -dz; iz = -1; }
    if( dc<0 ) { dc = -dc; ic = -1; }

    if( dx>dy )
    {   zrate = dz/dx;  dz -= dx*zrate;
        crate = dc/dx;  dc -= dx*crate;
        err = zerr = cerr = -(dx>>1);
        col = col1;

        if( col1 != col2 )
        {   mid = (x1+x2)>>1;
            while( x1!=mid ) XVectStep;
            col = col2;
        }
        while( x1!=x2 ) XVectStep;
    } else
    {   zrate = dz/dy;  dz -= dy*zrate;
        crate = dc/dy;  dc -= dy*crate;
        err = zerr = cerr = -(dy>>1);
        col = col1;

        if( col1 != col2 )
        {   mid = (y1+y2)>>1;
            while( y1!=mid ) YVectStep;
            col=col2;
        }
        while( y1!=y2 ) YVectStep;
    }
}

static void ClipVector(x1,y1,z1,x2,y2,z2,col)
    int x1,y1,z1,x2,y2,z2,col;
{
    register int code1,code2;
    register int delta,rest;
    register int temp;

    while( True )
    {   code1 = OutCode(x1,y1,z2);
        code2 = OutCode(x2,y2,z2);
        if( Reject(code1,code2) ) return;
        if( Accept(code1,code2) ) break;

        if( !code1 )
        {   temp=x1; x1=x2; x2=temp;
            temp=y1; y1=y2; y2=temp;
            temp=z1; z1=z2; z2=temp;
            code1 = code2;
        }

        if( code1 & BitAbove )
        {   delta = y2-y1;
            x1 += (y1*(x1-x2))/delta;  
            z1 += (y1*(z1-z2))/delta;
            y1 = 0;
        } else if( code1 & BitLeft )
        {   delta = x2-x1;
            y1 += (x1*(y1-y2))/delta;
            z1 += (x1*(z1-z2))/delta;
            x1 = 0;
        } else if( code1 & BitRight )
        {   delta = x2-x1;
            temp=XRange-1; rest=temp-x1;
            y1 += (rest*(y2-y1))/delta;
            z1 += (rest*(z2-z1))/delta;
            x1 = temp;
        } else if( code1 & BitBelow )
        {   delta = y2-y1;
            temp=YRange-1; rest=temp-y1;
            y1 += (rest*(x2-x1))/delta;
            z1 += (rest*(z2-z1))/delta;
            y1 = temp;
        } else /* SLAB */
        {   delta = z2-z1;
            rest = SlabValue-z1;
            x1+=(rest*(x2-x1))/delta;
            y1+=(rest*(y2-y1))/delta;
            z1 = SlabValue;
        }
    }
    DrawTwinVector(x1,y1,z1,x2,y2,z2,col,col);
}


void ClipTwinVector(x1,y1,z1,x2,y2,z2,col1,col2)
    int x1,y1,z1,x2,y2,z2,col1,col2;
{
    register int xmid,ymid,zmid;
    register int code1,code2;

    if( col1!=col2 )
    {   code1 = OutCode(x1,y1,z1);
        code2 = OutCode(x2,y2,z2);
        if( !Reject(code1,code2) )
        {   if( !Accept(code1,code2) )
            {  xmid = (x1+x2)/2;
               ymid = (y1+y2)/2;
               zmid = (z1+z2)/2;
               ClipVector(x1,y1,z1,xmid,ymid,zmid,col1);
               ClipVector(xmid,ymid,zmid,x2,y2,z2,col2);
            } else
               DrawTwinVector(x1,y1,z1,x2,y2,z2,col1,col2);
        }
    } else
        ClipVector(x1,y1,z1,x2,y2,z2,col1);
}


void ClipDashVector(x1,y1,z1,x2,y2,z2,col1,col2)
    int x1,y1,z1,x2,y2,z2,col1,col2;
{
    register Long offset;
    register Pixel __huge *fptr;
    register short __huge *dptr;
    register int ix,iy,iz,ic;
    register int dx,dy,dz,dc;
    register int crate, cerr;
    register int zrate, zerr;
    register int ystep,err;
    register int col, mid;
    register int c1, c2;


    if( (x1==x2) && (y1==y2) ) return;
    if( Reject(OutCode(x1,y1,z1),OutCode(x2,y2,z2)) )
        return;

    c1 = (ColourDepth*z1)/ImageSize;
    c2 = (ColourDepth*z2)/ImageSize;

    dz = (z2 - z1)<<1;  
    dc = (c2 - c1)<<1;
    dx = x2 - x1;  
    dy = y2 - y1;

    offset = (Long)y1*XRange + x1;
    fptr = FBuffer+offset;
    dptr = DBuffer+offset;


    ystep = XRange;
    ix = iy = iz = ic = 1;
    if( dy<0 ) { dy = -dy; iy = -1; ystep = -ystep; }
    if( dx<0 ) { dx = -dx; ix = -1; }
    if( dz<0 ) { dz = -dz; iz = -1; }
    if( dc<0 ) { dc = -dc; ic = -1; }


    if( dx>dy )
    {   if( x2<x1 )
        {   mid = col1;
            col1 = col2;
            col2 = mid;
        }
        zrate = dz/dx;  dz -= dx*zrate;
        crate = dc/dx;  dc -= dx*crate;
        err = zerr = cerr = -(dx>>1);
        mid = (x1+x2)/2;

        while( x1!=x2 )
        {   if( XValid(x1) && YValid(y1) && ZValid(z1) )
            {   col = (x1<mid)? col1 : col2;
                SETPIXEL(dptr,fptr,z1,Lut[col+c1]);
            }

            if( (x1+ix)==x2 ) 
                break;

            if( (err+=dy+dy)>0 )
                if( (err-=dx)>0 )
                {   err-=dx;
                    fptr+=ystep+ystep;
                    dptr+=ystep+ystep;
                    y1+=iy+iy;
                } else
                {   fptr+=ystep;
                    dptr+=ystep;
                    y1+=iy;
                }

            if( (zerr+=dz)>0 )
                if( (zerr-=dx)>0 )
                {   zerr -= dx;
                    z1 += iz+iz;
                } else z1 += iz;

            if( (cerr+=dc)>0 )
                if( (cerr-=dx)>0 )
                {   cerr -= dx;
                    c1 += ic+ic;
                } else c1 += ic;

            fptr+=ix+ix; dptr+=ix+ix; x1+=ix+ix;
            z1 += zrate;   c1 += crate;
        }
    } else
    {   if( y1>y2 )
        {   mid = col1;
            col1 = col2;
            col2 = mid;
        }

        zrate = dz/dy;  dz -= dy*zrate;
        crate = dc/dy;  dc -= dy*crate;
        err = zerr = cerr = -(dy>>1);
        mid = (y1+y2)/2;

        
        while( y1!=y2 )
        {   if( XValid(x1) && YValid(y1) && ZValid(z1) )
            {   col = (y1<mid)? col1 : col2;
                SETPIXEL(dptr,fptr,z1,Lut[col+c1]);
            }

            if( (y1+iy)==y2 ) 
                break;

            if( (err+=dx+dx)>0 )
                if( (err-=dy)>0 )
                {   err-=dy;
                    fptr+=ix+ix;
                    dptr+=ix+ix;
                    x1+=ix+ix;
                } else
                {   fptr+=ix;
                    dptr+=ix;
                    x1+=ix;
                }

            if( (zerr+=dz)>0 )
                if( (zerr-=dy)>0 )
                {   zerr -= dy;
                    z1 += iz+iz;
                } else z1 += iz;

            if( (cerr+=dc)>0 )
                if( (cerr-=dy)>0 )
                {   cerr -= dy;
                    c1 += ic+ic;
                } else c1 += ic;

            fptr+=ystep+ystep; dptr+=ystep+ystep; y1+=iy+iy;
            z1 += zrate;   c1 += crate;
        }
    }
}


/* SplineCount is either 1, 2, 3, 4, 5 or 9! */

void StrandRibbon( src, dst, col )
    Knot *src, *dst;  int col;
{
    register int hsx, hsy, hsz;
    register int hdx, hdy, hdz;
    register int qsx, qsy, qsz;
    register int qdx, qdy, qdz;

    if( SplineCount != 4 )
    {   if( SplineCount != 2 )
        {   ClipVector( src->px, src->py, src->pz,
                        dst->px, dst->py, dst->pz, col );
            if( SplineCount==1 ) return;
        }

        ClipVector( src->px+src->wx, src->py+src->wy, src->pz+src->wz,
                    dst->px+dst->wx, dst->py+dst->wy, dst->pz+dst->wz, col );
        ClipVector( src->px-src->wx, src->py-src->wy, src->pz-src->wz,
                    dst->px-dst->wx, dst->py-dst->wy, dst->pz-dst->wz, col );
        if( SplineCount<=3 ) return;

        hsx = src->wx/2;  hsy = src->wy/2;  hsz = src->wz/2;
        hdx = dst->wx/2;  hdy = dst->wy/2;  hdz = dst->wz/2;

        ClipVector( src->px+hsx, src->py+hsy, src->pz+hsz,
                    dst->px+hdx, dst->py+hdy, dst->pz+hdz, col );
        ClipVector( src->px-hsx, src->py-hsy, src->pz-hsz,
                    dst->px-hdx, dst->py-hdy, dst->pz-hdz, col );
        if( SplineCount==5 ) return;
    } else /* SplineCount == 4 */
    {   hsx = src->wx/2;  hsy = src->wy/2;  hsz = src->wz/2;
        hdx = dst->wx/2;  hdy = dst->wy/2;  hdz = dst->wz/2;
    }

    qsx = hsx/2;  qsy = hsy/2;  qsz = hsz/2;
    qdx = hdx/2;  qdy = hdy/2;  qdz = hdz/2;

    ClipVector( src->px+hsx+qsx, src->py+hsy+qsy, src->pz+hsz+qsz,
                dst->px+hdx+qdx, dst->py+hdy+qdy, dst->pz+hdz+qdz, col );
    ClipVector( src->px+hsx-qsx, src->py+hsy-qsy, src->pz+hsz-qsz,
                dst->px+hdx-qdx, dst->py+hdy-qdy, dst->pz+hdz-qdz, col );
    ClipVector( src->px-hsx+qsx, src->py-hsy+qsy, src->pz-hsz+qsz,
                dst->px-hdx+qdx, dst->py-hdy+qdy, dst->pz-hdz+qdz, col );
    ClipVector( src->px-hsx-qsx, src->py-hsy-qsy, src->pz-hsz-qsz,
                dst->px-hdx-qdx, dst->py-hdy-qdy, dst->pz-hdz-qdz, col );
}


#ifndef PIXUTILS  /* Unused Function */
static void OutLinePolygon( v1, v2, v3 )
    Vertex *v1, *v2, *v3;
{
    ClipLine( v1->x, v1->y, v1->z, v2->x, v2->y, v2->z, v1->inten );
    ClipLine( v2->x, v2->y, v2->z, v3->x, v3->y, v3->z, v2->inten );
    ClipLine( v3->x, v3->y, v3->z, v1->x, v1->y, v1->z, v3->inten );
}
#endif


static void RenderPolygon( v1, v2, v3 )
    Vertex *v1, *v2, *v3;
{
    register Pixel __huge *fptr;
    register short __huge *dptr;
    register Long offset;

    register int startx, stopx;
    register int startz, stopz;
    register Real ratio;

    register Vertex *vtmp;
    register int xmin, xmax;
    register int ymin, ymax;
    register int x,y,z;

    /* Sort increasing y! */
    if( v1->y > v2->y )
    {   vtmp = v2; 
        v2 = v1; 
        v1 = vtmp;
    }

    if( v2->y > v3->y )
    {   if( v1->y > v3->y )
        {   vtmp = v3; 
            v3 = v2;
            v2 = v1; 
            v1 = vtmp;
        } else
        {   vtmp = v3; 
            v3 = v2; 
            v2 = vtmp;
        }
    }

    /* Perform Polygon Y Culling! */
    if( (v3->y < 0) || (v1->y >= YRange) )
        return;

    xmin = MinFun(v1->x,v2->x);
    xmax = MaxFun(v1->x,v2->x);
    if( v3->x < xmin ) 
    {   xmin = v3->x;
    } else if( v3->x > xmax )
        xmax = v3->x;

    if( (xmax < 0) || (xmin >= XRange) )
        return;

    if( UseSlabPlane )  /* Perform Polygon Z Culling! */
    {   if( (v1->z>=SlabValue) && (v2->z>=SlabValue) && (v3->z>=SlabValue) )
            return;
    }


    if( v1->y == v3->y )
        return;

    ymax = MinFun(v3->y,YRange-1);
    ymin = MaxFun(v1->y,0);

    y = ymin;
    offset = (Long)y*XRange;
    while( y <= ymax )
    {   if( y < v2->y )
        {   ratio = (Real)(y - v1->y)/(v2->y - v1->y);
            startx = (int)(ratio*(v2->x - v1->x)) + v1->x;
            startz = (int)(ratio*(v2->z - v1->z)) + v1->z;
        } else if( y > v2->y )
        {   ratio = (Real)(y - v2->y)/(v3->y - v2->y);
            startx = (int)(ratio*(v3->x - v2->x)) + v2->x;
            startz = (int)(ratio*(v3->z - v2->z)) + v2->z;
        } else /* y == v2->y */
        {   startx = v2->x;
            startz = v2->z;
        }

        ratio = (float)(y - v1->y)/(v3->y - v1->y);
        stopx = (int)(ratio*(v3->x - v1->x)) + v1->x;
        stopz = (int)(ratio*(v3->z - v1->z)) + v1->z;

        if( startx < stopx )
        {   xmin = startx;
            xmax = stopx;
        } else
        {   xmax = startx;
            xmin = stopx;
        }

        if( xmax>=XRange )
            xmax = XRange-1;
        if( xmin<0 ) 
            xmin = 0;

        z = (startz + stopz)/2;
        fptr = FBuffer+offset+xmin;
        dptr = DBuffer+offset+xmin;

        for( x=xmin; x<=xmax; x++ )
            if( z > *dptr )
            {   *fptr++ = Lut[v1->inten];
                *dptr++ = z;
            } else 
            {   fptr++;
                dptr++;
            }
        offset += XRange;
        y++;
    }
}


void SolidRibbon( src, dst, col )
    Knot *src, *dst;  int col;
{
    auto Vertex v1, v2, v3;

    v1.x = src->px+src->wx;  
    v1.y = src->py+src->wy;  
    v1.z = src->pz+src->wz;
    v1.inten = src->inten+col;

    v2.x = dst->px+dst->wx;  
    v2.y = dst->py+dst->wy;  
    v2.z = dst->pz+dst->wz;
    v2.inten = dst->inten+col;

    v3.x = dst->px-dst->wx;
    v3.y = dst->py-dst->wy;  
    v3.z = dst->pz-dst->wz;
    v3.inten = dst->inten+col;

    /* OutLinePolygon(&v1,&v2,&v3); */
    RenderPolygon( &v1, &v2, &v3 );

    v2.x = src->px-src->wx;  
    v2.y = src->py-src->wy;  
    v2.z = src->pz-src->wz;
    v2.inten = src->inten+col;

    /* OutLinePolygon(&v1,&v2,&v3); */
    RenderPolygon( &v1, &v2, &v3 );
}



static int TestSphere( x, y, z, rad )
    register int x, y, z, rad;
{
    register int temp;

    ClipStatus = 0;

    if( UseSlabPlane )
    {   if( z-rad>=SlabValue )
            return( False );

        if( z+rad>=SlabValue )
        {   if( SlabMode )
            {   ClipStatus |= BitFront;
            } else return( False );
        } else if( SlabMode==SlabSection )
            return( False );
    }

    temp = y-rad;
    if( temp>=YRange ) return( False );
    if( temp<0 ) ClipStatus |= BitAbove;

    temp = y+rad;
    if( temp<0 ) return( False );
    if( temp>=YRange ) ClipStatus |= BitBelow;

    temp = x+rad;
    if( temp<0 ) return( False );
    if( temp>=XRange ) ClipStatus |= BitRight;

    temp = x-rad;
    if( temp>=XRange ) return( False );
    if( temp<0 ) ClipStatus |= BitLeft;

    return True;
}


#ifdef INVERT
#define CalcInten(dz)    inten = (dz)+(dz)+dx+dy
#else
#define CalcInten(dz)    inten = (dz)+(dz)+dx-dy
#endif

#define UpdateAcross(dz)    \
        depth = (dz)+z;                    \
        if( depth > *dptr )                \
        {   *dptr = depth;                 \
            fptr = fold+dx;                \
            CalcInten((dz));               \
            if( inten>0 )                  \
            {      inten = (int)((inten*ColConst[rad])>>ColBits); \
                   *fptr = Lut[col+inten]; \
            } else *fptr = Lut[col];       \
        }                                  \
        dptr++;  dx++;


#define UpdateLine  \
        dx = -wide;                   \
        dptr = dold-wide;             \
        tptr = LookUp[wide]+wide;     \
        while( dx<0 ) { UpdateAcross(*tptr); tptr--; }       \
        do { UpdateAcross(*tptr); tptr++; } while(dx<=wide); \
        dold += XRange;  fold += XRange;                     \
        dy++;


void DrawSphere(x,y,z,rad,col)
    int x,y,z,rad,col;
{
    register Pixel __huge *fptr, __huge *fold;
    register short __huge *dptr, __huge *dold;
    register char __far *tptr;

    register Long offset;
    register int depth,wide,inten;
    register int dx,dy;

    offset = (Long)(y-rad)*XRange + x;
    fold=FBuffer+offset;  
    dold=DBuffer+offset;

    dy = -rad;
    while( dy<0 ) 
    {   wide = LookUp[rad][-dy]; 
        UpdateLine; 
    }

    do { 
        wide = LookUp[rad][dy];  
        UpdateLine; 
    } while( dy<=rad );
}


void ClipSphere(x,y,z,rad,col)
    int x,y,z,rad,col;
{
    register Pixel __huge *fptr, __huge *fold;
    register short __huge *dptr, __huge *dold;

    register int lastx,lasty,dx,dy,dz;
    register int depth,wide,inten,side;
    register int crad,cwide,temp;
    register Long offset;


    /* Visibility Tests */
    if( !TestSphere(x,y,z,rad) )
        return;

    if( !ClipStatus )
    {   DrawSphere(x,y,z,rad,col);
        return;
    }

    if( ClipStatus&BitAbove )
    {   dy = -y;
        fold = FBuffer + x;
        dold = DBuffer + x;
    } else
    {   dy = -rad;
        offset = (Long)(y+dy)*XRange+x;
        fold = FBuffer + offset;
        dold = DBuffer + offset;
    }

    if( ClipStatus&BitBelow )
    {   lasty = (YRange-1)-y;
    } else lasty = rad;


    side = (XRange-1)-x;
    /* No Slab Plane Clipping */
    if( !(ClipStatus&BitFront) )
    {   while( dy<=lasty )
        {   wide = LookUp[rad][AbsFun(dy)];
            lastx = MinFun(wide,side);
            dx = - MinFun(wide,x);
            dptr = dold + dx;

            while( dx<=lastx )
            {   dz = LookUp[wide][AbsFun(dx)];
                UpdateAcross(dz);
            }
            dold += XRange;
            fold += XRange;
            dy++;
        }
        return;
    }


    dz = SlabValue-z;
    crad = LookUp[rad][AbsFun(dz)];

   
    if( (z>SlabValue) || (SlabMode==SlabSection) )
    {   if( crad<lasty ) lasty = crad;
        if( -crad>dy ) 
        {   dy = -crad;
            offset = (Long)(y+dy)*XRange+x;
            fold = FBuffer + offset;
            dold = DBuffer + offset;
        }
    }

    while( dy<=lasty )
    {   temp = AbsFun(dy);
        wide = LookUp[rad][temp];
        lastx = MinFun(wide,side);
        dx = - MinFun(x,wide);
        dptr = dold + dx;

        if( temp<=crad )
        {   cwide = LookUp[crad][temp];
            while( dx<=lastx )
            {   temp = AbsFun(dx);
                if( temp<=cwide )
                {    /* Slab Plane Clipping Modes */
                    switch( SlabMode )
                    {   case( SlabFinal ):
                                fold[dx] = Lut[col+SlabInten];
                                *dptr = SliceValue;
                                break;

                        case( SlabHollow ):
                                dz = LookUp[wide][temp];
                                depth = z - dz;
                                if( depth>*dptr )
                                {   *dptr = depth;
#ifdef INVERT
                                    inten = dz+dz-dx-dy;
#else
                                    inten = dz+dz-dx+dy;
#endif
                                    if( inten>0 )
                                    {   inten=(int)( (inten*ColConst[rad])
                                                     >>(ColBits+1));
                                        fold[dx] = Lut[col+inten];
                                    } else fold[dx] = Lut[col];
                                }
                                break;

                        case( SlabSection ):
                        case( SlabClose ):
                                dz = SlabValue-z;
                                depth = dx*dx+dy*dy+dz*dz+SliceValue;
                                if( (*dptr<SliceValue) || (depth<*dptr) )
                                {   fold[dx] = Lut[col+SlabInten];
                                    *dptr = depth;
                                }
                                break;
                    }
                    dptr++;  dx++;
                } else if( (z<SlabValue) && (SlabMode!=SlabSection) )
                {    dz = LookUp[wide][temp];
                     UpdateAcross(dz);
                } else
                {   dptr++;  dx++;
                }
            }
        } else /* Slabless ScanLine */
            while( dx<=lastx )
            {   dz = LookUp[wide][AbsFun(dx)];
                UpdateAcross(dz);
            }

        dold += XRange;
        fold += XRange;
        dy++;
    }
}


#if defined(__STDC__) || defined(IBMPC)
/* Function Prototypes */
static void DrawArcDn( short __huge*, Pixel __huge*, int, int );
static void DrawArcAc( short __huge*, Pixel __huge*, int, int );
#endif


static void DrawArcAc(dbase,fbase,z,c)
    register short __huge *dbase;
    register Pixel __huge *fbase;
    register int z,c;
{
    register ArcEntry __far *ptr;
    register short __huge *dptr;
    register short depth;

    for( ptr=ArcAc; ptr<ArcAcPtr; ptr++ )
    {   dptr = dbase+ptr->offset;  depth = ptr->dz+z;
        SETPIXEL(dptr,fbase+ptr->offset,depth,Lut[ptr->inten+c]);
    }
}

static void DrawArcDn(dbase,fbase,z,c)
    register short __huge *dbase;
    register Pixel __huge *fbase;
    register int z,c;
{
    register ArcEntry __far *ptr;
    register short __huge *dptr;
    register short depth;

    for( ptr=ArcDn; ptr<ArcDnPtr; ptr++ )
    {   dptr = dbase+ptr->offset;  depth = ptr->dz+z;
        SETPIXEL(dptr,fbase+ptr->offset,depth,Lut[ptr->inten+c]);
    }
}


static void DrawCylinderCaps( x1,y1,z1,c1,x2,y2,z2,c2,rad )
    int x1,y1,z1,c1, x2,y2,z2,c2, rad;
{
    register short __huge *dold, __huge *dptr;
    register Pixel __huge *fold;

    register Long offset,temp,end;
    register int inten,zrate,absx;
    register int wide,depth;
    register int ix,iy,ax,ay;
    register int lx,ly,lz;
    register int dx,dy,dz;

    lx = x2-x1;
    ly = y2-y1;
    lz = z2-z1;

    if( ly>0 ) { ay = ly; iy = 1; }
    else { ay = -ly; iy = -1; }
    if( lx>0 ) { ax = lx; ix = 1; }
    else { ax = -lx; ix = -1; }
    zrate = lz/MaxFun(ax,ay);


    end = (Long)ly*XRange+lx;
    temp = (Long)y1*XRange+x1;
    fold = FBuffer+temp;
    dold = DBuffer+temp;

    ArcAcPtr = ArcAc;
    ArcDnPtr = ArcDn;

    temp = (Long)-(rad*XRange);
    for( dy= -rad; dy<=rad; dy++ )
    {   wide = LookUp[rad][AbsFun(dy)];

        for( dx= -wide; dx<=wide; dx++ )
        {   absx = AbsFun(dx);
            dz = LookUp[wide][absx];
            CalcInten(dz);
            if( inten>0 )
            {   inten = (int)((inten*ColConst[rad])>>ColBits);
            } else inten = 0;
            offset = temp+dx;

            if( XValid(x1+dx) && YValid(y1+dy) )
            {   dptr = dold+offset; depth = z1+dz;
                SETPIXEL(dptr,fold+offset,depth,Lut[c1+inten]);
            }

            if( XValid(x2+dx) && YValid(y2+dy) )
            {   dptr = dold+(offset+end); depth = z2+dz;
                SETPIXEL(dptr,fold+(offset+end),depth,Lut[c2+inten]);
            }

/*
            k1 = AbsFun(dx+ix); 
            k2 = AbsFun(dx-ix);

            if( ((k1>wide)||(dz>=LookUp[wide][k1]-zrate)) &&
                ((k2>wide)||(dz>LookUp[wide][k2]+zrate)) )
 */
            {   ArcAcPtr->offset = offset; ArcAcPtr->inten = inten;
                ArcAcPtr->dx=dx; ArcAcPtr->dy=dy; ArcAcPtr->dz=dz;
                ArcAcPtr++;
            }

/*
            k1 = AbsFun(dy+iy);
            k2 = AbsFun(dy-iy);

            high = LookUp[rad][absx];
            if( ((k1>high)||(dz>=LookUp[LookUp[rad][k1]][absx]-zrate)) &&
                ((k2>high)||(dz>LookUp[LookUp[rad][k2]][absx]+zrate)) )
 */
            {   ArcDnPtr->offset = offset; ArcDnPtr->inten = inten;
                ArcDnPtr->dx=dx; ArcDnPtr->dy=dy; ArcDnPtr->dz=dz;
                ArcDnPtr++;
            }
        }
        temp += XRange;
    }
}


void DrawCylinder( x1,y1,z1,c1,x2,y2,z2,c2,rad )
    int x1,y1,z1,c1, x2,y2,z2,c2, rad;
{
    register short __huge *dbase;
    register Pixel __huge *fbase;

    register int zrate,zerr,ystep,err;
    register int ix,iy,ax,ay;
    register int lx,ly,lz;
    register int mid,tmp;
    register Long temp;


    /* Trivial Case */
    if( (x1==x2) && (y1==y2) )
    {   if( z1>z2 )
        {      DrawSphere(x1,y1,z1,rad,c1);
        } else DrawSphere(x2,y2,z2,rad,c2);
        return;
    }

    if( z1<z2 )
    {   tmp=x1; x1=x2; x2=tmp;
        tmp=y1; y1=y2; y2=tmp;
        tmp=z1; z1=z2; z2=tmp;
        tmp=c1; c1=c2; c2=tmp;
    }

    DrawCylinderCaps(x1,y1,z1,c1,x2,y2,z2,c2,rad);

    lx = x2-x1;
    ly = y2-y1;
    lz = z2-z1;

    if( ly>0 ) { ystep = XRange; ay = ly; iy = 1; }
    else {   ystep = -XRange; ay = -ly; iy = -1; }
    if( lx>0 ) { ax = lx; ix = 1; }
    else { ax = -lx; ix = -1; }
    zrate = lz/MaxFun(ax,ay);

    temp = (Long)y1*XRange+x1;
    fbase = FBuffer+temp;
    dbase = DBuffer+temp;

    if( ax>ay )
    {   lz -= ax*zrate;
        zerr = err = -(ax>>1);

        if( c1 != c2 )
        {   mid = (x1+x2)>>1;
            while( x1!=mid )
            {   z1 += zrate;  if( (zerr-=lz)>0 ) { zerr-=ax; z1--; }
                fbase+=ix; dbase+=ix; x1+=ix;
                if( (err+=ay)>0 )
                {   fbase+=ystep; dbase+=ystep; err-=ax;
                       DrawArcDn(dbase,fbase,z1,c1);
                } else DrawArcAc(dbase,fbase,z1,c1);
            }
        }

        while( x1!=x2 )
        {   z1 += zrate;  if( (zerr-=lz)>0 ) { zerr-=ax; z1--; }
            fbase+=ix; dbase+=ix; x1+=ix;
            if( (err+=ay)>0 )
            {   fbase+=ystep; dbase+=ystep; err-=ax;
                   DrawArcDn(dbase,fbase,z1,c2);
            } else DrawArcAc(dbase,fbase,z1,c2);
        }
    } else /*ay>=ax*/
    {   lz -= ay*zrate;
        zerr = err = -(ay>>1);

        if( c1 != c2 )
        {   mid = (y1+y2)>>1;
            while( y1!=mid )
            {   z1 += zrate;  if( (zerr-=lz)>0 ) { zerr-=ay; z1--; }
                fbase+=ystep; dbase+=ystep; y1+=iy;
                if( (err+=ax)>0 )
                {   fbase+=ix; dbase+=ix; err-=ay; 
                       DrawArcAc(dbase,fbase,z1,c1);
                } else DrawArcDn(dbase,fbase,z1,c1);
            }
        }

        while( y1!=y2 )
        {   z1 += zrate;  if( (zerr-=lz)>0 ) { zerr-=ay; z1--; }
            fbase+=ystep; dbase+=ystep; y1+=iy;
            if( (err+=ax)>0 )
            {   fbase+=ix; dbase+=ix; err-=ay; 
                   DrawArcAc(dbase,fbase,z1,c2);
            } else DrawArcDn(dbase,fbase,z1,c2);
        }
    }
}


static int TestCylinder( x1,y1,z1,x2,y2,z2,rad )
    int x1,y1,z1,x2,y2,z2,rad;
{
    register int tmp1, tmp2;

    if( UseSlabPlane )
        if( (z1+rad>SlabValue) || (z2+rad>SlabValue) )
            return(False);

    ClipStatus = False;

    tmp1 = x1+rad;  tmp2 = x2+rad;
    if( (tmp1<0) && (tmp2<0) )
        return( False );
    if( (tmp1>=XRange) || (tmp2>=XRange) )
        ClipStatus = True;

    tmp1 = x1-rad;  tmp2 = x2-rad;
    if( (tmp1>=XRange) && (tmp2>=XRange) )
        return( False );
    if( (tmp1<0) || (tmp2<0) )
        ClipStatus = True;

    tmp1 = y1+rad;  tmp2 = y2+rad;
    if( (tmp1<0) && (tmp2<0) )
        return( False );
    if( (tmp1>=YRange) || (tmp2>=YRange) )
        ClipStatus = True;

    tmp1 = y1-rad;  tmp2 = y2-rad;
    if( (tmp1>=YRange) && (tmp2>=YRange) )
        return( False );
    if( (tmp1<0) || (tmp2<0) )
        ClipStatus = True;

    return( True );
}

#if defined(__STDC__) || defined(IBMPC)
/* Function Prototypes */
static void ClipArcDn( short __huge*, Pixel __huge*, int, int, int, int );
static void ClipArcAc( short __huge*, Pixel __huge*, int, int, int, int );
#endif


static void ClipArcAc(dbase,fbase,x,y,z,c)
    register short __huge *dbase;
    register Pixel __huge *fbase;
    register int x,y,z,c;
{
    register ArcEntry __far *ptr;
    register short __huge *dptr;
    register short depth;
    register int temp;

    ptr = ArcAc;
    while( (temp=y+ptr->dy) < 0 )
        if( ++ptr == ArcAcPtr )
            return;

    while( (temp<YRange) && (ptr<ArcAcPtr) )
    {   temp = x+ptr->dx;
        if( XValid(temp) )
        {   dptr = dbase+ptr->offset;  depth = ptr->dz+z;
            SETPIXEL(dptr,fbase+ptr->offset,depth,Lut[ptr->inten+c]);
        }
        ptr++;
        temp = y+ptr->dy;
    }
}

static void ClipArcDn(dbase,fbase,x,y,z,c)
    register short __huge *dbase;
    register Pixel __huge *fbase;
    register int x,y,z,c;
{
    register ArcEntry __far *ptr;
    register short __huge *dptr;
    register short depth;
    register int temp;

    ptr = ArcDn;
    while( (temp=y+ptr->dy) < 0 )
        if( ++ptr == ArcDnPtr )
            return;

    while( (temp<YRange) && (ptr<ArcDnPtr) )
    {   temp = x+ptr->dx;
        if( XValid(temp) )
        {   dptr = dbase+ptr->offset;  depth = ptr->dz+z;
            SETPIXEL(dptr,fbase+ptr->offset,depth,Lut[ptr->inten+c]);
        }
        ptr++;
        temp = y+ptr->dy;
    }
}


void ClipCylinder( x1,y1,z1,c1,x2,y2,z2,c2,rad )
    int x1,y1,z1,c1, x2,y2,z2,c2, rad;
{
    register short __huge *dbase;
    register Pixel __huge *fbase;

    register int zrate,zerr,ystep,err;
    register int ix,iy,ax,ay;
    register int lx,ly,lz;
    register int mid,tmp;
    register Long temp;

    /* Visibility Tests */
    if( !TestCylinder(x1,y1,z1,x2,y2,z2,rad) )
        return;

    if( !ClipStatus )
    {   DrawCylinder(x1,y1,z1,c1,x2,y2,z2,c2,rad);
        return;
    }

    /* Trivial Case */
    if( (x1==x2) && (y1==y2) )
    {   if( z1>z2 )
        {      ClipSphere(x1,y1,z1,rad,c1);
        } else ClipSphere(x2,y2,z2,rad,c2);
        return;
    }

    if( z1<z2 )
    {   tmp=x1; x1=x2; x2=tmp;
        tmp=y1; y1=y2; y2=tmp;
        tmp=z1; z1=z2; z2=tmp;
        tmp=c1; c1=c2; c2=tmp;
    }

    DrawCylinderCaps(x1,y1,z1,c1,x2,y2,z2,c2,rad);

    lx = x2-x1;
    ly = y2-y1;
    lz = z2-z1;

    if( ly>0 ) { ystep = XRange; ay = ly; iy = 1; }
    else {   ystep = -XRange; ay = -ly; iy = -1; }
    if( lx>0 ) { ax = lx; ix = 1; }
    else { ax = -lx; ix = -1; }
    zrate = lz/MaxFun(ax,ay);

    temp = (Long)y1*XRange+x1;
    fbase = FBuffer+temp;
    dbase = DBuffer+temp;

    if( ax>ay )
    {   if( x2<x1 )
        {   tmp = c1;
            c1 = c2;
            c2 = tmp;
        }
        lz -= ax*zrate;
        zerr = err = -(ax>>1);
        mid = (x1+x2)/2;

        while( x1!=x2 )
        {   z1 += zrate;  if( (zerr-=lz)>0 ) { zerr-=ax; z1--; }
            fbase+=ix; dbase+=ix; x1+=ix;
            if( (err+=ay)>0 )
            {   fbase += ystep;  err -= ax;
                dbase += ystep;  y1 += iy;
                   ClipArcDn(dbase,fbase,x1,y1,z1,(x1<mid?c1:c2));
            } else ClipArcAc(dbase,fbase,x1,y1,z1,(x1<mid?c1:c2));
        }
    } else /*ay>=ax*/
    {   if( y2<y1 )
        {   tmp = c1;
            c1 = c2;
            c2 = tmp;
        }
        lz -= ay*zrate;
        zerr = err = -(ay>>1);
        mid = (y1+y2)/2;

        while( y1!=y2 )
        {   z1 += zrate;  if( (zerr-=lz)>0 ) { zerr-=ay; z1--; }
            fbase+=ystep; dbase+=ystep; y1+=iy;
            if( (err+=ax)>0 )
            {   fbase += ix;  err -= ay;
                dbase += ix;  x1 += ix; 
                   ClipArcAc(dbase,fbase,x1,y1,z1,(y1<mid?c1:c2));
            } else ClipArcDn(dbase,fbase,x1,y1,z1,(y1<mid?c1:c2));
        }
    }
}


void InitialisePixUtils()
{
#ifdef IBMPC
    ArcAc = (ArcEntry __far*)_fmalloc(ARCSIZE*sizeof(ArcEntry));
    ArcDn = (ArcEntry __far*)_fmalloc(ARCSIZE*sizeof(ArcEntry));
#endif
    SplineCount = 5;
}


