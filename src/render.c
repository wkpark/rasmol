/* render.c
 * RasMol2 Molecular Graphics
 * Roger Sayle, February 1994
 * Version 2.3
 */
#include "rasmol.h"

#ifdef IBMPC
#include <windows.h>
#include <malloc.h>
#endif
#ifndef sun386
#include <stdlib.h>
#endif

#include <stdio.h>
#include <math.h>

#define RENDER
#include "graphics.h"
#include "render.h"
#include "molecule.h"
#include "transfor.h"
#include "command.h"
#include "abstree.h"
#include "pixutils.h"


#define PoolSize        16
#define RootSix         2.4494897427831780
#define OneOverRootSix  0.4082482904638630
#define TwoOverRootSix  0.8164965809277260
#define ApproxZero      1.0E-3
#define INFINITY       200000
#define FUDGEFACTOR    1000


typedef struct _Item {
                struct _Item __far *list;
                Atom  __far *data;
               } Item;

typedef struct { Real h,l; } Interval;

static Atom __far * __far *YBucket;
static Atom __far * __far *IBuffer;
static int BuckY,ItemX;
static int FBufX,FBufY;
static int DBClear;

static Atom __far *SBuffer;
static Atom __far *Exclude;
static Real ShadowI, ShadowJ, ShadowK;
static int  ShadowX, ShadowY, ShadowZ;
static int deltax, deltay, deltaz;
static int xcord, ycord, zcord;
static int xflag, yflag, zflag;
static int xhash, yhash, zhash;
static int RayCount;

static Item __far *FreeItem;
static Real VoxRatio,IVoxRatio;
static int VoxelCount,InVoxCount;
static int ProbeCount;
static int VoxelsDone;


/* Macros for commonly used loops */
#define ForEachAtom  for(chain=Database->clist;chain;chain=chain->cnext) \
                     for(group=chain->glist;group;group=group->gnext)    \
                     for(aptr=group->alist;aptr;aptr=aptr->anext)
#define ForEachBond  for(bptr=Database->blist;bptr;bptr=bptr->bnext)
#define ForEachBack  for(chain=Database->clist;chain;chain=chain->cnext) \
                     for(bptr=chain->blist;bptr;bptr=bptr->bnext)




static void FatalRenderError(ptr)
    char *ptr;
{
    char buffer[80];

    sprintf(buffer,"Renderer Error: Unable to allocate %s!",ptr);
    RasMolFatalExit(buffer);
}


int isqrt( val )
    Card val;
{
    register Card side, left;
    register Card temp, result;
    register int i;

    result = side = left = 0;
    for( i=0; i<sizeof(Card)*4; i++ )
    {   left = (left<<2) + (val>>30);
        result <<= 1; side <<= 1;
        temp = side | 1;
        val <<= 2;

        if( left >= temp )
        {   side = temp+1;
            left -= temp;
            result |= 1;
        }
    }
    return( (int)result );
}


void ClearBuffers()
{
    register Long __huge *ptr;
    register Long __huge *end;
    register Long fill;

    if( !FBClear )
    {   FBClear = True;

#ifdef IBMPC
        FBuffer = (Pixel __huge*)GlobalLock(FBufHandle);
        fill = 0;
#else
        fill = Lut[5];
#ifdef EIGHTBIT
        fill |= fill<<8;
        fill |= fill<<16;
#endif
#endif

        ptr = (Long __huge*)FBuffer;
        end = (Long __huge*)(FBuffer+(Long)XRange*YRange);
        do { *ptr++=fill; *ptr++=fill;
             *ptr++=fill; *ptr++=fill;
        } while( ptr<end );
#ifdef IBMPC
        GlobalUnlock(FBufHandle);
#endif
    }

    if( !DBClear )
    {   DBClear = True;
#ifdef IBMPC
        DBuffer = (short __huge*)GlobalLock(DBufHandle);
#endif
        ptr = (Long __huge*)DBuffer;
        end = (Long __huge*)(DBuffer+(Long)XRange*YRange);
        do { *ptr++=0; *ptr++=0;
             *ptr++=0; *ptr++=0;
        } while( ptr<end );
#ifdef IBMPC
        GlobalUnlock(DBufHandle);
#endif
    }
}


void ReAllocBuffers()
{
    register Atom __far * __far *iptr;
    register int index,len;
    register Long temp;

    temp = (Long)XRange*YRange*sizeof(short)+32;
#ifdef IBMPC
    if( DBufHandle ) GlobalFree(DBufHandle);
    DBufHandle = GlobalAlloc(GMEM_MOVEABLE,temp);
    if( !DBufHandle ) FatalRenderError("depth buffer");
#else
    if( DBuffer ) free( DBuffer );
    DBuffer = (short*)malloc( temp );
    if( !DBuffer ) FatalRenderError("depth buffer");
#endif
    DBClear=False;

    if( YBucket && (BuckY<YRange) )
    {   _ffree(YBucket); 
        YBucket=(void __far*)0; 
    }

    if( !YBucket )
    {   len = YRange*sizeof(Atom __far*);
        YBucket = (Atom __far* __far*)_fmalloc( len );
        if( !YBucket ) FatalRenderError("Y buckets");
        BuckY = YRange;
    }

    if( IBuffer && (ItemX<XRange) )
    {   _ffree(IBuffer); 
        IBuffer=(void __far*)0; 
    }

    if( !IBuffer )
    {   len = (XRange+4)*sizeof(Atom __far*);
        IBuffer = (Atom __far* __far*)_fmalloc(len);
        if( !IBuffer ) FatalRenderError("item buffer");
        len = XRange>>2;  iptr = IBuffer;
        for( index=0; index<=len; index++ )
        {   *iptr++ = (void __far*)0;  *iptr++ = (void __far*)0;
            *iptr++ = (void __far*)0;  *iptr++ = (void __far*)0;
        }
        ItemX = XRange;
    }
}


void ReSizeScreen()
{
    register Real orig;
    register Long temp;

    if( Range != ZoomRange )
    {   orig = MaxZoom;
        MaxZoom = 0.236*(WorldSize+1000)/Range;
        ZoomRange = Range;  MaxZoom -= 1.0;

        /* Handle Change in MaxZoom */
        if( DialValue[3]>0.0 )
        {   DialValue[3] *= orig/MaxZoom;
            if( DialValue[3]>1.0 )
                DialValue[3] = 1.0;
        }
    }

#ifdef IBMPC
    if( !FBufHandle || (FBufX!=XRange) || (FBufY!=YRange) )
    {   temp = (Long)XRange*YRange*sizeof(Pixel)+16;
        if( FBufHandle ) GlobalFree(FBufHandle);
        FBufHandle = GlobalAlloc(GMEM_MOVEABLE,temp);
        if( !FBufHandle ) FatalRenderError("frame buffer");

        BucketFlag = False;
        FBufX=XRange;  FBufY=YRange;  FBClear = False;
        ReAllocBuffers();
        ClearBuffers();
    }
#else /* UNIX */
    if( !FBuffer || (FBufX!=XRange) || (FBufY!=YRange) )
    {   if( !Interactive )
        {   if( FBuffer ) free(FBuffer);
            temp = (Long)XRange*YRange*sizeof(Pixel);
            FBuffer = (Pixel*)malloc( temp+32 );
        } else FBuffer = CreateImage();
        if( !FBuffer ) FatalRenderError("frame buffer");

        BucketFlag = False;
        FBufX=XRange;  FBufY=YRange;  FBClear = False;
        ReAllocBuffers();
        ClearBuffers();
    }
#endif
}


static void DisplayWireFrame()
{
    register Bond __far *bptr;
    register Atom __far *s;
    register Atom __far *d;
    register int sc,dc;

    if( UseClipping )
    {   ForEachBond
           if( bptr->flag&WireFlag )
           {   s = bptr->srcatom; d = bptr->dstatom;
               if( !bptr->col ) 
               {   sc = s->col;  dc = d->col;
               } else sc = dc = bptr->col;
               ClipTwinVector(s->x,s->y,s->z,d->x,d->y,d->z,sc,dc);
           } else if( bptr->flag&CylinderFlag )
           {   s = bptr->srcatom; d = bptr->dstatom;
               if( !bptr->col ) 
               {   sc = s->col;  dc = d->col;
               } else sc = dc = bptr->col;

               if( bptr->irad>0 )
               {  ClipCylinder(s->x,s->y,s->z,sc,d->x,d->y,d->z,dc,bptr->irad);
               } else ClipTwinLine(s->x,s->y,s->z,d->x,d->y,d->z,sc,dc);
           }
    } else
        ForEachBond
           if( bptr->flag&WireFlag )
           {   s = bptr->srcatom; d = bptr->dstatom;
               if( !bptr->col ) 
               {   sc = s->col;  dc = d->col;
               } else sc = dc = bptr->col;
               DrawTwinVector(s->x,s->y,s->z,d->x,d->y,d->z,sc,dc);
           } else if( bptr->flag&CylinderFlag )
           {   s = bptr->srcatom; d = bptr->dstatom;
               if( !bptr->col ) 
               {   sc = s->col;  dc = d->col;
               } else sc = dc = bptr->col;

               if( bptr->irad>0 )
               {  DrawCylinder(s->x,s->y,s->z,sc,d->x,d->y,d->z,dc,bptr->irad);
               } else DrawTwinLine(s->x,s->y,s->z,d->x,d->y,d->z,sc,dc);
           }
}


static void DisplayBackBone()
{
    register Chain __far *chain;
    register Bond __far *bptr;
    register Atom __far *s;
    register Atom __far *d;
    register int sc,dc;

    ForEachBack
       if( bptr->flag&DrawBondFlag )
       {   s = bptr->srcatom; d = bptr->dstatom;
           if( !bptr->col ) 
           {   sc = s->col;  dc = d->col;
           } else sc = dc = bptr->col;

           if( bptr->flag&CylinderFlag )
           {   if( bptr->irad>0 )
               { ClipCylinder(s->x,s->y,s->z,sc,d->x,d->y,d->z,dc,bptr->irad);
               } else ClipTwinLine(s->x,s->y,s->z,d->x,d->y,d->z,sc,dc);
           } else ClipTwinVector(s->x,s->y,s->z,d->x,d->y,d->z,sc,dc);
       }
}


static void PrepareYBucket()
{
    register Atom __far * __far *temp;
    register Chain __far *chain;
    register Group __far *group;
    register Atom __far *aptr;
    register int scan;
    register int rad;

    temp = YBucket;
    for( scan=0; scan<BuckY; scan++ )
        *temp++ = (void __far*)0;

    if( UseClipping )
    {   ForEachAtom
            if( aptr->flag&SphereFlag )
            {   rad = aptr->irad;
                if( (aptr->x-rad>=XRange) || 
                    (aptr->x+rad<0) || (aptr->y+rad<0) )
                    continue;
                if( (scan=aptr->y-rad) > BuckY ) continue;

                if( scan>0 )
                {   aptr->bucket = YBucket[scan];
                    YBucket[scan] = aptr;
                } else
                {   aptr->bucket = *YBucket;
                    *YBucket = aptr;
                }
            }
    } else
        ForEachAtom
            if( aptr->flag&SphereFlag )
            {   scan = aptr->y-aptr->irad;
                aptr->bucket = YBucket[scan];
                YBucket[scan] = aptr;
            }
    BucketFlag = True;
}

#if defined(__STDC__) || defined(IBMPC)
/* Function Prototypes */
static void SqrInterval( Interval __far* );
static void VoxelInsert( Atom __far*, int );
static int AtomInter( Atom __far* );
#endif


static void SqrInterval( ival )
    Interval __far *ival;
{   register Real l,h;

    l = ival->l;
    h = ival->h;

    if( l>=0.0 )
    {   ival->l = l*l;
        ival->h = h*h;
    } else if( h<0.0 )
    {   ival->l = h*h;
        ival->h = l*l;
    } else
    {   ival->h = (-l>h)? l*l : h*h;
        ival->l = 0.0;
    }
}

static void VoxelInsert( ptr, ref )
    Atom __far *ptr;
    int ref;
{
    register int i;
    register Item __far *datum;

    if( !FreeItem )
    {   datum = (Item __far*)_fmalloc( PoolSize*sizeof(Item) );
        if( !datum ) FatalRenderError("voxel item");
        for( i=1; i<PoolSize; i++ )
        {   datum->list = FreeItem;
            FreeItem = datum++;
        }
    } else
    {   datum = FreeItem;
        FreeItem = datum->list;
    }
    datum->data = ptr;
    InVoxCount++;

    if( !HashTable[ref] ) VoxelCount++;
    datum->list = (Item __far*)HashTable[ref];
    HashTable[ref] = (void __far*)datum;
}


void ResetVoxelData()
{
    register Item __far *datum;
    register int i;

    if( VoxelsDone )
    {   for( i=0; i<VOXSIZE; i++ )
            if( datum = (Item __far*)HashTable[i] )
            {   while( datum->list ) datum = datum->list;
                datum->list = FreeItem;
                FreeItem = (Item __far*)HashTable[i];
                HashTable[i] = (void __far*)0;
            }
        VoxelsDone = False;
    } else for( i=0; i<VOXSIZE; i++ )
        HashTable[i] = (void __far*)0;
    VoxelsClean = True;
}


void CreateVoxelData()
{
    static Interval ix, iy, iz;
    register int lvx, lvy, lvz;
    register int hvx, hvy, hvz;
    register Long mx, my, mz;
    register int px, py, pz;
    register int i, rad;
    register Real rad2;

    register Chain __far *chain;
    register Group __far *group;
    register Atom __far *aptr;


    ResetVoxelData();
    ProbeCount = InVoxCount = VoxelCount = 0;
    VoxRatio = (Real)SideLen/VOXORDER;
    IVoxRatio = 1.0/VoxRatio;
    VoxelsDone = True;

    ForEachAtom
    if( aptr->flag&SphereFlag )
    {   mx = aptr->xorg+Offset;
        my = aptr->yorg+Offset;
        mz = aptr->zorg+Offset;
        rad = aptr->radius;  
        rad2 = (Long)rad*rad;

        lvx = (int)((mx-rad)*IVoxRatio);  hvx = (int)((mx+rad)*IVoxRatio);
        lvy = (int)((my-rad)*IVoxRatio);  hvy = (int)((my+rad)*IVoxRatio);
        lvz = (int)((mz-rad)*IVoxRatio);  hvz = (int)((mz+rad)*IVoxRatio);


        for( px=lvx; px<=hvx; px++ )
        {   ix.l=px*VoxRatio-mx;
            ix.h=ix.l+VoxRatio;  
            SqrInterval(&ix);
            i = VOXORDER2*px + VOXORDER*lvy;
       
            for( py=lvy; py<=hvy; py++ )
            {   iy.l=py*VoxRatio-my;
                iy.h=iy.l+VoxRatio;
                SqrInterval(&iy);
                
                for( pz=lvz; pz<=hvz; pz++ )
                {   iz.l=pz*VoxRatio-mz; 
                    iz.h=iz.l+VoxRatio;
                    SqrInterval(&iz);

                    if( ((ix.h+iy.h+iz.h)>rad2) && ((ix.l+iy.l+iz.l)<rad2) )
                        VoxelInsert( aptr, i+pz );
                } /*pz*/
                i += VOXORDER;
	    } /*py*/
        } /*px*/
    }
}


void ShadowTransform()
{
    ShadowI = OneOverRootSix*(DirX[0]-DirX[1]+DirX[2]+DirX[2]);
    ShadowK = OneOverRootSix*(DirZ[0]-DirZ[1]+DirZ[2]+DirZ[2]);
#ifdef INVERT
    ShadowJ = -OneOverRootSix*(DirY[0]-DirY[1]+DirY[2]+DirY[2]);
#else
    ShadowJ = OneOverRootSix*(DirY[0]-DirY[1]+DirY[2]+DirY[2]);
#endif

    if( ShadowI>ApproxZero )
    {   deltax =  (int)(FUDGEFACTOR/ShadowI); xhash =  VOXORDER2; xflag =  1;
    } else if( ShadowI<-ApproxZero )
    {   deltax = -(int)(FUDGEFACTOR/ShadowI); xhash = -VOXORDER2; xflag = -1;
    } else xflag = 0;

    if( ShadowJ>ApproxZero )
    {   deltay =  (int)(FUDGEFACTOR/ShadowJ); yhash =  VOXORDER; yflag =  1;
    } else if( ShadowJ<-ApproxZero )
    {   deltay = -(int)(FUDGEFACTOR/ShadowJ); yhash = -VOXORDER; yflag = -1;
    } else yflag = 0;

    if( ShadowK>ApproxZero )
    {   deltaz =  (int)(FUDGEFACTOR/ShadowK); zhash = zflag =  1;
    } else if( ShadowK<-ApproxZero )
    {   deltaz = -(int)(FUDGEFACTOR/ShadowK); zhash = zflag = -1;
    } else zflag = 0;
}


static int AtomInter( ptr )
    Atom __far *ptr;
{
    register int vx, vy, vz;
    register Long modv,rad2;
    register Real tca;

    if( ptr->mbox == RayCount )
        return( False );
    ptr->mbox = RayCount;

    vx = (int)ptr->xorg-ShadowX;
    vy = (int)ptr->yorg-ShadowY;
    vz = (int)ptr->zorg-ShadowZ;

    tca = vx*ShadowI + vy*ShadowJ + vz*ShadowK;
    if( tca<0.0 ) return( False );
    
    rad2 = ptr->radius+10;  rad2 = rad2*rad2;
    modv = (Long)vx*vx + (Long)vy*vy + (Long)vz*vz - rad2;
    return( modv<tca*tca );
}


static int ShadowRay()
{
    register Item __far * __far *ident;
    register Item __far *ptr;
    register Real ex, ey, ez;
    register Long dx, dy, dz;
    register int ref;

   
    RayCount++;
    if( SBuffer )
    {   if( (SBuffer!=Exclude) && AtomInter(SBuffer) )
            return( True );
        SBuffer = (void __far*)0;
    }

    ex = IVoxRatio*(ShadowX+Offset);  xcord = (int)ex;
    ey = IVoxRatio*(ShadowY+Offset);  ycord = (int)ey;
    ez = IVoxRatio*(ShadowZ+Offset);  zcord = (int)ez;

    ref = VOXORDER2*xcord+VOXORDER*ycord+zcord;
    ident = (Item __far* __far*)(HashTable+ref);

    if( xflag==1 ) 
    {   dx = (Long)(((xcord+1)-ex)*deltax);
    } else if( xflag == -1 )
    {   dx = (Long)((ex-xcord)*deltax); 
    } else dx = INFINITY;

    if( yflag==1 ) 
    {   dy = (Long)(((ycord+1)-ey)*deltay);
    } else if( yflag == -1 )
    {   dy = (Long)((ey-ycord)*deltay); 
    } else dy = INFINITY;

    if( zflag==1 ) 
    {   dz = (Long)(((zcord+1)-ez)*deltaz);
    } else if( zflag == -1 )
    {   dz = (Long)((ez-zcord)*deltaz); 
    } else dz = INFINITY;

    
    while( True )
    {   for( ptr = *ident; ptr; ptr = ptr->list )
            if( (ptr->data!=Exclude) && AtomInter(ptr->data) )
            {   SBuffer = ptr->data;
                return( True );
            }

        if( (dx<=dy) && (dx<=dz) )
        {   xcord += xflag;
            if( (xcord<0) || (xcord>=VOXORDER) ) return( False );
            ident += xhash; dx += deltax;
        } else if( dy<=dz  ) /*(dy<=dx)*/
        {   ycord += yflag;
            if( (ycord<0) || (ycord>=VOXORDER) ) return( False );
            ident += yhash; dy += deltay;
        } else /* (dz<=dx) && (dz<=dy) */
        {   zcord += zflag;
            if( (zcord<0) || (zcord>=VOXORDER) ) return( False );
            ident += zhash; dz += deltaz;
        }
    }
}



#define UpdateScanAcross \
        if( depth>*dptr )   \
        {   *dptr = depth;  \
            iptr[dx] = ptr; \
        } dptr++; dx++;


/* ScanLine for Shadows! */
static void ScanLine()
{
    static Atom __far *list;
    register Atom __far *ptr;
    register Atom __far * __far *iptr;
    register Atom __far * __far *prev;
    register short __huge *dbase;
    register short __huge *dptr;
    register Pixel __huge *fptr;
    register char __far *tptr;

    register int pos,depth,inten;
    register int lastx,wide,scan;
    register int dx,dy,dz;

    fptr = FBuffer;
    dbase = DBuffer;
    list = (void __far*)0;  

    wide = XRange>>2;  iptr = IBuffer;
    for( pos=0; pos<=wide; pos++ )
    {   *iptr++ = (void __far*)0;  *iptr++ = (void __far*)0;
        *iptr++ = (void __far*)0;  *iptr++ = (void __far*)0;
    }


    for( scan=0; scan<YRange; scan++ )
    {   for( ptr = YBucket[scan]; ptr; ptr = ptr->bucket )
        {    ptr->next = list; list = ptr; }

        prev = &list;
        for( ptr=list; ptr; ptr=ptr->next )
        {   dy = scan - ptr->y;
            wide = LookUp[ptr->irad][AbsFun(dy)];
            lastx = (XRange-1)-ptr->x;
            if( wide<lastx ) lastx=wide;
            dx = - MinFun(wide,ptr->x);

            iptr = IBuffer+ptr->x;
            tptr = LookUp[wide];

            dptr = dbase+ptr->x+dx;
            while( dx<=lastx )
            {   depth = tptr[AbsFun(dx)]+ptr->z;
                UpdateScanAcross;
            }

            /* Remove completed atoms */
            if( dy == ptr->irad )
            {   *prev = ptr->next;
            } else prev = &ptr->next;
        } /*ptr*/


        /* Process visible scanline */
        prev = (Atom __far* __far*)IBuffer;
        SBuffer = (void __far*)0;
        dptr = dbase; 

        for( pos=0; pos<XRange; pos++ )
        {   if( ptr = *prev )
            {   dz = *dptr-ptr->z;
#ifdef INVERT
                inten = (dz<<1)+(pos-ptr->x)+(scan-ptr->y);
#else
                inten = (dz<<1)+(pos-ptr->x)-(scan-ptr->y);
#endif
                if( inten>0 )
                {   inten = (int)( (inten*ColConst[ptr->irad])>>ColBits);
                    dz = *dptr-ImageRadius;
                    dy = scan-YOffset;
                    dx = pos-XOffset;

                    ShadowX = (int)(dx*InvX[0] + dy*InvX[1] + dz*InvX[2]);
                    ShadowY = (int)(dx*InvY[0] + dy*InvY[1] + dz*InvY[2]);
                    ShadowZ = (int)(dx*InvZ[0] + dy*InvZ[1] + dz*InvZ[2]);

                    Exclude = ptr;
                    if( ShadowRay() )
                    {   *fptr = Lut[ptr->col+(inten>>2)];
                    } else *fptr = Lut[ptr->col+inten];
                } else *fptr = Lut[ptr->col];
                *prev = (void __far*)0;
            }
            dptr++; fptr++; prev++;
        }
        dbase = dptr;
    } /*scan*/
}


static void DisplaySpaceFill()
{
    register Chain __far *chain;
    register Group __far *group;
    register Atom __far *aptr;

    if( UseShadow )
    {   if( !BucketFlag )
            PrepareYBucket();
        ScanLine();
    } else if( UseClipping )
    {   ForEachAtom
            if( aptr->flag&SphereFlag )
                ClipSphere(aptr->x,aptr->y,aptr->z,aptr->irad,aptr->col);
    } else 
        ForEachAtom
            if( aptr->flag&SphereFlag )
                DrawSphere(aptr->x,aptr->y,aptr->z,aptr->irad,aptr->col);
}


#if defined(__STDC__) || defined(IBMPC)
/* Function Prototype */
static void DisplayHBonds( HBond __far *, int );
#endif

static void DisplayHBonds( list, mode )
    HBond __far *list; 
    int mode;
{
    register HBond __far *ptr;
    register Atom __far *s;
    register Atom __far *d;
    register int sc,dc;

    for( ptr=list; ptr; ptr=ptr->hnext )
        if( ptr->flag & DrawBondFlag )
        {   if( mode )
            {   s = ptr->srcCA; d = ptr->dstCA;
                if( !s || !d ) continue;
            } else
            {   d = ptr->src;
                s = ptr->dst;
            }

            if( !ptr->col )
            {   sc = s->col;  dc = d->col;
            } else sc = dc = ptr->col;
            if( ptr->flag & CylinderFlag )
            {   if( ptr->irad>0 )
                {   ClipCylinder(s->x,s->y,s->z,sc,
                                 d->x,d->y,d->z,dc,ptr->irad);
                } else ClipTwinLine(s->x,s->y,s->z,d->x,d->y,d->z,sc,dc);
            } else ClipDashVector(s->x,s->y,s->z,d->x,d->y,d->z,sc,dc);
        }
}

#if defined(__STDC__) || defined(IBMPC)
/* Function Prototype */
static void DisplayRibbon( Chain __far * );
#endif


static void CalculateInten( ptr )
    Knot *ptr;
{
    register Real inten;
    register int size;

    size = isqrt( (Long)ptr->cx*ptr->cx + 
                  (Long)ptr->cy*ptr->cy + 
                  (Long)ptr->cz*ptr->cz );

    if( size )
    {   inten = (Real)(ptr->cx-ptr->cy+ptr->cz+ptr->cz)/(size*RootSix);
        if( ptr->cz < 0 ) inten = -inten;

        ptr->inten = (int)(ColourMask*inten);
        if( ptr->inten<0 ) ptr->inten = 0;
    } else ptr->inten = ColourMask;
}


static void DisplayRibbon( chain )
    Chain  __far *chain;
{
    register Group __far *group;
    register Atom __far *captr;
    register Atom __far *o1ptr;
    register Atom __far *o2ptr;
    register Atom __far *next;

    register int prev,wide;
    register int col1,col2;
    register int bx,by,bz;
    register int dx,dy,dz;
    register int size;

    static Knot mid1, mid2, mid3;
    static Knot knot1, knot2;

    prev = False;  
    group = chain->glist;
    if( IsAmino(group->refno) )
    {   captr = FindGroupAtom(group,1);
    } else captr = FindGroupAtom(group,7);

    while( group->gnext )
    {   if( IsAmino(group->gnext->refno) )
        {   next = FindGroupAtom(group->gnext,1);
            o1ptr = FindGroupAtom(group,3);
        } else /* Nucleic Acid */
        {   next = FindGroupAtom(group->gnext,7);
            o1ptr = FindGroupAtom(group->gnext,10);
        }

        /* When not to have a control point! */
        if( !next || !captr || !o1ptr || (next->flag&BreakFlag) ||
            !((group->flag|group->gnext->flag)&RibbonFlag) )
        {   group = group->gnext;
            captr = next;
            prev = False;
            continue;
        }

        knot2.tx = next->x - captr->x;
        knot2.ty = next->y - captr->y;
        knot2.tz = next->z - captr->z;

        if( IsAmino(group->refno) )
        {   bx = o1ptr->x - captr->x;
            by = o1ptr->y - captr->y;
            bz = o1ptr->z - captr->z;

        } else if( !FindGroupAtom(group,17) && 
                   (o2ptr=FindGroupAtom(group,8)) )
        {   /* Deoxyribonucleic Acid */
            o2ptr = FindGroupAtom(group,8);
            bx = (o1ptr->x + o2ptr->x)/2 - captr->x;
            by = (o1ptr->y + o2ptr->y)/2 - captr->y;
            bz = (o1ptr->z + o2ptr->z)/2 - captr->z;

        } else /* Ribonucleic Acid */
        {   bx = o1ptr->x - captr->x;
            by = o1ptr->y - captr->y;
            bz = o1ptr->z - captr->z;
        }

        /* c := a x b */
        knot2.cx = knot2.ty*bz - knot2.tz*by;
        knot2.cy = knot2.tz*bx - knot2.tx*bz;
        knot2.cz = knot2.tx*by - knot2.ty*bx;

        knot2.px = (captr->x + next->x)/2;
        knot2.py = (captr->y + next->y)/2;
        knot2.pz = (captr->z + next->z)/2;

        if( (group->flag&group->gnext->flag) & HelixFlag )
        {   size = isqrt((Long)knot2.cx*knot2.cx + 
                         (Long)knot2.cy*knot2.cy + 
                         (Long)knot2.cz*knot2.cz);

            if( size )
            {   wide = (int)(375*Scale);

#ifdef INVERT
                knot2.px += (int)(((Long)wide*knot2.cx)/size);
                knot2.py += (int)(((Long)wide*knot2.cy)/size);
                knot2.pz += (int)(((Long)wide*knot2.cz)/size);
#else
                knot2.px -= (int)(((Long)wide*knot2.cx)/size);
                knot2.py -= (int)(((Long)wide*knot2.cy)/size);
                knot2.pz -= (int)(((Long)wide*knot2.cz)/size);
#endif
            }
        }

        /* d := c x a */
        dx = (int)(((Long)knot2.cy*knot2.tz - 
                    (Long)knot2.cz*knot2.ty)/96);
        dy = (int)(((Long)knot2.cz*knot2.tx - 
                    (Long)knot2.cx*knot2.tz)/96);
        dz = (int)(((Long)knot2.cx*knot2.ty - 
                    (Long)knot2.cy*knot2.tx)/96);

        /* Average Ribbon Width */
        if( group->flag & RibbonFlag )
        {   if( group->gnext->flag & RibbonFlag )
            {   wide = (group->width+group->gnext->width)>>1;
            } else wide = group->width;
        } else wide = group->gnext->width;
        wide = (int)(wide*Scale);

        size = isqrt((Long)dx*dx + (Long)dy*dy + (Long)dz*dz);

        if( size )
        {   dx = (int)(((Long)wide*dx)/size);
            dy = (int)(((Long)wide*dy)/size);
            dz = (int)(((Long)wide*dz)/size);

            /* Handle Carbonyl Oxygen Flip */
            if( prev && ((knot1.wx*dx + knot1.wy*dy + knot1.wz*dz) < 0) )
            {   knot2.wx = -dx;
                knot2.wy = -dy;
                knot2.wz = -dz;
            } else
            {   knot2.wx = dx;
                knot2.wy = dy;
                knot2.wz = dz;
            }
        } else
        {   knot2.wx = 0;
            knot2.wy = 0;
            knot2.wz = 0;
        }


        if( RibbonMode )
            CalculateInten( &knot2 );
        if( !(col2 = group->col) )
            col2 = captr->col;
            

        if( prev )
        {   /* Approximate spline segment with straight line! */
            /* StrandRibbon( &knot1, &knot2, col1, col2 );   */

            /* Calculate Hermite Spline Points */
            mid1.px = (int)(((Long)54*knot1.px + (Long)9*knot1.tx +
                             (Long)10*knot2.px - (Long)3*knot2.tx)>>6);
            mid1.py = (int)(((Long)54*knot1.py + (Long)9*knot1.ty +
                             (Long)10*knot2.py - (Long)3*knot2.ty)>>6);
            mid1.pz = (int)(((Long)54*knot1.pz + (Long)9*knot1.tz +
                             (Long)10*knot2.pz - (Long)3*knot2.tz)>>6);

            mid2.px = (4*knot1.px + knot1.tx + 4*knot2.px - knot2.tx)>>3;
            mid2.py = (4*knot1.py + knot1.ty + 4*knot2.py - knot2.ty)>>3;
            mid2.pz = (4*knot1.pz + knot1.tz + 4*knot2.pz - knot2.tz)>>3;

            mid3.px = (int)(((Long)10*knot1.px + (Long)3*knot1.tx +
                             (Long)54*knot2.px - (Long)9*knot2.tx)>>6);
            mid3.py = (int)(((Long)10*knot1.py + (Long)3*knot1.ty +
                             (Long)54*knot2.py - (Long)9*knot2.ty)>>6);
            mid3.pz = (int)(((Long)10*knot1.pz + (Long)3*knot1.tz +
                             (Long)54*knot2.pz - (Long)9*knot2.tz)>>6);

            /* Calculate Hermite Spline Widths */
            mid1.wx = (27*knot1.wx + 5*knot2.wx)>>5;
            mid1.wy = (27*knot1.wy + 5*knot2.wy)>>5;
            mid1.wz = (27*knot1.wz + 5*knot2.wz)>>5;

            mid2.wx = (knot1.wx + knot2.wx)>>1;
            mid2.wy = (knot1.wy + knot2.wy)>>1;
            mid2.wz = (knot1.wz + knot2.wz)>>1;

            mid3.wx = (5*knot1.wx + 27*knot2.wx)>>5;
            mid3.wy = (5*knot1.wy + 27*knot2.wy)>>5;
            mid3.wz = (5*knot1.wz + 27*knot2.wz)>>5;

            /* Draw the Spline Segments */
            if( RibbonMode ) /* Solid! */
            {   mid1.cx = (27*knot1.cx + 5*knot2.cx)>>5;
                mid1.cy = (27*knot1.cy + 5*knot2.cy)>>5;
                mid1.cz = (27*knot1.cz + 5*knot2.cz)>>5;
                CalculateInten(&mid1);

                mid2.cx = (knot1.cx + knot2.cx)>>1;
                mid2.cy = (knot1.cy + knot2.cy)>>1;
                mid2.cz = (knot1.cz + knot2.cz)>>1;
                CalculateInten(&mid2);

                mid3.cx = (5*knot1.cx + 27*knot2.cx)>>5;
                mid3.cy = (5*knot1.cy + 27*knot2.cy)>>5;
                mid3.cz = (5*knot1.cz + 27*knot2.cz)>>5;
                CalculateInten(&mid3);

                SolidRibbon( &knot1, &mid1,  col1 );
                SolidRibbon( &mid1,  &mid2,  col1 );
                SolidRibbon( &mid2,  &mid3,  col2 );
                SolidRibbon( &mid3,  &knot2, col2 );
                
            } else /* Strands! */
            {   StrandRibbon( &knot1, &mid1,  col1 );
                StrandRibbon( &mid1,  &mid2,  col1 );
                StrandRibbon( &mid2,  &mid3,  col2 );
                StrandRibbon( &mid3,  &knot2, col2 );
            }
        } else prev = True;

        group = group->gnext;
        captr = next;

        knot1 = knot2;
        col1 = col2;
    }
}


static void DisplaySelected()
{
    register Atom __far *s, __far *d;
    register Chain __far *chain;
    register Group __far *group;
    register Bond __far *bptr;
    register Atom __far *aptr;
    register int irad,sc,dc;
    register int col;

    irad = (int)(Scale*20);

    if( irad>0 )
    {   ForEachBond
        {   s = bptr->srcatom;  
            col = (s->flag&SelectFlag)? 1 : 0;
            sc = Shade2Colour(col);

            d = bptr->dstatom;  
            col = (d->flag&SelectFlag)? 1 : 0;
            dc = Shade2Colour(col);
            ClipCylinder(s->x,s->y,s->z,sc,d->x,d->y,d->z,dc,irad);
        }
    } else ForEachBond
        {   s = bptr->srcatom;  
            col = (s->flag&SelectFlag)? 1 : 0;
            sc = Shade2Colour(col);

            d = bptr->dstatom;  
            col = (d->flag&SelectFlag)? 1 : 0;
            dc = Shade2Colour(col);
            ClipTwinLine(s->x,s->y,s->z,d->x,d->y,d->z,sc,dc);
        }


    irad = (int)(Scale*50);
    ForEachAtom
        if( aptr->flag&NonBondFlag )
        {   col = Shade2Colour( (aptr->flag&SelectFlag)? 1 : 0 );
            ClipSphere(aptr->x,aptr->y,aptr->z,irad,col);
        }
}


void DrawFrame()
{
    register Chain __far *chain;

    if( !Database ) 
        return;

    ClearBuffers();

    if( UseSlabPlane )
    {   SlabValue = (int)(DialValue[7]*ImageRadius)+ImageRadius;
        SlabInten = (int)(ColourMask*TwoOverRootSix);
        SliceValue = SlabValue+16;
        UseClipping = True;
    } else UseClipping = UseScreenClip;

#ifdef IBMPC
    /* Lock Buffers into Memory */
    FBuffer = (Pixel __huge*)GlobalLock(FBufHandle);
    DBuffer = (short __huge*)GlobalLock(DBufHandle);
#endif

    if( !DisplayMode )
    {   if( UseShadow && DrawAtoms )
            if( !VoxelsClean )
                CreateVoxelData();

        if( DrawAtoms ) 
            DisplaySpaceFill();

        if( !UseSlabPlane || (SlabMode != SlabSection) )
        {   if( DrawBonds ) 
                DisplayWireFrame();

            if( DrawRibbon )
                for( chain=Database->clist; chain; chain=chain->cnext )
                    if( chain->glist )
                        DisplayRibbon( chain );

            DisplayHBonds( Database->slist, SSBondMode );
            DisplayHBonds( Database->hlist, HBondMode );
            DisplayBackBone();
        }
    } else DisplaySelected();

#ifdef IBMPC
    /* Unlock Buffers */
    GlobalUnlock(FBufHandle);
    GlobalUnlock(DBufHandle);
#endif
    DBClear = False;
    FBClear = False;
}


/* Identified Atom Info */
static Long IdentDist;
static int IdentFound;
static int IdentDepth;

#if defined(__STDC__) || defined(IBMPC)
/* Function Prototype */
static void TestAtomProximity( Atom __far *, int, int );
#endif

static void TestAtomProximity( ptr, xpos, ypos )
    Atom __far *ptr;
    int xpos, ypos;
{
    register Long dist;
    register int dx,dy;

    if( UseSlabPlane && (ptr->z>SlabValue) )
        return;

    dx = ptr->x - xpos;
    dy = ptr->y - ypos;

    dist = (Long)dx*dx + (Long)dy*dy;

    if( IdentFound )
    {   if( dist==IdentDist )
        {   if( ptr->z<IdentDepth )
                return;
        } else if( dist>IdentDist ) 
            return;
    }

    IdentDepth = ptr->z;
    IdentFound = True;
    IdentDist = dist;
    QAtom = ptr;
}

void IdentifyAtom( xpos, ypos )
    int xpos, ypos;
{
    register int rad, wide, dpth;
    register int new, dx, dy, dz;
    register Chain __far *chain;
    register Group __far *group;
    register HBond __far *hptr;
    register Atom  __far *aptr;
    register Bond __far *bptr;
    register char *str;
    char buffer[40];

    /* Reset Search */
    QChain = (void __far*)0;
    QGroup = (void __far*)0;
    QAtom = (void __far*)0;
    IdentFound = False;

    if( !DisplayMode )
    {   if( !UseSlabPlane || (SlabMode != SlabSection) )
        {   if( DrawBonds )
                ForEachBond
                    if( bptr->flag&DrawBondFlag )
                    {   TestAtomProximity(bptr->srcatom,xpos,ypos);
                        TestAtomProximity(bptr->dstatom,xpos,ypos);
                    }

            ForEachBack
                if( bptr->flag&DrawBondFlag )
                {   TestAtomProximity(bptr->srcatom,xpos,ypos);
                    TestAtomProximity(bptr->dstatom,xpos,ypos);
                }

            for( hptr=Database->hlist; hptr; hptr=hptr->hnext )
                if( hptr->flag )
                {   if( HBondMode )
                    {   TestAtomProximity(hptr->srcCA,xpos,ypos);
                        TestAtomProximity(hptr->dstCA,xpos,ypos);
                    } else
                    {   TestAtomProximity(hptr->src,xpos,ypos);
                        TestAtomProximity(hptr->dst,xpos,ypos);
                    }
                }

            for( hptr=Database->slist; hptr; hptr=hptr->hnext )
                if( hptr->flag )
                {   if( HBondMode )
                    {   TestAtomProximity(hptr->srcCA,xpos,ypos);
                        TestAtomProximity(hptr->dstCA,xpos,ypos);
                    } else
                    {   TestAtomProximity(hptr->src,xpos,ypos);
                        TestAtomProximity(hptr->dst,xpos,ypos);
                    }
                }
        }

        ForEachAtom
        {   /* Identify bond! */
            if( aptr == QAtom )
            {   QChain = chain;
                QGroup = group;
            }

            if( aptr->flag & SphereFlag )
            {   dy = AbsFun(aptr->y-ypos);
                if( dy>aptr->irad ) continue;
                rad = LookUp[aptr->irad][dy];
                dx = AbsFun(aptr->x-xpos);
                if( dx>rad ) continue;

                new = False;
                dpth = aptr->z+LookUp[rad][dx];
                if( UseSlabPlane && (aptr->z+rad>=SlabValue) )
                {   dz = SlabValue-aptr->z;
                    if( SlabMode && (dz >= -rad) )
                    {   wide = LookUp[aptr->irad][AbsFun(dz)];
                        if( (dy<=wide) && (dx<=(int)(LookUp[wide][dy])) )
                        {   if( SlabMode == SlabFinal )
                            {   dpth = SliceValue;
                                new = True;
                            } else if( SlabMode == SlabHollow )
                            {   dpth = aptr->z-LookUp[rad][dx];
                                new = !IdentFound || (dpth>IdentDepth);
                            } else if( SlabMode != SlabHalf )
                            {   /* SlabClose, SlabSection */
                                dpth = dx*dx+dy*dy+dz*dz+SliceValue;
                                if( IdentFound )
                                {   new = (IdentDepth<SliceValue) 
                                          || (dpth<IdentDepth);
                                } else new=True;
                            }
                        } else if( (dz>0) && (SlabMode!=SlabSection) )
                            new = !IdentFound || (dpth>IdentDepth);
                    }
                } else if( !UseSlabPlane || (SlabMode != SlabSection) )
                    new = !IdentFound || (dpth>IdentDepth);

                if( new )
                {   IdentFound = True;
                    IdentDepth = dpth;
                    IdentDist = 0;

                    QChain = chain;
                    QGroup = group;
                    QAtom = aptr;
                }
            } 
        }
    } else /* Display Mode */
    {   ForEachAtom
        {   TestAtomProximity(aptr,xpos,ypos);
            /* Identify bond! */
            if( aptr == QAtom )
            {   QChain = chain;
                QGroup = group;
            }
        }
    }


    if( IdentFound && (IdentDist<50) )
    {   if( CommandActive )
            WriteChar('\n');
        CommandActive = False;

        WriteString("Atom: ");
        str = ElemDesc[QAtom->refno];
        if( str[0]!=' ' )   WriteChar(str[0]);
        WriteChar(str[1]);  WriteChar(str[2]);
        if( str[3]!=' ' )   WriteChar(str[3]);

        sprintf(buffer," %d  ",QAtom->serno);
        WriteString(buffer);

        str = Residue[QGroup->refno];
        if( QAtom->flag&HeteroFlag )
        {   WriteString("Hetero: ");
        } else WriteString("Group: ");

        if( str[0]!=' ' )  WriteChar(str[0]);
        WriteChar(str[1]); WriteChar(str[2]);

        sprintf(buffer," %d",QGroup->serno);
        WriteString(buffer);

        if( QChain->ident!=' ' )
        {   WriteString("  Chain: ");
            WriteChar(QChain->ident);
        }
        WriteChar('\n');
    }
}


void ResetRenderer()
{
    DrawAtoms = False;  MaxAtomRadius = 0;
    DrawBonds = False;  MaxBondRadius = 0;
    DrawRibbon = False;

    SlabMode = SlabClose;
    UseSlabPlane = False;
    UseShadow = False;

    SSBondMode = False;
    HBondMode = False;
    DisplayMode = 0;
    RibbonMode = 0;

    DrawBoundBox = False;
    DrawUnitCell = False;
    DrawAxes = False;

}


void InitialiseRenderer()
{
    register char __far *ptr;
    register int index,rad;
    register int maxval;
    register Real rad2;

    FBuffer = (void __huge*)0;  
    DBuffer = (void __huge*)0;
    IBuffer = (void __far*)0;   
    YBucket = (void __far*)0;

#ifdef IBMPC
    /* Allocate tables on FAR heaps */ 
    Array = (char __far*)_fmalloc(7261*sizeof(char));
    LookUp = (char __far* __far*)_fmalloc(120*sizeof(char __far*));
    HashTable = (void __far* __far*)_fmalloc(VOXSIZE*sizeof(void __far*));
    ColConst = (Card __far*)_fmalloc(120*sizeof(Card));
    
    if( !Array || !LookUp || !HashTable || !ColConst )
        FatalRenderError("tables");    

    FBufHandle = NULL;
    DBufHandle = NULL;
#endif

    ResetRenderer();
    ReSizeScreen();

    ptr = Array;
    for( rad=0; rad<120; rad++ )
    {   rad2 = rad*rad;
        LookUp[rad] = ptr;
        for( index=0; index<=rad; index++ )
#ifdef ISQRT 
            *ptr++ = isqrt( (Card)(rad2-index*index) );
#else 
            *ptr++ = sqrt( (double)(rad2-index*index) );
#endif 

        maxval = (int)(RootSix*rad)+2;
        ColConst[rad] = ((Card)ColourDepth<<ColBits)/maxval;
    }

    FreeItem = (void __far*)0;
    VoxelsClean = False;
    VoxelsDone = False;

    BucketFlag = False;
    RayCount = 0;
}
