/* transfor.c
 * RasMol2 Molecular Graphics
 * Roger Sayle, February 1994
 * Version 2.3
 */
#include "rasmol.h"

#ifdef IBMPC
#include <windows.h>
#endif
#include <stdio.h>
#include <math.h>

#define TRANSFORM
#include "transfor.h"
#include "molecule.h"
#include "command.h"
#include "abstree.h"
#include "render.h"
#include "graphics.h"


typedef struct {
		short col;
		short shade;
                unsigned char r;
                unsigned char g;
                unsigned char b;
              } ShadeRef;


static ShadeRef CPKShade[] = {
     { 0, 0, 200, 200, 200 },       /*  0 Light Grey   */
     { 0, 0, 143, 143, 255 },       /*  1 Sky Blue     */
     { 0, 0, 240,   0,   0 },       /*  2 Red          */
     { 0, 0, 255, 200,  50 },       /*  3 Yellow       */
     { 0, 0, 255, 255, 255 },       /*  4 White        */
     { 0, 0, 255, 192, 203 },       /*  5 Pink         */
     { 0, 0, 218, 165,  32 },       /*  6 Golden Rod   */
     { 0, 0,   0,   0, 255 },       /*  7 Blue         */
     { 0, 0, 255, 165,   0 },       /*  8 Orange       */
     { 0, 0, 128, 128, 144 },       /*  9 Dark Grey    */
     { 0, 0, 165,  42,  42 },       /* 10 Brown        */
     { 0, 0, 160,  32, 240 },       /* 11 Purple       */
     { 0, 0, 255,  20, 147 },       /* 12 Deep Pink    */
     { 0, 0,   0, 255,   0 },       /* 13 Green        */
     { 0, 0, 178,  34,  34 },       /* 14 Fire Brick   */
     { 0, 0,  34, 139,  34 } };     /* 15 Forest Green */


static int CPKIndex[] = { 12, /*Unknown*/
      9, /*Ag*/   9, /*Al*/   6, /*Au*/  13, /* B*/
     10, /*Br*/   0, /* C*/   9, /*Ca*/  13, /*Cl*/
      9, /*Cr*/  10, /*Cu*/   6, /* F*/   8, /*Fe*/
      4, /* H*/   5, /*He*/  11, /* I*/  12, /* K*/
     14, /*Li*/  15, /*Mg*/   9, /*Mn*/   1, /* N*/
      7, /*Na*/  10, /*Ni*/   2, /* O*/   8, /* P*/
      3, /* S*/   6, /*Si*/  10  /*Zn*/ };

static ShadeRef Shapely[] = {
     { 0, 0, 140, 255, 140 },    /* ALA */
     { 0, 0, 255, 255, 255 },    /* GLY */
     { 0, 0,  69,  94,  69 },    /* LEU */
     { 0, 0, 255, 112,  66 },    /* SER */
     { 0, 0, 255, 140, 255 },    /* VAL */
     { 0, 0, 184,  76,   0 },    /* THR */
     { 0, 0,  71,  71, 184 },    /* LYS */
     { 0, 0, 160,   0,  66 },    /* ASP */
     { 0, 0,   0,  76,   0 },    /* ILE */
     { 0, 0, 255, 124, 112 },    /* ASN */
     { 0, 0, 102,   0,   0 },    /* GLU */
     { 0, 0,  82,  82,  82 },    /* PRO */
     { 0, 0,   0,   0, 124 },    /* ARG */
     { 0, 0,  83,  76,  66 },    /* PHE */
     { 0, 0, 255,  76,  76 },    /* GLN */
     { 0, 0, 140, 112,  76 },    /* TYR */
     { 0, 0, 112, 112, 255 },    /* HIS */
     { 0, 0, 255, 255, 112 },    /* CYS */
     { 0, 0, 184, 160,  66 },    /* MET */
     { 0, 0,  79,  70,   0 },    /* TRP */

     { 0, 0, 255,   0, 255 },    /* ASX */
     { 0, 0, 255,   0, 255 },    /* GLX */
     { 0, 0, 255,   0, 255 },    /* PCA */
     { 0, 0, 255,   0, 255 },    /* HYP */

     { 0, 0, 160, 160, 255 },    /*   A */
     { 0, 0, 255, 140,  75 },    /*   C */
     { 0, 0, 255, 112, 112 },    /*   G */
     { 0, 0, 160, 255, 160 },    /*   T */

     { 0, 0, 184, 184, 184 },    /* 28 -> BackBone */
     { 0, 0,  94,   0,  94 },    /* 29 -> Special  */
     { 0, 0, 255,   0, 255 } };  /* 30 -> Default  */

     
static ShadeRef AminoShade[] = {
     { 0, 0, 230,  10,  10 },    /*  0  ASP, GLU      */
     { 0, 0,  20,  90, 255 },    /*  1  LYS, ARG      */
     { 0, 0, 130, 130, 210 },    /*  2  HIS           */
     { 0, 0, 250, 150,   0 },    /*  3  SER, THR      */
     { 0, 0,   0, 220, 220 },    /*  4  ASN, GLN      */
     { 0, 0, 230, 230,   0 },    /*  5  CYS, MET      */
     { 0, 0, 200, 200, 200 },    /*  6  ALA           */
     { 0, 0, 235, 235, 235 },    /*  7  GLY           */
     { 0, 0,  15, 130,  15 },    /*  8  LEU, VAL, ILE */
     { 0, 0,  50,  50, 170 },    /*  9  PHE, TYR      */
     { 0, 0, 180,  90, 180 },    /* 10  TRP           */
     { 0, 0, 220, 150, 130 },    /* 11  PRO, PCA, HYP */
     { 0, 0, 190, 160, 110 } };  /* 12  Others        */

static int AminoIndex[] = {
      6, /*ALA*/   7, /*GLY*/   8, /*LEU*/   3,  /*SER*/
      8, /*VAL*/   3, /*THR*/   1, /*LYS*/   0,  /*ASP*/
      8, /*ILE*/   4, /*ASN*/   0, /*GLU*/  11,  /*PRO*/
      1, /*ARG*/   9, /*PHE*/   4, /*GLN*/   9,  /*TYR*/
      2, /*HIS*/   5, /*CYS*/   5, /*MET*/  10,  /*TRP*/
      4, /*ASX*/   4, /*GLX*/  11, /*PCA*/  11   /*HYP*/
                          };

static ShadeRef HBondShade[] = {
     { 0, 0, 255, 255, 255 },    /* 0  Offset =  2   */
     { 0, 0, 255,   0, 255 },    /* 1  Offset =  3   */
     { 0, 0, 255,   0,   0 },    /* 2  Offset =  4   */
     { 0, 0, 255, 165,   0 },    /* 3  Offset =  5   */
     { 0, 0,   0, 255, 255 },    /* 4  Offset = -3   */
     { 0, 0,   0, 255,   0 },    /* 5  Offset = -4   */
     { 0, 0, 255, 255,   0 } };  /* 6  Others        */


static ShadeRef StructShade[] = {
     { 0, 0, 255, 255, 255 },    /* 0  Default     */
     { 0, 0, 240,   0, 128 },    /* 1  Alpha Helix */
     { 0, 0, 255, 255,   0 },    /* 2  Beta Sheet  */
     { 0, 0,  96, 128, 255 } };  /* 3  Turn        */


typedef struct { 
                Long refcount;
                unsigned char r;
                unsigned char g;
                unsigned char b;
              } ShadeDesc;


/* Macros for commonly used loops */
#define ForEachAtom  for(chain=Database->clist;chain;chain=chain->cnext) \
                     for(group=chain->glist;group;group=group->gnext)    \
                     for(ptr=group->alist;ptr;ptr=ptr->anext)
#define ForEachBond  for(bptr=Database->blist;bptr;bptr=bptr->bnext) 
#define ForEachBack  for(chain=Database->clist;chain;chain=chain->cnext) \
                     for(bptr=chain->blist;bptr;bptr=bptr->bnext)

#define MAXSHADE         32
#define MatchChar(a,b)   (((a)=='#')||((a)==(b)))
#define RootSix          2.44948974278


static ShadeDesc Shade[MAXSHADE];
static ShadeRef ScaleRef[MAXSHADE];
static int MaskColour[MAXMASK];
static int MaskShade[MAXMASK];
static int ScaleCount;
static int LastShade;

static Real LastRX,LastRY,LastRZ;
static Long OrigCX,OrigCY,OrigCZ;
static Real Zoom;



void DetermineClipping()
{
    register int temp;
    register int max;

    max = 0;
    if( DrawAtoms && (MaxAtomRadius>max) )  max = MaxAtomRadius;
    if( DrawBonds && (MaxBondRadius>max) )  max = MaxBondRadius;
       
    temp = ImageRadius + max;
    UseScreenClip = (XOffset<temp) || (XOffset+temp>=XRange) ||
                    (YOffset<temp) || (YOffset+temp>=YRange);
}



void SetRadiusValue( rad )
    int rad;
{
    register int irad,change;
    register Chain __far *chain;
    register Group __far *group;
    register Atom __far *ptr;

    if( !Database )
        return;

    irad = (int)(Scale*rad);
    MaxAtomRadius = 0;
    DrawAtoms = False;
    change = False;

    ForEachAtom
        if( ptr->flag & SelectFlag )
        {   if( irad>MaxAtomRadius )
                MaxAtomRadius = irad;
            ptr->flag |= SphereFlag;
            ptr->radius = rad;
            ptr->irad = irad;
            change = True;
        } else if( ptr->flag & SphereFlag )
        {   DrawAtoms = True;
            if( ptr->irad>MaxAtomRadius )
                MaxAtomRadius = ptr->irad;
        }

    if( change )
    {   DrawAtoms = True;
        DetermineClipping();
        VoxelsClean = False;
        BucketFlag = False;
    }
}

void SetRadiusTemperature()
{
    register int rad,irad,change;
    register Chain __far *chain;
    register Group __far *group;
    register Atom __far *ptr;

    if( !Database )
        return;

    MaxAtomRadius = 0;
    DrawAtoms = False;
    change = False;

    ForEachAtom
        if( (ptr->flag&SelectFlag) && (ptr->temp>0) )
        {   rad = (5*ptr->temp)>>1;
            if( rad>500 ) rad = 500;

            irad = (int)(Scale*rad);
            if( irad>MaxAtomRadius )
                MaxAtomRadius = irad;
            ptr->flag |= SphereFlag;
            ptr->radius = rad;
            ptr->irad = irad;
            change = True;
        } else if( ptr->flag & SphereFlag )
        {   DrawAtoms = True;
            if( ptr->irad>MaxAtomRadius )
                MaxAtomRadius = ptr->irad;
        }

    if( change )
    {   DrawAtoms = True;
        DetermineClipping();
        VoxelsClean = False;
        BucketFlag = False;
    }
}


void SetVanWaalRadius()
{
    register int rad,change;
    register Chain __far *chain;
    register Group __far *group;
    register Atom __far *ptr;

    if( !Database )
        return;

    MaxAtomRadius = 0;
    DrawAtoms = False;
    change = False;

    ForEachAtom
        if( ptr->flag&SelectFlag )
        {   rad = VanWaalRadius[ GetElemIdent(ptr) ];
            ptr->irad = (int)(Scale*rad);
            ptr->radius = rad;
            change = True;

            ptr->flag |=SphereFlag;
            if( ptr->irad>MaxAtomRadius )
                MaxAtomRadius = ptr->irad;
        } else if( ptr->flag&SphereFlag )
        {   DrawAtoms = True;
            if( ptr->irad>MaxAtomRadius )
                MaxAtomRadius = ptr->irad;
        }

    if( change )
    {   DrawAtoms = True;
        DetermineClipping();
        VoxelsClean = False;
        BucketFlag = False;
    }
}


void DisableSpacefill()
{
    register Chain __far *chain;
    register Group __far *group;
    register Atom __far *ptr;

    if( !Database || !DrawAtoms )
        return;

    MaxAtomRadius = 0;
    DrawAtoms = False;
    
    ForEachAtom
        if( !(ptr->flag&SelectFlag) )
        {   if( ptr->flag&SphereFlag )
            {   if( ptr->irad>MaxAtomRadius )
                    MaxAtomRadius = ptr->irad;
                DrawAtoms = True;
            }
        } else if( ptr->flag&SphereFlag )
            ptr->flag &= ~SphereFlag;

    DetermineClipping();
    VoxelsClean = False;
    BucketFlag = False;
}



void EnableWireFrame( depth, rad )
    int depth, rad;
{
    register Bond __far *bptr;
    register int flag, irad;

    if( !Database )
        return;

    DrawBonds = False;
    MaxBondRadius = 0;
    irad = (int)(Scale*rad);

    ForEachBond
    {   flag = ZoneBoth? bptr->dstatom->flag & bptr->srcatom->flag
                       : bptr->dstatom->flag | bptr->srcatom->flag;

        if( flag&SelectFlag )
        {   DrawBonds = True;
            bptr->flag &= ~DrawBondFlag;
            if( !depth )
            {   if( irad>MaxBondRadius )
                    MaxBondRadius = irad;
                bptr->flag |= CylinderFlag;
                bptr->radius = rad;
                bptr->irad = irad;
            } else bptr->flag |= WireFlag;
        } else if( bptr->flag&DrawBondFlag )
        {    DrawBonds = True;
             if( bptr->flag&CylinderFlag )
                 if( bptr->irad>MaxBondRadius )
                     MaxBondRadius = bptr->irad;
        }
    }
    DetermineClipping();
}


void DisableWireFrame()
{
    register Bond __far *bptr;
    register int flag;

    if( !Database || !DrawBonds )
        return;

    DrawBonds = False;
    MaxBondRadius = 0;

    ForEachBond
    {   flag = ZoneBoth? bptr->dstatom->flag & bptr->srcatom->flag
                       : bptr->dstatom->flag | bptr->srcatom->flag;

        if( flag&SelectFlag )
        {   bptr->flag &= ~DrawBondFlag;
        } else if( bptr->flag&DrawBondFlag )
        {   DrawBonds = True;
            if( bptr->flag&CylinderFlag )
                if( bptr->irad>MaxBondRadius )
                    MaxBondRadius = bptr->irad;
        }
    }
    DetermineClipping();
}


void EnableBackBone( depth, rad )
    int depth, rad;
{
    register Chain __far *chain;
    register Bond __far *bptr;
    register int flag,irad;

    if( !Database )
        return;

    irad = (int)(Scale*rad);

    ForEachBack
    {   flag = ZoneBoth? bptr->dstatom->flag & bptr->srcatom->flag
                       : bptr->dstatom->flag | bptr->srcatom->flag;

        if( flag&SelectFlag )
        {   bptr->flag &= ~DrawBondFlag;
            if( !depth )
            {   bptr->flag |= CylinderFlag;
                bptr->radius = rad;
                bptr->irad = irad;
            } else bptr->flag |= WireFlag;
        } 
    }
    DetermineClipping();
}


void DisableBackBone()
{
    register Chain __far *chain;
    register Bond __far *bptr;

    if( !Database )
        return;

    if( ZoneBoth )
    {   ForEachBack
            if( (bptr->dstatom->flag&bptr->srcatom->flag) & SelectFlag )
                bptr->flag &= ~DrawBondFlag;
    } else ForEachBack
        if( (bptr->dstatom->flag|bptr->srcatom->flag) & SelectFlag )
            bptr->flag &= ~DrawBondFlag;
    DetermineClipping();
}


void SetHBondStatus( hbonds, enable, rad )
    int hbonds, enable, rad;
{
    register HBond __far *list;
    register HBond __far *ptr;
    register Atom __far *src;
    register Atom __far *dst;
    register int flag, irad;

    if( !Database )
        return;

    if( hbonds )
    {   if( InfoHBondCount<0 )
            CalcHydrogenBonds();
        list = Database->hlist;
    } else 
    {   if( InfoSSBondCount<0 )
            FindDisulphideBridges();
        list = Database->slist;
    }

    irad = (int)(Scale*rad);
    for( ptr=list; ptr; ptr=ptr->hnext )
    {    if( !hbonds && SSBondMode )
         {   src = ptr->srcCA;
             dst = ptr->dstCA;
         } else
         {   src = ptr->src;
             dst = ptr->dst;
         }

         flag = ZoneBoth? src->flag&dst->flag : src->flag|dst->flag;
         if( flag & SelectFlag ) 
         {   ptr->flag &= ~DrawBondFlag;
             if( enable )
             {   if( rad )
                 {   ptr->flag |= CylinderFlag;
                     ptr->radius = rad;
                     ptr->irad = irad;
                 } else ptr->flag |= WireFlag;
             }
         }
    }
}


void SetRibbonStatus( enable, width )
    int enable, width;
{
    register Chain __far *chain;
    register Group __far *group;
    register Atom __far *ptr;

    if( !Database )
        return;

    /* Ribbons already disabled! */
    if( !enable && !DrawRibbon )
        return;

    if( InfoHelixCount<0 )
        DetermineStructure();

    DrawRibbon = False;
    for( chain=Database->clist; chain; chain=chain->cnext )
        for( group=chain->glist; group; group=group->gnext )
            if( enable )
            {   if( group->flag & RibbonFlag )
                    DrawRibbon = True;
                
                for( ptr=group->alist; ptr; ptr=ptr->anext )
                    if( IsAlphaCarbon(ptr->refno) )
                    {   if( ptr->flag&SelectFlag )
                        {   group->flag |= RibbonFlag;
                            if( !width )
                            {   if( group->flag & (HelixFlag|SheetFlag) )
                                {      group->width = 380;
                                } else group->width = 100;
                            } else group->width = width;
                            DrawRibbon = True;
                        }
                        break;

                    } else if( IsSugarPhosphate(ptr->refno) )
                    {   if( ptr->flag&SelectFlag )
                        {   group->width = width? width : 720;
                            group->flag |= RibbonFlag;
                            DrawRibbon = True;
                        }
                        break;
                    }


            } else  /* Disable Ribbon */
                if( group->flag & RibbonFlag )
                {   for( ptr=group->alist; ptr; ptr=ptr->anext )
                        if( IsAlphaCarbon(ptr->refno) ||
                            IsSugarPhosphate(ptr->refno) )
                        {   if( ptr->flag&SelectFlag )
                                group->flag &= ~RibbonFlag;
                            break;
                        }
                    if( group->flag & RibbonFlag ) 
                        DrawRibbon = True;
                }
}


void SelectZone( mask )
    int mask;
{
    register Bond __far *bptr;
    register Chain __far *chain;
    register Group __far *group;
    register Atom __far *ptr;

    if( !Database )
        return;

    SelectCount = 0;
    ForEachAtom
        if( ptr->flag & mask )
        {   ptr->flag |= SelectFlag;
            SelectCount++;
        } else ptr->flag &= ~SelectFlag;
    DisplaySelectCount();

    if( ZoneBoth )
    {   ForEachBond
           if( (bptr->srcatom->flag&bptr->dstatom->flag) & SelectFlag )
           {   bptr->flag |= SelectFlag;
           } else bptr->flag &= ~SelectFlag;
    } else
        ForEachBond
           if( (bptr->srcatom->flag|bptr->dstatom->flag) & SelectFlag )
           {   bptr->flag |= SelectFlag;
           } else bptr->flag &= ~SelectFlag;

}


void RestrictZone( mask )
    int mask;
{
    register Bond __far *bptr;
    register Chain __far *chain;
    register Group __far *group;
    register Atom __far *ptr;
    register int flag;

    if( !Database )
        return;

    DrawAtoms = False;   MaxAtomRadius = 0;
    DrawBonds = False;   MaxBondRadius = 0;
    
    SelectCount = 0;
    ForEachAtom
        if( ptr->flag & mask )
        {   ptr->flag |= SelectFlag;
            SelectCount++;

            if( ptr->flag & SphereFlag )
            {   DrawAtoms = True;
                if( ptr->irad>MaxAtomRadius )
                    MaxAtomRadius = ptr->irad;
            }
        } else ptr->flag &= ~(SelectFlag|SphereFlag);
    DisplaySelectCount();
    
    ForEachBond
    {   /* Ignore ZoneBoth setting! */
        flag = bptr->dstatom->flag & bptr->srcatom->flag;
        if( flag & SelectFlag )
        {   bptr->flag |= SelectFlag;
            if( bptr->flag&DrawBondFlag )
            {   DrawBonds = True;
                if( bptr->flag & CylinderFlag )
                    if( bptr->irad>MaxBondRadius )
                        MaxBondRadius = bptr->irad;
            } 
        } else bptr->flag &= ~(SelectFlag|DrawBondFlag);
    }

    DetermineClipping();
    VoxelsClean = False;
    BucketFlag = False;
}



int DefineShade( r, g, b )
    unsigned char r, g, b;
{
    register int d,dr,dg,db;
    register int dist,best;
    register int i;

    /* Already defined! */
    for( i=0; i<LastShade; i++ )
        if( Shade[i].refcount )
            if( (Shade[i].r==r)&&(Shade[i].g==g)&&(Shade[i].b==b) )
                return(i);

    /* Allocate request */
    for( i=0; i<LastShade; i++ )
         if( !Shade[i].refcount )
         {   Shade[i].r = r;
             Shade[i].g = g;
             Shade[i].b = b;
             Shade[i].refcount = 0;
             return(i);
         }

    if( CommandActive )
        WriteChar('\n');
    WriteString("Warning: Unable to allocate shade!\n");
    CommandActive = False;

    /* To avoid lint warning! */
    best = dist = 0;

    /* Nearest match */
    for( i=0; i<LastShade; i++ )
    {   dr = Shade[i].r - r;
        dg = Shade[i].g - g;
        db = Shade[i].b - b;
        d = dr*dr + dg*dg + db*db;
        if( !i || (d<dist) )
        {   dist = d;
            best = i;
        }
    }
    return( best );
}


void ScaleColourMap( count )
    int count;
{
    register Real hue;
    register int fract, sextant;
    register int i, r, g, b;

    ScaleCount=0;
    for( i=0; i<LastShade; i++ )
        if( !Shade[i].refcount )
            ScaleCount++;

    /* If there are no shades free! */
    if( !ScaleCount ) ScaleCount = LastShade;

    if( count && (count<ScaleCount) )
        ScaleCount = count;

    for( i=0; i<ScaleCount; i++ )
    {   sextant = (int)(hue = (4.0*i)/(ScaleCount-1));
        fract = (int)(255.0*(hue-sextant));

        switch( sextant )
        {   case(0): r = 0;         g = fract;     b = 255;         break;
            case(1): r = 0;         g = 255;       b = 255-fract;   break;
            case(2): r = fract;     g = 255;       b = 0;           break;
            case(3): r = 255;       g = 255-fract; b = 0;           break;
            default: r = 255;       g = 0;         b = 0;
        }
        ScaleRef[i].r = r;
        ScaleRef[i].g = g;
        ScaleRef[i].b = b;
        ScaleRef[i].shade = 0;
        ScaleRef[i].col = 0;
    }
}


static void SetLutEntry( i, r, g, b )
    int i, r, g, b;
{
    ULut[i] = True;
    RLut[i] = r;
    GLut[i] = g;
    BLut[i] = b;

#ifdef EIGHTBIT
    Lut[i] = i;
#else
    Lut[i] = ( ((r<<8)|g)<<8 ) | b;
#endif
}


static Real Power( x, y )
    Real x; int y;
{
    register Real result;

    result = x;
    while( y>1 )
    {   if( y&1 ) { result *= x; y--; }
        else { result *= result; y>>=1; }
    }
    return( result );
}


void DefineColourMap()
{
    register Real diffuse;
    register Real temp, inten;
    register int col, r, g, b;
    register int i, j, k;

#ifdef EIGHTBIT
    for( i=0; i<256; i++ )
        ULut[i] = False;
#endif

#ifdef IBMPC
    SetLutEntry(0,BackR,BackG,BackB);
#else
    SetLutEntry(5,BackR,BackG,BackB);
#endif

    diffuse = 1.0 - Ambient;
    for( i=0; i<ColourDepth; i++ )
    {   temp = (Real)i/ColourMask;
        inten = diffuse*temp + Ambient;

        if( DisplayMode )
        {   /* Unselected [0,96,255] */
            /* Selected   [255,0,0]  */
            r = b = (int)(255*inten);
            g = (int)(96*inten);

            SetLutEntry( Shade2Colour(0)+i, 0, g, b );
            SetLutEntry( Shade2Colour(1)+i, r, 0, 0 );

        } else 
        {   if( FakeSpecular )
            {   temp = Power(temp,SpecPower);
                inten *= 1.0 - temp;
                k = (int)(255*temp);
            }

            for( j=0; j<LastShade; j++ )
                if( Shade[j].refcount )
                {   col = Shade2Colour(j);
                    r = (int)(Shade[j].r*inten); 
                    g = (int)(Shade[j].g*inten);
                    b = (int)(Shade[j].b*inten);

                    if( FakeSpecular )
                    {   r += k;
                        g += k;
                        b += k;
                    }
                    SetLutEntry( col+i, r, g, b );
                }
        }
    }

    if( Interactive )
        AllocateColourMap();
}


void ResetColourMap()
{
    register int i;

#ifdef EIGHTBIT
    for( i=0; i<256; i++ )
        ULut[i] = False;
#endif

    BackR = BackG = BackB = 0;

    for( i=0; i<LastShade; i++ )
        Shade[i].refcount = 0;
    ScaleCount = 0;
}


void ColourBondNone()
{
    register Bond __far *bptr;

    if( Database )
        ForEachBond
            if( (bptr->flag&SelectFlag) && bptr->col )
            {   Shade[Colour2Shade(bptr->col)].refcount--;
                bptr->col = 0;
            }
}


void ColourBondAttrib( r, g, b )
    int r, g, b;
{
    register Bond __far *bptr;
    register int shade,col;

    if( Database )
    {   ForEachBond
            if( (bptr->flag&SelectFlag) && bptr->col )
                Shade[Colour2Shade(bptr->col)].refcount--;

        shade = DefineShade((Byte)r,(Byte)g,(Byte)b);
        col = Shade2Colour(shade);

        ForEachBond
            if( bptr->flag&SelectFlag )
            {   Shade[shade].refcount++;
                bptr->col = col;
            }
    }
}


void ColourBackNone()
{
    register Chain __far *chain;
    register Bond __far *bptr;
    register int flag;

    if( Database )
        ForEachBack
        {   flag = ZoneBoth? bptr->dstatom->flag & bptr->srcatom->flag
                           : bptr->dstatom->flag | bptr->srcatom->flag;

            if( flag&SelectFlag )
            {   bptr->flag |= SelectFlag;
                if( bptr->col )
                {   Shade[Colour2Shade(bptr->col)].refcount--;
                    bptr->col = 0;
                }
            } else bptr->flag &= ~SelectFlag;
        }
}


void ColourBackAttrib( r, g, b )
    int r, g, b;
{
    register int shade,col;
    register Chain __far *chain;
    register Bond __far *bptr;

    if( Database )
    {   ColourBackNone();
        shade = DefineShade((Byte)r,(Byte)g,(Byte)b);
        col = Shade2Colour(shade);

        ForEachBack
            if( bptr->flag&SelectFlag )
            {   Shade[shade].refcount++;
                bptr->col = col;
            }
    }
}


void ColourHBondNone( hbonds )
    int hbonds;
{
    register HBond __far *list;
    register HBond __far *ptr;
    register Atom __far *src;
    register Atom __far *dst;

    if( !Database )
        return;

    if( hbonds )
    {   if( InfoHBondCount<0 )
        {   CalcHydrogenBonds();
            return;
        } else list = Database->hlist;
    } else
        if( InfoSSBondCount<0 )
        {   FindDisulphideBridges();
            return;
        } else list = Database->slist;


    if( ZoneBoth )
    {   for( ptr=list; ptr; ptr=ptr->hnext )
        {   if( !hbonds && SSBondMode )
            {   src = ptr->srcCA;
                dst = ptr->dstCA;
            } else
            {   src = ptr->src;
                dst = ptr->dst;
            }

            if( (src->flag&dst->flag) & SelectFlag )
            {   ptr->flag |= SelectFlag;
                if( ptr->col )
                {   Shade[Colour2Shade(ptr->col)].refcount--;
                    ptr->col = 0;
                }
            } else ptr->flag &= ~SelectFlag;
        }
    } else
        for( ptr=list; ptr; ptr=ptr->hnext )
        {   if( !hbonds && SSBondMode )
            {   src = ptr->srcCA;
                dst = ptr->dstCA;
            } else
            {   src = ptr->src;
                dst = ptr->dst;
            }

            if( (src->flag|dst->flag) & SelectFlag )
            {   ptr->flag |= SelectFlag;
                if( ptr->col )
                {   Shade[Colour2Shade(ptr->col)].refcount--;
                    ptr->col = 0;
                }
            } else ptr->flag &= ~SelectFlag;
        }
}

void ColourHBondType()
{
    register HBond __far *ptr;
    register ShadeRef *ref;
    register int i;

    if( !Database ) return;
    for( i=0; i<7; i++ )
        HBondShade[i].col = 0;
    ColourHBondNone( True );

    for( ptr=Database->hlist; ptr; ptr=ptr->hnext )
        if( ptr->flag & SelectFlag )
        {   switch( ptr->offset )
            {   case(  2 ):  ref = HBondShade;     break;
                case(  3 ):  ref = HBondShade+1;   break;
                case(  4 ):  ref = HBondShade+2;   break;
                case(  5 ):  ref = HBondShade+3;   break;
                case( -3 ):  ref = HBondShade+4;   break;
                case( -4 ):  ref = HBondShade+5;   break;
                default:     ref = HBondShade+6;   break;
            }

            if( !ref->col )
            {   ref->shade = DefineShade( ref->r, ref->g, ref->b );
                ref->col = Shade2Colour(ref->shade);
            }
            Shade[ref->shade].refcount++;
            ptr->col = (Byte)ref->col;
        }
}


void ColourHBondAttrib( hbonds, r, g, b )
    int r, g, b;
{
    register HBond __far *list;
    register HBond __far *ptr;
    register int col,shade;

    if( !Database )
        return;

    ColourHBondNone( hbonds );
    shade = DefineShade((Byte)r,(Byte)g,(Byte)b);
    col = Shade2Colour(shade);

    list = hbonds? Database->hlist : Database->slist;
    for( ptr=list; ptr; ptr=ptr->hnext )
        if( ptr->flag & SelectFlag )
        {   Shade[shade].refcount++;
            ptr->col = col;
        }
}


void ColourRibbonNone()
{
    register Chain __far *chain;
    register Group __far *group;
    register Atom __far *aptr;

    if( !Database )
        return;

    if( InfoHelixCount<0 )
    {   DetermineStructure();
        return;
    }

    for( chain=Database->clist; chain; chain=chain->cnext )
        for( group=chain->glist; group; group=group->gnext )
            if( group->col && (aptr=FindGroupAtom(group,1) ) 
                && (aptr->flag&SelectFlag) )
            {   Shade[Colour2Shade(group->col)].refcount--;
                group->col = 0;
            }
}


void ColourRibbonAttrib( r, g, b )
    int r, g, b;
{
    register int shade, col;
    register Chain __far *chain;
    register Group __far *group;
    register Atom __far *aptr;

    if( Database )
    {   ColourRibbonNone();
        shade = DefineShade((Byte)r,(Byte)g,(Byte)b);
        col = Shade2Colour(shade);

        for( chain=Database->clist; chain; chain=chain->cnext )
            for( group=chain->glist; group; group=group->gnext )
                if( (aptr=FindGroupAtom(group,1)) 
                    && (aptr->flag&SelectFlag) )
                {   Shade[shade].refcount++;
                    group->col = col;
                }
    }
}


static void ResetColourAttrib()
{
    register Chain __far *chain;
    register Group __far *group;
    register Atom __far *ptr;

    ForEachAtom
        if( (ptr->flag&SelectFlag) && ptr->col )
            Shade[Colour2Shade(ptr->col)].refcount--;
}


void MonoColourAttrib( r, g, b )
    int r, g, b;
{
    register int shade,col;
    register Chain __far *chain;
    register Group __far *group;
    register Atom __far *ptr;

    if( Database )
    {   ResetColourAttrib();
        shade = DefineShade((Byte)r,(Byte)g,(Byte)b);
        col = Shade2Colour(shade);

        ForEachAtom
            if( ptr->flag&SelectFlag )
            {   Shade[shade].refcount++;
                ptr->col = col;
            }
    }
}


void ScaleColourAttrib( attr )
    int attr;
{
    register ShadeRef *ref;
    register int count, attrno, factor;
    register Chain __far *chain;
    register Group __far *group;
    register Atom __far *ptr;

    if( !Database ) return;

    switch( attr )
    {   case(ChainAttr):   attrno = InfoChainCount;   
                           break;

        case(GroupAttr):   factor = MinMainRes;
                           attrno = MaxMainRes;
                           if( HetaGroups && HetaGroupCount )
                           {   if( MinHetaRes < factor )
                                   factor = MinHetaRes;
                               if( MaxHetaRes > attrno )
                                   attrno = MaxHetaRes;
                           } 
                           attrno -= factor;
                           break;

        case(TempAttr):    factor = MinMainTemp;
                           attrno = MaxMainTemp;
                           if( HetaGroups && HetaGroupCount )
                           {   if( MinHetaTemp < factor )
                                   factor = MinHetaTemp;
                               if( MaxHetaTemp > attrno )
                                   attrno = MaxHetaTemp;
                           }
                           attrno -= factor;
                           break;

        default:           return;
    }

    if( attrno<2 )
    {   if( CommandActive )
            WriteChar('\n');
        WriteString("Warning: Only a single attribute value!\n");
        CommandActive = False;
        return;
    }

    ResetColourAttrib();
    ScaleColourMap(attrno);

    switch( attr )
    {    case(ChainAttr):
                 count = 0;
                 for( chain=Database->clist; chain; chain=chain->cnext )
                 {   ref = &(ScaleRef[(count*ScaleCount)/attrno]);
                     if( !(ref->col && Shade[ref->shade].refcount) )
                     {   ref->shade = DefineShade(ref->r,ref->g,ref->b);
                         ref->col = Shade2Colour(ref->shade);
                     }
                     for( group=chain->glist; group; group=group->gnext )
                         for( ptr=group->alist; ptr; ptr=ptr->anext )
                             if( ptr->flag&SelectFlag )
                             {   Shade[ref->shade].refcount++;
                                 ptr->col = ref->col;
                             }
                     count++;
                 }
                 break;


         case(GroupAttr):
                 factor++;
                 for( chain=Database->clist; chain; chain=chain->cnext )
                     for( group=chain->glist; group; group=group->gnext )
                     {   count = group->serno-factor;
                         if( count<0 )
                         {   ref = ScaleRef;
                         } else if( count<attrno )
                         {   count = (ScaleCount*count)/attrno;
                             ref = ScaleRef + count;
                         } else ref = ScaleRef + (ScaleCount-1);

                         if( !(ref->col && Shade[ref->shade].refcount) )
                         {   ref->shade = DefineShade(ref->r,ref->g,ref->b);
                             ref->col = Shade2Colour(ref->shade);
                         }

                         for( ptr=group->alist; ptr; ptr=ptr->anext )
                             if( ptr->flag&SelectFlag )
                             {   Shade[ref->shade].refcount++;
                                 ptr->col = ref->col;
                             }
                     }
                 break;


         case(TempAttr):
                 factor++;
                 ForEachAtom
                     if( ptr->flag&SelectFlag )
                     {   count = ptr->temp-factor;
                         if( count<0 )
                         {   ref = ScaleRef;
                         } else if( count<attrno )
                         {   count = (ScaleCount*count)/attrno;
                             ref = ScaleRef + count;
                         } else ref = ScaleRef + (ScaleCount-1);

                         if( !ref->col )
                         {   ref->shade = DefineShade(ref->r,ref->g,ref->b);
                             ref->col = Shade2Colour(ref->shade);
                         }
                         Shade[ref->shade].refcount++;
                         ptr->col = ref->col;
                     }
                 break;
    }
}



static int MatchNumber( len, value, mask )
    int len, value;
    char *mask;
{
    register char digit, template;
    register int result;
    register int i;

    result = True;
    for( i=0; i<len; i++ )
    {   digit = (value%10) + '0';
        template = mask[len-i];
        if( template==' ' )
        {   if( value ) result = False;
        } else if( !MatchChar(template,digit) )
            result = False;
        value/=10;
    }
    return( result );
}


void UserMaskAttrib( fields )
    int fields;
{
    register MaskDesc *mptr;
    register char *temp, *name;
    register int shade, change;
    register int i, rad, match;

    register Chain __far *chain;
    register Group __far *group;
    register Atom __far *ptr;


    if( !Database ) return;

    if( !MaskCount )
    {   if( CommandActive )
            WriteChar('\n');
        WriteString("Warning: No user supplied colour records!\n");
        CommandActive = False;
        return;
    }

    change = False;
    ResetColourAttrib();
    if( fields&MaskColourFlag )
        for( i=0; i<MaskCount; i++ )
            MaskShade[i] = -1;

    if( fields&MaskRadiusFlag )
    {   MaxAtomRadius = 0;
        DrawAtoms = False;
    }


    ForEachAtom
    if( ptr->flag&SelectFlag )
    {   for( i=0; i<MaskCount; i++ )
        {   mptr = UserMask+i;
            temp = mptr->mask;
            match = True;

            if( !MatchChar(temp[13],chain->ident) ) match=False;
            if( !MatchChar(temp[9],ptr->altl) )     match=False;

            /* Atom Name */
            if( match )
            {   name = ElemDesc[ptr->refno];
                if( !MatchChar(temp[5],name[0]) ) match=False;
                if( !MatchChar(temp[6],name[1]) ) match=False;
	        if( !MatchChar(temp[7],name[2]) ) match=False;
	        if( !MatchChar(temp[8],name[3]) ) match=False;
            }

            /* Group Name */
            if( match )
            {   name = Residue[group->refno];
                if( !MatchChar(temp[10],name[0]) ) match=False;
                if( !MatchChar(temp[11],name[1]) ) match=False;
                if( !MatchChar(temp[12],name[2]) ) match=False;
            }


            if( match && (mptr->flags&SerNoFlag) )
                match = MatchNumber(4,ptr->serno,&temp[0]);
            if( match && (mptr->flags&ResNoFlag) )
                match = MatchNumber(3,group->serno,&temp[14]);
            if( match ) break;
        }

        if( fields&MaskColourFlag )
        {   if( match )
            {   if( MaskShade[i] == -1 )
                {   MaskShade[i] = DefineShade(mptr->r,mptr->g,mptr->b);
                    MaskColour[i] = Shade2Colour(MaskShade[i]);
                }
                Shade[MaskShade[i]].refcount++;
                ptr->col = MaskColour[i];
            } else
            {   shade = DefineShade(255,255,255);
                ptr->col = Shade2Colour(shade);
                Shade[shade].refcount++;
            }
        }

        if( fields&MaskRadiusFlag )
        {   rad = match? mptr->radius : 375;
            ptr->irad = (int)(Scale*rad);
            ptr->flag |= SphereFlag;
            ptr->radius = rad;

            if( ptr->irad>MaxAtomRadius )
                MaxAtomRadius = ptr->irad;
            change = True;
        }
    } else if( ptr->flag&SphereFlag )
    {   DrawAtoms = True;
        if( ptr->irad>MaxAtomRadius )
            MaxAtomRadius = ptr->irad;
    }

    if( change )
    {   DrawAtoms = True;
        DetermineClipping();
        VoxelsClean = False;
        BucketFlag = False;
    }
}


void CPKColourAttrib()
{
    register ShadeRef *ref;
    register Chain __far *chain;
    register Group __far *group;
    register Atom __far *ptr;
    register int i;

    if( !Database ) return;
    for( i=0; i<7; i++ )
        CPKShade[i].col = 0;
    ResetColourAttrib();


    ForEachAtom
        if( ptr->flag&SelectFlag )
        {   i = GetElemIdent( ptr );
            ref = CPKShade + CPKIndex[i];

            if( !ref->col )
            {   ref->shade = DefineShade( ref->r, ref->g, ref->b );
                ref->col = Shade2Colour(ref->shade);
            }
            Shade[ref->shade].refcount++;
            ptr->col = ref->col;
        }
}


void AminoColourAttrib()
{
    register ShadeRef *ref;
    register Chain __far *chain;
    register Group __far *group;
    register Atom __far *ptr;
    register int i;

    if( !Database ) return;
    for( i=0; i<13; i++ )
        AminoShade[i].col = 0;
    ResetColourAttrib();

    ForEachAtom
        if( ptr->flag&SelectFlag )
        {   if( IsAmino(group->refno) )
            {   ref = AminoShade + AminoIndex[group->refno];
            } else ref = AminoShade+12;

            if( !ref->col )
            {   ref->shade = DefineShade( ref->r, ref->g, ref->b );
                ref->col = Shade2Colour(ref->shade);
            }
            Shade[ref->shade].refcount++;
            ptr->col = ref->col;
        }
}


void ShapelyColourAttrib()
{
    register ShadeRef *ref;
    register Chain __far *chain;
    register Group __far *group;
    register Atom __far *ptr;
    register int i;

    if( !Database ) return;
    for( i=0; i<30; i++ )
        Shapely[i].col = 0;
    ResetColourAttrib();

    ForEachAtom
        if( ptr->flag&SelectFlag )
        {   if( IsAminoNucleo(group->refno) )
            {   ref = Shapely + group->refno;
            } else ref = Shapely+30;

/*  Original Colour Scheme
 *
 *  ref = &(Shapely[26]);
 *  if( IsNucleo(group->refno) )
 *  {   ref = Shapely + group->refno;
 *  } else if( IsShapelyBackbone(ptr->refno) )
 *  {   ref = &(Shapely[24]);
 *  } else if( IsShapelySpecial(ptr->refno) )
 *  {   ref = &(Shapely[25]);
 *  } else if( IsAmino(group->refno) )
 *      ref = Shapely + group->refno;
 */

            if( !ref->col )
            {   ref->shade = DefineShade( ref->r, ref->g, ref->b );
                ref->col = Shade2Colour(ref->shade);
            }
            Shade[ref->shade].refcount++;
            ptr->col = ref->col;
        }
}


void StructColourAttrib()
{
    register ShadeRef *ref;
    register Chain __far *chain;
    register Group __far *group;
    register Atom __far *ptr;
    register int i;

    if( !Database )
        return;

    if( InfoHelixCount<0 )
        DetermineStructure();

    for( i=0; i<30; i++ )
        StructShade[i].col = 0;
    ResetColourAttrib();

    ForEachAtom
        if( ptr->flag&SelectFlag )
        {   if( group->flag & HelixFlag )
            {   ref = StructShade+1;
            } else if( group->flag & SheetFlag )
            {   ref = StructShade+2;
            } else if( group->flag & TurnFlag )
            {   ref = StructShade+3;
            } else ref = StructShade;

            if( !ref->col )
            {   ref->shade = DefineShade( ref->r, ref->g, ref->b );
                ref->col = Shade2Colour(ref->shade);
            }
            Shade[ref->shade].refcount++;
            ptr->col = ref->col;
        }
}


void InitialTransform()
{
    register Card dist;
    register Long x, y, z;
    register Long dx, dy, dz;
    register Card ax, ay, az;
    register Chain __far *chain;
    register Group __far *group;
    register Atom __far *ptr;


    dx = MaxX-MinX;   OrigCX = (dx>>1)+MinX;
    dy = MaxY-MinY;   OrigCY = (dy>>1)+MinY;
    dz = MaxZ-MinZ;   OrigCZ = (dz>>1)+MinZ;

    MaxX -= OrigCX;   MinX -= OrigCX;
    MaxY -= OrigCY;   MinY -= OrigCY;
    MaxZ -= OrigCZ;   MinZ -= OrigCZ;

    WorldRadius = 0;
    SideLen = MaxFun(dx,dy);
    if( dz>SideLen ) SideLen = dz;
    SideLen += 1000;  Offset = SideLen>>1;
    XOffset = WRange;  YOffset = HRange;

    ForEachAtom
    {   x = ptr->xorg-OrigCX;   ptr->xorg = x;  ax = (Card)AbsFun(x);
        y = ptr->yorg-OrigCY;   ptr->yorg = y;  ay = (Card)AbsFun(y);
        z = ptr->zorg-OrigCZ;   ptr->zorg = z;  az = (Card)AbsFun(z);
        dist = ax*ax + ay*ay + az*az;
        if( dist>WorldRadius )
	    WorldRadius = dist;
    }

    WorldRadius = (Card)sqrt((double)WorldRadius);
    WorldSize = WorldRadius<<1;
    DScale = 1.0/(WorldSize+1000);

    /* MaxZoom*DScale*Range*500 == 118 */
    MaxZoom = 0.236*(WorldSize+1000)/Range;
    ZoomRange = Range;
    MaxZoom -= 1.0;

}


void ReviseInvMatrix()
{
    InvX[0] = IScale*( DirX[0] = RotY[1]*RotZ[2]-RotY[2]*RotZ[1] );
    InvX[1] = IScale*( DirX[1] = RotX[2]*RotZ[1]-RotX[1]*RotZ[2] );
    InvX[2] = IScale*( DirX[2] = RotX[1]*RotY[2]-RotX[2]*RotY[1] );

    InvY[0] = IScale*( DirY[0] = RotY[2]*RotZ[0]-RotY[0]*RotZ[2] );
    InvY[1] = IScale*( DirY[1] = RotX[0]*RotZ[2]-RotX[2]*RotZ[0] );
    InvY[2] = IScale*( DirY[2] = RotX[2]*RotY[0]-RotX[0]*RotY[2] );

    InvZ[0] = IScale*( DirZ[0] = RotY[0]*RotZ[1]-RotY[1]*RotZ[0] );
    InvZ[1] = IScale*( DirZ[1] = RotX[1]*RotZ[0]-RotX[0]*RotZ[1] );
    InvZ[2] = IScale*( DirZ[2] = RotX[0]*RotY[1]-RotX[1]*RotY[0] );

    ShadowTransform();
}


void PrepareTransform()
{
    register Real theta, temp;
    register Real cost, sint;
    register Real x, y, z;
    register Real ncost;

    if( (ReDrawFlag&RFRotateX) && (DialValue[0]!=LastRX) )
    {   theta = PI*(DialValue[0]-LastRX);
        cost = cos(theta);  sint = sin(theta);
        LastRX = DialValue[0];

        y=RotY[0]; z=RotZ[0];
        RotY[0]=cost*y+sint*z; 
        RotZ[0]=cost*z-sint*y;

        y=RotY[1]; z=RotZ[1];
        RotY[1]=cost*y+sint*z;
        RotZ[1]=cost*z-sint*y;

        y=RotY[2]; z=RotZ[2];
        RotY[2]=cost*y+sint*z;
        RotZ[2]=cost*z-sint*y;

        if( UseShadow )
        {   y=DirX[1]; z=DirX[2];
            DirX[1]=cost*y+sint*z;
            DirX[2]=cost*z-sint*y;

            y=DirY[1]; z=DirY[2];
            DirY[1]=cost*y+sint*z;
            DirY[2]=cost*z-sint*y;

            y=DirZ[1]; z=DirZ[2];
            DirZ[1]=cost*y+sint*z;
            DirZ[2]=cost*z-sint*y;
        }

        if( CenX || CenY || CenZ )
        {   ncost = 1.0-cost;
            temp =  CenX*(ncost*RotY[0] + sint*RotZ[0]);
            temp += CenY*(ncost*RotY[1] + sint*RotZ[1]);
            temp += CenZ*(ncost*RotY[2] + sint*RotZ[2]);
            temp = DialValue[5] - (Scale*temp)/YRange;

            if( temp < -1.0 )
            {   DialValue[5] = -1.0;
            } else if( temp > 1.0 )
            {   DialValue[5] = 1.0;
            } else DialValue[5] = temp;
        }
    }

    if( (ReDrawFlag&RFRotateY) && (DialValue[1]!=LastRY) )
    {   theta = PI*(DialValue[1]-LastRY);
        cost = cos(theta);  sint = sin(theta);
        LastRY = DialValue[1];

        x=RotX[0]; z=RotZ[0];
        RotX[0]=cost*x+sint*z;
        RotZ[0]=cost*z-sint*x;

        x=RotX[1]; z=RotZ[1];
        RotX[1]=cost*x+sint*z;
        RotZ[1]=cost*z-sint*x;

        x=RotX[2]; z=RotZ[2];
        RotX[2]=cost*x+sint*z;
        RotZ[2]=cost*z-sint*x;

        if( UseShadow )
        {   x=DirX[0]; z=DirX[2];
            DirX[0]=cost*x+sint*z;
            DirX[2]=cost*z-sint*x;

            x=DirY[0]; z=DirY[2];
            DirY[0]=cost*x+sint*z;
            DirY[2]=cost*z-sint*x;

            x=DirZ[0]; z=DirZ[2];
            DirZ[0]=cost*x+sint*z;
            DirZ[2]=cost*z-sint*x;
        }

        if( CenX || CenY || CenZ )
        {   ncost = 1.0-cost;
            temp =  CenX*(ncost*RotX[0] + sint*RotZ[0]);
            temp += CenY*(ncost*RotX[1] + sint*RotZ[1]);
            temp += CenZ*(ncost*RotX[2] + sint*RotZ[2]);
            temp = DialValue[4] - (Scale*temp)/XRange;

            if( temp < -1.0 )
            {   DialValue[4] = -1.0;
            } else if( temp > 1.0 )
            {   DialValue[4] = 1.0;
            } else DialValue[4] = temp;
        }
    }

    if( (ReDrawFlag&RFRotateZ) && (DialValue[2]!=LastRZ) )
    {   theta = PI*(DialValue[2]-LastRZ);
        cost = cos(theta);  sint = sin(theta);
        LastRZ = DialValue[2];

        x=RotX[0]; y=RotY[0];
        RotX[0]=cost*x-sint*y;
        RotY[0]=cost*y+sint*x;

        x=RotX[1]; y=RotY[1];
        RotX[1]=cost*x-sint*y;
        RotY[1]=cost*y+sint*x;

        x=RotX[2]; y=RotY[2];
        RotX[2]=cost*x-sint*y;
        RotY[2]=cost*y+sint*x;

        if( UseShadow )
        {   x=DirX[0]; y=DirX[1];
            DirX[0]=cost*x-sint*y;
            DirX[1]=cost*y+sint*x;

            x=DirY[0]; y=DirY[1];
            DirY[0]=cost*x-sint*y;
            DirY[1]=cost*y+sint*x;

            x=DirZ[0]; y=DirZ[1];
            DirZ[0]=cost*x-sint*y;
            DirZ[1]=cost*y+sint*x;
        }

        if( CenX || CenY || CenZ )
        {   ncost = 1.0-cost;
            temp =  CenX*(ncost*RotX[0] - sint*RotY[0]);
            temp += CenY*(ncost*RotX[1] - sint*RotY[1]);
            temp += CenZ*(ncost*RotX[2] - sint*RotY[2]);
            temp = DialValue[4] - (Scale*temp)/XRange;

            if( temp < -1.0 )
            {   DialValue[4] = -1.0;
            } else if( temp > 1.0 )
            {   DialValue[4] = 1.0;
            } else DialValue[4] = temp;

            ncost = 1.0-cost;
            temp =  CenX*(ncost*RotY[0] + sint*RotX[0]);
            temp += CenY*(ncost*RotY[1] + sint*RotX[1]);
            temp += CenZ*(ncost*RotY[2] + sint*RotX[2]);
            temp = DialValue[5] - (Scale*temp)/YRange;

            if( temp < -1.0 )
            {   DialValue[5] = -1.0;
            } else if( temp > 1.0 )
            {   DialValue[5] = 1.0;
            } else DialValue[5] = temp;
        }
    }
}


void ApplyTransform()
{
    register int temp;
    register Real x, y, z;
    register int oldx,oldy;
    register Chain __far *chain;
    register Group __far *group;
    register HBond __far *hptr;
    register Bond __far *bptr;
    register Atom __far *ptr;


    if( ReDrawFlag & RFMagnify )
    {   if( DialValue[3] <= 0.0 )
        {   Zoom = DialValue[3]+1.0;
            if( Zoom<0.1 ) Zoom=0.1;
        } else Zoom = (DialValue[3]*MaxZoom) + 1.0;

        Scale = Zoom*DScale*Range;
        ImageSize = (int)(Scale*WorldSize);
        ImageRadius = ImageSize>>1;
        IScale = 1.0/Scale;

        MaxAtomRadius = 0;
        MaxBondRadius = 0;
    }

    if( ReDrawFlag & RFRotate )
    {   PrepareTransform();
        if( UseShadow )
            ShadowTransform();
    }

    if( ReDrawFlag & (RFRotate|RFMagnify) )
    {   MatX[0] = Scale*RotX[0]; 
        MatX[1] = Scale*RotX[1];
        MatX[2] = Scale*RotX[2];

        MatY[0] = Scale*RotY[0];
        MatY[1] = Scale*RotY[1];
        MatY[2] = Scale*RotY[2];

        MatZ[0] = Scale*RotZ[0];
        MatZ[1] = Scale*RotZ[1];
        MatZ[2] = Scale*RotZ[2];

        if( UseShadow )
        {   InvX[0] = IScale*DirX[0]; 
            InvX[1] = IScale*DirX[1];
            InvX[2] = IScale*DirX[2];

            InvY[0] = IScale*DirY[0];
            InvY[1] = IScale*DirY[1];
            InvY[2] = IScale*DirY[2];

            InvZ[0] = IScale*DirZ[0];
            InvZ[1] = IScale*DirZ[1];
            InvZ[2] = IScale*DirZ[2];
        }
    }

    oldx = XOffset;
    oldy = YOffset;
    XOffset = WRange + (int)(DialValue[4]*XRange);
    YOffset = HRange + (int)(DialValue[5]*YRange);

    /* Zoom dependent Translation! */
    /* XOffset = WRange + (int)(DialValue[4]*ImageSize); */
    /* YOffset = HRange + (int)(DialValue[5]*ImageSize); */


    switch( ReDrawFlag )
    {   case(RFTransX):
                if( temp = XOffset-oldx ) 
                    ForEachAtom ptr->x += temp;
                break;

        case(RFTransY):
                if( temp = YOffset-oldy ) 
                    ForEachAtom ptr->y += temp;
                break;

        case(RFRotateX):
            ForEachAtom
            {   x = ptr->xorg; y = ptr->yorg; z = ptr->zorg;
                ptr->z = (int)(x*MatZ[0]+y*MatZ[1]+z*MatZ[2])+ImageRadius;
                ptr->y = (int)(x*MatY[0]+y*MatY[1]+z*MatY[2])+YOffset;
            }
            break;

        case(RFRotateY):
            ForEachAtom
            {   x = ptr->xorg; y = ptr->yorg; z = ptr->zorg;
                ptr->z = (int)(x*MatZ[0]+y*MatZ[1]+z*MatZ[2])+ImageRadius;
                ptr->x = (int)(x*MatX[0]+y*MatX[1]+z*MatX[2])+XOffset;
            }
            break;

        case(RFRotateZ):
            ForEachAtom
            {   x = ptr->xorg; y = ptr->yorg; z = ptr->zorg;
                ptr->x = (int)(x*MatX[0]+y*MatX[1]+z*MatX[2])+XOffset;
                ptr->y = (int)(x*MatY[0]+y*MatY[1]+z*MatY[2])+YOffset;
            }
            break;

        default:
            if( DrawAtoms && (ReDrawFlag&RFMagnify) )
            {   ForEachAtom 
                {   x = ptr->xorg; y = ptr->yorg; z = ptr->zorg;
                    ptr->z = (int)(x*MatZ[0]+y*MatZ[1]+z*MatZ[2])+ImageRadius;
                    ptr->x = (int)(x*MatX[0]+y*MatX[1]+z*MatX[2])+XOffset;
                    ptr->y = (int)(x*MatY[0]+y*MatY[1]+z*MatY[2])+YOffset;
                    if( ptr->flag&SphereFlag )
                    {   ptr->irad = (int)(Scale*ptr->radius);
                        if( ptr->irad>MaxAtomRadius )
                            MaxAtomRadius = ptr->irad;
                    }
                }
            } else
                ForEachAtom 
                {   x = ptr->xorg; y = ptr->yorg; z = ptr->zorg;
                    ptr->z = (int)(x*MatZ[0]+y*MatZ[1]+z*MatZ[2])+ImageRadius;
                    ptr->x = (int)(x*MatX[0]+y*MatX[1]+z*MatX[2])+XOffset;
                    ptr->y = (int)(x*MatY[0]+y*MatY[1]+z*MatY[2])+YOffset;
                }

            if( ReDrawFlag & RFMagnify )
            {   if( DrawBonds )
                    ForEachBond
                        if( bptr->flag&CylinderFlag )
                        {   bptr->irad = (int)(Scale*bptr->radius);
                            if( bptr->irad>MaxBondRadius )
                            MaxBondRadius = bptr->irad;
                        }

                for( hptr=Database->hlist; hptr; hptr=hptr->hnext )
                    if( hptr->flag&CylinderFlag )
                        hptr->irad = (int)(Scale*hptr->radius);

                for( hptr=Database->slist; hptr; hptr=hptr->hnext )
                    if( hptr->flag&CylinderFlag )
                        hptr->irad = (int)(Scale*hptr->radius);

                ForEachBack
                    if( bptr->flag&CylinderFlag )
                        bptr->irad = (int)(Scale*bptr->radius);
            }
    }

    DetermineClipping();
    if( UseScreenClip || ReDrawFlag!=RFRotateY )
        BucketFlag = False;
}


void ResetTransform()
{
    RotX[0] = 1.0;  RotX[1] = 0.0;  RotX[2] = 0.0;
    RotY[0] = 0.0;  RotY[1] = 1.0;  RotY[2] = 0.0;
    RotZ[0] = 0.0;  RotZ[1] = 0.0;  RotZ[2] = 1.0;

    DirX[0] = 1.0;  DirX[1] = 0.0;  DirX[2] = 0.0;
    DirY[0] = 0.0;  DirY[1] = 1.0;  DirY[2] = 0.0;
    DirZ[0] = 0.0;  DirZ[1] = 0.0;  DirZ[2] = 1.0;

    LastRX = LastRY = LastRZ = 0.0;
    CenX = CenY = CenZ = 0;
}


void InitialiseTransform()
{
    ColourDepth = DefaultColDepth;
    ColourMask = ColourDepth-1;
    Ambient = DefaultAmbient;

    LastShade = Colour2Shade(LutSize);

    ResetColourMap();
    ResetTransform();

    ZoneBoth = False;
    HetaGroups = True;
    Hydrogens = True;
}

