/* outfile.c
 * RasMol2 Molecular Graphics
 * Roger Sayle, February 1994
 * Version 2.3
 */
#define OUTFILE
#include "rasmol.h"

#ifdef IBMPC
#include <windows.h>
#include <malloc.h>
#endif
#ifndef sun386
#include <stdlib.h>
#endif

#include <stdio.h>
#include <ctype.h>
#include <math.h>

#include "outfile.h"
#include "molecule.h"
#include "command.h"
#include "transfor.h"
#include "render.h"
#include "graphics.h"
#include "pixutils.h"


#ifdef EIGHTBIT
#define RComp(x)   (RLut[LutInv[x]])
#define GComp(x)   (GLut[LutInv[x]])
#define BComp(x)   (BLut[LutInv[x]])
#else
#define RComp(x)   (((x)>>16)&0xff)
#define GComp(x)   (((x)>>8)&0xff)
#define BComp(x)   ((x)&0xff)
#endif


/* Sun rasterfile.h macro defns */
#define RAS_MAGIC       0x59a66a95
#define RAS_RLE         0x80
#define RT_STANDARD     1
#define RT_BYTE_ENCODED 2
#define RMT_NONE        0
#define RMT_EQUAL_RGB   1

/* Standard A4 size page: 8.267x11.811 inches */
/* U.S. Normal size page: 8.500x11.000 inches */
#define PAGEHIGH  (11.811*72.0)
#define PAGEWIDE  (8.267*72.0)
#define BORDER    0.90

/* Compression Ratio   0<x<127 */
#define EPSFCompRatio  32

#define Round(x)  ((int)(x))
#define Rad2Deg   (180.0/PI)
#define Deg2Rad   (PI/180.0)


#define PSLine    0x00
#define PSStick   0x01
#define PSSphere  0x02

typedef void __far* PSItemPtr;


#ifdef IBMPC
static short __far *ABranch;
static short __far *DBranch;
static short __far *Hash;
static Byte __far *Node;
#else
static short ABranch[4096];
static short DBranch[4096];
static short Hash[256];
static Byte Node[4096];
#endif


typedef struct {
	Byte len;
	Byte ch;
        } BMPPacket;
        

static BMPPacket BMPBuffer[10];
static int BMPCount,BMPTotal;        
static int BMPPad;

static int GIFClrCode; 
static int GIFEOFCode;

static Card RLEFileSize;
static int RLEEncode;
static int RLEOutput;
static int RLELength;
static int RLEPixel;
static int RLEChar;


static Byte Buffer[256];
static Byte LutInv[256];
static int LineLength;
static FILE *OutFile;
static Card BitBuffer;
static int BitBufLen;
static int PacketLen;
static int CodeSize;

static Real LineWidth;
static int VectSolid;
static int VectCol;

/* Macros for commonly used loops */
#define ForEachAtom  for(chain=Database->clist;chain;chain=chain->cnext) \
                     for(group=chain->glist;group;group=group->gnext)    \
                     for(aptr=group->alist;aptr;aptr=aptr->anext)
#define ForEachBond  for(bptr=Database->blist;bptr;bptr=bptr->bnext)
#define ForEachBack  for(chain=Database->clist;chain;chain=chain->cnext) \
                     for(bptr=chain->blist;bptr;bptr=bptr->bnext)



static void FatalOutputError( ptr )
    char *ptr;
{
    if( CommandActive ) WriteChar('\n');
    WriteString("Output Error: Unable to create file `");
    WriteString( ptr );  WriteString("'!\n");
    CommandActive = False;
}



static void WritePPMWord( i )
    int i;
{
    if( i>99 )
    {   fputc((i/100)+'0',OutFile); i %= 100;
        fputc((i/10) +'0',OutFile); i %= 10;
    } else if( i>9 )
    {   fputc((i/10)+'0',OutFile);  i %= 10;
    }
    putc(i+'0',OutFile);
}


int WritePPMFile( name, raw )
    char *name;  int raw;
{
    register Pixel __huge *ptr;
    register int i, c;
    register int x,y;

#ifdef IBMPC
    OutFile = fopen(name, (raw?"wb":"w") );
#else
    OutFile = fopen(name,"w");
#endif

    if( !OutFile ) 
    {   FatalOutputError(name);
        return( False );
    }
    fputc('P',OutFile); fputc((raw?'6':'3'),OutFile);
    fprintf(OutFile," %d %d 255\n",XRange,YRange);

#ifdef EIGHTBIT
    for( i=0; i<256; i++ )
        if( ULut[i] )
            LutInv[Lut[i]] = i;
#endif

#ifdef IBMPC
    FBuffer = (Pixel __huge*)GlobalLock(FBufHandle);
#endif

#ifndef INVERT
    ptr = FBuffer;
#endif

    c = 0;

    if( !raw )
    {   c = 0;
        for( y=YRange-1; y>=0; y-- )
        {
#ifdef INVERT
            ptr = FBuffer + (Long)y*XRange;
#endif
            for( x=0; x<XRange; x++ )
            {   i = *ptr++; c++;
                WritePPMWord((int)RComp(i));  fputc(' ',OutFile);
                WritePPMWord((int)GComp(i));  fputc(' ',OutFile);
                WritePPMWord((int)BComp(i));  
                if( c==5 )
                { c=0; fputc('\n',OutFile);
                } else fputc(' ',OutFile);
            }
        }
    } else
        for( y=YRange-1; y>=0; y-- )
        {
#ifdef INVERT
            ptr = FBuffer + (Long)y*XRange;
#endif
            for( x=0; x<XRange; x++ )
            {   i = *ptr++;
                fputc((int)RComp(i),OutFile);
                fputc((int)GComp(i),OutFile);
                fputc((int)BComp(i),OutFile);
            }
        }

    fclose(OutFile);
#ifdef IBMPC
    GlobalUnlock(FBufHandle); 
#endif
    return( True );
}



#ifdef EIGHTBIT
static int CompactColourMap()
{
    register Pixel __huge *ptr;
    register Long pos, count;
    register int i, cols;

    for( i=0; i<256; i++ )
    {   if( ULut[i] )
            LutInv[Lut[i]] = i;
        Buffer[i] = 0;
        Node[i] = 5;
    }

#ifdef IBMPC
    ptr = (Pixel __huge*)GlobalLock(FBufHandle);    
#else
    ptr = FBuffer;
#endif

    cols = 0;
    count = (Long)XRange*YRange;
    for( pos=0; pos<count; pos++ )
    {   i = LutInv[*ptr++];
        if( !Buffer[i] ) 
        {   Node[cols++] = i;
            Buffer[i] = cols;
        }
    }

    for( i=0; i<256; i++ )
        LutInv[i] = Buffer[LutInv[i]]-1;
#ifdef IBMPC
    GlobalUnlock(FBufHandle);
#endif
    return( cols );
}
#endif



static void WriteGIFCode( code )
    int code;
{
    register int max;

    max = (code==GIFEOFCode)? 0 : 7;
    BitBuffer |= ((Card)code<<BitBufLen);
    BitBufLen += CodeSize;

    while( BitBufLen > max )
    {    
#ifdef IBMPC
         Buffer[PacketLen++]=(Byte)(BitBuffer&0xff);
#else
         Buffer[PacketLen++]=BitBuffer;
#endif
         BitBuffer >>= 8;
         BitBufLen -= 8;

        if( PacketLen==255 )
        {   fputc(0xff,OutFile);
            fwrite((char*)Buffer,1,255,OutFile);
            PacketLen = 0;
        }
    }
}

int WriteGIFFile( name )
    char *name;
{
    register int i,j,cols;
    register int pref,next,last;
    register int isize, ilast;
    register Pixel __huge *ptr;
    register short __far *prev;
    register int x,y,init;

#ifdef EIGHTBIT
    cols = CompactColourMap();
    if( cols<2 ) return( False );

    for( isize=0; isize<8; isize++ )
        if( (1<<isize)>=cols ) break;
    cols = 1<<isize;

#ifdef IBMPC
    OutFile = fopen(name,"wb");
#else
    OutFile = fopen(name,"w");
#endif
    if( !OutFile ) 
    {    FatalOutputError(name);
         return( False );
    }
    fwrite("GIF87a",1,6,OutFile);
    fputc(XRange&0xff,OutFile);  fputc((XRange>>8)&0xff,OutFile);
    fputc(YRange&0xff,OutFile);  fputc((YRange>>8)&0xff,OutFile);
    fputc(0xf0|(isize-1),OutFile); 
    fputc(0x00,OutFile); 
    fputc(0x00,OutFile);

    for( j=0; j<cols; j++ )
    {   i = Node[j];
        fputc((int)RLut[i],OutFile);
        fputc((int)GLut[i],OutFile);
        fputc((int)BLut[i],OutFile);
    }

    fputc(',',OutFile);
    fputc(0x00,OutFile);  fputc(0x00,OutFile);
    fputc(0x00,OutFile);  fputc(0x00,OutFile);
    fputc(XRange&0xff,OutFile);  fputc((XRange>>8)&0xff,OutFile);
    fputc(YRange&0xff,OutFile);  fputc((YRange>>8)&0xff,OutFile);
    fputc(0x00,OutFile);  fputc(isize,OutFile);

    PacketLen=0;
    BitBuffer=0;
    BitBufLen=0;

    GIFClrCode = (1<<isize);
    GIFEOFCode = GIFClrCode+1;
    ilast = (GIFClrCode<<1)-GIFEOFCode;
    isize++;

    CodeSize = isize;
    last = ilast;
    next = 1;  
   
    WriteGIFCode(GIFClrCode);
    for( i=0; i<cols; i++ )
        Hash[i]=0;

#ifdef IBMPC
    FBuffer = (Pixel __huge*)GlobalLock(FBufHandle);    
#endif

#ifndef INVERT
    ptr = FBuffer;
#endif

    /* Avoid Warnings! */
    prev = (short __far*)0; 
    pref = 0;

    init = False;
    for( y=YRange-1; y>=0; y-- )
    {   
#ifdef INVERT
        ptr = FBuffer + (Long)y*XRange;
#endif
        for( x=0; x<XRange; x++ )
        {   if( !init )
            {   pref = LutInv[*ptr++];
                prev = Hash+pref;
                init = True;
                continue;
            }

            i = LutInv[*ptr++];

            while( *prev && (Node[*prev] != (Byte)i) )
                prev = ABranch+*prev;

            if( *prev )
            {   pref = *prev+GIFEOFCode;
                prev = DBranch+*prev;
            } else
            {   WriteGIFCode(pref);
                if( next==last )
                {   if( CodeSize==12 )
                    {   WriteGIFCode(GIFClrCode);
                        pref = i;  prev = Hash+i;
                        for( i=0; i<cols; i++ )
                            Hash[i] = 0;
                        CodeSize = isize;
                        last = ilast;
                        next = 1; 
                        continue;
                    }
                    last = (last<<1)+GIFEOFCode;
                    CodeSize++;
                }
                *prev = next;
                ABranch[next] = 0;
                DBranch[next] = 0;
                Node[next] = i;
                prev = Hash+i;
                pref = i;
                next++;
            }
        }
    }


    WriteGIFCode(pref);
    WriteGIFCode(GIFEOFCode);
    if( PacketLen )
    {   fputc(PacketLen,OutFile);
        fwrite((char*)Buffer,1,PacketLen,OutFile);
    }

    fputc(0x00,OutFile);
    fputc(';',OutFile);
    fclose(OutFile);

#ifdef IBMPC
    GlobalUnlock(FBufHandle);
#endif
    return( True );
#else
    if( CommandActive ) WriteChar('\n');
    WriteString("Output Error: 24 bit GIF files unsupported!\n");
    CommandActive = False;
    return( False );
#endif
}

/* Function Prototype */
#if defined(__STDC__) || defined(IBMPC)
static void WriteRastInt( Card );
#endif

static void WriteRastInt( val )
    Card val;
{
    fputc((int)((val>>24)&0xff),OutFile);
    fputc((int)((val>>16)&0xff),OutFile);
    fputc((int)((val>>8) &0xff),OutFile);
    fputc((int)(val&0xff),OutFile);
}


static void FlushRastRLE()
{
    if( RLEChar==RAS_RLE )
    {   if( RLEOutput )
        {   fputc(RAS_RLE,OutFile);
            fputc(RLELength-1,OutFile);
            if( RLELength!=1 )
                fputc(RAS_RLE,OutFile);
        } else RLEFileSize += (RLELength>1)? 3 : 2;
    } else
        if( RLEOutput )
        {   if( RLELength>2 )
            {   fputc(RAS_RLE,OutFile);
                fputc(RLELength-1,OutFile);
            } else if( RLELength==2 )
                fputc(RLEChar,OutFile);
            fputc(RLEChar,OutFile);
        } else RLEFileSize += MinFun(RLELength,3);
}


static void WriteRastRLECode( val )
    int val;
{
    if( RLEEncode )
    {   if( !RLELength )
        {   RLELength = 1;
            RLEChar = val;
        } else
            if( (RLEChar!=val) || (RLELength==256) )
            {   FlushRastRLE();
                RLELength = 1;
                RLEChar = val;
            } else
                RLELength++;
    } else
        fputc(val,OutFile);
}

static void WriteRastRLEPad()
{
    if( RLEEncode )
    {   if( !RLELength || (RLELength==256) )
        {   WriteRastRLECode(0x00);
        } else RLELength++;
    } else fputc(0x00,OutFile);
}


static void WriteRastData( output )
    int output;
{
    register Pixel __huge *ptr;
    register int x,y,pad;
#ifndef EIGHTBIT
    register int i;
#endif

#ifdef IBMPC
    FBuffer = (Pixel __huge*)GlobalLock(FBufHandle);
#endif

#ifndef INVERT
    ptr = FBuffer;
#endif

    pad = XRange%2;

    RLEOutput = output;
    RLEFileSize = 0;
    RLELength = 0;

    for( y=YRange-1; y>=0; y-- )
    {   
#ifdef INVERT
        ptr = FBuffer + (Long)y*XRange;
#endif
        for( x=0; x<XRange; x++ )
#ifndef EIGHTBIT
        {   i = *ptr++;
            WriteRastRLECode((int)BComp(i));
            WriteRastRLECode((int)GComp(i));
            WriteRastRLECode((int)RComp(i));
        }
#else
            WriteRastRLECode((int)LutInv[*ptr++]);
#endif
        if( pad ) WriteRastRLEPad();
    }

    if( RLEEncode && RLELength )
        FlushRastRLE();
#ifdef IBMPC
    GlobalUnlock(FBufHandle);
#endif
}



int WriteRastFile( name, encode )
    char *name;
    int encode;
{
    register int i,size,cols;

#ifdef IBMPC
    OutFile = fopen(name,"wb");
#else
    OutFile = fopen(name,"w");
#endif

    if( !OutFile )
    {   FatalOutputError(name);
        return(False);
    }
    WriteRastInt( RAS_MAGIC );
    WriteRastInt(XRange);  
    WriteRastInt(YRange);
    RLEEncode = encode;


#ifdef EIGHTBIT
    WriteRastInt(8);

    if( encode )
    {   WriteRastData(False);
        WriteRastInt(RLEFileSize);
        WriteRastInt(RT_BYTE_ENCODED);
    } else
    {   size = (XRange%2)? XRange+1 : XRange;
        WriteRastInt(size*YRange);
        WriteRastInt(RT_STANDARD);
    }

    cols = CompactColourMap();
    WriteRastInt(RMT_EQUAL_RGB);
    WriteRastInt(cols*3);

    for( i=0; i<cols; i++ ) fputc((int)RLut[Node[i]],OutFile);
    for( i=0; i<cols; i++ ) fputc((int)GLut[Node[i]],OutFile);
    for( i=0; i<cols; i++ ) fputc((int)BLut[Node[i]],OutFile);
#else
    WriteRastInt(24);

    if( encode )
    {   WriteRastData(False);
        WriteRastInt(RLEFileSize);
        WriteRastInt(RT_BYTE_ENCODED);
    } else
    {   size = XRange*3;
        if( size&1 ) size++;
        WriteRastInt(size*YRange);
        WriteRastInt(RT_STANDARD);
    }
    WriteRastInt(RMT_NONE);
    WriteRastInt(0);
#endif

    WriteRastData(True);
    fclose( OutFile );
    return( True );
}


static void OutputEPSFByte( val )
    int val;
{
    register int i;

    i = val/16;  fputc( (i>9)? (i-10)+'A' : i+'0', OutFile );
    i = val%16;  fputc( (i>9)? (i-10)+'A' : i+'0', OutFile );
    if( (LineLength+=2) > 72 )
    {   fputc('\n',OutFile);
        LineLength = 0;
    }
}

static void EncodeEPSFPixel( val, col )
    int val, col;
{
    register int r, g, b;
    register int i;

    r = RComp(val);
    g = GComp(val);
    b = BComp(val);

    if( col )
    {   OutputEPSFByte( r );
        OutputEPSFByte( g );
        OutputEPSFByte( b );
    } else
    {   i = (20*r + 32*g + 12*b)>>6;
        OutputEPSFByte( i );
    }
}

int WriteEPSFFile( name, col, comp )
    char *name;  int col, comp;
{
    register int xpos, ypos;
    register int xsize, ysize;
    register int rotpage;
    register int x, y, i, j;

    register Pixel __huge *ptr;
    int RLEBuffer[128];

#ifdef EIGHTBIT
    for( i=0; i<256; i++ )
        if( ULut[i] )
            LutInv[Lut[i]] = i;
#endif

    OutFile = fopen(name,"w");
    if( !OutFile )
    {   FatalOutputError(name);
        return(False);
    }

    if( XRange <= YRange )
    {   rotpage = False; 
        xsize = XRange; 
        ysize = YRange;
    } else
    {   rotpage = True;  
        xsize = YRange; 
        ysize = XRange;
    }

    if( xsize > (int)(BORDER*PAGEWIDE) )
    {   ysize = (int)(ysize*(BORDER*PAGEWIDE)/xsize);
        xsize = (int)(BORDER*PAGEWIDE);
    }
    if( ysize > (int)(BORDER*PAGEHIGH) )
    {   xsize = (int)(xsize*(BORDER*PAGEHIGH)/ysize);
        ysize = (int)(BORDER*PAGEHIGH);
    }

    xpos = (int)(PAGEWIDE-xsize)/2;
    ypos = (int)(PAGEHIGH-ysize)/2;

    fputs("%!PS-Adobe-2.0 EPSF-2.0\n",OutFile);
    fputs("%%Creator: RasMol Version 2.3\n",OutFile);
    fprintf(OutFile,"%%%%Title: %s\n",name);
    fprintf(OutFile,"%%%%BoundingBox: %d %d ",xpos,ypos);
    fprintf(OutFile,"%d %d\n",xpos+xsize,ypos+ysize);

    fputs("%%Pages: 1\n",OutFile);
    fputs("%%EndComments\n",OutFile);
    fputs("%%EndProlog\n",OutFile);
    fputs("%%Page: 1 1\n",OutFile);

    fputs("gsave\n",OutFile);
    fputs("10 dict begin\n",OutFile);
    fprintf(OutFile,"%d %d translate\n",xpos,ypos);
    fprintf(OutFile,"%d %d scale\n",xsize,ysize);
    if( rotpage )
    {   fputs("0.5 0.5 translate\n",OutFile);
        fputs("90 rotate\n",OutFile);
        fputs("-0.5 -0.5 translate\n",OutFile);
    }
    fputc('\n',OutFile);

    if( comp )
    {   fputs("/rlecount 0 def\n",OutFile);
        fputs("/rlebyte 1 string def\n",OutFile);
        fprintf(OutFile,"/pixbuf %d string def\n", col?3:1 );
        fputs("/reppixel { pixbuf } def\n",OutFile);
        fputs("/getpixel { \n",OutFile);
        fputs("  currentfile pixbuf readhexstring pop\n",OutFile);
        fputs("} def\n\n",OutFile);

        if( col )
        {   fputs("/colorimage where {\n",OutFile);
            fputs("  pop\n",OutFile);
            fputs("} {\n",OutFile);
            fputs("  /bytebuf 1 string def\n",OutFile);
            fputs("  /colorimage { pop pop image } def\n",OutFile);
            fputs("  /reppixel { bytebuf } def\n",OutFile);
            fputs("  /getpixel {\n",OutFile);
            fputs("    currentfile pixbuf readhexstring pop pop\n",OutFile);
            fputs("    bytebuf 0\n",OutFile);
            fputs("    pixbuf 0 get 20 mul\n",OutFile);
            fputs("    pixbuf 1 get 32 mul\n",OutFile);
            fputs("    pixbuf 2 get 12 mul\n",OutFile);
            fputs("    add add -6 bitshift put bytebuf\n",OutFile);
            fputs("  } def\n",OutFile);
            fputs("} ifelse\n\n",OutFile);
        }

        fputs("/rledecode {\n",OutFile);
        fputs("  rlecount 0 eq {\n",OutFile);
        fputs("    currentfile rlebyte readhexstring pop\n",OutFile);
        fprintf(OutFile,"    0 get dup %d gt {\n",EPSFCompRatio);
        fprintf(OutFile,"      /rlecount exch %d sub def\n",EPSFCompRatio);
        fputs("      /rleflag true def\n",OutFile);
        fputs("    } {\n",OutFile);
        fputs("      /rlecount exch def\n",OutFile);
        fputs("      /rleflag false def\n",OutFile);
        fputs("    } ifelse getpixel\n",OutFile);
        fputs("  } {\n",OutFile);
        fputs("    /rlecount rlecount 1 sub def\n",OutFile);
        fputs("    rleflag { reppixel } { getpixel } ifelse\n",OutFile);
        fputs("  } ifelse\n",OutFile);
        fputs("} def\n",OutFile);
    } else if( col )
    {   fprintf(OutFile,"/scanbuf %d string def\n",XRange*3);
        fputs("/colorimage where {\n",OutFile);
        fputs("  pop\n",OutFile);
        fputs("} {\n",OutFile);
        fputs("  /pixbuf 3 string def\n",OutFile);
        fputs("  /bytebuf 1 string def\n",OutFile);
        fputs("  /colorimage {\n",OutFile);
        fputs("    pop pop pop {\n",OutFile);
        fputs("      currentfile pixbuf readhexstring pop pop\n",OutFile);
        fputs("      bytebuf 0\n",OutFile);
        fputs("      pixbuf 0 get 20 mul\n",OutFile);
        fputs("      pixbuf 1 get 32 mul\n",OutFile);
        fputs("      pixbuf 2 get 12 mul\n",OutFile);
        fputs("      add add -6 bitshift put bytebuf\n",OutFile);
        fputs("    } image\n",OutFile);
        fputs("  } def\n",OutFile);
        fputs("} ifelse\n\n",OutFile);
    } else fprintf(OutFile,"/scanbuf %d string def\n\n",XRange);

    fprintf(OutFile,"%d %d %d\n",XRange,YRange,8);
    fprintf(OutFile,"[%d 0 0 -%d 0 %d]\n",XRange,YRange,YRange);

    if( !comp )
    {   fputs("{ currentfile scanbuf readhexstring pop }\n",OutFile);
    } else fputs("{ rledecode }\n",OutFile);
    if( col ) fputs("false 3 color",OutFile);
    fputs("image\n",OutFile);

#ifdef IBMPC
    FBuffer = (Pixel __huge*)GlobalLock(FBufHandle);
#endif

#ifndef INVERT
    ptr = FBuffer; 
#endif

    RLELength = 0;
    LineLength = 0;
    for( y=YRange-1; y>=0; y-- )
    {
#ifdef INVERT
        ptr = FBuffer + (Long)y*XRange;
#endif
        for( x=0; x<XRange; x++ )
        {   i = *ptr++;
            if( comp )
            {   if( RLELength )
                {   if( RLEEncode )
                    {   if( (RLEPixel!=i) || (RLELength==256-EPSFCompRatio) )
                        {   OutputEPSFByte(RLELength+EPSFCompRatio-1);
                            EncodeEPSFPixel(RLEPixel,col);
                            RLEEncode = False;
                            RLEBuffer[0] = i;
                            RLELength = 1;
                        } else RLELength++;
                    } else if( RLEBuffer[RLELength-1] == i )
                    {   if( RLELength>1 )
                        {   OutputEPSFByte(RLELength-2);
                            for( j=0; j<RLELength-1; j++ )
                                EncodeEPSFPixel(RLEBuffer[j],col);
                        }
                        RLEEncode = True;
                        RLELength = 2;
                        RLEPixel = i;
                    } else if( RLELength == EPSFCompRatio+1 )
                    {   OutputEPSFByte(EPSFCompRatio);
                        for( j=0; j<RLELength; j++ )
                             EncodeEPSFPixel(RLEBuffer[j],col);
                        RLEEncode = False;
                        RLEBuffer[0] = i;
                        RLELength = 1;
                    } else RLEBuffer[RLELength++] = i;
                } else
                {   RLEEncode = False;
                    RLEBuffer[0] = i;
                    RLELength = 1;
                }
            } else EncodeEPSFPixel( i, col );
        }
    }

    if( comp && RLELength )
    {   if( RLEEncode )
        {   OutputEPSFByte(RLELength+EPSFCompRatio-1);
            EncodeEPSFPixel(RLEPixel,col);
        } else
        {   OutputEPSFByte(RLELength-1);
            for( j=0; j<RLELength; j++ )
                EncodeEPSFPixel(RLEBuffer[j],col);
        }
    }

    if( LineLength ) 
        fputc('\n',OutFile);
    fputs("end\n",OutFile);
    fputs("grestore\n",OutFile);
    fputs("showpage\n",OutFile);
    fputs("%%Trailer\n",OutFile);
    fputs("%%EOF\n",OutFile);
    fclose( OutFile );
#ifdef IBMPC
    GlobalUnlock(FBufHandle);
#endif
    return( True );
}



/* Flush AbsMode buffer */
static void FlushBMPBuffer()
{
    if( RLEOutput )
    {   fputc(0x00,OutFile);
        fputc(PacketLen,OutFile);
        fwrite((char*)Buffer,1,PacketLen,OutFile);
        if( PacketLen%2 ) fputc(0x00,OutFile);    
    } else
        RLEFileSize += (PacketLen%2)+PacketLen+2;
     
    PacketLen = 0;
}

/* Flush RLEMode buffer */
static void FlushBMPPackets()
{
    register int i;
    
    if( PacketLen )
        FlushBMPBuffer();

    if( RLEOutput )
    {   for( i=0; i<BMPCount; i++ )
        {   fputc(BMPBuffer[i].len,OutFile);
            fputc(BMPBuffer[i].ch,OutFile);
        }
    } else RLEFileSize += (BMPCount<<1);
    BMPCount = BMPTotal = 0;
}

static void ProcessBMPPacket()
{
    register int cost;
    register int i,j;
    
    BMPBuffer[BMPCount].len = RLELength;
    BMPBuffer[BMPCount].ch = RLEChar;
    BMPTotal += RLELength;
    RLELength = 0;
    BMPCount++;	


    /* RLEMode is more efficient */
    if( BMPTotal > BMPCount+5 )
    {   FlushBMPPackets();
        return;
    }
    
    /* Flush AbsMode buffer */
    if( PacketLen+BMPTotal>255 )
        FlushBMPBuffer();
    

    /* Cannot leave RLEMode */
    if( PacketLen+BMPTotal<3 )
        return;
        
    /* Determine AbsMode cost */
    if( PacketLen )
    {   cost = BMPTotal - (PacketLen%2);
        cost += (cost%2);
    } else cost = (BMPTotal%2)+BMPTotal+2;    

    /* Put RLE Packets into AbsMode buffer */
    if( cost <= (int)(BMPCount<<1) )
    {   for( i=0; i<BMPCount; i++ )
            for( j=0; j<(int)BMPBuffer[i].len; j++ )
                Buffer[PacketLen++] = BMPBuffer[i].ch;
        BMPCount = BMPTotal = 0;
    }
}
 
/* Collect pixels into RLE Packets */
static void WriteBMPCode( val )
    int val;
{
    if( !RLELength )
    {   RLELength = 1;
        RLEChar = val;
    } else
        if( (RLEChar!=val) || (RLELength==255) )
        {   ProcessBMPPacket();
            RLELength = 1;
            RLEChar = val;
        } else
            RLELength++;
}

static void WriteBMPData( output )
    int output;
{
    register Pixel __huge *ptr;
    register int x,y;
    
    RLEOutput = output;
    RLEFileSize = 0;   BMPCount = 0;
    RLELength = 0;     BMPTotal = 0;
    PacketLen = 0; 

#ifdef INVERT
    ptr = FBuffer;
#endif

    for( y=YRange-1; y>=0; y-- )
    {
#ifndef INVERT
        ptr = FBuffer + (Long)y*XRange;
#endif
        for( x=0; x<XRange; x++ )
            WriteBMPCode(*ptr++);

        for( x=0; x<BMPPad; x++ )
            WriteBMPCode(0x00);

        /* Flush RLE codes */
        ProcessBMPPacket();
        FlushBMPPackets();
        
        if( RLEOutput )
        {   /* End of line code */
            fputc(0x00,OutFile);
            fputc((y?0x00:0x01),OutFile);
        } else RLEFileSize += 2;
    }
}


#if defined(__STDC__) || defined(IBMPC)
static void WriteBMPWord( Card );
#endif

static void WriteBMPWord( val )
    Card val;
{
    fputc((int)(val&0xff),OutFile);
    fputc((int)((val>>8) &0xff),OutFile);
    fputc((int)((val>>16)&0xff),OutFile);
    fputc((int)((val>>24)&0xff),OutFile);
}
    
int WriteBMPFile( name )
    char *name;
{
    register Pixel __huge *ptr;
    register int x,y,i,raw;
    register Card size;

#ifdef EIGHTBIT
#ifdef IBMPC
    OutFile = fopen(name,"wb");
#else
    OutFile = fopen(name,"w");
#endif
    if( !OutFile )
    {    FatalOutputError(name);
         return( False );
    }

#ifdef IBMPC
    FBuffer = (Pixel __huge*)GlobalLock(FBufHandle);
#endif

    /* zero-pad scanlines to long */
    if( (BMPPad=XRange%4) ) 
         BMPPad = 4-BMPPad; 
   
    WriteBMPData(False);
    size = (Long)(XRange+BMPPad)*YRange;
    if( RLEFileSize<size )
    {   size = RLEFileSize;
        raw = False;
    } else raw = True;


    fputc('B',OutFile); 
    fputc('M',OutFile);
    WriteBMPWord(size+1078);
    WriteBMPWord((Card)0);
    WriteBMPWord((Card)1078);

    WriteBMPWord((Card)40);
    WriteBMPWord((Card)XRange);
    WriteBMPWord((Card)YRange);
    fputc(0x01,OutFile);  fputc(0x00,OutFile);
    fputc(0x08,OutFile);  fputc(0x00,OutFile);
    WriteBMPWord(raw? (Card)0 : (Card)1);
    WriteBMPWord(size);
    
    WriteBMPWord((Card)0);
    WriteBMPWord((Card)0);
    WriteBMPWord((Card)256);
    WriteBMPWord((Card)256);

    for( i=0; i<256; i++ )
    {   fputc((int)BLut[i],OutFile);
        fputc((int)GLut[i],OutFile);
        fputc((int)RLut[i],OutFile);
        fputc(0x00,OutFile);
    }

    if( raw )
    {   
#ifdef INVERT
        ptr = FBuffer;
#endif
        for( y=YRange-1; y>=0; y-- )
        {
#ifndef INVERT
            ptr = FBuffer + (Long)y*XRange;
#endif
            for( x=0; x<XRange; x++ )
                fputc((int)(*ptr++),OutFile);
            for( x=0; x<BMPPad; x++ )
                fputc(0x00,OutFile);
            fputc(0x00,OutFile);
            fputc((y?0x00:0x01),OutFile);
        }
    } else
        WriteBMPData(True);

#ifdef IBMPC
    GlobalUnlock(FBufHandle);
#endif
    fclose(OutFile);
    return(True);
#else
    if( CommandActive ) WriteChar('\n');
    WriteString("Output Error: 24 bit BMP files unsupported!\n");
    CommandActive = False;
    return( False );
#endif
}


static void MolScriptSegment( ptr, src, dst, chain )
    char *ptr;  int src, dst;  char chain;
{   
    putchar(' ');
    if( (chain!=' ') && !isdigit(chain) ) 
    {   fprintf(OutFile,"%s from %c%d to %c%d;\n",ptr,chain,src,chain,dst);
    } else fprintf(OutFile,"%s from %d to %d\n",ptr,src,dst);
    fputs(";\n",OutFile);
}


int WriteMolScriptFile( name )
    char *name;
{
    register Real psi, phi, theta;
    register Chain __far *chain;
    register Group __far *group;
    register Group __far *next;
    register Group __far *prev;
    register int flag,len;
    register char *ptr;

    if( !Database )
        return(False);

    OutFile = fopen(name,"w");
    if( !OutFile )
    {   FatalOutputError(name);
        return(False);
    }
    fprintf(OutFile,"! File: %s\n",name);
    fputs("! Creator: RasMol Version 2.3\n",OutFile);
    fputs("! Version: MolScript v1.3\n\n",OutFile);

    fputs("plot\n",OutFile);
    if( BackR || BackG || BackB )
    {   fputs("  background rgb ",OutFile);
        fprintf(OutFile,"%#g ",  BackR/255.0);
        fprintf(OutFile,"%#g ",  BackG/255.0);
        fprintf(OutFile,"%#g;\n",BackB/255.0);
    }
    fputc('\n',OutFile);

    fprintf(OutFile,"  read mol \"%s\";\n",InfoFileName);
    fputs("  transform atom *\n",OutFile);
    fputs("    by centre position atom *\n",OutFile);
    fputs("    by rotation x 180.0",OutFile);

    phi = Rad2Deg*asin(RotX[2]);
    if( (int)phi == 90 )
    {   theta = -Rad2Deg*atan2(RotY[0],RotY[1]);
        psi = 0;
    } else if( (int)phi == -90 )
    {   theta = Rad2Deg*atan2(RotY[0],RotY[1]);
        psi = 0;
    } else /* General Case! */
    {   theta = Rad2Deg*atan2(RotY[2],RotZ[2]);
        psi =  -Rad2Deg*atan2(RotX[1],RotX[0]);
    }

#ifdef INVERT
    if( (int)psi )   fprintf(OutFile,"\n    by rotation z %#g",psi);
    if( (int)phi )   fprintf(OutFile,"\n    by rotation y %#g",phi);
    if( (int)theta ) fprintf(OutFile,"\n    by rotation x %#g",-theta);
#else
    if( (int)psi )   fprintf(OutFile,"\n    by rotation z %#g",-psi);
    if( (int)phi )   fprintf(OutFile,"\n    by rotation y %#g",phi);
    if( (int)theta ) fprintf(OutFile,"\n    by rotation x %#g",theta);
#endif
    fputs(";\n\n",OutFile);

    /* fputs("  trace amino-acids;\n",OutFile); */

    if( Database->clist )
    {   if( InfoHelixCount<0 )
            DetermineStructure();

        for( chain=Database->clist; chain; chain=chain->cnext )
        {   prev = chain->glist;
            for( group=prev; group && group->gnext; group=next )
            {   next = group->gnext;
                flag = group->flag & next->flag;

                if( flag&HelixFlag )
                {   flag = HelixFlag;
                    ptr = "helix";
                } else if( flag&SheetFlag )
                {   flag = SheetFlag;
                    ptr = "strand";
                } else if( flag&TurnFlag )
                {   fputs("  turn residue ",OutFile);
                    if( chain->ident != ' ' )
                        fputc(chain->ident,OutFile);
                    fprintf(OutFile,"%d;\n",group->serno);
                    continue;
                } else continue;

                len = 2;  /* Determine Structure Length */
                while( next->gnext && (next->gnext->flag&flag) )
                {   next = next->gnext;
                    len++;
                }

                if( len>2 )
                {   if( prev ) /* MolScript coil or turn? */
                       MolScriptSegment("coil",prev->serno,group->serno,
                                              chain->ident);
                    MolScriptSegment(ptr,group->serno,next->serno,
                                         chain->ident);
	            prev = next;
                } 
            }

            if( prev != group )  /* C-terminal coil/turn */
                MolScriptSegment("coil",prev->serno,group->serno,chain->ident);
        }
    }

    fputs("end_plot\n",OutFile);
    fclose(OutFile);
    return( True );
}



int WriteScriptFile( name )
    char *name;
{
    register int theta,phi,psi;
    register char *ptr;
    register int temp;

    OutFile = fopen(name,"w");
    if( !OutFile )
    {   FatalOutputError(name);
        return(False);
    }

    fprintf(OutFile,"# File: %s\n",name);
    fputs("# Creator: RasMol Version 2.3\n\n",OutFile);
    fprintf(OutFile,"background [%d,%d,%d]\n",BackR,BackG,BackB);
    if( !Database )
    {   /* No Molecule! */
        fclose(OutFile);
        return(True);
    }

    /* Colour Details */
    fprintf(OutFile,"set ambient %d\n", (int)(100*Ambient) );
    fputs("set specular ",OutFile);
    if( FakeSpecular )
    {   fprintf(OutFile,"on\nset specpower %d\n",SpecPower);
    } else fputs("off\n",OutFile);
    putc('\n',OutFile);

    /* Transformation */
    fputs("reset\n",OutFile);
    if( UseSlabPlane )
    {   if( (temp = (int)(50.0*DialValue[7])) )
        {   fprintf(OutFile,"slab %d\n",temp+50);
        } else fputs("slab on\n",OutFile);

        fputs("set slabmode ",OutFile);
        switch( SlabMode )
        {   default:            
            case(SlabClose):    ptr = "solid";    break;
            case(SlabReject):   ptr = "reject";   break;
            case(SlabHalf):     ptr = "half";     break;
            case(SlabHollow):   ptr = "hollow";   break;
            case(SlabSection):  ptr = "section";
        }
        fputs(ptr,OutFile);
        putc('\n',OutFile);
    } else fputs("slab off\n",OutFile);

    phi = Round(Rad2Deg*asin(RotX[2]));
    if( phi == 90 )
    {   theta = -Round(Rad2Deg*atan2(RotY[0],RotY[1]));
        psi = 0;
    } else if( phi == -90 )
    {   theta = Round(Rad2Deg*atan2(RotY[0],RotY[1]));
        psi = 0;
    } else /* General Case! */
    {   theta = Round(Rad2Deg*atan2(RotY[2],RotZ[2]));
        psi =  Round(-Rad2Deg*atan2(RotX[1],RotX[0]));
    }

#ifdef INVERT
    if( psi )   fprintf(OutFile,"rotate z %d\n",-psi);
    if( phi )   fprintf(OutFile,"rotate y %d\n",phi);
    if( theta ) fprintf(OutFile,"rotate x %d\n",-theta);
#else
    if( psi )   fprintf(OutFile,"rotate z %d\n",psi);
    if( phi )   fprintf(OutFile,"rotate y %d\n",phi);
    if( theta ) fprintf(OutFile,"rotate x %d\n",theta);
#endif

    if( (temp = (int)(100.0*DialValue[4])) )
        fprintf(OutFile,"translate x %d\n",temp);
    if( (temp = (int)(100.0*DialValue[5])) )
#ifdef INVERT
        fprintf(OutFile,"translate y %d\n",-temp);
#else
        fprintf(OutFile,"translate y %d\n",temp);
#endif

    if( DialValue[3] != 0.0 )
    {   if( DialValue[3]<0.0 )
        {   temp = (int)(100*DialValue[3]);
        } else temp = (int)(100*MaxZoom*DialValue[3]);
        fprintf(OutFile,"zoom %d\n",temp+100);
    }
    putc('\n',OutFile);

    /* Rendering */
    if( DrawAtoms )
    {   fprintf(OutFile,"set shadows %s\n", UseShadow? "on" : "off" );
    }

    if( UseSlabPlane && (SlabMode==SlabSection) )
    {   /* Section Mode Slabbing! */
        fclose(OutFile);
        return(True);
    }

    if( Database->slist )
    {   fputs("set ssbond ",OutFile);
        fputs(SSBondMode?"backbone":"sidechain",OutFile);
        fputs("\nssbonds on\n",OutFile);
    }

    if( Database->hlist )
    {   fputs("set hbond ",OutFile);
        fputs(HBondMode?"backbone":"sidechain",OutFile);
        fputs("\nhbonds on\n",OutFile);
    }

    if( DrawRibbon )
    {   if( !RibbonMode )
        {   fprintf(OutFile,"set strands %d\n",SplineCount);
            fputs("set ribbons strands\n",OutFile);
        } else fputs("set ribbons solid\n",OutFile);
        
    }

    fprintf(OutFile,"set axes %s\n", DrawAxes? "on":"off" );
    fprintf(OutFile,"set boundingbox %s\n", DrawBoundBox? "on":"off" );
    fprintf(OutFile,"set unitcell %s\n", DrawUnitCell? "on":"off" );

    fclose(OutFile);
    return( True );
}


#if defined(__STDC__) || defined(IBMPC)
static int FindDepth( PSItemPtr, int );
static void DepthSort( PSItemPtr __far*, char __far*, int );
#endif

static int FindDepth( item, type )
     PSItemPtr item;  int type;
{
    register Atom __far *atom;
    register Bond __far *bond;
    register int result;

    if( type==PSSphere )
    {   atom = (Atom __far*)item;
        result = atom->z;
    } else /* PSLine or PSStick */
    {   bond = (Bond __far*)item;
        result = bond->srcatom->z;
        if( result < bond->dstatom->z )
            result = bond->dstatom->z;
    }
    return( result );
}


static void DepthSort( data, type, count )
    PSItemPtr __far *data;
    char __far *type;
    int count;
{
    register char ttmp;
    register void __far *dtmp;
    register int i, j, k;
    register int depth;
    register int temp;

    for( i=1; i<count; i++ )
    {   dtmp = data[i];  
        ttmp = type[i];

        j = i-1;
        depth = FindDepth(dtmp,ttmp);
        temp = FindDepth(data[j],type[j]);
        while( (depth<temp) || ((depth==temp)&&(ttmp<type[j])) )
            if( j-- ) 
            {   temp = FindDepth(data[j],type[j]);
            } else break;
        j++;

        if( j != i )
        {   for( k=i; k>j; k-- )
            {    data[k] = data[k-1];
                 type[k] = type[k-1];
            }
            data[j] = dtmp;
            type[j] = ttmp;
        }
    }
}

#if defined(__STDC__) || defined(IBMPC)
static void WriteVectSphere( Atom __far* );
static void WriteVectStick( Bond __far* );
static void WriteVectWire( Bond __far* );
#endif


static void WriteVectColour( col )
    int col;
{
    if( col != VectCol )
    {   fprintf(OutFile,"%g ",RLut[col]/255.0);
        fprintf(OutFile,"%g ",GLut[col]/255.0);
        fprintf(OutFile,"%g ",BLut[col]/255.0);
        fputs("setrgbcolor\n",OutFile);
        VectCol = col;
    }
}


static void WriteVectSphere( ptr )
    Atom __far *ptr;
{
    register Real temp;
    register int i;

    temp = ptr->z/(2.0*ImageSize)+1.0;
    if( temp != LineWidth )
    {   fprintf(OutFile,"%g setlinewidth\n",temp);
        LineWidth = temp;
    }

    i = ptr->col + ColourMask;
    fprintf(OutFile,"%g ",RLut[i]/255.0);
    fprintf(OutFile,"%g ",GLut[i]/255.0);
    fprintf(OutFile,"%g ",BLut[i]/255.0);
    fprintf(OutFile,"%g ",Scale*ptr->radius);

    fprintf(OutFile,"%d %d ",ptr->x,ptr->y);
    fputs("Sphere\n",OutFile);
}

static void WriteVectWire( ptr )
    Bond __far *ptr;
{
    register Atom __far *src;
    register Atom __far *dst;
    register Real radius;
    register Real temp;
    register Real dist;

    register Real midx, midy;
    register Real endx, endy;
    register int col1, col2;
    register int dx, dy, dz;
    register Long dist2;


    src = ptr->srcatom;  
    dst = ptr->dstatom;
    if( src->z > dst->z )
    {   src = ptr->dstatom;
        dst = ptr->srcatom;
    }

    if( !ptr->col )
    {   col1 = src->col;
        col2 = dst->col;
    } else col1 = col2 = ptr->col;

    dx = dst->x - src->x;  
    dy = dst->y - src->y;
    dist2 = dx*dx + dy*dy;
    dist = sqrt( (double)dist2 );

    if( dst->flag & SphereFlag )
    {   radius = dst->radius*Scale;
        if( dist <= radius ) return;

        /* Test for second half obscured! */
        if( (col1!=col2) && (0.5*dist < radius) )
            col2 = col1;
    }

    if( src->flag & SphereFlag )
    {   radius = src->radius*Scale;
        if( dist <= radius ) return;

        /* Test for first half obscured! */
        if( (col1!=col2) && (0.5*dist < radius) )
            col1 = col2;
    }


    WriteVectColour( col1+ColourMask );
    temp = (src->z+dst->z)/(2.0*ImageSize)+1.0;
    if( temp != LineWidth )
    {   fprintf(OutFile,"%g setlinewidth\n",temp);
        LineWidth = temp;
    }


    if( src->flag & SphereFlag )
    {   dz = dst->z - src->z;
        dist = sqrt( (double)(dist2 + dz*dz) );
        endx = src->x + (radius*dx)/dist;
        endy = src->y + (radius*dy)/dist;
        fprintf(OutFile,"%g %g ",endx,endy);
    } else
        fprintf(OutFile,"%d %d ",src->x,src->y);

    if( col1 != col2 )
    {   midx = 0.5*(src->x + dst->x);
        midy = 0.5*(src->y + dst->y);
        fprintf(OutFile,"%g %g Wire\n",midx,midy);

        WriteVectColour( col2+ColourMask );
        fprintf(OutFile,"%g %g ",midx,midy);
    } 
    fprintf(OutFile,"%d %d Wire\n",dst->x,dst->y);
}


static void WriteVectStick( ptr )
    Bond __far *ptr;
{
    register Atom __far *src;
    register Atom __far *dst;
    register Real midx, midy;
    register Real relx, rely;
    register Real endx, endy;
    register Real radius, angle;
    register Real ratio, dist;
    register Real temp;

    register Long dist2;
    register int dx, dy, dz;
    register int col1, col2;
    register int i;

    if( !ptr->radius )
    {   WriteVectWire(ptr);
        return;
    }

    src = ptr->srcatom;  
    dst = ptr->dstatom;
    if( src->z > dst->z )
    {   src = ptr->dstatom;
        dst = ptr->srcatom;
    }

    if( !ptr->col )
    {   col1 = src->col;
        col2 = dst->col;
    } else col1 = col2 = ptr->col;

    dx = dst->x - src->x;  
    dy = dst->y - src->y;
    dist2 = dx*dx + dy*dy;
    dist = sqrt( (double)dist2 );

    if( dst->flag & SphereFlag )
    {   radius = dst->radius*Scale;
        if( dist <= radius ) return;

        /* Test for nearest half obscured! */
        if( (col1!=col2) && (0.5*dist < radius) )
            col2 = col1;
    }

    if( src->flag & SphereFlag )
    {   radius = src->radius*Scale;
        if( dist <= radius ) return;

        /* Test for furthest half obscured! */
        if( (col1!=col2) && (0.5*dist < radius) )
            col1 = col2;
    }

    temp = (src->z+dst->z)/(2.0*ImageSize)+1.0;
    if( temp != LineWidth )
    {   fprintf(OutFile,"%g setlinewidth\n",temp);
        LineWidth = temp;
    }

    radius = ptr->radius*Scale;
    angle = Rad2Deg*atan2((double)dy,(double)dx);

    if( col1 != col2 )
    {   midx = 0.5*(src->x + dst->x);
        midy = 0.5*(src->y + dst->y);
        relx = (radius*dx)/dist;
        rely = (radius*dy)/dist;

        fprintf(OutFile,"%g %g moveto\n",midx+rely,midy-relx);
        fprintf(OutFile,"%g %g lineto\n",midx-rely,midy+relx);

        dz = dst->z - src->z;
        dist = sqrt( (double)(dist2 + dz*dz) );
        ratio = dz/dist;

        if( (src->flag&SphereFlag) && (src->radius>ptr->radius) )
        {   temp = (Scale*src->radius)/dist;
            endx = src->x + temp*dx;
            endy = src->y + temp*dy;

            fprintf(OutFile,"%g %g %g ",radius,ratio,angle);
            fprintf(OutFile,"%g %g StickEnd\n",endx,endy);
        } else
        {   fprintf(OutFile,"%d %d %g ",src->x,src->y,radius);
            fprintf(OutFile,"%g %g arc\n",angle+90,angle-90);
        }
        fputs("closepath ",OutFile);

        i = col1 + ColourMask;
        fprintf(OutFile,"%g ",RLut[i]/255.0);
        fprintf(OutFile,"%g ",GLut[i]/255.0);
        fprintf(OutFile,"%g ",BLut[i]/255.0);
        fputs("setrgbcolor fill\n",OutFile);

        fprintf(OutFile,"%d %d %g ",dst->x,dst->y,radius);
        fprintf(OutFile,"%g %g arc\n",angle-90,angle+90);
        fprintf(OutFile,"%g %g %g ",radius,ratio,angle);
        fprintf(OutFile,"%g %g StickEnd\n",midx,midy);
        fputs("closepath ",OutFile);

        i = col2 + ColourMask;
        fprintf(OutFile,"%g ",RLut[i]/255.0);
        fprintf(OutFile,"%g ",GLut[i]/255.0);
        fprintf(OutFile,"%g ",BLut[i]/255.0);
        fputs("setrgbcolor fill\n",OutFile);

        fprintf(OutFile,"%d %d %g ",dst->x,dst->y,radius);
        fprintf(OutFile,"%g %g arc\n",angle-90,angle+90);
        if( (src->flag&SphereFlag) && (src->radius>ptr->radius) )
        {   fprintf(OutFile,"%g %g %g ",radius,ratio,angle);
            fprintf(OutFile,"%g %g StickEnd\n",endx,endy);
        } else
        {   fprintf(OutFile,"%d %d %g ",src->x,src->y,radius);
            fprintf(OutFile,"%g %g arc\n",angle+90,angle-90);
        }
        fputs("closepath 0 setgray stroke\n",OutFile);

    } else /* col1 == col2! */
    {   fprintf(OutFile,"%d %d %g ",dst->x,dst->y,radius);
        fprintf(OutFile,"%g %g arc\n",angle-90,angle+90);

        if( (src->flag&SphereFlag) && (src->radius>ptr->radius) )
        {   dz = dst->z - src->z;
            dist = sqrt( (double)(dist2 + dz*dz) );
            temp = (Scale*src->radius)/dist;
            endx = src->x + temp*dx;
            endy = src->y + temp*dy;
            ratio = dz/dist;

            fprintf(OutFile,"%g %g %g ",radius,ratio,angle);
            fprintf(OutFile,"%g %g StickEnd\n",endx,endy);
        } else
        {   fprintf(OutFile,"%d %d %g ",src->x,src->y,radius);
            fprintf(OutFile,"%g %g arc\n",angle+90,angle-90);
        }

        i = col1 + ColourMask;
        fprintf(OutFile,"%g ",RLut[i]/255.0);
        fprintf(OutFile,"%g ",GLut[i]/255.0);
        fprintf(OutFile,"%g ",BLut[i]/255.0);
        fputs("Stick\n",OutFile);
    }
    VectCol = 0;
}


int WriteVectPSFile( name )
    char *name;
{
    register Real ambi;
    register Real temp, inten;
    register int xsize, ysize;
    register int xpos, ypos;
    register Long count;
    register int i;

    PSItemPtr __far *data;
    char __far *type;

    register Chain __far *chain;
    register Group __far *group;
    register Bond __far *bptr;
    register Atom __far *aptr;


    /* Determine the number of objects to draw! */
    count = 0;
    ForEachAtom if( aptr->flag&SphereFlag ) count++;
    ForEachBond if( bptr->flag&DrawBondFlag ) count++;
    ForEachBack if( bptr->flag&DrawBondFlag ) count++;

    if( !count ) 
        return( True );

#ifdef IBMPC
    if( count > 16383 )
    {   if( CommandActive ) WriteChar('\n');
        WriteString("Output Error: Too many PostScript objects!\n");
        CommandActive = False;
        return( False );
    }
#endif

    /* Allocate arrays for objects! */
    data = (PSItemPtr __far*)_fmalloc((size_t)count*sizeof(PSItemPtr));
    type = (char __far*)_fmalloc((size_t)count*sizeof(char));
    if( !data || !type )
    {   if( CommandActive ) WriteChar('\n');
        WriteString("Output Error: Not enough memory to create PostScript!\n");
        CommandActive = False;

        if( data ) _ffree( data );
        if( type ) _ffree( type );
        return( False );
    }

    OutFile = fopen(name,"w");
    if( !OutFile )
    {   FatalOutputError(name);
        return(False);
    }

    /* Determine the size of the image */
    ysize = (int)(YRange*(BORDER*PAGEWIDE)/XRange);
    if( ysize > (int)(BORDER*PAGEHIGH) )
    {   xsize = (int)(XRange*(BORDER*PAGEHIGH)/YRange);
        ysize = (int)(BORDER*PAGEHIGH);
    } else xsize = (int)(BORDER*PAGEWIDE);

    xpos = (int)(PAGEWIDE-xsize)/2;
    ypos = (int)(PAGEHIGH-ysize)/2;

    fputs("%!PS-Adobe-2.0 EPSF-2.0\n",OutFile);
    fputs("%%Creator: RasMol Version 2.3\n",OutFile);
    fprintf(OutFile,"%%%%Title: %s\n",name);
    fprintf(OutFile,"%%%%BoundingBox: %d %d ",xpos,ypos);
    fprintf(OutFile,"%d %d\n",xpos+xsize,ypos+ysize);

    fputs("%%Pages: 1\n",OutFile);
    fputs("%%EndComments\n",OutFile);
    fputs("%%EndProlog\n",OutFile);
    fputs("%%BeginSetup\n",OutFile);

    fputs("1 setlinecap 1 setlinejoin [] 0 setdash\n",OutFile);
    fputs("1 setlinewidth 0 setgray\n",OutFile);
    fputs("%%EndSetup\n",OutFile);
    fputs("%%Page: 1 1\n",OutFile);

    fputs("gsave\n",OutFile);
    fputs("5 dict begin\n\n",OutFile);
    fputs("/Inten {\n  dup 4 index mul exch\n",OutFile);
    fputs("  dup 4 index mul exch\n",OutFile);
    fputs("  3 index mul setrgbcolor\n} def\n\n",OutFile);

    fputs("/Wire {\n  moveto lineto stroke\n} def\n\n",OutFile);

    fputs("/Stick {\n  closepath gsave setrgbcolor fill\n",OutFile);
    fputs("  grestore 0 setgray stroke\n} def\n\n",OutFile);
    fputs("/StickEnd {\n  matrix currentmatrix 6 1 roll\n",OutFile);
    fputs("  translate rotate 1 scale\n",OutFile);
    fputs("  0 0 3 2 roll 90 -90 arc\n  setmatrix\n} def\n\n",OutFile);

    ambi = 0.5*Ambient;
    fputs("/Sphere {\n  gsave translate",OutFile);
    fputs(" 0 0 2 index 0 360 arc gsave\n",OutFile);
    fputs("  45 rotate dup -0.81649658092 mul scale\n",OutFile);
    fprintf(OutFile,"  gsave %g Inten fill grestore clip\n",ambi);
    inten = (1.0-ambi)/31;
    for( i=0; i<31; i++ )
    {   temp = (Real)(i+1)/32;
        fprintf(OutFile,"  newpath 0 %g ",(Real)i/32);
        fprintf(OutFile,"%g 0 360 arc ",sqrt(1.0-temp*temp));
        fprintf(OutFile,"%g Inten fill\n",i*inten+ambi);
    }
    fputs("  grestore 0 setgray stroke\n",OutFile);
    fputs("  pop pop pop grestore\n} def\n\n",OutFile);

    fprintf(OutFile,"%d %d translate\n",xpos,ypos+ysize);
    fprintf(OutFile,"%g ",(Real)xsize/XRange);
    fprintf(OutFile,"%g ",(Real)-ysize/YRange);
    fputs("scale\n\n",OutFile);

    fputs("newpath 0 0 moveto 0 ",OutFile);
    fprintf(OutFile,"%d rlineto %d 0 rlineto 0 %d",YRange,XRange,-YRange);
    fputs(" rlineto closepath\ngsave ",OutFile);
    if( BackR || BackG || BackB )
    {   fprintf(OutFile,"%g %g %g",BackR/255.0,BackG/255.0,BackB/255.0);
        fputs(" setrgbcolor fill grestore\ngsave ",OutFile);
    }
    fputs("clip newpath\n\n",OutFile);

    LineWidth = 1.0;
    VectSolid = True;
    VectCol = 0;

    i = 0;
    ForEachAtom
        if( aptr->flag&SphereFlag )
        {   type[i] = PSSphere; data[i] = aptr; i++;
        }

    ForEachBond
       if( bptr->flag&CylinderFlag )
       {   type[i] = PSStick; data[i] = bptr; i++;
       } else if( bptr->flag&WireFlag )
       {   type[i] = PSLine;  data[i] = bptr; i++;
       }

    ForEachBack
       if( bptr->flag&CylinderFlag )
       {   type[i] = PSStick; data[i] = bptr; i++;
       } else if( bptr->flag&WireFlag )
       {   type[i] = PSLine;  data[i] = bptr; i++;
       }

    if( count>1 )
        DepthSort(data,type,(int)count);

    for( i=0; i<count; i++ )
        if( type[i]==PSSphere )
        {   WriteVectSphere( data[i] );
        } else if( type[i]==PSStick )
        {   WriteVectStick( data[i] );
        } else /* PSWire */
            WriteVectWire( data[i] );
 
     

    fputs("grestore stroke end grestore\n",OutFile);
    fputs("showpage\n",OutFile);
    fputs("%%Trailer\n",OutFile);
    fputs("%%EOF\n",OutFile);

    fclose( OutFile );
    _ffree( data );
    _ffree( type );
    return(True);
}


void InitialiseOutFile()
{
#ifdef IBMPC
    /* Allocate Tables on FAR Heap */
    ABranch = (short __far*)_fmalloc(4096*sizeof(short));
    DBranch = (short __far*)_fmalloc(4096*sizeof(short));
    Hash = (short __far*)_fmalloc(256*sizeof(short));
    Node = (Byte __far*)_fmalloc(4096*sizeof(Byte));
#endif
}

