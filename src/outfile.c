/* outfile.c
 * RasMol2 Molecular Graphics
 * Roger Sayle, October 1994
 * Version 2.5
 */
#define OUTFILE
#include "rasmol.h"

#ifdef IBMPC
#include <windows.h>
#include <malloc.h>
#endif
#ifdef APPLEMAC
#ifdef __CONDITIONALMACROS__
#include <Printing.h>
#else
#include <PrintTraps.h>
#endif
#include <Types.h>
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
#include "abstree.h"
#include "transfor.h"
#include "render.h"
#include "graphics.h"
#include "pixutils.h"
#include "script.h"


#ifdef EIGHTBIT
#define RComp(x)   (RLut[LutInv[x]])
#define GComp(x)   (GLut[LutInv[x]])
#define BComp(x)   (BLut[LutInv[x]])
#else
#define RComp(x)   (((x)>>16)&0xff)
#define GComp(x)   (((x)>>8)&0xff)
#define BComp(x)   ((x)&0xff)
#endif


#ifdef INVERT
#define InvertY(y) (y)
#else
#define InvertY(y) (-(y))
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

#define Round(x)       ((int)(x))


#define PSLine    0x00
#define PSStick   0x01
#define PSSphere  0x02
#define PSRibbon  0x03

typedef void __far* PSItemPtr;


#if defined(IBMPC) || defined(APPLEMAC)
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


/* Apple PICT macros */
#define PICTcliprgn         0x0001
#define PICTpicversion      0x0011
#define PICTpackbitsrect    0x0098
#define PICTdirectbitsrect  0x009a
#define PICTendofpict       0x00ff
#define PICTheaderop        0x0c00

static Byte PICTBuffer[128];
static Byte PICTPixel;

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



#ifdef APPLEMAC
/* External RasMac Function Declaration! */
void SetFileInfo( char*, OSType, OSType, short );
#endif

static void FatalOutputError( ptr )
    char *ptr;
{
    if( CommandActive ) WriteChar('\n');
    WriteString("Output Error: Unable to create file `");
    WriteString( ptr );  WriteString("'!\n");
    CommandActive = False;
}


/*  Integer Output Routines  */
#ifdef FUNCPROTO
static void WriteLSBLong( Card );
static void WriteMSBLong( Card );
#endif

static void WriteLSBShort( val )
    int val;
{
    putc( val&0xff, OutFile );
    putc( (val>>8)&0xff, OutFile );
}

static void WriteMSBShort( val )
    int val;
{
    putc( (val>>8)&0xff, OutFile );
    putc( val&0xff, OutFile );
}

static void WriteLSBLong( val )
    Card val;
{
    putc((int)(val&0xff),OutFile);
    putc((int)((val>>8) &0xff),OutFile);
    putc((int)((val>>16)&0xff),OutFile);
    putc((int)((val>>24)&0xff),OutFile);
}

static void WriteMSBLong( val )
    Card val;
{
    putc((int)((val>>24)&0xff),OutFile);
    putc((int)((val>>16)&0xff),OutFile);
    putc((int)((val>>8) &0xff),OutFile);
    putc((int)(val&0xff),OutFile);
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

#if defined(IBMPC) || defined(APPLEMAC)
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
#ifdef APPLEMAC
    /* Avoid ANSI trigraph problems! */
    SetFileInfo(name,'\?\?\?\?','\?\?\?\?',134);
#endif
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

#if defined(IBMPC) || defined(APPLEMAC)
    OutFile = fopen(name,"wb");
#else
    OutFile = fopen(name,"w");
#endif
    if( !OutFile ) 
    {    FatalOutputError(name);
	 return( False );
    }
    fwrite("GIF87a",1,6,OutFile);
    WriteLSBShort(XRange);
    WriteLSBShort(YRange);
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
    WriteLSBShort(XRange);
    WriteLSBShort(YRange);
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

#ifdef APPLEMAC
    /* Avoid ANSI trigraph problems! */
    SetFileInfo(name,'\?\?\?\?','GIFf',134);
#endif
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

#if defined(IBMPC) || defined(APPLEMAC)
    OutFile = fopen(name,"wb");
#else
    OutFile = fopen(name,"w");
#endif

    if( !OutFile )
    {   FatalOutputError(name);
	return(False);
    }
    WriteMSBLong( RAS_MAGIC );
    WriteMSBLong(XRange);  
    WriteMSBLong(YRange);
    RLEEncode = encode;


#ifdef EIGHTBIT
    WriteMSBLong(8);

    if( encode )
    {   WriteRastData(False);
	WriteMSBLong(RLEFileSize);
	WriteMSBLong(RT_BYTE_ENCODED);
    } else
    {   size = (XRange%2)? XRange+1 : XRange;
	WriteMSBLong(size*YRange);
	WriteMSBLong(RT_STANDARD);
    }

    cols = CompactColourMap();
    WriteMSBLong(RMT_EQUAL_RGB);
    WriteMSBLong(cols*3);

    for( i=0; i<cols; i++ ) fputc((int)RLut[Node[i]],OutFile);
    for( i=0; i<cols; i++ ) fputc((int)GLut[Node[i]],OutFile);
    for( i=0; i<cols; i++ ) fputc((int)BLut[Node[i]],OutFile);
#else
    WriteMSBLong(24);

    if( encode )
    {   WriteRastData(False);
	WriteMSBLong(RLEFileSize);
	WriteMSBLong(RT_BYTE_ENCODED);
    } else
    {   size = XRange*3;
	if( size&1 ) size++;
	WriteMSBLong(size*YRange);
	WriteMSBLong(RT_STANDARD);
    }
    WriteMSBLong(RMT_NONE);
    WriteMSBLong(0);
#endif

    WriteRastData(True);
    fclose( OutFile );
#ifdef APPLEMAC
    /* Avoid ANSI trigraph problems! */
    SetFileInfo(name,'\?\?\?\?','\?\?\?\?',134);
#endif
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

int WriteEPSFFile( name, col, compr )
    char *name;  int col, compr;
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
    fputs("%%Creator: RasMol Version 2.5\n",OutFile);
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

    if( compr )
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

    if( !compr )
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
	    if( compr )
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

    if( compr && RLELength )
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
#ifdef APPLEMAC
    /* Avoid ANSI trigraph problems! */
    SetFileInfo(name,'\?\?\?\?','EPSF',134);
#endif
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
	    WriteBMPCode(LutInv[*ptr++]);

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


    
int WriteBMPFile( name )
    char *name;
{
    register Pixel __huge *ptr;
    register int x,y,i,raw;
    register Card size;

#if defined(IBMPC) || defined(APPLEMAC)
    OutFile = fopen(name,"wb");
#else
    OutFile = fopen(name,"w");
#endif
    if( !OutFile )
    {    FatalOutputError(name);
	 return( False );
    }

#ifdef EIGHTBIT
    for( i=0; i<256; i++ )
	if( ULut[i] )
	    LutInv[Lut[i]] = i;

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
    WriteLSBLong(size+1078);
    WriteLSBLong((Card)0);
    WriteLSBLong((Card)1078);

    WriteLSBLong((Card)40);
    WriteLSBLong((Card)XRange);
    WriteLSBLong((Card)YRange);
    fputc(0x01,OutFile);  fputc(0x00,OutFile);
    fputc(0x08,OutFile);  fputc(0x00,OutFile);
    WriteLSBLong(raw? (Card)0 : (Card)1);
    WriteLSBLong(size);
    
    WriteLSBLong((Card)0);
    WriteLSBLong((Card)0);
    WriteLSBLong((Card)256);
    WriteLSBLong((Card)256);

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
		fputc((int)LutInv[*ptr++],OutFile);
	    for( x=0; x<BMPPad; x++ )
		fputc(0x00,OutFile);
	    fputc(0x00,OutFile);
	    fputc((y?0x00:0x01),OutFile);
	}
    } else
	WriteBMPData(True);
    fclose(OutFile);

#ifdef APPLEMAC
    /* Avoid ANSI trigraph problems! */
    SetFileInfo(name,'\?\?\?\?','\?\?\?\?',134);
#endif
#ifdef IBMPC
    GlobalUnlock(FBufHandle);
#endif
    return(True);
#else
    if( CommandActive ) WriteChar('\n');
    WriteString("Output Error: 24 bit BMP files unsupported!\n");
    CommandActive = False;
    return( False );
#endif
}



#ifdef FUNCPROTO
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

#ifdef FUNCPROTO
static int ClipVectSphere( Atom __far* );
static int ClipVectBond( Bond __far* );

static void WriteVectSphere( PSItemPtr __far*, char __far*, int );
static void WriteVectStick( Bond __far* );
static void WriteVectWire( Bond __far* );
#endif


static int ClipVectSphere( ptr )
    Atom __far *ptr;
{
    register int rad;

    rad = ptr->irad;

    if( ptr->x + rad < 0 )  return( True );
    if( ptr->y + rad < 0 )  return( True );
    if( ptr->x - rad >= XRange )  return( True );
    if( ptr->y - rad >= YRange )  return( True );
    return( False );
}


static int ClipVectBond( ptr )
    Bond __far *ptr;
{
    register Atom __far *src;
    register Atom __far *dst;

    src = ptr->srcatom;
    dst = ptr->dstatom;

    if( (src->x<0) && (dst->x<0) )  return( True );
    if( (src->y<0) && (dst->y<0) )  return( True );
    if( (src->x>=XRange) && (dst->x>=XRange) )  return( True );
    if( (src->y>=YRange) && (dst->y>=YRange) )  return( True );
    return( False );
}



static void WriteVectColour( col )
    int col;
{
    if( col != VectCol )
    {   fprintf(OutFile,"%g ",(Real)RLut[col]/255.0);
	fprintf(OutFile,"%g ",(Real)GLut[col]/255.0);
	fprintf(OutFile,"%g ",(Real)BLut[col]/255.0);
	fputs("setrgbcolor\n",OutFile);
	VectCol = col;
    }
}


#define MAXSECT 5
typedef struct {
	/* Ellipse */
	Real ephi,epsi;
	Real etheta;
	Real ex,ey;
	Real erad;

	/* Sphere */
	Real sphi,spsi;
	int sx,sy;
	Real srad;
    } SphereSect;


static int VectClipContain( x, y )
    SphereSect *x; SphereSect *y;
{
    if( x->erad != 0.0 )
    {   if( y->erad != 0.0 )
	    /* Simple segment containment test! */
	    return( ((x->sphi+x->spsi)>=(y->sphi+y->spsi)) &&
		    ((x->sphi-x->spsi)<=(y->sphi-y->spsi)) );
    } else if( y->erad == 0.0 )
	return( x->srad >= y->srad );
    return( False );
}


static void WriteVectSphere( data, type, index )
    PSItemPtr __far*data; 
    char __far *type;
    int index;
{
    register int ecount, count;
    register Atom __far *atm;
    register Atom __far *ptr;
    register Long dist2,dist3;
    register int dx, dy, dz;
    register int i,j,k;

    register Real b,d,f,g,x;
    register Real radf,radb;
    register Real phi1,phi2;
    register Real temp,psi;
    register Real theta;

    register SphereSect *sptr;
    SphereSect sect[MAXSECT];

    ptr = (Atom __far*)data[index];
    radf = ptr->radius*Scale;

    count = 0;
    ecount = 0;
    sptr = sect;
    for( i=index-1; i>=0; i-- )
    {   if( type[i] != PSSphere )
	    continue;

	atm = (Atom __far*)data[i];
	/* Atom can't intersect visibly! */
	if( atm->z + atm->irad < ptr->z )
	    continue;

	dx = atm->x - ptr->x; 
	dy = atm->y - ptr->y; 
	dz = atm->z - ptr->z;

	dist2 = (Long)dx*dx + (Long)dy*dy;
	dist3 = dist2 + dz*dz;

	radb = atm->radius*Scale;  
	temp = radf + radb;

	/* Atoms don't intersect! */
	if( dist3 > temp*temp ) continue;


	d = sqrt( (double)dist3 );
	f = (temp*(radf-radb)+dist3)/(2.0*d);
	theta = -dz/d;

	if( f>0 )
	{   temp = radf*radf;
	    /* Intersection not visible! */
	    if( theta*temp > temp-f*f )
		continue;
	} else if( f < -radf )
	    return;

	x = sqrt( (radf-f)*(radf+f) );

	if( dx || dy )
	{   g = sqrt( (double)dist2 );
	    psi = Rad2Deg*atan2(dy,dx);
	    b = (f*(dz*dz))/(d*g);

	    if( AbsFun(b)>x )
		continue;

	    phi1 = b + (f*g)/d;
	    phi1 = Rad2Deg*acos(phi1/radf);
	    if( phi1!=phi1 ) continue;

	    phi2 = (d*b)/g;
	    if( AbsFun(phi2) < x )
	    {   phi2 = Rad2Deg*acos(phi2/x);
		if( phi2!=phi2 ) continue;
		if( phi2 > 90.0 ) 
		    phi2 = 180.0-phi1;
	    } else phi2 = 90.0;

	    sptr->erad = x;
	    sptr->etheta = -theta;
	    sptr->ephi = psi;
	    sptr->epsi = phi2;

	    temp = f/d;
	    sptr->ex = ptr->x+temp*dx;
	    sptr->ey = ptr->y+temp*dy;

	    sptr->srad = radf;
	    sptr->sphi = psi;
	    sptr->spsi = phi1;
	    sptr->sx = ptr->x;
	    sptr->sy = ptr->y;

	} else
	{   x = sqrt( (radf-g)*(radf+g) );

	    sptr->srad = x;
	    sptr->erad = 0.0;
	    sptr->sx = ptr->x;
	    sptr->sy = ptr->y;
	    sptr->sphi = 180;
	    sptr->spsi = -180;
	}

	/* Optimize Segments */
	j = 0;
	while( j<count )
	    if( VectClipContain(sptr,sect+j) )
	    {   /* Delete Segment sect[j] */
		for( k=j; k<count; k++ )
		    sect[k] = sect[k+1];
		count--;  sptr--;
	    } else if( VectClipContain(sect+j,sptr) )
	    {   break;  /* Exclude Segment */
	    } else j++;
	   

	if( j==count )
	{   count++;  sptr++;
	    if( sptr->erad != 0.0 )
		ecount++;
	    if( count==MAXSECT )
		break;
	}
    }

    if( UseOutLine )
    {   temp = (ptr->z-ZOffset)/ImageSize + 1.0;
	if( temp != LineWidth )
	{   fprintf(OutFile,"%g setlinewidth\n",temp);
	    LineWidth = temp;
	}
    }

    if( count )
    {   fputs("gsave\n",OutFile);
	fprintf(OutFile,"%%%% %d %d\n",count,ecount);

	sptr = sect;
	for( i=0; i<count; i++ )
	{   if( sptr->erad != 0.0 )
	    {   fprintf(OutFile,"%g %g %g %g %g %g ClipEllips\n",
			    sptr->erad,sptr->epsi,sptr->etheta,
			    sptr->ephi,sptr->ex,sptr->ey);
	    }

	    if( (i==count-1) || (sptr->erad==0.0) )
	    {   fprintf(OutFile,"%g %g %g %d %d ClipSphere\n",sptr->srad,
				sptr->sphi+sptr->spsi,sptr->sphi-sptr->spsi,
				sptr->sx, sptr->sy );
	    } else fprintf(OutFile,"%g %g %g %d %d ClipBox\n",
				    sptr->srad+sptr->srad+2,
				    sptr->srad+1, sptr->ephi,
				    sptr->sx, sptr->sy );
	    sptr++;
	}

	i = ptr->col + ColourMask;
	fprintf(OutFile,"%g ",(Real)RLut[i]/255.0);
	fprintf(OutFile,"%g ",(Real)GLut[i]/255.0);
	fprintf(OutFile,"%g ",(Real)BLut[i]/255.0);
	fprintf(OutFile,"%g Shade\n",radf);
	fputs("grestore\n\n",OutFile);
    } else
    {   i = ptr->col + ColourMask;
	fprintf(OutFile,"%g ",(Real)RLut[i]/255.0);
	fprintf(OutFile,"%g ",(Real)GLut[i]/255.0);
	fprintf(OutFile,"%g ",(Real)BLut[i]/255.0);
	fprintf(OutFile,"%g %d %d ",radf,ptr->x,ptr->y);
	fputs("Sphere\n\n",OutFile);
    }
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
    temp = ((src->z-ZOffset)+(dst->z-ZOffset))/ImageSize + 1.0;
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
    register Real dist, dist3;
    register Real temp, ratio;

    register Long dist2;
    register int dx, dy, dz;
    register int col1, col2;
    register int i, inten;

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
    dz = dst->z - src->z;
    dist2 = dx*dx + dy*dy;
    dist3 = sqrt( (double)(dist2 + dz*dz) );
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

    temp = ((src->z-ZOffset)+(dst->z-ZOffset))/ImageSize + 1.0;
    if( temp != LineWidth )
    {   fprintf(OutFile,"%g setlinewidth\n",temp);
	LineWidth = temp;
    }

    radius = ptr->radius*Scale;
    angle = Rad2Deg*atan2((double)dy,(double)dx);
    inten = (int)((dist/dist3)*ColourMask);

    if( col1 != col2 )
    {   midx = 0.5*(src->x + dst->x);
	midy = 0.5*(src->y + dst->y);
	relx = (radius*dx)/dist;
	rely = (radius*dy)/dist;

	fprintf(OutFile,"%g %g moveto\n",midx+rely,midy-relx);
	fprintf(OutFile,"%g %g lineto\n",midx-rely,midy+relx);

	ratio = dz/dist3;

	if( (src->flag&SphereFlag) && (src->radius>ptr->radius) )
	{   temp = (Scale*src->radius)/dist3;
	    endx = src->x + temp*dx;
	    endy = src->y + temp*dy;

	    fprintf(OutFile,"%g %g %g ",radius,ratio,angle);
	    fprintf(OutFile,"%g %g StickEnd\n",endx,endy);
	} else
	{   fprintf(OutFile,"%d %d %g ",src->x,src->y,radius);
	    fprintf(OutFile,"%g %g arc\n",angle+90,angle-90);
	}
	fputs("closepath ",OutFile);

	i = col1 + inten;
	fprintf(OutFile,"%g ",(Real)RLut[i]/255.0);
	fprintf(OutFile,"%g ",(Real)GLut[i]/255.0);
	fprintf(OutFile,"%g ",(Real)BLut[i]/255.0);
	fputs("setrgbcolor fill\n",OutFile);

	fprintf(OutFile,"%d %d %g ",dst->x,dst->y,radius);
	fprintf(OutFile,"%g %g arc\n",angle-90,angle+90);
	fprintf(OutFile,"%g %g %g ",radius,ratio,angle);
	fprintf(OutFile,"%g %g StickEnd\n",midx,midy);
	fputs("closepath ",OutFile);

	i = col2 + inten;
	fprintf(OutFile,"%g ",(Real)RLut[i]/255.0);
	fprintf(OutFile,"%g ",(Real)GLut[i]/255.0);
	fprintf(OutFile,"%g ",(Real)BLut[i]/255.0);
	fputs("setrgbcolor fill\n",OutFile);

	if( UseOutLine )
	{   fprintf(OutFile,"%d %d %g ",dst->x,dst->y,radius);
	    fprintf(OutFile,"%g %g arc\n",angle-90,angle+90);
	    if( (src->flag&SphereFlag) && (src->radius>ptr->radius) )
	    {   fprintf(OutFile,"%g %g %g ",radius,ratio,angle);
		fprintf(OutFile,"%g %g StickEnd\n",endx,endy);
	    } else
	    {   fprintf(OutFile,"%d %d %g ",src->x,src->y,radius);
		fprintf(OutFile,"%g %g arc\n",angle+90,angle-90);
	    }
	    fputs("closepath 0 setgray stroke\n",OutFile);
	}
    } else /* col1 == col2! */
    {   fprintf(OutFile,"%d %d %g ",dst->x,dst->y,radius);
	fprintf(OutFile,"%g %g arc\n",angle-90,angle+90);

	if( (src->flag&SphereFlag) && (src->radius>ptr->radius) )
	{   temp = (Scale*src->radius)/dist3;
	    endx = src->x + temp*dx;
	    endy = src->y + temp*dy;
	    ratio = dz/dist3;

	    fprintf(OutFile,"%g %g %g ",radius,ratio,angle);
	    fprintf(OutFile,"%g %g StickEnd\n",endx,endy);
	} else
	{   fprintf(OutFile,"%d %d %g ",src->x,src->y,radius);
	    fprintf(OutFile,"%g %g arc\n",angle+90,angle-90);
	}

	i = col1 + inten;
	fprintf(OutFile,"%g ",(Real)RLut[i]/255.0);
	fprintf(OutFile,"%g ",(Real)GLut[i]/255.0);
	fprintf(OutFile,"%g ",(Real)BLut[i]/255.0);
	fputs("Stick\n",OutFile);
    }
    VectCol = 0;
}


static void WriteVectDots()
{
    register DotStruct __far *ptr;
    register Real x,y,z;
    register Real xi,yi;
    register int inten;
    register int temp;
    register int zi;
    register int i;


    if( LineWidth != 1.0 )
    {   fputs("1 setlinewidth\n",OutFile);
	LineWidth = 1.0;
    }

    temp = SlabValue - ZOffset;
    for( ptr=DotPtr; ptr; ptr=ptr->next )
	for( i=0; i<ptr->count; i++ )
	{   x = ptr->xpos[i];  
	    y = ptr->ypos[i];  
	    z = ptr->zpos[i];


	    xi = (x*MatX[0]+y*MatX[1]+z*MatX[2]) + XOffset;
	    if( (xi<0.0) || (xi>=XRange) ) continue;
	    yi = (x*MatY[0]+y*MatY[1]+z*MatY[2]) + YOffset;
	    if( (yi<0.0) || (yi>=YRange) ) continue;

	    zi = (int)(x*MatZ[0]+y*MatZ[1]+z*MatZ[2]);
	    if( UseSlabPlane && (zi>=temp) ) continue;

	    inten = (ColourDepth*(zi+ImageRadius))/ImageSize;
	    WriteVectColour( ptr->col[i]+inten );
	    fprintf(OutFile,"%g %g Dot\n",xi,yi);
	}
}


static void WriteVectLabels()
{
    register Chain __far *chain;
    register Group __far *group;
    register Atom __far *aptr;
    register Label *label;
    auto char buffer[80];

    fputs("/Times-Roman",OutFile); /* Courier or Courier-Bold? */
    fprintf(OutFile," findfont %d scalefont setfont\n",FontSize<<1);

    if( UseLabelCol )
    {   if( BackR || BackG || BackB )
	{   fprintf(OutFile,"%g %g %g setrgbcolor\n",
		    LabR/250.0, LabG/250.0, LabB/250.0);
	} else fputs("0 setgray\n",OutFile);
    } else VectCol = 0;

    ForEachAtom
	if( aptr->label )
	{   if( !UseLabelCol && (aptr->col!=VectCol) )
		 WriteVectColour( aptr->col );

	    label = (Label*)aptr->label;
	    FormatLabel(chain,group,aptr,label->label,buffer);
	    fprintf(OutFile,"(%s) %d %d Label\n",buffer,aptr->x,aptr->y);
	}
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
    ForEachAtom 
	if( aptr->flag&SphereFlag ) 
	    if( !UseClipping || !ClipVectSphere(aptr) )
		count++;
    ForEachBond 
	if( bptr->flag&DrawBondFlag ) 
	    if( !UseClipping || !ClipVectBond(bptr) )
		count++;
    ForEachBack 
	if( bptr->flag&DrawBondFlag )
	    if( !UseClipping || !ClipVectBond(bptr) )
		count++;

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
    fputs("%%Creator: RasMol Version 2.5\n",OutFile);
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
    fputs("14 dict begin\n\n",OutFile);
    fputs("/handleerror { showpage } def\n\n",OutFile);
    fputs("/Inten {\n  dup 4 index mul exch\n",OutFile);
    fputs("  dup 4 index mul exch\n",OutFile);
    fputs("  3 index mul setrgbcolor\n} def\n\n",OutFile);

    fputs("/Dot {\n  moveto 0 0 rlineto stroke\n} def\n\n",OutFile);
    fputs("/Wire {\n  moveto lineto stroke\n} def\n\n",OutFile);
#ifdef INVERT
    fputs("/Label {\n  moveto show\n} def\n\n",OutFile);
#else
    fputs("/Label {\n  moveto 1 -1 scale\n",OutFile);
    fputs("  show mtrx setmatrix\n} def\n\n",OutFile);
#endif

    if( UseOutLine )
    {   fputs("/Stick {\n  closepath gsave setrgbcolor fill\n",OutFile);
	fputs("  grestore 0 setgray stroke\n} def\n\n",OutFile);
    } else
	fputs("/Stick {\n  closepath setrgbcolor fill\n} def\n\n",OutFile);

    fputs("/StickEnd {\n  matrix currentmatrix 6 1 roll\n",OutFile);
    fputs("  translate rotate 1 scale\n",OutFile);
    fputs("  0 0 3 2 roll 90 -90 arc\n  setmatrix\n} def\n\n",OutFile);

    ambi = 0.5*Ambient;
    if( UseOutLine )
    {   fputs("/Shade {\n  closepath gsave clip\n",OutFile);
    } else fputs("/Shade {\n  closepath clip\n",OutFile);
    fputs("  45 rotate dup -0.81649658092 mul scale\n",OutFile);
    fprintf(OutFile,"  %g Inten fill\n",ambi);
    inten = (1.0-ambi)/31;
    for( i=0; i<31; i++ )
    {   temp = (Real)(i+1)/32;
	fprintf(OutFile,"  0 %g ",(Real)i/32);
	fprintf(OutFile,"%g 0 360 arc ",sqrt(1.0-temp*temp));
	fprintf(OutFile,"%g Inten fill\n",i*inten+ambi);
    }
    if( UseOutLine )
    {   fputs("  grestore 0 setgray stroke pop pop pop\n} def\n\n",OutFile);
    } else fputs("  pop pop pop\n} def\n\n",OutFile);

    fputs("/ClipSphere {\n  translate 0 0 5 2 roll arc\n} def\n\n",OutFile);
    fputs("/ClipBox {\n  translate rotate\n  dup lineto dup neg ",OutFile);
    fputs("dup\n  0 rlineto 0 exch rlineto 0 rlineto closepath\n",OutFile);
    fputs("  clip newpath mtrx setmatrix\n} def\n\n",OutFile);
    fputs("/ClipEllips {\n  translate rotate 1 scale\n",OutFile);
    fputs("  0 0 4 2 roll dup neg arc\n",OutFile);
    fputs("  reversepath mtrx setmatrix\n} def\n\n",OutFile);

    fputs("/Sphere {\n  gsave\n",OutFile);
    fputs("  translate 0 0 2 index 0 360 arc\n",OutFile);
    if( UseOutLine )
    {   fputs("  gsave Shade grestore\n",OutFile);
	fputs("  0 setgray stroke\n  ",OutFile);
	fputs("  grestore\n} def\n\n",OutFile);
    } else
	fputs("  Shade grestore\n} def\n\n",OutFile);

#ifdef INVERT
    fprintf(OutFile,"%d %d translate\n",xpos,ypos);
    fprintf(OutFile,"%g ",(Real)xsize/XRange);
    fprintf(OutFile,"%g ",(Real)ysize/YRange);
#else
    fprintf(OutFile,"%d %d translate\n",xpos,ypos+ysize);
    fprintf(OutFile,"%g ",(Real)xsize/XRange);
    fprintf(OutFile,"%g ",(Real)-ysize/YRange);
#endif
    fputs("scale\n/mtrx matrix currentmatrix def\n\n",OutFile);

    fputs("newpath 0 0 moveto 0 ",OutFile);
    fprintf(OutFile,"%d rlineto %d 0 rlineto 0 %d",YRange,XRange,-YRange);
    fputs(" rlineto\nclosepath clip ",OutFile);
    if( BackR || BackG || BackB )
    {   fprintf(OutFile,"%g %g %g",BackR/255.0,BackG/255.0,BackB/255.0);
	fputs(" setrgbcolor fill\n\n",OutFile);
    } else fputs("newpath\n\n",OutFile);

    LineWidth = 1.0;
    VectSolid = True;
    VectCol = 0;

    i = 0;
    ForEachAtom
	if( aptr->flag&SphereFlag )
	    if( !UseClipping || !ClipVectSphere(aptr) )
	    {   type[i] = PSSphere; 
		data[i++] = aptr;
	    }

    ForEachBond
       if( bptr->flag&DrawBondFlag )
	   if( !UseClipping || !ClipVectBond(bptr) )
	   {   type[i] = (bptr->flag&CylinderFlag)? PSStick:PSLine;
	       data[i++] = bptr;
	   } 

    ForEachBack
       if( bptr->flag&DrawBondFlag )
	   if( !UseClipping || !ClipVectBond(bptr) )
	   {   type[i] = (bptr->flag&CylinderFlag)? PSStick:PSLine;
	       data[i++] = bptr; 
	   } 


    if( count>1 )
	DepthSort(data,type,(int)count);

    for( i=0; i<count; i++ )
	if( type[i]==PSSphere )
	{   WriteVectSphere( data, type, i );
	} else if( type[i]==PSStick )
	{   WriteVectStick( data[i] );
	} else /* PSWire */
	    WriteVectWire( data[i] );
 
    if( DrawDots )
	WriteVectDots();
    if( DrawLabels )
	WriteVectLabels();

    fputs("newpath 0 0 moveto 0 ",OutFile);
    fprintf(OutFile,"%d rlineto %d 0 rlineto 0 %d",YRange,XRange,-YRange);
    fputs(" rlineto\nclosepath 0 setgray 1 setlinewidth stroke\n",OutFile);
    fputs("end grestore\nshowpage\n",OutFile);
    fputs("%%Trailer\n",OutFile);
    fputs("%%EOF\n",OutFile);

    fclose( OutFile );
#ifdef APPLEMAC
    /* Avoid ANSI trigraph problems! */
    SetFileInfo(name,'\?\?\?\?','EPSF',134);
#endif
    _ffree( data );
    _ffree( type );
    return(True);
}


static void WritePICTPixel( i )
    int i;
{
#ifdef EIGHTBIT
    putc( (int)LutInv[i], OutFile );
#else
    putc( i, OutFile );
#endif
}


static int FlushPICTBuffer()
{
    register int result;
    register int i;

    if( RLEOutput )
    {   putc(PacketLen-1,OutFile);
	for( i=0; i<PacketLen; i++ )
	    WritePICTPixel( PICTBuffer[i] );
    }
    result = PacketLen+1;
    PacketLen = 0;
    return( result );
}


static int FlushPICTPacket()
{
    register int result;
    register int i;

    if( RLELength>2 )
    {   if( PacketLen )
	{   result = FlushPICTBuffer();
	} else result = 0;

	if( RLEOutput )
	{   putc(257-RLELength,OutFile);
	    WritePICTPixel(PICTPixel);
	}
	result += 2;
    } else 
    {   result = 0;
	for( i=0; i<RLELength; i++ )
	{   PICTBuffer[PacketLen++] = PICTPixel;
	    if( PacketLen==128 ) 
		result  += FlushPICTBuffer();
	}
    }
    return( result );
}


static int WritePICTCode( val )
    int val;
{
    register int result;

    result = 0;
    if( !RLELength )
    {   PICTPixel = val;
	    RLELength = 1;
    } else if( (val!=PICTPixel) || (RLELength==128) )
    {   result = FlushPICTPacket();
	    PICTPixel = val;
	    RLELength = 1;
    } else RLELength++;
    return( result );
}

static Card WritePICTData()
{
#ifndef EIGHTBIT
    register Pixel data;
#endif
    register Pixel __huge *ptr;
    register Pixel __huge *tmp;
    register int rowbytes;
    register Card result;
    register int size;
    register int x,y;

#ifdef IBMPC
    FBuffer = (Pixel __huge*)GlobalLock(FBufHandle);
#endif

#ifdef EIGHTBIT
    rowbytes = XRange;
#else
    rowbytes = XRange*3;
#endif

    result = 0;
#ifndef INVERT
    ptr = FBuffer;
#endif
    for( y=YRange-1; y>=0; y-- )
    {
#ifdef INVERT
	ptr = FBuffer + (Long)y*XRange;
#endif

	RLEOutput = False;
	PacketLen = 0;
	RLELength = 0;

	size = 0;
#ifdef EIGHTBIT
	tmp = ptr;
	for( x=0; x<XRange; x++ )
	    size += WritePICTCode( *tmp++ );
#else
	tmp = ptr;
	for( x=0; x<XRange; x++ )
	    size += WritePICTCode( (int)RComp(*tmp++) );
	tmp = ptr;
	for( x=0; x<XRange; x++ )
	    size += WritePICTCode( (int)GComp(*tmp++) );
	tmp = ptr;
	for( x=0; x<XRange; x++ )
	    size += WritePICTCode( (int)BComp(*tmp++) );
#endif
	size += FlushPICTPacket();
	if( PacketLen )
	    size += FlushPICTBuffer();

#ifdef EIGHTBIT
	if( rowbytes > 250 )
	{   WriteMSBShort(size);
	    result += (size+2);
	} else
	{   putc(size,OutFile);
	    result += (size+1);
	}
#else
	WriteMSBShort(size);
	result += (size+2);
#endif

	RLEOutput = True;
	PacketLen = 0;
	RLELength = 0;

#ifdef EIGHTBIT
	for( x=0; x<XRange; x++ )
	    (void)WritePICTCode( *ptr++ );
#else
	tmp = ptr;
	for( x=0; x<XRange; x++ )
	    (void)WritePICTCode( (int)RComp(*tmp++) );
	tmp = ptr;
	for( x=0; x<XRange; x++ )
	    (void)WritePICTCode( (int)GComp(*tmp++) );
	for( x=0; x<XRange; x++ )
	    (void)WritePICTCode( (int)BComp(*ptr++) );
#endif
	(void)FlushPICTPacket();
	if( PacketLen )
	    (void)FlushPICTBuffer();
    }

#ifdef IBMPC
    GlobalUnlock(FBufHandle);
#endif
    return(result);
}


int WritePICTFile( name )
    char *name;
{
#ifdef EIGHTBIT
    register int j,r,g,b;
#endif
    register Card size;
    register int cols;
    register int i;

#if defined(IBMPC) || defined(APPLEMAC)
    OutFile = fopen(name,"wb");
#else
    OutFile = fopen(name,"w");
#endif
    if( !OutFile )
    {    FatalOutputError(name);
	 return( False );
    }

#ifdef EIGHTBIT
    cols = CompactColourMap();
#endif

    /* Write out header */
    for( i=0; i<512; i++ )
	putc( 0, OutFile );

    size = 0;
    WriteMSBShort( size );  /* picSize */
    WriteMSBShort(0);       /* picFrame.top    */
    WriteMSBShort(0);       /* picFrame.left   */
    WriteMSBShort(YRange);  /* picFrame.bottom */
    WriteMSBShort(XRange);  /* picFrame.right  */

    WriteMSBShort(PICTpicversion);
    WriteMSBShort(0x02FF);

    WriteMSBShort(PICTheaderop);
    WriteMSBLong((Card)0xffffffff);
    WriteMSBShort(0);      WriteMSBShort(0);
    WriteMSBShort(0);      WriteMSBShort(0);
    WriteMSBShort(XRange); WriteMSBShort(0);
    WriteMSBShort(YRange); WriteMSBShort(0);
    WriteMSBLong(0);

    WriteMSBShort(PICTcliprgn);
    WriteMSBShort(10);      /* rgnSize */
    WriteMSBShort(0);       /* rgnBBox.top    */
    WriteMSBShort(0);       /* rgnBBox.left   */
    WriteMSBShort(YRange);  /* rgnBBox.bottom */
    WriteMSBShort(XRange);  /* rgnBBox.right  */

#ifdef EIGHTBIT
    WriteMSBShort(PICTpackbitsrect);
#else
    WriteMSBShort(PICTdirectbitsrect);
    WriteMSBShort(0x0000);  /* baseAddr      */
    WriteMSBShort(0x00ff);
#endif
    i = (XRange*sizeof(Pixel)) | 0x8000;
    WriteMSBShort( i );     /* rowBytes      */
    WriteMSBShort(0);       /* bounds.top    */
    WriteMSBShort(0);       /* bounds.left   */
    WriteMSBShort(YRange);  /* bounds.bottom */
    WriteMSBShort(XRange);  /* bounds.right  */
    WriteMSBShort(0);       /* pmVersion     */
#ifdef EIGHTBIT
    WriteMSBShort(0);       /* packType      */
#else
    WriteMSBShort(4);       /* packType      */
#endif
    WriteMSBLong(0);        /* packSize      */
    WriteMSBLong(72);       /* hRes          */
    WriteMSBLong(72);       /* vRes          */

#ifdef EIGHTBIT
    WriteMSBShort(0);       /* pixelType     */
    WriteMSBShort(8);       /* pixelSize     */
    WriteMSBShort(1);       /* cmpCount      */
    WriteMSBShort(8);       /* cmpSize       */
#else
    WriteMSBShort(16);      /* pixelType     */
    WriteMSBShort(32);      /* pixelSize     */
    WriteMSBShort(3);       /* cmpCount      */
    WriteMSBShort(8);       /* cmpSize       */
#endif

    WriteMSBLong(0);        /* planeBytes    */
    WriteMSBLong(0);        /* pmTable       */
    WriteMSBLong(0);        /* pmReserved    */

#ifdef EIGHTBIT
    WriteMSBLong(0);        /* ctSeed        */
    WriteMSBShort(0);       /* ctFlags       */
    WriteMSBShort(cols-1);  /* ctSize        */

    for( i=0; i<cols; i++ )
    {    WriteMSBShort(i);  /* value */
	 j=Node[i]; r=RLut[j]; g=GLut[j]; b=BLut[j];
	 WriteMSBShort( (r<<8)|r );  /* rgb.red */
	 WriteMSBShort( (g<<8)|g );  /* rgb.green */
	 WriteMSBShort( (b<<8)|b );  /* rgb.blue  */
    }
#endif

    WriteMSBShort(0);       /* srcRect.top    */
    WriteMSBShort(0);       /* srcRect.left   */
    WriteMSBShort(YRange);  /* srcRect.bottom */
    WriteMSBShort(XRange);  /* srcRect.right  */
    WriteMSBShort(0);       /* dstRect.top    */
    WriteMSBShort(0);       /* dstRect.left   */
    WriteMSBShort(YRange);  /* dstRect.bottom */
    WriteMSBShort(XRange);  /* dstRect.right  */
    WriteMSBShort(0);       /* mode (srcCopy) */

    size = WritePICTData();
    if( size & 0x01 ) putc(0,OutFile);
    WriteMSBShort(PICTendofpict);
    fclose(OutFile);
#ifdef APPLEMAC
    SetFileInfo(name,'ttxt','PICT',134);
#endif
    return( True );
}


int WriteIRISFile( name, encode )
    char *name;  int encode;
{
    return( True );
}


void InitialiseOutFile()
{
#if defined(IBMPC) || defined(APPLEMAC)
    /* Allocate Tables on FAR Heap */
    ABranch = (short __far*)_fmalloc(4096*sizeof(short));
    DBranch = (short __far*)_fmalloc(4096*sizeof(short));
    Hash = (short __far*)_fmalloc(256*sizeof(short));
    Node = (Byte __far*)_fmalloc(4096*sizeof(Byte));
#endif

    KinemageFlag = False;
    UseOutLine = False;
}

