/* rasmol.c
 * RasMol2 Molecular Graphics
 * Roger Sayle, October 1994
 * Version 2.5
 */

#ifndef sun386
#include <stdlib.h>
#endif
#include <signal.h>
#include <stdio.h>
#include <math.h>

#define RASMOL
#include "rasmol.h"
#include "graphics.h"
#include "molecule.h"
#include "transfor.h"
#include "command.h"
#include "abstree.h"
#include "render.h"
#include "pixutils.h"
#include "outfile.h"

#ifdef TERMIOS
#include <sys/types.h>
#include <sys/time.h>

#ifdef esv
#include <sysv/sys/termio.h>
#else
#ifdef __FreeBSD__
#include <sys/ioctl.h>
#include <sys/termios.h>
#define TCSETAW TIOCSETAW
#define TCGETA  TIOCGETA
#else
#include <sys/termio.h>
#endif /* __FreeBSD__ */
#endif /* esv */

#ifdef esv
#include <sysv/unistd.h>
#else
#include <unistd.h>
#endif

#if defined(_SEQUENT_) || defined(_AIX)
#include <sys/select.h>
#endif
#endif /* TERMIOS */

#ifdef VMS
#include <tt2def.h>
#include <iodef.h>
#endif


#define TwoPi    6.2832


#define MainMenuSize  6
static char *MainMenu[] = { "Load","Display","Colours",
                            "Options","Export","Quit" };
#define DispMenuSize  7
static char *DispMenu[] = { "Wireframe","Backbone","Sticks","Spacefill",
                            "Ball & Stick", "Ribbons", "Cancel" };

#define ColrMenuSize  8
static char *ColrMenu[] = { "Mono","CPK","Shapely","Group", "Chain",
                            "Temperature","Structure", "Cancel" };
#define OptnMenuSize  6
static char *OptnMenu[] = { "Slab Mode","Hydrogen","Het Atoms",
                            "Specular","Shadow","Cancel" };
#define SaveMenuSize  5
static char *SaveMenu[] = { "GIF","EPSF","PPM","Rast","Cancel" };

#define ConfMenuSize  2
static char *ConfMenu[] = { "Confirm", "Cancel" };


#ifdef VMS
static struct {
        unsigned short size; 
        unsigned short type;
        char *string;
        } StdInDesc = { 10, 0, "SYS$INPUT:" };

/* Character Terminator Mask! */
static int StdInMask[9] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };

static short StdInBlck[4];
static int StdInMode[3];
static int StdInOrig[3];
static short StdInChan;
static int StdInStatus;
static char StdInChar;
static int StdInFlag;
#endif

#ifdef TERMIOS
#ifdef __FreeBSD__
static struct termios OrigTerm;
static struct termios IntrTerm;
#else
static struct termio OrigTerm;
static struct termio IntrTerm;
#endif

static struct fd_set WaitSet;
static struct timeval TimeOut;
static int SocketNo,FileNo;
#endif

static char *FileNamePtr;
static char *ScriptNamePtr;
static int FileFormat;
static int ProfCount;
static int LexState;


/* Function Prototype */
#ifdef FUNCPROTO
static int HandleEvents( int );
#endif
int ProcessCommand();
void RasMolExit();


void WriteChar( ch )
    char ch;
{   putc(ch,stderr);
}

void WriteString( ptr )
    char *ptr;
{   fputs(ptr,stderr);
}


static void ResetTerminal()
{
#ifdef TERMIOS
    if( isatty(FileNo) )
        ioctl(FileNo, TCSETAW, &OrigTerm);
#endif

#ifdef VMS
    StdInFlag = False;
    if( StdInStatus & 0x01 )
        sys$cancel(StdInChan);

    sys$qiow( 0, StdInChan, IO$_SETMODE, 0, 0, 0,
              StdInOrig, 12, 0, 0, 0, 0 );
#endif
}


void RasMolExit()
{
    WriteChar('\n');
    if( CommandActive )
        WriteChar('\n');

    if( Interactive )
        CloseDisplay();
    ResetTerminal();
    exit(0);
}


void RasMolFatalExit( msg )
    char *msg;
{
    WriteChar('\n');
    WriteString(msg);
    WriteChar('\n');
    WriteChar(0x07);

    if( Interactive )
        CloseDisplay();
    ResetTerminal();
    exit(1);
}


#ifdef VMS
static int StdInASTEvent()
{
    register int ch;
    register int i;

    if( !StdInFlag )
        return( False );

    if( StdInBlck[0] & 0x01 )
    {   if( StdInBlck[1] )
        {   if( (StdInChar==0x03) || (StdInChar==0x1a) )
                RasMolFatalExit("*** Quit ***");

            if( LexState == 0 )
            {   if( StdInChar == 0x1b )
                {   LexState = 1;  ch = 0;
                } else ch = StdInChar;

            } else if( LexState == 1 )
            {   if( StdInChar=='[' )
                {   LexState = 2;  ch = 0;
                } else if( StdInChar=='O' )
                {   LexState = 3;  ch = 0;
                } else if( StdInChar != 0x1b )
                {   LexState = 0;  ch = StdInChar;
                } else ch = 0;

            } else /* LexState == 2 or 3 */
            {   LexState = 0;
                switch( StdInChar )
                {   case('A'): ch = 0x10;  break;
                    case('B'): ch = 0x0e;  break;
                    case('C'): ch = 0x06;  break;
                    case('D'): ch = 0x02;  break;
                    default:   if( LexState==2 )
                               {      ProcessCharacter('[');
                               } else ProcessCharacter('O');

                               if( StdInChar == 0x1b )
                               {   LexState = 1;  ch = 0;
                               } else ch = StdInChar;
                }
            }
        } else ch = '\n';

        if( ch && ProcessCharacter(ch) )
        {   if( ProcessCommand() )
                RasMolExit();

            if( !CommandActive )
                ResetCommandLine(0);
        }
    }

    /* Queue an Asynchronous I/O Request */
    StdInStatus = sys$qio( 0, StdInChan, IO$_READVBLK|IO$M_NOECHO, 
                           StdInBlck, StdInASTEvent, 0, &StdInChar, 
                           1, 0, StdInMask, 0, 0);
    return( True );
}
#endif


static void InitTerminal(socket)
    int socket;
{
#ifdef SIGTTIN
    signal(SIGTTIN,SIG_IGN);
#endif
#ifdef SIGTTOU
    signal(SIGTTOU,SIG_IGN);
#endif

#ifdef TERMIOS
    FileNo = fileno(stdin);
    SocketNo = socket;
    FD_ZERO(&WaitSet);

    if( isatty(FileNo) )
    {   ioctl(FileNo, TCGETA, &OrigTerm);
        IntrTerm = OrigTerm;

        IntrTerm.c_iflag |= IGNBRK|IGNPAR;
        IntrTerm.c_iflag &= ~(BRKINT|PARMRK|INPCK|IXON|IXOFF);
        IntrTerm.c_lflag &= ~(ICANON|ECHO|ECHOE|ECHOK|ECHONL|NOFLSH);
        IntrTerm.c_lflag |= ISIG;

        IntrTerm.c_cc[VMIN] = 1;
        IntrTerm.c_cc[VTIME] = 0;

#ifdef VSUSP /* Disable ^Z */
        IntrTerm.c_cc[VSUSP] = 0;
#endif

        ioctl(FileNo, TCSETAW, &IntrTerm);
    }
#endif /* TERMIOS */

#ifdef VMS
    /* Associate "SYS$INPUT" with channel StdInChan! */
    StdInStatus = sys$assign(&StdInDesc, &StdInChan, 0, 0, 0);
    if( StdInStatus & 0x01 )
    {   StdInFlag = True;

        /* Determine Original Terminal Mode */
        sys$qiow( 0, StdInChan, IO$_SENSEMODE, 0, 0, 0,
                  StdInOrig, 12, 0, 0, 0, 0 );

        StdInMode[0] = StdInOrig[0];
        StdInMode[1] = StdInOrig[1];
        StdInMode[2] = StdInOrig[2];

        /* Select "RAW" Terminal Mode */
        StdInMode[2] |= TT2$M_PASTHRU;

        /* Set Current Terminal Mode */
        sys$qiow( 0, StdInChan, IO$_SETMODE, 0, 0, 0,
                  StdInMode, 12, 0, 0, 0, 0 );

        if( socket )
        {   /* Queue an Asynchronous I/O Request */
            StdInStatus = sys$qio( 0, StdInChan, IO$_READVBLK|IO$M_NOECHO, 
                                   StdInBlck, StdInASTEvent, 0, &StdInChar, 
                                   1, 0, StdInMask, 0, 0);
        } else StdInStatus = False;
    } else StdInFlag = False;

#else /* !VMS */
    setbuf(stdin,(char*)NULL);
#endif
}


static int FetchCharacter()
{
#ifdef TERMIOS
    register int status;
    register int width;

    if( SocketNo )
        do {
            if( !CommandActive )
                ResetCommandLine(0);
            HandleEvents(False);

            /* To avoid slow response time */
            if( !CommandActive )
                ResetCommandLine(0);

            FD_SET(FileNo,&WaitSet);
            FD_SET(SocketNo,&WaitSet);
            TimeOut.tv_usec = 0;
            TimeOut.tv_sec = 1;

            width = MaxFun(FileNo,SocketNo)+1;
#ifdef __hpux
            status = select( width, (int*)&WaitSet, (int*)NULL, 
                                    (int*)NULL, &TimeOut );
#else
            status = select( width, &WaitSet, (struct fd_set*)NULL, 
                                    (struct fd_set*)NULL, &TimeOut );
#endif
        } while( (status<1) || !FD_ISSET(FileNo,&WaitSet) );
#endif /* TERMIOS */

    if( !CommandActive )
        ResetCommandLine(0);

#ifdef VMS
    sys$qiow( 0, StdInChan, IO$_READVBLK|IO$M_NOECHO, StdInBlck,
              0, 0, &StdInChar, 1, 0, StdInMask, 0, 0);
    return( StdInChar );
#else
    return( getc(stdin) );
#endif
}

static int ReadCharacter()
{
    register int tmp;
    register int ch;

    if( LexState )
    {   ch = LexState;
        LexState = 0;
        return( ch );
    }

    ch = FetchCharacter();
    if( ch!=0x1b ) return( ch );

    ch = FetchCharacter();
    if( (ch!='[') && (ch!='O') ) 
        return( ch );

    switch( tmp=FetchCharacter() )
    {   case('A'): return( 0x10 );
        case('B'): return( 0x0e );
        case('C'): return( 0x06 );
        case('D'): return( 0x02 );
    }
    LexState = tmp;
    return(ch);
}


void RasMolSignalExit( i )
    int i;
{
    RasMolFatalExit("*** Quit ***");
}


static void LoadInitFile()
{
    register char *src,*dst;
    register FILE *initrc;
    register char *fname;
    char fnamebuf[128];

#ifdef VMS
    fname = "RASMOL.INI";
#else
    fname = ".rasmolrc";
#endif

    initrc = fopen(fname,"r");
    if( !initrc && (src=(char*)getenv("HOME")) )
    {   dst = fnamebuf; 
        while( *src )
            *dst++ = *src++;
#ifndef VMS
        *dst++ = '/';
#endif

        src = fname; fname = fnamebuf;
        while( (*dst++ = *src++) );
        initrc = fopen(fname,"r");
    }

    if( !initrc && (src=(char*)getenv("RASMOLPATH")) )
    {   dst = fnamebuf; 
        while( *src )
            *dst++ = *src++;
#ifndef VMS
        *dst++ = '/';
#endif

        src = "rasmolrc"; fname = fnamebuf;
        while( (*dst++ = *src++) );
        initrc = fopen(fname,"r");
    }

    if( initrc )
        LoadScriptFile(initrc,fname);
}


static void HandleMenu( hand )
     int hand;
{
    register int menu;
    register int item;
    register int mask;

    menu = hand>>8;
    item = hand&0xff;
    switch( menu )
    {   case(0):  /* File Menu */
                  switch( item )
                  {   case(1):  /* Open */
                                if( !Database )
                                    ResetCommandLine(2);
                                break;

                      case(2):  /* Save As */
                                break;
                      case(3):  /* Close */
                                ZapDatabase();
                                break;

                      case(5):  /* Exit */
                                RasMolExit();
                                break;
                  } 
                  break;

        case(1):  /* Display Menu */
                  switch( item )
                  {   case(1):  /* Wireframe */
                                DisableSpacefill();
                                EnableWireFrame(True,0);
                                SetRibbonStatus(False,0,0);
                                DisableBackBone();
                                ReDrawFlag |= RFRefresh;
                                break;

                      case(2):  /* Backbone */
                                DisableSpacefill();
                                DisableWireFrame();
                                SetRibbonStatus(False,0,0);
                                EnableBackBone(False,80);
                                ReDrawFlag |= RFRefresh;
                                break;

                      case(3):  /* Sticks */
                                DisableSpacefill();
                                if( MainAtomCount<256 )
                                {   EnableWireFrame(False,40);
                                } else EnableWireFrame(False,80);
                                SetRibbonStatus(False,0,0);
                                ReDrawFlag |= RFRefresh;
                                DisableBackBone();
                                break;

                      case(4):  /* Spheres */
                                SetVanWaalRadius();
                                DisableWireFrame();
                                SetRibbonStatus(False,0,0);
                                DisableBackBone();
                                ReDrawFlag |= RFRefresh;
                                break;

                      case(5):  /* Ball & Stick */
                                SetRadiusValue(120);
                                EnableWireFrame(False,40);
                                SetRibbonStatus(False,0,0);
                                DisableBackBone();
                                ReDrawFlag |= RFRefresh;
                                break;

                      case(6):  /* Ribbons */
                                DisableSpacefill();
                                DisableWireFrame();
                                SetRibbonStatus(True,RibbonFlag,0);
                                DisableBackBone();
                                ReDrawFlag |= RFRefresh;
                                break;

                      case(7):  /* Strands */
                                DisableSpacefill();
                                DisableWireFrame();
                                SetRibbonStatus(True,StrandFlag,0);
                                DisableBackBone();
                                ReDrawFlag |= RFRefresh;
                  }
                  break;

        case(2):  /* Colours Menu */
                  switch( item )
                  {   case(1):  /* Monochrome */
                                MonoColourAttrib(255,255,255);
                                ReDrawFlag |= RFColour;  break;
                      case(2):  /* CPK */
                                CPKColourAttrib();
                                ReDrawFlag |= RFColour;  break;
                      case(3):  /* Shapely */
                                ShapelyColourAttrib();
                                ReDrawFlag |= RFColour;  break;
                      case(4):  /* Group */
                                ScaleColourAttrib( GroupAttr );
                                ReDrawFlag |= RFColour;  break;
                      case(5):  /* Chain */
                                ScaleColourAttrib( ChainAttr );
                                ReDrawFlag |= RFColour;  break;
                      case(6):  /* Temperature */
                                ScaleColourAttrib( TempAttr );
                                ReDrawFlag |= RFColour;  break;
                      case(7):  /* Structure */
                                StructColourAttrib();
                                ReDrawFlag |= RFColour;  break;
                      case(8):  /* User */
                                UserMaskAttrib(MaskColourFlag);
                                ReDrawFlag |= RFColour;  break;
                  }
                  break;

        case(3):  /* Option Menu */
                  switch( item )
                  {   case(1):  /* Slabbing */
                                ReDrawFlag |= RFRefresh;
                                UseSlabPlane = !UseSlabPlane;
                                if( UseSlabPlane )
                                    UseShadow = False;
                                break;

                      case(2):  /* Hydrogens */
                                mask = NormAtomFlag;
                                if( HetaGroups )
                                    mask |= HeteroFlag;
                                Hydrogens = !Hydrogens;
                                ReDrawFlag |= RFRefresh;
                                      
                                if( Hydrogens )
                                {   SelectZone(mask|HydrogenFlag);
                                } else RestrictZone(mask);
                                break;

                      case(3):  /* Hetero Atoms */
                                mask = NormAtomFlag;
                                if( Hydrogens )
                                    mask |= HydrogenFlag;
                                HetaGroups = !HetaGroups;
                                ReDrawFlag |= RFRefresh;
                                
                                if( HetaGroups )
                                {   SelectZone(mask|HeteroFlag);
                                } else RestrictZone(mask);
                                break;

                      case(4):  /* Specular */
                                FakeSpecular = !FakeSpecular;
                                ReDrawFlag |= RFColour;
                                break;

                      case(5):  /* Shadows */
                                ReDrawFlag |= RFRefresh;
                                UseShadow = !UseShadow;
                                if( UseShadow )
                                {   ReviseInvMatrix();
                                    VoxelsClean = False;
                                    UseSlabPlane = False;
                                    ReAllocBuffers();
                                }
                                break;
                  }
                  break;

        case(4):  /* Export Menu */
                  ResetCommandLine(3);
                  StateOption = item;
                  break;

        case(5):  /* Help Menu */
                  break;
    }
}


void RefreshScreen()
{
    if( !UseSlabPlane )
    {   ReDrawFlag &= ~(RFTransZ|RFSlab|RFPoint);
    } else ReDrawFlag &= ~(RFTransZ|RFPoint);

    if( ReDrawFlag )
    {   if( ReDrawFlag & RFReSize )
            ReSizeScreen();

        if( ReDrawFlag & RFColour )
        {   if( Interactive ) 
                ClearImage();
            DefineColourMap();
        }

        if( Database )
        {   if( Interactive )
                BeginWait();
            if( ReDrawFlag & RFApply ) 
                ApplyTransform();
            DrawFrame();
            if( Interactive )
            {   TransferImage();
                EndWait();
            }
        } else if( Interactive )
        {   ClearBuffers();
            TransferImage();
        }
        ReDrawFlag = 0;
    }
}


void AdviseUpdate( item )
    int item;
{
}


int ProcessCommand()
{
    switch(CurState)
    {    case(1):  /* RasMol Prompt */
                   return( ExecuteCommand() );

         case(2):  /* PDB Filename */
                   if( *CurLine && FetchFile(FormatPDB,False,CurLine) )
                   {   ReDrawFlag |= RFRefresh | RFColour;
                       if( InfoBondCount < 1 )
                       {   EnableBackBone(False,80);
                       } else EnableWireFrame(True,0);
                       CPKColourAttrib();
                   }
                   ResetCommandLine(1);
                   break;

         case(3):  /* Output Filename */
                   if( *CurLine ) switch( StateOption )
                   {   case(1):   WriteGIFFile(CurLine);            break;
                       case(2):   WriteEPSFFile(CurLine,True,True); break;
                       case(3):   WritePPMFile(CurLine,True);       break;
                       case(4):   WriteRastFile(CurLine,True);      break;
                       case(5):   WriteBMPFile(CurLine);            break;
                       case(6):   WritePICTFile(CurLine);           break;
                   }
                   ResetCommandLine(1);
                   break;
    }
    return( False );
}


static int HandleEvents( wait )
    int wait;
{
    register int result;

    result = FetchEvent( wait );
    while( ReDrawFlag || result )
    {   if( !result )
        {   if( ReDrawFlag&RFPoint )
            {   if( Database )
                    IdentifyAtom(PointX,PointY);
                ReDrawFlag &= ~RFPoint;
            }
            if( ReDrawFlag )
                RefreshScreen();
        } else HandleMenu( result );
        result = FetchEvent( False );
    }
    return( True );
}


static void ProfileExecution()
{
    register int start,stop;
    register Real delta;
    register int i;

    delta = TwoPi/ProfCount;

    start = time((long*)NULL);
    for( i=0; i<ProfCount; i++ )
    {   DrawFrame();
        if( Interactive )
            TransferImage();
        ReDrawFlag = RFRotateY;
        DialValue[1] += delta;
        /* UpdateScrollBars(); */
        ApplyTransform();
    }

    stop = time((long*)NULL);
    fprintf(stderr,"Execution of %d frames\n",ProfCount);
    fprintf(stderr,"Duration = %ld seconds\n",stop-start);
    RasMolExit();
}


static void DisplayUsage()
{
    fputs("usage: rasmol [-nodisplay] [-script scriptfile] ",stderr);
    fputs("[[-format] file]\n",stderr);
    fputs("    formats: -pdb -mdl -mol2 -xyz -alchemy -charmm\n",stderr);
    putc('\n',stderr);
    exit(1);
}

    
static void ProcessOptions(argc,argv)
    int argc;  char *argv[];
{
    register char *ptr;
    register int i;

    for( i=1; i<argc; i++ )
    {   ptr = argv[i];
#ifdef VMS
        if( (*ptr=='/') || (*ptr=='-') )
#else
        if( *ptr == '-' )
#endif
        {   ptr++;
            if( !strcmp(ptr,"pdb") )
            {   FileFormat = FormatPDB;
            } else if( !strcmp(ptr,"mdl") )
            {   FileFormat = FormatMDL;
            } else if( !strcmp(ptr,"charmm") )
            {   FileFormat = FormatCharmm;
            } else if( !strcmp(ptr,"alchemy") )
            {   FileFormat = FormatAlchemy;
            } else if( !strcmp(ptr,"mol2") || !strcmp(ptr,"sybyl") )
            {   FileFormat = FormatMol2;
            } else if( !strcmp(ptr,"xyz") )
            {   FileFormat = FormatXYZ;

            } else if( !strcmp(ptr,"nodisplay") )
            {   Interactive = False;
            } else if( !strcmp(ptr,"profile") )
            {   ProfCount = 10;
            } else if( !strcmp(ptr,"script") )
            {   if( i != argc-1 )
                {   ScriptNamePtr = argv[++i];
                } else DisplayUsage();
         

            } else DisplayUsage();
        } else
            if( !FileNamePtr )
            {   FileNamePtr = ptr;
            } else DisplayUsage();
    }
}


int main(argc,argv)
int argc; char *argv[];
{
    register FILE *fp;
    register int temp;
    register char ch;

    FileNamePtr = NULL;
    ScriptNamePtr = NULL;
    FileFormat = FormatPDB;
    Interactive = True;
    ProfCount = 0;

    ProcessOptions(argc,argv);
    ReDrawFlag = 0;
    
    temp = Interactive;
    setbuf(stderr,(char *)NULL);
    Interactive = OpenDisplay(InitialWide,InitialHigh);
    InitTerminal(Interactive);
    signal(SIGINT,RasMolSignalExit);

    fputs("RasMol Molecular Renderer\n",stderr);
    fputs("Roger Sayle, October 1994\n",stderr);
    fputs("Version 2.5.1\n",            stderr);

#ifdef EIGHTBIT
    fputs("[8bit version]\n\n",  stderr);
#else
    fputs("[24bit version]\n\n", stderr);
#endif

    if( !Interactive )
    {   if( !temp )
        {   fputs("Display window disabled!\n",stderr);
        } else fputs("No suitable display detected!\n",stderr);
    }

    InitialiseCommand();
    InitialiseTransform();
    InitialiseDatabase();
    InitialiseRenderer();
    InitialisePixUtils();
    InitialiseAbstree();
    InitialiseOutFile();

    if( ProfCount )
    {   UseShadow = True;
        if( !FileNamePtr && !FetchFile(FileFormat,True,FileNamePtr) )
            RasMolFatalExit("Profile Error: Invalid PDB file name!");

        SetVanWaalRadius();
        CPKColourAttrib();
        DefineColourMap();
        if( Interactive ) 
            FetchEvent(False);
        ProfileExecution();
    }

    if( FileNamePtr && FetchFile(FileFormat,True,FileNamePtr) )
    {   ReDrawFlag |= RFRefresh | RFColour;
        if( InfoBondCount < 1 )
        {   EnableBackBone(False,80);
        } else EnableWireFrame(True,0);
        CPKColourAttrib();
        RefreshScreen();
    }

    LexState = 0;
    ResetCommandLine(1);

    LoadInitFile();
    if( ScriptNamePtr )
    {   if( (fp=fopen(ScriptNamePtr,"r")) )
        {   LoadScriptFile(fp,ScriptNamePtr);
            fclose(fp);
        } else fprintf(stderr,"Error: File '%s' not found!\n",ScriptNamePtr);
    }

#ifndef TERMIOS
    if( Interactive )
    {   while( HandleEvents(True) )
            if( !CommandActive )
                ResetCommandLine(0);
    } else
#endif
    while( True )
    {   ch = ReadCharacter();
        if( ProcessCharacter(ch) )
            if( ProcessCommand() )
                break;
        RefreshScreen();
    }
    RasMolExit();
    return( 0 );
}
