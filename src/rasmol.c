/* rasmol.c
 * RasMol2 Molecular Graphics
 * Roger Sayle, February 1994
 * Version 2.3
 */

#ifndef sun386
#include <stdlib.h>
#endif
#include <signal.h>
#include <stdio.h>
#include <math.h>

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
#include <sysv/unistd.h>
#else
#include <sys/termio.h>
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
static struct fd_set WaitSet;
static struct timeval TimeOut;
static struct termio OrigTerm;
static struct termio IntrTerm;
static int SocketNo,FileNo;
#endif

static char fnamebuf[128];
static char *FileNamePtr;
static int ProfCount;
static int LexState;
static int State;


#ifdef __STDC__
/* Function Prototype */
int HandleEvents( int );
#endif


void WriteChar( ch )
    char ch;
{   putc(ch,stderr);
}

void WriteString( ptr )
    char *ptr;
{   fputs(ptr,stderr);
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
            {   if( StdInChar == '[' )
                {   LexState = 2;  ch = 0;
                } else if( StdInChar != 0x1b )
                {   LexState = 0;  ch = StdInChar;
                } else ch = 0;

            } else /* LexState == 2 */
            {   LexState = 0;
                switch( StdInChar )
                {   case('A'): ch = 0x10;  break;
                    case('B'): ch = 0x0e;  break;
                    case('C'): ch = 0x06;  break;
                    case('D'): ch = 0x02;  break;
                    default:   ProcessCharacter('[');
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


void InitTerminal(socket)
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

        /* Queue an Asynchronous I/O Request */
        StdInStatus = sys$qio( 0, StdInChan, IO$_READVBLK|IO$M_NOECHO, 
                               StdInBlck, StdInASTEvent, 0, &StdInChar, 
                               1, 0, StdInMask, 0, 0);
    } else StdInFlag = False;
#endif

    setbuf(stdin,(char*)NULL);
}


void ResetTerminal()
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


static int FetchCharacter()
{
#ifdef TERMIOS
    register int status;
    register int width;

    if( SocketNo )
        do {
            if( !CommandActive )
                ResetCommandLine(0);
            if( HandleEvents(False) )
                RasMolExit();

            /* To avoid slow response time */
            if( !CommandActive )
                ResetCommandLine(0);

            FD_SET(FileNo,&WaitSet);
            FD_SET(SocketNo,&WaitSet);
            TimeOut.tv_usec = 0;
            TimeOut.tv_sec = 1;

            width = MaxFun(FileNo,SocketNo)+1;
            status = select( width, &WaitSet,
                     (struct fd_set*)NULL, (struct fd_set*)NULL, &TimeOut );
        } while( (status<1) || !FD_ISSET(FileNo,&WaitSet) );
#endif /* TERMIOS */

    if( !CommandActive )
        ResetCommandLine(0);
    return( getc(stdin) );
}

static int ReadCharacter()
{
    register int ch;

    if( LexState )
    {   ch = LexState;
        LexState = 0;
        return( ch );
    }

    ch = FetchCharacter();
    if( ch!=0x1b ) return( ch );
    ch = FetchCharacter();
    if( ch!='[' ) return( ch );

    switch( ch=FetchCharacter() )
    {   case('A'): return( 0x10 );
        case('B'): return( 0x0e );
        case('C'): return( 0x06 );
        case('D'): return( 0x02 );
    }
    LexState = ch;
    return('[');
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


void RasMolSignalExit()
{
    RasMolFatalExit("*** Quit ***");
}


static void LoadInitFile()
{
    register char *src,*dst;
    register FILE *initrc;
    register char *fname;

    fname = ".rasmolrc";
    initrc = fopen(fname,"r");
    if( !initrc && (src=(char*)getenv("HOME")) )
    {   dst = fnamebuf; 
        while( *src )
            *dst++ = *src++;
        *dst++ = '/';

        src = fname; fname = fnamebuf;
        while( (*dst++ = *src++) );
        initrc = fopen(fname,"r");
    }

    if( !initrc && (src=(char*)getenv("RASMOLPATH")) )
    {   dst = fnamebuf; 
        while( *src )
            *dst++ = *src++;
        *dst++ = '/';

        src = "rasmolrc"; fname = fnamebuf;
        while( (*dst++ = *src++) );
        initrc = fopen(fname,"r");
    }

    if( initrc )
        LoadScriptFile(initrc,fname);
}


static void HandleMenu( option )
     int option;
{
    register int mask;

    switch( State )
    {   /* Main Menu */
        case(1):  if( !option )
                  {   if( !Database )
                          ResetCommandLine(2);
                      break;
                  } else State = option+1;

                  switch(State)
                  {   case(1): NewMenu(MainMenuSize,MainMenu); break;
                      case(2): NewMenu(DispMenuSize,DispMenu); break;
                      case(3): NewMenu(ColrMenuSize,ColrMenu); break;
                      case(4): NewMenu(OptnMenuSize,OptnMenu); break;
                      case(5): NewMenu(SaveMenuSize,SaveMenu); break;
                      case(6): NewMenu(ConfMenuSize,ConfMenu); break;
                  }
                  break;

        /* Display Menu */
        case(2):  switch( option )
                  {   case(0):  DisableSpacefill();
                                EnableWireFrame(True,0);
                                SetRibbonStatus(False,0);
                                DisableBackBone();
                                break;

                      case(1):  DisableSpacefill();
                                DisableWireFrame();
                                SetRibbonStatus(False,0);
                                EnableBackBone(False,80);
                                break;

                      case(2):  DisableSpacefill();
                                EnableWireFrame(False,80);
                                SetRibbonStatus(False,0);
                                DisableBackBone();
                                break;

                      case(3):  SetVanWaalRadius();
                                DisableWireFrame();
                                SetRibbonStatus(False,0);
                                DisableBackBone();
                                break;

                      case(4):  SetRadiusValue(120);
                                EnableWireFrame(False,40);
                                SetRibbonStatus(False,0);
                                DisableBackBone();
                                break;

                      case(5):  DisableSpacefill();
                                DisableWireFrame();
                                SetRibbonStatus(True,0);
                                DisableBackBone();
                  }

                  if( option != 6 )
                      ReDrawFlag |= RFRefresh;
                  NewMenu(MainMenuSize,MainMenu);
                  State = 1;  break;

        /* Colour Menu */
        case(3):  if( option != 7 )
                  {   switch( option )
                      {   case(0): MonoColourAttrib(255,255,255);  break;
                          case(1): CPKColourAttrib();              break;
                          case(2): ShapelyColourAttrib();          break;
                          case(3): ScaleColourAttrib(GroupAttr);   break;
                          case(4): ScaleColourAttrib(ChainAttr);   break;
                          case(5): ScaleColourAttrib(TempAttr);    break;
                          case(6): StructColourAttrib();           break;
                      }
                      ReDrawFlag |= RFColour;
                  }

                  NewMenu(MainMenuSize,MainMenu);
                  State = 1;  break;

        /* Option Menu */
        case(4):  switch( option )
                  {   case(0):  UseSlabPlane = !UseSlabPlane;
                                if( UseSlabPlane )
                                    UseShadow = False;
                                break;

                      case(1):  mask = NormAtomFlag;
                                if( HetaGroups )
                                    mask |= HeteroFlag;
                                CommandActive = False;
                                Hydrogens = !Hydrogens;

                                if( Hydrogens )
                                {   fputs("\nHydrogens included!\n",stderr);
                                    SelectZone(mask|HydrogenFlag);
                                } else 
                                {   fputs("\nHydrogens removed!\n",stderr);
                                    RestrictZone(mask);
                                } break;

                      case(2):  mask = NormAtomFlag;
                                if( Hydrogens )
                                    mask |= HydrogenFlag;
                                CommandActive = False;
                                HetaGroups = !HetaGroups;

                                if( HetaGroups )
                                {   fputs("\nHETA atoms selected!\n",stderr);
                                    SelectZone(mask|HeteroFlag);
                                } else 
                                {   fputs("\nHETA atoms removed!\n",stderr);
                                    RestrictZone(mask);
                                } break;

                      case(3):  FakeSpecular = !FakeSpecular;
                                ReDrawFlag |= RFColour;
                                break;

                      case(4):  UseShadow = !UseShadow;
                                if( UseShadow )
                                {   ReviseInvMatrix();
                                    VoxelsClean = False;
                                    UseSlabPlane = False;
                                    ReAllocBuffers();
                                }
                                break;

                  }

                  if( option != 5 )
                      ReDrawFlag |= RFRefresh;
                  NewMenu(MainMenuSize,MainMenu);
                  State = 1;  break;

        /* Save Menu */
        case(5):  if( option==4 )
                  {   NewMenu(MainMenuSize,MainMenu);
                      State = 1;
                  } else
                  {   ResetCommandLine(3);
                      StateOption = option;
                  }
                  break;

        /* Quit Menu */
        case(6):  if( !option )
                      RasMolExit();
                  NewMenu(MainMenuSize,MainMenu);
                  State = 1;  break;
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


static int ProcessCommand()
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
                   NewMenu(MainMenuSize,MainMenu);
                   ResetCommandLine(1);
                   State = 1;
                   break;

         case(3):  /* Output Filename */
                   if( *CurLine ) switch( StateOption )
                   {   case(0):   WriteGIFFile(CurLine);            break;
                       case(1):   WriteEPSFFile(CurLine,True,True); break;
                       case(2):   WritePPMFile(CurLine,True);       break;
                       case(3):   WriteRastFile(CurLine,True);      break;
                   }
                   NewMenu(MainMenuSize,MainMenu);
                   ResetCommandLine(1);
                   State = 1;
                   break;
    }
    return( False );
}


int HandleEvents( wait )
    int wait;
{
    register int result;

    result = FetchEvent( wait );
    while( ReDrawFlag || result )
    {   if( result )
        {   if( result>0 )
            {   if( ProcessCharacter(result) )
                    if( ProcessCommand() )
                        return( True );
            } else if ( CurState==1 )
                HandleMenu( result+ButMax );
        } else
        {   if( ReDrawFlag&RFPoint )
            {   if( Database )
                    IdentifyAtom(PointX,PointY);
                ReDrawFlag &= ~RFPoint;
            }
            if( ReDrawFlag )
                RefreshScreen();
        }
        result = FetchEvent( False );
    }
    return( False );
}


#ifdef PROFILE
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
#endif /* PROFILE */

static void DisplayUsage()
{
    fputs("usage: rasmol [file]\n",stderr);
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
        if( *ptr == '-' )
        {
        } else
            if( !FileNamePtr )
            {   FileNamePtr = ptr;
            } else DisplayUsage();
    }
}


int main(argc,argv)
int argc; char *argv[];
{
    register char ch;

    FakeSpecular = False; SpecPower = 8;
    FileNamePtr = NULL;   ProfCount = 0;

    ProcessOptions(argc,argv);
    ReDrawFlag = 0;
#ifdef PROFILE
    ProfCount = 10;
#endif
    
    setbuf(stderr,(char *)NULL);
    if( (Interactive=OpenDisplay(CanvWidth,CanvHeight)) )
        NewMenu(MainMenuSize,MainMenu);
    InitTerminal(Interactive);
    signal(SIGINT,RasMolSignalExit);


    fputs("RasMol Molecular Renderer\n",stderr);
    fputs("Roger Sayle, February 1994\n",stderr);
    fputs("Version 2.3\n",              stderr);

#ifdef EIGHTBIT
    fputs("[8bit version]\n\n",  stderr);
#else
    fputs("[24bit version]\n\n", stderr);
#endif

    if( !Interactive )
        fputs("No suitable display detected!\n",stderr);

    InitialiseCommand();
    InitialiseTransform();
    InitialiseDatabase();
    InitialiseRenderer();
    InitialisePixUtils();
    InitialiseAbstree();
    InitialiseOutFile();

    LoadInitFile();

#ifdef PROFILE
    if( ProfCount )
    {   UseShadow = True;
        if( !FileNamePtr && !FetchFile(FormatPDB,True,FileNamePtr) )
            RasMolFatalExit("Profile Error: Invalid PDB file name!");

        SetVanWaalRadius();
        CPKColourAttrib();
        DefineColourMap();
        if( Interactive ) 
            FetchEvent(False);
        ProfileExecution();
    }
#endif

    if( FileNamePtr && FetchFile(FormatPDB,True,FileNamePtr) )
    {   ReDrawFlag |= RFRefresh | RFColour;
        if( InfoBondCount < 1 )
        {   EnableBackBone(False,80);
        } else EnableWireFrame(True,0);
        CPKColourAttrib();
        RefreshScreen();
    }

    State = 1;
    LexState = 0;
    ResetCommandLine(1);

#ifndef TERMIOS
    if( Interactive )
    {   while( !HandleEvents(True) )
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
