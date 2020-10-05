/* command.c
 * RasMol2 Molecular Graphics
 * Roger Sayle, October 1994
 * Version 2.5
 */
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

#include <string.h>
#include <ctype.h>
#include <stdio.h>

#if !defined(IBMPC) && !defined(VMS) && !defined(APPLEMAC)
#include <pwd.h>
#endif

#define COMMAND
#include "command.h"
#include "tokens.h"
#include "molecule.h"
#include "abstree.h"
#include "transfor.h"
#include "render.h"
#include "graphics.h"
#include "pixutils.h"
#include "outfile.h"
#include "script.h"


/* Macros for commonly used loops */
#define ForEachAtom  for(chain=Database->clist;chain;chain=chain->cnext) \
		     for(group=chain->glist;group;group=group->gnext)    \
		     for(ptr=group->alist;ptr;ptr=ptr->anext)
#define ForEachBond  for(bptr=Database->blist;bptr;bptr=bptr->bnext)


#define IsIdentChar(x)  ((isalnum(x))||((x)=='_')||((x)=='$'))


#ifndef VMS
#ifdef IBMPC
#define DirChar  '\\'
#else
#define DirChar  '/'
#endif
#endif


#define ErrSyntax        0
#define ErrBigNum        1
#define ErrBadOpt        2
#define ErrParam         3
#define ErrFilNam        4
#define ErrBadLoad       5
#define ErrNotNum        6
#define ErrNotSep        7
#define ErrNotBrac       8
#define ErrNoCol         9
#define ErrColour       10
#define ErrBadArg       11
#define ErrBadExpr      12
#define ErrParen        13
#define ErrScript       14
#define ErrFunc         15
#define ErrSetName      16
#define ErrBadSet       17
#define ErrInScrpt      18


static char *ErrorMsg[] = {
	"Invalid command syntax",            /* ErrSyntax  */
	"Parameter value too large",         /* ErrBigNum  */
	"Invalid parameter setting",         /* ErrBadOpt  */
	"Invalid parameter name",            /* ErrParam   */
	"Filename string expected",          /* ErrFilNam  */
	"Molecule database loaded",          /* ErrBadLoad */
	"Integer value expected",            /* ErrNotNum  */
	"Comma separator missing",           /* ErrNotSep  */
	"Close bracket ']' expected",        /* ErrNotBrac */
	"No colour specified",               /* ErrNoCol   */
	"Unknown or incorrect colour",       /* ErrColour  */
	"Invalid command argument",          /* ErrBadArg  */
	"Syntax error in expression",        /* ErrBadExpr */
	"Close parenthesis ')' expected",    /* ErrParen   */
	"Script command stack too deep",     /* ErrScript  */
	"Open parenthesis '(' expected",     /* ErrFunc    */
	"Invalid or missing atom set name",  /* ErrSetName */
	"Not enough memory to define set",   /* ErrBadSet  */
	"Disabled command in script file"    /* ErrInScrpt */
    };


typedef struct _HlpEntry {
		struct _HlpEntry __far *next;
		struct _HlpEntry __far *info;
		char __far *keyword;
		Long fpos;
		} HlpEntry;

#define HelpPool   16
static char *HelpFileName;
static char HelpFileBuf[80];
static HlpEntry __far *FreeInfo;
static HlpEntry __far *HelpInfo;


#define STACKSIZE  10
static char *NameStack[STACKSIZE];
static int LineStack[STACKSIZE];

#define HISTSIZE    4096
#define HISTMASK    4095
static char HistBuff[HISTSIZE];
static int MinHist,MaxHist;
static int CurHist;

static char *CurPrompt;
static int CurPos,MaxPos;
static int FileDepth;

static int TokenLength;
static Long TokenValue;
static char TokenIdent[128];
static char *TokenStart;
static char *TokenPtr;
static int CurToken;


static int RVal, GVal, BVal;

#ifdef FUNCPROTO
/* Forward Declarations */
int ExecuteCommand();
int ProcessLine();
int ExecuteIPCCommand( char __huge* );
#endif



static void UpdateLine()
{
    register int i;

    for( i=CurPos; i<MaxPos; i++ )
	WriteChar(CurLine[i]);
    WriteChar(' ');
    for( i=MaxPos+1; i>CurPos; i-- )
	WriteChar(0x08);
}

static void CopyHistory()
{
    register int i;

    for( i=CurPos; i>0; i-- )
	WriteChar(0x08);
    for( i=0; i<MaxPos; i++ )
	WriteChar(' ');
    WriteChar(0x0D);
    WriteString(CurPrompt);

    CurPos = 0;
    if( (i=CurHist) != MaxHist )
	while( HistBuff[i] )
	{   CurLine[CurPos++] = HistBuff[i];
	    WriteChar(HistBuff[i]);
	    i = (i+1) & HISTMASK;
	}
    CurLine[CurPos] = 0;
    MaxPos = CurPos;
}


int ProcessCharacter( ch )
    int ch;
{
    register int i;

    if( !ch ) return( False );
    if( (ch>=' ') && (ch<='~') )
    {   if( MaxPos<MAXLINELEN )
	{   for( i=MaxPos; i>CurPos; i-- )
		CurLine[i] = CurLine[i-1];
	    CurLine[CurPos++] = ch;
	    CurLine[++MaxPos] = 0;

	    WriteChar(ch);
	    if( CurPos<MaxPos )
		UpdateLine();
	} else 
	    WriteChar(0x07);
	
    } else
	switch( ch )
	{    case( 0x7f ):  /* DEL and ^H */
	     case( 0x08 ):  if( CurPos>0 )
			    {   for( i=CurPos; i<=MaxPos; i++ )
				    CurLine[i-1] = CurLine[i];
				CurPos--; MaxPos--;
				WriteChar(0x08);
				UpdateLine();
			    }
			    break;

	     case( 0x04 ):  if( CurPos<MaxPos ) /* ^D */
			    {   for( i=CurPos; i<MaxPos; i++ )
				    CurLine[i] = CurLine[i+1];
				MaxPos--; UpdateLine();
			    }
			    break;

	     case( 0x0d ):  /* ^M and ^J */
	     case( 0x0a ):  WriteChar('\n');
			    if( MaxPos )
				for( i=0; i<=MaxPos; i++ )
				{    HistBuff[MaxHist] = CurLine[i];
				     MaxHist=(MaxHist+1)&HISTMASK;
				     if( MaxHist==MinHist )
				     {   while( HistBuff[MinHist] )
					     MinHist=(MinHist+1)&HISTMASK;
					 MinHist=(MinHist+1)&HISTMASK;
				     }
				}
			    CommandActive = False;
			    return( True );
			    break;

	     case( 0x02 ):  if( CurPos>0 )  /* ^B */
			    {    WriteChar(0x08);
				 CurPos--;
			    }
			    break;

	     case( 0x06 ):  if( CurPos<MaxPos )  /* ^F */
				WriteChar(CurLine[CurPos++]);
			    break;

	     case( 0x01 ):  while( CurPos>0 )   /* ^A */
			    {    WriteChar(0x08);
				 CurPos--;
			    }
			    break;

	     case( 0x05 ):  while( CurPos<MaxPos )  /* ^E */
				WriteChar(CurLine[CurPos++]);
			    break;

	     case( 0x0c ):  WriteChar('\n');    /* ^L */
			    WriteString(CurPrompt);
			    for( i=0; i<MaxPos; i++ )
				WriteChar(CurLine[i]);
			    for( i=CurPos; i<MaxPos; i++ )
				WriteChar(0x08);
			    break;

	     case( 0x10 ):  if( CurHist != MinHist ) /* ^P */
			    {   CurHist -= 2;
				if( CurHist<0 )
				    CurHist += HISTSIZE;
				while( HistBuff[CurHist] )
				    CurHist=CurHist?CurHist-1:HISTMASK;
				CurHist = (CurHist+1)&HISTMASK;
				CopyHistory();
			    }
			    break;

	     case( 0x0e ):  if( CurHist != MaxHist ) /* ^N */
			    {   while( HistBuff[CurHist] )
				    CurHist = (CurHist+1)&HISTMASK;
				CurHist = (CurHist+1)&HISTMASK;
				CopyHistory();
			    }
			    break;
	}
    return( False );
}


void ResetCommandLine( state )
     int state;
{
    if( state )
    {   EnableMenus(state==1);
	switch( CurState=state )
	{   case(1):   CurPrompt="RasMol> ";            break;
	    case(2):   CurPrompt="PDB file name:";      break;
	    case(3):   CurPrompt="Image file name:";    break;
	    case(4):   CurPrompt="Molecule file name:"; break;
	}
    }

    if( CommandActive )
	WriteChar('\n');
    CommandActive = True;
    WriteString(CurPrompt);

    CurHist = MaxHist;
    CurPos = MaxPos = 0;
    CurLine[0] = 0;
}



static void CommandError( error )
    register char *error;
{
    register char *ptr;
    char buffer[40];

    if( TokenPtr )
    {   if( FileDepth > -1 )
	{   if( CommandActive )
		WriteChar('\n');
	    CommandActive=False;
	    
	    WriteString(CurLine);
	    WriteChar('\n');
	} else WriteString("        ");

	for( ptr=CurLine; ptr<TokenStart; ptr++ )
	    WriteChar(' ');
	WriteString("^\n");
    }

    if( FileDepth > -1 )
    {   if( LineStack[FileDepth] )
	{   if( NameStack[FileDepth] )
	    {   WriteChar('"');
		WriteString(NameStack[FileDepth]);
		WriteString("\",");
	    }
	    sprintf(buffer,"line %d: ",LineStack[FileDepth]);
	    WriteString(buffer);
	} else
	{   WriteString(NameStack[FileDepth]);
	    WriteString(": ");
	}
    }

    if( error )
    {   WriteString(error);
	WriteString("!\n");
    }
    CommandActive = False;
    CurToken = 0;
}


static char *ProcessFileName( name )
    char *name;
{
#if !defined(IBMPC) && !defined(VMS) && !defined(APPLEMAC)
    register struct passwd *entry;
    register char *temp;
    char username[64];
#endif
    register char *ptr;


    while( *name==' ' )
	name++;

#if defined(IBMPC) || defined(VMS)
    ptr = DataFileName;
    while( *name && (*name!=' ') )
    {   *ptr++ = ToUpper(*name);
	name++;
    }
#else
#ifdef APPLEMAC
    ptr = DataFileName;
    while( *name )
	*ptr++ = *name++;

#else /* UNIX */
    /* Perform filename globbing */
    if( *name=='~' )
    {   ptr = username;  name++;
	while( *name && (*name!=' ') && (*name!='/') )
	    *ptr++ = *name++;
	*ptr = '\0';

	ptr = DataFileName;
	if( *username )
	{   if( entry=getpwnam(username) )
	    {   temp = entry->pw_dir;
		endpwent();
	    } else /* Unknown user! */
	    {   temp = username;
		*ptr++ = '~';
	    }

	} else if( !(temp=(char*)getenv("HOME")) )
	    temp = ".";

	while( *temp )
	    *ptr++ = *temp++;
    } else ptr = DataFileName;

    while( *name && (*name!=' ') )
	*ptr++ = *name++;
#endif
#endif
    *ptr = '\0';
    return ptr;
}


/* Filename extensions! */
#define MaxFileExt  4
static char *FileExt[] = { "", ".Z", ".gz", ".z" };

static FILE *OpenDataFile( begin, end )
    char *begin, *end;
{
    register FILE *fp;
#if !defined(IBMPC) && !defined(APPLEMAC) && !defined(VMS)
    register char *src, *dst;
    register int i;
    
    for( i=0; i<MaxFileExt; i++ )
    {   dst = end; src = FileExt[i];
	while( *dst++ = *src++ );
	if( fp=fopen(begin,"r") )
	    break;
    }
#else
    fp = fopen(begin,"r");
#endif
    *end = '\0';
    return fp;
}


int FetchFile( format, info, name )
    int format, info;
    char *name;
{
#ifndef APPLEMAC
#if !defined(IBMPC) && !defined(VMS)
    register int ch, comp;
#endif /* UNIX */
    register char *src,*dst;
    char buffer[128];
#endif /* APPLEMAC */

    register int done;
    register FILE *fp;

    DataFileFormat = 0;
    name = ProcessFileName(name);
    fp = OpenDataFile(DataFileName,name);

#ifndef APPLEMAC
    /* Try using a default file path! */
    if( !fp )
    {   switch( format )
	{   case(FormatMol2):    src = (char*)getenv("RASMOLMOL2PATH"); break;
	    case(FormatAlchemy): src = (char*)getenv("RASMOLMOLPATH");  break;
	    case(FormatMDL):     src = (char*)getenv("RASMOLMDLPATH");  break;
	    case(FormatPDB):     src = (char*)getenv("RASMOLPDBPATH");  break;
	    case(FormatXYZ):     src = (char*)getenv("RASMOLXYZPATH");  break;
	    default:             src = NULL;
	}

	if( src && *src )
	{   dst = buffer;
	    while( *src ) *dst++ = *src++;
#ifndef VMS
	    if( *(dst-1) != DirChar )
		*dst++ = DirChar;
#endif

	    src = DataFileName;
	    while( *dst = *src++ )
#ifndef VMS
		if( *dst == DirChar )
		{   break;
		} else dst++;
#else
		dst++;
#endif

	    if( !(*dst) && (fp=OpenDataFile(buffer,dst)) )
	    {   dst = DataFileName;  src=buffer;
		while( *dst++ = *src++ );
	    }
	}
    }
#endif /* APPLEMAC */


    if( !fp )
    {   *name = '\0';
	if( CommandActive )
	    WriteChar('\n');
	WriteString("Error: File '");
	WriteString(DataFileName);
	WriteString("' not found!\n\n");
	CommandActive=False;
	return( False );
    }

#if !defined(IBMPC) && !defined(VMS) && !defined(APPLEMAC)
    done = getc(fp);
    if( done == 0x1f )
    {   done = getc(fp);
	fclose(fp);

	if( done == 0x9d )
	{   sprintf(buffer,"uncompress -c %s 2> /dev/null\n",DataFileName);
	} else if( done == 0x8b )
	{   sprintf(buffer,"gzip -cdq %s 2> /dev/null\n",DataFileName);
	} else /* bad magic number! */
	{   if( CommandActive )
		WriteChar('\n');
	    WriteString("Error: Unrecognised compression format!\n\n");
	    CommandActive=False;
	    return( False );
	}
   
	comp = True;
	if( !(fp=popen(buffer,"r")) )
	{   if( CommandActive )
		WriteChar('\n');
	    WriteString("Error: Unable to decompress file!\n\n");
	    CommandActive=False;
	    return( False );
	}
    } else /* Uncompressed! */
    {   ungetc(done,fp);
	comp = False;
    }
#endif

    done = False;
    switch( format )
    {   case(FormatAlchemy): done = LoadAlchemyMolecule(fp); break;
	case(FormatCharmm):  done = LoadCharmmMolecule(fp);  break;
	case(FormatMol2):    done = LoadMol2Molecule(fp);    break;
	case(FormatPDB):     done = LoadPDBMolecule(fp);     break;
	case(FormatMDL):     done = LoadMDLMolecule(fp);     break;
	case(FormatXYZ):     done = LoadXYZMolecule(fp);     break;
    }

#if !defined(IBMPC) && !defined(VMS) && !defined(APPLEMAC)
    if( comp )
    {   if( pclose(fp) )
	{   if( CommandActive )
		WriteChar('\n');
	    WriteString("Error: Unable to decompress file!\n\n");
	    CommandActive=False;
	    return(False);
	}
    } else fclose(fp);
#else /* !defined(UNIX) */
    fclose(fp);
#endif

    if( !done ) 
    {   return( False );
    } else if( !Database ) 
	return( True );

    if( info )
	DescribeMolecule();
    DataFileFormat = format;
    AdviseUpdate(AdvName);
    AdviseUpdate(AdvClass);
    AdviseUpdate(AdvIdent);

#if !defined(IBMPC) && !defined(APPLEMAC)        
    if( Interactive ) 
       FetchEvent(False);
#endif

    ReDrawFlag |= RFInitial;
    if( InfoBondCount<=MainAtomCount )
    {   if( MainAtomCount+HetaAtomCount > 255 )
        {   CreateMoleculeBonds(info,False);
        } else CreateMoleculeBonds(info,True);
    }
    InitialTransform();

    VoxelsClean = False;
    ApplyTransform();
    return( True );
}


void LoadScriptFile( fileptr, name )
    FILE *fileptr;
    char *name;
{
    register char *ptr;
    register int ch,len;

    if( fileptr )
    {   len = 1;
	for( ptr=name; *ptr; ptr++ )
	    len++;

	FileDepth++;
	ptr = (char*)malloc( len );
	NameStack[FileDepth] = ptr;
	while( *ptr++ = *name++ );
	LineStack[FileDepth] = 0;

	do {
	    len = 0;
	    ch = getc(fileptr);
	    while( (ch!='\n') && (ch!=EOF) )
	    {   if( len<MAXBUFFLEN )
		    CurLine[len++] = ch;
		ch = getc(fileptr);
	    }

	    LineStack[FileDepth]++;
	    if( len<MAXBUFFLEN )
	    {   CurLine[len] = '\0';
		if( ExecuteCommand() )
		    break;
	    } else CommandError("Script command line too long");
	} while( ch!=EOF );
	free(NameStack[FileDepth]);
	fclose( fileptr );
	FileDepth--;
    } else
    {   CommandError( (char*)NULL );
	WriteString("Cannot open script file '");
	WriteString(name);  WriteString("'\n");
    }
}

#ifdef FUNCPROTO
/* Function Prototypes */
static int PrefixString( char __far*, char __far* );
#endif


static int PrefixString( str1, str2 )
    register char __far *str1, __far *str2;
{
    while( *str1 == *str2++ )
	if( *str1++ == '\0' )
	    return( True );
    return( *str1 == '\0' );
}


static HlpEntry __far *EnterHelpInfo( text )
    register char *text;
{
    register HlpEntry __far * __far *tmp;
    register HlpEntry __far *ptr;
    register int res,len,i;
    register char ch;

    char keyword[32];

    ptr = (void __far*)0;
    while( *text && (*text!='\n') )
    {   while( *text && (*text!='\n') && (*text==' ') )
	    text++;

	len = 0;
	while( *text && (*text!='\n') && (*text!=' ') )
	    if( len<31 )
	    {   ch = *text++;
		keyword[len++] = ToUpper(ch);
	    } else text++;
	keyword[len]='\0';

	if( ptr )
	{   tmp = &ptr->info;
	    ptr = (void __far*)0;
	} else tmp = &HelpInfo;

	while( *tmp )
	{   res = _fstrcmp(keyword,(*tmp)->keyword);
	    if( res==0 ) /* Exact Match */
	    {   ptr = *tmp;
		break;
	    } else if( res<0 )
		break;
	    tmp = &(*tmp)->next;
	}

	if( !ptr )
	{   if( !FreeInfo )
	    {   ptr = (HlpEntry __far*)_fmalloc(HelpPool*sizeof(HlpEntry));
		if( !ptr ) 
		    RasMolFatalExit("Command Error: Insufficient memory!");
		for( i=1; i<HelpPool; i++ )
		{   ptr->next = FreeInfo;
		    FreeInfo = ptr++;
		}
	    } else
	    {   ptr = FreeInfo;
		FreeInfo = ptr->next;
	    }

	    ptr->keyword = (char __far*)_fmalloc(len+1);
	    for( i=0; i<=len; i++ )
		ptr->keyword[i] = keyword[i];

	    ptr->info = (void __far*)0;
	    ptr->next = *tmp;
	    ptr->fpos = 0;
	    *tmp = ptr;
	}
    }
    return( ptr );
}

static void InitHelpFile()
{
    register char *src,*dst;
    register HlpEntry __far *fix;
    register HlpEntry __far *ptr;
    register FILE *fp;
    register Long pos;

    char buffer[82];


    HelpFileName = "rasmol.hlp";
    fp=fopen(HelpFileName,"r");

    if( !fp && (src=(char*)getenv("RASMOLPATH")) )
    {   HelpFileName = dst = HelpFileBuf; 
	while( *src )
	    *dst++ = *src++;
#ifndef VMS
	if( (dst!=HelpFileBuf) && (*(dst-1)!=DirChar) )
	    *dst++ = DirChar;
#endif

	src = "rasmol.hlp"; 
	while( *dst++ = *src++ );
	fp = fopen(HelpFileName,"r");
    }

#ifdef RASMOLDIR
    if( !fp )
    {   src = RASMOLDIR;
	HelpFileName = dst = HelpFileBuf;
	while( *src )
	    *dst++ = *src++;
#ifndef VMS
	if( (dst!=HelpFileBuf) && (*(dst-1)!=DirChar) )
	    *dst++ = DirChar;
#endif

	src = "rasmol.hlp"; 
	while( *dst++ = *src++ );
	fp = fopen(HelpFileName,"r");
    }
#endif

    if( !fp )
    {   if( CommandActive )
	    WriteChar('\n');
	CommandActive = False;
	
	WriteString("Unable to find RasMol help file!\n");
	HelpFileName = NULL;
	return;
    }

    pos = 0;
    fgets(buffer,80,fp);
    while( !feof(fp) )
    {    fix = (void __far*)0;
	 while( *buffer=='?' )
	 {   if( ptr = EnterHelpInfo(buffer+1) )
	     {   ptr->info = fix;
		 fix = ptr;
	     }

	     pos = ftell(fp);
	     if( !fgets(buffer,80,fp) )
		 break;
	 }

	 while( fix )
	 {   ptr = fix->info;
	     fix->info = (void __far*)0;
	     fix->fpos = pos;
	     fix = ptr;
	 }

	 while( fgets(buffer,80,fp) )
	     if( *buffer=='?' )
		 break;
    }
    fclose(fp);
}

static void FindHelpInfo()
{
    register HlpEntry __far * __far *tmp;
    register HlpEntry __far *ptr;
    register int res,len;
    register Long pos;
    register FILE *fp;
    register char ch;

    char keyword[32];
    char buffer[82];

    while( *TokenPtr && (*TokenPtr==' ') )
	TokenPtr++;

    if( *TokenPtr )
    {   ptr = NULL;
	do {
	    len = 0;
	    while( *TokenPtr && (*TokenPtr!=' ') )
		if( len<31 )
		{   ch = *TokenPtr++;
		    keyword[len++] = ToUpper(ch);
		} else TokenPtr++;
	    keyword[len]='\0';

	    if( ptr )
	    {   tmp = &ptr->info;
		ptr = (void __far*)0;
	    } else tmp = &HelpInfo;

	    while( *tmp )
	    {   res = _fstrcmp(keyword,(*tmp)->keyword);
		if( res<0 )
		{   if( PrefixString(keyword,(*tmp)->keyword) )
		    {   ptr = *tmp;
			if( ptr->next && 
			    PrefixString(keyword,ptr->next->keyword) )
			{   if( CommandActive ) WriteChar('\n');
			    WriteString("Ambiguous help topic requested!\n");
			    CommandActive = False;
			    return;
			} else break;
		    } else break;
		} else if( res==0 ) 
		{   ptr = *tmp;
		    break;
		}
		tmp = &(*tmp)->next;
	    }

	    while( *TokenPtr && (*TokenPtr==' ') )
		TokenPtr++;
	} while( *TokenPtr && ptr );

	if( !ptr || !ptr->fpos )
	{   if( CommandActive )
		WriteChar('\n');
	    WriteString("No available help on requested topic!\n");
	    CommandActive=False;
	    return;
	} else pos=ptr->fpos;
    } else pos=0;


    if( !(fp=fopen(HelpFileName,"r")) )
	RasMolFatalExit("Command Error: Unable to reopen help file!");

    if( CommandActive )
	WriteChar('\n');
    CommandActive = False;

    fseek(fp,pos,0);
    while( fgets(buffer,80,fp) )
	if( *buffer!='?' )
	{   WriteString(buffer);
	} else break;
    fclose(fp);
}


static int LookUpKeyword()
{
    register int mid,res;
    register int lo, hi;

    if( TokenLength>MAXKEYLEN )
	return( IdentTok );

    lo = KeyLen[TokenLength-1];
    hi = KeyLen[TokenLength]-1;

    while( hi>=lo )
    {   mid = (hi+lo)>>1;
	res = _fstrcmp(TokenIdent,Keyword[mid].ident);
	if( !res ) return( Keyword[mid].token );

	if( res>0 )
	{      lo = mid+1;
	} else hi = mid-1;
    }
    return( IdentTok );
}


static int FetchToken()
{
    register char ch;

    CurToken = 0;
    while( True )
    {    ch = *TokenPtr++;
	 if( !ch || (ch=='#') ) 
	     return(0);
	 if( isspace(ch) )
	     continue;

	 TokenStart = TokenPtr-1;
	 if( isalpha(ch) )
	 {   TokenLength = 1;
	     *TokenIdent = ToUpper(ch);
	     while( IsIdentChar(*TokenPtr) && (TokenLength<32) )
	     {   ch = *TokenPtr++;
		 TokenIdent[TokenLength++] = ToUpper(ch);
	     }
	     if( TokenLength==32 )
	     {   CommandError("Identifier too long");
		 return(0);
	     } else TokenIdent[TokenLength] = '\0';
	     return( CurToken = LookUpKeyword() );

	 } else if( isdigit(ch) )
	 {   TokenValue = ch-'0';
	     while( isdigit(*TokenPtr) )
		 TokenValue = 10*TokenValue + (*TokenPtr++)-'0';
	     return( CurToken = NumberTok );

	 } else if( (ch=='\'') || (ch=='\"') || (ch=='`') )
	 {   TokenLength = 0;
	     while( *TokenPtr && (TokenLength<128) && (*TokenPtr!=ch) )
		 TokenIdent[TokenLength++] = *TokenPtr++;

	     if( ch != *TokenPtr )
	     {   if( *TokenPtr )
		 {   CommandError("String constant unterminated");
		 } else CommandError("String constant too long");
		 return( 0 );
	     } else TokenPtr++;

	     TokenIdent[TokenLength]='\0';
	     return( CurToken = StringTok );
	 } else if( ispunct(ch) )
	     return( CurToken = ch );
    }
}


static int NextIf( token, error )
    int token, error;
{
    if( FetchToken()!=token )
    {   CommandError(ErrorMsg[error]);
	return( True );
    } else return( False );
}


static void FetchFloat( value, scale )
    Long value;  int scale;
{
    register int count;
    register int mant;

    if( !value && !isdigit(*TokenPtr) )
    {   CommandError("Invalid floating point number");
	TokenValue = 0;
	return;
    }

    mant = 0;
    count = 1;
    while( isdigit(*TokenPtr) )
    {   if( count < scale )
	{   mant = 10*mant + (*TokenPtr-'0');
	    count *= 10;
	}
	TokenPtr++;
    }

    mant = (scale*mant)/count;
    TokenValue = value*scale + mant;
}


static int ParseColour()
{
    switch( CurToken )
    {   case(BlueTok):        RVal=0;   GVal=0;   BVal=255; break;
	case(BlackTok):       RVal=0;   GVal=0;   BVal=0;   break;
	case(CyanTok):        RVal=0;   GVal=255; BVal=255; break;
	case(GreenTok):       RVal=0;   GVal=255; BVal=0;   break;
	case(GreenblueTok):   RVal=46;  GVal=139; BVal=87;  break;
	case(MagentaTok):     RVal=255; GVal=0;   BVal=255; break;
	case(OrangeTok):      RVal=255; GVal=165; BVal=0;   break;
	case(PurpleTok):      RVal=160; GVal=32;  BVal=240; break;
	case(RedTok):         RVal=255; GVal=0;   BVal=0;   break;
	case(RedorangeTok):   RVal=255; GVal=69;  BVal=0;   break;
	case(VioletTok):      RVal=238; GVal=130; BVal=238; break;
	case(WhiteTok):       RVal=255; GVal=255; BVal=255; break; 
	case(YellowTok):      RVal=255; GVal=255; BVal=0;   break;

	case('['):    RVal = GVal = BVal = 0;

		      if( NextIf(NumberTok,ErrNotNum) ) { return(False);
		      } else if( TokenValue>255 )
		      {   CommandError(ErrorMsg[ErrBigNum]); return( False );
		      } else RVal = (int)TokenValue;

		      if( NextIf(',',ErrNotSep) ) return(False);
		      if( NextIf(NumberTok,ErrNotNum) ) { return(False);
		      } else if( TokenValue>255 )
		      {   CommandError(ErrorMsg[ErrBigNum]); return( False );
		      } else GVal = (int)TokenValue;

		      if( NextIf(',',ErrNotSep) ) return(False);
		      if( NextIf(NumberTok,ErrNotNum) ) { return(False);
		      } else if( TokenValue>255 )
		      {   CommandError(ErrorMsg[ErrBigNum]); return( False );
		      } else BVal = (int)TokenValue;

		      return( !NextIf(']',ErrNotBrac) );

	case(IdentTok): if( Interactive )
			return( LookUpColour(TokenIdent,&RVal,&GVal,&BVal) );
		      
	default:  return(False);
    }
    return( True );
}


void DisplaySelectCount()
{
    char buffer[40];

    if( FileDepth == -1 )
    {   if( CommandActive )
	   WriteChar('\n');
	CommandActive=False;
    
	if( SelectCount==0 )
	{   WriteString("No atoms selected!\n");
	} else if( SelectCount>1 )
	{   sprintf(buffer,"%ld atoms selected!\n",SelectCount);
	    WriteString(buffer);
	} else WriteString("1 atom selected!\n");
    }

    if( DisplayMode )
	ReDrawFlag |= RFRefresh;
    AdviseUpdate(AdvSelectCount);
}


static void SelectZoneExpr( expr )
    Expr *expr;
{
    register Bond __far *bptr;

    if( !Database )
	return;

    SelectCount = 0;
    for( QChain=Database->clist; QChain; QChain=QChain->cnext )
	for( QGroup=QChain->glist; QGroup; QGroup=QGroup->gnext )
	    for( QAtom=QGroup->alist; QAtom; QAtom=QAtom->anext )
		if( EvaluateExpr(expr) )
		{   QAtom->flag |= SelectFlag;
		    SelectCount++;
		} else QAtom->flag &= ~SelectFlag;
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


static void RestrictZoneExpr( expr )
    Expr *expr;
{
    register Bond __far *bptr;
    register int flag;

    if( !Database )
	return;

    DrawAtoms = False;   MaxAtomRadius = 0;
    DrawBonds = False;   MaxBondRadius = 0;

    SelectCount = 0;
    for( QChain=Database->clist; QChain; QChain=QChain->cnext )
	for( QGroup=QChain->glist; QGroup; QGroup=QGroup->gnext )
	    for( QAtom=QGroup->alist; QAtom; QAtom=QAtom->anext )
		if( EvaluateExpr(expr) )
		{   QAtom->flag |= SelectFlag;
		    SelectCount++;

		    if( QAtom->flag & SphereFlag )
		    {   DrawAtoms = True;
			if( QAtom->irad>MaxAtomRadius )
			    MaxAtomRadius = QAtom->irad;
		    }
		}  else QAtom->flag &= ~(SelectFlag|SphereFlag);
    DisplaySelectCount();

    ForEachBond
    {   /* Ignore ZoneBoth setting! */
	flag = bptr->dstatom->flag & bptr->srcatom->flag;
	if( flag & SelectFlag )
	{   bptr->flag |= SelectFlag;
	    if( bptr->flag & CylinderFlag )
	    {   DrawBonds = True;
		if( bptr->irad>MaxBondRadius )
		    MaxBondRadius = bptr->irad;
	    } else if( bptr->flag&WireFlag )
		DrawBonds = True;
	} else bptr->flag &= ~(SelectFlag|DrawBondFlag);
    }

    DetermineClipping();
    VoxelsClean = False;
    BucketFlag = False;
}


static void CentreZoneExpr( expr )
    Expr *expr;
{
    register Real x, y, z;
    register Long count;

    if( !Database )
	return;

    count = 0;
    x = y = z = 0.0;
    for( QChain=Database->clist; QChain; QChain=QChain->cnext )
	for( QGroup=QChain->glist; QGroup; QGroup=QGroup->gnext )
	    for( QAtom=QGroup->alist; QAtom; QAtom=QAtom->anext )
		if( EvaluateExpr(expr) )
		{   x += (Real)QAtom->xorg;
		    y += (Real)QAtom->yorg;
		    z += (Real)QAtom->zorg;
		    count++;
		}

    if( count )
    {   CenX = (Long)(x/count);
	CenY = (Long)(y/count);
	CenZ = (Long)(z/count);
    } else
    {   if( CommandActive ) WriteChar('\n');
	WriteString("No Atoms to Centre!\n");
	CommandActive = False;
    }
}


static Expr *ParseRange( neg )
    int neg;
{
    register Expr *tmp1,*tmp2;
    register char ch;

    tmp1 = AllocateNode();
    tmp1->type = OpLftProp|OpRgtVal;
    tmp1->rgt.val = neg? -(int)TokenValue : (int)TokenValue;
    tmp1->lft.val = PropResId;

    if( *TokenPtr == '-' )
    {   TokenPtr++;
	neg = (*TokenPtr=='-');
	if( neg ) TokenPtr++;
	FetchToken();

	if( CurToken != NumberTok )
	{   CommandError(ErrorMsg[ErrNotNum]);
	    DeAllocateExpr( tmp1 );
	    return( (Expr*)NULL );
	}

	tmp1->type |= OpMoreEq;
	tmp2 = AllocateNode();
	tmp2->rgt.ptr = tmp1;
	tmp2->type = OpAnd;

	tmp1 = AllocateNode();
	tmp1->type = OpLftProp|OpRgtVal|OpLessEq;
	tmp1->rgt.val = neg? -(int)TokenValue : (int)TokenValue;
	tmp1->lft.val = PropResId;
	tmp2->lft.ptr = tmp1;
	tmp1 = tmp2;
    } else tmp1->type |= OpEqual;

    if( *TokenPtr == ':' )
	TokenPtr++;

    ch = *TokenPtr;
    if( isalnum(ch) )
    {   ch = ToUpper(ch);
	TokenPtr++;

	tmp2 = AllocateNode();
	tmp2->type = OpAnd;
	tmp2->rgt.ptr = tmp1;

	tmp1 = AllocateNode();
	tmp1->type = OpEqual | OpLftProp | OpRgtVal;
	tmp1->lft.val = PropChain;               
	tmp1->rgt.val = ch;

	tmp2->lft.ptr = tmp1;
	tmp1 = tmp2;
    } else if( (ch=='?') || (ch=='%') || (ch=='*') )
	TokenPtr++;

    FetchToken();
    return( tmp1 );
}


static Expr *ParseExpression( level )
    int level;
{
    register Expr *tmp1,*tmp2;
    register int done, pred;
    register int neg;

    switch( level )
    {    case(0): /* Disjunctions */
		  tmp1 = ParseExpression(1);
		  while( (CurToken==OrTok) || (CurToken=='|') ||
			 (CurToken==',') )
		  {   if( CurToken=='|' )
		      {   if( FetchToken()=='|' )
			      FetchToken();
		      } else FetchToken();

		      tmp2 = AllocateNode();
		      tmp2->type = OpOr;
		      tmp2->lft.ptr = tmp1;
		      tmp2->rgt.ptr = NULL;
		      if( !(tmp1=ParseExpression(1)) )
		      {   DeAllocateExpr(tmp2);
			  return( tmp1 );
		      }
		      tmp2->rgt.ptr = tmp1;
		      tmp1 = tmp2;
		  }
		  return( tmp1 );

	 case(1): /* Conjunctions */
		  tmp1 = ParseExpression(2);
		  while( (CurToken==AndTok) || (CurToken=='&') )
		  {   if( CurToken=='&' )
		      {   if( FetchToken()=='&' )
			      FetchToken();
		      } else FetchToken();

		      tmp2 = AllocateNode();
		      tmp2->type = OpAnd;
		      tmp2->lft.ptr = tmp1;
		      tmp2->rgt.ptr = NULL;
		      if( !(tmp1=ParseExpression(2)) )
		      {   DeAllocateExpr(tmp2);
			  return( tmp1 );
		      }
		      tmp2->rgt.ptr = tmp1;
		      tmp1 = tmp2;
		  }
		  return( tmp1 );

	 case(2): /* Primitives */
		  if( IsPredTok(CurToken) || (CurToken==BackboneTok) )
		  {   switch( CurToken )
		      {   case(HelixTok):    if( InfoHelixCount<0 )
						 DetermineStructure();
					     pred = PredHelix;
					     break;
			  case(SheetTok):    if( InfoLadderCount<0 )
						 DetermineStructure();
					     pred = PredSheet;
					     break;
			  case(TurnTok):     if( InfoTurnCount<0 )
						 DetermineStructure();
					     pred = PredTurn;
					     break;
			  case(CystineTok):  if( InfoSSBondCount<0 )
						 FindDisulphideBridges();
					     pred = PredCystine;     
					     break;
			  case(BackboneTok): pred = PredMainChain;   break;
			  case(SelectedTok): pred = PropSelect;      break;
			  default:  pred = PredAbsChr(PredTokOrd(CurToken));
		      }

		      tmp1 = AllocateNode();
		      tmp1->type = OpConst|OpLftProp|OpRgtVal;
		      tmp1->lft.val = pred;
		      FetchToken();
		      return( tmp1 );

		  } else if( IsPropTok(CurToken) )
		  {   tmp1 = AllocateNode();
		      tmp1->type = OpLftProp|OpRgtVal;
		      switch( CurToken )
		      {   case(TemperatureTok): pred = PropTemp;    break;
			  case(RadiusTok):      pred = PropRad;     break;
			  case(AtomNoTok):      pred = PropIdent;   break;
			  case(ElemNoTok):      pred = PropElemNo;  break;
			  case(ResNoTok):       pred = PropResId;   break;
		      }
		      tmp1->lft.val = pred;

		      FetchToken();
		      if( CurToken=='=' )
		      {   tmp1->type |= OpEqual;
			  if( FetchToken()=='=' )
			      FetchToken();
		      } else if( CurToken=='<' )
		      {   FetchToken();
			  if( CurToken=='>' )
			  {   tmp1->type |= OpNotEq;
			      FetchToken();
			  } else if( CurToken=='=' )
			  {   tmp1->type |= OpLessEq;
			      FetchToken();
			  } else tmp1->type |= OpLess;
		      } else if( CurToken=='>' )
		      {   if( FetchToken()=='=' )
			  {   tmp1->type |= OpMoreEq;
			      FetchToken();
			  } else tmp1->type |= OpMore;
		      } else if( (CurToken=='!') || (CurToken=='/') )
		      {   if( NextIf('=',ErrBadExpr) )
			  {   DeAllocateExpr( tmp1 );
			      return( (Expr*)NULL );
			  } else tmp1->type |= OpNotEq;
			  FetchToken();
		      } else
		      {   CommandError(ErrorMsg[ErrBadExpr]);
			  DeAllocateExpr( tmp1 );
			  return( (Expr*)NULL );
		      }


		      if( CurToken == '-' )
		      {   FetchToken();
			  neg = True;
		      } else neg = False;

		      if( CurToken!=NumberTok )
		      {   CommandError(ErrorMsg[ErrNotNum]);
			  DeAllocateExpr( tmp1 );
			  return( (Expr*)NULL );
		      } 

		      if( neg )
		      {     tmp1->rgt.val = -(int)TokenValue; 
		      } else tmp1->rgt.val = (int)TokenValue;
		      FetchToken();
		      return( tmp1 );
		      
		  } else switch( CurToken )
		  {   case('('):    FetchToken();
				    if( !(tmp1=ParseExpression(0)) )
					return( (Expr*)NULL );

				    if( CurToken!=')' )
				    {   CommandError(ErrorMsg[ErrParen]);
					DeAllocateExpr( tmp1 );
					return( (Expr*)NULL );
				    }
				    FetchToken();
				    return(tmp1);

		      case('!'): case('~'):
		      case(NotTok): FetchToken();
				    if( !(tmp1=ParseExpression(2)) )
					return( (Expr*)NULL );

				    tmp2 = AllocateNode();
				    tmp2->type = OpNot | OpRgtVal;
				    tmp2->lft.ptr = tmp1;
				    return( tmp2 );

		      case('-'):    if( NextIf(NumberTok,ErrNotNum) )
					return( (Expr*)NULL );
				    return( ParseRange(True) );

		      case(NumberTok):
				    return( ParseRange(False) );

		      case(WithinTok):
				    if( NextIf('(',ErrFunc) )
					return( (Expr*)NULL );

				    FetchToken();
				    if( CurToken==NumberTok )
				    {   if( *TokenPtr=='.' )
					{   TokenPtr++;
					    FetchFloat(TokenValue,250);
					}
				    } else if( CurToken!='.' )
				    {   CommandError(ErrorMsg[ErrNotNum]);
					return( (Expr*)NULL );
				    } else FetchFloat(0,250);

				    if( TokenValue>10000 )
				    {   CommandError(ErrorMsg[ErrBigNum]);
					return( (Expr*)NULL );
				    } else pred = (int)TokenValue;
				    if( NextIf(',',ErrNotSep) )
					return( (Expr*)NULL );

				    FetchToken();
				    if( !(tmp1=ParseExpression(0)) )
					return( (Expr*)NULL );

				    if( CurToken!=')' )
				    {   CommandError(ErrorMsg[ErrParen]);
					DeAllocateExpr( tmp1 );
					return( (Expr*)NULL );
				    }

				    FetchToken();
				    if( !pred )
					return( tmp1 );

				    tmp2 = AllocateNode();
				    tmp2->type = OpWithin;
				    tmp2->lft.limit = (Long)pred*pred;
				    tmp2->rgt.set = BuildAtomSet(tmp1);
				    DeAllocateExpr(tmp1);
				    return( tmp2 );

		      default:      if( CurToken==IdentTok )
				    {   tmp1 = LookUpSetExpr(TokenIdent);
					if( !tmp1 ) 
					    tmp1 = LookUpElement(TokenIdent);

					if( tmp1 )
					{   FetchToken();
					    return(tmp1);
					}
				    }

				    TokenPtr = TokenStart;
				    done = ParsePrimitiveExpr(&TokenPtr);
				    FetchToken();

				    if( !done )
				    {   CommandError(ErrorMsg[ErrBadExpr]);
					DeAllocateExpr( QueryExpr );
					return( (Expr*)NULL );
				    } else return( QueryExpr );
		  }
    }
    return( (Expr*)NULL );
}

static void ExecuteSetCommand()
{
    register int option;

    switch( FetchToken() )
    {   case(SlabTok):
        case(SlabModeTok):
	    option = -1;
	    FetchToken();
	    if( CurToken==RejectTok )
	    {   option = SlabReject;
	    } else if( CurToken==HalfTok )
	    {   option = SlabHalf;
	    } else if( CurToken==HollowTok )
	    {   option = SlabHollow;
	    } else if( CurToken==SolidTok )
	    {   option = SlabClose;
	    } else if( CurToken==SectionTok )
		option = SlabSection;

	    if( option != -1 )
	    {   if( UseSlabPlane && (SlabMode!=option) )
		    ReDrawFlag |= RFRefresh;
		SlabMode = option;
	    } else CommandError(ErrorMsg[ErrBadOpt]);
	    break;

	case(ShadowTok):
	    FetchToken();
	    if( CurToken==TrueTok )
	    {   UseShadow = True;
		ReviseInvMatrix();
		VoxelsClean = False;
		UseSlabPlane = False;
		ReDrawFlag |= RFRefresh;
		ReAllocBuffers();
	    } else if( CurToken==FalseTok )
	    {   ReDrawFlag |= RFRefresh;
		UseShadow = False;
	    } else CommandError(ErrorMsg[ErrBadOpt]);
	    break;
				  
	case(SpecularTok):
	    FetchToken();
	    if( CurToken==TrueTok )
	    {   FakeSpecular = True;
		ReDrawFlag |= RFColour;
	    } else if( CurToken==FalseTok )
	    {   FakeSpecular = False;
		ReDrawFlag |= RFColour;
	    } else CommandError(ErrorMsg[ErrBadOpt]);
	    break;

	case(SpecPowerTok):
	    FetchToken();
	    if( !CurToken )
	    {   SpecPower = 8;
		ReDrawFlag |= RFColour;
	    } else if( CurToken==NumberTok )
	    {   if( TokenValue<=100 )
		{   ReDrawFlag |= RFColour;
		    SpecPower = (int)TokenValue;
		} else 
		    CommandError(ErrorMsg[ErrBigNum]);
	    } else CommandError(ErrorMsg[ErrNotNum]);
	    break;

	case(AmbientTok):
	    FetchToken();
	    if( !CurToken )
	    {   ReDrawFlag |= RFColour;
		Ambient = DefaultAmbient;
	    } else if( CurToken==NumberTok )
	    {   if( TokenValue<=100 )
		{   Ambient = TokenValue/100.0;
		    ReDrawFlag |= RFColour;
		} else
		    CommandError(ErrorMsg[ErrBigNum]); 
	    } else CommandError(ErrorMsg[ErrNotNum]);
	    break;

	case(HeteroTok):
	    FetchToken();
	    if( CurToken==TrueTok )
	    {   HetaGroups = True;
	    } else if( CurToken==FalseTok )
	    {   HetaGroups = False;
	    } else CommandError(ErrorMsg[ErrBadOpt]);
	    break;
				  
	case(HydrogenTok):
	    FetchToken();
	    if( CurToken==TrueTok )
	    {   Hydrogens = True;
	    } else if( CurToken==FalseTok )
	    {   Hydrogens = False;
	    } else CommandError(ErrorMsg[ErrBadOpt]);
	    break;
				  

	case(BackgroundTok):
	    FetchToken();
	    if( !CurToken )
	    {   CommandError(ErrorMsg[ErrNoCol]);
	    } else if( ParseColour() )
	    {   ReDrawFlag |= RFColour;
		BackR = RVal;
		BackG = GVal;
		BackB = BVal;
#ifndef IBMPC
		FBClear = False;
#endif
	    } else if( CurToken )
		CommandError(ErrorMsg[ErrColour]);
	    break;

	case(BondModeTok):
	    FetchToken();
	    if( !CurToken || (CurToken==AndTok) )
	    {   ZoneBoth = True;
	    } else if( CurToken==OrTok )
	    {   ZoneBoth = False;
	    } else CommandError(ErrorMsg[ErrBadOpt]);
	    break;
	    
	case(HBondTok):
	    FetchToken();
	    if( (CurToken==BackboneTok) || (CurToken==MainChainTok) )
	    {   ReDrawFlag |= RFRefresh;
		HBondMode = True;
	    } else if( !CurToken || (CurToken==SidechainTok) )
	    {   ReDrawFlag |= RFRefresh;
		HBondMode = False;
	    } else CommandError(ErrorMsg[ErrBadOpt]);
	    break;

	case(SSBondTok):
	    FetchToken();
	    if( (CurToken==BackboneTok) || (CurToken==MainChainTok) )
	    {   ReDrawFlag |= RFRefresh;
		SSBondMode = True;
	    } else if( !CurToken || (CurToken==SidechainTok) )
	    {   ReDrawFlag |= RFRefresh;
		SSBondMode = False;
	    } else CommandError(ErrorMsg[ErrBadOpt]);
	    break;

	case(HourGlassTok):
	    FetchToken();
	    if( CurToken==TrueTok )
	    {   UseHourGlass = True;
	    } else if( CurToken==FalseTok )
	    {   UseHourGlass = False;
	    } else CommandError(ErrorMsg[ErrBadOpt]);
	    break;

	case(StrandsTok):
	    FetchToken();
	    if( !CurToken )
	    {   ReDrawFlag |= RFRefresh;
		SplineCount = 5;
	    } else if( CurToken==NumberTok )
	    {   if( (TokenValue>0) && (TokenValue<=5) )
		{   SplineCount = (int)TokenValue;
		    ReDrawFlag |= RFRefresh;
		} else if( TokenValue==9 )
		{   ReDrawFlag |= RFRefresh;
		    SplineCount = 9;
		} else CommandError(ErrorMsg[ErrBadOpt]);
	    } else CommandError(ErrorMsg[ErrNotNum]);
	    break;

	case(MouseTok):
	    FetchToken();
	    if( !CurToken || (CurToken==RasMolTok) )
	    {   if( Interactive )
		    SetMouseMode( MMRasMol );
	    } else if( CurToken==InsightTok )
	    {   if( Interactive )
		    SetMouseMode( MMInsight );
	    } else if( CurToken==QuantaTok )
	    {   if( Interactive )
		    SetMouseMode( MMQuanta );
	    } else CommandError(ErrorMsg[ErrBadOpt]);
	    break;

	case(DisplayTok):
	    FetchToken();
	    if( !CurToken || (CurToken==NormalTok) )
	    {   ReDrawFlag |= RFRefresh | RFColour;
		DisplayMode = 0;
	    } else if( CurToken==SelectedTok )
	    {   ReDrawFlag |= RFRefresh | RFColour;
		DisplayMode = 1;
	    } else CommandError(ErrorMsg[ErrBadOpt]);
	    break;

	case(AxesTok):
	    FetchToken();
	    if( !CurToken || (CurToken==FalseTok) )
	    {   ReDrawFlag |= RFRefresh;
		DrawAxes = False;
	    } else if( CurToken == TrueTok )
	    {   ReDrawFlag |= RFRefresh;
		DrawAxes = True;
	    } else CommandError(ErrorMsg[ErrBadOpt]);
	    break;

	case(BoundBoxTok):
	    FetchToken();
	    if( !CurToken || (CurToken==FalseTok) )
	    {   ReDrawFlag |= RFRefresh;
		DrawBoundBox = False;
	    } else if( CurToken == TrueTok )
	    {   ReDrawFlag |= RFRefresh;
		DrawBoundBox = True;
	    } else CommandError(ErrorMsg[ErrBadOpt]);
	    break;

	case(UnitCellTok):
	    FetchToken();
	    if( !CurToken || (CurToken==FalseTok) )
	    {   ReDrawFlag |= RFRefresh;
		DrawUnitCell = False;
	    } else if( CurToken == TrueTok )
	    {   ReDrawFlag |= RFRefresh;
		DrawUnitCell = True;
	    } else CommandError(ErrorMsg[ErrBadOpt]);
	    break;

	case(VectPSTok):
	    FetchToken();
	    if( !CurToken || (CurToken==FalseTok) )
	    {   UseOutLine = False;
	    } else if( CurToken == TrueTok )
	    {   UseOutLine = True;
	    } else CommandError(ErrorMsg[ErrBadOpt]);
	    break;

	case(KinemageTok):
	    FetchToken();
	    if( !CurToken || (CurToken==FalseTok) )
	    {   KinemageFlag = False;
	    } else if( CurToken == TrueTok )
	    {   KinemageFlag = True;
	    } else CommandError(ErrorMsg[ErrBadOpt]);
	    break;

	case(MenusTok):
	    FetchToken();
	    if( !CurToken || (CurToken==TrueTok) )
	    {   EnableMenus(True);
	    } else if( CurToken == FalseTok )
	    {   EnableMenus(False);
	    } else CommandError(ErrorMsg[ErrBadOpt]);
	    break;

	case(RadiusTok):
	    FetchToken();
	    if( !CurToken )
	    {   ProbeRadius = SolventDots? 300 : 0;
	    } else if( CurToken==NumberTok )
	    {   if( *TokenPtr=='.' )
		{   TokenPtr++;
		    FetchFloat(TokenValue,250);
		}

		if( TokenValue>500 )
		{   CommandError(ErrorMsg[ErrBigNum]);
		} else ProbeRadius = (int)TokenValue;
	    } else if( CurToken=='.' )
	    {   FetchFloat(0,250);
		if( TokenValue>500 )
		{   CommandError(ErrorMsg[ErrBigNum]);
		} else ProbeRadius = (int)TokenValue;

	    } else CommandError(ErrorMsg[ErrNotNum]);
	    break;

	case(SolventTok):
	    FetchToken();
	    if( !CurToken || (CurToken==FalseTok) )
	    {   SolventDots = False;
		ProbeRadius = 0;
	    } else if( CurToken == TrueTok )
	    {   SolventDots = True;
		ProbeRadius = 300;
	    } else CommandError(ErrorMsg[ErrBadOpt]);
	    break;

	case(FontSizeTok):
	    FetchToken();
	    if( CurToken==NumberTok )
	    {   if( TokenValue<=32 )
		{   if( DrawLabels )
			ReDrawFlag |= RFRefresh;
		    SetFontSize((int)TokenValue);
		} else CommandError(ErrorMsg[ErrBigNum]);
	    } else if( !CurToken )
	    {   if( DrawLabels )
		    ReDrawFlag |= RFRefresh;
		SetFontSize(8);
	    } else CommandError(ErrorMsg[ErrBadOpt]);
	    break;

	default:
	    CommandError(ErrorMsg[ErrParam]);
    }
}


static void ExecuteColourCommand()
{
    register int flag;

    flag = 0;
    switch( FetchToken() )
    {   case(AtomTok):
	    FetchToken();
	default:
	    if( !CurToken )
	    {   CommandError(ErrorMsg[ErrNoCol]);
	    } else switch( CurToken )
	    {   case(CPKTok):         CPKColourAttrib(); 
				      ReDrawFlag |= RFColour; break;

		case(AminoTok):       AminoColourAttrib();
				      ReDrawFlag |= RFColour; break;

		case(ShapelyTok):     ShapelyColourAttrib();
				      ReDrawFlag |= RFColour; break;
		
		case(UserTok):        UserMaskAttrib(MaskColourFlag);
				      ReDrawFlag |= RFColour; break;

		case(GroupTok):       ScaleColourAttrib(GroupAttr);
				      ReDrawFlag |= RFColour; break;

		case(ChainTok):       ScaleColourAttrib(ChainAttr);
				      ReDrawFlag |= RFColour; break;

		case(ChargeTok):      ScaleColourAttrib(ChargeAttr);
				      ReDrawFlag |= RFColour; break;

		case(TemperatureTok): ScaleColourAttrib(TempAttr);
				      ReDrawFlag |= RFColour; break;

		case(StructureTok):   StructColourAttrib();
				      ReDrawFlag |= RFColour; break;

		default:  if( ParseColour() )
			  {   MonoColourAttrib(RVal,GVal,BVal);
			      ReDrawFlag |= RFColour;
			  } else CommandError(ErrorMsg[ErrColour]);
	    }
	    break;

	case(BondTok):    
	    FetchToken();
	    if( !CurToken )
	    {   CommandError(ErrorMsg[ErrNoCol]);
	    } else if( CurToken==NoneTok )
	    {   ColourBondNone();
		ReDrawFlag |= RFColour;
	    } else if( ParseColour() )
	    {   ColourBondAttrib(RVal,GVal,BVal);
		ReDrawFlag |= RFColour;
	    } else CommandError(ErrorMsg[ErrColour]);
	    break;

	case(TraceTok):
	case(BackboneTok):
	    FetchToken();
	    if( !CurToken )
	    {   CommandError(ErrorMsg[ErrNoCol]);
	    } else if( CurToken==NoneTok )
	    {   ColourBackNone();
		ReDrawFlag |= RFColour;
	    } else if( ParseColour() )
	    {   ColourBackAttrib(RVal,GVal,BVal);
		ReDrawFlag |= RFColour;
	    } else CommandError(ErrorMsg[ErrColour]);
	    break;

	case(SSBondTok):
	    FetchToken();
	    if( !CurToken )
	    {   CommandError(ErrorMsg[ErrNoCol]);
	    } else if( CurToken==NoneTok )
	    {   ReDrawFlag |= RFColour;
		ColourHBondNone( False );
	    } else if( ParseColour() )
	    {   ReDrawFlag |= RFColour;
		ColourHBondAttrib(False,RVal,GVal,BVal);
	    } else CommandError(ErrorMsg[ErrColour]);
	    break;

	case(HBondTok):
	    FetchToken();
	    if( !CurToken )
	    {   CommandError(ErrorMsg[ErrNoCol]);
	    } else if( CurToken==NoneTok )
	    {   ReDrawFlag |= RFColour;
		ColourHBondNone( True );
	    } else if( CurToken==TypeTok )
	    {   ReDrawFlag |= RFColour;
		ColourHBondType();
	    } else if( ParseColour() )
	    {   ReDrawFlag |= RFColour;
		ColourHBondAttrib(True,RVal,GVal,BVal);
	    } else CommandError(ErrorMsg[ErrColour]);
	    break;

	case(DotsTok):
	    FetchToken();
	    if( !CurToken )
	    {   CommandError(ErrorMsg[ErrNoCol]);
	    } else if( CurToken==PotentialTok )
	    {   ReDrawFlag |= RFColour;
		ColourDotsPotential();
	    } else if( ParseColour() )
	    {   ReDrawFlag |= RFColour;
		ColourDotsAttrib(RVal,GVal,BVal);
	    } else CommandError(ErrorMsg[ErrColour]);
	    break;

	case(AxesTok):
	case(BoundBoxTok):
	case(UnitCellTok):
	    FetchToken();
	    if( !CurToken )
	    {   CommandError(ErrorMsg[ErrNoCol]);
	    } else if( ParseColour() )
	    {   BoxR = RVal;  BoxG = GVal;  BoxB = BVal;
		ReDrawFlag |= RFColour;
	    } else CommandError(ErrorMsg[ErrColour]);
	    break;

	case(LabelTok):
	    FetchToken();
	    if( !CurToken )
	    {   CommandError(ErrorMsg[ErrNoCol]);
	    } else if( CurToken==NoneTok )
	    {   ReDrawFlag |= RFColour;
		UseLabelCol = False;
	    } else if( ParseColour() )
	    {   LabR = RVal;  LabG = GVal;  LabB = BVal;
		ReDrawFlag |= RFColour;
		UseLabelCol = True;
	    } else CommandError(ErrorMsg[ErrColour]);
	    break;

	case(RibbonTok):   flag = RibColBoth;     break;
	case(Ribbon1Tok):  flag = RibColInside;   break;
	case(Ribbon2Tok):  flag = RibColOutside;  break;
    }

    if( flag )
    {   FetchToken();
	if( !CurToken )
	{   CommandError(ErrorMsg[ErrNoCol]);
	} else if( CurToken==NoneTok )
	{   ReDrawFlag |= RFColour;
	    ColourRibbonNone(flag);
	} else if( ParseColour() )
	{   ReDrawFlag |= RFColour;
	    ColourRibbonAttrib(flag,RVal,GVal,BVal);
	} else CommandError(ErrorMsg[ErrColour]);
    }
}


static void ExecuteShowCommand()
{
    register Chain __far *chn;
    register Group __far *grp;
    register int chain,count;
    register Real temp;
    register char *str;
    char buffer[40];

    switch( FetchToken() )
    {   case(InfoTok):
		DescribeMolecule();
		break;

	case(SequenceTok):
		if( CommandActive )
		    WriteChar('\n');
		CommandActive = False;
		if( !Database )
		    return;

		for( chn=Database->clist; chn; chn=chn->cnext )
		{   chain = (InfoChainCount<2);  count = 0;
		    for( grp=chn->glist; grp; grp=grp->gnext )
			if( grp->alist && !(grp->alist->flag&HeteroFlag) )
			{   if( !chain )
			    {   WriteString("Chain ");
				WriteChar(chn->ident);
				WriteString(":\n");
				chain = True;
			    }

			    if( count == 10 )
			    {   WriteChar('\n');
				count = 1;
			    } else count++;

			    str = Residue[grp->refno];
			    WriteChar(str[0]);
			    WriteChar(str[1]);
			    WriteChar(str[2]);

			    sprintf(buffer,"%-3d ",grp->serno);
			    WriteString(buffer);
			}
		    WriteChar('\n');
		}

		WriteChar('\n');
		break;

	case(SymmetryTok):
		if( CommandActive )
		    WriteChar('\n');
		CommandActive = False;

		if( *InfoSpaceGroup )
		{   sprintf(buffer,"Space Group ...... %s\n",InfoSpaceGroup);
		    WriteString(buffer);

		    sprintf(buffer,"Unit cell A ...... %g\n",InfoCellA);
		    WriteString(buffer);
		    sprintf(buffer,"Unit cell B ...... %g\n",InfoCellB);
		    WriteString(buffer);
		    sprintf(buffer,"Unit cell C ...... %g\n",InfoCellC);
		    WriteString(buffer);

		    temp = Rad2Deg*InfoCellAlpha;
		    sprintf(buffer,"Unit cell alpha .. %g\n",temp);
		    WriteString(buffer);
		    temp = Rad2Deg*InfoCellBeta;
		    sprintf(buffer,"Unit cell beta ... %g\n",temp);
		    WriteString(buffer);
		    temp = Rad2Deg*InfoCellGamma;
		    sprintf(buffer,"Unit cell gamma .. %g\n",temp);
		    WriteString(buffer);

		} else WriteString("No crystal symmetry data!\n");
		WriteChar('\n');
		break;

	default:
	    CommandError(ErrorMsg[ErrBadArg]);
    }
}

void ZapDatabase()
{
    register int i;

    for( i=0; i<8; i++ )
	DialValue[i] = 0.0;
    SelectCount = 0;

    DestroyDatabase();
    ResetSymbolTable();
    ResetTransform();
    ResetRenderer();

    ZoneBoth = True;
    HetaGroups = True;    
    Hydrogens = True;

    ResetColourMap();
    DefineColourMap();
    ClearBuffers();
    ReDrawFlag = 0;

    if( Interactive )
    {   UpdateScrollBars();
	ClearImage();
    }
    AdviseUpdate(AdvName);
    AdviseUpdate(AdvClass);
    AdviseUpdate(AdvIdent);
}


static void WriteImageFile( name, type )
    char *name;  int type;
{
    if( !type )
#ifdef EIGHTBIT
	type = GIFTok;
#else
	type = PPMTok;
#endif


    switch( type )
    {   case(GIFTok):     WriteGIFFile(name);             break;
	case(BMPTok):     WriteBMPFile(name);             break;
	case(PPMTok):     WritePPMFile(name,True);        break;
	case(SUNTok):     WriteRastFile(name,False);      break;
	case(SUNRLETok):  WriteRastFile(name,True);       break;
	case(PICTTok):    WritePICTFile(name);            break;
	case(IRISTok):    WriteIRISFile(name);            break;
	case(EPSFTok):    WriteEPSFFile(name,True,True);  break;
	case(MonoPSTok):  WriteEPSFFile(name,False,True); break;
	case(VectPSTok):  WriteVectPSFile(name);          break;

	case(RasMolTok):
	case(ScriptTok):     WriteScriptFile(name);     break;
	case(POVRayTok):     WritePOVRayFile(name);     break;
	case(KinemageTok):   WriteKinemageFile(name);   break;
	case(MolScriptTok):  WriteMolScriptFile(name);  break;
    }
}


static void ExecuteConnect( flag )
    int flag;
{
    register Bond __far *bptr;

    if( Database )
    {   ForEachBond
	    if( bptr->col )
		Shade[Colour2Shade(bptr->col)].refcount--;
	ReDrawFlag |= RFRefresh|RFColour;
	CreateMoleculeBonds(True,flag);
        EnableWireFrame(True,0);
    }
}

#ifdef IBMPC
/* Avoid Optimizer Warning */
#pragma optimize("g",off)
#endif

int ExecuteCommand()
{
    register char *param;
    register int option;
    register int i,done;
    register Long temp;
    FILE *script;

    TokenPtr = CurLine;
    if( !FetchToken() )
    {   TokenPtr = NULL;
	return( False );
    }

    switch( CurToken )
    {   case(LoadTok):    if( !Database )
			  {   FetchToken();
			      option = FormatPDB;
			      if( !*TokenPtr || *TokenPtr==' ' )
			      {   if( CurToken==PDBTok )
				  {   FetchToken();  option = FormatPDB;
				  } else if( CurToken==MDLTok )
				  {   FetchToken();  option = FormatMDL;
				  } else if( CurToken==XYZTok )
				  {   FetchToken();  option = FormatXYZ;
				  } else if( CurToken==AlchemyTok )
				  {   FetchToken();  option = FormatAlchemy;
				  } else if( CurToken==CharmmTok )
				  {   FetchToken();  option = FormatCharmm;
				  } else if( CurToken==Mol2Tok )
				  {   FetchToken();  option = FormatMol2;
				  }
			      }

			      done = (FileDepth == -1);
			      if( !CurToken )
			      {   CommandError(ErrorMsg[ErrFilNam]);
				  break;
			      } else if( CurToken==StringTok )
			      {      FetchFile(option,done,TokenIdent);
			      } else FetchFile(option,done,TokenStart);
			      CurToken = 0;

			      if( Database )
			      {   ReDrawFlag |= RFRefresh | RFColour;
				  if( InfoBondCount < 1 )
				  {   EnableBackBone(False,80);
				  } else EnableWireFrame(True,0);
				  CPKColourAttrib();
			      }
			  } else CommandError(ErrorMsg[ErrBadLoad]);
			  break;

	case(SelectTok):  FetchToken();
			  if( !CurToken )
			  {   option = NormAtomFlag;
			      if( HetaGroups ) option |= HeteroFlag;
			      if( Hydrogens )  option |= HydrogenFlag;
			      SelectZone(option);
			  } else if( CurToken==AllTok )
			  {   SelectZone(AllAtomFlag);
			  } else if( CurToken==NoneTok )
			  {   SelectZone(0x00);
			  } else
			      if( QueryExpr=ParseExpression(0) )
			      {   if( !CurToken )
				  {   SelectZoneExpr(QueryExpr);
				  } else CommandError(ErrorMsg[ErrSyntax]);
				  DeAllocateExpr(QueryExpr);
			      }
			  break;

	case(RestrictTok):
			  FetchToken();
			  if( !CurToken )
			  {   option = NormAtomFlag;
			      if( HetaGroups ) option |= HeteroFlag;
			      if( Hydrogens )  option |= HydrogenFlag;
			      RestrictZone(option);
			      ReDrawFlag |= RFRefresh;
			  } else if( CurToken==AllTok )
			  {   RestrictZone(AllAtomFlag);
			      ReDrawFlag |= RFRefresh;
			  } else if( CurToken==NoneTok )
			  {   RestrictZone(0x00);
			      ReDrawFlag |= RFRefresh;
			  } else
			      if( QueryExpr=ParseExpression(0) )
			      {   if( !CurToken )
				  {   RestrictZoneExpr(QueryExpr);
				      ReDrawFlag |= RFRefresh;
				  } else CommandError(ErrorMsg[ErrSyntax]);
				  DeAllocateExpr(QueryExpr);
			      } 
			  break;


	case(ColourTok):  ExecuteColourCommand();
			  break;


	case(WireframeTok):
			  FetchToken();
			  if( CurToken==FalseTok )
			  {   ReDrawFlag |= RFRefresh;
			      DisableWireFrame();
			  } else if( (CurToken==TrueTok) || !CurToken )
			  {   ReDrawFlag |= RFRefresh;
			      EnableWireFrame(True,0);
			  } else if( CurToken==NumberTok )
			  {   if( *TokenPtr=='.' )
			      {   TokenPtr++;
				  FetchFloat(TokenValue,250);
			      }

			      if( TokenValue<500 )
			      {   EnableWireFrame(False,(int)TokenValue);
				  ReDrawFlag |= RFRefresh;
			      } else CommandError(ErrorMsg[ErrBigNum]);
			  } else if( CurToken=='.' )
			  {   FetchFloat(0,250);
			      if( TokenValue<500 )
			      {   EnableWireFrame(False,(int)TokenValue);
				  ReDrawFlag |= RFRefresh;
			      } else CommandError(ErrorMsg[ErrBigNum]);
			  } else CommandError(ErrorMsg[ErrBadArg]);
			  break;

	case(TraceTok):
	case(BackboneTok):
			  FetchToken();
			  if( CurToken==FalseTok )
			  {   ReDrawFlag |= RFRefresh;
			      DisableBackBone();
			  } else if( (CurToken==TrueTok) || !CurToken )
			  {   ReDrawFlag |= RFRefresh;
			      EnableBackBone(True,0);
			  } else if( CurToken==NumberTok )
			  {   if( *TokenPtr=='.' )
			      {   TokenPtr++;
				  FetchFloat(TokenValue,250);
			      }


			      if( TokenValue<500 )
			      {   EnableBackBone(False,(int)TokenValue);
				  ReDrawFlag |= RFRefresh;
			      } else CommandError(ErrorMsg[ErrBigNum]);
			  } else if( CurToken=='.' )
			  {   FetchFloat(0,250);
			      if( TokenValue<500 )
			      {   EnableBackBone(False,(int)TokenValue);
				  ReDrawFlag |= RFRefresh;
			      } else CommandError(ErrorMsg[ErrBigNum]);
			  } else CommandError(ErrorMsg[ErrBadArg]);
			  break;

	case(CPKTok):
	case(SpacefillTok):
			  FetchToken();
			  if( CurToken==FalseTok )
			  {   ReDrawFlag |= RFRefresh;
			      DisableSpacefill();
			  } else if( CurToken==NumberTok )
			  {   if( *TokenPtr=='.' )
			      {   TokenPtr++;
				  FetchFloat(TokenValue,250);
			      }

			      if( TokenValue<=500 )
			      {   SetRadiusValue(MaxFun((int)TokenValue,1));
				  ReDrawFlag |= RFRefresh;
			      } else CommandError(ErrorMsg[ErrBigNum]);
			  } else if( CurToken=='.' )
			  {   FetchFloat(0,250);
			      if( TokenValue<=500 )
			      {   SetRadiusValue(MaxFun((int)TokenValue,1));
				  ReDrawFlag |= RFRefresh;
			      } else CommandError(ErrorMsg[ErrBigNum]);
			  } else if( CurToken==UserTok )
			  {   UserMaskAttrib(MaskRadiusFlag);
			      ReDrawFlag |= RFRefresh;
			  } else if( CurToken==TemperatureTok )
			  {   ReDrawFlag |= RFRefresh;
			      SetRadiusTemperature();
			  } else if( (CurToken==TrueTok) || !CurToken )
			  {   ReDrawFlag |= RFRefresh;
			      SetVanWaalRadius();
			  } else CommandError(ErrorMsg[ErrBadArg]);
			  break;

	case(SSBondTok):  FetchToken();
			  if( CurToken==NumberTok )
			  {   if( *TokenPtr=='.' )
			      {   TokenPtr++;
				  FetchFloat(TokenValue,250);
			      }

			      if( TokenValue<=500 )
			      {   SetHBondStatus(False,True,(int)TokenValue);
				  ReDrawFlag |= RFRefresh;
			      } else CommandError(ErrorMsg[ErrBigNum]);
			  } else if( CurToken=='.' )
			  {   FetchFloat(0,250);
			      if( TokenValue<=500 )
			      {   SetHBondStatus(False,True,(int)TokenValue);
				  ReDrawFlag |= RFRefresh;
			      } else CommandError(ErrorMsg[ErrBigNum]);
			  } else if( CurToken==FalseTok )
			  {   ReDrawFlag |= RFRefresh;
			      SetHBondStatus(False,False,0);
			  } else if( (CurToken==TrueTok) || !CurToken )
			  {   ReDrawFlag |= RFRefresh;
			      SetHBondStatus(False,True,0);
			  } else CommandError(ErrorMsg[ErrBadArg]);
			  break;

	case(HBondTok):   FetchToken();
			  if( CurToken==NumberTok )
			  {   if( *TokenPtr=='.' )
			      {   TokenPtr++;
				  FetchFloat(TokenValue,250);
			      }

			      if( TokenValue<=500 )
			      {   SetHBondStatus(True,True,(int)TokenValue);
				  ReDrawFlag |= RFRefresh;
			      } else CommandError(ErrorMsg[ErrBigNum]);
			  } else if( CurToken=='.' )
			  {   FetchFloat(0,250);
			      if( TokenValue<=500 )
			      {   SetHBondStatus(True,True,(int)TokenValue);
				  ReDrawFlag |= RFRefresh;
			      } else CommandError(ErrorMsg[ErrBigNum]);
			  } else if( CurToken==FalseTok )
			  {   ReDrawFlag |= RFRefresh;
			      SetHBondStatus(True,False,0);
			  } else if( (CurToken==TrueTok) || !CurToken )
			  {   ReDrawFlag |= RFRefresh;
			      SetHBondStatus(True,True,0);
			  } else CommandError(ErrorMsg[ErrBadArg]);
			  break;

	case(RibbonTok):  FetchToken();
			  if( CurToken==NumberTok )
			  {   if( *TokenPtr=='.' )
			      {   TokenPtr++;
				  FetchFloat(TokenValue,250);
			      }

			      if( TokenValue<=1000 )
			      {   SetRibbonStatus(True,RibbonFlag,
						  (int)TokenValue);
				  ReDrawFlag |= RFRefresh;
			      } else CommandError(ErrorMsg[ErrBigNum]);
			  } else if( CurToken=='.' )
			  {   FetchFloat(0,250);
			      if( TokenValue<=1000 )
			      {   SetRibbonStatus(True,RibbonFlag,
						  (int)TokenValue);
				  ReDrawFlag |= RFRefresh;
			      } else CommandError(ErrorMsg[ErrBigNum]);
			  } else if( CurToken==FalseTok )
			  {   ReDrawFlag |= RFRefresh;
			      SetRibbonStatus(False,RibbonFlag,0);
			  } else if( (CurToken==TrueTok) || !CurToken )
			  {   ReDrawFlag |= RFRefresh;
			      SetRibbonStatus(True,RibbonFlag,0);
			  } else CommandError(ErrorMsg[ErrBadArg]);
			  break;

	case(StrandsTok): FetchToken();
			  if( CurToken==NumberTok )
			  {   if( *TokenPtr=='.' )
			      {   TokenPtr++;
				  FetchFloat(TokenValue,250);
			      }

			      if( TokenValue<=1000 )
			      {   SetRibbonStatus(True,StrandFlag,
						  (int)TokenValue);
				  ReDrawFlag |= RFRefresh;
			      } else CommandError(ErrorMsg[ErrBigNum]);
			  } else if( CurToken=='.' )
			  {   FetchFloat(0,250);
			      if( TokenValue<=1000 )
			      {   SetRibbonStatus(True,StrandFlag,
						  (int)TokenValue);
				  ReDrawFlag |= RFRefresh;
			      } else CommandError(ErrorMsg[ErrBigNum]);
			  } else if( CurToken==FalseTok )
			  {   ReDrawFlag |= RFRefresh;
			      SetRibbonStatus(False,StrandFlag,0);
			  } else if( (CurToken==TrueTok) || !CurToken )
			  {   ReDrawFlag |= RFRefresh;
			      SetRibbonStatus(True,StrandFlag,0);
			  } else CommandError(ErrorMsg[ErrBadArg]);
			  break;

	case(DotsTok):    FetchToken();
			  if( CurToken==NumberTok )
			  {   if( TokenValue<=1000 )
			      {   if( TokenValue )
				  {   CalculateSurface((int)TokenValue);
				  } else CalculateSurface(1);
				  ReDrawFlag |= RFRefresh;
			      } else CommandError(ErrorMsg[ErrBigNum]);
			  } else if( CurToken==FalseTok )
			  {   ReDrawFlag |= RFRefresh;
			      DeleteSurface();
			  } else if( (CurToken==TrueTok) || !CurToken )
			  {   ReDrawFlag |= RFRefresh;
			      CalculateSurface(100);
			  } else CommandError(ErrorMsg[ErrBadArg]);
			  break;

	case(SlabTok):    FetchToken();
			  if( (CurToken==NumberTok) || (CurToken=='.') )
			  {   if( CurToken==NumberTok )
			      {   if( *TokenPtr=='.' )
				  {   TokenPtr++;
				      FetchFloat(TokenValue,100);
				  } else TokenValue *= 100;
			      } else FetchFloat(0,100);

			      if( TokenValue<=10000 )
			      {   DialValue[7] = (TokenValue-5000)/5000.0;
				  /* UpdateScrollBars(); */
				  ReDrawFlag |= RFSlab;
				  UseSlabPlane = True;
				  UseShadow = False;
			      } else CommandError(ErrorMsg[ErrBigNum]);

			  } else if( CurToken==FalseTok )
			  {   if( UseSlabPlane )
			      {   ReDrawFlag |= RFRefresh;
				  UseSlabPlane = False;
			      }
			  } else if( !CurToken || (CurToken==TrueTok) )
			  {   if( !UseSlabPlane )
			      {   ReDrawFlag |= RFRefresh;
				  UseSlabPlane = True;
				  UseShadow = False;
			      }
			  } else CommandError(ErrorMsg[ErrSyntax]);
			  break;

	case(ZoomTok):    FetchToken();
			  if( (CurToken==NumberTok) || (CurToken=='.') )
			  {   if( CurToken==NumberTok )
			      {   if( *TokenPtr=='.' )
				  {   TokenPtr++;
				      FetchFloat(TokenValue,100);
				  } else TokenValue *= 100;
			      } else FetchFloat(0,100);

			      if( TokenValue<=10000 )
			      {   DialValue[3] = (TokenValue-10000)/10000.0;
				  ReDrawFlag |= RFZoom;
			      } else if( Database )
			      {   /* Magnification */
				  TokenValue -= 10000;
				  temp = (Long)(MaxZoom*10000);
				  if( TokenValue<=temp )
				  {   DialValue[3] = (Real)TokenValue/temp;
				      ReDrawFlag |= RFZoom;
				  } else CommandError(ErrorMsg[ErrBigNum]);
			      }
			  } else if( CurToken==TrueTok )
			  {   ReDrawFlag |= RFZoom;
			      DialValue[3] = 0.5;
			  } else if( !CurToken || (CurToken==FalseTok) )
			  {   ReDrawFlag |= RFZoom;
			      DialValue[3] = 0.0;
			  } else CommandError(ErrorMsg[ErrSyntax]);
			  /* UpdateScrollBars(); */
			  break;

	case(RotateTok):  FetchToken();
			  if( CurToken==XTok )
			  {   option = 0;
			  } else if( CurToken==YTok )
			  {   option = 1;
			  } else if( CurToken==ZTok )
			  {   option = 2;
			  } else
			  {   CommandError(ErrorMsg[ErrSyntax]);
			      break;
			  }

			  FetchToken();
			  if( done=(CurToken=='-') )
			      FetchToken();
#ifdef INVERT
			  if( option != 1 )
			      done = !done;
#endif
			  if( (CurToken==NumberTok) || (CurToken=='.') )
			  {   if( CurToken==NumberTok )
			      {   if( *TokenPtr=='.' )
				  {   TokenPtr++;
				      FetchFloat(TokenValue,100);
				  } else TokenValue *= 100;
			      } else FetchFloat(0,100);

			      if( TokenValue )
			      {   if( ReDrawFlag & RFRotate )
				      PrepareTransform();

				  ReDrawFlag |= (1<<option);
				  if( done ) TokenValue = -TokenValue;
				  DialValue[option] += TokenValue/18000.0;

				  while( DialValue[option]<-1.0 )
				      DialValue[option] += 2.0;
				  while( DialValue[option]>1.0 )
				      DialValue[option] -= 2.0;
				  if( Interactive )
				      UpdateScrollBars();
			      }
			  } else CommandError(ErrorMsg[ErrNotNum]);
			  break;

	case(TranslateTok):
			  FetchToken();
			  if( CurToken==XTok )
			  {   option = 4;
			  } else if( CurToken==YTok )
			  {   option = 5;
			  } else if( CurToken==ZTok )
			  {   option = 6;
			  } else
			  {   CommandError(ErrorMsg[ErrSyntax]);
			      break;
			  }

			  FetchToken();
			  if( done=(CurToken=='-') )
			      FetchToken();
#ifdef INVERT
			  if( option == 5 )
			      done = !done;
#endif

			  if( (CurToken==NumberTok) || (CurToken=='.') )
			  {   if( CurToken==NumberTok )
			      {   if( *TokenPtr=='.' )
				  {   TokenPtr++;
				      FetchFloat(TokenValue,100);
				  } else TokenValue *= 100;
			      } else FetchFloat(0,100);

			      if( TokenValue<=10000 )
			      {   ReDrawFlag |= (1<<option);
				  if( done ) TokenValue = -TokenValue;
				  DialValue[option] = TokenValue/10000.0;
				  /* UpdateScrollBars(); */
			      } else CommandError(ErrorMsg[ErrBigNum]);
			  } else CommandError(ErrorMsg[ErrNotNum]);
			  break;

	case(StereoTok):  FetchToken();
			  if( !CurToken || (CurToken==TrueTok) )
			  {   ReDrawFlag |= RFRefresh;
			      SetStereoMode(True);
			  } else if( CurToken==FalseTok )
			  {   ReDrawFlag |= RFRefresh;
			      SetStereoMode(False);
			  } else CommandError(ErrorMsg[ErrSyntax]);
			  break;

	case(CentreTok):  FetchToken();
			  if( !CurToken || (CurToken==AllTok) )
			  {   CenX = CenY = CenZ = 0;
			  } else
			      if( QueryExpr=ParseExpression(0) )
			      {   if( !CurToken )
				  {   CentreZoneExpr(QueryExpr);
				  } else CommandError(ErrorMsg[ErrSyntax]);
				  DeAllocateExpr(QueryExpr);
			      }
			  break;

	case(ResizeTok):  FetchToken();
			  break;

	case(ResetTok):   for( i=0; i<8; i++ )
			      DialValue[i] = 0.0;
			  ReDrawFlag |= RFDials;
			  ResetTransform();

			  /* ReDrawFlag |= RFRefresh|RFColour; */
			  /* DisplayMode = 0;                  */

			  if( Interactive )
			      UpdateScrollBars();
			  break;

	case('?'):
	case(HelpTok):    if( !HelpFileName )
			      InitHelpFile();
			  if( HelpInfo )
			      FindHelpInfo();
			  CurToken=0;
			  break;

	case(SetTok):     ExecuteSetCommand();
			  break;

	case(LabelTok):   FetchToken();
			  if( !CurToken || (CurToken==TrueTok) )
			  {   if( InfoChainCount>1 )
			      {   DefineLabels("%n%r:%c.%a");
			      } else if( MainGroupCount>1 )
			      {   DefineLabels("%n%r.%a");
			      } else DefineLabels("%e%i");
			  } else if( CurToken==FalseTok )
			  {   DeleteLabels();
			  } else if( CurToken==StringTok )
			  {   DefineLabels(TokenIdent);
			  } else DefineLabels(TokenStart);
			  ReDrawFlag |= RFRefresh;
			  break;

	case(EchoTok):    FetchToken();
			  if( CommandActive )
			      WriteChar('\n');
			  CommandActive = False;

			  if( CurToken==StringTok )
			  {   WriteString(TokenIdent);
			  } else if( CurToken )
			      WriteString(TokenStart);
			  WriteChar('\n');
			  CurToken = 0;
			  break;

	case(DefineTok):  FetchToken();
			  if( CurToken != IdentTok ) 
			  {   CommandError(ErrorMsg[ErrSetName]);
			      break;
			  }

			  if( (param = (char*)malloc(TokenLength+1)) )
			  {   for( i=0; i<=TokenLength; i++ )
				  param[i] = TokenIdent[i];

			      if( FetchToken() )
			      {   if( QueryExpr=ParseExpression(0) )
				  {   done = DefineSetExpr(param,QueryExpr);
				  } else done = True;
			      } else done = DefineSetExpr(param,(Expr*)NULL);
			  } else done = False;

			  if( !done )
			      CommandError(ErrorMsg[ErrBadSet]);
			  break;

	case(BackgroundTok):
			  FetchToken();
			  if( !CurToken )
			  {   CommandError(ErrorMsg[ErrNoCol]);
			  } else if( ParseColour() )
			  {   ReDrawFlag |= RFColour;
			      BackR = RVal;
			      BackG = GVal;
			      BackB = BVal;
#ifndef IBMPC
			      FBClear = False;
#endif
			  } else if( CurToken )
			      CommandError(ErrorMsg[ErrColour]);
			  break;

	case(WriteTok):
	case(SaveTok):    i = CurToken; /* Save keyword! */
			  if( (FileDepth!=-1) && LineStack[FileDepth] )
			  {   CommandError(ErrorMsg[ErrInScrpt]);
			      break;
			  }

			  option = FetchToken();
			  if( (option==RasMolTok) || (option==ScriptTok)
			      || IsMoleculeFormat(option)
			      || IsImageFormat(option) )
			  {   if( !*TokenPtr || *TokenPtr==' ' )
				  FetchToken();
			  } else if( i==SaveTok )
			  {   option = PDBTok;
			  } else option = 0;

			  if( !CurToken )
			  {   CommandError(ErrorMsg[ErrFilNam]);
			      break;
			  } else if( CurToken==StringTok )
			  {      ProcessFileName(TokenIdent);
			  } else ProcessFileName(TokenStart);
			  param = DataFileName;
			  CurToken = 0;

			  if( !IsMoleculeFormat(option) )
			  {   if( ReDrawFlag ) RefreshScreen();
			      WriteImageFile( param, option );

			  } else switch(option)
			  {   case(PDBTok):     SavePDBMolecule(param); break;
			      case(XYZTok):     SaveXYZMolecule(param); break;
			      case(CIFTok):     SaveCIFMolecule(param); break;
			      case(AlchemyTok): SaveAlchemyMolecule(param);
						break;
			  } break;

	case(SourceTok):
	case(ScriptTok):  FetchToken();
			  if( FileDepth<STACKSIZE )
			  {   if( !CurToken )
			      {   CommandError(ErrorMsg[ErrFilNam]);
				  break;
			      } else if( CurToken==StringTok )
			      {      ProcessFileName(TokenIdent);
			      } else ProcessFileName(TokenStart);
			      CurToken = 0;

			      script = fopen(DataFileName,"r");
			      LoadScriptFile(script,DataFileName);
			  } else CommandError(ErrorMsg[ErrScript]);
			  break;

	case(PrintTok):   if( !PrintImage() )
			  {   if( CommandActive )
				  WriteChar('\n');
			      WriteString("Warning: No suitable printer!\n");
			      CommandActive = False;
			  }
			  break;

	case(ClipboardTok):
			  if( !ClipboardImage() )
			  {   if( CommandActive )
				  WriteChar('\n');
			      WriteString("Unable to copy to clipboard!\n");
			      CommandActive = False;
			  }
			  break;

	case(RenumTok):   FetchToken();
			  if( CurToken )
			  {   if( done = (CurToken=='-') )
				  FetchToken();

			      if( CurToken==NumberTok )
			      {   if( done )
				  {     RenumberMolecule(-(int)TokenValue);
				  } else RenumberMolecule((int)TokenValue); 
			      } else CommandError(ErrorMsg[ErrNotNum]);
			  } else RenumberMolecule(1);
			  break;

	case(StructureTok):
			  DetermineStructure();
			  break;

	case(ConnectTok): FetchToken();
                          if( !CurToken )
                          {   if( MainAtomCount+HetaAtomCount > 255 )
                              {   ExecuteConnect(False);
                              } else ExecuteConnect(True);
                          } else if( CurToken==TrueTok )
                          {   ExecuteConnect(True);
                          } else if( CurToken==FalseTok )
                          {   ExecuteConnect(False);
			  } else CommandError(ErrorMsg[ErrSyntax]);
			  break;

	case(ShowTok):    ExecuteShowCommand();
			  break;

	case(ZapTok):     ZapDatabase();
			  break;

	case(QuitTok):    return( True );
	default:          CommandError("Unrecognised command");
			  break;
    }

    if( CurToken )
	if( FetchToken() )
	    CommandError("Warning: Ignoring rest of command");
    TokenPtr = NULL;
    return( False );
}

#ifdef IBMPC
/* Avoid Optimizer Warning! */
#pragma optimize("g",)
#endif


int ExecuteIPCCommand( ptr )
    char __huge *ptr;
{
    register char *src,*dst;
    register int len,depth;
    register int result;
    auto char buffer[256];

    result = 0;
    FileDepth = 0;
    *LineStack = 0;
#ifdef IBMPC
    *NameStack = "DDE Error";
#else
    *NameStack = "IPC Error";
#endif

    /* Save command line */
    src=CurLine;  dst=buffer;
    while( *dst++ = *src++ );
   
    while( *ptr && (*ptr==' ') )
	ptr++;

    if( *ptr=='[' )
    {   depth = 0;
	while( *ptr++ == '[' )
	{   dst = CurLine;
	    depth=0;  len=1;

	
	    do {
		if( *ptr==']' )
		{   if( !depth )
		    {   *dst = '\0';
			ptr++; break;
		    } else depth--;
		} else if( *ptr=='[' )
		    depth++;

		if( len<255 )
		{   *dst++ = *ptr;
		    len++;
		}
	    } while( *ptr++ );

	    if( len==255 )
	    {   if( CommandActive )
		    WriteChar('\n');
		WriteString("Warning: Remote command too long!\n");
		CommandActive = False;
	    } else if( ExecuteCommand() )
		result = 2;

	    while( *ptr && ((*ptr==' ')||(*ptr==';')) )
		ptr++;
	}
    } else if( *ptr )
    {   dst = CurLine;
	len = 0;

	while( True )
	{   if( len==255 )
	    {   if( CommandActive )
		    WriteChar('\n');
		WriteString("Warning: Remote command too long!\n");
		CommandActive = False;
		break;
	    }

	    if( !(*dst++ = *ptr++) )
	    {   if( ExecuteCommand() )
		    result = 2;
		break;
	    } else len++;
	}
    }

    FileDepth = -1;
    if( CommandActive )
    {   src=buffer; dst=CurLine;
	while( *dst++ = *src++ );
	if( !result ) result = 1;
    }

    return( result );
}


void InitialiseCommand()
{
    MaxHist = MinHist = 1;
    HistBuff[0] = 0;

    HelpFileName = NULL;
    FreeInfo = (void __far*)0;
    HelpInfo = (void __far*)0;

    CommandActive = False;
    SelectCount = 0;
    TokenPtr = NULL;
    FileDepth = -1;
}
