/* command.h
 * RasMol2 Molecular Graphics
 * Roger Sayle, October 1994
 * Version 2.5
 */

#define MAXBUFFLEN   256
#define MAXLINELEN   256

#define FormatPDB        1
#define FormatXYZ        2
#define FormatAlchemy    3
#define FormatCharmm     4
#define FormatMol2       5
#define FormatCIF        6
#define FormatMDL        7


#ifdef COMMAND
int DataFileFormat;
char DataFileName[256];
char CurLine[MAXBUFFLEN];
int CurState,StateOption;
int CommandActive;
Long SelectCount;
int Interactive;

#else
extern int DataFileFormat;
extern char DataFileName[256];
extern char CurLine[MAXBUFFLEN];
extern int CurState,StateOption;
extern int CommandActive;
extern Long SelectCount;
extern int Interactive;

#ifdef FUNCPROTO
int ProcessCharacter( int );
int FetchFile( int, int, char* );
void LoadScriptFile( FILE*, char* );
void DisplaySelectCount();
void ResetCommandLine( int );
void InitialiseCommand();
int ExecuteIPCCommand( char __huge* );
int ExecuteCommand();
void ZapDatabase();

#else /* non-ANSI C compiler */
int ProcessCharacter();
int FetchFile();
void LoadScriptFile();
void DisplaySelectCount();
void ResetCommandLine();
void InitialiseCommand();
int ExecuteIPCCommand();
int ExecuteCommand();
void ZapDatabase();

#endif
#endif

