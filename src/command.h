/* command.h
 * RasMol2 Molecular Graphics
 * Roger Sayle, February 1994
 * Version 2.3
 */

#define MAXBUFFLEN   256
#define MAXLINELEN   256

#define FormatPDB        1
#define FormatXYZ        2
#define FormatAlchemy    3
#define FormatSybyl      4
#define FormatCIF        5
#define FormatMDL        6


#ifdef COMMAND
char DataFileName[80];
char CurLine[MAXBUFFLEN];
int CurState,StateOption;
int CommandActive;
Long SelectCount;
int Interactive;

#else
extern char DataFileName[80];
extern char CurLine[MAXBUFFLEN];
extern int CurState,StateOption;
extern int CommandActive;
extern Long SelectCount;
extern int Interactive;

#ifdef __STDC__
int ProcessCharacter( char );
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

