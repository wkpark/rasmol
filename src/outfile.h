/* outfile.h
 * RasMol2 Molecular Graphics
 * Roger Sayle, February 1994
 * Version 2.3
 */

#ifdef OUTFILE

#else

#ifdef __STDC__
int WriteVectPSFile( char* );
int WriteEPSFFile( char*, int, int );
int WriteRastFile( char*, int );
int WritePPMFile( char*, int );
int WriteGIFFile( char* );
int WriteBMPFile( char* );

int WriteMolScriptFile( char* );
int WriteScriptFile( char* );
void InitialiseOutFile();

#else /* non-ANSI C compiler */
int WriteVectPSFile();
int WriteEPSFFile();
int WriteRastFile();
int WritePPMFile();
int WriteGIFFile();
int WriteBMPFile();

int WriteMolScriptFile();
int WriteScriptFile();
void InitialiseOutFile();

#endif
#endif

