/* outfile.h
 * RasMol2 Molecular Graphics
 * Roger Sayle, October 1994
 * Version 2.5
 */

#ifdef OUTFILE
int UseOutLine;

#else
extern int UseOutLine;

#ifdef FUNCPROTO
int WriteVectPSFile( char* );
int WriteEPSFFile( char*, int, int );
int WriteRastFile( char*, int );
int WritePICTFile( char* );
int WriteIRISFile( char* );
int WritePPMFile( char*, int );
int WriteGIFFile( char* );
int WriteBMPFile( char* );
void InitialiseOutFile();

#else /* non-ANSI C compiler */
int WriteVectPSFile();
int WriteEPSFFile();
int WriteRastFile();
int WritePICTFile();
int WriteIRISFile();
int WritePPMFile();
int WriteGIFFile();
int WriteBMPFile();
void InitialiseOutFile();

#endif
#endif

