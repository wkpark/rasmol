/* rasmol.h
 * RasMol2 Molecular Graphics
 * Roger Sayle, February 1994
 * Version 2.3
 */

#ifndef True
#define True  1
#define False 0
#endif

#define Real double
#define Byte unsigned char

#ifndef PI   /* Avoid Linux Warnings! */
#define PI   3.14159265358979323846
#endif

#ifdef __STDC__
#define Char signed char
#else
#define Char char
#endif

#ifdef _LONGLONG
#define Long int
#define Card unsigned int
#else
#define Long long
#define Card unsigned long
#endif


#define AbsFun(a)    (((a)<0)? -(a) : (a))
#define MinFun(a,b)  (((a)<(b))? (a) : (b) )
#define MaxFun(a,b)  (((a)>(b))? (a) : (b) )

#define EIGHTBIT
/* #define DIALBOX */
/* #define PROFILE */
#define TERMIOS
#define MITSHM
/* #define INVERT */
#define ISQRT
/* #define IBMPC */


#ifdef EIGHTBIT
#define Pixel    Byte
#else
#define Pixel    Long
#endif

#ifndef IBMPC
#define _fmalloc  malloc
#define _ffree    free
#define __huge
#define __far
#endif


#define ItemCount       6
#define AdvPickAtom	0
#define AdvPickNumber   1
#define AdvSelectCount  2
#define AdvName		3
#define AdvIdent        4
#define AdvClass	5


#ifndef RASMOL
#ifdef __STDC__
void WriteChar( int );
void WriteString( char* );
void RasMolFatalExit( char* );
void AdviseUpdate( int );
void RefreshScreen();
void RasMolExit();

#else /* non-ANSI C compiler */
void WriteChar();
void WriteString();
void RasMolFatalExit();
void RefreshScreen();
void AdviseUpdate();
void RasMolExit();

#endif
#endif
