/***************************************************************************
 *                             RasMol 2.7.2.1                              *
 *                                                                         *
 *                                 RasMol                                  *
 *                 Molecular Graphics Visualisation Tool                   *
 *                              14 April 2001                              *
 *                                                                         *
 *                   Based on RasMol 2.6 by Roger Sayle                    *
 * Biomolecular Structures Group, Glaxo Wellcome Research & Development,   *
 *                      Stevenage, Hertfordshire, UK                       *
 *         Version 2.6, August 1995, Version 2.6.4, December 1998          *
 *                   Copyright (C) Roger Sayle 1992-1999                   *
 *                                                                         *
 *                          and Based on Mods by                           *
 *Author             Version, Date             Copyright                   *
 *Arne Mueller       RasMol 2.6x1   May 98     (C) Arne Mueller 1998       *
 *Gary Grossman and  RasMol 2.5-ucb Nov 95     (C) UC Regents/ModularCHEM  *
 *Marco Molinaro     RasMol 2.6-ucb Nov 96         Consortium 1995, 1996   *
 *                                                                         *
 *Philippe Valadon   RasTop 1.3     Aug 00     (C) Philippe Valadon 2000   *
 *                                                                         *
 *Herbert J.         RasMol 2.7.0   Mar 99     (C) Herbert J. Bernstein    * 
 *Bernstein          RasMol 2.7.1   Jun 99         1998-2001               *
 *                   RasMol 2.7.1.1 Jan 01                                 *
 *                   RasMol 2.7.2   Aug 00                                 *
 *                   RasMol 2.7.2.1 Apr 01                                 *
 *                                                                         *
 *                    and Incorporating Translations by                    *
 *  Author                               Item                      Language*
 *  Isabel Serván Martínez,                                                *
 *  José Miguel Fernández Fernández      2.6   Manual              Spanish *
 *  José Miguel Fernández Fernández      2.7.1 Manual              Spanish *
 *  Fernando Gabriel Ranea               2.7.1 menus and messages  Spanish *
 *  Jean-Pierre Demailly                 2.7.1 menus and messages  French  *
 *  Giuseppe Martini, Giovanni Paolella, 2.7.1 menus and messages          *
 *  A. Davassi, M. Masullo, C. Liotto    2.7.1 help file           Italian *
 *                                                                         *
 *                             This Release by                             *
 * Herbert J. Bernstein, Bernstein + Sons, P.O. Box 177, Bellport, NY, USA *
 *                       yaya@bernstein-plus-sons.com                      *
 *               Copyright(C) Herbert J. Bernstein 1998-2001               *
 *              Copyright (C) Herbert J. Bernstein 1998-2001               *
 *                                                                         *
 * Please read the file NOTICE for important notices which apply to this   *
 * package. If you are not going to make changes to RasMol, you are not    *
 * only permitted to freely make copies and distribute them, you are       *
 * encouraged to do so, provided you do the following:                     *
 *   * 1. Either include the complete documentation, especially the file   *
 *     NOTICE, with what you distribute or provide a clear indication      *
 *     where people can get a copy of the documentation; and               *
 *   * 2. Please give credit where credit is due citing the version and    *
 *     original authors properly; and                                      *
 *   * 3. Please do not give anyone the impression that the original       *
 *     authors are providing a warranty of any kind.                       *
 *                                                                         *
 * If you would like to use major pieces of RasMol in some other program,  *
 * make modifications to RasMol, or in some other way make what a lawyer   *
 * would call a "derived work", you are not only permitted to do so, you   *
 * are encouraged to do so. In addition to the things we discussed above,  *
 * please do the following:                                                *
 *   * 4. Please explain in your documentation how what you did differs    *
 *     from this version of RasMol; and                                    *
 *   * 5. Please make your modified source code available.                 *
 *                                                                         *
 * This version of RasMol is not in the public domain, but it is given     *
 * freely to the community in the hopes of advancing science. If you make  *
 * changes, please make them in a responsible manner, and please offer us  *
 * the opportunity to include those changes in future versions of RasMol.  *
 ***************************************************************************/

/* vector.c
 $Log: vector.c,v $
 Revision 1.1  2001/01/31 02:13:45  yaya
 Initial revision

 Revision 1.3  2000/08/26 18:12:49  yaya
 Updates to header comments in all files

 Revision 1.2  2000/08/26 17:31:12  yaya
 Fix for world rot, remove refs to toolbar

 Revision 1.1  2000/08/09 01:18:23  yaya
 Rough cut with ucb

 */
/* Original file header from before 2.7.2 modifications:
 */
/**********************************************************************
  Copyright (c) 1995 UC Regents, ModularCHEM Consortium

  vector.c
  Vector Library for Rasmol
  
  Author:      Gary Grossman (garyg@cory.EECS.Berkeley.EDU)
  Last Update: November 14, 1995
 **********************************************************************/

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

#include <ctype.h>
#include <stdio.h>

#if !defined(IBMPC) && !defined(VMS) && !defined(APPLEMAC)
#include <pwd.h>
#endif

#include <math.h>

#include "command.h"
#include "tokens.h"
#include "molecule.h"
#include "abstree.h"
#include "transfor.h"
#include "vector.h"
#include "multiple.h"

#define PI 3.14159265358979323846

Real Length( Real u[3] )
{
    return sqrt(sqr(u[0])+sqr(u[1])+sqr(u[2]));
}

Real DotProduct( Real u[3], Real v[3])
{
    return u[0]*v[0] + u[1]*v[1] + u[2]*v[2];
}

void AbsoluteVector( Real v[3] )
{
    int i;
    for (i=0; i<3; i++) v[i] = fabs(v[i]);
}

void CrossProduct( Real u[3], Real v[3], Real n[3] )
{
    n[0] = u[1]*v[2] - u[2]*v[1];
    n[1] = u[2]*v[0] - u[0]*v[2];
    n[2] = u[0]*v[1] - u[1]*v[0];
} 

Real RadToDegrees(Real rad)
{
    return rad * (180.0 / PI);
}

/* Subtract two atoms and get a vector */
void SubtractAtoms( RAtom __far *a, RAtom __far *b, Real c[3] )
{
    c[0] = (double)(a->xorg - b->xorg)
           + (double)((a->xtrl - b->xtrl)/40.);
    c[1] = (double)(a->yorg - b->yorg)
           + (double)((a->ytrl - b->ytrl)/40.);
    c[2] = (double)(a->zorg - b->zorg)
           + (double)((a->ztrl - b->ztrl)/40.);
}

void NormalizeVector(Real u[3])
{
    Real n = Length(u);
    if (n) {
	u[0] /= n;
	u[1] /= n;
	u[2] /= n;
    }
}

void MultMatrix( Real A[4][4], Real B[4][4], Real C[4][4])
{
    int i, j, k;
    
    for (i=0; i<4; i++)
	for (j=0; j<4; j++) {
	    C[i][j] = 0;
	    for (k=0; k<4; k++)
		C[i][j] += A[i][k] * B[k][j];
	}
}

void IdentityMatrix(Real M[4][4])
{
    int i, j;
    for (i=0; i<4; i++) {
	for (j=0; j<4; j++) M[i][j] = 0;
	M[i][i] = 1;
    }
}

void PrintVector(Real v[3])
{
    int i; char buf[40];
    
    for (i=0; i<3; i++) {
	sprintf(buf, "%25.15g ", v[i]);
	WriteString(buf);
    }
    WriteString("\n");
}

void PrintMatrix(Real M[4][4])
{
    int i, j;
    char buf[20];
    
    WriteString("\n");
    for (i=0; i<4; i++) {
	for (j=0; j<4; j++) {
	    sprintf(buf, "%5.3f ", M[i][j]);
	    WriteString(buf);
	}
	WriteString("\n");
    }
}

void RotXMatrix(Real A[4][4], Real theta)
{
    Real cost=cos(theta), sint=sin(theta);
    
    IdentityMatrix(A);
    A[1][1] = cost;
    A[1][2] = sint;
    A[2][1] = -sint;
    A[2][2] = cost;
}

void RotYMatrix(Real A[4][4], Real theta)
{
    Real cost=cos(theta), sint=sin(theta);
    
    IdentityMatrix(A);
    A[0][0] = cost;
    A[0][2] = -sint;
    A[2][0] = sint;
    A[2][2] = cost;
}

void RotZMatrix(Real A[4][4], Real theta)
{
    Real cost=cos(theta), sint=sin(theta);
    
    IdentityMatrix(A);
    A[0][0] = cost;
    A[0][1] = sint;
    A[1][0] = -sint;
    A[1][1] = cost;
}

void RotateAxisMatrix(Real matrix[4][4], Real angle, Real rx, Real ry, Real rz)
{
    Real s=sin(angle), c=cos(angle);
    
    IdentityMatrix(matrix);
    
    matrix[0][0] = sqr(rx)+c*(1-sqr(rx));
    matrix[0][1] = rx*ry*(1-c)-rz*s;
    matrix[0][2] = rz*rx*(1-c)+ry*s;
    
    matrix[1][0] = rx*ry*(1-c)+rz*s;
    matrix[1][1] = sqr(ry)+c*(1-sqr(ry));
    matrix[1][2] = ry*rz*(1-c)-rx*s;
    
    matrix[2][0] = rz*rx*(1-c)-ry*s;
    matrix[2][1] = ry*rz*(1-c)+rx*s;
    matrix[2][2] = sqr(rz)+c*(1-sqr(rz));
}

