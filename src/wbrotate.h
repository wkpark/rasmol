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

/* wbrotate.h
 $Log: wbrotate.h,v $
 Revision 1.1  2001/01/31 02:13:45  yaya
 Initial revision

 Revision 1.7  2000/08/27 00:55:05  yaya
 create rotation bond database

 Revision 1.6  2000/08/26 18:13:03  yaya
 Updates to header comments in all files

 Revision 1.5  2000/08/26 17:31:26  yaya
 Fix for world rot, remove refs to toolbar

 Revision 1.4  2000/08/26 03:14:31  yaya
 Mods for mac compilations

 Revision 1.2  2000/08/21 21:08:10  yaya
 semi-final ucb mods

 Revision 1.1  2000/08/09 01:18:45  yaya
 Rough cut with ucb

 */

/* Original header for this routine:
 */
/**********************************************************************
  Copyright (c) 1995 UC Regents, ModularCHEM Consortium

  wbrotate.h
  World Rotate/Bond Rotate
  
  Author:      Gary Grossman (garyg@cory.EECS.Berkeley.EDU)
  Last Update: November 14, 1995
 **********************************************************************/

#ifndef WBROTATE_H
#define WBROTATE_H

/*========================*/
/* Bond Rotation Database */
/*========================*/

typedef struct _BondRot {
        struct _BondRot __far *brnext;     /* Next bond for rotation */
        RAtom __far *BSrcAtom;             /* First atom in the bond */
        RAtom __far *BDstAtom;             /* Last atom in the bond  */
        Real BRotValue;                    /* Angle of rotation      */    
    } BondRot;

void InitialiseWBRotate( void );
void WorldRotate( void );
void BondRotate( void );
void CreateBondAxis( Long, Long );
void SetBondAxis( RAtom __far *, RAtom __far * );
int RemoveBond(  Long , Long );
void ResetBondsSel( void );

#ifdef WBROTATE
Real WLastRX, WLastRY, WLastRZ;
Real WTransX, WTransY, WTransZ;
Real WLastTX, WLastTY, WLastTZ;
BondRot *BondSelected;
BondRot *BondsSelected;
RAtom __far *BSrcAtom;
RAtom __far *BDstAtom;
Real BAxis[3];
Real BRotValue, BLastRot;
Real WRotValue[3];
Real WLRotX[3],WLRotY[3],WLRotZ[3];
Real WIRotX[3],WIRotY[3],WIRotZ[3];
#else
extern Real WLastRX, WLastRY, WLastRZ;
extern Real WTransX, WTransY, WTransZ;
extern Real WLastTX, WLastTY, WLastTZ;
extern BondRot *BondSelected;
extern BondRot *BondsSelected;
#ifndef GRAPHICS
extern RAtom __far *BSrcAtom;
extern RAtom __far *BDstAtom;
#endif
extern Real BAxis[3];
extern Real BRotValue, BLastRot;
extern Real WRotValue[3];
extern Real WLRotX[3],WLRotY[3],WLRotZ[3];
extern Real WIRotX[3],WIRotY[3],WIRotZ[3];
#endif

#endif


