/***************************************************************************
 *                            RasMol 2.7.1.1                               *
 *                                                                         *
 *                                RasMol                                   *
 *                 Molecular Graphics Visualisation Tool                   *
 *                            20 January 2001                              *
 *                                                                         *
 *                   Based on RasMol 2.6 by Roger Sayle                    *
 * Biomolecular Structures Group, Glaxo Wellcome Research & Development,   *
 *                      Stevenage, Hertfordshire, UK                       *
 *         Version 2.6, August 1995, Version 2.6.4, December 1998          *
 *                   Copyright (C) Roger Sayle 1992-1999                   *
 *                                                                         *
 *                  and Based on Mods by Arne Mueller                      *
 *                      Version 2.6x1, May 1998                            *
 *                   Copyright (C) Arne Mueller 1998                       *
 *                                                                         *
 *       Version 2.7.0, 2.7.1, 2.7.1.1 Mods by Herbert J. Bernstein        *
 *           Bernstein + Sons, P.O. Box 177, Bellport, NY, USA             *
 *                      yaya@bernstein-plus-sons.com                       *
 *           2.7.0 March 1999, 2.7.1 June 1999, 2.7.1.1 Jan 2001           *
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

/* langsel.h
 */


typedef enum {
       English     =    0,
       Spanish     =    1
} language;

typedef enum {
       ErrSyntax   =    0,
       ErrBigNum   =    1,
       ErrBadOpt   =    2,
       ErrParam    =    3,
       ErrFilNam   =    4,
       ErrBadLoad  =    5,
       ErrNotNum   =    6,
       ErrNotSep   =    7,
       ErrNotBrac  =    8,
       ErrNoCol    =    9,
       ErrColour   =   10,
       ErrBadArg   =   11,
       ErrBadExpr  =   12,
       ErrParen    =   13,
       ErrScript   =   14,
       ErrFunc     =   15,
       ErrSetName  =   16,
       ErrBadSet   =   17,
       ErrInScrpt  =   18,
       ErrOutScrpt =   19,
       ErrBadMolDB =   20,

       StrErrFile  =   30,
       StrNotFnd   =   31,
       StrCFmt     =   32,
       StrDcmp     =   33,
       StrSLong    =   34,
       StrSMem     =   35,
       StrHFil     =   36,
       StrHTop     =   37,
       StrHNone    =   38,
       StrHROpn    =   39,
       StrCTerm    =   40,
       StrCLong    =   41,
       StrFNum     =   42,
       StrCent     =   43,
       StrCClip    =   44,
       StrDFile    =   45,
       StrNPrint   =   46,
       StrUCell    =   47,
       StrSGroup   =   48,
       StrSymm     =   49,
       StrUnrec    =   50,
       StrIgnore   =   51,
       StrRCLong   =   52,
       StrSFile    =   53,
       StrILong    =   54,
       StrMolNam   =   55,
       StrClass    =   56,
       StrSecSt    =   57,
       StrNoAsmt   =   58,
       StrPDBRec   =   59,
       StrCalc     =   60,
       StrDBCode   =   61,
       StrExpTec   =   62,
       StrNumChn   =   63,
       StrNumGrp   =   64,
       StrNumAtm   =   65,
       StrNumBnd   =   66,
       StrNumBrg   =   67,
       StrNumHbd   =   68,
       StrNumHel   =   69,
       StrNumStrnd =   70,
       StrNumTrn   =   71,
       StrMalloc   =   72,
       StrXSRes    =   73,
       StrXSAtyp   =   74,

       StrMOpen    =   80,
       StrMInfo    =   81,
       StrMSaveAs  =   82,
       StrMClose   =   83,
       StrMPrint   =   84,
       StrMPSetup  =   85,
       StrMExit    =   86,
       StrMEmpty   =   87,

       StrMWirefr  =   90,
       StrMBackbn  =   91,
       StrMSticks  =   92,
       StrMSpacefl =   93,
       StrMBallStk =   94,
       StrMRibbons =   95,
       StrMStrands =   96,
       StrMCartoon =   97,

       StrMMonochr =   98,
       StrMCPK     =   99,
       StrMShapely =  100,
       StrMGroup   =  101,
       StrMChain   =  102,
       StrMTemp    =  103,
       StrMStruct  =  104,
       StrMUser    =  105,
       StrMModel   =  106,
       StrMAlt     =  107,

       StrMSlab    =  108,
       StrMHydr    =  109,
       StrMHet     =  110,
       StrMSpec    =  111,
       StrMShad    =  112,
       StrMStereo  =  113,
       StrMLabel   =  114,
       
       StrMPOff    =  115,
       StrMPIdent  =  116,
       StrMPDist   =  117,
       StrMPMon    =  118,
       StrMPAng    =  119,
       StrMPTrsn   =  120,
       StrMPLabl   =  121,
       StrMPCent   =  122,
       StrMPCoord  =  123,
       StrMPBond   =  124,
       StrMRBond   =  125,
       StrMRMol    =  126,
       StrMRAll    =  127,

       StrMGIF     =  128,
       StrMPostscr =  129,
       StrMPPM     =  130,
       StrMIRGB    =  131,
       StrMSRast   =  132,
       StrMBMP     =  133,
       StrMPICT    =  134,

       StrMAbout   =  135,
       StrMUserM   =  136,
       
       StrMUndo    =  137,
       StrMCut     =  138,
       StrMCopy    =  139,
       StrMPaste   =  140,
#ifdef APPLEMAC
       StrMClear   =  141,
#else
       StrMDelete  =  141,
#endif
       StrMSelAll  =  142,

       StrMFile    =  143,
       StrMEdit    =  144,
       StrMDisplay =  145,
       StrMColour  =  146,
       StrMOpt     =  147,
       StrMExport  =  148,
#ifdef APPLEMAC
       StrMWindow  =  149,
       StrMHelp    =  150,
       StrMMainWin =  151,
       StrMCmndLin =  152,
#else
       StrMHelp    =  149,
#endif
       StrMSettings=  153,

       StrPrmtPDB  =  154,
       StrPrmtImg  =  155,
       StrPrmtMol  =  156


} strflag;

#define MaxStrFlag     160
typedef struct {
      char *    msg;
      strflag   msgno;
      language  lang;
      int       aux;
} langstr;


#ifdef LANGSEL

char * MsgStrs[MaxStrFlag];
int    MsgLens[MaxStrFlag];
int    MsgAuxl[MaxStrFlag];

#else

extern char * MsgStrs[MaxStrFlag];
extern int    MsgLens[MaxStrFlag];
extern int    MsgAuxl[MaxStrFlag];

#endif

void SwitchLang( language );
