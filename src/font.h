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

/* font.h
 $Log: font.h,v $
 Revision 1.1  2001/01/31 02:13:45  yaya
 Initial revision

 Revision 1.3  2000/08/26 18:12:55  yaya
 Updates to header comments in all files

 Revision 1.2  2000/08/09 01:18:32  yaya
 Rough cut with ucb

 */

static char *VectFont[97] = {
    /*  32 ' ' */  "",
    /*  33 '!' */  "ChdhdicichCkdkewbwck",
    /*  34 '"' */  "CvbvbwcwcubtFvevewfwfuet",
    /*  35 '#' */  "AlmlAqmqDidtJijt",
    /*  36 '$' */  "AjcihijjkljnhodobparbtduiuktFgfw",
    /*  37 '%' */  "AgmwIhlhmimllmimhlhiihCqfqgrgufvcvbubrcq",
    /*  38 '&' */  "Kgihbsbtcvewfwgugsfqboamakbichegfghhkk",
    /*  39 ''' */  "Evdvdweweucs",
    /*  40 '(' */  "Egchbiakasbucvew",
    /*  41 ')' */  "Agchdiekesducvaw",
    /*  42 '*' */  "AplpCkjuCujk",
    /*  43 '+' */  "CokoGkgs",
    /*  44 ',' */  "Dgcgchdhdgdfce",
    /*  45 '-' */  "Coko",
    /*  46 '.' */  "Ggfgfhghgg",
    /*  47 '/' */  "Agmw",
    /*  48 '0' */  "EwcvbtbjchegigkhljltkviwewChkv",
    /*  49 '1' */  "DtgwggDgjg",
    /*  50 '2' */  "Kgagaibkhojpkrktjvhwdwbvat",
    /*  51 '3' */  "AtbvdwhwjvktkrjphoeoHojnklkjjhhgdgbhaj",
    /*  52 '4' */  "Klalhwhg",
    /*  53 '5' */  "Kwawaocphpjokmkjjhhgcgah",
    /*  54 '6' */  "Jvhwdwbvatajbhdghgjhkjkmjohpdpboan",
    /*  55 '7' */  "Awkwdg",
    /*  56 '8' */  "DwbvatarbpdohojpkrktjvhwdwDobnalajbhdghgjhkjkljnho",
    /*  57 '9' */  "Bhdghgjhkjktjvhwdwbvataqbodnhnjokq",
    /*  58 ':' */  "GififjgjgiGpfpfqgqgp",
    /*  59 ';' */  "DicicjdjdidhcgCpdpdqcqcp",
    /*  60 '<' */  "Lgbolw",
    /*  61 '=' */  "AsisAnjn",
    /*  62 '>' */  "Bglobw",
    /*  63 '?' */  "ArasbucvewgwivjuktkrjphognflfjFhfgggghfh",
    /*  64 '@' */  "Knjlhkflenepfrhsjrkqkkljmjnkomosnumvkwewcvbuasakbichegkg",
    /*  65 'A' */  "AgcmhwmmogCmmm",
    /*  66 'B' */  "AgigkhlimklmkniokplqmslukviwawBwbgboio",
    /*  67 'C' */  "Milhjgegchbiakasbucvewjwlvmu",
    /*  68 'D' */  "AgigkhlimkmslukviwawBwbg",
    /*  69 'E' */  "KwawagkgAoho",
    /*  70 'F' */  "KwawagAoho",
    /*  71 'G' */  "MilhjgegchbiakasbucvewjwlvmuMgmmHmom",
    /*  72 'H' */  "AgawKgkwAoko",
    /*  73 'I' */  "AgggAwgwDgdw",
    /*  74 'J' */  "AmakbichegghhiikiwDwnw",
    /*  75 'K' */  "AgawAmkwCokg",
    /*  76 'L' */  "Kgagaw",
    /*  77 'M' */  "Agawhmowog",
    /*  78 'N' */  "Agawmgmw",
    /*  79 'O' */  "Akbichegigkhlimkmslukviwewcvbuasak",
    /*  80 'P' */  "Agawiwkvlumsmrlpkoinan",
    /*  81 'Q' */  "AkbichegigkhlimkmslukviwewcvbuasakHlmg",
    /*  82 'R' */  "AgawiwkvlumsmrlpkoinanInmg",
    /*  83 'S' */  "Ajchegigkhlimklmkncpbqasbucvewiwkvmt",
    /*  84 'T' */  "HghwAwow",
    /*  85 'U' */  "Awakbichegigkhlimkmw",
    /*  86 'V' */  "Awggmw",
    /*  87 'W' */  "Aweghokgow",
    /*  88 'X' */  "AwkgAgkw",
    /*  89 'Y' */  "AwfokwFgfo",
    /*  90 'Z' */  "Awkwagkg",
    /*  91 '[' */  "Igdgdwiw",
    /*  92 '\' */  "Awmg",
    /*  93 ']' */  "Dgigiwdw",
    /*  94 '^' */  "Cqgwkq",
    /*  95 '_' */  "Bglg",
    /*  96 '`' */  "Cvdvdwcwcues",
    /*  97 'a' */  "BqdrirjqkokhlgKhigdgbhajbldmimkk",
    /*  98 'b' */  "BichegggihjikkkmjoipgqeqcpboAwbwbgag",
    /*  99 'c' */  "Jqhrercqbpanakbicheghgjh",
    /* 100 'd' */  "KiggegchbiakanbpcqerhrjqkoLgkgkwlw",
    /* 101 'e' */  "Kiihggegchbiakanbpcqerhrjqkokmam",
    /* 102 'f' */  "BgdgCgcrdtfuhtirAnfn",
    /* 103 'g' */  "KiihggegchbiakambodphpjokmLpkpkdjbhaeacbbc",
    /* 104 'h' */  "AwbwbgagBncpeqhqjpknkg",
    /* 105 'i' */  "AgegCgcpbpCrcsbsbrcr",
    /* 106 'j' */  "GshshrgrgsGphphegcebdbbcae",
    /* 107 'k' */  "AgbgbwawBjhpIgem",
    /* 108 'l' */  "AgegCgcubu",
    /* 109 'm' */  "AqbqbgBncpeqfqhpiniginjplqmqoppnpg",
    /* 110 'n' */  "AqbqbgagBncpeqhqjpknkg",
    /* 111 'o' */  "Akbichegggihjikkkmjoipgqeqcpboamak",
    /* 112 'p' */  "BichegggihjikkkmjoipgqeqcpboAqbqbaaa",
    /* 113 'q' */  "KiggegchbiakanbpcqerhrjqkoLpkpkdlbna",
    /* 114 'r' */  "AqbqbgBncpeqhqjo",
    /* 115 's' */  "Ahcgggihjjikbmanbpdqhqjp",
    /* 116 't' */  "AqiqEuehfghgih",
    /* 117 'u' */  "AqakbichegggihjikkKqkglg",
    /* 118 'v' */  "Arfgkr",
    /* 119 'w' */  "Areghmkgor",
    /* 120 'x' */  "ArkgAgkr",
    /* 121 'y' */  "AqakbichegggihjikkKqkdjbhaeacbbc",
    /* 122 'z' */  "Arkragkg",
    /* 123 '{' */  "Egchbjblcnaocpbrbtcvew",
    /* 124 '|' */  "Gggw",
    /* 125 '}' */  "Agchdjdlcneocpdrdtcvaw",
    /* 126 '~' */  "Arbtcueuftfsgrirjsku",
    /* 127 'Angstroms' */  "AgdkhqlkogDklkHwjuhrfuhw", 
    /* 128 'Degrees' */  "Aqbocnemgminjokqksjuivgwewcvbuasaq" 

        };

