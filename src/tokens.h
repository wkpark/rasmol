/* tokens.h
 * RasMol2 Molecular Graphics
 * Roger Sayle, October 1994
 * Version 2.5
 */

/* Lexeme Tokens */
#define IdentTok       256
#define NumberTok      257
#define FloatTok       258
#define StringTok      259

/* Command Tokens */
#define BackboneTok    260
#define CentreTok      261
#define ClipboardTok   262
#define ColourTok      263
#define ConnectTok     264
#define DefineTok      265
#define DisplayTok     266
#define EchoTok        267
#define HelpTok        268
#define LabelTok       269
#define LoadTok        270
#define PrintTok       271
#define QuitTok        272
#define RenumTok       273
#define ResetTok       274
#define ResizeTok      275
#define RestrictTok    276
#define RotateTok      277
#define SaveTok        278
#define ScriptTok      279
#define SelectTok      280
#define SetTok         281
#define ShowTok        282
#define SlabTok        283
#define SourceTok      284
#define SpacefillTok   285
#define StructureTok   286
#define SymmetryTok    287
#define TraceTok       288
#define TranslateTok   289
#define WaitTok        290
#define WireframeTok   291
#define WriteTok       292
#define ZapTok         293
#define ZoomTok        294

/* Predicate Tokens */
#define IsPredTok(x)   (((x)>=300) && ((x)<=338))
#define PredTokOrd(x)  ((x)-300)
#define PredTokChr(x)  ((x)+300)

#define AlphaTok       300
#define AminoTok       301
#define ATTok          302
#define BondedTok      303
#define CGTok          304
#define CystineTok     305
#define DNATok         306
#define HelixTok       307
#define HeteroTok      308
#define HydrogenTok    309
#define IonTok         310
#define LigandTok      311
#define MainChainTok   312
#define NucleicTok     313
#define ProteinTok     314
#define PurineTok      315
#define PyrimidineTok  316
#define RNATok         317
#define SelectedTok    318
#define SheetTok       319
#define SidechainTok   320
#define SolventTok     321
#define TurnTok        322
#define WaterTok       323

#define AcidicTok      324
#define AcyclicTok     325
#define AliphaticTok   326
#define AromaticTok    327
#define BasicTok       328
#define BuriedTok      329
#define ChargedTok     330
#define CyclicTok      331
#define HydrophobicTok 332
#define LargeTok       333
#define MediumTok      334
#define NeutralTok     335
#define PolarTok       336
#define SmallTok       337
#define SurfaceTok     338


/* Property Tokens */
#define IsPropTok(x)   (((x)>=340) && ((x)<=344))
#define TemperatureTok 340
#define RadiusTok      341
#define AtomNoTok      342
#define ElemNoTok      343
#define ResNoTok       344

/* File Format Tokens */
#define IsMoleculeFormat(x)  (((x)>=350) && ((x)<=356))
#define PDBTok         350
#define AlchemyTok     351
#define CharmmTok      352
#define Mol2Tok        353
#define XYZTok         354
#define CIFTok         355
#define MDLTok         356

/* Raster Tokens */
#define IsImageFormat(x) (((x)>=360) && ((x)<=372))
#define GIFTok         360
#define PPMTok         361
#define SUNTok         362
#define SUNRLETok      363
#define EPSFTok        364
#define PICTTok        365
#define IRISTok        366
#define BMPTok         367
#define MonoPSTok      368
#define VectPSTok      369
#define KinemageTok    370
#define MolScriptTok   371
#define POVRayTok      372

/* Feature Tokens */
#define AtomTok        380
#define BondTok        381
#define DotsTok        382
#define HBondTok       383
#define RibbonTok      384
#define SSBondTok      385
#define Ribbon1Tok     386
#define Ribbon2Tok     387

/* Expression Tokens */
#define TrueTok        390
#define FalseTok       391
#define AllTok         392
#define NoneTok        393
#define AndTok         394
#define OrTok          395
#define NotTok         396
#define WithinTok      397
#define XorTok         398

/* Colour Tokens */
#define BlueTok        400
#define BlueTintTok    401
#define BlackTok       402
#define BrownTok       403
#define CyanTok        404
#define GoldTok        405
#define GrayTok        406
#define GreenTok       407
#define GreenblueTok   408
#define GreenTintTok   409
#define HotPinkTok     410
#define MagentaTok     411
#define OrangeTok      412
#define PinkTok        413
#define PinkTintTok    414
#define PurpleTok      415
#define RedTok         416
#define RedorangeTok   417
#define SeaTok         418
#define SkyTok         419
#define VioletTok      420
#define WhiteTok       421
#define YellowTok      422
#define YellowTintTok  423

#define CPKTok         424
#define ShapelyTok     425
#define UserTok        426
#define GroupTok       427
#define ChainTok       428
#define TypeTok        429
#define PotentialTok   430
#define ChargeTok      431

/* Variable Tokens */
#define AmbientTok     440
#define AxesTok        441
#define BackgroundTok  442
#define BondModeTok    443
#define BoundBoxTok    444
#define FontSizeTok    445
#define HourGlassTok   446
#define MenusTok       447
#define MouseTok       448
#define ShadowTok      449
#define SlabModeTok    450
#define SpecularTok    451
#define SpecPowerTok   452
#define StrandsTok     453
#define UnitCellTok    454

/* SlabMode Tokens */
#define RejectTok      460
#define HalfTok        461
#define HollowTok      462
#define SolidTok       463
#define SectionTok     464

/* MouseMode Tokens */
#define RasMolTok      465
#define InsightTok     466
#define QuantaTok      467

/* Information Tokens */
#define InfoTok        470
#define SequenceTok    471
#define VersionTok     472

/* Display Mode Tokens */
#define NormalTok      475
#define StereoTok      476
#define MonoTok        477

/* Axis Tokens */
#define XTok           480
#define YTok           481
#define ZTok           482


typedef struct {
                char *ident;
                int token;
               } KeywordEntry;


#define MAXKEYLEN 11
static int KeyLen[MAXKEYLEN+1] = {
        0, 3, 8, 27, 57, 88, 134, 162, 182, 197, 200, 205 };

static KeywordEntry Keyword[] = {
            { "X",  XTok },
            { "Y",  YTok },
            { "Z",  ZTok },  /* 3 */

            { "AT", ATTok   },
            { "CG", CGTok   },
            { "ON", TrueTok },
            { "OR", OrTok   },
            { "PS", EPSFTok },  /* 8 */

            { "ALL", AllTok   },
            { "AND", AndTok   },
            { "BMP", BMPTok   },
            { "CIF", CIFTok   },
            { "CPK", CPKTok   },
            { "DNA", DNATok   },
            { "GIF", GIFTok   },
            { "ION", IonTok   },
            { "MDL", MDLTok   },
            { "NOT", NotTok   },
            { "OFF", FalseTok },
            { "PDB", PDBTok   },
            { "PPM", PPMTok   },
            { "RED", RedTok   },
            { "RNA", RNATok   },
            { "SET", SetTok   },
            { "SUN", SUNTok   },
            { "XYZ", XYZTok   },
            { "ZAP", ZapTok   }, /* 27 */

            { "ATOM", AtomTok },
            { "AXES", AxesTok },
            { "AXIS", AxesTok },
            { "BLUE", BlueTok },
            { "BOND", BondTok },
            { "CYAN", CyanTok },
            { "DOTS", DotsTok },
            { "ECHO", EchoTok },
            { "EPSF", EPSFTok },
            { "EXIT", QuitTok },
            { "HALF", HalfTok },
            { "HELP", HelpTok },
            { "INFO", InfoTok },
            { "IONS", IonTok  },
            { "IRIS", IRISTok },
            { "LOAD", LoadTok },
            { "MOL2", Mol2Tok },
            { "MONO", MonoTok },
            { "NONE", NoneTok },
            { "PICT", PICTTok },
            { "QUIT", QuitTok },
            { "SAVE", SaveTok },
            { "SHOW", ShowTok },
            { "SLAB", SlabTok },
            { "TRUE", TrueTok },
            { "TURN", TurnTok },
            { "TYPE", TypeTok },
            { "USER", UserTok },
            { "WAIT", WaitTok },
            { "ZOOM", ZoomTok }, /* 57 */

            { "ALPHA", AlphaTok  },
            { "AMINO", AminoTok  },
            { "ATOMS", AtomTok   },
            { "BASIC", BasicTok  },
            { "BLACK", BlackTok  },
            { "BONDS", BondTok   },
            { "CHAIN", ChainTok  },
            { "COLOR", ColourTok },
            { "FALSE", FalseTok  },
            { "GREEN", GreenTok  },
            { "GROUP", GroupTok  },
            { "HBOND", HBondTok  },
            { "HELIX", HelixTok  },
            { "LABEL", LabelTok  },
            { "LARGE", LargeTok  },
            { "MENUS", MenusTok  },
            { "MOUSE", MouseTok  },
            { "PAUSE", WaitTok   },
            { "POLAR", PolarTok  },
            { "PRINT", PrintTok  },
            { "RENUM", RenumTok  },
            { "RESET", ResetTok  },
            { "RESNO", ResNoTok  },
            { "SHEET", SheetTok  },
            { "SMALL", SmallTok  },
            { "SOLID", SolidTok  },
            { "TRACE", TraceTok  },
            { "TURNS", TurnTok   },
            { "WATER", WaterTok  },
            { "WHITE", WhiteTok  },
            { "WRITE", WriteTok  },  /* 88 */

            { "ACIDIC", AcidicTok },
            { "ATOMNO", AtomNoTok },
            { "BONDED", BondedTok },
            { "BURIED", BuriedTok },
            { "CENTER", CentreTok },
            { "CENTRE", CentreTok },
            { "CHARGE", ChargeTok },
            { "CHARMM", CharmmTok },
            { "COLORS", ColourTok },
            { "COLOUR", ColourTok },
            { "CYCLIC", CyclicTok },
            { "DEFINE", DefineTok },
            { "ELEMNO", ElemNoTok },
            { "HBONDS", HBondTok  },
            { "HETERO", HeteroTok },
            { "HOLLOW", HollowTok },
            { "LABELS", LabelTok  },
            { "LIGAND", LigandTok },
            { "MEDIUM", MediumTok },
            { "MONOPS", MonoPSTok },
            { "NORMAL", NormalTok },
            { "ORANGE", OrangeTok },
            { "POVRAY", POVRayTok },
            { "PURINE", PurineTok },
            { "PURPLE", PurpleTok },
            { "QUANTA", QuantaTok },
            { "RADIUS", RadiusTok },
            { "RASMOL", RasMolTok },
            { "RASWIN", RasMolTok },
            { "REJECT", RejectTok },
            { "RESIZE", ResizeTok },
            { "RIBBON", RibbonTok },
            { "ROTATE", RotateTok },
            { "SCRIPT", ScriptTok },
            { "SELECT", SelectTok },
            { "SHADOW", ShadowTok },
            { "SHEETS", SheetTok  },
            { "SOURCE", SourceTok },
            { "SSBOND", SSBondTok },
            { "STEREO", StereoTok },
            { "SUNRLE", SUNRLETok },
            { "VECTPS", VectPSTok },
            { "VIOLET", VioletTok },
            { "WATERS", WaterTok  },
            { "WITHIN", WithinTok },
            { "YELLOW", YellowTok },  /* 134 */

            { "ACYCLIC", AcyclicTok },
            { "ALCHEMY", AlchemyTok },
            { "AMBIENT", AmbientTok },
            { "CHARGED", ChargedTok },
            { "CHARGES", ChargeTok  },
            { "COLOURS", ColourTok  },
            { "CONNECT", ConnectTok },
            { "CYSTINE", CystineTok },
            { "DISPLAY", DisplayTok },
            { "HELICES", HelixTok   },
            { "INSIGHT", InsightTok },
            { "LIGANDS", LigandTok  },
            { "MAGENTA", MagentaTok },
            { "NEUTRAL", NeutralTok },
            { "NUCLEIC", NucleicTok },
            { "PROTEIN", ProteinTok },
            { "PURINES", PurineTok  },
            { "RESIDUE", GroupTok   },
            { "RIBBON1", Ribbon1Tok },
            { "RIBBON2", Ribbon2Tok },
            { "RIBBONS", RibbonTok  },
            { "SECTION", SectionTok },
            { "SHADOWS", ShadowTok  },
            { "SHAPELY", ShapelyTok },
            { "SOLVENT", SolventTok },
            { "SSBONDS", SSBondTok  },
            { "STRANDS", StrandsTok },
            { "SURFACE", SurfaceTok },  /* 162 */

            { "AROMATIC", AromaticTok },
            { "BACKBONE", BackboneTok },
            { "BONDMODE", BondModeTok },
            { "BOUNDBOX", BoundBoxTok },
            { "FONTSIZE", FontSizeTok },
            { "HYDROGEN", HydrogenTok },
            { "KINEMAGE", KinemageTok },
            { "NEGATIVE", AcidicTok   },
            { "POSITIVE", BasicTok    },
            { "RENUMBER", RenumTok    },
            { "RESTRICT", RestrictTok },
            { "RIBBONS1", Ribbon1Tok },
            { "RIBBONS2", Ribbon2Tok },
            { "SELECTED", SelectedTok },
            { "SEQUENCE", SequenceTok },
            { "SLABMODE", SlabModeTok },
            { "SOLVENTS", SolventTok  },
            { "SPECULAR", SpecularTok }, 
            { "SYMMETRY", SymmetryTok },
            { "UNITCELL", UnitCellTok },  /* 182 */

            { "ALIPHATIC", AliphaticTok },
            { "CLIPBOARD", ClipboardTok },
            { "GREENBLUE", GreenblueTok },
            { "HOURGLASS", HourGlassTok },
            { "MAINCHAIN", MainChainTok },
            { "MOLSCRIPT", MolScriptTok },
            { "MOUSEMODE", MouseTok     },
            { "POTENTIAL", PotentialTok },
            { "REDORANGE", RedorangeTok },
            { "SIDECHAIN", SidechainTok },
            { "SPACEFILL", SpacefillTok },
            { "SPECPOWER", SpecPowerTok },
            { "STRUCTURE", StructureTok },
            { "TRANSLATE", TranslateTok },
            { "WIREFRAME", WireframeTok },  /* 197 */

            { "BACKGROUND", BackgroundTok },
            { "MONOCHROME", MonoTok       },
            { "PYRIMIDINE", PyrimidineTok },  /* 200 */

            { "BOUNDINGBOX", BoundBoxTok    },
            { "HYDROPHOBIC", HydrophobicTok },
            { "INFORMATION", InfoTok        },
            { "PYRIMIDINES", PyrimidineTok, },
            { "TEMPERATURE", TemperatureTok }  /* 205 */
                };

