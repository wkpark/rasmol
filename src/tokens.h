/* tokens.h
 * RasMol2 Molecular Graphics
 * Roger Sayle, February 1994
 * Version 2.3
 */

/* Lexeme Tokens */
#define IdentTok   256
#define NumberTok  257
#define StringTok  258

/* Command Tokens */
#define CentreTok     260
#define ColourTok     261
#define DefineTok     262
#define DisplayTok    263
#define EchoTok       264
#define HelpTok       265
#define LoadTok       266
#define QuitTok       267
#define RenumTok      268
#define ResetTok      269
#define ResizeTok     270
#define RestrictTok   271
#define RotateTok     272
#define SaveTok       273
#define ScriptTok     274
#define SelectTok     275
#define SetTok        276
#define ShowTok       277
#define SlabTok       278
#define SpacefillTok  279
#define StructureTok  280
#define SymmetryTok   281
#define TraceTok      282
#define TranslateTok  283
#define WaitTok       284
#define WireframeTok  285
#define WriteTok      286
#define ZapTok        287
#define ZoomTok       288

/* Feature Tokens */
#define AtomTok      290
#define BondTok      291
#define SSBondTok    292
#define HBondTok     293
#define RibbonTok    294

/* Expression Tokens */
#define TrueTok      300
#define FalseTok     301
#define AllTok       302
#define NoneTok      303
#define AndTok       304
#define OrTok        305
#define NotTok       306
#define WithinTok    307

/* Axis Tokens */
#define XTok  310
#define YTok  311
#define ZTok  312

/* Raster Tokens */
#define IsImageFormat(x) (((x)>=320) && ((x)<=328))

#define GIFTok        320
#define PPMTok        321
#define SUNTok        322
#define SUNRLETok     323
#define EPSFTok       324
#define BMPTok        325
#define MonoPSTok     326
#define VectPSTok     327
#define MolScriptTok  328

/* SlabMode Tokens */
#define RejectTok     330
#define HalfTok       331
#define HollowTok     332
#define SolidTok      333
#define SectionTok    334

/* MouseMode Tokens */
#define RasMolTok     335
#define InsightTok    336
#define QuantaTok     337

/* Colour Tokens */
#define BlueTok       340
#define BlackTok      341
#define CyanTok       342
#define GreenTok      343
#define GreenblueTok  344
#define MagentaTok    345
#define OrangeTok     346
#define PurpleTok     347
#define RedTok        348
#define RedorangeTok  349
#define VioletTok     350
#define WhiteTok      351
#define YellowTok     352

#define CPKTok        353
#define ShapelyTok    354
#define UserTok       355
#define GroupTok      356
#define ChainTok      357
#define TypeTok       358

/* Variable Tokens */
#define AmbientTok     360
#define AxesTok        361
#define BackgroundTok  362
#define BondModeTok    363
#define BoundBoxTok    364
#define HourGlassTok   365
#define MouseTok       366
#define ShadowTok      367
#define SlabModeTok    368
#define SpecularTok    369
#define SpecPowerTok   370
#define StrandsTok     371
#define UnitCellTok    372

/* Predicate Tokens */
#define IsPredTok(x)   (((x)>=380) && ((x)<=418))
#define PredTokOrd(x)  ((x)-380)
#define PredTokChr(x)  ((x)+380)

#define AlphaTok       380
#define AminoTok       381
#define ATTok          382
#define BackboneTok    383
#define BondedTok      384
#define CGTok          385
#define CystineTok     386
#define DNATok         387
#define HelixTok       388
#define HeteroTok      389
#define HydrogenTok    390
#define IonTok         391
#define LigandTok      392
#define NucleicTok     393
#define ProteinTok     394
#define PurineTok      395
#define PyrimidineTok  396
#define RNATok         397
#define SelectedTok    398
#define SheetTok       399
#define SidechainTok   400
#define SolventTok     401
#define TurnTok        402
#define WaterTok       403

#define AcidicTok      404
#define AcyclicTok     405
#define AliphaticTok   406
#define AromaticTok    407
#define BasicTok       408
#define BuriedTok      409
#define ChargedTok     410
#define CyclicTok      411
#define HydrophobicTok 412
#define LargeTok       413
#define MediumTok      414
#define NeutralTok     415
#define PolarTok       416
#define SmallTok       417
#define SurfaceTok     418


/* Property Tokens */
#define IsPropTok(x)   (((x)>=420) && ((x)<=423))
#define TemperatureTok 420
#define RadiusTok      421
#define AtomNoTok      422
#define ResNoTok       423

/* Information Tokens */
#define InfoTok        430
#define SequenceTok    431
#define VersionTok     432

/* Display Mode Tokens */
#define NormalTok      435
#define StereoTok      436
#define MonoTok        437

/* File Format Tokens */
#define IsMoleculeFormat(x)  (((x)>=440) && ((x)<=442))

#define PDBTok         440
#define XYZTok         441
#define AlchemyTok     442
#define CIFTok         443
#define MDLTok         444

