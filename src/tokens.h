/***************************************************************************
 *                            RasMol 2.7.1.1                               *
 *                                                                         *
 *                                RasMol                                   *
 *                 Molecular Graphics Visualisation Tool                   *
 *                            17 January 2001                              *
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


/* tokens.h
 */

/* Lexeme Tokens */
#define IdentTok       256
#define NumberTok      257
#define FloatTok       258
#define StringTok      259

/* Command Tokens */
#define AdviseTok      260
#define BackboneTok    261
#define CartoonTok     262
#define CentreTok      263
#define ClipboardTok   264
#define ColourTok      265
#define ConnectTok     266
#define DashTok        267
#define DefineTok      268
#define DelayTok       269
#define DisplayTok     270
#define EchoTok        271
#define ExitTok        272
#define HelpTok        273
#define LabelTok       274
#define LoadTok        275
#define LoopTok        276
#define MonitorTok     277
#define MoveTok        278
#define PrintTok       279
#define QuitTok        280
#define RefreshTok     281
#define RenumTok       282
#define ResetTok       283
#define ResizeTok      284
#define RestoreTok     285
#define RestrictTok    286
#define RotateTok      287
#define SaveTok        288
#define ScriptTok      289
#define SelectTok      290
#define SetTok         291
#define ShowTok        292
#define SlabTok        293
#define SourceTok      294
#define SpacefillTok   295
#define StarTok        296
#define StructureTok   297
#define SymmetryTok    298
#define TitleTok       299
#define TraceTok       300
#define TranslateTok   301
#define ViewTok        302
#define WaitTok        303
#define WireframeTok   304
#define WriteTok       305
#define ZapTok         306
#define ZoomTok        307

/* Predicate Tokens */
#define IsPredTok(x)   (((x)>=310) && ((x)<=349))
#define PredTokOrd(x)  ((x)-310)
#define PredTokChr(x)  ((x)+310)

#define AlphaTok       310
#define AminoTok       311
#define ATTok          312
#define BondedTok      313
#define CGTok          314
#define CystineTok     315
#define DNATok         316
#define HelixTok       317
#define HeteroTok      318
#define HydrogenTok    319
#define IonTok         320
#define LigandTok      321
#define MainChainTok   322
#define NucleicTok     323
#define ProteinTok     324
#define PurineTok      325
#define PyrimidineTok  326
#define RNATok         327
#define SelectedTok    328
#define SheetTok       329
#define SidechainTok   330
#define SolventTok     331
#define TurnTok        332
#define WaterTok       333

#define AcidicTok      334
#define AcyclicTok     335
#define AliphaticTok   336
#define AromaticTok    337
#define BasicTok       338
#define BuriedTok      339
#define ChargedTok     340
#define CisBondedTok   341
#define CyclicTok      342
#define HydrophobicTok 343
#define LargeTok       344
#define MediumTok      345
#define NeutralTok     346
#define PolarTok       347
#define SmallTok       348
#define SurfaceTok     349


/* Property Tokens */
#define IsPropTok(x)   (((x)>=350) && ((x)<=356))
#define TemperatureTok 350
#define RadiusTok      351
#define AtomNoTok      352
#define ElemNoTok      353
#define ModelTok       354
#define ResNoTok       355
#define AltlTok        356

/* File Format Tokens */
/* Warning! Tokens are related to Format values */
#define IsMoleculeToken(x)  (((x)>=360) && ((x)<=375))

#define PDBTok         360
#define MacroModelTok  361
#define GaussianTok    362
#define AlchemyTok     363
#define NMRPDBTok      364
#define CharmmTok      365
#define BiosymTok      366
#define MOPACTok       367
#define SHELXTok       368
#define Mol2Tok        369
#define FDATTok        370
#define MMDBTok        371
#define MDLTok         372
#define XYZTok         373
#define CIFTok         374
#define CEXTok         375

/* Raster Tokens */
#define IsImageToken(x) (((x)>=380) && ((x)<=398) || ((x) == PhiPsiTok) )
#define GIFTok         380
#define PPMTok         381
#define SUNTok         382
#define SUNRLETok      383
#define EPSFTok        384
#define PICTTok        385
#define IRISTok        386
#define BMPTok         387
#define MonoPSTok      388
#define JPEGTok        389
#define PNGTok         390
#define VectPSTok      391
#define KinemageTok    392
#define MolScriptTok   393
#define POVRayTok      394
#define POVRay2Tok     394
#define POVRay3Tok     395
#define VRMLTok        396
#define RamachanTok    397  /* ok, this isn't a real image format ... */
#define RamPrintTok    398

/* Feature Tokens */
#define AtomTok        400
#define BondTok        401
#define DotsTok        402
#define HBondTok       403
#define RibbonTok      404
#define SSBondTok      405
#define Ribbon1Tok     406
#define Ribbon2Tok     407

/* Expression Tokens */
#define TrueTok        410
#define FalseTok       411
#define AllTok         412
#define NoneTok        413
#define AndTok         414
#define OrTok          415
#define NotTok         416
#define WithinTok      417
#define XorTok         418

/* Colour Tokens */
/* Warning! Tokens are related to colour values */
#define IsColourToken(x) (((x)>=420) && ((x)<=443))
#define Token2Colour(x)  ((x)-420)

#define BlackTok       420
#define BlueTok        421
#define BlueTintTok    422
#define BrownTok       423
#define CyanTok        424
#define GoldTok        425
#define GrayTok        426
#define GreenTok       427
#define GreenBlueTok   428
#define GreenTintTok   429
#define HotPinkTok     430
#define MagentaTok     431
#define OrangeTok      432
#define PinkTok        433
#define PinkTintTok    434
#define PurpleTok      435
#define RedTok         436
#define RedOrangeTok   437
#define SeaGreenTok    438
#define SkyBlueTok     439
#define VioletTok      440
#define WhiteTok       441
#define YellowTok      442
#define YellowTintTok  443

#define CPKTok         444
#define ShapelyTok     445
#define ResidueTok     446
#define UserTok        447
#define GroupTok       448
#define ChainTok       449
#define TypeTok        450
#define PotentialTok   451
#define ChargeTok      452

/* Variable Tokens */
#define AmbientTok     460
#define AxesTok        461
#define BackFadeTok    462
#define BackgroundTok  463
#define BondModeTok    464
#define BoundBoxTok    465
#define CisAngleTok    466
#define DepthCueTok    467
#define FontSizeTok    468
#define FontStrokeTok  469
#define HourGlassTok   470
#define MenusTok       471
#define MouseTok       472
#define PickingTok     473
#define ShadowTok      474
#define SlabModeTok    475
#define SpecularTok    476
#define SpecPowerTok   477
#define StrandsTok     478
#define TransparentTok 479
#define UnitCellTok    480

/* SlabMode Tokens */
#define RejectTok      481
#define HalfTok        482
#define HollowTok      483
#define SolidTok       484
#define SectionTok     485

/* MouseMode Tokens */
#define RasMolTok      486
#define InsightTok     487
#define QuantaTok      488
#define SybylTok       489

/* Information Tokens */
#define InfoTok        490
#define SequenceTok    491
#define VersionTok     492
#define PhiPsiTok      493

/* Display Mode Tokens */
#define NormalTok      495
#define StereoTok      496
#define MonoTok        497
#define HardwareTok    498

/* Axis Tokens */
#define XTok           500
#define YTok           501
#define ZTok           502

/* Picking Tokens */
#define IdentifyTok    505
#define CoordTok       506
#define DistanceTok    507
#define AngleTok       508
#define TorsionTok     509
#define OriginTok      510

/* Misc Tokens */
#define InLineTok      511
#define VDWTok         512
#define HeaderTok      513
#define CIFDataTok     514
#define FSTok          515
#define PSTok          EPSFTok

/* Language Tokens */
#define EnglishTok     600
#define SpanishTok     601

int LookUpKeyword( char *ptr );

