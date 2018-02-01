#define ADD 257
#define AND 258
#define VAR 259
#define QUIT 260
#define DEL 261
#define LIST 262
#define VARS 263
#define DEPENDS 264
#define ON 265
#define IN 266
#define IS 267
#define STRUCTURE 268
#define OF 269
#define FAM 270
#define FOR 271
#define RULE 272
#define ALL 273
#define DESCRIPTIONS 274
#define SET 275
#define CURRENT 276
#define SHOW 277
#define OPT 278
#define NODE 279
#define DERIVE 280
#define SEL 281
#define LOAD 282
#define UNDEF 283
#define TEST 284
#define EVALUATE 285
#define EXPECT 286
#define COMMENT 287
#define LEARN 288
#define STAT 289
#define GA 290
#define POPULATION 291
#define POINTS 292
#define PLOT 293
#define MAN 294
#define WAIT 295
#define DEBUG 296
#define ITERATIONS 297
#define LEX 298
#define PRINT 299
#define FREQ 300
#define ERROR 301
#define MAX 302
#define YECHO 303
#define SAVE 304
#define IMP 305
#define WEIGHT 306
#define MIN 307
#define MULT 308
#define INFORMATIVITY 309
#define GAABS 310
#define YES 311
#define NO 312
#define PREFER 313
#define MUTATION 314
#define POLICY 315
#define CUSTOM 316
#define OR 317
#define LSQR 318
#define RESET 319
#define CROSSOVER 320
#define USAGE 321
#define NORM 322
#define WIDTH 323
#define TABLE 324
#define SPLIT 325
#define DECOMPOSE 326
#define REDO 327
#define GAMMA 328
#define TS 329
#define USING 330
#define CL 331
#define SIGM 332
#define KNN 333
#define COLOR 334
#define COMPUTE 335
#define ONE 336
#define COMPOSE 337
#define HEURISTIC 338
#define COPY 339
#define TO 340
#define QUALITATIVE 341
#define FROM 342
#define QUANTITATIVE 343
#define COMPARE 344
#define STEP 345
#define REPEAT 346
#define YLOG 347
#define GR 348
#define GINI 349
#define RELIEFF 350
#define FUZZY 351
#define CRISP 352
#define PREFIX 353
#define JOIN 354
#define INTERVAL 355
#define INSTANCE 356
#define TREE 357
#define DM 358
#define PJOIN 359
#define FJOIN 360
#define SURE 361
#define NEED 362
#define REDUNDANT 363
#define VECT 364
#define SAME 365
#define GROUP 366
#define ANALYZE 367
#define NAMES 368
#define WRITE 369
#define DIFFERENT 370
#define SOFT 371
#define MNUMBER 372
#define DOT 373
#define COUNT 374
#define LOCAL 375
#define GLOBAL 376
#define SDTIC 377
#define DFC 378
#define DTIC 379
#define CM 380
#define MERGE 381
#define PAIRS 382
#define LAPLACE 383
#define CLUSTER 384
#define NOISE 385
#define SEED 386
#define DONTCARE 387
#define DONTKNOW 388
#define APRIORY 389
#define DISTRIBUTION 390
#define CV 391
#define SAMPLE 392
#define SORT 393
#define MDL 394
#define LEAF 395
#define NOD 396
#define DUPLICATE 397
#define EXPAND 398
#define NUM 399
#define INUM 400
#define STRING 401
#define FNAME 402
#define ID 403
#ifdef YYSTYPE
#undef  YYSTYPE_IS_DECLARED
#define YYSTYPE_IS_DECLARED 1
#endif
#ifndef YYSTYPE_IS_DECLARED
#define YYSTYPE_IS_DECLARED 1
typedef union {
  int intV;
  double doubV;
  char chV;
  Str255 strV;
  char *chpV;
  list_of_vars *varL;
  list_of_str *strL;
  list_of_num *numL;
  list_of_q *strQ;
  fa_type faV;
  ga_error_type gaeV;
  ga_policy_type gapV;
  eval_method_type evalV;
  coloring_type colV;
  er_type erV;
} YYSTYPE;
#endif /* !YYSTYPE_IS_DECLARED */
extern YYSTYPE yylval;
