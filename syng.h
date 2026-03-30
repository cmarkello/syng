/*  File: syng.h
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2023
 *-------------------------------------------------------------------
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Mar 12 22:48 2026 (rd109)
 * Created: Mon May 29 08:19:18 2023 (rd109)
 *-------------------------------------------------------------------
 */

#include "utils.h"
#include "array.h"
#include "hash.h"
#include "dict.h"
#include "ONElib.h"
#include "syncmerset.h"

#define SYNG_VERSION  "2.1"

typedef struct {
  int   fixedLen ;
  Array node ;   // of Node - only defined in syngbwt.c
  Array status ; // of U8 bitFlags
  Array length ; // of I32, if fixedLen == 0
  // the next annoying set of properties are to manage start counts
  Hash  startHash ;      // for nodes with starts - good if relatively few starts compared to nodes
  Array startHashCount ; // of I32, indexed by hash value
} SyngBWT ;

typedef struct {
  SyngBWT *sb ;
  I32      lastNode, thisNode ;
  U32      lastOff, jLast, jMax ;
} SyngBWTpath ;

// in syngbwt.c

SyngBWT       *syngBWTcreate (int fixedLen, I64 max) ;
void           syngBWTdestroy (SyngBWT *sb) ;
void           syngBWTwrite (OneFile *of, SyngBWT *sb) ;
SyngBWT       *syngBWTread  (OneFile *of) ;
SyngBWTpath   *syngBWTpathStartNew (SyngBWT *sb,   I32 startNode) ;
void           syngBWTpathAdd (SyngBWTpath *sbp,   I32 nextNode, U32 offset) ; // add to a new path
void           syngBWTpathFinish (SyngBWTpath *sbp) ; // use when creating a new path
SyngBWTpath   *syngBWTpathStartOld (SyngBWT *sb,   I32 startNode, U32 count) ; // follow existing path
bool           syngBWTpathNext (SyngBWTpath *sbp,  I32 *nextNode, U32 *nextPos) ;
SyngBWTpath   *syngBWTmatchStart (SyngBWT *sb,     I32 startNode, U32 *high) ;
bool           syngBWTmatchNext (SyngBWTpath *sbp, I32 nextNode, U32 nextOff, U32 *low, U32 *high) ;
void           syngBWTpathDestroy (SyngBWTpath *sbp) ;
void           syngBWTstat (SyngBWT *sb) ;

// Annotation-aware spliced pangenome graph structures
// Inspired by vg rna (https://github.com/vgteam/vg/wiki/Transcriptomic-analyses)

// A splice edge connects an exon donor site to an exon acceptor site in the syncmer graph.
// These are additional edges layered on top of the genomic pangenome graph.
typedef struct {
  I32 donorNode ;      // syncmer node at the donor (end of upstream exon)
  I32 acceptorNode ;   // syncmer node at the acceptor (start of downstream exon)
  U32 donorOffset ;    // offset within the donor syncmer node
  U32 acceptorOffset ; // offset within the acceptor syncmer node
  I32 transcriptIdx ;  // index into transcript list
  I32 iSample ;        // which sample/haplotype this junction belongs to
} SyncSpliceEdge ;

// A transcript path through the syncmer graph, representing a single transcript
// across exons connected by splice junctions. Analogous to GBWT paths in vg.
typedef struct {
  char *transcriptId ;
  char *geneId ;
  I32   iSample ;       // sample/haplotype index
  I32   nNodes ;         // number of syncmer nodes in this transcript path
  I32  *nodes ;          // array of syncmer node IDs (negative if reverse strand)
  U32  *offsets ;        // array of offsets between consecutive nodes
  I32   nExons ;         // number of exons
  I32  *exonBoundaries ; // indices into nodes[] where each exon starts
} SyncTranscriptPath ;

// The spliced pangenome graph extends SyngBWT with transcript annotation.
// It layers splice junction edges and transcript paths onto the base syncmer graph.
typedef struct {
  SyngBWT    *baseBWT ;       // the underlying genomic syncmer graph
  SyncmerSet *sms ;           // the syncmer set used for node lookups
  Array       spliceEdges ;   // of SyncSpliceEdge
  Array       txPaths ;       // of SyncTranscriptPath
  DICT       *txDict ;        // transcript_id -> index into txPaths
  DICT       *geneDict ;      // gene_id -> index (for lookup)
  int         nSamples ;      // number of samples with annotations
  char      **sampleNames ;   // sample names parallel to annotation files
} SplicedSyngBWT ;

SplicedSyngBWT *splicedSyngBWTcreate (SyngBWT *baseBWT, SyncmerSet *sms) ;
void            splicedSyngBWTdestroy (SplicedSyngBWT *ssbwt) ;
void            splicedSyngBWTaddTranscriptPath (SplicedSyngBWT *ssbwt,
                    SyncTranscriptPath *txPath) ;
void            splicedSyngBWTwrite (OneFile *of, SplicedSyngBWT *ssbwt) ;

static char *syngSchemaText =
  "1 3 def 1 0               schema for syng\n"
  ".\n"
  "P 3 seq                   SEQUENCE\n"
  "S 4 path                  contains P (path) objects = syncmer sequences\n"
  "S 3 gfa                   sequence graph - contains V (vertex) objects, probably with E lines\n"
  "S 4 gbwt                  gbwt: a gfa with B, C, Z lines\n"
  ".\n"
  "D h 3 3 INT 3 INT 3 INT   k, w, seed for the seqhash: for syncs k = |smer|, w+k = |syncmer|\n"
  ".\n"
  "O S 1 3 DNA               sequence of the node\n" // for the future, for general GFA
  ".\n"
  "O V 1 3 INT               graph node (vertex): length\n"
  "D K 1 3 INT               coverage count of node\n"
  "D E 3 3 INT 3 INT 3 INT   edge +: adjacent node (- if reversed), offset, count\n"
  "D e 3 3 INT 3 INT 3 INT   edge -: adjacent node (- if reversed), offset, count\n"
  "D B 1 8 INT_LIST          GBWT +: list of node indices from opposite-signed E lines\n"
  "D b 1 8 INT_LIST          GBWT -: list of node indices from opposite-signed E lines\n"
  "D C 1 8 INT_LIST          GBWT +: list of run-length counts\n"
  "D c 1 8 INT_LIST          GBWT -: list of run-length counts\n"
  ".\n"
  "O P 3 3 INT 3 INT 3 INT   path: length in bp, source file number, sequence number in file\n"
  "D Z 4 3 INT 3 INT 3 INT 3 INT   GBWT path: starting node, pos, count, then length in nodes\n"
  "D z 1 8 INT_LIST          alternative explicit list of node ids (-ve if reversed)\n"
  "D o 1 8 INT_LIST          if z, then offsets of the nodes from start of sequence, 1:1 with z\n"
  "D X 1 3 DNA               prefix before first node - required to fully reconstruct\n"
  "D Y 1 3 DNA               suffix after last node - required to fully reconstruct\n"
  ".\n"
  "P 5 khash                 KMER HASH\n"
  "S 7 syncset               SYNCMER SET\n"
  "D h 3 3 INT 3 INT 3 INT   k, w, seed for the seqhash: for syncs k = |smer|, w+k = |syncmer|\n"
  "O t 3 3 INT 3 INT 3 INT   max, len, dim for KmerHash table\n"
  "D S 1 3 DNA               packed sequences aligned to 64-bit boundaries\n" 
  "D L 1 8 INT_LIST          locations in the table\n"
  "D C 1 8 INT_LIST          kmer counts\n"
  "D M 1 6 STRING            maximum count in any input - (1..127)\n"
  ".\n"
  "P 3 map                   MAP\n"
  "P 3 ref                   REFERENCE INFORMATION\n"
  "O N 1 6 STRING            name of sequence, e.g. chr1\n"
  "O I 1 8 INT_LIST          list of indexes into name table for each sync: 0 if missing, -ve if RC\n"
  "D P 1 8 INT_LIST          positions in sequence\n"
  ".\n"
  "P 6 sgtxom                SPLICED PANGENOME TRANSCRIPTOME\n"
  "S 4 stxn                  spliced transcriptome: annotations over a gbwt\n"
  "D h 3 3 INT 3 INT 3 INT   k, w, seed for the seqhash\n"
  ".\n"
  "O T 1 3 INT               transcript: sample index\n"
  "D i 1 6 STRING            transcript_id\n"
  "D g 1 6 STRING            gene_id for this transcript\n"
  "D n 1 8 INT_LIST          node list for transcript path through syncmer graph\n"
  "D f 1 8 INT_LIST          offset list for transcript path (1:1 with n)\n"
  "D x 1 8 INT_LIST          exon boundary indices into node list\n"
  ".\n"
  "O J 4 3 INT 3 INT 3 INT 3 INT   splice junction: donor_node, acceptor_node, donor_off, acceptor_off\n"
  "D a 1 3 INT               sample index for this splice junction\n"
  ;


/****************** end of file ********************/
