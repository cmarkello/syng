/*  File: gffio.h
 *  Author: cmarkello
 *  Copyright (C) 2026
 *-------------------------------------------------------------------
 * Description: GFF/GTF annotation file parser for gene annotation
 *   in spliced pangenome graph construction.
 *   Supports GFF3 and GTF formats, extracting gene, transcript,
 *   and exon features with their attributes.
 * Exported functions: see below
 * HISTORY:
 * Created: Mar 27 2026
 *-------------------------------------------------------------------
 */

#ifndef GFFIO_DEFINED
#define GFFIO_DEFINED

#include "utils.h"
#include "array.h"
#include "dict.h"

typedef enum { GFF_UNKNOWN=0, GFF_GFF3, GFF_GTF } GffFormat ;

typedef struct {
  char *seqid ;       // chromosome/contig name
  char *source ;      // annotation source
  char *type ;        // feature type: gene, mRNA, exon, CDS, etc.
  I64   start ;       // 1-based start position
  I64   end ;         // 1-based end position (inclusive)
  float score ;       // score or -1 if '.'
  char  strand ;      // '+', '-', or '.'
  int   phase ;       // 0, 1, 2, or -1 if '.'
  char *geneId ;      // gene_id attribute
  char *transcriptId ;// transcript_id attribute
  char *geneName ;    // gene_name attribute (optional)
  int   exonNumber ;  // exon_number attribute or -1
} GffFeature ;

typedef struct {
  char *transcriptId ;
  char *geneId ;
  char *seqid ;
  char  strand ;
  I64   txStart ;     // transcript start (min of exon starts)
  I64   txEnd ;       // transcript end (max of exon ends)
  int   nExon ;
  Array exons ;       // of GffExon, sorted by start position
} GffTranscript ;

typedef struct {
  I64 start ;
  I64 end ;
} GffExon ;

typedef struct {
  GffFormat format ;
  char     *fileName ;
  int       iSample ;    // index of sample this annotation belongs to
  Array     features ;   // of GffFeature - all parsed features
  Array     transcripts ;// of GffTranscript - assembled transcripts
  DICT     *geneDict ;   // gene_id -> index
  DICT     *txDict ;     // transcript_id -> index
  DICT     *seqDict ;    // seqid (chromosome) -> index
} GffFile ;

GffFile    *gffFileOpen (char *fileName, int iSample) ;
void        gffFileClose (GffFile *gf) ;
int         gffTranscriptCount (GffFile *gf) ;
GffTranscript *gffTranscript (GffFile *gf, int i) ;

// splice junction representation for graph integration
typedef struct {
  char *seqid ;       // chromosome
  I64   donorEnd ;    // end of upstream exon (donor site)
  I64   acceptorStart ; // start of downstream exon (acceptor site)
  char  strand ;
  char *transcriptId ;
  char *geneId ;
  int   iSample ;     // which sample this junction comes from
} SpliceJunction ;

// extract all splice junctions from a GffFile
Array gffExtractSpliceJunctions (GffFile *gf) ; // returns Array of SpliceJunction

// comparison function for sorting exons by start position
int gffExonCompare (const void *a, const void *b) ;

#endif // GFFIO_DEFINED

/******* end of file ********/
