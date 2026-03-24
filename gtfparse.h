/* File: gtfparse.h
 * Author: (syngpanrna extension)
 *-------------------------------------------------------------------
 * Description: GTF/GFF3 transcript annotation parser for syngpanrna.
 *   Parses exon records into per-transcript exon lists, sorted and
 *   validated, ready for spliced-sequence extraction.
 *-------------------------------------------------------------------
 */

#ifndef GTFPARSE_DEFINED
#define GTFPARSE_DEFINED

#include "utils.h"
#include "array.h"

/* ── Single exon interval (0-based, half-open: [start, end) ) ──────────── */
typedef struct {
  I64  start ;   /* 0-based genomic start (inclusive) */
  I64  end ;     /* 0-based genomic end (exclusive)   */
} Exon ;

/* ── One transcript model ───────────────────────────────────────────────── */
typedef struct {
  char  *txId ;       /* transcript_id attribute value (owned)    */
  char  *geneId ;     /* gene_id attribute value (owned)          */
  char  *chrom ;      /* seqname / chromosome (owned)             */
  char   strand ;     /* '+' or '-'                               */
  Array  exons ;      /* Array of Exon, sorted by start coord     */
  I64    splicedLen ; /* sum of (e.end - e.start) over all exons  */
} TranscriptModel ;

/* ── Parse a GTF or GFF3 annotation file ───────────────────────────────── *
 *
 * Returns an Array of TranscriptModel* (caller owns all memory).
 * Only lines whose feature column (col 3) equals featureType are parsed.
 * transcriptTag is the attribute key used to extract transcript IDs
 * (GTF default: "transcript_id"; GFF3 common alternative: "Parent").
 *
 * Exons within each transcript are sorted by start coordinate.
 * Transcripts are sorted by (chrom, txId).
 * Returns NULL and prints an error on failure.
 *
 * NOTE: Array is already a pointer typedef in syng's array.h
 *   (typedef struct ArrayStruct *Array), so we return Array not Array*.
 */
Array gtfParse (const char *fname,
                const char *featureType,    /* NULL -> "exon"          */
                const char *transcriptTag   /* NULL -> "transcript_id" */
                ) ;

/* ── Parse a BED intron database (alternative/supplementary input) ─────── *
 *
 * Each BED interval [start, end] is treated as an intron.  A synthetic
 * two-exon TranscriptModel is created representing the flanking exonic
 * regions on each side (FLANK=128 bp), giving the syncmer iterator enough
 * context to generate junction-spanning syncmers.
 *
 * Strand is read from column 6 if present, otherwise defaults to '+'.
 * Returns an Array of TranscriptModel* or NULL on failure.
 */
Array intronBedParse (const char *fname) ;

/* ── Free a single TranscriptModel and its contents ────────────────────── */
void transcriptModelDestroy (TranscriptModel *tm) ;

/* ── Free the entire array returned by gtfParse / intronBedParse ────────── */
void transcriptModelArrayDestroy (Array models) ;

#endif /* GTFPARSE_DEFINED */

