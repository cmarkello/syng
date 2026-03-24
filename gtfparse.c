/* File: gtfparse.c
 * Author: Charles Markello
 *-------------------------------------------------------------------
 * Description: GTF/GFF3 transcript annotation parser.
 *   Supports both GTF (tab-separated, attribute format key "value";)
 *   and GFF3 (tab-separated, attribute format key=value;).
 *   Only exon (or user-specified feature) lines are processed.
 *-------------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "gtfparse.h"
#include "dict.h"   /* for string interning of txId / geneId */

/* ── internal helpers ─────────────────────────────────────────────────── */

static int exonCmp (const void *a, const void *b)
{
  const Exon *ea = (const Exon *)a ;
  const Exon *eb = (const Exon *)b ;
  if (ea->start < eb->start) return -1 ;
  if (ea->start > eb->start) return  1 ;
  return 0 ;
}

static int txModelCmp (const void *a, const void *b)
{
  const TranscriptModel *ta = *(const TranscriptModel **)a ;
  const TranscriptModel *tb = *(const TranscriptModel **)b ;
  int c = strcmp (ta->chrom, tb->chrom) ;
  if (c) return c ;
  return strcmp (ta->txId, tb->txId) ;
}

/* Strip leading/trailing whitespace and surrounding quotes from a string.
   Returns a malloc'd copy. */
static char *stripQuotes (const char *s)
{
  while (*s == ' ' || *s == '\t') ++s ;
  if (*s == '"') ++s ;
  size_t len = strlen (s) ;
  while (len > 0 && (s[len-1] == '"' || s[len-1] == ';' ||
                     s[len-1] == ' ' || s[len-1] == '\t' || s[len-1] == '\n'))
    --len ;
  char *out = (char *)malloc (len + 1) ;
  memcpy (out, s, len) ;
  out[len] = '\0' ;
  return out ;
}

/* Extract attribute value from a GTF attribute string (key "value"; format).
   Returns malloc'd string or NULL if key not found. */
static char *gtfAttr (const char *attrs, const char *key)
{
  size_t klen = strlen (key) ;
  const char *p = attrs ;
  while (*p)
    {
      while (*p == ' ' || *p == '\t') ++p ;  /* skip leading space */
      if (strncmp (p, key, klen) == 0 && (p[klen] == ' ' || p[klen] == '\t'))
        {
          p += klen ;
          while (*p == ' ' || *p == '\t') ++p ;
          return stripQuotes (p) ;
        }
      /* skip to next attribute (past next ';') */
      while (*p && *p != ';') ++p ;
      if (*p == ';') ++p ;
    }
  return NULL ;
}

/* Extract attribute value from a GFF3 attribute string (key=value; format).
   Returns malloc'd string or NULL if key not found. */
static char *gff3Attr (const char *attrs, const char *key)
{
  size_t klen = strlen (key) ;
  const char *p = attrs ;
  while (*p)
    {
      while (*p == ';') ++p ;
      if (strncmp (p, key, klen) == 0 && p[klen] == '=')
        {
          p += klen + 1 ;
          return stripQuotes (p) ;
        }
      while (*p && *p != ';') ++p ;
    }
  return NULL ;
}

/* ── Dict-backed transcript model table ──────────────────────────────────
   Maps txId string → index into models array. */

typedef struct {
  DICT  *txDict ;   /* txId → index in models */
  Array models ;   /* Array of TranscriptModel* */
} TxTable ;

static TxTable *txTableCreate (void)
{
  TxTable *tt = (TxTable *)malloc (sizeof(TxTable)) ;
  tt->txDict  = dictCreate (4096) ;
  tt->models  = arrayCreate (4096, TranscriptModel*) ;
  return tt ;
}

/* Get or create a TranscriptModel for this txId.  geneId/chrom/strand are
   set only when the model is first created; subsequent calls are no-ops for
   those fields (all exons of a transcript should agree). */
static TranscriptModel *txTableGet (TxTable      *tt,
                                    char   *txId,
                                    const char   *geneId,
                                    const char   *chrom,
                                    char          strand)
{
  U64 idx ;
  if (!dictFind (tt->txDict, txId, &idx))
    {
      /* new transcript */
      TranscriptModel *tm = (TranscriptModel *)calloc (1, sizeof(TranscriptModel)) ;
      tm->txId    = strdup (txId) ;
      tm->geneId  = geneId ? strdup (geneId) : strdup ("") ;
      tm->chrom   = strdup (chrom) ;
      tm->strand  = strand ;
      tm->exons   = arrayCreate (8, Exon) ;
      tm->splicedLen = 0 ;

      idx = arrayMax (tt->models) ;
      array(tt->models, idx, TranscriptModel*) = tm ;
      dictAdd (tt->txDict, txId, &idx) ;
    }
  return arr(tt->models, idx, TranscriptModel*) ;
}

/* ── Main parser ──────────────────────────────────────────────────────── */

Array gtfParse (const char *fname,
                 const char *featureType,
                 const char *transcriptTag)
{
  if (!featureType)    featureType   = "exon" ;
  if (!transcriptTag)  transcriptTag = "transcript_id" ;

  FILE *f = fopen (fname, "r") ;
  if (!f) { fprintf (stderr, "gtfParse: cannot open %s\n", fname) ; return NULL ; }

  /* Detect format: GFF3 uses key=value attributes, GTF uses key "value" */
  bool isGFF3 = (strstr (fname, ".gff") || strstr (fname, ".gff3")) ? true : false ;
  /* Also detect by peeking at the first ##gff-version pragma */

  TxTable *tt  = txTableCreate () ;
  char    *line = NULL ;
  size_t   cap  = 0 ;
  I64      lineNo = 0 ;
  I64      nExon  = 0 ;

  while (getline (&line, &cap, f) > 0)
    {
      ++lineNo ;
      if (line[0] == '#')
        {
          /* GFF3 version pragma overrides extension-based detection */
          if (strncmp (line, "##gff-version", 13) == 0) isGFF3 = true ;
          continue ;
        }

      /* Split into 9 tab-separated fields */
      char *fields[9] ;
      int   nf = 0 ;
      char *p  = line ;
      while (nf < 9)
        {
          fields[nf++] = p ;
          char *tab = strchr (p, '\t') ;
          if (!tab) break ;
          *tab = '\0' ;
          p = tab + 1 ;
        }
      if (nf < 9) continue ; /* skip malformed lines */

      /* col 3 (0-indexed col 2): feature type */
      if (strcmp (fields[2], featureType) != 0) continue ;

      char *chrom  = fields[0] ;
      /* col 4,5: start and end (1-based, inclusive in GTF/GFF3) */
      I64   start  = atoll (fields[3]) - 1 ;  /* convert to 0-based */
      I64   end    = atoll (fields[4]) ;       /* half-open end */
      char  strand = fields[6][0] ;
      if (strand != '+' && strand != '-') strand = '+' ;
      char *attrs  = fields[8] ;
      /* strip trailing newline from attrs */
      size_t alen = strlen (attrs) ;
      while (alen > 0 && (attrs[alen-1] == '\n' || attrs[alen-1] == '\r')) attrs[--alen] = '\0' ;

      char *txId = isGFF3 ? gff3Attr (attrs, transcriptTag)
                           : gtfAttr  (attrs, transcriptTag) ;
      if (!txId)
        {
          fprintf (stderr, "gtfParse line %lld: no %s attribute, skipping\n",
                   lineNo, transcriptTag) ;
          continue ;
        }

      /* Try both gene_id (GTF) and gene_id / Parent (GFF3) */
      char *geneId = isGFF3 ? gff3Attr (attrs, "gene_id")
                            : gtfAttr  (attrs, "gene_id") ;
      if (!geneId && isGFF3)
        geneId = gff3Attr (attrs, "Parent") ;

      TranscriptModel *tm = txTableGet (tt, txId, geneId, chrom, strand) ;

      Exon *e  = arrayp (tm->exons, arrayMax(tm->exons), Exon) ;
      e->start = start ;
      e->end   = end ;
      tm->splicedLen += (end - start) ;
      ++nExon ;

      free (txId) ;
      if (geneId) free (geneId) ;
    }
  free (line) ;
  fclose (f) ;

  /* Sort exons within each transcript and compute final splicedLen */
  for (I64 i = 0 ; i < arrayMax(tt->models) ; ++i)
    {
      TranscriptModel *tm = arr(tt->models, i, TranscriptModel*) ;
      arraySort (tm->exons, exonCmp) ;
      /* Recompute splicedLen from sorted exons (deduplication safety) */
      tm->splicedLen = 0 ;
      for (int j = 0 ; j < arrayMax(tm->exons) ; ++j)
        {
          Exon *e = arrp(tm->exons, j, Exon) ;
          tm->splicedLen += e->end - e->start ;
        }
    }

  /* Sort transcript models by (chrom, txId) */
  qsort (arrp(tt->models, 0, TranscriptModel*),
         arrayMax(tt->models),
         sizeof(TranscriptModel*),
         txModelCmp) ;

  fprintf (stdout, "gtfParse: %s  %lld transcripts  %lld exon records\n",
           fname, arrayMax(tt->models), nExon) ;

  Array result = tt->models ;
  dictDestroy (tt->txDict) ;
  free (tt) ;
  return result ;
}

/* ── BED intron database ──────────────────────────────────────────────── */

Array intronBedParse (const char *fname)
{
  FILE *f = fopen (fname, "r") ;
  if (!f) { fprintf (stderr, "intronBedParse: cannot open %s\n", fname) ; return NULL ; }

  Array models = arrayCreate (1024, TranscriptModel*) ;
  char  *line   = NULL ;
  size_t cap    = 0 ;
  I64    lineNo = 0 ;

  while (getline (&line, &cap, f) > 0)
    {
      ++lineNo ;
      if (line[0] == '#' || line[0] == '\n') continue ;

      char *fields[7] = { NULL } ;
      int   nf = 0 ;
      char *p  = line ;
      while (nf < 7)
        {
          fields[nf++] = p ;
          char *tab = strchr (p, '\t') ;
          if (!tab) break ;
          *tab = '\0' ;
          p = tab + 1 ;
        }
      if (nf < 3) continue ;

      char *chrom  = fields[0] ;
      I64   iStart = atoll (fields[1]) ; /* intron start (0-based) */
      I64   iEnd   = atoll (fields[2]) ; /* intron end   (0-based, half-open) */
      char  strand = (nf >= 6 && fields[5]) ? fields[5][0] : '+' ;
      if (strand != '+' && strand != '-') strand = '+' ;

      /* Build a synthetic transcript with two 1-base pseudo-exons flanking
         the intron boundary.  This gives the syncmer iterator enough sequence
         context (w+k bases on each side) to generate junction syncmers.
         We extend by FLANK bases to ensure at least one syncmer is generated
         spanning the junction. */
      const I64 FLANK = 128 ; /* must be > w+k = 63 by default */

      TranscriptModel *tm = (TranscriptModel *)calloc (1, sizeof(TranscriptModel)) ;
      char txIdBuf[256] ;
      snprintf (txIdBuf, sizeof(txIdBuf), "intron:%s:%lld-%lld:%c",
                chrom, iStart, iEnd, strand) ;
      tm->txId   = strdup (txIdBuf) ;
      tm->geneId = strdup ("intron_db") ;
      tm->chrom  = strdup (chrom) ;
      tm->strand = strand ;
      tm->exons  = arrayCreate (2, Exon) ;

      /* donor exon: the region just before the intron */
      Exon *e1   = arrayp (tm->exons, 0, Exon) ;
      e1->start  = (iStart >= FLANK) ? iStart - FLANK : 0 ;
      e1->end    = iStart ;

      /* acceptor exon: the region just after the intron */
      Exon *e2   = arrayp (tm->exons, 1, Exon) ;
      e2->start  = iEnd ;
      e2->end    = iEnd + FLANK ;

      tm->splicedLen = (e1->end - e1->start) + (e2->end - e2->start) ;

      array(models, arrayMax(models), TranscriptModel*) = tm ;
    }
  free (line) ;
  fclose (f) ;

  fprintf (stdout, "intronBedParse: %s  %lld intron junctions\n",
           fname, arrayMax(models)) ;
  return models ;
}

/* ── Cleanup ──────────────────────────────────────────────────────────── */

void transcriptModelDestroy (TranscriptModel *tm)
{
  if (!tm) return ;
  free (tm->txId) ;
  free (tm->geneId) ;
  free (tm->chrom) ;
  arrayDestroy (tm->exons) ;
  free (tm) ;
}

void transcriptModelArrayDestroy (Array models)
{
  if (!models) return ;
  for (I64 i = 0 ; i < arrayMax(models) ; ++i)
    transcriptModelDestroy (arr(models, i, TranscriptModel*)) ;
  arrayDestroy (models) ;
}
