/*  File: gffio.c
 *  Author: cmarkello
 *  Copyright (C) 2026
 *-------------------------------------------------------------------
 * Description: GFF/GTF annotation file parser implementation.
 *   Parses GFF3 and GTF files, extracts exon features, assembles
 *   them into transcripts, and generates splice junctions for
 *   spliced pangenome graph construction.
 * HISTORY:
 * Created: Mar 27 2026
 *-------------------------------------------------------------------
 */

#include "gffio.h"
#include <ctype.h>

static GffFormat detectFormat (char *fileName)
{
  int len = strlen (fileName) ;
  if (len > 4 && !strcmp (fileName + len - 4, ".gtf")) return GFF_GTF ;
  if (len > 7 && !strcmp (fileName + len - 7, ".gtf.gz")) return GFF_GTF ;
  if (len > 5 && !strcmp (fileName + len - 5, ".gff3")) return GFF_GFF3 ;
  if (len > 4 && !strcmp (fileName + len - 4, ".gff")) return GFF_GFF3 ;
  if (len > 7 && !strcmp (fileName + len - 7, ".gff.gz")) return GFF_GFF3 ;
  if (len > 8 && !strcmp (fileName + len - 8, ".gff3.gz")) return GFF_GFF3 ;
  return GFF_UNKNOWN ;
}

static char *strdupTrim (char *s)
{
  while (*s && isspace(*s)) ++s ;
  if (*s == '"') ++s ;  // strip leading quote
  int len = strlen(s) ;
  while (len > 0 && (isspace(s[len-1]) || s[len-1] == '"' || s[len-1] == ';'))
    --len ;
  char *r = new(len+1, char) ;
  memcpy (r, s, len) ;
  r[len] = 0 ;
  return r ;
}

// parse GTF attributes: key "value"; key "value"; ...
static void parseGtfAttributes (char *attrs, GffFeature *feat)
{
  char *p = attrs ;
  while (p && *p)
    { while (*p && isspace(*p)) ++p ;
      if (!*p) break ;
      char *key = p ;
      while (*p && !isspace(*p)) ++p ;
      if (!*p) break ;
      *p++ = 0 ;
      while (*p && isspace(*p)) ++p ;
      char *val = p ;
      // find end of value - skip to semicolon or end
      if (*val == '"')
        { ++val ; p = strchr(val, '"') ;
          if (p) { *p++ = 0 ; while (*p && (*p == ';' || isspace(*p))) ++p ; }
        }
      else
        { p = strchr(val, ';') ;
          if (p) *p++ = 0 ;
          else p = val + strlen(val) ;
        }
      if (!strcmp(key, "gene_id"))
        feat->geneId = strdupTrim(val) ;
      else if (!strcmp(key, "transcript_id"))
        feat->transcriptId = strdupTrim(val) ;
      else if (!strcmp(key, "gene_name"))
        feat->geneName = strdupTrim(val) ;
      else if (!strcmp(key, "exon_number"))
        feat->exonNumber = atoi(val) ;
    }
}

// parse GFF3 attributes: key=value;key=value;...
static void parseGff3Attributes (char *attrs, GffFeature *feat)
{
  char *p = attrs ;
  while (p && *p)
    { char *key = p ;
      char *eq = strchr(p, '=') ;
      if (!eq) break ;
      *eq = 0 ;
      char *val = eq + 1 ;
      p = strchr(val, ';') ;
      if (p) *p++ = 0 ;
      else p = val + strlen(val) ;

      if (!strcmp(key, "gene_id") || !strcmp(key, "gene_ID"))
        feat->geneId = strdupTrim(val) ;
      else if (!strcmp(key, "transcript_id") || !strcmp(key, "transcript_ID"))
        feat->transcriptId = strdupTrim(val) ;
      else if (!strcmp(key, "ID"))
        { // ID on a transcript line is the transcript_id; on an exon it's the exon ID
          if (!feat->transcriptId) feat->transcriptId = strdupTrim(val) ;
        }
      else if (!strcmp(key, "Parent"))
        { // Parent on an exon points to transcript_id; on a transcript points to gene_id
          // We store it in transcriptId for exons (overwritten if transcript_id also present)
          // and in geneId for transcripts
          char *trimmed = strdupTrim(val) ;
          if (!feat->transcriptId) feat->transcriptId = trimmed ;
          else if (!feat->geneId) feat->geneId = trimmed ;
          else free(trimmed) ;
        }
      else if (!strcmp(key, "Name") || !strcmp(key, "gene_name"))
        feat->geneName = strdupTrim(val) ;
      else if (!strcmp(key, "exon_number"))
        feat->exonNumber = atoi(val) ;
    }
}

static void parseLine (char *line, GffFormat format, Array features)
{
  if (!line || line[0] == '#' || line[0] == '\n' || line[0] == '\0')
    return ;

  // GFF/GTF has 9 fields, separated by tabs (standard) or spaces (some GFF3 files)
  // detect delimiter: if no tab found before first space, use space
  char delim = '\t' ;
  if (!strchr(line, '\t'))
    delim = ' ' ;

  char *fields[9] ;
  int nField = 0 ;
  char *p = line ;
  while (nField < 9)
    { fields[nField++] = p ;
      if (nField < 9)
        { p = strchr(p, delim) ;
          if (!p) break ;
          *p++ = 0 ;
          // for space delimiter, skip consecutive spaces
          if (delim == ' ')
            while (*p == ' ') ++p ;
        }
    }
  if (nField < 9) return ;

  // only keep exon features (and optionally gene, mRNA/transcript for metadata)
  bool isExon = !strcmp(fields[2], "exon") ;
  bool isTranscript = !strcmp(fields[2], "transcript") || !strcmp(fields[2], "mRNA") ;
  bool isGene = !strcmp(fields[2], "gene") ;
  if (!isExon && !isTranscript && !isGene) return ;

  GffFeature *feat = arrayp(features, arrayMax(features), GffFeature) ;
  memset(feat, 0, sizeof(GffFeature)) ;

  feat->seqid = strdupTrim(fields[0]) ;
  feat->source = strdupTrim(fields[1]) ;
  feat->type = strdupTrim(fields[2]) ;
  feat->start = atoll(fields[3]) ;
  feat->end = atoll(fields[4]) ;
  feat->score = (fields[5][0] == '.') ? -1.0f : atof(fields[5]) ;
  feat->strand = fields[6][0] ;
  feat->phase = (fields[7][0] == '.') ? -1 : atoi(fields[7]) ;
  feat->exonNumber = -1 ;

  // strip trailing newline from attributes
  int alen = strlen(fields[8]) ;
  while (alen > 0 && (fields[8][alen-1] == '\n' || fields[8][alen-1] == '\r'))
    fields[8][--alen] = 0 ;

  // parse attributes
  char *attrCopy = strdupTrim(fields[8]) ;
  if (format == GFF_GTF)
    parseGtfAttributes(attrCopy, feat) ;
  else
    parseGff3Attributes(attrCopy, feat) ;
  free(attrCopy) ;
}

static void assembleTranscripts (GffFile *gf)
{
  int i ;
  gf->transcripts = arrayCreate(1024, GffTranscript) ;
  gf->txDict = dictCreate(1024) ;
  gf->geneDict = dictCreate(1024) ;
  gf->seqDict = dictCreate(256) ;

  for (i = 0 ; i < arrayMax(gf->features) ; ++i)
    { GffFeature *feat = arrp(gf->features, i, GffFeature) ;
      if (strcmp(feat->type, "exon") != 0) continue ;
      if (!feat->transcriptId) continue ;

      U64 txIdx ;
      if (dictAdd(gf->txDict, feat->transcriptId, &txIdx))
        { // new transcript
          GffTranscript *tx = arrayp(gf->transcripts, txIdx, GffTranscript) ;
          memset(tx, 0, sizeof(GffTranscript)) ;
          tx->transcriptId = strdup(feat->transcriptId) ;
          tx->geneId = feat->geneId ? strdup(feat->geneId) : 0 ;
          tx->seqid = strdup(feat->seqid) ;
          tx->strand = feat->strand ;
          tx->txStart = feat->start ;
          tx->txEnd = feat->end ;
          tx->exons = arrayCreate(16, GffExon) ;
          tx->nExon = 0 ;
          if (feat->geneId) dictAdd(gf->geneDict, feat->geneId, 0) ;
          dictAdd(gf->seqDict, feat->seqid, 0) ;
        }

      GffTranscript *tx = arrp(gf->transcripts, txIdx, GffTranscript) ;
      GffExon *exon = arrayp(tx->exons, tx->nExon, GffExon) ;
      exon->start = feat->start ;
      exon->end = feat->end ;
      tx->nExon++ ;
      if (feat->start < tx->txStart) tx->txStart = feat->start ;
      if (feat->end > tx->txEnd) tx->txEnd = feat->end ;
    }

  // sort exons within each transcript by start position
  for (i = 0 ; i < arrayMax(gf->transcripts) ; ++i)
    { GffTranscript *tx = arrp(gf->transcripts, i, GffTranscript) ;
      if (tx->nExon > 1)
        qsort(arrp(tx->exons, 0, GffExon), tx->nExon, sizeof(GffExon), gffExonCompare) ;
    }
}

int gffExonCompare (const void *a, const void *b)
{
  const GffExon *ea = (const GffExon *)a ;
  const GffExon *eb = (const GffExon *)b ;
  if (ea->start < eb->start) return -1 ;
  if (ea->start > eb->start) return 1 ;
  if (ea->end < eb->end) return -1 ;
  if (ea->end > eb->end) return 1 ;
  return 0 ;
}

GffFile *gffFileOpen (char *fileName, int iSample)
{
  GffFormat format = detectFormat(fileName) ;
  if (format == GFF_UNKNOWN)
    { warn("gffFileOpen: cannot determine format of %s (expected .gff, .gff3, or .gtf)", fileName) ;
      format = GFF_GTF ; // default guess
    }

  FILE *f = fzopen(fileName, "r") ;
  if (!f)
    { warn("gffFileOpen: cannot open %s", fileName) ;
      return 0 ;
    }

  GffFile *gf = new0(1, GffFile) ;
  gf->format = format ;
  gf->fileName = strdup(fileName) ;
  gf->iSample = iSample ;
  gf->features = arrayCreate(65536, GffFeature) ;

  char line[65536] ;
  while (fgets(line, sizeof(line), f))
    { char *lineCopy = new(strlen(line)+1, char) ;
      strcpy(lineCopy, line) ;
      parseLine(lineCopy, format, gf->features) ;
      newFree(lineCopy, strlen(line)+1, char) ;
    }
  fclose(f) ;

  fprintf(stdout, "  parsed %lld features from %s (%s format)\n",
          (long long)arrayMax(gf->features), fileName,
          format == GFF_GTF ? "GTF" : "GFF3") ;

  assembleTranscripts(gf) ;

  fprintf(stdout, "  assembled %lld transcripts from %lld genes on %lld sequences\n",
          (long long)arrayMax(gf->transcripts),
          (long long)dictMax(gf->geneDict),
          (long long)dictMax(gf->seqDict)) ;

  return gf ;
}

void gffFileClose (GffFile *gf)
{
  if (!gf) return ;
  int i ;
  for (i = 0 ; i < arrayMax(gf->features) ; ++i)
    { GffFeature *f = arrp(gf->features, i, GffFeature) ;
      if (f->seqid) free(f->seqid) ;
      if (f->source) free(f->source) ;
      if (f->type) free(f->type) ;
      if (f->geneId) free(f->geneId) ;
      if (f->transcriptId) free(f->transcriptId) ;
      if (f->geneName) free(f->geneName) ;
    }
  arrayDestroy(gf->features) ;
  for (i = 0 ; i < arrayMax(gf->transcripts) ; ++i)
    { GffTranscript *tx = arrp(gf->transcripts, i, GffTranscript) ;
      if (tx->transcriptId) free(tx->transcriptId) ;
      if (tx->geneId) free(tx->geneId) ;
      if (tx->seqid) free(tx->seqid) ;
      if (tx->exons) arrayDestroy(tx->exons) ;
    }
  arrayDestroy(gf->transcripts) ;
  if (gf->geneDict) dictDestroy(gf->geneDict) ;
  if (gf->txDict) dictDestroy(gf->txDict) ;
  if (gf->seqDict) dictDestroy(gf->seqDict) ;
  free(gf->fileName) ;
  newFree(gf, 1, GffFile) ;
}

int gffTranscriptCount (GffFile *gf)
{ return arrayMax(gf->transcripts) ; }

GffTranscript *gffTranscript (GffFile *gf, int i)
{ return arrp(gf->transcripts, i, GffTranscript) ; }

Array gffExtractSpliceJunctions (GffFile *gf)
{
  Array junctions = arrayCreate(4096, SpliceJunction) ;
  int i, j ;

  for (i = 0 ; i < arrayMax(gf->transcripts) ; ++i)
    { GffTranscript *tx = arrp(gf->transcripts, i, GffTranscript) ;
      for (j = 0 ; j < tx->nExon - 1 ; ++j)
        { GffExon *e1 = arrp(tx->exons, j, GffExon) ;
          GffExon *e2 = arrp(tx->exons, j+1, GffExon) ;
          SpliceJunction *sj = arrayp(junctions, arrayMax(junctions), SpliceJunction) ;
          sj->seqid = tx->seqid ;         // shared pointer, not owned
          sj->donorEnd = e1->end ;
          sj->acceptorStart = e2->start ;
          sj->strand = tx->strand ;
          sj->transcriptId = tx->transcriptId ; // shared pointer
          sj->geneId = tx->geneId ;             // shared pointer
          sj->iSample = gf->iSample ;
        }
    }

  return junctions ;
}

/******* end of file ********/
