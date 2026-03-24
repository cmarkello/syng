/* File: syngpanrna.c
 * Author: (syngpanrna extension)
 *-------------------------------------------------------------------
 * Description:
 *   Splice-junction injector for syng pangenome graphs.
 *
 * INCLUDE ORDER: syng.h transitively includes syncmerset.h, kmerhash.h,
 *   ONElib.h, seqio.h and utils.h.  We include syng.h last and do NOT
 *   repeat those headers here to avoid the "typedef redefinition" errors
 *   caused by including headers that lack include guards more than once.
 *-------------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <dirent.h>
#include <sys/stat.h>

/* syng.h pulls in: syncmerset.h -> kmerhash.h, ONElib.h, seqio.h,
   seqhash.h, utils.h, array.h, dict.h — include it first so all the
   syng types are defined before anything else we add. */
#include "syng.h"

/* Our own parser — comes after syng.h so Array/Dict are already defined. */
#include "gtfparse.h"
#include "seqhash.h"
#include "dict.h"

/* ── Version ──────────────────────────────────────────────────────────── */
#define SYNGPANRNA_VERSION "0.1"

/* ── SyncPos (mirrors the definition in syng.c / syngpath2gbwt.c) ──────
   Not exported from syng.h, so we define it locally.               ── */
typedef struct {
  I32 sync ;   /* signed: negative = reverse-strand */
  U32 pos ;    /* offset within the sequence        */
} SyncPos ;

/* ── CoordIndex: maps syncmer id -> list of (hapId, seqOffset) ────────── */

typedef struct {
  I32  hapId ;
  I64  seqOffset ;
} SyncOcc ;

typedef struct {
  I64   nSyncs ;
  Array occs ;   /* Array of Array-of-SyncOcc, one per syncmer id */
} CoordIndex ;

static CoordIndex *coordIndexCreate (I64 nSyncs)
{
  CoordIndex *ci = (CoordIndex *)calloc (1, sizeof(CoordIndex)) ;
  ci->nSyncs = nSyncs ;
  /* Store as a flat Array of Array (each inner Array holds SyncOcc items) */
  ci->occs   = arrayCreate (nSyncs + 1, Array) ;
  for (I64 i = 0 ; i <= nSyncs ; ++i)
    array(ci->occs, i, Array) = arrayCreate (2, SyncOcc) ;
  return ci ;
}

static void coordIndexAdd (CoordIndex *ci, I64 syncId, I32 hapId, I64 offset)
{
  if (syncId < 0) syncId = -syncId ;
  if (syncId <= 0 || syncId >= ci->nSyncs) return ;
  Array inner = arr(ci->occs, syncId, Array) ;
  SyncOcc *o  = arrayp (inner, arrayMax(inner), SyncOcc) ;
  o->hapId     = hapId ;
  o->seqOffset = offset ;
}

static void coordIndexDestroy (CoordIndex *ci)
{
  for (I64 i = 0 ; i < ci->nSyncs ; ++i)
    arrayDestroy (arr(ci->occs, i, Array)) ;
  arrayDestroy (ci->occs) ;
  free (ci) ;
}

/* ── HapInfo: metadata parsed from FASTA filenames ───────────────────── */

typedef struct {
  char  *sample ;
  char  *haplotype ;
  char  *chrom ;
  char  *fastaPath ;
  char  *seqName ;
  I32    sourceId ;
} HapInfo ;

static HapInfo *hapInfoFromFilename (const char *dir, const char *fname, I32 sourceId)
{
  const char *SUFFIX = "_hprc_r2_v1.0.1." ;
  char *base = strdup (fname) ;
  char *dot  = strstr (base, ".fa.gz") ;
  if (!dot) { free(base) ; return NULL ; }
  *dot = '\0' ;

  char *mid = strstr (base, SUFFIX) ;
  if (!mid) { free(base) ; return NULL ; }
  *mid = '\0' ;
  char *chrom = mid + strlen(SUFFIX) ;

  char *under = strrchr (base, '_') ;
  if (!under) { free(base) ; return NULL ; }
  *under = '\0' ;
  char *hap    = under + 1 ;
  char *sample = base ;

  if (strcmp(hap,"mat") != 0 && strcmp(hap,"pat") != 0) { free(base) ; return NULL ; }

  HapInfo *hi   = (HapInfo *)calloc (1, sizeof(HapInfo)) ;
  hi->sample    = strdup (sample) ;
  hi->haplotype = strdup (hap) ;
  hi->chrom     = strdup (chrom) ;
  hi->sourceId  = sourceId ;

  size_t plen   = strlen(dir) + 1 + strlen(sample) + 1 + strlen(fname) + 1 ;
  hi->fastaPath = (char *)malloc (plen) ;
  snprintf (hi->fastaPath, plen, "%s/%s/%s", dir, sample, fname) ;

  free (base) ;
  return hi ;
}

static void hapInfoDestroy (HapInfo *hi)
{
  if (!hi) return ;
  free (hi->sample) ; free (hi->haplotype) ; free (hi->chrom) ;
  free (hi->fastaPath) ; free (hi->seqName) ;
  free (hi) ;
}

/* ── Scan fasta_dir for all matching .fa.gz files ─────────────────────── */

static Array collectHaplotypes (const char *fastaDir, const char *targetChrom)
{
  Array haps = arrayCreate (32, HapInfo*) ;
  DIR  *d    = opendir (fastaDir) ;
  if (!d) die ("syngpanrna: cannot open fasta_dir %s", fastaDir) ;

  struct dirent *e ;
  I32 sourceId = 0 ;

  while ((e = readdir(d)) != NULL)
    {
      if (e->d_type != DT_DIR) continue ;
      const char *sample = e->d_name ;
      if (sample[0] == '.') continue ;

      char subdir[4096] ;
      snprintf (subdir, sizeof(subdir), "%s/%s", fastaDir, sample) ;
      DIR *sd = opendir (subdir) ;
      if (!sd) continue ;

      struct dirent *se ;
      while ((se = readdir(sd)) != NULL)
        {
          const char *fn = se->d_name ;
          char suffix[256] ;
          snprintf (suffix, sizeof(suffix), ".%s.fa.gz", targetChrom) ;
          size_t fnlen  = strlen (fn) ;
          size_t suflen = strlen (suffix) ;
          if (fnlen < suflen) continue ;
          if (strcmp (fn + fnlen - suflen, suffix) != 0) continue ;
          if (!strstr (fn, "_hprc_r2_v1.0.1.")) continue ;

          HapInfo *hi = hapInfoFromFilename (fastaDir, fn, ++sourceId) ;
          if (hi)
            array(haps, arrayMax(haps), HapInfo*) = hi ;
        }
      closedir (sd) ;
    }
  closedir (d) ;

  fprintf (stdout, "collectHaplotypes: found %lld haplotype FASTA(s) for chrom %s\n",
           arrayMax(haps), targetChrom) ;
  return haps ;
}

/* ── Extract a region via samtools faidx ────────────────────────────────
   Returns malloc'd char[] in dna2index4Conv encoding, length (end-start).
   Coordinates are 0-based half-open [start, end).                      */

static char *fastaExtractRegion (const char *fastaPath,
                                 const char *seqName,
                                 I64 start, I64 end)
{
  if (end <= start) return NULL ;
  I64  len = end - start ;
  char cmd[4096] ;
  snprintf (cmd, sizeof(cmd),
            "samtools faidx %s '%s:%lld-%lld' 2>/dev/null",
            fastaPath, seqName, start+1, end) ;

  FILE *p = popen (cmd, "r") ;
  if (!p) { fprintf (stderr, "fastaExtractRegion: popen failed\n") ; return NULL ; }

  char  *buf  = (char *)malloc (len + 1) ;
  I64    pos  = 0 ;
  char  *line = NULL ;
  size_t cap  = 0 ;

  while (getline (&line, &cap, p) > 0)
    {
      if (line[0] == '>') continue ;
      for (char *c = line ; *c ; ++c)
        {
          if (*c == '\n' || *c == '\r') continue ;
          if (pos < len)
            buf[pos++] = (char)dna2index4Conv[(unsigned char)*c] ;
        }
    }
  free (line) ;
  pclose (p) ;

  if (pos != len)
    {
      fprintf (stderr,
               "fastaExtractRegion: expected %lld bases, got %lld (%s:%lld-%lld)\n",
               len, pos, seqName, start+1, end) ;
      free (buf) ;
      return NULL ;
    }
  buf[len] = '\0' ;
  return buf ;
}

/* ── Reverse-complement in dna2index4Conv encoding in-place ────────────── */
static void reverseComplementIndex4 (char *seq, I64 len)
{
  static const char comp[5] = { 3, 2, 1, 0, 4 } ;
  I64 i = 0, j = len - 1 ;
  while (i <= j)
    {
      char tmp = comp[(unsigned char)seq[i]] ;
      seq[i]   = comp[(unsigned char)seq[j]] ;
      seq[j]   = tmp ;
      ++i ; --j ;
    }
}

/* ── Lazy-load the sequence name from the .fai sidecar ─────────────────── */
static bool hapInfoEnsureSeqName (HapInfo *hi)
{
  if (hi->seqName) return true ;
  char fai[4096] ;
  snprintf (fai, sizeof(fai), "%s.fai", hi->fastaPath) ;
  FILE *f = fopen (fai, "r") ;
  if (!f) { fprintf (stderr, "hapInfoEnsureSeqName: cannot open %s\n", fai) ; return false ; }
  char  *line = NULL ;
  size_t cap  = 0 ;
  if (getline (&line, &cap, f) > 0)
    {
      char *tab = strchr (line, '\t') ;
      if (tab) *tab = '\0' ;
      size_t l = strlen(line) ;
      while (l > 0 && (line[l-1]=='\n'||line[l-1]=='\r')) line[--l]='\0' ;
      hi->seqName = strdup (line) ;
    }
  free (line) ;
  fclose (f) ;
  return hi->seqName != NULL ;
}

/* ── Build the spliced sequence for one (transcript, haplotype) pair ────── */
static char *buildSplicedSeq (TranscriptModel *tm, HapInfo *hi, I64 *outLen)
{
  if (tm->splicedLen <= 0) return NULL ;
  char *out = (char *)malloc (tm->splicedLen + 1) ;
  I64   pos = 0 ;

  for (int i = 0 ; i < arrayMax(tm->exons) ; ++i)
    {
      Exon *e   = arrp(tm->exons, i, Exon) ;
      I64   len = e->end - e->start ;
      char *seg = fastaExtractRegion (hi->fastaPath, hi->seqName, e->start, e->end) ;
      if (!seg) { free(out) ; return NULL ; }
      memcpy (out + pos, seg, len) ;
      pos += len ;
      free (seg) ;
    }

  if (tm->strand == '-')
    reverseComplementIndex4 (out, tm->splicedLen) ;

  *outLen = tm->splicedLen ;
  return out ;
}

/* ── Process one (transcript, haplotype) pair ───────────────────────────── */
static I64 processTxHap (TranscriptModel *tm,
                          HapInfo         *hi,
                          SyncmerSet      *sms,
                          Seqhash         *sh,
                          OneFile         *ofPath,
                          OneFile         *ofTxPath,
                          I32              txIdx,
                          I32              geneIdx,
                          I32              hapId,
                          bool             doCollapse,
                          DICT            *collapseDict,
                          I64             *nNovel)
{
  if (!hapInfoEnsureSeqName (hi)) return -1 ;

  I64   seqLen ;
  char *seq = buildSplicedSeq (tm, hi, &seqLen) ;
  if (!seq) return -1 ;

  SeqhashIterator *sit = syncmerIterator (sh, seq, (int)seqLen) ;
  Array syncPos = arrayCreate (256, SyncPos) ;
  int   pos ;
  while (syncmerNext (sit, NULL, &pos, NULL))
    {
      I64  syncId = 0 ;
      bool added  = kmerHashAdd (sms->kh, seq + pos, &syncId) ;
      if (added)
        {
          ++(*nNovel) ;
          I64 i = (syncId < 0) ? -syncId : syncId ;
          array(sms->count, i, I64) = 1 ;
          array(sms->thisCount, i, char) = 1 ;
        }
      else
        syncmerCount (sms, (I32)syncId) ;

      SyncPos *sp = arrayp (syncPos, arrayMax(syncPos), SyncPos) ;
      sp->sync = (I32)syncId ;
      sp->pos  = (U32)pos ;
    }
  seqhashIteratorDestroy (sit) ;

  I64 nSync = arrayMax (syncPos) ;

  /* Optional collapse */
  if (doCollapse && collapseDict && nSync > 0)
    {
      char  *key = (char *)malloc (nSync * 14 + 1) ;
      char  *kp  = key ;
      for (I64 j = 0 ; j < nSync ; ++j)
        kp += sprintf (kp, "%d,", arrp(syncPos, j, SyncPos)->sync) ;
      U64 dummy ;
      if (dictFind (collapseDict, key, &dummy))
        { free (key) ; free (seq) ; arrayDestroy (syncPos) ; return 0 ; }
      dictAdd (collapseDict, key, &dummy) ;
      free (key) ;
    }

  if (nSync == 0)
    { free (seq) ; arrayDestroy (syncPos) ; return 0 ; }

  /* Reusable buffer for z/o lines */
  static I64   *zbuf     = NULL ;
  static size_t zbufSize = 0 ;
  if (!zbuf || (size_t)nSync > zbufSize)
    { free(zbuf) ; zbufSize = (size_t)nSync * 2 ; zbuf = (I64 *)malloc(zbufSize * sizeof(I64)) ; }

  OneFile *ofs[2] = { ofPath, ofTxPath } ;
  for (int fi = 0 ; fi < 2 ; ++fi)
    {
      OneFile *of = ofs[fi] ;
      if (!of) continue ;

      oneInt(of,0) = seqLen ; oneInt(of,1) = hi->sourceId ; oneInt(of,2) = txIdx ;
      oneWriteLine (of, 'P', 0, 0) ;

      oneInt(of,0) = txIdx ; oneInt(of,1) = geneIdx ; oneInt(of,2) = hapId ;
      oneWriteLine (of, 'T', 0, 0) ;

      for (I64 j = 0 ; j < nSync ; ++j) zbuf[j] = arrp(syncPos, j, SyncPos)->sync ;
      oneWriteLine (of, 'z', nSync, zbuf) ;

      for (I64 j = 0 ; j < nSync ; ++j) zbuf[j] = arrp(syncPos, j, SyncPos)->pos ;
      oneWriteLine (of, 'o', nSync, zbuf) ;

      SyncPos *sp0    = arrp(syncPos, 0, SyncPos) ;
      if (sp0->pos > 0)
        oneWriteLine (of, 'X', sp0->pos, seq) ;

      SyncPos *spLast = arrp(syncPos, nSync-1, SyncPos) ;
      I64 endOff = spLast->pos + sms->kh->len ;
      if (endOff < seqLen)
        oneWriteLine (of, 'Y', seqLen - endOff, seq + endOff) ;
    }

  free (seq) ;
  arrayDestroy (syncPos) ;
  return nSync ;
}

/* ── Copy all paths from an existing .1path file to an output .1path ──── */
static I64 copyPaths (OneFile *ofIn, OneFile *ofOut)
{
  I64 nCopied = 0 ;
  if (!oneGoto (ofIn, 'P', 1)) return 0 ;
  oneReadLine (ofIn) ;
  while (ofIn->lineType == 'P')
    {
      oneInt(ofOut,0) = oneInt(ofIn,0) ;
      oneInt(ofOut,1) = oneInt(ofIn,1) ;
      oneInt(ofOut,2) = oneInt(ofIn,2) ;
      oneWriteLine (ofOut, 'P', 0, 0) ;
      ++nCopied ;

      while (oneReadLine (ofIn) && ofIn->lineType != 'P')
        {
          switch (ofIn->lineType)
            {
            case 'z': oneWriteLine (ofOut, 'z', oneLen(ofIn), oneIntList(ofIn)) ; break ;
            case 'o': oneWriteLine (ofOut, 'o', oneLen(ofIn), oneIntList(ofIn)) ; break ;
            case 'X': oneWriteLine (ofOut, 'X', oneLen(ofIn), oneDNAchar(ofIn)) ; break ;
            case 'Y': oneWriteLine (ofOut, 'Y', oneLen(ofIn), oneDNAchar(ofIn)) ; break ;
            case 'Z':
              oneInt(ofOut,0) = oneInt(ofIn,0) ; oneInt(ofOut,1) = oneInt(ofIn,1) ;
              oneInt(ofOut,2) = oneInt(ofIn,2) ; oneInt(ofOut,3) = oneInt(ofIn,3) ;
              oneWriteLine (ofOut, 'Z', 0, 0) ;
              break ;
            default: break ;
            }
        }
    }
  return nCopied ;
}

/* ── Build GBWT from a .1path file ──────────────────────────────────────── */
static void buildGBWT (const char *pathFile, const char *gbwtFile,
                       OneSchema *schema, int syncLen)
{
  OneFile *ofIn  = oneFileOpenRead (pathFile, schema, "path", 1) ;
  if (!ofIn) die ("buildGBWT: cannot open %s", pathFile) ;

  OneFile *ofOut = oneFileOpenWriteNew (gbwtFile, schema, "gbwt", true, 1) ;
  if (!ofOut) die ("buildGBWT: cannot open %s for write", gbwtFile) ;
  oneAddProvenance (ofOut, "syngpanrna", SYNGPANRNA_VERSION, "") ;
  oneAddReference  (ofOut, pathFile, 1) ;

  if (oneGoto (ofIn, 'h', 1))
    {
      oneReadLine (ofIn) ;
      oneInt(ofOut,0) = oneInt(ofIn,0) ; oneInt(ofOut,1) = oneInt(ofIn,1) ;
      oneInt(ofOut,2) = oneInt(ofIn,2) ;
      oneWriteLine (ofOut, 'h', 0, 0) ;
    }

  SyngBWT *gbwt  = syngBWTcreate (syncLen, 0) ;
  Array    syncPos = arrayCreate (1024, SyncPos) ;
  I64      totSync = 0 ;
  int      pathCount = 0 ;

  if (!oneGoto (ofIn, 'P', 1)) die ("buildGBWT: no P records in %s", pathFile) ;
  oneReadLine (ofIn) ;

  while (ofIn->lineType == 'P')
    {
      ++pathCount ;
      I64 len      = oneInt(ofIn,0) ;
      I64 source   = oneInt(ofIn,1) ;
      I64 inSource = oneInt(ofIn,2) ;
      I64 nSync    = 0 ;
      SyncPos *sp  = NULL ;
      char *dnaX = NULL, *dnaY = NULL ;
      I64 dnaXlen = 0, dnaYlen = 0 ;

      while (oneReadLine (ofIn) && ofIn->lineType != 'P')
        {
          I64 *iList ;
          switch (ofIn->lineType)
            {
            case 'z':
              nSync = oneLen(ofIn) ; totSync += nSync ;
              arrayp(syncPos, nSync, SyncPos)->pos = 0 ;
              sp    = arrp(syncPos, 0, SyncPos) ;
              iList = oneIntList(ofIn) ;
              for (I64 j = 0 ; j < nSync ; ++j) sp[j].sync = (I32)iList[j] ;
              break ;
            case 'o':
              if (sp && oneLen(ofIn) == nSync)
                { iList = oneIntList(ofIn) ;
                  for (I64 j = 0 ; j < nSync ; ++j) sp[j].pos = (U32)iList[j] ; }
              break ;
            case 'X': dnaXlen = oneLen(ofIn) ; dnaX = oneDNAchar(ofIn) ; break ;
            case 'Y': dnaYlen = oneLen(ofIn) ; dnaY = oneDNAchar(ofIn) ; break ;
            default:  break ;
            }
        }

      if (!nSync || !sp) continue ;

      oneInt(ofOut,0) = len ; oneInt(ofOut,1) = source ; oneInt(ofOut,2) = inSource ;
      oneWriteLine (ofOut, 'P', 0, 0) ;

      SyngBWTpath *sbp = syngBWTpathStartNew (gbwt, sp[0].sync) ;
      U32 j0 = sbp->jLast ;
      oneInt(ofOut,0) = sp[0].sync ; oneInt(ofOut,1) = sp[0].pos ;
      oneInt(ofOut,2) = j0          ; oneInt(ofOut,3) = nSync ;
      oneWriteLine (ofOut, 'Z', 0, 0) ;

      if (dnaXlen) oneWriteLine (ofOut, 'X', dnaXlen, dnaX) ;
      if (dnaYlen) oneWriteLine (ofOut, 'Y', dnaYlen, dnaY) ;

      for (I64 k = 1 ; k < nSync ; ++k)
        syngBWTpathAdd (sbp, sp[k].sync, sp[k].pos - sp[k-1].pos) ;
      syngBWTpathFinish (sbp) ;

      sbp = syngBWTpathStartNew (gbwt, -sp[nSync-1].sync) ;
      for (I64 k = nSync-2 ; k >= 0 ; --k)
        syngBWTpathAdd (sbp, -sp[k].sync, sp[k+1].pos - sp[k].pos) ;
      syngBWTpathFinish (sbp) ;
    }

  oneFileClose (ofIn) ;
  arrayDestroy (syncPos) ;
  fprintf (stdout, "buildGBWT: %d paths, %lld syncmers\n", pathCount, totSync) ;

  syngBWTwrite (ofOut, gbwt) ;
  oneFileClose (ofOut) ;
  syngBWTdestroy (gbwt) ;
  fprintf (stdout, "buildGBWT: wrote %s\n", gbwtFile) ;
}

/* ── Write tx_info.tsv ────────────────────────────────────────────────── */
static void writeTxInfo (const char *fname, Array txModels, Array haps)
{
  FILE *f = fopen (fname, "w") ;
  if (!f) { fprintf (stderr, "writeTxInfo: cannot open %s\n", fname) ; return ; }
  fprintf (f, "path_index\ttx_id\tgene_id\tsample\thaplotype\tchrom\tstrand\tspliced_len\n") ;
  I64 pidx = 1 ;
  for (I64 ti = 0 ; ti < arrayMax(txModels) ; ++ti)
    {
      TranscriptModel *tm  = arr(txModels, ti, TranscriptModel*) ;
      for (I64 hi = 0 ; hi < arrayMax(haps) ; ++hi)
        {
          HapInfo *hap = arr(haps, hi, HapInfo*) ;
          if (strcmp (hap->chrom, tm->chrom) != 0) continue ;
          fprintf (f, "%lld\t%s\t%s\t%s\t%s\t%s\t%c\t%lld\n",
                   pidx++, tm->txId, tm->geneId,
                   hap->sample, hap->haplotype, hap->chrom,
                   tm->strand, tm->splicedLen) ;
        }
    }
  fclose (f) ;
  fprintf (stdout, "writeTxInfo: wrote %s\n", fname) ;
}

/* ═══════════════════════════════════════════════════════════════════════ */
/*                              main                                      */
/* ═══════════════════════════════════════════════════════════════════════ */

static char usage[] =
  "Usage: syngpanrna -K <.1khash> -p <.1path> -f <fasta_dir> -o <prefix>\n"
  "                  [options] <annotation.[gtf|gff3]>\n"
  "\n"
  "Required:\n"
  "  -K <.1khash>        existing syncmer hash table from syng\n"
  "  -p <.1path>         existing path file from syng (-outputEnds)\n"
  "  -f <fasta_dir>      directory with per-sample FASTA subdirs\n"
  "  -o <prefix>         output file prefix\n"
  "  <annotation>        GTF or GFF3 transcript annotation\n"
  "\n"
  "Options:\n"
  "  --feature-type <s>   GTF/GFF3 feature to parse [exon]\n"
  "  --transcript-tag <s> attribute key for transcript ID [transcript_id]\n"
  "  --introns <bed>      additional intron BED database\n"
  "  --chrom <id>         restrict to this chromosome\n"
  "  --collapse           collapse identical transcript paths across haplotypes\n"
  "  --write-info <tsv>   write per-path metadata TSV\n" ;

int main (int argc, char *argv[])
{
  timeUpdate (0) ;
  storeCommandLine (argc, argv) ;

  char *khashFile   = NULL ;
  char *pathFile    = NULL ;
  char *fastaDir    = NULL ;
  char *outPrefix   = NULL ;
  char *annotFile   = NULL ;
  char *intronBed   = NULL ;
  char *chromId     = NULL ;
  char *featureType = NULL ;
  char *txTag       = NULL ;
  char *infoFile    = NULL ;
  bool  doCollapse  = false ;

  argc-- ; ++argv ;
  if (!argc) { fprintf (stderr, "%s", usage) ; exit (0) ; }

  while (argc > 0)
    {
      if      (!strcmp(*argv,"-K") && argc>1)            { khashFile   = argv[1] ; argc-=2 ; argv+=2 ; }
      else if (!strcmp(*argv,"-p") && argc>1)            { pathFile    = argv[1] ; argc-=2 ; argv+=2 ; }
      else if (!strcmp(*argv,"-f") && argc>1)            { fastaDir    = argv[1] ; argc-=2 ; argv+=2 ; }
      else if (!strcmp(*argv,"-o") && argc>1)            { outPrefix   = argv[1] ; argc-=2 ; argv+=2 ; }
      else if (!strcmp(*argv,"--feature-type") && argc>1){ featureType = argv[1] ; argc-=2 ; argv+=2 ; }
      else if (!strcmp(*argv,"--transcript-tag") && argc>1){ txTag     = argv[1] ; argc-=2 ; argv+=2 ; }
      else if (!strcmp(*argv,"--introns") && argc>1)     { intronBed   = argv[1] ; argc-=2 ; argv+=2 ; }
      else if (!strcmp(*argv,"--chrom") && argc>1)       { chromId     = argv[1] ; argc-=2 ; argv+=2 ; }
      else if (!strcmp(*argv,"--collapse"))               { doCollapse  = true ; --argc ; ++argv ; }
      else if (!strcmp(*argv,"--write-info") && argc>1)  { infoFile    = argv[1] ; argc-=2 ; argv+=2 ; }
      else if (**argv != '-')                             { annotFile   = *argv  ; --argc ; ++argv ; }
      else die ("syngpanrna: unknown argument %s\n%s", *argv, usage) ;
    }

  if (!khashFile || !pathFile || !fastaDir || !outPrefix || !annotFile)
    die ("syngpanrna: missing required arguments\n%s", usage) ;

  /* ── 1. Read syncmer set ── */
  fprintf (stdout, "Reading syncmer hash: %s\n", khashFile) ;
  SyncmerSet *sms = syncmerSetRead (khashFile) ;
  if (!sms) die ("syngpanrna: cannot read %s", khashFile) ;
  SyncmerParams params = sms->params ;
  fprintf (stdout, "  k=%d  w=%d  seed=%d  nSyncmers=%lld\n",
           params.k, params.w, params.seed, kmerHashMax(sms->kh)) ;
  int syncLen = params.w + params.k ;
  Seqhash *sh = seqhashCreate (params.k, params.w + 1, params.seed) ;

  /* ── 2. Collect haplotypes ── */
  if (!chromId) chromId = "" ;
  Array haps = collectHaplotypes (fastaDir, chromId) ;
  if (arrayMax(haps) == 0)
    die ("syngpanrna: no haplotype FASTAs found in %s for chrom '%s'", fastaDir, chromId) ;

  /* ── 3. Parse annotation ── */
  fprintf (stdout, "Parsing annotation: %s\n", annotFile) ;
  Array txModels = gtfParse (annotFile, featureType, txTag) ;
  if (!txModels || arrayMax(txModels) == 0)
    die ("syngpanrna: no transcripts parsed from %s", annotFile) ;

  if (intronBed)
    {
      fprintf (stdout, "Parsing intron BED: %s\n", intronBed) ;
      Array introns = intronBedParse (intronBed) ;
      if (introns)
        {
          for (I64 i = 0 ; i < arrayMax(introns) ; ++i)
            array(txModels, arrayMax(txModels), TranscriptModel*) =
              arr(introns, i, TranscriptModel*) ;
          arrayDestroy (introns) ;
        }
    }

  /* ── 4. Build string tables ── */
  DICT *txDict   = dictCreate (arrayMax(txModels) * 2) ;
  DICT *geneDict = dictCreate (arrayMax(txModels) * 2) ;
  for (I64 i = 0 ; i < arrayMax(txModels) ; ++i)
    {
      TranscriptModel *tm = arr(txModels, i, TranscriptModel*) ;
      U64 dummy ;
      dictAdd (txDict,   tm->txId,   &dummy) ;
      dictAdd (geneDict, tm->geneId, &dummy) ;
    }

  /* ── 5. Open output files ── */
  OneSchema *schema = oneSchemaCreateFromText (syngSchemaText) ;

  char splicedPathFile[4096], splicedKhashFile[4096] ;
  char splicedGbwtFile[4096], txPathFile[4096] ;
  snprintf (splicedPathFile,  sizeof(splicedPathFile),  "%s.spliced.1path",  outPrefix) ;
  snprintf (splicedKhashFile, sizeof(splicedKhashFile), "%s.spliced.1khash", outPrefix) ;
  snprintf (splicedGbwtFile,  sizeof(splicedGbwtFile),  "%s.spliced.1gbwt",  outPrefix) ;
  snprintf (txPathFile,       sizeof(txPathFile),        "%s.tx.1path",       outPrefix) ;

  OneFile *ofPath   = oneFileOpenWriteNew (splicedPathFile, schema, "path", true, 1) ;
  if (!ofPath)   die ("syngpanrna: cannot open %s for write", splicedPathFile) ;
  OneFile *ofTxPath = oneFileOpenWriteNew (txPathFile, schema, "path", true, 1) ;
  if (!ofTxPath) die ("syngpanrna: cannot open %s for write", txPathFile) ;

  oneAddProvenance (ofPath,   "syngpanrna", SYNGPANRNA_VERSION, getCommandLine()) ;
  oneAddProvenance (ofTxPath, "syngpanrna", SYNGPANRNA_VERSION, getCommandLine()) ;
  syncmerParamsWrite (ofPath, params) ; syncmerParamsWrite (ofTxPath, params) ;

  for (I64 i = 1 ; i <= dictMax(txDict) ; ++i)
    { char *s = dictName (txDict, i) ;
      oneWriteLine (ofPath, 'N', strlen(s), s) ;
      oneWriteLine (ofTxPath, 'N', strlen(s), s) ; }
  for (I64 i = 1 ; i <= dictMax(geneDict) ; ++i)
    { char *s = dictName (geneDict, i) ;
      oneWriteLine (ofPath, 'G', strlen(s), s) ;
      oneWriteLine (ofTxPath, 'G', strlen(s), s) ; }

  /* ── 6. Copy genomic paths ── */
  fprintf (stdout, "Copying genomic paths from %s\n", pathFile) ;
  OneFile *ofPathIn = oneFileOpenRead (pathFile, schema, "path", 1) ;
  if (!ofPathIn) die ("syngpanrna: cannot open %s", pathFile) ;
  I64 nGenomic = copyPaths (ofPathIn, ofPath) ;
  oneFileClose (ofPathIn) ;
  fprintf (stdout, "  copied %lld genomic paths\n", nGenomic) ;

  /* ── 7. Generate transcript paths ── */
  fprintf (stdout, "Processing %lld transcripts x %lld haplotypes\n",
           arrayMax(txModels), arrayMax(haps)) ;

  DICT *collapseDict = doCollapse ? dictCreate (65536) : NULL ;
  I64  nTxPaths = 0, nNovel = 0, nFailed = 0 ;

  for (I64 ti = 0 ; ti < arrayMax(txModels) ; ++ti)
    {
      TranscriptModel *tm = arr(txModels, ti, TranscriptModel*) ;
      U64 txIdx, geneIdx ;
      dictFind (txDict,   tm->txId,   &txIdx) ;
      dictFind (geneDict, tm->geneId, &geneIdx) ;

      for (I64 hi = 0 ; hi < arrayMax(haps) ; ++hi)
        {
          HapInfo *hap = arr(haps, hi, HapInfo*) ;
          if (strcmp(hap->chrom, tm->chrom) != 0) continue ;

          I64 nSync = processTxHap (tm, hap, sms, sh,
                                    ofPath, ofTxPath,
                                    (I32)txIdx, (I32)geneIdx, (I32)hi,
                                    doCollapse, collapseDict, &nNovel) ;
          if (nSync < 0) { ++nFailed ; continue ; }
          if (nSync > 0) ++nTxPaths ;
        }

      if ((ti % 1000) == 0 && ti > 0)
        fprintf (stdout, "  %lld/%lld transcripts, %lld novel syncmers\n",
                 ti, arrayMax(txModels), nNovel) ;
    }

  oneFileClose (ofPath) ;
  oneFileClose (ofTxPath) ;
  fprintf (stdout, "Transcript paths: %lld written, %lld failed, %lld novel syncmers\n",
           nTxPaths, nFailed, nNovel) ;

  /* ── 8. Write augmented .1khash ── */
  fprintf (stdout, "Writing augmented syncmer hash: %s\n", splicedKhashFile) ;
  OneFile *ofKhash = oneFileOpenWriteNew (splicedKhashFile, schema, "khash", true, 1) ;
  if (!ofKhash) die ("syngpanrna: cannot open %s for write", splicedKhashFile) ;
  oneAddProvenance (ofKhash, "syngpanrna", SYNGPANRNA_VERSION, getCommandLine()) ;
  syncmerSetWrite (sms, ofKhash) ;
  oneFileClose (ofKhash) ;

  /* ── 9. Build GBWT ── */
  fprintf (stdout, "Building GBWT: %s\n", splicedGbwtFile) ;
  buildGBWT (splicedPathFile, splicedGbwtFile, schema, syncLen) ;

  /* ── 10. tx_info.tsv ── */
  if (infoFile) writeTxInfo (infoFile, txModels, haps) ;

  /* ── Cleanup ── */
  syncmerSetDestroy (sms) ;
  seqhashDestroy (sh) ;
  transcriptModelArrayDestroy (txModels) ;
  for (I64 i = 0 ; i < arrayMax(haps) ; ++i)
    hapInfoDestroy (arr(haps, i, HapInfo*)) ;
  arrayDestroy (haps) ;
  dictDestroy (txDict) ; dictDestroy (geneDict) ;
  if (collapseDict) dictDestroy (collapseDict) ;

  fprintf (stdout, "\nOutputs:\n  %s\n  %s\n  %s\n  %s\n",
           splicedKhashFile, splicedPathFile, splicedGbwtFile, txPathFile) ;
  if (infoFile) fprintf (stdout, "  %s\n", infoFile) ;

  fprintf (stdout, "total: ") ; timeTotal (stdout) ;
  return 0 ;
}

