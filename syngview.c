/*  File: syngview.c
 *  Author: Charles Markello
 *  Copyright (C) 2026
 *-------------------------------------------------------------------
 * Description: exports syng pangenome graph as rGFA and/or BED
 * HISTORY:
 * Created: Mar 31 2026
 *-------------------------------------------------------------------
 */

#include <pthread.h>

#include "seqio.h"
#include "seqhash.h"
#include "syng.h"

extern int pathCount ;

/**************** reference coordinate mapping ********************/

typedef struct {
  char *sampleId ;     // e.g. "CHM13"
  char *seqName ;      // e.g. "chr19"
  I64   seqLen ;       // length of the reference sequence
  SyncmerSet *sms ;    // syncmer set (shared with main)
  KmerHash   *kh ;     // kmer hash (shared with main)
  I64  *nodeRefPos ;   // 1-based node index -> reference position (-1 if unmapped)
  bool *nodeOnRef ;    // 1-based node index -> true if on reference
  I64   nMapped ;      // number of nodes mapped to reference
} RefMap ;

static RefMap *refMapCreate (char *sampleId, char *refFasta, SyncmerSet *sms, int nThread)
{
  RefMap *rm = new0 (1, RefMap) ;
  rm->sms = sms ;
  rm->kh = sms->kh ;

  // extract sample id (everything before first '_' or the whole string)
  rm->sampleId = strdup (sampleId) ;

  // read the reference fasta to find syncmers and their positions
  SeqIO *sio = seqIOopenRead (refFasta, dna2index4Conv, 0) ;
  if (!sio) die ("failed to open reference fasta %s", refFasta) ;

  I64 nodeMax = kmerHashMax (rm->kh) + 1 ;
  rm->nodeRefPos = new (nodeMax, I64) ;
  rm->nodeOnRef = new0 (nodeMax, bool) ;
  for (I64 i = 0 ; i < nodeMax ; ++i) rm->nodeRefPos[i] = -1 ;

  Seqhash *sh = seqhashCreate (sms->params.k, sms->params.w + 1, sms->params.seed) ;
  U64 *uBuf = new (rm->kh->plen, U64) ;

  rm->nMapped = 0 ;
  while (seqIOread (sio))
    { if (!rm->seqName) rm->seqName = strdup (sqioId(sio)) ;
      rm->seqLen = sio->seqLen ;
      char *seq = sqioSeq(sio) ;
      int pos ;
      SeqhashIterator *sit = syncmerIterator (sh, seq, sio->seqLen) ;
      while (syncmerNext (sit, 0, &pos, 0))
        { I64 sync = 0 ;
          kmerHashFindThreadSafe (rm->kh, seq + pos, &sync, uBuf) ;
          if (sync && (sync > 2 || sync < -2))
            { I64 absSync = sync > 0 ? sync : -sync ;
              if (absSync < nodeMax && rm->nodeRefPos[absSync] < 0)
                { rm->nodeRefPos[absSync] = pos ;
                  rm->nodeOnRef[absSync] = true ;
                  ++rm->nMapped ;
                }
            }
        }
      seqhashIteratorDestroy (sit) ;
      break ; // only process first sequence in reference fasta
    }

  newFree (uBuf, rm->kh->plen, U64) ;
  seqhashDestroy (sh) ;
  seqIOclose (sio) ;

  fprintf (stderr, "reference %s %s: mapped %'lld of %'lld nodes\n",
           rm->sampleId, rm->seqName, rm->nMapped, kmerHashMax(rm->kh)) ;
  return rm ;
}

static void refMapDestroy (RefMap *rm)
{
  I64 nodeMax = kmerHashMax (rm->kh) + 1 ;
  newFree (rm->nodeRefPos, nodeMax, I64) ;
  newFree (rm->nodeOnRef, nodeMax, bool) ;
  free (rm->sampleId) ;
  if (rm->seqName) free (rm->seqName) ;
  newFree (rm, 1, RefMap) ;
}

/**************** region parsing ********************/

typedef struct {
  char *chrom ;
  I64   start ;
  I64   end ;
} Region ;

static Region parseRegion (char *regionStr)
{
  Region r = { 0, 0, 0 } ;
  char *s = strdup (regionStr) ;
  char *colon = strchr (s, ':') ;
  if (!colon) die ("invalid region format '%s', expected chr:start-end", regionStr) ;
  *colon = 0 ;
  r.chrom = strdup (s) ;
  char *dash = strchr (colon + 1, '-') ;
  if (!dash) die ("invalid region format '%s', expected chr:start-end", regionStr) ;
  *dash = 0 ;
  r.start = atoll (colon + 1) ;
  r.end = atoll (dash + 1) ;
  free (s) ;
  if (r.end <= r.start) die ("region end %lld <= start %lld", r.end, r.start) ;
  return r ;
}

/**************** node filtering by region ********************/

// returns a bool array indexed by node (1-based), true if node is in region
static bool *filterNodesByRegion (RefMap *rm, Region *region, I64 nodeMax, int syncLen)
{
  bool *inRegion = new0 (nodeMax, bool) ;
  I64 nIn = 0 ;
  for (I64 i = 1 ; i < nodeMax ; ++i)
    if (rm->nodeOnRef[i] && rm->nodeRefPos[i] >= 0)
      { I64 pos = rm->nodeRefPos[i] ;
        if (pos + syncLen > region->start && pos < region->end)
          { inRegion[i] = true ; ++nIn ; }
      }
  fprintf (stderr, "region %s:%lld-%lld: %'lld nodes in region\n",
           region->chrom, region->start, region->end, nIn) ;
  return inRegion ;
}

/**************** GFA output ********************/

// derive sample name and haplotype from source filename
// e.g. "HG01123_mat_hprc_r2_v1.0.1.chr19.fa.gz" -> sample "HG01123", hap 1 (mat)
// e.g. "CHM13_chm13v2.0.chr19.fa.gz" -> sample "CHM13", hap 0
static void parseSampleFromFilename (const char *filename, char *sampleBuf, int *hapIndex, char *seqBuf)
{
  // extract basename
  const char *base = strrchr (filename, '/') ;
  base = base ? base + 1 : filename ;

  // take up to first '_' as sample name
  const char *p = base ;
  int i = 0 ;
  while (*p && *p != '_' && i < 127) sampleBuf[i++] = *p++ ;
  sampleBuf[i] = 0 ;

  // check for _mat or _pat to determine haplotype
  *hapIndex = 0 ;
  if (strstr (base, "_mat")) *hapIndex = 1 ;
  else if (strstr (base, "_pat")) *hapIndex = 2 ;

  // try to extract chromosome from filename (e.g. ".chr19.")
  const char *chr = strstr (base, ".chr") ;
  if (!chr) chr = strstr (base, "_chr") ;
  if (chr)
    { chr = (*chr == '.') ? chr + 1 : chr + 1 ;
      i = 0 ;
      while (*chr && *chr != '.' && *chr != '_' && i < 127) seqBuf[i++] = *chr++ ;
      seqBuf[i] = 0 ;
    }
  else
    strcpy (seqBuf, "unknown") ;
}

static void writeGFA (FILE *f, SyngBWT *sb, KmerHash *kh, RefMap *rm,
                      bool *inRegion, I64 nodeMax, int syncLen,
                      OneFile *ofGbwt, SyncmerSet *sms)
{
  I64 i, nSeg = 0, nLink = 0 ;

  // header
  fprintf (f, "H\tVN:Z:1.1\n") ;

  // S-lines: one per existing node
  char *seqBuf = new (kh->len + 1, char) ;
  for (i = 1 ; i < nodeMax ; ++i)
    { if (!syngBWTnodeExists (sb, i)) continue ;
      if (inRegion && !inRegion[i]) continue ;
      kmerHashSeq (kh, i, seqBuf) ;
      // kmerHashSeq returns lowercase acgt - uppercase for GFA
      for (int j = 0 ; j < kh->len ; ++j)
        if (seqBuf[j] >= 'a' && seqBuf[j] <= 'z') seqBuf[j] -= 32 ;
      seqBuf[kh->len] = 0 ;

      fprintf (f, "S\ts%lld\t%s", i, seqBuf) ;
      // rGFA tags
      if (rm && rm->nodeOnRef[i] && rm->nodeRefPos[i] >= 0)
        fprintf (f, "\tSN:Z:%s\tSO:i:%lld\tSR:i:0", rm->seqName, rm->nodeRefPos[i]) ;
      else if (rm)
        fprintf (f, "\tSN:Z:%s\tSO:i:0\tSR:i:1", rm->seqName ? rm->seqName : "unknown") ;
      fprintf (f, "\n") ;
      ++nSeg ;
    }
  newFree (seqBuf, kh->len + 1, char) ;

  // L-lines: edges from the out-side of each node
  // For each node i in + orientation, emit out-edges
  for (i = 1 ; i < nodeMax ; ++i)
    { if (!syngBWTnodeExists (sb, i)) continue ;
      if (inRegion && !inRegion[i]) continue ;
      int nOut = syngBWTnodeOutDegree (sb, i) ;
      for (int j = 0 ; j < nOut ; ++j)
        { I32 adj ; U32 offset ; I64 count ;
          if (!syngBWTnodeOutEdge (sb, i, j, &adj, &offset, &count)) continue ;
          if (adj == 0) continue ; // path terminator
          I64 adjAbs = adj > 0 ? adj : -adj ;
          if (inRegion && !inRegion[adjAbs]) continue ;
          if (!syngBWTnodeExists (sb, adjAbs)) continue ;
          char fromOrient = '+' ;
          char toOrient = adj > 0 ? '+' : '-' ;
          fprintf (f, "L\ts%lld\t%c\ts%lld\t%c\t0M\tRC:i:%lld\n",
                   i, fromOrient, adjAbs, toOrient, count) ;
          ++nLink ;
        }
      // Also emit edges from the in-side (these represent edges arriving at -i, i.e. from -i's out)
      // In rGFA, we also need L-lines for edges leaving node i in - orientation
      int nIn = syngBWTnodeInDegree (sb, i) ;
      for (int j = 0 ; j < nIn ; ++j)
        { I32 adj ; U32 offset ; I64 count ;
          if (!syngBWTnodeInEdge (sb, i, j, &adj, &offset, &count)) continue ;
          if (adj == 0) continue ; // path start
          // in-edge with adj means: path comes from adj to +i
          // In the reverse orientation of i: -i has out-edges to -adj
          I64 adjAbs = adj > 0 ? adj : -adj ;
          if (inRegion && !inRegion[adjAbs]) continue ;
          if (!syngBWTnodeExists (sb, adjAbs)) continue ;
          // L-line: from i in - orientation to -adj (flipped)
          char toOrient = adj > 0 ? '-' : '+' ;
          fprintf (f, "L\ts%lld\t-\ts%lld\t%c\t0M\tRC:i:%lld\n",
                   i, adjAbs, toOrient, count) ;
          ++nLink ;
        }
    }

  // W-lines: haplotype paths
  I64 nPath = syngBWTpathCount (sb) ;
  I64 nWalk = 0 ;

  // collect source reference filenames for sample name derivation
  int nRef = ofGbwt ? oneReferenceCount (ofGbwt) : 0 ;
  char **refSample = 0, **refSeq = 0 ;
  int   *refHap = 0 ;
  if (nRef > 0)
    { refSample = new0 (nRef + 1, char*) ; // 1-based indexing
      refSeq = new0 (nRef + 1, char*) ;
      refHap = new0 (nRef + 1, int) ;
      for (int r = 0 ; r < nRef ; ++r)
        { refSample[r+1] = new (128, char) ;
          refSeq[r+1] = new (128, char) ;
          parseSampleFromFilename (ofGbwt->reference[r].filename,
                                  refSample[r+1], &refHap[r+1], refSeq[r+1]) ;
        }
    }

  if (nPath > 0 && ofGbwt)
    { I64 np ;
      if (oneStats (ofGbwt, 'P', &np, 0, 0) && np && oneGoto (ofGbwt, 'P', 1))
        { // position to first P line
          oneReadLine (ofGbwt) ;
          for (I64 pathIdx = 0 ; pathIdx < np ; ++pathIdx)
            { if (ofGbwt->lineType != 'P') break ;
              I64 len = oneInt(ofGbwt, 0) ;
              I64 source = oneInt(ofGbwt, 1) ;
              I64 inSource = oneInt(ofGbwt, 2) ;
              char *pathName = 0 ;
              I64 nSync = 0 ;
              I32 startSync = 0 ;
              U32 startPos = 0, startCount = 0 ;

              while (oneReadLine (ofGbwt) && ofGbwt->lineType != 'P' && ofGbwt->lineType != 'V')
                switch (ofGbwt->lineType)
                  { case 'I': pathName = strdup (oneString(ofGbwt)) ; break ;
                    case 'Z':
                      startSync = oneInt(ofGbwt, 0) ;
                      startPos = oneInt(ofGbwt, 1) ;
                      startCount = oneInt(ofGbwt, 2) ;
                      nSync = oneInt(ofGbwt, 3) ;
                      break ;
                  }

              if (nSync > 0)
                { char sampleBuf[128], seqBuf[128] ;
                  int hapIdx = 0 ;

                  if (pathName)
                    { char *p = strdup (pathName) ;
                      char *hash1 = strchr (p, '#') ;
                      if (hash1)
                        { *hash1 = 0 ;
                          strncpy (sampleBuf, p, 127) ; sampleBuf[127] = 0 ;
                          char *hash2 = strchr (hash1 + 1, '#') ;
                          if (hash2)
                            { *hash2 = 0 ;
                              hapIdx = atoi (hash1 + 1) ;
                              strncpy (seqBuf, hash2 + 1, 127) ; seqBuf[127] = 0 ;
                            }
                          else
                            { hapIdx = atoi (hash1 + 1) ;
                              strncpy (seqBuf, hash1 + 1, 127) ; seqBuf[127] = 0 ;
                            }
                        }
                      else
                        { strncpy (sampleBuf, p, 127) ; sampleBuf[127] = 0 ;
                          strcpy (seqBuf, "unknown") ;
                        }
                      free (p) ;
                    }
                  else if (refSample && source >= 1 && source <= nRef)
                    { strcpy (sampleBuf, refSample[source]) ;
                      strcpy (seqBuf, refSeq[source]) ;
                      hapIdx = refHap[source] ;
                    }
                  else
                    { sprintf (sampleBuf, "sample%lld", source) ;
                      sprintf (seqBuf, "seq%lld", inSource) ;
                    }

                  // build walk string
                  SyngBWTpath *sbp = syngBWTpathStartOld (sb, startSync, startCount) ;
                  I64 walkCap = nSync * 12, walkLen = 0 ;
                  char *walk = new (walkCap, char) ;

                  I64 absStart = startSync > 0 ? startSync : -startSync ;
                  bool anyInRegion = (!inRegion || (absStart < nodeMax && inRegion[absStart])) ;

                  char seg[32] ;
                  sprintf (seg, "%c" "s%lld", startSync > 0 ? '>' : '<', absStart) ;
                  walkLen = strlen (seg) ;
                  if (walkLen >= walkCap)
                    { walkCap *= 2 ; walk = newResize (walk, walkCap/2, walkCap, char) ; }
                  strcpy (walk, seg) ;

                  I32 nextSync ; U32 nextOff ;
                  for (I64 k = 1 ; k < nSync ; ++k)
                    { if (!syngBWTpathNext (sbp, &nextSync, &nextOff))
                        break ;
                      I64 absSync = nextSync > 0 ? nextSync : -nextSync ;
                      if (inRegion && absSync < nodeMax && inRegion[absSync]) anyInRegion = true ;
                      sprintf (seg, "%c" "s%lld", nextSync > 0 ? '>' : '<', absSync) ;
                      I64 segLen = strlen (seg) ;
                      while (walkLen + segLen + 1 >= walkCap)
                        { I64 oldCap = walkCap ; walkCap *= 2 ;
                          walk = newResize (walk, oldCap, walkCap, char) ;
                        }
                      strcpy (walk + walkLen, seg) ;
                      walkLen += segLen ;
                    }
                  syngBWTpathDestroy (sbp) ;

                  if (!inRegion || anyInRegion)
                    { fprintf (f, "W\t%s\t%d\t%s\t0\t%lld\t%s\n",
                               sampleBuf, hapIdx, seqBuf, len, walk) ;
                      ++nWalk ;
                    }
                  newFree (walk, walkCap, char) ;
                  if (pathName) free (pathName) ;
                }
              if (ofGbwt->lineType == 'V') break ;
            }
        }
    }

  // cleanup reference info
  if (refSample)
    { for (int r = 1 ; r <= nRef ; ++r)
        { newFree (refSample[r], 128, char) ; newFree (refSeq[r], 128, char) ; }
      newFree (refSample, nRef + 1, char*) ;
      newFree (refSeq, nRef + 1, char*) ;
      newFree (refHap, nRef + 1, int) ;
    }

  fprintf (stderr, "wrote GFA: %'lld segments, %'lld links, %'lld walks\n", nSeg, nLink, nWalk) ;
}

/**************** BED output ********************/

static void writeBED (FILE *f, SyngBWT *sb, KmerHash *kh, RefMap *rm,
                      bool *inRegion, I64 nodeMax, int syncLen,
                      OneFile *ofGbwt)
{
  if (!rm) { fprintf (stderr, "BED output requires -reference\n") ; return ; }

  I64 nPath = syngBWTpathCount (sb) ;
  I64 nBed = 0 ;

  // collect source reference filenames
  int nRef = ofGbwt ? oneReferenceCount (ofGbwt) : 0 ;
  char **refSample = 0, **refSeq = 0 ;
  int   *refHap = 0 ;
  if (nRef > 0)
    { refSample = new0 (nRef + 1, char*) ;
      refSeq = new0 (nRef + 1, char*) ;
      refHap = new0 (nRef + 1, int) ;
      for (int r = 0 ; r < nRef ; ++r)
        { refSample[r+1] = new (128, char) ;
          refSeq[r+1] = new (128, char) ;
          parseSampleFromFilename (ofGbwt->reference[r].filename,
                                  refSample[r+1], &refHap[r+1], refSeq[r+1]) ;
        }
    }

  if (nPath > 0 && ofGbwt)
    { I64 np ;
      if (oneStats (ofGbwt, 'P', &np, 0, 0) && np && oneGoto (ofGbwt, 'P', 1))
        { oneReadLine (ofGbwt) ;
          for (I64 pathIdx = 0 ; pathIdx < np ; ++pathIdx)
            { if (ofGbwt->lineType != 'P') break ;
              I64 len = oneInt(ofGbwt, 0) ;
              I64 source = oneInt(ofGbwt, 1) ;
              I64 inSource = oneInt(ofGbwt, 2) ;
              char *pathName = 0 ;
              I64 nSync = 0 ;
              I32 startSync = 0 ;
              U32 startPos = 0, startCount = 0 ;

              while (oneReadLine (ofGbwt) && ofGbwt->lineType != 'P' && ofGbwt->lineType != 'V')
                switch (ofGbwt->lineType)
                  { case 'I': pathName = strdup (oneString(ofGbwt)) ; break ;
                    case 'Z':
                      startSync = oneInt(ofGbwt, 0) ;
                      startPos = oneInt(ofGbwt, 1) ;
                      startCount = oneInt(ofGbwt, 2) ;
                      nSync = oneInt(ofGbwt, 3) ;
                      break ;
                  }

              if (nSync > 0)
                { char labelBuf[256] ;
                  if (pathName)
                    strncpy (labelBuf, pathName, 255) ;
                  else if (refSample && source >= 1 && source <= nRef)
                    sprintf (labelBuf, "%s_hap%d", refSample[source], refHap[source]) ;
                  else
                    sprintf (labelBuf, "sample%lld_seq%lld", source, inSource) ;
                  labelBuf[255] = 0 ;

                  SyngBWTpath *sbp = syngBWTpathStartOld (sb, startSync, startCount) ;
                  I64 absStart = startSync > 0 ? startSync : -startSync ;
                  I64 hapPos = 0 ;

                  if (absStart < nodeMax && rm->nodeOnRef[absStart] && rm->nodeRefPos[absStart] >= 0)
                    { I64 refPos = rm->nodeRefPos[absStart] ;
                      if (!inRegion || inRegion[absStart])
                        { fprintf (f, "%s\t%lld\t%lld\t%s\t%lld\t%c\n",
                                   rm->seqName, refPos, refPos + syncLen,
                                   labelBuf, hapPos, startSync > 0 ? '+' : '-') ;
                          ++nBed ;
                        }
                    }

                  I32 nextSync ; U32 nextOff ;
                  for (I64 k = 1 ; k < nSync ; ++k)
                    { if (!syngBWTpathNext (sbp, &nextSync, &nextOff)) break ;
                      hapPos += nextOff ;
                      I64 absSync = nextSync > 0 ? nextSync : -nextSync ;
                      if (absSync < nodeMax && rm->nodeOnRef[absSync] && rm->nodeRefPos[absSync] >= 0)
                        { I64 refPos = rm->nodeRefPos[absSync] ;
                          if (!inRegion || inRegion[absSync])
                            { fprintf (f, "%s\t%lld\t%lld\t%s\t%lld\t%c\n",
                                       rm->seqName, refPos, refPos + syncLen,
                                       labelBuf, hapPos, nextSync > 0 ? '+' : '-') ;
                              ++nBed ;
                            }
                        }
                    }
                  syngBWTpathDestroy (sbp) ;
                }
              if (pathName) free (pathName) ;
              if (ofGbwt->lineType == 'V') break ;
            }
        }
    }

  if (refSample)
    { for (int r = 1 ; r <= nRef ; ++r)
        { newFree (refSample[r], 128, char) ; newFree (refSeq[r], 128, char) ; }
      newFree (refSample, nRef + 1, char*) ;
      newFree (refSeq, nRef + 1, char*) ;
      newFree (refHap, nRef + 1, int) ;
    }
  fprintf (stderr, "wrote BED: %'lld entries\n", nBed) ;
}

/**************** main program ********************/

static char usage[] =
  "Usage: syngview [options]\n"
  "Required inputs:\n"
  "  -readK <.1khash file>    : syncmer hash table (provides node sequences)\n"
  "  -readGBWT <.1gbwt file>  : graph BWT with paths (provides haplotype paths)\n"
  "Optional inputs:\n"
  "  -reference <sample_id> <ref.fa> : reference sample for coordinate mapping\n"
  "Options:\n"
  "  -o <outfile prefix>      : [syngView] output prefix (.gfa, .bed)\n"
  "  -region <chr>:<start>-<end> : only output nodes in this genomic coordinate range\n"
  "                             (requires -reference)\n"
  "  -T <threads>             : [8] number of threads\n"
  "  -gfa                     : [default] write GFA output\n"
  "  -bed                     : write BED output for haplotype coordinates\n" ;

int main (int argc, char *argv[])
{
  char *khashFile = 0, *gbwtFile = 0 ;
  char *outPrefix = "syngView" ;
  char *refSampleId = 0, *refFasta = 0 ;
  char *regionStr = 0 ;
  int   nThread = 8 ;
  bool  doGFA = false, doBED = false ;

  timeUpdate (0) ;
  pathCount = 1 ; // avoid PATH_DEBUG==0 match in syngbwt3.c
  storeCommandLine (argc, argv) ;
  argc-- ; ++argv ;
  if (!argc) { fprintf (stderr, "%s", usage) ; exit (0) ; }

  while (argc > 0 && **argv == '-')
    { if (!strcmp (*argv, "-readK") && argc > 1)
        { khashFile = argv[1] ; argc -= 2 ; argv += 2 ; }
      else if (!strcmp (*argv, "-readGBWT") && argc > 1)
        { gbwtFile = argv[1] ; argc -= 2 ; argv += 2 ; }
      else if (!strcmp (*argv, "-reference") && argc > 2)
        { refSampleId = argv[1] ; refFasta = argv[2] ; argc -= 3 ; argv += 3 ; }
      else if (!strcmp (*argv, "-o") && argc > 1)
        { outPrefix = argv[1] ; argc -= 2 ; argv += 2 ; }
      else if (!strcmp (*argv, "-region") && argc > 1)
        { regionStr = argv[1] ; argc -= 2 ; argv += 2 ; }
      else if (!strcmp (*argv, "-T") && argc > 1)
        { nThread = atoi (argv[1]) ; argc -= 2 ; argv += 2 ; }
      else if (!strcmp (*argv, "-gfa"))
        { doGFA = true ; argc-- ; argv++ ; }
      else if (!strcmp (*argv, "-bed"))
        { doBED = true ; argc-- ; argv++ ; }
      else die ("unknown parameter %s\n%s", *argv, usage) ;
    }

  if (!khashFile) die ("must provide -readK <.1khash file>\n%s", usage) ;
  if (!gbwtFile) die ("must provide -readGBWT <.1gbwt file>\n%s", usage) ;
  if (!doGFA && !doBED) doGFA = true ; // default to GFA

  if (regionStr && !refFasta)
    die ("-region requires -reference\n%s", usage) ;

  // read the syncmer hash
  SyncmerSet *sms = syncmerSetRead (khashFile) ;
  KmerHash *kh = sms->kh ;
  int syncLen = kh->len ;
  fprintf (stderr, "read %'lld syncmers of length %d from %s\n",
           kmerHashMax(kh), syncLen, khashFile) ;
  timeUpdate (stderr) ;

  // read the GBWT
  OneSchema *schema = oneSchemaCreateFromText (syngSchemaText) ;
  OneFile *ofGbwt = oneFileOpenRead (gbwtFile, schema, "gbwt", nThread) ;
  if (!ofGbwt) die ("failed to open GBWT file %s", gbwtFile) ;

  // read syncmer params from GBWT and check
  if (oneGoto (ofGbwt, 'h', 1))
    { oneReadLine (ofGbwt) ;
      syncmerParamsCheck (ofGbwt, sms->params) ;
    }

  SyngBWT *sb = syngBWTread (ofGbwt) ;
  I64 nodeMax = syngBWTnodeMax (sb) ;
  fprintf (stderr, "GBWT: %'lld nodes, %'lld paths\n", nodeMax - 1, syngBWTpathCount (sb)) ;
  timeUpdate (stderr) ;

  // build reference map if requested
  RefMap *rm = 0 ;
  if (refFasta)
    { rm = refMapCreate (refSampleId, refFasta, sms, nThread) ;
      timeUpdate (stderr) ;
    }

  // parse region and filter nodes
  Region region = { 0, 0, 0 } ;
  bool *inRegion = 0 ;
  if (regionStr)
    { region = parseRegion (regionStr) ;
      inRegion = filterNodesByRegion (rm, &region, nodeMax, syncLen) ;
    }

  // write GFA
  if (doGFA)
    { char *fname = fnameTag (outPrefix, "gfa") ;
      FILE *f = fopen (fname, "w") ;
      if (!f) die ("failed to open GFA file %s for writing", fname) ;
      // re-open the gbwt file for path traversal (need fresh read position)
      OneFile *ofGbwt2 = oneFileOpenRead (gbwtFile, schema, 0, 1) ;
      writeGFA (f, sb, kh, rm, inRegion, nodeMax, syncLen, ofGbwt2, sms) ;
      fclose (f) ;
      if (ofGbwt2) oneFileClose (ofGbwt2) ;
      fprintf (stderr, "wrote GFA to %s\n", fname) ;
      timeUpdate (stderr) ;
    }

  // write BED
  if (doBED)
    { char *fname = fnameTag (outPrefix, "bed") ;
      FILE *f = fopen (fname, "w") ;
      if (!f) die ("failed to open BED file %s for writing", fname) ;
      OneFile *ofGbwt2 = oneFileOpenRead (gbwtFile, schema, 0, 1) ;
      writeBED (f, sb, kh, rm, inRegion, nodeMax, syncLen, ofGbwt2) ;
      fclose (f) ;
      if (ofGbwt2) oneFileClose (ofGbwt2) ;
      fprintf (stderr, "wrote BED to %s\n", fname) ;
      timeUpdate (stderr) ;
    }

  // cleanup
  if (inRegion) newFree (inRegion, nodeMax, bool) ;
  if (rm) refMapDestroy (rm) ;
  if (region.chrom) free (region.chrom) ;
  syngBWTdestroy (sb) ;
  oneFileClose (ofGbwt) ;
  oneSchemaDestroy (schema) ;
  syncmerSetDestroy (sms) ;

  fprintf (stderr, "total: ") ; timeTotal (stderr) ;
  return 0 ;
}

/*********************** end of file **********************/
