/*  File: syngview.c
 *  Author: cmarkello
 *  Copyright (C) 2026
 *-------------------------------------------------------------------
 * Description: Visualization exporter for spliced syncmer pangenome graphs.
 *   Reads a syncmer kmer hash, a spliced transcriptome annotation file,
 *   and optionally a GBWT with genomic sample paths, then outputs:
 *
 *   1. GFA (Graphical Fragment Assembly) format for tools like Bandage,
 *      odgi viz, or vg view. The GFA includes:
 *        - S (segment) lines for each syncmer node with its DNA sequence
 *        - L (link) lines for edges between adjacent syncmers
 *        - W (walk) lines for genomic sample paths
 *        - W (walk) lines for transcript annotation paths
 *
 *   2. BED format with transcript coordinates mapped back to reference
 *      positions, suitable for IGV or UCSC genome browser.
 *
 *   This allows researchers to see how transcript annotations from
 *   multiple samples/haplotypes thread through the pangenome graph
 *   alongside the underlying genomic variation.
 *
 * Usage:
 *   syngview -readK <.1khash> -readTx <.1sgtxom> [-readGBWT <.1gbwt>]
 *            [-region <node_from>:<node_to>] [-o <prefix>]
 *
 * HISTORY:
 * Created: Mar 30 2026
 *-------------------------------------------------------------------
 */

#include "seqio.h"
#include "seqhash.h"
#include "syng.h"

static OneSchema *schema ;

/**************** GFA output ******************/

// Write GFA header
static void gfaWriteHeader (FILE *f)
{
  fprintf(f, "H\tVN:Z:1.1\n") ;
}

// Write a segment line: S <id> <sequence>
static void gfaWriteSegment (FILE *f, I64 nodeId, KmerHash *kh)
{
  char *seq = kmerHashSeq(kh, nodeId, 0) ;
  fprintf(f, "S\t%lld\t%s\n", nodeId, seq) ;
}

// Write a link line: L <from> <from_orient> <to> <to_orient> <cigar>
static void gfaWriteLink (FILE *f, I32 from, I32 to, U32 offset, int syncmerLen)
{
  char fromOrient = (from >= 0) ? '+' : '-' ;
  char toOrient = (to >= 0) ? '+' : '-' ;
  I64 fromAbs = (from >= 0) ? from : -from ;
  I64 toAbs = (to >= 0) ? to : -to ;
  // overlap = syncmerLen - offset (if offset < syncmerLen)
  int overlap = (offset < syncmerLen) ? syncmerLen - offset : 0 ;
  fprintf(f, "L\t%lld\t%c\t%lld\t%c\t%dM\n", fromAbs, fromOrient, toAbs, toOrient, overlap) ;
}

// Write a walk line (GFA 1.1): W <sample> <hap> <seq> <start> <end> <walk>
// walk is a string of >node or <node segments
static void gfaWriteWalk (FILE *f, char *sampleName, int hapIdx,
                           char *seqName, I64 start, I64 end,
                           int nNodes, I32 *nodes)
{
  fprintf(f, "W\t%s\t%d\t%s\t%lld\t%lld\t", sampleName, hapIdx, seqName, start, end) ;
  int i ;
  for (i = 0 ; i < nNodes ; ++i)
    { char orient = (nodes[i] >= 0) ? '>' : '<' ;
      I64 id = (nodes[i] >= 0) ? nodes[i] : -nodes[i] ;
      fprintf(f, "%c%lld", orient, id) ;
    }
  fprintf(f, "\n") ;
}

/**************** BED output ******************/

// Write a BED12 line for a transcript with exon structure
static void bedWriteTranscript (FILE *f, char *chrom, I64 txStart, I64 txEnd,
                                 char *name, char strand, int nExons,
                                 I64 *exonStarts, I64 *exonEnds)
{
  fprintf(f, "%s\t%lld\t%lld\t%s\t0\t%c\t%lld\t%lld\t0,0,0\t%d\t",
          chrom, txStart, txEnd, name, strand, txStart, txEnd, nExons) ;
  // block sizes
  int i ;
  for (i = 0 ; i < nExons ; ++i)
    { fprintf(f, "%lld", exonEnds[i] - exonStarts[i]) ;
      if (i < nExons - 1) fprintf(f, ",") ;
    }
  fprintf(f, "\t") ;
  // block starts (relative to txStart)
  for (i = 0 ; i < nExons ; ++i)
    { fprintf(f, "%lld", exonStarts[i] - txStart) ;
      if (i < nExons - 1) fprintf(f, ",") ;
    }
  fprintf(f, "\n") ;
}

/**************** region filtering ******************/

typedef struct {
  I64 nodeFrom ;
  I64 nodeTo ;
  bool isActive ;
} Region ;

static bool nodeInRegion (I64 node, Region *region)
{
  if (!region->isActive) return true ;
  I64 absNode = (node >= 0) ? node : -node ;
  return (absNode >= region->nodeFrom && absNode <= region->nodeTo) ;
}

/**************** main ******************/

static char usage[] =
  "Usage: syngview [options]\n"
  "Export syncmer pangenome graph with transcript annotations for visualization.\n"
  "\n"
  "Required inputs:\n"
  "  -readK <.1khash file>    : syncmer hash table (provides node sequences)\n"
  "  -readTx <.1sgtxom file>  : spliced transcriptome annotation from syngannotate\n"
  "\n"
  "Optional inputs:\n"
  "  -readPaths <.1path file> : genomic sample paths to include in GFA\n"
  "\n"
  "Options:\n"
  "  -o <outfile prefix>      : [syngView] output prefix (.gfa, .bed)\n"
  "  -region <from>:<to>      : only output nodes in this ID range\n"
  "  -gfa                     : [default] write GFA output\n"
  "  -bed                     : write BED output for transcript coordinates\n" ;

int main (int argc, char *argv[])
{
  char       *outPrefix = "syngView" ;
  SyncmerSet *sms = 0 ;
  SyncmerParams params = syncmerParamsDefault() ;
  char       *txFile = 0 ;
  char       *pathFile = 0 ;
  bool        doGfa = true ;
  bool        doBed = false ;
  Region      region = { 0, 0, false } ;

  timeUpdate(0) ;
  schema = oneSchemaCreateFromText(syngSchemaText) ;

  storeCommandLine(argc, argv) ;
  --argc ; ++argv ;
  if (!argc) { fprintf(stderr, "%s", usage) ; exit(0) ; }

  while (argc > 0 && **argv == '-')
    if (!strcmp(*argv, "-readK") && argc > 1)
      { sms = syncmerSetRead(argv[1]) ;
        if (!sms) die("failed to read syncmer set from %s", argv[1]) ;
        fprintf(stdout, "read %llu syncmers (len %d)\n", kmerHashMax(sms->kh), sms->kh->len) ;
        argc -= 2 ; argv += 2 ;
      }
    else if (!strcmp(*argv, "-readTx") && argc > 1)
      { txFile = argv[1] ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp(*argv, "-readPaths") && argc > 1)
      { pathFile = argv[1] ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp(*argv, "-o") && argc > 1)
      { outPrefix = argv[1] ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp(*argv, "-region") && argc > 1)
      { if (sscanf(argv[1], "%lld:%lld", &region.nodeFrom, &region.nodeTo) == 2)
          region.isActive = true ;
        else
          die("-region requires format <from>:<to>, got %s", argv[1]) ;
        argc -= 2 ; argv += 2 ;
      }
    else if (!strcmp(*argv, "-gfa")) { doGfa = true ; --argc ; ++argv ; }
    else if (!strcmp(*argv, "-bed")) { doBed = true ; --argc ; ++argv ; }
    else die("unknown parameter %s\n%s", *argv, usage) ;

  if (!sms) die("must provide -readK <.1khash file>") ;
  if (!txFile) die("must provide -readTx <.1sgtxom file>") ;

  params = sms->params ;
  int syncmerLen = params.k + params.w ;
  KmerHash *kh = sms->kh ;

  // =========== GFA output ===========
  if (doGfa)
    { char *gfaName = fnameTag(outPrefix, "gfa") ;
      FILE *gfaFile = fopen(gfaName, "w") ;
      if (!gfaFile) die("failed to open %s for writing", gfaName) ;
      gfaWriteHeader(gfaFile) ;

      // collect all nodes referenced by transcripts and paths, then write segments
      // first pass: gather referenced node IDs from the transcriptome
      Hash usedNodes = hashCreate(65536) ;
      I64 nTx = 0, nJunctions = 0 ;

      OneFile *ofTx = oneFileOpenRead(txFile, schema, "sgtxom", 1) ;
      if (!ofTx) die("failed to open transcriptome file %s", txFile) ;

      // read transcript node lists to find all referenced nodes
      while (oneReadLine(ofTx))
        { if (ofTx->lineType == 'n')
            { I64 len = oneLen(ofTx) ;
              I64 *nodeList = oneIntList(ofTx) ;
              int i ;
              for (i = 0 ; i < len ; ++i)
                { I64 absNode = (nodeList[i] >= 0) ? nodeList[i] : -nodeList[i] ;
                  if (nodeInRegion(absNode, &region))
                    hashAdd(usedNodes, hashInt(absNode), 0) ;
                }
            }
          else if (ofTx->lineType == 'T') ++nTx ;
          else if (ofTx->lineType == 'J') ++nJunctions ;
        }
      oneFileClose(ofTx) ;

      // also gather nodes from genomic paths if provided
      if (pathFile)
        { OneFile *ofP = oneFileOpenRead(pathFile, schema, "path", 1) ;
          if (ofP)
            { while (oneReadLine(ofP))
                if (ofP->lineType == 'z')
                  { I64 len = oneLen(ofP) ;
                    I64 *nodeList = oneIntList(ofP) ;
                    int i ;
                    for (i = 0 ; i < len ; ++i)
                      { I64 absNode = (nodeList[i] >= 0) ? nodeList[i] : -nodeList[i] ;
                        if (nodeInRegion(absNode, &region))
                          hashAdd(usedNodes, hashInt(absNode), 0) ;
                      }
                  }
              oneFileClose(ofP) ;
            }
        }

      // if no region filter and no paths, write ALL segments (for completeness)
      if (!region.isActive && !pathFile)
        { I64 nodeId ;
          for (nodeId = 1 ; nodeId <= kmerHashMax(kh) ; ++nodeId)
            { gfaWriteSegment(gfaFile, nodeId, kh) ;
              hashAdd(usedNodes, hashInt(nodeId), 0) ;
            }
          fprintf(stdout, "wrote %llu segments to GFA\n", kmerHashMax(kh)) ;
        }
      else
        { // write only referenced/filtered segments
          I64 nSeg = 0 ;
          I64 nodeId ;
          for (nodeId = 1 ; nodeId <= kmerHashMax(kh) ; ++nodeId)
            if (hashFind(usedNodes, hashInt(nodeId), 0))
              { gfaWriteSegment(gfaFile, nodeId, kh) ; ++nSeg ; }
          fprintf(stdout, "wrote %lld segments to GFA\n", nSeg) ;
        }

      // second pass: write links and walks from the transcriptome
      ofTx = oneFileOpenRead(txFile, schema, "sgtxom", 1) ;
      if (!ofTx) die("failed to reopen transcriptome file %s", txFile) ;

      I64 nLinks = 0, nWalks = 0 ;
      // track written links to avoid duplicates
      Hash writtenLinks = hashCreate(65536) ;

      // state for current transcript
      char *curTxId = 0, *curGeneId = 0 ;
      I32   curSample = 0 ;
      I64  *curNodes = 0, *curOffsets = 0 ;
      I64   curNNodes = 0 ;

      while (oneReadLine(ofTx))
        switch (ofTx->lineType)
          {
          case 'T':
            curSample = oneInt(ofTx, 0) ;
            break ;
          case 'i':
            if (curTxId) free(curTxId) ;
            curTxId = strndup(oneString(ofTx), oneLen(ofTx)) ;
            break ;
          case 'g':
            if (curGeneId) free(curGeneId) ;
            curGeneId = strndup(oneString(ofTx), oneLen(ofTx)) ;
            break ;
          case 'n':
            curNNodes = oneLen(ofTx) ;
            curNodes = oneIntList(ofTx) ;
            // write links between consecutive nodes
            { int i ;
              for (i = 0 ; i < curNNodes - 1 ; ++i)
                { I32 from = curNodes[i], to = curNodes[i+1] ;
                  I64 fromAbs = (from >= 0) ? from : -from ;
                  I64 toAbs = (to >= 0) ? to : -to ;
                  if (!nodeInRegion(fromAbs, &region) || !nodeInRegion(toAbs, &region))
                    continue ;
                  // dedup links using hash of (from, to) pair
                  if (hashAdd(writtenLinks, hashInt2(from, to), 0))
                    { gfaWriteLink(gfaFile, from, to, 0, syncmerLen) ;
                      ++nLinks ;
                    }
                }
            }
            break ;
          case 'f':
            curOffsets = oneIntList(ofTx) ;
            // now we have a complete transcript - write walk
            if (curNNodes > 0 && curTxId)
              { // check if any node is in region
                bool inRegion = false ;
                int i ;
                for (i = 0 ; i < curNNodes ; ++i)
                  if (nodeInRegion(curNodes[i], &region)) { inRegion = true ; break ; }
                if (inRegion)
                  { // compute bp length
                    I64 bpLen = 0 ;
                    for (i = 0 ; i < curNNodes ; ++i) bpLen += curOffsets[i] ;
                    bpLen += syncmerLen ;
                    // create walk label: sample_txId
                    char walkName[512] ;
                    snprintf(walkName, sizeof(walkName), "tx_%s",
                             curTxId ? curTxId : "unknown") ;
                    I32 *nodes32 = new(curNNodes, I32) ;
                    for (i = 0 ; i < curNNodes ; ++i) nodes32[i] = curNodes[i] ;
                    char sampleBuf[64] ;
                    snprintf(sampleBuf, sizeof(sampleBuf), "sample%d", curSample) ;
                    gfaWriteWalk(gfaFile, sampleBuf, curSample,
                                 curGeneId ? curGeneId : "unknown",
                                 0, bpLen, curNNodes, nodes32) ;
                    newFree(nodes32, curNNodes, I32) ;
                    ++nWalks ;
                  }
              }
            curNNodes = 0 ;
            break ;
          case 'J':
            { // splice junction - write as a link
              I32 donor = oneInt(ofTx, 0) ;
              I32 acceptor = oneInt(ofTx, 1) ;
              I64 donorAbs = (donor >= 0) ? donor : -donor ;
              I64 accAbs = (acceptor >= 0) ? acceptor : -acceptor ;
              if (nodeInRegion(donorAbs, &region) && nodeInRegion(accAbs, &region))
                { if (hashAdd(writtenLinks, hashInt2(donor, acceptor), 0))
                    { // splice junctions get 0* overlap to mark them as spliced
                      fprintf(gfaFile, "L\t%lld\t%c\t%lld\t%c\t0M\tSJ:Z:splice\n",
                              donorAbs, (donor >= 0) ? '+' : '-',
                              accAbs, (acceptor >= 0) ? '+' : '-') ;
                      ++nLinks ;
                    }
                }
            }
            break ;
          }

      oneFileClose(ofTx) ;

      // add genomic sample paths as walks
      if (pathFile)
        { OneFile *ofP = oneFileOpenRead(pathFile, schema, "path", 1) ;
          if (ofP)
            { int pathIdx = 0 ;
              I64 pathLen = 0 ;
              I32 pathSource = 0 ;
              while (oneReadLine(ofP))
                { if (ofP->lineType == 'P')
                    { pathLen = oneInt(ofP, 0) ;
                      pathSource = oneInt(ofP, 1) ;
                      ++pathIdx ;
                    }
                  else if (ofP->lineType == 'z')
                    { I64 len = oneLen(ofP) ;
                      I64 *nodeList = oneIntList(ofP) ;
                      // write links
                      int i ;
                      for (i = 0 ; i < len - 1 ; ++i)
                        { I32 from = nodeList[i], to = nodeList[i+1] ;
                          I64 fromAbs = (from >= 0) ? from : -from ;
                          I64 toAbs = (to >= 0) ? to : -to ;
                          if (!nodeInRegion(fromAbs, &region) || !nodeInRegion(toAbs, &region))
                            continue ;
                          if (hashAdd(writtenLinks, hashInt2(from, to), 0))
                            { gfaWriteLink(gfaFile, from, to, 0, syncmerLen) ;
                              ++nLinks ;
                            }
                        }
                      // write walk
                      I32 *nodes32 = new(len, I32) ;
                      for (i = 0 ; i < len ; ++i) nodes32[i] = nodeList[i] ;
                      char sampleBuf[64] ;
                      snprintf(sampleBuf, sizeof(sampleBuf), "genome%d", pathSource) ;
                      gfaWriteWalk(gfaFile, sampleBuf, pathSource,
                                   "chr", 0, pathLen, len, nodes32) ;
                      newFree(nodes32, len, I32) ;
                      ++nWalks ;
                    }
                }
              oneFileClose(ofP) ;
            }
        }

      hashDestroy(writtenLinks) ;
      hashDestroy(usedNodes) ;
      fclose(gfaFile) ;
      fprintf(stdout, "wrote GFA: %lld links, %lld walks (transcripts + paths) to %s\n",
              nLinks, nWalks, gfaName) ;
      if (curTxId) free(curTxId) ;
      if (curGeneId) free(curGeneId) ;
      timeUpdate(stdout) ;
    }

  // =========== BED output ===========
  if (doBed)
    { char *bedName = fnameTag(outPrefix, "bed") ;
      FILE *bedFile = fopen(bedName, "w") ;
      if (!bedFile) die("failed to open %s for writing", bedName) ;

      // write BED header
      fprintf(bedFile, "#chrom\tchromStart\tchromEnd\tname\tscore\tstrand\n") ;

      // read transcript annotations and compute approximate genomic positions
      OneFile *ofTx = oneFileOpenRead(txFile, schema, "sgtxom", 1) ;
      if (!ofTx) die("failed to open transcriptome file %s", txFile) ;

      char *curTxId = 0, *curGeneId = 0 ;
      I32   curSample = 0 ;
      I64  *curNodes = 0, *curOffsets = 0 ;
      I64   curNNodes = 0 ;
      I64  *curExonBounds = 0 ;
      I64   curNExons = 0 ;
      I64   nBed = 0 ;

      while (oneReadLine(ofTx))
        switch (ofTx->lineType)
          {
          case 'T':
            curSample = oneInt(ofTx, 0) ;
            break ;
          case 'i':
            if (curTxId) free(curTxId) ;
            curTxId = strndup(oneString(ofTx), oneLen(ofTx)) ;
            break ;
          case 'g':
            if (curGeneId) free(curGeneId) ;
            curGeneId = strndup(oneString(ofTx), oneLen(ofTx)) ;
            break ;
          case 'n':
            curNNodes = oneLen(ofTx) ;
            curNodes = oneIntList(ofTx) ;
            break ;
          case 'f':
            curOffsets = oneIntList(ofTx) ;
            break ;
          case 'x':
            curNExons = oneLen(ofTx) ;
            curExonBounds = oneIntList(ofTx) ;
            // now we have a complete transcript with exon info - write BED
            if (curNNodes > 0 && curTxId)
              { // compute cumulative positions from offsets
                I64 pos = 0 ;
                I64 txStart = 0, txEnd = 0 ;
                int i ;
                for (i = 0 ; i < curNNodes ; ++i)
                  { pos += curOffsets[i] ; }
                txEnd = pos + syncmerLen ;
                // write BED with exon blocks
                if (curNExons > 0 && curExonBounds)
                  { I64 *exStarts = new(curNExons, I64) ;
                    I64 *exEnds = new(curNExons, I64) ;
                    for (i = 0 ; i < curNExons ; ++i)
                      { I64 exStartIdx = curExonBounds[i] ;
                        I64 exEndIdx = (i+1 < curNExons) ? curExonBounds[i+1] - 1
                                                          : curNNodes - 1 ;
                        // compute positions
                        I64 sPos = 0 ;
                        int j ;
                        for (j = 0 ; j <= exStartIdx && j < curNNodes ; ++j)
                          sPos += curOffsets[j] ;
                        exStarts[i] = sPos ;
                        I64 ePos = 0 ;
                        for (j = 0 ; j <= exEndIdx && j < curNNodes ; ++j)
                          ePos += curOffsets[j] ;
                        exEnds[i] = ePos + syncmerLen ;
                      }
                    char nameBuf[512] ;
                    snprintf(nameBuf, sizeof(nameBuf), "%s|%s|s%d",
                             curTxId, curGeneId ? curGeneId : ".", curSample) ;
                    bedWriteTranscript(bedFile, curGeneId ? curGeneId : ".",
                                       exStarts[0], exEnds[curNExons-1],
                                       nameBuf, '+', curNExons,
                                       exStarts, exEnds) ;
                    newFree(exStarts, curNExons, I64) ;
                    newFree(exEnds, curNExons, I64) ;
                  }
                else // single-exon transcript
                  fprintf(bedFile, "%s\t%lld\t%lld\t%s|%s|s%d\t0\t+\n",
                          curGeneId ? curGeneId : ".",
                          txStart, txEnd, curTxId,
                          curGeneId ? curGeneId : ".", curSample) ;
                ++nBed ;
              }
            curNNodes = 0 ;
            curNExons = 0 ;
            break ;
          }

      oneFileClose(ofTx) ;
      fclose(bedFile) ;
      if (curTxId) free(curTxId) ;
      if (curGeneId) free(curGeneId) ;
      fprintf(stdout, "wrote %lld transcript entries to %s\n", nBed, bedName) ;
      timeUpdate(stdout) ;
    }

  if (sms) syncmerSetDestroy(sms) ;
  fprintf(stdout, "total: ") ; timeTotal(stdout) ;
  return 0 ;
}

/*********************** end of file **********************/
