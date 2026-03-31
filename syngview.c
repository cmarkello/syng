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
 *   syngview -readK <.1khash> -readTx <.1sgtxom> [-readPaths <.1path>]
 *            [-reference <sample_id> <ref.fa>] [-region chr:start-end]
 *            [-gfa] [-bed] [-o <prefix>]
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

/************ mapping genomic coordinates to syncmer nodes ************/

// SyncmerCoord: maps a genomic position to a syncmer node
typedef struct {
  I64 genomePos ;    // position in the reference genome sequence
  I32 syncNode ;     // syncmer node ID (negative if reverse complement)
} SyncmerCoord ;

// Build a coordinate map for a reference sequence: scan it for syncmers
// and record each syncmer's genomic position and node ID.
static Array buildSyncmerCoordMap (char *seq, I64 seqLen, Seqhash *sh, KmerHash *kh)
{
  Array coordMap = arrayCreate(seqLen / 50, SyncmerCoord) ;
  U64 *uBuf = new(kh->plen, U64) ;

  // convert sequence to index form
  char *seqIdx = new(seqLen, char) ;
  I64 i ;
  for (i = 0 ; i < seqLen ; ++i)
    seqIdx[i] = dna2index4Conv[(unsigned char)seq[i]] ;

  SeqhashIterator *sit = syncmerIterator(sh, seqIdx, seqLen) ;
  int pos ;
  while (syncmerNext(sit, 0, &pos, 0))
    { I64 sync = 0 ;
      kmerHashFindThreadSafe(kh, seqIdx + pos, &sync, uBuf) ;
      if (sync && (sync > 2 || sync < -2)) // skip poly-X
        { SyncmerCoord *sc = arrayp(coordMap, arrayMax(coordMap), SyncmerCoord) ;
          sc->genomePos = pos ;
          sc->syncNode = sync ;
        }
    }
  seqhashIteratorDestroy(sit) ;
  newFree(uBuf, kh->plen, U64) ;
  newFree(seqIdx, seqLen, char) ;

  return coordMap ;
}

/**************** region filtering ******************/

typedef struct {
  char  chrName[256] ;  // chromosome name from -region
  I64   startPos ;      // start coordinate (0-based)
  I64   endPos ;        // end coordinate
  bool  isActive ;
  Hash  nodeHash ;      // hash of node IDs in the region (built from coord map)
} Region ;

static bool nodeInRegion (I64 node, Region *region)
{
  if (!region->isActive) return true ;
  I64 absNode = (node >= 0) ? node : -node ;
  return hashFind(region->nodeHash, hashInt(absNode), 0) ;
}

// Populate region->nodeHash with all syncmer nodes whose genomic position
// falls within [region->startPos, region->endPos] on the matching chromosome.
static void regionBuildNodeHash (Region *region, Array coordMap, int syncmerLen)
{
  region->nodeHash = hashCreate(4096) ;
  I64 i ;
  for (i = 0 ; i < arrayMax(coordMap) ; ++i)
    { SyncmerCoord *sc = arrp(coordMap, i, SyncmerCoord) ;
      // a syncmer at genomePos covers [genomePos, genomePos + syncmerLen)
      if (sc->genomePos + syncmerLen > region->startPos && sc->genomePos < region->endPos)
        { I64 absNode = (sc->syncNode >= 0) ? sc->syncNode : -sc->syncNode ;
          hashAdd(region->nodeHash, hashInt(absNode), 0) ;
        }
    }
}

// Node-to-position map: hash maps node ID -> index into posArray
typedef struct {
  Hash  hash ;       // node ID -> index
  Array posArray ;   // parallel array of I64 genomic positions
} NodePosMap ;

static NodePosMap *buildNodePosMap (Array coordMap)
{
  NodePosMap *npm = new0(1, NodePosMap) ;
  npm->hash = hashCreate(4096) ;
  npm->posArray = arrayCreate(4096, I64) ;
  I64 i ;
  for (i = 0 ; i < arrayMax(coordMap) ; ++i)
    { SyncmerCoord *sc = arrp(coordMap, i, SyncmerCoord) ;
      I64 absNode = (sc->syncNode >= 0) ? sc->syncNode : -sc->syncNode ;
      int idx ;
      if (hashAdd(npm->hash, hashInt(absNode), &idx))
        array(npm->posArray, idx, I64) = sc->genomePos ;  // store first occurrence
    }
  return npm ;
}

static bool nodePosMapFind (NodePosMap *npm, I64 absNode, I64 *posOut)
{
  int idx ;
  if (hashFind(npm->hash, hashInt(absNode), &idx))
    { *posOut = arr(npm->posArray, idx, I64) ;
      return true ;
    }
  return false ;
}

static void nodePosMapDestroy (NodePosMap *npm)
{
  if (!npm) return ;
  hashDestroy(npm->hash) ;
  arrayDestroy(npm->posArray) ;
  newFree(npm, 1, NodePosMap) ;
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
  "  -reference <sample_id> <ref.fa> : reference sample for coordinate mapping\n"
  "\n"
  "Options:\n"
  "  -o <outfile prefix>      : [syngView] output prefix (.gfa, .bed)\n"
  "  -region <chr>:<start>-<end> : only output nodes in this genomic coordinate range\n"
  "                             (requires -reference)\n"
  "  -gfa                     : [default] write GFA output\n"
  "  -bed                     : write BED output for transcript coordinates\n" ;

int main (int argc, char *argv[])
{
  char       *outPrefix = "syngView" ;
  SyncmerSet *sms = 0 ;
  SyncmerParams params = syncmerParamsDefault() ;
  char       *txFile = 0 ;
  char       *pathFile = 0 ;
  char       *refSampleId = 0 ;
  char       *refFasta = 0 ;
  bool        doGfa = true ;
  bool        doBed = false ;
  Region      region ;
  memset(&region, 0, sizeof(Region)) ;

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
    else if (!strcmp(*argv, "-reference") && argc > 2)
      { refSampleId = argv[1] ; refFasta = argv[2] ;
        argc -= 3 ; argv += 3 ;
      }
    else if (!strcmp(*argv, "-region") && argc > 1)
      { // parse chr:start-end format
        char *arg = argv[1] ;
        char *colon = strchr(arg, ':') ;
        if (!colon) die("-region requires format <chr>:<start>-<end>, got %s", arg) ;
        int chrLen = colon - arg ;
        if (chrLen >= (int)sizeof(region.chrName))
          die("-region chromosome name too long") ;
        strncpy(region.chrName, arg, chrLen) ;
        region.chrName[chrLen] = '\0' ;
        char *rest = colon + 1 ;
        char *dash = strchr(rest, '-') ;
        if (!dash)
          die("-region requires format <chr>:<start>-<end>, got %s", arg) ;
        region.startPos = atoll(rest) ;
        region.endPos = atoll(dash + 1) ;
        region.isActive = true ;
        argc -= 2 ; argv += 2 ;
      }
    else if (!strcmp(*argv, "-gfa")) { doGfa = true ; --argc ; ++argv ; }
    else if (!strcmp(*argv, "-bed")) { doBed = true ; --argc ; ++argv ; }
    else die("unknown parameter %s\n%s", *argv, usage) ;

  if (!sms) die("must provide -readK <.1khash file>") ;
  if (!txFile) die("must provide -readTx <.1sgtxom file>") ;
  if (region.isActive && !refFasta)
    die("-region requires -reference <sample_id> <ref.fa>") ;

  params = sms->params ;
  int syncmerLen = params.k + params.w ;
  KmerHash *kh = sms->kh ;

  // Build coordinate map from reference if provided
  DICT *refSeqDict = 0 ;          // chr name -> index
  Array refCoordMaps = 0 ;        // Array of (Array of SyncmerCoord), one per chr
  NodePosMap *nodePosMap = 0 ;    // node ID -> genomic position (for BED output)

  if (refFasta)
    { Seqhash *sh = seqhashCreate(params.k, params.w+1, params.seed) ;
      fprintf(stdout, "scanning reference %s for syncmer coordinates...\n", refFasta) ;
      SeqIO *sio = seqIOopenRead(refFasta, dna2indexConv, 0) ;
      if (!sio) die("failed to open reference FASTA %s", refFasta) ;

      refSeqDict = dictCreate(64) ;
      refCoordMaps = arrayCreate(64, Array) ;

      while (seqIOread(sio))
        { char *seqName = strdup(sqioId(sio)) ;
          U64 seqIdx ;
          dictAdd(refSeqDict, seqName, &seqIdx) ;

          char *seqIdx4 = new(sio->seqLen, char) ;
          I64 j ;
          for (j = 0 ; j < sio->seqLen ; ++j)
            seqIdx4[j] = sqioSeq(sio)[j] ;

          Array cmap = buildSyncmerCoordMap(seqIdx4, sio->seqLen, sh, kh) ;
          array(refCoordMaps, seqIdx, Array) = cmap ;
          newFree(seqIdx4, sio->seqLen, char) ;
          free(seqName) ;
        }
      seqIOclose(sio) ;
      seqhashDestroy(sh) ;
      fprintf(stdout, "built coordinate maps for %lld reference sequences\n",
              (long long)dictMax(refSeqDict)) ;

      // if region filter is active, build the node hash
      if (region.isActive)
        { U64 chrIdx ;
          if (!dictFind(refSeqDict, region.chrName, &chrIdx))
            die("chromosome %s not found in reference %s", region.chrName, refFasta) ;
          Array cmap = arr(refCoordMaps, chrIdx, Array) ;
          regionBuildNodeHash(&region, cmap, syncmerLen) ;
          fprintf(stdout, "region %s:%lld-%lld contains %d nodes\n",
                  region.chrName, region.startPos, region.endPos,
                  hashCount(region.nodeHash)) ;
        }

      // build node-to-position map (for BED and walk coordinate output)
      // Use the region chromosome if active, otherwise first chr
      if (region.isActive)
        { U64 chrIdx ;
          dictFind(refSeqDict, region.chrName, &chrIdx) ;
          nodePosMap = buildNodePosMap(arr(refCoordMaps, chrIdx, Array)) ;
        }
      else if (arrayMax(refCoordMaps) > 0)
        nodePosMap = buildNodePosMap(arr(refCoordMaps, 0, Array)) ;

      timeUpdate(stdout) ;
    }

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
              { // find the sub-range of nodes that are in the region
                int clipFirst = 0, clipLast = curNNodes - 1 ;
                if (region.isActive)
                  { clipFirst = -1 ; clipLast = -1 ;
                    int i ;
                    for (i = 0 ; i < curNNodes ; ++i)
                      if (nodeInRegion(curNodes[i], &region))
                        { if (clipFirst < 0) clipFirst = i ;
                          clipLast = i ;
                        }
                  }
                if (clipFirst >= 0)
                  { int clipN = clipLast - clipFirst + 1 ;
                    I64 walkStart = 0, walkEnd = 0 ;
                    char *walkChr = curGeneId ? curGeneId : "unknown" ;
                    if (nodePosMap)
                      { I64 fAbs = (curNodes[clipFirst] >= 0) ?
                                    curNodes[clipFirst] : -curNodes[clipFirst] ;
                        I64 lAbs = (curNodes[clipLast] >= 0) ?
                                    curNodes[clipLast] : -curNodes[clipLast] ;
                        I64 pos ;
                        if (nodePosMapFind(nodePosMap, fAbs, &pos))
                          walkStart = pos ;
                        if (nodePosMapFind(nodePosMap, lAbs, &pos))
                          walkEnd = pos + syncmerLen ;
                        if (walkStart > walkEnd)
                          { I64 tmp = walkStart ; walkStart = walkEnd ; walkEnd = tmp ; }
                        if (region.isActive) walkChr = region.chrName ;
                        // skip walk if coordinates fall outside the region window
                        if (region.isActive &&
                            (walkEnd <= region.startPos || walkStart >= region.endPos))
                          goto skipTxWalk ;
                      }
                    else
                      { int i ;
                        I64 bpLen = 0 ;
                        for (i = clipFirst ; i <= clipLast ; ++i) bpLen += curOffsets[i] ;
                        walkEnd = bpLen + syncmerLen ;
                      }
                    { I32 *nodes32 = new(clipN, I32) ;
                      int i ;
                      for (i = 0 ; i < clipN ; ++i)
                        nodes32[i] = curNodes[clipFirst + i] ;
                      char sampleBuf[64] ;
                      snprintf(sampleBuf, sizeof(sampleBuf), "sample%d", curSample) ;
                      gfaWriteWalk(gfaFile, sampleBuf, curSample,
                                   walkChr, walkStart, walkEnd,
                                   clipN, nodes32) ;
                      newFree(nodes32, clipN, I32) ;
                      ++nWalks ;
                    }
                  skipTxWalk: ;
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
                      // write walk - clip to region if active
                      int cFirst = 0, cLast = len - 1 ;
                      if (region.isActive)
                        { cFirst = -1 ; cLast = -1 ;
                          for (i = 0 ; i < len ; ++i)
                            if (nodeInRegion(nodeList[i], &region))
                              { if (cFirst < 0) cFirst = i ;
                                cLast = i ;
                              }
                        }
                      if (cFirst >= 0)
                        { int cN = cLast - cFirst + 1 ;
                          I64 wStart = 0, wEnd = pathLen ;
                          char *wChr = "chr" ;
                          bool skipPath = false ;
                          if (nodePosMap && cN > 0)
                            { I64 fAbs = (nodeList[cFirst] >= 0) ?
                                          nodeList[cFirst] : -nodeList[cFirst] ;
                              I64 lAbs = (nodeList[cLast] >= 0) ?
                                          nodeList[cLast] : -nodeList[cLast] ;
                              I64 pos ;
                              if (nodePosMapFind(nodePosMap, fAbs, &pos))
                                wStart = pos ;
                              if (nodePosMapFind(nodePosMap, lAbs, &pos))
                                wEnd = pos + syncmerLen ;
                              if (wStart > wEnd)
                                { I64 tmp = wStart ; wStart = wEnd ; wEnd = tmp ; }
                              if (region.isActive) wChr = region.chrName ;
                              if (region.isActive &&
                                  (wEnd <= region.startPos || wStart >= region.endPos))
                                skipPath = true ;
                            }
                          if (!skipPath)
                            { I32 *nodes32 = new(cN, I32) ;
                              for (i = 0 ; i < cN ; ++i)
                                nodes32[i] = nodeList[cFirst + i] ;
                              char sampleBuf[64] ;
                              snprintf(sampleBuf, sizeof(sampleBuf), "genome%d", pathSource) ;
                              gfaWriteWalk(gfaFile, sampleBuf, pathSource,
                                           wChr, wStart, wEnd, cN, nodes32) ;
                              newFree(nodes32, cN, I32) ;
                              ++nWalks ;
                            }
                        }
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
              { // region filter: skip transcripts with no nodes in region
                if (region.isActive)
                  { bool inRegion = false ;
                    int i ;
                    for (i = 0 ; i < curNNodes ; ++i)
                      if (nodeInRegion(curNodes[i], &region)) { inRegion = true ; break ; }
                    if (!inRegion) { curNNodes = 0 ; curNExons = 0 ; break ; }
                  }
                // determine chromosome name
                char *bedChr = curGeneId ? curGeneId : "." ;
                if (region.isActive) bedChr = region.chrName ;
                else if (refSampleId) bedChr = "ref" ;

                if (nodePosMap)
                  { // use real genomic coordinates from reference
                    I64 firstPos = 0, lastPos = 0 ;
                    bool hasFirst = false, hasLast = false ;
                    I64 firstAbs = (curNodes[0] >= 0) ? curNodes[0] : -curNodes[0] ;
                    I64 lastAbs = (curNodes[curNNodes-1] >= 0) ?
                                   curNodes[curNNodes-1] : -curNodes[curNNodes-1] ;
                    I64 pos ;
                    if (nodePosMapFind(nodePosMap, firstAbs, &pos))
                      { firstPos = pos ; hasFirst = true ; }
                    if (nodePosMapFind(nodePosMap, lastAbs, &pos))
                      { lastPos = pos ; hasLast = true ; }

                    // detect strand: if first node maps after last, it's minus strand
                    bool isReversed = hasFirst && hasLast && firstPos > lastPos ;
                    char strand = isReversed ? '-' : '+' ;
                    I64 txStart = isReversed ? lastPos : firstPos ;
                    I64 txEnd = (isReversed ? firstPos : lastPos) + syncmerLen ;

                    if (curNExons > 0 && curExonBounds)
                      { I64 *exStarts = new(curNExons, I64) ;
                        I64 *exEnds = new(curNExons, I64) ;
                        int i ;
                        for (i = 0 ; i < curNExons ; ++i)
                          { // for reversed transcripts, iterate exons in reverse
                            int ei = isReversed ? (curNExons - 1 - i) : i ;
                            I64 exStartIdx = curExonBounds[ei] ;
                            I64 exEndIdx = (ei+1 < curNExons) ? curExonBounds[ei+1] - 1
                                                               : curNNodes - 1 ;
                            // map exon boundary nodes to genomic positions
                            I64 absS = (curNodes[exStartIdx] >= 0) ?
                                        curNodes[exStartIdx] : -curNodes[exStartIdx] ;
                            I64 absE = (curNodes[exEndIdx] >= 0) ?
                                        curNodes[exEndIdx] : -curNodes[exEndIdx] ;
                            I64 posS = txStart, posE = txEnd ;
                            if (nodePosMapFind(nodePosMap, absS, &pos)) posS = pos ;
                            if (nodePosMapFind(nodePosMap, absE, &pos)) posE = pos + syncmerLen ;
                            // for reversed exons, the start/end may be swapped
                            if (posS > posE)
                              { I64 tmp = posS ; posS = posE ; posE = tmp ; }
                            exStarts[i] = posS ;
                            exEnds[i] = posE ;
                          }
                        // ensure monotonic ordering and valid block sizes
                        for (i = 1 ; i < curNExons ; ++i)
                          { if (exStarts[i] < exStarts[i-1])
                              exStarts[i] = exStarts[i-1] ;
                            if (exEnds[i] < exStarts[i])
                              exEnds[i] = exStarts[i] + syncmerLen ;
                          }
                        // ensure txStart/txEnd encompass all exons
                        if (exStarts[0] < txStart) txStart = exStarts[0] ;
                        if (exEnds[curNExons-1] > txEnd) txEnd = exEnds[curNExons-1] ;

                        char nameBuf[512] ;
                        snprintf(nameBuf, sizeof(nameBuf), "%s|%s|s%d",
                                 curTxId, curGeneId ? curGeneId : ".", curSample) ;
                        bedWriteTranscript(bedFile, bedChr,
                                           txStart, txEnd,
                                           nameBuf, strand, curNExons,
                                           exStarts, exEnds) ;
                        newFree(exStarts, curNExons, I64) ;
                        newFree(exEnds, curNExons, I64) ;
                      }
                    else
                      { if (txStart > txEnd)
                          { I64 tmp = txStart ; txStart = txEnd ; txEnd = tmp ; }
                        fprintf(bedFile, "%s\t%lld\t%lld\t%s|%s|s%d\t0\t%c\n",
                                bedChr, txStart, txEnd, curTxId,
                                curGeneId ? curGeneId : ".", curSample, strand) ;
                      }
                  }
                else
                  { // no reference: use offset-based positions (original behavior)
                    I64 pos = 0 ;
                    I64 txStart = 0, txEnd = 0 ;
                    int i ;
                    for (i = 0 ; i < curNNodes ; ++i)
                      { pos += curOffsets[i] ; }
                    txEnd = pos + syncmerLen ;
                    if (curNExons > 0 && curExonBounds)
                      { I64 *exStarts = new(curNExons, I64) ;
                        I64 *exEnds = new(curNExons, I64) ;
                        for (i = 0 ; i < curNExons ; ++i)
                          { I64 exStartIdx = curExonBounds[i] ;
                            I64 exEndIdx = (i+1 < curNExons) ? curExonBounds[i+1] - 1
                                                              : curNNodes - 1 ;
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
                        bedWriteTranscript(bedFile, bedChr,
                                           exStarts[0], exEnds[curNExons-1],
                                           nameBuf, '+', curNExons,
                                           exStarts, exEnds) ;
                        newFree(exStarts, curNExons, I64) ;
                        newFree(exEnds, curNExons, I64) ;
                      }
                    else
                      fprintf(bedFile, "%s\t%lld\t%lld\t%s|%s|s%d\t0\t+\n",
                              bedChr, txStart, txEnd, curTxId,
                              curGeneId ? curGeneId : ".", curSample) ;
                  }
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

  // cleanup
  if (region.isActive && region.nodeHash) hashDestroy(region.nodeHash) ;
  if (nodePosMap) nodePosMapDestroy(nodePosMap) ;
  if (refCoordMaps)
    { int c ;
      for (c = 0 ; c < arrayMax(refCoordMaps) ; ++c)
        { Array cmap = arr(refCoordMaps, c, Array) ;
          if (cmap) arrayDestroy(cmap) ;
        }
      arrayDestroy(refCoordMaps) ;
    }
  if (refSeqDict) dictDestroy(refSeqDict) ;
  if (sms) syncmerSetDestroy(sms) ;
  fprintf(stdout, "total: ") ; timeTotal(stdout) ;
  return 0 ;
}

/*********************** end of file **********************/
