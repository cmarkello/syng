/*  File: syngannotate.c
 *  Author: cmarkello
 *  Copyright (C) 2026
 *-------------------------------------------------------------------
 * Description: Gene annotation tool for spliced pangenome and
 *   pantranscriptome construction using syncmer graphs and GBWT.
 *
 *   Implements algorithms inspired by vg rna
 *   (https://github.com/vgteam/vg/wiki/Transcriptomic-analyses)
 *   adapted to the syncmer-based pangenome framework of syng.
 *
 *   The tool:
 *   1. Reads an existing syncmer GBWT pangenome graph (.1gbwt + .1khash)
 *   2. Reads per-sample GFF/GTF annotation files
 *   3. For each sample's reference sequence, maps exon boundaries to
 *      syncmer node coordinates by scanning the reference through the
 *      kmer hash table
 *   4. Creates splice junction edges between consecutive exon boundaries
 *   5. Builds transcript paths through the syncmer graph, connecting
 *      exons via the mapped syncmer nodes
 *   6. Stores transcript paths in the GBWT for efficient compression,
 *      analogous to haplotype-specific transcripts (HSTs) in vg rna
 *   7. Outputs a spliced pangenome transcriptome (.1sgtxom) file
 *
 * HISTORY:
 * Created: Mar 27 2026
 *-------------------------------------------------------------------
 */

#include <pthread.h>

#include "seqio.h"
#include "seqhash.h"
#include "syng.h"
#include "gffio.h"

static OneSchema *schema ;

/**************** spliced syngbwt implementation ******************/

SplicedSyngBWT *splicedSyngBWTcreate (SyngBWT *baseBWT, SyncmerSet *sms)
{
  SplicedSyngBWT *ssbwt = new0(1, SplicedSyngBWT) ;
  ssbwt->baseBWT = baseBWT ;
  ssbwt->sms = sms ;
  ssbwt->spliceEdges = arrayCreate(4096, SyncSpliceEdge) ;
  ssbwt->txPaths = arrayCreate(1024, SyncTranscriptPath) ;
  ssbwt->txDict = dictCreate(1024) ;
  ssbwt->geneDict = dictCreate(512) ;
  ssbwt->nSamples = 0 ;
  ssbwt->sampleNames = 0 ;
  return ssbwt ;
}

void splicedSyngBWTdestroy (SplicedSyngBWT *ssbwt)
{
  if (!ssbwt) return ;
  int i ;
  if (ssbwt->spliceEdges) arrayDestroy(ssbwt->spliceEdges) ;
  for (i = 0 ; i < arrayMax(ssbwt->txPaths) ; ++i)
    { SyncTranscriptPath *tp = arrp(ssbwt->txPaths, i, SyncTranscriptPath) ;
      if (tp->transcriptId) free(tp->transcriptId) ;
      if (tp->geneId) free(tp->geneId) ;
      if (tp->nodes) newFree(tp->nodes, tp->nNodes, I32) ;
      if (tp->offsets) newFree(tp->offsets, tp->nNodes, U32) ;
      if (tp->exonBoundaries) newFree(tp->exonBoundaries, tp->nExons, I32) ;
    }
  if (ssbwt->txPaths) arrayDestroy(ssbwt->txPaths) ;
  if (ssbwt->txDict) dictDestroy(ssbwt->txDict) ;
  if (ssbwt->geneDict) dictDestroy(ssbwt->geneDict) ;
  for (i = 0 ; i < ssbwt->nSamples ; ++i)
    if (ssbwt->sampleNames[i]) free(ssbwt->sampleNames[i]) ;
  if (ssbwt->sampleNames) newFree(ssbwt->sampleNames, ssbwt->nSamples, char*) ;
  // do not destroy baseBWT or sms - they are owned externally
  newFree(ssbwt, 1, SplicedSyngBWT) ;
}

void splicedSyngBWTaddTranscriptPath (SplicedSyngBWT *ssbwt, SyncTranscriptPath *txPath)
{
  U64 idx ;
  if (!dictAdd(ssbwt->txDict, txPath->transcriptId, &idx))
    return ; // duplicate transcript, skip

  SyncTranscriptPath *tp = arrayp(ssbwt->txPaths, idx, SyncTranscriptPath) ;
  *tp = *txPath ; // copy the struct
  tp->transcriptId = strdup(txPath->transcriptId) ;
  tp->geneId = txPath->geneId ? strdup(txPath->geneId) : 0 ;

  // copy node and offset arrays
  tp->nodes = new(tp->nNodes, I32) ;
  memcpy(tp->nodes, txPath->nodes, tp->nNodes * sizeof(I32)) ;
  tp->offsets = new(tp->nNodes, U32) ;
  memcpy(tp->offsets, txPath->offsets, tp->nNodes * sizeof(U32)) ;
  if (tp->nExons > 0 && txPath->exonBoundaries)
    { tp->exonBoundaries = new(tp->nExons, I32) ;
      memcpy(tp->exonBoundaries, txPath->exonBoundaries, tp->nExons * sizeof(I32)) ;
    }

  if (tp->geneId) dictAdd(ssbwt->geneDict, tp->geneId, 0) ;
}

/**************** write spliced transcriptome ******************/

void splicedSyngBWTwrite (OneFile *of, SplicedSyngBWT *ssbwt)
{
  int i, j ;

  // write splice junctions
  for (i = 0 ; i < arrayMax(ssbwt->spliceEdges) ; ++i)
    { SyncSpliceEdge *se = arrp(ssbwt->spliceEdges, i, SyncSpliceEdge) ;
      oneInt(of, 0) = se->donorNode ;
      oneInt(of, 1) = se->acceptorNode ;
      oneInt(of, 2) = se->donorOffset ;
      oneInt(of, 3) = se->acceptorOffset ;
      oneWriteLine(of, 'J', 0, 0) ;
      oneInt(of, 0) = se->iSample ;
      oneWriteLine(of, 'a', 0, 0) ;
    }

  // write transcript paths
  for (i = 0 ; i < arrayMax(ssbwt->txPaths) ; ++i)
    { SyncTranscriptPath *tp = arrp(ssbwt->txPaths, i, SyncTranscriptPath) ;
      if (tp->nNodes == 0) continue ;

      oneInt(of, 0) = tp->iSample ;
      oneWriteLine(of, 'T', 0, 0) ;
      if (tp->transcriptId)
        oneWriteLine(of, 'i', strlen(tp->transcriptId), tp->transcriptId) ;
      if (tp->geneId)
        oneWriteLine(of, 'g', strlen(tp->geneId), tp->geneId) ;

      // write node list
      I64 *nodeList = new(tp->nNodes, I64) ;
      for (j = 0 ; j < tp->nNodes ; ++j) nodeList[j] = tp->nodes[j] ;
      oneWriteLine(of, 'n', tp->nNodes, nodeList) ;
      newFree(nodeList, tp->nNodes, I64) ;

      // write offset list
      I64 *offList = new(tp->nNodes, I64) ;
      for (j = 0 ; j < tp->nNodes ; ++j) offList[j] = tp->offsets[j] ;
      oneWriteLine(of, 'f', tp->nNodes, offList) ;
      newFree(offList, tp->nNodes, I64) ;

      // write exon boundaries
      if (tp->nExons > 0 && tp->exonBoundaries)
        { I64 *exonList = new(tp->nExons, I64) ;
          for (j = 0 ; j < tp->nExons ; ++j) exonList[j] = tp->exonBoundaries[j] ;
          oneWriteLine(of, 'x', tp->nExons, exonList) ;
          newFree(exonList, tp->nExons, I64) ;
        }

    }
}

/************ mapping genomic coordinates to syncmer nodes ************/

// SyncmerCoord: maps a genomic position range to a syncmer node
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

// Find the syncmer node closest to a given genomic coordinate.
// Uses binary search on the sorted coordMap.
static int findClosestSyncmer (Array coordMap, I64 targetPos, I32 *node, U32 *offset)
{
  if (arrayMax(coordMap) == 0) return 0 ;

  I64 lo = 0, hi = arrayMax(coordMap) - 1 ;
  while (lo < hi)
    { I64 mid = (lo + hi) / 2 ;
      if (arrp(coordMap, mid, SyncmerCoord)->genomePos < targetPos) lo = mid + 1 ;
      else hi = mid ;
    }

  // check both lo and lo-1 for closest
  SyncmerCoord *best = arrp(coordMap, lo, SyncmerCoord) ;
  if (lo > 0)
    { SyncmerCoord *prev = arrp(coordMap, lo-1, SyncmerCoord) ;
      I64 dPrev = targetPos - prev->genomePos ;
      I64 dCurr = best->genomePos - targetPos ;
      if (dPrev >= 0 && dPrev < dCurr) // prev is actually closer and before target
        best = prev ;
    }

  *node = best->syncNode ;
  *offset = (targetPos >= best->genomePos) ? (U32)(targetPos - best->genomePos) : 0 ;
  return 1 ;
}

// Binary search: find first index in coordMap where genomePos >= targetPos
static I64 coordMapLowerBound (Array coordMap, I64 targetPos)
{
  I64 lo = 0, hi = arrayMax(coordMap) ;
  while (lo < hi)
    { I64 mid = (lo + hi) / 2 ;
      if (arrp(coordMap, mid, SyncmerCoord)->genomePos < targetPos) lo = mid + 1 ;
      else hi = mid ;
    }
  return lo ;
}

/**************** main program ********************/

static char usage[] =
  "Usage: syngannotate [options] <input files>\n"
  "Build a spliced pangenome transcriptome from syncmer graphs and gene annotations.\n"
  "\n"
  "Required inputs:\n"
  "  -readK <.1khash file>    : syncmer hash table\n"
  "  -readGBWT <.1gbwt file>  : existing syncmer GBWT graph\n"
  "\n"
  "Annotation inputs (repeatable, one per sample):\n"
  "  -annotate <sample_name> <ref.fa> <annotation.gff/gtf>\n"
  "     Provide a sample name, its reference FASTA, and GFF/GTF annotation.\n"
  "     The reference FASTA is scanned for syncmers to map annotation\n"
  "     coordinates to graph nodes. Can be specified multiple times\n"
  "     for multiple samples/haplotypes.\n"
  "\n"
  "Options:\n"
  "  -w <window length>       : [55] syncmer length = w + k\n"
  "  -k <smer length>         : [8] must be under 32\n"
  "  -seed <seed>             : [7] for the hashing function\n"
  "  -o <outfile prefix>      : [syngAnnotate] output file prefix\n"
  "  -writeGBWT               : write transcript paths as .txpaths.1path (can be fed to syngpath2gbwt)\n"
  "  -T <threads>             : [8] number of threads\n" ;

// per-sample annotation specification
typedef struct {
  char *sampleName ;
  char *refFasta ;
  char *gffFile ;
} AnnotSpec ;

int main (int argc, char *argv[])
{
  char *outPrefix = "syngAnnotate" ;
  int   nThread = 8 ;
  SyncmerSet *sms = 0 ;
  SyngBWT    *gbwt = 0 ;
  SyncmerParams params = syncmerParamsDefault() ;
  bool  writeGBWT = false ;
  I64   i ;

  int maxAnnot = 32 ;
  int nAnnot = 0 ;
  AnnotSpec *annots = new0(maxAnnot, AnnotSpec) ;

  timeUpdate(0) ;
  schema = oneSchemaCreateFromText(syngSchemaText) ;

  storeCommandLine(argc, argv) ;
  --argc ; ++argv ;
  if (!argc) { fprintf(stderr, "%s", usage) ; exit(0) ; }

  while (argc > 0 && **argv == '-')
    if (!strcmp(*argv, "-w") && argc > 1)
      { params.w = atoi(argv[1]) ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp(*argv, "-k") && argc > 1)
      { params.k = atoi(argv[1]) ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp(*argv, "-seed") && argc > 1)
      { params.seed = atoi(argv[1]) ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp(*argv, "-T") && argc > 1)
      { nThread = atoi(argv[1]) ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp(*argv, "-o") && argc > 1)
      { outPrefix = argv[1] ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp(*argv, "-readK") && argc > 1)
      { sms = syncmerSetRead(argv[1]) ;
        if (!sms) die("failed to read syncmer set from %s", argv[1]) ;
        fprintf(stdout, "read syncmer set with %llu syncmers\n", kmerHashMax(sms->kh)) ;
        timeUpdate(stdout) ;
        argc -= 2 ; argv += 2 ;
      }
    else if (!strcmp(*argv, "-readGBWT") && argc > 1)
      { OneFile *of = oneFileOpenRead(argv[1], schema, "gbwt", nThread) ;
        if (!of) die("failed to open GBWT file %s", argv[1]) ;
        oneReadLine(of) ;
        if (of->lineType == 'h') syncmerParamsCheck(of, params) ;
        gbwt = syngBWTread(of) ;
        oneFileClose(of) ;
        fprintf(stdout, "read GBWT graph\n") ;
        timeUpdate(stdout) ;
        argc -= 2 ; argv += 2 ;
      }
    else if (!strcmp(*argv, "-writeGBWT"))
      { writeGBWT = true ; --argc ; ++argv ; }
    else if (!strcmp(*argv, "-annotate") && argc > 3)
      { if (nAnnot >= maxAnnot)
          { AnnotSpec *newA = new0(maxAnnot*2, AnnotSpec) ;
            memcpy(newA, annots, maxAnnot * sizeof(AnnotSpec)) ;
            newFree(annots, maxAnnot, AnnotSpec) ;
            annots = newA ;
            maxAnnot *= 2 ;
          }
        annots[nAnnot].sampleName = argv[1] ;
        annots[nAnnot].refFasta = argv[2] ;
        annots[nAnnot].gffFile = argv[3] ;
        ++nAnnot ;
        argc -= 4 ; argv += 4 ;
      }
    else
      die("unknown parameter %s\n%s", *argv, usage) ;

  if (!sms) die("must provide -readK <.1khash file>") ;
  if (!nAnnot) die("must provide at least one -annotate <sample> <ref.fa> <annotation.gff/gtf>") ;

  fprintf(stdout, "k, w, seed are %d %d %d\n", params.k, params.w, params.seed) ;
  Seqhash *sh = seqhashCreate(params.k, params.w+1, params.seed) ;

  // create spliced pangenome structure
  SplicedSyngBWT *ssbwt = splicedSyngBWTcreate(gbwt, sms) ;
  ssbwt->nSamples = nAnnot ;
  ssbwt->sampleNames = new0(nAnnot, char*) ;

  // process each sample's annotation
  I64 totalTranscripts = 0, totalSpliceJunctions = 0 ;

  for (i = 0 ; i < nAnnot ; ++i)
    { fprintf(stdout, "\nProcessing sample %d: %s\n", (int)(i+1), annots[i].sampleName) ;
      ssbwt->sampleNames[i] = strdup(annots[i].sampleName) ;

      // parse the GFF/GTF annotation
      GffFile *gf = gffFileOpen(annots[i].gffFile, i) ;
      if (!gf) die("failed to parse annotation file %s", annots[i].gffFile) ;

      // read the reference FASTA and build syncmer coordinate map
      fprintf(stdout, "  scanning reference %s for syncmer coordinates...\n", annots[i].refFasta) ;
      SeqIO *sio = seqIOopenRead(annots[i].refFasta, dna2indexConv, 0) ;
      if (!sio) die("failed to open reference FASTA %s", annots[i].refFasta) ;

      // build coordinate maps per chromosome and process transcripts
      DICT *seqNameDict = dictCreate(64) ;
      Array coordMaps = arrayCreate(64, Array) ;  // Array of Array of SyncmerCoord

      while (seqIOread(sio))
        { // get the chromosome name from the sequence
          char *seqName = strdup(sqioId(sio)) ;
          U64 seqIdx ;
          dictAdd(seqNameDict, seqName, &seqIdx) ;

          // convert to index4 form for syncmer scanning
          char *seqIdx4 = new(sio->seqLen, char) ;
          I64 j ;
          for (j = 0 ; j < sio->seqLen ; ++j)
            seqIdx4[j] = sqioSeq(sio)[j] ;

          Array cmap = buildSyncmerCoordMap(seqIdx4, sio->seqLen, sh, sms->kh) ;
          array(coordMaps, seqIdx, Array) = cmap ;
          newFree(seqIdx4, sio->seqLen, char) ;
          free(seqName) ;
        }
      seqIOclose(sio) ;
      fprintf(stdout, "  built syncmer coordinate maps for %lld sequences\n",
              (long long)dictMax(seqNameDict)) ;
      timeUpdate(stdout) ;

      // now build transcript paths using the coordinate maps
      I64 splicesBefore = arrayMax(ssbwt->spliceEdges) ;
      I64 txBefore = arrayMax(ssbwt->txPaths) ;

      int t ;
      for (t = 0 ; t < gffTranscriptCount(gf) ; ++t)
        { GffTranscript *tx = gffTranscript(gf, t) ;
          U64 seqIdx ;
          if (dictFind(seqNameDict, tx->seqid, &seqIdx) &&
              seqIdx < arrayMax(coordMaps))
            { Array cmap = arr(coordMaps, seqIdx, Array) ;
              if (cmap)
                { // build paths for just this transcript using a temp GffFile view
                  // We process transcripts one at a time against their chromosome's coordMap
                  Array pathNodes = arrayCreate(64, I32) ;
                  Array pathOffsets = arrayCreate(64, U32) ;
                  Array exonBounds = arrayCreate(tx->nExon, I32) ;
                  I64 lastPos = 0 ;
                  int j ;

                  for (j = 0 ; j < tx->nExon ; ++j)
                    { GffExon *exon = arrp(tx->exons, j, GffExon) ;
                      array(exonBounds, j, I32) = arrayMax(pathNodes) ;

                      // use binary search to find syncmers within this exon
                      I64 kStart = coordMapLowerBound(cmap, exon->start) ;
                      I64 k ;
                      for (k = kStart ; k < arrayMax(cmap) ; ++k)
                        { SyncmerCoord *sc = arrp(cmap, k, SyncmerCoord) ;
                          if (sc->genomePos > exon->end) break ;
                          I32 nd = (tx->strand == '-') ? -sc->syncNode : sc->syncNode ;
                          U32 off = (arrayMax(pathNodes) == 0) ? 0
                            : (U32)(sc->genomePos - lastPos) ;
                          arrayp(pathNodes, arrayMax(pathNodes), I32)[0] = nd ;
                          arrayp(pathOffsets, arrayMax(pathOffsets), U32)[0] = off ;
                          lastPos = sc->genomePos ;
                        }

                      // add splice junction edge between consecutive exons
                      if (j > 0)
                        { I32 donorNode = 0, acceptorNode = 0 ;
                          U32 donorOff = 0, acceptorOff = 0 ;
                          GffExon *prevExon = arrp(tx->exons, j-1, GffExon) ;

                          if (findClosestSyncmer(cmap, prevExon->end, &donorNode, &donorOff) &&
                              findClosestSyncmer(cmap, exon->start, &acceptorNode, &acceptorOff))
                            { SyncSpliceEdge *se = arrayp(ssbwt->spliceEdges,
                                                           arrayMax(ssbwt->spliceEdges),
                                                           SyncSpliceEdge) ;
                              se->donorNode = (tx->strand == '-') ? -donorNode : donorNode ;
                              se->acceptorNode = (tx->strand == '-') ? -acceptorNode : acceptorNode ;
                              se->donorOffset = donorOff ;
                              se->acceptorOffset = acceptorOff ;
                              se->transcriptIdx = t ;
                              se->iSample = i ;
                            }
                        }
                    }

                  // reverse path for minus strand transcripts
                  if (tx->strand == '-' && arrayMax(pathNodes) > 1)
                    { int nN = arrayMax(pathNodes) ;
                      I32 *nds = arrp(pathNodes, 0, I32) ;
                      U32 *ofs = arrp(pathOffsets, 0, U32) ;
                      for (j = 0 ; j < nN/2 ; ++j)
                        { I32 tmpN = nds[j] ; nds[j] = nds[nN-1-j] ; nds[nN-1-j] = tmpN ;
                          U32 tmpO = ofs[j] ; ofs[j] = ofs[nN-1-j] ; ofs[nN-1-j] = tmpO ;
                        }
                    }

                  if (arrayMax(pathNodes) > 0)
                    { SyncTranscriptPath txPath ;
                      memset(&txPath, 0, sizeof(txPath)) ;
                      txPath.transcriptId = tx->transcriptId ;
                      txPath.geneId = tx->geneId ;
                      txPath.iSample = i ;
                      txPath.nNodes = arrayMax(pathNodes) ;
                      txPath.nodes = arrp(pathNodes, 0, I32) ;
                      txPath.offsets = arrp(pathOffsets, 0, U32) ;
                      txPath.nExons = tx->nExon ;
                      txPath.exonBoundaries = arrp(exonBounds, 0, I32) ;

                      splicedSyngBWTaddTranscriptPath(ssbwt, &txPath) ;
                    }

                  arrayDestroy(pathNodes) ;
                  arrayDestroy(pathOffsets) ;
                  arrayDestroy(exonBounds) ;
                }
            }
        }

      I64 newSplices = arrayMax(ssbwt->spliceEdges) - splicesBefore ;
      I64 newTx = arrayMax(ssbwt->txPaths) - txBefore ;
      totalTranscripts += newTx ;
      totalSpliceJunctions += newSplices ;
      fprintf(stdout, "  added %lld transcript paths, %lld splice junctions for sample %s\n",
              (long long)newTx, (long long)newSplices, annots[i].sampleName) ;

      // clean up coordinate maps
      int c ;
      for (c = 0 ; c < arrayMax(coordMaps) ; ++c)
        { Array cmap = arr(coordMaps, c, Array) ;
          if (cmap) arrayDestroy(cmap) ;
        }
      arrayDestroy(coordMaps) ;
      dictDestroy(seqNameDict) ;
      gffFileClose(gf) ;
      timeUpdate(stdout) ;
    }

  fprintf(stdout, "\nTotal: %lld transcript paths, %lld splice junctions from %d samples\n",
          (long long)totalTranscripts, (long long)totalSpliceJunctions, nAnnot) ;

  // write the spliced transcriptome output
  char *fname = fnameTag(outPrefix, "1sgtxom") ;
  OneFile *ofOut = oneFileOpenWriteNew(fname, schema, "sgtxom", true, 1) ;
  if (!ofOut) die("failed to open output file %s", fname) ;
  oneAddProvenance(ofOut, "syngannotate", SYNG_VERSION, getCommandLine()) ;
  syncmerParamsWrite(ofOut, params) ;
  splicedSyngBWTwrite(ofOut, ssbwt) ;
  oneFileClose(ofOut) ;
  fprintf(stdout, "wrote spliced transcriptome to %s\n", fname) ;
  timeUpdate(stdout) ;

  // optionally write transcript paths as a separate GBWT-compatible path file
  if (writeGBWT)
    { fname = fnameTag(outPrefix, "txpaths.1path") ;
      ofOut = oneFileOpenWriteNew(fname, schema, "path", true, 1) ;
      if (!ofOut) die("failed to open transcript path output file %s", fname) ;
      oneAddProvenance(ofOut, "syngannotate", SYNG_VERSION, getCommandLine()) ;
      syncmerParamsWrite(ofOut, params) ;
      // write transcript paths as path objects with z/o lines
      int ti ;
      for (ti = 0 ; ti < arrayMax(ssbwt->txPaths) ; ++ti)
        { SyncTranscriptPath *tp = arrp(ssbwt->txPaths, ti, SyncTranscriptPath) ;
          if (tp->nNodes == 0) continue ;
          // compute total bp length from offsets
          I64 bpLen = 0 ;
          int j ;
          for (j = 0 ; j < tp->nNodes ; ++j) bpLen += tp->offsets[j] ;
          bpLen += sms->kh->len ; // add final syncmer length
          oneInt(ofOut, 0) = bpLen ;
          oneInt(ofOut, 1) = tp->iSample + 1 ;
          oneInt(ofOut, 2) = ti + 1 ;
          oneWriteLine(ofOut, 'P', 0, 0) ;
          I64 *nodeList = new(tp->nNodes, I64) ;
          I64 *offList = new(tp->nNodes, I64) ;
          for (j = 0 ; j < tp->nNodes ; ++j)
            { nodeList[j] = tp->nodes[j] ;
              offList[j] = tp->offsets[j] ;
            }
          oneWriteLine(ofOut, 'z', tp->nNodes, nodeList) ;
          oneWriteLine(ofOut, 'o', tp->nNodes, offList) ;
          newFree(nodeList, tp->nNodes, I64) ;
          newFree(offList, tp->nNodes, I64) ;
        }
      oneFileClose(ofOut) ;
      fprintf(stdout, "wrote %lld transcript paths to %s\n",
              (long long)arrayMax(ssbwt->txPaths), fname) ;
      timeUpdate(stdout) ;
    }

  // cleanup
  splicedSyngBWTdestroy(ssbwt) ;
  if (gbwt) syngBWTdestroy(gbwt) ;
  if (sms) syncmerSetDestroy(sms) ;
  seqhashDestroy(sh) ;
  newFree(annots, maxAnnot, AnnotSpec) ;

  fprintf(stdout, "total: ") ; timeTotal(stdout) ;
  return 0 ;
}

/*********************** end of file **********************/
