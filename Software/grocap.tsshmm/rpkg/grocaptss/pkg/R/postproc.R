#

subpeak.split <- function(predictions, peaks, max.dist = 30) {
  chroms = NULL
  starts = NULL
  ends = NULL
  strands = NULL

  #
  peaks.chrom = lapply(levels(peaks[,1]), function(chrom) {
    list(plus = peaks[peaks[,1] == chrom & peaks[,6] == '+', 2:3],
         minus = peaks[peaks[,1] == chrom & peaks[,6] == '-', 2:3])
  })
  names(peaks.chrom) <- levels(peaks[,1])

  #
  foreach.bed(predictions, function(i, chrom, start, end, strand) {
    sname = "plus"
    if (strand == '-')
      sname = "minus"
    
    pchr = peaks.chrom[[chrom]][[sname]]

    #
    N = dim(pchr)[1]
    idxs = which(pchr[,1] >= start & pchr[,2] <= end)
    if (length(idxs) == 1) {
      # just trim ?
      #
      chroms <<- c(chroms, chrom)
      starts <<- c(starts, start)
      ends <<- c(ends, end)
      strands <<- c(strands, strand)
      
    } else {
      #
      # group peaks
      subset = pchr[idxs,]
      M = length(idxs)

      dists = subset[2:M,1] - subset[1:(M-1), 2]

      gidx = 1
      gIDs = vector(mode="integer", length=M)
      gIDs[1] = gidx

      for (i in 1:(M-1)) {
        if (dists[i] > max.dist) 
          gidx = gidx + 1
        gIDs[i + 1] = gidx
      }

      # add results
      chroms <<- c(chroms, rep(chrom, gidx))
      strands <<- c(strands, rep(strand, gidx))

      starts <<- c(starts, sapply(1:gidx, function(id)
                                  subset[which(gIDs == id)[1], 1]))
      ends <<- c(ends, sapply(1:gidx, function(id) {
        idxs = which(gIDs == id)
        subset[idxs[length(idxs)], 2]
      }))
    }
  })

  return(data.frame(chroms, starts, ends, "M2", 0, strands))
}

#' Trim edges
#'
#' @export
edge.trim2 <- function(predictions, bwSet, bwBck, scale.factor, step = 10) {
  result = NULL
  for (chrom in levels(predictions[,1])) {
    subset = predictions[predictions[,1] == chrom,]

    if (!is.null(dim(subset)[1]) && dim(subset)[1] > 0) {
      tap.reads.plus = chromStepSum.bigWig(bwSet$GROcap.plus, chrom, step, 0)
      tap.reads.minus = abs(chromStepSum.bigWig(bwSet$GROcap.minus, chrom, step, 0))
      notap.reads.plus = chromStepSum.bigWig(bwBck$GROcap.plus, chrom, step, 0)
      notap.reads.minus = abs(chromStepSum.bigWig(bwBck$GROcap.minus, chrom, step, 0))

      foreach.bed(subset, function(i, chrom, start, end, strand) {
        istart = (start + 1) %/% step + 1
        iend = end %/% step
        range = istart:iend

#        print(range)
        
        mask = NULL
        if (strand == '+')
          mask = tap.reads.plus[range] < notap.reads.plus[range] * scale.factor
        else
          mask = tap.reads.minus[range] < notap.reads.minus[range] * scale.factor
#        print(mask)

        # edge adjustment
        N = (end - start) %/% step
        n.left.drop = 0
        for (j in 1:N) {
          if (mask[j])
            n.left.drop = n.left.drop + 1
          else
            break
        }
        
        n.right.drop = 0
        for (j in N:1) {
          if (mask[j])
            n.right.drop = n.right.drop + 1
          else
            break
        }

#        print(n.left.drop)
#        print(n.right.drop)
        
        # set new values
        subset[i, 2] <<- start + n.left.drop * step
        subset[i, 3] <<- end - n.right.drop * step
      })

      result = rbind(result, subset)
    }
  }

  return(result)
}


#' Edge extension
#'
#' edge extension: for predictions smaller than X bp (100 bp) find overlapping (up to) X bp window
#'                 that maximally increases the number of stesp with TAP+ > TAP- (scaled).
#'                 Don't extend if increase is zero; don't add useless edges (TAP- (scaled) > TAP+ or
#'                 no data.), but these can be fixed on a subsequent edge-trim pass.
#'
#'                 Secundary criterion to resolve ties: total read count (TAP+)
#' @export
edge.expand <- function(predictions, bwSet, bwBck, scale.factor, step = 10, thresh = 100) {
  result = NULL
  
  for (chrom in levels(predictions[,1])) {
    subset = predictions[predictions[,1] == chrom,]

    if (!is.null(dim(subset)[1]) && dim(subset)[1] > 0) {
      tap.reads.plus = chromStepSum.bigWig(bwSet$GROcap.plus, chrom, step, 0)
      tap.reads.minus = abs(chromStepSum.bigWig(bwSet$GROcap.minus, chrom, step, 0))
      notap.reads.plus = chromStepSum.bigWig(bwBck$GROcap.plus, chrom, step, 0)
      notap.reads.minus = abs(chromStepSum.bigWig(bwBck$GROcap.minus, chrom, step, 0))

      foreach.bed(subset, function(i, chrom, start, end, strand) {
        if (end - start < thresh) {
          istart = (start + 1) %/% step + 1
          iend = end %/% step
          range = istart:iend

          iextra = (thresh - (end - start)) %/% step
          
          mask = NULL
          if (strand == '+')
            mask = tap.reads.plus[range] > notap.reads.plus[range] * scale.factor
          else
            mask = tap.reads.minus[range] > notap.reads.minus[range] * scale.factor

          # base score
          base.score = sum(mask)
          base.reads = sum(tap.reads.plus[range])
          if (strand == '-')
            base.reads = sum(tap.reads.minus[range])
          
          # scan/score alternatives
          scores = NULL
          reads = NULL
          for (before in 0:iextra) {
            after = iextra - before

            erange = (istart - before):(iend + after)
            #print(erange)
            emask = NULL
            if (strand == '+')
              emask = tap.reads.plus[erange] > notap.reads.plus[erange] * scale.factor
            else
              emask = tap.reads.minus[erange] > notap.reads.minus[erange] * scale.factor

            # score
            score = sum(emask) - base.score

            # reads
            reads.i = sum(tap.reads.plus[erange])
            if (strand == '-')
              reads.i = sum(tap.reads.minus[erange])

            scores = c(scores, score)
            reads = c(reads, reads.i - base.reads)
          }

          # pick solution
          if (any(scores > 0)) {
            ord = order(scores, reads, decreasing = TRUE)

            before = ord[1] - 1
            after = iextra - before

            subset[i, 2] <<- start - before * step
            subset[i, 3] <<- end + after * step
          }
        }
      })

      result = rbind(result, subset)
    }
  }

  return(result)
}

#' Write BED file
#'
#' @export
write.bed <- function(table, filename) {
  # make sure start/end are formated as integers
  M = dim(table)[2]
  bed = data.frame(table[,1], as.integer(table[,2]), as.integer(table[,3]),
    table[,4:M])
  # write file
  write.table(bed, file=filename, sep='\t', quote=F, col.names=F, row.names=F)
}

#' Filter/split
#'
#' @export
filter.split <- function(preds, peaks, scale.factor, thresh = 100) {
  sizes = preds[,3] - preds[,2]

  preds.small = preds[sizes < 100,]
  preds.large.m2 = preds[sizes >= 100 & preds[,4] == 'M2',]
  preds.large.m1 = preds[sizes >= 100 & preds[,4] == 'M1',]

  # filter/split large/M2
  preds.large.m2.split = subpeak.split(preds.large.m2, peaks, max.dist=20)
  colnames(preds.large.m2.split) <- paste("V", 1:6, sep='')

  # recombine
  return(rbind(preds.small[,1:6], preds.large.m1[,1:6], preds.large.m2.split))
}

#' Create pairs
#'
#' @export
create.pairs <- function(preds, thresh = 150) {
  pair.plus.150 = NULL
  pair.minus.150 = NULL
  N = dim(preds)[1]

  # sort
  ord = order(preds[,1], preds[,2])
  preds = preds[ord,]

  # look for pairs
  for (i in 2:N) {
    d = abs(preds[i, 2] - preds[i - 1, 3])
    
    if (preds[i, 1] == preds[i - 1, 1] &&
      preds[i - 1, 6] == '-' &&
      preds[i, 6] == '+' &&
      d < thresh) {

      pair.plus.150 = c(pair.plus.150, i)
      pair.minus.150 = c(pair.minus.150, i - 1)
    }
  }

  # create result
  bed.plus = preds[pair.plus.150,]
  bed.minus = preds[pair.minus.150,]

  return(list(plus = bed.plus, minus = bed.minus))
}

#' Select preds that have no neighbor up or down stream (up to thresh)
#'
#' @export
create.single <- function(preds, thresh = 250) {
  N = dim(preds)[1]

  # sort
  ord = order(preds[,1], preds[,2])
  preds = preds[ord,]

  # distances
  rlow = 1:(N-2)
  rmid = 2:(N-1)
  rhigh = 3:N
  dist.prev = preds[rmid, 2] - preds[rlow, 3]
  dist.next = preds[rhigh, 2] - preds[rmid, 3]

  # handle edges
  idxs = NULL
  if (preds[2, 2] - preds[1, 3] >= thresh)
    idxs = c(idxs, 1)
  if (preds[N, 2] - preds[N - 1, 3] >= thresh)
    idxs = c(idxs, N)

  # handle mid
  idxs.mid = 1 + which(dist.prev >= thresh & dist.next >= thresh)
  idxs = c(idxs, idxs.mid)

  return(preds[idxs,])
}
