#
# GROcap TSS HMM
#

# port of HMM2 to new framework
#

#' Create TSS HMM
#'
#' @export
make.qhmm2 <- function(free.params.vec = c(
                         (1 - 0.9)/2, (1 - 0.9)/2,
                         0.9, 0.01,
                         0.45, 0.45)) {
  full.params.vec <- function(a, b) {
    c(1 - (a + b), a, b)
  }

  # Data: (Y, X)
  #
  # Y: 0 :: no peak
  #    1 :: peak
  #
  # X: 1 :: [no signap TAP+ = 0],
  #    2 :: [enriched TAP+ > TAP-],
  #    3 :: [depleated TAP- > TAP+ > 0]
  #    
  
  # states: B, M1.[1..3], M2.1, M2.2, M2.3 :: total 7
  #

  B = 1
  M1.1 = 2
  M1.2 = 3
  M1.3 = 4
  M2.1 = 5
  M2.2 = 6
  M2.3 = 7
  n.states = 7
  
  # valid transitions
  tm = matrix(data=0, nrow=n.states, ncol=n.states)

  tm[B,B] = 1
  tm[B, M1.1] = 2
  tm[B, M2.1] = 3

  tm[M1.1, M1.2] = 1

  tm[M1.2, M1.2] = 1
  tm[M1.2, M1.3] = 2

  tm[M1.3, B] = 1
  
  tm[M2.1, M2.1] = 1
  tm[M2.1, M2.2] = 2

  tm[M2.2, M2.1] = 1
  tm[M2.2, M2.2] = 2
  tm[M2.2, M2.3] = 3

  tm[M2.3, M2.3] = 1
  tm[M2.3, B] = 2

  #
  egrp = new.emission.groups(n.states, 2)
  egrp = add.emission.groups(egrp, group=c(2, 2:4))
  egrp = add.emission.groups(egrp, group=c(2, 5:7))
  
  #
  hmm <- new.qhmm(list(c(1,1), NULL),
                  tm, rep("discrete", n.states),
                  rep(list(c("discrete", "discrete")), n.states),
                  emission.groups = egrp)

  # set initial probabilities
  set.initial.probs.qhmm(hmm, c(1, rep(0, n.states - 1)))
  
  # set transition parameters
  set.transition.params.qhmm(hmm, B, c(0.99, 0.005, 0.005))
  set.transition.params.qhmm(hmm, M1.1, 1)
  set.transition.params.qhmm(hmm, M1.2, c(0.5, 0.5))
  set.transition.params.qhmm(hmm, M1.3, 1)
  set.transition.params.qhmm(hmm, M2.1, c(0.5, 0.5))
  set.transition.params.qhmm(hmm, M2.2, c(0.45, 0.1, 0.45))
  set.transition.params.qhmm(hmm, M2.3, c(0.5, 0.5))

  # set emission parameters & options

  # . Y
  Y.params = c(1, 0)
  Y.params.1 = c(0, 1)
  set.emission.params.qhmm(hmm, 1:n.states, Y.params, slot = 1, fixed = c(T, T))
  set.emission.params.qhmm(hmm, M2.2, Y.params.1, slot = 1, fixed = c(T, T))

  for (i in 1:n.states)
    set.emission.option.qhmm(hmm, i, "offset", 0, slot = 1)
  
  # . X
  X.B.params = full.params.vec(free.params.vec[1], free.params.vec[2])
  X.M1.params = full.params.vec(free.params.vec[3], free.params.vec[4])
  X.M2.params = full.params.vec(free.params.vec[5], free.params.vec[6])
  
  set.emission.params.qhmm(hmm, B, X.B.params, slot = 2)
  set.emission.params.qhmm(hmm, M1.1:M1.3, X.M1.params, slot = 2)
  set.emission.params.qhmm(hmm, M2.1:M2.3, X.M2.params, slot = 2)

  return(hmm)
}

#' Decode TSS HMM predictions
#'
decode.data.qhmm <- function (hmm, data, start = 0, step = 10,
                              m1.range = 2:4, m2.range = NA,
                              m2.start = 5, m2.mid = c(5, 6),
                              m2.end = 7) {
  
  path = viterbi.qhmm(hmm, data)
  m1.blocks = path.blocks.qhmm(path, m1.range)

  m2.blocks = NULL
  if (is.na(m2.range))
    m2.blocks = path.blocks2.qhmm(path, m2.start, m2.mid, m2.end)
  else
    m2.blocks = path.blocks.qhmm(path, m2.range)
  
  starts = as.integer((c(m1.blocks[1, ], m2.blocks[1, ]) - 1) * step + start)
  ends = as.integer(c(m1.blocks[2, ], m2.blocks[2, ]) * step + start)
  
  m1.dim = dim(m1.blocks)[2]
  if (is.null(m1.blocks))
    m1.dim = 0
  m2.dim = dim(m2.blocks)[2]
  if (is.null(m2.blocks))
    m2.dim = 0
  types = c(rep("M1", m1.dim), rep("M2", m2.dim))
  data.frame(starts, ends, types)
}

#' Process chromosome with TSS HMM
#'
process.chromosome.qhmm <- function (bwSet, bwBck, chrom, scale.factor, hmm = make.qhmm2(), enable.EM = FALSE, log2.thresh = 1, free.param.lst = NULL) {
  decode <- function(hmm, chrom, start, strand, data) {
    tmp = decode.data.qhmm(hmm, data, start = start)
    if (dim(tmp)[1] == 0) 
      return(NULL)
    else {
      bed = cbind(chrom, tmp, 0, strand)
      colnames(bed) <- c("chrom", "start", "end", "type", 
                         "score", "strand")
      return(bed)
    }
  }
  decode.peaks <- function(hmm, chrom, start, strand, data) {
    tmp = decode.data.qhmm(hmm, data, start = start, m2.range = 6)
    if (dim(tmp)[1] == 0) 
      return(NULL)
    else {
      bed = cbind(chrom, tmp, 0, strand)
      colnames(bed) <- c("chrom", "start", "end", "type", 
                         "score", "strand")
      return(bed[bed[,4] == 'M2',])
    }
  }
  combine.strands <- function(plus.lst, minus.lst) {
    if (length(plus.lst) == 1)
      return(rbind(plus.lst[[1]], minus.lst[[1]]))
    return(lapply(1:length(plus.lst), function(idx) {
      rbind(plus.lst[[idx]], minus.lst[[idx]])
    }))
  }

  em.trace.plus = NULL
  em.trace.minus = NULL
  em.params.plus = NULL
  em.params.minus = NULL
  start.params = collect.params.qhmm(hmm)

  N = max(1, length(free.param.lst))
  bed.plus.lst = vector(mode="list", length=N)
  bed.minus.lst = vector(mode="list", length=N)
  bed.peaks.plus.lst = vector(mode="list", length=N)
  bed.peaks.minus.lst = vector(mode="list", length=N)

  cat(chrom, "\n")
  cat(" * data (+)\n")
  tap.plus = chromStepSum.bigWig(bwSet$GROcap.plus, chrom, 
    step = 10, defaultValue = 0)
  notap.plus = chromStepSum.bigWig(bwBck$GROcap.plus, chrom, 
    step = 10, defaultValue = 0)
  res.plus = sequence.to.data(tap.plus, notap.plus, scale.factor, log2.thresh = log2.thresh)
  tap.plus = NULL
  notap.plus = NULL
  gc()

  cat(" * process (+)\n")
  if (is.null(free.param.lst)) {
    if (enable.EM) {
      em.trace.plus = em.qhmm(hmm, list(res.plus))
      em.params.plus = collect.params.qhmm(hmm)
    }
    bed.plus = decode(hmm, chrom, 0, "+", res.plus)
    bed.peaks.plus = decode.peaks(hmm, chrom, 0, "+", res.plus)
    res.plus = NULL
    if (enable.EM)
      restore.params.qhmm(hmm, start.params)
    gc()

    # store in list
    bed.plus.lst[[1]] = bed.plus
    bed.peaks.plus.lst[[1]] = bed.peaks.plus
  } else {
    for (i in 1:N) {
      # make new HMM
      free.params.i = free.param.lst[[i]]
      hmm = make.qhmm2(free.params.i)

      # parse data
      bed.plus = decode(hmm, chrom, 0, "+", res.plus)
      bed.peaks.plus = decode.peaks(hmm, chrom, 0, "+", res.plus)

      # store results
      bed.plus.lst[[i]] = bed.plus
      bed.peaks.plus.lst[[i]] = bed.peaks.plus
    }
    # clean up
    res.plus = NULL
    gc()
  }
  
  cat(" * data (-)\n")
  tap.minus = chromStepSum.bigWig(bwSet$GROcap.minus, chrom, 
    step = 10, defaultValue = 0)
  notap.minus = chromStepSum.bigWig(bwBck$GROcap.minus, chrom, 
    step = 10, defaultValue = 0)
  res.minus = sequence.to.data(tap.minus, notap.minus, scale.factor, log2.thresh = log2.thresh)
  tap.minus = NULL
  notap.minus = NULL
  gc()
  cat(" * process (-)\n")
  if (is.null(free.param.lst)) {
    if (enable.EM) {
      em.trace.minus = em.qhmm(hmm, list(res.minus))
      em.params.minus = collect.params.qhmm(hmm)
    }
    bed.minus = decode(hmm, chrom, 0, "-", res.minus)
    bed.peaks.minus = decode.peaks(hmm, chrom, 0, "-", res.minus)
    res.minus = NULL
    if (enable.EM)
      restore.params.qhmm(hmm, start.params)
    gc()

    # store in list
    bed.minus.lst[[1]] = bed.minus
    bed.peaks.minus.lst[[1]] = bed.peaks.minus
  } else {
    for (i in 1:N) {
      # make new HMM
      free.params.i = free.param.lst[[i]]
      hmm = make.qhmm2(free.params.i)

      # parse data
      bed.minus = decode(hmm, chrom, 0, "-", res.minus)
      bed.peaks.minus = decode.peaks(hmm, chrom, 0, "-", res.minus)

      # store results
      bed.minus.lst[[i]] = bed.minus
      bed.peaks.minus.lst[[i]] = bed.peaks.minus
    }
    # clean up
    res.minus = NULL
    gc()
  }

  #
  # prepare final result
  #
  if (length(bed.plus.lst)>0 && length(bed.minus.lst)>0){
      preds = combine.strands(bed.plus.lst, bed.minus.lst)
      peaks = combine.strands(bed.peaks.plus.lst, bed.peaks.minus.lst)
  }else{
      cat(paste("Failed to run viterbi on", chrom, "\n"))
      preds = data.frame(chrom=factor(), start=integer(), end=integer(), type=factor(), score=numeric(), strand=factor())
      peaks = data.frame(chrom=factor(), start=integer(), end=integer(), type=factor(), score=numeric(), strand=factor())
  }
  
  if (!enable.EM)
    return(list(preds = preds, peaks = peaks))
  else
    return(list(preds = preds, peaks = peaks,
  	        em = list(trace.plus = em.trace.plus,
	                  trace.minus = em.trace.minus,
		  	  params.plus = em.params.plus,
			  params.minus = em.params.minus)))
}

#' Process entire genome with TSS HMM
#'
#' @export
process.genome.qhmm <- function (bwSet, bwBck, chroms, scale.factor, hmm = make.qhmm2(), enable.EM = FALSE, log2.thresh = 1, free.param.lst = NULL) {
  combine.chroms <- function(curSet, new.lst) {
    if (is.null(curSet))
      return(new.lst)
    else {
      N = length(curSet)
      return(lapply(1:N, function(idx) {
        rbind(curSet[[idx]], new.lst[[idx]])
      }))
    }
  }
  
  bed = NULL
  peaks = NULL
  
  em.data = NULL

  for (chrom in chroms) {
    res.chrom = process.chromosome.qhmm(bwSet, bwBck, chrom, scale.factor, hmm = hmm,
      enable.EM = enable.EM, log2.thresh = log2.thresh, free.param.lst = free.param.lst)

    if (is.null(free.param.lst)) {
      bed.chrom = res.chrom$preds
      bed.peaks.chrom = res.chrom$peaks
    
      bed = rbind(bed, bed.chrom)
      peaks = rbind(peaks, bed.peaks.chrom)
      if (enable.EM) {
        em.data = c(em.data, list(res.chrom$em))
      }
    } else {
      bed.chrom.lst = res.chrom$preds
      bed.peaks.chrom.lst = res.chrom$peaks

      bed = combine.chroms(bed, bed.chrom.lst)
      peaks = combine.chroms(peaks, bed.peaks.chrom.lst)
    }
  }
  if (!enable.EM)
    return(list(preds = bed, peaks = peaks))
  else {
    names(em.data) <- chroms
    return(list(preds = bed, peaks = peaks, em = em.data))
  }
}
