library(grocaptss)
GROcapTSSHMM <- function(prefix_output, treatment_plus, treatment_minus, control_plus, control_minus, pairing_distance, CAGE_plus, CAGE_minus){
    bwSet_treatment = list(GROcap.plus=load.bigWig(treatment_plus),
                        GROcap.minus=load.bigWig(treatment_minus))
    bwSet_contorl = list(GROcap.plus = load.bigWig(control_plus),
                        GROcap.minus = load.bigWig(control_minus))
    cov_treatment = compute.normalization(bwSet_treatment)
    cov_control = compute.normalization(bwSet_contorl)
    scale.factor = (cov_treatment$GROcap.plus - cov_treatment$GROcap.minus) / (cov_control$GROcap.plus - cov_control$GROcap.minus)
    res = process.genome.qhmm(bwSet_treatment, bwSet_contorl, bwSet_treatment$GROcap.plus$chroms, scale.factor)
    preds = res$preds
    peaks = res$peaks
    peak.col.plus = do.call("paste", c(as.list(col2rgb(gcol$cage)), list(sep=',')))
    peak.col.minus = do.call("paste", c(as.list(col2rgb(gcol$groseq)), list(sep=',')))
    write.track.bed(preds, filename=paste(prefix_output, "new_hmm2b.bed", sep = "."), "hmm2b.preds")
    write.track.bed(peaks, filename=paste(prefix_output, "new_hmm2bp.bed", sep = "."), "hmm2bp.preds.peaks", color.plus = peak.col.plus, color.minus = peak.col.minus)
    preds = read.table(paste(prefix_output, "new_hmm2b.bed", sep = "."), skip=1, colClasses = c("factor", "integer", "integer", "factor", "integer", "factor", "integer", "integer", "factor"))
    peaks = read.table(paste(prefix_output, "new_hmm2bp.bed", sep = "."), skip=1)
    preds.filtered = filter.split(preds, peaks, scale.factor, 100)
    preds.filtered.post1 = edge.expand(preds.filtered, bwSet_treatment, bwSet_contorl, scale.factor)
    preds.filtered.post2b = edge.trim2(preds.filtered.post1, bwSet_treatment, bwSet_contorl, scale.factor)
    preds.filtered.post2a = edge.trim2(preds.filtered, bwSet_treatment, bwSet_contorl, scale.factor)
    pairs.150.b = create.pairs(preds.filtered.post2b, pairing_distance)
    pairs.150.a = create.pairs(preds.filtered.post2a, pairing_distance)
    write.bed(pairs.150.a$minus, paste(prefix_output, "new_hmm2bp.post2.pair_minus.bed", sep = "."))
    write.bed(pairs.150.a$plus, paste(prefix_output, "new_hmm2bp.post2.pair_plus.bed", sep = "."))

    pre_result <- cbind(pairs.150.a$plus, pairs.150.a$minus)
    colnames(pre_result) <- c("Chromosome", "Start_P1", "End_P1", 
                              "Status_P1", "Score_P1", "Strand_P1", 
                              "Chromosome_P2", "Start_P2", "End_P2", 
                              "Status_P2", "Score_P2", "Strand_P2")
    for(i in 1:nrow(pairs.150.a$plus)){
        row <- pre_result[i, ]
        start <- min(row$Start_P1, row$End_P1, row$Start_P2, row$End_P2)
        end <- max(row$Start_P1, row$End_P1, row$Start_P2, row$End_P2)
        # do stuff with row
        pre_result[i, "Start"] <- start
        pre_result[i, "End"] <- end
    }

    ord = order(pre_result[, 1], pre_result[, 2])
    pre_result = pre_result[ord, ]
    write.table(pre_result[, c("Chromosome", "Start", "End")], 
                paste(prefix_output, "bed", sep = "."), 
                sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0){
    stop("At least five args must be supplied", call.=FALSE)
}else if (length(args)==5){
    GROcapTSSHMM(args[[1]], args[[2]], args[[3]], args[[4]], args[[5]], 
        300,
        "/local/storage/ly349/projects/trivial/GROcap_hg38/data/cage_k562/k562_plus.bw", 
        "/local/storage/ly349/projects/trivial/GROcap_hg38/data/cage_k562/k562_minus.bw")
}else if (length(args)==6){
    GROcapTSSHMM(args[[1]], args[[2]], args[[3]], args[[4]], args[[5]], args[[6]],
        "/local/storage/ly349/projects/trivial/GROcap_hg38/data/cage_k562/k562_plus.bw", 
        "/local/storage/ly349/projects/trivial/GROcap_hg38/data/cage_k562/k562_minus.bw")
}else if (length(args)==8){
    GROcapTSSHMM(args[[1]], args[[2]], args[[3]], args[[4]], args[[5]], args[[6]], args[[7]], args[[8]])
}