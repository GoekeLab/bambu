findPolyATail = function(seqs, softClipPrime, side, polyAPattern, polyAPatternLong, mat, searchLengthUnaligned = 20){
    names(seqs) = c(1:length(seqs))
    if(side == "start"){
        polyASeqPrime <- GRanges(seqnames = c(1:length(seqs)),
            IRanges(rep(1,length(seqs)), width=softClipPrime))
    } else{
        polyASeqPrime <- GRanges(seqnames = c(1:length(seqs)),
            IRanges(end = lengths(seqs), width=softClipPrime))
    }
    polyASeqPrime = BSgenome::getSeq(seqs, polyASeqPrime)
    paPrimeTable <- barcodeAlignmentExtended(
        pattern=polyAPattern, 
        subject=polyASeqPrime,
        type='overlap',
        gapOpening=-2,
        gapExtension= -3,
        substitutionMatrix=mat)

    extendedSet <- which(paPrimeTable[,'score']> (0.95*searchLengthUnaligned))
    paPrimeTableLong <- barcodeAlignmentExtended(
        pattern=polyAPatternLong, 
        subject=polyASeqPrime[extendedSet], 
        type='overlap',
        gapOpening=-2,
        gapExtension= -3,
        substitutionMatrix=mat)
    paPrimeTable[extendedSet,] <- paPrimeTableLong
    return(paPrimeTable[,'score'])
}
