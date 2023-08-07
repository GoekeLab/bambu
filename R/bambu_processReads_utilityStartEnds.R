findPolyATail = function(reads, softClipPrime, side, polyAPattern, polyAPatternLong, mat, searchLengthUnaligned = 20){
    if(side == "start"){
        polyASeqPrime <- GRanges(seqnames = c(1:length(reads)),
            IRanges(rep(1,length(reads)), width=softClipPrime))
    } else{
        polyASeqPrime <- GRanges(seqnames = c(1:length(reads)),
            IRanges(end = lengths(x), width=softClipPrime))
    }
    x = mcols(reads)$seq
    names(x) = c(1:length(reads))
    polyASeqPrime = BSgenome::getSeq(x, polyASeqPrime)
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
