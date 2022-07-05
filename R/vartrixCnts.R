#' Get counts at variants using vartrix
#'
#' Rather than using alleleCounter, with all the cruft that requires, use the dedicated 10X tool vartrix to get counts.  Unlike alleleCount, it will work with indels and other things just as well as single base pair events.
#'
#' Only BAMs and barcodes are allowed to be not length 1, in which case their lengths must match (if files) or barcodes must contain cellID_barcode formatted barcodes.   All BAMs are run with the same parameters in series.  The results are made non-ambiguous by the cellID, which is a concatenation of the BAM file label and the barcode, separated by an underscore.
#'
#' @param loci GRanges object that contains REF and ALT.
#' @param BAMs Cellranger 10X BAM files, or at least the least something that has the cellranger tags for UMI and CB.  Must be uniquely named by a label to form cell IDs.
#' @param barcodes The barcodes (or barcode files) to consider.
#' @param refGenome fasta reference genome to use.
#' @param outputs Path to save output for each BAM to.  Temp files used if NULL.
#' @param skipIfExists Should we load from file if possible.
#' @param mapq mapq threshold to use.
#' @param nBasesPad How many bases on either side of variant to use in mapping.
#' @param autoChr Automatically try and strip/add 'chr' to chromosome names as needed.
#' @param nParallel How many threads.
#' @param noOutputCheck Suppresses output warning, used for recursion.
#' @param ... Currently does nothing.
#' @return Loci, expanded to have one entry per cell/snp combination with reads, with refCount and altCount columns added.
#' @importFrom utils write.table read.delim read.table
#' @importFrom GenomeInfoDb seqnames seqinfo 
#' @importFrom GenomicRanges GRangesList
#' @importFrom BiocGenerics start
#' @importFrom methods as
#' @importFrom Rsamtools BamFile
#' @export
vartrixCnts = function(loci,BAMs,barcodes,refGenome,outputs=NULL,skipIfExists=TRUE,mapq=30,nBasesPad=500,autoChr=TRUE,nParallel=1,noOutputCheck=FALSE,...){
  #Check that things exist and are sensible.
  if(is.null(names(BAMs)) || any(duplicated(names(BAMs))))
     stop("BAM files lack unique names")
  if(is.null(outputs) && !noOutputCheck)
    warning("No output files specified.  These calculations are time consuming, it is a good idea to save the results.")
  #Store things to delete at the end
  trashHeap=c()
  if(length(BAMs)>1){
    #Check length of outputs if not null
    if(!is.null(outputs) && length(outputs)!=length(BAMs))
      stop("Length of outputs must match the length of BAMs")
    #Is barcodes just a vector of files?
    if(length(barcodes)==length(BAMs) && all(file.exists(barcodes))){
      bcodeFiles = barcodes
    }else{
      #Assume it's a vector of cellIDs, and try splitting
      grp = gsub('_[ACGT]+(-[0-9]+)?$','',barcodes)
      if(!all(grp %in% names(BAMs))|| !all(names(BAMs) %in% grp))
        stop("Could not resolve barcodes by cellID.  Check naming of BAMs and barcodes.")
      tmp = split(gsub('.*_([ACGT]+(-[0-9]+)?)$','\\1',barcodes),grp)
      #Make sure order matches BAM files
      tmp = tmp[names(BAMs)]
      bcodeFiles = sapply(BAMs,function(e) tempfile())
      trashHeap = append(trashHeap,bcodeFiles) 
      for(i in seq_along(tmp))
        write.table(tmp[[i]],bcodeFiles[i],sep='\t',row.names=FALSE,col.names=FALSE,quote=FALSE)
    }
    #Execute one at a time
    outs = list()
    for(i in seq_along(BAMs))
      outs[[i]] = vartrixCnts(loci,BAMs[i],bcodeFiles[i],refGenome,outputs[i],skipIfExists,mapq,nBasesPad,autoChr,nParallel,noOutputCheck=TRUE,...)
    #Merge and return
    unlink(trashHeap)
    return(unlist(GRangesList(outs)))
  }
  if(!file.exists(BAMs))
    stop("BAM files not found.")
  if(!file.exists(refGenome))
    stop("Reference genome not found.")
  #Does it already exist?
  if(is.null(outputs) || !file.exists(outputs) || !skipIfExists){
    #Nope, so calculate it.
    #If barcodes isn't already a file, make it one
    if(length(barcodes)==1 && is.character(barcodes) && file.exists(barcodes)){
      bcodesFile=barcodes
    }else{
      bcodesFile = tempfile()
      trashHeap = append(trashHeap,bcodesFile)
      #Make sure you trash the experiment ID if present
      write.table(gsub('.*_([ACGT]+(-[0-9]+)?)$','\\1',barcodes),bcodesFile,row.names=FALSE,col.names=FALSE,sep='\t',quote=FALSE)
    }
    bcodes = read.table(bcodesFile,sep='\t',header=FALSE)[,1]
    #Construct the VCF input
    tmpVCF = tempfile()
    trashHeap = append(trashHeap,tmpVCF)
    tmp = data.frame(paste0('chr',as.character(seqnames(loci))),
                     start(loci),
                     '.',
                     loci$REF,
                     loci$ALT,
                     '.',
                     '.',
                     '.')
    #Does this need syncing with BAM file chr status
    tmpGenome = refGenome
    if(autoChr){
      hasChr = seqinfo(BamFile(BAMs))
      hasChr = any(grepl('^chr',seqnames(hasChr)))
      vcfHasChr = any(grepl('^chr',tmp[,1]))
      if(hasChr!=vcfHasChr){
        #Need to add it in?
        if(hasChr){
          #Yes
          tmp[,1] = paste0('chr',tmp[,1])
        }else{
          #No, take it out
          tmp[,1] = gsub('^chr','',tmp[,1])
        }
      }
      #Now need to check (and correct) the genome reference
      genomeHasChr = system(sprintf("grep '>' %s",refGenome),intern=TRUE)
      genomeHasChr = any(grepl('^>chr',genomeHasChr))
      if(hasChr!=genomeHasChr){
        tmpGenome=tempfile()
        trashHeap = append(trashHeap,tmpGenome)
        #Do we need to add it?
        if(hasChr){
          #Yes
          system(sprintf("sed 's/^>/>chr/g' %s > %s",refGenome,tmpGenome))
        }else{
          #Nope, remove it
          system(sprintf("sed 's/^>chr/>/g' %s > %s",refGenome,tmpGenome))
        }
        #Need to index modified version
        system(sprintf('samtools faidx %s',tmpGenome))
        trashHeap = append(trashHeap,paste0(tmpGenome,'.fai'))
      }
    }
    #Actually write the vcf
    writeChar(paste0('##fileformat=VCFv4.3\n',
                     paste(sprintf('##contig=<ID=%s>',unique(tmp[,1])),collapse='\n'),
                     '\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'),tmpVCF,eos=NULL)
    write.table(tmp,tmpVCF,append=TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE,sep='\t')
    #Construct the command
    tmpRefMtx = tempfile()
    tmpAltMtx = tempfile()
    trashHeap = append(trashHeap,c(tmpRefMtx,tmpAltMtx))
    cmd = sprintf('vartrix --umi --mapq %d -b %s -c %s --scoring-method coverage --threads %d --ref-matrix %s --out-matrix %s -v %s --fasta %s -p %d',mapq,BAMs,bcodesFile,nParallel,tmpRefMtx,tmpAltMtx,tmpVCF,tmpGenome,nBasesPad)
    #message(cmd)
    system(cmd)
    #Load them in
    ref = read.delim(tmpRefMtx,sep=' ',header=FALSE,comment.char='%')[-1,]
    ref = data.frame(snpIdx = ref[,1],
                     barcode = bcodes[ref[,2]],
                     refCnt = ref[,3])
    alt = read.delim(tmpAltMtx,sep=' ',header=FALSE,comment.char='%')[-1,]
    alt = data.frame(snpIdx = alt[,1],
                     barcode = bcodes[alt[,2]],
                     altCnt = alt[,3])
    cnts = merge(ref,alt,by=c('snpIdx','barcode'),all=TRUE)
    #Save if outputs is not null
    if(!is.null(outputs))
      write.table(cnts,outputs,row.names=FALSE,col.names=TRUE,sep='\t',quote=FALSE)
  }else{
    cnts = read.table(outputs,sep='\t',header=TRUE)
  }
  out = loci[cnts$snpIdx]
  out$barcode = cnts$barcode
  out$refCount = cnts$refCnt
  out$altCount = cnts$altCnt
  out$Tot = cnts$refCnt+cnts$altCnt
  out$cellID = paste0(names(BAMs),'_',out$barcode)
  #Take out the trash
  unlink(trashHeap)
  return(out)
}
