#' Load genomic features and genome sequence
#'
#' This function loads the annotation created by the \code{prepare_annotation_files function}
#' @keywords RiboseQC
#' @author Lorenzo Calviello, \email{calviello.l.bio@@gmail.com}
#' @param path Full path to the *Rannot R file made by the \code{prepare_annotation_files} function (its value)
#' @param envir The environment to put \code{GTF_annotation} and \code{genome_seq} in; defaults to the environment of the caller
#' @details Reads the *Rannot files of all versions of RiboseQC and ORFquant: with the name of the BSgenome package in
#' \code{genome_package} (ORFquant, with \code{forge_BSgenome=TRUE}) or in \code{genome} (RiboseQC before 1.2.0), or with an
#' \code{FaFile_Circ} object in \code{genome} (all other files, and all files of RiboseQC 1.2.0 and later).
#' @return invisibly, a list with the annotated features (\code{GTF_annotation}) and the genome sequence (\code{genome_seq}).
#' The function also assigns the two objects to \code{GTF_annotation} and \code{genome_seq} in \code{envir}.
#' @seealso \code{\link{prepare_annotation_files}}
#' @examples
#' annot_file <- prepare_annotation_files(annotation_directory = tempdir(),
#'   gtf_file = system.file("ext_data", "simp_arab.gtf.gz", package = "RiboseQC"),
#'   genome_seq = system.file("ext_data", "simp_arab.fa.gz", package = "RiboseQC"))
#' load_annotation(annot_file)
#' GTF_annotation$genetic_codes
#' genome_seq
#' @export

load_annotation<-function(path,envir=parent.frame()){

    #a file that save() did not make, e.g. the GTF file, stops load() with "bad restore file magic number". The files
    #of save() start with RDA, RDB or RDX, the version and a newline
    con<-gzfile(path,"rb")
    magic<-readBin(con,"raw",5)
    close(con)
    if(!identical(magic[c(1,2,5)],charToRaw("RD\n"))){
        stop("The annotation file ",path," is not a file of prepare_annotation_files(): give the *_Rannot file that it makes, not the GTF file",call. = FALSE)
    }
    GTF_annotation<-get(load(path))
	  genome_package <- if(is.character(GTF_annotation$genome)) GTF_annotation$genome else GTF_annotation$genome_package
	  haspackage <- isTRUE(is.character(genome_package))
    where <- if(identical(envir,globalenv())) 'the global workspace' else 'the environment envir'
    if(haspackage){
      library(genome_package,character.only = TRUE)
      genome_sequence<-get(genome_package)
	    message(paste0('assigning genome package ',genome_package,' to ',where,' as genome_seq'))
      assign('genome_seq',genome_sequence,envir = envir)
    }else{
      genome_sequence<-GTF_annotation$genome
      message(paste0('assigning FaFile_Circ object ',GTF_annotation$genome$path,' to ',where,' as genome_seq'))
      assign('genome_seq',genome_sequence,envir = envir)
    }

	  message(paste0('assigning GTF_annotation object GTF_annotation to ',where))
    assign('GTF_annotation',GTF_annotation,envir = envir)
    invisible(list(GTF_annotation=GTF_annotation,genome_seq=genome_sequence))
}


#' Prepare comprehensive sets of annotated genomic features
#'
#' This function processes a gtf file and a FASTA file of the genome to create a comprehensive set of genomic regions of interest in genomic and transcriptomic space (e.g. introns, UTRs, start/stop codons).
#'    In addition, by linking genome sequence and annotation, it extracts additional info, such as gene and transcript biotypes, genetic codes for different organelles, or chromosomes and transcripts lengths.
#' @keywords RiboseQC
#' @author Lorenzo Calviello, \email{calviello.l.bio@@gmail.com}
#' @param annotation_directory The target directory which will contain the output files
#' @param twobit_file Not used. Earlier versions forged a \code{BSgenome} package from this twobit file
#' @param gtf_file Full path to the annotation file in GTF format
#' @param scientific_name Not used. Earlier versions named the forged \code{BSgenome} package with it
#' @param annotation_name Not used. Earlier versions named the forged \code{BSgenome} package with it
#' @param export_bed_tables_TxDb Export coordinates and info about different genomic regions in the annotation_directory? It defaults to \code{TRUE}
#' @param forge_BSgenome Not used. Earlier versions forged and installed a \code{BSgenome} package from \code{twobit_file} when it was \code{TRUE};
#' this failed with current BSgenomeForge versions (ohlerlab/RiboseQC#21, #23). It defaults to \code{FALSE}
#' @param genome_seq FASTA file of the genome sequence, uncompressed or compressed with \code{bgzip} (a path or an \code{FaFile}). Required; \code{twoBitToFa} of the UCSC tools
#' (http://hgdownload.soe.ucsc.edu/admin/exe/) makes one from a twobit file. The annotation refers to the FASTA file by its path:
#' a relative path is stored as an absolute path, the path of an \code{FaFile} as it is. If the FASTA file has no index (.fai, and .gzi
#' with \code{bgzip}), it is made next to the file
#' @param circ_chroms Chromosomes to make circular in the genome sequence - defaults to "chrM", "MT", "MtDNA", "mit", "Mito",
#' "mitochondrion", "dmel_mitochondrion_genome", "Pltd", "ChrC", "Pt", "chloroplast", "Chloro", "2micron", "2-micron", "2uM", "Mt",
#' "NC_001879.2", "NC_006581.1", "ChrM" and "mitochondrion_genome". \code{RiboseQC_analysis} analyzes each circular chromosome
#' as its own compartment
#' @param create_TxDb Create a \code{TxDb} object and a *Rannot object? It defaults to \code{TRUE}
#' @param annot_file specify an exact file name for the *Rannot R file (made with \code{save}) created by this function, defaults to annotation_directory/basename(gtf)_Rannot
#' @param genetic_codes The genetic codes of some chromosomes, which replace the codes that their names give (see Details): a named
#' vector, with the chromosome names as names and ids of \code{\link[Biostrings]{GENETIC_CODE_TABLE}} as values, e.g.
#' \code{c(chrM = "3")} for the mitochondrion of yeast with UCSC names. It defaults to \code{NULL}: the names give all codes
#' @details This function uses the \code{makeTxDbFromGFF} function to  create a TxDb object and extract
#' genomic regions and other info to a *Rannot R file; the \code{pmapToTranscripts} and \code{pmapFromTranscripts} functions are used to
#' map features to genomic or transcript-level coordinates. GTF file must contain "exon" and "CDS" lines,
#' where each line contains "transcript_id" and "gene_id" values. The CDS should include the stop codon, or the file should have
#' "stop_codon" lines, which are added to the CDS (as in GENCODE, Ensembl and RefSeq files). If most CDS end before their stop codon
#' (e.g. CDS lines made by factR, or a GENCODE file without its "stop_codon" lines), \code{stop_in_gtf} is \code{NA}, and
#' ORFquant's \code{annotate_ORFs} uses the codon after each CDS as its stop codon. It does so for all CDS of the file, so an ORF at the
#' stop codon of a CDS that includes it gets a wrong category, e.g. "C_truncation" instead of "ORF_annotated". If \code{stop_in_gtf} is
#' \code{"*"}, ORFquant 1.3.1 and later keep the last codon of a CDS that ends without a stop codon (e.g. with the GENCODE tag
#' "cds_end_NF"); earlier versions remove it, so an ORF at its end gets a wrong category too. Biotypes and gene names are read
#' from "gene_biotype" or "gene_type", "transcript_biotype" or "transcript_type", and "gene_name", "gene_symbol", "gene" (NCBI) or
#' "ref_gene_name" (StringTie) values, on any line of the transcript or (for genes) of the gene. Missing biotypes are "no_type", apart
#' from those of transcripts with CDS lines and of their genes, which are "protein_coding" (without biotypes, e.g. in UCSC's GTFs,
#' transcripts with a CDS that GENCODE calls "nonsense_mediated_decay" are "protein_coding" too). The transcript biotype "mRNA" (NCBI)
#' is read as "protein_coding"; if no gene has a name, all are "no_name". The CDS of transcripts with a CDS shorter than 3 nt are not used:
#' these transcripts are not in \code{cds_txs}, \code{cds_txs_coords} and \code{start_stop_codons}.
#' The order of the lines in the GTF file does not change the results if the lines agree: all lines of a transcript give the same
#' values, and all transcripts of a gene the same gene name and biotype (as in GENCODE and Ensembl files). If they do not agree, the
#' order chooses the value: a transcript gets the value of its first line with one, and the read counts of the genes and the top 50
#' positions of \code{RiboseQC_analysis} get the gene name and biotype of the first transcript of the gene in the file. Then this
#' function gives a warning, with some of these transcripts and genes.
#' The genetic code of a chromosome comes from its name: "chrM", "MT", "MtDNA", "mit" and "mitochondrion" get code 2 (vertebrate
#' mitochondria), "Mito" gets code 3 (yeast mitochondria), "dmel_mitochondrion_genome" gets code 5 (invertebrate mitochondria), and all
#' other chromosomes get code 1 (standard). So the mitochondrion of e.g. yeast or Drosophila with UCSC names ("chrM") gets code 2: give
#' its code with \code{genetic_codes}. \code{RiboseQC_analysis} uses the code of each circular chromosome for the codons of its CDS,
#' and code 1 for the codons of the other chromosomes. ORFquant uses the code of each chromosome to find and translate its ORFs.
#' The genome sequence is read from the FASTA file \code{genome_seq}.\cr\cr
#' The resulting GTF_annotation object (obtained after runnning \code{load_annotation}) contains:\cr\cr
#' \code{txs}: annotated transcript boundaries.\cr
#' \code{txs_gene}: GRangesList including transcript grouped by gene.\cr
#' \code{seqinfo}: indicating chromosomes and chromosome lengths.\cr
#' \code{start_stop_codons}: the set of annotated start and stop codon, with respective transcript and gene_ids.
#' reprentative_mostcommon,reprentative_boundaries and reprentative_5len represent the most common start/stop codon,
#' the most upstream/downstream start/stop codons and the start/stop codons residing on transcripts with the longest 5'UTRs\cr
#' \code{cds_txs}: GRangesList including CDS grouped by transcript.\cr
#' \code{introns_txs}: GRangesList including introns grouped by transcript.\cr
#' \code{cds_genes}: GRangesList including CDS grouped by gene.\cr
#' \code{exons_txs}: GRangesList including exons grouped by transcript.\cr
#' \code{exons_bins}: the list of exonic bins with associated transcripts and genes.\cr
#' \code{junctions}: the list of annotated splice junctions, with associated transcripts and genes.\cr
#' \code{genes}: annotated genes coordinates, without genes with exons on both strands or on more than one chromosome.\cr
#' \code{threeutrs}: collapsed set of 3'UTR regions, with correspinding gene_ids. This set does not overlap CDS region.\cr
#' \code{fiveutrs}: collapsed set of 5'UTR regions, with correspinding gene_ids. This set does not overlap CDS region.\cr
#' \code{ncIsof}: collapsed set of exonic regions of protein_coding genes, with correspinding gene_ids. This set does not overlap CDS region.\cr
#' \code{ncRNAs}: collapsed set of exonic regions of non_coding genes, with correspinding gene_ids. This set does not overlap CDS region.\cr
#' \code{introns}: collapsed set of intronic regions, with correspinding gene_ids. This set does not overlap exonic region.\cr
#' \code{intergenicRegions}: set of intergenic regions, defined as regions with no gene of \code{genes} on either strand.\cr
#' \code{trann}: DataFrame object including the mapping between gene_id, gene_biotype, gene_name, transcript_id and transcript_biotype, for each transcript.\cr
#' \code{cds_txs_coords}: transcript-level coordinates of ORF boundaries, for each annotated coding transcript (the seqnames are the transcript ids).
#' Additional columns are gene_id, upstr_stasto and mostfreq_stasto (how many of the start and stop codons of the transcript are the most
#' upstream/downstream ones or the most common ones of the gene), lentx (the length of the transcript), and the reprentative_mostcommon,
#' reprentative_boundaries and reprentative_5len of the \code{start_stop_codons} object, for the transcript.\cr
#' \code{genetic_codes}: an object containing the list of genetic code ids used for each chromosome/organelle. see GENETIC_CODE_TABLE for more info.\cr
#' \code{genome_package}: \code{NULL}. In annotations made by ORFquant with \code{forge_BSgenome=TRUE}, the name of the forged
#' BSgenome package, which \code{load_annotation} loads.\cr
#' \code{stop_in_gtf}: \code{"*"} if most CDS end with a stop codon (of the genetic code of their chromosome), \code{NA} if most end
#' before it. The CDS on circular chromosomes count only if all CDS are on them.\cr
#' \code{genome}: the genome sequence, an FaFile_Circ object. Loaded with \code{load_annotation} function.\cr
#' @return The path of the *Rannot file (\code{annot_file}). The function creates this file and a TxDb file,
#' \code{<annotation_directory>/<basename(gtf_file)>_TxDb}. With \code{export_bed_tables_TxDb = TRUE}, it also creates these tab-separated files
#' in \code{annotation_directory}: fiveutrs_similbed.bed, threeutrs_similbed.bed, ncIsof_similbed.bed, ncRNAs_similbed.bed, introns_similbed.bed,
#' cds_txs_coords_similbed.bed, table_gene_tx_IDs, seqinfo and genetic_codes. With \code{create_TxDb = FALSE}, it creates no file.
#' @seealso \code{\link{load_annotation}}, \code{\link[txdbmaker]{makeTxDbFromGFF}}.
#' @examples
#' annot_file <- prepare_annotation_files(annotation_directory = tempdir(),
#'   gtf_file = system.file("ext_data", "simp_arab.gtf.gz", package = "RiboseQC"),
#'   genome_seq = system.file("ext_data", "simp_arab.fa.gz", package = "RiboseQC"))
#' annot_file
#' @export


prepare_annotation_files<-function(annotation_directory,twobit_file=NULL,gtf_file,scientific_name="Homo.sapiens",
                                   annotation_name="genc25",export_bed_tables_TxDb=TRUE,forge_BSgenome=FALSE,genome_seq=NULL,circ_chroms=DEFAULT_CIRC_SEQS,create_TxDb=TRUE,annot_file=NULL,genetic_codes=NULL){
    #forging a BSgenome package from twobit_file is removed: it fails with current BSgenomeForge (ohlerlab/RiboseQC#21, #23)
    if(is.null(genome_seq)){
        stop("Please give the genome sequence as a FASTA file (genome_seq). prepare_annotation_files() no longer forges a BSgenome package ",
             "from a twobit file (forge_BSgenome, twobit_file); twoBitToFa of the UCSC tools makes a FASTA file from it",call. = FALSE)
    }
    if(forge_BSgenome){
        message('fasta file passed - cancelling BSgenome creation')
    }

    DEFAULT_CIRC_SEQS <- unique(c("chrM","MT","MtDNA","mit","Mito","mitochondrion",
                                  "dmel_mitochondrion_genome","Pltd","ChrC","Pt","chloroplast",
                                  "Chloro","2micron","2-micron","2uM",
                                  "Mt", "NC_001879.2", "NC_006581.1","ChrM","mitochondrion_genome"))

    filestotest <- c(gtf_file)
    if(is.character(genome_seq)) filestotest <- c(filestotest,genome_seq)
    check_files_exist(filestotest)
    if(create_TxDb) check_gtf(gtf_file)

    if(!dir.exists(annotation_directory)){dir.create(path = annotation_directory,recursive = TRUE)}
    annotation_directory<-normalizePath(annotation_directory)
    gtf_file<-normalizePath(gtf_file)
    #the annotation stores the FASTA path: make it absolute, but keep the file name, so that the .fai and .gzi next to it are used
    if(is.character(genome_seq)) genome_seq<-paste(normalizePath(dirname(genome_seq)),basename(genome_seq),sep="/")


    if(!is(genome_seq,'FaFile')){
        genome_seq <- Rsamtools::FaFile(genome_seq)
    }
    if(!is(genome_seq,'FaFile_Circ')){
        genome_seq <- FaFile_Circ(genome_seq,circularRanges=circ_chroms)
    }
    seqinfo_genome<-seqinfo(genome_seq)
    seqinfo_genome@is_circular[which(seqnames(seqinfo_genome)%in%circ_chroms)]<-TRUE
    check_genetic_codes(genetic_codes,seqnames(seqinfo_genome))
    genome <- genome_seq
    pkgnm<-NULL
    #Create the TxDb from GTF and BSGenome info

    if(is.null(annot_file)){
        annot_file <- paste(annotation_directory,"/",basename(gtf_file),"_Rannot",sep="")
    }

    if(create_TxDb){
        message("Creating the TxDb object ... ",date())

        annotation<-tryCatch(txdbmaker::makeTxDbFromGFF(file=gtf_file,format="gtf",chrominfo = seqinfo_genome),
                             error=function(e){
                                 stop("Reading the GTF file ",gtf_file," failed: ",conditionMessage(e),"\nThe GTF file needs exon and CDS lines with transcript_id and gene_id, ",
                                      "and its chromosome names must be those of the genome sequence (e.g. ",paste(head(seqnames(seqinfo_genome),3),collapse=", "),")",call. = FALSE)
                             })
        if(length(GenomicFeatures::cds(annotation))==0){
            stop("The GTF file has no CDS lines: the annotation needs the coding sequences (CDS) of the transcripts")
        }

        saveDb(annotation, file=paste(annotation_directory,"/",basename(gtf_file),"_TxDb",sep=""))
        message("Creating the TxDb object --- Done! ",date())
        message("Extracting genomic regions ... ",date())


        genes<-genes(annotation)
        exons_ge<-exonsBy(annotation,by="gene")
        exons_ge<-GenomicRanges::reduce(exons_ge)

        cds_gen<-cdsBy(annotation,"gene")
        cds_ge<-reduce(cds_gen)

        #define regions not overlapping CDS ( or exons when defining introns)

        threeutrs<-reduce(GenomicRanges::setdiff(unlist(threeUTRsByTranscript(annotation)),unlist(cds_ge),ignore.strand=FALSE))

        fiveutrs<-reduce(GenomicRanges::setdiff(unlist(fiveUTRsByTranscript(annotation)),unlist(cds_ge),ignore.strand=FALSE))

        introns<-reduce(GenomicRanges::setdiff(unlist(intronsByTranscript(annotation)),unlist(exons_ge),ignore.strand=FALSE))

        nc_exons<-reduce(GenomicRanges::setdiff(unlist(exons_ge),reduce(c(unlist(cds_ge),fiveutrs,threeutrs)),ignore.strand=FALSE))

        #assign gene ids (mutiple when overlapping multiple genes; none for regions of genes that genes() drops,
        #with exons on both strands or on several chromosomes, e.g. UCSC's PAR genes)
        ov<-findOverlaps(threeutrs,genes)
        ov<-split(subjectHits(ov),factor(queryHits(ov),levels=seq_along(threeutrs)))
        threeutrs$gene_id<-CharacterList(lapply(ov,FUN = function(x){names(genes)[x]}))
        ov<-findOverlaps(fiveutrs,genes)
        ov<-split(subjectHits(ov),factor(queryHits(ov),levels=seq_along(fiveutrs)))
        fiveutrs$gene_id<-CharacterList(lapply(ov,FUN = function(x){names(genes)[x]}))
        ov<-findOverlaps(introns,genes)
        ov<-split(subjectHits(ov),factor(queryHits(ov),levels=seq_along(introns)))
        introns$gene_id<-CharacterList(lapply(ov,FUN = function(x){names(genes)[x]}))
        ov<-findOverlaps(nc_exons,genes)
        ov<-split(subjectHits(ov),factor(queryHits(ov),levels=seq_along(nc_exons)))
        nc_exons$gene_id<-CharacterList(lapply(ov,FUN = function(x){names(genes)[x]}))

        intergenicRegions<-genes
        strand(intergenicRegions)<-"*"
        intergenicRegions <- gaps(reduce(intergenicRegions))
        intergenicRegions<-intergenicRegions[strand(intergenicRegions)=="*"]

        cds_tx<-cdsBy(annotation,"tx",use.names=TRUE)

        #filter out abnormally short cds (I"m looking at you maize annotation)
        cds_tx <- cds_tx[sum(width(cds_tx))>=3]
        txs_gene<-transcriptsBy(annotation,by="gene")

        exons_tx<-exonsBy(annotation,"tx",use.names=TRUE)

        transcripts_db<-transcripts(annotation)
        intron_names_tx<-intronsByTranscript(annotation,use.names=TRUE)


        #define exonic bins, including regions overlapping multiple genes
        nsns<-exonicParts(annotation, linked.to.single.gene.only=F)



        #define tx_coordinates of ORF boundaries

        exsss_cds<-exons_tx[names(cds_tx)]
        chunks<-seq(1,length(cds_tx),by = 20000)
        if(length(chunks)==1 || chunks[length(chunks)]<length(cds_tx)){chunks<-c(chunks,length(cds_tx))}
        mapp<-GRangesList()
        for(i in seq_len(length(chunks)-1)){
            if(i!=(length(chunks)-1)){
                mapp<-suppressWarnings(c(mapp,pmapToTranscripts(cds_tx[chunks[i]:(chunks[i+1]-1)],transcripts = exsss_cds[chunks[i]:(chunks[i+1]-1)])))
            }
            if(i==(length(chunks)-1)){
                mapp<-suppressWarnings(c(mapp,pmapToTranscripts(cds_tx[chunks[i]:(chunks[i+1])],transcripts = exsss_cds[chunks[i]:(chunks[i+1])])))
            }
        }
        cds_txscoords<-unlist(mapp)


        #extract biotypes and ids, one row per transcript: each value comes from the first line of the transcript
        #with it (e.g. gffread writes biotypes on transcript lines only), and for genes from a line of the gene
        #(e.g. NCBI's gene lines, with transcript_id ""). gene_type and transcript_type are GENCODE's names,
        #gene NCBI's and ref_gene_name StringTie's

        message("Extracting ids and biotypes ... ",date())

        gtf_ids<-data.frame(unique(mcols(import.gff2(gtf_file,colnames=c("gene_id","gene_biotype","gene_type","gene_name","gene_symbol","gene","ref_gene_name","transcript_id","transcript_biotype","transcript_type")))),stringsAsFactors=FALSE)
        gtf_ids$transcript_id[gtf_ids$transcript_id%in%""]<-NA
        first_value<-function(cols,by,ids){
            res<-rep(NA_character_,length(ids))
            for(cl in cols){
                ok<-!is.na(gtf_ids[,cl]) & !is.na(gtf_ids[,by])
                res[is.na(res)]<-gtf_ids[ok,cl][match(ids[is.na(res)],gtf_ids[ok,by])]
            }
            res
        }
        txs_ids<-unique(gtf_ids$transcript_id[!is.na(gtf_ids$transcript_id)])
        trann<-data.frame(gene_id=first_value("gene_id","transcript_id",txs_ids),stringsAsFactors=FALSE)
        gene_cols<-list(gene_biotype=c("gene_biotype","gene_type"),gene_name=c("gene_name","gene_symbol","gene","ref_gene_name"))
        for(cl in names(gene_cols)){
            trann[,cl]<-first_value(gene_cols[[cl]],"transcript_id",txs_ids)
            miss<-is.na(trann[,cl])
            trann[miss,cl]<-first_value(gene_cols[[cl]],"gene_id",trann$gene_id[miss])
        }
        trann$transcript_id<-txs_ids
        trann$transcript_biotype<-first_value(c("transcript_biotype","transcript_type"),"transcript_id",txs_ids)

        trann$gene_biotype[is.na(trann$gene_biotype)]<-"no_type"
        trann$transcript_biotype[is.na(trann$transcript_biotype)]<-"no_type"
        #NCBI's coding transcripts are "mRNA"
        trann$transcript_biotype[trann$transcript_biotype=="mRNA"]<-"protein_coding"
        #without biotypes (e.g. UCSC's GTFs), transcripts with a CDS and their genes are taken as protein_coding
        tx_cds<-trann$transcript_id%in%names(cds_tx)
        trann$transcript_biotype[trann$transcript_biotype=="no_type" & tx_cds]<-"protein_coding"
        trann$gene_biotype[trann$gene_biotype=="no_type" & trann$gene_id%in%trann$gene_id[tx_cds]]<-"protein_coding"
        if(all(is.na(trann$gene_name))){trann$gene_name<-"no_name"}
        check_gtf_values(gtf_ids,trann)

        trann<-DataFrame(trann)



        #introns and transcript_ids/gene_ids
        unq_intr<-sort(unique(unlist(intron_names_tx)))
        names(unq_intr)<-NULL
        all_intr<-unlist(intron_names_tx)

        ov<-findOverlaps(unq_intr,all_intr,type="equal")
        ov<-split(subjectHits(ov),queryHits(ov))
        a_nam<-CharacterList(lapply(ov,FUN = function(x){unique(names(all_intr)[x])}))

        unq_intr$type=rep("J",length(unq_intr))
        unq_intr$tx_name<-a_nam


        mat_genes<-match(unq_intr$tx_name,trann$transcript_id)
        g<-unlist(apply(cbind(seq_along(mat_genes),Y = elementNROWS(mat_genes)),FUN =function(x) rep(x[1],x[2]),MARGIN = 1))
        g2<-split(trann[unlist(mat_genes),"gene_id"],g)
        unq_intr$gene_id<-CharacterList(lapply(g2,unique))


        #filter ncRNA and ncIsof regions
        ncrnas<-nc_exons[!nc_exons%over%genes[names(genes)%in%trann$gene_id[trann$gene_biotype=="protein_coding"]]]
        ncisof<-nc_exons[nc_exons%over%genes[names(genes)%in%trann$gene_id[trann$gene_biotype=="protein_coding"]]]


        # define genetic codes to use
        # IMPORTANT : modify if needed (e.g. different organelles or species) check ids of GENETIC_CODE_TABLE for more info

        ifs<-seqinfo(annotation)
        # close the SQLite connection now; a finalizer run during S4 method lookup fails (lcalviell/ORFquant#3)
        annotation$finalize()
        translations<-as.data.frame(ifs)
        translations$genetic_code<-"1"

        #insert new codes for chromosome name

        #Mammalian mito
        translations$genetic_code[rownames(translations)%in%c("chrM","MT","MtDNA","mit","mitochondrion")]<-"2"

        #Yeast mito
        translations$genetic_code[rownames(translations)%in%c("Mito")]<-"3"

        #Drosophila mito
        translations$genetic_code[rownames(translations)%in%c("dmel_mitochondrion_genome")]<-"5"

        #the codes given by the user replace those of the names
        if(!is.null(genetic_codes)){translations[names(genetic_codes),"genetic_code"]<-as.character(genetic_codes)}

        circs<-ifs@seqnames[which(ifs@is_circular)]


        #define start and stop codons (genome space)

        #the CDS of the chromosomes that are not circular (of the circular ones if all CDS are on them), each
        #translated with the genetic code of its chromosome
        tocheck<-as.character(runValue(seqnames(cds_tx)))
        lin<-!tocheck%in%circs
        if(!any(lin)){lin<-!lin}
        cd<-translations[tocheck[lin],"genetic_code"]
        tocheck<-cds_tx[lin]
        #a CDS past the end of a chromosome that is not circular stops extractTranscriptSeqs() with "Out of bounds
        #ranges detected": the GTF file is not of the genome sequence
        cds_ranges<-unlist(tocheck)
        cds_chrs<-as.character(seqnames(cds_ranges))
        past_end<-which(end(cds_ranges)>seqlengths(ifs)[cds_chrs] & !cds_chrs%in%circs)
        if(length(past_end)>0){
            chr<-cds_chrs[past_end[1]]
            stop("The GTF file ",gtf_file," is not of the genome sequence ",genome$path,": the CDS of ",length(unique(names(cds_ranges)[past_end])),
                 " transcripts end after the end of their chromosome, e.g. ",names(cds_ranges)[past_end[1]]," at ",end(cds_ranges)[past_end[1]],
                 " on ",chr,", which has ",seqlengths(ifs)[[chr]]," nt in the genome sequence. The GTF file and the FASTA file must be of the same genome assembly",
                 call. = FALSE)
        }
        seqcds<-extractTranscriptSeqs(genome,transcripts = tocheck)
        trslend<-character(length(seqcds))
        for(code in unique(cd)){
            trsl<-suppressWarnings(translate(seqcds[cd==code],genetic.code = getGeneticCode(code),if.fuzzy.codon = "solve"))
            trslend[cd==code]<-as.character(narrow(trsl,end = width(trsl),width = 1))
        }
        stop_inannot<-NA
        if(names(sort(table(trslend),decreasing = TRUE)[1])=="*"){stop_inannot<-"*"}

        cds_txscoords$gene_id<-trann$gene_id[match(as.vector(seqnames(cds_txscoords)),trann$transcript_id)]
        cds_cc<-cds_txscoords
        strand(cds_cc)<-"*"
        sta_cc<-resize(cds_cc,width = 1,"start")
        sta_cc<-unlist(pmapFromTranscripts(sta_cc,exons_tx[seqnames(sta_cc)],ignore.strand=FALSE))
        sta_cc$gene_id<-trann$gene_id[match(names(sta_cc),trann$transcript_id)]
        sta_cc<-sta_cc[sta_cc$hit]
        strand(sta_cc)<-structure(as.vector(strand(transcripts_db)),names=transcripts_db$tx_name)[names(sta_cc)]
        sta_cc$type<-"start_codon"
        mcols(sta_cc)<-mcols(sta_cc)[,c("exon_rank","type","gene_id")]

        sto_cc<-resize(cds_cc,width = 1,"end")
        #stop codon is the 1st nt, e.g. U of the UAA
        #To-do: update with regards to different organelles, and different annotations
        sto_cc<-shift(sto_cc,-2)

        if(is.na(stop_inannot)){sto_cc<-resize(trim(shift(sto_cc,3)),width = 1,fix = "end")}

        sto_cc<-unlist(pmapFromTranscripts(sto_cc,exons_tx[seqnames(sto_cc)],ignore.strand=FALSE))
        sto_cc<-sto_cc[sto_cc$hit]
        sto_cc$gene_id<-trann$gene_id[match(names(sto_cc),trann$transcript_id)]
        strand(sto_cc)<-structure(as.vector(strand(transcripts_db)),names=transcripts_db$tx_name)[names(sto_cc)]
        sto_cc$type<-"stop_codon"
        mcols(sto_cc)<-mcols(sto_cc)[,c("exon_rank","type","gene_id")]


        #define most common, most upstream/downstream

        message("Defining most common start/stop codons ... ",date())

        start_stop_cc<-sort(c(sta_cc,sto_cc))
        start_stop_cc$transcript_id<-names(start_stop_cc)
        start_stop_cc$most_up_downstream<-FALSE
        start_stop_cc$most_frequent<-FALSE

        df<-cbind.DataFrame(start(start_stop_cc),start_stop_cc$type,start_stop_cc$gene_id)
        colnames(df)<-c("start_pos","type","gene_id")
        upst<-by(df$start_pos,INDICES = df$gene_id,function(x){x==min(x) | x==max(x)})
        start_stop_cc$most_up_downstream[unlist(split(seq_len(nrow(df)),df$gene_id))]<-unlist(upst)

        mostfr<-by(df[,c("start_pos","type")],INDICES = df$gene_id,function(x){
            mfreq<-table(x)
            x$start_pos%in%as.numeric(names(which(mfreq[,1]==max(mfreq[,1])))) | x$start_pos%in%as.numeric(names(which(mfreq[,2]==max(mfreq[,2]))))
        })

        start_stop_cc$most_frequent[unlist(split(seq_len(nrow(df)),df$gene_id))]<-unlist(mostfr)

        names(start_stop_cc)<-NULL



        #define transcripts as containing frequent start/stop codons or most upstream ones, in relation with 5'UTR length

        mostupstr_tx<-sum(LogicalList(split(start_stop_cc$most_up_downstream,start_stop_cc$transcript_id)))[as.character(seqnames(cds_txscoords))]
        cds_txscoords$upstr_stasto<-mostupstr_tx
        mostfreq_tx<-sum(LogicalList(split(start_stop_cc$most_frequent,start_stop_cc$transcript_id)))[as.character(seqnames(cds_txscoords))]
        cds_txscoords$mostfreq_stasto<-mostfreq_tx
        cds_txscoords$lentx<-sum(width(exons_tx[as.character(seqnames(cds_txscoords))]))
        df<-cbind.DataFrame(as.character(seqnames(cds_txscoords)),width(cds_txscoords),start(cds_txscoords),cds_txscoords$mostfreq_stasto,cds_txscoords$gene_id)
        colnames(df)<-c("txid","cdslen","utr5len","var","gene_id")
        repres_freq<-by(df[,c("txid","cdslen","utr5len","var")],df$gene_id,function(x){
            x<-x[order(x$var,x$utr5len,x$cdslen,decreasing = TRUE),]
            x<-x[x$var==max(x$var),]
            ok<-x$txid[which(x$cdslen==max(x$cdslen) & x$utr5len==max(x$utr5len) & x$var==max(x$var))][1]
            if(length(ok)==0 | is.na(ok[1])){ok<-x$txid[1]}
            ok
        })

        df<-cbind.DataFrame(as.character(seqnames(cds_txscoords)),width(cds_txscoords),start(cds_txscoords),cds_txscoords$upstr_stasto,cds_txscoords$gene_id)
        colnames(df)<-c("txid","cdslen","utr5len","var","gene_id")
        repres_upstr<-by(df[,c("txid","cdslen","utr5len","var")],df$gene_id,function(x){
            x<-x[order(x$var,x$utr5len,x$utr5len,decreasing = TRUE),]
            x<-x[x$var==max(x$var),]
            ok<-x$txid[which(x$cdslen==max(x$cdslen) & x$utr5len==max(x$utr5len) & x$var==max(x$var))][1]
            if(length(ok)==0 | is.na(ok[1])){ok<-x$txid[1]}
            ok
        })
        df<-cbind.DataFrame(as.character(seqnames(cds_txscoords)),width(cds_txscoords),start(cds_txscoords),cds_txscoords$upstr_stasto,cds_txscoords$gene_id)
        colnames(df)<-c("txid","cdslen","utr5len","var","gene_id")
        repres_len5<-by(df[,c("txid","cdslen","utr5len","var")],df$gene_id,function(x){
            x<-x[order(x$utr5len,x$var,x$cdslen,decreasing = TRUE),]
            ok<-x$txid[which(x$utr5len==max(x$utr5len) & x$var==max(x$var))][1]
            if(length(ok)==0 | is.na(ok[1])){ok<-x$txid[1]}
            ok
        })

        cds_txscoords$reprentative_mostcommon<-as.character(seqnames(cds_txscoords))%in%unlist(repres_freq)
        cds_txscoords$reprentative_boundaries<-as.character(seqnames(cds_txscoords))%in%unlist(repres_upstr)
        cds_txscoords$reprentative_5len<-as.character(seqnames(cds_txscoords))%in%unlist(repres_len5)
        unq_stst<-start_stop_cc
        mcols(unq_stst)<-NULL
        unq_stst<-sort(unique(unq_stst))
        ov<-findOverlaps(unq_stst,start_stop_cc,type="equal")
        ov<-split(subjectHits(ov),queryHits(ov))
        unq_stst$type<-CharacterList(lapply(ov,FUN = function(x){unique(start_stop_cc$type[x])}))
        unq_stst$transcript_id<-CharacterList(lapply(ov,FUN = function(x){start_stop_cc$transcript_id[x]}))
        unq_stst$gene_id<-CharacterList(lapply(ov,FUN = function(x){unique(start_stop_cc$gene_id[x])}))

        unq_stst$reprentative_mostcommon<-sum(!is.na(match(unq_stst$transcript_id,unlist(as(repres_freq,"CharacterList")))))>0
        unq_stst$reprentative_boundaries<-sum(!is.na(match(unq_stst$transcript_id,unlist(as(repres_upstr,"CharacterList")))))>0
        unq_stst$reprentative_5len<-sum(!is.na(match(unq_stst$transcript_id,unlist(as(repres_len5,"CharacterList")))))>0


        #put in a list
        GTF_annotation<-list(transcripts_db,txs_gene,ifs,unq_stst,cds_tx,intron_names_tx,cds_gen,exons_tx,nsns,unq_intr,genes,threeutrs,fiveutrs,ncisof,ncrnas,introns,intergenicRegions,trann,cds_txscoords,translations,pkgnm,stop_inannot,genome)
        names(GTF_annotation)<-c("txs","txs_gene","seqinfo","start_stop_codons","cds_txs","introns_txs","cds_genes","exons_txs","exons_bins","junctions","genes","threeutrs","fiveutrs","ncIsof","ncRNAs","introns","intergenicRegions","trann","cds_txs_coords","genetic_codes","genome_package","stop_in_gtf","genome")

        txs_all<-unique(GTF_annotation$trann$transcript_id)
        txs_exss<-unique(names(GTF_annotation$exons_txs))

        txs_notok<-txs_all[!txs_all%in%txs_exss]
        if(length(txs_notok)>0){
            examples<-withr::with_seed(666,txs_notok[sample(seq_along(txs_notok),size = min(3,length(txs_notok)),replace = FALSE)])
            message(paste("Warning: ",length(txs_notok)," txs with incorrect/unspecified exon boundaries - e.g. trans-splicing events, examples: "
                          ,paste(examples,collapse=", ")," - ",date(),sep = ""))
        }

        #Save as a RData object
        save(GTF_annotation,file=annot_file)
        message("Rannot object created!   ",date())

        #create tables and bed files (with colnames, so with header)
        if(export_bed_tables_TxDb==TRUE){
            message("Exporting annotation tables ... ",date())
            for(bed_file in c("fiveutrs","threeutrs","ncIsof","ncRNAs","introns","cds_txs_coords")){
                bf<-GTF_annotation[[bed_file]]
                bf_t<-data.frame(chromosome=seqnames(bf),start=start(bf),end=end(bf),name=rep(".",length(bf)),score=width(bf),strand=strand(bf))
                meccole<-mcols(bf)
                for(mecc in names(meccole)){
                    if(is(meccole[,mecc],"CharacterList") | is(meccole[,mecc],"NumericList") | is(meccole[,mecc],"IntegerList")){
                        meccole[,mecc]<-paste(meccole[,mecc],collapse=";")
                    }
                }
                bf_t<-cbind.data.frame(bf_t,meccole)
                write.table(bf_t,file = paste(annotation_directory,"/",bed_file,"_similbed.bed",sep=""),sep="\t",quote = FALSE,row.names = FALSE)

            }

            write.table(GTF_annotation$trann,file = paste(annotation_directory,"/table_gene_tx_IDs",sep=""),sep="\t",quote = FALSE,row.names = FALSE)
            seqi<-as.data.frame(GTF_annotation$seqinfo)
            seqi$chromosome<-rownames(seqi)
            write.table(seqi,file = paste(annotation_directory,"/seqinfo",sep=""),sep="\t",quote = FALSE,row.names = FALSE)

            gen_cod<-as.data.frame(GTF_annotation$genetic_codes)
            gen_cod$chromosome<-rownames(gen_cod)
            write.table(gen_cod,file = paste(annotation_directory,"/genetic_codes",sep=""),sep="\t",quote = FALSE,row.names = FALSE)
            message("Exporting annotation tables --- Done! ",date())

        }

    }

    return(annot_file)
}
