# kmeans() picks random start centres. With a fixed seed and fixed kinds of
# random numbers, the same data give the same clusters, and so the same offsets,
# in each run. with_seed() restores the random numbers of the session after the call.
kmeans_seeded<-function(...){
    withr::with_seed(666,kmeans(...),.rng_kind="Mersenne-Twister",.rng_normal_kind="Inversion",.rng_sample_kind="Rejection")
}


#' Filter read lengths for P-sites position calculation
#'
#' This function selects a subset of readlenghts to be used in the P-sites calculation step
#' @keywords RiboseQC
#' @author Lorenzo Calviello, \email{calviello.l.bio@@gmail.com}
#' @param summary_data for each compartment, a matrix with the \code{cutoff} and the \code{frame_preference} of each read length (rows),
#' from the output of the \code{calc_cutoffs_from_profiles} function: \code{res_all$selection_cutoffs$results_cutoffs}
#' @param choice Method used to select readlengths, defaults to "max_coverage"
#' @param nt_signals Profiles of 5'ends of each transcript, for the compartments of \code{summary_data}, in the same order:
#' \code{res_all$profiles_fivepr$five_prime_subcodon}. Only the windows of the CDS are used
#' @details Three different methods are available to choose readlengths: the "max_coverage" method
#' selects all read lenghts with more in-frame signal compared to out-of-frame signal, on the codons of the windows of the CDS
#' (the positions with signal after the cutoff, on average over the transcripts);
#' the "max_inframe" method starts with the most accurate read length and
#' progressively selects read lengths which add in-frame signals in codons not covered by previous read lengths;
#' the "all" method selects all available read lengths
#' @return a list object containing different compartments (those with read lengths with a cutoff). Each sub-list contains \code{final_choice},
#' the set of chosen read lengths with cutoffs, and \code{data}, the complete stats for each selection method
#' @seealso \code{\link{calc_cutoffs_from_profiles}}
#' @examples
#' \dontrun{
#' # The 5' profiles of each transcript are only in the *_results_RiboseQC_all file
#' # (RiboseQC_analysis() with write_tmp_files = TRUE).
#' load("my_sample_results_RiboseQC_all")
#' choose_readlengths(res_all$selection_cutoffs$results_cutoffs, choice = "max_coverage",
#'   nt_signals = res_all$profiles_fivepr[["five_prime_subcodon"]])
#' }
#' @export

choose_readlengths<-function(summary_data,choice="max_coverage",nt_signals){

    if(!choice%in%c("max_coverage","max_inframe","all")){
        stop("choice must be one of 'max_coverage','max_inframe' or 'all'")
    }
    stopifnot(names(summary_data) == names(nt_signals))
    choices_res<-list()
    for(i in names(summary_data)){
        choices_res[[i]]<-c()
        data<-DataFrame(summary_data[[i]])
        sign<-nt_signals[[i]]
        data<-data[grep(rownames(data),pattern = "all",invert = TRUE),]
        if(dim(data)[1]==0){next}
        data<-data[order(data$frame_preference,decreasing = TRUE),]
        mat_fr<-matrix(0,nrow = nrow(sign[[1]]),ncol = 99)
        net_all<-c()
        net_all_cov<-c()
        gains_cod_all<-c()
        gains_sign_all<-c()
        for(cnt in seq_along(rownames(data))){
            rl<-rownames(data)[cnt]
            cutf<-data$cutoff[cnt]
            sig<-as.matrix(sign[[rl]])
            cdss<-grep(colnames(sig),pattern = "CDS")
            sig_ctf<-cbind(matrix(0,ncol = cutf,nrow = nrow(sig)),sig[,seq_len(length(sig[1,])-cutf)])[,cdss]
            sig_ctf<-as.matrix(sig_ctf)
            sig_ctf_masked<-sig_ctf
            sig_ctf_masked[mat_fr>0]<-NA

            codon_cov<-t(apply(as.matrix(sig_ctf),1,function(y){
                fr1<-y[seq(1,99,by=3)]
                fr2<-y[seq(2,99,by=3)]
                fr3<-y[seq(3,99,by=3)]
                fr1<-fr1[!is.na(fr1)]
                fr2<-fr2[!is.na(fr2)]
                fr3<-fr3[!is.na(fr3)]
                codongain<-c(sum(fr1>0),sum(fr2>0),sum(fr3>0))
                codongain<-codongain/sum(codongain)
                codongain

            }))
            codon_cov<-codon_cov[!apply(codon_cov,1,FUN = function(y){sum(is.na(y))>0}),]


            codon_gain<-t(apply(sig_ctf_masked,1,function(y){
                fr1<-y[seq(1,99,by=3)]
                fr2<-y[seq(2,99,by=3)]
                fr3<-y[seq(3,99,by=3)]
                fr1<-fr1[!is.na(fr1)]
                fr2<-fr2[!is.na(fr2)]
                fr3<-fr3[!is.na(fr3)]
                codongain<-c(sum(fr1>0),sum(fr2>0),sum(fr3>0))
                codongain<-codongain/sum(codongain)
                codongain

            }))
            codon_gain<-codon_gain[!apply(codon_gain,1,FUN = function(y){sum(is.na(y))>0}),]


            av_codoncov<-matrix(colMeans(codon_cov),ncol=3)
            av_codongain<-matrix(colMeans(codon_gain),ncol=3)
            colnames(av_codongain)<-c("gain_in_frame","gain_off_1","gain_off_2")
            rownames(av_codongain)<-rl
            colnames(av_codoncov)<-c("in_frame","off_1","off_2")
            rownames(av_codoncov)<-rl



            gain_cod<-av_codongain[1]-(av_codongain[2]+av_codongain[3])
            gain_sign<-av_codoncov[1]-(av_codoncov[2]+av_codoncov[3])
            gains_cod_all<-c(gains_cod_all,gain_cod)
            gains_sign_all<-c(gains_sign_all,gain_sign)

            if(gain_cod<=0 & gain_sign<=0){next}
            if(gain_cod>0){
                mat_fr[sig_ctf>1]<-1
                net_all<-rbind(net_all,av_codongain)

            }

            if(gain_sign>0){
                net_all_cov<-rbind(net_all_cov,av_codoncov)

            }


        }

        data$read_length<-rownames(data)
        data$gain_codons<-gains_sign_all
        data$gain_new_codons<-gains_cod_all


        data$all<-TRUE
        data$max_inframe<-FALSE
        if(length(net_all)>0){
            mtc<-match(rownames(net_all),rownames(data))
            data$max_inframe[mtc]<-TRUE
        }

        data$max_coverage<-FALSE
        if(length(net_all_cov)>0){
            mtc<-match(rownames(net_all_cov),rownames(data))
            data$max_coverage[mtc]<-TRUE
        }

        final_choice<-data[,c("read_length","cutoff")][data[,choice],]
        choices<-list(final_choice=final_choice,data=data)
        choices_res[[i]]<-choices
    }
    choices_res
}


#' Offset spliced reads on plus strand
#'
#' This function calculates P-sites positions for spliced reads on the plus strand
#' @keywords RiboseQC
#' @author Lorenzo Calviello, \email{calviello.l.bio@@gmail.com}
#' @param x a \code{GAlignments} object with a cigar string
#' @param cutoff number representing the offset value
#' @return a \code{GRanges} object with offset reads: the P-site (1 nt) of each read, in the order of \code{x}
#' @export

get_ps_fromspliceplus<-function(x,cutoff){
    rang<-cigarillo::cigars_as_ranges_along_ref(cigar(x), lmmpos=start(x),ops="M")
    cs<-lapply(rang,function(x){cumsum(x@width)})
    rangok<-lapply(which(IntegerList(cs)>cutoff),"[[",1)
    rangok<-unlist(rangok)
    gr<-as(x,"GRanges")
    #do the splitaslist thingy
    ones<-rangok==1
    mores<-rangok>1
    psmores<-GRanges()
    psones<-GRanges()

    if(sum(ones)>0){
        psones<-shift(resize(gr[ones],width=1,fix="start"),shift=cutoff)
    }

    if(sum(mores)>0){
        rangmore<-rang[mores]
        rangok<-rangok[mores]
        cms<-cumsum(width(rangmore))
        shft<-c()
        for(i in seq_along(rangok)){shft<-c(shft,cutoff-cms[[i]][rangok[i]-1])}
        stt<-start(rangmore)
        stok<-c()
        for(i in seq_along(shft)){stok<-c(stok,stt[[i]][rangok[i]]+shft[i])}
        psmores<-GRanges(IRanges(start=stok,width = 1),seqnames = seqnames(x[mores]),strand=strand(x[mores]),seqlengths=seqlengths(x[mores]))

    }
    # in the order of x: psites_of_chunk() gives each P-site the mcols of its read
    ps<-c(psones,psmores)[order(c(which(ones),which(mores)))]
    return(ps)
}

#' Offset spliced reads on minus strand
#'
#' This function calculates P-sites positions for spliced reads on the minus strand
#' @keywords RiboseQC
#' @author Lorenzo Calviello, \email{calviello.l.bio@@gmail.com}
#' @param x a \code{GAlignments} object with a cigar string
#' @param cutoff number representing the offset value
#' @return a \code{GRanges} object with offset reads: the P-site (1 nt) of each read, in the order of \code{x}
#' @export

get_ps_fromsplicemin<-function(x,cutoff){
    rang<-cigarillo::cigars_as_ranges_along_ref(cigar(x), lmmpos=start(x),ops="M")
    rang<-endoapply(rang,rev)
    cs<-lapply(rang,function(x){cumsum(x@width)})
    rangok<-lapply(which(IntegerList(cs)>cutoff),"[[",1)
    rangok<-unlist(rangok)
    gr<-as(x,"GRanges")
    #do the splitaslist thingy
    ones<-rangok==1
    mores<-rangok>1
    psmores<-GRanges()
    psones<-GRanges()


    if(sum(ones)>0){
        psones<-shift(resize(gr[ones],width=1,fix="start"),shift=-cutoff)
    }


    if(sum(mores)>0){
        rangmore<-rang[mores]
        rangok<-rangok[mores]
        cms<-cumsum(width(rangmore))
        shft<-c()
        for(i in seq_along(rangok)){shft<-c(shft,cutoff-cms[[i]][rangok[i]-1])}
        #start?
        stt<-end(rangmore)
        stok<-c()
        for(i in seq_along(shft)){stok<-c(stok,stt[[i]][rangok[i]]-shft[i])}
        psmores<-GRanges(IRanges(start=stok,width = 1),seqnames = seqnames(x[mores]),strand=strand(x[mores]),seqlengths=seqlengths(x[mores]))

    }
    # in the order of x: psites_of_chunk() gives each P-site the mcols of its read
    ps<-c(psones,psmores)[order(c(which(ones),which(mores)))]
    return(ps)
}

#' Calculate offsets from 5'end profiles
#'
#' This function calculates cutoffs and frame resolution for Ribo-seq reads,
#'  for each read length and compartment.
#' @details Four methods are used and combined in the final choice:
#'  the position of maximum coverage from 24 nt before the start codon to its first nt (so the cutoff is at most 24) is calculated for each transcript, and
#'  the most frequent one is stored in the "*_tab" objects. Such frequency values are also
#'  subjected to k-means clustering (centers=3) and the first value belonging to the highest cluster
#'  is selected, output as "*km_tab" objects. Analysis of aggregate plots, instead of frequencies, is performed
#'  again using kmeans (centers=3) using the same analysis above and stored in the "*km_meta" objects,
#'  and by simply calculating the maximum value in the profile, stored in the "*meta" objects.
#'  for each method, all positions ("absolute_") or only in-frame positions ("in_frame_") are considered. The final choice takes the most
#'  frequent cutoff chosen in all methods applied to in-frame positions.
#' @param reads_profile Profile of 5'ends around start and stop codon, as a DataFrame object with
#'   tx_ids as rows and positions as columns: the columns of \code{five_prime_subcodon} in \code{\link{RiboseQC_analysis}}
#'   (\code{CDS_1} is the first nt of the start codon)
#' @param length_max Maximum cutoff to use (the cutoff is at most 24 in any case)
#' @keywords RiboseQC
#' @author Lorenzo Calviello, \email{calviello.l.bio@@gmail.com}
#' @return a list with a \code{final_cutoff} object, the frame analysis containing the \% of transcripts
#' displaying the max frame and the average \% of frame precision per transcript in \code{frames_res},
#' all the calculated cutoffs in \code{cutoffs}, data used for the frame analysis in \code{frames},
#' and profiles around start codons in  \code{profiles_start}. \code{NULL} if fewer than 5 transcripts have 5'ends
#' in the CDS, or from 24 nt before the start codon to its first nt.
#' @seealso \code{\link{RiboseQC_analysis}}
#' @examples
#' \dontrun{
#' # The 5' profiles of each transcript are only in the *_results_RiboseQC_all file
#' # (RiboseQC_analysis() with write_tmp_files = TRUE).
#' load("my_sample_results_RiboseQC_all")
#' calc_cutoffs_from_profiles(res_all$profiles_fivepr[["five_prime_subcodon"]][["nucl"]][["28"]],
#'   length_max = 28)
#' }
#' @export

calc_cutoffs_from_profiles<-function(reads_profile,length_max){
    profok<-as.matrix(reads_profile[,grep("CDS",x = colnames(reads_profile))])
    sums<-rowSums(profok)
    profok<-profok[sums>0,,drop=FALSE]
    #what if no reads?
    if((dim(profok)[1])<5){return(c())}

    frames<-t(apply(profok,1,function(y){
        sumy<-sum(y)
        fr0<-sum(y[seq(1,length(y),by=3)])/sumy
        fr1<-sum(y[seq(2,length(y),by=3)])/sumy
        fr2<-sum(y[seq(3,length(y),by=3)])/sumy
        c(fr0,fr1,fr2)
    }))

    tbfra<-table(apply(frames,1,which.max))
    tbfra<-(tbfra/sum(tbfra))*100
    fra<-as.numeric(names(which.max(tbfra)))-1
    frameans<-colMeans(frames)
    frameans<-(frameans/sum(frameans))*100
    names(frameans)<-seq_len(3)

    infra<-seq(fra,-24,by=-3)

    #just offset
    profok<-as.matrix(reads_profile)
    Atg<-which(colnames(profok)=="CDS_1")
    colnames(profok)<-as.character(seq(from=-(Atg-1),length.out = length(profok[1,]),by=1))

    profok<-profok[,colnames(profok)%in%as.character(-24:0)]
    profok<-profok[,abs(as.numeric(colnames(profok)))<=length_max]
    sums<-rowSums(profok)
    profok<-profok[sums>0,,drop=FALSE]
    #what if no reads?
    if((dim(profok)[1])<5){return(c())}


    metaprof<-colSums(profok)
    names(metaprof)<-as.numeric(names(metaprof))
    metaprof_fra<-metaprof[names(metaprof)%in%infra]
    km_meta<-NA
    km_meta_fra<-NA
    if(length(unique(metaprof))>3){
        kms<-kmeans_seeded(metaprof,centers=3)$cluster
        km_meta<-abs(as.numeric(names(metaprof[which(kms==kms[which.max(metaprof)])[1]])))
    }
    if(length(unique(metaprof_fra))>3){
        kms<-kmeans_seeded(metaprof_fra,centers=3)$cluster
        km_meta_fra<-abs(as.numeric(names(metaprof_fra[which(kms==kms[which.max(metaprof_fra)])[1]])))

    }
    tbcut<-table(apply(profok,1,which.max))
    names(tbcut)<-abs(as.numeric(colnames(profok)[as.numeric(names(tbcut))]))
    cutoff<-abs(as.numeric(names(tbcut)[which.max(tbcut)]))
    cutoff_km<-NA
    cutoff_km_fra<-NA
    if(length(cutoff)==0){cutoff=NA}
    tbcut_fra<-tbcut[as.numeric(names(tbcut))%in%(-infra)]
    if(length(tbcut_fra)==0){tbcut_fra=NA}
    cutoff_fra<-abs(as.numeric(names(tbcut_fra)[which.max(tbcut_fra)]))
    if(length(cutoff_fra)==0){cutoff_fra=NA}
    if(length(cutoff)>0 & length(unique(tbcut))>4){
        kms<-kmeans_seeded(tbcut,centers=3)$cluster
        cutoff_km<-abs(as.numeric(names(tbcut[which(kms==kms[which.max(tbcut)])[1]])))
        if(length(unique(tbcut_fra))>3){
            kms<-kmeans_seeded(tbcut_fra,centers=3)$cluster
            cutoff_km_fra<-abs(as.numeric(names(tbcut_fra[which(kms==kms[which.max(tbcut_fra)])[1]])))
        }
    }
    cutoff_meta<-abs(as.numeric(names(metaprof)[which.max(metaprof)]))
    cutoff_meta_fra<-abs(as.numeric(names(metaprof_fra)[which.max(metaprof_fra)]))
    if(length(cutoff_meta)==0){cutoff_meta=NA}
    if(length(cutoff_meta_fra)==0){cutoff_meta_fra=NA}

    cutoffs<-c(cutoff,cutoff_fra,cutoff_km,cutoff_km_fra,km_meta,km_meta_fra,cutoff_meta,cutoff_meta_fra)
    names(cutoffs)<-c("absolute_tab","in_frame_tab","absolute_km_tab","in_frame_km_tab","absolute_km_meta","in_frame_km_meta","absolute_meta","in_frame_meta")
    cutt<-cutoffs[!is.na(cutoffs)]
    final_cutoff<-as.numeric(names(sort(table(cutt[grep(names(cutt),pattern = "frame")]),decreasing = TRUE)[1]))
    frames_res<-c(sort(tbfra,decreasing = TRUE)[1],as.numeric(frameans[which.max(frameans)]))
    names(frames_res)<-c("max_frame_in_pctORF","mean_pct_max_frame")

    results<-list(final_cutoff=final_cutoff,frames_res=frames_res,cutoffs=cutoffs,frames=frames,profiles_start=profok)
    return(results)

}
