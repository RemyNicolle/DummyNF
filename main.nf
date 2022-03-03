

// ~/nextflow  -c nextflowGC.config run . -work-dir gs://pipebuff/test3p -resume


// ~/nextflow  -c nextflow.config run . -resume


Channel.fromList(file(params.sampleList).readLines())
.map {
  [it ,  file(params.sampleInputDir + "/" + it + params.samPsuffix )] }
.set { samples_ch }


/*
* Step 1. Builds the genome index required by the mapping process
*/

process buildIndex {
  // publishDir  "${params.outputDir}", mode: 'copy'

  label 'big'

  input:
  val fa from params.refFasta
  val gtf from params.refGTF
  val info from params.ensemblRefID

  output:
  path "StarIndexed" into index_ch


  """
  wget -nv -O /tmp/genome.fa.gz $fa
  wget -nv -O /tmp/trascriptome.gtf.gz $gtf


  gunzip /tmp/genome.fa.gz
  gunzip /tmp/trascriptome.gtf.gz



  STAR --runThreadN $task.cpus --runMode genomeGenerate --genomeDir StarIndexed \
  --genomeFastaFiles /tmp/genome.fa  --sjdbOverhang 100 --sjdbGTFfile /tmp/trascriptome.gtf

  mv /tmp/trascriptome.gtf StarIndexed/ref.gtf

  echo -e "refGT\t"$info >StarIndexed/description.txt


  grep -P "\tgene\t" StarIndexed/ref.gtf >ref.GeneLvlOnly.gtf


  R -e 'library(parallel) ;gtf=read.delim("ref.GeneLvlOnly.gtf",sep="\\t",as.is=T,header=F) ; gtf.getmeta=function(apieceofgtf){strsplit((apieceofgtf)[,9],";| ")} ; gtf.getmetavalue=function(apieceofgtf,field){ metadata=gtf.getmeta(apieceofgtf)   ;    unlist(mclapply(metadata,function(x){if(field%in% x){return(x[which(x==field)+1])}else{return(NA)} }))}  ; geneTab =gtf[which(gtf[,3] =="gene"),];geneTab[["GeneID"]]=gtf.getmetavalue(geneTab,"gene_id")  ; geneTab[["GeneName"]]=gtf.getmetavalue(geneTab,"gene_name");geneTab=unique(geneTab)  ; rownames(geneTab)=geneTab[["GeneID"]];geneTab=unique(geneTab[,- c(9,6,2,3,8)])  ; saveRDS(geneTab,file="StarIndexed/refGeneID.rds")'


  """

  // indexSTAR $gtf $fa $params.outputDir $task.cpus $info
}



/*
* Step 2. Maps each read
*/
process doSTAR {
  tag "STAR $sample"
  label 'big'


  input:
  path index from index_ch
  // set sample from samples_ch
  tuple val(sample), file(fq) from samples_ch


  // val sample from samples_ch
  // path inputd from params.sampleInputDir
  // val suffix from params.samPsuffix

  output:
  // tuple val(sample), file("${sample}StarOutAligned.sortedByCoord.out.bam") ,file("${sample}CountSummary.txt" ) into aligned_ch
  file "${sample}StarOutAligned.sortedByCoord.out.bam" into bam_ch
  file "${sample}CountSummary.txt" into starlog_ch

 // gunzip -c ${inputd}/${sample}${suffix} > /tmp/this.fq
 // echo ${sample}
 // echo ${fq}

  """
  gunzip -c ${fq} > /tmp/this.fq

  STAR --genomeDir $index \
  --readFilesIn /tmp/this.fq \
  --outFileNamePrefix $sample'StarOut' \
  --runThreadN $task.cpus --sjdbGTFfile $index/ref.gtf \
  --twopassMode None --outFilterType BySJout  --seedSearchStartLmax 12 \
  --alignEndsType Local --outSAMtype BAM SortedByCoordinate \
  --alignIntronMax 1000000 \
  --alignMatesGapMax 1000000 \
  --limitOutSJcollapsed 1000000 \
  --limitSjdbInsertNsj 1000000 \
  --outFilterMultimapNmax 100 --winAnchorMultimapNmax 50 \
  --alignSJoverhangMin 15 \
  --alignSJDBoverhangMin 1 \
  --alignIntronMin 20 \
  --outFilterMatchNminOverLread 0 \
  --outFilterScoreMinOverLread   ${params.filterScoreMinOverLread} \
  --outFilterMismatchNmax ${params.filterMismatchNmax} \
  --outFilterMismatchNoverLmax ${params.filterMismatchNoverLmax}

  sed "s/ | */;/" ${sample}StarOutLog.final.out | sed "s/^  *//" | grep -v ":\$" | \
  sed "s/% ?//g"|sed "s/^%/Prop/" | sed "s/ /_/g" | tr -cd "[:alnum:]._;[:space:]" | \
  grep -v "^Start"|grep -v "^Mapp"|grep -v "^Fin" | grep -v "^\$" | tr -d ";" > "${sample}CountSummary.txt"
  """




  // countLexoFWD "/tmp/StarOut" $sample $index $task.cpus
}

/*
* Step 2. Maps each read
*/
process CountNagreg {
  tag "FC "
  label 'small'

  publishDir "${params.outputDir}", mode: 'copy'


  input:
  path index from index_ch
  file allbams from bam_ch.collect()
  file starlogs from starlog_ch.collect()
  val info from params.ensemblRefID

  output:
  file "FCcounts.tsv.gz" into countout
  file "RNAseqSummary.tsv" into seqsumout
  file "referenceGeneAnnotation_${info}.rds" into annotout

  // subaaCountSummary.txt subabCountSummary.txt subacCountSummary.txt

  """

  featureCounts -O -F GTF -T $task.cpus -t exon -g gene_id -s 1 -a $index/ref.gtf \
  -o FCcounts_all $allbams
  echo $starlogs
  cat $starlogs >tmpstarlogs
  sed "s/StarOutAligned.sortedByCoord.out.bam//g" FCcounts_all >fctmp

  R -e 'x=read.delim("fctmp",header=T,as.is=T,comment.char="#");rownames(x)=x[,"Geneid"]; write.table(x[, -(1:6) ],file="FCcounts.tsv",sep="\\t",quote=F);starsum=do.call(rbind,lapply(strsplit("$starlogs"," ")[[1]],function(y){s=gsub("CountSummary.txt","",y);ysum=data.frame(t(read.delim(y,header=F)[,2]),row.names=s);colnames(ysum)=read.delim(y,header=F)[,1];ysum}));fcsum=t(read.delim("FCcounts_all.summary",as.is=T,header=T)[c(1,9,12),-1]);rownames(fcsum)=gsub("StarOutAligned.sortedByCoord.out.bam","",rownames(fcsum));colnames(fcsum)=c("Assigned","Unassigned_MultiMapping","Unassigned_NoFeatures");write.table(data.frame(starsum,fcsum[rownames(starsum),],nGenes=apply(x[,rownames(starsum)]>0,2,sum)),file="RNAseqSummary.tsv",sep="\\t",quote=F)'


  gzip FCcounts.tsv

  cp $index/refGeneID.rds referenceGeneAnnotation_${info}.rds


  """

}


/*
* completion handler
*/
workflow.onComplete {
  log.info ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}

// workflow.onComplete {
// 	log.info ( workflow.success ? "\nDone! Open the following report in your browser --> $params.outdir/multiqc_report.html\n" : "Oops .. something went wrong" )
// }
