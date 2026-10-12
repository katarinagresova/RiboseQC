

testthat::test_that("prepare_annotation_files works as expected",{
  
  arab_fasta <- system.file(package='RiboseQC',"ext_data","simp_arab.fa.gz")
  tmp_dir = tempdir()
  # An uncompressed copy: Rsamtools reads only uncompressed or bgzip FASTA files
  Biostrings::writeXStringSet(Biostrings::readDNAStringSet(arab_fasta),file.path(tmp_dir,'simp_arab.fa'))
  arab_fasta <- file.path(tmp_dir,'simp_arab.fa')
  
  arab_gtf <- system.file(package='RiboseQC',"ext_data","simp_arab.gtf.gz")
 
   
  annotation <- prepare_annotation_files(annotation_directory = tempdir(),
                                         genome_seq = arab_fasta,
                                         gtf_file = arab_gtf,
                                         scientific_name = "arabidopsis.test",
                                         annotation_name = "araport11_custom",export_bed_tables_TxDb = TRUE,
                                         forge_BSgenome = FALSE,
                                         create_TxDb = TRUE)
  
  expect_true(file.exists(annotation))
  
})

