#!/usr/bin/Rscript


# Author: David Fournier - March/June 2019
# a script to extract 1:1 orthologues of mouse to fill putative missing data from biomaRt
# to launch with launch_Rscript.sh 
# or directly from command line with: 
# srun -p short -A jgu-cbdm R --no-save --slave -f ./blast.R &

# LIBRARIES
print(paste0("time started: ", Sys.time()))

#library("biomaRt")  # not used
library("foreach")
library("doParallel")
#library("data.table")  # not used
 library("stringr")

# GLOBAL VARIABLES

#wdir="/project/jgu-cbdm/andradeLab/scratch/dfournie/prot_evol/orthologues_finding/blastp" # location of blast software suite
wdir="/home/dfournie/dfournie/prot_evol/orthologues_finding/blastp" # location of blast software suite
dbdir=wdir

#args=commandArgs(TRUE)
#dbdir=args[1] # location of blast db on the ramdisk
args=c("",1)
#inputf=args[3] # file chunk to process



# FUNCTIONS 

get_match_id <- function(blast_output_file) {
tmp=readLines(blast_output_file)
rec.id.tmp="" # variable initiation
if(length(tmp)>0){ # retrieving Ensembl gene sequence corres. to the protein in blastp output
tmp2=str_match(tmp, "gene:([A-Z0-9]+)") # pattern b. "gene:" and a dot (Ensembl gene has a version with it, like gene:ENSMUSG00000020681.14) or a space (if gene with no version like gene:ENSMUSG00000020681)
if(length(tmp2)>1){
rec.id.tmp=tmp2[,2] 
}
}
return(rec.id.tmp)
}

ortho.blast.CDS <- function(mouse_gene_id,mouse_prot_id, species_id) {  # script to find 1:1 orthologue sequence by doing a reciprocal blast on Ensembl CDS protein database

 seq_name="no_ortho" # in case no orthologue present even on transcriptome sequencing data
 seq_id="" # intermediary variable
# step1: tblastn against against Ensembl CDS database
cmd=paste0("system(\"",wdir, "/ncbi-blast-2.8.1+/bin/tblastn -db ",dbdir,"/db/nt/",species_id,"_CDS -query ",dbdir,"/",mouse_prot_id,".fasta -out ",dbdir,"/tmp_",mouse_prot_id,".txt -num_threads 2 -outfmt '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq' -max_target_seqs 1\")") # tblastn to query target species CDS from Ensembl    --- qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore
eval(parse(text=cmd))

if(length(readLines(paste0(dbdir,"/tmp_",mouse_prot_id,".txt")))>0){ # homologue found
# step2: blastp of the result nucleotide sequence against mouse protein db
tmp=read.table(paste0(dbdir,"/tmp_",mouse_prot_id,".txt"),h=F,stringsAsFactors = FALSE) # score and evalue don't matter here - could be low
tar_seq=gsub('-', '', tmp[1,13])
seq_id=tmp[1,2]  # for final output table
cmd=paste0("system(\"echo '>target_seq\n",tar_seq,"' > ",dbdir,"/target_",mouse_prot_id,".fasta\")")
eval(parse(text=cmd))
cmd=paste0("system(\"",wdir,"/ncbi-blast-2.8.1+/bin/blastp -db ",dbdir,"/db/nr/Mmusculus_Ensembl -query ",dbdir,"/target_",mouse_prot_id,".fasta -out ",dbdir,"/tmp_",mouse_prot_id,".txt -num_threads 2 -outfmt '6 stitle' -max_target_seqs 1\")") # blastp possible as we get the predicted protein sequence
eval(parse(text=cmd))

rec.id=get_match_id(paste0(dbdir,"/tmp_",mouse_prot_id,".txt"))
if( rec.id==mouse_gene_id ){ seq_name=seq_id }  # 1:1 orthology detected

}
  return(seq_name)  # value "no_ortho" by default
}

ortho.blast <- function(mouse_gene_id, mouse_prot_id,species_id) {  # script to find 1:1 orthologue sequence by doing a reciprocal blast
# gene id needed when doing recriprocal blast to avoid ambiguity b. multiple protein ids
presence=0; seq_name="no_ortho" # values by default

if(species_id %in% c("Dnovemcinctus","Mdomestica","Ggallus","Xtropicalis","Drerio","Cmilii","Cintestinalis","Bfloridae","Dmelanogaster","Celegans")){ # choice 1
#seq_id="" # for final output table - empty by default
# choice 1 step1: blastp for all but olatipes, pmarinus and eburgeri
cmd=paste0("system(\"",wdir, "/ncbi-blast-2.8.1+/bin/blastp -db ",dbdir,"/db/nr/",species_id," -query ",dbdir,"/",mouse_prot_id,".fasta -out ",dbdir,"/tmp_",mouse_prot_id,".txt -num_threads 2 -outfmt '6 sgi sseq' -max_target_seqs 1\")")  # parameters left aside: qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore
eval(parse(text=cmd))

if(length(readLines(paste0(dbdir,"/tmp_",mouse_prot_id,".txt")))>0){ # homologue found
# choice 1 step2: blastp back to mouse proteome
tmp=read.table(paste0(dbdir,"/tmp_",mouse_prot_id,".txt"),h=F,stringsAsFactors = FALSE) # score and evalue don't matter here - could be low
tar_seq=gsub('-', '', tmp[1,2])
seq_id=tmp[1,1]  # for final output table
cmd=paste0("system(\"echo '>target_seq\n",tar_seq,"' > ",dbdir,"/target_",mouse_prot_id,".fasta\")")
eval(parse(text=cmd))
cmd=paste0("system(\"",wdir,"/ncbi-blast-2.8.1+/bin/blastp -db ",dbdir,"/db/nr/Mmusculus_Ensembl -query ",dbdir,"/target_",mouse_prot_id,".fasta -out ",dbdir,"/tmp_",mouse_prot_id,".txt -num_threads 2 -outfmt '6  stitle' -max_target_seqs 1\")") # blastp againts mouse Ensembl protein sequences 
eval(parse(text=cmd))

rec.id=get_match_id(paste0(dbdir,"/tmp_",mouse_prot_id,".txt"))
if( !(species_id %in% c("Cmilii","Bfloridae")) && ((length(tmp)==0) || (rec.id!=mouse_gene_id)) ){ 
# case where no orthologue in transcriptome data, we try CDS predicted from Ensembl CDS data 
seq_id=ortho.blast.CDS(mouse_gene_id,mouse_prot_id,species_id)
if(seq_id!="no_ortho"){ rec.id=mouse_gene_id } # orthologue was found by ortho.blast.CDS
}
if( rec.id==mouse_gene_id ){ presence=3; seq_name=seq_id }  # 1:1 orthology detected

} 
}else{ # Olatipes, Pmarinus & Eburgeri
# choice 2 step1: tblastn for olatipes, pmarinus and eburgeri against Ensembl cds database
cmd=paste0("system(\"",wdir, "/ncbi-blast-2.8.1+/bin/tblastn -db ",dbdir,"/db/nt/",species_id,"_CDS -query ",dbdir,"/",mouse_prot_id,".fasta -out ",dbdir,"/tmp_",mouse_prot_id,".txt -num_threads 2 -outfmt '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq' -max_target_seqs 1\")") # tblastn to query target species transcriptome (downloaded from NCBI)
eval(parse(text=cmd))  # OPTIMIZATION STEP: ONLY GET stitle and SEQUENCE AND THEN LOAD IT ONCE WITH readLines - no read.table and gsub after - also just display 1 line to go faster 

if(FALSE){  # beta-test of new solution if no orthologue found
cmd=paste0("system(\"",wdir, "/ncbi-blast-2.8.1+/bin/tblastn -db ",dbdir,"/db2/nt/",species_id,"_WGS -query ",dbdir,"/",mouse_prot_id,".fasta -out ",dbdir,"/tmp_",mouse_prot_id,".txt -num_threads 2 -outfmt '6 sqaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq' -max_target_seqs 100\")") # tblastn to query target species transcriptome (downloaded from NCBI)  # back-up: accver sstart send evalue bitscore sseq
eval(parse(text=cmd))  # OPTIMIZATION STEP: ONLY GET stitle and SEQUENCE AND THEN LOAD IT ONCE WITH readLines - no read.table and gsub after - also just display 1 line to go faster 
if(length(readLines(paste0(dbdir,"/tmp_",mouse_prot_id,".txt")))>0){ # homologue found
# choice 2 step2: blastp of the result nucleotide sequence against mouse protein db
tmp=read.table(paste0(dbdir,"/tmp_",mouse_prot_id,".txt"),h=F,stringsAsFactors = FALSE) # score and evalue don't matter here - could be low
tar_seq=gsub('-', '', tmp[1,12])
seq_id=paste0(tmp[1,1],"_",tmp[1,2],"_",tmp[1,3]) # for final output table
cmd=paste0("system(\"echo '>target_seq\n",tar_seq,"' > ",dbdir,"/target_",mouse_prot_id,".fasta\")")
eval(parse(text=cmd))
cmd=paste0("system(\"",wdir,"/ncbi-blast-2.8.1+/bin/blastp -db ",dbdir,"/db/nr/Mmusculus_Ensembl -query ",dbdir,"/target_",mouse_prot_id,".fasta -out ",dbdir,"/tmp_",mouse_prot_id,".txt -num_threads 2 -outfmt '6 stitle' -max_target_seqs 1\")") # IMPORTANT !!! new version -- information left out for recollection: sframe sgi sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq
eval(parse(text=cmd))
rec.id=get_match_id(paste0(dbdir,"/tmp_",mouse_prot_id,".txt"))
if( rec.id==mouse_gene_id ){ presence=3; seq_name=seq_id }  # 1:1 orthology detected
}
} # end of if(FALSE)

if(length(readLines(paste0(dbdir,"/tmp_",mouse_prot_id,".txt")))>0){ # homologue found
# choice 2 step2: blastp of the result nucleotide sequence against mouse protein db
tmp=read.table(paste0(dbdir,"/tmp_",mouse_prot_id,".txt"),h=F,stringsAsFactors = FALSE) # score and evalue don't matter here - could be low
tar_seq=gsub('-', '', tmp[1,13])
seq_id=tmp[1,2]  # for final output table
cmd=paste0("system(\"echo '>target_seq\n",tar_seq,"' > ",dbdir,"/target_",mouse_prot_id,".fasta\")")
eval(parse(text=cmd))
cmd=paste0("system(\"",wdir,"/ncbi-blast-2.8.1+/bin/blastp -db ",dbdir,"/db/nr/Mmusculus_Ensembl -query ",dbdir,"/target_",mouse_prot_id,".fasta -out ",dbdir,"/tmp_",mouse_prot_id,".txt -num_threads 2 -outfmt '6 stitle' -max_target_seqs 1\")") # IMPORTANT !!! new version -- information left out for recollection: sframe sgi sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq
eval(parse(text=cmd))

rec.id=get_match_id(paste0(dbdir,"/tmp_",mouse_prot_id,".txt"))
if( rec.id==mouse_gene_id ){ presence=3; seq_name=seq_id }  # 1:1 orthology detected
}

} 
system(paste0("rm ",dbdir,"/target_",mouse_prot_id,".fasta"))
return(list(presence,seq_name))

}


# MAIN SCRIPT

# recall variables

#res=read.table("orthos_1to1_absence_presence.txt",h=T,row.names=1,stringsAsFactors = FALSE) # 1st col: gene name
#res2=read.table("orthos_1to1_gene_names.txt",h=T,row.names=1,stringsAsFactors = FALSE)
#res=read.table(paste0("tmp1/",inputf),h=F,row.names=1,stringsAsFactors = FALSE) # first column are gene names
#res2=read.table(paste0("tmp2/",inputf),h=F,row.names=1,stringsAsFactors = FALSE)


#for debug. ace1: aaa3 &aaa4; ace2: aaa1 & aaa2
res=read.table("aaa3",h=F,row.names=1,stringsAsFactors = FALSE) # debug 29/05/2019
res2=read.table("aaa4",h=F,row.names=1,stringsAsFactors = FALSE) # debug 29/05/2019
#res=read.table("aaa1",h=F,row.names=1,stringsAsFactors = FALSE) # debug 29/05/2019
#res2=read.table("aaa2",h=F,row.names=1,stringsAsFactors = FALSE) # debug 29/05/2019

species.id.list=c("Dnovemcinctus","Mdomestica","Ggallus","Xtropicalis","Drerio","Olatipes","Cmilii","Pmarinus","Eburgeri","Cintestinalis","Bfloridae","Dmelanogaster","Celegans")
#for recollection, associated NCBI taxonomy ids: 9361, 13616, 9031,8364,7955,8090,7868,7757,7764,7719,7739,7227,6239

prot.ids=res2[,14] # prot ids
gene.ids=row.names(res2) # gene ids
res2=res2[,1:13] # removes last column

cores=detectCores(); cl=makeCluster(as.numeric(args[2]))  # ncores threads 
registerDoParallel(cl)

res_table=foreach(i=1:nrow(res2), .combine='rbind', .packages=c("stringr")) %dopar% { # chunks table

tab.01=as.numeric(res[i,]) # boolean values of 1:1 orthologues for gene i
tab.names=as.character(res2[i,])  # Ensembl gene ids of 1:1 orthologues for gene i 

token=0
prot.id=prot.ids[i]; gene.id=gene.ids[i]

for(j in 1:13){

system(paste0("echo 'inside the j loop number: ", j,"'>> Rreport.txt"))  

#if( j==13 || (tab.01[j]==0 && ( (j==1) || (sum(tab.01[(j+1):13])>0) || (tab.01[j-1]!=0) )) ){ 
system(paste0("echo 'strange case at species: ", j,"'>> Rreport.txt"))  
    # A/ j==11 every time -- C.elegans has to be done systematically (biomaRt 1:1 orthologues are not reliable)
    # B/ 1st species already equal to zero; 
    # C/ if no 1:1 orthologues according to biomaRt but more "ancestral" species have it, recip. blast to test for 1:1 orthology. 
    # D/ An orthologue present to the left, & no orthologues to the right. --> tries to extend the orthology to the right
if(token==0){ #&& length(list.files(pattern = "\\.fasta$"))>0){
system(paste0("rm ",dbdir,"/",prot.id,".fasta"))
system(paste0("rm ",dbdir,"/*",prot.id,".fasta"))
}
token=1

if(!file.exists(paste0(dbdir,"/",prot.id,".fasta"))){
cmd=paste0("system(\"",wdir,"/ncbi-blast-2.8.1+/bin/blastdbcmd -entry ",prot.id," -db ",dbdir,"/db/nr/Mmusculus_Ensembl > ",dbdir,"/",prot.id,".fasta\")") # retrieves sequence from a local mouse UniProt proteome 
eval(parse(text=cmd))
}

system(paste0("echo 'protein id passed to function: ", prot.id,"'>> Rreport.txt"))  

out=ortho.blast(gene.id, prot.id, species.id.list[j]) # reciprocal blast function
#mouse_gene_id=gene.id; mouse_prot_id=prot.id; species_id=species.id.list[j]  # for debug

tab.names[j]=out[[2]] # sequences found with BLAST -- could find an id or name in the BLAST output

#} # end of if( j==11 || 

} # end of for(j in 1:11){


tab.names # variable passed to the parfor output  

}

write.table(res_table,"aaa2_out",row.names=T,sep="\t",col.names=F,quote=F) #debug
#write.table(res_table,paste0("tmp2/",inputf,"_2"),row.names=T,sep="\t",col.names=F,quote=F)


print(paste0("time finished: ", Sys.time()))







