#!/usr/bin/Rscript


# Author: David Fournier - March/April 2019
# a script to extract 1:1 orthologues of mouse to fill putative missing data from biomaRt


# LIBRARIES

library("biomaRt")
library("foreach")
library("doParallel")
library("data.table")


# GLOBAL VARIABLES

#blast.dir="/project/jgu-cbdm/andradeLab/scratch/dfournie/prot_evol/orthologues_finding/blastp/ncbi-blast-2.8.1+/bin"
if(TRUE){ # old version
#wdir="/project/jgu-cbdm/andradeLab/scratch/dfournie/prot_evol/orthologues_finding/blastp" # location of blast software suite
wdir="/home/dfournie/dfournie/prot_evol/orthologues_finding/blastp" # location of blast software suite
dbdir=wdir
}
args=commandArgs(TRUE)
dbdir=args[1] # location of blast db on the ramdisk
inputf=args[3] # file chunk to process



# FUNCTIONS 

ortho.blast.WGS <- function(mouse_Ensembl_id, species_id) {  # script to find 1:1 ortologue sequence by doing a reciprocal blast on WGS data

 seq_id="no_ortho" # in case no orthologue present even on shotgun sequencing data
 seq_id2="" # intermediary variable
# step1: tblastn agains whole genome shotgun sequencing data
cmd=paste0("system(\"",wdir, "/ncbi-blast-2.8.1+/bin/tblastn -db ",dbdir,"/db/nt/",species_id,"_CDS -query ",dbdir,"/",mouse_Ensembl_id,".fasta -out ",dbdir,"/tmp_",mouse_Ensembl_id,".txt -num_threads 2 -outfmt '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq' -max_target_seqs 1\")") # tblastn to query target species whole genome shotgun sequencing
eval(parse(text=cmd))

if(length(readLines(paste0(dbdir,"/tmp_",mouse_Ensembl_id,".txt")))>0){ # homologue found
# step2: blastp of the result nucleotide sequence against mouse protein db
tmp=read.table(paste0(dbdir,"/tmp_",mouse_Ensembl_id,".txt"),h=F) # score and evalue don't matter here - could be low
tar_seq=gsub('-', '', tmp[1,13])
seq_id2=tmp[1,2]  # for final output table
cmd=paste0("system(\"echo '>target_seq\n",tar_seq,"' > ",dbdir,"/target_",mouse_Ensembl_id,".fasta\")")
eval(parse(text=cmd))
cmd=paste0("system(\"",wdir,"/ncbi-blast-2.8.1+/bin/blastp -db ",dbdir,"/db/nr/Mmusculus -query ",dbdir,"/target_",mouse_Ensembl_id,".fasta -out ",dbdir,"/tmp_",mouse_Ensembl_id,".txt -num_threads 2 -outfmt '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq' -max_target_seqs 1\")") # blastp possible as we get the predicted protein sequence
eval(parse(text=cmd))

if(length(readLines(paste0(dbdir,"/tmp_",mouse_Ensembl_id,".txt")))>0){
seq_id=seq_id2 # storage of id now that we know that orthologue 1:1 exists
tmp=read.table(paste0(dbdir,"/tmp_",mouse_Ensembl_id,".txt"),h=F) # score and evalue don't matter here - could be low
tar_seq=gsub('-', '', tmp[1,13])
cmd=paste0("system(\"echo '>target_seq\n",tar_seq,"' > ",dbdir,"/target_",mouse_Ensembl_id,".fasta\")")
eval(parse(text=cmd))
cmd=paste0("system(\"",wdir,"/ncbi-blast-2.8.1+/bin/blastp -query ",dbdir,"/",mouse_Ensembl_id,".fasta -subject ",dbdir,"/target_",mouse_Ensembl_id,".fasta -out ",dbdir,"/final_result_",mouse_Ensembl_id,".txt -outfmt 6 -max_target_seqs 1\")") # retrieves sequence and & pipe output to blast2p to check if pipeline input & output are same (signing 1:1 orthology)
eval(parse(text=cmd))
}
}
  return(seq_id)

}

ortho.blast <- function(mouse_Ensembl_id, species_id) {  # script to find 1:1 ortologue sequence by doing a reciprocal blast

presence=0; seq_name="no_ortho" # values by default

if(species_id %in% c("Dnovemcinctus","Mdomestica","Ggallus","Xtropicalis","Drerio","Cintestinalis","Dmelanogaster","Celegans")){
seq_id="" # for final output table - empty by default
# choice 1 step1: blastp for all but olatipes, pmarinus and eburgeri
cmd=paste0("system(\"",wdir, "/ncbi-blast-2.8.1+/bin/blastp -db ",dbdir,"/db/nr/",species_id," -query ",dbdir,"/",mouse_Ensembl_id,".fasta -out ",dbdir,"/tmp_",mouse_Ensembl_id,".txt -num_threads 2 -outfmt '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq' -max_target_seqs 1\")") # retrieves sequence from a local mouse UniProt proteome & pipe output to blastp to query target species proteome   # TO DO: PRINT ONLY SEQ IN OUTPUT AND PIPE TO NEXT COMMAND
eval(parse(text=cmd))

if(length(readLines(paste0(dbdir,"/tmp_",mouse_Ensembl_id,".txt")))>0){ # homologue found
# choice 1 step2: blastp back to mouse proteome
tmp=read.table(paste0(dbdir,"/tmp_",mouse_Ensembl_id,".txt"),h=F) # score and evalue don't matter here - could be low
#tar_seq=tmp[1,2] # sequence id with top-score
tar_seq=gsub('-', '', tmp[1,13])
seq_id=tmp[1,2]  # for final output table
cmd=paste0("system(\"echo '>target_seq\n",tar_seq,"' > ",dbdir,"/target_",mouse_Ensembl_id,".fasta\")")
eval(parse(text=cmd))
cmd=paste0("system(\"",wdir,"/ncbi-blast-2.8.1+/bin/blastp -db ",dbdir,"/db/nr/Mmusculus -query ",dbdir,"/target_",mouse_Ensembl_id,".fasta -out ",dbdir,"/tmp_",mouse_Ensembl_id,".txt -num_threads 2 -outfmt '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq' -max_target_seqs 1\")")
eval(parse(text=cmd))

if(length(readLines(paste0(dbdir,"/tmp_",mouse_Ensembl_id,".txt")))>0){
tmp=read.table(paste0(dbdir,"/tmp_",mouse_Ensembl_id,".txt"),h=F) # score and evalue don't matter here - could be low
tar_seq=gsub('-', '', tmp[1,13])
cmd=paste0("system(\"echo '>target_seq\n",tar_seq,"' > ",dbdir,"/target_",mouse_Ensembl_id,".fasta\")")
eval(parse(text=cmd))

# to change: find the external gene id and compare it in two sequences (avoids having alignments with paralogues)
cmd=paste0("system(\"",wdir,"/ncbi-blast-2.8.1+/bin/blastp -query ",dbdir,"/",mouse_Ensembl_id,".fasta -subject ",dbdir,"/target_",mouse_Ensembl_id,".fasta -out ",dbdir,"/final_result_",mouse_Ensembl_id,".txt -outfmt 6 -max_target_seqs 1\")") # retrieves sequence
eval(parse(text=cmd)) # only seq every time? no need to select the best score?

}


#if(( file.info(paste0(dbdir,"/tmp_",mouse_Ensembl_id,".txt"))$size == 0 || (file.info(paste0(dbdir,"/final_result_",mouse_Ensembl_id,".txt"))$size == 0)) && (species_id %in% c("Ggallus", "Drerio","Cintestinalis") )){ # case where no orthologue in transcriptome data, we try CDS predicted from WGS data  -- for Ggallus, Drerio & Cintestinalis -- OLD 
if( (file.info(paste0(dbdir,"/tmp_",mouse_Ensembl_id,".txt"))$size == 0) || (file.info(paste0(dbdir,"/final_result_",mouse_Ensembl_id,".txt"))$size == 0) ){ # case where no orthologue in transcriptome data, we try CDS predicted from WGS data  -- for any species
print(paste0("trying to do on WGS data of ",species_id,"..."))
seq_id=ortho.blast.WGS(mouse_Ensembl_id,species_id)
# convert transcript into gene id: to do after getting the final table
}

if(file.exists(paste0(dbdir,"/final_result_",mouse_Ensembl_id,".txt"))){
if(file.info(paste0(dbdir,"/final_result_",mouse_Ensembl_id,".txt"))$size> 0){
reciprocal.search=read.table(paste0(dbdir,"/final_result_",mouse_Ensembl_id,".txt"),h=F)
if(reciprocal.search[1,11]<1e-4){ presence=3; seq_name=seq_id } # test if evalue=0 alignment initial mouse query/final mouse sequence from reciprocal  BLAST  -- VALUE 3 FOR DEBUGGING, TO CHANGE TO 1 IN FUTURE
} # if output file is empty, the 1:1 orthology is considered not to exist
}

}

}else{ # Olatipes, Pmarinus & Eburgeri
# choice 2 step1: tblastn for olatipes, pmarinus and eburgeri
system(paste0("echo 'about to tblastn against ", dbdir,"/",mouse_Ensembl_id,".fasta' >> Rreport.txt"))  
#print(paste0("about to tblastn against ",dbdir,"/",mouse_Ensembl_id,".fasta"))
system(paste0("cat ", dbdir,"/",mouse_Ensembl_id,".fasta >> Rreport.txt"))  
#paste0("system(\"cat ",dbdir,"/",mouse_Ensembl_id,".fasta\")")
system(paste0("echo 'current species: ", species_id,"' >> Rreport.txt"))  
#print(paste0("species being ", species_id))

cmd=paste0("system(\"",wdir, "/ncbi-blast-2.8.1+/bin/tblastn -db ",dbdir,"/db/nt/",species_id,"_CDS -query ",dbdir,"/",mouse_Ensembl_id,".fasta -out ",dbdir,"/tmp_",mouse_Ensembl_id,".txt -num_threads 2 -outfmt '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq' -max_target_seqs 1\")") # tblastn to query target species transcriptome
eval(parse(text=cmd))

if(length(readLines(paste0(dbdir,"/tmp_",mouse_Ensembl_id,".txt")))>0){ # homologue found
# choice 2 step2: blastp of the result nucleotide sequence against mouse protein db
tmp=read.table(paste0(dbdir,"/tmp_",mouse_Ensembl_id,".txt"),h=F) # score and evalue don't matter here - could be low
#tar_seq=tmp[1,2] # sequence id with top-score
tar_seq=gsub('-', '', tmp[1,13])
seq_id=tmp[1,2]  # for final output table
cmd=paste0("system(\"echo '>target_seq\n",tar_seq,"' > ",dbdir,"/target_",mouse_Ensembl_id,".fasta\")")
eval(parse(text=cmd))
cmd=paste0("system(\"",wdir,"/ncbi-blast-2.8.1+/bin/blastp -db ",dbdir,"/db/nr/Mmusculus -query ",dbdir,"/target_",mouse_Ensembl_id,".fasta -out ",dbdir,"/tmp_",mouse_Ensembl_id,".txt -num_threads 2 -outfmt '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq' -max_target_seqs 1\")") # blastp possible as we get the predicted protein sequence
eval(parse(text=cmd))
#cmd=paste0("system(\"",wdir,"/ncbi-blast-2.8.1+/bin/blastx -db ",dbdir,"/db/nr/Mmusculus -query ",dbdir,"/target_",mouse_Ensembl_id,".fasta -out ",dbdir,"/tmp_",mouse_Ensembl_id,".txt -num_threads 2 -outfmt '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq' -max_target_seqs 1\")") # blastx: translates a nucleotide query in six frames and searches it against a protein database.

if(length(readLines(paste0(dbdir,"/tmp_",mouse_Ensembl_id,".txt")))>0){
tmp=read.table(paste0(dbdir,"/tmp_",mouse_Ensembl_id,".txt"),h=F) # score and evalue don't matter here - could be low
tar_seq=gsub('-', '', tmp[1,13])
cmd=paste0("system(\"echo '>target_seq\n",tar_seq,"' > ",dbdir,"/target_",mouse_Ensembl_id,".fasta\")")
eval(parse(text=cmd))
cmd=paste0("system(\"",wdir,"/ncbi-blast-2.8.1+/bin/blastp -query ",dbdir,"/",mouse_Ensembl_id,".fasta -subject ",dbdir,"/target_",mouse_Ensembl_id,".fasta -out ",dbdir,"/final_result_",mouse_Ensembl_id,".txt -outfmt 6 -max_target_seqs 1\")") # retrieves sequence and & pipe output to blast2p to check if pipeline input & output are same (signing 1:1 orthology)
eval(parse(text=cmd))
#cmd=paste0("system(\"",wdir,"/ncbi-blast-2.8.1+/bin/blastdbcmd -entry ",tar_seq," -db ",dbdir,"/db/nr/Mmusculus -outfmt '%f' | ",wdir,"/ncbi-blast-2.8.1+/bin/blastp -query ",dbdir,"/",mouse_Ensembl_id,".fasta -subject - -out ",dbdir,"/final_result_",mouse_Ensembl_id,".txt -outfmt 6\")") # retrieves sequence and & pipe output to blast2p to check if pipeline input & output are same (signing 1:1 orthology) # alternative: take full sequence of target - slower but should give better evalues

if(file.exists(paste0(dbdir,"/final_result_",mouse_Ensembl_id,".txt"))){
if(file.info(paste0(dbdir,"/final_result_",mouse_Ensembl_id,".txt"))$size> 0){
reciprocal.search=read.table(paste0(dbdir,"/final_result_",mouse_Ensembl_id,".txt"),h=F)
if(reciprocal.search[1,11]<1e-4){ presence=3; seq_name=seq_id } # test if evalue=0 alignment initial mouse query/final mouse sequence from reciprocal BLAST  -- VALUE 2 FOR DEBUGGING, TO CHANGE TO 1 IN FUTURE
} # if output file is empty, the 1:1 orthology is considered not to exist
}

}
}

}

system(paste0("rm ",dbdir,"/target_",mouse_Ensembl_id,".fasta"))

  return(list(presence,seq_name))

}



# MAIN SCRIPT

# recall variables

#res=read.table("orthos_1to1_absence_presence.txt",h=T,row.names=1,stringsAsFactors = FALSE) # first column are gene names
#res2=read.table("orthos_1to1_gene_names.txt",h=T,row.names=1,stringsAsFactors = FALSE)
res=read.table(paste0("tmp1/",inputf),h=F,row.names=1,stringsAsFactors = FALSE) # first column are gene names
res2=read.table(paste0("tmp2/",inputf),h=F,row.names=1,stringsAsFactors = FALSE)

#res=read.table("aaa1",h=F,row.names=1,stringsAsFactors = FALSE) # debug 29/05/2019
#res2=read.table("aaa2",h=F,row.names=1,stringsAsFactors = FALSE) # debug 29/05/2019


#res2=read.table("test2.txt",h=T,row.names=1,stringsAsFactors = FALSE) # debug


#mouse.genes=colnames(res2)

species.id.list=c("Dnovemcinctus","Mdomestica","Ggallus","Xtropicalis","Drerio","Olatipes","Pmarinus","Eburgeri","Cintestinalis","Dmelanogaster","Celegans")
#for recollection, associated NCBI taxonomy ids: 9361, 13616, 9031,8364,7955,8090,7757,7764,7719,7227,6239

prot.ids=res2[,12] # prot ids
res2=res2[,1:11] # removes last column

cores=detectCores(); cl=makeCluster(as.numeric(args[2]))  # ncores threads 
registerDoParallel(cl)

res_table=foreach(i=1:nrow(res2), .combine='rbind') %dopar% { # chunks tabl

tmp=as.numeric(res[i,])
tmp2=res2[i,]

token=0
for(j in 1:10){

system(paste0("echo 'inside the j loop number: ", j,"'>> Rreport.txt"))  

id=as.character(prot.ids[i])
if(tmp[j]==0 && sum(tmp[(j+1):11])>0){ # if no 1:1 orthologues according to biomaRt but more "ancestral" species have it, we perform recip. blast to test for 1:1 orthology
if(token==0){ #&& length(list.files(pattern = "\\.fasta$"))>0){
system(paste0("rm ",dbdir,"/",id,".fasta"))
system(paste0("rm ",dbdir,"/*",id,".fasta"))
}
token=1
system(paste0("echo 'strange case at species: ", j,"'>> Rreport.txt"))  

if(!file.exists(paste0(dbdir,"/",id,".fasta"))){
cmd=paste0("system(\"",wdir,"/ncbi-blast-2.8.1+/bin/blastdbcmd -entry ",id," -db ",dbdir,"/db/nr/Mmusculus_Ensembl > ",dbdir,"/",id,".fasta\")") # retrieves sequence from a local mouse UniProt proteome 
eval(parse(text=cmd))
}

system(paste0("echo 'protein id passed to function: ", id,"'>> Rreport.txt"))  

out=ortho.blast(id, species.id.list[j]) # reciprocal blast function

tmp2[j]=as.character(out[[2]]) # sequences found with BLAST -- could find an id or name in the BLAST output

}

}

tmp2 # variable passed to the parfor output  

}

write.table(res_table,paste0("tmp2/",inputf,"_2"),row.names=T,sep="\t",col.names=F,quote=F) #debug










