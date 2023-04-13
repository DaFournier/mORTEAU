#!/usr/bin/Rscript

# Author: David Fournier - March/June 2019
# a script to compute 1:1 orthologues of all Ensembl mouse proteins
# to launch with main_script.sh (for 1 chunk of data - run the entire biomaRt proteome w. 11ortho_find.sh)
# or directly from command line with: 
# module load lang/R; > ../Rreport.txt; srun -e err.out -p short -A jgu-cbdm --mem=30000M -c 30 --time=01:00:00 Rscript ../../../blast.R xab 15 & 
# > ../Rreport.txt; srun -p short -A jgu-cbdm R --no-save --slave -f ./blast.R &

# LIBRARIES

#library("biomaRt")  # not used
library("foreach")
library("doParallel") # local installation in ~/R/x86_64-pc-linux-gnu-library/3.5
#library("data.table")  # not used
 library("stringr")

# GLOBAL VARIABLES

args=commandArgs(TRUE)
inputf=args[1]
ncores=args[2]

wdir="/project/jgu-cbdm/andradeLab/scratch/dfournie/prot_evol/orthologues_finding/blastp" # location of blast software suite
dbdir=wdir  # location of the sequences data bases

#system("> ../Rreport.txt")
system(paste0("echo 'BEGINNING OF SCRIPT: ", Sys.time(),"' >> ../Rreport.txt"))  

# FUNCTIONS 

get_match_id <- function(blast_output_file) {
tmp=readLines(blast_output_file)[1]
rec.id.tmp="" # variable initiation
if(length(tmp)>0){ # retrieving Ensembl gene sequence corres. to the protein in blastp output
tmp2=str_match(tmp, "gene:([A-Z0-9]+)") # pattern b. "gene:" and a dot (Ensembl gene has a version with it, like gene:ENSMUSG00000020681.14) or a space (if gene with no version like gene:ENSMUSG00000020681)
if(length(tmp2)>1){
rec.id.tmp=tmp2[,2] 
}
}
return(rec.id.tmp)
}


ortho.blast.PEP <- function(mouse_gene_id,mouse_prot_id, species_id) {  # script to find 1:1 orthologue sequence by doing a reciprocal blast (blastp,blastp) on protein seqs computed from predicted cds of Ensembl/NCBI genomes

 seq_name="no_ortho" # in case no orthologue present even on transcriptome sequencing data
 seq_id="" # intermediary variable
# step1: blastp against Ensembl/NCBI predicted protein sequences from WGS genomes
cmd=paste0("system(\"",wdir, "/ncbi-blast-2.8.1+/bin/blastp -db ",dbdir,"/db/nr/",species_id,"_PEP -query ../",mouse_prot_id,".fasta -out ../tmp_",mouse_prot_id,".txt -num_threads 2 -outfmt '6 evalue saccver stitle sseq' -max_target_seqs 1\")") # blastp to query target species protein seqs from Ensembl    --- qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore
eval(parse(text=cmd))

if(length(readLines(paste0("../tmp_",mouse_prot_id,".txt")))>0){ # homologue found
tmp=read.table(paste0("../tmp_",mouse_prot_id,".txt"),h=F,stringsAsFactors = FALSE) 
if(tmp[1,1]<=0.01){ # threshold to avoid bogus hits
# step2: blastp of the result nucleotide sequence against mouse protein db
tar_seq=gsub('-', '', tmp[1,ncol(tmp)])
if(species_id=="Bfloridae"){
seq_id=tmp[1,2]  # for final table - sequence id for Branchiostoma
}else{
seq_id=tmp[1,4]  # other species: coordinate in the WGS data
}

cmd=paste0("system(\"echo '>target_seq\n",tar_seq,"' > ../target_",mouse_prot_id,".fasta\")")
eval(parse(text=cmd))
cmd=paste0("system(\"",wdir,"/ncbi-blast-2.8.1+/bin/blastp -db ",dbdir,"/db/nr/Mmusculus_Ensembl -query ../target_",mouse_prot_id,".fasta -out ../tmp_",mouse_prot_id,".txt -num_threads 2 -outfmt '6 evalue stitle' -max_target_seqs 10\")") # blastp possible as we get the predicted protein sequence
eval(parse(text=cmd))

if(length(readLines(paste0("../tmp_",mouse_prot_id,".txt")))>0){ # gets a result from reciprocal blast
tmp2=readLines(paste0("../tmp_",mouse_prot_id,".txt")); good=str_detect(tmp2,"GRCm38:[0-9XYMT]+:") # keeps only canonical chromosomes
tmp=read.table(paste0("../tmp_",mouse_prot_id,".txt"),sep="",fill=T,h=F,stringsAsFactors = FALSE); tmp=tmp[which(good==TRUE),]
if(tmp[1,1]<=0.01){ # threshold to avoid bogus hits
write.table(tmp,paste0("../tmp_",mouse_prot_id,".txt"),row.names=F,sep="\t",col.names=F,quote=F)  
rec.id=get_match_id(paste0("../tmp_",mouse_prot_id,".txt"))
if( rec.id==mouse_gene_id ){ seq_name=seq_id }  # 1:1 orthology detected
}
}

} ## end of if(tmp[1,1]<=0.01)
}  ## end of if(length(readLines))

  return(seq_name)  # value "no_ortho" by default
}

ortho.blast <- function(mouse_gene_id, mouse_prot_id,species_id) {  # script to find 1:1 orthologue sequence by doing a reciprocal blast (blastp,blastp) 
# mouse gene id needed when doing recriprocal blast to avoid ambiguity b. multiple protein ids
presence=0; seq_name="no_ortho"; rec.id="" # values by default

# choice 1 step1: blastp 
cmd=paste0("system(\"",wdir, "/ncbi-blast-2.8.1+/bin/blastp -db ",dbdir,"/db/nr/",species_id," -query ../",mouse_prot_id,".fasta -out ../tmp_",mouse_prot_id,".txt -num_threads 2 -outfmt '6 evalue sgi sseq' -max_target_seqs 1\")")  # parameters left aside: qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore
eval(parse(text=cmd))

if(length(readLines(paste0("../tmp_",mouse_prot_id,".txt")))>0){ # homologue found
tmp=read.table(paste0("../tmp_",mouse_prot_id,".txt"),h=F,stringsAsFactors = FALSE) 
if(tmp[1,1]<=0.01){ # threshold to avoid bogus hits
# choice 1 step2: blastp back against mouse proteome
tar_seq=gsub('-', '', tmp[1,3])
seq_id=tmp[1,2]  # for final output table
cmd=paste0("system(\"echo '>target_seq\n",tar_seq,"' > ../target_",mouse_prot_id,".fasta\")")
eval(parse(text=cmd))
cmd=paste0("system(\"",wdir,"/ncbi-blast-2.8.1+/bin/blastp -db ",dbdir,"/db/nr/Mmusculus_Ensembl -query ../target_",mouse_prot_id,".fasta -out ../tmp_",mouse_prot_id,".txt -num_threads 2 -outfmt '6 evalue stitle' -max_target_seqs 10\")") # blastp againts mouse Ensembl protein sequences 
eval(parse(text=cmd))

if(length(readLines(paste0("../tmp_",mouse_prot_id,".txt")))>0){ # reciprocal blast target detected  # NEW !! (+ next 2 lines)
tmp2=readLines(paste0("../tmp_",mouse_prot_id,".txt")); good=str_detect(tmp2,"GRCm38:[0-9XYMT]+:") # keeps only canonical chromosomes
tmp=read.table(paste0("../tmp_",mouse_prot_id,".txt"),sep="",fill=T,h=F,stringsAsFactors = FALSE); tmp=tmp[which(good==TRUE),]
if(tmp[1,1]<=0.01){ # threshold to avoid bogus hits
write.table(tmp,paste0("../tmp_",mouse_prot_id,".txt"),row.names=F,sep="\t",col.names=F,quote=F) 
rec.id=get_match_id(paste0("../tmp_",mouse_prot_id,".txt"))
if( (species_id!="Cmilii") && ((length(tmp)==0) || (rec.id!=mouse_gene_id)) ){  # no predicted protein sequences for shark
# case where no orthologue in transcriptome data, we try prot seqs computed from predicted Ensembl/NCBI CDS data 
seq_id=ortho.blast.PEP(mouse_gene_id,mouse_prot_id,species_id)
if(seq_id!="no_ortho"){ rec.id=mouse_gene_id } # orthologue was found by ortho.blast.PEP
}
}
} # end of if(length(readLines(paste0("../tmp_",mouse_prot_id,".txt")))>0)  # NEW !!!

}else{ ## evalue>0.01
if( species_id!="Cmilii"){  # no predicted protein sequences for shark
seq_id=ortho.blast.PEP(mouse_gene_id,mouse_prot_id,species_id)
if(seq_id!="no_ortho"){ rec.id=mouse_gene_id } # orthologue was found by ortho.blast.PEP
}
}

}else{  # no homologue found with the initial blastp query
if( species_id!="Cmilii"){  # no predicted protein sequences for shark
seq_id=ortho.blast.PEP(mouse_gene_id,mouse_prot_id,species_id)
if(seq_id!="no_ortho"){ rec.id=mouse_gene_id } # orthologue was found by ortho.blast.PEP
}

} ## end of if(length(readLines))

if( rec.id==mouse_gene_id ){ presence=3; seq_name=seq_id }  # 1:1 orthology detected

# not needed because target_ files are only created when blast produces a significant output (w. an aligned sequence associated)
#if(file.exists(paste0("../target_",mouse_prot_id,".fasta"))){
#	system(paste0("rm ../target_",mouse_prot_id,".fasta"))
#}

return(list(presence,seq_name))

}

ortho.blast.WGS <- function(mouse_gene_id,mouse_prot_id, species_id) {  # script to find 1:1 orthologue sequence by doing a reciprocal blast (tblastn,blastp) on whole-genome sequencing data from Ensembl or NCBI repository

 seq_name="no_ortho" # in case no orthologue present even on transcriptome sequencing data
 seq_id="" # intermediary variable
# step1: blastp against Ensembl/NCBI predicted protein sequences from WGS genomes
cmd=paste0("system(\"",wdir, "/ncbi-blast-2.8.1+/bin/tblastn -db ",dbdir,"/db2/nt/",species_id,"_WGS -query ../",mouse_prot_id,".fasta -out ../tmp_",mouse_prot_id,".txt -num_threads 3 -outfmt '6 evalue saccver stitle sseq' -max_target_seqs 1\")") # tblastn to query WGS data from Ensembl or NCBI of target species protein seqs 
eval(parse(text=cmd))

if(length(readLines(paste0("../tmp_",mouse_prot_id,".txt")))>0){ # homologue found
tmp=read.table(paste0("../tmp_",mouse_prot_id,".txt"),h=F,stringsAsFactors = FALSE) 
if(tmp[1,1]<=0.01){ # threshold to avoid bogus hits
# step2: blastp of the result nucleotide sequence against mouse protein db
tar_seq=gsub('-', '', tmp[1,ncol(tmp)])
if(species_id %in% c("Bfloridae","Cmilii","Dnovemcinctus","Drerio")){
seq_id=tmp[1,2]  # for final table - sequence ids for some species
}else{
seq_id=tmp[1,4]  # other species: coordinate in the WGS data
}

cmd=paste0("system(\"echo '>target_seq\n",tar_seq,"' > ../target_",mouse_prot_id,".fasta\")")
eval(parse(text=cmd))
cmd=paste0("system(\"",wdir,"/ncbi-blast-2.8.1+/bin/blastp -db ",dbdir,"/db/nr/Mmusculus_Ensembl -query ../target_",mouse_prot_id,".fasta -out ../tmp_",mouse_prot_id,".txt -num_threads 2 -outfmt '6 evalue stitle' -max_target_seqs 10\")") # blastp possible as we get the predicted protein sequence
eval(parse(text=cmd))

if(length(readLines(paste0("../tmp_",mouse_prot_id,".txt")))>0){ # gets a result from reciprocal blast
tmp2=readLines(paste0("../tmp_",mouse_prot_id,".txt")); good=str_detect(tmp2,"GRCm38:[0-9XYMT]+:") # keeps only canonical chromosomes
tmp=read.table(paste0("../tmp_",mouse_prot_id,".txt"),sep="",fill=T,h=F,stringsAsFactors = FALSE); tmp=tmp[which(good==TRUE),]
if(tmp[1,1]<=0.01){ # threshold to avoid bogus hits
write.table(tmp,paste0("../tmp_",mouse_prot_id,".txt"),row.names=F,sep="\t",col.names=F,quote=F)  
rec.id=get_match_id(paste0("../tmp_",mouse_prot_id,".txt"))
if( rec.id==mouse_gene_id ){ seq_name=seq_id }  # 1:1 orthology detected
}
}

} ## end of if(tmp[1,1]<=0.01)
}  ## end of if(length(readLines))

  return(seq_name)  # value "no_ortho" by default
}


# MAIN SCRIPT

# recall variables

#res=read.table("orthos_1to1_absence_presence.txt",h=T,row.names=1,stringsAsFactors = FALSE) # 1st col: gene name
#res2=read.table("orthos_1to1_gene_names.txt",h=T,row.names=1,stringsAsFactors = FALSE)

res=read.table(inputf,h=F,row.names=1,stringsAsFactors = FALSE) # we are in the "presence" folder (see main_script.sh)
res2=read.table(paste0("../names/",inputf),h=F,row.names=1,stringsAsFactors = FALSE) # "names" folder

if(FALSE){ # debug
res=read.table("xaa_1_5_presence",h=F,row.names=1,stringsAsFactors = FALSE) 
res2=read.table("xaa_1_5_names",h=F,row.names=1,stringsAsFactors = FALSE) 
res=read.table("xaa_1_presence",h=F,row.names=1,stringsAsFactors = FALSE) 
res2=read.table("xaa_1_names",h=F,row.names=1,stringsAsFactors = FALSE)

res=read.table("aaa1",h=F,row.names=1,stringsAsFactors = FALSE) 
res2=read.table("aaa2",h=F,row.names=1,stringsAsFactors = FALSE)
res=read.table("aaa3",h=F,row.names=1,stringsAsFactors = FALSE)
res2=read.table("aaa4",h=F,row.names=1,stringsAsFactors = FALSE)
res=read.table("aaa5",h=F,row.names=1,stringsAsFactors = FALSE)
res2=read.table("aaa6",h=F,row.names=1,stringsAsFactors = FALSE)
}

species.id.list=c("Dnovemcinctus","Mdomestica","Ggallus","Xtropicalis","Drerio","Olatipes","Cmilii","Pmarinus","Eburgeri","Cintestinalis","Bfloridae","Dmelanogaster","Celegans")
#for recollection, associated NCBI taxonomy ids: 9361, 13616, 9031,8364,7955,8090,7868,7757,7764,7719,7739,7227,6239

prot.ids=res2[,14] # prot ids
gene.ids=row.names(res2) # gene ids
res2=res2[,1:13] # removes last column

cores=detectCores(); cl=makeCluster(as.numeric(ncores))  # ncores threads 
registerDoParallel(cl)

#res_table=res2
#foreach(i=1:nrow(res2), .packages=c("stringr")) %dopar% { # chunks table  # new method - does not work apparently 

#system("> ../tmp_table_progress.txt")  # DEBUG - TO DELETE

system(paste0("echo 'number of ncores: ",ncores,"'")) # debug

res_table=res2
n=as.integer((nrow(res_table)))/as.integer(ncores)   # number of chunk of size "ncores" -- important: size of a chunk has to be a multiple of number of cores allowed

#system("echo 'BEFORE for(k in ...) loop '") # debug

for(k in 0:(n-1)){ # subprocess (to accelerate the script)
interval=(10*k+1):(10*(k+1))

res_table_intermediate=foreach(i=interval, .combine='rbind',.packages=c("stringr")) %dopar% { # too slow...ramdisk?

tab.01=as.numeric(res[i,]) # boolean values of 1:1 orthologues for gene i
tab.names=as.character(res2[i,])  # Ensembl gene ids of 1:1 orthologues for gene i 

#token=0
prot.id=prot.ids[i]; gene.id=gene.ids[i]
system(paste0("echo 'gene: ", gene.id," AT TIME: ", Sys.time(),"'>> ../Rreport.txt")) # debug

for(j in 1:13){

#system(paste0("echo 'species: ", j,"'>> ../Rreport.txt")) # DEBUG - TO DELETE

if(!file.exists(paste0("../",prot.id,".fasta"))){
cmd=paste0("system(\"",wdir,"/ncbi-blast-2.8.1+/bin/blastdbcmd -entry ",prot.id," -db ",dbdir,"/db/nr/Mmusculus_Ensembl > ../",prot.id,".fasta\")") # retrieves sequence from a local mouse UniProt proteome 
eval(parse(text=cmd))
}

out=ortho.blast(gene.id, prot.id, species.id.list[j]) # reciprocal blast function
#mouse_gene_id=gene.id; mouse_prot_id=prot.id; species_id=species.id.list[j]  # for debug

tab.names[j]=out[[2]] # sequences found with BLAST -- could find an id or name in the BLAST output

#write.table(tab.names,"../tmp_table.txt",row.names=T,sep="\t",col.names=F,quote=F)   # DEBUG - TO DELETE
#system("cat ../tmp_table.txt >> ../tmp_table_progress.txt")   # DEBUG - TO DELETE

#} # end of if( j==13 || 

} # end of for(j in 1:13)

# optional tests using tblastn against WGS data


tab.names.tmp=tab.names

if(TRUE){ # set to FALSE for debug

tab.01=rep(0,13)

for(j in 1:13){ if(tab.names.tmp[j]!="no_ortho"){ tab.01[j]=1} }

for(j in 1:13){ 

token=0; token2=0
if(j!=13){ token=sum(tab.01[(j+1):13]) }; if(j==9){if((tab.01[j-4]!=0)||(tab.01[j-2]!=0)){ token2=1 }}

if( (tab.01[j]==0) && ( (j==1) || (tab.01[j-1]!=0) || (token>0) || (token2==1) ) ){
#system(paste0("echo 'ortho.WGS performing for ", species.id.list[j]," ...' >> ../Rreport.txt"))  
    # A/ 1st species already equal to zero; 
    # B/ if no 1:1 orthos according to biomaRt but more "ancestral" sp have it, recip. blast to test for 1:1 orthology. 
    # C/ An orthologue present to the left, & no orthologues to the right. --> tries to extend the orthology to the right
    # D/ same as C for E.burgeri, goes for the tblastn if sequence found in C.milii or D.rerio
    seq_id=ortho.blast.WGS(gene.id, prot.id, species.id.list[j]) # blast against WGS data
    if(seq_id!="no_ortho"){ tab.names.tmp[j]=seq_id; tab.01[j]=1 }
} ## end of if( (tab.01[j]==0

} # end of for(j in 1:13)

tab.names=tab.names.tmp

} # end of if(FALSE){ # toggled off for debug

system(paste0("rm ../",prot.id,".fasta"))
system(paste0("rm ../target_",prot.id,".fasta"))
system(paste0("rm ../tmp_",prot.id,".txt"))

 tab.names # rbind method: very slow......

} # end of for each ... do par loop

res_table[interval,]=res_table_intermediate

} # end of for(k in 0:(n-1))

system(paste0("echo 'time final table is built: ", Sys.time(),"' >> ../Rreport.txt"))

row.names(res_table)=gene.ids
write.table(res_table,paste0("../",inputf,"_out.txt"),row.names=T,sep="\t",col.names=F,quote=F)

stopCluster(cl)

system(paste0("echo 'time script is FINISHED: ", Sys.time(),"' >> ../Rreport.txt"))


####################################
## OLD STUFF
#####################################

#}else{ # evalue>0.01 -- apparently duplicated
#if( species_id!="Cmilii"){  # no predicted protein sequences for shark 
#seq_id=ortho.blast.PEP(mouse_gene_id,mouse_prot_id,species_id) 
#if(seq_id!="no_ortho"){ rec.id=mouse_gene_id } # orthologue was found by ortho.blast.PEP
#}

#} # end of if(tmp[1,1]<=0.01)







