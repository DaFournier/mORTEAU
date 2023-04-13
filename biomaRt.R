#!/usr/bin/Rscript


# Author: David Fournier - March 2019
# a script to extract 1:1 orthologues of mouse from biomaRt database
# first step of the pipeline
# mouse model is chosen because gene expression data for embryonic stages are available to test molecular evo-devo Ho

# to be launched from script launch_Rscript.sh



# LIBRARIES

library("biomaRt")


# MAIN SCRIPT

mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl") # mm10 

# initiating table with absence/presence of orthologues (nrow=nb of species)

# - retrieving info about mouse genome and proteome 
#proteins.mouse=getBM( attributes=c("ensembl_gene_id","transcript_length","ensembl_peptide_id_version","external_gene_name","gene_biotype") ,values="*",mart=mouse) # retrieving gene-protein correspondance OLD

proteins.mouse=getBM( attributes=c("ensembl_gene_id","ensembl_transcript_id","cds_length","gene_biotype","external_gene_name","refseq_peptide") ,values="*",mart=mouse) # retrieving gene-protein correspondance  
prot.trans.conversion=getBM( attributes=c("ensembl_transcript_id", "ensembl_peptide_id_version") ,values="*",mart=mouse)
prot.trans.conversion=subset(prot.trans.conversion, ensembl_transcript_id !="" & ensembl_peptide_id_version !="")

proteins.mouse=subset(proteins.mouse,gene_biotype=="protein_coding")
proteins.mouse=subset(proteins.mouse,!is.na(cds_length))
proteins.mouse=proteins.mouse[order(proteins.mouse$ensembl_gene_id),]
genes.unique.id=sort(unique(proteins.mouse$ensembl_gene_id)) #22376 unique spliced forms

# - selecting the longest cds for each gene (selecting on total transcript length with transcript_length leads to not selecting longest protein sometimes)
df=as.data.frame(matrix(0,nrow=length(genes.unique.id),ncol=5)) #22376 mouse genes 
colnames(df)=colnames(proteins.mouse)
n=1
for(current.id in genes.unique.id){
tmp=subset(proteins.mouse, ensembl_gene_id==current.id)
df[n,]=tmp[(which(tmp$cds_length==(max(tmp$cds_length))))[1], ] # selection of longest alternative spliced protein
n=n+1
}
paste0("nb of gene without prot sequence: ", nrow(subset(df, cds_length==""))) # test - checks if all gene seq has a prot seq
# converting transcript ids into peptide/protein sequence ids: 
proteins.mouse=merge(df,prot.trans.conversion, by="ensembl_transcript_id")
nrow(subset(proteins.mouse,ensembl_peptide_id_version=="")) # test - checks if all transcripts have their proteins sequence version associated 
proteins.mouse=proteins.mouse[order(proteins.mouse$ensembl_gene_id),] # sorting step

# retrieving orthologues from the Ensembl database

res=as.data.frame(matrix(0,nrow=nrow(df),ncol=13)) #22376 mouse genes (prot),13 species, 2 (cmilii & bfloridae) with no info about 1:1 orthology in Ensembl db
rownames(res)=df$ensembl_gene_id
species=c("dnovemcinctus","mdomestica","ggallus","xtropicalis","drerio","olatipes","cmilii","pmarinus","eburgeri","cintestinalis","bfloridae","dmelanogaster","celegans"); colnames(res)=species
res2=res

for(i in 1:13){  # all but cmilii & blforidae

sp=species[i]
print(paste0("working on species ... ",sp))

if(i in c(1:6,8:10,12:13)){
ensembl_gene=paste0(sp,"_homolog_ensembl_gene")
ensembl_prot=paste0(sp,"_homolog_ensembl_peptide")
orthology_type=paste0(sp,"_homolog_orthology_type")
tmp.bm=getBM( attributes=c("ensembl_gene_id",ensembl_gene,ensembl_prot,orthology_type),values="*",mart=mouse) 
cmd=paste0("tmp.bm.1to1=subset(tmp.bm,",orthology_type,"==\"ortholog_one2one\")") # keeps only 1to1 associations
eval(parse(text=cmd)) # dim 

# tests
test.dupl=(nrow(tmp.bm.1to1)==length(table(tmp.bm.1to1$ensembl_gene_id))) # tests if duplicates, normally none should be there (TRUE)
print(paste0("no duplicates for ",sp,": ",test.dupl)) # output of test for the log file
cmd=paste0("tmp=nrow(subset(tmp.bm.1to1,", ensembl_gene," != \"\" ))") # tests if number of gene sequence orthologues
eval(parse(text=cmd)); print(paste0("number of genes for ", sp," is: ",tmp))
cmd=paste0("tmp2=nrow(subset(tmp.bm.1to1,", ensembl_prot," != \"\" ))") # tests if number of protein sequence orthologues
eval(parse(text=cmd)); print(paste0("number of genes for ", sp," is: ",tmp2))

tmp.bm.1to1=subset(tmp.bm.1to1,ensembl_gene_id %in% rownames(res)) # removes non-coding sequences
j=tmp.bm.1to1[,1]
res[j,i]=rep(1,length(j))  #orthologues presence/absence
res2[j,i]=tmp.bm.1to1[,2]  #orthologues names

}else{ # Cmilii or Bfloridae
res[j,i]=0; res2[j,i]=0
}

} # end of for(i in 1:13)

#res=res[1:nrow(genes.mouse),]
#res2=res2[1:nrow(genes.mouse),]
write.table(res,"orthos_1to1_absence_presence.txt",row.names=T,sep="\t",col.names=T,quote=F)
write.table(res2,"orthos_1to1_gene_names.txt",row.names=T,sep="\t",col.names=T,quote=F)
write.table(proteins.mouse,"orthos_1to1_metadata.txt",row.names=F,sep="\t",col.names=T,quote=F)



