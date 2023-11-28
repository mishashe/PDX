library(readxl)

cores <- 40
wd <- "/home/misha/Documents/Development/PDX/"
system(paste0("mdir -r ",wd,"commands/"))


# import samples info----
sheets <- excel_sheets(paste0(wd,"data/Samples.xlsx"))
SamplesInfo <- read_excel(paste0(wd,"data/Samples.xlsx"),sheet="Oncoprint")
names <- SamplesInfo$WES; names <- names[!is.na(names)]
IDs <- SamplesInfo$`Sample ID...6`; IDs <- IDs[!is.na(IDs)]


# rename ---
out_dir <- paste0("/home/m.sheinman/PDX/data/WES/fastq/renamed/")
commandfile <- paste0(wd,"commands/WES/rename.txt")
cat("cd ",wd,"\n",file = commandfile,append=FALSE)
cat("mkdir ",out_dir,"\n",file = commandfile,append=TRUE)
cat("rm ",out_dir,"*.fastq \n",file = commandfile,append=TRUE)
for (i in 1:length(names))
{
  command <- paste0("gzip -d /home/m.sheinman/PDX/data/WES/fastq/raw/*",names[i],"*_R1.fastq.gz --stdout >> ", "/home/m.sheinman/PDX/data/WES/fastq/renamed/",IDs[i],"_R1.fastq \n")
  cat(command,file = commandfile,append=TRUE)
  command <- paste0("gzip -d /home/m.sheinman/PDX/data/WES/fastq/raw/*",names[i],"*_R2.fastq.gz --stdout >> ", "/home/m.sheinman/PDX/data/WES/fastq/renamed/",IDs[i],"_R2.fastq \n")
  cat(command,file = commandfile,append=TRUE)
  command <- paste0("gzip /home/m.sheinman/PDX/data/WES/fastq/renamed/",IDs[i],"_R1.fastq \n")
  cat(command,file = commandfile,append=TRUE)
  command <- paste0("gzip /home/m.sheinman/PDX/data/WES/fastq/renamed/",IDs[i],"_R2.fastq \n")
  cat(command,file = commandfile,append=TRUE)
}




paired <- TRUE
# cutadapt ----
commandfile <- paste0(wd,"commands/WES/cutadapt.txt")
cat("cd ",wd,"\n",file = commandfile,append=FALSE)
input_dir <- "/home/m.sheinman/PDX/data/WES/fastq/renamed/"
out_dir <- paste0("/home/m.sheinman/PDX/data/WES/fastq/trimmed/")
cat("cd ",wd,"\n",file = commandfile,append=FALSE)
cat("mkdir ",out_dir,"\n",file = commandfile,append=TRUE)
for (i in 1:length(IDs))
{
  file1 <- paste0(input_dir,IDs[i],"_R1.fastq.gz \n")
  file2 <- paste0(input_dir,IDs[i],"_R2.fastq.gz \n")
  command <- paste0("cutadapt --cores=",1," -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -q 20 --minimum-length 60 -o ",out_dir,IDs[i],"_R1.fastq.gz -p ",out_dir,IDs[i],"_R2.fastq.gz ",input_dir,IDs[i],"_R1.fastq.gz ",input_dir,IDs[i],"_R2.fastq.gz > ",out_dir,IDs[i],".log & \n")
  cat(command,file = commandfile,append=TRUE)
}
cat("\n \n",file = commandfile,append=TRUE)
cat("\n \n",file = commandfile,append=TRUE)




# qc ----
commandfile <- paste0(wd,"commands/WES/qc.txt")
input_dir <- "/home/m.sheinman/PDX/data/WES/fastq/trimmed/"
out_dir <- "/home/m.sheinman/PDX/data/WES/fastq/"
cat("cd ",wd,"\n",file = commandfile,append=FALSE)
cat("mkdir -r ",out_dir,"\n",file = commandfile,append=TRUE)
for (i in 1:length(IDs))
{
  command <- paste0("fastqc ",input_dir,IDs[i],"_R1.fastq.gz ",input_dir,IDs[i],"_R2.fastq.gz --outdir ",out_dir," & \n")
  cat(command,file = commandfile,append=TRUE)
}
cat("echo 123 \n \n",file = commandfile,append=TRUE)
command <- paste0("multiqc ",out_dir," -o ",out_dir," \n")
cat(command,file = commandfile,append=TRUE)
cat("\n \n",file = commandfile,append=TRUE)






# bwa ----
index.host <- paste0("/home/m.sheinman/PDX/data/external/GCF_000001635.27_GRCm39_genomic/Mus_musculus.GRCm39.dna.toplevel.fa")
index.graft <- paste0("/home/m.sheinman/PDX/data/external/Homo_sapiens_assembly38/Homo_sapiens_assembly38.fasta")
commandfile <- paste0(wd,"commands/WES/bwa.txt")
input_dir <- "/home/m.sheinman/PDX/data/WES/fastq/trimmed/"
out_dir <- "/home/m.sheinman/PDX/data/WES/bam/aligned/"
cat("conda activate PDX3.5"," \n",file = commandfile,append=FALSE)
cat("cd ",wd," \n",file = commandfile,append=TRUE)
cat("mkdir -p ",out_dir,"\n",file = commandfile,append=TRUE)
for (i in 1:length(IDs))
{
  command <- paste0("bwa mem -t ",2," ",index.host," ",input_dir,IDs[i],"_R1.fastq.gz ",input_dir,IDs[i],"_R2.fastq.gz | samtools view -Sb -> ",out_dir,IDs[i],".host.bam & \n")
  cat(command,file = commandfile,append=TRUE)
  command <- paste0("bwa mem -t ",2," ",index.graft," ",input_dir,IDs[i],"_R1.fastq.gz ",input_dir,IDs[i],"_R2.fastq.gz | samtools view -Sb -> ",out_dir,IDs[i],".graft.bam & \n")
  cat(command,file = commandfile,append=TRUE)
}
cat("echo 123 \n \n",file = commandfile,append=TRUE)
cat("\n \n",file = commandfile,append=TRUE)


# sambamba sort ----
commandfile <- paste0(wd,"commands/WES/sambamba_sort.txt")
input_dir <- "/home/m.sheinman/PDX/data/WES/bam/aligned/"
out_dir <- "/home/m.sheinman/PDX/data/WES/bam/sorted/"
cat("cd ",wd," \n",file = commandfile,append=FALSE)
cat("conda activate PDX3.5"," \n",file = commandfile,append=FALSE)
cat("mkdir ",out_dir,"\n",file = commandfile,append=TRUE)
for (i in 1:length(IDs))
{
  command <- paste0("sambamba sort -t ",1," ",input_dir,IDs[i],".host.bam --out=",out_dir,IDs[i],".host.bam & \n")
  cat(command,file = commandfile,append=TRUE)
  command <- paste0("sambamba sort -t ",1," ",input_dir,IDs[i],".graft.bam --out=",out_dir,IDs[i],".graft.bam & \n")
  cat(command,file = commandfile,append=TRUE)
}
cat("echo 123 \n \n",file = commandfile,append=TRUE)
cat("\n \n",file = commandfile,append=TRUE)
                                      
                                      
# disambiguate ----
commandfile <- paste0(wd,"commands/WES/disambiguate.txt")
input_dir <- "/home/m.sheinman/PDX/data/WES/bam/sorted/"
out_dir <- "/home/m.sheinman/PDX/data/WES/bam/disambiguated/"
cat("cd ",wd," \n",file = commandfile,append=FALSE)
cat("conda activate PDX3.5"," \n",file = commandfile,append=TRUE)
cat("mkdir ",out_dir,"\n",file = commandfile,append=TRUE)
for (i in 1:length(IDs))
{
  unit <- SamplesInfo$unit[i]
  command <- paste0("ngs_disambiguate -a bwa -o ",out_dir," -s ",IDs[i]," ",input_dir,IDs[i],".graft.bam ", input_dir,IDs[i],".host.bam & \n")
  cat(command,file = commandfile,append=TRUE)
}
cat("echo 123 \n \n",file = commandfile,append=TRUE)
cat("\n \n",file = commandfile,append=TRUE)



# mark duplicates picard----
commandfile <- paste0(wd,"commands/WES/MarkDuplicates.txt")
input_dir <- "/home/m.sheinman/PDX/data/WES/bam/disambiguated/"
out_dir <- "/home/m.sheinman/PDX/data/WES/bam/markedDups/"
cat("cd ",wd," \n",file = commandfile,append=FALSE)
cat("conda activate PDX3.5"," \n",file = commandfile,append=TRUE)
cat("mkdir ",out_dir,"\n",file = commandfile,append=TRUE)
for (i in 1:length(IDs))
{
  command <- paste0("/home/m.sheinman/Software/jre1.8.0_391/bin/java -XX:ParallelGCThreads=40 -jar /home/m.sheinman/Software/picard.jar MarkDuplicates I=",input_dir,IDs[i],".disambiguatedSpeciesA.bam O=",out_dir,IDs[i],".bam M=",out_dir,IDs[i],".bam.txt  \n")
  cat(command,file = commandfile,append=TRUE)
}
cat("echo 123 \n \n",file = commandfile,append=TRUE)
cat("\n \n",file = commandfile,append=TRUE)

# GATK add groups  ----
commandfile <- paste0(wd,"commands/WES/addGrous.txt")
input_dir <- "/home/m.sheinman/PDX/data/WES/bam/markedDups/"
out_dir <- "/home/m.sheinman/PDX/data/WES/bam/addedGroups/"
cat("cd ",wd," \n",file = commandfile,append=FALSE)
cat("conda activate gatk"," \n",file = commandfile,append=TRUE)
cat("mkdir ",out_dir,"\n",file = commandfile,append=TRUE)
for (i in 1:length(IDs))
{
  command <- paste0("/home/m.sheinman/Software/jre1.8.0_391/bin/java -XX:ParallelGCThreads=1 -jar /home/m.sheinman/Software/picard.jar AddOrReplaceReadGroups I=",input_dir,IDs[i],".bam  O=",out_dir,IDs[i],".bam   RGID=4   RGLB=lib1   RGPL=ILLUMINA   RGPU=unit1  RGSM=20 & \n")
  cat(command,file = commandfile,append=TRUE)
}
cat("echo 123 \n \n",file = commandfile,append=TRUE)
cat("\n \n",file = commandfile,append=TRUE)


# run to index vcf files
# gatk IndexFeatureFile --feature-file /home/m.sheinman/PDX/data/external/Homo_sapiens_assembly38.dbsnp138.vcf
# gatk IndexFeatureFile --feature-file /home/m.sheinman/PDX/data/external/Mills_and_1000G_gold_standard.indels.hg38.vcf
# data is here https://hub.docker.com/r/adgh456/gatk-library
# GATK make tablr for recalibration  ----
commandfile <- paste0(wd,"commands/WES/recalibrateTable.txt")
input_dir <- "/home/m.sheinman/PDX/data/WES/bam/addedGroups/"
out_dir <- "/home/m.sheinman/PDX/data/WES/bam/addedGroups/"
cat("cd ",wd," \n",file = commandfile,append=FALSE)
cat("conda activate gatk"," \n",file = commandfile,append=TRUE)
cat("mkdir ",out_dir,"\n",file = commandfile,append=TRUE)
for (i in 1:length(IDs))
{
  command <- paste0("gatk BaseRecalibrator -I ",input_dir,IDs[i],".bam -R /home/m.sheinman/PDX/data/external/Homo_sapiens_assembly38/Homo_sapiens_assembly38.fasta --known-sites /home/m.sheinman/PDX/data/external/Homo_sapiens_assembly38.dbsnp138.vcf  --known-sites /home/m.sheinman/PDX/data/external/Mills_and_1000G_gold_standard.indels.hg38.vcf  -O ",out_dir,IDs[i],".txt & \n")
  cat(command,file = commandfile,append=TRUE)
}
cat("echo 123 \n \n",file = commandfile,append=TRUE)
cat("\n \n",file = commandfile,append=TRUE)

# GATK recalibration  ----
commandfile <- paste0(wd,"commands/WES/recalibrate.txt")
input_dir <- "/home/m.sheinman/PDX/data/WES/bam/addedGroups/"
out_dir <- "/home/m.sheinman/PDX/data/WES/bam/recalibrated/"
cat("cd ",wd," \n",file = commandfile,append=FALSE)
cat("conda activate gatk"," \n",file = commandfile,append=TRUE)
cat("mkdir ",out_dir,"\n",file = commandfile,append=TRUE)
for (i in 1:length(IDs))
{
  command <- paste0("gatk ApplyBQSR --java-options '-Xmx6G -Dsamjdk.compression_level=5' -R /home/m.sheinman/PDX/data/external/Homo_sapiens_assembly38/Homo_sapiens_assembly38.fasta  -I ",input_dir,IDs[i],".bam --bqsr-recal-file ",input_dir,IDs[i],".txt -O ",out_dir,IDs[i],".bam & \n")
  cat(command,file = commandfile,append=TRUE)
}
cat("echo 123 \n \n",file = commandfile,append=TRUE)
cat("\n \n",file = commandfile,append=TRUE)



# Mutect2  ----
commandfile <- paste0(wd,"commands/WES/Mutect2.txt")
input_dir <- "/home/m.sheinman/PDX/data/WES/bam/recalibrated/"
out_dir <- "/home/m.sheinman/PDX/data/WES/vcf/mutect2/"
cat("cd ",wd," \n",file = commandfile,append=FALSE)
cat("conda activate gatk"," \n",file = commandfile,append=TRUE)
cat("mkdir ",out_dir,"\n",file = commandfile,append=TRUE)
for (i in 1:length(IDs))
{

  command <- paste0("gatk Mutect2 -R /home/m.sheinman/PDX/data/external/Homo_sapiens_assembly38/Homo_sapiens_assembly38.fasta  -I ",input_dir,IDs[i],".bam -tumor 20 -O ",out_dir,IDs[i],".vcf.gz & \n")
  cat(command,file = commandfile,append=TRUE)
}
cat("echo 123 \n \n",file = commandfile,append=TRUE)
cat("\n \n",file = commandfile,append=TRUE)


# filter  ----
commandfile <- paste0(wd,"commands/WES/filter.txt")
input_dir <- "/home/m.sheinman/PDX/data/WES/vcf/mutect2/"
out_dir <- "/home/m.sheinman/PDX/data/WES/vcf/filtered/"
cat("cd ",wd," \n",file = commandfile,append=FALSE)
cat("conda activate gatk"," \n",file = commandfile,append=TRUE)
cat("mkdir ",out_dir,"\n",file = commandfile,append=TRUE)
for (i in 1:length(IDs))
{
  
  command <- paste0("gatk FilterMutectCalls -R /home/m.sheinman/PDX/data/external/Homo_sapiens_assembly38/Homo_sapiens_assembly38.fasta  -V ",input_dir,IDs[i],".vcf.gz  -O ",out_dir,IDs[i],".vcf.gz & \n")
  cat(command,file = commandfile,append=TRUE)
}
cat("echo 123 \n \n",file = commandfile,append=TRUE)
cat("\n \n",file = commandfile,append=TRUE)


# download annovar databases
# /home/m.sheinman/Software/annovar/annotate_variation.pl -buildver hg38 -downdb -webfrom annovar cytoBand  /home/m.sheinman/Software/annovar/humandb/
# hrcr1 gnomad_genome gnomad_exome esp6500siv2_aa esp6500siv2_ea exac03nontcga gme kaviar_20150923 1000g2015aug cytoBand clinvar_20221231
#/home/m.sheinman/Software/annovar/table_annovar.pl -buildver hg38 -out myanno -remove -protocol esp6500siv2_aa,esp6500siv2_ea,exac03nontcga,gme,kaviar_20150923,ALL.sites.2015_08 -operation 6 -nastring . -vcfinput -polish -outfile /home/m.sheinman/PDX/data/WES/annotation/annovar/IDC001.csv --xreffile /home/m.sheinman/Software/annovar/example/gene_xref.txt  /home/m.sheinman/PDX/data/WES/vcf/filtered/IDC001.vcf.gz humandb/
#/home/m.sheinman/Software/annovar/table_annovar.pl /home/m.sheinman/PDX/data/WES/vcf/filtered/IDC001.vcf.gz humandb/ -buildver hg38 -out myanno -remove -protocol refGene,exac03,avsnp147,dbnsfp30a -operation g,r,f,f,f -nastring . -vcfinput -polish
  
# annovar  ----
commandfile <- paste0(wd,"commands/WES/annovar.txt")
input_dir <- "/home/m.sheinman/PDX/data/WES/vcf/filtered/"
out_dir <- "/home/m.sheinman/PDX/data/WES/annotation/annovar/"
cat("cd ",wd," \n",file = commandfile,append=FALSE)
cat("mkdir ",out_dir,"\n",file = commandfile,append=TRUE)
for (i in 1:length(IDs))
{
  command <- paste0("/home/m.sheinman/Software/annovar/table_annovar.pl -buildver hg38 -out myanno -remove -protocol refGene,cytoBand,hrcr1,gnomad_genome,gnomad_exome,esp6500siv2_aa,esp6500siv2_ea,exac03nontcga,gme,kaviar_20150923,ALL.sites.2015_08,clinvar_20221231 -operation gx,r,f,f,f,f,f,f,f,f,f,f -nastring . -vcfinput -polish -outfile ",out_dir,IDs[i],".csv --xreffile /home/m.sheinman/Software/annovar/example/gene_xref.txt ",input_dir,IDs[i],".vcf.gz humandb/ & \n")
  cat(command,file = commandfile,append=TRUE)
}
cat("echo 123 \n \n",file = commandfile,append=TRUE)
cat("\n \n",file = commandfile,append=TRUE)




# import samples info----
datall <- data.frame()
library(readxl)
library(DelayedMatrixStats)
library(stringr)
sheets <- excel_sheets("/home/m.sheinman/PDX/data/WES/Samples.xlsx")
SamplesInfo <- read_excel("/home/m.sheinman/PDX/data/WES/Samples.xlsx",sheet="Oncoprint")
names <- SamplesInfo$WES; names <- names[!is.na(names)]
IDs <- SamplesInfo$`Sample ID...6`; IDs <- IDs[!is.na(IDs)]
for (i in 1:length(IDs))
{
  print(i)
  dat <- read.table(file = paste0("/home/m.sheinman/PDX/data/WES/annotation/annovar/",IDs[i],".csv.hg38_multianno.txt") ,header=TRUE, row.names=NULL,sep="\t",skip = 0,fill=NA)
  dat <- dat[!is.na(dat$Otherinfo10),]
  dat <- dat[dat$Otherinfo10=="PASS",]
  dat$TLOD <- as.numeric(sapply(1:nrow(dat),function(j){strsplit(strsplit(dat$Otherinfo11[j],"TLOD=")[[1]][2],",")[[1]][1]}))
  dat <- dat[dat$TLOD>10,]
  dat$DP <- as.numeric(sapply(1:nrow(dat),function(j){str_match(dat$Otherinfo11[j], "DP=\\s*(.*?)\\s*;")[2]}))
  dat <- dat[dat$DP>15,]
  dat$AF <- sapply(1:nrow(dat),function(j){str_match(dat$Otherinfo13[j], ":\\s*(.*?)\\s*:")[2]})
  dat$AF <- as.numeric(sapply(1:nrow(dat),function(j){strsplit(dat$AF[j],",")[[1]][2]}))/dat$DP
  dat <- dat[dat$AF>0.2,]
  dat[dat=="."] <- "0"
  dat[,c(13:54,57)] <- sapply(dat[,c(13:54,57)], as.numeric)
  dat <- dat[rowProds(dat[,c(13:54,57)]<0.01)==1,]
  dat$ID <- IDs[i]
  datall <- rbind(datall,dat)
  print(nrow(dat))
  print(table(dat$ExonicFunc.refGene))
  print(table(is.na(dat$Otherinfo10)))
}
write.table(datall,file = "/home/m.sheinman/PDX/data/WES/annotation/all.csv.hg38_multianno.txt", row.names = FALSE);
panel_genes <- as.data.frame(read_excel("/home/m.sheinman/PDX/data/WES/Gene_panel.xlsx"))[,1]
datpanel <- datall[datall$Gene.refGene %in% panel_genes,]
print(table(datpanel$ExonicFunc.refGene))
write.table(dat_panel,file = "/home/m.sheinman/PDX/data/WES/annotation/all.csv.hg38_multianno.panel.txt", row.names = FALSE);

####################################################################################################################################################################################################################################################################################

