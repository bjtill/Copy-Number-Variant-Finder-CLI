#!/bin/bash
#BT February 24, 2023
#Command line version

##################################################################################################################################

# Help
[ "$1" = "-h" -o "$1" = "--help" ] && echo "
Copy Number Variation Finder (CNVF) command line version 1.5
ABOUT: 
This program aids in the identification of copy number variants between control and test samples by parsing coverage data from a BAM (or SAM) alignment file and plotting the results.  

PREREQUISITES:
1) BAM files that are position sorted and have PCR duplicates removed, 2) An index of each BAM file, 3)the following tools installed: samtools, bash, datamash, awk, sed, zenity, R, and the R package ggplot2. This program was built to run on Ubuntu 20.04 and higher. See the readme file for information on using with other opperating systems.  

TO RUN:
1) Provide permission for this program to run on your computer (open a terminal window and type chmod +x Copy_Number_Variant_Hunter_CLI_V1_5.sh).  Check to make sure that the name is an exact match to the .sh file you are using as the version may change.

2) Test to see if this works by launching the program without any arguments (./CNV_Finder_CLI_V1_3.sh). You should should see some information on what you need to do to run the program.  

3) Run the program with the arguments.  Prior to doing this, you should collect the information needed.  It may be helpful to create a text file that contains all the information you need. If this seems like a lot of work, try running the GUI version, where you don't need to do any of this.  The following information is needed: 

	a: A file that lists the paths of all the bam files. It should be a plain text file and look like:
	
/home/brad/Documents/43803_rg_dedup.bam
/home/brad/Documents/43802_rg_dedup.bam
/home/brad/Documents/43801_rg_dedup.bam
/home/brad/Documents/43804_rg_dedup.bam 

	It is easiest to put all the bam files in the same directory.  Next, open a terminal window and drag a bam file into the window.  The path should appear.  Copy this into the text file and repeat for all bams.  The text file containing BAM names and paths should be in the same directory as this program.  

	b: The ploidy of the samples being evaluated. 
	
	c: The bin size for your analysis.  This is a number that defines the non-overlapping window used to calculate the average coverage for that region.  It should not contain any symblols.  For example: 100000
	
	
	d: The name of the sample to be used as a control.  Note that the name should match the exact name of the bam file.  For example from section a above: 43801_rg_dedup.bam (and not 43801).  If there is no control, type none in all lower case.
	
4) Once you have collected all the information in section 3, you are ready to run the program.  The command structure is	
	
LICENSE:  
MIT License 

Copyright (c) 2023 Bradley John Till

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the *Software*), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED *AS IS*, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
Version Information:  Version 1.5, February 24, 2023
" && exit

############################################################################################################################

helpFunction()
{
   echo ""
   echo "Usage: $0 -a BAM list -b Sample ploidy -c Bin size -d Control sample"
   echo ""
   echo -e "\t-a A file containing the paths to the selected BAM file.  
   	For example /home/brad/Documents/43801_rg_dedup.bam \n"
   echo -e "\t-b The ploidy of your samples (e.g. 2)"
   echo -e "\t-c The bin size for your analysis. Symbols not allowed (e.g. 100000 NOT 100,000)"
   echo -e "\t-d The name of the BAM file of your control sample (e.g. 43801_rg_dedup.bam).  If no control sample, type 
   	none"
   echo ""
   echo "For more detailed help and to view the license, type $0 -h"
   exit 1 # Exit script after printing help
}

while getopts "a:b:c:d:e:" opt
do
   case "$opt" in
      a ) parameterA="$OPTARG" ;;
      b ) parameterB="$OPTARG" ;;
      c ) parameterC="$OPTARG" ;;
      d ) parameterD="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Print helpFunction in case parameters are empty
if [ -z "$parameterA" ] || [ -z "$parameterB" ] || [ -z "$parameterC" ] || [ -z "$parameterD" ] 
then
   echo "Some or all of the parameters are empty";
   helpFunction
fi
echo $(head -1000 "$parameterA") > bp1
tr ' ' '\t' < bp1 | datamash transpose > bp3
mkdir CNVF_Analysis
mv bp3 ./CNVF_Analysis/
rm bp1
cd CNVF_Analysis
# Store parameters for later use

echo "$parameterB" > ploidy
echo "$parameterC" > binsize
echo "$parameterD" > answer.control


############################################################################################################################




log=CNVFt.log
now=$(date)  
echo "Copy Number Variant Finder (CNVF) command line version 1.0
-h for help
Script Started $now."  > $log

OR='\033[0;33m'
NC='\033[0m'
printf "${OR}Copy Number Variant Finder (CNVF) command line version 1.0
-h for help
Script Started $now.${NC}\n" 

############################################################################################################################
##Determine if there is a control sample (used later)
for file in *.control
do
  
   if [ -f "$file" ]
   then
       newname=`head -1 $file`
       if [ -f "$newname" ]
       then
              echo "Cannot rename $file to $newname - file already exists"
       else
              mv "$file" "$newname".controlsample
       fi
   fi
done

#Get coverage 
echo "Generating coverage data from the BAM files. This may take a long time." >> $log
printf "${OR}Generating coverage data from the BAM files. This may take a long time.${NC}\n"
samtools depth -a -H -f bp3 -o STdepth2
tail -n +2 STdepth2 > OT
head -1 STdepth2 > h2
rm STdepth2 
#Split the table into info columns and data colums for processing
awk '{print $1, $2}' OT > OTA
cut -f3- OT > OTB
#Get info columns with binsize
a=$(head -1 binsize)
awk -v var=$a 'NR % var == 0' OTA > OT1
awk 'BEGIN{print "Chromosome", "Position"}1' OT1 > IC2  #Info Columns with header
rm OT1
#Take the mean of every n lines with binsize variable 
a=$(head -1 binsize)
awk -v N=$a '{ for (i = 1; i <= NF; i++) sum[i] += $i } NR % N == 0 { for (i = 1; i <= NF; i++) {printf("%.6f%s", sum[i]/N, (i == NF) ? "\n" : " ") 
sum[i] = 0}}' OTB > OTC  #NOTE the newline here is intentional 
cut -f3- h2 > h3
rm OTB
awk -F'/' '{print $NF}' bp3 > bp4
datamash transpose < bp4 > bp4t
###################################################################################################
echo "Processing coverage data." >> $log
printf "${OR}Processing coverage data.${NC}\n"


[[ -f none.controlsample ]] ||
echo "yescontrol" > YES_control 
if [ -f "YES_control" ]  
then 
b=$(head -1 *.controlsample)
awk -v b="$b" '{for (i=1;i<=NF;i++) { if ($i == b) { print i } }}' bp4t > h4
c=$(head -1 h4)
awk -v c="$c" '{ for( i=1;i<=NF;i++ ) { printf "%f%s",  $i/($c+0.0000000000001), OFS }; printf "%s", ORS  }' OTC > OTG
p=$(head -1 ploidy)
awk -v p="$p" '{ for( i=1;i<=NF;i++ ) { printf "%f%s",  $i*p, OFS }; printf "%s", ORS  }' OTG > OTH
rm OTG
datamash transpose < bp4 > bp4t
awk '{for (i=1;i<=NF;i++) printf $i "_controlcompare "}' bp4t | awk '{print $0}' > cc1
cat cc1 OTH > OTI #add header
rm OTH
mkdir ControlCompare
mv OTI ./ControlCompare/
cp IC2 ./ControlCompare/
cp ploidy ./ControlCompare/
cd ControlCompare
#split each sample column to process for plotting
awk '$1=$1' OTI | awk -F '[ ]' '{for(i=1; i<=NF; i++)  print $i >> ("column" i ".txt"); close("column" i ".txt")}'
#Combine data and sample name and fill the empty columns
for i in *.txt; do 
tail -n +2 $i > ${i%.*}.tail
head -1 $i > ${i%.*}.hed
paste ${i%.*}.tail ${i%.*}.hed > ${i%.*}.sail
awk '$2==""{$2=p}{p=$2}1' ${i%.*}.sail > ${i%.*}.bail ; done 
tail -n +2 IC2 > ICUP
for i in *.bail; do 
paste ICUP ${i%.*}.bail > ${i%.*}.boil; done 
cat *.boil > controldata
sed 's/_controlcompare //g' controldata > controldata2
#then split by chromosome
awk '{print > ($1".bb")}' controldata2
#Split by chromosome
for i in *.bb; do 
#Prepare coverage groups
a=$(head -1 ploidy)
 awk -v v="$a" '{if($3 <= (v+0.7) && $3 >= (v-0.7)) print $0, "3"; else if ($3 > (v+0.7) && $3 <= (v+1.7)) print $0, "4"; else if ($3 > (v+1.7) && $3 <= (v+2.7)) print $0, "5"; else if ($3 > (v+2.7)) print $0, "6"; else if ($3 < (v-0.7) && $3 >= (v-1.7)) print $0, "2"; else if ($3 < (v-1.7) && $3 >= (v-2.7)) print $0, "1"; else if ($3 < (v-2.7)) print $3, "0"; else print $0, "ERROR"}' $i > ${i%.*}.nf3
 
 done 
 
#Add header and prepare for plotting
for i in *.nf3; do 
awk 'BEGIN{print "Chromosome", "Position", "CovBMean", "Sample", "CovGp"}1' $i > ${i%.*}.nf4; done 
for i in *.nf4; do 
tr ' ' ',' < $i > ${i%.*}.nf5
tr '\t' ',' < ${i%.*}.nf5 > ${i%.*}.nf6
done

printf 'library(ggplot2) \nfile_list=list.files(full.names=F, pattern="\\\.nf6") \nfilenames <- gsub("\\\.nf6$","", list.files(pattern="\\\_tmp$")) \nfor (i in file_list){ \ng<-read.csv(i) \np <- ggplot(g, aes(x=Position, y=Sample, fill=factor(CovGp))) + geom_tile() + scale_fill_manual(values =c("0"="firebrick1","1"="firebrick", "2"="dodgerblue4", "3"="darkolivegreen", "4"="darkblue", "5"="darkorange2", "6"="darkorange"), name = "Coverage Groups") + theme(axis.text.x = element_text(angle = 90, size =8, vjust = 0.5, hjust=1), axis.text.y = element_text(size = 8)) + xlab("Position")\np2 <- p + labs(title= sub("\\\.nf6$","",i)) \nggsave(plot = p2, filename= paste0(i, "AFbins.jpeg")) \n#p3 <- ggplotly(p2) \n#htmlwidgets::saveWidget(p3, file = paste0(i, ".html")) \n}' > AFbin.R

Rscript AFbin.R
for i in *nf6AFbins*; do mv $i ${i%.nf6*}_ControlCompare_Coverage.jpeg; done
mkdir CovPlot
cp *.nf6 ./CovPlot
cd CovPlot 

printf 'library(ggplot2) \nfile_list=list.files(full.names=F, pattern="\\\\.nf6") \nfilenames <- gsub("\\\.nf6$","", list.files(pattern="\\\_tmp$")) \nfor (i in file_list){ \ng<-read.csv(i) \np <- ggplot(g, aes(x=Position, y=CovBMean, color=Sample)) + geom_point(size=1, alpha = 0.5) + theme(axis.text.x = element_text(angle = 90, size =8), axis.text.y = element_text(size = 8)) + labs(title="_Coverage_Variation") + ylab("Coverage Compared to Mean of All Samples") \np2 <- p + labs(title= sub("\\\.nf6$","",i))\nggsave(plot = p2, filename= paste0(i, "Covbins.jpeg"), width=10, height=5, units=c("in"))  \n}' > BFbin.R

Rscript BFbin.R
for i in *nf6Covbins*; do mv $i ${i%.nf6*}_ControlCompare_CoverageGroups.jpeg; done
cp *.jpeg ..
cd ..
rm -r CovPlot
#Prepare a data table for keeping 
cat *.nf3 | awk 'BEGIN{print "Chromosome", "Position", "CovBMean", "Sample", "CovGp"}1'| tr ' ' ',' | tr '\t' ',' > AllData_ControlCompare.csv
rm *.sail *.hed *.boil *.bail *.tail *.txt controldata controldata2 IC2 ICUP OTI ploidy *.bb *.nf3 *.nf4 *.nf5 *.nf6 AFbin.R
cd ..
fi
#done

echo "Generating plots." >> $log
printf "${OR}Generating plots.${NC}\n"


#Process the data compared to the mean of all samples  
awk '{sum=0; for (i=1;i<=NF;i++)sum+=$i; print $0,sum/(NF)}' OTC > CT1
awk 'NR==1 {print NF}' CT1 > meancol
d=$(head -1 meancol)
awk -v d="$d" '{ for( i=1;i<=NF;i++ ) { printf "%f%s",  $i/($d+0.0000000000001), OFS }; printf "%s", ORS  }' CT1 > CT2   
#remove terminal column
awk 'NF{NF-=1};1' CT2 > CT3
#Setting values to user selected ploidy. 
p=$(head -1 ploidy)
awk -v p="$p" '{ for( i=1;i<=NF;i++ ) { printf "%f%s",  $i*p, OFS }; printf "%s", ORS  }' CT3 > CT4
datamash transpose < bp4 > bp4t
awk '{for (i=1;i<=NF;i++) printf $i "_meancompare "}' bp4t | awk '{print $0}' > mc1
cat mc1 CT4 > ct5
mkdir MeanCompare
mv ct5 ./MeanCompare/
cp IC2 ./MeanCompare/
cp ploidy ./MeanCompare/
cd MeanCompare

#Split tables by sample, add sample name and fill empty rows in table. 
awk '$1=$1' ct5 | awk -F '[ ]' '{for(i=1; i<=NF; i++)  print $i >> ("column" i ".txt"); close("column" i ".txt")}'
for i in *.txt; do 
tail -n +2 $i > ${i%.*}.tail
head -1 $i > ${i%.*}.hed
paste ${i%.*}.tail ${i%.*}.hed > ${i%.*}.sail
awk '$2==""{$2=p}{p=$2}1' ${i%.*}.sail > ${i%.*}.bail ; done 
tail -n +2 IC2 > ICUP
for i in *.bail; do 
paste ICUP ${i%.*}.bail > ${i%.*}.boil; done 
cat *.boil > meandata
sed 's/_meancompare//g' meandata > meandata2
#Split data by chromosome
awk '{print > ($1".bb")}' meandata2

 
for i in *.bb; do 
#Prepare coverage groups
a=$(head -1 ploidy)
 awk -v v="$a" '{if($3 <= (v+0.7) && $3 >= (v-0.7)) print $0, "3"; else if ($3 > (v+0.7) && $3 <= (v+1.7)) print $0, "4"; else if ($3 > (v+1.7) && $3 <= (v+2.7)) print $0, "5"; else if ($3 > (v+2.7)) print $0, "6"; else if ($3 < (v-0.7) && $3 >= (v-1.7)) print $0, "2"; else if ($3 < (v-1.7) && $3 >= (v-2.7)) print $0, "1"; else if ($3 < (v-2.7)) print $3, "0"; else print $0, "ERROR"}' $i > ${i%.*}.nf3
 
 done 
 
#Format for plotting
for i in *.nf3; do 
awk 'BEGIN{print "Chromosome", "Position", "CovBMean", "Sample", "CovGp"}1' $i > ${i%.*}.nf4; done 
for i in *.nf4; do 
tr ' ' ',' < $i > ${i%.*}.nf5
tr '\t' ',' < ${i%.*}.nf5 > ${i%.*}.nf6
done
 
printf 'library(ggplot2) \nfile_list=list.files(full.names=F, pattern="\\\.nf6") \nfilenames <- gsub("\\\.nf6$","", list.files(pattern="\\\_tmp$")) \nfor (i in file_list){ \ng<-read.csv(i) \np <- ggplot(g, aes(x=Position, y=Sample, fill=factor(CovGp))) + geom_tile() + scale_fill_manual(values =c("0"="firebrick1","1"="firebrick", "2"="dodgerblue4", "3"="darkolivegreen", "4"="darkblue", "5"="darkorange2", "6"="darkorange"), name = "Coverage Groups") + theme(axis.text.x = element_text(angle = 90, size =8, vjust = 0.5, hjust=1), axis.text.y = element_text(size = 8)) + xlab("Position")\np2 <- p + labs(title= sub("\\\.nf6$","",i)) \nggsave(plot = p2, filename= paste0(i, "AFbins.jpeg")) \n#p3 <- ggplotly(p2) \n#htmlwidgets::saveWidget(p3, file = paste0(i, ".html")) \n}' > AFbin.R

Rscript AFbin.R
#rename plots 


mkdir CovPlot
cp *.nf6 ./CovPlot
cd CovPlot 

printf 'library(ggplot2) \nfile_list=list.files(full.names=F, pattern="\\\\.nf6") \nfilenames <- gsub("\\\.nf6$","", list.files(pattern="\\\_tmp$")) \nfor (i in file_list){ \ng<-read.csv(i) \np <- ggplot(g, aes(x=Position, y=CovBMean, color=Sample)) + geom_point(size=1, alpha = 0.5) + theme(axis.text.x = element_text(angle = 90, size =8), axis.text.y = element_text(size = 8)) + labs(title="_Coverage_Variation") + ylab("Coverage Compared to Mean of All Samples") \np2 <- p + labs(title= sub("\\\.nf6$","",i))\nggsave(plot = p2, filename= paste0(i, "Covbins.jpeg"), width=10, height=5, units=c("in"))  \n}' > BFbin.R

Rscript BFbin.R
cp *.jpeg ..
cd .. 
rm -r CovPlot
for i in *nf6AFbins*; do mv $i ${i%.nf6*}_MeanCompare_CoverageGroups.jpeg; done
#Prepare a data table for keeping
cat *.nf3 | awk 'BEGIN{print "Chromosome", "Position", "CovBMean", "Sample", "CovGp"}1'| tr ' ' ',' | tr '\t' ',' > AllData_MeanCompare.csv
rm *.sail *.hed *.boil *.bail *.tail *.txt IC2 ICUP ploidy *.bb *.nf3 *.nf4 *.nf5 *.nf6 AFbin.R ct5 meandata meandata2
for i in *nf6Covbins*; do mv $i ${i%.nf6*}_MeanCompare_Coverage.jpeg; done
cd ..

echo "Final processing steps. Program almost finished." >> $log
printf "${OR}Final processing steps. Program almost finished.${NC}\n"



now=$(date)  
printf "${OR}Script Finished $now. The logfile is named CNVF.log.${NC}\n" 
echo "Script Finished $now." >> log
#Collect info on user parameters
awk 'BEGIN{print "BAM files selected and their path:"}1' bp3 > B4L
awk 'BEGIN{print "Bin size selected (in base pairs):"}1' binsize > BSL
 #awk 'BEGIN{print "Chromosomes selected:"}1' CHRL > CSL
awk 'BEGIN{print "Ploidy selected:"}1' ploidy > PSL

cat CNVFt.log B4L BSL PSL > CNVF.log

rm binsize bp3 bp4 bp4t cc1 CNVFt.log CT1 CT2 CT3 CT4 h2 h3 h4 IC2 mc1 meancol OT OTA OTC ploidy B4L BSL PSL *.controlsample *_control log

#End of program
