# Copy-Number-Variant-Finder-CLI
A command line interface version of the Copy-Number-Variant-Finder tool

GUI version is available at:  
_________________________________________________________________________________________________________________________________________________

Use at your own risk.

I cannot provide support. All information obtained/inferred with this script is without any implied warranty of fitness for any purpose or use whatsoever.

ABOUT: 

This program identifies putative copy number variants from BAM files and plots the results in both scatter and tile format to facilitate the identification of rare events such as those induced when treating cells with ionizing radiation. User selected parameters include sample ploidy, the bin size (the number of consecutive non-overlapping base pairs from which to derive mean coverage for a region), choice of chromosomes/contigs to evaluate, and the choice of a control sample to compare data to. 

PREREQUISITES:

    1. A set of BAM files that have had PCR duplicates removed. 
    2. The following tools installed: samtools, zenity, awk, tr, datamash, R, and the R package ggplot. This program was built to run on Ubuntu 20.04 and higher and has been tested on macOS Mojave. See installation notes about running on other systems. 
    
INSTALLATION:

Note: A GUI version of this tool exists. You may want to try it first. <URL>

Linux/Ubuntu:
Most tools that this program requires can be installed in the Linux command line by typing the name of the tool. Either version information (if already installed) or installation instructions (if not installed) will appear in the terminal. Follow the installation instructions to install the tool. To install ggplot2 first launch R by typing R in a terminal window. Next, type install.packages("ggplot2") to install the package. Next, downlaod the .sh file from this page and provide it permission to execute using chmod +x .

macOS:  Install homebrew from the terminal window. Next, install other tools using brew install from the terminal (for example brew install samtools). The tools are: samtools, datamash, and R.  If running an older version of macOS (e.g. Mojave), you may need to install R by visiting cran.r-project.org to select and install the correct version.  Once installed, open a terminal window and launch R by typing R.  Install ggplot2 by typing install.packages("ggplot2") from the R terminal. 

Windows: 
NOT TESTED. In theory you can install Linux bash shell on Windows (https://itsfoss.com/install-bash-on-windows/) and install the dependencies from the command line. If you try this and it works, please let me know. I don't have a Windows machine for testing.

EXAMPLE DATA:

Example data can be found in https://www.ncbi.nlm.nih.gov/sra/?term=SRP130725. Using SRA toolkit download the four “Non treated Grande Naine” samples (e.g. prefetch SRX3579227) and the four “Gamma irradiated Novaria” samples.  These represent four “wild-type” replicates and four replicates of the mutant banana variety Novaria.  A large deletion in chromosome 5 was previously reported in Novaria by evaluating coverage variation from BAM files (https://pubmed.ncbi.nlm.nih.gov/29476650/).  Convert the downloaded .sra files into .sam files using sra-dump (e.g. sam-dump SRR6489404.sra > N1.sam, making sure to correctly name the files).  This banana data was originally deposited as BAM data that had PCR duplicates removed.  Therefore the .sam files generated can be directly used to test this program.  When using the test data, select 3 for ploidy, a bin size of 100000 and GN1 as the control.  Test the chromosome selection feature by selecting chromosome 5.  

TO RUN:

Collect the command line argument information needed prior to running the program. It may be helpful to create a text file that contains all the information you need. If this seems like a lot of work, try running the GUI version, where you don't need to do any of this. The following information is needed:
a: A file that lists the paths of all the bam files. It should be a plain text file and look like:
/home/brad/Documents/43803_rg_dedup.bam
/home/brad/Documents/43802_rg_dedup.bam
/home/brad/Documents/43801_rg_dedup.bam
/home/brad/Documents/43804_rg_dedup.bam 

It is easiest to put all the bam files in the same directory.  Next, open a terminal window and drag a bam file into the window.  The path should appear.  Copy this into the text file and repeat for all bams.  The text file containing BAM names and paths should be in the same directory as this program.  

b: The ploidy of the samples being evaluated. The ploidy parameter is used to convert data for plotting so that it is easier to identify single and multi-copy number events (e.g. if ploidy is set to 2, a single copy deletion would appear near 1 on the y-axis of plots). 
	
c: The bin size for your analysis.  This is a number that defines the non-overlapping window used to calculate the average coverage for that region.  It should not contain any symblols.  For example: 100000

Bin size is the size of consecutive non-overlapping base pairs for coverage analysis.  For example, if the bin size is set to 100000, then the average coverage every 100000 base pairs is calculated and subsequently plotted.  Larger bin sizes are better for lower coverage data.  For example, the banana data described above has between 2 and 3.5x coverage per sample and a large ~3.8 Mbp deletion was discovered on chromosome 5 when using a bin size of 100000.  Lowering the bin size may enable recovery of smaller indels.  

d: The name of the sample to be used as a control.  Note that the name should match the exact name of the bam file.  For example from section a above: 43801_rg_dedup.bam (and not 43801).  If there is no control, type none in all lower case.
	
	
Once you have collected all the information you are ready to run the program.  From a terminal window give permission for the program to run using chmod +x. The command structure to run the program is	
./CLITest_CNVT_1b.sh -a mybamfiles.txt -b 3 -c 100000 -d GN1.bam

Launching the program without arguments will result in a text output that explains the variables required. 
 
OUTPUT:  

A new directory is created named CNVF_Analysis that contains the results in sub-directories. When no control sample is selected, the program creates a directory named MeanCompare where the coverage of every sample is compared to the average, or mean, of all samples for each bin of base pairs.  Two plots are generated for each selected chromosome: one ending in _Coverage.jpeg and the other in _CoverageGroups.jpeg.  The former is a scatter plot of coverage values relative to the mean of all samples.  The latter is the same data in a tile plot where specific coverage ranges are grouped and color coded for ease of interpretation.  Groups range from 0 to 6, with 3 indicating coverage near the selected ploidy (no CNV observed).  Group 0 indicates an approximate 3x or lower reduction in copy, 1 indicates an approximate 2x reduction in copy, 2 indicates an approximate 1x reduction in copy, 4 indicates a ~1x increase in copy number, 5 a ~2x increase and 6 a ~3x or higher increase.  In addition to the plots, a .csv file is provided containing the data used to make the plots.  This can be used to automate the identification of candidate variants.  

When a control sample is selected, the mean comparison data is produced along with data compared to the selected control sample.  This data appears in a directory titled “ControlCompare”.  It is similar to data in the MeanCompare directory except that data is plotted relative to the values of the control sample.  As such, the control sample appears as a solid line on the scatter plot or as a solid green bar on the tile plot.  
