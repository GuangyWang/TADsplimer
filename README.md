# TADsplimer - a pipeline for analyzing TAD splitting and merging

Introduction
----------

This pipeline is designed for automated end-to-end analyzing TAD splitting and merging. 

If you have any questions or issues, you can ask in "Issues" of this Git or send email to Guangyu.Wang@childrens.harvard.edu or guangywang@gmail.com.

Required packages for executing TADsplimer
----------

When executing TADsplimer, user should install following packages / libraries in the system:
1. R environment (we have used 3.4.1)
2. Python version 3.5.5
3. The R packages: bcp, cluster, mass, essentials, spatstat, hicrep, changepoint, optparse
4. The python packages: numpy, cv2, PIL, scipy, imagehash

User should include the PATH of above mentioned libraries / packages inside their SYSTEM PATH variable. 

Installation
----------

The docker can be directly downloaded from dockerhub (https://hub.docker.com/r/guangywang/tadsplimer) with the following command.

	docker pull guangywang/tadsplimer:v1.0.3		 

Conda setup (Contribute by Jimin Tan, Thanks):

	conda create --prefix ./tadsplimer --file tadsplimer_pkgs.txt
	conda activate ./tadsplimer

Execution
----------
In general, TADsplimer can be executed by following command line options:

	docker run -v /<path>/:/data/ -t guangywang/tadsplimer:v1.0.3 python3 /bin/TADsplimer.py  <command>  <path> [optional arguments]		 

Please make sure output filder is empty or doesn't exit.

TADsplimer involves following command options:

split_TADs: 
	
split TAD detection using two contact maps as input files

	-h, --help            show this help message and exit
	-c, --contact_maps CONTACT_MAP
		paths to two contact maps. paths must be separated by
		the comma ','. (default: None)
	--contact_maps_aliases
                A set of short aliases for two contact maps. Aliases
		must be separated by the comma ','. (default: None)
	-u, --up_cutoff UP_CUTOFF
		paths for up cutoff of two contact maps,paths must be 
		separated by the comma ','. (default: None)
	-d, --down_cutoff DOWN_CUTOFF
		paths for down cutoff of two contact maps, paths must 
		be separated by the comma ','. (default: None)
	-j, --adjust_quality ADJUST_QUALITY
		set as 0 to auto optimize up_cutoff and down_cutoff, 
		set as 1 not auto optimize up_cutoff and down_cutoff
		(default: 0)
	-o, --output OUTPUT
		path to output files (default: None)
	-d, --split_direction DIRECTION
		set as 0: output TADs split in both two contact maps, 
		set as 1: output TADs split in contact map1, set as 2: 
		output TADs split in contact map2 (default: 0)

The following commands can be used to detect TAD split and merge events.

	docker run -v /<path>/:/data/ -t guangywang/tadsplimer:v1.0.3 python3 /bin/TADsplimer.py split_TADs -c /data/simulation_merge.txt,/data/simulation_split.txt --contact_maps_aliases merge,split -o /data/output

TAD_calculator:

topological domain identification

	-h, --help            show this help message and exit
	-c, --contact_map CONTACT_MAP
		path to Hi-C contact map (default: None)
	-u, --up UP_CUTOFF
		up cutoff for Hi-C contact map detection (default: 0)
	-d, --down DOWN_CUTOFF
		down cutoff for Hi-C contact map detection (default: 0)
	-o, --TAD_output OUTPUT
		path for the output file of TADs (default: None)
	-p, --TAD_plot PLOT
		Set to 1 to plot the contact map and TAD, else set to 
		0 to cancel this analysis. (default: 1)
	--sub_map SUB_MAP
		Set to 1 to output sub contact maps and TADs, else set 
		to 0 to cancel this analysis. (default: 1)

The following commands can be used to detect TADs.

	docker run -v /<path>/:/data/ -t guangywang/tadsplimer:v1.0.3 python3 /bin/TADsplimer.py TAD_calculator -c /data/simulation_merge.txt -u 1.7 -d 0.2 -o /data/output
	 

TAD_similarity:

calculating four similarity scores for given TADs

	-h, --help            show this help message and exit
	-c, --contact_maps CONTACT_MAP
		paths to Hi-C contact maps in two conditions. paths must 
		be separated by the comma ','. (default: None)
	-t, --TAD TAD
		input files of TADs for two compared Hi-C contact maps. 
		Paths must be separated by the comma ','. (default: None)
	-o, --output OUTPUT
		path to output files (default: None)
		 
The following commands can be used to calculate four similarity scores for given TADs.

	docker run -v /<path>/:/data/ -t guangywang/tadsplimer:v1.0.3 python3 /bin/TADsplimer.py TAD_similarity -c /data/simulation_merge.txt -t /data/tad.txt -o /data/output

Input
----------

The input are two Hi-C matrices (full TSV matrices) to be compared.  The Hi-C matrices 
should have the dimension N * N. 10Kb is the default resolution for the input matrices.

Output
----------

The main outputs are a split TAD file and a merged TAD file.
Columns and explaination are as follwing:

### Split TAD file
Files of *.all.split.txt
column | explaination
------ | ------------
1th | Row number
2th | Start position of split TAD  
3th | End position of split TAD

### Merged TAD file
Files of *.all.merge.txt
column | explaination
------ | ------------
1th | P value for TAD matching
2th | Start position of merged TAD  
3th | End position of merged TAD
4th | Laplacian matrix similarity score
5th | Corner split ratio
6th | Stratum-adjusted correlation coeffient
7th | Image hashing similarity score

In the output files, the "A->B.all.merge.txt" is the coordinate of merged TADs comparing A sample to B sample. The file "A->B.all.split.txt" is the coordinate of split TADs comparing A sample to B sample.

## Update
08/17/2020
1.	remove python dependency rpy2
2.	fix bugs for small size input file

Simulation
----------
In general, simulation of TADs can be executed by following command line options:

	Rscript simulation.R -i <path of reference TAD> -o <path of output>		 

Reference TADs can be downloaded from ./reference folder.

This source code is released under an open source licence compliant with MIT license which is approved by Open Source Licenses (OSI) (https://opensource.org/licenses ).


How to Cite?
----------
Please cite the following required publication.

Wang, G., Meng, Q., Xia, B. et al. TADsplimer reveals splits and mergers of topologically associating domains for epigenetic regulation of transcription. Genome Biol 21, 84 (2020). https://doi.org/10.1186/s13059-020-01992-7
