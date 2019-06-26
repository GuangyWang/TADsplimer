# TADsplimer - a pipeline for analysising TAD splitting and merging

Introduction
----------

This pipeline is designed for automated end-to-end analysising TAD splitting and merging. 

Required packages for executing TADsplimer
----------

When executing TADsplimer, user should install following packages / libraries in the system:
1. R environment (we have used 3.4.1)
2. Python version 3.5.5
3. The R packages: bcp, cluster, mass, essentials, spatstat, hicrep, changepoint
4. The python packages: numpy, rpy2, cv2, PIL, scipy, imagehash

User should include the PATH of above mentioned libraries / packages inside their SYSTEM PATH variable. 

Installation
----------

The docker can be directly downloaded from dockerhub (https://hub.docker.com/r/guangywang/tadsplimer) with the following command.

	docker push guangywang/tadsplimer		 


Execution
----------
In general, TADsplimer can be executed by following command line options:

	docker run -v /<path>/:/data/ -t tadsplimer:v1 python3 /bin/TADsplimer.py  <command>  <path> [optional arguments]		 

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

	docker run -v /<path>/:/data/ -t tadsplimer:v1 python3 /bin/TADsplimer.py split_TADs -c /data/simulation_merge.txt,/data/simulation_split.txt --contact_maps_aliases merge,split -o /data/output

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

	docker run -v /<path>/:/data/ -t tadsplimer:v1 python3 /bin/TADsplimer.py TAD_calculator -c /data/simulation_merge.txt -u 1.7 -d 0.2 -o /data/output
	 

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

	docker run -v /<path>/:/data/ -t tadsplimer:v1 python3 /bin/TADsplimer.py TAD_similarity -c /data/simulation_merge.txt -t /data/tad.txt -o /data/output


