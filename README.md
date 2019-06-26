# TADsplimer - a pipeline for analysising TAD splitting and merging

Introduction
----------

This pipeline is designed for automated end-to-end analysising TAD splitting and merging. 

Required packages for executing TADsplimer
----------

When executing TADsplimer, user should install following packages / libraries in the system:
1. R environment (we have used 4.1.11)
2. Python version 3.5
3. The R packages: bcp, cluster, mass, essentials, spatstat, hicrep, changepoint
4. The python packages: numpy, rpy2, cv2, PIL, scipy, imagehash

User should include the PATH of above mentioned libraries / packages inside their SYSTEM PATH variable. 

Installation
----------

The docker can be directly downloaded from dockerhub(https://hub.docker.com/r/guangywang/tadsplimer) with the command.

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
  	--up_cutoff
		paths for up cutoff of two contact maps,paths must be 
		separated by the comma ','. (default: None)
  	--down_cutoff
		paths for down cutoff of two contact maps, paths must 
		be separated by the comma ','. (default: None)
  	-j, --adjust_quality ADJUST_QUALITY
		set as 1 to normalize sequence quality for two Hi-C 
		contact maps, set as 0 not to normalize sequence quality 
		for two Hi-C contact maps (default: 0)
  	-o, --output OUTPUT
		path to output files (default: None)
  	-d, --split_direction DIRECTION
		set as 0: output TADs split in both two contact maps, 
		set as 1: output TADs split in contact map1, set as 2: 
		output TADs split in contact map2 (default: 0)
		 
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
	-p PLOT, --TAD_plot PLOT
		Set to 1 to plot the contact map and TAD, else set to 
		0 to cancel this analysis. (default: 1)
	--sub_map SUB_MAP
		Set to 1 to output sub contact maps and TADs, else set 
		to 0 to cancel this analysis. (default: 1)
		 

split_TADs_alternate:

split TAD detection using two contact maps and detected TADs for these two contact maps as input files

	  -h, --help            show this help message and exit
	  -c, --contact_maps CONTACT_MAP
		paths to Hi-C contact maps in two conditions. paths must 
		be separated by the comma ','. (default: None)
	  --contact_maps_aliases ALIASES
		A set of short aliases for the contact map. Paths must be 
		separated by the comma ','. (default: None)
	  -t, --TAD TAD
		input files of TADs for two compared Hi-C contact maps. 
		Paths must be separated by the comma ','. (default: None)
	  -u, --up_cutoff UP
	  	up cutoff for two compared Hi-C contact maps, paths
		must be separated by the comma ','. (default: 0,0)
	  -j, --adjust_quality ADJUST_QUALITY
	  	set as 1 to normalize sequence quality for two Hi-C contact 
		maps, set as 0 not to normalize sequence quality for two 
		Hi-C contact maps (default: 0)
	  -o, --output OUTPUT
	  	path to output files (default: None)
	  -d, --split_direction DIRECTION
	  	set as 0: output TADs split in both two contact maps,
		set as 1: output TADs split in contact map 1, set as 2: 
		output TADs split in contact map 2 (default: 0)
		 
TAD_similarity:

calculating four similarity scores for given regions

	docker run -v /<path>/:/data/ -t tadsplimer:v1 python3 /bin/TADsplimer.py  <command>  <path> [optional arguments]		 



