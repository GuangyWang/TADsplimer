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

TADsplimer involves following command line options:

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
		 
TAD_calculator: topological domain identification

	docker run -v /<path>/:/data/ -t tadsplimer:v1 python3 /bin/TADsplimer.py  <command>  <path> [optional arguments]		 

split_TADs_alternate: split TAD detection using two contact maps and detected TADs for these two contact maps as input files

	docker run -v /<path>/:/data/ -t tadsplimer:v1 python3 /bin/TADsplimer.py  <command>  <path> [optional arguments]		 

TAD_similarity: calculating four similarity scores for given regions

	docker run -v /<path>/:/data/ -t tadsplimer:v1 python3 /bin/TADsplimer.py  <command>  <path> [optional arguments]		 



