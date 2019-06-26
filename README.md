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
In general, TADsplimer involves following command line options:

	docker run -v /<path>/:/data/ -t tadsplimer:v1 python3 /bin/TADsplimer.py  <command>  <path> [optional arguments]		 





