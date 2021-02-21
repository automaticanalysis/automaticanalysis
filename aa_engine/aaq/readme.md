# Queue processors in automaticanalysis
-- Work in progress --
### aaq.m
Base class of all queue processors; cannot be chosen as queue processor.

###aaq_localsingle.m
Simplest queue processor: runs all tasks serially in the local Matlab session, that is, without explicit parallelisation. Useful for development, debugging, and estimating baseline performance.

### aaq_qsub.m
Queue processor running independent __batch__ jobs on a parallel cluster initiated with __parcluster__. 

Very robust queue processor with safety nets like repeated submission of jobs to catch hickups in the cluster.
The correct order of execution of modules is ensured via 'job done' flags, namely, files written into each module's working directory.

### aaq_parpool.m
aaq_parpool starts and reserves a pool of Matlab workers with __parpool__ prior to any module computation. Parallelisation is achieved via __parfeval__. 

The advantages compared to other queue processors are in terms of the availability of workers: 
i. Short latency (because Matlab processes don't need to be started up for each task)
ii. Workers are reserved, so their availability is guaranteed 

Compared to other queue processors with parallelisation, aaq_parpool features an alternative, faster computation of job dependencies, which rather than relying on 'job done' flags in the form of files written into each module's working directory, keeps track of jobs via variables within the client session's workspace.

aaq_parpool lends itself to analyses with many tasks of little complexity, and/or computations on individual multicore machines, or other environments in which blocking CPUs is not an issue.

aaq_parpool has been developed from aaq_matlab_pct (see below). The streamcache/real-time code in aaq_matlab_pct has not been ported.

### aaq_matlab_pct.m
aaq_matlab_pct starts and reserves a pool of Matlab workers with __parpool__ prior to any module computation. Parallelisation is achieved via __spmd__. 
aaq_matlab_pct is outdated and deprecated; its designated successor is aaq_parpool.

Scenarios for usage are the same as for aaq_parpool.

### aaq_qsub_debug.m, aaq_qsub_monitor_jobs.m
Outdated variants of aaq_qsub using createJob/createTask, and intended for debugging. 

### aaq_condor.m
-- to come --
