function copy_globals_to_workers(aacache_from_client)
% COPY_GLOBALS_TO_WORKERS - assign input args to global variables.
%
% copy_globals_to_workers(aacache_from_client) assigns input arguments to
% global variable aacache. Currently, only input argument
% aacache_from_client is assigned to aacache. This seemingly nonsensical
% action is necessary for aa to work with the spmd-based explicit
% parallelism as implemented in the matlab_pct queuer (aaq_matlab_pct). The
% reason is that workers instantiated by the spmd construct each possess
% versions of global variables (note the oxymoron) distinct from that of
% the client. These are initially empty, and will not be automatically
% synchronized among each other or with the client.
%
% copy_globals_to_workers must be called by parfevalOnAll from the client
% before the spmd ... end block, with aacache_from_client corresponding to
% global variable aacache. 
%
% The whole idea will fail if any of the workers changes aacache AND code
% executed by any of the other workers within the same spmd block relies on
% such changes, or if any of the changes need to be present also in the
% global variable known to the client. As of Dec 28, 2020, it appears that
% aacache is only read, but not changed, within the spmd code block, so the
% approach seems safe.
% 
%  HH Dec 2020

global aacache
aacache = aacache_from_client;

