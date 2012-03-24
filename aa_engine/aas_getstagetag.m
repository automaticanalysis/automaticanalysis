% Automatic Analysis
%  Stages are processing stages (e.g., aamod_realign)
%  From aa_ver_3 onwards, it is possible to have multiple instances of a
%  single module in a task list - so realign could occur twice in theory
%  The done flags and dependencies for the mulitple instances are
%  differentiated by the addtion of _0002 _0003 etc
%  The stage tag is the stage name, with this index (e.g.,
%  aamod_realign_0002)
%
%  RC 15/2/2010
%  From version 4.0, even when only one stage, index is added, as this way
%  if a later stage is added that repeats one already present, previous one
%  will be identified
%
%  RC 15/2/2010
%  Changed to 5 digit tag to allow for highly branched designs
%
% function [stagetag]=aas_getstagetag(aap,stage)
%  stage = stage number

function [stagetag]=aas_getstagetag(aap,stage)

% allow full path of module to be provided
[stagepath stagename]=fileparts(aap.tasklist.main.module(stage).name);

 stagetag=sprintf('%s_%05d',stagename, aap.tasklist.main.module(stage).index);

