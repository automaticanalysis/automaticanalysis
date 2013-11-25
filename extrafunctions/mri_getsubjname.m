% wrapper for aamod_evaluatesubjectnames
% usefull in pre-initialisation phase when subjectsnames are not evaluated yet
%
% Tibor Auer MRC CBU Cambridge 2012-2013

function strSubj = mri_getsubjname(aap,subj)

aap = aamod_evaluatesubjectnames(aap,'doit',subj);
strSubj = basename(aas_getsubjpath(aap,subj));