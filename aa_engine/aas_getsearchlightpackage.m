function [searchlights]=aas_getsearchlightpackage(aap,packagenumber)
packagesize=aap.options.searchlight.N/aap.options.searchlight.Npackage;
searchlights=max(1,ceil(((packagenumber-1)*packagesize+1))):min(floor(((packagenumber)*packagesize+1)),aap.options.searchlight.N);