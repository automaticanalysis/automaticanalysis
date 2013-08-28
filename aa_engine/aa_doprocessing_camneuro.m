function aa_doprocessing_camneuro(aap,bucket_foranalyzed,bucket_forrawdicom)

global cc

if (~exist('bucket_foranalyzed','var'))
    bucket_foranalyzed='imaging-analysis';
end;

if (~exist('bucket_forrawdicom','var'))
    bucket_forrawdicom='imaging-rawdata';
end;

import camneuro.aacc.*;


if (~exist('cc','var') || isempty(cc))
    cc=Aacc();
end;

cc.bucketname=bucket_foranalyzed;
cc.rawdicombucketname=bucket_forrawdicom;
  
if (isstruct(aap))
    aap_file=[tempname '.mat'];
    save(aap_file,'aap');
    deletetemp=true;
else
    aap_file=aap;
    deletetemp=false;
    load aap_file;
end;

cc.uploadAap(aap_file,aap.directory_conventions.analysisid);

if (deletetemp)
    delete(aap_file);
end;