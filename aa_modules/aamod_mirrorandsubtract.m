% AA module
% [aap,resp]=aamod_mirrorandsubtract(aap,task,i,j)
% Mirror and subtract images from one another, so that we can obtain
% hemispheric differences...
% Modified for aa by Rhodri Cusack Mar 2006
% @@@ THIS IS NOT YET TRANSFORMED TO AA4 @@@

function [aap,resp]=aamod_mirrorandsubtract(aap,task,i,j)

resp='';

switch task
    case 'domain'
        resp='session';   % this module needs to be run once per subject

    case 'description'
        resp='SPM5 mirror and subtract';

    case 'summary'
        subjpath=aas_getsubjpath(i);
        resp=sprintf('Mirror and subtract %s\n',subjpath);

    case 'report'
        
    case 'doit'

        imgs= aas_getimages(aap,i,j,'swar');
        for i=1:size(imgs,1)
            V=spm_vol(imgs(i,:));
            Y=spm_read_vols(V);
            Ym=flipdim(Y,1);
            Y=Y-Ym;
            [pth fle ext]=fileparts(V.fname);
            V.fname=fullfile(pth,['m' fle ext]);
            
            % djm: files were not writing consistently; hoping this might
            % help...
            V=rmfield(V,{'private','pinfo'});
            aas_log(aap,false,sprintf('%g of %g',i,size(imgs,1)))
            
            spm_write_vol(V,Y);
        end;
        
    case 'checkrequirements'

    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end;



