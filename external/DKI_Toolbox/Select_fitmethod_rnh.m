function [fitmethod, sel]=Select_fitmethod_rnh(DWI_name,Mask_name,bval_name,bvec_name)
% Rafael Neto Henriques
% 2 July 2013 IBEB
% Last update 04/04/2015, Rafael NH, MRC CBU
% Output: sel =
%               1 if Diffusion and Kurtosis Tensors are computed
%               2 if only Diffusion Tensor is computed
%               3 if the invariant metrics Mentrics MD and MK are computed

hfig=figure;
set(hfig,'color',[0.8 0.8 0.8],'name','Fit method Selection',...
    'Units','normalized','position',[0.5,0.35,0.3,0.3],'Toolbar',...
    'none','NumberTitle','off','MenuBar','none');

upind=uibuttongroup('FontSize',12,...
    'BackgroundColor',[0.93 0.93 0.93],'Units','Normalized',...
    'Position',[.025 .2042 .95 0.7708],'HighlightColor','b',...
    'ForegroundColor','b');

% Botton after Method's Selected
uicontrol('Units','Normalized','Position',[0.025,0.025,0.95,0.925/6],...
    'String','Process','Callback',{@Call_Fit});

% Select Fit Method
uicontrol('Style','Radio', 'Units','normalized','Parent',upind, ...
    'HandleVisibility','off', 'Position',[0.05 19/25 0.9 3/25], ...
    'String','linear DTI (DT)', 'Tag','LDTI','FontWeight','bold')

uicontrol('Style','Radio', 'Units','normalized','Parent',upind, ...
    'HandleVisibility','off', 'Position',[0.05 15/25 0.9 3/25], ...
    'String','non-linear DTI (DT)', 'Tag','NDTI','FontWeight','bold')

uicontrol('Style','Radio', 'Units','normalized','Parent',upind, ...
    'HandleVisibility','off', 'Position',[0.05 11/25 0.9 3/25], ...
    'String','OLS DKI (DT/KT)', 'Tag','ODKI','FontWeight','bold')

uicontrol('Style','Radio', 'Units','normalized','Parent',upind, ...
    'HandleVisibility','off', 'Position',[0.05 7/25 0.9 3/25], ...
    'String','CLS DKI (DT/KT)', 'Tag','CDKI','FontWeight','bold')

uicontrol('Style','Radio', 'Units','normalized','Parent',upind, ...
    'HandleVisibility','off', 'Position',[0.05 3/25 0.9 3/25], ...
    'String','DLS DKI (MD/MK)', 'Tag','DDKI','FontWeight','bold')

% uicontrol('Style','Radio', 'Units','normalized','Parent',upind, ...
%     'HandleVisibility','off', 'Position',[0.05 1/25 0.9 3/25], ...
%     'String','RLS DKI (DT/KT)', 'Tag','RDKI','FontWeight','bold')

archstr = computer('arch');
st_comp=archstr(1:3);
if strcmp(st_comp,'win')
    sc='\';
else
    sc='/';
end

V=load_untouch_nii(DWI_name);
DWI=double(V.img);
V.hdr.dime.datatype=16;
Vm=load_untouch_nii(Mask_name);
Mask=double(Vm.img);
bval=load(bval_name);
bvec=load(bvec_name);

    function Call_Fit(source, eventdata)
        % Resume and close
        switch get(get(upind,'SelectedObject'),'Tag')
            
            case 'ODKI'
                
                fitmethod='ols';
                fun_olsDKI(DWI_name,bval_name,bvec_name,Mask_name)
                sel=1;
                
            case 'CDKI',
                
                fitmethod='cls';
                fun_clsDKI(DWI_name,bval_name,bvec_name,Mask_name)
                sel=1;
                
            case 'DDKI',
                
                fitmethod='ldf';
                fun_ldsDKI(DWI_name,bval_name,Mask_name)
                sel=3;
                
            case 'NDTI',
                
                fitmethod='nls';
                fun_nlsDTI(DWI_name,bval_name,bvec_name,Mask_name)
                sel=2;
                
            case 'LDTI'
                
                fitmethod='ols';
                fun_olsDTI(DWI_name,bval_name,bvec_name,Mask_name)
                sel=2;
                
        end
        uiresume(hfig)
        close(hfig)
    end
uiwait(hfig)
end