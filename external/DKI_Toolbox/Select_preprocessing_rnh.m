function DWI_name=Select_preprocessing_rnh(DWI_name)
% 1 July 2013

archstr = computer('arch');
st_comp=archstr(1:3);
if strcmp(st_comp,'win')
    sc='\';
else
    sc='/';
end

hfig=figure;
set(hfig,'color',[0.8 0.8 0.8],'name','Preprocessing DWI',...
'Units','normalized','position',[0.5,0.5,0.3,0.3],'Toolbar',...
'none','NumberTitle','off','MenuBar','none');

upind=uipanel('FontSize',12,...
    'BackgroundColor',[0.93 0.93 0.93],'Units','Normalized',...
    'Position',[.025 .2042 .95 0.7708],'HighlightColor','b',...
    'ForegroundColor','b');

% Botton after Method's Selected
uicontrol('Units','Normalized','Position',[0.025,0.025,0.95,0.925/6],...
   'String','Process','Callback',{@Call_Process});

% Select Erode Brain Mask
EM_sel=uicontrol('Parent',upind,'Units','Normalized','Style','Checkbox',...
    'String','Eroded Brain Mask','position',[0.05  14/27 0.6 1/9],...
    'FontWeight','bold','Callback',{@check_Mask});

% Erode Brain Mask parameters
ttDWI=uicontrol('Parent',upind,'Units','Normalized','Style','Text',...
    'String','DWI Threshold:','position',[0.35  10/27 0.3 1/9],'visible', 'off');
tDWI=uicontrol('Parent',upind,'Units','Normalized','Style','Edit',...
    'String','100','position',[0.65  10/27 0.3 1/9],'visible', 'off');
tB0vol=uicontrol('Parent',upind,'Units','Normalized','Style','Text',...
    'String','B0 vol:','position',[0.05  10/27 0.15 1/9],'visible', 'off');
B0vol=uicontrol('Parent',upind,'Units','Normalized','Style','Edit',...
    'String','1','position',[0.20  10/27 0.15 1/9],'visible', 'off');

V=load_untouch_nii(DWI_name);
DWI=double(V.img);
originalR=V.hdr.dime.pixdim(2:4);
V.hdr.dime.datatype=16;

% Downsampling
DS_sel=uicontrol('Parent',upind,'Units','Normalized','Style','Checkbox',...
    'String','Downsampling','position',[0.05  23/27 0.6 1/9],...
    'FontWeight','bold','Callback',{@check_Down});
% Final resolution
trs=uicontrol('Parent',upind,'Units','Normalized','Style','Text',...
    'String','Final Resolution (mm):','position',[0.35  19/27 0.3 1/9],...
    'visible', 'off');
rX=uicontrol('Parent',upind,'Units','Normalized','Style','Edit',...
    'String',num2str(max(originalR)),'position',[0.65  19/27 0.1 1/9],...
    'visible', 'off');
rY=uicontrol('Parent',upind,'Units','Normalized','Style','Edit',...
    'String',num2str(max(originalR)),'position',[0.75  19/27 0.1 1/9],...
    'visible', 'off');
rZ=uicontrol('Parent',upind,'Units','Normalized','Style','Edit',...
    'String',num2str(max(originalR)),'position',[0.85  19/27 0.1 1/9],...
    'visible', 'off');
tX=uicontrol('Parent',upind,'Units','Normalized','Style','Text',...
    'String','X','position',[0.65  23/27 0.1 1/9],...
    'visible', 'off');
tY=uicontrol('Parent',upind,'Units','Normalized','Style','Text',...
    'String','Y','position',[0.75  23/27 0.1 1/9],...
    'visible', 'off');
tZ=uicontrol('Parent',upind,'Units','Normalized','Style','Text',...
    'String','Z','position',[0.85  23/27 0.1 1/9],...
    'visible', 'off');

% Filter
GF_sel=uicontrol('Parent',upind,'Units','Normalized','Style','Checkbox',...
    'String','Filter','position',[0.05  5/27 0.6 1/9],...
    'FontWeight','bold','Callback',{@check_Filt});
% FWHM
tpfwhm=uicontrol('Parent',upind,'Units','Normalized','Style','Text',...
    'String','FWHM (normalized):','position',[0.35  1/27 0.3 1/9],'visible', 'off');
pfwhm=uicontrol('Parent',upind,'Units','Normalized','Style','Edit',...
    'String','1.5','position',[0.65  1/27 0.3 1/9],'visible', 'off');

    function check_Mask(source, eventdata)
        selMask=get(source,'Value');
        if selMask == 1
            %set([tDWI, B0vol, ttDWI, tB0vol],'visible','on') % This is not
            %working propertly. The variables left are useful for
            %downsampling proposes
            set([B0vol, tB0vol],'visible','on')
        else
            %set([tDWI, B0vol, ttDWI, tB0vol],'visible','off')
            %working propertly. The variables left are useful for
            %downsampling proposes
            set([B0vol, tB0vol],'visible','off')
        end
    end
    function check_Down(source, eventdata)
        selMask=get(source,'Value');
        if selMask == 1
            set([trs, rX, rY, rZ, tX, tY, tZ],'visible','on')
        else
            set([trs, rX, rY, rZ, tX, tY, tZ],'visible','off')
        end
    end
    function check_Filt(source, eventdata)
        selMask=get(source,'Value');
        if selMask == 1
            set([tpfwhm, pfwhm],'visible','on')
        else
            set([tpfwhm, pfwhm],'visible','off')
        end
    end
    function Call_Process(source, eventdata)
        % Downsampling
        sel=get(DS_sel,'value');
        if sel==1
            % read par
            rx=str2double(get(rX,'String'));
            ry=str2double(get(rY,'String'));
            rz=str2double(get(rZ,'String'));
            resfinal=[rx ry rz];
            S=size(DWI);
            Nvol=S(4);
            % process
            for vol=Nvol:-1:1
                DWId(:,:,:,vol)=resampling_for_camcan_rh(squeeze(DWI(:,:,:,vol)),originalR,resfinal);
            end
            DWI=DWId;
            V.hdr.dime.pixdim(2:4)=resfinal;
            S=size(DWId);
            V.hdr.dime.dim(2:4)=S(1:3);
            [pth,nnname,daa]=fileparts(DWI_name);
            DWI_name=[pth,sc,nnname,'_downS.nii'];
        end
        % Erode Mask
        sel=get(EM_sel,'value');
        if sel==1
            % read par
            v=str2double(get(B0vol,'String'));
            Mask=squeeze(DWI(:,:,:,v));
           % dwiT=str2double(get(B0vol,'String'));
            Mask=Mask>0;
            %SE=strel('disk', 4, 4);
            %Mask=imerode(Mask,SE);
            %Mask=imdilate(Mask,SE);
            %Mask=imclose(Mask,SE);
            %Mask=double(Mask>0);
            Masknii=V;
            Masknii.hdr.dime.dim(5)=1;
            Masknii.hdr.dime.pixdim(5)=0;
            Masknii.hdr.dime.xyz_t=0;
            Masknii.img=Mask;
            [pth,nnname,ext_ext]=fileparts(DWI_name);
            save_untouch_nii(Masknii,[pth,sc, nnname,'_EMask.nii'])
            DWI_name=[pth,sc,nnname,'_EM.nii'];
            S=size(DWI);
            for v=1:S(4)
                DWIvol=squeeze(DWI(:,:,:,v));
                DWIvol(Mask==0)=0;
                DWI(:,:,:,v)=DWIvol;
            end
        end        
        % Filter
        sel=get(GF_sel,'value');
        if sel==1
            % read par
            ssd=get(pfwhm,'String');
            fwhm=str2double(ssd);
            sd=fwhm/(sqrt(8*log(2))); % Convert FWHM to sd
            S=size(DWI);
            Nvol=S(4);
            % Filtro gaussiano
            size_k=[7 7 7];
            for im=Nvol:-1:1
                A=squeeze(DWI(:,:,:,im));
                W = smooth3(A,'gaussian',size_k,sd);
                DWI(:,:,:,im)=W;
            end
            [pth,nnname,ext_ext]=fileparts(DWI_name);
            DWI_name=[pth,sc,nnname,'_GF',num2str(fwhm*100),'.nii'];
        end
        V.img=DWI;
        save_untouch_nii(V,DWI_name)
        % Resume and close
        uiresume(hfig)
        close(hfig)
    end
uiwait(hfig)
end