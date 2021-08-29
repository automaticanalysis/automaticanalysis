function varargout = BrainNet_LoadFiles(varargin)
%BrainNet Viewer, a graph-based brain network mapping tool, by Mingrui Xia
%Function to load files for graph drawing
%-----------------------------------------------------------
%	Copyright(c) 2017
%	State Key Laboratory of Cognitive Neuroscience and Learning, Beijing Normal University
%	Written by Mingrui Xia
%	Mail to Author:  <a href="mingruixia@gmail.com">Mingrui Xia</a>
%   Version 1.6;
%   Date 20110531;
%   Last edited 20170330
%-----------------------------------------------------------
%


% BrainNet_LoadFiles MATLAB code for BrainNet_LoadFiles.fig
%      BrainNet_LoadFiles, by itself, creates a new BrainNet_LoadFiles or raises the existing
%      singleton*.
%
%      H = BrainNet_LoadFiles returns the handle to a new BrainNet_LoadFiles or the handle to
%      the existing singleton*.
%
%      BrainNet_LoadFiles('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BrainNet_LoadFiles.M with the given input arguments.
%
%      BrainNet_LoadFiles('Property','Value',...) creates a new BrainNet_LoadFiles or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before BrainNet_LoadFiles_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to BrainNet_LoadFiles_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help BrainNet_LoadFiles

% Last Modified by GUIDE v2.5 26-Oct-2011 19:54:43

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @BrainNet_LoadFiles_OpeningFcn, ...
    'gui_OutputFcn',  @BrainNet_LoadFiles_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin
    if ischar(varargin{1}), gui_State.gui_Callback = str2func(varargin{1}); end
    if isstruct(varargin{1})
        OK_button_Callback([], [], [], varargin{1});
        return
    end
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before BrainNet_LoadFiles is made visible.
function BrainNet_LoadFiles_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to BrainNet_LoadFiles (see VARARGIN)

% Choose default command line output for BrainNet_LoadFiles
handles.output = hObject;
guidata(hObject, handles);
h_NV=findobj('Tag','NV_fig');
h_NV=guihandles(h_NV);
setappdata(handles.LF_fig,'h_NV',h_NV);
movegui(handles.LF_fig,'center');
global File
set(handles.MF_edit,'string',File.MF);
set(handles.NI_edit,'string',File.NI);
set(handles.NT_edit,'string',File.NT);
set(handles.VF_edit,'string',File.VF);
global FLAG
FLAG.LF=0;


% Update handles structure


% UIWAIT makes BrainNet_LoadFiles wait for user response (see UIRESUME)
% uiwait(handles.LF_fig);


% --- Outputs from this function are returned to the command line.
function varargout = BrainNet_LoadFiles_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in OK_button.
function OK_button_Callback(hObject, eventdata, handles, varargin)
% hObject    handle to OK_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global File
global FLAG
global EC
File.MF=[];
File.NI=[];
File.NT=[];
File.VF=[];
if nargin >= 4
    File.MF = varargin{1}.surf;
    File.NI = varargin{1}.node;
    File.NT = varargin{1}.edge;
    File.VF = varargin{1}.map;
else
    File.MF=get(handles.MF_edit,'string');
    File.NI=get(handles.NI_edit,'string');
    File.NT=get(handles.NT_edit,'string');
    File.VF=get(handles.VF_edit,'String');
end
if isempty(File.MF)
    mf=0;
else
    mf=1;
end
if isempty(File.NI)
    ni=0;
else
    ni=2;
end
if isempty(File.NT)
    nt=0;
else
    nt=4;
end
if isempty(File.VF)
    vf=0;
else
    vf=8;
end

FLAG.Loadfile=mf+ni+nt+vf;
global surf
if FLAG.Loadfile==0
    h=msgbox('Please select file!','Error','error');
    uiwait(h);
    return;
else
    if FLAG.Loadfile==4||FLAG.Loadfile==5|| FLAG.Loadfile == 13 %%% Edited by Mingrui Xia 20120219, add volume, node and edge mode
        h=msgbox('Please select node file!','Error','error');
        uiwait(h);
        return;
        %     elseif FLAG.Loadfile>9 %%%Edited by Mingrui Xia 20111116, add volume
        %     and & mode
    elseif FLAG.Loadfile == 8 || FLAG.Loadfile == 10 || FLAG.Loadfile == 12 || FLAG.Loadfile == 14
        h=msgbox('Please select surf file!','Error','error');
        uiwait(h);
        return;
    else
        surf=[];
        switch FLAG.Loadfile
            case 1
                [surf.vertex_number surf.coord surf.ntri surf.tri]=MF_load(File.MF);
                EC.msh.alpha = 1;
            case 2
                [surf.nsph surf.sphere surf.label]=NI_load(File.NI);
            case 3
                [surf.vertex_number surf.coord surf.ntri surf.tri]=MF_load(File.MF);
                [surf.nsph surf.sphere surf.label]=NI_load(File.NI);
            case 6
                [surf.nsph surf.sphere surf.label]=NI_load(File.NI);
                surf.net=NT_load(File.NT);
            case 7
                [surf.vertex_number surf.coord surf.ntri surf.tri]=MF_load(File.MF);
                [surf.nsph surf.sphere surf.label]=NI_load(File.NI);
                surf.net=NT_load(File.NT);
            case 9
                [surf.vertex_number surf.coord surf.ntri surf.tri]=MF_load(File.MF);
                [path, fname,ext]=fileparts(File.VF);
                switch ext
                    case '.txt'
                        [surf.T]=load(File.VF);
                        surf.test = 'No';
                        FLAG.MAP=1;
                        
                        % test code
                        %                         facevalue = zeros(size(surf.tri));
                        %                         facevalue(:) = surf.T(surf.tri(:));
                        %                         facevalue_judge = 1 - all(diff(facevalue,[],2)==0,2);
                        %                         face_remove = surf.tri(find(facevalue_judge),:);
                        %                         surf.T(unique(face_remove)) = 0;
                        
                        
                        
                        
                        
                    case '.gz' % Added by Mingrui, 20120611, support for .nii.gz files
                        tmp_folder = tempdir;
                        gunzip(File.VF,tmp_folder);
                        filename = [tmp_folder,fname];
                        [surf.hdr,surf.vol]=VF_load(filename);
                        FLAG.MAP=2;
                        FLAG.ROICheck = 0;
                        delete(filename);
                    case '.annot'
                        [surf.T] = VF_load_a(File.VF);
                        surf.test = 'No';
                        FLAG.MAP = 1;
                    otherwise
                        [surf.hdr,surf.vol]=VF_load(File.VF);
                        FLAG.MAP=2;
                        FLAG.ROICheck = 0;
                end
                
                % Added by Mingrui, 20140924,  statistic for SPM or REST nifti files
                if FLAG.MAP == 2
                    if isempty(strfind(surf.hdr.descrip,'REST'))&&isempty(strfind(surf.hdr.descrip,'SPM'))&&isempty(strfind(surf.hdr.descrip,'DPABI'))
                        surf.test = 'No';
                    else
                        % Modified from REST rest_sliceviwer.m
                        if ~isempty(strfind(surf.hdr.descrip,'{T_['))
                            surf.test = 'T';
                            Tstart = strfind(surf.hdr.descrip,'{T_[')+length('{T_[');
                            Tend = strfind(surf.hdr.descrip,']}')-1;
                            surf.df = str2num(surf.hdr.descrip(Tstart:Tend));
                        elseif ~isempty(strfind(surf.hdr.descrip,'{F_['))
                            surf.test = 'F';
                            Tstart = strfind(surf.hdr.descrip,'{F_[')+length('{F_[');
                            Tend = strfind(surf.hdr.descrip,']}')-1;
                            surf.df = str2num(surf.hdr.descrip(Tstart:Tend));
                        elseif ~isempty(strfind(surf.hdr.descrip,'{R_['))
                            surf.test = 'R';
                            Tstart = strfind(surf.hdr.descrip,'{R_[')+length('{R_[');
                            Tend = strfind(surf.hdr.descrip,']}')-1;
                            surf.df = str2num(surf.hdr.descrip(Tstart:Tend));
                        elseif ~isempty(strfind(surf.hdr.descrip,'{Z_['))
                            surf.test = 'Z';
                        else
                            surf.test = 'No';
                        end
                    end
                end
                
                EC.msh.alpha = 1;
            case 11
                [surf.vertex_number surf.coord surf.ntri surf.tri]=MF_load(File.MF);
                [path, fname,ext]=fileparts(File.VF);
                switch ext
                    case '.txt'
                        [surf.T]=load(File.VF);
                        surf.test = 'No';
                        FLAG.MAP=1;
                    case '.gz' % Added by Mingrui, 20120611, support for .nii.gz files
                        tmp_folder = tempdir;
                        gunzip(File.VF,tmp_folder);
                        filename = [tmp_folder,fname];
                        [surf.hdr surf.vol]=VF_load(filename);
                        FLAG.MAP=2;
                        delete(filename);
                        FLAG.ROICheck = 0;
                    case '.annot'
                        [surf.T] = VF_load_a(File.VF);
                        surf.test = 'No';
                        FLAG.MAP = 1;
                    otherwise
                        [surf.hdr surf.vol]=VF_load(File.VF);
                        FLAG.MAP=2;
                        FLAG.ROICheck = 0;
                end
                
                % Added by Mingrui, 20140924,  statistic for SPM or REST nifti files
                if FLAG.MAP == 2
                    if isempty(strfind(surf.hdr.descrip,'REST'))&&isempty(strfind(surf.hdr.descrip,'SPM'))&&isempty(strfind(surf.hdr.descrip,'DPABI'))
                        surf.test = 'No';
                    else
                        % Modified from REST rest_sliceviwer.m
                        if ~isempty(strfind(surf.hdr.descrip,'{T_['))
                            surf.test = 'T';
                            Tstart = strfind(surf.hdr.descrip,'{T_[')+length('{T_[');
                            Tend = strfind(surf.hdr.descrip,']}')-1;
                            surf.df = str2num(surf.hdr.descrip(Tstart:Tend));
                        elseif ~isempty(strfind(surf.hdr.descrip,'{F_['))
                            surf.test = 'F';
                            Tstart = strfind(surf.hdr.descrip,'{F_[')+length('{F_[');
                            Tend = strfind(surf.hdr.descrip,']}')-1;
                            surf.df = str2num(surf.hdr.descrip(Tstart:Tend));
                        elseif ~isempty(strfind(surf.hdr.descrip,'{R_['))
                            surf.test = 'R';
                            Tstart = strfind(surf.hdr.descrip,'{R_[')+length('{R_[');
                            Tend = strfind(surf.hdr.descrip,']}')-1;
                            surf.df = str2num(surf.hdr.descrip(Tstart:Tend));
                        elseif ~isempty(strfind(surf.hdr.descrip,'{Z_['))
                            surf.test = 'Z';
                        else % Added by Mingrui 20190723, fix a bug when header description includes toolbox name but not a statistical map
                            surf.test = 'No';
                            
                        end
                    end
                end
                
                [surf.nsph surf.sphere surf.label]=NI_load(File.NI);
            case 15 %%% Added by Mingrui Xia, 20100210, add volume, node and edge mode.
                [surf.vertex_number surf.coord surf.ntri surf.tri]=MF_load(File.MF);
                [path, fname,ext]=fileparts(File.VF);
                switch ext
                    case '.txt'
                        [surf.T]=load(File.VF);
                        surf.test = 'No';
                        FLAG.MAP=1;
                    case '.gz' % Added by Mingrui, 20120611, support for .nii.gz files
                        tmp_folder = tempdir;
                        gunzip(File.VF,tmp_folder);
                        filename = [tmp_folder,fname];
                        [surf.hdr surf.vol]=VF_load(filename);
                        FLAG.MAP=2;
                        delete(filename);
                        FLAG.ROICheck = 0;
                    case '.annot'
                        [surf.T] = VF_load_a(File.VF);
                        surf.test = 'No';
                        FLAG.MAP = 1;
                    otherwise
                        [surf.hdr surf.vol]=VF_load(File.VF);
                        FLAG.MAP=2;
                        FLAG.ROICheck = 0;
                end
                
                % Added by Mingrui, 20140924,  statistic for SPM or REST nifti files
                if FLAG.MAP == 2
                    if isempty(strfind(surf.hdr.descrip,'REST'))&&isempty(strfind(surf.hdr.descrip,'SPM'))&&isempty(strfind(surf.hdr.descrip,'DPABI'))
                        surf.test = 'No';
                    else
                        % Modified from REST rest_sliceviwer.m
                        if ~isempty(strfind(surf.hdr.descrip,'{T_['))
                            surf.test = 'T';
                            Tstart = strfind(surf.hdr.descrip,'{T_[')+length('{T_[');
                            Tend = strfind(surf.hdr.descrip,']}')-1;
                            surf.df = str2num(surf.hdr.descrip(Tstart:Tend));
                        elseif ~isempty(strfind(surf.hdr.descrip,'{F_['))
                            surf.test = 'F';
                            Tstart = strfind(surf.hdr.descrip,'{F_[')+length('{F_[');
                            Tend = strfind(surf.hdr.descrip,']}')-1;
                            surf.df = str2num(surf.hdr.descrip(Tstart:Tend));
                        elseif ~isempty(strfind(surf.hdr.descrip,'{R_['))
                            surf.test = 'R';
                            Tstart = strfind(surf.hdr.descrip,'{R_[')+length('{R_[');
                            Tend = strfind(surf.hdr.descrip,']}')-1;
                            surf.df = str2num(surf.hdr.descrip(Tstart:Tend));
                        elseif ~isempty(strfind(surf.hdr.descrip,'{Z_['))
                            surf.test = 'Z';
                        else % Added by Mingrui 20190723, fix a bug when header description includes toolbox name but not a statistical map
                            surf.test = 'No';
                        end
                    end
                end
                
                [surf.nsph surf.sphere surf.label]=NI_load(File.NI);
                surf.net=NT_load(File.NT);
        end
        FLAG.LF=1;
        close(findobj('Tag','LF_fig'));
    end
end

function t = VF_load_a(filename)

global EC

fp = fopen(filename, 'r', 'b');

verbosity = 1;

if(fp < 0)
    return;
end

A = fread(fp, 1, 'int');
tmp = fread(fp, 2*A, 'int');
label = tmp(2:2:end);

bool = fread(fp, 1, 'int');
if(isempty(bool)) %means no colortable
    colortable = struct([]);
    fclose(fp);
    return;
end

if(bool)
    %Read colortable
    numEntries = fread(fp, 1, 'int');
    if(numEntries > 0)
        colortable.numEntries = numEntries;
        len = fread(fp, 1, 'int');
        colortable.orig_tab = fread(fp, len, '*char')';
        colortable.orig_tab = colortable.orig_tab(1:end-1);
        colortable.struct_names = cell(numEntries,1);
        colortable.table = zeros(numEntries,5);
        for i = 1:numEntries
            len = fread(fp, 1, 'int');
            colortable.struct_names{i} = fread(fp, len, '*char')';
            colortable.struct_names{i} = colortable.struct_names{i}(1:end-1);
            colortable.table(i,1) = fread(fp, 1, 'int');
            colortable.table(i,2) = fread(fp, 1, 'int');
            colortable.table(i,3) = fread(fp, 1, 'int');
            colortable.table(i,4) = fread(fp, 1, 'int');
            colortable.table(i,5) = colortable.table(i,1) + colortable.table(i,2)*2^8 + colortable.table(i,3)*2^16 + colortable.table(i,4)*2^24;
        end
        if verbosity
            disp(['colortable with ' num2str(colortable.numEntries) ' entries read (originally ' colortable.orig_tab ')']);
        end
    else
        version = -numEntries;
        if verbosity
            if(version~=2)
                disp(['Error! Does not handle version ' num2str(version)]);
            else
                disp(['Reading from version ' num2str(version)]);
            end
        end
        numEntries = fread(fp, 1, 'int');
        colortable.numEntries = numEntries;
        len = fread(fp, 1, 'int');
        colortable.orig_tab = fread(fp, len, '*char')';
        colortable.orig_tab = colortable.orig_tab(1:end-1);
        colortable.struct_names = cell(numEntries,1);
        colortable.table = zeros(numEntries,5);
        numEntriesToRead = fread(fp, 1, 'int');
        for i = 1:numEntriesToRead
            structure = fread(fp, 1, 'int')+1;
            if (structure < 0)
                if verbosity, disp(['Error! Read entry, index ' num2str(structure)]); end;
            end
            if(~isempty(colortable.struct_names{structure}))
                if verbosity, disp(['Error! Duplicate Structure ' num2str(structure)]); end;
            end
            len = fread(fp, 1, 'int');
            colortable.struct_names{structure} = fread(fp, len, '*char')';
            colortable.struct_names{structure} = colortable.struct_names{structure}(1:end-1);
            colortable.table(structure,1) = fread(fp, 1, 'int');
            colortable.table(structure,2) = fread(fp, 1, 'int');
            colortable.table(structure,3) = fread(fp, 1, 'int');
            colortable.table(structure,4) = fread(fp, 1, 'int');
            colortable.table(structure,5) = colortable.table(structure,1) + colortable.table(structure,2)*2^8 + colortable.table(structure,3)*2^16 + colortable.table(structure,4)*2^24;
        end
    end
end
fclose(fp);

tmp = sortrows(colortable.table,5);
tmp_t = label;
for i = 1:size(tmp,1)
    tmp_t(tmp_t==tmp(i,5)) = i - 1;
end
t = tmp_t;
EC.vol.CM_annot = tmp(:,1:3)/255;
EC.vol.color_map = 25;




function [hdr,vol]=VF_load(filename)
% YAN Chao-Gan 111028. Add the path of BrainNet SPM files every time.
% [BrainNetPath, fileN, extn] = fileparts(which('BrainNet.m'));
[BrainNetPath] = fileparts(which('BrainNet.m')); %%% Edited by Mingrui Xia, 20111112, remove two unused var.
BrainNet_SPMPath = fullfile(BrainNetPath, 'BrainNet_spm12_files'); % Edited by Mingrui Xia 20181020, use updated spm12 files
% rmpath(BrainNet_SPMPath); %%% Edited by Mingrui Xia, 20111116, clear
% warning
if ~exist('BrainNet_spm_file.m','file') %%% Edited by Mingrui Xia, 111103, check if SPM is installed.
    %    hdr=spm_vol(filename); %%% Edited by Mingrui Xia, 111026, integrated SPM NIFTI into BrainNet Viewer.
    %    vol=spm_read_vols(hdr);
    %else
    addpath(BrainNet_SPMPath);
end
hdr=BrainNet_spm_vol(filename); %%% Edited by Mingrui Xia, 111026, integrated SPM NIFTI into BrainNet Viewer.
vol=BrainNet_spm_read_vols(hdr);
%    rmpath(BrainNet_SPMPath);
%end
vol(isnan(vol)) = 0; % Added by Mingrui, 20170605, replace NaN to 0;
% Nii = nifti(filename);
% hdr.descrip = Nii.descrip;
% hdr.dim = Nii.dat.dim;
% hdr.mat = Nii.mat;
% vol = double(Nii.dat);






% --- Executes on button press in MF_button.
function MF_button_Callback(hObject, eventdata, handles)
% hObject    handle to MF_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[BrainSurfPath] = fileparts(which('BrainMesh_ICBM152_smoothed.nv')); %%% Edited by Mingrui Xia, 20120412, add default template path to ICBM152.

% Edited by Mingrui Xia, 20120918, add support for BYU '*.g' file.
if ispc % Edited by Mingrui, 20170321, fix the bug that couldn't locate the surface template in Linux
    [filename,pathname]=uigetfile({'*.nv','NetViewer Files (*.nv)';'*.mesh',...
        'BrainVISA Mesh (*.mesh)';'*.pial','FreeSurfer Mesh (*.pial)';...
        '*.g','BYU file (*.g)';'*.obj','Objective Files (*.obj)';...
        '*.gii','GIfTI Files (*.gii)';...
        '*.mz3','Surf Ice Files (*.mz3)';...
        '*.*','All Files (*.*)'},...
        'Select brain template',[BrainSurfPath,'\BrainMesh_ICBM152_smoothed.nv']);
else
    [filename,pathname]=uigetfile({'*.nv','NetViewer Files (*.nv)';'*.mesh',...
        'BrainVISA Mesh (*.mesh)';'*.pial','FreeSurfer Mesh (*.pial)';...
        '*.g','BYU file (*.g)';'*.obj','Objective Files (*.obj)';...
        '*.gii','GIfTI Files (*.gii)';'*.mz3','Surf Ice Files (*.mz3)';...
        '*.*','All Files (*.*)'},...
        'Select brain template',[BrainSurfPath,'/BrainMesh_ICBM152_smoothed.nv']);
end

if isequal(filename,0)||isequal(pathname,0)
    return;
else
    fpath=fullfile(pathname,filename);
    set(handles.MF_edit,'string',fpath);
end


function MF_edit_Callback(hObject, eventdata, handles)
% hObject    handle to MF_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MF_edit as text
%        str2double(get(hObject,'String')) returns contents of MF_edit as a double


% --- Executes during object creation, after setting all properties.
function MF_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MF_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function [vertex, faces, vertex_number, face_number] =loadpial(filename)
fid = fopen(filename, 'rb', 'b') ;
b1 = fread(fid, 1, 'uchar') ;
b2 = fread(fid, 1, 'uchar') ;
b3 = fread(fid, 1, 'uchar') ;
% magic = bitshift(b1, 16) + bitshift(b2,8) + b3 ;
fgets(fid);
fgets(fid);
v = fread(fid, 1, 'int32') ;
t = fread(fid, 1, 'int32') ;
vertex= fread(fid, [3 v], 'float32') ;
faces= fread(fid, [3 t], 'int32')' + 1 ;
fclose(fid) ;
vertex_number=size(vertex,2);
face_number=size(faces,1);


function [vertex_number coord ntri tri]=MF_load(MF)
[path, fname, ext] = fileparts(MF);
switch ext
    case '.pial'
        [coord, tri, vertex_number, ntri] =loadpial(MF);
    case '.nv'
        fid=fopen(MF);
        %         vertex_number=fscanf(fid,'%f',1);
        %         coord=fscanf(fid,'%f',[3,vertex_number]);
        %         ntri=fscanf(fid,'%f',1);
        %         tri=fscanf(fid,'%d',[3,ntri])';
        % modified by Mingrui 20141030, add support for comments
        data = textscan(fid,'%f','CommentStyle','#');
        vertex_number = data{1}(1);
        coord  = reshape(data{1}(2:1+3*vertex_number),[3,vertex_number]);
        ntri = data{1}(3*vertex_number+2);
        tri = reshape(data{1}(3*vertex_number+3:end),[3,ntri])';
        fclose(fid);
    case '.mesh' %modified from BrainVISA loadmesh.m, Copyright (C) 2003 Denis Schwartz & Guillaume Flandin
        fid=fopen(MF,'r');
        [file_format, COUNT]=fread(fid, 5, 'uchar');  %- 'ascii' or 'binar'
        switch char(file_format)'
            case 'binar'
                [byte_swapping, COUNT]=fread(fid, 1, 'uint32'); %- 'ABCD' or 'DCBA'
                ff = strcmp(dec2hex(byte_swapping),'41424344');
                if ~ff
                    [fn, pm, mf]=fopen(1); %- machine format
                    fclose(fid);
                    if strmatch(mf,'ieee-le');
                        fid=fopen(filename,'r','ieee-be');
                    else
                        fid=fopen(filename,'r','ieee-le');
                    end
                    [file_format, COUNT]=fread(fid, 5, 'uchar');
                    [byte_swapping, COUNT]=fread(fid, 1, 'uint32');
                end
                [arg_size, COUNT]= fread(fid, 1, 'uint32'); %- length('VOID')
                [VOID, COUNT]= fread(fid, arg_size, 'uchar'); %- VOID
                [polygon_dimension, COUNT]=fread(fid, 1, 'uint32'); %- 3 for triangles
                [mesh_time, COUNT]=fread(fid, 1, 'uint32'); %- number of meshes
                vertex=cell(1,mesh_time);
                normals=cell(1,mesh_time);
                faces=cell(1,mesh_time);
                for i=1:mesh_time
                    [mesh_step, COUNT]=fread(fid, 1, 'uint32'); %- [0 ... mesh_time-1]
                    %- Get vertices
                    [vertex_number, COUNT]=fread(fid, 1, 'uint32');
                    [vtx, COUNT]=fread(fid, 3*vertex_number, 'float32');
                    vertex{i}=reshape(vtx, 3, vertex_number)';
                    %- Get normals
                    [normal_number, COUNT]=fread(fid, 1, 'uint32');
                    [nrml, COUNT]=fread(fid, 3*normal_number, 'float32');
                    normal{i}=reshape(nrml, 3, normal_number)';
                    [arg_size, COUNT]=fread(fid, 1, 'uint32'); %- no data ('VOID')
                    %- Get faces
                    [faces_number, COUNT]  = fread(fid, 1, 'uint32');
                    [fcs, COUNT] = fread(fid, polygon_dimension*faces_number, 'uint32');
                    faces{i} = reshape(fcs, polygon_dimension, faces_number)';
                end
            case 'ascii'
                VOID = fscanf(fid,'%s',1);
                polygon_dimension = fscanf(fid,'%d',1);
                mesh_time = fscanf(fid,'%d',1);
                for i=1:mesh_time
                    mesh_step = fscanf(fid,'\n%d',1);
                    vertex_number = fscanf(fid,'\n%d\n',1);
                    vtx = fscanf(fid,'(%f ,%f ,%f) ',3*vertex_number);
                    vertex{i} = reshape(vtx, 3, vertex_number)';
                    normal_number = fscanf(fid,'\n%d\n',1);
                    nrml = fscanf(fid,'(%f ,%f ,%f) ',3*normal_number);
                    normal{i} = reshape(nrml, 3, normal_number)';
                    arg_size = fscanf(fid,'\n%d\n',1);
                    faces_number = fscanf(fid,'\n%d\n',1);
                    fcs = fscanf(fid,'(%d ,%d ,%d) ',polygon_dimension*faces_number);
                    faces{i} = reshape(fcs, polygon_dimension, faces_number)';
                end
        end
        if mesh_time == 1
            vertex = vertex{1};
            normal = normal{1};
            faces = faces{1};
        end
        fclose(fid);
        coord=vertex';
        tri=faces+1;
        ntri=faces_number;
        %         coord(1,:)=91-coord(1,:);
        %         coord(2,:)=91-coord(2,:);
        %         coord(3,:)=109-coord(3,:);
        coord(1,:)=-coord(1,:);
        coord(2,:)=-coord(2,:);
        coord(3,:)=-coord(3,:);
        
        % Edited by Mingrui Xia, 20120918, add support for BYU '*.g' file.
    case '.g'
        fid=fopen(MF);
        fscanf(fid,'%f',1);
        vertex_number = fscanf(fid,'%f',1);
        ntri = fscanf(fid,'%f',1);
        fscanf(fid,'%f',3);
        coord = fscanf(fid,'%f',[3,vertex_number]);
        tri = fscanf(fid,'%d',[3,ntri])';
        tri(:,3) = -tri(:,3);
        fclose(fid);
        
        %    Edited by Mingrui Xia, 20131215, add support for '*.obj'file.
    case '.obj' % Modified from SurfStat SurfStatReadSurf1.m
        fid=fopen(MF);
        FirstChar=fscanf(fid,'%1s',1);
        if FirstChar=='P' % ASCII
            fscanf(fid,'%f',5);
            vertex_number=fscanf(fid,'%f',1);
            coord=fscanf(fid,'%f',[3,vertex_number]);
            fscanf(fid,'%f',[3,vertex_number]);
            ntri=fscanf(fid,'%f',1);
            ind=fscanf(fid,'%f',1);
            if ind==0
                fscanf(fid,'%f',4);
            else
                fscanf(fid,'%f',[4,vertex_number]);
            end
            fscanf(fid,'%f',ntri);
            tri=fscanf(fid,'%f',[3,ntri])'+1;
            fclose(fid);
        else
            fclose(fid);
            fid=fopen(filename,'r','b');
            FirstChar=fread(fid,1);
            if FirstChar==uint8(112) % binary
                fread(fid,5,'float');
                vertex_number=fread(fid,1,'int');
                coord=fread(fid,[3,vertex_number],'float');
                
                fread(fid,[3,vertex_number],'float');
                
                ntri=fread(fid,1,'int');
                ind=fread(fid,1,'int');
                if ind==0
                    uint8(fread(fid,4,'uint8'));
                else
                    uint8(fread(fid,[4,vertex_number],'uint8'));
                end
                fread(fid,ntri,'int');
                tri=fread(fid,[3,ntri],'int')'+1;
                fclose(fid);
            end
        end
    case '.gii' % Added by Mingrui, 20170309, load gii by using gifti toolbox
        g = gifti(MF);
        vertex_number = size(g.vertices,1);
        coord = g.vertices';
        ntri = size(g.faces,1);
        tri = g.faces;
    case '.mz3'% Added by Mingrui, 20170330, support mz3 from Surf_Ice, modified from readMz3.m, https://github.com/bonilhamusclab/MRIcroS/blob/development/%2BfileUtils/%2Bmz3/readMz3.m
        %Decode gzip data
        % http://undocumentedmatlab.com/blog/savezip-utility
        % http://www.mathworks.com/matlabcentral/fileexchange/39526-byte-encoding-utilities/content/encoder/gzipdecode.m
        streamCopier = com.mathworks.mlwidgets.io.InterruptibleStreamCopier.getInterruptibleStreamCopier;
        baos = java.io.ByteArrayOutputStream;
        fis  = java.io.FileInputStream(MF);
        zis  = java.util.zip.GZIPInputStream(fis);
        streamCopier.copyStream(zis,baos);
        fis.close;
        data = baos.toByteArray;
        %mz3 ALWAYS little endian
        machine = 'ieee-le';
        magic = typecast(data(1:2),'uint16');
        if magic ~= 23117, fprintf('Signature is not MZ3\n'); return; end;
        %attr reports attributes and version
        attr = typecast(data(3:4),'uint16');
        if (attr == 0) || (attr > 7), fprintf('This file uses unsupported features\n'); end;
        isFace = bitand(attr,1);
        isVert = bitand(attr,2);
        isRGBA = bitand(attr,4);
        isSCALAR = bitand(attr,8);
        %read attributes
        nFace = typecast(data(5:8),'uint32');
        nVert = typecast(data(9:12),'uint32');
        nSkip = typecast(data(13:16),'uint32');
        hdrSz = 16+nSkip; %header size in bytes
        %read faces
        if isFace
            facebytes = nFace * 3 * 4; %each face has 3 indices, each 4 byte int
            faces = typecast(data(hdrSz+1:hdrSz+facebytes),'int32');
            faces = double(faces')+1; %matlab indices arrays from 1 not 0
            %faces = reshape(faces,3,nFace)';
            faces = reshape(faces,3, nFace)';
            hdrSz = hdrSz + facebytes;
        end;
        %read vertices
        if isVert
            vertbytes = nVert * 3 * 4; %each vertex has 3 values (x,y,z), each 4 byte float
            vertices = typecast(data(hdrSz+1:hdrSz+vertbytes),'single');
            vertices = double(vertices); %matlab wants doubles
            %vertices = reshape(vertices,nVert,3);
            vertices = reshape(vertices,3,nVert);
            hdrSz = hdrSz + vertbytes;
        end
        vertex_number = nVert;
        coord = vertices;
        ntri = nFace;
        tri = faces;
end


function [nsph, sphere, label]=NI_load(NI)
% fid=fopen(NI);
% i=0;
% while ~feof(fid)
%     curr=fscanf(fid,'%f',5);
%     if ~isempty(curr)
%         i=i+1;
%         textscan(fid,'%s',1);
%     end
% end
% nsph=i;
% fclose(fid);
% sphere=zeros(nsph,5);
% label=cell(nsph,1);
% fid=fopen(NI);
% i=0;
% while ~feof(fid)
%     curr=fscanf(fid,'%f',5);
%     if ~isempty(curr)
%         i=i+1;
%         sphere(i,1:5)=curr;
%         label{i}=textscan(fid,'%s',1);
%     end
% end
% fclose(fid);

% modified by Mingrui 20141028, add support for comments in node file according to Chris Rorden's suggestion
fid = fopen(NI);
data = textscan(fid,'%f %f %f %f %f %s','CommentStyle','#');
fclose(fid);
sphere = [cell2mat(data(1)) cell2mat(data(2)) cell2mat(data(3)) cell2mat(data(4)) cell2mat(data(5))];
label = data{6};
nsph = size(sphere,1);


function net=NT_load(NT)
% net=load(NT);
% modified by Mingrui 20141030, add support for comments in edge file
fid = fopen(NT);
data = textscan(fid,'%f','CommentStyle','#');
net = reshape(data{1},[sqrt(length(data{1})),sqrt(length(data{1}))])';
fclose(fid);


% --- Executes on button press in Cancel_button.
function Cancel_button_Callback(hObject, eventdata, handles)
% hObject    handle to Cancel_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global FLAG
FLAG.LF=0;
close(findobj('Tag','LF_fig'));


function NI_edit_Callback(hObject, eventdata, handles)
% hObject    handle to NI_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NI_edit as text
%        str2double(get(hObject,'String')) returns contents of NI_edit as a double


% --- Executes during object creation, after setting all properties.
function NI_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NI_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function NT_edit_Callback(hObject, eventdata, handles)
% hObject    handle to NT_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NT_edit as text
%        str2double(get(hObject,'String')) returns contents of NT_edit as a double


% --- Executes during object creation, after setting all properties.
function NT_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NT_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in NI_button.
function NI_button_Callback(hObject, eventdata, handles)
% hObject    handle to NI_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename,pathname]=uigetfile({'*.node','Node files (*.node)';'*.*','All Files (*.*)'});
if isequal(filename,0)||isequal(pathname,0)
    return;
else
    fpath=fullfile(pathname,filename);
    set(handles.NI_edit,'string',fpath);
end

% --- Executes on button press in NT_button.
function NT_button_Callback(hObject, eventdata, handles)
% hObject    handle to NT_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename,pathname]=uigetfile({'*.edge','Edge files (*.edge)';'*.*','All Files (*.*)'});
if isequal(filename,0)||isequal(pathname,0)
    return;
else
    fpath=fullfile(pathname,filename);
    set(handles.NT_edit,'string',fpath);
end


% --- Executes on button press in Reset_button.
function Reset_button_Callback(hObject, eventdata, handles)
% hObject    handle to Reset_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.MF_edit,'string','');
set(handles.NI_edit,'string','');
set(handles.NT_edit,'string','');
set(handles.VF_edit,'string','');



function NL_edit_Callback(hObject, eventdata, handles)
% hObject    handle to NL_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NL_edit as text
%        str2double(get(hObject,'String')) returns contents of NL_edit as a double


% --- Executes during object creation, after setting all properties.
function NL_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NL_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in NL_button.
function NL_button_Callback(hObject, eventdata, handles)
% hObject    handle to NL_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename,pathname]=uigetfile({'*.txt','Text files (*.txt)';'*.*','All Files (*.*)'});
if isequal(filename,0)||isequal(pathname,0)
    return;
else
    fpath=fullfile(pathname,filename);
    set(handles.NL_edit,'string',fpath);
end



function VF_edit_Callback(hObject, eventdata, handles)
% hObject    handle to VF_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of VF_edit as text
%        str2double(get(hObject,'String')) returns contents of VF_edit as a double


% --- Executes during object creation, after setting all properties.
function VF_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to VF_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in VF_pushbutton.
function VF_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to VF_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename,pathname]=uigetfile({'*.nii;*.hdr;*.img;*.nii.gz','NIFTI files (*.nii,*.hdr,*.img,*.nii.gz)';'*.txt','Text files (*.txt)';'*.annot','Annot files (*.annot)';'*.*','All Files (*.*)'});
if isequal(filename,0)||isequal(pathname,0)
    return;
else
    fpath=fullfile(pathname,filename);
    set(handles.VF_edit,'string',fpath);
end
