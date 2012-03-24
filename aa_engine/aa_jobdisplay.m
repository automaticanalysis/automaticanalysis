function varargout = aa_jobdisplay(varargin)
% AA_JOBDISPLAY M-file for aa_jobdisplay.fig
%      AA_JOBDISPLAY, by itself, creates a new AA_JOBDISPLAY or raises the existing
%      singleton*.
%
%      H = AA_JOBDISPLAY returns the handle to a new AA_JOBDISPLAY or the handle to
%      the existing singleton*.
%
%      AA_JOBDISPLAY('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in AA_JOBDISPLAY.M with the given input arguments.
%
%      AA_JOBDISPLAY('Property','Value',...) creates a new AA_JOBDISPLAY or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before aa_jobdisplay_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to aa_jobdisplay_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help aa_jobdisplay

% Last Modified by GUIDE v2.5 26-Aug-2008 10:03:33

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @aa_jobdisplay_OpeningFcn, ...
                   'gui_OutputFcn',  @aa_jobdisplay_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before aa_jobdisplay is made visible.
function aa_jobdisplay_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to aa_jobdisplay (see VARARGIN)

% Choose default command line output for aa_jobdisplay
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes aa_jobdisplay wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = aa_jobdisplay_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
