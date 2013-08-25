function varargout = aa_parallel_status(varargin)
% AA_PARALLEL_STATUS M-file for aa_parallel_status.fig
%      AA_PARALLEL_STATUS, by itself, creates a new AA_PARALLEL_STATUS or raises the existing
%      singleton*.
%
%      H = AA_PARALLEL_STATUS returns the handle to a new AA_PARALLEL_STATUS or the handle to
%      the existing singleton*.
%
%      AA_PARALLEL_STATUS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in AA_PARALLEL_STATUS.M with the given input arguments.
%
%      AA_PARALLEL_STATUS('Property','Value',...) creates a new AA_PARALLEL_STATUS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before aa_parallel_status_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to aa_parallel_status_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help aa_parallel_status

% Last Modified by GUIDE v2.5 26-Aug-2008 12:05:41

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @aa_parallel_status_OpeningFcn, ...
                   'gui_OutputFcn',  @aa_parallel_status_OutputFcn, ...
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


% --- Executes just before aa_parallel_status is made visible.
function aa_parallel_status_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to aa_parallel_status (see VARARGIN)

% Choose default command line output for aa_parallel_status
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes aa_parallel_status wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = aa_parallel_status_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
