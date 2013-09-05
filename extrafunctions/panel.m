
% Panel is an alternative to Matlab's "subplot" function.
% 
% INSTALLATION. To install panel, place the file "panel.m"
% on your Matlab path.
% 
% DOCUMENTATION. Scan the introductory information in the
% folder "docs". Learn to use panel by working through the
% demonstration scripts in the folder "demo". Reference
% information is available through "doc panel" or "help
% panel". For the change log, use "edit panel" to view the
% file "panel.m".



% CHANGE LOG
% 
% ############################################################
% 22/05/2011
% First Public Release Version 2.0
% ############################################################
% 
% 23/05/2011
% Incorporated an LP solver, since the one we were using
% "linprog()" is not available to users who do not have the
% Optimisation Toolbox installed.
% 
% 21/06/2011
% Added -opdf option, and changed PageSize to be equal to
% PaperPosition.
%
% 12/07/2011
% Made some linprog optimisations, inspired by "Ian" on
% Matlab central. Tested against subplot using
% demopanel2(N=9). Subplot is faster, by about 20%, but
% panel is better :). For my money, 20% isn't much of a hit
% for the extra functionality. NB: Using Jeff Stuart's
% linprog (unoptimised), panel is much slower (especially
% for large N problems); we will probably have to offer a
% faster solver at some point (optimise Jeff's?).
%
% NOTES: You will see a noticeable delay, also, on resize.
% That's the price of using physical units for the layout,
% because we have to recalculate everything when the
% physical canvas size changes. I suppose in the future, we
% could offer an option so that physical units are only used
% during export; that would make resizes fast, and the user
% may not care so much about layout on screen, if they are
% aiming for print figures. Or, you could have the ability
% to turn off auto-refresh on resize().
%
% ############################################################
% 20/07/2011
% Release Version 2.1
% ############################################################
%
% 05/10/2011
% Tidied in-file documentation (panel.m).
%
% 11/12/2011
% Added flag "no-manage-font" to constructor, as requested
% by Matlab Central user Mukhtar Ullah.
%
% ############################################################
% 13/12/2011
% Release Version 2.2
% ############################################################
%
% 21/01/2012
% Fixed bug in explicit height export option "-hX" which
% wasn't working right at all.
%
% 25/01/12
% Fixed bug in tick label display during print. _Think_ I've
% got it right, this time! Some notes below, search for
% "25/01/12".
%
% 25/01/12
% Fixed DPI bug in smoothed export figures. Bug was flagged
% up by Jesper at Matlab Central.
%
% ############################################################
% 26/01/2012
% Release Version 2.3
% ############################################################
%
% 09/03/12
% Fixed bug whereby re-positioning never got done if only
% one panel was created in an existing figure window.
%
% ############################################################
% 13/03/2012
% Release Version 2.4
% ############################################################
%
% 15/03/12
% NB: On 2008b, and possibly later versions, the fact that
% the resizeCallback() and closeCallback() are private makes
% things not work. You can fix this by removing the "Access
% = Private" modifier on that section of "methods". It works
% fine in later versions, they must have changed the access
% rules I guess.
%
% 19/07/12
% Modified so that more than one object can be managed by
% one axis. Just use p.select([h1 h2 ...]). Added function
% "getAllManagedAxes()" which returns only objects from the
% "object list" (h_object), as it now is, which represent
% axes. Suggested by Brendan Sullivan @ Matlab Central.
%
% 19/07/12
% Added support for zlabel() call (not applicable to parent
% panels, since they are implicitly 2D for axis labelling).
%
% 19/07/12
% Fixed another export bug - how did this one not get
% noticed? XLimMode (etc.) was not getting locked during
% export, so that axes without manual limits might get
% re-dimensioned during export, which is bad news. Added
% locking of limits as well as ticks, in storeAxisState().
% Hope this has no side effects!
%
% ############################################################
% 19/07/12
% Release Version 2.5
%
% NB: Owing to the introduction of management of multiple
% objects by each panel, this release should be considered
% possibly flaky. Revert to 2.4 if you have problems with
% 2.5.
% ############################################################
%
% 23/07/12
% Improved documentation for figure export in demopanelA.
%
% 24/07/12
% Added support for export to SVG, using "plot2svg" (Matlab
% Central File Exchange) as the renderer. Along the way,
% tidied the behaviour of export() a little, and improved
% reporting to the user. Changed default DPI for EPS to 600,
% since otherwise the output files are pretty shoddy, and
% the filesize is relatively unaffected.
%
% 24/07/12
% Updated documentation, particularly HTML pages and
% associated figures. Bit nicer, now.
%
% ############################################################
% 24/07/12
% Release Version 2.6
% ############################################################
%
% 22/09/12
% Added demopanelH, which illustrates how to do insets. Kudos
% to Ann Hickox for the idea.
%
% 20/03/13
% Added panel.plot() to work around poor rendering of dashed
% lines, etc. Added demopanelI to illustrate its use.
%
% 20/03/13
% Renamed setCallback to addCallback, so we can have more
% than one. Added "userdata" argument to addCallback(), and
% "event" field (and "userdata" field) to "data" passed when
% callback is fired.
%
% ############################################################
% 21/03/13
% Release Version 2.7
% ############################################################
%
% 21/03/13
% Fixed bug in panel.plot() which did not handle solid lines
% correctly.
%
% 12/04/13
% Added back setCallback() with appropriate semantics, for
% the use of legacy code (or, really, future code, these
% semantics might be useful to someone). Also added the
% function clearCallbacks().
%
% 12/04/13
% Removed panel.plot() because it just seemed to be too hard
% to manage. Instead, we'll let the user plot things in the
% usual way, but during export (when things are well under
% our control), we'll fix up any dashed lines that the user
% has requested using the call fixdash(). Thus, we apply the
% fix only where it's needed, when printing to an image
% file, and save all the faffing with resize callbacks.
%
% ############################################################
% 12/04/13
% Release Version 2.8
% ############################################################



classdef (Sealed = true) panel < handle
	
	% panel is an alternative to subplot
	%
	% for detailed reference information, use "doc panel". to
	% learn how to use panel, work through the demos (in the
	% "demo" folder), starting with "demopanel1".
	
	
	
	
	
	
	%% ---- PROPERTIES ----
	
	properties (Constant = true, Hidden = true)
		
		PANEL_TYPE_UNCOMMITTED = 0;
		PANEL_TYPE_PARENT = 1;
		PANEL_TYPE_OBJECT = 2;
		
	end
	
	properties (Constant = true)
		
		RENDER_MODE_NORMAL = 0;
		RENDER_MODE_PREPRINT = 1;
		RENDER_MODE_POSTPRINT = 2;
		
	end
	
	properties
		
		% these properties are only here for documentation. they
		% are actually stored in "prop". it's just some subsref
		% madness.
		
		% font name to use for axis text (inherited)
		%
		% access: read/write
		% default: normal
		fontname
		
		% font size to use for axis text (inherited)
		%
		% access: read/write
		% default: normal
		fontsize
		
		% font weight to use for axis text (inherited)
		%
		% access: read/write
		% default: normal
		fontweight
		
		% boolean indicating whether this panel aligns its children
		%
		% if align is true, this parent panel will force its
		% children to be aligned on their outer edges. if not,
		% this alignment will only happen if they share the same
		% margin values. align can usually be left at false, but
		% may need to be set to true for particular panels in
		% some layouts where mixed margins are used.
		%
		% NB: it's possible that there are no circumstances
		% where you would need this, so don't use it unless you
		% feel you have to. if you come up with a scenario in
		% which you definitely do need it, please let me know so
		% that i can add a demo.
		%
		% access: read/write
		% default: false
		%
		% see also: margin
		align
		
		% the units that are used when reading/writing margins
		%
		% units can be set to any of 'mm', 'cm', 'in' or 'pt'.
		% it only affects the read/write interface; values
		% stored already are not re-interpreted.
		%
		% access: read/write
		% default: mm
		units
		
		% the panel's margin vector in the form [left bottom right top]
		%
		% the margin is key to the layout process. layout
		% renders all objects as large as possible whilst not
		% violating margin constraints. margins are respected
		% between objects and between an object and the edges of
		% the figure.
		%
		% access: read/write
		% default: [12 10 2 2] (mm)
		%
		% see also: align, marginleft, marginbottom, marginright, margintop
		margin
		
		% one element of the margin vector
		%
		% access: read/write
		% default: see margin
		%
		% see also: margin
		marginleft
		
		% one element of the margin vector
		%
		% access: read/write
		% default: see margin
		%
		% see also: margin
		marginbottom
		
		% one element of the margin vector
		%
		% access: read/write
		% default: see margin
		%
		% see also: margin
		marginright
		
		% one element of the margin vector
		%
		% access: read/write
		% default: see margin
		%
		% see also: margin
		margintop
		
		% return position of panel
		%
		% return the panel's position in normalized
		% coordinates (normalized to the figure window that
		% is associated with the panel). note that parent
		% panels have positions too, even though nothing is
		% usually rendered. uncommitted panels, too.
		%
		% access: read only
		position
		
		% return handle of associated figure
		%
		% access: read only
		figure
		
		% return handle of associated axis
		%
		% if the panel is not an axis panel, empty is returned.
		% object includes axis, but axis does not include
		% object.
		%
		% access: read only
		%
		% see also: object
		axis
		
		% return handle of associated object
		%
		% if the panel is not an object panel, empty is
		% returned. object includes axis, but axis does not
		% include object.
		%
		% access: read only
		%
		% see also: axis
		object
		
		% access properties of panel's children
		%
		% if the panel is a parent panel, "children" gives
		% access to some properties of its children (direct
		% descendants). "children" can be abbreviated "ch".
		% properties that can be accessed are as follows.
		%
		% axis: read-only, returns an array
		% object: read-only, returns an array
		%
		% margin: write-only
		% fontname: write-only
		% fontsize: write-only
		% fontweight: write-only
		% align: write-only
		%
		% EXAMPLE:
		%   h = p.ch.axis;
		%   p.ch.margin = 3;
		%
		% see also: descendants
		children
		
		% access properties of panel's descendants
		%
		% if the panel is a parent panel, "descendants" gives
		% access to some properties of its descendants
		% (children, grandchildren, etc.). "descendants" can be
		% abbreviated "de". properties that can be accessed are
		% as follows.
		%
		% axis: read-only, returns an array
		% object: read-only, returns an array
		%
		% margin: write-only
		% fontname: write-only
		% fontsize: write-only
		% fontweight: write-only
		% align: write-only
		%
		% EXAMPLE:
		%   h = p.de.axis;
		%   p.de.margin = 3;
		%
		% see also: children
		descendants
		
	end
	
	properties (Access = private)
		
		% associated figure window
		h_figure
		
		% parent graphics object
		h_parent
		
		% this is empty for the root PANEL, populated for all others
		parent
		
		% this is always the root panel associated with this
		m_root
		
		% packing position within parent
		%
		% empty:            relative positioning mode (stretch)
		% scalar fraction:  relative positioning mode
		% 1x4 dimension:    absolute positioning mode
		packpos
		
		% packing dimension of children
		%
		% 1 : horizontal
		% 2 : vertical
		packdim
		
		% panel type
		m_panelType
		
		% fixdash lines
		m_fixdash
		m_fixdash_restore
		
		% associated managed graphics object (usually, an axis)
		h_object
		
		% show axis (only the root has this extra axis, if show() is active)
		h_showAxis
		
		% children (only a parent panel has non-empty, here)
		m_children
		
		% callback (any functions listed in this cell array are called when events occur)
		m_callback
		
		% local properties (actual properties is this overlaid on inherited/default properties)
		%
		% see getPropertyValue()
		prop
		
		% state
		%
		% private state information used during various processing
		state
		
		% rendering context for this panel
		%
		% this is the most recent rendering context passed to
		% the PANEL - if it is asked to renderPanel() again,
		% without being passed a new context, it uses this one.
		% thus, we can re-render some things about a child
		% without having to re-render the whole tree.
		m_context
		
	end
	
	
	
	
	
	
	
	%% ---- PUBLIC CTOR/DTOR ----
	
	methods
		
		function p = panel(varargin)
			
			% create a new panel
			%
			% p = panel(...)
			%   create a new panel. initially, the panel is an
			%   "uncommitted panel". calling pack() or select() on
			%   the panel will commit it as a "parent panel" or an
			%   "object panel", respectively. the following
			%   arguments may be passed, in any order.
			%
			% h_parent
			%   a handle to a graphics object that will act as the
			%   parent of the new panel. this is usually a figure
			%   handle, but may be a handle to any graphics
			%   object, in principle. currently, an error is
			%   raised unless it's a figure or a uipanel.
			%
			% 'add'
			%   existing panels attached to the figure should not
			%   be destroyed before the new one is created.
			%
			% 'defer'
			%   the panel should be created with render()ing
			%   deferred; it will not be rendered until you call
			%   refresh() or export() on the panel. if you create
			%   a complex layout, this will improve performance.
			%
			% 'no-manage-font'
			%   by default, the panel will manage fonts of titles
			%   and axis labels. this prevents the user from
			%   setting individual fonts on those items. pass this
			%   flag to disable font management for this panel.
			%
			% see also: panel (overview), pack(), select()

			% PRIVATE DOCUMENTATION
			%
			% h_parent
			%   can also be a handle to an existing panel. this is
			%   used internally when pack()ing child panels into a
			%   parent panel.
			
			% default condition
			passed_h_parent = [];
			add = false;
			
			% default state
			p.state = [];
			p.state.name = '';
			p.state.defer = 0;
			p.state.manage_font = 1;
			
			% peel off args
			while ~isempty(varargin)
				
				% get arg
				arg = varargin{1};
				varargin = varargin(2:end);
				
				% handle text
				if ischar(arg)
					
					switch arg
						
						case 'add'
							add = true;
							continue;
							
						case 'defer'
							p.state.defer = 1;
							continue;
							
						case 'no-manage-font'
							p.state.manage_font = 0;
							continue;
							
						otherwise
							error('panel:InvalidArgument', ['unrecognised text argument "' arg '"']);
							
					end
					
				end
				
				% handle parent
				if isscalar(arg) && (ishandle(arg) || isa(arg, 'panel'))
					passed_h_parent = arg;
					continue;
				end
				
				% error
				error('panel:InvalidArgument', 'unrecognised argument to panel constructor');
				
			end
			
			% no callbacks
			p.m_callback = {};
			
			% no fixdash
			p.m_fixdash = {};
			
			% debug output
			panel.debugmsg('creating new panel...');
			
			% attach to current figure if no parent supplied
			if nargin < 1 || isempty(passed_h_parent)
				passed_h_parent = gcf;
				
				% this might cause a figure to be created - if so,
				% give it time to display now so we don't get a (or
				% two, in fact!) resize event(s) later
				drawnow
			end
			
			% see what our parent is
			if ishandle(passed_h_parent)
				
				parentType = get(passed_h_parent, 'type');
				
				% parent is a graphics object - become a root panel
				p.state.name = 'root';
				p.parent = [];
				p.m_root = p;
				
				switch parentType
					
					case 'uipanel'
						p.h_parent = passed_h_parent;
						p.h_figure = getParentFigure(passed_h_parent);
						
					case 'figure'
						p.h_parent = passed_h_parent;
						p.h_figure = passed_h_parent;
						
					otherwise
						error('panel:InvalidArgument', ['panel() cannot be attached to an object of type "' parentType '"']);
						
				end
				
			elseif isa(passed_h_parent, 'panel')
				
				% parent is a panel - become its child
				indexInParent = int2str(length(passed_h_parent.m_children)+1);
				if passed_h_parent.isRoot()
					p.state.name = ['(' indexInParent ')'];
				else
					p.state.name = [passed_h_parent.state.name(1:end-1) ',' indexInParent ')'];
				end
				p.h_parent = passed_h_parent.h_parent;
				p.h_figure = passed_h_parent.h_figure;
				p.parent = passed_h_parent;
				p.m_root = passed_h_parent.m_root;
				
			else
				
				% error
				error('panel:InvalidArgument', 'argument to panel() constructor must be a figure handle or a panel object');
				
			end
			
			% default state
			p.packpos = [];
			p.packdim = 2;
			p.m_panelType = p.PANEL_TYPE_UNCOMMITTED;
			p.prop = panel.getPropertyInitialState();
			
			% if we are root
			if p.isRoot()
				
				% lay in callbacks
				addHandleCallback(p.h_figure, 'CloseRequestFcn', @panel.closeCallback);
				addHandleCallback(p.h_parent, 'ResizeFcn', @panel.resizeCallback);
				
				% register for callbacks
				if add
					panel.callbackDispatcher('registerNoClear', p);
				else
					panel.callbackDispatcher('register', p);
				end
				
				% lock class in memory (prevent persistent from being cleared by clear all)
				panel.lockClass();
				
			end
			
			% debug output
			panel.debugmsg(['created "' p.state.name '"!']);
			
		end
		
		function delete(p)
			
			% destroy a panel
			%
			% delete(p)
			%   destroy the passed panel, deleting all associated
			%   graphics objects.
			%
			% NB: you won't usually have to call this explicitly.
			% it is called automatically for all attached panels
			% when you close the associated figure.
			
			% debug output
			panel.debugmsg(['deleting "' p.state.name '"...']);
			
			% delete managed graphics objects
			for n = 1:length(p.h_object)
				h = p.h_object(n);
				if ishandle(h)
					delete(h);
				end
			end
			
			% delete associated show axis
			if ~isempty(p.h_showAxis) && ishandle(p.h_showAxis)
				delete(p.h_showAxis);
			end
			
			% delete all children (child will remove itself from
			% "m_children" on delete())
			while ~isempty(p.m_children)
				delete(p.m_children(end));
			end
			
			% unregister...
			if p.isRoot()
				
				% ...for callbacks
				panel.callbackDispatcher('unregister', p);
				
			else
				
				% ...from parent
				p.parent.removeChild(p);
				
			end
			
			% debug output
			panel.debugmsg(['deleted "' p.state.name '"!']);
			
		end
		
	end
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	%% ---- PUBLIC DISPLAY ----

	methods (Hidden = true)
		
		function disp(p)
			
			display(p);
			
		end
		
		function display(p, indent)

			% default
			if nargin < 2
				indent = '';
			end
			
			% handle non-scalar (should not exist!)
			nels = numel(p);
			if nels > 1
				sz = size(p);
				sz = sprintf('%dx', sz);
				disp([sz(1:end-1) ' array of panel objects']);
				return
			end
			
			% header
			header = indent;
			if p.isObject()
				header = [header 'Object ' p.state.name ': '];
			elseif p.isParent()
				header = [header 'Parent ' p.state.name ': '];
			else
				header = [header 'Uncommitted ' p.state.name ': '];
			end
			if p.isRoot()
				pp = ['attached to Figure ' num2str(p.h_figure)];
			else
				if isempty(p.packpos)
					pp = 'stretch';
				else
					pp = sprintf('%.3f ', p.packpos);
					pp = pp(1:end-1);
				end
			end
			header = [header '[' pp];
			if p.isParent()
				edges = {'hori' 'vert'};
				header = [header ', ' edges{p.packdim}];
			end
			header = [header ']'];

			% margin
			header = rpad(header, 60);
			header = [header '[ margin ' sprintf('%.3g ', p.getPropertyValue('margin')) p.getPropertyValue('units') ']'];
			
% 			% index
% 			if isfield(p.state, 'index')
% 				header = [header ' (' int2str(p.state.index) ')'];
% 			end

			% display
			disp(header);
			
			% children
			for c = 1:length(p.m_children)
				p.m_children(c).display([indent '  ']);
			end
						
		end
			
	end
	
	
	
	
	
	
	
	
	
	
	%% ---- PUBLIC METHODS ----

	methods
		
		function xlabel(p, text)
			
			% apply an xlabel to the panel (or group)
			%
			% p.xlabel(...)
			%   behaves just like xlabel() at the prompt (you can
			%   use that as an alternative) when called on an axis
			%   panel. when called on a parent panel, however, the
			%   group of objects within that parent have a label
			%   applied. when called on a non-axis object panel,
			%   an error is raised.
			
			h = get(p.getOrCreateAxis(), 'xlabel');
			set(h, 'string', text);
			if p.isParent()
				set(h, 'visible', 'on');
			end
			
		end
		
		function ylabel(p, text)
			
			% apply a ylabel to the panel (or group)
			%
			% p.ylabel(...)
			%   behaves just like ylabel() at the prompt (you can
			%   use that as an alternative) when called on an axis
			%   panel. when called on a parent panel, however, the
			%   group of objects within that parent have a label
			%   applied. when called on a non-axis object panel,
			%   an error is raised.
			
			h = get(p.getOrCreateAxis(), 'ylabel');
			set(h, 'string', text);
			if p.isParent()
				set(h, 'visible', 'on');
			end
			
		end
		
		function zlabel(p, text)
			
			% apply a zlabel to the panel (or group)
			%
			% p.zlabel(...)
			%   behaves just like zlabel() at the prompt (you can
			%   use that as an alternative) when called on an axis
			%   panel. when called on a parent panel, however,
			%   this method raises an error, since parent panels
			%   are assumed to be 2D, with respect to axes.
			
			if p.isParent()
				error('panel:ZLabelOnParentAxis', 'can only call zlabel() on an object panel');
			end
			
			h = get(p.getOrCreateAxis(), 'zlabel');
			set(h, 'string', text);
			
		end
		
		function title(p, text)
			
			% apply a title to the panel (or group)
			%
			% p.title(...)
			%   behaves just like title() at the prompt (you can
			%   use that as an alternative) when called on an axis
			%   panel. when called on a parent panel, however, the
			%   group of objects within that parent have a title
			%   applied. when called on a non-axis object panel,
			%   an error is raised.
			
			h = title(p.getOrCreateAxis(), text);
			if p.isParent()
				set(h, 'visible', 'on');
			end
			
		end
		
% 		function p = xlabels(p, text)
% 			
% 			if p.isParent()
% 				for c = 1:length(p.m_children)
% 					p.m_children(c).xlabels(text);
% 				end
% 			elseif p.isObject()
% 				p.xlabel(text);
% 			end
% 			
% 		end
% 		
% 		function p = ylabels(p, text)
% 			
% 			if p.isParent()
% 				for c = 1:length(p.m_children)
% 					p.m_children(c).ylabels(text);
% 				end
% 			elseif p.isObject()
% 				p.ylabel(text);
% 			end
% 			
% 		end
		
		function hold(p, state)
			
			% set the hold on/off state of the associated axis
			% 
			% p.hold('on' / 'off')
			%   behaves just like the matlab "hold on", "hold off"
			%   command. however, matlab's "hold off" sets the
			%   axis state to a state that is unsuitable for
			%   panel, so you should use this function provided by
			%   panel instead, if you want your panel properties
			%   (e.g. fontname) to be rendered correctly.
			%
			%   NB: if you prefer you can use the Matlab "hold"
			%   command, and call refresh() on the root panel when
			%   you've finished setting properties to cause fonts
			%   to render correctly.

			% because the matlab "hold off" command sets an axis's
			% nextplot state to "replace", we lose control over
			% the axis properties (such as fontname). we set
			% nextplot to "replacechildren" when we create an
			% axis, but if the user does a "hold on, hold off"
			% cycle, we lose that. therefore, we offer this
			% alternative.
			
			% check
			if ~p.isObject()
				error('panel:HoldWhenNotObjectPanel', 'can only call hold() on an object panel');
			end
			
			% check
			h_axes = p.getAllManagedAxes();
			if isempty(h_axes)
				error('panel:HoldWhenNoAxes', 'can only call hold() on a panel that manages one or more axes');
			end
			
			% switch
			switch state
				case {'on' true 1}
					set(h_axes, 'nextplot', 'add');
				case {'off' false 0}
					set(h_axes, 'nextplot', 'replacechildren');
				otherwise
					error('panel:InvalidArgument', 'argument to hold() must be ''on'', ''off'', or boolean');
			end
			
		end
		
		function fixdash(p, hs, linestyle)
			
			% pass dashed lines to be fixed up during export
			%
			% p.fixdash(h, linestyle)
			%   add the lines specified as handles in "h" to the
			%   list of lines to be "fixed up" during export.
			%   panel will attempt to get the lines to look right
			%   during export to all formats where they would
			%   usually get mussed up. see demopanelI for an
			%   example of how it works.
			%
			%   the above is the usual usage of fixdash(), but
			%   you can get more control over linestyle by
			%   specifying the additional argument, "linestyle".
			%   if "linestyle" is supplied, it is used as the
			%   linestyle; if not, the current linestyle of the
			%   line (-, --, -., :) is used. "linestyle" can
			%   either be a text string or a series of numbers, as
			%   described below.
			%
			%     '-' solid
			%     '--' dashed, equal to [2 0.75]
			%     '-.' dash-dot, equal to [2 0.75 0.5 0.75]
			%     ':', '.' dotted, equal to [0.5 0.5]
			%
			%   a number series should be 1xN, where N is a
			%   multiple of 2, as in the examples above, and
			%   specifies the lengths of any number of dash
			%   components that are used before being repeated.
			%   for instance, '-.' generates a 2 unit segment
			%   (dash), a 0.75 unit gap, then a 0.5 unit segment
			%   (dot) and a final 0.75 unit gap. at present, the
			%   units are always millimetres. this system is
			%   extensible, so that the following examples are
			%   also valid:
			%
			%     '--..' dash-dash-dot-dot
			%     '-..-.' dash-dot-dot-dash-dot
			%     [2 1 4 1 6 1] 2 dash, 4 dash, 6 dash

			% default
			if nargin < 3
				linestyle = [];
			end
			
			% bubble up to root
			if ~p.isRoot()
				p.m_root.fixdash(hs, linestyle);
				return
			end
			
			% for each passed handle
			for h = (hs(:)')
				
				% check it's still a handle
				if ~ishandle(h)
					continue
				end
				
				% check it's a line
				if ~isequal(get(h, 'type'), 'line')
					continue
				end
				
				% update if in list
				found = false;
				for i = 1:length(p.m_fixdash)
					if h == p.m_fixdash{i}.h
						p.m_fixdash{i}.linestyle = linestyle;
						found = true;
						break
					end
				end
				
				% else add to list
				if ~found
					p.m_fixdash{end+1} = struct('h', h, 'linestyle', linestyle);
				end
				
			end
			
		end
		
		function show(p)
			
			% highlight the outline of the panel
			%
			% p.show()
			%   make the outline of the panel "p" show up in red
			%   in the figure window. this is useful for
			%   understanding a complex layout.
			%
			% see also: identify()

			r = p.getObjectPosition();
			
			if ~isempty(r)
				h = p.getShowAxis();
				delete(get(h, 'children'));
				xdata = [r(1) r(1)+r(3) r(1)+r(3) r(1) r(1)];
				ydata = [r(2) r(2) r(2)+r(4) r(2)+r(4) r(2)];
				plot(h, xdata, ydata, 'r-', 'linewidth', 5);
				axis([0 1 0 1])
			end
			
		end
		
		function export(p, varargin)
			
			% to export the root panel to an image file
			%
			% p.export(...)
			%
			% export the figure containing panel "p" to an image
			% file. you must supply at least the filename to
			% export to (you can, optionally, omit the file
			% extension). further arguments must be option
			% strings starting with the '-' character. "p" should
			% be the root panel (the first panel that was
			% created).
			%
			% if you are targeting a print publication, you may
			% find it easiest to size your output using the "paper
			% model". if you prefer, you can use the "explicit
			% sizing model", instead. these two sizing models are
			% described below. underneath these are listed the
			% options unrelated to sizing (which apply regardless
			% of which sizing model you use).
			%
			% NB: if you pass 'defer' to the constructor, calling
			% export() both exports the panel and releases the
			% defer mode. future changes to properties (e.g.
			% margins) will be rendered immediately.
			%
			%
			%
			% PAPER SIZING MODEL:
			%
			% using the paper sizing model, you specify your
			% target as a region of a piece of paper, and the
			% actual size in millimeters is calculated for you.
			% this is usually very convenient, but if you find it
			% unsuitable, the explicit sizing model (next section)
			% is provided as an alternative.
			%
			% to specify the region, you specify the type (size)
			% of paper, the orientation, the number of columns,
			% and the aspect ratio of the figure (or the fraction
			% of a column to fill). usually, the remaining options
			% can be left as defaults.
			%
			% -pX
			%   X is the paper type, A2-A6, letter (default is
			%   A4).
			%
			% -l
			%   specify landscape mode (default is portrait).
			%
			% -mX
			%   X is the paper margins in mm. you can provide a
			%   scalar (same margins all round) or a
			%   comma-separated list of four values, specifying
			%   the left, bottom, right, top margins separately
			%   (default is 20mm all round).
			%
			% -iX
			%   X is the inter-column space in mm (default is
			%   5mm).
			%
			% -cX
			%   X is the number of columns (default is 1).
			%
			% NB: the following two options represent two ways to
			% specify the height of the figure relative to the
			% space defined by the above options. if you supply
			% both, whichever comes second will be used.
			%
			% -aX
			%   X is the aspect ratio of the resulting image file.
			%   X can be one of the following strings - s
			%   (square), g (landscape golden ratio), gp (portrait
			%   golden ratio), h (half-height), d (double-height)
			%   - or a number greater than zero, to specify the
			%   aspect ratio explicitly. note that, if using the
			%   numeric form, the ratio is expressed as the
			%   quotient of width over height, in the usual way.
			%   ratios greater than 10 or less than 0.1 are
			%   disallowed, since these can cause a very large
			%   figure file to be created accidentally. default is
			%   to use the landscape golden ratio.
			%
			% -fX
			%   X is the fraction of the column (or page, if there
			%   are not columns) to fill. X can be one of the
			%   following strings - a (all), tt (two thirds), h
			%   (half), t (third), q (quarter) - or a fraction
			%   between 0 and 1, to specify the fraction of the
			%   space to fill as a number. default is to use
			%   aspect ratio, not fill fraction.
			%
			%
			%
			% EXPLICIT SIZING MODEL:
			%
			% these two options override any output of the paper
			% model, so you can override just one, or both (in
			% which case all paper model options are ignored).
			%
			% -wX
			%   X is explicit width in mm (default is to use the
			%   width produced by the paper model).
			%
			% -hX
			%   X is explicit height in mm (default is to use the
			%   height produced by the paper model).
			%
			%
			%
			% NON-SIZING OPTIONS:
			%
			% finally, a few options are provided to control how
			% the prepared figure is exported. note that DPI below
			% 150 is only recommended for sizing drafts, since
			% font and line sizes are not rendered even vaguely
			% accurately in some cases. at the other end, DPI
			% above 600 is unlikely to be useful except when
			% submitting camera-ready copy.
			%
			% -rX
			%   X is the resolution (DPI) at which to produce the
			%   output file. X can be one of the following strings
			%   - d (draft, 75DPI), n (normal, 150DPI), h (high,
			%   300DPI), p (publication quality, 600DPI), x
			%   (extremely high quality, 1200DPI) - or just
			%   the DPI as a number (must be in 75-2400). the
			%   default depends on the output format (see below).
			%
			% -rX/S
			%   X is the DPI, S is the smoothing factor, which can
			%   be 2 or 4. the output file is produced at S times
			%   the specified DPI, and then reduced in size to the
			%   specified DPI by averaging. thus, the hard edges
			%   produced by the renderer are smoothed - the effect
			%   is much like "anti-aliasing".
			%
			% NB: the DPI setting might be expected to have no
			% effect on vector formats. this is true for SVG, but
			% not for EPS, where the DPI affects the numerical
			% precision used as well as the size of some image
			% elements, but has little effect on file size. for
			% this reason, the default DPI is 150 for bitmap
			% formats but 600 for vector formats.
			%
			% -s
			%   print sideways (default is to print upright)
			%
			% -oX
			%   X is the output format - choose from most of the
			%   built-in image device drivers supported by "print"
			%   (try "help print"). this includes "png", "jpg",
			%   "tif", "eps" and "pdf". note that "eps"/"ps"
			%   resolve to "epsc2"/"psc2", for convenience. to use
			%   the "eps"/"ps" devices, use "-oeps!"/"-ops!". you
			%   may also specify "svg", if you have the tool
			%   "plot2svg" on your path (available at Matlab
			%   Central File Exchange).
			%
			%
			%
			% EXAMPLES:
			%
			% default export of 'myfig', creates 'myfig.png' at a
			% size of 170x105mm (1004x620px). this size comes
			% from: A4 (210mm wide), minus two 20mm margins
			% (170mm), and using the golden aspect ratio to give a
			% height of 105mm, and finally 150DPI to give the
			% pixel size.
			%
			% p.export('myfig')
			%
			% when producing the final camera-ready image for a
			% square figure that will sit in one of the two
			% columns of a letter-size paper journal with default
			% margins and inter-column space, we might use this:
			%
			% p.export('myfig', '-pletter', '-c2', '-as', '-rp');
			
			% check
			if ~p.isRoot()
				error('panel:ExportWhenNotRoot', 'cannot export() this panel - it is not the root panel');
			end
			
			% parameters
			pars = [];
			pars.filename = '';
			pars.fmt = 'png';
			pars.ext = 'png';
			pars.dpi = [];
			pars.smooth = 1;
			pars.paper = 'A4';
			pars.landscape = false;
			pars.fill = -1.618;
			pars.cols = 1;
			pars.intercolumnspacing = 5;
			pars.margin = 20;
			pars.sideways = false;
			pars.width = 0;
			pars.height = 0;
			invalid = false;
			
			% interpret args
			for a = 1:length(varargin)
				
				% extract
				arg = varargin{a};
				
				% all arguments must be non-empty strings
				if ~isstring(arg)
					error('panel:InvalidArgument', ...
						'all arguments to export() must be non-empty strings');
				end
				
				% is filename?
				if arg(1) ~= '-'
					
					% error if already set
					if ~isempty(pars.filename)
						error('panel:InvalidArgument', ...
							['at argument "' arg '", the filename is already set ("' pars.filename '")']);
					end
					
					% ok, continue
					pars.filename = arg;
					continue
					
				end

				% split off option key and option value
				if length(arg) < 2
					error('panel:InvalidArgument', ...
						['at argument "' arg '", no option specified']);
				end
				key = arg(2);
				val = arg(3:end);
				
				% switch on option key
				switch key

					case 'p'
						pars.paper = validate_par(val, arg, {'A2' 'A3' 'A4' 'A5' 'A6' 'letter'});

					case 'l'
						pars.landscape = true;
						validate_par(val, arg, 'empty');

					case 'm'
						pars.margin = validate_par(str2num(val), arg, 'dimension', 'nonneg');

					case 'i'
						pars.intercolumnspacing = validate_par(str2num(val), arg, 'scalar', 'nonneg');

					case 'c'
						pars.cols = validate_par(str2num(val), arg, 'scalar', 'integer');

					case 'f'
						switch val
							case 'a', pars.fill = 1;      % all
							case 'w', pars.fill = 1;      % whole (legacy, not documented)
							case 'tt', pars.fill = 2/3;   % two thirds
							case 'h', pars.fill = 1/2;    % half
							case 't', pars.fill = 1/3;    % third
							case 'q', pars.fill = 1/4;    % quarter
							otherwise
								pars.fill = validate_par(str2num(val), arg, 'scalar', [0 1]);
						end

					case 'a'
						switch val
							case 's', pars.fill = -1;         % square
							case 'g', pars.fill = -1.618;     % golden ratio (landscape)
							case 'gp', pars.fill = -1/1.618;  % golden ratio (portrait)
							case 'h', pars.fill = -2;         % half height
							case 'd', pars.fill = -0.5;       % double height
							otherwise
								pars.fill = -validate_par(str2num(val), arg, 'scalar', [0.1 10]);
						end

					case 'w'
						pars.width = validate_par(str2num(val), arg, 'scalar', 'nonneg', [10 Inf]);

					case 'h'
						pars.height = validate_par(str2num(val), arg, 'scalar', 'nonneg', [10 Inf]);

					case 'r'
						% peel off smoothing ("/...")
						if any(val == '/')
							f = find(val == '/', 1);
							switch val(f+1:end)
								case '2', pars.smooth = 2;
								case '4', pars.smooth = 4;
								otherwise, error('panel:InvalidArgument', ...
										['invalid argument "' arg '", part after / must be "2" or "4"']);
							end
							val = val(1:end-2);
						end

						switch val
							case 'd', pars.dpi = 75;      % draft
							case 'n', pars.dpi = 150;     % normal
							case 'h', pars.dpi = 300;     % high
							case 'p', pars.dpi = 600;     % publication quality
							case 'x', pars.dpi = 1200;    % extremely high quality
							otherwise
								pars.dpi = validate_par(str2num(val), arg, 'scalar', [75 2400]);
						end

					case 's'
						pars.sideways = true;
						validate_par(val, arg, 'empty');

					case 'o'
						fmts = {
							'png' 'png' 'png'
							'tif' 'tiff' 'tif'
							'tiff' 'tiff' 'tif'
							'jpg' 'jpeg' 'jpg'
							'jpeg' 'jpeg' 'jpg'
							'ps' 'psc2' 'ps'
							'ps!' 'psc' 'ps'
							'psc' 'psc' 'ps'
							'ps2' 'ps2' 'ps'
							'psc2' 'psc2' 'ps'
							'eps' 'epsc2' 'eps'
							'eps!' 'eps' 'eps'
							'epsc' 'epsc' 'eps'
							'eps2' 'eps2' 'eps'
							'epsc2' 'epsc2' 'eps'
							'pdf' 'pdf' 'pdf'
							'svg' 'svg' 'svg'
							};
						validate_par(val, arg, fmts(:, 1)');
						index = isin(fmts(:, 1), val);
						pars.fmt = fmts{index, 2};
						pars.ext = fmts{index, 3};

					otherwise
						error('panel:InvalidArgument', ...
							['invalid argument "' argtext '", option is not recognised']);

				end
				
			end
			
			% extract
			is_bitmap = ismember(pars.fmt, {'png' 'jpeg' 'tiff'});
			
			% default DPI
			if isempty(pars.dpi)
				if is_bitmap
					pars.dpi = 150;
				else
					pars.dpi = 600;
				end
			end

			% validate
			if isequal(pars.fmt, 'svg') && isempty(which('plot2svg'))
				error('panel:Plot2SVGMissing', 'export to SVG requires plot2svg (Matlab Central File Exchange)');
			end
			
			% validate
			if ~is_bitmap && pars.smooth ~= 1
				pars.smooth = 1;
				warning('panel:NoSmoothVectorFormat', 'requested smoothing will not be performed (chosen export format is not a bitmap format)');
			end
			
			% validate
			if isempty(pars.filename)
				error('panel:InvalidArgument', 'filename not supplied');
			end
			
			% make sure filename has extension
			if ~any(pars.filename == '.')
				pars.filename = [pars.filename '.' pars.ext];
			end
			
			
			
%%%% GET TARGET DIMENSIONS (BEGIN)
			
			% get space for figure
			switch pars.paper
				case 'A0', sz = [841 1189];
				case 'A1', sz = [594 841];
				case 'A2', sz = [420 594];
				case 'A3', sz = [297 420];
				case 'A4', sz = [210 297];
				case 'A5', sz = [148 210];
				case 'A6', sz = [105 148];
				case 'letter', sz = [216 279];
				otherwise
					error(['unrecognised paper size "' pars.paper '"'])
			end
			
			% orientation of paper
			if pars.landscape
				sz = sz([2 1]);
			end
			
			% paper margins (scalar or quad)
			if isscalar(pars.margin)
				pars.margin = pars.margin * [1 1 1 1];
			end
			sz = sz - pars.margin(1:2) - pars.margin(3:4);
			
			% divide by columns
			w = (sz(1) + pars.intercolumnspacing) / pars.cols - pars.intercolumnspacing;
			sz(1) = w;
			
			% explicit measurement overrides automatic
			if pars.width
				sz(1) = pars.width;
			end
			
			% apply fill / aspect ratio
			if pars.fill > 0
				% fill fraction
				sz(2) = sz(2) * pars.fill;
			elseif pars.fill < 0
				% aspect ratio
				sz(2) = sz(1) * (-1 / pars.fill);
			end
			
			% explicit measurement overrides automatic
			if pars.height
				sz(2) = pars.height;
			end
			
%%%% GET TARGET DIMENSIONS (END)

			
			
			% orientation of figure is upright, unless printing
			% sideways, in which case the printing space is rotated too
			if pars.sideways
				set(p.h_figure, 'PaperOrientation', 'landscape')
				sz = fliplr(sz);
			else
				set(p.h_figure, 'PaperOrientation', 'portrait')
			end
			
			% report export size
			msg = ['exporting to ' int2str(sz(1)) 'x' int2str(sz(2)) 'mm'];
			if is_bitmap
				psz = sz / 25.4 * pars.dpi;
				msg = [msg ' (' int2str(psz(1)) 'x' int2str(psz(2)) 'px @ ' int2str(pars.dpi) 'DPI)'];
			else
				msg = [msg ' (vector format @ ' int2str(pars.dpi) 'DPI)'];
			end
			disp(msg);
			
			% if we are in defer state, we need to do a clean
			% render first so that axes get positioned so that
			% axis ticks get set correctly (if they are in
			% automatic mode), since the RENDER_MODE_PREPRINT
			% render will store the tick states.
			if p.state.defer
				p.state.defer = 0;
				p.renderAll();
			end

			% enable rendering
			p.state.defer = 0;
			
			% do a pre-print render
			context.size_in_mm = sz;
			context.mode = panel.RENDER_MODE_PREPRINT;
			p.renderAll(context);
			
			% need also to disable the warning that we should set
			% PaperPositionMode to auto during this operation -
			% we're setting it explicitly.
			w = warning('off', 'MATLAB:Print:CustomResizeFcnInPrint');
			
			% handle smoothing
			pars.write_dpi = pars.dpi;
			if pars.smooth > 1
				pars.write_dpi = pars.write_dpi * pars.smooth;
				print_filename = [pars.filename '-temp'];
			else
				print_filename = pars.filename;
			end

			% disable rendering so we don't get automatic
			% rendering during any figure resize operations.
			p.state.defer = 1;
			
			% set size of figure now. it's important we do this
			% after the pre-print render, because in SVG export
			% mode the on-screen figure size is changed and that
			% would otherwise affect ticks and limits.
			switch pars.fmt
				
				case 'svg'
					% plot2svg (our current SVG export mechanism) uses
					% 'Units' and 'Position' (i.e. on-screen position)
					% rather than the Paper- prefixed ones used by the
					% Matlab export functions.
					
					% store old on-screen position
					svg_units = get(p.h_figure, 'Units');
					svg_pos = get(p.h_figure, 'Position');
					
					% update on-screen position
					set(p.h_figure, 'Units', 'centimeters');
					pos = get(p.h_figure, 'Position');
					pos(3:4) = sz / 10;
					set(p.h_figure, 'Position', pos);
					
				otherwise
					set(p.h_figure, ...
						'PaperUnits', 'centimeters', ...
						'PaperPosition', [0 0 sz] / 10, ...
						'PaperSize', sz / 10 ... % * 1.5 / 10 ... % CHANGED 21/06/2011 so that -opdf works correctly - why was this * 1.5, anyway? presumably was spurious...
						);
					
			end
			
			% do fixdash (not for SVG, since plot2svg does a nice
			% job of dashed lines without our meddling...)
			if ~isequal(pars.fmt, 'svg')
				p.do_fixdash(context);
			end
			
			% do the export
			switch pars.fmt
				case 'svg'
					plot2svg(print_filename, p.h_figure);
				otherwise
					print(p.h_figure, '-loose', ['-d' pars.fmt], ['-r' int2str(pars.write_dpi)], print_filename)
			end

			% undo fixdash
			if ~isequal(pars.fmt, 'svg')
				p.do_fixdash([]);
			end
			
			% set on-screen figure size back to what it was, if it
			% was changed.
			switch pars.fmt
				case 'svg'
					set(p.h_figure, 'Units', svg_units);
					set(p.h_figure, 'Position', svg_pos);
			end
			
			% enable rendering
			p.state.defer = 0;
			
			% enable warnings
			warning(w);
			
			% do a post-print render
			context.size_in_mm = [];
			context.mode = panel.RENDER_MODE_POSTPRINT;
			p.renderAll(context);
			
			% handle smoothing
			if pars.smooth > 1
				psz = sz * pars.smooth / 25.4 * pars.dpi;
				msg = [' (reducing from ' int2str(psz(1)) 'x' int2str(psz(2)) 'px)'];
				disp(['smoothing by factor ' int2str(pars.smooth) msg]);
				im1 = imread(print_filename);
				delete(print_filename);
				sz = size(im1);
				sz = [sz(1)-mod(sz(1),pars.smooth) sz(2)-mod(sz(2),pars.smooth)] / pars.smooth;
				im = zeros(sz(1), sz(2), 3);
				mm = 1:pars.smooth:(sz(1) * pars.smooth);
				nn = 1:pars.smooth:(sz(2) * pars.smooth);
				for m = 0:pars.smooth-1
					for n = 0:pars.smooth-1
						im = im + double(im1(mm+m, nn+n, :));
					end
				end
				im = uint8(im / (pars.smooth^2));
				
				% set the DPI correctly in the new file
				switch pars.fmt
					case 'png'
						dpm = pars.dpi / 25.4 * 1000; 
						imwrite(im, pars.filename, ... 
							'XResolution', dpm, ... 
							'YResolution', dpm, ... 
							'ResolutionUnit', 'meter');
					case 'tiff'
						imwrite(im, pars.filename, ... 
							'Resolution', pars.dpi * [1 1]);
					otherwise
						imwrite(im, pars.filename);
				end
			end
			
		end

		function clearCallbacks(p)
			
			% clear all callback functions for the panel
			%
			% p.clearCallbacks()
			p.m_callback = {};
			
		end
		
		function setCallback(p, func, userdata)
			
			% set the callback function for the panel
			%
			% p.setCallback(myCallbackFunction, userdata)
			%
			% NB: this function clears all current callbacks, then
			%   calls addCallback(myCallbackFunction, userdata).
			p.clearCallbacks();
			p.addCallback(func, userdata);
			
		end
		
		function addCallback(p, func, userdata)
			
			% attach a callback function to the resize event
			%
			% p.addCallback(myCallbackFunction, userdata)
			%   register myCallbackFunction() to be called when
			%   the panel is updated (usually, resized).
			%   myCallbackFunction() should accept one argument,
			%   "data", which will have the following fields.
			%
			% "userdata": the userdata passed to this function, if
			%     any was supplied, else empty.
			%
			% "panel": a reference to the panel on which the
			%     callback was set. this object can be queried in
			%     the usual way.
			%
			% "event": name of event (currently only
			%	    "render-complete").
			%
			% "context": the rendering context for the panel.
			
			invalid = ~isscalar(func) || ~isa(func, 'function_handle');
			if invalid
				error('panel:InvalidArgument', 'argument to callback() must be a function handle');
			end
			if nargin == 2
				p.m_callback{end+1} = {func []};
			else
				p.m_callback{end+1} = {func userdata};
			end
			
		end
		
		function identify(p)

			% add annotations to help identify individual panels
			%
			% p.identify()
			%   when creating a complex layout, it can become
			%   confusing as to which panel is which. this
			%   function adds a text label to each axis panel
			%   indicating how to reference the axis panel through
			%   the root panel. for instance, if "(2, 3)" is
			%   indicated, you can find that panel at p(2, 3).
			%
			% see also: show()
			
			if p.isObject()
				
				% get managed axes
				h_axes = p.getAllManagedAxes();
			
				% if no axes, ignore
				if isempty(h_axes)
					return
				end
				
				% mark first axis
				h_axes = h_axes(1);
				cla(h_axes);
				text(0.5, 0.5, p.state.name, 'fontsize', 12, 'hori', 'center', 'parent', h_axes);
				axis(h_axes, [0 1 0 1]);
				grid(h_axes, 'off')

			else
				
				% recurse
				for c = 1:length(p.m_children)
					p.m_children(c).identify();
				end
				
			end
			
		end
		
		function repack(p, packpos)
			
			% set a new packing position for an existing panel
			%
			% p.repack(packpos)
			%   packpos can be anything you would pass to pack().
			%   however, a panel cannot have its positioning mode
			%   changed by repack, so if you originally packed the
			%   panel in relative mode, packpos must be a scalar,
			%   or if you originally packed the panel in absolute
			%   mode, packpos must be a 1x4 row vector.
			
			% must match current packing mode if there is more
			% than one panel packed (if not, it's free for all)
			if p.isRoot()
				
				% root can only accept absolute positioning
				if ~isofsize(packpos, [1 4])
					error('panel:InvalidArgument', 'root panel can only use absolute positioning mode');
				end
				
			else
				
				% do we have siblings?
				nsiblings = numel(p.parent.m_children);
				
				% check
				if nsiblings > 1
					if isofsize(packpos, [1 4])
						if ~isofsize(p.packpos, [1 4])
							error('panel:InvalidArgument', 'repack() cannot change the packing mode - this panel''s siblings use relative positioning');
						end
					elseif (isscalar(packpos) && (packpos == -1 || (packpos > 0 && packpos <= 100)))
						if ~isscalar(p.packpos) && ~isempty(p.packpos)
							error('panel:InvalidArgument', 'repack() cannot change the packing mode - this panel''s siblings use absolute positioning');
						end
					else
						error('panel:InvalidArgument', 'repack() only accepts -1 (stretch), positive scalar <= 1 (relative positioning), or 1x4 (absolute positioning)');
					end
				end
				
			end
			
			% update the packpos
			p.packpos = packpos;
			
			% and renderAll
			p.renderAll();
			
		end
		
		function pack(p, varargin)
			
			% add (pack) child panel(s) into an existing panel
			%
			% p.pack(...)
			%   add children to the panel "p" and, in doing so,
			%   commit the panel as a "parent panel" (if it is not
			%   already committed). after this, it can no longer
			%   be associated with an object (though it can still
			%   have xlabel/ylabel and title). newly created child
			%   panels begin as "uncommitted panels". pack can be
			%   followed by any number of arguments, drawn from
			%   the following. 
			%
			% 'h'/'v'
			%   switch to horizontal or vertical packing
			%   direction.
			%
			% small integer (1 to 32)
			%   pack this many panels along the packing direction.
			%
			% 1xN row vector (without 'abs')
			%   pack N new panels along the packing dimension,
			%   with the relative size of each given by the
			%   elements of the vector. -1 can be passed for any
			%   elements to mark those panel as 'stretchable', so
			%   that they fill available space left over by other
			%   panels packed alongside. the sum of the vector
			%   (apart from any -1 entries) should not come to
			%   more than 1, or a warning will be generated during
			%   rendering. an example would be [1/4 1/4 -1], to
			%   pack 3 panels, at 25, 25 and 50% relative sizes.
			%   though, NB, you can use percentages instead of
			%   fractions if you prefer, in which case they should
			%   not sum to over 100. so that same pack() would be
			%   [25 25 -1].
			%
			% 'abs'
			%   the next argument will be an absolute position, as
			%   described below. you should avoid using absolute
			%   positioning mode, in general, since this does not
			%   take advantage of panel's automatic layout.
			%   however, on occasion, you may need to take manual
			%   control of the position of one or more panels.
			%
			% 1x4 row vector (after 'abs')
			%   pack 1 new panel using absolute positioning. the
			%   argument indicates the [left bottom width height]
			%   of the new panel, in normalised coordinates.
			%   panels using absolute positioning mode are ignored
			%   for the sake of layout, much like items using
			%   'position:absolute' in CSS.
			%
			% see also: panel (overview), panel/panel(), select()
			
			% check m_panelType
			if p.isObject()
				error('panel:PackWhenObjectPanel', 'cannot pack() into this panel - it is already committed as an object panel');
			end
			
			% check abs mode
			%
			% we don't currently support packing children into an
			% abs panel. there's no fundamental reason why we
			% can't, it's just that abs panels aren't currently
			% included in the sizing process, so we can't apply
			% the appropriate constraints on their children. we
			% should add this in future (TODO). there are two
			% distinct ways of approaching it. probably the first
			% is preferred, to just bring the edges of an abs
			% panel into the linprog process. the alternative is
			% to do a separate linprog process within the abs
			% panel, probably much hairier.
			if isofsize(p.packpos, [1 4])
				error('panel:PackWhenAbsPanel', 'cannot pack() into this panel - it uses absolute positioning mode');
			end
			
			% if no arguments, simulate an argument of [], to pack
			% a single panel with unspecified size
			if isempty(varargin)
				varargin = {[]};
			end
			
			% state
			norender = false;
			absolute = false;
			
			% handle arguments one by one
			while ~isempty(varargin)
				
				% extract
				arg = varargin{1};
				varargin = varargin(2:end);
				
				% is char?
				if ischar(arg)
					
					% handle string arguments
					switch arg
						case {'abs' 'absolute'}
							absolute = true;
						case 'h'
							p.packdim = 1;
						case 'v'
							p.packdim = 2;
						case 'norender'
							norender = true;
						otherwise
							error('panel:InvalidArgument', ['pack() did not recognise the argument "' arg '"']);
					end
					
				else
					
					% handle numeric arguments
					if isnumeric(arg) && isscalar(arg) && arg >= 1 && arg <= 32 && arg == round(arg)
						
						% error if absolute
						if absolute
							error('panel:InvalidArgument', 'after "abs", pack() expects a [1x4] numeric argument to specify the absolute position');
						end
						
						% treat as "number of panels to pack"
						p.pack('norender', -ones(1, arg), varargin{:});
						
						% and we're done
						break
						
					elseif isfloat(arg)
						
						% if [] is the argument, convert to -1 (stretchy)
						if isempty(arg)
							arg = -1;
						end
						
						% handle absolute
						if absolute
							
							% error if absolute
							if ~isofsize(arg, [1 4])
								error('panel:InvalidArgument', 'after "abs", pack() expects a [1x4] numeric argument to specify the absolute position');
							end
							
							% error if non-positive width or height
							if any(arg(3:4) <= 0)
								error('panel:InvalidArgument', 'absolute position must have non-zero width and height');
							end
							
							% error if any panels are already packed
							for c = 1:length(p.m_children)
								if ~isofsize(p.m_children(c).packpos, [1 4])
									error('panel:IllegalPack', 'all panels pack()ed into a single parent panel must use the same positioning mode (relative or absolute)');
								end
							end
							
							% commit as parent
							p.commitAsParent();

							% create
							child = panel(p);
							
							% store passed packpos
							child.packpos = arg;
							
							% store
							if isempty(p.m_children)
								p.m_children = child;
							else
								p.m_children(end+1) = child;
							end
							
							% recurse
							if ~isempty(varargin)
								child.pack('norender', varargin{:});
							end
							
						else
							
							% error if not in range
							if size(arg, 1) ~= 1
								error('panel:InvalidArgument', 'argument to pack() must be a row vector');
							end
							
							% error if not in range
							if ~all(arg == -1 | (arg > 0 & arg <= 100)) || any(isnan(arg))
								error('panel:InvalidArgument', 'argument to pack() must contain only -1 (stretch) and non-zero values no larger than 100');
							end
							
							% error if any panels are already absolute
							for c = 1:length(p.m_children)
								if isofsize(p.m_children(c).packpos, [1 4])
									error('panel:IllegalPack', 'all panels pack()ed into a single parent panel must use the same positioning mode (relative or absolute)');
								end
							end
							
							% treat as "widths of panels to pack"
							for i = 1:length(arg)
								
								% commit as parent
								p.commitAsParent();

								% create
								child = panel(p);
								
								% store passed packpos
								child.packpos = arg(i);
								if child.packpos == -1
									child.packpos = [];
								end
								
								% store
								if isempty(p.m_children)
									p.m_children = child;
								else
									p.m_children(end+1) = child;
								end
								
								% recurse
								if ~isempty(varargin)
									subpackdim = flipdim(p.packdim);
									edges = 'hv';
									child.pack('norender', edges(subpackdim), varargin{:});
								end
								
							end
							
						end
						
						% and we're done
						break
						
					else
						
						error('panel:InvalidArgument', 'invalid numerical argument passed to pack()');
						
					end
					
				end
				
			end
			
			% this must generate a renderAll(), since the addition
			% of new panels may affect the layout. any recursive
			% call passes 'norender', so that only the root call
			% actually bothers doing a render.
			if ~norender
				p.renderAll();
			end
			
		end
		
		function h_out = select(p, h_object)
			
			% create or select an axis or object panel
			%
			% h = p.select(h)
			%   this call will return the handle of the object
			%   associated with the panel. if the panel is not yet
			%   committed, this will involve first committing it
			%   as an "object panel". if a list of objects ("h")
			%   is passed, these are the objects associated with
			%   the panel. if not, a new axis is created by the
			%   panel.
			%
			%   if the object list includes axes, then the "object
			%   panel" is also known as an "axis panel". in this
			%   case, the call to select() will make the (first)
			%   axis current, unless an output argument is
			%   requested, in which case the handle of the axes
			%   are returned but no axis is made current.
			%
			%   the passed objects can be user-created axes (e.g.
			%   a colorbar) or any graphics object that is to have
			%   its position managed (e.g. a uipanel). your
			%   mileage may vary with different types of graphics
			%   object, please let me know.
			%
			% see also: panel (overview), panel/panel(), pack()
			
			% handle "all" and "data"
			if nargin == 2 && isstring(h_object) && (strcmp(h_object, 'all') || strcmp(h_object, 'data'))
				
				% collect
				h_out = [];
				
				% commit all uncommitted panels as axis panels by
				% selecting them once
				if p.isParent()

					% recurse
					for c = 1:length(p.m_children)
						h_out = [h_out p.m_children(c).select(h_object)];
					end

				elseif p.isUncommitted()

					% select in an axis
					h_out = p.select();
					
					% plot some data
					if strcmp(h_object, 'data')
						plot(h_out, randn(100, 1), 'k-');
					end

				end
				
				% ok
				return
				
			end
			
			% check m_panelType
			if p.isParent()
				error('panel:SelectWhenParent', 'cannot select() this panel - it is already committed as a parent panel');
			end
			
			% commit as object
			p.commitAsObject();

			% assume not a new object
			newObject = false;
			
			% use passed graphics object
			if nargin >= 2
				
				% validate
				if ~all(ishandle(h_object))
					error('panel:InvalidArgument', 'argument to select() must be a list of handles to graphics objects');
				end
				
				% validate
				if ~isempty(p.h_object)
					error('panel:SelectWithObjectWhenObject', 'cannot select() new objects into this panel - it is already managing objects');
				end
				
				% store
				p.h_object = h_object;
				newObject = true;
				
				% make sure it has the correct parent - this doesn't
				% seem to affect axes, so we do it for all
 				set(p.h_object, 'parent', p.h_parent);
				
			end
			
			% create new axis if necessary
			if isempty(p.h_object)
				% 'NextPlot', 'replacechildren'
				%   make sure fonts etc. don't get changed when user
				%   plots into it
				p.h_object = axes( ...
					'Parent', p.h_parent, ...
					'NextPlot', 'replacechildren' ...
					);
				newObject = true;
			end
			
			% if wrapped objects include an axis, and no output args, make it current
			h_axes = p.getAllManagedAxes();
			if ~isempty(h_axes) && ~nargout
				set(p.h_figure, 'CurrentAxes', h_axes(1));
				
				% 12/07/11: this call is slow, because it implies "drawnow"
% 				figure(p.h_figure);

				% 12/07/11: this call is fast, because it doesn't
				set(0, 'CurrentFigure', p.h_figure);
				
			end
			
			% and return object list
			if nargout
				h_out = p.h_object;
			end
			
			% this must generate a renderPanel(), since the axis
			% will need positioning appropriately
			if newObject
				% 09/03/12 mitch
				% if there isn't a context yet, we'll have to
				% renderAll(), in fact, to generate a context first.
				% this will happen, for instance, if a single panel
				% is generated in a window that was already open
				% (no resize event will fire, and since pack() is
				% not called, it will not call renderAll() either).
				% nonetheless, we have to reposition this object, so
				% this forces us to renderAll() now and generate
				% that context we need.
				if isempty(p.m_context)
					p.renderAll();
				else
					p.renderPanel();
				end
			end
			
		end
		
	end
	
	
	
	
	
	
	
	
	
	
	
	
	
	%% ---- HIDDEN OVERLOADS ----
	
	methods (Hidden = true)
		
		function out = vertcat(p, q)
			error('panel2:MethodNotImplemented', 'concatenation is not supported by panel (use a cell array instead)');
		end
		
		function out = horzcat(p, q)
			error('panel2:MethodNotImplemented', 'concatenation is not supported by panel (use a cell array instead)');
		end
		
		function out = cat(dim, p, q)
			error('panel2:MethodNotImplemented', 'concatenation is not supported by panel (use a cell array instead)');
		end
		
		function out = ge(p, q)
			error('panel2:MethodNotImplemented', 'inequality operators are not supported by panel');
		end
		
		function out = le(p, q)
			error('panel2:MethodNotImplemented', 'inequality operators are not supported by panel');
		end
		
		function out = lt(p, q)
			error('panel2:MethodNotImplemented', 'inequality operators are not supported by panel');
		end
		
		function out = gt(p, q)
			error('panel2:MethodNotImplemented', 'inequality operators are not supported by panel');
		end
		
		function out = eq(p, q)
			out = eq@handle(p, q);
		end
		
		function out = ne(p, q)
			out = ne@handle(p, q);
		end
		
	end
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	%% ---- PUBLIC HIDDEN GET/SET ----
	
	methods (Hidden = true)
		
		function p = descend(p, indices)
			
			while ~isempty(indices)

				% validate
				if numel(p) > 1
					error('panel:InvalidIndexing', 'you can only use () on a single (scalar) panel');
				end

				% validate
				if ~p.isParent()
					error('panel:InvalidIndexing', 'you can only use () on a parent panel');
				end
				
				% extract
				index = indices{1};
				indices = indices(2:end);

				% only accept numeric
				if ~isnumeric(index) || ~isscalar(index)
					error('panel:InvalidIndexing', 'you can only use () with scalar indices');
				end
					
				% do the reference
				p = p.m_children(index);
					
			end
			
		end

		function p_out = subsasgn(p, refs, value)
			
			% output is always subject
			p_out = p;
			
			% handle () indexing
			if strcmp(refs(1).type, '()')
				p = p.descend(refs(1).subs);
				refs = refs(2:end);
			end
				
			% is that it?
			if isempty(refs)
				error('panel:InvalidIndexing', 'you cannot assign to a child panel');
			end
			
 			% next ref must be .
 			if ~strcmp(refs(1).type, '.')
				panel.error('InvalidIndexing');
			end
			
			% either one (.X) or two (.ch.X)
			switch numel(refs)
				
				case 2
			
					% validate
					if ~strcmp(refs(2).type, '.')
						panel.error('InvalidIndexing');
					end
					
					% validate
					switch refs(2).subs
						case {'fontname' 'fontsize' 'fontweight'}
						case {'margin' 'marginleft' 'marginbottom' 'marginright' 'margintop' 'align'}
						otherwise
							panel.error('InvalidIndexing');
					end
					
					% avoid re-rendering whilst setting descendant properties
					p.defer();
					
					% recurse
					switch refs(1).subs
						case {'children' 'ch'}
							cs = p.m_children;
							for c = 1:length(cs)
								subsasgn(cs(c), refs(2:end), value);
							end
						case {'descendants' 'de'}
							cs = p.getPanels('*');
							for c = 1:length(cs)
								if cs{c} ~= p
									subsasgn(cs{c}, refs(2:end), value);
								end
							end
					end
					
					% release for rendering
					p.undefer();

					% mark for appropriate updates
					refs(1).subs = refs(2).subs;
					
				case 1

					% delegate
					p.setPropertyValue(refs(1).subs, value);
					
				otherwise
					panel.error('InvalidIndexing');
	
			end
			
			% update rendering
			switch refs(1).subs
				case {'fontname' 'fontsize' 'fontweight'}
					p.renderPanel('recurse');
				case {'margin' 'marginleft' 'marginbottom' 'marginright' 'margintop' 'align'}
					p.renderAll();
			end

		end
		
		function out = subsref(p, refs)
			
			% handle () indexing
			if strcmp(refs(1).type, '()')
				p = p.descend(refs(1).subs);
				refs = refs(2:end);
			end
				
			% is that it?
			if isempty(refs)
				out = p;
				return
			end

 			% next ref must be .
 			if ~strcmp(refs(1).type, '.')
				panel.error('InvalidIndexing');
			end
			
			% switch on "fieldname"
			switch refs(1).subs
				
				case { ...
						'fontname' 'fontsize' 'fontweight' ...
						'margin' 'marginleft' 'marginbottom' 'marginright' 'margintop' ...
 						'units'  'align' ...
						}

					% delegate this property get
					out = p.getPropertyValue(refs(1).subs);
					
				case 'position'
					out = p.getObjectPosition();
					
				case 'figure'
					out = p.h_figure;
					
				case 'packpos'
					out = p.packpos;
					
				case 'axis'
					if p.isObject()
						out = p.getAllManagedAxes();
					else
						out = [];
					end
					
				case 'object'
					if p.isObject()
						h = p.h_object;
						ih = ishandle(h);
						out = h(ih);
					else
						out = [];
					end
					
 				case {'ch' 'children' 'de' 'descendants'}
					
					% validate
					if length(refs) < 2 || ~strcmp(refs(2).type, '.')
						panel.error('InvalidIndexing');
					end
					
					% validate
					switch refs(2).subs
						case {'axis' 'object'}
						otherwise
							panel.error('InvalidIndexing');
					end
					
					% handle
					out = [];
					if strcmp(refs(1).subs(1:2), 'ch')
						cs = p.m_children;
						for c = 1:length(cs)
							out = cat(2, out, subsref(cs(c), refs(2)));
						end
					else
						cs = p.getPanels('*');
						for c = 1:length(cs)
							if cs{c} ~= p
								out = cat(2, out, subsref(cs{c}, refs(2)));
							end
						end
					end
					
					% we've used two
					refs = refs(2:end);
					
				case { ...
						'addCallback' 'setCallback' 'clearCallbacks' ...
						'xlabel' 'ylabel' 'zlabel' 'title' 'hold' ...
						'refresh' 'export' ...
						'pack' 'repack' ...
						'identify' 'show' ...
						}

					% validate
					if length(refs) ~= 2 || ~strcmp(refs(2).type, '()')
						error('panel:InvalidIndexing', ['"' refs(1).subs '" is a function (try "help panel/' refs(1).subs '")']);
					end
					
					% delegate this function call with no output
					builtin('subsref', p, refs);
					return
					
				case { ...
						'select' 'fixdash' ...
						}
					
					% validate
					if length(refs) ~= 2 || ~strcmp(refs(2).type, '()')
						error('panel:InvalidIndexing', ['"' refs(1).subs '" is a function (try "help panel/' refs(1).subs '")']);
					end
					
					% delegate this function call with output
					if nargout
						out = builtin('subsref', p, refs);
					else
						builtin('subsref', p, refs);
					end
					return
					
				otherwise
					panel.error('InvalidIndexing');
							
			end
			
			% continue
			if length(refs) > 1
				out = subsref(out, refs(2:end));
			end
			
		end
		
	end
	
	
	
	
	
	
	
	
	
	
	
	
	
	%% ---- UTILITY METHODS ----
	
	methods (Access = private)
		
		function b = ismanagefont(p)
			
			% ask root
			b = p.m_root.state.manage_font;
			
		end
		
		function b = isdefer(p)
			
			% ask root
			b = p.m_root.state.defer ~= 0;
			
		end
		
		function defer(p)
			
			% increment
			p.m_root.state.defer = p.m_root.state.defer + 1;
			
		end

		function undefer(p)
			
			% decrement
			p.m_root.state.defer = p.m_root.state.defer - 1;
			
		end

		function packing_size = autoPack(p)
			
			packing_size = [];
			for c = 1:length(p.m_children)
				pos = p.m_children(c).packpos;
				if isempty(pos)
					pos = -1;
				end
				if numel(pos) == 4
					% this is an absolute positioning mode panel, so
					% we just return empty which means don't do
					% relative sizing
					packing_size = [];
					return
				end
				packing_size(c) = pos;
			end
			used = sum(packing_size(packing_size ~= -1));
			if used > 1 & used <= 100
				% interpret as percentages
				used = used / 100;
				packing_size(packing_size ~= -1) = packing_size(packing_size ~= -1) / 100;
			end
			over = 1 - used;
			packing_size(packing_size == -1) = over / length(packing_size(packing_size == -1));
			
		end
		
		function cs = getPanels(p, panelTypes, edgespec, all)
			
			% return all the panels that match the specification.
			%
			% panelTypes "*": return all panels
			% panelTypes "s": return all sizeable panels (parent,
			%		object and uncommitted)
			% panelTypes "p": return only physical panels (object
			%   or uncommitted)
			% panelTypes "o": return only object panels
			%
			% if edgespec/all is specified, only panels matching
			% the edgespec are returned (all of them if "all" is
			% true, or any of them - the first one, in fact - if
			% "all" is false).
			
			cs = {};
			
			% do not include any that use absolute positioning -
			% they stand outside of the sizing model
			skip = (numel(p.packpos) == 4) && ~any(panelTypes == '*');
			
			if p.isParent()
				
				% return if appropriate type
 				if any(panelTypes == '*s') && ~skip
					cs = {p};
 				end
				
				% if edgespec was supplied
				if nargin == 4

					% if we are perpendicular to the specified edge
					if p.packdim ~= edgespec(1)

						if all
							
							% return all matching
							for c = 1:length(p.m_children)
								ppp = p.m_children(c).getPanels(panelTypes, edgespec, all);
								cs = cat(2, cs, ppp);
							end
							
						else
							
							% return only the first one
							cs = cat(2, cs, p.m_children(1).getPanels(panelTypes, edgespec, all));
							
						end

					else

						% if we are parallel to the specified edge
						if edgespec(2) == 2
							
							% use last
							ppp = p.m_children(end).getPanels(panelTypes, edgespec, all);
							cs = cat(2, cs, ppp);
							
						else
							
							% use first
							cs = cat(2, cs, p.m_children(1).getPanels(panelTypes, edgespec, all));
							
						end

					end
					
				else
					
					% else, return all
					for c = 1:length(p.m_children)
						ppp = p.m_children(c).getPanels(panelTypes);
						cs = cat(2, cs, ppp);
					end
					
				end
				
			elseif p.isObject()
				
				% return if appropriate type
				if any(panelTypes == '*spo') && ~skip
					cs = {p};
				end
				
			else
				
				% return if appropriate type
				if any(panelTypes == '*sp') && ~skip
					cs = {p};
				end
				
			end
			
		end
		
		function commitAsParent(p)
			
			if p.isUncommitted()
				p.m_panelType = p.PANEL_TYPE_PARENT;
			elseif p.isObject()
				error('panel:AlreadyCommitted', 'cannot make this panel a parent panel, it is already an object panel');
			end

		end
		
		function commitAsObject(p)
			
			if p.isUncommitted()
				p.m_panelType = p.PANEL_TYPE_OBJECT;
			elseif p.isParent()
				error('panel:AlreadyCommitted', 'cannot make this panel an object panel, it is already a parent panel');
			end

		end
		
		function b = isRoot(p)
			
			b = isempty(p.parent);
			
		end
		
		function b = isParent(p)
			
			b = p.m_panelType == p.PANEL_TYPE_PARENT;
			
		end
		
		function b = isObject(p)
			
			b = p.m_panelType == p.PANEL_TYPE_OBJECT;
			
		end
		
		function b = isUncommitted(p)
			
			b = p.m_panelType == p.PANEL_TYPE_UNCOMMITTED;
			
		end

		function h_axes = getAllManagedAxes(p)
			
			h_axes = [];
			for n = 1:length(p.h_object)
				h = p.h_object(n);
				if isaxis(h)
					h_axes = [h_axes h];
				end
			end
			
		end
		
		function h_object = getOrCreateAxis(p)
			
			switch p.m_panelType
				
				case p.PANEL_TYPE_PARENT
					
					% create if not present
					if isempty(p.h_object)
						
						% 'Visible', 'off'
						%   this is the hidden axis of a parent panel,
						%   used for displaying a parent panel's xlabel,
						%   ylabel and title, but not as a plotting axis
						%
						% 'NextPlot', 'replacechildren'
						%   make sure fonts etc. don't get changed when user
						%   plots into it
						p.h_object = axes( ...
							'Parent', p.h_parent, ...
							'Visible', 'off', ...
							'NextPlot', 'replacechildren' ...
							);
						
						% make sure it's unitary, to help us in
						% positioning labels and title
						axis(p.h_object, [0 1 0 1]);
						
						% refresh this axis position
						p.renderPanel();
						
					end
					
					% ok
					h_object = p.h_object;
					
				case p.PANEL_TYPE_OBJECT
					
					% ok
					h_object = p.getAllManagedAxes();
					if isempty(h_object)
						error('panel:ManagedObjectNotAnAxis', 'this object panel does not manage an axis');
					end
					
				case p.PANEL_TYPE_UNCOMMITTED
					
					panel.error('PanelUncommitted');
					
			end
			
		end
		
		function removeChild(p, child)
			
			% if not a parent, fail but warn (shouldn't happen)
			if ~p.isParent()
				warning('panel:NotParentOnRemoveChild', 'i am not a parent (in removeChild())');
				return
			end
			
			% remove from children
			for c = 1:length(p.m_children)
				if p.m_children(c) == child
					p.m_children = p.m_children([1:c-1 c+1:end]);
					return
				end
			end
			
			% warn
			warning('panel:ChildAbsentOnRemoveChild', 'child not found (in removeChild())');
			
		end
		
		function h = getShowAxis(p)
			
			if p.isRoot()
				if isempty(p.h_showAxis)
					
					% create
					p.h_showAxis = axes( ...
						'Parent', p.h_parent, ...
						'units', 'normalized', ...
						'position', [0 0 1 1] ...
						);
					
					% move to bottom
					c = get(p.h_parent, 'children');
					c = [c(2:end); c(1)];
					set(p.h_parent, 'children', c);
					
					% finalise axis
					set(p.h_showAxis, ...
						'xtick', [], 'ytick', [], ...
						'color', 'none', 'box', 'off' ...
						);
					axis(p.h_showAxis, [0 1 0 1]);
					
					% hold
					hold(p.h_showAxis, 'on');
					
				end
				
				% return it
				h = p.h_showAxis;
				
			else
				h = p.parent.getShowAxis();
			end
			
		end
		
		function fireCallbacks(p, event)
		
			% for each attached callback
			for c = 1:length(p.m_callback)
				
				% extract
				callback = p.m_callback{c};
				func = callback{1};
				userdata = callback{2};
				
				% fire
				data = [];
				data.panel = p;
				data.event = event;
				data.context = p.m_context;
				data.userdata = userdata;
				func(data);
				
			end
				
		end
		
	end
	
	
	
	

		
	
	
	
	
	
		
	%% ---- RENDER METHODS ----
	
	methods

		function refresh(p)
			
			% re-render all panels
			%
			% p.refresh()
			%   re-render the panel from scratch.
			%
			% NB: if you pass 'defer' to the constructor, calling
			% refresh() both renders the panel and releases the
			% defer mode. future changes to properties (e.g.
			% margins) will be rendered immediately.
			
			% bubble up to root
			if ~p.isRoot()
				p.m_root.refresh();
				return
			end
			
			% release defer
			p.state.defer = 0;

			% debug output
			panel.debugmsg(['refresh "' p.state.name '"...']);
			
			% call renderAll
			p.renderAll();
			
		end
		
	end
		
	methods (Access = private)
		
		function do_fixdash(p, context)
			
			% if context is [], this is _after_ the render for
			% export, so we need to restore
			if isempty(context)
				
				% restore lines we changed to their original state
				for r = 1:length(p.m_fixdash_restore)
					
					% get
					restore = p.m_fixdash_restore{r};
					
					% if empty, no change was made
					if ~isempty(restore)
						set(restore.h_line, ...
							'xdata', restore.xdata, 'ydata', restore.ydata);
						delete([restore.h_supp restore.h_mark]);
					end
					
				end
				
			else
				
% 				% get handles to objects that still exist
% 				h_lines = p.m_fixdash(ishandle(p.m_fixdash));
				
				% no restores
				p.m_fixdash_restore = {};
				
				% for each line
				for i = 1:length(p.m_fixdash)
					
					% get
					fix = p.m_fixdash{i};
					
					% final check
					if ~ishandle(fix.h) || ~isequal(get(fix.h, 'type'), 'line')
						continue
					end
					
					% apply dashstyle
					p.m_fixdash_restore{end+1} = dashstyle_line(fix, context);

				end
				
			end

		end

		function p = renderAll(p, varargin)
			
			% this function renders everything from scratch. it
			% starts by calculating the sizes of everything (using
			% linprog), then starts a recursive renderPanel().
			
			% bubble up to root
			if ~p.isRoot()
				p.m_root.renderAll(varargin{:});
				return
			end
			
			% skip if disabled
			if p.isdefer()
				return
			end

			% debug output
			panel.debugmsg(['renderAll "' p.state.name '"...']);
			
			% once we reach the root, we need to create the
			% rendering context, as follows
			if nargin == 2
				
				% supplied
				context = varargin{1};
				
			else
				
				% not supplied (use figure window)
				context = [];
				context.mode = panel.RENDER_MODE_NORMAL;
				context.size_in_mm = [];
				
			end
			
			% but root itself may have a packpos
			if ~isempty(p.packpos)
				if isscalar(p.packpos)
					% this should never happen, because it should be
					% caught when the packpos is set in repack()
					warning('panel:RootPanelCannotUseRelativeMode', 'the root panel uses relative positioning mode - this is ignored');
				else
					context.rect = p.packpos;
				end
			end
			
			% get context (whole parent) size in its units
			pp = get(p.h_figure, 'position');
			context_size = pp(3:4);
			
			% defaults, in case this fails for any reason
			screen_size = [1280 1024];
			if ismac
				screen_dpi = 72;
			else
				screen_dpi = 96;
			end
				
			% get screen DPI
			try
				local_screen_dpi = get(0, 'ScreenPixelsPerInch');
				if ~isempty(local_screen_dpi)
					screen_dpi = local_screen_dpi;
				end
			end
			
			% get screen size
			try
				local_screen_size = get(0, 'ScreenSize');
				if ~isempty(local_screen_size)
					screen_size = local_screen_size;
				end
			end
			
			% if not given a context size, use the size on screen
			if isempty(context.size_in_mm)
				
				% get figure width and height on screen
				switch get(p.h_figure,'Units')
					
					case 'points'
						points_per_inch = 72;
						context.size_in_mm = context_size / points_per_inch * 25.4;
						
					case 'inches'
						context.size_in_mm = context_size * 25.4;
						
					case 'centimeters'
						context.size_in_mm = context_size * 10.0;
						
					case 'pixels'
						context.size_in_mm = context_size / screen_dpi * 25.4;
						
					case 'characters'
						context_size = context_size .* [5 13]; % convert to pixels (based on empirical measurement)
						context.size_in_mm = context_size / screen_dpi * 25.4;
						
					case 'normalized'
						context_size = context_size .* screen_size(3:4); % convert to pixels (based on screen size)
						context.size_in_mm = context_size / screen_dpi * 25.4;
						
					otherwise
						error('panel:CaseNotCoded', ['case not coded, (Parent Units are ' get(p.h_figure, 'Units') ')']);
						
				end
				
			end
			
			% that's the figure size, now we need the size of our
			% parent, if it's not the figure too
			if p.h_parent ~= p.h_figure
				pos = get(p.h_parent, 'position');
				context.size_in_mm = context.size_in_mm .* pos(3:4);
			end
			
			
			
			% LINEAR PROGRAMMING APPROACH TO LAYOUT
			%
			% parent panels represent constraints, whereas object
			% panels are the subjects of constraints. uncommitted
			% panels, for the sake of sizing, behave like object
			% panels (it's just that they won't actually display
			% anything in the rectangle sized up for them).
			%
			% we go through first and get all the constraints
			% placed by all the parent panels of the layout. we
			% then maximize the size of all panels subject to
			% those constraints.
			%
			% but before we do that, we have to index all the
			% non-parent panels so that we can refer to them by
			% index in the parameters of the linear programming
			% operation.
			%
			% CONSTRAINT TYPES:
			%
			% constraint 1: margin constraints
			%   each margin element must be respected. an object
			%   must have clear space around itself of the amount
			%   specified by its margin. these are inequality
			%   constraints.
			%
			% constraint 2: relative size constraints
			%   when we pack() panels, we ask them to be of
			%   certain sizes relative to one another. these are
			%   equality constraints.
			%
			% constraint 3: alignment constraints
			%   if two panels are packed alongside, then the axes
			%   within them on their outside edges, perpendicular
			%   to the packing dimension, must align. these are
			%   equality constraints.
			
			
% 			t0 = clock();
% 			t_linprog = 0;
% 			r_linprog = [];
			
			% index all sizeable panels (keep references to them
			% as well so that we can give them their rect when
			% we've finished computing the layout). numPanels is the
			% number of these panels. numVars is twice this,
			% because, for each dimension, there are two variables
			% to find.
			allPanels = p.getPanels('*');
			
			% separate sizeables from those with absolute
			% positioning mode
			abs_mode = logical([]);
			for i = 1:length(allPanels)
				if isofsize(allPanels{i}.packpos, [1 4])
					abs_mode(i) = true;
				else
					abs_mode(i) = false;
				end
			end
			absPanels = allPanels(abs_mode);
			allPanels = allPanels(~abs_mode);
			
			% index the sizeable ones
			numPanels = length(allPanels);
			numVars = numPanels * 2;
			for i = 1:numPanels
				allPanels{i}.state.index = i;
			end
			
			% create LP problem
			lp = linprog_create(numVars);
			p.state.lp = repmat(lp, 1, 2);
			p.state.lp(1).size_in_mm = context.size_in_mm(1);
			p.state.lp(2).size_in_mm = context.size_in_mm(2);
			
			% add constraints
			p.addConstraints(p, [1 1]);
			p.addConstraintsFigureEdge();
			
			% add maximization
			for dim = 1:2
				for n = 1:numPanels
					i2 = n * 2;
					i1 = i2 - 1;
					p.state.lp(dim) = linprog_addmaxim(p.state.lp(dim), i1, i2);
				end
			end
			
			% data for panels
			xxyy = [];
			
			% solve lp problems
			p.state.lp(1) = linprog_solve(p.state.lp(1));
			p.state.lp(2) = linprog_solve(p.state.lp(2));
			opt_success = ~isempty(p.state.lp(1).x) && ~isempty(p.state.lp(2).x);
			
			% if successful
			if opt_success
			
				% distribute solution amongst panels
				for dim = 1:2
					xxyy(:, dim*2+[-1 0]) = reshape(p.state.lp(dim).x, 2, numPanels)' / context.size_in_mm(dim);
				end

				% because we pack top-bottom (zero to one), but the
				% actual object positions in matlab are +ve is upwards,
				% we have to flip the y data, here.
				yy = xxyy(:, [3 4]);
				yy = fliplr(1 - yy);
				xxyy(:, [3 4]) = yy;

				% pass context to each panel
				for n = 1:numPanels
					context.rect = [xxyy(n, [1 3]) xxyy(n, [2 4])-xxyy(n, [1 3])];
					allPanels{n}.m_context = context;
				end

				% pass context to abs panels
				for n = 1:length(absPanels)
					context.rect = [];
					absPanels{n}.m_context = context;
				end

				% clear h_showAxis
				if ~isempty(p.h_showAxis)
					delete(p.h_showAxis);
					p.h_showAxis = [];
				end

	% 			% timing report
	% 			tf = clock();
	% 			t_total = etime(tf, t0);
	%  			disp(['LINPROG: ' sprintf('%dms (%dms)', round([t_total*1000 t_linprog*1000])), ' constraints ' r_linprog]);
	
			end
			
			% now we can begin rendering into that context
			p.renderPanel('recurse');
			
		end
		
		function renderPanel(p, varargin)
			
			% skip if disabled
			if p.isdefer()
				return
			end
			
			% debug output
			panel.debugmsg(['renderPanel "' p.state.name '"...']);
			
			% defaults
			recurse = false;
			
			% handle arguments
			while ~isempty(varargin)
				
				% get
				arg = varargin{1};
				varargin = varargin(2:end);
				
				% handle
				switch arg
					case 'recurse'
						recurse = true;
					otherwise
						panel.error('InternalError');
				end
				
			end
			
			% recurse
			if recurse
				pp = p.getPanels('*');
			else
				pp = {p};
			end
			
			% who do we have to split the render cycle into two?
			% because the "group labels" are positioned with
			% respect to the axes in their group depending on
			% whether those axes have tick labels, and what those
			% tick labels are. if those tick labels are in
			% automatic mode (as they usually are), we don't know
			% this until those axes have been positioned. since an
			% axis group may contain many of these nested deep, we
			% have to position all axes (cycle 1) first, then
			% (cycle 2) position any group labels.
			
			% render cycle 1
			for pi = 1:length(pp)
				pp{pi}.renderPanel1();
			end
			
			% render cycle 2
			for pi = 1:length(pp)
				pp{pi}.renderPanel2();
			end
			
			% callbacks
			for pi = 1:length(pp)
				fireCallbacks(pp{pi}, 'render-complete');
			end
			
		end
		
		function r = getObjectPosition(p)
			
			% get packed position
			r = p.m_context.rect;
			
			% if empty, must be absolute position
			if isempty(r)
				r = p.packpos;
				pp = getObjectPosition(p.parent);
				r = panel.getRectangleOfRectangle(pp, r);
			end
			
		end
		
		function renderPanel1(p)
			
			% if no context yet, skip this call
			if isempty(p.m_context)
				return
			end
			
			% if no managed objects, skip this call
			if isempty(p.h_object)
				return
			end

			% debug output
			panel.debugmsg(['renderPanel1 "' p.state.name '"...']);
			
			% handle RENDER_MODE
			switch p.m_context.mode

				case panel.RENDER_MODE_PREPRINT

					% if in RENDER_MODE_PREPRINT, store current axis
					% layout (ticks and ticklabels) and lock them into
					% manual mode so they don't get changed during the
					% print operation
					h_axes = p.getAllManagedAxes();
					for n = 1:length(h_axes)
						p.state.store{n} = storeAxisState(h_axes(n));
					end

				case panel.RENDER_MODE_POSTPRINT

					% if in RENDER_MODE_POSTPRINT, restore axis
					% layout, leaving it as it was before we ran
					% export
					h_axes = p.getAllManagedAxes();
					for n = 1:length(h_axes)
						restoreAxisState(h_axes(n), p.state.store{n});
					end

			end
			
			% position it
			try
				set(p.h_object, 'position', p.getObjectPosition(), 'units', 'normalized');
			catch err
				if strcmp(err.identifier, 'MATLAB:hg:set_chck:DimensionsOutsideRange')
					w = warning('query', 'backtrace');
					warning off backtrace
					warning('panel:PanelZeroSize', 'a panel had zero size, and the managed object was hidden');
					set(p.h_object, 'position', [-0.3 -0.3 0.2 0.2]);
					if strcmp(w.state, 'on')
						warning on backtrace
					end
				elseif strcmp(err.identifier, 'MATLAB:class:InvalidHandle')
					% this will happen if the user deletes the managed
					% objects manually. an obvious way that this
					% happens is if the user select()s some panels so
					% that axes get created, then calls clf. it would
					% be nice if we could clear the panels attached to
					% a figure in response to a clf call, but there
					% doesn't seem any obvious way to pick up the clf
					% call, only the delete(objects) that follows, and
					% this is indistinguishable from a call by the
					% user to delete(my_axis), for instance. how are
					% we to respond if the user deletes the axis the
					% panel is managing? it's not clear. so, we'll
					% just fail silently, for now, and these panels
					% will either never be used again (and will be
					% destroyed when the figure is closed) or will be
					% destroyed when the user creates a new panel on
					% this figure. either way, i think, no real harm
					% done.
% 					w = warning('query', 'backtrace');
% 					warning off backtrace
% 					warning('panel:PanelObjectDestroyed', 'the object managed by a panel has been destroyed');
% 					if strcmp(w.state, 'on')
% 						warning on backtrace
% 					end
					panel.debugmsg('***WARNING*** the object managed by a panel has been destroyed');
					return
				else
					rethrow(err)
				end
			end

			% if managing fonts
			if p.ismanagefont()
				
				% apply properties to objects
				h = p.h_object;
				
				% get those which are axes
				h_axes = p.getAllManagedAxes();

				% and labels/title objects, for any that are axes
				for n = 1:length(h_axes)
					h = [h ...
						get(h_axes(n), 'xlabel') ...
						get(h_axes(n), 'ylabel') ...
						get(h_axes(n), 'zlabel') ...
						get(h_axes(n), 'title') ...
						];
				end

				% apply font properties
				set(h, ...
					'fontname', p.getPropertyValue('fontname'), ...
					'fontsize', p.getPropertyValue('fontsize'), ...
					'fontweight', p.getPropertyValue('fontweight') ...
					);
				
			end

		end
			
		function renderPanel2(p)
			
			% if no context yet, skip this call
			if isempty(p.m_context)
				return
			end
			
			% if no object, skip this call
			if isempty(p.h_object)
				return
			end

			% if not a parent, skip this call
			if ~p.isParent()
				return
			end

			% if not an axis, skip this call - NB: this is not a
			% displayed and managed object, rather it is the
			% invisible axis used to display parent labels/titles.
			% we checked above if this panel is a parent. thus,
			% the member h_object must be scalar, if it is
			% non-empty.
			if ~isaxis(p.h_object)
				return
			end

			% debug output
			panel.debugmsg(['renderPanel2 "' p.state.name '"...']);
			
			% matlab moves x/ylabels around depending on
			% whether the axis in question has any x/yticks,
			% so that the label is always "near" the axis.
			% we try to do the same, but it's hack-o-rama.

			% calibration offsets - i measured these
			% empirically, what a load of shit
			font_fudge = [2 1/3];
			nofont_fudge = [2 0];

			% do xlabel
			cs = p.getPanels('o', [2 2], true);
			y = 0;
			for c = 1:length(cs)
				ch = cs{c};
				h_axes = ch.getAllManagedAxes();
				for h_axis = h_axes
					% only if there are some tick labels, and they're
					% at the bottom...
					if ~isempty(get(h_axis, 'xticklabel')) && ~isempty(get(h_axis, 'xtick')) ...
							&& strcmp(get(h_axis, 'xaxislocation'), 'bottom')
						fontoffset_mm = get(h_axis, 'fontsize') * font_fudge(2) + font_fudge(1);
						y = max(y, fontoffset_mm);
					end
				end
			end
			y = max(y, get(p.h_object, 'fontsize') * nofont_fudge(2) + nofont_fudge(1));

			% convert and lay in
			axisheight_mm = p.m_context.size_in_mm(2) * p.m_context.rect(4);
			y = y / axisheight_mm;
			set(get(p.h_object, 'xlabel'), ...
				'VerticalAlignment', 'Cap', ...
				'Units', 'Normalized', ...
				'Position', [0.5 -y 1]);

			% calibration offsets - i measured these
			% empirically, what a load of shit
			font_fudge = [3 1/6];
			nofont_fudge = [2 0];

			% do ylabel
			cs = p.getPanels('o', [1 1], true);
			x = 0;
			for c = 1:length(cs)
				ch = cs{c};
				h_axes = ch.getAllManagedAxes();
				for h_axis = h_axes
					% only if there are some tick labels, and they're
					% at the left...
					if ~isempty(get(h_axis, 'yticklabel')) && ~isempty(get(h_axis, 'ytick')) ...
							&& strcmp(get(h_axis, 'yaxislocation'), 'left')
						yt = get(h_axis, 'yticklabel');
						if ischar(yt)
							ml = size(yt, 2);
						else
							ml = 0;
							for i = 1:length(yt)
								ml = max(ml, length(yt{i}));
							end
						end
						fontoffset_mm = get(h_axis, 'fontsize') * ml * font_fudge(2) + font_fudge(1);
						x = max(x, fontoffset_mm);
					end
				end
			end
			x = max(x, get(p.h_object, 'fontsize') * nofont_fudge(2) + nofont_fudge(1));

			% convert and lay in
			axisheight_mm = p.m_context.size_in_mm(1) * p.m_context.rect(3);
			x = x / axisheight_mm;
			set(get(p.h_object, 'ylabel'), ...
				'VerticalAlignment', 'Bottom', ...
				'Units', 'Normalized', ...
				'Position', [-x 0.5 1]);

			% calibration offsets - made up based on the
			% ones i measured for the labels
			nofont_fudge = [2 0];

			% get y position
			y = max(y, get(p.h_object, 'fontsize') * nofont_fudge(2) + nofont_fudge(1));

			% convert and lay in
			axisheight_mm = p.m_context.size_in_mm(2) * p.m_context.rect(4);
			y = y / axisheight_mm;
			set(get(p.h_object, 'title'), ...
				'VerticalAlignment', 'Bottom', ...
				'Position', [0.5 1+y 1]);

		end
		
	end
	
	
	
	
	
	
	
	
	
	
	%% ---- CONSTRAINT METHODS ----
	
	methods (Access = private)

		function addConstraint_Alignment(p, placer, p1, p2, edgespec)
			
			% "placer" is the calling panel, being some parent panel,
			% and is placing a constraint between "p1" and "p2",
			% that their edges, specified by edgespec, must be
			% aligned.
% 			nc = rpad(int2str(size(p.state.lin(1).A, 1)+size(p.state.lin(2).A, 1)+1), 4);
% 			disp([nc 'ALIGN   ' lpad(placer.state.name) ' --- ' ...
% 				lpad(p1.state.name) ' : ' rpad(p2.state.name) ...
% 				'align on edge ' edgestr(edgespec) ...
% 				]);
			
			% add
			dim = edgespec(1);
			p1 = p1.state.index;
			p2 = p2.state.index;
			i1 = (p1-1) * 2 + edgespec(2);
			i2 = (p2-1) * 2 + edgespec(2);
			p.state.lp(dim) = linprog_equality(p.state.lp(dim), i1, i2);
			
		end
		
		function addConstraint_RelativeSize(p, placer, p1, p2, dim, ratio)
			
			% "placer" is the calling panel, being some parent panel,
			% and is placing a constraint between "f1/f2" and "o1/o2",
			% that they have the specified size "ratio" along the
			% specified "dim"ension.
% 			nc = rpad(int2str(size(p.state.lin(1).A, 1)+size(p.state.lin(2).A, 1)+1), 4);
% 			disp([nc 'RELSIZE ' lpad(placer.state.name) ' --- ' ...
% 				lpad(p1.state.name) ' : ' rpad(p2.state.name) ...
% 				sprintf('%i  %.3f', dim, ratio) ...
% 				]);
			
			% add
			p1 = p1.state.index;
			p2 = p2.state.index;
			i1 = p1*2-1;
			i2 = p1*2;
			i3 = p2*2-1;
			i4 = p2*2;
			p.state.lp(dim) = linprog_addrelsize(p.state.lp(dim), i1, i2, i3, i4, ratio);
			
		end
		
		function addConstraint_EdgeMargin(p, placer, dim, p1, e1, p2, e2, margin)
			
			% "placer" is the calling panel, being some parent panel,
			% and is placing a constraint between the edges p1/e1
			% and p2/e2, with at least the specified margin
			% between the edges. e1 and e1 are edgespec's. either
			% p1 or p2 can be empty, in which case the figure is
			% assumed to be that side of the constraint.
% 			if isempty(p1) p1name = 'figure'; else p1name = p1.state.name; end
% 			if isempty(p2) p2name = 'figure'; else p2name = p2.state.name; end
% 			nc = rpad(int2str(size(p.state.lin(1).A, 1)+size(p.state.lin(2).A, 1)+1), 4);
% 			disp([nc 'EDGE    ' lpad(placer.state.name) ' --- ' ...
% 				lpad([p1name ' ' edgestr([dim e1])]) ' : ' rpad([p2name ' ' edgestr([dim e2])]) ...
% 				sprintf('[margin = %i]', margin) ...
% 				]);
			
			% add
			if isempty(p1)
				i1 = {0};
			else
 				p1 = p1.state.index;
				i1 = (p1-1) * 2 + e1;
			end
			if isempty(p2)
				i2 = {p.state.lp(dim).size_in_mm};
			else
 				p2 = p2.state.index;
				i2 = (p2-1) * 2 + e2;
			end
			p.state.lp(dim) = linprog_addmargin(p.state.lp(dim), i1, i2, margin);

			
		end
		
		function addConstraintsFigureEdge(p)
			
			% ---- FIGURE EDGE CONSTRAINTS ----

			% here, we only use the top-level margin for this. an
			% alternative is to use the margin of each individual
			% panel (see below), but this seems to be unhelpful in
			% making layouts. doing it this way is part of the
			% general policy of "margins apply only at the level
			% that they are set" - thus, the margin setting only
			% of the top-level container (the root panel) is
			% applied wrt the figure edge.
			margin = p.getPropertyValue('margin', 'mm');
			
			% for each dim
			for dim = 1:2
				
				% for each edge
				for edge = 1:2
					
					% get all children on edge
					cs = p.getPanels('s', [dim edge], true);

					% place margins with figure edges
					for c = 1:length(cs)
						ch = cs{c};
						
						% use the margin of each individual panel
% 						margin = ch.getPropertyValue('margin', 'mm');

						% extract the margin for the appropriate edge
						m = margin(edgeindex([dim edge]));
						
						if edge == 1
							p.addConstraint_EdgeMargin(p, dim, [], 2, ch, 1, m);
						else
							p.addConstraint_EdgeMargin(p, dim, ch, 2, [], 1,  m);
						end
					end
					
				end

			end
			
		end
		
		function addConstraints(p, root, xy)
			
			% parameters
			packdim = p.packdim;
			
			% as part of the policy of "margins apply only at the
			% level that they are set", we take the margins only
			% of the parent containers when applying margins.
			
			% ---- MARGIN CONSTRAINTS ----
			
			% only multi-child parents place internal margin constraints
			if p.isParent() && length(p.m_children) > 1
				
				% for each child pair
				for c = 1:(length(p.m_children)-1)
					
					% extract local pair of children, left then right
					% (or top then bottom)
					p1 = p.m_children(c);
					p2 = p.m_children(c+1);
					
					% and get the appropriate margins from them
					m1 = p1.getPropertyValue('margin', 'mm');
					m2 = p2.getPropertyValue('margin', 'mm');
					m1 = m1(edgeindex([packdim 2]));
					m2 = m2(edgeindex([packdim 1]));
					
					% get all children on the appropriate edge of each
					% of them
					p1c = p1.getPanels('s', [packdim 2], true);
					p2c = p2.getPanels('s', [packdim 1], true);
					
					% for each pair, get the maximum margin between
					% them (they each specify a margin) and apply that
					for i1 = 1:length(p1c)
						for i2 = 1:length(p2c)
							
							% get individual children
							c1 = p1c{i1};
							c2 = p2c{i2};
							
% 							% this alternative model uses the margins of
% 							% the individual children, not just of the
% 							% packing parents
% 	 						m1 = c1.getPropertyValue('margin', 'mm');
% 	 						m2 = c2.getPropertyValue('margin', 'mm');
% 							m1 = m1(edgeindex([packdim 2]));
% 							m2 = m2(edgeindex([packdim 1]));
							
							% both margins give a constraint - the tighter
							% constraint is given simply by the larger margin
							m = max(m1, m2);
							
							% apply the margin
		 					root.addConstraint_EdgeMargin(p, packdim, c1, 2, c2, 1, m);
							
						end
					end
					
				end
				
			end
			
			
			
			% ---- ALIGNMENT CONSTRAINTS ----
			
			if p.isParent() && p.getPropertyValue('align')
				
				% align on opposite to packdim
				aligndim = flipdim(packdim);
				
				% collect items on extreme edges
 				e1 = {};
				e2 = {};
				
				% for each child, get items on edge
				for c = 1:length(p.m_children)

					% get extreme edges
					ch = p.m_children(c);
					e1 = cat(2, e1, ch.getPanels('p', [aligndim 1], true));
					e2 = cat(2, e2, ch.getPanels('p', [aligndim 2], true));

				end
				
				% for each, align with the first
				for e = 2:length(e1)
					root.addConstraint_Alignment(p, e1{1}, e1{e}, [aligndim 1]);
				end
				
				% for each, align with the first
				for e = 2:length(e2)
					root.addConstraint_Alignment(p, e2{1}, e2{e}, [aligndim 2]);
				end

			end
			
			
			
			% ---- CONTAINMENT CONSTRAINTS ----
			
			if p.isParent()
				
				% for each dim
				for containdim = 1:2
					
					% in packing dim, contain first and last
					if containdim == packdim
						
						% get extreme edges of first and last
						o1 = p.m_children(1).getPanels('s', [containdim 1], false);
						o2 = p.m_children(end).getPanels('s', [containdim 2], false);

						% add a constraint for the edge
						if ~isempty(o1)
							root.addConstraint_EdgeMargin(p, containdim, p, 1, o1{1}, 1, 0);
						end
						
						% add a constraint for each edge
						if ~isempty(o2)
							root.addConstraint_EdgeMargin(p, containdim, o2{1}, 2, p, 2, 0);
						end
						
					else
						
						% across packing-dim, contain all
						for c = 1:length(p.m_children)

							% get extreme edges of child
							ch = p.m_children(c);
							o1 = ch.getPanels('s', [containdim 1], false);
							o2 = ch.getPanels('s', [containdim 2], false);

							% add a constraint for the edge
							if ~isempty(o1)
								root.addConstraint_EdgeMargin(p, containdim, p, 1, o1{1}, 1, 0);
							end
							
							% add a constraint for the edge
							if ~isempty(o2)
								root.addConstraint_EdgeMargin(p, containdim, o2{1}, 2, p, 2, 0);
							end

						end
						
					end
					
				end
				
			end
			
			
			
			% ---- SIZING CONSTRAINTS ----
			
			if p.isParent()

				% auto-pack to get the size of stretch panels
				packing_size = p.autoPack();
				
				% if that returned empty, this panel contains
				% absolute positioning mode panels, and we can skip
				% it
				if ~isempty(packing_size)

					% between the first and each other, place a relative
					% size relationship. we size against the object edges
					% that are exposed by each child. if the child is an
					% object, this will just be the two edges of its object;
					% if the child is a parent, this will be the edges
					% of its two extreme object panel children,
					% grandchildren, etc.

					% get first child
					ch = p.m_children(1);

					% for each other child
					for c = 2:length(p.m_children)

						% get other child
						och = p.m_children(c);

						% add a sizing constraint for the child
						d = packing_size(c) / packing_size(1);
						root.addConstraint_RelativeSize(p, ch, och, packdim, d);

					end

				end
				
			end
			
			
			
			% ---- RECURSE ----
			
			if p.isParent()
				
				% recurse
				for c = 1:length(p.m_children)
					
					if isempty(packing_size)
						
						% absolute packing
						xy_ = panel.getRectangleOfRectangle([0 0 xy], p.m_children(c).packpos);
						p.m_children(c).addConstraints(root, xy_(3:4));
						
					else
						
						% relative packing
						xy_ = xy;
						xy_(packdim) = xy_(2) * packing_size(c);
						p.m_children(c).addConstraints(root, xy_);
					
					end
					
				end
				
			end
			
		end

	end
	
	
	
	
	
	
	%% ---- PROPERTY METHODS ----
	
	methods (Access = private)
		
		function value = getPropertyValue(p, key, units)

			value = p.prop.(key);
			
			if isempty(value)

				% inherit
				if ~isempty(p.parent)
					switch key
						case {'fontname' 'fontsize' 'fontweight' 'margin' 'units'}
							if nargin == 3
								value = p.parent.getPropertyValue(key, units);
							else
								value = p.parent.getPropertyValue(key);
							end
							return
					end
				end
				
				% default
				if isempty(value)
					value = panel.getPropertyDefault(key);
				end
				
			end
			
			% translate dimensions
			switch key
				case {'margin'}
					if nargin < 3
						units = p.getPropertyValue('units');
					end
					value = panel.resolveUnits(value, units);
			end
			
		end
		
		function setPropertyValue(p, key, value)
			
			% root properties
			switch key
				case 'units'
					if ~isempty(p.parent)
						p.parent.setPropertyValue(key, value);
						return
					end
			end
			
			% value validation
			switch key
				case 'units'
					invalid = ~( (isstring(value) && isin({'mm' 'in' 'cm' 'pt'}, value)) || isempty(value) );
				case 'fontname'
					invalid = ~( isstring(value) || isempty(value) );
				case 'fontsize'
					invalid = ~( (isnumeric(value) && isscalar(value) && value >= 4 && value <= 60) || isempty(value) );
				case 'fontweight'
					invalid = ~( (isstring(value) && isin({'normal' 'bold'}, value)) || isempty(value) );
				case 'margin'
					invalid = ~( (isdimension(value)) || isempty(value) );
				case {'marginleft' 'marginbottom' 'marginright' 'margintop'}
					invalid = ~isscalardimension(value);
				case 'align'
					invalid = ~( (isscalar(value) && (isnumeric(value) || islogical(value))) || isempty(value) );
					value = logical(value);
				otherwise
					error('panel:UnrecognisedProperty', ['unrecognised property "' key '"']);
			end
			
			% value validation
			if invalid
				error('panel:InvalidValueForProperty', ['invalid value for property "' key '"']);
			end
			
			% marginX properties
			switch key
				case {'marginleft' 'marginbottom' 'marginright' 'margintop'}
					index = isin({'left' 'bottom' 'right' 'top'}, key(7:end));
					element = value;
					value = p.getPropertyValue('margin');
					value(index) = element;
					key = 'margin';
			end
			
			% translate dimensions
			switch key
				case {'margin'}
					if isscalar(value)
						value = value * [1 1 1 1];
					end
					if ~isempty(value)
						units = p.getPropertyValue('units');
						value = {panel.resolveUnits({value units}, 'mm') 'mm'};
					end
			end
			
			% lay in
			p.prop.(key) = value;
			
		end
		
	end	
		
	methods (Static = true, Access = private)
	
		function prop = getPropertyInitialState()
			
			prop = panel.getPropertyDefaults();
			for key = fieldnames(prop)'
				prop.(key{1}) = [];
			end
			
		end
		
		function value = getPropertyDefault(key)
			
			persistent defprop
			
			if isempty(defprop)
				defprop = panel.getPropertyDefaults();
			end
			
			value = defprop.(key);
			
		end
		
		function defprop = getPropertyDefaults()
			
			% root properties
			defprop.units = 'mm';
			
			% inherited properties
			defprop.fontname = get(0, 'defaultAxesFontName');
			defprop.fontsize = get(0, 'defaultAxesFontSize');
			defprop.fontweight = 'normal';
			defprop.margin = {[15 15 5 5] 'mm'};
			
			% not inherited properties
			defprop.align = false;
			
		end
		
	end	
		
	
	

	
	
	
	%% ---- STATIC PUBLIC METHODS ----
	
	methods (Static = true)
	
		function p = recover(h_figure)
			
			% get a handle to the root panel associated with a figure
			%
			% p = recover(h_fig)
			%   if you have not got a handle to the root panel of
			%   the figure h_fig, this call will retrieve it. if
			%   h_fig is not supplied, gcf is used.
			
			if nargin < 1
				h_figure = gcf;
			end
			
			p = panel.callbackDispatcher('recover', h_figure);
			
		end
		
		function panic()
			
			% call delete on all children of the global workspace,
			% to recover from bugs that leave us with uncloseable
			% figures. call this as "panel.panic()".
			%
			% NB: if you have to call panic(), something has gone
			% wrong. if you are able to reproduce the problem,
			% please contact me to report the bug.
			delete(allchild(0));
			
		end
		
	end
	
	
	
	
	
	
	%% ---- STATIC PRIVATE METHODS ----
	
	methods (Static = true, Access = private)
		
		function error(id)

			switch id
				case 'PanelUncommitted'
					throwAsCaller(MException('panel:PanelUncommitted', 'this action cannot be performed on an uncommitted panel'));
				case 'InvalidIndexing'
					throwAsCaller(MException('panel:InvalidIndexing', 'you cannot index a panel object in this way'));
				case 'InternalError'
					throwAsCaller(MException('panel:InternalError', 'an internal error occurred'));
				otherwise
					throwAsCaller(MException('panel:UnknownError', ['an unknown error was generated with id "' id '"']));
			end
				
		end
		
		function lockClass()
			
			persistent hasLocked
			if isempty(hasLocked)
				
				% only lock if not in debug mode
				if ~panel.isDebug()
					% in production code, must mlock() file at this point,
					% to avoid persistent variables being cleared by user
					if strcmp(getenv('USERDOMAIN'), 'BERGEN')
						% do nothing
					else
						mlock
					end
				end
				
				% mark that we've handled this
				hasLocked = true;
				
			end
			
		end
		
		function debugmsg(msg)
			
			% display, if in debug mode
			if panel.isDebug()
				disp(msg);
			end
			
		end
		
		function state = isDebug()
			
			% persistent
			persistent debug
			
			% create
			if isempty(debug)
				try
					debug = panel_debug_state();
				catch
	 				debug = false;
				end
			end
			
			% ok
			state = debug;
			
		end
		
		function r = getFractionOfRectangle(r, dim, range)
			
			switch dim
				case 1
					r = [r(1)+range(1)*r(3) r(2) range(2)*r(3) r(4)];
				case 2
					r = [r(1) r(2)+(1-sum(range))*r(4) r(3) range(2)*r(4)];
				otherwise
					error('panel:CaseNotCoded', ['case not coded, dim = ' dim ' (internal error)']);
			end
			
		end
		
		function r = getRectangleOfRectangle(r, s)
			
			w = r(3);
			h = r(4);
			r = [r(1)+s(1)*w r(2)+s(2)*h s(3)*w s(4)*h];
			
		end
		
		function a = getUnionRect(a, b)
			
			if isempty(a)
				a = b;
			end
			if ~isempty(b)
				d = a(1) - b(1);
				if d > 0
					a(1) = a(1) - d;
					a(3) = a(3) + d;
				end
				d = a(2) - b(2);
				if d > 0
					a(2) = a(2) - d;
					a(4) = a(4) + d;
				end
				d = b(1) + b(3) - (a(1) + a(3));
				if d > 0
					a(3) = a(3) + d;
				end
				d = b(2) + b(4) - (a(2) + a(4));
				if d > 0
					a(4) = a(4) + d;
				end
			end
			
		end
		
		function r = reduceRectangle(r, margin)
			
			r(1:2) = r(1:2) + margin(1:2);
			r(3:4) = r(3:4) - margin(1:2) - margin(3:4);
			
		end
		
		function v = normaliseDimension(v, space_size_in_mm)
			
			v = v ./ [space_size_in_mm space_size_in_mm];
			
		end
		
		function v = resolveUnits(d, units)
			
			% first, convert into mm
			v = d{1};
			switch d{2}
				case 'mm'
					% ok
				case 'cm'
					v = v * 10.0;
				case 'in'
					v = v * 25.4;
				case 'pt'
					v = v / 72.0 * 25.4;
				otherwise
					error('panel:CaseNotCoded', ['case not coded, storage units = ' units ' (internal error)']);
			end
			
			% then, convert to specified units
			switch units
				case 'mm'
					% ok
				case 'cm'
					v = v / 10.0;
				case 'in'
					v = v / 25.4;
				case 'pt'
					v = v / 25.4 * 72.0;
				otherwise
					error('panel:CaseNotCoded', ['case not coded, requested units = ' units ' (internal error)']);
			end
			
		end
		
		function resizeCallback(obj, evt)
			
			panel.callbackDispatcher('resize', obj);
			
		end
		
		function closeCallback(obj, evt)
			
			panel.callbackDispatcher('delete', obj);
			delete(obj);
			
		end
		
		function out = callbackDispatcher(op, data)
			
			% debug output
			panel.debugmsg(['callbackDispatcher(' op ')...'])
			
			% persistent store
			persistent registeredPanels
			
			% switch on operation
			switch op
				
				case {'register' 'registerNoClear'}
					
					% if a root panel is already attached to this
					% figure, we could throw an error and refuse to
					% create the new object, we could delete the
					% existing panel, or we could allow multiple
					% panels to be attached to the same figure.
					%
					% we should allow multiple panels, because they
					% may have different parents within the same
					% figure (e.g. uipanels). but by default we don't,
					% unless the panel.add() static constructor is
					% used.
					
					if strcmp(op, 'register')
						
						argument_h_figure = data.h_figure;
						i = 0;
						while i < length(registeredPanels)
							i = i + 1;
							if registeredPanels(i).h_figure == argument_h_figure
								delete(registeredPanels(i));
								i = 0;
							end
						end
						
					end
					
					% register the new panel
					if isempty(registeredPanels)
						registeredPanels = data;
					else
						registeredPanels(end+1) = data;
					end
					
					% debug output
					panel.debugmsg(['panel registered (' int2str(length(registeredPanels)) ' now registered)']);
					
				case 'unregister'
					
					% debug output
					panel.debugmsg(['on unregister, ' int2str(length(registeredPanels)) ' registered']);
					
					for r = 1:length(registeredPanels)
						if registeredPanels(r) == data
							registeredPanels = registeredPanels([1:r-1 r+1:end]);

							% debug output
							panel.debugmsg(['panel unregistered (' int2str(length(registeredPanels)) ' now registered)']);
							
							return
						end
					end
					
					% warn
					warning('panel:AbsentOnCallbacksUnregister', 'panel was absent from the callbacks register when it tried to unregister itself');
					
				case 'resize'
					
					argument_h_parent = data;
					for r = 1:length(registeredPanels)
						if registeredPanels(r).h_parent == argument_h_parent
							registeredPanels(r).renderAll();
						end
					end
					
				case 'recover'
					
					argument_h_figure = data;
					out = [];
					for r = 1:length(registeredPanels)
						if registeredPanels(r).h_figure == argument_h_figure
							if isempty(out)
								out = registeredPanels(r);
							else
								out(end+1) = registeredPanels(r);
							end
						end
					end
					
				case 'delete'
					
					argument_h_figure = data;
					i = 0;
					while i < length(registeredPanels)
						i = i + 1;
						if registeredPanels(i).h_figure == argument_h_figure
							delete(registeredPanels(i));
							i = 0;
						end
					end
					
			end
			
		end
		
	end
	
	
	
	
end

















% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% HELPERS
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function restore = dashstyle_line(fix, context)

% get axis size in mm
h_line = fix.h;
h_axis = get(h_line, 'parent');
u = get(h_axis, 'units');
set(h_axis, 'units', 'norm');
pos = get(h_axis, 'position');
set(h_axis, 'units', u);
axis_in_mm = pos(3:4) .* context.size_in_mm;

% recover data
xdata = get(h_line, 'xdata');
ydata = get(h_line, 'ydata');
zdata = get(h_line, 'zdata');
linestyle = get(h_line, 'linestyle');
marker = get(h_line, 'marker');

% empty restore
restore = [];

% do not handle 3D
if ~isempty(zdata)
	warning('panel:NoFixdash3D', 'panel cannot fixdash() a 3D line - no action taken');
	return
end

% get range of axis
ax = axis(h_axis);

% get scale in each dimension (mm per unit)
sc = axis_in_mm ./ (ax([2 4]) - ax([1 3]));

% create empty line
data = NaN;

% override linestyle
if ~isempty(fix.linestyle)
	linestyle = fix.linestyle;
end

% transcribe linestyle
linestyle = dashstyle_parse_linestyle(linestyle);
if isempty(linestyle)
	return
end

% scale
scale = 1;
dashes = linestyle * scale;

% store for restore
restore.h_line = h_line;
restore.xdata = xdata;
restore.ydata = ydata;

% create another, separate, line to overlay on the original
% line and render the fixed-up dashes.
restore.h_supp = copyobj(h_line, h_axis);

% if the original line has markers, we'll have to create yet
% another separate line instance to represent them, because
% they shouldn't be "dashed", as it were. note that we don't
% currently attempt to get the z-order right for these
% new lines.
if ~isequal(marker, 'none')
	restore.h_mark = copyobj(h_line, h_axis);
	set(restore.h_mark, 'linestyle', 'none');
	set(restore.h_supp, 'marker', 'none');
else
	restore.h_mark = [];
end

% hide the original line. this line remains in existence so
% that if there is a legend, it doesn't get messed up.
set(h_line, 'xdata', NaN, 'ydata', NaN);

% extract pattern length
patlen = sum(dashes);

% position within pattern is initially zero
pos = 0;

% linedata
line_xy = complex(xdata, ydata);

% for each line segment
while length(line_xy) > 1
	
	% get line segment
	xy = line_xy(1:2);
	line_xy = line_xy(2:end);
	
	% any NaNs, and we're outta here
	if any(isnan(xy))
		continue
	end
	
	% get start etc.
	O = xy(1);
	V = xy(2) - xy(1);
	
	% get mm length of this line segment
	d = sqrt(sum(([real(V) imag(V)] .* sc) .^ 2));
	
	% and mm unit vector
	U = V / d;
	
	% generate a long-enough pattern for this segment
	n = ceil((pos + d) / patlen);
	pat = [0 cumsum(repmat(dashes, [1 n]))] - pos;
	pos = d - (pat(end) - patlen);
	pat = [pat(1:2:end-1); pat(2:2:end)];
	
	% trim spurious segments
	pat = pat(:, any(pat >= 0) & any(pat <= d));
	
	% skip if that's it
	if isempty(pat)
		continue
	end
	
	% and reduce ones that are oversized
	pat(1) = max(pat(1), 0);
	pat(end) = min(pat(end), d);

	% finally, add these segments to the line data
	seg = [O + pat * U; NaN(1, size(pat, 2))];
	data = [data seg(:).'];
	
end

% update line
set(restore.h_supp, 'xdata', real(data), 'ydata', imag(data), ...
	'linestyle', '-');

end


function linestyle = dashstyle_parse_linestyle(linestyle)

if isequal(linestyle, 'none') || isequal(linestyle, '-')
	linestyle = [];
	return
end

while 1

	% if numbers
	if isnumeric(linestyle)
		if ~isa(linestyle, 'double') || ~isrow(linestyle) || mod(length(linestyle), 2) ~= 0
			break
		end
		% no need to parse
		return
	end

	% else, must be char
	if ~ischar(linestyle) || ~isrow(linestyle)
		break
	end
	
	% translate matlab non-standard codes into codes we can
	% easily parse
	switch linestyle
		case ':'
			linestyle = '.';
		case '--'
			linestyle = '-';
	end
	
	% must be only - and .
	if any(linestyle ~= '.' & linestyle ~= '-')
		break
	end
	
	% transcribe
	c = linestyle;
	linestyle = [];
	for l = c
		switch l
			case '-'
				linestyle = [linestyle 2 0.75];
			case '.'
				linestyle = [linestyle 0.5 0.75];
		end
	end
	return

end

warning('panel:BadFixdashLinestyle', 'unusable linestyle in fixdash()');
linestyle = [];

end



% MISCELLANEOUS

function index = isin(list, value)

for i = 1:length(list)
	if strcmp(value, list{i})
		index = i;
		return
	end
end

index = 0;

end

function dim = flipdim(dim)

% this function, used between arguments in a recursive call,
% causes the dim to be switched with each recurse, so that
% we build a grid, rather than a long, long row.
dim = 3 - dim;

end



% STRING PADDING FUNCTIONS

function s = rpad(s, l)

if nargin < 2
	l = 16;
end

if length(s) < l
	s = [s repmat(' ', 1, l - length(s))];
end

end

function s = lpad(s, l)

if nargin < 2
	l = 16;
end

if length(s) < l
	s = [repmat(' ', 1, l - length(s)) s];
end

end



% HANDLE GRAPHICS HELPERS

function h = getParentFigure(h)

if strcmp(get(h, 'type'), 'figure')
	return
else
	h = getParentFigure(get(h, 'parent'));
end

end

function addHandleCallback(h, name, func)

% % get current list of callbacks
% callbacks = get(h, name);
%
% % if empty, turn into a cell
% if isempty(callbacks)
% 	callbacks = {};
% elseif iscell(callbacks)
% 	% only add ourselves once
% 	for c = 1:length(callbacks)
% 		if callbacks{c} == func
% 			return
% 		end
% 	end
% else
% 	callbacks = {callbacks};
% end
%
% % and add ours (this is friendly, in case someone else has a
% % callback attached)
% callbacks{end+1} = func;
%
% % lay in
% set(h, name, callbacks);

% the above isn't as simple as i thought - for now, we'll
% just stamp on any existing callbacks
set(h, name, func);

end

function store = storeAxisState(h)

% LOCK TICKS AND LIMITS
%
% (LOCK TICKS)
%
% lock state so that the ticks and labels do not change when
% the figure is resized for printing. this is what the user
% will expect, which is why we go through this palaver.
%
% however, for fuck's sake. the following code illustrates
% an idiosyncrasy of matlab (i would call this an
% inconsistency, myself, but there you go).
%
% figure
% axis([0 1 0 1])
% set(gca, 'ytick', [-1 0 1 2])
% get(gca, 'yticklabel')
% set(gca, 'yticklabelmode', 'manual')
%
% now, resize the figure window. at least in R2011b, the
% tick labels change on the first resize event. presumably,
% this is because matlab treats the ticklabel value
% differently depending on if the ticklabelmode is auto or
% manual. if it's manual, the value is used as documented,
% and [0 1] is used to label [-1 0 1 2], cyclically.
% however, if the ticklabelmode is auto, and the ticks
% extend outside the figure, then the ticklabels are set
% sensibly, but the _value_ of ticklabel is not consistent
% with what it would need to be to get this tick labelling
% were the mode manual. and, in a final bizarre twist, this
% doesn't become evident until the resize event. i think
% this is a bug, no other way of looking at it; at best it's
% an inconsistency that is either tedious or impossible to
% work around in the general case.
%
% in any case, we have to lock the ticks to manual as we go
% through the print cycle, so that the ticks do not get
% changed if they were in automatic mode. but we mustn't fix
% the tick labels to manual, since if we do we may encounter
% this inconsistency and end up with the wrong tick labels
% in the print out. i can't, at time of writing, think of a
% case where we'd have to fix the tick labels to manual too.
% the possible cases are:
%
% ticks auto, labels auto: in this case, fixing the ticks to
%		manual should be enough.
%
% ticks manual, labels auto: leave as is.
%
% ticks manual, labels manual: leave as is.
%
% the only other case is ticks auto, labels manual, which is
% a risky case to use, but in any case we can also fix the
% ticks to manual in that case. thus, our preferred solution
% is to always switch the ticks to manual, if they're not
% already, and otherwise leave things be.
%
% (LOCK LIMITS)
%
% the other thing that may get modified, if the user hasn't
% fixed it, is the axis limits. so we lock them too, any
% that are set to auto, and mark them for unlocking when the
% print is complete.

store = '';

% manual-ise ticks on any axis where they are currently
% automatic, and indicate that we need to switch them back
% afterwards.
if strcmp(get(h, 'XTickMode'), 'auto')
	store = [store 'X'];
	set(h, 'XTickMode', 'manual');
end
if strcmp(get(h, 'YTickMode'), 'auto')
	store = [store 'Y'];
	set(h, 'YTickMode', 'manual');
end
if strcmp(get(h, 'ZTickMode'), 'auto')
	store = [store 'Z'];
	set(h, 'ZTickMode', 'manual');
end

% manual-ise limits on any axis where they are currently
% automatic, and indicate that we need to switch them back
% afterwards.
if strcmp(get(h, 'XLimMode'), 'auto')
	store = [store 'x'];
	set(h, 'XLimMode', 'manual');
end
if strcmp(get(h, 'YLimMode'), 'auto')
	store = [store 'y'];
	set(h, 'YLimMode', 'manual');
end
if strcmp(get(h, 'ZLimMode'), 'auto')
	store = [store 'z'];
	set(h, 'ZLimMode', 'manual');
end

% % OLD CODE OBSOLETED 25/01/12 - see notes above
% 
% % store current state
% store.XTick = get(h, 'XTick');
% store.XTickMode = get(h, 'XTickMode');
% store.XTickLabel = get(h, 'XTickLabel');
% store.XTickLabelMode = get(h, 'XTickLabelMode');
% store.YTickMode = get(h, 'YTickMode');
% store.YTick = get(h, 'YTick');
% store.YTickLabel = get(h, 'YTickLabel');
% store.YTickLabelMode = get(h, 'YTickLabelMode');
% store.ZTick = get(h, 'ZTick');
% store.ZTickMode = get(h, 'ZTickMode');
% store.ZTickLabel = get(h, 'ZTickLabel');
% store.ZTickLabelMode = get(h, 'ZTickLabelMode');
% 
% % lock state to manual
% set(h, 'XTickLabelMode', 'manual');
% set(h, 'XTickMode', 'manual');
% set(h, 'YTickLabelMode', 'manual');
% set(h, 'YTickMode', 'manual');
% set(h, 'ZTickLabelMode', 'manual');
% set(h, 'ZTickMode', 'manual');

end

function restoreAxisState(h, store)

% unmanualise
for item = store
	switch item
		case {'X' 'Y' 'Z'}
			set(h, [item 'TickMode'], 'auto');
		case {'x' 'y' 'z'}
			set(h, [upper(item) 'TickMode'], 'auto');
	end
end

% % OLD CODE OBSOLETED 25/01/12 - see notes above
% 
% % restore passed state
% set(h, 'XTick', store.XTick);
% set(h, 'XTickMode', store.XTickMode);
% set(h, 'XTickLabel', store.XTickLabel);
% set(h, 'XTickLabelMode', store.XTickLabelMode);
% set(h, 'YTick', store.YTick);
% set(h, 'YTickMode', store.YTickMode);
% set(h, 'YTickLabel', store.YTickLabel);
% set(h, 'YTickLabelMode', store.YTickLabelMode);
% set(h, 'ZTick', store.ZTick);
% set(h, 'ZTickMode', store.ZTickMode);
% set(h, 'ZTickLabel', store.ZTickLabel);
% set(h, 'ZTickLabelMode', store.ZTickLabelMode);

end



% DIM AND EDGE HANDLING

% we describe each edge of a panel in terms of "dim" (1 or
% 2, horizontal or vertical) and "edge" (1 or 2, former or
% latter). together, [dim edge] is an "edgespec".

function s = edgestr(edgespec)

s = 'lbrt';
s = s(edgeindex(edgespec));

end

function i = edgeindex(edgespec)

% edge indices. margins are stored as [l b r t] but
% dims are packed left to right and top to bottom, so
% relationship between 'dim' and 'end' and index into
% margin is non-trivial. we call the index into the margin
% the "edgeindex". an "edgespec" is just [dim end], in a
% single array.
i = [1 3; 4 2];
i = i(edgespec(1), edgespec(2));

end



% VARIABLE TYPE HELPERS

function val = validate_par(val, argtext, varargin)

% this helper validates arguments to some functions in the
% main body

for n = 1:length(varargin)
	
	% get validation constraint
	arg = varargin{n};
	
	% handle string list
	if iscell(arg)
		% string list
		if ~isin(arg, val)
			error('panel:InvalidArgument', ...
				['invalid argument "' argtext '", "' val '" is not a recognised data value for this option']);
		end
		continue;
	end
	
	% handle strings
	if isstring(arg)
		switch arg
			case 'empty'
				if ~isempty(val)
					error('panel:InvalidArgument', ...
						['invalid argument "' argtext '", option does not expect any data']);
				end
			case 'dimension'
				if ~isdimension(val)
					error('panel:InvalidArgument', ...
						['invalid argument "' argtext '", option expects a dimension']);
				end
			case 'scalar'
				if ~(isnumeric(val) && isscalar(val) && ~isnan(val))
					error('panel:InvalidArgument', ...
						['invalid argument "' argtext '", option expects a scalar value']);
				end
			case 'nonneg'
				if any(val(:) < 0)
					error('panel:InvalidArgument', ...
						['invalid argument "' argtext '", option expects non-negative values only']);
				end
			case 'integer'
				if any(val(:) ~= round(val(:)))
					error('panel:InvalidArgument', ...
						['invalid argument "' argtext '", option expects integer values only']);
				end
		end
		continue;
	end
	
	% handle numeric range
	if isnumeric(arg) && isofsize(arg, [1 2])
		if any(val(:) < arg(1)) || any(val(:) > arg(2))
			error('panel:InvalidArgument', ...
				['invalid argument "' argtext '", option data must be between ' num2str(arg(1)) ' and ' num2str(arg(2))]);
		end
		continue;
	end
	
	% not recognised
	arg
	error('panel:InternalError', 'internal error - bad argument to validate_par (above)');
	
end

end

function b = checkpar(value, mn, mx)

b = isscalar(value) && isnumeric(value) && ~isnan(value);
if b
	if nargin >= 2
		b = b && value >= mn;
	end
	if nargin >= 3
		b = b && value <= mx;
	end
end

end

function b = isintegral(v)

b = all(all(v == round(v)));

end

function b = isstring(value)

sz = size(value);
b = ischar(value) && length(sz) == 2 && sz(1) == 1 && sz(2) >= 1;

end

function b = isdimension(value)

b = isa(value, 'double') && (isscalar(value) || isofsize(value, [1 4]));

end

function b = isscalardimension(value)

b = isa(value, 'double') && isscalar(value);

end

function b = isofsize(value, siz)

sz = size(value);
b = length(sz) == length(siz) && all(sz == siz);

end

function b = isaxis(h)

b = ishandle(h) && strcmp(get(h, 'type'), 'axes');

end






%% LINPROG PROBLEM
%
% NB: This could really replace the functions above
%
% addConstraint_Alignment
% addConstraint_RelativeSize
% addConstraint_EdgeMargin
%
% in a later version, so long as everything is working
% correctly then we can get rid of all lines with the old
% implementation tag "lin" in them.
		
function lp = linprog_create(N)

lp = [];
lp.N = N;
lp.A = zeros(0, N);
lp.b = [];
lp.c = [];
lp.x = [];
lp.xo = [];

% collect constraints
lp.margins = {};
lp.relsizes = {};

end

function lp = linprog_addrelsize(lp, i1, i2, i3, i4, relsize)

% set (x(i4) - x(i3)) - (relsize * (x(i2) - x(i1))) = 0

% store
lp.relsizes{end+1} = [i1 i2 i3 i4 relsize];

end

function lp = linprog_processrelsizes(lp)

for m = 1:length(lp.relsizes)

	relsize = lp.relsizes{m};
	i1 = relsize(1);
	i2 = relsize(2);
	i3 = relsize(3);
	i4 = relsize(4);
	relsize = relsize(5);
	
	row = size(lp.A, 1) + 1;
	lp.A(row, i4) = 1;
	lp.A(row, i3) = -1;
	lp.A(row, i2) = -relsize;
	lp.A(row, i1) = relsize;
	lp.b(row, 1) = 0;

end

end

function lp = linprog_addmaxim(lp, i1, i2)

% maximise x(i2) - x(i1)

if (length(lp.c) >= i1 && lp.c(i1)) || (length(lp.c) >= i2 && lp.c(i2))
	error('assumption that maxims are only introduced once has broken');
end

lp.c(i2, 1) = 1;
lp.c(i1, 1) = -1;

end

function lp = linprog_addmargin(lp, i1, i2, margin)

% set x(i2) - x(i1) >= margin

% for convenience of the caller, we allow these to build up
% one by one into this cell array.
%
% for efficiency of constructing the A matrix, we actually
% don't modify the A and b and c matrix until at the end, in
% linprog_processmargins, avoiding lots of spurious copies
% of a large lp object.
lp.margins{end+1} = {i1 i2 margin};

end

function lp = linprog_processmargins(lp)

% number of margins we will add
N = length(lp.margins);

% count up indices
row = size(lp.A, 1);
i3 = size(lp.A, 2);

% pre-allocate expanded A array
sz = size(lp.A) + N;
lp.A(sz(1), sz(2)) = 0;

% for each margin
for m = 1:length(lp.margins)
	
	% extract
	margin = lp.margins{m};
	i1 = margin{1};
	i2 = margin{2};
	margin = margin{3};

	% introduce a slack variable
	i3 = i3 + 1;
	row = row + 1;

	% then x(i2) - x(i1) - x(i3) = margin
	bval = margin;

	% if i1 or i2 is a cell, it's a constant, so treat it thus
	if iscell(i2)
		bval = bval - i2{1};
		i2 = [];
	end
	if iscell(i1)
		bval = bval + i1{1};
		i1 = [];
	end

	% add the constraint
	lp.b(row, 1) = bval;
	lp.A(row, i3) = -1;
	lp.A(row, i2) = 1;
	lp.A(row, i1) = -1;

	% augment c
	lp.c(i3, 1) = 0;

end

end

function lp = linprog_equality(lp, i1, i2)

% set x(i2) - x(i1) = 0

row = size(lp.A, 1) + 1;
lp.A(row, i2) = 1;
lp.A(row, i1) = -1;
lp.b(row, 1) = 0;

end

function lp = linprog_solve(lp)

lp = linprog_processrelsizes(lp);
lp = linprog_processmargins(lp);

persistent have_opt_toolbox have_warned
if isempty(have_opt_toolbox)
	have_opt_toolbox = ~isempty(which('linprog'));
end

if have_opt_toolbox

	% find using linprog() (Optimization Toolbox)
	NN = size(lp.A, 2);
	opt = optimset('Display', 'off');
	x = linprog(-lp.c, -eye(NN), zeros(NN, 1), lp.A, lp.b, [], [], [], opt);
	
else
	
% 	c1 = clock();
	
	% make sure all constraints have non-negative b
	f = find(lp.b < 0);
	lp.A(f, :) = -lp.A(f, :);
	lp.b(f) = -lp.b(f);
	
	% find using Jeff Stuart's implementation
	x = js_linprog_stub(lp.A, lp.b, lp.c');
	
% 	c2 = clock();
% 	
% 	T = etime(c2, c1);
% 	if T > 1 && isempty(have_warned)
% 		disp(['Panel warning: rendering is slow because optimization toolbox is not installed.']);
% 		have_warned = true;
% 	end

end

if isempty(x)
	lp.x = [];
else
	lp.x = x(1:lp.N);
	if any(isnan(lp.x))
		% something wrong - let's give up rather than bother the
		% user
% 		error('some NaNs in lp.x');
		lp.x = [];
	end
end

end







%% LINPROG SOLVER
%
% this LP solver is due to Jeff Stuart (was at USM). i
% incorporated it here in version 2.1 so that users no
% longer need the Optimization Toolbox to run this. thanks
% to "Daniel" at Matlab Central for pointing out this
% problem.



function x = js_linprog_stub(A, b, c)

% this is a stub written by ben mitch to interface the LP
% problem that panel has to solve with Jeff Stuart's code
% for solving an LP problem in standard form. i've hacked
% Jeff's code around to remove all the messaging and what
% have you, and to simply raise an error if something goes
% awry, on the grounds that the LP problem solved by panel
% either works, or does not have a feasible solution (which
% means we can't do the resizing according to the user's
% constraints at the current window size).


try

	[zmax,PHIiter,PHIIiter,xbasic,ibasic] = js_linprog(A,b,c);

	x = NaN(size(A, 2), 1);
	x(ibasic) = xbasic;
	
catch err
	
 	switch err.identifier
		
		case 'panel:InfeasibleLP'
			% no solution - figure window too small, probably
			x = [];
			return
			
		otherwise
			rethrow(err);
			
	end
	
end

end




function [zmax,PHIiter,PHIIiter,xbasic,ibasic] = js_linprog(A,b,c);
%
%LINPROG uses the two phase simplex method to solve the linear
%program maximize cx subject to the constraints Ax = b and x >= 0 ,
%where A is m x n , rank(A) = m , and b >= 0 .    The output vector
%is [zmax,PHIiter,PHIIiter,xbasic,ibasic], where zmax is the maximal
%objective value, PHIiter and PHIIiter are the phase I and phase II
%iteration counts, respectively, where xbasic is the vector of basic
%x variables at optimality, and where ibasic is the indices of the
%optimal basis columns in A (and hence the indices for the entries
%in xbasic).  LINPROG detects infeasibility and unboundedness, and
%provides appropriate output messages in such cases.  LINPROG also
%contains a heuristic check for cycling, terminating the algorithm
%when m Phase II iterations occur without a change in the objective
%value.  See also PHASEI and PHASEII.
%
%Written for MATLAB version 5.0
%
%Written by Jeff Stuart, Department of Mathematics,
%University of Southern Mississippi, Hattiesburg, MS 39406.
%December, 1993.  Revised, October, 1997.
%jeffrey.stuart@usm.edu
%
[m,n]=size(A);
if max(size(b) ~=[m 1]);
	disp('The dimensions of b do not match the dimensions of A.')
	return
end
if min(b) < 0;
	disp('The RHS vector b must be nonnegative.')
	return
end
if max(size(c) ~=[1 n]);
	disp('The dimensions of c do not match the dimensions of A.')
	return
end
if rank(A) ~=m;
	disp('A does not have full row rank.')
	return
end
PHIiter=0;
PHIIiter=0;
tol=eps;
xbasic=zeros(1,n);
[wmax,ibasic,PHIiter]=js_linprog_phase1(A,b);
if wmax < -tol
	error('panel:InfeasibleLP', 'LP is infeasible');
	%      disp('The original LP is infeasible.  Infeasibility was')
	%      disp('detected during Phase I.  The total number of phase')
	%      disp('one iterations performed was: '), disp(PHIiter)
else
% 	disp('Phase I completed.  Original LP is feasible.')
% 	disp('The total number of Phase I iterations was: '),disp(PHIiter)
% 	disp('Starting Phase II.')
	[zmax,xbasic,ibasic,ienter,PHIIiter,PCOL,OPTEST,CYCTEST]=js_linprog_phase2(A,b,c,ibasic);
	xbasic=xbasic';
% 	if CYCTEST==1;
% 		return
% 	end
% 	if OPTEST == 0;
% 		disp('The orginal LP is unbounded. An unbounded ray was')
% 		disp('detected during Phase II.  The output objective')
% 		disp('value is for the last basic solution found.')
% 		disp('The number of Phase II iterations was: '),disp(PHIIiter)
% 		disp('Last objective value is '),disp(zmax)
% 		disp('The last basic solution, xbasic is '),disp(xbasic)
% 		disp('The column indices for the last basis: '),disp(ibasic)
% 		disp('The index of the unbounded entering variable: '),disp(ienter)
% 		disp('The unbounded ray column is: '),disp(PCOL)
% 	else
% 		disp('The original LP has an optimal solution.')
% 		disp('The number of Phase II iterations was: '),disp(PHIIiter)
% 		disp('The optimal objective value is '),disp(zmax)
% 		disp('The indices for the basic columns: '),disp(ibasic)
% 		disp('The optimal, basic solution is '),disp(xbasic)
% 	end
end

end









function [wmax,ibasic,PHIiter] = js_linprog_phase1(A,b)
%PHASEI performs Phase I of the simplex method on the constraints
%Ax = b and x >= 0 (where A is m x n, rank(A) = m , and b >= 0)
%to determine whether there is a feasible point.  The function
%output is wmax, the artificial objective value; ibasic, the
%indices of the basic variables at optimality; and PHIiter, the
%number of Phase I iterations performed.  If  wmax < 0 ,then the
%original LP is infeasible.  If wmax = 0 , the original LP is
%feasible, and ibasic is the index set for a feasible basis.
%To allow for round-off error, the tests are  wmax < -tol  for
%infeasibility, and  wmax >= -tol  for feasibility, where "tol"
% is a preset tolerance (see the initialization value in the fourth
%line below).  At the expense of additional computation, an adaptive
%choice for "tol" based on A and b could be selected.
%
%See also PHASEII and LINPROG.
%Written for MATLAB version 5.0 .
%Written by Jeff Stuart, Department of Mathematics, University of
%Southern Mississippi, Hattiesburg, MS 39406.  October, 1997.
%jeff.stuart@usm.edu

[m,n]=size(A);
A=[A,eye(m)];
PHIiter=0;
tol=0.0000001;
ztol=0.0000001;
X=zeros(1,n+m);
J=[zeros(1,n),ones(1,m)];
c=-J;
K=[1:n+m];
J=logical(J);
ibasic=K(J);
inon=K(~J);
B=eye(m);
xbasic=b;
w=-sum(xbasic);
X(ibasic)=b;
Cred= ones(1,m)*A(:,inon);
loop =1;
while loop ==1;
	if max(Cred) > ztol ;
		PHIiter=PHIiter + 1;
		[Maxcost,j]=max(Cred);
		ienter=inon(j);
		PCOL=B\A(:,ienter);
		J(ienter)=1;
		TESTROWS=find(PCOL > ztol);
		TESTCOL=PCOL(TESTROWS);
		[minrat,j]=min(xbasic(TESTROWS)./TESTCOL);
		iexit=ibasic(TESTROWS(j));
		J(iexit)=0;
		if minrat > 0;
			xbasic=xbasic - minrat*PCOL;
		end
		X(ibasic)=xbasic;
		X(ienter)=minrat;
		X(iexit)=0;
		w=w + Maxcost*minrat;
		ibasic=K(J);
		inon=K(~J);
		B=A(:,ibasic);
		xbasic=X(ibasic)';
		Cred=c(inon) - (c(ibasic)/B)*A(:,inon);
	elseif Cred <= ztol;
		loop = 0;
	end
end
wmax=-sum(X(n+1:n+m));
if wmax >= -tol;
	X=X(1:n);
	last=ibasic(m);
	K=K(1:last);
	ibasic=K(J(1:last));
	while last > n;
		J=J(1:last);
		K=[1:last];
		inon=K(~J);
		B=A(:,ibasic);
		inon=inon(inon <= n);
		j=find(([zeros(1,m-1),1]/B)*A(:,inon));
		ienter=inon(j(1));
		J(ienter)=1;
		J(last)=0;
		ibasic=K(J);
		last=ibasic(m);
		PHIiter=PHIiter+1;
	end
end

end




function [z,xbasic,ibasic,ienter,iter,PCOL,OPTEST,CYCTEST] = js_linprog_phase2(A,b,c,ibasic);
%PHASEII performs phase II of the simplex method starting with the basic
%columns specified by the vector ibasic.
%
%See also PHASEI and LINPROG.
%Written for Matlab version 5.0.
%
%Written by Jeff Stuart, Department of Mathematics, University of Southern Mississippi,
%Hattiesburg, MS 39406. October, 1993. Revised October, 1997.
%jeffrey.stuart@usm.edu

[m,n]=size(A);
PCOL=[];
ienter=[];
iter=0;
cycle=0;
CYCTEST=0;
X=zeros(1,n);
J=X;
J(ibasic)=ones(1,m);
K=[1:n];
inon=K(~J);
B=A(:,ibasic);
xbasic=B\b;
z=c(ibasic)*xbasic;
OPTEST=0;
if m<n;
	X(ibasic)=xbasic;
	Cred=c(inon) - (c(ibasic)/B)*A(:,inon);
	OPTEST=1;
	loop =1;
	while loop ==1;
		if max(Cred) > 0;
			iter=iter + 1;
			[Maxcost,j]=max(Cred);
			ienter=inon(j);
			PCOL=B\A(:,ienter);
			if PCOL <= 0 , OPTEST = 0;
				loop = 0;
			else
				J(ienter)=1;
				TESTROWS=find(PCOL > 0);
				TESTCOL=PCOL(TESTROWS);
				[minrat,j]=min(xbasic(TESTROWS)./TESTCOL);
				if minrat <=0, cycle = cycle+1;
					if cycle > m;
						error('panel:ExcessiveCycling', 'excessive cycling');
% 						disp('Algorithm terminated due to excessive cycling.')
% 						disp('Restart algorithm from phase II using a perturbed')
% 						disp(' RHS vector b and the current basis.')
% 						disp(ibasic)
% 						CYCTEST=1;
% 						break
					end
				else
					cycle = 0;
				end
				iexit=ibasic(TESTROWS(j));
				J(iexit)=0;
				xbasic=xbasic - minrat*PCOL;
				X(ibasic)=xbasic;
				X(ienter)=minrat;
				X(iexit)=0;
				z=z + Maxcost*minrat;
				J=logical(J);
				ibasic=K(J);
				inon=K(~J);
				B=A(:,ibasic);
				xbasic=X(ibasic)';
				Cred=c(inon) - (c(ibasic)/B)*A(:,inon);
			end
		else
			loop = 0;
		end
	end
end

end






