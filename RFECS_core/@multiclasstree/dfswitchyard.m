function varargout = dfswitchyard(action,varargin)
% DFSWITCHYARD switchyard for Distribution Fitting.
% Helper function for the Distribution Fitting tool

%   Copyright 2003-2007 The MathWorks, Inc.
%   $Revision: 1.1.6.3 $ $Date: 2007/05/23 19:15:24 $

% Calls from Java prefer the if/else version.
% [varargout{1:max(nargout,1)}]=feval(action,varargin{:});
if nargout==0
	feval(action,varargin{:});
else    
	[varargout{1:nargout}]=feval(action,varargin{:});
end

% The following lines list functions that are called via this
% function from other Statistics Toolbox functions.  These lines
% insure that the compiler will include the functions listed.
%#function statgetargs
%#function mgrp2idx
%#function dfgetdistributions
%#function dfhistbins
