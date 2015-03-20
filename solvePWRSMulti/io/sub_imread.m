%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   This file is part of the HCI-Correspondece Estimation Benchmark.
%
%   More information on this benchmark can be found unter:
%       http://hci.iwr.uni-heidelberg.de/Benchmarks/
%
%    Copyright (C) 2011  <name of author>
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [X, map, alpha, warnCode, warnMsg] = sub_imread(varargin);
% This function replaces the imread() function coming with Matlab and contains
% workarounds for Matlab bugs and shortcomings.  You can use it in the same
% syntax as imread().
%
% Currently implemented improvements over imread():
% (1) Correct support of 12 bit PGM images
%
% IN: see documentation of imread()
% OUT:   X : image, see documentation of imread()
%      map : colormap, see documentation of imread()
%    alpha : see documentation of imread()
% warnCode : 0 if okay, unique warning code number otherwise
%  warnMsg : string, warning message
%
%------------------------------------------------------------------------------
% Details:
% ========
%
% (1) Correct support of 12 bit PGM images
%     ------------------------------------
% Matlab's imread function has a bug that occurs when reading 2-byte-PGM files
% with a max value that does not fully exploit the two bytes.  For example
% when the header of the PGM file reads like this (bit-depth = 12):
%
% P5
% # 
% 1024 512
% 4095
%
% Matlab will not understand that the maximum value that can occur is 2^12-1 =
% 4095.  Instead, Matlab will scale all values with 65535/4095 and round them,
% which is really nonsense!
% The workaround function provided here reverts this scaling.  Before doing
% this, it checks whether the file type is PGM with 2 bytes per pixel.
%
%------------------------------------------------------------------------------

%==============================================================================

%------------------------------------------------------------------------------
% Settings

Verbose = 1; % be talkative when a workaround or bugfix is applied

%------------------------------------------------------------------------------
if nargin < 1,
    error('ERROR (sub_imread): too few arguments');
end

%------------------------------------------------------------------------------
fileName = varargin{1};

% initialize warning message
warnCode = 0;
warnMsg = '';

if ~exist(fileName, 'file'),
    error(['ERROR (sub_imread): file not found: ', fileName]);
end

[X, map, alpha] = imread(varargin{:});

% get file information
fi = imfinfo(fileName);

if strcmp(fi.Format, 'PGM') & (fi.BitDepth > 8) & (fi.BitDepth < 16),
    warnCode = 1;
    warnMsg = sprintf(['INFO (sub_imread): workaround (1) for PGM with ', ...
                      'bit-depth %g is applied'], fi.BitDepth);

    % the input and output data format is uint16
    X = uint16(double(X)*((2^fi.BitDepth-1)/(2^16-1)));
end

if (nargout < 5) & (warnCode == 1),
    fprintf('%s\n', warnMsg);
end
