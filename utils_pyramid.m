% UTILS_PYRAMID solves pyramid of any image

%   Author: Menandro Roxas
%   Contact: menandro.roxas@gmail.com
%   $Date: 2020-11-16 (Monday, 16 Nov 2020) $

% Copyright 2020, Menandro Roxas.
%
%                         All Rights Reserved
%
% Permission to use, copy, modify, and distribute this software and its
% documentation for any purpose other than its incorporation into a
% commercial product is hereby granted without fee, provided that the
% above copyright notice appear in all copies and that both that
% copyright notice and this permission notice appear in supporting
% documentation, and that the name of the author not be used in
% advertising or publicity pertaining to distribution of the software
% without specific, written prior permission.
%
% THE AUTHOR DISCLAIM ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
% INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR ANY
% PARTICULAR PURPOSE.  IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
% ANY SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
% WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
% ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
% OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE. 

function varargout = utils_pyramid(varargin) %src,npyr,pyscale
    switch nargin
        case 3
            src = varargin{1};
            npyr = varargin{2};
            pyscale = varargin{3};
            minWidth = 16;
            mode = 'bicubic';
        case 4
            src = varargin{1};
            npyr = varargin{2};
            pyscale = varargin{3};
            minWidth = varargin{4};
            mode = 'bicubic';
        case 5
            src = varargin{1};
            npyr = varargin{2};
            pyscale = varargin{3};
            minWidth = varargin{4};
            mode = varargin{5};
    end
    [h,w,~]=size(src);
    cnt = 0;
    sizes={};
    while (w>=minWidth && h>=minWidth && cnt<npyr)
        cnt = cnt+1;
        sizes{cnt} = [h,w];
        w = int32(floor(pyscale*w));
        h = int32(floor(pyscale*h));
    end
    py = {};
    py{1} = src;
    mask = (1/256)*[1 4 6 4 1]'*[1 4 6 4 1];
    for i=2:cnt
        %py{i-1} = imfilter(py{i-1},mask,'corr');
        py{i} = imresize(py{i-1},sizes{i},mode);
    end
    py{1} = src;
    nLevel = cnt;
    varargout{1}=py;
    varargout{2}=nLevel;
    varargout{3}=sizes;
        
%---------------------------------------
function [src,npyr,pyscale,minWidth] = parse_inputs(varargin)
    narginchk(3,4);
    validateattributes(varargin{1},{'uint8','double'},{'real','nonsparse'},...
        mfilename,'I,X or RGB',1);
    switch nargin
        case 3
            src = varargin{1};
            npyr = varargin{2};
            pyscale = varargin{3};
        case 4
            src = varargin{1};
            npyr = varargin{2};
            pyscale = varargin{3};
            minWidth = varargin{4};
    end