function info = adiosopen(varargin)
%ADIOSOPEN Open an ADIOS BP file and provide information on it.
%   
%   ADIOSOPEN opens an ADIOS BP file and returns a structure that
%   contains the file handler and information on all groups,
%   variables and attributes.
%
%   If you have opened the file with ADIOSOPEN, you can supply the
%   file handler and the group handler to read in a variable/attribute
%   with ADIOSREAD. When finished reading all data, you need to close 
%   the file with ADIOSCLOSE. 
%
%   Note: the adios group is not the same as HDF5 groups. Each variable has
%   a path, which defines a logical hierarchy of the variables within one 
%   adios group. This logical hierarchy is what is similar to HDF5 groups.
%
%   INFO = ADIOSOPEN(FILE) 
%      Open FILE and return an information structure. 
%
%   The returned INFO structure is the following
%
%     Name        file path
%     FileHandler int64 file handler
%     TimeStart   First timestep in file
%     TimeEnd     Last timestep in file
%     Groups      Adios groups in the file. Usually 1 group is in a file.
%                 This is a structure array of 
%
%        Name          group name
%        GroupHandler  int64 group handler
%        Variables     structure array of variables
%           Name          path of variable
%           Type          Matlab type class of data
%           Dims          Array of dimensions
%           Timedim       The time dimension, 0 if there is no time varying 
%                         part of the variable
%           GlobalMin     global minimum  of the variable (1-by-1 mxArray)
%           GlobalMax     global maximum of the variable
%           
%        Attributes  structure array of attributes
%           Name          path of attribute
%           Type          Matlab type class of data
%           Value         attribute value
%
%
%   INFO = ADIOSOPEN(FILE, 'Verbose', LEVEL)
%      To get logging from the adiosopen code, set Verbose to 1 or higher.
%      Higher values cause more and more details to be printed.
%
%   See also ADIOSREAD, ADIOSCLOSE, ADIOS.

%   Copyright 2009 UT-BATTELLE, LLC
%   $Revision: 1.0 $  $Date: 2009/08/05 12:53:41 $
%   Author: Norbert Podhorszki <pnorbert@ornl.gov>

%
% Process arguments.
%

checkArgCounts(varargin{:});
[args, msg] = parse_inputs(varargin{:});
if (~isempty(msg))
    error('MATLAB:adiosopen:inputParsing', '%s', msg);
end

if (isnumeric(args.File))
    fn=sprintf('File handler=%lld',args.File);
else
    fn=sprintf('File name=%s',args.File);
end
verbose=sprintf('%d ', args.Verbose);

input = sprintf('adiosopenc:\n  %s \n  Verbose=%s', fn, verbose);
if (args.Verbose > 0) 
    CallArguments = input
end


if (nargout == 1)
    info = adiosopenc(args.File, args.Verbose);
else 
    adiosopenc(args.File, args.Verbose);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION:   checkArgCounts %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function checkArgCounts(varargin)

if (nargin < 1)
    error('MATLAB:adiosopen:notEnoughInputs', ...
          'ADIOSOPEN requires at least one input argument.')
end


if (nargout > 1)
    error('MATLAB:adiosopen:tooManyOutputs', ...
          'ADIOSOPEN requires one or zero output arguments.')
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION:   parse_inputs   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [args, msg] = parse_inputs(varargin)

args.File    = '';
args.Verbose = 0;  % verbosity, default is off, i.e. 0

msg = '';

% Arg 1: file name
if ischar(varargin{1})
    args.File = varargin{1};
    %
    % Verify existence of filename
    %
    fid = fopen(args.File);    
    if (fid == -1)
        % Look for filename with extensions.
        fid = fopen([args.File '.bp']);
    end

    if (fid == -1)
        error('MATLAB:hdf5read:fileOpen', ...
              'Couldn''t open file (%s).', ...
              args.File)
    else
        % Get full filename 
        args.File = fopen(fid);
        fclose(fid);
    end
    varargin = {varargin{2:end}};
else
    msg = 'FILE input argument to ADIOSOPEN must be a string ';
    return
end


% Parse optional arguments based on their number.
if (length(varargin) > 0)
    
    paramStrings = {'verbose'};
    
    % For each pair
    for k = 1:2:length(varargin)
        param = lower(varargin{k});
            
        if (~ischar(param))
            msg = 'Parameter name must be a string.';
            return
        end
        
        idx = strmatch(param, paramStrings);
        
        if (isempty(idx))
            msg = sprintf('Unrecognized parameter name "%s".', param);
            return
        elseif (length(idx) > 1)
            msg = sprintf('Ambiguous parameter name "%s".', param);
            return
        end

        switch (paramStrings{idx})
        % VERBOSE
        case 'verbose'
            if (k == length(varargin))
                msg = 'No value specified for Verbose option.';
                return
            end
        
            args.Verbose = varargin{k+1};
            if ((~isnumeric(args.Verbose)) || ...
                (~isempty(find(rem(args.Verbose, 1) ~= 0))))
                
                msg = sprintf('''VERBOSE'' must be an integer.');
                return
            end
            if (args.Verbose < 0)
                msg = sprintf('''VERBOSE'' must be greater or equal to zero.');
                return
            end
        end
    end
end
