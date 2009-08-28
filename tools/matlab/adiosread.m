function [data, attributes] = adiosread(varargin)
%ADIOSREAD Read data from an ADIOS BP file.
%   
%   ADIOSREAD can read in data with different call options.
%   The simplest way is to give the name of the file and the path
%   to a variable/attribute. If you have more than one adios groups
%   stored in one file, you should supply the group name as well.
%   In these cases, the library opens the file, reads in the data
%   and closes the file. 
%
%   If you have opened the file with ADIOSOPEN, you can supply the
%   file handler and the group handler to read in a variable/attribute.
%   When finished reading all data, you need to close the file with 
%   ADIOSCLOSE. 
%
%   Note: the adios group is not the same as HDF5 groups. Each variable has
%   a path, which defines a logical hierarchy of the variables within one 
%   adios group. This logical hierarchy is what is similar to HDF5 groups.
%
%   DATA = ADIOSREAD(FILE, VARPATH) 
%      Read complete DATA from FILE with path VARPATH.
%      Use only if you have only one adios group in the file or you know that
%      the VARPATH is in the first adios group in the file.
%
%   DATA = ADIOSREAD(FILE, GROUP, VARPATH) 
%      Read complete DATA from an adios GROUP of FILE with path VARPATH.
%      Use this form if you have more than one ADIOS group in the file.
%      GROUP is either a group name (string) or an index of the groups
%      (integer, starting from 1) or a GROUP handler from ADIOSOPEN (int64).
% 
%   ATTR = ADIOSREAD(FILE, ATTRPATH) 
%   ATTR = ADIOSREAD(FILE, GROUP, ATTRPATH) 
%      Read in an attribute similarly as you read in variables.
%      ATTR is just a variable as DATA variables in ADIOS, although only
%      either strings or 1-by-1 scalars and do not vary over time.
% 
%   DATA,ATTRS = ADIOSREAD(FILE, VARPATH) 
%   DATA,ATTRS = ADIOSREAD(FILE, GROUP, VARPATH) 
%      Read in a variable and all attributes 'belonging' to that variable.
%      An attribute belongs to a variable if ATTRPATH=VARPATH/NAME.
%      ATTRS is an array of structure where a structure contains the Name 
%      and Value of an attribute.
% 
%   DATA... = ADIOSREAD( ..., 'Time', TIME)
%      ADIOS files can contain several timesteps of data. In this case you 
%      can specify the timestep you want to read in (default is the first
%      timestep). Time indices start from 1 but since an ADIOS file can be 
%      a split part containing only an interval, it may start from a larger 
%      number (see ADIOSOPEN structure TimeOffsets field).
%      A negative index is interpreted as counting from the last time.
%      -1 refers to the last timestep in the file, -2 the previous one, etc.
%
%   DATA... = ADIOSREAD( ..., 'Slice', SLICEDEF)
%      You can read a portion of a variable. 
%      A slice is defined as an N-by-2 array of integers, where N is the 
%      number of dimensions of the variable (or less). A tuple describes the 
%      "start" and "count" values. The "start" values start from 1.
%          E.g. [1 10; 5 2] reads the first 10 values in the first dimension
%      and 2 values from the 5th value in the second dimension resulting in
%      a 10-by-2 array. 
%          You can use negative numbers to index from the end of the array
%      as in python. -1 refers to the last element of the array, -2 the one
%      before and so on. 
%          E.g. [-1 1] reads in the last value of a 1D array. 
%               [1 -1] reads in the complete dimension.
%      If the slice definition has less rows than the number of dimensions
%      of the variable, [1 -1] rows are added automatically to read those
%      dimensions completely.
%          If the slice definition has more rows than the number of dimensions
%      of the variable, the extra slice definitions will be ignored.
%
%   DATA... = ADIOSREAD( ..., 'Verbose', VALUE)
%      To get logging from the adiosread code, set Verbose to 1 or higher.
%      Higher values cause more and more details to be printed. 
%
%   Please read the file adioscopyright.txt for more information.
%
%   See also ADIOSOPEN, ADIOSCLOSE, ADIOS.

%   Copyright 2009 Oak Ridge National Laboratory
%   $Revision: 1.0 $  $Date: 2009/08/05 12:53:41 $
%   Author: Norbert Podhorszki <pnorbert@ornl.gov>

%
% Process arguments.
%

checkArgCounts(varargin{:});
[args, msg] = parse_inputs(varargin{:});
if (~isempty(msg))
    error('MATLAB:adiosread:inputParsing', '%s', msg);
end

if (isnumeric(args.File))
    fn=sprintf('File handler=%lld',args.File);
else
    fn=sprintf('File name=%s',args.File);
end
if (isnumeric(args.Group))
    if (isa(args.Group, 'int64'))
        gn=sprintf('Group handler=%lld',args.Group);
    else
        gn=sprintf('Group index=%d',args.Group);
    end
else
    gn=sprintf('Group name=%s',args.Group);
end
offsets=sprintf('%d ', args.Offsets);
counts=sprintf('%d ', args.Counts);
verbose=sprintf('%d ', args.Verbose);

input = sprintf('adiosreadc:\n  %s \n  %s \n  Var=%s\n  Time=%d  Offsets=[%s]  Counts=[%s] \n  Verbose=%s', ...
 fn, gn, args.Path, args.Time, offsets, counts, verbose);
if (args.Verbose > 0) 
    CallArguments = input
end


if (nargout == 1)
    data = adiosreadc(args.File, args.Group, args.Path, args.Time, args.Offsets, args.Counts, args.Verbose);
elseif (nargout == 2)
    [data, attributes] = adiosreadc(args.File, args.Group, args.Path, args.Time, args.Offsets, args.Counts, args.Verbose);
else 
    adiosreadc(args.File, args.Group, args.Path, args.Time, args.Offsets, args.Counts, args.Verbose);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION:   checkArgCounts %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function checkArgCounts(varargin)

if (nargin < 3)
    error('MATLAB:adiosread:notEnoughInputs', ...
          'ADIOSREAD requires at least two input arguments.')
end


if (nargout > 2)
    error('MATLAB:adiosread:tooManyOutputs', ...
          'ADIOSREAD requires two or fewer output arguments.')
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION:   parse_inputs   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [args, msg] = parse_inputs(varargin)

args.File    = '';
args.Group   = '';
args.Path    = '';
args.Time    = int32(-1); % time, default is 'last timestep', i.e. -1
args.Offsets = []; % start positions in each dimension for slicing
args.Counts  = []; % counts in each dimension for slicing
args.Verbose = 0;  % verbosity, default is off, i.e. 0

msg = '';

% Arg 1: file name or int64 file handler
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
elseif (isa(varargin{1}, 'int64'))
    args.File = varargin{1};
else
    msg = ['FILE input argument to ADIOSREAD must be a string ' ...
           'or an int64 handler from ADIOSOPEN'];
    return
end

% Arg 2 and maybe 3: group and varpath
if (rem(nargin, 2) == 0)
    % even number of arguments: FILE, VARPATH, ...
    args.Path = varargin{2};
    varargin = {varargin{3:end}};
else
    % odd number of arguments: FILE, GROUP, VARPATH, ...
    args.Group = varargin{2};
    args.Path = varargin{3};
    varargin = {varargin{4:end}};
    % check type of Group
    if ((~ischar(args.Group)) && ...
        (~isnumeric(args.Group)) && ...
        (~isa(args.Group, 'int64'))) 
        msg = ['GROUP input argument to ADIOSREAD must be a string ' ...
           'or an number or an int64 handler from ADIOSOPEN'];
        return
    end
    if (isnumeric(args.Group) && (~isa(args.Group, 'int64')))
        % convert group index to int32
        args.Group = int32(args.Group);
    end
end
% check type of Path
if ~ischar(args.Path)
    msg = ['VARPATH input argument to ADIOSREAD must be string.'];
    return
end

% Parse optional arguments based on their number.
if (length(varargin) > 0)
    
    paramStrings = {'time', 'slice', 'verbose'};
    
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
        % TIME
        case 'time'
            if (k == length(varargin))
                msg = 'No time value specified for Time option.';
                return
            end
        
            time = varargin{k+1};
            if ((~isnumeric(time)) || ...
                (~isempty(find(rem(time, 1) ~= 0))))
                
                msg = sprintf('''TIME'' must be an integer.');
                return
            end
            if (time == 0)
                msg = sprintf('''TIME'' in adios files start from 1, not from 0 .');
                return
            end
            args.Time = int32(fix(time));

        % SLICE
        case 'slice'
            if (k == length(varargin))
                msg = 'No slicing value specified for Slice option.';
                return
            end
        
            slices = varargin{k+1};
            if ((~isnumeric(slices)) || ...
                (~isempty(slices) && size(slices, 2) ~= 2) || ...
                (~isempty(find(rem(slices, 1) ~= 0))))

                msg = 'Slice values must be n-by-2 array of integers.';
                return
            end

            if (~isempty(slices))
                args.Offsets = int32(fix(slices(:,1)));
                args.Counts = int32(fix(slices(:,2)));
            else
                args.Offsets = int32([]);
                args.Counts = int32([]);
            end
            
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
