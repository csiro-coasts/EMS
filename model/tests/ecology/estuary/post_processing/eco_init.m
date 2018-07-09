function out = eco_init
% ECO_INIT  Initialises the ecology post process structure
%

s.name       = '';
s.start_file = '';
s.filesRange = [];
s.files = {}; % cell array or sting with range
s.vars  = [];
s.t     = [];

% Substruct of regions
s.region.name   = '';
s.region.Irange = -1;
s.region.Jrange = -1;
s.region.Zrange = {}; % cell array of depth-ranges

% What function to perform, eg. mean, depth integration etc...
s.function = 'mean';

% Internal states. i.e. calculated values
s.region.vols       = {};
s.region.total_vols = {};
s.region.data       = {};
s.region.area       ={};

% Assign output var
out = s;

% EOF
