function out = eco_post_process(s)

% ECO_POST_PROCESS : Ecology post process SHOC output files
% 
% Input : initialised structure, see eco_init
%
% xxx porosity for sediments??
%
out = [];

%%%%%%%%%%%%%%%
% SETUP FILES %
%%%%%%%%%%%%%%%

% Quiet benign, ncmex warnings
ncquiet;

% Construct the fullfile cell array
if ~iscell(s.files)
  s.files = getFiles(s.start_file, s.filesRange);
end

% Open first file
fid = ncmex('open', s.files{1} , 'nowrite');
if fid <= 0
  error(['Unable to open file : ', s.files{1}]);
end

% Loop over each region to compute volumes
for i=1:length(s.region)
  %
end

% Get time
s.t = ncmex('varget', fid, 't', 0, [-1], 0)'; % make col vector

% Done with this file
ncmex('close', fid);

%%%%%%%%%%%%%
% MAIN-LOOP %
%%%%%%%%%%%%%

% Loop over each file
for i=1:length(s.files)
  fid = ncmex('open', s.files{i} , 'nowrite');
  if fid <= 0
    error(['Unable to open file : ', s.files{i}]);
    %
    % optionally make this a warning and continue
    %
  else
    disp(['Starting to process file: ', s.files{i}]);
  end
  
  % over each variable
  for l=1:length(s.vars)
    var = s.vars{l};

    disp(['Var : ', var]);
    
    % Over each region
    for r=1:length(s.region);
      reg = s.region(r);

      % Setup J,I dimensions
      [begin, count] = getBegCntJI(reg);
      
      % Compute volumes
      % if isempty(findstr(var ,'_sed'))
        reg = computeVolumes(reg, fid, begin, count);
      % end
      disp(['Region ', num2str(r), ' out of ', ...
            num2str(length(s.region))]);
      
      % Make 4-D
      % Assume all time go all the way to the bottom of the column
      begin = [ 0,  0, begin];
      count = [-1, -1, count];
      
      % Get the variable in this region
      v = ncmex('varget', fid, var, begin, count, 0);
      porosity=ncmex('varget', fid, 'porosity_sed', begin, count, 0);
      if v<0
        error(['Variable ',var,' not found in file ',s.files{i},'!']);
      end
      varid = ncmex('varid',fid,var);
      
      % Do function
      switch s.function,
       case 'mean'
        % Each cell contains data for each file
        % Each column in data contains the timeseries means for
        %      each variable
        for z=1:length(reg.Zrange)
          k_z   = reg.Zrange{z};
          vdata = v(:,:,k_z,:);
          data  = calc_sums(reg.vols{z}, vdata) / reg.total_vols{z};
          % Fill the region struct
          reg.data{i}(:,l,z) = data;                    
        end

       case 'sum'
        % Each cell contains data for each file
        % Each column in data contains the timeseries means for
        %      each variable
        
        
        for z=1:length(reg.Zrange)
          if isempty(or(findstr(var ,'MA'),findstr(var ,'SG')) )
              
              if isempty(findstr(var ,'_sed'))
                k_z   = reg.Zrange{z};
                vdata = v(:,:,k_z,:);
                data  = calc_sums(reg.vols{z}, vdata);
             
              else
                  
                data = calc_sums_sed(fid, v, varid,porosity);
                
              end
              % Fill the region struct
             % reg.data{i}(:,l,z) = data;  
              
          else
             k_z   = reg.Zrange{z};
             vdata = v(:,:,:);
             data = calc_sums_epi(reg.area{z}, vdata);
            %data = calc_sums_sed_epi(fid, v, varid);

          end   
        reg.data{i}(:,l,z) = data;      
        end
        
       otherwise
        error([s.function, ' not supported!']);
      end
      
      % Replace reg
      s.region(r) = reg;
      
    end
  end
  % We're done with this file
  ncmex('close', fid);
end

% Assign output
out = s;

% end eco_post_process

%%%%%%%%%%%%%%%%%%%%%%%%
% Local help functions %
%%%%%%%%%%%%%%%%%%%%%%%%

%
% Expand into a file cell array
%
function files = getFiles(start_file, file_range)

% sanity checking
if isempty(start_file)
  error('start_file not supplied');
end
if isempty(file_range)
  error('file_range seems to be empty');
end

len   = length(file_range);
files = cell(1, len);

% Loop over and 
for i=1:len
  files{i} = [start_file, num2str(file_range(i)),'.nc'];
end

% end getFiles

%
% Computes the volumes of each region using the netcdf filename
% provided.
%
function s = computeVolumes(s, fid, begin, count)

% Extract the width and hiegts of cells using the specified 
% I's and J's.
h1   = ncmex('varget', fid, 'h1acell', begin, count, 0);
h2   = ncmex('varget', fid, 'h2acell', begin, count, 0);
botz = ncmex('varget', fid, 'botz',    begin, count, 0);

% Calculate areas
areas = h1 .* h2;

% Zero out the land cells
areas(find(isnan(areas))) = 0;

% Get the layer depths and the differences
% Z  = ncmex('varget', fid, 'z_centre', 0, -1, 0);
% dZ = diff([Z 0]'); % Add zero to diff against the surface layer
% Use actual layerfaces
Z  = ncmex('varget', fid, 'z_grid', 0, -1, 0);
Z(end) = [];
dZ = diff([Z 0]'); % Add zero to diff against the surface layer

% Preallocate memory
dZZ = zeros([size(h1) length(dZ)]);

% Iterate over each column and clamp 
len_dZ = length(dZ);
for i=1:size(h1,1)
  for j=1:size(h1,2)
    bot_k  = find(Z > botz(i,j));
    col_dz = [zeros(len_dZ - length(bot_k), 1); dZ(bot_k)];
    dZZ(i,j,:) = col_dz;
  end
end

if isempty(s.Zrange)
  s.Zrange = {(1:length(Z))};
end

% Divvy up the depths
for z=1:length(s.Zrange)
  % Note s.Zrange must be indicies
  k_z = s.Zrange{z};
  
  % Expand the areas
  area = repmat(areas, [1 1 length(dZ(k_z))]);
  
  % Calculate the volume on a per cell basis and cache
  vol = area .* dZZ(:,:,k_z);
  
  % Cache stuff
  s.vols{z}       = vol;
  s.total_vols{z} = sum(sum(sum(vol)));
  s.area{z}=area;
end

% end computeVolumes

%
% Gets the begin and counts given the range
%
function [begin, count] = getBegCntJI(s)

begin = [s.Jrange(1) s.Irange(1)];
count = [s.Jrange(end) - s.Jrange(1)+1, ...
         s.Irange(end) - s.Irange(1)+1];

% override if defaults
if s.Jrange == -1
  begin(1) = 0;
  count(1) = -1;
end
if s.Irange == -1
  begin(2) = 0;
  count(2) = -1;
end

% end getBegCntJI

%%%%%%%%%%%%%
% FUNCTIONS %
%%%%%%%%%%%%%

% Calculate the total sums
function out = calc_sums(vol, v)

out = [];

t_len  = size(v, 4);
vol_4d = repmat(vol, [1 1 1 t_len]);

% Replace nan's with 0
v(find(isnan(v))) = 0;

% xxx
v(find(v<0))=0;


% Multiply out
% i.e. multiply each cell's concentration with its volume
mass = v .* vol_4d;

% Calculate mean
% The first 3 dimenstions are in space
mn = sum(sum(sum(mass)));

% Make 1-D column vector
mn = shiftdim(mn);

% Assign output
out = mn;

% end calc_sums


% Calculate the total sums in sed
function out = calc_sums_sed(fid, v,varid,porosity)
% This works a little different as we have a full 4D volumes for sediments
% 
out = [];

begin = [0, 0];
count = [-1, -1];

% Get the dissolved attribute to know if we need to apply porosity
% in the mass calculation
dissol = ncmex('attget', fid, varid, 'dissol');

% Replace nan's with 0
porosity(find(isnan(porosity))) = 0;
porosity(find(porosity<0))=0;

if dissol==1
    v=v.*porosity;
end
% Extract the width and hiegts of cells using the specified 
% I's and J's.
h1   = ncmex('varget', fid, 'h1acell',  begin, count, 0);
h2   = ncmex('varget', fid, 'h2acell', begin, count, 0);

% Compute the areas as a 2D plane
areas = h1 .* h2;

% Zero out the land cells
areas(find(isnan(areas))) = 0;

% Add the more dimensions to these vector
begin = [ 0,  0, 0, 0];
count = [-1, -1, -1, -1];

% Replace nan's with 0
v(find(isnan(v))) = 0;
v(find(v<0))=0;


z_grid_sed = ncmex('varget', fid, 'z_grid_sed', begin, count, 0);
z_grid_sed(find(isnan(z_grid_sed))) = 0;
dZ_sed = diff(z_grid_sed, [], 3);

ar = repmat(areas, [1 1 size(dZ_sed, 3)]);

% Do this to save memory
% Maybe we should stride this?
rec_len = size(dZ_sed, 4);
for i=1:rec_len
  dz = dZ_sed(:,:,:,i);
  % 3D volumes
  vols = ar .* dz;
  % Tracer sums
  tr = sum(sum(sum((vols .* v(:,:,:,i)))));
  % Assign output
  out(i,1) = tr;
end

% end calc_sums_sed



% Calculate the total sums
function out = calc_sums_epi(vol, v)

out = [];

t_len  = size(v, 3);
vol_4d = repmat(vol, [1 1 1 t_len]);

% Replace nan's with 0
v(find(isnan(v))) = 0;

% xxx
v(find(v<0))=0;


% Multiply out
% i.e. multiply each cell's concentration with its volume
mass = v .* squeeze(vol_4d(:,:,1,:));

% Calculate mean
% The first 3 dimenstions are in space
mn = sum(sum((mass)));

% Make 1-D column vector
mn = shiftdim(mn);

% Assign output
out = mn;

% end calc_sums_epi












function out = calc_sums_sed_epi(fid, v,varid)
% This works a little different as we have a full 4D volumes for sediments
% 
out = [];

begin = [0, 0];
count = [-1, -1];

% Get the dissolved attribute to know if we need to apply porosity
% in the mass calculation

% Replace nan's with 0

% Extract the width and hiegts of cells using the specified 
% I's and J's.
h1   = ncmex('varget', fid, 'h1acell', begin, count, 0);
h2   = ncmex('varget', fid, 'h2acell', begin, count, 0);

% Compute the areas as a 2D plane
areas = h1 .* h2;

% Zero out the land cells
areas(find(isnan(areas))) = 0;

% Add the more dimensions to these vector
begin = [ 0,  0, 0, 0];
count = [-1, -1, -1, -1];

% Replace nan's with 0
v(find(isnan(v))) = 0;
v(find(v<0))=0;


%z_grid_sed = ncmex('varget', fid, 'z_grid_sed', begin, count, 0);
%z_grid_sed(find(isnan(z_grid_sed))) = 0;
%dZ_sed = diff(z_grid_sed, [], 3);

ar = areas;%repmat(areas, [1 1 size(dZ_sed, 3)]);
size(ar)
% Do this to save memory
% Maybe we should stride this?
%rec_len = size(dZ_sed, 4);
%for i=1:rec_len
%  dz = dZ_sed(:,:,:,i);
  % 3D volumes
%  vols = ar .* dz;
  % Tracer sums
  [a,b,c]=size(v)
  for i=1:c
  tr(i,:,:) = sum(sum(((ar.* v(:,:,i)))));
  end
  'tr'
  size(tr)
  % Assign output
  out = tr;
%end

% end calc_sums_sed

