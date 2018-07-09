
%	plotts
%
%------------------------------------------------------------------------------
%
% 	Plots time series from ascii ts files, or netcdf SHOC/BM files.
%
%------------------------------------------------------------------------------
%
%   Usage :
%
%	Plotts should be called after defining the necessary input parameters.
%	The following parameters are compulsory :-
%
%		inpfile  : An array of input file names.
%		var_name : An array of plot variable names.
%		file_no  : A vector, of length equal to var_name, specifying
%			   which input files to take the variables from.
%		cell_no  : Only compulsory if netcdf files are being used.
%			   An array of cell indices of height equal to the
%			   the no. of variables from netcdf files, and width
%			   equal to 2 for BM files and equal to 3 for SHOC
%			   files. 
%
%	All other plot parameters will default if not specified.
%	An example of a script which uses plotts is given below.
%
%------------- Example script --------------------------------------------------
%
%	header = 'ANM  (P3 bot)' ;
%	plot_start = [1990 01 01 00 00 00] ;
%	plot_end   = [1990 01 25 24 00 00] ;
%	landscape = 0 ;
%	nplotperpage = 4 ;
%
%	inpfile = ['g00/run61/P3_bot.ts              '
%		   'g00/run62/P3_bot.ts              '
%		   'g00/run63/P3_bot.ts              '
%		   '/mgk/andrewar/anm99/run61_surf.nc'] ;
%
%	var_name = ['eta             '
%                   'u v             '
%                   'salt            '
%                   'salt            '
%                   'BW_bottom_tracer'
%                   'eta             '] ;
%	cell_no   = [60 15 15] ;
%	file_no   = [  1   1  1  3  2   4] ;
%	plot_no   = [  1   2  3  3  4   1] ;
%	stick     = [  0   1  0  0  0   0] ;
%	decim     = [  1   1  1  1  1   7] ;
%	scalefac  = [  1   1  1  1  1   1] ;
%	fill2nan  = [  0   0  0  0  0   0] ;
%	miss2nan  = [  1   1  1  1  1   1] ;
%	linewidth = [ .5  .5 .5 .5 .5  .5] ;
%	linemark  = [  1   1  1  1  1   5] ;
%	linetype  = ['b- ' ; 'b- ' ; 'c- ' ; 'm- ' ; 'b- '; 'r^ '] ;
%	font1 = 13 ; font2 = 11 ; font3 = 10 ;
%
%	yclip     = [  0   1  0  0 ] ;
%	yaxis_min = [-.5 -.02  5  5] ;
%	yaxis_max = [ .5 .02 30 30 ] ;
%	ytick_min = [-.5 -.02  5  5] ;
%	ytick_max = [ .5 .02 30 30] ;
%	nytick    = [  3   3  6  6] ;
%	zero_line = [  1   1  1  1] ;
%	plotbot    = [.2 .4 .6 .8] ;
%	plotheight = [.15 .15 .15 .15] ;
%	y_label  = ['Sea-Level       '
%                    'Current         '
%                    'Salinity        '
%                    'BW Bottom Tracer'] ;
%	y_units  = ['(m)      '
%                    '(m s^-^1)'
%                    '(PSU)    '
%                    '         '] ;
%	ptitle = ['          Run61            '
%		  '  Run 61 - current sticks  '
%		  'Runs 61 (blue) and 63 (red)'
%		  '         Run 62            '] ;
%	subtitle = ['                               '
%		    '                               '
%		    '90 cumecs (Neap & Spring Tides)'
%		    '                               '] ;
%	plotts ; pause ;
%	print -dpsc plot1.ps ;
%
%	inpfile = ['g00/run61/Q3_bot.ts              '
%		   'g00/run62/Q3_bot.ts              '
%		   'g00/run63/Q3_bot.ts              '
%		   '/mgk/andrewar/anm99/run61_surf.nc'] ;
%	header = 'ANM  (Q3 bot)' ;
%	plotts ; pause ;
%	print -dpsc plot2.ps ;
%
%------------------------------------------------------------------------------

	clf ;
	global maxfigl maxfigr ;
%
%...... Initialise.
%
	numnc = 0 ;
	labside = ['left '] ;
	month = ['Jan';'Feb';'Mar';'Apr';'May';'Jun';...
	         'Jul';'Aug';'Sep';'Oct';'Nov';'Dec'] ;
%
%...... Set paper orientation.
%
	set(gcf,'Color',[1 1 1]) ;
	set(gcf,'PaperUnits','inches') ;
	if (~exist('landscape') | ~landscape),
	    set(gcf,'PaperOrientation','portrait') ;
	    set(gcf,'PaperPosition',[.25 0.5 8.0 11.5]) ;
%	    set(gcf,'PaperPosition',[0.5 0.5 6.0 2.0]) ;   % gipp ts
%	    set(gcf,'PaperPosition',[.25 1.0 8.0 6.0]) ;   % gipp talk
	else
	    set(gcf,'PaperOrientation','landscape') ;
	    set(gcf,'PaperPosition',[0.0 0.5 11.0 8.0]) ;
	end
%
%...... Check there is at least 1 input file listed.
%
	if (~exist('inpfile')),
	    error('** Need an array of input file names (inpfile) **') ;
	else
	    [ninp,n] = size(inpfile) ;
	end
%
%...... Ensure at least 1 file has been requested.
%
	if (~exist('file_no')),
	    error('** Need a vector of file nos. (file_no) **') ;
	end
%
%...... Ensure all file nos are legal.
%
	if (max(file_no) > ninp),
	    error('** Illegal file_no specified **') ;
	end
%
%...... Sort file nos. into ascending order, to allow more efficient
%...... reading/use of ascii data in ts files.
%
	[dum,new_order] = sort(file_no) ;
%
%...... Determine how many separate files being used to plot.
%
	nfile = sum(diff(new_order)) + 1 ;
%
%...... Determine how many individual variables being plotted.
%
	nvar = length(file_no) ;
%
%...... Determine whether its to be a 'single period' plot, or
%...... a 'multiple period' plot.
%
	if (exist('plot_start') & exist('plot_end')),
	    [nperiod,n] = size(plot_start) ; [m,n] = size(plot_end) ;
	    if (nperiod ~= m),
	        error('** Mismatch in no. of plot start & end times **') ;
	    else
		if (nperiod>1),
	    	    mult_period = 1 ;
		    if (nperiod < nvar),
	    		error('** Insufficient plot starts specified **') ;
		    end
		else
	    	    mult_period = 0 ;
		end
	    end
	else
	    if (nfile>1),
	    	 mult_period = 1 ;
	    else
	    	 mult_period = 0 ;
	    end
	end
%
%...... Assign defaults for most variable-related parameters that do not exist.
%
	if (~exist('var_name')),
	    error('** Need an array of variable names (var_name)**') ;
	else
	    [m,n] = size(var_name) ;
	    if (m<nvar),error('** Insufficient variable names specified **');end
	end
	if (~exist('plot_no')),
	    disp('** Warning : default plot nos. being used **') ;
	    plot_no = [1:1:nvar] ;
	else
	    if (length(plot_no)<nvar),
	        error('** Insufficient plot nos. specified **') ;
	    end
	end
	if (~exist('stick') | length(stick)<nvar), stick = 0*ones(1,nvar) ; end 
	if (~exist('decim') | length(decim)<nvar), decim = ones(1,nvar) ; end 
	if (~exist('scalefac') | length(scalefac)<nvar),
	    scalefac = 1*ones(1,nvar) ;
	end 
	if (~exist('fill2nan') | length(fill2nan)<nvar),
	    fill2nan = 0*ones(1,nvar) ;
	end 
	if (~exist('miss2nan') | length(miss2nan)<nvar),
	    miss2nan = 0*ones(1,nvar) ;
	end 
	if (~exist('linewidth') | length(linewidth)<nvar),
	    linewidth = .5*ones(1,nvar) ;
	end 
	if (~exist('linemark') | length(linemark)<nvar),
	    linemark = ones(1,nvar) ;
	end 
	if (~exist('linetype')), linetype = ['k- '] ; end
	if (~exist('font1')), font1 = 13 ; end
	if (~exist('font2')), font2 = 11 ; end
	if (~exist('font3')), font3 = 10 ; end
%
%...... Determine how many separate plot boxes being plotted.
%...... & how many variables per plot box.
%
	nplot = max(plot_no) ; done_plot = 0 * ones(nplot,1) ;
	for i = 1:nplot,
	    ind = find(plot_no == i) ;
	    if (isempty(ind)),
	        error('** Missing plot no. **') ;
	    else
		nplotperbox(i) = length(ind) ;
	    end
	end
	if (~exist('nplotperpage')), nplotperpage = nplot ; end
%
%...... Assign defaults for most plot-box-related parameters that do not exist.
%
	if (exist('plotbot')),
	   m = length(plotbot) ;
	   if (m<nplot), error('** Insufficient plot bots. specified **') ; end
	end
	if (exist('plotheight')),
	    m = length(plotheight) ;
	    if (m<nplot),error('** Insufficient plot heights specified **') ; end
	end
	if (~exist('yclip') | length(yclip)<nplot),
	    yclip = ones(1,nplot) ;
	end 
	if (~exist('zero_line') | length(zero_line)<nplot),
	    zero_line = ones(1,nplot) ;
	end 
	if (~exist('nytick') | length(nytick)<nplot),
	    nytick = 2*ones(1,nplot) ;
	end 
	if (exist('y_label')),
	    [m,n] = size(y_label) ;
	    if(m<nplot),y_label(m+1:nplot,:)=repmat(y_label(m,:),nplot-m,1) ; end
	end
	if (exist('y_units')),
	    [m,n] = size(y_units) ;
	    if(m<nplot),y_units(m+1:nplot,:)=repmat(y_units(m,:),nplot-m,1) ; end
	end
	if (exist('ptitle')),
	    [m,n] = size(ptitle) ;
	    if (m<nplot), error('** Insufficient titles specified **') ; end
	end
	if (exist('subtitle')),
	    [m,n] = size(subtitle) ;
	    if (m<nplot), error('** Insufficient subtitles specified **') ; end
	end
%
%...... If plot box positions and heights not given, calculate
%...... a default (in normalised page units).
%
	plotwidth = 0.7 ;
	plotleft = 0.7*(1.-plotwidth) ;
	if (~exist('plotbot') | ~exist('plotheight')),
	    botmargin = 0.05 ; topmargin = 0.02 ;
	    headergap = 0.05 ;
	    titlegap = 0.04 ;
	    if (nplotperpage < 5), titlegap = 0.05 ; end
	    if (nplotperpage < 4), titlegap = 0.06 ; end
	    if (nplotperpage < 3), titlegap = 0.07 ; end
	    dategap = 0 ; if (mult_period), dategap = 0.04 ; end
	    plotgap = (titlegap + dategap) * ones(1,nplotperpage) ;
	    plotgap(1) = plotgap(1) + headergap ;
	    totgap = sum(plotgap) ;
	    pheight = ((1.0-botmargin-topmargin-totgap)/nplotperpage) ;
	    plotheight = pheight * ones(1,nplotperpage) ;
	    cumgap = cumsum(plotgap) ;
	    plotbot = 1.0 - topmargin - cumgap -[1:1:nplotperpage]*pheight ;
	end
%
%---------------- Loop over the variables in the sorted order ------------------
%
	iplot = 0 ;
	for ivar = new_order,
	    iplot = iplot + 1 ;
%
%.......... Check whether input file is ascii ts or netcdf.
%
	    inpstr = deblank(inpfile(file_no(ivar),:)) ;
	    len = length(inpstr) ;
	    if (inpstr(len-1:len) == 'ts' | inpstr(len-1:len) == 'es'),
		ts = 1 ; nc = 0 ;
	    elseif (inpstr(len-1:len) == 'nc'),
		nc = 1 ; ts = 0 ;
	    else
		error('** Unrecognizeable file suffix **') ;
	    end
%
%.......... Find all variables & constants in the var_name string.
%
	    [vars] = split_maths(deblank(var_name(ivar,:))) ;
	    [numvar,n] = size(vars) ;
%
%.......... Set some defaults.
%
	    longname = [] ; unitstring = [] ;
%
%.......... Open and read from an ascii ts file.
%
	    if (ts),
		if (iplot==1 | file_no(ivar)~=file_no(new_order(iplot-1))),
	            [head,alldata] = hdrload(inpstr) ;
		end
	        [mh,nh] = size(head) ;
	        head_row = reshape(head',1,mh*nh) ;
%
%.............. Extract time info., data, longname, units.
%
	        s = 'COLUMN1.units' ;
	        ind = findstr(head_row,s) + length(s) + 1 ;
		timestring = fliplr(deblank(fliplr(deblank(head_row(ind(1):ind(1)+nh-18))))) ;
		times = alldata(:,1) ;
		for jvar = 1:numvar,
	            shortname = deblank(vars(jvar,:)) ;
%.................. Check this is a variable & not a constant.
	    	    if (isempty(str2num(shortname))),
	            	ind = findstr(head_row,['name ' shortname]) ;
		    	if (isempty(ind)),
	            	    ind = findstr(head_row,['name  ' shortname]) ;
		    	    if (isempty(ind)),
			        error('** Variable not found **') ;
			    end
			end
		        s = head_row(ind(1):mh*nh) ;
		        ind = findstr(s,'COLUMN') + 6 ; ind2 = findstr(s,'.') ;
		        colno = str2num(s(ind(1):ind2(1)-1)) ;
			eval([shortname ' = alldata(:,colno) ;']) ;
		        ndigit = length((ind(1):ind2(1)-1)) ;
	    		if (fill2nan(ivar)),
		    	    ind = findstr(s,'fill_value ') + 11 ;
		    	    if (~isempty(ind)),
		                fillval = str2num(deblank(s(ind(1):ind(1)+nh-22-ndigit))) ;
				eval([shortname ' = change(' shortname ',''=='',fillval,NaN) ;']) ;
			    else
		    		disp ('** Warning : No fill_value available **') ;
		    	    end
		    	end
	    		if (miss2nan(ivar)),
		    	    ind = findstr(s,'missing_value') + 14 ;
		    	    if (~isempty(ind)),
		                missval = str2num(deblank(s(ind(1):ind(1)+nh-25-ndigit))) ;
				eval([shortname ' = change(' shortname ',''=='',missval,NaN) ;']) ;
			    else
		    		disp ('** Warning : No missing_value available **') ;
		    	    end
		    	end
			if (numvar==1),
		    	    ind = findstr(s,'long') + 9 ;
		    	    longname = deblank(s(ind(1):ind(1)+nh-21-ndigit)) ;
			end
		    	ind = findstr(s,'units') + 6 ;
		    	unitstring = deblank(s(ind(1):ind(1)+nh-17-ndigit)) ;
		    end
		end
	    end
%
%......... Open and read from a netcdf file.
%
	    if (nc),
		numnc = numnc + 1 ;
		if (~exist('cell_no')),
		    error('** Need a vector cell_no for netcdf files **') ;
		else
		    [m,nc] = size(cell_no) ;
		    if (m<numnc),
		        error('** Need cell indices for each netcdf file **') ;
		    end
		end
	        [cdfid] = ncmex('open',inpstr,'nowrite') ; ncquiet ;
	        [dum] = ncmex('varget',cdfid,'nominal_dz',[0 0],[-1 -1],0) ;
	        [m,nn] = size(dum) ;
	        if (nn>0),
		    shoc = 0 ; bm = 1 ;
		    if (nc<2),
		     error('** Need (i,j) cell indices for bm plotting **');
		    end
	        else
		    shoc = 1 ; bm = 0 ;
		    if (nc~=3),
		     error('** Need (i,j,k) cell indices for shoc plotting **');
		    end
	        end
%
%.............. Read in time info., data, longname, units.
%
	        timestring = ncmex('attget',cdfid,'t','units') ;
	        times = ncmex('varget',cdfid,'t',[0],[-1],0) ;
	        times = times' ;
	        i = cell_no(numnc,1) ;
	        j = cell_no(numnc,2) ;
	        if (shoc), k = cell_no(numnc,3) ; end
		for jvar = 1:numvar,
	            shortname = deblank(vars(jvar,:)) ;
%.................. Check a variable & not a constant.
	    	    if (isempty(str2num(shortname))),
		    	if (shoc),
	                    str = [shortname ' = ncmex(''varget'',cdfid,''' shortname ''',[0 k j i],[-1 1 1 1],0) ;'] ;
		    	else
	                    str = [shortname ' = ncmex(''varget'',cdfid,''' shortname ''',[0 i j],[-1 1 1],0) ;'] ;
		    	end
			eval(str) ; eval([shortname ' = squeeze(' shortname ') ;']) ;
	    		if (fill2nan(ivar)),
			    fillval = [] ;
		    	    [dtype,len,stat] = ncmex('attinq',cdfid,shortname,'_FillValue') ;
		    	    if (stat>0), fillval = ncmex('attget',cdfid,shortname,'_FillValue') ; end
		    	    [dtype,len,stat] = ncmex('attinq',cdfid,shortname,'_FillValueWC') ;
		    	    if (stat>0), fillval = ncmex('attget',cdfid,shortname,'_FillValueWC') ; end
			    if (~isempty(fillval)),
				eval([shortname ' = change(' shortname ',''=='',fillval,NaN) ;']) ;
			    else
		    		disp ('** Warning : No fill_value available **') ;
		    	    end
		    	end
	    		if (miss2nan(ivar)),
		    	    disp ('** Warning : No missing_value available **') ;
		    	end
			if (numvar==1),
	            	    longname = ncmex('attget',cdfid,shortname,'long_name') ;
			end
	            	unitstring = ncmex('attget',cdfid,shortname,'units') ;
		    end
		end
	        ncmex('close',cdfid) ;
	    end
%
%.......... Put the appropriate data into an array called 'data',
%.......... evaluating any maths. expression involved.
%
	    if (~stick(ivar)),
		eval(['data = ' deblank(var_name(ivar,:)) ' ;']) ;
	    else
		ind = findstr(' ',deblank(var_name(ivar,:))) ;
		if (numvar~=2 | isempty(ind)),
	    	    error('** Stick variables not specified correctly **') ;
		else
		    eval(['data = [' deblank(vars(1,:)) ' ' deblank(vars(2,:)) '] ;']) ;
	        end
	    end
%
%.......... Decimate data if required.
%
	    if (decim(ivar)>1),
		[m,n] = size(data) ;
	        decdata = data(1:decim(ivar):m,:) ;
		dectimes = times(1:decim(ivar):m) ;
	        clear data times ;
		data = decdata ; times = dectimes ;
		clear decdata ; dectimes ;
	    end
%
%.......... Scale data if required.
%
	    if (scalefac(ivar)~=1),
	        data = scalefac(ivar) * data ;
	    end
%
%.......... Put elapsed times into days.
%
	    if (timestring(1:4) == 'seco'), times = times/(3600*24) ; end
	    if (timestring(1:4) == 'minu'), times = times/(60*24) ; end
	    if (timestring(1:4) == 'hour'), times = times/24 ; end
%
%.......... Get start & end time/dates of data.
%
	    time_orig = time2num(timestring) ;
	    data_start = t2dat(time_orig,times(1)) ;
	    data_end   = t2dat(time_orig,times(length(times))) ;
%
%.......... Get start and end times of this plot.
%
	    if (exist('plot_start')),
		if (mult_period),
		    pstart = plot_start(ivar,:) ; pend = plot_end(ivar,:) ;
		else
		    pstart = plot_start(1,:) ; pend = plot_end(1,:) ;
		end
	    else
		pstart = [data_start(1:3) 0 0 0] ;
		if (data_end(4:6)==[0 0 0]),
		    pend = data_end ;
		else
		    pend = t2dat([data_end(1:3) 0 0 0],0.9999) ;
		end
	    end
%
%.......... Put elapsed times in days from start of plot.
%
	    times = times - dat2t(time_orig,pstart) ;
%
%.......... Get plot x-range (in days).
%
	    xmin = 0. ; xmax = dat2t(pstart,pend) ;
%
%.......... Truncate the data & times arrays to just that within range.
%
	    indt = find(times>=xmin & times<=xmax) ;
	    times = times(indt,:) ;
	    data = data(indt,:) ;
%
%.......... Test whether variable is a direction.
%
	    idir = 0 ;
	    ind = findstr(shortname,'Dir') ; ind2 = findstr(shortname,'dir') ;
	    if (~isempty(ind) | ~isempty(ind2)), idir = 1 ; end
%
%.......... If y-axis limits not given, calculate defaults.
%
	    if (~exist('yaxis_min') | ~exist('yaxis_max') | ...
	    length(yaxis_min)<nplot | length(yaxis_max)<nplot),
		if (idir),
		    yaxmin = 0 ; yaxmax = 360 ;
		    ytmin = 0 ; ytmax = 360 ;
		    ntick = 5 ; 
		else
		    datmin = min(min(data)) ;
		    datmax = max(max(data)) ;
		    datrange = datmax - datmin ;
		    logr = log10(datrange/2) ;
		    lexp = floor(logr) ;
		    lmant = logr - lexp ;
		    a = 10^(lmant) ;
		    b = [1 2 2.5 5.0 10.0] ;
		    c = abs(b-a) ;
		    [cmin ind] = min(c) ;
		    dy = 10^(lexp) * b(ind) ;
		    yaxmin = floor(datmin/dy)*dy ;
		    yaxmax = ceil(datmax/dy)*dy ;
		    ytmin = yaxmin ; ytmax = yaxmax ;
		    ntick = (ytmax-ytmin)/dy + 1 ;
		end
		yclip = ones(1,nplot) ; 
	    else
		yaxmin = yaxis_min(plot_no(ivar)) ; yaxmax = yaxis_max(plot_no(ivar)) ;
	        if (~exist('ytick_min') | ~exist('ytick_max') | ...
		length(ytick_min)<nplot | length(ytick_max)<nplot),
		    ytmin = yaxmin ; ytmax = yaxmax ;
		else
		    ytmin = ytick_min(plot_no(ivar)) ; ytmax = ytick_max(plot_no(ivar)) ;
		end
		ntick = nytick(plot_no(ivar)) ; 
	    end
%
%.......... For multi-period plots, ensure same horizontal scale.
%
	    if (ivar==1), plotwidth0 = plotwidth ; xmax0 = xmax ; end
	    if (mult_period), plotwidth = plotwidth0 * xmax / xmax0 ; end
%
%.......... Define, but not draw, axes and get scales.
%
	    if (yclip(plot_no(ivar))),
	        axes ('position',[plotleft plotbot(plot_no(ivar))...
		      plotwidth plotheight(plot_no(ivar))]) ;
		ylo = yaxmin ; yhi = yaxmax ;
	    else
		pbot = plotbot(plot_no(ivar)) - 1 ;
		pheight = plotheight(plot_no(ivar)) + 2 ;
	        axes ('position',[plotleft pbot plotwidth pheight]) ;
		scale = (yaxis_max(plot_no(ivar))-yaxis_min(plot_no(ivar)))/...
			   plotheight(plot_no(ivar));
		ylo = yaxmin - 1*scale ;
		yhi = yaxmax + 1*scale ;
	    end
	    axis ([xmin xmax ylo yhi]); axis ('off'); hold on ;
	    ppos = get(gcf,'paperposition') ;
	    pos = get(gca,'position') ;
	    hscale = (xmax - xmin) / (pos(3)*ppos(3)) ;   %units/inch
	    vscale = (yhi - ylo) / (pos(4)*ppos(4)) ;     %units/inch
%
%.......... Plot the data (trace or stick).
%
	    [m,n] = size(linetype) ; ltype = linetype(mod(ivar-1,m)+1,:) ;
	    if (~stick(ivar)),
		if (idir), ltype = [ltype(1) '. '] ; end
		plot_handle = plot (times,data,ltype,'LineWidth',linewidth(ivar),...
		                    'MarkerSize',linemark(ivar)) ;
	    else
		drawstick (times,data(:,1),data(:,2),hscale,vscale,ltype) ;
	    end
%
%.......... Count no. of plots done in this plot box.
%
	    done_plot(plot_no(ivar)) = done_plot(plot_no(ivar)) + 1 ;
	    if (~isempty(plot_handle)),
	        hand(done_plot(plot_no(ivar))) = plot_handle ;
	    end
%
%..................... If this last plot in this box, plot the box .............................
%
%
	    if (done_plot(plot_no(ivar)) == nplotperbox(plot_no(ivar))),
%
%.......... Draw the y-axes lines (& x-axis lines if clipping turned on
%.......... or its a multi-period plot).
%
	    if (yclip(plot_no(ivar)) | mult_period),
		xbox = [xmin xmax xmax xmin xmin] ;
		dyy = yaxmax-yaxmin ;
		ybox = [yaxmin yaxmin yaxmax-0.001*dyy yaxmax-0.001*dyy yaxmin] ;
		plot (xbox,ybox,'k-') ;
	    else
		plot([xmin xmin],[yaxmin yaxmax],'k-') ;
		plot([xmax xmax],[yaxmin yaxmax],'k-') ;
	    end
%
%.......... Set y-axes tick lengths & positions.
%
	    yticklen = 0.20 * hscale / 2.54 ;
	    ytrange = ytmax - ytmin ; dyt = ytrange/(ntick-1) ;
	    yticks = [ytmin:dyt:ytmax] ;
%
%.......... Tick & number the y-axes.
%
	    itick = 1 ;
	    if (labside == 'left '), ilabl = 1 ; ilabr = 0 ; end
	    if (labside == 'right'), ilabl = 0 ; ilabr = 1 ; end
	    if (~idir),
	        yleft (yaxmin,yaxmax,xmin,...
		       yticks,itick,yticklen,ilabl,font3) ;
	        yright (yaxmin,yaxmax,xmax,...
		        yticks,itick,yticklen,ilabr,font3);
	    else
	        yleft_dir (yaxmin,yaxmax,xmin,...
		           yticks,itick,yticklen,ilabl,font3) ;
	        yright_dir (yaxmin,yaxmax,xmax,...
			    yticks,itick,yticklen,ilabr,font3);
	    end
%
%.......... Tidy up variable names etc.
%
	    ind = findstr('_',longname) ; longname(ind) = blanks(length(ind)) ;
	    ind = findstr('_',unitstring);unitstring(ind) = blanks(length(ind));
%
%.......... Label the y-axes.
%
	    if (exist('y_label')),
		labname = deblank(y_label(plot_no(ivar),:)) ;
	    else
		if (~isempty(longname)),
		    labname = longname ;
		else
		    labname = ' ' ;
		end
	    end
	    if (exist('y_units')),
		unitname = deblank(y_units(plot_no(ivar),:)) ;
	    else
		if (~isempty(unitstring)),
		    unitname = ['(' unitstring ')'] ;
		else
		    unitname = ' ' ;
		end
	    end
	    yaxlab = 0.5 * (yaxmax + yaxmin);
	    fonthdat3 = font3 * hscale / 72. ;
	    if (labside == 'left '),
	        xaxlab1 = -yticklen - 1.0*(maxfigl)*fonthdat3 ;
	        xaxlab2 = xaxlab1 - 1.4*fonthdat3 ;
%	        xaxlab1 = -yticklen - 1.2*(maxfigl)*fonthdat3 ;
%	        xaxlab2 = xaxlab1 - 1.8*fonthdat3 ;
	        vlabel (xaxlab2,yaxlab,labname,1,font2) ;
	        if (~idir), vlabel (xaxlab1,yaxlab,unitname,1,font2) ; end
	    else
	    	xaxlab1 = xmax + yticklen + (maxfigr)*fonthdat3 ;
	        xaxlab2 = xaxlab1 + 1.6*fonthdat3 ;
	        vlabel (xaxlab2,yaxlab,labname,-1,font2) ;
	        if (~idir), vlabel (xaxlab1,yaxlab,unitname,-1,font2) ; end
	    end
%
%.......... Determine x-axis tick lengths & what time resolution to label.
%
	    xticklen = min([0.075*(yaxmax-yaxmin),0.20*vscale/2.54]) ;
	    tickday = 0 ; tickmon = 0 ; tickyear = 0 ;
	    ilabday = 0 ; ilabmon = 0 ; ilabyear = 0 ;
	    if (xmax<=92),
		tickday = xticklen ;
		tickmon = 2*xticklen ;
%		tickyear = yaxmax - yaxmin ;
		tickyear = 2*xticklen ;
	        if (mult_period | plot_no(ivar)==nplot),
		    ilabday = 1 ; ilabmon = 1 ; ilabyear = 1 ;
		end
	    end
	    if (xmax>92 & xmax<=1200),
		tickmon = xticklen ;
		tickyear =  2*xticklen ;
	        if (mult_period | plot_no(ivar)==nplot),
		    ilabmon=1 ; ilabyear=1 ;
		end
	    end
	    if (xmax>1200),
		tickyear = xticklen ;
	        if (mult_period | plot_no(ivar)==nplot),
		    ilabyear = 1 ;
		end
	    end
%
%.......... Tick & label the bottom & top x-axes (& draw x-axis lines if
%.......... clipping turned off).
%
	    fontvdat3 = font3 * vscale / 72. ;
	    if (yclip(plot_no(ivar)) | mult_period),
	    	ydate = yaxmin ;
		dates (0,pstart,pend,ydate,tickday,tickmon,tickyear,...
	               ilabday,ilabmon,ilabyear,font3,fontvdat3) ;
	    else
	    	if (plot_no(ivar)==nplot),
	    	    ydate = yaxmin - 0.2*(yaxmax-yaxmin) ;
	   	    plot([xmin xmax],[ydate ydate],'k-') ;
	   	    plot([xmin xmin],[ydate ydate+xticklen],'k-') ;
	    	    plot([xmax xmax],[ydate ydate+xticklen],'k-') ;
		    dates (0,pstart,pend,ydate,tickday,tickmon,tickyear,...
			   ilabday,ilabmon,ilabyear,font3,fontvdat3);
		end
		if (plot_no(ivar)==1),
	    	    ydate = yaxmax + 0.25*(yaxmax-yaxmin) ;
	   	    plot([xmin xmax],[ydate ydate],'k-') ;
	   	    plot([xmin xmin],[ydate ydate-xticklen],'k-') ;
	    	    plot([xmax xmax],[ydate ydate-xticklen],'k-') ;
		    dates (0,pstart,pend,ydate,-tickday,-tickmon,-tickyear,...
		           0,0,0,0,0);
		end
	    end
%
%.......... Plot a zero line.
%
	    if (zero_line(plot_no(ivar))),
		if (yaxmin<=0 & yaxmax>=0),
	            hz = plot ([0 xmax],[0 0],'k-') ; set(hz,'Linewidth',0.3) ;
		end
	    end
%
%.......... Plot header at top of page.
%
	    if (exist('header') & plot_no(ivar)==1),
		yheader = yaxmax + 0.02*(yaxmax-yaxmin) + 3.*font1*vscale/72. ;
%		yheader = yaxmax + 0.2*(yaxmax-yaxmin) + 3.*font1*vscale/72. ;
%		yheader = yaxmax + 0.1*(yaxmax-yaxmin) ;
	        hlabel2 (0.5*xmax,yheader,header,font1) ;
	    end
%
%.......... Put a title on top,centre of plot.
%
	    xtitle = xmin + 0.5*(xmax-xmin) ;
	    ytitle = yaxmax + 0.025*nplotperpage*(yaxmax-yaxmin) ;
	    if (exist('ptitle')),
		hlabel2 (xtitle,ytitle,ptitle(plot_no(ivar),:),font2) ;
	    else
		if (~isempty(longname)),
		    hlabel2 (xtitle,ytitle,longname,font2) ;
		end
	    end
%
%.......... Put a sub-title in top, left corner of plot.
%
	    if (exist('subtitle')),
		xsubtitle = xmin + 0.02*(xmax-xmin) ;
	        ysubtitle = yaxmax - .15*(yaxmax-yaxmin) ;
		hlabell (xsubtitle,ysubtitle,subtitle(plot_no(ivar),:),font3) ;
	    end
%
%.......... Plot a legend in top, right corner of plot.
%
	    if (exist('plegend')),
		[m,n] = size(plegend) ;
		if (m>=plot_no(ivar)),
		    eval(['legend(gca,hand,' plegend(plot_no(ivar),:) ');']) ;
	        end
	    end
%
	    end
	end
