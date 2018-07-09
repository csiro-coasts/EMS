clear all
close all
% this routine uses the plotts.m fucntion to plot time serie of the model
% output


% define some parameters for the plot
header = 'Estuary test case'
nplotperpage = 2 ;
inpfile = ['../out/loc1.ts']
plot_start   = [2000 01 10 01 00 00] ;
plot_end   = [2000 5 01 01 00 00] ;
plotbot    = [0.05 0.30 0.55 0.80 ] ;
plotheight = [.15 .15 0.15 0.14 ] ;
file_no   = [  1   1  1  1] ;
linetype  = ['b- ' ; 'g- ' ; 'c- ' ; 'm- '] ;
% then plot each variable

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
var_name = ['salt' 
            'temp'
            'TN  '
            'TP  ']
plotts;
print -dpsc ts.ps ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure 
var_name = ['NO3  ' 
            'NH4  '
            'DIP  '
            'Chl_a']
plotts
print -dpsc nutrient.ps ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure 
var_name = ['DetR_N ' 
            'DetPL_N'
            'DOR_N  '
            'DIN    ']
plotts
print -dpsc det.ps ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure 
var_name = ['PhyS_N' 
            'PhyL_N'
            'ZooS_N' 
            'ZooL_N']
plotts

print -dpsc phy_zoo.ps ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure 
var_name = ['PhyS_N' 
            'PhyL_N'
            'ZooS_N' 
            'ZooL_N']
plotts

print -dpsc phy_zoo.ps ;
