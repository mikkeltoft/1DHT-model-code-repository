% This MATLAB script saves the figures (simulation data) created 
% by the script Code_1DHT_model.m
% It must be run after conse

fpath='C:\Documents (2)\MATLAB\ARTICLE\For repository\Results\';   % << Set absolute directory (~\) and name folder ('Name') according to scenario
for i=11:10:121
    if ismember((i-1)/10,col_incl)
    Fpath=[fpath 'Temperature profiles'];
    mkdir(Fpath);
    FigH = figure(i);
    name = ['TempProfiles' num2str(i)];
    saveas(FigH,fullfile(Fpath,name), 'fig');
    end
end
for i=12:10:122
    if ismember((i-2)/10,col_incl)
    Fpath=[fpath 'Volumetric heat capacity'];
    mkdir(Fpath);
    FigH = figure(i);
    name = ['VolHeatCap' num2str(i)];
    saveas(FigH,fullfile(Fpath,name), 'fig');
    end
end
for i=13:10:123
    if ismember((i-3)/10,col_incl)
    Fpath=[fpath 'Equivalent thermal diffusivity'];
    mkdir(Fpath);
    FigH = figure(i);
    name = ['EquivThermDiff' num2str(i)];
    saveas(FigH,fullfile(Fpath,name), 'fig');
    end
end
for i=14:10:124
    if ismember((i-4)/10,col_incl)
    Fpath=[fpath 'Ice and water fractions'];
    mkdir(Fpath);
    FigH = figure(i);
    name = ['FractionIceWater' num2str(i)];
    saveas(FigH,fullfile(Fpath,name), 'fig');
    end
end
for i=16:10:126
    if ismember((i-6)/10,col_incl)
    Fpath=[fpath 'Aggradation rate'];
    mkdir(Fpath);
    FigH = figure(i);
    name = ['AggrRate' num2str(i)];
    saveas(FigH,fullfile(Fpath,name), 'fig');
    end
end
for i=17:10:127
    if ismember((i-7)/10,col_incl)
    Fpath=[fpath 'FF and PF depths'];
    mkdir(Fpath);
    FigH = figure(i);
    name = ['FFandPFdepths' num2str(i)];
    saveas(FigH,fullfile(Fpath,name), 'fig');
    end
end
%%
Fpath=[fpath 'Temperature curve used for simulation'];
    mkdir(Fpath);
    FigH = figure(15);
    name = ['Temperature during simulation'];
    saveas(FigH,fullfile(Fpath,name), 'fig');

Fpath=[fpath 'FF and PF depths\Final FF & PF'];
    mkdir(Fpath);
    FigH = figure(18);
    name = ['FinalPF&FFdepth'];
    saveas(FigH,fullfile(Fpath,name), 'fig');

Fpath=[fpath 'Aggradation rate\Final aggrRate'];
    mkdir(Fpath);
    FigH = figure(19);
    name = ['FinalAggrRate'];
    saveas(FigH,fullfile(Fpath,name), 'fig');

close all force
toc
