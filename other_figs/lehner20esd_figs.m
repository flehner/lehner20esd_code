close all
clear all

% ---------------------------
% lehner20esd_figs.m
%
% Description: script to create Figures 1-2 and 5-9 and Figures S1-S5 in Lehner et al., 2020,
% https://doi.org/10.5194/esd-11-491-2020
%
% The full model data is too large to share on github, but is available via:
% - Large Ensembles: http://www.cesm.ucar.edu/projects/community-projects/MMLEA/
% - CMIP5: https://esgf-node.llnl.gov/search/cmip5/
% - CMIP6: https://esgf-node.llnl.gov/search/cmip6/
% Note: the models were regridded to a common 2.5x2.5° grid (indicated by the "g025" in the file names)
%
% Global observations are provided in the repo:
%   Mildly pre-processed observational data; for original data check
%   the institutional sources:
%   - Land_and_Ocean_LatLong1.185001-201812_ts.nc (Berkeley Earth: http://berkeleyearth.org/archive/land-and-ocean-data/)
%   - precip.mon.mean.197901-201812_ts.nc (GPCP: http://gpcp.umd.edu/)
%
% Special functions and landmasks are also provided in the repo.
%
% Author: Flavio Lehner, May 2020, flehner@ucar.edu or flavio.lehner@cornell.edu
%
% ------------------------------------------------------------------------------

% -- which plot?
plot1   = 1; % [Fig. 1+2]
plot5   = 0; % [Fig. 5]
plot7   = 0; % [Fig. 7]
plot8   = 0; % [Fig. 8]
plot9   = 0; % [Fig. 9]

plotS1  = 0; % [Fig. S1]
plotS2  = 0; % [Fig. S2]
plotS3  = 0; % [Fig. S3]
plotS4  = 0; % [Fig. S4]
plotS5  = 0; % [Fig. S5]

% -- which variable?
v = 1; % 1=tas, 2=pr

% -- which season?
s = 3; % 1=DJF, 2=JJA, 3=annual

% -- which region?
r = 1; % e.g., 1=global

% ------------------------------------------------------------------------------

% -- region to plot --
regions     = {'global','uk','sahel','europe','southern_ocean',...
               'se_asia','us_southwest','nino34','north_america','kolkata',...
               'dallas','alaska','seattle','upper_colorado','sydney',...
               'NH','SH','arctic','AA','zurich',...
               'southern_europe','northern_europe','denver','anchorage','india',...
               'sahara'};

load_control = 0; % 0=no, 1=yes
if plotS1 == 1
  load_control = 1; % 0=no, 1=yes
end

pathin      = '~/Dropbox/work/';
pathout_fig = '~/Dropbox/publication/lehner19_revisiting_hawkins_sutton/fig/hawkins_plots/'
vars        = {'tas','pr','psl'};
comp        = 'Amon';
seasons     = {'DJF','JJA','annual','JJAS'};

% -- CLIVAR parameters
models      = {'cesm_lens','canesm2_lens','csiro_mk36_lens','gfdl_cm3_lens','gfdl_esm2m_lens','ec_earth_lens','mpi_lens'};
model_names = {'CESM1-CAM5','CanESM2','CSIRO-Mk3-6-0','GFDL-CM3','GFDL-ESM2M','EC-EARTH','MPI-ESM'};
ensmem0     = [35,50,30,20,30,16,100];
start0      = [1920,1950,1850,1920,1950,1860,1850];
ende0       = [2100,2100,2100,2100,2100,2100,2099];

% -- CMIP5 parameters
scen_cmip5 = {'rcp26','rcp45','rcp85'};
scen_cmip6 = {'ssp126','ssp245','ssp370','ssp585'};
% -- selection for scenario uncertainy w/ full CMIP5
models_cmip5 = {'bcc-csm1-1-m','bcc-csm1-1','BNU-ESM',...
                'CanESM2','CCSM4','CESM1-CAM5',...
                'CNRM-CM5','CSIRO-Mk3-6-0','EC-EARTH','FGOALS-g2',...
                'FIO-ESM','GFDL-CM3','GFDL-ESM2G','GFDL-ESM2M',...
                'GISS-E2-H','GISS-E2-R','HadGEM2-AO',...
                'HadGEM2-ES','IPSL-CM5A-LR','IPSL-CM5A-MR',...
                'MIROC5','MIROC-ESM-CHEM','MIROC-ESM','MPI-ESM-LR','MPI-ESM-MR',...
                'MRI-CGCM3','NorESM1-ME','NorESM1-M'};
% -- how to exclude certain models from above list:
% models_cmip5(ismember(models_cmip5,'FIO-ESM')) = [];
% models_cmip5(ismember(models_cmip5,'GFDL-CM3')) = [];
ensmem_cmip5  = ones(length(models_cmip5),1);
% -- selection for scenario uncertainy w/ only 7 LEs
models_le_cmip5  = {'CESM1-CAM5','CanESM2','CSIRO-Mk3-6-0','GFDL-CM3','GFDL-ESM2M','EC-EARTH','MPI-ESM-LR'};
models_le_cmip5_id = find(contains(models_cmip5,models_le_cmip5));
% -- CMIP6 parameters
constrain = 0; % weight models according to trend performance

% ------------------------------------------------------------------------------
% -- selection used for Lehner et al. 2020, Earth System Dynamics (fewer models
%    than full CMIP6, since compiled in November 2019) --------------
models_cmip6_paper = ...
               {'BCC-CSM2-MR','CAMS-CSM1-0','CESM2','CESM2-WACCM','CNRM-CM6-1',...
                'CNRM-ESM2-1','CanESM5','EC-Earth3','EC-Earth3-Veg','FGOALS-f3-L',...
                'FGOALS-g3','GFDL-ESM4','INM-CM4-8','INM-CM5-0','IPSL-CM6A-LR',...
                'MCM-UA-1-0','MIROC-ES2L','MIROC6','MPI-ESM1-2-HR','MRI-ESM2-0',...
                'UKESM1-0-LL'};
models_cmip6 = models_cmip6_paper;
% ------------------------------------------------------------------------------
% % -- activate this section if you want more CMIP6 models than are in the paper
% command = 'ls -1 /Users/flehner/Dropbox/work/cmip6-ng/tas/tas*ssp*nc | awk -F''[_]'' ''{print $3}'' | sort | uniq > /Users/flehner/Dropbox/work/cmip6-ng/tas/tmp.txt';
% system(command);
% models_cmip6 = table2cell(readtable('/Users/flehner/Dropbox/work/cmip6-ng/tas/tmp.txt'))'
models_cmip6 = {'ACCESS-CM2','ACCESS-ESM1-5','AWI-CM-1-1-MR','BCC-CSM2-MR','CAMS-CSM1-0',...
                'CESM2','CESM2-WACCM','CNRM-CM6-1','CNRM-CM6-1-HR','CNRM-ESM2-1',...
                'CanESM5','CanESM5-CanOE','EC-Earth3','EC-Earth3-Veg','FGOALS-f3-L',...
                'FGOALS-g3','FIO-ESM-2-0','GFDL-CM4','GFDL-ESM4','GISS-E2-1-G',...
                'HadGEM3-GC31-LL','INM-CM4-8','INM-CM5-0','IPSL-CM6A-LR','KACE-1-0-G',...
                'MCM-UA-1-0','MIROC-ES2L','MIROC6','MPI-ESM1-2-HR','MPI-ESM1-2-LR',...
                'MRI-ESM2-0','NESM3','NorESM2-LM','NorESM2-MM','UKESM1-0-LL'};
% -- how to exclude certain models from above list:
models_cmip6(ismember(models_cmip6,'GFDL-CM4')) = [];
models_cmip6(ismember(models_cmip6,'NorESM2-LM')) = []; % pr: ssp585 r1i1p1f1 exists, but the historical r1i1p1f1 only starts in 1950 (?)
models_cmip6(ismember(models_cmip6,'FIO-ESM-2-0')) = []; % no ssp370
models_cmip6(ismember(models_cmip6,'HadGEM3-GC31-LL')) = []; % no ssp370
models_cmip6(ismember(models_cmip6,'NESM3')) = []; % no ssp370
% ------------------------------------------------------------------------------
ensmem_cmip6  = ones(length(models_cmip6),1);
% ------------------------------------------------------------------------------

start_cmip5  = 1870;
ende_cmip5   = 2100;
start_cmip6  = 1850;
ende_cmip6   = 2100;
% -- general paramaters
start       = 1950; % 1950
ende        = 2099;
time        = start:ende;
refstart    = 1995; % 1986 1995 1971 1950 2001 1979 1991
refende     = 2014; % 2005 2014 2000 1979 2030 2000 2020
reflength   = refende-refstart+1;
lati        = [-88.75:2.5:88.75];
loni        = [1.25:2.5:358.75];
wl          = 10; % 1 10 % runnig mean length in years
% -- area --
filein  = [pathin 'cmip5-ng/area_g025.nc'];
area    = ncread([filein],'AREA');
filein  = [pathin 'cmip5-ng/landfrac_g025.nc'];
land    = ncread([filein],'LANDFRAC');
id = isnan(land);
land(id) = 0;
landarea = land.*area;


% -- load CLIVAR data
region = regions{r}
clear('lon2','var0','seas0')
if strcmp(region,'upper_colorado')==1
  lat         = [38.75 41.25];
  lon         = [248.75 251.25 253.75];
    seas0     = 3;
    var0      = 2;
end
if strcmp(region,'us_southwest')==1
  lat         = [32.5 42.1];
  lon         = [-125 -103]+360;
end
if strcmp(region,'europe')==1
  lat         = [36.4 70.1];
  lon         = [0 23.4];
  lon2        = [350 360];
end
if strcmp(region,'southern_europe')==1
  lat         = [36.4 50.1];
  lon         = [0 23.4];
  lon2        = [350 360];
    % seas0     = 2;
    % var0      = 1;
end
if strcmp(region,'northern_europe')==1
  lat         = [50.1 70.1];
  lon         = [0 23.4];
  lon2        = [350 360];
end
if strcmp(region,'se_asia')==1
  lat         = [7.9 30.1];
  lon         = [93.2 122.1];
end
if strcmp(region,'india')==1
  lat         = [6.2 30.1];
  lon         = [69.9 87.3];
    seas0     = 2;
    var0      = 2;
end
if strcmp(region,'nino34')==1
  lat         = [-5 5];
  lon         = [-170 -120]+360;
  % lon         = [-190 -140]+360;
    seas0        = 1;
    var0        = 1;
end
if strcmp(region,'uk')==1
  lat         = [50.2:2.5:58.8];
  lon         = [-10.9:2.5:0]+360;
    seas0       = 3;
    var0        = 1;
end
if strcmp(region,'sydney')==1
  lat         = -33.8;
  lon         = 150.0;
end
if strcmp(region,'dallas')==1
  lat         = 32.7;
  lon         = -96.8+360;
end
if strcmp(region,'seattle')==1
  lat         = 49.2;
  lon         = -121.4+360;
    var0      = 2;
    seas0      = 1;
end
if strcmp(region,'sahara')==1
  lat         = 22.5;
  lon         = 9.6;
    var0      = 2;
    seas0     = 3;
end
if strcmp(region,'denver')==1
  lat         = 39.9;
  lon         = -100.1+360;
end
if strcmp(region,'anchorage')==1
  lat         = 61.2;
  lon         = -149.1+360;
end
if strcmp(region,'north_america')==1
  lat         = [22.1 70.5];
  lon         = [-124.4+360 -53.7+360];
end
if strcmp(region,'kolkata')==1
  lat         = 22.6;
  lon         = 88.4;
end
if strcmp(region,'southern_ocean')==1
  lat         = [-70 -60];
  lon         = [0 360];
    seas0      = 3;
    var0      = 1;
end
if strcmp(region,'sahel')==1
  lat         = [10 20];
  lon         = [0 40];
    var0      = 2;
    seas0     = 2;
end
if strcmp(region,'NH')==1
  lat         = [0 90];
  lon         = [0 360];
end
if strcmp(region,'SH')==1
  lat         = [-90 0];
  lon         = [0 360];
end
if strcmp(region,'alaska')==1
  lat         = [51.6 73.9];
  lon         = [191.2 255.9];
end
if strcmp(region,'zurich')==1
  lat         = 47.5;
  lon         = 8.5;
end
if strcmp(region,'arctic')==1
  lat         = [75 90];
  lon         = [0 360];
end
if strcmp(region,'AA')==1
  lat1        = [0 90];
  lat2        = [75 90];
  lon1        = [0 360];
  lon2        = [0 360];
end
% --
% ncells      = (length(lat)*length(lon));
% -- find grid cells to plot in models
clear('tmp','tmp0')
if strcmp(region,'global')==1
  ii      = 1:length(lati);
  jj      = 1:length(loni);
  weights = area(jj,ii);
  weights = weights/sum(sum(weights));
else if strcmp(region,'AA')==1
  for x = 1:length(lat1)
    tmp0    = find(abs(lati-lat1(x))==min(abs(lati-lat1(x))));
    tmp1(x)  = tmp0(1);
  end
  for x = 1:length(lat2)
    tmp0    = find(abs(lati-lat2(x))==min(abs(lati-lat2(x))));
    tmp2(x)  = tmp0(1);
  end
  ii1 = tmp1(1):tmp1(end);
  ii2 = tmp2(1):tmp2(end);
  for x = 1:length(lon1)
    tmp0    = find(abs(loni-lon1(x))==min(abs(loni-lon1(x))));
    tmp(x)  = tmp0(1);
  end
  jj1 = tmp(1):tmp(end);
  jj2 = jj1;
  % -- cut out area --
  weights1 = area(jj1,ii1);
  weights2 = area(jj2,ii2);
  weights1 = weights1/sum(sum(weights1));
  weights2 = weights2/sum(sum(weights2));
  clear('tmp1','tmp2')
else
  for x = 1:length(lat)
    tmp0    = find(abs(lati-lat(x))==min(abs(lati-lat(x))));
    tmp(x)  = tmp0(1);
  end
  ii = tmp(1):tmp(end);
  for x = 1:length(lon)
    tmp0    = find(abs(loni-lon(x))==min(abs(loni-lon(x))));
    tmp(x)  = tmp0(1);
  end
  jj = tmp(1):tmp(end);
  if exist('lon2')==1
    clear('tmp')
    for x = 1:length(lon2)
      tmp0    = find(abs(loni-lon2(x))==min(abs(loni-lon2(x))));
      tmp(x)  = tmp0(1);
    end
    jj = [jj tmp(1):tmp(end)];
  end
  % -- cut out area --
  if strcmp(region,'nino34')==1 ||...
     strcmp(region,'southern_ocean')==1 ||...
     strcmp(region,'arctic')==1
    weights = area(jj,ii);
    weights = weights/sum(sum(weights));
  else
    weights = landarea(jj,ii);
    weights = weights/sum(sum(weights));
  end
end
end
clear('tmp','tmp0')

% -- variable --
if r == 1
  vari = vars{v}
else
  if exist('var0')==0
    var0 = v;
  end
  vari = vars{var0}
end
if strcmp(vari,'tas')==1
  units       = 'K';
  f           = 1;
  var_name    = 'Temperature'
  var_name2   = 'temperature';
elseif strcmp(vari,'pr')==1
  units       = '%';
  f           = 86400;
  var_name    = 'Precipitation'
  var_name2   = 'precipitation';
elseif strcmp(vari,'psl')==1
  units       = 'unitless';
  f           = 1;
  var_name    = 'Sea level pressure'
end
% -- load global mean observations
start_trend = 1981;
ende_trend  = 2014;
time_trend  = start_trend:ende_trend;
if strcmp(vari,'tas')==1 && strcmp(region,'global')==1
  obs_name = 'BEST';
  time_obs = 1850:2018;
  obs_raw = ncread([pathin 'observations/temperature/best/Land_and_Ocean_LatLong1.185001-201812_ts.nc'],'TEMPERATURE');
  obs = rm(am(obs_raw),wl);
  obs = obs-nanmean(obs(refstart-1850+1:refende-1850+1));
  if wl == 1
    tmp = polyfit(time_trend',obs(start_trend-1850+1:ende_trend-1850+1),1);
  else
    tmp = polyfit(time_trend,obs(start_trend-1850+1:ende_trend-1850+1),1);
  end
  obs_trend   = tmp(1);
else
  obs_name = 'GPCP';
  time_obs = 1986:2018;
  obs_raw  = ncread([pathin 'observations/precipitation/gpcp/precip.mon.mean.197901-201812_ts.nc'],'PRECIP');
  obs_raw  = obs_raw(85:end); % limit to 1986 onwards, based on Angie's assessment
  obs  = rm(am(obs_raw),wl);
  ref  = nanmean(obs(refstart-1986+1:refende-1986+1));
  obs  = ((obs-ref)/ref)*100;
end
tmp1                = rm(am(obs_raw),wl)'; % do smoothing
if wl > 1
  tmp1 = tmp1';
end
idx = ~isnan(tmp1);
fit                 = rm(polyval(polyfit(time_obs(idx),tmp1(idx),4),time_obs),10); % do polyfit and do smoothing
idx = ~isnan(fit);
obs_ts_em_hs        = fit-nanmean(fit(refstart-time_obs(1)+1:refende-time_obs(1)+1));  % do anomalies
obs_residual        = am(obs_raw)-fit;


% -- load GMST (need always)
time_tas_obs = 1850:2018;
tas_obs_raw = ncread([pathin 'observations/temperature/best/Land_and_Ocean_LatLong1.185001-201812_ts.nc'],'TEMPERATURE');
tas_obs = rm(am(tas_obs_raw),wl);
ref_pi = nanmean(tas_obs(1850-1850+1:1899-1850+1)); % preindustrial reference
ref = nanmean(tas_obs(refstart-1850+1:refende-1850+1));
tas_obs_wtd = ref-ref_pi; % wtd = warming-to-date
tas_obs = tas_obs-nanmean(tas_obs(refstart-1850+1:refende-1850+1));


% -- season --
if r == 1
  seas = seasons{s}
else
  if exist('seas0')==0
    seas0 = s;
  end
  seas  = seasons{seas0}
end



% -- load CMIP5 data --
start_cmip5 = 1870;
ende_cmip5  = 2100;
cmip5_ts_em  = NaN(length(models_cmip5),length(time));
scen        = scen_cmip5;
pathin_tmp  = [pathin 'cmip5-ng/' vari '/'];
for sc = 1:3 % scenarios
  eval(['cmip5_' scen{sc} '_ts_residual = cmip5_ts_em;'])
  clear('ne')
  for m = 1:length(models_cmip5)
    ['scenario = ' scen{sc} ' / model = ' models_cmip5{m}]
    files  = dir([pathin_tmp vari '_mon_' models_cmip5{m} '_' scen{sc} '_r*i1p1_g025.nc']);
    if sc == 3
      ne(m) = length(files);
    else
      ne(m) = 1;
    end
    for e = 1:ne(m) %1:length(ensmem_cmip6(m))
      % if strcmp(models_cmip5{m},'EC-EARTH')==1
      %   filein  = [pathin 'cmip5-ng/' vari '/' vari '_mon_' models_cmip5{m} '_' scen{sc} '_r8i1p1_g025.nc'];
      % else
      %   filein  = [pathin 'cmip5-ng/' vari '/' vari '_mon_' models_cmip5{m} '_' scen{sc} '_r' num2str(e) 'i1p1_g025.nc'];
      % end
      tmp0    = ncread([pathin_tmp files(e).name],[vari])*f;
      for i = 1:length(tmp0)
        if strcmp(region,'AA')==1
          tmp1(i)       = nansum(nansum(squeeze(tmp0(jj1,ii1,i)).*weights1)); % do area-weighted mean
          tmp2(i)       = nansum(nansum(squeeze(tmp0(jj2,ii2,i)).*weights2)); % do area-weighted mean
          tmp(i)     = tmp1(i)-tmp2(i);
        else
          tmp(i)     = sum(sum(squeeze(tmp0(jj,ii,i)).*weights));
        end
      end
      tmp     = tmp((start-start_cmip5)*12+1:(ende-start_cmip5+1)*12);
      if strcmp(seas,'annual')==1
        tmp    = am(tmp);
      else
        tmp    = seasmean(tmp,seas);
      end
      eval(['cmip5_' scen{sc} '_ts_raw(m,e,:) = tmp;'])
      ref = nanmean(tmp(refstart-start+1:refende-start+1));
      if strcmp(vari,'pr')==1
        eval(['cmip5_' scen{sc} '_ts(m,:)   = ((tmp-ref)/ref)*100;'])
      else
        eval(['cmip5_' scen{sc} '_ts(m,:)   = tmp-ref;'])
      end
      eval(['cmip5_' scen{sc} '_ts_all(m,e,:) = cmip5_' scen{sc} '_ts(m,:);'])
      eval(['cmip5_' scen{sc} '_ts_all(cmip5_' scen{sc} '_ts_all==0) = NaN;'])
    end
      % % -- ensemble mean for each model
      % eval(['cmip5_' scen{sc} '_ts_em(m,:)  = nanmean(cmip5_' scen{sc} '_ts(m,:,:),2);'])
    % -- H&S style calculation --
    eval(['tmp1 = squeeze(cmip5_' scen{sc} '_ts_raw(m,1,:));'])
    idx = ~isnan(tmp1);
    fit = NaN(size(tmp1));
    % fit(idx) = rm(polyval(polyfit(time(idx),tmp1(idx)',4),time(idx)),10); % polynomial fit and then 10-yr running mean
    fit(idx) = polyval(polyfit(time(idx),tmp1(idx)',4),time(idx)); % polynomial fit
    idx = ~isnan(fit);
    ref = nanmean(fit(refstart-start+1:refende-start+1));
    if strcmp(vari,'pr')==1
      eval(['cmip5_' scen{sc} '_ts_em(m,:)   = ((fit-ref)/ref)*100;'])
    else
      eval(['cmip5_' scen{sc} '_ts_em(m,:)   = fit-ref;'])
    end
    tmp1 = tmp1';
    residual  = tmp1-fit' + nanmean(fit);
    if strcmp(vari,'pr')==1
      ref       = nanmean(residual);
      residual  = (residual/ref)*100;
    end
    tmp2 = rm(residual,wl);
    eval(['cmip5_' scen{sc} '_ts_residual(m,:)  = tmp2;'])
    eval(['cmip5_' scen{sc} '_ts_var1(m)        = nanvar(tmp2);'])
    eval(['cmip5_' scen{sc} '_ts_var2(m)        = nanvar(tmp2(1:refende-start+1));'])
    eval(['cmip5_' scen{sc} '_ts_var3(m)        = nanvar(tmp2(refende-start+1+1:end));'])
  end
end
% -- calculates weights based on trend performance:
if strcmp(region,'global')==1 & strcmp(seas,'annual')==1 & strcmp(vari,'tas')==1
  clear('w')
  for m = 1:length(models_cmip5)
    for e = 1:ne(m)
      tmp = polyfit(time_trend',squeeze(cmip5_rcp85_ts_raw(m,e,start_trend-start+1:ende_trend-start+1)),1);
      trends1_cmip5(m,e) = tmp(1);
      tmp = polyfit(1:ende-start_trend+1,squeeze(cmip5_rcp85_ts_raw(m,e,start_trend-start+1:end))',1);
      trends2_cmip5(m,e) = tmp(1);
      w(m,e) = 1/(obs_trend+abs(trends1_cmip5(m,e)-obs_trend));
    end
  end
  weights_cmip5 = (w(:,1)./sum(w(:,1)))./(1/length(models_cmip5));
  for i = 1:100
    clear('wi')
    for m = 1:length(models_cmip5)
      ni = randi(ne(m));
      wi(m) = w(m,ni);
    end
    weights_cmip5_bs(:,i) = (wi./sum(wi))./(1/length(models_cmip5));
  end
  % -- write trends and weights into txt file:
  fid = fopen([pathin 'time_series/trends1_cmip5.txt'], 'w');
  fprintf(fid, '%f \n', trends1_cmip5);
  fclose(fid);
  fid = fopen([pathin 'time_series/trends2_cmip5.txt'], 'w');
  fprintf(fid, '%f \n', trends2_cmip5);
  fclose(fid);
  fid = fopen([pathin 'time_series/weights_cmip5.txt'], 'w');
  fprintf(fid, '%f \n', weights_cmip5);
  fclose(fid);
  fid = fopen([pathin 'time_series/weights_cmip5_bs.txt'], 'w');
  fprintf(fid, '%f \n', weights_cmip5_bs);
  fclose(fid);
  weights_cmip5     = (repmat(weights_cmip5',length(time),1))';
  weights_cmip5_bs  = repmat(weights_cmip5_bs',[1,1,length(time)]);
else
  weights_cmip5     = ones(length(models_cmip5),length(time));
  weights_cmip5_bs  = ones(100,length(models_cmip5),length(time));
  % -- GMST trends
  trends1_cmip5     = load([pathin 'time_series/trends1_cmip5.txt'])';
  trends2_cmip5     = load([pathin 'time_series/trends2_cmip5.txt'])';
  % -- specific vari trends
  for m = 1:length(models_cmip5)
    tmp = polyfit(time_trend',squeeze(cmip5_rcp85_ts_raw(m,1,start_trend-start+1:ende_trend-start+1)),1);
    trends1_cmip5_vari(m) = tmp(1);
    tmp = polyfit(1:ende-start_trend+1,squeeze(cmip5_rcp85_ts_raw(m,1,start_trend-start+1:end))',1);
    trends2_cmip5_vari(m) = tmp(1);
  end
end





% -- SMILE models
le_ts        = NaN(length(models),max(ensmem0),length(time));
le_raw_ts    = le_ts;
le_ts_em_hs  = NaN(length(models),max(ensmem0),length(time));
for m = 1:length(models)
  ['SMILE model = ' models{m} ]
  filein  = [pathin models{m} '/' comp '/' vari '/' vari '_' comp '_' model_names{m} '_' seas '_g025.nc'];
  tmp0    = ncread([filein],'data')*f;
  for e = 1:ensmem0(m)
    for i = 1:length(tmp0)
      if strcmp(region,'AA')==1
        tmp1       = nansum(nansum(squeeze(tmp0(jj1,ii1,e,i)).*weights1)); % do area-weighted mean
        tmp2       = nansum(nansum(squeeze(tmp0(jj2,ii2,e,i)).*weights2)); % do area-weighted mean
        tmp(i)     = tmp1-tmp2;
      else
        tmp(i)     = nansum(nansum(squeeze(tmp0(jj,ii,e,i)).*weights)); % do area-weighted mean
      end
    end
    tmp              = tmp(start-start0(m)+1:ende-start0(m)+1); % cut to common length
    if strcmp(vari,'pr')==1
      ref = nanmean(tmp);
      le_raw_ts(m,e,:) = (tmp/ref)*100; % raw at seas resolution
      ref              = nanmean(tmp(refstart-start+1:refende-start+1));
      le_ts_anom(m,e,:) = ((tmp-ref)/ref)*100; % raw at seas resolution
    else
      le_raw_ts(m,e,:) = tmp; % raw at seas resolution
      ref              = nanmean(tmp(refstart-start+1:refende-start+1));
      le_ts_anom(m,e,:) = tmp-ref; % raw at seas resolution
    end
    tmp              = rm(tmp,wl); % do smoothing %
    ref              = nanmean(tmp(refstart-start+1:refende-start+1));
    le_ts(m,e,:)     = tmp-ref; % do anomalies
    if strcmp(vari,'pr')==1
      le_ts(m,e,:)     = ((tmp-ref)/ref)*100; % do anomalies
    end
    % -- H&S style calculation
    tmp1                 = squeeze(le_raw_ts(m,e,:)); % pick one member
    tmp1                 = rm(tmp1,wl); % do smoothing
    idx = ~isnan(tmp1);
    fit                  = rm(polyval(polyfit(time(idx),tmp1(idx),4),time),10); % do polyfit and do smoothing
    idx = ~isnan(fit);
    le_ts_em_hs(m,e,:)        = fit-nanmean(fit(refstart-start+1:refende-start+1));  % do anomalies
    le_ts_residual_hs(m,e,:)  = tmp1-fit; % do residual
    % -- ensure there's no unwanted zeros:
    le_ts_residual_hs(le_ts_residual_hs==0) = NaN;
    le_ts_var_hs1(m,e)        = nanvar(squeeze(le_ts_residual_hs(m,e,idx))); % do variance
    le_ts_var_hs2(m,e)        = nanvar(squeeze(le_ts_residual_hs(m,e,1:refende-start))); % do variance
    le_ts_var_hs3(m,e)        = nanvar(squeeze(le_ts_residual_hs(m,e,refende-start+1:end))); % do variance
  end
  % -- LE-style forced response calculation
  % le_ts_em(m,:)       = rm(nanmean(le_ts(m,:,:),2),10); % do EM with additional running mean as in HS09 (not really necessary)
  le_ts_em(m,:)       = nanmean(le_ts(m,:,:),2); % do EM
  % -- LE-style int var uncertainty calculation: variance across ensemble members
  le_ts_tmp(m,:,:)    = rm(squeeze(le_raw_ts(m,:,:)),wl,2); % do smoothing
  le_ts_var(m,:)      = nanvar(squeeze(le_ts_tmp(m,:,:)),[],1); % do smoothing and do variance
  % -- LE style int var uncertainty calculation: variance across time but after removing the "good" forced response (=ensemble mean)
  % -- residual (already-smoothed) after removing ensemble mean:
  le_ts_residual(m,:,:) = squeeze(le_ts(m,:,:)) - repmat(le_ts_em(m,:),[100,1]);
  le_ts_var1(m,:) = le_ts_em(m,:)*0+nanmean(nanvar(le_ts_residual(m,:,:),[],3));
  le_ts_var2(m,:) = le_ts_em(m,:)*0+nanmean(nanvar(le_ts_residual(m,:,1:refende-start+1),[],3));
  le_ts_var3(m,:) = le_ts_em(m,:)*0+nanmean(nanvar(le_ts_residual(m,:,refende-start+1+1:end),[],3));
end
% -- ensure no unwanted zeros:
le_ts_var_hs1(le_ts_var_hs1==0) = NaN;
le_ts_var_hs2(le_ts_var_hs2==0) = NaN;
le_ts_var_hs3(le_ts_var_hs3==0) = NaN;
le_ts_var1(le_ts_var1==0) = NaN;
le_ts_var2(le_ts_var2==0) = NaN;
le_ts_var3(le_ts_var3==0) = NaN;
clear('tmp0','tmp')


% -- load CMIP6 data --
cmip6_ts_em  = NaN(length(models_cmip6),length(time));
scen        = scen_cmip6;
pathin_tmp  = [pathin 'cmip6-ng/' vari '/'];
for sc = 1:length(scen)
  fid = fopen(['/Users/flehner/Dropbox/tmp.txt'], 'w');
  clear('ne')
  for m = 1:length(models_cmip6)
    ['scenario = ' scen{sc} ' / model = ' models_cmip6{m} ]
    files   = dir([pathin 'cmip6-ng/' vari '/' vari '_mon_' models_cmip6{m} '_' scen{sc} '_r*_g025.nc']);
    if sc == 4
      ne(m) = length(files);
    else
      ne(m) = 1;
    end
    for e = 1:ne(m) %1:length(ensmem_cmip6(m))
      % ['scenario = ' scen{sc} ' / model = ' models_cmip6{m} ' / ensmem = ' num2str(e)]
      % if strcmp(models_cmip6{m},'CAMS-CSM1-0')==1
      %   files   = dir([pathin 'cmip6-ng/' vari '/' vari '_mon_' models_cmip6{m} '_' scen{sc} '_r2i1p1f1_g025.nc']);
      % end
      tmp0    = ncread([pathin_tmp files(e).name],[vari])*f; % yyy
      fprintf(fid, '%s\n', files(e).name);
      % filein  = [pathin 'cmip6-ng/' vari '/' vari '_mon_' models_cmip6{m} '_' scen{sc} '_r' num2str(e) 'i1p1_g025.nc'];
      % tmp0    = ncread([filein],[vari])*f;
      for i = 1:length(tmp0)
        if strcmp(region,'AA')==1
          tmp1(i)       = nansum(nansum(squeeze(tmp0(jj1,ii1,i)).*weights1)); % do area-weighted mean
          tmp2(i)       = nansum(nansum(squeeze(tmp0(jj2,ii2,i)).*weights2)); % do area-weighted mean
          tmp(i)     = tmp1(i)-tmp2(i);
        else
          tmp(i)     = sum(sum(squeeze(tmp0(jj,ii,i)).*weights));
        end
      end
      tmp     = tmp((start-start_cmip6)*12+1:(ende-start_cmip6+1)*12);
      if strcmp(seas,'annual')==1
        tmp    = am(tmp);
      else
        tmp    = seasmean(tmp,seas);
      end
      eval(['cmip6_' scen{sc} '_ts_raw(m,e,:) = tmp;'])
      ref = nanmean(tmp(refstart-start+1:refende-start+1));
      if strcmp(vari,'pr')==1
        eval(['cmip6_' scen{sc} '_ts(m,:)   = ((tmp-ref)/ref)*100;'])
      else
        eval(['cmip6_' scen{sc} '_ts(m,:)   = tmp-ref;'])
      end
    end
    % % -- ensemble mean for each model
    % eval(['cmip6_' scen{sc} '_ts_em(m,:)  = nanmean(cmip6_' scen{sc} '_ts(m,:,:),2);'])
    % -- H&S style calculation --
    eval(['tmp1 = squeeze(cmip6_' scen{sc} '_ts_raw(m,1,:));'])
    idx = ~isnan(tmp1);
    fit = NaN(size(tmp1));
    % fit(idx) = rm(polyval(polyfit(time(idx),tmp1(idx)',4),time(idx)),10); % polynomial fit and then 10-yr running mean
    fit(idx) = polyval(polyfit(time(idx),tmp1(idx)',4),time(idx)); % polynomial fit
    idx = ~isnan(fit);
    ref = nanmean(fit(refstart-start+1:refende-start+1));
    if strcmp(vari,'pr')==1
      eval(['cmip6_' scen{sc} '_ts_em(m,:)   = ((fit-ref)/ref)*100;'])
    else
      eval(['cmip6_' scen{sc} '_ts_em(m,:)   = fit-ref;'])
    end
    tmp1 = tmp1';
    residual  = tmp1-fit' + nanmean(fit);
    if strcmp(vari,'pr')==1
      ref       = nanmean(residual);
      residual  = (residual/ref)*100;
    end
    tmp2 = rm(residual,wl);
    eval(['cmip6_' scen{sc} '_ts_residual(m,:)  = tmp2;'])
    eval(['cmip6_' scen{sc} '_ts_var1(m)        = nanvar(tmp2);'])
    eval(['cmip6_' scen{sc} '_ts_var2(m)        = nanvar(tmp2(1:refende-start+1));'])
    eval(['cmip6_' scen{sc} '_ts_var3(m)        = nanvar(tmp2(refende-start+1+1:end));'])
  end
  fclose(fid);
end
% -- calculates weights based on trend performance:
if strcmp(region,'global')==1 & strcmp(seas,'annual')==1 & strcmp(vari,'tas')==1
  clear('w')
  for m = 1:length(models_cmip6)
    for e = 1:ne(m)
      tmp = polyfit(time_trend',squeeze(cmip6_ssp585_ts_raw(m,e,start_trend-start+1:ende_trend-start+1)),1);
      trends1_cmip6(m,e) = tmp(1);
      tmp = polyfit(1:ende-start_trend+1,squeeze(cmip6_ssp585_ts_raw(m,e,start_trend-start+1:end))',1);
      trends2_cmip6(m,e) = tmp(1);
      w(m,e) = 1/(obs_trend+abs(trends1_cmip6(m,e)-obs_trend));
    end
  end
  weights_cmip6 = (w(:,1)./sum(w(:,1)))./(1/length(models_cmip6));
  for i = 1:100
    clear('wi')
    for m = 1:length(models_cmip6)
      ni = randi(ne(m));
      wi(m) = w(m,ni);
    end
    weights_cmip6_bs(:,i) = (wi./sum(wi))./(1/length(models_cmip6));
  end
  % -- write weights into txt file:
  fid = fopen([pathin 'time_series/trends1_cmip6.txt'], 'w');
  fprintf(fid, '%f \n', trends1_cmip6);
  fclose(fid);
  fid = fopen([pathin 'time_series/trends2_cmip6.txt'], 'w');
  fprintf(fid, '%f \n', trends2_cmip6);
  fclose(fid);
  fid = fopen([pathin 'time_series/weights_cmip6.txt'], 'w');
  fprintf(fid, '%f \n', weights_cmip6);
  fclose(fid);
  fid = fopen([pathin 'time_series/weights_cmip6_bs.txt'], 'w');
  fprintf(fid, '%f \n', weights_cmip6_bs);
  fclose(fid);
  weights_cmip6     = (repmat(weights_cmip6',length(time),1))';
  weights_cmip6_bs  = repmat(weights_cmip6_bs',[1,1,length(time)]);
else
  weights_cmip6     = ones(length(models_cmip6),length(time));
  weights_cmip6_bs  = ones(100,length(models_cmip6),length(time));
  % -- GMST trends
  % trends1_cmip6     = load([pathin 'time_series/trends1_cmip6.txt'])';
  % trends2_cmip6     = load([pathin 'time_series/trends2_cmip6.txt'])';
  % weights_cmip6     = load([pathin 'time_series/weights_cmip6.txt'])';
  % weights_cmip6_bs  = load([pathin 'time_series/weights_cmip6_bs.txt'])';
  % weights_cmip6     = (repmat(weights_cmip6,length(time),1))';
  % weights_cmip6_bs  = reshape(weights_cmip6_bs,[length(models_cmip6),100]);
  % weights_cmip6_bs  = repmat(weights_cmip6_bs',[1,1,length(time)]);
  % -- specific vari trends
  for m = 1:length(models_cmip6)
    tmp = polyfit(time_trend',squeeze(cmip6_ssp585_ts_raw(m,1,start_trend-start+1:ende_trend-start+1)),1);
    trends1_cmip6_vari(m) = tmp(1);
    tmp = polyfit(1:ende-start_trend+1,squeeze(cmip6_ssp585_ts_raw(m,1,start_trend-start+1:end))',1);
    trends2_cmip6_vari(m) = tmp(1);
  end
end




if load_control == 1
  clear('tmp','tmp1','tmp2')
  % -- load CMIP5 piControl data --
  cmip5_piControl_ts_raw  = NaN(length(models_cmip5),252);
  cmip5_piControl_ts      = NaN(length(models_cmip5),252);
  for m = 1:length(models_cmip5)
    ['model = ' models_cmip5{m}]
    filein  = dir([pathin 'cmip5-ng/' vari '/' vari '_mon_' models_cmip5{m} '_piControl_*_g025.nc']);
    tmp0    = ncread([filein(1).folder '/' filein(1).name],[vari])*f;
    tmp0    = tmp0(:,:,end-252*12+1:end);
    for i = 1:length(tmp0)
      if strcmp(region,'AA')==1
        tmp1(i)    = nansum(nansum(squeeze(tmp0(jj1,ii1,i)).*weights1)); % do area-weighted mean
        tmp2(i)    = nansum(nansum(squeeze(tmp0(jj2,ii2,i)).*weights2)); % do area-weighted mean
        tmp(i)     = tmp1(i)-tmp2(i);
      else
        tmp(i)     = sum(sum(squeeze(tmp0(jj,ii,i)).*weights));
      end
    end
    if strcmp(seas,'annual')==1
      tmp    = am(tmp);
    else
      tmp    = seasmean(tmp,seas);
    end
    cmip5_piControl_ts_raw(m,:) = tmp;
    if wl > 1
      tmp = rm(tmp,wl);
    else
      tmp = tmp;
    end
    ref = nanmean(tmp);
    if strcmp(vari,'pr')==1
      cmip5_piControl_ts(m,:)   = ((tmp-ref)/ref)*100;
    else
      cmip5_piControl_ts(m,:)   = tmp-ref;
    end
    % cmip5_piControl_ts_var(m)   = nanvar(cmip5_piControl_ts(m,:));
    idx = ~isnan(cmip5_piControl_ts(m,:));
    cmip5_piControl_ts_var(m)   = nanvar(detrend(cmip5_piControl_ts(m,idx)));
  end

  % -- load CMIP6 piControl data --
  clear('tmp','tmp1','tmp2')
  cmip6_piControl_ts_raw  = NaN(length(models_cmip6),252);
  cmip6_piControl_ts      = NaN(length(models_cmip6),252);
  for m = 1:length(models_cmip6)
    ['model = ' models_cmip6{m}]
    filein  = dir([pathin 'cmip6-ng/' vari '/' vari '_mon_' models_cmip6{m} '_piControl_*_g025.nc']);
    tmp0    = ncread([filein(1).folder '/' filein(1).name],[vari])*f;
    tmp0    = tmp0(:,:,end-252*12+1:end);
    for i = 1:length(tmp0)
      if strcmp(region,'AA')==1
        tmp1(i)    = nansum(nansum(squeeze(tmp0(jj1,ii1,i)).*weights1)); % do area-weighted mean
        tmp2(i)    = nansum(nansum(squeeze(tmp0(jj2,ii2,i)).*weights2)); % do area-weighted mean
        tmp(i)     = tmp1(i)-tmp2(i);
      else
        tmp(i)     = sum(sum(squeeze(tmp0(jj,ii,i)).*weights));
      end
    end
    if strcmp(seas,'annual')==1
      tmp    = am(tmp);
    else
      tmp    = seasmean(tmp,seas);
    end
    cmip6_piControl_ts_raw(m,:) = tmp;
    if wl > 1
      tmp = rm(tmp,wl);
    else
      tmp = tmp;
    end
    ref = nanmean(tmp);
    if strcmp(vari,'pr')==1
      cmip6_piControl_ts(m,:)   = ((tmp-ref)/ref)*100;
    else
      cmip6_piControl_ts(m,:)   = tmp-ref;
    end
    % cmip6_piControl_ts_var(m)   = nanvar(cmip6_piControl_ts(m,:));
    idx = ~isnan(cmip6_piControl_ts(m,:));
    cmip6_piControl_ts_var(m)   = nanvar(detrend(cmip6_piControl_ts(m,idx)));
  end
end



% -- MPI-LE ensemble mean
scen = {'rcp26','rcp45','rcp85'};
for sc = 1:length(scen)
  filein  = [pathin 'mpi_lens/' comp '/' vari '/' vari '_' comp '_MPI-ESM_historical_' scen{sc} '_ensmean_185001-209912_g025.nc'];
  tmp0    = ncread([filein],[vari])*f;
  for i = 1:length(tmp0)
    if strcmp(region,'AA')==1
      tmp1(i)    = nansum(nansum(squeeze(tmp0(jj1,ii1,i)).*weights1)); % do area-weighted mean
      tmp2(i)    = nansum(nansum(squeeze(tmp0(jj2,ii2,i)).*weights2)); % do area-weighted mean
      tmp(i)     = tmp1(i)-tmp2(i);
    else
      tmp(i)     = sum(sum(squeeze(tmp0(jj,ii,i)).*weights));
    end
  end
  tmp     = tmp((start-1850)*12+1:(ende-1850+1)*12);
  if strcmp(seas,'annual')==1
    tmp    = am(tmp);
  else
    tmp    = seasmean(tmp,seas);
  end
  if wl > 1
    tmp = rm(tmp,wl);
  end
  ref = nanmean(tmp(refstart-start+1:refende-start+1));
  if strcmp(vari,'pr')==1
    tmp   = ((tmp-ref)/ref)*100;
  else
    tmp   = tmp-ref;
  end
  mpi_ensmean(sc,:) = tmp;
  idx = ~isnan(tmp);
  mpi_ensmean2(sc,:) = polyval(polyfit(time(idx),tmp(idx),4),time);
end
scen_u_mpi = nanvar(mpi_ensmean2,1);
% -- CanESM5-LE ensemble mean
scen = {'ssp126','ssp245','ssp370','ssp585'};
for sc = 1:length(scen)
  filein  = [pathin 'cmip6-ng/' vari '/' vari '_mon_CanESM5_' scen{sc} '_ensmean_g025.nc'];
  tmp0    = ncread([filein],[vari])*f;
  for i = 1:length(tmp0)
    if strcmp(region,'AA')==1
      tmp1(i)    = nansum(nansum(squeeze(tmp0(jj1,ii1,i)).*weights1)); % do area-weighted mean
      tmp2(i)    = nansum(nansum(squeeze(tmp0(jj2,ii2,i)).*weights2)); % do area-weighted mean
      tmp(i)     = tmp1(i)-tmp2(i);
    else
      tmp(i)     = sum(sum(squeeze(tmp0(jj,ii,i)).*weights));
    end
  end
  tmp     = tmp((start-1850)*12+1:(ende-1850+1)*12);
  if strcmp(seas,'annual')==1
    tmp    = am(tmp);
  else
    tmp    = seasmean(tmp,seas);
  end
  if wl > 1
    tmp = rm(tmp,wl);
  end
  ref = nanmean(tmp(refstart-start+1:refende-start+1));
  if strcmp(vari,'pr')==1
    tmp   = ((tmp-ref)/ref)*100;
  else
    tmp   = tmp-ref;
  end
  canesm5_ensmean(sc,:) = tmp;
  idx = ~isnan(tmp);
  canesm5_ensmean2(sc,:) = polyval(polyfit(time(idx),tmp(idx),4),time);
end
scen_u_canesm5 = nanvar(canesm5_ensmean2,1);


% -- write times series into txt file:
writematrix(le_ts_em,[pathin 'time_series/' vari '_' region '_le_' num2str(wl) 'yr.txt']);
writematrix(cmip5_rcp26_ts_em,[pathin 'time_series/' vari '_' region '_cmip5_rcp26_em_' num2str(wl) 'yr.txt']);
writematrix(cmip5_rcp45_ts_em,[pathin 'time_series/' vari '_' region '_cmip5_rcp45_em_' num2str(wl) 'yr.txt']);
writematrix(cmip5_rcp85_ts_em,[pathin 'time_series/' vari '_' region '_cmip5_rcp85_em_' num2str(wl) 'yr.txt']);
writematrix(cmip6_ssp126_ts_em,[pathin 'time_series/' vari '_' region '_cmip6_ssp126_em_' num2str(wl) 'yr.txt']);
writematrix(cmip6_ssp245_ts_em,[pathin 'time_series/' vari '_' region '_cmip6_ssp245_em_' num2str(wl) 'yr.txt']);
writematrix(cmip6_ssp370_ts_em,[pathin 'time_series/' vari '_' region '_cmip6_ssp370_em_' num2str(wl) 'yr.txt']);
writematrix(cmip6_ssp585_ts_em,[pathin 'time_series/' vari '_' region '_cmip6_ssp585_em_' num2str(wl) 'yr.txt']);

writematrix(cmip5_rcp26_ts,[pathin 'time_series/' vari '_' region '_cmip5_rcp26_' num2str(wl) 'yr.txt']);
writematrix(cmip5_rcp45_ts,[pathin 'time_series/' vari '_' region '_cmip5_rcp45_' num2str(wl) 'yr.txt']);
writematrix(cmip5_rcp85_ts,[pathin 'time_series/' vari '_' region '_cmip5_rcp85_' num2str(wl) 'yr.txt']);
writematrix(cmip6_ssp126_ts,[pathin 'time_series/' vari '_' region '_cmip6_ssp126_' num2str(wl) 'yr.txt']);
writematrix(cmip6_ssp245_ts,[pathin 'time_series/' vari '_' region '_cmip6_ssp245_' num2str(wl) 'yr.txt']);
writematrix(cmip6_ssp370_ts,[pathin 'time_series/' vari '_' region '_cmip6_ssp370_' num2str(wl) 'yr.txt']);
writematrix(cmip6_ssp585_ts,[pathin 'time_series/' vari '_' region '_cmip6_ssp585_' num2str(wl) 'yr.txt']);



%% -- UNCERTAINTIES --
% -- scenario uncertainty as in Brekke and Barsugli (2013) --
clear('tmp')
for m = 1:length(models_cmip5)
  tmp(m,:) = nanvar([cmip5_rcp85_ts_em(m,:); cmip5_rcp45_ts_em(m,:); cmip5_rcp26_ts_em(m,:)]);
end
scen_u_cmip5_brekke = nanmean(tmp);
clear('tmp')
for m = 1:length(models_cmip6)
  tmp(m,:) = nanvar([cmip6_ssp585_ts_em(m,:); cmip6_ssp370_ts_em(m,:); cmip6_ssp245_ts_em(m,:); cmip6_ssp126_ts_em(m,:)]);
end
scen_u_cmip6_brekke = nanmean(tmp);
% -- normal HS09 uncertainties --
scen_u_cmip5      = nanvar([nanmean(cmip5_rcp85_ts_em,1); nanmean(cmip5_rcp45_ts_em,1); nanmean(cmip5_rcp26_ts_em,1)]);
% le_cmip5_rcp85_ts_em = [cmip5_rcp85_ts_em(models_le_cmip5_id,:); mpi_ensmean(3,:)]
% le_cmip5_rcp45_ts_em = [cmip5_rcp45_ts_em(models_le_cmip5_id,:); mpi_ensmean(2,:)]
% le_cmip5_rcp26_ts_em = [cmip5_rcp26_ts_em(models_le_cmip5_id,:); mpi_ensmean(1,:)]
le_cmip5_rcp85_ts_em = [cmip5_rcp85_ts_em(models_le_cmip5_id,:)];
le_cmip5_rcp45_ts_em = [cmip5_rcp45_ts_em(models_le_cmip5_id,:)];
le_cmip5_rcp26_ts_em = [cmip5_rcp26_ts_em(models_le_cmip5_id,:)];
scen_u_le_cmip5   = nanvar([nanmean(le_cmip5_rcp85_ts_em,1); nanmean(le_cmip5_rcp45_ts_em,1); nanmean(le_cmip5_rcp26_ts_em,1)]);
% scen_u_cmip5_w    = nanvar([nanmean(weights_cmip5.*cmip5_rcp85_ts_em,1); nanmean(weights_cmip5.*cmip5_rcp45_ts_em,1); nanmean(weights_cmip5.*cmip5_rcp26_ts_em,1)]);
scen_u_cmip5_w    = scen_u_cmip5;
% % -- CMIP6 scen uncertainty with 3 scenarios:
% scen_u_cmip6      = nanvar([nanmean(cmip6_ssp585_ts_em,1); nanmean(cmip6_ssp245_ts_em,1); nanmean(cmip6_ssp126_ts_em,1)]); % using 3 scenarios
% -- CMIP6 scen uncertainty with 4 scenarios:
scen_u_cmip6      = nanvar([nanmean(cmip6_ssp585_ts_em,1); nanmean(cmip6_ssp370_ts_em,1); nanmean(cmip6_ssp245_ts_em,1); nanmean(cmip6_ssp126_ts_em,1)]); % using 4 scenarios

% scen_u_cmip6_w    = nanvar([nanmean(weights_cmip6.*cmip6_ssp585_ts_em,1); nanmean(weights_cmip6.*cmip6_ssp370_ts_em,1); nanmean(weights_cmip6.*cmip6_ssp245_ts_em,1); nanmean(weights_cmip6.*cmip6_ssp126_ts_em,1)]);
scen_u_cmip6_w    = scen_u_cmip6;
[a,i] = min(scen_u_cmip5);
scen_u_cmip5(1:i) = 0;
scen_u_cmip5_brekke(1:i) = 0;
scen_u_le_cmip5(1:i) = 0;
[a,i] = min(scen_u_cmip6);
scen_u_cmip6(1:i) = 0;
scen_u_cmip6_brekke(1:i) = 0;

model_u_le        = nanvar(le_ts_em,1);
model_u_le_hs     = nanvar(squeeze(le_ts_em_hs(:,1,:)),1);
% -- treat ensemble member as individual models:
for m = 1:length(models)
  model_u_le_hs_all(m,:) = nanvar(squeeze(le_ts_em_hs(m,:,:)),1);
end
% -- subsample MPI (100) as if it were EC-EARTH (16) to calculate the forced response/ensemble mean;
% -- then calculate the variance across these estimates as another pseudo model uncertainty
m = 7; % MPI model id
ensmem_subsample = 30; % 16 30
clear('tmp1')
for i = 1:100
  idx = randi(100,[1,ensmem_subsample]);
  tmp0 = rm(nanmean(squeeze(le_ts_anom(m,idx,:))),wl);
  tmp1(i,:) = tmp0 - nanmean(tmp0(refstart-start+1:refende-start+1));
end
model_u_le_mpi_subsample = nanvar(tmp1,1);

model_u_cmip5     = nanvar(cmip5_rcp85_ts_em,1);
model_u_le_cmip5  = nanvar(le_cmip5_rcp85_ts_em,1);
model_u_cmip5_w   = nanvar(weights_cmip5.*cmip5_rcp85_ts_em,1);
model_u_cmip6     = nanvar(cmip6_ssp585_ts_em,1);
model_u_cmip6_w   = nanvar(weights_cmip6.*cmip6_ssp585_ts_em,1);
for i = 1:100
  model_u_cmip5_w_bs(i,:) = nanvar(squeeze(weights_cmip5_bs(i,:,:)).*cmip5_rcp85_ts_em,1);
  model_u_cmip6_w_bs(i,:) = nanvar(squeeze(weights_cmip6_bs(i,:,:)).*cmip6_ssp585_ts_em,1);
end

% -- what HS09 did, namely calculate the mean across model uncertainties from the different scenarios):
model_u_cmip5_hs  = nanmean([nanvar(cmip5_rcp85_ts_em,1); nanvar(cmip5_rcp45_ts_em,1); nanvar(cmip5_rcp26_ts_em,1)]);
model_u_cmip6_hs  = nanmean([nanvar(cmip6_ssp585_ts_em,1); nanvar(cmip6_ssp370_ts_em,1); nanvar(cmip6_ssp245_ts_em,1); nanvar(cmip6_ssp126_ts_em,1)]);

if wl == 1
  int_u_le_mean     = nanmean(le_ts_var,1);
  int_u_le_max      = max(le_ts_var);
  int_u_le_min      = min(le_ts_var);
else
  int_u_le_mean     = rm(nanmean(le_ts_var,1),10);
  int_u_le_max      = rm(max(le_ts_var),10);
  int_u_le_min      = rm(min(le_ts_var),10);
end
% int_u_le_mean     = nanmean(le_ts_var1,1);
% int_u_le_max      = max(le_ts_var1);
% int_u_le_min      = min(le_ts_var1);
int_u_le_fixed    = int_u_le_mean*0 + nanmean(int_u_le_mean(1:refende-start));
int_u_le_hs_mean  = nanmean(le_ts_var_hs1(:,1)); % use just first ensemble member
int_u_le_hs_max   = max(le_ts_var_hs1(:,1)); % use just first ensemble member
int_u_le_hs_min   = min(le_ts_var_hs1(:,1)); % use just first ensemble member
int_u_cmip5_mean  = nanmean(cmip5_rcp85_ts_var1);
int_u_le_cmip5_mean  = nanmean(cmip5_rcp85_ts_var1(models_le_cmip5_id));
int_u_cmip5_w_mean= nanmean(weights_cmip5(:,1)'.*cmip5_rcp85_ts_var1);
int_u_cmip5_max   = max(cmip5_rcp85_ts_var1);
int_u_cmip5_min   = min(cmip5_rcp85_ts_var1);
int_u_cmip6_mean  = nanmean(cmip6_ssp585_ts_var1);
int_u_cmip6_w_mean= nanmean(weights_cmip6(:,1)'.*cmip6_ssp585_ts_var1);
int_u_cmip6_max   = max(cmip6_ssp585_ts_var1);
int_u_cmip6_min   = min(cmip6_ssp585_ts_var1);

total_u_le_mean       = scen_u_cmip5 + model_u_le + int_u_le_mean;
total_u_le_mean_brekke= scen_u_cmip5_brekke + model_u_le + int_u_le_mean;
total_u_le_cmip5_mean = scen_u_le_cmip5 + model_u_le_cmip5 + int_u_le_cmip5_mean;
total_u_le_mean_mpi   = scen_u_mpi + model_u_le + int_u_le_mean;
total_u_le_mean_canesm5 = scen_u_canesm5 + model_u_le + int_u_le_mean;
total_u_le_mpi_subsample = scen_u_cmip5 + model_u_le_mpi_subsample + int_u_le_mean;
total_u_le_fixed      = scen_u_cmip5 + model_u_le + int_u_le_fixed;
total_u_le_max        = scen_u_cmip5 + model_u_le + int_u_le_max;
total_u_le_min        = scen_u_cmip5 + model_u_le + int_u_le_min;
total_u_le_hs_mean    = scen_u_cmip5 + model_u_le_hs + int_u_le_hs_mean;
total_u_le_hs_all_mean= scen_u_cmip5 + model_u_le_hs_all + int_u_le_hs_mean;
total_u_le_mpi_subsample = scen_u_cmip5 + model_u_le_mpi_subsample + int_u_le_mean;

total_u_cmip5_mean        = scen_u_cmip5 + model_u_cmip5 + int_u_cmip5_mean;
total_u_cmip5_mean_brekke = scen_u_cmip5_brekke + model_u_cmip5 + int_u_cmip5_mean;
total_u_cmip5_mean_mpi    = scen_u_mpi + model_u_cmip5 + int_u_cmip5_mean;
% total_u_cmip5_w_mean    = scen_u_cmip5_w + model_u_cmip5_w + int_u_cmip5_w_mean;
total_u_cmip5_w_mean      = scen_u_cmip5 + model_u_cmip5_w + int_u_cmip5_mean;
total_u_cmip5_w_bs_mean   = scen_u_cmip5 + model_u_cmip5_w_bs + int_u_cmip5_mean;
total_u_cmip6_mean        = scen_u_cmip6 + model_u_cmip6 + int_u_cmip6_mean;
total_u_cmip6_mean_brekke = scen_u_cmip6_brekke + model_u_cmip6 + int_u_cmip6_mean;
total_u_cmip6_mean_canesm5 = scen_u_canesm5 + model_u_cmip6 + int_u_cmip6_mean;
% total_u_cmip6_w_mean    = scen_u_cmip6_w + model_u_cmip6_w + int_u_cmip6_w_mean;
total_u_cmip6_w_mean      = scen_u_cmip6 + model_u_cmip6_w + int_u_cmip6_mean;
total_u_cmip6_w_bs_mean   = scen_u_cmip6 + model_u_cmip6_w_bs + int_u_cmip6_mean;
total_u_cmip5_mean_hs     = scen_u_cmip5 + model_u_cmip5_hs + int_u_cmip5_mean;
total_u_cmip6_mean_hs     = scen_u_cmip6 + model_u_cmip6_hs + int_u_cmip6_mean;

% -- fractional uncertainties
scen_u_frac_le        = (scen_u_cmip5./total_u_le_mean)*100;
scen_u_frac_le_brekke = (scen_u_cmip5_brekke./total_u_le_mean_brekke)*100;
scen_u_frac_le_cmip5  = (scen_u_le_cmip5./total_u_le_cmip5_mean)*100;
scen_u_frac_le_mpi    = (scen_u_mpi./total_u_le_mean_mpi)*100;
scen_u_frac_le_canesm5 = (scen_u_canesm5./total_u_le_mean_canesm5)*100;
scen_u_frac_le_mpi_subsample = (scen_u_le_cmip5./total_u_le_mpi_subsample)*100;
scen_u_frac_le_fixed  = (scen_u_cmip5./total_u_le_fixed)*100;
scen_u_frac_le_max    = (scen_u_cmip5./total_u_le_max)*100;
scen_u_frac_le_min    = (scen_u_cmip5./total_u_le_min)*100;
scen_u_frac_le_hs     = (scen_u_cmip5./total_u_le_hs_mean)*100;
scen_u_frac_le_hs_all = (scen_u_cmip5./total_u_le_hs_all_mean)*100;
scen_u_frac_cmip5     = (scen_u_cmip5./total_u_cmip5_mean)*100;
scen_u_frac_cmip5_brekke = (scen_u_cmip5_brekke./total_u_cmip5_mean_brekke)*100;
scen_u_frac_cmip5_mpi = (scen_u_mpi./total_u_cmip5_mean_mpi)*100;
% scen_u_frac_cmip5_w   = (scen_u_cmip5_w./total_u_cmip5_w_mean)*100;
scen_u_frac_cmip5_w   = (scen_u_cmip5./total_u_cmip5_w_mean)*100;
scen_u_frac_cmip5_w_bs= (scen_u_cmip5./total_u_cmip5_w_bs_mean)*100;
scen_u_frac_cmip6     = (scen_u_cmip6./total_u_cmip6_mean)*100;
scen_u_frac_cmip6_brekke = (scen_u_cmip6_brekke./total_u_cmip6_mean_brekke)*100;
scen_u_frac_cmip6_canesm5 = (scen_u_canesm5./total_u_cmip6_mean_canesm5)*100;
% scen_u_frac_cmip6_w   = (scen_u_cmip6_w./total_u_cmip6_w_mean)*100;
scen_u_frac_cmip6_w   = (scen_u_cmip6./total_u_cmip6_w_mean)*100;
scen_u_frac_cmip6_w_bs= (scen_u_cmip6./total_u_cmip6_w_bs_mean)*100;
scen_u_frac_cmip5_hs  = (scen_u_cmip5./total_u_cmip5_mean_hs)*100;
scen_u_frac_cmip6_hs  = (scen_u_cmip6./total_u_cmip6_mean_hs)*100;
model_u_frac_le       = (model_u_le./total_u_le_mean)*100;
model_u_frac_le_brekke= (model_u_le./total_u_le_mean_brekke)*100;
model_u_frac_le_cmip5 = (model_u_le_cmip5./total_u_le_cmip5_mean)*100;
model_u_frac_le_mpi   = (model_u_le./total_u_le_mean_mpi)*100;
model_u_frac_le_canesm5 = (model_u_le./total_u_le_mean_canesm5)*100;
model_u_frac_le_mpi_subsample = (model_u_le_mpi_subsample./total_u_le_mpi_subsample)*100;
model_u_frac_le_fixed = (model_u_le./total_u_le_fixed)*100;
model_u_frac_le_max   = (model_u_le./total_u_le_max)*100;
model_u_frac_le_min   = (model_u_le./total_u_le_min)*100;
model_u_frac_le_hs    = (model_u_le_hs./total_u_le_hs_mean)*100;
model_u_frac_le_hs_all= (model_u_le_hs_all./total_u_le_hs_all_mean)*100;
model_u_frac_cmip5    = (model_u_cmip5./total_u_cmip5_mean)*100;
model_u_frac_cmip5_brekke = (model_u_cmip5./total_u_cmip5_mean_brekke)*100;
model_u_frac_cmip5_mpi = (model_u_cmip5./total_u_cmip5_mean_mpi)*100;
model_u_frac_cmip5_w  = (model_u_cmip5_w./total_u_cmip5_w_mean)*100;
model_u_frac_cmip5_w_bs  = (model_u_cmip5_w_bs./total_u_cmip5_w_bs_mean)*100;
model_u_frac_cmip6    = (model_u_cmip6./total_u_cmip6_mean)*100;
model_u_frac_cmip6_brekke = (model_u_cmip6./total_u_cmip6_mean_brekke)*100;
model_u_frac_cmip6_canesm5 = (model_u_cmip6./total_u_cmip6_mean_canesm5)*100;
model_u_frac_cmip6_w  = (model_u_cmip6_w./total_u_cmip6_w_mean)*100;
model_u_frac_cmip6_w_bs  = (model_u_cmip6_w_bs./total_u_cmip6_w_bs_mean)*100;
model_u_frac_cmip5_hs = (model_u_cmip5_hs./total_u_cmip5_mean_hs)*100;
model_u_frac_cmip6_hs = (model_u_cmip6_hs./total_u_cmip6_mean_hs)*100;
int_u_frac_le_mean        = (int_u_le_mean./total_u_le_mean)*100;
int_u_frac_le_mean_brekke = (int_u_le_mean./total_u_le_mean_brekke)*100;
int_u_frac_le_cmip5_mean  = (int_u_le_mean./total_u_le_cmip5_mean)*100;
int_u_frac_le_mean_mpi    = (int_u_le_mean./total_u_le_mean_mpi)*100;
int_u_frac_le_mean_canesm5= (int_u_le_mean./total_u_le_mean_canesm5)*100;
int_u_frac_le_mpi_subsample  = (int_u_le_mean./total_u_le_mpi_subsample)*100;
int_u_frac_le_fixed       = (int_u_le_fixed./total_u_le_fixed)*100;
int_u_frac_le_max         = (int_u_le_max./total_u_le_max)*100;
int_u_frac_le_min         = (int_u_le_min./total_u_le_min)*100;
int_u_frac_le_hs_mean     = (int_u_le_hs_mean./total_u_le_hs_mean)*100;
int_u_frac_le_hs_all_mean = (int_u_le_hs_mean./total_u_le_hs_all_mean)*100;
int_u_frac_cmip5_mean     = (int_u_cmip5_mean./total_u_cmip5_mean)*100;
int_u_frac_cmip5_mean_brekke = (int_u_cmip5_mean./total_u_cmip5_mean_brekke)*100;
int_u_frac_cmip5_mean_mpi = (int_u_cmip5_mean./total_u_cmip5_mean_mpi)*100;
% int_u_frac_cmip5_w_mean   = (int_u_cmip5_w_mean./total_u_cmip5_w_mean)*100;
int_u_frac_cmip5_w_mean   = (int_u_cmip5_mean./total_u_cmip5_w_mean)*100;
int_u_frac_cmip6_mean     = (int_u_cmip6_mean./total_u_cmip6_mean)*100;
int_u_frac_cmip6_mean_brekke = (int_u_cmip6_mean./total_u_cmip6_mean_brekke)*100;
int_u_frac_cmip6_mean_canesm5 = (int_u_cmip6_mean./total_u_cmip6_mean_canesm5)*100;
% int_u_frac_cmip6_w_mean   = (int_u_cmip6_w_mean./total_u_cmip6_w_mean)*100;
int_u_frac_cmip6_w_mean   = (int_u_cmip6_mean./total_u_cmip6_w_mean)*100;
int_u_frac_cmip6_w_bs_mean   = (int_u_cmip6_mean./total_u_cmip6_w_bs_mean)*100;

% -- load global mean temperature change (forced response, NOTE: added a warming-to-date "tas_obs_wtd" estimate from observations):
% offset = tas_obs_wtd;
offset = 0;
tic
tasg_le_ts_em            = load([pathin 'time_series/tas_global_le_' num2str(wl) 'yr.txt'])' + offset;
tasg_cmip5_rcp26_ts_em   = load([pathin 'time_series/tas_global_cmip5_rcp26_em_' num2str(wl) 'yr.txt'])' + offset;
tasg_cmip5_rcp45_ts_em   = load([pathin 'time_series/tas_global_cmip5_rcp45_em_' num2str(wl) 'yr.txt'])' + offset;
tasg_cmip5_rcp85_ts_em   = load([pathin 'time_series/tas_global_cmip5_rcp85_em_' num2str(wl) 'yr.txt'])' + offset;
tasg_cmip6_ssp126_ts_em  = load([pathin 'time_series/tas_global_cmip6_ssp126_em_' num2str(wl) 'yr.txt'])' + offset;
tasg_cmip6_ssp245_ts_em  = load([pathin 'time_series/tas_global_cmip6_ssp245_em_' num2str(wl) 'yr.txt'])' + offset;
tasg_cmip6_ssp370_ts_em  = load([pathin 'time_series/tas_global_cmip6_ssp370_em_' num2str(wl) 'yr.txt'])' + offset;
tasg_cmip6_ssp585_ts_em  = load([pathin 'time_series/tas_global_cmip6_ssp585_em_' num2str(wl) 'yr.txt'])' + offset;
% -- sort variable according to tasg:
tmp   = [tasg_le_ts_em(:); tasg_cmip5_rcp85_ts_em(:); tasg_cmip6_ssp585_ts_em(:)] + offset;
range = [min(tmp):(max(tmp)-min(tmp))/100:max(tmp)] + offset;
clear('le_ts_em_vs_tasg','le_var_vs_tasg','cmip5_rcp85_ts_em_vs_tasg','cmip6_ssp585_ts_em_vs_tasg')
le_ts_em_vs_tasg            = NaN(length(models),length(range));
le_var_vs_tasg              = NaN(length(models),length(range));
cmip5_rcp85_ts_em_vs_tasg   = NaN(length(models_cmip5),length(range));
cmip6_ssp585_ts_em_vs_tasg  = NaN(length(models_cmip6),length(range));
for i = 1:length(range)
  for m = 1:length(models)
    [null,idx] = min(abs(tasg_le_ts_em(:,m)-range(i)));
    le_ts_em_vs_tasg(m,i) = le_ts_em(m,idx);
    le_var_vs_tasg(m,i) = nanvar(le_ts_tmp(m,:,idx));
    % le_var_vs_tasg(m,i) = int_u_le_mean(idx);
  end
  for s = 1:length(scen_cmip5)
    for m = 1:length(models_cmip5)
      eval(['tmp = islocalmin((abs(tasg_cmip5_' scen_cmip5{s} '_ts_em(:,m)-range(i))));']);
      tmp = find(tmp==1);
      if isempty(tmp)==1
        eval(['[null,idx] = min(abs(tasg_cmip5_' scen_cmip5{s} '_ts_em(:,m)-range(i)));']);
      else
        idx = tmp(1);
      end
      eval(['cmip5_' scen_cmip5{s} '_ts_em_vs_tasg(m,i) = cmip5_' scen_cmip5{s} '_ts_em(m,idx);']);
    end
  end
  for s = 1:length(scen_cmip6)
    for m = 1:length(models_cmip6)
      eval(['[null,idx] = min(abs(tasg_cmip6_' scen_cmip6{s} '_ts_em(:,m)-range(i)));']);
      eval(['cmip6_' scen_cmip6{s} '_ts_em_vs_tasg(m,i) = cmip6_' scen_cmip6{s} '_ts_em(m,idx);']);
    end
  end
end


model_u_le_tasg           = nanvar(le_ts_em_vs_tasg);
model_u_cmip5_rcp26_tasg  = nanvar(cmip5_rcp26_ts_em_vs_tasg);
model_u_cmip5_rcp45_tasg  = nanvar(cmip5_rcp45_ts_em_vs_tasg);
model_u_cmip5_rcp85_tasg  = nanvar(cmip5_rcp85_ts_em_vs_tasg);
model_u_cmip6_ssp126_tasg = nanvar(cmip6_ssp126_ts_em_vs_tasg);
model_u_cmip6_ssp245_tasg = nanvar(cmip6_ssp245_ts_em_vs_tasg);
model_u_cmip6_ssp370_tasg = nanvar(cmip6_ssp370_ts_em_vs_tasg);
model_u_cmip6_ssp585_tasg = nanvar(cmip6_ssp585_ts_em_vs_tasg);
int_u_le_tasg_mean        = nanmean(le_var_vs_tasg);
int_u_cmip5_tasg_mean     = int_u_cmip5_mean;
int_u_cmip6_tasg_mean     = int_u_cmip6_mean;
total_u_le_tasg_mean      = model_u_le_tasg + int_u_le_tasg_mean;
total_u_cmip5_rcp26_tasg_mean   = model_u_cmip5_rcp26_tasg + int_u_cmip5_tasg_mean;
total_u_cmip5_rcp45_tasg_mean   = model_u_cmip5_rcp45_tasg + int_u_cmip5_tasg_mean;
total_u_cmip5_rcp85_tasg_mean   = model_u_cmip5_rcp85_tasg + int_u_cmip5_tasg_mean;
total_u_cmip6_ssp126_tasg_mean   = model_u_cmip6_ssp126_tasg + int_u_cmip6_tasg_mean;
total_u_cmip6_ssp245_tasg_mean   = model_u_cmip6_ssp245_tasg + int_u_cmip6_tasg_mean;
total_u_cmip6_ssp370_tasg_mean   = model_u_cmip6_ssp370_tasg + int_u_cmip6_tasg_mean;
total_u_cmip6_ssp585_tasg_mean   = model_u_cmip6_ssp585_tasg + int_u_cmip6_tasg_mean;

model_u_frac_le_tasg      = (model_u_le_tasg./total_u_le_tasg_mean)*100;
model_u_frac_cmip5_rcp26_tasg   = (model_u_cmip5_rcp26_tasg./total_u_cmip5_rcp26_tasg_mean)*100;
model_u_frac_cmip5_rcp45_tasg   = (model_u_cmip5_rcp45_tasg./total_u_cmip5_rcp45_tasg_mean)*100;
model_u_frac_cmip5_rcp85_tasg   = (model_u_cmip5_rcp85_tasg./total_u_cmip5_rcp85_tasg_mean)*100;
model_u_frac_cmip6_ssp126_tasg   = (model_u_cmip6_ssp126_tasg./total_u_cmip6_ssp126_tasg_mean)*100;
model_u_frac_cmip6_ssp126_tasg   = (model_u_cmip6_ssp126_tasg./total_u_cmip6_ssp126_tasg_mean)*100;
model_u_frac_cmip6_ssp126_tasg   = (model_u_cmip6_ssp126_tasg./total_u_cmip6_ssp126_tasg_mean)*100;
model_u_frac_cmip6_ssp126_tasg   = (model_u_cmip6_ssp126_tasg./total_u_cmip6_ssp126_tasg_mean)*100;
int_u_frac_le_tasg        = (int_u_le_tasg_mean./total_u_le_tasg_mean)*100;
int_u_frac_cmip5_rcp26_tasg     = (int_u_cmip5_tasg_mean./total_u_cmip5_rcp26_tasg_mean)*100;
int_u_frac_cmip5_rcp45_tasg     = (int_u_cmip5_tasg_mean./total_u_cmip5_rcp45_tasg_mean)*100;
int_u_frac_cmip5_rcp85_tasg     = (int_u_cmip5_tasg_mean./total_u_cmip5_rcp85_tasg_mean)*100;
int_u_frac_cmip6_ssp126_tasg     = (int_u_cmip6_tasg_mean./total_u_cmip6_ssp126_tasg_mean)*100;
int_u_frac_cmip6_ssp245_tasg     = (int_u_cmip6_tasg_mean./total_u_cmip6_ssp245_tasg_mean)*100;
int_u_frac_cmip6_ssp370_tasg     = (int_u_cmip6_tasg_mean./total_u_cmip6_ssp370_tasg_mean)*100;
int_u_frac_cmip6_ssp585_tasg     = (int_u_cmip6_tasg_mean./total_u_cmip6_ssp585_tasg_mean)*100;
toc





% --- PLOTTING ---------------------------------------------------------------
close all

cols = [255 0 0;
        255 160 16;
        255 224 32;
        0 192 0;
        80 208 255;
        0 32 255;
        160 32 255]/255;

hs09_cols = ...
[53 74 161;...
 255 110 4;...
 0 127 60]/255;
hs09_cols_light = ...
[164 180 245;...
 252 210 179;...
 172 232 200]/255;


rcp_cols = [217,37,42;...
           232,126,63;...
           133,177,212;...
           51,74,141]/255;
rcp_cols_light = ...
[235,144,115;...
243,183,136;...
189,211,227;...
138,141,185]/255;

if wl > 1
  xlim = [refende+1 2099-round(wl/2)];
  xlim0 = [1950+round(wl/2) 2099-round(wl/2)];
else
  xlim = [refende+1 2099-round(wl/2)];
  xlim0 = [1950+round(wl/2) 2099-round(wl/2)];
end
if strcmp(vari,'tas')==1
  yincr = 1;
  if refstart == 1995
    ylim = [-1.25 5.5];
  else
    ylim = [0 5.5];
  end
else
  yincr = 2;
  ylim = [-2.9 10];
end

tl = 0.015; % tick length
sh = 0.05; % horizontal space between panels
sv = 0.04; % vertical space between panels



% ----------------------------------
if plotS2 == 1;
% -- plot --------
close all


  figure1 = figure;
  set(figure1, 'units', 'centimeters', 'pos', [10 10 30 8]) % on home screen

  % -- SMILEs:
  subplot(1,3,1)
  subaxis(1,3,1, 'sh', sh, 'sv', sv)
  hold on
  % title([region ' ' vari ' ' seas ' ' num2str(wl) '-yr means (' units ')' char(10) ['relative to ' num2str(refstart) '-' num2str(refende)]],'Interpreter','none')%,'FontSize',10)
  title('SMILEs')
  idx = find(~isnan(total_u_le_mean));
  x = time(idx);
  y0        = x*0;
  ym_le     = model_u_frac_le(idx); % model u
  yms_le    = model_u_frac_le(idx) + scen_u_frac_le(idx); % model u plus scen u
  ym_le_hs  = model_u_frac_le_hs(idx); % model u
  yms_le_hs = model_u_frac_le_hs(idx) + scen_u_frac_le_hs(idx); % model u plus scen u
  ytop      = time(idx)*0+100;
  le_fixed  = model_u_frac_le_fixed(idx) + scen_u_frac_le_fixed(idx);
  le_max1   = model_u_frac_le_max(idx);
  le_min1   = model_u_frac_le_min(idx);
  le_max2   = model_u_frac_le_max(idx) + scen_u_frac_le_max(idx);
  le_min2   = model_u_frac_le_min(idx) + scen_u_frac_le_min(idx);
  ym_cmip5  = model_u_frac_cmip5(idx); % model u
  yms_cmip5 = model_u_frac_cmip5(idx) + scen_u_frac_cmip5(idx); % model u plus scen u
  ym_mpi    = model_u_frac_le_mpi(idx);
  yms_mpi   = model_u_frac_le_mpi(idx)+scen_u_frac_le_mpi(idx);
  ym_le_cmip5    = model_u_frac_le_cmip5(idx);
  yms_le_cmip5   = model_u_frac_le_cmip5(idx)+scen_u_frac_le_cmip5(idx);
  h2 = patch([x fliplr(x)],[y0 fliplr(ym_le)],[53 74 161]/255,'LineWidth',.5,'Edgecolor','none')
  h1 = patch([x fliplr(x)],[yms_le fliplr(ytop)],[255 110 4]/255,'LineWidth',.5,'Edgecolor','none')
  h3 = patch([x fliplr(x)],[ym_le fliplr(yms_le)],[0 127 60]/255,'LineWidth',.5,'Edgecolor','none');
  patch([x fliplr(x)],[ym_le fliplr(yms_le)],[0 127 60]/255,'LineWidth',.5)%,'Edgecolor','none')
  hold on
  h5 = plot(x,[ym_mpi; yms_mpi],'k:','LineWidth',2)
  hp = patch([x fliplr(x)],[ym_mpi fliplr(yms_mpi)],[1 1 1],'LineWidth',.5,'Edgecolor','none','FaceAlpha',.1);
  hatchfill(hp,'single',65,5,[0 127 60]/255);
  h6 = plot(x,[ym_le_cmip5; yms_le_cmip5],'k--','LineWidth',2)
  hp2 = patch([x fliplr(x)],[ym_le_cmip5 fliplr(yms_le_cmip5)],[1 1 1],'LineWidth',.5,'Edgecolor','none','FaceAlpha',.1);
  hatchfill(hp2,'single',25,10,[0 127 60]/255);
  set(gca,'Layer','top','XLim',xlim,'YLim',[0 100],'TickLength',[tl tl])
  if strcmp(region,'global')==1 & strcmp(seas,'annual')==1 % & strcmp(vari,'tas')==1
    % legend([h1(1) h3(1) h2(1)],'Int. variability','Scenario (CMIP5)','Model','Location','SouthWest','FontSize',8)
    legend([h1(1) h3(1) h2(1) h5(1) h6(1)],'Int. variability','Scenario (CMIP5)','Model','Scenario (MPI-LE)','Scenario (SMILEs)','Location','SouthWest','FontSize',8)
    % legend([h1(1) h3(1) h2(1) h5(1)],'Int. variability','Scenario (CMIP5)','Model','Scenario (MPI-LE)','Location','SouthWest','FontSize',8)
  end
  ylabel({'Fractional contribution','to total uncertainty (%)'})
  xlabel('Time (Year)')
  box on
  text(2085,.1*100,'(a)')

  % -- CMIP5:
  subplot(1,3,2)
  subaxis(1,3,2, 'sh', sh, 'sv', sv)
  hold on
  title('CMIP5')%,'FontSize',10)
  idx       = find(~isnan(total_u_cmip5_mean));
  x         = time(idx);
  y0        = x*0;
  % ym_cmip5  = model_u_frac_cmip5_hs(idx); % model u
  % yms_cmip5 = model_u_frac_cmip5_hs(idx) + scen_u_frac_cmip5_hs(idx); % model u plus scen u
  if constrain == 0
    ym_cmip5  = model_u_frac_cmip5(idx); % model u
    yms_cmip5 = model_u_frac_cmip5(idx) + scen_u_frac_cmip5(idx); % model u plus scen u
  else
    ym_cmip5  = model_u_frac_cmip5_w(idx); % model u
    yms_cmip5 = model_u_frac_cmip5_w(idx) + scen_u_frac_cmip5_w(idx); % model u plus scen u
  end
  ym_cmip5_mpi  = model_u_frac_cmip5_mpi(idx); % model u
  yms_cmip5_mpi = model_u_frac_cmip5_mpi(idx) + scen_u_frac_cmip5_mpi(idx); % model u plus scen u
  ytop      = time(idx)*0+100;
  h2 = patch([x fliplr(x)],[y0 fliplr(ym_cmip5)],[53 74 161]/255,'LineWidth',.5,'Edgecolor','none')
  h1 = patch([x fliplr(x)],[yms_cmip5 fliplr(ytop)],[255 110 4]/255,'LineWidth',.5,'Edgecolor','none')
  h3 = patch([x fliplr(x)],[ym_cmip5 fliplr(yms_cmip5)],[0 127 60]/255,'LineWidth',.5,'Edgecolor','none')
  patch([x fliplr(x)],[ym_cmip5 fliplr(yms_cmip5)],[0 127 60]/255,'LineWidth',.5)%,'Edgecolor','none')
  hold on
  set(gca,'Layer','top','XLim',xlim,'YLim',[0 100],'TickLength',[tl tl])
  h5 = plot(x,ym_cmip5_mpi,'k:','LineWidth',2)
  plot(x,yms_cmip5_mpi,'k:','LineWidth',2)
  hp = patch([x fliplr(x)],[ym_cmip5_mpi fliplr(yms_cmip5_mpi)],[1 1 1],'LineWidth',.5,'Edgecolor','none','FaceAlpha',.1);
  hatchfill(hp,'single',49,8,[0 127 60]/255);
  if strcmp(region,'global')==1 & strcmp(seas,'annual')==1 % & strcmp(vari,'tas')==1
    % legend([h1(1) h3(1) h2(1)],'Int. variability','Scenario','Model','Location','SouthWest','FontSize',8)
    legend([h1(1) h3(1) h2(1) h5(1)],'Int. variability','Scenario','Model','Scenario (MPI-LE)','Location','SouthWest','FontSize',8)
  end
  % ylabel({'Fractional contribution','to total uncertainty (%)'})
  xlabel('Time (Year)')
  box on
  text(2085,.1*100,'(b)')


  % -- CMIP6:
  subplot(1,3,3)
  subaxis(1,3,3, 'sh', sh, 'sv', sv)
  hold on
  title('CMIP6')%,'FontSize',10)
  idx       = find(~isnan(total_u_cmip6_mean));
  x         = time(idx);
  y0        = x*0;
  if constrain == 0
    % ym_cmip6  = model_u_frac_cmip6_hs(idx); % model u
    % yms_cmip6 = model_u_frac_cmip6_hs(idx) + scen_u_frac_cmip6_hs(idx); % model u plus scen u
    ym_cmip6  = model_u_frac_cmip6(idx); % model u
    yms_cmip6 = model_u_frac_cmip6(idx) + scen_u_frac_cmip6(idx); % model u plus scen u
  else
    ym_cmip6  = model_u_frac_cmip6_w(idx); % model u
    yms_cmip6 = model_u_frac_cmip6_w(idx) + scen_u_frac_cmip6_w(idx); % model u plus scen u
  end
  ym_cmip6_canesm5  = model_u_frac_cmip6_canesm5(idx); % model u
  yms_cmip6_canesm5 = model_u_frac_cmip6_canesm5(idx) + scen_u_frac_cmip6_canesm5(idx); % model u plus scen u
  ytop      = time(idx)*0+100;
  h2 = patch([x fliplr(x)],[y0 fliplr(ym_cmip6)],[53 74 161]/255,'LineWidth',.5,'Edgecolor','none')
  h1 = patch([x fliplr(x)],[yms_cmip6 fliplr(ytop)],[255 110 4]/255,'LineWidth',.5,'Edgecolor','none')
  h3 = patch([x fliplr(x)],[ym_cmip6 fliplr(yms_cmip6)],[0 127 60]/255,'LineWidth',.5,'Edgecolor','none')
  patch([x fliplr(x)],[ym_cmip6 fliplr(yms_cmip6)],[0 127 60]/255,'LineWidth',.5)%,'Edgecolor','none')
  hold on
  set(gca,'Layer','top','XLim',xlim,'YLim',[0 100],'TickLength',[tl tl])
  h5 = plot(x,ym_cmip6_canesm5,'k:','LineWidth',2)
  plot(x,yms_cmip6_canesm5,'k:','LineWidth',2)
  hp = patch([x fliplr(x)],[ym_cmip6_canesm5 fliplr(yms_cmip6_canesm5)],[1 1 1],'LineWidth',.5,'Edgecolor','none','FaceAlpha',.1);
  hatchfill(hp,'single',49,8,[0 127 60]/255);
  if strcmp(region,'global')==1 & strcmp(seas,'annual')==1 % & strcmp(vari,'tas')==1
    % legend([h1(1) h3(1) h2(1)],'Int. variability','Scenario','Model','Location','SouthWest','FontSize',8)
    legend([h1(1) h3(1) h2(1) h5(1)],'Int. variability','Scenario','Model','Scenario (CanESM5-LE)','Location','SouthWest','FontSize',8)
  end
  % ylabel({'Fractional contribution','to total uncertainty (%)'})
  xlabel('Time (Year)')
  box on
  text(2085,.1*100,'(c)')

  tightfig

  % return
  set(gcf,'PaperPositionMode','auto');
  set(gcf,'renderer','Painters')
  fileo = [pathout_fig vari '/hawkins_plots_3x3panels_' region '_' vari '_' seas '_' num2str(wl) 'yr_ref' num2str(refstart) '-' num2str(refende) '_additional_scen_u'];
  print('-r300','-loose', '-depsc', ['' fileo '.eps'])
  save2pdf(['' fileo '.pdf'])
  saveas(gcf,fileo,'jpg')

end








% ----------------------------------
if plot1 == 1;
% -- plot --------
close all


  figure1 = figure;
  set(figure1, 'units', 'centimeters', 'pos', [10 10 25 20])

  % -- SMILEs:
  subplot(3,3,1)
  subaxis(3,3,1, 'sh', sh, 'sv', sv)
  hold on
  % title([region ' ' vari ' ' seas ' ' num2str(wl) '-yr means (' units ')' char(10) ['relative to ' num2str(refstart) '-' num2str(refende)]],'Interpreter','none')%,'FontSize',10)
  title([region ' decadal mean ' seas char(10) var_name2 ' (' units ') ' ['relative to ' num2str(refstart) '-' num2str(refende)]],'Interpreter','none')%,'FontSize',10)
  plot(time,le_ts_em','Color',rcp_cols_light(1,:))
  h1 = plot(time,nanmean(le_ts_em),'Color',rcp_cols(1,:),'LineWidth',3)
  if strcmp(region,'global')==1 & strcmp(seas,'annual')==1 % & strcmp(vari,'tas')==1
    h0 = plot(time_obs,obs,'k','LineWidth',2)
  end
  set(gca,'YLim',ylim,'YTick',[ceil(ylim(1)):yincr:round(ylim(end))],'XLim',xlim0,'TickLength',[tl tl])
  hline(0,'k')
  box on
  % ylabel([var_name ' change (' units ') relative to ' num2str(refstart) '-' num2str(refende)])
  % ylabel(['Change (' units ') relative ' num2str(refstart) '-' num2str(refende)])
  ylabel('SMILEs (7)','FontWeight','bold')
  if strcmp(region,'global')==1 & strcmp(seas,'annual')==1 % & strcmp(vari,'tas')==1
    legend([h1 h0],'RCP8.5',['Observations' char(10) '(' obs_name ')'],'Location','NorthWest')
  else
    legend([h1],'RCP8.5','Location','NorthWest')
  end
  legend boxoff
  text(2080,ylim(1)+abs(.1*(ylim(2)-ylim(1))),'(a)')

  subplot(3,3,2)
  subaxis(3,3,2, 'sh', sh, 'sv', sv)
  hold on
  % title(['Sources of uncertainty for' char(10) region ' ' vari ' ' seas ' ' num2str(wl) '-yr means (' units ')'],'Interpreter','none')%,'FontSize',10)
  title(['Sources of uncertainty (' units ')'],'Interpreter','none')%,'FontSize',10)
  % % -- as in HS11 Fig. 2
  % -- new way
  tmp   = nanmean(le_ts_em);
  idx   = ~isnan(tmp);
  i     = sqrt(int_u_le_mean(idx));
  imin  = sqrt(int_u_le_min(idx));
  imax  = sqrt(int_u_le_max(idx));
  m     = sqrt(model_u_le(idx));
  s     = sqrt(scen_u_cmip5(idx));
  h3 = patch([time(idx) fliplr(time(idx))],[tmp(idx)-i-m fliplr(tmp(idx)+i+m)],hs09_cols(2,:),'LineWidth',.1,'Edgecolor','none')
  h2 = patch([time(idx) fliplr(time(idx))],[tmp(idx)-m fliplr(tmp(idx)+m)],hs09_cols(1,:),'LineWidth',.1,'Edgecolor','none')
  hold on
  set(gca,'YLim',ylim,'YTick',[ceil(ylim(1)):yincr:round(ylim(end))],'XLim',xlim0,'TickLength',[tl tl],'Layer','top')
  hline(0,'k')
  box on
  % ylabel(['Change (' units ') relative to ' num2str(refstart) '-' num2str(refende)])
  % legend([h3(1) h2(1) h22(1)],'Internal variability','Model',['Int. var. range'],'Location','NorthWest','Interpreter','none')
  legend([h3(1) h2(1)],'Internal variability','Model','Location','NorthWest','Interpreter','none')
  legend boxoff
  text(2080,ylim(1)+abs(.1*(ylim(2)-ylim(1))),'(d)')

  subplot(3,3,3)
  subaxis(3,3,3, 'sh', sh, 'sv', sv)
  hold on
  title({'Fractional contribution','to total uncertainty (%)'})%,'FontSize',10)
  idx = find(~isnan(total_u_le_mean));
  x = time(idx);
  y0        = x*0;
  ym_le     = model_u_frac_le(idx); % model u
  yms_le    = model_u_frac_le(idx) + scen_u_frac_le(idx); % model u plus scen u
  ym_le_hs  = model_u_frac_le_hs(idx); % model u
  yms_le_hs = model_u_frac_le_hs(idx) + scen_u_frac_le_hs(idx); % model u plus scen u
  ytop      = time(idx)*0+100;
  le_fixed  = model_u_frac_le_fixed(idx) + scen_u_frac_le_fixed(idx);
  le_max1   = model_u_frac_le_max(idx);
  le_min1   = model_u_frac_le_min(idx);
  le_max2   = model_u_frac_le_max(idx) + scen_u_frac_le_max(idx);
  le_min2   = model_u_frac_le_min(idx) + scen_u_frac_le_min(idx);
  ym_cmip5  = model_u_frac_cmip5(idx); % model u
  yms_cmip5 = model_u_frac_cmip5(idx) + scen_u_frac_cmip5(idx); % model u plus scen u
  h2 = patch([x fliplr(x)],[y0 fliplr(ym_le)],[53 74 161]/255,'LineWidth',.5,'Edgecolor','none')
  h1 = patch([x fliplr(x)],[yms_le fliplr(ytop)],[255 110 4]/255,'LineWidth',.5,'Edgecolor','none')
  h3 = patch([x fliplr(x)],[ym_le fliplr(yms_le)],[0 127 60]/255,'LineWidth',.5,'Edgecolor','none');
  patch([x fliplr(x)],[ym_le fliplr(yms_le)],[0 127 60]/255,'LineWidth',.5)%,'Edgecolor','none')
  hold on
  set(gca,'Layer','top','XLim',xlim,'YLim',[0 100],'TickLength',[tl tl])
  % jbfill((refstart:refende),(refstart:refende)*0+100,(refstart:refende)*0,[1 1 1],'none',1,.8)
  % hold on
  % if strcmp(region,'global')==1 & strcmp(seas,'annual')==1 % & strcmp(vari,'tas')==1
  %   legend([h1(1) h3(1) h2(1) h4(1)],'Int. variability','Scenario (CMIP5)','Model','Int. var. range','Int. var. fixed','SMILEs HS09','CMIP5 HS09','Location','SouthWest','FontSize',8)
  % end
  if strcmp(region,'global')==1 & strcmp(seas,'annual')==1 % & strcmp(vari,'tas')==1
    legend([h1(1) h3(1) h2(1)],'Int. variability','Scenario (CMIP5)','Model','Location','SouthWest','FontSize',8)
  end
  % ylabel({'Fractional contribution','to total uncertainty (%)'})
  % xlabel('Time (Year)')
  box on
  text(2085,.1*100,'(g)')


  % -- CMIP5:
  subplot(3,3,4)
  % subaxis(3,3,4, 'sh', 0.03, 'sv', 0.01, 'padding', 0, 'margin', 0);
  subaxis(3,3,4, 'sh', sh, 'sv', sv)
  hold on
  % title([region ' ' vari ' ' seas ' ' num2str(wl) '-yr means'],'Interpreter','none')%,'FontSize',10)
  plot(time,cmip5_rcp85_ts_em','Color',rcp_cols_light(1,:))
  plot(time,cmip5_rcp45_ts_em','Color',rcp_cols_light(3,:))
  plot(time,cmip5_rcp26_ts_em','Color',rcp_cols_light(4,:))
  h3 = plot(time,nanmean(cmip5_rcp85_ts_em),'Color',rcp_cols(1,:),'LineWidth',3)
  h2 = plot(time,nanmean(cmip5_rcp45_ts_em),'Color',rcp_cols(3,:),'LineWidth',3)
  h1 = plot(time,nanmean(cmip5_rcp26_ts_em),'Color',rcp_cols(4,:),'LineWidth',3)
  if strcmp(region,'global')==1 & strcmp(seas,'annual')==1 % & strcmp(vari,'tas')==1
    h0 = plot(time_obs,obs,'k','LineWidth',2)
  end
  set(gca,'YLim',ylim,'YTick',[ceil(ylim(1)):yincr:round(ylim(end))],'XLim',xlim0,'TickLength',[tl tl],'Layer','top')
  hline(0,'k')
  box on
  % ylabel([var_name ' change (' units ') relative to ' num2str(refstart) '-' num2str(refende)])
  % ylabel(['Change (' units ') relative ' num2str(refstart) '-' num2str(refende)])
  % if strcmp(region,'global')==1 & strcmp(seas,'annual')==1 % & strcmp(vari,'tas')==1
    legend([h3 h2 h1],'RCP8.5','RCP4.5','RCP2.6','Location','NorthWest')
  % else
  %   legend([h3 h2 h1 h0],'RCP8.5','RCP4.5','RCP2.6',['Observations' char(10) '(' obs_name ')'],'Location','NorthWest')
  % end
  legend boxoff
  ylabel(['CMIP5 (' num2str(length(models_cmip5)) ')'],'FontWeight','bold')
  text(2080,ylim(1)+abs(.1*(ylim(2)-ylim(1))),'(b)')

  subplot(3,3,5)
  subaxis(3,3,5, 'sh', sh, 'sv', sv)
  hold on
  % title('Sources of uncertainty')
  % % -- as in HS11 Fig. 2
  % -- new way
  if constrain == 0
    tmp = nanmean([nanmean(cmip5_rcp85_ts_em,1); nanmean(cmip5_rcp45_ts_em,1); nanmean(cmip5_rcp26_ts_em,1)]);
    i = sqrt(int_u_cmip5_mean);
    % m = sqrt(model_u_cmip5_hs(idx));
    m = sqrt(model_u_cmip5(idx));
    s = sqrt(scen_u_cmip5(idx));
  else
    tmp = nanmean([nanmean(weights_cmip5.*cmip5_rcp85_ts_em,1); nanmean(weights_cmip5.*cmip5_rcp45_ts_em,1); nanmean(weights_cmip5.*cmip5_rcp26_ts_em,1)]);
    i = sqrt(int_u_cmip5_w_mean);
    % m = sqrt(model_u_cmip5_hs(idx));
    m = sqrt(model_u_cmip5_w(idx));
    s = sqrt(scen_u_cmip5_w(idx));
  end
  h3 = patch([time(idx) fliplr(time(idx))],[tmp(idx)-(s+m+i) fliplr(tmp(idx)+(s+m+i))],hs09_cols(2,:),'LineWidth',.1,'Edgecolor','none')
  h2 = patch([time(idx) fliplr(time(idx))],[tmp(idx)-(s+m) fliplr(tmp(idx)+(s+m))],hs09_cols(1,:),'LineWidth',.1,'Edgecolor','none')
  h1 = patch([time(idx) fliplr(time(idx))],[tmp(idx)-s fliplr(tmp(idx)+s)],hs09_cols(3,:),'LineWidth',.1,'Edgecolor','none')
  hold on
  set(gca,'YLim',ylim,'YTick',[ceil(ylim(1)):yincr:round(ylim(end))],'XLim',xlim0,'TickLength',[tl tl],'Layer','top')
  hline(0,'k')
  box on
  % ylabel(['Change (' units ') relative to ' num2str(refstart) '-' num2str(refende)])
  legend([h3(1) h2(1) h1(1)],'Internal variability','Model','Scenario','Location','NorthWest','Interpreter','none')
  legend boxoff
  text(2080,ylim(1)+abs(.1*(ylim(2)-ylim(1))),'(e)')

  subplot(3,3,6)
  subaxis(3,3,6, 'sh', sh, 'sv', sv)
  hold on
  % title({'Fractional contribution','to total uncertainty (%)'})%,'FontSize',10)
  idx       = find(~isnan(total_u_cmip5_mean));
  x         = time(idx);
  y0        = x*0;
  % ym_cmip5  = model_u_frac_cmip5_hs(idx); % model u
  % yms_cmip5 = model_u_frac_cmip5_hs(idx) + scen_u_frac_cmip5_hs(idx); % model u plus scen u
  if constrain == 0
    % ym_cmip5  = model_u_frac_cmip5(idx); % model u
    % yms_cmip5 = model_u_frac_cmip5(idx) + scen_u_frac_cmip5(idx); % model u plus scen u
    ym_cmip5  = model_u_frac_cmip5_brekke(idx); % model u
    yms_cmip5 = model_u_frac_cmip5_brekke(idx) + scen_u_frac_cmip5_brekke(idx); % model u plus scen u
  else
    ym_cmip5  = model_u_frac_cmip5_w(idx); % model u
    yms_cmip5 = model_u_frac_cmip5_w(idx) + scen_u_frac_cmip5_w(idx); % model u plus scen u
  end
  ytop      = time(idx)*0+100;
  h2 = patch([x fliplr(x)],[y0 fliplr(ym_cmip5)],[53 74 161]/255,'LineWidth',.5,'Edgecolor','none')
  h1 = patch([x fliplr(x)],[yms_cmip5 fliplr(ytop)],[255 110 4]/255,'LineWidth',.5,'Edgecolor','none')
  h3 = patch([x fliplr(x)],[ym_cmip5 fliplr(yms_cmip5)],[0 127 60]/255,'LineWidth',.5,'Edgecolor','none')
  patch([x fliplr(x)],[ym_cmip5 fliplr(yms_cmip5)],[0 127 60]/255,'LineWidth',.5)%,'Edgecolor','none')
  hold on
  set(gca,'Layer','top','XLim',xlim,'YLim',[0 100],'TickLength',[tl tl])
  % jbfill((refstart:refende),(refstart:refende)*0+100,(refstart:refende)*0,[1 1 1],'none',1,.8)
  % hold on
  if strcmp(region,'global')==1 & strcmp(seas,'annual')==1 % & strcmp(vari,'tas')==1
    legend([h1(1) h3(1) h2(1)],'Int. variability','Scenario','Model','Location','SouthWest','FontSize',8)
  end
  % ylabel({'Fractional contribution','to total uncertainty (%)'})
  % xlabel('Time (Year)')
  box on
  text(2085,.1*100,'(h)')


  % -- CMIP6:
  subplot(3,3,7)
  subaxis(3,3,7, 'sh', sh, 'sv', sv)
  hold on
  % title([region ' ' vari ' ' seas ' ' num2str(wl) '-yr means'],'Interpreter','none')%,'FontSize',10)
  plot(time,cmip6_ssp585_ts_em','Color',rcp_cols_light(1,:))
  plot(time,cmip6_ssp370_ts_em','Color',rcp_cols_light(2,:))
  plot(time,cmip6_ssp245_ts_em','Color',rcp_cols_light(3,:))
  plot(time,cmip6_ssp126_ts_em','Color',rcp_cols_light(4,:))
  h4 = plot(time,nanmean(cmip6_ssp585_ts_em),'Color',rcp_cols(1,:),'LineWidth',3)
  h3 = plot(time,nanmean(cmip6_ssp370_ts_em),'Color',rcp_cols(2,:),'LineWidth',3)
  h2 = plot(time,nanmean(cmip6_ssp245_ts_em),'Color',rcp_cols(3,:),'LineWidth',3)
  h1 = plot(time,nanmean(cmip6_ssp126_ts_em),'Color',rcp_cols(4,:),'LineWidth',3)
  if strcmp(region,'global')==1 & strcmp(seas,'annual')==1 % & strcmp(vari,'tas')==1
    h0 = plot(time_obs,obs,'k','LineWidth',2)
  end
  set(gca,'YLim',ylim,'YTick',[ceil(ylim(1)):yincr:round(ylim(end))],'XLim',xlim0,'TickLength',[tl tl],'Layer','top')
  hline(0,'k')
  box on
  % ylabel([var_name ' change (' units ') relative to ' num2str(refstart) '-' num2str(refende)])
  % ylabel(['Change (' units ') relative ' num2str(refstart) '-' num2str(refende)])
  if strcmp(region,'global')==1 & strcmp(seas,'annual')==1 % & strcmp(vari,'tas')==1
    legend([h4 h3 h2 h1 h0],'SSP5-8.5','SSP3-7.0','SSP2-4.5','SSP1-2.6',['Observations' char(10) '(' obs_name ')'],'Location','NorthWest')
  else
    legend([h4 h3 h2 h1],'SSP5-8.5','SSP3-7.0','SSP2-4.5','SSP1-2.6','Location','NorthWest')
  end
  legend boxoff
  ylabel(['CMIP6 (' num2str(length(models_cmip6)) ')'],'FontWeight','bold')
  xlabel('Time (Year)')
  text(2080,ylim(1)+abs(.1*(ylim(2)-ylim(1))),'(c)')

  subplot(3,3,8)
  subaxis(3,3,8, 'sh', sh, 'sv', sv)
  hold on
  % title('Sources of uncertainty')
  % -- new way
  if constrain == 0
    tmp = nanmean([nanmean(cmip6_ssp585_ts_em,1); nanmean(cmip6_ssp370_ts_em,1); nanmean(cmip6_ssp245_ts_em,1); nanmean(cmip6_ssp126_ts_em,1)]);
    i = sqrt(int_u_cmip6_mean);
    % m = sqrt(model_u_cmip6_hs(idx));
    m = sqrt(model_u_cmip6(idx));
    s = sqrt(scen_u_cmip6(idx));
  else
    tmp = nanmean([nanmean(weights_cmip6.*cmip6_ssp585_ts_em,1); nanmean(weights_cmip6.*cmip6_ssp370_ts_em,1); nanmean(weights_cmip6.*cmip6_ssp245_ts_em,1); nanmean(weights_cmip6.*cmip6_ssp126_ts_em,1)]);
    i = sqrt(int_u_cmip6_w_mean);
    m = sqrt(model_u_cmip6_w(idx));
    s = sqrt(scen_u_cmip6_w(idx));
  end
  h3 = patch([time(idx) fliplr(time(idx))],[tmp(idx)-(i+m+s) fliplr(tmp(idx)+(i+m+s))],hs09_cols(2,:),'LineWidth',.1,'Edgecolor','none')
  h2 = patch([time(idx) fliplr(time(idx))],[tmp(idx)-(s+m) fliplr(tmp(idx)+(s+m))],hs09_cols(1,:),'LineWidth',.1,'Edgecolor','none')
  h1 = patch([time(idx) fliplr(time(idx))],[tmp(idx)-(s) fliplr(tmp(idx)+(s))],hs09_cols(3,:),'LineWidth',.1,'Edgecolor','none')
  hold on
  set(gca,'YLim',ylim,'YTick',[ceil(ylim(1)):yincr:round(ylim(end))],'XLim',xlim0,'TickLength',[tl tl],'Layer','top')
  hline(0,'k')
  box on
  % ylabel(['Change (' units ') relative to ' num2str(refstart) '-' num2str(refende)])
  legend([h3(1) h2(1) h1(1)],'Internal variability','Model','Scenario','Location','NorthWest','Interpreter','none')
  legend boxoff
  xlabel('Time (Year)')
  text(2080,ylim(1)+abs(.1*(ylim(2)-ylim(1))),'(f)')

  subplot(3,3,9)
  subaxis(3,3,9, 'sh', sh, 'sv', sv)
  hold on
  % title({'Fractional contribution','to total uncertainty (%)'})%,'FontSize',10)
  idx       = find(~isnan(total_u_cmip6_mean));
  x         = time(idx);
  y0        = x*0;
  if constrain == 0
    % ym_cmip6  = model_u_frac_cmip6_hs(idx); % model u
    % yms_cmip6 = model_u_frac_cmip6_hs(idx) + scen_u_frac_cmip6_hs(idx); % model u plus scen u
    ym_cmip6  = model_u_frac_cmip6(idx); % model u
    yms_cmip6 = model_u_frac_cmip6(idx) + scen_u_frac_cmip6(idx); % model u plus scen u
  else
    ym_cmip6  = model_u_frac_cmip6_w(idx); % model u
    yms_cmip6 = model_u_frac_cmip6_w(idx) + scen_u_frac_cmip6_w(idx); % model u plus scen u
  end
  ytop      = time(idx)*0+100;
  h2 = patch([x fliplr(x)],[y0 fliplr(ym_cmip6)],[53 74 161]/255,'LineWidth',.5,'Edgecolor','none')
  h1 = patch([x fliplr(x)],[yms_cmip6 fliplr(ytop)],[255 110 4]/255,'LineWidth',.5,'Edgecolor','none')
  h3 = patch([x fliplr(x)],[ym_cmip6 fliplr(yms_cmip6)],[0 127 60]/255,'LineWidth',.5,'Edgecolor','none')
  patch([x fliplr(x)],[ym_cmip6 fliplr(yms_cmip6)],[0 127 60]/255,'LineWidth',.5)%,'Edgecolor','none')
  hold on
  set(gca,'Layer','top','XLim',xlim,'YLim',[0 100],'TickLength',[tl tl])
  if strcmp(region,'global')==1 & strcmp(seas,'annual')==1 % & strcmp(vari,'tas')==1
    legend([h1(1) h3(1) h2(1)],'Int. variability','Scenario','Model','Location','SouthWest','FontSize',8)
  end
  % ylabel({'Fractional contribution','to total uncertainty (%)'})
  xlabel('Time (Year)')
  box on
  text(2085,.1*100,'(i)')

  tightfig

  return
  set(gcf,'PaperPositionMode','auto');
  set(gcf,'renderer','Painters')
  if constrain == 0
    fileo = [pathout_fig vari '/hawkins_plots_3x3panels_' region '_' vari '_' seas '_' num2str(wl) 'yr_ref' num2str(refstart) '-' num2str(refende) '_reversed_uncertainties'];
  else
    fileo = [pathout_fig vari '/hawkins_plots_3x3panels_' region '_' vari '_' seas '_' num2str(wl) 'yr_ref' num2str(refstart) '-' num2str(refende) '_reversed_uncertainties_constrained'];
  end
  print('-r300','-loose', '-depsc', ['' fileo '.eps'])
  save2pdf(['' fileo '.pdf'])
  saveas(gcf,fileo,'jpg')

end








% ----------------------------------
if plotS3 == 1;
% -- plot --------
close all


  figure1 = figure;
  set(figure1, 'units', 'centimeters', 'pos', [10 10 30 8]) % on home screen

  subplot(1,3,1)
  subaxis(1,3,1, 'sh', sh, 'sv', sv)
  hold on
  title('SMILEs')
  idx = find(~isnan(total_u_le_mean));
  x = time(idx);
  y0        = x*0;
  ym_le     = model_u_frac_le(idx); % model u
  yms_le    = model_u_frac_le(idx) + scen_u_frac_le(idx); % model u plus scen u
  ym_le2    = model_u_frac_le_brekke(idx); % model u
  yms_le2   = model_u_frac_le_brekke(idx) + scen_u_frac_le_brekke(idx); % model u plus scen u
  ytop      = time(idx)*0+100;
  % le_fixed  = model_u_frac_le_fixed(idx) + scen_u_frac_le_fixed(idx);
  % le_max1   = model_u_frac_le_max(idx);
  % le_min1   = model_u_frac_le_min(idx);
  % le_max2   = model_u_frac_le_max(idx) + scen_u_frac_le_max(idx);
  % le_min2   = model_u_frac_le_min(idx) + scen_u_frac_le_min(idx);
  h2 = patch([x fliplr(x)],[y0 fliplr(ym_le)],[53 74 161]/255,'LineWidth',.1,'Edgecolor','none')
  h1 = patch([x fliplr(x)],[yms_le fliplr(ytop)],[255 110 4]/255,'LineWidth',.1,'Edgecolor','none')
  h3 = patch([x fliplr(x)],[ym_le fliplr(yms_le)],[0 127 60]/255,'LineWidth',.1,'Edgecolor','none');
  patch([x fliplr(x)],[ym_le fliplr(yms_le)],[0 127 60]/255,'LineWidth',.1)%,'Edgecolor','none')
  h4 = plot(x,[ym_le2; yms_le2],'k:','LineWidth',2)
  hp = patch([x fliplr(x)],[ym_le2 fliplr(yms_le2)],[1 1 1],'LineWidth',.5,'Edgecolor','none','FaceAlpha',.1);
  hatchfill(hp,'single',65,5,[0 127 60]/255);
  hold on
  set(gca,'Layer','top','XLim',xlim,'YLim',[0 100],'TickLength',[tl tl])
  % if strcmp(region,'global')==1 & strcmp(seas,'annual')==1 % & strcmp(vari,'tas')==1
  legend([h1(1) h3(1) h2(1) h4(1)],'Int. variability','Scenario (CMIP5; HS09)','Model', 'Scenario (CMIP5; Brekke&Barsugli)','Location','SouthWest','FontSize',8)
  % end
  ylabel({'Fractional contribution','to total uncertainty (%)'})
  xlabel('Time (Year)')
  box on
  if strcmp(region,'global')==1 & strcmp(seas,'annual')==1 % & strcmp(vari,'tas')==1
    text(2085,.1*100,'(a)')
  else
    text(2085,.1*100,'(d)')
  end


  % -- CMIP5:
  subplot(1,3,2)
  subaxis(1,3,2, 'sh', sh, 'sv', sv)
  hold on
  title('CMIP5')
  idx       = find(~isnan(total_u_cmip5_mean));
  x         = time(idx);
  y0        = x*0;
  ym_cmip5  = model_u_frac_cmip5(idx); % model u
  yms_cmip5 = model_u_frac_cmip5(idx) + scen_u_frac_cmip5(idx); % model u plus scen u
  ym_cmip52 = model_u_frac_cmip5_brekke(idx); % model u
  yms_cmip52= model_u_frac_cmip5_brekke(idx) + scen_u_frac_cmip5_brekke(idx); % model u plus scen u
  ytop      = time(idx)*0+100;
  h2 = patch([x fliplr(x)],[y0 fliplr(ym_cmip5)],[53 74 161]/255,'LineWidth',.1,'Edgecolor','none')
  h1 = patch([x fliplr(x)],[yms_cmip5 fliplr(ytop)],[255 110 4]/255,'LineWidth',.1,'Edgecolor','none')
  h3 = patch([x fliplr(x)],[ym_cmip5 fliplr(yms_cmip5)],[0 127 60]/255,'LineWidth',.1,'Edgecolor','none')
  patch([x fliplr(x)],[ym_cmip5 fliplr(yms_cmip5)],[0 127 60]/255,'LineWidth',.1)%,'Edgecolor','none')
  h4 = plot(x,[ym_cmip52; yms_cmip52],'k:','LineWidth',2)
  hp = patch([x fliplr(x)],[ym_cmip52 fliplr(yms_cmip52)],[1 1 1],'LineWidth',.5,'Edgecolor','none','FaceAlpha',.1);
  hatchfill(hp,'single',65,5,[0 127 60]/255);
  hold on
  set(gca,'Layer','top','XLim',xlim,'YLim',[0 100],'TickLength',[tl tl])
  % jbfill((refstart:refende),(refstart:refende)*0+100,(refstart:refende)*0,[1 1 1],'none',1,.8)
  % hold on
  legend([h1(1) h3(1) h2(1) h4(1)],'Int. variability','Scenario (HS09)','Model','Scenario (Brekke&Barsugli)','Location','SouthWest','FontSize',8)
  % ylabel({'Fractional contribution','to total uncertainty (%)'})
  xlabel('Time (Year)')
  box on
  if strcmp(region,'global')==1 & strcmp(seas,'annual')==1 % & strcmp(vari,'tas')==1
    text(2085,.1*100,'(b)')
  else
    text(2085,.1*100,'(e)')
  end


  % -- CMIP6:
  subplot(1,3,3)
  subaxis(1,3,3, 'sh', sh, 'sv', sv)
  hold on
  title('CMIP6')
  idx       = find(~isnan(total_u_cmip6_mean));
  x         = time(idx);
  y0        = x*0;
  ym_cmip6  = model_u_frac_cmip6(idx); % model u
  yms_cmip6 = model_u_frac_cmip6(idx) + scen_u_frac_cmip6(idx); % model u plus scen u
  ym_cmip62 = model_u_frac_cmip6_brekke(idx); % model u
  yms_cmip62= model_u_frac_cmip6_brekke(idx) + scen_u_frac_cmip6_brekke(idx); % model u plus scen u
  ytop      = time(idx)*0+100;
  h2 = patch([x fliplr(x)],[y0 fliplr(ym_cmip6)],[53 74 161]/255,'LineWidth',.1,'Edgecolor','none')
  h1 = patch([x fliplr(x)],[yms_cmip6 fliplr(ytop)],[255 110 4]/255,'LineWidth',.1,'Edgecolor','none')
  h3 = patch([x fliplr(x)],[ym_cmip6 fliplr(yms_cmip6)],[0 127 60]/255,'LineWidth',.1,'Edgecolor','none')
  patch([x fliplr(x)],[ym_cmip6 fliplr(yms_cmip6)],[0 127 60]/255,'LineWidth',.1)%,'Edgecolor','none')
  h4 = plot(x,[ym_cmip62; yms_cmip62],'k:','LineWidth',2)
  hp = patch([x fliplr(x)],[ym_cmip62 fliplr(yms_cmip62)],[1 1 1],'LineWidth',.5,'Edgecolor','none','FaceAlpha',.1);
  hatchfill(hp,'single',65,5,[0 127 60]/255);
  hold on
  set(gca,'Layer','top','XLim',xlim,'YLim',[0 100],'TickLength',[tl tl])
  legend([h1(1) h3(1) h2(1) h4(1)],'Int. variability','Scenario (HS09)','Model','Scenario (Brekke&Barsugli)','Location','SouthWest','FontSize',8)
  % ylabel({'Fractional contribution','to total uncertainty (%)'})
  xlabel('Time (Year)')
  box on
  if strcmp(region,'global')==1 & strcmp(seas,'annual')==1 % & strcmp(vari,'tas')==1
    text(2085,.1*100,'(c)')
  else
    text(2085,.1*100,'(f)')
  end

  tightfig

  % return
  set(gcf,'PaperPositionMode','auto');
  set(gcf,'renderer','Painters')
  fileo = [pathout_fig vari '/hawkins_plots_3x1panels_' region '_' vari '_' seas '_' num2str(wl) 'yr_ref' num2str(refstart) '-' num2str(refende) '_reversed_uncertainties_brekke'];
  % print('-r300','-loose', '-depsc', ['' fileo '.eps'])
  save2pdf(['' fileo '.pdf'])
  % saveas(gcf,fileo,'jpg')
  return

end









% -----------------------------

if plot7 == 1;
  % -- plot 2 --------

  figure1 = figure;
  % set(figure1, 'units', 'centimeters', 'pos', [10 10 10 8]) % office screen
  set(figure1, 'units', 'centimeters', 'pos', [10 10 9 8]) % home screen

  hold on
  % title({'Fractional contribution','to total uncertainty (%)'})%,'FontSize',10)
  title([region ' ' vari ' ' seas],'Interpreter','none')%,'FontSize',10)
  idx = find(~isnan(total_u_le_mean));
  x = time(idx);
  y0        = x*0;
  ym_le     = model_u_frac_le(idx); % model u
  yms_le    = model_u_frac_le(idx) + scen_u_frac_le(idx); % model u plus scen u
  ym_le_hs  = model_u_frac_le_hs(idx); % model u
  yms_le_hs = model_u_frac_le_hs(idx) + scen_u_frac_le_hs(idx); % model u plus scen u
  ytop      = time(idx)*0+100;
  le_fixed  = model_u_frac_le_fixed(idx) + scen_u_frac_le_fixed(idx);
  le_max1   = model_u_frac_le_max(idx);
  le_min1   = model_u_frac_le_min(idx);
  le_max2   = model_u_frac_le_max(idx) + scen_u_frac_le_max(idx);
  le_min2   = model_u_frac_le_min(idx) + scen_u_frac_le_min(idx);
  ym_cmip5  = model_u_frac_cmip5(idx); % model u
  yms_cmip5 = model_u_frac_cmip5(idx) + scen_u_frac_cmip5(idx); % model u plus scen u
  ym_cmip6  = model_u_frac_cmip6(idx); % model u
  yms_cmip6 = model_u_frac_cmip6(idx) + scen_u_frac_cmip6(idx); % model u plus scen u
  h2 = patch([x fliplr(x)],[y0 fliplr(ym_le)],[53 74 161]/255,'LineWidth',1.5,'Edgecolor','none')
  h1 = patch([x fliplr(x)],[yms_le fliplr(ytop)],[255 110 4]/255,'LineWidth',2,'Edgecolor','none')
  h3 = patch([x fliplr(x)],[ym_le fliplr(yms_le)],[0 127 60]/255,'LineWidth',2,'Edgecolor','none')
  patch([x fliplr(x)],[ym_le fliplr(yms_le)],[0 127 60]/255,'LineWidth',2)%,'Edgecolor','none')
  h4 = jbfill(time(idx),le_max1,le_min1,[.9 .9 .9],'none',1,.3)
  h4 = jbfill(time(idx),le_max2,le_min2,[.9 .9 .9],'none',1,.5)
  hold on
  plot(time(idx),yms_le,'k','LineWidth',2)
  h5 = plot(time(idx),le_fixed,'k--','LineWidth',2)
  h6 = plot(time(idx),ym_le_hs,'Color','k','LineWidth',1,'LineStyle',':')
  plot(time(idx),yms_le_hs,'Color','k','LineWidth',1,'LineStyle',':')
  h7 = plot(time(idx),ym_cmip5,'Color','k','LineWidth',1,'LineStyle','-')
  plot(time(idx),yms_cmip5,'Color','k','LineWidth',1,'LineStyle','-')

  set(gca,'Layer','top','XLim',xlim,'YLim',[0 100],'TickLength',[.04 .04])
  if strcmp(region,'global')==1 & strcmp(seas,'annual')==1 & strcmp(vari,'tas')==1
    lgd = legend([h1(1) h3(1) h2(1) h4(1) h5(1) h6(1) h7(1)],'Int. variability','Scenario (CMIP5)','Model','Int. var. range','Int. var. fixed','SMILEs HS09','CMIP5 HS09',['CMIP5 HS09' char(10) 'model range'],'Location','NorthEast','FontSize',9)
  end
  % ylabel({'Fractional contribution','to total uncertainty (%)'})
  box on

  tightfig

  % return
  set(gcf,'PaperPositionMode','auto');
  set(gcf,'renderer','Painters')
  fileo = [pathout_fig vari '/hawkins_plots_incl_scenario_from_' num2str(length(models_cmip5)) 'models_' region '_' vari '_' seas '_' num2str(wl) 'yr_ref' num2str(refstart) '-' num2str(refende) '_reduced'];
  print('-r300','-loose', '-depsc', ['' fileo '.eps'])
  save2pdf(['' fileo '.pdf'])
  saveas(gcf,fileo,'jpg')
end






if plot5 == 1;
% -- plot 5 --------
close all
  figure5 = figure;
  if r == 1 || r == 2 || r == 3
    set(figure5, 'units', 'centimeters', 'pos', [10 10 19 8])
  else
    set(figure5, 'units', 'centimeters', 'pos', [10 10 9 7.8])
  end

  if r == 1 || r == 2 || r == 3
    if r == 3
      ylim = [-12 15];
    end
    subplot(1,2,1)
    subaxis(1,2,1, 'sh', sh)%, 'sv', sv)
    hold on
    % title(['Sources of uncertainty for' char(10) region ' ' vari ' ' seas ' ' num2str(wl) '-yr means (' units ')'],'Interpreter','none')%,'FontSize',10)
    % -- new way
    tmp   = nanmean(le_ts_em);
    idx   = ~isnan(tmp);
    i     = sqrt(int_u_le_mean(idx));
    imin  = sqrt(int_u_le_min(idx));
    imax  = sqrt(int_u_le_max(idx));
    % m     = sqrt(model_u_le(idx));
    m     = sqrt(nanmean(model_u_le_hs_all(:,idx),1));
    mmax  = sqrt(max(model_u_le_hs_all(:,idx)));
    mmin  = sqrt(min(model_u_le_hs_all(:,idx)));
    s     = sqrt(scen_u_cmip5(idx));
    % h1 = patch([time(idx) fliplr(time(idx))],[tmp(idx)-(i+m+s) fliplr(tmp(idx)+(i+m+s))],hs09_cols(3,:),'LineWidth',.1,'Edgecolor','none')
    % h22 = jbfill(time(idx),tmp(idx)+i+mmax,tmp(idx)-i-mmax,[.9 .9 .9],[.4 .4 .4],1,.5)
    % hold on
    % h2 = patch([time(idx) fliplr(time(idx))],[tmp(idx)-(i+m) fliplr(tmp(idx)+(i+m))],hs09_cols(1,:),'LineWidth',.1,'Edgecolor','none')
    h2 = patch([time(idx) fliplr(time(idx))],[tmp(idx)-(i+m) fliplr(tmp(idx)+(i+m))],[235, 82, 235]/255,'LineWidth',.1,'Edgecolor','none')
    h3 = patch([time(idx) fliplr(time(idx))],[tmp(idx)-i fliplr(tmp(idx)+i)],hs09_cols(2,:),'LineWidth',.1,'Edgecolor','none')
    h22 = jbfill(time(idx),tmp(idx)+i+m+(mmax-m),tmp(idx)+i+m+(mmin-m),[.9 .9 .9],[.4 .4 .4],1,.5)
    h22 = jbfill(time(idx),tmp(idx)-i-m+(mmax-m),tmp(idx)-i-m+(mmin-m),[.9 .9 .9],[.4 .4 .4],1,.5)
    hold on
    set(gca,'YLim',ylim,'YTick',[ceil(ylim(1)):yincr:round(ylim(end))],'XLim',xlim0,'TickLength',[tl tl],'Layer','top')
    hline(0,'k')
    box on
    % ylabel(['Change (' units ') relative to ' num2str(refstart) '-' num2str(refende)])
    % legend([h3(1) h2(1) h22(1) h1(1)],'Internal variability','Potential method bias','Range method bias','Scenario','Location','NorthWest','Interpreter','none')
    legend([h3(1) h2(1) h22(1)],'Internal variability','Potential method bias',['Range of potential' char(10) 'method bias'],'Location','NorthWest','Interpreter','none','FontSize',8)
    legend boxoff
    % text(2080,ylim(1)+abs(.1*(ylim(2)-ylim(1))),'(a)')
  end

  if r == 1 || r == 2 || r == 3
    subplot(1,2,2)
    subaxis(1,2,2, 'sh', sh)%, 'sv', sv)
  end
  hold on
  idx = find(~isnan(total_u_le_mean));
  x = time(idx);
  y0        = x*0;
  ym_le     = nanmean(model_u_frac_le_hs_all(:,idx),1); % model u
  yms_le    = nanmean(model_u_frac_le_hs_all(:,idx) + scen_u_frac_le_hs_all(:,idx),1); % model u plus scen u
  ytop      = time(idx)*0+100;
  h2 = patch([x fliplr(x)],[y0 fliplr(ym_le)],[235, 82, 235]/255,'LineWidth',1.5,'Edgecolor','none')
  h1 = patch([x fliplr(x)],[yms_le fliplr(ytop)],[255 110 4]/255,'LineWidth',2,'Edgecolor','none')
  h3 = patch([x fliplr(x)],[ym_le fliplr(yms_le)],[0 127 60]/255,'LineWidth',2,'Edgecolor','none')
  patch([x fliplr(x)],[ym_le fliplr(yms_le)],[0 127 60]/255,'LineWidth',2)%,'Edgecolor','none')
  plot(time(idx),ym_le,'Color',[145, 2, 156]/255,'LineWidth',2)
  h4 = jbfill(time(idx),max(model_u_frac_le_hs_all(:,idx),[],1),min(model_u_frac_le_hs_all(:,idx),[],1),[.9 .9 .9],'none',1,.5)
  hold on
  set(gca,'Layer','top','XLim',xlim,'YLim',[0 100],'TickLength',[.04 .04])
  if strcmp(region,'sahel')==1
    lgd = legend([h1(1) h3(1) h2(1) h4(1)],'Internal variability','Scenario (CMIP5)','Potential method bias',['Range of potential' char(10) 'method bias'],'Location','NorthWest','FontSize',8)
  else
    lgd = legend([h1(1) h3(1) h2(1) h4(1)],'Internal variability','Scenario (CMIP5)','Potential method bias',['Range of potential' char(10) 'method bias'],'Location','NorthEast','FontSize',8)
  end
  box on

  tightfig

  % return
  set(gcf,'PaperPositionMode','auto');
  set(gcf,'renderer','Painters')
  fileo = [pathout_fig vari '/hawkins_plots_incl_scenario_from_1le_only_' region '_' vari '_' seas '_' num2str(wl) 'yr_ref' num2str(refstart) '-' num2str(refende)];
  print('-r300','-loose', '-depsc', ['' fileo '.eps'])
  save2pdf(['' fileo '.pdf'])
  saveas(gcf,fileo,'jpg')
  % return

end




if plotS4 == 1;
% -- plot 5.2 --------
close all
  figure5 = figure;
  if r == 1 || r == 2 || r == 3
    set(figure5, 'units', 'centimeters', 'pos', [10 10 19 8])
  else
    set(figure5, 'units', 'centimeters', 'pos', [10 10 9 7.8])
  end

  if r == 1 || r == 2 || r == 3
    if r == 3
      ylim = [-12 15];
    end
    subplot(1,2,1)
    subaxis(1,2,1, 'sh', sh)%, 'sv', sv)
    hold on
    % title(['Sources of uncertainty for' char(10) region ' ' vari ' ' seas ' ' num2str(wl) '-yr means (' units ')'],'Interpreter','none')%,'FontSize',10)
    % -- new way
    tmp   = nanmean(le_ts_em);
    idx   = ~isnan(tmp);
    i     = sqrt(int_u_le_mean(idx));
    imin  = sqrt(int_u_le_min(idx));
    imax  = sqrt(int_u_le_max(idx));
    % m     = sqrt(model_u_le(idx));
    m     = sqrt(nanmean(model_u_le_mpi_subsample(:,idx),1));
    mmax  = sqrt(max(model_u_le_mpi_subsample(:,idx)));
    mmin  = sqrt(min(model_u_le_mpi_subsample(:,idx)));
    s     = sqrt(scen_u_cmip5(idx));
    h2 = patch([time(idx) fliplr(time(idx))],[tmp(idx)-(i+m) fliplr(tmp(idx)+(i+m))],[235, 82, 235]/255,'LineWidth',.1,'Edgecolor','none')
    h3 = patch([time(idx) fliplr(time(idx))],[tmp(idx)-i fliplr(tmp(idx)+i)],hs09_cols(2,:),'LineWidth',.1,'Edgecolor','none')
    hold on
    set(gca,'YLim',ylim,'YTick',[ceil(ylim(1)):yincr:round(ylim(end))],'XLim',xlim0,'TickLength',[tl tl],'Layer','top')
    hline(0,'k')
    box on
    legend([h3(1) h2(1)],'Internal variability','Potential method bias','Location','NorthWest','Interpreter','none','FontSize',8)
    legend boxoff
  end

  if r == 1 || r == 2 || r == 3
    subplot(1,2,2)
    subaxis(1,2,2, 'sh', sh)%, 'sv', sv)
  end
  hold on
  idx = find(~isnan(total_u_le_mean));
  x = time(idx);
  y0        = x*0;
  ym_le     = nanmean(model_u_frac_le_mpi_subsample(:,idx),1); % model u
  yms_le    = nanmean(model_u_frac_le_mpi_subsample(:,idx) + scen_u_frac_le_mpi_subsample(:,idx),1); % model u plus scen u
  ytop      = time(idx)*0+100;
  h2 = patch([x fliplr(x)],[y0 fliplr(ym_le)],[235, 82, 235]/255,'LineWidth',1.5,'Edgecolor','none')
  h1 = patch([x fliplr(x)],[yms_le fliplr(ytop)],[255 110 4]/255,'LineWidth',2,'Edgecolor','none')
  h3 = patch([x fliplr(x)],[ym_le fliplr(yms_le)],[0 127 60]/255,'LineWidth',2,'Edgecolor','none')
  patch([x fliplr(x)],[ym_le fliplr(yms_le)],[0 127 60]/255,'LineWidth',2)%,'Edgecolor','none')
  plot(time(idx),ym_le,'Color',[145, 2, 156]/255,'LineWidth',2)
  set(gca,'Layer','top','XLim',xlim,'YLim',[0 100],'TickLength',[.04 .04])
  if strcmp(region,'sahel')==1
    lgd = legend([h1(1) h3(1) h2(1)],'Internal variability','Scenario (CMIP5)','Potential method bias','Location','NorthWest','FontSize',8)
  else
    lgd = legend([h1(1) h3(1) h2(1)],'Internal variability','Scenario (CMIP5)','Potential method bias','Location','NorthEast','FontSize',8)
  end
  box on

  tightfig

  % return
  set(gcf,'PaperPositionMode','auto');
  set(gcf,'renderer','Painters')
  fileo = [pathout_fig vari '/hawkins_plots_incl_scenario_from_subsampled_mpi_n' num2str(ensmem_subsample) '_' region '_' vari '_' seas '_' num2str(wl) 'yr_ref' num2str(refstart) '-' num2str(refende)];
  print('-r300','-loose', '-depsc', ['' fileo '.eps'])
  save2pdf(['' fileo '.pdf'])
  saveas(gcf,fileo,'jpg')
  % return

end





if plot8 == 1;
% -- plot --------
close all

  xlim_tasg_le    = [0 min(tasg_le_ts_em(end,:))];
  xlim_tasg_cmip5 = [0 min(tasg_cmip5_rcp85_ts_em(end-wl/2,:))];
  xlim_tasg_cmip6 = [0 min(tasg_cmip6_ssp585_ts_em(end-wl/2,:))];
  xlim_tasg       = [0 max([xlim_tasg_le xlim_tasg_cmip5 xlim_tasg_cmip6])];


  figure1 = figure;
  set(figure1, 'units', 'centimeters', 'pos', [10 10 28 7])

  % -- SMILEs:
  subplot(1,3,1)
  subaxis(1,3,1, 'sh', sh, 'sv', sv)
  hold on
  title([region ' ' vari ' ' seas ' ' num2str(wl) '-yr means (' units ')' char(10) ['relative to ' num2str(refstart) '-' num2str(refende)]],'Interpreter','none')%,'FontSize',10)
  plot(tasg_cmip5_rcp85_ts_em(1:end-ceil(wl/2)+1,:),cmip5_rcp85_ts_em(:,1:end-ceil(wl/2)+1)','Color',[204 204 252]/255)
  plot(tasg_cmip6_ssp585_ts_em(1:end-ceil(wl/2)+1,:),cmip6_ssp585_ts_em(:,1:end-ceil(wl/2)+1)','Color',[252 204 204]/255)
  plot(tasg_le_ts_em,le_ts_em','Color',[.5 .5 .5])
  h3 = plot(nanmean(tasg_cmip6_ssp585_ts_em(1:end-ceil(wl/2)+1,:)'),nanmean(cmip6_ssp585_ts_em(:,1:end-ceil(wl/2)+1)),'Color','r','LineWidth',3)
  h2 = plot(nanmean(tasg_cmip5_rcp85_ts_em(1:end-ceil(wl/2)+1,:)'),nanmean(cmip5_rcp85_ts_em(:,1:end-ceil(wl/2)+1)),'Color','b','LineWidth',3)
  h1 = plot(nanmean(tasg_le_ts_em'),nanmean(le_ts_em),'Color','k','LineWidth',3)
  t1 = tasg_cmip6_ssp585_ts_em(end-ceil(wl/2)+1,:);
  t2 = cmip6_ssp585_ts_em(:,end-ceil(wl/2)+1);
  plot([nanmean(t1) nanmean(t1)],[nanmean(t2) nanmean(t2)],'or','MarkerSize',6)
  t1 = tasg_cmip5_rcp85_ts_em(end-ceil(wl/2)+1,:);
  t2 = cmip5_rcp85_ts_em(:,end-ceil(wl/2)+1);
  plot([nanmean(t1) nanmean(t1)],[nanmean(t2) nanmean(t2)],'ob','MarkerSize',6)
  t1 = tasg_le_ts_em(end-ceil(wl/2)+1,:);
  t2 = le_ts_em(:,end-ceil(wl/2)+1);
  plot([nanmean(t1) nanmean(t1)],[nanmean(t2) nanmean(t2)],'ok','MarkerSize',6)
  set(gca,'YLim',ylim,'YTick',[ceil(ylim(1)):yincr:round(ylim(end))],'TickLength',[tl tl])%,'XLim',xlim_tasg)
  hline(0,'k')
  box on
  xlabel(['Global temperature change' char(10) 'from ' num2str(refstart) '-' num2str(refende) ' (\circC)'])
  legend([h1 h2 h3],'SMILEs','CMIP5','CMIP6',['Observations' char(10) '(' obs_name ')'],'Location','NorthWest','FontSize',8)
  legend boxoff

  subplot(1,3,2)
  subaxis(1,3,2, 'sh', sh, 'sv', sv)
  title(['Absolute uncertainty' char(10) '(standard deviation; ' units ')'])
  hold on
  h22 = plot(range,range*0+int_u_cmip5_tasg_mean.^.5,'b--','LineWidth',2)
  h12 = plot(range,int_u_le_tasg_mean.^.5,'k--','LineWidth',2)
  h32 = plot(range,range*0+int_u_cmip6_tasg_mean.^.5,'r--','LineWidth',2)
  h11 = plot(range,model_u_le_tasg.^.5,'k','LineWidth',2)
  h21 = plot(range,model_u_cmip5_rcp85_tasg.^.5,'b','LineWidth',2)
  h31 = plot(range,model_u_cmip6_ssp585_tasg.^.5,'r','LineWidth',2)
  set(gca,'XLim',xlim_tasg,'TickLength',[tl tl],'Layer','top')
  box on
  xlabel(['Global temperature change' char(10) 'from ' num2str(refstart) '-' num2str(refende) ' (\circC)'])
  % ylabel(['(' units '^2) relative to ' num2str(refstart) '-' num2str(refende)])
  legend([h11 h21 h31 h12 h22 h32],'Model SMILEs','Model CMIP5','Model CMIP6',...
  'Int. var. SMILEs','Int. var. CMIP5','Int. var. CMIP6',...
  'Location','NorthWest','Interpreter','none','FontSize',8)
  legend boxoff

  subplot(1,3,3)
  subaxis(1,3,3, 'sh', sh, 'sv', sv)
  title({'Fractional contribution','to total uncertainty (%)'})%,'FontSize',10)
  hold on
  idx       = find(~isnan(total_u_le_tasg_mean));
  x         = range(idx);
  y0        = x*0;
  ym_le     = model_u_frac_le_tasg(idx); % model u
  ytop      = range(idx)*0+100;
  h2 = patch([x fliplr(x)],[y0 fliplr(ym_le)],[53 74 161]/255,'LineWidth',.5,'Edgecolor','none')
  h1 = patch([x fliplr(x)],[ym_le fliplr(ytop)],[255 110 4]/255,'LineWidth',.5,'Edgecolor','none')
  h3 = plot(x,model_u_frac_cmip5_tasg(idx),'k','LineWidth',2)
  h4 = plot(x,model_u_frac_cmip6_tasg(idx),'k:','LineWidth',2)
  set(gca,'Layer','top','XLim',xlim_tasg,'YLim',[0 100],'TickLength',[tl tl])
  xlabel(['Global temperature change' char(10) 'from ' num2str(refstart) '-' num2str(refende) ' (\circC)'])
  legend([h1(1) h2(1) h3(1) h4(1)],'Int. var. SMILEs','Model SMILEs','Dividing line CMIP5','Dividing line CMIP6','Location','SouthEast','FontSize',8)
  box on

  tightfig

  return
  set(gcf,'PaperPositionMode','auto');
  set(gcf,'renderer','Painters')
  fileo = [pathout_fig vari '/hawkins_plots_vs_tasg_' region '_' vari '_' seas '_' num2str(wl) 'yr_ref' num2str(refstart) '-' num2str(refende)];
  print('-r300','-loose', '-depsc', ['' fileo '.eps'])
  save2pdf(['' fileo '.pdf'])
  saveas(gcf,fileo,'jpg')

end





% ----------------------------------
if plot9 == 1;
% -- plot --------
close all


  figure1 = figure;
  set(figure1, 'units', 'centimeters', 'pos', [10 10 17 15])

  % -- CMIP5:
  subplot(2,2,1)
  subaxis(2,2,1, 'sh', sh, 'sv', sv)
  hold on
  % title([region ' ' vari ' ' seas ' ' num2str(wl) '-yr means (' units ')' char(10) ['relative to ' num2str(refstart) '-' num2str(refende)]],'Interpreter','none')%,'FontSize',10)
  title(['Global annual decadal mean' char(10) 'temperature (' units '), ' ['relative to ' num2str(refstart) '-' num2str(refende)]],'Interpreter','none')%,'FontSize',10)
  % -- new way
  tmp1 = nanmean([nanmean(cmip5_rcp85_ts_em,1); nanmean(cmip5_rcp45_ts_em,1); nanmean(cmip5_rcp26_ts_em,1)]);
  idx   = ~isnan(tmp1);
  i1 = sqrt(int_u_cmip5_mean);
  m1 = sqrt(model_u_cmip5(idx));
  s1 = sqrt(scen_u_cmip5(idx));
  tmp2 = nanmean([nanmean(weights_cmip5.*cmip5_rcp85_ts_em,1); nanmean(weights_cmip5.*cmip5_rcp45_ts_em,1); nanmean(weights_cmip5.*cmip5_rcp26_ts_em,1)]);
  i2 = sqrt(int_u_cmip5_w_mean);
  m2 = sqrt(model_u_cmip5_w(idx));
  s2 = sqrt(scen_u_cmip5_w(idx));
  x = time(idx);
  % h1 = patch([x fliplr(x)],[tmp1(idx)-(i1+m1+s1) fliplr(tmp1(idx)+(i1+m1+s1))],hs09_cols(3,:),'LineWidth',.1,'Edgecolor','none')
  % h2 = patch([x fliplr(x)],[tmp1(idx)-(i1+m1) fliplr(tmp1(idx)+(i1+m1))],hs09_cols(1,:),'LineWidth',.1,'Edgecolor','none')
  % h3 = patch([x fliplr(x)],[tmp1(idx)-i1 fliplr(tmp1(idx)+i1)],hs09_cols(2,:),'LineWidth',.1,'Edgecolor','none')
  % plot(x,tmp2(idx),'k:','LineWidth',2)
  % h4 = plot(x,tmp2(idx),'k:','LineWidth',2)
  % plot(x,[tmp2(idx)-(i2+m2); tmp2(idx)+(i2+m2)],'k:','LineWidth',2)
  % plot(x,[tmp2(idx)-(i2+m2+s2); tmp2(idx)+(i2+m2+s2)],'k:','LineWidth',2)
  h3 = patch([x fliplr(x)],[tmp1(idx)-(i1+m1+s1) fliplr(tmp1(idx)+i1+m1+s1)],hs09_cols(2,:),'LineWidth',.1,'Edgecolor','none')
  h2 = patch([x fliplr(x)],[tmp1(idx)-(m1+s1) fliplr(tmp1(idx)+(m1+s1))],hs09_cols(1,:),'LineWidth',.1,'Edgecolor','none')
  h1 = patch([x fliplr(x)],[tmp1(idx)-(s1) fliplr(tmp1(idx)+(s1))],hs09_cols(3,:),'LineWidth',.1,'Edgecolor','none')
  h4 = plot(x,[tmp2(idx)-(s2); tmp2(idx)+(s2)],'k:','LineWidth',2)
  plot(x,[tmp2(idx)-(m2+s2); tmp2(idx)+(m2+s2)],'k:','LineWidth',2)
  plot(x,[tmp2(idx)-(i2+m2+s2); tmp2(idx)+(i2+m2+s2)],'k:','LineWidth',2)
  hold on
  set(gca,'YLim',ylim,'YTick',[ceil(ylim(1)):yincr:round(ylim(end))],'XLim',xlim0,'TickLength',[tl tl],'Layer','top')
  hline(0,'k')
  box on
  ylabel(['CMIP5 (' num2str(length(models_cmip5)) ')'])
  legend([h3(1) h2(1) h1(1) h4(1)],'Internal variability','Model','Scenario','Constrained','Location','NorthWest','Interpreter','none')
  legend boxoff
  text(2080,ylim(1)+abs(.1*(ylim(2)-ylim(1))),'(a)')

  subplot(2,2,2)
  subaxis(2,2,2, 'sh', sh, 'sv', sv)
  hold on
  title({'Fractional contribution','to total uncertainty (%)'})%,'FontSize',10)
  idx       = find(~isnan(total_u_cmip5_mean));
  x         = time(idx);
  y0        = x*0;
  ym1_cmip5  = model_u_frac_cmip5(idx); % model u
  yms1_cmip5 = model_u_frac_cmip5(idx) + scen_u_frac_cmip5(idx); % model u plus scen u
  ym2_cmip5  = model_u_frac_cmip5_w(idx); % model u
  yms2_cmip5 = model_u_frac_cmip5_w(idx) + scen_u_frac_cmip5_w(idx); % model u plus scen u
  % -- find bootstrap member with highest/lowest summed fraction:
  [tmp,ni]            = max(nansum(model_u_frac_cmip5_w_bs(:,idx),2));
  ym3_cmip5_max       = model_u_frac_cmip5_w_bs(ni,idx); % model u
  yms3_cmip5_max      = model_u_frac_cmip5_w_bs(ni,idx) + scen_u_frac_cmip5_w_bs(ni,idx); % model u plus scen u
  [tmp,ni]            = min(nansum(model_u_frac_cmip5_w_bs(:,idx),2));
  ym3_cmip5_min       = model_u_frac_cmip5_w_bs(ni,idx); % model u
  yms3_cmip5_min      = model_u_frac_cmip5_w_bs(ni,idx) + scen_u_frac_cmip5_w_bs(ni,idx); % model u plus scen u
  ytop      = time(idx)*0+100;
  h2 = patch([x fliplr(x)],[y0 fliplr(ym1_cmip5)],[53 74 161]/255,'LineWidth',.5,'Edgecolor','none')
  h1 = patch([x fliplr(x)],[yms1_cmip5 fliplr(ytop)],[255 110 4]/255,'LineWidth',.5,'Edgecolor','none')
  h3 = patch([x fliplr(x)],[ym1_cmip5 fliplr(yms1_cmip5)],[0 127 60]/255,'LineWidth',.5,'Edgecolor','none')
  patch([x fliplr(x)],[ym1_cmip5 fliplr(yms1_cmip5)],[0 127 60]/255,'LineWidth',.5)%,'Edgecolor','none')
  h4 = plot(x,ym2_cmip5,'k:','LineWidth',2)
  plot(x,yms2_cmip5,'k:','LineWidth',2)
  hold on
  set(gca,'Layer','top','XLim',xlim,'YLim',[0 100],'TickLength',[tl tl])
  if strcmp(region,'global')==1 & strcmp(seas,'annual')==1 % & strcmp(vari,'tas')==1
    legend([h1(1) h3(1) h2(1) h4(1)],'Int. variability','Scenario','Model','Constrained','Location','NorthEast','FontSize',8)
  end
  box on
  text(2085,.1*100,'(b)')


  % -- CMIP6:
  subplot(2,2,3)
  subaxis(2,2,3, 'sh', sh, 'sv', sv)
  hold on
  % -- new way
  tmp1 = nanmean([nanmean(cmip6_ssp585_ts_em,1); nanmean(cmip6_ssp370_ts_em,1); nanmean(cmip6_ssp245_ts_em,1); nanmean(cmip6_ssp126_ts_em,1)]);
  i1 = sqrt(int_u_cmip6_mean);
  m1 = sqrt(model_u_cmip6(idx));
  s1 = sqrt(scen_u_cmip6(idx));
  tmp2 = nanmean([nanmean(weights_cmip6.*cmip6_ssp585_ts_em,1); nanmean(weights_cmip6.*cmip6_ssp370_ts_em,1); nanmean(weights_cmip6.*cmip6_ssp245_ts_em,1); nanmean(weights_cmip6.*cmip6_ssp126_ts_em,1)]);
  i2 = sqrt(int_u_cmip6_w_mean);
  m2 = sqrt(model_u_cmip6_w(idx));
  s2 = sqrt(scen_u_cmip6_w(idx));
  h3 = patch([x fliplr(x)],[tmp1(idx)-(i1+m1+s1) fliplr(tmp1(idx)+i1+m1+s1)],hs09_cols(2,:),'LineWidth',.1,'Edgecolor','none')
  h2 = patch([x fliplr(x)],[tmp1(idx)-(m1+s1) fliplr(tmp1(idx)+(m1+s1))],hs09_cols(1,:),'LineWidth',.1,'Edgecolor','none')
  h1 = patch([x fliplr(x)],[tmp1(idx)-(s1) fliplr(tmp1(idx)+(s1))],hs09_cols(3,:),'LineWidth',.1,'Edgecolor','none')
  h4 = plot(x,[tmp2(idx)-(s2); tmp2(idx)+(s2)],'k:','LineWidth',2)
  plot(x,[tmp2(idx)-(m2+s2); tmp2(idx)+(m2+s2)],'k:','LineWidth',2)
  plot(x,[tmp2(idx)-(i2+m2+s2); tmp2(idx)+(i2+m2+s2)],'k:','LineWidth',2)
  hold on
  set(gca,'YLim',ylim,'YTick',[ceil(ylim(1)):yincr:round(ylim(end))],'XLim',xlim0,'TickLength',[tl tl],'Layer','top')
  hline(0,'k')
  box on
  ylabel(['CMIP6 (' num2str(length(models_cmip6)) ')'])
  xlabel('Time (Year)')
  text(2080,ylim(1)+abs(.1*(ylim(2)-ylim(1))),'(c)')

  subplot(2,2,4)
  subaxis(2,2,4, 'sh', sh, 'sv', sv)
  hold on
  % title({'Fractional contribution','to total uncertainty (%)'})%,'FontSize',10)
  idx       = find(~isnan(total_u_cmip6_mean));
  x         = time(idx);
  y0        = x*0;
  ym1_cmip6  = model_u_frac_cmip6(idx); % model u
  yms1_cmip6 = model_u_frac_cmip6(idx) + scen_u_frac_cmip6(idx); % model u plus scen u
  ym2_cmip6  = model_u_frac_cmip6_w(idx); % model u
  yms2_cmip6 = model_u_frac_cmip6_w(idx) + scen_u_frac_cmip6_w(idx); % model u plus scen u
  % -- find bootstrap member with highest/lowest summed fraction:
  [tmp,ni]            = max(nansum(model_u_frac_cmip6_w_bs(:,idx),2));
  ym3_cmip6_max       = model_u_frac_cmip6_w_bs(ni,idx); % model u
  yms3_cmip6_max      = model_u_frac_cmip6_w_bs(ni,idx) + scen_u_frac_cmip6_w_bs(ni,idx); % model u plus scen u
  [tmp,ni]            = min(nansum(model_u_frac_cmip6_w_bs(:,idx),2));
  ym3_cmip6_min       = model_u_frac_cmip6_w_bs(ni,idx); % model u
  yms3_cmip6_min      = model_u_frac_cmip6_w_bs(ni,idx) + scen_u_frac_cmip6_w_bs(ni,idx); % model u plus scen u

  ytop      = time(idx)*0+100;
  h2 = patch([x fliplr(x)],[y0 fliplr(ym1_cmip6)],[53 74 161]/255,'LineWidth',.5,'Edgecolor','none')
  h1 = patch([x fliplr(x)],[yms1_cmip6 fliplr(ytop)],[255 110 4]/255,'LineWidth',.5,'Edgecolor','none')
  h3 = patch([x fliplr(x)],[ym1_cmip6 fliplr(yms1_cmip6)],[0 127 60]/255,'LineWidth',.5,'Edgecolor','none')
  patch([x fliplr(x)],[ym1_cmip6 fliplr(yms1_cmip6)],[0 127 60]/255,'LineWidth',.5)%,'Edgecolor','none')
  plot(x,ym2_cmip6,'k:','LineWidth',2)
  plot(x,yms2_cmip6,'k:','LineWidth',2)
  plot(x,ym3_cmip6_max,':','Color',[.7 .7 .7],'LineWidth',1)
  plot(x,ym3_cmip6_min,':','Color',[.7 .7 .7],'LineWidth',1)
  hold on
  set(gca,'Layer','top','XLim',xlim,'YLim',[0 100],'TickLength',[tl tl])
  xlabel('Time (Year)')
  box on
  text(2085,.1*100,'(d)')

  tightfig

  return
  set(gcf,'PaperPositionMode','auto');
  set(gcf,'renderer','Painters')
  fileo = [pathout_fig vari '/hawkins_plots_3x3panels_' region '_' vari '_' seas '_' num2str(wl) 'yr_ref' num2str(refstart) '-' num2str(refende) '_constrained'];
  print('-r300','-loose', '-depsc', ['' fileo '.eps'])
  save2pdf(['' fileo '.pdf'])
  saveas(gcf,fileo,'jpg')

end







% ----------------------------------
if plotS1 == 1;
% -- plot --------
close all

  % -- across-member variance from SMILEs:
  tmp01 = [zeros(size(nanmean(le_ts_var,2))) nanmean(le_ts_var(:,1:refende-start+1),2) nanmean(le_ts_var,2) nanmean(le_ts_var(:,refende-start+1+1:end),2)].^.5;
  % -- across-time variance from HS09 residuals from SMILEs:
  tmp02 = [cmip5_piControl_ts_var([models_le_cmip5_id])' nanmean(le_ts_var_hs2,2) nanmean(le_ts_var_hs1,2) nanmean(le_ts_var_hs3,2)].^.5;
  tmp0  = [tmp02 tmp01];

  % -- across-time variance from HS09 residuals from CMIP5
  tmp1  = [cmip5_piControl_ts_var; cmip5_rcp85_ts_var2; cmip5_rcp85_ts_var1; cmip5_rcp85_ts_var3;].^.5;
  % -- across-time variance from HS09 residuals from CMIP6
  tmp2  = [cmip6_piControl_ts_var; cmip6_ssp585_ts_var2; cmip6_ssp585_ts_var1; cmip6_ssp585_ts_var3].^.5;

  % -- observations
  obs_std = nanstd(rm(obs_residual,wl));

  xlim = [0 length(models_cmip5)+3];
  ylim = [0 1.1*max([tmp0(:); tmp1(:); tmp2(:)])];
  cols = [0 115 195;
          212 0 11;
          244 116 0;
          252 210 3;
          0 195 115;
          23 125 0;
          116 244 0;
          210 252 3]/255;

  offset = .05;
  tl1 = .0;
  tl2 = .01;
  pp1 = 10;
  pp2 = 90;

  figure1 = figure;
  set(figure1, 'units', 'centimeters', 'pos', [10 10 35 25])
  hold on

  subplot(3,1,1)
  title(['Variability of ' region ' ' num2str(wl) '-yr mean ' seas ' ' vari ],'Interpreter','none')%,'FontSize',10)
  % title(['Variability of ' region ' decadal mean annual ' var_name2],'Interpreter','none')%,'FontSize',10)
  hold on
  plot([0 11.6],[obs_std obs_std],'k','LineWidth',2)
  text(9.5,double(obs_std),['Observations ' char(10) '1950-2014 (fit residual)'],'FontSize',10)
  b = bar([tmp0]);
  for i = [1:4]
    set(b(i),'FaceColor',cols(i,:))
    plot([i/6 i/6]+7.5,[prctile(tmp0(:,i),pp1) prctile(tmp0(:,i),pp2)],'Color',cols(i,:),'LineWidth',4)
    plot([i/6 i/6]+7.5,[nanmean(tmp0(:,i)) nanmean(tmp0(:,i))],'.','Color',cols(i,:),'MarkerSize',25)
  end
  for i = 6:8
    set(b(i),'FaceColor',cols(i,:))
    plot([i/6 i/6]+7.5,[prctile(tmp0(:,i),pp1) prctile(tmp0(:,i),pp2)],'Color',cols(i,:),'LineWidth',4)
    plot([i/6 i/6]+7.5,[nanmean(tmp0(:,i)) nanmean(tmp0(:,i))],'.','Color',cols(i,:),'MarkerSize',25)
  end
  set(gca,'XLim',xlim/2,'YLim',ylim,'XTick',[1:1:length(models) 7.5+5/6],'XTickLabel',[model_names 'Multi-model'],'FontSize',10,'YGrid','on')
  h = gca;
  h.XRuler.TickLength = [tl1,tl1];
  h.YRuler.TickLength = [tl2,tl2];
  xtickangle(30)
  legend([b(1:4) b(6:8)],'piControl','1950-2014 (fit residual)','1950-2099 (fit residual)','2015-2099 (fit residual)','1950-2014 (across member)','1950-2099 (across member)','2015-2099 (across member)','Location','NorthEast','FontSize',10)
  legend boxoff
  ylabel(['Standard deviation (' units ')'],'FontSize',10)
  text(.5,.85*ylim(2),'(a) SMILEs','FontWeight','bold','FontSize',12)
  box on

  subplot(3,1,2)
  hold on
  plot([0 xlim(2)],[obs_std obs_std],'k','LineWidth',2)
  b = bar([tmp1]');
  for i = 1:4
    set(b(i),'FaceColor',cols(i,:))
    plot([i/3 i/3]+28.5,[prctile(tmp1(i,:),pp1) prctile(tmp1(i,:),pp2)],'Color',cols(i,:),'LineWidth',4)
    plot([i/3 i/3]+28.5,[nanmean(tmp1(i,:)) nanmean(tmp1(i,:))],'.','Color',cols(i,:),'MarkerSize',25)
  end
  set(gca,'XLim',xlim,'YLim',ylim,'XTick',[1:1:length(models_cmip5) 28.5+2.5/3],'XTickLabel',[models_cmip5 'Multi-model'],'FontSize',10,'YGrid','on')
  h = gca;
  h.XRuler.TickLength = [tl1,tl1];
  h.YRuler.TickLength = [tl2,tl2];
  xtickangle(30)
  legend([b(1:4)],'piControl','1950-2014 (fit residual)','1950-2099 (fit residual)','2015-2099 (fit residual)','Location','NorthEast','FontSize',10)
  legend boxoff
  ylabel(['Standard deviation (' units ')'],'FontSize',10)
  text(1,.85*ylim(2),'(b) CMIP5','FontWeight','bold','FontSize',12)
  box on

  subplot(3,1,3)
  hold on
  plot([0 xlim(2)],[obs_std obs_std],'k','LineWidth',2)
  b = bar([tmp2]');
  for i = 1:4
    set(b(i),'FaceColor',cols(i,:))
    plot([i/3 i/3]+21.5,[prctile(tmp2(i,:),pp1) prctile(tmp2(i,:),pp2)],'Color',cols(i,:),'LineWidth',4)
    plot([i/3 i/3]+21.5,[nanmean(tmp2(i,:)) nanmean(tmp2(i,:))],'.','Color',cols(i,:),'MarkerSize',25)
  end
  set(gca,'XLim',xlim,'YLim',ylim,'XTick',[1:1:length(models_cmip6) 21.5+2.5/3],'XTickLabel',[models_cmip6 'Multi-model'],'FontSize',10,'YGrid','on')
  h = gca;
  h.XRuler.TickLength = [tl1,tl1];
  h.YRuler.TickLength = [tl2,tl2];
  xtickangle(30)
  legend([b(1:4)],'piControl','1950-2014 (fit residual)','1950-2099 (fit residual)','2015-2099 (fit residual)','Location','NorthEast','FontSize',10)
  legend boxoff
  ylabel(['Standard deviation (' units ')'],'FontSize',10)
  text(1,.85*ylim(2),'(c) CMIP6','FontWeight','bold','FontSize',12)
  box on

  tightfig

  return
  set(gcf,'PaperPositionMode','auto');
  set(gcf,'renderer','Painters')
  fileo = [pathout_fig vari '/variability_by_model_' region '_' vari '_' seas '_' num2str(wl) 'yr'];
  print('-r300','-loose', '-depsc', ['' fileo '.eps'])
  save2pdf(['' fileo '.pdf'])
  saveas(gcf,fileo,'jpg')
  return

end







if plotS5 == 1;
% -- plot --------
  close all

  m1 = 4;
  m2 = 5;

  figure1 = figure;
  set(figure1, 'units', 'centimeters', 'pos', [10 10 11 8])
  hold on
  h1 = plot(time,cmip5_rcp85_ts_em,'Color',[1 .5 .5])
  h2 = plot(time,le_ts_em,'Color',[.4 .4 .4],'LineWidth',2)
  h3 = plot(time,le_ts_em(m1,:),'Color','b','LineWidth',3)
  h4 = plot(time,le_ts_em(m2,:),'Color','r','LineWidth',3)
  box on
  legend([h1(1) h2(1) h3(1) h4(1)],'CMIP5','SMILEs','GFDL-CM3','GFDL-ESM2M','location','northwest')
  legend boxoff
  set(gca,'XLim',[1960 2090])
  hline(0,'k')
  xlabel('Time (Year)')
  ylabel([ region ' ' vari ' ' seas ' change (' units ')'],'Interpreter','none')

  tightfig

  % return
  set(gcf,'PaperPositionMode','auto');
  set(gcf,'renderer','Painters')
  fileo = [pathout_fig vari '/hawkins_plots_model_independence_' region '_' vari '_' seas];
  print('-r300','-loose', '-depsc', ['' fileo '.eps'])
  % save2pdf(['' fileo '.pdf'])
  % saveas(gcf,fileo,'jpg')
  return

end
