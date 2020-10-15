; load "~asphilli/ncl/asphilli_util.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin

  plot1   = "yes" ; normal maps with fraction of uncertainty
  plot2   = "no" ; maps with varying reference periods
  plot3   = "no" ; maps with difference between "SMILEs" and "SMILEs with HS09 approach"
  plot4   = "no" ; maps with individual SMILEs masquerading as different models to estimate a moch "model uncertainty" --> method error

  constrain = "no" ; yes=only use constrained CMIP6 ensemble (Tokarska et al)

  vars    = (/"tas","pr"/)    ; tas or pr
  seasons = (/"DJF","JJA","annual"/)    ; DJF JJA annual
  comp    = "Amon"
  domain  = "global"    ; global or NA
  ; pathin  = "/project/yampa03/clivar_wg_le/"
  pathin  = "/net/meso/climphys/flehner/"
  pathout = "/net/meso/climphys/flehner/"
  start_trend = 1950
  ende_trend  = 2099 ; MPI ends in 2099
  refstart    = 1995; 1986 1990 1971 1991
  refende     = 2014; 2014 2019 2000 2020
  trend_time  = fspan(start_trend,ende_trend,ende_trend-start_trend+1)
  ntime       = dimsizes(trend_time)
  models      = (/"canesm2","cesm","csiro_mk36","ec_earth","gfdl_cm3","gfdl_esm2m","mpi"/)
  model_names = (/"CanESM2","CESM1-CAM5","CSIRO-Mk-3-6-0","EC-EARTH","GFDL-CM3","GFDL-ESM2M","MPI"/)
  start0      = (/1950,1920,1850,1860,1920,1950,1850/)
  ende0       = (/2100,2100,2100,2100,2100,2100,2099/)
  ; ensmem0     = (/50,40,30,16,20,30,100/)
  ensmem0     = (/50,35,30,16,20,30,100/)

  ; --> regridding of CLIVAR LEs with ~/scripts/clivar_le_atm_regrid.sh and cmip5 ad hoc with
  ;     cdo remapcon2,/project/cas/flehner/gridfile_lamber_equal_area.txt filein.nc fileout.nc
  grid    = "g025" ; g025 gea(equal-area)

  ; -- method for model uncertainty --
  mu      = 1 ; 1=just use rcp85 or ssp585;  2=hs09 (mean across scenarios)

  ; -- recalc ensemble statistics? --
  recalc1 = "no" ; LEs
  recalc2 = "no" ; cmip5
  recalc3 = "no" ; cmip6

  ; -- regrid simulations? (only in case of gea [18x10]) --
  regrid1 = "no" ; LEs
  regrid2 = "no" ; cmip5
  regrid3 = "no" ; cmip6

  ; -- CMIP5 --
  ; -- consistent set of models for rcp26, rcp45, and rcp85:
  cmip5_models = (/"BNU-ESM","CanESM2","CCSM4","CESM1-CAM5","CNRM-CM5",\
                   "CSIRO-Mk3-6-0","CanESM2","FGOALS-g2","FIO-ESM","GFDL-CM3",\
                   "GFDL-ESM2G","GFDL-ESM2M","GISS-E2-H","GISS-E2-R","HadGEM2-AO",\
                   "HadGEM2-ES","IPSL-CM5A-LR","IPSL-CM5A-MR","MIROC-ESM-CHEM","MIROC-ESM",\
                   "MIROC5","MPI-ESM-LR","MPI-ESM-MR","MRI-CGCM3","NorESM1-ME",\
                   "NorESM1-M","bcc-csm1-1-m","bcc-csm1-1"/)

  ; -- CMIP6 --
  if constrain .eq. "no"
    ending = ""
    ; ; -- consistent set of models for ssp126, ssp245, and ssp585 (status 12/04/2019):
    ; cmip6_models = (/"BCC-CSM2-MR","CAMS-CSM1-0","CESM2","CESM2-WACCM","CNRM-CM6-1",\
    ;                 "CNRM-ESM2-1","CanESM5","EC-Earth3","EC-Earth3-Veg","FGOALS-f3-L",\
    ;                 "FGOALS-g3","GFDL-ESM4","INM-CM4-8","INM-CM5-0","IPSL-CM6A-LR",\
    ;                 "MCM-UA-1-0","MIROC-ES2L","MIROC6","MPI-ESM1-2-HR","MRI-ESM2-0",\
    ;                 "NESM3","UKESM1-0-LL"/)
    ; -- consistent set of models for ssp126, ssp245, ssp370, and ssp585 (status 12/12/2019):
    cmip6_models = (/"BCC-CSM2-MR","CAMS-CSM1-0","CESM2","CESM2-WACCM","CNRM-CM6-1",\
                    "CNRM-ESM2-1","CanESM5","EC-Earth3","EC-Earth3-Veg","FGOALS-f3-L",\
                    "FGOALS-g3","GFDL-ESM4","INM-CM4-8","INM-CM5-0","IPSL-CM6A-LR",\
                    "MCM-UA-1-0","MIROC-ES2L","MIROC6","MPI-ESM1-2-HR","MRI-ESM2-0",\
                    "UKESM1-0-LL"/)
  else
    ending = "_wo_highsens_models"
    ; -- same as above but WITHOUT high-sensitivity models (according to Kasia):
    cmip6_models = (/"BCC-CSM2-MR","CAMS-CSM1-0","CESM2","CESM2-WACCM","CNRM-CM6-1",\
                    "CNRM-ESM2-1","EC-Earth3","FGOALS-f3-L",\
                    "FGOALS-g3","GFDL-ESM4","INM-CM4-8","INM-CM5-0","IPSL-CM6A-LR",\
                    "MCM-UA-1-0","MIROC-ES2L","MIROC6","MPI-ESM1-2-HR","MRI-ESM2-0"/)
  end if

  if grid .eq. "g025"
    nlat = 72
    nlon = 144
  else
    nlat = 10
    nlon = 18
    ; -- midpoints of 18x10 grid cells with approx. equal area:
    lat_new   = (/-71.5500,-45.3795,-30.5306,-17.7528,-5.8381,5.8381,17.7528,30.5306,45.3795,71.5500/)
    lon_new   = fspan(10,350,18)
    ; -- bounds to those midpoints:
    lat_bnds  = (/-90.0,-53.1,-36.8,-23.5,-11.5,0.0,11.5,23.5,36.8,53.1,90.0/)
    lon_bnds  = fspan(0,360,19)
    nbnd      = 2
    lat = lat_new
    lon = lon_new
    lat_bnds1  = new((/nlat,nbnd/),double)
    lon_bnds1  = new((/nlon,nbnd/),double)
    lat_bnds1(:,0) = lat_bnds(0:nlat-1)
    lat_bnds1(:,1) = lat_bnds(1:nlat)
    lon_bnds1(:,0) = lon_bnds(0:nlon-1)
    lon_bnds1(:,1) = lon_bnds(1:nlon)
  end if

  ; -- which calculation of uncertainties to plot
  if constrain .eq. "no"
    tags    = (/"cmip5","cmip6","le","le_hs09"/)
  else
    tags    = (/"cmip6"/)
  end if
;
;-----------------------------------------------------------
; -- regrid resources
if regrid1 .eq. "yes" .or. regrid2 .eq. "yes" .or. regrid3 .eq. "yes"
   Opt                	= True
  Opt@DstGridLat    	= lat_new
  Opt@DstGridLon    	= lon_new
  ; Opt@Overwrite      	= True
  Opt@ForceOverwrite 	= True
  Opt@InterpMethod  	= "conserve" ; conserve bilinear patch
  Opt@Pole	          = "none"
  Opt@SkipSrcGrid	    = False
  Opt@SkipDstGrid	    = False
  Opt@SkipWgtGen	    = False

  ; del_FillValue(lat_bnds1)
  ; del_FillValue(lon_bnds1)
  lat@bounds  = "lat_bnds"
  lon@bounds  = "lon_bnds"

  var_x_bnds      = (/ "lon_bnds" /)
  var_x_bnds_type = (/ "double" /)
  var_y_bnds      = (/ "lat_bnds" /)
  var_y_bnds_type = (/ "double" /)
  var_t_bnds      = (/ "time_bnds" /)
  var_t_bnds_type = (/ "double" /)
  var_t           = (/ "time" /)
  var_t_type      = (/ "double" /)
  dimNames = (/ "lon", "lat", "bnds", "time" /)
  dimSizes = (/  nlon,  nlat,   nbnd, -1     /)
  dimUnlim = (/ False, False,  False, True    /)
end if

;-----------------------------------------------------------

do v = 0,1 ; full: 0,1     ; tas pr
  var = vars(v)
  print_clock("variable = "+var)
  do s = 2,2 ; full: 0,2      ; DJF JJA annual
    seas = seasons(s)

    ; -- LEs --
    if recalc1 .eq. "yes"
      le_ensmean_all      = new((/ende_trend-start_trend+1,dimsizes(models),nlat,nlon/),float) ; ensemble mean from each model
      le_intvar_all       = new((/ende_trend-start_trend+1,dimsizes(models),max(ensmem0),nlat,nlon/),float) ; variance across ensemble from each model
      le_ensmean_hs09_all = le_intvar_all
      le_intvar_hs09_all  = le_intvar_all
      do m = 0,dimsizes(models)-1
        print_clock("model = "+models(m)+" ("+m+"/"+(dimsizes(models)-1)+")")
        start   = start0(m)
        ende    = ende0(m)
        ensmem  = ensmem0(m)
        ; fil    = systemfunc("ls "+pathin+models(m)+"_lens/"+comp+"/"+var+"/"+var+"_"+comp+"_*_historical_rcp85_r2i1p1_*g025.nc")
        ; print(fil)
        ; a = addfile(fil(0),"r")
        ; dummy   = a->$var$
        ; dims    = dimsizes(dummy)
        ; ntime   = dims(0)
        ; nlat    = dims(1)
        ; nlon    = dims(2)
        data    = new((/ende_trend-start_trend+1,ensmem,nlat,nlon/),float)
        do e = 1,ensmem
          print_clock("  ensmem = "+e+"/"+ensmem)
          if regrid1 .eq. "yes"
            fil    = systemfunc("ls "+pathin+models(m)+"_lens/"+comp+"/"+var+"/"+var+"_"+comp+"_*_historical_rcp85_r"+e+"i1p1_*_g025.nc")
            a = addfile(fil(0),"r")
            raw       = a->$var$
            tmp0      = ESMF_regrid(raw,Opt)
            tmp0&lat  = lat
            tmp0&lon  = lon
            print("  save regrided data...")
            ofile1 = pathin+models(m)+"_lens/"+comp+"/"+var+"/"+var+"_"+comp+"_"+model_names(m)+"_historical_rcp85_r"+e+"i1p1_"+start+"01-"+ende+"12_"+grid+".nc"
            system("rm -rf "+ofile1)
            o = addfile(ofile1,"c")
            filedimdef(o, dimNames, dimSizes, dimUnlim )
            filevardef(o, var_x_bnds, var_x_bnds_type, (/"lon","bnds"/) )
            filevardef(o, var_y_bnds, var_y_bnds_type, (/"lat","bnds"/) )
            filevardef(o, var_t_bnds, var_t_bnds_type, (/"time","bnds"/) )
            filevardef(o, var_t, var_t_type, (/"time"/))
            o->$var$    = tmp0
            o->lat      = lat
            o->lat_bnds = lat_bnds1
            o->lon      = lon
            o->lon_bnds = lon_bnds1
            delete(raw)
            raw = tmp0
            delete(tmp0)
          else
            fil    = systemfunc("ls "+pathin+models(m)+"_lens/"+comp+"/"+var+"/"+var+"_"+comp+"_*_historical_rcp85_r"+e+"i1p1_*_"+grid+".nc")
            a = addfile(fil(0),"r")
            raw  = a->$var$
          end if
          if seas .eq. "annual"
            raw_seas = month_to_annual(raw,1)
          else
            raw_seas = month_to_season(raw,seas)
          end if
          raw_seas_cut = raw_seas(start_trend-start:ende_trend-start,:,:) ; cut to common length
          raw_seas_cut!0 = "time"
          ; if e .eq. 1
            ; -- Hawkins&Sutton calculation --
            ; -- temporarily save data:
            ofile1  = pathin+"tmp999.nc"
            system("rm -rf "+ofile1)
            o = addfile(ofile1,"c")
            o->tmp  = raw_seas_cut
            ; -- call matlab to fit data:
            system("cd /home/flehner/scripts/matlab ; echo 'polyfit_3d_netcdf2' |matlab -display /dev/null")
            ; -- load fitted data:
            o = addfile(ofile1,"r")
            fit     = o->tmp
            copy_VarCoords(raw_seas_cut,fit)
            ; raw_seas_cut&time = trend_time
            ; x                 = conform_dims(dimsizes(raw_seas_cut({lat|:},{lon|:},{time|:})),trend_time,(/2/))
            ; c                 = lspoly_n(x,raw_seas_cut({lat|:},{lon|:},{time|:}),1.,5,2)
            ; tmp99             = raw_seas_cut({lat|:},{lon|:},{time|:})
            ; do i = 0,dimsizes(trend_time)-1
            ;   tmp99(:,:,i) = c(:,:,0) + c(:,:,1)*x(:,:,i) + c(:,:,2)*x(:,:,i)^2 + c(:,:,3)*x(:,:,i)^3 + c(:,:,4)*x(:,:,i)^4
            ; end do
            ; fit = tmp99({time|:},{lat|:},{lon|:})
            ; delete(tmp99)
            ; -- additional smoothing of polynomial to get forced response:
            fit_smoothed  = runave_n_Wrap(fit,10,0,0)
            ; -- calculate anomalies to reference period:
            ref   = conform_dims(dimsizes(raw_seas_cut),dim_avg_n_Wrap(fit_smoothed(refstart-start_trend:refende-start_trend,:,:),0),(/1,2/))
            if var .eq. "pr"
              le_ensmean_hs09_all(:,m,e-1,:,:) = tofloat( ((fit_smoothed - ref) / ref) * 100)
            else
              le_ensmean_hs09_all(:,m,e-1,:,:) = tofloat(fit_smoothed - ref )
            end if
            ; -- subtract forced response to get residual (still at seas resolution, so noisy):
            residual = (raw_seas_cut-fit_smoothed) + conform_dims(dimsizes(raw_seas_cut),dim_avg_n_Wrap(fit_smoothed,0),(/1,2/)) ; do residual, add back long-term mean from fit
            if var .eq. "pr"
              ref_tmp = conform_dims(dimsizes(residual),dim_avg_n_Wrap(residual(refstart-start_trend:refende-start_trend,:,:),0),(/1,2/))
              le_intvar_hs09_all(:,m,e-1,:,:) = tofloat( (residual / ref_tmp) * 100 )
              delete([/ref_tmp/])
            else
              le_intvar_hs09_all(:,m,e-1,:,:) = tofloat( residual )
            end if
            copy_VarCoords(raw_seas_cut,le_ensmean_hs09_all(:,0,0,:,:))
            copy_VarCoords(le_ensmean_hs09_all,le_intvar_hs09_all)
            ; delete([/raw,raw_seas,raw_seas_cut,fit,fit_smoothed/])
            delete([/fit,fit_smoothed/])
          ; else
            ; -- SMILE calculation --
            data(:,e-1,:,:) = raw_seas_cut
            copy_VarCoords(raw_seas_cut,data(:,0,:,:))
            ref_tmp = conform_dims(dimsizes(data(:,0,:,:)),dim_avg_n_Wrap(data(refstart-start_trend:refende-start_trend,e-1,:,:),0),(/1,2/))
            if var .eq. "pr"
              le_intvar_all(:,m,e-1,:,:)   = (data(:,e-1,:,:) / ref_tmp) * 100
            else
              le_intvar_all(:,m,e-1,:,:)   = data(:,e-1,:,:)
            end if
            delete([/raw,raw_seas,ref_tmp/])
          ; end if
        end do
        forced             = runave_n_Wrap(dim_avg_n_Wrap(data,1),10,0,0)
        ; -- calculate anomalies to reference period:
        ref             = conform_dims(dimsizes(forced),dim_avg_n_Wrap(forced(refstart-start_trend:refende-start_trend,:,:),0),(/1,2/))
        if var .eq. "pr"
          le_ensmean_all(:,m,:,:) = tofloat( ((forced - ref) / ref) * 100)
        else
          le_ensmean_all(:,m,:,:) = tofloat( forced - ref )
        end if
        copy_VarCoords(raw_seas_cut,le_ensmean_all(:,0,:,:))
        copy_VarCoords(raw_seas_cut,le_intvar_all(:,0,0,:,:))
        delete([/data,fil,forced,raw_seas_cut/])
      end do
      ; -- save data --
      print_clock("save SMILEs data...")
      setfileoption("nc","format","netcdf4")
      ofile1 = pathin+"/"+var+"_"+comp+"_all_lens_forced_and_intvar_"+grid+"_"+seas+"_ref"+refstart+"-"+refende+".nc"
      system("rm -rf "+ofile1)
      o = addfile(ofile1,"c")
      o->ensmean_all      = le_ensmean_all
      o->intvar_all       = le_intvar_all
      o->ensmean_hs09_all = le_ensmean_hs09_all
      o->intvar_hs09_all  = le_intvar_hs09_all
    else
      ; -- load data --
      print_clock("load SMILEs data...")
      ofile1 = pathin+"/"+var+"_"+comp+"_all_lens_forced_and_intvar_"+grid+"_"+seas+"_ref"+refstart+"-"+refende+".nc"
      o = addfile(ofile1,"r")
      le_ensmean_all        = o->ensmean_all
      le_intvar_all         = o->intvar_all
      le_ensmean_hs09_all   = o->ensmean_hs09_all
      le_intvar_hs09_all    = o->intvar_hs09_all
    end if


    ; -- CMIP5 --
    scen  = (/"rcp26","rcp45","rcp85"/)
    if recalc2 .eq. "yes"
      start = 1870
      ende  = 2100
      do sc = 0,2 ; 0,2 ; scenarios
        cmip5_ensmean_all = new((/ende_trend-start_trend+1,dimsizes(cmip5_models),nlat,nlon/),float) ; ensemble mean estimated via 4th order poly from each CMIP5 model
        cmip5_intvar_all  = new((/ende_trend-start_trend+1,dimsizes(cmip5_models),nlat,nlon/),float) ; int var estimated from residual after subtracting the 4th order poly from each CMIP5 model
        do m = 0,dimsizes(cmip5_models)-1
          print_clock("scenario = "+scen(sc)+" | model = "+cmip5_models(m)+" ("+m+"/"+(dimsizes(cmip5_models)-1)+")")
          if regrid2 .eq. "yes" ; .and. m .eq. 1
            fil    = systemfunc("ls "+pathin+"cmip5-ng/"+var+"/"+var+"_mon_"+cmip5_models(m)+"_"+scen(sc)+"_r1i1p1_native.nc")
            a = addfile(fil(0),"r")
            tmp0    = a->$var$
            raw   = ESMF_regrid(tmp0,Opt)
            print("  save regrided data...")
            ofile1 = pathin+"cmip5-ng/"+var+"/"+var+"_mon_"+cmip5_models(m)+"_"+scen(sc)+"_r1i1p1_"+grid+".nc"
            system("rm -rf "+ofile1)
            o = addfile(ofile1,"c")
            filedimdef(o, dimNames, dimSizes, dimUnlim )
            filevardef(o, var_x_bnds, var_x_bnds_type, (/"lon","bnds"/) )
            filevardef(o, var_y_bnds, var_y_bnds_type, (/"lat","bnds"/) )
            filevardef(o, var_t_bnds, var_t_bnds_type, (/"time","bnds"/) )
            filevardef(o, var_t, var_t_type, (/"time"/))
            o->$var$    = raw
            o->lat      = lat
            o->lat_bnds = lat_bnds1
            o->lon      = lon
            o->lon_bnds = lon_bnds1
            delete(tmp0)
          else
            fil   = systemfunc("ls "+pathin+"cmip5-ng/"+var+"/"+var+"_mon_"+cmip5_models(m)+"_"+scen(sc)+"_r1i1p1_"+grid+".nc")
            a     = addfile(fil(0),"r")
            raw   = a->$var$
          end if
          if seas .eq. "annual"
            raw_seas = month_to_annual(raw,1)
          else
            raw_seas = month_to_season(raw,seas)
          end if
          raw_seas_cut = raw_seas(start_trend-start:ende_trend-start,:,:) ; cut to common length
          raw_seas_cut!0 = "time"
          ; -- Hawkins&Sutton calculation --
          ; -- temporarily save data:
          ofile1  = pathin+"tmp999.nc"
          system("rm -rf "+ofile1)
          o = addfile(ofile1,"c")
          o->tmp  = raw_seas_cut
          ; -- call matlab to fit data:
          system("cd /home/flehner/scripts/matlab ; echo 'polyfit_3d_netcdf2' |matlab -display /dev/null")
          ; -- load fitted data:
          o = addfile(ofile1,"r")
          fit     = o->tmp
          copy_VarCoords(raw_seas_cut,fit)
          ; -- additional smoothing of polynomial to get forced response:
          fit_smoothed  = runave_n_Wrap(fit,10,0,0)
          ; -- calculate anomalies to reference period:
          ref_fit   = conform_dims(dimsizes(raw_seas_cut),dim_avg_n_Wrap(fit_smoothed(refstart-start_trend:refende-start_trend,:,:),0),(/1,2/))
          if var .eq. "pr"
            ref_fit   = where(ref_fit .eq. 0,.001,ref_fit)
            cmip5_ensmean_all(:,m,:,:) = tofloat( ((fit_smoothed - ref_fit) / ref_fit) * 100 )
          else
            cmip5_ensmean_all(:,m,:,:) = tofloat( fit_smoothed - ref_fit )
          end if
          ; -- subtract forced response to get residual (still at seas resolution, so noisy):
          residual = (raw_seas_cut-fit_smoothed) + conform_dims(dimsizes(raw_seas_cut),dim_avg_n_Wrap(fit_smoothed,0),(/1,2/)) ; do residual, add back long-term mean from fit
          if var .eq. "pr"
            ref_tmp = conform_dims(dimsizes(residual),dim_avg_n_Wrap(residual(refstart-start_trend:refende-start_trend,:,:),0),(/1,2/))
            ref_tmp = where(ref_tmp .eq. 0,.001,ref_tmp)
            cmip5_intvar_all(:,m,:,:) = tofloat( (residual / ref_tmp) * 100 )
            delete([/ref_tmp/])
          else
            cmip5_intvar_all(:,m,:,:) = tofloat( residual )
          end if
          copy_VarCoords(raw_seas_cut,cmip5_ensmean_all(:,0,:,:))
          copy_VarCoords(cmip5_ensmean_all,cmip5_intvar_all)
          delete([/raw,raw_seas,raw_seas_cut,fit,fit_smoothed/])
        end do
        ; -- save data --
        print_clock("save CMIP5 data...")
        setfileoption("nc","format","netcdf4")
        ofile1 = pathin+"cmip5-ng/"+var+"/"+var+"_mon_all_models_forced_and_intvar_"+scen(sc)+"_r1i1p1_"+grid+"_"+seas+"_ref"+refstart+"-"+refende+"_1.nc"
        system("rm -rf "+ofile1)
        o = addfile(ofile1,"c")
        o->ensmean_all  = cmip5_ensmean_all
        ofile1 = pathin+"cmip5-ng/"+var+"/"+var+"_mon_all_models_forced_and_intvar_"+scen(sc)+"_r1i1p1_"+grid+"_"+seas+"_ref"+refstart+"-"+refende+"_2.nc"
        system("rm -rf "+ofile1)
        o = addfile(ofile1,"c")
        o->intvar_all   = cmip5_intvar_all
        if sc .eq. 0
          cmip5_ensmean_rcp26_all  = cmip5_ensmean_all
          cmip5_intvar_rcp26_all   = cmip5_intvar_all
        else if sc .eq. 1
          cmip5_ensmean_rcp45_all  = cmip5_ensmean_all
          cmip5_intvar_rcp45_all   = cmip5_intvar_all
        else if sc .eq. 2
          cmip5_ensmean_rcp85_all  = cmip5_ensmean_all
          cmip5_intvar_rcp85_all   = cmip5_intvar_all
        end if
        end if
        end if
      end do
    else
      ; -- load data --
      print_clock("load CMIP5 data...")
      ; -- rcp26
      sc = 0
      ofile1 = pathin+"cmip5-ng/"+var+"/"+var+"_mon_all_models_forced_and_intvar_"+scen(sc)+"_r1i1p1_"+grid+"_"+seas+"_ref"+refstart+"-"+refende+"_1.nc"
      o = addfile(ofile1,"r")
      cmip5_ensmean_rcp26_all  = o->ensmean_all
      ofile1 = pathin+"cmip5-ng/"+var+"/"+var+"_mon_all_models_forced_and_intvar_"+scen(sc)+"_r1i1p1_"+grid+"_"+seas+"_ref"+refstart+"-"+refende+"_2.nc"
      o = addfile(ofile1,"r")
      cmip5_intvar_rcp26_all   = o->intvar_all
      ; -- rcp45
      sc = 1
      ofile1 = pathin+"cmip5-ng/"+var+"/"+var+"_mon_all_models_forced_and_intvar_"+scen(sc)+"_r1i1p1_"+grid+"_"+seas+"_ref"+refstart+"-"+refende+"_1.nc"
      o = addfile(ofile1,"r")
      cmip5_ensmean_rcp45_all  = o->ensmean_all
      ofile1 = pathin+"cmip5-ng/"+var+"/"+var+"_mon_all_models_forced_and_intvar_"+scen(sc)+"_r1i1p1_"+grid+"_"+seas+"_ref"+refstart+"-"+refende+"_2.nc"
      o = addfile(ofile1,"r")
      cmip5_intvar_rcp45_all   = o->intvar_all
      ; -- rcp85
      sc = 2
      ofile1 = pathin+"cmip5-ng/"+var+"/"+var+"_mon_all_models_forced_and_intvar_"+scen(sc)+"_r1i1p1_"+grid+"_"+seas+"_ref"+refstart+"-"+refende+"_1.nc"
      o = addfile(ofile1,"r")
      cmip5_ensmean_rcp85_all  = o->ensmean_all
      ofile1 = pathin+"cmip5-ng/"+var+"/"+var+"_mon_all_models_forced_and_intvar_"+scen(sc)+"_r1i1p1_"+grid+"_"+seas+"_ref"+refstart+"-"+refende+"_2.nc"
      o = addfile(ofile1,"r")
      cmip5_intvar_rcp85_all   = o->intvar_all
    end if



    ; -- CMIP6 --
    delete(scen)
    scen  = (/"ssp126","ssp245","ssp370","ssp585"/)
    if recalc3 .eq. "yes"
      start = 1850
      ende  = 2100
      do sc = 0,3 ; 0,3 ; scenarios
        cmip6_ensmean_all = new((/ende_trend-start_trend+1,dimsizes(cmip6_models),nlat,nlon/),float) ; ensemble mean estimated via 4th order poly from each cmip6 model
        cmip6_intvar_all  = new((/ende_trend-start_trend+1,dimsizes(cmip6_models),nlat,nlon/),float) ; int var estimated from residual after subtracting the 4th order poly from each cmip6 model
        do m = 0,dimsizes(cmip6_models)-1
          print_clock("scenario = "+scen(sc)+" | model = "+cmip6_models(m)+" ("+m+"/"+(dimsizes(cmip6_models)-1)+")")
          if regrid3 .eq. "yes" ; .and. m .eq. 1
            fil   = systemfunc("ls "+pathin+"cmip6-ng/"+var+"/"+var+"_mon_"+cmip6_models(m)+"_"+scen(sc)+"_*_g025.nc")
            ; fil   = systemfunc("ls /net/cfc/cmip6/Next_Generation/"+var+"/mon/g025/"+var+"_mon_"+cmip6_models(m)+"_"+scen(sc)+"_*_g025.nc")
            ; fil   = systemfunc("ls /net/meso/climphys/flehner/scratch/"+var+"/"+var+"_mon_"+cmip6_models(m)+"_"+scen(sc)+"_*_g025.nc")
            fil0  = fil(0)
            ripf  = systemfunc("echo "+fil0+" | awk -F'[_.]' '{print $5}'")
            a = addfile(fil0,"r")
            tmp0    = a->$var$
            raw   = ESMF_regrid(tmp0,Opt)
            print("  save regrided data...")
            ofile1 = pathin+"cmip6-ng/"+var+"/"+var+"_mon_"+cmip6_models(m)+"_"+scen(sc)+"_"+ripf+"_"+grid+".nc"
            system("rm -rf "+ofile1)
            o = addfile(ofile1,"c")
            filedimdef(o, dimNames, dimSizes, dimUnlim )
            filevardef(o, var_x_bnds, var_x_bnds_type, (/"lon","bnds"/) )
            filevardef(o, var_y_bnds, var_y_bnds_type, (/"lat","bnds"/) )
            filevardef(o, var_t_bnds, var_t_bnds_type, (/"time","bnds"/) )
            filevardef(o, var_t, var_t_type, (/"time"/))
            o->$var$    = raw
            o->lat      = lat
            o->lat_bnds = lat_bnds1
            o->lon      = lon
            o->lon_bnds = lon_bnds1
            delete([/tmp0,fil,fil0,ripf/])
          else
            fil    = systemfunc("ls "+pathin+"cmip6-ng/"+var+"/"+var+"_mon_"+cmip6_models(m)+"_"+scen(sc)+"_*_"+grid+".nc")
            a = addfile(fil(0),"r")
            raw  = a->$var$
            delete(fil)
          end if
          if seas .eq. "annual"
            raw_seas = month_to_annual(raw,1)
          else
            raw_seas = month_to_season(raw,seas)
          end if
          raw_seas_cut = raw_seas(start_trend-start:ende_trend-start,:,:) ; cut to common length
          raw_seas_cut!0 = "time"
          ; -- Hawkins&Sutton calculation --
          ; -- temporarily save data:
          ofile1  = pathin+"tmp.nc"
          system("rm -rf "+ofile1)
          o = addfile(ofile1,"c")
          o->tmp  = raw_seas_cut
          ; -- call matlab to fit data:
          system("cd /home/flehner/scripts/matlab ; echo 'polyfit_3d_netcdf2' |matlab -display /dev/null")
          ; -- load fitted data:
          o = addfile(ofile1,"r")
          fit     = o->tmp
          copy_VarCoords(raw_seas_cut,fit)
          ; -- additional smoothing of polynomial to get forced response:
          fit_smoothed  = runave_n_Wrap(fit,10,0,0)
          ; -- calculate anomalies to reference period:
          ref_fit   = conform_dims(dimsizes(raw_seas_cut),dim_avg_n_Wrap(fit_smoothed(refstart-start_trend:refende-start_trend,:,:),0),(/1,2/))
          if var .eq. "pr"
            ref_fit   = where(ref_fit .eq. 0,.001,ref_fit)
            cmip6_ensmean_all(:,m,:,:) = tofloat( ((fit_smoothed - ref_fit) / ref_fit) * 100 )
          else
            cmip6_ensmean_all(:,m,:,:) = tofloat( fit_smoothed - ref_fit )
          end if
          ; -- subtract forced response to get residual (still at seas resolution, so noisy):
          residual = (raw_seas_cut-fit_smoothed) + conform_dims(dimsizes(raw_seas_cut),dim_avg_n_Wrap(fit_smoothed,0),(/1,2/)) ; do residual, add back long-term mean from fit
          if var .eq. "pr"
            ref_tmp = conform_dims(dimsizes(residual),dim_avg_n_Wrap(residual(refstart-start_trend:refende-start_trend,:,:),0),(/1,2/))
            ref_tmp = where(ref_tmp .eq. 0,.001,ref_tmp)
            cmip6_intvar_all(:,m,:,:) = tofloat( (residual / ref_tmp) * 100 )
            delete([/ref_tmp/])
          else
            cmip6_intvar_all(:,m,:,:) = tofloat( residual )
          end if
          copy_VarCoords(raw_seas_cut,cmip6_ensmean_all(:,0,:,:))
          copy_VarCoords(cmip6_ensmean_all,cmip6_intvar_all)
          delete([/raw,raw_seas,raw_seas_cut,fit,fit_smoothed/])

        end do
        ; -- save data --
        print_clock("save CMIP6 data...")
        setfileoption("nc","format","netcdf4")
        ofile1 = pathin+"cmip6-ng/"+var+"/"+var+"_mon_all_models_forced_and_intvar_"+scen(sc)+"_r1i1p1f1_"+grid+"_"+seas+"_ref"+refstart+"-"+refende+"_1"+ending+".nc"
        system("rm -rf "+ofile1)
        o = addfile(ofile1,"c")
        o->ensmean_all  = cmip6_ensmean_all
        ofile1 = pathin+"cmip6-ng/"+var+"/"+var+"_mon_all_models_forced_and_intvar_"+scen(sc)+"_r1i1p1f1_"+grid+"_"+seas+"_ref"+refstart+"-"+refende+"_2"+ending+".nc"
        system("rm -rf "+ofile1)
        o = addfile(ofile1,"c")
        o->intvar_all   = cmip6_intvar_all
        if sc .eq. 0
          cmip6_ensmean_ssp126_all  = cmip6_ensmean_all
          cmip6_intvar_ssp126_all   = cmip6_intvar_all
        else if sc .eq. 1
          cmip6_ensmean_ssp245_all  = cmip6_ensmean_all
          cmip6_intvar_ssp245_all   = cmip6_intvar_all
        else if sc .eq. 2
          cmip6_ensmean_ssp370_all  = cmip6_ensmean_all
          cmip6_intvar_ssp370_all   = cmip6_intvar_all
        else if sc .eq. 3
          cmip6_ensmean_ssp585_all  = cmip6_ensmean_all
          cmip6_intvar_ssp585_all   = cmip6_intvar_all
        end if
        end if
        end if
        end if
      end do
    else
      ; -- load data --
      print_clock("load CMIP6 data...")
      ; -- ssp126
      sc = 0
      ofile1 = pathin+"cmip6-ng/"+var+"/"+var+"_mon_all_models_forced_and_intvar_"+scen(sc)+"_r1i1p1f1_"+grid+"_"+seas+"_ref"+refstart+"-"+refende+"_1"+ending+".nc"
      o = addfile(ofile1,"r")
      cmip6_ensmean_ssp126_all  = o->ensmean_all
      ofile1 = pathin+"cmip6-ng/"+var+"/"+var+"_mon_all_models_forced_and_intvar_"+scen(sc)+"_r1i1p1f1_"+grid+"_"+seas+"_ref"+refstart+"-"+refende+"_2"+ending+".nc"
      o = addfile(ofile1,"r")
      cmip6_intvar_ssp126_all   = o->intvar_all
      ; -- ssp245
      sc = 1
      ofile1 = pathin+"cmip6-ng/"+var+"/"+var+"_mon_all_models_forced_and_intvar_"+scen(sc)+"_r1i1p1f1_"+grid+"_"+seas+"_ref"+refstart+"-"+refende+"_1"+ending+".nc"
      o = addfile(ofile1,"r")
      cmip6_ensmean_ssp245_all  = o->ensmean_all
      ofile1 = pathin+"cmip6-ng/"+var+"/"+var+"_mon_all_models_forced_and_intvar_"+scen(sc)+"_r1i1p1f1_"+grid+"_"+seas+"_ref"+refstart+"-"+refende+"_2"+ending+".nc"
      o = addfile(ofile1,"r")
      cmip6_intvar_ssp245_all   = o->intvar_all
      ; -- ssp370
      sc = 2
      ofile1 = pathin+"cmip6-ng/"+var+"/"+var+"_mon_all_models_forced_and_intvar_"+scen(sc)+"_r1i1p1f1_"+grid+"_"+seas+"_ref"+refstart+"-"+refende+"_1"+ending+".nc"
      o = addfile(ofile1,"r")
      cmip6_ensmean_ssp370_all  = o->ensmean_all
      ofile1 = pathin+"cmip6-ng/"+var+"/"+var+"_mon_all_models_forced_and_intvar_"+scen(sc)+"_r1i1p1f1_"+grid+"_"+seas+"_ref"+refstart+"-"+refende+"_2"+ending+".nc"
      o = addfile(ofile1,"r")
      cmip6_intvar_ssp370_all   = o->intvar_all
      ; -- ssp585
      sc = 3
      ofile1 = pathin+"cmip6-ng/"+var+"/"+var+"_mon_all_models_forced_and_intvar_"+scen(sc)+"_r1i1p1f1_"+grid+"_"+seas+"_ref"+refstart+"-"+refende+"_1"+ending+".nc"
      o = addfile(ofile1,"r")
      cmip6_ensmean_ssp585_all  = o->ensmean_all
      ofile1 = pathin+"cmip6-ng/"+var+"/"+var+"_mon_all_models_forced_and_intvar_"+scen(sc)+"_r1i1p1f1_"+grid+"_"+seas+"_ref"+refstart+"-"+refende+"_2"+ending+".nc"
      o = addfile(ofile1,"r")
      cmip6_intvar_ssp585_all   = o->intvar_all
    end if
    delete(scen)




; -- UNCERTAINTIES --
print_clock("calculate uncertainties...")
rmwl = 10 ; running mean window length in years

; -- scenario uncertainty --
; -- like Froelicher et al. 2015:
; scen_u_cmip5              = abs(dim_avg_n_Wrap(cmip5_ensmean_rcp85_all,1)-dim_avg_n_Wrap(cmip5_ensmean_rcp26_all,1))
; -- like HS09:
tmp               = new((/ntime,3,nlat,nlon/),float)
tmp(:,0,:,:)      = dim_avg_n_Wrap(cmip5_ensmean_rcp26_all,1)
tmp(:,1,:,:)      = dim_avg_n_Wrap(cmip5_ensmean_rcp45_all,1)
tmp(:,2,:,:)      = dim_avg_n_Wrap(cmip5_ensmean_rcp85_all,1)
scen_u_cmip5            = dim_variance_n_Wrap(tmp,1)
delete(tmp)

tmp               = new((/ntime,4,nlat,nlon/),float)
tmp(:,0,:,:)      = dim_avg_n_Wrap(cmip6_ensmean_ssp126_all,1)
tmp(:,1,:,:)      = dim_avg_n_Wrap(cmip6_ensmean_ssp245_all,1)
tmp(:,2,:,:)      = dim_avg_n_Wrap(cmip6_ensmean_ssp370_all,1)
tmp(:,3,:,:)      = dim_avg_n_Wrap(cmip6_ensmean_ssp585_all,1)
scen_u_cmip6      = dim_variance_n_Wrap(tmp,1)
delete(tmp)

model_u_le          = dim_variance_n_Wrap(le_ensmean_all,1)
model_u_le_hs09     = dim_variance_n_Wrap(le_ensmean_hs09_all(:,:,0,:,:),1)
; model_u_le_hs09_all = dim_variance_n_Wrap(le_ensmean_hs09_all,2) ; NOTE: this treats individual ensemble members as forced response to calculate a mock "model uncertainty"
model_u_le_hs09_all = dim_avg_n_Wrap(dim_variance_n_Wrap(le_ensmean_hs09_all,2),1) ; NOTE: this treats individual ensemble members as forced response to calculate a mock "model uncertainty"
if mu .eq. 1
  model_u_cmip5       = dim_variance_n_Wrap(cmip5_ensmean_rcp85_all,1)
  model_u_cmip6       = dim_variance_n_Wrap(cmip6_ensmean_ssp585_all,1)
else
  delete(tmp)
  dims                = dimsizes(dim_variance_n_Wrap(cmip5_ensmean_rcp26_all,1))
  tmp                 = new((/3,dims(0),dims(1),dims(2)/),float)
  tmp(0,:,:,:)        = dim_variance_n_Wrap(cmip5_ensmean_rcp26_all,1)
  tmp(1,:,:,:)        = dim_variance_n_Wrap(cmip5_ensmean_rcp45_all,1)
  tmp(2,:,:,:)        = dim_variance_n_Wrap(cmip5_ensmean_rcp85_all,1)
  model_u_cmip5       = dim_avg_n_Wrap(tmp,0)
  delete(tmp)
  dims                = dimsizes(dim_variance_n_Wrap(cmip6_ensmean_ssp126_all,1))
  tmp                 = new((/4,dims(0),dims(1),dims(2)/),float)
  tmp(0,:,:,:)        = dim_variance_n_Wrap(cmip6_ensmean_ssp126_all,1)
  tmp(1,:,:,:)        = dim_variance_n_Wrap(cmip6_ensmean_ssp245_all,1)
  tmp(2,:,:,:)        = dim_variance_n_Wrap(cmip6_ensmean_ssp370_all,1)
  tmp(3,:,:,:)        = dim_variance_n_Wrap(cmip6_ensmean_ssp585_all,1)
  model_u_cmip6       = dim_avg_n_Wrap(tmp,0)
  delete(tmp)
end if

int_u_mean_le       = dim_avg_n_Wrap(dim_variance_n_Wrap(runave_n_Wrap(le_intvar_all,rmwl,0,0),2),1) ; first smooth, then variance across members, then mean across SMILEs
int_u_mean_le_hs09  = conform_dims(dimsizes(scen_u_cmip5),dim_avg_n_Wrap(dim_variance_n_Wrap(runave_n_Wrap(le_intvar_hs09_all(:,:,0,:,:),rmwl,0,0),0),0),(/1,2/)); first smooth, then variance across time, then mean across SMILEs, and broadcasting
int_u_mean_cmip5    = conform_dims(dimsizes(scen_u_cmip5),dim_avg_n_Wrap(dim_variance_n_Wrap(runave_n_Wrap(cmip5_intvar_rcp85_all,rmwl,0,0),0),0),(/1,2/)); first smooth, then variance across time, then mean across CMIP5 models, and broadcasting
int_u_mean_cmip6    = conform_dims(dimsizes(scen_u_cmip5),dim_avg_n_Wrap(dim_variance_n_Wrap(runave_n_Wrap(cmip6_intvar_ssp585_all,rmwl,0,0),0),0),(/1,2/)); first smooth, then variance across time, then mean across CMIP6 models, and broadcasting

total_u_le          = model_u_le + int_u_mean_le + scen_u_cmip5
total_u_le_hs09     = model_u_le_hs09 + int_u_mean_le_hs09 + scen_u_cmip5
; total_u_le_hs09_all = model_u_le_hs09_all
; do m = 0,dimsizes(models)-1
;   total_u_le_hs09_all(:,m,:,:) = model_u_le_hs09_all(:,m,:,:) + int_u_mean_le_hs09 + scen_u_cmip5
; end do
total_u_le_hs09_all = model_u_le_hs09_all + int_u_mean_le + scen_u_cmip5
total_u_cmip5       = model_u_cmip5 + int_u_mean_cmip5 + scen_u_cmip5
total_u_cmip6       = model_u_cmip6 + int_u_mean_cmip6 + scen_u_cmip6

copy_VarCoords(cmip5_ensmean_rcp85_all(:,0,:,:),scen_u_cmip5)
copy_VarCoords(cmip5_ensmean_rcp85_all(:,0,:,:),scen_u_cmip6)
copy_VarCoords(cmip5_ensmean_rcp85_all(:,0,:,:),int_u_mean_cmip5)
copy_VarCoords(cmip6_ensmean_ssp585_all(:,0,:,:),int_u_mean_cmip6)


; -- fractional uncertainty --
int_u_mean_le_frac      = (int_u_mean_le/total_u_le)*100
model_u_le_frac         = (model_u_le/total_u_le)*100
scen_u_le_frac          = (scen_u_cmip5/total_u_le)*100
copy_VarCoords(model_u_le,int_u_mean_le_frac)
copy_VarCoords(model_u_le,model_u_le_frac)
copy_VarCoords(model_u_le,scen_u_le_frac)

int_u_mean_le_hs09_frac = (int_u_mean_le_hs09/total_u_le_hs09)*100
model_u_le_hs09_frac    = (model_u_le_hs09/total_u_le_hs09)*100
scen_u_le_hs09_frac     = (scen_u_cmip5/total_u_le_hs09)*100
copy_VarCoords(model_u_le,int_u_mean_le_hs09_frac)
copy_VarCoords(model_u_le,model_u_le_hs09_frac)
copy_VarCoords(model_u_le,scen_u_le_hs09_frac)

int_u_mean_le_hs09_all_frac = (int_u_mean_le/total_u_le_hs09_all)*100
model_u_le_hs09_all_frac    = (model_u_le_hs09_all/total_u_le_hs09_all)*100
scen_u_le_hs09_all_frac     = (scen_u_cmip5/total_u_le_hs09_all)*100
copy_VarCoords(model_u_le,int_u_mean_le_hs09_all_frac)
copy_VarCoords(model_u_le,model_u_le_hs09_all_frac)
copy_VarCoords(model_u_le,scen_u_le_hs09_all_frac)

int_u_mean_cmip5_frac   = (int_u_mean_cmip5/total_u_cmip5)*100
model_u_cmip5_frac      = (model_u_cmip5/total_u_cmip5)*100
scen_u_cmip5_frac       = (scen_u_cmip5/total_u_cmip5)*100
copy_VarCoords(model_u_cmip5,int_u_mean_cmip5_frac)
copy_VarCoords(model_u_cmip5,model_u_cmip5_frac)
copy_VarCoords(model_u_cmip5,scen_u_cmip5_frac)

int_u_mean_cmip6_frac   = (int_u_mean_cmip6/total_u_cmip6)*100
model_u_cmip6_frac      = (model_u_cmip6/total_u_cmip6)*100
scen_u_cmip6_frac       = (scen_u_cmip6/total_u_cmip6)*100
copy_VarCoords(model_u_cmip6,int_u_mean_cmip6_frac)
copy_VarCoords(model_u_cmip6,model_u_cmip6_frac)
copy_VarCoords(model_u_cmip6,scen_u_cmip6_frac)








;========================================================================================
  print_clock("plotting...")

if plot1 .eq. "yes"

do t = 0,dimsizes(tags)-1
  tag = tags(t)
  ; -- which uncertainties to plot? --
  if tag .eq. "le"
    int_u_frac = int_u_mean_le_frac
    model_u_frac = model_u_le_frac
    scen_u_frac = scen_u_le_frac
  else if tag .eq. "le_hs09"
    int_u_frac = int_u_mean_le_hs09_frac
    model_u_frac = model_u_le_hs09_frac
    scen_u_frac = scen_u_le_hs09_frac
  else if tag .eq. "cmip5"
    int_u_frac = int_u_mean_cmip5_frac
    model_u_frac = model_u_cmip5_frac
    scen_u_frac = scen_u_cmip5_frac
  else if tag .eq. "cmip6"
    int_u_frac = int_u_mean_cmip6_frac
    model_u_frac = model_u_cmip6_frac
    scen_u_frac = scen_u_cmip6_frac
  end if
  end if
  end if
  end if

  wks_type = "eps" ; eps png
  if (wks_type.eq."png") then
     wks_type@wkWidth = 1500
     wks_type@wkHeight = 1500
  end if
  ofile = "/home/flehner/publication/lehner19_revisiting_hawkins_sutton/fig/"+get_script_prefix_name()+"_"+var+"_"+seas+"_"+tag+"_"+grid+"_ref"+refstart+"_"+refende+ending
  system("rm -rf "+ofile+"*")
  wks = gsn_open_wks(wks_type,ofile)

  res = True
  proj = "CylindricalEquidistant" ; Robinson CylindricalEquidistant
  res@mpProjection      = proj
  ; -- user-defined aspect ratio (only used with CylindricalEquidistant):
  if proj .eq. "CylindricalEquidistant"
    res@mpShapeMode       = "FreeAspect"
    res@vpWidthF          = 0.45
    res@vpHeightF         = 0.4
    res@tmXBOn            = False
    res@tmYLOn            = False
    res@tmXTOn            = False
    res@tmYROn            = False
  end if

  res@mpFillOn = True
  res@mpFillColors = (/-1,0,-1,-1/)
  res@cnFillMode    = "RasterFill"; RasterFill CellFill
  res@gsnAddCyclic         = True
  res@mpGeophysicalLineColor = "gray35"
  res@mpUSStateLineColor = res@mpGeophysicalLineColor
  res@mpNationalLineColor = "gray20"
  res@mpPerimOn    = False
  res@mpGridLatSpacingF =  90            ; change latitude  line spacing
  res@mpGridLonSpacingF = 180.           ; change longitude line spacing
  res@mpGridLineColor   = "transparent"  ; trick ncl into drawing perimeter
  res@mpGridAndLimbOn   = True           ; turn on lat/lon lines
  res@mpOutlineOn = True
  res@mpDataSetName = "Earth..4"
  res@mpGeophysicalLineThicknessF = 1.5
  res@mpPerimDrawOrder = "PostDraw"
  res@mpOutlineDrawOrder = "Draw"
  res@mpFillDrawOrder = "PostDraw"
  res@cnFillDrawOrder = "PreDraw"

  if grid .eq. "gea"
    lati = new(11,typeof(lat_bnds1))
    lati(:9) = lat_bnds1(:,0)
    lati(10) = lat_bnds1(9,1)
    loni = new(19,typeof(lon_bnds1))
    loni(:17) = lon_bnds1(:,0)
    loni(18) = lon_bnds1(17,1)
    ; res@sfXArray = loni
    res@sfYArray = lati
  end if

  if (domain.eq."NA") then
    res@mpDataBaseVersion     = "MediumRes"
    res@mpOutlineBoundarySets = "GeophysicalAndUSStates"
    res@mpOutlineSpecifiers = (/"Canada : Provinces","Mexico : States"/)
    res@mpLimitMode = "LatLon"
    res@mpMinLonF = 200.
    res@mpMaxLonF = 298.
    res@mpCenterLonF = (res@mpMinLonF+res@mpMaxLonF)/2.0
    res@mpMinLatF = 24.
    res@mpMaxLatF = 72.
    res@mpPerimOn    = True
    res@gsnLeftStringOrthogonalPosF = -.96
    res@gsnLeftStringParallelPosF = 0.025
    res@gsnLeftStringFontHeightF = 0.022
    res@gsnLeftStringFontColor = "gray25"
  else
    if var .eq. "tas"
      res@mpCenterLonF = 0. ; 0. 210.
    end if
    if var .eq. "pr"
      res@mpCenterLonF = 0.;210.
    end if
    res@mpFillOn = False
    res@gsnCenterStringFontHeightF = 0.03
    res@gsnCenterStringFontColor = "gray25"
    tw = "gsnCenterString"
  end if

  res@gsnFrame     =  False
  res@gsnDraw = False
  res@cnLineLabelsOn = False
  res@cnFillOn        = True
  res@cnLinesOn       = False
  res@tiMainOn = False
  res@cnLevelSelectionMode = "ExplicitLevels"
  res@lbLabelBarOn = False
  res@tiXAxisSide = "Top"


  panres = True
  panres@gsnMaximize = False
  panres@gsnPaperOrientation = "portrait"
  panres@gsnPanelLabelBar = True
  panres@lbTitleString = "Fractional uncertainty (%)"
  ; panres@gsnPanelTop = .90
  panres@gsnPanelBottom = .05
  panres@pmLabelBarWidthF = 0.62
  panres@pmLabelBarHeightF = 0.062
  panres@lbBoxEndCapStyle = "TriangleBothEnds"
  panres@lbBoxLineColor = "gray30"
  panres@lbTitlePosition = "bottom"
  panres@lbTitleOn = True

  panres@lbLabelFontHeightF = 0.02;145
  panres@lbTitleFontHeightF = panres@lbLabelFontHeightF
  panres@pmLabelBarOrthogonalPosF = -0.02;13


  res@cnLevels = (/10,20,30,40,50,60,70,80,90/)
  ; gsn_define_colormap(wks,"BlueDarkRed18")
  ; res@cnFillColors = (/2,4,5,7,9,12,13,14,16,18/)
  hs09_cols = new((/12,3/),float)
  hs09_cols(2,:) = (/25.,25.,152./)/255
  hs09_cols(3,:) = (/79.,86.,184./)/255
  hs09_cols(4,:) = (/117,146.,207./)/255
  hs09_cols(5,:) = (/188.,216.,234./)/255
  hs09_cols(6,:) = (/233.,241.,249./)/255
  hs09_cols(7,:) = (/249.,224.,150./)/255
  hs09_cols(8,:) = (/234.,153.,87./)/255
  hs09_cols(9,:) = (/226.,84.,45./)/255
  hs09_cols(10,:) = (/206.,37.,22./)/255
  hs09_cols(11,:) = (/154.,21.,20./)/255
  hs09_cols(0,:) = (/1,1,1/)
  hs09_cols(1,:) = (/0,0,0/)
  gsn_define_colormap(wks,hs09_cols)
  res@cnFillColors = (/2,3,4,5,6,7,8,9,10,11/)

  if grid .eq. "gea"
    gw = 1
  else
    lat = int_u_frac&lat
    rad = 4.*atan(1.)/180.
    gw  = sqrt( cos(rad*lat) )   ; gw for gaussian grid
  end if

  res@gsnRightStringOrthogonalPosF = -0.95
  res@gsnRightStringParallelPosF = 0.94
  res@gsnRightStringFontHeightF = 0.03
  res@tiXAxisFontHeightF = 0.005;0.03
  res@tiYAxisFontHeightF = 0.005;0.03

  ; ----------------------------------------------------------------------------

  nplot   = 9
  plot    = new(nplot,graphic)
  labels  = new((/nplot/),string)
  decades = (/"1st decade","4th decade","8th decade"/)
  if refstart .eq. 1971
    ; yyyy = (/2004,2034,2074/); 1st, 4th, 8th decade
    yyyy = (/2004,2034,2084/); 1st, 4th, 9th decade
  else if refstart .eq. 1991
    yyyy = (/2024,2054,2094/); 1st, 4th, 8th decade
  else if refstart .eq. 1995
    yyyy = (/2019,2049,2089/); 1st, 4th, 8th decade
  end if
  end if
  end if
  do i = 0,dimsizes(yyyy)-1
    res@tiYAxisString = "" ; decades(i)+" ("+(yyyy(i)-4)+"s)"
    ; ------
    if i .eq. 0
      res@tiXAxisString = "" ; "Internal variability"
    else
      res@tiXAxisString = ""
    end if
    frac = decimalPlaces(wgt_areaave(int_u_frac(yyyy(i)-start_trend,:,:),gw,1,1),0,True)
    labels(i*3) = frac+"%"
    plot(i*3) = gsn_csm_contour_map(wks,int_u_frac(yyyy(i)-start_trend,:,:),res)
    ; ------
    res@tiYAxisString = ""
    if i .eq. 0
      res@tiXAxisString = "" ; "Model uncertainty"
    else
      res@tiXAxisString = ""
    end if
    frac = decimalPlaces(wgt_areaave(model_u_frac(yyyy(i)-start_trend,:,:),gw,1,1),0,True)
    labels(i*3+1) = frac+"%"
    plot(i*3+1) = gsn_csm_contour_map(wks,model_u_frac(yyyy(i)-start_trend,:,:),res)
    ; ------
    if i .eq. 0
      res@tiXAxisString = "" ; "Scenario uncertainty"
    else
      res@tiXAxisString = ""
    end if
    frac = decimalPlaces(wgt_areaave(scen_u_frac(yyyy(i)-start_trend,:,:),gw,1,1),0,True)
    labels(i*3+2) = frac+"%"
    plot(i*3+2) = gsn_csm_contour_map(wks,scen_u_frac(yyyy(i)-start_trend,:,:),res)
  end do
  ; ----------------------------------------------------------------------------

  panres@gsnPanelFigureStrings = labels
  panres@gsnPanelFigureStringsFontHeightF = 0.02
  panres@gsnPanelFigureStringsPerimOn = False
  panres@amJust          = "BottomRight"

  panres@gsnPanelMainString = "" ; "Fractional uncertainty "+var+" "+seas
  panres@gsnPanelXWhiteSpacePercent = 0.
  panres@gsnPanelYWhiteSpacePercent = 0.

  gsn_panel(wks,plot,(/3,3/),panres)
  delete([/plot,decades,frac,labels/])
  delete(wks)

end do ; end of tags loop (cmip5, le, le_hs09)

end if ; end of plot1 loop





; ==============================================================================





if plot2 .eq. "yes"

  decades = (/"decade1","decade4"/)
  do d = 0,dimsizes(decades)-1
    if refstart .eq. 1971 .and. d .eq. 0
      yyyy = 2004
    else if refstart .eq. 1971 .and. d .eq. 1
      yyyy = 2034
    else if refstart .eq. 1991 .and. d .eq. 0
      yyyy = 2024
    else if refstart .eq. 1991 .and. d .eq. 1
      yyyy = 2054
    else if refstart .eq. 1995 .and. d .eq. 0
      yyyy = 2019
    else if refstart .eq. 1995 .and. d .eq. 1
      yyyy = 2049
    end if
    end if
    end if
    end if
    end if
    end if

    wks_type = "eps" ; eps png
    if (wks_type.eq."png") then
       wks_type@wkWidth = 1500
       wks_type@wkHeight = 1500
    end if
    if constrain .eq. "no"
      ofile = "/home/flehner/publication/lehner19_revisiting_hawkins_sutton/fig/"+get_script_prefix_name()+"_"+var+"_"+seas+"_"+grid+"_ref"+refstart+"_"+refende+"_"+decades(d)
    else
      ofile = "/home/flehner/publication/lehner19_revisiting_hawkins_sutton/fig/"+get_script_prefix_name()+"_"+var+"_"+seas+"_"+grid+"_ref"+refstart+"_"+refende+"_"+decades(d)+"_wo_highsens_models"
    end if
    system("rm -rf "+ofile+"*")
    wks = gsn_open_wks(wks_type,ofile)

    res = True
    proj = "CylindricalEquidistant" ; Robinson CylindricalEquidistant
    res@mpProjection      = proj
    ; -- user-defined aspect ratio (only used with CylindricalEquidistant):
    if proj .eq. "CylindricalEquidistant"
      res@mpShapeMode       = "FreeAspect"
      res@vpWidthF          = 0.45
      res@vpHeightF         = 0.4
      res@tmXBOn            = False
      res@tmYLOn            = False
      res@tmXTOn            = False
      res@tmYROn            = False
    end if

    res@mpFillOn = True
    res@mpFillColors = (/-1,0,-1,-1/)
    res@cnFillMode    = "RasterFill"; RasterFill CellFill
    res@gsnAddCyclic         = True
    res@mpGeophysicalLineColor = "gray35"
    res@mpUSStateLineColor = res@mpGeophysicalLineColor
    res@mpNationalLineColor = "gray20"
    res@mpPerimOn    = False
    res@mpGridLatSpacingF =  90            ; change latitude  line spacing
    res@mpGridLonSpacingF = 180.           ; change longitude line spacing
    res@mpGridLineColor   = "transparent"  ; trick ncl into drawing perimeter
    res@mpGridAndLimbOn   = True           ; turn on lat/lon lines
    res@mpOutlineOn = True
    res@mpDataSetName = "Earth..4"
    res@mpGeophysicalLineThicknessF = 1.5
    res@mpPerimDrawOrder = "PostDraw"
    res@mpOutlineDrawOrder = "Draw"
    res@mpFillDrawOrder = "PostDraw"
    res@cnFillDrawOrder = "PreDraw"

    if grid .eq. "gea"
      lati = new(11,typeof(lat_bnds1))
      lati(:9) = lat_bnds1(:,0)
      lati(10) = lat_bnds1(9,1)
      loni = new(19,typeof(lon_bnds1))
      loni(:17) = lon_bnds1(:,0)
      loni(18) = lon_bnds1(17,1)
      ; res@sfXArray = loni
      res@sfYArray = lati
    end if

    if (domain.eq."NA") then
      res@mpDataBaseVersion     = "MediumRes"
      res@mpOutlineBoundarySets = "GeophysicalAndUSStates"
      res@mpOutlineSpecifiers = (/"Canada : Provinces","Mexico : States"/)
      res@mpLimitMode = "LatLon"
      res@mpMinLonF = 200.
      res@mpMaxLonF = 298.
      res@mpCenterLonF = (res@mpMinLonF+res@mpMaxLonF)/2.0
      res@mpMinLatF = 24.
      res@mpMaxLatF = 72.
      res@mpPerimOn    = True
      res@gsnLeftStringOrthogonalPosF = -.96
      res@gsnLeftStringParallelPosF = 0.025
      res@gsnLeftStringFontHeightF = 0.022
      res@gsnLeftStringFontColor = "gray25"
    else
      if var .eq. "tas"
        res@mpCenterLonF = 0. ; 0. 210.
      end if
      if var .eq. "pr"
        res@mpCenterLonF = 0.;210.
      end if
      res@mpFillOn = False
      res@gsnCenterStringFontHeightF = 0.03
      res@gsnCenterStringFontColor = "gray25"
      tw = "gsnCenterString"
    end if

    res@gsnFrame     =  False
    res@gsnDraw = False
    res@cnLineLabelsOn = False
    res@cnFillOn        = True
    res@cnLinesOn       = False
    res@tiMainOn = False
    res@cnLevelSelectionMode = "ExplicitLevels"
    res@lbLabelBarOn = False
    res@tiXAxisSide = "Top"


    panres = True
    panres@gsnMaximize = False
    panres@gsnPaperOrientation = "portrait"
    panres@gsnPanelLabelBar = True
    panres@lbTitleString = "Fractional uncertainty (%)"
    ; panres@gsnPanelTop = .90
    panres@gsnPanelBottom = .05
    panres@pmLabelBarWidthF = 0.62
    panres@pmLabelBarHeightF = 0.062
    panres@lbBoxEndCapStyle = "TriangleBothEnds"
    panres@lbBoxLineColor = "gray30"
    panres@lbTitlePosition = "bottom"
    panres@lbTitleOn = True

    panres@lbLabelFontHeightF = 0.02;145
    panres@lbTitleFontHeightF = panres@lbLabelFontHeightF
    panres@pmLabelBarOrthogonalPosF = -0.02;13

    res@cnLevels = (/10,20,30,40,50,60,70,80,90/)
    ; gsn_define_colormap(wks,"BlueDarkRed18")
    ; res@cnFillColors = (/2,4,5,7,9,12,13,14,16,18/)
    hs09_cols = new((/12,3/),float)
    hs09_cols(2,:) = (/25.,25.,152./)/255
    hs09_cols(3,:) = (/79.,86.,184./)/255
    hs09_cols(4,:) = (/117,146.,207./)/255
    hs09_cols(5,:) = (/188.,216.,234./)/255
    hs09_cols(6,:) = (/233.,241.,249./)/255
    hs09_cols(7,:) = (/249.,224.,150./)/255
    hs09_cols(8,:) = (/234.,153.,87./)/255
    hs09_cols(9,:) = (/226.,84.,45./)/255
    hs09_cols(10,:) = (/206.,37.,22./)/255
    hs09_cols(11,:) = (/154.,21.,20./)/255
    hs09_cols(0,:) = (/1,1,1/)
    hs09_cols(1,:) = (/0,0,0/)
    gsn_define_colormap(wks,hs09_cols)
    res@cnFillColors = (/2,3,4,5,6,7,8,9,10,11/)

    if grid .eq. "gea"
      gw = 1
    else
      lat = int_u_mean_le_frac&lat
      rad = 4.*atan(1.)/180.
      gw  = sqrt( cos(rad*lat) )   ; gw for gaussian grid
    end if

    res@gsnRightStringOrthogonalPosF = -0.95
    res@gsnRightStringParallelPosF = 0.94
    res@gsnRightStringFontHeightF = 0.03
    res@tiXAxisFontHeightF = 0.005;0.03
    res@tiYAxisFontHeightF = 0.005;0.03

    ; ----------------------------------------------------------------------------

    nplot   = 9
    plot    = new(nplot,graphic)
    labels  = new((/nplot/),string)
    do i = 0,2
      if i .eq. 2
        int_u_frac    = int_u_mean_le_frac
        model_u_frac  = model_u_le_frac
        scen_u_frac   = scen_u_le_frac
      else if i .eq. 3
        int_u_frac    = int_u_mean_le_hs09_frac
        model_u_frac  = model_u_le_hs09_frac
        scen_u_frac   = scen_u_le_hs09_frac
      else if i .eq. 0
        int_u_frac    = int_u_mean_cmip5_frac
        model_u_frac  = model_u_cmip5_frac
        scen_u_frac   = scen_u_cmip5_frac
      else if i .eq. 1
        int_u_frac    = int_u_mean_cmip6_frac
        model_u_frac  = model_u_cmip6_frac
        scen_u_frac   = scen_u_cmip6_frac
      end if
      end if
      end if
      end if
      res@tiYAxisString = ""
      ; ------
      if i .eq. 0
        res@tiXAxisString = "" ; "Internal variability"
      else
        res@tiXAxisString = ""
      end if
      frac        = decimalPlaces(wgt_areaave(int_u_frac(yyyy-start_trend,:,:),gw,1,1),0,True)
      labels(i*3) = frac+"%"
      plot(i*3)   = gsn_csm_contour_map(wks,int_u_frac(yyyy-start_trend,:,:),res)
      ; ------
      res@tiYAxisString = ""
      if i .eq. 0
        res@tiXAxisString = "" ; "Model uncertainty"
      else
        res@tiXAxisString = ""
      end if
      frac = decimalPlaces(wgt_areaave(model_u_frac(yyyy-start_trend,:,:),gw,1,1),0,True)
      labels(i*3+1) = frac+"%"
      plot(i*3+1) = gsn_csm_contour_map(wks,model_u_frac(yyyy-start_trend,:,:),res)
      ; ------
      if i .eq. 0
        res@tiXAxisString = "" ; "Scenario uncertainty"
      else
        res@tiXAxisString = ""
      end if
      frac = decimalPlaces(wgt_areaave(scen_u_frac(yyyy-start_trend,:,:),gw,1,1),0,True)
      labels(i*3+2) = frac+"%"
      plot(i*3+2) = gsn_csm_contour_map(wks,scen_u_frac(yyyy-start_trend,:,:),res)
    end do
    ; ----------------------------------------------------------------------------

    panres@gsnPanelFigureStrings = labels
    panres@gsnPanelFigureStringsFontHeightF = 0.02
    panres@gsnPanelFigureStringsPerimOn = False
    panres@amJust          = "BottomRight"

    panres@gsnPanelMainString = "" ; "Fractional uncertainty "+var+" "+seas
    panres@gsnPanelXWhiteSpacePercent = 0.
    panres@gsnPanelYWhiteSpacePercent = 0.

    gsn_panel(wks,plot,(/3,3/),panres)
    delete([/plot/])
    delete(wks)

  end do ; end of decades loop

end if ; end of plot2 loop





if plot3 .eq. "yes"

do t = 0,1
  if t .eq. 0
    tag_tmp     = "le"
    int_diff    = int_u_mean_le_hs09_frac-int_u_mean_le_frac
    model_diff  = model_u_le_hs09_frac-model_u_le_frac
    scen_diff   = scen_u_le_hs09_frac-scen_u_le_frac
  else
    tag_tmp     = "cmip5"
    int_diff    = int_u_mean_cmip5_frac-int_u_mean_le_frac
    model_diff  = model_u_cmip5_frac-model_u_le_frac
    scen_diff   = scen_u_cmip5_frac-scen_u_le_frac
  end if
  copy_VarCoords(int_u_mean_le_frac,int_diff)
  copy_VarCoords(int_u_mean_le_frac,model_diff)
  copy_VarCoords(int_u_mean_le_frac,scen_diff)

  wks_type = "eps" ; eps png
  if (wks_type.eq."png") then
     wks_type@wkWidth = 1500
     wks_type@wkHeight = 1500
  end if
  ofile = "/home/flehner/publication/lehner19_revisiting_hawkins_sutton/fig/"+get_script_prefix_name()+"_"+var+"_"+seas+"_le_"+tag_tmp+"_hs09_diff_"+grid+"_ref"+refstart+"_"+refende
  system("rm -rf "+ofile+"*")
  wks = gsn_open_wks(wks_type,ofile)

  res = True
  proj = "CylindricalEquidistant" ; Robinson CylindricalEquidistant
  res@mpProjection      = proj
  ; -- user-defined aspect ratio (only used with CylindricalEquidistant):
  if proj .eq. "CylindricalEquidistant"
    res@mpShapeMode       = "FreeAspect"
    res@vpWidthF          = 0.45
    res@vpHeightF         = 0.4
    res@tmXBOn            = False
    res@tmYLOn            = False
    res@tmXTOn            = False
    res@tmYROn            = False
  end if

  res@mpFillOn = True
  res@mpFillColors = (/-1,0,-1,-1/)
  res@cnFillMode    = "RasterFill"; RasterFill CellFill
  res@gsnAddCyclic         = True
  res@mpGeophysicalLineColor = "gray35"
  res@mpUSStateLineColor = res@mpGeophysicalLineColor
  res@mpNationalLineColor = "gray20"
  res@mpPerimOn    = False
  res@mpGridLatSpacingF =  90            ; change latitude  line spacing
  res@mpGridLonSpacingF = 180.           ; change longitude line spacing
  res@mpGridLineColor   = "transparent"  ; trick ncl into drawing perimeter
  res@mpGridAndLimbOn   = True           ; turn on lat/lon lines
  res@mpOutlineOn = True
  res@mpDataSetName = "Earth..4"
  res@mpGeophysicalLineThicknessF = 1.5
  res@mpPerimDrawOrder = "PostDraw"
  res@mpOutlineDrawOrder = "Draw"
  res@mpFillDrawOrder = "PostDraw"
  res@cnFillDrawOrder = "PreDraw"

  if grid .eq. "gea"
    lati = new(11,typeof(lat_bnds1))
    lati(:9) = lat_bnds1(:,0)
    lati(10) = lat_bnds1(9,1)
    loni = new(19,typeof(lon_bnds1))
    loni(:17) = lon_bnds1(:,0)
    loni(18) = lon_bnds1(17,1)
    ; res@sfXArray = loni
    res@sfYArray = lati
  end if

  if (domain.eq."NA") then
    res@mpDataBaseVersion     = "MediumRes"
    res@mpOutlineBoundarySets = "GeophysicalAndUSStates"
    res@mpOutlineSpecifiers = (/"Canada : Provinces","Mexico : States"/)
    res@mpLimitMode = "LatLon"
    res@mpMinLonF = 200.
    res@mpMaxLonF = 298.
    res@mpCenterLonF = (res@mpMinLonF+res@mpMaxLonF)/2.0
    res@mpMinLatF = 24.
    res@mpMaxLatF = 72.
    res@mpPerimOn    = True
    res@gsnLeftStringOrthogonalPosF = -.96
    res@gsnLeftStringParallelPosF = 0.025
    res@gsnLeftStringFontHeightF = 0.022
    res@gsnLeftStringFontColor = "gray25"
  else
    if var .eq. "tas"
      res@mpCenterLonF = 0. ; 0. 210.
    end if
    if var .eq. "pr"
      res@mpCenterLonF = 0.;210.
    end if
    res@mpFillOn = False
    res@gsnCenterStringFontHeightF = 0.03
    res@gsnCenterStringFontColor = "gray25"
    tw = "gsnCenterString"
  end if

  res@gsnFrame     =  False
  res@gsnDraw = False
  res@cnLineLabelsOn = False
  res@cnFillOn        = True
  res@cnLinesOn       = False
  res@tiMainOn = False
  res@cnLevelSelectionMode = "ExplicitLevels"
  res@lbLabelBarOn = False
  res@tiXAxisSide = "Top"


  panres = True
  panres@gsnMaximize = False
  panres@gsnPaperOrientation = "portrait"
  panres@gsnPanelLabelBar = True
  panres@lbTitleString = "Fractional uncertainy (%)";"Method error (percentage points)"
  panres@gsnPanelBottom = .05
  panres@pmLabelBarWidthF = 0.62
  panres@pmLabelBarHeightF = 0.062
  panres@lbBoxEndCapStyle = "TriangleBothEnds"
  panres@lbBoxLineColor = "gray30"
  panres@lbTitlePosition = "bottom"
  panres@lbTitleOn = True

  panres@lbLabelFontHeightF = 0.02;145
  panres@lbTitleFontHeightF = panres@lbLabelFontHeightF
  panres@pmLabelBarOrthogonalPosF = -0.02;13


  ; res@cnLevels = (/10,20,30,40,50,60,70,80,90/)
  res@cnLevels = (/-10,-7.5,-5,-2.5,0,2.5,5,7.5,10/)*2
  ; gsn_define_colormap(wks,"BlueDarkRed18")
  ; res@cnFillColors = (/2,4,5,7,9,12,13,14,16,18/)
  hs09_cols = new((/12,3/),float)
  hs09_cols(2,:) = (/25.,25.,152./)/255
  hs09_cols(3,:) = (/79.,86.,184./)/255
  hs09_cols(4,:) = (/117,146.,207./)/255
  hs09_cols(5,:) = (/188.,216.,234./)/255
  hs09_cols(6,:) = (/233.,241.,249./)/255
  hs09_cols(7,:) = (/249.,224.,150./)/255
  hs09_cols(8,:) = (/234.,153.,87./)/255
  hs09_cols(9,:) = (/226.,84.,45./)/255
  hs09_cols(10,:) = (/206.,37.,22./)/255
  hs09_cols(11,:) = (/154.,21.,20./)/255
  hs09_cols(0,:) = (/1,1,1/)
  hs09_cols(1,:) = (/0,0,0/)
  gsn_define_colormap(wks,hs09_cols)
  res@cnFillColors = (/2,3,4,5,6,7,8,9,10,11/)

  if grid .eq. "gea"
    gw = 1
  else
    lat = int_diff&lat
    rad = 4.*atan(1.)/180.
    gw  = sqrt( cos(rad*lat) )   ; gw for gaussian grid
  end if

  res@gsnRightStringOrthogonalPosF = -0.95
  res@gsnRightStringParallelPosF = 0.94
  res@gsnRightStringFontHeightF = 0.03
  res@tiXAxisFontHeightF = 0.005;0.03
  res@tiYAxisFontHeightF = 0.005;0.03

  ; ----------------------------------------------------------------------------

  nplot   = 9
  plot    = new(nplot,graphic)
  labels  = new((/nplot/),string)
  decades = (/"1st decade","4th decade","8th decade"/)
  if refstart .eq. 1971
    ; yyyy = (/2004,2034,2074/); 1st, 4th, 8th decade
    yyyy = (/2004,2034,2084/); 1st, 4th, 9th decade
  else if refstart .eq. 1991
    yyyy = (/2024,2054,2094/); 1st, 4th, 8th decade
  else if refstart .eq. 1995
    yyyy = (/2019,2049,2089/); 1st, 4th, 8th decade
  end if
  end if
  end if
  do i = 0,dimsizes(yyyy)-1
    res@tiYAxisString = "" ; decades(i)+" ("+(yyyy(i)-4)+"s)"
    ; ------
    if i .eq. 0
      res@tiXAxisString = "" ; "Internal variability"
    else
      res@tiXAxisString = ""
    end if
    frac = decimalPlaces(wgt_areaave(int_diff(yyyy(i)-start_trend,:,:),gw,1,1),0,True)
    labels(i*3) = frac+"%"
    plot(i*3) = gsn_csm_contour_map(wks,int_diff(yyyy(i)-start_trend,:,:),res)
    ; ------
    res@tiYAxisString = ""
    if i .eq. 0
      res@tiXAxisString = "" ; "Model uncertainty"
    else
      res@tiXAxisString = ""
    end if
    frac = decimalPlaces(wgt_areaave(model_diff(yyyy(i)-start_trend,:,:),gw,1,1),0,True)
    labels(i*3+1) = frac+"%"
    plot(i*3+1) = gsn_csm_contour_map(wks,model_diff(yyyy(i)-start_trend,:,:),res)
    ; ------
    if i .eq. 0
      res@tiXAxisString = "" ; "Scenario uncertainty"
    else
      res@tiXAxisString = ""
    end if
    frac = decimalPlaces(wgt_areaave(scen_diff(yyyy(i)-start_trend,:,:),gw,1,1),0,True)
    labels(i*3+2) = frac+"%"
    plot(i*3+2) = gsn_csm_contour_map(wks,scen_diff(yyyy(i)-start_trend,:,:),res)
  end do
  ; ----------------------------------------------------------------------------

  panres@gsnPanelFigureStrings = labels
  panres@gsnPanelFigureStringsFontHeightF = 0.02
  panres@gsnPanelFigureStringsPerimOn = False
  panres@amJust          = "BottomRight"

  panres@gsnPanelMainString = "" ; "Fractional uncertainty "+var+" "+seas
  panres@gsnPanelXWhiteSpacePercent = 0.
  panres@gsnPanelYWhiteSpacePercent = 0.

  gsn_panel(wks,plot,(/3,3/),panres)
  delete([/plot,decades,frac,labels/])
  delete(wks)

end do ; end of tag loop (le_hs09 and cmip5)

end if ; end of plot3 loop





if plot4 .eq. "yes"

  ; xxx
  int_u_frac = int_u_mean_le_hs09_all_frac
  model_u_frac = model_u_le_hs09_all_frac
  scen_u_frac = scen_u_le_hs09_all_frac

  wks_type = "eps" ; eps png
  if (wks_type.eq."png") then
     wks_type@wkWidth = 1500
     wks_type@wkHeight = 1500
  end if
  ofile = "/home/flehner/publication/lehner19_revisiting_hawkins_sutton/fig/"+get_script_prefix_name()+"_"+var+"_"+seas+"_le_method_error_"+grid+"_ref"+refstart+"_"+refende
  system("rm -rf "+ofile+"*")
  wks = gsn_open_wks(wks_type,ofile)

  res = True
  proj = "CylindricalEquidistant" ; Robinson CylindricalEquidistant
  res@mpProjection      = proj
  ; -- user-defined aspect ratio (only used with CylindricalEquidistant):
  if proj .eq. "CylindricalEquidistant"
    res@mpShapeMode       = "FreeAspect"
    res@vpWidthF          = 0.45
    res@vpHeightF         = 0.4
    res@tmXBOn            = False
    res@tmYLOn            = False
    res@tmXTOn            = False
    res@tmYROn            = False
  end if

  res@mpFillOn = True
  res@mpFillColors = (/-1,0,-1,-1/)
  res@cnFillMode    = "RasterFill"; RasterFill CellFill
  res@gsnAddCyclic         = True
  res@mpGeophysicalLineColor = "gray35"
  res@mpUSStateLineColor = res@mpGeophysicalLineColor
  res@mpNationalLineColor = "gray20"
  res@mpPerimOn    = False
  res@mpGridLatSpacingF =  90            ; change latitude  line spacing
  res@mpGridLonSpacingF = 180.           ; change longitude line spacing
  res@mpGridLineColor   = "transparent"  ; trick ncl into drawing perimeter
  res@mpGridAndLimbOn   = True           ; turn on lat/lon lines
  res@mpOutlineOn = True
  res@mpDataSetName = "Earth..4"
  res@mpGeophysicalLineThicknessF = 1.5
  res@mpPerimDrawOrder = "PostDraw"
  res@mpOutlineDrawOrder = "Draw"
  res@mpFillDrawOrder = "PostDraw"
  res@cnFillDrawOrder = "PreDraw"

  if grid .eq. "gea"
    lati = new(11,typeof(lat_bnds1))
    lati(:9) = lat_bnds1(:,0)
    lati(10) = lat_bnds1(9,1)
    loni = new(19,typeof(lon_bnds1))
    loni(:17) = lon_bnds1(:,0)
    loni(18) = lon_bnds1(17,1)
    ; res@sfXArray = loni
    res@sfYArray = lati
  end if

  if (domain.eq."NA") then
    res@mpDataBaseVersion     = "MediumRes"
    res@mpOutlineBoundarySets = "GeophysicalAndUSStates"
    res@mpOutlineSpecifiers = (/"Canada : Provinces","Mexico : States"/)
    res@mpLimitMode = "LatLon"
    res@mpMinLonF = 200.
    res@mpMaxLonF = 298.
    res@mpCenterLonF = (res@mpMinLonF+res@mpMaxLonF)/2.0
    res@mpMinLatF = 24.
    res@mpMaxLatF = 72.
    res@mpPerimOn    = True
    res@gsnLeftStringOrthogonalPosF = -.96
    res@gsnLeftStringParallelPosF = 0.025
    res@gsnLeftStringFontHeightF = 0.022
    res@gsnLeftStringFontColor = "gray25"
  else
    if var .eq. "tas"
      res@mpCenterLonF = 0. ; 0. 210.
    end if
    if var .eq. "pr"
      res@mpCenterLonF = 0.;210.
    end if
    res@mpFillOn = False
    res@gsnCenterStringFontHeightF = 0.03
    res@gsnCenterStringFontColor = "gray25"
    tw = "gsnCenterString"
  end if

  res@gsnFrame     =  False
  res@gsnDraw = False
  res@cnLineLabelsOn = False
  res@cnFillOn        = True
  res@cnLinesOn       = False
  res@tiMainOn = False
  res@cnLevelSelectionMode = "ExplicitLevels"
  res@lbLabelBarOn = False
  res@tiXAxisSide = "Top"


  panres = True
  panres@gsnMaximize = False
  panres@gsnPaperOrientation = "portrait"
  panres@gsnPanelLabelBar = True
  panres@lbTitleString = "Fractional uncertainy (%)";"Method error (percentage points)"
  panres@gsnPanelBottom = .05
  panres@pmLabelBarWidthF = 0.62
  panres@pmLabelBarHeightF = 0.062
  panres@lbBoxEndCapStyle = "TriangleBothEnds"
  panres@lbBoxLineColor = "gray30"
  panres@lbTitlePosition = "bottom"
  panres@lbTitleOn = True

  panres@lbLabelFontHeightF = 0.02;145
  panres@lbTitleFontHeightF = panres@lbLabelFontHeightF
  panres@pmLabelBarOrthogonalPosF = -0.02;13

  res@cnLevels = (/10,20,30,40,50,60,70,80,90/)
  ; gsn_define_colormap(wks,"BlueDarkRed18")
  ; res@cnFillColors = (/2,4,5,7,9,12,13,14,16,18/)
  hs09_cols = new((/12,3/),float)
  hs09_cols(2,:) = (/25.,25.,152./)/255
  hs09_cols(3,:) = (/79.,86.,184./)/255
  hs09_cols(4,:) = (/117,146.,207./)/255
  hs09_cols(5,:) = (/188.,216.,234./)/255
  hs09_cols(6,:) = (/233.,241.,249./)/255
  hs09_cols(7,:) = (/249.,224.,150./)/255
  hs09_cols(8,:) = (/234.,153.,87./)/255
  hs09_cols(9,:) = (/226.,84.,45./)/255
  hs09_cols(10,:) = (/206.,37.,22./)/255
  hs09_cols(11,:) = (/154.,21.,20./)/255
  hs09_cols(0,:) = (/1,1,1/)
  hs09_cols(1,:) = (/0,0,0/)
  gsn_define_colormap(wks,hs09_cols)
  res@cnFillColors = (/2,3,4,5,6,7,8,9,10,11/)

  if grid .eq. "gea"
    gw = 1
  else
    lat = int_u_frac&lat
    rad = 4.*atan(1.)/180.
    gw  = sqrt( cos(rad*lat) )   ; gw for gaussian grid
  end if

  res@gsnRightStringOrthogonalPosF = -0.95
  res@gsnRightStringParallelPosF = 0.94
  res@gsnRightStringFontHeightF = 0.03
  res@tiXAxisFontHeightF = 0.005;0.03
  res@tiYAxisFontHeightF = 0.005;0.03

  ; ----------------------------------------------------------------------------

  nplot   = 9
  plot    = new(nplot,graphic)
  labels  = new((/nplot/),string)
  decades = (/"1st decade","4th decade","8th decade"/)
  if refstart .eq. 1971
    ; yyyy = (/2004,2034,2074/); 1st, 4th, 8th decade
    yyyy = (/2004,2034,2084/); 1st, 4th, 9th decade
  else if refstart .eq. 1991
    yyyy = (/2024,2054,2094/); 1st, 4th, 8th decade
  else if refstart .eq. 1995
    yyyy = (/2019,2049,2089/); 1st, 4th, 8th decade
  end if
  end if
  end if
  do i = 0,dimsizes(yyyy)-1
    res@tiYAxisString = "" ; decades(i)+" ("+(yyyy(i)-4)+"s)"
    ; ------
    if i .eq. 0
      res@tiXAxisString = "" ; "Internal variability"
    else
      res@tiXAxisString = ""
    end if
    frac = decimalPlaces(wgt_areaave(int_u_frac(yyyy(i)-start_trend,:,:),gw,1,1),0,True)
    labels(i*3) = frac+"%"
    ; res@gsnRightString = frac+"%"
    ; if proj .eq. "CylindricalEquidistant" .and. frac .lt. 15
    ;   res@gsnRightStringFontColor = "white"
    ; else
    ;   res@gsnRightStringFontColor = "black"
    ; end if
    ; ; - text string version (not worked out):
    ; txres                       = True
    ; txres@txPerimOn             = True
    ; txres@txBackgroundFillColor = "White"
    ; txres@txFontHeightF         = 0.03
    ; txid1 = gsn_create_text(wks, frac+"%", txres)
    plot(i*3) = gsn_csm_contour_map(wks,int_u_frac(yyyy(i)-start_trend,:,:),res)
    ; ------
    res@tiYAxisString = ""
    if i .eq. 0
      res@tiXAxisString = "" ; "Model uncertainty"
    else
      res@tiXAxisString = ""
    end if
    frac = decimalPlaces(wgt_areaave(model_u_frac(yyyy(i)-start_trend,:,:),gw,1,1),0,True)
    labels(i*3+1) = frac+"%"
    plot(i*3+1) = gsn_csm_contour_map(wks,model_u_frac(yyyy(i)-start_trend,:,:),res)
    ; ------
    if i .eq. 0
      res@tiXAxisString = "" ; "Scenario uncertainty"
    else
      res@tiXAxisString = ""
    end if
    frac = decimalPlaces(wgt_areaave(scen_u_frac(yyyy(i)-start_trend,:,:),gw,1,1),0,True)
    labels(i*3+2) = frac+"%"
    plot(i*3+2) = gsn_csm_contour_map(wks,scen_u_frac(yyyy(i)-start_trend,:,:),res)
  end do
  ; ----------------------------------------------------------------------------

  panres@gsnPanelFigureStrings = labels
  panres@gsnPanelFigureStringsFontHeightF = 0.02
  panres@gsnPanelFigureStringsPerimOn = False
  panres@amJust          = "BottomRight"

  panres@gsnPanelMainString = "" ; "Fractional uncertainty "+var+" "+seas
  panres@gsnPanelXWhiteSpacePercent = 0.
  panres@gsnPanelYWhiteSpacePercent = 0.

  gsn_panel(wks,plot,(/3,3/),panres)
  delete([/plot,decades,frac,labels/])
  delete(wks)

end if ; end of plot4 loop



end do ; end of season loop
delete([/cmip6_ensmean_ssp126_all,cmip6_ensmean_ssp245_all,cmip6_ensmean_ssp585_all/])
delete([/cmip6_intvar_ssp126_all,cmip6_intvar_ssp245_all,cmip6_intvar_ssp585_all/])
end do ; end of variable loop

print_clock("...done")

end
