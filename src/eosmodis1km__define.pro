;+
; EosModis1km Class
;
;   MODIS Level 1/2 Producs 1km resolution data class
;
;-----------------------------------
; Copyright and contact information
;  Copyright:
;   Name            (initials)  year       contact-information
;   Souichiro Hioki (SH)        2015-2017  s.hioki@tamu.edu
;   Justus Milligan (JM)        2015       justus.milligan@tamu.edu
;
;-----------------------------------
; Revision History
;  2015.Aug.21 SH coded
;  2015.Sep.07 SH adjusted for the new specification for chunk system
;  2015.Sep.21 SH migrated to Platyco
;  2015.Sep.25 SH backported getAngles from EosMIL1B2E
;  2015.Sep.26 SH backported error handling from EosMIL1B2E
;  2015.Oct.17 SH integrated into getMOD02ReflectivityProduct() function
;                 implemented getMOD06CloudPhaseInfrared
;                 implemented loading program for MOD06CloudPhaseInfrared
;  2015.Oct.18 SH implemented loading program for MOD02RadianceB31
;  2015.Oct.23 SH implemented inverse Planck function
;  2015.Oct.24 SH fixed bugs in infrared data processing
;                 implemented congridext() to interpolate geometric data
;  2015.Oct.27 SH modified loadPosGeom() to be compatible to both MOD02 and MOD06
;                 added getLatLonRange()
;  2015.Nov.01 SH fixed bugs in getLatLon()
;  2015.Nov.03 SH fixed bugs about date line treatment in loadPosGeom()
;  2015.Nov.04 SH improved MOD06 file treatment in loadMOD06()
;  2015.Nov.09 SH changed the name of varible from MOD06EffectiveRadius to MOD06CloudEffectiveRadius
;                 fixed bugs in __invPlanck()
;  2015.Nov.12 JM added writeImage() and writePng()
;  2015.Nov.17 JM added bands 1, 3, 4, and added visible mode to writepng()
;  2015.Dec.01 SH added getMOD02ReflectivityProductVariance()
;  2015.Dec.01 JM added band 31 to writepng()
;  2016.Apr.12 SH corrected the sign of relative azimuth angle
;  2016.May.16 SH added cloud emissivity reader and getMOD06CloudEmissivity()
;  2016.Jun.14 SH added cloud pressure color shading in writePng()
;  2016.Jun.21 SH added cloud mask SPI and reflectivity with atmospheric correction
;                 rewrote data extraction program for readibility
;  2016.Jun.25 SH improved accuracy of viewing angle interpolations in loadPosGeom()
;                 added getScatteringAngle() and getGlitterAngle()
;  2016.Jul.27 SH added latitude and longitude grid overlay in writePng()
;                 added plotReverse variable to correctly plot both Aqua and Terra MODIS
;  2016.Sep.09 SH added QualityAssurance1km
;  2017.Feb.23 SH added gamma correction to writeImage(), writePng()
;  2017.Mar.15 SH added getSolarAzimuthAngle() and getSpaceCraftAltitude()
;                 added a new variable: MOD03SpacecraftPosition, reader is not coded yet
;                 corrected a bug in solar azimuth direction
;  2017.Apr.02 SH added SoG2 and L5oG2 modes in getMOD02ReflectivityProduct(), getMOD06ReflectivityAC()
;  2017.May.09 SH added RadianceB29, RadianceB32, CloudTopHeightInfrared, CloudTopTemperatureInfrared
;
;-----------------------------------
; This program is free software: you can redistribute it and/or
; modify it as long as proper credits are given for all contributors.
; When you redistribute the program with modifications, please add
; your name and your initials to "Copyright and contact information"
; section and briefly explain the modification in the "Revision History"
; section with your initials.
;
; This program is disclosed in the hope that it will be useful,
; but WITHOUT ANY WARRANTY; without even the implied warranty of
; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
;-
;-----------------------------------
; constructor and destructor
function EosModis1km::init, filename
  void = self->IntfData::init()
  if n_elements(filename) eq 0 then filename =''
  self.filename                         = filename ; To use MOD03/MYD03 geolocation data, set this filename
  self.recordDimension                  = [0UL,0UL]
  self.plotReverse                      = 0 ; 0 when ascending, 1 when descending
  self.position                         = obj_new('EosPosition', self)
  self.grid                             = self.position
  self.geometry                         = obj_new('EosGeometry', self)
  self.MOD03SpacecraftPosition          = ptr_new(/allocate)
  self.MOD02ReflectivityProductB01      = ptr_new(/allocate)
  self.MOD02ReflectivityProductB02      = ptr_new(/allocate)
  self.MOD02ReflectivityProductB03      = ptr_new(/allocate)
  self.MOD02ReflectivityProductB04      = ptr_new(/allocate)
  self.MOD02ReflectivityProductB07      = ptr_new(/allocate)
  self.MOD02RadianceB29                 = ptr_new(/allocate)
  self.MOD02RadianceB31                 = ptr_new(/allocate)
  self.MOD02RadianceB32                 = ptr_new(/allocate)
  self.MOD06CloudPhaseInfrared          = ptr_new(/allocate)
  self.MOD06CloudTopPressureInfrared    = ptr_new(/allocate)
  self.MOD06CloudTopHeightInfrared      = ptr_new(/allocate)
  self.MOD06CloudTopTemperatureInfrared = ptr_new(/allocate)
  self.MOD06CloudPhaseOpticalProperties = ptr_new(/allocate)
  self.MOD06CloudOpticalThickness       = ptr_new(/allocate)
  self.MOD06CloudEffectiveRadius        = ptr_new(/allocate)
  self.MOD06CloudEmissivity11           = ptr_new(/allocate)
  self.MOD06CloudEmissivity12           = ptr_new(/allocate)
  self.MOD06CloudEmissivity13           = ptr_new(/allocate)
  self.MOD06CloudMaskSPIB01             = ptr_new(/allocate)
  self.MOD06CloudMaskSPIB02             = ptr_new(/allocate)
  self.MOD06ReflectivityB02AC           = ptr_new(/allocate)
  self.MOD06ReflectivityB07AC           = ptr_new(/allocate)
  self.MOD06CloudMask1km                = ptr_new(/allocate)
  self.MOD06QualityAssurance1km         = ptr_new(/allocate)
  return, 1
end
pro EosModis1km::cleanup
  ptr_free, self.MOD03SpacecraftPosition
  ptr_free, self.MOD02ReflectivityProductB01
  ptr_free, self.MOD02ReflectivityProductB02
  ptr_free, self.MOD02ReflectivityProductB03
  ptr_free, self.MOD02ReflectivityProductB04
  ptr_free, self.MOD02ReflectivityProductB07
  ptr_free, self.MOD02RadianceB29
  ptr_free, self.MOD02RadianceB31
  ptr_free, self.MOD02RadianceB32
  ptr_free, self.MOD06CloudPhaseInfrared
  ptr_free, self.MOD06CloudTopPressureInfrared
  ptr_free, self.MOD06CloudTopTemperatureInfrared
  ptr_free, self.MOD06CloudTopHeightInfrared
  ptr_free, self.MOD06CloudPhaseOpticalProperties
  ptr_free, self.MOD06CloudOpticalThickness
  ptr_free, self.MOD06CloudEffectiveRadius
  ptr_free, self.MOD06CloudEmissivity11
  ptr_free, self.MOD06CloudEmissivity12
  ptr_free, self.MOD06CloudEmissivity13
  ptr_free, self.MOD06CloudMaskSPIB01
  ptr_free, self.MOD06CloudMaskSPIB02
  ptr_free, self.MOD06ReflectivityB02AC
  ptr_free, self.MOD06ReflectivityB07AC
  ptr_free, self.MOD06CloudMask1km
  ptr_free, self.MOD06QualityAssurance1km
  obj_destroy, self.position
  self->IntfData::cleanup
  return
end
;-----------------------------------
; override of IntfData class methods
function EosModis1km::createSelector
  obj = obj_new('EosSelector', self)
  return, obj
end
pro EosModis1km::collocatePosition, target
  obj = obj_new('EosLink', self)
  obj->collocatePosition, target
end
pro EosModis1km::collocateGrid, target
  message, 'Not yet implemented'
;  obj = obj_new('EosLink', self)
;  obj->collocateGrid, target
end
;-----------------------------------
; methods
;+
; loadPosGeom
;    Procedure : load position and geometry information from
;                MOD03/MYD03 data file or interpolate from the
;                given 5km data in a given MOD/MYD file
;    Args:
;      MODfilename  (character): file name of the product that
;                                triggered this loading function
;      sizx              (long): data size x
;      sizy              (long): data size y
;-
pro EosModis1km::loadPosGeom, MODfilename, sizx, sizy
  if self.filename ne '' then begin
    fname = self.filename
    intp  = 'NO'
  endif else begin
    if n_elements(MODfilename) eq 0 then message, 'MOD03 file and MODIS product file name are not given.'
    fname = MODfilename
    intp  = 'YES'
  endelse

  fid = hdf_sd_Start(fname, /read)
  did_lat = hdf_sd_NameToIndex(fid,'Latitude')
  did_lon = hdf_sd_NameToIndex(fid,'Longitude')
  if strmid(file_basename(fname),3,2) eq '06' then begin
    did_vza = hdf_sd_NameToIndex(fid,'Sensor_Zenith')
    did_vaa = hdf_sd_NameToIndex(fid,'Sensor_Azimuth')
    did_sza = hdf_sd_NameToIndex(fid,'Solar_Zenith')
    did_saa = hdf_sd_NameToIndex(fid,'Solar_Azimuth')
  endif else begin
    did_vza = hdf_sd_NameToIndex(fid,'SensorZenith')
    did_vaa = hdf_sd_NameToIndex(fid,'SensorAzimuth')
    did_sza = hdf_sd_NameToIndex(fid,'SolarZenith')
    did_saa = hdf_sd_NameToIndex(fid,'SolarAzimuth')
  endelse
  dat_lat = hdf_sd_Select(fid, did_lat)
  dat_lon = hdf_sd_Select(fid, did_lon)
  dat_vza = hdf_sd_Select(fid, did_vza)
  dat_vaa = hdf_sd_Select(fid, did_vaa)
  dat_sza = hdf_sd_Select(fid, did_sza)
  dat_saa = hdf_sd_Select(fid, did_saa)
  aid_vza = hdf_sd_AttrFind(dat_vza, 'scale_factor')
  aid_vaa = hdf_sd_AttrFind(dat_vaa, 'scale_factor')
  aid_sza = hdf_sd_AttrFind(dat_sza, 'scale_factor')
  aid_saa = hdf_sd_AttrFind(dat_saa, 'scale_factor')
  hdf_sd_GetData, dat_lat, lat
  hdf_sd_GetData, dat_lon, lon
  hdf_sd_GetData, dat_vza, vza
  hdf_sd_GetData, dat_vaa, vaa ; clockwize from the north
  hdf_sd_GetData, dat_sza, sza
  hdf_sd_GetData, dat_saa, saa ; clockwize from the north
  hdf_sd_AttrInfo, dat_vza, aid_vza, data=scl_vza
  hdf_sd_AttrInfo, dat_vaa, aid_vaa, data=scl_vaa
  hdf_sd_AttrInfo, dat_sza, aid_sza, data=scl_sza
  hdf_sd_AttrInfo, dat_saa, aid_saa, data=scl_saa
  hdf_sd_end, fid

  vza = vza * scl_vza[0]
  sza = sza * scl_sza[0]
  vaa = vaa * scl_vaa[0]
  saa = saa * scl_saa[0]

  if (intp eq 'YES') then begin
    sizd = size(lat)
    if strmid(file_basename(fname),3,2) eq '06' then begin
      ; MOD06 in Collection 6 fails to carry one element
      lat = [lat, 2.0*lat[(sizd(1)-1),*]-lat[(sizd(1)-2),*]]
      lon = [lon, 2.0*lon[(sizd(1)-1),*]-lon[(sizd(1)-2),*]]
      vza = [vza, 2.0*vza[(sizd(1)-1),*]-vza[(sizd(1)-2),*]]
      sza = [sza, 2.0*sza[(sizd(1)-1),*]-sza[(sizd(1)-2),*]]
      vaa = [vaa, 2.0*vaa[(sizd(1)-1),*]-vaa[(sizd(1)-2),*]]
      saa = [saa, 2.0*saa[(sizd(1)-1),*]-saa[(sizd(1)-2),*]]
    endif
    sizd = size(lat)
    ; check if near date line
    mn = min(lon, max=mx)
    if (min(abs(lon)) ge 90) and (mx - mn ge 350) then ndl = 1 else ndl = 0
    lon = lon + 360.0 * (fix(-(lon-360)/360) > 0) * ndl ; [-360:0] --> [0:360]
        lat = congridext(lat,5*sizd(1),5*sizd(2)) & lat = lat[0:sizx-1,0:sizy-1]
    lon = congridext(lon,5*sizd(1),5*sizd(2)) & lon = lon[0:sizx-1,0:sizy-1]
    lon = lon - 360.0 * float( fix(lon/180.0)/2 + fix(lon/180.0) mod 2 ) ; remap to [-180:180]

    vac = congridext(tan(!DPI/180.0*vza)*cos(!DPI/180.0*vaa), sizx, sizy)
    vas = congridext(tan(!DPI/180.0*vza)*sin(!DPI/180.0*vaa), sizx, sizy)
    vza = 180.0 / !DPI * atan(sqrt(vas^2+vac^2), 1.0) & vza = vza[0:sizx-1,0:sizy-1]
    vaa = 180.0 / !DPI * atan(vas,vac)                & vaa = vaa[0:sizx-1,0:sizy-1]

    sza = congridext(sza,5*sizd(1),5*sizd(2)) & sza = sza[0:sizx-1,0:sizy-1]
    saa = congridext(saa,5*sizd(1),5*sizd(2)) & saa = saa[0:sizx-1,0:sizy-1]
  endif

  raa = -(vaa-saa)             ; azimuth angles are measured clockwise from the local north
  raa = temporary(raa) + 180.0 ; solar azimuth angle is reported as the direction to the sun
  raa = raa - 360.0 * float( fix(raa/180.0)/2 + fix(raa/180.0) mod 2 ) ; remap to [-180:180]

  self.geometry->setGeometry, sza, 270.0-saa, vza, raa
  self.position->setPosition, lat, lon
  self.recordDimension = [sizx,sizy]

  if lat[0] le lat[n_elements(lat)-1UL] then begin
    self.plotReverse = 0
  endif else begin
    self.plotReverse = 1
  endelse

end
;+
; loadMOD06
;    Procedure : load MODIS Level 2 cloud product
;    Args:
;      filenameMOD06 (string): name of MOD06 file
;-
pro EosModis1km::loadMOD06, filenameMOD06

  ; dataset pointer
  nDset     = 14
  sds       = strarr(nDset)
  dset      = ptrarr(nDset)
  proc      = bytarr(nDset) ; post-processing mode
  brf       = bytarr(nDset)

  ; pointers for special datasets
  ptr_Atm_Corr_Reff  = ptr_new(/allocate)
  ptr_Cloud_Mask_SPI = ptr_new(/allocate)

  ; dataset registrations
  ; Note: to stop loading a dataset, replace the sds name with blank string (e.g., sds[0]='')
  ; Post-processing modes
  ;   0 : No post processing
  ;   1 : Convert to float and add offset then scale, fill value with NaN
  dset[0]  = self.MOD06CloudPhaseInfrared          & proc[0]  = 0 & sds[0]  = 'Cloud_Phase_Infrared_1km'
  dset[1]  = self.MOD06CloudTopPressureInfrared    & proc[1]  = 1 & sds[1]  = 'cloud_top_pressure_1km'
  dset[2]  = self.MOD06CloudTopHeightInfrared      & proc[2]  = 1 & sds[2]  = 'cloud_top_height_1km'
  dset[3]  = self.MOD06CloudTopTemperatureInfrared & proc[3]  = 1 & sds[3]  = 'cloud_top_temperature_1km'
  dset[4]  = self.MOD06CloudPhaseOpticalProperties & proc[4]  = 0 & sds[4]  = 'Cloud_Phase_Optical_Properties'
  dset[5]  = self.MOD06CloudOpticalThickness       & proc[5]  = 1 & sds[5]  = 'Cloud_Optical_Thickness'
  dset[6]  = self.MOD06CloudEffectiveRadius        & proc[6]  = 1 & sds[6]  = 'Cloud_Effective_Radius'
  dset[7]  = self.MOD06CloudEmissivity11           & proc[7]  = 1 & sds[7]  = 'cloud_emiss11_1km'
  dset[8]  = self.MOD06CloudEmissivity12           & proc[8]  = 1 & sds[8]  = 'cloud_emiss12_1km'
  dset[9]  = self.MOD06CloudEmissivity13           & proc[9]  = 1 & sds[9]  = 'cloud_emiss13_1km'
  dset[10] = ptr_Cloud_Mask_SPI                    & proc[10] = 1 & sds[10] = 'Cloud_Mask_SPI'
  dset[11] = ptr_Atm_Corr_Reff                     & proc[11] = 1 & sds[11] = 'Atm_Corr_Refl'
  dset[12] = self.MOD06CloudMask1km                & proc[12] = 0 & sds[12] = 'Cloud_Mask_1km'
  dset[13] = self.MOD06QualityAssurance1km         & proc[13] = 0 & sds[13] = 'Quality_Assurance_1km'

  fid = hdf_sd_Start(filenameMOD06, /read)
  for i=0, nDset-1 do begin
    if sds[i] ne '' then begin
      did = hdf_sd_NameToIndex(fid, sds[i])
      dat = hdf_sd_Select(fid, did)
      aid_scl = hdf_sd_AttrFind(dat, 'scale_factor')
      aid_ofs = hdf_sd_AttrFind(dat, 'add_offset')
      aid_fil = hdf_sd_AttrFind(dat, '_FillValue')
      hdf_sd_AttrInfo, dat, aid_scl, data=scl
      hdf_sd_AttrInfo, dat, aid_ofs, data=ofs
      hdf_sd_AttrInfo, dat, aid_fil, data=fil
      delvar, data
      hdf_sd_GetData, dat, data
      ; post processing
      case proc[i] of
        1 : begin
          invalid = where(data eq fil[0], cnt)
          data = (float(data) - ofs[0]) * scl[0]
          if (cnt ne 0) then data[invalid] = !Values.f_NaN
        end
        else :
      endcase
      *(dset[i]) = data
      siz = size(data)
      sizx = siz[siz[0]-1]
      sizy = siz[siz[0]  ]
    endif
  endfor
  hdf_sd_end, fid

  ; load geometric parameters and position data
  if total(self.recordDimension) eq 0 then self.loadPosGeom, filenameMOD06, sizx, sizy

  ; split dataset
  tgt = where(sds eq 'Atm_Corr_Refl', cnt)
  if cnt eq 1 then begin
    *self.MOD06ReflectivityB02AC = reform((*dset[tgt[0]])[1,*,*])
    *self.MOD06ReflectivityB07AC = reform((*dset[tgt[0]])[4,*,*])
  endif

  ; split dataset
  tgt = where(sds eq 'Cloud_Mask_SPI', cnt)
  if cnt eq 1 then begin
    *self.MOD06CloudMaskSPIB01 = reform((*dset[tgt[0]])[0,*,*])
    *self.MOD06CloudMaskSPIB02 = reform((*dset[tgt[0]])[1,*,*])
  endif

end
;+
; loadMOD021km
;    Procedure : load MODIS Level 1B radiance product
;    Args:
;      filenameMOD021km (string): name of MOD021km file
;-
pro EosModis1km::loadMOD021km, filenameMOD021km

  fid = hdf_sd_Start(filenameMOD021km, /read)
  did_emssv  = hdf_sd_NameToIndex(fid,'EV_1KM_Emissive')
  did_rp250  = hdf_sd_NameToIndex(fid,'EV_250_Aggr1km_RefSB')
  did_rp500  = hdf_sd_NameToIndex(fid,'EV_500_Aggr1km_RefSB')
  dat_emssv  = hdf_sd_Select(fid, did_emssv)
  dat_rp250  = hdf_sd_Select(fid, did_rp250)
  dat_rp500  = hdf_sd_Select(fid, did_rp500)
  aid_emssvs = hdf_sd_AttrFind(dat_emssv, 'radiance_scales')
  aid_rp250s = hdf_sd_AttrFind(dat_rp250, 'reflectance_scales')
  aid_rp500s = hdf_sd_AttrFind(dat_rp500, 'reflectance_scales')
  aid_emssvo = hdf_sd_AttrFind(dat_emssv, 'radiance_offsets')
  aid_rp250o = hdf_sd_AttrFind(dat_rp250, 'reflectance_offsets')
  aid_rp500o = hdf_sd_AttrFind(dat_rp500, 'reflectance_offsets')
  aid_emssvf = hdf_sd_AttrFind(dat_emssv, '_FillValue')
  aid_rp250f = hdf_sd_AttrFind(dat_rp250, '_FillValue')
  aid_rp500f = hdf_sd_AttrFind(dat_rp500, '_FillValue')
  hdf_sd_GetData, dat_emssv, emssv
  hdf_sd_GetData, dat_rp250, rp250
  hdf_sd_GetData, dat_rp500, rp500
  hdf_sd_AttrInfo, dat_emssv, aid_emssvs, data=scl_emssv
  hdf_sd_AttrInfo, dat_rp250, aid_rp250s, data=scl_rp250
  hdf_sd_AttrInfo, dat_rp500, aid_rp500s, data=scl_rp500
  hdf_sd_AttrInfo, dat_emssv, aid_emssvo, data=ofst_emssv
  hdf_sd_AttrInfo, dat_rp250, aid_rp250o, data=ofst_rp250
  hdf_sd_AttrInfo, dat_rp500, aid_rp500o, data=ofst_rp500
  hdf_sd_AttrInfo, dat_emssv, aid_emssvf, data=fill_emssv
  hdf_sd_AttrInfo, dat_rp250, aid_rp250f, data=fill_rp250
  hdf_sd_AttrInfo, dat_rp500, aid_rp500f, data=fill_rp500
  hdf_sd_end, fid

  siz = size(emssv)
  if total(self.recordDimension) eq 0 then self.loadPosGeom, filenameMOD021km, siz[1], siz[2]

  *self.MOD02RadianceB29            = (float(emssv[*,*, 8]) - ofst_emssv[ 8]) * scl_emssv[ 8]
  *self.MOD02RadianceB31            = (float(emssv[*,*,10]) - ofst_emssv[10]) * scl_emssv[10]
  *self.MOD02RadianceB32            = (float(emssv[*,*,11]) - ofst_emssv[11]) * scl_emssv[11]
  *self.MOD02ReflectivityProductB01 = (float(rp250[*,*, 0]) - ofst_rp250[ 0]) * scl_rp250[ 0]
  *self.MOD02ReflectivityProductB02 = (float(rp250[*,*, 1]) - ofst_rp250[ 1]) * scl_rp250[ 1]
  *self.MOD02ReflectivityProductB03 = (float(rp500[*,*, 0]) - ofst_rp500[ 0]) * scl_rp500[ 0]
  *self.MOD02ReflectivityProductB04 = (float(rp500[*,*, 1]) - ofst_rp500[ 1]) * scl_rp500[ 1]
  *self.MOD02ReflectivityProductB07 = (float(rp500[*,*, 4]) - ofst_rp500[ 4]) * scl_rp500[ 4]

  invalid_emssv = where((max(emssv,dimension=3) eq fill_emssv[0]) or (max(emssv,dimension=3) ge 32768), cnt_emssv)
  invalid_rp250 = where((max(rp250,dimension=3) eq fill_rp250[0]) or (max(rp250,dimension=3) ge 32768), cnt_rp250)
  invalid_rp500 = where((max(rp500,dimension=3) eq fill_rp500[0]) or (max(rp500,dimension=3) ge 32768), cnt_rp500)

  if cnt_emssv ne 0 then (*self.MOD02RadianceB29           )[invalid_emssv] = !Values.f_NaN
  if cnt_emssv ne 0 then (*self.MOD02RadianceB31           )[invalid_emssv] = !Values.f_NaN
  if cnt_emssv ne 0 then (*self.MOD02RadianceB32           )[invalid_emssv] = !Values.f_NaN
  if cnt_rp250 ne 0 then (*self.MOD02ReflectivityProductB01)[invalid_rp250] = !Values.f_NaN
  if cnt_rp250 ne 0 then (*self.MOD02ReflectivityProductB02)[invalid_rp250] = !Values.f_NaN
  if cnt_rp500 ne 0 then (*self.MOD02ReflectivityProductB03)[invalid_rp500] = !Values.f_NaN
  if cnt_rp500 ne 0 then (*self.MOD02ReflectivityProductB04)[invalid_rp500] = !Values.f_NaN
  if cnt_rp500 ne 0 then (*self.MOD02ReflectivityProductB07)[invalid_rp500] = !Values.f_NaN

end
;+
; getRecordDimension
;    Function  : return the record dimension array
;    Args:
;      none
;-
function EosModis1km::getRecordDimension
  return, self.recordDimension
end
;+
; getLatLonRange
;    Function  : return the latitude and longitude range
;                returned four element floating point array contains:
;                [latmin, lonmin, latmax, lonmax]
;    Args:
;      none
;-
function EosModis1km::getLatLonRange
  return, self.position->getLatLonRange()
end
;+
; getAngles
;    Function : get solar    zenith  angle (element 0),
;                   viewing  zenith  angle (element 1), and
;                   relative azimuth angle (element 2)
;    Args:
;      sel (intfSelector): selector with pixel of interest flagged
;-
function EosModis1km::getAngles, sel
  ret = self.geometry->getMeanAngles(sel)
  return, ret
end
;+
; getSolarAzimuthAngle
;    Function : get solar azimuth angles
;    Args:
;      sel (intfSelector): selector with pixel of interest flagged
;-
function EosModis1km::getSolarAzimuthAngle, sel
  ret = self.geometry->getMeanSolarAzimuthAngle(sel)
  return, ret
end
;+
; getScatteringAngle
;    Function : get scattering angle
;    Args:
;      sel (intfSelector): selector with pixel of interest flagged (optional)
;-
function EosModis1km::getScatteringAngle, sel
  ret = self.geometry->getMeanScatteringAngle(sel)
  return, ret
end
;+
; getGlitterAngle
;    Function : get glitter angle
;    Args:
;      sel (intfSelector): selector with pixel of interest flagged (optional)
;-
function EosModis1km::getGlitterAngle, sel
  ret = self.geometry->getMeanGlitterAngle(sel)
  return, ret
end
;+
; getLatLon
;    Function  : return the position of specified pixel
;                for invalid pixel id, return [lat,lon]=[0,0]
;                if no pixel id is given, return null pointer
;    Args:
;      sel (intfSelector): selector storing the point of interest
;-
function EosModis1km::getLatLon, sel
  return, self.position->getMeanLatLon(sel)
end
;+
; getSpaceCraftAltitude
;    Function  : return the space craft altitude
;                if MOD03 is loaded, compute the height using the orbital data (not implemented yet)
;                if MOD03 is not loaded, return the constatn value
;    Args:
;      sel (intfSelector): selector storing the point of interest
;-
function EosModis1km::getSpacecraftAltitude, sel
  sela = sel->reformForData(self)
  cntC = sela->countChunk()
  if cntC ne 0 then pixid = sela->getSelected1dIndex(chunk=ulindgen(cntC)+1UL)
  valid = where(pixid ne !NULL, cntV) ; valid chunks
  if cntV eq 0 then begin
    ret = !NULL
  endif else begin
    ret    = fltarr(cntC) ; array to retrun
    if n_elements(*self.MOD03SpacecraftPosition) ne 0 then begin
      ; take the norm of MOD04SpacecraftPosition (L) and subtract the ellipsoid altitude
      ; L - a^2 * (cos(lat)^2 + (1-f)^2*sin(lat)^2)
      ; where a = ellipsoid semi major axis, f = ellipsoid flattening
      ret[*] = 705000.0
    endif else begin
      ret[*] = 705000.0
    endelse
  endelse
  return, ret
end
;+
; getMOD06CloudPhaseInfrared
;    Function  : return the IR Cloud Phase
;                  0 -- cloud free
;                  1 -- water cloud
;                  2 -- ice cloud
;                If multiple pixels are provided in a chunk, and contradicts
;                to each other, 3 -- undetermined is assigned
;
;    Args:
;      sel (intfSelector): selector storing the point of interest
;-
function EosModis1km::getMOD06CloudPhaseInfrared, sel
  return, self->extractData(sel, self.MOD06CloudPhaseInfrared, mode='flag', fill=3)
end
;+
; getMOD06CloudPhaseOpticalProperties
;    Function  : return the Cloud Phase
;                  0 -- cloud mask undetermined
;                  1 -- clear sky
;                  2 -- liquid water cloud
;                  3 -- ice cloud
;                  4 -- undetermined phase cloud
;                If multiple pixels are provided in a chunk, and contradicts
;                to each other, 4 -- undetermined is assigned
;
;    Args:
;      sel (intfSelector): selector storing the point of interest
;-
function EosModis1km::getMOD06CloudPhaseOpticalProperties, sel
  return, self->extractData(sel, self.MOD06CloudPhaseOpticalProperties, mode='flag', fill=4)
end

;function EosModis1km::getMOD02RadianceB29, sel
;  return, self->extractData(sel, *self.MOD02RadianceB29, mode='median', fill=-999)
;end


;+
; getMOD06CloudTopPressureInfrared
;    Function  : return the Cloud Top Pressure
;                If multiple pixels are provided in a chunk, the median is reported
;
;    Args:
;      sel (intfSelector): selector storing the point of interest
;-
function EosModis1km::getMOD06CloudTopPressureInfrared, sel
  return, self->extractData(sel, self.MOD06CloudTopPressureInfrared, mode='median', fill=-999)
end
;+
; getMOD06CloudTopHeightInfrared
;    Function  : return the Cloud Top Height
;                If multiple pixels are provided in a chunk, the median is reported
;
;    Args:
;      sel (intfSelector): selector storing the point of interest
;-
function EosModis1km::getMOD06CloudTopHeightInfrared, sel
  return, self->extractData(sel, self.MOD06CloudTopHeightInfrared, mode='median', fill=-999)
end
;+
; getMOD06CloudTopTemperatureInfrared
;    Function  : return the Cloud Top Temperature
;                If multiple pixels are provided in a chunk, the median is reported
;
;    Args:
;      sel (intfSelector): selector storing the point of interest
;-
function EosModis1km::getMOD06CloudTopTemperatureInfrared, sel
  return, self->extractData(sel, self.MOD06CloudTopTemperatureInfrared, mode='median', fill=-999)
end
;+
; getMOD06CloudOpticalThickness
;    Function  : return the Cloud Optical Thickness
;                If multiple pixels are provided in a chunk, the mean is reported
;
;    Args:
;      sel (intfSelector): selector storing the point of interest
;-
function EosModis1km::getMOD06CloudOpticalThickness, sel
  return, self->extractData(sel, self.MOD06CloudOpticalThickness, mode='mean', fill=-999)
end
;+
; getMOD06CloudEffectiveRadius
;    Function  : return the Cloud Effective Radius
;                If multiple pixels are provided in a chunk, the mean is reported
;
;    Args:
;      sel (intfSelector): selector storing the point of interest
;-
function EosModis1km::getMOD06CloudEffectiveRadius, sel
  return, self->extractData(sel, self.MOD06CloudEffectiveRadius, mode='mean', fill=-999)
end
;+
; getMOD06SubpixelHeterogeneityIndex
;    Function  : return the cloud subpixel heterogeneity index
;                If multiple pixels are provided in a chunk, the mean is reported
;
;    Function  : return the variance of reflectivity products from the specified band
;    Args:
;      sel (intfSelector): selector storing the point of interest
;      band     (integer): band of the data (1,2)
;-
function EosModis1km::getMOD06SubpixelHeterogeneityIndex, sel, band
  if n_elements(sel)  eq 0 then message, 'Selector not given'
  if n_elements(band) eq 0 then message, 'Band not given'
  case band of
    1 : ptr = self.MOD06CloudMaskSPIB01
    2 : ptr = self.MOD06CloudMaskSPIB02
    else : message, 'Improper band number'
  endcase
  return, self->extractData(sel, ptr, mode='mean', fill=-999)
end
;+
; getMOD06ReflectivityAC
;    Function  : return the reflectivity from the specified band with atmospheric
;                correction applied. Note that you should multiply with cos(sza).
;                If multiple pixels are provided in a chunk, the mean is reported
;
;    Args:
;      sel (intfSelector): selector storing the point of interest
;      band     (integer): band of the data (2,7)
;      filter   (char)   : filter to be applied (None, Sobel, Laplace, Laplace5, L5oG1)
;-
function EosModis1km::getMOD06ReflectivityAC, sel, band, filter=filter
  if n_elements(sel)  eq 0 then message, 'Selector not given'
  if n_elements(band) eq 0 then message, 'Band not given'
  if n_elements(filter)  eq 0 then filter='None'
  case band of
    2 : ptr = self.MOD06ReflectivityB02AC
    7 : ptr = self.MOD06ReflectivityB07AC
    else : message, 'Improper band number'
  endcase
  if filter ne 'None' then dataArr = ptr_new(/allocate)
  case filter  of
    'Sobel'   : *dataArr = sobel(*ptr)
    'SoG2'    : *dataArr = sobel(gauss_smooth(*ptr, 2.0))
    'Laplace' : *dataArr = laplacian(*ptr)
    'Laplace5': *dataArr = laplacian(*ptr, kernel_size=5)
    'L5oG1'   : *dataArr = gauss_smooth(laplacian(*ptr, kernel_size=5), 1.0)
    'L5oG2'   : *dataArr = gauss_smooth(laplacian(*ptr, kernel_size=5), 2.0)
    'L3oG2'   : *dataArr = laplacian(gauss_smooth(*ptr, 2.0))
    'None'    :  dataArr = ptr
  endcase
  return, self->extractData(sel, dataArr, mode='mean', fill=-999)
end
;+
; getMOD06CloudEmissivity
;    Function  : return the Cloud Emissivity
;                If multiple pixels are provided in a chunk, the mean is reported
;
;    Args:
;      sel (intfSelector): selector storing the point of interest
;      band     (integer): band of the data (11,12,13)
;-
function EosModis1km::getMOD06CloudEmissivity, sel, band
  if n_elements(sel)  eq 0 then message, 'Selector not given'
  if n_elements(band) eq 0 then message, 'Band number is not given'
  sela = sel->reformForData(self)
  cntC = sela->countChunk()
  case band of
   11 : ptr = self.MOD06CloudEmissivity11
   12 : ptr = self.MOD06CloudEmissivity12
   13 : ptr = self.MOD06CloudEmissivity13
   else : message, 'Improper band number'
  endcase
  return, self->extractData(sel, ptr, mode='mean', fill=-999)
end
;+
; getMOD06CloudMask1km
;    Function  : return the cloud mask flags
;    Args:
;      sel  (intfSelector): selector storing the point of interest
;      flagname   (string): flag names
;    Description:
;      This function returns value of flags. Please see QA plan document by
;      MODIS science team for the meanings of values.
;      You can extract single flag data by specifying single element character,
;        (example) obj->getMOD06CloudMask1km, sel, 'CloudMaskCloudiness'
;      or multiple camera data by specifying an array of string.
;        (example) obj->getMOD06CloudMask1km, sel, ['CloudMaskCloudiness','SnowIce']
;-
function EosModis1km::getMOD06CloudMask1km, sel, flagnames
  if n_elements(sel)       eq 0 then message, 'Selector not given'
  nFlg = n_elements(flagnames)
  if nFlg eq 0 then message, 'Flag name not given'
  sela = sel->reformForData(self)
  cntC = sela->countChunk()
  if cntC ne 0 then pixid = sela->getSelected1dIndex(chunk=ulindgen(cntC)+1UL)
  valid  = where(pixid ne !NULL, cntV) ; valid chunks
  if cntV eq 0 then begin
    ret = !NULL
  endif else begin
    ret    = bytarr(cntC, nFlg) ; array to retrun
    ret[*] = 255B ; missing value: 
    for iFlg=0, nFlg-1 do begin
      case flagnames[iFlg] of
       'CloudMaskStatus'     : begin & slice = 0 & mask  = '00000001'BB & end
       'CloudMaskCloudiness' : begin & slice = 0 & mask  = '00000110'BB & end
       'DayNight'            : begin & slice = 0 & mask  = '00001000'BB & end
       'Sunglint'            : begin & slice = 0 & mask  = '00010000'BB & end
       'SnowIce'             : begin & slice = 0 & mask  = '00100000'BB & end
       'SurfaceType'         : begin & slice = 0 & mask  = '11000000'BB & end
       'HeavyAerosol'        : begin & slice = 1 & mask  = '00000001'BB & end
       'ThinCirrus'          : begin & slice = 1 & mask  = '00000010'BB & end
       'Shadow'              : begin & slice = 1 & mask  = '00000100'BB & end
       else : message, 'Improper flag name'
      endcase
      dat = (reform((*self.MOD06CloudMask1km)[slice,*,*]) and mask) / (mask and (mask xor byte(2*mask)))
      for i=0UL, cntV-1UL do begin
        val = dat[pixid[valid[i]]-1UL]
        if product(val eq val[0]) eq 1 $
          then ret[valid[i],iFlg] = val[0]  $
          else ret[valid[i],iFlg] = 255 ; do not agree with each other: 255 -- undetermined
      endfor
    endfor
  endelse
  return, ret
end
;+
; getMOD06QualityAssurance1km
;    Function  : return the quality assurance flags
;    Args:
;      sel  (intfSelector): selector storing the point of interest
;      flagname   (string): flag names
;    Description:
;      This function returns value of flags. Please see QA plan document by
;      MODIS science team for the meanings of values.
;      You can extract single flag data by specifying single element character,
;        (example) obj->getMOD06QualityAssurance1km, sel, 'CloudMaskCloudiness'
;      or multiple camera data by specifying an array of string.
;        (example) obj->getMOD06QualityAssurance1km, sel, ['CloudMaskCloudiness','SnowIce']
;-
function EosModis1km::getMOD06QualityAssurance1km, sel, flagnames
  if n_elements(sel)       eq 0 then message, 'Selector not given'
  nFlg = n_elements(flagnames)
  if nFlg eq 0 then message, 'Flag name not given'
  sela = sel->reformForData(self)
  cntC = sela->countChunk()
  if cntC ne 0 then pixid = sela->getSelected1dIndex(chunk=ulindgen(cntC)+1UL)
  valid  = where(pixid ne !NULL, cntV) ; valid chunks
  if cntV eq 0 then begin
    ret = !NULL
  endif else begin
    ret    = bytarr(cntC, nFlg) ; array to retrun
    ret[*] = 255 ; missing value: 255
    for iFlg=0, nFlg-1 do begin
      case flagnames[iFlg] of
       'PrimaryCloudOpticalThicknessUsefulness'     : begin & slice = 0 & mask  = '00000001'BB & end
       'PrimaryCloudOpticalThicknessConfidence'     : begin & slice = 0 & mask  = '00000110'BB & end
       'PrimaryCloudEffectiveRadiusUsefulness'      : begin & slice = 0 & mask  = '00010001'BB & end
       'PrimaryCloudEffectiveRadiusConfidence'      : begin & slice = 0 & mask  = '01100000'BB & end
       'PrimaryCloudWaterPathUsefulness'            : begin & slice = 1 & mask  = '00000001'BB & end
       'PrimaryCloudWaterPathConfidence'            : begin & slice = 1 & mask  = '00000110'BB & end
       '1.6-2.1CloudRetrievalPhase'                 : begin & slice = 1 & mask  = '00111000'BB & end
       '1.6-2.1CloudRetrievalOutcome'               : begin & slice = 1 & mask  = '01000000'BB & end
       '1.6-2.1CloudRetrievalPhase&Outcome'         : begin & slice = 1 & mask  = '01111000'BB & end
       'PrimaryCloudRetrievalPhase'                 : begin & slice = 2 & mask  = '00000111'BB & end
       'PrimaryCloudRetrievalOutcome'               : begin & slice = 2 & mask  = '00001000'BB & end
       'PrimaryCloudRetrievalPhase&Outcome'         : begin & slice = 2 & mask  = '00001111'BB & end
       'RayleighCorrection'                         : begin & slice = 2 & mask  = '00010000'BB & end
       'AtmosphericWaterVaporCorrection'            : begin & slice = 2 & mask  = '00100000'BB & end
       'BandUsedForPriparyOpticalThicknessRetrieval': begin & slice = 2 & mask  = '11000000'BB & end
       '1.6-2.1CloudOpticalThicknessUsefulnessFlag' : begin & slice = 3 & mask  = '00000001'BB & end
       '1.6-2.1CloudOpticalThicknessConfidenceFlag' : begin & slice = 3 & mask  = '00000110'BB & end
       '1.6-2.1CloudEffectiveRadiusUsefulnessFlag'  : begin & slice = 3 & mask  = '00001000'BB & end
       '1.6-2.1CloudEffectiveRadiusConfidenceFlag'  : begin & slice = 3 & mask  = '00110000'BB & end
       'ClearSkyRestoralType'                       : begin & slice = 3 & mask  = '11000000'BB & end
       '1.6-2.1CloudWaterPathUsefulnessFlag'        : begin & slice = 4 & mask  = '00000001'BB & end
       '1.6-2.1CloudWaterPathConfidenceFlag'        : begin & slice = 4 & mask  = '00000110'BB & end
       'PrimaryCloudRetrievalMultilayerCloud&Phase' : begin & slice = 4 & mask  = '00111000'BB & end
       'PhaseDifferenceMultilayerTest'              : begin & slice = 5 & mask  = '00000001'BB & end
       'DeltaPrecipitableWaterMultilayerTest'       : begin & slice = 5 & mask  = '00000010'BB & end
       'DeltaPrecipitableWaterAt900mbTest'          : begin & slice = 5 & mask  = '00000100'BB & end
       'TauDifferenceVIS-NIRMultilayerTest'         : begin & slice = 5 & mask  = '00001000'BB & end
       'Pavolonis-HeidingerMultilayerTest'          : begin & slice = 5 & mask  = '00010000'BB & end
       'VISWIR-1.6CloudRetrievalPhase&Outcome'      : begin & slice = 6 & mask  = '00001111'BB & end
       'VISWIR-1.6PCLCloudRetrievalPhase&Outcome'   : begin & slice = 6 & mask  = '11110000'BB & end
       'VISWIR-3.7CloudRetrievalPhase&Outcome'      : begin & slice = 7 & mask  = '00001111'BB & end
       'VISWIR-3.7PCLCloudRetrievalPhase&Outcome'   : begin & slice = 7 & mask  = '11110000'BB & end
       'VISWIR-2.1CloudRetrievalPhase&Outcome'      : begin & slice = 8 & mask  = '00001111'BB & end
       'VISWIR-2.1PCLCloudRetrievalPhase&Outcome'   : begin & slice = 8 & mask  = '11110000'BB & end
       else : message, 'Improper flag name'
      endcase
      dat = (reform((*self.MOD06QualityAssurance1km)[slice,*,*]) and mask) / (mask and (mask xor byte(2*mask)))
      for i=0UL, cntV-1UL do begin
        val = dat[pixid[valid[i]]-1UL]
        if product(val eq val[0]) eq 1 $
          then ret[valid[i],iFlg] = val[0]  $
          else ret[valid[i],iFlg] = 255 ; do not agree with each other: 255 -- undetermined
      endfor
    endfor
  endelse
  return, ret
end
;+
; getMOD02ReflectivityProductB02
;    Function  : return the reflectivity product from Band 02
;    Args:
;      sel (intfSelector): selector storing the point of interest
;    Note:
;      THIS FUNCTION IS KEPT FOR BACKWARD COMPATIBILITY.
;      PLEASE USE getMOD02ReflectivityProduct() for new programs.
;-
function EosModis1km::getMOD02ReflectivityProductB02, sel
  return, self->getMOD02ReflectivityProduct(sel, 2)
end
;+
; getMOD02ReflectivityProductB07
;    Function  : return the reflectivity product from Band 07
;    Args:
;      sel (intfSelector): selector storing the point of interest
;    Note:
;      THIS FUNCTION IS KEPT FOR BACKWARD COMPATIBILITY.
;      PLEASE USE getMOD02ReflectivityProduct() for new programs.
;-
function EosModis1km::getMOD02ReflectivityProductB07, sel
  return, self->getMOD02ReflectivityProduct(sel, 7)
end
;+
; getMOD02ReflectivityProduct
;    Function  : return the reflectivity product from the specified band
;    Args:
;      sel (intfSelector): selector storing the point of interest
;      band     (integer): band number
;      filter      (char): filter to be applied (None, Sobel, Laplace, Laplace5, L5oG1)
;-
function EosModis1km::getMOD02ReflectivityProduct, sel, band, filter=filter
  if n_elements(sel)  eq 0 then message, 'Selector not given'
  if n_elements(band) eq 0 then message, 'Band not given'
  if n_elements(filter)  eq 0 then filter='None'
  case band of
    1 : ptr = self.MOD02ReflectivityProductB01
    2 : ptr = self.MOD02ReflectivityProductB02
    3 : ptr = self.MOD02ReflectivityProductB03
    4 : ptr = self.MOD02ReflectivityProductB04
    7 : ptr = self.MOD02ReflectivityProductB07
   else : message, 'Improper band number'
  endcase
  help, *ptr
  if filter ne 'None' then dataArr = ptr_new(/allocate)
  case filter  of
    'Sobel'   : *dataArr = sobel(*ptr)
    'SoG2'    : *dataArr = sobel(gauss_smooth(*ptr, 2.0))
    'Laplace' : *dataArr = laplacian(*ptr)
    'Laplace5': *dataArr = laplacian(*ptr, kernel_size=5)
    'L5oG1'   : *dataArr = gauss_smooth(laplacian(*ptr, kernel_size=5), 1.0)
    'L5oG2'   : *dataArr = gauss_smooth(laplacian(*ptr, kernel_size=5), 2.0)
    'L3oG2'   : *dataArr = laplacian(gauss_smooth(*ptr, 2.0))
    'None'    :  dataArr = ptr
  endcase
  return, self->extractData(sel, dataArr, mode='mean', fill=-999)
end
;+
; getMOD02ReflectivityProductVariance
;    Function  : return the variance of reflectivity products from the specified band
;    Args:
;      sel (intfSelector): selector storing the point of interest
;      band      (integer): band number
;-
function EosModis1km::getMOD02ReflectivityProductVariance, sel, band
  if n_elements(sel)  eq 0 then message, 'Selector not given'
  if n_elements(band) eq 0 then message, 'Band not given'
  case band of
    2 : ptr = self.MOD02ReflectivityProductB02
    7 : ptr = self.MOD02ReflectivityProductB07
   else : message, 'Improper band number'
  endcase
  return, self->extractData(sel, ptr, mode='variance', fill=-999)
end
;+
; getMOD02BrightnessTemperature
;    Function  : return the reflectivity product from the specified band
;    Args:
;      sel  (intfSelector): selector storing the point of interest
;      band      (integer): band number
;
;    The center wavelength from Guenther et al. (1996, JAOT, 13, 274-285)
;    --------------------
;      27     6.715
;      28     7.325
;      29     8.550
;      30     9.730
;      31    11.030
;      32    12.020
;      33    13.335
;      34    13.635
;      35    13.935
;      36    14.235
;-
function EosModis1km::getMOD02BrightnessTemperature, sel, band
  if n_elements(sel)  eq 0 then message, 'Selector not given'
  if n_elements(band) eq 0 then message, 'Band not given'
  case band of
   29 : begin
        ptr = self.MOD02RadianceB29
        wvl =  8.550
     end
   31 : begin
        ptr = self.MOD02RadianceB31
        wvl = 11.030
     end
   32 : begin
        ptr = self.MOD02RadianceB32
        wvl = 12.020
     end
   else : message, 'Improper band number'
  endcase
;  return, self->extractData(sel, ptr, mode='meanPlanck', fill=-999, wavelength=wvl)
   return, self->extractData(sel, ptr, mode='median', fill=-999)
end

;return, self->extractData(sel, *self.MOD02RadianceB29, mode='median', fill=-999)

;+
; writeImage
;    Function  : write out data of a channel with image-compatible 2d array
;    Args:
;      band  (integer)        : band to be retireved (MOD02ReflectivityProductB02)
;-
function EosModis1km::writeImage, band=band, gammac=gammac, mincap=mincap, maxcap=maxcap, cbaxis=cbaxis
  if n_elements(band) eq 0 then message, 'Band is not provided'
  if n_elements(mincap)  eq 0 then mincap=0.0
  if n_elements(maxcap)  eq 0 then maxcap=0.6
  if n_elements(gammac)  eq 0 then gammac='None'
  case band of
    1 : dataArr = self.MOD02ReflectivityProductB01
    2 : dataArr = self.MOD02ReflectivityProductB02
    3 : dataArr = self.MOD02ReflectivityProductB03
    4 : dataArr = self.MOD02ReflectivityProductB04
    7 : dataArr = self.MOD02ReflectivityProductB07
    29: dataArr = self.MOD02RadianceB29
    31: dataArr = self.MOD02RadianceB31
    32: dataArr = self.MOD02RadianceB32
  endcase
  if *dataArr eq !NULL then message, 'Radiance data not available in the specified band'
  img = *dataArr
  ; assign 0 to the invalid values
  invalid = where(~ finite(img), cnt)
  if cnt ne 0 then img[invalid] = 0.0

  ; apply caps
  img = (img - mincap) / (maxcap - mincap) > 0.0 < 1.0

  ;gamma correction
  case gammac of
    'sRGB' : img = 12.92*img*(img le 0.0031308) + (1.055*img^0.416666-0.055)*(img gt 0.0031308)
    'None' : img = img
  endcase

  ;convert to byte
  img = bytscl(img)

  ; record the axis bounds for colorbar
  cbaxis = obj_new('AxisScale')
  cbaxis->setAxis, [mincap, maxcap]

  ; reverse image if needed
  if self.plotReverse eq 1 then begin
    img = reverse(img, 2)
  endif else begin
    img = reverse(img, 1)
  endelse

  return, img
end
;+
; writePng
;    Procedure : plot the data
;    Args:
;      band     (string)        : plot mode (visible,pressure,1,2,3,4,7)
;      filename (string)        : filename of output
;      sel      (intfSelector): color the selected point (optional)
;-
pro EosModis1km::writePng, band=band, filename=filename, sel=sel, mincap=mincap, maxcap=maxcap, $
                           grid=grid, colorbar=colorbar, cbarea=cbarea, gammac=gammac
  if n_elements(band)     eq 0 then band='2'
  if n_elements(filename) eq 0 then filename='modis1km.png'
  if (n_elements(gammac)  eq 0) and (band eq 'visible') then gammac='sRGB'
  if n_elements(gammac)   eq 0 then gammac='None'
  ; extract image
  img = obj_new('ImagePlate')
  case band of
    'visible'     : begin
      if (n_elements(*self.MOD02ReflectivityProductB01) eq 0) or $
         (n_elements(*self.MOD02ReflectivityProductB03) eq 0) or $
         (n_elements(*self.MOD02ReflectivityProductB04) eq 0) then message, 'Missing band to write png image'
      img->addToR, self->writeImage(band=1, gammac=gammac, maxcap=maxcap, mincap=mincap)
      img->addToG, self->writeImage(band=4, gammac=gammac, maxcap=maxcap, mincap=mincap)
      img->addToB, self->writeImage(band=3, gammac=gammac, maxcap=maxcap, mincap=mincap)
    end
    'pressure'   : begin
      dataArr = self.MOD06CloudTopPressureInfrared
      if n_elements(*dataArr) eq 0 then message, 'Cloud top pressure data is not given'
      ; reverse image if needed
      if self.plotReverse eq 1 then begin
        dataPlt = reverse(float(*dataArr),2)
      endif else begin
        dataPlt = reverse(float(*dataArr),1)
      endelse
      blt  = self->writeImage(band=2, gammac=gammac, maxcap=maxcap, mincap=mincap)
      ccld = applyColorShading(blt, float(dataPlt)/750*255, mode='EsGnAm')
      cbaxis = obj_new('AxisScale')
      cbaxis->setAxis, [0.0, 750.0]
      cbshading = 'EsGnAm'
      img->addToR, ccld[*,*,0]
      img->addToG, ccld[*,*,1]
      img->addToB, ccld[*,*,2]
    end
    '1'     : begin
      if n_elements(*self.MOD02ReflectivityProductB01) eq 0 then message, 'Missing band to write png image'
      img->addToR, self->writeImage(band=1, gammac=gammac, maxcap=maxcap, mincap=mincap)
      img->addToG, self->writeImage(band=1, gammac=gammac, maxcap=maxcap, mincap=mincap)
      img->addToB, self->writeImage(band=1, gammac=gammac, maxcap=maxcap, mincap=mincap, cbaxis=cbaxis)
    end
    '2' : begin
      if (n_elements(*self.MOD02ReflectivityProductB02)   eq 0) then message, 'Missing band to write png image'
      img->addToR, self->writeImage(band=2, gammac=gammac, maxcap=maxcap, mincap=mincap)
      img->addToG, self->writeImage(band=2, gammac=gammac, maxcap=maxcap, mincap=mincap)
      img->addToB, self->writeImage(band=2, gammac=gammac, maxcap=maxcap, mincap=mincap, cbaxis=cbaxis)
    end
    '3'     : begin
      if n_elements(*self.MOD02ReflectivityProductB03) eq 0 then message, 'Missing band to write png image'
      img->addToR, self->writeImage(band=3, gammac=gammac, maxcap=maxcap, mincap=mincap)
      img->addToG, self->writeImage(band=3, gammac=gammac, maxcap=maxcap, mincap=mincap)
      img->addToB, self->writeImage(band=3, gammac=gammac, maxcap=maxcap, mincap=mincap, cbaxis=cbaxis)
    end
    '4'     : begin
      if n_elements(*self.MOD02ReflectivityProductB04) eq 0 then message, 'Missing band to write png image'
      img->addToR, self->writeImage(band=4, gammac=gammac, maxcap=maxcap, mincap=mincap)
      img->addToG, self->writeImage(band=4, gammac=gammac, maxcap=maxcap, mincap=mincap)
      img->addToB, self->writeImage(band=4, gammac=gammac, maxcap=maxcap, mincap=mincap, cbaxis=cbaxis)
    end
    '7'     : begin
      if n_elements(*self.MOD02ReflectivityProductB07) eq 0 then message, 'Missing band to write png image'
      img->addToR, self->writeImage(band=7, gammac=gammac, maxcap=maxcap, mincap=mincap)
      img->addToG, self->writeImage(band=7, gammac=gammac, maxcap=maxcap, mincap=mincap)
      img->addToB, self->writeImage(band=7, gammac=gammac, maxcap=maxcap, mincap=mincap, cbaxis=cbaxis)
    end
    '29'    : begin
      if n_elements(*self.MOD02RadianceB29) eq 0 then message, 'Missing band to write png image'
      img->addToR, self->writeImage(band=29, gammac=gammac, maxcap=maxcap, mincap=mincap)
      img->addToG, self->writeImage(band=29, gammac=gammac, maxcap=maxcap, mincap=mincap)
      img->addToB, self->writeImage(band=29, gammac=gammac, maxcap=maxcap, mincap=mincap, cbaxis=cbaxis)
    end
    '31'    : begin
      if n_elements(*self.MOD02RadianceB31) eq 0 then message, 'Missing band to write png image'
      img->addToR, self->writeImage(band=31, gammac=gammac, maxcap=maxcap, mincap=mincap)
      img->addToG, self->writeImage(band=31, gammac=gammac, maxcap=maxcap, mincap=mincap)
      img->addToB, self->writeImage(band=31, gammac=gammac, maxcap=maxcap, mincap=mincap, cbaxis=cbaxis)
    end
    '32'    : begin
      if n_elements(*self.MOD02RadianceB32) eq 0 then message, 'Missing band to write png image'
      img->addToR, self->writeImage(band=32, gammac=gammac, maxcap=maxcap, mincap=mincap)
      img->addToG, self->writeImage(band=32, gammac=gammac, maxcap=maxcap, mincap=mincap)
      img->addToB, self->writeImage(band=32, gammac=gammac, maxcap=maxcap, mincap=mincap, cbaxis=cbaxis)
    end
    else      : begin
      if (n_elements(*self.MOD02ReflectivityProductB01) eq 0) or $
         (n_elements(*self.MOD02ReflectivityProductB03) eq 0) or $
         (n_elements(*self.MOD02ReflectivityProductB04) eq 0) then message, 'Missing band to write png image'
      img->addToR, self->writeImage(band=1, gammac=gammac, maxcap=maxcap, mincap=mincap)
      img->addToG, self->writeImage(band=4, gammac=gammac, maxcap=maxcap, mincap=mincap)
      img->addToB, self->writeImage(band=3, gammac=gammac, maxcap=maxcap, mincap=mincap)
    end
  endcase

  ; selection coloring
  if n_elements(sel) ne 0 then begin
    sela   = sel->reformForData(self)
    lincol = sela->getSelectedLinCol()
    if lincol ne !NULL then begin
      if self.plotReverse eq 1 then begin
        img->drawPoint, lincol[*,1], self.recordDimension[1]-lincol[*,0]-1L, color='ffff00'
      endif else begin
        img->drawPoint, self.recordDimension[0]-lincol[*,1]-1L, lincol[*,0], color='ffff00'
      endelse
    endif
  endif

  ; adding Lat/Lon grid to the image
  cex = (norm(self.recordDimension) / 1200.0) > 0.8
  if n_elements(grid) ne 0 then begin
    range = self.position->getLatLonRange()
    idxlatb = fix(range[0] / grid[0]) - (range[0] le 0) ; bottom
    idxlatt = fix(range[2] / grid[0]) + (range[2] ge 0) ; top
    idxlonl = fix(range[1] / grid[1]) - (range[1] le 0) ; left
    idxlonr = fix(range[3] / grid[1]) + (range[3] ge 0) ; right
    if idxlonl ge idxlonr then idxlonr = idxlonr + fix(360.0 / grid[1]); at the date line
    lsLat = grid[0] * (idxlatb + lindgen(idxlatt - idxlatb + 1)) ; latitude list
    lsLon = grid[1] * (idxlonl + lindgen(idxlonr - idxlonl + 1)) ; longitude list
    lsLon = lsLon - 360.0 * float( fix(lsLon/180.0)/2 + fix(lsLon/180.0) mod 2 ) ; rescale to [-180.0:180.0]
    ret = self.position->getLatLonGridPoint(lsLat, lsLon)
    if self.plotReverse eq 1 then begin
      lincol = [[self.recordDimension[1]-long(ret[*,2])-1UL], [long(ret[*,3])]]
    endif else begin
      lincol = [[long(ret[*,2])], [self.recordDimension[0]-long(ret[*,3])-1UL]]
    endelse
    img->plotLatLonGrid, lsLat, lsLon, float(ret[*,0:1]), lincol, charsize=cex, thick=cex
  endif

  ; adding colorbar
  if n_elements(colorbar) ne 0 then begin
    img->addColorBar, cbaxis, shading=cbshading, side=colorbar, charsize=cex, thick=cex, grapharea=cbarea
  endif


  ; plot
  img->writePng, filename
end

;-----------------------------------
; definition
pro EosModis1km__define
  void={EosModis1km,                                $
     inherits IntfData                            , $
     filename                         : ''        , $
     plotReverse                      : 0         , $
     recordDimension                  : ulonarr(2), $
     geometry                         : obj_new() , $
     MOD03SpacecraftPosition          : ptr_new() , $
     MOD02ReflectivityProductB01      : ptr_new() , $
     MOD02ReflectivityProductB02      : ptr_new() , $
     MOD02ReflectivityProductB03      : ptr_new() , $
     MOD02ReflectivityProductB04      : ptr_new() , $
     MOD02ReflectivityProductB07      : ptr_new() , $
     MOD02RadianceB29                 : ptr_new() , $
     MOD02RadianceB31                 : ptr_new() , $
     MOD02RadianceB32                 : ptr_new() , $
     MOD06CloudPhaseInfrared          : ptr_new() , $
     MOD06CloudTopPressureInfrared    : ptr_new() , $
     MOD06CloudTopHeightInfrared      : ptr_new() , $
     MOD06CloudTopTemperatureInfrared : ptr_new() , $
     MOD06CloudPhaseOpticalProperties : ptr_new() , $
     MOD06CloudOpticalThickness       : ptr_new() , $
     MOD06CloudEffectiveRadius        : ptr_new() , $
     MOD06CloudEmissivity11           : ptr_new() , $
     MOD06CloudEmissivity12           : ptr_new() , $
     MOD06CloudEmissivity13           : ptr_new() , $
     MOD06CloudMaskSPIB01             : ptr_new() , $
     MOD06CloudMaskSPIB02             : ptr_new() , $
     MOD06ReflectivityB02AC           : ptr_new() , $
     MOD06ReflectivityB07AC           : ptr_new() , $
     MOD06CloudMask1km                : ptr_new() , $
     MOD06QualityAssurance1km         : ptr_new()   $
  }
end
