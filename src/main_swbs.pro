pro main_swbs, mod02=mod02fname, mod06=mod06fname, mod03=mod03fname
  !PATH = !PATH+":~/Documents/Platyco/"
  ;----
  ; 1. Initial setup
  ;--
  ; 1.1 Filenames of Satellite Data
  if n_elements(mod06fname)   eq 0 then $
    mod06fname   = 'MYD06_L2.A2012245.0915.006.2014095023356.hdf'
  if n_elements(mod02fname)   eq 0 then $
    mod02fname   = 'MYD021KM.A2012245.0915.006.2012245193439.hdf'
  if n_elements(mod03fname)   eq 0 then $
    mod03fname   = strjoin(strsplit(string(format='(2A)', $
                   strmid(mod02fname, 0, strlen(mod02fname)-22), '*.hdf'),$
                   'D021[Kk][Mm]',/regex,/extract),'D03')
  print, mod06fname
  print, mod02fname
  print, mod03fname

  ;--
  ; 1.2 Filenames of prt groups
  prtGrpReffname = ['Bispb02r050/prtgrp/Bispb02r050.prtgrp', 'Bispb02r050/prtgrp/Bispb02r050.prtgrp']
  prtGrpAbsfname = ['Bispb07r050/prtgrp/Bispb07r050.prtgrp', 'Bispb07r050/prtgrp/Bispb07r050.prtgrp']
  otspalfname    = ['Bispb02r050/prtgrp/Bispb02r050.otspal', $
                    'Bispb02r050/prtgrp/Bispb02r050.otspal']
  nMdl = n_elements(prtGrpReffname)

  ;--
  ; 1.3 Constants
  outpath = 'out/' ; output path
  err = [6.4e-5, 1.0e-6] ; variance of channel (Ref, Abs)

  ;----
  ; 2. File name completion
  ;--
  ; 2.1 Set paths
  dataPath    = ['', '/cir0/pub/Spaceborne/MYD/2009253/']
  prtGrpPath  = ['', '/cir2/ywang/ScaData/BiSpect/', '/cir2/ywang/ScaData/BiSpect/']
  nDataPath   = n_elements(dataPath)
  nPrtGrpPath = n_elements(prtGrpPath)
  ;--
  ; 2.2 Complete filenames
  for i=nDataPath-1, 0, -1 do begin
    f = (file_search(dataPath[i]+mod06fname  ))[0] & if f ne '' then mod06Cname   = f
    f = (file_search(dataPath[i]+mod02fname  ))[0] & if f ne '' then mod02Cname   = f
    f = (file_search(dataPath[i]+mod03fname  ))[0] & if f ne '' then mod03Cname   = f
  endfor
  prtgrpRefCname = strarr(nMdl)
  prtgrpAbsCname = strarr(nMdl)
  otspalCname    = strarr(nMdl)
  for j=0, nMdl-1 do begin
    for i=nPrtGrpPath-1, 0, -1 do begin
      f = (file_search(prtGrpPath[i]+prtGrpReffname[j]))[0] & if f ne '' then prtGrpRefCname[j] = f
      f = (file_search(prtGrpPath[i]+prtGrpAbsfname[j]))[0] & if f ne '' then prtGrpAbsCname[j] = f
      f = (file_search(prtGrpPath[i]+otspalfname[j]))[0]    & if f ne '' then otspalCname[j] = f
print, 'prtGrpPath[i]+otspalfname[j] = ', prtGrpPath[i]+otspalfname[j]
    endfor
  endfor
  if n_elements(mod06Cname    ) eq 0 then message, 'MOD06_L2 file not found: '+mod06fname
  if n_elements(mod02Cname    ) eq 0 then message, 'MOD021km file not found: '+mod02fname
  if n_elements(mod03Cname    ) eq 0 then mod03Cname = !NULL
  for j=0, nMdl-1 do begin
    if n_elements(prtGrpRefCname[j]) eq 0 then message, 'prtgrp file not found: '+prtGrpReffname[j]
    if n_elements(prtGrpAbsCname[j]) eq 0 then message, 'prtgrp file not found: '+prtGrpAbsfname[j]
    if n_elements(otspalCname[j]   ) eq 0 then message, 'otspal file not found: '+otspalfname[j]
  endfor

print, mod06Cname
print, mod02Cname
print, mod03Cname
print, prtGrpRefCname[0]
print, prtGrpAbsCname[0]
print, otspalCname[0]
;stop

  ;--
  ; 2.3 Define data identifier
  modisid   = strmid(mod06Cname, 34,13,/reverse_offset)
  matchupid = modisid
  print, systime()+' ['+matchupid+'] started '
  ;--
  ; 2.4 Write empty output files
  openw, ui,  outpath+matchupid+'_obs.inv',  /get_lun
  close, ui

  ;----
  ; 3. Load data and register collocation
  ;--
  ; 3.1 MODIS
  mod1k = obj_new('EosModis1km', mod03Cname)
  mod1k->loadMOD021km, mod02Cname
  mod1k->loadMOD06,    mod06Cname
  ; 3.2 Load RT table for retrieval model
  mdl = objarr(nMdl)
  for j=0, nMdl-1 do begin
    mdl[j] = obj_new('SwbsRetModelDLM')
    mdl[j]->setPrtGrpRfl, prtGrpRefCname[j]
    mdl[j]->setPrtGrpAbs, prtGrpAbsCname[j]
    print, 'prtGrpRefCname[j] = ', prtGrpRefCname[j]
    print, 'prtGrpAbsCname[j] = ', prtGrpAbsCname[j]
  endfor
  ; 3.3 Load spherical albedo --> optical thickness conversion table
  spalot  = objarr(nMdl)
  maxspal = fltarr(nMdl)
  for j=0, nMdl-1 do begin
    d = read_ascii(otspalCname[j], comment_symbol='#')
    x = reform(d.Field1[0,*])
    y = reform(d.Field1[1,*])
    spalot[j]  = obj_new('Spline1D', alog(x+1.0), y)
    maxspal[j] = max(y)
  endfor

  ;----
  ; 4 Select pixels of interest
  print, systime()+' ['+matchupid+'] loaded all data'
  sp = mod1k->createSelector()
  sp->selectAll
  ;--
;rad29  = mod1k->MOD02RadianceB29(sp) ; !!!!!!!!!!!!!!!!!!!!!!!!1


  ; 4.1 Filter by cloud phase
  idxice = where(mod1k->getMOD06CloudPhaseOpticalProperties(sp) eq 3, cnt)
  if cnt ne 0 then sp->sieveChunk, idxice+1UL $
              else message, systime()+' ['+matchupid+'] aborted, no valid pixel (no ice clouds)'
  ;--
  ; 4.2 Filter by validity of MODIS retrievals
  idxval = where(mod1k->getMOD06CloudOpticalThickness(sp)    ge 0, cnt)
  if cnt ne 0 then sp->sieveChunk, idxval+1UL $
              else message, systime()+' ['+matchupid+'] aborted, no valid pixel (no successful MOD06 tau retrievals)'
  ; ot < 5 for SplitWindow retrievals, split window can only identify ot < 5
  idxval = where(mod1k->getMOD06CloudOpticalThickness(sp)    le 5, cnt)
  if cnt ne 0 then sp->sieveChunk, idxval+1UL $
              else message, systime()+' ['+matchupid+'] aborted, no valid pixel (no successful MOD06 tau retrievals)'
  idxval = where(mod1k->getMOD06CloudEffectiveRadius(sp)     ge 0, cnt)
  if cnt ne 0 then sp->sieveChunk, idxval+1UL $
              else message, systime()+' ['+matchupid+'] aborted, no valid pixel (no successful MOD06 reff retrievals)'
  idxval = where(mod1k->getMOD06CloudTopPressureInfrared(sp) ge 0, cnt)
  if cnt ne 0 then sp->sieveChunk, idxval+1UL $
              else message, systime()+' ['+matchupid+'] aborted, no valid pixel (no successful MOD06 ctp retrievals)'
  ;--
  ; 4.2 Filter by data availability
  idxval = where(mod1k->getMOD02BrightnessTemperature(sp,31) ge 0, cnt)
  if cnt ne 0 then sp->sieveChunk, idxval+1UL $
              else message, systime()+' ['+matchupid+'] aborted, no valid pixel (no IR brightness temperatures)'

  ;----
  ; 5 Extract non-directional data for recording and initial conidtion
  print, systime()+' ['+matchupid+'] retrieval candidates  '+string(format='(I8)',sp->countChunk())+' pixels'
  pixid  = sp->getSelected1dIndex(chunk=ulindgen(sp->countChunk())+1UL)
  angles = mod1k->getAngles(sp)
  radRef = mod1k->getMOD06ReflectivityAC(sp, 2)
  radAbs = mod1k->getMOD06ReflectivityAC(sp, 7)
  sz     = cos(!DPI/180.0*angles[*,0])
  radRef = radRef * sz
  radAbs = radAbs * sz
  latlon = mod1k->getLatLon(sp)
  ot06   = mod1k->getMOD06CloudOpticalThickness(sp)
  reff06 = mod1k->getMOD06CloudEffectiveRadius(sp)
;  bt29   = mod1k->getMOD02BrightnessTemperature(sp, 29)
;  bt31   = mod1k->getMOD02BrightnessTemperature(sp, 31)
;  bt32   = mod1k->getMOD02BrightnessTemperature(sp, 32)
  cth    = mod1k->getMOD06CloudTopHeightInfrared(sp)

  rad29  = mod1k->MOD02RadianceB29(sp)
  rad31  = mod1k->MOD02RadianceB31(sp)
  rad32  = mod1k->MOD02RadianceB32(sp)

  ;----
  ; 6 retrieval
  ;--
  ; 6.1 Initialize counters
  cntVR    = 0UL
  cntCRobs = 0UL
  ;--
  ; 6.2 Write headers
  header_c1 = '#PIXID     Lat        Lon       '
  header_c2 = 'Rad29    Rad31    Rad32     CTH     SZA     VZA     SCT  MOD06OT MOD06Reff  '
  header_c = header_c1 + header_c2
  header_d = 'OT     OtSd      Reff   ReffSd    Chisq     Corr12   '
  invfmt = '(I7,X,F9.5,X,F10.5,X,3(F8.3,X),F8.1,X,3(F7.2,X),2(F6.2,X),3(2(F8.4,X,E8.2,X),E11.4,X,F6.3,X))' 
  openw,  ui, outpath+matchupid+'_obs.inv',  /get_lun
  printf, ui, format='(A,$)', header_c
  for i=0, nMdl-1 do printf, ui, format='(A,$)', header_d
  printf, ui, ''

  ;--
  ; 6.3 Retrieval Loop
  retres = fltarr(6, nMdl)
  t = systime(1)
  spd = mod1k->createSelector()
  for i=0UL, sp->countChunk()-1UL do begin

    ; 6.3.1 calculate angles
    sza = angles[i,0]
    vza = angles[i,1]
    raa = angles[i,2]
    sz = cos(!DPI/180.0*sza)
    vz = cos(!DPI/180.0*vza)
    ra = cos(!DPI/180.0*raa)
    ga = 180.0/!DPI*acos( sz*vz + sqrt((1.0-sz*sz)*(1.0-vz*vz)) * ra) ; glitter angles
    sa = 180.0/!DPI*acos(-sz*vz + sqrt((1.0-sz*sz)*(1.0-vz*vz)) * ra) ; scattering angles

    ; 6.3.2 check the validity
    if (ga ge 35.0) and (sa ge 60) and (sa le 160) and $
       (angles[i,0] lt 81) and (angles[i,0] ge 0)  and $
       (latlon[i,1] lt 179.5) and (latlon[i,1] ge -179.5) and $
       (angles[i,1] lt 60) and (angles[i,1] ge 0) then begin
; (latlon[i,1] lt 179.5) and (latlon[i,1] ge -179.5) and $ for SplitWindow retrievals
; Kuo's code has an error, lon can only support [-179.5, 179.5]
  
      ; 6.3.3 prepare angles and variance data to input 
      rad = [radRef[i], radAbs[i]]
      varcov = diag_matrix(err)

      retres[*] = -999.0
      for j=0, nMdl-1 do begin
        ; 6.3.4 build a model and invert observation
        res = [replicate(0.0,4), -1.0, 0.0]
        res = mdl[j]->InvertMLE(sza, vza, raa, rad, varcov)
        if product(finite(res)) eq 0 then res = [replicate(0.0, 4), -1.0, 0.0]
        if res[4] ge 0 then begin
          otv   = spalot[j]->invertX(res[0] < maxspal[j], 1.0)
          ot    = exp(otv[*,0]) - 1.0
          oterr = otv[*,1] / otv[*,0] * (ot + 1.0) * res[1]
          retres[0,j] = [ot, oterr, 0.5*res[2:3], res[4:5]]
        endif
      endfor

      if total(retres[0,*] le -10.0) eq 0 then begin
        cntCRobs = cntCRobs + 1
        printf, ui, format=invfmt, pixid[i], latlon[i,0], latlon[i,1], $
                    rad29[i], rad31[i], rad32[i], cth[i], sza, vza, sa, ot06[i], reff06[i], retres
      endif
      cntVR = cntVR + 1

    endif 
  endfor
  print, systime()+' ['+matchupid+'] retrieval time '+string(format='(F10.5)',systime(1)-t)+' seconds'

  ;--
  ; 6.4 message
  print, systime()+' ['+matchupid+'] retrieval valid       (obs) '+string(format='(I8)',cntVR)+' pixels'
  print, systime()+' ['+matchupid+'] retrieval convergence (obs) '+string(format='(I8)',cntCRobs)+' pixels'
  close, ui
end
