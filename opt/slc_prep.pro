pro slc_prep
;-- based on LC_PREP.pro from FBEYE
;   take in list of short (and long) cadence files from Kepler
;   output single giant LC, stichted together
;   output files for each month also

;-- run on laptop in dir /Users/james/data/kepler17_data/

dir = "/Users/james/data/kepler17_data/"

fluxmed =  27283.5761719d0 ; found from median of LLC data
BJDREF = 2454833d0

spawn, 'ls ' + dir + '*slc.fits > slc.lis'
spawn, 'ls ' + dir + '*llc.fits > llc.lis'

;-- start w/ LLC data
readcol,'llc.lis',lfile,f='(A)',/silent

time = [-99]
flux = [-99]
ferr = [-99]
for n=0l,n_elements(lfile)-1 do begin
   t = MRDFITS( lfile[n], 1, /silent, hdr)
   time = [time, t.time]
   ferr = [ferr, real(t.pdcsap_flux_err,-99.)]

   ftmp = real(t.pdcsap_flux, -999.-fluxmed)
;   x = WHERE(ftmp gt 0)
;   fit = POLY_FIT((t.time)[x], ftmp[x],1)
;   ftmp = ftmp - POLY(t.time, fit) + fluxmed
   ftmp = ftmp - median(ftmp) + fluxmed

   flux = [flux, ftmp]
endfor

;-- remove bad data points
remove,where(flux lt -90), time, flux, ferr

s = sort(time)
time = time[s]
flux = flux[s]
ferr = ferr[s]

plot,time,flux,/xsty,/ysty,psym=3,xtitle='Time',ytitle='Flux'


;; tmp = LNP_TEST(time, flux, /double, wk1 = freq, wk2 = pwr, $
;;                jmax = jj, ofac = 16, hifac = 2)
;; plot,1./freq, pwr,xrange=[10,15],xtitle='Period (days)',ytitle='Power'
;; p_rot = 1./freq[jj]
;--- this is the rotation period I will use, have put in the param file
;; p_rot = 12.25817669188d0 

readcol,'~/Dropbox/python/kepler17/kep17.params',params,f='(D)',/silent
p_rot = params[3]
t_dur = params[2]
t_0 = params[1]
p_orb = params[0]

t_0 = t_0 - p_orb * 100.


set_plot,'ps'
device,filename='kepler17_diffrot_whole.eps',$
       /encap,/color,/inch,xsize=9,ysize=6
cubehelix,rot=-0.5,start=0.3

xbinx = p_rot*2.
ybiny = 0.05

pixel_plus, time, (time mod p_rot)/p_rot, flux/1d4,$
            xbin=xbinx,ybin=ybiny,lvl=[findgen(20)*.01+2.64],$
            /xsty,/ysty,xrange=[100,1615],yrange=[-0.5,1.5],$
            xtitle='!7BJD - 2454833.11567 (days)',$
            ticklen=-.02,ytitle='!7Phase',font=0,$
            position=posgen(1,5,2,ysp=4)

pixel_plus, time, (time mod p_rot)/p_rot-1, flux/1d4,$
            xbin=xbinx,ybin=ybiny,lvl=[findgen(20)*.01+2.64],$
            /xsty,/ysty,xrange=[100,1615],yrange=[-0.5,1.5],/overplot

pixel_plus, time, (time mod p_rot)/p_rot+1, flux/1d4,$
            xbin=xbinx,ybin=ybiny,lvl=[findgen(20)*.01+2.64],$
            /xsty,/ysty,xrange=[100,1615],yrange=[-0.5,1.5],/overplot


plot,time,(flux-fluxmed)/fluxmed,psym=8,symsize=0.2,/noerase,position=posgen(1,5,1),xtickname=replicate(' ',8),ytitle='!7Flux',font=0,charsize=1,xrange=[100,1610],/xsty,/ysty

device,/close
set_plot,'X'




n=0
plot,time,flux,/xsty,/ysty,xrange=[-1., 1.]*t_dur + t_0 - BJDREF + p_orb*n,psym=-4



;-- now the SLC data
readcol,'slc.lis',sfile,f='(A)',/silent

time_s = [-99]
flux_s = [-99]
ferr_s = [-99]
for n=0l,n_elements(sfile)-1 do begin
   ts = MRDFITS( sfile[n], 1, /silent, hdr)
   time_s = [time_s, ts.time]
   ferr_s = [ferr_s, real(ts.pdcsap_flux_err,-99.)]

   ftmp = real(ts.pdcsap_flux, -999.-fluxmed)
;   x = WHERE(ftmp gt 0)
;   fit = POLY_FIT((t.time)[x], ftmp[x],1)
;   ftmp = ftmp - POLY(t.time, fit) + fluxmed
   ftmp = ftmp - median(ftmp) + fluxmed

   flux_s = [flux_s, ftmp]
endfor

;-- remove bad data points
remove,where(flux_s lt -90), time_s, flux_s, ferr_s

medflux = [-99]
medtime = [-99]
;-- stitch in SLC data in the transit regions
isdone = 0.
n=0
while isdone lt 1 do begin
   if (t_dur + t_0 - BJDREF + p_orb*n) gt max(time) then $
      isdone = 1.
   xl = where(time ge (-1.*t_dur + t_0 - BJDREF + p_orb*n) and $
              time le (1. * t_dur + t_0 - BJDREF + p_orb*n))
   xs = where(time_s ge (-0.8 * t_dur + t_0 - BJDREF + p_orb*n) and $
              time_s le ( 0.8 * t_dur + t_0 - BJDREF + p_orb*n))

   if (xl[0] ne -1) and (xs[0] ne -1) then begin
      xs1 = where(time_s ge (-2. * t_dur + t_0 - BJDREF + p_orb*n) and $
                  time_s le (-0.8 * t_dur + t_0 - BJDREF + p_orb*n))
      xs2 = where(time_s le (2. * t_dur + t_0 - BJDREF + p_orb*n) and $
                  time_s ge (0.8 * t_dur + t_0 - BJDREF + p_orb*n))

      xl1 = where(time ge (-2. * t_dur + t_0 - BJDREF + p_orb*n) and $
                  time le (-0.8 * t_dur + t_0 - BJDREF + p_orb*n))
      xl2 = where(time le (2. * t_dur + t_0 - BJDREF + p_orb*n) and $
                  time ge (0.8 * t_dur + t_0 - BJDREF + p_orb*n))
      

      fit_s = poly_fit(time_s[[xs1,xs2]], flux_s[[xs1,xs2]], 1)
      fit_l = poly_fit(time[[xl1,xl2]], flux[[xl1,xl2]], 1)

      ;; ms = median(flux_s[xs])
      ;; ml = median(flux[xl])

      ;; ml = interpol(flux[[xl1,xl2]], time[[xl1,xl2]], time_s[xs])
      ;; ms = interpol(flux_s[[xs1,xs2]], time_s[[xs1,xs2]], time_s[xs])

      ;remove, xl, time, flux, ferr
      time = [time, time_s[xs]]
      ferr = [ferr, ferr_s[xs]]
      flux = [flux, flux_s[xs] - POLY(time_s[xs],fit_s)+POLY(time_s[xs],fit_l)]



      ;; hide in this loop: gather data to
      ;;     create median model of all transits stacked

      medflux = [medflux, flux_s[xs] - POLY(time_s[xs],fit_s) + fluxmed]
      medtime = [medtime, time_s[xs]]

   endif
   n = n+1
endwhile


s = sort(time)
time = time[s]
flux = flux[s]
ferr = ferr[s]



;; n=900
;; n++ & plot,time,flux,/xsty,/ysty,xrange=[-1.5, 1.5]*t_dur + t_0 - BJDREF + p_orb*n,psym=4

set_plot,'ps'
device,filename='kepler17_datasample.eps',/encap,/inch,xsize=10,ysize=5
plot,time, (flux - median(flux))/median(flux), $
     psym=8, symsize=0.5, xrange=[406,408],/ysty,font=0,charsize=1.5,$
     xtitle='!7Time (BJD - 2454833 days)', ytitle='!7Relative Flux'
device,/close


;-- output datafiles for STSP, (or STSPy)
forprint,textout='kepler17_whole.dat',time,flux,ferr,f='(D,D,D)',/nocomm

;-- shift slightly to the left of starting at 868, avoid missing transit
x868 = where(time ge 868-p_rot*2./3. and time le 868+p_rot/3.)
forprint,textout='kepler17_868.dat',time[x868],flux[x868],ferr[x868],f='(D,D,D)',/nocomm


loadct,39,/silent

device,filename='kepler17_869.eps',/encap,/inch,xsize=10,ysize=5
ploterror,time[x868],smooth(flux[x868],1),ferr[x868],$
          /ysty,psym=3,xrange=[869.5,870]
oplot,time[x868],smooth(flux[x868],5),color=50,thick=5
device,/close





remove,0,medflux,medtime

phz = ((medtime- t_0 + BJDREF - p_orb/2.) mod p_orb) / p_orb-0.5
s1 = sort(phz)
s2 = sort(medtime[s1])

medbin, (phz[s1]), (smooth(medflux,5,/edge))[s1], $
        xbin, ybin, .0005, min(phz), max(phz), std=stdbin

set_plot,'X'

cubehelix,rot=0,start=0
contour_plus, ((medtime- t_0 + BJDREF - p_orb/2.) mod p_orb) / p_orb-0.5, smooth(medflux,5), psym=3,/xsty,/ysty,1,xbin=.0002,ybin=.001d4,/pix,/rev

loadct,39,/silent
oplot,xbin,ybin,color=250
oplot,xbin,ybin+stdbin,color=250,linestyle=2
oplot,xbin,ybin-stdbin,color=250,linestyle=2

model = (interpol(ybin,xbin,phz[s1]))[s2]

;; ok = where(medtime ge 400 and medtime le 420)
;; forprint, textout='kepler17_mediantransit.dat',medtime[ok],model[ok],model[ok]*1d-5,$
forprint, textout='kepler17_mediantransit.dat',medtime,model,model*1d-5,$
          f='(D,D,D)',/silent,/nocomm

ybin = (ybin - fluxmed)/fluxmed
stdbin = stdbin/fluxmed

set_plot,'ps'
cubehelix,rot=0,start=0


device,filename='kepler17_alltransits.eps',/encap,/color,/inch,xsize=7,ysize=5

contour_plus,((medtime- t_0 + BJDREF - p_orb/2.) mod p_orb) / p_orb-0.5, $
             (smooth(medflux,5)-fluxmed)/fluxmed,font=0,$
             /xsty,/ysty,1,xbin=.0004,ybin=.0002,/pix,/rev,$
             xtitle='!7Orbital Phase (P = 1.486 days)',ytitle='!7Relative Flux'
loadct,39,/silent
oplot,xbin,ybin,color=250,thick=3
oplot,xbin,ybin+stdbin,color=250,linestyle=2,thick=3
oplot,xbin,ybin-stdbin,color=250,linestyle=2,thick=3
device,/close




t_leslie = [554.73444, 568.10582, 579.9915, 591.87721, $
            603.76275, 615.64890, 627.53428]
set_plot,'X'




t_pk = 568.103 ;-- from Leslie
;plot, time, flux, xrange=[-t_dur,t_dur] + t_pk,/xsty,/ysty


set_plot,'ps'
device,filename='a_spot_evol.eps',/encap,/color,/inch,xsize=4,ysize=10
plot,[0],/nodata,xsty=5,ysty=5
for n=-4,4 do begin

   ;; plot,time, flux, xrange=[-0.5,0.5]*p_rot + t_pk + (8.*p_orb)*n,$
   ;;      /noerase,position=posgen(3,12,2,n+5,xspa=2),xsty=5,ysty=5

   nn = where(time ge -t_dur+ t_pk + (8.*p_orb)*n and $
              time le t_dur+ t_pk + (8.*p_orb)*n)
   plot, time, flux, xrange=[-t_dur,t_dur] + t_pk + (8.*p_orb)*n,$
         /noerase,position=posgen(1,9,1,n+5),xsty=5,ysty=5
   ;; oplot,xbin*p_orb + t_pk + (8.*p_orb)*n, $
   ;;       ybin*fluxmed+fluxmed,color=250

endfor
device,/close

t_pk = 869.707

device,filename='b_spot_evol.eps',/encap,/color,/inch,xsize=4,ysize=10
plot,[0],/nodata,xsty=5,ysty=5
for n=-4,4 do $
   plot, time, flux, xrange=[-t_dur,t_dur] + t_pk + (8.*p_orb)*n,$
         /noerase,position=posgen(1,9,1,n+5),xsty=5,ysty=5
device,/close


loadct,0,/silent
device,filename='bump.eps',/encap,/color,/inch,xsize=9,ysize=5
plot, time, (flux-fluxmed)/fluxmed, xrange=[-t_dur,t_dur] + t_pk,$
      thick=3,xtitle='!7Time',ytitle='!7Flux',font=0,ysty=5,xsty=5,/nodata
oploterror,time,(flux-fluxmed)/fluxmed,ferr/fluxmed,color=150,thick=3,/nohat
oplot,time,(smooth(flux,5)-fluxmed)/fluxmed,thick=5
device,/close






n=10
t = MRDFITS( lfile[n], 1, /silent, hdr)

ttmp = t.time
ftmp = real(t.pdcsap_flux,-99)
remove,where(ftmp lt 0),ttmp,ftmp



device,filename='llc.eps',/encap,/color,/inch,xsize=7,ysize=5
plot,ttmp, (ftmp-median(ftmp))/median(ftmp),psym=8,symsize=0.5,$
     xtitle='!7Time',ytitle='!7Relative Flux',font=0,charsize=1.5,$
     yrange=[-.035,.02],/ysty,/xsty

device,/close


set_plot,'X'
stop
end

