pro plot_harmonic,infile,lmin,lmax,rootname


if n_elements(rootname) then print,'Using ',rootname,' as root name' $
  else rootname = ''

set_plot,'ps'
loadct,40
;!p.charsize=1.0
;!p.thick=2
;!p.charthick=2

W = 0
B = 50
G = 150
R = 250

nalm  = (lmax+1)^2 - lmin^2
nalm2 = nalm - (nalm - (lmax+1-lmin))/2
print,'nondegenerate alms for one mode:',nalm2

;read in the complex matrix
m  = dcomplexarr(3*nalm2,3*nalm2)
m2 = dcomplexarr(3*nalm2,3*nalm2)
openr,unit,infile,/get
readu,unit,m
close,unit


; rearrange to George's style
for i = 0, 3*nalm2-1 do begin
    i2lm,i,ell,emm
    i_new  = lm2j(ell,emm,lmin,lmax)
    i_new += floor(i/nalm2)*nalm2
    for j = 0, 3*nalm2-1 do begin
        i2lm,j,ell,emm
        j_new  = lm2j(ell,emm,lmin,lmax)
        j_new += floor(j/nalm2)*nalm2
        m2[i_new,j_new] = m[i,j]
    endfor
endfor

; normalize the matrix to appear as correlation factor
for i = 0, 3*nalm2-1 do begin
    isigma = m2[i,i]^0.5
    for j = 0, 3*nalm2-1 do begin
        if i eq j then continue
        jsigma = m2[j,j]^0.5
        m2[i,j] /= isigma*jsigma
    endfor
endfor

; mask out the uninteresting diagonal
for i = 0, 3*nalm2-1 do begin
    m2[i,i] = 0
endfor

; get tick boundaries
nticks = 0
case lmax of
    (lmax < 20): begin
        ticks    = dblarr(3*(lmax+1)+1)
        for isig = 0,2 do begin
            for emm = 0, lmax do begin
                ell = max([lmin,emm])
                ticks[nticks++] = lm2j(ell,emm,lmin,lmax)+isig*nalm2
            endfor
        endfor
        ticks[nticks] = 3*nalm2-1
    end
    else: begin
        nticks   = 0
        interval = floor(3*(lmax-lmin+1)/60.0) ; maximum of 60 ticks
        nskipped = 0
        ticks    = dblarr(3*(lmax+1)/(interval+1)+1)
        for i = 0, 3*nalm2-1 do begin
            if i mod nalm2 eq 0 then begin
                emm = -1
                ell = lmax
            endif
            ell++
            if ell gt lmax then begin
                emm++
                ell = max([lmin, emm])
                if nskipped eq interval then begin
                    nskipped = 0
                    ticks[nticks++] = i
                endif else nskipped++
            endif
        endfor
        ticks[nticks] = 3*nalm2-1
    end
endcase

;print,nticks
;print,ticks

; Make a lego plot of the covariance matrix
device,/encapsulated,/color,filename=rootname+'_harmonic2.eps',bits=8
surface,abs(m2),xstyle=1,ystyle=1,/lego,$
  xtickformat='j2lm',ytickformat='j2lm',$
  xticks=nticks,xtickv=ticks,yticks=nticks,ytickv=ticks,$
  xtickunits=['numeric'],ytickunits=['numeric'],$
  title='Madam 1.25s',/save,$
  xtitle='m-mode',ytitle='m-mode'

; first plot contour over lego plot
ncontours = 50
contour,abs(m2),findgen(3*nalm2),findgen(3*nalm2),$
  /t3d,/noerase,$
  levels=max(m2)*findgen(ncontours)/ncontours,$
  zvalue=1.0,/noclip,$
  xtickformat='j2lm',ytickformat='j2lm',$
  xticks=nticks,xtickv=ticks,yticks=nticks,ytickv=ticks,$
  xtickunits=['numeric'],ytickunits=['numeric'],$
  c_colors=[float(R-B)*findgen(ncontours)/ncontours+B],$
  xticklen=1.0,yticklen=1.0,xtitle='m-mode',ytitle='m-mode'
device,/close

; then plot contour separately to its own file
device,/encapsulated,/color,filename=rootname+'_harmonic2_contour.eps',bits=8
contour,abs(m2),findgen(3*nalm2),findgen(3*nalm2),$
  levels=max(m2)*findgen(ncontours)/ncontours,$
  xtickformat='j2lm',ytickformat='j2lm',$
  xticks=nticks,xtickv=ticks,yticks=nticks,ytickv=ticks,$
  xtickunits=['numeric'],ytickunits=['numeric'],$
  c_colors=[float(R-B)*findgen(ncontours)/ncontours+B],$
  xticklen=1.0,yticklen=1.0,xtitle='m-mode',ytitle='m-mode'

plot,[0,1],/noerase,/nodata,position=[0.7,0.25,0.95,0.3],$
  xticks=1,yticks=1,xtickname=[' ',' '],ytickname=[' ',' ']
nlines=100
for i=0,nlines-1 do begin
    oplot,[1,1]*float(i)/nlines,[0,1],color=fix(float(R-B)*i/nlines+B),$
      thick=8
endfor
xyouts,0.1,0.3,'0.0'
xyouts,0.7,0.3,string(max(abs(m2)),format='(f4.2)')

device,/close

set_plot,'X'
!p.multi = 0

end



function ind2lm,axis,tick_index,value;,level
level=1
i = value mod scope_varfetch('nalm2',level=-1)
ell = 1
m   = 1
for j = 0, i do begin
    m = m + 1
    if m gt ell then begin
        ell++
        m = 0
    endif
endfor

case level of 
    0 : return, string(m, format='(i)')
    1 : return, string(ell, format='(i)')
endcase

end

pro i2lm,index,ell,m
i   = index mod scope_varfetch('nalm2',level=-1)
ell = 1
m   = 1
for j = 0, i do begin
    m = m + 1
    if m gt ell then begin
        ell++
        m = 0
    endif
endfor
end

function j2lm,axis,tick_index,index,level
j   = index mod scope_varfetch('nalm2',level=-1)
ell = scope_varfetch('lmin',level=-1) - 1
m   = 0
for i = 0, j do begin
    ell = ell + 1
    if ell gt scope_varfetch('lmax',level=-1) then begin
        m++
        ell = max([scope_varfetch('lmin',level=-1),m])
    endif
endfor
case level of 
    0 : if   m gt 10 then return,'' else return, string(m, format='(i)')
    1 : if ell gt 10 then return,'' else return, string(ell, format='(i)')
endcase
end

function lm2j,ell,m,lmin,lmax
; returns index in Efstathiou ordering, assumes lmin = 2
j = m*(lmax+2)-m*(m+1)/2+ell-m-3
if m eq 0 then j++
return,j
end

function lm2i,ell,m,lmin,lmax
; returns index in Healpix ordering
return,ell*(ell+1)/2+lmin*(lmin+1)/2+m
end
