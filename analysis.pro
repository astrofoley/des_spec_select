PRO analysis, card, mag

;readcol, 'out_temp', format='X,D,X,X,X,I,X,D,X,D,X,A,X,D,X,D,X,D,D,D,D,D,D,D,D,D,D,D,X,D', peak_r, cid, z, phase, name, current_r, z_err, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, ia_prob
;readcol, 'out_temp', format='X,D,X,X,X,X,X,X,X,X,X,X,X,X,X,D,X,D,X,X,X,A,X,D,D,D,D,D,D,D,D,D,D,D,X,D,X,X,X,I,X,D,X,D', z, current_r, phase, name, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, peak_r, cid, ia_prob, z_err
readcol, 'out_temp3', format='X,D,X,X,X,X,X,X,X,X,X,X,X,X,X,D,X,D,X,X,X,A,X,D,D,D,D,D,D,D,D,D,D,D,X,D,X,X,X,L,X,D,X,D', z, current_r, phase, name, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, peak_r, cid, ia_prob, z_err





plothist, current_r, bin=0.1
plothist, peak_r, bin=0.1
plothist, z, bin=0.1

print, median(z)

print, n_elements(where(z le 0.2))
print, n_elements(where(peak_r le 20.5))
print, n_elements(where(current_r le 19))
print, n_elements(where((z le 0.2) or (peak_r le 20.5) or (current_r le 19)))

plothist, ia_prob, bin=0.1
print, median(ia_prob)
print, n_elements(where(ia_prob gt 0.8))

;plothist, p1, bin=1
;plothist, p2, bin=1
;plothist, p3, bin=1
;plothist, p4, bin=1
;plothist, p5, bin=1
;plothist, p6, bin=1
;plothist, p7, bin=1
;plothist, p8, bin=1
;plothist, p9, bin=1
;plothist, p10, bin=1
;plothist, p11, bin=1

print, n_elements(where(p6 gt 3))
print, n_elements(where((p6 gt 3) and (peak_r gt 20.5)))
print, n_elements(where((p6 gt 3) and (ia_prob gt 0.7)))

plot, ia_prob, p9, psym=3

case mag of
  20.0: ii = reverse(sort(p1))
  20.5: ii = reverse(sort(p2))
  21.0: ii = reverse(sort(p3))
  21.5: ii = reverse(sort(p4))
  22.0: ii = reverse(sort(p5))
  22.5: ii = reverse(sort(p6))
  23.0: ii = reverse(sort(p7))
  23.5: ii = reverse(sort(p8))
  24.0: ii = reverse(sort(p9))
  24.5: ii = reverse(sort(p10))
  25.0: ii = reverse(sort(p11))
  else: print, 'wrong limiting magnitude'
endcase

case mag of
  20.0: pp = p1
  20.5: pp = p2
  21.0: pp = p3
  21.5: pp = p4
  22.0: pp = p5
  22.5: pp = p6
  23.0: pp = p7
  23.5: pp = p8
  24.0: pp = p9
  24.5: pp = p10
  25.0: pp = p11
  else: print, 'wrong limiting magnitude'
endcase

;for i = 0, n_elements(ii) - 1 do $
;for i = 0, 250 do begin
;  if (((card eq 'N') or (card eq 'n')) and (strmid(name[ii[i]],5,1) ne 'E')) then $
;  print, name[ii[i]], pp[ii[i]], z[ii[i]], z_err[ii[i]], ia_prob[ii[i]], current_r[ii[i]], '  http://dessne.cosmology.illinois.edu/SNWG/web/display/examineCand_Y2.php?Name='+name[ii[i]]
;endfor

;https://portal-auth.nersc.gov/atc2/web/targets/'+name[ii[i]]

;for i = 0, n_elements(ii)-1 do begin
for i = 0, n_elements(ii)-1 do begin
  if ((pp[ii[i]] gt 1d-15) and (strmid(name[ii[i]],5,1) eq 'C1')) then $
;  if (pp[ii[i]] gt 1d-15) then $
  print, name[ii[i]], pp[ii[i]], z[ii[i]], z_err[ii[i]], ia_prob[ii[i]], current_r[ii[i]], '  http://dessne.cosmology.illinois.edu/SNWG/web/display/examineCand_Y2.php?Name='+name[ii[i]]
endfor

;for i = 0, 150 do print, cid[ii[i]]

;ii = where((name eq 'DES14X3qhy') or (name eq 'DES14C3rad') or (name eq 'DES14C3plr'))
;i = 0

;print, ' '
;  print, name[ii[i]], pp[ii[i]], z[ii[i]], z_err[ii[i]], ia_prob[ii[i]], current_r[ii[i]], '  http://dessne.cosmology.illinois.edu/SNWG/web/display/examineCand_Y2.php?Name='+name[ii[i]]


stop


END
