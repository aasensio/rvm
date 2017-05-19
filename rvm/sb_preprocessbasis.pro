; Simply normalises the basis vectors to unit length, returning the
; original lengths so that the weights can be rescaled appropriately
; before returning to the user.

pro sb_preprocessbasis, BASIS, scales

 	scales = sqrt(total(BASIS^2, 1))

; Work-around divide-by-zero
 	ind = where(Scales eq 0, count)
 	if (count ne 0) then begin
		scales[ind] = 1.d0		
 	endif

 	for i = 0, n_elements(BASIS[0,*])-1 do begin
		BASIS[*,i] = BASIS[*,i] / scales[i]
	endfor	
end