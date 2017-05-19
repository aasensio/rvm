@routines
@sparsebayes
@sb_initialization
@sb_preprocessbasis
@sb_fullstatistics

pro bayes, tb, pc

	y = tb
	s = size(pc)

;Basis has to be in the format [Ndata, Nbasis], since we have 
; s[1] PCs of s[2] stars, the stars are the number of different data.

	basis = (transpose(pc))[*,0:15]
	
	basis = [[basis],[replicate(1.d0,234)]]

	noise=100

	outputs=transpose(y)
	Nfine = s[2]
	Basisfine=basis

	Basisfine_normalized = basisfine
	sb_preprocessbasis, Basisfine_normalized, scales

; Number of basis functions that we consider. It can be much larger
; than the number of measured points. In principle, it can even be an
; infinite number of functions if you wish
	M = n_elements(Basis[0,*])

;***********************
; Sparse Bayes section
;***********************

; Set the noise level of the observations. Another option is to infer the noise from the observations
	settings = create_struct('noiseStdDev', noise)

; Some options.
; fixednoise : 0 -> estimate it from the observations / 1 -> fixed noise given by seetings.noiseStdDev
; diagnosticLevel : verbosity level
; iterations : maximum number of iterations
	OPTIONS = create_struct('fixednoise', 1, 'diagnosticLevel', 2, 'iterations', 1000)

; Do the work
	SparseBayes, BASIS, outputs, options, settings, parameter, hyperparameter, diagnostic, SIGMA

; In the output, parameter.relevant contains the elements of the basis that
; are computed to contribute to the signal. The value of the associated coefficient is
; given by parameter.value

; We compute the most probable fit (maximum a-posteriori)
	w_infer = fltarr(M)
	w_infer[parameter.relevant] = parameter.value
	y_pred = transpose(BASIS) ## w_infer

	y_pred_fine = transpose(BASISfine) ## w_infer
  
  
	location= where(w_infer ne 0.,count)
	print, w_infer[parameter.relevant], parameter.relevant, 'N of relevant componentes=' ,count
	print, stddev(y_pred_fine-Tb)
	
; Compute 1-sigma error of the maximum a-posteriori fit
	error = fltarr(Nfine)
	for i = 0, Nfine-1 do begin
		f = Basisfine_normalized[i,parameter.relevant]
		error[i] = sqrt(settings.noiseStdDev^2 + transpose(f)##(SIGMA ## f))
	endfor

	plot, tb, y_pred, psym=3
	stop
end


pro dobayes
	restore,'datos.idl'
	bayes, tb, pc
end
