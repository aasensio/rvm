@rvm/routines
@rvm/sparsebayes
@rvm/sb_initialization
@rvm/sb_preprocessbasis
@rvm/sb_fullstatistics

pro rvm, x, y, noises, nPoly

; Number of observed points	
	N = n_elements(x)

; Measured intensity. For some stupid reasons that I have not fixed,
; it has to be in dimensions [1,Ndata]
	Outputs = transpose(y)

; Generate the Gaussians using a fast routine the distances. In any case,
; it can be done using loops
	Basis = fltarr(nPoly,N)
	for i = 0, nPoly-1 do begin
		Basis[i,*] = x^i
	endfor

; Transpose the basis matrix. It has to be in dimensions [Ndata, Nbasis]
	Basis = transpose(Basis)
	
; Number of basis functions that we consider. It can be much larger
; than the number of measured points. In principle, it can even be an
; infinite number of functions if you wish
	M = n_elements(Basis[0,*])

;***********************
; Sparse Bayes section
;***********************

; Set the noise level of the observations. Another option is to infer the noise from the observations
	noise = mean(noises)
	settings = create_struct('noiseStdDev', noise)

; Some options.
; fixednoise : 0 -> estimate it from the observations / 1 -> fixed noise given by seetings.noiseStdDev
; diagnosticLevel : verbosity level
; iterations : maximum number of iterations
	OPTIONS = create_struct('fixednoise', 1, 'diagnosticLevel', 2, 'iterations', 200)
	
; Do the work
	SparseBayes, BASIS, Outputs, options, settings, parameter, hyperparameter, diagnostic, SIGMA

; In the output, parameter.relevant contains the elements of the basis that
; are computed to contribute to the signal. The value of the associated coefficient is
; given by parameter.value

; We compute the most probable fit (maximum a-posteriori)
; w_infer contains the weights of the basis functions and should be a sparse vector by construction
	w_infer = fltarr(M)
	w_infer[parameter.relevant] = parameter.value
	
	y_pred = transpose(BASIS) ## w_infer
	
	print, 'RVM: ', w_infer

end
