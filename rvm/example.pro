@routines
@sparsebayes
@sb_initialization
@sb_preprocessbasis
@sb_fullstatistics

pro example

; Read the data
	nlines = file_lines('stenflo_mg.dat')
	openr,2,'stenflo_mg.dat'
	dat = dblarr(2,nlines)
	readf,2,dat
	close,2

; Number of observed points	
	N = n_elements(dat[0,*])

; X is the wavelength in our case
	X = reform(dat[0,*])

; Additionally, generate a fine axis so that the final plots are visually better
	Xfine = findgen(100) / 99.d0 * (max(X)-min(X)) + min(X)
	Nfine = n_elements(Xfine)

; Measured intensity. For some stupid reasons that I have not fixed,
; it has to be in dimensions [1,Ndata]
	Outputs = transpose(reform(dat[1,*]))

; The basis will consist of Gaussians centered at the observed points
; with a fixed width. We set the width to 2 A
	basisWidth = 2.d0

; Generate the Gaussians using a fast routine the distances. In any case,
; it can be done using loops
	Basis = exp(-distSquared(X,X) / basisWidth^2)
	Basisfine = exp(-distSquared(X,Xfine) / basisWidth^2)

; Add a constant bias to the set of basis functions
	Basis = [Basis,replicate(1.d0,1,N)]
	Basisfine = [Basisfine,replicate(1.d0,1,Nfine)]

; Transpose the basis matrix. It has to be in dimensions [Ndata, Nbasis]
	Basis = transpose(Basis)
	Basisfine = transpose(Basisfine)

; Compute a normalized version of the basis vectors of the prediction
	Basisfine_normalized = Basisfine
	sb_preprocessbasis, Basisfine_normalized, scales

; Number of basis functions that we consider. It can be much larger
; than the number of measured points. In principle, it can even be an
; infinite number of functions if you wish
	M = n_elements(Basis[0,*])

;***********************
; Sparse Bayes section
;***********************

; Set the noise level of the observations. Another option is to infer the noise from the observations
	noise = 0.009d0
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

	y_pred_fine = transpose(BASISfine) ## w_infer

; Compute 1-sigma error of the maximum a-posteriori fit
	error = fltarr(Nfine)
   for i = 0, Nfine-1 do begin
		f = Basisfine_normalized[i,parameter.relevant]
		error[i] = sqrt(settings.noiseStdDev^2 + transpose(f)##(SIGMA ## f))
	endfor

; Do some plots
	cwpal
	plot, X, outputs, psym=8, xtit='Wavelength [A]', ytit='Q/I', yran=[-0.05,0.05], $
		tit='N='+strtrim(string(n_elements(parameter.relevant),FORMAT='(I1)'),2)
	oplot, X, y_pred, col=2, thick=3
	oplot, Xfine, y_pred_fine, col=3, thick=3
	errplot, X, outputs-settings.noiseStdDev, outputs+settings.noiseStdDev

	oplot, Xfine, y_pred_fine-error, col=2, thick=2,line=1
	oplot, Xfine, y_pred_fine+error, col=2, thick=2,line=1
	
	stop

end
