@routines
@sparsebayes
@sb_initialization
@sb_preprocessbasis
@sb_fullstatistics

pro stenflo_mg

; Random seed
; 	seed = 1

;**********************************************
; Compute first the case at mu=0.1
;**********************************************
	dat = ddread('stenflo_mg.dat')
	N = n_elements(dat[0,*])
	X = reform(dat[0,*])
	Outputs = transpose(reform(dat[1,*]))

	nwidth = 11
	widths = findgen(nwidth) / (nwidth-1.d0) * (5.d0-1.d0) + 1.d0

	lnZ_noise = fltarr(nwidth)
	lnZ_signal = fltarr(nwidth)

	!p.multi = [0,3,4]

	abre_ps,'stenflo_mu0.1.ps',/todo, /color

	for i = 0, nwidth-1 do begin

; Basis width (data will be in [0,1])
		basisWidth = widths[i]

		iterations = 500

; Define the basis
; Locate basis functions at data points
		C = X

; The basis functions are Gaussians centered on the location of the points
; and a constant bias
		Basis = exp(-distSquared(X,C) / basisWidth^2)
		Basis = [basis,replicate(1.d0,1,N)]
 		Basis = transpose(Basis)

		M = n_elements(Basis[0,*])

; Sparse Bayes section
		settings = create_struct('noiseStdDev', 0.018/2.d0);0.018)

		OPTIONS = create_struct('fixednoise', 1, 'diagnosticLevel', 2)

		SparseBayes, BASIS, Outputs, options, settings, parameter, hyperparameter, diagnostic

; Compute the prediction
		w_infer = fltarr(M)
		w_infer[parameter.relevant] = parameter.value
		y_pred = transpose(BASIS) ## w_infer

 		cwpal
 		plot, X, outputs, psym=8, xtit='Wavelength [A]', ytit='Q/I', yran=[-0.05,0.05], tit='Width='+strtrim(string(widths[i],FORMAT='(F3.1)'),2)+$
			' - N='+strtrim(string(n_elements(parameter.relevant),FORMAT='(I1)'),2)
 		oplot, X, y_pred, col=2, thick=3
 		errplot, X, outputs-settings.noiseStdDev, outputs+settings.noiseStdDev
 
 		for j = 0, M-1 do begin
 			oplot, X, basis[*,j]*w_infer[j], col=3
 		endfor

		lnZ_noise[i] = -N/2.d0*alog(2.d0*!DPI) - N * alog(settings.noiseStdDev) - 0.5*total(Outputs^2) / settings.noiseStdDev^2
		lnZ_signal[i] = max(diagnostic.logmarginal)
		
	endfor

	plot, widths, lnZ_signal - lnZ_noise, xtit='Basis width [A]',ytit='ln (Z!dsignal!n/Z!dnoise!n)'

	!p.multi = 0

	cierra_ps


;**********************************************
; Compute then the case at mu=1
;**********************************************
	dat = ddread('stenflo_mg_center.dat')
	N = n_elements(dat[0,*])
	X = reform(dat[0,*])
	Outputs = transpose(reform(dat[1,*]))

	nwidth = 11
	widths = findgen(nwidth) / (nwidth-1.d0) * (5.d0-1.d0) + 1.d0

	lnZ_noise = fltarr(nwidth)
	lnZ_signal = fltarr(nwidth)

	!p.multi = [0,3,4]

	abre_ps,'stenflo_mu1.0.ps',/todo, /color

	for i = 0, nwidth-1 do begin

; Basis width (data will be in [0,1])
		basisWidth = widths[i]

		iterations = 500

; Define the basis
; Locate basis functions at data points
		C = X

; The basis functions are Gaussians centered on the location of the points
; and a constant bias
		Basis = exp(-distSquared(X,C) / basisWidth^2)
		Basis = [basis,transpose(replicate(1.d0,N))]
		Basis = transpose(Basis)
		

		M = n_elements(Basis[0,*])

; Sparse Bayes section
		settings = create_struct('noiseStdDev', 0.018/2.d0)

		OPTIONS = create_struct('fixednoise', 1, 'diagnosticLevel', 2)

		SparseBayes, BASIS, Outputs, options, settings, parameter, hyperparameter, diagnostic

	; Compute the prediction
		w_infer = fltarr(M)
		w_infer[parameter.relevant] = parameter.value
		y_pred = transpose(BASIS) ## w_infer

 		cwpal
 		plot, X, outputs, psym=8, xtit='Wavelength [A]', ytit='Q/I', yran=[-0.05,0.05], tit='Width='+strtrim(string(widths[i],FORMAT='(F3.1)'),2)+$
			' - N='+strtrim(string(n_elements(parameter.relevant),FORMAT='(I1)'),2)
 		oplot, X, y_pred, col=2, thick=3
 		errplot, X, outputs-settings.noiseStdDev, outputs+settings.noiseStdDev

 		for j = 0, M-1 do begin
 			oplot, X, basis[*,j]*w_infer[j], col=3
 		endfor

		lnZ_noise[i] = -5.d0*alog(2.d0*!DPI) - 10.d0 * alog(settings.noiseStdDev) - 0.5*total(Outputs^2) / settings.noiseStdDev^2
		lnZ_signal[i] = max(diagnostic.logmarginal)
	endfor

	plot, widths, lnZ_signal - lnZ_noise, xtit='Basis width [A]',ytit='ln (Z!dsignal!n/Z!dnoise!n)'

	!p.multi = 0

	cierra_ps

	stop

end
