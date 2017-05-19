@routines
@sparsebayes
@sb_initialization
@sb_preprocessbasis
@sb_fullstatistics
; Do an example of 

pro SparseBayesDemo, noiseToSignal

; Random seed
	seed = 1

; Set default S/N in case it does not exist
	if (not keyword_set(noiseToSignal)) then begin
		noiseToSignal = 0.2
	endif

; Number of points
	N = 100

; Number of basis functions
	M = 100

; Basis width (data will be in [0,1])
	basisWidth = 0.05

; Define probability of a basis function NOT being used by the generative
; model. i.e. if pSparse=0.90, only 10% of basis functions (on average) will
; be used to synthesise the data.
	pSparse = 0.90
	iterations = 500

; Synthetic data generation
	X = findgen(N) / N

; Define the basis
; Locate basis functions at data points
	C = X

	Basis = exp(-distSquared(X,C) / basisWidth^2)
 	Basis = [basis,transpose(replicate(1.d0,100))]
 	Basis = transpose(Basis)

; Temporarily read data from matlab
	w = ddread('../SB2_Release_200/w.txt')
	Outputs = ddread('../SB2_Release_200/output.txt')
	Outputs = Outputs

	M = n_elements(Basis[0,*])

; Sparse Bayes section
	settings = create_struct('noiseStdDev', 1)

	OPTIONS = create_struct('fixednoise', 1, 'diagnosticLevel', 2)
	
	SparseBayes, BASIS, Outputs, options, settings, parameter, hyperparameter, diagnostic

; Compute the prediction
	w_infer = fltarr(M)
	w_infer[parameter.relevant] = parameter.value
	y_pred = transpose(BASIS) ## w_infer

	cwpal
	plot, X, outputs, psym=8
	oplot, X, y_pred, col=2, thick=3
	
	stop

end