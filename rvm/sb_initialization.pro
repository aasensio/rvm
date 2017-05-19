pro sb_initialization, BASIS, scales, settings, Targets, Alpha, beta, mu, PHI, Used

; A "reasonable" initial value for the noise in the Gaussian case
	GAUSSIAN_SNR_INIT	= 0.01;

; "Reasonable" initial alpha bounds
	INIT_ALPHA_MAX	= 1e3;
	INIT_ALPHA_MIN	= 1e-3;

; BASIS PREPROCESSING:
; Scale basis vectors to unit norm. This eases some calculations and
; will improve numerical robustness later.
	sb_preprocessbasis, BASIS, scales

; Compute beta from the noise standard deviation
	beta = 1.d0 / settings.noiseStdDev^2

; Initialize basis (phi), mu and alpha

; First compute linearised output for use in heuristic initialization
	TargetsPseudoLinear = Targets

; 1) the starting basis, PHI

; Take account of "free basis": it needs to be included from the outset

; Set initial basis to be the largest projection with the targets
	proj = BASIS ## TargetsPseudoLinear

	foo = max(abs(proj), Used)
	print, format='(A,I2)', 'Initialising with maximally aligned basis vector ', Used

	PHI	= BASIS[*,Used]
	M	= n_elements(Used)
 
; 2) the most probable weights
;  mu will be calculated analytically later in the Gaussian case
; 	mu = []

; 3) the hyperparameters

; Exact for single basis function case (diag irrelevant),
; heuristic in the multiple case	
	p = total(PHI*PHI)*beta
	q = total(PHI*Targets)*beta
	Alpha	= [p^2 /(q^2-p)]
	
; Its possible that one of the basis functions could already be irrelevant (alpha<0), so trap that
	ind = where(Alpha lt 0, count)
	if (count eq n_elements(Alpha)) then begin
		print, 'WARNING: no relevant basis function at  initialisation'
	endif
	if (count ne 0) then begin
		Alpha[ind] = INIT_ALPHA_MAX
	endif

	print, format='(A,F9.3)', 'Initial alpha = ', Alpha

end