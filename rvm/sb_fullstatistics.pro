; OUTPUT ARGUMENTS:
;
;	SIGMA			Posterior covariance matrix for relevant bases
;	MU				Posterior mean
;	S_IN			S-factors for in-model (relevant) basis vectors
;	Q_IN			Q-factors for in-model (relevant) basis vectors
;	S_OUT			S-factors for all basis vectors
;	Q_OUT			Q-factors for all basis vectors
;	FACTOR			Q^2-S relevance factors
;	LOGML			Log-marginal-likelihood
;	GAMMA			"Well-determinedness" factors
;	BETABASIS_PHI	Cached value of BASIS'*B*PHI matrix
;	BETA			Inverse noise variance (vector of beta
;					approximations in non-Gaussian case)
;
; INPUT ARGUMENTS:
;
;	LIKELIHOOD		LIKELIHOOD structure
;	BASIS			Full NxM basis matrix
;	PHI				Current relevant basis matrix
;	TARGETS			N-vector with target output values
;	USED			Relevant basis vector indices
;	ALPHA			Hyperparameters
;	BETA			Inverse noise variance (ignored in non-Gauss case)
;	MU				Current posterior mean (for non-Gauss)
;	BASIS_PHI		Cached value of BASIS'*PHI (for Gauss only)
;	BASIS_TARGETS	Cached value of BASIS'*TARGETS (for Gauss only)
; This function computes the posterior, and other, statistics for the
; SPARSEBAYES algorithm in "long-hand" fashion. i.e. it does not use the
; more effecient "iterative" updates.
; It is required on every iteration (unless approximating) in the
; non-Gaussian case, and on iterations in the Gaussian case where the noise
; estimate (BETA) is updated.

pro sb_fullstatistics, BASIS, PHI, Targets, Used, Alpha, beta, mu, BASIS_PHI, BASIS_Targets, sigma_pos, mu_pos, S_in, Q_in,$
	S_out, Q_out, factor, logml, gamm, betabasis_phi, beta_noise, SIGMA

	N = N_elements(BASIS[*,0])	
	
	U = PHI ## transpose(PHI) * beta + diag_matrix(Alpha)

; 	if (n_elements(Alpha) eq 12) then begin
; 		stop
; 	endif

	if (n_elements(U) eq 1) then begin
		U = sqrt(U)
	endif else begin

; Cholesky decomposition and zero out the lower triangular portion
		la_choldc, U, status=var, /upper
		U = zero_lower_triangular(U)
	endelse
	
	Ui = invert(U)
	SIGMA = Ui ## transpose(Ui)

; Posterior mean
	Mu_pos = reform((SIGMA ## (PHI ## Targets)) * beta)

; Data error and likelihood
	y = transpose(PHI) ## Mu_pos
	e = Targets - y
	ED = total(e*e)

	dataLikely = 0.5d0 * (N*alog(beta) - beta*ED)

; Compute log marginal likelihood
	logdetHOver2 = total(alog(diag_matrix(U)))
	logML = dataLikely - 0.5d0 * Mu_pos^2 ## transpose(Alpha) + 0.5d0 * total(alog(Alpha)) - logdetHOver2 - N/2.d0 * alog(2.d0*!DPI)

; Well-determinedness factors
	DiagC = total(Ui^2, 2)
	Gamm = 1.d0 - Alpha * DiagC

; Compute Q & S values
; Q: "quality" factor - related to how well the basis function contributes
; to reducing the error
; S: "sparsity factor" - related to how orthogonal a given basis function
; is to the currently used set of basis functions

	betaBASIS_PHI = beta * BASIS_PHI

	S_in = beta - total( (betaBASIS_PHI ## Ui)^2,1)	
	Q_in = beta * (BASIS_Targets - BASIS_PHI##Mu_pos)

	S_out = S_in
	Q_out = Q_in

	S_out[Used] = (Alpha * S_in[Used]) / (Alpha - S_in[Used])
	Q_out[Used] = (Alpha * Q_in[Used]) / (Alpha - S_in[Used])

; Relevance factor
	Factor = Q_out * Q_out - S_out

end