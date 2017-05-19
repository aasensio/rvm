pro SparseBayes, BASIS_IN, Targets, options, settings, parameter, hyperparameter, diagnostic, SIGMA

	BASIS = BASIS_IN
	sb_initialization, BASIS, BasisScales, settings, Targets, Alpha, beta, mu, PHI, Used
	
; Cache some quantities
 	BASIS_PHI = BASIS ## PHI
 	BASIS_Targets = BASIS ## Targets

; Full computation
;
; Initialise with a full explicit computation of the statistics
;
; NOTE: The AISTATS paper uses "S/Q" (upper case) to denote the key
; "sparsity/quality" Factors for "included" basis functions, and "s/q"
; (lower case) for the factors calculated when the relevant basis
; functions are "excluded".
;
; Here, for greater clarity:
;
;	All S/Q are denoted by vectors S_in, Q_in
;	All s/q are denoted by vectors S_out, Q_out

	sb_fullstatistics, BASIS, PHI, Targets, Used, Alpha, beta, mu, BASIS_PHI, BASIS_Targets, sigma_pos, mu, S_in, Q_in,$
		S_out, Q_out, factor, logml, gamm, BASIS_B_PHI, beta_noise, SIGMA
	

	N = n_elements(BASIS[*,0])
	M_full = n_elements(BASIS[0,*])
	M = n_elements(PHI[0,*])

	addCount = 0
	deleteCount = 0
	updateCount = 0

	maxLogLike = -1.d10

	if (n_tags(OPTIONS) eq 0) then begin
		OPTIONS = create_struct('monitor', 1)
	endif else begin
		OPTIONS = create_struct('monitor', 1, OPTIONS)
	endelse

	CONTROLS = create_struct('BetaUpdateStart', 10, 'BetaUpdateFrequency', 5, 'ZeroFactor', 1.d-12, $
		'PriorityAddition', 0, 'PriorityDeletion', 1, 'AlignmentMax', 0.999d0, 'MinDeltaLogAlpha', 1.d-3)

	logMarginalLog = [logMl]
	count_loop = 0

; Keep lists of functions that are near identical
; 	Aligned_out = NULL
; 	Aligned_in = NULL
	alignDeferCount = 0

	loop = 0
	LAST_ITERATION = 0

	ACTION_REESTIMATE	= 0
	ACTION_ADD			= 1
	ACTION_DELETE		= -1
	ACTION_TERMINATE		= 10
	ACTION_NOISE_ONLY		= 11
	ACTION_ALIGNMENT_SKIP	= 12

; 	freeBasis = NULL

	while (not LAST_ITERATION) do begin
		loop = loop + 1

;*************************
; DECISION PHASE
;*************************

; Compute change in likelihood
		DeltaML = fltarr(M_full)
		Action = replicate(ACTION_REESTIMATE, M_full)
		UsedFactor = Factor[Used]

; Re-estimation
		iu = where(UsedFactor gt CONTROLS.ZeroFactor, count)
		if (count ne 0) then begin
			index = Used[iu]
			NewAlpha = S_out[index]^2 / Factor[index]
			Delta = 1.d0 / NewAlpha - 1.d0 / Alpha[iu]

; Quick computation of change in log-likelidoog given all re-estimations
			DeltaML[index] = 0.5d0 * (Delta * Q_in[index]^2 / (Delta * S_in[index] + 1.d0) - alog(1.d0+S_in[index]*Delta))
		endif


; Deletion
		iu = where(UsedFactor le CONTROLS.ZeroFactor, count)
		if (count ne 0) then begin
			index = Used[iu]
		endif
		
; Any to delete
		anyToDelete = count ne 0 and M gt 1
		if (anyToDelete) then begin
			DeltaML[index] = -0.5d0 * (Q_out[index]^2 / (S_out[index]+Alpha[iu]) - alog(1.d0+S_out[index] / Alpha[iu]))
			Action[index] = ACTION_DELETE			
		endif

; Addition
		index = where(Factor gt CONTROLS.ZeroFactor, count)
		GoodFactor = Factor * 0
		if (count ne 0) then begin
			GoodFactor[index] = 1
		endif
		GoodFactor[Used] = 0
		if (n_elements(Aligned_out) ne 0) then GoodFactor[Aligned_out] = 0
		index = where(GoodFactor gt CONTROLS.ZeroFactor, anyToAdd)
		
		if (anyToAdd ne 0) then begin

; Quick computation of change in log-likelidoog
			quot = Q_in[index]^2 / S_in[index]
			DeltaML[index] = 0.5d0 * (quot - 1.d0 - alog(quot))			
			Action[index] = ACTION_ADD
		endif
		
		if ( (anyToAdd and CONTROLS.PriorityAddition) or (anyToDelete and CONTROLS.PriorityDeletion) ) then begin

; We won't perform re-estimation this iteration
			ind = where(Action eq ACTION_REESTIMATE, count)
			if (count ne 0) then DeltaML[ind] = 0

; We should enforce add if preferred and delete
			if (anyToAdd and CONTROLS.PriorityAddition and not CONTROLS.PriorityDeletion) then begin
				ind = where(Action eq ACTION_DELETE, count)
				if (count ne 0) then DeltaML[ind] = 0
			endif
			if (anyToDelete and CONTROLS.PriorityDeletion and not CONTROLS.PriorityAddition) then begin
				ind = where(Action eq ACTION_ADD, count)
				if (count ne 0) then DeltaML[ind] = 0				
			endif			
		endif

; Finally, choose the action that results in the greatest change in likelihood
		deltaLogMarginal = max(DeltaML, nu)
		selectedAction = Action[nu]

		anyWorthwhileAction = deltaLogMarginal gt 0

; If basis nu is already in the model, find its index
		if (selectedAction eq ACTION_REESTIMATE or selectedAction eq ACTION_DELETE) then begin
			j = where(Used eq nu)
			j = j[0]
		endif

; Get the individual basis vector for update and compute its optimal alpha
		PPhi = BASIS[*,nu]		
		newAlpha = S_out[nu]^2 / Factor[nu]

; Termination conditions				
		if (not anyWorthwhileAction or $
			(selectedAction eq ACTION_REESTIMATE and not anyToDelete)) then begin
			if (n_elements(j) eq 0) then begin
				change = 0.d0			
			endif else begin
				change = abs(alog(newAlpha) - alog(Alpha[j]))
; 				print, 'change ', Factor
			endelse
			if (change lt CONTROLS.MinDeltaLogAlpha) then begin
				selectedAction = ACTION_TERMINATE			
			endif
		endif

; Alignment checks for addition
		if (selectedAction eq ACTION_ADD) then begin

; Basic test for correlated basis vectors
			p = Pphi ## transpose(Phi)			
			findAligned = where(p gt CONTROLS.AlignmentMax, numAligned)			
			if (numAligned gt 0) then begin
; The added basis function is effectively indistinguishable from one present already
				selectedAction = ACTION_ALIGNMENT_SKIP
				alignDeferCount = alignDeferCount + 1
; Take note not to try this again
; 				Aligned_out = [Aligned_out, replicate(nu, numAligned)]
				push, Aligned_out, replicate(nu, numAligned)
				
; 				Aligned_in = [Aligned_in, Used[findAligned]]
				push, Aligned_in, Used[findAligned]
			endif
		endif

; Alignment checks for deletion
		if (selectedAction eq ACTION_DELETE) then begin
			numAligned = 0
			if (n_elements(Aligned_in) ne 0) then begin
				findAligned = where(Aligned_in eq nu, numAligned)			
			endif
			if (numAligned gt 0) then begin
				reinstated = Aligned_out[findAligned]
; 				Aligned_in = delete_column_array(Aligned_in, findAligned)
; 				Aligned_out = delete_column_array(Aligned_out, findAligned)
				print, 'Alignment reinstatement of ', reinstated
			endif
		endif

;*************************
; ACTION PHASE
;*************************
		UPDATE_REQUIRED = 0

		switch(selectedAction) of
			ACTION_REESTIMATE: begin
			
; Basis function nu is already in the model and we're reeestimatig its alpha				
				oldAlpha = Alpha[j]
				Alpha[j] = newAlpha
				s_j = SIGMA[*,j]
				deltaInv = 1.d0 / (NewAlpha - oldAlpha)
				kappa = 1.d0 / (SIGMA[j,j] + deltaInv)
				tmp = kappa * s_j
				SIGMANEW = SIGMA - tmp ## s_j
				deltaMu = -Mu[j]*tmp
				Mu = Mu + deltaMu

				S_in = S_in + kappa*(BASIS_B_PHI ## s_j)^2
				Q_in = Q_in - BASIS_B_PHI ## deltaMu
				
				updateCount = updateCount + 1
				UPDATE_REQUIRED = 1				
				break
			end

			ACTION_ADD: begin
				BASIS_Pphi = BASIS ## PPhi
				BASIS_PHI = [BASIS_PHI, BASIS_Pphi]
				B_Pphi = beta*Pphi
				BASIS_B_PpHI = beta * BASIS_Pphi

				tmp = (B_Pphi ## transpose(PHI)) ## SIGMA

				Alpha = [Alpha, newAlpha]
				
				PHI = [[PHI], [PPhi]]

				s_ii = 1.d0 / (newAlpha + S_in[nu])
				s_i = -s_ii * tmp

				TAU = reform(-s_i ## tmp)
				SIGMANEW = [[SIGMA+TAU,transpose(s_i)],[s_i,s_ii]]
				mu_i = s_ii * Q_in[nu]
				deltaMu = [-mu_i*tmp, mu_i]
				Mu = [Mu, 0] + deltaMu

				mCi = BASIS_B_PpHI - BASIS_B_PHI ## tmp
				S_in = S_in - s_ii * mCi^2
				Q_in = Q_in - mu_i * mCi

				Used = [Used, nu]
				addCount = addCount + 1

				UPDATE_REQUIRED = 1
				break
			end

			ACTION_DELETE: begin				
				BASIS_PHI = delete_column_array(BASIS_PHI, j)
				PHI = delete_row_array(PHI, j)
				Alpha = delete_column_array(Alpha, j)

				s_jj = SIGMA[j,j]
				s_j = SIGMA[*,j]
				tmp = s_j / s_jj[0]
				SIGMANEW = SIGMA- tmp ## s_j

				SIGMANEW = delete_column_array(SIGMANEW, j)
				SIGMANEW = delete_row_array(SIGMANEW, j)

				deltaMu = -Mu[j]*tmp
				mu_j = Mu[j]
				Mu = Mu + deltaMu
				
				Mu = delete_column_array(Mu, j)

				jPm = BASIS_B_PHI ## s_j
				S_in = S_in + jPm^2 / s_jj
				Q_in = Q_in + jPm * mu_j / s_jj

				Used = delete_column_array(Used, j)
				deleteCount = deleteCount + 1

				UPDATE_REQUIRED = 1				
				break
			end
		endswitch

		M = n_elements(Used)

		if (UPDATE_REQUIRED) then begin
			S_out = S_in
			Q_out = Q_in
			
			tmp = Alpha / (Alpha - S_in[Used])
			S_out[Used] = tmp * S_in[Used]
			Q_out[Used] = tmp * Q_in[Used]

			Factor = (Q_out * Q_out - S_out)
			SIGMA = SIGMANEW
			Gamm = 1.d0 - Alpha * diag_matrix(SIGMA)
			BASIS_B_PHI = beta * BASIS_PHI

			logMl = logMl + deltaLogMarginal
			count_loop = count_loop + 1
			logMarginalLog = [logMarginalLog, logML]

		endif

; Something went wrong. Recompute statistics
		if (total(Gamm) lt 0) then begin
			print, 'Recomputing full statistics'
			sb_fullstatistics, BASIS, PHI, Targets, Used, Alpha, beta, mu, BASIS_PHI, BASIS_Targets, sigma_pos, mu, S_in, Q_in,$
					S_out, Q_out, factor, logml, gamm, BASIS_B_PHI, beta_noise, SIGMA
		endif
			
; Recompute statistics if needed
		if (not OPTIONS.fixedNoise and (selectedAction eq ACTION_TERMINATE or loop le CONTROLS.BetaUpdateStart or $
			loop mod CONTROLS.BetaUpdateFrequency eq 0)) then begin
; Gaussian noise estimate			
			betaZ1 = beta
			y = transpose(PHI) ## Mu
			e = (Targets - y)
			beta = (N- total(Gamm)) / total(e^2)
; Work-around zero-noise issue
			beta = min([beta, 1.d6 / variance(Targets)])
			deltaLogBeta = alog(beta) - alog(betaZ1)

; Full re-computation of statistics after beta update
			if (abs(deltaLogBeta) gt 1.d-6) then begin
				sb_fullstatistics, BASIS, PHI, Targets, Used, Alpha, beta, mu, BASIS_PHI, BASIS_Targets, sigma_pos, mu, S_in, Q_in,$
					S_out, Q_out, factor, logml, gamm, BASIS_B_PHI, beta_noise, SIGMA

				count_loop = count_loop + 1
				logMarginalLog = [logMarginalLog, logML]

				if (selectedAction eq ACTION_TERMINATE) then begin
					selectedAction = ACTION_NOISE_ONLY
					print, 'Noise update. Termination deferred'
				endif
			endif
		endif


; End of cycle processing
		if (selectedAction eq ACTION_TERMINATE) then begin
			print, format='(A,I3,A,E12.4,A)', 'Stopping at iteration ', loop, ' (Max_delta_ml=', deltaLogMarginal,' )'
			print, format='(A,E12.4,A,F8.3,A,I3,A,F8.3)','L=', logMl[0]/N, ' - Gamma=', total(Gamm), ' (M=', M,') -  s=', sqrt(1.d0/beta)
			break
		endif

; Check for natural termination
		ITERATION_LIMIT = loop eq OPTIONS.iterations
		LAST_ITERATION = ITERATION_LIMIT		

; Print some info
 		if (loop / OPTIONS.monitor eq (1.d0*loop) / OPTIONS.monitor) then begin
			print, format='(I4,A,E12.4,A,F8.3,A,I3,A,F8.3)', loop, 'L=', logMl[0]/N, ' Gamma=', total(Gamm), ' (M=', M,')   s=', sqrt(1.d0/beta)
		endif

	endwhile

; Final computation of the statistics
	print, 'Final computation of statistics'
	sb_fullstatistics, BASIS, PHI, Targets, Used, Alpha, beta, mu, BASIS_PHI, BASIS_Targets, sigma_pos, mu, S_in, Q_in,$
		S_out, Q_out, factor, logml, gamm, BASIS_B_PHI, beta_noise, SIGMA
		
; Output values
	index = sort(Used)
	
	parameter = create_struct('Relevant', Used[index], 'Value', Mu[index] / BasisScales[Used[index]])
 	hyperparameter = create_struct('Alpha', Alpha[index] / BasisScales[Used[index]]^2, 'beta', beta)
 	diagnostic = create_struct('Gamma', Gamm[index], 'logMarginal', logMarginalLog[0:count_loop], 'iterations', loop, $
		'S_Factor', S_out, 'Q_Factor', Q_out, 'M_full', M_full)

; Summary
	if (OPTIONS.diagnosticLevel gt 1) then begin
		print, 'Action summary'
		print, '=============='
		tot = addCount + deleteCount + updateCount + alignDeferCount
		print, format='(A,I3,A,I3,A)', 'Added       -> ', addCount, ' (', addCount/(tot*1.d0)*100.d0, '%)'
		print, format='(A,I3,A,I3,A)', 'Deleted     -> ', deleteCount, ' (', deleteCount/(tot*1.d0)*100.d0, '%)'
		print, format='(A,I3,A,I3,A)', 'Reestimated -> ', updateCount, ' (', updateCount/(tot*1.d0)*100.d0, '%)'
		print, format='(A,I3,A)', 'Total of ', count, ' likelihood updates'
	endif

end