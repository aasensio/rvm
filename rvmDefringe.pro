@rvm/routines
@rvm/sparsebayes
@rvm/sb_initialization
@rvm/sb_preprocessbasis
@rvm/sb_fullstatistics

pro rvmDefringe

  noise = 1e-3

; Generate mock spectrum with fringes
  lambda = findgen(500)
  line = 1.0 - 0.4*voigt(0.1, (lambda - 300)/20.0) - 0.5 * voigt(0.1, (lambda - 100)/10.0)
  fringes = 0.01 * sin(5*lambda)

  spectrum = line + fringes + noise * randomn(seed,500)

; Number of observed points 
  N = n_elements(spectrum)

; X is the wavelength in our case
  X = lambda

; Measured intensity. For some stupid reasons that I have not fixed,
; it has to be in dimensions [1,Ndata]
  Outputs = transpose(spectrum)
  
  
; The basis will consist of Gaussians centered at the observed points
; with several widths and sines and cosines of different frequencies
  left = [80, 250]
  right = [130, 380]
  NGauss = total(right-left)


  NFreqs = 201
  NWidths = 50
  NBasis = 2 * NFreqs + NWidths * NGauss + 1
  Basis = dblarr(NBasis,N)

  omegas = findgen(NFreqs) / (NFreqs-1.d0) * (6.d0 - 4.d0) + 4.d0
  widths = findgen(NWidths) / (NWidths-1.d0) * (30.d0 - 5.d0) + 5.d0


; The fringes
  loop = 0
  for i = 0, NFreqs-1 do begin    
    Basis[loop,*] = sin(omegas[i] * lambda)
    loop = loop + 1    
  endfor
;  for i = 0, NFreqs-1 do begin
;    Basis[loop,*] = cos(omegas[i] * lambda)
;    loop = loop + 1    
;  endfor

; Now the lines
  for i = 0, 1 do begin
    for k = 0, NWidths-1 do begin
      for j = left[i], right[i] do begin
        Basis[loop,*] = voigt(0.0, (lambda - lambda[j])/widths[k])
        loop = loop + 1
      endfor
    endfor
  endfor

; Constant bias
  Basis[loop,*] = 1.d0

  
  
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

  ; Plot the original spectrum together with the inferred fringes and the substraction
  cwpal
  plot, spectrum
  oplot, y_pred, col=2

  w_infer_fringes = fltarr(M)
  w_infer_fringes[parameter.relevant] = parameter.value
  ind = where(parameter.relevant gt 2*NFreqs)
  w_infer_fringes[parameter.relevant[ind]] = 0.d0
  y_pred_fringes = transpose(BASIS) ## w_infer_fringes

  oplot, y_pred_fringes+1, col=3
  oplot, spectrum - y_pred_fringes, col=4
  
  stop  
  
end