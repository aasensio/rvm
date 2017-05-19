; Support function to compute basis
function distSquared, X, Y
	nx = n_elements(X)
	ny = n_elements(Y)
	
	return, replicate(1.d0,ny) # X^2 + Y^2 # replicate(1.d0,nx) - 2.d0 * X#Y
end

; Support function to compute basis
function distSquared, X, Y
	nx = n_elements(X)
	ny = n_elements(Y)

	res = fltarr(nx,ny)
	for i = 0, ny-1 do begin
		res[*,i] = (x-y[i])^2
	endfor
end