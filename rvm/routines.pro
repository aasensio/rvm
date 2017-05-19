function delete_column_array, X, which
	n = n_elements(X[*,0])
	ind = replicate(1,n)
	ind[which] = 0
	keeprow = where(ind eq 1,count)
	if (count ne -1) then begin
		out = X[keeprow,*]
	endif else begin
		out = X
	endelse
	return, out 
end

function delete_row_array, X, which
	n = n_elements(X[0,*])
	ind = replicate(1,n)
	ind[which] = 0
	keeprow = where(ind eq 1,count)
	if (count ne -1) then begin
		out = X[*,keeprow]
	endif else begin
		out = X
	endelse
	return, out
end

; Support function to compute basis
; function distSquared, X, Y
; 	nx = n_elements(X)
; 	ny = n_elements(Y)
; 
; 	return, replicate(1.d0,ny) # X^2 + Y^2 # replicate(1.d0,nx) - 2.d0 * X#Y
; end

function distSquared, X, Y
	nx = n_elements(X)
	ny = n_elements(Y)

	res = fltarr(nx,ny)
	for i = 0, ny-1 do begin
		res[*,i] = (x-y[i])^2
	endfor
	return, res
end

; Zero out the lower triangular of a matrix
function zero_lower_triangular, X
	n = n_elements(X[*,0])
	for i = 1, n-1 do X[0:i-1,i] = 0.d0
	return, X
end

; Append a new element on an array if it exists. If not, create the array
pro push, array, append
	if n_elements(array) eq 0 then array=[append] else array=[array,append]
end