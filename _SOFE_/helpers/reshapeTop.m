function R = reshapeTop(slots, varargin)
	% Creates matrix with column i containing the next slots(i)
	% entries from varargin
	% Columnwise, the matrix is filled from top to bottom,
	% unspecified entries being filled with zero.
	rows = max(slots);
	cols = numel(slots);
	R = repmat((1:rows)',1,cols); 
	R = bsxfun(@le,R,slots(:)');
	[ii,jj] = find(R>0);
	if nargin > 1
		if sum(slots) ~= numel(varargin{1})
		  error('! Incompatible data !');
		end
		R = full(sparse(ii, jj, varargin{1}, rows, cols));
	else
		R = full(sparse(ii, jj, ii, rows, cols));
	end
end
