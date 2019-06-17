tmp = pwd;
cd([SOFE.getCorePath() '/helpers/stenglib/']);
%
if ~exist('OCTAVE_VERSION', 'builtin') % Matlab
  mex CFLAGS='$CFLAGS -std=c99 -o2 -march=native' tprod.c -lmwblas
  mex CFLAGS='$CFLAGS -std=c99 -o2 -march=native' fsparse.c -lmwblas
else % Octave
  mex -std=c99 tprod.c
  mex -std=c99 fsparse.c
%  delete('fsparse.o');
%  delete('tprod.o');
end
%
cd(tmp);