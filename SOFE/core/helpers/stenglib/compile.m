tmp = pwd;
cd([SOFE.getCorePath() '/helpers/stenglib/']);
%
if ~exist('OCTAVE_VERSION', 'builtin') % Matlab
  mex tprod.c -lmwblas
  mex fsparse.c -lmwblas
else % Octave
  mex -std=c99 tprod.c
  mex -std=c99 fsparse.c
%  delete('fsparse.o');
%  delete('tprod.o');
end
%
cd(tmp);