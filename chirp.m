function   [x, nn] = chirp( T, W, p )
%CHIRP    generate a sampled chirp signal
%-----       exp(j(W/T)pi*t^2)   -T/2 <= t < +T/2
%
%   Usage:   X = chirp( T, W, <P> )
%
%      X :  N=pTW samples of a "chirp" signal
%      T :  time duration from -T/2 to +T/2
%      W :  swept bandwidth from -W/2 to +W/2
%
%   optional (default is P = 1)
%      P :  samples at P times the Nyquist rate (W)
%            i.e., sampling interval is 1/(PW)

%---------------------------------------------------------------
% copyright 1994, by C.S. Burrus, J.H. McClellan, A.V. Oppenheim,
% T.W. Parks, R.W. Schafer, & H.W. Schussler.  For use with the book
% "Computer-Based Exercises for Signal Processing Using MATLAB"
% (Prentice-Hall, 1994).
%---------------------------------------------------------------

if nargin < 3
   p = 1;   end    
J = sqrt(-1);
%--------------
delta_t = 1/(p*W);
N = round( p*T*W );   %--- same as T/delta_t
nn = [0:N-1]';
x = exp( J*pi*W/T * (delta_t*nn - T/2).^2 );
%plot(nn, x)