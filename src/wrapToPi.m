function lambda = wrapToPi(lambda)
%wrapToPi Wrap angle in radians to [-pi pi]
%
%   lambdaWrapped = wrapToPi(LAMBDA) wraps angles in LAMBDA, in radians,
%   to the interval [-pi pi] such that pi maps to pi and -pi maps to
%   -pi.  (In general, odd, positive multiples of pi map to pi and odd,
%   negative multiples of pi map to -pi.)
%
%   See also wrapTo2Pi, wrapTo180, wrapTo360.

% Copyright 2007 The MathWorks, Inc.
% $Revision: 1.1.6.2 $  $Date: 2007/08/20 16:36:02 $

lambda = wrapTo2Pi(lambda + pi) - pi;
