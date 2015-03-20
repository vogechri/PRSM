function [X, fX, i] = minimizeVektor(X, f, length, varargin)

% Minimize a differentiable multivariate function. With reparameterization.
%
% Usage: [X, fX, i] = minimize(X, f, length, P1, P2, P3, ... )
%
% where the starting point is given by "X" (D by 1), and the function named in
% the string "f", must return a function value and a vector of partial
% derivatives of f wrt X, the "length" gives the length of the run: if it is
% positive, it gives the maximum number of line searches, if negative its
% absolute gives the maximum allowed number of function evaluations. You can
% (optionally) give "length" a second component, which will indicate the
% reduction in function value to be expected in the first line-search (defaults
% to 1.0). The parameters P1, P2, P3, ... are passed on to the function f.
%
% The function returns when either its length is up, or if no further progress
% can be made (ie, we are at a (local) minimum, or so close that due to
% numerical problems, we cannot get any closer). NOTE: If the function
% terminates within a few iterations, it could be an indication that the
% function values and derivatives are not consistent (ie, there may be a bug in
% the implementation of your "f" function). The function returns the found
% solution "X", a vector of function values "fX" indicating the progress made
% and "i" the number of iterations (line searches or function evaluations,
% depending on the sign of "length") used.
%
% The Polack-Ribiere flavour of conjugate gradients is used to compute search
% directions, and a line search using quadratic and cubic polynomial
% approximations and the Wolfe-Powell stopping criteria is used together with
% the slope ratio method for guessing initial step sizes. Additionally a bunch
% of checks are made to make sure that exploration is taking place and that
% extrapolation will not be unboundedly large.
%
% See also: checkgrad 
%
% Copyright (C) 2001 - 2006 by Carl Edward Rasmussen (2006-09-08).

% this for per pixel optimization
%MAX = 5;                    % max 20 function evaluations per line search, since parallel: inner + outer: 10 each, 8 seems ok
%SIG = 0.051; RHO = SIG/2; % SIG and RHO are the constants controlling the Wolfe-

INT = 0.1;    % don't reevaluate within 0.1 of the limit of the current bracket
EXT = 3.0;                  % extrapolate maximum 3 times the current step-size
%MAX = 5;                    % max 20 function evaluations per line search, since parallel: inner + outer: 10 each, 8 seems ok
MAX = 20;                   % max 20 function evaluations per line search, since parallel: inner + outer: 10 each, 8 seems ok
RATIO = 10;                                       % maximum allowed slope ratio
%SIG = 0.051; RHO = SIG/2; % SIG and RHO are the constants controlling the Wolfe-
SIG = 0.1; RHO = SIG/2; % SIG and RHO are the constants controlling the Wolfe-
% Powell conditions. SIG is the maximum allowed absolute ratio between
% previous and new slopes (derivatives in the search direction), thus setting
% SIG to low (positive) values forces higher precision in the line-searches.
% RHO is the minimum allowed fraction of the expected (from the slope at the
% initial point in the linesearch). Constants must satisfy 0 < RHO < SIG < 1.
% Tuning of SIG (depending on the nature of the function to be optimized) may
% speed up the minimization; it is probably not worth playing much with RHO.

% The code falls naturally into 3 parts, after the initial line search is
% started in the direction of steepest descent. 1) we first enter a while loop
% which uses point 1 (p1) and (p2) to compute an extrapolation (p3), until we
% have extrapolated far enough (Wolfe-Powell conditions). 2) if necessary, we
% enter the second loop which takes p2, p3 and p4 chooses the subinterval
% containing a (local) minimum, and interpolates it, unil an acceptable point
% is found (Wolfe-Powell conditions). Note, that points are always maintained
% in order p0 <= p1 <= p2 < p3 < p4. 3) compute a new search direction using
% conjugate gradients (Polack-Ribiere flavour), or revert to steepest if there
% was a problem in the previous line-search. Return the best value so far, if
% two consecutive line-searches fail, or whenever we run out of function
% evaluations or line-searches. During extrapolation, the "f" function may fail
% either with an error or returning Nan or Inf, and minimize should handle this
% gracefully.

if max(size(length)) == 2, red=length(2); length=length(1); else red=1; end
if length>0, S='Linesearch'; else S='Function evaluation'; end 

[split, N] = size(X);
optimize3D = 1;   % 3D entities are optimized or 2d entities

ratio = RATIO * ones(1,N);

i = 0;                                            % zero the run length counter
ls_failed = 0;                                    % no previous line search has failed
[f0, df0]  = feval(f, X, varargin{:});             % get function value and gradient
fX = f0;
i = i + (length<0);                               % count epochs?!
s = -df0; d0 = -sum( s.^2, 1 );%diag(-s'*s);        % initial search direction (steepest) and slope
temp = (d0 == 0);
d0 = -max(-d0,0.0001);
df0(:,temp) = 0.0001;
x3 = red./(1-d0);                                  % initial step is red/(|s|+1) SHOULD BE HALF A PIXEL (or less 0.25 pixel ?)!! 
%x3 = red./sqrt(-d0);

%if any(any(isnan(df0)))
%  a =1;
%end
  %x3 =0.25;
% THIS WORKS GREAT FOR SMALL s:
%
% s*s ~ 0 -> x3 = red, in reallity to reach about red, x3 * s = red -> x3
% big lol

%if f0 < 0.0001
%  fX = f0;
%  return;
%end
%checkgrad(f, X, 0.0001, varargin{:})
%checkgrad(f, X, 0.000000000001, varargin{:});
%checkgradVector(f, X, 0.05, pcg0, varargin{:});
%feval(@EnergyInDataMinimize, X, varargin{:});


x4 = zeros(1,N); d4 = zeros(1,N); f4 = zeros(1,N);
tttemp = zeros(1,N);
ttemp4D = zeros(size(X));

while i < abs(length)                                      % while not finished
  i = i + (length>0);                                      % count iterations?!

  X0 = X; F0 = f0; dF0 = df0;                   % make a copy of current values
  if length>0, M = MAX; else M = min(MAX, -length-i); end
  allImproved = false(1,N);

  startx3 = x3;
  
  %  alreadyImproved = zeros(1,N);
  while 1                             % keep extrapolating as long as necessary
    x2 = zeros(1,N); f2 = f0; d2 = d0; f3 = f0; df3 = df0;
    success = 0;
    while ~success && M > 0
      try
        ttemp4D =  X + s .* repmat( x3, split, 1);% moved from below
        M = M - 1; i = i + (length<0);                         % count epochs?!
        [f3 df3] = feval(f, ttemp4D, varargin{:});
        allImproved = (isnan(f3) | isinf(f3) | any(isnan(df3)+isinf(df3)));
        if any(allImproved), error(''), end
        success = 1;
      catch                            % catch any error which occured in f
        ttemp = (x2+x3)/2;
        x3(allImproved) = ttemp(allImproved);   % bisect fails and try again
      end
    end
    improved = f3 < F0; % keep best values
%    ttemp4D =  X + pcg .* s * diag(x3); %moved up

    X0(:,improved)   = ttemp4D(:,improved);
    F0(improved)     = f3(improved); 
    dF0(:,improved)  = df3(:,improved);

%    alreadyImproved = improved | alreadyImproved;
     % by keeping x3 constant when allImproved is true d3 can only vary
     % when allImproved is false for that pixel
     d3 = sum( df3.*s, 1 );%diag(df3'*s);                   % new slope
    % note that the already improved items do not get updated but remain
    % constant, that saves conditionally updating x1,x2, d3, etc.
    
    % when d3 is positive, the derivative has switched sign. So we have
    % moved too far and can now be sure to find the minimum along the line
    % f(X+x3*s) in between 0 and x3.
    allImproved = (d3 > SIG*d0) | (f3 > f0+(x3.*d0*RHO));
    if all(allImproved) || M == 0 % are we done extrapolating?
      break
    end

    % this can go through completely (already for solved directions) BUT x3
    % must only be altered for unsolved directions, save because d1,f1,x1
    % are not used below and x2,f2, etc. get overwritten in this loop again
    x1 = x2; f1 = f2; d1 = d2;                        % move point 2 to point 1
    x2 = x3; f2 = f3; d2 = d3;                        % move point 3 to point 2
    A = 6*(f1-f2)+3*(d2+d1).*(x2-x1);                 % make cubic extrapolation
    B = 3*(f2-f1)-(2*d1+d2).*(x2-x1);
    ttemp = B.^2-A.*d1.*(x2-x1);
    temp  = x1-d1.*(x2-x1).^2./(B+sqrt(max(ttemp,0)));

    x3(~allImproved) = temp(~allImproved) ; % num. error possible, ok!
    temp = ((ttemp < 0) | isnan(x3) | isinf(x3) | (x3 < 0));

    % ~allImproved geschnitten temp = temp
    if any(temp & ~allImproved) % num prob | wrong sign?
      x3(temp & ~allImproved) = x2(temp & ~allImproved)*EXT; % extrapolate maximum amount
    end

    temp = (x3 > x2*EXT) & ~temp;
%    if any(temp && ~alreadyImproved) % new point beyond extrapolation limit?
    x3(temp & ~allImproved) = x2(temp & ~allImproved)*EXT;         % extrapolate maximum amount

    ttemp = x2+INT*(x2-x1);
    temp = (x3 < ttemp) & ~temp;         % new point too close to previous point?
    x3(temp & ~allImproved) = ttemp(temp & ~allImproved);
  end                                                       % end extrapolation

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  M = max(M,MAX);
  allImproved = (abs(d3) > -SIG*d0) | (f3 > f0+(x3.*d0*RHO));
  while any(allImproved) && M > 0  % keep interpolating
    temp = (d3 > 0) | (f3 > f0+(x3.*d0*RHO));
    x4(temp) = x3(temp); f4(temp) = f3(temp); d4(temp) = d3(temp);       % move point 3 to point 4
    x2(~temp) = x3(~temp); f2(~temp) = f3(~temp); d2(~temp) = d3(~temp); % move point 3 to point 2

    % alter x3 only if loop is entered, that is when ~allImproved
    temp = f4 > f0;
    if any(temp)
      ttemp = x2-(0.5*d2.*(x4-x2).^2)./((f4-f2)-d2.*(x4-x2));
      x3(temp & allImproved) = ttemp(temp & allImproved);% quadratic interpolation
    end
    if any(~temp)
      A = 6*(f2-f4)./(x4-x2)+3.*(d4+d2);         % cubic interpolation
      B = 3*(f4-f2)-(2*d2+d4).*(x4-x2);
      tttemp = B.^2-A.*d2.*(x4-x2).^2;     
      ttemp = x2+(sqrt(tttemp)-B)./A;
      x3(~temp & allImproved) = ttemp (~temp & allImproved); % num. error possible, ok!
    end
    
    temp = (tttemp < 0) | isnan(x3) | isinf(x3);
    if any(temp)
      x3(temp) = (x2(temp)+x4(temp))/2;               % if we had a numerical problem then bisect
    end


    % probably only active for alterd x3 anyway but still check ~allImproved
    ttemp = max(min(x3, x4-INT*(x4-x2)),x2+INT*(x4-x2));
    x3(allImproved) = ttemp(allImproved);  % don't accept too close
    ttemp4D = X + s .* repmat( x3, split, 1);%diag(x3);
    [f3 df3] = feval(f, ttemp4D, varargin{:});
    
    % with update here too ?SEE BELOW
%    improved = (f3 < F0) & ~allImproved; % keep best values and keep the
%    values of the improved pixels by the extrapolation part
    improved = (f3 < F0); % keep best values
    X0(:,improved) = ttemp4D(:,improved);
    F0(improved) = f3(improved);
    dF0(:,improved) = df3(:,improved);
    %    if f3 < F0, X0 = X+x3*s; F0 = f3; dF0 = df3; pcg0 = pcgx;end % keep best values

    M = M - 1; i = i + (length<0);                             % count epochs?!
    ttemp = sum(df3.*s,1);%diag(df3'*s); % MAYBE A LITTLE LARGE A MATRIX : N^4 only N^2 
    % needed here that is N = 128 pixel, double : 8 byte 2GB :( 64 bit
    % worse :-O
    
    % since d3 is used below, only update if improved when ~allImproved is
    % true. So F0, etc correspond to d3 here. update for the others to
    % steer the next iteration.

    % with update here too ?
    % all need a new d3, especially those which did improve and those for
    % which these loop is about [~allImproved]
    d3( allImproved ) = ttemp( allImproved ); % new slope: yields direction of the update
    % here we update the improved pixels in general - disregarding whether the
    % extrapolation already worked for them. SEE ABOVE
    
    % observe line:
    %ttemp = s*diag(diag(df3'*df3-df0'*df3)./diag(df0'*df0)) - df3;
    % here it is assumed that df3 contains the correct derivative for 
    % for the new point X (df0) [why not use df0 here?] this also
    % effectively states that df3 must be corresponding to d3 or allImproved
    % respectively
%    d3(improved | ~allImproved) = ttemp(improved | ~allImproved); % new slope: yields direction of the update
    % x3 and therefore d3 remain constant for pixels with allImproved true
    % therefore allImproved remains true for those pixels
    allImproved = (abs(d3) > -SIG*d0) | (f3 > f0+(x3.*d0*RHO));
%    disp([x3 d3 d0 f3' f0' F0' allImproved']);
  end                                                       % end interpolation
%d3' when it is 0 the minimum is found
%idSucc = (abs(d3) < 0.01);
%sumSucc = sum(idSucc);
%fprintf('Iteration %d: %d\n', i, sumSucc);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % must just invert, that means the following row is unnecessary 
  % unless the former if clause was entered:
%  allImproved = (abs(d3) < -SIG*d0)' & (f3 < f0+(x3.*d0*RHO)'); 
%  if abs(d3) < -SIG*d0 && f3 < f0+x3*RHO*d0          % if line search succeeded
  if any(~allImproved)
    ttemp4D = X + s .* repmat( x3, split, 1);%diag(x3);
    X(:,~allImproved) = ttemp4D(:,~allImproved); 
    f0(~allImproved) = f3(~allImproved);
    fX = [fX' f0']';                     % update variables
%    fprintf('%s %6i;  Value %4.6e\r', S, i, f0);
%    ttemp4D = s*diag(sum(df3.*df3-df0.*df3)./sum(df0.*df0)) - df3;

    ttemp = sum(df0.*df0);
    ttemp (ttemp == 0) = 0.00001;
    ttemp4D = s.* repmat( sum(df3.*df3-df0.*df3)./ttemp, split, 1) - df3;
    %    ttemp4D = s*diag(diag(df3'*df3-df0'*df3)./diag(df0'*df0)) - df3;
    s(:,~allImproved) = ttemp4D(:,~allImproved);   % Polack-Ribiere CG direction

    df0(:,~allImproved) = df3(:,~allImproved);     % swap derivatives
    d3(~allImproved) = d0(~allImproved); 
    ttemp = sum(df0.*s,1);% diag(df0'*s);
    d0(~allImproved) = ttemp(~allImproved);
    temp = d0 > 0; %(if d0 > 0)  % new slope must be negative
    s(:,temp & ~allImproved) = -df0(:,temp  & ~allImproved); 

    ttemp = -sum(s.*s, 1);%-diag(s'*s);
    d0(temp & ~allImproved) = ttemp(temp & ~allImproved); % otherwise use steepest direction
%  end
    ttemp =  x3 .* min(ratio, d3./(d0-realmin));
    x3(~allImproved) =ttemp(~allImproved);      % slope ratio but max RATIO

%    temp = (d0 == 0);
%    d0 = -max(-d0,0.0001);
%    df0(:,temp) = 0.0001;

%  if any(any(isnan(s)))
%    a =1;
%  end

    ls_failed = 0;                              % this line search did not fail
  end
  if any(allImproved)
    X(:,allImproved)   = X0(:,allImproved);
    f0(allImproved)    = F0(allImproved); 
    df0(:,allImproved) = dF0(:,allImproved);    % restore best point so far

    if ls_failed || i > abs(length)      % line search failed twice in a row
      break;                             % or we ran out of time, so we give up
    end

    s(:,allImproved)  = -df0(:,allImproved); 

    ttemp = -sum(s.*s,1);%diag(-s'*s);
    d0(allImproved) = ttemp(allImproved);       % try steepest

    temp = (d0 == 0);
    d0 = -max(-d0,0.0001);

    ttemp = 1./(1-d0);%red./sqrt(-d0+0.0001);
    %ttemp = 1./sqrt(-d0+0.0001);
    x3(allImproved) = ttemp(allImproved);
    if all(allImproved)
      ls_failed = 1;                            % this line search failed
    end
  end

end
%fprintf('\n');