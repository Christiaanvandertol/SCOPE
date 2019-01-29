function [b, err2, fcounter] = fixedp_brent_ari(func, x0, corner, tolFn, verbose)
% Find a fixed point of func(x) using Brent's method, as described by Brent 1971
%   func is a single-argument function, f(x) that returns a value the same size as x:
%        The goal is to find f(x) = x (or for Brent, f(x) - x = 0).
%   x0 is the initial guess (or 2 x n matrix if we want to generalize)
%   tol is the tolerance in x (or if two-valued, x, f(x)? )
%   corner (optional) is a known "edge" in the function that could slow down the algorithm
%      if specified and the first two points include the corner, the corner will be substituted as a starting point.
%
%  Written by: Ari Kornfeld, 2016-10

tolx_1 = 0; % we don't want to converge in x
subsetLimit = 0; % never use subsets of x
accelBisection = false; % should we divide the interval by ever-increasing fractions (1/2, 1/3, 1/4) when consecutive calls make no improvement?
accelLimit = 3; % but reset if more than accelLimit bisections don't help.
rotatePrev = true; % do we use Brent's original method of wiping out "xprev" in certain instances (forcing secant) or alway preserve a third value (allowing more inverse-quad).
iter_limit = 100;

if nargin < 5
    verbose = false;
end

if nargin < 4
    tolFn = eps; % default to FP precision: 2^-52 = 2.2204e-16
end
if nargin < 3
    corner = []; % default, no corners
end
track_fcount = (nargout > 2) || verbose;

recompute_b = false; % ensure that we return after a call on func(s); this allows us to skip the re-call in the main body (since func sets globals)

%  We keep track of the best value of x (b), the past two iterations (c, d),
%    and one contrapoint (a), i.e. a point that stradles the fixed point (i.e. err1 has the opposite sign of err2)
%  err1, err2, ... are the corresponding "y" values for y = f(x) - x (i.e. distance from the fixed point)
% We start by finding bracketing points a, b
a = x0;
[err1, b] = func(a);  % guess the second point by computing x_1 = f(x_0)

% special case: func may return a vector even though 'a' was scalar (because other variables not visible here were nonscalar).
% If so, expand a:
if length(a) == 1 && length(b) > 1  % anything else is an error
    a = repmat(a, size(b));  % or a * ones(size(b))?
end

err2 = func(b);
err2(isnan(err2)) = 0;  % if isnan, we're done, i.e. we shouldn't try
if track_fcount
    % count each item separately: count only the "necessary" iterations
    fcounter = 2*ones(size(b));
end

err_outside_tol = abs(err2) > tolFn;
if ~any(err_outside_tol)
    return  % we're done!
end
% ELSE
recompute_b = true;  % we'll be messing with it

% Now confirm that the two first guesses bracket zero.
%  NOTE: the algorithm may still succeed w/o bracketting, though it's not guaranteed.
not_bracketting_zero = (sign(err1) == sign(err2)) & err_outside_tol;
if any(not_bracketting_zero)
    %warning( 'Not all initial guesses bracket zero. Will fix it now.' );
    
    % first try a simple secant extrapolation
    x1 = b - err2.*(b - a)./(err2 - err1);
    err_x1 = func(x1);
    if track_fcount
        fcounter = fcounter + not_bracketting_zero; % count only the ones that needed this fcall
    end
    
    % since sign(err1) == sign(err2), compare the new value to either of those
    use_x1 = (sign(err_x1) ~= sign(err1)) & not_bracketting_zero;
    if any(use_x1)
        % save the better of the two original points into 'a'
        swap_to_a = (abs(err2) < abs(err1)  & use_x1);
        a(swap_to_a) = b(swap_to_a);      err1(swap_to_a) = err2(swap_to_a);
        % then put the new contrapoint into 'b'
        b(use_x1) = x1(use_x1);           err2(use_x1) = err_x1(use_x1);
    end
    
    % recompute a_too_high and iterate if necessary
    err_outside_tol = min(abs(err1), abs(err2)) > tolFn;
    not_bracketting_zero = (sign(err1) == sign(err2)) & err_outside_tol;
   % make 'a' the lower value, to make the rest simpler
   if any(not_bracketting_zero)
       swap_to_a = (err2 < err1  & not_bracketting_zero);
       [a(swap_to_a), b(swap_to_a)] = deal(b(swap_to_a), a(swap_to_a));
       [err1(swap_to_a), err2(swap_to_a)] = deal( err2(swap_to_a), err1(swap_to_a) );
   end
   % if both values > 0, need to find a negative value:
   both_positive = err1 > 0 & not_bracketting_zero;
   ntries=1; % up to 10 tries
   while any(both_positive)
        % err1 < err2, so assuming a monotonic fn, we can decrease err1 by increasing distance from 'b'
       diffab = b(both_positive) - a(both_positive);  
       a(both_positive) = a(both_positive) - diffab; % walk out in steps that double with each iteration
       % recompute the new a values (note, it might be smarter to shift a[n-1] into b when we're done
       err1 = func(a);
       if track_fcount
           fcounter = fcounter + both_positive;
       end
       
       % for a severely ill-behaved function err1 can go in the wrong direction as we move apart, so fix it now
       swap_to_a = (err2 < err1  & not_bracketting_zero);
       [a(swap_to_a), b(swap_to_a)] = deal(b(swap_to_a), a(swap_to_a));
       [err1(swap_to_a), err2(swap_to_a)] = deal( err2(swap_to_a), err1(swap_to_a) );
       err_outside_tol = min(abs(err1), abs(err2)) > tolFn;
       not_bracketting_zero = (sign(err1) == sign(err2) & err_outside_tol);
       both_positive = not_bracketting_zero;
       if any(both_positive) && ntries > 10
            error('Couldn''t find contrapoint in 10 tries!')
        end
        ntries = ntries + 1;
   end
    
%    ntries=1; % if using while loop for b
    both_negative = err2 < 0 & not_bracketting_zero;
    if any(both_negative)
        %  for Ci, A(0) -> Ci >= 0 so  f(0) - 0 >= 0
        b(both_negative) = 0; % just go to zero and be done NOT GENERAL!!!!
        err2 = func(b);  % for now don't calculate on subsets
        if track_fcount
            fcounter = fcounter + both_negative;
        end
    end
    
    recompute_b = true; % we can no longer be certain that s (b) is the best)
end

if ~isempty(corner)
    % special case: guesses that bracket a "corner" may result in very slow convergence
    %  if we replace an endpoint with corner, the remaining iterations should be much simplified
    bracket_corner = sign(corner - a) ~= sign(corner - b); %( a < corner & b > corner) | ( a > corner & b < corner);
    if any(bracket_corner)
        % replace one endpoint with corner:
        x1 = b; % initialize it for the function call; we'll only use the bracket_corner elements
        if length(corner) == 1 % as above, make sure value is same size as x
            corner = repmat(corner, size(b));
        end
        x1(bracket_corner) = corner(bracket_corner);
        
        errCorner = func(x1);
        if track_fcount
            fcounter = fcounter + bracket_corner;
        end
        
        % sort the result into a b so that the sign is preserved
        save_into_b = bracket_corner & (sign(errCorner) == sign(err2)); % error on corner and b have same sign
        save_into_a = bracket_corner & ~save_into_b;
        [a(save_into_a),b(save_into_b)] = deal(x1(save_into_a), x1(save_into_b));
        [err1(save_into_a),err2(save_into_b)] = deal(errCorner(save_into_a), errCorner(save_into_b));

        recompute_b = true; % we can no longer be certain that s (i.e. matches the last call to computeA)
    end
end


% Brent tolerance:
% BRENT: 2 * eps * abs(b)+ toler %NUMERICAL RECIPES: 2 * eps * abs(b)+ 0.5*toler
tolx = 2* max(1, abs(b)).*tolx_1; % MATLAB's criterea (and also Wilkins 2013, except for the factor of 4)
err_outside_tol =  0.5 .* abs(a - b) > tolx & min(abs(err1), abs(err2)) > tolFn ;


% make sure that 'b' is the best 
err1_is_best = abs(err2) > abs(err1); % and therefore needs to be swapped
if any(err1_is_best)
    [a(err1_is_best), b(err1_is_best)] = deal(b(err1_is_best), a(err1_is_best));
    [err1(err1_is_best), err2(err1_is_best)] = deal(err2(err1_is_best), err1(err1_is_best));
    recompute_b = true;
end

% initialize the search array with the current best guess
%[s, err_s] = deal(b, err2); 
ab_gap =  (a - b);  % do this only after we've sorted a and b

% since b is the current best guess, 'a' therefore stands in as the "previous best"
[c, err3] = deal(a, err1);  % the "previous" value of b
best_is_unchanged = abs(err2) == abs(err1);

% initialize additional vectors to nan
%[p, q] = deal(nan(size(b)));  % components of the new increment
[xstep, xstep1] = deal(3*ab_gap); % the step_size of one- or two-steps ago; initialize to prevent triggering on the first round
q = ones(size(b)); % needed for first iteration; don't use nan: we need p/q to be a valid number
p = 0 .* q;
%--------------------------------------------------------
%  MAIN LOOP
%  note: Stage 2013 says we should add   abs(a - b)<= 0.5(a+b)*eps as a stopping condition
%   to avoid machine accuracy issues when abs(x_root) >> 1.
%     fzero.m uses: 0.5*abs(a-b) <= 2.0*tol*max(abs(b),1.0) [adjusting for the fact that c in fzero is our a]
%used_bisection = true(size(b));  % 'mflag': did we use bisection on the previous round (starts true)
counter = 0;
accel_bi = zeros(size(b));
if verbose
    fprintf('init''l ')
    fprintf('%d: a: %9.5g (%9.3g), b:  %9.5g (%9.3g),c: %9.5g (%9.3g), s: %9.3g\n', fcounter, a, err1, b, err2, c, err3, err2); % a, b, c, s);
end
while any( err_outside_tol )
    % 0. Setup
    %*** NOTE: See 2013 Wilkins about an objection to the magnitude of tol in these tests 
    xstep2 = xstep1; % Advance the record of previous step-sizes; Brent calls this "e"
    xstep1 = xstep;
    
    %ab_gap =  (a - b); % done at end of loop since it's needed for the exit test
    p = 0.*p; % clear p, xstep (this is a bit faster than zeros() or even p(:)=0 in R2015b )
    xstep = 0.*xstep;
%    [p, xstep] = deal(zeros(size(b))); % don't use nan, so p/q is a valid number
    %q = ones(size(b)); % don't worry about q: if not used, p/q = 0
    
    %use_bisection = (abs(xstep2) < tol) | (abs_err_s >= abs_err_s1); % if the latest err is larger than the previous iteration give up on interpolation
    %abs_err_s1 = abs_err_s;
    use_bisection = (abs(xstep2) < tolx) | best_is_unchanged ; % if the latest err is larger than the previous iteration give up on interpolation
    
    r2 = err2 ./ err1; % f(b)/f(c) in Brent - see my comments in secant
    
    % the new guess will be stored in 's' using either secant, inverse quadratic, or bisection:
    try_interp = ~use_bisection & err_outside_tol;
    %err3_close_enough = abs(ab_gap)*20 >= abs(err3 - err1); % if 3rd best is too far, just use two 
    quad_is_safe = (err1 ~= err3 & err2 ~= err3); % (note err1 never equals err2 -- since they have opp. signs -- unless they're zero, in which case try_interp is false)

    %-----
    % 1. Inverse quadratic Method
    % when three values are available, use inverse quadratic method (note: fzero has a more clever algorithm?)
    use_quad = try_interp & quad_is_safe; % see prev note about (err1 ~= err2) 
    if any(use_quad)
        % this way is 3x faster than with subsetting!
           % (defining intermediate difference variables doesn't make much if any speed difference.)
%         s1 =  a .* err2 .* err3 ./ ((err1 - err2).*(err1 - err3)) ...
%              + b .* err1 .* err3 ./ ((err2 - err1).*(err2 - err3))  ...
%              + c .* err1 .* err2 ./ ((err3 - err1).*(err3 - err2));
%         s(use_quad) = s1(use_quad);
        r1 = err3 ./ err1; % Brent Section 4, not ALGOL, notation, swapping a/c from Brent's notation
        r3 = err2 ./ err3;
        p = r3 .* ( ab_gap .* r1 .* (r1 - r2) - (b - c).*(r2 - 1) );
        q = (r1 - 1) .* (r2 - 1) .* (r3 - 1);
%         p(~use_quad) = 0; % so xstep = 0
%         q(~use_quad) = 1; % so xstep is not nan
%         p(use_quad) = p1(use_quad); % don't bother; we'll overwrite what needs overwriting
%         q(use_quad) = q1(use_quad);
        if verbose, fprintf('**quad '), end
    end
    
    %-----
    % 2. Secant method
    % I've found no case in which doing secant helps when inv. quad failed (it actually makes things worse)
    %s_test = (quad_is_safe & abs(pq1) >= abs(ab_gap) ) ; % if inverse-quad went too far
    use_secant = try_interp & ~quad_is_safe;
    % secant:  b - f(b) * (b - a) / ( f(b) - (fa) ): derivation: using point-slope of a line solve for y=0:
    %      (y - y1) = m(x - x1):  m = (y2-y1)/(x2-x1);  b = y1 - m x1 = y2 - m x2;  
    %      y = 0 => -y1 = m x - m x1 => -mx = y1 - m x1 => x = -y1/m + x1
    %      0 = mx + b => x = -b/m  = -(y1 - m x1)/m =  -y1/m + x1 OR -(y2 - m x2)/m = x2 - y2/m
    %       x = x1 - y1 (x2 - x1)/(y2 - y1);   OR x2 - y2 (x2 - x1)/(y2 - y1)
    %s(use_secant)	 = b(use_secant) - err2(use_secant) .* (b(use_secant) - a(use_secant)) ./ (err2(use_secant) - err1(use_secant));
    if any(use_secant)
        % NOTE: We only take a secant when a = c, 
        %   so it doesn't matter whether we use err2/err1 or err2/err3
        % p /q = ((a - b) * err2/err1) / ((err1 - err2)/err1)
        %      = err2 (a - b) / (err1 - err2)  -- compare to secant formula
        p1 = ab_gap .* r2;
        p(use_secant) = p1(use_secant);
        q(use_secant) = 1 - r2(use_secant);
        if verbose, fprintf('secant '), end
    end
    
    if any(try_interp)
        %pq1 = p ./ q;    % divide now (before subsetting, since it's faster); 
        %xstep(try_interp) = pq1(try_interp);  % though only for the interpolated points, so far.
        % note: this does not need subsetting because all(p(~try_interp)) = 0; 
        %  (it should be impossible for q=0 since that would imply err1 = err2, but then err isn't out of tolerance
        xstep = p ./ q;    % divide now (before subsetting, since it's faster); 
        xstep(~try_interp) = 0; % to clear nan's
    end
    
    %-----
    % 3. Override any of the above with bisection , depending on values of s
    % 3a: s is NOT between (3 * a + b)/4 and b
    %	 Brent says: new point should be up to 3/4 way between b and a
    %bi_test1a = sign(s - a) == sign(s - b);  % s is outside (a ... b); i.e. just a safety-check (and Dekker's test)
    bi_test1 = ( abs(p) >= 0.75 .* abs(ab_gap .* q) - 0.5*abs(tolx.*q) ) ; % i.e, toler

    % second test done previously
    % third test: |the current xstep| > 1/2 |xstep[n-2]|, i.e. we're moving too far away from the best. but here it's different
    %bi_test3 = abs(xstep) >= 0.5 * abs(xstep2); % or the "safer":
    bi_test3 = abs(p) >= 0.5 * abs(xstep2 .* q); % need to think about p and q here
    
    % update the use_bisection flag with the (p/q)-dependent criteria
    use_bisection = ( use_bisection | bi_test1 | bi_test3)  & err_outside_tol; % 
        
    % now accept any qualifying interpolated point:
    
    if any(use_bisection)  % it's a pretty rare event
       % s(use_bisection) = (a(use_bisection) + b(use_bisection))*0.5;
        m = -ab_gap./(2 + accel_bi); %(1+ 2.^conseqs);
        %s(use_bisection) =  b(use_bisection) - m;  % b + (a-b)/2 = (b + a)/2
        %  set xstep1 so it makes its way into xstep2 (as per Brent 1971)
        [xstep(use_bisection), xstep1(use_bisection)] =  deal(m(use_bisection)); 
        if verbose, fprintf('bisect '), end
    end

    % xstep (d in Brent) is fully updated, now compute the new guess (note xstep=0 if err is within tol)
    s = b - xstep; 
    xstep_too_small = abs(xstep) < tolx & err_outside_tol; % prevent underflow, but it converges faster if we test against tol rather than eps
    if any(xstep_too_small)
        s2 = b + sign(ab_gap) .* tolx;
        s(xstep_too_small) = s2(xstep_too_small);
    end
    %s(err_outside_tol) = s1(err_outside_tol); % preserve values that have already converged
    
    % Quick sanity check
    if ~all(use_secant | use_quad | use_bisection |  ~err_outside_tol)
        error('Somehow, we didn''t update idx: %d\n', find(~(use_secant | use_quad | use_bisection |  ~err_outside_tol)));
    end
    
    %-----
    % compute error-value for s (how far from the fixed point)
    %  this values should always be "better" than either a or b (but not necessarily better than both)
    if subsetLimit <= 0 || mean(err_outside_tol(:)) > subsetLimit
         err_s = func(s);
    else
        err_s(err_outside_tol) = func(s(err_outside_tol));
    end
    if track_fcount
        fcounter = fcounter + err_outside_tol;
    end
    
    
    %fprintf('Laggards: %d\n', sum(err_outside_tol));
    counter = counter + 1;
    if (counter > iter_limit)
        error('iteration limit exceeded');
    end
    %-----------------
    %  Now reorganize a, b, c so b is in best
    %    Also: set conseqs, err_increased, ab_gap for next round
    if all(abs(err_s) < tolFn)
        % converged in y; s hold the full answer; no need to sort, etc.
        b = s;
        err2 = err_s;
        err_outside_tol=false;  % we're done!
        recompute_b = false;
    else
        % first, test that our new guess was an improvement
        %  if the new guess is larger than the old guess then most likely b is close to zero
        %  and we're just whittling away at a. (Alternatively, the function is not monotonic on [a b])
        %  Either way, we will try bisecting on the next round.
        best_is_unchanged = abs(err_s) > abs(err2); % strictly-greater performs what we really mean and in one extreme case (Heaviside) saves us from incorrect behavior
        if accelBisection
            % Furthermore, if we're just chipping away at the far side, let "bisection" divide by increasing integers.
            accel_bi = accel_bi + best_is_unchanged; 
            accel_bi(~best_is_unchanged | (accel_bi >= accelLimit)) = 0; % reset anything that behaved properly; or after 3 tries
            %      limit the number of consecutive forced-bisections (3 may be best for well-behaved fns (but may prevent convergence in singularities);
            %                                                        Inf is probably safest, though with the current method, it's never needed)
            best_is_unchanged =  accel_bi > 0; %best_is_unchanged & ; 
        end
        
        % update previous  states: (Brent's, a, fa)
        % first store prev. round's best in c: (prev round's result is in s right now)
%        d= c; err4 = err3;
        c = b;   err3 = err2; 

        %swap a,b as needed so that all b matches the sign of s
%         s_a_sign_match = (sign(err_s) == sign(err1) ) & err_outside_tol;
%         if any(s_a_sign_match)
%             % s belong in a, so move b into its place. but if we moved b in, b is now in a and c, 
%             % so lets swap 'a' into 'c' to retain three points
%             c(s_a_sign_match) = a(s_a_sign_match);      err3(s_a_sign_match)= err1(s_a_sign_match);
%             a(s_a_sign_match) = b(s_a_sign_match);      err1(s_a_sign_match)= err2(s_a_sign_match);
%             xstep1(s_a_sign_match) = xstep(s_a_sign_match); % makes a very tiny improvement in some cases
%         end
% 
%         % and now we can just copy all of s into b
%         b = s; err2 = err_s;
% 
%         % finally swap a, b if necessary, so that b is the best answer
%         err1_is_best = (abs(err2) > abs(err1)) & err_outside_tol;
%         if any(err1_is_best)
%             [a(err1_is_best), b(err1_is_best)] = deal(b(err1_is_best), a(err1_is_best));
%             [err1(err1_is_best), err2(err1_is_best)] = deal(err2(err1_is_best), err1(err1_is_best));
%         % Adding this makes it much closer to Brent's (i.e. erasing the improvement)
%         %    err_increased = ( abs(err2) >= abs(err1) );
%         end
%         
%         % very slight improvement, though makes Stage's hyperbolic case (singularity at 0) worse
%         d_closer_than_c = abs(c - b) > abs(d - b) & err_outside_tol;
%         if (d_closer_than_c)
%             c(d_closer_than_c) = d(d_closer_than_c);      err3(d_closer_than_c)= err4(d_closer_than_c);
%         end

        s_b_sign_match = (sign(err_s) == sign(err2) );
        err_s_is_best = ( abs(err_s) <= abs(err2) ) & err_outside_tol;
        a_into_b = (s_b_sign_match & ~err_s_is_best) & err_outside_tol;
        if any(a_into_b)
            % move a into b, because we're going to move s into a (note, b has already been moved into c
            b(a_into_b) = a(a_into_b);      err2(a_into_b)= err1(a_into_b);
        end

        b_into_a = (~s_b_sign_match & err_s_is_best); %  & err_outside_tol is redundant
        if any(b_into_a)
            % move b into a because we're going to move s into b; we need to save a into c (since prev-b isn't being lost, but a will)
            c(b_into_a) = a(b_into_a);      err3(b_into_a)= err1(b_into_a);
            a(b_into_a) = b(b_into_a);      err1(b_into_a)= err2(b_into_a);
        end
        % now copy s into a or b
        if any(err_s_is_best)
            b(err_s_is_best) = s(err_s_is_best);      err2(err_s_is_best)= err_s(err_s_is_best);
        end
        err_s_not_best = ~err_s_is_best & err_outside_tol;
        if any(err_s_not_best)
            a(err_s_not_best) = s(err_s_not_best);      err1(err_s_not_best)= err_s(err_s_not_best);
            xstep1(err_s_not_best) = xstep(err_s_not_best); % As per Brent; makes a very tiny improvement in some cases (& a bit worse in others?)
        end
        
        %-----
        % 0. Test if it's time to exit
        ab_gap =  (a - b); % i.e. 2 * m in Brent ( m = 0.5*(a - b) translating his paralance )
        %clear xstep so it can be partially filled w/o carryover of previous step

        %tol = eps(b) + 0.5*toler;  % BRENT: 2 * eps * abs(b)+ toler
        tolx =  2*max(1,abs(b)).*tolx_1;  % w/o max(abs(b),1.0) it fails on test f2_2
        err_outside_tol =  (0.5 .* abs(ab_gap) > tolx   &   abs(err2) > tolFn) ;
        recompute_b = true;
    end
    
    if verbose
        fprintf('%d: a: %9.5g (%9.3g), b:  %9.5g (%9.3g),c: %9.5g (%9.3g), s: %9.3g\n', fcounter, a, err1, b, err2, c, err3, err_s); % a, b, c, s);
    end
end
if recompute_b
    % should never be needed if TolX = 0 and TolFn > 0 (and at least one iter)
    err2 = func(b);
    if track_fcount
        fcounter = fcounter + 1;
    end
end
if verbose, fprintf('Brent p/q preserving C. iterations: %d; fcalls: %d; xval: %0.5g\n', counter, fcounter, b), end

