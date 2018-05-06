%Function my KSTest

function [H, pValue, KSstatistic] = kstest_local(x, CDF, alpha)

%KSTEST Single sample Kolmogorov-Smirnov goodness-of-fit hypothesis test.
%   H = KSTEST(X,CDF,ALPHA,TAIL) performs a Kolmogorov-Smirnov (K-S) test
%   to determine if a random sample X could have the hypothesized, continuous
%   cumulative distribution function CDF. CDF is optional: if omitted or
%   unspecified (i.e., set to an empty matrix []), the hypothetical c.d.f is
%   assumed to be a standard normal, N(0,1). ALPHA and TAIL are optional
%   scalar inputs: ALPHA is the desired significance level (default = 0.05);
%   TAIL indicates the type of test (default = 'unequal'). H indicates the
%   result of the hypothesis test:
%      H = 0 => Do not reject the null hypothesis at significance level ALPHA.
%      H = 1 => Reject the null hypothesis at significance level ALPHA.
%
%   Let S(x) be the empirical c.d.f. estimated from the sample vector X, F(x)
%   be the corresponding true (but unknown) population c.d.f., and CDF be the
%   known input c.d.f. specified under the null hypothesis.
%
%   The one-sample K-S test tests the null hypothesis that F(x) = CDF for
%   all x against the alternative specified by TAIL:
%       'unequal' -- "F(x) not equal to CDF" (two-sided test)
%       'larger'  -- "F(x) > CDF" (one-sided test)
%       'smaller' -- "F(x) < CDF" (one-sided test)
%
%   For TAIL = 'unequal', 'larger', and 'smaller', the test statistics are
%   max|S(x) - CDF|, max[S(x) - CDF], and max[S(x) - CDF], respectively.
%
%   X is a vector representing a random sample from some underlying
%   distribution. Missing observations in X, indicated by NaNs
%   (Not-a-Number), are ignored.
%
%   CDF is the c.d.f. under the null hypothesis. If specified, it must be an
%   explicit 2-column matrix of paired (x,y) values. Column 1 contains the
%   x-axis data, column 2 the corresponding y-axis c.d.f data. Since the K-S
%   test statistic will occur at one of the observations in X, the calculation
%   is most efficient when CDF is only specified at the observations in X.
%   When column 1 of CDF represents x-axis points independent of X, CDF is
%   're-sampled' at the observations found in the vector X via interpolation.
%   In this case, the interval along the x-axis (the column 1 spread of CDF)
%   must span the observations in X for successful interpolation.
%
%   [H,P] = KSTEST(...) also returns the asymptotic P-value P.
%
%   [H,P,KSSTAT] = KSTEST(...) also returns the K-S test statistic KSSTAT
%   defined above for the test type indicated by TAIL.
%
%   [H,P,KSSTAT,CV] = KSTEST(...) returns the critical value of the test CV.
%
%   The decision to reject the null hypothesis is based on comparing the
%   p-value P with ALPHA, not by comparing the statistic KSSTAT with the
%   critical value CV.  CV is computed separately using an approximate
%   formula or by interpolation in a table.  The formula and table cover
%   the range 0.01<=ALPHA<=0.2 for two-sided tests and 0.005<=ALPHA<=0.1
%   for one-sided tests.  CV is returned as NaN if ALPHA is outside this
%   range.  Since CV is approximate, a comparison of KSSTAT with CV may
%   occasionally lead to a different conclusion than a comparison of P with
%   ALPHA.  
%
%   See also KSTEST2, LILLIETEST, CDFPLOT.

% Copyright 1993-2007 The MathWorks, Inc.
% $Revision: 1.5.2.5 $   $ Date: 1998/01/30 13:45:34 $

% References:
%   Massey, F.J., "The Kolmogorov-Smirnov Test for Goodness of Fit,"
%         Journal of the American Statistical Association, 46 (March 1956), 68-77.
%   Miller, L.H., "Table of Percentage Points of Kolmogorov Statistics,"
%         Journal of the American Statistical Association, (March 1951), 111-121.
%   Marsaglia, G., W.W. Tsang, and J. Wang (2003), "Evaluating Kolmogorov`s
%         Distribution," Journal of Statistical Software, vol. 8, issue 18.

% Get sample cdf, display error message if any
[sampleCDF,x,n,emsg,eid] = cdfcalc(x);

% Check & scrub the hypothesized CDF specified under the null hypothesis.
% If CDF has been specified, remove any rows with NaN's in them and sort
% x-axis data found in column 1 of CDF. If CDF has not been specified, then
% allow the convenience of x ~ N(0,1) under the null hypothesis.
if (nargin >= 2) && ~isempty(CDF)

   CDF  =  CDF(~isnan(sum(CDF,2)),:);

   [xCDF,i] =  sort(CDF(:,1));    % Sort the theoretical CDF.
   yCDF = CDF(i,2);
   
   ydiff = diff(yCDF); 
   % Remove duplicates, but it's an error if they are not consistent
   dups = find(diff(xCDF) == 0);
   if ~isempty(dups)
      if ~all(ydiff(dups) == 0)
         error('stats:kstest:BadCdf',...
               'CDF must not have duplicate X values.');
      end
      xCDF(dups) = [];
      yCDF(dups) = [];
   end

end

% If CDF's x-axis values have been specified by the data sample X, then just
% assign column 2 to the null CDF; if not, then we interpolate subject to the
% check that the x-axis interval of CDF must bound the observations in X. Note
% that both X and CDF have been sorted and have had NaN's removed.
if isequal(x,xCDF)
   nullCDF  =  yCDF;       % CDF has been specified at the observations in X.
else
    
    if (x(1) < xCDF(1)) || (x(end) > xCDF(end))
      fprintf('%f %f %f %f',x(1),xCDF(1),x(end),xCDF(end));
     error('stats:kstest:BadCdf',...
           'Hypothesized CDF matrix must span the observations interval in X.');
   end
   nullCDF  =  interp1(xCDF, yCDF, x, 'linear');
end

% Compute the test statistic of interest.

%  2-sided test: T = max|S(x) - CDF(x)|.
      delta1    =  sampleCDF(1:end-1) - nullCDF;   % Vertical difference at jumps approaching from the LEFT.
      delta2    =  sampleCDF(2:end)   - nullCDF;   % Vertical difference at jumps approaching from the RIGHT.
      deltaCDF  =  abs([delta1 ; delta2]);
      alpha1    =  alpha / 2;

KSstatistic   =  max(deltaCDF);



    s = n*KSstatistic^2;

    % For d values that are in the far tail of the distribution (i.e.
    % p-values > .999), the following lines will speed up the computation
    % significantly, and provide accuracy up to 7 digits.
    if (s > 7.24) ||((s > 3.76) && (n > 99))
        pValue = 2*exp(-(2.000071+.331/sqrt(n)+1.409/n)*s);
    else
        % Express d as d = (k-h)/n, where k is a +ve integer and 0 < h < 1.
        k = ceil(KSstatistic*n);
        h = k - KSstatistic*n;
        m = 2*k-1;

        % Create the H matrix, which describes the CDF, as described in Marsaglia,
        % et al. 
        if m > 1
            c = 1./gamma((1:m)' + 1);

            r = zeros(1,m);
            r(1) = 1; 
            r(2) = 1;

            T = toeplitz(c,r);

            T(:,1) = T(:,1) - (h.^(1:m)')./gamma((1:m)' + 1);

            T(m,:) = fliplr(T(:,1)');
            T(m,1) = (1 - 2*h^m + max(0,2*h-1)^m)/gamma(m+1);
         else
             T = (1 - 2*h^m + max(0,2*h-1)^m)/gamma(m+1);
         end

        % Scaling before raising the matrix to a power
        if ~isscalar(T)
            lmax = max(eig(T));
            T = (T./lmax)^n;
        else
            lmax = 1;
        end

        % Pr(Dn < d) = n!/n * tkk ,  where tkk is the kth element of Tn = T^n.
        % p-value = Pr(Dn > d) = 1-Pr(Dn < d)
        pValue = (1 - exp(gammaln(n+1) + n*log(lmax) - n*log(n)) * T(k,k));
    end

    % "H = 0" implies that we "Do not reject the null hypothesis at the
    % significance level of alpha," and "H = 1" implies that we "Reject null
    % hypothesis at significance level of alpha."

    H  =  (pValue < alpha1);

end