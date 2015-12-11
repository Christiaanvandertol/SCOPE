 function [M] = aggreg(atmfile,SCOPEspec)

% Aggregate MODTRAN data over SCOPE bands by averaging (over rectangular band
% passes)

% Read .atm file with MODTRAN data
s   = importdata(atmfile);
wlM = s.data(:,2);
T   = s.data(:,3:20);

% Extract 6 relevant columns from T

%  1: <Eso*costts/pi>
%  3: <rdd>
%  4: <tss>
%  5: <tsd>
% 12: <tssrdd>
% 16: <La(b)>

U     = [T(:,1) T(:,3) T(:,4) T(:,5) T(:,12) T(:,16)];

nwM   = length(wlM);

nreg  = SCOPEspec.nreg;
streg = SCOPEspec.start;
enreg = SCOPEspec.end;
width = SCOPEspec.res;

% Nr. of bands in each region

nwreg = int32((enreg-streg)./width)+1;

off   = int32(zeros(nreg,1));

for i=2:nreg
    off(i) = off(i-1)+nwreg(i-1);
end

nwS = sum(nwreg);
n   = zeros(nwS,1);    % Count of MODTRAN data contributing to a band
S   = zeros(nwS,6);    % Intialize sums

%k   = int32(0);
j   = int32(zeros(nreg,1));  % Band index within regions

for iwl = 1:nwM
    w   = wlM(iwl);    % MODTRAN wavelength
    for r = 1:nreg
        j(r) = int32(round(w-streg(r))./(width(r)))+1;
        if j(r)>0 && j(r)<=nwreg(r)                 % test if index is in valid range
            k      = j(r)+off(r);                   % SCOPE band index
            S(k,:) = S(k,:)+U(iwl,:);               % Accumulate from contributing MODTRAN data
            n(k)   = n(k)+1;                        % Increment count
        end
    end
end

M = zeros(size(S,1),6);
for i = 1:6
    M(:,i) = S(:,i)./n;      % Calculate averages per SCOPE band
end

end



