function [vnew]=count(nvars,v,vmax,id)

% nvars = number of digits
% v     = current vector of digits
% vmax  = maximum values of digits
% id    = starting digit
% vnew  = new vector of digits

i=id;

% starting at id, set digits which are at its maximum equal to 1
% first digit that is not at its maximum is incremented

while v(i)==vmax(i)
    v(i)=1;
    i=rem(i,nvars)+1;
end
v(i)=rem(v(i),vmax(i))+1;
vnew=v;
end