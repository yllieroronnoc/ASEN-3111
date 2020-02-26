function [xb,yb] = NACA_Airfoils(m, p, t, c, N)
%{
Auuthor: Connor O'Reilly
Collaborators:
Last Revision: 11/6/2019

Inputs:
    m:  maximum camber
    p: location of maximum camber
    t: thickness
    c: chord length
    N: number of employed panels
Outputs: 
    xb and yb are double vectors containing the locations of the 
    boundary points.
%}

%initialization
m = m/100;
p = p/10;
t = t/100;
xloc = linspace(0,c,N);


%set value for easy reading
frac_xc = xloc./c; %x/c


%thickness distribution of the airfoil from the mean camber line
y_t = (t*c/0.2).*( 0.2969 .* sqrt(frac_xc) - 0.1260 .*(frac_xc) ...
                   - 0.3516 .*(frac_xc).^2 + 0.2843 .*(frac_xc).^3 ...
                    - 0.1036 .* (frac_xc).^4);

%mean camber line
if(m ~= 0) && (p ~= 0)
    for i = 1:length(xloc)
        if (0 <= xloc(i)) && (xloc(i) <= p*c)
            y_c(i) = xloc(i) * (m/p^2) * (2*p - (xloc(i)/c));
        else 
            y_c(i) = (m/(1-p)^2) * (c-xloc(i)) * (1 + (xloc(i)/c) - 2*p); 
        end
    end
else
    y_c = zeros(1, length(xloc));
end 
%calculate zeta
zeti = atan(diff(y_c));
%append to end of vector
zeti = [zeti,0];

%calculate coordinates of upper and lower surface

%upper
xu = xloc - (y_t.*sin(zeti));
yu = y_c + y_t.*cos(zeti);

%lower
xl =  xloc + (y_t.*sin(zeti));
yl = y_c - y_t.*cos(zeti);

%combine both upper and lower coords to get full surface
xb = [flip(xl),xu(2:end)];
yb = [flip(yl),yu(2:end)];

end

