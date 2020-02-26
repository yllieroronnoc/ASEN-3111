function cl = Vortex_Panel(xb,yb,V_inf,alpha,plot_cp)
%{
Author: Connor T. O'Reilly
Collaborators:
Last Revision: 11/7/2019

Inputs:
    x: vector containing boundary points in panel method
    y: vector containing the y-location of coordinates
    V_inf: free-stream flow speed
    alpha: angle of attack (deg)
Outputs: 
    cl: section coefficient of lift
    %}
    %obtain number of vortex panels
    m = length(xb)-1;

    %initialize for optimization
    x = zeros(1,m);
    y = zeros(1,m);
    s = zeros(1,m);
    sine = zeros(1,m);
    cosine = zeros(1,m);
    theta = zeros(1,m);
    v = zeros(1,m);
    cp = zeros(1,m);
    rhs = zeros(1,m);
    cn1 = zeros(m,m);
    cn2 = zeros(m,m);
    ct1 = zeros(m,m);
    ct2 = zeros(m,m);
    an = zeros(m+1,m+1);
    at = zeros(m,m+1);
    
    %initialize other vals for cl
    rho = 1.225; %kg/m^3
    dyn = 0.5 * rho * V_inf^2;

    mp1 = m + 1;

    %convert AOA to radians
    alpha = deg2rad(alpha);


    for i = 1:m
        ip1 = i + 1;
        x(i) = 0.5*(xb(i) + xb(ip1));
        y(i) = 0.5*(yb(i) + yb(ip1));
        s(i) = sqrt( (xb(ip1) - xb(i))^2 + (yb(ip1) - yb(i))^2 );
        theta(i) = atan2( (yb(ip1)-yb(i)), (xb(ip1)-xb(i)) );
        sine(i) = sin(theta(i));
        cosine(i) = cos(theta(i));
        rhs(i) = sin(theta(i) - alpha);
    end

    for i = 1:m
        for j = 1:m
            if(i == j)
                cn1(i,j) = -1.0;
                cn2(i,j) = 1.0;
                ct1(i,j) = 0.5*pi;
                ct2(i,j) = 0.5*pi;
            else
                a = -(x(i)-xb(j)) * cosine(j) - (y(i) - yb(j) ) * sine(j);
                b = ( x(i) - xb(j) )^2 + ( y(i) - yb(j) )^2;
                c = sin( theta(i) - theta(j) );
                d = cos( theta(i) - theta(j) );
                e = ( x(i) - xb(j) ) * sine(j) - ( y(i) - yb(j) ) * cosine(j);
                f = log( 1.0 + ( s(j)^2 + 2*a*s(j) )/b );
                g = atan2( e*s(j), b + a*s(j) );
                p = (x(i) - xb(j)) * sin( theta(i) - 2*theta(j) ) + ...
                    (y(i) - yb(j)) * cos( theta(i) - 2*theta(j) );
                q = (x(i) - xb(j)) * cos( theta(i) - 2*theta(j) ) - ...
                    (y(i) - yb(j)) * sin( theta(i) - 2*theta(j) );
                cn2(i,j) = d + 0.5*q*f/s(j) - (a*c + d*e)*g/s(j);
                cn1(i,j) = 0.5*d*f + c*g - cn2(i,j);
                ct2(i,j) = c + 0.5*p*f/s(j) + (a*d - c*e)*g/s(j);
                ct1(i,j) = 0.5*c*f - d*g - ct2(i,j);
            end
        end
    end

    for i = 1:m
        an(i,1) = cn1(i,1);
        an(i,mp1) = cn2(i,m);
        at(i,1) = ct1(i,1);
        at(i,mp1) = ct2(i,m);
        for j = 2:m
            an(i,j) = cn1(i,j) + cn2(i,j-1);
            at(i,j) = ct1(i,j) + ct2(i,j-1);
        end
    end

    an(mp1,1) = 1.0;
    an(mp1,mp1) = 1.0;

    for j = 2:m
        an(mp1,j) = 0;
    end

    rhs(mp1) = 0.0;
    gama = an\rhs';
   % [~, ~, gama, ~] = cramer(an, rhs, gama, mp1, m);
    
    for i = 1:m
        v(i) = cos( theta(i) - alpha);
        for j = 1:mp1
            v(i) = v(i) + at(i,j) * gama(j);
            cp(i) = 1.0 - v(i)^2;
        end
    end  
    %add value to cp for closed loop 
%     cp(end+1) = cp(1);
    
    %caclulate cl using kutta
    
    %infinitesimal length ds along airfoil
    for i = 1:length(xb)-1
        ds(i,1) = sqrt((xb(i) - xb(i+1))^2 + (yb(i) - yb(i+1))^2);
    end
    %get arrays to be same size
    ds(end+1) = 0;
    
    %calculate circ from vortex sheet strength
    GAMA =  sum(2 *  pi * V_inf * (gama.*ds));
    
    %kutta
    lift = rho * V_inf * GAMA;
    
    %calc chord
    chord = max(xb) - min(xb);
    
    %calc cl
    cl = lift/(dyn*chord);
    if(plot_cp == true)
        %get panel length
        %plot cp against x/c
        yeet = num2str(alpha);
        figure('name',yeet)
        plot(xb(2:end),cp)
        set(gca, 'YDir','reverse')
        xlabel('x/c')
        ylabel('c_p')
        title("C_p vs x/c for an AOA of " + alpha * (180/pi) + "deg")
    end
    
end


