% Returns the density as a function of altitude
function rho = getDensity(altitude)
    a = 1.2;
    b = 2.9e-5;
    rho = a*exp(-b*altitude^1.15);
end
