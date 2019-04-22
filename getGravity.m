% Returns the acceleration due to gravity as a function of altitude
function g = getGravity(altitude, g0)
    Re = 3.67e6; % Radius of earth [m]
    g = g0*(Re/(Re + altitude))^2;
end