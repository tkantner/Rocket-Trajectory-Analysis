% Returns the linearly interpolated thrust given the thrust profile and the
% current time
function T = getThrust(time, THRUST)
    %Loop through the time
    for i = 1:size(THRUST,1)
        if(THRUST(i, 1) == time) %If exact time, return it
            T = THRUST(i, 1);
            return;
        %Else we need to linearly interpolate
        elseif(time < THRUST(i + 1, 1) && time > THRUST(i, 1)) 
            y0 = THRUST(i, 2);
            x0 = THRUST(i, 1);
            y1 = THRUST(i + 1, 2);
            x1 = THRUST(i + 1, 1);
            T = (y0*(x1 - time) + y1*(time - x0))/(x1 - x0);
            return;
        end
    end
    T = 0;
    return;
end