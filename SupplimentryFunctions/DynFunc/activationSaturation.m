function actSat = activationSaturation(act)
%activationSaturation: saturate the reflex generated activation between 0
%and 1

delta = 1e-4;

% first make sure it is more than 0
act1 = (sqrt(act.^2 + delta) + act)/2;   % 0 < act1 < inf

% then transfer act1 to [-inf, 1]
act2 = 1 - act1;   % -inf < act2 < 1

% then regulate it to be more than 0
act3 = (sqrt(act2.^2 + delta) + act2)/2;   % 0 < act3 < 1

% transfer it back to the original order of act
actSat = 1 - act3;

end

