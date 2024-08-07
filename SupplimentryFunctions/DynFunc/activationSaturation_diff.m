function dactSat_dact = activationSaturation_diff(act)
%activationSaturation: saturate the reflex generated activation between 0
%and 1

delta = 1e-4;

% first make sure it is more than 0
act1 = (sqrt(act.^2 + delta) + act)/2;   % 0 < act1 < inf
dact1_dact = (act./sqrt(act.^2 + delta) + 1)/2;

% then transfer act1 to [-inf, 1]
act2 = 1 - act1;   % -inf < act2 < 1
dact2_dact = - dact1_dact;

% then regulate it to be more than 0
% act3 = (sqrt(act2.^2 + delta) + act2)/2;   % 0 < act3 < 1
dact3_dact = ((act2./sqrt(act2.^2 + delta) + 1)/2).*dact2_dact;

% transfer it back to the original order of act
% actSat = 1 - act3;
dactSat_dact = -dact3_dact;

end

