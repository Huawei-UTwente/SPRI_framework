% test the reflex control parameter plot

tPhase = 100;
M = 4;
nOpt = 10;
gains_fse = zeros(nOpt, tPhase, M);
gains_lce = zeros(nOpt, tPhase, M);
res = zeros(nOpt, 3);

fse_rand = 5*rand(1, M);
lce_rand = 5*rand(1, M);
res_rand = rand(1, 3);



for opt = 1:nOpt

gains_fse(opt, :, :) = 3*sin(((1:tPhase)')*ones(1, M)*0.1 + fse_rand + rand(1, M));
gains_lce(opt, :, :) = 3*sin(((1:tPhase)')*ones(1, M)*0.1 + lce_rand + rand(1, M));
res(opt, :) = 0.6*res_rand + 0.3*rand(1, 3);

end

hostMus = "Sol";
refMus = ["TA", "Sol", "mGas", "lGas"];

fig = reflexParPlot(gains_fse, gains_lce, res, tPhase, M, hostMus, refMus);