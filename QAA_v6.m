function [a, bb, adg443, adg, aph] = QAA_v6_v2(Rrs,wl)
% input Rrs, wavelength (412, 443, 490, 555, 670).
% references: QAA_v6 pdf file
% output: a(λ) bbp(λ) ag443 adg(λ) aph(λ)

idx412 = find(abs(wl-412)==min(abs(wl-412))); 
idx443 = find(abs(wl-443)==min(abs(wl-443))); 
idx490 = find(abs(wl-490)==min(abs(wl-490))); 
idx55x = find(abs(wl-550)==min(abs(wl-550))); 
idx670 = find(abs(wl-670)==min(abs(wl-670))); 


%% aw & bbw

load('a_bb_water.mat', 'Lee_Zhang')
aw = interp1(Lee_Zhang.wl', Lee_Zhang.aw', wl, 'spline');
bbw = interp1(Lee_Zhang.wl', Lee_Zhang.bbw', wl, 'spline');

bbw55x = bbw(:, idx55x);
[nsample, nband] = size(Rrs);
repbbw = repmat(bbw, nsample, 1);

%% step1
rrs = Rrs ./ (0.52 + 1.7 .* Rrs);
g0 = 0.0895; 
g1 = 0.1247;


u = (-g0 + (g0.^2 + 4*g1.*rrs).^0.5) ./ (2.*g1);

%% step 2-5: get bbp(λ)

bbp = zeros(nsample, nband);
% 
idx_clear_water = Rrs(:, idx670) < 0.0015;  % 清水条件
idx_turbid_water = ~idx_clear_water;           % 浑水条件

% Step 2: 预先计算所需的rrs和u值
rrs443 = rrs(:, idx443);
rrs490 = rrs(:, idx490);
rrs55x = rrs(:, idx55x);
rrs670 = rrs(:, idx670);
u55x = u(:, idx55x);
u670 = u(:, idx670);

% Step 3: 处理清水情况 (Rrs670 < 0.0015)
kappa = log10((rrs443(idx_clear_water) + rrs490(idx_clear_water)) ./ ...
              (rrs55x(idx_clear_water) + 5 * (rrs670(idx_clear_water) ./ rrs490(idx_clear_water)) .* rrs670(idx_clear_water)));

h0 = -1.146; 
h1 = -1.366;
h2 = -0.469;
a555 = aw(:, idx55x) + 10 .^ (h0 + h1 .* kappa + h2 .* kappa .^ 2);

bbp555 = (u55x(idx_clear_water) .* a555) ./ (1 - u55x(idx_clear_water)) - bbw55x;

yita_clear = 2 * (1 - 1.2 .* exp(-0.9 .* (rrs443(idx_clear_water) ./ rrs55x(idx_clear_water))));
bbp(idx_clear_water, :) = bbp555 .* (wl(idx55x) ./ wl) .^ yita_clear;

% Step 4: 处理浑水情况 (Rrs670 >= 0.0015)
a670 = aw(:, idx670) + 0.39 * (Rrs(idx_turbid_water, idx670) ./ ...
       (Rrs(idx_turbid_water, idx443) + Rrs(idx_turbid_water, idx490))) .^ 1.14;

bbp670 = (u670(idx_turbid_water) .* a670) ./ (1 - u670(idx_turbid_water)) - bbw(:, idx670);

yita_turbid = 2 * (1 - 1.2 .* exp(-0.9 .* (rrs443(idx_turbid_water) ./ rrs55x(idx_turbid_water))));
bbp(idx_turbid_water, :) = bbp670 .* (670 ./ wl) .^ yita_turbid;

% Step 5: bb的计算
bb = repbbw + bbp;

%% step 6 get a(λ)
a = (1 - u) .* (repbbw + bbp) ./ u;

%% step 7 & 8
epsilon1 = 0.74 + 0.2 ./ (0.8+rrs(:, idx443)./rrs(:, idx55x));
S = 0.015 + 0.002./(0.6+rrs(:, idx443)./rrs(:, idx55x));
epsilon2 = exp(S*(442.5-415.5));

%% step 9 & 10
a412 = a(:, idx412);
a443 = a(:, idx443);
adg443 =  (a412 - epsilon1 .* a443) ./ (epsilon2 - epsilon1) - (aw(:, idx412) - epsilon1.*aw(:, idx443)) ./ (epsilon2-epsilon1);
adg = adg443 .* exp(-S.*(wl - 443));
aph = a - adg - aw;

end

