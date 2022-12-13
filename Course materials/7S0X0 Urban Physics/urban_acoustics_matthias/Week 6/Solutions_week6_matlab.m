clear all; close all

freq = [63, 125, 250, 500, 1000, 2000, 4000];  % octave band center frequencies

%% Question 1a: rolling noise and propagation noise

% Define coefficients for the lightweight vehicles (m=1) and the heavy duty
% vehicles (m=3), for each octave band
A_R = [79.7, 85.7, 84.5,  90.2,  97.3, 93.9, 84.1; ...
       87.0, 91.7, 94.1, 100.7, 100.8, 94.3, 87.1];
B_R = [30.0, 41.5, 38.9, 25.7, 32.5, 37.2, 39.0; ...
       30.0, 33.5, 31.3, 25.4, 31.8, 37.1, 38.6];
A_P = [ 94.5,  89.2,  88.0,  85.9,  84.2, 86.9, 83.3; ...
       104.4, 100.6, 101.7, 101.0, 100.1, 95.9, 91.3];
B_P = [-1.3, 7.2, 7.7, 8.0, 8.0, 8.0, 8.0; ...
        0.0, 3.0, 4.6, 5.0, 5.0, 5.0, 5.0,];
% speed of the lightweight and heavy duty vehicles
v = kron(ones(size(freq)), [120; 80]);
vref = 70; % reference speed

% Compute sound power associated with 1 vehicle
L_WR = A_R + B_R.*log10(v/vref) + 0; % rolling noise
L_WP = A_P + B_P.*(v-vref)/vref + 0; % propulsion noise
L_W = 10*log10( 10.^(L_WR/10) + 10.^(L_WP/10) ); % total sound power

% Compute sound power associated with the flow of vehicles
Q = kron(ones(size(freq)), [3600; 400]); % number of vehicles per hour
% sound power per meter, for each  vehicule type
L_WEQ = L_W + 10*log10(Q./(1000*v));
% total sound power per meter, with both types of vehicles combined
L_WEQ_total = 10*log10( 10.^(L_WEQ(1,:)/10) + 10.^(L_WEQ(2,:)/10) );


%% Question 1b: compute the noise level Lp at the closest façade
% Sound level formula: Lp = L_WEQ_total -Adiv -Aatm -Aground -Adif

d = sqrt(4^2 + 110^2); % direct path length between source and receiver

% there is no barrier (so no diffraction), thus
Adif=0;

% geometrical divergence
Adiv = 10*log10(d)+8;

% atmospheric absorption
alpha = [0 0.381 1.13 2.36 4.08 8.75 26.4]; % coefficients
Aatm = alpha*d/1000;

% ground effects
dp = 110; zs = 0; zr = 4; % geometry
Gs = 0; Gpath=0.7; % ground properties
if dp <= 30*(zs+zr) % here this condition is true
    Gprime_path = Gpath*dp/(30*(zs+zr)) + Gs*(1-dp/(30*(zs+zr)));
else
    Gprime_path = Gpath;
end
% from Gprime_path and figure 1, we find the corresponding (approximate) values of
% Agroud for each octave band
Aground = [-1.5 -1.5 -1.5 -1.5 -1.5 4.5 8.5];

% compute noise level
Lp = L_WEQ_total -Adiv -Aatm -Aground -Adif;


%% Question 1c: compute the corresponding A-weighted noise level

% A-weighting correction factor for each octave band
WA = [-26.2, -16.1, -8.6, -3.2, 0, 1.2, 1];

% Logarithmic sum with corrections:
LpA = 10*log10( sum( 10.^((Lp + WA)/10) ) );


%% Question 1d: compute LpA in case of snow

% With the snow, Gprime_path is 0.9. The new values for Aground are thus (approximately)
Aground_snow = [-0.5, -0.5, -0.5, 0, 4, 8, 12];

% Re-compute corresponding sound level:
Lp_snow = L_WEQ_total -Adiv -Aatm -Aground_snow -Adif;

% Apply A correction:
LpA_snow = 10*log10( sum( 10.^((Lp_snow + WA)/10) ) );

% The sound reduction due to the snow is thus
sound_reduction_snow = LpA_snow - LpA;


%% Question 1e: compute LpA with an additional 4m-high barrier at 10m
%% from the road

% Define height of the barrier
height_barrier = 4;

% L_WEQ_total, Adiv and Aatm are the same as before
% Adif is not zero anymore
% the ground effects are now taken into account with Adif, so
Aground_barrier = 0;

% Formula to use: Adif = Ddif_SR + Dground_SO + Dground_OR = ?

% computation of Ddif_SR
SO = sqrt(height_barrier^2+10^2);
OR = 100;
SR = d;
delta = SO + OR - SR; % path difference
lambda = 340./freq; % wavelength associated with each frequency band

for ii=1:length(freq) % loop over the octave bands
    if 40*delta/lambda(ii)>=-2
        Ch = min(freq(ii)*height_barrier/250,1);
        Ddif_SR(ii) = Ch*10*log10(3+40*delta/lambda(ii));
    else
        Ddif_SR(ii) = 0;
    end
end
Ddif_SR = max(Ddif_SR,0); % diffraction cannot increase sound level (see
                          % Cnossos report)

% To compute the ground effects we will need to compute additional
% diffraction terms, between the image-source Sp and the receiver R (Ddif_SpR) and
% between the source S and the image-receiver Rp (Ddif_SRp).
%
% However, here the source is located at zs=0 meter, so the source S and
% the image-source Sp are located at the same position (S==Sp).
% We thus have directly:
Ddif_SpR = Ddif_SR;

% Let's compute the remaining diffraction term, Ddif_SRp
SO = sqrt(height_barrier^2+10^2); % same as before
ORp = sqrt(100^2 + (height_barrier+zr)^2); % distance between top of
                                           % barrier and image-receiver
SRp = sqrt(110^2 + zr^2); % distance between source and image-receiver
delta = SO + ORp - SRp; % path difference (or "detour")
for ii=1:length(freq) % loop over the octave bands (same formula as before)
    if 40*delta/lambda(ii)>=-2
        Ch = min(freq(ii)*height_barrier/250,1);
        Ddif_SRp(ii) = Ch*10*log10(3+40*delta/lambda(ii));
    else
        Ddif_SRp(ii) = 0;
    end
end
Ddif_SRp = max(Ddif_SRp,0); % diffraction cannot increase sound level (see
                            % Cnossos report)

% also, according to Cnossos report, the attenuation due to diffraction
% on a horizontal edge cannot exceed 25 dB, so
Ddif_SR = min(Ddif_SR,25);


% Let us now compute the ground attenuation on the source side
% (Aground_SO) and on the receiver side (Aground_OR).

% On the source side of the barrier, Gs_SO and Gpath_SO are both 0 (dense asphalt), so Gprime_path_SO is simply
% 0. From figure 2, we find that
Aground_SO = -3;

% On the receiver side, since the virtual source is the top of
% the barrier O, we have Gs=Gpath=Gprime_path (see Cnossos report, page
% 95). Therefore, with Gpath=0.7, we directly get from figure 2 (very approximately)
Aground_OR = [-0.9, -0.75, -0.9, -0.9, -0.9, -0.9, -0.9];

% we can now compute the terms Dground_SO and Dground_OR as
Dground_SO = -20*log10(1+(10.^(-Aground_SO/20)-1).*10.^(-(Ddif_SpR-Ddif_SR)/20));
Dground_OR = -20*log10(1+(10.^(-Aground_OR/20)-1).*10.^(-(Ddif_SRp-Ddif_SR)/20));

% Therefore, the attenuation due to the barrier (per octave band) is
Adif_barrier = Ddif_SR + Dground_SO + Dground_OR;

% and the sound level is
Lp_barrier = L_WEQ_total -Adiv -Aatm -Aground_barrier -Adif_barrier;

% We can now compute the A-weighted sound level:
LpA_barrier = 10*log10( sum( 10.^((Lp_barrier + WA)/10) ) );


%% Question 1f: how much would the A-weighted sound level reduce if the
%% speed limit is decreased to 100 km/h?

% We need to recalculate the corresponding sound power (see question 1a)

v = kron(ones(size(freq)), [100; 80]); % speed of the lightweight and heavy duty vehicles
L_WR = A_R + B_R.*log10(v/vref) + 0; % rolling noise
L_WP = A_P + B_P.*(v-vref)/vref + 0; % propulsion noise
L_W = 10*log10( 10.^(L_WR/10) + 10.^(L_WP/10) ); % total sound power
Q = kron(ones(size(freq)), [3600; 400]); % number of vehicles per hour
L_WEQ = L_W + 10*log10(Q./(1000*v)); %sound power per meter, for each vehicule type

% total sound power per meter, with both types of vehicles combined
L_WEQ_total_100 = 10*log10( 10.^(L_WEQ(1,:)/10) + 10.^(L_WEQ(2,:)/10) );

% calculate the difference with the previous sound power, for each octave
% band:
sound_reduction_100 = L_WEQ_total_100-L_WEQ_total;

% modify the sound levels accordingly (just to avoid using the full
% formula again)
Lp_100 = Lp + sound_reduction_100;
Lp_snow_100 = Lp_snow + sound_reduction_100;
Lp_barrier_100 = Lp_barrier + sound_reduction_100;

% compute the A-weighted sound levels
LpA_100 = 10*log10( sum( 10.^((Lp_100 + WA)/10) ) );
LpA_snow_100 = 10*log10( sum( 10.^((Lp_snow_100 + WA)/10) ) );
LpA_barrier_100 = 10*log10( sum( 10.^((Lp_barrier_100 + WA)/10) ) );