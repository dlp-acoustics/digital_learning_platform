clear all; close all

freq = [63, 125, 250, 500, 1000, 2000, 4000];  % octave band
lambda = 340./freq; % wavelength

% The beginning of the script corresponds to the exercises of
% week 6. The exercises for week 7 starts at line 194

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



%% Question 2a: compute the A-weighted sound level with a car parking
%% space of 10 meter width in front of the residential area

% sound level: Lp_parking = L_WEQ_total -Adiv -Aatm -Aground_parking;
% Aground_parking = ?   (the other terms are the same as before, without
% the barrier)

dp = 110; zs = 0; zr = 4; % geometry
Gs = 0; % ground properties near the source
Gpath= (100*0.7 + 10*0.3)/dp; % ground properties along the propagation path
Gprime_path =  Gpath*dp/(30*(zs+zr)) + Gs*(1-dp/(30*(zs+zr)));
% From Gprime_path, we find
Aground_parking = [-1.7, -1.7, -1.7, -1.7, -1.7, 4.3, 1.0];

% sound level per otave band
Lp_parking = L_WEQ_total -Adiv -Aatm -Aground_parking;

% A-weighted sound level
LpA_parking = 10*log10( sum( 10.^((Lp_parking + WA)/10) ) );


%% Question 2b: A-weighted sound level in favorable conditions (without barrier)
% ground effects
dp = 110; zs = 0; zr = 4; % geometry (same as before)
Gs = 0; Gpath=0.7; % ground properties (same as before)
Gprime_path =  Gpath*dp/(30*(zs+zr)) + Gs*(1-dp/(30*(zs+zr))); % (same as before)
% from Gprime_path and figure 1, we find the corresponding values of
% Agroud for each octave band
AgroundF = [-1 -1 -1 -1 -1 5 0.5];

% sound level per octave band
Lp_F = L_WEQ_total -Adiv -Aatm -AgroundF;

% A-weighted sound level
LpA_F = 10*log10( sum( 10.^((Lp_F + WA)/10) ) );

%% Question 2c: A-weighted sound level in favorable conditions with a
%% barrier

% computation of DdifF_SR
gamma = max(1000,SR);
SOF = 2*gamma*asin(SO/2/gamma);
ORF = 2*gamma*asin(OR/2/gamma);
SRF = 2*gamma*asin(SR/2/gamma);
deltaF = SOF + ORF - SRF; % path length difference

for ii=1:length(freq) % loop over the octave bands
    if 40*deltaF/lambda(ii)>=-2
        Ch = min(freq(ii)*height_barrier/250,1);
        DdifF_SR(ii) = Ch*10*log10(3+40*deltaF/lambda(ii));
    else
        DdifF_SR(ii) = 0;
    end
end
% The values of Ddif should be bound:
DdifF_SR = max(DdifF_SR,0); % diffraction cannot increase sound level (see
                            % Cnossos report)

% computation of DdifF_SpR
DdifF_SpR = DdifF_SR; % since the source is located at zs=0, those two terms
                        % are the same (same as question 1e)

% computation of DdifF_SRp. This is the same as before, just replace R with
% Rp in the calculations for DdifF_SR.
gamma = max(1000,SRp);
ORpF = 2*gamma*asin(ORp/2/gamma);
SRpF = 2*gamma*asin(SRp/2/gamma);
delta = SOF + ORpF - SRpF; % path length difference
for ii=1:length(freq) % loop over the octave bands
    if 40*delta/lambda(ii)>=-2
        Ch = min(freq(ii)*height_barrier/250,1);
        DdifF_SRp(ii) = Ch*10*log10(3+40*delta/lambda(ii));
    else
        DdifF_SRp(ii) = 0;
    end
end
DdifF_SRp = max(DdifF_SRp,0);

% also, according to Cnossos report, the attenuation due to diffraction
% on a horizontal edge cannot exceed 25 dB, so
DdifF_SR = min(DdifF_SR,25);

% computation of ground effects on the source side AgroundF_SO: we first
% need to compute Gprime_path (it is the same as for the homogeneous
% conditions). On the source side of the barrier, Gs and Gpath are
% both 0 (dense asphalt), so Gprime_path is simply 0. From figure 2,
% we find that
AgroundF_SO = -3;

% computation of of ground effects on the receiver side AgroundF_OR: since
% the virtual (or secondary) source is here the top of the barrier O, we
% have Gs=Gpath=Gprime_path (see Cnossos report, page 95). Therefore, with
% Gpath=0.7, we directly get from figure 2 (approximately)
AgroundF_OR = [-1 -0.75 -1 -1 -1 -1 -1];

% we can now compute
DgroundF_SO = -20*log10(1+(10.^(-AgroundF_SO/20)-1).*10.^(-(DdifF_SpR-DdifF_SR)/20));
DgroundF_OR = -20*log10(1+(10.^(-AgroundF_OR/20)-1).*10.^(-(DdifF_SRp-DdifF_SR)/20));

% and finally we can get the attenuation due to the barrier
Adif_barrierF = DdifF_SR + DgroundF_SO + DgroundF_OR;

% sound level per octave band
Lp_barrierF = Lp_barrier + Adif_barrier - Adif_barrierF;

% A-weighted sound level
LpA_barrierF = 10*log10( sum( 10.^((Lp_barrierF + WA)/10) ) );


%% Question 2e: A-weighted sound level with 2 barriers in homogeneous
%% conditions

% The total sound level corresponds to the contributions of the sound
% levels associated with each path.

% For path 1 (see slides), we can use the sound level computed previously
% (for week 6), i.e.,
Lp_path1 = Lp_barrier;

% For path 2:

% we compute sound power of the image source Sp (Op corresponds to the top of
% the second barrier)

%geometry
SpR = sqrt(130^2 + 4^2); % distance between image-source and receiver
OOp = 20; % distance between the two barriers
SpO = sqrt(height_barrier^2+30^2); % distance between image-source and
                                   % top of first barrier
SOp = sqrt(height_barrier^2+10^2); % distance between source and top of
                                   % the second barrier

deltap = SpO - SOp - OOp; % detour

% retrodiffraction coefficient
for ii=1:length(freq) % loop over the octave bands
    if 40*deltap/lambda(ii)>=-2
        Ch = min(freq(ii)*height_barrier/250,1);
        Dretrodif(ii) = Ch*10*log10(3+40*deltap/lambda(ii));
    else
        Dretrodif(ii) = 0;
    end
end

% the attenuation due to diffraction shall be bounded:
Dretrodif = max(Dretrodif,0);
Dretrodif = min(Dretrodif,25);

% Compute the sound power of the image-source:
alpha_barrier = 0; % absorption coefficient of the barrier
L_WEQ_image = L_WEQ_total + 10*log10(1-alpha_barrier) - Dretrodif;


% now we need to compute the other terms

% geometrical divergence
Adiv_image = 10*log10(SpR)+8;

% atmospheric absorption
alpha = [0 0.381 1.13 2.36 4.08 8.75 26.4]; % coefficients
Aatm_image = alpha*SpR/1000;

% diffraction by the first barrier
delta = SpO + OR - SpR; % path difference
for ii=1:length(freq) % loop over the octave bands
    if 40*delta/lambda(ii)>=-2
        Ch = min(freq(ii)*height_barrier/250,1);
        Ddif_image_SR(ii) = Ch*10*log10(3+40*delta/lambda(ii));
    else
        Ddif_image_SR(ii) = 0;
    end
end
Ddif_image_SR = max(Ddif_image_SR,0);

Ddif_image_SpR = Ddif_image_SR; % because the image source is located at zs=0

Ddif_image_SR = min(Ddif_image_SR,25); % we put this after because this
                                       % bound doesn't apply to
                                       % Ddif_image_SpR (see Cnossos report)

% Compute the remaining diffraction term
SpRp = SpR; % distance between image of image-source and
            % image-receiver
delta = SpO + ORp - SpRp; % path difference
for ii=1:length(freq) % loop over the octave bands
    if 40*delta/lambda(ii)>=-2
        Ch = min(freq(ii)*height_barrier/250,1);
        Ddif_image_SRp(ii) = Ch*10*log10(3+40*delta/lambda(ii));
    else
        Ddif_image_SRp(ii) = 0;
    end
end
Ddif_image_SRp = max(Ddif_image_SRp,0);
Ddif_image_SRp = min(Ddif_image_SRp,25);


% ground

% Aground_OR is the same as previously
% So we just need to estimate Aground_SpO, the ground attenuation on the
% source side of the barrier.
% Here, Gpath=0 and Gs=0 so Gprime_path=0, and we find from figure 4
Aground_SpO = -3;
% for all octave bands

% we can thus compute the remaining terms
Dground_SpO = -20*log10(1+(10.^(-Aground_SpO/20)-1).*10.^(-(Ddif_image_SpR-Ddif_image_SR)/20));
Dground_OR = -20*log10(1+(10.^(-Aground_OR/20)-1).*10.^(-(Ddif_image_SRp-Ddif_image_SR)/20));
Adif_image = Ddif_image_SR + Dground_SpO + Dground_OR;

% sound level associated with path 2:
Lp_path2 = L_WEQ_image -Adiv_image -Aatm_image -Adif_image;


% the total sound level is obtained by combining the 2 sound paths
% (logarithmic sum):
Lp_2barriers = 10*log10( 10.^(Lp_path1/10) + 10.^(Lp_path2/10) );

% A-weighted sound level
LpA_2barriers = 10*log10( sum( 10.^((Lp_2barriers + WA)/10) ) );


%% Question 2f: same as question e with a different absorption
%% coefficient

% with an absorption coefficient of the barriers
alpha_barrier = 0.9;
% the sound power of the source image is now
L_WEQ_image_new = L_WEQ_total + 10*log10(1-alpha_barrier) - Dretrodif;

% We need to recompute the sound level associated with path 2:
Lp_path2_new = L_WEQ_image_new -Adiv_image -Aatm_image -Adif_image;

% the total sound level is then
Lp_2barriers_new = 10*log10( 10.^(Lp_path1/10) + 10.^(Lp_path2_new/10) );

% A-weighted sound level
LpA_2barriers_new = 10*log10( sum( 10.^((Lp_2barriers_new + WA)/10) ) );