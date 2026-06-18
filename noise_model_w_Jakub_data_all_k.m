

%% COMPUTING STEADY AND UNSTEADY LOADING NOISE WITH HANSON FORMULA - MODEL INPUTS
% 
% Emma Vella - 11/2024
%
%
%

%%

clear all
close all
clc

%% reading Jakub data - test case

precision = 'double';
addpath('/home/daep/e.vella/Documents/phD/MATLAB/postpro_exp_crossflow/easyh5-master/');
addpath('/home/daep/e.vella/Documents/phD/2026/Jakub_comp_Hanson/propeller_noise_tests-version1')

file_loading = loadh5('test_loading_prediction.h5');
file_acou = loadh5('test_propagator.h5');





%% INPUT PARAMETERS

B = double(file_loading.inputs.B);
% RPM = 8000;
omega = file_loading.inputs.Omega_rad_p_s;
f0 = omega./(2.*pi);
c0 = file_loading.inputs.sos_m_p_s;


% Hanson params
R = 1.62; % distance rotor-observer [m]
Theta = (0:5:180)*pi/180; 
phi = 90*pi/180; % azimuth angle
phi = pi-phi;
pref = 2e-5;


%% INPUT DATA FILES 


r = file_loading.inputs.radius_inner_m;
dr = abs(file_loading.inputs.radius_outer_m(1)-file_loading.inputs.radius_outer_m(2));

Fz_tour_bis = file_loading.inputs.F_z_prime_N_p_m.*dr;
Fy_tour_bis = file_loading.inputs.F_phi_prime_N_p_m.*dr;

% rotor
Fz_k_1 = file_acou.inputs.blade_loading_harmonics_N_p_m.real(2,2:end,:)...
       + 1i.*file_acou.inputs.blade_loading_harmonics_N_p_m.imag(2,2:end,:);
Fphi_k_1 = file_acou.inputs.blade_loading_harmonics_N_p_m.real(3,2:end,:)...
         + 1i.*file_acou.inputs.blade_loading_harmonics_N_p_m.imag(3,2:end,:);

Fz_k_1 = Fz_k_1.*dr;
Fphi_k_1 = Fphi_k_1.*dr;

% beam
Fk_x_1 = file_acou.inputs.beam_loading_harmonics_N_p_m.real(2,3:2:end,:)...
       + 1i.*file_acou.inputs.beam_loading_harmonics_N_p_m.imag(2,3:2:end,:);
Fk_phi_1 = file_acou.inputs.beam_loading_harmonics_N_p_m.real(3,3:2:end,:)...
         + 1i.*file_acou.inputs.beam_loading_harmonics_N_p_m.imag(3,3:2:end,:);

% Fk_x_1 = file_acou.inputs.beam_loading_harmonics_N_p_m.real(2,2:end,:)...
%        + 1i.*file_acou.inputs.beam_loading_harmonics_N_p_m.imag(2,2:end,:);
% Fk_phi_1 = file_acou.inputs.beam_loading_harmonics_N_p_m.real(3,2:end,:)...
%          + 1i.*file_acou.inputs.beam_loading_harmonics_N_p_m.imag(3,2:end,:);

Fk_x_1 = Fk_x_1.*dr;
Fk_phi_1 = Fk_phi_1.*dr;


%% HANSON - FAR-FIELD FORMULATION 

K = -size(Fz_k_1,2):size(Fz_k_1,2);

for m = 1:10 % sound harmonics

    % unsteady loading noise 
    for k = 1:size(K,2) % loading harmonics

        for idx_r = 1:length(r)
    
            n = m*B-K(k); 
            km = m.*B.*omega./c0;

            if K(k) < 0

                Fz_k(1,k,idx_r) = conj(Fz_k_1(1,abs(K(k)),idx_r));
                Fphi_k(1,k,idx_r) = conj(Fphi_k_1(1,abs(K(k)),idx_r));

            elseif K(k) == 0

                Fz_k(1,k,idx_r) = 0;
                Fphi_k(1,k,idx_r) = 0;

            elseif K(k) > 0

                Fz_k(1,k,idx_r) = Fz_k_1(1,abs(K(k)),idx_r);
                Fphi_k(1,k,idx_r) = Fphi_k_1(1,abs(K(k)),idx_r);

            end


            Pm_r_unst_rotor(k,idx_r,:) = 1i.*km.*B./(4*pi*R).*exp(1i.*n.*(phi-pi/2)).*exp(1i.*km.*R)...
                                    .*(Fz_k(1,k,idx_r).*cos(Theta)-n./(km.*r(idx_r)).*Fphi_k(1,k,idx_r)).*besselj(n,km.*r(idx_r).*sin(Theta)); % pas de multiplication par dr

            Pm_r_unst_beam(k,idx_r,:) = 1i.*km./(4*pi*R).*exp(1i.*n.*(phi-pi/2)).*exp(1i.*km.*R)...
                                    .*(Fk_x_1(1,m,idx_r).*cos(Theta)-n./(km.*r(idx_r)).*Fk_phi_1(1,m,idx_r)).*besselj(n,km.*r(idx_r).*sin(Theta)); % pas de multiplication par dr



%             % formulation 1.44 for beam noise
%             km = m.*omega./c0;
%             Pm_r_unst_beam(k,idx_r,:) = 1i.*km./(4*pi*R).*exp(-1i.*km.*r(idx_r).*sin(Theta).*cos(phi)).*exp(1i.*km.*R)...
%                                         .*(Fk_x_1(1,m,idx_r).*cos(Theta)-Fk_phi_1(1,m,idx_r).*sin(Theta).*sin(phi));



        end

        Pm_rotor(k,:) = sum(Pm_r_unst_rotor(k,:,:),2);
        Pm_beam(k,:) = sum(Pm_r_unst_beam(k,:,:),2);

    end

    % steady loading noise 
    for idx_r = 1:length(r)

        n = m*B;
        Pm_r_st(idx_r,:) = 1i.*km.*B./(4*pi*R).*exp(1i.*n.*(phi-pi/2)).*exp(1i.*m.*B.*omega.*R./c0)...
                                .*(Fz_tour_bis(idx_r).*cos(Theta)-n./(km.*r(idx_r)).*Fy_tour_bis(idx_r)).*besselj(n,km.*r(idx_r).*sin(Theta));
    
    end

    Pm_r_unst(m,:) = sum(Pm_rotor,1);
    Pm_r_steady(m,:) = sum(Pm_r_st,1); 
    Pm_b(m,:) = sum(Pm_beam,1); % beam noise 

end




%

p_f_rotor_unst = Pm_r_unst;
p_ff_rotor_unst = abs(p_f_rotor_unst).*sqrt(2);
SPL_rotor_unst = 20.*log10(p_ff_rotor_unst./pref);

p_f_rotor_st = Pm_r_steady;
p_ff_rotor_st = abs(p_f_rotor_st).*sqrt(2);
SPL_rotor_st = 20.*log10(p_ff_rotor_st./pref);

p_f_beam = Pm_b;
p_ff_beam = abs(p_f_beam).*sqrt(2);
SPL_beam = 20.*log10(p_ff_beam./pref);

%


for m = 1:10
    BPF(m) = m*B*f0;
end


idx_angle = 19; % 19 pour 90° (plan rotor) et 27 pour 130° (-40°) et 9 pour 40°
figure;
hold on
plot(BPF,SPL_rotor_st(:,idx_angle),'o')
plot(BPF,SPL_rotor_unst(:,idx_angle),'o')
plot(BPF,SPL_beam(:,idx_angle),'o')
xlabel('f (Hz)') 
ylabel('SPL (dB)')
legend('rotor st','rotor unst','beam')




% save('sound_harmonics_model_Jakub_data.mat','BPF','p_f_rotor_unst', ...
%     'p_ff_rotor_unst','SPL_rotor_unst','p_f_rotor_st','p_ff_rotor_st','SPL_rotor_st','p_f_beam','p_ff_beam','SPL_beam')



