% Random walk residence time simulation on a line and comparison with Holzer & Hall (2000) theory
% TWNH Mar '16, Feb '23, Jul '23, Jan '24, Jun '24

%% Housekeeping
clearvars ;
close all ;
more off ;
clc
fprintf(1,'\n Random walk residence time simulation on a line and comparison with Holzer & Hall (2000) theory.\n twnh Mar ''16, Feb ''23, Jul ''23, Jan ''24, Jun ''24\n\n')
fprintf(1,' Residence time means: At a given moment, in equilibrium, of all the particles in the domain, what is the distribution over ages when they have all eventually died?\n\n')

%% Parameters
no_particles = 2^10 ;           % Number of particles released each timestep
exp_bias     = 0.00 ;           % Explicit bias in each step
delta        = 1.0 ;            % Step length
Delta_t      = 1.0 ;            % Time step
no_steps_1   = 8e4 ;            % Number of steps to reach equilibrium
no_steps_2   = 2e5 ;            % Number of steps to allow particles to return to the origin
wall_L       = 2e8 ;            % x-position of reflecting wall
scale_H      = 32 ;             % Exponential density decay scale along the line. This allows to simluate the Holzer & Hall (2000) Appendix B case.
rand_FLAG    = 'discrete' ;

% Derived parameters
k_theory     = delta^2/(24*Delta_t) ;
u_theory     = exp_bias*delta/Delta_t - k_theory/scale_H ;
bias         = u_theory*Delta_t/delta ;

fprintf(1,' Particle random walk on a line with steps of length [%g], timestep [%g], and bias of [%g].\n',delta,Delta_t,bias) ;
fprintf(1,' Particles start at the origin and take steps drawn from a [%s] pdf.\n',rand_FLAG)
fprintf(1,' Particles returning to the origin are removed.\n')
fprintf(1,' Particles reflect off a wall at x = [%g].\n\n',wall_L) ;
fprintf(1,' [%d] particles are released at the origin every timestep until a steady state is reached.\n',no_particles) ;
fprintf(1,' Then particles step until all have returned to the origin.\n The residence time of the steady state population of particles is computed.\n\n') ;
fprintf(1,' Two absorbing barriers guarantees that all particles will eventually be absorbed, see Cox & Miller (1977, p. 29).\n')
fprintf(1,' Therefore, a reflecting plus an absorbing barrier should also guarantee that all particles will be absorbed (although it might take a long time!).\n')
fprintf(1,' For full theory, see Weesakul (1961), although it looks unwieldy.\n\n')

%% Start computing.
fprintf(1,' Stepping particles to establish equilibrium....\n') ;
tic
new_particles  = zeros(no_particles,2) ;     % First column is position, second column is age. Both start at zero.
particles      = [] ;
for tt = 1:no_steps_1
    particles = [particles;new_particles] ;

    % Make random steps and increase age of particles.
    particles = step_and_age_particles(particles,delta,Delta_t,bias,rand_FLAG) ;

    % Find particles reaching far wall, and reflect them back towards the
    % origin.
    wall_inds              = find(particles(:,1) >= wall_L) ;
    particles(wall_inds,1) = 2*wall_L - particles(wall_inds,1) ;

    % Find particles returning to the source point and remove them.
    fld_inds              = find(particles(:,1) <= 0.0) ;
    particles(fld_inds,:) = [] ;
    particles_left        = size(particles,1) ;

    msg = sprintf(' Step [%6.0d] with [%7.0d] particles left: Mean/std/max particle age = [%8.2f]/[%8.2f]/[%6.0d] steps. Mean/std/max particle distance = [%8.2f]/[%8.2f]/[%8.2f]',tt,particles_left,mean(particles(:,2)),std(particles(:,2)),max(particles(:,2)),mean(particles(:,1)),std(particles(:,1)),max(particles(:,1))) ;
    % dispProgressMsg(msg,false)
    dispProgressMsg(msg,false)
    if(max(particles(:,2)) < tt)
        fprintf(1,' Last particle from first release has died.')
        break
    end
end % tt
particles_in_domain = particles ;
dispProgressMsg('test',true) ;
fprintf(1,' Done in [%g]s and [%d] steps.\n\n',toc,tt)

%% Plot statistical equilibrium distribution
figure
subplot(2,1,1)
histogram(particles(:,1))
hold on
yyaxis right
x_vals = 0:1:max(particles(:,1)) ;
y_vals = erfc(x_vals./(6*sqrt(k_theory*tt))) ;              % This is a guessed functional form. Need to justify theoretically.
y_vals2 = exp((u_theory/k_theory).*x_vals) ;
plot(x_vals,y_vals,'r-','LineWidth',2)
hold on
plot(x_vals,y_vals2,'b-','LineWidth',2)
legend('particles','theory (guess): erfc','theory: exp')
title('Spatial distribution of particles at equilibrium')
xlabel('Distance')

%% Display output
% Bin the age of the particles_in_domain
bin_edges  = 1:max(particles_in_domain(:,2)) ;
calR       = histcounts(particles_in_domain(:,2),bin_edges,'Normalization','probability') ;
calR_times = (bin_edges(2:end) + bin_edges(1:end-1))./2 ;

% Add theory. See Holzer & Hall (2000) Appendix B and Julia theory
% notebooks (testing_HH00_theory_v0_dimensional.ipynb).
% The theory assumes no explicit advection and exponential
% density decay. The exponential density decay amounts to a uniform
% advection, however.
calR_theory = sqrt(k_theory./(4*pi*(scale_H^2).*calR_times)).*exp(-k_theory.*calR_times./(4*(scale_H^2))) ;
% Normalize, but trapz(calR_times,calR_theory) is already close to one.
calR_theory = calR_theory./sum(calR_theory) ;     

% Add power law fits
fprintf(1,' Age distribution at statistical equilibrium:\n')
fit_inds = find(calR_times(calR_times < 5e2)) ;
fit_power_law(calR_times(fit_inds),calR(fit_inds)) ;

% Plot
figure
set(gcf,'PaperPosition', [0 0 7 3]);
set(gcf,'Papersize',[7 3]) ;
loglog(calR_times,calR,'ko') ;
hold on
slope_theory = -1/2 ;
fit_theory   = log10(calR_theory(1)) + log10(calR_times(fit_inds)).*slope_theory ;
loglog(calR_times(fit_inds),10.^fit_theory,'r-','linewidth',2) ;

grid on
ylabel('Frequency')
xlabel('Age [steps]')
set(gca,'XLim',[1e0 1e4])
tmp_txt = sprintf('%3.1f power law',slope_theory) ;
legend('Simulation',tmp_txt,'interpreter','latex')
title('Age distribution at statistical equilibrium','Interpreter','latex') 

%% Continue stepping until particles have left the domain:
fprintf(1,' Stepping [%d] particles until they die....\n',size(particles_in_domain,1)) ;
tic
dead_particles = [] ;
for tt = 1:no_steps_2
    % Make random steps and increase age of particles.
    particles = step_and_age_particles(particles,delta,Delta_t,bias,rand_FLAG) ;

    % Find particles reaching far wall, and reflect them back towards the
    % origin.
    wall_inds              = find(particles(:,1) >= wall_L) ;
    particles(wall_inds,1) = 2*wall_L - particles(wall_inds,1) ;

    % Find particles returning to the source point, record them, and remove them.
    fld_inds              = find(particles(:,1) <= 0.0) ;
    dead_particles        = [dead_particles; particles(fld_inds,:)] ;
    particles(fld_inds,:) = [] ;
    particles_left        = numel(particles(:,1)) ;
    msg = sprintf(' Step [%6.0d] with [%7.0d] particles left. Mean/std/max particle age = [%8.2f]/[%8.2f]/[%6.0d] steps. Mean/std/max particle distance = [%8.2f]/[%8.2f]/[%8.2f]',tt,particles_left,mean(particles(:,2)),std(particles(:,2)),max(particles(:,2)),mean(particles(:,1)),std(particles(:,1)),max(particles(:,1))) ;
    dispProgressMsg(msg,false)
    if(particles_left < 1)
        break
    end
end % tt
fprintf(1,' Done in [%g]s and [%d] steps.\n\n',toc,tt)

%% Make residence time plot
% Bin the age at death of the dead_particles
bin_edges  = 1:max(dead_particles(:,2)) ;
calR       = histcounts(dead_particles(:,2),bin_edges,'Normalization','probability') ;
calR_times = (bin_edges(2:end) + bin_edges(1:end-1))./2 ;

% Add theory. See Holzer & Hall (2000) Appendix B and Julia theory
% notebooks (testing_HH00_theory_v0_dimensional.ipynb).
% The theory assumes no explicit advection and exponential
% density decay. The exponential density decay amounts to a uniform
% advection, however.
calR_theory = sqrt(k_theory./(4*pi*(scale_H^2).*calR_times)).*exp(-k_theory.*calR_times./(4*(scale_H^2))) ;
% Normalize, but trapz(calR_times,calR_theory) is already close to one.
calR_theory = calR_theory./sum(calR_theory) ;     

% Add power law fits
fprintf(1,' R, residence-time distribution:\n')
fit_inds = find(calR_times(calR_times < 5e2)) ;
calR_fit = fit_power_law(calR_times(fit_inds),calR(fit_inds)) ;

% Plot
figure
set(gcf,'PaperPosition', [0 0 7 3]);
set(gcf,'Papersize',[7 3]) ;
loglog(calR_times,calR,'ko') ;
hold on
loglog(calR_times(fit_inds),calR_fit(fit_inds),'r-','linewidth',2) ;
slope_theory = -1/2 ;
fit_theory   = log10(calR_theory(1)) + log10(calR_times(fit_inds)).*slope_theory ;
loglog(calR_times(fit_inds),10.^fit_theory,'b-','linewidth',2) ;

grid on
ylabel('Frequency')
xlabel('Residence time [steps]')
set(gca,'XLim',[1e0 1e4])
tmp_txt = sprintf('%3.1f power law',slope_theory) ;
legend('Simulation',tmp_txt,'Fit','interpreter','latex')
title('${\mathcal R}$, Residence time distribution','Interpreter','latex') 

%% Function definitions
function fit = fit_power_law(x,y)

inds = y > 0 ;
coeffs = polyfit(log10(x(inds)),log10(y(inds)),1) ;
fit = 10.^polyval(coeffs,log10(x)) ;
fprintf(1,' Coefficients of power law fit: pre-factor = [%4.3g], exponent = [%4.3g].\n\n',coeffs(2),coeffs(1))

end

function particles = step_and_age_particles(particles,delta,Delta_t,bias,rand_FLAG)

particles_left = size(particles,1) ;
random_nos = rand(particles_left,1) - 0.5 + bias ;
switch rand_FLAG
    case 'uniform'    % Uniform step pdf.
        steps = delta.*random_nos ;
    case 'discrete'   % Discrete pdf (steps right or left 1 space).
        steps = delta.*(random_nos >= 0) - (random_nos < 0) ;
end % switch
particles(:,1) = particles(:,1) + steps ;       % Move particles
particles(:,2) = particles(:,2) + Delta_t ;     % Age  particles

end

function dispProgressMsg(msg,resetFLAG)
% From https://www.mathworks.com/matlabcentral/answers/292577-how-to-print-text-on-the-same-line#answer_747418

ASCII_BKSP_CHAR = 8;
persistent prevMsgLen ;

if(resetFLAG)
    prevMsgLen = 0 ;
    return
else
    if isempty(prevMsgLen)
        prevMsgLen = 0;
    end

    disp([ char(repmat(ASCII_BKSP_CHAR,1,prevMsgLen)) msg]);

    prevMsgLen = numel(msg)+1;
end % if

end