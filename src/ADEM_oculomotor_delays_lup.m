% Eye movements with sensory-motor delays under active inference:
%__________________________________________________________________________
% This demo illustrates tracking eye movements. The focus
% here is on oculomotor delays and their compensation in generalised
% coordinates of motion. This is illustrates using a "sweep" of motion to
% examine the effects of motor delays, sensory delays and their interaction
% under active inference. We then move on to smooth pursuit of sine
% wave motion, in which the trajectory entrains following under compensated
% dynamics. This entrainment can be destroyed by rectifying the sine wave
% which calls for a more realistic (hierarchical) model of motion that
% registers its phase and anticipates the next onset of motion (cf
% movement behind an occluder). These simulations depend delicately on
% delays and precisions (gains) that are chosen to make pursuit initiation
% under uncertainty relatively slow. The dependency on visual uncertainty
% (contrast) is illustrated by changing the precision of the generalised
% motion of the sensory target.
%__________________________________________________________________________
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: ADEM_occulomotor_delays.m 4663 2012-02-27 11:56:23Z karl $
%close all; 
clear all
spm_jobman('initcfg');
system('rm -f ../../../lup_figures_dev/*');

switch_section1 = 1;
switch_section2 = 1;
% definition of time
%-------------------
N       = 96;                                 % length of data sequence
%N       = 64;                                 % length of data sequence
dt      = 16;                                 % time step (ms)
% TODO control to see if dt has some importance at all
% HACK : dt should be the same in function spm_DEM_qU(qU,pU)
%N       = 128;                                 % length of data sequence
%dt      = 8;                                 % time step (ms)
%N       = 256;                                 % length of data sequence
%dt      = 4;                                 % time step (ms)
%N       = 512;                                 % length of data sequence
%dt      = 2;                                 % time step (ms)

t_m     = 16;                               % time constant for sensory (ms)
t_s     = 32;                               % time constant for sensory (ms)
t_a     = 64; % 100 would be more realistic                               % time constant for action (ms)
t_o     = 1024;                            % damping time constant (ms)

M1V = exp(4);% error precision
M1W = exp(4);% error precision
M2V = exp(0);
G1V = exp(16);% error precision (errors)
G1W = exp(16);% error precision (motion)
G2V  = exp(16);
G1U = exp(2); % motor gain for the sweep
% motor gain after fig 4 (sinusoids)
G1U_bis = exp(4);
% only in the anticipatory model
M2V_bis = exp(16);%
M2W = exp(4);% 
M3V  = exp(4);


Ta=32;
Ty=32;

t_sweep = 80;
t_ramp = 1024;

t_period = 512;
%t_period = 512;
U       =  2*pi/(t_period/dt); % 2*pi/32; %

t_start = 256;

t       = [1:N]*dt  - t_start;                           % PST

amplitude = 1;


% trajectories
%-------------
% sigmoid (potentially not causal as it starts before t_start...)
% C_sweep    = [amplitude; 0]*spm_phi(((1:N) - t_start/dt)/(t_sweep/dt/2));    % simple motion (sweep)
% pure ramp
C_ramp    = [amplitude; 0]*cumsum( (1:N > t_start/dt).*(1:N < ((t_start+t_ramp)/dt))/(t_ramp/dt));
%C_ramp    = C_ramp - 0.02 * [amplitude; 0]*( (1:N > t_start/dt).*(1:N < ((t_start+t_ramp)/dt)));
% simple motion (ramp)
%C_sweep    = C_sweep.*([amplitude; 0]*(1:N < (t_start/dt+t_sweep/dt)));
%% ramp returns to zero
%C_sweep    = C_sweep.*([amplitude; 0]*(1-(1-exp(- (((1:N) - (t_start+t_sweep)/dt)/(30/dt)))).*((1:N) > (t_start/dt+t_sweep/dt))));
% first-order
% C_sweep    = [amplitude; 0]*(1-exp(- (((1:N) - t_start/dt)/(t_sweep/2/dt)).*((1:N) > (t_start/dt))));
% gaussian-like
C_sweep    = [amplitude; 0]*(1-exp(- .5 * (((1:N) - t_start/dt).^2/(t_sweep/dt)^2).*((1:N) > (t_start/dt))));

C_sin       = [amplitude; 0]*((sin(2*pi*((1:N) - t_start/dt)/(t_period/dt))).*((1:N) > (t_start/dt)));

% hidden causes and states
%==========================================================================
% x    - hidden states:
%   x.o(1) - oculomotor angle
%   x.o(2) - oculomotor angle
%   x.x(1) - target location (visual) - extrinsic coordinates (Cartesian)
%   x.x(2) - target location (visual) - extrinsic coordinates (Cartesian)
%
% v    - causal states: force on target
%
% g    - sensations:
%   g(1) - oculomotor angle (proprioception)
%   g(2) - oculomotor angle (proprioception)
%   g(3) - target location (visual) - intrinsic coordinates (polar)
%   g(4) - target location (visual) - intrinsic coordinates (polar)
%--------------------------------------------------------------------------
 
% parameters mapping from (unstable) point attractors to visual space
%--------------------------------------------------------------------------
x.o = [0;0];                                  % oculomotor angle
x.x = [0;0];                                  % target location
 
% Recognition model
%==========================================================================
M(1).E.s = 1/2; %1/2;                               % smoothness - s.d. (bins)
M(1).E.n = 4; % 4;                                % order of
M(1).E.d = 1; % 1;                                % generalised motion
 
% level 1: Displacement dynamics and mapping to sensory/proprioception
%--------------------------------------------------------------------------
M(1).f   = ['[(x.x - x.o)/', num2str(t_s/dt), '; (v - x.x)/', num2str(t_m/dt), ']'];        % extrinsic coordinates
%M(1).f   = '[(x.x - x.o)/2; v - x.x]';        % extrinsic coordinates
M(1).g   = '[x.o; x.x - x.o]';                % intrinsic coordinate
%M(1).g  = '[x.o; tanh(1*(x.x - x.o))]';                 % intrinsic coordinates
%M(1).g  = '[exp((x(1) - 1)*4); 0]';          
M(1).x   = x;                                 % hidden states
M(1).V   = M1V;                            % error precision
M(1).W   = M1W;                            % error precision

% TODO: implement the fact that precision decreases when the object is farther away

% level 2:
%--------------------------------------------------------------------------
%M(2).v  = '[v/4]';                             % inputs
M(2).v  = [0; 0];                             % inputs
M(2).V  = M2V;
 

% generative model
%==========================================================================
 
% first level
%--------------------------------------------------------------------------
G(1).f  = ['[a/', num2str(t_a/dt), ' - x.o/', num2str(t_o/dt), '; (v - x.x)/', num2str(t_m/dt), ']'];          % extrinsic coordinates
%G(1).f  = '[a/4 - x.o/32; v - x.x]';          % extrinsic coordinates
G(1).g  = '[x.o; x.x - x.o]';                 % intrinsic coordinates
%G(1).g  = '[x.o; tanh(1.*(x.x - x.o))]';                 % intrinsic coordinates
G(1).x  = x;                                  % hidden states
G(1).V  = G1V;                            % error precision (errors)
G(1).W  = G1W;                            % error precision (motion)
G(1).U  = [1 1 0 0]*G1U;         % error precision (action)
 
% second level
%--------------------------------------------------------------------------
G(2).v  = [0; 0];                             % exogenous forces
G(2).a  = [0; 0];                             % action forces
G(2).V  = G2V;

DEM.G   = G;
DEM.M   = M;


% Solve or integrate
%==========================================================================
% Note : the 4 diferrent sweep stimulations produce similar results

DEM.C   = C_sweep;
%DEM.C   = C_sin;

DEM.U   = zeros(2,N);

% std sweep = t_sweep / units of angular dispacement = 4deg
disp(['peak speed = ', num2str(4/sqrt(2*pi)/(t_sweep/1000))])


if switch_section1,
% note : In terms of the code; Ta and Ty represent the difference between
% the true delay and the expected delay - it would be redundant to specify 
% both separately because of the nature of the temporal delay operators 
% (only the difference is needed).
    
% compensated delays
%--------------------------------------------------------------------------
DEM.M(1).Ta = 0;
DEM.M(1).Ty = 0;
DEM         = spm_ADEM(DEM);
 
 
 
% Compare latencies graphically
%--------------------------------------------------------------------------
PF   = spm_platform('fonts');
spm_figure('GetWin','Figure 3', 'FontName', spm_platform('Font','Times'), 'Visible', 'on'); clf

spm_DEM_qU(DEM.qU,DEM.pU)
h1 = subplot(2,2,1); title({'Pursuit initiation:', 'prediction and error'},'FontSize',14);
h2 = subplot(2,2,2); title({'Ramp motion', 'hidden states'},'FontSize',14);
h3 = subplot(2,2,3); 
h4 = subplot(2,2,4); %hold on; plot(t,DEM.pU.v{1},'g',t,DEM.qU.v{1},'--g')
spm_axis([h1, h2, h3, h4],[min(t), max(t), -0.19, 1.19])%AXIS([XMIN XMAX YMIN YMAX])

system('rm -f ../../../lup_figures_dev/figure3.ps')
spm_figure('Print', 'Graphics', '../../../lup_figures_dev/figure3.ps');


spm_figure('GetWin','Figure 3 bis', 'FontName', spm_platform('Font','Times'), 'Visible', 'on'); clf

% TODO step ramp stimulus from Rashbass/ Barnes

DEM.M(1).Ta = 0;
DEM.M(1).Ty = 0;

amp = (0.25:0.25:1);
for i = 1:length(amp)
    DEM.C      = amp(i) * C_ramp; %C_sweep;
    DEM.U   = zeros(2,N);
    DEM        = spm_ADEM(DEM);
    spm_figure('GetWin','Figure 3 bis');
    subplot(2,1,1)
    plot(t,DEM.qU.v{1}(1:2,:), 'Color', [1-amp(i),0,0]), hold on
    plot(t,DEM.C, '--', 'Color', [1-amp(i),0,0]), hold on
    subplot(2,1,2)
    plot(t(1:(N-1)),conv2(full(DEM.qU.v{1}(1:2,:)),[1/2*dt,-1/2*dt],'valid'), 'Color', [1-amp(i),0,0]), hold on
    plot(t(1:(N-1)),conv2(DEM.C,[1/2*dt,-1/2*dt],'valid'), '--', 'Color', [1-amp(i),0,0]), hold on
end
 
% Graphics
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 3 bis');
h1 = subplot(2,1,1), hold off
xlabel('Time (ms)','FontSize',12), spm_axis tight
ylabel('prediction and error','FontSize',12), spm_axis tight
 
h2 = subplot(2,1,2)
xlabel('Time (ms)','FontSize',12), spm_axis tight
ylabel('prediction and error (speed)','FontSize',12)

spm_axis(h1,[min(t), max(t), -.19, 1.19])%AXIS([XMIN XMAX YMIN YMAX])
spm_axis(h2,[min(t), max(t), -.019, .19])%AXIS([XMIN XMAX YMIN YMAX])

system('rm -f  ../../../lup_figures_dev/figure3bis.ps');
spm_figure('Print', 'Graphics', '../../../lup_figures_dev/figure3bis.ps');

% uncompensated sensory delays
%--------------------------------------------------------------------------
DEM.C   = C_sweep;
DEM.M(1).Ta = 0;
DEM.M(1).Ty = 0;
DEM        = spm_ADEM(DEM);
DEM.M(1).Ta = 0;
DEM.M(1).Ty = Ty/dt;
DEMy        = spm_ADEM(DEM);
 
% uncompensated motor delays
%--------------------------------------------------------------------------
DEM.M(1).Ta = Ta/dt;
DEM.M(1).Ty = 0;
DEMa        = spm_ADEM(DEM);
 
% uncompensated delays
%--------------------------------------------------------------------------
DEM.M(1).Ta = Ta/dt;
DEM.M(1).Ty = Ty/dt;
DEMay       = spm_ADEM(DEM);

spm_figure('GetWin','Figure 4', 'FontName', spm_platform('Font','Times')); clf
h1 = subplot(3,2,1); 
%disp(DEM.pU.v{1})
plot(t,DEM.pU.v{1},'b',t,DEM.qU.v{1},'--b'), hold on
subplot(3,2,1); plot(t,DEMy.pU.v{1},'r',t,DEMy.qU.v{1},'--r'), hold on
p = plot(0,0,'k',0,0,'k--',0,0,'b',0,0,'r'); plot(t,C_sweep,'k--')
legend(p,'sensory','estimation','control','delayed','Location','East')
h2 = subplot(3,2,2); plot(t,DEM.qU.a{2},'b'), hold on
subplot(3,2,2); plot(t,DEMy.qU.a{2},'r'), hold on
plot(t,C_sweep,'k--')
p2 = plot(0,0,'b',0,0,'r');
legend(p2,'control','delayed','Location','East')

h3 = subplot(3,2,3); plot(t,DEM.pU.v{1},'b',t,DEM.qU.v{1},'--b'), hold on, plot(t,C_sweep,'k--')
h4 = subplot(3,2,4); plot(t,DEM.qU.a{2},'b'), hold on
subplot(3,2,3); plot(t,DEMa.pU.v{1},'r',t,DEMa.qU.v{1},'--r'), hold on
subplot(3,2,4); 
plot(t,DEMa.qU.a{2},'r'), hold on, 
 
h5 = subplot(3,2,5); plot(t,DEM.pU.v{1},'b',t,DEM.qU.v{1},'--b'), hold on, plot(t,C_sweep,'k--')
h6 = subplot(3,2,6); plot(t,DEM.qU.a{2},'b'), hold on
subplot(3,2,5); plot(t,DEMay.pU.v{1},'r',t,DEMay.qU.v{1},'--r'), hold on
subplot(3,2,6); plot(t,DEMay.qU.a{2},'r'), % hold on, plot(t,C_sweep,'k--')


% labels and titles
%--------------------------------------------------------------------------
subplot(3,2,1); title('Sensory input and predictions','FontSize',14), spm_axis tight
subplot(3,2,2); title('Action','FontSize',14), spm_axis tight
 
subplot(3,2,5); xlabel('Time (ms)','FontSize',12), spm_axis tight
subplot(3,2,6); xlabel('Time (ms)','FontSize',12), spm_axis tight
subplot(3,2,4); spm_axis tight
 
% subplot(3,2,1); ylabel('Sensory delays','FontSize',14), spm_axis tight
% subplot(3,2,3); ylabel('Motor delays','FontSize',14), spm_axis tight
% subplot(3,2,5); ylabel('Sensorimotor delays','FontSize',14), spm_axis tight
subplot(3,2,1); ylabel('Position','FontSize',12)
subplot(3,2,3); ylabel('Position','FontSize',12)
subplot(3,2,5); ylabel('Position','FontSize',12)

spm_axis([h1, h2, h3, h4, h5, h6],[min(t), max(t), -0.19*amplitude, 1.19*amplitude])%AXIS([XMIN XMAX YMIN YMAX])

subplot(3,2,1); text(max(t)*1.18, -0.46, 'Sensory delays','FontSize',14, 'HorizontalAlignment', 'center'), %spm_axis tight
subplot(3,2,3); text(max(t)*1.18, -0.46, 'Motor delays','FontSize',14, 'HorizontalAlignment', 'center'), %spm_axis tight
subplot(3,2,5); text(max(t)*1.18, -0.46, 'Sensory and motor delays','FontSize',14, 'HorizontalAlignment', 'center'), %spm_axis tight
% set(hax,'XTick',[],'YTick',[])
system('rm -f ../../../lup_figures_dev/figure4.ps');
spm_figure('Print', 'Graphics', '../../../lup_figures_dev/figure4.ps');

spm_figure('GetWin','Figure 4 bis', 'FontName', spm_platform('Font','Times')); clf

% oscillations in velocity to step ramp stimulus from Goldreich: The natural delay (typically ~90ms) because of the computation in the visual inputs to drive the pursuit system leads to oscillations (Goldreich et al 1992). Changing the target viewing conditions (e.g. contrast and size of the stimulus) changes the natural delay. Electronically controlled feed-back delay is imposed and the responses for different delays show oscillations with greater amplitudes and lower frequency (tuned to the total delay).

delay = (0.25:0.5:2) * Ty/dt;
for i = 1:length(delay)

    DEM.M(1).Ta = delay(i); % 
    DEM.M(1).Ty = 0; % delay(i);
    DEM.C      = C_ramp; %C_sweep;
    DEM.U      = zeros(2,N);
    DEM        = spm_ADEM(DEM);
    spm_figure('GetWin','Figure 4 bis');
    subplot(2,1,1)
    plot(t,DEM.qU.v{1}(1:2,:), 'Color', [1-amp(i),0,0]), hold on
    plot(t,DEM.C, '--', 'Color', [1-amp(i),0,0]), hold on
    subplot(2,1,2)
    plot(t(1:(N-1)),conv2(full(DEM.qU.v{1}(1:2,:)),[1/2*dt,-1/2*dt],'valid'), 'Color', [1-amp(i),0,0]), hold on
    plot(t(1:(N-1)),conv2(DEM.C,[1/2*dt,-1/2*dt],'valid'), '--', 'Color', [1-amp(i),0,0]), hold on
    
end
 
% Graphics
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 4 bis');
h1 = subplot(2,1,1), hold off
xlabel('Time (ms)','FontSize',12), spm_axis tight
ylabel('prediction and error','FontSize',12), spm_axis tight
 
h2 = subplot(2,1,2)
xlabel('Time (ms)','FontSize',12), spm_axis tight
ylabel('prediction and error (speed)','FontSize',12)

spm_axis(h1,[min(t), max(t), -.19, 1.19])%AXIS([XMIN XMAX YMIN YMAX])
spm_axis(h2,[min(t), max(t), -.019, .19])%AXIS([XMIN XMAX YMIN YMAX])

system('rm -f  ../../../lup_figures_dev/figure3bis.ps')
spm_figure('Print', 'Graphics', '../../../lup_figures_dev/figure4bis.ps');

% Contrast-dependent latencies (compensated for delays)
%==========================================================================
DEM.C   = C_sweep;
DEM.U   = zeros(2,N);

DEM.M(1).Ta = 0;
DEM.M(1).Ty = 0;
%DEM.M(1).Ta = Ta/dt;
%DEM.M(1).Ty = Ty/dt;
DEM.G(1).U  = [1 1 0 0]*G1U;
 
% change prior precision (contrast)
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 5', 'FontName', spm_platform('Font','Times'), 'Visible', 'on'); clf
 
V = (2:8);                         
for i = 1:length(V)
    DEM.M(1).V = diag(exp([4 4 V(i) V(i)])); % error precision
    DEM        = spm_ADEM(DEM);
    pd         = sqrt(sum(DEM.pU.v{1}(3:4,:).^2));
    qd         = sqrt(sum(DEM.qU.v{1}(3:4,:).^2));
    pD(i)      = max(pd);
    qD(i)      = max(qd);
    
    spm_figure('GetWin','Figure 5');
    subplot(2,1,1)
    plot(t,pd,'k',t,qd,':k', 'Color', [(i-1)/(length(V)-1),0,0]), hold on
    
end
 
% Graphics
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 5');
subplot(2,1,1), hold off
xlabel('Time (ms)','FontSize',12), spm_axis tight
ylabel('Spatial lag','FontSize',12), spm_axis tight
title('Distance from target','FontSize',14)
axis square
 
subplot(2,1,2)
%plot(V,pD,'-k',V,pD,'ok',V,qD,':k',V,qD,'xk'), hold on 
plot(V,pD,'-k',V,qD,':k'), hold on 
for i = 1:length(V)
    plot(V(i),pD(i),'o', 'Color', [(i-1)/(length(V)-1),0,0]), hold on 
    plot(V(i),qD(i),'x', 'Color', [(i-1)/(length(V)-1),0,0]), hold on 
end
xlabel('log-Precision','FontSize',12)
ylabel('Spatial lag','FontSize',12)
title('Peak lag','FontSize',14)
axis square
 
system('rm -f ../../../lup_figures_dev/figure5.ps')
spm_figure('Print', 'Graphics', '../../../lup_figures_dev/figure5.ps');
%spm_figure('Close','Figure 5');



end

if switch_section2, 
% Oculomotor latencies (sinusoidal movement)
%==========================================================================
DEM.M(1).Ta = 0;
DEM.M(1).Ty = 0;
%DEM.M(1).Ta = Ta/dt;
%DEM.M(1).Ty = Ty/dt;
DEM.G(1).U  = [1 1 0 0]*G1U_bis;
 
% Pursuit initiation with (first order) generative model and sine wave
%--------------------------------------------------------------------------
DEM.C   = C_sin;
DEM.U   = zeros(2,N);
DEM     = spm_ADEM(DEM);
 
spm_figure('GetWin','Figure 6', 'FontName', spm_platform('Font','Times'));
spm_DEM_qU(DEM.qU,DEM.pU)
h1 = subplot(2,2,1); title({'Pursuit initiation:', 'prediction and error'},'FontSize',14)
h2 = subplot(2,2,2); title({'Periodic motion:', 'hidden states'},'FontSize',14)
h3 = subplot(2,2,3); 
h4 = subplot(2,2,4); %hold on; plot(t,DEM.pU.v{1},'g',t,DEM.qU.v{1},'--g')
spm_axis([h1, h2, h3, h4],[min(t), max(t), -1.19, 1.19])%AXIS([XMIN XMAX YMIN YMAX])

system('rm -f ../../../lup_figures_dev/figure6.ps')
spm_figure('Print', 'Graphics', '../../../lup_figures_dev/figure6.ps');

% smooth pursuit following with (second order) generative model
%--------------------------------------------------------------------------
mx.o = [0;0];                                  % oculomotor angle
mx.x = [0;0];                                  % target location
mx.v = [0;0];                                  % oculomotor velocity
 
% level 1: Displacement dynamics and mapping to sensory/proprioception
%--------------------------------------------------------------------------
M(1).f   = '[x.v; v - x.x; (v - x.o)/4 - x.v/2]';
M(1).g   = '[x.o; x.x - x.o]';                 % intrinsic coordinate
M(1).x   = mx;                                 % hidden states
M(1).V   = M1V;                             % error precision
M(1).W   = M1W;                             % error precision
 
% level 2:
%--------------------------------------------------------------------------
M(2).v  = [0; 0];                              % inputs
M(2).V  = M2V;
 
 
% generative model
%==========================================================================
gx.o = [0;0];                                  % oculomotor angle
gx.x = [0;0];                                  % target location
gx.v = [0;0];                                  % oculomotor velocity
 
% first level
%--------------------------------------------------------------------------
G(1).f  = ['[x.v; v - x.x; a/', num2str(t_a/dt), ' - x.v]'];
G(1).g  = '[x.o; x.x - x.o]';                 % intrinsic coordinates
G(1).x  = gx;                                 % hidden states
G(1).V  = G1V;                            % error precision (errors)
G(1).W  = G1W;                            % error precision (motion)
G(1).U  = [1 1 0 0]*G1U_bis;                   % motor gain
 
% second level
%--------------------------------------------------------------------------
G(2).v  = [0; 0];                             % exogenous forces
G(2).a  = [0; 0];                             % action forces
G(2).V  = G2V;
 
% Sine wave cause
%--------------------------------------------------------------------------
DEM.M   = M;
DEM.G   = G;
DEM.C   = C_sin;
DEM.U   = zeros(2,N);
DEM     = spm_ADEM(DEM);
 
spm_figure('GetWin','Figure 7', 'FontName', spm_platform('Font','Times'));
spm_DEM_qU(DEM.qU,DEM.pU)
h1 = subplot(2,2,1); title({'Smooth pursuit:', 'prediction and error'},'FontSize',14)
h2 = subplot(2,2,2); title({'Periodic motion:', 'hidden states'},'FontSize',14)
h3 = subplot(2,2,3); 
h4 = subplot(2,2,4); %hold on; plot(t,DEM.pU.v{1},'g',t,DEM.qU.v{1},'--g')
spm_axis([h1, h2, h3, h4],[min(t), max(t), -1.19, 1.19])%AXIS([XMIN XMAX YMIN YMAX])

system('rm -f ../../../lup_figures_dev/figure7.ps')
spm_figure('Print', 'Graphics', '../../../lup_figures_dev/figure7.ps');

% % Control with the sweep cause
% %--------------------------------------------------------------------------
% DEM.M   = M;
% DEM.G   = G;
% DEM.C   = C_sweep;
% DEM.U   = zeros(2,N);
% DEM     = spm_ADEM(DEM);
%  
% spm_figure('GetWin','Figure 7 bis');
% spm_DEM_qU(DEM.qU,DEM.pU)
% h1 = subplot(2,2,1); title({'Smooth pursuit:', 'prediction and error'},'FontSize',16)
% h2 = subplot(2,2,2); title({'Sweep motion:', 'hidden states'},'FontSize',16)
% h3 = subplot(2,2,3); 
% h4 = subplot(2,2,4); hold on; plot(t,DEM.pU.v{1},'g',t,DEM.qU.v{1},'--g')
% spm_axis([h1, h2, h3, h4],[min(t), max(t), -.19, 1.19])%AXIS([XMIN XMAX YMIN YMAX])
% 
% system('rm -f ../../../lup_figures_dev/figure7bis.ps')
% spm_figure('Print', 'Graphics',
% '../../../lup_figures_dev/figure7bis.ps');

  
% rectified sinusoidal movement
%--------------------------------------------------------------------------
%DEM.C      = exp((C_sin - 1)*4);
DEM.C      = C_sin.*(C_sin>0);
DEM        = spm_ADEM(DEM);
 
% Graphics
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 8', 'FontName', spm_platform('Font','Times'));
spm_DEM_qU(DEM.qU,DEM.pU)
h1 = subplot(2,2,1); title({'Smooth pursuit:', 'prediction and error'},'FontSize',14)
h2 = subplot(2,2,2); title({'Aperiodic motion:', 'hidden states'},'FontSize',14)
h3 = subplot(2,2,3); 
h4 = subplot(2,2,4); %hold on; plot(t,DEM.pU.v{1},'g',t,DEM.qU.v{1},'--g')
spm_axis([h1, h2, h3, h4],[min(t), max(t), -1.19, 1.19])%AXIS([XMIN XMAX YMIN YMAX])

system('rm -f ../../../lup_figures_dev/figure8.ps')
spm_figure('Print', 'Graphics', '../../../lup_figures_dev/figure8.ps');
% 
% % Smooth pursuit with a hierarchical generative model
% %==========================================================================
% % Endow the model with internal dynamics (a simple integrator) so that is
% % recognises and remembers the smooth trajectory.
%  
%  
% % level 1: Displacement dynamics and mapping to sensory/proprioception
% %--------------------------------------------------------------------------
% 
% M(1).f   = '[x.v; v - x.x; (v - x.o)/4 - x.v/2]';
% M(1).g   = '[x.o; x.x - x.o]';                 % intrinsic coordinate
% M(1).x   = mx;                                 % hidden states
% M(1).V   = M1V;                             % error precision
% M(1).W   = M1W;                             % error precision
%  
% % level 2: With hidden (memory) states
% %--------------------------------------------------------------------------
% 
% M(2).f  = '[1 0; 0 1]*x*v'; 
% M(2).g  = '[x(1); 0]';          
% M(2).x  = [0; 0];                              % hidden states
% M(2).v  = [1;  0];                           % hidden causes
% M(2).V   = M2V;                             % error precision
% M(2).W   = M2W;                             % error precision
%  
% % level 3: Encoding frequency of memory states (U)
% %--------------------------------------------------------------------------
% U       = 1;
% M(3).v  = U;
% M(3).V  = M3V;
% 
%  
% 
% % Control with the sweep cause
% %--------------------------------------------------------------------------
% DEM.M   = M;
% DEM.G   = G;
% DEM.C   = C_sweep;
% DEM.U   = zeros(1,N) + U;
% DEM     = spm_ADEM(DEM);
%  
% spm_figure('GetWin','Figure 8 bis');
% spm_DEM_qU(DEM.qU,DEM.pU)
% %h1 = subplot(2,2,1); title({'Smooth pursuit:', 'prediction and error'},'FontSize',16)
% %h2 = subplot(2,2,2); title({'Sweep motion:', 'hidden states'},'FontSize',16)
% h1 = subplot(3,2,1); 
% h2 = subplot(3,2,2); 
% h3 = subplot(3,2,3);
% h4 = subplot(3,2,4); 
% h5 = subplot(3,2,5); 
% h6 = subplot(3,2,6); 
% spm_axis([h1, h2, h3, h4, h6],[min(t), max(t), -.19, 1.19])%AXIS([XMIN XMAX YMIN YMAX])
% spm_axis([h5],[min(t), max(t), 0, 1.1])%AXIS([XMIN XMAX YMIN YMAX])
% %hold on; plot(t,DEM.pU.v{1},'g',t,DEM.qU.v{1},'--g')
% 
% system('rm -f ../../../lup_figures_dev/figure8bis.ps')
% spm_figure('Print', 'Graphics', '../../../lup_figures_dev/figure8bis.ps')
% %spm_figure('Print', 'Graphics', '../../../lup_figures_dev/figures_lup.ps')

 
% Smooth pursuit with a hierarchical generative model
%==========================================================================
% Endow the model with internal dynamics (a simple oscillator) so that is
% recognises and remembers the trajectory to anticipate jumps in rectified
% sinusoidal motion.
 
 
% level 1: Displacement dynamics and mapping to sensory/proprioception
%--------------------------------------------------------------------------
%M(1).Ta = Ta/dt;
%M(1).Ty = Ty/dt;

M(1).f   = '[x.v; v - x.x; (v - x.o)/4 - x.v/2]';
M(1).g   = '[x.o; x.x - x.o]';                 % intrinsic coordinate
M(1).x   = mx;                                 % hidden states
M(1).V   = M1V;                             % error precision
M(1).W   = M1W;                             % error precision
 
% level 2: With hidden (memory) states
%--------------------------------------------------------------------------
M(2).f  = '[0 1; -1 0]*x*v'; 
%M(2).f  = '[0 -1; 1 0]*x*v'; 
%M(2).g  = '[exp((x(1) - 1)*4); 0]';          
M(2).g  = '[x(1).*(x(1)>0); 0]';          
%M(2).x  = [0; 0];                              % hidden states
M(2).x  = [0; -1]; % giving the phase and amplitude                              % hidden states
M(2).V   = M2V_bis;                             % error precision
M(2).W   = M2W;                             % error precision
 
% level 3: Encoding frequency of memory states (U)
%--------------------------------------------------------------------------
M(3).v  = U;
M(3).V  = M3V;
 
% give the model prior beliefs about angular motion (but not its phase)
%--------------------------------------------------------------------------
DEM.G   = G;
DEM.M   = M;
%DEM.C      = exp((C_sin - 1)*4);
DEM.C      = C_sin.*(C_sin>0);
DEM.U   = zeros(1,N) + U;
DEM     = spm_ADEM(DEM);
 
% Graphics
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 9', 'FontName', spm_platform('Font','Times'));
spm_DEM_qU(DEM.qU,DEM.pU)
subplot(3,2,1), title({'Anticipatory model:', 'prediction and error'},'FontSize',14)
subplot(3,2,2), title({'Aperiodic motion:', 'hidden states'},'FontSize',14)
h1 = subplot(3,2,1); 
h2 = subplot(3,2,2); 
h3 = subplot(3,2,3);
h4 = subplot(3,2,4); 
h5 = subplot(3,2,5); 
h6 = subplot(3,2,6); 
spm_axis([h1, h2, h3, h4, h6],[min(t), max(t), -1.19, 1.19])%AXIS([XMIN XMAX YMIN YMAX])
spm_axis([h5],[min(t), max(t), 0, 1.])%AXIS([XMIN XMAX YMIN YMAX])

system('rm -f ../../../lup_figures_dev/figure9.ps');
spm_figure('Print', 'Graphics', '../../../lup_figures_dev/figure9.ps');

% % Control with the sweep cause
% %--------------------------------------------------------------------------
% % level 3: Encoding frequency of memory states (U)
% %--------------------------------------------------------------------------
% U       = 0; % 2*pi/32;
% M(3).v  = U;
% M(3).V  = M3V;
% 
% DEM.M   = M;
% DEM.G   = G;
% DEM.C   = C_sweep;
% DEM.U   = zeros(1,N) + U;
% DEM     = spm_ADEM(DEM);
% 
% spm_figure('GetWin','Figure 9 bis');
% spm_DEM_qU(DEM.qU,DEM.pU)
% subplot(3,2,1), title({'Anticipatory model:', 'prediction and
% error'},'FontSize',14)
% subplot(3,2,2), title({'Aperiodic motion:', 'hidden states'},'FontSize',14)
% h1 = subplot(3,2,1);
% h2 = subplot(3,2,2);
% h3 = subplot(3,2,3);
% h4 = subplot(3,2,4);
% h5 = subplot(3,2,5);
% h6 = subplot(3,2,6);
% spm_axis([h1, h2, h3, h4, h6],[min(t), max(t), -1.19, 1.19])%AXIS([XMIN XMAX YMIN YMAX])
% spm_axis([h5],[min(t), max(t), 0, 1.])%AXIS([XMIN XMAX YMIN YMAX])
% 
% system('rm -f ../../../lup_figures_dev/lup_figures_dev/figure9bis.ps')
% spm_figure('Print', 'Graphics', '../../../lup_figures_dev/lup_figures_dev/figure9bis.ps')


% TODO replicate Anna + Bennett, S. J., Orban de Xivry, J.-J., Lefèvre, P.,
% & Barnes, G. R. (2010). Oculomotor prediction of accelerative target motion during occlusion: long-term and short-term effects. Experimental brain research. Experimentelle Hirnforschung. Expérimentation cérébrale, 204(4), 493?504. doi:10.1007/s00221-010-2313-4


 
% create movie in extrinsic and intrinsic coordinates
%==========================================================================
% TODO : show eye image in movie (as for mountain car)
%spm_figure('GetWin','Figure 10'); clf
%spm_dem_pursuit_movie(DEM,0)

end

system('sh ../../../pdf_cnv.sh ../../../*.ps');
disp('Goodbye')
