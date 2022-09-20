%% Codes for QLM
% Notions
%   M - matter sites, and G - gauge sites
% 
% State configuration, start and end both with matter sites
%    lattice-site idx 'k':    1 2 3 4 5 6 7 8 9   ...
%     QLM state configs. :  | M G M G M G M G M >
%
% 1. for matter sites, maps between atom cofiguration (AN: atom number)
%    current site index 'j' = (k+1)/2, (1,2,3,...)
%      AN  |  Odd   |  Even
%    ----- | ------ | ------
%      1   |   +e   |   -e
%      0   |    0   |    0
%
% 2. for gauge sites, lattice sites index/2
%    current site index 'j' = k/2, (1,2,3,...)
%    and '<-' = -0.5, '->' = 0.5 
%       AN  |   Odd  |  Even
%     ----- | ------ | ------
%        2  |   <-   |   ->
%        0  |   ->   |   <-
%

close all
clc

addpath(genpath('libs'))


%% 
init_atom_configS = [0,2,0,0,0,2,0,0,1,0,0,0,0,2,0];
% init_atom_configS = [0,2,0,0,0,2,0,0,0,2,0,0,0,2,0];
% init_atom_configS = [1,0,1,0,1,0,1,0,0,0,1,0,1,0,1];

L = length(init_atom_configS);
if ~mod(L,2)
    error('Error! Initial sites number should be odd!')
end

init_atom_matterS = init_atom_configS(1:2:end);
init_atom_gaugeS = init_atom_configS(2:2:end);

n_matter = length(init_atom_matterS);
n_gauge = length(init_atom_gaugeS);


init_state_matterS = init_atom_matterS .* ((-1).^((1:n_matter)-1));
init_state_gaugeS = (init_atom_gaugeS-1)/2 .* ((-1).^(1:n_gauge));

init_Gl = zeros(1,n_matter);
init_Gl(1) = 0;
init_Gl(end) = 0;
init_Gl(2:end-1) = init_state_gaugeS(2:end) ...
    - init_state_gaugeS(1:end-1) - init_state_matterS(2:end-1);

gauge_edge = zeros(1,2);
gauge_edge(1) = - init_Gl(1) - init_state_matterS(1) ...
    + init_state_gaugeS(1);
gauge_edge(2) = init_Gl(end) + init_state_matterS(end) ...
    + init_state_gaugeS(end);


%% basis
basis_qlm = boson_basis_1dqlm(L,init_Gl,gauge_edge);
ns = basis_qlm.ns;
fprintf('Total basis number is %d.\n',ns)


%% initial state
init_idxC = [base2dec(strrep(num2str(abs(init_state_matterS)),' ',''),2),...
    base2dec(strrep(num2str(init_state_gaugeS+0.5),' ',''),2)];
state_idxlt = basis_qlm.state_idxlt;

[~,state_idx_init] = ismember(init_idxC,...
                state_idxlt(:,1:2),'rows');

phi_init = zeros(ns,1);
phi_init(state_idx_init) = 1;


%% Hamiltonian parameters
w0 = 25 * 2 * pi;
w = ones(1,n_gauge) * w0;    % gauge-matter coupling
m = -6 * abs(w0);     % matter field
h = 0.0 * 2 * pi;     % background-electronic field 
% h = (rand(1,n_gauge)-0.5) * 20.0 * 2 * pi; % background-electronic field 
% h = ones(1,n_gauge) * 4.0 * 2 * pi; % background-electronic field 

% w = 1;
% m = abs(4 * w);
% h = 0.0 * w;


%% generate hamiltonian
tic
ham_elems = hamiltonian_1dqlm_elements_generate(basis_qlm);
toc

%
ham_all = hamiltonian_1dqlm(basis_qlm,ham_elems,w,m,h);
ham = ham_all.ham;
fprintf('Number of nonzero hamiltonian elements: %d.\n',nnz(ham))


%% for ground state check
% [V,D] = eig(full(ham));
% [V,D] = eigs(ham,3,'sr')
% % 
% psi_gd = V(:,1)
% probl = psi_gd .* conj(psi_gd)
% 
% basis_qlm.matter_atom_stateS(7,:)
% basis_qlm.gauge_atom_stateS(7,:)



%% ramp time list
nt = 201;
T = 120 * 1e-03;
% nt = 201;
% T = 3/abs(w)*(2*pi);

tl = linspace(0,T,nt);
dt = tl(2) - tl(1);


%% for two-point correlation statistics
init_matter_occuP = init_atom_matterS;
init_matter_occuP(init_matter_occuP==0) = -1;
init_matter_occuP = [0 0 0 0 0 -1 0 0];


%%
density_gauge_Mt = [];
density_matter_Mt = [];
density_all_Mt = [];
twoPts_corrMt = [];

matterF_Mt = [];
gaugeF_Mt = [];

psic = phi_init;
probl = abs(conj(psic).*psic);

stat_nC = site_stat_1dqlm_FCN(psic,basis_qlm,init_matter_occuP);
density_gauge_Mt = cat(1,density_gauge_Mt,stat_nC.density_guage_ltC);
density_matter_Mt = cat(1,density_matter_Mt,stat_nC.density_matter_ltC);
density_all_Mt = cat(1,density_all_Mt,stat_nC.density_all);

matterF_Mt = cat(1,matterF_Mt,stat_nC.matter_fieldC);
gaugeF_Mt = cat(1,gaugeF_Mt,stat_nC.gauge_fieldC);

twoPts_corrMt = cat(1,twoPts_corrMt,stat_nC.two_pts_correlation_ltC);

tStart = tic; 
for kk = 2:nt
    fprintf('Current process: %04d / %04d.\n',kk,nt)
    
    tSC = tic;
    psic = expv(-1i*dt,ham,psic,1.0e-8,50);
    % psic = expm(-1i*ham*dt)*psic;
    tEC = toc(tSC);
    fprintf('time for evolution: %.6f seconds.\n',tEC)
    
    stat_nC = site_stat_1dqlm_FCN(psic,basis_qlm,init_matter_occuP);
    density_gauge_Mt = cat(1,density_gauge_Mt,stat_nC.density_guage_ltC);
    density_matter_Mt = cat(1,density_matter_Mt,stat_nC.density_matter_ltC);
    density_all_Mt = cat(1,density_all_Mt,stat_nC.density_all);
    
    matterF_Mt = cat(1,matterF_Mt,stat_nC.matter_fieldC);
    gaugeF_Mt = cat(1,gaugeF_Mt,stat_nC.gauge_fieldC);
    twoPts_corrMt = cat(1,twoPts_corrMt,stat_nC.two_pts_correlation_ltC);
    
    % fprintf('\n')
end
tEnd = toc(tStart);
fprintf('\nElapsed time is %.6f seconds.\n',tEnd)


%%
x = 0:n_matter-1;
y = tl*1000;

figure('Color','w','Position',[120 120 560 420])
imagesc(x,y,twoPts_corrMt)
colormap(hot)
clb = colorbar;
clb.Title.String = '';
clb.Title.FontSize = 14;
ax = gca;
ax.FontSize = 14;
xlabel('\itr','FontSize',16)
ylabel('Evolution time (ms)','FontSize',16)

% export_fig('-r300','TwoPts_correlation_Conf.png')

%
x = 1:L;
y = tl*1000;

figure('Color','w','Position',[720 120 560 420])
imagesc(x,y,density_all_Mt)
colormap(hot)
clb = colorbar;
clb.Title.String = '\rho';
clb.Title.FontSize = 14;
ax = gca;
ax.FontSize = 14;
xlabel('Sites','FontSize',16)
ylabel('Evolution time (ms)','FontSize',16)

% export_fig('-r300','Density_profile_Conf.png')

% figure('Color','w','Position',[120 120 560 420])
% imagesc(x,y,gaugeF_Mt)
% ax = gca;
% ax.FontSize = 14;
% xlabel('Sites','FontSize',16)
% ylabel('Evolution time (ms)','FontSize',16)

% figure('Color','w','Position',[720 120 560 420])
% plot(tl*1000,-mean(gaugeF_Mt,2),'LineWidth',2)
% ax = gca;
% ax.FontSize = 14;
% xlabel('Evolution time (ms)','FontSize',16)
% ylabel('\langle \itE \rangle','FontSize',16)


return

