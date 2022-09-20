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
% Sample number
n_sample = 50;
r_amp = 10.0;

% parameters
w0 = 25 * 2 * pi;
w = ones(1,n_gauge) * w0;    % gauge-matter coupling
m = -4 * abs(w0);     % matter field
% h = 0.0 * 2 * pi;     % background-electronic field 
h_mat = r_amp * randi([-1,1],n_sample,n_gauge) * 2 * pi; % background-electronic field 
% h = (rand(1,n_gauge)-0.5) * 20.0 * 2 * pi; % background-electronic field 


%% ramp time list
nt = 201;
T = 120 * 1e-03;
tl = linspace(0,T,nt);
dt = tl(2) - tl(1);


%% generate hamiltonian
tic
ham_elems = hamiltonian_1dqlm_elements_generate(basis_qlm);
toc


%%
sample_data_mt = [];
rand_sel_mt = rand(n_sample,nt);

tStart = tic; 
for tt = 1:n_sample
    fprintf('Current sample loop: %04d / %04d.\n',tt,n_sample)
    
    tSC = tic;
    
    h = h_mat(tt,:);
    ham_all = hamiltonian_1dqlm(basis_qlm,ham_elems,w,m,h);
    ham = ham_all.ham;
    fprintf('Number of nonzero hamiltonian elements: %d.\n',nnz(ham))
    
    % 
    sample_mt_curr = [];
    
    %
    psic = phi_init;
    probl = abs(conj(psic).*psic);
    
    B = cumsum(probl);
    r_valC = rand_sel_mt(tt,1);
    idx_cur = find(r_valC < B);
    idx_cur = idx_cur(1);
    sample_mt_curr = cat(1,sample_mt_curr,...
        basis_qlm.atom_occupationS(idx_cur,:));
    
    for kk = 2:nt
        psic = expv(-1i*dt,ham,psic,1.0e-8,50);
        
        probl = abs(conj(psic).*psic);
    
        B = cumsum(probl);
        r_valC = rand_sel_mt(tt,kk);
        idx_cur = find(r_valC < B);
        idx_cur = idx_cur(1);
        sample_mt_curr = cat(1,sample_mt_curr,...
            basis_qlm.atom_occupationS(idx_cur,:));
        
    end
    
    tEC = toc(tSC);    
    fprintf('time for evolution: %.6f seconds.\n',tEC)
    
    sample_data_mt = cat(3,sample_data_mt,...
        sample_mt_curr);
    
end

tEnd = toc(tStart);
fprintf('\nElapsed time is %.6f seconds.\n',tEnd)



%%
density_sample = mean(sample_data_mt,3);

gauge_Mt = zeros(size(sample_data_mt));
gauge_Mt(:,[2,6,14],:) = 2;
A = sample_data_mt-gauge_Mt;
A(:,9,:) = 0;
density_sample2 = mean(A,3);

%
x = 1:L;
y = tl*1000;

figure('Color','w','Position',[120 120 560 420])
imagesc(x,y,density_sample2)
colormap(hot)
clb = colorbar;
clb.Title.String = '\rho';
clb.Title.FontSize = 14;
ax = gca;
ax.FontSize = 14;
xlabel('Sites','FontSize',16)
ylabel('Evolution time (ms)','FontSize',16)


%
figure('Color','w','Position',[120 120 560 420])
imagesc(x,y,density_sample)
colormap(hot)
clb = colorbar;
clb.Title.String = '\rho';
clb.Title.FontSize = 14;
ax = gca;
ax.FontSize = 14;
xlabel('Sites','FontSize',16)
ylabel('Evolution time (ms)','FontSize',16)

% export_fig('-r300',sprintf('Sample%03d_densityProfile_Amp_%03dHz.png',n_sample,r_amp))




