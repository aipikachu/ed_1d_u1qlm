function [basis] = boson_basis_1dqlm(L,init_Gl,gauge_edge)
% [basis] = boson_basis_1dqlm(L,init_Gl,gauge_edge)
%
%


%%
basis = struct();


%%
n_matter = length(init_Gl);
n_gauge = L - n_matter;


%%
tic
gauge_state_raw = [];
for kk = 0:n_gauge
    basis_cur = boson_basis_1d(n_gauge,kk,1);
    gauge_state_raw = cat(1,gauge_state_raw,...
        basis_cur.state-0.5);
end
toc


%%
gauge_raw_num = size(gauge_state_raw,1);
gauge_stateS_RAW = [repmat(gauge_edge(1),gauge_raw_num,1),...
    gauge_state_raw,repmat(gauge_edge(2),gauge_raw_num,1)];


matter_stateS_RAW = gauge_stateS_RAW(:,2:end) ...
    - gauge_stateS_RAW(:,1:end-1) ...
    - repmat(init_Gl,gauge_raw_num,1);
matter_stateS_RAW_Od = matter_stateS_RAW(:,1:2:end);
matter_stateS_RAW_En = matter_stateS_RAW(:,2:2:end);


% for matter sites: 0/1 at 'odd'-sites, and -1/0 at 'even'-sites
idx_Gl_Od = sum(matter_stateS_RAW_Od >= 0,2) ...
    == size(matter_stateS_RAW_Od,2);
idx_Gl_En = sum(matter_stateS_RAW_En <= 0,2) ...
    == size(matter_stateS_RAW_En,2);

idx_Gl = idx_Gl_Od & idx_Gl_En;

gauge_stateS = gauge_state_raw(idx_Gl,:);
gauge_stateS_Ep = gauge_stateS_RAW(idx_Gl,:);
matter_stateS = matter_stateS_RAW(idx_Gl,:);

ns = sum(idx_Gl);


%%
% matter_site_idx, gauge_site_idx,gauge_expand_site_idx
state_idxlt = zeros(ns,3);  
for kk = 1:ns
    matterS_cur = abs(matter_stateS(kk,:));
    gaugeS_cur = gauge_stateS(kk,:)+0.5;
    gaugeSE_cur = gauge_stateS_Ep(kk,:)+0.5;
    
    state_idxlt(kk,:) = [base2dec(strrep(num2str(matterS_cur),' ',''),2),...
        base2dec(strrep(num2str(gaugeS_cur),' ',''),2),...
        base2dec(strrep(num2str(gaugeSE_cur),' ',''),2)];
end


%% in atom number case
matter_atom_stateS = abs(matter_stateS);
gauge_atom_stateS = zeros(size(gauge_stateS));
atom_occupationS = zeros(ns,L);

gauge_atom_stateS(:,1:2:end) = -2*(gauge_stateS(:,1:2:end)-0.5);
gauge_atom_stateS(:,2:2:end) = 2*(gauge_stateS(:,2:2:end)+0.5);

atom_occupationS(:,1:2:end) = matter_atom_stateS;
atom_occupationS(:,2:2:end) = gauge_atom_stateS;


%%
basis.matter_stateS = matter_stateS;
basis.gauge_stateS = gauge_stateS;
basis.gauge_stateS_Ep = gauge_stateS_Ep;
basis.ns = ns;
basis.matter_atom_stateS = matter_atom_stateS;
basis.gauge_atom_stateS = gauge_atom_stateS;
basis.atom_occupationS = atom_occupationS;
basis.state_idxlt = state_idxlt;

basis.n_matter = n_matter;
basis.n_gauge = n_gauge;

basis.L = L;
basis.init_Gl = init_Gl;
basis.gauge_edge = gauge_edge;

