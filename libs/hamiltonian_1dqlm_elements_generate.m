function [ham_elems] = hamiltonian_1dqlm_elements_generate(basis_qlm)



%%
ham_elems = struct();


%%
ns = basis_qlm.ns;
matter_stateS = basis_qlm.matter_stateS;
gauge_stateS = basis_qlm.gauge_stateS;
state_idxlt = basis_qlm.state_idxlt;

% init_Gl = basis_qlm.init_Gl;
% gauge_edge = basis_qlm.gauge_edge;

n_matter = basis_qlm.n_matter;
n_gauge = basis_qlm.n_gauge;


%% diag-terms,
i_lt = (1:ns)';
j_lt = i_lt;

% 01: mass term
fprintf('\ngenerating hamiltonian elements of mass term.\n')
tVal_mass = tic;
matterS_dotM = matter_stateS.*matter_stateS;
ham_elems_mass = struct();
for kk = 1:n_matter
    field_cur = ['matterS_',num2str(kk)];
    ham_elems_mass.(field_cur) = sparse(i_lt,j_lt,...
        matterS_dotM(:,kk),ns,ns);
end
ham_elems_mass.sumAllSite = sparse(i_lt,j_lt,...
    sum(matterS_dotM,2),ns,ns);
tD_mass = toc(tVal_mass);
fprintf('elapsed time is %.6f seconds.\n',tD_mass)


% 02: background electronic term
fprintf('\ngenerating hamiltonian elements of electronic term.\n')
tVal_elec = tic;
ham_elems_elec = struct();
for kk = 1:n_gauge
    field_cur = ['gaugeS_',num2str(kk)];
    ham_elems_elec.(field_cur) = sparse(i_lt,j_lt,...
        gauge_stateS(:,kk),ns,ns);
end
ham_elems_elec.sumAllSite = sparse(i_lt,j_lt,...
    sum(gauge_stateS,2),ns,ns);
tD_elec = toc(tVal_elec);
fprintf('elapsed time is %.6f seconds.\n',tD_elec)


%% off-diag term, interactions
fprintf('\ngenerating hamiltonian elements of interaction term.\n')
tVal_inter = tic;

seq_stat_all = (1:ns)';

i_idx = [];
j_idx = [];
v_idx = [];

idx_i_all = struct();
idx_j_all = struct();
idx_k_all = struct();
for kk = 1:n_gauge
    field_cur = ['mgm_',num2str(kk),'_',num2str(kk+1)];
    idx_i_all.(field_cur) = [];
    idx_j_all.(field_cur) = [];
    idx_k_all.(field_cur) = [];
end


%% Method 01
vC = -1;
for jj = 1:n_gauge
    field_cur = ['mgm_',num2str(jj),'_',num2str(jj+1)];
    
    % 01: '++-' terms
    gauge_nxt = gauge_stateS;
    matter_nxt = matter_stateS;
    
    gauge_nxt(:,jj) = gauge_nxt(:,jj) + 1;
    matter_nxt(:,jj) = matter_nxt(:,jj) + 1;
    matter_nxt(:,jj+1) = matter_nxt(:,jj+1) - 1;
    
    matter_nxt_od = matter_nxt(:,1:2:end);
    matter_nxt_en = matter_nxt(:,2:2:end);
    idx_sel = (sum(abs(matter_nxt) <= 1,2) == n_matter) ...
        & (sum(matter_nxt_od >= 0,2) == size(matter_nxt_od,2)) ...
        & (sum(matter_nxt_en <= 0,2) == size(matter_nxt_en,2)) ...
        & (sum(abs(gauge_nxt) < 1,2) == n_gauge);
    
    seq_stat_cur = seq_stat_all(idx_sel);
    gauge_nxt = gauge_nxt(idx_sel,:);
    matter_nxt = matter_nxt(idx_sel,:);
    idx_mt_nxt = abs(matter_nxt) * (2.^(n_matter-1:-1:0))';
    idx_gt_nxt = (gauge_nxt+0.5) * (2.^(n_gauge-1:-1:0))';
    
    [~,seq_stat_nxt] = ismember([idx_mt_nxt,idx_gt_nxt],...
        state_idxlt(:,1:2),'rows');
    
    i_idx = cat(1,i_idx,seq_stat_nxt);
    j_idx = cat(1,j_idx,seq_stat_cur);
    v_idx = cat(1,v_idx,vC*ones(size(seq_stat_nxt)));
    
    idx_i_all.(field_cur) = cat(1,...
        idx_i_all.(field_cur),seq_stat_nxt);
    idx_j_all.(field_cur) = cat(1,...
        idx_j_all.(field_cur),seq_stat_cur);
    idx_k_all.(field_cur) = cat(1,...
        idx_k_all.(field_cur),vC*ones(size(seq_stat_nxt)));
    

    % 02: '--+' terms
    gauge_nxt = gauge_stateS;
    matter_nxt = matter_stateS;
    
    gauge_nxt(:,jj) = gauge_nxt(:,jj) - 1;
    matter_nxt(:,jj) = matter_nxt(:,jj) - 1;
    matter_nxt(:,jj+1) = matter_nxt(:,jj+1) + 1;
    
    matter_nxt_od = matter_nxt(:,1:2:end);
    matter_nxt_en = matter_nxt(:,2:2:end);
    idx_sel = (sum(abs(matter_nxt) <= 1,2) == n_matter) ...
        & (sum(matter_nxt_od >= 0,2) == size(matter_nxt_od,2)) ...
        & (sum(matter_nxt_en <= 0,2) == size(matter_nxt_en,2)) ...
        & (sum(abs(gauge_nxt) < 1,2) == n_gauge);
    
    seq_stat_cur = seq_stat_all(idx_sel);
    gauge_nxt = gauge_nxt(idx_sel,:);
    matter_nxt = matter_nxt(idx_sel,:);
    idx_mt_nxt = abs(matter_nxt) * (2.^(n_matter-1:-1:0))';
    idx_gt_nxt = (gauge_nxt+0.5) * (2.^(n_gauge-1:-1:0))';
    
    [~,seq_stat_nxt] = ismember([idx_mt_nxt,idx_gt_nxt],...
        state_idxlt(:,1:2),'rows');
    
    i_idx = cat(1,i_idx,seq_stat_nxt);
    j_idx = cat(1,j_idx,seq_stat_cur);
    v_idx = cat(1,v_idx,vC*ones(size(seq_stat_nxt)));
    
    idx_i_all.(field_cur) = cat(1,...
        idx_i_all.(field_cur),seq_stat_nxt);
    idx_j_all.(field_cur) = cat(1,...
        idx_j_all.(field_cur),seq_stat_cur);
    idx_k_all.(field_cur) = cat(1,...
        idx_k_all.(field_cur),vC*ones(size(seq_stat_nxt)));
    
end


%%
for kk = 1:n_gauge
    field_cur = ['mgm_',num2str(kk),'_',num2str(kk+1)];
    i_cur = idx_i_all.(field_cur);
    j_cur = idx_j_all.(field_cur);
    k_cur = idx_k_all.(field_cur);
    
	ham_elems_interaction.(field_cur) = sparse(i_cur,j_cur,k_cur,ns,ns);

end
ham_elems_interaction.sumAllSite = sparse(i_idx,j_idx,v_idx,ns,ns);

tD_inter = toc(tVal_inter);
fprintf('elapsed time is %.6f seconds.\n',tD_inter)


%% results
ham_elems.ham_elems_mass = ham_elems_mass;
ham_elems.ham_elems_interaction = ham_elems_interaction;
ham_elems.ham_elems_elec = ham_elems_elec;

ham_elems.basis_qlm = basis_qlm;


