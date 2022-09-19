function [ham_elems] = hamiltonian_1dqlm_elements_generate(basis_qlm)


%%
ham_elems = struct();


%% Input parameter check
TAG_s = true;
try
    cache_dir = [pwd,filesep,'cache'];
    cache_info = load([cache_dir,filesep,'ham_elems_cache.mat']);
    if (cache_info.basis_qlm.L ~= basis_qlm.L) ...
            || sum(cache_info.basis_qlm.init_Gl ~= basis_qlm.init_Gl) ...
            || sum(cache_info.basis_qlm.gauge_edge ~= basis_qlm.gauge_edge)
        TAG_s = false; 
    end
catch
    TAG_s = false;    
end

%
if TAG_s
    ham_elems = cache_info.ham_elems;
    fprintf('Using the cache data!\n')
    return
else
    fprintf('Re-generate the hamiltonian elements!\n')
end


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
% 01: mass term
matterS_dotM = matter_stateS.*matter_stateS;
ham_elems_mass = struct();
for kk = 1:n_matter
    field_cur = ['matterS_',num2str(kk)];
    ham_elems_mass.(field_cur) = sparse(diag(matterS_dotM(:,kk)));
end
ham_elems_mass.sumAllSite = sparse(diag(sum(matterS_dotM,2)));


% 02: background electronic term
ham_elems_elec = struct();
for kk = 1:n_gauge
    field_cur = ['gaugeS_',num2str(kk)];
    ham_elems_elec.(field_cur) = sparse(diag(gauge_stateS(:,kk)));
end
ham_elems_elec.sumAllSite = sparse(diag(sum(gauge_stateS,2)));


%% off-diag term, interactions
% ham_elems_interaction = sparse(ns,ns);
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

% matter-gauge field couping 
for kk = 1:ns
    state_idx_cur = kk;
    
    gauge_cur = gauge_stateS(kk,:);
    matter_cur = matter_stateS(kk,:);
    
    for jj = 1:n_gauge
        field_cur = ['mgm_',num2str(jj),'_',num2str(jj+1)];
        
        gauge_nxt = gauge_cur;
        matter_nxt = matter_cur;
        
        % '++-' terms
        if gauge_cur(jj) < 0
            gauge_nxt(jj) = gauge_nxt(jj) + 1;
            matter_nxt(jj) = matter_nxt(jj) + 1;
            matter_nxt(jj+1) = matter_nxt(jj+1) - 1;
            
            vC = -1;

        end

        % '--+' terms
        if gauge_cur(jj) > 0
            gauge_nxt(jj) = gauge_nxt(jj) - 1;
            matter_nxt(jj) = matter_nxt(jj) - 1;
            matter_nxt(jj+1) = matter_nxt(jj+1) + 1;
            
            vC = -1;

        end
        
        % bypass
        if sum(abs(matter_nxt) > 1) || sum(matter_nxt(1:2:end) < 0) ...
                || sum(matter_nxt(2:2:end) > 0)
            continue
        end

        % idx_nxt
        idxlt_nxt = [base2dec(strrep(num2str(abs(matter_nxt)),' ',''),2),...
            base2dec(strrep(num2str(gauge_nxt+0.5),' ',''),2)];

        try
            [~,state_idx_nxt] = ismember(idxlt_nxt,...
                state_idxlt(:,1:2),'rows');
            
            i_idx = cat(1,i_idx,state_idx_nxt);
            j_idx = cat(1,j_idx,state_idx_cur);
            v_idx = cat(1,v_idx,vC);
            
            idx_i_all.(field_cur) = cat(1,...
                idx_i_all.(field_cur),state_idx_nxt);
            idx_j_all.(field_cur) = cat(1,...
                idx_j_all.(field_cur),state_idx_cur);
            idx_k_all.(field_cur) = cat(1,...
                idx_k_all.(field_cur),vC);
            
        catch
            continue
        end

    end
    
end

%
for kk = 1:n_gauge
    field_cur = ['mgm_',num2str(kk),'_',num2str(kk+1)];
    i_cur = idx_i_all.(field_cur);
    j_cur = idx_j_all.(field_cur);
    k_cur = idx_k_all.(field_cur);
    
    ham_elems_interaction.(field_cur) = sparse(i_cur,j_cur,k_cur,ns,ns);

end
ham_elems_interaction.sumAllSite = sparse(i_idx,j_idx,v_idx,ns,ns);


%% results
ham_elems.ham_elems_mass = ham_elems_mass;
ham_elems.ham_elems_interaction = ham_elems_interaction;
ham_elems.ham_elems_elec = ham_elems_elec;

ham_elems.basis_qlm = basis_qlm;


%% write cache info. into disk
try
    cache_dir = [pwd,filesep,'cache'];
    if ~exist(cache_dir,'dir')
        mkdir(cache_dir)
    end
    
    save([cache_dir,filesep,'ham_elems_cache.mat'],...
        'basis_qlm','ham_elems');
    
catch
    
end
