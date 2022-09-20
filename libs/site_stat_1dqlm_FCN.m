function stat_n = site_stat_1dqlm_FCN(psic,basis_qlm,init_matter_occuP)


%%
if nargin<3
    % for two-point correlation statistics
    init_matter_occuP = zeros(1,basis_qlm.n_matter);
end
    

%%
stat_n = struct();


%%
ns = basis_qlm.ns;
L = basis_qlm.L;

n_matter = basis_qlm.n_matter;
n_gauge = basis_qlm.n_gauge;


matter_atom_stateS = basis_qlm.matter_atom_stateS;
gauge_atom_stateS = basis_qlm.gauge_atom_stateS;

matter_stateS = basis_qlm.matter_stateS;
gauge_stateS = basis_qlm.gauge_stateS;


%% probability list
prob_ltC = abs(conj(psic).*psic);


%%
matter_fieldC = sum(matter_stateS .* repmat(prob_ltC,1,n_matter),1);
gauge_fieldC = sum(gauge_stateS .* repmat(prob_ltC,1,n_gauge),1);


%% density profile
density_matter_ltC = sum(matter_atom_stateS .* repmat(prob_ltC,1,n_matter),1);
density_guage_ltC = sum(gauge_atom_stateS .* repmat(prob_ltC,1,n_gauge),1);

density_all = zeros(1,L);
density_all(1:2:end) = density_matter_ltC;
density_all(2:2:end) = density_guage_ltC;


%% two-point correlations
% init_matter_occuP
% matter_atom_stateS
two_pts_correlation_ltC = zeros(1,n_matter);

for kk = 1:n_matter
    r = kk - 1;
    G_cur = zeros(ns,1);
    try
        for jj = 1:n_matter
            G_cur = G_cur + (matter_atom_stateS(:,jj) - init_matter_occuP(jj)) ...
                .* (matter_atom_stateS(:,jj+r) - init_matter_occuP(jj+r));
        end        
    catch
        
    end    
    two_pts_correlation_ltC(kk) = sum(G_cur .* prob_ltC,1,'omitnan');   

end



%%
stat_n.ns = ns;
stat_n.L = L;
stat_n.prob_ltC = prob_ltC;
stat_n.matter_fieldC = matter_fieldC;
stat_n.gauge_fieldC = gauge_fieldC;


% density profile
stat_n.density_matter_ltC = density_matter_ltC;
stat_n.density_guage_ltC = density_guage_ltC;
stat_n.density_all = density_all;
stat_n.two_pts_correlation_ltC = two_pts_correlation_ltC;


