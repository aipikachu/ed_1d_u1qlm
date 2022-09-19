function [ham_mat] = hamiltonian_1dqlm(basis_qlm,ham_elems,w,m,h)

%
if nargin < 5
    h = 0;
end


%%
ham_mat = struct();


%%
ns = basis_qlm.ns;
L = basis_qlm.L;

n_matter = basis_qlm.n_matter;
n_gauge = basis_qlm.n_gauge;


%% hamiltonian generate
ham = sparse(ns,ns);

% 01: matter-gauge coupling
len_w = length(w);
if len_w == 1
    ham = ham + w * ham_elems.ham_elems_interaction.sumAllSite;
elseif len_w == n_gauge
    for kk = 1:n_gauge
        field_cur = ['mgm_',num2str(kk),'_',num2str(kk+1)];
        ham = ham + w(kk) * ham_elems.ham_elems_interaction.(field_cur);
    end
else
    error('Error! Invalid input parameter "w"!')
end

% 02: mass term
len_m = length(m);
if len_m == 1
    ham = ham + m * ham_elems.ham_elems_mass.sumAllSite;
elseif len_m == n_matter
    for kk = 1:n_matter
        field_cur = ['matterS_',num2str(kk)];
        ham = ham + m(kk) * ham_elems.ham_elems_mass.(field_cur);
    end
else
    error('Error! Invalid input parameter "m"!')
end

% 03: background electronic term
len_h = length(h);
if len_h == 1
    ham = ham + h * ham_elems.ham_elems_elec.sumAllSite;
elseif len_h == n_gauge
    for kk = 1:n_gauge
        field_cur = ['gaugeS_',num2str(kk)];
        ham = ham + h(kk) * ham_elems.ham_elems_elec.(field_cur);
    end
else
    error('Error! Invalid input parameter "h"!')
end


%%
ham_mat.ham_elems = ham_elems;
ham_mat.ham = ham;
ham_mat.ns = ns;
ham_mat.L = L;

