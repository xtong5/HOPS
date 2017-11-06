%% tool functions
% compute_errors_2d_polar.m 
%                 ---compute error between real Nuemann data and estimeted
%                 data by DNO solvers
% diff_bessel.m   ---nth derivative of bessel function of first kind 
%                 ---OR hankel function of first kind 
% diff_besselh.m  ---nth derivative of hankel function of first kind
% diff_besselj.m  ---nth derivative of bessel function of first kind
% dno_fe_helmholtz_polar_exterior.m 
%                 ---Dirichlet to Nuemann operater using exterior data
% dno_fe_helmholtz_polar_interior.m
%                 ---Dirichlet to Nuemann operater using exterior data
% field_fe_helmholtz_polar_exterior.m 
%                 ---Single layer Field expansion (exterior)
% field_fe_helmholtz_polar_interior.m 
%                 ---Single layer Field expansion (interior)
% field_fe_helmholtz_twolayer_polar.m 
%                 ---Two layer Field expansion (NOT in use)
% make_plots_polar.m 
%                 ---plot figures using errors
% padesum.m       ---Sums a truncated Taylor series via Pade approximation
% taylorsum.m     ---Sums a truncated Taylor series
% ri_perm.m       ---generate refractive index
% drude.m         ---data used in ri_perm.m
% twolayer_dno_fe_helmholtz_polar.m
%                 ---Two layer Field expansion using DNO


%% test files
% plot_f.m        ---plot the functions 
% test_helmholtz_polar_exterior.m
%                 ---test Helmholtz DNO solvers in polar (exterior)
% test_helmholtz_polar_interior.m
%                 ---test Helmholtz DNO solvers in polar (interior)
% test_helmholtz_polar_twolayer.m
%                 ---test Helmholtz DNO solvers in polar (two layers)

%% functions generating data
% FE_app_pade.m   ---data given by different lambda and epsilon (padesum)
% FE_app_taylor.m ---data given by different lambda and epsilon (taylorsum)

%% data files
% make_data_all_eps.m
%                 ---generate data using index of real materials 

%% data folder
% data_material(pade)
%                 ---includes the data and plot files
% data_paper
%                 ---includes the data and plot files for paper






