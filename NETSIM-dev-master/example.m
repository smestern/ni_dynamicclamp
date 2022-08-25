% Oct 11 2021
% Examples for running NETSIM and loading outputs into MATLAB
% Gabriel Benigno

%% Important: the 'analysis' and 'src' folders must be added to the path

%% EXAMPLE 1: basic 1D network
addpath('./analysis')
addpath('./src')
clearvars; clc
paramfile = './tests/test1.parameters';
outputdir = '~/tmpdir';
jobid = 2;

runNetsim( paramfile , outputdir , jobid );
[params, si, st, v, lz, ge, gi] = nsbin2mat( paramfile , outputdir , jobid );
rmpath('./analysis')
rmpath('./src')

%% EXAMPLE 2: basic 2D network
clearvars; clc
paramfile = './tests/test7.parameters';
outputdir = '~/tmpdir';
jobid = 1;

netsim( paramfile , outputdir , jobid );
[params, si, st, v, lz, ge, gi] = nsbin2mat(paramfile, outputdir, jobid);
rmpath('./analysis')

%% EXAMPLE 3: sparse network
clearvars; clc
addpath('./src')
addpath('./analysis')

outputdir = '../sparse';
paramfile = sprintf('%s/sparse.parameters',outputdir);

makeParameterFile(paramfile, ...
                    'N' , 450000 , ...
                    'K' , 3000 , ...
                    'T' , 1.2 , ...
                    'ge' , 1e-9 , ...
                    'gi' , 10e-9 , ...
                    'connector' , 4 , ...
                    'L' , 4e-3 , ...
                    'p_release' , 0.1 , ...
                    'sigma_space' , .4e-3 );
runNetsim( paramfile , outputdir , 'release_probability' , true );
[params, si, st, v, lz, ge, gi, ge_indiv, gi_indiv, v_indiv] = nsbin2mat(paramfile, outputdir, 0);

rmpath('./src')
rmpath('./analysis')

%% EXAMPLE 4: dense network
clearvars; clc
addpath('./src')
addpath('./analysis')

outputdir = '../dense';
paramfile = sprintf('%s/dense.parameters',outputdir);

makeParameterFile(paramfile, ...
                    'N' , 12500 , ...
                    'K' , 100 , ...
                    'T' , 1.2 , ...
                    'ge' , 9.3157e-09 , ...
                    'gi' , 1.1754e-06 , ...
                    'connector' , 4 , ...
                    'L' , 4e-3 , ...
                    'p_release' , 0.1 , ...
                    'sigma_space' , .4e-3 );
runNetsim( paramfile , outputdir , 'release_probability' , true );
[params, si, st, v, lz, ge, gi, ge_indiv, gi_indiv, v_indiv] = nsbin2mat(paramfile, outputdir, 0);

rmpath('./src')
rmpath('./analysis')