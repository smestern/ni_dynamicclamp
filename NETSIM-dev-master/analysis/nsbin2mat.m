function [params, si, st, v, lz, ge, gi, ge_indiv, gi_indiv, v_indiv] = nsbin2mat(paramfile, outputdir, jobid)
% Gabriel Benigno, 2021

arguments
    paramfile
    outputdir
    jobid = 0
end

disp('CONVERTING BINARY DATA TO MATLAB FORMAT')

% import parameters
params = importParams( paramfile );

% import spike, v, ge, gi, individual ge data
[si, st] = load_netsim_spikes(fullfile(outputdir, sprintf('%08dspk.bin', jobid ) ) , params );
v = load_netsim_data(fullfile(outputdir, sprintf('%08dvms.bin', jobid ) ) , params );
ge = load_netsim_data(fullfile(outputdir, sprintf('%08dge.bin', jobid ) ) , params );
gi = load_netsim_data(fullfile(outputdir, sprintf('%08dgi.bin', jobid ) ) , params );
ge_indiv = load_netsim_individual(fullfile(outputdir, sprintf('%08dindividualge.bin',jobid)));
ge_indiv = load_netsim_individual(fullfile(outputdir, sprintf('%08dindividualgi.bin',jobid)));
v_indiv = load_netsim_individual(fullfile(outputdir, sprintf('%08dindividualvms.bin',jobid)));

% calculate lfp based on v, ge, gi
lz = nslfp(params,v,ge,gi);

% remove initial transient of spikes, v, ge, gi, lfp
it = 0.2;
[params,v,ge,gi,lz,st,si,ge_indiv] = nsrminit(it,params,v,ge,gi,lz,st,si,ge_indiv);

% save data into mat files
save( sprintf( '%s/params.mat' , outputdir ) , 'params' )
save( sprintf( '%s/ge.mat' , outputdir ) , 'ge' )
save( sprintf( '%s/gi.mat' , outputdir ) , 'gi' )
save( sprintf( '%s/v.mat' , outputdir ) , 'v' )
save( sprintf( '%s/lz.mat' , outputdir ) , 'lz' )
save( sprintf( '%s/ge_indiv.mat' , outputdir ) , 'ge_indiv' )
save( sprintf( '%s/gi_indiv.mat' , outputdir ) , 'gi_indiv' )
save( sprintf( '%s/v_indiv.mat' , outputdir ) , 'v_indiv' )

% save note that first 0.2 s have already been removed from the data
fid = fopen(sprintf('%s/NOTE_INITIAL_TRANSIENT_ALREADY_REMOVED',outputdir),'w');
fclose(fid);

disp('DONE')

%{
% notes
notes = ["p struct contains parameters from parameter file";
             "all times in seconds; all positions and distances in meters";
             "st and si form ordered pairs of spike time and neuron index, respectively";
              "lz is z-scored LFP"; 
              "time series data (v, ge, gi, lz) are downsampled in time by a factor of vars.record_downsample_factor";
              "time series data are downsampled in space by a factor of (vars.bin_size)^dim, where dim is the spatial dimension";
              "initial transient of simulation (the first 0.2 seconds of data) has been removed from time series data and spike data (st and si)";
              "t is time vector corresponding to time series data, excluding initial transient and set as beginning at 0 seconds";
              "st has been adjusted such time=0 starts just after the initial transient";
              "data are of exc. neurons only--no inh. were saved";
              "to reshape time series data in 2D space (if applicable), use reshape(xx, sqrt(nbins), sqrt(nbins), []), where nbins = Ne / vars.bin_size^2, where Ne = floor(0.8*vars.N)";
              "neuron indexing in 2D space is linear in the matlab sense. see sub2ind function documentation"];

% save
save(fullfile(outputdir,'p.mat'),'p')
save(fullfile(outputdir,'spikes.mat'),'st','si')
save(fullfile(outputdir,'v.mat'),'t','v')
save(fullfile(outputdir,'ge.mat'),'t','ge')
save(fullfile(outputdir,'gi.mat'),'t','gi')
save(fullfile(outputdir,'lz.mat'),'t','lz')
save(fullfile(outputdir,'notes.mat'),'notes')

disp('SELECTED DATA OUTPUT TO WORKSPACE AND ALL DATA SAVED IN OUTPUT DIRECTORY')
disp(' ')
%}

end