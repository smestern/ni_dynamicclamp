function params = importParams( paramfile )

% open and read file
fid = fopen( paramfile, 'r' );
if fid == -1, error( 'Bad filepath.' ); end

% parse and save each line of parameter file
line = fgets( fid );
while line ~= -1
    
    if ( line(1) ~= '#' ) && ~( any( (line == '(') | (line == '[') | (line == '{') ) )
        if line(1) == '$'
            
            s = textscan( line, '$%s = %s' );
            if ~isempty( s ) && ( ~isempty(s{1}) ) %#ok<ALIGN>
                eval( strcat( 'params.', s{1}{1}, ' = ''', char(s{2}), ''';') );
            end
            
        elseif line(1) == '@'
            s = textscan( line, '%s = %f : %f : %f' );
            if ~isempty( s ) && ( ~isempty(s{1}) ) %#ok<ALIGN>
                eval( strcat( 'params.', s{1}{1}(2:end), ' = ','[',num2str( linspace(s{2},s{3},s{4})),']',';') );
            end
            
        else
            
            s = textscan( line, '%s = %f' );
            if ~isempty( s ) && ( ~isempty(s{1}) ) %#ok<ALIGN>
                eval( strcat( 'params.', s{1}{1}, ' = ', num2str(s{2}), ';') );
            end
            
        end
        
    end
    
    line = fgets( fid );
    
end

% close file
fclose( fid );

if isfield(params,'output_file_code')
    params = rmfield(params,'output_file_code');
end

if params.connector == 4, params.dim = 2; else, params.dim = 1; end
params.Ne = floor(0.8*params.N);
params.Ni = params.N - params.Ne;

params.bin_size_x = params.bin_size ^ (1/params.dim) ; % number of neurons in one dimension of one bin (if 1D network, same as bin_size)
params.nbins = params.N  /  ( params.bin_size ^ params.dim ) ; % number of bins (exc+inh)
params.nbins_exc = 0.8 * params.nbins; % number of bins (exc)
params.nbins_inh = params.nbins - params.nbins_exc; % number of bins (inh)

DT = params.dt*params.record_downsample_factor;
params.t = DT : DT : params.T;

end
