function [ data ] = load_netsim_individual( filepath )

assert( ~isempty( regexp( filepath, '\.bin?$', 'match', 'once' ) ) )

% open and read file
fid = fopen( filepath, 'rb' );
sz = fread( fid, 1, 'int' );
data = fread( fid, 'double' );
data = reshape(data, [], sz)';

% close file
fclose( fid );

end