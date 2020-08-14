function para = Para(varargin)


p = inputParser;
defultPath = 'ParaSet.txt';
p.addOptional('parafile',defultPath);

p.parse(varargin{:});

parafile = p.Results.parafile;

fid = fopen(parafile);
if fid < 0 
	error(['Please Check if',parafile,'exist']);
end

C = textscan(fid,'%s%s');
para = containers.Map(C{1},C{2});

fclose(fid);
end
