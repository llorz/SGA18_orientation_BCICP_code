function T = read_map_with_NaN(filename)
tmp = strsplit(filename,'.');
if length(tmp) == 1
    filename = [filename,'.map'];
end
if exist(filename,'file')
    fid = fopen(filename);
    T = cell2mat(textscan(fid,'%d\n'));
    fclose(fid);
    T = double(T);
    T(T == 0) = nan;
else
    error('File not found.')
end
end
