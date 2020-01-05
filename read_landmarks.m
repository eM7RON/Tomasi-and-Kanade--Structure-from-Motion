function shape = read_landmarks(filename_lndm)
% Read landmarks from the file with name filename_lndm
%
% (OUTPUT):
% shape: Nlandmarks x 2 array whose columns contain the X,Y coordinated of the
% 2D landmarks (in pixels)


fid = fopen(filename_lndm);
tline = fgetl(fid);
start = 1;
while ~isstrprop(tline(1), 'digit')
    start = start + 1;
    tline = fgetl(fid);
end
end_line = start-1;
while isstrprop(tline(1), 'digit')
    end_line = end_line + 1;
    tline = fgetl(fid);
end

fclose(fid);

% read shape with dlmream:
shape =  dlmread(filename_lndm, ' ', [start-1 0 end_line-1 1]);

