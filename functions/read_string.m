function str = read_string(fid) 
    str = "";
    ch = fread(fid, 1, 'char=>char')';
    while (ch ~= 0)
        str = str + ch;
        ch = fread(fid, 1, 'char=>char')';
    end
end
