function concatBinFiles(srcA, srcB, dst, chunkBytes)
% concatBinFiles  Stream-copy srcA then srcB into dst (binary append).
%
% PURPOSE:
%   Concatenate two binary files (typically Kilosort-style int16 .bin)
%   into a single output without loading either fully into memory.
%   Chunks default to 64 MiB which keeps memory bounded even for
%   tens-of-GB recordings.
%
% USAGE:
%   concatBinFiles('A.bin', 'B.bin', 'merged.bin');
%   concatBinFiles('A.bin', 'B.bin', 'merged.bin', 16*1024*1024);
%
% Last modified 09.06.2026 (Jesus)

    if nargin < 4 || isempty(chunkBytes), chunkBytes = 64 * 1024 * 1024; end
    assert(isfile(srcA), 'NGL:concatBin', 'Source A not found: %s', srcA);
    assert(isfile(srcB), 'NGL:concatBin', 'Source B not found: %s', srcB);

    fOut = fopen(dst, 'w');
    if fOut < 0, error('NGL:concatBin', 'Cannot open %s for writing.', dst); end
    cleanupOut = onCleanup(@() fclose(fOut));

    for src = {srcA, srcB}
        fIn = fopen(src{1}, 'r');
        if fIn < 0, error('NGL:concatBin', 'Cannot open %s for reading.', src{1}); end
        cleanupIn = onCleanup(@() fclose(fIn));
        while ~feof(fIn)
            buf = fread(fIn, chunkBytes, '*uint8');
            if isempty(buf), break, end
            fwrite(fOut, buf, 'uint8');
        end
        clear cleanupIn
    end
end
