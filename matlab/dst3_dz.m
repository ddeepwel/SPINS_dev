function dd = dst3_dz(f)
% Takes the derivative along the third (z) dimension of
% three-dimensional array f via the DST (sine transform).
% This does not normalize for length, so the result should
% be multiplied by pi/Lz.

[Nx Ny Nz] = size(f);
persistent ls
% Create an array for the wavenumbers ls.  This is
% declared 'persistent' so it is -not- regenerated
% at each call.  In MATLAB, 'persistent' variables
% stay around after execution like 'global' ones, but
% do not show up in the global namespace.

% Since DCT and DST transforms have implied symmetry
% -and- MATLAB doesn't provide good built-in DCT/DST
% operators (only in the JPEG toolbox, which is itself
% highly specialized), this code fakes the real-transform
% with a full FFT on a copied matrix.  The upshot is
% that the wavenumber vector is doubled.
if (~exist('ls') || size(ls,1) ~= Nx || size(ls,3) ~= 2*Nz || size(ls,2) ~= Ny)
    ls = [0:(Nz-1) -Nz:-1];
    ls = permute(ls,[1 3 2]);
    ls = repmat(ls,[Nx Ny 1]);
end
% For the Fourier transform, f gets doubled.  A DCT doubles
% with even symmetry, and a DST doubles with odd symmetry.

dd = real(ifft(1i.*ls.*fft(cat(3,f,-flipdim(f,3)),[],3),[],3));

% Now, dd contains the derivative of f and its symmetric 
% extension, so give it a haircut.

dd = dd(1:Nx,1:Ny,1:Nz);
