%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   inv_f_trans_on_vector
%
%   Description:
%   Makes use of the FFT algorithm in order to conveniently calculate
%   the inverse Fourier transform of the function which is sampled in the
%   arguments.
%
%   Notes:
%   - FFT/IFFT are fastest when N = 2^p (for N even - when N is odd, this
%     function will be fastest for N = 2^p + 1).
%   - See associated notes for an explanation of the code.
%   - FFT is used for this algorithm due to sign convention.
%
%   Arguments:
%   xilist - A list of equally spaced temporal frequencies with length N. 
% If N is
%           odd, then the values should be centered about zero. If N is
%           even, then the values should be in the form -W:dw:W-dw.
%   fhatlist - A list of the values of the function being inverse Fourier 
%           transformed evaluated at the points in xilist.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




function [tlist,flist]=inv_f_trans_on_vector_3(xilist,fhatlist)


s1 = size(xilist);
s2 = size(fhatlist);

if s1(2) > 1
    if s1(1) > 1
        error('xilist must be an Nx1 vector.')
    else
        xilist = xilist';
    end
end

N=length(xilist);

if s2(2) > 1
    if (s2(1) > 1) || (s2(2) ~= N)
        error('fhatlist must be an Nx1 vector,and size must match xilist.')
    else
        fhatlist = fhatlist';
    end
end

jlist=(0:N-1)';

if(mod(N,2)==0)

    duration=1/(xilist(2)-xilist(1));
    dt=duration/N;

    tlist=-duration/2 + dt*jlist;

    gtildelist=fhatlist.*exp(-pi*1i*jlist);
    ftildelist=fft(gtildelist);

    flist=1/(dt*N) * exp(pi*1i*N/2) * exp(-pi*1i*jlist) .*ftildelist;
else

    duration=1/(xilist(2)-xilist(1));
    dt=duration/(N-1);
    
    tlist = -duration/2 + dt*jlist;
    
    gtildelist = fhatlist.*exp(-pi*1i*jlist);
    ftildelist = fft(gtildelist(1:end-1)) + gtildelist(end);
    ftildelist(end+1) = sum(gtildelist);
    
    flist =...
        1./(dt*(N-1)).* exp(pi*1i*(N-1)/2).*exp(-pi*1i*jlist).*ftildelist;
end
    
    

end

    