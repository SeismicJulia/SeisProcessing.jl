"""
    SeisRadon_fx(In,Par,itype);
RADON_FX: Operators for f-x, tau-p forward and adjoint Radon transforms
          to compute linear or  parabolic Radon transform.
IN    Radon coefficients if itype =  1 (Forward transform)  with In(np,nt)
      CMP gather         if itype = -1 (Adjoint transform)  with In(nh,nt)

OUT   CMP gather         if itype =  1 (Forward transform)  with Out(nh,nt)
      Radon coeffcients  if itype = -1 (Adjoint transform)  with Out(np,nt)
# Parameters
- `Par.h`:vector containing the nh offsets
- `Par.p`:vector containing the np curvatures normalized by far offset (s) if itype='parab'
- `Par.p`:vector containing the np dips in s/m if itype='linear'
- `Par.dt`:sampling interval
- `Par.f`:frequency corners of BP operator - this acts like a zero phase wavelet.
- `Par.transf `:'parab', 'linear'
itype =  1  :  forward
itype = -1  :  adjoint
"""
function SeisRadon_fx(In,Par,itype);
    h = Par[1][:h];
    p = Par[2][:p];
    dt = Par[3][:dt];
    f  = Par[4][:f];
    f1 = f[1];
    f2 = f[2];
    f3 = f[3];
    f4 = f[4];
    np = length(p);
    nh = length(h);
    if itype == 1;
       m = In;
       (nt, np) = size(m);
       d = zeros(nt,nh);
       m = SeisBPFilter(m,dt,f1,f2,f3,f4);
    end;
    if itype == -1;
       d = In;
       (nt, nh) = size(d);
       m = zeros(nt,np);
    end;
    if cmp(Par[5][:transf],"parab")==0 ;
      hmax = maximum(abs.(h));
      h=(h ./hmax) .^2;
    end
    nfft = 2*nextpow(2,nt);
    ilow  = floor(Int,f[1]*dt*nfft)+1;
    ihigh = floor(Int,f[4]*dt*nfft)+1;
     A = h*p';

     if itype == 1;
         m = cat( m, zeros(nfft-size(m)[1], size(m)[2]),dims=1)
         M = fft(m,1);
         D = zeros(ComplexF64,nfft,nh);

        for ifreq = ilow:ihigh
             w = 2*pi*(ifreq-1)/(nfft*dt);
             L = exp.(im*w*A');
             x = M[ifreq,:];
             y = L'*x;
             D[ifreq,:] =y;
         end

         for k=Int(nfft/2+2):nfft
             D[k,:] = conj(D[nfft-k+2,:]);
         end

         d = ifft(D,1);
         d = real(d[1:nt,:]);

         Out = d;

     else
         d = cat(d,zeros(nfft-size(d)[1], size(d)[2]),dims=1)
         D = fft(d,1);
         M = zeros(ComplexF64,nfft,np);

         for ifreq = ilow:ihigh
             w = 2*pi*(ifreq-1)/(nfft*dt);
             L = exp.(im*w*A');
             y = D[ifreq,:];
             x = L*y;
             M[ifreq,:]=x;
         end

         for k = Int(nfft/2+2):nfft
             M[k,:] = conj(M[nfft-k+2,:]);
         end

         m = ifft(M,1);
         m = real( m[1:nt,:]);

         Out = m;

     end
     if itype == -1;
         m = SeisBPFilter(m,dt,f1,f2,f3,f4);
         Out = m;
     end;
     return Out
  end
