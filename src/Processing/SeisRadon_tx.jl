"""
    SeisRadon_tx(In,Par,itype);
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
- `Par.transf `:'parab', 'linear', 'hyperb'
itype =  1  :  forward
itype = -1  :  adjoint

"""


function SeisRadon_tx(In,Par,itype);
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
   tmax = (nt-1)*dt;
   tau = collect(0:1:nt-1)*dt;


     for itau = 1:1:nt
       for ih = 1:nh
        for ip = 1:np
           if cmp(Par[5][:transf],"parab")==0 ;
            hmax = maximum(abs.(h));
            t=tau[itau]+ (h[ih]/hmax)^2*p[ip];
         elseif cmp(Par[5][:transf],"linear")==0;
            t=tau[itau] + (h[ih]*p[ip]);
         elseif cmp(Par[5][:transf],"hyperb")==0;
            t=sqrt(tau[itau]^2 + (h[ih]^2/p[ip]^2));
         end
          if t>0 && t<tmax
            it    = t/dt+1;
            it1   = floor(Int,it);
            it2   = it1 + 1;
            alpha = it  - it1;
               if itype == 1;
                    d[it1,ih] = d[it1,ih]+ (1-alpha) * m[itau,ip];
                    d[it2,ih] = d[it2,ih] +     alpha * m[itau,ip];
                 else
                    m[itau,ip] = m[itau,ip] + (1-alpha) * d[it1,ih] + alpha*d[it2,ih];

          end
        end
       end
      end
   end
   if itype == 1;
      Out = d;
   end;

   if itype == -1;
      m = SeisBPFilter(m,dt,f1,f2,f3,f4);
      Out = m;
   end;
   return Out
end
