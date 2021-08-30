function SeisBPFilter(d,dt,f1,f2,f3,f4)
    Ndims = ndims(d);
    N     = size(d);
    nt    = N[1];
    x     = reshape(d,nt,prod(N[2:end]));
    nf=4*nextpow(2,nt)

    i1 = floor(Int,nf*f1*dt)+1;
    i2 = floor(Int,nf*f2*dt)+1;
    i3 = floor(Int,nf*f3*dt)+1;
    i4 = floor(Int,nf*f4*dt)+1;
    up=collect(1:1:(i2-i1))/(i2-i1);
    down=collect((i4-i3):-1:1)/(i4-i3);
    aux1 = append!(vec(zeros(1,i1)),up,vec(ones(1,i3-i2)), down,vec(zeros(1,Int(nf/2)+1-i4)));
    aux2 = reverse(aux1[2:Int(nf/2)]);
    c = 0;
    F=append!(aux1,aux2);

    Phase = (pi/180.)*append!([0],-c*ones(1,Int(nf/2-1)),[0],c*ones(1,Int(nf/2-1)));
    Transfer = F.*exp.(-im .*Phase);

    xp=zeros(nf,size(x)[2]);
    xp[1:size(x)[1],:]=x;
    X = fft(xp,1);
    Y=zeros(ComplexF64,size(X));
    for k = 1:prod(N[2:end])
         Y[:,k] = Transfer .*X[:,k];
    end

    o = ifft(Y,1);

    o = real(o[1:nt,:]);

    o = reshape(o, (nt, N[end]));
    return o
end
