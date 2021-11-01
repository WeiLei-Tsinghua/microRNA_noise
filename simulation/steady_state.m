function  sol = steady_state(type, para)
    para = num2cell(para);
    [KR, GR, KT, GT, Kon, Koff, Alpha, GC, Kp, Gp, ...
         KT2, GT2, Kon2, Koff2, Alpha2, GC2, ...
         KonDouble, KoffDouble, KR2, GR2, ...
         KonB, KoffB, KonA_CB, KoffA_CB, KonB_CA, KoffB_CA, ...
         GCB, AlphaB, GCAB, AlphaAB] = deal(para{:});


    function ds = eqn(t, x)
        x = num2cell(x);
        [R, T, C, p, T2, C2, CB, CAB, R2] = deal(x{:});
        dR = KR - GR*R - Kon*R*T + Koff*C + GC*(1 - Alpha)*C;
        if type == 1
            dR = dR - Kon2*R*T2 + Koff2*C2 + GC2*(1 - Alpha2)*C2;
        elseif type == 2
            dR = dR - Kon*R*T + Koff*CB - KonDouble*R*(C + CB) + 2*KoffDouble*CAB ...
                + GCB*(1 - AlphaB)*CB + GCAB*(1-AlphaAB)*CAB*2;
        elseif type == 3
            dR = dR - KonA_CB*R*CB + KoffA_CB*CAB + GCAB * (1 - AlphaAB) * CAB;
        end
        
        dT = KT - GT*T - Kon*R*T + Koff*C;
        if type == 2
            dT = dT - Kon*R*T + Koff*CB;
        elseif type == 3
            dT = dT - KonB*R2*T + KoffB*CB;
        end
        
        dC = Kon*R*T - Koff*C - GC*C;
        if type == 2
            dC = dC - KonDouble*R*C + KoffDouble*CAB;
        elseif type == 3
            dC = dC - KonB_CA*R2*C + KoffB_CA*CAB;
        end
        
        dp = Kp*T - Gp*p;
        
        dT2 = 0;
        dC2 = 0;
        if type == 1
            dT2 = KT2 - GT2*T2 - Kon2*R*T2 + Koff2*C2;
            dC2 = Kon2*R*T2 - Koff2*C2 -GC2*C2;
        end
        
        dCB = 0;
        if type == 2
            dCB = Kon*R*T - Koff*CB - KonDouble*R*CB + KoffDouble*CAB - GCB*CB;
        elseif type == 3
            dCB = KonB*R2*T - KoffB*CB - KonA_CB*R*CB + KoffA_CB*CAB - GCB*CB;
        end
        
        dCAB = 0;
        if type == 2
            dCAB = KonDouble*R*(C + CB) - 2*KoffDouble*CAB - GCAB*CAB;
        elseif type == 3
            dCAB = KonA_CB*R*CB + KonB_CA*R2*C - (KoffA_CB + KoffB_CA)*CAB - GCAB*CAB;
        end
        
        dR2 = 0;
        if type == 3
            dR2 = KR2 - GR2*R2 - KonB*R2*T + KoffB*CB - KonB_CA*R2*C + KoffB_CA*CAB ...
                + (1 - AlphaB)*GCB*CB + (1 - AlphaAB)*GCAB*CAB;
        end
          
        ds=[dR, dT, dC, dp, dT2, dC2, dCB, dCAB, dR2]';
    end


    ODEFUN=@eqn;
    x0=zeros(9,1);
    Tend=1e12;
    [t,species]=ode23s(ODEFUN,[0,Tend],x0);
    sol = species(end,:);

end