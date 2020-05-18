function ψ_m(ζ,stab_type_stable) # momentum correction function
  FT = eltype(ζ)

    if(ζ>=0) # stable conditions
        if(stab_type_stable==1) # Grachev default value - Grachev et al SHEBA 2007
          x  = (1+ζ)^FT(1/3);
          am = FT(5.0);
          bm = am/FT(6.5);
          ah = bh = FT(5.0);
          ch = FT(3.0);
          Bm = ((FT(1.0)-bm)/bm)^FT(1/3);
          Bh = sqrt(5);
          ψ = -FT(3)*am/bm*(x-1) + am*Bm/(2bm)*( FT(2)*log((x+Bm)/(1+Bm))
                                                   - log( (x*x-x*Bm+Bm*Bm)/(1-Bm+Bm*Bm) )
                                                   +2*sqrt(FT(3.0))*( atan((2x-Bm)/(sqrt(FT(3.0))*Bm)) -  atan((2-Bm)/(sqrt(FT(3.0))*Bm)))  );
        elseif(stab_type_stable==2) # Webb correction as in NoahMP, quite numerically unstable at night - lots of oscillations
          alpha = FT(5.0);
          """ should be the correct one but not continuously differntialble
          if(ζ>=0.0 && ζ<=1)
            ψ = -alpha * ζ;
          else
            ψ = -alpha; # should be z less anyways but here just made to create a plateau in the value
          end
          """
          # Pierre's new one:
          ψ = -alpha*tanh(ζ);
        elseif(stab_type_stable==3) # Beljaars-Holtslag, 1991, Journal of Applied Meteorology
          # Constants
          a = FT(1.0);
          b = FT(0.667);
          c = FT(5);
          d = FT(0.35);
          # Stability correction for momentum
          ψ = -(a * ζ
                   + b * (ζ - c / d) * exp(-d * ζ)
                   + b * c / d);
        elseif(stab_type_stable==4) #Cheng and Brutsaert (2005) as described by Jimenez et al. (2012)
          a = FT(6.1);
          b = FT(2.5);
          c = FT(5.3);
          d = FT(1.1);
          ψ = -a * log(ζ + (1 + ζ^b) ^ (1 / b));
        end

    else #  (ζ<0) unstable
        x   =   (1-16ζ)^FT(0.25);
        ψ   =   FT(2.0)*log((1+x)/2) + log((1+x*x)/2) - 2*atan(x) + FT(π)/2;
    end
end




function ψ_h(ζ,stab_type_stable) # momentum correction function

if(ζ>=0.0) # stable conditions
      if(stab_type_stable==1) # Grachev default value - Grachev et al SHEBA 2007
        # issue when x ~ 0
        x  = (1.0+ζ)^(1.0/3.0);
        am = 5.0;
        bm = am/6.5;
        ah = bh = 5.0;
        ch = 3.0;
        Bm = ((1.0-bm)/bm)^(1.0/3.0);
        Bh = sqrt(5);
        ψ = -bh/2*log(1+ch*ζ+ζ*ζ) + (-ah/Bh + bh*ch/(2.0*Bh))*(log((2.0*ζ+ch-Bh)/(2.0*ζ+ch+Bh)) - log((ch-Bh)/(ch+Bh)) );
      elseif(stab_type_stable==2) # Webb correction as in NoahMP
        alpha = 5.0;
        #if(ζ>=0.0 && ζ<=1)
        #  ψ = -alpha * ζ;
        #else
        #  ψ = -alpha;
        #end
        # Pierre - I didn't like this so I changed it to be C1
        ψ = -alpha*tanh(ζ);
      elseif(stab_type_stable==3) # Beljaars-Holtslag, 1991, Journal of Applied Meteorology
        # Constants
        a = 1.0;
        b = 0.667;
        c = 5;
        d = 0.35;
        # Stability correction for momentum
        ψ = (- (1. + (2. / 3.) * a * ζ) ^ (3. / 2.)
                  - b * (ζ - (c / d)) * exp(-d * ζ)
                  - (b * c) / d + 1);
      elseif(stab_type_stable==4) #Cheng and Brutsaert (2005) as described by Jimenez et al. (2012)
        a = 6.1;
        b = 2.5;
        c = 5.3;
        d = 1.1;
        ψ =  -c * log(ζ + (1 + ζ^d) ^ (1. / d));
      end

    else # (ζ<0) unstable
        x   = (1.0-16.0*ζ)^0.25;
        ψ = 2.0*log((1.0+x*x)/2.0);
    end

    return ψ
end