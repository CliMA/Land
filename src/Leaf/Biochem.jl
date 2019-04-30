
"""
    Calculates net assimilation rate A, fluorescence F using biochemical model
"""
function biochem_out = biochemical(flux::fluxes, leaf::leaf_params,T,Cs)
    # Adjust rates to leaf Temperature:
    setLeafT!(leaf, T)

    # Calculate potential electron transport rate:
    Je_pot = 0.5*leaf.po0 * flux.APAR; # potential electron transport rate

    # Actual Max Je given Jmax
    Je = minimum(quadratic(leaf.θ_j, -(Je_pot + leaf.jmax), Je_pot * leaf.jmax))
    # C3: Rubisco-limited photosynthesis
    ac = flux.vcmax * max(ci_val - flux.cp, 0) / (ci_val + flux.kc * (1 + atmos.o2air / flux.ko));

    # C3: RuBP regeneration-limited photosynthesis
    flux.aj = flux.je * max(ci_val - flux.cp, 0) / (4 * ci_val + 8 * flux.cp);



# Compute Assimilation.
function computeA!(Ci, leaf::leaf_params, bio::bioOut)
    if leaf.C3
        wc = (leaf.Vcmax .* (Ci - Gamma_star))./(Ci + MM_const); % This is Rubisco Limiting Step

            % Compute the Potential Electron Transport Rate
            thetaPSII = 0.7;            % Curvature parameter
            f = 0.15;
            phiPSII = (1-f);            % Quantum Yeild of PSII
            IPSII = 0.5*phiPSII.*Q;     % light utilized in electron transport by PSII

            % Je is the potential electron transport rate computed as the
            % smaller root of the quadratic
            Je = sel_root(thetaPSII, -(IPSII+Jmax), IPSII.*Jmax, -1);

            wj = Je.*(Ci-Gamma_star)./(4.*Ci + 8.*Gamma_star);      % This is light limiting step

            ws = 3*TPU; % This is the product limiting step (doesn't change with iteration of A-Ci)

            thetacj = 0.98;
            thetais = 0.95;

            CO2_per_electron = (Ci-Gamma_star)./(Ci + 2.*Gamma_star) .* effcon; % same as previous version


        else % Compute A for C4

            wc = Vcmax;                 % Rubisco Limiting Step
            alpha = 0.05;               % Quantum Yeild [mol mol^-1]
            wj = alpha .* Q;            % Light limiting rate Q is already in mumols m-2 s-1
            ws = Ke .* Ci * 1E-6;       % Product limitation

            thetacj = 0.80;
            thetais = 0.95;
            CO2_per_electron = effcon; % same as previous version


            % Compute Je for C4
            f = 0.15;
            phiPSII = (1-f);            % Quantum Yeild of PSII
            IPSII = 0.5*phiPSII.*Q;     % light utilized in electron transport by PSII
            Je=IPSII;

        end


        % find the smoothed minimum of we, wc and then ws

        wi      = sel_root(thetacj,-(wc+wj),wc.*wj, -1 ); % i.e. sign(Gamma_star - Ci)
        Ag      = sel_root(thetais,-(wi+ws),wi.*ws, -1);
        A       = Ag - Rd;


        % Compute the Photorespiration (rate of oxygenation vo) from net
        % assimilation
        % The formulations are provided in Sharkey 1988 and Bernachhi et al
        % 2001 Papers
        if strcmpi('C3', Type)
            phi = 2.*Gamma_star./Ci;
            vo = (A + Rd)./(1./phi - 0.5);
        end


        if nargout > 1
            A_out.A = A;
            A_out.Ag = Ag;
            A_out.Vc = wc;
            A_out.Vs = ws;
            A_out.Ve = wj;
            A_out.Je = Je;
            A_out.CO2_per_electron = CO2_per_electron;
            if strcmpi('C3', Type)
                A_out.vo = vo;
            end
        end
        fcount = fcount + 1; % # of times we called computeA

        end

struct flux
         APAR
         gbc
         ac
         aj
         ap
         ag
         an
         cs
         ci
end

struct atmos
         o2air              # Atmospheric O2 (mmol/mol)
         co2air             # Atmospheric CO2 (μmol/mol)
         eair               # Vapor pressure of air (Pa)
end
