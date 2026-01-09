using LinearAlgebra
include("./lib-robotique.jl"); 

function zMouvement(θ,robot,pd,Ż)
    global T = MGD(θ,robot);
    global p0 = T[1:3,4];
    global dt = 0.05;

    vec_erreur = pd - p0;
    v = Ż + vec_erreur;
    J = Jacobian(θ,robot,p0);
    # Jz ici est J(1:3,1:7)
    q̇ = pinv(J[1:3,1:7])*v;
    qc = θ + q̇*dt;
    return qc;
end

function eloignementButees(θ,robot,pd,Ż,θmin,θmax)
    global T = MGD(θ,robot);
    global p0 = T[1:3,4];
    global dt = 0.05;
    α = -2;

    vec_erreur = pd - p0;
    v = Ż + vec_erreur;
    J = Jacobian(θ,robot,p0);
    # Calcule ∇q
    θmoy = zeros(7);
    ∇Φ = zeros(7);
    for i=1:7
        θmoy = (θmin[i]+θmax[i])/2; # Toujours égale à 0
        ∇Φ[i] = 2*(θ[i]-θmoy) / (θmax[i]-θmin[i])^2;
    end
    # Jz ici est J(1:3,1:7)
    q̇ = pinv(J[1:3,1:7])*v + (I - pinv(J[1:3,1:7])*J[1:3,1:7])*α*∇Φ;
    qc = θ + q̇*dt;
    return qc;
end

function chapeau(S)
    chap = [0 -S[3] S[2];
            S[3] 0 -S[1];
            -S[2] S[1] 0];
    return chap;
end

function MCI(θ,robot,pd,pPointd,adesiree,ωd)
    Kp = 1.0;
    K0 = 1.0;

    global T = MGD(θ,robot);
    global p0 = T[1:3,4];
    Ae = T[1:3,1:3];
    global dt = 0.05;

    # Position
    vec_erreur = pd .- p0;
    v = pPointd + Kp*vec_erreur;
    J = Jacobian(θ,robot,p0);

    # Orientation
    global a = adesiree*(Ae'); # matrice de rotation
    # Erreur d'orientation
    global ϵ0 = 0.5*(cross(Ae[1:3,1],adesiree[1:3,1])+cross(Ae[1:3,2],adesiree[1:3,2])+cross(Ae[1:3,3],adesiree[1:3,3]))
    if (norm(ϵ0)<0.01) # Check si orientation ok
        ωe = ωd; # Pas de correction
    else
        se = chapeau(Ae[1:3,1]);
        sd = chapeau(adesiree[1:3,1]);
        ne = chapeau(Ae[1:3,2]);
        nd = chapeau(adesiree[1:3,2]);
        ae = chapeau(Ae[1:3,3]);
        ad = chapeau(adesiree[1:3,3]);
        L = -0.5*(se*sd+ne*nd+ae*ad);
        ωe = K0*ϵ0 + (L')*ωd;
    end

    q̇ = pinv(J)*[v;ωe];
    qc = θ + q̇*dt;
    return qc;
end

function MCI_eloignementButees(θ,robot,pd,pPointd,adesiree,ωd,θmin,θmax)
    Kp = 1.0;
    K0 = 1.0;
    α = -1.0;

    global T = MGD(θ,robot);
    global p0 = T[1:3,4];
    Ae = T[1:3,1:3];
    global dt = 0.05;

    # Position
    vec_erreur = pd .- p0;
    v = pPointd + Kp*vec_erreur;
    J = Jacobian(θ,robot,p0);

    # Orientation
    global a = adesiree*(Ae'); # matrice de rotation
    # Erreur d'orientation
    global ϵ0 = 0.5*(cross(Ae[1:3,1],adesiree[1:3,1])+cross(Ae[1:3,2],adesiree[1:3,2])+cross(Ae[1:3,3],adesiree[1:3,3]))
    if (norm(ϵ0)<0.01) # Check si orientation ok
        ωe = ωd; # Pas de correction
    else
        se = chapeau(Ae[1:3,1]);
        sd = chapeau(adesiree[1:3,1]);
        ne = chapeau(Ae[1:3,2]);
        nd = chapeau(adesiree[1:3,2]);
        ae = chapeau(Ae[1:3,3]);
        ad = chapeau(adesiree[1:3,3]);
        L = -0.5*(se*sd+ne*nd+ae*ad);
        ωe = K0*ϵ0 + (L')*ωd;
    end

    # Calcule ∇q eloignement butees
    θmoy = zeros(7);
    ∇Φ = zeros(7);
    for i=1:7
        θmoy = (θmin[i]+θmax[i])/2; # Toujours égale à 0
        ∇Φ[i] = 2*(θ[i]-θmoy) / (θmax[i]-θmin[i])^2;
    end

    q̇ = pinv(J)*[v;ωe] + (I - pinv(J[1:3,1:7])*J[1:3,1:7])*α*∇Φ;
    qc = θ + q̇*dt;
    return qc;
end