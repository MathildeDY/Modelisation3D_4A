#connexion avec VREP
PATHCSIM="./";
include("./lib-robotique.jl"); 
include("./lib-CSim.jl");
include("./fonctions.jl");

using Plots

# Variables
global it_max = 250;

global clientID=startsimulation(simx_opmode_oneshot) # On lance une instance de connexion avec VREP

if clientID==0 println("Connexion établie")
    else println("Erreur de connexion")
end

global rob=CreateRobotKukaLwr();

global angle = [0.143019, −0.109465, −0.011994, −1.1788, −0.154233, 0.93555, 0.264868];
global pos_init = [−0.3668; −0.0379; 0.8634];
global pos_desiree = [−0.3668; −0.0379; 0.5];
global vitesse_desiree = [0.0; 0.0; 0.0];
global orient_desiree = RTL2R(pi/12,pi/6,0);
global ωd = zeros(3);
global θmin = [-2.97, -2.09, -2.97, -2.09, -2.97, -2.09, -3.05];
global θmax = [ 2.97,  2.09,  2.97,  2.09,  2.97,  2.09,  3.05];
global erreur = 1;
global i = 1;

# Inits
global ϵ0 = 1;
mgd = MGD(angle,rob);
global a = mgd[1:3,1:3];
global p0 = [−0.3668; −0.0379; 0.8634];
global traj = zeros(7,it_max);
global traj_z = zeros(it_max);
global X = zeros(it_max);
global Y = zeros(it_max);
global Z = zeros(it_max);
global error_orient = zeros(it_max);

setjointposition(clientID,angle,7,0,objectname_kuka);

while (norm(p0 - pos_desiree) > 0.001 || norm(ϵ0) > 0.001) && (i < it_max)
    # q = zMouvement(angle,rob,pos_desiree,vitesse_desiree);
    # q = eloignementButees(angle,rob,pos_desiree,vitesse_desiree,θmin,θmax);
    # q = MCI(angle,rob,pos_desiree,vitesse_desiree,orient_desiree,ωd);
    q = MCI_eloignementButees(angle,rob,pos_desiree,vitesse_desiree,orient_desiree,ωd,θmin,θmax)

    setjointposition(clientID,q,7,0,objectname_kuka);
    sleep(0.1)
    # Variables pour les plots
    global traj[:,i] = vec(q); 
    global traj_z[i] = p0[3,1];
    global X[i] = p0[1]
    global Y[i] = p0[2]
    global Z[i] = p0[3]
    global error_orient[i] = norm(ϵ0);
    # Variables de boucle
    global i = i+1;
    # Mise à jour variable globale
    global angle = q;
end

mgd = MGD(angle,rob);
global a = mgd[1:3,1:3];

stopsimulation(clientID,simx_opmode_oneshot) # Arrêt de la simulation

# Plots
# Articulations
p1 = plot(
    traj[1,1:i-1], label="q1", xlabel="Itération", ylabel="Angle (rad)", title="Mouvement des articulations")

for j = 2:7
    plot!(traj[j,1:i-1], label="q$j")
end

display(p1)

# Mouvement effecteur Z
p2 = plot(
    traj_z[1:i-1],
    xlabel="Itération",
    ylabel="Position z (m)",
    title="Mouvement de l’effecteur selon l’axe z",
    legend=false
)

display(p2)

# Mouvement effecteur 3D
p3 = plot(
    X[1:i-1],
    Y[1:i-1],
    Z[1:i-1],
    xlabel = "X (m)",
    ylabel = "Y (m)",
    zlabel = "Z (m)",
    title = "Trajectoire 3D de l’effecteur",
    legend = false,
    lw = 2,
    aspect_ratio = :equal
)

display(p3)

# Erreur d’orientation
p4 = plot(
    error_orient[1:i-1],
    xlabel="Itération",
    ylabel="Erreur d’orientation",
    title="Évolution de l’erreur d’orientation",
    legend=false
)
