# Authors: Matthew J Simpson and Pascal R Buenzli, 2023-2024

# compile and run with
# 	julia [-i] spring-diffusion-continuum.jl
# or, from the interpreter: include("spring-diffusion-continuum.jl")

using DifferentialEquations

# spring/cell model parameters
R = 1 # circle radius 1
L=2*π*R
N = 8 # number of cells (to set up initial density)
k1 = 1 # cell spring constant k*
eta1 = 1 # cell drag coefficient eta*
a1 = L/N # cell resting length a*

# discretisation parameters
dt = 0.001
t_end = 2
dt_frame = 0.001
incr_frame = 0
every_nth_dt = round(Int, dt_frame/dt)
dt_frame = every_nth_dt*dt
T = 0:dt:t_end
N_frames = length(T)
ds=L/1000
Ns=round(Int, L/ds)+1 # discretisation of interface
ds=L/Ns
s=LinRange(-L/2,L/2,Ns)

datadir="../datadir-continuum"

###### choose restoring force law #####
# Hookean springs / nonlinear diffusion
function Diffusivity(eta1,k1,a1,q)
    return (k1/eta1)/q^2
end

# # nonlinear springs / linear diffusion
# function Diffusivity(eta1,k1,a1,q)
#     return k1*a1^2/eta1
# end

# # nonlinear springs 2 / porous medium diffusion
# function Diffusivity(eta1,k1,a1,q)
#     return (k1*a1^3/eta1)*q
# end
########################################


function diff!(dq,q,p,t)
    dx,Ns,eta1,k1=p
    for i in 1:Ns
        # periodic BCs
        i==1 ? im1=Ns : im1=i-1
        i==Ns ? ip1=1 : ip1=i+1

        Di = Diffusivity(eta1,k1,a1,q[i])
        Dip = Diffusivity(eta1,k1,a1,q[ip1])
        Dim = Diffusivity(eta1,k1,a1,q[im1])
        dq[i]=((Di+Dip)*(q[ip1]-q[i])-(Di+Dim)*(q[i]-q[im1]))/(2*dx^2)
    end
end

function pdesolver(L,dx,Ns,T,q0,eta1,k1)
    p=(dx,Ns,eta1,k1)
    tspan=(0.0,maximum(T))
    prob=ODEProblem(diff!,q0,tspan,p)
    sol=solve(prob,saveat=T);

    alg = sol.alg;
    println("Algorithm used to solve the ODEs: $alg");
    # println("Algorithm stats: $sol.stats");
    # println("Algorithm stats: $sol.DEStats");
    
    for i in 1:length(sol[:,])
        qc[i,:]=sol[:,i]
    end
    
    return qc
end



# save metadata
mkpath(datadir)
open(joinpath(datadir, "metadata.dat"), "w") do metadata
    println(metadata, "dt=$dt")
    println(metadata, "t_end=$t_end")
    println(metadata, "dt_frame=$dt_frame")
    println(metadata, "every_nth_dt=$every_nth_dt")
    println(metadata, "N_frames=$N_frames")
    println(metadata, "# T=$T")
    println(metadata, "L=$L")
    println(metadata, "ds=$ds")
    println(metadata, "Ns=$Ns")
    println(metadata, "k1=$k1")
    println(metadata, "a1=$a1")
    println(metadata, "eta1=$eta1")
end

# initialisation:

# Circle, curved springs, N=8, m=8, k=8, a=0.0981748, eta=0.125 (post-scaling with m)
# NB: scaling with m does not affect continuum model, but it is used to setup the discrete model
q0 = ones(Ns)/(L/N) # uniform density
i0=round(Int, 2.35619/ds) # low index of extended spring
i1=round(Int, 3.53429/ds) # boundary extended / compressed spring
i2=round(Int, 3.92699/ds) # upper index compressed spring
for i in i0:i1 # for s in [2.35619, 3.53429]
    q0[i]= 1.0/(8*0.147129)
end
for i in i1+1:i2 # for s in [3.53429, 3.92699] 
    q0[i]=1.0/(8*0.0490825)
end

# # Circle, straight springs
# q0 = ones(Ns)/(2*R*sin(π/N)) # uniform density regular polygon edge

qc=zeros(length(T),Ns)

# solve
qc=pdesolver(L,ds,Ns,T,q0,eta1,k1)

# save data 
for t = 1:every_nth_dt:length(T)
    global incr_frame
    open(joinpath(datadir, "data_$incr_frame.dat"), "w") do fp
        println(fp, "# i, s_i, q(s_i, t)")
        for i in 1:Ns
            println(fp, i, "\t", i*ds, "\t", qc[t,i])
        end
        incr_frame = incr_frame + 1
    end
end
