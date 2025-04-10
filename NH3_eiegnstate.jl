#This is a code to calculate the eigenstates of a 1-Dimensional box
using QuantumOptics
using Plots

#Mass of particle and length of box
m = 1 #mass
l = 60 #length
n = 1 #Quantum number

#Defining a Basis to represent everything
xmin = -l
xmax = l
Npoints = 1000
b_position = PositionBasis(xmin, xmax, Npoints) #Position Basis
b_momentum = MomentumBasis(b_position) #Momentum Basis

#Defining Operators 
x = position(b_position) #Position Operator in Position Basis
p = momentum(b_momentum) #Momentum Operator in Momentum Basis
px = momentum(b_position) #Momentum Operator in Position Basis

#Defining the Potential 
function V(x)
    #v = 0.0105567-3.5627*10^(-5)*x^(2)+3.12338*10^(-8)*x^(4)
    v = 0.0105567-3.5627*10^(-5)*x^(2)+3.12338*10^(-8)*x^(4)
    return v
end

#Plotting the Potential
PE = potentialoperator(b_position, V)
PE = dense(PE)
ptsx = samplepoints(b_position)
#plot(ptsx , V)


#Defining the Hamiltonian
H = LazySum(px^2/2m , PE)
# A more feasible way is to use Fourier transforms(FT), so defining FT operators
Txp = transform(b_position, b_momentum)
Tpx = transform(b_momentum, b_position)
Hkin = LazyProduct(Txp, p^2/2m, Tpx) #A faster method to do Hkin * psi
H = dense(LazySum(Hkin, PE))
n= 1
#finding eigenstates
E,ψ = eigenstates(H) #Computes all  eigenstates
En,ψn = eigenstates(H,n) #Computes only the required eigenstate
ψn
E[4]-E[3]
E[2]-E[1]
E[6]-E[5]
#Plotting the States
#plot!(ptsx, 500*real(ψ[n].data)) #Can be plotted from all eigenstates 
n=1
Ep = [0.1*(ψ[n].data[i]) + real(E[n]) for i in 1:length(ψ[n])]
plot(ptsx, [real(Ep) , V], title="Ammonia v=0 Eigenstate", label=["v=0" "V(x)"])
#plot!(ptsx,V)
n=2
Ep = [0.1*(ψ[n].data[i]) + real(E[n]) for i in 1:length(ψ[n])]
plot!(ptsx, real(Ep), label = "v=1")

xlabel!("x")
ylabel!("ψ(x)")
savefig("NH3_Eigenstates_psi_v=0.png")
print("The Energy of the corresponding state is  : $(real(E[n])) ")
