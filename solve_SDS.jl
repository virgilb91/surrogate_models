using ITensors
using LinearAlgebra


include("dmrg_obs.jl") #observer file for monitoring DMRG

function solve_S_D_S(Δ,ϕ,L,ϵ_d,U,Γ,Π)

  #import surrogate model parameters
  ξ=readdlm("xi_L="*string(L)*"_B=10_wc=10.dat",Float64)
  γ=readdlm("gamma_L="*string(L)*"_B=10_wc=10.dat",Float64)
  ξ=[(ξ...)...];  γ=[(γ...)...]
  p = sortperm(ξ); γ=γ[p]; ξ=ξ[p]

  #effective hopping amplitudes; note that Γ_L = Γ_R = Γ/2
  t=sqrt.(γ*Γ/2.0)   
  t₁=t  # hopping L
  t₂=t  # hopping R

  Δ₁=Δ*exp(im*ϕ/2)  # gap L
  Δ₂=Δ*exp(-im*ϕ/2) # gap R

  # site indices for the type Electron (spin-1/2 fermion) with parity and S_z conservation
  sites = siteinds("Electron", 2*L+1,conserve_nfparity=true,conserve_sz=true)

  #################################
  #####   S-D-S HAMILTONIAN   #####
  #################################
  ampo = OpSum()

  for i in 1:L
    ampo += ξ[i], "Ntot", i

    ampo += -Δ₁, "Cdagup", i, "Cdagdn",i
    ampo += -conj(Δ₁), "Cdn", i, "Cup", i

    ampo += t₁[i], "Cdagup", i, "Cup", L+1
    ampo += t₁[i], "Cdagdn", i, "Cdn", L+1
    ampo += conj(t₁[i]), "Cdagup", L+1, "Cup", i
    ampo += conj(t₁[i]), "Cdagdn", L+1, "Cdn", i
  end

  ampo += ϵ_d, "Ntot", L+1
  ampo += U, "Nup",L+1, "Ndn",L+1

  for i in L+2:2*L+1
    ii=i-L-1
    ampo += ξ[ii], "Ntot", i

    ampo += -Δ₂, "Cdagup", i, "Cdagdn",i
    ampo += -conj(Δ₂), "Cdn", i, "Cup", i

    ampo += t₂[ii], "Cdagup", i, "Cup", L+1
    ampo += t₂[ii], "Cdagdn", i, "Cdn", L+1
    ampo += conj(t₂[ii]), "Cdagup", L+1, "Cup", i
    ampo += conj(t₂[ii]), "Cdagdn", L+1, "Cdn", i
  end

  H = MPO(ampo, sites)


  #################################
  ######    DOT S_z OPERATOR  #####
  #################################
  
  ampo = OpSum()
  ampo += +1/2, "Nup",L+1
  ampo += -1/2, "Ndn",L+1
  Sz_dot = MPO(ampo, sites)
 


  #input state for dmrg: defines the parity and S_z sector
  state = ["Emp" for n in 1:2*L+1]
  if(Π==1)
    state[L+1]= "Up"
  end
  ψ₀₀ = randomMPS(sites,state)


  #dmrg parameters; see ITensor documentation for details
  cutoff = 1E-10

  maxdim1 = [100,200,300,400]; maxdim2 = [500 for j=1:100]
  maxdim = vcat(maxdim1,maxdim2)

  nsweeps = length(maxdim)
  noise1 = [1E-6,1E-7,1E-8]; noise2=  [1E-10 for j=1:100]
  noise = vcat(noise1,noise2)

  #energy tolerance: dmrg stops when |E_new-E_old|<etol
  etol = 1E-4
  obs = DemoObserver(etol)

  #DMRG call: outputs the lowest energy and associated wavefunction in the chosen parity and S_z sector
  E₀, ψ₀ = dmrg(H, ψ₀₀; nsweeps, cutoff, maxdim, noise, observer=obs, outputlevel=1)

  #return the energy and Sz_dot expectation value
  return E₀, real(inner(ψ₀,Sz_dot,ψ₀))


end