#################################### 
####   S-D-S: ⟨Sz_dot⟩ vs Γ      ####
####################################

#import the S-D-S solver function
include("solve_SDS.jl")

using DelimitedFiles #for writing the output data to disk

let
  π=4.0*atan(1.0)

  #SC gap Δ = Δ_L = Δ_R
  Δ=1.0
  
  #phase bias
  ϕ=0.0

  #Coulomb interaction strength
  U=5.0

  #dot level position
  ϵ_d=-U/2.0  # ph symmetric point here

  #number of effective levels
  L=2


  #grid for Γ scan
  n=11
  Γ= [(i-1)*10/(n-1) for i in 1:n]
  Γ[1]=0.001

  #loop over the parity Π=even/odd
  for Π in 0:1

    Sz_dot=fill(0.0,n)
    energy=fill(0.0,n)
    
    
    #do the Γ scan in parallel
    Threads.@threads for i in 1:n
      energy[i],Sz_dot[i]=solve_S_D_S(Δ,ϕ,L,ϵ_d,U,Γ[i],Π)
    end
  
     
    # write the data
      open("SDS_Szdot_L="*string(L)*"_parity="*string(Π)*".txt", "w") do g1  # "w" for writing, "a" for appending
      open("SDS_energy_L="*string(L)*"_parity="*string(Π)*".txt", "w") do g2  # "w" for writing, "a" for appending
        for i in 1:n
          write(g1,"$(Γ[i]) $(Sz_dot[i])"); write(g2,"$(Γ[i]) $(energy[i])")
          if(i<n)
            write(g1,"\n"); write(g2,"\n")
          end
        end
      end
      end
    end

end


