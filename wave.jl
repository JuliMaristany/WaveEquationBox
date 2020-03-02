using Plots

#Useful functions

function hat(i, N, x)
    @assert 0 <= i < N
    h = 1 / (N-1) #Domain then goes from 0 to 1
    xim1 = h * (i-1)
    xi = h * i
    xip1 = h * (i+1)
    if x <= xim1
        return 0.
    elseif x <= xi
        return (x - xim1) / h
    elseif x <= xip1
        return -(x - xip1) / h
    else
        return 0.
    end
end

function sample(f, N)
    h = 1 / (N-1)
    return coeffs = Float64[f(h*i) for i in 0:N-1]
end

function evaluate(coeffs, x)
    N = length(coeffs)
    sum(coeffs[i+1] * hat(i, N, x) for i in 0:N-1)
end

function deriv(f)
    N = length(f)
    h = 1 / (N-1)
    deriv_coeffs = similar(f)
    for i in 1:N
        if i == 1
            deriv_coeffs[i] = (f[2] - f[1]) / h
        elseif i == N
            deriv_coeffs[i] = (f[N] - f[N-1]) / h
        else
            deriv_coeffs[i] = (f[i+1] - f[i-1]) / (2*h)
        end
    end
    return deriv_coeffs
end

function deriv2(f)
    N = length(f)
    h = 1 / (N-1)
    deriv2_coeffs = similar(f)
    for i in 1:N
        if i == 1
            deriv2_coeffs[i] = (f[3] - 2*f[2] + f[1]) / h^2
        elseif i == N
            deriv2_coeffs[i] = (f[N] - 2*f[N-1] + f[N-2]) / h^2
        else
            deriv2_coeffs[i] = (f[i+1] - 2*f[i] + f[i-1]) / h^2
        end
    end
    return deriv2_coeffs
end

function rk2step(f, y0, h)
    k0 = f(y0)
    y1 = y0 + (h/2)*k0
    k1 = f(y1)
    y2 = y0 + h*k1
    y2
end

function rk2integrate(f, y0, h, n)
    res = [y0]
    for i in 1:n
        y1 = rk2step(f, y0, h)
        push!(res, y1)
        y0 = y1
    end
    res
end

#Solving now the wave equation in a box.
#My state vector will be an array where the first index
#indicates position and the second index indicates time step

L=100
T=100

#Create the array

u = zeros(L,T)

#Boundary Conditions

for n=1:T
   u[1,n]=0.
   u[L,n]=0.
end

#Initial Conditions

f(x)=sin(2*pi*x)

#Discretization Parameter.
#It has a big effect on stability

s = 1

#Finite difference scheme

fi=sample(f,L)
u[:,1]=fi

for j=2:L-1
  u[j,2]=u[j,1]-1/2*s*(u[j+1,1]-2*u[j,1]+u[j-1,1])
end

for n=2:T-1
    for j=2:L-1
        u[j,n+1]=s*(u[j+1,n]+u[j-1,n]) + 2*(1.0-s)*u[j,n]-u[j,n-1]
    end
end

#Plot the solution

g(x,j)=evaluate(u[:,j],x)

u5(x)=g(x,5)
u10(x)=g(x,10)
u25(x)=g(x,25)
u27(x)=g(x,27)
u30(x)=g(x,30)
u50(x)=g(x,50)

labels = ["t=5" "t=10" "t=25" "t=27" "t=30" "t=50"]

plot([u5,u10,u25,u27,u30,u50],0,1,xlabel="x",ylabel="y",plot_title=`Hello`,label=labels)

#Energy of the solution (Conserved)- Using to check for stability

ud=similar(u)
udt=similar(u)

for n=1:T
   ud[:,n]=deriv(u[:,n])
end

for i=1:L
   udt[i,:]=deriv(u[i,:])
end

h(x,j)=evaluate(ud[:,j],x)
k(x,j)=evaluate(udt[:,j],x)

E=[]
ind=[]

for j=1:T

     p(x)=(h(x,j)^2+k(x,j)^2)

     Energy=rk2integrate(p,0.001,0.001,9000)[9000]

     push!(E,Energy)

     push!(ind,j)

end

plot(ind,E,ylabel="Energy",xlabel="Time step",label="Energy")

## Check for convergence

E=[]
ind=[]

for i=1:10

    T=10*i

    L=10*i

    u = zeros(L,T)

    for n=1:T
       u[1,n]=0.
       u[L,n]=0.
    end

    fi=sample(f,L)
    u[:,1]=fi

    for j=2:L-1
      u[j,2]=u[j,1]-1/2*s*(u[j+1,1]-2*u[j,1]+u[j-1,1])
    end

    for n=2:T-1
        for j=2:L-1
            u[j,n+1]=s*(u[j+1,n]+u[j-1,n]) + 2*(1.0-s)*u[j,n]-u[j,n-1]
        end
    end

    ud=similar(u)
    udt=similar(u)

    for n=1:T
       ud[:,n]=deriv(u[:,n])
    end

    for i=1:L
       udt[i,:]=deriv(u[i,:])
    end

    #Let's take a particular j, that is, a particular time step

    j=10

    g(x,j)=evaluate(u[:,j],x)
    h(x,j)=evaluate(ud[:,j],x)
    k(x,j)=evaluate(udt[:,j],x)

    p(x)=(h(x,j)^2+k(x,j)^2)

    Energy=rk2integrate(p,0.001,0.001,9000)[9000]

    push!(E,Energy)

    push!(ind,L)

end

plot(ind,ylabel="Energy",xlabel="Resolution",label="Energy")
