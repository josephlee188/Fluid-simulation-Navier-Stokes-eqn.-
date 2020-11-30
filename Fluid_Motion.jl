function diffuse!(b, x, x0, diff::Float64, dt::Float64)
    a = dt * diff * (N-2)^2
    lin_solve!(b, x, x0, a, 1+4*a)
    set_bnd!(b, x)
end

function lin_solve!(b, x, x0, a, c)
    # Via Gauss-Seidel method, this function
    # solves linear system of 'Ax = b' with 'A' in a specific form; 
    # matrix of Laplacian, ∇² 
    @inbounds for k in 1:iter

        @. x[2:end-1,2:end-1] = (x0[2:end-1,2:end-1] + a * (x[3:end,2:end-1] +
                                                            x[1:end-2,2:end-1] +
                                                            x[2:end-1,1:end-2] + 
                                                            x[2:end-1,3:end])) / c
        set_bnd!(b, x)
    end
end

function project!(velx, vely, p, div)

    @. div[2:end-1,2:end-1] = -0.5 * (velx[3:end,2:end-1] - velx[1:end-2,2:end-1] + 
                                    vely[2:end-1,3:end] - vely[2:end-1,1:end-2]) / N
    @. p[2:end-1,2:end-1] = 0

    set_bnd!(0, div)
    set_bnd!(0, p)
    lin_solve!(0, p, div, 1, 4)

    @. velx[2:end-1,2:end-1] -= 0.5 * (p[3:end,2:end-1] - p[1:end-2,2:end-1]) * N
    @. vely[2:end-1,2:end-1] -= 0.5 * (p[2:end-1,3:end] - p[2:end-1,1:end-2]) * N

    set_bnd!(1, velx)
    set_bnd!(2, vely)
end

inflate(f, is, js) = [f(i,j) for i in is, j in js]

function rev!(d, dtx, dty, velx, vely, d0, aux)

    Cx = @view aux[:,:,1]
    Cy = @view aux[:,:,2]
    CI0 = @view aux[:,:,3]
    CI1 = @view aux[:,:,4]
    CJ0 = @view aux[:,:,5]
    CJ1 = @view aux[:,:,6]
    CS0 = @view aux[:,:,7]
    CS1 = @view aux[:,:,8]
    CT0 = @view aux[:,:,9]
    CT1 = @view aux[:,:,10]
    mgridx = @view aux[:,:,11]
    mgridy = @view aux[:,:,12]

    @. Cx = mgridx - dtx * velx[2:N-1,2:N-1]
    @. Cy = mgridy - dty * vely[2:N-1,2:N-1]
    
    Cx[findall(Cx .< 1.5)] .= 1.5
    Cy[findall(Cy .< 1.5)] .= 1.5
    Cx[findall(Cx .> N-0.5)] .= N-0.5
    Cy[findall(Cy .> N-0.5)] .= N-0.5

    @. CI0 = floor(Cx)
    @. CI1 = CI0 + 1.0
    @. CJ0 = floor(Cy)
    @. CJ1 = CJ0 + 1.0

    @. CS1 = Cx - CI0
    @. CS0 = 1.0 - CS1
    @. CT1 = Cy - CJ0
    @. CT0 = 1.0 - CT1

    A = Int.(CI0 .+ N .* (CJ0 .- 1))
    B = Int.(CI0 .+ N .* (CJ1 .- 1))
    C = Int.(CI1 .+ N .* (CJ0 .- 1))
    D = Int.(CI1 .+ N .* (CJ1 .- 1))

    @. d[2:end-1,2:end-1] = 
    CS0 * CT0 * d0[A] + 
    CS0 * CT1 * d0[B] + 
    CS1 * CT0 * d0[C] + 
    CS1 * CT1 * d0[D]
end

function advect!(b, d, d0, velx, vely, dt, aux)
    dtx = dt * (N-2)
    dty = dt * (N-2)

    # Original version
    # for j in 2:N-1
    #     for i in 2:N-1
    #         x = i - dtx * velx[i,j]
    #         y = j - dty * vely[i,j]
    #         if x < 1.5 x = 1.5 end
    #         if y < 1.5 y = 1.5 end
    #         if x > N-0.5 x = N-0.5 end
    #         if y > N-0.5 y = N-0.5 end
    #         i0 = floor(x)
    #         i1 = i0 + 1.0
    #         j0 = floor(y)
    #         j1 = j0 + 1.0

    #         s1 = x - i0
    #         s0 = 1.0 - s1
    #         t1 = y - j0
    #         t0 = 1.0 - t1

    #         i0 = Int(i0)
    #         i1 = Int(i1)
    #         j0 = Int(j0)
    #         j1 = Int(j1)

    #         d[i,j] = 
    #         s0 * (t0 * d0[i0,j0] + t1 * d0[i0,j1]) + 
    #         s1 * (t0 * d0[i1,j0] + t1 * d0[i1,j1])
    #     end
    # end

    rev!(d, dtx, dty, velx, vely, d0, aux)

    set_bnd!(b, d)
end

function set_bnd!(b, x)

    a = b==2 ? -1 : 1
    @. x[2:end-1,1] = a * x[2:end-1,2]
    @. x[2:end-1,N] = a * x[2:end-1,N-1]

    a = b==1 ? -1 : 1
    @. x[1,2:end-1] = a * x[2,2:end-1]
    @. x[N,2:end-1] = a * x[N-1,2:end-1]

    x[bnd_ind1] .= 0.5 .*(x[bnd_ind2] .+ x[bnd_ind3])

end






