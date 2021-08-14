

module TwoDerrors

import Printf
import JSON
import Random
import Statistics
import JLD
import PyPlot

const C = 330

# Sensor positions
sensorA=(10,-5,0,"A")
sensorB=(40,-5,0,"B")
sensorC=(-5,45,-90,"C")
sensorsAB=(sensorA,sensorB)

# For constructing file names.
sensor2string(s) = Printf.@sprintf("%+04d%+04d%+04d",s[1],s[2],s[3])
sensors2string(s) = join(map(sensor2string,s))

# Angle from microphone array to a point on the modelled area.
theta(x,y,sensorPos) = atan.(x-sensorPos[1], y-sensorPos[2])

# Equation 8
function phi(theta, ifParams)
    f,d,C=ifParams
    2*pi*f*d.*(sin.(theta)./C)
end

# Equation 8
function phi(theta,f,d)
    return 2*pi*f*d/C * sin(deg2rad(theta))
end

function errorMapFileName(sensors,d,f) 
    fmt=Printf.@sprintf("jld/%s_%02d_%05d",sensors2string(sensors),
                        convert(Int64,d*1000),convert(Int64,f))
    return "errormap_$fmt.jld"
end

coeffsFileName(froot,f,d) = return Printf.@sprintf("json/%s_%05d_%02d_%02d.json",
                                                   froot,
                                                   Int(floor(f)),
                                                   Int(d*1000),
                                                   5)

# Reads pre-calculated model coefficients from json files.
function makeCoeffs(froot)
    dict=Dict()
    for d in 1:200
        for f in [2000,4000,6000,8000]#,4000,6000,8000]
            di=Dict()
            open(coeffsFileName(froot,f,d/1000.0), "r") do f
                dicttxt = read(f,String)  # file information to string
                di=JSON.parse(dicttxt)
            end
            dict[(f,d)]=di
        end
    end
    return dict
end

BiasCoeffs=makeCoeffs("biascoeffs")
ErrorBarCoeffs=makeCoeffs("polycoeffs")

# Equation 9
# Takes radians returns degrees
function biasModel(theta,f,d)
    c=BiasCoeffs[(f,d)]
    alpha=c["alpha"]
    beta=c["beta"]
    gamma=c["gamma"]
    model = alpha/(cos(gamma*theta)^beta) - alpha
    return model*sign(theta)
end

# Equation 10
# Takes radians and returns degrees
function errorBarModel(theta,f,d)
    c=ErrorBarCoeffs[(f,d)]["coeffs"]
    return c[1]*theta^4 + c[2]*theta^3 + c[3]*theta^2 + c[4]*theta + c[1]
end

# Purturb the angle
# Takes radians returns degrees
function newAngleRad(theta,f,d)
    bias=biasModel(theta,f,d)
    err=errorBarModel(theta,f,d)
    #println("Bias: ",bias," error: ",err)
    return rad2deg(theta)-bias + err*Random.rand() - err/2
end

# Purturb the angle
# Takes degrees returns degrees
function newAngle(theta,f,d)
    return newAngleRad(deg2rad(theta),f,d)
end

# Generates and saves an error map
function errorMapN(sensors;f=2000.0,d=0.020,M=1000)
    N=length(sensors)
    GZ=51
    ifparams=(f,d,330)
    posx=Array{Float64}(undef,GZ,GZ)
    posy=Array{Float64}(undef,GZ,GZ)
    errx=Array{Float64}(undef,GZ,GZ)
    erry=Array{Float64}(undef,GZ,GZ)
    th_=zeros(N)
    th=Array{Float64}(undef,N)
    phi_=Array{Float64}(undef,N)
    slopeErr=zeros(N)
    for x = 0:GZ-1
        for y = 0:GZ-1
            for i = 1:N
                th_[i] = theta(x,y,sensors[i]) 
            end
            xp=zeros(M)
            yp=zeros(M)
            for test = 1:M
                th=newAngle.(th_,Int(floor(f)),Int(floor(d*1000)))
                if false
                    println("x: ",x," y: ",y)
                    println("Angle:")
                    println(map(x -> Printf.@sprintf("%.1f",x), rad2deg.(th_)))
                    println("Error:")
                    println(map(x -> Printf.@sprintf("%.1f",x), slopeErr))
                    println("New angle:")
                    println(map(x -> Printf.@sprintf("%.1f",x), rad2deg.(th)))
                    println()
                end

                A=hcat(tan.(th),ones(N))
                b=[sensors[i][1]*tan(th[i]) + sensors[i][2] for i in 1:N]
                if N > 2
                    m=GLM.lm(A,b)
                    xp[test],yp[test]=GLM.coef(m)
                else
                    xp[test],yp[test]=A\b
                end
            end
            errx[x+1,y+1]=Statistics.std(xp)
            erry[x+1,y+1]=Statistics.std(yp)
            posx[x+1,y+1]=Statistics.mean(xp)
            posy[x+1,y+1]=Statistics.mean(yp)
        end
    end
    fname=errorMapFileName(sensors,d,f)
    JLD.save(fname,"errx",errx,"erry",erry,"posx",posx,"posy",posy)
end

# Load a pre-calculated map
function loadMap(sensors,f,d)
    fname=errorMapFileName(sensors,d,f)
    return JLD.load(fname)#,"errx",errx,"erry",erry,"posx",posx,"posy",posy)
end

function plotErrorMapN(sensors,f,d;colour=false,levelsx=false,levelsy=false,PLOTSENSORS=false,contcolours="none",cmap="gist_rainbow",vmin=0,vmax=15)
    println("cmap=",cmap)
    emap=loadMap(sensors,f,d)
    errx=emap["errx"]
    erry=emap["erry"]
    posx=emap["posx"]
    posy=emap["posy"]
    fmt=Printf.@sprintf("%s_%02d_%04d",sensors2string(sensors),
                        convert(Int64,d*1000),convert(Int64,f))

    fint=convert(Int64, f)
    fig1=PyPlot.figure(1)
    fig1.clf()
    if colour
        if levelsx==false
            PyPlot.contourf(errx',cmap=cmap,vmin=vmin,vmax=vmax)
        else
            PyPlot.contourf(errx',cmap=cmap,levels=levelsx,vmin=vmin,vmax=vmax)
        end
    end

    if colour
        if contcolours == "none"
            contcolours="white"
        end
    else
        contcolours="black"
    end
    
    if levelsx==false
        CS = PyPlot.contour(errx',colors=contcolours)
    else
        CS = PyPlot.contour(errx',colors=contcolours,levels=levelsx)
    end
    
    PyPlot.clabel(CS)
    PyPlot.xlabel("Distance (meters)")
    PyPlot.ylabel("Distance (meters)")
    PyPlot.savefig("mc_error_x_$fmt.png")
    
    fig2=PyPlot.figure(2)
    fig2.clf()
    if colour
        if levelsy==false
            PyPlot.contourf(erry',cmap=cmap)
        else
            PyPlot.contourf(erry',cmap=cmap,levels=levelsy)
        end
    end

    if levelsy==false
        CS = PyPlot.contour(erry',colors=contcolours)
    else
        CS = PyPlot.contour(erry',colors=contcolours,levels=levelsy)
    end

    if PLOTSENSORS
        for s=sensors
            x,y=s[1],s[2]
            PyPlot.plot(x,y,"ko")
            PyPlot.annotate(s[4],xy=(x+2,y))
        end
        PyPlot.ylim([-7,50])
        PyPlot.xlim([-7,50])
    end


    
    
    PyPlot.clabel(CS)
    PyPlot.xlabel("Distance (meters)")
    PyPlot.ylabel("Distance (meters)")
    PyPlot.savefig("mc_error_y_$fmt.png")

    (errx,erry,posx,posy)
end

function paperPlot(sensors,f,d;M=600,colour=false,contcolours="black",levs=[0,1,1.5,2,3,4,5,6,8,10,15],cmap="gist_rainbow",vmin=0,vmax=15,calculate=true)
    if calculate
        errorMapN(sensors,f=f,d=d,M=M)#,levelsx=mylevs,levelsy=mylevs)#,contcolours="black")
    end
    plotErrorMapN(sensors,f,d;colour=colour,levelsx=levs,levelsy=levs,vmin=vmin,vmax=vmax,cmap=cmap,contcolours=contcolours)
    fig=PyPlot.figure(1)
    froot="$(sensors[1][4])$(sensors[2][4])$(Int(f))-$(Int(d*1000))"
    println(froot)
    PyPlot.savefig(froot*"x.png")
    fig=PyPlot.figure(2)
    PyPlot.savefig(froot*"y.png")
end

# For interactive testing
function ppRun(f,d;M=1000)
    paperPlot(TwoDerrors.sensorsAB,f,d;M=M,levs=[0,1,1.5,2,3,4,5,6,8,10,15],colour=true,cmap="jet",calculate=true)
end

# Plots for paper
function makePlots(;M=1000)
    for (f,d) in zip((2000,8000),(0.008,0.012,0.020))
        paperPlot(TwoDerrors.sensorsAB,f,d;M=M,levs=[0,1,1.5,2,3,4,5,6,8,10,15],colour=true,cmap="jet",calculate=true)
    end
    paperPlot(TwoDerrors.sensorsAB,2000,0.060;M=M,levs=[0,1,1.5,2,3,4,5,6,8,10,15],colour=true,cmap="jet",calculate=true)
end

# To generate the area coverage.
percentUnder(a,thresh) = (sum(a.<thresh) / length(a)) * 100

# To generate the area coverage.
function percentUnder(sensors,f,d,thresh)
    emap=loadMap(sensors,f,d)
    errx=emap["errx"]
    erry=emap["erry"]
    percentUnder(errx,thresh),percentUnder(erry,thresh)
end

end
