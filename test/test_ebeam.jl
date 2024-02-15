using Revise, EPIC
using Test
using BenchmarkTools
using Plots
using DelimitedFiles
theme(:ggplot2)

begin
    ebeam=BunchedBeam(ELECTRON, 1.72e11, 10e9,  50000, [20e-9/0.8, 1.3e-9, 1.36e-4])
    opIPe=optics4DUC(0.45*0.8,0.0,0.056,0.0)
    vbase=3.42*8.5e6
    ϕs=10.0
    vact=vbase/cos(ϕs*π/180.0)
    mainRFe=AccelCavity(591e6, vact, 7560.0, π-ϕs*π/180.0)
    tunex, tuney=50.08, 44.14
    αc=3.42/tunex/tunex
    lmap=LongitudinalRFMap(αc, mainRFe)
    tunez=get_synchrotron_tune(ebeam, lmap)
    initilize_6DGaussiandist!(ebeam, opIPe, lmap)

    get_emittance!(ebeam)
    eqbs=1.0.*ebeam.beamsize


    ezcrab1=easyCrabCavity(394.0e6, 12.5e-3*0.0)
    #track!(ebeam, ezcrab1)
    lb=LorentzBoost(12.5e-3)
    #track!(ebeam, lb)
    #scatter(ebeam.dist.z, ebeam.dist.x; markersize=0.1, color=:oslo)
    opIPp=optics4DUC(0.8, 0.0, 0.072, 0.0)
    pstrong=StrongGaussianBeam(PROTON, 0.688e11, 275e9,  opIPp, [95e-6, 8.5e-6, 0.06], 9)
    initilize_zslice!(pstrong, :gaussian, :evennpar, 5.0)
    #lumi=track!(ebeam, pstrong) 
    ilb=InvLorentzBoost(12.5e-3)
    #track!(ebeam, ilb)
    ezcrab2=easyCrabCavity(394.0e6, 12.5e-3*0.0, π)
    crab_ratio=0.33
    overcrab=1.0
    pcrab1st = easyCrabCavity(197.0e6, overcrab*12.5e-3*(1+crab_ratio))
    pcrab2nd = easyCrabCavity(197.0e6*2.0, -overcrab*12.5e-3*crab_ratio)
    crab_crossing_setup!(pstrong, 12.5e-3, pcrab1st, pcrab2nd)

    #track!(ebeam, ezcrab2)
    #scatter(ebeam.dist.x, ebeam.dist.px; markersize=0.1, markercolor=:red)
    RLCwake = LongitudinalRLCWake(180e9, 5.5e3, 3.0)

    externalwake=readdlm("test/example_wake.txt", ' ', Float64, '\n')
    times=externalwake[1:5000,1]
    wakes=-externalwake[1:5000,2]
    
    arbWake=LongitudinalWake(times, wakes)



    oneturn=TransferMap4DChrom(opIPe, tunex, tuney, 2.0, 2.0)
    oneturndamp=OneTurnRadiation(4000.0, 4000.0, 2000.0)
    eqbs


    line1=Lattice([ezcrab1, lb])
    line2=Lattice([ilb, ezcrab2, oneturn, lmap])
    longiline=Lattice([lmap, RLCwake])
end


plot(pstrong.zslice_center, pstrong.xoffsets)

# externalwake=readdlm("test/example_wake.txt", ' ', Float64, '\n')
# times=externalwake[1:5000,1]
# wakes=-externalwake[1:5000,2]
# plot(times.*3e8, wakes, legend=false, xlabel="z [m]", ylabel="W [V/C]", title="External Wakefield")
# arbWake=LongitudinalWake(times, wakes)

##plot wake field
# wake1=RLCwake.wakefield.(ebeam.dist.z./3e8)
# wake2=arbWake.wakefield.(ebeam.dist.z./3e8)
# plot(ebeam.dist.z, wake1, legend=false, markersize=0.2, seriestype=:scatter, xlabel="z [m]", ylabel="W [V/C]", title="RLC Wakefield")
# plot!(ebeam.dist.z, wake2, legend=false, markersize=0.2, seriestype=:scatter, xlabel="z [m]", ylabel="W [V/C]", title="RLC Wakefield", xlim=(-0.01, 0.005))
# old_dp=ebeam.dist.dp.*1.0

# track!(ebeam, RLCwake)
# plot(ebeam.dist.z[1:100000], (ebeam.dist.dp.-old_dp)[1:100000], markersize=1, seriestype=:scatter, legend=false, xlabel="z [m]", ylabel="dp/p" )



turns=20000
coor_rec = Array{Float64, 2}(undef, 6, turns)
size_rec = Array{Float64, 2}(undef, 6, turns)
z2nd_rec = Array{Float64, 2}(undef, 6, turns)
lumis = Vector{Float64}(undef, turns)
formfactors = Vector{Float64}(undef, turns)
formfactor2s = Vector{Float64}(undef, turns)

function generate_noise(turns, centerfreq, bandwidth, nnoise)
    freqs = randn(nnoise) .* bandwidth .+ centerfreq
    phase = rand(Float64, nnoise) .* 2.0 .* π
    noiseraw = zeros(turns)
    for i in 1:nnoise
        noiseraw .+= sin.((2.0 * π * freqs[i]/78000) .* (1:turns) .+ phase[i])
    end
    ave = mean(noiseraw)
    noiseraw .-= ave
    size = sqrt(sum(noiseraw .* noiseraw) / turns)
    noiseraw .*= 1.0 / size
    return noiseraw
end

noise = generate_noise(turns, 120.0, 10.0, 1000)
fftnoise = fft(noise)
fftfreqs = fftfreq(turns, 1.0)
plot(fftfreqs, abs.(fftnoise)*1e-6/1.35, xlabel="Frequency [kHz]", ylabel="Amplitude", title="Noise Spectrum")


for i in 1:turns
    track!(ebeam, line1)
    lumis[i]=track!(ebeam, pstrong)
    track!(ebeam, line2)
    track!(ebeam, oneturndamp, eqbs)
    
    #track!(ebeam, longiline)
    get_emittance!(ebeam)
    coor_rec[:,i] .= ebeam.centroid
    size_rec[:,i] .= ebeam.beamsize
    z2nd_rec[:,i] .= ebeam.moment2nd[5,:]
    formfactors[i] = ebeam.num_macro*sum(ebeam.dist.x .^ 4)/sum(ebeam.dist.x .^ 2)^2
    formfactor2s[i] = ebeam.num_macro^2*sum(ebeam.dist.x .^ 6)/sum(ebeam.dist.x .^ 2)^3
    #ebeam.dist.x .+= 1e-6*noise[i]
end

size_rec


using FFTW

plot(@view size_rec[1,:])

plot(formfactor2s, label="Nominal paramter, no eCC", xlabel="Turns", ylabel="Form Factor, <x^6>/<x^2>^3")
plot(1e3.*z2nd_rec[1,:]./z2nd_rec[5,:])

plot(1e3.*z2nd_rec[2,:]./z2nd_rec[5,:])
plot(lumis)
fftx=fft(coor_rec[1,:])
ffty=fft(coor_rec[3,:])
fftz=fft(coor_rec[5,:])

plot!(fftfreqs, abs.(fftx), xlim=(0, 0.2), xlabel="Frequency [kHz]", ylabel="Amplitude", title="x")

get_emittance!(ebeam)
ebeam.emittance
ebeam.centroid
ebeam.beamsize

eqbs
ebeam.beamsize
colorcode= ebeam.dist.z .^ 2 ./ ebeam.beamsize[5]^2 .+ ebeam.dist.dp .^ 2 ./ ebeam.beamsize[6]^2
colorcode
plot(ebeam.dist.z, ebeam.dist.x, seriestype=:scatter, markersize=1, xlabel="z", ylabel="x", label="Nominal paramter, no eCC",)
savefig("test.png")

plot(histogram(ebeam.dist.x, nbins=100))
shapefactor=ebeam.num_macro*sum(ebeam.dist.x .^ 4)/sum(ebeam.dist.x .^ 2)^2
shapefactor2=ebeam.num_macro^2*sum(ebeam.dist.x .^ 6)/sum(ebeam.dist.x .^ 2)^3

dist08=copy(ebeam.dist)
formfac2_08 = copy(formfactor2s)
plot!(formfac2_08, label="0.8x e β_x, no eCC")
plot!(dist08.z, dist08.x, seriestype=:scatter, markersize=1, xlabel="z", ylabel="x", label="0.8x e β_x, no eCC",)


b_range = range(-10, 10, length=101)
plot()
histogram!(ebeam.dist.x ./ebeam.beamsize[1], bins=b_range, yaxis=(:log10, (1, 5000)), label="Nominal paramter, no eCC")
histogram!(dist08.x ./ std(dist08.x) , bins=b_range, yaxis=(:log10, (1, 5000)), label="0.8x e β_x, no eCC")

using Statistics
std(dist08.x)
ebeam.beamsize[1]
ebeam.dist.x