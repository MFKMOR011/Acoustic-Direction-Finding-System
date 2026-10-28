using LibSerialPort
using GLMakie
using FFTW
using Statistics
using WAV
using Dates
using LinearAlgebra


const PORT_NAME       = "COM3"
const FS              = 44118                 # Teensy runs at 44117.647 Hz
const BLOCK_SAMPLES   = 128
const N_MICS          = 3
const MAGIC           = UInt8[0x53, 0x54, 0x41, 0x52, 0x54]      # "START"
const BODY_BYTES      = 8 + N_MICS * BLOCK_SAMPLES * 2
const PKT_BYTES       = 5 + BODY_BYTES                           # 781

const BASELINE        = 0.017                 # m
const C               = 343.0
const FMAX            = 10088.0
const FMIN            = 2000.0
const FCUTOFF         = 2000.0
const FRAME_LEN       = 1024
const HOP             = FRAME_LEN ÷ 2
const FRAME_GATE      = 0.05                   # only use frames with >=5 % of the peak power

const SAVE_RECORDINGS = false
const OUTPUT_DIR      = "recordings/live"

# DOA Algorithm: MUSIC
# A 3x3 spatial covariance matrix is formed for every significant frequency bin across the frames,
# its noise subspace is found, and the pseudo-spectra of the bins are fused (weighted by bin energy).

const ANGLE_STEP    = 0.1                 # deg
const ANGLES = collect(0:ANGLE_STEP:360-ANGLE_STEP)
const THETAS = ANGLES .* (pi / 180)
const TAU21  = BASELINE .* sin.(pi/6 .- THETAS) ./ C     # delay of mic 2 relative to mic 1
const TAU31  = BASELINE .* cos.(THETAS) ./ C             # delay of mic 3 relative to mic 1
const WIN    = 0.5 .- 0.5 .* cos.(2 * pi .* (0:FRAME_LEN-1) ./ FRAME_LEN)
const FREQS  = (0:FRAME_LEN÷2) .* (FS / FRAME_LEN)
const INBAND = findall(f -> FMIN <= f <= FMAX, FREQS)

const BIN_FRACTION  = 0.3                 # keep bins with more than 30 % of the peak bin energy
const EIG_FRACTION  = 0.01                # eigenvalues above 1 % of the largest are counted as sources

function estimate_doa(raw::Matrix{Int16})
    x = Float64.(raw)
    x .-= mean(x, dims = 1)
    N = size(x, 1)

    # high-pass in frequency domain
    sc = round(Int, FCUTOFF / (FS / N))
    H  = zeros(N); H[sc+1:N-sc+1] .= 1.0
    sig = [real(ifft(H .* fft(x[:, m]))) for m in 1:N_MICS]

    Ns = length(sig[1])
    Nf = 1 + div(Ns - FRAME_LEN, HOP)
    nb = length(FREQS)

    # spectra of all frames, and the in-band power of mic 1 for the frame gate
    X  = Array{ComplexF64}(undef, N_MICS, nb, Nf)
    pw = zeros(Nf)
    for fi in 1:Nf
        idx = (fi-1)*HOP .+ (1:FRAME_LEN)
        for m in 1:N_MICS
            X[m, :, fi] = rfft(WIN .* sig[m][idx])
        end
        pw[fi] = sum(abs2, @view X[1, INBAND, fi])
    end
    keep = findall(p -> p >= FRAME_GATE * maximum(pw), pw)
    isempty(keep) && return NaN, 0
    Xk = X[:, :, keep]

    # energy of each bin over the kept frames and all microphones
    energy = [sum(abs2, @view Xk[:, b, :]) for b in 1:nb]
    thr = BIN_FRACTION * maximum(energy[INBAND])

    fused = zeros(length(ANGLES)); used = 0
    for b in INBAND
        energy[b] > thr || continue
        Sb = Xk[:, b, :]                              # 3 x number of frames
        R  = (Sb * Sb') / size(Sb, 2)                 # spatial covariance matrix
        E  = eigen(Hermitian(R))                      # eigenvalues in ascending order
        ev = max.(E.values, eps())
        ns = count(>(EIG_FRACTION * maximum(ev)), ev) # number of sources
        ns >= N_MICS && continue
        En = E.vectors[:, 1:N_MICS-ns]                # noise subspace
        P  = En * En'
        f  = FREQS[b]
        a2 = cis.(-2 * pi * f .* TAU21)
        a3 = cis.(-2 * pi * f .* TAU31)
        # a' * P * a for a = [1, a2, a3], for all angles at once
        den  = real(P[1,1] + P[2,2] + P[3,3]) .+ 2 .* real.(P[1,2] .* a2 .+ P[1,3] .* a3 .+ P[2,3] .* a3 .* conj.(a2))
        spec = 1 ./ max.(den, floatmin(Float64))
        fused .+= energy[b] .* (spec ./ maximum(spec))
        used += 1
    end
    used == 0 && return NaN, 0

    aoa = ANGLES[argmax(fused)]
    return round(aoa, digits = 2), length(keep)
end


function find_magic(buf)
    i = findfirst(==(0x53), buf)
    while i !== nothing
        i + 4 > length(buf) && return nothing
        @views(buf[i:i+4]) == MAGIC && return i
        i = findnext(==(0x53), buf, i + 1)
    end
    return nothing
end

function pull!(sp, buf)
    n = bytesavailable(sp)
    if n > 0
        try
            append!(buf, read(sp, n))
        catch e
            e isa LibSerialPort.Timeout || rethrow()
        end
        return true
    end
    return false
end


function receive_recording!(sp, buf, running)
    s = nothing
    while running[]
        s = find_magic(buf)
        s !== nothing && length(buf) >= s + 12 && break
        if s === nothing && length(buf) > 4
            deleteat!(buf, 1:length(buf)-4)
        end
        pull!(sp, buf) || sleep(0.005)
    end
    running[] || return nothing

    total = Int(reinterpret(UInt32, buf[s+9:s+12])[1])
    (total == 0 || total > 10_000) && (deleteat!(buf, 1:s+4); return nothing)
    need = s - 1 + total * PKT_BYTES

    t_last = time()
    while length(buf) < need && running[]
        if pull!(sp, buf)
            t_last = time()
        else
            time() - t_last > 1.0 && break       # stream stopped early
            sleep(0.002)
        end
    end

    raw = zeros(Int16, total * BLOCK_SAMPLES, N_MICS)
    got = 0
    for k in 0:total-1
        p = s + k * PKT_BYTES
        p + PKT_BYTES - 1 > length(buf) && break
        @views(buf[p:p+4]) == MAGIC || continue
        c = Int(reinterpret(UInt32, buf[p+5:p+8])[1])
        c < total || continue
        samples = reinterpret(Int16, buf[p+13:p+PKT_BYTES-1])
        raw[c*BLOCK_SAMPLES+1:(c+1)*BLOCK_SAMPLES, :] = reshape(samples, BLOCK_SAMPLES, N_MICS)
        got += 1
    end
    deleteat!(buf, 1:min(length(buf), need))
    return raw, total - got, total
end

# ───────────── GLMakie UI ─────────────
# equilateral triangle of side BASELINE, centred on its centroid (same geometry as the delay equations)
const mic_x = BASELINE .* [1/2, 0.0, -1/2]
const mic_y = BASELINE .* [-sqrt(3)/6, sqrt(3)/3, -sqrt(3)/6]

fig = Figure(size = (1200, 560))
ax  = Axis(fig[1, 1], aspect = DataAspect(), limits = (-1.3, 1.3, -1.3, 1.3),
           title = "MUSIC: mics and source (mic layout scaled ×1/baseline)")
scatter!(ax, mic_x ./ BASELINE, mic_y ./ BASELINE, markersize = 16, color = :black, label = "MICS")
text!(ax, mic_x ./ BASELINE, mic_y ./ BASELINE; text = ["1", "2", "3"], offset = (8, 8))
doa_pts  = Observable([Point2f(0, 0), Point2f(0, 0)])
doa_line = lines!(ax, doa_pts, color = :red, linewidth = 4, label = "DOA")
status   = Observable("Waiting for trigger…")
aoa_txt  = Observable("Angle of Arrival: —")
Label(fig[2, 1], aoa_txt, fontsize = 22, tellwidth = false)
Label(fig[3, 1], status,  fontsize = 14, tellwidth = false)

ax2 = Axis(fig[1:2, 2], title = "Last recording (mean removed)", xlabel = "time (s)", ylabel = "mic (offset)")
wf  = [Observable(Point2f[]) for _ in 1:N_MICS]
for m in 1:N_MICS
    lines!(ax2, wf[m], label = "Mic $m")
end
axislegend(ax2, position = :rt)

function update_ui!(raw, aoa, nframes, ndrop, total)
    if isnan(aoa)
        aoa_txt[] = "Angle of Arrival: no usable frames"
    else
        θ = deg2rad(aoa)
        doa_pts[] = [Point2f(0, 0), Point2f(1.2cos(θ), 1.2sin(θ))]
        aoa_txt[] = "Angle of Arrival: $aoa CCW  ($(nframes) frames)"
    end
    ds = 20
    t  = (0:ds:size(raw, 1)-1) ./ FS
    x  = Float64.(raw); x .-= mean(x, dims = 1)
    sc = max(maximum(abs, x), 1.0)
    for m in 1:N_MICS
        wf[m][] = Point2f.(t, x[1:ds:end, m] ./ sc .+ (m - 1) .* 2)
    end
    autolimits!(ax2)
    status[] = ndrop == 0 ? "Done. Waiting for next trigger…" :
                            "Done ($ndrop of $total blocks dropped). Waiting for next trigger…"
end


SAVE_RECORDINGS && mkpath(OUTPUT_DIR)
running = Ref(true)
display(fig)

listener = @async try
    LibSerialPort.open(PORT_NAME, 115200) do sp
        set_read_timeout(sp, 0.5)
        sp_flush(sp, SP_BUF_INPUT)
        buf = UInt8[]
        while running[]
            res = receive_recording!(sp, buf, running)
            res === nothing && continue
            raw, ndrop, total = res
            status[] = "Processing…"
            t0 = time()
            aoa, nframes = estimate_doa(raw)
            println("AOA = $aoa ($nframes frames, DOA took $(round(time()-t0, digits=2)) s)")
            update_ui!(raw, aoa, nframes, ndrop, total)
            if SAVE_RECORDINGS
                fn = joinpath(OUTPUT_DIR, Dates.format(now(), "yyyymmdd_HHMMSS") * "_$(aoa)deg.wav")
                wavwrite(raw, fn; Fs = FS)
            end
            sp_flush(sp, SP_BUF_INPUT); empty!(buf)
        end
    end
catch e
    @error "listener stopped" exception = (e, catch_backtrace())
    status[] = "Listener error: $(sprint(showerror, e))"
end

# keep the script alive until the window is closed
while isopen(fig.scene)
    sleep(0.1)
end
running[] = false