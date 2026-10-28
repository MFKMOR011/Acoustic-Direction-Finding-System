# Acoustic-Direction-Finding-System

Final-year project (EEE4022S), Mechatronics Engineering, University of Cape Town.
Author: Morena Mofokeng (MFKMOR011). Supervisor: Assoc. Prof. A. J. Wilkinson.

A low-cost system that estimates the direction (0 to 360 degrees) of a single sound source in the horizontal plane. Three SPH0645 I2S microphones on a 17 mm equilateral triangle are sampled by a Teensy 4.0 at about 44.1 kHz and streamed over USB to a laptop. Julia programs estimate the angle and display it in real time. Three algorithms are implemented and compared: Trig + FFT, MUSIC, and a look-up table with a cost function (LUT + cost, the chosen algorithm).

## Repository layout

| Folder | Contents |
|---|---|
| `firmware/` | Teensy 4.0 sketch: records the three channels and sends them over USB serial |
| `realtime/` | Julia programs that receive the data, estimate the angle and display it: `realtime_trig.jl`, `realtime_music.jl`, `realtime_lut.jl` |
| `simulation/` | Signal simulation scripts for Trig + FFT, MUSIC and LUT + cost |
| `montecarlo/` | `MonteCarlo.ipynb` and the result tables it writes to `results/` |
| `recordings/` | Recorded data; the file name gives the distance, true angle and signal, for example `180cm_90_d_man_sin_4000hz.wav` |
| `cad/` | CAD files of the enclosure and the changes made to the turntable |
| `docs/` | SPH0645 microphone datasheet and the Teensy 4.0 processor manual |
| `report/` | Report pdf|

## Hardware

- Teensy 4.0 (600 MHz) with three SPH0645 I2S microphone boards
- USB cable that supports data transfer
- Breadboard, green LED and resistor (recording indicator)
- 3D-printed enclosure; turntable (Big Simi, "Manual Turning Table", MakerWorld) with the changes in `cad/`

Microphone geometry: equilateral triangle, side 17 mm. Microphone 1 is at the origin, microphone 3 at (-b, 0) and microphone 2 at (-b/2, 0.866b), with b = 17 mm. Angles are measured counter-clockwise from the +x axis.

## Setup

### 1. Firmware

1. Install the Arduino IDE with Teensyduino.
2. Open the sketch in `firmware/`, select **Teensy 4.0** and upload.
3. Connect the board by USB. The green LED is the status indicator.

The firmware sends packets that begin with the marker `START` followed by a counter and 128 samples from each of the three microphones (781 bytes per packet). The host opens the port at 115200 baud.

### 2. Julia

1. Install Julia.
2. From the repository folder, install the packages:

```julia
using Pkg
Pkg.add(["LibSerialPort", "GLMakie", "FFTW", "WAV", "Statistics", "Dates"])
```

MUSIC also uses `LinearAlgebra`, which ships with Julia.

## Running the real-time program

1. Plug in the Teensy and find its serial port (Windows: Device Manager, for example `COM3`).
2. Open the program you want to run in `realtime/` and set the port at the top:

```julia
const PORT_NAME = "COM3"
```

3. Run it:

```julia
include("realtime/realtime_lut.jl")     # LUT + cost function (chosen algorithm)
include("realtime/realtime_trig.jl")    # Trig + FFT
include("realtime/realtime_music.jl")   # MUSIC
```

A polar display opens and shows the estimated direction of the source. Each program works on 1024-sample frames, uses only frames with at least 5 % of the peak power, and uses the 2 to 10 kHz band.

Optional settings at the top of each program:

| Setting | Meaning |
|---|---|
| `PORT_NAME` | Serial port of the Teensy |
| `SAVE_RECORDINGS` | Set to `true` to save the received audio as WAV files in `recordings/live` |
| `ANGLE_STEP` | Angle resolution of the look-up table (LUT program only) |

## Simulations and Monte Carlo analysis

1. Run the scripts in `simulation/` for the known-angle simulation of the three algorithms.
2. Open `montecarlo/MonteCarlo.ipynb` in Jupyter (IJulia) and run all cells. It reproduces the experiments on noise, microphone position error, gain and delay mismatch, and frame and observation length. The tables are written to `montecarlo/results/`.

## Recordings

The files in `recordings/` are the real-array measurements used in Chapter 5. The file name gives the distance from the source, the true angle and the signal (sine, speech or clap). Source positions were set on the turntable.

## Report

The report source is in `report/`. 

## Licence and credits

- Turntable model: Big Simi, "Manual Turning Table (Fully 3D printable)", MakerWorld, 2025. No licence is stated on the model page.
