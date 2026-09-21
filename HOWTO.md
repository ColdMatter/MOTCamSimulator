# How to use MOTCamSimulator

## 1. Just look at the results

Open `CaF_MOT_Camera_Comparison.nb` in Mathematica / Wolfram (14.x). Every cell is already evaluated, so nothing needs to run. When the front end asks whether to enable dynamic content, click **Enable** — that activates the dashboard in Section 6, which carries its own definitions and works without evaluating anything else.

Sections, in order of usefulness:

| Section | What it answers |
|---|---|
| 2 Primary benchmark | Integrated SNR vs molecule number for every camera and binning; where each saturates |
| 3 SNR definitions | Peak-pixel vs top-hat vs matched-filter SNR; precision on number and on cloud width |
| 4 Visibility | How many molecules you need before the cloud shows in a single raw frame, vs binning |
| 5 Resolution | Object-space pixel size and resolution elements across the cloud, vs binning |
| 6 Dashboard | Interactive synthetic image, cross-section fit and metrics card |
| 7 Thresholds | Table + CSV of molecules for SNR = 3 / 10 per camera, readout mode and binning |

## 2. Use the dashboard (Section 6)

Controls, top to bottom: camera, readout mode (only cameras with more than one mode), binning (only the bins the camera supports), log10 molecules, scattering rate per molecule, exposure, background multiplier, noise seed.

- The image and grey points are one random realisation; the **noise seed** picks which one. Keep it fixed while you change other controls so you see only the effect of the change; step it to "take another shot".
- The metrics card is analytic and does not depend on the seed, except the two rows marked "from frame" and the fitted width.
- Sliders update on release (the frame generation takes ~0.1–0.3 s).

## 3. Change the physics

All inputs are top-level parameters in `CaFCameraNoise.wl`, Section 1. The ones most likely to need editing:

| Parameter | Section | Default | Meaning |
|---|---|---|---|
| `cloudFWHM` | 1.2 | 2.8 mm | measured MOT FWHM (object space); everything else about the cloud derives from it |
| `scatteringRate` | 1.2 | 1e6 /s | photons scattered per molecule per second — **assumed, replace with the measured value** |
| `baselineScatterRate`, `backgroundMultiplierDefault` | 1.2 | 39 /s/px, 4 | laser-scatter background per 6.5 µm pixel — **unverified** |
| `exposureTimeNominal` | 1.2 | 20 ms | exposure used for all benchmark plots and tables |
| `imagingLensRadius`, `imagingLensEFL` | 1.1 | 20 mm, 28.6 mm | second lens (Comar 29 AF 40); sets the stop and the magnification |
| `viewportClearRadius` | 1.1 | Infinity | set to the viewport semi-diameter if it clips before the lenses |
| `workingPointMolecules` | 3 | 1e4 | the molecule number used in the summary tables |

After editing, regenerate (step 5). The verification suite will tell you if a change violates an assumption (e.g. the synthetic frame is too small to hold the reference annulus).

## 4. Add a camera

Append an entry to `cameraDatabase` in Section 1.3. All keys are required:

```mathematica
"Vendor Model (type)" -> <|
   "ShortName" -> "Model", "Type" -> "sCMOS",
   "PixelsX" -> 2048, "PixelsY" -> 2048, "PixelPitch" -> 6.5,     (* um *)
   "QE" -> 0.80,                                                     (* at 606 nm, not the peak *)
   "ReadNoiseRMS" -> 1.0,                                            (* e- rms, primary mode *)
   "ReadNoiseModes" -> <|"Mode name" -> 1.0|>,                       (* every quoted mode *)
   "DarkCurrent" -> 0.1,                                             (* e-/pixel/s at the operating temperature *)
   "FullWell" -> 30000,                                              (* e- *)
   "BinningType" -> "Software",                                      (* label only *)
   "HardwareBinLimit" -> 1,      (* 1 = digital binning only, Infinity = on-chip for all bins (CCD/EMCCD), 2 = on-chip 2x2 then digital *)
   "SupportedBins" -> {1, 2, 4, 8},
   "ENF" -> 1.0,                                                     (* Sqrt[2] for EMCCD with gain *)
   "Notes" -> "Where the numbers came from."|>
```

Read the QE off the vendor's curve at 606 nm — the headline "peak QE" is usually at 500–550 nm and can be 10 points higher. Regenerate; the new camera appears in every plot, table and the dashboard automatically.

## 5. Regenerate the notebook

```bash
./regenerate.sh
```

or, explicitly, `wolframscript -file GenerateCaFComparisonNotebook.wls`. On this Mac `wolframscript` is not on the PATH; the script looks in `/Applications/Wolfram.app` and `/Applications/Mathematica.app` for you. Takes ~30 s and produces:

- `CaF_MOT_Camera_Comparison.nb` — rebuilt from scratch (never edit it by hand)
- `CaF_MOT_Camera_Thresholds.csv`
- `verification/*.png` — every key figure, for checking without opening Mathematica
- a PASS/FAIL list for the 41 verification tests; the exit code is 1 if any fail

## 6. One-off calculations without the notebook

Load the model into any kernel and call the functions directly:

```mathematica
Get["CaFCameraNoise.wl"];
cam = cameraDatabase["Thorlabs Zelux CS165MU (Compact CMOS)"];

integratedSNRMolecules[cam, 4, 10^4]              (* top-hat SNR, 4x4 bin, 1e4 molecules, nominal conditions *)
peakPixelSNR[cam, 4, 10^4]                        (* single-shot visibility *)
visibilityMolecules[cam, 4, 5]                    (* molecules needed for peak-pixel SNR = 5 *)
minimumMolecules[cam, 1, 3]                       (* molecules for integrated SNR = 3 *)
snrEstimators[cam, 1, 10^4]                       (* peak / top-hat / matched-filter SNR, dN/N, dSigma/Sigma *)
noiseFloorBreakdown[cam, 1, 0.02, 4., 4.]          (* read vs scatter vs dark variance, dominant term *)
```

Every function has a long form with explicit `tExp` (s), background multiplier, read noise (e⁻) and scattering rate (/s) — see Section 1.4/1.5 of the `.wl`.
