# MOTCamSimulator

Analytical noise model, benchmark plots and an interactive dashboard for imaging CaF MOT fluorescence (606 nm) with seven candidate cameras: Hamamatsu ORCA-Quest 2 and ORCA-R2, Thorlabs Zelux CS165MU, Teledyne Kinetix, Andor iXon Ultra 888, Andor CB2 High Res and CB2 High Speed 7.1F.

## Quick start

Open `CaF_MOT_Camera_Comparison.nb` in Mathematica / Wolfram 14.x. Everything is pre-evaluated; click **Enable** when asked about dynamic content to activate the dashboard.

| Section | What it answers |
|---|---|
| 1 Noise engine | Optics, cloud, camera database, SNR functions |
| 2 Primary benchmark | Integrated SNR vs molecule number for every camera and binning; where each saturates |
| 3 SNR definitions | Peak-pixel vs top-hat vs matched-filter SNR; precision on number and on cloud width |
| 4 Visibility | Molecules needed before the cloud shows in a single raw frame, vs binning |
| 5 Resolution | Object-space pixel size and resolution elements across the cloud, vs binning |
| 6 Dashboard | Interactive synthetic image, cross-section fit and metrics card |
| 7 Thresholds | Table + CSV of molecules for SNR = 3 / 10 per camera, readout mode and binning |
| 8 Tests | The verification suite |

**Dashboard controls**: camera, readout mode, binning (only the bins the camera supports), log10 molecules, scattering rate, exposure, background multiplier, noise seed. The image is one random realisation chosen by the seed — keep it fixed while changing other controls, step it to "take another shot". The metrics card is analytic except the rows marked "from frame" and the fitted width. Sliders update on release.

## Files

| File | Purpose |
|---|---|
| `CaFCameraNoise.wl` | The model: physics, camera database, plots, dashboard and tests, in Wolfram package-cell format. **Edit this.** |
| `GenerateCaFComparisonNotebook.wls` | Builds and evaluates the notebook from the `.wl`, runs the tests, exports the CSV and figures. |
| `regenerate.sh` | Finds `wolframscript` and runs the generator. |
| `CaF_MOT_Camera_Comparison.nb` | Generated notebook — never edit by hand, it is overwritten. |
| `CaF_MOT_Camera_Thresholds.csv` | Molecules for SNR = 3 / 10 and saturation molecule number per camera, mode and binning. |
| `verification/*.png` | Every key figure, for checking without opening Mathematica. |

## Change the physics

All inputs are top-level parameters in `CaFCameraNoise.wl`, Section 1:

| Parameter | Section | Default | Meaning |
|---|---|---|---|
| `cloudFWHM` | 1.2 | 2.8 mm | measured MOT FWHM (object space); σ, aperture and fields of view derive from it |
| `scatteringRate` | 1.2 | 1e6 /s | photons scattered per molecule per second — **assumed, replace with the measured value** |
| `baselineScatterRate`, `backgroundMultiplierDefault` | 1.2 | 39 /s/px, 4 | laser-scatter background per 6.5 µm pixel — **unverified** |
| `exposureTimeNominal` | 1.2 | 20 ms | exposure for all benchmark plots and tables |
| `imagingLensRadius`, `imagingLensEFL` | 1.1 | 20 mm, 28.6 mm | second lens (Comar 29 AF 40); sets the limiting stop and the magnification |
| `viewportClearRadius` | 1.1 | Infinity | set to the viewport semi-diameter if it clips before the lenses |
| `workingPointMolecules` | 3 | 1e4 | molecule number used in the summary tables |

Then regenerate:

```bash
./regenerate.sh
```

(~30 s; equivalent to `wolframscript -file GenerateCaFComparisonNotebook.wls`, the script also looks inside `/Applications/Wolfram.app` and `/Applications/Mathematica.app` if `wolframscript` is not on your PATH). It rebuilds the notebook, CSV and figures and prints PASS/FAIL for the 41 verification tests; exit code 1 if any fail.

## Add a camera

Append an entry to `cameraDatabase` (Section 1.3). All keys are required:

```mathematica
"Vendor Model (type)" -> <|
   "ShortName" -> "Model", "Type" -> "sCMOS",
   "PixelsX" -> 2048, "PixelsY" -> 2048, "PixelPitch" -> 6.5,     (* um *)
   "QE" -> 0.80,                                                     (* at 606 nm, not the peak *)
   "ReadNoiseRMS" -> 1.0,                                            (* e- rms, primary mode *)
   "ReadNoiseModes" -> <|"Mode name" -> 1.0|>,                       (* every quoted mode *)
   "DarkCurrent" -> 0.1,                                             (* e-/pixel/s at operating temperature *)
   "FullWell" -> 30000,                                              (* e-, only used for saturation markers *)
   "BinningType" -> "Software",                                      (* label only *)
   "HardwareBinLimit" -> 1,      (* 1 = digital binning only, Infinity = on-chip at every bin (CCD/EMCCD), 2 = on-chip 2x2 then digital *)
   "SupportedBins" -> {1, 2, 4, 8},
   "ENF" -> 1.0,                                                     (* Sqrt[2] for EMCCD with gain *)
   "Notes" -> "Where the numbers came from."|>
```

Read the QE off the vendor's curve at 606 nm — the headline peak is usually at 500–550 nm and can be 10 points higher. Regenerate; the camera appears in every plot, table and the dashboard.

## One-off calculations

```mathematica
Get["CaFCameraNoise.wl"];
cam = cameraDatabase["Thorlabs Zelux CS165MU (Compact CMOS)"];

integratedSNRMolecules[cam, 4, 10^4]        (* top-hat SNR, 4x4 bin, 1e4 molecules, nominal conditions *)
peakPixelSNR[cam, 4, 10^4]                  (* single-shot visibility *)
visibilityMolecules[cam, 4, 5]              (* molecules needed for peak-pixel SNR = 5 *)
minimumMolecules[cam, 1, 3]                 (* molecules for integrated SNR = 3 *)
snrEstimators[cam, 1, 10^4]                 (* peak / top-hat / matched-filter SNR, dN/N, dSigma/Sigma *)
noiseFloorBreakdown[cam, 1, 0.02, 4., 4.]    (* read vs scatter vs dark variance, dominant term *)
```

Every function has a long form with explicit exposure (s), background multiplier, read noise (e⁻) and scattering rate (/s); see Sections 1.4–1.5 of the `.wl`.

## How the model works

- **Signal**: `nSig = nMol × scatteringRate × tExp × etaGeom × T × fAperture` photons arrive inside the 2σ integration aperture (86.5 % of the cloud's fluorescence); signal electrons are `nSig × QE`. All axes are in molecules (10 to 10⁶).
- **Optics**: LA1401-A collector (2", EFL 59.8 mm) at its BFL of 49.1 mm, Comar 29 AF 40 imaging lens (40 mm, EFL 28.6 mm). The Comar is the limiting stop, giving a geometric efficiency of 3.69 % (5.59 % would be the unvignetted 2" lens). Magnification M = f2/f1 = 0.478: a pixel of pitch p samples p/M at the MOT. Cloud geometry is kept in object space; background and dark electrons use the physical sensor pixel area and do not scale with M.
- **Noise per superpixel**: shot noise on signal, background and dark electrons, plus read noise `sigmaRead × b / Min[b, HardwareBinLimit]`; a factor 2 on all signal-independent terms accounts for subtracting an identical off-target reference. Software binning sums `b²` independent reads, so it leaves the aperture-integrated read noise unchanged — only on-chip binning lowers it. It does raise the peak-pixel SNR in proportion to b, which is what matters for seeing the cloud in a single shot.
- **Three SNRs**: peak-pixel (visibility), top-hat over the 2σ aperture (the conventional integrated SNR), and the matched filter (inverse-variance, profile-weighted), whose variance is the Cramér–Rao bound on the molecule number. δN/N and δσ/σ follow from the Fisher information.
- **Regimes**: the log-log slope is 1 where a signal-independent term (read noise, laser scatter or dark current) dominates and ½ once signal shot noise dominates; the regime table in Section 2 says which term sets the floor for each camera. Filled markers show where the peak superpixel reaches full well.
