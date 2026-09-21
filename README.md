# CaF MOT Camera Comparison

Analytical noise model, benchmark plots and an interactive dashboard for imaging CaF MOT fluorescence (606 nm) with five candidate cameras.

## Files

| File | Purpose |
|---|---|
| `CaFCameraNoise.wl` | All physics, camera database, plots, dashboard and verification tests, in Wolfram package-cell format. `Get`-able standalone, and the front end can open it directly as a notebook. |
| `GenerateCaFComparisonNotebook.wls` | Generator: parses the `.wl`, evaluates every Input cell in the kernel, runs the `VerificationTest` suite, writes the notebook through the headless front end, rasterises verification figures and exports the threshold CSV. Exit code 1 on any failed test. |
| `CaF_MOT_Camera_Comparison.nb` | Generated, fully evaluated notebook (Input + Output cells; the `Manipulate` dashboard carries its own definitions via `SaveDefinitions`). |
| `CaF_MOT_Camera_Thresholds.csv` | Minimum molecule number (and equivalent photons) for SNR = 3 and 10, plus the saturation molecule number, per camera, readout mode and binning. |
| `verification/*.png` | Rasterised figures produced by the last generator run. |

## Regenerating

Quick start for users is in [HOWTO.md](HOWTO.md). The one-liner is `./regenerate.sh`, which locates `wolframscript` and runs the generator.

`wolframscript` is not on `$PATH` on this machine; use the Wolfram.app binary directly
(`/Applications/Mathematica.app` reports a licence problem, `/Applications/Wolfram.app` 14.3 is licensed):

```bash
/Applications/Wolfram.app/Contents/MacOS/wolframscript -file "GenerateCaFComparisonNotebook.wls"
```

Or add it to your path once (`~/.zshrc`):

```bash
export PATH="/Applications/Wolfram.app/Contents/MacOS:$PATH"
```

To edit the model, change `CaFCameraNoise.wl` and re-run the generator; never edit the `.nb` by hand (it is overwritten).

## Notes on the model

- **Cloud size**: the primary input is the measured **FWHM = 2.8 mm** (σ = 1.189 mm, 1/e² diameter 4.76 mm), which supersedes the 1.5 mm 1/e² diameter of the original brief. The 2σ integration aperture is therefore 2.378 mm in radius and covers 10× more pixels than the original spec, which raises the integrated read-noise floor by √10 and is the single biggest lever on the low-molecule SNR. The synthetic-frame and estimator fields of view scale with the cloud automatically.
- **Optical relay**: Thorlabs LA1401-A collector (2", EFL 59.8 mm) at BFL 49.1 mm + Comar 29 AF 40 achromat (40 mm dia, EFL 28.6 mm). The Comar is the limiting stop, so the geometric efficiency is 3.69 %, not the 5.59 % of the unvignetted 2" lens (`viewportClearRadius` models any tighter upstream stop). The relay is **not** 1:1: M = f2/f1 = 0.478, so a pixel of pitch p samples p/M at the MOT and the cloud image is demagnified 2.1×. Cloud geometry is kept in object space; background and dark electrons use the physical sensor pixel area (they do not scale with M).
- **Three SNR definitions** are computed (Section 1.5, compared in Section 3): peak-pixel, flat top-hat over the 2σ aperture, and the matched filter (inverse-variance, profile-weighted), whose variance is the Cramér–Rao bound on molecule number. Section 3 also gives δN/N and δσ/σ, the figures of merit for number and for time-of-flight thermometry.
- **Visibility (Section 4)** is kept separate from precision: `peakPixelSNR` and `visibilityMolecules` give the SNR of the brightest superpixel and the molecule number at which the cloud becomes unmistakable in a raw single shot (peak-pixel SNR = 5). This is the one figure of merit that software binning improves — in proportion to b, for every camera.
- All benchmark axes and tables are in **number of molecules in the MOT** (10 to 10^6). Internally `nSig = nMol * scatteringRate * tExp * etaGeom * T * fAperture` is the number of signal photons arriving inside the 0.75 mm integration aperture (2 sigma of the cloud, 86.5 % of the fluorescence), about 920 photons per molecule at 20 ms; signal electrons are `nSig * QE`. The scattering rate (default 1e6 photons/s per molecule) is a dashboard control.
- Filled markers on the benchmark curves show where the peak superpixel reaches the full-well capacity at 20 ms  — with the 2.8 mm cloud there is ample headroom: saturation at 20 ms arrives at 1.1–3.9 × 10⁶ molecules unbinned.
- **Andor CB2 High Res and High Speed 7.1F** were added from `cb2-specifications.pdf` (QE at 606 nm read off the brochure curves: 62 % / 68 %, validated against the quoted 74 % peaks). They bin 2×2 in the charge domain on chip and digitally beyond, so the read-noise rule was generalised to `sigmaReadEff = sigmaRead * b / Min[b, HardwareBinLimit]` with `HardwareBinLimit` = 1 (digital only), ∞ (CCD/EMCCD) or 2 (CB2).
- Software binning (qCMOS, sCMOS, CMOS) sums `b^2` independent reads, so the aperture-integrated read-noise variance is independent of `b`; the 1x1/2x2/4x4 SNR curves of those cameras coincide. Only hardware binning (ORCA-R2, iXon) lowers the floor.
- The log-log slope is 1 wherever a signal-independent term (read noise **or** laser scatter **or** dark current) dominates, and 1/2 once the signal shot noise dominates; the regime table in Section 2 reports which floor term dominates for each camera.
- Kinetix and iXon full-well values and the default 1e6 photons/s per molecule scattering rate are not part of the brief; the scattering rate scales every molecule axis linearly, the full-well values only affect the saturation markers.
