(* ::Package:: *)

(* ::Title:: *)
(*CaF MOT Fluorescence Imaging: Camera Comparison*)


(* ::Text:: *)
(*Analytical noise model and interactive dashboard for imaging the fluorescence of a CaF magneto-optical trap (A2Pi1/2 -> X2Sigma+ (0,0) band at 606 nm) with five candidate cameras. The notebook is generated and evaluated automatically by GenerateCaFComparisonNotebook.wls from the source file CaFCameraNoise.wl. Every user-defined symbol is lowerCamelCase; underscores appear only in Blank[] patterns.*)


(* ::Text:: *)
(*Units convention: lengths in mm unless stated, pixel pitches in micrometres, times in seconds, charges in electrons, photon fluxes in photons per second per square micrometre at the sensor.*)


(* ::Section:: *)
(*1. Analytical Noise Engine*)


(* ::Subsection:: *)
(*1.1 Optical relay geometry (Thorlabs LA1401-A collector + Comar 29 AF 40 imaging lens)*)


(* ::Text:: *)
(*The MOT sits at the back focal length of the plano-convex singlet (flat side facing the MOT), so the fluorescence leaves the collector collimated with a beam diameter equal to the illuminated part of that lens (up to 50.8 mm). The second lens, a Comar 29 AF 40 achromat (40 mm diameter, EFL 28.6 mm), focuses the beam onto the sensor.*)


(* ::Text:: *)
(*Two consequences of the real relay, both of which matter for the SNR:*)


(* ::Text:: *)
(*(i) Vignetting. The collimated beam is clipped by whichever element has the smallest clear aperture; here the 40 mm Comar, not the 2 inch collector, is the stop. The accepted solid angle is that of a disc of the limiting radius at the working distance, so the geometric efficiency drops from 5.59% (unvignetted 2 inch) to 3.69%. Any smaller restriction upstream (a viewport, a re-entrant window) would reduce it further: set viewportClearRadius to model it.*)


(* ::Text:: *)
(*(ii) Magnification. With two lenses in an infinity relay the magnification is M = f2/f1 = 28.6/59.8 = 0.478, so the MOT image is demagnified. This has no effect on the total number of collected photons, but it concentrates them onto 1/M^2 = 4.4 times fewer pixels, which lowers the integrated read-noise and dark-current floor by the same factor, raises the peak electrons per pixel (earlier saturation) and coarsens the object-space sampling to pixelPitch/M.*)


(* ::Input:: *)
collectorRadius = 25.4;           (* mm, LA1401-A clear semi-diameter (2 inch diameter) *)
collectorEFL = 59.8;              (* mm *)
collectorBFL = 49.1;              (* mm *)
collectorCenterThickness = 16.3;  (* mm *)
workingDistance = collectorBFL;   (* mm, MOT to flat collector face *)
imagingLensRadius = 20.0;         (* mm, Comar 29 AF 40 achromat: 40 mm diameter *)
imagingLensEFL = 28.6;            (* mm *)
viewportClearRadius = Infinity;   (* mm, set to the vacuum viewport semi-diameter if it clips before the lenses *)
transmissionOptics = 0.95;        (* AR coated lens surfaces and viewports *)
wavelengthNm = 606;               (* nm, CaF A-X (0,0) *)
limitingApertureRadius = Min[collectorRadius, imagingLensRadius, viewportClearRadius];
limitingElement = First[Keys[TakeSmallest[<|"Collector (2 in)" -> collectorRadius, "Comar 29 AF 40" -> imagingLensRadius, "Viewport" -> viewportClearRadius|>, 1]]];
magnification = imagingLensEFL/collectorEFL;   (* infinity relay: M = f2/f1 = 0.478 *)
solidAngleCollection[d0_, r_] := 2 Pi (1 - d0/Sqrt[d0^2 + r^2]);
solidAngleMOT = solidAngleCollection[workingDistance, limitingApertureRadius];
geometricEfficiency = solidAngleMOT/(4 Pi);
geometricEfficiencyUnvignetted = solidAngleCollection[workingDistance, collectorRadius]/(4 Pi);   (* 2 inch collector with no downstream clipping, for reference *)
vignettingFactor = geometricEfficiency/geometricEfficiencyUnvignetted;
numericalApertureCollection = Sin[ArcTan[limitingApertureRadius/workingDistance]];
imagingLensFNumber = imagingLensEFL/(2 imagingLensRadius);
Grid[{{"Limiting aperture", Row[{limitingElement, ", radius ", limitingApertureRadius, " mm"}]},
      {"Solid angle (sr)", solidAngleMOT},
      {"Geometric collection efficiency (vignetted)", geometricEfficiency},
      {"Geometric collection efficiency (2 in, unvignetted)", geometricEfficiencyUnvignetted},
      {"Vignetting factor", vignettingFactor},
      {"Collection NA", numericalApertureCollection},
      {"Magnification M = f2/f1", magnification},
      {"Imaging lens f-number", imagingLensFNumber},
      {"Total optical throughput (geometry x transmission)", geometricEfficiency transmissionOptics}},
     Alignment -> Left, Dividers -> All]


(* ::Subsection:: *)
(*1.2 Molecular cloud, integration aperture and laser-scatter background*)


(* ::Text:: *)
(*The cloud is a 2D Gaussian of measured FWHM 2.8 mm in object space at the MOT (sigma = FWHM/(2 Sqrt[2 Log[2]]) = 1.189 mm, 1/e^2 diameter 4.76 mm). The FWHM is the primary input; sigma and the 1/e^2 radius derive from it. The integration aperture is a disc of radius 2 sigma = 2.378 mm at the MOT, which contains a fraction 1 - Exp[-2] = 86.5% of the cloud fluorescence. Signal photon numbers quoted in this notebook (nSig) are the photons arriving on the sensor inside that aperture. All cloud geometry is kept in object space and pixels are projected back through the magnification: a superpixel of physical pitch b x pixelPitch on the chip covers b x pixelPitch / M at the MOT. Background and dark electrons, in contrast, are properties of the physical sensor pixel and do not scale with M.*)


(* ::Text:: *)
(*The laser-scatter background is specified as 39 photons/s per 6.5 um x 6.5 um pixel at baseline; the active condition is 4x baseline, i.e. 3.69 photons/(s um^2). The multiplier is a free parameter (0.5 to 10).*)


(* ::Text:: *)
(*Molecule number to signal photons: each molecule scatters scatteringRate photons/s, of which the fraction geometricEfficiency x transmissionOptics reaches the sensor and apertureSignalFraction of those land inside the integration aperture, so nSig = nMol x scatteringRate x tExp x geometricEfficiency x transmissionOptics x apertureSignalFraction. The scattering rate (default 1e6 photons/s per molecule, a typical CaF MOT value) is an assumed parameter and is exposed as a dashboard control; all benchmark axes are in molecules.*)


(* ::Input:: *)
cloudFWHM = 2.8;                                     (* mm, measured FWHM of the MOT cloud in object space *)
sigmaCloud = cloudFWHM/(2 Sqrt[2 Log[2]]);           (* mm, ~1.189 mm *)
cloudWaist = 2 sigmaCloud;                           (* mm, 1/e^2 radius, ~2.378 mm *)
apertureRadius = 2 sigmaCloud;                       (* mm *)
apertureAreaMicron = Pi (1000 apertureRadius)^2;     (* um^2 *)
apertureSignalFraction = 1 - Exp[-apertureRadius^2/(2 sigmaCloud^2)];
exposureTimeNominal = 0.020;                         (* s *)
exposureTimeRange = {0.001, 0.200};                  (* s *)
baselineScatterRate = 39.;                           (* photons/s per baseline pixel *)
baselinePixelPitch = 6.5;                            (* um *)
baselinePixelArea = baselinePixelPitch^2;            (* um^2 *)
backgroundMultiplierDefault = 4.0;
backgroundMultiplierRange = {0.5, 10.0};
backgroundFluxDensity[multiplier_] := multiplier baselineScatterRate/baselinePixelArea;  (* photons/(s um^2) *)
scatteringRate = 1.0*^6;   (* photons/s scattered per molecule in the MOT: assumed typical CaF value, replace with the measured rate *)
scatteringRateRange = {1.0*^5, 5.0*^6};
moleculeNumberRange = {10, 10^6};
photonsPerMolecule[tExp_, rate_] := rate tExp geometricEfficiency transmissionOptics apertureSignalFraction;   (* signal photons arriving in the aperture per molecule *)
photonsPerMolecule[tExp_] := photonsPerMolecule[tExp, scatteringRate];
signalPhotonsFromMolecules[nMol_, tExp_, rate_] := nMol photonsPerMolecule[tExp, rate];
signalPhotonsFromMolecules[nMol_, tExp_] := nMol photonsPerMolecule[tExp, scatteringRate];
moleculesFromSignalPhotons[nSig_, tExp_, rate_] := nSig/photonsPerMolecule[tExp, rate];
moleculesFromSignalPhotons[nSig_, tExp_] := nSig/photonsPerMolecule[tExp, scatteringRate];
Grid[{{"Cloud FWHM (mm)", cloudFWHM},
      {"Cloud sigma (mm)", sigmaCloud},
      {"Cloud 1/e^2 diameter (mm)", 2 cloudWaist},
      {"Aperture radius (mm)", apertureRadius},
      {"Aperture area (um^2)", apertureAreaMicron},
      {"Fraction of cloud signal inside aperture", apertureSignalFraction},
      {"Background flux density at 4x baseline (photons/(s um^2))", backgroundFluxDensity[backgroundMultiplierDefault]},
      {"Background flux density at 4x baseline (photons/(s m^2))", 10.^12 backgroundFluxDensity[backgroundMultiplierDefault]},
      {"Signal photons in aperture per molecule at 20 ms (scattering rate 1e6 /s)", photonsPerMolecule[exposureTimeNominal]}},
     Alignment -> Left, Dividers -> All]


(* ::Subsection:: *)
(*1.3 Camera specification database*)


(* ::Text:: *)
(*ReadNoiseRMS is the primary (lowest-noise) readout mode used for the benchmark curves; ReadNoiseModes lists all quoted modes. QE values are at 606 nm, read off each vendor's published curve where no number is quoted at that wavelength (the CB2 curves give 62% High Res and 68% High Speed). FullWell is only used for the saturation markers. HardwareBinLimit encodes how binning is implemented: 1 = digital only, Infinity = on-chip for every bin (CCD serial register, EMCCD), 2 = on-chip 2x2 charge-domain binning with digital summing beyond (Andor CB2).*)


(* ::Input:: *)
cameraDatabase = <|
  "Hamamatsu ORCA-Quest 2 (qCMOS C15550-22UP)" -> <|
     "ShortName" -> "ORCA-Quest 2", "Type" -> "qCMOS",
     "PixelsX" -> 4096, "PixelsY" -> 2304, "PixelPitch" -> 4.6,
     "QE" -> 0.72, "ReadNoiseRMS" -> 0.30,
     "ReadNoiseModes" -> <|"Ultra Quiet" -> 0.30, "Standard" -> 0.43|>,
     "DarkCurrent" -> 0.006, "FullWell" -> 7000,
     "BinningType" -> "Software", "HardwareBinLimit" -> 1, "SupportedBins" -> {1, 2, 4, 8}, "ENF" -> 1.0,
     "Notes" -> "Photon-number resolving qCMOS; digital binning only."|>,
  "Hamamatsu ORCA-R2 (Cooled CCD C10600-10B)" -> <|
     "ShortName" -> "ORCA-R2", "Type" -> "Interline CCD",
     "PixelsX" -> 1344, "PixelsY" -> 1024, "PixelPitch" -> 6.45,
     "QE" -> 0.68, "ReadNoiseRMS" -> 6.0,
     "ReadNoiseModes" -> <|"Normal Scan" -> 6.0, "Fast Scan" -> 10.0|>,
     "DarkCurrent" -> 0.0005, "FullWell" -> 18000,
     "BinningType" -> "Hardware", "HardwareBinLimit" -> Infinity, "SupportedBins" -> {1, 2, 4, 8}, "ENF" -> 1.0,
     "Notes" -> "On-chip binning in the serial register: one read per superpixel."|>,
  "Thorlabs Zelux CS165MU (Compact CMOS)" -> <|
     "ShortName" -> "Zelux CS165MU", "Type" -> "Standard CMOS",
     "PixelsX" -> 1440, "PixelsY" -> 1080, "PixelPitch" -> 3.45,
     "QE" -> 0.68, "ReadNoiseRMS" -> 4.0,
     "ReadNoiseModes" -> <|"Standard" -> 4.0|>,
     "DarkCurrent" -> 3.0, "FullWell" -> 11000,
     "BinningType" -> "Software", "HardwareBinLimit" -> 1, "SupportedBins" -> {1, 2, 4, 8, 16}, "ENF" -> 1.0,
     "Notes" -> "Uncooled; dark current quoted at ambient."|>,
  "Teledyne Photometrics Kinetix (sCMOS Benchmark)" -> <|
     "ShortName" -> "Kinetix", "Type" -> "sCMOS",
     "PixelsX" -> 3200, "PixelsY" -> 3200, "PixelPitch" -> 6.5,
     "QE" -> 0.91, "ReadNoiseRMS" -> 1.2,
     "ReadNoiseModes" -> <|"Sensitivity" -> 1.2|>,
     "DarkCurrent" -> 0.15, "FullWell" -> 15000,
     "BinningType" -> "Software", "HardwareBinLimit" -> 1, "SupportedBins" -> {1, 2, 4, 8}, "ENF" -> 1.0,
     "Notes" -> "FullWell is a nominal data-sheet value."|>,
  "Andor iXon Ultra 888 (EMCCD Benchmark)" -> <|
     "ShortName" -> "iXon Ultra 888", "Type" -> "EMCCD",
     "PixelsX" -> 1024, "PixelsY" -> 1024, "PixelPitch" -> 13.0,
     "QE" -> 0.92, "ReadNoiseRMS" -> 0.15,
     "ReadNoiseModes" -> <|"EM Gain 300x" -> 0.15|>,
     "DarkCurrent" -> 0.0005, "FullWell" -> 80000,
     "BinningType" -> "Hardware", "HardwareBinLimit" -> Infinity, "SupportedBins" -> {1, 2, 4, 8}, "ENF" -> Sqrt[2.],
     "Notes" -> "Effective read noise 0.15 e- at EM gain 300; excess noise factor Sqrt[2]. FullWell is the conventional-mode value."|>,
  "Andor CB2 High Res (24.5 MP BSI sCMOS)" -> <|
     "ShortName" -> "CB2 High Res", "Type" -> "sCMOS (BSI, global shutter)",
     "PixelsX" -> 5328, "PixelsY" -> 4608, "PixelPitch" -> 2.74,
     "QE" -> 0.62, "ReadNoiseRMS" -> 1.4,
     "ReadNoiseModes" -> <|"12-bit, 24 dB gain" -> 1.4|>,
     "DarkCurrent" -> 0.043, "FullWell" -> 9500,
     "BinningType" -> "On-chip 2\[Times]2 + software", "HardwareBinLimit" -> 2, "SupportedBins" -> {1, 2, 4, 8}, "ENF" -> 1.0,
     "Notes" -> "Andor CB2 brochure: QE 62% at 606 nm read from the published curve (peak 74% at 500 nm); read noise 1.4 e- in 12-bit at 24 dB analogue gain; dark current 0.043 e-/p/s at +20 C (lower when cooled to -5 C air / -40 C liquid); full well 9.5 ke- at lowest gain (16-bit HDR mode combines both). 2x2 charge-domain (floating-diffusion) binning on chip, further binning is digital."|>,
  "Andor CB2 High Speed 7.1F (sCMOS)" -> <|
     "ShortName" -> "CB2 HS 7.1F", "Type" -> "sCMOS (FSI, global shutter)",
     "PixelsX" -> 3216, "PixelsY" -> 2208, "PixelPitch" -> 4.5,
     "QE" -> 0.68, "ReadNoiseRMS" -> 1.4,
     "ReadNoiseModes" -> <|"12-bit, 24 dB gain" -> 1.4|>,
     "DarkCurrent" -> 0.26, "FullWell" -> 23000,
     "BinningType" -> "On-chip 2\[Times]2 + software", "HardwareBinLimit" -> 2, "SupportedBins" -> {1, 2, 4, 8}, "ENF" -> 1.0,
     "Notes" -> "Andor CB2 brochure: QE 68% at 606 nm read from the published curve (peak 74% at 550 nm); read noise 1.4 e- in 12-bit at 24 dB analogue gain; dark current 0.26 e-/p/s at +20 C; full well 23 ke- at lowest gain (14-bit HDR mode combines both). 2x2 charge-domain binning on chip gives 9 um superpixels at the same 1.4 e-, which is why the 1.7F/0.5F models (9 um, 2.6 e-) are not listed separately."|>
|>;
cameraNames = Keys[cameraDatabase];
cameraShortName[name_String] := cameraDatabase[name]["ShortName"];
cameraColors = AssociationThread[cameraNames, ColorData[97] /@ Range[Length[cameraNames]]];
Dataset[KeyDrop[#, {"Notes", "ReadNoiseModes"}] & /@ cameraDatabase]


(* ::Subsection:: *)
(*1.4 Superpixel noise terms and integrated-aperture SNR*)


(* ::Text:: *)
(*For a b x b superpixel: backgroundElectrons = flux x area x tExp x QE, darkElectrons = darkCurrent x b^2 x tExp. Hardware (on-chip) binning reads the superpixel once, so the effective read noise is unchanged; software binning sums b^2 independently read pixels, so the effective read-noise variance is b^2 sigmaRead^2 (effective rms = b sigmaRead). The general rule, sigmaReadEff = sigmaRead x b / Min[b, hardwareBinLimit], also covers sensors such as the Andor CB2 that bin 2x2 in the charge domain on chip and digitally beyond that.*)


(* ::Text:: *)
(*The aperture-integrated SNR uses signal electrons S = nSig QE and the variance ENF^2 (S + 2 nSuper B + 2 nSuper D) + 2 nSuper sigmaReadEff^2. The factor 2 on the background, dark and read terms accounts for subtracting an identical off-target reference aperture. The minimum signal for a target SNR follows from the quadratic S^2 - SNR^2 ENF^2 S - SNR^2 V0 = 0 in closed form.*)


(* ::Input:: *)
sensorSuperpixelPitch[cam_Association, binFactor_] := binFactor cam["PixelPitch"];                    (* um on the chip *)
sensorSuperpixelArea[cam_Association, binFactor_] := sensorSuperpixelPitch[cam, binFactor]^2;         (* um^2 on the chip *)
superpixelPitch[cam_Association, binFactor_] := sensorSuperpixelPitch[cam, binFactor]/magnification;  (* um projected back to the MOT *)
superpixelArea[cam_Association, binFactor_] := superpixelPitch[cam, binFactor]^2;                     (* um^2 at the MOT *)
superpixelCount[cam_Association, binFactor_] := Round[apertureAreaMicron/superpixelArea[cam, binFactor]];
backgroundElectrons[cam_Association, binFactor_, tExp_, multiplier_] :=
  backgroundFluxDensity[multiplier] sensorSuperpixelArea[cam, binFactor] tExp cam["QE"];   (* scatter is collected per unit sensor area, so no magnification factor *)
darkElectrons[cam_Association, binFactor_, tExp_] := cam["DarkCurrent"] binFactor^2 tExp;
effectiveReadNoise[cam_Association, binFactor_, readNoiseRMS_] :=
  readNoiseRMS binFactor/Min[binFactor, cam["HardwareBinLimit"]];   (* one read per on-chip superpixel up to HardwareBinLimit (1 = digital only, Infinity = CCD serial-register binning, 2 = CB2 charge-domain 2x2), digital summing of b^2/limit^2 reads beyond *)
effectiveReadNoise[cam_Association, binFactor_] := effectiveReadNoise[cam, binFactor, cam["ReadNoiseRMS"]];
noiseFloorVariance[cam_Association, binFactor_, tExp_, multiplier_, readNoiseRMS_] :=
  Module[{nSuper, bgElectrons, dkElectrons, readEff},
    nSuper = superpixelCount[cam, binFactor];
    bgElectrons = backgroundElectrons[cam, binFactor, tExp, multiplier];
    dkElectrons = darkElectrons[cam, binFactor, tExp];
    readEff = effectiveReadNoise[cam, binFactor, readNoiseRMS];
    cam["ENF"]^2 (2 nSuper bgElectrons + 2 nSuper dkElectrons) + 2 nSuper readEff^2];
integratedSNR[cam_Association, binFactor_, nSig_, tExp_, multiplier_, readNoiseRMS_] :=
  Module[{signalElectrons, floorVariance},
    signalElectrons = nSig cam["QE"];
    floorVariance = noiseFloorVariance[cam, binFactor, tExp, multiplier, readNoiseRMS];
    signalElectrons/Sqrt[cam["ENF"]^2 signalElectrons + floorVariance]];
integratedSNR[cam_Association, binFactor_, nSig_] :=
  integratedSNR[cam, binFactor, nSig, exposureTimeNominal, backgroundMultiplierDefault, cam["ReadNoiseRMS"]];
minimumSignalPhotons[cam_Association, binFactor_, targetSNR_, tExp_, multiplier_, readNoiseRMS_] :=
  Module[{enfSquared, floorVariance, signalElectrons},
    enfSquared = cam["ENF"]^2;
    floorVariance = noiseFloorVariance[cam, binFactor, tExp, multiplier, readNoiseRMS];
    signalElectrons = (targetSNR^2 enfSquared + Sqrt[targetSNR^4 enfSquared^2 + 4 targetSNR^2 floorVariance])/2;
    signalElectrons/cam["QE"]];
minimumSignalPhotons[cam_Association, binFactor_, targetSNR_] :=
  minimumSignalPhotons[cam, binFactor, targetSNR, exposureTimeNominal, backgroundMultiplierDefault, cam["ReadNoiseRMS"]];
noiseFloorBreakdown[cam_Association, binFactor_, tExp_, multiplier_, readNoiseRMS_] :=
  Module[{nSuper, readVariance, scatterVariance, darkVariance, total},
    nSuper = superpixelCount[cam, binFactor];
    readVariance = 2 nSuper effectiveReadNoise[cam, binFactor, readNoiseRMS]^2;
    scatterVariance = cam["ENF"]^2 2 nSuper backgroundElectrons[cam, binFactor, tExp, multiplier];
    darkVariance = cam["ENF"]^2 2 nSuper darkElectrons[cam, binFactor, tExp];
    total = readVariance + scatterVariance + darkVariance;
    <|"Superpixels" -> nSuper, "ReadVariance" -> readVariance, "ScatterVariance" -> scatterVariance,
      "DarkVariance" -> darkVariance, "FloorVariance" -> total,
      "DominantFloorTerm" -> First[Keys[TakeLargest[<|"Read noise" -> readVariance, "Laser scatter" -> scatterVariance, "Dark current" -> darkVariance|>, 1]]],
      "ShotNoiseCrossoverPhotons" -> total/(cam["ENF"]^2 cam["QE"])|>];
localLogSlope[cam_Association, binFactor_, nSig_] :=
  Module[{ratio = 1.01},
    (Log[integratedSNR[cam, binFactor, nSig ratio]] - Log[integratedSNR[cam, binFactor, nSig/ratio]])/(2 Log[ratio])];
integratedSNRMolecules[cam_Association, binFactor_, nMol_, tExp_, multiplier_, readNoiseRMS_, rate_] :=
  integratedSNR[cam, binFactor, signalPhotonsFromMolecules[nMol, tExp, rate], tExp, multiplier, readNoiseRMS];
integratedSNRMolecules[cam_Association, binFactor_, nMol_] :=
  integratedSNRMolecules[cam, binFactor, nMol, exposureTimeNominal, backgroundMultiplierDefault, cam["ReadNoiseRMS"], scatteringRate];
minimumMolecules[cam_Association, binFactor_, targetSNR_, tExp_, multiplier_, readNoiseRMS_, rate_] :=
  moleculesFromSignalPhotons[minimumSignalPhotons[cam, binFactor, targetSNR, tExp, multiplier, readNoiseRMS], tExp, rate];
minimumMolecules[cam_Association, binFactor_, targetSNR_] :=
  minimumMolecules[cam, binFactor, targetSNR, exposureTimeNominal, backgroundMultiplierDefault, cam["ReadNoiseRMS"], scatteringRate];
peakSuperpixelFraction[cam_Association, binFactor_] :=
  Module[{pitchMm = superpixelPitch[cam, binFactor]/1000.}, Erf[pitchMm/(2 Sqrt[2.] sigmaCloud)]^2];   (* fraction of the cloud in the superpixel centred on the cloud, using its object-space footprint *)
peakSuperpixelElectrons[cam_Association, binFactor_, nSig_] := nSig cam["QE"]/apertureSignalFraction peakSuperpixelFraction[cam, binFactor];
saturationMolecules[cam_Association, binFactor_, tExp_, rate_] :=
  cam["FullWell"]/peakSuperpixelElectrons[cam, binFactor, photonsPerMolecule[tExp, rate]];   (* molecule number at which the peak superpixel mean reaches full well *)
saturationMolecules[cam_Association, binFactor_] := saturationMolecules[cam, binFactor, exposureTimeNominal, scatteringRate];
(* single-shot visibility: signal to noise of the brightest superpixel, which is what decides whether the cloud shows up in a raw frame *)
peakPixelSNR[cam_Association, binFactor_, nMol_, tExp_, multiplier_, readNoiseRMS_, rate_] :=
  Module[{peakElectrons, floorVariance},
    peakElectrons = peakSuperpixelElectrons[cam, binFactor, signalPhotonsFromMolecules[nMol, tExp, rate]];
    floorVariance = cam["ENF"]^2 2 (backgroundElectrons[cam, binFactor, tExp, multiplier] + darkElectrons[cam, binFactor, tExp]) + 2 effectiveReadNoise[cam, binFactor, readNoiseRMS]^2;
    peakElectrons/Sqrt[cam["ENF"]^2 peakElectrons + floorVariance]];
peakPixelSNR[cam_Association, binFactor_, nMol_] :=
  peakPixelSNR[cam, binFactor, nMol, exposureTimeNominal, backgroundMultiplierDefault, cam["ReadNoiseRMS"], scatteringRate];
visibilityMolecules[cam_Association, binFactor_, targetSNR_, tExp_, multiplier_, readNoiseRMS_, rate_] :=
  Module[{enfSquared, floorVariance, peakElectrons, electronsPerMolecule},
    enfSquared = cam["ENF"]^2;
    floorVariance = enfSquared 2 (backgroundElectrons[cam, binFactor, tExp, multiplier] + darkElectrons[cam, binFactor, tExp]) + 2 effectiveReadNoise[cam, binFactor, readNoiseRMS]^2;
    peakElectrons = (targetSNR^2 enfSquared + Sqrt[targetSNR^4 enfSquared^2 + 4 targetSNR^2 floorVariance])/2;
    electronsPerMolecule = peakSuperpixelElectrons[cam, binFactor, photonsPerMolecule[tExp, rate]];
    peakElectrons/electronsPerMolecule];
visibilityMolecules[cam_Association, binFactor_, targetSNR_] :=
  visibilityMolecules[cam, binFactor, targetSNR, exposureTimeNominal, backgroundMultiplierDefault, cam["ReadNoiseRMS"], scatteringRate];
roundSignificant[x_?NumericQ, digits_Integer] := If[x == 0, 0, N[Round[x, 10^(Floor[Log10[Abs[x]]] - digits + 1)]]];
Column[{Row[{"Example: ORCA-Quest 2, 1x1, 10^4 photons, nominal conditions -> SNR = ",
             integratedSNR[cameraDatabase[cameraNames[[1]]], 1, 10^4]}],
        Row[{"Example: ORCA-Quest 2, 1x1, 100 molecules, nominal conditions -> SNR = ",
             integratedSNRMolecules[cameraDatabase[cameraNames[[1]]], 1, 100]}],
        Row[{"Example: ORCA-Quest 2, 1x1 saturates (peak superpixel = full well) at ",
             roundSignificant[saturationMolecules[cameraDatabase[cameraNames[[1]]], 1], 3], " molecules"}],
        Row[{"Example: ORCA-R2 effective read noise at 4x4 (hardware) = ",
             effectiveReadNoise[cameraDatabase[cameraNames[[2]]], 4], " e-"}],
        Row[{"Example: Zelux effective read noise at 4x4 (software) = ",
             effectiveReadNoise[cameraDatabase[cameraNames[[3]]], 4], " e-"}]}]


(* ::Subsection:: *)
(*1.5 How the SNR is defined: per pixel, top-hat aperture, or optimally weighted*)


(* ::Text:: *)
(*"The SNR" is not a property of the camera alone: it depends on the estimator applied to the image. Three definitions are computed here, all for a signal frame minus an identical off-target reference frame.*)


(* ::Text:: *)
(*peakPixelSNR: signal to noise of the single brightest superpixel. This sets whether the cloud is visible by eye in a single shot, but it is NOT the precision of a molecule-number measurement.*)


(* ::Text:: *)
(*topHatSNR: the flat sum over the 2 sigma aperture minus a reference aperture, i.e. the definition used in the benchmark sections. It is the natural estimator of the molecule number and is what most groups quote.*)


(* ::Text:: *)
(*matchedFilterSNR: the same measurement made with the optimal (inverse-variance, profile-weighted) estimator, N proportional to Sum[wi di] with wi proportional to pi/sigmai^2. Its variance is the Cramer-Rao bound 1/Sum[pi^2/sigmai^2], so no unbiased estimator of the molecule number can do better. In the shot-noise limit it coincides with the top-hat sum; when read noise or scatter dominates it wins, because pixels in the wings carry full noise but little signal and the flat sum weights them equally.*)


(* ::Text:: *)
(*Also computed: widthPrecision, the fractional uncertainty on the cloud width from the Fisher information, delta sigma / sigma = 1/(sigma Sqrt[Sum[(dsi/dsigma)^2/sigmai^2]]). This is the figure of merit for time-of-flight thermometry and behaves differently from the number precision: it depends on resolving the profile, so it rewards small pixels and penalises coarse binning.*)


(* ::Input:: *)
estimatorFOV = 6 sigmaCloud;   (* mm at the MOT: +-3 sigma, which holds 99% of a 2D Gaussian; the truncated wings change the estimators by well under 1% *)
gaussianFractions1D[centers_List, pitchMm_, sigmaMm_] :=
  (Erf[(centers + pitchMm/2)/(Sqrt[2.] sigmaMm)] - Erf[(centers - pitchMm/2)/(Sqrt[2.] sigmaMm)])/2;
estimatorCoordinates[cam_Association, binFactor_] :=
  Module[{pitchMm, count},
    pitchMm = superpixelPitch[cam, binFactor]/1000.;
    count = Max[Round[estimatorFOV/pitchMm], 9];
    {(Range[count] - (count + 1)/2.) pitchMm, pitchMm}];
signalFractionGrid[cam_Association, binFactor_, sigmaMm_] :=
  Module[{coords, pitchMm, fractions},
    {coords, pitchMm} = estimatorCoordinates[cam, binFactor];
    fractions = gaussianFractions1D[coords, pitchMm, sigmaMm];
    Outer[Times, fractions, fractions]];
snrEstimators[cam_Association, binFactor_, nMol_, tExp_, multiplier_, readNoiseRMS_, rate_] :=
  Module[{coords, pitchMm, fractionGrid, totalElectrons, signalGrid, varianceGrid, floorPerPixel, apertureMask,
          peakSignal, peakVariance, topHat, matched, fisherWidth, derivativeGrid, delta, nSig},
    nSig = signalPhotonsFromMolecules[nMol, tExp, rate];
    totalElectrons = nSig cam["QE"]/apertureSignalFraction;   (* electrons from the whole cloud, not only the aperture *)
    {coords, pitchMm} = estimatorCoordinates[cam, binFactor];
    fractionGrid = signalFractionGrid[cam, binFactor, sigmaCloud];
    signalGrid = totalElectrons fractionGrid;
    floorPerPixel = cam["ENF"]^2 2 (backgroundElectrons[cam, binFactor, tExp, multiplier] + darkElectrons[cam, binFactor, tExp]) + 2 effectiveReadNoise[cam, binFactor, readNoiseRMS]^2;
    varianceGrid = cam["ENF"]^2 signalGrid + floorPerPixel;
    apertureMask = UnitStep[apertureRadius - Sqrt[Outer[Plus, coords^2, coords^2]]];
    topHat = Total[signalGrid apertureMask, 2]/Sqrt[Total[varianceGrid apertureMask, 2]];
    matched = totalElectrons Sqrt[Total[fractionGrid^2/varianceGrid, 2]];
    peakSignal = Max[signalGrid];
    peakVariance = cam["ENF"]^2 peakSignal + floorPerPixel;
    delta = 0.01 sigmaCloud;
    derivativeGrid = totalElectrons (signalFractionGrid[cam, binFactor, sigmaCloud + delta] - signalFractionGrid[cam, binFactor, sigmaCloud - delta])/(2 delta);
    fisherWidth = Total[derivativeGrid^2/varianceGrid, 2];
    <|"PeakPixelElectrons" -> peakSignal, "PeakPixelSNR" -> peakSignal/Sqrt[peakVariance],
      "TopHatSNR" -> topHat, "MatchedFilterSNR" -> matched, "MatchedGain" -> matched/topHat,
      "NumberPrecision" -> 1/matched, "WidthPrecision" -> 1/(sigmaCloud Sqrt[fisherWidth]),
      "PixelsInAperture" -> Total[apertureMask, 2]|>];
snrEstimators[cam_Association, binFactor_, nMol_] :=
  snrEstimators[cam, binFactor, nMol, exposureTimeNominal, backgroundMultiplierDefault, cam["ReadNoiseRMS"], scatteringRate];


(* ::Section:: *)
(*2. Primary Benchmark: Integrated SNR vs Number of Molecules in the MOT (log-log)*)


(* ::Text:: *)
(*All five cameras at 1x1, 2x2 and 4x4 binning, 20 ms exposure, 4x baseline scatter, 10 to 10^6 molecules in the MOT (about 607 signal photons in the aperture per molecule at 20 ms and the assumed 1e6 /s scattering rate). Red dashed lines mark SNR = 3 (detection threshold) and SNR = 10 (quantitative sizing limit). Grey guide lines show the two asymptotic regimes: where the signal-independent noise floor (read noise, laser scatter, dark current) dominates, SNR is proportional to nMol (slope 1 in log-log); once the signal shot noise dominates, SNR is proportional to Sqrt[nMol] (slope 1/2). Filled markers show where the peak superpixel reaches the full-well capacity at 20 ms: beyond that point the curve is only reachable with a shorter exposure or coarser binning.*)


(* ::Text:: *)
(*Note on binning: binning sums pixels, it does not average them. Over a fixed aperture the total signal, background and dark electrons are the same electrons however they are grouped, so their shot-noise variances are untouched by binning; only the read noise depends on the grouping, because it is paid once per readout. Hardware binning performs one read per superpixel, so the integrated read-noise variance 2 nSuper sigmaRead^2 falls as 1/b^2. Software binning still reads every physical pixel and adds afterwards, giving 2 nSuper (b sigmaRead)^2 = 2 apertureArea sigmaRead^2 / pixelPitch^2, which is exactly independent of b. Hence the 1x1, 2x2 and 4x4 curves of the qCMOS, sCMOS and CMOS cameras coincide to within the rounding of nSuper, while the ORCA-R2 and iXon separate. Binning a software-binned camera is still worth doing for per-pixel visibility, data volume and frame rate, but it adds no information about the molecule number: see Section 3, where the matched-filter SNR of the Zelux is identical at every binning while its peak-pixel SNR grows as b.*)


(* ::Text:: *)
(*Caveat on the laser-scatter background: it is specified as a flux per sensor pixel area and implemented that way, so the scatter collected inside the aperture scales as M^2 and demagnification reduces it. That is the right model for stray light reaching the sensor from outside the imaged region. If instead the scatter is light from the imaged MOT region (chamber walls, viewports in the field of view), it is imaged like the signal, its total inside the aperture is independent of M, and only the read-noise term benefits from demagnification. The distinction does not affect the read-noise-limited cameras (for the Zelux the scatter is 3.6% of the noise-floor variance) but it matters for the ORCA-Quest 2 and the iXon, which are scatter-dominated. Measure the background per pixel per frame to settle it.*)


(* ::Input:: *)
moleculeMin = First[moleculeNumberRange];
moleculeMax = Last[moleculeNumberRange];
benchmarkBins = {1, 2, 4};
binDashing = <|1 -> Dashing[{}], 2 -> Dashing[{0.025, 0.012}], 4 -> Dashing[{0.005, 0.010}]|>;
binLabel[binFactor_Integer] := ToString[binFactor] <> "\[Times]" <> ToString[binFactor];
benchmarkCurveKeys = Flatten[Table[{name, binFactor}, {name, cameraNames}, {binFactor, benchmarkBins}], 1];
benchmarkCurveStyles = (Directive[cameraColors[#1], binDashing[#2], AbsoluteThickness[1.8]] &) @@@ benchmarkCurveKeys;
saturationMarker[name_String, binFactor_] :=
  Module[{cam = cameraDatabase[name], nMolSat},
    nMolSat = saturationMolecules[cam, binFactor];
    If[moleculeMin <= nMolSat <= moleculeMax,
      {EdgeForm[Directive[Black, AbsoluteThickness[0.8]]], cameraColors[name], Disk[{Log[nMolSat], Log[integratedSNRMolecules[cam, binFactor, nMolSat]]}, Offset[{4.5, 4.5}]]},
      {}]];
saturationMarkers = saturationMarker @@@ benchmarkCurveKeys;
asymptoteAnnotations = {
   Directive[Gray, AbsoluteThickness[1], Dashing[{0.01, 0.006}]],
   Line[{{Log[12.], Log[0.6]}, {Log[300.], Log[15.]}}],
   Line[{{Log[2.*^4], Log[3.*^4]}, {Log[6.*^5], Log[3.*^4 Sqrt[6.*^5/2.*^4]]}}],
   Text[Style["Noise-floor limited\n(read noise + scatter + dark)\nslope \[TildeTilde] 1", 10, Gray, FontFamily -> "Helvetica"], {Log[120.], Log[1.6]}, {-1, 0}],
   Text[Style["Signal shot-noise limited\nslope \[TildeTilde] 1/2", 10, Gray, FontFamily -> "Helvetica"], {Log[1.5*^4], Log[1.5*^5]}, {-1, 0}]};
benchmarkPlot = LogLogPlot[
   Evaluate[(integratedSNRMolecules[cameraDatabase[#1], #2, nMol] &) @@@ benchmarkCurveKeys], {nMol, moleculeMin, moleculeMax},
   PlotStyle -> benchmarkCurveStyles,
   PlotRange -> {{moleculeMin, moleculeMax}, {0.5, 3*^5}},
   Frame -> True, Axes -> False,
   FrameLabel -> {"Number of molecules in the MOT, nMol", "Integrated aperture SNR"},
   PlotLabel -> Style["CaF MOT fluorescence: integrated SNR at 20 ms exposure, 4\[Times] baseline scatter, 10^6 photons/s per molecule", 13],
   GridLines -> {Automatic, {{3, Directive[Red, Dashed, AbsoluteThickness[1.5]]}, {10, Directive[Darker[Red], Dashed, AbsoluteThickness[1.5]]}}},
   Epilog -> Join[asymptoteAnnotations, saturationMarkers,
                  {Text[Style["SNR = 3 (detection)", 10, Red], {Log[1.1*^4], Log[3.]}, {-1, -1.2}],
                   Text[Style["SNR = 10 (sizing)", 10, Darker[Red]], {Log[1.1*^4], Log[10.]}, {-1, -1.2}]}],
   ImageSize -> 720, AspectRatio -> 0.65, PerformanceGoal -> "Quality"];
cameraLegend = LineLegend[Values[cameraColors], cameraShortName /@ cameraNames, LegendLabel -> "Camera", LegendLayout -> "Column"];
binningLegend = LineLegend[(Directive[Black, binDashing[#], AbsoluteThickness[1.8]] &) /@ benchmarkBins, ("b = " <> binLabel[#] &) /@ benchmarkBins, LegendLabel -> "Binning"];
saturationLegend = PointLegend[{Directive[Gray, EdgeForm[Black]]}, {"peak superpixel = full well"}, LegendMarkers -> {Graphics[{EdgeForm[Black], Gray, Disk[]}]}, LegendMarkerSize -> 10, LegendLabel -> "Saturation (20 ms)"];
benchmarkFigure = Legended[benchmarkPlot, Placed[Column[{cameraLegend, binningLegend, saturationLegend}, Spacings -> 1.5], Right]]


(* ::Text:: *)
(*The same curves split by camera, which makes the hardware-versus-software binning behaviour visible.*)


(* ::Input:: *)
benchmarkPanel[name_String] := LogLogPlot[
   Evaluate[integratedSNRMolecules[cameraDatabase[name], #, nMol] & /@ benchmarkBins], {nMol, moleculeMin, moleculeMax},
   PlotStyle -> (Directive[cameraColors[name], binDashing[#], AbsoluteThickness[1.8]] & /@ benchmarkBins),
   PlotRange -> {{moleculeMin, moleculeMax}, {0.5, 3*^5}}, Frame -> True, Axes -> False,
   FrameLabel -> {"nMol", "SNR"}, PlotLabel -> Style[cameraShortName[name] <> " (" <> cameraDatabase[name]["BinningType"] <> " binning)", 11],
   GridLines -> {None, {{3, Directive[Red, Dashed]}, {10, Directive[Darker[Red], Dashed]}}},
   Epilog -> (saturationMarker[name, #] & /@ benchmarkBins),
   PlotLegends -> Placed[LineLegend[binLabel /@ benchmarkBins], {0.25, 0.8}], ImageSize -> 300];
benchmarkGrid = Multicolumn[benchmarkPanel /@ cameraNames, 3, Appearance -> "Horizontal"]


(* ::Text:: *)
(*Asymptotic diagnostics at nominal conditions: which term sets the noise floor, the molecule number at which signal shot noise overtakes the floor, the molecule number at which the peak superpixel saturates (20 ms), and the local log-log slope at three molecule numbers.*)


(* ::Input:: *)
regimeDataset = Dataset[Flatten[Table[
    Module[{cam = cameraDatabase[name], breakdown},
      breakdown = noiseFloorBreakdown[cam, binFactor, exposureTimeNominal, backgroundMultiplierDefault, cam["ReadNoiseRMS"]];
      <|"Camera" -> cameraShortName[name], "Bin" -> binLabel[binFactor],
        "nSuper" -> breakdown["Superpixels"],
        "Read var (e-^2)" -> roundSignificant[breakdown["ReadVariance"], 3],
        "Scatter var (e-^2)" -> roundSignificant[breakdown["ScatterVariance"], 3],
        "Dark var (e-^2)" -> roundSignificant[breakdown["DarkVariance"], 3],
        "Dominant floor" -> breakdown["DominantFloorTerm"],
        "Shot-noise crossover nMol" -> roundSignificant[moleculesFromSignalPhotons[breakdown["ShotNoiseCrossoverPhotons"], exposureTimeNominal], 3],
        "Saturation nMol" -> roundSignificant[saturationMolecules[cam, binFactor], 3],
        "Slope @10 mol" -> roundSignificant[localLogSlope[cam, binFactor, signalPhotonsFromMolecules[10., exposureTimeNominal]], 3],
        "Slope @1e3 mol" -> roundSignificant[localLogSlope[cam, binFactor, signalPhotonsFromMolecules[10.^3, exposureTimeNominal]], 3],
        "Slope @1e5 mol" -> roundSignificant[localLogSlope[cam, binFactor, signalPhotonsFromMolecules[10.^5, exposureTimeNominal]], 3]|>],
    {name, cameraNames}, {binFactor, benchmarkBins}], 1]]


(* ::Section:: *)
(*3. Choosing the SNR Definition: Per Pixel vs Top-Hat vs Optimal Weighting*)


(* ::Text:: *)
(*At the working point (10^4 molecules, 20 ms, 4x scatter) the three definitions differ by orders of magnitude, so quoting "the SNR" without saying which one is meaningless. The table gives all three, the gain from optimal weighting, and the resulting fractional precision on the molecule number and on the cloud width.*)


(* ::Input:: *)
workingPointMolecules = 10^4;
estimatorDataset = Dataset[Flatten[Table[
    Module[{cam = cameraDatabase[name], est},
      est = snrEstimators[cam, binFactor, workingPointMolecules];
      <|"Camera" -> cameraShortName[name], "Bin" -> binLabel[binFactor],
        "Object pitch (um)" -> roundSignificant[superpixelPitch[cam, binFactor], 3],
        "Peak e-/pixel" -> roundSignificant[est["PeakPixelElectrons"], 3],
        "Peak-pixel SNR" -> roundSignificant[est["PeakPixelSNR"], 3],
        "Top-hat SNR" -> roundSignificant[est["TopHatSNR"], 4],
        "Matched-filter SNR" -> roundSignificant[est["MatchedFilterSNR"], 4],
        "Optimal gain" -> roundSignificant[est["MatchedGain"], 3],
        "dN/N (%)" -> roundSignificant[100 est["NumberPrecision"], 3],
        "dWidth/Width (%)" -> roundSignificant[100 est["WidthPrecision"], 3]|>],
    {name, cameraNames}, {binFactor, benchmarkBins}], 1]]


(* ::Text:: *)
(*How the three definitions separate with molecule number, at 1x1 binning. The gain from optimal weighting is largest exactly where the measurement is hardest (read-noise limited, few molecules) and vanishes once signal shot noise dominates, where the flat sum is already optimal.*)


(* ::Input:: *)
estimatorSampleMolecules = 10.^Range[1, 6, 0.25];
estimatorPanel[name_String, binFactor_] :=
  Module[{samples, series, labels},
    samples = snrEstimators[cameraDatabase[name], binFactor, #] & /@ estimatorSampleMolecules;
    labels = {"TopHatSNR", "MatchedFilterSNR", "PeakPixelSNR"};
    series = Table[Transpose[{estimatorSampleMolecules, #[key] & /@ samples}], {key, labels}];
    ListLogLogPlot[series, Joined -> True,
      PlotStyle -> {Directive[cameraColors[name], AbsoluteThickness[1.8]],
                    Directive[Black, AbsoluteThickness[1.6], Dashing[{0.02, 0.012}]],
                    Directive[Gray, AbsoluteThickness[1.2], Dashing[{0.004, 0.008}]]},
      Frame -> True, Axes -> False, PlotRange -> {{moleculeMin, moleculeMax}, {0.05, 3*^5}},
      FrameLabel -> {"Molecules in the MOT", "SNR"},
      PlotLabel -> Style[cameraShortName[name] <> ", " <> binLabel[binFactor], 11],
      GridLines -> {None, {{3, Directive[Red, Dashed]}, {10, Directive[Darker[Red], Dashed]}}},
      PlotLegends -> Placed[LineLegend[{"top-hat aperture", "matched filter (optimal)", "peak pixel"}], Below],
      ImageSize -> 330]];
estimatorComparisonPlot = Grid[{{estimatorPanel[cameraNames[[3]], 1], estimatorPanel[cameraNames[[2]], 1], estimatorPanel[cameraNames[[1]], 1]}}, Spacings -> {1, 0}]


(* ::Section:: *)
(*4. Single-Shot Visibility: Peak-Pixel SNR vs Binning*)


(* ::Text:: *)
(*This section is about seeing the cloud in a raw single shot, not about measuring the molecule number. The relevant quantity is the signal to noise of the brightest superpixel: the cloud is marginal at peakPixelSNR ~ 3 and unmistakable at ~5, whereas the aperture-integrated SNR of Section 2 can be in the hundreds while no individual pixel stands out.*)


(* ::Text:: *)
(*This is the one place where software binning genuinely helps. Binning b x b multiplies the signal per superpixel by b^2 and its noise by b, so peakPixelSNR grows in proportion to b for every camera, hardware-binned or not. It buys visibility, not precision: the integrated SNR is unchanged (software) or improved separately (hardware), as Section 3 shows.*)


(* ::Input:: *)
visibilityThreshold = 5;      (* peak-pixel SNR at which the cloud is unmistakable in a raw frame *)
visibilityThresholdMarginal = 3;
visibilityPanel[name_String] :=
  Module[{cam = cameraDatabase[name], bins},
    cam = cameraDatabase[name]; bins = cam["SupportedBins"];
    LogLogPlot[Evaluate[peakPixelSNR[cam, #, nMol] & /@ bins], {nMol, moleculeMin, moleculeMax},
      PlotStyle -> (Directive[cameraColors[name], AbsoluteThickness[1.6], Opacity[0.35 + 0.65 (Position[bins, #][[1, 1]] - 1)/(Length[bins] - 1)]] & /@ bins),
      PlotRange -> {{moleculeMin, moleculeMax}, {0.02, 2000}}, Frame -> True, Axes -> False,
      FrameLabel -> {"Molecules in the MOT", "Peak-pixel SNR"},
      PlotLabel -> Style[cameraShortName[name] <> " (" <> cam["BinningType"] <> " binning)", 11],
      GridLines -> {{{workingPointMolecules, Directive[Blue, Dashed, Opacity[0.5]]}},
                    {{visibilityThreshold, Directive[Red, Dashed, AbsoluteThickness[1.5]]}, {visibilityThresholdMarginal, Directive[Orange, Dashed]}}},
      PlotLegends -> Placed[LineLegend[binLabel /@ bins, LegendLabel -> "bin"], {0.83, 0.26}],
      ImageSize -> 330]];
visibilityGrid = Multicolumn[visibilityPanel /@ cameraNames, 3, Appearance -> "Horizontal"]


(* ::Text:: *)
(*The actionable summary: how many molecules are needed before the cloud is visible in one shot, as a function of binning. Lower is better. The blue line is the working point of 10^4 molecules; any camera/binning combination whose curve lies below it shows the MOT in a raw frame.*)


(* ::Input:: *)
visibilityThresholdPlot = Module[{series, plotOptions},
   series = Table[Table[{binFactor, visibilityMolecules[cameraDatabase[name], binFactor, visibilityThreshold]}, {binFactor, cameraDatabase[name]["SupportedBins"]}], {name, cameraNames}];
   ListPlot[series, PlotStyle -> Values[cameraColors], PlotMarkers -> {Automatic, 9}, Joined -> True,
     Frame -> True, Axes -> False, ScalingFunctions -> {"Log2", "Log"},
     FrameTicks -> {{Automatic, None}, {{1, 2, 4, 8, 16}, None}},
     FrameLabel -> {"Binning factor b", "Molecules needed for peak-pixel SNR = 5"},
     PlotLabel -> Style["Single-shot visibility threshold (20 ms, 4\[Times] baseline scatter)", 12],
     GridLines -> {None, {{workingPointMolecules, Directive[Blue, Dashed]}}},
     Epilog -> {Text[Style["working point, 10^4 molecules", 9, Blue], {Log2[1.1], Log[workingPointMolecules]}, {-1, -1.3}]},
     PlotLegends -> Placed[LineLegend[cameraShortName /@ cameraNames], Right], ImageSize -> 560]]


(* ::Input:: *)
visibilityDataset = Dataset[Flatten[Table[
    Module[{cam = cameraDatabase[name]},
      <|"Camera" -> cameraShortName[name], "Binning" -> cam["BinningType"], "Bin" -> binLabel[binFactor],
        "Object pitch (um)" -> roundSignificant[superpixelPitch[cam, binFactor], 3],
        "Elements across FWHM" -> Round[resolutionElementsAcrossFWHM[cam, binFactor]],
        "Peak SNR at 10^4 mol" -> roundSignificant[peakPixelSNR[cam, binFactor, workingPointMolecules], 3],
        "Molecules for peak SNR = 3" -> roundSignificant[visibilityMolecules[cam, binFactor, visibilityThresholdMarginal], 3],
        "Molecules for peak SNR = 5" -> roundSignificant[visibilityMolecules[cam, binFactor, visibilityThreshold], 3],
        "Visible at 10^4?" -> If[peakPixelSNR[cam, binFactor, workingPointMolecules] >= visibilityThreshold, "yes", "no"]|>],
    {name, cameraNames}, {binFactor, cameraDatabase[name]["SupportedBins"]}], 1]]


(* ::Section:: *)
(*5. Spatial Resolution vs Binning Factor*)


(* ::Text:: *)
(*With magnification M = 0.478 a superpixel of physical pitch b x pixelPitch samples b x pixelPitch / M at the MOT, so the demagnification coarsens the object-space sampling by a factor 2.1 relative to a 1:1 relay. The number of resolution elements across the MOT FWHM (2.8 mm) is FWHM / object-space superpixel pitch; at least ~5 elements are needed to fit a Gaussian width reliably.*)


(* ::Input:: *)
resolutionElementsAcrossFWHM[cam_Association, binFactor_] := 1000 cloudFWHM/superpixelPitch[cam, binFactor];
resolutionDataset = Dataset[Flatten[Table[
    Module[{cam = cameraDatabase[name]},
      <|"Camera" -> cameraShortName[name], "Binning" -> cam["BinningType"], "Bin" -> binLabel[binFactor],
        "Sensor pitch (um)" -> sensorSuperpixelPitch[cam, binFactor],
        "Object pitch (um)" -> roundSignificant[superpixelPitch[cam, binFactor], 3],
        "Elements across FWHM" -> roundSignificant[resolutionElementsAcrossFWHM[cam, binFactor], 4],
        "Superpixels in aperture" -> superpixelCount[cam, binFactor]|>],
    {name, cameraNames}, {binFactor, cameraDatabase[name]["SupportedBins"]}], 1]]


(* ::Input:: *)
resolutionPlots = Module[{superpixelSeries, elementSeries, plotOptions},
   superpixelSeries = Table[Table[{binFactor, superpixelPitch[cameraDatabase[name], binFactor]}, {binFactor, cameraDatabase[name]["SupportedBins"]}], {name, cameraNames}];
   elementSeries = Table[Table[{binFactor, resolutionElementsAcrossFWHM[cameraDatabase[name], binFactor]}, {binFactor, cameraDatabase[name]["SupportedBins"]}], {name, cameraNames}];
   plotOptions = {PlotStyle -> Values[cameraColors], PlotMarkers -> {Automatic, 9}, Joined -> True, Frame -> True, Axes -> False,
                  ScalingFunctions -> {"Log2", "Log"}, FrameTicks -> {{Automatic, None}, {{1, 2, 4, 8, 16}, None}}, ImageSize -> 400, GridLines -> Automatic};
   Grid[{{ListPlot[superpixelSeries, FrameLabel -> {"Binning factor b", "Object-space superpixel size (um)"}, PlotLabel -> "Object-space sampling (M = 0.478)", PlotRange -> {{0.8, 20}, {2, 150}}, Evaluate[Sequence @@ plotOptions]],
        ListPlot[elementSeries, FrameLabel -> {"Binning factor b", "Resolution elements across FWHM"}, PlotLabel -> "Sampling of the 2.8 mm MOT FWHM", PlotRange -> {{0.8, 20}, {3, 400}},
                 Epilog -> {Directive[Red, Dashed], Line[{{-0.2, Log[5]}, {4.2, Log[5]}}], Text[Style["5 elements", 9, Red], {3.4, Log[5]}, {0, -1}]},
                 PlotLegends -> Placed[LineLegend[cameraShortName /@ cameraNames], Below], Evaluate[Sequence @@ plotOptions]]}}, Alignment -> Top, Spacings -> {2, 0}]]


(* ::Section:: *)
(*6. Interactive MOT Dashboard*)


(* ::Text:: *)
(*Controls: camera, readout mode, binning, number of molecules (10 to 10^6), scattering rate per molecule, exposure, background multiplier and noise seed. A synthetic frame is generated for the selected camera and binning over a 2.5 mm x 2.5 mm field of view: the expected electrons per superpixel are the exact Gaussian pixel integrals (via Erf) plus background and dark electrons; the shot noise is Poisson (excess-noise factor applied as a variance multiplier for the EMCCD), and Gaussian read noise with the effective (binned) rms is added. The centre row is fitted with a Gaussian plus offset and the metrics card compares the fit with the analytic model.*)


(* ::Input:: *)
syntheticFrameFOV = 2.6 cloudFWHM;   (* mm at the MOT: 1.5 aperture diameters, leaving room for the off-target reference annulus *)
poissonSample = Compile[{{mu, _Real}},
   Module[{limit, k, p},
     If[mu <= 0., 0.,
       If[mu > 50.,
         Max[0., Round[mu + Sqrt[mu] Sqrt[-2. Log[RandomReal[]]] Cos[2. Pi RandomReal[]]]],
         limit = Exp[-mu]; k = 0; p = 1.;
         While[p > limit, k++; p *= RandomReal[]];
         N[k - 1]]]],
   RuntimeAttributes -> {Listable}, Parallelization -> True];
frameDimensions[cam_Association, binFactor_] :=
  Module[{fovPixels},
    fovPixels = Floor[1000 syntheticFrameFOV/superpixelPitch[cam, binFactor]];
    {Min[fovPixels, Floor[cam["PixelsY"]/binFactor]], Min[fovPixels, Floor[cam["PixelsX"]/binFactor]]}];   (* {rows, columns} *)
frameCoordinates[cam_Association, binFactor_] :=
  Module[{rows, columns, pitchMm},
    {rows, columns} = frameDimensions[cam, binFactor];
    pitchMm = superpixelPitch[cam, binFactor]/1000.;
    {(Range[columns] - (columns + 1)/2.) pitchMm, (Range[rows] - (rows + 1)/2.) pitchMm}];   (* {xs, ys} in mm *)
gaussianPixelFractions[centers_List, pitchMm_] :=
  Module[{scale = Sqrt[2.] sigmaCloud},
    (Erf[(centers + pitchMm/2)/scale] - Erf[(centers - pitchMm/2)/scale])/2];
signalMeanFrame[cam_Association, binFactor_, nSig_] :=
  Module[{xs, ys, pitchMm, totalElectrons, fractionsX, fractionsY},
    {xs, ys} = frameCoordinates[cam, binFactor];
    pitchMm = superpixelPitch[cam, binFactor]/1000.;
    totalElectrons = nSig cam["QE"]/apertureSignalFraction;
    fractionsX = gaussianPixelFractions[xs, pitchMm];
    fractionsY = gaussianPixelFractions[ys, pitchMm];
    totalElectrons Outer[Times, fractionsY, fractionsX]];
syntheticFrame[cam_Association, binFactor_, nSig_, tExp_, multiplier_, readNoiseRMS_] :=
  Module[{meanFrame, shotFrame, readFrame, readEff},
    meanFrame = signalMeanFrame[cam, binFactor, nSig] + backgroundElectrons[cam, binFactor, tExp, multiplier] + darkElectrons[cam, binFactor, tExp];
    shotFrame = poissonSample[meanFrame];
    readEff = effectiveReadNoise[cam, binFactor, readNoiseRMS];
    readFrame = RandomVariate[NormalDistribution[0., readEff], Dimensions[meanFrame]];
    meanFrame + cam["ENF"] (shotFrame - meanFrame) + readFrame];
frameImage[frame_List, xs_List, ys_List, title_String] :=
  Module[{flat = Flatten[frame], low, high},
    {low, high} = Quantile[flat, {0.002, 0.999}];
    If[high <= low, high = low + 1.];
    ArrayPlot[frame, DataRange -> {{First[xs], Last[xs]}, {First[ys], Last[ys]}}, DataReversed -> True,
      ColorFunction -> (Blend[{Black, RGBColor[0.10, 0.15, 0.55], RGBColor[0.85, 0.35, 0.15], RGBColor[1., 0.95, 0.75]}, #] &),
      ClippingStyle -> Automatic, PlotRange -> {All, All, {low, high}},
      Frame -> True, FrameLabel -> {{"x (mm)", None}, {"y (mm)", None}}, PlotLabel -> Style[title, 11],
      PlotLegends -> BarLegend[Automatic, LegendLabel -> "e-/superpixel"],
      Epilog -> {Directive[White, Dashed, AbsoluteThickness[1]], Circle[{0, 0}, apertureRadius]},
      ImageSize -> 360]];
ClearAll[fitAmplitude, fitCentre, fitWidth, fitOffset, fitVariable];   (* formal fit-parameter symbols, never assigned: NonlinearModelFit caches its model symbols, so Module-local ones would persist as Global` temporaries *)
fitCrossSection[xs_List, row_List] :=
  Module[{data, amplitudeGuess, offsetGuess, pitchMm, fit, parameters},
    data = Transpose[{xs, row}];
    offsetGuess = Median[row];
    amplitudeGuess = Max[Max[row] - offsetGuess, 1.];
    pitchMm = If[Length[xs] > 1, xs[[2]] - xs[[1]], sigmaCloud];
    fit = Quiet[Check[NonlinearModelFit[data, fitAmplitude Exp[-(fitVariable - fitCentre)^2/(2 fitWidth^2)] + fitOffset,
        {{fitAmplitude, amplitudeGuess}, {fitCentre, 0.}, {fitWidth, sigmaCloud}, {fitOffset, offsetGuess}}, fitVariable, MaxIterations -> 300], $Failed]];
    If[fit === $Failed, Return[$Failed]];
    parameters = fit["BestFitParameters"][[All, 2]];   (* {amplitude, centre, width, offset} *)
    (* reject fits that locked onto a single noise spike or diverged: the width must lie between one superpixel and the field of view *)
    If[parameters[[1]] > 0 && pitchMm <= Abs[parameters[[3]]] <= (Last[xs] - First[xs]) && First[xs] <= parameters[[2]] <= Last[xs], fit, $Failed]];
crossSectionPlot[frame_List, meanFrame_List, xs_List, ys_List, fit_] :=
  Module[{rowIndex, row, expectedRow, graphicsList, t},
    rowIndex = First[Ordering[Abs[ys], 1]];
    row = frame[[rowIndex]];
    expectedRow = meanFrame[[rowIndex]];
    graphicsList = {ListPlot[Transpose[{xs, row}], PlotStyle -> Directive[Gray, PointSize[0.006]], Joined -> False],
                    ListLinePlot[Transpose[{xs, expectedRow}], PlotStyle -> Directive[RGBColor[0.1, 0.5, 0.9], AbsoluteThickness[1.5]]]};
    If[fit =!= $Failed, AppendTo[graphicsList, Quiet[Plot[fit[t], {t, First[xs], Last[xs]}, PlotStyle -> Directive[Red, AbsoluteThickness[2]], PlotPoints -> 100], General::munfl]]];
    Show[graphicsList,
         Frame -> True, Axes -> False, FrameLabel -> {"x (mm)", "e-/superpixel"},
         PlotLabel -> Style["Centre row: data (grey), expected mean (blue), Gaussian fit (red)", 10],
         PlotRange -> All, ImageSize -> 420, AspectRatio -> 0.7]];
metricsCard[cam_Association, name_String, mode_String, binFactor_, nMol_, rate_, tExp_, multiplier_, readNoiseRMS_, frame_List, meanFrame_List, xs_List, ys_List, fit_] :=
  Module[{nSig, signalElectrons, nSuper, bgElectrons, dkElectrons, readEff, breakdown, snr, estimators, minThree, minTen, peakElectrons, fullWellFraction, nMolSat,
          radiusFrame, apertureMask, outerMask, empiricalSignal, empiricalBackground, fitParameters, fitErrors, fitSigma, fitSigmaError, fitFWHM, rows},
    nSig = signalPhotonsFromMolecules[nMol, tExp, rate];
    signalElectrons = nSig cam["QE"];
    nSuper = superpixelCount[cam, binFactor];
    bgElectrons = backgroundElectrons[cam, binFactor, tExp, multiplier];
    dkElectrons = darkElectrons[cam, binFactor, tExp];
    readEff = effectiveReadNoise[cam, binFactor, readNoiseRMS];
    breakdown = noiseFloorBreakdown[cam, binFactor, tExp, multiplier, readNoiseRMS];
    snr = integratedSNR[cam, binFactor, nSig, tExp, multiplier, readNoiseRMS];
    estimators = snrEstimators[cam, binFactor, nMol, tExp, multiplier, readNoiseRMS, rate];
    minThree = minimumSignalPhotons[cam, binFactor, 3, tExp, multiplier, readNoiseRMS];
    minTen = minimumSignalPhotons[cam, binFactor, 10, tExp, multiplier, readNoiseRMS];
    peakElectrons = Max[meanFrame];
    fullWellFraction = peakElectrons/cam["FullWell"];
    nMolSat = saturationMolecules[cam, binFactor, tExp, rate];
    radiusFrame = Sqrt[Outer[Plus, ys^2, xs^2]];
    apertureMask = UnitStep[apertureRadius - radiusFrame];
    outerMask = UnitStep[radiusFrame - 1.4 apertureRadius];
    empiricalBackground = Total[frame outerMask, 2]/Max[Total[outerMask, 2], 1];
    empiricalSignal = Total[frame apertureMask, 2] - empiricalBackground Total[apertureMask, 2];
    {fitSigma, fitSigmaError, fitFWHM} = If[fit === $Failed, {Style["no reliable fit (peak not resolved above the noise)", Italic, Gray], "", ""},
       fitParameters = fit["BestFitParameters"][[All, 2]];
       fitErrors = Quiet[Check[fit["ParameterErrors"], ConstantArray[Indeterminate, 4]]];
       {Abs[fitParameters[[3]]], fitErrors[[3]], 2 Sqrt[2 Log[2]] Abs[fitParameters[[3]]]}];
    rows = {
      {"Camera", cameraShortName[name] <> "  (" <> cam["Type"] <> ", " <> cam["BinningType"] <> " binning)"},
      {"Readout mode / read noise", mode <> ": " <> ToString[readNoiseRMS] <> " e- rms per pixel"},
      {"Binning / superpixel", binLabel[binFactor] <> "  ->  " <> ToString[NumberForm[sensorSuperpixelPitch[cam, binFactor], {5, 2}]] <> " um on chip = " <>
            ToString[NumberForm[superpixelPitch[cam, binFactor], {5, 2}]] <> " um at the MOT (M = " <> ToString[NumberForm[magnification, {4, 3}]] <> "), frame " <> ToString[Dimensions[frame]] <> " superpixels"},
      {"Superpixels in aperture", nSuper},
      {"Exposure / background multiplier", ToString[NumberForm[1000 tExp, {5, 1}]] <> " ms / " <> ToString[NumberForm[multiplier, {4, 1}]] <> "\[Times] baseline"},
      {"Molecules / scattering rate", Row[{roundSignificant[nMol, 4], "  /  ", ScientificForm[rate, 2], " photons/s per molecule"}]},
      {"Signal photons in aperture (per molecule / total)", Row[{roundSignificant[photonsPerMolecule[tExp, rate], 4], "  /  ", roundSignificant[nSig, 4]}]},
      {"Signal electrons in aperture (model)", roundSignificant[signalElectrons, 4]},
      {"Signal electrons in aperture (from frame)", roundSignificant[empiricalSignal, 4]},
      {"Background e-/superpixel", roundSignificant[bgElectrons, 4]},
      {"Dark e-/superpixel", roundSignificant[dkElectrons, 4]},
      {"Effective read noise (e- rms/superpixel)", roundSignificant[readEff, 4]},
      {"Noise-floor variance (e-^2)", Row[{roundSignificant[breakdown["FloorVariance"], 4], "   dominant: ", breakdown["DominantFloorTerm"]}]},
      {"Integrated aperture SNR (top-hat, analytic)", Style[roundSignificant[snr, 4], Bold, If[snr >= 10, Darker[Green], If[snr >= 3, Darker[Orange], Red]]]},
      {"SNR by definition (peak pixel / top-hat / matched filter)",
        Row[{roundSignificant[estimators["PeakPixelSNR"], 3], "  /  ", roundSignificant[estimators["TopHatSNR"], 4], "  /  ",
             Style[roundSignificant[estimators["MatchedFilterSNR"], 4], Bold], "   (optimal weighting gains ", roundSignificant[estimators["MatchedGain"], 3], "\[Times])"}]},
      {"Precision on molecule number / cloud width", Row[{roundSignificant[100 estimators["NumberPrecision"], 3], "%  /  ", roundSignificant[100 estimators["WidthPrecision"], 3], "%"}]},
      {"Molecules for SNR = 3 / SNR = 10", Row[{roundSignificant[moleculesFromSignalPhotons[minThree, tExp, rate], 4], "  /  ", roundSignificant[moleculesFromSignalPhotons[minTen, tExp, rate], 4],
            "   (", roundSignificant[minThree, 4], " / ", roundSignificant[minTen, 4], " photons)"}]},
      {"Peak superpixel mean (e-) / full well", Row[{roundSignificant[peakElectrons, 4], "  (", roundSignificant[100 fullWellFraction, 3], "% of ", cam["FullWell"], " e-)",
            If[fullWellFraction >= 1, Style["  SATURATED", Bold, Red], ""]}]},
      {"Molecules at full-well saturation (this exposure)", Style[roundSignificant[nMolSat, 3], If[nMol >= nMolSat, Red, Black]]},
      {"Fitted sigma (mm) / true sigma", Row[{fitSigma, If[fitSigmaError === "", "", Row[{" \[PlusMinus] ", fitSigmaError}]], "   /   ", sigmaCloud}]},
      {"Fitted FWHM (mm) / true FWHM", Row[{fitFWHM, "   /   ", roundSignificant[cloudFWHM, 4]}]}};
    Panel[Grid[rows, Alignment -> {{Left, Left}}, Dividers -> {False, {False, {LightGray}, False}}, Spacings -> {2, 0.6}, BaseStyle -> {FontSize -> 11}],
      Style["Performance metrics", Bold, 12]]];
dashboardView[name_String, modeRequested_, binRequested_, nMol_, rate_, tExp_, multiplier_, seed_] :=
  Module[{cam, bins, binFactor, modes, mode, readNoiseRMS, nSig, frame, meanFrame, xs, ys, fit, rowIndex, image, cut, card},
    cam = cameraDatabase[name];
    nSig = signalPhotonsFromMolecules[nMol, tExp, rate];
    bins = cam["SupportedBins"];
    binFactor = If[MemberQ[bins, binRequested], binRequested, First[bins]];
    modes = cam["ReadNoiseModes"];
    mode = If[KeyExistsQ[modes, modeRequested], modeRequested, First[Keys[modes]]];
    readNoiseRMS = modes[mode];
    SeedRandom[Round[seed]];
    {xs, ys} = frameCoordinates[cam, binFactor];
    meanFrame = signalMeanFrame[cam, binFactor, nSig] + backgroundElectrons[cam, binFactor, tExp, multiplier] + darkElectrons[cam, binFactor, tExp];
    frame = syntheticFrame[cam, binFactor, nSig, tExp, multiplier, readNoiseRMS];
    rowIndex = First[Ordering[Abs[ys], 1]];
    fit = fitCrossSection[xs, frame[[rowIndex]]];
    image = frameImage[frame, xs, ys, cameraShortName[name] <> ", " <> binLabel[binFactor] <> ", " <> ToString[NumberForm[1000 tExp, {5, 1}]] <> " ms, " <> ToString[roundSignificant[nMol, 3]] <> " molecules"];
    cut = crossSectionPlot[frame, meanFrame, xs, ys, fit];
    card = metricsCard[cam, name, mode, binFactor, nMol, rate, tExp, multiplier, readNoiseRMS, frame, meanFrame, xs, ys, fit];
    Column[{Grid[{{image, cut}}, Alignment -> Top, Spacings -> {2, 0}], card}, Spacings -> 1.5]];
motDashboard = Manipulate[
   dashboardView[selectedCamera, readoutMode, binFactor, 10.^logMolecules, 10.^6 scatteringRateMHz, tExpMs/1000., backgroundMultiplier, noiseSeed],
   {{selectedCamera, cameraNames[[1]], "Camera"}, cameraNames, ControlType -> PopupMenu},
   {{readoutMode, "Ultra Quiet", "Readout mode"}, Keys[cameraDatabase[selectedCamera]["ReadNoiseModes"]], ControlType -> PopupMenu},
   {{binFactor, 1, "Binning"}, (# -> binLabel[#] &) /@ cameraDatabase[selectedCamera]["SupportedBins"], ControlType -> SetterBar},
   Delimiter,
   {{logMolecules, 2., "log10 molecules in the MOT"}, 1., 6., 0.05, Appearance -> "Labeled"},
   {{scatteringRateMHz, 1., "Scattering rate (10^6 photons/s per molecule)"}, 0.1, 5., 0.1, Appearance -> "Labeled"},
   {{tExpMs, 20, "Exposure tExp (ms)"}, 1, 200, 1, Appearance -> "Labeled"},
   {{backgroundMultiplier, backgroundMultiplierDefault, "Background multiplier"}, 0.5, 10., 0.1, Appearance -> "Labeled"},
   {{noiseSeed, 1, "Noise realisation seed"}, 1, 50, 1, Appearance -> "Labeled"},
   ControlPlacement -> Top, ContinuousAction -> False, SynchronousUpdating -> False, SaveDefinitions -> True,
   TrackedSymbols :> {selectedCamera, readoutMode, binFactor, logMolecules, scatteringRateMHz, tExpMs, backgroundMultiplier, noiseSeed},
   FrameMargins -> 8]


(* ::Section:: *)
(*7. Evaluation Dataset: Minimum Number of Molecules for SNR = 3 and SNR = 10*)


(* ::Text:: *)
(*Closed-form minimum number of molecules in the MOT (20 ms exposure, 4x baseline scatter, 1e6 photons/s per molecule) for every camera, quoted readout mode and benchmark binning, together with the equivalent signal photons in the aperture and the molecule number at which the peak superpixel saturates.*)


(* ::Input:: *)
thresholdTable = Dataset[Flatten[Table[
    Module[{cam = cameraDatabase[name], readNoiseRMS = cameraDatabase[name]["ReadNoiseModes"][mode], minThree, minTen},
      minThree = minimumSignalPhotons[cam, binFactor, 3, exposureTimeNominal, backgroundMultiplierDefault, readNoiseRMS];
      minTen = minimumSignalPhotons[cam, binFactor, 10, exposureTimeNominal, backgroundMultiplierDefault, readNoiseRMS];
      <|"Camera" -> cameraShortName[name], "Readout mode" -> mode, "Read noise (e-)" -> readNoiseRMS,
        "Bin" -> binLabel[binFactor], "nSuper" -> superpixelCount[cam, binFactor],
        "Eff. read noise (e-)" -> effectiveReadNoise[cam, binFactor, readNoiseRMS],
        "Bg e-/superpixel" -> roundSignificant[backgroundElectrons[cam, binFactor, exposureTimeNominal, backgroundMultiplierDefault], 4],
        "nMol (SNR=3)" -> roundSignificant[moleculesFromSignalPhotons[minThree, exposureTimeNominal], 4],
        "nMol (SNR=10)" -> roundSignificant[moleculesFromSignalPhotons[minTen, exposureTimeNominal], 4],
        "nSig (SNR=3)" -> roundSignificant[minThree, 4], "nSig (SNR=10)" -> roundSignificant[minTen, 4],
        "Saturation nMol" -> roundSignificant[saturationMolecules[cam, binFactor], 3],
        "Primary mode" -> (readNoiseRMS == cam["ReadNoiseRMS"])|>],
    {name, cameraNames}, {mode, Keys[cameraDatabase[name]["ReadNoiseModes"]]}, {binFactor, benchmarkBins}], 2]]


(* ::Text:: *)
(*Ranking at 1x1 binning in the primary readout mode (lower is better).*)


(* ::Input:: *)
thresholdTable[Select[#Bin == "1\[Times]1" && #["Primary mode"] &]][SortBy[#["nMol (SNR=3)"] &]][All, {"Camera", "Readout mode", "nMol (SNR=3)", "nMol (SNR=10)", "Saturation nMol", "nSig (SNR=3)"}]


(* ::Section:: *)
(*8. Verification Tests*)


(* ::Text:: *)
(*Automated assertions run by the generator script: collection geometry, read-noise binning rules, synthetic frame dimensions against the ROI calculation, consistency of the closed-form threshold with the SNR function, physical monotonicity, and absence of global-scope leakage from the dashboard code.*)


(* ::Input:: *)
verificationTests = {
   VerificationTest[Abs[geometricEfficiencyUnvignetted - 0.0559] < 0.001, True, TestID -> "GeometricEfficiencyUnvignetted"],
   VerificationTest[Abs[solidAngleCollection[workingDistance, collectorRadius] - 0.7025] < 0.001, True, TestID -> "SolidAngleUnvignetted"],
   VerificationTest[Abs[geometricEfficiency - 0.0369] < 0.001, True, TestID -> "GeometricEfficiencyVignetted"],
   VerificationTest[limitingElement, "Comar 29 AF 40", TestID -> "LimitingApertureIsImagingLens"],
   VerificationTest[Abs[magnification - 28.6/59.8] < 10^-9 && Abs[magnification - 0.4783] < 0.001, True, TestID -> "Magnification"],
   VerificationTest[Abs[superpixelPitch[cameraDatabase["Thorlabs Zelux CS165MU (Compact CMOS)"], 1] - 3.45/magnification] < 10^-9, True, TestID -> "ObjectSpacePitch"],
   VerificationTest[Module[{cam = cameraDatabase["Teledyne Photometrics Kinetix (sCMOS Benchmark)"]},
       Abs[backgroundElectrons[cam, 1, 0.02, 4.] - backgroundFluxDensity[4.] 6.5^2 0.02 cam["QE"]] < 10^-9], True, TestID -> "BackgroundUsesSensorAreaNotObjectArea"],
   VerificationTest[Abs[backgroundFluxDensity[4.0] - 3.6923] < 0.001, True, TestID -> "BackgroundFluxDensity"],
   VerificationTest[Abs[apertureAreaMicron - 1.7766*^7] < 20000, True, TestID -> "ApertureArea"],
   VerificationTest[Abs[cloudFWHM - 2.8] < 10^-9 && Abs[sigmaCloud - 1.18903] < 10^-4, True, TestID -> "CloudSize"],
   VerificationTest[Abs[apertureSignalFraction - (1 - Exp[-2.])] < 10^-9, True, TestID -> "ApertureSignalFraction"],
   VerificationTest[effectiveReadNoise[cameraDatabase["Hamamatsu ORCA-R2 (Cooled CCD C10600-10B)"], 4], 6.0, TestID -> "ReadNoiseHardwareBinningORCAR2"],
   VerificationTest[effectiveReadNoise[cameraDatabase["Thorlabs Zelux CS165MU (Compact CMOS)"], 4], 16.0, TestID -> "ReadNoiseSoftwareBinningZelux"],
   VerificationTest[effectiveReadNoise[cameraDatabase["Hamamatsu ORCA-Quest 2 (qCMOS C15550-22UP)"], 2, 0.43], 0.86, TestID -> "ReadNoiseSoftwareBinningQuestStandardMode"],
   VerificationTest[effectiveReadNoise[cameraDatabase["Andor CB2 High Speed 7.1F (sCMOS)"], #] & /@ {1, 2, 4, 8}, {1.4, 1.4, 2.8, 5.6}, TestID -> "ReadNoiseHybridBinningCB2"],
   VerificationTest[Abs[cameraDatabase["Andor CB2 High Res (24.5 MP BSI sCMOS)"]["QE"] - 0.62] < 0.03 && Abs[cameraDatabase["Andor CB2 High Speed 7.1F (sCMOS)"]["QE"] - 0.68] < 0.03, True, TestID -> "CB2QEAt606nm"],
   VerificationTest[And @@ (KeyExistsQ[#, "HardwareBinLimit"] & /@ Values[cameraDatabase]), True, TestID -> "AllCamerasHaveHardwareBinLimit"],
   VerificationTest[superpixelCount[cameraDatabase["Hamamatsu ORCA-Quest 2 (qCMOS C15550-22UP)"], 1], Round[apertureAreaMicron/(4.6/magnification)^2], TestID -> "SuperpixelCountQuest"],
   VerificationTest[Dimensions[syntheticFrame[cameraDatabase[#], 1, 10.^4, 0.02, 4., 1.]] == frameDimensions[cameraDatabase[#], 1] & /@ cameraNames, ConstantArray[True, Length[cameraNames]], TestID -> "FrameDimensionsBin1"],
   VerificationTest[Dimensions[syntheticFrame[cameraDatabase[#], 4, 10.^4, 0.02, 4., 1.]] == frameDimensions[cameraDatabase[#], 4] & /@ cameraNames, ConstantArray[True, Length[cameraNames]], TestID -> "FrameDimensionsBin4"],
   VerificationTest[frameDimensions[cameraDatabase["Thorlabs Zelux CS165MU (Compact CMOS)"], 1],
       {Min[Floor[1000 syntheticFrameFOV magnification/3.45], 1080], Min[Floor[1000 syntheticFrameFOV magnification/3.45], 1440]}, TestID -> "FrameDimensionsZeluxROI"],
   VerificationTest[frameDimensions[cameraDatabase["Andor iXon Ultra 888 (EMCCD Benchmark)"], 8],
       {Min[Floor[1000 syntheticFrameFOV magnification/104.], 128], Min[Floor[1000 syntheticFrameFOV magnification/104.], 128]}, TestID -> "FrameDimensionsIXonBin8"],
   VerificationTest[Module[{cam = cameraDatabase["Thorlabs Zelux CS165MU (Compact CMOS)"], halfWidth},
       (* the synthetic frame must be wide enough to hold the off-target reference annulus used by the metrics card *)
       halfWidth = Last[First[frameCoordinates[cam, 1]]];
       halfWidth > 1.4 apertureRadius], True, TestID -> "FrameHoldsReferenceAnnulus"],
   VerificationTest[Module[{cam = cameraDatabase["Andor iXon Ultra 888 (EMCCD Benchmark)"], frame, xs, ys, mask},
       {xs, ys} = frameCoordinates[cam, 1]; frame = signalMeanFrame[cam, 1, 10.^6];
       mask = UnitStep[apertureRadius - Sqrt[Outer[Plus, ys^2, xs^2]]];
       Abs[Total[frame mask, 2]/(10.^6 cam["QE"]) - 1] < 0.02], True, TestID -> "SyntheticFrameApertureSignal"],
   VerificationTest[And @@ Flatten[Table[Abs[integratedSNR[cameraDatabase[name], binFactor, minimumSignalPhotons[cameraDatabase[name], binFactor, target]] - target] < 10^-8,
       {name, cameraNames}, {binFactor, benchmarkBins}, {target, {3, 10}}]], True, TestID -> "ThresholdRoundTrip"],
   VerificationTest[And @@ Table[Less @@ (integratedSNR[cameraDatabase[name], 1, #] & /@ {10.^2, 10.^4, 10.^6}), {name, cameraNames}], True, TestID -> "SNRMonotonicInSignal"],
   VerificationTest[And @@ Table[integratedSNR[cameraDatabase[name], 1, 10.^4, 0.02, 10., cameraDatabase[name]["ReadNoiseRMS"]] < integratedSNR[cameraDatabase[name], 1, 10.^4, 0.02, 0.5, cameraDatabase[name]["ReadNoiseRMS"]], {name, cameraNames}], True, TestID -> "SNRDecreasesWithScatter"],
   (* the asymptotic slopes are defined relative to each camera's own shot-noise crossover, so they stay meaningful if the cloud size, exposure or background change *)
   VerificationTest[Module[{cam = cameraDatabase["Thorlabs Zelux CS165MU (Compact CMOS)"], crossover},
       crossover = noiseFloorBreakdown[cam, 1, exposureTimeNominal, backgroundMultiplierDefault, cam["ReadNoiseRMS"]]["ShotNoiseCrossoverPhotons"];
       Abs[localLogSlope[cam, 1, crossover/100] - 1] < 0.01], True, TestID -> "SlopeReadNoiseLimited"],
   VerificationTest[Module[{cam = cameraDatabase["Andor iXon Ultra 888 (EMCCD Benchmark)"], crossover},
       crossover = noiseFloorBreakdown[cam, 1, exposureTimeNominal, backgroundMultiplierDefault, cam["ReadNoiseRMS"]]["ShotNoiseCrossoverPhotons"];
       Abs[localLogSlope[cam, 1, 100 crossover] - 0.5] < 0.01], True, TestID -> "SlopeShotNoiseLimited"],
   VerificationTest[Module[{before, after},
       before = Names["Global`*"];
       dashboardView[cameraNames[[3]], "Standard", 2, 100., 1.*^6, 0.02, 4., 7];
       after = Names["Global`*"];
       Complement[after, before]], {}, TestID -> "NoGlobalScopeLeakage"],
   VerificationTest[ValueQ /@ {fitAmplitude, fitCentre, fitWidth, fitOffset, fitVariable}, {False, False, False, False, False}, TestID -> "FitParametersUnassigned"],
   VerificationTest[And @@ Flatten[Table[Module[{est = snrEstimators[cameraDatabase[name], binFactor, nMol]},
       est["MatchedFilterSNR"] >= est["TopHatSNR"] - 10^-6 && est["PeakPixelSNR"] <= est["TopHatSNR"]],
       {name, cameraNames}, {binFactor, benchmarkBins}, {nMol, {10., 10.^4}}]], True, TestID -> "MatchedFilterBeatsTopHat"],
   VerificationTest[Module[{cam = cameraDatabase["Hamamatsu ORCA-Quest 2 (qCMOS C15550-22UP)"], est, shotLimit},
       (* with no read noise, no background and no dark current the matched filter must reach the pure shot-noise limit Sqrt[S] *)
       est = snrEstimators[cam, 1, 10.^4, exposureTimeNominal, 0., 0., scatteringRate];
       shotLimit = Sqrt[signalPhotonsFromMolecules[10.^4, exposureTimeNominal] cam["QE"]/apertureSignalFraction];
       Abs[est["MatchedFilterSNR"]/shotLimit - 1] < 0.01], True, TestID -> "MatchedFilterReachesShotNoiseLimit"],
   VerificationTest[And @@ Flatten[Table[Abs[peakPixelSNR[cameraDatabase[name], binFactor, 10.^4]/snrEstimators[cameraDatabase[name], binFactor, 10.^4]["PeakPixelSNR"] - 1] < 0.02,
       {name, cameraNames}, {binFactor, benchmarkBins}]], True, TestID -> "PeakPixelSNRClosedFormMatchesGrid"],
   VerificationTest[And @@ Flatten[Table[Abs[peakPixelSNR[cameraDatabase[name], binFactor, visibilityMolecules[cameraDatabase[name], binFactor, target]] - target] < 10^-8,
       {name, cameraNames}, {binFactor, benchmarkBins}, {target, {3, 5}}]], True, TestID -> "VisibilityThresholdRoundTrip"],
   VerificationTest[Module[{cam = cameraDatabase["Thorlabs Zelux CS165MU (Compact CMOS)"]},
       (* software binning must leave the integrated SNR alone while the peak-pixel SNR grows in proportion to b *)
       And[Abs[integratedSNRMolecules[cam, 4, 10.^4]/integratedSNRMolecules[cam, 1, 10.^4] - 1] < 0.01,
           Abs[peakPixelSNR[cam, 4, 10.^4]/peakPixelSNR[cam, 1, 10.^4] - 4] < 0.05]], True, TestID -> "BinningTradesPrecisionForVisibility"],
   VerificationTest[Module[{cam = cameraDatabase["Thorlabs Zelux CS165MU (Compact CMOS)"], topHat, analytic},
       (* the top-hat estimator in snrEstimators must reproduce the closed-form integratedSNR of Section 1.4 *)
       topHat = snrEstimators[cam, 1, 10.^4]["TopHatSNR"];
       analytic = integratedSNRMolecules[cam, 1, 10.^4];
       Abs[topHat/analytic - 1] < 0.02], True, TestID -> "TopHatMatchesClosedForm"],
   VerificationTest[Abs[photonsPerMolecule[0.02, 1.*^6] - 1.*^6 0.02 geometricEfficiency 0.95 (1 - Exp[-2])] < 10^-6, True, TestID -> "PhotonsPerMoleculeConversion"],
   VerificationTest[And @@ Flatten[Table[Abs[integratedSNRMolecules[cameraDatabase[name], binFactor, minimumMolecules[cameraDatabase[name], binFactor, target]] - target] < 10^-8,
       {name, cameraNames}, {binFactor, benchmarkBins}, {target, {3, 10}}]], True, TestID -> "MoleculeThresholdRoundTrip"],
   VerificationTest[Module[{cam = cameraDatabase["Hamamatsu ORCA-Quest 2 (qCMOS C15550-22UP)"], nMolSat},
       nMolSat = saturationMolecules[cam, 1];
       Abs[Max[signalMeanFrame[cam, 1, signalPhotonsFromMolecules[nMolSat, exposureTimeNominal]]]/cam["FullWell"] - 1] < 0.02], True, TestID -> "SaturationMoleculesMatchesFrame"],
   VerificationTest[StringFreeQ[#, "_"] & /@ Select[Names["Global`*"], StringFreeQ[#, "$"] &], ConstantArray[True, Length[Select[Names["Global`*"], StringFreeQ[#, "$"] &]]], TestID -> "NoUnderscoresInSymbolNames"]};
verificationReport = TestReport[verificationTests]


(* ::Input:: *)
Dataset[Table[<|"TestID" -> result["TestID"], "Outcome" -> result["Outcome"]|>, {result, Values[verificationReport["TestResults"]]}]]
