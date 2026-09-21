#!/usr/bin/env wolframscript
(* ::Package:: *)

(* GenerateCaFComparisonNotebook.wls
   Builds and evaluates CaF_MOT_Camera_Comparison.nb from CaFCameraNoise.wl.

   Pipeline
     1. Parse CaFCameraNoise.wl (Wolfram package-cell format) into Title/Section/Text/Input chunks.
     2. Evaluate every Input chunk in this kernel, in order, capturing the result.
     3. Run the VerificationTest suite defined in the source (Section 6) and print a report.
     4. Reparse the Input source into box structures with the (headless) front end, append Output cells
        built from the kernel results, and save the notebook.
     5. Rasterise the key figures to verification/*.png and export the threshold table to CSV.
   Exit code is 1 if any verification test fails or any Input chunk raises a message.

   NotebookEvaluate is deliberately not used: under wolframscript it re-enters the calling kernel and deadlocks,
   so evaluation happens in the kernel and results are typeset with ToBoxes instead. *)

scriptDirectory = If[$InputFileName === "", Directory[], DirectoryName[$InputFileName]];
sourceFile = FileNameJoin[{scriptDirectory, "CaFCameraNoise.wl"}];
outputNotebook = FileNameJoin[{scriptDirectory, "CaF_MOT_Camera_Comparison.nb"}];
outputCSV = FileNameJoin[{scriptDirectory, "CaF_MOT_Camera_Thresholds.csv"}];
verificationDirectory = FileNameJoin[{scriptDirectory, "verification"}];
generatorStartTime = AbsoluteTime[];

logLine[args___] := Print["[", ToString[NumberForm[AbsoluteTime[] - generatorStartTime, {6, 1}]], " s] ", args];

logLine["Wolfram kernel ", $Version];
logLine["Source: ", sourceFile];
If[! FileExistsQ[sourceFile], Print["ERROR: source file not found"]; Exit[2]];

(* ---------- 1. Parse the package-format source into cell chunks ---------- *)

sourceText = Import[sourceFile, "Text"];
sourceLines = StringSplit[sourceText, "\n", All];
cellMarkerPattern = StartOfString ~~ "(* ::" ~~ style : (WordCharacter ..) ~~ ":: *)" ~~ EndOfString;

parseSourceChunks[lines_List] :=
  Module[{chunks = {}, currentStyle = None, buffer = {}, flush},
    flush[] := If[currentStyle =!= None && currentStyle =!= "Package",
       AppendTo[chunks, <|"Style" -> currentStyle, "Lines" -> buffer|>]];
    Do[
      If[StringMatchQ[line, cellMarkerPattern],
        flush[];
        currentStyle = First[StringCases[line, cellMarkerPattern :> style]];
        buffer = {},
        AppendTo[buffer, line]],
      {line, lines}];
    flush[];
    chunks];

stripCommentDelimiters[line_String] :=
  StringReplace[line, {StartOfString ~~ "(*" -> "", "*)" ~~ EndOfString -> ""}];

chunkContent[chunk_Association] :=
  Module[{lines = chunk["Lines"], style = chunk["Style"]},
    lines = Reverse[Drop[Reverse[lines], LengthWhile[Reverse[lines], StringTrim[#] === "" &]]];
    lines = Drop[lines, LengthWhile[lines, StringTrim[#] === "" &]];
    If[style === "Input",
      StringRiffle[lines, "\n"],
      StringRiffle[stripCommentDelimiters /@ Select[lines, StringTrim[#] =!= "" &], "\n"]]];

sourceChunks = parseSourceChunks[sourceLines];
sourceChunks = Map[Append[#, "Content" -> chunkContent[#]] &, sourceChunks];
inputChunkCount = Count[sourceChunks, KeyValuePattern["Style" -> "Input"]];
logLine["Parsed ", Length[sourceChunks], " cells (", inputChunkCount, " Input cells)"];

(* ---------- 2. Evaluate every Input chunk in the kernel ---------- *)

evaluationMessagesSeen = 0;
evaluateChunk[code_String, index_Integer] :=
  Module[{result, timing, messageCountBefore = Length[$MessageList]},
    {timing, result} = AbsoluteTiming[ToExpression[code, InputForm]];
    logLine["Evaluated Input ", index, " in ", ToString[NumberForm[timing, {5, 2}]], " s -> ", Head[result]];
    result];

inputCounter = 0;
sourceChunks = Map[
   If[#["Style"] === "Input",
     inputCounter++;
     Append[#, {"Index" -> inputCounter, "Result" -> evaluateChunk[#["Content"], inputCounter]}],
     #] &,
   sourceChunks];

(* ---------- 3. Verification report ---------- *)

If[! MatchQ[verificationReport, _TestReportObject],
  Print["ERROR: verificationReport was not produced by the source file"]; Exit[1]];

testOutcomes = Table[{result["TestID"], result["Outcome"]}, {result, Values[verificationReport["TestResults"]]}];
Print[""];
Print["=== Verification tests ==="];
Scan[Print[StringPadRight[#[[1]], 44], If[#[[2]] === "Success", "PASS", "FAIL  <-- " <> ToString[#[[2]]]]] &, testOutcomes];
allTestsPassed = TrueQ[verificationReport["AllTestsSucceeded"]];
Print["All tests passed: ", allTestsPassed, "  (", verificationReport["TestsSucceededCount"], "/", Length[testOutcomes], ")"];
Print[""];

(* Independent cross-checks of the headline assertions, evaluated here rather than inside the notebook *)
crossChecks = {
   {"etaGeom unvignetted (2 in) = 0.0559 +- 0.001", Abs[geometricEfficiencyUnvignetted - 0.0559] < 0.001},
   {"etaGeom vignetted by Comar 40 mm = 0.0369", Abs[geometricEfficiency - 0.0369] < 0.001},
   {"magnification f2/f1 = 0.478", Abs[magnification - 0.4783] < 0.001},
   {"matched filter >= top-hat for every camera/bin",
     And @@ Flatten[Table[Module[{est = snrEstimators[cameraDatabase[name], binFactor, 10.^4]}, est["MatchedFilterSNR"] >= est["TopHatSNR"] - 10^-6],
       {name, cameraNames}, {binFactor, benchmarkBins}]]},
   {"ORCA-R2 sigmaReadEff(4x4) = 6.0", effectiveReadNoise[cameraDatabase["Hamamatsu ORCA-R2 (Cooled CCD C10600-10B)"], 4] == 6.0},
   {"Zelux sigmaReadEff(4x4) = 16.0", effectiveReadNoise[cameraDatabase["Thorlabs Zelux CS165MU (Compact CMOS)"], 4] == 16.0},
   {"Frame dims == ROI calc (all cameras, bins 1,2,4,8)",
     And @@ Flatten[Table[Dimensions[syntheticFrame[cameraDatabase[name], binFactor, 10.^4, 0.02, 4., 1.]] === frameDimensions[cameraDatabase[name], binFactor],
       {name, cameraNames}, {binFactor, {1, 2, 4, 8}}]]}};
Print["=== Generator cross-checks ==="];
Scan[Print[StringPadRight[#[[1]], 56], If[TrueQ[#[[2]]], "PASS", "FAIL"]] &, crossChecks];
allCrossChecksPassed = And @@ (TrueQ[#[[2]]] & /@ crossChecks);
Print[""];

(* ---------- 4. Build the notebook with the headless front end ---------- *)

outputBoxes[result_] :=
  Quiet[Check[ToBoxes[result, StandardForm], Check[ToBoxes[Short[result, 10], StandardForm], ToBoxes["<< output could not be typeset >>"]]]];

buildNotebookCells[chunks_List] :=
  Module[{cells = {}, style, content, inputBoxes, result},
    Do[
      style = chunk["Style"];
      content = chunk["Content"];
      Switch[style,
        "Input",
          inputBoxes = MathLink`CallFrontEnd[FrontEnd`ReparseBoxStructurePacket[content]];
          AppendTo[cells, Cell[BoxData[inputBoxes], "Input", CellLabel -> "In[" <> ToString[chunk["Index"]] <> "]:="]];
          result = chunk["Result"];
          If[result =!= Null,
            AppendTo[cells, Cell[BoxData[outputBoxes[result]], "Output", CellLabel -> "Out[" <> ToString[chunk["Index"]] <> "]="]]],
        _,
          AppendTo[cells, Cell[content, style]]],
      {chunk, chunks}];
    cells];

notebookOptions = {
   WindowTitle -> "CaF MOT Camera Comparison",
   StyleDefinitions -> "Default.nb",
   WindowSize -> {1400, 900},
   TaggingRules -> {"GeneratedBy" -> "GenerateCaFComparisonNotebook.wls", "GeneratedOn" -> DateString[], "KernelVersion" -> $Version}};

logLine["Building notebook cells via the front end ..."];
notebookExpression = UsingFrontEnd[Notebook[buildNotebookCells[sourceChunks], Sequence @@ notebookOptions]];
logLine["Notebook expression built: ", Length[First[notebookExpression]], " cells, ", ByteCount[notebookExpression], " bytes"];

UsingFrontEnd[
  Module[{notebookObject},
    notebookObject = NotebookPut[notebookExpression];
    NotebookSave[notebookObject, outputNotebook];
    NotebookClose[notebookObject]]];
logLine["Saved ", outputNotebook, " (", FileByteCount[outputNotebook], " bytes)"];

(* ---------- 5. Verification images and CSV export ---------- *)

If[! DirectoryQ[verificationDirectory], CreateDirectory[verificationDirectory]];
rasterTargets = {
   {"benchmark_snr_loglog.png", benchmarkFigure},
   {"benchmark_per_camera.png", benchmarkGrid},
   {"regime_table.png", regimeDataset},
   {"resolution_vs_binning.png", resolutionPlots},
   {"estimator_table.png", estimatorDataset},
   {"estimator_comparison.png", estimatorComparisonPlot},
   {"visibility_per_camera.png", visibilityGrid},
   {"visibility_threshold.png", visibilityThresholdPlot},
   {"visibility_table.png", visibilityDataset},
   {"dashboard_initial_state.png", motDashboard},
   {"threshold_table.png", thresholdTable}};
UsingFrontEnd[
  Do[
    Module[{path = FileNameJoin[{verificationDirectory, target[[1]]}]},
      Export[path, Rasterize[Style[target[[2]], LightDark -> "Light"], "Image", ImageResolution -> 110, Background -> White]];   (* force light appearance: the headless front end follows the system dark mode *)
      logLine["Rasterised ", target[[1]]]],
    {target, rasterTargets}]];

Export[outputCSV, thresholdTable];
logLine["Exported ", outputCSV];

Print[""];
Print["Notebook:      ", outputNotebook];
Print["Threshold CSV: ", outputCSV];
Print["Verification:  ", verificationDirectory];
If[! (allTestsPassed && allCrossChecksPassed), Print["*** VERIFICATION FAILED ***"]; Exit[1]];
Print["*** ALL CHECKS PASSED ***"];
