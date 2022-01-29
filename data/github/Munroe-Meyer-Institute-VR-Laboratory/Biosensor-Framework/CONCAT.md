[![status](https://joss.theoj.org/papers/ce2098ed62af72faa9a4817dabdc34e3/status.svg)](https://joss.theoj.org/papers/ce2098ed62af72faa9a4817dabdc34e3) [![DOI](https://zenodo.org/badge/354086592.svg)](https://zenodo.org/badge/latestdoi/354086592)

## Biosensor Framework: A C# Library for Affective Computing
### Parties Involved
Institution: Munroe Meyer Institute in the University of Nebraska Medical Center<br>
Laboratory: Virtual Reality Laboratory<br>
Advisor: Dr. James Gehringer<br>
Developer: Walker Arce<br>

### Introduction
Developed in the Virtual Reality Laboratory at the Munroe Meyer Institute, this library provides an interface for interacting with biosensors and performing affective computing tasks on their output data streams.  Currently, it natively supports the [Empatica E4](https://www.empatica.com/research/e4/) and communicates through a TCP connection with their provided server.  

This library provides the following capabilities:
-	Connection to the [TCP client](https://github.com/empatica/ble-client-windows) that Empatica provides for connecting their devices to a user's PC
-	Parsing and packing of commands to communicate bidirectionally with the TCP client
-	Collection and delivery of biometric readings to an end user's software application
-	Extraction of literature supported feature vectors to allow characterization of the biological signal features
-	Microsoft.ML support for inferencing on the feature vectors natively in C# in both multi-class and binary classification tasks
-	Support for simulating a data stream for debugging and testing applications without an Empatica E4 device
-	Support for parsing open source Empatica E4 datasets and training machine learning models on them (i.e. [WESAD](https://archive.ics.uci.edu/ml/datasets/WESAD+%28Wearable+Stress+and+Affect+Detection%29))
-	Support for interfacing new sensors into the same pipeline
-	Contains models trained on WESAD dataset in the Microsoft.ML framework

### Motivation
For our lab to use biosensor data in our research projects, we needed an interface to pull in data from our Empatica E4, perform feature extraction, and return inferences from a trained model.  From our research, there was no existing all in one library to do this and major parts of this project did not exist.  Additionally, there was no Unity compatible way to process this data and perform actions from the output.  From a programming perspective, it is best if there was no string manipulation when interacting with the server, so there needed to be a true wrapper for the TCP client communications.  The development of the Biosensor Framework has allowed us to use these features reliably in all of our research projects.  There will be errors about the loaded dependencies, though the project should still run with those present.

### Dependencies

    .NETStandard 2.0
        LightGBM (>= 2.2.3)
        MathNet.Numerics (>= 4.15.0)
        Microsoft.ML (>= 1.5.5)
        Microsoft.ML.FastTree (>= 1.5.5)
        Microsoft.ML.LightGbm (>= 1.5.5)

### Installation
The library can be installed using the [NuGet package manager](https://www.nuget.org/packages/Biosensor-Framework/) and is listed publicly.  For Unity, downloading the package and importing it through the Unity UI is the best approach.  This includes all of the dependency libraries, a test script, a test scene, and an example output file.  This example performs the same tasks as Example2 and needs to be supplied the E4 Streaming Server exe path, an output path, and your API key for your device.

The following command can be used to install the package using the NuGet CLI:

`Install-Package Biosensor-Framework -Version 1.0.0`

### Unity Package
The included Unity package includes the required DLLs, a test script, and an example output file.  To install the package, clone this repository and open the target Unity project. 

Click on Assets -> Import Package -> Custom Package...

![image](https://user-images.githubusercontent.com/22334349/127775654-de8f4aeb-92c4-4d99-8c17-1b8f3b738e78.png)

Navigate to the [Unity Package](https://github.com/Munroe-Meyer-Institute-VR-Laboratory/Biosensor-Framework/tree/main/Unity%20Package) folder and click on the "BiosensorFramework_Unity.unitypackage" file.  When the file has loaded in, click "All" and "Import" to import the package into your project.

![image](https://user-images.githubusercontent.com/22334349/127775776-724fa58d-6365-4613-9268-a722b357a8b5.png)

### Testing Installation
Once the package has been installed, this repo can be cloned to get the test [datasets](https://github.com/Munroe-Meyer-Institute-VR-Laboratory/Biosensor-Framework/tree/main/Dataset) and the [trained models](https://github.com/Munroe-Meyer-Institute-VR-Laboratory/Biosensor-Framework/tree/main/Trained%20Models).  
Follow Example 1 to test the functionality of the Empatica E4 TCP server communication and window data gathering.  
Follow Example 2 to test the functionality of Microsoft.ML for affective computing tasks with a connected body-worn sensor.  
Follow Example 3 to demonstrate the functionality of Microsoft.ML for affective computing tasks using the WESAD dataset without a connected body-worn sensor. 

### Basic Usage
This repository has two example projects.  Example1 shows the basic usage to communicate with the server and pull data off of a biosensor.  Example2 expands Example1 by adding in Microsoft.ML inferencing.  Example2 will be used to explain the operation of the library.

```
        using System;
        using System.Diagnostics;
        using System.IO;
        using System.Threading;
        using Microsoft.ML;

        using MMIVR.BiosensorFramework.Biosensors.EmpaticaE4;
        using MMIVR.BiosensorFramework.MachineLearningUtilities;
        
        namespace Example2_ComputingWithMicrosoftML
        {
            class Program
            {
                public static ServerClient Device1;
                static ServerClient client;
                public static MLContext mlContext;
                public static ITransformer Model;
                public static DataViewSchema ModelSchema;
                // TODO: Fill these out with your own values
                public static string APIKey = "";
                public static string ServerPath = @"";
                public static string WesadDirectory = @"";
                public static string SessionOutputPath = @"";
                public static string ModelDir = @"";
```
For each project, a client connection is established using the ServerClient class.  This handles the TCP initialization, command communication, and error handling.  Each device will have its own ServerClient instance, so Device1 is represented by a device name and a TCP connection.  Additionally, the API key and the path to the server executable need to be defined.  The API key can be found by following Empatica's directions from the Installation section.  The path to the server is found in the installation directory for the server (i.e. @"{installation path}\EmpaticaBLEServer.exe").

Next the Microsoft.ML variables are created.
```
        static void Main(string[] args)
        {
            mlContext = new MLContext();
            Train.RunBenchmarks(WesadDirectory, out ITransformer RegModel, out ITransformer MultiModel, out Model, ModelDir);
```
In this example, the models are trained on the WESAD data using the RunBenchmarks function.  The best performing models for each class are output.  For this example, the binary classification model (BinModel) is used.
```
            Console.WriteLine("E4 Console Interface - Press ENTER to begin the client");
            Console.ReadKey();

            Console.WriteLine("Step 1 - Start Empatica server");
            Utilities.StartE4ServerGUI(ServerPath);
```
The E4 server is started through a Process call.  This will open the GUI for the server, if that is not desired, the command line variant of the command can be called.
```
            client = new ServerClient();
            Console.ReadKey();
            client.StartClient();
            Utilities.StartE4Server(ServerPath, APIKey);
```
The E4 server can also be started through the command line variant as shown.  
```
            Console.WriteLine("Step 2 - Getting connected E4 devices");
            Utilities.ListDiscoveredDevices(client);
            Console.ReadKey();

            Console.WriteLine("     Available Devices:");
            Utilities.AvailableDevices.ForEach(i => Console.WriteLine("     Device Name: {0}", i));
            Console.ReadKey();
```
Listing the connected devices will show all devices connected through the Bluetooth dongle.  These are managed internally.
```
            Device1 = new ServerClient();
            Device1.StartClient();
            Device1.DeviceName = Utilities.AvailableDevices[0];
            Utilities.ConnectDevice(Device1);
            Console.ReadKey();
```
To connect a device, it needs to be assigned one of the available devices.  Once that's done, the ConnectDevice function can be called and its TCP connection will be established.
```
            Utilities.SuspendStreaming(Device1);
```
Since there are configurations to be done on the device, the streaming needs to be suspended.
```
            Console.WriteLine("Step 3 - Adding biometric data streams");

            foreach (ServerClient.DeviceStreams stream in Enum.GetValues(typeof(ServerClient.DeviceStreams)))
            {
                Thread.Sleep(100);
                Console.WriteLine("     Adding new stream: " + stream.ToString());
                Utilities.SubscribeToStream(Device1, stream);
            }
```
Each available device stream is assigned, but any number of streams can be assigned.  
```
            Utilities.StartStreaming(Device1);
            var timer = Utilities.StartTimer(5);
            Utilities.WindowedDataReadyEvent += PullData;

            Console.WriteLine("ENTER to end program");
            Console.ReadKey();
        }
```
When the device has streaming restarted, it will begin collecting data as quickly as the server sends it.  The window size for this example is 5 seconds, but can be any length.  An internal timer will trigger an event at each expiration and then reset.  Ending the program can be done by pressing the 'Enter' key.
```
        private static void PullData()
        {
            var watch = Stopwatch.StartNew();
            var WindowData = Utilities.GrabWindow(Device1, @"C:\Readings.data");
            var pred = Predict.PredictWindow(mlContext, Model, WindowData);
            watch.Stop();
            Console.WriteLine("Time: {0} | Prediction: {1}", DateTime.Now, pred.Prediction);
            Console.WriteLine("Processing Time: {0}", watch.ElapsedMilliseconds);
        }
```
The data can be manipulated once the timer expires for the window size.  The data can be grabbed by the GrabWindow function, which also can be fed a filepath to write the readings to a data file.  Predictions can be done on the raw data using the PredictWindow function and a prediction is returned in a simple structure that contains the output bool.

### Class Documentation
[Biosensor Framework Documentation](https://github.com/Munroe-Meyer-Institute-VR-Laboratory/Biosensor-Framework/blob/main/BiosensorFramework/BiosensorFrameworkDoc.md)

### Paper Citation
Arce, W., Gehringer, J., (2021). Biosensor Framework: A C# Library for Affective Computing. Journal of Open Source Software, 6(64), 3455, https://doi.org/10.21105/joss.03455

### Contact
To address issues in the codebase, please create an issue in the repo and it will be addressed by a maintainer.  For collaboration inquiries, please address Dr. James Gehringer (james.gehringer@unmc.edu).  For technical questions, please address Walker Arce (walker.arce@unmc.edu).
### Biosensor Framework: A C# Library for Affective Computing
C# library that handles the full process of gathering biometric data from a body-worn sensor, transforming it into handcrafted feature vectors, and delivering an inferencing result in thirty-five lines of code.

### Basic Usage

```
        public static ServerClient Device1;
        static ServerClient client;
        public static MLContext mlContext;
        public static ITransformer Model;
        public static DataViewSchema ModelSchema;
        // TODO: Fill these out with your own values
        public static string APIKey = "";
        public static string ServerPath = @"";
        public static string WesadDirectory = @"";
        public static string SessionOutputPath = @"";
```
For each project, a client connection is established using the ServerClient class.  This handles the TCP initialization, command communication, and error handling.  Each device will have its own ServerClient instance, so Device1 is represented by a device name and a TCP connection.  Additionally, the API key and the path to the server executable need to be defined.  The API key can be found by following Empatica's directions from the Installation section.  The path to the server is found in the installation directory for the server (i.e. @"{installation path}\EmpaticaBLEServer.exe").

Next the Microsoft.ML variables are created.
```
        static void Main(string[] args)
        {
            mlContext = new MLContext();
            Train.RunBenchmarks(WesadDirectory, out ITransformer RegModel, out ITransformer MultiModel, out Model);
            Console.ReadKey();
```
In this example, the models are trained on the WESAD data using the RunBenchmarks function.  The best performing models for each class are output.  For this example, the binary classification model (BinModel) is used.
```
            Console.WriteLine("E4 Console Interface - Press ENTER to begin the client");

            Console.WriteLine("Step 1 - Start Empatica server");
            Utilities.StartE4ServerGUI(ServerPath);
```
The E4 server is started through a Process call.  This will open the GUI for the server, if that is not desired, the command line variant of the command can be called.
```
            client = new ServerClient();
            Console.ReadKey();
            client.StartClient();
            Utilities.StartE4Server(ServerPath, APIKey);
```
The E4 server can also be started through the command line variant as shown.  
```
            Console.WriteLine("Step 2 - Getting connected E4 devices");
            Utilities.ListDiscoveredDevices(client);
            Console.ReadKey();
            Console.WriteLine("     Available Devices:");
            Utilities.AvailableDevices.ForEach(i => Console.WriteLine("     Device Name: {0}", i));
            Console.ReadKey();
```
Listing the connected devices will show all devices connected through the Bluetooth dongle.  These are managed internally.
```
            Device1 = new ServerClient();
            Device1.StartClient();
            Device1.DeviceName = Utilities.AvailableDevices[0];
            Utilities.ConnectDevice(Device1);
            Console.ReadKey();
```
To connect a device, it needs to be assigned one of the available devices.  Once that's done, the ConnectDevice function can be called and its TCP connection will be established.
```
            Utilities.SuspendStreaming(Device1);
```
Since there are configurations to be done on the device, the streaming needs to be suspended.
```
            Console.WriteLine("Step 3 - Adding biometric data streams");

            foreach (ServerClient.DeviceStreams stream in Enum.GetValues(typeof(ServerClient.DeviceStreams)))
            {
                Thread.Sleep(100);
                Console.WriteLine("     Adding new stream: " + stream.ToString());
                Utilities.SubscribeToStream(Device1, stream);
            }
```
Each available device stream is assigned, but any number of streams can be assigned.  
```
            Utilities.StartStreaming(Device1);
            var timer = Utilities.StartTimer(5);
            Utilities.WindowedDataReadyEvent += PullData;

            Console.WriteLine("ENTER to end program");
            Console.ReadKey();
        }
```
When the device has streaming restarted, it will begin collecting data as quickly as the server sends it.  The window size for this example is 5 seconds, but can be any length.  An internal timer will trigger an event at each expiration and then reset.  Ending the program can be done by pressing the 'Enter' key.
```
        private static void PullData()
        {
            var watch = Stopwatch.StartNew();
            var WindowData = Utilities.GrabWindow(Device1, Path.Combine(SessionOutputPath, "Readings.data"));
            var pred = Predict.PredictWindow(mlContext, Model, WindowData);
            watch.Stop();
            Console.WriteLine("Time: {0} | Prediction: {1}", DateTime.Now, pred.Prediction);
            Console.WriteLine("Processing Time: {0}", watch.ElapsedMilliseconds);
        }
```
The data can be manipulated once the timer expires for the window size.  The data can be grabbed by the GrabWindow function, which also can be fed a filepath to write the readings to a data file.  Predictions can be done on the raw data using the PredictWindow function and a prediction is returned in a simple structure that contains the output bool.

### Class Documentation
[Biosensor Framework Documentation](https://github.com/Munroe-Meyer-Institute-VR-Laboratory/Biosensor-Framework/blob/main/BiosensorFramework/BiosensorFrameworkDoc.md)

### Contact
To address issues in the codebase, please create an issue in the repo and it will be addressed by a maintainer.  For collaboration inquiries, please address Dr. James Gehringer (james.gehringer@unmc.edu).  For technical questions, please address Walker Arce (walker.arce@unmc.edu).
<a name='assembly'></a>
# Biosensor Framework

## Contents

- [Analyze](#T-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-Analyze 'MMIVR.BiosensorFramework.DataProcessing.ThirdParty.DSPLib.DSP.Analyze')
  - [FindMaxAmplitude(inData)](#M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-Analyze-FindMaxAmplitude-System-Double[]- 'MMIVR.BiosensorFramework.DataProcessing.ThirdParty.DSPLib.DSP.Analyze.FindMaxAmplitude(System.Double[])')
  - [FindMaxFrequency(inData,fSpan)](#M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-Analyze-FindMaxFrequency-System-Double[],System-Double[]- 'MMIVR.BiosensorFramework.DataProcessing.ThirdParty.DSPLib.DSP.Analyze.FindMaxFrequency(System.Double[],System.Double[])')
  - [FindMaxPosition(inData)](#M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-Analyze-FindMaxPosition-System-Double[]- 'MMIVR.BiosensorFramework.DataProcessing.ThirdParty.DSPLib.DSP.Analyze.FindMaxPosition(System.Double[])')
  - [FindMean(inData,startBin,stopBin)](#M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-Analyze-FindMean-System-Double[],System-UInt32,System-UInt32- 'MMIVR.BiosensorFramework.DataProcessing.ThirdParty.DSPLib.DSP.Analyze.FindMean(System.Double[],System.UInt32,System.UInt32)')
  - [FindRms(inData,startBin,stopBin)](#M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-Analyze-FindRms-System-Double[],System-UInt32,System-UInt32- 'MMIVR.BiosensorFramework.DataProcessing.ThirdParty.DSPLib.DSP.Analyze.FindRms(System.Double[],System.UInt32,System.UInt32)')
  - [UnwrapPhaseDegrees(inPhaseDeg)](#M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-Analyze-UnwrapPhaseDegrees-System-Double[]- 'MMIVR.BiosensorFramework.DataProcessing.ThirdParty.DSPLib.DSP.Analyze.UnwrapPhaseDegrees(System.Double[])')
  - [UnwrapPhaseRadians(inPhaseRad)](#M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-Analyze-UnwrapPhaseRadians-System-Double[]- 'MMIVR.BiosensorFramework.DataProcessing.ThirdParty.DSPLib.DSP.Analyze.UnwrapPhaseRadians(System.Double[])')
- [ArrayExtensions](#T-MMIVR-BiosensorFramework-Extensions-ArrayExtensions 'MMIVR.BiosensorFramework.Extensions.ArrayExtensions')
  - [AllIndexesOf(array,value)](#M-MMIVR-BiosensorFramework-Extensions-ArrayExtensions-AllIndexesOf-System-Collections-Generic-List{MMIVR-BiosensorFramework-MachineLearningUtilities-ExtractedMultiFeatures},System-Int32- 'MMIVR.BiosensorFramework.Extensions.ArrayExtensions.AllIndexesOf(System.Collections.Generic.List{MMIVR.BiosensorFramework.MachineLearningUtilities.ExtractedMultiFeatures},System.Int32)')
  - [GetSubArray\`\`1(array,start,end)](#M-MMIVR-BiosensorFramework-Extensions-ArrayExtensions-GetSubArray``1-``0[],System-Int32,System-Int32- 'MMIVR.BiosensorFramework.Extensions.ArrayExtensions.GetSubArray``1(``0[],System.Int32,System.Int32)')
  - [Shuffle\`\`1(list)](#M-MMIVR-BiosensorFramework-Extensions-ArrayExtensions-Shuffle``1-System-Collections-Generic-IList{``0}- 'MMIVR.BiosensorFramework.Extensions.ArrayExtensions.Shuffle``1(System.Collections.Generic.IList{``0})')
  - [SplitAcc3D(array,x,y,z)](#M-MMIVR-BiosensorFramework-Extensions-ArrayExtensions-SplitAcc3D-System-Double[],System-Double[]@,System-Double[]@,System-Double[]@- 'MMIVR.BiosensorFramework.Extensions.ArrayExtensions.SplitAcc3D(System.Double[],System.Double[]@,System.Double[]@,System.Double[]@)')
  - [ToFloat(array)](#M-MMIVR-BiosensorFramework-Extensions-ArrayExtensions-ToFloat-System-Double[]- 'MMIVR.BiosensorFramework.Extensions.ArrayExtensions.ToFloat(System.Double[])')
- [BloodVolumePulse](#T-MMIVR-BiosensorFramework-DataProcessing-BloodVolumePulse 'MMIVR.BiosensorFramework.DataProcessing.BloodVolumePulse')
  - [GetHeartFeatures(Beats,SamplingRate)](#M-MMIVR-BiosensorFramework-DataProcessing-BloodVolumePulse-GetHeartFeatures-System-Int32[],System-Double- 'MMIVR.BiosensorFramework.DataProcessing.BloodVolumePulse.GetHeartFeatures(System.Int32[],System.Double)')
- [CollectedData](#T-MMIVR-BiosensorFramework-MachineLearningUtilities-DataImport-CollectedData 'MMIVR.BiosensorFramework.MachineLearningUtilities.DataImport.CollectedData')
- [ConvertComplex](#T-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-ConvertComplex 'MMIVR.BiosensorFramework.DataProcessing.ThirdParty.DSPLib.DSP.ConvertComplex')
  - [ToMagnitude(rawFFT)](#M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-ConvertComplex-ToMagnitude-System-Numerics-Complex[]- 'MMIVR.BiosensorFramework.DataProcessing.ThirdParty.DSPLib.DSP.ConvertComplex.ToMagnitude(System.Numerics.Complex[])')
  - [ToMagnitudeDBV(rawFFT)](#M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-ConvertComplex-ToMagnitudeDBV-System-Numerics-Complex[]- 'MMIVR.BiosensorFramework.DataProcessing.ThirdParty.DSPLib.DSP.ConvertComplex.ToMagnitudeDBV(System.Numerics.Complex[])')
  - [ToMagnitudeSquared(rawFFT)](#M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-ConvertComplex-ToMagnitudeSquared-System-Numerics-Complex[]- 'MMIVR.BiosensorFramework.DataProcessing.ThirdParty.DSPLib.DSP.ConvertComplex.ToMagnitudeSquared(System.Numerics.Complex[])')
  - [ToPhaseDegrees(rawFFT)](#M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-ConvertComplex-ToPhaseDegrees-System-Numerics-Complex[]- 'MMIVR.BiosensorFramework.DataProcessing.ThirdParty.DSPLib.DSP.ConvertComplex.ToPhaseDegrees(System.Numerics.Complex[])')
  - [ToPhaseRadians(rawFFT)](#M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-ConvertComplex-ToPhaseRadians-System-Numerics-Complex[]- 'MMIVR.BiosensorFramework.DataProcessing.ThirdParty.DSPLib.DSP.ConvertComplex.ToPhaseRadians(System.Numerics.Complex[])')
- [ConvertMagnitude](#T-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-ConvertMagnitude 'MMIVR.BiosensorFramework.DataProcessing.ThirdParty.DSPLib.DSP.ConvertMagnitude')
  - [ToMagnitudeDBV(magnitude)](#M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-ConvertMagnitude-ToMagnitudeDBV-System-Double[]- 'MMIVR.BiosensorFramework.DataProcessing.ThirdParty.DSPLib.DSP.ConvertMagnitude.ToMagnitudeDBV(System.Double[])')
  - [ToMagnitudeSquared(magnitude)](#M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-ConvertMagnitude-ToMagnitudeSquared-System-Double[]- 'MMIVR.BiosensorFramework.DataProcessing.ThirdParty.DSPLib.DSP.ConvertMagnitude.ToMagnitudeSquared(System.Double[])')
- [ConvertMagnitudeSquared](#T-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-ConvertMagnitudeSquared 'MMIVR.BiosensorFramework.DataProcessing.ThirdParty.DSPLib.DSP.ConvertMagnitudeSquared')
  - [ToMagnitude(magSquared)](#M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-ConvertMagnitudeSquared-ToMagnitude-System-Double[]- 'MMIVR.BiosensorFramework.DataProcessing.ThirdParty.DSPLib.DSP.ConvertMagnitudeSquared.ToMagnitude(System.Double[])')
  - [ToMagnitudeDBV(magSquared)](#M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-ConvertMagnitudeSquared-ToMagnitudeDBV-System-Double[]- 'MMIVR.BiosensorFramework.DataProcessing.ThirdParty.DSPLib.DSP.ConvertMagnitudeSquared.ToMagnitudeDBV(System.Double[])')
- [DFT](#T-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DFT 'MMIVR.BiosensorFramework.DataProcessing.ThirdParty.DSPLib.DFT')
  - [#ctor()](#M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DFT-#ctor 'MMIVR.BiosensorFramework.DataProcessing.ThirdParty.DSPLib.DFT.#ctor')
  - [IsUsingCached](#P-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DFT-IsUsingCached 'MMIVR.BiosensorFramework.DataProcessing.ThirdParty.DSPLib.DFT.IsUsingCached')
  - [Dft(timeSeries)](#M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DFT-Dft-System-Double[]- 'MMIVR.BiosensorFramework.DataProcessing.ThirdParty.DSPLib.DFT.Dft(System.Double[])')
  - [DftCached(timeSeries)](#M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DFT-DftCached-System-Double[]- 'MMIVR.BiosensorFramework.DataProcessing.ThirdParty.DSPLib.DFT.DftCached(System.Double[])')
  - [Execute(timeSeries)](#M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DFT-Execute-System-Double[]- 'MMIVR.BiosensorFramework.DataProcessing.ThirdParty.DSPLib.DFT.Execute(System.Double[])')
  - [FrequencySpan(samplingFrequencyHz)](#M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DFT-FrequencySpan-System-Double- 'MMIVR.BiosensorFramework.DataProcessing.ThirdParty.DSPLib.DFT.FrequencySpan(System.Double)')
  - [Initialize(inputDataLength,zeroPaddingLength,forceNoCache)](#M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DFT-Initialize-System-UInt32,System-UInt32,System-Boolean- 'MMIVR.BiosensorFramework.DataProcessing.ThirdParty.DSPLib.DFT.Initialize(System.UInt32,System.UInt32,System.Boolean)')
- [DataImport](#T-MMIVR-BiosensorFramework-MachineLearningUtilities-DataImport 'MMIVR.BiosensorFramework.MachineLearningUtilities.DataImport')
  - [ConvertRawToBin(mlContext,RawFeatures,TrainTestRatio)](#M-MMIVR-BiosensorFramework-MachineLearningUtilities-DataImport-ConvertRawToBin-Microsoft-ML-MLContext,System-Collections-Generic-List{System-Tuple{System-Double[],System-Int32}},System-Double- 'MMIVR.BiosensorFramework.MachineLearningUtilities.DataImport.ConvertRawToBin(Microsoft.ML.MLContext,System.Collections.Generic.List{System.Tuple{System.Double[],System.Int32}},System.Double)')
  - [ConvertRawToMulti(mlContext,RawFeatures,TrainTestRatio)](#M-MMIVR-BiosensorFramework-MachineLearningUtilities-DataImport-ConvertRawToMulti-Microsoft-ML-MLContext,System-Collections-Generic-List{System-Tuple{System-Double[],System-Int32}},System-Double- 'MMIVR.BiosensorFramework.MachineLearningUtilities.DataImport.ConvertRawToMulti(Microsoft.ML.MLContext,System.Collections.Generic.List{System.Tuple{System.Double[],System.Int32}},System.Double)')
  - [LoadCollectedDataset(WesadDirectory,DirectoryPath,SearchPattern,WindowSize)](#M-MMIVR-BiosensorFramework-MachineLearningUtilities-DataImport-LoadCollectedDataset-System-String,System-String,System-String,System-Int32- 'MMIVR.BiosensorFramework.MachineLearningUtilities.DataImport.LoadCollectedDataset(System.String,System.String,System.String,System.Int32)')
  - [LoadDataset(DirectoryPath)](#M-MMIVR-BiosensorFramework-MachineLearningUtilities-DataImport-LoadDataset-System-String- 'MMIVR.BiosensorFramework.MachineLearningUtilities.DataImport.LoadDataset(System.String)')
  - [LoadFile(Filepath,WindowSize)](#M-MMIVR-BiosensorFramework-MachineLearningUtilities-DataImport-LoadFile-System-String,System-Int32- 'MMIVR.BiosensorFramework.MachineLearningUtilities.DataImport.LoadFile(System.String,System.Int32)')
  - [SchmidtDatasetPipeline(DirectoryPath,mlContext,MultiClass,BinClass,RegClass,WindowSize,TrainTestRatio)](#M-MMIVR-BiosensorFramework-MachineLearningUtilities-DataImport-SchmidtDatasetPipeline-System-String,Microsoft-ML-MLContext,Microsoft-ML-DataOperationsCatalog-TrainTestData@,Microsoft-ML-DataOperationsCatalog-TrainTestData@,Microsoft-ML-DataOperationsCatalog-TrainTestData@,System-Int32,System-Double- 'MMIVR.BiosensorFramework.MachineLearningUtilities.DataImport.SchmidtDatasetPipeline(System.String,Microsoft.ML.MLContext,Microsoft.ML.DataOperationsCatalog.TrainTestData@,Microsoft.ML.DataOperationsCatalog.TrainTestData@,Microsoft.ML.DataOperationsCatalog.TrainTestData@,System.Int32,System.Double)')
- [DeviceStreams](#T-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-DeviceStreams 'MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient.DeviceStreams')
  - [ACC](#F-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-DeviceStreams-ACC 'MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient.DeviceStreams.ACC')
  - [BAT](#F-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-DeviceStreams-BAT 'MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient.DeviceStreams.BAT')
  - [BVP](#F-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-DeviceStreams-BVP 'MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient.DeviceStreams.BVP')
  - [GSR](#F-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-DeviceStreams-GSR 'MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient.DeviceStreams.GSR')
  - [HR](#F-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-DeviceStreams-HR 'MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient.DeviceStreams.HR')
  - [IBI](#F-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-DeviceStreams-IBI 'MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient.DeviceStreams.IBI')
  - [TAG](#F-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-DeviceStreams-TAG 'MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient.DeviceStreams.TAG')
  - [TMP](#F-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-DeviceStreams-TMP 'MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient.DeviceStreams.TMP')
- [ExtractedBinFeatures](#T-MMIVR-BiosensorFramework-MachineLearningUtilities-ExtractedBinFeatures 'MMIVR.BiosensorFramework.MachineLearningUtilities.ExtractedBinFeatures')
  - [Features](#F-MMIVR-BiosensorFramework-MachineLearningUtilities-ExtractedBinFeatures-Features 'MMIVR.BiosensorFramework.MachineLearningUtilities.ExtractedBinFeatures.Features')
  - [Label](#F-MMIVR-BiosensorFramework-MachineLearningUtilities-ExtractedBinFeatures-Label 'MMIVR.BiosensorFramework.MachineLearningUtilities.ExtractedBinFeatures.Label')
- [ExtractedFeatures](#T-MMIVR-BiosensorFramework-MachineLearningUtilities-ExtractedFeatures 'MMIVR.BiosensorFramework.MachineLearningUtilities.ExtractedFeatures')
  - [StressFeatures](#F-MMIVR-BiosensorFramework-MachineLearningUtilities-ExtractedFeatures-StressFeatures 'MMIVR.BiosensorFramework.MachineLearningUtilities.ExtractedFeatures.StressFeatures')
- [ExtractedMultiFeatures](#T-MMIVR-BiosensorFramework-MachineLearningUtilities-ExtractedMultiFeatures 'MMIVR.BiosensorFramework.MachineLearningUtilities.ExtractedMultiFeatures')
  - [Result](#F-MMIVR-BiosensorFramework-MachineLearningUtilities-ExtractedMultiFeatures-Result 'MMIVR.BiosensorFramework.MachineLearningUtilities.ExtractedMultiFeatures.Result')
  - [StressFeatures](#F-MMIVR-BiosensorFramework-MachineLearningUtilities-ExtractedMultiFeatures-StressFeatures 'MMIVR.BiosensorFramework.MachineLearningUtilities.ExtractedMultiFeatures.StressFeatures')
- [ExtractedRegFeatures](#T-MMIVR-BiosensorFramework-MachineLearningUtilities-ExtractedRegFeatures 'MMIVR.BiosensorFramework.MachineLearningUtilities.ExtractedRegFeatures')
  - [Result](#F-MMIVR-BiosensorFramework-MachineLearningUtilities-ExtractedRegFeatures-Result 'MMIVR.BiosensorFramework.MachineLearningUtilities.ExtractedRegFeatures.Result')
  - [StressFeatures](#F-MMIVR-BiosensorFramework-MachineLearningUtilities-ExtractedRegFeatures-StressFeatures 'MMIVR.BiosensorFramework.MachineLearningUtilities.ExtractedRegFeatures.StressFeatures')
- [FFT](#T-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-FFT 'MMIVR.BiosensorFramework.DataProcessing.ThirdParty.DSPLib.FFT')
  - [#ctor()](#M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-FFT-#ctor 'MMIVR.BiosensorFramework.DataProcessing.ThirdParty.DSPLib.FFT.#ctor')
  - [BitReverse()](#M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-FFT-BitReverse-System-UInt32,System-UInt32- 'MMIVR.BiosensorFramework.DataProcessing.ThirdParty.DSPLib.FFT.BitReverse(System.UInt32,System.UInt32)')
  - [Execute(timeSeries)](#M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-FFT-Execute-System-Double[]- 'MMIVR.BiosensorFramework.DataProcessing.ThirdParty.DSPLib.FFT.Execute(System.Double[])')
  - [FrequencySpan(samplingFrequencyHz)](#M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-FFT-FrequencySpan-System-Double- 'MMIVR.BiosensorFramework.DataProcessing.ThirdParty.DSPLib.FFT.FrequencySpan(System.Double)')
  - [Initialize(inputDataLength,zeroPaddingLength)](#M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-FFT-Initialize-System-UInt32,System-UInt32- 'MMIVR.BiosensorFramework.DataProcessing.ThirdParty.DSPLib.FFT.Initialize(System.UInt32,System.UInt32)')
- [FIRFilterImplementation](#T-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-FIRFilterImplementation 'MMIVR.BiosensorFramework.DataProcessing.ThirdParty.DSPLib.FIRFilterImplementation')
- [FeatureExtraction](#T-MMIVR-BiosensorFramework-DataProcessing-FeatureExtraction 'MMIVR.BiosensorFramework.DataProcessing.FeatureExtraction')
  - [FindLocalMaxima(array,Threshold,Samples)](#M-MMIVR-BiosensorFramework-DataProcessing-FeatureExtraction-FindLocalMaxima-System-Double[],System-Double,System-Int32- 'MMIVR.BiosensorFramework.DataProcessing.FeatureExtraction.FindLocalMaxima(System.Double[],System.Double,System.Int32)')
  - [PSI(Accelerations,AverageIndex,MaxIndex,WindowWidth)](#M-MMIVR-BiosensorFramework-DataProcessing-FeatureExtraction-PSI-System-Collections-Generic-List{System-Double},System-Double@,System-Double@,System-Int32- 'MMIVR.BiosensorFramework.DataProcessing.FeatureExtraction.PSI(System.Collections.Generic.List{System.Double},System.Double@,System.Double@,System.Int32)')
  - [PSIVector(Accelerations,AverageIndex,MaxIndex,WindowWidth)](#M-MMIVR-BiosensorFramework-DataProcessing-FeatureExtraction-PSIVector-System-Collections-Generic-List{System-Collections-Generic-List{System-Double}},System-Double@,System-Double@,System-Int32- 'MMIVR.BiosensorFramework.DataProcessing.FeatureExtraction.PSIVector(System.Collections.Generic.List{System.Collections.Generic.List{System.Double}},System.Double@,System.Double@,System.Int32)')
  - [RemoveBias(AccelerationVector,BiasX,BiasY,BiasZ)](#M-MMIVR-BiosensorFramework-DataProcessing-FeatureExtraction-RemoveBias-System-Numerics-Vector3,System-Double,System-Double,System-Double- 'MMIVR.BiosensorFramework.DataProcessing.FeatureExtraction.RemoveBias(System.Numerics.Vector3,System.Double,System.Double,System.Double)')
  - [RootMeanSquare(FirstTerm,Summation)](#M-MMIVR-BiosensorFramework-DataProcessing-FeatureExtraction-RootMeanSquare-System-Double,System-Double- 'MMIVR.BiosensorFramework.DataProcessing.FeatureExtraction.RootMeanSquare(System.Double,System.Double)')
  - [SignalAbsolute(AUC)](#M-MMIVR-BiosensorFramework-DataProcessing-FeatureExtraction-SignalAbsolute-System-Double- 'MMIVR.BiosensorFramework.DataProcessing.FeatureExtraction.SignalAbsolute(System.Double)')
  - [SignalCorrelation(Signal,TimeLabels)](#M-MMIVR-BiosensorFramework-DataProcessing-FeatureExtraction-SignalCorrelation-System-Double[],System-Double[]- 'MMIVR.BiosensorFramework.DataProcessing.FeatureExtraction.SignalCorrelation(System.Double[],System.Double[])')
  - [SignalDynamicRange(Signal)](#M-MMIVR-BiosensorFramework-DataProcessing-FeatureExtraction-SignalDynamicRange-System-Double[]- 'MMIVR.BiosensorFramework.DataProcessing.FeatureExtraction.SignalDynamicRange(System.Double[])')
  - [SignalFreqBandEnergies(Signal,Bands,SamplingFrequency,Order)](#M-MMIVR-BiosensorFramework-DataProcessing-FeatureExtraction-SignalFreqBandEnergies-System-Double[],System-Collections-Generic-List{System-Tuple{System-Double,System-Double}},System-Double,System-Int32- 'MMIVR.BiosensorFramework.DataProcessing.FeatureExtraction.SignalFreqBandEnergies(System.Double[],System.Collections.Generic.List{System.Tuple{System.Double,System.Double}},System.Double,System.Int32)')
  - [SignalFreqNormalize(Freq)](#M-MMIVR-BiosensorFramework-DataProcessing-FeatureExtraction-SignalFreqNormalize-System-Double[]- 'MMIVR.BiosensorFramework.DataProcessing.FeatureExtraction.SignalFreqNormalize(System.Double[])')
  - [SignalFreqRatio(LowFreq,HighFreq)](#M-MMIVR-BiosensorFramework-DataProcessing-FeatureExtraction-SignalFreqRatio-System-Double,System-Double- 'MMIVR.BiosensorFramework.DataProcessing.FeatureExtraction.SignalFreqRatio(System.Double,System.Double)')
  - [SignalFreqSummation(Freq)](#M-MMIVR-BiosensorFramework-DataProcessing-FeatureExtraction-SignalFreqSummation-System-Collections-Generic-List{System-Double[]}- 'MMIVR.BiosensorFramework.DataProcessing.FeatureExtraction.SignalFreqSummation(System.Collections.Generic.List{System.Double[]})')
  - [SignalIntegral(Signal)](#M-MMIVR-BiosensorFramework-DataProcessing-FeatureExtraction-SignalIntegral-System-Double[]- 'MMIVR.BiosensorFramework.DataProcessing.FeatureExtraction.SignalIntegral(System.Double[])')
  - [SignalMean(Signal)](#M-MMIVR-BiosensorFramework-DataProcessing-FeatureExtraction-SignalMean-System-Collections-Generic-List{System-Single}- 'MMIVR.BiosensorFramework.DataProcessing.FeatureExtraction.SignalMean(System.Collections.Generic.List{System.Single})')
  - [SignalMinMax(Signal)](#M-MMIVR-BiosensorFramework-DataProcessing-FeatureExtraction-SignalMinMax-System-Collections-Generic-List{System-Single}- 'MMIVR.BiosensorFramework.DataProcessing.FeatureExtraction.SignalMinMax(System.Collections.Generic.List{System.Single})')
  - [SignalPeakFrequency(Signal,SamplingRate)](#M-MMIVR-BiosensorFramework-DataProcessing-FeatureExtraction-SignalPeakFrequency-System-Double[],System-Double- 'MMIVR.BiosensorFramework.DataProcessing.FeatureExtraction.SignalPeakFrequency(System.Double[],System.Double)')
  - [SignalPercentile(Signal,percentile)](#M-MMIVR-BiosensorFramework-DataProcessing-FeatureExtraction-SignalPercentile-System-Double[],System-Int32- 'MMIVR.BiosensorFramework.DataProcessing.FeatureExtraction.SignalPercentile(System.Double[],System.Int32)')
  - [SignalRMS(Signal)](#M-MMIVR-BiosensorFramework-DataProcessing-FeatureExtraction-SignalRMS-System-Double[]- 'MMIVR.BiosensorFramework.DataProcessing.FeatureExtraction.SignalRMS(System.Double[])')
  - [SignalRelativePower(Energies,Resolution,Epsilon)](#M-MMIVR-BiosensorFramework-DataProcessing-FeatureExtraction-SignalRelativePower-System-Collections-Generic-List{System-Double[]},System-Double,System-Double- 'MMIVR.BiosensorFramework.DataProcessing.FeatureExtraction.SignalRelativePower(System.Collections.Generic.List{System.Double[]},System.Double,System.Double)')
  - [SignalSlope(Signal,TimeLabels)](#M-MMIVR-BiosensorFramework-DataProcessing-FeatureExtraction-SignalSlope-System-Double[],System-Double[]- 'MMIVR.BiosensorFramework.DataProcessing.FeatureExtraction.SignalSlope(System.Double[],System.Double[])')
  - [SignalStandardDeviation(Signal)](#M-MMIVR-BiosensorFramework-DataProcessing-FeatureExtraction-SignalStandardDeviation-System-Collections-Generic-List{System-Single}- 'MMIVR.BiosensorFramework.DataProcessing.FeatureExtraction.SignalStandardDeviation(System.Collections.Generic.List{System.Single})')
  - [SquareSummation(Value)](#M-MMIVR-BiosensorFramework-DataProcessing-FeatureExtraction-SquareSummation-System-Double- 'MMIVR.BiosensorFramework.DataProcessing.FeatureExtraction.SquareSummation(System.Double)')
  - [SquareSummation(Value)](#M-MMIVR-BiosensorFramework-DataProcessing-FeatureExtraction-SquareSummation-System-Numerics-Vector3- 'MMIVR.BiosensorFramework.DataProcessing.FeatureExtraction.SquareSummation(System.Numerics.Vector3)')
- [Generate](#T-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-Generate 'MMIVR.BiosensorFramework.DataProcessing.ThirdParty.DSPLib.DSP.Generate')
  - [LinSpace(startVal,stopVal,points)](#M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-Generate-LinSpace-System-Double,System-Double,System-UInt32- 'MMIVR.BiosensorFramework.DataProcessing.ThirdParty.DSPLib.DSP.Generate.LinSpace(System.Double,System.Double,System.UInt32)')
  - [NoisePsd(amplitudePsd (Vrms / rt-Hz)(Vrms/rt-Hz),samplingFrequencyHz,points)](#M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-Generate-NoisePsd-System-Double,System-Double,System-UInt32- 'MMIVR.BiosensorFramework.DataProcessing.ThirdParty.DSPLib.DSP.Generate.NoisePsd(System.Double,System.Double,System.UInt32)')
  - [NoiseRms(amplitudeVrms,points,dcV)](#M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-Generate-NoiseRms-System-Double,System-UInt32,System-Double- 'MMIVR.BiosensorFramework.DataProcessing.ThirdParty.DSPLib.DSP.Generate.NoiseRms(System.Double,System.UInt32,System.Double)')
  - [ToneCycles(amplitudeVrms,cycles,points,dcV,phaseDeg)](#M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-Generate-ToneCycles-System-Double,System-Double,System-UInt32,System-Double,System-Double- 'MMIVR.BiosensorFramework.DataProcessing.ThirdParty.DSPLib.DSP.Generate.ToneCycles(System.Double,System.Double,System.UInt32,System.Double,System.Double)')
  - [ToneSampling(amplitudeVrms,frequencyHz,samplingFrequencyHz,points,dcV,phaseDeg)](#M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-Generate-ToneSampling-System-Double,System-Double,System-Double,System-UInt32,System-Double,System-Double- 'MMIVR.BiosensorFramework.DataProcessing.ThirdParty.DSPLib.DSP.Generate.ToneSampling(System.Double,System.Double,System.Double,System.UInt32,System.Double,System.Double)')
- [HealeyStressDetection](#T-MMIVR-BiosensorFramework-DataProcessing-HealeyStressDetection 'MMIVR.BiosensorFramework.DataProcessing.HealeyStressDetection')
  - [ConstructPeaks(Peaks,ZeroCrossings)](#M-MMIVR-BiosensorFramework-DataProcessing-HealeyStressDetection-ConstructPeaks-System-Collections-Generic-List{System-Collections-Generic-List{System-Int32}},System-Tuple{System-Int32[],System-Int32[]}- 'MMIVR.BiosensorFramework.DataProcessing.HealeyStressDetection.ConstructPeaks(System.Collections.Generic.List{System.Collections.Generic.List{System.Int32}},System.Tuple{System.Int32[],System.Int32[]})')
  - [FindAbove(array,Threshold,SamplingRate)](#M-MMIVR-BiosensorFramework-DataProcessing-HealeyStressDetection-FindAbove-System-Double[],System-Double,System-Double- 'MMIVR.BiosensorFramework.DataProcessing.HealeyStressDetection.FindAbove(System.Double[],System.Double,System.Double)')
  - [FindNearestMinus(PeakStart,ZeroCrossings)](#M-MMIVR-BiosensorFramework-DataProcessing-HealeyStressDetection-FindNearestMinus-System-Int32,System-Int32[]- 'MMIVR.BiosensorFramework.DataProcessing.HealeyStressDetection.FindNearestMinus(System.Int32,System.Int32[])')
  - [FindNearestPlus(PeakStart,ZeroCrossings)](#M-MMIVR-BiosensorFramework-DataProcessing-HealeyStressDetection-FindNearestPlus-System-Int32,System-Int32[]- 'MMIVR.BiosensorFramework.DataProcessing.HealeyStressDetection.FindNearestPlus(System.Int32,System.Int32[])')
  - [FindZeroCrossings(array)](#M-MMIVR-BiosensorFramework-DataProcessing-HealeyStressDetection-FindZeroCrossings-System-Double[]- 'MMIVR.BiosensorFramework.DataProcessing.HealeyStressDetection.FindZeroCrossings(System.Double[])')
  - [FiniteDifference(Signal)](#M-MMIVR-BiosensorFramework-DataProcessing-HealeyStressDetection-FiniteDifference-System-Double[]@- 'MMIVR.BiosensorFramework.DataProcessing.HealeyStressDetection.FiniteDifference(System.Double[]@)')
  - [GetPeakFeatures(Signal,SamplingRate,CutoffFreq,Threshold)](#M-MMIVR-BiosensorFramework-DataProcessing-HealeyStressDetection-GetPeakFeatures-System-Double[],System-Double,System-Double,System-Double- 'MMIVR.BiosensorFramework.DataProcessing.HealeyStressDetection.GetPeakFeatures(System.Double[],System.Double,System.Double,System.Double)')
  - [MergeOverlappingPeaks(Peaks)](#M-MMIVR-BiosensorFramework-DataProcessing-HealeyStressDetection-MergeOverlappingPeaks-System-Collections-Generic-List{System-Int32[]}- 'MMIVR.BiosensorFramework.DataProcessing.HealeyStressDetection.MergeOverlappingPeaks(System.Collections.Generic.List{System.Int32[]})')
  - [SumFeatures(PeakFeatures)](#M-MMIVR-BiosensorFramework-DataProcessing-HealeyStressDetection-SumFeatures-System-Collections-Generic-List{System-Double[]}- 'MMIVR.BiosensorFramework.DataProcessing.HealeyStressDetection.SumFeatures(System.Collections.Generic.List{System.Double[]})')
- [Math](#T-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-Math 'MMIVR.BiosensorFramework.DataProcessing.ThirdParty.DSPLib.DSP.Math')
  - [Add()](#M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-Math-Add-System-Double[],System-Double[]- 'MMIVR.BiosensorFramework.DataProcessing.ThirdParty.DSPLib.DSP.Math.Add(System.Double[],System.Double[])')
  - [Add()](#M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-Math-Add-System-Double[],System-Double- 'MMIVR.BiosensorFramework.DataProcessing.ThirdParty.DSPLib.DSP.Math.Add(System.Double[],System.Double)')
  - [Divide()](#M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-Math-Divide-System-Double[],System-Double[]- 'MMIVR.BiosensorFramework.DataProcessing.ThirdParty.DSPLib.DSP.Math.Divide(System.Double[],System.Double[])')
  - [Divide()](#M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-Math-Divide-System-Double[],System-Double- 'MMIVR.BiosensorFramework.DataProcessing.ThirdParty.DSPLib.DSP.Math.Divide(System.Double[],System.Double)')
  - [Log10()](#M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-Math-Log10-System-Double[]- 'MMIVR.BiosensorFramework.DataProcessing.ThirdParty.DSPLib.DSP.Math.Log10(System.Double[])')
  - [Multiply()](#M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-Math-Multiply-System-Double[],System-Double[]- 'MMIVR.BiosensorFramework.DataProcessing.ThirdParty.DSPLib.DSP.Math.Multiply(System.Double[],System.Double[])')
  - [Multiply()](#M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-Math-Multiply-System-Double[],System-Double- 'MMIVR.BiosensorFramework.DataProcessing.ThirdParty.DSPLib.DSP.Math.Multiply(System.Double[],System.Double)')
  - [RemoveMean()](#M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-Math-RemoveMean-System-Double[]- 'MMIVR.BiosensorFramework.DataProcessing.ThirdParty.DSPLib.DSP.Math.RemoveMean(System.Double[])')
  - [Sqrt()](#M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-Math-Sqrt-System-Double[]- 'MMIVR.BiosensorFramework.DataProcessing.ThirdParty.DSPLib.DSP.Math.Sqrt(System.Double[])')
  - [Square()](#M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-Math-Square-System-Double[]- 'MMIVR.BiosensorFramework.DataProcessing.ThirdParty.DSPLib.DSP.Math.Square(System.Double[])')
  - [Subtract()](#M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-Math-Subtract-System-Double[],System-Double[]- 'MMIVR.BiosensorFramework.DataProcessing.ThirdParty.DSPLib.DSP.Math.Subtract(System.Double[],System.Double[])')
  - [Subtract()](#M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-Math-Subtract-System-Double[],System-Double- 'MMIVR.BiosensorFramework.DataProcessing.ThirdParty.DSPLib.DSP.Math.Subtract(System.Double[],System.Double)')
- [Predict](#T-MMIVR-BiosensorFramework-MachineLearningUtilities-Predict 'MMIVR.BiosensorFramework.MachineLearningUtilities.Predict')
  - [MakeBinPrediction(mlContext,LiveData,Model)](#M-MMIVR-BiosensorFramework-MachineLearningUtilities-Predict-MakeBinPrediction-Microsoft-ML-MLContext,MMIVR-BiosensorFramework-MachineLearningUtilities-ExtractedBinFeatures,Microsoft-ML-ITransformer- 'MMIVR.BiosensorFramework.MachineLearningUtilities.Predict.MakeBinPrediction(Microsoft.ML.MLContext,MMIVR.BiosensorFramework.MachineLearningUtilities.ExtractedBinFeatures,Microsoft.ML.ITransformer)')
  - [MakeMultiPrediction(mlContext,LiveData,Model)](#M-MMIVR-BiosensorFramework-MachineLearningUtilities-Predict-MakeMultiPrediction-Microsoft-ML-MLContext,MMIVR-BiosensorFramework-MachineLearningUtilities-ExtractedMultiFeatures,Microsoft-ML-ITransformer- 'MMIVR.BiosensorFramework.MachineLearningUtilities.Predict.MakeMultiPrediction(Microsoft.ML.MLContext,MMIVR.BiosensorFramework.MachineLearningUtilities.ExtractedMultiFeatures,Microsoft.ML.ITransformer)')
  - [MakeRegPrediction(mlContext,LiveData,Model)](#M-MMIVR-BiosensorFramework-MachineLearningUtilities-Predict-MakeRegPrediction-Microsoft-ML-MLContext,Microsoft-ML-IDataView,Microsoft-ML-ITransformer- 'MMIVR.BiosensorFramework.MachineLearningUtilities.Predict.MakeRegPrediction(Microsoft.ML.MLContext,Microsoft.ML.IDataView,Microsoft.ML.ITransformer)')
  - [PredictWindow(mlContext,Model,WindowReadings)](#M-MMIVR-BiosensorFramework-MachineLearningUtilities-Predict-PredictWindow-Microsoft-ML-MLContext,Microsoft-ML-ITransformer,System-Collections-Generic-List{System-Collections-Generic-List{System-Double}}- 'MMIVR.BiosensorFramework.MachineLearningUtilities.Predict.PredictWindow(Microsoft.ML.MLContext,Microsoft.ML.ITransformer,System.Collections.Generic.List{System.Collections.Generic.List{System.Double}})')
- [PredictionBinResult](#T-MMIVR-BiosensorFramework-MachineLearningUtilities-PredictionBinResult 'MMIVR.BiosensorFramework.MachineLearningUtilities.PredictionBinResult')
  - [Prediction](#F-MMIVR-BiosensorFramework-MachineLearningUtilities-PredictionBinResult-Prediction 'MMIVR.BiosensorFramework.MachineLearningUtilities.PredictionBinResult.Prediction')
  - [Probability](#F-MMIVR-BiosensorFramework-MachineLearningUtilities-PredictionBinResult-Probability 'MMIVR.BiosensorFramework.MachineLearningUtilities.PredictionBinResult.Probability')
  - [Score](#F-MMIVR-BiosensorFramework-MachineLearningUtilities-PredictionBinResult-Score 'MMIVR.BiosensorFramework.MachineLearningUtilities.PredictionBinResult.Score')
- [PredictionMultiResult](#T-MMIVR-BiosensorFramework-MachineLearningUtilities-PredictionMultiResult 'MMIVR.BiosensorFramework.MachineLearningUtilities.PredictionMultiResult')
  - [Result](#F-MMIVR-BiosensorFramework-MachineLearningUtilities-PredictionMultiResult-Result 'MMIVR.BiosensorFramework.MachineLearningUtilities.PredictionMultiResult.Result')
- [PredictionRegResult](#T-MMIVR-BiosensorFramework-MachineLearningUtilities-PredictionRegResult 'MMIVR.BiosensorFramework.MachineLearningUtilities.PredictionRegResult')
  - [Result](#F-MMIVR-BiosensorFramework-MachineLearningUtilities-PredictionRegResult-Result 'MMIVR.BiosensorFramework.MachineLearningUtilities.PredictionRegResult.Result')
- [ProcessFFT](#T-MMIVR-BiosensorFramework-DataProcessing-ProcessFFT 'MMIVR.BiosensorFramework.DataProcessing.ProcessFFT')
  - [ProcessSignal(Signal,SamplingRate,FreqSpan,MagSpectrum)](#M-MMIVR-BiosensorFramework-DataProcessing-ProcessFFT-ProcessSignal-System-Double[]@,System-Double,System-Double[]@,System-Double[]@- 'MMIVR.BiosensorFramework.DataProcessing.ProcessFFT.ProcessSignal(System.Double[]@,System.Double,System.Double[]@,System.Double[]@)')
- [ScaleFactor](#T-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-Window-ScaleFactor 'MMIVR.BiosensorFramework.DataProcessing.ThirdParty.DSPLib.DSP.Window.ScaleFactor')
  - [NENBW(windowCoefficients)](#M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-Window-ScaleFactor-NENBW-System-Double[]- 'MMIVR.BiosensorFramework.DataProcessing.ThirdParty.DSPLib.DSP.Window.ScaleFactor.NENBW(System.Double[])')
  - [Noise(windowCoefficients,samplingFrequencyHz)](#M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-Window-ScaleFactor-Noise-System-Double[],System-Double- 'MMIVR.BiosensorFramework.DataProcessing.ThirdParty.DSPLib.DSP.Window.ScaleFactor.Noise(System.Double[],System.Double)')
  - [Signal(windowCoefficients)](#M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-Window-ScaleFactor-Signal-System-Double[]- 'MMIVR.BiosensorFramework.DataProcessing.ThirdParty.DSPLib.DSP.Window.ScaleFactor.Signal(System.Double[])')
- [ServerClient](#T-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient 'MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient')
  - [ACCReadings3D](#F-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-ACCReadings3D 'MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient.ACCReadings3D')
  - [ACCReadingsX](#F-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-ACCReadingsX 'MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient.ACCReadingsX')
  - [ACCReadingsY](#F-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-ACCReadingsY 'MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient.ACCReadingsY')
  - [ACCReadingsZ](#F-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-ACCReadingsZ 'MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient.ACCReadingsZ')
  - [ACCTimestamps](#F-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-ACCTimestamps 'MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient.ACCTimestamps')
  - [BATReadings](#F-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-BATReadings 'MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient.BATReadings')
  - [BATTimestamps](#F-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-BATTimestamps 'MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient.BATTimestamps')
  - [BVPReadings](#F-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-BVPReadings 'MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient.BVPReadings')
  - [BVPTimestamps](#F-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-BVPTimestamps 'MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient.BVPTimestamps')
  - [Buffer](#F-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-Buffer 'MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient.Buffer')
  - [BufferSize](#F-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-BufferSize 'MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient.BufferSize')
  - [ConnectDone](#F-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-ConnectDone 'MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient.ConnectDone')
  - [DeviceName](#F-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-DeviceName 'MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient.DeviceName')
  - [GSRReadings](#F-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-GSRReadings 'MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient.GSRReadings')
  - [GSRTimestamps](#F-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-GSRTimestamps 'MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient.GSRTimestamps')
  - [HRReadings](#F-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-HRReadings 'MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient.HRReadings')
  - [HRTimestamps](#F-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-HRTimestamps 'MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient.HRTimestamps')
  - [IBIReadings](#F-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-IBIReadings 'MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient.IBIReadings')
  - [IBITimestamps](#F-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-IBITimestamps 'MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient.IBITimestamps')
  - [ReceiveDone](#F-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-ReceiveDone 'MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient.ReceiveDone')
  - [Sb](#F-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-Sb 'MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient.Sb')
  - [SendDone](#F-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-SendDone 'MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient.SendDone')
  - [SocketConnection](#F-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-SocketConnection 'MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient.SocketConnection')
  - [SubscribedStreams](#F-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-SubscribedStreams 'MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient.SubscribedStreams')
  - [TAGReadings](#F-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-TAGReadings 'MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient.TAGReadings')
  - [TAGTimestamps](#F-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-TAGTimestamps 'MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient.TAGTimestamps')
  - [TMPReadings](#F-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-TMPReadings 'MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient.TMPReadings')
  - [TMPTimestamps](#F-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-TMPTimestamps 'MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient.TMPTimestamps')
  - [ConnectCallback(ar)](#M-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-ConnectCallback-System-IAsyncResult- 'MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient.ConnectCallback(System.IAsyncResult)')
  - [Dispose()](#M-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-Dispose 'MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient.Dispose')
  - [HandleDataStream(DataPacket,response)](#M-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-HandleDataStream-System-String[],System-String- 'MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient.HandleDataStream(System.String[],System.String)')
  - [HandleDiscoverResponse(responses)](#M-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-HandleDiscoverResponse-System-String[]- 'MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient.HandleDiscoverResponse(System.String[])')
  - [HandleDiscoverResponseBTLE(responses)](#M-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-HandleDiscoverResponseBTLE-System-String[]- 'MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient.HandleDiscoverResponseBTLE(System.String[])')
  - [HandleErrorCodes(ErrorResponse)](#M-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-HandleErrorCodes-System-String[]- 'MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient.HandleErrorCodes(System.String[])')
  - [HandleResponseFromEmpaticaBLEServer(response)](#M-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-HandleResponseFromEmpaticaBLEServer-System-String- 'MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient.HandleResponseFromEmpaticaBLEServer(System.String)')
  - [Receive()](#M-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-Receive 'MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient.Receive')
  - [ReceiveCallback(ar)](#M-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-ReceiveCallback-System-IAsyncResult- 'MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient.ReceiveCallback(System.IAsyncResult)')
  - [Send(client,data)](#M-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-Send-System-Net-Sockets-Socket,System-String- 'MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient.Send(System.Net.Sockets.Socket,System.String)')
  - [SendCallback(ar)](#M-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-SendCallback-System-IAsyncResult- 'MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient.SendCallback(System.IAsyncResult)')
  - [StartClient()](#M-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-StartClient 'MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient.StartClient')
  - [TagStressEvent()](#M-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-TagStressEvent 'MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient.TagStressEvent')
- [SignalProcessing](#T-MMIVR-BiosensorFramework-InputPipeline-SignalProcessing 'MMIVR.BiosensorFramework.InputPipeline.SignalProcessing')
  - [ProcessAccSignal(AccSignal,SamplingRate)](#M-MMIVR-BiosensorFramework-InputPipeline-SignalProcessing-ProcessAccSignal-System-Double[],System-Double- 'MMIVR.BiosensorFramework.InputPipeline.SignalProcessing.ProcessAccSignal(System.Double[],System.Double)')
  - [ProcessEdaSignal(EdaSignal,SamplingRate)](#M-MMIVR-BiosensorFramework-InputPipeline-SignalProcessing-ProcessEdaSignal-System-Double[],System-Double- 'MMIVR.BiosensorFramework.InputPipeline.SignalProcessing.ProcessEdaSignal(System.Double[],System.Double)')
  - [ProcessPpgSignal(PpgSignal,SamplingRate,Threshold,PeakSize)](#M-MMIVR-BiosensorFramework-InputPipeline-SignalProcessing-ProcessPpgSignal-System-Double[],System-Double,System-Double,System-Int32- 'MMIVR.BiosensorFramework.InputPipeline.SignalProcessing.ProcessPpgSignal(System.Double[],System.Double,System.Double,System.Int32)')
  - [ProcessTmpSignal(TempSignal)](#M-MMIVR-BiosensorFramework-InputPipeline-SignalProcessing-ProcessTmpSignal-System-Double[]- 'MMIVR.BiosensorFramework.InputPipeline.SignalProcessing.ProcessTmpSignal(System.Double[])')
- [TarvainenDetrending](#T-MMIVR-BiosensorFramework-DataProcessing-TarvainenDetrending 'MMIVR.BiosensorFramework.DataProcessing.TarvainenDetrending')
  - [GetResidual(TrendedSignal,DetrendedSignal)](#M-MMIVR-BiosensorFramework-DataProcessing-TarvainenDetrending-GetResidual-System-Double[],System-Double[]- 'MMIVR.BiosensorFramework.DataProcessing.TarvainenDetrending.GetResidual(System.Double[],System.Double[])')
  - [RemoveTrend(InputSignal,Lambda,Filter)](#M-MMIVR-BiosensorFramework-DataProcessing-TarvainenDetrending-RemoveTrend-System-Double[],System-Int32,System-Double[]- 'MMIVR.BiosensorFramework.DataProcessing.TarvainenDetrending.RemoveTrend(System.Double[],System.Int32,System.Double[])')
- [Train](#T-MMIVR-BiosensorFramework-MachineLearningUtilities-Train 'MMIVR.BiosensorFramework.MachineLearningUtilities.Train')
  - [BuildAndTrainBinClassModels(mlContext,TrainingSet)](#M-MMIVR-BiosensorFramework-MachineLearningUtilities-Train-BuildAndTrainBinClassModels-Microsoft-ML-MLContext,Microsoft-ML-IDataView- 'MMIVR.BiosensorFramework.MachineLearningUtilities.Train.BuildAndTrainBinClassModels(Microsoft.ML.MLContext,Microsoft.ML.IDataView)')
  - [BuildAndTrainMultiClassModels(mlContext,TrainingSet)](#M-MMIVR-BiosensorFramework-MachineLearningUtilities-Train-BuildAndTrainMultiClassModels-Microsoft-ML-MLContext,Microsoft-ML-IDataView- 'MMIVR.BiosensorFramework.MachineLearningUtilities.Train.BuildAndTrainMultiClassModels(Microsoft.ML.MLContext,Microsoft.ML.IDataView)')
  - [BuildAndTrainRegressionModels(mlContext,TrainingSet)](#M-MMIVR-BiosensorFramework-MachineLearningUtilities-Train-BuildAndTrainRegressionModels-Microsoft-ML-MLContext,Microsoft-ML-IDataView- 'MMIVR.BiosensorFramework.MachineLearningUtilities.Train.BuildAndTrainRegressionModels(Microsoft.ML.MLContext,Microsoft.ML.IDataView)')
  - [MultiToBin(FeatureSet)](#M-MMIVR-BiosensorFramework-MachineLearningUtilities-Train-MultiToBin-System-Collections-Generic-List{MMIVR-BiosensorFramework-MachineLearningUtilities-ExtractedMultiFeatures}- 'MMIVR.BiosensorFramework.MachineLearningUtilities.Train.MultiToBin(System.Collections.Generic.List{MMIVR.BiosensorFramework.MachineLearningUtilities.ExtractedMultiFeatures})')
  - [MultiToReg(FeatureSet)](#M-MMIVR-BiosensorFramework-MachineLearningUtilities-Train-MultiToReg-System-Collections-Generic-List{MMIVR-BiosensorFramework-MachineLearningUtilities-ExtractedMultiFeatures}- 'MMIVR.BiosensorFramework.MachineLearningUtilities.Train.MultiToReg(System.Collections.Generic.List{MMIVR.BiosensorFramework.MachineLearningUtilities.ExtractedMultiFeatures})')
  - [PrintBinMetrics(metrics)](#M-MMIVR-BiosensorFramework-MachineLearningUtilities-Train-PrintBinMetrics-Microsoft-ML-Data-BinaryClassificationMetrics- 'MMIVR.BiosensorFramework.MachineLearningUtilities.Train.PrintBinMetrics(Microsoft.ML.Data.BinaryClassificationMetrics)')
  - [PrintMultiMetrics(metrics)](#M-MMIVR-BiosensorFramework-MachineLearningUtilities-Train-PrintMultiMetrics-Microsoft-ML-Data-MulticlassClassificationMetrics- 'MMIVR.BiosensorFramework.MachineLearningUtilities.Train.PrintMultiMetrics(Microsoft.ML.Data.MulticlassClassificationMetrics)')
  - [PrintRegMetrics(metrics)](#M-MMIVR-BiosensorFramework-MachineLearningUtilities-Train-PrintRegMetrics-Microsoft-ML-Data-RegressionMetrics- 'MMIVR.BiosensorFramework.MachineLearningUtilities.Train.PrintRegMetrics(Microsoft.ML.Data.RegressionMetrics)')
  - [RunBenchmarks(DirectoryPath,BestRegModel,BestMultiModel,BestBinModel)](#M-MMIVR-BiosensorFramework-MachineLearningUtilities-Train-RunBenchmarks-System-String,Microsoft-ML-ITransformer@,Microsoft-ML-ITransformer@,Microsoft-ML-ITransformer@- 'MMIVR.BiosensorFramework.MachineLearningUtilities.Train.RunBenchmarks(System.String,Microsoft.ML.ITransformer@,Microsoft.ML.ITransformer@,Microsoft.ML.ITransformer@)')
  - [TrimFeatureSet(FeatureSet,LabelsToRemove)](#M-MMIVR-BiosensorFramework-MachineLearningUtilities-Train-TrimFeatureSet-System-Collections-Generic-List{MMIVR-BiosensorFramework-MachineLearningUtilities-ExtractedMultiFeatures},System-Collections-Generic-List{System-Int32}- 'MMIVR.BiosensorFramework.MachineLearningUtilities.Train.TrimFeatureSet(System.Collections.Generic.List{MMIVR.BiosensorFramework.MachineLearningUtilities.ExtractedMultiFeatures},System.Collections.Generic.List{System.Int32})')
- [Type](#T-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-Window-Type 'MMIVR.BiosensorFramework.DataProcessing.ThirdParty.DSPLib.DSP.Window.Type')
- [Utilities](#T-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-Utilities 'MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.Utilities')
  - [AvailableDevices](#F-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-Utilities-AvailableDevices 'MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.Utilities.AvailableDevices')
  - [ClearReadings(E4Object)](#M-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-Utilities-ClearReadings-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient- 'MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.Utilities.ClearReadings(MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient)')
  - [ConnectDevice(E4Object)](#M-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-Utilities-ConnectDevice-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient- 'MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.Utilities.ConnectDevice(MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient)')
  - [ConnectDeviceBTLE(E4Object,Timeout)](#M-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-Utilities-ConnectDeviceBTLE-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient,System-Int32- 'MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.Utilities.ConnectDeviceBTLE(MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient,System.Int32)')
  - [DisconnectDevice(E4Object)](#M-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-Utilities-DisconnectDevice-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient- 'MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.Utilities.DisconnectDevice(MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient)')
  - [DisconnectDeviceBTLE(E4Object)](#M-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-Utilities-DisconnectDeviceBTLE-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient- 'MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.Utilities.DisconnectDeviceBTLE(MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient)')
  - [GrabWindow(E4Object,Filepath)](#M-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-Utilities-GrabWindow-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient,System-String- 'MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.Utilities.GrabWindow(MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient,System.String)')
  - [ListDiscoveredDevices(E4Object)](#M-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-Utilities-ListDiscoveredDevices-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient- 'MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.Utilities.ListDiscoveredDevices(MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient)')
  - [ListDiscoveredDevicesBTLE(E4Object)](#M-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-Utilities-ListDiscoveredDevicesBTLE-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient- 'MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.Utilities.ListDiscoveredDevicesBTLE(MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient)')
  - [SaveReadings(E4Object,Filepath)](#M-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-Utilities-SaveReadings-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient,System-String- 'MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.Utilities.SaveReadings(MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient,System.String)')
  - [StartE4Server(ServerPath,APIKey,IPaddress,Port)](#M-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-Utilities-StartE4Server-System-String,System-String,System-String,System-String- 'MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.Utilities.StartE4Server(System.String,System.String,System.String,System.String)')
  - [StartE4ServerGUI(ServerPath)](#M-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-Utilities-StartE4ServerGUI-System-String- 'MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.Utilities.StartE4ServerGUI(System.String)')
  - [StartStreaming(E4Object)](#M-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-Utilities-StartStreaming-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient- 'MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.Utilities.StartStreaming(MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient)')
  - [StartTimer(Seconds)](#M-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-Utilities-StartTimer-System-Single- 'MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.Utilities.StartTimer(System.Single)')
  - [SubscribeToStream(E4Object,Stream)](#M-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-Utilities-SubscribeToStream-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient,MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-DeviceStreams- 'MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.Utilities.SubscribeToStream(MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient,MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient.DeviceStreams)')
  - [SuspendStreaming(E4Object)](#M-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-Utilities-SuspendStreaming-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient- 'MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.Utilities.SuspendStreaming(MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient)')
  - [UnsubscribeToStream(E4Object,Stream)](#M-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-Utilities-UnsubscribeToStream-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient,MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-DeviceStreams- 'MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.Utilities.UnsubscribeToStream(MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient,MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient.DeviceStreams)')
- [Window](#T-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-Window 'MMIVR.BiosensorFramework.DataProcessing.ThirdParty.DSPLib.DSP.Window')
  - [Coefficients(windowName,points)](#M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-Window-Coefficients-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-Window-Type,System-UInt32- 'MMIVR.BiosensorFramework.DataProcessing.ThirdParty.DSPLib.DSP.Window.Coefficients(MMIVR.BiosensorFramework.DataProcessing.ThirdParty.DSPLib.DSP.Window.Type,System.UInt32)')
- [WindowedDataReady](#T-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-Utilities-WindowedDataReady 'MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.Utilities.WindowedDataReady')

<a name='T-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-Analyze'></a>
## Analyze `type`

##### Namespace

MMIVR.BiosensorFramework.DataProcessing.ThirdParty.DSPLib.DSP

##### Summary

DFT / FFT Output Analysis Functions

<a name='M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-Analyze-FindMaxAmplitude-System-Double[]-'></a>
### FindMaxAmplitude(inData) `method`

##### Summary

Finds the maximum value in an array.

##### Returns

Maximum value of input array

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| inData | [System.Double[]](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double[] 'System.Double[]') |  |

<a name='M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-Analyze-FindMaxFrequency-System-Double[],System-Double[]-'></a>
### FindMaxFrequency(inData,fSpan) `method`

##### Summary

Finds the maximum frequency from the given inData and fSpan arrays.

##### Returns

Maximum frequency from input arrays

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| inData | [System.Double[]](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double[] 'System.Double[]') |  |
| fSpan | [System.Double[]](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double[] 'System.Double[]') |  |

<a name='M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-Analyze-FindMaxPosition-System-Double[]-'></a>
### FindMaxPosition(inData) `method`

##### Summary

Finds the position in the inData array where the maximum value happens.

##### Returns

Position of maximum value in input array

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| inData | [System.Double[]](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double[] 'System.Double[]') |  |

<a name='M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-Analyze-FindMean-System-Double[],System-UInt32,System-UInt32-'></a>
### FindMean(inData,startBin,stopBin) `method`

##### Summary

Finds the mean of the input array.

##### Returns

Mean value of input array between start and stop bins.

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| inData | [System.Double[]](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double[] 'System.Double[]') | = of N data points, 0 based. |
| startBin | [System.UInt32](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.UInt32 'System.UInt32') | = Bin to start the counting at (0 based)."> |
| stopBin | [System.UInt32](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.UInt32 'System.UInt32') | = Bin FROM END to stop counting at (Max = N - 1)."> |

<a name='M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-Analyze-FindRms-System-Double[],System-UInt32,System-UInt32-'></a>
### FindRms(inData,startBin,stopBin) `method`

##### Summary

Find the RMS value of a[].

##### Returns

RMS value of input array between start and stop bins.

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| inData | [System.Double[]](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double[] 'System.Double[]') | = of N data points, 0 based. |
| startBin | [System.UInt32](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.UInt32 'System.UInt32') | = Bin to start the counting at (0 based)."> |
| stopBin | [System.UInt32](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.UInt32 'System.UInt32') | = Bin FROM END to stop counting at (Max = N - 1)."> |

<a name='M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-Analyze-UnwrapPhaseDegrees-System-Double[]-'></a>
### UnwrapPhaseDegrees(inPhaseDeg) `method`

##### Summary

Unwraps the phase so that it is continuous, without jumps.

##### Returns

Continuous Phase data

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| inPhaseDeg | [System.Double[]](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double[] 'System.Double[]') | Array of Phase Data from FT in Degrees |

<a name='M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-Analyze-UnwrapPhaseRadians-System-Double[]-'></a>
### UnwrapPhaseRadians(inPhaseRad) `method`

##### Summary

Unwraps the phase so that it is continuous, without jumps.

##### Returns

Continuous Phase data

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| inPhaseRad | [System.Double[]](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double[] 'System.Double[]') | Array of Phase Data from FT in Radians |

<a name='T-MMIVR-BiosensorFramework-Extensions-ArrayExtensions'></a>
## ArrayExtensions `type`

##### Namespace

MMIVR.BiosensorFramework.Extensions

##### Summary

Extension methods to support the array manipulations needed for feature extraction and machine learning.

<a name='M-MMIVR-BiosensorFramework-Extensions-ArrayExtensions-AllIndexesOf-System-Collections-Generic-List{MMIVR-BiosensorFramework-MachineLearningUtilities-ExtractedMultiFeatures},System-Int32-'></a>
### AllIndexesOf(array,value) `method`

##### Summary

Gets all indices from an ExtractedMultiFeatures class object that match target label.

##### Returns

A List of indices where the feature set matches the target label.

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| array | [System.Collections.Generic.List{MMIVR.BiosensorFramework.MachineLearningUtilities.ExtractedMultiFeatures}](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Collections.Generic.List 'System.Collections.Generic.List{MMIVR.BiosensorFramework.MachineLearningUtilities.ExtractedMultiFeatures}') | A list of labeled feature sets. |
| value | [System.Int32](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Int32 'System.Int32') | The class to search for. |

<a name='M-MMIVR-BiosensorFramework-Extensions-ArrayExtensions-GetSubArray``1-``0[],System-Int32,System-Int32-'></a>
### GetSubArray\`\`1(array,start,end) `method`

##### Summary

Generic extension to split an array from the start index to end index.

##### Returns

The sub array from start to end of original array.

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| array | [\`\`0[]](#T-``0[] '``0[]') | An array of values with more than one value. |
| start | [System.Int32](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Int32 'System.Int32') | The index to start the sub array. |
| end | [System.Int32](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Int32 'System.Int32') | The noninclusive index to end the sub array. |

##### Generic Types

| Name | Description |
| ---- | ----------- |
| T | Standard data types. |

<a name='M-MMIVR-BiosensorFramework-Extensions-ArrayExtensions-Shuffle``1-System-Collections-Generic-IList{``0}-'></a>
### Shuffle\`\`1(list) `method`

##### Summary

Generic implementation of Fisher-Yates shuffling algorithm.

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| list | [System.Collections.Generic.IList{\`\`0}](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Collections.Generic.IList 'System.Collections.Generic.IList{``0}') | List of values to be shuffled. |

##### Generic Types

| Name | Description |
| ---- | ----------- |
| T | Standard data types. |

<a name='M-MMIVR-BiosensorFramework-Extensions-ArrayExtensions-SplitAcc3D-System-Double[],System-Double[]@,System-Double[]@,System-Double[]@-'></a>
### SplitAcc3D(array,x,y,z) `method`

##### Summary

Splits the 3D accelerometer readings into three separate arrays for the input pipeline.

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| array | [System.Double[]](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double[] 'System.Double[]') | Array with 3D accelerometer readings |
| x | [System.Double[]@](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double[]@ 'System.Double[]@') | Array with only x axis readings |
| y | [System.Double[]@](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double[]@ 'System.Double[]@') | Array with only y axis readings |
| z | [System.Double[]@](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double[]@ 'System.Double[]@') | Array with only z axis readings |

<a name='M-MMIVR-BiosensorFramework-Extensions-ArrayExtensions-ToFloat-System-Double[]-'></a>
### ToFloat(array) `method`

##### Summary

Converts a double array to a float array.

##### Returns

The converted double array in float format.

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| array | [System.Double[]](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double[] 'System.Double[]') | Original double array. |

<a name='T-MMIVR-BiosensorFramework-DataProcessing-BloodVolumePulse'></a>
## BloodVolumePulse `type`

##### Namespace

MMIVR.BiosensorFramework.DataProcessing

<a name='M-MMIVR-BiosensorFramework-DataProcessing-BloodVolumePulse-GetHeartFeatures-System-Int32[],System-Double-'></a>
### GetHeartFeatures(Beats,SamplingRate) `method`

##### Summary

Computes the features of the PPG signal.  Returns the heartrate, RR intervals, NN50 metric, pNN50 metric.

##### Returns



##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| Beats | [System.Int32[]](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Int32[] 'System.Int32[]') | Indices of heart beats in PPG signal. |
| SamplingRate | [System.Double](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double 'System.Double') | The sensor sampling rate. |

<a name='T-MMIVR-BiosensorFramework-MachineLearningUtilities-DataImport-CollectedData'></a>
## CollectedData `type`

##### Namespace

MMIVR.BiosensorFramework.MachineLearningUtilities.DataImport

##### Summary

The organization of the imported data.

<a name='T-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-ConvertComplex'></a>
## ConvertComplex `type`

##### Namespace

MMIVR.BiosensorFramework.DataProcessing.ThirdParty.DSPLib.DSP

##### Summary

DFT / FFT Format Conversion Functions.

<a name='M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-ConvertComplex-ToMagnitude-System-Numerics-Complex[]-'></a>
### ToMagnitude(rawFFT) `method`

##### Summary

Convert Complex DFT/FFT Result to: Magnitude Vrms

##### Returns

double[] Magnitude Format (Vrms)

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| rawFFT | [System.Numerics.Complex[]](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Numerics.Complex[] 'System.Numerics.Complex[]') |  |

<a name='M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-ConvertComplex-ToMagnitudeDBV-System-Numerics-Complex[]-'></a>
### ToMagnitudeDBV(rawFFT) `method`

##### Summary

Convert Complex DFT/FFT Result to: Log Magnitude dBV

##### Returns

double[] Magnitude Format (dBV)

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| rawFFT | [System.Numerics.Complex[]](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Numerics.Complex[] 'System.Numerics.Complex[]') | Complex[] input array"> |

<a name='M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-ConvertComplex-ToMagnitudeSquared-System-Numerics-Complex[]-'></a>
### ToMagnitudeSquared(rawFFT) `method`

##### Summary

Convert Complex DFT/FFT Result to: Magnitude Squared V^2 rms

##### Returns

double[] MagSquared Format

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| rawFFT | [System.Numerics.Complex[]](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Numerics.Complex[] 'System.Numerics.Complex[]') |  |

<a name='M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-ConvertComplex-ToPhaseDegrees-System-Numerics-Complex[]-'></a>
### ToPhaseDegrees(rawFFT) `method`

##### Summary

Convert Complex DFT/FFT Result to: Phase in Degrees

##### Returns

double[] Phase (Degrees)

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| rawFFT | [System.Numerics.Complex[]](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Numerics.Complex[] 'System.Numerics.Complex[]') | Complex[] input array"> |

<a name='M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-ConvertComplex-ToPhaseRadians-System-Numerics-Complex[]-'></a>
### ToPhaseRadians(rawFFT) `method`

##### Summary

Convert Complex DFT/FFT Result to: Phase in Radians

##### Returns

double[] Phase (Degrees)

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| rawFFT | [System.Numerics.Complex[]](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Numerics.Complex[] 'System.Numerics.Complex[]') | Complex[] input array"> |

<a name='T-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-ConvertMagnitude'></a>
## ConvertMagnitude `type`

##### Namespace

MMIVR.BiosensorFramework.DataProcessing.ThirdParty.DSPLib.DSP

##### Summary

DFT / FFT Format Conversion Functions

<a name='M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-ConvertMagnitude-ToMagnitudeDBV-System-Double[]-'></a>
### ToMagnitudeDBV(magnitude) `method`

##### Summary

Convert Magnitude FT Result to: Magnitude dBVolts

##### Returns

double[] array

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| magnitude | [System.Double[]](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double[] 'System.Double[]') |  |

<a name='M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-ConvertMagnitude-ToMagnitudeSquared-System-Double[]-'></a>
### ToMagnitudeSquared(magnitude) `method`

##### Summary

Convert Magnitude FT Result to: Magnitude Squared Format

##### Returns



##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| magnitude | [System.Double[]](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double[] 'System.Double[]') |  |

<a name='T-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-ConvertMagnitudeSquared'></a>
## ConvertMagnitudeSquared `type`

##### Namespace

MMIVR.BiosensorFramework.DataProcessing.ThirdParty.DSPLib.DSP

##### Summary

DFT / FFT Format Conversion Functions

<a name='M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-ConvertMagnitudeSquared-ToMagnitude-System-Double[]-'></a>
### ToMagnitude(magSquared) `method`

##### Summary

Convert Magnitude Squared FFT Result to: Magnitude Vrms

##### Returns

double[] array

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| magSquared | [System.Double[]](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double[] 'System.Double[]') |  |

<a name='M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-ConvertMagnitudeSquared-ToMagnitudeDBV-System-Double[]-'></a>
### ToMagnitudeDBV(magSquared) `method`

##### Summary

Convert Magnitude Squared FFT Result to: Magnitude dBVolts

##### Returns

double[] array

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| magSquared | [System.Double[]](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double[] 'System.Double[]') |  |

<a name='T-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DFT'></a>
## DFT `type`

##### Namespace

MMIVR.BiosensorFramework.DataProcessing.ThirdParty.DSPLib

##### Summary

DFT Base Class

<a name='M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DFT-#ctor'></a>
### #ctor() `constructor`

##### Summary

DFT Class

##### Parameters

This constructor has no parameters.

<a name='P-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DFT-IsUsingCached'></a>
### IsUsingCached `property`

##### Summary

Read only Boolean property. True meas the currently defined DFT is using cached memory to speed up calculations.

<a name='M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DFT-Dft-System-Double[]-'></a>
### Dft(timeSeries) `method`

##### Summary

A brute force DFT - Uses Task / Parallel pattern

##### Returns

Complex[] result

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| timeSeries | [System.Double[]](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double[] 'System.Double[]') |  |

<a name='M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DFT-DftCached-System-Double[]-'></a>
### DftCached(timeSeries) `method`

##### Summary

DFT with Pre-calculated Sin/Cos arrays + Task / Parallel pattern.
DFT can only be so big before the computer runs out of memory and has to use
the brute force DFT.

##### Returns

Complex[] result

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| timeSeries | [System.Double[]](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double[] 'System.Double[]') |  |

<a name='M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DFT-Execute-System-Double[]-'></a>
### Execute(timeSeries) `method`

##### Summary

Execute the DFT.

##### Returns

Complex[] FFT Result

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| timeSeries | [System.Double[]](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double[] 'System.Double[]') |  |

<a name='M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DFT-FrequencySpan-System-Double-'></a>
### FrequencySpan(samplingFrequencyHz) `method`

##### Summary

Return the Frequency Array for the currently defined DFT.
Takes into account the total number of points and zero padding points that were defined.

##### Returns



##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| samplingFrequencyHz | [System.Double](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double 'System.Double') |  |

<a name='M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DFT-Initialize-System-UInt32,System-UInt32,System-Boolean-'></a>
### Initialize(inputDataLength,zeroPaddingLength,forceNoCache) `method`

##### Summary

Pre-Initializes the DFT.
Must call first and this anytime the FFT setup changes.

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| inputDataLength | [System.UInt32](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.UInt32 'System.UInt32') |  |
| zeroPaddingLength | [System.UInt32](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.UInt32 'System.UInt32') |  |
| forceNoCache | [System.Boolean](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Boolean 'System.Boolean') | True will force the DFT to not use pre-calculated caching. |

<a name='T-MMIVR-BiosensorFramework-MachineLearningUtilities-DataImport'></a>
## DataImport `type`

##### Namespace

MMIVR.BiosensorFramework.MachineLearningUtilities

##### Summary

Class for importing data from CSV files and txt files.

<a name='M-MMIVR-BiosensorFramework-MachineLearningUtilities-DataImport-ConvertRawToBin-Microsoft-ML-MLContext,System-Collections-Generic-List{System-Tuple{System-Double[],System-Int32}},System-Double-'></a>
### ConvertRawToBin(mlContext,RawFeatures,TrainTestRatio) `method`

##### Summary

Converts raw List of Tuples into a TrainTestData set for Microsoft.ML for binary classification training.

##### Returns

The TrainTestData to be used in training a Microsoft.ML model.

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| mlContext | [Microsoft.ML.MLContext](#T-Microsoft-ML-MLContext 'Microsoft.ML.MLContext') | The Microsoft.ML context. |
| RawFeatures | [System.Collections.Generic.List{System.Tuple{System.Double[],System.Int32}}](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Collections.Generic.List 'System.Collections.Generic.List{System.Tuple{System.Double[],System.Int32}}') | The data and tags loaded in. |
| TrainTestRatio | [System.Double](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double 'System.Double') | The percentage of the data that is train and test. Default is 0.1. |

<a name='M-MMIVR-BiosensorFramework-MachineLearningUtilities-DataImport-ConvertRawToMulti-Microsoft-ML-MLContext,System-Collections-Generic-List{System-Tuple{System-Double[],System-Int32}},System-Double-'></a>
### ConvertRawToMulti(mlContext,RawFeatures,TrainTestRatio) `method`

##### Summary

Converts raw List of Tuples into a TrainTestData set for Microsoft.ML multi-class training.

##### Returns

The TrainTestData to be used in training a Microsoft.ML model.

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| mlContext | [Microsoft.ML.MLContext](#T-Microsoft-ML-MLContext 'Microsoft.ML.MLContext') | The Microsoft.ML context. |
| RawFeatures | [System.Collections.Generic.List{System.Tuple{System.Double[],System.Int32}}](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Collections.Generic.List 'System.Collections.Generic.List{System.Tuple{System.Double[],System.Int32}}') | The data and tags loaded in. |
| TrainTestRatio | [System.Double](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double 'System.Double') | The percentage of the data that is train and test. Default is 0.1. |

<a name='M-MMIVR-BiosensorFramework-MachineLearningUtilities-DataImport-LoadCollectedDataset-System-String,System-String,System-String,System-Int32-'></a>
### LoadCollectedDataset(WesadDirectory,DirectoryPath,SearchPattern,WindowSize) `method`

##### Summary

Parses and transforms the data collected from a biosensor from a txt file.  To facilitate training the Fast Forest classifier used, the WESAD data is loaded in 
first to allow retraining.

##### Returns

List of Tuples with data and tags.

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| WesadDirectory | [System.String](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.String 'System.String') | The top directory of the WESAD dataset. |
| DirectoryPath | [System.String](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.String 'System.String') | The top directory with the data. |
| SearchPattern | [System.String](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.String 'System.String') | The search pattern to find the relevant files in the DirectoryPath. |
| WindowSize | [System.Int32](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Int32 'System.Int32') | The window size, in seconds, for collecting data. Defaults to 5 seconds. |

<a name='M-MMIVR-BiosensorFramework-MachineLearningUtilities-DataImport-LoadDataset-System-String-'></a>
### LoadDataset(DirectoryPath) `method`

##### Summary

Loads in a CSV datasets from the filesystem from the WESAD dataset.

##### Returns

List of Tuples with Subject ID, sensor data, and tags.

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| DirectoryPath | [System.String](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.String 'System.String') | The top directory with the data. |

<a name='M-MMIVR-BiosensorFramework-MachineLearningUtilities-DataImport-LoadFile-System-String,System-Int32-'></a>
### LoadFile(Filepath,WindowSize) `method`

##### Summary

Loads a TXT file of stored readings.

##### Returns

List of Tuples containing data and tags.

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| Filepath | [System.String](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.String 'System.String') | The absolute filepath to the file. |
| WindowSize | [System.Int32](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Int32 'System.Int32') | The window size, in seconds, for collecting data. Defaults to 5 seconds. |

<a name='M-MMIVR-BiosensorFramework-MachineLearningUtilities-DataImport-SchmidtDatasetPipeline-System-String,Microsoft-ML-MLContext,Microsoft-ML-DataOperationsCatalog-TrainTestData@,Microsoft-ML-DataOperationsCatalog-TrainTestData@,Microsoft-ML-DataOperationsCatalog-TrainTestData@,System-Int32,System-Double-'></a>
### SchmidtDatasetPipeline(DirectoryPath,mlContext,MultiClass,BinClass,RegClass,WindowSize,TrainTestRatio) `method`

##### Summary

Pipeline implemented to process the WESAD dataset.
Schmidt, Philip, et al. "Introducing wesad, a multimodal dataset for wearable stress and affect detection." 
    Proceedings of the 20th ACM International Conference on Multimodal Interaction. 2018.

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| DirectoryPath | [System.String](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.String 'System.String') | The top directory with the data. |
| mlContext | [Microsoft.ML.MLContext](#T-Microsoft-ML-MLContext 'Microsoft.ML.MLContext') | The Microsoft.ML context. |
| MultiClass | [Microsoft.ML.DataOperationsCatalog.TrainTestData@](#T-Microsoft-ML-DataOperationsCatalog-TrainTestData@ 'Microsoft.ML.DataOperationsCatalog.TrainTestData@') | The returned TrainTestData for the multi-class models. |
| BinClass | [Microsoft.ML.DataOperationsCatalog.TrainTestData@](#T-Microsoft-ML-DataOperationsCatalog-TrainTestData@ 'Microsoft.ML.DataOperationsCatalog.TrainTestData@') | The returned TrainTestData for the binary models. |
| RegClass | [Microsoft.ML.DataOperationsCatalog.TrainTestData@](#T-Microsoft-ML-DataOperationsCatalog-TrainTestData@ 'Microsoft.ML.DataOperationsCatalog.TrainTestData@') | The returned TrainTestData for the regression models. |
| WindowSize | [System.Int32](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Int32 'System.Int32') | The window size, in seconds, for collecting data. Defaults to 5 seconds. |
| TrainTestRatio | [System.Double](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double 'System.Double') | The percentage of the data that is train and test. Default is 0.1. |

<a name='T-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-DeviceStreams'></a>
## DeviceStreams `type`

##### Namespace

MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient

##### Summary

The applicable device data streams for the E4 server.

<a name='F-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-DeviceStreams-ACC'></a>
### ACC `constants`

##### Summary

Accelerometer data stream.

<a name='F-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-DeviceStreams-BAT'></a>
### BAT `constants`

##### Summary

Battery data stream.

<a name='F-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-DeviceStreams-BVP'></a>
### BVP `constants`

##### Summary

Blood volume pulse data stream.

<a name='F-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-DeviceStreams-GSR'></a>
### GSR `constants`

##### Summary

Galvanic skin response data stream.

<a name='F-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-DeviceStreams-HR'></a>
### HR `constants`

##### Summary

Heart rate data stream.

<a name='F-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-DeviceStreams-IBI'></a>
### IBI `constants`

##### Summary

Inter-beat interval data stream.

<a name='F-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-DeviceStreams-TAG'></a>
### TAG `constants`

##### Summary

Tags data stream. This is the onboard E4 button.

<a name='F-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-DeviceStreams-TMP'></a>
### TMP `constants`

##### Summary

Temperature data stream.

<a name='T-MMIVR-BiosensorFramework-MachineLearningUtilities-ExtractedBinFeatures'></a>
## ExtractedBinFeatures `type`

##### Namespace

MMIVR.BiosensorFramework.MachineLearningUtilities

##### Summary

Class to hold extracted features for training a binary model.

<a name='F-MMIVR-BiosensorFramework-MachineLearningUtilities-ExtractedBinFeatures-Features'></a>
### Features `constants`

##### Summary

Array to store the hand crafted features.

<a name='F-MMIVR-BiosensorFramework-MachineLearningUtilities-ExtractedBinFeatures-Label'></a>
### Label `constants`

##### Summary

The label for the stress features.

<a name='T-MMIVR-BiosensorFramework-MachineLearningUtilities-ExtractedFeatures'></a>
## ExtractedFeatures `type`

##### Namespace

MMIVR.BiosensorFramework.MachineLearningUtilities

##### Summary

Class to hold extracted features for inferencing.

<a name='F-MMIVR-BiosensorFramework-MachineLearningUtilities-ExtractedFeatures-StressFeatures'></a>
### StressFeatures `constants`

##### Summary

Array to store the hand crafted features.

<a name='T-MMIVR-BiosensorFramework-MachineLearningUtilities-ExtractedMultiFeatures'></a>
## ExtractedMultiFeatures `type`

##### Namespace

MMIVR.BiosensorFramework.MachineLearningUtilities

##### Summary

Class to hold extracted features for training a multi-class model.

<a name='F-MMIVR-BiosensorFramework-MachineLearningUtilities-ExtractedMultiFeatures-Result'></a>
### Result `constants`

##### Summary

The label for the stress features.

<a name='F-MMIVR-BiosensorFramework-MachineLearningUtilities-ExtractedMultiFeatures-StressFeatures'></a>
### StressFeatures `constants`

##### Summary

Array to store the hand crafted features.

<a name='T-MMIVR-BiosensorFramework-MachineLearningUtilities-ExtractedRegFeatures'></a>
## ExtractedRegFeatures `type`

##### Namespace

MMIVR.BiosensorFramework.MachineLearningUtilities

##### Summary

Class to hold extracted features for training a regression model.

<a name='F-MMIVR-BiosensorFramework-MachineLearningUtilities-ExtractedRegFeatures-Result'></a>
### Result `constants`

##### Summary

The label for the stress features.

<a name='F-MMIVR-BiosensorFramework-MachineLearningUtilities-ExtractedRegFeatures-StressFeatures'></a>
### StressFeatures `constants`

##### Summary

Array to store the hand crafted features.

<a name='T-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-FFT'></a>
## FFT `type`

##### Namespace

MMIVR.BiosensorFramework.DataProcessing.ThirdParty.DSPLib

##### Summary

FFT Base Class

<a name='M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-FFT-#ctor'></a>
### #ctor() `constructor`

##### Summary

FFT Class

##### Parameters

This constructor has no parameters.

<a name='M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-FFT-BitReverse-System-UInt32,System-UInt32-'></a>
### BitReverse() `method`

##### Parameters

This method has no parameters.

<a name='M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-FFT-Execute-System-Double[]-'></a>
### Execute(timeSeries) `method`

##### Summary

Executes a FFT of the input time series.

##### Returns

Complex[] Spectrum

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| timeSeries | [System.Double[]](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double[] 'System.Double[]') |  |

<a name='M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-FFT-FrequencySpan-System-Double-'></a>
### FrequencySpan(samplingFrequencyHz) `method`

##### Summary

Return the Frequency Array for the currently defined FFT.
Takes into account the total number of points and zero padding points that were defined.

##### Returns



##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| samplingFrequencyHz | [System.Double](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double 'System.Double') |  |

<a name='M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-FFT-Initialize-System-UInt32,System-UInt32-'></a>
### Initialize(inputDataLength,zeroPaddingLength) `method`

##### Summary

Initialize the FFT. Must call first and this anytime the FFT setup changes.

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| inputDataLength | [System.UInt32](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.UInt32 'System.UInt32') |  |
| zeroPaddingLength | [System.UInt32](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.UInt32 'System.UInt32') |  |

<a name='T-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-FIRFilterImplementation'></a>
## FIRFilterImplementation `type`

##### Namespace

MMIVR.BiosensorFramework.DataProcessing.ThirdParty.DSPLib

##### Summary

Source: https://www.codeproject.com/Tips/5070936/Lowpass-Highpass-and-Bandpass-Butterworth-Filters
Modified to process samples in an array, rather than individually.

<a name='T-MMIVR-BiosensorFramework-DataProcessing-FeatureExtraction'></a>
## FeatureExtraction `type`

##### Namespace

MMIVR.BiosensorFramework.DataProcessing

##### Summary

Class with methods for extracted features from time-series data.

<a name='M-MMIVR-BiosensorFramework-DataProcessing-FeatureExtraction-FindLocalMaxima-System-Double[],System-Double,System-Int32-'></a>
### FindLocalMaxima(array,Threshold,Samples) `method`

##### Summary

Finds the local maxima of the signal.

##### Returns



##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| array | [System.Double[]](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double[] 'System.Double[]') | Array of data. |
| Threshold | [System.Double](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double 'System.Double') | The threshold for the local maxima. |
| Samples | [System.Int32](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Int32 'System.Int32') | The number of samples required to classify a maxima. |

<a name='M-MMIVR-BiosensorFramework-DataProcessing-FeatureExtraction-PSI-System-Collections-Generic-List{System-Double},System-Double@,System-Double@,System-Int32-'></a>
### PSI(Accelerations,AverageIndex,MaxIndex,WindowWidth) `method`

##### Summary

Calculates the physical stillness index, specified in the paper referenced below.  
Chang, Kang-Ming, et al. "A wireless accelerometer-based body posture stability detection system and its application for 
    meditation practitioners." Sensors 12.12 (2012): 17620-17632.
Calculates a metric of physical stillness from acceleration readings. 
Uses a sliding window, root mean square filter to get the absolute magnitude of the overall acceleration in the readings.
This method works with the absolute magnitude of accelerations and not three-axis readings, use PSIVector for those purposes.

##### Returns

The output filtered list of RMS values, can be discarded in favor of 'out' variables

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| Accelerations | [System.Collections.Generic.List{System.Double}](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Collections.Generic.List 'System.Collections.Generic.List{System.Double}') | List of acceleration readings |
| AverageIndex | [System.Double@](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double@ 'System.Double@') | Output PSI, the average of the filtered output |
| MaxIndex | [System.Double@](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double@ 'System.Double@') | Output PSI, the max value of the filtered output |
| WindowWidth | [System.Int32](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Int32 'System.Int32') | The number of values to incorporate in the window on both sides |

<a name='M-MMIVR-BiosensorFramework-DataProcessing-FeatureExtraction-PSIVector-System-Collections-Generic-List{System-Collections-Generic-List{System-Double}},System-Double@,System-Double@,System-Int32-'></a>
### PSIVector(Accelerations,AverageIndex,MaxIndex,WindowWidth) `method`

##### Summary

Performs the calculation of the physical stillness index using a 3D acceleration vector

##### Returns

The output filtered list of RMS values, can be discarded in favor of 'out' variables

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| Accelerations | [System.Collections.Generic.List{System.Collections.Generic.List{System.Double}}](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Collections.Generic.List 'System.Collections.Generic.List{System.Collections.Generic.List{System.Double}}') | List of acceleration readings |
| AverageIndex | [System.Double@](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double@ 'System.Double@') | Output PSI, the average of the filtered output |
| MaxIndex | [System.Double@](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double@ 'System.Double@') | Output PSI, the max value of the filtered output |
| WindowWidth | [System.Int32](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Int32 'System.Int32') | The number of values to incorporate in the window on both sides |

<a name='M-MMIVR-BiosensorFramework-DataProcessing-FeatureExtraction-RemoveBias-System-Numerics-Vector3,System-Double,System-Double,System-Double-'></a>
### RemoveBias(AccelerationVector,BiasX,BiasY,BiasZ) `method`

##### Summary

Removes mean from each axis of the 3D vector

##### Returns

Vector3 with subtracted bias.

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| AccelerationVector | [System.Numerics.Vector3](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Numerics.Vector3 'System.Numerics.Vector3') | 3D acceleration vector |
| BiasX | [System.Double](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double 'System.Double') | Mean of X axis |
| BiasY | [System.Double](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double 'System.Double') | Mean of Y axis |
| BiasZ | [System.Double](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double 'System.Double') | Mean of Z axis |

<a name='M-MMIVR-BiosensorFramework-DataProcessing-FeatureExtraction-RootMeanSquare-System-Double,System-Double-'></a>
### RootMeanSquare(FirstTerm,Summation) `method`

##### Summary

Calculates the root mean square of the input variables

##### Returns

Root mean square

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| FirstTerm | [System.Double](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double 'System.Double') | Calculated first term of form - 1/2n |
| Summation | [System.Double](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double 'System.Double') | The square sum of the input values |

<a name='M-MMIVR-BiosensorFramework-DataProcessing-FeatureExtraction-SignalAbsolute-System-Double-'></a>
### SignalAbsolute(AUC) `method`

##### Summary

Returns the absolute value of the area under the curve.

##### Returns

The absolute value of the AUC.

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| AUC | [System.Double](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double 'System.Double') | The area under the curve of the signal. |

<a name='M-MMIVR-BiosensorFramework-DataProcessing-FeatureExtraction-SignalCorrelation-System-Double[],System-Double[]-'></a>
### SignalCorrelation(Signal,TimeLabels) `method`

##### Summary

Computes the Pearson correlation coefficient of the signal.

##### Returns

The correlation coefficient of the data.

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| Signal | [System.Double[]](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double[] 'System.Double[]') | Array of data. |
| TimeLabels | [System.Double[]](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double[] 'System.Double[]') | The time labels for each data point. |

<a name='M-MMIVR-BiosensorFramework-DataProcessing-FeatureExtraction-SignalDynamicRange-System-Double[]-'></a>
### SignalDynamicRange(Signal) `method`

##### Summary

Returns the dynamic range of the signal (max / min).

##### Returns

The dynamic range of the signal.

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| Signal | [System.Double[]](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double[] 'System.Double[]') | Array of data. |

<a name='M-MMIVR-BiosensorFramework-DataProcessing-FeatureExtraction-SignalFreqBandEnergies-System-Double[],System-Collections-Generic-List{System-Tuple{System-Double,System-Double}},System-Double,System-Int32-'></a>
### SignalFreqBandEnergies(Signal,Bands,SamplingFrequency,Order) `method`

##### Summary

Computes the frequency band energies of the signal.

##### Returns

List of frequency band energies.

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| Signal | [System.Double[]](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double[] 'System.Double[]') | Array of data. |
| Bands | [System.Collections.Generic.List{System.Tuple{System.Double,System.Double}}](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Collections.Generic.List 'System.Collections.Generic.List{System.Tuple{System.Double,System.Double}}') | The frequency bands to compute. |
| SamplingFrequency | [System.Double](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double 'System.Double') | The sampling rate of the signal. |
| Order | [System.Int32](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Int32 'System.Int32') | The filter order. Defaults to 5. |

<a name='M-MMIVR-BiosensorFramework-DataProcessing-FeatureExtraction-SignalFreqNormalize-System-Double[]-'></a>
### SignalFreqNormalize(Freq) `method`

##### Summary

Normalizes the frequency components.

##### Returns

The normalized frequency components.

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| Freq | [System.Double[]](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double[] 'System.Double[]') | The signal band energies. |

<a name='M-MMIVR-BiosensorFramework-DataProcessing-FeatureExtraction-SignalFreqRatio-System-Double,System-Double-'></a>
### SignalFreqRatio(LowFreq,HighFreq) `method`

##### Summary

Returns the frequency ratio of the signal.

##### Returns

Returns the ratio of the freqencies (low / high).

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| LowFreq | [System.Double](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double 'System.Double') | The low frequency value. |
| HighFreq | [System.Double](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double 'System.Double') | The high frequency value. |

<a name='M-MMIVR-BiosensorFramework-DataProcessing-FeatureExtraction-SignalFreqSummation-System-Collections-Generic-List{System-Double[]}-'></a>
### SignalFreqSummation(Freq) `method`

##### Summary

Sums the frequency components of the signal band energies.

##### Returns

List of doubles with the summed signal band energies.

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| Freq | [System.Collections.Generic.List{System.Double[]}](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Collections.Generic.List 'System.Collections.Generic.List{System.Double[]}') | The signal band energies. |

<a name='M-MMIVR-BiosensorFramework-DataProcessing-FeatureExtraction-SignalIntegral-System-Double[]-'></a>
### SignalIntegral(Signal) `method`

##### Summary

Computes the integral of the signal using Simpson's Rule.

##### Returns

The integral of the data.

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| Signal | [System.Double[]](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double[] 'System.Double[]') | Array of data. |

<a name='M-MMIVR-BiosensorFramework-DataProcessing-FeatureExtraction-SignalMean-System-Collections-Generic-List{System-Single}-'></a>
### SignalMean(Signal) `method`

##### Summary

Returns the mean of the data.

##### Returns

Mean of the signal.

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| Signal | [System.Collections.Generic.List{System.Single}](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Collections.Generic.List 'System.Collections.Generic.List{System.Single}') | List of data. |

<a name='M-MMIVR-BiosensorFramework-DataProcessing-FeatureExtraction-SignalMinMax-System-Collections-Generic-List{System-Single}-'></a>
### SignalMinMax(Signal) `method`

##### Summary

Gets the signal min and max from the array of data.

##### Returns

Tuple with max and min.

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| Signal | [System.Collections.Generic.List{System.Single}](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Collections.Generic.List 'System.Collections.Generic.List{System.Single}') | List of data. |

<a name='M-MMIVR-BiosensorFramework-DataProcessing-FeatureExtraction-SignalPeakFrequency-System-Double[],System-Double-'></a>
### SignalPeakFrequency(Signal,SamplingRate) `method`

##### Summary

Computes the peak frequency of the signal.

##### Returns

Peak frequency component of the signal.

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| Signal | [System.Double[]](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double[] 'System.Double[]') | Array of data. |
| SamplingRate | [System.Double](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double 'System.Double') | The sampling rate of the signal. |

<a name='M-MMIVR-BiosensorFramework-DataProcessing-FeatureExtraction-SignalPercentile-System-Double[],System-Int32-'></a>
### SignalPercentile(Signal,percentile) `method`

##### Summary

Calculates the signal p-percentile.

##### Returns

The signal percentile.

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| Signal | [System.Double[]](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double[] 'System.Double[]') | Array of data. |
| percentile | [System.Int32](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Int32 'System.Int32') | The percentile to calculate. |

<a name='M-MMIVR-BiosensorFramework-DataProcessing-FeatureExtraction-SignalRMS-System-Double[]-'></a>
### SignalRMS(Signal) `method`

##### Summary

Computes the root mean square of the signal.

##### Returns

The root mean square of the signal.

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| Signal | [System.Double[]](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double[] 'System.Double[]') | Array of data. |

<a name='M-MMIVR-BiosensorFramework-DataProcessing-FeatureExtraction-SignalRelativePower-System-Collections-Generic-List{System-Double[]},System-Double,System-Double-'></a>
### SignalRelativePower(Energies,Resolution,Epsilon) `method`

##### Summary

Computes the signal relative power from the frequency band energies.

##### Returns

List of doubles for the signal relative power.

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| Energies | [System.Collections.Generic.List{System.Double[]}](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Collections.Generic.List 'System.Collections.Generic.List{System.Double[]}') | The frequency band energies. |
| Resolution | [System.Double](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double 'System.Double') | The resolution of the relative power. |
| Epsilon | [System.Double](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double 'System.Double') | Default 1e-7. |

<a name='M-MMIVR-BiosensorFramework-DataProcessing-FeatureExtraction-SignalSlope-System-Double[],System-Double[]-'></a>
### SignalSlope(Signal,TimeLabels) `method`

##### Summary

The slope of the data.

##### Returns

The slope of the signal.

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| Signal | [System.Double[]](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double[] 'System.Double[]') | Array of data. |
| TimeLabels | [System.Double[]](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double[] 'System.Double[]') | Array of time labels for each data point. |

<a name='M-MMIVR-BiosensorFramework-DataProcessing-FeatureExtraction-SignalStandardDeviation-System-Collections-Generic-List{System-Single}-'></a>
### SignalStandardDeviation(Signal) `method`

##### Summary

Returns the standard deviation of the data.

##### Returns

Standard deviation of the signal.

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| Signal | [System.Collections.Generic.List{System.Single}](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Collections.Generic.List 'System.Collections.Generic.List{System.Single}') | List of data. |

<a name='M-MMIVR-BiosensorFramework-DataProcessing-FeatureExtraction-SquareSummation-System-Double-'></a>
### SquareSummation(Value) `method`

##### Summary

Performs square of input values

##### Returns



##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| Value | [System.Double](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double 'System.Double') |  |

<a name='M-MMIVR-BiosensorFramework-DataProcessing-FeatureExtraction-SquareSummation-System-Numerics-Vector3-'></a>
### SquareSummation(Value) `method`

##### Summary

Performs square summation of Vector3 variables

##### Returns

The squared Vector3.

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| Value | [System.Numerics.Vector3](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Numerics.Vector3 'System.Numerics.Vector3') | Vector3 with x, y, and z values |

<a name='T-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-Generate'></a>
## Generate `type`

##### Namespace

MMIVR.BiosensorFramework.DataProcessing.ThirdParty.DSPLib.DSP

<a name='M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-Generate-LinSpace-System-Double,System-Double,System-UInt32-'></a>
### LinSpace(startVal,stopVal,points) `method`

##### Summary

Generate linearly spaced array. Like the Octave function of the same name.
EX: DSP.Generate.LinSpace(1, 10, 10) -> Returns array: 1, 2, 3, 4....10.

##### Returns

double[] array

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| startVal | [System.Double](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double 'System.Double') | Any value |
| stopVal | [System.Double](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double 'System.Double') | Any value > startVal |
| points | [System.UInt32](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.UInt32 'System.UInt32') | Number of points to generate |

<a name='M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-Generate-NoisePsd-System-Double,System-Double,System-UInt32-'></a>
### NoisePsd(amplitudePsd (Vrms / rt-Hz)(Vrms/rt-Hz),samplingFrequencyHz,points) `method`

##### Summary

Generates a normal distribution noise signal of the specified power spectral density (Vrms / rt-Hz).

##### Returns

double[] array

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| amplitudePsd (Vrms / rt-Hz)(Vrms/rt-Hz) | [System.Double](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double 'System.Double') |  |
| samplingFrequencyHz | [System.Double](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double 'System.Double') |  |
| points | [System.UInt32](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.UInt32 'System.UInt32') |  |

<a name='M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-Generate-NoiseRms-System-Double,System-UInt32,System-Double-'></a>
### NoiseRms(amplitudeVrms,points,dcV) `method`

##### Summary

Generates a normal distribution noise signal of the specified Volts RMS.

##### Returns

double[] array

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| amplitudeVrms | [System.Double](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double 'System.Double') |  |
| points | [System.UInt32](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.UInt32 'System.UInt32') |  |
| dcV | [System.Double](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double 'System.Double') |  |

<a name='M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-Generate-ToneCycles-System-Double,System-Double,System-UInt32,System-Double,System-Double-'></a>
### ToneCycles(amplitudeVrms,cycles,points,dcV,phaseDeg) `method`

##### Summary

Generates a Sine Wave Tone using Number of Cycles Terms.

##### Returns

double[] array

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| amplitudeVrms | [System.Double](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double 'System.Double') |  |
| cycles | [System.Double](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double 'System.Double') |  |
| points | [System.UInt32](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.UInt32 'System.UInt32') |  |
| dcV | [System.Double](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double 'System.Double') | [Optional] DC Voltage offset |
| phaseDeg | [System.Double](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double 'System.Double') | [Optional] Phase of signal in degrees |

<a name='M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-Generate-ToneSampling-System-Double,System-Double,System-Double,System-UInt32,System-Double,System-Double-'></a>
### ToneSampling(amplitudeVrms,frequencyHz,samplingFrequencyHz,points,dcV,phaseDeg) `method`

##### Summary

Generates a Sine Wave Tone using Sampling Terms.

##### Returns

double[] array

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| amplitudeVrms | [System.Double](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double 'System.Double') |  |
| frequencyHz | [System.Double](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double 'System.Double') |  |
| samplingFrequencyHz | [System.Double](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double 'System.Double') |  |
| points | [System.UInt32](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.UInt32 'System.UInt32') |  |
| dcV | [System.Double](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double 'System.Double') | [Optional] DC Voltage offset |
| phaseDeg | [System.Double](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double 'System.Double') | [Optional] Phase of signal in degrees |

<a name='T-MMIVR-BiosensorFramework-DataProcessing-HealeyStressDetection'></a>
## HealeyStressDetection `type`

##### Namespace

MMIVR.BiosensorFramework.DataProcessing

<a name='M-MMIVR-BiosensorFramework-DataProcessing-HealeyStressDetection-ConstructPeaks-System-Collections-Generic-List{System-Collections-Generic-List{System-Int32}},System-Tuple{System-Int32[],System-Int32[]}-'></a>
### ConstructPeaks(Peaks,ZeroCrossings) `method`

##### Summary

Takes the zero crossings and peaks and reconstructs the index range for the peaks.

##### Returns

List of peak indices.

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| Peaks | [System.Collections.Generic.List{System.Collections.Generic.List{System.Int32}}](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Collections.Generic.List 'System.Collections.Generic.List{System.Collections.Generic.List{System.Int32}}') | List of List of indices that represent the indices above threshold |
| ZeroCrossings | [System.Tuple{System.Int32[],System.Int32[]}](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Tuple 'System.Tuple{System.Int32[],System.Int32[]}') | List of indices where zero crossings happen |

<a name='M-MMIVR-BiosensorFramework-DataProcessing-HealeyStressDetection-FindAbove-System-Double[],System-Double,System-Double-'></a>
### FindAbove(array,Threshold,SamplingRate) `method`

##### Summary

Returns indices above a threshold, removes zero elements.  Intended for time-series data collection.

##### Returns

List of List of indices above the threshold.

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| array | [System.Double[]](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double[] 'System.Double[]') | An array of data. |
| Threshold | [System.Double](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double 'System.Double') | The threshold for inclusion in return data. |
| SamplingRate | [System.Double](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double 'System.Double') | The sampling rate of the data. |

<a name='M-MMIVR-BiosensorFramework-DataProcessing-HealeyStressDetection-FindNearestMinus-System-Int32,System-Int32[]-'></a>
### FindNearestMinus(PeakStart,ZeroCrossings) `method`

##### Summary

Finds the nearest index for zero crossings that go from above zero to below zero

##### Returns

int, the index in ZeroCrossings that is closest to PeakStart

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| PeakStart | [System.Int32](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Int32 'System.Int32') | Peak starting index |
| ZeroCrossings | [System.Int32[]](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Int32[] 'System.Int32[]') | Array of indices that represent the zero crossing points |

<a name='M-MMIVR-BiosensorFramework-DataProcessing-HealeyStressDetection-FindNearestPlus-System-Int32,System-Int32[]-'></a>
### FindNearestPlus(PeakStart,ZeroCrossings) `method`

##### Summary

Finds the nearest index for zero crossings that go from below zero to above zero

##### Returns

int, the index in ZeroCrossings that is closest to PeakStart

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| PeakStart | [System.Int32](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Int32 'System.Int32') | Peak starting index |
| ZeroCrossings | [System.Int32[]](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Int32[] 'System.Int32[]') | Array of indices that represent the zero crossing points |

<a name='M-MMIVR-BiosensorFramework-DataProcessing-HealeyStressDetection-FindZeroCrossings-System-Double[]-'></a>
### FindZeroCrossings(array) `method`

##### Summary

Find zero crossings for two cases:
    1. Goes from below zero to above zero
    2. Goes from above zero to below zero

##### Returns

A tuple of the form: Item1 is case one and Item2 is case twos

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| array | [System.Double[]](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double[] 'System.Double[]') | Array to find the indices of zero crossings in |

<a name='M-MMIVR-BiosensorFramework-DataProcessing-HealeyStressDetection-FiniteDifference-System-Double[]@-'></a>
### FiniteDifference(Signal) `method`

##### Summary

Calculates the differences between adjacent elements, equating to the approximate first derivative
Algorithm: Y = [X(1)-X(0), X(2)-X(1), ..., X(m + 1)-X(m)]
Source: https://www.mathworks.com/help/matlab/ref/diff.html?s_tid=srchtitle

##### Returns



##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| Signal | [System.Double[]@](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double[]@ 'System.Double[]@') |  |

<a name='M-MMIVR-BiosensorFramework-DataProcessing-HealeyStressDetection-GetPeakFeatures-System-Double[],System-Double,System-Double,System-Double-'></a>
### GetPeakFeatures(Signal,SamplingRate,CutoffFreq,Threshold) `method`

##### Summary

Uses the approximate first derivative of the GSR signal to find the magnitude, duration, and amplitude of the peaks in the 
signal.

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| Signal | [System.Double[]](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double[] 'System.Double[]') |  |
| SamplingRate | [System.Double](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double 'System.Double') |  |
| CutoffFreq | [System.Double](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double 'System.Double') |  |
| Threshold | [System.Double](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double 'System.Double') |  |

<a name='M-MMIVR-BiosensorFramework-DataProcessing-HealeyStressDetection-MergeOverlappingPeaks-System-Collections-Generic-List{System-Int32[]}-'></a>
### MergeOverlappingPeaks(Peaks) `method`

##### Summary

If reconstructed peaks are overlapping in their indices, which can happen if a peak does not zero cross before the next peak,
then remove the copies.

##### Returns

List with any overlapping peaks removed

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| Peaks | [System.Collections.Generic.List{System.Int32[]}](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Collections.Generic.List 'System.Collections.Generic.List{System.Int32[]}') | List of indices for the reconstructed peaks |

<a name='M-MMIVR-BiosensorFramework-DataProcessing-HealeyStressDetection-SumFeatures-System-Collections-Generic-List{System-Double[]}-'></a>
### SumFeatures(PeakFeatures) `method`

##### Summary



##### Returns



##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| PeakFeatures | [System.Collections.Generic.List{System.Double[]}](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Collections.Generic.List 'System.Collections.Generic.List{System.Double[]}') | Output of the GetPeakFeatures function |

<a name='T-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-Math'></a>
## Math `type`

##### Namespace

MMIVR.BiosensorFramework.DataProcessing.ThirdParty.DSPLib.DSP

##### Summary

Double[] Array Math Operations (All Static)

<a name='M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-Math-Add-System-Double[],System-Double[]-'></a>
### Add() `method`

##### Summary

result[] = a[] + b[]

##### Parameters

This method has no parameters.

<a name='M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-Math-Add-System-Double[],System-Double-'></a>
### Add() `method`

##### Summary

result[] = a[] + b

##### Parameters

This method has no parameters.

<a name='M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-Math-Divide-System-Double[],System-Double[]-'></a>
### Divide() `method`

##### Summary

result[] = a[] / b[]

##### Parameters

This method has no parameters.

<a name='M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-Math-Divide-System-Double[],System-Double-'></a>
### Divide() `method`

##### Summary

result[] = a[] / b

##### Parameters

This method has no parameters.

<a name='M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-Math-Log10-System-Double[]-'></a>
### Log10() `method`

##### Summary

Log10 a[].

##### Parameters

This method has no parameters.

<a name='M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-Math-Multiply-System-Double[],System-Double[]-'></a>
### Multiply() `method`

##### Summary

result[] = a[] * b[]

##### Parameters

This method has no parameters.

<a name='M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-Math-Multiply-System-Double[],System-Double-'></a>
### Multiply() `method`

##### Summary

result[] = a[] * b

##### Parameters

This method has no parameters.

<a name='M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-Math-RemoveMean-System-Double[]-'></a>
### RemoveMean() `method`

##### Summary

Removes mean value from a[].

##### Parameters

This method has no parameters.

<a name='M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-Math-Sqrt-System-Double[]-'></a>
### Sqrt() `method`

##### Summary

Square root of a[].

##### Parameters

This method has no parameters.

<a name='M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-Math-Square-System-Double[]-'></a>
### Square() `method`

##### Summary

Squares a[].

##### Parameters

This method has no parameters.

<a name='M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-Math-Subtract-System-Double[],System-Double[]-'></a>
### Subtract() `method`

##### Summary

result[] = a[] - b[]

##### Parameters

This method has no parameters.

<a name='M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-Math-Subtract-System-Double[],System-Double-'></a>
### Subtract() `method`

##### Summary

result[] = a[] - b

##### Parameters

This method has no parameters.

<a name='T-MMIVR-BiosensorFramework-MachineLearningUtilities-Predict'></a>
## Predict `type`

##### Namespace

MMIVR.BiosensorFramework.MachineLearningUtilities

##### Summary

Class for performing predictions on Microsoft.ML models.

<a name='M-MMIVR-BiosensorFramework-MachineLearningUtilities-Predict-MakeBinPrediction-Microsoft-ML-MLContext,MMIVR-BiosensorFramework-MachineLearningUtilities-ExtractedBinFeatures,Microsoft-ML-ITransformer-'></a>
### MakeBinPrediction(mlContext,LiveData,Model) `method`

##### Summary

Performs a prediction on a dataset for binary classification models.

##### Returns



##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| mlContext | [Microsoft.ML.MLContext](#T-Microsoft-ML-MLContext 'Microsoft.ML.MLContext') | Microsoft ML context for operations to be performed in. |
| LiveData | [MMIVR.BiosensorFramework.MachineLearningUtilities.ExtractedBinFeatures](#T-MMIVR-BiosensorFramework-MachineLearningUtilities-ExtractedBinFeatures 'MMIVR.BiosensorFramework.MachineLearningUtilities.ExtractedBinFeatures') | The extracted features packaged in the ExtractedMultiFeatures class. |
| Model | [Microsoft.ML.ITransformer](#T-Microsoft-ML-ITransformer 'Microsoft.ML.ITransformer') | The model to perform inferencing with. |

<a name='M-MMIVR-BiosensorFramework-MachineLearningUtilities-Predict-MakeMultiPrediction-Microsoft-ML-MLContext,MMIVR-BiosensorFramework-MachineLearningUtilities-ExtractedMultiFeatures,Microsoft-ML-ITransformer-'></a>
### MakeMultiPrediction(mlContext,LiveData,Model) `method`

##### Summary

Performs a prediction on a dataset for multi-class models.

##### Returns



##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| mlContext | [Microsoft.ML.MLContext](#T-Microsoft-ML-MLContext 'Microsoft.ML.MLContext') | Microsoft ML context for operations to be performed in. |
| LiveData | [MMIVR.BiosensorFramework.MachineLearningUtilities.ExtractedMultiFeatures](#T-MMIVR-BiosensorFramework-MachineLearningUtilities-ExtractedMultiFeatures 'MMIVR.BiosensorFramework.MachineLearningUtilities.ExtractedMultiFeatures') | The extracted features packaged in the ExtractedMultiFeatures class. |
| Model | [Microsoft.ML.ITransformer](#T-Microsoft-ML-ITransformer 'Microsoft.ML.ITransformer') | The model to perform inferencing with. |

<a name='M-MMIVR-BiosensorFramework-MachineLearningUtilities-Predict-MakeRegPrediction-Microsoft-ML-MLContext,Microsoft-ML-IDataView,Microsoft-ML-ITransformer-'></a>
### MakeRegPrediction(mlContext,LiveData,Model) `method`

##### Summary

Performs a prediction on a dataset for regression classification models.

##### Returns



##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| mlContext | [Microsoft.ML.MLContext](#T-Microsoft-ML-MLContext 'Microsoft.ML.MLContext') | Microsoft ML context for operations to be performed in. |
| LiveData | [Microsoft.ML.IDataView](#T-Microsoft-ML-IDataView 'Microsoft.ML.IDataView') | The extracted features packaged in the ExtractedMultiFeatures class. |
| Model | [Microsoft.ML.ITransformer](#T-Microsoft-ML-ITransformer 'Microsoft.ML.ITransformer') | The model to perform inferencing with. |

<a name='M-MMIVR-BiosensorFramework-MachineLearningUtilities-Predict-PredictWindow-Microsoft-ML-MLContext,Microsoft-ML-ITransformer,System-Collections-Generic-List{System-Collections-Generic-List{System-Double}}-'></a>
### PredictWindow(mlContext,Model,WindowReadings) `method`

##### Summary

Takes the readings from the windowed data, extracts the features, and runs it through a prediction pipeline.

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| mlContext | [Microsoft.ML.MLContext](#T-Microsoft-ML-MLContext 'Microsoft.ML.MLContext') | Microsoft ML context for operations to be performed in. |
| Model | [Microsoft.ML.ITransformer](#T-Microsoft-ML-ITransformer 'Microsoft.ML.ITransformer') | The loaded model for operations to be performed on. |
| WindowReadings | [System.Collections.Generic.List{System.Collections.Generic.List{System.Double}}](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Collections.Generic.List 'System.Collections.Generic.List{System.Collections.Generic.List{System.Double}}') | Packaged List of List of sensor readings. |

<a name='T-MMIVR-BiosensorFramework-MachineLearningUtilities-PredictionBinResult'></a>
## PredictionBinResult `type`

##### Namespace

MMIVR.BiosensorFramework.MachineLearningUtilities

##### Summary

Class to hold binary prediction result.

<a name='F-MMIVR-BiosensorFramework-MachineLearningUtilities-PredictionBinResult-Prediction'></a>
### Prediction `constants`

##### Summary

The label of the inferenced data.

<a name='F-MMIVR-BiosensorFramework-MachineLearningUtilities-PredictionBinResult-Probability'></a>
### Probability `constants`

##### Summary

The probability of the predicted label.

<a name='F-MMIVR-BiosensorFramework-MachineLearningUtilities-PredictionBinResult-Score'></a>
### Score `constants`

##### Summary

The score of the predicted label.

<a name='T-MMIVR-BiosensorFramework-MachineLearningUtilities-PredictionMultiResult'></a>
## PredictionMultiResult `type`

##### Namespace

MMIVR.BiosensorFramework.MachineLearningUtilities

##### Summary

Class to hold multi-class prediction result.

<a name='F-MMIVR-BiosensorFramework-MachineLearningUtilities-PredictionMultiResult-Result'></a>
### Result `constants`

##### Summary

The label of the inferenced data.

<a name='T-MMIVR-BiosensorFramework-MachineLearningUtilities-PredictionRegResult'></a>
## PredictionRegResult `type`

##### Namespace

MMIVR.BiosensorFramework.MachineLearningUtilities

##### Summary

Class to hold regression prediction result.

<a name='F-MMIVR-BiosensorFramework-MachineLearningUtilities-PredictionRegResult-Result'></a>
### Result `constants`

##### Summary

The label of the inferenced data.

<a name='T-MMIVR-BiosensorFramework-DataProcessing-ProcessFFT'></a>
## ProcessFFT `type`

##### Namespace

MMIVR.BiosensorFramework.DataProcessing

##### Summary

Class to process a signal through FFT and return the frequency and magnitude spectrums.

<a name='M-MMIVR-BiosensorFramework-DataProcessing-ProcessFFT-ProcessSignal-System-Double[]@,System-Double,System-Double[]@,System-Double[]@-'></a>
### ProcessSignal(Signal,SamplingRate,FreqSpan,MagSpectrum) `method`

##### Summary

Processes the input signal through the FFT and returns the frequency span and magnitude spectrum.

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| Signal | [System.Double[]@](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double[]@ 'System.Double[]@') | An array containing the original sensor signals. |
| SamplingRate | [System.Double](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double 'System.Double') | The sampling rate of the input signal. |
| FreqSpan | [System.Double[]@](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double[]@ 'System.Double[]@') | The frequency span of the input signal. |
| MagSpectrum | [System.Double[]@](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double[]@ 'System.Double[]@') | The magnitude spectrum of the input signal. |

<a name='T-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-Window-ScaleFactor'></a>
## ScaleFactor `type`

##### Namespace

MMIVR.BiosensorFramework.DataProcessing.ThirdParty.DSPLib.DSP.Window

<a name='M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-Window-ScaleFactor-NENBW-System-Double[]-'></a>
### NENBW(windowCoefficients) `method`

##### Summary

Calculate Normalized, Equivalent Noise BandWidth from window coefficient array.

##### Returns

double NENBW

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| windowCoefficients | [System.Double[]](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double[] 'System.Double[]') |  |

<a name='M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-Window-ScaleFactor-Noise-System-Double[],System-Double-'></a>
### Noise(windowCoefficients,samplingFrequencyHz) `method`

##### Summary

Calculate Noise scale factor from window coefficient array.
 Takes into account the bin width in Hz for the final result also.
 Designed to be applied to the "Magnitude" result.

##### Returns

double scaleFactor

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| windowCoefficients | [System.Double[]](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double[] 'System.Double[]') |  |
| samplingFrequencyHz | [System.Double](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double 'System.Double') |  |

<a name='M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-Window-ScaleFactor-Signal-System-Double[]-'></a>
### Signal(windowCoefficients) `method`

##### Summary

Calculate Signal scale factor from window coefficient array.
Designed to be applied to the "Magnitude" result.

##### Returns

double scaleFactor

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| windowCoefficients | [System.Double[]](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double[] 'System.Double[]') |  |

<a name='T-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient'></a>
## ServerClient `type`

##### Namespace

MMIVR.BiosensorFramework.Biosensors.EmpaticaE4

##### Summary

Class to instantiate a client to communicate with TCP server.

<a name='F-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-ACCReadings3D'></a>
### ACCReadings3D `constants`

##### Summary

List for 3D accelerometer readings.

<a name='F-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-ACCReadingsX'></a>
### ACCReadingsX `constants`

##### Summary

List for X accelerometer readings.

<a name='F-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-ACCReadingsY'></a>
### ACCReadingsY `constants`

##### Summary

List for Y accelerometer readings.

<a name='F-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-ACCReadingsZ'></a>
### ACCReadingsZ `constants`

##### Summary

List for Z accelerometer readings.

<a name='F-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-ACCTimestamps'></a>
### ACCTimestamps `constants`

##### Summary

List for accelerometer timestamps.

<a name='F-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-BATReadings'></a>
### BATReadings `constants`

##### Summary

List for battery readings.

<a name='F-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-BATTimestamps'></a>
### BATTimestamps `constants`

##### Summary

List for battery timestamps.

<a name='F-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-BVPReadings'></a>
### BVPReadings `constants`

##### Summary

List for blood volume pulse readings.

<a name='F-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-BVPTimestamps'></a>
### BVPTimestamps `constants`

##### Summary

List for blood volume pulse timestamps.

<a name='F-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-Buffer'></a>
### Buffer `constants`

##### Summary

Receive buffer.

<a name='F-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-BufferSize'></a>
### BufferSize `constants`

##### Summary

Size of receive buffer.

<a name='F-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-ConnectDone'></a>
### ConnectDone `constants`

##### Summary

ManualResetEvent instances signal completion.

<a name='F-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-DeviceName'></a>
### DeviceName `constants`

##### Summary

Name of the device for connection

<a name='F-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-GSRReadings'></a>
### GSRReadings `constants`

##### Summary

List for galvanic skin response readings.

<a name='F-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-GSRTimestamps'></a>
### GSRTimestamps `constants`

##### Summary

List for galvanic skin response timestamps.

<a name='F-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-HRReadings'></a>
### HRReadings `constants`

##### Summary

List for heart rate readings.

<a name='F-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-HRTimestamps'></a>
### HRTimestamps `constants`

##### Summary

List for heart rate timestamps.

<a name='F-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-IBIReadings'></a>
### IBIReadings `constants`

##### Summary

List for inter-beat interval readings.

<a name='F-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-IBITimestamps'></a>
### IBITimestamps `constants`

##### Summary

List for inter-beat interval timestamps.

<a name='F-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-ReceiveDone'></a>
### ReceiveDone `constants`

##### Summary

ManualResetEvent signal completion of receive.

<a name='F-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-Sb'></a>
### Sb `constants`

##### Summary

Received data string.

<a name='F-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-SendDone'></a>
### SendDone `constants`

##### Summary

ManualResetEvent signal completion of send.

<a name='F-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-SocketConnection'></a>
### SocketConnection `constants`

##### Summary

Connection to be used in communications

<a name='F-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-SubscribedStreams'></a>
### SubscribedStreams `constants`

##### Summary

List of subscribed streams for synchronizing readings

<a name='F-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-TAGReadings'></a>
### TAGReadings `constants`

##### Summary

List for tag readings.

<a name='F-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-TAGTimestamps'></a>
### TAGTimestamps `constants`

##### Summary

List for tag timestamps.

<a name='F-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-TMPReadings'></a>
### TMPReadings `constants`

##### Summary

List for temperature readings.

<a name='F-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-TMPTimestamps'></a>
### TMPTimestamps `constants`

##### Summary

List for temperature timestamps.

<a name='M-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-ConnectCallback-System-IAsyncResult-'></a>
### ConnectCallback(ar) `method`

##### Summary

The callback for connect events to the TCP server.

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| ar | [System.IAsyncResult](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.IAsyncResult 'System.IAsyncResult') | The async result. |

<a name='M-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-Dispose'></a>
### Dispose() `method`

##### Summary

Closes the TCP client.

##### Parameters

This method has no parameters.

<a name='M-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-HandleDataStream-System-String[],System-String-'></a>
### HandleDataStream(DataPacket,response) `method`

##### Summary

Takes the TCP data and places it into the corresponding list.

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| DataPacket | [System.String[]](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.String[] 'System.String[]') | The data packet from TCP server. |
| response | [System.String](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.String 'System.String') | Unused. |

<a name='M-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-HandleDiscoverResponse-System-String[]-'></a>
### HandleDiscoverResponse(responses) `method`

##### Summary

Handles the device discovered response.

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| responses | [System.String[]](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.String[] 'System.String[]') | The responses from the TCP server. |

<a name='M-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-HandleDiscoverResponseBTLE-System-String[]-'></a>
### HandleDiscoverResponseBTLE(responses) `method`

##### Summary

Handles the discover device from BTLE connection.

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| responses | [System.String[]](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.String[] 'System.String[]') | The responses from the TCP server. |

<a name='M-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-HandleErrorCodes-System-String[]-'></a>
### HandleErrorCodes(ErrorResponse) `method`

##### Summary

Handles the received error codes.

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| ErrorResponse | [System.String[]](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.String[] 'System.String[]') | TCP error response to be parsed. |

<a name='M-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-HandleResponseFromEmpaticaBLEServer-System-String-'></a>
### HandleResponseFromEmpaticaBLEServer(response) `method`

##### Summary

Handles the response from the E4 TCP server.

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| response | [System.String](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.String 'System.String') | The response from the server. |

<a name='M-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-Receive'></a>
### Receive() `method`

##### Summary

Starts the receive process.

##### Parameters

This method has no parameters.

<a name='M-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-ReceiveCallback-System-IAsyncResult-'></a>
### ReceiveCallback(ar) `method`

##### Summary

TCP receive callback.

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| ar | [System.IAsyncResult](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.IAsyncResult 'System.IAsyncResult') | Async result. |

<a name='M-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-Send-System-Net-Sockets-Socket,System-String-'></a>
### Send(client,data) `method`

##### Summary

Send data to TCP server.

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| client | [System.Net.Sockets.Socket](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Net.Sockets.Socket 'System.Net.Sockets.Socket') | The TCP client. |
| data | [System.String](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.String 'System.String') | The data to send. |

<a name='M-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-SendCallback-System-IAsyncResult-'></a>
### SendCallback(ar) `method`

##### Summary

The callback for send completion.

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| ar | [System.IAsyncResult](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.IAsyncResult 'System.IAsyncResult') | The async result. |

<a name='M-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-StartClient'></a>
### StartClient() `method`

##### Summary

Starts the TCP client.

##### Parameters

This method has no parameters.

<a name='M-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-TagStressEvent'></a>
### TagStressEvent() `method`

##### Summary

Adds tag event to tag list.

##### Parameters

This method has no parameters.

<a name='T-MMIVR-BiosensorFramework-InputPipeline-SignalProcessing'></a>
## SignalProcessing `type`

##### Namespace

MMIVR.BiosensorFramework.InputPipeline

##### Summary

Public methods to compute features of biosensor signals.

<a name='M-MMIVR-BiosensorFramework-InputPipeline-SignalProcessing-ProcessAccSignal-System-Double[],System-Double-'></a>
### ProcessAccSignal(AccSignal,SamplingRate) `method`

##### Summary

Method to compute features for 3-axis accelerometer signals.

##### Returns



##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| AccSignal | [System.Double[]](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double[] 'System.Double[]') | An array of accelerometer data, packaged in [X, Y, Z] format in a single array. |
| SamplingRate | [System.Double](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double 'System.Double') | The sampling rate of the sensor. Defaults to 32 Hz. |

<a name='M-MMIVR-BiosensorFramework-InputPipeline-SignalProcessing-ProcessEdaSignal-System-Double[],System-Double-'></a>
### ProcessEdaSignal(EdaSignal,SamplingRate) `method`

##### Summary

Method to compute features of EDA signal.

##### Returns



##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| EdaSignal | [System.Double[]](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double[] 'System.Double[]') | An array containing the EDA signal readings. |
| SamplingRate | [System.Double](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double 'System.Double') | The sampling rate of the sensor. Defaults to 4 Hz. |

<a name='M-MMIVR-BiosensorFramework-InputPipeline-SignalProcessing-ProcessPpgSignal-System-Double[],System-Double,System-Double,System-Int32-'></a>
### ProcessPpgSignal(PpgSignal,SamplingRate,Threshold,PeakSize) `method`

##### Summary

Method to compute the features for the PPG signal.

##### Returns



##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| PpgSignal | [System.Double[]](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double[] 'System.Double[]') | An array containing readings from the PPG sensor. |
| SamplingRate | [System.Double](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double 'System.Double') | The sampling rate of the sensor. Defaults to 64 Hz. |
| Threshold | [System.Double](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double 'System.Double') | The signal threshold for feature computation. Defaults to 3.5. |
| PeakSize | [System.Int32](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Int32 'System.Int32') | The peak size needed for feature computation. Defaults to 3. |

<a name='M-MMIVR-BiosensorFramework-InputPipeline-SignalProcessing-ProcessTmpSignal-System-Double[]-'></a>
### ProcessTmpSignal(TempSignal) `method`

##### Summary

Method to compute the features of the temperature signal.

##### Returns



##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| TempSignal | [System.Double[]](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double[] 'System.Double[]') | An array containing the temperature sensor readings. |

<a name='T-MMIVR-BiosensorFramework-DataProcessing-TarvainenDetrending'></a>
## TarvainenDetrending `type`

##### Namespace

MMIVR.BiosensorFramework.DataProcessing

<a name='M-MMIVR-BiosensorFramework-DataProcessing-TarvainenDetrending-GetResidual-System-Double[],System-Double[]-'></a>
### GetResidual(TrendedSignal,DetrendedSignal) `method`

##### Summary

Removes the stationary signal from the EDA signal and returns the residual signal

##### Returns



##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| TrendedSignal | [System.Double[]](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double[] 'System.Double[]') | Original EDA signal |
| DetrendedSignal | [System.Double[]](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double[] 'System.Double[]') | Output of RemoveTrend function |

<a name='M-MMIVR-BiosensorFramework-DataProcessing-TarvainenDetrending-RemoveTrend-System-Double[],System-Int32,System-Double[]-'></a>
### RemoveTrend(InputSignal,Lambda,Filter) `method`

##### Summary

A time-varying finite-impulse-response high-pass filter for detrending
If using this in a published study, cite:
    Tarvainen, Mika P., Perttu O. Ranta-Aho, and Pasi A.Karjalainen. "An advanced detrending method with
        application to HRV analysis." IEEE Transactions on Biomedical Engineering 49.2 (2002): 172-175.

##### Returns

Filtered signal without a trend, double[]

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| InputSignal | [System.Double[]](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double[] 'System.Double[]') | Original EDA signal |
| Lambda | [System.Int32](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Int32 'System.Int32') | The 'smoothing' factor of the filter |
| Filter | [System.Double[]](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Double[] 'System.Double[]') | Coefficients for the 0, 1, 2 diagonals in the second derivative matrix.  Can be null, will default to standard value |

<a name='T-MMIVR-BiosensorFramework-MachineLearningUtilities-Train'></a>
## Train `type`

##### Namespace

MMIVR.BiosensorFramework.MachineLearningUtilities

##### Summary

Class to perform training of Microsoft.ML models.

<a name='M-MMIVR-BiosensorFramework-MachineLearningUtilities-Train-BuildAndTrainBinClassModels-Microsoft-ML-MLContext,Microsoft-ML-IDataView-'></a>
### BuildAndTrainBinClassModels(mlContext,TrainingSet) `method`

##### Summary

Trains binary classification models built-in to Microsoft.ML on the provided TrainingSet data.

##### Returns

List of models that can be used in performance benchmarks.

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| mlContext | [Microsoft.ML.MLContext](#T-Microsoft-ML-MLContext 'Microsoft.ML.MLContext') | The Microsoft.ML context to perform operations in. |
| TrainingSet | [Microsoft.ML.IDataView](#T-Microsoft-ML-IDataView 'Microsoft.ML.IDataView') | The time-series dataset to train the models on. |

<a name='M-MMIVR-BiosensorFramework-MachineLearningUtilities-Train-BuildAndTrainMultiClassModels-Microsoft-ML-MLContext,Microsoft-ML-IDataView-'></a>
### BuildAndTrainMultiClassModels(mlContext,TrainingSet) `method`

##### Summary

Trains multi-class models built-in to Microsoft.ML on the TrainingSet provided.

##### Returns

List of models that can be used in performance benchmarks.

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| mlContext | [Microsoft.ML.MLContext](#T-Microsoft-ML-MLContext 'Microsoft.ML.MLContext') | The Microsoft.ML context to perform operations in. |
| TrainingSet | [Microsoft.ML.IDataView](#T-Microsoft-ML-IDataView 'Microsoft.ML.IDataView') | The time-series dataset to train the models on. |

<a name='M-MMIVR-BiosensorFramework-MachineLearningUtilities-Train-BuildAndTrainRegressionModels-Microsoft-ML-MLContext,Microsoft-ML-IDataView-'></a>
### BuildAndTrainRegressionModels(mlContext,TrainingSet) `method`

##### Summary

Train regression classification models built-in to Microsoft.ML on the provided TrainingSet data.

##### Returns

List of models that can be used in performance benchmarks.

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| mlContext | [Microsoft.ML.MLContext](#T-Microsoft-ML-MLContext 'Microsoft.ML.MLContext') | The Microsoft.ML context to perform operations in. |
| TrainingSet | [Microsoft.ML.IDataView](#T-Microsoft-ML-IDataView 'Microsoft.ML.IDataView') | The time-series dataset to train the models on. |

<a name='M-MMIVR-BiosensorFramework-MachineLearningUtilities-Train-MultiToBin-System-Collections-Generic-List{MMIVR-BiosensorFramework-MachineLearningUtilities-ExtractedMultiFeatures}-'></a>
### MultiToBin(FeatureSet) `method`

##### Summary

Converts multi-class feature set to binary class representation.

##### Returns

Binary class representation of the input data.

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| FeatureSet | [System.Collections.Generic.List{MMIVR.BiosensorFramework.MachineLearningUtilities.ExtractedMultiFeatures}](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Collections.Generic.List 'System.Collections.Generic.List{MMIVR.BiosensorFramework.MachineLearningUtilities.ExtractedMultiFeatures}') | The feature set to convert. |

<a name='M-MMIVR-BiosensorFramework-MachineLearningUtilities-Train-MultiToReg-System-Collections-Generic-List{MMIVR-BiosensorFramework-MachineLearningUtilities-ExtractedMultiFeatures}-'></a>
### MultiToReg(FeatureSet) `method`

##### Summary

Converts multi-class feature dataset to regression class feature dataset.

##### Returns

Regression class representation of the input data.

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| FeatureSet | [System.Collections.Generic.List{MMIVR.BiosensorFramework.MachineLearningUtilities.ExtractedMultiFeatures}](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Collections.Generic.List 'System.Collections.Generic.List{MMIVR.BiosensorFramework.MachineLearningUtilities.ExtractedMultiFeatures}') | the feature set to convert. |

<a name='M-MMIVR-BiosensorFramework-MachineLearningUtilities-Train-PrintBinMetrics-Microsoft-ML-Data-BinaryClassificationMetrics-'></a>
### PrintBinMetrics(metrics) `method`

##### Summary

Prints the performance metrics of the binary classification test to Console.

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| metrics | [Microsoft.ML.Data.BinaryClassificationMetrics](#T-Microsoft-ML-Data-BinaryClassificationMetrics 'Microsoft.ML.Data.BinaryClassificationMetrics') | The metrics from the test set. |

<a name='M-MMIVR-BiosensorFramework-MachineLearningUtilities-Train-PrintMultiMetrics-Microsoft-ML-Data-MulticlassClassificationMetrics-'></a>
### PrintMultiMetrics(metrics) `method`

##### Summary

Prints the performance of the multi-class classification test to Console.

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| metrics | [Microsoft.ML.Data.MulticlassClassificationMetrics](#T-Microsoft-ML-Data-MulticlassClassificationMetrics 'Microsoft.ML.Data.MulticlassClassificationMetrics') | The metrics from the test set. |

<a name='M-MMIVR-BiosensorFramework-MachineLearningUtilities-Train-PrintRegMetrics-Microsoft-ML-Data-RegressionMetrics-'></a>
### PrintRegMetrics(metrics) `method`

##### Summary

Prints the performance metrics of the regression classification test to Console.

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| metrics | [Microsoft.ML.Data.RegressionMetrics](#T-Microsoft-ML-Data-RegressionMetrics 'Microsoft.ML.Data.RegressionMetrics') | The metrics from the test set. |

<a name='M-MMIVR-BiosensorFramework-MachineLearningUtilities-Train-RunBenchmarks-System-String,Microsoft-ML-ITransformer@,Microsoft-ML-ITransformer@,Microsoft-ML-ITransformer@-'></a>
### RunBenchmarks(DirectoryPath,BestRegModel,BestMultiModel,BestBinModel) `method`

##### Summary

Runs regression, multi-class, and binary classification tasks on the WESAD dataset and compares performance.  Returns the best performing model of each category.

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| DirectoryPath | [System.String](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.String 'System.String') | Path to directory with WESAD dataset. |
| BestRegModel | [Microsoft.ML.ITransformer@](#T-Microsoft-ML-ITransformer@ 'Microsoft.ML.ITransformer@') | The best regression ITransformer model. |
| BestMultiModel | [Microsoft.ML.ITransformer@](#T-Microsoft-ML-ITransformer@ 'Microsoft.ML.ITransformer@') | The best multi-class ITransformer model. |
| BestBinModel | [Microsoft.ML.ITransformer@](#T-Microsoft-ML-ITransformer@ 'Microsoft.ML.ITransformer@') | The best binary ITransformer model. |

<a name='M-MMIVR-BiosensorFramework-MachineLearningUtilities-Train-TrimFeatureSet-System-Collections-Generic-List{MMIVR-BiosensorFramework-MachineLearningUtilities-ExtractedMultiFeatures},System-Collections-Generic-List{System-Int32}-'></a>
### TrimFeatureSet(FeatureSet,LabelsToRemove) `method`

##### Summary

Removes specified labels from the dataset.

##### Returns

Trimmed list.

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| FeatureSet | [System.Collections.Generic.List{MMIVR.BiosensorFramework.MachineLearningUtilities.ExtractedMultiFeatures}](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Collections.Generic.List 'System.Collections.Generic.List{MMIVR.BiosensorFramework.MachineLearningUtilities.ExtractedMultiFeatures}') | The data to remove samples from. |
| LabelsToRemove | [System.Collections.Generic.List{System.Int32}](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Collections.Generic.List 'System.Collections.Generic.List{System.Int32}') | List of labels to remove from dataset. |

<a name='T-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-Window-Type'></a>
## Type `type`

##### Namespace

MMIVR.BiosensorFramework.DataProcessing.ThirdParty.DSPLib.DSP.Window

##### Summary

ENUM Types for included Windows.

<a name='T-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-Utilities'></a>
## Utilities `type`

##### Namespace

MMIVR.BiosensorFramework.Biosensors.EmpaticaE4

##### Summary

Utilities for sending commands to TCP server.

<a name='F-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-Utilities-AvailableDevices'></a>
### AvailableDevices `constants`

##### Summary

List of available devices on TCP server.

<a name='M-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-Utilities-ClearReadings-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-'></a>
### ClearReadings(E4Object) `method`

##### Summary

Clears the data lists.

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| E4Object | [MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient](#T-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient 'MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient') | The TCP client handling the responses. |

<a name='M-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-Utilities-ConnectDevice-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-'></a>
### ConnectDevice(E4Object) `method`

##### Summary

Connect device from TCP server.

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| E4Object | [MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient](#T-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient 'MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient') | The TCP client handling resposnes. |

<a name='M-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-Utilities-ConnectDeviceBTLE-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient,System-Int32-'></a>
### ConnectDeviceBTLE(E4Object,Timeout) `method`

##### Summary

Connect device over BTLE.

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| E4Object | [MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient](#T-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient 'MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient') | The TCP client handling resposnes. |
| Timeout | [System.Int32](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Int32 'System.Int32') | The timeout for the client. Defaults to -1. |

<a name='M-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-Utilities-DisconnectDevice-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-'></a>
### DisconnectDevice(E4Object) `method`

##### Summary

Disconnect device from TCP server.

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| E4Object | [MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient](#T-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient 'MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient') | The TCP client handling resposnes. |

<a name='M-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-Utilities-DisconnectDeviceBTLE-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-'></a>
### DisconnectDeviceBTLE(E4Object) `method`

##### Summary

Disconnect device connected over BTLE.

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| E4Object | [MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient](#T-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient 'MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient') | The TCP client handling resposnes. |

<a name='M-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-Utilities-GrabWindow-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient,System-String-'></a>
### GrabWindow(E4Object,Filepath) `method`

##### Summary

Grab data window.

##### Returns

List of List of data from TCP server.

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| E4Object | [MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient](#T-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient 'MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient') | The TCP client handling the responses. |
| Filepath | [System.String](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.String 'System.String') | The filepath to save out the data. Defaults to null. |

<a name='M-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-Utilities-ListDiscoveredDevices-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-'></a>
### ListDiscoveredDevices(E4Object) `method`

##### Summary

Lists the discovered device on the TCP server to the AvailableDevices list.

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| E4Object | [MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient](#T-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient 'MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient') | The TCP client handling resposnes. |

<a name='M-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-Utilities-ListDiscoveredDevicesBTLE-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-'></a>
### ListDiscoveredDevicesBTLE(E4Object) `method`

##### Summary

List devices discovered over BTLE.

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| E4Object | [MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient](#T-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient 'MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient') | The TCP client handling resposnes. |

<a name='M-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-Utilities-SaveReadings-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient,System-String-'></a>
### SaveReadings(E4Object,Filepath) `method`

##### Summary

Saves the readings from the lists of data to file.

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| E4Object | [MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient](#T-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient 'MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient') | The TCP client handling the responses. |
| Filepath | [System.String](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.String 'System.String') | The filepath to the file to save data to. |

<a name='M-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-Utilities-StartE4Server-System-String,System-String,System-String,System-String-'></a>
### StartE4Server(ServerPath,APIKey,IPaddress,Port) `method`

##### Summary

Starts the E4 server.

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| ServerPath | [System.String](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.String 'System.String') | The path to the exe. |
| APIKey | [System.String](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.String 'System.String') | The API key to complete connection. |
| IPaddress | [System.String](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.String 'System.String') | The IP address to the TCP server. Defaults to 127.0.0.1 |
| Port | [System.String](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.String 'System.String') | The port to connect through. Defaults to 8000. |

<a name='M-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-Utilities-StartE4ServerGUI-System-String-'></a>
### StartE4ServerGUI(ServerPath) `method`

##### Summary

Starts the E4 server GUI with a process call.

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| ServerPath | [System.String](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.String 'System.String') | The path to the server exe. |

<a name='M-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-Utilities-StartStreaming-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-'></a>
### StartStreaming(E4Object) `method`

##### Summary

Starts the data streaming on the E4 device.

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| E4Object | [MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient](#T-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient 'MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient') | The TCP client handling resposnes. |

<a name='M-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-Utilities-StartTimer-System-Single-'></a>
### StartTimer(Seconds) `method`

##### Summary

Starts the windowed data reading timer.

##### Returns

Returns the timer instance.

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| Seconds | [System.Single](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.Single 'System.Single') | The number of seconds for the timer. |

<a name='M-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-Utilities-SubscribeToStream-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient,MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-DeviceStreams-'></a>
### SubscribeToStream(E4Object,Stream) `method`

##### Summary

Subscribes device to data stream.

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| E4Object | [MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient](#T-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient 'MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient') | The TCP client handling resposnes. |
| Stream | [MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient.DeviceStreams](#T-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-DeviceStreams 'MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient.DeviceStreams') | The stream to unsubscribe from. |

<a name='M-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-Utilities-SuspendStreaming-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-'></a>
### SuspendStreaming(E4Object) `method`

##### Summary

Suspends the data streaming on the E4 device.

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| E4Object | [MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient](#T-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient 'MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient') | The TCP client handling resposnes. |

<a name='M-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-Utilities-UnsubscribeToStream-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient,MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-DeviceStreams-'></a>
### UnsubscribeToStream(E4Object,Stream) `method`

##### Summary

Unsubscribes device from data streams.

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| E4Object | [MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient](#T-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient 'MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient') | The TCP client handling resposnes. |
| Stream | [MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient.DeviceStreams](#T-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-ServerClient-DeviceStreams 'MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.ServerClient.DeviceStreams') | The stream to unsubscribe from. |

<a name='T-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-Window'></a>
## Window `type`

##### Namespace

MMIVR.BiosensorFramework.DataProcessing.ThirdParty.DSPLib.DSP

<a name='M-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-Window-Coefficients-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-Window-Type,System-UInt32-'></a>
### Coefficients(windowName,points) `method`

##### Summary

Calculates a set of Windows coefficients for a given number of points and a window type to use.

##### Returns

double[] array of the calculated window coefficients

##### Parameters

| Name | Type | Description |
| ---- | ---- | ----------- |
| windowName | [MMIVR.BiosensorFramework.DataProcessing.ThirdParty.DSPLib.DSP.Window.Type](#T-MMIVR-BiosensorFramework-DataProcessing-ThirdParty-DSPLib-DSP-Window-Type 'MMIVR.BiosensorFramework.DataProcessing.ThirdParty.DSPLib.DSP.Window.Type') |  |
| points | [System.UInt32](http://msdn.microsoft.com/query/dev14.query?appId=Dev14IDEF1&l=EN-US&k=k:System.UInt32 'System.UInt32') |  |

<a name='T-MMIVR-BiosensorFramework-Biosensors-EmpaticaE4-Utilities-WindowedDataReady'></a>
## WindowedDataReady `type`

##### Namespace

MMIVR.BiosensorFramework.Biosensors.EmpaticaE4.Utilities

##### Summary

Delegate to alert that windowed data is ready. Triggered by stopwatch.
---
title: 'Biosensor Framework: A C# Library for Affective Computing'
tags:
  - C#
  - Affective Computing
  - Biometrics
  - Machine Learning
  - Unity
  - Microsoft.ML
authors:
  - name: Walker Arce
    orcid: 0000-0003-0819-0710
    affiliation: 1
  - name: James Gehringer, PhD
    orcid: 0000-0003-3457-2288
    affiliation: 1
affiliations:
 - name: Virtual Reality Laboratory in the Munroe Meyer Institute at the University of Nebraska Medical Center
   index: 1
date: 15 May 2021
bibliography: paper.bib
---

# Summary

Developed in the Virtual Reality Laboratory at the Munroe Meyer Institute, the Biosensor Framework library provides an interface for interacting with biosensors and performing affective computing tasks on their output data streams.  Currently, it natively supports the Empatica E4 biosensor and communicates through a TCP connection with their provided server.  

This library provides the following capabilities:
- Connection to the [TCP client](https://github.com/empatica/ble-client-windows) that Empatica provides for connecting their devices to a user's PC
- Parsing and packing of commands to communicate bidirectionally with the TCP client
- Collection and delivery of biometric readings to an end user's software application
- Extraction of literature supported feature vectors to allow characterization of the biological signal features
- Microsoft.ML support for inferencing on the feature vectors natively in C# in both multi-class and binary classification tasks
- Support for simulating a data stream for debugging and testing applications without an Empatica E4 device
- Support for parsing open source Empatica E4 datasets and training machine learning models on them (i.e., WESAD [@schmidt2018introducing])
- Support for interfacing new sensors into the same pipeline
- Contains models trained on WESAD dataset in the Microsoft.ML framework.

# Background

Biosensor Framework is a C# library that handles the full process of gathering biometric data from a body-worn sensor, transforming it into handcrafted feature vectors, and delivering an inferencing result in thirty-five lines of code.  Purpose built to handle the Empatica E4 biosensor, the Empatica provided TCP client is wrapped in a static C# class to provide convenient function calls and error handling to simplify code structure.  With the introduction of Microsoft.ML, the entire pipeline is handled in C# managed code and is fully compatible with the Unity game engine, allowing for affective computing to be easily introduced into existing virtual reality therapy tools.  The decoupled structure of the library also allows for new biosensors to be introduced into the pipeline, requiring only a single code file to be added to retrieve data in a list of doubles.  

Generally, affective computing is concerned with the physiological changes of the subject and how that relates to their internal state.  Most implementations of affective computing look at electrodermal activity (EDA), like a lie detector test, to measure the onset of stress before the subject may even be aware of it [@gjoreski2016continuous; @healey2000wearable; @healey2005detecting; @kutt2018towards].   There is a growing body of literature that incorporates the use of multiple body sensors, primarily worn on the wrist, such as photoplethysmography (PPG), three axis acceleration (ACC), and skin temperature (TMP) [@choi2011development; @gjoreski2016continuous; @kaczor2020objective; @kleiman2019using; @kutt2018towards; @van2019heartpy; @zhou2017indexing].  As these sensors are miniaturized, they can be easily incorporated into wrist-worn packages and can be fused with other sensor modalities, such as respiratory sensors, eye trackers, or body trackers, to estimate a subject’s emotional state more accurately.

Biosensor Framework incorporates the algorithms used in the WESAD dataset for computing feature vectors for all sensors on the Empatica E4: EDA, PPG, ACC, and TMP [@schmidt2018introducing].  Additionally, the E4 provides estimations on the heart rate and inter-beat interval, as well as a button for manually tagging the data stream, which has been implemented in a function to mark relevant events.  To detrend the EDA signal, a time-varying, finite-impulse-response, high-pass filter is used [@tarvainen2002advanced].  For detecting major responses in the EDA, the method proposed by Healey is used [@healey2000wearable].  To perform the Fast Fourier Transform, the [DSPLib](https://www.codeproject.com/Articles/1107480/DSPLib-FFT-DFT-Fourier-Transform-Library-for-NET-6) project was incorporated to allow real-time computation of the signal spectra.  Finally, a stillness metric was added to compute a unitless measure of movement from the accelerometer signals [@chang2012wireless].

# Statement of need

In the application of virtual reality therapy tools for patients with autism spectrum disorder, the continual management of stress is of utmost importance [@kildahl2019identification; @bishop2015relationship].  These tools can introduce stressful situations such as navigating an airport, handling an airplane ride, or crossing a busy intersection, which can be difficult or impossible to simulate in a controlled clinical environment [@poyade2017using; @miller2020virtual; @saiano2015natural].  This stress can be monitored using biosensors and literature supported machine learning models, which gives the clinician an additional layer of information for tuning the environment for its intended therapeutic effect.  Additionally, the library can be used in applications for monitoring and managing ambulatory stress for at-risk patients, allowing another layer of information for medical personnel to assess the state of their patients [@kleiman2019using].

The most common tool to implement virtual reality therapy tools is the Unity game engine, which runs a C# interpreter based on the Mono compiler.  To reasonably use common open-source machine learning tools, such as Python, MATLAB, or R, would require the use of additional TCP servers to act as the intermediary between the virtual environment and the biosensor readings.  This adds an additional layer to the technology stack and, for Unity developers, requires the compilation and management of multiple projects. The Biosensor Framework library removes this additional layer and simplifies the machine learning down to simple function calls.

# Affect Classification Model Performance

We split the WESAD dataset into the binary classification and multi-class organization as specified in [@schmidt2018introducing] and calculated the feature vectors accordingly.  We used a train-test split of 0.3 and trained using the standard Fit method built into the Microsoft.ML models.  The binary classification task was split into _stress vs. non-stress_ and the multi-class classification task was split into _baseline vs. stress vs. amusement_.  The regression tasks also used the latter dataset.  The metrics were calculated from the test fit and were automatically generated by Microsoft.ML.

Regression Model Metrics

|                            |     Fast Forest    |     Fast Tree    |     Fast Tree Tweedie    |     LBFGS    |
|----------------------------|--------------------|------------------|--------------------------|--------------|
|     Loss Function          |     19.54          |     11.34        |     12.27                |     37.11    |
|     Mean Absolute Error    |     28.33          |     16.64        |     16.63                |     45.89    |
|     Mean Squared Error     |     19.54          |     11.14        |     12.27                |     37.11    |
|     RMS Error              |     44.20          |     33.67        |     35.03                |     60.92    |
|     R Squared              |     65.17          |     **79.79**        |     78.12                |     26.11    |

Multi-Class Classification Metrics

|                          |     Fast Forest    |     Fast Tree    |     LBFGS Regression    |     LBFGS Max Entropy    |
|--------------------------|--------------------|------------------|-------------------------|--------------------------|
|     Log Loss             |     31.30          |     22.27        |     45.11               |     41.51                |
|     LL Reduction         |     68.51          |     77.59        |     51.82               |     55.67                |
|     Macro Acc            |     84.59          |     **93.78**        |     68.80               |     71.97                |
|     Micro Acc            |     89.06          |     **95.40**        |     84.76               |     85.40                |
|     Class 0 Precision    |     93.72          |     95.90        |     83.09               |     84.69                |
|     Class 0 Recall       |     93.36          |     97.82        |     95.85               |     94.54                |
|     Class 0 Loss         |     22.03          |     9.92         |     23.13               |     23.59                |
|     Class 1 Precision    |     83.98          |     98.34        |     91.46               |     90.56                |
|     Class 1 Recall       |     92.97          |     94.73        |     91.09               |     91.30                |
|     Class 1 Loss         |     25.96          |     27.68        |     42.37               |     34.16                |
|     Class 2 Precision    |     83.21          |     88.53        |     63.64               |     85.40                |
|     Class 2 Recall       |     67.26          |     88.79        |     19.44               |     30.09                |
|     Class 2 Loss         |     71.23          |     52.41        |     152.03              |     140.35               |

Binary Classification Metrics

|                           |     Fast Forest    |     Fast Tree    |
|---------------------------|--------------------|------------------|
|     Accuracy              |     92.35          |     **96.85**        |
|     AUPRC                 |     97.95          |     99.05        |
|     AURC                  |     97.39          |     98.88        |
|     F1 Score              |     92.86          |     97.06        |
|     Negative Precision    |     92.34          |     97.07        |
|     Negative Recall       |     91.19          |     96.17        |
|     Positive Precision    |     92.36          |     96.67        |
|     Positive Recall       |     93.37          |     97.46        |

For our evaluation of the machine learning framework, we focused on the binary classification task of _stress v. non-stress_.  Microsoft.ML has built in models that allow for quick training, evaluation, and deployment of machine learning applications.  The ones that we evaluated were the Fast Forest and Fast Tree models.  The models evaluated performed very well in the binary classification task and the Fast Tree outperformed the Fast Forest.  For all three classification tasks the Fast Tree outperformed all other models.  

# Acknowledgements

We would like to thank the Integrated Center for Autism Spectrum Disorder and the Munroe Meyer Guild at the Munroe Meyer Institute at the University of Nebraska Medical Center for supporting our work.

# References
