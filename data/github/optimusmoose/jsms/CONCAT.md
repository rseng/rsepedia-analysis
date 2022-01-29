# JS-MS
*JavaScript Mass Spectrometry*, a visualization GUI for Mass Spectrometry signal data. Presents a robust, 3-dimensional graph view that allows the creation, editing, and saving of XIC and envelope annotations in addition to panning, rotating, and zooming actions. Interfaces via HTTP to the MsDataServer API for access to MS data stored in the MzTree format.

## Citations
If you are publishing any work related to the execution of this software, *please cite the following papers* which describe it:
> @article{rosen2017js,
  title={JS-MS: a cross-platform, modular javascript viewer for mass spectrometry signals},
  author={Rosen, Jebediah and Handy, Kyle and Gillan, Andr{\'e} and Smith, Rob},
  journal={BMC bioinformatics},
  volume={18},
  number={1},
  pages={469},
  year={2017},
  publisher={BioMed Central}
}

> @article{handy2017fast,
  title={Fast, axis-agnostic, dynamically summarized storage and retrieval for mass spectrometry data},
  author={Handy, Kyle and Rosen, Jebediah and Gillan, Andr{\'e} and Smith, Rob},
  journal={PloS one},
  volume={12},
  number={11},
  pages={e0188059},
  year={2017},
  publisher={Public Library of Science}
}

> @article{gutierrez2019xnet,
  title={XNet: A Bayesian approach to Extracted Ion Chromatogram Clustering for Precursor Mass Spectrometry Data},
  author={Gutierrez, Matthew and Handy, Kyle and Smith, Rob},
  journal={Journal of proteome research},
  year={2019},
  publisher={ACS Publications}
}

## Dependencies
* Run: Java 8 JRE or higher, WebGL (included in major browsers)
* Build (optional): Java 8 JDK, Apache Maven

## Usage: GUI
- Navigate to /msDataServer/target in a file manager
- Double click msDataServer-*version*.jar
- You may need to make the .jar file executable first. This procedure varies depending on your system configuration.
- A server window will open. You will now be able to open a mass spectrometry data file and start the HTTP server and launch JS-MS.
- The msDataServer jar file is self-contained can be run independently from the project folder as a standalone package
- Enter desired port number in the MsDataServer window
- Click **Start Server**, JS-MS will open in your default browser and MsDataServer will be activated
- Click **Open File**, a file chooser dialog will appear over the MsDataServer window
- Select an *mzML* or *MzTree* file
- Begin interacting with the graph view
- **NOTE**: JSMS' file storage and retrieval system incurs overhead costs that amortize well on files of reasonable or large size, but which causes notable delays even if loading small files. Please be patient! This is a one-time cost when building a tree with a new raw file. Subsequent file loads are very fast once an mzTree file is built.

## Usage: Terminal/Command-line start
- `java -jar msDataServer/target/msDataServer-<version>.jar`
- Consider expanding the maximum memory available to the application.
  This can be done by passing the `-Xmx` flag. For example, add `-Xmx8g` to give
  the application a maximum 8 gigabytes of memory.

#### See *instruction.html* for instructions on interacting with JS-MS
*instruction.html* is accessible by clicking the **?** button in the top-right of the graph view.

#### Creating ground truth data
- *Isotopic Trace Mode.* When a user enters isotopic trace mode, they are given the option to create a new trace or select an existing trace to edit. Each time a new trace is created, the trace is given an ID and color. Users select the points belonging to the trace by clicking and dragging a rectangle over the desired points to highlight them in the given color. The same procedure is used to edit an existing trace, only the control key is depressed while drawing the rectangle.

- *Isotopic Envelope Mode.* After the user has identified isotopic traces, they can group them together with isotopic envelope mode. Similar to isotopic trace mode, this mode creates a new envelope ID and color for each new envelope created. The user then selects all isotopic traces that belong to the same group. Isotopic traces can be grouped by clicking each trace or simply by dragging a line across all traces in an envelope. To help the user distinguish which isotopic traces belong together, the ruler tool shows m/z intervals corresponding to specific charge states. The ruler will appear wherever the mouse is placed when users select a number from the keypad. The ruler moves with the graph as the user zooms or pans, and will remain present until the user hits the tilde key. The m/z distance displayed is 1/z, where z is the number selected and the charge state of a hypothetical compound at the given mass. Users can also toggle between 2-D and 3-D mode while in either isotopic trace or isotopic envelope mode to ensure peak alignment. Isotopic traces can be added to existing isotopic envelopes at any time following this procedure and they can be removed in the same way while depressing the control key.

- *Mark as Noise Button.* When all distinguishable points in the current view have been annotated the user can mark all other points in the view as noise. When a point is marked as noise it will be colored gray in the view and given an ID of -1 when exported to .csv. To prevent users from marking unseen points as noise, the graph view must be displaying a number of points below the point threshold to ensure that the user is viewing every point within the (m/z, RT) coordinates and none are hidden.

#### Exporting ground truth data
- In the server window (where you start the server and launch the browser), you can export segmented data to a .csv that will feature all points in the file in the form: m/z, RT, intensity, isotopic trace ID, isotopic envelope ID. Isotopic IDs are unique, meaning they are not dependent on the isotopic envelope ID. A trace/envelope ID of zero means that point has not yet been assigned to a trace/envelope, and an ID of -1 means it has been assigned to noise using the "mark as noise" button.

#### Creating and editing bookmarks
- Bookmarks can be created in the program by first displaying the bookmarks bar by clicking on the bookmarks button in the lefthand bar. Next, you enter a label, m/z value, and RT value in the respective fields on the list. Existing entries can be edited by clicking the edit button beside the entry, or deleted by clicking the delete button beside the entry.

- To create a bookmark list outside of JS-MS, create a .tsv file that contains, for each line/bookmark, a text label, a float value for m/z, and a float value for RT. In JS-MS, you can open this bookmark list using the open button in the bookmark list. You can also export a list using the export button in the bookmarks list view.

- NOTE: Re-opening a bookmark list in JS-MS will overwrite any changes you have made to the list. Exporting a bookmark list will export the current version of the list as displayed. This is important to know if you edit or remove entries to the bookmark list while in JS-MS.

#### Inspecting data using bookmarks
- In the righthand panel, locate the "Jump To" radio button and ensure "Next Bookmark" is selected. Create or load a bookmark list as described above. Click on the bookmarks button in the lefthand bar to display the bookmark list. You can now navigate to a specific bookmark by clicking on the name of the bookmark, or the next bookmark in visited order by clicking the jump button (right arrow) on the lefthand bar.

#### Inspecting data using "Jump to Window"
- In the righthand panel, locate the "Jump to Window" area. Enter the m/z and RT window you would like to graph, and click the check button. The graph will respond by plotting the described area.

## Build (optional)
- Run `mvn package` from within the project root directory

## Modules
- msDataServer: mzTree core and HTTP API. Includes a GUI with buttons to perform most actions provided by other modules
- tracesegmentation: creates segmentation data, extensible with different file format accessors
- xnet: trace clustering (traces -> envelopes) and data types used by mzTree storage
- correspondence: matching of finished envelope data from several files to highlight similarities and differences between data sets

## License
This work is published under the MIT license.
# Mass Spec - Trace Segmentation

The trace segmentation algorithm is derived from Rob Smith's
original Ruby code and lives in the `edu.umt.ms.traceSeg` namespace.
The "meat" of the code is in the `TraceSegmenter` Java class.

Trace segmentation depends on a `PointDatabaseConnection`; a reference
`HttpPointDatabaseConnection` has been provided with support for the
`MsHttpServer` API. Additionally, a `TraceParametersProvider` is used
to obtain the tuning parameters. This plugin may be removed in the
future if good enough defaults are determined. If desired, the two
interfaces can be reimplemented by consumers of this library to work
efficiently with their own data structures and other needs.

Trace segmentation works against a running MS HTTP server with the
reference implementation, and msViz has support for running
direct-access trace segmentation (as opposed to slow HTTP) in several
harnesses and as an option in the GUI (_Cluster_ tab, _Trace_ button).

## Metrics
In addition to trace segmentation itself, metrics have been implemented
in `edu.umt.ms.traceSeg.metrics`; Purity, NMI, SSE, and NTCD are
implemented at this time. These all run against CSV exports of
trace-segmented mass spec data.

## Usage
The project is packaged as a maven artifact. For code usage, see
`src/main/java/edu/umt/ms/traceSeg/TraceSegmentation.java` for an
example. This example uses an `HttpPointDatabaseConnection`
(connecting to a running MS HTTP server) and uses
`ConsoleTraceParametersProvider` to interactively request parameters
from the user at runtime.

## Algorithm
TraceSegmenter.java implements a probability-based single-pass assignment
of points to traces. The point with the highest intensity is assigned
first, followed by the next-highest point, until a point with a lower
intensity than the limit is reached. For each point, all "nearby" existing
traces, that is, traces close enough that the point could be included
in them, are checked. The trace with the highest probability score is
selected, as long as it is within the probability threshold. If no trace
scores the point within the probability threshold, a new trace is created.
The point is assigned to the trace and the trace's statistics are updated
for use when inspecting future points.

The probability score is a function of distance in m/z, distance in RT,
and intensity. For m/z and RT, closer points are higher-scoring. For
intensity, less-intense points are more likely to be associated with any
trace and more-intense points are less likely to be associated with any
trace. Additionally, no point can be added to an existing trace unless
its intensity is strictly lower than the intensities of all other points
currently in the trace.
# msDataServer Web API

###HTTP GET /*
(.html, .js, etc.)
		
Serves the contents of the "jsms" folder (i.e. the user interface files).

---

###HTTP GET /api/v2/filestatus
Checks the status of the server's open-file process. Returns additional file data if a file has been loaded and data model is ready.. 

####Server response:
		
	HTTP 200 (OK): File has been loaded and ata model is ready, begin querying for points.
		Payload: { "mzmin" : float, "mzmax" : float, "rtmin" : float, "rtmax" : float, "intmin" : float, "intmax" : float, "pointcount": integer, "progress": float }
	HTTP 204 (No Content): No file has been selected, open a file before continuing.
	HTTP 400 (Bad Request): Malformed request (usually a missing parameter)
	HTTP 403 (Forbidden): The file is currently being selected
	HTTP 406 (Not Acceptable): The selected file is of the wrong file format, reselect file before continuing.
	HTTP 409 (Conflict): The server is selecting a file or processing the selected file. Continue checking file status.

---

###HTTP GET /api/v2/getpoints

Queries the database for points to display on the graph. Requests a specific number of points within the given bounds. Server determines the detail level that would provide the same or more points than requested and samples this set to the requested number of points.

####URL parameters:
	mzmin (double): mz lower bound (0 for global mz minimum)
	mzmax (double): mz upper bound (0 for global mz maximum)
	rtmin (float): rt lower bound (0 for global rt minimum)
	rtmax (float): rt upper bound (0 for global rt maximum)
	numpoints (int): the number of points to be returned (0 for no limit)

####Server response:
	HTTP 200 (OK): Query successfully serviced, returning points
		Payload: [[<pointId>,<traceId>,<mz>,<rt>,<intensity>], ... ]
	HTTP 204 (No Content): No file has been selected, open a file before continuing.
	HTTP 400 (Bad Request): Malformed request, missing parameter or invalid query range (i.e. mzmin > mzmax).
	HTTP 406 (Not Acceptable): The previously selected file is of the wrong file format, reselect file before continuing.
	HTTP 409 (Conflict): The server is selecting a file or processing the selected file. Continue checking file status.

---

###HTTP GET /api/v2/gethighestuntraced
		
Queries the database for the highest intensity point that has not yet been assigned to a trace.
	
####Server response:
	HTTP 200 (OK): Query serviced, returning a single point
		[<pointId>, <traceId>, <mz>, <rt>, <intensity>]
	HTTP 204 (No Content): No file has been selected, open a file before continuing.
	HTTP 400 (Bad Request): Malformed request, missing parameter or invalid query range (i.e. mzmin > mzmax).
	HTTP 406 (Not Acceptable): The previously selected file is of the wrong file format, reselect file before continuing.
	HTTP 409 (Conflict): The server is selecting a file or processing the selected file. Continue checking file status.

---

###HTTP GET /api/v2/gettracemap
Returns the current Trace Map containing traceID->envelopeID mappings.

####Server response:
	HTTP 200 (OK): Query successfully serviced, returning Trace Map
		Payload: { "<traceID>":<envID>, ... }
	HTTP 204 (No Content): No file has been selected, open a file before continuing.
	HTTP 406 (Not Acceptable): The previously selected file is of the wrong file format, reselect file before continuing.
 	HTTP 409 (Conflict): The server is selecting a file or processing the selected file. Continue checking file status.

---

###HTTP GET /api/v2/getnextids
Returns the next trace ID and envelope ID. "Next" is defined as the ID after the largest used ID for that type.

####Server response:
	HTTP 200 (OK): Query successfully serviced, returning next IDs
		Payload: { "nextTrace":<traceId>, "nextEnvelope":<envelopeId> }
	HTTP 204 (No Content): No file has been selected, open a file before continuing.
	HTTP 406 (Not Acceptable): The previouslyj selected file is of the wrong file format, reselect file before continuing.
 	HTTP 409 (Conflict): The server is selecting a file or processing the selected file. Continue checking file status.

---

### HTTP POST /api/v2/updatesegmentation
Sends segmentation modifications to the server to update the mzTree data model. 
####POST Parameters (JSON):
	[ <action>, <action>, ... ]
	Action format one of following:
		{ "type": "undo" }
		{ "type": "redo" }
		{ "type": "set-trace", "trace": <traceid>, "points": [<pointid>, ...] }
		{ "type": "set-envelope", "envelope": <envelopeid>, "traces": [<traceid>, ...] }
		{ "type": "rectangle", "bounds": [<lower_mz>, <upper_mz>, <lower_rt>, <upper_rt>], "id": <traceID>, "isAdd": <boolean> }

	
####Server response:
	HTTP 200 (OK): Successfully updated segmentation data.
	HTTP 500 (INTERNAL_SERVER_ERROR): Failed to update segmentation data.

---

### HTTP POST /api/v2/getenvelopeinfo
Requests that the server recompile trace data, and returns envelope statistics.
####POST Parameters (urlencoded):
	id1,id2,id3,...
	where each id is an envelope ID

####Server response:
	HTTP 200 (OK): Successful query, returning envelope information
		Payload: { "env_id": [ env_mz, env_rt, env_intensity, [ [trace1_mz, trace1_intensity], [trace2_mz, trace2_intensity], ... ] ], ... }
	HTTP 400 (Bad Request): Malformed request (a given ID was not parseable as an integer)
	HTTP 500 (INTERNAL_SERVER_ERROR): An error occurred.
