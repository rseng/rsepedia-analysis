# Via Appia Changelog

This changelog is updated manually when deploying major features to the application

### September 2021

*   technical: using composition API
*   feat: Navigate between stories
# Via Appia Viewer

> [Beginner Tutorial](./TUTORIAL.MD)

![Netlify Status](https://api.netlify.com/api/v1/badges/ff9d22c2-1548-448b-a6c8-f54573e6df3e/deploy-status)

Previous project (2015-2016): https://github.com/Via-Appia/PattyVis

## Build Setup

Enable `eslint --fix on file save` in your code editor.

```bash
# install dependencies
$ yarn install

# serve with hot reload at localhost:3000
$ yarn dev

# generate static project
$ yarn generate
```
ok


Archeological pointcloud viewer for use in a museum setting. Based on storylines (narratives) that are defined by an artist or expert, a user should be able to view multiple _storylines_ , all consisting of multiple _pages_.

## high resolution pointcloud locally:

You need to download and place the pointclouds data into the `/static/pointclouds/highres`.  
The structure looks like:

![img.png](img.png)

## Upload the PointCloud data in cloud storage

*   You need to have installed locally [gsutils](https://cloud.google.com/storage/docs/gsutil_install)
*   log in you google account to get writing permissions
*   Navigate to the root where the data folder is placed and start the copy of the files to the google cloud storage: 

```shell
gsutil -m cp -r ./data gs://via-appia-20540.appspot.com
```

*    Access the cloud storage dashboard [here](https://console.cloud.google.com/storage/browser/via-appia-20540.appspot.com)

# Enable settings locally

You can change local setting by creating a `.env` file and enabling the settings you want to have:

```shell
LOCAL_POINTCLUDS = true
POINTS_BUDGET = 1000000
```

## Converting big LAS files to LAZ and viewing it in Potree Desktop

To visualize point clouds, which are usually in .las format, in potree they need to be covertered using the potree converter tooling. Using the the [PotreeConverter 2](https://github.com/potree/PotreeConverter/releases/tag/2.0) is a good option, since it generated only 3 files. However, it does not compress these files which results in  a converted point cloud which will be as big as the original LAS file. It therefore makes sense to use the .laz conversion from the [PotreeConverter 1.7](https://github.com/potree/PotreeConverter/releases/tag/1.7). Converting a LAS file (in our case 21GB) to LAZ for potree is easiest using windows. It is possible to compile the PotreeConverter and LasZIP for Linux and Mac but you have to do it yourself. 

### What you need:

*   Your LAS file(s), (in our case has been a 21GB size files).
*   Lastools
    *   Download here https://rapidlasso.com/lastools/     
*   PotreeConverter
    *   Download here [PotreeConverter\_1.7\_windows\_x64.zip](https://github.com/potree/PotreeConverter/releases/tag/1.7)
*   Potree Desktop (optional, use it to test your converted files)
    *   Web: [https://github.com/potree/PotreeDesktop/releases/tag/1.8](https://github.com/potree/PotreeDesktop/releases/tag/1.8)
    *   File: [**PotreeDesktop\_1.8\_windows\_x64.zip**](https://github.com/potree/PotreeDesktop/releases/download/1.8/PotreeDesktop_1.8_windows_x64.zip)

First we need to merge the various .las files (if any) into one file. To do so we use lastools. Browse to the bin folder and run: 

`las2las -i C:\\...\\001.las C:\\...\\002.las C:\\...\\003.las ... -merged -o C:\\...\\merged.las`   

Next we need to convert this .las to using the potree converter 1.7 and specifically indicate that it need to be compressed as .laz .

Navigate to the folder where PotreeConverter is, and run the command (replace the \<names>):

`./PotreeConverter.exe .\<fileName>.las -o ./<outputDirectory> --output-format LAZ`

Once the process has finished, you can drag and drop the new output directory to PotreeDesktop 1.8 to test it. In our repository the files need to be placed in `\static\pointclouds\highres\` .
# STATIC

**This directory is not required, you can delete it if you don't want to use it.**

This directory contains your static files.
Each file inside this directory is mapped to `/`.
Thus you'd want to delete this README.md before deploying to production.

Example: `/static/robots.txt` is mapped as `/robots.txt`.

More information about the usage of this directory in [the documentation](https://nuxtjs.org/guide/assets#static).


# About

* Potree is a free open-source WebGL based point cloud renderer for large point clouds. It is based on the [TU Wien Scanopy project](https://www.cg.tuwien.ac.at/research/projects/Scanopy/) and research projects [Harvest4D](https://harvest4d.org/), [GCD Doctoral College](https://gcd.tuwien.ac.at/) and [Superhumans](https://www.cg.tuwien.ac.at/research/projects/Superhumans/).
* Newest information and work in progress is usually available on [twitter](https://twitter.com/m_schuetz)
* Contact: Markus Schütz (mschuetz@potree.org)
* References: 
    * [Potree: Rendering Large Point Clouds in Web Browsers](https://www.cg.tuwien.ac.at/research/publications/2016/SCHUETZ-2016-POT/SCHUETZ-2016-POT-thesis.pdf) (2016)
    * [Fast Out-of-Core Octree Generation for Massive Point Clouds](https://www.cg.tuwien.ac.at/research/publications/2020/SCHUETZ-2020-MPC/) (2020)
    
<a href="http://potree.org/wp/demo/" target="_blank"> ![](./docs/images/potree_screens.png) </a>

# Getting Started

### Install on your PC

Install [node.js](http://nodejs.org/)

Install dependencies, as specified in package.json, and create a build in ./build/potree.

```bash
npm install
```

### Run on your PC

Use the `npm start` command to 

* create ./build/potree 
* watch for changes to the source code and automatically create a new build on change
* start a web server at localhost:1234. 

Go to http://localhost:1234/examples/ to test the examples.

### Deploy to a server

* Simply upload the Potree folderm with all your point clouds, the build directory, and your html files to a web server.
* It is not required to install node.js on your webserver. All you need is to host your files online. 

### Convert Point Clouds to Potree Format

Download [PotreeConverter](https://github.com/potree/PotreeConverter) and run it like this:

    ./PotreeConverter.exe C:/pointclouds/data.las -o C:/pointclouds/data_converted

Copy the converted directory into &lt;potreeDirectory&gt;/pointclouds/data_converted. Then, duplicate and rename one of the examples and modify the path in the html file to your own point cloud.

# Downloads

* [Potree](https://github.com/potree/potree/releases)
* [PotreeConverter ](https://github.com/potree/PotreeConverter/releases) - Convert your point cloud to the Potree format.
* [PotreeDesktop ](https://github.com/potree/PotreeDesktop/releases) - Desktop version of Potree. Allows drag&drop of point clouds into the viewer.

# Examples

<table>
	<tr>
		<td style="padding: 0px">
			<a href="http://potree.org/potree/examples/viewer.html" target="_blank">
				<img src="examples/thumbnails/viewer.png" width="100%" />
			</a>
		</td><td>
			<a href="http://potree.org/potree/examples/ca13.html" target="_blank">
				<img src="examples/thumbnails/ca13.png" width="100%" />
			</a>
		</td><td>
			<a href="http://potree.org/potree/examples/cesium_retz.html" target="_blank">
				<img src="examples/thumbnails/cesium_retz.png" width="100%" />
			</a>
		</td><td>
			<a href="http://potree.org/potree/examples/classifications.html" target="_blank">
				<img src="examples/thumbnails/classifications.jpg" width="100%" />
			</a>
		</td><td>
			<a href="http://potree.org/potree/examples/features_sorvilier.html" target="_blank">
				<img src="examples/thumbnails/features_sorvilier.png" width="100%" />
			</a>
		</td><td>
			<a href="http://potree.org/potree/examples/toolbar.html" target="_blank">
				<img src="examples/thumbnails/toolbar.jpg" width="100%" />
			</a>
		</td>
	</tr>
	<tr>
		<th>Basic Viewer</th><th>CA13 (18 billion Points)</th><th>Retz (Potree + Cesium)</th><th>Classifications</th><th>Various Features</th><th>Toolbar</th>
	</tr>
</table>

<details>
<summary>More Examples</summary>


<table>
	<tr>
		<td>
			<a href="http://potree.org/potree/examples/load_project.html" target="_blank">
				<img src="examples/thumbnails/load_project.jpg" width="100%" />
			</a>
		</td><td>
			<a href="http://potree.org/potree/examples/matcap.html" target="_blank">
				<img src="examples/thumbnails/matcap.jpg" width="100%" />
			</a>
		</td><td>
			<a href="https://potree.org/potree/examples/vr_heidentor.html" target="_blank">
				<img src="examples/thumbnails/heidentor.jpg" width="100%" />
			</a>
		</td><td>
			<a href="http://potree.org/potree/examples/heidentor.html" target="_blank">
				<img src="examples/thumbnails/heidentor.png" width="100%" />
			</a>
		</td><td>
			<a href="http://potree.org/potree/examples/lion.html" target="_blank">
				<img src="examples/thumbnails/lion.png" width="100%" />
			</a>
		</td><td>
			<a href="http://potree.org/potree/examples/lion_las.html" target="_blank">
				<img src="examples/thumbnails/lion_las.png" width="100%" />
			</a>
		</td>
	</tr>
	<tr>
		<th>Load Project</th><th>Matcap</th><th>Virtual Reality</th><th>Heidentor</th><th>Lion</th><th>Lion LAS</th>
	</tr><tr>
		<td>
			<a href="http://potree.org/potree/examples/lion_laz.html" target="_blank">
				<img src="examples/thumbnails/lion_las.png" width="100%" />
			</a>
		</td><td>
			<a href="http://potree.org/potree/examples/ept.html" target="_blank">
				<img src="examples/thumbnails/lion.png" width="100%" />
			</a>
		</td><td>
			<a href="http://potree.org/potree/examples/ept_binary.html" target="_blank">
				<img src="examples/thumbnails/lion_las.png" width="100%" />
			</a>
		</td><td>
			<a href="http://potree.org/potree/examples/ept_zstandard.html" target="_blank">
				<img src="examples/thumbnails/lion_las.png" width="100%" />
			</a>
		</td><td>
			<a href="http://potree.org/potree/examples/clipping_volume.html" target="_blank">
				<img src="examples/thumbnails/clipping_volume.png" width="100%" />
			</a>
		</td><td>
			<a href="http://potree.org/potree/examples/oriented_images.html" target="_blank">
				<img src="examples/thumbnails/oriented_images.jpg" width="100%" />
			</a>
		</td>
	</tr>
	<tr>
		<th>Lion LAZ</th><th>EPT</th><th>EPT Binary</th><th>EPT zstandard</th><th>Clipping Volume</th><th>Oriented Images</th>
	</tr><tr>
		<td>
			<a href="http://potree.org/potree/examples/elevation_profile.html" target="_blank">
				<img src="examples/thumbnails/elevation_profile.png" width="100%" />
			</a>
		</td><td>
			<a href="http://potree.org/potree/examples/measurements.html" target="_blank">
				<img src="examples/thumbnails/measurements.png" width="100%" />
			</a>
		</td><td>
			<a href="http://potree.org/potree/examples/meshes.html" target="_blank">
				<img src="examples/thumbnails/meshes.png" width="100%" />
			</a>
		</td><td>
			<a href="http://potree.org/potree/examples/multiple_pointclouds.html" target="_blank">
				<img src="examples/thumbnails/multiple_point_clouds.png" width="100%" />
			</a>
		</td><td>
			<a href="http://potree.org/potree/examples/camera_animation.html" target="_blank">
				<img src="examples/thumbnails/camera_animation.jpg" width="100%" />
			</a>
		</td><td>
			<a href="http://potree.org/potree/examples/features_ca13.html" target="_blank">
				<img src="examples/thumbnails/features_ca13.png" width="100%" />
			</a>
		</td>
	</tr>
	<tr>
		<th>Elevation Profile</th><th>Measurements</th><th>Meshes</th><th>Multiple Point Clouds</th><th>Camera Animation</th><th>Features (CA13)</th>
	</tr><tr>
		<td>
			<a href="http://potree.org/potree/examples/annotations.html" target="_blank">
				<img src="examples/thumbnails/annotations.png" width="100%" />
			</a>
		</td><td>
			<a href="http://potree.org/potree/examples/annotation_hierarchy.html" target="_blank">
				<img src="examples/thumbnails/annotation_hierarchy.png" width="100%" />
			</a>
		</td><td>
			<a href="http://potree.org/potree/examples/animation_paths.html" target="_blank">
				<img src="examples/thumbnails/animation_paths.png" width="100%" />
			</a>
		</td><td>
			<a href="http://potree.org/potree/examples/shapefiles.html" target="_blank">
				<img src="examples/thumbnails/shapefiles.png" width="100%" />
			</a>
		</td><td>
			<a href="http://potree.org/potree/examples/cesium_ca13.html" target="_blank">
				<img src="examples/thumbnails/cesium_ca13.png" width="100%" />
			</a>
		</td><td>
			<a href="http://potree.org/potree/examples/geopackage.html" target="_blank">
				<img src="examples/thumbnails/geopackage.jpg" width="100%" />
			</a>
		</td>
	</tr>
	<tr>
		<th>Annotations</th><th>Hierarchical Annotations</th><th>Animation Path</th><th>Shapefiles</th><th>Cesium CA13</th><th>Geopackage</th>
	</tr><tr>
		<td>
			<a href="http://potree.org/potree/examples/cesium_sorvilier.html" target="_blank">
				<img src="examples/thumbnails/cesium_sorvilier.png" width="100%" />
			</a>
		</td><td>
			<a href="http://potree.org/potree/examples/custom_sidebar_section.html" target="_blank">
				<img src="examples/thumbnails/custom_sidebar_section.png" width="100%" />
			</a>
		</td><td>
			<a href="http://potree.org/potree/examples/embedded_iframe.html" target="_blank">
				<img src="examples/thumbnails/embedded_iframe.png" width="100%" />
			</a>
		</td><td>
			<a href="http://potree.org/potree/examples/gradient_colors.html" target="_blank">
				<img src="examples/thumbnails/gradient_colors.png" width="100%" />
			</a>
		</td>
	</tr>
	<tr>
		<th>Cesium Sorvilier</th><th>Custom Sidebar Section</th><th>Embedded Iframe</th><th>Gradient Colors</th>
	</tr>
</table>
</details>

# VR

<table>
	<tr>
		<td>
			<a href="https://potree.org/potree/examples/vr_heidentor.html" target="_blank">
				<img src="examples/thumbnails/heidentor.jpg" width="100%" />
			</a>
		</td><td>
			<a href="https://potree.org/potree/examples/vr_eclepens.html" target="_blank">
				<img src="examples/thumbnails/eclepens.jpg" width="100%" />
			</a>
		</td><td>
			<a href="https://potree.org/potree/examples/vr_morro_bay.html" target="_blank">
				<img src="examples/thumbnails/ca13.png" width="100%" />
			</a>
		</td><td>
			<a href="https://potree.org/potree/examples/vr_lion.html" target="_blank">
				<img src="examples/thumbnails/lion.png" width="100%" />
			</a>
		</td><td>
			<a href="https://potree.org/potree/examples/vr_dechen_cave.html" target="_blank">
				<img src="examples/thumbnails/dechen_cave.jpg" width="100%" />
			</a>
		</td>
	</tr>
	<tr>
		<th>Heidentor</th><th>Eclepens</th><th>Morro Bay</th><th>Lion</th><th>Dechen Cave</th>
	</tr>
</table>

# Showcase

<table>
	<tr>
		<td>
			<a href="http://potree.org/potree/examples/showcase/matterhorn.html" target="_blank">
				<img src="examples/thumbnails/matterhorn.jpg" width="100%" />
			</a>
		</td><td>
			<a href="http://potree.org/potree/examples/showcase/retz.html" target="_blank">
				<img src="examples/thumbnails/retz.jpg" width="100%" />
			</a>
		</td><td>
			<a href="http://potree.org/potree/examples/showcase/lake_tahoe.html" target="_blank">
				<img src="examples/thumbnails/lake_tahoe.jpg" width="100%" />
			</a>
		</td><td>
			<a href="http://potree.org/potree/examples/showcase/sorvilier.html" target="_blank">
				<img src="examples/thumbnails/vol_total.jpg" width="100%" />
			</a>
		</td><td>
			<a href="http://potree.org/potree/examples/showcase/grab_15.html" target="_blank">
				<img src="examples/thumbnails/grab_15.jpg" width="100%" />
			</a>
		</td><td>
			<a href="http://potree.org/potree/examples/showcase/tern_auscover_chowilla.html" target="_blank">
				<img src="examples/thumbnails/chowilla.jpg" width="100%" />
			</a>
		</td>
	</tr>
	<tr>
		<th>Matterhorn</th><th>Retz</th><th>Lake Tahoe</th><th>Sorvilier</th><th>Grave</th><th>Chowilla</th>
	</tr>
</table>

<details>
<summary>More</summary>

<table>
	<tr>
		<td>
			<a href="http://potree.org/potree/examples/showcase/chiller.html" target="_blank">
				<img src="examples/thumbnails/chiller.jpg" width="100%" />
			</a>
		</td><td>
			<a href="http://potree.org/potree/examples/showcase/cooler_tower.html" target="_blank">
				<img src="examples/thumbnails/cooler_tower.jpg" width="100%" />
			</a>
		</td><td>
			<a href="http://potree.org/potree/examples/showcase/dechen_cave.html" target="_blank">
				<img src="examples/thumbnails/dechen_cave.jpg" width="100%" />
			</a>
		</td><td>
			<a href="http://potree.org/potree/examples/showcase/doverMillRuins.html" target="_blank">
				<img src="examples/thumbnails/DoverMillRuins.jpg" width="100%" />
			</a>
		</td><td>
			<a href="http://potree.org/potree/examples/showcase/eclepens.html" target="_blank">
				<img src="examples/thumbnails/eclepens.jpg" width="100%" />
			</a>
		</td><td>
			<a href="http://potree.org/potree/examples/showcase/heidentor.html" target="_blank">
				<img src="examples/thumbnails/heidentor.jpg" width="100%" />
			</a>
		</td>
	</tr>
	<tr>
		<th>Chiller</th><th>Cooler</th><th>Dechen Cave</th><th>Ruins</th><th>Eclepens</th><th>Heidentor</th>
	</tr><tr>
		<td>
			<a href="http://potree.org/potree/examples/showcase/land_building.html" target="_blank">
				<img src="examples/thumbnails/land_building.jpg" width="100%" />
			</a>
		</td><td>
			<a href="http://potree.org/potree/examples/showcase/LDHI_module.html" target="_blank">
				<img src="examples/thumbnails/LDHI_module.jpg" width="100%" />
			</a>
		</td><td>
			<a href="http://potree.org/potree/examples/showcase/lion_head_simone_garagnani.html" target="_blank">
				<img src="examples/thumbnails/lion_head.jpg" width="100%" />
			</a>
		</td><td>
			<a href="http://potree.org/potree/examples/showcase/overpass.html" target="_blank">
				<img src="examples/thumbnails/overpass.jpg" width="100%" />
			</a>
		</td><td>
			<a href="http://potree.org/potree/examples/showcase/pielach.html" target="_blank">
				<img src="examples/thumbnails/pielach.jpg" width="100%" />
			</a>
		</td><td>
			<a href="http://potree.org/potree/examples/showcase/pompei.html" target="_blank">
				<img src="examples/thumbnails/pompei.jpg" width="100%" />
			</a>
		</td>
	</tr>
	<tr>
		<th>Building</th><th>LDHI</th><th>Lion Head</th><th>Overpass</th><th>Pielach</th><th>pompei</th>
	</tr><tr>
		<td>
			<a href="http://potree.org/potree/examples/showcase/santorini.html" target="_blank">
				<img src="examples/thumbnails/santorini.jpg" width="100%" />
			</a>
		</td><td>
			<a href="http://potree.org/potree/examples/showcase/skatepark.html" target="_blank">
				<img src="examples/thumbnails/skatepark.jpg" width="100%" />
			</a>
		</td><td>
			<a href="http://potree.org/potree/examples/showcase/subsea_equipment.html" target="_blank">
				<img src="examples/thumbnails/subsea_equipment.jpg" width="100%" />
			</a>
		</td><td>
			<a href="http://potree.org/potree/examples/showcase/subsea_manifold.html" target="_blank">
				<img src="examples/thumbnails/subseamanifold.jpg" width="100%" />
			</a>
		</td><td>
			<a href="http://potree.org/potree/examples/showcase/westend_palais.html" target="_blank">
				<img src="examples/thumbnails/westend_palais.jpg" width="100%" />
			</a>
		</td><td>
			<a href="http://potree.org/potree/examples/showcase/whitby.html" target="_blank">
				<img src="examples/thumbnails/whitby.jpg" width="100%" />
			</a>
		</td>
	</tr>
	<tr>
		<th>Santorini</th><th>Skatepark</th><th>Subsea Eq.</th><th>Subsea Man.</th><th>Westend Palais</th><th>Whitby</th>
	</tr>
</table>

</details>

# Funding

Potree is funded by a combination of research projects, companies and institutions. 

Research projects who's funding contributes to Potree:

<table>
	<tr>
		<th>Project Name</th>
		<th>Funding Agency</th>
	</tr>
	<tr>
		<td><a href="https://projekte.ffg.at/projekt/3851914">LargeClouds2BIM</a></td>
		<td><a href="https://www.ffg.at/">FFG</a></td>
	</tr>
	<tr>
		<td><a href="https://harvest4d.org/">Harvest4D</a></td>
		<td><a href="https://ec.europa.eu/transport/themes/research/fp7_en">EU 7th Framework Program 323567</a></td>
	</tr>
	<tr>
		<td><a href="https://gcd.tuwien.ac.at/">GCD Doctoral College</a></td>
		<td><a href="https://www.tuwien.at/en/">TU Wien</a></td>
	</tr>
	<tr>
		<td><a href="https://www.cg.tuwien.ac.at/research/projects/Superhumans/">Superhumans</a></td>
		<td><a href="https://www.fwf.ac.at/">FWF</a></td>
	</tr>
</table>

We would like to thank our sponsors for their financial contributions that keep this project up and running!

<table>
	<tr>
		<th>
			Diamond<br>
			€ 15,000+
		</th>
		<td>
			<a href="http://www.ne.ch/autorites/DDTE/SGRF/SITN/Pages/accueil.aspx">
				<img src="docs/sponsors/sitn_logo.png" height="80px"/> &nbsp;
			</a> &nbsp;
			<a href="http://www.synth3d.co">
				<img src="docs/sponsors/synth.png" height="120"/>
			</a> &nbsp;
			<a href="http://www.geocue.com">
				<img src="docs/sponsors/geocue.png" height="120px"/>
			</a> &nbsp;
			<a href="http://rapidlasso.com">
				<img src="./docs/sponsors/rapidlasso_square_256x2561.png" width="150" height="150"/>
			</a> &nbsp;
		</td>
	</tr>
	<tr>
		<th>
			Gold<br>
			€ 10,000+
		</th>
		<td>
			<a href="https://www.bart.gov">
				<img src="docs/sponsors/bart.png" height="100"/>
			</a>
		</td>
	</tr>
	<tr>
		<th>
			Silver<br>
			€ 5,000+
		</th>
		<td>
			<a href="https://biology.anu.edu.au/research/facilities/australian-plant-phenomics-facility-anu">
				<img src="docs/sponsors/APPF full logo.png" height="70"/> &nbsp;
			</a>
			<a href="https://www.limit-addict.fr/">
				<img src="docs/sponsors/limitaddict.png" height="45"/>
			</a>
			<a href="http://georepublic.info">
				<img src="docs/sponsors/georepublic.png" height="45"/>
			</a>
		</td>
	</tr>
	<tr>
		<th>
			Bronze<br>
			€ 1,000+
		</th>
		<td>
			<a href="https://www.eventart.at/">
				<img src="docs/sponsors/eventart.png" height="55"/> &nbsp;
			</a>
			<a href="https://www.geodelta.com/">
				<img src="docs/sponsors/geodelta.png" height="35"/> &nbsp;
			</a>
			<a href="https://www.e-cassini.fr/">
				<img src="docs/sponsors/e_cassini.jpg" height="70"/> &nbsp;
			</a>
			<a href="https://www.sogelink.fr/">
				<img src="docs/sponsors/SOGELINK_SO-EASY.png" height="40"/> &nbsp;
			</a>
			<b>Data-viewer</b>
			<a href="http://www.helimap.com/">
				<img src="docs/sponsors/helimap.gif" height="60"/> &nbsp;
			</a>
			<a href="http://www.vevey.ch/">
				<img src="docs/sponsors/vevey.png" height="60"/> &nbsp;
			</a>
			<a href="https://www.yverdon-les-bains.ch/">
				<img src="docs/sponsors/Logo-YLB.png" height="60"/> &nbsp;
			</a>
			<a href="http://archpro.lbg.ac.at">
				<img src="docs/sponsors/archpro_EN_small.png" height="60"/> 
			</a> &nbsp;
			<br>
			<a href="http://www.kts.co.jp">
				<img src="docs/sponsors/kts.png" height="32"/> &nbsp;
			</a>
			<a href="http://veesus.com">
				<img src="docs/sponsors/veesus_small.png" height="40"/> &nbsp;
			</a>
			<a href="http://www.sigeom.ch">
				<img src="docs/sponsors/logo_sigeom.png" height="40"/> &nbsp;
			</a>
		</td>
	</tr>
</table>



# Credits

* The multi-res-octree algorithms used by this viewer were developed at the Vienna University of Technology by Michael Wimmer and Claus Scheiblauer as part of the [Scanopy Project](http://www.cg.tuwien.ac.at/research/projects/Scanopy/).
* [Three.js](https://github.com/mrdoob/three.js), the WebGL 3D rendering library on which potree is built.
* [plas.io](http://plas.io/) point cloud viewer. LAS and LAZ support have been taken from the laslaz.js implementation of plas.io. Thanks to [Uday Verma](https://twitter.com/udaykverma) and [Howard Butler](https://twitter.com/howardbutler) for this!
* [Harvest4D](https://harvest4d.org/) Potree currently runs as Master Thesis under the Harvest4D Project
* Christian Boucheny (EDL developer) and Daniel Girardeau-Montaut ([CloudCompare](http://www.danielgm.net/cc/)). The EDL shader was adapted from the CloudCompare source code!
* [Martin Isenburg](http://rapidlasso.com/), [Georepublic](http://georepublic.de/en/),
[Veesus](http://veesus.com/), [Sigeom Sa](http://www.sigeom.ch/), [SITN](http://www.ne.ch/sitn), [LBI ArchPro](http://archpro.lbg.ac.at/),  [Pix4D](http://pix4d.com/) as well as all the contributers to potree and PotreeConverter and many more for their support.
three.js
========

[![NPM package][npm]][npm-url]
[![Build Size][build-size]][build-size-url]
[![Build Status][build-status]][build-status-url]
[![Dependencies][dependencies]][dependencies-url]
[![Dev Dependencies][dev-dependencies]][dev-dependencies-url]
[![Language Grade][lgtm]][lgtm-url]

#### JavaScript 3D library ####

The aim of the project is to create an easy to use, lightweight, 3D library with a default WebGL renderer. The library also provides Canvas 2D, SVG and CSS3D renderers in the examples.

[Examples](http://threejs.org/examples/) &mdash;
[Documentation](http://threejs.org/docs/) &mdash;
[Wiki](https://github.com/mrdoob/three.js/wiki) &mdash;
[Migrating](https://github.com/mrdoob/three.js/wiki/Migration-Guide) &mdash;
[Questions](http://stackoverflow.com/questions/tagged/three.js) &mdash;
[Forum](https://discourse.threejs.org/) &mdash;
[Gitter](https://gitter.im/mrdoob/three.js) &mdash;
[Slack](https://join.slack.com/t/threejs/shared_invite/enQtMzYxMzczODM2OTgxLTQ1YmY4YTQxOTFjNDAzYmQ4NjU2YzRhNzliY2RiNDEyYjU2MjhhODgyYWQ5Y2MyZTU3MWNkOGVmOGRhOTQzYTk)

### Usage ###

Download the [minified library](http://threejs.org/build/three.min.js) and include it in your HTML, or install and import it as a [module](http://threejs.org/docs/#manual/introduction/Import-via-modules),
Alternatively, see [how to build the library yourself](https://github.com/mrdoob/three.js/wiki/Build-instructions).

```html
<script src="js/three.min.js"></script>
```

This code creates a scene, a camera, and a geometric cube, and it adds the cube to the scene. It then creates a `WebGL` renderer for the scene and camera, and it adds that viewport to the `document.body` element. Finally, it animates the cube within the scene for the camera.

```javascript
var camera, scene, renderer;
var geometry, material, mesh;

init();
animate();

function init() {

	camera = new THREE.PerspectiveCamera( 70, window.innerWidth / window.innerHeight, 0.01, 10 );
	camera.position.z = 1;

	scene = new THREE.Scene();

	geometry = new THREE.BoxGeometry( 0.2, 0.2, 0.2 );
	material = new THREE.MeshNormalMaterial();

	mesh = new THREE.Mesh( geometry, material );
	scene.add( mesh );

	renderer = new THREE.WebGLRenderer( { antialias: true } );
	renderer.setSize( window.innerWidth, window.innerHeight );
	document.body.appendChild( renderer.domElement );

}

function animate() {

	requestAnimationFrame( animate );

	mesh.rotation.x += 0.01;
	mesh.rotation.y += 0.02;

	renderer.render( scene, camera );

}
```

If everything went well you should see [this](https://jsfiddle.net/f2Lommf5/).

### Change log ###

[Releases](https://github.com/mrdoob/three.js/releases)


[npm]: https://img.shields.io/npm/v/three.svg
[npm-url]: https://www.npmjs.com/package/three
[build-size]: https://badgen.net/bundlephobia/minzip/three
[build-size-url]: https://bundlephobia.com/result?p=three
[build-status]: https://travis-ci.org/mrdoob/three.js.svg?branch=dev
[build-status-url]: https://travis-ci.org/mrdoob/three.js
[dependencies]: https://img.shields.io/david/mrdoob/three.js.svg
[dependencies-url]: https://david-dm.org/mrdoob/three.js
[dev-dependencies]: https://img.shields.io/david/dev/mrdoob/three.js.svg
[dev-dependencies-url]: https://david-dm.org/mrdoob/three.js#info=devDependencies
[lgtm]: https://img.shields.io/lgtm/grade/javascript/g/mrdoob/three.js.svg?label=code%20quality
[lgtm-url]: https://lgtm.com/projects/g/mrdoob/three.js/
# JSON5 – JSON for Humans

[![Build Status](https://travis-ci.org/json5/json5.svg)][Build Status]
[![Coverage
Status](https://coveralls.io/repos/github/json5/json5/badge.svg)][Coverage
Status]

The JSON5 Data Interchange Format (JSON5) is a superset of [JSON] that aims to
alleviate some of the limitations of JSON by expanding its syntax to include
some productions from [ECMAScript 5.1].

This JavaScript library is the official reference implementation for JSON5
parsing and serialization libraries.

[Build Status]: https://travis-ci.org/json5/json5

[Coverage Status]: https://coveralls.io/github/json5/json5

[JSON]: https://tools.ietf.org/html/rfc7159

[ECMAScript 5.1]: https://www.ecma-international.org/ecma-262/5.1/

## Summary of Features
The following ECMAScript 5.1 features, which are not supported in JSON, have
been extended to JSON5.

### Objects
- Object keys may be an ECMAScript 5.1 _[IdentifierName]_.
- Objects may have a single trailing comma.

### Arrays
- Arrays may have a single trailing comma.

### Strings
- Strings may be single quoted.
- Strings may span multiple lines by escaping new line characters.
- Strings may include character escapes.

### Numbers
- Numbers may be hexadecimal.
- Numbers may have a leading or trailing decimal point.
- Numbers may be [IEEE 754] positive infinity, negative infinity, and NaN.
- Numbers may begin with an explicit plus sign.

### Comments
- Single and multi-line comments are allowed.

### White Space
- Additional white space characters are allowed.

[IdentifierName]: https://www.ecma-international.org/ecma-262/5.1/#sec-7.6

[IEEE 754]: http://ieeexplore.ieee.org/servlet/opac?punumber=4610933

## Short Example
```js
{
  // comments
  unquoted: 'and you can quote me on that',
  singleQuotes: 'I can use "double quotes" here',
  lineBreaks: "Look, Mom! \
No \\n's!",
  hexadecimal: 0xdecaf,
  leadingDecimalPoint: .8675309, andTrailing: 8675309.,
  positiveSign: +1,
  trailingComma: 'in objects', andIn: ['arrays',],
  "backwardsCompatible": "with JSON",
}
```

## Specification
For a detailed explanation of the JSON5 format, please read the [official
specification](https://json5.github.io/json5-spec/).

## Installation
### Node.js
```sh
npm install json5
```

```js
const JSON5 = require('json5')
```

### Browsers
```html
<script src="https://unpkg.com/json5@^2.0.0/dist/index.min.js"></script>
```

This will create a global `JSON5` variable.

## API
The JSON5 API is compatible with the [JSON API].

[JSON API]:
https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Global_Objects/JSON

### JSON5.parse()
Parses a JSON5 string, constructing the JavaScript value or object described by
the string. An optional reviver function can be provided to perform a
transformation on the resulting object before it is returned.

#### Syntax
    JSON5.parse(text[, reviver])

#### Parameters
- `text`: The string to parse as JSON5.
- `reviver`: If a function, this prescribes how the value originally produced by
  parsing is transformed, before being returned.

#### Return value
The object corresponding to the given JSON5 text.

### JSON5.stringify()
Converts a JavaScript value to a JSON5 string, optionally replacing values if a
replacer function is specified, or optionally including only the specified
properties if a replacer array is specified.

#### Syntax
    JSON5.stringify(value[, replacer[, space]])
    JSON5.stringify(value[, options])

#### Parameters
- `value`: The value to convert to a JSON5 string.
- `replacer`: A function that alters the behavior of the stringification
  process, or an array of String and Number objects that serve as a whitelist
  for selecting/filtering the properties of the value object to be included in
  the JSON5 string. If this value is null or not provided, all properties of the
  object are included in the resulting JSON5 string.
- `space`: A String or Number object that's used to insert white space into the
  output JSON5 string for readability purposes. If this is a Number, it
  indicates the number of space characters to use as white space; this number is
  capped at 10 (if it is greater, the value is just 10). Values less than 1
  indicate that no space should be used. If this is a String, the string (or the
  first 10 characters of the string, if it's longer than that) is used as white
  space. If this parameter is not provided (or is null), no white space is used.
  If white space is used, trailing commas will be used in objects and arrays.
- `options`: An object with the following properties:
  - `replacer`: Same as the `replacer` parameter.
  - `space`: Same as the `space` parameter.
  - `quote`: A String representing the quote character to use when serializing
    strings.

#### Return value
A JSON5 string representing the value.

### Node.js `require()` JSON5 files
When using Node.js, you can `require()` JSON5 files by adding the following
statement.

```js
require('json5/lib/register')
```

Then you can load a JSON5 file with a Node.js `require()` statement. For
example:

```js
const config = require('./config.json5')
```

## CLI
Since JSON is more widely used than JSON5, this package includes a CLI for
converting JSON5 to JSON and for validating the syntax of JSON5 documents.

### Installation
```sh
npm install --global json5
```

### Usage
```sh
json5 [options] <file>
```

If `<file>` is not provided, then STDIN is used.

#### Options:
- `-s`, `--space`: The number of spaces to indent or `t` for tabs
- `-o`, `--out-file [file]`: Output to the specified file, otherwise STDOUT
- `-v`, `--validate`: Validate JSON5 but do not output JSON
- `-V`, `--version`: Output the version number
- `-h`, `--help`: Output usage information

## Contributing
### Development
```sh
git clone https://github.com/json5/json5
cd json5
npm install
```

When contributing code, please write relevant tests and run `npm test` and `npm
run lint` before submitting pull requests. Please use an editor that supports
[EditorConfig](http://editorconfig.org/).

### Issues
To report bugs or request features regarding the JSON5 data format, please
submit an issue to the [official specification
repository](https://github.com/json5/json5-spec).

To report bugs or request features regarding the JavaScript implementation of
JSON5, please submit an issue to this repository.

## License
MIT. See [LICENSE.md](./LICENSE.md) for details.

## Credits
[Assem Kishore](https://github.com/aseemk) founded this project.

[Michael Bolin](http://bolinfest.com/) independently arrived at and published
some of these same ideas with awesome explanations and detail. Recommended
reading: [Suggested Improvements to JSON](http://bolinfest.com/essays/json.html)

[Douglas Crockford](http://www.crockford.com/) of course designed and built
JSON, but his state machine diagrams on the [JSON website](http://json.org/), as
cheesy as it may sound, gave us motivation and confidence that building a new
parser to implement these ideas was within reach! The original
implementation of JSON5 was also modeled directly off of Doug’s open-source
[json_parse.js] parser. We’re grateful for that clean and well-documented
code.

[json_parse.js]:
https://github.com/douglascrockford/JSON-js/blob/03157639c7a7cddd2e9f032537f346f1a87c0f6d/json_parse.js

[Max Nanasy](https://github.com/MaxNanasy) has been an early and prolific
supporter, contributing multiple patches and ideas.

[Andrew Eisenberg](https://github.com/aeisenberg) contributed the original
`stringify` method.

[Jordan Tucker](https://github.com/jordanbtucker) has aligned JSON5 more closely
with ES5, wrote the official JSON5 specification, completely rewrote the
codebase from the ground up, and is actively maintaining this project.
# jstree

[jsTree](http://www.jstree.com/) is jquery plugin, that provides interactive trees. It is absolutely free, [open source](https://github.com/vakata/jstree) and distributed under the MIT license.

jsTree is easily extendable, themable and configurable, it supports HTML & JSON data sources, AJAX & async callback loading.

jsTree functions properly in either box-model (content-box or border-box), can be loaded as an AMD module, and has a built in mobile theme for responsive design, that can easily be customized. It uses jQuery's event system, so binding callbacks on various events in the tree is familiar and easy.

You also get:
 * drag & drop support
 * keyboard navigation
 * inline edit, create and delete
 * tri-state checkboxes
 * fuzzy searching
 * customizable node types

_Aside from this readme you can find a lot more info on [jstree.com](http://www.jstree.com) & [the discussion group](https://groups.google.com/forum/#!forum/jstree)_.

---

<!-- MarkdownTOC depth=0 autolink=true bracket=round -->

- [Getting Started](#getting-started)
  - [Include all neccessary files](#include-all-neccessary-files)
  - [Populating a tree using HTML](#populating-a-tree-using-html)
  - [Populating a tree using an array \(or JSON\)](#populating-a-tree-using-an-array-or-json)
    - [The required JSON format](#the-required-json-format)
  - [Populating the tree using AJAX](#populating-the-tree-using-ajax)
  - [Populating the tree using AJAX and lazy loading nodes](#populating-the-tree-using-ajax-and-lazy-loading-nodes)
  - [Populating the tree using a callback function](#populating-the-tree-using-a-callback-function)
- [Working with events](#working-with-events)
- [Interacting with the tree using the API](#interacting-with-the-tree-using-the-api)
- [More on configuration](#more-on-configuration)
- [Plugins](#plugins)
  - [checkbox](#checkbox)
  - [contextmenu](#contextmenu)
  - [dnd](#dnd)
  - [massload](#massload)
  - [search](#search)
  - [sort](#sort)
  - [state](#state)
  - [types](#types)
  - [unique](#unique)
  - [wholerow](#wholerow)
  - [More plugins](#more-plugins)
- [PHP demos moved to new repository](#php-demos-moved-to-new-repository)
- [License & Contributing](#license--contributing)

<!-- /MarkdownTOC -->


---

## Getting Started

### Include all neccessary files
To get started you need 3 things in your page:
 1. jQuery (anything above 1.9.1 will work)
 2. A jstree theme (there is only one theme supplied by default)
 3. The jstree source file

```html
<script src="//cdnjs.cloudflare.com/ajax/libs/jquery/3.1.1/jquery.min.js"></script>

<link rel="stylesheet" href="//cdnjs.cloudflare.com/ajax/libs/jstree/3.3.8/themes/default/style.min.css" />
<script src="//cdnjs.cloudflare.com/ajax/libs/jstree/3.3.8/jstree.min.js"></script>
```

_If you decide to host jstree yourself - the files are located in the `dist` folder. You can safely ignore the `dist/libs` folder._

---

### Populating a tree using HTML

Now we are all set to create a tree, inline HTML is the easiest option (suitable for menus). All you need to do is select a node (using a jQuery selector) and invoke the `.jstree()` function to let jstree know you want to render a tree inside the selected node. `$.jstree.create(element)` can be used too.

```html
<div id="container">
  <ul>
    <li>Root node
      <ul>
        <li>Child node 1</li>
        <li>Child node 2</li>
      </ul>
    </li>
  </ul>
</div>
<script>
$(function() {
  $('#container').jstree();
});
</script>
```

[view result](http://jsfiddle.net/vakata/2kwkh2uL/)

_You can add a few options when rendering a node using a data-attribute (note the quotes):_
```html
<li data-jstree='{ "selected" : true, "opened" : true }'>Root node ...
```

---

### Populating a tree using an array (or JSON)

Building trees from HTML is easy, but it is not very flexible, inline JS data is a better option:

```html
<div id="container"></div>
<script>
$(function() {
  $('#container').jstree({
    'core' : {
      'data' : [
        { "text" : "Root node", "children" : [
            { "text" : "Child node 1" },
            { "text" : "Child node 2" }
          ]
        }
      ]
    }
  });
});
</script>
```

[view result](http://jsfiddle.net/vakata/2kwkh2uL/4478/)

Unlike the previous simple HTML example, this time the `.jstree()` function accepts a config object.

For now it is important to note that jstree will try to parse any data you specify in the  `core.data` key and use it to create a tree. As seen in the previous example, if this key is missing jstree will try to parse the inline HTML of the container.

#### The required JSON format

The data you use must be in a specific format, each branch of the tree is represented by an object, which must at least have a `text` key. The `children` key can be used to add children to the branch, it should be an array of objects.

_Keep in mind, you can use a simple string instead of an object if all you need is node with the given text, the above data can be written as:_

```js
[ { "text" : "Root node", "children" : [ "Child node 1", "Child node 2" ] } ]
```

There are other available options for each node, only set them if you need them like:

 * `id` - makes if possible to identify a node later (will also be used as a DOM ID of the `LI` node). _Make sure you do not repeat the same ID in a tree instance (that would defeat its purpose of being a unique identifier and may cause problems for jstree)_.
 * `icon` - a string which will be used for the node's icon - this can either be a path to a file, or a className (or list of classNames), which you can style in your CSS (font icons also work).
 * `data` - this can be anything you want - it is metadata you want attached to the node - you will be able to access and modify it any time later - it has no effect on the visuals of the node.
 * `state` - an object specifyng a few options about the node:
   - `selected` - if the node should be initially selected
   - `opened` - if the node should be initially opened
   - `disabled` - if the node should be disabled
   - `checked` - __checkbox plugin specific__ - if the node should be checked (only used when `tie_selection` is `false`, which you should only do if you really know what you are doing)
   - `undetermined` - __checkbox plugin specific__ - if the node should be rendered in undetermined state (only used with lazy loading and when the node is not yet loaded, otherwise this state is automatically calculated).
 * `type` - __types plugin specific__ - the type of the nodes (should be defined in the types config), if not set `"default"` is assumed.
 * `li_attr` - object of values which will be used to add HTML attributes on the resulting `LI` DOM node.
 * `a_attr` - object of values which will be used to add HTML attributes on the resulting `A` node.

Here is a new demo with some of those properties set:

```html
<div id="container"></div>
<script>
$(function() {
  $('#container').jstree({
    'core' : {
      'data' : [
          {
              "text" : "Root node",
              "state" : {"opened" : true },
              "children" : [
                  {
                    "text" : "Child node 1",
                    "state" : { "selected" : true },
                    "icon" : "glyphicon glyphicon-flash"
                  },
                  { "text" : "Child node 2", "state" : { "disabled" : true } }
              ]
        }
      ]
    }
  });
});
</script>
```

[view result](http://jsfiddle.net/vakata/2kwkh2uL/4479/)

---

### Populating the tree using AJAX

Building off of the previous example, let's see how to have jstree make AJAX requests for you.

```html
<div id="container"></div>
<script>
$(function() {
  $('#container').jstree({
    'core' : {
      'data' : {
        "url" : "//www.jstree.com/fiddle/",
        "dataType" : "json" // needed only if you do not supply JSON headers
      }
    }
  });
});
</script>
```

The server response is:
```json
[{
  "id":1,"text":"Root node","children":[
    {"id":2,"text":"Child node 1"},
    {"id":3,"text":"Child node 2"}
  ]
}]
```

[view result](http://jsfiddle.net/vakata/2kwkh2uL/4480/)

Instead of a JS array, you can set `core.data` to a [jQuery AJAX config](http://api.jquery.com/jQuery.ajax/). 
jsTree will hit that URL, and provided you return properly formatted JSON it will be displayed.

_If you cannot provide proper JSON headers, set `core.data.dataType` to `"json"`._

The ids in the server response make it possible to identify nodes later (which we will see in the next few demos), but they are not required.

__WHEN USING IDS MAKE SURE THEY ARE UNIQUE INSIDE A PARTICULAR TREE__

---

### Populating the tree using AJAX and lazy loading nodes

Lazy loading means nodes will be loaded when they are needed. Imagine you have a huge amount of nodes you want to show, but loading them with a single request is way too much traffic. Lazy loading makes it possible to load nodes on the fly - jstree will perform AJAX requests as the user browses the tree.

Here we take our previous example, and lazy load the "Child node 1" node.

```html
<div id="container"></div>
<script>
$(function() {
  $('#container').jstree({
    'core' : {
      'data' : {
        "url" : "//www.jstree.com/fiddle/?lazy",
        "data" : function (node) {
          return { "id" : node.id };
        }
      }
    }
  });
});
</script>
```

The initial server response is:
```json
[{
  "id":1,"text":"Root node","children":[
    {"id":2,"text":"Child node 1","children":true},
    {"id":3,"text":"Child node 2"}
  ]
}]
```

[view result](http://jsfiddle.net/vakata/2kwkh2uL/4481/)

Now to focus on what is different. First off the `"data"` config option of the data object. If you check with jQuery, it is supposed to be a string or an object. But jstree makes it possible to set a function.

Each time jstree needs to make an AJAX call this function will be called and will receive a single parameter - the node that is being loaded. The return value of this function will be used as the actual `"data"` of the AJAX call. To understand better open up the demo and see the requests go off in the console.

You will notice that the first request goes off to:
`http://www.jstree.com/fiddle?lazy&id=#`
`#` is the special ID that the function receives when jstree needs to load the root nodes.

Now go ahead and open the root node - two children will be shown, but no request will be made - that is because we loaded those children along with the first request.

Onto the next difference - "Child node 1" appears closed - that is because in the data we supplied `true` as the `"children"` property of this node (you can see it in the server response). This special value indicated to jstree, that it has to lazy load the "Child node 1" node.

Proceed and open this node - you will see a next request fire off to:
`http://www.jstree.com/fiddle?lazy&id=2`
ID is set to `2` because the node being loaded has an ID of `2`, and we have configured jstree to send the node ID along with the AJAX request (the `data` function).

The server response is:
```json
["Child node 3","Child node 4"]
```

_You can also set `"url"` to a function and it works exactly as with `"data"` - each time a request has to be made, jstree will invoke your function and the request will go off to whatever you return in this function. This is useful when dealing with URLs like: `http://example.com/get_children/1`._

### Populating the tree using a callback function

Sometimes you may not want jsTree to make AJAX calls for you - you might want to make them yourself, or use some other method of populating the tree. In that case you can use a callback function.

```html
<div id="container"></div>
<script>
$(function() {
  $('#container').jstree({
    'core' : {
      'data' : function (node, cb) {
        if(node.id === "#") {
          cb([{"text" : "Root", "id" : "1", "children" : true}]);
        }
        else {
          cb(["Child"]);
        }
      }
    }
  });
});
</script>
```

[view result](http://jsfiddle.net/vakata/2kwkh2uL/4482/)

As you can see your function will receive two arguments - the node whose children need to be loaded and a callback function to call with the data once you have it. The data follows the same familiar JSON format and lazy loading works just as with AJAX (as you can see in the above example).

---

## Working with events

jstree provides a lot of events to let you know something happened with the tree. The events are the same regardless of how you populate the tree.
Let's use the most basic event `changed` - it fires when selection on the tree changes:

```html
<div id="container"></div>
<script>
$(function() {
  $('#container').jstree({
    'core' : {
      'data' : [
        {"id" : 1, "text" : "Node 1"},
        {"id" : 2, "text" : "Node 2"},
      ]
    }
  });
  $('#container').on("changed.jstree", function (e, data) {
    console.log("The selected nodes are:");
    console.log(data.selected);
  });
});
</script>
```

[view result](http://jsfiddle.net/vakata/2kwkh2uL/4483/)

All jstree events fire in a special `".jstree"` namespace - this is why we listen for `"changed.jstree"`. The handler itself receives one additional parameter - it will be populated with all you need to know about the event that happened. In this case `data.selected` is an array of selected node IDs (please note, that if you have not specified IDs they will be autogenerated).

Let's extend this a bit and log out the text of the node instead of the ID.

```js
$('#container').on("changed.jstree", function (e, data) {
  console.log(data.instance.get_selected(true)[0].text);
  console.log(data.instance.get_node(data.selected[0]).text);
});
```

The two rows above achieve exactly the same thing - get the text of the first selected node.

In the `data` argument object you will always get an `instance` key - that is a reference to the tree instance, so that you can easily invoke methods.

__All available functions and events are documented in the API docs__

---

## Interacting with the tree using the API

We scratched the surface on interacting with the tree in the previous example. Let's move on to obtaining an instance and calling a method on this instance:

```html
<button>Select node 1</button>
<div id="container"></div>
<script>
$(function() {
  $('#container').jstree({
    'core' : {
      'data' : [
        {"id" : 1, "text" : "Node 1"},
        {"id" : 2, "text" : "Node 2"},
      ]
    }
  });
  $('button').on("click", function () {
    var instance = $('#container').jstree(true);
    instance.deselect_all();
    instance.select_node('1');
  });
});
</script>
```

[view result](http://jsfiddle.net/vakata/2kwkh2uL/4484/)

The above example shows how to obtain a reference to a jstree instance (again with a selector, but this time instead of a config, we pass a boolean `true`), and call a couple of methods - the latter one is selecting a node by its ID.

Methods can also be invoked like this:

```js
$('#container').jstree("select_node", "1");
```

__All available functions and events are documented in the API docs__

## More on configuration

We already covered the config object in general (when we specified inline & AJAX data sources).

```js
$("#tree").jstree({ /* config object goes here */ });
```

Each key in the config object corresponds to a plugin, and the value of that key is the configuration for that plugin. There are also two special keys `"core"` and `"plugins"`:
 * `"core"` stores the core configuration options
 * `"plugins"` is an array of plugin names (strings) you want active on the instance

When configuring you only need to set values that you want to be different from the defaults.

__All config options and defaults are documented in the API docs__

```js
$("#tree").jstree({
  "core" : { // core options go here
    "multiple" : false, // no multiselection
    "themes" : {
      "dots" : false // no connecting dots between dots
    }
  },
  "plugins" : ["state"] // activate the state plugin on this instance
});
```

[view result](http://jsfiddle.net/vakata/2kwkh2uL/4485/)

We will cover all plugins further down.

__Keep in mind by default all modifications to the structure are prevented - that means drag'n'drop, create, rename, delete will not work unless you enable them.__

```js
$("#tree").jstree({
  "core" : {
    "check_callback" : true, // enable all modifications
  },
  "plugins" : ["dnd","contextmenu"]
});
```

[view result](http://jsfiddle.net/vakata/2kwkh2uL/4486/)

`"core.check_callback"` can also be set to a function, that will be invoked every time a modification is about to happen (or when jstree needs to check if a modification is possible). If you return `true` the operation will be allowed, a value of `false` means it will not be allowed. The possible operation you can expect are `create_node`, `rename_node`, `delete_node`, `move_node` and `copy_node`. The `more` parameter will contain various information provided by the plugin that is invoking the check. For example the DND plugin will provide an object containing information about the move or copy operation that is being checked - is it a multi tree operation, which node is currently hovered, where the insert arrow is pointing - before, after or inside, etc.

```js
$("#tree").jstree({
  "core" : {
    "check_callback" : function (operation, node, parent, position, more) {
      if(operation === "copy_node" || operation === "move_node") {
        if(parent.id === "#") {
          return false; // prevent moving a child above or below the root
        }
      },
      return true; // allow everything else
    }
  },
  "plugins" : ["dnd","contextmenu"]
});
```

[view result](http://jsfiddle.net/vakata/2kwkh2uL/4487/)

The `more` parameter you receive contains other information related to the check being performed.

__For example__: `move_node` & `copy_node` checks will fire repeatedly while the user drags a node, if the check was triggered by the `dnd` plugin `more` will contain a `dnd` key, which will be set to `true`.
You can check for `more.dnd` and only perform a certain action if `dnd` triggered the check.
If you only want to perform an operation when a node is really about to be dropped check for `more.core`.

## Plugins

jsTree comes with a few plugin bundled, but they will only modify your tree if you activate them using the `"plugins"` config option. Here is a brief description of each plugin. You can read more on the available config options for each plugin in the API docs.

### checkbox
Renders a checkbox icon in front of each node, making multiselection easy. It also has a "tri-state" option, meaning a node with some of its children checked will get a "square" icon.

_Keep in mind that if any sort of cascade is enabled, disabled nodes may be checked too (not by themselves, but for example when a parent of a disabled node is checked and selection is configured to cascade down)._

```js
$("#tree").jstree({
  "plugins" : ["checkbox"]
});
```

[view result](http://jsfiddle.net/vakata/2kwkh2uL/4488/)

### contextmenu
Makes it possible to right click nodes and shows a list of configurable actions in a menu.

```js
$("#tree").jstree({
  "core" : { "check_callback" : true }, // so that modifying operations work
  "plugins" : ["contextmenu"]
});
```

[view result](http://jsfiddle.net/vakata/2kwkh2uL/4489/)

### dnd
Makes it possible to drag and drop tree nodes and rearrange the tree.

```js
$("#tree").jstree({
  "core" : { "check_callback" : true }, // so that operations work
  "plugins" : ["dnd"]
});
```

[view result](http://jsfiddle.net/vakata/2kwkh2uL/4490/)

### massload
Makes it possible to load multiple nodes in a single go (for a lazy loaded tree).

```js
$("#tree").jstree({
  "core" : {
    "data" : { .. AJAX config .. }
  },
  "massload" : {
    "url" : "/some/path",
    "data" : function (nodes) {
      return { "ids" : nodes.join(",") };
    }
  },
  "plugins" : [ "massload", "state" ]
});
```

### search
Adds the possibility to search for items in the tree and show only matching nodes. It also has AJAX / callback hooks, so that search will work on lazy loaded trees too.

```html
<form id="s">
  <input type="search" id="q" />
  <button type="submit">Search</button>
</form>
<script>
$("#container").jstree({
  "plugins" : ["search"]
});
$("#s").submit(function(e) {
  e.preventDefault();
  $("#container").jstree(true).search($("#q").val());
});
</script>
```

[view result](http://jsfiddle.net/vakata/2kwkh2uL/4491/)

### sort
Automatically arranges all sibling nodes according to a comparison config option function, which defaults to alphabetical order.

```js
$("#tree").jstree({
  "plugins" : ["sort"]
});
```

[view result](http://jsfiddle.net/vakata/2kwkh2uL/4492/)

### state
Saves all opened and selected nodes in the user's browser, so when returning to the same tree the previous state will be restored.

```js
$("#tree").jstree({
  // the key is important if you have multiple trees in the same domain
  "state" : { "key" : "state_demo" },
  "plugins" : ["state"]
});
```

[view result](http://jsfiddle.net/vakata/2kwkh2uL/4493/)

### types
Makes it possible to add a "type" for a node, which means to easily control nesting rules and icon for groups of nodes instead of individually. To set a node type add a type property to the node structure.

```js
$("#tree").jstree({
  "types" : {
    "default" : {
      "icon" : "glyphicon glyphicon-flash"
    },
    "demo" : {
      "icon" : "glyphicon glyphicon-ok"
    }
  },
  "plugins" : ["types"]
});
```

[view result](http://jsfiddle.net/vakata/2kwkh2uL/4494/)

### unique
Enforces that no nodes with the same name can coexist as siblings - prevents renaming and moving nodes to a parent, which already contains a node with the same name.

```js
$("#tree").jstree({
  "plugins" : ["unique"]
});
```

[view result](http://jsfiddle.net/vakata/2kwkh2uL/4495/)

### wholerow
Makes each node appear block level which makes selection easier. May cause slow down for large trees in old browsers.

```js
$("#tree").jstree({
  "plugins" : ["wholerow"]
});
```

[view result](http://jsfiddle.net/vakata/2kwkh2uL/4496/)

### More plugins
If you create your own plugin (or download a 3rd party one) you must include its source on the page and list its name in the `"plugins"` config array.

```js
// conditional select
(function ($, undefined) {
  "use strict";
  $.jstree.defaults.conditionalselect = function () { return true; };
  $.jstree.plugins.conditionalselect = function (options, parent) {
    this.activate_node = function (obj, e) {
      if(this.settings.conditionalselect.call(this, this.get_node(obj))) {
        parent.activate_node.call(this, obj, e);
      }
    };
  };
})(jQuery);
$("#tree").jstree({
  "conditionalselect" : function (node) {
    return node.text === "Root node" ? false : true;
  },
  "plugins" : ["conditionalselect"]
});
```

[view result](http://jsfiddle.net/vakata/2kwkh2uL/4497/)

As seen here when creating a plugin you can define a default config, add your own functions to jstree, or override existing ones while maintaining the ability to call the overridden function.

## PHP demos moved to new repository
https://github.com/vakata/jstree-php-demos

## License & Contributing

_Please do NOT edit files in the "dist" subdirectory as they are generated via grunt. You'll find source code in the "src" subdirectory!_

If you want to you can always [donate a small amount][paypal] to help the development of jstree.

[paypal]: https://www.paypal.com/cgi-bin/webscr?cmd=_xclick&business=paypal@vakata.com&currency_code=USD&amount=&return=http://jstree.com/donation&item_name=Buy+me+a+coffee+for+jsTree

Copyright (c) 2014 Ivan Bozhanov (http://vakata.com)

Licensed under the [MIT license](http://www.opensource.org/licenses/mit-license.php).
# Spectrum
## The No Hassle Colorpicker

See the demo and docs: http://bgrins.github.io/spectrum.

I wanted a colorpicker that didn't require images, and that had an API that made sense to me as a developer who has worked with color in a number of applications.  I had tried a number of existing plugins, but decided to try and make a smaller, simpler one.

I started using canvas, then switched to CSS gradients, since it turned out to be easier to manage, and provided better cross browser support.

### Basic Usage

Head over to the [docs](http://bgrins.github.io/spectrum) for more information.  There is a visual demo of the different options hosted at: http://bgrins.github.io/spectrum.

    <script src='spectrum.js'></script>
    <link rel='stylesheet' href='spectrum.css' />

    <input id='colorpicker' />

    <script>
    $("#colorpicker").spectrum({
        color: "#f00"
    });
    </script>

### npm

Spectrum is registered as package with npm.  It can be installed with:

    npm install spectrum-colorpicker

### Bower

Spectrum is registered as a package with [Bower](http://bower.io/), so it can be pulled down using:

    bower install spectrum

### Using spectrum with a CDN

CDN provided by [cdnjs](https://cdnjs.com/libraries/spectrum)

    <script src="https://cdnjs.cloudflare.com/ajax/libs/spectrum/1.8.0/spectrum.min.js"></script>
    <link rel="stylesheet" type="text/css" href="https://cdnjs.cloudflare.com/ajax/libs/spectrum/1.8.0/spectrum.min.css">

### Continuous Integration

[![Build Status](https://secure.travis-ci.org/bgrins/spectrum.png?branch=master)](http://travis-ci.org/bgrins/spectrum)

Visit https://travis-ci.org/bgrins/spectrum to view the status of the automated tests.

### Building Spectrum Locally

If you'd like to download and use the plugin, head over to http://bgrins.github.io/spectrum/ and click the 'Download Zip' button.

If you'd like to run the development version, spectrum uses Grunt to automate the testing, linting, and building.  Head over to http://gruntjs.com/getting-started for more information.  First, clone the repository, then run:

    npm install -g grunt-cli
    npm install

    # runs jshint and the unit test suite
    grunt

    # runs jshint, the unit test suite, and builds a minified version of the file.
    grunt build

### Internationalization

If you are able to translate the text in the UI to another language, please do!  You can do so by either [filing a pull request](https://github.com/bgrins/spectrum/pulls) or [opening an issue]( https://github.com/bgrins/spectrum/issues) with the translation.  The existing languages are listed at: https://github.com/bgrins/spectrum/tree/master/i18n.

For an example, see the [Dutch translation](i18n/jquery.spectrum-nl.js).
# PROJ4JS [![Build Status](https://travis-ci.org/proj4js/proj4js.svg)](https://travis-ci.org/proj4js/proj4js)

Proj4js is a JavaScript library to transform point coordinates from one coordinate system to another, including datum transformations.
Originally a port of [PROJ](https://proj.org/) ([then known as PROJ.4](https://proj.org/faq.html#what-happened-to-proj-4)) and GCTCP C ([Archive](https://web.archive.org/web/20130523091752/http://edcftp.cr.usgs.gov/pub/software/gctpc/)) it is
a part of the [MetaCRS](https://trac.osgeo.org/metacrs/wiki) group of projects.

## Installing

Depending on your preferences

```bash
npm install proj4
bower install proj4
component install proj4js/proj4js
```

or just manually grab the file `proj4.js` from the [latest release](https://github.com/proj4js/proj4js/releases)'s `dist/` folder.

If you do not want to download anything, Proj4js is also hosted on [cdnjs](https://www.cdnjs.com/libraries/proj4js) for direct use in your browser applications.

## Using

The basic signature is:

```javascript
proj4(fromProjection[, toProjection, coordinates])
```

Projections can be proj or wkt strings.

Coordinates may an object of the form `{x:x,y:y}` or an array of the form `[x,y]`.

When all 3 arguments  are given, the result is that the coordinates are transformed from projection1 to projection 2. And returned in the same format that they were given in.

```javascript
var firstProjection = 'PROJCS["NAD83 / Massachusetts Mainland",GEOGCS["NAD83",DATUM["North_American_Datum_1983",SPHEROID["GRS 1980",6378137,298.257222101,AUTHORITY["EPSG","7019"]],AUTHORITY["EPSG","6269"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4269"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Lambert_Conformal_Conic_2SP"],PARAMETER["standard_parallel_1",42.68333333333333],PARAMETER["standard_parallel_2",41.71666666666667],PARAMETER["latitude_of_origin",41],PARAMETER["central_meridian",-71.5],PARAMETER["false_easting",200000],PARAMETER["false_northing",750000],AUTHORITY["EPSG","26986"],AXIS["X",EAST],AXIS["Y",NORTH]]';
var secondProjection = "+proj=gnom +lat_0=90 +lon_0=0 +x_0=6300000 +y_0=6300000 +ellps=WGS84 +datum=WGS84 +units=m +no_defs";
//I'm not going to redefine those two in latter examples.
proj4(firstProjection,secondProjection,[2,5]);
// [-2690666.2977344505, 3662659.885459918]
```

If only 1 projection is given then it is assumed that it is being projected *from* WGS84 (fromProjection is WGS84).

```javascript
proj4(firstProjection,[-71,41]);
// [242075.00535055372, 750123.32090043]
```

If no coordinates are given an object with two methods is returned, its methods are `forward` which projects from the first projection to the second and `inverse` which projects from the second to the first.

```javascript
proj4(firstProjection,secondProjection).forward([2,5]);
// [-2690666.2977344505, 3662659.885459918]
proj4(secondProjection,firstProjection).inverse([2,5]);
// [-2690666.2977344505, 3662659.885459918]
```

And as above if only one projection is given, it's assumed to be coming from wgs84:

```javascript
proj4(firstProjection).forward([-71,41]);
// [242075.00535055372, 750123.32090043]
proj4(firstProjection).inverse([242075.00535055372, 750123.32090043]);
//[-71, 40.99999999999986]
//the floating points to answer your question
```

## Named Projections

If you prefer to define a projection as a string and reference it that way, you may use the proj4.defs method which can be called 2 ways, with a name and projection:

```js
proj4.defs('WGS84', "+title=WGS 84 (long/lat) +proj=longlat +ellps=WGS84 +datum=WGS84 +units=degrees");
```

or with an array

```js
proj4.defs([
  [
    'EPSG:4326',
    '+title=WGS 84 (long/lat) +proj=longlat +ellps=WGS84 +datum=WGS84 +units=degrees'],
  [
    'EPSG:4269',
    '+title=NAD83 (long/lat) +proj=longlat +a=6378137.0 +b=6356752.31414036 +ellps=GRS80 +datum=NAD83 +units=degrees'
  ]
]);
```

you can then do

```js
proj4('EPSG:4326');
```

instead of writing out the whole proj definition, by default proj4 has the following projections predefined:

- 'EPSG:4326', which has the following alias
    - 'WGS84'
- 'EPSG:4269'
- 'EPSG:3857', which has the following aliases
    - 'EPSG:3785'
    - 'GOOGLE'
    - 'EPSG:900913'
    - 'EPSG:102113'

Defined projections can also be accessed through the proj4.defs function (`proj4.defs('EPSG:4326')`).

proj4.defs can also be used to define a named alias:

```javascript
proj4.defs('urn:x-ogc:def:crs:EPSG:4326', proj4.defs('EPSG:4326'));
```

## TypeScript

TypeScript implementation was added to the [DefinitelyTyped repository](https://github.com/DefinitelyTyped/DefinitelyTyped).

```bash
$ npm install --save @types/proj4
```

## Developing
To set up build tools make sure you have node and grunt-cli installed and then run `npm install`.

To do the complete build and browser tests run:

```bash
node_modules/.bin/grunt
```

To run node tests run:

```bash
npm test
```

To run node tests with coverage run:

```bash
npm test --coverage
```

To create a build with only default projections (latlon and Mercator) run:

```bash
node_modules/.bin/grunt build
```

To create a build with only custom projections include a comma separated list of projections codes (the file name in 'lib/projections' without the '.js') after a colon, e.g.:

```bash
node_modules/.bin/grunt build:tmerc
#includes transverse Mercator
node_modules/.bin/grunt build:lcc
#includes lambert conformal conic
node_modules/.bin/grunt build:omerc,moll
#includes oblique Mercator and Mollweide
```
# MeshLine
Mesh replacement for ```THREE.Line```

Instead of using GL_LINE, it uses a strip of triangles billboarded. Some examples:

[![Demo](screenshots/demo.jpg)](https://www.clicktorelease.com/code/THREE.MeshLine/demo/index.html)
[![Graph](screenshots/graph.jpg)](https://www.clicktorelease.com/code/THREE.MeshLine/demo/graph.html)
[![Spinner](screenshots/spinner.jpg)](https://www.clicktorelease.com/code/THREE.MeshLine/demo/spinner.html)
[![SVG](screenshots/svg.jpg)](https://www.clicktorelease.com/code/THREE.MeshLine/demo/svg.html)
[![Shape](screenshots/shape.jpg)](https://www.clicktorelease.com/code/THREE.MeshLine/demo/shape.html)
[![Shape](screenshots/birds.jpg)](https://www.clicktorelease.com/code/THREE.MeshLine/demo/birds.html)

* [Demo](https://www.clicktorelease.com/code/THREE.MeshLine/demo/index.html): play with the different settings of materials
* [Graph](https://www.clicktorelease.com/code/THREE.MeshLine/demo/graph.html): example of using ```MeshLine``` to plot graphs
* [Spinner](https://www.clicktorelease.com/code/THREE.MeshLine/demo/spinner.html): example of dynamic ```MeshLine``` with texture
* [SVG](https://www.clicktorelease.com/code/THREE.MeshLine/demo/svg.html): example of ```MeshLine``` rendering SVG Paths
* [Shape](https://www.clicktorelease.com/code/THREE.MeshLine/demo/shape.html): example of ```MeshLine``` created from a mesh
* [Birds](https://www.clicktorelease.com/code/THREE.MeshLine/demo/birds.html): example of ```MeshLine.advance()``` by @caramelcode (Jared Sprague) and @mwcz (Michael Clayton)

### How to use ####

* Include script
* Create and populate a geometry
* Create a MeshLine and assign the geometry
* Create a MeshLineMaterial
* Use MeshLine and MeshLineMaterial to create a THREE.Mesh

#### Include the script

Include script after THREE is included
```js
<script src="THREE.MeshLine.js"></script>
```
or use npm to install it
```
npm i three.meshline
```
and include it in your code (don't forget to require three.js)
```js
var THREE = require( 'three' );
var MeshLine = require( 'three.meshline' );
```

##### Create and populate a geometry #####

First, create the list of vertices that will define the line. ```MeshLine``` accepts ```THREE.Geometry``` (looking up the ```.vertices``` in it) and ```Array```/```Float32Array```. ```THREE.BufferGeometry``` coming soon, and may be others like ```Array``` of ```THREE.Vector3```.

```js
var geometry = new THREE.Geometry();
for( var j = 0; j < Math.PI; j += 2 * Math.PI / 100 ) {
	var v = new THREE.Vector3( Math.cos( j ), Math.sin( j ), 0 );
	geometry.vertices.push( v );
}
```

##### Create a MeshLine and assign the geometry #####

Once you have that, you can create a new ```MeshLine```, and call ```.setGeometry()``` passing the vertices.

```js
var line = new MeshLine();
line.setGeometry( geometry );
```

Note: ```.setGeometry``` accepts a second parameter, which is a function to define the width in each point along the line. By default that value is 1, making the line width 1 * lineWidth.

```js
line.setGeometry( geometry, function( p ) { return 2; } ); // makes width 2 * lineWidth
line.setGeometry( geometry, function( p ) { return 1 - p; } ); // makes width taper
line.setGeometry( geometry, function( p ) { return 2 + Math.sin( 50 * p ); } ); // makes width sinusoidal
```

##### Create a MeshLineMaterial #####

A ```MeshLine``` needs a ```MeshLineMaterial```:

```js
var material = new MeshLineMaterial();
```

By default it's a white material of width 1 unit.

```MeshLineMaterial``` has several attributes to control the appereance of the ```MeshLine```:

* ```map``` - a ```THREE.Texture``` to paint along the line (requires ```useMap``` set to true)
* ```useMap``` - tells the material to use ```map``` (0 - solid color, 1 use texture)
* ```alphaMap``` - a ```THREE.Texture``` to use as alpha along the line (requires ```useAlphaMap``` set to true)
* ```useAlphaMap``` - tells the material to use ```alphaMap``` (0 - no alpha, 1 modulate alpha)
* ```repeat``` - THREE.Vector2 to define the texture tiling (applies to map and alphaMap - MIGHT CHANGE IN THE FUTURE)
* ```color``` - ```THREE.Color``` to paint the line width, or tint the texture with
* ```opacity``` - alpha value from 0 to 1 (requires ```transparent``` set to ```true```)
* ```alphaTest``` - cutoff value from 0 to 1
* ```dashArray``` - THREE.Vector2 to define the dashing (NOT IMPLEMENTED YET)
* ```resolution``` - ```THREE.Vector2``` specifying the canvas size (REQUIRED)
* ```sizeAttenuation``` - makes the line width constant regardless distance (1 unit is 1px on screen) (0 - attenuate, 1 - don't attenuate)
* ```lineWidth``` - float defining width (if ```sizeAttenuation``` is true, it's world units; else is screen pixels)
* ```near``` - camera near clip plane distance  (REQUIRED if ```sizeAttenuation``` set to false)
* ```far``` - camera far clip plane distance  (REQUIRED if ```sizeAttenuation``` set to false)

If you're rendering transparent lines or using a texture with alpha map, you should set ```depthTest``` to ```false```, ```transparent``` to ```true``` and ```blending``` to an appropriate blending mode, or use ```alphaTest```.

##### Use MeshLine and MeshLineMaterial to create a THREE.Mesh #####

Finally, we create a mesh and add it to the scene:

```js
var mesh = new THREE.Mesh( line.geometry, material ); // this syntax could definitely be improved!
scene.add( mesh );
```

### TODO ###

* Better miters
* Proper sizes
* Support for dashArray

### Support ###

Tested successfully on

* Chrome OSX, Windows, Android
* Firefox OSX, Windows, Anroid
* Safari OSX, iOS
* Internet Explorer 11 (SVG and Shape demo won't work because they use Promises)
* Opera OSX, Windows

### References ###

* [Drawing lines is hard](http://mattdesl.svbtle.com/drawing-lines-is-hard)
* [WebGL rendering of solid trails](http://codeflow.org/entries/2012/aug/05/webgl-rendering-of-solid-trails/)
* [Drawing Antialiased Lines with OpenGL](https://www.mapbox.com/blog/drawing-antialiased-lines/)

#### License ####

MIT licensed

Copyright (C) 2015-2016 Jaume Sanchez Elias, http://www.clicktorelease.com
# GeoPackage JS

GeoPackage JS is an implementation of the OGC GeoPackage spec.  This library works in both the browser and Node 4+.

### Demo ###
[GeoPackage JS Demo Page](http://ngageoint.github.io/geopackage-js/)

Cloning this repository and opening the docs/index.html in your browser will run the demo locally.

### Installation ###

[![Build Status](https://travis-ci.org/ngageoint/geopackage-js.svg?branch=master)](https://travis-ci.org/ngageoint/geopackage-js)
[![NPM](https://img.shields.io/npm/v/@ngageoint/geopackage.svg)](https://www.npmjs.com/package/@ngageoint/geopackage)
[![Coverage Status](https://coveralls.io/repos/github/ngageoint/geopackage-js/badge.svg)](https://coveralls.io/github/ngageoint/geopackage-js)

```sh
$ npm install @ngageoint/geopackage
```

#### GeoPackage JS Library ####

The [GeoPackage Libraries](http://ngageoint.github.io/GeoPackage/) were developed at the [National Geospatial-Intelligence Agency (NGA)](http://www.nga.mil/) in collaboration with [BIT Systems](http://www.bit-sys.com/). The government has "unlimited rights" and is releasing this software to increase the impact of government investments by providing developers with the opportunity to take things in new directions. The software use, modification, and distribution rights are stipulated within the [MIT license](http://choosealicense.com/licenses/mit/).

### Pull Requests ###
If you'd like to contribute to this project, please make a pull request. We'll review the pull request and discuss the changes. All pull request contributions to this project will be released under the MIT license.

Software source code previously released under an open source license and then modified by NGA staff is considered a "joint work" (see 17 USC § 101); it is partially copyrighted, partially public domain, and as a whole is protected by the copyrights of the non-government authors and must be released according to the terms of the original open source license.

### About ###

[GeoPackage JS](https://github.com/ngageoint/geopackage-js) is a [GeoPackage Library](http://ngageoint.github.io/GeoPackage/) JavaScript implementation of the Open Geospatial Consortium [GeoPackage](http://www.geopackage.org/) [spec](http://www.geopackage.org/spec/).  It is listed as an [OGC GeoPackage Implementation](http://www.geopackage.org/#implementations_nga) by the National Geospatial-Intelligence Agency.

The GeoPackage JavaScript library currently provides the ability to read GeoPackage files.  This library works both in the browser and in Node.  In the browser tiles are rendered using HTML5 Canvas and GeoPackages are read using [sql.js](https://github.com/kripken/sql.js/).  In Node tiles are rendered  [PureImage](https://github.com/joshmarinacci/node-pureimage) and GeoPackages are read using [node-sqlite3](https://github.com/mapbox/node-sqlite3).

### Changelog

##### 2.1.0

- Implementation of the Feature Style Extension and Contents ID Extension

##### 2.0.8

- Checks for Electron when returning a tile creator

##### 2.0

- All new API utilizing Promises

##### 1.1.4

- Adds a method to retrieve tiles in EPSG:4326

##### 1.1.3

- Fixes issue #115

##### 1.1.2

- fix case where GeoPackage Zoom does not correspond to the web map zoom

##### 1.1.1

- fix more instances of proj4 bug for react
- fixed tile generation for images with different x and y pixel densities

##### 1.1.0

- accept pull request adding support for react
- fix bug with projected tiles that spanned the date line

##### 1.0.25

- ensure we use proj4 2.4.3 instead of 2.4.4

##### 1.0.22

- Fixed bug where querying for indexed features only returned the geometry instead of the entire feature

##### 1.0.19

- Remove dependency on Lwip

### Usage ###

View examples using [Bower](https://github.com/ngageoint/geopackage-js/tree/master/docs/bower) and [Browserify](https://github.com/ngageoint/geopackage-js/tree/master/docs)

View the latest [docs](http://ngageoint.github.io/geopackage-js/jsdoc/module-geoPackage-GeoPackage.html) (currently being updated).

#### Browser Usage ####
```javascript

// attach this method to a file input onchange event
window.loadGeoPackage = function(files) {
  var f = files[0];
  var r = new FileReader();
  r.onload = function() {
    var array = new Uint8Array(r.result);
    loadByteArray(array);
  }
  r.readAsArrayBuffer(f);
}

function loadByteArray(array, callback) {
  var db = new SQL.Database(array);
  GeoPackageConnection.connectWithDatabase(db, function(err, connection) {
    var geoPackage = new GeoPackage('', '', connection);

    // Now you can operate on the GeoPackage

    // get the tile table names
    geoPackage.getTileTables(function(err, tileTableNames) {
      // tileTableNames is an array of all tile table names

      // get the info for the first table
      geoPackage.getTileDaoWithTableName(tileTableNames[0], function(err, tileDao) {
        geoPackage.getInfoForTable(tileDao, function(err, info) {
          // do something with the tile table info
        });

        // draw a tile into a canvas for an XYZ tile
        var canvas = canvasFromSomewhere;
        var gpr = new GeoPackageTileRetriever(tileDao, 256, 256);
        var x = 0;
        var y = 0;
        var zoom = 0;

        console.time('Draw tile ' + x + ', ' + y + ' zoom: ' + zoom);
        gpr.drawTileIn(x, y, zoom, canvas, function() {
          console.timeEnd('Draw tile ' + x + ', ' + y + ' zoom: ' + zoom);
        });

        // or get a tile base64 data URL for an XYZ tile
        gpr.getTile(x, y, zoom, function(err, tileBase64DataURL) {
          console.log('got the base64 data url');
        });

        // or get a tile from a GeoPackage tile column and tile row
        tileDao.queryForTile(tileColumn, tileRow, zoom, function(err, tile) {
          var tileData = tile.getTileData();  // the raw bytes from the GeoPackage
        });

      });
    });

    // get the feature table names
    geoPackage.getFeatureTables(function(err, featureTableNames) {
      // featureTableNames is an array of all feature table names

      // get the info for the first table
      geoPackage.getFeatureDaoWithTableName(featureTableNames[0], function(err, featureDao) {
        geoPackage.getInfoForTable(featureDao, function(err, info) {
          // do something with the feature table info
        });

        // query for all features
        featureDao.queryForEach(function(err, row, rowDone) {
          var feature = featureDao.getFeatureRow(row);
          var geometry = currentRow.getGeometry();
          if (geometry) {
            var geom = geometry.geometry;
            var geoJson = geometry.geometry.toGeoJSON();

            geoJson.properties = {};
            for (var key in feature.values) {
              if(feature.values.hasOwnProperty(key) && key != feature.getGeometryColumn().name) {
                var column = info.columnMap[key];
                geoJson.properties[column.displayName] = currentRow.values[key];
              }
            }
          }
          rowDone();
        });
      });
    });
  });
}

```

#### NodeJS Usage ####

```javascript
var GeoPackageAPI = require('@ngageoint/geopackage')
  , GeoPackageManager = GeoPackageAPI.GeoPackageManager
  , GeoPackageConnection = GeoPackageAPI.GeoPackageConnection
  , GeoPackageTileRetriever = GeoPackageAPI.GeoPackageTileRetriever;

GeoPackageAPI.open(filename, function(err, geoPackage) {

  // Now you can operate on the GeoPackage

  // get the tile table names
  geoPackage.getTileTables(function(err, tileTableNames) {
    // tileTableNames is an array of all tile table names

    // get the info for the first table
    geoPackage.getTileDaoWithTableName(tileTableNames[0], function(err, tileDao) {
      geoPackage.getInfoForTable(tileDao, function(err, info) {
        // do something with the tile table info
      });

      // draw a tile into a canvas for an XYZ tile
      var canvas = canvasFromSomewhere;
      var gpr = new GeoPackageTileRetriever(tileDao, 256, 256);
      var x = 0;
      var y = 0;
      var zoom = 0;

      console.time('Draw tile ' + x + ', ' + y + ' zoom: ' + zoom);
      gpr.drawTileIn(x, y, zoom, canvas, function() {
        console.timeEnd('Draw tile ' + x + ', ' + y + ' zoom: ' + zoom);
      });

      // or get a tile base64 data URL for an XYZ tile
      gpr.getTile(x, y, zoom, function(err, tileBase64DataURL) {
        console.log('got the base64 data url');
      });

      // or get a tile from a GeoPackage tile column and tile row
      tileDao.queryForTile(tileColumn, tileRow, zoom, function(err, tile) {
        var tileData = tile.getTileData();  // the raw bytes from the GeoPackage
      });

    });
  });

  // get the feature table names
  geoPackage.getFeatureTables(function(err, featureTableNames) {
    // featureTableNames is an array of all feature table names

    // get the info for the first table
    geoPackage.getFeatureDaoWithTableName(featureTableNames[0], function(err, featureDao) {
      geoPackage.getInfoForTable(featureDao, function(err, info) {
        // do something with the feature table info
      });

      // query for all features
      featureDao.queryForEach(function(err, row, rowDone) {
        var feature = featureDao.getFeatureRow(row);
        var geometry = currentRow.getGeometry();
        if (geometry) {
          var geom = geometry.geometry;
          var geoJson = geometry.geometry.toGeoJSON();

          geoJson.properties = {};
          for (var key in feature.values) {
            if(feature.values.hasOwnProperty(key) && key != feature.getGeometryColumn().name) {
              var column = info.columnMap[key];
              geoJson.properties[column.displayName] = currentRow.values[key];
            }
          }
        }
        rowDone();
      });
    });
  });
});

```
# SQLite compiled to javascript
[![Build Status](https://travis-ci.org/kripken/sql.js.svg?branch=master)](http://travis-ci.org/kripken/sql.js) [![CDNJS version](https://img.shields.io/cdnjs/v/sql.js.svg)](https://cdnjs.com/libraries/sql.js)

For the impatients, try the demo here: http://kripken.github.io/sql.js/examples/GUI

*sql.js* is a port of [SQLite](http://sqlite.org/about.html) to Webassembly, by compiling the SQLite C code with [Emscripten](http://kripken.github.io/emscripten-site/docs/introducing_emscripten/about_emscripten.html). It uses a [virtual database file stored in memory](https://kripken.github.io/emscripten-site/docs/porting/files/file_systems_overview.html), and thus **doesn't persist the changes** made to the database. However, it allows you to **import** any existing sqlite file, and to **export** the created database as a [javascript typed array](https://developer.mozilla.org/en-US/docs/Web/JavaScript/Typed_arrays).

There are no C bindings or node-gyp compilation here, sql.js is a simple javascript file, that can be used like any traditional javascript library. If you are building a native application in javascript (using Electron for instance), or are working in node.js, you will likely prefer to use [a native binding of SQLite to javascript](https://www.npmjs.com/package/sqlite3).

SQLite is public domain, sql.js is MIT licensed.

Sql.js predates WebAssembly, and thus started as an [asm.js](https://en.wikipedia.org/wiki/Asm.js) project. It still supports asm.js for backwards compatability.

## Version of binaries
Sql.js was last built with:
Emscripten version 1.38.30 (2019-04-16) [Release History](https://emscripten.org/docs/introducing_emscripten/release_notes.html)
SqlLite version: 3.28.0 (2019-04-16) [Release History](https://www.sqlite.org/changes.html)

## Documentation
A [full documentation](http://kripken.github.io/sql.js/documentation/#http://kripken.github.io/sql.js/documentation/class/Database.html) generated from comments inside the source code, is available.

## Usage

```javascript
var initSqlJs = require('sql-wasm.js');
// or if you are in a browser:
//var initSqlJs = window.initSqlJs;

initSqlJs().then(function(SQL){

  // Create a database
  var db = new SQL.Database();
  // NOTE: You can also use new SQL.Database(data) where
  // data is an Uint8Array representing an SQLite database file

  // Execute some sql
  sqlstr = "CREATE TABLE hello (a int, b char);";
  sqlstr += "INSERT INTO hello VALUES (0, 'hello');"
  sqlstr += "INSERT INTO hello VALUES (1, 'world');"
  db.run(sqlstr); // Run the query without returning anything

  var res = db.exec("SELECT * FROM hello");
  /*
  [
    {columns:['a','b'], values:[[0,'hello'],[1,'world']]}
  ]
  */

  // Prepare an sql statement
  var stmt = db.prepare("SELECT * FROM hello WHERE a=:aval AND b=:bval");

  // Bind values to the parameters and fetch the results of the query
  var result = stmt.getAsObject({':aval' : 1, ':bval' : 'world'});
  console.log(result); // Will print {a:1, b:'world'}

  // Bind other values
  stmt.bind([0, 'hello']);
  while (stmt.step()) console.log(stmt.get()); // Will print [0, 'hello']

  // You can also use javascript functions inside your SQL code
  // Create the js function you need
  function add(a, b) {return a+b;}
  // Specifies the SQL function's name, the number of it's arguments, and the js function to use
  db.create_function("add_js", add);
  // Run a query in which the function is used
  db.run("INSERT INTO hello VALUES (add_js(7, 3), add_js('Hello ', 'world'));"); // Inserts 10 and 'Hello world'

  // free the memory used by the statement
  stmt.free();
  // You can not use your statement anymore once it has been freed.
  // But not freeing your statements causes memory leaks. You don't want that.

  // Export the database to an Uint8Array containing the SQLite database file
  var binaryArray = db.export();
});

```

## Demo
There are a few examples [available here](https://kripken.github.io/sql.js/index.html). The most full-featured is the [Sqlite Interpreter](https://kripken.github.io/sql.js/examples/GUI/index.html).

## Examples
The test files provide up to date example of the use of the api.
### Inside the browser
#### Example **HTML** file:
```html
<meta charset="utf8" />
<html>
  <script src='/dist/sql-wasm.js'></script>
  <script>
    config = {
      locateFile: url => `/dist/${filename}` 
    }
    // The `initSqlJs` function is globally provided by all of the main dist files if loaded in the browser.
    // We must specify this locateFile function if we are loading a wasm file from anywhere other than the current html page's folder.
    initSqlJs(config).then(function(SQL){
      //Create the database
      var db = new SQL.Database();
      // Run a query without reading the results
      db.run("CREATE TABLE test (col1, col2);");
      // Insert two rows: (1,111) and (2,222)
      db.run("INSERT INTO test VALUES (?,?), (?,?)", [1,111,2,222]);
  
      // Prepare a statement
      var stmt = db.prepare("SELECT * FROM test WHERE col1 BETWEEN $start AND $end");
      stmt.getAsObject({$start:1, $end:1}); // {col1:1, col2:111}
  
      // Bind new values
      stmt.bind({$start:1, $end:2});
      while(stmt.step()) { //
        var row = stmt.getAsObject();
        console.log('Here is a row: ' + JSON.stringify(row));
      }
    });
  </script>
  <body>
    Output is in Javscript console
  </body>
</html>
```

#### Creating a database from a file choosen by the user
`SQL.Database` constructor takes an array of integer representing a database file as an optional parameter.
The following code uses an HTML input as the source for loading a database:
```javascript
dbFileElm.onchange = () => {
  var f = dbFileElm.files[0];
  var r = new FileReader();
  r.onload = function() {
    var Uints = new Uint8Array(r.result);
    db = new SQL.Database(Uints);
  }
  r.readAsArrayBuffer(f);
}
```
See : http://kripken.github.io/sql.js/examples/GUI/gui.js

#### Loading a database from a server

```javascript
var xhr = new XMLHttpRequest();
// For example: https://github.com/lerocha/chinook-database/raw/master/ChinookDatabase/DataSources/Chinook_Sqlite.sqlite
xhr.open('GET', '/path/to/database.sqlite', true);
xhr.responseType = 'arraybuffer';

xhr.onload = e => {
  var uInt8Array = new Uint8Array(this.response);
  var db = new SQL.Database(uInt8Array);
  var contents = db.exec("SELECT * FROM my_table");
  // contents is now [{columns:['col1','col2',...], values:[[first row], [second row], ...]}]
};
xhr.send();
```
See: https://github.com/kripken/sql.js/wiki/Load-a-database-from-the-server


### Use from node.js

`sql.js` is [hosted on npm](https://www.npmjs.org/package/sql.js). To install it, you can simply run `npm install sql.js`.
Alternatively, you can simply download `sql-wasm.js` and `sql-wasm.wasm`, from the download link below.

#### read a database from the disk:
```javascript
var fs = require('fs');
var initSqlJs = require('sql-wasm.js');
var filebuffer = fs.readFileSync('test.sqlite');
 
initSqlJs().then(function(SQL){
  // Load the db
  var db = new SQL.Database(filebuffer);
});

```

#### write a database to the disk
You need to convert the result of `db.export` to a buffer
```javascript
var fs = require("fs");
// [...] (create the database)
var data = db.export();
var buffer = new Buffer(data);
fs.writeFileSync("filename.sqlite", buffer);
```

See : https://github.com/kripken/sql.js/blob/master/test/test_node_file.js

### Use as web worker
If you don't want to run CPU-intensive SQL queries in your main application thread,
you can use the *more limited* WebWorker API.

You will need to download [dist/worker.sql-wasm.js](dist/worker.sql-wasm.js) [dist/worker.sql-wasm.wasm](dist/worker.sql-wasm.wasm).

Example:
```html
<script>
  var worker = new Worker("/dist/worker.sql-wasm.js");
  worker.onmessage = () => {
    console.log("Database opened");
    worker.onmessage = event => {
      console.log(event.data); // The result of the query
    };
	
    worker.postMessage({
      id: 2,
      action: 'exec',
      sql: 'SELECT * FROM test'
    });
  };

  worker.onerror = e => console.log("Worker error: ", e);
  worker.postMessage({
    id:1,
    action:'open',
    buffer:buf, /*Optional. An ArrayBuffer representing an SQLite Database file*/
  });
</script>
```

See [examples/GUI/gui.js](examples/GUI/gui.js) for a full working example.

## Flavors/versions Targets/Downloads

This library includes both WebAssembly and asm.js versions of Sqlite. (WebAssembly is the newer, preferred way to compile to Javascript, and has superceded asm.js. It produces smaller, faster code.) Asm.js versions are included for compatibility.

## Upgrading from 0.x to 1.x

Version 1.0 of sql.js must be loaded asynchronously, whereas asm.js was able to be loaded synchronously. 

So in the past, you would:
```html
<script src='js/sql.js'></script>
<script>
  var db = new SQL.Database();
  //...
</script>
```
or:
```javascript
var SQL = require('sql.js');
var db = new QL.Database();
//...
```

Version 1.x:
```html
<script src='dist/sql-wasm.js'></script>
<script>
  initSqlJs({ locateFile: filename => `/dist/${filename}` }).then(function(SQL){
    var db = new SQL.Database();
    //...
  });
</script>
```
or:
```javascript
var initSqlJs = require('sql-wasm.js');
initSqlJs().then(function(SQL){
  var db = new SQL.Database();
  //...
});
```

`NOTHING` is now a reserved word in SQLite, whereas previously it was not. This could cause errors like `Error: near "nothing": syntax error`

### Downloading/Using: ###
Although asm.js files were distributed as a single Javascript file, WebAssembly libraries are most efficiently distributed as a pair of files, the `.js`  loader and the `.wasm` file, like [dist/sql-wasm.js]([dist/sql-wasm.js]) and [dist/sql-wasm.wasm]([dist/sql-wasm.wasm]). The `.js` file is reponsible for wrapping/loading the `.wasm` file. 




## Versions of sql.js included in `dist/`
 - `sql-wasm.js` : The Web Assembly version of Sql.js. Minified and suitable for production. Use this. If you use this, you will need to include/ship `sql-wasm.wasm` as well.
 - `sql-wasm-debug.js` : The Web Assembly, Debug version of Sql.js. Larger, with assertions turned on. Useful for local development. You will need to include/ship `sql-wasm-debug.wasm` if you use this.
 - `sql-asm.js` : The older asm.js version of Sql.js. Slower and larger. Provided for compatiblity reasons.
 - `sql-asm-memory-growth.js` : Asm.js doesn't allow for memory to grow by default, because it is slower and de-optimizes. If you are using sql-asm.js and you see this error (`Cannot enlarge memory arrays`), use this file.
 - `sql-asm-debug.js` : The _Debug_ asm.js version of Sql.js. Use this for local development.
 - `worker.*` - Web Worker versions of the above libraries. More limited API. See [examples/GUI/gui.js](examples/GUI/gui.js) for a good example of this.

## Compiling

- Install the EMSDK, [as described here](https://kripken.github.io/emscripten-site/docs/getting_started/downloads.html)
- Run `npm run rebuild`

# GeoPackage JS

GeoPackage JS is an implementation of the OGC GeoPackage spec.  This library works in both the browser and Node 4+.

### Demo ###
[GeoPackage JS Demo Page](http://ngageoint.github.io/geopackage-js/)

Cloning this repository and opening the docs/index.html in your browser will run the demo locally.

### Installation ###

[![Build Status](https://travis-ci.org/ngageoint/geopackage-js.svg?branch=master)](https://travis-ci.org/ngageoint/geopackage-js)
[![NPM](https://img.shields.io/npm/v/@ngageoint/geopackage.svg)](https://www.npmjs.com/package/@ngageoint/geopackage)
[![Coverage Status](https://coveralls.io/repos/github/ngageoint/geopackage-js/badge.svg)](https://coveralls.io/github/ngageoint/geopackage-js)

```sh
$ npm install @ngageoint/geopackage
```

#### GeoPackage JS Library ####

The [GeoPackage Libraries](http://ngageoint.github.io/GeoPackage/) were developed at the [National Geospatial-Intelligence Agency (NGA)](http://www.nga.mil/) in collaboration with [BIT Systems](http://www.bit-sys.com/). The government has "unlimited rights" and is releasing this software to increase the impact of government investments by providing developers with the opportunity to take things in new directions. The software use, modification, and distribution rights are stipulated within the [MIT license](http://choosealicense.com/licenses/mit/).

### Pull Requests ###
If you'd like to contribute to this project, please make a pull request. We'll review the pull request and discuss the changes. All pull request contributions to this project will be released under the MIT license.

Software source code previously released under an open source license and then modified by NGA staff is considered a "joint work" (see 17 USC § 101); it is partially copyrighted, partially public domain, and as a whole is protected by the copyrights of the non-government authors and must be released according to the terms of the original open source license.

### About ###

[GeoPackage JS](https://github.com/ngageoint/geopackage-js) is a [GeoPackage Library](http://ngageoint.github.io/GeoPackage/) JavaScript implementation of the Open Geospatial Consortium [GeoPackage](http://www.geopackage.org/) [spec](http://www.geopackage.org/spec/).  It is listed as an [OGC GeoPackage Implementation](http://www.geopackage.org/#implementations_nga) by the National Geospatial-Intelligence Agency.

The GeoPackage JavaScript library currently provides the ability to read GeoPackage files.  This library works both in the browser and in Node.  In the browser tiles are rendered using HTML5 Canvas and GeoPackages are read using [sql.js](https://github.com/kripken/sql.js/).  In Node tiles are rendered  [PureImage](https://github.com/joshmarinacci/node-pureimage) and GeoPackages are read using [node-sqlite3](https://github.com/mapbox/node-sqlite3).

### Changelog

##### 2.1.0

- Implementation of the Feature Style Extension and Contents ID Extension

##### 2.0.8

- Checks for Electron when returning a tile creator

##### 2.0

- All new API utilizing Promises

##### 1.1.4

- Adds a method to retrieve tiles in EPSG:4326

##### 1.1.3

- Fixes issue #115

##### 1.1.2

- fix case where GeoPackage Zoom does not correspond to the web map zoom

##### 1.1.1

- fix more instances of proj4 bug for react
- fixed tile generation for images with different x and y pixel densities

##### 1.1.0

- accept pull request adding support for react
- fix bug with projected tiles that spanned the date line

##### 1.0.25

- ensure we use proj4 2.4.3 instead of 2.4.4

##### 1.0.22

- Fixed bug where querying for indexed features only returned the geometry instead of the entire feature

##### 1.0.19

- Remove dependency on Lwip

### Usage ###

View examples using [Bower](https://github.com/ngageoint/geopackage-js/tree/master/docs/bower) and [Browserify](https://github.com/ngageoint/geopackage-js/tree/master/docs)

View the latest [docs](http://ngageoint.github.io/geopackage-js/jsdoc/module-geoPackage-GeoPackage.html) (currently being updated).

#### Browser Usage ####
```javascript

// attach this method to a file input onchange event
window.loadGeoPackage = function(files) {
  var f = files[0];
  var r = new FileReader();
  r.onload = function() {
    var array = new Uint8Array(r.result);
    loadByteArray(array);
  }
  r.readAsArrayBuffer(f);
}

function loadByteArray(array, callback) {
  var db = new SQL.Database(array);
  GeoPackageConnection.connectWithDatabase(db, function(err, connection) {
    var geoPackage = new GeoPackage('', '', connection);

    // Now you can operate on the GeoPackage

    // get the tile table names
    geoPackage.getTileTables(function(err, tileTableNames) {
      // tileTableNames is an array of all tile table names

      // get the info for the first table
      geoPackage.getTileDaoWithTableName(tileTableNames[0], function(err, tileDao) {
        geoPackage.getInfoForTable(tileDao, function(err, info) {
          // do something with the tile table info
        });

        // draw a tile into a canvas for an XYZ tile
        var canvas = canvasFromSomewhere;
        var gpr = new GeoPackageTileRetriever(tileDao, 256, 256);
        var x = 0;
        var y = 0;
        var zoom = 0;

        console.time('Draw tile ' + x + ', ' + y + ' zoom: ' + zoom);
        gpr.drawTileIn(x, y, zoom, canvas, function() {
          console.timeEnd('Draw tile ' + x + ', ' + y + ' zoom: ' + zoom);
        });

        // or get a tile base64 data URL for an XYZ tile
        gpr.getTile(x, y, zoom, function(err, tileBase64DataURL) {
          console.log('got the base64 data url');
        });

        // or get a tile from a GeoPackage tile column and tile row
        tileDao.queryForTile(tileColumn, tileRow, zoom, function(err, tile) {
          var tileData = tile.getTileData();  // the raw bytes from the GeoPackage
        });

      });
    });

    // get the feature table names
    geoPackage.getFeatureTables(function(err, featureTableNames) {
      // featureTableNames is an array of all feature table names

      // get the info for the first table
      geoPackage.getFeatureDaoWithTableName(featureTableNames[0], function(err, featureDao) {
        geoPackage.getInfoForTable(featureDao, function(err, info) {
          // do something with the feature table info
        });

        // query for all features
        featureDao.queryForEach(function(err, row, rowDone) {
          var feature = featureDao.getFeatureRow(row);
          var geometry = currentRow.getGeometry();
          if (geometry) {
            var geom = geometry.geometry;
            var geoJson = geometry.geometry.toGeoJSON();

            geoJson.properties = {};
            for (var key in feature.values) {
              if(feature.values.hasOwnProperty(key) && key != feature.getGeometryColumn().name) {
                var column = info.columnMap[key];
                geoJson.properties[column.displayName] = currentRow.values[key];
              }
            }
          }
          rowDone();
        });
      });
    });
  });
}

```

#### NodeJS Usage ####

```javascript
var GeoPackageAPI = require('@ngageoint/geopackage')
  , GeoPackageManager = GeoPackageAPI.GeoPackageManager
  , GeoPackageConnection = GeoPackageAPI.GeoPackageConnection
  , GeoPackageTileRetriever = GeoPackageAPI.GeoPackageTileRetriever;

GeoPackageAPI.open(filename, function(err, geoPackage) {

  // Now you can operate on the GeoPackage

  // get the tile table names
  geoPackage.getTileTables(function(err, tileTableNames) {
    // tileTableNames is an array of all tile table names

    // get the info for the first table
    geoPackage.getTileDaoWithTableName(tileTableNames[0], function(err, tileDao) {
      geoPackage.getInfoForTable(tileDao, function(err, info) {
        // do something with the tile table info
      });

      // draw a tile into a canvas for an XYZ tile
      var canvas = canvasFromSomewhere;
      var gpr = new GeoPackageTileRetriever(tileDao, 256, 256);
      var x = 0;
      var y = 0;
      var zoom = 0;

      console.time('Draw tile ' + x + ', ' + y + ' zoom: ' + zoom);
      gpr.drawTileIn(x, y, zoom, canvas, function() {
        console.timeEnd('Draw tile ' + x + ', ' + y + ' zoom: ' + zoom);
      });

      // or get a tile base64 data URL for an XYZ tile
      gpr.getTile(x, y, zoom, function(err, tileBase64DataURL) {
        console.log('got the base64 data url');
      });

      // or get a tile from a GeoPackage tile column and tile row
      tileDao.queryForTile(tileColumn, tileRow, zoom, function(err, tile) {
        var tileData = tile.getTileData();  // the raw bytes from the GeoPackage
      });

    });
  });

  // get the feature table names
  geoPackage.getFeatureTables(function(err, featureTableNames) {
    // featureTableNames is an array of all feature table names

    // get the info for the first table
    geoPackage.getFeatureDaoWithTableName(featureTableNames[0], function(err, featureDao) {
      geoPackage.getInfoForTable(featureDao, function(err, info) {
        // do something with the feature table info
      });

      // query for all features
      featureDao.queryForEach(function(err, row, rowDone) {
        var feature = featureDao.getFeatureRow(row);
        var geometry = currentRow.getGeometry();
        if (geometry) {
          var geom = geometry.geometry;
          var geoJson = geometry.geometry.toGeoJSON();

          geoJson.properties = {};
          for (var key in feature.values) {
            if(feature.values.hasOwnProperty(key) && key != feature.getGeometryColumn().name) {
              var column = info.columnMap[key];
              geoJson.properties[column.displayName] = currentRow.values[key];
            }
          }
        }
        rowDone();
      });
    });
  });
});

```
# SQLite compiled to javascript
[![Build Status](https://travis-ci.org/kripken/sql.js.svg?branch=master)](http://travis-ci.org/kripken/sql.js) [![CDNJS version](https://img.shields.io/cdnjs/v/sql.js.svg)](https://cdnjs.com/libraries/sql.js)

For the impatients, try the demo here: http://kripken.github.io/sql.js/examples/GUI

*sql.js* is a port of [SQLite](http://sqlite.org/about.html) to Webassembly, by compiling the SQLite C code with [Emscripten](http://kripken.github.io/emscripten-site/docs/introducing_emscripten/about_emscripten.html). It uses a [virtual database file stored in memory](https://kripken.github.io/emscripten-site/docs/porting/files/file_systems_overview.html), and thus **doesn't persist the changes** made to the database. However, it allows you to **import** any existing sqlite file, and to **export** the created database as a [javascript typed array](https://developer.mozilla.org/en-US/docs/Web/JavaScript/Typed_arrays).

There are no C bindings or node-gyp compilation here, sql.js is a simple javascript file, that can be used like any traditional javascript library. If you are building a native application in javascript (using Electron for instance), or are working in node.js, you will likely prefer to use [a native binding of SQLite to javascript](https://www.npmjs.com/package/sqlite3).

SQLite is public domain, sql.js is MIT licensed.

Sql.js predates WebAssembly, and thus started as an [asm.js](https://en.wikipedia.org/wiki/Asm.js) project. It still supports asm.js for backwards compatability.

## Version of binaries
Sql.js was last built with:
Emscripten version 1.38.30 (2019-04-16) [Release History](https://emscripten.org/docs/introducing_emscripten/release_notes.html)
SqlLite version: 3.28.0 (2019-04-16) [Release History](https://www.sqlite.org/changes.html)

## Documentation
A [full documentation](http://kripken.github.io/sql.js/documentation/#http://kripken.github.io/sql.js/documentation/class/Database.html) generated from comments inside the source code, is available.

## Usage

```javascript
var initSqlJs = require('sql-wasm.js');
// or if you are in a browser:
//var initSqlJs = window.initSqlJs;

initSqlJs().then(function(SQL){

  // Create a database
  var db = new SQL.Database();
  // NOTE: You can also use new SQL.Database(data) where
  // data is an Uint8Array representing an SQLite database file

  // Execute some sql
  sqlstr = "CREATE TABLE hello (a int, b char);";
  sqlstr += "INSERT INTO hello VALUES (0, 'hello');"
  sqlstr += "INSERT INTO hello VALUES (1, 'world');"
  db.run(sqlstr); // Run the query without returning anything

  var res = db.exec("SELECT * FROM hello");
  /*
  [
    {columns:['a','b'], values:[[0,'hello'],[1,'world']]}
  ]
  */

  // Prepare an sql statement
  var stmt = db.prepare("SELECT * FROM hello WHERE a=:aval AND b=:bval");

  // Bind values to the parameters and fetch the results of the query
  var result = stmt.getAsObject({':aval' : 1, ':bval' : 'world'});
  console.log(result); // Will print {a:1, b:'world'}

  // Bind other values
  stmt.bind([0, 'hello']);
  while (stmt.step()) console.log(stmt.get()); // Will print [0, 'hello']

  // You can also use javascript functions inside your SQL code
  // Create the js function you need
  function add(a, b) {return a+b;}
  // Specifies the SQL function's name, the number of it's arguments, and the js function to use
  db.create_function("add_js", add);
  // Run a query in which the function is used
  db.run("INSERT INTO hello VALUES (add_js(7, 3), add_js('Hello ', 'world'));"); // Inserts 10 and 'Hello world'

  // free the memory used by the statement
  stmt.free();
  // You can not use your statement anymore once it has been freed.
  // But not freeing your statements causes memory leaks. You don't want that.

  // Export the database to an Uint8Array containing the SQLite database file
  var binaryArray = db.export();
});

```

## Demo
There are a few examples [available here](https://kripken.github.io/sql.js/index.html). The most full-featured is the [Sqlite Interpreter](https://kripken.github.io/sql.js/examples/GUI/index.html).

## Examples
The test files provide up to date example of the use of the api.
### Inside the browser
#### Example **HTML** file:
```html
<meta charset="utf8" />
<html>
  <script src='/dist/sql-wasm.js'></script>
  <script>
    config = {
      locateFile: url => `/dist/${filename}` 
    }
    // The `initSqlJs` function is globally provided by all of the main dist files if loaded in the browser.
    // We must specify this locateFile function if we are loading a wasm file from anywhere other than the current html page's folder.
    initSqlJs(config).then(function(SQL){
      //Create the database
      var db = new SQL.Database();
      // Run a query without reading the results
      db.run("CREATE TABLE test (col1, col2);");
      // Insert two rows: (1,111) and (2,222)
      db.run("INSERT INTO test VALUES (?,?), (?,?)", [1,111,2,222]);
  
      // Prepare a statement
      var stmt = db.prepare("SELECT * FROM test WHERE col1 BETWEEN $start AND $end");
      stmt.getAsObject({$start:1, $end:1}); // {col1:1, col2:111}
  
      // Bind new values
      stmt.bind({$start:1, $end:2});
      while(stmt.step()) { //
        var row = stmt.getAsObject();
        console.log('Here is a row: ' + JSON.stringify(row));
      }
    });
  </script>
  <body>
    Output is in Javscript console
  </body>
</html>
```

#### Creating a database from a file choosen by the user
`SQL.Database` constructor takes an array of integer representing a database file as an optional parameter.
The following code uses an HTML input as the source for loading a database:
```javascript
dbFileElm.onchange = () => {
  var f = dbFileElm.files[0];
  var r = new FileReader();
  r.onload = function() {
    var Uints = new Uint8Array(r.result);
    db = new SQL.Database(Uints);
  }
  r.readAsArrayBuffer(f);
}
```
See : http://kripken.github.io/sql.js/examples/GUI/gui.js

#### Loading a database from a server

```javascript
var xhr = new XMLHttpRequest();
// For example: https://github.com/lerocha/chinook-database/raw/master/ChinookDatabase/DataSources/Chinook_Sqlite.sqlite
xhr.open('GET', '/path/to/database.sqlite', true);
xhr.responseType = 'arraybuffer';

xhr.onload = e => {
  var uInt8Array = new Uint8Array(this.response);
  var db = new SQL.Database(uInt8Array);
  var contents = db.exec("SELECT * FROM my_table");
  // contents is now [{columns:['col1','col2',...], values:[[first row], [second row], ...]}]
};
xhr.send();
```
See: https://github.com/kripken/sql.js/wiki/Load-a-database-from-the-server


### Use from node.js

`sql.js` is [hosted on npm](https://www.npmjs.org/package/sql.js). To install it, you can simply run `npm install sql.js`.
Alternatively, you can simply download `sql-wasm.js` and `sql-wasm.wasm`, from the download link below.

#### read a database from the disk:
```javascript
var fs = require('fs');
var initSqlJs = require('sql-wasm.js');
var filebuffer = fs.readFileSync('test.sqlite');
 
initSqlJs().then(function(SQL){
  // Load the db
  var db = new SQL.Database(filebuffer);
});

```

#### write a database to the disk
You need to convert the result of `db.export` to a buffer
```javascript
var fs = require("fs");
// [...] (create the database)
var data = db.export();
var buffer = new Buffer(data);
fs.writeFileSync("filename.sqlite", buffer);
```

See : https://github.com/kripken/sql.js/blob/master/test/test_node_file.js

### Use as web worker
If you don't want to run CPU-intensive SQL queries in your main application thread,
you can use the *more limited* WebWorker API.

You will need to download [dist/worker.sql-wasm.js](dist/worker.sql-wasm.js) [dist/worker.sql-wasm.wasm](dist/worker.sql-wasm.wasm).

Example:
```html
<script>
  var worker = new Worker("/dist/worker.sql-wasm.js");
  worker.onmessage = () => {
    console.log("Database opened");
    worker.onmessage = event => {
      console.log(event.data); // The result of the query
    };
	
    worker.postMessage({
      id: 2,
      action: 'exec',
      sql: 'SELECT * FROM test'
    });
  };

  worker.onerror = e => console.log("Worker error: ", e);
  worker.postMessage({
    id:1,
    action:'open',
    buffer:buf, /*Optional. An ArrayBuffer representing an SQLite Database file*/
  });
</script>
```

See [examples/GUI/gui.js](examples/GUI/gui.js) for a full working example.

## Flavors/versions Targets/Downloads

This library includes both WebAssembly and asm.js versions of Sqlite. (WebAssembly is the newer, preferred way to compile to Javascript, and has superceded asm.js. It produces smaller, faster code.) Asm.js versions are included for compatibility.

## Upgrading from 0.x to 1.x

Version 1.0 of sql.js must be loaded asynchronously, whereas asm.js was able to be loaded synchronously. 

So in the past, you would:
```html
<script src='js/sql.js'></script>
<script>
  var db = new SQL.Database();
  //...
</script>
```
or:
```javascript
var SQL = require('sql.js');
var db = new QL.Database();
//...
```

Version 1.x:
```html
<script src='dist/sql-wasm.js'></script>
<script>
  initSqlJs({ locateFile: filename => `/dist/${filename}` }).then(function(SQL){
    var db = new SQL.Database();
    //...
  });
</script>
```
or:
```javascript
var initSqlJs = require('sql-wasm.js');
initSqlJs().then(function(SQL){
  var db = new SQL.Database();
  //...
});
```

`NOTHING` is now a reserved word in SQLite, whereas previously it was not. This could cause errors like `Error: near "nothing": syntax error`

### Downloading/Using: ###
Although asm.js files were distributed as a single Javascript file, WebAssembly libraries are most efficiently distributed as a pair of files, the `.js`  loader and the `.wasm` file, like [dist/sql-wasm.js]([dist/sql-wasm.js]) and [dist/sql-wasm.wasm]([dist/sql-wasm.wasm]). The `.js` file is reponsible for wrapping/loading the `.wasm` file. 




## Versions of sql.js included in `dist/`
 - `sql-wasm.js` : The Web Assembly version of Sql.js. Minified and suitable for production. Use this. If you use this, you will need to include/ship `sql-wasm.wasm` as well.
 - `sql-wasm-debug.js` : The Web Assembly, Debug version of Sql.js. Larger, with assertions turned on. Useful for local development. You will need to include/ship `sql-wasm-debug.wasm` if you use this.
 - `sql-asm.js` : The older asm.js version of Sql.js. Slower and larger. Provided for compatiblity reasons.
 - `sql-asm-memory-growth.js` : Asm.js doesn't allow for memory to grow by default, because it is slower and de-optimizes. If you are using sql-asm.js and you see this error (`Cannot enlarge memory arrays`), use this file.
 - `sql-asm-debug.js` : The _Debug_ asm.js version of Sql.js. Use this for local development.
 - `worker.*` - Web Worker versions of the above libraries. More limited API. See [examples/GUI/gui.js](examples/GUI/gui.js) for a good example of this.

## Compiling

- Install the EMSDK, [as described here](https://kripken.github.io/emscripten-site/docs/getting_started/downloads.html)
- Run `npm run rebuild`

This folder contains the neccesary files to control the global scope for Potree viewer using composition API  
# LAYOUTS

**This directory is not required, you can delete it if you don't want to use it.**

This directory contains your Application Layouts.

More information about the usage of this directory in [the documentation](https://nuxtjs.org/guide/views#layouts).
# PAGES

This directory contains your Application Views and Routes.
The framework reads all the `*.vue` files inside this directory and creates the router of your application.

More information about the usage of this directory in [the documentation](https://nuxtjs.org/guide/routing).
---
# TODO: WORK IN PROGRESS
# position: position of the label in the map
# Title and Description will be visible when hovering the label.
#
# - title: string
#   position: [x,y,x]
#   cameraPosition: [x,y,x]
#   cameraTarget: [x,y,x]
#   description: srting
##

labels:
  - title: Piramide 
    position: [296249.19705492083, 4633726.448203472, 128.7690558114301]
    cameraPosition: [296260.6092322839, 4633696.311992192, 130.7733116269969]
    cameraTarget: [296255.42565172265, 4633711.747814068, 133.44087523067427]
    description: Click on the annotation label to move a predefined view. <br>Click on the icon to execute the specified action.<br>In this case, the action will bring you to another scene and point cloud.

  - title: Easter egg
    position: [296267.3540227735,4633689.354615031,127.652587177568]
    cameraPosition:  [296268.0535866889,4633686.065772681,129.3927963214762]
    cameraTarget:  [296267.33982642664,4633687.710605324,129.18163495696592]
    description: This is an example rabbit in the scene

---
# Storylines file structure
- Each folder contains a story line, 
- Each md file represents a story page
---
title: III Curiatii

mediaPath: /videos/c_07_gec1900-1080p.mp4
mediaPosition:  [296057.1795848146,4633948.971329112,130.48376955673638]
mediaRotation:  [0.34408049273351005,-0.643849383635389,-0.6027499849171594,0.3221165028292691]
mediaScale: 1
cameraFOV: 37.36

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[296054.1931038446,4633950.96750176,130.24664897386705],[296061.46748833,4633946.105281812,130.8242204724034]],
    [[296051.5405588594,4633954.793215068,129.95335636457568],[296069.124996751,4633940.556095265,130.56847420441596]],
    [[296047.2815013388,4633965.57257587,129.28770062308524],[296060.25982305116,4633947.090213253,130.788315106573]],
    [[296046.58604174503,4633985.16486595,130.32754795980398],[296057.4844652355,4633965.332025003,130.73887700145636]],
    [[296048.09380919614,4634010.3055092925,152.09152780893103],[296053.69075344194,4633989.114253503,146.4442895991443]],
    [[296013.62269847124,4634021.973557029,164.72647153868942],[296027.7679362976,4634006.352478239,156.46932006500464]]

]

# ViewPoints overrides in seconds
animationEntry:
---
---
title: III Curiatii

mediaPath: /videos/c_02_cl-1080p.mp4
mediaPosition:  [295976.2124383994,4633832.429925773,132.39068275360825]
mediaRotation:  [0.502784015144051,-0.5055979127976997,-0.4971531320061001,0.4943862336541762]
mediaScale: 1
cameraFOV: 32

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[295972.6130050703,4633832.450014424,132.33005144337795],[296000.5515015024,4633832.294088024,132.800666534611]],
    [[295896.71896302194,4633842.749220234,155.73519598554955],[295911.0923866786,4633844.079734429,152.24851785827786]],
    [[295961.21334609983,4633919.326354789,158.4163485137719],[295969.7455074782,4633890.99621932,149.00556220045902]],
    [[296100.84482554445,4633982.914183403,159.95245641293613],[296095.3139897102,4633953.26178987,152.5970050203098]]
]


# ViewPoints overrides in seconds
animationEntry:
---
---
title: III Curiatii

mediaPath: /videos/c_13_er1930-1080p.mp4
mediaPosition:  [296136.60576125083,4633858.933033093,128.85051979602719]
mediaRotation:  [-0.7042249803027177,-0.258333958148526,-0.2277507988077569,-0.620854505361265]
mediaScale: 1
cameraFOV: 36.54

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[296138.9153455968,4633856.208657591,128.39930114695883],[296133.85228836234,4633862.181017562,129.3884600185369]],
    [[296141.35309367836,4633851.489688265,127.74540032136855],[296126.83238734415,4633868.618249384,130.58228025906027]],
    [[296144.25275771745,4633845.4056430515,130.17600731063808],[296128.9493674441,4633861.924004673,127.88769315540343]],
    [[296156.0758040967,4633829.367951202,132.18812194731777],[296140.77241382335,4633845.886312824,129.8998077920831]],
    [[296175.955758152,4633807.6969803935,130.09379344173587],[296159.4084735207,4633823.135119273,129.7227874199425]]
]

# ViewPoints overrides in seconds
animationEntry:
---
---
title: III Curiatii

mediaPath: /videos/c_04_rmacp1850s-1080p.mp4
mediaPosition:  [296040.84433442773,4633987.852257391,130.53462988422825]
mediaRotation:  [0.27170561012487443,-0.6731846028580056,-0.6377613193782398,0.2574083359304938]
mediaScale: 1
cameraFOV: 32

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[296038.34905449883,4633990.439878435,130.34021945791963],[296054.03725706966,4633974.171113304,131.56250723192224]],
    [[296032.69593864185,4633996.369619758,131.00823515080447],[296048.4047284712,4633980.079505497,130.61771543984256]],
    [[296025.143152492,4634004.201906695,131.19599702783182],[296040.85194232134,4633987.911792434,130.8054773168699]],
    [[296021.617166337,4634004.881394485,131.78830584202422],[296037.32265479164,4633988.594703769,131.3947249627685]],
    [[296003.74484399194,4634023.415105633,131.0982210965927],[296019.4503324466,4634007.128414918,131.704640217337]]
]


# ViewPoints overrides in seconds
animationEntry:
---
---
title: III Curiatii

mediaPath: /videos/C_10_1905a-1080p.mp4
mediaPosition:  [296139.24153037503,4633857.854031051,129.78717226889756]
mediaRotation:  [-0.6642303061933994,-0.34533112335013494,-0.30582049889614643,-0.5882332343848559]
mediaScale: 1
cameraFOV: 35.28

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[296142.1666778832,4633855.801217414,129.35189278107325],[296135.0416876808,4633860.80140185,130.41213404323042]],
    [[296158.69743095664,4633841.17023041,128.2679580161413],[296131.2297026,4633855.076104793,132.28141152119053]],
    [[296152.21405617526,4633832.937157529,128.39709013424306],[296144.80712888454,4633854.323611106,128.18867295286697]],
    [[296158.533516294,4633826.486268551,127.95590259503611],[296144.5687625713,4633844.251278684,126.6611150949184]],
    [[296169.0936630589,4633813.052367885,128.93502090260966],[296155.1289093362,4633830.817378017,127.64023340249196]]
]


# ViewPoints overrides in seconds
animationEntry:
---
---
title: III Curiatii

mediaPath: /videos/C_11_1905d-1080p.mp4
mediaPosition:  [296067.9450247373,4633950.447018224,128.2258335169477]
mediaRotation:  [-0.23555199323674567,0.6775802680368137,0.6580758494300156,-0.22877153504384926]
mediaScale: 1
cameraFOV: 37.79

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[296065.7128652136,4633953.269501884,128.12071523666305],[296071.1498952068,4633946.394577332,128.3767593216298]],
    [[296062.2951367346,4633957.591093899,127.95976538681745],[296076.3290614996,4633939.845711956,128.62065992471307]],
    [[296059.29740369826,4633961.224311991,129.14667911055477],[296074.6188838128,4633945.123412672,124.86886489634765]],
    [[296054.56631680916,4633966.259749442,130.6638242131597],[296068.8344282612,4633950.226489427,123.47758910319797]],
    [[296047.12586743286,4633975.122295576,134.34999016699803],[296062.21462414734,4633959.123400329,128.9973106614901]],
    [[296037.6092733908,4633985.30484699,137.00971077831696],[296052.53950982186,4633969.150815633,131.6785085458364]]
]


# ViewPoints overrides in seconds
animationEntry:
---
---
title: III Curiatii

mediaPath: /videos/c_09_anonearly1900-1080p.mp4
mediaPosition:  [296089.3683794771,4633905.102576009,130.96776925540894]
mediaRotation:  [-0.31246392186229144,0.6900580493537798,0.5947025641098184,-0.2692861791806252]
mediaScale: 1
cameraFOV: 37.79

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[296086.6925269016,4633907.4514829535,130.4363097950157],[296093.21029167436,4633901.730082948,131.7308234691821]],
    [[296084.4166496632,4633914.562320727,130.76776117675664],[296096.1411726965,4633895.532788923,127.20421165770159]],
    [[296086.4750041632,4633924.4003962,129.24621816285838],[296095.661415272,4633903.754027389,130.52046776951408]],
    [[296081.0609766165,4633938.663689658,128.38728897802895],[296090.2473877253,4633918.017320847,129.66153858468465]],
    [[296073.5402249991,4633949.015720067,127.77668267358175],[296084.482180544,4633929.218915143,128.58075348412544]]
]

# ViewPoints overrides in seconds
animationEntry:
---
---
title: III Curiatii

mediaPath: /videos/c_03_lc-1080p.mp4
mediaPosition:  [296093.5433670494,4633912.739320675,128.06704018933536]
mediaRotation:  [0.11738557061262322,-0.7147115073713547,-0.6803827014798023,0.11174734256330754]
mediaScale: 1
cameraFOV: 34.87

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[296092.3932806417,4633916.146071383,127.88997843044434],[296092.9123759718,4633914.608423105,127.96989584496512]],
    [[296090.10257355927,4633922.931533701,127.53731221131076],[296097.3333408661,4633901.512776389,128.6505262870782]],
    [[296088.23120295594,4633926.674408419,129.38903397803682],[296098.10837503197,4633906.352155133,128.0722578528838]],
    [[296077.6397793022,4633940.131836061,130.74531853931487],[296089.2355255331,4633920.893563375,127.96792055104254]],
    [[296068.29729128105,4633953.138598328,133.31007244387015],[296080.52718462795,4633934.590383026,128.98805624642776]],
    [[296068.29729128105,4633953.138598328,133.31007244387015],[296080.52718462795,4633934.590383026,128.98805624642776]]
]

# ViewPoints overrides in seconds
animationEntry:
---
---
title: III Curiatii

mediaPath: /videos/c_08_c-gca1900-1080p.mp4
mediaPosition:  [296044.5999809695,4633976.48482735,130.48376955673703]
mediaRotation:  [0.3440804926561288,-0.6438493836728728,-0.6027499849582837,0.3221165027600516]
mediaScale: 1
cameraFOV: 37.36

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[296041.6135,4633978.481,130.24664897390358],[296048.88788448425,4633973.6187800495,130.82422047235252]],
    [[296037.9440840249,4633982.224982695,129.90793626208287],[296056.7205418804,4633969.674743126,131.3987492542583]],
    [[296029.3222680948,4633995.903933822,128.98542430602208],[296046.58753065986,4633981.270944224,128.71573898492568]],
    [[296038.243835136,4633996.372456005,141.68299544864328],[296050.22147078824,4633977.994473929,136.10873090607487]],
    [[295995.3428236987,4633994.197745317,149.71501689768402],[296014.01919273945,4633984.109624278,141.85948002404314]]
]

# ViewPoints overrides in seconds
animationEntry:
---

---
cameraFOV: 60

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
  [[296145.23893221654,4633909.824099188,159.40584846815443],[296064.9830625395,4633848.548224539,72.33093620244487]],
  [[296173.7993210386,4633931.630140443,190.39290724648816],[296096.1661795595,4633872.356735706,99.63203124424209]]
]

---
---
title: III Curiatii

mediaPath: /videos/c_14_anon1930-1080p.mp4
mediaPosition:  [296145.37817308365,4633846.381874894,128.2361403364562]
mediaRotation:  [0.6775370344582602,0.24604955359439326,0.2365877638869976,0.6514824741253263]
mediaScale: 1
cameraFOV: 37.36

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[296147.6864494795,4633843.622898327,128.09504326207346],[296142.62625952566,4633849.671110548,128.40435543959003]],
    [[296177.2649246794,4633806.603293696,133.7059685038686],[296156.1632534752,4633829.377610476,133.79093951936846]],
    [[296183.29578327097,4633794.394878338,136.7314102096456],[296162.5517123532,4633817.130917485,132.6432201264756]]
]

# ViewPoints overrides in seconds
animationEntry:
---
---
title: III Curiatii

mediaPath: /videos/c_01_gbp-1080p.mp4
mediaPosition:  [296219.26107785135,4633867.709464002,133.40150829467254]
mediaRotation:  [-0.4194067254941228,-0.5711270301585142,-0.5687468650954337,-0.4176588529849806]
mediaScale: 1
cameraFOV: 36.25

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[296222.696,4633868.78699999,133.38647407570735],[296206.9526068186,4633863.84829338,133.4553809126311]],
    [[296229.61941477866,4633869.724579454,135.32005973582022],[296208.2099419223,4633864.28994283,130.3818816332365]],
    [[296241.7075221199,4633855.652041152,138.82855831338284],[296221.89723228,4633865.524016445,134.09749854885348]],
    [[296225.68727121694,4633817.128920737,141.77270486968803],[296224.83744182426,4633839.234306215,136.98562548184827]],
    [[296174.44860012224,4633799.862267898,130.47716422287382],[296164.88698015624,4633820.333983268,129.14648782118329]],
    [[296178.23350603017,4633793.82167212,130.2544812067517],[296166.2217449762,4633812.992080812,130.96118706563266]]
]



# ViewPoints overrides in seconds
animationEntry:
---
---
title: III Curiatii

mediaPath: /videos/C_12_1905e-1080p.mp4
mediaPosition:  [296064.7540066245,4633947.595149822,129.2183269508244]
mediaRotation:  [0.25506772306293124,-0.6704776042566997,-0.6511776386532632,0.2477254967895885]
mediaScale: 1
cameraFOV: 37.79

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[296062.36224729905,4633950.283725183,129.11320867048457],[296073.3244775407,4633937.961088109,129.59500078870886]],
    [[296053.81820304604,4633966.539379014,130.403355769477],[296066.92555642535,4633948.335476715,127.386464209164]],
    [[296046.7716898693,4633976.325796892,132.02523667231333],[296059.8790432486,4633958.121894593,129.00834511200034]]
]


# ViewPoints overrides in seconds
animationEntry:
---
---
title: III Curiatii

mediaPath: /videos/c_15_tv1956-1080p.mp4
mediaPosition:  [296066.8459923548,4633954.588549161,128.38105519411164]
mediaRotation:  [0.21752443437896457,-0.690390597657515,-0.6580682798804119,0.20734049804475435]
mediaScale: 1
cameraFOV: 30.71

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[296064.7846909576,4633957.5349513665,128.20857158453396],[296068.444896936,4633952.303092027,128.5148468045497]],
    [[296061.2116978462,4633962.642149563,128.37735813943152],[296074.32042530016,4633944.192978825,128.65262133682424]],
    [[296054.68564004335,4633970.03442264,128.25805906967716],[296068.6034934804,4633952.187742498,128.53332226706988]],
    [[296043.37893886096,4633983.258215391,129.2304462834397],[296058.0648736345,4633966.082227923,127.96744171954509]],
    [[296040.21956430847,4633996.500916561,132.768040297091],[296054.7659646728,4633979.298577891,130.58442244805536]],
    [[296019.0447660271,4634004.897660223,134.84102150844217],[296037.04872827296,4633991.47630403,132.01065976311077]]
]

# ViewPoints overrides in seconds
animationEntry:
---
---
title: III Curiatii

mediaPath: /videos/c_06_da1890-1080p.mp4
mediaPosition:  [296032.791473029,4633968.792806959,130.82106321979043]
mediaRotation:  [0.35653648295147217,-0.6062043201722179,-0.6127846859246004,0.3604066969763776]
mediaScale: 1
cameraFOV: 36.95

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[296029.64535963617,4633970.542223991,130.8599292412002],[296031.06536839315,4633969.752618788,130.84238693300654]],
    [[296025.7946523094,4633976.549944257,130.9277870576215],[296045.57474337105,4633965.551094606,130.6834304986848]],
    [[296022.75907071476,4633991.314492212,130.82296247778547],[296035.187277969,4633972.450209712,129.42120771171926]],
    [[296017.30180922133,4634007.907157232,130.5294847760481],[296030.0042179403,4633989.339731343,128.04233484161207]],
    [[295977.0960663542,4634061.573542704,143.43102184995203],[295989.80312601355,4634044.071448778,136.76065834881584]]
]


# ViewPoints overrides in seconds
animationEntry:
---
---
title: III Curiatii

mediaPath: /videos/c_05_rmacp1862-1080p.mp4
mediaPosition:  [296063.61795906385,4633945.939122305,130.16304082206082]
mediaRotation:  [0.2993553585400927,-0.6354129956628881,-0.6439038041198369,0.3033555427152267]
mediaScale: 1
cameraFOV: 38.62

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[296060.84227188415,4633948.231125218,130.21082502070436],[296068.3544263634,4633942.0280198855,130.081501268193]],
    [[296046.2442431575,4633967.688935353,129.5694612217775],[296060.7151864137,4633950.2892494565,129.21010322984154]],
    [[296025.93610639294,4634004.632850104,148.1462969833519],[296046.2067208105,4633982.620488797,139.86889441451032]]

]


# ViewPoints overrides in seconds
animationEntry:
---
---
title: I Mausoleo Rotondo

mediaPath: /videos/mr_04_erf-1080p.mp4
mediaPosition:  [296009.6728270726,4634021.916551706,128.96809265104994]
mediaRotation:  [-0.17749404887464218,0.7113776189048389,0.6598002140657462,-0.16462509970889433]
mediaScale: 1
cameraFOV: 37.7

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[296007.986434268,4634025.085611041,128.6976444218245],[296015.7157346223,4634010.560755759,129.93719880575267]],
    [[296008.0257193365,4634025.128638689,128.82057641799273],[296015.49195854657,4634010.48812509,130.29048579787462]],
    [[295980.26934436645,4634054.82205554,131.57193854818877],[295993.3046320715,4634044.809546258,130.12939484503102]],
    [[295972.3111422813,4634065.561850869,132.55532631503362],[295984.43900167756,4634054.424956594,131.49014446686178]]
]


# ViewPoints overrides in seconds
animationEntry:
---
---
title: I Mausoleo Rotondo

mediaPath: /videos/mr_09_1905a-1080p.mp4
mediaPosition:  [296097.0021353372,4633921.222725828,128.27621883624673]
mediaRotation:  [0.6860624997349487,0.19728023952813836,0.1935276506988402,0.6730124827699271]
mediaScale: 1
cameraFOV: 46

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[296098.91405305604,4633918.173169824,128.20708986180063],[296096.6866271612,4633921.725969165,128.28762662666776]],
    [[296100.5226523999,4633912.501158718,128.09838166907707],[296091.7596961872,4633926.478290402,128.41522280200874]],
    [[296102.968569152,4633905.016277739,128.6352814177397],[296095.6357644681,4633919.786108885,129.2115308273205]],
    [[296110.8554894431,4633889.130339686,130.33043676130714],[296103.6295143943,4633903.963195351,130.4776692468356]],
    [[296117.55445004534,4633875.379293628,130.1939424317208],[296110.32847499655,4633890.212149293,130.34117491724928]]
]

# ViewPoints overrides in seconds
animationEntry:
---
---
title: I Mausoleo Rotondo

mediaPath: /videos/mr_15_bfr-1080p.mp4
mediaPosition:  [295982.8642286003,4634052.075568891,129.45580414927727]
mediaRotation:  [0.22403910365583485,-0.7096611540243712,-0.6369837879471411,0.2010949536771107]
mediaScale: 1
cameraFOV: 30

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[295980.809215013,4634055.005893413,129.06835434009085],[295981.3569164697,4634054.2249044115,129.171617319077]],
    [[295939.2624261372,4634105.899806225,130.45357127652252],[295949.85826453165,4634093.281273568,129.58705970358204]],
    [[295929.9394382801,4634115.525618284,132.04298182627085],[295941.22170392604,4634103.632816568,130.16634940832978]]
]

# ViewPoints overrides in seconds
animationEntry:
---
---
title: I Mausoleo Rotondo

mediaPath: /videos/mr_03_fllid'a-1080p.mp4
mediaPosition:  [296030.3479579428,4633989.401942355,129.3487837832159]
mediaRotation:  [0.6848898682181741,0.2659395982258716,0.2455509520438992,0.6323817901061737]
mediaScale: 1
cameraFOV: 52.62

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[296032.76968311507,4633986.753709207,129.0622381724703],[296021.6701094088,4633998.891444466,130.37557222172103]],
    [[296037.8362685251,4633981.21323765,128.46274495688309],[296026.7366948198,4633993.3509729095,129.77607900603635]],
    [[296053.13108154526,4633966.524055845,129.79123801426908],[296041.8111535942,4633978.521461204,130.20535857152004]],
    [[296064.6973013733,4633954.50364016,136.78241103919896],[296053.4432303969,4633966.14433475,133.60591675844748]]
]



# ViewPoints overrides in seconds
animationEntry:
---
---
title: I Mausoleo Rotondo

mediaPath: /videos/mr_13_fllia2-1080p.mp4
mediaPosition:  [296040.2432447413,4633999.123614633,131.71465940635144]
mediaRotation:  [-0.523427624965146,-0.5319466352519582,-0.47445155605959166,-0.4668533170280562]
mediaScale: 1
cameraFOV: 36

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[296043.8193558778,4633999.181351332,131.30466475675559],[296027.428846502,4633998.916724797,133.18380690074477]],
    [[296043.82021534175,4633999.128117627,131.30466475675559],[296042.2061256712,4633999.102057977,131.48971720720232]],
    [[296053.42145264463,4633999.566055508,139.77482148537152],[296039.86874020024,4633999.02307149,133.72902511468592]],
    [[296062.9711250356,4633956.911539005,135.30492483391703],[296051.9874058906,4633969.128696553,133.77228815563413]]
]

# ViewPoints overrides in seconds
animationEntry:
---
---
title: I Mausoleo Rotondo

mediaPath: /videos/mr_02_da2-1080p.mp4
mediaPosition:  [296042.48715267057,4633990.733943746,130.46400905933385]
mediaRotation:  [0.6209115935397026,0.3522940928127215,0.3455644810433871,0.60905078164178]
mediaScale: 1
cameraFOV: 31

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[296045.5768885639,4633988.887673911,130.39458424861004],[296044.3217822376,4633989.637661945,130.42278585617484]],
    [[296052.2043758284,4633982.624702855,130.22288512089855],[296038.04308631783,4633991.0867729345,130.5410821700314]],
    
    [[296061.25794632756,4633963.6478963215,130.42843686475624],[296051.2172902531,4633976.668951353,129.05462775856853]],
    [[296064.42448994843,4633955.506522902,131.27374431955624],[296053.57541442115,4633967.778000272,129.28417065455943]]
]

# ViewPoints overrides in seconds
animationEntry:
---
---
title: I Mausoleo Rotondo

mediaPath: /videos/mr_12_fllia1-1080p.mp4
mediaPosition:  [296013.22135282785,4633979.595465134,132.29563047113066]
mediaRotation:  [0.7442793156478824,0.09943023276462744,0.0874510715554438,0.6546099901429069]
mediaScale: 1
cameraFOV: 37

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[296014.158620369,4633976.15013989,131.83599636073052],[296013.2741512747,4633979.401382,132.26973821937037]],
    [[296016.65464327636,4633966.320056503,133.44331961789493],[296011.08020967635,4633981.612027617,130.73568052131066]],
    [[296027.5593892429,4633945.478441913,140.42032948541095],[296020.69630880177,4633959.641371881,135.46457105818214]],
    [[296038.66942258796,4633951.435088397,141.0980939157392],[296028.02891503234,4633962.766521643,135.5638171299014]],
    [[296073.8124128873,4633923.97395449,159.12897251896416],[296061.9871914993,4633933.382590545,152.50389088057034]],
    [[296094.50416474586,4633923.588061362,166.62794320380664],[296081.85089346196,4633931.471850488,159.55754046933592]]
]

# ViewPoints overrides in seconds
animationEntry:
---
---
title: I Mausoleo Rotondo

mediaPath: /videos/mr_19_en-1080p.mp4
mediaPosition:  [295980.21999169124,4633962.508446806,133.49356248609092]
mediaRotation:  [0.6844462457794682,-0.1788803696257808,-0.17871414178156828,0.6838102116290845]
mediaScale: 1
cameraFOV: 38

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[295978.45858047414,4633959.368794067,133.49021556609767],[295980.9920157008,4633963.884552098,133.49502943632254]],
    [[295963.40121988434,4633930.819013654,143.11703266636803],[295969.71979837416,4633943.621858239,139.0322760760358]],
    [[296024.1716123595,4633924.230682101,141.76404148392706],[296014.3806082896,4633937.103772519,138.49757819789656]],
    [[296051.28165102453,4633928.155989705,144.94414411815052],[296038.8151492209,4633938.5518509215,141.9840085919256]]
]

# ViewPoints overrides in seconds
animationEntry:
---
---
title: I Mausoleo Rotondo

mediaPath: /videos/mr_08_anon1900-1080p.mp4
mediaPosition:  [296040.31980402797,4633981.29109547,129.64225267784576]
mediaRotation:  [0.6961142301098672,0.26503801932881554,0.23741042953020644,0.6235512127316315]
mediaScale: 1
cameraFOV: 38

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[296042.69961683656,4633978.618886877,129.24754743187717],[296040.4538643838,4633981.140563784,129.62001801908463]],
    [[296048.0602199125,4633972.650827083,130.43431567070618],[296036.90397848294,4633984.582734382,128.10668482264404]],
    [[296053.4962844862,4633966.836806388,130.41514259530055],[296042.26902446884,4633978.844670022,128.99745311871214]],
    [[296066.04568705655,4633953.414874896,130.47217076404888],[296054.78872473526,4633965.454505962,131.23248852928788]],
    [[296082.72849214653,4633934.43679741,133.65211807501754],[296071.27805419563,4633946.316587659,133.5623718327131]]
]


# ViewPoints overrides in seconds
animationEntry:
---
---
title: I Mausoleo Rotondo

mediaPath: /videos/mr_10_1905g-1080p.mp4
mediaPosition:  [296046.6725875602,4633977.098428775,129.31704996722706]
mediaRotation:  [0.7083788608819993,0.2870507157804339,0.24217050143332672,0.5976242333280946]
mediaScale: 1
cameraFOV: 48

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[296049.142885441,4633974.550858838,128.7108191546104],[296046.26493446354,4633977.518833447,129.4170912901677]],
    [[296054.4928788518,4633967.861592617,129.45999242590372],[296042.58575715177,4633979.246831746,130.380200441229]],
    [[296064.2095849235,4633956.876861346,129.84212589537393],[296052.3024632235,4633968.262100475,130.7623339106992]],
    [[296068.68828700663,4633951.970300026,131.64174152789494],[296057.2333678569,4633963.803325762,130.63460240426085]],
    [[296078.3149349177,4633939.333375189,132.61568348066254],[296066.9755845459,4633951.254799324,131.37118368129387]]
]

# ViewPoints overrides in seconds
animationEntry:
---
---
title: I Mausoleo Rotondo

mediaPath: /videos/mr_11_s-tco1908-1080p.mp4
mediaPosition:  [296065.45426688873,4633943.081264223,130.28362654266823]
mediaRotation:  [0.6207631821208122,0.19884752139041117,0.23134483983959833,0.722213472624074]
mediaScale: 1
cameraFOV: 38

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[296067.52225605753,4633940.184551688,130.82443823438118],[296058.0439723671,4633953.461150807,128.34571798069678]],
    [[296072.9129399887,4633934.5176018225,130.5330012706072],[296061.85118996823,4633946.384015229,127.52037461528114]],
    [[296084.410593079,4633928.350027186,132.1633349611518],[296072.0466296449,4633938.674946866,128.58864031597997]],
    [[296093.0616147225,4633921.164116482,132.54819431367278],[296080.37879164465,4633931.555104481,130.6979599527207]],
    [[296103.4690919357,4633903.555752718,135.77018254629937],[296092.51713822153,4633915.333614937,132.0841723194061]]
]


# ViewPoints overrides in seconds
animationEntry:
---
---
title: I Mausoleo Rotondo

mediaPath: /videos/mr_14_cinec-1080p.mp4
mediaPosition:  [296104.168593022,4633896.786586147,128.83772192321757]
mediaRotation:  [0.6815823519865024,0.15869795407481999,0.16198869203961738,0.6957155456681328]
mediaScale: 1
cameraFOV: 36

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[296105.7584773481,4633893.557528967,128.91159721244443],[296098.4715075199,4633908.357374375,128.57300213681012]],
    [[296105.8175163313,4633893.586597872,128.91159721244443],[296105.09991674614,4633895.044043584,128.8782533546337]],
    [[296110.75493813935,4633886.038582075,128.92615571965214],[296101.7714903358,4633899.878448172,129.0021293971797]],
    [[296115.8771047575,4633874.219908515,129.32836683588192],[296107.71472862334,4633888.555447028,128.9849263334182]],
    [[296126.4304869763,4633863.162522158,136.2654908398502],[296116.77428594034,4633875.588042754,131.30421887943623]],
    [[296133.8807069317,4633836.428017561,145.068538364793],[296126.5334852799,4633850.206943192,139.73847809662036]]
]


# ViewPoints overrides in seconds
animationEntry:
---



[[296105.7584773481,4633893.557528967,128.91159721244443],[296098.4715075199,4633908.357374375,128.57300213681012]]
---
title: I Mausoleo Rotondo

mediaPath: /videos/mr_16_doc1164-1080p.mp4
mediaPosition:  [296010.58031073265,4633981.964168431,132.9535872228128]
mediaRotation:  [0.7765819254364801,0.027643577694282728,0.022390578001031832,0.6290111348095591]
mediaScale: 1
cameraFOV: 37

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[296010.8306996344,4633978.451578429,132.20591291183575],[296010.7176857698,4633980.036997624,132.54337820197594]],
    [[296011.4139566698,4633971.86418139,134.40846782573232],[296009.5749225757,4633988.102235678,132.1295324273154]],
    [[296023.2088996643,4633956.510872779,140.24952857946093],[296015.2805507621,4633970.058237101,135.16423570424183]],
    [[296042.61810364615,4633943.244924738,143.43465621752998],[296031.26182869857,4633954.945519635,140.90857261717477]],
    [[296069.1971344336,4633948.008029201,148.58278600052944],[296056.23160663055,4633957.186024958,144.12078166731726]]
]

# ViewPoints overrides in seconds
animationEntry:
---
---
title: I Mausoleo Rotondo

mediaPath: /videos/mr_01_da1-1080p.mp4
mediaPosition:  [296030.00320750225,4633986.913381339,130.38025388346102]
mediaRotation:  [0.7155481052417728,0.2335493785750627,0.20427969664465953,0.6258717140044568]
mediaScale: 1
cameraFOV: 44

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
  [[296032.10808358045,4633984.032431494,129.90106218079623],[296024.4309698231,4633994.540118819,131.64881790432258]],
  [[296038.3262287553,4633978.926977197,129.88155718266722],[296026.269332394,4633990.180985145,129.40348267635943]],
  [[296047.78487469134,4633971.969497356,130.22073510909095],[296033.8789104232,4633980.837732949,129.74266060278316]],
  [[296056.34083757777,4633964.862448437,130.54061476126336],[296042.55832396867,4633973.921353634,130.06254025495556]],
  [[296065.3951435988,4633955.3442426175,133.89494287443372],[296052.52342077193,4633965.224289522,130.90271550162615]]
]

# ViewPoints overrides in seconds
animationEntry:
---
---
cameraFOV: 60

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
  [[296042.24160671013,4634028.000173633,163.82310337277625],[295964.97084735223,4633950.327378622,87.83662798872214]],
  [[296073.7401086453,4634059.80243414,196.5547448611831],[295998.9193766951,4633984.13185226,116.22045027836168]]
]

---
---
title: I Mausoleo Rotondo

mediaPath: /videos/mr_18_jhl-1080p.mp4
mediaPosition:  [296003.2116087451,4634028.404700467,129.8215907952377]
mediaRotation:  [0.189317771527056,-0.7548728139313343,-0.609090546482942,0.15275641510761365]
mediaScale: 1
cameraFOV: 38

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[296001.5511207703,4634031.50693848,129.06073667167192],[296005.19123275334,4634024.706230534,130.72867661133486]],
    [[295997.70091201254,4634036.33334773,128.86362635662937],[296007.3624705686,4634022.96461058,129.2892474013709]],
    [[295991.9168770954,4634041.636038262,128.6615414845894],[296003.3426050852,4634029.73975369,129.08716252933093]],
    [[295980.6717140627,4634053.344321505,131.5038987524547],[295992.02098132856,4634041.527646611,129.5518928674254]],
    [[295966.35638606333,4634068.978565416,134.4797244603222],[295977.2855622522,4634057.360635615,130.25754315714]]
]

# ViewPoints overrides in seconds
animationEntry:
---
---
title: I Mausoleo Rotondo

mediaPath: /videos/mr_06_anon2-1080p.mp4
mediaPosition:  [296002.2499244519,4634016.404459633,130.00663866617558]
mediaRotation:  [0.1800624865495736,-0.7121361845468546,-0.6578536892714215,0.1663372451034737]
mediaScale: 1
cameraFOV: 43

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[296000.5441757484,4634019.561877925,129.72180346647133],[296005.16063536145,4634011.016601592,130.49268501139088]],
    [[295998.26169601333,4634026.046749059,129.18285175231927],[296006.07971090236,4634011.575248553,130.48834641735453]],
    [[295997.53539133974,4634032.400650265,131.4737465115422],[296005.2601623748,4634017.831399955,130.9126245442714]],
    [[295989.93344946334,4634041.434482504,131.62391191728943],[296000.06393662526,4634028.421751441,131.0834860191703]],
    [[295979.8789893025,4634055.353713452,132.08984289757586],[295990.3493174895,4634042.603373042,131.86390682810077]]
]

# ViewPoints overrides in seconds
animationEntry:
---
---
title: I Mausoleo Rotondo

mediaPath: /videos/mr_05_anon1-1080p.mp4
mediaPosition:  [295987.6396659968,4634039.563619223,129.32914461780814]
mediaRotation:  [0.2206033365679761,-0.7015248790932908,-0.6464323215599567,0.2032787877542191]
mediaScale: 1
cameraFOV: 45

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[295985.58615216793,4634042.505858741,129.03536310354733],[295991.1438076632,4634034.542944933,129.83045706601172]],
    [[295982.25148089736,4634048.754986823,128.45951053246972],[295991.66341928067,4634035.2697223695,129.8060091397197]],
    [[295979.21164092084,4634054.393820027,129.7000585373195],[295987.2284612537,4634039.975609586,130.00957337570367]],
    [[295968.4960547576,4634067.798572968,129.24309885786016],[295979.28270876815,4634055.3308225265,129.9163701125969]],
    [[295958.3450007393,4634081.942554186,129.82424192182546],[295969.1355956547,4634069.470248651,130.3288764516245]]

]

# ViewPoints overrides in seconds
animationEntry:
---
---
title: I Mausoleo Rotondo

mediaPath: /videos/mr_17_mh-1080p.mp4
mediaPosition:  [295983.46563318034,4634043.170043939,128.96309737686892]
mediaRotation:  [0.22063612588327808,-0.6640197934168958,-0.6779745730660185,0.22527292820432604]
mediaScale: 1
cameraFOV: 36

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[295981.31159934157,4634046.053537305,129.0379587569907],[295982.95807997027,4634043.849479093,128.98073690188917]],
    [[295978.8603013871,4634052.163665808,129.17028664375528],[295987.977479846,4634038.415212806,128.84302985954326]],
    [[295970.35514604556,4634068.632118452,129.450543452109],[295978.78766659845,4634054.45709999,128.99070792996488]],
    [[295961.91699912836,4634079.384549949,129.87861766782842],[295972.2386913286,4634066.5209946865,129.38707303973536]],
    [[295953.2684532336,4634090.162922903,130.29048291170562],[295963.59014543385,4634077.29936764,129.79893828361256]]
]

# ViewPoints overrides in seconds
animationEntry:
---
---
title: I Mausoleo Rotondo

mediaPath: /videos/mr_07_c-gca1900-1080p.mp4
mediaPosition:  [296005.989614073,4634041.921464165,130.8960383037013]
mediaRotation:  [0.08872389018247133,-0.7196990514115831,-0.6834201046864746,0.08425145232270059]
mediaScale: 1
cameraFOV: 31.79

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[296005.11646013247,4634045.4090122925,130.71],[296009.11841569335,4634029.424416706,131.56267555863084]],
    [[296004.2818564825,4634048.460495357,131.07936615260334],[296008.9528897919,4634032.751579811,129.16588237129267]],
    [[296005.04362361634,4634051.320861584,136.08846913063977],[296007.34066620347,4634037.431626258,131.36245205214476]],
    [[296003.33200849366,4634046.782836024,130.56238127932076],[296009.70874486276,4634031.564849828,130.55100423675984]],
    [[296000.58830107213,4634044.745659683,130.5566870017558],[296014.81669796444,4634036.396570196,130.86559898532175]],
    [[295993.61090716755,4634037.803141222,131.33706232070531],[296007.6925537786,4634046.037426472,128.8563819405876]],
    [[295978.9908722678,4634029.254031798,133.91258738067967],[295993.07251887885,4634037.488317047,131.43190700056195]],
    [[295962.4898378049,4634045.473702696,140.90768491899988],[295977.05515873287,4634040.83017228,134.69970203125675]],
    [[295933.0864964038,4634058.40603086,159.22376982702093],[295947.1069218871,4634054.252897819,151.5798732283207]]

]



# ViewPoints overrides in seconds
animationEntry:
---
---
title: I Mausoleo Rotondo

mediaPath: /videos/mr_20_anon1973-1080p.mp4
mediaPosition:  [296042.1682328262,4633978.837745597,129.50119198080546]
mediaRotation:  [0.7246838282206157,0.2624343128517174,0.2169480835705427,0.5990785504327701]
mediaScale: 1
cameraFOV: 32

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[296044.4321830815,4633976.121848601,128.8241153607911],[296039.4691646326,4633982.07562204,130.30839874936146]],
    [[296051.38462544145,4633966.716689559,130.3978214415481],[296040.44242838543,4633979.00806983,129.19769869297207]],
    [[296063.7679415044,4633953.737162919,132.37349915984862],[296052.3798276281,4633965.507849626,130.37055213625797]],
    [[296081.5638071404,4633933.671894817,134.1456634172921],[296069.8679071175,4633945.136803271,132.14271639370145]],
    [[296098.92645517405,4633911.017751159,136.05333750720698],[296088.14738160226,4633923.378464628,134.24382501093504]]
]

# ViewPoints overrides in seconds
animationEntry:
---
---
title: II Lataritio

mediaPath: /videos/l_09_uu-1080p.mp4
mediaPosition:  [295929.81155456527,4634129.472590153,127.4764834460679]
mediaRotation:  [0.277929048016758,-0.5397106771502349,-0.7064784552321997,0.3638076710725855]
mediaScale: 1
cameraFOV: 60

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[295926.98410182836,4634131.489895123,128.42305181658534],[295939.94326020585,4634122.243914012,124.08461345171374]],
    [[295921.2760280912,4634135.562439068,130.33398811120176],[295939.0526513784,4634122.879337405,124.38276950424134]],
    [[295916.57897872565,4634136.940673029,129.40745878653624],[295937.0060380178,4634127.895732862,125.77280443590261]],
    [[295908.756289069,4634140.732458989,131.09895828770888],[295928.3207862324,4634132.42452134,123.32110325166501]],
    [[295899.4947244857,4634151.79217599,136.82214518450198],[295916.85001777054,4634141.208161512,126.86931449392455]]
]


# ViewPoints overrides in seconds
animationEntry:
---
---
title: II Lataritio

mediaPath: /videos/l_03_lc-1080p.mp4
mediaPosition:  [295979.67624905723,4634052.450032415,129.30674310633924]
mediaRotation:  [0.353474446801508,-0.6576283060365246,-0.5859841165683883,0.3149657785973974]
mediaScale: 1
cameraFOV: 33.62

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[295976.6935711325,4634054.423029605,128.89332505755613],[295987.766762928,4634047.098277535,130.42813956366348]],
    [[295972.32113560295,4634058.029951607,128.241700141472],[295991.07368314086,4634045.625419659,130.84092198456528]],
    [[295969.1320310919,4634063.164987291,129.8283338487062],[295985.71681524586,4634047.841848732,128.26860134680385]],
    [[295963.27306370705,4634073.410534289,130.5917935798446],[295979.65996627894,4634057.896126825,128.84279633651613]],
    [[295955.9975810825,4634082.815708284,131.30252289149973],[295972.1822499686,4634067.064459871,129.80550162777445]]
]

# ViewPoints overrides in seconds
animationEntry:
---
---
title: II Lataritio

mediaPath: /videos/l_10_morpurgodoc1183-1080p.mp4
mediaPosition:  [295981.485204174,4634035.503504097,131.23335649760566]
mediaRotation:  [0.5938141622233302,-0.4256574377352556,-0.39779333362621166,0.5549422944470958]
mediaScale: 1
cameraFOV: 32

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[295978.0837036359,4634034.35,130.99],[296104.06520504673,4634077.072374005,140.00320361502355]],
    [[295957.9187101371,4634024.020414886,138.3080838156943],[295970.6923534253,4634030.700026342,134.73886430021997]],
    [[295932.3684001269,4634071.928069886,149.0009527164594],[295944.5419696395,4634065.293988815,143.67970533820517]],
    [[295938.0271846535,4634103.620598886,153.4085754516735],[295947.049181708,4634092.420185026,149.71035986656517]]
]


# ViewPoints overrides in seconds
animationEntry:
---
---
title: II Lataritio

mediaPath: /videos/l_12_mh-1080p.mp4
mediaPosition:  [295957.46904114314,4634088.831891755,128.35504726933203]
mediaRotation:  [-0.20795546868908163,0.6937913475962276,0.6604673724604887,-0.1979670148233108]
mediaScale: 1
cameraFOV: 32.68

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[295955.4912327944,4634091.834710975,128.17798551044103],[296028.7433938578,4633980.6191843,134.73582843233044]],
    [[295953.92244964116,4634094.137503885,128.0379348812164],[295967.3574567148,4634075.96150729,129.22963069318757]],
    [[295951.96749000857,4634098.122675536,127.80769925028675],[295965.4024970822,4634079.946678941,128.99939506225792]],
    [[295943.72865844925,4634102.332487901,128.28941194405226],[295964.6293318731,4634093.64739862,128.15719689101266]],
    [[295933.41528368404,4634113.759733533,129.7582910815163],[295950.0243458315,4634098.645857866,126.93053170036839]],
    [[295926.09954544366,4634120.416891969,131.00382483024347],[295942.7086075911,4634105.303016302,128.17606544909557]]
]


# ViewPoints overrides in seconds
animationEntry:
---
---
title: II Lataritio

mediaPath: /videos/l_05_doc1024-1080p.mp4
mediaPosition:  [296018.7584286861,4634004.507203135,128.81097010072668]
mediaRotation:  [0.6802352347910724,0.21057970406487023,0.20762544479924705,0.6706920964620062]
mediaScale: 1
cameraFOV: 33.62

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[296020.7922003481,4634001.5371549595,128.7601108385364],[296011.47074689745,4634015.149875762,128.99321579024175]],
    [[296030.62563977216,4633990.10334601,130.80232755670426],[296016.127638033,4634006.904321901,126.35018639969519]],
    [[296035.23218712496,4633985.532745132,132.0345026810969],[296019.32970558346,4634001.010815084,127.58176896045538]],
    [[296046.7739008096,4633974.299061175,135.3367989540107],[296030.9032731335,4633989.74612739,130.66799517521642]],
    [[296061.6528100691,4633960.150016829,137.43449672116932],[296045.11334161,4633975.224562749,134.04473941301308]],
    [[296068.6805593239,4633945.996521681,137.96333688377712],[296053.3989003267,4633962.570915734,135.9516792965021]]
]



# ViewPoints overrides in seconds
animationEntry:
---
---
title: II Lataritio

mediaPath: /videos/l_08_anonaftera-1080p.mp4
mediaPosition:  [295973.9765602476,4634060.4704073565,128.35202862906414]
mediaRotation:  [0.2881756230199118,-0.6611364258884115,-0.635012864665879,0.27678890579179055]
mediaScale: 1
cameraFOV: 32

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[295971.34142896585,4634062.918884863,128.20697337302033],[295983.2070010268,4634051.893784711,128.86013384563455]],
    [[295967.8699919203,4634066.144430238,128.01588225247306],[295984.43746111455,4634050.750481581,128.92786660121862]],
    [[295960.9292638001,4634072.761447543,127.83379860749977],[295977.29713773826,4634057.13374487,128.22625001330923]],
    [[295958.0449027033,4634080.625150065,128.2310342512706],[295972.84042212256,4634063.585372169,126.49178712934501]],
    [[295951.9076812977,4634087.9566779155,130.87117052155756],[295966.17801832355,4634070.820455389,126.99899220813285]],
    [[295946.1315245474,4634094.892849586,132.4384992302368],[295960.4018615733,4634077.75662706,128.5663209168121]]
]


# ViewPoints overrides in seconds
animationEntry:
---
---
title: II Lataritio

mediaPath: /videos/l_07_tha-1080p.mp4
mediaPosition:  [295972.7495712026,4634058.523449884,127.94878203210638]
mediaRotation:  [0.3133139071851656,-0.6359583407226416,-0.6326490093169527,0.3116835180755444]
mediaScale: 1
cameraFOV: 33.28

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[295969.8952358476,4634060.717172044,127.93],[295982.9776062247,4634050.662612143,128.0160843138209]],
    [[295966.3928454227,4634063.408962086,127.9069536507715],[295984.33848380425,4634049.616698852,128.0250391292878]],
    [[295962.9731094268,4634070.38924577,128.85682453019433],[295978.79336403764,4634054.538371249,125.57723199983398]],
    [[295960.8099248408,4634076.144473232,129.89549702745032],[295976.5809440181,4634060.063357204,127.6686006182648]],
    [[295942.22432647494,4634097.922161409,133.43705489057925],[295957.629847835,4634081.538572712,130.88081952715288]]
]

# ViewPoints overrides in seconds
animationEntry:
---
title: II Lataritio

mediaPath: /videos/l_13_en-1080p.mp4
mediaPosition:  [296007.2195179315,4634012.575564231,130.41950586925992]
mediaRotation:  [0.6883263438769712,-0.025198574066360124,-0.026522033031623468,0.7244780589871946]
mediaScale: 1
cameraFOV: 60

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[296006.9566342095,4634008.989899315,130.6036233778278],[296016.69306835823,4634141.792303621,123.78445639383216]],
    [[296006.0461035785,4633996.57050067,137.49863238330684],[296005.4291619688,4634010.753634286,133.1419588967056]],
    [[296034.63255971135,4633990.682315243,141.69768556775898],[296023.12304670183,4633997.600636831,135.35809278084878]],
    [[296051.72111732716,4633980.410452572,160.66775029824632],[296040.21160431765,4633987.328774159,154.32815751133612]]

]


# ViewPoints overrides in seconds
animationEntry:
---
---
title: II Lataritio

mediaPath: /videos/l_06_da-1080p.mp4
mediaPosition:  [295945.37963804876,4634094.750070662,128.30463065297377]
mediaRotation:  [0.2144932720551034,-0.6804183388348559,-0.6683105260287648,0.21067643726688132]
mediaScale: 1
cameraFOV: 32

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[295943.3154252435,4634097.698773491,128.24],[295954.9956417342,4634081.01372662,128.60570842732068]],
    [[295939.2243212669,4634103.542865649,128.11190723713457],[295952.20233959076,4634085.003924681,128.51824993417173]],
    [[295932.0798828267,4634113.565387277,129.74377608374508],[295945.1776379302,4634095.106691979,129.6352157781097]],
    [[295922.40201496606,4634126.294814495,129.82047755835427],[295936.1534809703,4634108.317834672,129.7119172527189]],
    [[295912.9755079834,4634136.212008635,129.88600094016064],[295927.577342914,4634118.918627313,129.77744063452528]]
]

# ViewPoints overrides in seconds
animationEntry:
---
---
cameraFOV: 60

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
  [[295951.2189581529,4633997.797692887,169.5321795877354],[296038.0927089406,4634075.278521221,104.51333276250504]],
  [[295901.36590624216,4633952.561407406,208.1466148864326],[295986.73964796076,4634030.660916376,141.89279460835658]]
]
---
---
title: II Lataritio

mediaPath: /videos/l_04_repropb-1080p.mp4
mediaPosition:  [295973.51953941095,4634065.265634941,129.1063337949326]
mediaRotation:  [-0.17966862357752497,0.7051976549584036,0.6646362084920034,-0.16933447228576914]
mediaScale: 1
cameraFOV: 33.62

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[295971.7999738832,4634068.42122073,128.89332505756167],[295978.18386090494,4634056.706108487,129.68411999505125]],
    [[295966.71502895455,4634073.931527652,128.4456244234578],[295980.9233208697,4634056.330518486,129.23004994067708]],
    [[295958.8575066606,4634081.090122916,129.25729254312415],[295973.28892724944,4634063.665793975,128.61443413057728]],
    [[295950.8653859389,4634090.739716284,129.61330753197473],[295965.2968065277,4634073.315387343,128.97044911942785]],
    [[295943.5366295166,4634099.588371332,129.93977246247974],[295957.9680501054,4634082.164042391,129.29691404993287]]
]

# ViewPoints overrides in seconds
animationEntry:
---
---
title: II Lataritio

mediaPath: /videos/l_01_gbp-1080p.mp4
mediaPosition:  [295988.4288916751,4634048.552979653,129.12995054961874]
mediaRotation:  [0.4955672507227679,-0.533420474969569,-0.5021949207591757,0.46655756178418306]
mediaScale: 1
cameraFOV: 57.35

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[295984.8451441457,4634048.817006568,128.91305556725783],[296001.2706536554,4634047.606883205,129.90715756974538]],
    [[295980.18676220835,4634052.8898631,128.59213978546688],[295996.2909194237,4634049.387280499,129.39091103420478]],
    [[295977.0078573169,4634057.889158933,130.0803407226552],[295990.24598302686,4634048.739422631,126.435570045962]],
    [[295970.98743230046,4634065.130705993,132.15312180644716],[295984.4661468205,4634056.415207314,128.33025146096168]],
    [[295967.71740319533,4634068.813520942,131.13522741662942],[295980.5398971743,4634058.962108401,127.85143097053324]]
]


# ViewPoints overrides in seconds
animationEntry:
---
---
title: II Lataritio

mediaPath: /videos/l_11_gb-1080p.mp4
mediaPosition:  [295947.45706993673,4634101.795880663,128.72348129287627]
mediaRotation:  [0.22663136741219683,-0.6794275062170797,-0.6620080237841688,0.22082088566505073]
mediaScale: 1
cameraFOV: 35

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[295945.2966122519,4634104.674019617,128.63],[295955.198709974,4634091.48254941,129.0584559256829]],
    [[295941.7629300993,4634108.349499583,129.09185173113482],[295955.1257791956,4634090.112566778,128.02651042763307]],
    [[295931.5269565254,4634116.372172767,129.28464548663663],[295949.38965473534,4634102.510446238,128.2536476649642]],
    [[295916.78284950904,4634129.517253618,130.1866311142116],[295933.71758921444,4634114.535997301,129.15563329253916]],
    [[295900.3778876949,4634154.1622176105,133.65064133606938],[295915.63707211526,4634137.6513495315,131.0362709877609]],
    [[295885.43938864185,4634168.31525069,136.78222645210545],[295901.84448059194,4634153.071565366,133.49750088188998]]

]




# ViewPoints overrides in seconds
animationEntry:
---
---
title: II Lataritio

mediaPath: /videos/l_02_cl-1080p.mp4
mediaPosition:  [295974.01251872926,4634055.690910948,129.0347577030788]
mediaRotation:  [0.3437805315965897,-0.6484006519130134,-0.6001249346090145,0.31818485751913383]
mediaScale: 1
cameraFOV: 34.46

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[295971.04163645505,4634057.705007552,128.75677679607523],[295976.31696449587,4634054.128622051,129.25038116614354]],
    [[295966.7748483915,4634062.034969717,129.15035768954476],[295985.16398109135,4634048.888736481,128.0089741207865]],
    [[295964.1865563574,4634068.069218763,129.46756604286844],[295981.7545450052,4634054.05608584,126.76862426019791]],
    [[295962.45772999,4634071.635438252,130.03536370800438],[295978.6719053416,4634056.442164829,125.72823397196467]],
    [[295954.55906454753,4634083.967917735,133.1150151404958],[295969.98353529646,4634068.107939822,128.3366044310981]],
    [[295943.0918609484,4634100.402849886,137.7802346741112],[295957.36485518375,4634084.102294504,131.23348463950245]]
]

# ViewPoints overrides in seconds
animationEntry:
---
---
title: VII Horatio

mediaPath: /videos/o_06_fllia-doc1205-1080p.mp4
mediaPosition:  [296281.3438851955,4633644.320760986,135.49342475781827]
mediaRotation:  [0.3174651248933075,-0.5560583403324293,-0.6670646348758321,0.3808408979103919]
mediaScale: 1
cameraFOV: 38

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[296278.29440068576,4633646.12092767,136.14153287750545],[296292.2712046889,4633637.870163702,133.17103732893918]],
    [[296258.14037813555,4633648.338218782,144.5677378103193],[296273.2841281333,4633647.053921963,138.1437656991585]],
    [[296270.65906772384,4633672.832560453,132.1328162421325],[296269.66325812205,4633656.417483199,133.4757987481838]],
    [[296264.84755485796,4633697.005328261,129.7608352914263],[296272.7697992155,4633682.547167,130.43137036686943]]
]


# ViewPoints overrides in seconds
animationEntry:
---
---
title: VII Horatio

mediaPath: /videos/o_13_wwm1931-1080p.mp4
mediaPosition:  [296269.42349810235,4633679.871358391,130.1693037420755]
mediaRotation:  [0.24750856045777517,-0.7909875628530889,-0.5340014141950473,0.1670948160462801]
mediaScale: 1
cameraFOV: 60
# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[296267.52025123546,4633682.614782166,128.82346669504082],[296276.2434660419,4633670.040756533,134.99188649394975]],
    [[296262.3023348634,4633687.168037223,128.57091769940814],[296278.25943624065,4633668.771661901,132.7394869746622]],
    [[296256.8167449558,4633694.097717676,129.22901712821138],[296271.82462282793,4633674.507589396,130.42197769140535]],
    [[296249.02193833486,4633710.870328794,128.98767153755634],[296261.9554544392,4633689.821540036,129.31222072037642]],
    [[296241.6678447517,4633722.047128253,130.00064957193146],[296255.209632592,4633701.4662306495,128.1353573406845]]
]


# ViewPoints overrides in seconds
animationEntry:
---
---
title: VII Horatio

mediaPath: /videos/Tumuli_NL_final.m4v
mediaPosition:  [296267.6764461834,4633681.292140067,129.4651425209769]
mediaRotation:  [0.042262894311928575,-0.7566631997791755,-0.6514221376413435,0.03638472832248394]
mediaScale: 1
cameraFOV: 40.12

playDT: 306  # 5.06 min

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[296267.28,4633684.83,128.93],[296269.52326544083,4633664.811245235,131.9580698203426]],
    [[296259.5259870479,4633698.672024895,130.28253846183435],[296269.6875584719,4633681.186404892,127.84377887524437]],
    [[296246.0600045914,4633736.76597715,147.23277992603101],[296251.01585910906,4633719.301816105,137.99186146337055]],
    [[296391.6478276091,4633758.370369164,217.82906967374362],[296360.4952625182,4633734.343272505,191.53184056680539]]
]
---
---
title: VII Horatio

mediaPath: /videos/o_08_stereoph1900-1080p.mp4
mediaPosition:  [296292.80817627144,4633660.533849068,130.00576565358114]
mediaRotation:  [0.6407602271715748,0.31802356891735517,0.3106590136824367,0.6259219744552478]
mediaScale: 1
cameraFOV: 64

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
[[296295.67461061006,4633658.357512136,129.92143450797545],[296282.53678655816,4633668.332389739,130.30795225866814]],
[[296302.4167907011,4633653.201294178,129.72255058771447],[296282.7443436664,4633668.1375766685,130.30131842236983]],
[[296302.93162869534,4633642.981207559,129.6893504711867],[296283.5692887321,4633657.229420387,123.98664045703154]],
[[296307.53123708995,4633640.561576929,133.89196454996554],[296288.7713091622,4633654.397427872,125.70251939163052]]
]


# ViewPoints overrides in seconds
animationEntry:
---
---
title: VII Horatio

mediaPath: /videos/o_10_anonbc-1080p.mp4
mediaPosition:  [296315.86224985344,4633628.70328375,129.6823531752794]
mediaRotation:  [0.6562879868903723,0.29063670872921493,0.2819409036612436,0.6366519525181032]
mediaScale: 1
cameraFOV: 32.95

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[296318.52674561786,4633626.284914256,129.5730311242569],[296306.3144733642,4633637.369107771,130.07409052477664]],
    [[296324.3296090237,4633620.790130591,129.93577931057712],[296305.23473813717,4633636.408942001,128.56781910550694]],
    [[296333.72710887034,4633614.796247588,131.2915550520482],[296313.1690005006,4633627.799989232,126.96725595747267]],
    [[296340.8536280812,4633605.203554902,132.5935176628926],[296321.1223574813,4633619.536473097,128.63408642929463]],
    [[296340.5940716387,4633598.812506522,132.2713779626357],[296322.8945239197,4633615.48860476,127.9048647930136]]
]


# ViewPoints overrides in seconds
animationEntry:
---
---
title: VII Horatio

mediaPath: /videos/o_11_cf-1080p.mp4
mediaPosition:  [296327.1235663774,4633608.811843919,129.6133526881316]
mediaRotation:  [-0.6590613333033295,-0.23715254936750252,-0.24165219217643388,-0.6715661138631518]
mediaScale: 1
cameraFOV: 38

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[296329.4169624474,4633606.0377211785,129.6810098444967],[296318.90556379326,4633618.752450403,129.37091454448998]],
    [[296336.7117719505,4633597.213821466,129.89621299831006],[296320.97210817767,4633616.252730515,129.4318793810409]],
    [[296346.7853243179,4633586.656869722,129.90906264067033],[296329.64244991605,4633604.448844315,129.92234711979125]],
    [[296352.6965035537,4633581.380918222,131.07007995273642],[296336.34581352293,4633599.623796364,127.86305457937794]]

]


# ViewPoints overrides in seconds
animationEntry:
---
---
title: VII Horatio

mediaPath: /videos/o_12_kh1925-1080p.mp4
mediaPosition:  [296271.040463389,4633667.625302507,131.06975683074054]
mediaRotation:  [-0.28820051604787067,0.6870330691442924,0.6151005467878766,-0.25802585489281926]
mediaScale: 1
cameraFOV: 38

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[296268.48774234095,4633670.132567377,130.67322420679434],[296280.1877138113,4633658.640936722,132.49066539988115]],
    [[296264.13433201046,4633678.533402552,131.0613010987852],[296278.73682032106,4633658.628278729,130.06829562672178]],
    [[296263.1284797532,4633686.67848325,130.7526287744177],[296276.812615471,4633666.222903524,128.57401295513702]],
    [[296257.242682809,4633699.160336925,131.88397082859098],[296270.18824151997,4633678.192266343,130.09925658541664]],
    [[296245.32207750273,4633717.2954250835,132.65491915883518],[296260.6697355129,4633697.944739289,131.99686745836425]]
]


# ViewPoints overrides in seconds
animationEntry:
---
---
title: VII Horatio

mediaPath: /videos/o_09_anon1905-1080p.mp4
mediaPosition:  [296328.4071638867,4633603.846947455,129.51731003503295]
mediaRotation:  [-0.6259149506804219,-0.16506657272623707,-0.19436831597560006,-0.7370240558000453]
mediaScale: 1
cameraFOV: 60

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[296330.1590395894,4633600.756486676,130.10039124283156],[296322.1296092855,4633614.921098581,127.42793570708798]],
    [[296335.4743813869,4633594.713636665,131.02945360161755],[296319.6203996532,4633612.222950145,123.78322121526078]],
    [[296341.28370725946,4633592.712839693,132.02001553083795],[296322.8768240659,4633607.486093988,124.71451843597308]],
    [[296346.2590418207,4633588.350019979,134.06823339987147],[296328.5132283828,4633603.911141301,126.7627363050066]]
]


# ViewPoints overrides in seconds
animationEntry:
---
---
title: VII Horatio

mediaPath: /videos/o_01_jhp-1080p.mp4
mediaPosition:  [296278.9448397351,4633685.012502358,129.44849388063759]
mediaRotation:  [0.15790910494069182,0.7171345215450977,0.6629287957412491,0.14597330016880003]
mediaScale: 1
cameraFOV: 42.3

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
[[296280.4522676311,4633688.269484362,129.1661299872504],[296273.54322310793,4633673.34165018,130.46029783194166]],
[[296271.30300208076,4633694.983466773,129.0091050014435],[296282.70014707115,4633683.054522674,128.7737781228023]],
[[296264.6142570416,4633695.8550837785,129.5484930388437],[296277.76742907707,4633685.912275758,128.92830430785907]],
[[296257.72575526783,4633700.779236812,129.12232092773263],[296270.8789273033,4633690.836428791,128.502132196748]]
]

# ViewPoints overrides in seconds
animationEntry:
---
---
title: VII Horatio

mediaPath: /videos/o_14_ld1957-1080p.mp4
mediaPosition:  [296265.3515622494,4633678.235480788,128.9749748310552]
mediaRotation:  [0.2952799363457504,-0.5660787510530326,-0.6823944116363039,0.3559529094774068]
mediaScale: 1
cameraFOV: 38

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[296262.45,4633680.26,129.64],[296275.74882697646,4633670.980953614,126.59196797566968]],
    [[296257.9735056914,4633685.786752902,130.1398672781183],[296277.8870369066,4633671.89240116,125.57577444712246]],
    [[296255.02562171506,4633690.31993821,130.07835481746898],[296272.92326609883,4633673.798383573,125.93767406954292]],
    [[296257.61986240785,4633698.458838994,131.11116690177337],[296267.9833223006,4633676.587036213,126.14572418829702]],
    [[296257.88932713185,4633706.264850471,137.74334939800153],[296268.20418151544,4633685.375562261,129.5169226954936]]

]


# ViewPoints overrides in seconds
animationEntry:
---
---
title: VII Horatio

mediaPath: /videos/O_03_PPM-1856-1080p.mp4
mediaPosition:  [296299.6866396207,4633633.6146786995,130.17635778273052]
mediaRotation:  [0.5796498791218397,0.42981574855569576,0.41234554762762754,0.556089551492577]
mediaScale: 1
cameraFOV: 69.46

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
[[296303.12847069517,4633632.569925356,130.02706175020884],[296287.35341160384,4633637.358378181,130.7113352325999]],
[[296323.42033147515,4633628.039965847,135.23674949007096],[296308.947751032,4633633.933243305,129.93923157002888]],
[[296318.1515853004,4633628.009722124,131.8449053039193],[296302.37652621145,4633632.798174958,132.52917878621597]],
[[296310.30374983023,4633629.850592081,131.15812854448626],[296294.33203538344,4633633.7338506505,129.7177672985143]]
]


# ViewPoints overrides in seconds
animationEntry:
---
---
title: VII Horatio

mediaPath: /videos/o_16_en8962-1080p.mp4
mediaPosition:  [296281.06573691807,4633646.365376869,136.50237163979122]
mediaRotation:  [-0.3342670528553986,0.47459220086404563,0.6657166365684768,-0.46888073115702045]
mediaScale: 1
cameraFOV: 33.55

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[296277.8613461288,4633647.511707183,137.6761716576228],[296292.5481372462,4633642.257693246,132.29625490922803]],
    [[296267.44607705175,4633651.670539702,136.79517651798182],[296282.83210607286,4633646.140288535,134.57336802504184]],
    [[296265.2601306981,4633654.839318328,133.6913365252595],[296279.4135753853,4633646.358237991,133.65556211605124]],
    [[296264.0058015141,4633663.019399114,132.95572940701572],[296273.2836575967,4633649.390903076,132.29581860867953]],
    [[296273.3496952584,4633681.632409446,130.95351421377381],[296275.42583555117,4633665.346385411,129.3088067583222]],
    [[296269.82496467384,4633691.337449143,130.0267903483712],[296277.38765323913,4633676.672873341,129.94960097381028]]
    
]

# ViewPoints overrides in seconds
animationEntry:
---
---
title: VII Horatio

mediaPath: /videos/o_05_fllia-doc1243-720p.mp4
mediaPosition:  [296320.30155777535,4633621.396739068,129.9413161361397]
mediaRotation:  [0.6555889766913671,0.22386753978154497,0.23304906370706432,0.6824767777574453]
mediaScale: 1
cameraFOV: 38

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[296322.50165309495,4633618.550923719,130.08593834770477],[296312.4178828801,4633631.5942440685,129.42308654469826]],
    [[296329.3851018951,4633613.136496061,130.29337200192836],[296316.35293686204,4633623.213315033,129.36042294034445]],
    [[296340.8945858969,4633604.666607679,130.99515391259342],[296326.264517631,4633612.210761569,129.8563630577775]],
    [[296355.6567980208,4633596.814889012,139.16150124893002],[296341.74783012335,4633603.83858324,133.7339482627857]]
]


# ViewPoints overrides in seconds
animationEntry:
---
---
cameraFOV: 60

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
  [[296325.76537415583,4633699.399426047,177.22194053436957],[296246.8796625461,4633633.276422189,92.47189249654215]],
  [[296355.32871531876,4633725.955291409,215.04293271724057],[296283.3627133376,4633662.6203386905,122.37494417350551]]
]
---
---
title: VII Horatio

mediaPath: /videos/o_07_av-1080p.mp4
mediaPosition:  [296314.7998573784,4633646.647266331,131.21341609172242]
mediaRotation:  [0.5627296968174571,0.41672620963112805,0.424869030889573,0.5737254231051304]
mediaScale: 1
cameraFOV: 37.78

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
[[296318.24270183995,4633645.597518769,131.2830728843111],[296302.46299805795,4633650.408861764,130.9638125849463]],
[[296324.5381053638,4633640.516764739,131.43012327785254],[296302.1029563091,4633650.827558057,130.54451208002268]],
[[296328.839183178,4633633.169648084,132.20281999362652],[296310.87908934744,4633649.910169299,129.44190380471747]],
[[296325.24801417347,4633624.669441907,131.64834734368048],[296313.85139066755,4633646.402928262,128.78387136124076]],
[[296324.32691261254,4633614.870199172,132.66862130087017],[296314.0321609845,4633637.045127768,129.1016464171916]]

]


# ViewPoints overrides in seconds
animationEntry:
---
---
title: VII Horatio

mediaPath: /videos/o_15_en710-1080p.mp4
mediaPosition:  [296314.9457049262,4633629.465386899,129.91177725275656]
mediaRotation:  [0.6601168672060141,0.2966568287278902,0.28287952251259624,0.6294597869143621]
mediaScale: 1
cameraFOV: 38

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[296317.6346679627,4633627.077874959,129.74070847547102],[296305.31025404565,4633638.020638016,130.52477370469649]],
    [[296326.581529129,4633618.040665917,131.00391569746031],[296314.5372482453,4633629.299243531,131.65945219364082]],
    [[296337.63283553894,4633603.141297237,131.5858641396188],[296326.2366270895,4633615.056321192,132.22511145384905]],
    [[296349.06780231773,4633591.139328745,132.40155306006994],[296337.7526146539,4633603.1404371075,132.83735830823247]],
    [[296364.2783148446,4633580.916955519,132.31511860539433],[296351.98793895676,4633591.919717688,131.94649252825286]]
]

# ViewPoints overrides in seconds
animationEntry:
---
---
title: VII Horatio

mediaPath: /videos/o_02_tha-1080p.mp4
mediaPosition:  [296304.8048670289,4633640.945327751,128.94492268554686]
mediaRotation:  [-0.6738408353456847,-0.27785121509842603,-0.2609890788716513,-0.6329470211626602]
mediaScale: 1
cameraFOV: 38.68

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
[[296307.3373244534,4633638.39660395,128.71983074794736],[296295.7302279245,4633650.078254704,129.75150212861175]],
[[296314.56229011144,4633631.18770575,130.0428827952662],[296302.9551935827,4633642.869356504,131.07455417591277]],
[[296320.84050314437,4633628.892547226,131.57404123397845],[296306.2692159243,4633636.539033352,130.36621729722248]],
[[296330.854276373,4633632.533385606,132.37307133287746],[296314.56556013523,4633633.982922357,130.17612651486849]]
]

# ViewPoints overrides in seconds
animationEntry:
---
---
title: VII Horatio

mediaPath: /videos/O_04_PPM-0131-1080p.mp4
mediaPosition:  [296306.44280968956,4633639.799743179,129.42487967355106]
mediaRotation:  [0.608257107039782,0.4525048656295728,0.3892397213454918,0.5232160907732496]
mediaScale: 1
cameraFOV: 36.57

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
[[296309.85212239623,4633638.776500506,128.88677066943458],[296294.2261058241,4633643.466362758,131.3531036049684]],
[[296317.16250540066,4633636.682352375,130.43624303524726],[296301.40807580965,4633641.435029823,129.22793622978466]],
[[296316.3059144765,4633627.544224055,130.6952234605514],[296306.43448221666,4633640.331988238,127.33696314678429]],
[[296322.4462672035,4633619.851459048,132.802043578451],[296311.22091946006,4633631.776108826,130.7910001362211]]
]

# ViewPoints overrides in seconds
animationEntry:
---
---
title: V Piramide

mediaPath: /videos/p_09_tha1313-1080p.mp4
mediaPosition:  [296129.86820310395,4633857.526811792,129.11260190597346]
mediaRotation:  [0.33983239311468716,-0.6736811109499271,-0.5859257897265207,0.295565008538849]
mediaScale: 1
cameraFOV: 38

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[296127.00092459226,4633859.645660715,128.61340938441583],[296130.71252690523,4633856.902877084,129.25959848327474]],
    [[296110.942786915,4633880.407340599,132.2999004177942],[296179.52874724293,4633810.449650386,118.49655375709992]],
    [[296085.6367242989,4633905.712243116,169.64246003064073],[296153.1571246125,4633842.760739255,134.05105016756013]]


    ]

# ViewPoints overrides in seconds
animationEntry:

---
---
title: V Piramide

mediaPath: /videos/p_05_jhp-1080p.mp4
mediaPosition:  [296272.63161998463,4633707.11997207,130.88051880923908]
mediaRotation:  [0.734488980519281,0.15809622664875408,0.13887096476501587,0.6451715862934502]
mediaScale: 1
cameraFOV: 37.8

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[296274.1004083686,4633703.866173665,130.41634575766832],[296267.3684616086,4633718.779416355,132.5438055773676]],
    [[296290.5863403461,4633686.358645248,131.81310893896355],[296227.4511858436,4633762.488275892,129.19082805690005]],
    [[296309.5149322304,4633689.445581479,147.33423038810864],[296222.74197712244,4633725.192786191,116.00953772902247]],
    [[296291.3962075654,4633650.329868766,158.09097378605674],[296252.6560938854,4633722.070896988,102.04569553081083]]
 
   
]
# ViewPoints overrides in seconds
animationEntry:
---
---
title: Piramide

mediaPath: /videos/Piramide_NL_final.m4v
mediaPosition:  [296255.264388564,4633689.545982254,129.91779272168478]
mediaRotation:  [-0.7646388723342585,-0.06252773282003134,-0.05227698700136337,-0.6392845956018647]
mediaScale: 1
cameraFOV: 43.1

playDT: 289  # 4.49 min

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[296255.84,4633686.05,129.28],[296253.2017809183,4633702.073251996,132.20321664105526]],
    [[296264.01584076684,4633680.246540324,129.37794680269718],[296254.9838066788,4633694.047281402,128.917614399112]],
    [[296280.9222167666,4633664.234170015,130.7784677493043],[296272.0517804476,4633678.095488321,129.58297022893225]],
    [[296284.0482965831,4633620.980983169,154.08785247249625],[296282.98450312024,4633636.730884951,149.28588023806142]],
    [[296172.7522053552,4633670.847710761,159.06321479976793],[296238.1949463322,4633726.515405634,135.67181358310688]]
]




# ViewPoints overrides in seconds
animationEntry:
---


#title 

<html>
---
title: V Piramide

mediaPath: /videos/p_16_cinec1955-1080p.mp4
mediaPosition:  [296327.24094382575,4633608.135185997,130.09782779923]
mediaRotation:  [0.6842438506693153,0.10685432700846685,0.11130499378310539,0.7127437856479273]
mediaScale: 1
cameraFOV: 37.36

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[296328.3376443344,4633604.709446645,130.24465383969974],[296327.65863277507,4633606.830460178,130.15374791058957]],
    [[296331.83914827625,4633597.218290736,130.601775630329],[296325.6932296471,4633609.07480818,130.0784179251639]],
    [[296340.1697282867,4633596.62262854,131.1052091117439],[296330.08016772167,4633605.303951589,129.89747457325373]],
    [[296342.20008419105,4633594.875656964,133.10262455848886],[296332.11052362603,4633603.556980014,131.89489001999868]]
    ]

# ViewPoints overrides in seconds
animationEntry:
---
---
title: V Piramide

mediaPath: /videos/p_07_da1890-1080p.mp4
mediaPosition:  [296389.483768546,4633554.8862472875,132.0313962924016]
mediaRotation:  [0.7163559011333975,0.10284848970174086,0.09807514889416333,0.6831088319197876]
mediaScale: 1
cameraFOV: 37.79

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[296390.495465194,4633551.435559722,131.8604424521591],[296389.18585554033,4633555.902366799,132.08173684901442]],
    [[296395.0951009992,4633552.784113026,133.05009245225105],[296380.3674781091,4633559.759539055,130.46348159470466]],
    [[296453.5065608116,4633529.2711146055,157.42921862879516],[296364.33311837015,4633560.881518374,128.48920671202546]],
    [[296399.3280654969,4633555.09817578,137.02853777549015],[296384.2963170496,4633560.532627328,132.93423991372936]],
    [[296382.50420694356,4633540.223987577,138.24109540555455],[296374.6808523131,4633553.994564593,133.61331721714018]],
    [[296346.55422760383,4633586.490866935,145.00558959989354],[296357.7211675596,4633575.054067764,140.91302600768287]]

    ]

# ViewPoints overrides in seconds
animationEntry:
---
---
title: V Piramide

mediaPath: /videos/p_18_en8961-1080p.mp4
mediaPosition:  [296267.68504012364,4633654.44555816,135.1295542051014]
mediaRotation:  [0.7349615609161011,0.09397642017011176,0.08517672554544556,0.6661417730882327]
mediaScale: 1
cameraFOV: 37.8

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[296268.5865034396,4633650.978153406,134.77675374738143],[296264.454796575,4633666.8704251945,136.39375584526474]],
    [[296279.17338599195,4633645.568643073,135.8123350758666],[296268.0528632074,4633657.644383421,134.15086944622473]],
    [[296293.9244549263,4633621.697087516,146.32365672299454],[296248.3195816809,4633708.111120746,130.78482953640213]],
    [[296314.8885013864,4633633.577392351,132.95593603306256],[296301.54619511944,4633643.255194937,132.19896800938778]]
]


# ViewPoints overrides in seconds
animationEntry:
---
---
title: V Piramide

mediaPath: /videos/p_10_tha1314-1080p.mp4
mediaPosition:  [296412.6487352082,4633500.714295922,132.22797250200384]
mediaRotation:  [0.7062976093082939,0.16812629423762127,0.15924065630082582,0.6689690946913734]
mediaScale: 1
cameraFOV: 37.5

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[296414.26832185406,4633497.50512106,132.03268860315728],[296412.1718176197,4633501.6592975,132.2854774995605]],
    [[296417.3883198871,4633491.322914454,132.41027652696835],[296411.37560446403,4633503.236976129,133.1352680014782]],
    [[296435.042038029,4633477.918515747,133.34182843525247],[296425.0749369685,4633486.77964893,132.46979622913406]],
    [[296465.17370685586,4633459.362694111,160.7743180009233],[296397.68768460985,4633527.068175406,135.27433060754205]]

    ]

# ViewPoints overrides in seconds
animationEntry:
---
---
title: V Piramide

mediaPath: /videos/p_15_am1950s-1080p.mp4
mediaPosition:  [296219.2543402447,4633746.389190704,129.3620121278183]
mediaRotation:  [0.32657019415962507,-0.7085747092771706,-0.5680872270890002,0.2618218708924996]
mediaScale: 1
cameraFOV: 60

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[296216.58284711756,4633748.671799666,128.5791834133654],[296218.23687636235,4633747.258544307,129.06386427374818]],
    [[296211.8488435159,4633753.594280476,130.65850737440846],[296222.01207913575,4633744.932203686,130.11051794408203]],
    [[296206.9104416072,4633761.813033074,130.71057226668984],[296215.57145475945,4633751.634579439,130.806268504914]],
    [[296208.3177200565,4633764.592384644,134.28582259160427],[296215.7773586605,4633754.495888771,129.69912574424808]]
    ]

# ViewPoints overrides in seconds
animationEntry:
---
---
title: V Piramide

mediaPath: /videos/p_11_anonls-1080p.mp4
mediaPosition:  [296226.4356814254,4633727.015552151,129.09895036517995]
mediaRotation:  [0.4786144806511472,-0.5436175952525063,-0.5175067201453011,0.45562581538732844]
mediaScale: 1
cameraFOV: 35


# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[296222.869,4633727.471,128.92188860628895], [296223.5619805521,4633727.382509758,128.95629043274812]],
    [[296218.8438690378,4633729.532915271,129.65371289466907],[296231.7952655059,4633726.358892655,128.75363791862674]],
    [[296166.8432677564,4633719.736192578,153.14314502155102],[296261.43383467384,4633733.542320218,127.63510231195939]]

    ]

# ViewPoints overrides in seconds
animationEntry:
---
---
title: V Piramide

mediaPath: /videos/p_01_gbp-1080p.mp4
mediaPosition:  [296260.3969488976,4633755.224059511,130.41739586745825]
mediaRotation:  [-0.09814850366851359,-0.7699304470663743,-0.6254730592184153,-0.07973349421517656]
mediaScale: 1
cameraFOV: 38

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[296261.2809540234,4633758.63502372,129.67992850839445],[296259.3430431802,4633751.157528021,131.2966001203134]],
    [[296265.5462703157,4633776.157408916,134.0583794793378],[296234.3944266787,4633646.514517946,133.81091898001787]],
    [[296286.75633808767,4633774.645358458,137.02080523405093],[296271.83619517897,4633757.731480442,135.12436851272284]],
    [[296304.9724547897,4633739.859775183,135.79971395123428],[296282.46901526186,4633740.219369482,133.4013170776563]],
    [[296259.0717495834,4633700.535316206,136.49979086495173],[296259.58274742286,4633723.139622957,135.46509107110163]],
    [[296263.01091832394,4633692.773882436,130.69833351647264],[296248.8906871611,4633710.460201531,131.013920229808]]
]

# ViewPoints overrides in seconds
animationEntry:
---
---
title: V Piramide

mediaPath: /videos/p_08_da1920-1080p.mp4
mediaPosition:  [296265.9736016267,4633688.041618925,128.45771918019872]
mediaRotation:  [0.7049477458791472,0.13284741898340935,0.1290233529754182,0.6846555434934682]
mediaScale: 1
cameraFOV: 34.46

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[296267.2833496213,4633684.689975994,128.35260089991408],[296265.5879218258,4633689.028572892,128.4886732244729]],
    [[296269.83262856473,4633678.166375434,128.84539984724069],[296263.8272471603,4633693.534136421,128.7138079625158]],
    [[296276.9484124847,4633674.310228343,129.08367759466827],[296266.1980738428,4633686.809562426,128.41519679826285]],
    [[296292.3206313309,4633655.369584191,142.04421675217975],[296282.33286741906,4633667.882500148,138.0540405318861]]
    ]

# ViewPoints overrides in seconds
animationEntry:
---
---
title: V Piramide

mediaPath: /videos/p_17_ppp-1080p.mp4
mediaPosition:  [296249.9980481078,4633700.032650234,129.34376131982805]
mediaRotation:  [0.7373411047215563,0.02950938420611745,0.026987832303948637,0.6743359314461187]
mediaScale: 1
cameraFOV: 38

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[296250.2845971362,4633696.458431949,129.0230538188568],[296250.10718306457,4633698.671374783,129.2216167798624]],
    [[296250.59968798823,4633692.52820231,128.67040212686243],[296249.6422560466,4633704.4705591565,129.74196606470767]],
    [[296254.82622632117,4633691.071958861,130.4190016634729],[296248.4433313988,4633701.260653038,130.05297483042438]],
    [[296266.5243191027,4633684.647038104,131.05134375879487],[296258.83256888133,4633693.894767671,131.01544605373317]]
    ]

# ViewPoints overrides in seconds
animationEntry:
---
---
title: V Piramide

mediaPath: /videos/p_04_doc1294-1080p.mp4
mediaPosition:  [296410.6149488505,4633507.258520561,132.1715514454357]
mediaRotation:  [0.6807695560508688,0.21483495526833904,0.21074844738893353,0.66782021942459]
mediaScale: 1
cameraFOV: 34

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[296412.6809330788,4633504.311160666,132.1024224709896],[296410.73133051704,4633507.092488967,132.167657250532]],
    [[296416.34105369105,4633499.089585124,132.610252811966],[296406.8719593148,4633512.59831798,132.9270939447205]],
    [[296421.7747064126,4633493.345333479,133.91748345675958],[296411.76193202863,4633506.339014809,132.1402399026535]],
    [[296432.371391311,4633481.323498957,135.33617127441607],[296420.91221582703,4633493.157692673,134.3932204126244]],
    [[296326.1997821649,4633634.827896374,142.96622360499416],[296386.60386071,4633558.086310122,127.13205208589153]]
]

# ViewPoints overrides in seconds
animationEntry:
---
---
title: V Piramide

mediaPath: /videos/p_19_en5859-1080p.mp4
mediaPosition:  [296228.1019324927,4633717.206641001,130.11219001507243]
mediaRotation:  [0.679075263382079,-0.34847180523646626,-0.2949716148520356,0.5748181747730154]
mediaScale: 1
cameraFOV: 37

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[296225.2174983432,4633715.136241312,129.5176441659206],[296228.5779271175,4633717.548302162,130.2103030638653]],
    [[296222.27156966156,4633713.021701851,129.63147292998656],[296235.4713210865,4633722.880081682,130.54191756316706]],
    [[296188.22234214994,4633686.809506215,138.89025198759708],[296216.7073255064,4633711.901999345,133.58030133358787]],
    [[296199.5428901988,4633608.971742024,164.6795435968508],[296214.1934869745,4633643.07502513,155.11141108700437]]
    ]

# ViewPoints overrides in seconds
animationEntry:
---
---
cameraFOV: 60

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
[[296180.60198688484,4633687.925886628,183.84771359496858],[296263.4499923041,4633753.861060878,102.81353208845232]],
[[296156.1421396154,4633669.513035882,218.1342591472983],[296228.17551974684,4633723.337095716,119.68660550515644]]
]

# ViewPoints overrides in seconds
animationEntry:

# Exlude from page view points menu
exclude: true
---

#title 

<html>
---
title: V Piramide

mediaPath: /videos/p_20_js-1080p.mp4
mediaPosition:  [296229.8254746211,4633718.337931467,129.6028928524787]
mediaRotation:  [0.6707996132284927,-0.38407425428221714,-0.3152409214888981,0.5505797013229725]
mediaScale: 1
cameraFOV: 37.5

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[296226.7804003909,4633716.550511819,128.90099974432715],[296228.1548047853,4633717.35726963,129.21780153611297]],
    [[296223.2410091114,4633714.472934381,130.32366598840733],[296237.1051352677,4633723.346178534,129.1828136248722]],
    [[296164.8909522803,4633695.788251793,140.47637383704418],[296200.6890300446,4633708.896327622,136.48901372116737]],
    [[296161.91783343063,4633710.335302223,145.1807972006394],[296199.2860675021,4633715.647512562,138.5010840359585]]
    ]

# ViewPoints overrides in seconds
animationEntry:
---
---
title: V Piramide

mediaPath: /videos/p_03_lc-1080p.mp4
mediaPosition:  [296225.9229882462,4633739.996877952,129.51940815510264]
mediaRotation:  [0.367179183978474,-0.7125826102813708,-0.5314315149689552,0.2738357450649768]
mediaScale: 1
cameraFOV: 56.52

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[296223.1131077501,4633741.99950086,128.49272754182056],[296225.76470118354,4633740.109690328,129.46157286325848]],
    [[296220.24550796876,4633744.50945153,129.7523040677012],[296233.63938660576,4633734.974536185,131.14492907758668]],
    [[296211.9160385982,4633752.941794123,128.76330820651427],[296225.3099172352,4633743.406878778,130.15593321639975]],
    [[296210.8050975578,4633759.7009269735,128.12330883281317],[296221.25695977156,4633747.09682167,130.15926211921146]],  
    [[296196.9022813706,4633772.666234221,143.2681648347127],[296207.68108440045,4633761.096691275,138.55508927143177]]
]



# ViewPoints overrides in seconds
animationEntry:
---
---
title: V Piramide

mediaPath: /videos/p_02_cl-1080p.mp4
mediaPosition:  [296280.5478984773,4633731.9264766555,132.01560618748866]
mediaRotation:  [0.5211164922501835,0.5365616331356262,0.47613167821252356,0.4624260377122162]
mediaScale: 1
cameraFOV: 38

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[296284.1208274853,4633732.03084898,131.58748853802086],[296267.7449028653,4633731.552475827,133.54969443141496]],
    [[296290.9070106481,4633732.229086811,130.77435041574967],[296274.53108602803,4633731.750713661,132.73655630917466]],
    [[296304.9369984471,4633731.350885622,129.16930210467007],[296288.56107382703,4633730.872512473,131.13150799809506]],
    [[296288.7609095558,4633719.230629459,130.8789841178707],[296279.33056099195,4633732.7691735905,130.71666108735906]],
    [[296270.3126982619,4633694.305162148,131.7248114587763],[296274.6320982583,4633710.229736301,131.6979469420832]],
    [[296277.4472505592,4633672.19475005,128.94914516272402],[296267.5576366572,4633685.380271203,129.71567553203505]],
    [[296271.0729524511,4633685.805003644,149.19390469141447],[296267.5977615166,4633701.033325689,143.87683259536482]]
]

# ViewPoints overrides in seconds
animationEntry:
---
---
title: V Piramide

mediaPath: /videos/p_06_anon741a-1080p.mp4
mediaPosition:  [296275.73523494444,4633679.405492418,128.8302519253152]
mediaRotation:  [-0.6765413323171713,-0.2533805698351245,-0.24251043310449974,-0.6475174146940006]
mediaScale: 1
cameraFOV: 37.80


# Pair of camera points and targets: [final point], ... , [entrance point]

cameraPath: [
    [[296278.09781891824,4633676.693793407,128.67250073529024],[296275.0395278921,4633680.204002925,128.87670471456474]],
    [[296281.48966916313,4633672.80073519,128.4460248011342],[296270.66115928313,4633685.229355658,129.1690510887317]],
    [[296308.21248647675,4633636.666559246,147.9405297230549],[296245.01757910114,4633707.009339188,118.83826615341259]],
    [[296329.05805837637,4633664.354821563,147.9405297230549],[296244.4031891765,4633708.404759425,121.82921324620742]]
    ]

# ViewPoints overrides in seconds
animationEntry:
---
---
title: V Piramide

mediaPath: /videos/p_12_cf1911-1080p.mp4
mediaPosition:  [296228.28400306916,4633736.400851762,129.2511144984435]
mediaRotation:  [0.5269840627127194,-0.5733266570238038,-0.46188710518731796,0.4245522865965815]
mediaScale: 1
cameraFOV: 36

 
# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[296224.778940207,4633736.69662705,128.484921760791],[296234.26508615044,4633735.896137113,130.558555363383]],
    [[296216.2473001603,4633739.185407818,129.35],[296313.5188318435,4633726.721185817,142.44721523640993]],[[296205.1414360284,4633765.775039359,147.53590011585013],[296270.0562206802,4633699.505066656,113.14104287640416]],
    [[296179.9163777299,4633777.488411805,157.42212858171678],[296253.3373124641,4633720.812268976,122.9863052298312]]

    ]

# ViewPoints overrides in seconds
animationEntry:
---
---
title: V Piramide

mediaPath: /videos/p_13_bdc1919-1938-1080p.mp4
mediaPosition:  [296444.4168692821,4633464.622316913,133.04165050154955]
mediaRotation:  [-0.7058849158827654,-0.22290816719440346,-0.20245971131340537,-0.6411306417744402]
mediaScale: 1
cameraFOV: 32

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[296446.47482017265,4633461.6887883255,132.6963272301455],[296445.2006596956,4633463.505054302,132.91013080837945]],
    [[296449.86177526566,4633456.962197441,135.55507299913518],[296441.5490049243,4633467.01842608,132.65794686975082]],
    [[296453.2761071372,4633452.83176913,135.02048877096098],[296444.7625299597,4633462.964742666,133.1590173201688]],
    [[296489.6976483115,4633417.304520369,154.06665320194477],[296420.2281258672,4633486.175786515,139.25520817250606]]
    
    ]

# ViewPoints overrides in seconds
animationEntry:
---
---
title: V Piramide

mediaPath: /videos/p_14_gb-1080p.mp4
mediaPosition:  [296251.0723864292,4633701.959845937,129.05533985330038]
mediaRotation:  [0.7574961606925863,0.022911060029451716,0.019724458258800245,0.6521392455685483]
mediaScale: 1
cameraFOV: 34

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[296251.2875397294,4633698.406350258,128.52019733232348],[296251.1543296354,4633700.606462889,128.851525671085]],
    [[296251.4199469527,4633696.219497924,128.19086594246775],[296250.62119032594,4633709.411850631,130.17758255141112]],
    [[296263.81043535325,4633686.548917581,129.09250804427694],[296254.06993061723,4633695.700111679,129.1304431159717]],
    [[296276.5741153363,4633676.950861242,129.04722874270257],[296267.9206320729,4633687.135989225,129.10739838190275]],
    [[296276.5741153363,4633676.950861242,131.56163074259638],[296268.0167220122,4633687.173005566,130.6118887203029]]
    ]

# ViewPoints overrides in seconds
animationEntry:
---
---
title: Casal Rotondo

mediaPath: /videos/cr_26_mgc-1080p.mp4
mediaPosition:  [296938.6415785213,4632814.261376878,141.84827250158813]
mediaRotation:  [0.4902833384537672,-0.5469759523397676,-0.5052827759600953,0.4529115497695916]
mediaScale: 1
cameraFOV: 27

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[296935.07424166315,4632814.652499173,141.56343730188388],[296946.9936059405,4632813.345661806,142.51514291289567]],
    [[296931.3204334786,4632815.183795392,142.67845089047543],[296962.2317035338,4632812.352763953,142.01822590398257]],
    [[296921.1463275563,4632810.456019402,142.89611449236503],[296951.9692335495,4632814.137263926,142.30100328416606]],
    [[296931.9561774807,4632786.702243705,144.16201490264905],[296937.77819878224,4632817.199018079,144.07216384627444]],
    [[296963.4006140729,4632783.454592137,144.035533257345],[296941.7356310198,4632805.445731369,140.72259542951483]],
    [[296991.7376807767,4632774.581124641,148.22909674616503],[296968.14340473595,4632794.259885966,143.75542450834527]]
]

# ViewPoints overrides in seconds
animationEntry:
---
---
title: Casal Rotondo

mediaPath: /videos/cr_05_fllia-doc1210-1080p.mp4
mediaPosition:  [297019.87674980354,4632737.526190996,138.8215683104031]
mediaRotation:  [0.6904881945537289,0.14467669514336814,0.14534220707145987,0.6936644360985517]
mediaScale: 1
cameraFOV: 33



# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[297021.32189172896,4632734.229026784,138.83809020126625],[297019.95815798827,4632737.3404541155,138.82249902670137]],
    [[297024.15081726253,4632727.774690459,138.87043249104664],[297017.52725010423,4632742.886693098,138.79470715791717]],
    [[297024.15081726253,4632727.774690459,141.36028249098106],[297016.98105967586,4632742.530521578,139.59675644431846]],
    [[297026.8600820387,4632725.514193867,142.02144446366157],[297018.94518198026,4632739.771389654,139.50415914700056]],
    [[297051.79608813475,4632693.035674304,142.84200482622578],[297033.8749770918,4632714.407104372,141.138734345702]]

]

# ViewPoints overrides in seconds
animationEntry:
---
---
title: Casal Rotondo

mediaPath: /videos/cr_07_da1890-1080p.mp4
mediaPosition:  [296992.80614868924,4632757.377810636,140.89388979648496]
mediaRotation:  [0.7178019860693331,-0.023519886912034965,-0.022788363937090814,0.6954766812651857]
mediaScale: 1
cameraFOV: 36

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[296992.5706002156,4632753.78732498,140.78018108722983],[296992.7326419592,4632756.257341224,140.8584051454395]],
    [[296991.13077652344,4632709.241409425,154.98623842666697],[296988.9124971244,4632736.407721615,148.8318828167722]],
    [[297016.2680613347,4632739.962280817,143.28097692473688],[297004.7234322613,4632751.736257562,142.69360040919153]]
]
# ViewPoints overrides in seconds
animationEntry:
---
---
title: Casal Rotondo

mediaPath: /videos/cr_02_cl-1080p.mp4
mediaPosition:  [296940.25140633545,4632826.750343467,140.72027158246507]
mediaRotation:  [0.380689948779622,-0.6420036086179064,-0.5724393453931962,0.33944031178627554]
mediaScale: 1
cameraFOV: 35

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[296937.11333090236,4632828.466007063,140.3091985942873],[296939.2721188875,4632827.28574418,140.59198957455848]],
    [[296932.0811813138,4632831.812906609,142.21091970844603],[296946.91576768836,4632824.999042007,139.8116914388798]],
    [[296931.8868545879,4632839.220447031,142.71093189641877],[296940.75511417293,4632825.948806095,138.53151581418015]],
    [[296934.23063188774,4632848.69832915,142.86152831102262],[296945.47197377984,4632836.844621213,140.54380515747292]]

]

# ViewPoints overrides in seconds
animationEntry:
---
---
title: Casal Rotondo

mediaPath: /videos/cr_21_eisner-1080p.mp4
mediaPosition:  [297042.0402477955,4632831.392324763,138.4429678854768]
mediaRotation:  [-0.33038338273322215,-0.6645971686853058,-0.6001281900051146,-0.29833467674813374]
mediaScale: 1
cameraFOV: 37.5

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[297044.8953700889,4632833.554331232,138.07690133675925],[297042.47051209875,4632831.718137149,138.38780199332618]],
    [[297066.97864064743,4632850.276618015,146.70605771555068],[297049.0915425064,4632836.731831519,143.72664752459883]],
    [[297063.5326200608,4632733.813554655,148.69419510064242],[297046.8394630481,4632758.9643535195,141.43244063852154]],
    [[297031.62895526626,4632726.193773426,142.28113532396532],[297061.1029200704,4632720.033477699,134.71189442812985]],
    [[297014.4037261628,4632744.996146135,145.60912590682233],[297037.2214220592,4632724.398956651,141.2426866498023]]
]

# ViewPoints overrides in seconds
animationEntry:
---
---
title: Casal Rotondo

mediaPath: /videos/cr_16_emigliorini-1080p.mp4
mediaPosition:  [296897.4718234203,4632892.744298877,139.73800978195376]
mediaRotation:  [0.3444127081651618,-0.6410650772826525,-0.6041923982105522,0.3246028328396219]
mediaScale: 1
cameraFOV: 37

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[296894.4753052426,4632894.7281098785,139.52500104458284],[296908.2093468905,4632885.63564279,140.5012910908663]],
    [[296887.2663067835,4632899.500745856,140.77969639948986],[296900.55461895117,4632889.744326621,140.08468312837303]],
    [[296885.55880507827,4632908.173292268,141.51172963182637],[296896.42998205806,4632895.879857564,139.79738594915565]],
    [[296882.5279003772,4632911.6007251125,142.06778879462277],[296893.48417852377,4632899.425043571,140.077123534015]]
]

# ViewPoints overrides in seconds
animationEntry:
---
---
title: Casal Rotondo

mediaPath: /videos/cr_08_ppm-1847-1080p.mp4
mediaPosition:  [297022.2361005747,4632729.858347014,139.35180480328142]
mediaRotation:  [0.7169534213120289,0.12474910230558983,0.11757380329640904,0.6757158085477527]
mediaScale: 1
cameraFOV: 28

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[297023.44994771836,4632726.475862366,139.1387960659105],[297019.79975762655,4632736.647416512,139.77933996527707]],
    [[297027.7593585449,4632719.184281134,139.4668934001318],[297019.7647945776,4632733.592547956,138.60736203711062]],
    [[297036.48127107334,4632715.165452657,140.10591447959794],[297025.4152387121,4632727.286040263,138.40759652789646]],
    [[297040.86617526057,4632710.23647381,142.1659182797483],[297030.13999433786,4632722.54116082,139.75890344978362]]

]

# ViewPoints overrides in seconds
animationEntry:
---
---
title: Casal Rotondo

mediaPath: /videos/cr_20_pdea2-1080p.mp4
mediaPosition:  [296988.41041277576,4632786.868594005,139.81255060742663]
mediaRotation:  [0.7077279743155134,0.22557484952108756,0.2033144892388494,0.6378873882593575]
mediaScale: 1
cameraFOV: 27

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[296990.48244823905,4632783.948359585,139.43985773919718],[296989.80067278136,4632784.909223525,139.56248732508965]],
    [[296994.9026866668,4632775.896996137,138.55998279576255],[296986.79014326265,4632790.199024034,139.93425989406063]],
    [[296993.770438936,4632771.888035351,138.44427233359963],[296986.8588106399,4632786.606968936,141.2429127978262]],
    [[296995.2133827239,4632768.121093521,143.2745514316362],[296989.491908517,4632783.565760057,142.2861140451723]],
    [[296999.73307219637,4632761.150520475,142.69859870496205],[296989.9462331358,4632774.415387591,141.9836990019136]]
]

# ViewPoints overrides in seconds
animationEntry:
---
---
title: Casal Rotondo

mediaPath: /videos/cr_24_en-1853-1080p.mp4
mediaPosition:  [296982.9840754576,4632793.396584391,139.1961522650109]
mediaRotation:  [0.26304208287991926,-0.6486129893057708,-0.6618633284881863,0.26841569838696994]
mediaScale: 1
cameraFOV: 40

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[296980.47706557513,4632795.9791417895,139.26894464472468],[296983.1057147102,4632793.271279599,139.19262040391982]],
    [[296974.4782927671,4632801.782346546,139.4355593091544],[296996.79272524954,4632780.196456326,139.165931986716]],
    [[296967.87411919114,4632806.530421087,143.5103660780395],[296993.3366474464,4632789.056662488,140.3033943380557]],
    [[296954.40687568067,4632824.651575352,146.43775084235645],[296981.64374177397,4632809.7949953,145.2577478629383]],
    [[296947.3031811753,4632830.036821914,144.53189439219955],[296970.4051005754,4632809.77989518,140.0692167448114]]
]

# ViewPoints overrides in seconds
animationEntry:
---
---
title: Casal Rotondo

mediaPath: /videos/cr_17_italynews-1080p.mp4
mediaPosition:  [296948.54899371124,4632815.700116255,141.95114097205806]
mediaRotation:  [0.44116586282317793,-0.5755912987286951,-0.5464735750191229,0.41884838559210674]
mediaScale: 1
cameraFOV: 30

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[296945.0773627095,4632816.634415798,141.7644258078706],[296951.9268028871,4632814.791066427,142.13281010525276]],
    [[296942.8441432829,4632817.261879921,142.85668223796264],[296958.8809208125,4632813.459505282,142.073283505665]],
    [[296935.9182597427,4632810.895433801,143.25389469319808],[296951.5558096746,4632816.150468802,142.93508554907604]],
    [[296948.77800104296,4632801.739053373,143.2673396235779],[296945.05171989277,4632817.774421008,142.15746864440183]],
    [[296968.99472884607,4632801.681995601,143.93350756527585],[296957.0088012776,4632812.988414044,143.0660732224691]]
]

# ViewPoints overrides in seconds
animationEntry:
---
---
title: Casal Rotondo

mediaPath: /videos/cr_13_bulwerc1890-1080p.mp4
mediaPosition:  [296840.60913562856,4632939.949020138,141.0024421187111]
mediaRotation:  [-0.34891744078243564,0.6533198908811648,0.5926589324480149,-0.3165203460860577]
mediaScale: 1
cameraFOV: 28

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[296837.6313734812,4632941.941666217,140.652736244459],[296851.27944999025,4632932.808705023,142.25555483481435]],
    [[296838.1211994843,4632946.6545256395,142.09308883496146],[296847.5171994076,4632933.427808177,139.0899020955388]],
    [[296795.03040742426,4632965.46301966,151.67603819120316],[296818.4227642658,4632950.434993112,148.89010316465684]]
]

# ViewPoints overrides in seconds
animationEntry:
---
---
title: Casal Rotondo

mediaPath: /videos/cr_22_en-709-1080p.mp4
mediaPosition:  [296941.6250700751,4632820.468759457,142.23520762403857]
mediaRotation:  [0.47488154263549487,-0.5603998258988842,-0.5176834163668135,0.4386837539907507]
mediaScale: 1
cameraFOV: 42

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[296938.08499856427,4632821.057624983,141.95037242433432],[296938.8492323832,4632820.930500164,142.01186290320615]],
    [[296928.83556078677,4632819.072675213,143.86270126422096],[296959.3533318016,4632823.593737778,140.3723004234844]],
    [[296930.1086223235,4632813.348461116,143.86270126422096],[296959.8571000068,4632821.521736466,140.3723004234844]],
    [[296953.2284948294,4632798.856412477,143.66547478011697],[296926.74535094836,4632814.646514234,140.02225791351054]],
    [[296963.5849768292,4632813.104025575,145.51949596322325],[296947.2231413561,4632813.141044832,143.38901021459324]]
]

# ViewPoints overrides in seconds
animationEntry:
---
---
title: Casal Rotondo

mediaPath: /videos/cr_10_tha-1080p.mp4
mediaPosition:  [296999.116136317,4632755.867208812,139.5253647007758]
mediaRotation:  [0.7575015593183931,0.022731866748785107,0.019570188208430577,0.6521438933224821]
mediaScale: 1
cameraFOV: 34

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[296999.32960836944,4632752.313611738,138.9902221797989],[296999.1974392038,4632754.513787148,139.32155051856043]],
    [[296999.8150838821,4632748.845666032,141.2522781457522],[296997.4008826179,4632765.056749652,139.34901775096296]],
    [[297002.2479055656,4632745.151626623,141.90106673592578],[296995.7694824099,4632760.19164277,139.8816738266538]],
    [[297017.48716844106,4632738.939378372,141.75163225987578],[297005.26866676664,4632749.943659025,140.38633956068483]]
]
# ViewPoints overrides in seconds
animationEntry:
---
---
title: Casal Rotondo

mediaPath: /videos/cr_01_gbp-1080p.mp4
mediaPosition:  [296999.88909307844,4632755.697524181,139.4887531447218]
mediaRotation:  [0.7506203763417982,0.031484881582427446,0.027658769019413307,0.6594033252483404]
mediaScale: 1
cameraFOV: 37.5

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[297000.1880548712,4632752.140070859,139.02491296525793],[297000.05311747047,4632753.745739263,139.23426877492443]],
    [[297000.71852076065,4632748.353910485,142.38558460351166],[296997.13660481834,4632764.1701520225,139.3415113747053]],
    [[297012.74727664556,4632746.154566575,142.04764628845808],[297001.90292390797,4632758.451862522,140.19655848057053]],
    [[297022.856829675,4632731.315122587,141.80269992873437],[297011.74477223965,4632743.370240487,139.94627541793565]]
]

# ViewPoints overrides in seconds
animationEntry:
---

---
title: Casal Rotondo

mediaPath: /videos/cr_12_anonblew1900-1080p.mp4
mediaPosition:  [296951.6757288243,4632812.281750991,140.66699605290427]
mediaRotation:  [0.6529803341290381,-0.3252515338186515,-0.3049555599825983,0.6122338029095127]
mediaScale: 1
cameraFOV: 32

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[296948.80825706257,4632810.11750753,140.43535853842755],[296961.9508359704,4632820.036956726,141.4970304797792]],
    [[296931.3460996595,4632796.937827619,147.8202381873429],[296951.6089739799,4632814.956819686,141.0721107703129]],
    [[296949.0877734908,4632828.747883824,141.9093395092039],[296959.0810252713,4632815.810355066,139.67209376937654]],
    [[296940.14574530587,4632836.364715805,141.96354364843995],[296954.4690298529,4632828.228500003,141.01722858128144]]
]

# ViewPoints overrides in seconds
animationEntry:
---
---
title: VII Casal Rotondo

mediaPath: /videos/cr_06_fllia1890-1080p.mp4
mediaPosition:  [296945.33494263235,4632822.7629878,140.42717086617222]
mediaRotation:  [0.4480481622342859,-0.5710478576985801,-0.5411712277450123,0.4246067484123368]
mediaScale: 1
cameraFOV: 37

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[296941.8433634861,4632823.618278658,140.2339026236323],[296957.84643457306,4632819.698195558,141.11971540194173]],
    [[296917.62457231816,4632838.862155331,146.70664014891025],[296941.5647346665,4632827.257723842,138.16198269355502]],
    [[296942.0076741443,4632835.074815828,141.71321109268627],[296952.9107751012,4632822.952919936,139.17706648045015]]
]

# ViewPoints overrides in seconds
animationEntry:
---
---
title: Casal Rotondo

mediaPath: /videos/cr_03_jgc-1080p.mp4
mediaPosition:  [296910.29734363017,4632873.137401671,138.89428568234854]
mediaRotation:  [0.24867812551193538,-0.6795761917684637,-0.6481387201983533,0.23717417408517116]
mediaScale: 1
cameraFOV: 20

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[296907.9763815532,4632875.884054894,138.7239004865566],[296918.6141244059,4632863.295227625,139.5048326339363]],
    [[296905.11731086404,4632879.267512034,138.51401195628165],[296915.7550537202,4632866.678684769,139.2949441039151]],
    [[296903.8238952424,4632885.203461009,140.17666423999302],[296913.81970036187,4632872.189995036,138.44935899701377]]
]

# ViewPoints overrides in seconds
animationEntry:
---
---
title: Casal Rotondo

mediaPath: /videos/cr_14_bs-1080p.mp4
mediaPosition:  [296946.24335131433,4632826.098027713,140.38710457486843]
mediaRotation:  [0.36627962705513023,-0.6433538844682531,-0.5842139109329754,0.33260956152275883]
mediaScale: 1
cameraFOV: 37

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[296943.1619579051,4632827.927030629,140.04103667670498],[296957.2850110306,4632819.544100596,141.62718120995416]],
    [[296938.76074378385,4632830.539431062,141.15659116872428],[296952.9515042663,4632822.288113024,139.48684041730806]],
    [[296936.8212888076,4632837.788719209,141.44046696506922],[296948.9179778256,4632826.59441934,140.6608739688652]],
    [[296936.3653954944,4632843.744913024,141.9552685092986],[296947.80377140077,4632831.956787983,140.3888507324399]]
]

# ViewPoints overrides in seconds
animationEntry:
---
---
cameraFOV: 60

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
  [[296910.8632555868,4632739.143639741,204.55082802089478],[296995.64606282354,4632807.08486307,127.26147851243842]],
  [[296859.0044384939,4632697.586258139,251.82614680356733],[296943.78724573064,4632765.527481468,174.53679729511094]]
]
---
---
title: Casal Rotondo

mediaPath: /videos/cr_15_1905-1080p.mp4
mediaPosition:  [296947.5430050464,4632811.474955375,141.67715795896848]
mediaRotation:  [0.44304323974298676,-0.5201681930755518,-0.5558647327872339,0.47344707990854784]
mediaScale: 1
cameraFOV: 50

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[296943.9966866325,4632812.046531879,141.91574967710847],[296960.2506460293,4632809.426806238,140.82220430230018]],
    [[296940.91759122594,4632813.963403654,142.9865716458081],[296957.17155062186,4632811.343678007,141.8930262710389]],
    [[296934.1423332583,4632818.020357948,142.92845803625053],[296949.4048663844,4632812.8489217535,139.3842613345593]],
    [[296942.1083607127,4632830.520968929,142.28484329491712],[296936.45364891493,4632815.020386914,142.20615611135236]],
    [[296945.28474015585,4632834.443087063,142.39849340998416],[296957.3078695127,4632823.212120793,141.1495910466104]]
]

# ViewPoints overrides in seconds
animationEntry:
---
---
title: Casal Rotondo

mediaPath: /videos/cr_04_anonrs1869-1080p.mp4
mediaPosition:  [296938.8482914974,4632836.330151673,140.6673264669091]
mediaRotation:  [0.322657681428531,-0.6461971424269202,-0.618761549317856,0.30895860373230244]
mediaScale: 1
cameraFOV: 25

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[296935.9733578946,4632838.491261037,140.5112396377971],[296949.15013690735,4632828.586176454,141.22663760456035]],
    [[296932.39354685205,4632843.9537535375,140.24461081560142],[296945.5703258646,4632834.0486689545,140.96000878234898]],
    [[296929.02994738315,4632851.679148688,141.55751482106473],[296940.1087432523,4632839.61393968,139.57231616668318]],
    [[296924.59937075025,4632856.634933943,145.8210246082357],[296935.4182836348,4632844.672376672,142.34275671577868]]

    ]

# ViewPoints overrides in seconds
animationEntry:
---
---
title: Casal Rotondo

mediaPath: /videos/cr_23_en-5851-1080p.mp4
mediaPosition:  [297009.4566497859,4632746.912722128,139.92593455589463]
mediaRotation:  [0.6964737685090627,0.17564793857501346,0.17013887976237105,0.6746294190489608]
mediaScale: 1
cameraFOV: 27

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[297011.16301042726,4632743.744886646,139.81125361435718],[297010.4698560251,4632745.031718431,139.85783907673988]],
    [[297015.1246254279,4632736.390201904,142.17473861416656],[296998.63487533707,4632762.088837383,136.5510806519757]],
    [[297023.15091658203,4632731.423439047,141.93751178150472],[297001.8972668678,4632753.974080365,140.0118424353924]],
    [[297044.2554480555,4632705.88179909,146.25491248976672],[297024.4702099611,4632728.617290341,138.7984358714868]]
]

# ViewPoints overrides in seconds
animationEntry:
---
---
title: Casal Rotondo

mediaPath: /videos/cr_09_ppm-0159-1080p.mp4
mediaPosition:  [296934.72303731187,4632816.67377751,142.76696397696765]
mediaRotation:  [0.45243002943162186,-0.5655961274846097,-0.5384295861974208,0.43069904777023715]
mediaScale: 1
cameraFOV: 25

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[296931.215172637,4632817.463419554,142.58990221807665],[296939.75952856085,4632815.540032167,143.0211841598578]],
    [[296924.4350279619,4632816.680167795,144.26897258359915],[296940.8415879849,4632816.859497718,142.5246592541178]],
    [[296914.87994132796,4632816.788384672,145.41877017778373],[296931.26394754136,4632816.03781444,143.6157214467935]],
    [[296912.51318480907,4632839.809527022,144.08578871128137],[296923.6876265672,4632827.751859346,142.67352045418983]],
    [[296929.28785643005,4632856.122564406,142.63299224381888],[296936.3369849388,4632841.232928367,141.70642551429566]]
]

# ViewPoints overrides in seconds
animationEntry:
---
---
title: Casal Rotondo

mediaPath: /videos/cr_18_gs-1080p.mp4
mediaPosition:  [296932.3752148676,4632852.946841958,139.59366356283024]
mediaRotation:  [0.2949741505626587,-0.6578125175629119,-0.6323514659366846,0.28355698856686234]
mediaScale: 1
cameraFOV: 29

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[296929.68922122166,4632855.3395944,139.4516287469394],[296938.66379749106,4632847.344810305,139.92620257553472]],
    [[296924.859282369,4632859.642228545,140.92102252985873],[296937.17795360775,4632848.668438373,140.6448555709424]],
    [[296915.7724560902,4632864.540355387,141.27595109592687],[296930.0094525602,4632856.207771076,140.92102827417796]],
    [[296908.846294283,4632877.874034729,142.49937349849347],[296920.51275466796,4632866.206351424,142.40520138860677]]
]

# ViewPoints overrides in seconds
animationEntry:
---
---
title: Casal Rotondo

mediaPath: /videos/cr_25_en-5855-1080p.mp4
mediaPosition:  [296955.05684507865,4632827.278809856,140.79108258427547]
mediaRotation:  [0.40193969704834226,-0.6693276661510459,-0.5356883722730923,0.3216876172203178]
mediaScale: 1
cameraFOV: 36

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[296951.9563174009,4632828.929424433,140.00228626215812],[296952.78266888397,4632828.48950325,140.21251596444185]],
    [[296946.47255002504,4632831.857398947,141.3291502247845],[296974.12773226126,4632817.7460664995,141.4649607762072]],
    [[296935.4335354091,4632843.608299404,142.94429884399392],[296959.63616468635,4632824.556412461,139.04358192998174]]
    
]

# ViewPoints overrides in seconds
animationEntry:
---
---
title: Casal Rotondo

mediaPath: /videos/cr_11_anonred1890-1080p.mp4
mediaPosition:  [296948.1795287538,4632823.433635699,140.86566347987733]
mediaRotation:  [0.43596857247208,-0.5669443343978027,-0.554057306734503,0.42605871234490433]
mediaScale: 1
cameraFOV: 38

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[296944.70118610136,4632824.357910084,140.78290326355574],[296960.64358992496,4632820.121652485,141.1622209216964]],
    [[296938.573951924,4632827.551041896,140.6278761366487],[296954.51635574887,4632823.314784303,141.00719379482368]],
    [[296938.21781125996,4632841.324638402,140.71607353722277],[296949.6532896911,4632829.431977375,140.50524951954952]],
    [[296932.796768572,4632848.540265019,140.4979105547467],[296944.5862749472,4632837.004905139,140.0585793820419]]
]

# ViewPoints overrides in seconds
animationEntry:
---
---
title: Casal Rotondo

mediaPath: /videos/cr_19_pdea1-1080p.mp4
mediaPosition:  [296970.1768384401,4632793.087797696,139.83685766520327]
mediaRotation:  [0.70035740356016,0.1419694741444417,0.13897543455484554,0.6855873425610436]
mediaScale: 1
cameraFOV: 40

# Pair of camera points and targets: [final point], ... , [entrance point]
cameraPath: [
    [[296971.57842607296,4632789.772731203,139.76013573016854],[296964.440711276,4632806.655014268,140.15084928821574]],
    [[296973.48505245,4632785.26313575,139.6557683245132],[296966.52802230365,4632800.201164031,138.8161001670592]],
    [[296983.5552766577,4632780.577079458,140.21759741081414],[296971.15558872465,4632791.388312938,138.94678278248662]],
    [[296988.1637632812,4632777.681967447,142.57494963117006],[296976.48876867833,4632789.0246057585,139.8751242388577]]
]

# ViewPoints overrides in seconds
animationEntry:
---
