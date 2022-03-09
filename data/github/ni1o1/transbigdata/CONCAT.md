# Contributor Covenant Code of Conduct

## Our Pledge

We as members, contributors, and leaders pledge to make participation in our
community a harassment-free experience for everyone, regardless of age, body
size, visible or invisible disability, ethnicity, sex characteristics, gender
identity and expression, level of experience, education, socio-economic status,
nationality, personal appearance, race, religion, or sexual identity
and orientation.

We pledge to act and interact in ways that contribute to an open, welcoming,
diverse, inclusive, and healthy community.

## Our Standards

Examples of behavior that contributes to a positive environment for our
community include:

* Demonstrating empathy and kindness toward other people
* Being respectful of differing opinions, viewpoints, and experiences
* Giving and gracefully accepting constructive feedback
* Accepting responsibility and apologizing to those affected by our mistakes,
  and learning from the experience
* Focusing on what is best not just for us as individuals, but for the
  overall community

Examples of unacceptable behavior include:

* The use of sexualized language or imagery, and sexual attention or
  advances of any kind
* Trolling, insulting or derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or email
  address, without their explicit permission
* Other conduct which could reasonably be considered inappropriate in a
  professional setting

## Enforcement Responsibilities

Community leaders are responsible for clarifying and enforcing our standards of
acceptable behavior and will take appropriate and fair corrective action in
response to any behavior that they deem inappropriate, threatening, offensive,
or harmful.

Community leaders have the right and responsibility to remove, edit, or reject
comments, commits, code, wiki edits, issues, and other contributions that are
not aligned to this Code of Conduct, and will communicate reasons for moderation
decisions when appropriate.

## Scope

This Code of Conduct applies within all community spaces, and also applies when
an individual is officially representing the community in public spaces.
Examples of representing our community include using an official e-mail address,
posting via an official social media account, or acting as an appointed
representative at an online or offline event.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported to the community leaders responsible for enforcement at
1810514@tongji.edu.cn.
All complaints will be reviewed and investigated promptly and fairly.

All community leaders are obligated to respect the privacy and security of the
reporter of any incident.

## Enforcement Guidelines

Community leaders will follow these Community Impact Guidelines in determining
the consequences for any action they deem in violation of this Code of Conduct:

### 1. Correction

**Community Impact**: Use of inappropriate language or other behavior deemed
unprofessional or unwelcome in the community.

**Consequence**: A private, written warning from community leaders, providing
clarity around the nature of the violation and an explanation of why the
behavior was inappropriate. A public apology may be requested.

### 2. Warning

**Community Impact**: A violation through a single incident or series
of actions.

**Consequence**: A warning with consequences for continued behavior. No
interaction with the people involved, including unsolicited interaction with
those enforcing the Code of Conduct, for a specified period of time. This
includes avoiding interactions in community spaces as well as external channels
like social media. Violating these terms may lead to a temporary or
permanent ban.

### 3. Temporary Ban

**Community Impact**: A serious violation of community standards, including
sustained inappropriate behavior.

**Consequence**: A temporary ban from any sort of interaction or public
communication with the community for a specified period of time. No public or
private interaction with the people involved, including unsolicited interaction
with those enforcing the Code of Conduct, is allowed during this period.
Violating these terms may lead to a permanent ban.

### 4. Permanent Ban

**Community Impact**: Demonstrating a pattern of violation of community
standards, including sustained inappropriate behavior,  harassment of an
individual, or aggression toward or disparagement of classes of individuals.

**Consequence**: A permanent ban from any sort of public interaction within
the community.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage],
version 2.0, available at
https://www.contributor-covenant.org/version/2/0/code_of_conduct.html.

Community Impact Guidelines were inspired by [Mozilla's code of conduct
enforcement ladder](https://github.com/mozilla/diversity).

[homepage]: https://www.contributor-covenant.org

For answers to common questions about this code of conduct, see the FAQ at
https://www.contributor-covenant.org/faq. Translations are available at
https://www.contributor-covenant.org/translations.
[English](README.md) 中文版

# TransBigData 针对交通时空大数据处理的Python包

<img src="https://github.com/ni1o1/transbigdata/raw/main/docs/source/_static/logo-wordmark-dark.png" style="width:550px">

[![Documentation Status](https://readthedocs.org/projects/transbigdata/badge/?version=latest)](https://transbigdata.readthedocs.io/en/latest/?badge=latest) [![PyPI version](https://badge.fury.io/py/transbigdata.svg)](https://badge.fury.io/py/transbigdata) [![Downloads](https://pepy.tech/badge/transbigdata)](https://pepy.tech/project/transbigdata) ![GitHub commit activity](https://img.shields.io/github/commit-activity/m/ni1o1/transbigdata) [![bilibili](https://img.shields.io/badge/bilibili-%E5%90%8C%E6%B5%8E%E5%B0%8F%E6%97%AD%E5%AD%A6%E9%95%BF-green.svg)](https://space.bilibili.com/3051484) [![status](https://joss.theoj.org/papers/d1055fe3105dfa2dcff4cb6c7688a79b/status.svg)](https://joss.theoj.org/papers/d1055fe3105dfa2dcff4cb6c7688a79b) [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/ni1o1/transbigdata/d7d6fa33ff16440ba1698b10dd3cf3f76ff00abd?urlpath=lab%2Ftree%2Fexample%2FExample%201-Taxi%20GPS%20data%20processing.ipynb) [![Tests](https://github.com/ni1o1/transbigdata/actions/workflows/tests.yml/badge.svg)](https://github.com/ni1o1/transbigdata/actions/workflows/tests.yml) [![codecov](https://codecov.io/gh/ni1o1/transbigdata/branch/main/graph/badge.svg?token=GLAVYYCD9L)](https://codecov.io/gh/ni1o1/transbigdata) [![DOI](https://zenodo.org/badge/419559811.svg)](https://zenodo.org/badge/latestdoi/419559811) [![Gitter](https://badges.gitter.im/transbigdata/community.svg)](https://gitter.im/transbigdata/community?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge)



`TransBigData`是一个为交通时空大数据处理、分析和可视化而开发的Python包。`TransBigData`为处理常见的交通时空大数据（如出租车GPS数据、共享单车数据和公交车GPS数据）提供了快速而简洁的方法。`TransBigData`为交通时空大数据分析的各个阶段提供了多种处理方法,代码简洁、高效、灵活、易用，可以用简洁的代码实现复杂的数据任务。

对于一些特定类型的数据，`TransBigData`还提供了针对特定需求的工具，如从出租车GPS数据中提取出租车行程的起点和终点信息（OD），从公交车GPS数据中识别到离站信息。该包的最新稳定版本可以通过pip安装，完整的文档可以查看：[TransBigData的说明文档](https://transbigdata.readthedocs.io/zh_CN/latest/)

**技术特点**

* 面向交通时空大数据分析不同阶段的处理需求提供不同处理功能。
* 代码简洁、高效、灵活、易用，通过简短的代码即可实现复杂的数据任务。


**主要功能**

目前，`TransBigData`主要提供以下方法:

* **数据质量分析**: 提供快速获取数据集一般信息的方法，包括数据量、时间段和采样间隔。
* **数据预处理**: 提供清洗多种类型的数据错误的方法。
* **数据栅格化**: 提供在研究区域内生成多种类型的地理网格（矩形网格、六角形网格）的方法。提供快速算法将GPS数据映射到生成的网格上。
* **数据聚合集计**: 提供将GPS数据和OD数据聚合到地理多边形的方法。
* **数据可视化**: 内置的可视化功能，利用可视化包keplergl，用简单的代码在Jupyter笔记本上交互式地可视化数据。
* **轨迹数据处理**: 提供处理轨迹数据的方法，包括从GPS点生成轨迹线型，轨迹增密等。
* **地图底图**: 提供在matplotlib上显示Mapbox地图底图的方法。

## 安装

### 用pypi安装

在安装`TransBigData`之前，请确保已经安装了可用的geopandas包：https://geopandas.org/index.html  
如果你已经安装了geopandas，则直接在命令提示符中运行下面代码即可安装：

    pip install -U transbigdata

### 用conda-forge安装

你也可以用conda-forge安装`TransBigData`，这种方式会自动解决环境依赖，不过国内可能需要更换conda源。运行下面代码即可安装：

    conda install -c conda-forge transbigdata

## 可视化示例

### 可视化轨迹(基于keplergl)

![gif](images/tbdexample1.gif)

### 可视化数据分布(基于keplergl)

![gif](images/tbdexample2.gif)

### 可视化OD(基于keplergl)

![gif](images/tbdexample3.gif)

## 使用示例

下面例子展示如何使用`TransBigData`工具快速处理出租车GPS数据，实现数据栅格化，数据聚合集计与数据可视化:

```python
import transbigdata as tbd
import pandas as pd
#读取出租车GPS数据 
data = pd.read_csv('TaxiData-Sample.csv',header = None) 
data.columns = ['VehicleNum','time','lon','lat','OpenStatus','Speed'] 
data
```

<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>VehicleNum</th>
      <th>time</th>
      <th>lon</th>
      <th>lat</th>
      <th>OpenStatus</th>
      <th>Speed</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>34745</td>
      <td>20:27:43</td>
      <td>113.806847</td>
      <td>22.623249</td>
      <td>1</td>
      <td>27</td>
    </tr>
    <tr>
      <th>1</th>
      <td>34745</td>
      <td>20:24:07</td>
      <td>113.809898</td>
      <td>22.627399</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>2</th>
      <td>34745</td>
      <td>20:24:27</td>
      <td>113.809898</td>
      <td>22.627399</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>3</th>
      <td>34745</td>
      <td>20:22:07</td>
      <td>113.811348</td>
      <td>22.628067</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>4</th>
      <td>34745</td>
      <td>20:10:06</td>
      <td>113.819885</td>
      <td>22.647800</td>
      <td>0</td>
      <td>54</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>544994</th>
      <td>28265</td>
      <td>21:35:13</td>
      <td>114.321503</td>
      <td>22.709499</td>
      <td>0</td>
      <td>18</td>
    </tr>
    <tr>
      <th>544995</th>
      <td>28265</td>
      <td>09:08:02</td>
      <td>114.322701</td>
      <td>22.681700</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>544996</th>
      <td>28265</td>
      <td>09:14:31</td>
      <td>114.336700</td>
      <td>22.690100</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>544997</th>
      <td>28265</td>
      <td>21:19:12</td>
      <td>114.352600</td>
      <td>22.728399</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>544998</th>
      <td>28265</td>
      <td>19:08:06</td>
      <td>114.137703</td>
      <td>22.621700</td>
      <td>0</td>
      <td>0</td>
    </tr>
  </tbody>
</table>
<p>544999 rows × 6 columns</p>
</div>


### 数据预处理

首先定义研究范围，并使用`tbd.clean_outofbounds`剔除研究范围外的数据

```python
#定义研究范围
bounds = [113.75, 22.4, 114.62, 22.86]
#剔除研究范围外的数据
data = tbd.clean_outofbounds(data,bounds = bounds,col = ['lon','lat'])
```

### 数据栅格化

以栅格形式表达数据分布是最基本的表达方法。GPS数据经过栅格化后，每个数据点都含有对应的栅格信息，采用栅格表达数据的分布时，其表示的分布情况与真实情况接近。如果要使用`TransBigData`工具进行栅格划分，首先需要确定栅格化的参数（可以理解为定义了一个栅格坐标系），参数可以帮助我们快速进行栅格化:

```python
#获取栅格化参数
params = tbd.grid_params(bounds,accuracy = 1000)
```

取得栅格化参数后，将GPS对应至栅格。使用`tbd.GPS_to_grids`方法,该方法会生成`LONCOL`列与`LATCOL`列，并由这两列共同指定一个栅格:

```python
#将GPS数据对应至栅格
data['LONCOL'],data['LATCOL'] = tbd.GPS_to_grids(data['lon'],data['lat'],params)
```

聚合集计栅格内数据量，并为栅格生成几何图形：

```python
#聚合集计栅格内数据量
grid_agg = data.groupby(['LONCOL','LATCOL'])['VehicleNum'].count().reset_index()
#生成栅格的几何图形
grid_agg['geometry'] = tbd.gridid_to_polygon(grid_agg['LONCOL'],grid_agg['LATCOL'],params)
#转换为GeoDataFrame
import geopandas as gpd
grid_agg = gpd.GeoDataFrame(grid_agg)
#绘制栅格
grid_agg.plot(column = 'VehicleNum',cmap = 'autumn_r')
```

![png](images/output_5_1.png)
    

### 数据可视化(在matplotlib中绘制地图底图)

对于一个正式的数据可视化图来说，我们还需要添加底图、色条、指北针和比例尺。 用`tbd.plot_map`加载地图底图，并用`tbd.plotscale`添加指北针和比例尺:

```python
import matplotlib.pyplot as plt
fig =plt.figure(1,(8,8),dpi=300)
ax =plt.subplot(111)
plt.sca(ax)
#加载地图底图
tbd.plot_map(plt,bounds,zoom = 11,style = 4)
#定义色条位置
cax = plt.axes([0.05, 0.33, 0.02, 0.3])
plt.title('Data count')
plt.sca(ax)
#绘制数据
grid_agg.plot(column = 'VehicleNum',cmap = 'autumn_r',ax = ax,cax = cax,legend = True)
#添加指北针和比例尺
tbd.plotscale(ax,bounds = bounds,textsize = 10,compasssize = 1,accuracy = 2000,rect = [0.06,0.03],zorder = 10)
plt.axis('off')
plt.xlim(bounds[0],bounds[2])
plt.ylim(bounds[1],bounds[3])
plt.show()
```

    
![png](images/output_7_0.png)
    


## 相关链接

* 小旭学长的b站： https://space.bilibili.com/3051484
* 小旭学长的七天入门交通时空大数据分析课程（零基础免费课）： https://www.lifangshuju.com/#/introduce/166  
* 小旭学长的交通时空大数据分析课程： https://www.lifangshuju.com/#/introduce/154  
* 小旭学长的数据可视化课程： https://www.lifangshuju.com/#/introduce/165  
* 本项目的github页面： https://github.com/ni1o1/transbigdata/  
* 有bug请在这个页面提交： https://github.com/ni1o1/transbigdata/issues  

## 引用信息
如果你想要引用`TransBigData`，请引用[这个DOI](https://doi.org/10.5281/zenodo.5912101)，引用信息在这个文件中[CITATION.cff](https://github.com/ni1o1/transbigdata/blob/main/CITATION.cff)。

## 介绍视频

* [Bilibili](https://www.bilibili.com/video/BV1nq4y1u7i1)
* [Youtube](https://www.youtube.com/watch?v=V_KHFv75W_w)
# Contributing to TransBigData

Whether you are a novice or experienced software developer, all contributions and suggestions are welcome!

## Getting Started

If you are looking to contribute to the *TransBigData* codebase, the best place to start is the [GitHub &#34;issues&#34; tab](https://github.com/ni1o1/transbigdata/issues). This is also a great place for filing bug reports and making suggestions for ways in which we can improve the code and documentation.

## Step-by-step Instructions of Contribute

The code is hosted on [GitHub](https://github.com/ni1o1/transbigdata),
so you will need to use [Git](http://git-scm.com/) to clone the project and make
changes to the codebase. 

1. Fork the [TransBigData repository](https://github.com/ni1o1/transbigdata).
2. Create a new branch from the `TransBigData` master branch.
3. Within your forked copy, the source code of `TransBigData` is located at the [src](https://github.com/ni1o1/transbigdata/tree/main/src) folder, you can make and test changes in the source code.
4. Before submitting your changes for review, make sure to check that your changes do not break any tests by running: ``pytest``. The tests are located in the [tests](https://github.com/ni1o1/transbigdata/tree/main/src/transbigdata/tests) folder.
5. When you are ready to submit your contribution, raise the Pull Request(PR). Once you finished your PR, the github [testing workflow](https://github.com/ni1o1/transbigdata/actions/workflows/tests.yml) will test your code. We will review your changes, and might ask you to make additional changes before it is finally ready to merge. However, once it's ready, we will merge it, and you will have successfully contributed to the codebase!

# 为TransBigData贡献代码

无论您是新手还是经验丰富的软件开发人员，欢迎您提供所有意见和建议！

## 开始

如果你想为*TransBigData*代码库做贡献，最好从[GitHub issues](https://github.com/ni1o1/transbigdata/issues)开始。你可以在这里提交BUG报告，并提出改进代码和文档的方法和建议。

## 如何贡献代码

代码托管在[GitHub](https://github.com/ni1o1/transbigdata)，所以你需要使用[Git](http://git-scm.com/)克隆项目并对代码做出更改。具体方法如下：
1. Fork [`TransBigData`仓库](https://github.com/ni1o1/transbigdata).
2. 以`TransBigData`的`main`分支为基础创建新分支。
3. 在您的分支仓库中，`TransBigData`的源代码位于[src](https://github.com/ni1o1/transbigdata/tree/main/src)文件夹，您可以在源代码中进行和测试更改，如果你使用的是jupyter notebook,可以在src文件夹下建立ipynb文件进行调试，这样修改transbigdata的源码时可以直接读取到。
4. 在提交更改以供审阅之前，请运行`pytest`来测试代码，确保您对代码的更改不会破坏任何测试结果。测试代码位于[tests](https://github.com/ni1o1/transbigdata/tree/main/src/transbigdata/tests)文件夹中
5. 当你准备好提交你的贡献时，提交Pull Request（PR）。完成PR后，github提供的[测试工作流](https://github.com/ni1o1/transbigdata/actions/workflows/tests.yml)将测试您的代码，并将测试结果做出分析。
6. test分两部分，一部分是旧的代码会test保证输出一致，另一部分是你增加的方法需要自己写个test文件，增加test，这样后面贡献的人要改你代码时也会test，确保不会更变你的程序功能。`TransBigData`的测试结果在[![codecov](https://codecov.io/gh/ni1o1/transbigdata/branch/main/graph/badge.svg?token=GLAVYYCD9L)](https://codecov.io/gh/ni1o1/transbigdata)这里可以看到，其中的百分比表示单元测试覆盖率，表明有多少比例的代码通过了测试。
7. 测试成功后，我们将检查您的更改，并可能要求您在最终准备合并之前进行其他更改。如果成功，我们将merge到`main`分支中，贡献就成功啦。

## 如何贡献代码的视频介绍
[bilibili](https://www.bilibili.com/video/BV1K44y1H7ML/)
[Youtube](https://www.youtube.com/watch?v=ocjzT-23pak)
---
title: 'TransBigData: A Python package for transportation spatio-temporal big data processing, analysis and visualization'
tags:
  - Python
  - transportation
  - spatio-temporal data
  - geospatial data
  - GIS
  - data quality analysis
  - data pre-processing
  - data visualization
  - taxi GPS data
  - bus GPS data
  - bike sharing data
authors:
  - name: Qing Yu
    orcid: 0000-0003-2513-2969
    affiliation: 1
  - name: Jian Yuan^[corresponding author]
    orcid: 0000-0002-7202-0946
    affiliation: 1
affiliations:
 - name: Key Laboratory of Road and Traffic Engineering of the Ministry of Education, Tongji University, 4800 Cao’an Road, Shanghai 201804, People’s Republic of China
   index: 1
date: 30 November 2021
bibliography: paper.bib
---
# Summary

In recent years, data generated in the field of transportation has begun to explode. Individual continuous tracking data, such as mobile phone data, IC smart card data, taxi GPS data, bus GPS data and bicycle sharing order data, also known as "spatio-temporal big data" or "Track &Trace data" [@harrison:2020], has great potential for applications in data-driven transportation research. These spatio-temporal big data typically require three aspects of information [@zhang:2021]: Who? When? Where? They are characterized by high data quality, large collection scope, and fine-grained spatio-temporal information, which can fully capture the daily activities of individuals and their travel behavior in the city in both temporal and spatial dimensions. The emergence of these data provides new ways and opportunities for potential transportation demand analysis and travel mechanism understanding in supporting urban transportation planning and management [@chen:2021; @zhang:2020]. However, processing with these multi-source spatio-temporal big data usually requires a series of similar processing procedure (e.g., data quality assessment, data preprocessing, data cleaning, data gridding, data aggregation, and data visualization). There is an urgent need for a one-size-fits-all tool that can adapt to the various processing demands of different transportation data in this field.

# **State of the art**

Typical processing for spatiotemporal data analysis involves multiple procedures, including data acquisition, data preprocessing, data analysis, data visualization, etc. Currently, there exists several open source packages in these domains:

- Data acquisition: In this aspect, existing software mainly provides data acquisition of basic geospatial data. `osmnx` [@boeing2017osmnx] is a Python package to download geospatial data from OpenStreetMap. Some packages provide tools to generate synthetic mobility data using standard mathematical models. `scikit-mobility` [@scikitmobility] is a library that allows for managing and generation of mobility data of various formats.
- Data preprocessing and data analysis: Spatio-temporal data in the field of transportation involves multiple source of data. Most existing tools are developed for specific type of data, such as trajectory data, air traffic data, etc. `MovingPandas` [@graser2019movingpandas] provides trajectory data structures and functions for data exploration and visualization. `PTRAIL` [@haidri2021ptrail] is a library for mobility data preprocessing especially in feature generation and trajectory interpolation. `traffic` [@Olive2019] is a toolbox for processing and analyzing air traffic data. `trackintel` [@axhausen2007definition] is a framework for spatio-temporal analysis of tracking data focusing on human mobility. `PySAL` [@pysal2007; @Rey2021] is a family of packages that allows for advanced geospatial data science, which supports the development of high-level applications.
- Data visualization: Apart from data processing, `MovingPandas` also support trajectory visualizations. `moveVis` [@2020moveVis] provides tools to visualize movement data and temporal changes of environmental data by creating video animations. `TraViA` [@Siebinga2021] provides tools to visualize and annotate movement data in an interactive approach.

The above-mentioned libraries provide preprocessing and geometric analysis functions from sepcific aspects. However, much remains to be done in dealing with the task of transforming raw spatio-temporal data into valuable insights, and these libraries provide no single solution. Thus, a library compatible with existing tools that provides and analysis framework for transportation spatio-temporal big data can effectively facilitate the research progress in this field.

# Statement of need

`TransBigData` is a Python package developed for transportation spatio-temporal big data processing, analysis, and visualization that provides fast and concise methods for processing taxi GPS data, bicycle sharing data, and bus GPS data, for example. Further, a variety of processing methods for each stage of transportation spatio-temporal big data analysis are included in the code base. `TransBigData` provides clean, efficient, flexible, and easy to use API, allowing complex data tasks to be achieved with concise code, and it has already been used in a number of scientific publications [@yu:2020-1; @yu:2020; @li2021taxi; @yu2021partitioning].

For some types of data, `TransBigData` also provides targeted tools for specific needs, such as extraction of origins and destinations (OD) of taxi trips from taxi GPS data and identification of arrival and departure information from bus GPS data.

Currently, `TransBigData` mainly provides the following methods:

- *Data Quality*: Provides methods to quickly obtain the general information of the dataset, including the data amount, the time period, and the sampling interval.
- *Data Preprocess*: Provides methods to fix multiple types of data error.
- *Data Gridding*: Provides methods to generate multiple types of geographic grids (rectangular and hexagonal) in the research area. Provides fast algorithms to map GPS data to the generated grids (\autoref{fig:fig1}).
- *Data Aggregating*: Provides methods to aggregate GPS and OD data into geographic polygons.
- *Trajectory Processing*: Provides quick methods to re-organize the data structure and implement data augmentation from various data formats, including generating trajectory linestrings from GPS points, and trajectory densification, etc.
- *Data Visualization*: Built-in visualization capabilities leverage the visualization package `keplergl` to interactively visualize data in Jupyter notebooks with simple code.
- *Basemap Loading*: Provides methods to display Mapbox basemaps in `matplotlib` figures (\autoref{fig:fig2}).

The target audience of `TransBigData` includes: 1) Data science researchers and data engineers in the field of transportation big data, smart transportation systems, and urban computing, particularly those who want to integrate innovative algorithms into intelligent trasnportation systems; 2) Government, enterprises, or other entities who expect efficient and reliable management decision support through transportation spatio-temporal data analysis.

The latest stable release of the software can be installed via `pip` and full documentation
can be found at https://transbigdata.readthedocs.io/en/latest/.

![TransBigData</code></code></code></code></code></code></code></code> generates rectangular grids and aggregates GPS data to the grids.\label{fig:fig1}](images/figure1.png){ width=100% }

![TransBigData</code></code></code></code></code></code></code></code> visualizes taxi trip ODs and displays basemaps with matplotlib</code></code></code></code></code></code></code></code>.\label{fig:fig2}](images/figure2.png){ width=100% }

# References
English [中文版](README-zh_CN.md)

# TransBigData

<img src="https://github.com/ni1o1/transbigdata/raw/main/docs/source/_static/logo-wordmark-dark.png" style="width:550px">

[![Documentation Status](https://readthedocs.org/projects/transbigdata/badge/?version=latest)](https://transbigdata.readthedocs.io/en/latest/?badge=latest) [![Downloads](https://pepy.tech/badge/transbigdata)](https://pepy.tech/project/transbigdata) [![Conda Downloads](https://img.shields.io/conda/dn/conda-forge/transbigdata.svg)](https://anaconda.org/conda-forge/transbigdata)  [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/ni1o1/transbigdata/d7d6fa33ff16440ba1698b10dd3cf3f76ff00abd?urlpath=lab%2Ftree%2Fexample%2FExample%201-Taxi%20GPS%20data%20processing.ipynb) [![Tests](https://github.com/ni1o1/transbigdata/actions/workflows/tests.yml/badge.svg)](https://github.com/ni1o1/transbigdata/actions/workflows/tests.yml) [![codecov](https://codecov.io/gh/ni1o1/transbigdata/branch/main/graph/badge.svg?token=GLAVYYCD9L)](https://codecov.io/gh/ni1o1/transbigdata)

## Introduction

`TransBigData` is a Python package developed for transportation spatio-temporal big data processing, analysis and visualization. `TransBigData` provides fast and concise methods for processing common transportation spatio-temporal big data such as Taxi GPS data, bicycle sharing data and bus GPS data. `TransBigData` provides a variety of processing methods for each stage of transportation spatio-temporal big data analysis. The code with `TransBigData` is clean, efficient, flexible, and easy to use, allowing complex data tasks to be achieved with concise code.

For some specific types of data, `TransBigData` also provides targeted tools for specific needs, such as extraction of Origin and Destination(OD) of taxi trips from taxi GPS data and identification of arrival and departure information from bus GPS data. The latest stable release of the software can be installed via pip and full documentation
can be found at https://transbigdata.readthedocs.io/en/latest/.

### Target Audience

The target audience of `TransBigData` includes:

- Data science researchers and data engineers in the field of transportation big data, smart transportation systems, and urban computing, particularly those who want to integrate innovative algorithms into intelligent trasnportation systems
- Government, enterprises, or other entities who expect efficient and reliable management decision support through transportation spatio-temporal data analysis.

### Technical Features

* Provide a variety of processing methods for each stage of transportation spatio-temporal big data analysis.
* The code with `TransBigData` is clean, efficient, flexible, and easy to use, allowing complex data tasks to be achieved with concise code.

### Main Functions

Currently, `TransBigData` mainly provides the following methods:

* **Data Quality**: Provides methods to quickly obtain the general information of the dataset, including the data amount the time period and the sampling interval.
* **Data Preprocess**: Provides methods to clean multiple types of data error.
* **Data Gridding**: Provides methods to generate multiple types of geographic grids (Rectangular grids, Hexagonal grids) in the research area. Provides fast algorithms to map GPS data to the generated grids.
* **Data Aggregating**: Provides methods to aggregate GPS data and OD data into geographic polygon.
* **Data Visualization**: Built-in visualization capabilities leverage the visualization package keplergl to interactively visualize data on Jupyter notebook with simple code.
* **Trajectory Processing**: Provides methods to process trajectory data, including generating trajectory linestring from GPS points, and trajectory densification, etc.
* **Basemap Loading**: Provides methods to display Mapbox basemap on matplotlib figures

## Installation

It is recommended to use `Python 3.7, 3.8, 3.9`

### Using pypi [![PyPI version](https://badge.fury.io/py/transbigdata.svg)](https://badge.fury.io/py/transbigdata)

`TransBigData` can be installed by using `pip install`. Before installing `TransBigData`, make sure that you have installed the available [geopandas package](https://geopandas.org/en/stable/getting_started/install.html). If you already have geopandas installed, run the following code directly from the command prompt to install `TransBigData`:

    pip install transbigdata

### Using conda-forge [![Conda Version](https://img.shields.io/conda/vn/conda-forge/transbigdata.svg)](https://anaconda.org/conda-forge/transbigdata)

You can also install `TransBigData` by `conda-forge`, this will automaticaly solve the dependency, it can be installed with:

    conda install -c conda-forge transbigdata

## Contributing to TransBigData [![GitHub contributors](https://img.shields.io/github/contributors/ni1o1/transbigdata.svg)](https://github.com/ni1o1/transbigdata/graphs/contributors) [![Join the chat at https://gitter.im/transbigdata/community](https://badges.gitter.im/transbigdata/community.svg)](https://gitter.im/transbigdata/community?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge) ![GitHub commit activity](https://img.shields.io/github/commit-activity/m/ni1o1/transbigdata)

All contributions, bug reports, bug fixes, documentation improvements, enhancements and ideas are welcome. A detailed overview on how to contribute can be found in the [contributing guide](https://github.com/ni1o1/transbigdata/blob/master/CONTRIBUTING.md) on GitHub.

## Examples

### Example of data visualization

#### Visualize trajectories (with keplergl)

![gif](https://github.com/ni1o1/transbigdata/raw/main/images/tbdexample1.gif)

#### Visualize data distribution (with keplergl)

![gif](https://github.com/ni1o1/transbigdata/raw/main/images/tbdexample2.gif)

#### Visualize OD (with keplergl)

![gif](https://github.com/ni1o1/transbigdata/raw/main/images/tbdexample3.gif)

### Example of taxi GPS data processing

The following example shows how to use the `TransBigData` to perform data gridding, data aggregating and data visualization for taxi GPS data.

#### Read the data

```python
import transbigdata as tbd
import pandas as pd
#Read taxi gps data  
data = pd.read_csv('TaxiData-Sample.csv',header = None) 
data.columns = ['VehicleNum','time','lon','lat','OpenStatus','Speed'] 
data
```

<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>VehicleNum</th>
      <th>time</th>
      <th>lon</th>
      <th>lat</th>
      <th>OpenStatus</th>
      <th>Speed</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>34745</td>
      <td>20:27:43</td>
      <td>113.806847</td>
      <td>22.623249</td>
      <td>1</td>
      <td>27</td>
    </tr>
    <tr>
      <th>1</th>
      <td>34745</td>
      <td>20:24:07</td>
      <td>113.809898</td>
      <td>22.627399</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>2</th>
      <td>34745</td>
      <td>20:24:27</td>
      <td>113.809898</td>
      <td>22.627399</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>3</th>
      <td>34745</td>
      <td>20:22:07</td>
      <td>113.811348</td>
      <td>22.628067</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>4</th>
      <td>34745</td>
      <td>20:10:06</td>
      <td>113.819885</td>
      <td>22.647800</td>
      <td>0</td>
      <td>54</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>544994</th>
      <td>28265</td>
      <td>21:35:13</td>
      <td>114.321503</td>
      <td>22.709499</td>
      <td>0</td>
      <td>18</td>
    </tr>
    <tr>
      <th>544995</th>
      <td>28265</td>
      <td>09:08:02</td>
      <td>114.322701</td>
      <td>22.681700</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>544996</th>
      <td>28265</td>
      <td>09:14:31</td>
      <td>114.336700</td>
      <td>22.690100</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>544997</th>
      <td>28265</td>
      <td>21:19:12</td>
      <td>114.352600</td>
      <td>22.728399</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>544998</th>
      <td>28265</td>
      <td>19:08:06</td>
      <td>114.137703</td>
      <td>22.621700</td>
      <td>0</td>
      <td>0</td>
    </tr>
  </tbody>
</table>
<p>544999 rows × 6 columns</p>
</div>

#### Data pre-processing

Define the study area and use the `tbd.clean_outofbounds` method to delete the data out of the study area

```python
#Define the study area
bounds = [113.75, 22.4, 114.62, 22.86]
#Delete the data out of the study area
data = tbd.clean_outofbounds(data,bounds = bounds,col = ['lon','lat'])
```

#### Data gridding

The most basic way to express the data distribution is in the form of geograpic grids. `TransBigData` provides methods to generate multiple types of geographic grids (Rectangular grids, Hexagonal grids) in the research area. For rectangular gridding, you need to determine the gridding parameters at first (which can be interpreted as defining a grid coordinate system):

```python
#Obtain the gridding parameters
params = tbd.grid_params(bounds,accuracy = 1000)
```

The next step is to map the GPS data to their corresponding grids. Using the `tbd.GPS_to_grids`, it will generate the `LONCOL` column and the `LATCOL` column. The two columns together can specify a grid:

```python
#Map the GPS data to grids
data['LONCOL'],data['LATCOL'] = tbd.GPS_to_grids(data['lon'],data['lat'],params)
```

Count the amount of data in each grids, generate the geometry of the grids and transform it into a GeoDataFrame:

```python
#Aggregate data into grids
grid_agg = data.groupby(['LONCOL','LATCOL'])['VehicleNum'].count().reset_index()
#Generate grid geometry
grid_agg['geometry'] = tbd.gridid_to_polygon(grid_agg['LONCOL'],grid_agg['LATCOL'],params)
#Change the type into GeoDataFrame
import geopandas as gpd
grid_agg = gpd.GeoDataFrame(grid_agg)
#Plot the grids
grid_agg.plot(column = 'VehicleNum',cmap = 'autumn_r')
```

![png](https://github.com/ni1o1/transbigdata/raw/main/images/output_5_1.png)

#### Data Visualization(with basemap)

For a geographical data visualization figure, we still have to add the basemap, the colorbar, the compass and the scale. Use `tbd.plot_map` to load the basemap and `tbd.plotscale` to add compass and scale in matplotlib figure:

```python
import matplotlib.pyplot as plt
fig =plt.figure(1,(8,8),dpi=300)
ax =plt.subplot(111)
plt.sca(ax)
#Load basemap
tbd.plot_map(plt,bounds,zoom = 11,style = 4)
#Define colorbar
cax = plt.axes([0.05, 0.33, 0.02, 0.3])
plt.title('Data count')
plt.sca(ax)
#Plot the data
grid_agg.plot(column = 'VehicleNum',cmap = 'autumn_r',ax = ax,cax = cax,legend = True)
#Add scale
tbd.plotscale(ax,bounds = bounds,textsize = 10,compasssize = 1,accuracy = 2000,rect = [0.06,0.03],zorder = 10)
plt.axis('off')
plt.xlim(bounds[0],bounds[2])
plt.ylim(bounds[1],bounds[3])
plt.show()
```

![png](https://github.com/ni1o1/transbigdata/raw/main/images/output_7_0.png)

## Citation information [![DOI](https://zenodo.org/badge/419559811.svg)](https://zenodo.org/badge/latestdoi/419559811) [![status](https://joss.theoj.org/papers/d1055fe3105dfa2dcff4cb6c7688a79b/status.svg)](https://joss.theoj.org/papers/d1055fe3105dfa2dcff4cb6c7688a79b)

Please cite [this](https://doi.org/10.21105/joss.04021) when using `TransBigData` in your research. Citation information can be found at [CITATION.cff](https://github.com/ni1o1/transbigdata/blob/main/CITATION.cff).

## Introducing Video (In Chinese) [![bilibili](https://img.shields.io/badge/bilibili-%E5%90%8C%E6%B5%8E%E5%B0%8F%E6%97%AD%E5%AD%A6%E9%95%BF-green.svg)](https://space.bilibili.com/3051484)

* [Bilibili](https://www.bilibili.com/video/BV1nq4y1u7i1)
* [Youtube](https://www.youtube.com/watch?v=V_KHFv75W_w)
---
name: Bug report
about: Create a report to help us improve
title: ''
labels: ''
assignees: ''

---

**Describe the bug**
A clear and concise description of what the bug is.

**To Reproduce**
Steps to reproduce the behavior:
1. Go to '...'
2. Click on '....'
3. Scroll down to '....'
4. See error

**Expected behavior**
A clear and concise description of what you expected to happen.

**Screenshots**
If applicable, add screenshots to help explain your problem.

**Desktop (please complete the following information):**
 - OS: [e.g. iOS]
 - Browser [e.g. chrome, safari]
 - Version [e.g. 22]

**Smartphone (please complete the following information):**
 - Device: [e.g. iPhone6]
 - OS: [e.g. iOS8.1]
 - Browser [e.g. stock browser, safari]
 - Version [e.g. 22]

**Additional context**
Add any other context about the problem here.
---
name: Feature request
about: Suggest an idea for this project
title: ''
labels: ''
assignees: ''

---

**Is your feature request related to a problem? Please describe.**
A clear and concise description of what the problem is. Ex. I'm always frustrated when [...]

**Describe the solution you'd like**
A clear and concise description of what you want to happen.

**Describe alternatives you've considered**
A clear and concise description of any alternative solutions or features you've considered.

**Additional context**
Add any other context or screenshots about the feature request here.
.. _gisprocess:


******************************
GIS处理
******************************

近邻匹配
================

| 下面的案例展示如何用TransBigData包进行点与点、点与线的近邻匹配。该方法使用的是KDTree算法，可查看wiki：https://en.wikipedia.org/wiki/K-d_tree，算法复杂度为o(log(n))


点与点匹配（DataFrame与DataFrame）
----------------------------------

| 导入TransBigData包

::

    import transbigdata as tbd

生成两个GeoDataFrame表，但它们只有经纬度列

::

    import pandas as pd
    import geopandas as gpd
    from shapely.geometry import LineString
    dfA = gpd.GeoDataFrame([[1,2],[2,4],[2,6],
                            [2,10],[24,6],[21,6],
                            [22,6]],columns = ['lon1','lat1'])
    dfB = gpd.GeoDataFrame([[1,3],[2,5],[2,2]],columns = ['lon','lat'])

使用tbd.ckdnearest进行点与点匹配，如果是DataFrame与DataFrame匹配（不含有地理信息），则需要指定前后两个表的经纬度列


.. function:: transbigdata.ckdnearest(dfA_origin,dfB_origin,Aname = ['lon','lat'],Bname = ['lon','lat'])

输入两个DataFrame，分别指定经纬度列名，为表A匹配表B中最近点，并计算距离

**输入**

dfA_origin : DataFrame
    表A
dfB_origin : DataFrame
    表B
Aname : List
    表A中经纬度列字段
Bname : List
    表B中经纬度列字段

**输出**

gdf : DataFrame
    为A匹配到B上最近点的表

::

    tbd.ckdnearest(dfA,dfB,Aname=['lon1','lat1'],Bname=['lon','lat'])
    #此时计算出的距离为经纬度换算实际距离




.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>lon1</th>
          <th>lat1</th>
          <th>index</th>
          <th>lon</th>
          <th>lat</th>
          <th>dist</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>0</th>
          <td>1</td>
          <td>2</td>
          <td>0</td>
          <td>1</td>
          <td>3</td>
          <td>1.111949e+05</td>
        </tr>
        <tr>
          <th>1</th>
          <td>2</td>
          <td>4</td>
          <td>1</td>
          <td>2</td>
          <td>5</td>
          <td>1.111949e+05</td>
        </tr>
        <tr>
          <th>2</th>
          <td>2</td>
          <td>6</td>
          <td>1</td>
          <td>2</td>
          <td>5</td>
          <td>1.111949e+05</td>
        </tr>
        <tr>
          <th>3</th>
          <td>2</td>
          <td>10</td>
          <td>1</td>
          <td>2</td>
          <td>5</td>
          <td>5.559746e+05</td>
        </tr>
        <tr>
          <th>4</th>
          <td>24</td>
          <td>6</td>
          <td>1</td>
          <td>2</td>
          <td>5</td>
          <td>2.437393e+06</td>
        </tr>
        <tr>
          <th>5</th>
          <td>21</td>
          <td>6</td>
          <td>1</td>
          <td>2</td>
          <td>5</td>
          <td>2.105798e+06</td>
        </tr>
        <tr>
          <th>6</th>
          <td>22</td>
          <td>6</td>
          <td>1</td>
          <td>2</td>
          <td>5</td>
          <td>2.216318e+06</td>
        </tr>
      </tbody>
    </table>
    </div>



点与点匹配（GeoDataFrame与GeoDataFrame）
----------------------------------------

将A表B表变为含有点信息的GeoDataFrame

::

    dfA['geometry'] = gpd.points_from_xy(dfA['lon1'],dfA['lat1'])
    dfB['geometry'] = gpd.points_from_xy(dfB['lon'],dfB['lat'])

使用tbd.ckdnearest_point进行点与点匹配

.. function:: transbigdata.ckdnearest_point(gdA, gdB)

输入两个GeoDataFrame，gdfA、gdfB均为点，该方法会为gdfA表连接上gdfB中最近的点，并添加距离字段dist

**输入**

gdA : GeoDataFrame
    表A，点要素
gdB : GeoDataFrame
    表B，点要素

**输出**

gdf : GeoDataFrame
    为A匹配到B上最近点的表


::

    tbd.ckdnearest_point(dfA,dfB)
    #此时计算出的距离为经纬度距离




.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>lon1</th>
          <th>lat1</th>
          <th>geometry_x</th>
          <th>dist</th>
          <th>index</th>
          <th>lon</th>
          <th>lat</th>
          <th>geometry_y</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>0</th>
          <td>1</td>
          <td>2</td>
          <td>POINT (1.00000 2.00000)</td>
          <td>1.000000</td>
          <td>0</td>
          <td>1</td>
          <td>3</td>
          <td>POINT (1.00000 3.00000)</td>
        </tr>
        <tr>
          <th>1</th>
          <td>2</td>
          <td>4</td>
          <td>POINT (2.00000 4.00000)</td>
          <td>1.000000</td>
          <td>1</td>
          <td>2</td>
          <td>5</td>
          <td>POINT (2.00000 5.00000)</td>
        </tr>
        <tr>
          <th>2</th>
          <td>2</td>
          <td>6</td>
          <td>POINT (2.00000 6.00000)</td>
          <td>1.000000</td>
          <td>1</td>
          <td>2</td>
          <td>5</td>
          <td>POINT (2.00000 5.00000)</td>
        </tr>
        <tr>
          <th>3</th>
          <td>2</td>
          <td>10</td>
          <td>POINT (2.00000 10.00000)</td>
          <td>5.000000</td>
          <td>1</td>
          <td>2</td>
          <td>5</td>
          <td>POINT (2.00000 5.00000)</td>
        </tr>
        <tr>
          <th>4</th>
          <td>24</td>
          <td>6</td>
          <td>POINT (24.00000 6.00000)</td>
          <td>22.022716</td>
          <td>1</td>
          <td>2</td>
          <td>5</td>
          <td>POINT (2.00000 5.00000)</td>
        </tr>
        <tr>
          <th>5</th>
          <td>21</td>
          <td>6</td>
          <td>POINT (21.00000 6.00000)</td>
          <td>19.026298</td>
          <td>1</td>
          <td>2</td>
          <td>5</td>
          <td>POINT (2.00000 5.00000)</td>
        </tr>
        <tr>
          <th>6</th>
          <td>22</td>
          <td>6</td>
          <td>POINT (22.00000 6.00000)</td>
          <td>20.024984</td>
          <td>1</td>
          <td>2</td>
          <td>5</td>
          <td>POINT (2.00000 5.00000)</td>
        </tr>
      </tbody>
    </table>
    </div>



点与线匹配（GeoDataFrame与GeoDataFrame）
----------------------------------------

将A表变为地理点，B表为线

::

    dfA['geometry'] = gpd.points_from_xy(dfA['lon1'],dfA['lat1'])
    dfB['geometry'] = [LineString([[1,1],[1.5,2.5],[3.2,4]]),
                      LineString([[1,0],[1.5,0],[4,0]]),
                       LineString([[1,-1],[1.5,-2],[4,-4]])]
    dfB.plot()


.. image:: _static/output_15_1.png



.. function:: transbigdata.ckdnearest_line(gdfA, gdfB)

输入两个GeoDataFrame，其中gdfA为点，gdfB为线，该方法会为gdfA表连接上gdfB中最近的线，并添加距离字段dist

**输入**

gdA : GeoDataFrame
    表A，点要素
gdB : GeoDataFrame
    表B，线要素

**输出**

gdf : GeoDataFrame
    为A匹配到B中最近的线

用tbd.ckdnearest_line可以实现点匹配线，其原理是将线中的折点提取，然后使用点匹配点。

::

    tbd.ckdnearest_line(dfA,dfB)
    #此时计算出的距离为经纬度距离




.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>lon1</th>
          <th>lat1</th>
          <th>geometry_x</th>
          <th>dist</th>
          <th>index</th>
          <th>lon</th>
          <th>lat</th>
          <th>geometry_y</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>0</th>
          <td>1</td>
          <td>2</td>
          <td>POINT (1.00000 2.00000)</td>
          <td>0.707107</td>
          <td>0</td>
          <td>1</td>
          <td>3</td>
          <td>LINESTRING (1.00000 1.00000, 1.50000 2.50000, ...</td>
        </tr>
        <tr>
          <th>1</th>
          <td>2</td>
          <td>4</td>
          <td>POINT (2.00000 4.00000)</td>
          <td>1.200000</td>
          <td>0</td>
          <td>1</td>
          <td>3</td>
          <td>LINESTRING (1.00000 1.00000, 1.50000 2.50000, ...</td>
        </tr>
        <tr>
          <th>2</th>
          <td>2</td>
          <td>6</td>
          <td>POINT (2.00000 6.00000)</td>
          <td>2.332381</td>
          <td>0</td>
          <td>1</td>
          <td>3</td>
          <td>LINESTRING (1.00000 1.00000, 1.50000 2.50000, ...</td>
        </tr>
        <tr>
          <th>3</th>
          <td>2</td>
          <td>10</td>
          <td>POINT (2.00000 10.00000)</td>
          <td>6.118823</td>
          <td>0</td>
          <td>1</td>
          <td>3</td>
          <td>LINESTRING (1.00000 1.00000, 1.50000 2.50000, ...</td>
        </tr>
        <tr>
          <th>4</th>
          <td>21</td>
          <td>6</td>
          <td>POINT (21.00000 6.00000)</td>
          <td>17.912007</td>
          <td>0</td>
          <td>1</td>
          <td>3</td>
          <td>LINESTRING (1.00000 1.00000, 1.50000 2.50000, ...</td>
        </tr>
        <tr>
          <th>5</th>
          <td>22</td>
          <td>6</td>
          <td>POINT (22.00000 6.00000)</td>
          <td>18.906084</td>
          <td>0</td>
          <td>1</td>
          <td>3</td>
          <td>LINESTRING (1.00000 1.00000, 1.50000 2.50000, ...</td>
        </tr>
        <tr>
          <th>6</th>
          <td>24</td>
          <td>6</td>
          <td>POINT (24.00000 6.00000)</td>
          <td>20.880613</td>
          <td>1</td>
          <td>2</td>
          <td>5</td>
          <td>LINESTRING (1.00000 0.00000, 1.50000 0.00000, ...</td>
        </tr>
      </tbody>
    </table>
    </div>






打断线
===============

在实际应用中，我们可能会需要把很长的线打断为很多子线段，每一条线段不要超过一定的最大长度，此时则可以使用TransBigData包中的splitline_with_length方法。


.. function:: transbigdata.splitline_with_length(Centerline,maxlength = 100)

输入线GeoDataFrame要素，打断为最大长度maxlength的小线段

**输入**

Centerline : GeoDataFrame
    线要素
maxlength : number
    打断的线段最大长度

**输出**

splitedline : GeoDataFrame
    打断后的线

下面演示如何将线打断为100米一段的线段

::

    #读取线要素
    import geopandas as gpd
    Centerline = gpd.read_file(r'test_lines.json')
    Centerline.plot()





.. image:: splitline/output_2_1.png


::

    #转换线为投影坐标系
    Centerline.crs = {'init':'epsg:4326'}
    Centerline = Centerline.to_crs(epsg = '4517')
    #计算线的长度
    Centerline['length'] = Centerline.length
    Centerline




.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>Id</th>
          <th>geometry</th>
          <th>length</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>0</th>
          <td>0</td>
          <td>LINESTRING (29554925.232 4882800.694, 29554987...</td>
          <td>285.503444</td>
        </tr>
        <tr>
          <th>1</th>
          <td>0</td>
          <td>LINESTRING (29554682.635 4882450.554, 29554773...</td>
          <td>185.482276</td>
        </tr>
        <tr>
          <th>2</th>
          <td>0</td>
          <td>LINESTRING (29554987.079 4882521.969, 29555040...</td>
          <td>291.399180</td>
        </tr>
        <tr>
          <th>3</th>
          <td>0</td>
          <td>LINESTRING (29554987.079 4882521.969, 29555073...</td>
          <td>248.881529</td>
        </tr>
        <tr>
          <th>4</th>
          <td>0</td>
          <td>LINESTRING (29554987.079 4882521.969, 29554969...</td>
          <td>207.571197</td>
        </tr>
        <tr>
          <th>5</th>
          <td>0</td>
          <td>LINESTRING (29554773.177 4882288.671, 29554828...</td>
          <td>406.251357</td>
        </tr>
        <tr>
          <th>6</th>
          <td>0</td>
          <td>LINESTRING (29554773.177 4882288.671, 29554926...</td>
          <td>158.114403</td>
        </tr>
        <tr>
          <th>7</th>
          <td>0</td>
          <td>LINESTRING (29555060.286 4882205.456, 29555082...</td>
          <td>107.426629</td>
        </tr>
        <tr>
          <th>8</th>
          <td>0</td>
          <td>LINESTRING (29555040.278 4882235.468, 29555060...</td>
          <td>36.069941</td>
        </tr>
        <tr>
          <th>9</th>
          <td>0</td>
          <td>LINESTRING (29555060.286 4882205.456, 29555095...</td>
          <td>176.695446</td>
        </tr>
      </tbody>
    </table>
    </div>



::

    #将线打断为最长100米的线段
    import transbigdata as tbd
    splitedline = tbd.splitline_with_length(Centerline,maxlength = 100)

::

    #打断后线型不变
    splitedline.plot()








.. image:: splitline/output_5_1.png


::

    #但内容已经变成一段一段了
    splitedline




.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>geometry</th>
          <th>id</th>
          <th>length</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>0</th>
          <td>LINESTRING (29554925.232 4882800.694, 29554927...</td>
          <td>0</td>
          <td>100.000000</td>
        </tr>
        <tr>
          <th>1</th>
          <td>LINESTRING (29554946.894 4882703.068, 29554949...</td>
          <td>0</td>
          <td>100.000000</td>
        </tr>
        <tr>
          <th>2</th>
          <td>LINESTRING (29554968.557 4882605.443, 29554970...</td>
          <td>0</td>
          <td>85.503444</td>
        </tr>
        <tr>
          <th>0</th>
          <td>LINESTRING (29554682.635 4882450.554, 29554688...</td>
          <td>1</td>
          <td>100.000000</td>
        </tr>
        <tr>
          <th>1</th>
          <td>LINESTRING (29554731.449 4882363.277, 29554736...</td>
          <td>1</td>
          <td>85.482276</td>
        </tr>
        <tr>
          <th>0</th>
          <td>LINESTRING (29554987.079 4882521.969, 29554989...</td>
          <td>2</td>
          <td>100.000000</td>
        </tr>
        <tr>
          <th>1</th>
          <td>LINESTRING (29555005.335 4882423.650, 29555007...</td>
          <td>2</td>
          <td>100.000000</td>
        </tr>
        <tr>
          <th>2</th>
          <td>LINESTRING (29555023.592 4882325.331, 29555025...</td>
          <td>2</td>
          <td>91.399180</td>
        </tr>
        <tr>
          <th>0</th>
          <td>LINESTRING (29554987.079 4882521.969, 29554993...</td>
          <td>3</td>
          <td>100.000000</td>
        </tr>
        <tr>
          <th>1</th>
          <td>LINESTRING (29555042.051 4882438.435, 29555048...</td>
          <td>3</td>
          <td>99.855617</td>
        </tr>
        <tr>
          <th>2</th>
          <td>LINESTRING (29555111.265 4882370.450, 29555116...</td>
          <td>3</td>
          <td>48.881529</td>
        </tr>
        <tr>
          <th>0</th>
          <td>LINESTRING (29554987.079 4882521.969, 29554985...</td>
          <td>4</td>
          <td>100.000000</td>
        </tr>
        <tr>
          <th>1</th>
          <td>LINESTRING (29554973.413 4882422.908, 29554971...</td>
          <td>4</td>
          <td>99.756943</td>
        </tr>
        <tr>
          <th>2</th>
          <td>LINESTRING (29554930.341 4882335.023, 29554929...</td>
          <td>4</td>
          <td>7.571197</td>
        </tr>
        <tr>
          <th>0</th>
          <td>LINESTRING (29554773.177 4882288.671, 29554777...</td>
          <td>5</td>
          <td>100.000000</td>
        </tr>
        <tr>
          <th>1</th>
          <td>LINESTRING (29554816.361 4882198.476, 29554821...</td>
          <td>5</td>
          <td>99.782969</td>
        </tr>
        <tr>
          <th>2</th>
          <td>LINESTRING (29554882.199 4882125.314, 29554891...</td>
          <td>5</td>
          <td>99.745378</td>
        </tr>
        <tr>
          <th>3</th>
          <td>LINESTRING (29554976.612 4882096.588, 29554987...</td>
          <td>5</td>
          <td>100.000000</td>
        </tr>
        <tr>
          <th>4</th>
          <td>LINESTRING (29555076.548 4882100.189, 29555077...</td>
          <td>5</td>
          <td>6.251357</td>
        </tr>
        <tr>
          <th>0</th>
          <td>LINESTRING (29554773.177 4882288.671, 29554783...</td>
          <td>6</td>
          <td>100.000000</td>
        </tr>
        <tr>
          <th>1</th>
          <td>LINESTRING (29554869.914 4882314.006, 29554876...</td>
          <td>6</td>
          <td>58.114403</td>
        </tr>
        <tr>
          <th>0</th>
          <td>LINESTRING (29555060.286 4882205.456, 29555062...</td>
          <td>7</td>
          <td>100.000000</td>
        </tr>
        <tr>
          <th>1</th>
          <td>LINESTRING (29555081.239 4882107.675, 29555081...</td>
          <td>7</td>
          <td>7.426629</td>
        </tr>
        <tr>
          <th>0</th>
          <td>LINESTRING (29555040.278 4882235.468, 29555042...</td>
          <td>8</td>
          <td>36.069941</td>
        </tr>
        <tr>
          <th>0</th>
          <td>LINESTRING (29555060.286 4882205.456, 29555064...</td>
          <td>9</td>
          <td>100.000000</td>
        </tr>
        <tr>
          <th>1</th>
          <td>LINESTRING (29555094.981 4882299.244, 29555100...</td>
          <td>9</td>
          <td>76.419694</td>
        </tr>
      </tbody>
    </table>
    </div>

面要素处理
========================

面合并
-------------------------

.. function:: transbigdata.merge_polygon(data,col)

输入多边形GeoDataFrame数据，以及分组列名col，对不同组别进行分组的多边形进行合并

**输入**

data : GeoDataFrame
    多边形数据
col : str
    分组列名

**输出**

data1 : GeoDataFrame
    合并后的面


对面取外边界构成新多边形
-------------------------


.. function:: transbigdata.polyon_exterior(data,minarea = 0)

输入多边形GeoDataFrame数据，对多边形取外边界构成新多边形

**输入**

data : GeoDataFrame
    多边形数据
minarea : number
    最小面积，小于这个面积的面全部剔除
    
**输出**

data1 : GeoDataFrame
    处理后的面


置信椭圆
========================


置信椭圆参数估计
-------------------------

.. function:: transbigdata.ellipse_params(data,col = ['lon','lat'],confidence = 95,epsg = None)

输入点数据，获取置信椭圆的参数

**输入**

data : DataFrame
    点数据
confidence : number
    置信度，可选99，95，90
epsg : number
    如果给了，则将原始坐标从wgs84，转换至给定epsg坐标系下进行置信椭圆参数估计
col: List
    以[经度，纬度]形式存储的列名

**输出**

params: List
    质心椭圆参数，分别为[pos,width,height,theta,area,alpha]
    对应[中心点坐标，短轴，长轴，角度，面积，方向性]


置信椭圆绘制
-------------------------

.. function:: transbigdata.ellipse_plot(ellip_params,ax,**kwargs)

输入置信椭圆的参数，绘制置信椭圆

**输入**

ellip_params : List
    
ax : matplotlib.axes._subplots.AxesSubplot
    画板

用法
-------------------------

::

    import pandas as pd
    import transbigdata as tbd
    import numpy as np
    #生成测试用数据
    data = np.random.uniform(1,10,(100,2))
    data[:,1:] = 0.5*data[:,0:1]+np.random.uniform(-2,2,(100,1))
    data = pd.DataFrame(data,columns = ['x','y'])
    
    #绘制数据分布
    import matplotlib.pyplot as plt
    plt.figure(1,(5,5))
    #绘制数据点
    plt.scatter(data['x'],data['y'],s = 0.5)
    #绘制坐标轴
    plt.plot([-10,10],[0,0],c = 'k')
    plt.plot([0,0],[-10,10],c = 'k')
    plt.xlim(-15,15)
    plt.ylim(-15,15)
    plt.show()



.. image:: gisprocess/output_1_0.png

输入数据与xy坐标所在列名，置信度，估计椭圆参数
分别代表[中心点坐标，短轴，长轴，角度，面积，扁率

::

    
    ellip_params = tbd.ellipse_params(data,confidence=95,col = ['x','y'])
    ellip_params


.. parsed-literal::

    [array([5.78928146, 2.88466235]),
     4.6981983145616875,
     14.04315715927693,
     -58.15524535916836,
     51.8186366184246,
     0.6654457212665993]

再用tbd.ellipse_plot绘制置信椭圆

::

    #绘制数据分布
    import matplotlib.pyplot as plt
    plt.figure(1,(5,5))
    ax = plt.subplot(111)
    #绘制数据点
    plt.scatter(data['x'],data['y'],s = 0.5)
    #获取置信椭圆参数并绘制椭圆
    #99%置信椭圆
    ellip_params = tbd.ellipse_params(data,confidence=99,col = ['x','y'])
    tbd.ellipse_plot(ellip_params,ax,fill = False,edgecolor = 'r',linewidth = 1)
    #95%置信椭圆
    ellip_params = tbd.ellipse_params(data,confidence=95,col = ['x','y'])
    tbd.ellipse_plot(ellip_params,ax,fill = False,edgecolor = 'b',linewidth = 1)
    #90%置信椭圆
    ellip_params = tbd.ellipse_params(data,confidence=90,col = ['x','y'])
    tbd.ellipse_plot(ellip_params,ax,fill = False,edgecolor = 'k',linewidth = 1)
    #绘制坐标轴
    plt.plot([-10,10],[0,0],c = 'k')
    plt.plot([0,0],[-10,10],c = 'k')
    plt.xlim(-15,15)
    plt.ylim(-15,15)
    plt.show()



.. image:: gisprocess/output_3_0.png



.. _traj:


******************************
轨迹处理
******************************

停留与出行识别
==================
.. function:: transbigdata.traj_stay_move(data,params,col = ['ID','dataTime','longitude','latitude'],activitytime = 1800)

输入轨迹数据与栅格化参数，识别活动与出行

**输入**

data : DataFrame
    轨迹数据集
params : List
    栅格化参数
col : List
    数据的列名[个体，时间，经度，纬度]顺序
activitytime : Number
    多长时间识别为停留
    
**输出**

stay : DataFrame
    个体停留信息
move : DataFrame
    个体移动信息

.. function:: transbigdata.plot_activity(data,col = ['stime','etime','LONCOL','LATCOL'])

输入个体的活动数据（单一个体），绘制活动图

**输入**

data : DataFrame
    活动数据集
col : List
    列名，分别为[活动开始时间，活动结束时间，活动所在栅格经度编号，活动所在栅格纬度编号]

轨迹线型生成
==================

.. function:: transbigdata.points_to_traj(traj_points,col = ['Lng','Lat','ID'],timecol = None)

输入轨迹点，生成轨迹线型的GeoDataFrame

**输入**

traj_points : DataFrame
    轨迹点数据
col : List
    列名，按[经度,纬度,轨迹编号]的顺序
timecol : str
    可选，时间列的列名，如果给了则输出带有[经度,纬度,高度,时间]的geojson，可放入kepler中可视化轨迹

**输出**

traj : GeoDataFrame或json
    生成的轨迹数据，如果timecol没定义则为GeoDataFrame，否则为json


轨迹增密与稀疏化
==================

.. function:: transbigdata.traj_densify(data,col = ['Vehicleid','Time','Lng','Lat'],timegap = 15)

轨迹点增密，确保每隔timegap秒会有一个轨迹点

**输入**

data : DataFrame
    数据
col : List
    列名，按[车辆ID,时间,经度,纬度]的顺序
timegap : number
    单位为秒，每隔多长时间插入一个轨迹点

**输出**

data1 : DataFrame
    处理后的数据

.. function:: transbigdata.traj_sparsify(data,col = ['Vehicleid','Time','Lng','Lat'],timegap = 15)

轨迹点稀疏化。轨迹数据采样间隔过高的时候，数据量太大，不便于分析。这个函数可以将采样间隔扩大，缩减数据量

**输入**

data : DataFrame
    数据
col : List
    列名，按[车辆ID,时间,经度,纬度]的顺序
timegap : number
    单位为秒，每隔多长时间一个轨迹点
method : str
    可选interpolate插值或subsample子采样

**输出**

data1 : DataFrame
    处理后的数据

使用方法

::

    import transbigdata as tbd
    import pandas as pd
    #读取数据    
    data = pd.read_csv('TaxiData-Sample.csv',header = None) 
    data.columns = ['Vehicleid','Time','Lng','Lat','OpenStatus','Speed']      
    data['Time'] = pd.to_datetime(data['Time'])
    #轨迹增密前的采样间隔
    tbd.data_summary(data,col = ['Vehicleid','Time','Lng','Lat'],show_sample_duration=True)

::

    数据量
    -----------------
    数据总量 : 544999 条
    个体总量 : 180 个
    个体数据量均值 : 3027.77 条
    个体数据量上四分位 : 4056.25 条
    个体数据量中位数 : 2600.5 条
    个体数据量下四分位 : 1595.75 条

    数据时间段
    -----------------
    开始时间 : 2021-11-12 00:00:00
    结束时间 : 2021-11-12 23:59:59

    个体采样间隔
    -----------------
    均值 : 28.0 秒
    上四分位 : 30.0 秒
    中位数 : 20.0 秒
    下四分位 : 15.0 秒

进行轨迹增密，设置15秒一条数据::
    
    data1 = tbd.traj_densify(data,timegap = 15)
    #轨迹增密后的采样间隔
    tbd.data_summary(data1,show_sample_duration=True)

::

    数据量
    -----------------
    数据总量 : 1526524 条
    个体总量 : 180 个
    个体数据量均值 : 8480.69 条
    个体数据量上四分位 : 9554.75 条
    个体数据量中位数 : 8175.0 条
    个体数据量下四分位 : 7193.5 条

    数据时间段
    -----------------
    开始时间 : 2021-11-12 00:00:00
    结束时间 : 2021-11-12 23:59:59

    个体采样间隔
    -----------------
    均值 : 9.99 秒
    上四分位 : 15.0 秒
    中位数 : 11.0 秒
    下四分位 : 6.0 秒

增密后的效果

.. image:: example-taxi/densify.png

::

    #两辆车的数据测试
    tmp = data.iloc[:10]
    tmp1 = data.iloc[-100:]
    tmp = tmp.append(tmp1)

    #增密前数据
    import geopandas as gpd
    tmp['geometry'] = gpd.points_from_xy(tmp['Lng'],tmp['Lat'])
    tmp = gpd.GeoDataFrame(tmp)
    tmp[tmp['Vehicleid']==36805].plot()

    #进行轨迹增密，设置5秒一条数据
    tmp1 = tbd.traj_densify(tmp,timegap = 1)
    import geopandas as gpd
    tmp1['geometry'] = gpd.points_from_xy(tmp1['Lng'],tmp1['Lat'])
    tmp1 = gpd.GeoDataFrame(tmp1)
    tmp1[tmp1['Vehicleid']==36805].plot()

    #轨迹稀疏化，20秒一条数据
    tmp2 = tbd.traj_sparsify(tmp1,timegap = 20)
    import geopandas as gpd
    tmp2['geometry'] = gpd.points_from_xy(tmp2['Lng'],tmp2['Lat'])
    tmp2 = gpd.GeoDataFrame(tmp2)
    tmp2[tmp2['Vehicleid']==36805].plot()

.. image:: example-taxi/sparsify.png.. _plot_map:


***************
底图加载
***************

使用前的设置
=============================

| TransBigData包提供了在matplotlib上绘制地图底图的功能，底图由mapbox提供，坐标系为WGS84。如果你要使用该功能，首先需要点击\ `这个链接 <https://account.mapbox.com/auth/signin/?route-to=%22https://account.mapbox.com/%22>`__\ 注册一个mapbox的账号，mapbox上注册成为开发者，并获取到一个mapbox token。 `这个链接 <https://docs.mapbox.com/help/getting-started/access-tokens/#how-access-tokens-work>`__\ 介绍了mapbox token的作用。
| 如果你已经得到了mapbox token，可以用以下代码为TransBigData设置mapbox token(只需要设置一次，后面重新打开python也不需要再重新设置了)：

::

    import transbigdata as tbd
    #用下面代码设置你的mapboxtoken
    tbd.set_mapboxtoken('pk.eyxxxxxxxxxx.xxxxxxxxx')#必须在里面设置你申请的token，直接复制此行代码无效！

另外还需要设置一个地图底图的存储位置，下一次显示同一个位置时，地图会从本地读取加载。

::

    #设置你的地图底图存储路径
    #如果你是linux或者mac系统，路径是这么写，注意最后有一个反斜杠
    tbd.set_imgsavepath(r'/Users/xxxx/xxxx/')
    #如果是windows系统，路径这么写，最后注意要两个斜杠以防转义
    tbd.set_imgsavepath(r'E:\pythonscript\xxx\\')

设置好后，下次绘制底图时，会在你设置的路径下创建一个tileimg文件夹，底图都放在里面  
尝试一下下面的代码，看看能否成功绘制底图

::

    #定义显示范围范围
    bounds = [113.6,22.4,114.8,22.9]
    #创建图框
    import matplotlib.pyplot as plt
    fig =plt.figure(1,(8,8),dpi=250)
    ax =plt.subplot(111)
    plt.sca(ax)
    #添加地图底图
    tbd.plot_map(plt,bounds,zoom = 11,style = 4)
    #添加比例尺和指北针
    tbd.plotscale(ax,bounds = bounds,textsize = 10,compasssize = 1,accuracy = 2000,rect = [0.06,0.03],zorder = 10)
    plt.axis('off')
    plt.xlim(bounds[0],bounds[2])
    plt.ylim(bounds[1],bounds[3])
    plt.show()

.. image:: plot_map/output_6_0.png



地图底图加载
=============================

TransBigData包的底图绘制功能由plot_map包提供。首先确保你的plot_map包在0.3.3版本以上::

    pip install -U plot-map

.. function:: transbigdata.plot_map(plt,bounds,zoom='auto',style=4,printlog = False,styleid = 'dark')

添加地图底图

**输入**

bounds : List
    底图的绘图边界，[lon1,lat1,lon2,lat2] (WGS84坐标系) 其中，lon1,lat1是左下角坐标，lon2,lat2是右上角坐标 
zoom : number
    底图的放大等级，默认为auto自动选取。越大越精细，加载的时间也就越久，一般单个城市大小的范围，这个参数选取12到16之间 
printlog : bool
    是否显示日志                                                
style : number
    地图底图的样式，可选1-10，对应分别如下（需要plot_map包在0.3.3版本以上）   

底图样式1：streets
----------------------------------------

.. image:: plot_map/1.png


底图样式2：outdoors
----------------------------------------

.. image:: plot_map/2.png


底图样式3：satellite
----------------------------------------

.. image:: plot_map/3.png


底图样式4：light
----------------------------------------

.. image:: plot_map/4.png


底图样式5：dark
----------------------------------------

.. image:: plot_map/5.png


底图样式6：light-ch（中文）
----------------------------------------

.. image:: plot_map/6.png


底图样式7：ice creem
----------------------------------------

.. image:: plot_map/7.png


底图样式8：night
----------------------------------------

.. image:: plot_map/8.png


底图样式9：terrain
----------------------------------------

.. image:: plot_map/9.png


底图样式10：basic blue
----------------------------------------

.. image:: plot_map/10.png

用法
----------------------------------------

::

    #设定显示范围
    bounds = [lon1,lat1,lon2,lat2]  
    tbd.plot_map(plt,bounds,zoom = 12,style = 4)  

指北针和比例尺
=============================

.. function:: transbigdata.plotscale(ax,bounds,textcolor = 'k',textsize = 8,compasssize = 1,accuracy = 'auto',rect=[0.1,0.1],unit = "KM",style = 1,**kwargs)

为底图添加指北针和比例尺

**输入**

bounds : List
    底图的绘图边界，[lon1,lat1,lon2,lat2] (WGS84坐标系) 其中，lon1,lat1是左下角坐标，lon2,lat2是右上角坐标 
textsize : number
    标注文字大小                                                 
compasssize : number
    标注的指北针大小                                             
accuracy : number
    标注比例尺的长度（米）                                         
unit : str
    'KM','km','M','m' 比例尺的单位                               
style : number
    1或2，比例尺样式                                             
rect : List
    比例尺在图中的大致位置，如[0.9,0.9]则在右上角                    


::

    tbd.plotscale(ax,bounds = bounds,textsize = 10,compasssize = 1,accuracy = 2000,rect = [0.06,0.03])  

******************************
公交地铁网络拓扑建模
******************************

.. function:: transbigdata.split_subwayline(line,stop)

用公交/地铁站点对公交/地铁线进行切分，得到断面

**输入**

line : GeoDataFrame
    公交/地铁线路
stop : GeoDataFrame
    公交/地铁站点

**输出**

metro_line_splited : GeoDataFrame
    生成的断面线型


.. function:: transbigdata.metro_network(stop,traveltime = 3,transfertime = 5,nxgraph = True)

输入站点信息，输出网络信息，该方法依赖于NetworkX

**输入**

stop : GeoDataFrame
    公交站点
traveltime : number
    每个轨道断面的出行时长
transfertime : number
    每个轨道换乘的时长
nxgraph : bool
    默认True，如果True则直接输出由NetworkX构建的网络G，如果为False，则输出网络的边edge1,edge2,和节点node
    
**输出**

G : networkx.classes.graph.Graph
    networkx构建的网络G，nxgraph参数为True时只输出这个
edge1 : DataFrame
    轨道断面的边，nxgraph参数为False时输出这个
edge2 : DataFrame
    轨道换乘的边，nxgraph参数为False时输出这个
node : List
    网络节点，nxgraph参数为False时输出这个.. _bikedata:


******************************
共享单车数据处理
******************************

.. function:: transbigdata.bikedata_to_od(data,col = ['BIKE_ID','DATA_TIME','LONGITUDE','LATITUDE','LOCK_STATUS'],startend = None)

输入共享单车订单数据（只在开关锁时产生数据），指定列名，提取其中的骑行与停车信息

**输入**

data : DataFrame
    共享单车订单数据，只在开关锁时产生数据
col : List
    列名，顺序不能变，分别为[单车ID,时间,经度,纬度,锁状态]，例如['BIKE_ID','DATA_TIME','LONGITUDE','LATITUDE','LOCK_STATUS']
startend : List
    传入的为[开始时间,结束时间]，如['2018-08-27 00:00:00','2018-08-28 00:00:00']。
    如传入，则考虑（观测时段开始时到单车第一次出现）与（单车最后一次出现到观测时段结束）的骑行与停车情况。
    
**输出**

move_data : DataFrame
    骑行订单数据
stop_data : DataFrame
    停车数据.. _CoordinatesConverter:


******************************
坐标距离
******************************

火星坐标系互转
=============================

坐标互转方法
--------------------------

TransBigData包提供了GCJ02,BD09,BD09mc,WGS94坐标系互转。

.. function:: transbigdata.gcj02tobd09(lng, lat)

.. function:: transbigdata.bd09togcj02(bd_lon, bd_lat)

.. function:: transbigdata.wgs84togcj02(lng, lat)

.. function:: transbigdata.gcj02towgs84(lng, lat)

.. function:: transbigdata.wgs84tobd09(lon,lat)

.. function:: transbigdata.bd09towgs84(lon,lat)

.. function:: transbigdata.bd09mctobd09(lon,lat)

坐标互转，基于numpy列运算::

  data['Lng'],data['Lat'] = tbd.wgs84tobd09(data['Lng'],data['Lat'])  
  data['Lng'],data['Lat'] = tbd.wgs84togcj02(data['Lng'],data['Lat'])  
  data['Lng'],data['Lat'] = tbd.gcj02tobd09(data['Lng'],data['Lat'])  
  data['Lng'],data['Lat'] = tbd.gcj02towgs84(data['Lng'],data['Lat'])  
  data['Lng'],data['Lat'] = tbd.bd09togcj02(data['Lng'],data['Lat'])  
  data['Lng'],data['Lat'] = tbd.bd09towgs84(data['Lng'],data['Lat'])  
  data['Lng'],data['Lat'] = tbd.bd09mctobd09(data['Lng'],data['Lat']) 

对地理要素转换坐标
--------------------------

.. function:: transbigdata.transform_shape(gdf,method)

输入地理要素的GeoDataFrame，对整体做坐标转换

**输入**

gdf : GeoDataFrame
    地理要素
method : function
    坐标转换函数

**输出**

gdf : GeoDataFrame
    转换后结果


::

    #读取线要素
    import geopandas as gpd
    Centerline = gpd.read_file(r'test_lines.json')
    Centerline.plot()


.. image:: transform_shape/output_0_1.png


::

    #整体进行坐标转换
    import transbigdata as tbd
    Centerline_transformed = tbd.transform_shape(Centerline,tbd.bd09towgs84)
    Centerline_transformed.plot()

.. image:: transform_shape/output_1_1.png


经纬度计算距离
=============================

.. function:: transbigdata.getdistance(lon1, lat1, lon2, lat2)

按经度1，纬度1，经度2，纬度2 （十进制度数）顺序输入起终点经纬度，为DataFrame的列，获取距离（米），基于numpy列运算::
    
  data['distance'] = tbd.getdistance(data['Lng1'],data['Lat1'], data['Lng2'],data['Lat2'])  

.. _preprocess:


******************************
数据预处理
******************************


各类数据通用的预处理
============================

.. function:: transbigdata.clean_same(data,col = ['VehicleNum','Time','Lng','Lat'])

删除信息与前后数据相同的数据以减少数据量
如：某个体连续n条数据除了时间以外其他信息都相同，则可以只保留首末两条数据

**输入**

data : DataFrame
    数据
col : List
    列名，按[个体ID,时间,经度,纬度]的顺序，可以传入更多列。会以时间排序，再判断除了时间以外其他列的信息

**输出**

data1 : DataFrame
    清洗后的数据

.. function:: transbigdata.clean_drift(data,col = ['VehicleNum','Time','Lng','Lat'],speedlimit = 80)

删除漂移数据。条件是，此数据与前后的速度都大于speedlimit，但前后数据之间的速度却小于speedlimit。
传入的数据中时间列如果为datetime格式则计算效率更快

**输入**

data : DataFrame
    数据
col : List
    列名，按[个体ID,时间,经度,纬度]的顺序

**输出**

data1 : DataFrame
    研究范围内的数据


.. function:: transbigdata.clean_outofbounds(data,bounds,col = ['Lng','Lat'])

输入研究范围的左下右上经纬度坐标，剔除超出研究范围的数据

**输入**

data : DataFrame
    数据
bounds : List    
    研究范围的左下右上经纬度坐标，顺序为[lon1,lat1,lon2,lat2]
col : List
    经纬度列名

**输出**

data1 : DataFrame
    研究范围内的数据


.. function:: transbigdata.clean_outofshape(data,shape,col = ['Lng','Lat'],accuracy=500)

输入研究范围的GeoDataFrame，剔除超出研究区域的数据

**输入**

data : DataFrame
    数据
shape : GeoDataFrame    
    研究范围的GeoDataFrame
col : List
    经纬度列名
accuracy : number
    计算原理是先栅格化后剔除，这里定义栅格大小，越小精度越高

**输出**

data1 : DataFrame
    研究范围内的数据

.. function:: transbigdata.id_reindex(data,col,new = False,timegap = None,timecol = None,suffix = '_new',sample = None)

对数据的ID列重新编号

**输入**

data : DataFrame
    数据 
col : str
    要重新编号的ID列名
new : bool
    False，相同ID的新编号相同；True，依据表中的顺序，ID再次出现则编号不同
timegap : number
    如果个体在一段时间内没出现（timegap为时间阈值），则编号为新的个体。此参数与timecol同时设定才有效果。
timecol : str
    时间字段名称，此参数与timegap同时设定才有效果。
suffix : str
    新编号列名的后缀，设置为False时替代原有列名
sample : int
    传入数值，对重新编号的个体进行抽样
    
**输出**

data1 : DataFrame
    重新编号的数据

.. function:: transbigdata.id_reindex_disgap(data,col = ['uid','lon','lat'],disgap=1000,suffix = '_new')

对数据的ID列重新编号，如果相邻两条记录超过距离，则编号为新id

**输入**

data : DataFrame
    数据 
col : str
    要重新编号的ID列名
disgap : number
    如果个体轨迹超过一定距离，则编号为新的个体。
suffix : str
    新编号列名的后缀
    
**输出**

data1 : DataFrame
    重新编号的数据

数据格式转换
==================

.. function:: transbigdata.dumpjson(data,path)

输入json数据，存储为文件。这个方法主要是解决numpy数值型无法兼容json包报错的问题

**输入**

data : json
    要储存的json数据
path : str
    保存的路径


轨迹数据清洗
==================
.. function:: transbigdata.clean_traj(data,col = ['uid','str_time','lon','lat'],tripgap = 1800,disgap = 50000,speedlimit = 80)

轨迹数据清洗组合拳

**输入**

data : DataFrame
    轨迹数据
col : List
    列名，以[个体id,时间,经度,纬度]排列
tripgap : number
    多长的时间视为新的出行
disgap : number
    多长距离视为新的出行
speedlimit : number
    车速限制

**输出**

data1 : DataFrame
    清洗后的数据


出租车数据的预处理
==================

.. function:: transbigdata.clean_taxi_status(data,col = ['VehicleNum','Time','OpenStatus'],timelimit = None)

删除出租车数据中载客状态瞬间变化的记录，这些记录的存在会影响出行订单判断。
判断条件为:如果对同一辆车，上一条记录与下一条记录的载客状态都与本条记录不同，则本条记录应该删去

**输入**

data : DataFrame
    数据
col : List
    列名，按[车辆ID,时间,载客状态]的顺序
timelimit : number
    可选，单位为秒，上一条记录与下一条记录的时间小于该时间阈值才予以删除

**输出**

data1 : DataFrame
    清洗后的数据.. _grids:


***************
数据栅格化
***************

方形栅格渔网的生成与对应
=============================

.. function:: transbigdata.rect_grids(location,accuracy = 500,params='auto')

生成研究范围内的方形栅格

**输入**

location : bounds(List) or shape(GeoDataFrame)
    在哪生成栅格
    如果是生成范围的边界bounds，则内容为[lon1,lat1,lon2,lat2] (WGS84坐标系) 其中，lon1,lat1是左下角坐标，lon2,lat2是右上角坐标 
    如果是面要素，则必须是GeoDataFrame
accuracy : number
    栅格大小（米）
params : List
    栅格参数(lonStart,latStart,deltaLon,deltaLat)，或(lonStart,latStart,deltaLon,deltaLat,theta)，其中，lonStart,latStart分别为栅格左下角坐标，deltaLon,deltaLat为单个栅格的经纬度长宽，theta为栅格的角度，不给则默认为0
    默认值为auto自动生成，当给定栅格参数时，栅格大小将从栅格参数中计算得到                   
    

**输出**

grid : GeoDataFrame
    栅格的GeoDataFrame，其中LONCOL与LATCOL为栅格的编号，HBLON与HBLAT为栅格的中心点坐标 
params : List
    栅格参数(lonStart,latStart,deltaLon,deltaLat)，或(lonStart,latStart,deltaLon,deltaLat,theta)，其中，lonStart,latStart分别为栅格左下角坐标，deltaLon,deltaLat为单个栅格的经纬度长宽，theta为栅格的角度，不给则默认为0


::

    #设定范围
    bounds = [lon1,lat1,lon2,lat2]
    grid,params = tbd.rect_grids(bounds,accuracy = 500)


.. function:: transbigdata.grid_params(bounds,accuracy = 500)

栅格参数获取

**输入**

bounds : List
    生成范围的边界，[lon1,lat1,lon2,lat2] (WGS84坐标系) 其中，lon1,lat1是左下角坐标，lon2,lat2是右上角坐标 
accuracy : number
    栅格大小（米）
                                           

**输出**

params : List
    栅格参数(lonStart,latStart,deltaLon,deltaLat)，或(lonStart,latStart,deltaLon,deltaLat,theta)，其中，lonStart,latStart分别为栅格左下角坐标，deltaLon,deltaLat为单个栅格的经纬度长宽，theta为栅格的角度，不给则默认为0


::

    bounds = [113.75194,22.447837,114.624187,22.864748]
    tbd.grid_params(bounds,accuracy = 500)


.. function:: transbigdata.grid_params_best(data,col = ['lon','lat'],accuracy = 500,gap = 10,sample = 10000)

获取最佳的栅格化参数，以基尼系数最大为标准

**输入**

data : DataFrame
    数据
col : List
    经纬度列
accuracy : number
    网格大小
gap : number
    精度,越大越精确，效果越好，计算量越大
sample : number
    抽样多少数据做测试

**输出**

params : List
    最佳的栅格参数(lonStart,latStart,deltaLon,deltaLat)，或(lonStart,latStart,deltaLon,deltaLat,theta)，其中，lonStart,latStart分别为栅格左下角坐标，deltaLon,deltaLat为单个栅格的经纬度长宽，theta为栅格的角度，不给则默认为0


.. function:: transbigdata.GPS_to_grids(lon,lat,params)

GPS数据对应栅格编号。输入数据的经纬度列与栅格参数，输出对应的栅格编号

**输入**

lon : Series
    经度列
lat : Series
    纬度列
params : List
    栅格参数(lonStart,latStart,deltaLon,deltaLat)，或(lonStart,latStart,deltaLon,deltaLat,theta)，其中，lonStart,latStart分别为栅格左下角坐标，deltaLon,deltaLat为单个栅格的经纬度长宽，theta为栅格的角度，不给则默认为0
                                           
**输出**

LONCOL : Series
    经度栅格编号列
LATCOL : Series
    纬度栅格编号列

::

    data['LONCOL'],data['LATCOL'] = tbd.GPS_to_grids(data['Lng'],data['Lat'],params)

.. function:: transbigdata.grids_centre(loncol,latcol,params)

栅格编号对应栅格中心点经纬度。输入数据的栅格编号与栅格参数，输出对应的栅格中心点

**输入**

LONCOL : Series
    经度栅格编号列
LATCOL : Series
    纬度栅格编号列
params : List
    栅格参数(lonStart,latStart,deltaLon,deltaLat)，或(lonStart,latStart,deltaLon,deltaLat,theta)，其中，lonStart,latStart分别为栅格左下角坐标，deltaLon,deltaLat为单个栅格的经纬度长宽，theta为栅格的角度，不给则默认为0
                                           
**输出**

HBLON : Series
    栅格中心点经度列
HBLAT : Series
    栅格中心点纬度列


::

    data['HBLON'],data['HBLAT'] = tbd.grids_centre(data['LONCOL'],data['LATCOL'],params)

.. function:: transbigdata.gridid_to_polygon(loncol,latcol,params)

栅格编号生成栅格的地理信息列。输入数据的栅格编号与栅格参数，输出对应的地理信息列

**输入**

LONCOL : Series
    经度栅格编号列
LATCOL : Series
    纬度栅格编号列
params : List
    栅格参数(lonStart,latStart,deltaLon,deltaLat)，或(lonStart,latStart,deltaLon,deltaLat,theta)，其中，lonStart,latStart分别为栅格左下角坐标，deltaLon,deltaLat为单个栅格的经纬度长宽，theta为栅格的角度，不给则默认为0
                                           
**输出**

geometry : Series
    栅格的矢量图形列

::

    data['geometry'] = tbd.gridid_to_polygon(data['LONCOL'],data['LATCOL'],params)

.. function:: transbigdata.gridid_sjoin_shape(data,shape,params,col = ['LONCOL','LATCOL'])

输入数据（带有栅格经纬度编号两列），矢量图形与栅格化参数，输出数据栅格并对应矢量图形。

**输入**

data : DataFrame
    数据,（带有栅格经纬度编号两列）
shape : GeoDataFrame
    矢量图形
params : List
    栅格化参数
col : List
    列名，[经度栅格编号，纬度栅格编号]

**输出**

data1 : DataFrame
    数据栅格并对应矢量图形


.. function:: transbigdata.regenerate_params(grid):

从栅格数据重新生成栅格参数  

**输入**  
grid : GeoDataFrame  
    transbigdata中生成的grid                 

**输出**  
params : List  
    栅格参数(lonStart,latStart,deltaLon,deltaLat)，或(lonStart,latStart,deltaLon,deltaLat,theta)，其中，lonStart,latStart分别为栅格左下角坐标，deltaLon,deltaLat为单个栅格的经纬度长宽，theta为栅格的角度，不给则默认为0  


geohash编码
==============

geohash是一种公共域地理编码系统，它的作用是将经纬度地理位置编码为字母和数字组成的字符串，字符串也可解码为经纬度。每个字符串代表一个网格编号，字符串的长度越长则精度越高。根据\ `wiki <https://en.wikipedia.org/wiki/Geohash>`__\ ，geohash字符串长度对应精度表格如下：

========================= ======== ======== ========= ========= ========
geohash length(precision) lat bits lng bits lat error lng error km error
========================= ======== ======== ========= ========= ========
1                         2        3        ±23       ±23       ±2500
2                         5        5        ±2.8      ±5.6      ±630
3                         7        8        ±0.70     ±0.70     ±78
4                         10       10       ±0.087    ±0.18     ±20
5                         12       13       ±0.022    ±0.022    ±2.4
6                         15       15       ±0.0027   ±0.0055   ±0.61
7                         17       18       ±0.00068  ±0.00068  ±0.076
8                         20       20       ±0.000085 ±0.00017  ±0.019
========================= ======== ======== ========= ========= ========

TransBigData包中也提供了geohash的处理功能，主要包括三个函数：


.. function:: transbigdata.geohash_encode(lon,lat,precision=12)

输入经纬度与精度，输出geohash编码

**输入**

lon : Series
    经度列
lat : Series
    纬度列
precision : number
    geohash精度                       

**输出**

lon : Series
    经度列
lat : Series
    纬度列


.. function:: transbigdata.geohash_decode(geohash)

输入经纬度与精度，输出geohash编码

**输入**

geohash : Series
    geohash编码列                    

**输出**

geohash : Series
    geohash编码列

.. function:: transbigdata.geohash_togrid(geohash)

输入geohash编码，输出geohash网格的地理信息图形Series列

**输入**

geohash : Series
    geohash编码列                    

**输出**

poly : Series
    geohash的栅格列

相比TransBigData包中提供的方形栅格处理方法，geohash更慢，也无法提供自由定义的栅格大小。下面的示例展示如何利用这三个函数对数据进行geohash编码集计，并可视化

::

    import transbigdata as tbd
    import pandas as pd
    import geopandas as gpd
    #读取数据    
    data = pd.read_csv('TaxiData-Sample.csv',header = None) 
    data.columns = ['VehicleNum','time','slon','slat','OpenStatus','Speed'] 

::

    #依据经纬度geohash编码，精确度选6时，栅格大小约为±0.61km
    data['geohash'] = tbd.geohash_encode(data['slon'],data['slat'],precision=6)
    data['geohash']




.. parsed-literal::

    0         ws0btw
    1         ws0btz
    2         ws0btz
    3         ws0btz
    4         ws0by4
               ...  
    544994    ws131q
    544995    ws1313
    544996    ws131f
    544997    ws1361
    544998    ws10tq
    Name: geohash, Length: 544999, dtype: object



::

    #基于geohash编码集计
    dataagg = data.groupby(['geohash'])['VehicleNum'].count().reset_index()
    #geohash编码解码为经纬度
    dataagg['lon_geohash'],dataagg['lat_geohash'] = tbd.geohash_decode(dataagg['geohash'])
    #geohash编码生成栅格矢量图形
    dataagg['geometry'] = tbd.geohash_togrid(dataagg['geohash'])
    #转换为GeoDataFrame
    dataagg = gpd.GeoDataFrame(dataagg)
    dataagg




.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>geohash</th>
          <th>VehicleNum</th>
          <th>lon_geohash</th>
          <th>lat_geohash</th>
          <th>geometry</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>0</th>
          <td>w3uf3x</td>
          <td>1</td>
          <td>108.</td>
          <td>10.28</td>
          <td>POLYGON ((107.99561 10.27771, 107.99561 10.283...</td>
        </tr>
        <tr>
          <th>1</th>
          <td>webzz6</td>
          <td>12</td>
          <td>113.9</td>
          <td>22.47</td>
          <td>POLYGON ((113.87329 22.46704, 113.87329 22.472...</td>
        </tr>
        <tr>
          <th>2</th>
          <td>webzz7</td>
          <td>21</td>
          <td>113.9</td>
          <td>22.48</td>
          <td>POLYGON ((113.87329 22.47253, 113.87329 22.478...</td>
        </tr>
        <tr>
          <th>3</th>
          <td>webzzd</td>
          <td>1</td>
          <td>113.9</td>
          <td>22.47</td>
          <td>POLYGON ((113.88428 22.46704, 113.88428 22.472...</td>
        </tr>
        <tr>
          <th>4</th>
          <td>webzzf</td>
          <td>2</td>
          <td>113.9</td>
          <td>22.47</td>
          <td>POLYGON ((113.89526 22.46704, 113.89526 22.472...</td>
        </tr>
        <tr>
          <th>...</th>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
        </tr>
        <tr>
          <th>2022</th>
          <td>ws1d9u</td>
          <td>1</td>
          <td>114.7</td>
          <td>22.96</td>
          <td>POLYGON ((114.68628 22.96143, 114.68628 22.966...</td>
        </tr>
        <tr>
          <th>2023</th>
          <td>ws1ddh</td>
          <td>6</td>
          <td>114.7</td>
          <td>22.96</td>
          <td>POLYGON ((114.69727 22.96143, 114.69727 22.966...</td>
        </tr>
        <tr>
          <th>2024</th>
          <td>ws1ddj</td>
          <td>2</td>
          <td>114.7</td>
          <td>22.97</td>
          <td>POLYGON ((114.69727 22.96692, 114.69727 22.972...</td>
        </tr>
        <tr>
          <th>2025</th>
          <td>ws1ddm</td>
          <td>4</td>
          <td>114.7</td>
          <td>22.97</td>
          <td>POLYGON ((114.70825 22.96692, 114.70825 22.972...</td>
        </tr>
        <tr>
          <th>2026</th>
          <td>ws1ddq</td>
          <td>7</td>
          <td>114.7</td>
          <td>22.98</td>
          <td>POLYGON ((114.70825 22.97241, 114.70825 22.977...</td>
        </tr>
      </tbody>
    </table>
    <p>2027 rows × 5 columns</p>
    </div>



::

    #设定绘图边界
    bounds = [113.6,22.4,114.8,22.9]
    #创建图框
    import matplotlib.pyplot as plt
    import plot_map
    fig =plt.figure(1,(8,8),dpi=280)
    ax =plt.subplot(111)
    plt.sca(ax)
    #添加地图底图
    tbd.plot_map(plt,bounds,zoom = 12,style = 4)
    #绘制colorbar
    cax = plt.axes([0.05, 0.33, 0.02, 0.3])
    plt.title('count')
    plt.sca(ax)
    #绘制geohash的栅格集计
    dataagg.plot(ax = ax,column = 'VehicleNum',cax = cax,legend = True)
    #添加比例尺和指北针
    tbd.plotscale(ax,bounds = bounds,textsize = 10,compasssize = 1,accuracy = 2000,rect = [0.06,0.03],zorder = 10)
    plt.axis('off')
    plt.xlim(bounds[0],bounds[2])
    plt.ylim(bounds[1],bounds[3])
    plt.show()



.. image:: geohash/output_9_0.png




六边形渔网生成
=============================

.. function:: transbigdata.hexagon_grids(bounds,accuracy = 500)

生成研究范围内的六边形渔网。

**输入**

bounds : List
    生成范围的边界，[lon1,lat1,lon2,lat2] (WGS84坐标系) 其中，lon1,lat1是左下角坐标，lon2,lat2是右上角坐标 
accuracy : number
    六边形的边长（米）
                                           
**输出**

hexagon : GeoDataFrame
    六边形渔网的矢量图形

::

    
    #设定范围
    bounds = [113.6,22.4,114.8,22.9]
    hexagon = tbd.hexagon_grids(bounds,accuracy = 5000)
    hexagon.plot()

.. image:: _static/WX20211021-201747@2x.png
   :height: 200

******************************
数据可视化
******************************

在jupyter中显示可视化的设置
--------------------------------------

| TransBigData包也依托于kepler.gl提供的可视化插件提供了一键数据整理与可视化的方法
| 使用此功能请先安装python的keplergl包

::

    pip install keplergl

如果要在jupyter notebook中显示可视化，则需要勾选jupyter-js-widgets（可能需要另外安装）和keplergl-jupyter两个插件

.. image:: _static/jupytersettings.png

数据点分布可视化
-------------------

.. function:: transbigdata.visualization_data(data,col =  ['lon','lat'],accuracy = 500,height = 500,maptype = 'point',zoom = 'auto')

输入数据点，集计并可视化

**输入**

data : DataFrame
    数据点分布
col : List
    列名，可输入不带权重的OD，按[经度，纬度]的顺序，此时会自动集计。
    也可输入带权重的OD，按[经度，纬度，数量]的顺序。
zoom : number
    地图缩放等级,默认'auto'自动选择
height : number
    地图图框高度
accuracy : number
    集计的栅格大小
maptype : str
    出图类型，'point'或者'heatmap'

**输出**

vmap : keplergl.keplergl.KeplerGl
    keplergl提供的可视化

使用方法::

    import transbigdata as tbd
    import pandas as pd
    #读取数据    
    data = pd.read_csv('TaxiData-Sample.csv',header = None) 
    data.columns = ['VehicleNum','Time','Lng','Lat','OpenStatus','Speed']
    #可视化数据点分布
    tbd.visualization_data(data,col = ['Lng','Lat'],accuracy=300)

.. image:: example-taxi/datavis.png


轨迹可视化
-------------------

.. function:: transbigdata.visualization_trip(trajdata,col = ['Lng','Lat','ID','Time'],zoom = 10,height=500)

输入轨迹数据与列名，生成kepler的可视化

**输入**

trajdata : DataFrame
    轨迹点数据
col : List
    列名，按[经度,纬度,轨迹编号,时间]的顺序
zoom : number
    地图缩放等级
height : number
    地图图框高度

**输出**

vmap : keplergl.keplergl.KeplerGl
    keplergl提供的可视化

使用方法

::

    import transbigdata as tbd
    import pandas as pd
    #读取数据    
    data = pd.read_csv('TaxiData-Sample.csv',header = None) 
    data.columns = ['VehicleNum','Time','Lng','Lat','OpenStatus','Speed']  
    #轨迹数据可视化
    tbd.visualization_trip(data,col = ['Lng', 'Lat', 'VehicleNum', 'Time'])

.. image:: example-taxi/kepler-traj.png

OD可视化
--------------------

.. function:: transbigdata.visualization_od(oddata,col = ['slon','slat','elon','elat'],zoom = 'auto',height=500,accuracy = 500,mincount = 0)

输入od数据与列名，生成kepler的可视化

**输入**

oddata : DataFrame
    od数据
col : List
    列名，可输入不带权重的OD，按[起点经度，起点纬度，终点经度，终点纬度]的顺序，此时会自动集计。
    也可输入带权重的OD，按[起点经度，起点纬度，终点经度，终点纬度，数量]的顺序。
zoom : number
    地图缩放等级,默认'auto'自动选择
height : number
    地图图框高度
accuracy : number
    集计的栅格大小
mincount : number
    最小的od数，少于这个的od就不显示了

**输出**

vmap : keplergl.keplergl.KeplerGl
    keplergl提供的可视化

使用方法

::

    import transbigdata as tbd
    import pandas as pd
    #读取数据    
    data = pd.read_csv('TaxiData-Sample.csv',header = None) 
    data.columns = ['VehicleNum','Time','Lng','Lat','OpenStatus','Speed']
    #提取OD
    oddata = tbd.taxigps_to_od(data,col = ['VehicleNum','Time','Lng','Lat','OpenStatus'])
    #OD可视化
    tbd.visualization_od(oddata)

.. image:: example-taxi/odvisualization.png
******************************
数据获取
******************************

获取公交线路
=============================

.. function:: transbigdata.getbusdata(city,keywords,accurate=True)

通过输入城市与关键词，获取公交线路的线型与站点

**输入**

city : str
    城市
keywords : str或List
    关键词，线路名称，支持单个关键词或多个关键词
accurate : bool
    是否精确匹配，默认True

**输出**

data : GeoDataFrame
    生成的公交线路
stop : GeoDataFrame
    生成的公交站点

获取行政区划
=============================

.. function:: transbigdata.getadmin(keyword,ak,subdistricts = False)

输入关键词与高德ak，抓取行政区划gis

**输入**

keywords : str
    关键词，可以是名称，如"深圳市"，或行政区划编号，如440500
ak : str
    高德ak
subdistricts : bool
    是否输出子行政区划的信息

**输出**

admin : GeoDataFrame
    行政区划信息
districts : DataFrame
    子行政区划的信息，利用这个可以进一步抓下一级的行政区划

获取等时圈
=============================

.. function:: transbigdata.get_isochrone_amap(lon,lat,reachtime,ak,mode=2)

获取高德地图等时圈，支持`公交`、`地铁`、`公交+地铁`三种模式

**输入**

lon : float
    起点经度(WGS84)
lat : float
    起点纬度(WGS84)
reachtime : number
    等时圈时间
ak : str
    高德地图ak
mode : int or str
    出行方式，0`公交`、1`地铁`、2`公交+地铁`

**输出**

isochrone : GeoDataFrame
    等时圈的GeoDataFrame(WGS84)

.. function:: transbigdata.get_isochrone_mapbox(lon,lat,reachtime,access_token='auto',mode = 'driving')

获取mapbox地图等时圈，支持驾车、步行、骑行

**输入**

lon : float
    起点经度(WGS84)
lat : float
    起点纬度(WGS84)
reachtime : number
    等时圈时间
access_token : str
    Mapbox的access token，如果设置为 `auto`则会自动读取已经保存的access token
mode : bool
    出行方式，取值为 `driving`， `walking` 或 `cycling`

**输出**

isochrone : GeoDataFrame
    等时圈的GeoDataFrame(WGS84).. transbigdata documentation master file, created by
   sphinx-quickstart on Thu Oct 21 14:41:25 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

TransBigData 为交通时空大数据而生
========================================

.. image:: _static/logo-wordmark-dark.png



**主要功能**

TransBigData工具针对时空大数据处理而开发，依托于GeoPandas。TransBigData集成了交通时空大数据处理过程中常用的方法。包括栅格化、数据质量分析、数据预处理、数据集计、轨迹分析、GIS处理、地图底图加载、坐标与距离计算、数据可视化等通用方法。TransBigData也针对出租车GPS数据、共享单车数据、公交GPS数据等多种常见交通时空大数据提供了快速简洁的处理方法。

**技术特点**

* 面向交通时空大数据分析不同阶段的处理需求提供不同处理功能。
* 代码简洁、高效、灵活、易用，通过简短的代码即可实现复杂的数据任务。


TransBigData简介
==============================

快速入门
---------------

| 在安装TransBigData之前，请确保已经安装了可用的geopandas包：https://geopandas.org/index.html
| 如果你已经安装了geopandas，则直接在命令提示符中运行下面代码即可安装

::

    pip install -U transbigdata

下面例子展示如何使用TransBigData工具快速地从出租车GPS数据中提取出行OD
::

    #导入TransBigData包
    import transbigdata as tbd
    #读取数据    
    import pandas as pd
    data = pd.read_csv('TaxiData-Sample.csv',header = None) 
    data.columns = ['VehicleNum','time','slon','slat','OpenStatus','Speed'] 
    data

.. image:: _static/WX20211021-192131@2x.png
   :height: 300


使用*tbd.taxigps_to_od*方法，传入对应的列名，即可提取出行OD
::

    #从GPS数据提取OD
    oddata = tbd.taxigps_to_od(data,col = ['VehicleNum','time','slon','slat','OpenStatus'])
    oddata

.. image:: _static/WX20211021-190104@2x.png
   :height: 300

对提取出的OD进行OD的栅格集计::

    #定义研究范围
   bounds = [113.6,22.4,114.8,22.9]
   #输入研究范围边界bounds与栅格宽度accuracy，获取栅格化参数
   params = tbd.grid_params(bounds = bounds,accuracy = 1500)
   #栅格化OD并集计
   od_gdf = tbd.odagg_grid(oddata,params)
   od_gdf.plot(column = 'count')

.. image:: _static/WX20211021-190524@2x.png
   :height: 300

使用示例
---------------
.. raw:: html
   :file: gallery/html/grid.html



相关链接
---------------

* 小旭学长的b站： https://space.bilibili.com/3051484
* 小旭学长的七天入门交通时空大数据分析课程（零基础免费课）： https://www.lifangshuju.com/#/introduce/166  
* 小旭学长的交通时空大数据分析课程： https://www.lifangshuju.com/#/introduce/154  
* 小旭学长的数据可视化课程： https://www.lifangshuju.com/#/introduce/165  
* 本项目的github页面： https://github.com/ni1o1/transbigdata/  
* 有bug请在这个页面提交： https://github.com/ni1o1/transbigdata/issues

安装
=========================

.. toctree::
   :caption: 安装
   :maxdepth: 2
   
   getting_started.rst


使用示例
=========================

.. toctree::
   :caption: 使用示例
   :maxdepth: 2

   example-taxi/example-taxi.rst
   example-busgps/example-busgps.rst
   metromodel/metromodel.rst
   Example-pNEUMA/Example-pNEUMA.rst
   example-bikesharing/example-bikesharing.rst
   
通用方法
=========================

.. toctree::
   :caption: 通用方法
   :maxdepth: 2
   
   quality.rst
   preprocess.rst
   grids.rst
   odprocess.rst
   visualization.rst
   getbusdata.rst
   traj.rst
   gisprocess.rst
   plot_map.rst
   CoordinatesConverter.rst


各类数据处理方法
=========================

.. toctree::
   :caption: 各类数据处理方法
   :maxdepth: 2

   taxigps.rst
   bikedata.rst
   busgps.rst
   metroline.rst



.. _quality:


******************************
数据质量分析
******************************



.. function:: transbigdata.data_summary(data,col = ['Vehicleid','Time'],show_sample_duration = False)

输入数据，打印数据概况

**输入**

data : DataFrame
    轨迹点数据
col : List
    列名，按[个体ID，时间]的顺序
show_sample_duration : bool
    是否输出个体采样间隔信息
roundnum : number
    小数点取位数
    
使用方法

::

    import transbigdata as tbd
    import pandas as pd

    # Read data
    data = pd.read_csv('TaxiData-Sample.csv',header = None)
    data.columns = ['Vehicleid','Time','Lng','Lat','OpenStatus','Speed']
    data['Time'] = pd.to_datetime(data['Time'])

    # The sampling interval
    tbd.data_summary(data,col = ['Vehicleid','Time'],show_sample_duration=True)

::

    Amount of data
    -----------------
    Total number of data items:  544999
    Total number of individuals:  180
    Data volume of individuals(Mean):  3027.7722
    Data volume of individuals(Upper quartile):  4056.25
    Data volume of individuals(Median):  2600.5
    Data volume of individuals(Lower quartile):  1595.75

    Data time period
    -----------------
    Start time:  2022-01-09 00:00:00
    End time:  2022-01-09 23:59:59

    Sampling interval
    -----------------
    Mean:  27.995 s
    Upper quartile:  30.0 s
    Median:  20.0 s
    Lower quartile:  15.0 s

.. function:: transbigdata.sample_duration(data,col = ['Vehicleid','Time']):

统计数据采样间隔

**输入**

data : DataFrame
    数据
col : List
    列名，按[个体ID,时间]的顺序

**输出**

sample_duration : DataFrame
    一列的数据表，列名为duration，内容是数据的采样间隔，单位秒.. _taxigps:


******************************
出租车GPS数据处理
******************************

.. function:: transbigdata.taxigps_to_od(data,col = ['VehicleNum','Stime','Lng','Lat','OpenStatus'])

出租车OD提取算法，输入出租车GPS数据,提取OD

**输入**

data : DataFrame
	出租车GPS数据（清洗好的）
col : List            
	数据中各列列名，需要按顺序[车辆id，时间，经度，纬度，载客状态]

.. function:: transbigdata.taxigps_traj_point(data,oddata,col=['Vehicleid', 'Time', 'Lng', 'Lat', 'OpenStatus'])

输入出租车数据与OD数据，提取载客与空载的行驶路径点

**输入**

data : DataFrame
    出租车GPS数据，字段名由col变量指定
oddata : DataFrame
    出租车OD数据
col : List
    列名，按[车辆ID,时间,经度,纬度,载客状态]的顺序

**输出**

data_deliver : DataFrame
    载客轨迹点
data_idle : DataFrame
    空载轨迹点
.. _getting_started:


******************************
安装、依赖与更新日志
******************************

安装
=============================

TransBigData依赖geopandas，在安装TransBigData前需要根据 `这个链接 <https://geopandas.org/en/stable/getting_started.html#installation>`_ 中的方法安装geopandas。
如果你已经安装了geopandas，则直接在命令提示符中运行下面代码即可安装::

  pip install -U transbigdata

在Python中运行下面代码::

  import transbigdata as tbd

依赖包
=============================
TransBigData依赖如下包

* `pandas`
* `shapely`
* `rtree`
* `geopandas`
* `scipy`
* `matplotlib`
* `plot_map`>=0.3.5
* `CoordinatesConverter`>=0.1.4
* `networkx` (optional)
* `igraph` (optional)
* `osmnx` (optional)
* `seaborn` (optional)
* `keplergl` (optional)
* `leuvenmapmatching` (optional)

版本更新
=============================

0.2.3 (2021-11-30)
------------------------
版本更新太快来不及写更新日志：

. 将visualization_dataagg合并至visualization_data中
. 增加置信椭圆绘制功能
. 增加英文文档页面
. 为plot_map地图底图增加10个style，现在需要输入mapboxtoken才能绘制

0.1.28 (2021-11-22)
------------------------
为odagg_grid、odagg_shape增加了权重count的支持
增加可视化方法visualization_dataagg，支持带权重的集计可视化
visualization_dataagg的热力图模式加入更新计划

0.1.25 (2021-11-12)
------------------------
添加轨迹增密函数traj_densify

0.1.24 (2021-11-11)
------------------------
添加id_reindex_disgap方法对数据的ID列重新编号，如果相邻两条记录超过距离，则编号为新id
添加clean_traj方法轨迹数据清洗组合拳，能够有效解决轨迹漂移的问题

0.1.23 (2021-11-10)
------------------------
添加visualization_data方法,添加对keplergl的支持

0.1.22 (2021-11-10)
------------------------
添加visualization_od方法,添加对keplergl的支持

0.1.20 (2021-11-09)
------------------------
添加visualization_trip方法,添加对keplergl的支持

0.1.19 (2021-11-09)
------------------------
添加transform_shape方法,可以将GeoDataFrame整体做坐标转换

0.1.18 (2021-11-09)
------------------------
添加splitline_with_length方法,可以将长线段全部打断

0.1.17 (2021-11-08)
------------------------
添加geohash编码方法

0.1.16 (2021-11-05)
------------------------
添加getadmin方法，输入关键词抓取行政区划信息

0.1.15 (2021-11-04)
------------------------
添加metro_network方法，输入站点信息，输出网络信息

0.1.13 (2021-11-03)
------------------------
增加公交地铁线网拓扑建模的模块，添加split_subwayline方法，用公交/地铁站点对公交/地铁线进行切分，得到断面

0.1.12 (2021-11-02)
------------------------
增加爬虫模块，用getbusdata通过输入城市与关键词，获取公交线路的线型与站点，且不依赖于ak

0.1.11 (2021-11-01)
------------------------
公交GPS模块增加busgps_onewaytime函数，
输入到离站信息表arrive_info与站点信息表stop，计算单程耗时

0.1.10 (2021-10-31)
------------------------
增加公交GPS数据到离站信息识别的方法

0.1.9 (2021-10-28)
------------------------
改进id_reindex方法，sample支持所有模式，同时也添加了timegap和timecol参数

0.1.7 (2021-10-27)
------------------------
增加数据质量分析部分：

* sample_duration采样间隔

0.1.6 (2021-10-25)
------------------------
修正taxigps_traj_point的Bug

0.1.5 (2021-10-25)
------------------------
增加轨迹处理部分：

* taxigps_traj_point  输入出租车数据与OD数据，提取载客与空载的行驶路径点
* points_to_traj 输入轨迹点，生成轨迹线型的GeoDataFrame


0.1.4 (2021-10-24)
------------------------
增加栅格化的gridid_sjoin_shape方法，输入数据（带有栅格经纬度编号两列），矢量图形与栅格化参数，输出数据栅格并对应矢量图形。


0.1.3 (2021-10-23)
------------------------
增加预处理的clean_same,clean_drift,clean_taxi_status方法
为预处理的id_reindex方法加入sample参数

0.1.2 (2021-10-23)
------------------------
更新数据预处理的clean_outofshape方法
增加共享单车数据处理功能，bikedata_to_od提取骑行订单数据与停车数据

0.1.1 (2021-10-22)
------------------------
加入数据预处理的clean_outofbounds，dataagg，id_reindex方法

0.1.0 (2021-10-21)
------------------------
最初版本发布.. _odprocess:


***************
数据聚合集计
***************

数据集计
==========

.. function:: transbigdata.dataagg(data,shape,col = ['Lng','Lat','count'],accuracy=500)

数据集计至小区

**输入**

data : DataFrame
    数据
shape : GeoDataFrame
	小区
col : List
    可传入经纬度两列，如['Lng','Lat']，此时每一列权重为1。也可以传入经纬度和计数列三列，如['Lng','Lat','count']
accuracy : number
    计算原理是先栅格化后集计，这里定义栅格大小，越小精度越高

**输出**

aggresult : GeoDataFrame
    小区，其中count列为统计结果
data1 : DataFrame
    数据，对应上了小区

OD集计
==========

.. function:: transbigdata.odagg_grid(oddata,params,col = ['slon','slat','elon','elat'],arrow = False)


OD集计与地理信息生成（栅格）。输入OD数据（每一行数据是一个出行），栅格化OD并集计后生成OD的GeoDataFrame

**输入**

oddata : DataFrame 
    OD数据
col : List
    起终点列名,['slon','slat','elon','elat']，此时每一列权重为1。
    也可以传入权重列，如['slon','slat','elon','elat','count']
params : List
    栅格参数(lonStart,latStart,deltaLon,deltaLat)，分别为栅格左下角坐标与单个栅格的经纬度长宽
arrow : bool
    生成的OD地理线型是否包含箭头

**输出**

oddata1 : GeoDataFrame 
    集计后生成OD的GeoDataFrame

.. function:: transbigdata.odagg_shape(oddata,shape,col = ['slon','slat','elon','elat'],params = None,round_accuracy = 6,arrow = False)

OD集计与地理信息生成（小区集计）。输入OD数据（每一行数据是一个出行），栅格化OD并集计后生成OD的GeoDataFrame

**输入**

oddata : DataFrame 
    OD数据
shape : GeoDataFrame
    集计小区的GeoDataFrame
col : List   
    起终点列名,['slon','slat','elon','elat']，此时每一列权重为1。
    也可以传入权重列，如['slon','slat','elon','elat','count']
params : List 
    栅格化参数，如果传入，则先栅格化后以栅格中心点匹配小区，如果不传入，则直接以经纬度匹配。在数据量大时，用栅格化进行匹配速度会极大提升
round_accuracy : number
    集计时经纬度取小数位数
arrow : bool       
    生成的OD地理线型是否包含箭头

**输出**

oddata1 : GeoDataFrame 
    集计后生成OD的GeoDataFrame


.. _busgps:


******************************
公交车GPS数据处理
******************************

.. function:: transbigdata.busgps_arriveinfo(data,line,stop,col = ['VehicleId','GPSDateTime','lon','lat','stopname'],stopbuffer = 200,mintime = 300,project_epsg = 2416,timegap = 1800,method = 'project',projectoutput = False)

输入公交GPS数据、公交线路与站点的GeoDataFrame，该方法能够识别公交的到离站信息

**输入**

data : DataFrame
    公交GPS数据，单一公交线路，且需要含有车辆ID、GPS时间、经纬度（wgs84）
line : GeoDataFrame
    公交线型的GeoDataFrame数据，单一公交线路
stop : GeoDataFrame
    公交站点的GeoDataFrame数据
col : List
    列名，按[车辆ID,时间,经度,纬度，站点名称字段]的顺序
stopbuffer : number
    米，站点的一定距离范围，车辆进入这一范围视为到站，离开则视为离站
mintime : number
    秒，短时间内公交再次到站则需要与前一次的到站数据结合一起计算到离站时间，该参数设置阈值
project_epsg : number
    匹配时会将数据转换为投影坐标系以计算距离，这里需要给定投影坐标系的epsg代号
timegap : number
    秒，清洗数据用，多长时间车辆不出现，就视为新的车辆
method : str
    公交运行图匹配方法，可选'project'或'dislimit'；
    project为直接匹配线路上最近点，匹配速度快；
    dislimit则需要考虑前面点位置，加上距离限制，匹配速度慢。
projectoutput : bool
    是否输出投影后的数据

**输出**

arrive_info : DataFrame
    公交到离站信息

.. function:: transbigdata.busgps_onewaytime(arrive_info,start,end,col = ['VehicleId','stopname','arrivetime','leavetime'])

输入到离站信息表arrive_info与站点信息表stop，计算单程耗时

**输入**

arrive_info : DataFrame
    公交到离站数据
stop : GeoDataFrame
    公交站点的GeoDataFrame数据
start : Str
    起点站名字
end : Str
    终点站名字
col : List
    字段列名[车辆ID,站点名称,到站时间,离站时间]


**输出**

onewaytime : DataFrame
    公交单程耗时共享单车数据社区发现
========================================

| 这个案例的Jupyter notebook: `点击这里 <https://github.com/ni1o1/transbigdata/blob/main/example/Example%205-community%20detection%20for%20bikesharing%20data.ipynb>`__.
| 对于共享单车的出行，每一次出行都可以被看作是一个从起点行动到终点的出行过程。当我们把起点和终点视为节点，把它们之间的出行视为边时，就可以构建一个网络。通过分析这个网络，我们可以得到关于城市的空间结构、共享单车需求的宏观出行特征等信息。
| 社区发现，也可以叫图分割，帮助我们揭示网络中节点之间的隐藏关系。在这个例子中，我们将介绍如何将TransBigData整合到共享单车数据的社区发现分析过程中。


::

    import pandas as pd
    import numpy as np
    import geopandas as gpd
    import transbigdata as tbd

数据预处理
-------------------------

在社区发现之前，我们首先需要对数据进行预处理。从共享单车订单中提取出行OD并剔除异常出行，并以清洗好的数据作为研究对象。

::

    #读取共享单车数据
    bikedata = pd.read_csv(r'data/bikedata-sample.csv')
    bikedata.head(5)




.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>BIKE_ID</th>
          <th>DATA_TIME</th>
          <th>LOCK_STATUS</th>
          <th>LONGITUDE</th>
          <th>LATITUDE</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>0</th>
          <td>5</td>
          <td>2018-09-01 0:00:36</td>
          <td>1</td>
          <td>121.363566</td>
          <td>31.259615</td>
        </tr>
        <tr>
          <th>1</th>
          <td>6</td>
          <td>2018-09-01 0:00:50</td>
          <td>0</td>
          <td>121.406226</td>
          <td>31.214436</td>
        </tr>
        <tr>
          <th>2</th>
          <td>6</td>
          <td>2018-09-01 0:03:01</td>
          <td>1</td>
          <td>121.409402</td>
          <td>31.215259</td>
        </tr>
        <tr>
          <th>3</th>
          <td>6</td>
          <td>2018-09-01 0:24:53</td>
          <td>0</td>
          <td>121.409228</td>
          <td>31.214427</td>
        </tr>
        <tr>
          <th>4</th>
          <td>6</td>
          <td>2018-09-01 0:26:38</td>
          <td>1</td>
          <td>121.409771</td>
          <td>31.214406</td>
        </tr>
      </tbody>
    </table>
    </div>


读取研究区域的边界，并用tbd.clean_outofshape方法剔除研究区域以外的数据

::

    #读取上海行政区划边界
    shanghai_admin = gpd.read_file(r'data/shanghai.json')
    #剔除研究范围外的数据
    bikedata = tbd.clean_outofshape(bikedata, shanghai_admin, col=['LONGITUDE', 'LATITUDE'], accuracy=500)

用tbd.bikedata_to_od方法从单车数据中识别出行OD信息

::

    #识别单车出行OD
    move_data,stop_data = tbd.bikedata_to_od(bikedata,
                       col = ['BIKE_ID','DATA_TIME','LONGITUDE','LATITUDE','LOCK_STATUS'])
    move_data.head(5)




.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>BIKE_ID</th>
          <th>stime</th>
          <th>slon</th>
          <th>slat</th>
          <th>etime</th>
          <th>elon</th>
          <th>elat</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>96</th>
          <td>6</td>
          <td>2018-09-01 0:00:50</td>
          <td>121.406226</td>
          <td>31.214436</td>
          <td>2018-09-01 0:03:01</td>
          <td>121.409402</td>
          <td>31.215259</td>
        </tr>
        <tr>
          <th>561</th>
          <td>6</td>
          <td>2018-09-01 0:24:53</td>
          <td>121.409228</td>
          <td>31.214427</td>
          <td>2018-09-01 0:26:38</td>
          <td>121.409771</td>
          <td>31.214406</td>
        </tr>
        <tr>
          <th>564</th>
          <td>6</td>
          <td>2018-09-01 0:50:16</td>
          <td>121.409727</td>
          <td>31.214403</td>
          <td>2018-09-01 0:52:14</td>
          <td>121.412610</td>
          <td>31.214905</td>
        </tr>
        <tr>
          <th>784</th>
          <td>6</td>
          <td>2018-09-01 0:53:38</td>
          <td>121.413333</td>
          <td>31.214951</td>
          <td>2018-09-01 0:55:38</td>
          <td>121.412656</td>
          <td>31.217051</td>
        </tr>
        <tr>
          <th>1028</th>
          <td>6</td>
          <td>2018-09-01 11:35:01</td>
          <td>121.419261</td>
          <td>31.213414</td>
          <td>2018-09-01 11:35:13</td>
          <td>121.419518</td>
          <td>31.213657</td>
        </tr>
      </tbody>
    </table>
    </div>

我们需要剔除过长与过短的共享单车出行。用tbd.getdistance获取起终点之间的直线距离，并筛选删除直线距离小于100米与大于10千米的出行

::

    #计算骑行直线距离
    move_data['distance'] = tbd.getdistance(move_data['slon'],move_data['slat'],move_data['elon'],move_data['elat'])
    #清洗骑行数据，删除过长与过短的出行
    move_data = move_data[(move_data['distance']>100)&(move_data['distance']<10000)]

接下来，我们以500米×500米的栅格为最小分析单元，用tbd.grid_params方法获取栅格划分参数，再将参数输入tbd.odagg_grid方法，对OD进行栅格集计

::

    # 获取栅格划分参数
    bounds = (120.85, 30.67, 122.24, 31.87)
    params = tbd.grid_params(bounds,accuracy = 500)
    #集计OD
    od_gdf = tbd.odagg_grid(move_data, params, col=['slon', 'slat', 'elon', 'elat'])
    od_gdf.head(5)




.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>SLONCOL</th>
          <th>SLATCOL</th>
          <th>ELONCOL</th>
          <th>ELATCOL</th>
          <th>count</th>
          <th>SHBLON</th>
          <th>SHBLAT</th>
          <th>EHBLON</th>
          <th>EHBLAT</th>
          <th>geometry</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>0</th>
          <td>26</td>
          <td>95</td>
          <td>26</td>
          <td>96</td>
          <td>1</td>
          <td>120.986782</td>
          <td>31.097177</td>
          <td>120.986782</td>
          <td>31.101674</td>
          <td>LINESTRING (120.98678 31.09718, 120.98678 31.1...</td>
        </tr>
        <tr>
          <th>40803</th>
          <td>117</td>
          <td>129</td>
          <td>116</td>
          <td>127</td>
          <td>1</td>
          <td>121.465519</td>
          <td>31.250062</td>
          <td>121.460258</td>
          <td>31.241069</td>
          <td>LINESTRING (121.46552 31.25006, 121.46026 31.2...</td>
        </tr>
        <tr>
          <th>40807</th>
          <td>117</td>
          <td>129</td>
          <td>117</td>
          <td>128</td>
          <td>1</td>
          <td>121.465519</td>
          <td>31.250062</td>
          <td>121.465519</td>
          <td>31.245565</td>
          <td>LINESTRING (121.46552 31.25006, 121.46552 31.2...</td>
        </tr>
        <tr>
          <th>40810</th>
          <td>117</td>
          <td>129</td>
          <td>117</td>
          <td>131</td>
          <td>1</td>
          <td>121.465519</td>
          <td>31.250062</td>
          <td>121.465519</td>
          <td>31.259055</td>
          <td>LINESTRING (121.46552 31.25006, 121.46552 31.2...</td>
        </tr>
        <tr>
          <th>40811</th>
          <td>117</td>
          <td>129</td>
          <td>118</td>
          <td>126</td>
          <td>1</td>
          <td>121.465519</td>
          <td>31.250062</td>
          <td>121.470780</td>
          <td>31.236572</td>
          <td>LINESTRING (121.46552 31.25006, 121.47078 31.2...</td>
        </tr>
      </tbody>
    </table>
    </div>

对OD集计的结果在地图上可视化，用tbd.plot_map加载地图底图，并用tbd.plotscale添加比例尺与指北针

::

    #创建图框
    import matplotlib.pyplot as plt
    import plot_map
    fig =plt.figure(1,(8,8),dpi=300)
    ax =plt.subplot(111)
    plt.sca(ax)
    #添加地图底图
    tbd.plot_map(plt,bounds,zoom = 11,style = 8)
    #绘制colorbar
    cax = plt.axes([0.05, 0.33, 0.02, 0.3])
    plt.title('Data count')
    plt.sca(ax)
    #绘制OD
    od_gdf.plot(ax = ax,column = 'count',cmap = 'Blues_r',linewidth = 0.5,vmax = 10,cax = cax,legend = True)
    #添加比例尺和指北针
    tbd.plotscale(ax,bounds = bounds,textsize = 10,compasssize = 1,textcolor = 'white',accuracy = 2000,rect = [0.06,0.03],zorder = 10)
    plt.axis('off')
    plt.xlim(bounds[0],bounds[2])
    plt.ylim(bounds[1],bounds[3])
    plt.show()



.. image:: output_7_0.png


提取节点信息
------------------------

使用igraph包构建网络。在Python中，igraph与networkx功能类似，都提供了网络分析的功能，仅在部分算法的支持上有所区别。
构建网络时，我们需要向igraph提供网络的节点与边的信息。以OD数据中出现过的每个栅格作为节点，构建节点的信息时，需要为节点创建从0开始的数字编号，代码如下

::

    #把起终点的经纬度栅格编号变为一个字段
    od_gdf['S'] = od_gdf['SLONCOL'].astype(str) + ',' + od_gdf['SLATCOL'].astype(str)
    od_gdf['E'] = od_gdf['ELONCOL'].astype(str) + ',' + od_gdf['ELATCOL'].astype(str)
    #提取节点集合
    node = set(od_gdf['S'])|set(od_gdf['E'])
    #把节点集合变成DataFrame
    node = pd.DataFrame(node)
    #重新编号节点
    node['id'] = range(len(node))
    node




.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>0</th>
          <th>id</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>0</th>
          <td>118,134</td>
          <td>0</td>
        </tr>
        <tr>
          <th>1</th>
          <td>109,102</td>
          <td>1</td>
        </tr>
        <tr>
          <th>2</th>
          <td>59,71</td>
          <td>2</td>
        </tr>
        <tr>
          <th>3</th>
          <td>93,78</td>
          <td>3</td>
        </tr>
        <tr>
          <th>4</th>
          <td>96,17</td>
          <td>4</td>
        </tr>
        <tr>
          <th>...</th>
          <td>...</td>
          <td>...</td>
        </tr>
        <tr>
          <th>9806</th>
          <td>94,97</td>
          <td>9806</td>
        </tr>
        <tr>
          <th>9807</th>
          <td>106,152</td>
          <td>9807</td>
        </tr>
        <tr>
          <th>9808</th>
          <td>124,134</td>
          <td>9808</td>
        </tr>
        <tr>
          <th>9809</th>
          <td>98,158</td>
          <td>9809</td>
        </tr>
        <tr>
          <th>9810</th>
          <td>152,86</td>
          <td>9810</td>
        </tr>
      </tbody>
    </table>
    <p>9811 rows × 2 columns</p>
    </div>



提取边信息
----------------

将新的编号连接到OD信息表上，以提取新ID之间的出行量构成边

::

    #把新编号连接到OD数据上
    node.columns = ['S','S_id']
    od_gdf = pd.merge(od_gdf,node,on = ['S'])
    node.columns = ['E','E_id']
    od_gdf = pd.merge(od_gdf,node,on = ['E'])
    #提取边信息
    edge = od_gdf[['S_id','E_id','count']]
    edge




.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>S_id</th>
          <th>E_id</th>
          <th>count</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>0</th>
          <td>8261</td>
          <td>7105</td>
          <td>1</td>
        </tr>
        <tr>
          <th>1</th>
          <td>9513</td>
          <td>2509</td>
          <td>1</td>
        </tr>
        <tr>
          <th>2</th>
          <td>118</td>
          <td>2509</td>
          <td>3</td>
        </tr>
        <tr>
          <th>3</th>
          <td>348</td>
          <td>2509</td>
          <td>1</td>
        </tr>
        <tr>
          <th>4</th>
          <td>1684</td>
          <td>2509</td>
          <td>1</td>
        </tr>
        <tr>
          <th>...</th>
          <td>...</td>
          <td>...</td>
          <td>...</td>
        </tr>
        <tr>
          <th>68468</th>
          <td>8024</td>
          <td>4490</td>
          <td>2</td>
        </tr>
        <tr>
          <th>68469</th>
          <td>4216</td>
          <td>3802</td>
          <td>2</td>
        </tr>
        <tr>
          <th>68470</th>
          <td>4786</td>
          <td>6654</td>
          <td>2</td>
        </tr>
        <tr>
          <th>68471</th>
          <td>6484</td>
          <td>602</td>
          <td>3</td>
        </tr>
        <tr>
          <th>68472</th>
          <td>7867</td>
          <td>8270</td>
          <td>3</td>
        </tr>
      </tbody>
    </table>
    <p>68473 rows × 3 columns</p>
    </div>



构建网络
--------

导入igraph包，创建网络，添加节点，并将边数据输入网络。同时，为每一条边添加相应的权重

::

    import igraph
    #创建网络
    g = igraph.Graph()
    #在网络中添加节点。
    g.add_vertices(len(node))
    #在网络中添加边。
    g.add_edges(edge[['S_id','E_id']].values)
    #提取边的权重。
    edge_weights = edge[['count']].values
    #给边添加权重。
    for i in range(len(edge_weights)):
        g.es[i]['weight'] = edge_weights[i]

社区发现
-------------

在构建好的网络上应用社区发现算法。其中，我们使用igraph包自带的g.community_multilevel方法实现Fast unfolding社区发现算法。前面我们介绍过，Fast unfolding算法将社区逐层迭代合并直至模块度最优，而在g.community_multilevel方法中可以设定return_levels返回迭代的中间结果。这里我们设定return_levels为False，只返回最终结果进行分析

::

    #社区发现
    g_clustered = g.community_multilevel(weights = edge_weights, return_levels=False)


社区发现的结果存储在g_clustered变量中，可以用内置方法直接计算模块度

::

    #模块度
    g_clustered.modularity




.. parsed-literal::

    0.8496561130926571

一般来说，模块度在0.5以上已经属于较高值。而这一结果的模块度达到0.84，表明网络的社区结构非常明显，社区划分结果也能够很好地划分网络。接下来，我们将社区划分结果赋值到节点信息表上，为后面的可视化做准备。代码如下

::

    #将结果赋值到节点上
    node['group'] = g_clustered.membership
    #重命名列
    node.columns = ['grid','node_id','group']
    node




.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>grid</th>
          <th>node_id</th>
          <th>group</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>0</th>
          <td>118,134</td>
          <td>0</td>
          <td>0</td>
        </tr>
        <tr>
          <th>1</th>
          <td>109,102</td>
          <td>1</td>
          <td>1</td>
        </tr>
        <tr>
          <th>2</th>
          <td>59,71</td>
          <td>2</td>
          <td>2</td>
        </tr>
        <tr>
          <th>3</th>
          <td>93,78</td>
          <td>3</td>
          <td>3</td>
        </tr>
        <tr>
          <th>4</th>
          <td>96,17</td>
          <td>4</td>
          <td>4</td>
        </tr>
        <tr>
          <th>...</th>
          <td>...</td>
          <td>...</td>
          <td>...</td>
        </tr>
        <tr>
          <th>9806</th>
          <td>94,97</td>
          <td>9806</td>
          <td>8</td>
        </tr>
        <tr>
          <th>9807</th>
          <td>106,152</td>
          <td>9807</td>
          <td>36</td>
        </tr>
        <tr>
          <th>9808</th>
          <td>124,134</td>
          <td>9808</td>
          <td>37</td>
        </tr>
        <tr>
          <th>9809</th>
          <td>98,158</td>
          <td>9809</td>
          <td>9</td>
        </tr>
        <tr>
          <th>9810</th>
          <td>152,86</td>
          <td>9810</td>
          <td>26</td>
        </tr>
      </tbody>
    </table>
    <p>9811 rows × 3 columns</p>
    </div>



社区可视化
-------------

在社区发现的结果中，可能会存在部分社区中只存在少量的节点，无法形成规模较大的社区。这些社区为离群点，在可视化之前应该删去，这里我们保留包含10个栅格以上的社区

::

    #统计每个社区的栅格数量
    group = node['group'].value_counts()
    #提取大于10个栅格的社区
    group = group[group>10]
    #只保留这些社区的栅格
    node = node[node['group'].apply(lambda r:r in group.index)]

将栅格编号复原，再用tbd.gridid_to_polygon方法从栅格编号生成栅格的地理几何图形

::

    #切分获取栅格编号
    node['LONCOL'] = node['grid'].apply(lambda r:r.split(',')[0]).astype(int)
    node['LATCOL'] = node['grid'].apply(lambda r:r.split(',')[1]).astype(int)
    #生成栅格地理图形
    node['geometry'] = tbd.gridid_to_polygon(node['LONCOL'],node['LATCOL'],params)
    #转为GeoDataFrame
    import geopandas as gpd
    node = gpd.GeoDataFrame(node)
    node




.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>grid</th>
          <th>node_id</th>
          <th>group</th>
          <th>LONCOL</th>
          <th>LATCOL</th>
          <th>geometry</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>0</th>
          <td>118,134</td>
          <td>0</td>
          <td>0</td>
          <td>118</td>
          <td>134</td>
          <td>POLYGON ((121.46815 31.27030, 121.47341 31.270...</td>
        </tr>
        <tr>
          <th>1</th>
          <td>109,102</td>
          <td>1</td>
          <td>1</td>
          <td>109</td>
          <td>102</td>
          <td>POLYGON ((121.42080 31.12641, 121.42606 31.126...</td>
        </tr>
        <tr>
          <th>3</th>
          <td>93,78</td>
          <td>3</td>
          <td>3</td>
          <td>93</td>
          <td>78</td>
          <td>POLYGON ((121.33663 31.01849, 121.34189 31.018...</td>
        </tr>
        <tr>
          <th>4</th>
          <td>96,17</td>
          <td>4</td>
          <td>4</td>
          <td>96</td>
          <td>17</td>
          <td>POLYGON ((121.35241 30.74419, 121.35767 30.744...</td>
        </tr>
        <tr>
          <th>5</th>
          <td>156,117</td>
          <td>5</td>
          <td>5</td>
          <td>156</td>
          <td>117</td>
          <td>POLYGON ((121.66806 31.19385, 121.67332 31.193...</td>
        </tr>
        <tr>
          <th>...</th>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
        </tr>
        <tr>
          <th>9806</th>
          <td>94,97</td>
          <td>9806</td>
          <td>8</td>
          <td>94</td>
          <td>97</td>
          <td>POLYGON ((121.34189 31.10392, 121.34715 31.103...</td>
        </tr>
        <tr>
          <th>9807</th>
          <td>106,152</td>
          <td>9807</td>
          <td>36</td>
          <td>106</td>
          <td>152</td>
          <td>POLYGON ((121.40502 31.35124, 121.41028 31.351...</td>
        </tr>
        <tr>
          <th>9808</th>
          <td>124,134</td>
          <td>9808</td>
          <td>37</td>
          <td>124</td>
          <td>134</td>
          <td>POLYGON ((121.49971 31.27030, 121.50498 31.270...</td>
        </tr>
        <tr>
          <th>9809</th>
          <td>98,158</td>
          <td>9809</td>
          <td>9</td>
          <td>98</td>
          <td>158</td>
          <td>POLYGON ((121.36293 31.37822, 121.36819 31.378...</td>
        </tr>
        <tr>
          <th>9810</th>
          <td>152,86</td>
          <td>9810</td>
          <td>26</td>
          <td>152</td>
          <td>86</td>
          <td>POLYGON ((121.64702 31.05446, 121.65228 31.054...</td>
        </tr>
      </tbody>
    </table>
    <p>8527 rows × 6 columns</p>
    </div>

在这一步中，我们将每一个节点复原为栅格，标记上节点所属的社区编号，生成了每个栅格的地理信息，并将其转换为GeoDataFrame，可以用如下代码绘制栅格，测试是否生成成功

::

    node.plot('group')




.. image:: output_22_1.png


这里我们将group字段的分组编号映射到颜色上进行初步可视化，不同分组的颜色不同。从结果的图中可以看到，相同颜色的栅格在地理空间上大多聚集在一起，表明共享单车的空间联系可以将地理空间上接近的区域紧密地联系在一起

前面的结果可视化的效果并不明显，我们并不能从图中清晰地看出分区的情况。接下来，我们可以对分区结果进行一定的调整与可视化。可视化的调整主要有以下思路:

* 比较合适的分区结果应该是每个区域都为空间上连续的区域，在初步的可视化结果中，有不少的栅格在空间上为孤立存在，这些点应该予以剔除。
* 在可视化结果中，我们可以将同一个组别的栅格合并，为每个分区形成面要素，这样在下一步可视化中就可以绘制出分区的边界。
* 在分区结果中，有些区域的内部可能会存在其他区域的“飞地”，即隶属于本分区，却被其他分区所包围，只能“飞”过其他分区的属地，才能到达自己的飞地。这种分区在共享单车的实际运营中也是难以管理的，应该避免这种情况的出现。

解决上述问题，我们可以使用TransBigData所提供的两个GIS处理方法，tbd.merge_polygon和tbd.polyon_exterior。其中tbd.merge_polygon能够将同一个组别的面要素进行合并，而tbd.polyon_exterior则可以对多边形取外边界后再构成新的多边形，以此剔除飞地。同时，也可以设定最小面积，对小于此面积的面要素进行剔除。代码如下


::

    #以group字段为分组，将同一组别的面要素合并
    node_community = tbd.merge_polygon(node,'group')
    #输入多边形GeoDataFrame数据，对多边形取外边界构成新多边形
    #设定最小面积minarea，小于该面积的面全部剔除，避免大量离群点出现
    node_community = tbd.polyon_exterior(node_community,minarea = 0.000100)


处理好社区的面要素后，接下来需要对面要素进行可视化。我们希望对不同的面赋予不同的颜色。在可视化章节中我们提到，在显示的要素没有数值上的大小区别时，颜色的选择上需要保持它们各自的颜色具有相同的亮度与饱和度。而用seaborn的调色盘方法即可快速地生成同一亮度与饱和度下的多种颜色

::

    #生成调色盘
    import seaborn as sns
    ## l: 亮度
    ## s: 饱和度
    cmap = sns.hls_palette(n_colors=len(node_community), l=.7, s=0.8)
    sns.palplot(cmap)



.. image:: output_24_0.png

对社区结果进行可视化

::

    #创建图框
    import matplotlib.pyplot as plt
    import plot_map
    fig =plt.figure(1,(8,8),dpi=300)
    ax =plt.subplot(111)
    plt.sca(ax)
    #添加地图底图
    tbd.plot_map(plt,bounds,zoom = 10,style = 6)
    #设定colormap
    from matplotlib.colors import ListedColormap 
    #打乱社区的排列顺序
    node_community = node_community.sample(frac=1)
    #绘制社区
    node_community.plot(cmap = ListedColormap(cmap),ax = ax,edgecolor = '#333',alpha = 0.8)
    #添加比例尺和指北针
    tbd.plotscale(ax,bounds = bounds,textsize = 10,compasssize = 1,textcolor = 'k'
                  ,accuracy = 2000,rect = [0.06,0.03],zorder = 10)
    plt.axis('off')
    plt.xlim(bounds[0],bounds[2])
    plt.ylim(bounds[1],bounds[3])
    plt.show()



.. image:: output_25_0.png

至此，我们就已经成功地可视化出共享单车社区，并绘制出每一个社区的边界。在用社区发现模型进行分区时，并没有往模型中输入任何地理空间信息，模型对研究区域的分割也仅仅依靠共享单车出行需求所构成的网络联系。
公交GPS的到离站信息匹配
=======================

| 这个案例的Jupyter notebook: `点击这里 <https://github.com/ni1o1/transbigdata/blob/main/example/Example%202-Identifying%20arrival%20and%20departure%20information%20from%20Bus%20GPS%20data.ipynb>`__.
| 可以点击 `这个链接 <https://mybinder.org/v2/gh/ni1o1/transbigdata/9507de936806c34a4befd74aa9227b012569a6aa?urlpath=lab%2Ftree%2Fexample%2FExample%202-Identifying%20arrival%20and%20departure%20information%20from%20Bus%20GPS%20data.ipynb>`__ 在线编辑器中尝试
下面的案例展示如何用TransBigData包处理公交GPS数据，以内置方法计算公交车辆的到离站信息、统计公交单程耗时与运营车速

::

    import transbigdata as tbd
    import pandas as pd
    import geopandas as gpd

读取数据
--------

读取GPS数据

::

    BUS_GPS= pd.read_csv(r'busgps.csv',header = None)
    BUS_GPS.columns = ['GPSDateTime', 'LineId', 'LineName', 'NextLevel', 'PrevLevel',
           'Strlatlon', 'ToDir', 'VehicleId', 'VehicleNo', 'unknow']
    #时间转换为datetime格式
    BUS_GPS['GPSDateTime'] = pd.to_datetime(BUS_GPS['GPSDateTime'])

经纬度坐标转换

::

    #切分经纬度的字符串
    BUS_GPS['lon'] = BUS_GPS['Strlatlon'].apply(lambda r:r.split(',')[0])
    BUS_GPS['lat'] = BUS_GPS['Strlatlon'].apply(lambda r:r.split(',')[1])
    #坐标系转换
    BUS_GPS['lon'],BUS_GPS['lat'] = tbd.gcj02towgs84(BUS_GPS['lon'].astype(float),BUS_GPS['lat'].astype(float))
    BUS_GPS.head(5)




.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>GPSDateTime</th>
          <th>LineId</th>
          <th>LineName</th>
          <th>NextLevel</th>
          <th>PrevLevel</th>
          <th>Strlatlon</th>
          <th>ToDir</th>
          <th>VehicleId</th>
          <th>VehicleNo</th>
          <th>unknow</th>
          <th>lon</th>
          <th>lat</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>0</th>
          <td>2019-01-16 23:59:59</td>
          <td>7100</td>
          <td>71</td>
          <td>2</td>
          <td>1</td>
          <td>121.335413,31.173188</td>
          <td>1</td>
          <td>沪D-R7103</td>
          <td>Z5A-0021</td>
          <td>1</td>
          <td>121.330858</td>
          <td>31.175129</td>
        </tr>
        <tr>
          <th>1</th>
          <td>2019-01-17 00:00:00</td>
          <td>7100</td>
          <td>71</td>
          <td>2</td>
          <td>1</td>
          <td>121.334616,31.172271</td>
          <td>1</td>
          <td>沪D-R1273</td>
          <td>Z5A-0002</td>
          <td>1</td>
          <td>121.330063</td>
          <td>31.174214</td>
        </tr>
        <tr>
          <th>2</th>
          <td>2019-01-17 00:00:00</td>
          <td>7100</td>
          <td>71</td>
          <td>24</td>
          <td>23</td>
          <td>121.339955,31.173025</td>
          <td>0</td>
          <td>沪D-R5257</td>
          <td>Z5A-0020</td>
          <td>1</td>
          <td>121.335390</td>
          <td>31.174958</td>
        </tr>
        <tr>
          <th>3</th>
          <td>2019-01-17 00:00:01</td>
          <td>7100</td>
          <td>71</td>
          <td>14</td>
          <td>13</td>
          <td>121.409491,31.20433</td>
          <td>0</td>
          <td>沪D-R5192</td>
          <td>Z5A-0013</td>
          <td>1</td>
          <td>121.404843</td>
          <td>31.206179</td>
        </tr>
        <tr>
          <th>4</th>
          <td>2019-01-17 00:00:03</td>
          <td>7100</td>
          <td>71</td>
          <td>15</td>
          <td>14</td>
          <td>121.398615,31.200253</td>
          <td>0</td>
          <td>沪D-T0951</td>
          <td>Z5A-0022</td>
          <td>1</td>
          <td>121.393966</td>
          <td>31.202103</td>
        </tr>
      </tbody>
    </table>
    </div>



读取公交线数据

::

    shp = r'busline.json'
    linegdf = gpd.GeoDataFrame.from_file(shp,encoding = 'gbk')
    line = linegdf.iloc[:1].copy()
    line.plot()









.. image:: output_8_1.png


读取公交站点数据

::

    shp = r'busstop.json'
    stop = gpd.GeoDataFrame.from_file(shp,encoding = 'gbk')
    stop = stop[stop['linename'] == '71路(延安东路外滩-申昆路枢纽站)']
    stop.plot()









.. image:: output_10_1.png


到离站信息匹配
--------------

::

    arriveinfo = tbd.busgps_arriveinfo(BUS_GPS,line,stop)



数据清洗中...

运行位置匹配中......

匹配到离站信息.........................................................................................................................................................

::

    arriveinfo




.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>arrivetime</th>
          <th>leavetime</th>
          <th>stopname</th>
          <th>VehicleId</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>0</th>
          <td>2019-01-17 07:19:42</td>
          <td>2019-01-17 07:31:14</td>
          <td>延安东路外滩</td>
          <td>1</td>
        </tr>
        <tr>
          <th>1</th>
          <td>2019-01-17 09:53:08</td>
          <td>2019-01-17 10:09:34</td>
          <td>延安东路外滩</td>
          <td>1</td>
        </tr>
        <tr>
          <th>0</th>
          <td>2019-01-17 07:13:23</td>
          <td>2019-01-17 07:15:45</td>
          <td>西藏中路</td>
          <td>1</td>
        </tr>
        <tr>
          <th>1</th>
          <td>2019-01-17 07:34:24</td>
          <td>2019-01-17 07:35:38</td>
          <td>西藏中路</td>
          <td>1</td>
        </tr>
        <tr>
          <th>2</th>
          <td>2019-01-17 09:47:03</td>
          <td>2019-01-17 09:50:22</td>
          <td>西藏中路</td>
          <td>1</td>
        </tr>
        <tr>
          <th>...</th>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
        </tr>
        <tr>
          <th>2</th>
          <td>2019-01-17 16:35:52</td>
          <td>2019-01-17 16:36:49</td>
          <td>吴宝路</td>
          <td>148</td>
        </tr>
        <tr>
          <th>3</th>
          <td>2019-01-17 19:21:09</td>
          <td>2019-01-17 19:23:44</td>
          <td>吴宝路</td>
          <td>148</td>
        </tr>
        <tr>
          <th>0</th>
          <td>2019-01-17 13:36:26</td>
          <td>2019-01-17 13:45:04</td>
          <td>申昆路枢纽站</td>
          <td>148</td>
        </tr>
        <tr>
          <th>1</th>
          <td>2019-01-17 15:52:26</td>
          <td>2019-01-17 16:32:46</td>
          <td>申昆路枢纽站</td>
          <td>148</td>
        </tr>
        <tr>
          <th>2</th>
          <td>2019-01-17 19:24:54</td>
          <td>2019-01-17 19:25:55</td>
          <td>申昆路枢纽站</td>
          <td>148</td>
        </tr>
      </tbody>
    </table>
    <p>8984 rows × 4 columns</p>
    </div>



单程耗时
--------

::

    onewaytime = tbd.busgps_onewaytime(arriveinfo,
                                       start = '延安东路外滩',
                                       end = '申昆路枢纽站',col = ['VehicleId','stopname'])

::

    ## 绘制耗时分布箱型图
    import numpy as np
    import matplotlib.pyplot as plt
    import seaborn as sns
    plt.rcParams['font.sans-serif']=['SimHei']
    plt.rcParams['font.serif'] = ['SimHei']
    plt.rcParams['axes.unicode_minus']=False
    fig     = plt.figure(1,(8,4),dpi = 250)    
    ax1      = plt.subplot(111)
    
    sns.boxplot(x = 'shour',y = onewaytime['duration']/60,hue = '方向',data = onewaytime)
    
    plt.ylabel('始发站至终点站耗时（分钟）')
    plt.xlabel('小时')
    plt.ylim(0)
    plt.show()




.. image:: output_16_0.png


运营车速
--------

::

    #转换坐标系为投影坐标系，方便后面计算距离
    line.crs = {'init':'epsg:4326'}
    line_2416 = line.to_crs(epsg = 2416)
    #公交线路数据里面的geometry
    lineshp = line_2416['geometry'].iloc[0]
    linename = line_2416['name'].iloc[0]
    lineshp


.. parsed-literal::

    /opt/anaconda3/lib/python3.8/site-packages/pyproj/crs/crs.py:53: FutureWarning: '+init=<authority>:<code>' syntax is deprecated. '<authority>:<code>' is the preferred initialization method. When making the change, be mindful of axis order changes: https://pyproj4.github.io/pyproj/stable/gotchas.html#axis-order-changes-in-proj-6
      return _prepare_from_string(" ".join(pjargs))




.. image:: output_18_1.png



::

    #筛选去掉车速过快的
    #车速单位转换为km/h
    onewaytime['speed'] = (lineshp.length/onewaytime['duration'])*3.6
    onewaytime = onewaytime[onewaytime['speed']<=60]

::

    ## 车速分布
    import numpy as np
    import matplotlib.pyplot as plt
    import seaborn as sns
    plt.rcParams['font.sans-serif']=['SimHei']
    plt.rcParams['font.serif'] = ['SimHei']
    plt.rcParams['axes.unicode_minus']=False
    fig     = plt.figure(1,(8,4),dpi = 250)    
    ax1      = plt.subplot(111)
    sns.boxplot(x = 'shour',y = 'speed',hue = '方向',data = onewaytime)
    plt.ylabel('运营速度（km/h）')
    plt.xlabel('小时')
    plt.ylim(0)
    plt.show()



.. image:: output_20_0.png

pNEUMA轨迹数据处理
====================================

| 这个案例的Jupyter notebook: `点击这里 <https://github.com/ni1o1/transbigdata/blob/main/example/Example%204-pNEUMA%20trajectory%20dataset%20processing.ipynb>`__.
| 在这个例子中，我们将示例如何将``TransBigData``融入到雅典pNEUMA轨迹数据集的处理与可视化中。
| 请注意，样本数据已经被经过了一定的处理。原始版本的数据集可以在这里下载。
  `website <https://open-traffic.epfl.ch/>`__

::

    import transbigdata as tbd
    import pandas as pd
    import geopandas as gpd
    import matplotlib.pyplot as plt

读取数据
-------------

轨迹数据
~~~~~~~~~~~~~~~~~~~

::

    # 读取数据
    data = pd.read_csv('data/pNEUMA_tbd_sample.csv')
    # 将时间戳转换为时间格式
    data['time'] = pd.to_datetime(data['time'], unit='s')
    data.head()




.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>track_id</th>
          <th>lon</th>
          <th>lat</th>
          <th>speed</th>
          <th>time</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>0</th>
          <td>128</td>
          <td>23.730362</td>
          <td>37.990046</td>
          <td>12.5845</td>
          <td>1970-01-01 00:00:00.000</td>
        </tr>
        <tr>
          <th>1</th>
          <td>128</td>
          <td>23.730364</td>
          <td>37.990045</td>
          <td>12.4935</td>
          <td>1970-01-01 00:00:00.040</td>
        </tr>
        <tr>
          <th>2</th>
          <td>128</td>
          <td>23.730366</td>
          <td>37.990045</td>
          <td>12.3965</td>
          <td>1970-01-01 00:00:00.080</td>
        </tr>
        <tr>
          <th>3</th>
          <td>128</td>
          <td>23.730367</td>
          <td>37.990045</td>
          <td>12.2949</td>
          <td>1970-01-01 00:00:00.120</td>
        </tr>
        <tr>
          <th>4</th>
          <td>128</td>
          <td>23.730369</td>
          <td>37.990044</td>
          <td>12.1910</td>
          <td>1970-01-01 00:00:00.160</td>
        </tr>
      </tbody>
    </table>
    </div>





::

    # 输出数据大小信息
    data.info()


.. parsed-literal::

    <class 'pandas.core.frame.DataFrame'>
    RangeIndex: 581244 entries, 0 to 581243
    Data columns (total 5 columns):
     #   Column    Non-Null Count   Dtype         
    ---  ------    --------------   -----         
     0   track_id  581244 non-null  int64         
     1   lon       581244 non-null  float64       
     2   lat       581244 non-null  float64       
     3   speed     581244 non-null  float64       
     4   time      581244 non-null  datetime64[ns]
    dtypes: datetime64[ns](1), float64(3), int64(1)
    memory usage: 22.2 MB


OSM路网数据获取
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

你可以直接从``data``文件夹加载道路数据，或者使用OSMNX下载路网 `OSMNX <https://osmnx.readthedocs.io/en/stable/>`__

::

    # 从OSMNX中获取路网数据
    # OSM Graph
    import osmnx as ox
    bounds = [23.723577, 37.975462, 23.738471, 37.993053]
    north, south, east, west = bounds[3], bounds[1], bounds[2], bounds[0]
    G = ox.graph_from_bbox(north, south, east, west, network_type='drive')
    
    # 获取点和边
    nodes, edges = ox.graph_to_gdfs(G, nodes=True, edges=True)
    
    # 存储路网数据
    filepath = "data/pNEUMA_network.graphml"
    ox.save_graphml(G, filepath)

如果你没有OSMNX,可以运行下面代码读取已经现成的数据

::

    # 读取OSM数据
    import osmnx as ox
    filepath = "data/pNEUMA_network.graphml"
    G = ox.load_graphml(filepath)
    # 获取点和边
    nodes, edges = ox.graph_to_gdfs(G, nodes=True, edges=True)

地图底图加载
~~~~~~~~~~~~~~~~~~~~~

将地图底图加载并可视化

::

    # 可视化地图底图 tbd.plot_map
    bounds = [23.723577, 37.975462, 23.738471, 37.993053]
    
    fig = plt.figure(1, (12, 8), dpi=100)
    ax = plt.subplot(121)
    plt.sca(ax)
    tbd.plot_map(plt, bounds, zoom=18, style=1) # the map
    edges.plot(ax=ax, lw=1, color='grey') # edges
    nodes.plot(ax=ax, markersize = 8, color='red') # nodes
    plt.axis('off');
    
    ax = plt.subplot(122)
    plt.sca(ax)
    tbd.plot_map(plt, bounds, zoom=18, style=5) # the map
    edges.plot(ax=ax, lw=1, color='grey') # edges
    nodes.plot(ax=ax, markersize = 8, color='red') # nodes
    plt.axis('off');



.. image:: output_11_0.png


数据清洗
-------------

数据稀疏化
~~~~~~~~~~~~~~

| 数据集的采样间隔为 :math:`0.04` 秒, 非常小，不便于处理。
| 然而，一些宏观层面的研究不需要如此高的采样间隔。在这种情况下，数据可以使用``tbd.traj_sparsify``进行稀疏化。

::

    # 原始数据
    data.info()


.. parsed-literal::

    <class 'pandas.core.frame.DataFrame'>
    RangeIndex: 581244 entries, 0 to 581243
    Data columns (total 5 columns):
     #   Column    Non-Null Count   Dtype         
    ---  ------    --------------   -----         
     0   track_id  581244 non-null  int64         
     1   lon       581244 non-null  float64       
     2   lat       581244 non-null  float64       
     3   speed     581244 non-null  float64       
     4   time      581244 non-null  datetime64[ns]
    dtypes: datetime64[ns](1), float64(3), int64(1)
    memory usage: 22.2 MB


::

    #轨迹稀疏化
    data_sparsify = tbd.traj_sparsify(data, col=['track_id', 'time', 'lon', 'lat'],timegap=0.4,method='subsample')
    data_sparsify.info()


.. parsed-literal::

    <class 'pandas.core.frame.DataFrame'>
    Int64Index: 23293 entries, 0 to 581229
    Data columns (total 5 columns):
     #   Column    Non-Null Count  Dtype         
    ---  ------    --------------  -----         
     0   track_id  23293 non-null  int64         
     1   lon       23293 non-null  float64       
     2   lat       23293 non-null  float64       
     3   speed     23293 non-null  float64       
     4   time      23293 non-null  datetime64[ns]
    dtypes: datetime64[ns](1), float64(3), int64(1)
    memory usage: 1.1 MB


冗余数据剔除
~~~~~~~~~~~~~

在车辆停止运行时，位置没有发生移动，但仍然会产生大量GPS点，这些静止的GPS点除第一和最后一个点外的都可以删除。

::

    #用 tbd.clean_same 删除冗余数据
    data_sparsify_clean = tbd.clean_same(data_sparsify, col=['track_id', 'time', 'lon', 'lat'])
    data_sparsify_clean.info()


.. parsed-literal::

    <class 'pandas.core.frame.DataFrame'>
    Int64Index: 10674 entries, 0 to 581229
    Data columns (total 5 columns):
     #   Column    Non-Null Count  Dtype         
    ---  ------    --------------  -----         
     0   track_id  10674 non-null  int64         
     1   lon       10674 non-null  float64       
     2   lat       10674 non-null  float64       
     3   speed     10674 non-null  float64       
     4   time      10674 non-null  datetime64[ns]
    dtypes: datetime64[ns](1), float64(3), int64(1)
    memory usage: 500.3 KB


::

    data_sparsify_clean.head()




.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>track_id</th>
          <th>lon</th>
          <th>lat</th>
          <th>speed</th>
          <th>time</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>0</th>
          <td>128</td>
          <td>23.730362</td>
          <td>37.990046</td>
          <td>12.5845</td>
          <td>1970-01-01 00:00:00</td>
        </tr>
        <tr>
          <th>25</th>
          <td>128</td>
          <td>23.730399</td>
          <td>37.990040</td>
          <td>10.6835</td>
          <td>1970-01-01 00:00:01</td>
        </tr>
        <tr>
          <th>50</th>
          <td>128</td>
          <td>23.730429</td>
          <td>37.990036</td>
          <td>7.8580</td>
          <td>1970-01-01 00:00:02</td>
        </tr>
        <tr>
          <th>75</th>
          <td>128</td>
          <td>23.730443</td>
          <td>37.990033</td>
          <td>1.2661</td>
          <td>1970-01-01 00:00:03</td>
        </tr>
        <tr>
          <th>1775</th>
          <td>128</td>
          <td>23.730443</td>
          <td>37.990033</td>
          <td>0.0027</td>
          <td>1970-01-01 00:01:11</td>
        </tr>
      </tbody>
    </table>
    </div>



数据可视化
------------------

::

    gdf_data = gpd.GeoDataFrame(data_sparsify_clean, 
                                geometry=gpd.points_from_xy(data_sparsify_clean['lon'], 
                                                            data_sparsify_clean['lat']), 
                                crs=4326)
    gdf_data.head()




.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>track_id</th>
          <th>lon</th>
          <th>lat</th>
          <th>speed</th>
          <th>time</th>
          <th>geometry</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>0</th>
          <td>128</td>
          <td>23.730362</td>
          <td>37.990046</td>
          <td>12.5845</td>
          <td>1970-01-01 00:00:00</td>
          <td>POINT (23.73036 37.99005)</td>
        </tr>
        <tr>
          <th>25</th>
          <td>128</td>
          <td>23.730399</td>
          <td>37.990040</td>
          <td>10.6835</td>
          <td>1970-01-01 00:00:01</td>
          <td>POINT (23.73040 37.99004)</td>
        </tr>
        <tr>
          <th>50</th>
          <td>128</td>
          <td>23.730429</td>
          <td>37.990036</td>
          <td>7.8580</td>
          <td>1970-01-01 00:00:02</td>
          <td>POINT (23.73043 37.99004)</td>
        </tr>
        <tr>
          <th>75</th>
          <td>128</td>
          <td>23.730443</td>
          <td>37.990033</td>
          <td>1.2661</td>
          <td>1970-01-01 00:00:03</td>
          <td>POINT (23.73044 37.99003)</td>
        </tr>
        <tr>
          <th>1775</th>
          <td>128</td>
          <td>23.730443</td>
          <td>37.990033</td>
          <td>0.0027</td>
          <td>1970-01-01 00:01:11</td>
          <td>POINT (23.73044 37.99003)</td>
        </tr>
      </tbody>
    </table>
    </div>



::

    # 获取有最多数据点的车辆
    gdf_count = gdf_data.groupby('track_id')['lon'].count().rename('count').sort_values(ascending=False).reset_index()
    print(list(gdf_count.iloc[:20]['track_id']))


.. parsed-literal::

    [2138, 3290, 1442, 3197, 4408, 1767, 5002, 5022, 2140, 347, 2584, 4750, 4542, 2431, 4905, 4997, 1329, 4263, 1215, 3400]


可视化车辆

::

    fig = plt.figure(1, (6, 8), dpi=100)
    
    ax = plt.subplot(111)
    plt.sca(ax)
    
    # map
    tbd.plot_map(plt, bounds, zoom=18, style=4) # the map
    edges.plot(ax=ax, lw=1, color='grey') # edges
    # nodes.plot(ax=ax, markersize = 6, color='red') # nodes
    
    # trajectory
    gdf_data.plot(column='speed', ax=ax, markersize=0.5)
    
    plt.axis('off');



.. image:: output_22_0.png


可视化单辆车，并显示最短路径

::

    # select a vehicle
    tmpgdf_data = gdf_data[gdf_data['track_id']==2138]
    
    # the origin / destination location
    # o_point = [tmpgdf_data.iloc[0]['lon'], tmpgdf_data.iloc[0]['lat']]
    # d_point = [tmpgdf_data.iloc[-1]['lon'], tmpgdf_data.iloc[-1]['lat']]
    
    # get the nearest node of each point on the map
    tmpgdf_data = tbd.ckdnearest_point(tmpgdf_data, nodes)
    
    # extract the o/d node
    o_index, d_index = tmpgdf_data.iloc[0]['index'], tmpgdf_data.iloc[-1]['index']
    o_node_id, d_node_id = list(nodes[nodes['index']==o_index].index)[0], \
                           list(nodes[nodes['index']==d_index].index)[0]
    print(o_node_id, d_node_id)
    
    tmpgdf_data.head()


.. parsed-literal::

    250691723 358465943


.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>track_id</th>
          <th>lon</th>
          <th>lat</th>
          <th>speed</th>
          <th>time</th>
          <th>geometry_x</th>
          <th>dist</th>
          <th>index</th>
          <th>y</th>
          <th>x</th>
          <th>street_count</th>
          <th>highway</th>
          <th>geometry_y</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>0</th>
          <td>2138</td>
          <td>23.735287</td>
          <td>37.977435</td>
          <td>42.1006</td>
          <td>1970-01-01 00:01:35.560</td>
          <td>POINT (23.73529 37.97743)</td>
          <td>0.000779</td>
          <td>145</td>
          <td>37.978086</td>
          <td>23.734859</td>
          <td>4</td>
          <td>NaN</td>
          <td>POINT (23.73486 37.97809)</td>
        </tr>
        <tr>
          <th>1</th>
          <td>2138</td>
          <td>23.735254</td>
          <td>37.977473</td>
          <td>41.8663</td>
          <td>1970-01-01 00:01:36.000</td>
          <td>POINT (23.73525 37.97747)</td>
          <td>0.000729</td>
          <td>145</td>
          <td>37.978086</td>
          <td>23.734859</td>
          <td>4</td>
          <td>NaN</td>
          <td>POINT (23.73486 37.97809)</td>
        </tr>
        <tr>
          <th>2</th>
          <td>2138</td>
          <td>23.735181</td>
          <td>37.977558</td>
          <td>39.9012</td>
          <td>1970-01-01 00:01:37.000</td>
          <td>POINT (23.73518 37.97756)</td>
          <td>0.000618</td>
          <td>145</td>
          <td>37.978086</td>
          <td>23.734859</td>
          <td>4</td>
          <td>NaN</td>
          <td>POINT (23.73486 37.97809)</td>
        </tr>
        <tr>
          <th>3</th>
          <td>2138</td>
          <td>23.735111</td>
          <td>37.977638</td>
          <td>37.7748</td>
          <td>1970-01-01 00:01:38.000</td>
          <td>POINT (23.73511 37.97764)</td>
          <td>0.000514</td>
          <td>145</td>
          <td>37.978086</td>
          <td>23.734859</td>
          <td>4</td>
          <td>NaN</td>
          <td>POINT (23.73486 37.97809)</td>
        </tr>
        <tr>
          <th>4</th>
          <td>2138</td>
          <td>23.735047</td>
          <td>37.977712</td>
          <td>33.8450</td>
          <td>1970-01-01 00:01:39.000</td>
          <td>POINT (23.73505 37.97771)</td>
          <td>0.000418</td>
          <td>145</td>
          <td>37.978086</td>
          <td>23.734859</td>
          <td>4</td>
          <td>NaN</td>
          <td>POINT (23.73486 37.97809)</td>
        </tr>
      </tbody>
    </table>
    </div>






.. parsed-literal::

    250691723 358465943


::

    fig = plt.figure(1, (6, 8), dpi=100)
    
    ax = plt.subplot(111)
    plt.sca(ax)
    
    # map
    tbd.plot_map(plt, bounds, zoom=18, style=4) # the map
    edges.plot(ax=ax, lw=1, color='grey') # edges
    # nodes.plot(ax=ax, markersize = 6, color='red') # nodes
    
    # trajectory
    gdf_data[gdf_data['track_id']==2138].plot(ax=ax, markersize=5, color='red')
    
    
    plt.axis('off');



.. image:: output_27_0.png


我们可以将轨迹数据与最短路径做比对.

::

    # the shortest path (optional)
    # ax = plt.subplot(122)
    # plt.sca(ax)
    route = ox.shortest_path(G, o_node_id, d_node_id, weight="length")
    plt, ax = ox.plot_graph_route(G, route, route_color="green", route_linewidth=8, node_size=0)



.. image:: output_29_0.png

地铁网络拓扑建模
================

这个案例的Jupyter notebook: `点击这里 <https://github.com/ni1o1/transbigdata/blob/main/example/Example%203-Modeling%20for%20subway%20network%20topology.ipynb>`__.

| 可以点击 `这个链接 <https://mybinder.org/v2/gh/ni1o1/transbigdata/9507de936806c34a4befd74aa9227b012569a6aa?urlpath=lab%2Ftree%2Fexample%2FExample%203-Modeling%20for%20subway%20network%20topology.ipynb>`__ 在线编辑器中尝试

下面的案例展示如何用TransBigData包抓取地铁线路，并构建地铁线网的拓扑网络模型

爬取地铁线路
------------

首先爬取地铁线路使用tbd.getbusdata方法，输入城市跟公交或地铁线路名称的关键词，即可获取到线路数据，坐标系为wgs84。

::

    import transbigdata as tbd
    line,stop = tbd.getbusdata('厦门',['1号线','2号线','3号线'])



获取城市id: 厦门成功
1号线成功
2号线成功
3号线成功


::

    line.plot()








.. image:: output_5_1.png


::

    stop.plot()








.. image:: output_6_1.png


轨道断面信息获取
----------------

tbd.split_subwayline方法可以用轨道站点切分轨道线路，得到轨道断面信息（这一步骤主要在地铁客流可视化中有用）

::

    metroline_splited = tbd.split_subwayline(line,stop)
    metroline_splited.plot(column = 'o_project')





.. image:: output_9_1.png


轨道网络拓扑模型构建
--------------------

同时我们也可以直接使用站点数据，构建地铁网络的拓扑结构模型，方便后续地铁出行路径的识别。这一功能依赖于networkx包。

::

    #构建拓扑模型
    import networkx as nx
    G = tbd.metro_network(stop)
    nx.draw(G)



.. image:: output_12_0.png

出租车数据处理
==============

| 这个案例的Jupyter notebook: `点击这里 <https://github.com/ni1o1/transbigdata/blob/main/example/Example%201-Taxi%20GPS%20data%20processing.ipynb>`__.
| 可以点击 `这个链接 <https://mybinder.org/v2/gh/ni1o1/transbigdata/d7d6fa33ff16440ba1698b10dd3cf3f76ff00abd?urlpath=lab%2Ftree%2Fexample%2FExample%201-Taxi%20GPS%20data%20processing.ipynb>`__ 在线编辑器中尝试
| 使用示例中的样例数据集在github仓库中，链接为：https://github.com/ni1o1/transbigdata/tree/main/example
| 下面我们介绍如何使用TransBigData包，调用其中的函数实现对出租车GPS数据的快速处理。
| 首先我们引入TransBigData包，并读取数据:

::

    import transbigdata as tbd
    import pandas as pd
    import geopandas as gpd
    #读取数据    
    data = pd.read_csv('TaxiData-Sample.csv',header = None) 
    data.columns = ['VehicleNum','Time','Lng','Lat','OpenStatus','Speed']    
    data




.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>VehicleNum</th>
          <th>Time</th>
          <th>Lng</th>
          <th>Lat</th>
          <th>OpenStatus</th>
          <th>Speed</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>0</th>
          <td>34745</td>
          <td>20:27:43</td>
          <td>113.806847</td>
          <td>22.623249</td>
          <td>1</td>
          <td>27</td>
        </tr>
        <tr>
          <th>1</th>
          <td>34745</td>
          <td>20:24:07</td>
          <td>113.809898</td>
          <td>22.627399</td>
          <td>0</td>
          <td>0</td>
        </tr>
        <tr>
          <th>2</th>
          <td>34745</td>
          <td>20:24:27</td>
          <td>113.809898</td>
          <td>22.627399</td>
          <td>0</td>
          <td>0</td>
        </tr>
        <tr>
          <th>3</th>
          <td>34745</td>
          <td>20:22:07</td>
          <td>113.811348</td>
          <td>22.628067</td>
          <td>0</td>
          <td>0</td>
        </tr>
        <tr>
          <th>4</th>
          <td>34745</td>
          <td>20:10:06</td>
          <td>113.819885</td>
          <td>22.647800</td>
          <td>0</td>
          <td>54</td>
        </tr>
        <tr>
          <th>...</th>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
        </tr>
        <tr>
          <th>544994</th>
          <td>28265</td>
          <td>21:35:13</td>
          <td>114.321503</td>
          <td>22.709499</td>
          <td>0</td>
          <td>18</td>
        </tr>
        <tr>
          <th>544995</th>
          <td>28265</td>
          <td>09:08:02</td>
          <td>114.322701</td>
          <td>22.681700</td>
          <td>0</td>
          <td>0</td>
        </tr>
        <tr>
          <th>544996</th>
          <td>28265</td>
          <td>09:14:31</td>
          <td>114.336700</td>
          <td>22.690100</td>
          <td>0</td>
          <td>0</td>
        </tr>
        <tr>
          <th>544997</th>
          <td>28265</td>
          <td>21:19:12</td>
          <td>114.352600</td>
          <td>22.728399</td>
          <td>0</td>
          <td>0</td>
        </tr>
        <tr>
          <th>544998</th>
          <td>28265</td>
          <td>19:08:06</td>
          <td>114.137703</td>
          <td>22.621700</td>
          <td>0</td>
          <td>0</td>
        </tr>
      </tbody>
    </table>
    <p>544999 rows × 6 columns</p>
    </div>



::

    #读取区域信息
    import geopandas as gpd
    sz = gpd.read_file(r'sz/sz.shp')
    sz.crs = None
    sz.plot()





.. image:: output_3_1.png


数据预处理
----------------

TransBigData包也集成了数据预处理的常用方法。其中，tbd.clean_outofshape方法输入数据和研究范围区域信息，筛选剔除研究范围外的数据。而tbd.clean_taxi_status方法则可以剔除的载客状态瞬间变化的记录。在使用预处理的方法时，需要传入相应的列，代码如下：

::

    #数据预处理
    #剔除研究范围外的数据
    data = tbd.clean_outofshape(data, sz, col=['Lng', 'Lat'], accuracy=500)
    #剔除出租车数据中载客状态瞬间变化的记录
    data = tbd.clean_taxi_status(data, col=['VehicleNum', 'Time', 'OpenStatus'])

数据栅格化
----------------------------

以栅格形式表达数据分布是最基本的表达方法。GPS数据经过栅格化后，每个数据点都含有对应的栅格信息，采用栅格表达数据的分布时，其表示的分布情况与真实情况接近。如果要使用TransBigData工具进行栅格划分，首先需要确定栅格化的参数（可以理解为定义了一个栅格坐标系），参数可以帮助我们快速进行栅格化:

::

    #栅格化
    #定义范围，获取栅格化参数
    bounds = [113.6,22.4,114.8,22.9]
    params = tbd.grid_params(bounds,accuracy = 500)
    params

(113.6, 22.4, 0.004872390756896538, 0.004496605206422906)



取得栅格化参数后，将GPS对应至栅格，由LONCOL与LATCOL两列共同指定一个栅格:

::

    #将GPS栅格化
    data['LONCOL'],data['LATCOL'] = tbd.GPS_to_grids(data['Lng'],data['Lat'],params)

统计每个栅格的数据量:

::

    #集计栅格数据量
    datatest = data.groupby(['LONCOL','LATCOL'])['VehicleNum'].count().reset_index()

生成栅格的地理图形，并将它转化为GeoDataFrame:

::

    #生成栅格地理图形
    datatest['geometry'] = tbd.gridid_to_polygon(datatest['LONCOL'],datatest['LATCOL'],params)
    #转为GeoDataFrame
    import geopandas as gpd
    datatest = gpd.GeoDataFrame(datatest)


绘制栅格测试是否成功:

::

    #绘制
    datatest.plot(column = 'VehicleNum')



.. image:: output_17_1.png


出行OD提取与集计
----------------------

使用tbd.taxigps_to_od方法，传入对应的列名，即可提取出行OD:

::

    #从GPS数据提取OD
    oddata = tbd.taxigps_to_od(data,col = ['VehicleNum','Time','Lng','Lat','OpenStatus'])
    oddata




.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>VehicleNum</th>
          <th>stime</th>
          <th>slon</th>
          <th>slat</th>
          <th>etime</th>
          <th>elon</th>
          <th>elat</th>
          <th>ID</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>427075</th>
          <td>22396</td>
          <td>00:19:41</td>
          <td>114.013016</td>
          <td>22.664818</td>
          <td>00:23:01</td>
          <td>114.021400</td>
          <td>22.663918</td>
          <td>0</td>
        </tr>
        <tr>
          <th>131301</th>
          <td>22396</td>
          <td>00:41:51</td>
          <td>114.021767</td>
          <td>22.640200</td>
          <td>00:43:44</td>
          <td>114.026070</td>
          <td>22.640266</td>
          <td>1</td>
        </tr>
        <tr>
          <th>417417</th>
          <td>22396</td>
          <td>00:45:44</td>
          <td>114.028099</td>
          <td>22.645082</td>
          <td>00:47:44</td>
          <td>114.030380</td>
          <td>22.650017</td>
          <td>2</td>
        </tr>
        <tr>
          <th>376160</th>
          <td>22396</td>
          <td>01:08:26</td>
          <td>114.034897</td>
          <td>22.616301</td>
          <td>01:16:34</td>
          <td>114.035614</td>
          <td>22.646717</td>
          <td>3</td>
        </tr>
        <tr>
          <th>21768</th>
          <td>22396</td>
          <td>01:26:06</td>
          <td>114.046021</td>
          <td>22.641251</td>
          <td>01:34:48</td>
          <td>114.066048</td>
          <td>22.636183</td>
          <td>4</td>
        </tr>
        <tr>
          <th>...</th>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
        </tr>
        <tr>
          <th>57666</th>
          <td>36805</td>
          <td>22:37:42</td>
          <td>114.113403</td>
          <td>22.534767</td>
          <td>22:48:01</td>
          <td>114.114365</td>
          <td>22.550632</td>
          <td>5332</td>
        </tr>
        <tr>
          <th>175519</th>
          <td>36805</td>
          <td>22:49:12</td>
          <td>114.114365</td>
          <td>22.550632</td>
          <td>22:50:40</td>
          <td>114.115501</td>
          <td>22.557983</td>
          <td>5333</td>
        </tr>
        <tr>
          <th>212092</th>
          <td>36805</td>
          <td>22:52:07</td>
          <td>114.115402</td>
          <td>22.558083</td>
          <td>23:03:27</td>
          <td>114.118484</td>
          <td>22.547867</td>
          <td>5334</td>
        </tr>
        <tr>
          <th>119041</th>
          <td>36805</td>
          <td>23:03:45</td>
          <td>114.118484</td>
          <td>22.547867</td>
          <td>23:20:09</td>
          <td>114.133286</td>
          <td>22.617750</td>
          <td>5335</td>
        </tr>
        <tr>
          <th>224103</th>
          <td>36805</td>
          <td>23:36:19</td>
          <td>114.112968</td>
          <td>22.549601</td>
          <td>23:43:12</td>
          <td>114.089485</td>
          <td>22.538918</td>
          <td>5336</td>
        </tr>
      </tbody>
    </table>
    <p>5337 rows × 8 columns</p>
    </div>



对提取出的OD进行OD的栅格集计,并生成GeoDataFrame

::

    #栅格化OD并集计
    od_gdf = tbd.odagg_grid(oddata,params)
    od_gdf.plot(column = 'count')



.. image:: output_22_1.png


出行OD小区集计
--------------------------------

TransBigData包也提供了将OD直接集计到小区的方法

::

    #OD集计到小区（在不传入栅格化参数时，直接用经纬度匹配）
    od_gdf = tbd.odagg_shape(oddata,sz,round_accuracy=6)
    od_gdf.plot(column = 'count')





.. image:: output_25_1.png


::

    #OD集计到小区（传入栅格化参数时，先栅格化后匹配，可加快匹配速度，数据量大时建议使用）
    od_gdf = tbd.odagg_shape(oddata,sz,params = params)
    od_gdf.plot(column = 'count')




.. image:: output_26_1.png


基于matplotlib的地图绘制
------------------------------

tbd中提供了地图底图加载和比例尺指北针的功能。使用这个方法之前首先需要设置mapboxtoken和底图存放位置，详情看：\ `这个链接 <https://transbigdata.readthedocs.io/zh_CN/latest/plot_map.html>`__\ 。plot_map方法添加地图底图，plotscale添加比例尺和指北针:

::

    #创建图框
    import matplotlib.pyplot as plt
    import plot_map
    fig =plt.figure(1,(8,8),dpi=80)
    ax =plt.subplot(111)
    plt.sca(ax)
    #添加地图底图
    tbd.plot_map(plt,bounds,zoom = 12,style = 4)
    #绘制colorbar
    cax = plt.axes([0.05, 0.33, 0.02, 0.3])
    plt.title('count')
    plt.sca(ax)
    #绘制OD
    od_gdf.plot(ax = ax,vmax = 100,column = 'count',cax = cax,legend = True)
    #绘制小区底图
    sz.plot(ax = ax,edgecolor = (0,0,0,1),facecolor = (0,0,0,0.2),linewidths=0.5)
    #添加比例尺和指北针
    tbd.plotscale(ax,bounds = bounds,textsize = 10,compasssize = 1,accuracy = 2000,rect = [0.06,0.03],zorder = 10)
    plt.axis('off')
    plt.xlim(bounds[0],bounds[2])
    plt.ylim(bounds[1],bounds[3])
    plt.show()



.. image:: output_29_0.png


出租车轨迹的提取
----------------

使用tbd.taxigps_traj_point方法，输入数据和OD数据，可以提取出轨迹点

::

    data_deliver,data_idle = tbd.taxigps_traj_point(data,oddata,col=['VehicleNum', 'Time', 'Lng', 'Lat', 'OpenStatus'])

::

    data_deliver




.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>VehicleNum</th>
          <th>Time</th>
          <th>Lng</th>
          <th>Lat</th>
          <th>OpenStatus</th>
          <th>Speed</th>
          <th>LONCOL</th>
          <th>LATCOL</th>
          <th>ID</th>
          <th>flag</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>427075</th>
          <td>22396</td>
          <td>00:19:41</td>
          <td>114.013016</td>
          <td>22.664818</td>
          <td>1</td>
          <td>63.0</td>
          <td>85.0</td>
          <td>59.0</td>
          <td>0.0</td>
          <td>1.0</td>
        </tr>
        <tr>
          <th>427085</th>
          <td>22396</td>
          <td>00:19:49</td>
          <td>114.014030</td>
          <td>22.665483</td>
          <td>1</td>
          <td>55.0</td>
          <td>85.0</td>
          <td>59.0</td>
          <td>0.0</td>
          <td>1.0</td>
        </tr>
        <tr>
          <th>416622</th>
          <td>22396</td>
          <td>00:21:01</td>
          <td>114.018898</td>
          <td>22.662500</td>
          <td>1</td>
          <td>1.0</td>
          <td>86.0</td>
          <td>58.0</td>
          <td>0.0</td>
          <td>1.0</td>
        </tr>
        <tr>
          <th>427480</th>
          <td>22396</td>
          <td>00:21:41</td>
          <td>114.019348</td>
          <td>22.662300</td>
          <td>1</td>
          <td>7.0</td>
          <td>86.0</td>
          <td>58.0</td>
          <td>0.0</td>
          <td>1.0</td>
        </tr>
        <tr>
          <th>416623</th>
          <td>22396</td>
          <td>00:22:21</td>
          <td>114.020615</td>
          <td>22.663366</td>
          <td>1</td>
          <td>0.0</td>
          <td>86.0</td>
          <td>59.0</td>
          <td>0.0</td>
          <td>1.0</td>
        </tr>
        <tr>
          <th>...</th>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
        </tr>
        <tr>
          <th>170960</th>
          <td>36805</td>
          <td>23:42:31</td>
          <td>114.092766</td>
          <td>22.538317</td>
          <td>1</td>
          <td>66.0</td>
          <td>101.0</td>
          <td>31.0</td>
          <td>5336.0</td>
          <td>1.0</td>
        </tr>
        <tr>
          <th>170958</th>
          <td>36805</td>
          <td>23:42:37</td>
          <td>114.091721</td>
          <td>22.538349</td>
          <td>1</td>
          <td>65.0</td>
          <td>101.0</td>
          <td>31.0</td>
          <td>5336.0</td>
          <td>1.0</td>
        </tr>
        <tr>
          <th>170974</th>
          <td>36805</td>
          <td>23:42:43</td>
          <td>114.090752</td>
          <td>22.538300</td>
          <td>1</td>
          <td>60.0</td>
          <td>101.0</td>
          <td>31.0</td>
          <td>5336.0</td>
          <td>1.0</td>
        </tr>
        <tr>
          <th>170973</th>
          <td>36805</td>
          <td>23:42:49</td>
          <td>114.089813</td>
          <td>22.538099</td>
          <td>1</td>
          <td>62.0</td>
          <td>101.0</td>
          <td>31.0</td>
          <td>5336.0</td>
          <td>1.0</td>
        </tr>
        <tr>
          <th>253064</th>
          <td>36805</td>
          <td>23:42:55</td>
          <td>114.089500</td>
          <td>22.538067</td>
          <td>1</td>
          <td>51.0</td>
          <td>100.0</td>
          <td>31.0</td>
          <td>5336.0</td>
          <td>1.0</td>
        </tr>
      </tbody>
    </table>
    <p>190492 rows × 10 columns</p>
    </div>



::

    data_idle




.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>VehicleNum</th>
          <th>Time</th>
          <th>Lng</th>
          <th>Lat</th>
          <th>OpenStatus</th>
          <th>Speed</th>
          <th>LONCOL</th>
          <th>LATCOL</th>
          <th>ID</th>
          <th>flag</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>416628</th>
          <td>22396</td>
          <td>00:23:01</td>
          <td>114.021400</td>
          <td>22.663918</td>
          <td>0</td>
          <td>25.0</td>
          <td>86.0</td>
          <td>59.0</td>
          <td>0.0</td>
          <td>0.0</td>
        </tr>
        <tr>
          <th>401744</th>
          <td>22396</td>
          <td>00:25:01</td>
          <td>114.027115</td>
          <td>22.662100</td>
          <td>0</td>
          <td>25.0</td>
          <td>88.0</td>
          <td>58.0</td>
          <td>0.0</td>
          <td>0.0</td>
        </tr>
        <tr>
          <th>394630</th>
          <td>22396</td>
          <td>00:25:41</td>
          <td>114.024551</td>
          <td>22.659834</td>
          <td>0</td>
          <td>21.0</td>
          <td>87.0</td>
          <td>58.0</td>
          <td>0.0</td>
          <td>0.0</td>
        </tr>
        <tr>
          <th>394671</th>
          <td>22396</td>
          <td>00:26:21</td>
          <td>114.022797</td>
          <td>22.658367</td>
          <td>0</td>
          <td>0.0</td>
          <td>87.0</td>
          <td>57.0</td>
          <td>0.0</td>
          <td>0.0</td>
        </tr>
        <tr>
          <th>394672</th>
          <td>22396</td>
          <td>00:26:29</td>
          <td>114.022797</td>
          <td>22.658367</td>
          <td>0</td>
          <td>0.0</td>
          <td>87.0</td>
          <td>57.0</td>
          <td>0.0</td>
          <td>0.0</td>
        </tr>
        <tr>
          <th>...</th>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
        </tr>
        <tr>
          <th>64411</th>
          <td>36805</td>
          <td>23:53:09</td>
          <td>114.120354</td>
          <td>22.544300</td>
          <td>1</td>
          <td>2.0</td>
          <td>107.0</td>
          <td>32.0</td>
          <td>5336.0</td>
          <td>0.0</td>
        </tr>
        <tr>
          <th>64405</th>
          <td>36805</td>
          <td>23:53:15</td>
          <td>114.120354</td>
          <td>22.544300</td>
          <td>1</td>
          <td>1.0</td>
          <td>107.0</td>
          <td>32.0</td>
          <td>5336.0</td>
          <td>0.0</td>
        </tr>
        <tr>
          <th>64390</th>
          <td>36805</td>
          <td>23:53:21</td>
          <td>114.120354</td>
          <td>22.544300</td>
          <td>1</td>
          <td>0.0</td>
          <td>107.0</td>
          <td>32.0</td>
          <td>5336.0</td>
          <td>0.0</td>
        </tr>
        <tr>
          <th>64406</th>
          <td>36805</td>
          <td>23:53:27</td>
          <td>114.120354</td>
          <td>22.544300</td>
          <td>1</td>
          <td>0.0</td>
          <td>107.0</td>
          <td>32.0</td>
          <td>5336.0</td>
          <td>0.0</td>
        </tr>
        <tr>
          <th>64393</th>
          <td>36805</td>
          <td>23:53:33</td>
          <td>114.120354</td>
          <td>22.544300</td>
          <td>1</td>
          <td>0.0</td>
          <td>107.0</td>
          <td>32.0</td>
          <td>5336.0</td>
          <td>0.0</td>
        </tr>
      </tbody>
    </table>
    <p>312779 rows × 10 columns</p>
    </div>



对轨迹点生成载客与空载的轨迹

::

    traj_deliver = tbd.points_to_traj(data_deliver)
    traj_deliver.plot()




.. image:: output_36_1.png


::

    traj_idle = tbd.points_to_traj(data_idle)
    traj_idle.plot()

.. image:: output_37_1.png

轨迹可视化
------------------

| TransBigData包也依托于kepler.gl提供的可视化插件提供了一键数据整理与可视化的方法
| 使用此功能请先安装python的keplergl包


::

    pip install keplergl

将轨迹数据进行可视化：

::

    tbd.visualization_trip(data_deliver)

.. image:: kepler-traj.png
