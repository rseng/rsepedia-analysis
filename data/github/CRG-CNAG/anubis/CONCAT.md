----

<p align="center">
  <img src=".logo/anubis.png"/>
</p>

ANUBIS is a Python framework to analyze Tn-seq data.

## Requirements
Specific libraries are required by ANUBIS. We provide a [requirements](./requirements.txt) file to install everything at once. To do so, you will need first to have [pip](https://pip.pypa.io/en/stable/installing/) installed and then run:

```bash
pip3 --version                      # Check if installed
sudo apt-get install python3-pip    # if you need to install pip, you can check installation with the previous command
pip3 install -r requirements.txt
```

## Installation & Help

Download this repository and run:

```bash
python3 setup.py install
```

You may require to call it using sudo. Once installed, `anubis` will be integrated in your python distribution.

In the case you need to install the package in a specific directory of your system, you can call the argument *--install-lib* followed by a directory path:

```bash
python3 setup.py install --install-lib /custom/path/
```

## Examples

A short [manual](./Manual.ipynb) with examples can be visited from this same repository. We are currently writing a more in detail documentation. 

## Contact

This project has been fully developed at [Centre for Genomic Regulation](http://www.crg.eu/) at the group of [Design of Biological Systems](http://www.crg.eu/en/luis_serrano).

If you experience any problem at any step involving the program, you can use the 'Issues' page of this repository or contact:

[Miravet-Verde, Samuel](mailto:samuel.miravet@crg.eu)       
[Lluch-Senar, Maria](mailto:maria.lluch@crg.eu)       
[Serrano, Luis](mailto:luis.serrano@crg.eu)

## License

ANUBIS is under a common GNU GENERAL PUBLIC LICENSE. Plese, check [LICENSE](./LICENSE) for further information.

###### [2020] - Centre de Regulació Genòmica (CRG) - All Rights Reserved*

