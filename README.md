McSplicer
=========

This repository provides a software that implements 
a probabilistic model for alternative splicing in Python. 
The model described in [here](https://github.com/shimlab/Probsplicing).


Using this software
-------------------

## Installation<a name="installation"></a>

```shell
git clone https://github.com/shimlab/McSplicer.git
```

You can execute the script:

```shell
python2 ./python_code/McSplicer.py --help
```


### Dependencies<a name="dependencies"></a>

McSplicer was implemnted and tested on python 2.7, and requires only few standard packages:
- numpy>=1.13.1
- pandas>=0.20.3

## Usage <a name="usage"></a>

Execute McSplicer script with `--help` option for a complete list of options.  
Sample data and usage examples can be found at `examples` subfolder.

### Output: ###

The output of EM algorithm is the likelihood of the parameters:


&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&Theta;=(&Pi;,p<sub>1</sub>,...,p<sub>Ms</sub>,q<sub>1</sub>,q<sub>2</sub>,...,q<sub>Me</sub>)
 
Where <i>Ms</i> and <i>Me</i> are the total number of start and end sites, respectively.


Authors
-------
* Israa Alqassem (alqassem.isra@gmail.com)
* Yash Kumar Sonthalia (yashsonthalia@iitj.ac.in)
* Heejung Shim (hjshim@gmail.com)
* Stefan Canzar (canzar@genzentrum.lmu.de)




&copy; 2017 McSplicer





