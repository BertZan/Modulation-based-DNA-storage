# A  robust and efficient DNA storage architecture based on modulation encoding and decoding 


This repository provides code for reproducing the results in the paper:

**``A  robust and efficient DNA storage architecture based on modulation encoding and decoding ''**, by 
Xiangzhen Zan, Ranze Xie,Xiangyu Yao,Peng Xu,Wenbin Liu 

Code by: Xiangzhen Zan (xiangcheng2436@163.com)



The aim of the code is to recover the data stored on DNA sequences from millions of noisy reads. The overall procedure comprises of i) indel detection of noise reads, ii)clustering the aligned reads and obataining the consensus sequences, iii) decoding the consensus sequences and recovering the data. Details about "setup and installation," is provided below.

### List of contents
[Setup and Installation](#Setup-and-Installation) <br>
[Data dependencies](#Data-dependencies) <br>
[Running the code and reproducing the results](#Running-the-code-and-reproducing-the-results) <br>
[License](#License)

# Setup and Installation
It takes less than 10 minutes to setup and install everything, provided that you specifically follow the instructions step by step.

### OS requirements
The code has been tested on the following operating system:

	Linux: Ubuntu 20.04.1 LTS

### Python dependencies
To reproduce the results, the following softwares are required. Assuming the experiment is being performed in a docker container or a linux machine, the following libraries and packages need to be installed.
        
        apt-get install python3.8       # --> or any other system-specific command for installing python3 on your system.		
		pip install numpy

# Running the code and reproducing the results

To reproduce the result, run the Modulation-based DNA storage_demo.py, for example by running the following code in the command line:

`python "Modulation-based DNA storage_demo.py"`

# License

This project is licensed under the GNU General Public License, version 3
(GPLv3), for more information read the LICENSE file or refer to:

  http://www.gnu.org/licenses/

# Citation

If you use our software, please cite:

Zan, X., R. Xie, X. Yao, P. Xu, and W. Liu (2022) A  robust and efficient DNA storage architecture based on modulation encoding and decoding. bioRxiv. 2022.2005.2025.490755.

https://www.biorxiv.org/content/10.1101/2022.05.25.490755v2
