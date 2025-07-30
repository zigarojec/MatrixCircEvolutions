# MatrixCircEvolutions
An evolutionary algorithm for automatic analog circuit topology synthesis. 

INTRODUCTION

Let us imagine the typical work-flow of an analog circuit designer. He or she has to:
  1. Figure out desired circuit specifications,
  2. choose a suitable topology to perform a task,
  3. find the physical sizes (numerical parameters) of choosen electrical building-blocks,
  4. evaluate and confirm the design. 
  
It is often a time consuming and ambiguous task, specially for un-experienced designer.

This software is build to help you with owerall analog circuit design. This software basically joins steps 2 and 3 together and makes possible for a designer to jump from the 1st step directly to the 4th step. 

INPUTS-OUTPUTS

There are two main inputs to this process. The first is a fixed bank of electrical building blocks (SPICE models) you want to be used in the final design. There can be more than needed, but certainly not less than needed. The second input is a high-level definition of the desired circuit behaviour. That is, for example, how much gain do we want, the bandwidth, dumping, allowed THD, etc.
The output is simply a netlist (or a set of netlists) that represent the resulting circuit (topology+parameters). 

IDEA

The main idea behing this software is to take the building-blocks bank and connect terminals together using special binary-connection matrix. Altering the setup of the binary matrix provides changes in the topology. We modify those matrices utilizing the artificial evolution process. Please read REFERENCES to fully understand the main core of this procedure. 

WARNING

Do not expect miracles. They might come, but do not count on it. When designing a circuit, firstly make sure that you know, what you are doing. This sofware can indeed find novel topologies and even return an optimized solution - but only when high-level circuit definition is defined well. WYWIWYG!

INSTALLATIONS

Before using this software you will have to follow the installations guide for PyOPUS, a framework for analog circuit optimization in Python. 
http://fides.fe.uni-lj.si/pyopus/quickstart.html

PREPARATION

Put whole project under a new folder _MAIN_work

Create _MAIN_data folder in the same folder as _MAIN_work 

USAGE

Set (choose) your scorefunctions.

Set globalVars.py

Run main.py for single-objective search.

Run main_moea.py for multi-objective search.

Run main_comb.py for search based on globalVars.py setting (single or multiobjective - automatically).



REFERENCES

[1] Ž. Rojec, Á. Bűrmen and I. Fajfar, "An evolution-driven analog circuit topology synthesis," 2016 IEEE Symposium Series on Computational Intelligence (SSCI), Athens, 2016, pp. 1-6.
doi: 10.1109/SSCI.2016.7850184
http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=7850184&isnumber=7849361

[2] Ž. Rojec, J. Olenšek, I. Fajfar, "Analog Circuit Topology Representation for Automated Synthesis and Optimization", Informacije MIDEM, Journal of Microelectronics, Electronic Components and MaterialsVol. 48, No. 1(2018), 29 – 40.
http://www.midem-drustvo.si/Journal%20papers/MIDEM_48(2018)1p29.pdf

[3] Ž. Rojec, Á. Bűrmen and I. Fajfar, "Analog Circuit Topology Synthesis by Means of Evolutionary Computation", Engineering Applications of Artificial Intelligence, ISSN 0952-1976. [Print ed.], Apr. 2019, vol. 80, str. 48-65, ilustr. https://www.sciencedirect.com/science/article/pii/S0952197619300119, doi: 10.1016/j.engappai.2019.01.012.

[4] Rojec, Žiga, Iztok Fajfar, and Árpád Burmen. 2022. "Evolutionary Synthesis of Failure-Resilient Analog Circuits" Mathematics 10, no. 1: 156. https://doi.org/10.3390/math10010156 

[5] Rojec, Žiga. "Towards smaller single-point failure-resilient analog circuits by use of a genetic algorithm." Informacije MIDEM 53, no. 2 (2023): 103-117.
